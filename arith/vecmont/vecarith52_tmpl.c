/*
Copyright (c) 2026, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * Constant-S AVX-512 IFMA Montgomery mul/sqr for AVX-ECM (8 numbers, one
 * per lane, same modulus).  Drop-in replacements for vecmulmod52 /
 * vecsqrmod52 and the fixed416..1040 cores, for ANY word count S (not only
 * multiples of 4 words), generated from one always-inline template.
 *
 * Algorithm: column-wise (product scanning) Montgomery with lazy column
 * accumulation.  For column k the lo halves of the products a_i*b_(k-i)
 * and m_i*n_(k-i) accumulate in 64-bit registers, as do their hi halves;
 * the hi sums are added into column k+1.  52-bit words leave 12 spare bits,
 * so no per-product normalization is needed: one shift and two or three
 * adds per column.  m_k = acc_k * (-1/N) mod 2^52 depends only on the low
 * 52 bits of acc_k, so acc_k never needs normalizing before it is used.
 *
 * IFMA count (8 numbers): mul 4S^2 + S, sqr 3S^2 + 2S (the CIOS minimum).
 *
 * Output: almost-Montgomery (AMM) exactly like vecmulmod52 with USE_AMM:
 * c = (a*b + m*N) / R, minus N if that is >= R, so c < R and the result is
 * bit-identical to the existing kernels.  The conditional subtraction is a
 * rarely taken slow path behind a mask test.
 *
 * Inputs: a, b < R with normalized 52-bit words.  c may alias a or b.
 * s (scratch, >= S words) must not alias a, b or c.
 */

#include "vecmont.h"

#if defined(USE_AVX512F) && defined(IFMA)
#include <immintrin.h>

// Fully unroll the column loop for S <= this; above it the column loop is
// kept and only the inner dot-product loops are unrolled.
#ifndef AE_TMPL_FULL_MAX
#define AE_TMPL_FULL_MAX 24
#endif

// The low-column carry is computed without the m_k*n_0 lo product (exact;
// saves S IFMAs, 1-3% measured).  Define AE_TMPL_NO_CARRY_TRICK for the
// plain version.

// Keep m_i in registers (array scalarized by gcc) for S <= this, otherwise
// spill them to the scratch s and reload as memory operands.
#ifndef AE_TMPL_MREG_MAX
#define AE_TMPL_MREG_MAX 0
#endif

#if defined(_MSC_VER) && !defined(__clang__)
#define AE_INLINE static __forceinline
#define AE_NOINLINE __declspec(noinline)
#define AE_UNUSED __pragma(warning(suppress: 4100 4101))
#define AE_UNLIKELY(x) (x)
#define AE_UNROLL_ALL
#else
#define AE_INLINE static inline __attribute__((always_inline))
#define AE_NOINLINE __attribute__((noinline))
#define AE_UNUSED __attribute__((unused))
#define AE_UNLIKELY(x) __builtin_expect(!!(x), 0)
#if defined(__clang__)
#define AE_UNROLL_ALL _Pragma("unroll")
#else
#define AE_UNROLL_ALL _Pragma("GCC unroll 256")
#endif
#endif
#define AE_LD(p, i) _mm512_load_si512((const void *)((p) + (size_t)(i) * 8))
#define AE_ST(p, i, v) _mm512_store_si512((void *)((p) + (size_t)(i) * 8), (v))
#define AE_LO(acc, x, y) _mm512_madd52lo_epu64((acc), (x), (y))
#define AE_HI(acc, x, y) _mm512_madd52hi_epu64((acc), (x), (y))

#ifdef AE_TMPL_COUNT_SLOW
uint64_t ae_tmpl_slow_count = 0;
uint64_t ae_tmpl_calls = 0;
#endif

// conditional final subtraction c -= N in the lanes of tm (slow path)
static AE_NOINLINE void ae_tmpl_final_sub(uint64_t *C,
    const uint64_t *N, __mmask8 tm, int S)
{
    const __m512i mask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i br = _mm512_setzero_si512();
    int j;
    for (j = 0; j < S; j++)
    {
        __m512i t = _mm512_add_epi64(_mm512_sub_epi64(AE_LD(C, j), AE_LD(N, j)), br);
        br = _mm512_srai_epi64(t, 52);
        _mm512_mask_store_epi64((void *)(C + (size_t)j * 8), tm,
            _mm512_and_si512(t, mask));
    }
}

/*
 * One column of the Montgomery product.  Returns the column accumulator
 * (lo sums + previous hi sums + carry); *hk receives this column's hi sums.
 * mprev is m_(k-1) in a register (avoids a store-forward on the critical
 * path); M holds m_0..m_(k-2) (or Mr does, when register-resident).
 */
AE_INLINE __m512i ae_col(const uint64_t *A, const uint64_t *B,
    const uint64_t *N, const uint64_t *M, const __m512i *Mr, __m512i mprev,
    __m512i hiprev, __m512i carry, __m512i *hk, const int k, const int S,
    const int SQR, const int MREG, const int FULL)
{
    const __m512i zero = _mm512_setzero_si512();
    __m512i la0 = zero, ha0 = zero, la1 = zero, ha1 = zero;
    __m512i lm0 = zero, hm0 = zero, lm1 = zero, hm1 = zero;
    const int ilo = (k - (S - 1) > 0) ? k - (S - 1) : 0;
    const int ihi = (k < S - 1) ? k : S - 1;
    const int mhi = (k - 1 < S - 1) ? k - 1 : S - 1;    // last m index
    int i;

    if (!SQR)
    {
        if (FULL)
        {
            AE_UNROLL_ALL
            for (i = ilo; i <= ihi; i++)
            {
                __m512i x = AE_LD(A, i), y = AE_LD(B, k - i);
                if ((i - ilo) & 1) { la1 = AE_LO(la1, x, y); ha1 = AE_HI(ha1, x, y); }
                else { la0 = AE_LO(la0, x, y); ha0 = AE_HI(ha0, x, y); }
            }
        }
        else
        {
            for (i = ilo; i + 1 <= ihi; i += 2)
            {
                __m512i x0 = AE_LD(A, i), y0 = AE_LD(B, k - i);
                __m512i x1 = AE_LD(A, i + 1), y1 = AE_LD(B, k - i - 1);
                la0 = AE_LO(la0, x0, y0); ha0 = AE_HI(ha0, x0, y0);
                la1 = AE_LO(la1, x1, y1); ha1 = AE_HI(ha1, x1, y1);
            }
            if (i <= ihi)
            {
                __m512i x = AE_LD(A, i), y = AE_LD(B, k - i);
                la0 = AE_LO(la0, x, y); ha0 = AE_HI(ha0, x, y);
            }
        }
    }
    else
    {
        // off-diagonal products i < k - i, doubled, plus the diagonal
        if (FULL)
        {
            AE_UNROLL_ALL
            for (i = ilo; 2 * i < k; i++)
            {
                __m512i x = AE_LD(A, i), y = AE_LD(A, k - i);
                if ((i - ilo) & 1) { la1 = AE_LO(la1, x, y); ha1 = AE_HI(ha1, x, y); }
                else { la0 = AE_LO(la0, x, y); ha0 = AE_HI(ha0, x, y); }
            }
        }
        else
        {
            for (i = ilo; 2 * (i + 1) < k; i += 2)
            {
                __m512i x0 = AE_LD(A, i), y0 = AE_LD(A, k - i);
                __m512i x1 = AE_LD(A, i + 1), y1 = AE_LD(A, k - i - 1);
                la0 = AE_LO(la0, x0, y0); ha0 = AE_HI(ha0, x0, y0);
                la1 = AE_LO(la1, x1, y1); ha1 = AE_HI(ha1, x1, y1);
            }
            if (2 * i < k)
            {
                __m512i x = AE_LD(A, i), y = AE_LD(A, k - i);
                la0 = AE_LO(la0, x, y); ha0 = AE_HI(ha0, x, y);
            }
        }
        la0 = _mm512_add_epi64(la0, la1);
        ha0 = _mm512_add_epi64(ha0, ha1);
        la0 = _mm512_add_epi64(la0, la0);
        ha0 = _mm512_add_epi64(ha0, ha0);
        la1 = ha1 = zero;
        if (((k & 1) == 0) && (k / 2 <= S - 1))
        {
            __m512i x = AE_LD(A, k / 2);
            la1 = AE_LO(la1, x, x);
            ha1 = AE_HI(ha1, x, x);
        }
    }

    // reduction products m_i * n_(k-i), i = ilo..mhi; m_(k-1) is last
    // (it is the newest, on the critical path)
    if (FULL || MREG)
    {
        AE_UNROLL_ALL
        for (i = ilo; i <= mhi; i++)
        {
            __m512i x = (i == k - 1) ? mprev : (MREG ? Mr[i] : AE_LD(M, i));
            __m512i y = AE_LD(N, k - i);
            if ((i - ilo) & 1) { lm1 = AE_LO(lm1, x, y); hm1 = AE_HI(hm1, x, y); }
            else { lm0 = AE_LO(lm0, x, y); hm0 = AE_HI(hm0, x, y); }
        }
    }
    else
    {
        int iend = (mhi == k - 1) ? mhi - 1 : mhi;      // loop over stored m
        for (i = ilo; i + 1 <= iend; i += 2)
        {
            __m512i x0 = AE_LD(M, i), y0 = AE_LD(N, k - i);
            __m512i x1 = AE_LD(M, i + 1), y1 = AE_LD(N, k - i - 1);
            lm0 = AE_LO(lm0, x0, y0); hm0 = AE_HI(hm0, x0, y0);
            lm1 = AE_LO(lm1, x1, y1); hm1 = AE_HI(hm1, x1, y1);
        }
        if (i <= iend)
        {
            __m512i x = AE_LD(M, i), y = AE_LD(N, k - i);
            lm0 = AE_LO(lm0, x, y); hm0 = AE_HI(hm0, x, y);
        }
        if (mhi == k - 1 && mhi >= ilo)
        {
            __m512i y = AE_LD(N, 1);
            lm1 = AE_LO(lm1, mprev, y); hm1 = AE_HI(hm1, mprev, y);
        }
    }

    *hk = _mm512_add_epi64(_mm512_add_epi64(ha0, ha1), _mm512_add_epi64(hm0, hm1));
    return _mm512_add_epi64(
        _mm512_add_epi64(_mm512_add_epi64(la0, la1), _mm512_add_epi64(hiprev, carry)),
        _mm512_add_epi64(lm0, lm1));
}

AE_INLINE void ae_montmul_t(uint64_t *C, const uint64_t *A, const uint64_t *B,
    const uint64_t *N, const uint64_t *rho, uint64_t *M,
    const int S, const int SQR, const int FULL, const int MREG)
{
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64(0xfffffffffffffULL);
    const __m512i vrho = _mm512_load_si512((const void *)rho);
    const __m512i n0 = AE_LD(N, 0);
    AE_UNUSED const __m512i one = _mm512_set1_epi64(1);
    __m512i Mr[(AE_TMPL_MREG_MAX > 0) ? AE_TMPL_MREG_MAX : 1];
    __m512i carry = zero, hiprev = zero, mprev = zero, hk, acc;
    __mmask8 tm;
    int k;

#define AE_COLUMN_BODY                                                      \
    acc = ae_col(A, B, N, M, Mr, mprev, hiprev, carry, &hk, k, S, SQR,      \
        MREG, FULL);                                                        \
    if (k < S)                                                              \
    {                                                                       \
        __m512i mk = AE_LO(zero, acc, vrho);                                \
        if (MREG) Mr[k] = mk; else AE_ST(M, k, mk);                         \
        hk = AE_HI(hk, mk, n0);                                             \
        mprev = mk;                                                         \
        AE_LOW_CARRY                                                        \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        AE_ST(C, k - S, _mm512_and_si512(acc, mask));                       \
        carry = _mm512_srli_epi64(acc, 52);                                 \
    }                                                                       \
    hiprev = hk;

#ifndef AE_TMPL_NO_CARRY_TRICK
    // acc + lo(m*n0) == 0 mod 2^52, so its carry is (acc >> 52) + 1
    // unless acc's low 52 bits are already 0: no IFMA needed
#define AE_LOW_CARRY                                                        \
        carry = _mm512_add_epi64(_mm512_srli_epi64(acc, 52),                \
            _mm512_min_epu64(_mm512_and_si512(acc, mask), one));
#else
#define AE_LOW_CARRY                                                        \
        acc = AE_LO(acc, mk, n0);                                           \
        carry = _mm512_srli_epi64(acc, 52);
#endif

    if (FULL)
    {
        AE_UNROLL_ALL
        for (k = 0; k < 2 * S; k++)
        {
            AE_COLUMN_BODY
        }
    }
    else
    {
        for (k = 0; k < S; k++)
        {
            AE_COLUMN_BODY
        }
        for (k = S; k < 2 * S; k++)
        {
            AE_COLUMN_BODY
        }
    }
#undef AE_COLUMN_BODY
#undef AE_LOW_CARRY

    // carry out of column 2S-1 is the bit at R: result >= R, subtract N
    tm = _mm512_test_epi64_mask(carry, carry);
#ifdef AE_TMPL_COUNT_SLOW
    ae_tmpl_calls++;
    if (tm) ae_tmpl_slow_count++;
#endif
    if (AE_UNLIKELY(tm != 0))
        ae_tmpl_final_sub(C, N, tm, S);
}

#define AE_FULL(S) ((S) <= AE_TMPL_FULL_MAX)
#define AE_MREG(S) ((S) <= AE_TMPL_MREG_MAX)
#define AE_TMPL_INST(S)                                                     \
void vecmulmod52_s##S(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,    \
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)                   \
{                                                                           \
    c->signmask = a->signmask ^ b->signmask;                                \
    ae_montmul_t(c->data, a->data, b->data, n->data, mdata->vrho, s->data,  \
        S, 0, AE_FULL(S), AE_MREG(S));                                      \
    c->size = S;                                                            \
}                                                                           \
void vecsqrmod52_s##S(vec_bignum_t *a, vec_bignum_t *c,                     \
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)                   \
{                                                                           \
    c->signmask = 0;                                                        \
    ae_montmul_t(c->data, a->data, a->data, n->data, mdata->vrho, s->data,  \
        S, 1, AE_FULL(S), AE_MREG(S));                                      \
    c->size = S;                                                            \
}

AE_TMPL_INST(2)  AE_TMPL_INST(3)  AE_TMPL_INST(4)  AE_TMPL_INST(5)
AE_TMPL_INST(6)  AE_TMPL_INST(7)  AE_TMPL_INST(8)  AE_TMPL_INST(9)
AE_TMPL_INST(10) AE_TMPL_INST(11) AE_TMPL_INST(12) AE_TMPL_INST(13)
AE_TMPL_INST(14) AE_TMPL_INST(15) AE_TMPL_INST(16) AE_TMPL_INST(17)
AE_TMPL_INST(18) AE_TMPL_INST(19) AE_TMPL_INST(20) AE_TMPL_INST(21)
AE_TMPL_INST(22) AE_TMPL_INST(23) AE_TMPL_INST(24) AE_TMPL_INST(25)
AE_TMPL_INST(26) AE_TMPL_INST(27) AE_TMPL_INST(28) AE_TMPL_INST(29)
AE_TMPL_INST(30) AE_TMPL_INST(31) AE_TMPL_INST(32)

// general S (loop form) for anything larger
void vecmulmod52_sN(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    int S = (int)mdata->NWORDS;
    c->signmask = a->signmask ^ b->signmask;
    ae_montmul_t(c->data, a->data, b->data, n->data, mdata->vrho, s->data,
        S, 0, 0, 0);
    c->size = S;
}

void vecsqrmod52_sN(vec_bignum_t *a, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    int S = (int)mdata->NWORDS;
    c->signmask = 0;
    ae_montmul_t(c->data, a->data, a->data, n->data, mdata->vrho, s->data,
        S, 1, 0, 0);
    c->size = S;
}

#define AE_P(S) { vecmulmod52_s##S, vecsqrmod52_s##S }
static const struct { vecmont_mul_fn mul; vecmont_sqr_fn sqr; } ae_tmpl_tab[33] = {
    { 0, 0 }, { 0, 0 }, AE_P(2), AE_P(3), AE_P(4), AE_P(5), AE_P(6), AE_P(7),
    AE_P(8), AE_P(9), AE_P(10), AE_P(11), AE_P(12), AE_P(13), AE_P(14),
    AE_P(15), AE_P(16), AE_P(17), AE_P(18), AE_P(19), AE_P(20), AE_P(21),
    AE_P(22), AE_P(23), AE_P(24), AE_P(25), AE_P(26), AE_P(27), AE_P(28),
    AE_P(29), AE_P(30), AE_P(31), AE_P(32)
};

// pick the mul/sqr pair for NWORDS = S; returns 0 if S is unsupported
int vecmulmod52_tmpl_get(int S, vecmont_mul_fn *mul, vecmont_sqr_fn *sqr)
{
    if (S < 2)
        return 0;
    if (S <= 32)
    {
        *mul = ae_tmpl_tab[S].mul;
        *sqr = ae_tmpl_tab[S].sqr;
    }
    else
    {
        *mul = vecmulmod52_sN;
        *sqr = vecsqrmod52_sN;
    }
    return 1;
}

#endif
