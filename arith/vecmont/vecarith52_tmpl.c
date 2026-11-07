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
 *
 * Without IFMA (AVX-512F only: Skylake-X, Cascade Lake) the same template
 * runs on FP-emulated 52x52 products (AE_FP).  Operands are converted to
 * double once per call (A, B, N into mdata->fpwork, m_k into s as doubles),
 * and one product x*y < 2^104 yields both halves:
 *   t  = fma_rz(x, y, 2^104)          bits(t)  = BIAS_HI + floor(x*y / 2^52)
 *   t' = fma(x, y, 2^104 + 2^52 - t)  bits(t') = BIAS_LO + (x*y mod 2^52)
 * Both are exact.  The raw bit patterns are accumulated with integer adds
 * and the biases are removed once per column (count * bias, mod 2^64), so a
 * product costs 2 FMA + 1 FP sub + 2 integer adds.  Results are identical
 * to the IFMA version.
 */

#include "vecmont.h"

#if defined(USE_AVX512F)
#include <immintrin.h>

#ifdef IFMA
#define AE_FP 0
#else
#define AE_FP 1
#ifdef AE_TMPL_NO_CARRY_TRICK
#error "AE_TMPL_NO_CARRY_TRICK is IFMA-only"
#endif
#endif

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
#if !AE_FP
// operands: 52-bit integers
typedef __m512i ae_v;
#define AE_LDV(p, i) AE_LD(p, i)
#define AE_LO(acc, x, y) _mm512_madd52lo_epu64((acc), (x), (y))
#define AE_HI(acc, x, y) _mm512_madd52hi_epu64((acc), (x), (y))
// lo and hi halves of x*y added to l and h
#define AE_MAC(l, h, x, y) do { (l) = AE_LO((l), (x), (y)); (h) = AE_HI((h), (x), (y)); } while (0)
#else
// operands: the same values as doubles (exact, < 2^52), stored as raw bits
typedef __m512d ae_v;
#define AE_LDV(p, i) _mm512_castsi512_pd(AE_LD(p, i))
#define AE_RZ (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)
#define AE_BIAS_HI 0x4670000000000000ULL    // bits of 2^104
#define AE_BIAS_LO 0x4330000000000000ULL    // bits of 2^52
// needs ae_v dhi = 2^104, dlo = 2^104 + 2^52 in scope; adds biased bits
#define AE_MAC(l, h, x, y) do {                                             \
        __m512d _t = _mm512_fmadd_round_pd((x), (y), dhi, AE_RZ);           \
        (h) = _mm512_add_epi64((h), _mm512_castpd_si512(_t));               \
        _t = _mm512_fmadd_round_pd((x), (y), _mm512_sub_pd(dlo, _t), AE_RZ); \
        (l) = _mm512_add_epi64((l), _mm512_castpd_si512(_t));               \
    } while (0)
#endif

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
 * *np receives the number of products summed (FP: for the bias removal).
 */
AE_INLINE __m512i ae_col(const uint64_t *A, const uint64_t *B,
    const uint64_t *N, const uint64_t *M, const ae_v *Mr, ae_v mprev,
    __m512i hiprev, __m512i carry, __m512i *hk, int *np, const int k,
    const int S, const int SQR, const int MREG, const int FULL)
{
    const __m512i zero = _mm512_setzero_si512();
#if AE_FP
    const __m512d dhi = _mm512_castsi512_pd(_mm512_set1_epi64((long long)AE_BIAS_HI));
    const __m512d dlo = _mm512_castsi512_pd(_mm512_set1_epi64((long long)(AE_BIAS_HI + 1)));
#endif
    __m512i la0 = zero, ha0 = zero, la1 = zero, ha1 = zero;
    __m512i lm0 = zero, hm0 = zero, lm1 = zero, hm1 = zero;
    const int ilo = (k - (S - 1) > 0) ? k - (S - 1) : 0;
    const int ihi = (k < S - 1) ? k : S - 1;
    const int mhi = (k - 1 < S - 1) ? k - 1 : S - 1;    // last m index
    int i, n = 0;

    if (!SQR)
    {
        if (FULL)
        {
            AE_UNROLL_ALL
            for (i = ilo; i <= ihi; i++)
            {
                ae_v x = AE_LDV(A, i), y = AE_LDV(B, k - i);
                if ((i - ilo) & 1) AE_MAC(la1, ha1, x, y);
                else AE_MAC(la0, ha0, x, y);
            }
        }
        else
        {
            for (i = ilo; i + 1 <= ihi; i += 2)
            {
                ae_v x0 = AE_LDV(A, i), y0 = AE_LDV(B, k - i);
                ae_v x1 = AE_LDV(A, i + 1), y1 = AE_LDV(B, k - i - 1);
                AE_MAC(la0, ha0, x0, y0);
                AE_MAC(la1, ha1, x1, y1);
            }
            if (i <= ihi)
            {
                ae_v x = AE_LDV(A, i), y = AE_LDV(B, k - i);
                AE_MAC(la0, ha0, x, y);
            }
        }
        n = ihi - ilo + 1;
    }
    else
    {
        // off-diagonal products i < k - i, doubled, plus the diagonal
        int nd;
        if (FULL)
        {
            AE_UNROLL_ALL
            for (i = ilo; 2 * i < k; i++)
            {
                ae_v x = AE_LDV(A, i), y = AE_LDV(A, k - i);
                if ((i - ilo) & 1) AE_MAC(la1, ha1, x, y);
                else AE_MAC(la0, ha0, x, y);
            }
        }
        else
        {
            for (i = ilo; 2 * (i + 1) < k; i += 2)
            {
                ae_v x0 = AE_LDV(A, i), y0 = AE_LDV(A, k - i);
                ae_v x1 = AE_LDV(A, i + 1), y1 = AE_LDV(A, k - i - 1);
                AE_MAC(la0, ha0, x0, y0);
                AE_MAC(la1, ha1, x1, y1);
            }
            if (2 * i < k)
            {
                ae_v x = AE_LDV(A, i), y = AE_LDV(A, k - i);
                AE_MAC(la0, ha0, x, y);
            }
        }
        nd = (k + 1) / 2 - ilo;                 // i = ilo .. ceil(k/2) - 1
        n = (nd > 0) ? 2 * nd : 0;
        la0 = _mm512_add_epi64(la0, la1);
        ha0 = _mm512_add_epi64(ha0, ha1);
        la0 = _mm512_add_epi64(la0, la0);
        ha0 = _mm512_add_epi64(ha0, ha0);
        la1 = ha1 = zero;
        if (((k & 1) == 0) && (k / 2 <= S - 1))
        {
            ae_v x = AE_LDV(A, k / 2);
            AE_MAC(la1, ha1, x, x);
            n++;
        }
    }

    // reduction products m_i * n_(k-i), i = ilo..mhi; m_(k-1) is last
    // (it is the newest, on the critical path)
    if (FULL || MREG)
    {
        AE_UNROLL_ALL
        for (i = ilo; i <= mhi; i++)
        {
            ae_v x = (i == k - 1) ? mprev : (MREG ? Mr[i] : AE_LDV(M, i));
            ae_v y = AE_LDV(N, k - i);
            if ((i - ilo) & 1) AE_MAC(lm1, hm1, x, y);
            else AE_MAC(lm0, hm0, x, y);
        }
    }
    else
    {
        int iend = (mhi == k - 1) ? mhi - 1 : mhi;      // loop over stored m
        for (i = ilo; i + 1 <= iend; i += 2)
        {
            ae_v x0 = AE_LDV(M, i), y0 = AE_LDV(N, k - i);
            ae_v x1 = AE_LDV(M, i + 1), y1 = AE_LDV(N, k - i - 1);
            AE_MAC(lm0, hm0, x0, y0);
            AE_MAC(lm1, hm1, x1, y1);
        }
        if (i <= iend)
        {
            ae_v x = AE_LDV(M, i), y = AE_LDV(N, k - i);
            AE_MAC(lm0, hm0, x, y);
        }
        if (mhi == k - 1 && mhi >= ilo)
        {
            ae_v y = AE_LDV(N, 1);
            AE_MAC(lm1, hm1, mprev, y);
        }
    }
    if (mhi >= ilo)
        n += mhi - ilo + 1;
    *np = n;

    *hk = _mm512_add_epi64(_mm512_add_epi64(ha0, ha1), _mm512_add_epi64(hm0, hm1));
    return _mm512_add_epi64(
        _mm512_add_epi64(_mm512_add_epi64(la0, la1), _mm512_add_epi64(hiprev, carry)),
        _mm512_add_epi64(lm0, lm1));
}

AE_INLINE void ae_montmul_t(uint64_t *C, const uint64_t *A, const uint64_t *B,
    const uint64_t *N, const uint64_t *rho, uint64_t *M, uint64_t *W,
    const int S, const int SQR, const int FULL, const int MREG)
{
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64(0xfffffffffffffULL);
    AE_UNUSED const __m512i one = _mm512_set1_epi64(1);
    ae_v Mr[(AE_TMPL_MREG_MAX > 0) ? AE_TMPL_MREG_MAX : 1];
    __m512i carry = zero, hiprev = zero, hk, acc;
    const uint64_t *Ao = A, *Bo = B, *No = N;   // product operands
    __mmask8 tm;
    int k, np;
#if !AE_FP
    const __m512i vrho = _mm512_load_si512((const void *)rho);
    const __m512i n0 = AE_LD(N, 0);
    __m512i mprev = zero;
    (void)W;
#else
    const __m512d dhi = _mm512_castsi512_pd(_mm512_set1_epi64((long long)AE_BIAS_HI));
    const __m512d dlo = _mm512_castsi512_pd(_mm512_set1_epi64((long long)(AE_BIAS_HI + 1)));
    const __m512d d52 = _mm512_castsi512_pd(_mm512_set1_epi64((long long)AE_BIAS_LO));
    const __m512d rhod = _mm512_cvtepu64_pd(_mm512_load_si512((const void *)rho));
    __m512d mprev = _mm512_setzero_pd(), n0;
    int j;

    // convert the operands to double once (W: 3S vectors; B shares A's
    // slot when squaring).  A may alias C, so this must happen first.
    for (j = 0; j < S; j++)
    {
        AE_ST(W, j, _mm512_castpd_si512(_mm512_cvtepu64_pd(AE_LD(A, j))));
        AE_ST(W + (size_t)2 * (size_t)S * 8, j, _mm512_castpd_si512(_mm512_cvtepu64_pd(AE_LD(N, j))));
    }
    if (!SQR)
    {
        for (j = 0; j < S; j++)
            AE_ST(W + (size_t)S * 8, j, _mm512_castpd_si512(_mm512_cvtepu64_pd(AE_LD(B, j))));
    }
    Ao = W;
    Bo = SQR ? W : W + (size_t)S * 8;
    No = W + (size_t)2 * (size_t)S * 8;
    n0 = AE_LDV(No, 0);
#endif

#if !AE_FP
    // m_k = acc * rho mod 2^52 (IFMA uses the low 52 bits of acc)
#define AE_MK_STEP                                                          \
        __m512i mk = AE_LO(zero, acc, vrho);                                \
        if (MREG) Mr[k] = mk; else AE_ST(M, k, mk);                         \
        hk = AE_HI(hk, mk, n0);                                             \
        mprev = mk;
#define AE_FIX_BIAS
#else
#define AE_MK_STEP                                                          \
        __m512d a52 = _mm512_cvtepu64_pd(_mm512_and_si512(acc, mask));      \
        __m512d mt = _mm512_fmadd_round_pd(a52, rhod, dhi, AE_RZ);          \
        mt = _mm512_fmadd_round_pd(a52, rhod, _mm512_sub_pd(dlo, mt), AE_RZ); \
        __m512d mk = _mm512_sub_pd(mt, d52);        /* exact, < 2^52 */      \
        if (MREG) Mr[k] = mk; else AE_ST(M, k, _mm512_castpd_si512(mk));    \
        hk = _mm512_add_epi64(hk, _mm512_castpd_si512(                      \
            _mm512_fmadd_round_pd(mk, n0, dhi, AE_RZ)));                    \
        hk = _mm512_sub_epi64(hk, _mm512_set1_epi64((long long)AE_BIAS_HI)); \
        mprev = mk;
    // remove the biases of the np products summed in this column
#define AE_FIX_BIAS                                                         \
    acc = _mm512_sub_epi64(acc,                                             \
        _mm512_set1_epi64((long long)((unsigned long long)(unsigned)np * AE_BIAS_LO)));         \
    hk = _mm512_sub_epi64(hk,                                               \
        _mm512_set1_epi64((long long)((unsigned long long)(unsigned)np * AE_BIAS_HI)));
#endif

#define AE_COLUMN_BODY                                                      \
    acc = ae_col(Ao, Bo, No, M, Mr, mprev, hiprev, carry, &hk, &np, k, S,   \
        SQR, MREG, FULL);                                                   \
    AE_FIX_BIAS                                                             \
    if (k < S)                                                              \
    {                                                                       \
        AE_MK_STEP                                                          \
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
#undef AE_MK_STEP
#undef AE_FIX_BIAS

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
        mdata->fpwork, S, 0, AE_FULL(S), AE_MREG(S));                       \
    c->size = S;                                                            \
}                                                                           \
void vecsqrmod52_s##S(vec_bignum_t *a, vec_bignum_t *c,                     \
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)                   \
{                                                                           \
    c->signmask = 0;                                                        \
    ae_montmul_t(c->data, a->data, a->data, n->data, mdata->vrho, s->data,  \
        mdata->fpwork, S, 1, AE_FULL(S), AE_MREG(S));                       \
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
        mdata->fpwork, S, 0, 0, 0);
    c->size = S;
}

void vecsqrmod52_sN(vec_bignum_t *a, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    int S = (int)mdata->NWORDS;
    c->signmask = 0;
    ae_montmul_t(c->data, a->data, a->data, n->data, mdata->vrho, s->data,
        mdata->fpwork, S, 1, 0, 0);
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

#elif defined(USE_AVX2)

/*
 * AVX2 + FMA version (no AVX-512).  The same column algorithm with
 * FP-emulated products, on the same 8-lane vec_bignum_t layout: each call
 * processes lanes 0-3 and 4-7 as two 4-lane halves (stride 8 words).
 * AVX2 has no embedded rounding, so the call runs with MXCSR set to round
 * toward zero (restored on return; this file is built with -frounding-math
 * so the compiler keeps the FP ops between the MXCSR writes).  uint64 ->
 * double for values < 2^52 is bits(x | 2^52's bits) - 2^52 (no vcvtuqq2pd);
 * no vpminuq: the low-column carry uses vpcmpeqq.  Bit-identical to the
 * AVX-512 cores.  Measured in arith/avx2mont (arith/avx2_experiments.md):
 * full unrolling only pays up to S ~ 12 with 16 ymm registers.
 */
#include <immintrin.h>
// the MXCSR switch: keep the FP ops in order (gcc/clang: -frounding-math)
#if defined(_MSC_VER) && !defined(__clang__)
#pragma fenv_access(on)
#elif defined(__clang__)
#pragma STDC FENV_ACCESS ON
#endif

#ifndef AE2_FULL_MAX
#define AE2_FULL_MAX 12
#endif

#if defined(_MSC_VER) && !defined(__clang__)
#define AE2_INLINE static __forceinline
#define AE2_NOINLINE __declspec(noinline)
#define AE2_UNROLL_ALL
#else
#define AE2_INLINE static inline __attribute__((always_inline))
#define AE2_NOINLINE __attribute__((noinline))
#if defined(__clang__)
#define AE2_UNROLL_ALL _Pragma("unroll")
#else
#define AE2_UNROLL_ALL _Pragma("GCC unroll 256")
#endif
#endif

#define V2_LD(p, i) _mm256_load_si256((const __m256i *)(const void *)((p) + (size_t)(i) * 8))
#define V2_ST(p, i, v) _mm256_store_si256((__m256i *)(void *)((p) + (size_t)(i) * 8), (v))
#define V2_LDD(p, i) _mm256_castsi256_pd(V2_LD(p, i))
#define V2_BC(x) _mm256_set1_epi64x((long long)(x))
#define V2_BIAS_HI 0x4670000000000000ULL   // bits of 2^104
#define V2_BIAS_LO 0x4330000000000000ULL   // bits of 2^52
#define V2_MASK52 0xfffffffffffffULL

// both halves of x*y (MXCSR = round toward zero) as biased bit patterns
#define V2_MAC(l, h, x, y) do {                                             \
        __m256d _t = _mm256_fmadd_pd((x), (y), dhi);                        \
        (h) = _mm256_add_epi64((h), _mm256_castpd_si256(_t));               \
        _t = _mm256_fmadd_pd((x), (y), _mm256_sub_pd(dlo, _t));             \
        (l) = _mm256_add_epi64((l), _mm256_castpd_si256(_t));               \
    } while (0)

AE2_INLINE __m256d v2_cvt52(__m256i x, __m256d d52)
{
    return _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(x, V2_BC(V2_BIAS_LO))), d52);
}

// conditional final subtraction c -= N in the lanes of tm (rare); C and N
// point at this half's first lane
static AE2_NOINLINE void ae2_final_sub(uint64_t *C, const uint64_t *N, int tm, int S)
{
    int l, j;
    for (l = 0; l < 4; l++)
    {
        int64_t br = 0;
        if (!(tm & (1 << l)))
            continue;
        for (j = 0; j < S; j++)
        {
            int64_t t = (int64_t)C[j * 8 + l] - (int64_t)N[j * 8 + l] + br;
            C[j * 8 + l] = (uint64_t)t & V2_MASK52;
            br = t >> 52;
        }
    }
}

AE2_INLINE __m256i ae2_col(const uint64_t *A, const uint64_t *B, const uint64_t *N,
    const uint64_t *M, __m256d mprev, __m256i hiprev, __m256i carry,
    __m256i *hk, int *np, const int k, const int S, const int SQR, const int FULL)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256d dhi = _mm256_castsi256_pd(V2_BC(V2_BIAS_HI));
    const __m256d dlo = _mm256_castsi256_pd(V2_BC(V2_BIAS_HI + 1));
    __m256i la0 = zero, ha0 = zero, la1 = zero, ha1 = zero;
    __m256i lm0 = zero, hm0 = zero, lm1 = zero, hm1 = zero;
    const int ilo = (k - (S - 1) > 0) ? k - (S - 1) : 0;
    const int ihi = (k < S - 1) ? k : S - 1;
    const int mhi = (k - 1 < S - 1) ? k - 1 : S - 1;
    int i, n = 0;

    if (!SQR)
    {
        if (FULL)
        {
            AE2_UNROLL_ALL
            for (i = ilo; i <= ihi; i++)
            {
                __m256d x = V2_LDD(A, i), y = V2_LDD(B, k - i);
                if ((i - ilo) & 1) V2_MAC(la1, ha1, x, y); else V2_MAC(la0, ha0, x, y);
            }
        }
        else
        {
            for (i = ilo; i + 1 <= ihi; i += 2)
            {
                V2_MAC(la0, ha0, V2_LDD(A, i), V2_LDD(B, k - i));
                V2_MAC(la1, ha1, V2_LDD(A, i + 1), V2_LDD(B, k - i - 1));
            }
            if (i <= ihi)
                V2_MAC(la0, ha0, V2_LDD(A, i), V2_LDD(B, k - i));
        }
        n = ihi - ilo + 1;
    }
    else
    {
        int nd;
        if (FULL)
        {
            AE2_UNROLL_ALL
            for (i = ilo; 2 * i < k; i++)
            {
                __m256d x = V2_LDD(A, i), y = V2_LDD(A, k - i);
                if ((i - ilo) & 1) V2_MAC(la1, ha1, x, y); else V2_MAC(la0, ha0, x, y);
            }
        }
        else
        {
            for (i = ilo; 2 * (i + 1) < k; i += 2)
            {
                V2_MAC(la0, ha0, V2_LDD(A, i), V2_LDD(A, k - i));
                V2_MAC(la1, ha1, V2_LDD(A, i + 1), V2_LDD(A, k - i - 1));
            }
            if (2 * i < k)
                V2_MAC(la0, ha0, V2_LDD(A, i), V2_LDD(A, k - i));
        }
        nd = (k + 1) / 2 - ilo;
        n = (nd > 0) ? 2 * nd : 0;
        la0 = _mm256_add_epi64(la0, la1);
        ha0 = _mm256_add_epi64(ha0, ha1);
        la0 = _mm256_add_epi64(la0, la0);
        ha0 = _mm256_add_epi64(ha0, ha0);
        la1 = ha1 = zero;
        if (((k & 1) == 0) && (k / 2 <= S - 1))
        {
            __m256d x = V2_LDD(A, k / 2);
            V2_MAC(la1, ha1, x, x);
            n++;
        }
    }

    if (FULL)
    {
        AE2_UNROLL_ALL
        for (i = ilo; i <= mhi; i++)
        {
            __m256d x = (i == k - 1) ? mprev : V2_LDD(M, i);
            __m256d y = V2_LDD(N, k - i);
            if ((i - ilo) & 1) V2_MAC(lm1, hm1, x, y); else V2_MAC(lm0, hm0, x, y);
        }
    }
    else
    {
        int iend = (mhi == k - 1) ? mhi - 1 : mhi;
        for (i = ilo; i + 1 <= iend; i += 2)
        {
            V2_MAC(lm0, hm0, V2_LDD(M, i), V2_LDD(N, k - i));
            V2_MAC(lm1, hm1, V2_LDD(M, i + 1), V2_LDD(N, k - i - 1));
        }
        if (i <= iend)
            V2_MAC(lm0, hm0, V2_LDD(M, i), V2_LDD(N, k - i));
        if (mhi == k - 1 && mhi >= ilo)
            V2_MAC(lm1, hm1, mprev, V2_LDD(N, 1));
    }
    if (mhi >= ilo)
        n += mhi - ilo + 1;
    *np = n;

    *hk = _mm256_add_epi64(_mm256_add_epi64(ha0, ha1), _mm256_add_epi64(hm0, hm1));
    return _mm256_add_epi64(
        _mm256_add_epi64(_mm256_add_epi64(la0, la1), _mm256_add_epi64(hiprev, carry)),
        _mm256_add_epi64(lm0, lm1));
}

// one 4-lane half: all pointers at this half's first lane (stride 8);
// MXCSR is round-toward-zero (set by the caller)
AE2_INLINE void ae2_montmul_t(uint64_t *C, const uint64_t *A, const uint64_t *B,
    const uint64_t *N, const uint64_t *rho, uint64_t *M, uint64_t *W,
    const int S, const int SQR, const int FULL)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i mask = V2_BC(V2_MASK52);
    const __m256i one = V2_BC(1);
    const __m256d dhi = _mm256_castsi256_pd(V2_BC(V2_BIAS_HI));
    const __m256d dlo = _mm256_castsi256_pd(V2_BC(V2_BIAS_HI + 1));
    const __m256d d52 = _mm256_castsi256_pd(V2_BC(V2_BIAS_LO));
    const uint64_t *Ao, *Bo, *No;
    __m256d rhod, n0, mprev = _mm256_setzero_pd();
    __m256i carry = zero, hiprev = zero, hk, acc;
    int j, k, np, tm;

    rhod = v2_cvt52(_mm256_load_si256((const __m256i *)(const void *)rho), d52);
    for (j = 0; j < S; j++)
    {
        V2_ST(W, j, _mm256_castpd_si256(v2_cvt52(V2_LD(A, j), d52)));
        V2_ST(W + (size_t)2 * (size_t)S * 8, j, _mm256_castpd_si256(v2_cvt52(V2_LD(N, j), d52)));
    }
    if (!SQR)
        for (j = 0; j < S; j++)
            V2_ST(W + (size_t)S * 8, j, _mm256_castpd_si256(v2_cvt52(V2_LD(B, j), d52)));
    Ao = W;
    Bo = SQR ? W : W + (size_t)S * 8;
    No = W + (size_t)2 * (size_t)S * 8;
    n0 = V2_LDD(No, 0);

#define AE2_COLUMN_BODY                                                     \
    acc = ae2_col(Ao, Bo, No, M, mprev, hiprev, carry, &hk, &np, k, S, SQR, FULL); \
    acc = _mm256_sub_epi64(acc, V2_BC((unsigned long long)(unsigned)np * V2_BIAS_LO)); \
    hk = _mm256_sub_epi64(hk, V2_BC((unsigned long long)(unsigned)np * V2_BIAS_HI)); \
    if (k < S)                                                              \
    {                                                                       \
        __m256d a52 = v2_cvt52(_mm256_and_si256(acc, mask), d52);           \
        __m256d mt = _mm256_fmadd_pd(a52, rhod, dhi);                       \
        mt = _mm256_fmadd_pd(a52, rhod, _mm256_sub_pd(dlo, mt));            \
        __m256d mk = _mm256_sub_pd(mt, d52);                                \
        V2_ST(M, k, _mm256_castpd_si256(mk));                               \
        hk = _mm256_add_epi64(hk, _mm256_sub_epi64(_mm256_castpd_si256(     \
            _mm256_fmadd_pd(mk, n0, dhi)), V2_BC(V2_BIAS_HI)));             \
        mprev = mk;                                                         \
        carry = _mm256_add_epi64(_mm256_add_epi64(_mm256_srli_epi64(acc, 52), one), \
            _mm256_cmpeq_epi64(_mm256_and_si256(acc, mask), zero));         \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        V2_ST(C, k - S, _mm256_and_si256(acc, mask));                       \
        carry = _mm256_srli_epi64(acc, 52);                                 \
    }                                                                       \
    hiprev = hk;

    if (FULL)
    {
        AE2_UNROLL_ALL
        for (k = 0; k < 2 * S; k++)
        {
            AE2_COLUMN_BODY
        }
    }
    else
    {
        for (k = 0; k < S; k++)
        {
            AE2_COLUMN_BODY
        }
        for (k = S; k < 2 * S; k++)
        {
            AE2_COLUMN_BODY
        }
    }
#undef AE2_COLUMN_BODY

    tm = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(carry, zero)));
    if (tm)
        ae2_final_sub(C, N, tm, S);
}

// both halves, one MXCSR switch
AE2_INLINE void ae2_montmul(vec_bignum_t *c, vec_bignum_t *a, vec_bignum_t *b,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata, const int S,
    const int SQR, const int FULL)
{
    unsigned int csr = _mm_getcsr();
    int h;
    _mm_setcsr((csr & ~0x6000u) | 0x6000u);     // round toward zero
    for (h = 0; h < 8; h += 4)
        ae2_montmul_t(c->data + h, a->data + h, b->data + h, n->data + h,
            mdata->vrho + h, s->data + h, mdata->fpwork + h, S, SQR, FULL);
    _mm_setcsr(csr);
}

#define AE2_INST(S)                                                         \
void vecmulmod52_s##S(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,    \
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)                   \
{                                                                           \
    c->signmask = a->signmask ^ b->signmask;                                \
    ae2_montmul(c, a, b, n, s, mdata, S, 0, (S) <= AE2_FULL_MAX);           \
    c->size = S;                                                            \
}                                                                           \
void vecsqrmod52_s##S(vec_bignum_t *a, vec_bignum_t *c,                     \
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)                   \
{                                                                           \
    c->signmask = 0;                                                        \
    ae2_montmul(c, a, a, n, s, mdata, S, 1, (S) <= AE2_FULL_MAX);           \
    c->size = S;                                                            \
}

AE2_INST(2)  AE2_INST(3)  AE2_INST(4)  AE2_INST(5)
AE2_INST(6)  AE2_INST(7)  AE2_INST(8)  AE2_INST(9)
AE2_INST(10) AE2_INST(11) AE2_INST(12) AE2_INST(13)
AE2_INST(14) AE2_INST(15) AE2_INST(16) AE2_INST(17)
AE2_INST(18) AE2_INST(19) AE2_INST(20) AE2_INST(21)
AE2_INST(22) AE2_INST(23) AE2_INST(24) AE2_INST(25)
AE2_INST(26) AE2_INST(27) AE2_INST(28) AE2_INST(29)
AE2_INST(30) AE2_INST(31) AE2_INST(32)

void vecmulmod52_sN(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    int S = (int)mdata->NWORDS;
    c->signmask = a->signmask ^ b->signmask;
    ae2_montmul(c, a, b, n, s, mdata, S, 0, 0);
    c->size = S;
}

void vecsqrmod52_sN(vec_bignum_t *a, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    int S = (int)mdata->NWORDS;
    c->signmask = 0;
    ae2_montmul(c, a, a, n, s, mdata, S, 1, 0);
    c->size = S;
}

#define AE2_P(S) { vecmulmod52_s##S, vecsqrmod52_s##S }
static const struct { vecmont_mul_fn mul; vecmont_sqr_fn sqr; } ae2_tmpl_tab[33] = {
    { 0, 0 }, { 0, 0 }, AE2_P(2), AE2_P(3), AE2_P(4), AE2_P(5), AE2_P(6), AE2_P(7),
    AE2_P(8), AE2_P(9), AE2_P(10), AE2_P(11), AE2_P(12), AE2_P(13), AE2_P(14),
    AE2_P(15), AE2_P(16), AE2_P(17), AE2_P(18), AE2_P(19), AE2_P(20), AE2_P(21),
    AE2_P(22), AE2_P(23), AE2_P(24), AE2_P(25), AE2_P(26), AE2_P(27), AE2_P(28),
    AE2_P(29), AE2_P(30), AE2_P(31), AE2_P(32)
};

int vecmulmod52_tmpl_get(int S, vecmont_mul_fn *mul, vecmont_sqr_fn *sqr)
{
    if (S < 2)
        return 0;
    if (S <= 32)
    {
        *mul = ae2_tmpl_tab[S].mul;
        *sqr = ae2_tmpl_tab[S].sqr;
    }
    else
    {
        *mul = vecmulmod52_sN;
        *sqr = vecsqrmod52_sN;
    }
    return 1;
}

#endif
