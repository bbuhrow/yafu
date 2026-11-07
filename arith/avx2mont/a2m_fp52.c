/*
 * AVX2 + FMA Montgomery product with 52-bit digits: the AVX-512 template
 * (arith/vecmont/vecarith52_tmpl.c, AE_FP path) on 4 lanes.
 *
 * AVX2 has no embedded rounding, so the call runs with MXCSR set to round
 * toward zero (restored on return).  Products x*y < 2^104:
 *   t  = fma(x, y, 2^104)           bits(t)  = BIAS_HI + floor(x*y / 2^52)
 *   t' = fma(x, y, 2^104 + 2^52 - t)  bits(t') = BIAS_LO + (x*y mod 2^52)
 * Both exact under RZ.  Every other FP operation here is exact, so the
 * rounding mode only matters for the first FMA.  Build with
 * -frounding-math so the compiler keeps the FP ops between the MXCSR
 * writes.  -DA2M_CALLER_RZ: the caller sets RZ once (experiment).
 * uint64 -> double (< 2^52): bits(x | BIAS_LO) - 2^52 (no
 * vcvtuqq2pd in AVX2).
 */

#include "a2m.h"
#include <immintrin.h>
#include <string.h>

#define LD(p, i) _mm256_load_si256((const __m256i *)((p) + (size_t)(i) * 4))
#define ST(p, i, v) _mm256_store_si256((__m256i *)((p) + (size_t)(i) * 4), (v))
#define LDD(p, i) _mm256_castsi256_pd(LD(p, i))
#define BC(x) _mm256_set1_epi64x((long long)(x))
#define BIAS_HI 0x4670000000000000ULL   // bits of 2^104
#define BIAS_LO 0x4330000000000000ULL   // bits of 2^52
#define MASK52 0xfffffffffffffULL

#define MAC(l, h, x, y) do {                                                \
        __m256d _t = _mm256_fmadd_pd((x), (y), dhi);                        \
        (h) = _mm256_add_epi64((h), _mm256_castpd_si256(_t));               \
        _t = _mm256_fmadd_pd((x), (y), _mm256_sub_pd(dlo, _t));             \
        (l) = _mm256_add_epi64((l), _mm256_castpd_si256(_t));               \
    } while (0)

A2M_INLINE __m256d cvt52(__m256i x, __m256d d52)
{
    return _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(x, BC(BIAS_LO))), d52);
}

// conditional final subtraction c -= N in the lanes of tm (rare)
static A2M_NOINLINE void final_sub(uint64_t *C, const uint64_t *N, int tm, int S)
{
    int l, j;
    for (l = 0; l < 4; l++)
    {
        int64_t br = 0;
        if (!(tm & (1 << l)))
            continue;
        for (j = 0; j < S; j++)
        {
            int64_t t = (int64_t)C[j * 4 + l] - (int64_t)N[j * 4 + l] + br;
            C[j * 4 + l] = (uint64_t)t & MASK52;
            br = t >> 52;
        }
    }
}

A2M_INLINE __m256i col(const uint64_t *A, const uint64_t *B, const uint64_t *N,
    const uint64_t *M, __m256d mprev, __m256i hiprev, __m256i carry,
    __m256i *hk, int *np, const int k, const int S, const int SQR, const int FULL)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256d dhi = _mm256_castsi256_pd(BC(BIAS_HI));
    const __m256d dlo = _mm256_castsi256_pd(BC(BIAS_HI + 1));
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
            A2M_UNROLL_ALL
            for (i = ilo; i <= ihi; i++)
            {
                __m256d x = LDD(A, i), y = LDD(B, k - i);
                if ((i - ilo) & 1) MAC(la1, ha1, x, y); else MAC(la0, ha0, x, y);
            }
        }
        else
        {
            for (i = ilo; i + 1 <= ihi; i += 2)
            {
                MAC(la0, ha0, LDD(A, i), LDD(B, k - i));
                MAC(la1, ha1, LDD(A, i + 1), LDD(B, k - i - 1));
            }
            if (i <= ihi)
                MAC(la0, ha0, LDD(A, i), LDD(B, k - i));
        }
        n = ihi - ilo + 1;
    }
    else
    {
        int nd;
        if (FULL)
        {
            A2M_UNROLL_ALL
            for (i = ilo; 2 * i < k; i++)
            {
                __m256d x = LDD(A, i), y = LDD(A, k - i);
                if ((i - ilo) & 1) MAC(la1, ha1, x, y); else MAC(la0, ha0, x, y);
            }
        }
        else
        {
            for (i = ilo; 2 * (i + 1) < k; i += 2)
            {
                MAC(la0, ha0, LDD(A, i), LDD(A, k - i));
                MAC(la1, ha1, LDD(A, i + 1), LDD(A, k - i - 1));
            }
            if (2 * i < k)
                MAC(la0, ha0, LDD(A, i), LDD(A, k - i));
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
            __m256d x = LDD(A, k / 2);
            MAC(la1, ha1, x, x);
            n++;
        }
    }

    if (FULL)
    {
        A2M_UNROLL_ALL
        for (i = ilo; i <= mhi; i++)
        {
            __m256d x = (i == k - 1) ? mprev : LDD(M, i);
            __m256d y = LDD(N, k - i);
            if ((i - ilo) & 1) MAC(lm1, hm1, x, y); else MAC(lm0, hm0, x, y);
        }
    }
    else
    {
        int iend = (mhi == k - 1) ? mhi - 1 : mhi;
        for (i = ilo; i + 1 <= iend; i += 2)
        {
            MAC(lm0, hm0, LDD(M, i), LDD(N, k - i));
            MAC(lm1, hm1, LDD(M, i + 1), LDD(N, k - i - 1));
        }
        if (i <= iend)
            MAC(lm0, hm0, LDD(M, i), LDD(N, k - i));
        if (mhi == k - 1 && mhi >= ilo)
            MAC(lm1, hm1, mprev, LDD(N, 1));
    }
    if (mhi >= ilo)
        n += mhi - ilo + 1;
    *np = n;

    *hk = _mm256_add_epi64(_mm256_add_epi64(ha0, ha1), _mm256_add_epi64(hm0, hm1));
    return _mm256_add_epi64(
        _mm256_add_epi64(_mm256_add_epi64(la0, la1), _mm256_add_epi64(hiprev, carry)),
        _mm256_add_epi64(lm0, lm1));
}

A2M_INLINE void montmul_t(a2m_ctx_t *ctx, uint64_t *C, const uint64_t *A,
    const uint64_t *B, const int S, const int SQR, const int FULL)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i mask = BC(MASK52);
    const __m256i one = BC(1);
    const __m256d dhi = _mm256_castsi256_pd(BC(BIAS_HI));
    const __m256d dlo = _mm256_castsi256_pd(BC(BIAS_HI + 1));
    const __m256d d52 = _mm256_castsi256_pd(BC(BIAS_LO));
    uint64_t *W = ctx->W, *M = ctx->M;
    const uint64_t *N = ctx->N;
    const uint64_t *Ao, *Bo, *No;
    __m256d rhod, n0, mprev = _mm256_setzero_pd();
    __m256i carry = zero, hiprev = zero, hk, acc;
#ifndef A2M_CALLER_RZ
    unsigned int csr = _mm_getcsr();
#endif
    int j, k, np, tm;

#ifndef A2M_CALLER_RZ
    _mm_setcsr((csr & ~0x6000u) | 0x6000u);     // round toward zero
#endif

    rhod = cvt52(_mm256_loadu_si256((const __m256i *)ctx->rho), d52);
    for (j = 0; j < S; j++)
    {
        ST(W, j, _mm256_castpd_si256(cvt52(LD(A, j), d52)));
        ST(W + (size_t)2 * (size_t)S * 4, j, _mm256_castpd_si256(cvt52(LD(N, j), d52)));
    }
    if (!SQR)
        for (j = 0; j < S; j++)
            ST(W + (size_t)S * 4, j, _mm256_castpd_si256(cvt52(LD(B, j), d52)));
    Ao = W;
    Bo = SQR ? W : W + (size_t)S * 4;
    No = W + (size_t)2 * (size_t)S * 4;
    n0 = LDD(No, 0);

#define COLUMN_BODY                                                         \
    acc = col(Ao, Bo, No, M, mprev, hiprev, carry, &hk, &np, k, S, SQR, FULL); \
    acc = _mm256_sub_epi64(acc, BC((uint64_t)(unsigned)np * BIAS_LO));      \
    hk = _mm256_sub_epi64(hk, BC((uint64_t)(unsigned)np * BIAS_HI));        \
    if (k < S)                                                              \
    {                                                                       \
        __m256d a52 = cvt52(_mm256_and_si256(acc, mask), d52);              \
        __m256d mt = _mm256_fmadd_pd(a52, rhod, dhi);                       \
        mt = _mm256_fmadd_pd(a52, rhod, _mm256_sub_pd(dlo, mt));            \
        __m256d mk = _mm256_sub_pd(mt, d52);                                \
        ST(M, k, _mm256_castpd_si256(mk));                                  \
        hk = _mm256_add_epi64(hk, _mm256_sub_epi64(_mm256_castpd_si256(     \
            _mm256_fmadd_pd(mk, n0, dhi)), BC(BIAS_HI)));                   \
        mprev = mk;                                                         \
        /* acc + lo(m*n0) == 0 mod 2^52: carry (acc >> 52) + (lo52 != 0) */ \
        carry = _mm256_add_epi64(_mm256_add_epi64(_mm256_srli_epi64(acc, 52), one), \
            _mm256_cmpeq_epi64(_mm256_and_si256(acc, mask), zero));         \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        ST(C, k - S, _mm256_and_si256(acc, mask));                          \
        carry = _mm256_srli_epi64(acc, 52);                                 \
    }                                                                       \
    hiprev = hk;

    if (FULL)
    {
        A2M_UNROLL_ALL
        for (k = 0; k < 2 * S; k++)
        {
            COLUMN_BODY
        }
    }
    else
    {
        for (k = 0; k < S; k++)
        {
            COLUMN_BODY
        }
        for (k = S; k < 2 * S; k++)
        {
            COLUMN_BODY
        }
    }
#undef COLUMN_BODY

#ifndef A2M_CALLER_RZ
    _mm_setcsr(csr);
#endif
    tm = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(carry, zero)));
    if (tm)
        final_sub(C, N, tm, S);
}

#define INST(S)                                                             \
static void mul_s##S(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a, const uint64_t *b) \
{ montmul_t(ctx, c, a, b, S, 0, (S) <= A2M_FULL_MAX); }                    \
static void sqr_s##S(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a)       \
{ montmul_t(ctx, c, a, a, S, 1, (S) <= A2M_FULL_MAX); }

INST(2) INST(3) INST(4) INST(5) INST(6) INST(7) INST(8) INST(9)
INST(10) INST(11) INST(12) INST(13) INST(14) INST(15) INST(16) INST(17)
INST(18) INST(19) INST(20) INST(21) INST(22) INST(23) INST(24) INST(25)
INST(26) INST(27) INST(28) INST(29) INST(30) INST(31) INST(32)

static void mul_sN(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a, const uint64_t *b)
{ montmul_t(ctx, c, a, b, ctx->S, 0, 0); }
static void sqr_sN(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a)
{ montmul_t(ctx, c, a, a, ctx->S, 1, 0); }

#define P(S) { mul_s##S, sqr_s##S }
static const struct { a2m_mul_fn mul; a2m_sqr_fn sqr; } tab[33] = {
    { 0, 0 }, { 0, 0 }, P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9),
    P(10), P(11), P(12), P(13), P(14), P(15), P(16), P(17), P(18), P(19),
    P(20), P(21), P(22), P(23), P(24), P(25), P(26), P(27), P(28), P(29),
    P(30), P(31), P(32)
};

int a2m_get_fp52(int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr)
{
    if (S < 2)
        return 0;
    if (S <= 32) { *mul = tab[S].mul; *sqr = tab[S].sqr; }
    else { *mul = mul_sN; *sqr = sqr_sN; }
    return 1;
}
