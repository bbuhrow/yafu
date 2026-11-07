/*
 * AVX2 Montgomery product with integer digits and vpmuludq (32x32 -> 64,
 * 4 lanes): column-wise (product scanning), lazy column sums.
 *
 * D28: 28-bit digits.  A product is < 2^56, so a 64-bit lane holds the sum
 *   of up to 255 of them.  The a*b and m*n terms of a column go to separate
 *   accumulators (<= S terms each, S <= 255), split into low digit / high
 *   part before they are combined: 1 multiply + 1 add per product.
 * D32: 32-bit digits.  A product is < 2^64: the accumulator sums the full
 *   products mod 2^64 and a second one sums their high halves (p >> 32);
 *   the sum of the low halves is then (sum - (hisum << 32)) mod 2^64, exact
 *   since it is < 2^64.  1 multiply + 1 shift + 2 adds per product.
 *
 * m_k = (column low digit) * rho mod 2^D, computed with vpmuludq too.
 */

#include "a2m.h"
#include <immintrin.h>

#define LD(p, i) _mm256_load_si256((const __m256i *)((p) + (size_t)(i) * 4))
#define ST(p, i, v) _mm256_store_si256((__m256i *)((p) + (size_t)(i) * 4), (v))
#define BC(x) _mm256_set1_epi64x((long long)(x))
#define MUL(x, y) _mm256_mul_epu32((x), (y))

static A2M_NOINLINE void final_sub(uint64_t *C, const uint64_t *N, int tm,
    int S, const int D)
{
    const uint64_t mask = (D == 32) ? 0xffffffffULL : ((1ULL << D) - 1);
    int l, j;
    for (l = 0; l < 4; l++)
    {
        int64_t br = 0;
        if (!(tm & (1 << l)))
            continue;
        for (j = 0; j < S; j++)
        {
            int64_t t = (int64_t)C[j * 4 + l] - (int64_t)N[j * 4 + l] + br;
            C[j * 4 + l] = (uint64_t)t & mask;
            br = t >> D;
        }
    }
}

// accumulate one product into (l, h): D28: l += p; D32: l += p, h += p >> 32
#define ACC(l, h, x, y) do {                                                \
        __m256i _p = MUL((x), (y));                                         \
        (l) = _mm256_add_epi64((l), _p);                                    \
        if (D == 32) (h) = _mm256_add_epi64((h), _mm256_srli_epi64(_p, 32)); \
    } while (0)

// column sums -> (lo, hi): value = lo + hi * 2^D, lo < a few * 2^D
#define SPLIT(lo, hi, l, h) do {                                            \
        if (D == 32)                                                        \
        {                                                                   \
            __m256i _lo = _mm256_sub_epi64((l), _mm256_slli_epi64((h), 32)); /* sum of low halves */ \
            (lo) = _mm256_add_epi64((lo), _mm256_and_si256(_lo, dmask));    \
            (hi) = _mm256_add_epi64(_mm256_add_epi64((hi), (h)), _mm256_srli_epi64(_lo, 32)); \
        }                                                                   \
        else                                                                \
        {                                                                   \
            (lo) = _mm256_add_epi64((lo), _mm256_and_si256((l), dmask));    \
            (hi) = _mm256_add_epi64((hi), _mm256_srli_epi64((l), D));       \
        }                                                                   \
    } while (0)

A2M_INLINE void montmul_t(a2m_ctx_t *ctx, uint64_t *C, const uint64_t *A,
    const uint64_t *B, const int S, const int SQR, const int FULL, const int D)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i dmask = BC((D == 32) ? 0xffffffffULL : ((1ULL << D) - 1));
    const __m256i rho = _mm256_loadu_si256((const __m256i *)ctx->rho);
    const uint64_t *N = ctx->N;
    uint64_t *M = ctx->M;
    const __m256i n0 = LD(N, 0);
    __m256i carry = zero, mprev = zero;
    int k, tm;

#define COLUMN_BODY                                                         \
    {                                                                       \
        __m256i la0 = zero, la1 = zero, ha0 = zero, ha1 = zero;             \
        __m256i lm0 = zero, lm1 = zero, hm0 = zero, hm1 = zero;             \
        __m256i lo, hi;                                                     \
        const int ilo = (k - (S - 1) > 0) ? k - (S - 1) : 0;                \
        const int ihi = (k < S - 1) ? k : S - 1;                            \
        const int mhi = (k - 1 < S - 1) ? k - 1 : S - 1;                    \
        int i;                                                              \
        if (!SQR)                                                           \
        {                                                                   \
            if (FULL)                                                       \
            {                                                               \
                A2M_UNROLL_ALL                                              \
                for (i = ilo; i <= ihi; i++)                                \
                {                                                           \
                    if ((i - ilo) & 1) ACC(la1, ha1, LD(A, i), LD(B, k - i)); \
                    else ACC(la0, ha0, LD(A, i), LD(B, k - i));             \
                }                                                           \
            }                                                               \
            else                                                            \
            {                                                               \
                for (i = ilo; i + 1 <= ihi; i += 2)                         \
                {                                                           \
                    ACC(la0, ha0, LD(A, i), LD(B, k - i));                  \
                    ACC(la1, ha1, LD(A, i + 1), LD(B, k - i - 1));          \
                }                                                           \
                if (i <= ihi) ACC(la0, ha0, LD(A, i), LD(B, k - i));        \
            }                                                               \
        }                                                                   \
        else                                                                \
        {                                                                   \
            if (FULL)                                                       \
            {                                                               \
                A2M_UNROLL_ALL                                              \
                for (i = ilo; 2 * i < k; i++)                               \
                {                                                           \
                    if ((i - ilo) & 1) ACC(la1, ha1, LD(A, i), LD(A, k - i)); \
                    else ACC(la0, ha0, LD(A, i), LD(A, k - i));             \
                }                                                           \
            }                                                               \
            else                                                            \
            {                                                               \
                for (i = ilo; 2 * (i + 1) < k; i += 2)                      \
                {                                                           \
                    ACC(la0, ha0, LD(A, i), LD(A, k - i));                  \
                    ACC(la1, ha1, LD(A, i + 1), LD(A, k - i - 1));          \
                }                                                           \
                if (2 * i < k) ACC(la0, ha0, LD(A, i), LD(A, k - i));       \
            }                                                               \
            /* double the off-diagonal sums: < S products each, fits */     \
            la0 = _mm256_add_epi64(la0, la1); ha0 = _mm256_add_epi64(ha0, ha1); \
            la1 = ha1 = zero;                                               \
            if (D == 32)                                                    \
            {                                                               \
                /* doubling mod 2^64 loses the top bit: split first */      \
                __m256i _lo = _mm256_sub_epi64(la0, _mm256_slli_epi64(ha0, 32)); \
                la0 = _mm256_add_epi64(_lo, _lo);                           \
                ha0 = _mm256_add_epi64(ha0, ha0);                           \
                la0 = _mm256_add_epi64(la0, _mm256_slli_epi64(ha0, 32));    \
            }                                                               \
            else                                                            \
                la0 = _mm256_add_epi64(la0, la0);                           \
            if (((k & 1) == 0) && (k / 2 <= S - 1))                         \
                ACC(la1, ha1, LD(A, k / 2), LD(A, k / 2));                  \
        }                                                                   \
        if (FULL)                                                           \
        {                                                                   \
            A2M_UNROLL_ALL                                                  \
            for (i = ilo; i <= mhi; i++)                                    \
            {                                                               \
                __m256i x = (i == k - 1) ? mprev : LD(M, i);                \
                if ((i - ilo) & 1) ACC(lm1, hm1, x, LD(N, k - i));          \
                else ACC(lm0, hm0, x, LD(N, k - i));                        \
            }                                                               \
        }                                                                   \
        else                                                                \
        {                                                                   \
            int iend = (mhi == k - 1) ? mhi - 1 : mhi;                      \
            for (i = ilo; i + 1 <= iend; i += 2)                            \
            {                                                               \
                ACC(lm0, hm0, LD(M, i), LD(N, k - i));                      \
                ACC(lm1, hm1, LD(M, i + 1), LD(N, k - i - 1));              \
            }                                                               \
            if (i <= iend) ACC(lm0, hm0, LD(M, i), LD(N, k - i));           \
            if (mhi == k - 1 && mhi >= ilo) ACC(lm1, hm1, mprev, LD(N, 1));  \
        }                                                                   \
        /* combine: carry + the four sums, as lo + hi * 2^D */             \
        lo = _mm256_and_si256(carry, dmask);                                \
        hi = _mm256_srli_epi64(carry, D);                                   \
        SPLIT(lo, hi, la0, ha0);                                            \
        SPLIT(lo, hi, la1, ha1);                                            \
        SPLIT(lo, hi, lm0, hm0);                                            \
        SPLIT(lo, hi, lm1, hm1);                                            \
        if (k < S)                                                          \
        {                                                                   \
            __m256i mk = _mm256_and_si256(MUL(lo, rho), dmask);             \
            __m256i p = MUL(mk, n0);                                        \
            ST(M, k, mk);                                                   \
            mprev = mk;                                                     \
            /* lo + p == 0 mod 2^D */                                       \
            lo = _mm256_add_epi64(lo, _mm256_and_si256(p, dmask));          \
            hi = _mm256_add_epi64(hi, _mm256_srli_epi64(p, D));             \
            carry = _mm256_add_epi64(hi, _mm256_srli_epi64(lo, D));         \
        }                                                                   \
        else                                                                \
        {                                                                   \
            ST(C, k - S, _mm256_and_si256(lo, dmask));                      \
            carry = _mm256_add_epi64(hi, _mm256_srli_epi64(lo, D));         \
        }                                                                   \
    }

    if (FULL)
    {
        A2M_UNROLL_ALL
        for (k = 0; k < 2 * S; k++)
            COLUMN_BODY
    }
    else
    {
        for (k = 0; k < S; k++)
            COLUMN_BODY
        for (k = S; k < 2 * S; k++)
            COLUMN_BODY
    }
#undef COLUMN_BODY

    tm = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(carry, zero)));
    if (tm)
        final_sub(C, N, tm, S, D);
}

#define INST(S, D)                                                          \
static void mul_d##D##_s##S(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a, const uint64_t *b) \
{ montmul_t(ctx, c, a, b, S, 0, (S) <= A2M_FULL_MAX, D); }                 \
static void sqr_d##D##_s##S(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a) \
{ montmul_t(ctx, c, a, a, S, 1, (S) <= A2M_FULL_MAX, D); }

#define INSTD(D)                                                            \
INST(2, D) INST(3, D) INST(4, D) INST(5, D) INST(6, D) INST(7, D) INST(8, D) \
INST(9, D) INST(10, D) INST(11, D) INST(12, D) INST(13, D) INST(14, D)     \
INST(15, D) INST(16, D) INST(17, D) INST(18, D) INST(19, D) INST(20, D)    \
INST(21, D) INST(22, D) INST(23, D) INST(24, D) INST(25, D) INST(26, D)    \
INST(27, D) INST(28, D) INST(29, D) INST(30, D) INST(31, D) INST(32, D)    \
static void mul_d##D##_sN(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a, const uint64_t *b) \
{ montmul_t(ctx, c, a, b, ctx->S, 0, 0, D); }                              \
static void sqr_d##D##_sN(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a)  \
{ montmul_t(ctx, c, a, a, ctx->S, 1, 0, D); }

INSTD(28)
INSTD(32)

#define P(S, D) { mul_d##D##_s##S, sqr_d##D##_s##S }
#define TAB(D) { { 0, 0 }, { 0, 0 }, P(2, D), P(3, D), P(4, D), P(5, D),     \
    P(6, D), P(7, D), P(8, D), P(9, D), P(10, D), P(11, D), P(12, D),         \
    P(13, D), P(14, D), P(15, D), P(16, D), P(17, D), P(18, D), P(19, D),     \
    P(20, D), P(21, D), P(22, D), P(23, D), P(24, D), P(25, D), P(26, D),     \
    P(27, D), P(28, D), P(29, D), P(30, D), P(31, D), P(32, D) }

typedef struct { a2m_mul_fn mul; a2m_sqr_fn sqr; } pair_t;
static const pair_t tab28[33] = TAB(28);
static const pair_t tab32[33] = TAB(32);

int a2m_get_d28(int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr)
{
    if (S < 2 || S > 255)
        return 0;
    if (S <= 32) { *mul = tab28[S].mul; *sqr = tab28[S].sqr; }
    else { *mul = mul_d28_sN; *sqr = sqr_d28_sN; }
    return 1;
}

int a2m_get_d32(int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr)
{
    if (S < 2 || S > (1 << 28))
        return 0;
    if (S <= 32) { *mul = tab32[S].mul; *sqr = tab32[S].sqr; }
    else { *mul = mul_d32_sN; *sqr = sqr_d32_sN; }
    return 1;
}

int a2m_get(int variant, int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr)
{
    switch (variant)
    {
    case A2M_FP52: return a2m_get_fp52(S, mul, sqr);
    case A2M_D28: return a2m_get_d28(S, mul, sqr);
    case A2M_D32: return a2m_get_d32(S, mul, sqr);
    }
    return 0;
}
