/*
 * AVX2 vector Montgomery multiplication experiments (handoff 4.3).
 *
 * 4 numbers per vector (one per 64-bit lane), all modulo the same N, one
 * digit per lane-word: a number of S digits is uint64_t[S * 4], digit j of
 * lane l at [j * 4 + l], 32-byte aligned.  Three digit formats:
 *
 *   FP52  52-bit digits, products via FP64 FMA (two FMAs give the exact hi
 *         and lo halves; MXCSR round-toward-zero is set for the call)
 *   D28   28-bit digits, vpmuludq 32x32 -> 64: 56-bit products summed whole
 *   D32   32-bit digits, vpmuludq: 64-bit products, high halves split off
 *
 * Every kernel computes the same almost-Montgomery product as AVX-ECM's
 * cores: c = (a*b + m*N) / R, R = 2^(dbits * S), m = -a*b/N mod R, minus N
 * if that is >= R, so c < R for inputs a, b < R.  c may alias a or b.
 */

#ifndef A2M_H
#define A2M_H

#include <stdint.h>
#include <stddef.h>

#define A2M_VL 4

enum { A2M_FP52 = 0, A2M_D28 = 1, A2M_D32 = 2, A2M_NVARIANT = 3 };

typedef struct
{
    int variant;
    int dbits;          // 52, 28, 32
    int S;              // digits per number
    uint64_t *N;        // S * 4 words, normalized digits (all lanes equal)
    uint64_t rho[4];    // -1/N mod 2^dbits, per lane
    uint64_t *M;        // scratch, S * 4 words
    uint64_t *W;        // FP52 scratch: 3 * S * 4 words
} a2m_ctx_t;

typedef void (*a2m_mul_fn)(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a, const uint64_t *b);
typedef void (*a2m_sqr_fn)(a2m_ctx_t *ctx, uint64_t *c, const uint64_t *a);

// 0 if (variant, S) is not supported
int a2m_get(int variant, int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr);

int a2m_get_fp52(int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr);
int a2m_get_d28(int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr);
int a2m_get_d32(int S, a2m_mul_fn *mul, a2m_sqr_fn *sqr);

#if defined(_MSC_VER) && !defined(__clang__)
#define A2M_INLINE static __forceinline
#define A2M_NOINLINE __declspec(noinline)
#define A2M_UNROLL_ALL
#else
#define A2M_INLINE static inline __attribute__((always_inline))
#define A2M_NOINLINE __attribute__((noinline))
#if defined(__clang__)
#define A2M_UNROLL_ALL _Pragma("unroll")
#else
#define A2M_UNROLL_ALL _Pragma("GCC unroll 256")
#endif
#endif

// fully unroll the column loop for S <= this.  Much lower than the
// AVX-512 template's 24: with 16 ymm registers and the code size of a full
// unroll, the loop form is faster from S ~ 13 (measured: 12 is best for
// D28 at S = 10, 8-12 for FP52; 24 cost FP52 25% at S = 20)
#ifndef A2M_FULL_MAX
#define A2M_FULL_MAX 12
#endif

#endif
