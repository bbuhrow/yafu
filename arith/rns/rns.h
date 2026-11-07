/*
 * RNS Montgomery modular multiplication (small regime: 52-bit channels).
 *
 * Bajard-style two-base RNS Montgomery with Gandino/Yao precomputation:
 *   base B = {b_1..b_nB}, base C = {c_1..c_nC}, all b_i, c_j = 2^52 - d (d odd,
 *   small), pairwise coprime, plus a redundant channel m_r = 2^52 used for an
 *   exact (Shenoy-Kumaresan) second base extension.  The first extension
 *   B -> C is inexact by alpha*M_B (alpha <= nB); the bounds absorb it.
 *
 * Representation of an integer X (0 <= X < M_B):
 *   B channels:  |X|_{b_i}
 *   C channels:  |X * C_j^{-1}|_{c_j}      (pre-scaled, as in Yao Alg. 3)
 *   r channel:   X mod 2^52
 * Residues are kept partially reduced: any value < 2^52 congruent to the
 * true residue.  Montgomery form of a is X = a*M_B mod N (any representative
 * X < k*N, see below).
 *
 * Bounds: mulmod outputs are < k*N with k = nB + 2.  mulmod inputs may be
 * as large as h*k*N, where h = RNS_HEADROOM.  submod(x, y) represents
 * X - Y + RNS_SUBMULT*k*N, so Y must be < RNS_SUBMULT*k*N.  Callers are
 * responsible for keeping chains of add/sub within h*k*N between mulmods.
 *
 * Setup is split: rns_fixed_t depends only on the size class (max bits of
 * N); rns_ctx_t holds the N-dependent constants.
 */

#ifndef RNS_H
#define RNS_H

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>

// The AVX-512 vector kernels (*_avx512) are compiled for any AVX-512F
// target: the compiler's own macros (gcc/clang -march=...) or yafu's build
// defines (USE_AVX512F [+ IFMA]; MSVC has no IFMA macro).  With IFMA
// (RNS_HAVE_IFMA) the 52-bit products are madd52lo/hi; without it they are
// emulated exactly with FP64 FMAs (needs AVX512DQ for the conversions, which
// every AVX-512 CPU with AVX512BW has).  Otherwise only the scalar
// references and the setup code are built.
#if (defined(__AVX512F__) && defined(__AVX512IFMA__)) || \
    (defined(USE_AVX512F) && defined(IFMA))
#define RNS_HAVE_AVX512 1
#define RNS_HAVE_IFMA 1
#elif (defined(__AVX512F__) && defined(__AVX512DQ__)) || defined(USE_AVX512F)
#define RNS_HAVE_AVX512 1
#endif

// compiler shims for the kernels
#if defined(_MSC_VER) && !defined(__clang__)
#define RNS_INLINE static __forceinline
#define RNS_UNROLL(n)
#define RNS_LIKELY(x) (x)
#else
#define RNS_INLINE static inline __attribute__((always_inline))
#define RNS_PRAGMA_(x) _Pragma(#x)
#if defined(__clang__)
#define RNS_UNROLL(n) RNS_PRAGMA_(unroll n)
#else
#define RNS_UNROLL(n) RNS_PRAGMA_(GCC unroll n)
#endif
#define RNS_LIKELY(x) __builtin_expect(!!(x), 1)
#endif

#define RNS_WBITS       52
#define RNS_MASK52      0xfffffffffffffULL
#define RNS_HEADROOM    8       // max mulmod input = HEADROOM * k * N
#define RNS_SUBMULT     2       // submod adds SUBMULT * k * N

typedef struct
{
    int maxbits;            // N < 2^maxbits
    int nB, nC;             // channels per base
    int nch;                // total residues per number: nB + nC + 1
    int k;                  // output bound multiple: outputs < k*N

    uint64_t *b, *db;       // b_i = 2^52 - db_i
    uint64_t *c, *dc;       // c_j = 2^52 - dc_j

    // first extension B -> C u {r}: BtoC[j*nB + i] = |B_i|_{c_j},
    // row j = nC is |B_i|_r (mod 2^52)
    uint64_t *BtoC;

    // second extension C -> B u {r}: CtoB[i*nC + j] = |C_j|_{b_i},
    // row i = nB is |C_j|_r
    uint64_t *CtoB;
    uint64_t *negMC_b;      // |-M_C|_{b_i}
    uint64_t MCinv_r;       // |M_C^{-1}|_r

    uint64_t *CjMBinv_c;    // |C_j * M_B^{-1}|_{c_j}    (multiplies X'Y')
    uint64_t MBinv_r;       // |M_B^{-1}|_r

    mpz_t MB, MC;
    mpz_t *Bi;              // B_i = M_B / b_i, for CRT out
    uint64_t *Biinv_b;      // |B_i^{-1}|_{b_i}
    uint64_t *Cjinv_c;      // |C_j^{-1}|_{c_j}, for conversion in
} rns_fixed_t;

typedef struct
{
    rns_fixed_t *fx;
    mpz_t N;
    uint64_t *f1_b;         // |-N^{-1} * B_i^{-1}|_{b_i}
    uint64_t *g2_c;         // |N * M_B^{-1} * C_j^{-1}|_{c_j}
    uint64_t NMBinv_r;      // |N * M_B^{-1}|_r
    uint64_t *BtoCN;        // BtoC rows folded with g2_c / NMBinv_r (Gandino)
    uint64_t *subconst;     // representation of RNS_SUBMULT*k*N
    uint64_t *one;          // Montgomery form of 1 (M_B mod N)
    uint64_t *vscratch;     // nB * RNS_VL words, 64-byte aligned
    mpz_t MBmodN, MBinvmodN;
} rns_ctx_t;

// Vector layout: RNS_VL independent numbers modulo the same N (e.g. ECM
// curves), channel-major: v[ch * RNS_VL + lane], 64-byte aligned.  The
// vector routines use ctx->vscratch, so one ctx per thread.
#define RNS_VL 8

// fixed (size-class) setup
rns_fixed_t *rns_fixed_init(int maxbits);
void rns_fixed_free(rns_fixed_t *fx);

// N-dependent setup.  N must be odd, < 2^maxbits, coprime to M_B*M_C
// (the moduli are not prime); returns NULL otherwise.
rns_ctx_t *rns_ctx_init(rns_fixed_t *fx, mpz_t N);
void rns_ctx_free(rns_ctx_t *ctx);

// residue vectors are uint64_t[fx->nch]
uint64_t *rns_alloc(rns_ctx_t *ctx);

// raw conversion of an integer 0 <= X < M_B (no Montgomery scaling)
void rns_from_int(rns_ctx_t *ctx, uint64_t *x, mpz_t X);
void rns_to_int(rns_ctx_t *ctx, mpz_t X, uint64_t *x);

// a (any integer) -> Montgomery form; and back to a mod N in [0, N)
void rns_to_monty(rns_ctx_t *ctx, uint64_t *x, mpz_t a);
void rns_from_monty(rns_ctx_t *ctx, mpz_t a, uint64_t *x);

// scalar reference arithmetic
void rns_mulmod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rns_sqrmod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x);
void rns_addmod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rns_submod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);

// 64-byte aligned, zeroed; free with rns_vfree
uint64_t *rns_aligned_alloc(size_t words);

// vector residues: uint64_t[fx->nch * RNS_VL], free with rns_vfree
uint64_t *rns_valloc(rns_ctx_t *ctx);
void rns_vfree(uint64_t *v);
void rns_vset_lane(rns_ctx_t *ctx, uint64_t *v, int lane, uint64_t *x);
void rns_vget_lane(rns_ctx_t *ctx, uint64_t *x, uint64_t *v, int lane);

#ifdef RNS_HAVE_AVX512
// AVX-512 kernels (rns_avx512.c; IFMA or FP-emulated); results are bit-identical to the
// scalar reference, lane by lane.  z may alias x or y.
void rns_vmulmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rns_vsqrmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x);
void rns_vaddmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rns_vsubmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);

// d = b^e in Montgomery form, all lanes share the exponent e (> 0).
// Sliding window, the same window rule as AVX-ECM's vecslidingmodexp.
void rns_vpowm_avx512(rns_ctx_t *ctx, uint64_t *d, uint64_t *b, mpz_t e);
#endif

#endif
