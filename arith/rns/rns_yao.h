/*
 * RNS Montgomery modular multiplication, Yao regime: n + n multi-word
 * channels (Yao, Fan, Cheung, Verbauwhede, IEEE TC 2014) plus a redundant
 * channel r = 2^52.  See arith/rns_design.md §0b/§0c.
 *
 * Channels: b_i = 2^w - db_i, c_j = 2^w - dc_j, w = 52*s, db < dc (all of
 * B's offsets below C's) so the extension constants have uniform sign:
 *   |B_i|_{c_j} == prod_{k != i} (dc_j - db_k)                  (> 0)
 *   |C_j|_{b_i} == (-1)^(n-1) prod_{k != j} (dc_k - db_i)
 * Every difference is even, so a product of n-1 of them has >= n-1 trailing
 * zeros (Yao §3).  When a constant exceeds 52 bits, the common power of two
 * over each output channel is factored out: Bt[i][j] = |B_i|~ >> zB[j] and
 * Ct[j][i] = |C_j|~ >> zC[i], and the summed extension is shifted back.
 *
 * Two vector kernels share this setup:
 *  - rnsy_v*_avx512  (rns_yao_avx512.c): one number per call, n = 8, the 8
 *    channels of each base in the 8 lanes.  Lowest latency for one curve.
 *  - rnsyb_v*_avx512 (rns_yaob_avx512.c): 8 independent numbers mod the
 *    same N per call, one per lane, channels as the loop index, any n the
 *    constants allow.  Highest throughput for 8 curves.
 *
 * Single-number representation of X (0 <= X < M_B), radix 2^52:
 *   x[(0*s + t)*n + i]  word t of |X|_{b_i}
 *   x[(1*s + t)*n + j]  word t of |X * C_j^{-1}|_{c_j}   (pre-scaled)
 *   x[2*n*s]            X mod 2^52
 * (for n = 8 this is also the rnsy_v* vector layout).  Batch layout:
 *   v[((base*n + ch)*s + t)*8 + lane],  r at v[2*n*s*8 + lane].
 * Channel residues are lazily reduced (any value < 2^w in the class).
 * Bounds as in rns.h: outputs < k*N, k = n + 2; mulmod inputs < h*k*N.
 */

#ifndef RNS_YAO_H
#define RNS_YAO_H

#include "rns.h"

#define RNSY_NCH 8          // SIMD lanes (and n for the one-number kernel)
#define RNSY_MAXN 16        // max channels per base
#define RNSY_MAXS 10        // max words per channel (vector kernels)

typedef struct
{
    int maxbits;
    int n;                  // channels per base
    int s, w;               // words per channel, bits per channel
    int k;                  // output bound multiple, n + 2
    int nwords;             // words per number: 2*n*s + 1
    int rword;              // index of the r word: 2*n*s
    uint64_t db[RNSY_MAXN], dc[RNSY_MAXN];
    uint64_t Bt[RNSY_MAXN][RNSY_MAXN];  // Bt[i][j], first extension, >> zB[j]
    uint64_t Ct[RNSY_MAXN][RNSY_MAXN];  // Ct[j][i], second extension (negated), >> zC[i]
    int zB[RNSY_MAXN], zC[RNSY_MAXN];   // per-output shifts (0 if all fit)
    int shifted;                        // any nonzero zB/zC
    int cneg;                           // n even: |C_j|_{b_i} == -Ct, so U = K - S + ...
    int vB, vC;                         // max bits of the unshifted constants
    int vBt, vCt;                       // max bits of the stored constants
    uint64_t Br[RNSY_MAXN];             // |B_i|_r
    uint64_t Cr[RNSY_MAXN];             // |C_j|_r
    uint64_t MBinv_r, MCinv_r;
    mpz_t b[RNSY_MAXN], c[RNSY_MAXN];
    mpz_t negMC_b[RNSY_MAXN];           // |-M_C|_{b_i}
    mpz_t K_b[RNSY_MAXN];               // multiple of b_i > max S_i
    mpz_t g1_c[RNSY_MAXN];              // |C_j M_B^{-1}|_{c_j}
    mpz_t Cjinv_c[RNSY_MAXN];           // |C_j^{-1}|_{c_j}
    mpz_t Bi[RNSY_MAXN], Biinv_b[RNSY_MAXN];
    mpz_t MB, MC, two_w;
} rnsy_fixed_t;

typedef struct
{
    rnsy_fixed_t *fx;
    mpz_t N;
    mpz_t f1_b[RNSY_MAXN];              // |-N^{-1} B_i^{-1}|_{b_i}
    mpz_t g2_c[RNSY_MAXN];              // |N M_B^{-1} C_j^{-1}|_{c_j}
    uint64_t BrN[RNSY_MAXN];            // |B_i|_r * N M_B^{-1} mod 2^52
    uint64_t *subconst;                 // representation of SUBMULT*k*N
    uint64_t *one;                      // Montgomery 1
    mpz_t MBmodN, MBinvmodN;

    // rotated constants F_i = f * 2^(52 i) mod m (i < s), so a product by a
    // per-N constant is sum_i a_i F_i: s+1 columns, no high half to fold.
    mpz_t F1[RNSY_MAXN][RNSY_MAXS], G1[RNSY_MAXN][RNSY_MAXS], G2[RNSY_MAXN][RNSY_MAXS];

    // one-number kernel (n = 8, unshifted constants only): digit-major x 8
    // lanes, lane = channel, in one aligned block.  NULL when unavailable.
    uint64_t *vblock;
    uint64_t *vdb, *vdc;                // 1 vector each
    uint64_t *vf1, *vg1, *vg2, *vnegMC; // s vectors each
    uint64_t *vK, *vEb, *vEc;           // s+1 vectors each
    uint64_t *vBrN, *vCr;               // 1 vector each (lane = channel)
    uint64_t *vbt;                      // 8 vectors: vbt[i] lane j = Bt[i][j]
    uint64_t *vct;                      // 8 vectors: vct[j] lane i = Ct[j][i]
    uint64_t *vscratch;                 // 2*s vectors
    uint64_t *vF1, *vG1, *vG2;          // s*s vectors each, (i*s + digit)*8 + lane

    // batch kernel: per-channel scalar constants (broadcast in the kernel),
    // the alpha' table and scratch, in one aligned block
    uint64_t *bblock;
    uint64_t *bF1, *bG1, *bG2;          // [ch][rotation i][digit], n*s*s each
    uint64_t *bK;                       // [ch][digit], s+2 digits
    uint64_t *bEb, *bEc;                // [ch][digit], s+1 digits
    uint64_t *balpha;                   // [i][digit][16]: alpha * |-M_C|_{b_i}, alpha < 16
    uint64_t *bscratch;                 // n*s vectors
} rnsy_ctx_t;

// fixed setup for n channels per base; NULL if the extension constants
// can't be made to fit 52 bits (even after the trailing-zero shift)
rnsy_fixed_t *rnsy_fixed_init_n(int maxbits, int n);
// n = 8 (required by the one-number kernel)
rnsy_fixed_t *rnsy_fixed_init(int maxbits);
// cost-model choice of n for the batch kernel
int rnsyb_choose_n(int maxbits);
void rnsy_fixed_free(rnsy_fixed_t *fx);
rnsy_ctx_t *rnsy_ctx_init(rnsy_fixed_t *fx, mpz_t N);  // NULL if gcd(N, M_B M_C) != 1
void rnsy_ctx_free(rnsy_ctx_t *ctx);

uint64_t *rnsy_alloc(rnsy_ctx_t *ctx);          // 64-byte aligned, free with rns_vfree
void rnsy_from_int(rnsy_ctx_t *ctx, uint64_t *x, mpz_t X);
void rnsy_to_int(rnsy_ctx_t *ctx, mpz_t X, uint64_t *x);
void rnsy_to_monty(rnsy_ctx_t *ctx, uint64_t *x, mpz_t a);
void rnsy_from_monty(rnsy_ctx_t *ctx, mpz_t a, uint64_t *x);

// scalar reference (GMP per channel); with RNSY_COLUMN_FOLD = 0 it defines
// the exact residues both vector kernels reproduce
void rnsy_mulmod_ref(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rnsy_addmod_ref(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rnsy_submod_ref(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);

#ifdef RNS_HAVE_AVX512
// One-number kernels (rns_yao_avx512.c): n = 8, z may alias x or y,
// specialized for s = 1..10.
void rnsy_vmulmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rnsy_vsqrmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x);
void rnsy_vaddmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rnsy_vsubmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
// d = b^e in Montgomery form (sliding window, as rns_vpowm_avx512)
void rnsy_vpowm_avx512(rnsy_ctx_t *ctx, uint64_t *d, uint64_t *b, mpz_t e);
// 1 if the vector kernels were built bit-exact with the reference
// (RNSY_COLUMN_FOLD = 0), else 0: results are then valid representations
// (congruent, consistent across channels, same bounds) but other values
int rnsy_vec_exact(void);

// Batch kernels (rns_yaob_avx512.c): 8 numbers mod the same N per call,
// any n, s = 1..10; z may alias x or y.  Lane l of a batch vector holds
// number l; rnsyb_vset_lane / vget_lane convert from / to the
// single-number representation.
uint64_t *rnsyb_alloc(rnsy_ctx_t *ctx);         // (2*n*s + 1) * 8 words
void rnsyb_vset_lane(rnsy_ctx_t *ctx, uint64_t *v, int lane, uint64_t *x);
void rnsyb_vget_lane(rnsy_ctx_t *ctx, uint64_t *x, uint64_t *v, int lane);
void rnsyb_vmulmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rnsyb_vsqrmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x);
void rnsyb_vaddmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
void rnsyb_vsubmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y);
// all 8 lanes share the exponent e
void rnsyb_vpowm_avx512(rnsy_ctx_t *ctx, uint64_t *d, uint64_t *b, mpz_t e);
#endif

#endif
