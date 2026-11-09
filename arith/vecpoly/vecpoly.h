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
 * vecpoly: polynomials over Z/NZ, 8 independent polynomials per vector (one
 * per lane, all modulo the same N), for AVX-ECM's FFT stage 2.  Design:
 * arith/avxecm_stage2_fft_design.md.
 *
 * Coefficients: radix 2^52, n words, Montgomery form wrt R = 2^(52n), the
 * vec_bignum_t layout (word w of lane l at c[w*8 + l]); a polynomial is an
 * array of coefficients with stride 8n u64.  Invariant: coefficients are in
 * [0, 2N).  Monic polynomials are stored without their leading 1.
 *
 * Products of long polynomials go through NTTs modulo K primes p_i in
 * (2^49, 2^50) and the CRT; NTT-form layout: prime i, index t, lane l at
 * v[(i*L + t)*8 + l].
 */

#ifndef VECPOLY_H
#define VECPOLY_H

#include <stdint.h>
#include <stddef.h>
#include <gmp.h>
#include "../rns/rns.h"
#include "../vecmont/vecmont.h"

#if defined(RNS_HAVE_AVX512) || defined(RNS_HAVE_AVX2)
#define VP_HAVE_KERNELS 1
#endif

#define VP_LANES        8
#define VP_PRIME_BITS   50          // primes in (2^49, 2^50)
#define VP_MAX_LOGL     20          // p = 1 mod 2^20
#define VP_MASK52       0xfffffffffffffULL
#define VP_GRP          4           // primes per group for the second operand
#define VP_SCHOOL_CAP   16          // largest schoolbook degree (w->sm sizing)

// shared, read-only after vp_ctx_init (one per N and maximum length)
typedef struct
{
    int n;                  // words per coefficient
    int K;                  // NTT primes
    int logLmax;
    uint32_t Lmax;          // longest transform supported
    mpz_t N;
    uint64_t n0;            // -N^-1 mod 2^52
    uint64_t *Nw;           // N, n words
    uint64_t *N2w;          // 2N, n words
    uint64_t *p;            // [K]
    uint64_t *pneg;         // [K] -p^-1 mod 2^52
    double *pinvd;          // [K] 1/p
    uint64_t *cvt;          // [K][n] 2^(52k + 104) mod p_i
    uint64_t *crt;          // [K+1][n] C_i = (P/p_i) R^-1 2^104 mod N;
                            //          C_K = -P R^-1 2^104 mod N
    uint64_t *yscale;       // [logLmax+1][K][2] Shoup pair of
                            //   L^-1 2^52 (P/p_i)^-1 mod p_i, L = 2^l
    uint64_t *tw;           // [K][Lmax/2][2] forward twiddles w_Lmax^j and
                            //   Shoup w'; level m reads every (Lmax/2m)-th
    uint64_t *itw;          // [K][Lmax/2][2] inverse (w_Lmax^-j)
    int school_max;         // degrees <= this use schoolbook products
} vp_ctx_t;

// per thread
typedef struct
{
    vp_ctx_t *ctx;
    vec_monty_t md;         // for the template mul: vrho, NWORDS, fpwork
    vecmont_mul_fn mul;
    vecmont_sqr_fn sqr;
    vec_bignum_t vn;        // N broadcast to all lanes
    vec_bignum_t vs;        // template scratch
    uint64_t *ntt;          // NTT buffer, K * Lmax vectors
    uint64_t *grp;          // second operand, VP_GRP primes at a time:
                            //   VP_GRP * Lmax vectors
    uint64_t *tmp[2];       // coefficient scratch: Lmax/2 and Lmax/4
    uint64_t *sm;           // schoolbook scratch, 4 * VP_SCHOOL_CAP
    uint64_t *one;          // Montgomery 1 (R mod N), one coefficient
} vp_work_t;

// coefficient j of polynomial a
#define VP_COEF(ctx, a, j) ((a) + (size_t)(j) * 8 * (size_t)(ctx)->n)
// NTT residue vector, prime i, index t, length L
#define VP_RES(v, L, i, t) ((v) + ((size_t)(i) * (L) + (t)) * 8)

// ---- setup (vp_setup.c) ----
// N odd, n words with N < 2^(52n - 3); transforms up to 2^logLmax.
vp_ctx_t *vp_ctx_init(mpz_t N, int n, int logLmax);
void vp_ctx_free(vp_ctx_t *ctx);
vp_work_t *vp_work_init(vp_ctx_t *ctx);
void vp_work_free(vp_work_t *w);
// NTT primes needed for n words and transforms up to 2^logLmax
int vp_nprimes(int n, int logLmax);
// bytes owned by a context (shared by all threads) / by one work area
size_t vp_ctx_bytes(int n, int K, int logLmax);
size_t vp_work_bytes(int n, int K, int logLmax);
uint64_t *vp_alloc(size_t ncoef, int n);    // ncoef coefficients, zeroed
void vp_free(uint64_t *a);
// lane values <-> integers (no Montgomery conversion)
void vp_get(vp_ctx_t *ctx, mpz_t x, const uint64_t *c, int lane);
void vp_set(vp_ctx_t *ctx, uint64_t *c, mpz_t x, int lane);    // x < 2^(52n)

// ---- kernels (vp_kern.c) ----
// The _r forms work on primes i0..i1-1 only, with v holding just those
// (prime i at v + (i - i0) L 8); the plain forms on all K primes.
// residues of coefficients a[0..len-1] (< 2^(52n)) into v, zero-padded to L
void vpk_to_ntt_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L,
    const uint64_t *a, uint32_t len, int i0, int i1);
void vpk_fwd_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L, int i0, int i1);
void vpk_inv_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L, int i0, int i1);
// c = a * b * 2^-52, pointwise (c may alias a or b)
void vpk_pmul_r(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t L, int i0, int i1);
#define vpk_to_ntt(ctx, v, L, a, len) vpk_to_ntt_r((ctx), (v), (L), (a), (len), 0, (ctx)->K)
#define vpk_fwd(ctx, v, L) vpk_fwd_r((ctx), (v), (L), 0, (ctx)->K)
#define vpk_inv(ctx, v, L) vpk_inv_r((ctx), (v), (L), 0, (ctx)->K)
#define vpk_pmul(ctx, c, a, b, L) vpk_pmul_r((ctx), (c), (a), (b), (L), 0, (ctx)->K)
// CRT of indices start..start+count-1 of an inverse transform of length L
// (all values Montgomery products) into coefficients a[0..count-1]
void vpk_from_ntt(const vp_ctx_t *ctx, uint64_t *a, const uint64_t *v,
    uint32_t L, uint32_t start, uint32_t count);
// coefficient arithmetic mod 2N on [0, 2N): c = a + b, a - b, -a
void vpk_add(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t cnt);
void vpk_sub(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t cnt);
void vpk_neg(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a, uint32_t cnt);

// ---- polynomials (vp_poly.c) ----
// c[0..cnt) = a[0..cnt) * b[0..cnt), coefficient-wise (template mul)
void vp_cmul(vp_work_t *w, uint64_t *c, const uint64_t *a, const uint64_t *b,
    uint32_t cnt);
// full product c[0..la+lb-1) = a * b  (c must not alias a or b)
void vp_mul(vp_work_t *w, uint64_t *c, const uint64_t *a, uint32_t la,
    const uint64_t *b, uint32_t lb);
// monic product: a, b monic of degree d (d coefficients each); c (2d
// coefficients) = a*b, monic of degree 2d.  c may be a when b = a + d
// (in place); otherwise c must not overlap a or b.  d <= Lmax/4.
void vp_mul_monic(vp_work_t *w, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t d);
// product tree of the m roots (m a power of 2).  tree == NULL: only the
// root polynomial, written to out (m coefficients, monic; roots may be
// out).  Otherwise tree[s] (m coefficients each, s = 0..log2(m)-1)
// receives level s (nodes of degree 2^s), the root goes to out, and roots
// may be tree[0].  tree[s] may be NULL for any s >= 1 (level not kept; out
// is the scratch for it); vp_eval_tree then rebuilds it from the nearest
// kept level below.
void vp_from_roots(vp_work_t *w, uint64_t *out, uint64_t **tree,
    const uint64_t *roots, uint32_t m);
// inv[0..m) = floor(X^(2m-1) / F), F monic of degree m (f: m coefficients);
// scratch: m coefficients (not a w->tmp buffer)
void vp_recip(vp_work_t *w, uint64_t *inv, const uint64_t *f, uint32_t m,
    uint64_t *scratch);

// reduction modulo a fixed monic F of degree m (2m <= Lmax)
#define VP_MODF_KEEP_INV  1     // keep the length-2m transform of inv
#define VP_MODF_KEEP_F    2     // keep the length-m transform of f
typedef struct
{
    uint32_t m;
    uint64_t *f;            // F without its leading 1 (m coefficients, not owned)
    uint64_t *inv;          // floor(X^(2m-1)/F) (m coefficients, owned)
    uint64_t *ntt_inv;      // length-2m transform of inv, or NULL (recomputed)
    uint64_t *ntt_f;        // length-m transform of f, or NULL (recomputed)
} vp_modf_t;

// scratch: m coefficients for the reciprocal, or NULL (allocated briefly).
// Without VP_MODF_KEEP_F, f must stay valid while vp_mulmod_monic is used.
vp_modf_t *vp_modf_init(vp_work_t *w, uint64_t *f, uint32_t m, int flags,
    uint64_t *scratch);
void vp_modf_free(vp_modf_t *mf);
// bytes a vp_modf_t for (n, K, m, flags) owns
size_t vp_modf_bytes(int n, int K, uint32_t m, int flags);
// h = h * g mod F, h of degree < m, g monic of degree m (uses w->tmp[0])
void vp_mulmod_monic(vp_work_t *w, vp_modf_t *mf, uint64_t *h, const uint64_t *g);
// h = g mod F for g monic of degree m (g - F); h may be mf->f
void vp_mod_monic(vp_work_t *w, vp_modf_t *mf, uint64_t *h, const uint64_t *g);
// vals[k] = H(root_k), the roots of F, given the product tree of F
// (tree[s], from vp_from_roots).  h of degree < m.  vals may alias h.
// scratch: m coefficients, needed only when the tree has NULL levels.
void vp_eval_tree(vp_work_t *w, vp_modf_t *mf, uint64_t **tree,
    uint64_t *vals, const uint64_t *h, uint64_t *scratch);
// prod = product of a[0..cnt) (one coefficient)
void vp_prod(vp_work_t *w, uint64_t *prod, const uint64_t *a, uint32_t cnt);

#endif
