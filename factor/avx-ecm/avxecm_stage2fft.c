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
 * AVX-ECM stage 2, FFT continuation (arith/avxecm_stage2_fft_design.md):
 * baby roots x([k]Q), k in S1 = {1 <= k < D/2, gcd(k, D) = 1}, give F;
 * giant roots x([jD]Q) in blocks of dF give G_i; H = prod G_i mod F;
 * the lane result is prod_k H(x([k]Q)).  With a Brent-Suyama polynomial f
 * (design doc section 14) the roots are x([f(k)]Q) and x([f(jD)]Q).  One curve per lane, polynomial
 * arithmetic from arith/vecpoly.  Roots are computed in vecpoly's
 * radix-2^52 Montgomery domain, independent of the curve backend (radix,
 * RNS or Mersenne).
 *
 * Memory (design doc section 11): per thread, the vecpoly work area (one
 * NTT buffer of K * 2dF vectors), G and H (dF coefficients each), the kept
 * levels of F's product tree (dF coefficients each), 1/F and optionally
 * the transforms of 1/F and F.  setup picks dF, the tree levels kept and
 * whether the transforms are kept, for the least estimated time within
 * the memory limit.
 *
 * Layout of this file: cost / memory model and the automatic choice of B2
 * and of the Brent-Suyama polynomial (sections 11, 13.3, 14.9); setup and
 * per-thread state; x-only Montgomery point arithmetic and batch
 * normalization for the plain roots; the Brent-Suyama tables (affine
 * points, finite differences); avxecm_s2fft_run, the per-thread driver.
 */

#include <sys/time.h>	//for gettimeofday using gcc
#include "avx_ecm.h"
#include "avxecm_stage2fft.h"
#include "../../arith/vecpoly/vecpoly.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef VP_HAVE_KERNELS

// D with phi(D)/2 <= 2^t, densest (GMP-ECM bestd.c, lpo2), t = 1..15
static const uint32_t s2_D[16] = { 0, 12, 30, 60, 120, 240, 510, 1020, 2310,
    4620, 9240, 19110, 39270, 79170, 158340, 324870 };

#define S2_CHUNK 256        // roots normalized per batch inversion

typedef struct
{
    vp_work_t *w;
    uint64_t *X, *Z, *pre;                  // S2_CHUNK coefficients each
    uint64_t *G, *H, *F, **tree;            // dF coefficients each; F is
                                            // H (until block 0) unless lean
    vp_modf_t *mf;
    uint64_t *pt;                           // point scratch: 16 coefficients
    // Brent-Suyama (design doc section 14): tables of M progressions,
    // entry (e, r) at e M + r, e = 0..deg; affine, scaled (X = g x, Y = g^2 y)
    uint64_t *bsx, *bsy;
    uint64_t *bsT;                          // their integer multipliers mod S2_TP
    uint64_t *bden, *binv, *bpre, *blam, *bx3, *bt, *bsa;  // scratch, bsQ each
    uint64_t *ix, *iy, *iT;                 // 2-D setup table, (deg+1)(deg+2)/2
    uint64_t *bci;                          // chain inversion scratch, 16
    uint64_t *bc;                           // constants and temporaries, 20
    double tinit;                           // seconds in the table setups (-v -v)
} s2_thread_t;

struct avxecm_s2fft_s
{
    vp_ctx_t *ctx;
    int t;
    int kstride;                            // tree levels kept: s % kstride == 0
    int lean;                               // no stored transforms of 1/F, F
    int bs;                                 // Brent-Suyama: 0 off, S > 1 x^S,
                                            //   S < -1 Dickson D_|S|(x, -1)
    uint32_t bsM;                           // progressions side by side
    uint32_t bsQ;                           // Brent-Suyama scratch coefficients
    uint32_t dF, D, nbaby, chunk;
    uint32_t *baby;                         // the k in S1
    uint64_t j0;
    uint32_t blocks;
    uint64_t B1, B2;
    double mem, ctxmem;                     // bytes per thread, shared
    s2_thread_t *thr;
    int threads;
    mpz_t N, R;                             // arithmetic modulus, 2^(52n)
};

// Euclid; used for S1 membership (gcd(k, D) = 1)
static uint32_t s2_gcd(uint32_t a, uint32_t b)
{
    while (b) { uint32_t t = a % b; a = b; b = t; }
    return a;
}

// ---- memory and cost model (design doc section 11) ----

// Cost model constants per ISA, in template mulmods (one 8-lane mulmod of
// this build).  Fit with arith/vecpoly/ecm_tests/fit_s2.py to vp_test -b
// (n = 4..40 words, t = 8..13) and yafu stage 1 / stage 2 timings.
//   tc = a t + b + c t^2 / n: G's product tree, mulmods per root
//   mm, mi, ev: H G mod F, 1/F (+ stored transforms), evaluation, per tree
//   c1: stage 1 mulmods per unit of B1 (PRAC, curve types 0/1); c1e the same
//   for the Edwards stage 1 (curve type 2, arith/avxecm_edwards_plan.md F);
//   scale: s2_cost vs measured stage 2
typedef struct
{
    double a, b, c;
    double mm, mm_lean, mi, mi_lean, ev, ev_lean;
    double c1, c1e, scale;
} s2_model_t;

#if defined(USE_AVX512IFMA) || defined(IFMA)
static const s2_model_t s2m = { 6.819, -36.97, 3.336, 0.817, 0.939, 0.785, 0.566, 1.842, 1.906, 15.5, 12.0, 0.90 };
#elif defined(USE_AVX512F)
static const s2_model_t s2m = { 3.445, -7.58, 2.687, 0.721, 0.837, 0.762, 0.549, 1.916, 1.870, 14.2, 11.4, 0.92 };
#else
static const s2_model_t s2m = { 4.219, -14.68, 1.402, 0.669, 0.734, 0.679, 0.492, 1.828, 1.808, 14.5, 11.9, 1.0 };
#endif

// G's product tree (from dF = 2^t roots), in template mulmods per root
static double s2_tcost(int t, int n)
{
    double c = s2m.a * t + s2m.b + s2m.c * t * t / n;
    return c < 8.0 ? 8.0 : c;
}

// number of F's tree levels stored when every k-th one is kept
static int s2_levels(int t, int k)
{
    return (t + k - 1) / k;                 // s = 0, k, 2k, ... < t
}

// ---- Brent-Suyama cost terms (design doc section 14), template mulmods ----
#define S2_NCH      4       // interleaved chains in a batch inversion
#define S2_BS_MAXM  128     // progressions side by side, at most

// one batch inversion's fixed part: 8 lanes x mpz_invert, chain totals
static double s2_bs_inv(int n)
{
    return 8.0 * (4.0 + 160.0 / n) + 3.0 * S2_NCH;
}

// one batched affine addition (6 mulmods + add/sub)
static double s2_bs_add(int n)
{
    return 6.0 + 12.0 / n;
}

// table setup for one family of C progressions with arguments of ~lg bits:
// 1-D: C (deg+1) ladders; 2-D: (deg+1)(deg+2)/2 ladders, then C - 1 steps
// of deg (deg+1)/2 additions.  A ladder: ~10 mulmods per bit of an
// argument^deg, + y-recovery and normalization.
static double s2_bs_famcost(int deg, double C, double lg, int n, int *twod)
{
    double E = deg + 1.0, lad = 10.0 * deg * lg + 20.0;
    double c1 = C * E * lad;
    double c2 = E * (deg + 2.0) / 2.0 * lad
        + (C - 1.0) * (deg * E / 2.0 * s2_bs_add(n) + s2_bs_inv(n));
    if (twod != NULL)
        *twod = c2 < c1;
    return c2 < c1 ? c2 : c1;
}

// all table setups: giant (one family of M) and baby (two families of M/2)
static double s2_bs_init_cost(int t, int n, double blocks, int deg, uint32_t M)
{
    double dF = (double)(1u << t), D = s2_D[t];
    return s2_bs_famcost(deg, M, log2((blocks * dF + 1.0) * D), n, NULL)
        + 2.0 * s2_bs_famcost(deg, M / 2, log2(D / 2.0), n, NULL);
}

// entries of the 2-D setup triangle: Delta_sigma^j Delta_delta^e f, j + e <= deg
static uint32_t s2_bs_tri(int deg)
{
    return (uint32_t)((deg + 1) * (deg + 2) / 2);
}

// progressions side by side: one batch inversion per M roots against the
// table setups (cheap per progression with the 2-D setup)
static uint32_t s2_bs_m(int t, int n, double blocks, int deg)
{
    double dF = (double)(1u << t), roots = dF * (blocks + 1.6), best = 1e300;
    uint32_t M, bM = 2;
    for (M = 2; (M <= S2_BS_MAXM) && (M <= (uint32_t)dF); M *= 2)
    {
        double c = roots * s2_bs_inv(n) / M + s2_bs_init_cost(t, n, blocks, deg, M);
        if (c < best)
        {
            best = c;
            bM = M;
        }
    }
    return bM;
}

// Brent-Suyama scratch coefficients per array
static uint32_t s2_bs_q(int deg, uint32_t M)
{
    uint32_t a = (uint32_t)(deg + 1) * M, b = s2_bs_tri(deg);
    return a > b ? a : b;
}

// bytes per thread (excluding the shared context)
static double s2_mem(int t, int n, int k, int lean, int bs, uint32_t bsM, double *ctxb)
{
    int K = vp_nprimes(n, t + 1), deg = bs > 0 ? bs : -bs;
    double dF = (double)(1u << t), cb = 64.0 * n;
    uint32_t M = bsM ? bsM : (bs ? s2_bs_m(t, n, 8.0, deg) : 0);
    double b = (double)vp_work_bytes(n, K, t + 1)
        + (2 + s2_levels(t, k) + (lean ? 1 : 0)) * dF * cb
        + (3.0 * (dF < S2_CHUNK ? dF : S2_CHUNK) + 16.0) * cb
        + (bs ? (2.0 * (deg + 1) * M + 7.0 * s2_bs_q(deg, M) + 2.0 * s2_bs_tri(deg) + 36.0) * cb
            + 8.0 * ((double)(deg + 1) * M + s2_bs_tri(deg)) : 0.0)
        + (double)vp_modf_bytes(n, K, 1u << t,
            lean ? 0 : VP_MODF_KEEP_INV | VP_MODF_KEEP_F);
    if (ctxb != NULL)
        *ctxb = (double)vp_ctx_bytes(n, K, t + 1);
    return b;
}

// estimated mulmods for the whole stage 2 (blocks blocks)
static double s2_cost(int t, int n, int k, int lean, double blocks, int bs, uint32_t bsM)
{
    double dF = (double)(1u << t), T = s2_tcost(t, n) * dF;
    double groot = 10.0 * dF, broot = 17.0 * dF, init = 0;
    int deg = bs > 0 ? bs : -bs;
    double rebuilt = 0;                     // tree levels rebuilt in the evaluation
    int s;
    for (s = 1; s < t; s++)
        rebuilt += s % k;
    rebuilt /= t;
    // per block: giant roots, G's tree, H G mod F; once: baby roots, F's
    // tree, 1/F, evaluation (+ rebuilt tree levels), in units of a tree
    if (deg > 1)
    {
        // Brent-Suyama (design doc section 14): deg affine additions per
        // root, one batch inversion per M roots; D/6 ~ 1.6 dF baby
        // candidates; table setups
        double per;
        if (bsM == 0)
            bsM = s2_bs_m(t, n, blocks, deg);
        per = deg * s2_bs_add(n) + 1.0 + s2_bs_inv(n) / bsM;
        groot = dF * per;
        broot = 1.6 * dF * per;
        init = s2_bs_init_cost(t, n, blocks, deg, bsM);
    }
    return blocks * (groot + T * (1.0 + (lean ? s2m.mm_lean : s2m.mm)))
        + init + broot + T * (1.0 + (lean ? s2m.mi_lean : s2m.mi) +
            (lean ? s2m.ev_lean : s2m.ev) + rebuilt);
}

// ---- automatic B2 (design doc section 11.3) ----
// Choose (t, k, lean, blocks) for the least predicted time per factor,
// (T1 + T2) / P(B1, B2), at the factor size B1 is usually chosen for.
// T1 is predicted (stage 2 is set up before stage 1 runs): S2_C1 template
// mulmods per unit of B1 (per ISA and stage-1 arithmetic: s2m.c1 or s2m.c1e,
// chosen by avxecm_s2fft_setup).

static double s2_c1 = 0;
#define S2_C1           s2_c1
#define S2_COST_SCALE   s2m.scale
#define S2_RHO_H        0.001
#define S2_RHO_MAX      40.0

// Dickman rho on [0, S2_RHO_MAX], step S2_RHO_H (rho'(u) = -rho(u-1)/u)
static double *s2_rho_table(int *len)
{
    int m = (int)(S2_RHO_MAX / S2_RHO_H) + 2, k1 = (int)(1.0 / S2_RHO_H + 0.5), i;
    double *r = (double *)malloc((size_t)m * sizeof(double));
    for (i = 0; i <= k1 && i < m; i++)
        r[i] = 1.0;
    for (; i < m; i++)
    {
        double u = i * S2_RHO_H;
        r[i] = r[i - 1] - 0.5 * S2_RHO_H * (r[i - 1 - k1] / (u - S2_RHO_H) + r[i - k1] / u);
        if (r[i] < 0) r[i] = 0;
    }
    *len = m;
    return r;
}

// Dickman rho(u) by linear interpolation in the table (0 beyond it)
static double s2_rho(const double *r, int len, double u)
{
    double x, f;
    int i;
    if (u <= 1.0) return 1.0;
    x = u / S2_RHO_H;
    i = (int)x;
    if (i + 1 >= len) return 0.0;
    f = x - i;
    return r[i] * (1 - f) + r[i + 1] * f;
}

// Brent-Suyama: for a prime q > B2 the expected number of non-trivial
// pairs with f(jD) = +-f(k) mod q is pairs r / q, r the number of
// non-trivial x with f(x) = +-f(y); r depends on q's residue class, so the
// catch chance 1 - exp(-pairs r / q) is averaged over r's distribution
// (counted over random primes, arith/vecpoly/ecm_tests/bs_model.py).
typedef struct { int S, cnt; double r[6], w[6]; } s2_rhist_t;
static const s2_rhist_t s2_rhist[] = {
    {  2, 1, { 2 }, { 0.500 } },
    {  3, 1, { 4 }, { 0.467 } },
    {  6, 3, { 2, 4, 10 }, { 0.273, 0.240, 0.227 } },
    { -3, 1, { 4 }, { 0.485 } },
    { -6, 3, { 2, 4, 10 }, { 0.251, 0.262, 0.223 } },
    { -12, 5, { 2, 4, 6, 10, 22 }, { 0.147, 0.262, 0.104, 0.102, 0.122 } },
};

// r distribution of polynomial bs, NULL if not tabulated (no q > B2 term)
static const s2_rhist_t *s2_rhist_get(int bs)
{
    int i;
    for (i = 0; i < (int)(sizeof(s2_rhist) / sizeof(s2_rhist[0])); i++)
        if (s2_rhist[i].S == bs)
            return &s2_rhist[i];
    return NULL;
}

// probability that a group order of size p/23 (Suyama torsion) is
// B1-smooth but for one prime q in (B1, B2]; lnp = ln p.  With a
// Brent-Suyama polynomial bs (and its r distribution known) also one q
// above B2, caught with the chance above (pairs = nbaby x giant steps).
static double s2_prob(const double *r, int len, double lnp, double B1, double B2,
    double pairs, int bs)
{
    double L = lnp - log(23.0), a = log(B1), b = log(B2), h, sum = 0, p;
    const s2_rhist_t *rh = s2_rhist_get(bs);
    int i, j, steps = 100;
    h = (b - a) / steps;
    for (i = 0; i <= steps; i++)
    {
        double l = a + i * h, wt = (i == 0 || i == steps) ? 0.5 : 1.0;
        sum += wt * s2_rho(r, len, (L - l) / a) / l;
    }
    p = s2_rho(r, len, L / a) + sum * h;
    if ((rh != NULL) && (L > b))
    {
        h = (L - b) / steps;
        sum = 0;
        for (i = 0; i <= steps; i++)
        {
            double l = b + i * h, wt = (i == 0 || i == steps) ? 0.5 : 1.0;
            double E1 = pairs / exp(l), c = 0;
            for (j = 0; j < rh->cnt; j++)
                c += rh->w[j] * (1.0 - exp(-E1 * rh->r[j]));
            sum += wt * s2_rho(r, len, (L - l) / a) / l * c;
        }
        p += sum * h;
    }
    return p;
}

// baby steps for D: k < D/2 coprime to D
static double s2_nbaby(uint32_t D)
{
    uint32_t k, c = 0;
    for (k = 1; k < D / 2; k++)
        if (s2_gcd(k, D) == 1)
            c++;
    return (double)c;
}

// polynomials the automatic choice considers (design doc section 14.8)
static const int s2_bs_cand[] = { 0, -3, -6, -12 };
#define S2_BS_NCAND ((int)(sizeof(s2_bs_cand) / sizeof(s2_bs_cand[0])))

// digits of the factors B1 is usually chosen for (GMP-ECM's table,
// interpolated in log B1)
static double s2_target_digits(double B1)
{
    static const double b1[] = { 11e3, 5e4, 25e4, 1e6, 3e6, 11e6, 43e6, 11e7, 26e7, 85e7, 29e8 };
    int i, cnt = (int)(sizeof(b1) / sizeof(b1[0]));
    if (B1 <= b1[0]) return 20.0;
    for (i = 1; i < cnt; i++)
        if (B1 <= b1[i])
            return 20.0 + 5.0 * (i - 1 + log(B1 / b1[i - 1]) / log(b1[i] / b1[i - 1]));
    return 70.0;
}

// largest prime covered by blocks blocks of dF = 2^t giant steps from
// j0 = round(B1/D): the last giant step's jD plus D/2 - 1
static uint64_t s2_covered_B2(int t, uint64_t B1, uint64_t blocks)
{
    uint64_t D = s2_D[t], j0 = (B1 + D / 2) / D;
    if (j0 < 1) j0 = 1;
    return (j0 + blocks * (1u << t) - 1) * D + D / 2 - 1;
}

// least-cost configuration (t, k, lean) for a fixed B2 within the memory
// limit; returns its cost (1e300: none fits)
static double s2_fixed_cfg(int n, uint64_t B1, uint64_t B2, int threads, double cap,
    int bs, uint32_t bsM, int *bt, int *bk, int *bl, double *bmem, double *bctx)
{
    int t, k, lean;
    double best = 1e300;
    for (t = 3; t <= 15; t++)
    {
        double dF = (double)(1u << t), D = s2_D[t];
        double blocks = ceil((double)(B2 - B1) / (dF * D));
        if (B1 < D)                 // need j0 >= 1 and D's primes < B1
            break;
        for (lean = 0; lean <= 1; lean++)
        {
            for (k = 1; k <= 4; k++)
            {
                double ctxb, mem = s2_mem(t, n, k, lean, bs, bsM, &ctxb), cost;
                if (mem + ctxb / threads > cap)
                    continue;
                cost = s2_cost(t, n, k, lean, blocks, bs, bsM);
                if (cost < best)
                {
                    best = cost;
                    *bt = t; *bk = k; *bl = lean;
                    *bmem = mem; *bctx = ctxb;
                }
            }
        }
    }
    return best;
}

// Automatic B2: the configuration and block count with the least predicted
// time per factor.  *bs = AVXECM_BS_AUTO also chooses the polynomial among
// s2_bs_cand.  Returns the chosen B2 (0: none fits).
static uint64_t s2_auto_B2(int n, uint64_t B1, int threads, double cap, int *bs, uint32_t bsM,
    int *bt, int *bk, int *bl, double *ratio, double *pdig)
{
    int len, t, k, lean, c, ncand = (*bs == AVXECM_BS_AUTO) ? S2_BS_NCAND : 1, bbs = 0;
    double *r = s2_rho_table(&len);
    double T1 = S2_C1 * (double)B1, best = 0;
    double lnp;
    uint64_t bB2 = 0;

    *pdig = s2_target_digits((double)B1);
    lnp = *pdig * log(10.0);
    for (c = 0; c < ncand; c++)
    {
        int cs = (*bs == AVXECM_BS_AUTO) ? s2_bs_cand[c] : *bs;
        for (t = 3; t <= 15; t++)
        {
            double nb;
            if (B1 < s2_D[t])
                break;
            nb = s2_nbaby(s2_D[t]);
            for (lean = 0; lean <= 1; lean++)
            {
                for (k = 1; k <= 4; k++)
                {
                    double ctxb, mem = s2_mem(t, n, k, lean, cs, bsM, &ctxb);
                    uint64_t blocks = 1;
                    if (mem + ctxb / threads > cap)
                        continue;
                    while (1)
                    {
                        double T2 = S2_COST_SCALE * s2_cost(t, n, k, lean, (double)blocks, cs, bsM);
                        uint64_t B2c = s2_covered_B2(t, B1, blocks);
                        double e;
                        if ((T2 > 4.0 * T1) || (B2c > (uint64_t)1e15))
                            break;
                        e = s2_prob(r, len, lnp, (double)B1, (double)B2c,
                            nb * (double)blocks * (double)(1u << t), cs) / (T1 + T2);
                        if (e > best)
                        {
                            best = e;
                            bB2 = B2c;
                            bbs = cs;
                            *bt = t; *bk = k; *bl = lean;
                            *ratio = T2 / T1;
                        }
                        blocks = blocks < 8 ? blocks + 1 : blocks + blocks / 8;
                    }
                }
            }
        }
    }
    free(r);
    *bs = bbs;
    return bB2;
}

// fixed B2, bs = AVXECM_BS_AUTO: the polynomial with the least predicted
// time per factor at this B2
static int s2_auto_bs(int n, uint64_t B1, uint64_t B2, int threads, double cap, uint32_t bsM,
    double *pdig)
{
    int len, c, bbs = 0;
    double *r = s2_rho_table(&len);
    double T1 = S2_C1 * (double)B1, best = 0, lnp;

    *pdig = s2_target_digits((double)B1);
    lnp = *pdig * log(10.0);
    for (c = 0; c < S2_BS_NCAND; c++)
    {
        int bt = 0, bk = 1, bl = 0;
        double bmem, bctx, cost, e, dF, blocks;
        cost = s2_fixed_cfg(n, B1, B2, threads, cap, s2_bs_cand[c], bsM, &bt, &bk, &bl, &bmem, &bctx);
        if (bt == 0)
            continue;
        dF = (double)(1u << bt);
        blocks = ceil((double)(B2 - B1) / (dF * s2_D[bt]));
        e = s2_prob(r, len, lnp, (double)B1, (double)B2, s2_nbaby(s2_D[bt]) * blocks * dF,
            s2_bs_cand[c]) / (T1 + S2_COST_SCALE * cost);
        if (e > best)
        {
            best = e;
            bbs = s2_bs_cand[c];
        }
    }
    free(r);
    return bbs;
}

// Shared setup (see avxecm_stage2fft.h).  Order: resolve the polynomial
// option; automatic B2 (with the polynomial when bs is AVXECM_BS_AUTO) or,
// for a fixed B2, the automatic polynomial; then the least-cost (t, k,
// lean) within the memory limit; then D, S1, j0, the block count and the
// covered B2, and the vecpoly context (NTT primes, twiddles, CRT tables)
// shared read-only by all threads.  Per-thread state is allocated lazily
// by the first avxecm_s2fft_run of each thread.
avxecm_s2fft_t *avxecm_s2fft_setup(mpz_t N, uint64_t B1, uint64_t B2,
    int threads, double mem_mb, int auto_b2, int bs, int edw_stage1, int verbose)
{
    avxecm_s2fft_t *s;
    int n, t, bt = 0, bk = 1, bl = 0, bs_auto = (bs == AVXECM_BS_AUTO);
    uint32_t bsM = 0;                       // 0: s2_bs_m per configuration
    double bmem = 0, bctx = 0, cap = mem_mb * 1048576.0;
    double aratio = 0, adig = 0;
    uint32_t j;

    if (mpz_even_p(N))
        return NULL;
    s2_c1 = edw_stage1 ? s2m.c1e : s2m.c1;
    if (threads < 1) threads = 1;
    // Brent-Suyama: |S| <= 1 is the plain stage 2; |S| <= 32
    if (!bs_auto)
    {
        if ((bs >= -1) && (bs <= 1))
            bs = 0;
        if (bs > 32) bs = 32;
        if (bs < -32) bs = -32;
    }
    {
        // testing aid: AVXECM_S2FFT_BSM = progressions side by side
        const char *e = getenv("AVXECM_S2FFT_BSM");
        if ((e != NULL) && (atoi(e) >= 2) && (atoi(e) <= 256))
            bsM = (uint32_t)atoi(e);
        while (bsM & (bsM - 1)) bsM &= bsM - 1;    // power of 2
    }
    n = (int)((mpz_sizeinbase(N, 2) + 3 + 51) / 52);
    if (n < 2) n = 2;

    if (auto_b2)
    {
        uint64_t aB2 = s2_auto_B2(n, B1, threads, cap, &bs, bsM, &bt, &bk, &bl, &aratio, &adig);
        if (aB2 > B1)
        {
            B2 = aB2;
            bmem = s2_mem(bt, n, bk, bl, bs, bsM, &bctx);
        }
        else
        {
            bt = 0;
            auto_b2 = 0;
            if (bs_auto)
                bs = AVXECM_BS_AUTO;
        }
    }
    if (B2 <= B1)
        return NULL;
    if (bs == AVXECM_BS_AUTO)
        bs = s2_auto_bs(n, B1, B2, threads, cap, bsM, &adig);

    if (bt == 0)                    // fixed B2 (or auto found nothing)
        s2_fixed_cfg(n, B1, B2, threads, cap, bs, bsM, &bt, &bk, &bl, &bmem, &bctx);
    // testing aid: AVXECM_S2FFT_CFG="t,k,lean" forces a configuration
    {
        const char *e = getenv("AVXECM_S2FFT_CFG");
        int ft, fk, fl;
        if ((e != NULL) && (sscanf(e, "%d,%d,%d", &ft, &fk, &fl) == 3) &&
            (ft >= 1) && (ft <= 15) && (fk >= 1) && (B1 >= s2_D[ft]))
        {
            bt = ft; bk = fk; bl = (fl != 0);
            bmem = s2_mem(bt, n, bk, bl, bs, bsM, &bctx);
        }
    }
    if (bt == 0)
        return NULL;

    s = (avxecm_s2fft_t *)calloc(1, sizeof(avxecm_s2fft_t));
    s->t = t = bt;
    s->kstride = bk;
    s->lean = bl;
    s->bs = bs;
    s->bsM = bsM;
    s->D = s2_D[t];
    s->dF = 1u << t;
    s->j0 = (B1 + s->D / 2) / s->D;
    if (s->j0 < 1) s->j0 = 1;
    if (bs && (s->bsM == 0))
        s->bsM = s2_bs_m(t, n, (double)((B2 + s->D / 2 - s->j0 * s->D) / ((uint64_t)s->dF * s->D) + 1),
            bs > 0 ? bs : -bs);
    while (s->bsM > (1u << t)) s->bsM >>= 1;   // M divides dF (both powers of 2)
    if (s->bsM < 2) s->bsM = 2;
    s->bsQ = bs ? s2_bs_q(bs > 0 ? bs : -bs, s->bsM) : 0;
    s->mem = bmem;
    s->ctxmem = bctx;
    s->dF = 1u << t;
    s->chunk = s->dF < S2_CHUNK ? s->dF : S2_CHUNK;
    s->D = s2_D[t];
    s->B1 = B1;
    s->threads = threads;
    mpz_init_set(s->N, N);
    mpz_init(s->R);
    mpz_set_ui(s->R, 1);
    mpz_mul_2exp(s->R, s->R, (mp_bitcnt_t)(52 * n));

    s->baby = (uint32_t *)malloc(s->dF * sizeof(uint32_t));
    s->nbaby = 0;
    for (j = 1; j < s->D / 2; j++)
        if (s2_gcd(j, s->D) == 1)
            s->baby[s->nbaby++] = j;

    // giant steps j = j0 .. j0 + blocks*dF - 1 cover (j0 D - D/2, ...)
    s->j0 = (B1 + s->D / 2) / s->D;
    if (s->j0 < 1) s->j0 = 1;
    s->blocks = (uint32_t)((B2 + s->D / 2 - s->j0 * s->D) / ((uint64_t)s->dF * s->D) + 1);
    s->B2 = (s->j0 + (uint64_t)s->blocks * s->dF - 1) * s->D + s->D / 2 - 1;

    s->ctx = vp_ctx_init(N, n, t + 1);
    if (s->ctx == NULL)
    {
        avxecm_s2fft_free(s);
        return NULL;
    }
    s->thr = (s2_thread_t *)calloc((size_t)threads, sizeof(s2_thread_t));

    if (verbose > 0)
    {
        printf("ecm: FFT stage 2: B2 = %" PRIu64 ", D = %u, dF = %u (%u baby roots), "
            "%u block%s, %d NTT primes\n",
            s->B2, s->D, s->dF, s->nbaby, s->blocks, s->blocks > 1 ? "s" : "",
            s->ctx->K);
        printf("ecm: FFT stage 2: %1.1f MB/thread + %1.1f MB shared "
            "(limit %g MB/thread; %d of %d tree levels kept%s)\n",
            s->mem / 1048576.0, s->ctxmem / 1048576.0, mem_mb,
            s2_levels(t, s->kstride), t, s->lean ? ", 1/F and F transforms recomputed" : "");
        if (s->bs)
            printf("ecm: FFT stage 2: Brent-Suyama with %s of degree %d, %u progressions%s\n",
                s->bs > 0 ? "x^S" : "Dickson D_S(x, -1)", s->bs > 0 ? s->bs : -s->bs, s->bsM,
                bs_auto ? " (chosen automatically)" : "");
        else if (bs_auto)
            printf("ecm: FFT stage 2: no Brent-Suyama (chosen automatically)\n");
        if (auto_b2)
            printf("ecm: FFT stage 2: B2 chosen for the memory limit (target p%1.0f, "
                "predicted stage 2 / stage 1 time %1.2f)\n", adig, aratio);
    }
    return s;
}

uint64_t avxecm_s2fft_B2(avxecm_s2fft_t *s)
{
    return s->B2;
}

// frees one thread's work area, polynomials and tables
static void s2_thread_free(avxecm_s2fft_t *s, s2_thread_t *th)
{
    int i;
    if (th->w == NULL) return;
    vp_free(th->X); vp_free(th->Z); vp_free(th->pre);
    vp_free(th->G); vp_free(th->H); vp_free(th->pt);
    if (s->lean) vp_free(th->F);
    if (s->bs)
    {
        vp_free(th->bsx); vp_free(th->bsy); free(th->bsT);
        vp_free(th->bden); vp_free(th->binv); vp_free(th->bpre);
        vp_free(th->blam); vp_free(th->bx3); vp_free(th->bt); vp_free(th->bsa);
        vp_free(th->ix); vp_free(th->iy); free(th->iT);
        vp_free(th->bci); vp_free(th->bc);
    }
    for (i = 0; i < s->t; i++)
        if (th->tree[i] != NULL) vp_free(th->tree[i]);
    free(th->tree);
    vp_work_free(th->w);
    th->w = NULL;
}

void avxecm_s2fft_free(avxecm_s2fft_t *s)
{
    int i;
    if (s == NULL) return;
    if (s->thr != NULL)
    {
        for (i = 0; i < s->threads; i++)
            s2_thread_free(s, &s->thr[i]);
        free(s->thr);
    }
    vp_ctx_free(s->ctx);
    free(s->baby);
    mpz_clear(s->N);
    mpz_clear(s->R);
    free(s);
}

// one thread's state: the vecpoly work area (NTT buffers), the root batch
// (X, Z, prefix products), G, H (F lives in H unless lean), the kept tree
// levels and, with Brent-Suyama, its tables and scratch (sizes: s2_mem)
static void s2_thread_init(avxecm_s2fft_t *s, s2_thread_t *th)
{
    int n = s->ctx->n, i;
    th->w = vp_work_init(s->ctx);
    th->X = vp_alloc(s->chunk, n);
    th->Z = vp_alloc(s->chunk, n);
    th->pre = vp_alloc(s->chunk, n);
    th->G = vp_alloc(s->dF, n);
    th->H = vp_alloc(s->dF, n);
    th->F = s->lean ? vp_alloc(s->dF, n) : th->H;
    th->pt = vp_alloc(16, n);
    if (s->bs)
    {
        int deg = s->bs > 0 ? s->bs : -s->bs;
        size_t e = (size_t)s->bsM * (size_t)(deg + 1), q = s->bsQ, tri = s2_bs_tri(deg);
        th->bsx = vp_alloc(e, n);
        th->bsy = vp_alloc(e, n);
        th->bsT = (uint64_t *)calloc(e, sizeof(uint64_t));
        th->bden = vp_alloc(q, n);
        th->binv = vp_alloc(q, n);
        th->bpre = vp_alloc(q, n);
        th->blam = vp_alloc(q, n);
        th->bx3 = vp_alloc(q, n);
        th->bt = vp_alloc(q, n);
        th->bsa = vp_alloc(q, n);
        th->ix = vp_alloc(tri, n);
        th->iy = vp_alloc(tri, n);
        th->iT = (uint64_t *)calloc(tri, sizeof(uint64_t));
        th->bci = vp_alloc(16, n);
        th->bc = vp_alloc(20, n);
    }
    th->tree = (uint64_t **)malloc((size_t)s->t * sizeof(uint64_t *));
    for (i = 0; i < s->t; i++)
        th->tree[i] = (i % s->kstride == 0) ? vp_alloc(s->dF, n) : NULL;
}

// ---- x-only Montgomery curve arithmetic in the vecpoly domain ----
// points are (X, Z) pairs of single coefficients; scratch from th->pt

#define PC(i) VP_COEF(ctx, th->pt, (i))

typedef struct { uint64_t *X, *Z; } s2_pt;

// one coefficient (8 lanes x n words)
static void s2_copy(vp_ctx_t *ctx, uint64_t *d, const uint64_t *a)
{
    memcpy(d, a, (size_t)8 * (size_t)ctx->n * sizeof(uint64_t));
}

// R = 2P.  a24 = (A+2)/4.  R may alias P.
static void s2_dbl(s2_thread_t *th, s2_pt R, s2_pt P, const uint64_t *a24)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    uint64_t *t1 = PC(0), *t2 = PC(1), *t3 = PC(2), *t4 = PC(3);
    vpk_add(ctx, t1, P.X, P.Z, 1);
    vpk_sub(ctx, t2, P.X, P.Z, 1);
    vp_cmul(w, t1, t1, t1, 1);              // (x+z)^2
    vp_cmul(w, t2, t2, t2, 1);              // (x-z)^2
    vpk_sub(ctx, t3, t1, t2, 1);            // 4xz
    vp_cmul(w, t4, t3, a24, 1);
    vpk_add(ctx, t4, t4, t2, 1);
    vp_cmul(w, R.X, t1, t2, 1);
    vp_cmul(w, R.Z, t3, t4, 1);
}

// R = P + Q, D = P - Q.  R may alias P or Q, not D.
static void s2_add(s2_thread_t *th, s2_pt R, s2_pt P, s2_pt Q, s2_pt D)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    uint64_t *t1 = PC(0), *t2 = PC(1), *t3 = PC(2), *t4 = PC(3);
    vpk_sub(ctx, t1, P.X, P.Z, 1);
    vpk_add(ctx, t2, Q.X, Q.Z, 1);
    vp_cmul(w, t1, t1, t2, 1);              // u
    vpk_add(ctx, t3, P.X, P.Z, 1);
    vpk_sub(ctx, t4, Q.X, Q.Z, 1);
    vp_cmul(w, t3, t3, t4, 1);              // v
    vpk_add(ctx, t2, t1, t3, 1);
    vpk_sub(ctx, t4, t1, t3, 1);
    vp_cmul(w, t2, t2, t2, 1);
    vp_cmul(w, t4, t4, t4, 1);
    vp_cmul(w, R.X, D.Z, t2, 1);
    vp_cmul(w, R.Z, D.X, t4, 1);
}

// R = [k]P, k >= 1 (Montgomery ladder).  Uses PC(4..7).
static void s2_mul(s2_thread_t *th, s2_pt R, s2_pt P, uint64_t k, const uint64_t *a24)
{
    vp_ctx_t *ctx = th->w->ctx;
    s2_pt R0 = { PC(4), PC(5) }, R1 = { PC(6), PC(7) };
    int b = 63;
    s2_copy(ctx, R0.X, P.X);
    s2_copy(ctx, R0.Z, P.Z);
    if (k > 1)
    {
        s2_dbl(th, R1, P, a24);
        while (((k >> b) & 1) == 0) b--;
        for (b = b - 1; b >= 0; b--)
        {
            if ((k >> b) & 1)
            {
                s2_add(th, R0, R1, R0, P);
                s2_dbl(th, R1, R1, a24);
            }
            else
            {
                s2_add(th, R1, R1, R0, P);
                s2_dbl(th, R0, R0, a24);
            }
        }
    }
    s2_copy(ctx, R.X, R0.X);
    s2_copy(ctx, R.Z, R0.Z);
}

// dest[i] = X[i]/Z[i], i < m <= chunk (Montgomery batch inversion per
// lane).  A lane whose Z product is not invertible has found a factor: it
// is folded into found[lane] and the lane continues with garbage.
static void s2_normalize(avxecm_s2fft_t *s, s2_thread_t *th, uint32_t m,
    uint64_t *dest, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    uint64_t *I = PC(8), *T = PC(9);
    uint32_t i;
    int l;

    s2_copy(ctx, VP_COEF(ctx, th->pre, 0), VP_COEF(ctx, th->Z, 0));
    for (i = 1; i < m; i++)
        vp_cmul(w, VP_COEF(ctx, th->pre, i), VP_COEF(ctx, th->pre, i - 1),
            VP_COEF(ctx, th->Z, i), 1);
    for (l = 0; l < 8; l++)
    {
        // c = c_p R; want c_p^-1 R = c^-1 R^2
        vp_get(ctx, t1, VP_COEF(ctx, th->pre, m - 1), l);
        mpz_mod(t1, t1, s->N);
        if (mpz_invert(t2, t1, s->N) == 0)
        {
            mpz_mul(found[l], found[l], t1);
            mpz_mod(found[l], found[l], s->N);
            mpz_set_ui(t2, 1);
        }
        mpz_mul(t2, t2, s->R);
        mpz_mul(t2, t2, s->R);
        mpz_mod(t2, t2, s->N);
        vp_set(ctx, I, t2, l);
    }
    for (i = m - 1; i > 0; i--)
    {
        vp_cmul(w, T, I, VP_COEF(ctx, th->pre, i - 1), 1);         // 1/Z_i
        vp_cmul(w, I, I, VP_COEF(ctx, th->Z, i), 1);
        vp_cmul(w, VP_COEF(ctx, dest, i), VP_COEF(ctx, th->X, i), T, 1);
    }
    vp_cmul(w, VP_COEF(ctx, dest, 0), VP_COEF(ctx, th->X, 0), I, 1);
}

// ---- Brent-Suyama roots (design doc section 14) ----
// Roots x([f(i)]Q) for f(x) = x^S (S > 1) or the Dickson polynomial
// D_|S|(x, -1) (S < -1), i over arithmetic progressions, by finite
// differences of affine points.  Q gets a y-coordinate by moving to the
// curve g y^2 = x^3 + A x^2 + x, g = x0^3 + A x0^2 + x0, Q = (x0, 1): this
// is the stage-1 curve or its quadratic twist, with the same x-coordinates
// for every multiple of Q, which is all stage 2 uses.  M progressions run
// side by side so each step's M S affine additions share one batch
// inversion (one mpz_invert per lane).

// f(x) for the configured polynomial, exact
static void s2_fpoly(int S, mpz_t r, const mpz_t x)
{
    if (S > 1)
        mpz_pow_ui(r, x, (unsigned long)S);
    else
    {
        // D_0 = 2, D_1 = x, D_n = x D_{n-1} + D_{n-2}   (alpha = -1)
        mpz_t d0, d1, t;
        int i;
        mpz_init_set_ui(d0, 2);
        mpz_init_set(d1, x);
        mpz_init(t);
        for (i = 1; i < -S; i++)
        {
            mpz_mul(t, x, d1);
            mpz_add(t, t, d0);
            mpz_swap(d0, d1);
            mpz_swap(d1, t);
        }
        mpz_set(r, d1);
        mpz_clear(d0);
        mpz_clear(d1);
        mpz_clear(t);
    }
}

// out[e] = Delta^e f(i0), step delta, e = 0..deg
static void s2_findiff(int S, int deg, mpz_t *out, uint64_t i0, uint64_t delta)
{
    mpz_t x, d;
    int t, e;
    mpz_init(x);
    mpz_init(d);
    for (t = 0; t <= deg; t++)
    {
        mpz_set_ui(d, (unsigned long)delta);
        mpz_mul_ui(d, d, (unsigned long)t);
        mpz_add_ui(x, d, (unsigned long)i0);
        s2_fpoly(S, out[t], x);
    }
    for (e = 1; e <= deg; e++)
        for (t = deg; t >= e; t--)
            mpz_sub(out[t], out[t], out[t - 1]);
    mpz_clear(x);
    mpz_clear(d);
}

// out[i] = 1 / in[i], i < Q (Montgomery's batch inversion per lane, over
// S2_NCH interleaved chains so the products are vector calls; pre: Q
// coefficients; out must not overlap in).  A lane whose product is not
// invertible is folded into found[lane] (a stage-2 find when it shares a
// factor with N) and continues with garbage.
static void s2_inv_batch(avxecm_s2fft_t *s, s2_thread_t *th, uint32_t Q,
    const uint64_t *in, uint64_t *out, uint64_t *pre, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    uint32_t c = Q < S2_NCH ? Q : S2_NCH, i, j, ri, rows = (Q + c - 1) / c;
    uint64_t *tot = VP_COEF(ctx, th->bci, 0), *tp = VP_COEF(ctx, th->bci, S2_NCH);
    uint64_t *I = VP_COEF(ctx, th->bci, 2 * S2_NCH), *x = VP_COEF(ctx, th->bci, 3 * S2_NCH);
    int l;

    memcpy(pre, in, (size_t)c * 8 * (size_t)ctx->n * sizeof(uint64_t));
    for (i = c; i < Q; i += c)
        vp_cmul(w, VP_COEF(ctx, pre, i), VP_COEF(ctx, pre, i - c), VP_COEF(ctx, in, i),
            Q - i < c ? Q - i : c);
    // chain totals and their product
    for (j = 0; j < c; j++)
        s2_copy(ctx, VP_COEF(ctx, tot, j), VP_COEF(ctx, pre, j + c * ((Q - 1 - j) / c)));
    s2_copy(ctx, tp, tot);
    for (j = 1; j < c; j++)
        vp_cmul(w, VP_COEF(ctx, tp, j), VP_COEF(ctx, tp, j - 1), VP_COEF(ctx, tot, j), 1);
    for (l = 0; l < 8; l++)
    {
        // v = v_p R; want v_p^-1 R = v^-1 R^2
        vp_get(ctx, t1, VP_COEF(ctx, tp, c - 1), l);
        mpz_mod(t1, t1, s->N);
        if (mpz_invert(t2, t1, s->N) == 0)
        {
            mpz_mul(found[l], found[l], t1);
            mpz_mod(found[l], found[l], s->N);
            mpz_set_ui(t2, 1);
        }
        mpz_mul(t2, t2, s->R);
        mpz_mul(t2, t2, s->R);
        mpz_mod(t2, t2, s->N);
        vp_set(ctx, x, t2, l);
    }
    for (j = c - 1; j > 0; j--)
    {
        vp_cmul(w, VP_COEF(ctx, I, j), x, VP_COEF(ctx, tp, j - 1), 1);
        vp_cmul(w, x, x, VP_COEF(ctx, tot, j), 1);
    }
    s2_copy(ctx, I, x);
    // I[j] = 1 / (chain j's product so far); walk the rows back
    for (ri = rows - 1; ri >= 1; ri--)
    {
        uint32_t m;
        i = ri * c;
        m = Q - i < c ? Q - i : c;
        vp_cmul(w, VP_COEF(ctx, out, i), I, VP_COEF(ctx, pre, i - c), m);
        vp_cmul(w, I, I, VP_COEF(ctx, in, i), m);
    }
    memcpy(out, I, (size_t)c * 8 * (size_t)ctx->n * sizeof(uint64_t));
}

// x-only ladder with an mpz scalar c >= 1: R0 = [c]P, R1 = [c+1]P
static void s2_mul2(s2_thread_t *th, s2_pt R0o, s2_pt R1o, s2_pt P, const mpz_t c,
    const uint64_t *a24)
{
    vp_ctx_t *ctx = th->w->ctx;
    s2_pt R0 = { PC(4), PC(5) }, R1 = { PC(6), PC(7) };
    long b = (long)mpz_sizeinbase(c, 2) - 1;
    s2_copy(ctx, R0.X, P.X);
    s2_copy(ctx, R0.Z, P.Z);
    s2_dbl(th, R1, P, a24);
    for (b = b - 1; b >= 0; b--)
    {
        if (mpz_tstbit(c, (mp_bitcnt_t)b))
        {
            s2_add(th, R0, R1, R0, P);
            s2_dbl(th, R1, R1, a24);
        }
        else
        {
            s2_add(th, R1, R1, R0, P);
            s2_dbl(th, R0, R0, a24);
        }
    }
    s2_copy(ctx, R0o.X, R0.X); s2_copy(ctx, R0o.Z, R0.Z);
    s2_copy(ctx, R1o.X, R1.X); s2_copy(ctx, R1o.Z, R1.Z);
}

// BS curve constants, in th->bc: x0, A, 2A, g, 2g; the scaled curve
// Y^2 = X^3 + a2 X^2 + a4 X (X = g x, Y = g^2 y): a2 = g A, a4 = g^2, 2 a2
#define BC_X0(th)  VP_COEF(ctx, (th)->bc, 0)
#define BC_A(th)   VP_COEF(ctx, (th)->bc, 1)
#define BC_A2(th)  VP_COEF(ctx, (th)->bc, 2)
#define BC_G(th)   VP_COEF(ctx, (th)->bc, 3)
#define BC_G2(th)  VP_COEF(ctx, (th)->bc, 4)
#define BC_W2(th)  VP_COEF(ctx, (th)->bc, 5)
#define BC_W4(th)  VP_COEF(ctx, (th)->bc, 6)
#define BC_W22(th) VP_COEF(ctx, (th)->bc, 7)
#define BC_T(th, i) VP_COEF(ctx, (th)->bc, 8 + (i))     // 12 temporaries

// The tables' integer multipliers are tracked modulo the prime S2_TP to
// detect equal points (not mod 2^64: high differences such as
// deg! delta^deg are divisible by large powers of 2).
#define S2_TP 0x1FFFFFFFFFFFFFFFull     // 2^61 - 1

static uint64_t s2_tp_add(uint64_t a, uint64_t b)
{
    uint64_t c = a + b;                 // < 2^62
    return c >= S2_TP ? c - S2_TP : c;
}

// affine points [c_i]Q, i < cnt, scaled, into X[], Y[] (contiguous); the
// multipliers mod S2_TP into T[]
static void s2_bs_points(avxecm_s2fft_t *s, s2_thread_t *th, uint32_t cnt, mpz_t *c,
    uint64_t *X, uint64_t *Y, uint64_t *T, s2_pt Q, const uint64_t *a24,
    mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    s2_pt R0 = { BC_T(th, 0), BC_T(th, 1) }, R1 = { BC_T(th, 2), BC_T(th, 3) };
    uint64_t *u = BC_T(th, 4), *v = BC_T(th, 5), *z11 = BC_T(th, 6), *z12 = BC_T(th, 7);
    uint64_t *t3 = BC_T(th, 8);
    uint32_t q;

    for (q = 0; q < cnt; q++)
    {
        uint64_t *X1 = R0.X, *Z1 = R0.Z, *X2 = R1.X, *Z2 = R1.Z;
        mpz_set_ui(t1, 1);
        mpz_mul_2exp(t1, t1, 61);
        mpz_sub_ui(t1, t1, 1);
        mpz_mod(t2, c[q], t1);
        T[q] = (uint64_t)mpz_getlimbn(t2, 0);
        s2_mul2(th, R0, R1, Q, c[q], a24);
        // Okeya-Sakurai: y = [Z2((x0 X1 + Z1)(x0 Z1 + X1 + 2A Z1) - 2A Z1^2)
        //   - (x0 Z1 - X1)^2 X2] / (2 g Z1^2 Z2),  x = X1 2g Z1 Z2 / (same)
        vp_cmul(w, u, BC_X0(th), X1, 1);
        vpk_add(ctx, u, u, Z1, 1);                      // x0 X1 + Z1
        vp_cmul(w, v, BC_X0(th), Z1, 1);                // x0 Z1
        vp_cmul(w, t3, BC_A2(th), Z1, 1);
        vpk_add(ctx, t3, t3, X1, 1);
        vpk_add(ctx, t3, t3, v, 1);                     // x0 Z1 + X1 + 2A Z1
        vp_cmul(w, u, u, t3, 1);
        vp_csqr(w, z11, Z1, 1);                         // Z1^2
        vp_cmul(w, t3, BC_A2(th), z11, 1);
        vpk_sub(ctx, u, u, t3, 1);
        vp_cmul(w, u, u, Z2, 1);                        // Z2 (...)
        vpk_sub(ctx, v, v, X1, 1);
        vp_csqr(w, v, v, 1);
        vp_cmul(w, v, v, X2, 1);                        // (x0 Z1 - X1)^2 X2
        vpk_sub(ctx, u, u, v, 1);
        vp_cmul(w, VP_COEF(ctx, Y, q), u, BC_W4(th), 1);    // g^2 y numerator
        vp_cmul(w, z12, Z1, Z2, 1);
        vp_cmul(w, z12, z12, BC_G2(th), 1);             // 2g Z1 Z2
        vp_cmul(w, VP_COEF(ctx, th->bden, q), z12, Z1, 1);  // W = 2g Z1^2 Z2
        vp_cmul(w, z12, z12, BC_G(th), 1);
        vp_cmul(w, VP_COEF(ctx, X, q), X1, z12, 1);     // g x numerator
    }
    s2_inv_batch(s, th, cnt, th->bden, th->binv, th->bpre, found, t1, t2);
    vp_cmul(w, X, X, th->binv, cnt);
    vp_cmul(w, Y, Y, th->binv, cnt);
}

// (X1, Y1) <- (X1, Y1) + (X2, Y2) elementwise over nr rows of cnt[i]
// points (X2 may be X1 at a higher offset: rows are read before they are
// written), sharing one batch inversion; multipliers T1 += T2.  The
// multipliers are positive integers (positive coefficients, positive
// arguments); equal ones mean equal points: a doubling (2 Y denominator).
static void s2_aff_add(avxecm_s2fft_t *s, s2_thread_t *th, int nr,
    uint64_t **X1, uint64_t **Y1, uint64_t **X2, uint64_t **Y2,
    uint64_t **T1, uint64_t **T2, const uint32_t *cnt, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    uint32_t off = 0, q, m;
    int i;

    for (i = 0; i < nr; i++)
    {
        m = cnt[i];
        vpk_sub(ctx, VP_COEF(ctx, th->bden, off), X2[i], X1[i], m);
        for (q = 0; q < m; q++)
            if (T1[i][q] == T2[i][q])
                vpk_add(ctx, VP_COEF(ctx, th->bden, off + q), VP_COEF(ctx, Y1[i], q),
                    VP_COEF(ctx, Y1[i], q), 1);
        off += m;
    }
    s2_inv_batch(s, th, off, th->bden, th->binv, th->bpre, found, t1, t2);
    off = 0;
    for (i = 0; i < nr; i++)
    {
        uint64_t *lam = VP_COEF(ctx, th->blam, off), *x3 = VP_COEF(ctx, th->bx3, off);
        uint64_t *tt = VP_COEF(ctx, th->bt, off);
        m = cnt[i];
        vpk_sub(ctx, lam, Y2[i], Y1[i], m);             // (Y2 - Y1) / (X2 - X1)
        for (q = 0; q < m; q++)
            if (T1[i][q] == T2[i][q])
            {
                // (3 X^2 + 2 a2 X + a4) / (2 Y)
                uint64_t *l = VP_COEF(ctx, lam, q), *x = VP_COEF(ctx, X1[i], q);
                uint64_t *u = VP_COEF(ctx, tt, q);
                vp_csqr(w, l, x, 1);
                vpk_add(ctx, u, l, l, 1);
                vpk_add(ctx, l, l, u, 1);
                vp_cmul(w, u, BC_W22(th), x, 1);
                vpk_add(ctx, l, l, u, 1);
                vpk_add(ctx, l, l, BC_W4(th), 1);
            }
        vp_cmul(w, lam, lam, VP_COEF(ctx, th->binv, off), m);
        vp_csqr(w, x3, lam, m);
        vpk_add(ctx, tt, X1[i], X2[i], m);
        vpk_sub(ctx, x3, x3, tt, m);
        vpk_sub(ctx, x3, x3, th->bsa, m);               // lam^2 - a2 - X1 - X2
        vpk_sub(ctx, tt, X1[i], x3, m);
        vp_cmul(w, tt, tt, lam, m);
        vpk_sub(ctx, Y1[i], tt, Y1[i], m);              // lam (X1 - X3) - Y1
        memcpy(X1[i], x3, (size_t)m * 8 * (size_t)ctx->n * sizeof(uint64_t));
        for (q = 0; q < m; q++)
            T1[i][q] = s2_tp_add(T1[i][q], T2[i][q]);
        off += m;
    }
}

// Tables of M progressions in nfam families of C = M / nfam: progression
// r = f C + r' starts at base[f] + r' rstep, step delta; entry (e, r) =
// [Delta_delta^e f(i_r)]Q at e M + r in bsx, bsy, bsT.  1-D setup: a
// ladder per entry.  2-D setup (cheaper for C > ~deg/2): g_e = Delta_delta^e f
// is a polynomial of degree deg - e, so the entries of row e over r' are a
// finite-difference table in steps of rstep; one ladder per entry of that
// triangle (row e: deg - e + 1), then C - 1 steps of all rows at once.
static void s2_bs_init(avxecm_s2fft_t *s, s2_thread_t *th, int nfam, const uint64_t *base,
    uint64_t rstep, uint64_t delta, s2_pt Q, const uint64_t *a24, mpz_t *found,
    mpz_t t1, mpz_t t2)
{
    vp_ctx_t *ctx = th->w->ctx;
    int S = s->bs, twod, f;
    uint32_t deg = (uint32_t)(S > 0 ? S : -S), e, k;
    uint32_t M = s->bsM, C = M / (uint32_t)nfam, E = deg + 1, tri = s2_bs_tri((int)deg), r, nc;
    double lg = log2((double)base[nfam - 1] + (double)C * (double)rstep + (double)deg * (double)delta);
    mpz_t *c, *g;
    size_t cb = (size_t)8 * (size_t)ctx->n * sizeof(uint64_t);
    struct timeval tv0, tv1;

    gettimeofday(&tv0, NULL);
    s2_bs_famcost((int)deg, C, lg, ctx->n, &twod);
    nc = twod ? tri : E * M;
    c = (mpz_t *)malloc((size_t)nc * sizeof(mpz_t));
    g = (mpz_t *)malloc((size_t)E * E * sizeof(mpz_t));
    for (r = 0; r < nc; r++) mpz_init(c[r]);
    for (r = 0; r < E * E; r++) mpz_init(g[r]);

    if (!twod)
    {
        for (f = 0; f < nfam; f++)
            for (r = 0; r < C; r++)
            {
                s2_findiff(S, (int)deg, g, base[f] + r * rstep, delta);
                for (e = 0; e <= deg; e++)
                    mpz_set(c[e * M + (uint32_t)f * C + r], g[e]);
            }
        s2_bs_points(s, th, E * M, c, th->bsx, th->bsy, th->bsT, Q, a24, found, t1, t2);
    }
    else
    {
        uint64_t *X1[64], *Y1[64], *X2[64], *Y2[64], *T1[64], *T2[64];
        uint32_t cnt[64], off[64];
        for (e = 0, r = 0; e <= deg; e++)
        {
            off[e] = r;
            r += deg - e + 1;
        }
        for (e = 0; e < deg; e++)
        {
            X1[e] = VP_COEF(ctx, th->ix, off[e]);  X2[e] = VP_COEF(ctx, th->ix, off[e] + 1);
            Y1[e] = VP_COEF(ctx, th->iy, off[e]);  Y2[e] = VP_COEF(ctx, th->iy, off[e] + 1);
            T1[e] = th->iT + off[e];               T2[e] = th->iT + off[e] + 1;
            cnt[e] = deg - e;
        }
        for (f = 0; f < nfam; f++)
        {
            // g[k E + e] = Delta_delta^e f(base + k rstep), then differences
            // over k: entry (e, j) of the triangle = Delta_rstep^j of row e
            for (k = 0; k <= deg; k++)
                s2_findiff(S, (int)(deg - k), g + k * E, base[f] + (uint64_t)k * rstep, delta);
            for (e = 0; e <= deg; e++)
            {
                uint32_t L = deg - e + 1, j;
                for (k = 0; k < L; k++)
                    mpz_set(c[off[e] + k], g[k * E + e]);
                for (j = 1; j < L; j++)
                    for (k = L - 1; k >= j; k--)
                        mpz_sub(c[off[e] + k], c[off[e] + k], c[off[e] + k - 1]);
            }
            s2_bs_points(s, th, tri, c, th->ix, th->iy, th->iT, Q, a24, found, t1, t2);
            for (r = 0; r < C; r++)
            {
                uint32_t idx = (uint32_t)f * C + r;
                for (e = 0; e <= deg; e++)
                {
                    memcpy(VP_COEF(ctx, th->bsx, e * M + idx), VP_COEF(ctx, th->ix, off[e]), cb);
                    memcpy(VP_COEF(ctx, th->bsy, e * M + idx), VP_COEF(ctx, th->iy, off[e]), cb);
                    th->bsT[e * M + idx] = th->iT[off[e]];
                }
                if (r + 1 < C)
                    s2_aff_add(s, th, (int)deg, X1, Y1, X2, Y2, T1, T2, cnt, found, t1, t2);
            }
        }
    }
    for (r = 0; r < nc; r++) mpz_clear(c[r]);
    for (r = 0; r < E * E; r++) mpz_clear(g[r]);
    free(c);
    free(g);
    gettimeofday(&tv1, NULL);
    th->tinit += ytools_difftime(&tv0, &tv1);
}

// one step of all M tables: (e, r) <- (e, r) + (e + 1, r), e < deg
static void s2_bs_step(avxecm_s2fft_t *s, s2_thread_t *th, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_ctx_t *ctx = th->w->ctx;
    uint32_t M = s->bsM, Q = (uint32_t)(s->bs > 0 ? s->bs : -s->bs) * M;
    uint64_t *X1 = th->bsx, *Y1 = th->bsy, *X2 = VP_COEF(ctx, th->bsx, M);
    uint64_t *Y2 = VP_COEF(ctx, th->bsy, M), *T1 = th->bsT, *T2 = th->bsT + M;
    s2_aff_add(s, th, 1, &X1, &Y1, &X2, &Y2, &T1, &T2, &Q, found, t1, t2);
}

// baby roots X([f(k)]Q), k in S1, into dest (nbaby coefficients):
// families k = c + 6 (r' + h m), c in {1, 5}, r' < h = M/2
static void s2_bs_baby(avxecm_s2fft_t *s, s2_thread_t *th, uint64_t *dest, s2_pt Q,
    const uint64_t *a24, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_ctx_t *ctx = th->w->ctx;
    uint32_t M = s->bsM, h = M / 2, r, idx = 0;
    uint64_t base[2] = { 1, 5 }, m, delta = 6 * (uint64_t)h;
    s2_bs_init(s, th, 2, base, 6, delta, Q, a24, found, t1, t2);
    for (m = 0; ; m++)
    {
        int any = 0;
        for (r = 0; (r < M) && (idx < s->nbaby); r++)
        {
            uint64_t k = base[r / h] + 6 * (uint64_t)(r % h) + m * delta;
            if (k < s->D / 2)
            {
                any = 1;
                if (s2_gcd((uint32_t)k, s->D) == 1)
                    s2_copy(ctx, VP_COEF(ctx, dest, idx++), VP_COEF(ctx, th->bsx, r));
            }
        }
        if (!any || idx == s->nbaby)
            break;
        s2_bs_step(s, th, found, t1, t2);
    }
}

// Stage 2 for the 8 curves of one thread (called from
// vec_ecm_stage2_fft_work_fcn after stage 1).  Steps:
//  1. Q and (A+2)/4 from the backend's Montgomery form into vecpoly's
//     (x R_b^-1 R mod N per lane); with Brent-Suyama also affine x0 and the
//     scaled-curve constants.
//  2. Baby roots (x-only xADD chain + batch normalization, or the
//     Brent-Suyama tables), padded to dF; F's product tree; 1/F (+ its
//     transforms unless lean).
//  3. Per block: dF giant roots, G = their product (in place), H = G mod F
//     (block 0) or H G mod F.
//  4. H evaluated at the baby roots down F's tree; the product of the
//     values times any factors folded in by failed inversions is the
//     lane's result, written to work->stg2acc for the caller's gcd.
int avxecm_s2fft_run(avxecm_s2fft_t *s, int tid, vec_monty_t *mdata, ecm_work *work,
    ecm_pt *P, int verbose)
{
    s2_thread_t *th = &s->thr[tid];
    vp_ctx_t *ctx = s->ctx;
    vp_work_t *w;
    mpz_t found[8], t1, t2, Rbi;
    uint64_t *a24;
    s2_pt Q, DQ;
    uint32_t i, b, k, c;
    int l;
    struct timeval t0, t1v;
    double tb = 0, tg = 0, tf = 0, te = 0;

    if (th->w == NULL)
        s2_thread_init(s, th);
    w = th->w;
    th->tinit = 0;

    for (l = 0; l < 8; l++) mpz_init_set_ui(found[l], 1);
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(Rbi);

    // point scratch: PC(0..9) used by the helpers above
    Q.X = PC(10); Q.Z = PC(11);
    a24 = PC(12);

    // stage 1 output and (A+2)/4 from the backend (Montgomery wrt R_b):
    // into this domain, x R_b^-1 R mod N
    for (l = 0; l < 8; l++)
    {
        vec_lane_get(mdata, t1, mdata->one, l);
        mpz_mod(t1, t1, s->N);
        if (mpz_invert(Rbi, t1, s->N) == 0)
            mpz_set_ui(Rbi, 1);
        mpz_mul(Rbi, Rbi, s->R);
        mpz_mod(Rbi, Rbi, s->N);
        vec_lane_get(mdata, t1, P->X, l);
        mpz_mul(t1, t1, Rbi); mpz_mod(t1, t1, s->N);
        vp_set(ctx, Q.X, t1, l);
        vec_lane_get(mdata, t1, P->Z, l);
        mpz_mul(t1, t1, Rbi); mpz_mod(t1, t1, s->N);
        vp_set(ctx, Q.Z, t1, l);
        vec_lane_get(mdata, t1, work->s, l);
        mpz_mul(t1, t1, Rbi); mpz_mod(t1, t1, s->N);
        vp_set(ctx, a24, t1, l);
    }

    if (s->bs)
    {
        // affine Q = (x0, 1) on g y^2 = x^3 + A x^2 + x; constants in th->bc
        uint64_t *two = BC_T(th, 0);
        for (l = 0; l < 8; l++)
        {
            vp_get(ctx, t1, Q.X, l);
            vp_get(ctx, t2, Q.Z, l);
            if (mpz_invert(Rbi, t2, s->N) == 0)
            {
                // stage 1 already reached the identity mod a factor
                mpz_mul(found[l], found[l], t2);
                mpz_mod(found[l], found[l], s->N);
                mpz_set_ui(t1, 2);
                mpz_set_ui(Rbi, 1);
            }
            mpz_mul(t1, t1, Rbi);
            mpz_mul(t1, t1, s->R);
            mpz_mod(t1, t1, s->N);                  // x0 R
            vp_set(ctx, BC_X0(th), t1, l);
        }
        s2_copy(ctx, Q.X, BC_X0(th));
        s2_copy(ctx, Q.Z, w->one);
        vpk_add(ctx, two, w->one, w->one, 1);
        vpk_add(ctx, BC_A(th), a24, a24, 1);
        vpk_add(ctx, BC_A(th), BC_A(th), BC_A(th), 1);
        vpk_sub(ctx, BC_A(th), BC_A(th), two, 1);   // A = 4 (A+2)/4 - 2
        vpk_add(ctx, BC_A2(th), BC_A(th), BC_A(th), 1);
        vpk_add(ctx, BC_G(th), BC_X0(th), BC_A(th), 1);
        vp_cmul(w, BC_G(th), BC_G(th), BC_X0(th), 1);
        vpk_add(ctx, BC_G(th), BC_G(th), w->one, 1);
        vp_cmul(w, BC_G(th), BC_G(th), BC_X0(th), 1);   // g = x0 (x0 (x0 + A) + 1)
        vpk_add(ctx, BC_G2(th), BC_G(th), BC_G(th), 1);
        vp_cmul(w, BC_W2(th), BC_G(th), BC_A(th), 1);   // scaled curve: a2 = g A
        vp_csqr(w, BC_W4(th), BC_G(th), 1);             //   a4 = g^2
        vpk_add(ctx, BC_W22(th), BC_W2(th), BC_W2(th), 1);
        for (i = 0; i < s->bsQ; i++)
            s2_copy(ctx, VP_COEF(ctx, th->bsa, i), BC_W2(th));
    }

    gettimeofday(&t0, NULL);

    // ---- baby roots: [k]Q (or [f(k)]Q), k in S1, padded to dF, into tree[0] ----
    if (s->bs)
    {
        uint32_t idx;
        s2_bs_baby(s, th, th->tree[0], Q, a24, found, t1, t2);
        for (idx = s->nbaby; idx < s->dF; idx++)
            s2_copy(ctx, VP_COEF(ctx, th->tree[0], idx), VP_COEF(ctx, th->tree[0], 0));
    }
    else
    {
        s2_pt Pm = { PC(13), PC(14) };      // [k-2]Q
        s2_pt Pk = { PC(4), PC(5) };        // [k]Q  (ladder scratch is free here)
        s2_pt P2 = { PC(6), PC(7) };        // [2]Q
        uint64_t *r0 = th->tree[0];
        uint32_t idx = 0;
        c = 0;
        s2_dbl(th, P2, Q, a24);
        s2_copy(ctx, Pk.X, Q.X);            // k = 1
        s2_copy(ctx, Pk.Z, Q.Z);
        s2_copy(ctx, Pm.X, Q.X);            // [-1]Q = [1]Q for the first step
        s2_copy(ctx, Pm.Z, Q.Z);
        for (k = 1; k < s->D / 2; k += 2)
        {
            if ((idx < s->nbaby) && (s->baby[idx] == k))
            {
                s2_copy(ctx, VP_COEF(ctx, th->X, c), Pk.X);
                s2_copy(ctx, VP_COEF(ctx, th->Z, c), Pk.Z);
                idx++;
                if (++c == s->chunk)
                {
                    s2_normalize(s, th, c, VP_COEF(ctx, r0, idx - c), found, t1, t2);
                    c = 0;
                }
            }
            // [k+2] = [k] + [2], difference [k-2] ([-1] = [1] when k = 1)
            {
                uint64_t *nx = PC(8), *nz = PC(9);
                s2_pt Nn = { nx, nz };
                s2_add(th, Nn, Pk, P2, Pm);     // uses PC(0..3)
                s2_copy(ctx, Pm.X, Pk.X);
                s2_copy(ctx, Pm.Z, Pk.Z);
                s2_copy(ctx, Pk.X, nx);
                s2_copy(ctx, Pk.Z, nz);
            }
        }
        if (c > 0)
            s2_normalize(s, th, c, VP_COEF(ctx, r0, idx - c), found, t1, t2);
        // pad with copies of the first root
        for (; idx < s->dF; idx++)
            s2_copy(ctx, VP_COEF(ctx, r0, idx), VP_COEF(ctx, r0, 0));
    }
    vp_from_roots(w, th->F, th->tree, th->tree[0], s->dF);
    th->mf = vp_modf_init(w, th->F, s->dF,
        s->lean ? 0 : VP_MODF_KEEP_INV | VP_MODF_KEEP_F, th->G);
    gettimeofday(&t1v, NULL);
    tb = ytools_difftime(&t0, &t1v);

    if (verbose > 1)
    {
        printf("FFT stage 2 (thread %d): F + 1/F %1.3f s\n", tid, tb);
    }

    // ---- giant roots [jD]Q (or [f(jD)]Q), j = j0, j0 + 1, ... ----
    DQ.X = PC(15); DQ.Z = PC(13);           // PC(13) is free after the baby steps
    if (s->bs)
    {
        // M progressions j = j0 + r + M m, r < M: each step gives M roots
        uint64_t base = s->j0 * s->D;
        uint32_t M = s->bsM, st;
        gettimeofday(&t0, NULL);
        s2_bs_init(s, th, 1, &base, s->D, (uint64_t)M * s->D, Q, a24, found, t1, t2);
        for (b = 0; b < s->blocks; b++)
        {
            if (b > 0)
                gettimeofday(&t0, NULL);
            for (st = 0; st < s->dF / M; st++)
            {
                memcpy(VP_COEF(ctx, th->G, st * M), th->bsx,
                    (size_t)M * 8 * (size_t)ctx->n * sizeof(uint64_t));
                if ((b + 1 < s->blocks) || (st + 1 < s->dF / M))
                    s2_bs_step(s, th, found, t1, t2);
            }
            vp_from_roots(w, th->G, NULL, th->G, s->dF);
            gettimeofday(&t1v, NULL);
            tg += ytools_difftime(&t0, &t1v);
            t0 = t1v;
            if (b == 0)
                vp_mod_monic(w, th->mf, th->H, th->G);      // F may be H
            else
                vp_mulmod_monic(w, th->mf, th->H, th->G);
            gettimeofday(&t1v, NULL);
            tf += ytools_difftime(&t0, &t1v);
        }
    }
    else
    {
        uint64_t *pX = vp_alloc(3, ctx->n), *pZ = vp_alloc(3, ctx->n);
        s2_pt A = { VP_COEF(ctx, pX, 0), VP_COEF(ctx, pZ, 0) };     // [(j-1)D]Q
        s2_pt B = { VP_COEF(ctx, pX, 1), VP_COEF(ctx, pZ, 1) };     // [jD]Q
        s2_pt C = { VP_COEF(ctx, pX, 2), VP_COEF(ctx, pZ, 2) };

        s2_mul(th, DQ, Q, s->D, a24);
        s2_mul(th, A, Q, s->j0 * s->D, a24);                        // j0
        s2_mul(th, B, Q, (s->j0 + 1) * s->D, a24);                  // j0 + 1

        for (b = 0; b < s->blocks; b++)
        {
            gettimeofday(&t0, NULL);
            c = 0;
            for (i = 0; i < s->dF; i++)
            {
                // root i holds [(j0 + b dF + i) D]Q = A; advance A <- B <- B + DQ
                s2_copy(ctx, VP_COEF(ctx, th->X, c), A.X);
                s2_copy(ctx, VP_COEF(ctx, th->Z, c), A.Z);
                s2_add(th, C, B, DQ, A);
                s2_copy(ctx, A.X, B.X); s2_copy(ctx, A.Z, B.Z);
                s2_copy(ctx, B.X, C.X); s2_copy(ctx, B.Z, C.Z);
                if (++c == s->chunk)
                {
                    s2_normalize(s, th, c, VP_COEF(ctx, th->G, i + 1 - c), found, t1, t2);
                    c = 0;
                }
            }
            vp_from_roots(w, th->G, NULL, th->G, s->dF);
            gettimeofday(&t1v, NULL);
            tg += ytools_difftime(&t0, &t1v);
            t0 = t1v;
            if (b == 0)
                vp_mod_monic(w, th->mf, th->H, th->G);      // F may be H
            else
                vp_mulmod_monic(w, th->mf, th->H, th->G);
            gettimeofday(&t1v, NULL);
            tf += ytools_difftime(&t0, &t1v);
        }
        vp_free(pX);
        vp_free(pZ);
    }

    if (verbose > 1)
    {
        printf("FFT stage 2 (thread %d): H*G mod F %1.3f s\n", tid, tf);
    }

    // ---- evaluate H at the baby roots, multiply up ----
    gettimeofday(&t0, NULL);
    vp_eval_tree(w, th->mf, th->tree, th->H, th->H, th->G);
    vp_prod(w, th->G, th->H, s->dF);
    gettimeofday(&t1v, NULL);
    te = ytools_difftime(&t0, &t1v);

    if (verbose > 1)
    {
        printf("FFT stage 2 (thread %d): eval %1.3f s\n", tid, te);
    }

    if (verbose > 1)
    {
        printf("FFT stage 2 (thread %d): giant roots + G %1.3f s\n", tid, tg);
        if (s->bs)
            printf("FFT stage 2 (thread %d): Brent-Suyama table setups %1.3f s\n", tid, th->tinit);
    }

    // lane result = found * prod (as a plain integer), into the backend
    // accumulator (only its gcd with N is used)
    {
        mpz_t Ri;
        mpz_init(Ri);
        mpz_invert(Ri, s->R, s->N);
        for (l = 0; l < 8; l++)
        {
            vp_get(ctx, t1, th->G, l);
            mpz_mul(t1, t1, Ri);
            mpz_mul(t1, t1, found[l]);
            mpz_mod(t1, t1, s->N);
            vec_lane_set(mdata, work->stg2acc, t1, l);
        }
        mpz_clear(Ri);
    }

    // if (verbose > 1)
    // {
    //     printf("FFT stage 2 (thread %d): F + 1/F %1.3f s, giant roots + G %1.3f s, "
    //         "H*G mod F %1.3f s, eval %1.3f s", tid, tb, tg, tf, te);
    //     if (s->bs)
    //         printf(" (Brent-Suyama table setups %1.3f s)", th->tinit);
    //     printf("\n");
    // }

    vp_modf_free(th->mf);
    th->mf = NULL;
    for (l = 0; l < 8; l++) mpz_clear(found[l]);
    mpz_clear(t1);
    mpz_clear(t2);
    mpz_clear(Rbi);
    return 0;
}

#else

avxecm_s2fft_t *avxecm_s2fft_setup(mpz_t N, uint64_t B1, uint64_t B2,
    int threads, double mem_mb, int auto_b2, int bs, int edw_stage1, int verbose)
{
    return NULL;
}
uint64_t avxecm_s2fft_B2(avxecm_s2fft_t *s) { return 0; }
void avxecm_s2fft_free(avxecm_s2fft_t *s) { }
int avxecm_s2fft_run(avxecm_s2fft_t *s, int tid, vec_monty_t *mdata, ecm_work *work,
    ecm_pt *P, int verbose)
{
    return -1;
}

#endif
