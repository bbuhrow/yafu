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
 * the lane result is prod_k H(x([k]Q)).  One curve per lane, polynomial
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
 */

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
    uint64_t *bsx, *bsy, *bsd, *bsi, *bsp;  // Brent-Suyama tables, M (S+1) each
    uint64_t *bc;                           // Brent-Suyama constants, 16
    mpz_t *bsT;                             // the tables' integer multipliers
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
//   c1: stage 1 mulmods per unit of B1; scale: s2_cost vs measured stage 2
typedef struct
{
    double a, b, c;
    double mm, mm_lean, mi, mi_lean, ev, ev_lean;
    double c1, scale;
} s2_model_t;

#if defined(USE_AVX512IFMA) || defined(IFMA)
static const s2_model_t s2m = { 6.819, -36.97, 3.336, 0.817, 0.939, 0.785, 0.566, 1.842, 1.906, 15.5, 0.90 };
#elif defined(USE_AVX512F)
static const s2_model_t s2m = { 3.445, -7.58, 2.687, 0.721, 0.837, 0.762, 0.549, 1.916, 1.870, 14.2, 0.92 };
#else
static const s2_model_t s2m = { 4.219, -14.68, 1.402, 0.669, 0.734, 0.679, 0.492, 1.828, 1.808, 14.5, 1.0 };
#endif

// G's product tree (from dF = 2^t roots), in template mulmods per root
static double s2_tcost(int t, int n)
{
    double c = s2m.a * t + s2m.b + s2m.c * t * t / n;
    return c < 8.0 ? 8.0 : c;
}

static int s2_levels(int t, int k)
{
    return (t + k - 1) / k;                 // s = 0, k, 2k, ... < t
}

// bytes per thread (excluding the shared context)
static double s2_mem(int t, int n, int k, int lean, int bs, uint32_t bsM, double *ctxb)
{
    int K = vp_nprimes(n, t + 1);
    double dF = (double)(1u << t), cb = 64.0 * n;
    double b = (double)vp_work_bytes(n, K, t + 1)
        + (2 + s2_levels(t, k) + (lean ? 1 : 0)) * dF * cb
        + (3.0 * (dF < S2_CHUNK ? dF : S2_CHUNK) + 16.0) * cb
        + (bs ? (5.0 * (bsM ? bsM : 64) * ((bs > 0 ? bs : -bs) + 1) + 16.0) * cb : 0.0)
        + (double)vp_modf_bytes(n, K, 1u << t,
            lean ? 0 : VP_MODF_KEEP_INV | VP_MODF_KEEP_F);
    if (ctxb != NULL)
        *ctxb = (double)vp_ctx_bytes(n, K, t + 1);
    return b;
}

// Brent-Suyama progressions side by side: balances the table setup
// (M (deg+1) ladders over ~deg log2(argument) bits, ~10 mulmods per bit)
// against one batch inversion (8 lanes x ~12 mulmods) per M roots
static uint32_t s2_bs_m(int t, double blocks, int deg)
{
    double dF = (double)(1u << t), D = s2_D[t];
    double roots = dF * (blocks + 1.6);
    double bits = deg * (log((blocks * dF + 1.0) * D) + log(D / 2.0)) / log(2.0);
    double mopt = sqrt(roots * 8.0 * 12.0 / ((deg + 1.0) * 10.0 * bits));
    uint32_t M = 2;
    while ((M * 2 <= 64) && (M * 2 <= (uint32_t)dF) && (M * 1.41 < mopt))
        M *= 2;
    return M;
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
        // Brent-Suyama (fit to yafu timings, design doc section 14): deg
        // affine additions per root; D/6 ~ 1.6 dF baby candidates; tables:
        // 2 M (deg + 1) ladders of ~10 mulmods per bit
        double D = s2_D[t], per;
        if (bsM == 0)
            bsM = s2_bs_m(t, blocks, deg);
        // 7 mulmods + ~6 add/sub (O(n)) per affine addition; one batch
        // inversion (8 lanes x mpz_invert ~ 4 + 160/n mulmods) per M roots
        per = deg * (7.0 + 12.0 / n) + 1.0 + 8.0 * (4.0 + 160.0 / n) / bsM;
        groot = dF * per;
        broot = 1.6 * dF * per;
        init = bsM * (deg + 1.0) * 10.0 * deg *
            (log((blocks * dF + 1.0) * D) + log(D / 2.0)) / log(2.0);
    }
    return blocks * (groot + T * (1.0 + (lean ? s2m.mm_lean : s2m.mm)))
        + init + broot + T * (1.0 + (lean ? s2m.mi_lean : s2m.mi) +
            (lean ? s2m.ev_lean : s2m.ev) + rebuilt);
}

// ---- automatic B2 (design doc section 11.3) ----
// Choose (t, k, lean, blocks) for the least predicted time per factor,
// (T1 + T2) / P(B1, B2), at the factor size B1 is usually chosen for.
// T1 is predicted (stage 2 is set up before stage 1 runs): S2_C1 template
// mulmods per unit of B1 (per ISA, s2m.c1).

#define S2_C1           s2m.c1
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

// probability that a group order of size p/23 (Suyama torsion) is
// B1-smooth but for one prime in (B1, B2]; lnp = ln p
static double s2_prob(const double *r, int len, double lnp, double B1, double B2)
{
    double L = lnp - log(23.0), a = log(B1), b = log(B2), h, sum = 0;
    int i, steps = 100;
    h = (b - a) / steps;
    for (i = 0; i <= steps; i++)
    {
        double l = a + i * h, wt = (i == 0 || i == steps) ? 0.5 : 1.0;
        sum += wt * s2_rho(r, len, (L - l) / a) / l;
    }
    return s2_rho(r, len, L / a) + sum * h;
}

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

static uint64_t s2_covered_B2(int t, uint64_t B1, uint64_t blocks)
{
    uint64_t D = s2_D[t], j0 = (B1 + D / 2) / D;
    if (j0 < 1) j0 = 1;
    return (j0 + blocks * (1u << t) - 1) * D + D / 2 - 1;
}

// returns the chosen B2 (0: none fits) and the configuration
static uint64_t s2_auto_B2(int n, uint64_t B1, int threads, double cap, int bs, uint32_t bsM,
    int *bt, int *bk, int *bl, double *ratio, double *pdig)
{
    int len, t, k, lean;
    double *r = s2_rho_table(&len);
    double T1 = S2_C1 * (double)B1, best = 0;
    double lnp;
    uint64_t bB2 = 0;

    *pdig = s2_target_digits((double)B1);
    lnp = *pdig * log(10.0);
    for (t = 3; t <= 15; t++)
    {
        if (B1 < s2_D[t])
            break;
        for (lean = 0; lean <= 1; lean++)
        {
            for (k = 1; k <= 4; k++)
            {
                double ctxb, mem = s2_mem(t, n, k, lean, bs, bsM, &ctxb);
                uint64_t blocks = 1;
                if (mem + ctxb / threads > cap)
                    continue;
                while (1)
                {
                    double T2 = S2_COST_SCALE * s2_cost(t, n, k, lean, (double)blocks, bs, bsM);
                    uint64_t B2c = s2_covered_B2(t, B1, blocks);
                    double e;
                    if ((T2 > 4.0 * T1) || (B2c > (uint64_t)1e15))
                        break;
                    e = s2_prob(r, len, lnp, (double)B1, (double)B2c) / (T1 + T2);
                    if (e > best)
                    {
                        best = e;
                        bB2 = B2c;
                        *bt = t; *bk = k; *bl = lean;
                        *ratio = T2 / T1;
                    }
                    blocks = blocks < 8 ? blocks + 1 : blocks + blocks / 8;
                }
            }
        }
    }
    free(r);
    return bB2;
}

avxecm_s2fft_t *avxecm_s2fft_setup(mpz_t N, uint64_t B1, uint64_t B2,
    int threads, double mem_mb, int auto_b2, int bs, int verbose)
{
    avxecm_s2fft_t *s;
    int n, t, k, lean, bt = 0, bk = 1, bl = 0;
    uint32_t bsM = 0;                       // 0: s2_bs_m per configuration
    double best = 1e300, bmem = 0, bctx = 0, cap = mem_mb * 1048576.0;
    double aratio = 0, adig = 0;
    uint32_t j;

    if (mpz_even_p(N))
        return NULL;
    if (threads < 1) threads = 1;
    // Brent-Suyama: |S| <= 1 is the plain stage 2; |S| <= 32
    if ((bs >= -1) && (bs <= 1))
        bs = 0;
    if (bs > 32) bs = 32;
    if (bs < -32) bs = -32;
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
        uint64_t aB2 = s2_auto_B2(n, B1, threads, cap, bs, bsM, &bt, &bk, &bl, &aratio, &adig);
        if (aB2 > B1)
        {
            B2 = aB2;
            bmem = s2_mem(bt, n, bk, bl, bs, bsM, &bctx);
        }
        else
        {
            bt = 0;
            auto_b2 = 0;
        }
    }
    if (B2 <= B1)
        return NULL;

    for (t = (bt > 0) ? 16 : 3; t <= 15; t++)       // skipped when auto chose
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
                    bt = t; bk = k; bl = lean;
                    bmem = mem; bctx = ctxb;
                }
            }
        }
    }
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
        s->bsM = s2_bs_m(t, (double)((B2 + s->D / 2 - s->j0 * s->D) / ((uint64_t)s->dF * s->D) + 1),
            bs > 0 ? bs : -bs);
    while (s->bsM > (1u << t)) s->bsM >>= 1;   // M divides dF (both powers of 2)
    if (s->bsM < 2) s->bsM = 2;
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
            printf("ecm: FFT stage 2: Brent-Suyama with %s of degree %d, %u progressions\n",
                s->bs > 0 ? "x^S" : "Dickson D_S(x, -1)", s->bs > 0 ? s->bs : -s->bs, s->bsM);
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

static void s2_thread_free(avxecm_s2fft_t *s, s2_thread_t *th)
{
    int i;
    if (th->w == NULL) return;
    vp_free(th->X); vp_free(th->Z); vp_free(th->pre);
    vp_free(th->G); vp_free(th->H); vp_free(th->pt);
    if (s->lean) vp_free(th->F);
    if (s->bs)
    {
        vp_free(th->bsx); vp_free(th->bsy); vp_free(th->bsd);
        vp_free(th->bsi); vp_free(th->bsp); vp_free(th->bc);
        for (i = 0; i < (int)s->bsM * ((s->bs > 0 ? s->bs : -s->bs) + 1); i++)
            mpz_clear(th->bsT[i]);
        free(th->bsT);
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
        size_t q = (size_t)s->bsM * (size_t)((s->bs > 0 ? s->bs : -s->bs) + 1);
        th->bsx = vp_alloc(q, n);
        th->bsy = vp_alloc(q, n);
        th->bsd = vp_alloc(q, n);
        th->bsi = vp_alloc(q, n);
        th->bsp = vp_alloc(q, n);
        th->bc = vp_alloc(16, n);
        th->bsT = (mpz_t *)malloc(q * sizeof(mpz_t));
        for (i = 0; i < (int)q; i++)
            mpz_init(th->bsT[i]);
    }
    th->tree = (uint64_t **)malloc((size_t)s->t * sizeof(uint64_t *));
    for (i = 0; i < s->t; i++)
        th->tree[i] = (i % s->kstride == 0) ? vp_alloc(s->dF, n) : NULL;
}

// ---- x-only Montgomery curve arithmetic in the vecpoly domain ----
// points are (X, Z) pairs of single coefficients; scratch from th->pt

#define PC(i) VP_COEF(ctx, th->pt, (i))

typedef struct { uint64_t *X, *Z; } s2_pt;

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

// out[i] = 1 / in[i], i < m (Montgomery's batch inversion per lane; pre:
// m coefficients).  A lane whose product is not invertible is folded into
// found[lane] (a stage-2 find when it shares a factor with N).
static void s2_inv_batch(avxecm_s2fft_t *s, s2_thread_t *th, uint32_t m,
    const uint64_t *in, uint64_t *out, uint64_t *pre, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    uint64_t *I = PC(8), *T = PC(9);
    uint32_t i;
    int l;

    s2_copy(ctx, VP_COEF(ctx, pre, 0), VP_COEF(ctx, in, 0));
    for (i = 1; i < m; i++)
        vp_cmul(w, VP_COEF(ctx, pre, i), VP_COEF(ctx, pre, i - 1), VP_COEF(ctx, in, i), 1);
    for (l = 0; l < 8; l++)
    {
        vp_get(ctx, t1, VP_COEF(ctx, pre, m - 1), l);
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
        vp_cmul(w, T, I, VP_COEF(ctx, pre, i - 1), 1);
        vp_cmul(w, I, I, VP_COEF(ctx, in, i), 1);
        s2_copy(ctx, VP_COEF(ctx, out, i), T);
    }
    s2_copy(ctx, VP_COEF(ctx, out, 0), I);
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

// BS curve constants, in th->bc: x0, A, 2A, g, 2g
#define BC_X0(th)  VP_COEF(ctx, (th)->bc, 0)
#define BC_A(th)   VP_COEF(ctx, (th)->bc, 1)
#define BC_A2(th)  VP_COEF(ctx, (th)->bc, 2)
#define BC_G(th)   VP_COEF(ctx, (th)->bc, 3)
#define BC_G2(th)  VP_COEF(ctx, (th)->bc, 4)
#define BC_T(th, i) VP_COEF(ctx, (th)->bc, 5 + (i))     // 11 temporaries

// Tables of the M progressions i_r(m) = i0[r] + m delta: entry (r, e) =
// [Delta^e f(i_r(0))]Q, affine, in th->bsx / th->bsy at r (S+1) + e
static void s2_bs_init(avxecm_s2fft_t *s, s2_thread_t *th, const uint64_t *i0,
    uint64_t delta, s2_pt Q, const uint64_t *a24, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    int S = s->bs, deg = S > 0 ? S : -S, e;
    uint32_t r, M = s->bsM, E = (uint32_t)deg + 1;
    mpz_t c[64];
    s2_pt R0 = { BC_T(th, 0), BC_T(th, 1) }, R1 = { BC_T(th, 2), BC_T(th, 3) };
    uint64_t *u = BC_T(th, 4), *v = BC_T(th, 5), *z11 = BC_T(th, 6), *z12 = BC_T(th, 7);
    uint64_t *t3 = BC_T(th, 8);

    for (e = 0; e <= deg; e++) mpz_init(c[e]);
    for (r = 0; r < M; r++)
    {
        s2_findiff(S, deg, c, i0[r], delta);
        for (e = 0; e <= deg; e++)
        {
            uint32_t q = r * E + (uint32_t)e;
            uint64_t *X1 = R0.X, *Z1 = R0.Z, *X2 = R1.X, *Z2 = R1.Z;
            mpz_set(th->bsT[q], c[e]);
            s2_mul2(th, R0, R1, Q, c[e], a24);
            // Okeya-Sakurai: y = [Z2((x0 X1 + Z1)(x0 Z1 + X1 + 2A Z1) - 2A Z1^2)
            //   - (x0 Z1 - X1)^2 X2] / (2 g Z1^2 Z2),  x = X1 2g Z1 Z2 / (same)
            vp_cmul(w, u, BC_X0(th), X1, 1);
            vpk_add(ctx, u, u, Z1, 1);                      // x0 X1 + Z1
            vp_cmul(w, v, BC_X0(th), Z1, 1);                // x0 Z1
            vp_cmul(w, t3, BC_A2(th), Z1, 1);
            vpk_add(ctx, t3, t3, X1, 1);
            vpk_add(ctx, t3, t3, v, 1);                     // x0 Z1 + X1 + 2A Z1
            vp_cmul(w, u, u, t3, 1);
            vp_cmul(w, z11, Z1, Z1, 1);                     // Z1^2
            vp_cmul(w, t3, BC_A2(th), z11, 1);
            vpk_sub(ctx, u, u, t3, 1);
            vp_cmul(w, u, u, Z2, 1);                        // Z2 (...)
            vpk_sub(ctx, v, v, X1, 1);
            vp_cmul(w, v, v, v, 1);
            vp_cmul(w, v, v, X2, 1);                        // (x0 Z1 - X1)^2 X2
            vpk_sub(ctx, VP_COEF(ctx, th->bsy, q), u, v, 1);
            vp_cmul(w, z12, Z1, Z2, 1);
            vp_cmul(w, z12, z12, BC_G2(th), 1);             // 2g Z1 Z2
            vp_cmul(w, VP_COEF(ctx, th->bsx, q), X1, z12, 1);
            vp_cmul(w, VP_COEF(ctx, th->bsd, q), z12, Z1, 1);   // W = 2g Z1^2 Z2
        }
    }
    for (e = 0; e <= deg; e++) mpz_clear(c[e]);

    s2_inv_batch(s, th, M * E, th->bsd, th->bsi, th->bsp, found, t1, t2);
    vp_cmul(w, th->bsx, th->bsx, th->bsi, M * E);
    vp_cmul(w, th->bsy, th->bsy, th->bsi, M * E);
}

// one step of all M tables: (r, e) <- (r, e) + (r, e+1), e < deg
static void s2_bs_step(avxecm_s2fft_t *s, s2_thread_t *th, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_work_t *w = th->w;
    vp_ctx_t *ctx = w->ctx;
    uint32_t deg = (uint32_t)(s->bs > 0 ? s->bs : -s->bs), e;
    uint32_t r, M = s->bsM, E = deg + 1, q = 0;
    uint64_t *lam = BC_T(th, 0), *x3 = BC_T(th, 1);

    // the multipliers T(r, e) are positive integers (positive coefficients,
    // positive arguments); when T(r, e) = T(r, e+1) the two points are equal
    // and the addition is a doubling (denominator 2 g y)
    for (r = 0; r < M; r++)
        for (e = 0; e < deg; e++, q++)
        {
            if (mpz_cmp(th->bsT[r * E + e], th->bsT[r * E + e + 1]) == 0)
                vp_cmul(w, VP_COEF(ctx, th->bsd, q), BC_G2(th),
                    VP_COEF(ctx, th->bsy, r * E + e), 1);
            else
                vpk_sub(ctx, VP_COEF(ctx, th->bsd, q),
                    VP_COEF(ctx, th->bsx, r * E + e + 1), VP_COEF(ctx, th->bsx, r * E + e), 1);
        }
    s2_inv_batch(s, th, q, th->bsd, th->bsi, th->bsp, found, t1, t2);
    q = 0;
    for (r = 0; r < M; r++)
        for (e = 0; e < deg; e++, q++)
        {
            uint64_t *x1 = VP_COEF(ctx, th->bsx, r * E + e), *y1 = VP_COEF(ctx, th->bsy, r * E + e);
            uint64_t *x2 = VP_COEF(ctx, th->bsx, r * E + e + 1), *y2 = VP_COEF(ctx, th->bsy, r * E + e + 1);
            if (mpz_cmp(th->bsT[r * E + e], th->bsT[r * E + e + 1]) == 0)
            {
                // (3 x1^2 + 2 A x1 + 1) / (2 g y1)
                vpk_add(ctx, lam, x1, x1, 1);
                vpk_add(ctx, lam, lam, x1, 1);
                vpk_add(ctx, lam, lam, BC_A2(th), 1);
                vp_cmul(w, lam, lam, x1, 1);
                vpk_add(ctx, lam, lam, w->one, 1);
            }
            else
                vpk_sub(ctx, lam, y2, y1, 1);                  // (y2 - y1) / (x2 - x1)
            vp_cmul(w, lam, lam, VP_COEF(ctx, th->bsi, q), 1);
            mpz_add(th->bsT[r * E + e], th->bsT[r * E + e], th->bsT[r * E + e + 1]);
            vp_cmul(w, x3, lam, lam, 1);
            vp_cmul(w, x3, x3, BC_G(th), 1);
            vpk_sub(ctx, x3, x3, BC_A(th), 1);
            vpk_sub(ctx, x3, x3, x1, 1);
            vpk_sub(ctx, x3, x3, x2, 1);                       // g lam^2 - A - x1 - x2
            vpk_sub(ctx, x1, x1, x3, 1);
            vp_cmul(w, x1, x1, lam, 1);
            vpk_sub(ctx, y1, x1, y1, 1);                       // lam (x1 - x3) - y1
            s2_copy(ctx, x1, x3);
        }
}

// baby roots x([f(k)]Q), k in S1, into dest (nbaby coefficients); progressions
// k = c + 6 (r + (M/2) m), c in {1, 5}
static void s2_bs_baby(avxecm_s2fft_t *s, s2_thread_t *th, uint64_t *dest, s2_pt Q,
    const uint64_t *a24, mpz_t *found, mpz_t t1, mpz_t t2)
{
    vp_ctx_t *ctx = th->w->ctx;
    uint32_t M = s->bsM, h = M / 2, r, idx = 0, E = (uint32_t)(s->bs > 0 ? s->bs : -s->bs) + 1;
    uint64_t i0[256], m, delta = 6 * (uint64_t)h;
    for (r = 0; r < M; r++)
        i0[r] = ((r & 1) ? 5 : 1) + 6 * (uint64_t)(r >> 1);
    s2_bs_init(s, th, i0, delta, Q, a24, found, t1, t2);
    for (m = 0; ; m++)
    {
        int any = 0;
        for (r = 0; r < M; r++)
        {
            uint64_t k = i0[r] + m * delta;
            if (k < s->D / 2)
            {
                any = 1;
                if (s2_gcd((uint32_t)k, s->D) == 1)
                    s2_copy(ctx, VP_COEF(ctx, dest, idx++), VP_COEF(ctx, th->bsx, r * E));
            }
        }
        if (!any || idx == s->nbaby)
            break;
        s2_bs_step(s, th, found, t1, t2);
    }
}

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

    // ---- giant roots [jD]Q (or [f(jD)]Q), j = j0, j0 + 1, ... ----
    DQ.X = PC(15); DQ.Z = PC(13);           // PC(13) is free after the baby steps
    if (s->bs)
    {
        // M progressions j = j0 + r + M m, r < M: each step gives M roots
        uint64_t i0[256];
        uint32_t M = s->bsM, r, E = (uint32_t)(s->bs > 0 ? s->bs : -s->bs) + 1, st;
        for (r = 0; r < M; r++)
            i0[r] = (s->j0 + r) * s->D;
        gettimeofday(&t0, NULL);
        s2_bs_init(s, th, i0, (uint64_t)M * s->D, Q, a24, found, t1, t2);
        for (b = 0; b < s->blocks; b++)
        {
            if (b > 0)
                gettimeofday(&t0, NULL);
            for (st = 0; st < s->dF / M; st++)
            {
                for (r = 0; r < M; r++)
                    s2_copy(ctx, VP_COEF(ctx, th->G, st * M + r), VP_COEF(ctx, th->bsx, r * E));
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

    // ---- evaluate H at the baby roots, multiply up ----
    gettimeofday(&t0, NULL);
    vp_eval_tree(w, th->mf, th->tree, th->H, th->H, th->G);
    vp_prod(w, th->G, th->H, s->dF);
    gettimeofday(&t1v, NULL);
    te = ytools_difftime(&t0, &t1v);

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

    if (verbose > 1)
    {
        printf("FFT stage 2 (thread %d): F + 1/F %1.3f s, giant roots + G %1.3f s, "
            "H*G mod F %1.3f s, eval %1.3f s\n", tid, tb, tg, tf, te);
    }

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
    int threads, double mem_mb, int auto_b2, int bs, int verbose)
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
