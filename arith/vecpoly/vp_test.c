/*
 * vecpoly tests (every lane against GMP) and timings.
 *   vp_test [bits ...]        default: 200 400 520 1024 2048
 * Options: -l logL (largest test degree 2^logL, default 9), -b (benchmarks
 * only), -s school_max.
 */

#include "vecpoly.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef VP_HAVE_KERNELS

static gmp_randstate_t rs;
static mpz_t gN, gR, gRi;
static int fails = 0;

static double now(void)
{
    struct timespec ts;
#ifdef _WIN32
    timespec_get(&ts, TIME_UTC);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// plain value x (lane) -> Montgomery coefficient
static void set_mont(vp_ctx_t *ctx, uint64_t *c, mpz_t x, int lane)
{
    mpz_t t;
    mpz_init(t);
    mpz_mul(t, x, gR);
    mpz_mod(t, t, gN);
    vp_set(ctx, c, t, lane);
    mpz_clear(t);
}

// coefficient -> plain value mod N (also checks the [0, 2N) invariant)
static void get_plain(vp_ctx_t *ctx, mpz_t x, const uint64_t *c, int lane)
{
    mpz_t t;
    mpz_init(t);
    vp_get(ctx, x, c, lane);
    mpz_mul_2exp(t, gN, 1);
    if (mpz_cmp(x, t) >= 0)
    {
        if (fails < 20) printf("  coefficient out of range [0, 2N)\n");
        fails++;
    }
    mpz_mul(x, x, gRi);
    mpz_mod(x, x, gN);
    mpz_clear(t);
}

// random polynomial: plain values in P[lane][j], Montgomery coefficients in a
static void rand_poly(vp_ctx_t *ctx, uint64_t *a, mpz_t **P, uint32_t len, int extreme)
{
    int l;
    uint32_t j;
    mpz_t t;
    mpz_init(t);
    for (l = 0; l < 8; l++)
    {
        for (j = 0; j < len; j++)
        {
            if (extreme)
            {
                // largest allowed stored value: 2N - 1
                mpz_mul_2exp(t, gN, 1);
                mpz_sub_ui(t, t, 1);
                vp_set(ctx, VP_COEF(ctx, a, j), t, l);
                mpz_mul(P[l][j], t, gRi);
                mpz_mod(P[l][j], P[l][j], gN);
            }
            else
            {
                mpz_urandomm(P[l][j], rs, gN);
                set_mont(ctx, VP_COEF(ctx, a, j), P[l][j], l);
            }
        }
    }
    mpz_clear(t);
}

static mpz_t **plist(uint32_t len)
{
    mpz_t **P = (mpz_t **)malloc(8 * sizeof(mpz_t *));
    int l;
    uint32_t j;
    for (l = 0; l < 8; l++)
    {
        P[l] = (mpz_t *)malloc(len * sizeof(mpz_t));
        for (j = 0; j < len; j++) mpz_init(P[l][j]);
    }
    return P;
}

static void plist_free(mpz_t **P, uint32_t len)
{
    int l;
    uint32_t j;
    for (l = 0; l < 8; l++)
    {
        for (j = 0; j < len; j++) mpz_clear(P[l][j]);
        free(P[l]);
    }
    free(P);
}

static void vp_copy_test(vp_ctx_t *ctx, uint64_t *d, const uint64_t *a, uint32_t cnt)
{
    memcpy(d, a, (size_t)cnt * 8 * (size_t)ctx->n * sizeof(uint64_t));
}

static void check(vp_ctx_t *ctx, const char *what, const uint64_t *a, mpz_t **P,
    uint32_t len)
{
    int l;
    uint32_t j;
    mpz_t x;
    int bad = 0;
    mpz_init(x);
    for (l = 0; l < 8; l++)
        for (j = 0; j < len; j++)
        {
            get_plain(ctx, x, VP_COEF(ctx, a, j), l);
            if (mpz_cmp(x, P[l][j]) != 0)
            {
                if (bad < 3)
                    gmp_printf("  %s: lane %d coef %u: got %Zx want %Zx\n",
                        what, l, j, x, P[l][j]);
                bad++;
            }
        }
    if (bad) fails++;
    mpz_clear(x);
}

// plain reference: C = A * B (full)
static void ref_mul(mpz_t *C, mpz_t *A, uint32_t la, mpz_t *B, uint32_t lb)
{
    uint32_t i, j;
    for (i = 0; i < la + lb - 1; i++) mpz_set_ui(C[i], 0);
    for (i = 0; i < la; i++)
        for (j = 0; j < lb; j++)
            mpz_addmul(C[i + j], A[i], B[j]);
    for (i = 0; i < la + lb - 1; i++) mpz_mod(C[i], C[i], gN);
}

static void test_mul(vp_ctx_t *ctx, vp_work_t *w, uint32_t la, uint32_t lb, int extreme)
{
    uint64_t *a = vp_alloc(la, ctx->n), *b = vp_alloc(lb, ctx->n);
    uint64_t *c = vp_alloc(la + lb, ctx->n);
    mpz_t **A = plist(la), **B = plist(lb), **C = plist(la + lb);
    char what[64];
    int l;
    rand_poly(ctx, a, A, la, extreme);
    rand_poly(ctx, b, B, lb, extreme);
    vp_mul(w, c, a, la, b, lb);
    for (l = 0; l < 8; l++) ref_mul(C[l], A[l], la, B[l], lb);
    sprintf(what, "mul %ux%u%s", la, lb, extreme ? " (max)" : "");
    check(ctx, what, c, C, la + lb - 1);
    vp_free(a); vp_free(b); vp_free(c);
    plist_free(A, la); plist_free(B, lb); plist_free(C, la + lb);
}

// lean: only tree levels 0, 3, 6, ... kept, no stored transforms in the vp_modf_t
static void test_tree(vp_ctx_t *ctx, vp_work_t *w, uint32_t m, int lean)
{
    uint64_t *r = vp_alloc(m, ctx->n), *f = vp_alloc(m, ctx->n);
    uint64_t *h = vp_alloc(m, ctx->n), *g = vp_alloc(m, ctx->n);
    uint64_t *vals = vp_alloc(m, ctx->n), *scr = vp_alloc(m, ctx->n);
    uint64_t **tree;
    mpz_t **Rt = plist(m), **F = plist(m + 1), **H = plist(m), **G = plist(m + 1);
    mpz_t **T = plist(2 * m + 1), **V = plist(m);
    vp_modf_t *mf;
    int l, s, t = 0;
    uint32_t i, j;
    char what[64];
    mpz_t x;
    mpz_init(x);

    while ((1u << t) < m) t++;
    tree = (uint64_t **)malloc((t + 1) * sizeof(uint64_t *));
    for (s = 0; s < t; s++) tree[s] = (lean && (s % 3 != 0)) ? NULL : vp_alloc(m, ctx->n);

    rand_poly(ctx, r, Rt, m, 0);
    vp_from_roots(w, f, tree, r, m);
    // reference: F = prod (X - r), low coefficient first, F[l][m] = 1
    for (l = 0; l < 8; l++)
    {
        mpz_set_ui(F[l][0], 1);
        for (i = 1; i <= m; i++) mpz_set_ui(F[l][i], 0);
        for (i = 0; i < m; i++)
        {
            // F = F * (X - r_i): degree i -> i+1
            for (j = i + 1; j > 0; j--)
            {
                mpz_mul(x, F[l][j], Rt[l][i]);
                mpz_sub(F[l][j], F[l][j - 1], x);
                mpz_mod(F[l][j], F[l][j], gN);
            }
            mpz_mul(F[l][0], F[l][0], Rt[l][i]);
            mpz_neg(F[l][0], F[l][0]);
            mpz_mod(F[l][0], F[l][0], gN);
        }
        for (i = 0; i < m; i++) mpz_set(V[l][i], F[l][i]);
    }
    sprintf(what, "from_roots %u%s", m, lean ? " (lean)" : "");
    check(ctx, what, f, V, m);

    // the same without a tree, in place
    vp_copy_test(ctx, g, r, m);
    vp_from_roots(w, g, NULL, g, m);
    sprintf(what, "from_roots %u (no tree)", m);
    check(ctx, what, g, V, m);

    mf = vp_modf_init(w, f, m, lean ? 0 : VP_MODF_KEEP_INV | VP_MODF_KEEP_F,
        lean ? NULL : scr);
    // reciprocal: inv * F = X^(2m-1) + (degree < m)
    for (l = 0; l < 8; l++)
    {
        for (i = 0; i < m; i++)
            get_plain(ctx, V[l][i], VP_COEF(ctx, mf->inv, i), l);
        ref_mul(T[l], V[l], m, F[l], m + 1);
        for (i = m; i < 2 * m - 1; i++)
            if (mpz_sgn(T[l][i]) != 0) { fails++; printf("  recip %u: lane %d coef %u nonzero\n", m, l, i); break; }
        if (mpz_cmp_ui(T[l][2 * m - 1], 1) != 0) { fails++; printf("  recip %u: lane %d top != 1\n", m, l); }
    }

    // h = random, g = random monic; h * G mod F
    rand_poly(ctx, h, H, m, 0);
    rand_poly(ctx, g, G, m, 0);
    for (l = 0; l < 8; l++)
    {
        mpz_set_ui(G[l][m], 1);
        ref_mul(T[l], H[l], m, G[l], m + 1);     // degree <= 2m - 1
        // reduce mod monic F (low first): for k = 2m-1 down to m
        for (i = 2 * m - 1; i >= m; i--)
        {
            mpz_set(x, T[l][i]);
            for (j = 0; j <= m; j++)
                mpz_submul(T[l][i - m + j], x, F[l][j]);
            for (j = 0; j <= m; j++)
                mpz_mod(T[l][i - m + j], T[l][i - m + j], gN);
        }
    }
    vp_mulmod_monic(w, mf, h, g);
    sprintf(what, "mulmod %u%s", m, lean ? " (lean)" : "");
    check(ctx, what, h, T, m);

    // evaluation at the roots: H = the reduced product
    for (l = 0; l < 8; l++)
        for (i = 0; i < m; i++)
        {
            // Horner
            mpz_set_ui(x, 0);
            for (j = m; j > 0; j--)
            {
                mpz_mul(x, x, Rt[l][i]);
                mpz_add(x, x, T[l][j - 1]);
                mpz_mod(x, x, gN);
            }
            mpz_set(V[l][i], x);
        }
    vp_eval_tree(w, mf, tree, vals, h, scr);
    sprintf(what, "eval %u%s", m, lean ? " (lean)" : "");
    check(ctx, what, vals, V, m);
    // in place
    vp_eval_tree(w, mf, tree, h, h, scr);
    sprintf(what, "eval %u in place", m);
    check(ctx, what, h, V, m);

    vp_modf_free(mf);
    for (s = 0; s < t; s++) if (tree[s] != NULL) vp_free(tree[s]);
    free(tree);
    vp_free(r); vp_free(f); vp_free(h); vp_free(g); vp_free(vals); vp_free(scr);
    plist_free(Rt, m); plist_free(F, m + 1); plist_free(H, m);
    plist_free(G, m + 1); plist_free(T, 2 * m + 1); plist_free(V, m);
    mpz_clear(x);
}

static void bench(vp_ctx_t *ctx, vp_work_t *w, uint32_t m, int lean)
{
    uint64_t *r = vp_alloc(m, ctx->n), *f = vp_alloc(m, ctx->n), *scr = vp_alloc(m, ctx->n);
    uint64_t *h = vp_alloc(m, ctx->n), *g = vp_alloc(m, ctx->n);
    uint64_t *vals = vp_alloc(m, ctx->n);
    uint64_t **tree;
    mpz_t **Rt = plist(m);
    vp_modf_t *mf;
    double t0, tr, tt, tm, te, tn, tc;
    int s, t = 0, it, reps;
    vec_bignum_t A, C;

    while ((1u << t) < m) t++;
    tree = (uint64_t **)malloc((t + 1) * sizeof(uint64_t *));
    for (s = 0; s < t; s++) tree[s] = (lean && (s % 3 != 0)) ? NULL : vp_alloc(m, ctx->n);
    rand_poly(ctx, r, Rt, m, 0);
    rand_poly(ctx, h, Rt, m, 0);

    reps = m >= 1024 ? 1 : 4096 / m;
    vp_from_roots(w, f, tree, r, m);    // warm up
    t0 = now();
    for (it = 0; it < reps; it++) vp_from_roots(w, f, tree, r, m);
    tt = (now() - t0) / reps;
    t0 = now();
    for (it = 0; it < reps; it++) vp_from_roots(w, g, NULL, r, m);
    tr = (now() - t0) / reps;
    t0 = now();
    mf = vp_modf_init(w, f, m, lean ? 0 : VP_MODF_KEEP_INV | VP_MODF_KEEP_F, scr);
    tn = now() - t0;
    t0 = now();
    for (it = 0; it < reps; it++) vp_mulmod_monic(w, mf, h, g);
    tm = (now() - t0) / reps;
    t0 = now();
    for (it = 0; it < reps; it++) vp_eval_tree(w, mf, tree, vals, h, scr);
    te = (now() - t0) / reps;
    // one template mulmod, for scale
    A.data = (base_t *)h; A.size = ctx->n;
    C.data = (base_t *)vals; C.size = ctx->n;
    t0 = now();
    for (it = 0; it < 200000; it++) w->mul(&A, &A, &C, &w->vn, &w->vs, &w->md);
    tc = (now() - t0) / 200000;
    printf("  %s m=%5u: tree %8.2f ms, roots->poly %8.2f ms, modf_init %8.2f ms, "
        "mulmod %8.2f ms, eval %8.2f ms   (mulmod = %.0f ns; tree = %.1f mulmods/root)\n",
        lean ? "lean" : "full", m, tt * 1e3, tr * 1e3, tn * 1e3, tm * 1e3, te * 1e3, tc * 1e9,
        tr / tc / m);

    vp_modf_free(mf);
    for (s = 0; s < t; s++) if (tree[s] != NULL) vp_free(tree[s]);
    free(tree);
    vp_free(r); vp_free(f); vp_free(h); vp_free(g); vp_free(vals); vp_free(scr);
    plist_free(Rt, m);
}

int main(int argc, char **argv)
{
    int bitsl[16] = { 200, 400, 520, 1024, 2048 }, nb = 5, i, logL = 9;
    int benchonly = 0, school = -1;
    uint32_t m;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-l") && i + 1 < argc) logL = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-b")) benchonly = 1;
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) school = atoi(argv[++i]);
        else
        {
            if (nb == 5 && bitsl[0] == 200) nb = 0;
            bitsl[nb++] = atoi(argv[i]);
        }
    }

    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 12345);
    mpz_init(gN); mpz_init(gR); mpz_init(gRi);

    for (i = 0; i < nb; i++)
    {
        int bits = bitsl[i], n = (bits + 3 + 51) / 52;
        vp_ctx_t *ctx;
        vp_work_t *w;

        mpz_urandomb(gN, rs, bits);
        mpz_setbit(gN, bits - 1);
        mpz_setbit(gN, 0);
        mpz_set_ui(gR, 1);
        mpz_mul_2exp(gR, gR, 52 * n);
        mpz_invert(gRi, gR, gN);

        ctx = vp_ctx_init(gN, n, logL + 1 > 9 ? logL + 1 : 9);
        if (ctx == NULL) { printf("ctx init failed at %d bits\n", bits); return 1; }
        if (school >= 0) ctx->school_max = school;
        w = vp_work_init(ctx);
        printf("%d bits: n = %d words, K = %d primes, Lmax = %u, school_max = %d\n",
            bits, n, ctx->K, ctx->Lmax, ctx->school_max);

        if (!benchonly)
        {
            int sm = ctx->school_max;
            uint32_t sizes[] = { 1, 2, 3, 5, 8, 9, 16, 17, 31, 64, 100, 128 };
            unsigned si;
            int f0 = fails;
            for (si = 0; si < sizeof(sizes) / sizeof(sizes[0]); si++)
            {
                test_mul(ctx, w, sizes[si], sizes[si], 0);
                test_mul(ctx, w, sizes[si], sizes[si] + 3, 0);
            }
            ctx->school_max = 0;     // force the NTT path
            for (si = 0; si < sizeof(sizes) / sizeof(sizes[0]); si++)
                test_mul(ctx, w, sizes[si], sizes[si], 0);
            test_mul(ctx, w, 64, 64, 1);
            test_mul(ctx, w, 1u << logL, 1u << logL, 1);
            ctx->school_max = sm;
            for (m = 1; m <= (1u << logL); m <<= 1)
            {
                test_tree(ctx, w, m, 0);
                test_tree(ctx, w, m, 1);
            }
            ctx->school_max = 0;
            test_tree(ctx, w, 64, 0);
            test_tree(ctx, w, 64, 1);
            ctx->school_max = sm;
            printf("  tests: %s\n", fails == f0 ? "pass" : "FAIL");
        }
        for (m = 64; m <= (1u << logL); m <<= 1)
        {
            bench(ctx, w, m, 0);
            bench(ctx, w, m, 1);
        }
        printf("  memory: ctx %.1f MB, work %.1f MB (K = %d, logLmax = %d)\n",
            vp_ctx_bytes(n, ctx->K, ctx->logLmax) / 1048576.0,
            vp_work_bytes(n, ctx->K, ctx->logLmax) / 1048576.0, ctx->K, ctx->logLmax);

        vp_work_free(w);
        vp_ctx_free(ctx);
    }
    printf("%s\n", fails ? "FAILED" : "all passed");
    return fails ? 1 : 0;
}

#else
int main(void) { printf("vecpoly kernels not built for this target\n"); return 0; }
#endif
