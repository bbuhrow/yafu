/*
 * AVX2 Montgomery experiments: bit-exact tests and benchmarks (handoff 4.3).
 *
 *   a2m_bench [-t] [-b] [-s seed] [-v variant]
 *     -t  correctness only   -b  benchmarks only (default: both)
 *     -v  0 = FP52, 1 = D28, 2 = D32 (default: all)
 *
 * Correctness: every lane, bit for bit, against a GMP model of the
 * almost-Montgomery product (U = (a*b + m*N) / R, m = -a*b/N mod R,
 * U -= N if U >= R): random inputs in [0, R), worst case R - 1, aliasing,
 * long dependent chains, and powm against mpz_powm.
 * Benchmarks: mul / sqr per number, powm (full-size exponent) per number
 * vs mpz_powm, 256 .. 4096 bits.
 */

#include "a2m.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>

static const char *vname[A2M_NVARIANT] = { "FP52", "D28", "D32" };
static const int vdbits[A2M_NVARIANT] = { 52, 28, 32 };

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static uint64_t *a2m_valloc(int S)
{
    void *p;
    if (posix_memalign(&p, 64, (size_t)S * A2M_VL * sizeof(uint64_t) + 64) != 0)
    {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    memset(p, 0, (size_t)S * A2M_VL * sizeof(uint64_t) + 64);
    return (uint64_t *)p;
}

typedef struct
{
    a2m_ctx_t c;
    a2m_mul_fn mul;
    a2m_sqr_fn sqr;
    mpz_t N, R, Ninv;           // Ninv = N^-1 mod R
} bench_t;

static void set_lane(uint64_t *v, int lane, mpz_t x, int S, int D)
{
    mpz_t t;
    int j;
    mpz_init_set(t, x);
    for (j = 0; j < S; j++)
    {
        v[j * A2M_VL + lane] = mpz_get_ui(t) & ((D == 64) ? ~0ULL : ((1ULL << D) - 1));
        mpz_tdiv_q_2exp(t, t, (mp_bitcnt_t)D);
    }
    mpz_clear(t);
}

static void get_lane(mpz_t x, uint64_t *v, int lane, int S, int D)
{
    int j;
    mpz_set_ui(x, 0);
    for (j = S - 1; j >= 0; j--)
    {
        mpz_mul_2exp(x, x, (mp_bitcnt_t)D);
        mpz_add_ui(x, x, (unsigned long)v[j * A2M_VL + lane]);
    }
}

static int bench_init(bench_t *b, int variant, mpz_t N)
{
    int D = vdbits[variant];
    int bits = (int)mpz_sizeinbase(N, 2);
    int S = (bits + 3 + D - 1) / D, l;
    mpz_t t;

    if (S < 2)
        S = 2;
    if (!a2m_get(variant, S, &b->mul, &b->sqr))
        return 0;
    b->c.variant = variant;
    b->c.dbits = D;
    b->c.S = S;
    b->c.N = a2m_valloc(S);
    b->c.M = a2m_valloc(S);
    b->c.W = a2m_valloc(3 * S);
    mpz_init_set(b->N, N);
    mpz_init(b->R);
    mpz_init(b->Ninv);
    mpz_init(t);
    mpz_setbit(b->R, (mp_bitcnt_t)(D * S));
    mpz_invert(b->Ninv, N, b->R);
    // rho = -1/N mod 2^D
    mpz_set_ui(t, 1);
    mpz_mul_2exp(t, t, (mp_bitcnt_t)D);
    mpz_invert(t, N, t);
    mpz_neg(t, t);
    mpz_fdiv_r_2exp(t, t, (mp_bitcnt_t)D);
    for (l = 0; l < A2M_VL; l++)
    {
        set_lane(b->c.N, l, N, S, D);
        b->c.rho[l] = mpz_get_ui(t);
    }
    mpz_clear(t);
    return 1;
}

static void bench_free(bench_t *b)
{
    free(b->c.N);
    free(b->c.M);
    free(b->c.W);
    mpz_clear(b->N);
    mpz_clear(b->R);
    mpz_clear(b->Ninv);
}

// GMP model of the AMM product
static void amm_model(bench_t *b, mpz_t u, mpz_t x, mpz_t y)
{
    mpz_t t, m;
    mpz_init(t);
    mpz_init(m);
    mpz_mul(t, x, y);
    mpz_mul(m, t, b->Ninv);
    mpz_neg(m, m);
    mpz_mod(m, m, b->R);
    mpz_addmul(t, m, b->N);
    mpz_tdiv_q_2exp(u, t, (mp_bitcnt_t)(b->c.dbits * b->c.S));
    if (mpz_cmp(u, b->R) >= 0)
        mpz_sub(u, u, b->N);
    mpz_clear(t);
    mpz_clear(m);
}

static int check_vec(bench_t *b, uint64_t *v, mpz_t *expect, const char *what)
{
    mpz_t g;
    int l, bad = 0;
    mpz_init(g);
    for (l = 0; l < A2M_VL; l++)
    {
        get_lane(g, v, l, b->c.S, b->c.dbits);
        if (mpz_cmp(g, expect[l]) != 0)
        {
            if (bad == 0)
                gmp_printf("  %s lane %d: got %Zx\n  expected %Zx\n", what, l, g, expect[l]);
            bad++;
        }
    }
    mpz_clear(g);
    return bad;
}

// left-to-right sliding window powm, all lanes share e; d, b in Montgomery form
static void vpowm(bench_t *b, uint64_t *d, uint64_t *base, mpz_t e, uint64_t *one_m)
{
    int S = b->c.S;
    int bit = (int)mpz_sizeinbase(e, 2) - 1;
    int k = (bit > 1000) ? 6 : (bit > 200) ? 5 : 4;
    int ntab = 1 << (k - 1), i, j;
    uint64_t **g = (uint64_t **)malloc((size_t)ntab * sizeof(uint64_t *));
    uint64_t *b2 = a2m_valloc(S);
    for (i = 0; i < ntab; i++)
        g[i] = a2m_valloc(S);
    memcpy(g[0], base, (size_t)S * A2M_VL * sizeof(uint64_t));
    b->sqr(&b->c, b2, base);
    for (i = 1; i < ntab; i++)
        b->mul(&b->c, g[i], g[i - 1], b2);
    memcpy(d, one_m, (size_t)S * A2M_VL * sizeof(uint64_t));
    while (bit >= 0)
    {
        if (mpz_tstbit(e, (mp_bitcnt_t)bit) == 0)
        {
            b->sqr(&b->c, d, d);
            bit--;
        }
        else
        {
            int a = bit - k + 1, thisk, bstr = 1;
            if (a < 0) a = 0;
            a = (int)mpz_scan1(e, (mp_bitcnt_t)a);
            thisk = bit - a + 1;
            for (j = 1; j < thisk; j++)
                bstr = (bstr << 1) | mpz_tstbit(e, (mp_bitcnt_t)(bit - j));
            for (j = 0; j < thisk; j++)
                b->sqr(&b->c, d, d);
            b->mul(&b->c, d, d, g[bstr >> 1]);
            bit -= thisk;
        }
    }
    for (i = 0; i < ntab; i++)
        free(g[i]);
    free(g);
    free(b2);
}

static int test_size(int variant, int bits, gmp_randstate_t rs)
{
    bench_t b;
    mpz_t N, x[4], y[4], u[4], t, e;
    uint64_t *va, *vb, *vc, *one_m;
    int fails = 0, l, it, S;

    mpz_init(N); mpz_init(t); mpz_init(e);
    for (l = 0; l < 4; l++) { mpz_init(x[l]); mpz_init(y[l]); mpz_init(u[l]); }
    mpz_urandomb(N, rs, (mp_bitcnt_t)bits);
    mpz_setbit(N, (mp_bitcnt_t)(bits - 1));
    mpz_setbit(N, 0);
    if (!bench_init(&b, variant, N))
    {
        printf("  %s %d bits: unsupported\n", vname[variant], bits);
        return 0;
    }
    S = b.c.S;
    va = a2m_valloc(S); vb = a2m_valloc(S); vc = a2m_valloc(S); one_m = a2m_valloc(S);

    // random inputs in [0, R), mul and sqr, aliasing
    for (it = 0; it < 200; it++)
    {
        for (l = 0; l < 4; l++)
        {
            mpz_urandomm(x[l], rs, b.R);
            mpz_urandomm(y[l], rs, b.R);
            if (it == 0) { mpz_sub_ui(x[l], b.R, 1); mpz_sub_ui(y[l], b.R, 1); }
            set_lane(va, l, x[l], S, b.c.dbits);
            set_lane(vb, l, y[l], S, b.c.dbits);
            amm_model(&b, u[l], x[l], y[l]);
        }
        b.mul(&b.c, vc, va, vb);
        fails += check_vec(&b, vc, u, "mul");
        b.mul(&b.c, va, va, vb);                // c == a
        fails += check_vec(&b, va, u, "mul c==a");
        for (l = 0; l < 4; l++)
        {
            set_lane(va, l, x[l], S, b.c.dbits);
            amm_model(&b, u[l], x[l], x[l]);
        }
        b.sqr(&b.c, vc, va);
        fails += check_vec(&b, vc, u, "sqr");
        b.sqr(&b.c, va, va);
        fails += check_vec(&b, va, u, "sqr c==a");
    }

    // dependent chain of squarings and multiplies
    for (l = 0; l < 4; l++)
    {
        mpz_urandomm(x[l], rs, b.R);
        mpz_urandomm(y[l], rs, b.R);
        set_lane(va, l, x[l], S, b.c.dbits);
        set_lane(vb, l, y[l], S, b.c.dbits);
    }
    for (it = 0; it < 500; it++)
    {
        b.sqr(&b.c, va, va);
        b.mul(&b.c, va, va, vb);
        for (l = 0; l < 4; l++)
        {
            amm_model(&b, x[l], x[l], x[l]);
            amm_model(&b, x[l], x[l], y[l]);
        }
    }
    fails += check_vec(&b, va, x, "chain");

    // powm vs mpz_powm: base in Montgomery form, result out of it
    mpz_urandomb(e, rs, (mp_bitcnt_t)bits);
    mpz_setbit(e, (mp_bitcnt_t)(bits - 1));
    for (l = 0; l < 4; l++)
    {
        mpz_urandomm(x[l], rs, N);
        mpz_mul_2exp(t, x[l], (mp_bitcnt_t)(b.c.dbits * S));
        mpz_mod(t, t, N);
        set_lane(va, l, t, S, b.c.dbits);
        mpz_mul_2exp(t, t, 0);
        mpz_set_ui(t, 1);
        mpz_mul_2exp(t, t, (mp_bitcnt_t)(b.c.dbits * S));
        mpz_mod(t, t, N);
        set_lane(one_m, l, t, S, b.c.dbits);
        mpz_powm(u[l], x[l], e, N);
    }
    vpowm(&b, vc, va, e, one_m);
    {
        mpz_t g, Rinv;
        mpz_init(g);
        mpz_init(Rinv);
        mpz_invert(Rinv, b.R, N);
        for (l = 0; l < 4; l++)
        {
            get_lane(g, vc, l, S, b.c.dbits);
            mpz_mul(g, g, Rinv);
            mpz_mod(g, g, N);
            if (mpz_cmp(g, u[l]) != 0)
            {
                printf("  powm lane %d wrong\n", l);
                fails++;
            }
        }
        mpz_clear(g);
        mpz_clear(Rinv);
    }

    printf("  %-4s %4d bits (S = %3d): %s\n", vname[variant], bits, S, fails ? "FAILED" : "ok");
    free(va); free(vb); free(vc); free(one_m);
    bench_free(&b);
    mpz_clear(N); mpz_clear(t); mpz_clear(e);
    for (l = 0; l < 4; l++) { mpz_clear(x[l]); mpz_clear(y[l]); mpz_clear(u[l]); }
    return fails;
}

static void warmup(double sec)
{
    double t0 = now_sec();
    volatile double x = 1.0;
    while (now_sec() - t0 < sec)
        x = x * 1.0000001 + 1e-9;
}

// per-number times: mul ns, sqr ns, powm us
static void bench_size(int variant, int bits, gmp_randstate_t rs,
    double *tmul, double *tsqr, double *tpowm, int *S_out)
{
    bench_t b;
    mpz_t N, x, e;
    uint64_t *va, *vb, *one_m;
    int S, r, l;
    long it, iters;
    double best;

    mpz_init(N); mpz_init(x); mpz_init(e);
    mpz_urandomb(N, rs, (mp_bitcnt_t)bits);
    mpz_setbit(N, (mp_bitcnt_t)(bits - 1));
    mpz_setbit(N, 0);
    *tmul = *tsqr = *tpowm = 0;
    if (!bench_init(&b, variant, N))
        return;
    S = *S_out = b.c.S;
    va = a2m_valloc(S); vb = a2m_valloc(S); one_m = a2m_valloc(S);
    for (l = 0; l < 4; l++)
    {
        mpz_urandomm(x, rs, N);
        set_lane(va, l, x, S, b.c.dbits);
        mpz_urandomm(x, rs, N);
        set_lane(vb, l, x, S, b.c.dbits);
        mpz_set_ui(x, 1);
        mpz_mul_2exp(x, x, (mp_bitcnt_t)(b.c.dbits * S));
        mpz_mod(x, x, N);
        set_lane(one_m, l, x, S, b.c.dbits);
    }
    iters = (long)(4e8 / ((double)S * S * 2 + 50));
    if (iters < 50) iters = 50;

    warmup(0.2);
    for (best = 1e30, r = 0; r < 3; r++)
    {
        double t0 = now_sec();
        for (it = 0; it < iters; it++)
            b.mul(&b.c, va, va, vb);
        t0 = (now_sec() - t0) / (double)iters;
        if (t0 < best) best = t0;
    }
    *tmul = best * 1e9 / A2M_VL;
    for (best = 1e30, r = 0; r < 3; r++)
    {
        double t0 = now_sec();
        for (it = 0; it < iters; it++)
            b.sqr(&b.c, va, va);
        t0 = (now_sec() - t0) / (double)iters;
        if (t0 < best) best = t0;
    }
    *tsqr = best * 1e9 / A2M_VL;

    mpz_urandomb(e, rs, (mp_bitcnt_t)bits);
    mpz_setbit(e, (mp_bitcnt_t)(bits - 1));
    for (best = 1e30, r = 0; r < 3; r++)
    {
        int reps = (bits <= 1024) ? 20 : (bits <= 2048) ? 4 : 1, q;
        double t0 = now_sec();
        for (q = 0; q < reps; q++)
            vpowm(&b, vb, va, e, one_m);
        t0 = (now_sec() - t0) / reps;
        if (t0 < best) best = t0;
    }
    *tpowm = best * 1e6 / A2M_VL;

    free(va); free(vb); free(one_m);
    bench_free(&b);
    mpz_clear(N); mpz_clear(x); mpz_clear(e);
}

static double gmp_powm_us(int bits, gmp_randstate_t rs)
{
    mpz_t N, x, e, d;
    int r, q, reps = (bits <= 1024) ? 20 : (bits <= 2048) ? 4 : 1;
    double best = 1e30;
    mpz_init(N); mpz_init(x); mpz_init(e); mpz_init(d);
    mpz_urandomb(N, rs, (mp_bitcnt_t)bits);
    mpz_setbit(N, (mp_bitcnt_t)(bits - 1));
    mpz_setbit(N, 0);
    mpz_urandomm(x, rs, N);
    mpz_urandomb(e, rs, (mp_bitcnt_t)bits);
    mpz_setbit(e, (mp_bitcnt_t)(bits - 1));
    warmup(0.2);
    for (r = 0; r < 3; r++)
    {
        double t0 = now_sec();
        for (q = 0; q < reps; q++)
            mpz_powm(d, x, e, N);
        t0 = (now_sec() - t0) / reps;
        if (t0 < best) best = t0;
    }
    mpz_clear(N); mpz_clear(x); mpz_clear(e); mpz_clear(d);
    return best * 1e6;
}

int main(int argc, char **argv)
{
    static const int sizes[] = { 256, 384, 512, 640, 768, 1024, 1280, 1536,
        1792, 2048, 2560, 3072, 3584, 4096 };
    const int nsizes = (int)(sizeof(sizes) / sizeof(sizes[0]));
    int do_t = 1, do_b = 1, vsel = -1, i, v, fails = 0;
    unsigned long seed = (unsigned long)time(NULL);
    gmp_randstate_t rs;

    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-t") == 0) do_b = 0;
        else if (strcmp(argv[i], "-b") == 0) do_t = 0;
        else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) seed = strtoul(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "-v") == 0 && i + 1 < argc) vsel = atoi(argv[++i]);
    }
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, seed);
#ifdef A2M_CALLER_RZ
    // experiment: RZ for the whole run instead of per kernel call
    _mm_setcsr((_mm_getcsr() & ~0x6000u) | 0x6000u);
    printf("MXCSR: round toward zero for the whole run (A2M_CALLER_RZ)\n");
#endif
    printf("seed %lu\n", seed);

    if (do_t)
    {
        static const int tsizes[] = { 100, 200, 256, 300, 512, 700, 1000, 1024,
            1500, 2048, 3000, 4096, 5000 };
        printf("correctness (bit-exact AMM vs GMP, all lanes, R-1, aliasing, chains, powm)\n");
        for (v = 0; v < A2M_NVARIANT; v++)
        {
            if (vsel >= 0 && v != vsel) continue;
            for (i = 0; i < (int)(sizeof(tsizes) / sizeof(tsizes[0])); i++)
                fails += test_size(v, tsizes[i], rs);
        }
        printf("%s: %d failures\n", fails ? "TESTS FAILED" : "TESTS PASSED", fails);
    }

    if (do_b)
    {
        double tm[A2M_NVARIANT][16], ts[A2M_NVARIANT][16], tp[A2M_NVARIANT][16], tg[16];
        int Ss[A2M_NVARIANT][16];
        printf("\nper-number times (4 numbers per call); mul/sqr ns, powm us (full-size exponent)\n");
        printf(" bits |");
        for (v = 0; v < A2M_NVARIANT; v++)
            if (vsel < 0 || v == vsel) printf(" %-4s  S  mul    sqr     powm |", vname[v]);
        printf("  GMP powm | best/GMP\n");
        for (i = 0; i < nsizes; i++)
        {
            double best = 1e30;
            printf(" %4d |", sizes[i]);
            for (v = 0; v < A2M_NVARIANT; v++)
            {
                if (vsel >= 0 && v != vsel) continue;
                bench_size(v, sizes[i], rs, &tm[v][i], &ts[v][i], &tp[v][i], &Ss[v][i]);
                printf("  %3d %6.1f %6.1f %8.2f |", Ss[v][i], tm[v][i], ts[v][i], tp[v][i]);
                if (tp[v][i] > 0 && tp[v][i] < best) best = tp[v][i];
                fflush(stdout);
            }
            tg[i] = gmp_powm_us(sizes[i], rs);
            printf("  %8.2f | %6.2fx\n", tg[i], tg[i] / best);
            fflush(stdout);
        }
    }
    gmp_randclear(rs);
    return fails ? 1 : 0;
}
