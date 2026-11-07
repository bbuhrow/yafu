/*
 * AVX-ECM vector Montgomery arithmetic: correctness tests and benchmark.
 *
 *   ae_bench [-s seed] [-t] [-b] [-p] [-w S1,S2,...] [-spare k]
 *     -t  run correctness tests only     -b  benchmarks only (default both)
 *     -p  include powm timing            -w  word counts to run
 *     -neg  negative self-test: a corrupting kernel must FAIL
 *     -a  add/sub-mod range test on AMM outputs
 *     -spare k  N has 52*S - k bits (default: k = 1 and a mid-size N)
 *
 * Kernels compared per word count S:
 *   old : what yafu uses today (fixed416/624/832/1040 when S = 8/12/16/20,
 *         else the BFIPS vecmulmod52/vecsqrmod52 when S % 4 == 0)
 *   new : vecarith52_tmpl.c constant-S template (any S)
 *
 * Correctness is checked bit for bit against a GMP model of the AMM
 * (almost-Montgomery) product the kernels compute:
 *   U = (a*b + m*N) / R, m = -a*b/N mod R;  U -= N if U >= R
 * in all 8 lanes, with random inputs in [0, R), worst-case inputs (R-1),
 * aliasing (c == a, c == b), long dependent chains and powm vs mpz_powm.
 */

#include "avx_ecm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void *xmalloc_align(size_t len)
{
    void *p;
    if (posix_memalign(&p, 64, (len + 63) & ~(size_t)63) != 0)
    {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    return p;
}

typedef void (*mul_fn)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
typedef void (*sqr_fn)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_monty_t *);

int vecmulmod52_tmpl_get(int S, mul_fn *mul, sqr_fn *sqr);

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// ---------------------------------------------------------------------------
// IFMA calibration: 1 IFMA per cycle on p0 -> ns per cycle
// ---------------------------------------------------------------------------
static double calibrate_cycle_ns(long iters)
{
    long i;
    double t0, best = 1e9;
    int r;
    for (r = 0; r < 3; r++)
    {
        t0 = now_sec();
        for (i = 0; i < iters; i++)
        {
            __asm__ volatile(
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm10\n\t"
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm11\n\t"
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm12\n\t"
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm13\n\t"
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm14\n\t"
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm15\n\t"
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm16\n\t"
                "vpmadd52luq %%zmm1, %%zmm2, %%zmm17\n\t"
                ::: "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15",
                "xmm16", "xmm17");
        }
        t0 = (now_sec() - t0) / (iters * 8.0) * 1e9;
        if (t0 < best)
            best = t0;
    }
    return best;
}

// ---------------------------------------------------------------------------
// setup (mirrors avx_ecm_main.c / rns_bench.c)
// ---------------------------------------------------------------------------
typedef struct
{
    int words;
    vec_monty_t *md;
    vec_bignum_t *n, *s;
    mpz_t N, R, Rinv, nhat;
    mul_fn mul;
    sqr_fn sqr;
    const char *core;
} ae_t;

static void vec_insert(vec_bignum_t *v, mpz_t src, int lane, int words)
{
    mpz_t t;
    int i;
    mpz_init_set(t, src);
    for (i = 0; i < words; i++)
    {
        v->data[lane + i * VECLEN] = mpz_get_ui(t) & VEC_MAXDIGIT;
        mpz_tdiv_q_2exp(t, t, DIGITBITS);
    }
    v->size = words;
    mpz_clear(t);
}

static void vec_extract(mpz_t dest, vec_bignum_t *v, int lane, int words)
{
    int j;
    mpz_set_ui(dest, 0);
    for (j = words - 1; j >= 0; j--)
    {
        mpz_mul_2exp(dest, dest, DIGITBITS);
        mpz_add_ui(dest, dest, v->data[lane + j * VECLEN]);
    }
}

// returns 0 if no "old" kernel exists for this size
static int pick_old(ae_t *ae)
{
    int maxbits = ae->words * DIGITBITS;
    if (maxbits == 416)
    {
        ae->mul = &vecmulmod52_fixed416_bfips;
        ae->sqr = &vecsqrmod52_fixed416_bfips;
        ae->core = "fixed416";
    }
    else if (maxbits == 624)
    {
        ae->mul = &vecmulmod52_fixed624_bfips;
        ae->sqr = &vecsqrmod52_fixed624_bfips;
        ae->core = "fixed624";
    }
    else if (maxbits == 832)
    {
        ae->mul = &vecmulmod52_fixed832_bfips;
        ae->sqr = &vecsqrmod52_fixed832_bfips;
        ae->core = "fixed832";
    }
    else if (maxbits == 1040)
    {
        ae->mul = &vecmulmod52_fixed1040_bfips;
        ae->sqr = &vecsqrmod52_fixed1040_bfips;
        ae->core = "fixed1040";
    }
    else if ((ae->words % 4) == 0)
    {
        ae->mul = &vecmulmod52;
        ae->sqr = &vecsqrmod52;
        ae->core = "bfips";
    }
    else
        return 0;
    return 1;
}

static int pick_new(ae_t *ae)
{
    ae->core = "tmpl";
    return vecmulmod52_tmpl_get(ae->words, &ae->mul, &ae->sqr);
}

// negative self-test: a kernel that corrupts one word every 50th call
static mul_fn neg_mul_real;
static sqr_fn neg_sqr_real;
static long neg_calls = 0;
static void neg_corrupt(vec_bignum_t *c, vec_monty_t *md)
{
    if ((++neg_calls % 50) == 0)
        c->data[5 + (md->NWORDS / 2) * VECLEN] ^= 1ULL << (neg_calls % 52);
}
static void neg_mul(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *md)
{
    neg_mul_real(a, b, c, n, s, md);
    neg_corrupt(c, md);
}
static void neg_sqr(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n,
    vec_bignum_t *s, vec_monty_t *md)
{
    neg_sqr_real(a, c, n, s, md);
    neg_corrupt(c, md);
}
static int pick_neg(ae_t *ae)
{
    if (!pick_new(ae))
        return 0;
    neg_mul_real = ae->mul;
    neg_sqr_real = ae->sqr;
    ae->mul = neg_mul;
    ae->sqr = neg_sqr;
    ae->core = "corrupt";
    return 1;
}

static void ae_init(ae_t *ae, mpz_t N, int words)
{
    vec_monty_t *md;
    int j;

    ae->words = words;
    md = ae->md = vec_monty_alloc(words);
    md->NBLOCKS = words / BLOCKWORDS;
    md->MAXBITS = words * DIGITBITS;
    md->NWORDS = words;
    md->isMersenne = 0;
    md->nbits = mpz_sizeinbase(N, 2);
    md->use_vnbits = 0;

    mpz_init_set(ae->N, N);
    mpz_init(ae->R);
    mpz_init(ae->Rinv);
    mpz_init(ae->nhat);
    mpz_set_ui(ae->R, 1);
    mpz_mul_2exp(ae->R, ae->R, words * DIGITBITS);
    mpz_invert(md->nhat, N, ae->R);
    mpz_sub(md->nhat, ae->R, md->nhat);
    mpz_set(ae->nhat, md->nhat);
    mpz_mul(md->rhat, ae->R, ae->R);
    mpz_tdiv_r(md->rhat, md->rhat, N);
    ae->n = vecInit(words);
    ae->s = vecInit(words);
    for (j = 0; j < VECLEN; j++)
    {
        vec_insert(ae->n, N, j, words);
        vec_insert(md->n, N, j, words);
        vec_insert(md->vrhat, md->rhat, j, words);
        vec_insert(md->vnhat, md->nhat, j, words);
        md->vrho[j] = mpz_get_ui(md->nhat) & VEC_MAXDIGIT;
    }
    mpz_tdiv_r(md->rhat, ae->R, N);         // R mod N = Montgomery 1
    for (j = 0; j < VECLEN; j++)
        vec_insert(md->one, md->rhat, j, words);
    mpz_invert(ae->Rinv, ae->R, N);
}

static void ae_free(ae_t *ae)
{
    vecFree(ae->n);
    vecFree(ae->s);
    vec_monty_free(ae->md);
    free(ae->md);
    mpz_clear(ae->N);
    mpz_clear(ae->R);
    mpz_clear(ae->Rinv);
    mpz_clear(ae->nhat);
}

// GMP model of the AMM product
static void ref_amm(ae_t *ae, mpz_t u, mpz_t a, mpz_t b)
{
    mpz_t t, m;
    int bits = ae->words * DIGITBITS;
    mpz_init(t);
    mpz_init(m);
    mpz_mul(t, a, b);
    mpz_tdiv_r_2exp(m, t, bits);
    mpz_mul(m, m, ae->nhat);
    mpz_tdiv_r_2exp(m, m, bits);
    mpz_addmul(t, m, ae->N);
    mpz_tdiv_q_2exp(u, t, bits);
    if (mpz_cmp(u, ae->R) >= 0)
        mpz_sub(u, u, ae->N);
    mpz_clear(t);
    mpz_clear(m);
}

// ---------------------------------------------------------------------------
// correctness
// ---------------------------------------------------------------------------
static int check_lanes(ae_t *ae, vec_bignum_t *v, mpz_t *ref, const char *what)
{
    mpz_t t;
    int l, j, fails = 0;
    mpz_init(t);
    for (l = 0; l < VECLEN; l++)
    {
        vec_extract(t, v, l, ae->words);
        for (j = 0; j < ae->words; j++)
            if (v->data[l + j * VECLEN] > VEC_MAXDIGIT)
                break;
        if (mpz_cmp(t, ref[l]) != 0 || j < ae->words)
        {
            if (fails < 3)
                printf("  FAIL %s S=%d %s lane %d%s\n", ae->core, ae->words,
                    what, l, j < ae->words ? " (unnormalized word)" : "");
            fails++;
        }
    }
    mpz_clear(t);
    return fails;
}

static void set_inputs(ae_t *ae, vec_bignum_t *va, mpz_t *a, gmp_randstate_t rs,
    int special)
{
    int l;
    for (l = 0; l < VECLEN; l++)
    {
        mpz_urandomb(a[l], rs, ae->words * DIGITBITS);   // [0, R)
        if (special)
        {
            // worst cases: R-1, N-1, 0, 1, random high-weight
            switch (l)
            {
            case 0: mpz_sub_ui(a[l], ae->R, 1); break;
            case 1: mpz_sub_ui(a[l], ae->R, 1); break;
            case 2: mpz_sub_ui(a[l], ae->N, 1); break;
            case 3: mpz_set_ui(a[l], 0); break;
            case 4: mpz_set_ui(a[l], 1); break;
            case 5: mpz_rrandomb(a[l], rs, ae->words * DIGITBITS); break;
            default: break;
            }
        }
        vec_insert(va, a[l], l, ae->words);
    }
}

static int test_kernel(ae_t *ae, gmp_randstate_t rs, int iters, long *nslow)
{
    int S = ae->words, l, it, fails = 0;
    vec_bignum_t *va = vecInit(S), *vb = vecInit(S), *vc = vecInit(S);
    mpz_t a[VECLEN], b[VECLEN], r[VECLEN], x[VECLEN];

    for (l = 0; l < VECLEN; l++)
    {
        mpz_init(a[l]); mpz_init(b[l]); mpz_init(r[l]); mpz_init(x[l]);
    }

    for (it = 0; it < iters && fails == 0; it++)
    {
        int special = (it % 4) == 0;
        set_inputs(ae, va, a, rs, special);
        set_inputs(ae, vb, b, rs, special && (it % 8) == 0);

        // mul, c distinct
        for (l = 0; l < VECLEN; l++)
        {
            // the model itself: U*R == a*b (mod N) and U < R
            mpz_t t;
            ref_amm(ae, r[l], a[l], b[l]);
            mpz_init(t);
            mpz_mul(x[l], a[l], b[l]);
            mpz_mod(x[l], x[l], ae->N);
            mpz_mul(t, r[l], ae->R);
            mpz_mod(t, t, ae->N);
            if (mpz_cmp(t, x[l]) != 0 || mpz_cmp(r[l], ae->R) >= 0)
            {
                printf("  ref model inconsistent\n");
                fails++;
            }
            mpz_clear(t);
        }
        ae->mul(va, vb, vc, ae->n, ae->s, ae->md);
        fails += check_lanes(ae, vc, r, "mul");

        // sqr
        for (l = 0; l < VECLEN; l++)
            ref_amm(ae, r[l], a[l], a[l]);
        ae->sqr(va, vc, ae->n, ae->s, ae->md);
        fails += check_lanes(ae, vc, r, "sqr");

        // aliasing: c == a (mul), c == b (mul), c == a (sqr)
        for (l = 0; l < VECLEN; l++)
            ref_amm(ae, r[l], a[l], b[l]);
        vecCopyn(va, vc, S);
        ae->mul(vc, vb, vc, ae->n, ae->s, ae->md);
        fails += check_lanes(ae, vc, r, "mul c==a");
        vecCopyn(vb, vc, S);
        ae->mul(va, vc, vc, ae->n, ae->s, ae->md);
        fails += check_lanes(ae, vc, r, "mul c==b");
        for (l = 0; l < VECLEN; l++)
            ref_amm(ae, r[l], a[l], a[l]);
        vecCopyn(va, vc, S);
        ae->sqr(vc, vc, ae->n, ae->s, ae->md);
        fails += check_lanes(ae, vc, r, "sqr c==a");
    }

    // long dependent chains: x = x^2, then x = x*b, checked every step
    set_inputs(ae, va, a, rs, 0);
    set_inputs(ae, vb, b, rs, 0);
    for (l = 0; l < VECLEN; l++)
        mpz_set(x[l], a[l]);
    for (it = 0; it < 4 * iters && fails == 0; it++)
    {
        int domul = it & 1;
        for (l = 0; l < VECLEN; l++)
        {
            mpz_t u;
            mpz_init(u);
            ref_amm(ae, u, x[l], domul ? b[l] : x[l]);
            if (mpz_cmp(u, ae->R) >= 0)
                printf("  ref output >= R?\n");
            mpz_set(x[l], u);
            mpz_clear(u);
        }
        if (domul)
            ae->mul(va, vb, va, ae->n, ae->s, ae->md);
        else
            ae->sqr(va, va, ae->n, ae->s, ae->md);
        fails += check_lanes(ae, va, x, domul ? "chain mul" : "chain sqr");
    }

    // how often the final subtraction fires on random [0, R) inputs
    if (nslow)
    {
        // (a*b + m*N)/R >= R before the subtraction
        long cnt = 0;
        mpz_t a1, b1, u, m;
        mpz_init(a1); mpz_init(b1); mpz_init(u); mpz_init(m);
        for (it = 0; it < 2000; it++)
        {
            mpz_urandomb(a1, rs, S * DIGITBITS);
            mpz_urandomb(b1, rs, S * DIGITBITS);
            mpz_mul(u, a1, b1);
            mpz_tdiv_r_2exp(m, u, S * DIGITBITS);
            mpz_mul(m, m, ae->nhat);
            mpz_tdiv_r_2exp(m, m, S * DIGITBITS);
            mpz_addmul(u, m, ae->N);
            mpz_tdiv_q_2exp(u, u, S * DIGITBITS);
            if (mpz_cmp(u, ae->R) >= 0)
                cnt++;
        }
        *nslow = cnt;
        mpz_clear(a1); mpz_clear(b1); mpz_clear(u); mpz_clear(m);
    }

    for (l = 0; l < VECLEN; l++)
    {
        mpz_clear(a[l]); mpz_clear(b[l]); mpz_clear(r[l]); mpz_clear(x[l]);
    }
    vecFree(va); vecFree(vb); vecFree(vc);
    return fails;
}

// sliding-window powm (vecslidingmodexp), checked against mpz_powm
static void ae_powm(ae_t *ae, vec_bignum_t *d, vec_bignum_t *b, mpz_t e)
{
    vec_monty_t *md = ae->md;
    vec_bignum_t **g = md->g;
    int bit = (int)mpz_sizeinbase(e, 2) - 1;
    int k = get_winsize(bit);
    int i, j;

    vecCopyn(b, g[1], md->NWORDS);
    ae->sqr(b, g[2], ae->n, ae->s, md);
    for (i = 1; i < (1 << (k - 1)); i++)
        ae->mul(g[2 * i - 1], g[2], g[2 * i + 1], ae->n, ae->s, md);

    vecCopyn(md->one, d, md->NWORDS);
    while (bit >= 0)
    {
        if (mpz_tstbit(e, bit) == 0)
        {
            ae->sqr(d, d, ae->n, ae->s, md);
            bit--;
        }
        else
        {
            int a = bit - k + 1, thisk, bstr = 1;
            if (a < 0)
                a = 0;
            a = (int)mpz_scan1(e, a);
            thisk = bit - a + 1;
            for (j = 1; j < thisk; j++)
                bstr = (bstr << 1) | mpz_tstbit(e, bit - j);
            for (j = 0; j < thisk; j++)
                ae->sqr(d, d, ae->n, ae->s, md);
            ae->mul(d, g[bstr], d, ae->n, ae->s, md);
            bit -= thisk;
        }
    }
}

static int test_powm(ae_t *ae, gmp_randstate_t rs, double *us)
{
    int S = ae->words, l, fails = 0, reps;
    vec_bignum_t *vx = vecInit(S), *vd = vecInit(S);
    mpz_t a[VECLEN], r, t, e;
    double t0, el;
    int bits = (int)mpz_sizeinbase(ae->N, 2);

    mpz_init(r); mpz_init(t); mpz_init(e);
    for (l = 0; l < VECLEN; l++)
    {
        mpz_init(a[l]);
        mpz_urandomm(a[l], rs, ae->N);
        mpz_mul_2exp(t, a[l], S * DIGITBITS);
        mpz_mod(t, t, ae->N);
        vec_insert(vx, t, l, S);
    }
    mpz_urandomb(e, rs, bits);
    mpz_setbit(e, bits - 1);
    ae_powm(ae, vd, vx, e);
    for (l = 0; l < VECLEN; l++)
    {
        mpz_powm(r, a[l], e, ae->N);
        vec_extract(t, vd, l, S);
        mpz_mul(t, t, ae->Rinv);
        mpz_mod(t, t, ae->N);
        if (mpz_cmp(t, r) != 0)
        {
            printf("  FAIL %s S=%d powm lane %d\n", ae->core, S, l);
            fails++;
        }
    }
    if (us)
    {
        t0 = now_sec();
        for (reps = 0; (el = now_sec() - t0) < 0.4 || reps < 2; reps++)
            ae_powm(ae, vd, vx, e);
        *us = el / reps / VECLEN * 1e6;
    }
    for (l = 0; l < VECLEN; l++)
        mpz_clear(a[l]);
    mpz_clear(r); mpz_clear(t); mpz_clear(e);
    vecFree(vx); vecFree(vd);
    return fails;
}

// ---------------------------------------------------------------------------
// timing
// ---------------------------------------------------------------------------
typedef struct { double sqr, mul, sqr_ind, mul_ind; } tm_t;

static void time_kernel(ae_t *ae, gmp_randstate_t rs, tm_t *tm)
{
    int S = ae->words, l, i, tr, reps = 20000;
    vec_bignum_t *x = vecInit(S), *y = vecInit(S), *d = vecInit(S);
    vec_bignum_t *d2 = vecInit(S), *d3 = vecInit(S);
    mpz_t t;
    double t0, el;

    mpz_init(t);
    for (l = 0; l < VECLEN; l++)
    {
        mpz_urandomm(t, rs, ae->N);
        vec_insert(x, t, l, S);
        mpz_urandomm(t, rs, ae->N);
        vec_insert(y, t, l, S);
    }
    if (S > 40)
        reps = 5000;

    vecCopyn(x, d, S);
    t0 = now_sec();
    while (now_sec() - t0 < 0.15)
    {
        ae->sqr(d, d, ae->n, ae->s, ae->md);
        ae->mul(d, x, d, ae->n, ae->s, ae->md);
    }
    // min over several trials (timing noise on this machine is ~10%)
    tm->sqr = tm->mul = tm->sqr_ind = tm->mul_ind = 1e9;
    reps /= 4;
    for (tr = 0; tr < 7; tr++)
    {
        // dependent chains
        t0 = now_sec();
        for (i = 0; i < reps; i++)
            ae->sqr(d, d, ae->n, ae->s, ae->md);
        el = (now_sec() - t0) / reps;
        if (el < tm->sqr) tm->sqr = el;
        t0 = now_sec();
        for (i = 0; i < reps; i++)
            ae->mul(d, x, d, ae->n, ae->s, ae->md);
        el = (now_sec() - t0) / reps;
        if (el < tm->mul) tm->mul = el;
        // independent calls (no data dependency between consecutive calls)
        t0 = now_sec();
        for (i = 0; i < reps; i++)
        {
            ae->sqr(x, d2, ae->n, ae->s, ae->md);
            ae->sqr(y, d3, ae->n, ae->s, ae->md);
        }
        el = (now_sec() - t0) / reps / 2;
        if (el < tm->sqr_ind) tm->sqr_ind = el;
        t0 = now_sec();
        for (i = 0; i < reps; i++)
        {
            ae->mul(x, y, d2, ae->n, ae->s, ae->md);
            ae->mul(y, x, d3, ae->n, ae->s, ae->md);
        }
        el = (now_sec() - t0) / reps / 2;
        if (el < tm->mul_ind) tm->mul_ind = el;
    }

    mpz_clear(t);
    vecFree(x); vecFree(y); vecFree(d); vecFree(d2); vecFree(d3);
}

static void make_N(mpz_t N, int bits, gmp_randstate_t rs);

// vecaddmod52 / vecsubmod52 on AMM-product outputs (as ECM feeds them):
// count lanes where the result is not == a+-b mod N, or not < R
static void test_addsub(int S, gmp_randstate_t rs)
{
    int spares[] = { 1, 2, 3, 4, 6, 8, 16, 26 }, si, it, l;
    for (si = 0; si < (int)(sizeof(spares) / sizeof(spares[0])); si++)
    {
        ae_t ae;
        mpz_t N, t, u;
        long bad_add = 0, bad_sub = 0, tot = 0;
        vec_bignum_t *va = vecInit(S), *vb = vecInit(S), *vc = vecInit(S);
        mpz_t a[VECLEN], b[VECLEN];
        mpz_init(N); mpz_init(t); mpz_init(u);
        for (l = 0; l < VECLEN; l++) { mpz_init(a[l]); mpz_init(b[l]); }
        make_N(N, S * DIGITBITS - spares[si], rs);
        ae_init(&ae, N, S);
        pick_new(&ae);
        // seed values: AMM products of random residues, then keep
        // mixing (x = x*y, y = x+y, ...) like a curve-arithmetic chain
        for (l = 0; l < VECLEN; l++)
        {
            mpz_urandomm(a[l], rs, N);      // residues < N, as ECM starts
            mpz_urandomm(b[l], rs, N);
            vec_insert(va, a[l], l, S);
            vec_insert(vb, b[l], l, S);
        }
        for (it = 0; it < 4000; it++)
        {
            ae.mul(va, vb, va, ae.n, ae.s, ae.md);      // a = a*b (AMM)
            ae.sqr(vb, vb, ae.n, ae.s, ae.md);          // b = b^2
            for (l = 0; l < VECLEN; l++)
            {
                vec_extract(a[l], va, l, S);
                vec_extract(b[l], vb, l, S);
            }
            vecaddmod52(va, vb, vc, ae.md);
            for (l = 0; l < VECLEN; l++)
            {
                vec_extract(t, vc, l, S);
                mpz_add(u, a[l], b[l]);
                mpz_sub(u, u, t);
                if (!mpz_divisible_p(u, N) || mpz_cmp(t, ae.R) >= 0)
                    bad_add++;
            }
            vecsubmod52(va, vb, vc, ae.md);
            for (l = 0; l < VECLEN; l++)
            {
                vec_extract(t, vc, l, S);
                mpz_sub(u, a[l], b[l]);
                mpz_sub(u, u, t);
                if (!mpz_divisible_p(u, N) || mpz_cmp(t, ae.R) >= 0)
                    bad_sub++;
            }
            tot += VECLEN;
            // feed sums back in, as curve arithmetic does
            if (it & 1)
                vecaddmod52(va, vb, vb, ae.md);
        }
        printf("  S=%2d nbits %4d (spare %2d): addmod wrong %ld/%ld, submod wrong %ld/%ld\n",
            S, S * DIGITBITS - spares[si], spares[si], bad_add, tot, bad_sub, tot);
        for (l = 0; l < VECLEN; l++) { mpz_clear(a[l]); mpz_clear(b[l]); }
        vecFree(va); vecFree(vb); vecFree(vc);
        ae_free(&ae);
        mpz_clear(N); mpz_clear(t); mpz_clear(u);
    }
}

static void make_N(mpz_t N, int bits, gmp_randstate_t rs)
{
    mpz_urandomb(N, rs, bits);
    mpz_setbit(N, bits - 1);
    mpz_setbit(N, 0);
}

int main(int argc, char **argv)
{
    int def_words[] = { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20,
        22, 24, 28, 30, 32, 36, 40, 48, 60, 80 };
    int words[64], nw = 0, w, i, do_test = 1, do_bench = 1, do_powm = 0, do_addsub = 0, do_neg = 0;
    int spare = -1;
    unsigned long seed = (unsigned long)time(NULL);
    gmp_randstate_t rs;
    int fails = 0;
    double cyc;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-s") && i + 1 < argc)
            seed = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-t"))
            do_bench = 0;
        else if (!strcmp(argv[i], "-b"))
            do_test = 0;
        else if (!strcmp(argv[i], "-p"))
            do_powm = 1;
        else if (!strcmp(argv[i], "-neg"))
            do_neg = 1, do_bench = 0;
        else if (!strcmp(argv[i], "-a"))
            do_addsub = 1, do_test = do_bench = 0;
        else if (!strcmp(argv[i], "-spare") && i + 1 < argc)
            spare = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-w") && i + 1 < argc)
        {
            char *p = argv[++i];
            while (*p && nw < 64)
            {
                words[nw++] = (int)strtol(p, &p, 10);
                if (*p == ',')
                    p++;
            }
        }
    }
    if (nw == 0)
    {
        nw = (int)(sizeof(def_words) / sizeof(def_words[0]));
        memcpy(words, def_words, sizeof(def_words));
    }

    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, seed);
    cyc = calibrate_cycle_ns(20000000);
    printf("seed %lu; IFMA-calibrated cycle %.4f ns (%.2f GHz)\n", seed, cyc, 1.0 / cyc);

    if (do_addsub)
    {
        printf("addmod/submod on AMM outputs (result == a+-b mod N and < R)\n");
        for (w = 0; w < nw; w++)
            test_addsub(words[w], rs);
    }

    if (do_test)
    {
        printf("correctness (bit-exact AMM vs GMP, all lanes, worst cases, aliasing, chains, powm)\n");
        for (w = 0; w < nw; w++)
        {
            int S = words[w], v, k;
            for (v = 0; v < 2; v++)
            {
                for (k = 0; k < 3; k++)
                {
                    // N sizes: 1 spare bit, 2 spare bits, mid-size
                    int sp = (spare >= 0) ? spare : (k == 0 ? 1 : k == 1 ? 2 : 26);
                    ae_t ae;
                    mpz_t N;
                    long nslow = 0;
                    int f;
                    mpz_init(N);
                    make_N(N, S * DIGITBITS - sp, rs);
                    ae_init(&ae, N, S);
                    if ((v == 0 ? (do_neg ? pick_neg(&ae) : pick_old(&ae)) : pick_new(&ae)) == 0)
                    {
                        ae_free(&ae);
                        mpz_clear(N);
                        break;
                    }
                    f = test_kernel(&ae, rs, 200, &nslow);
                    f += test_powm(&ae, rs, NULL);
                    printf("  S=%2d %-9s nbits %4d: %s (final-sub on random [0,R) inputs: %ld/2000)\n",
                        S, ae.core, S * DIGITBITS - sp, f ? "FAIL" : "ok", nslow);
                    fails += f;
                    ae_free(&ae);
                    mpz_clear(N);
                    if (spare >= 0)
                        break;
                }
            }
        }
        printf("%s: %d failures\n", fails ? "TESTS FAILED" : "TESTS PASSED", fails);
    }

    if (do_bench)
    {
        printf("\ntimes per 8-lane call: ns (cycles); dep = dependent chain, ind = independent calls\n");
        printf("floor = IFMA count (mul 4S^2+S, sqr 3S^2+2S) in cycles; eff = floor / dep cycles\n");
        printf("  S  bits | core      |  sqr dep        ind  floor  eff |  mul dep        ind  floor  eff | per-number sqr/mul ns%s\n",
            do_powm ? " | powm us" : "");
        for (w = 0; w < nw; w++)
        {
            int S = words[w], v;
            for (v = 0; v < 2; v++)
            {
                ae_t ae;
                mpz_t N;
                tm_t tm;
                double pw = 0;
                double fs = 3.0 * S * S + 2 * S, fm = 4.0 * S * S + S;
                mpz_init(N);
                make_N(N, S * DIGITBITS - 1, rs);
                ae_init(&ae, N, S);
                if ((v == 0 ? pick_old(&ae) : pick_new(&ae)) == 0)
                {
                    ae_free(&ae);
                    mpz_clear(N);
                    continue;
                }
                time_kernel(&ae, rs, &tm);
                cyc = calibrate_cycle_ns(2000000);   // clock drifts; re-measure per row
                if (do_powm)
                    fails += test_powm(&ae, rs, &pw);
                printf(" %2d %5d | %-9s | %6.0f (%5.0f) %5.0f %5.0f %3.0f%% | %6.0f (%5.0f) %5.0f %5.0f %3.0f%% | %6.1f %6.1f",
                    S, S * DIGITBITS, ae.core,
                    tm.sqr * 1e9, tm.sqr * 1e9 / cyc, tm.sqr_ind * 1e9 / cyc, fs, 100 * fs / (tm.sqr * 1e9 / cyc),
                    tm.mul * 1e9, tm.mul * 1e9 / cyc, tm.mul_ind * 1e9 / cyc, fm, 100 * fm / (tm.mul * 1e9 / cyc),
                    tm.sqr * 1e9 / VECLEN, tm.mul * 1e9 / VECLEN);
                if (do_powm)
                    printf(" | %8.2f", pw);
                printf("\n");
                fflush(stdout);
                ae_free(&ae);
                mpz_clear(N);
            }
        }
    }
    gmp_randclear(rs);
    return fails != 0;
}
