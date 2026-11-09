/*
 * Benchmark: AVX-ECM vecarith52 (old fixed/BFIPS cores and the constant-S
 * template cores), the RNS kernels (small regime, Yao one-number, Yao
 * batch) and GMP.
 *   rns_bench [seed]
 * Every vector engine processes 8 numbers modulo the same N; times are
 * reported per number.  Modexp uses one sliding-window rule for all vector
 * engines (AVX-ECM's vecslidingmodexp), with one full-size exponent shared
 * by all lanes, and every engine's result is checked against mpz_powm in
 * every lane.
 */

#include "rns.h"
#include "rns_yao.h"
#include "avx_ecm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// AVX-ECM's allocator lives in ytools; a local stand-in avoids linking it
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

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// ---------------------------------------------------------------------------
// AVX-ECM glue (mirrors test_vecarith / assign_mul_ptrs in avx_ecm_main.c)
// ---------------------------------------------------------------------------

typedef void (*ae_mul_t)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
typedef void (*ae_sqr_t)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_monty_t *);

typedef struct
{
    int words;
    vec_monty_t *mdata;
    vec_bignum_t *n, *s;
    mpz_t R, Rinv;
    ae_mul_t mul;
    ae_sqr_t sqr;
    char core[16];
} avxecm_t;

static void vec_insert(vec_bignum_t *v, mpz_t src, int lane)
{
    mpz_t t;
    int i;
    mpz_init_set(t, src);
    for (i = 0; i < (int)v->WORDS_ALLOC; i++)
    {
        v->data[lane + i * VECLEN] = mpz_get_ui(t) & VEC_MAXDIGIT;
        mpz_tdiv_q_2exp(t, t, DIGITBITS);
    }
    v->size = v->WORDS_ALLOC;
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

// tmpl = 0: the old cores (fixed416/624/832/1040 BFIPS, else BFIPS), with
// NWORDS a multiple of 4.  tmpl = 1: the constant-S template cores with
// the exact word count.  Both keep >= AVXECM_SPARE_BITS spare bits.
static void avxecm_init(avxecm_t *ae, mpz_t N, int tmpl)
{
    vec_monty_t *md;
    int nbits = (int)mpz_sizeinbase(N, 2), j;
    int words = (nbits + AVXECM_SPARE_BITS + DIGITBITS - 1) / DIGITBITS;

    if (!tmpl)
        words = (words + BLOCKWORDS - 1) / BLOCKWORDS * BLOCKWORDS;
    else if (words < 2)
        words = 2;
    ae->words = words;

    if (tmpl)
    {
        if (!vecmulmod52_tmpl_get(words, &ae->mul, &ae->sqr))
        {
            fprintf(stderr, "no template core for %d words\n", words);
            exit(1);
        }
        snprintf(ae->core, sizeof(ae->core), "tmpl%s", words > 32 ? "-sN" : "");
    }
    else
    {
        switch (words * DIGITBITS)
        {
        case 416: ae->mul = vecmulmod52_fixed416_bfips; ae->sqr = vecsqrmod52_fixed416_bfips;
            strcpy(ae->core, "fixed416"); break;
        case 624: ae->mul = vecmulmod52_fixed624_bfips; ae->sqr = vecsqrmod52_fixed624_bfips;
            strcpy(ae->core, "fixed624"); break;
        case 832: ae->mul = vecmulmod52_fixed832_bfips; ae->sqr = vecsqrmod52_fixed832_bfips;
            strcpy(ae->core, "fixed832"); break;
        case 1040: ae->mul = vecmulmod52_fixed1040_bfips; ae->sqr = vecsqrmod52_fixed1040_bfips;
            strcpy(ae->core, "fixed1040"); break;
        default: ae->mul = vecmulmod52; ae->sqr = vecsqrmod52;
            strcpy(ae->core, "bfips"); break;
        }
    }

    md = ae->mdata = vec_monty_alloc(words);
    md->NBLOCKS = (words + BLOCKWORDS - 1) / BLOCKWORDS;
    md->MAXBITS = words * DIGITBITS;
    md->NWORDS = words;
    md->isMersenne = 0;
    md->nbits = nbits;
    md->use_vnbits = 0;

    mpz_init(ae->R);
    mpz_init(ae->Rinv);
    mpz_set_ui(ae->R, 1);
    mpz_mul_2exp(ae->R, ae->R, words * DIGITBITS);
    mpz_invert(md->nhat, N, ae->R);
    mpz_sub(md->nhat, ae->R, md->nhat);
    mpz_mul(md->rhat, ae->R, ae->R);
    mpz_tdiv_r(md->rhat, md->rhat, N);
    ae->n = vecInit(words);
    ae->s = vecInit(words);
    for (j = 0; j < VECLEN; j++)
    {
        vec_insert(ae->n, N, j);
        vec_insert(md->n, N, j);
        vec_insert(md->vrhat, md->rhat, j);
        vec_insert(md->vnhat, md->nhat, j);
        md->vrho[j] = mpz_get_ui(md->nhat) & VEC_MAXDIGIT;
    }
    mpz_tdiv_r(md->rhat, ae->R, N);         // R mod N = Montgomery 1
    for (j = 0; j < VECLEN; j++)
        vec_insert(md->one, md->rhat, j);
    mpz_invert(ae->Rinv, ae->R, N);
}

static void avxecm_free(avxecm_t *ae)
{
    vecFree(ae->n);
    vecFree(ae->s);
    vec_monty_free(ae->mdata);
    free(ae->mdata);
    mpz_clear(ae->R);
    mpz_clear(ae->Rinv);
}

// vecslidingmodexp (avxppm1.c) restated so the bench doesn't link P-1
static void avxecm_powm(avxecm_t *ae, vec_bignum_t *d, vec_bignum_t *b, mpz_t e)
{
    vec_monty_t *md = ae->mdata;
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

// ---------------------------------------------------------------------------
// engines and timing
// ---------------------------------------------------------------------------

enum { E_AEOLD, E_AETMPL, E_SMALL, E_YAO1, E_YAO8, E_GMP, NENG };
// AVX2 builds: the old AVX-ECM cores, the small-regime RNS kernel and the
// one-number Yao kernel are AVX-512 only (shown as n/a)
#ifdef RNS_HAVE_AVX512
static const int have_eng[NENG] = { 1, 1, 1, 1, 1, 1 };
#else
static const int have_eng[NENG] = { 0, 1, 0, 0, 1, 1 };
#endif
static const char *ename[NENG] = { "AE-old", "AE-tmpl", "RNS-sm", "Yao1", "Yao8", "GMP" };

typedef struct
{
    int bits;
    double sqr[NENG], mul[NENG], powm[NENG];   // ns, ns, us per number
    char aeold[16], aetmpl[16];
    int wold, wtmpl, smallnB, smallnC, yao1s, yao8n, yao8s;
} result_t;

// ns per call of op over reps iterations, after a 0.2 s warm-up
#define TIME_NS(out, reps, op) do {                                     \
        double _t0 = now_sec();                                         \
        int _i;                                                         \
        while (now_sec() - _t0 < 0.2) { op; }                           \
        _t0 = now_sec();                                                \
        for (_i = 0; _i < (reps); _i++) { op; }                         \
        (out) = (now_sec() - _t0) / (reps) * 1e9;                       \
    } while (0)

// us per call, repeated for at least 0.5 s
#define TIME_US(out, op) do {                                           \
        double _t0 = now_sec(), _el;                                    \
        int _r;                                                         \
        for (_r = 0; (_el = now_sec() - _t0) < 0.5 || _r < 2; _r++) { op; } \
        (out) = _el / _r * 1e6;                                         \
    } while (0)

static int check_lane(mpz_t got, mpz_t expect, const char *eng, int bits, int lane)
{
    if (mpz_cmp(got, expect) == 0)
        return 0;
    printf("FAIL %s powm %d bits lane %d\n", eng, bits, lane);
    return 1;
}

static int run_size(int bits, gmp_randstate_t rs, result_t *res)
{
    rns_fixed_t *fx = rns_fixed_init(bits);
    rnsy_fixed_t *yfx = rnsy_fixed_init(bits);
    rnsy_fixed_t *bfx = rnsy_fixed_init_n(bits, rnsyb_choose_n(bits));
    rns_ctx_t *ctx;
    rnsy_ctx_t *yctx, *bctx;
    avxecm_t ae[2];
    mpz_t N, e, t, a[VECLEN], r[VECLEN];
    uint64_t *vx, *vd, *tmp, *yx, *yd, *bx, *bd, *bt;
    vec_bignum_t *ax[2], *ad[2];
    const int reps = 20000;
    int l, k, fails = 0;

    mpz_init(N); mpz_init(e); mpz_init(t);
    for (l = 0; l < VECLEN; l++)
    {
        mpz_init(a[l]);
        mpz_init(r[l]);
    }
    // N must be coprime to every RNS setup's moduli
    for (;;)
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        if ((ctx = rns_ctx_init(fx, N)) == NULL)
            continue;
        if ((yctx = rnsy_ctx_init(yfx, N)) != NULL)
        {
            if ((bctx = rnsy_ctx_init(bfx, N)) != NULL)
                break;
            rnsy_ctx_free(yctx);
        }
        rns_ctx_free(ctx);
    }
    for (k = 0; k < 2; k++)
    {
        avxecm_init(&ae[k], N, k);
        ax[k] = vecInit(ae[k].words);
        ad[k] = vecInit(ae[k].words);
    }

    vx = rns_valloc(ctx); vd = rns_valloc(ctx); tmp = rns_alloc(ctx);
    yx = rnsy_alloc(yctx); yd = rnsy_alloc(yctx);
    bx = rnsyb_alloc(bctx); bd = rnsyb_alloc(bctx); bt = rnsy_alloc(bctx);

    for (l = 0; l < VECLEN; l++)
    {
        mpz_urandomm(a[l], rs, N);
        rns_to_monty(ctx, tmp, a[l]);
        rns_vset_lane(ctx, vx, l, tmp);
        rnsy_to_monty(bctx, bt, a[l]);
        rnsyb_vset_lane(bctx, bx, l, bt);
        for (k = 0; k < 2; k++)
        {
            mpz_mul_2exp(t, a[l], ae[k].words * DIGITBITS);
            mpz_mod(t, t, N);
            vec_insert(ax[k], t, l);
        }
    }
    rnsy_to_monty(yctx, yx, a[0]);

    res->bits = bits;
    strcpy(res->aeold, ae[0].core);
    strcpy(res->aetmpl, ae[1].core);
    res->wold = ae[0].words;
    res->wtmpl = ae[1].words;
    res->smallnB = fx->nB;
    res->smallnC = fx->nC;
    res->yao1s = yfx->s;
    res->yao8n = bfx->n;
    res->yao8s = bfx->s;

    // sqrmod / mulmod, ns per number
    for (k = 0; k < 2; k++)
    {
        vec_monty_t *md = ae[k].mdata;
        if (!have_eng[k])
            continue;
        vecCopyn(ax[k], ad[k], ae[k].words);
        TIME_NS(res->sqr[k], reps, ae[k].sqr(ad[k], ad[k], ae[k].n, ae[k].s, md));
        TIME_NS(res->mul[k], reps, ae[k].mul(ad[k], ax[k], ad[k], ae[k].n, ae[k].s, md));
        res->sqr[k] /= VECLEN;
        res->mul[k] /= VECLEN;
    }
    memcpy(vd, vx, (size_t)fx->nch * VECLEN * sizeof(uint64_t));
#ifdef RNS_HAVE_AVX512
    TIME_NS(res->sqr[E_SMALL], reps, rns_vsqrmod_avx512(ctx, vd, vd));
    TIME_NS(res->mul[E_SMALL], reps, rns_vmulmod_avx512(ctx, vd, vd, vx));
    res->sqr[E_SMALL] /= VECLEN;
    res->mul[E_SMALL] /= VECLEN;
    memcpy(yd, yx, (size_t)yfx->nwords * sizeof(uint64_t));
    TIME_NS(res->sqr[E_YAO1], reps, rnsy_vsqrmod_avx512(yctx, yd, yd));
    TIME_NS(res->mul[E_YAO1], reps, rnsy_vmulmod_avx512(yctx, yd, yd, yx));
#endif
    memcpy(bd, bx, (size_t)bfx->nwords * VECLEN * sizeof(uint64_t));
    TIME_NS(res->sqr[E_YAO8], reps, rnsyb_vsqrmod(bctx, bd, bd));
    TIME_NS(res->mul[E_YAO8], reps, rnsyb_vmulmod(bctx, bd, bd, bx));
    res->sqr[E_YAO8] /= VECLEN;
    res->mul[E_YAO8] /= VECLEN;
    {
        // GMP: mul + tdiv_r and sqr + tdiv_r, for reference
        mpz_t q, x;
        mpz_init(q);
        mpz_init_set(x, a[0]);
        TIME_NS(res->sqr[E_GMP], reps, (mpz_mul(q, x, x), mpz_tdiv_r(x, q, N)));
        TIME_NS(res->mul[E_GMP], reps, (mpz_mul(q, x, a[1]), mpz_tdiv_r(x, q, N)));
        mpz_clear(q);
        mpz_clear(x);
    }

    // modexp with a full-size exponent shared by all lanes, us per number
    mpz_urandomb(e, rs, bits);
    mpz_setbit(e, bits - 1);
    for (k = 0; k < 2; k++)
    {
        if (!have_eng[k])
            continue;
        TIME_US(res->powm[k], avxecm_powm(&ae[k], ad[k], ax[k], e));
        res->powm[k] /= VECLEN;
    }
#ifdef RNS_HAVE_AVX512
    TIME_US(res->powm[E_SMALL], rns_vpowm_avx512(ctx, vd, vx, e));
    res->powm[E_SMALL] /= VECLEN;
    TIME_US(res->powm[E_YAO1], rnsy_vpowm_avx512(yctx, yd, yx, e));
#endif
    TIME_US(res->powm[E_YAO8], rnsyb_vpowm(bctx, bd, bx, e));
    res->powm[E_YAO8] /= VECLEN;
    TIME_US(res->powm[E_GMP], for (l = 0; l < VECLEN; l++) mpz_powm(r[l], a[l], e, N));
    res->powm[E_GMP] /= VECLEN;

    // every engine, every lane, against mpz_powm
    for (l = 0; l < VECLEN; l++)
    {
        for (k = 0; k < 2; k++)
        {
            if (!have_eng[k])
                continue;
            vec_extract(t, ad[k], l, ae[k].words);
            mpz_mul(t, t, ae[k].Rinv);
            mpz_mod(t, t, N);
            fails += check_lane(t, r[l], ename[k], bits, l);
        }
#ifdef RNS_HAVE_AVX512
        rns_vget_lane(ctx, tmp, vd, l);
        rns_from_monty(ctx, t, tmp);
        fails += check_lane(t, r[l], ename[E_SMALL], bits, l);
#endif
        rnsyb_vget_lane(bctx, bt, bd, l);
        rnsy_from_monty(bctx, t, bt);
        fails += check_lane(t, r[l], ename[E_YAO8], bits, l);
    }
#ifdef RNS_HAVE_AVX512
    rnsy_from_monty(yctx, t, yd);
    fails += check_lane(t, r[0], ename[E_YAO1], bits, 0);
#endif

    printf("  %5d bits done\n", bits);
    fflush(stdout);

    for (l = 0; l < VECLEN; l++)
    {
        mpz_clear(a[l]);
        mpz_clear(r[l]);
    }
    for (k = 0; k < 2; k++)
    {
        vecFree(ax[k]);
        vecFree(ad[k]);
        avxecm_free(&ae[k]);
    }
    rns_vfree(vx); rns_vfree(vd); free(tmp);
    rns_vfree(yx); rns_vfree(yd);
    rns_vfree(bx); rns_vfree(bd); rns_vfree(bt);
    rns_ctx_free(ctx);
    rnsy_ctx_free(yctx);
    rnsy_ctx_free(bctx);
    rns_fixed_free(fx);
    rnsy_fixed_free(yfx);
    rnsy_fixed_free(bfx);
    mpz_clear(N); mpz_clear(e); mpz_clear(t);
    return fails;
}

static int best_vec(const double *v)
{
    int k, best = 0;
    for (k = 0; k < E_GMP; k++)
        if (have_eng[k] && (!have_eng[best] || v[k] < v[best]))
            best = k;
    return best;
}

int main(int argc, char **argv)
{
    int sizes[] = { 256, 384, 512, 640, 768, 1024, 1280, 1536, 1792, 2048,
        2560, 3072, 3584, 4096 };
    const int ns = (int)(sizeof(sizes) / sizeof(sizes[0]));
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : (unsigned long)time(NULL);
    result_t res[sizeof(sizes) / sizeof(sizes[0])];
    gmp_randstate_t rs;
    int s, k, fails = 0;

    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, seed);
    printf("seed %lu\n", seed);
    for (s = 0; s < ns; s++)
        fails += run_size(sizes[s], rs, &res[s]);

    printf("\nConfigurations (words per number for AVX-ECM; RNS channel layout):\n");
    printf(" bits | AE-old        | AE-tmpl      | RNS-sm nB+nC | Yao1 s | Yao8 n/s\n");
    for (s = 0; s < ns; s++)
        printf("%5d | %-9s %3d | %-8s %3d | %5d+%-5d  | %4d   | %3d/%-2d\n", res[s].bits,
            res[s].aeold, res[s].wold, res[s].aetmpl, res[s].wtmpl,
            res[s].smallnB, res[s].smallnC, res[s].yao1s, res[s].yao8n, res[s].yao8s);

    printf("\nsqrmod / mulmod, ns per number (GMP: mpz_mul + mpz_tdiv_r)\n");
    printf(" bits |");
    for (k = 0; k < NENG; k++)
        printf(" %15s |", ename[k]);
    printf(" fastest sqr\n");
    for (s = 0; s < ns; s++)
    {
        printf("%5d |", res[s].bits);
        for (k = 0; k < NENG; k++)
            if (have_eng[k])
                printf(" %7.1f %7.1f |", res[s].sqr[k], res[s].mul[k]);
            else
                printf("     n/a     n/a |");
        printf(" %s\n", ename[best_vec(res[s].sqr)]);
    }

    printf("\npowm (full-size exponent), us per number; ratios are speedups\n");
    printf(" bits |");
    for (k = 0; k < NENG; k++)
        printf(" %9s", ename[k]);
    printf(" | tmpl/old  Yao8/tmpl  best/GMP  fastest\n");
    for (s = 0; s < ns; s++)
    {
        int b = best_vec(res[s].powm);
        printf("%5d |", res[s].bits);
        for (k = 0; k < NENG; k++)
            if (have_eng[k])
                printf(" %9.1f", res[s].powm[k]);
            else
                printf("       n/a");
        printf(" | %7.2fx %9.2fx %8.2fx  %s\n",
            have_eng[E_AEOLD] ? res[s].powm[E_AEOLD] / res[s].powm[E_AETMPL] : 0.0,
            res[s].powm[E_AETMPL] / res[s].powm[E_YAO8],
            res[s].powm[E_GMP] / res[s].powm[b], ename[b]);
    }

    printf("\n%s: %d failures\n", fails ? "FAILED" : "PASSED", fails);
    gmp_randclear(rs);
    return fails != 0;
}
