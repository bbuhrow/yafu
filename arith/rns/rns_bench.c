/*
 * Benchmark: vector RNS (rns_avx512.c) vs AVX-ECM vecarith52 vs GMP.
 *   rns_bench [seed]
 * Both vector engines process 8 numbers modulo the same N; times are
 * reported per number.  Modexp uses the same sliding-window rule for both
 * vector engines (AVX-ECM's vecslidingmodexp), with one full-size exponent
 * shared by all lanes, and is checked against mpz_powm in every lane.
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

// pointers vec_common.c declares; set by pick_avxecm (as assign_mul_ptrs)
extern void (*vecmulmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
extern void (*vecsqrmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_monty_t *);

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// ---------------------------------------------------------------------------
// AVX-ECM glue (mirrors test_vecarith / assign_mul_ptrs in avx_ecm_main.c)
// ---------------------------------------------------------------------------

typedef struct
{
    int words;
    vec_monty_t *mdata;
    vec_bignum_t *n, *s;
    mpz_t R, Rinv;
    const char *core;
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

static void pick_avxecm(avxecm_t *ae)
{
    int maxbits = ae->words * DIGITBITS;
    if (maxbits == 416)
    {
        vecmulmod_ptr = &vecmulmod52_fixed416_bfips;
        vecsqrmod_ptr = &vecsqrmod52_fixed416_bfips;
        ae->core = "fixed416";
    }
    else if (maxbits == 624)
    {
        vecmulmod_ptr = &vecmulmod52_fixed624_bfips;
        vecsqrmod_ptr = &vecsqrmod52_fixed624_bfips;
        ae->core = "fixed624";
    }
    else if (maxbits == 832)
    {
        vecmulmod_ptr = &vecmulmod52_fixed832_bfips;
        vecsqrmod_ptr = &vecsqrmod52_fixed832_bfips;
        ae->core = "fixed832";
    }
    else if (maxbits == 1040)
    {
        vecmulmod_ptr = &vecmulmod52_fixed1040_bfips;
        vecsqrmod_ptr = &vecsqrmod52_fixed1040_bfips;
        ae->core = "fixed1040";
    }
    else
    {
        vecmulmod_ptr = &vecmulmod52;
        vecsqrmod_ptr = &vecsqrmod52;
        ae->core = "bfips";
    }
}

static void avxecm_init(avxecm_t *ae, mpz_t N)
{
    vec_monty_t *md;
    int words = 0, j;

    while (words * DIGITBITS <= (int)mpz_sizeinbase(N, 2))
        words += 4;
    ae->words = words;
    md = ae->mdata = vec_monty_alloc(words);
    md->NBLOCKS = words / BLOCKWORDS;
    md->MAXBITS = words * DIGITBITS;
    md->NWORDS = words;
    md->isMersenne = 0;
    md->nbits = mpz_sizeinbase(N, 2);
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
    pick_avxecm(ae);
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
    vecsqrmod_ptr(b, g[2], ae->n, ae->s, md);
    for (i = 1; i < (1 << (k - 1)); i++)
        vecmulmod_ptr(g[2 * i - 1], g[2], g[2 * i + 1], ae->n, ae->s, md);

    vecCopyn(md->one, d, md->NWORDS);
    while (bit >= 0)
    {
        if (mpz_tstbit(e, bit) == 0)
        {
            vecsqrmod_ptr(d, d, ae->n, ae->s, md);
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
                vecsqrmod_ptr(d, d, ae->n, ae->s, md);
            vecmulmod_ptr(d, g[bstr], d, ae->n, ae->s, md);
            bit -= thisk;
        }
    }
}

// ---------------------------------------------------------------------------

static rns_ctx_t *rns_setup(rns_fixed_t *fx, mpz_t N)
{
    return rns_ctx_init(fx, N);
}

static void random_N(mpz_t N, int bits, gmp_randstate_t rs, rns_fixed_t *fx,
    rns_ctx_t **ctx)
{
    do
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        *ctx = rns_setup(fx, N);
    } while (*ctx == NULL);
}

static int run_size(int bits, gmp_randstate_t rs)
{
    rns_fixed_t *fx = rns_fixed_init(bits);
    rns_ctx_t *ctx;
    avxecm_t ae;
    mpz_t N, e, t, a[VECLEN], r[VECLEN];
    uint64_t *vx, *vd, *tmp;
    vec_bignum_t *ax, *ad;
    double t0, el, sq_rns, sq_ae, mul_rns, mul_ae, pw_rns, pw_ae, pw_gmp;
    double sq_y, mul_y, pw_y, sq_b, mul_b, pw_b;
    rnsy_fixed_t *yfx = rnsy_fixed_init(bits);
    rnsy_fixed_t *bfx = rnsy_fixed_init_n(bits, rnsyb_choose_n(bits));
    rnsy_ctx_t *yctx, *bctx;
    uint64_t *yx, *yd, *bx, *bd, *bt;
    int i, l, reps, fails = 0;

    mpz_init(N); mpz_init(e); mpz_init(t);
    for (l = 0; l < VECLEN; l++)
    {
        mpz_init(a[l]);
        mpz_init(r[l]);
    }
    for (;;)
    {
        random_N(N, bits, rs, fx, &ctx);
        if ((yctx = rnsy_ctx_init(yfx, N)) != NULL)
        {
            if ((bctx = rnsy_ctx_init(bfx, N)) != NULL)
                break;
            rnsy_ctx_free(yctx);
        }
        rns_ctx_free(ctx);
    }
    avxecm_init(&ae, N);
    yx = rnsy_alloc(yctx);
    yd = rnsy_alloc(yctx);
    bx = rnsyb_alloc(bctx);
    bd = rnsyb_alloc(bctx);
    bt = rnsy_alloc(bctx);

    vx = rns_valloc(ctx);
    vd = rns_valloc(ctx);
    tmp = rns_alloc(ctx);
    ax = vecInit(ae.words);
    ad = vecInit(ae.words);

    for (l = 0; l < VECLEN; l++)
    {
        mpz_urandomm(a[l], rs, N);
        rns_to_monty(ctx, tmp, a[l]);
        rns_vset_lane(ctx, vx, l, tmp);
        rnsy_to_monty(bctx, bt, a[l]);
        rnsyb_vset_lane(bctx, bx, l, bt);
        mpz_mul_2exp(t, a[l], ae.words * DIGITBITS);
        mpz_mod(t, t, N);
        vec_insert(ax, t, l);
    }

    // sqrmod / mulmod throughput, ns per number.  Warm up first so the
    // first timed loop doesn't run while the core is still clocking up.
    reps = 20000;
    memcpy(vd, vx, (size_t)fx->nch * VECLEN * sizeof(uint64_t));
    t0 = now_sec();
    while (now_sec() - t0 < 0.2)
        rns_vsqrmod_avx512(ctx, vd, vd);
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        rns_vsqrmod_avx512(ctx, vd, vd);
    sq_rns = (now_sec() - t0) / reps / VECLEN * 1e9;
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        rns_vmulmod_avx512(ctx, vd, vd, vx);
    mul_rns = (now_sec() - t0) / reps / VECLEN * 1e9;

    // Yao: one number per call
    rnsy_to_monty(yctx, yx, a[0]);
    memcpy(yd, yx, (size_t)yfx->nwords * sizeof(uint64_t));
    t0 = now_sec();
    while (now_sec() - t0 < 0.2)
    {
        rnsy_vsqrmod_avx512(yctx, yd, yd);
        rnsy_vmulmod_avx512(yctx, yd, yd, yx);
    }
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        rnsy_vsqrmod_avx512(yctx, yd, yd);
    sq_y = (now_sec() - t0) / reps * 1e9;
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        rnsy_vmulmod_avx512(yctx, yd, yd, yx);
    mul_y = (now_sec() - t0) / reps * 1e9;

    // Yao batch: 8 numbers per call
    memcpy(bd, bx, (size_t)bfx->nwords * VECLEN * sizeof(uint64_t));
    t0 = now_sec();
    while (now_sec() - t0 < 0.2)
    {
        rnsyb_vsqrmod_avx512(bctx, bd, bd);
        rnsyb_vmulmod_avx512(bctx, bd, bd, bx);
    }
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        rnsyb_vsqrmod_avx512(bctx, bd, bd);
    sq_b = (now_sec() - t0) / reps / VECLEN * 1e9;
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        rnsyb_vmulmod_avx512(bctx, bd, bd, bx);
    mul_b = (now_sec() - t0) / reps / VECLEN * 1e9;

    vecCopyn(ax, ad, ae.words);
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        vecsqrmod_ptr(ad, ad, ae.n, ae.s, ae.mdata);
    sq_ae = (now_sec() - t0) / reps / VECLEN * 1e9;
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        vecmulmod_ptr(ad, ax, ad, ae.n, ae.s, ae.mdata);
    mul_ae = (now_sec() - t0) / reps / VECLEN * 1e9;

    // modexp with a full-size exponent shared by all lanes
    mpz_urandomb(e, rs, bits);
    mpz_setbit(e, bits - 1);

    t0 = now_sec();
    for (reps = 0; (el = now_sec() - t0) < 0.5 || reps < 2; reps++)
        rns_vpowm_avx512(ctx, vd, vx, e);
    pw_rns = el / reps / VECLEN * 1e6;

    // Yao powm: 8 numbers, one call each, to match the per-number accounting
    t0 = now_sec();
    for (reps = 0; (el = now_sec() - t0) < 0.5 || reps < 2; reps++)
        for (l = 0; l < VECLEN; l++)
            rnsy_vpowm_avx512(yctx, yd, yx, e);
    pw_y = el / reps / VECLEN * 1e6;

    t0 = now_sec();
    for (reps = 0; (el = now_sec() - t0) < 0.5 || reps < 2; reps++)
        rnsyb_vpowm_avx512(bctx, bd, bx, e);
    pw_b = el / reps / VECLEN * 1e6;

    t0 = now_sec();
    for (reps = 0; (el = now_sec() - t0) < 0.5 || reps < 2; reps++)
        avxecm_powm(&ae, ad, ax, e);
    pw_ae = el / reps / VECLEN * 1e6;

    t0 = now_sec();
    for (reps = 0; (el = now_sec() - t0) < 0.5 || reps < 2; reps++)
        for (l = 0; l < VECLEN; l++)
            mpz_powm(r[l], a[l], e, N);
    pw_gmp = el / reps / VECLEN * 1e6;

    // check every lane of both engines against mpz_powm
    for (l = 0; l < VECLEN; l++)
    {
        rns_vget_lane(ctx, tmp, vd, l);
        rns_from_monty(ctx, t, tmp);
        if (mpz_cmp(t, r[l]) != 0)
        {
            printf("FAIL rns powm %d bits lane %d\n", bits, l);
            fails++;
        }
        vec_extract(t, ad, l, ae.words);
        mpz_mul(t, t, ae.Rinv);
        mpz_mod(t, t, N);
        if (mpz_cmp(t, r[l]) != 0)
        {
            printf("FAIL avx-ecm powm %d bits lane %d\n", bits, l);
            fails++;
        }
    }

    // Yao powm (lane 0's base) and every batch lane against mpz_powm
    rnsy_from_monty(yctx, t, yd);
    if (mpz_cmp(t, r[0]) != 0)
    {
        printf("FAIL yao powm %d bits\n", bits);
        fails++;
    }
    for (l = 0; l < VECLEN; l++)
    {
        rnsyb_vget_lane(bctx, bt, bd, l);
        rnsy_from_monty(bctx, t, bt);
        if (mpz_cmp(t, r[l]) != 0)
        {
            printf("FAIL yao batch powm %d bits lane %d\n", bits, l);
            fails++;
        }
    }

    printf("%5d | %3d+%-3d %2d %2d/%-2d %4d %-9s | %6.1f %6.1f %6.1f %6.1f | %6.1f %6.1f %6.1f %6.1f |"
        " %7.1f %7.1f %7.1f %7.1f %8.1f | %5.2f %5.2f %5.2f\n",
        bits, fx->nB, fx->nC, yfx->s, bfx->n, bfx->s, ae.words, ae.core,
        sq_rns, sq_y, sq_b, sq_ae, mul_rns, mul_y, mul_b, mul_ae,
        pw_rns, pw_y, pw_b, pw_ae, pw_gmp, pw_ae / pw_y, pw_ae / pw_b, pw_gmp / pw_b);
    fflush(stdout);
    rns_vfree(yx); rns_vfree(yd);
    rns_vfree(bx); rns_vfree(bd); rns_vfree(bt);
    rnsy_ctx_free(yctx);
    rnsy_ctx_free(bctx);
    rnsy_fixed_free(yfx);
    rnsy_fixed_free(bfx);

    for (l = 0; l < VECLEN; l++)
    {
        mpz_clear(a[l]);
        mpz_clear(r[l]);
    }
    rns_vfree(vx); rns_vfree(vd); free(tmp);
    vecFree(ax); vecFree(ad);
    avxecm_free(&ae);
    rns_ctx_free(ctx);
    rns_fixed_free(fx);
    mpz_clear(N); mpz_clear(e); mpz_clear(t);
    return fails;
}

int main(int argc, char **argv)
{
    int sizes[] = { 256, 400, 512, 600, 800, 1000, 1024, 1536, 2048, 3072, 4096 };
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : (unsigned long)time(NULL);
    gmp_randstate_t rs;
    int s, fails = 0;

    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, seed);
    printf("seed %lu; times per number (vector times / 8)\n", seed);
    printf("Yao1 = one number per call (n = 8); Yao8 = batch, 8 numbers per call (n/s from the cost model)\n");
    printf("      | small Yao1 Yao8  AE             |  sqrmod ns                  |  mulmod ns                  |"
        "  powm us                                 | AE/Y1 AE/Y8 GMP/Y8\n");
    printf(" bits | nB+nC  s   n/s  words core      |  small   Yao1   Yao8 AVXECM |  small   Yao1   Yao8 AVXECM |"
        "   small    Yao1    Yao8  AVXECM      GMP |\n");
    for (s = 0; s < (int)(sizeof(sizes) / sizeof(sizes[0])); s++)
        fails += run_size(sizes[s], rs);
    printf("%s: %d failures\n", fails ? "FAILED" : "PASSED", fails);
    gmp_randclear(rs);
    return fails != 0;
}
