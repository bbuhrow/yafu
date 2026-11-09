/*
 * Above 16400 bits (handoff 4.7): per-number mulmod / sqrmod timings of
 *   AE-kara   AVX-ECM's karatsuba REDC (veckmul_redc / vecksqr_redc)
 *   AE-tmplN  the constant-S template, loop form (no karatsuba)
 *   Yao8      the batch Yao RNS kernel, large-s loop form, schoolbook
 *             channel products
 *   Yao8-kK   the same with karatsuba channel products down to K words
 *   GMP       mpz_mul + mpz_tdiv_r (per number)
 * at 16k .. 64k bits.  Every engine's mul and sqr result is checked
 * against GMP in every lane.
 *
 *   big_bench [seed] [bits ...]
 */

#include "rns.h"
#include "rns_yao.h"
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

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

// best of 3 runs of at least 0.3 s each (after a 0.2 s warm-up), ns per call
#define TIME_NS(out, op) do {                                           \
        double _t0 = now_sec(), _best = 1e30;                           \
        int _r, _n;                                                     \
        while (now_sec() - _t0 < 0.2) { op; }                           \
        for (_r = 0; _r < 3; _r++) {                                    \
            _t0 = now_sec();                                            \
            for (_n = 0; now_sec() - _t0 < 0.3 || _n < 3; _n++) { op; } \
            double _el = (now_sec() - _t0) / _n * 1e9;                  \
            if (_el < _best) _best = _el;                               \
        }                                                               \
        (out) = _best;                                                  \
    } while (0)

typedef void (*ae_mul_t)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
typedef void (*ae_sqr_t)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *,
    vec_bignum_t *, vec_monty_t *);

typedef struct
{
    int words;
    vec_monty_t *md;
    vec_bignum_t *n, *s;
    mpz_t R, Rinv;
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

static void ae_init(ae_t *ae, mpz_t N, int words)
{
    vec_monty_t *md;
    int j;
    ae->words = words;
    md = ae->md = vec_monty_alloc(words);
    md->NBLOCKS = (words + BLOCKWORDS - 1) / BLOCKWORDS;
    md->MAXBITS = words * DIGITBITS;
    md->NWORDS = words;
    md->isMersenne = 0;
    md->nbits = (int)mpz_sizeinbase(N, 2);
    md->use_vnbits = 0;
    mpz_init(ae->R);
    mpz_init(ae->Rinv);
    mpz_setbit(ae->R, (mp_bitcnt_t)(words * DIGITBITS));
    mpz_invert(md->nhat, N, ae->R);
    mpz_sub(md->nhat, ae->R, md->nhat);
    ae->n = vecInit(words);
    ae->s = vecInit(words);
    for (j = 0; j < VECLEN; j++)
    {
        vec_insert(ae->n, N, j, words);
        vec_insert(md->n, N, j, words);
        vec_insert(md->vnhat, md->nhat, j, words);
        md->vrho[j] = mpz_get_ui(md->nhat) & VEC_MAXDIGIT;
    }
    mpz_invert(ae->Rinv, ae->R, N);
}

static void ae_free(ae_t *ae)
{
    vecFree(ae->n);
    vecFree(ae->s);
    vec_monty_free(ae->md);
    free(ae->md);
    mpz_clear(ae->R);
    mpz_clear(ae->Rinv);
}

enum { E_KARA, E_TMPL, E_YAO, E_YAOK1, E_YAOK2, E_YAOK3, E_GMP, NENG };
static const char *ename[NENG] = { "AE-kara", "AE-tmplN", "Yao8", "Yao8-k16", "Yao8-k24", "Yao8-k32", "GMP" };
static const int kara_k[NENG] = { 0, 0, 0, 16, 24, 32, 0 };

int main(int argc, char **argv)
{
    int defsizes[] = { 16384, 24576, 32768, 49152, 65536 };
    int *sizes = defsizes, nsizes = 5, si, e, l, fails = 0;
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : 1;
    double tsqr[16][NENG], tmul[16][NENG];
    int yn[16], ys[16], kw[16], tw[16];
    gmp_randstate_t rs;

    if (argc > 2)
    {
        nsizes = argc - 2;
        sizes = (int *)malloc((size_t)nsizes * sizeof(int));
        for (si = 0; si < nsizes; si++)
            sizes[si] = atoi(argv[si + 2]);
    }
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, seed);
    printf("seed %lu\n", seed);

    for (si = 0; si < nsizes; si++)
    {
        int bits = sizes[si];
        int nY = rnsyb_choose_n(bits);
        rnsy_fixed_t *fx;
        rnsy_ctx_t *ctx = NULL;
        mpz_t N, a[VECLEN], b[VECLEN], t, u;
        ae_t ak, at;
        vec_bignum_t *va[2], *vb[2], *vc[2];
        uint64_t *ya, *yb, *yc, *y1;
        int wt = (bits + AVXECM_SPARE_BITS + DIGITBITS - 1) / DIGITBITS;
        int wk = (wt + BLOCKWORDS - 1) / BLOCKWORDS * BLOCKWORDS;
        ae_mul_t tmul_fn;
        ae_sqr_t tsqr_fn;

        mpz_init(N); mpz_init(t); mpz_init(u);
        for (l = 0; l < VECLEN; l++) { mpz_init(b[l]); mpz_init(a[l]); }
        fx = rnsy_fixed_init_n(bits, nY);
        for (;;)
        {
            mpz_urandomb(N, rs, (mp_bitcnt_t)bits);
            mpz_setbit(N, (mp_bitcnt_t)(bits - 1));
            mpz_setbit(N, 0);
            if ((ctx = rnsy_ctx_init(fx, N)) != NULL)
                break;
        }
        yn[si] = fx->n; ys[si] = fx->s; kw[si] = wk; tw[si] = wt;
        printf("%6d bits: Yao n = %d, s = %d; karatsuba REDC %d words; template %d words\n",
            bits, fx->n, fx->s, wk, wt);
        fflush(stdout);

        ae_init(&ak, N, wk);
        ae_init(&at, N, wt);
        vecmulmod52_tmpl_get(wt, &tmul_fn, &tsqr_fn);
        for (e = 0; e < 2; e++)
        {
            int w = e ? wt : wk;
            va[e] = vecInit(w); vb[e] = vecInit(w); vc[e] = vecInit(w);
        }
        ya = rnsyb_alloc(ctx); yb = rnsyb_alloc(ctx); yc = rnsyb_alloc(ctx);
        y1 = rnsy_alloc(ctx);
        for (l = 0; l < VECLEN; l++)
        {
            mpz_urandomm(a[l], rs, N);
            mpz_urandomm(b[l], rs, N);
            vec_insert(va[0], a[l], l, wk); vec_insert(vb[0], b[l], l, wk);
            vec_insert(va[1], a[l], l, wt); vec_insert(vb[1], b[l], l, wt);
            rnsy_to_monty(ctx, y1, a[l]); rnsyb_vset_lane(ctx, ya, l, y1);
            rnsy_to_monty(ctx, y1, b[l]); rnsyb_vset_lane(ctx, yb, l, y1);
        }

        // correctness: one mul and one sqr per engine, every lane
        for (e = 0; e < NENG; e++)
        {
            int sq;
            if (e == E_GMP)
                continue;
            for (sq = 0; sq < 2; sq++)
            {
                if (e == E_KARA || e == E_TMPL)
                {
                    int k = (e == E_TMPL);
                    ae_t *ae = k ? &at : &ak;
                    if (k)
                    {
                        if (sq) tsqr_fn(va[1], vc[1], at.n, at.s, at.md);
                        else tmul_fn(va[1], vb[1], vc[1], at.n, at.s, at.md);
                    }
                    else
                    {
                        if (sq) vecksqr_redc(va[0], vc[0], ak.n, ak.s, ak.md);
                        else veckmul_redc(va[0], vb[0], vc[0], ak.n, ak.s, ak.md);
                    }
                    for (l = 0; l < VECLEN; l++)
                    {
                        // Montgomery: out == a*b/R mod N
                        vec_extract(t, vc[k], l, ae->words);
                        mpz_mod(t, t, N);
                        mpz_mul(u, a[l], sq ? a[l] : b[l]);
                        mpz_mul(u, u, ae->Rinv);
                        mpz_mod(u, u, N);
                        if (mpz_cmp(t, u)) { printf("FAIL %s %s lane %d\n", ename[e], sq ? "sqr" : "mul", l); fails++; }
                    }
                }
                else
                {
                    rnsyb_kara_min = kara_k[e];
                    if (sq) rnsyb_vsqrmod(ctx, yc, ya);
                    else rnsyb_vmulmod(ctx, yc, ya, yb);
                    for (l = 0; l < VECLEN; l++)
                    {
                        rnsyb_vget_lane(ctx, y1, yc, l);
                        rnsy_from_monty(ctx, t, y1);
                        mpz_mul(u, a[l], sq ? a[l] : b[l]);
                        mpz_mod(u, u, N);
                        if (mpz_cmp(t, u)) { printf("FAIL %s %s lane %d\n", ename[e], sq ? "sqr" : "mul", l); fails++; }
                    }
                }
            }
        }

        // timings, ns per number
        TIME_NS(tsqr[si][E_KARA], vecksqr_redc(vc[0], vc[0], ak.n, ak.s, ak.md));
        TIME_NS(tmul[si][E_KARA], veckmul_redc(vc[0], vb[0], vc[0], ak.n, ak.s, ak.md));
        TIME_NS(tsqr[si][E_TMPL], tsqr_fn(vc[1], vc[1], at.n, at.s, at.md));
        TIME_NS(tmul[si][E_TMPL], tmul_fn(vc[1], vb[1], vc[1], at.n, at.s, at.md));
        for (e = E_YAO; e <= E_YAOK3; e++)
        {
            rnsyb_kara_min = kara_k[e];
            memcpy(yc, ya, (size_t)fx->nwords * RNSY_NCH * sizeof(uint64_t));
            TIME_NS(tsqr[si][e], rnsyb_vsqrmod(ctx, yc, yc));
            TIME_NS(tmul[si][e], rnsyb_vmulmod(ctx, yc, yc, yb));
        }
        rnsyb_kara_min = 0;
        for (e = 0; e < E_GMP; e++)
        {
            tsqr[si][e] /= VECLEN;
            tmul[si][e] /= VECLEN;
        }
        mpz_set(t, a[0]);
        TIME_NS(tsqr[si][E_GMP], (mpz_mul(u, t, t), mpz_tdiv_r(t, u, N)));
        TIME_NS(tmul[si][E_GMP], (mpz_mul(u, t, b[1]), mpz_tdiv_r(t, u, N)));

        for (e = 0; e < 2; e++) { vecFree(va[e]); vecFree(vb[e]); vecFree(vc[e]); }
        rns_vfree(ya); rns_vfree(yb); rns_vfree(yc); rns_vfree(y1);
        ae_free(&ak); ae_free(&at);
        rnsy_ctx_free(ctx);
        rnsy_fixed_free(fx);
        for (l = 0; l < VECLEN; l++) { mpz_clear(a[l]); mpz_clear(b[l]); }
        mpz_clear(N); mpz_clear(t); mpz_clear(u);
    }

    printf("\nsqrmod / mulmod, us per number\n  bits |");
    for (e = 0; e < NENG; e++)
        printf(" %17s |", ename[e]);
    printf(" fastest sqr  (vs AE-kara)\n");
    for (si = 0; si < nsizes; si++)
    {
        int best = 0;
        printf(" %5d |", sizes[si]);
        for (e = 0; e < NENG; e++)
        {
            printf(" %8.1f %8.1f |", tsqr[si][e] / 1000, tmul[si][e] / 1000);
            if (tsqr[si][e] < tsqr[si][best]) best = e;
        }
        printf(" %s (%.2fx)\n", ename[best], tsqr[si][E_KARA] / tsqr[si][best]);
    }
    printf("\nconfigurations:\n");
    for (si = 0; si < nsizes; si++)
        printf(" %5d bits: Yao n = %d, s = %d; AE-kara %d words, AE-tmplN %d words\n",
            sizes[si], yn[si], ys[si], kw[si], tw[si]);
    printf("\n%s: %d failures\n", fails ? "FAILED" : "PASSED", fails);
    gmp_randclear(rs);
    return fails != 0;
}
