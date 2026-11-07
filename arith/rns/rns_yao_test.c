/*
 * Test the Yao-regime scalar reference against GMP.
 *   rns_yao_test [iterations] [seed]
 */

#include "rns_yao.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

static int check(rnsy_ctx_t *ctx, uint64_t *x, mpz_t expect, const char *what,
    int bits, int iter)
{
    mpz_t a, raw, bound;
    int ok;

    mpz_init(a); mpz_init(raw); mpz_init(bound);
    rnsy_from_monty(ctx, a, x);
    ok = (mpz_cmp(a, expect) == 0);
    if (!ok)
        gmp_printf("FAIL %s: %d bits, iter %d\n  got %Zx\n  exp %Zx\n",
            what, bits, iter, a, expect);
    // represented integer must stay < h*k*N (outputs of mulmod < k*N)
    rnsy_to_int(ctx, raw, x);
    mpz_mul_ui(bound, ctx->N, (unsigned long)RNS_HEADROOM * ctx->fx->k);
    if (mpz_cmp(raw, bound) >= 0)
    {
        printf("FAIL bound %s: %d bits, iter %d\n", what, bits, iter);
        ok = 0;
    }
    mpz_clear(a); mpz_clear(raw); mpz_clear(bound);
    return ok;
}

static int check_out(rnsy_ctx_t *ctx, uint64_t *x, const char *what, int bits)
{
    mpz_t raw, bound;
    int ok;
    mpz_init(raw); mpz_init(bound);
    rnsy_to_int(ctx, raw, x);
    mpz_mul_ui(bound, ctx->N, ctx->fx->k);
    ok = (mpz_cmp(raw, bound) < 0);
    if (!ok)
        printf("FAIL output bound %s: %d bits\n", what, bits);
    mpz_clear(raw); mpz_clear(bound);
    return ok;
}

static int test_size(int bits, int iters, gmp_randstate_t rs)
{
    rnsy_fixed_t *fx = rnsy_fixed_init(bits);
    rnsy_ctx_t *ctx;
    mpz_t N, a, b, e, t;
    uint64_t *x, *y, *z, *w;
    int it, fails = 0;

    mpz_init(N); mpz_init(a); mpz_init(b); mpz_init(e); mpz_init(t);
    do
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        ctx = rnsy_ctx_init(fx, N);
    } while (ctx == NULL);

    printf("%5d bits: s = %d (w = %d, capacity %d), db = %llu..%llu, "
        "dc = %llu..%llu, v = %d, v' = %d\n", bits, fx->s, fx->w, 416 * fx->s - 14,
        (unsigned long long)fx->db[0], (unsigned long long)fx->db[RNSY_NCH - 1],
        (unsigned long long)fx->dc[0], (unsigned long long)fx->dc[RNSY_NCH - 1],
        fx->vB, fx->vC);

    x = rnsy_alloc(ctx); y = rnsy_alloc(ctx); z = rnsy_alloc(ctx); w = rnsy_alloc(ctx);

    for (it = 0; it < iters; it++)
    {
        mpz_urandomm(a, rs, N);
        mpz_urandomm(b, rs, N);
        rnsy_to_monty(ctx, x, a);
        rnsy_to_monty(ctx, y, b);

        rnsy_mulmod_ref(ctx, z, x, y);
        mpz_mul(e, a, b); mpz_mod(e, e, N);
        fails += !check(ctx, z, e, "mul", bits, it);
        fails += !check_out(ctx, z, "mul", bits);

        rnsy_mulmod_ref(ctx, w, z, z);
        mpz_mul(t, e, e); mpz_mod(t, t, N);
        fails += !check(ctx, w, t, "sqr", bits, it);

        // (z+w)(z-w)
        rnsy_addmod_ref(ctx, x, z, w);
        rnsy_submod_ref(ctx, y, z, w);
        rnsy_mulmod_ref(ctx, z, x, y);
        mpz_add(a, e, t);
        mpz_sub(b, e, t);
        mpz_mul(e, a, b); mpz_mod(e, e, N);
        fails += !check(ctx, z, e, "add/sub/mul", bits, it);
        fails += !check_out(ctx, z, "add/sub/mul", bits);
    }

    // worst case: X = h*k*N - 1, raw
    mpz_mul_ui(t, N, (unsigned long)RNS_HEADROOM * fx->k);
    mpz_sub_ui(t, t, 1);
    rnsy_from_int(ctx, x, t);
    rnsy_mulmod_ref(ctx, z, x, x);
    mpz_mul(e, t, t);
    mpz_mul(e, e, ctx->MBinvmodN);
    mpz_mul(e, e, ctx->MBinvmodN);
    mpz_mod(e, e, N);
    fails += !check(ctx, z, e, "worst-case", bits, 0);
    fails += !check_out(ctx, z, "worst-case", bits);

    // squaring chain
    mpz_urandomm(a, rs, N);
    rnsy_to_monty(ctx, x, a);
    for (it = 0; it < iters * 5; it++)
        rnsy_mulmod_ref(ctx, x, x, x);
    mpz_ui_pow_ui(t, 2, iters * 5);
    mpz_powm(e, a, t, N);
    fails += !check(ctx, x, e, "sqr chain", bits, it);

    rns_vfree(x); rns_vfree(y); rns_vfree(z); rns_vfree(w);
    rnsy_ctx_free(ctx);
    rnsy_fixed_free(fx);
    mpz_clear(N); mpz_clear(a); mpz_clear(b); mpz_clear(e); mpz_clear(t);
    return fails;
}

#ifdef RNS_HAVE_AVX512
#define HAVE_VEC 1

// v is a valid representation: digits < 2^52, every channel encodes the same
// integer X (C channels pre-scaled by C_j^{-1}, r = X mod 2^52), X < bm*N
static int quiet = 0;

static int vvalid(rnsy_ctx_t *ctx, uint64_t *v, int bm, const char *what,
    int bits, int iter)
{
    rnsy_fixed_t *fx = ctx->fx;
    mpz_t X, t, u;
    int j, k, bad = 0;

    for (k = 0; k < fx->nwords; k++)
        if (v[k] >> 52)
        {
            if (!quiet) printf("FAIL vec %s: %d bits, iter %d: word %d not normalized\n",
                what, bits, iter, k);
            return 1;
        }

    mpz_init(X); mpz_init(t); mpz_init(u);
    rnsy_to_int(ctx, X, v);
    mpz_mul_ui(t, ctx->N, (unsigned long)bm);
    if (mpz_cmp(X, t) >= 0)
    {
        if (!quiet) printf("FAIL vec %s: %d bits, iter %d: value >= %d*N\n", what, bits, iter, bm);
        bad = 1;
    }
    for (j = 0; j < fx->n && !bad; j++)
    {
        mpz_set_ui(u, 0);
        for (k = fx->s - 1; k >= 0; k--)
        {
            mpz_mul_2exp(u, u, 52);
            mpz_add_ui(u, u, v[(fx->s + k) * fx->n + j]);
        }
        mpz_mod(u, u, fx->c[j]);
        mpz_mod(t, X, fx->c[j]);
        mpz_mul(t, t, fx->Cjinv_c[j]);
        mpz_mod(t, t, fx->c[j]);
        if (mpz_cmp(t, u) != 0)
        {
            if (!quiet) printf("FAIL vec %s: %d bits, iter %d: C channel %d inconsistent\n",
                what, bits, iter, j);
            bad = 1;
        }
    }
    mpz_fdiv_r_2exp(t, X, 52);
    if (!bad && mpz_get_ui(t) != v[fx->rword])
    {
        if (!quiet) printf("FAIL vec %s: %d bits, iter %d: r channel inconsistent\n", what, bits, iter);
        bad = 1;
    }
    mpz_clear(X); mpz_clear(t); mpz_clear(u);
    return bad;
}

// exact kernels: word-for-word equal to the reference.  Column-fold
// kernels: a valid representation (bound bm*N) of the same value mod N.
static int vcmp(rnsy_ctx_t *ctx, uint64_t *v, uint64_t *r, int bm, const char *what,
    int bits, int iter)
{
    int k;
    if (rnsy_vec_exact())
    {
        for (k = 0; k < ctx->fx->nwords; k++)
        {
            if (v[k] != r[k])
            {
                printf("FAIL vec %s: %d bits, iter %d, word %d: %013llx != %013llx\n",
                    what, bits, iter, k, (unsigned long long)v[k], (unsigned long long)r[k]);
                return 1;
            }
        }
        return 0;
    }
    else
    {
        mpz_t a, b;
        int bad = vvalid(ctx, v, bm, what, bits, iter);
        mpz_init(a); mpz_init(b);
        rnsy_from_monty(ctx, a, v);
        rnsy_from_monty(ctx, b, r);
        if (!bad && mpz_cmp(a, b) != 0)
        {
            printf("FAIL vec %s: %d bits, iter %d: wrong value mod N\n", what, bits, iter);
            bad = 1;
        }
        mpz_clear(a); mpz_clear(b);
        return bad;
    }
}

static int test_vec(int bits, int iters, gmp_randstate_t rs)
{
    rnsy_fixed_t *fx = rnsy_fixed_init(bits);
    rnsy_ctx_t *ctx;
    mpz_t N, a, e, t;
    uint64_t *x, *y, *z, *w, *rz, *rw, *rx, *ry;
    size_t bytes;
    int it, fails = 0;

    mpz_init(N); mpz_init(a); mpz_init(e); mpz_init(t);
    do
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        ctx = rnsy_ctx_init(fx, N);
    } while (ctx == NULL);
    bytes = (size_t)fx->nwords * sizeof(uint64_t);

    x = rnsy_alloc(ctx); y = rnsy_alloc(ctx); z = rnsy_alloc(ctx); w = rnsy_alloc(ctx);
    rx = rnsy_alloc(ctx); ry = rnsy_alloc(ctx); rz = rnsy_alloc(ctx); rw = rnsy_alloc(ctx);

    for (it = 0; it < iters; it++)
    {
        mpz_urandomm(a, rs, N);
        rnsy_to_monty(ctx, x, a);
        mpz_urandomm(a, rs, N);
        rnsy_to_monty(ctx, y, a);

        rnsy_vmulmod_avx512(ctx, z, x, y);
        rnsy_mulmod_ref(ctx, rz, x, y);
        fails += vcmp(ctx, z, rz, fx->k, "mul", bits, it);

        rnsy_vsqrmod_avx512(ctx, w, z);
        rnsy_mulmod_ref(ctx, rw, rz, rz);
        fails += vcmp(ctx, w, rw, fx->k, "sqr", bits, it);

        rnsy_vaddmod_avx512(ctx, x, z, w);
        rnsy_addmod_ref(ctx, rx, rz, rw);
        fails += vcmp(ctx, x, rx, 2 * fx->k, "add", bits, it);

        rnsy_vsubmod_avx512(ctx, y, z, w);
        rnsy_submod_ref(ctx, ry, rz, rw);
        fails += vcmp(ctx, y, ry, 3 * fx->k, "sub", bits, it);

        // in place, inputs from add/sub (up to 4kN)
        memcpy(rz, x, bytes);
        rnsy_vmulmod_avx512(ctx, x, x, y);
        rnsy_mulmod_ref(ctx, rw, rz, y);
        fails += vcmp(ctx, x, rw, fx->k, "mul in place", bits, it);
        memcpy(rz, y, bytes);
        rnsy_vsqrmod_avx512(ctx, y, y);
        rnsy_mulmod_ref(ctx, rw, rz, rz);
        fails += vcmp(ctx, y, rw, fx->k, "sqr in place", bits, it);
    }

    // the validity check must reject corrupted representations
    if (!rnsy_vec_exact())
    {
        quiet = 1;
        memcpy(rz, z, bytes);
        rz[(fx->s + 0) * RNSY_NCH + 3] ^= 1;        // a C channel
        if (!vvalid(ctx, rz, fx->k, "self-test C", bits, 0))
        {
            printf("FAIL self-test: corrupted C channel accepted\n");
            fails++;
        }
        memcpy(rz, z, bytes);
        rz[16 * fx->s] ^= 1;                        // r channel
        if (!vvalid(ctx, rz, fx->k, "self-test r", bits, 0))
        {
            printf("FAIL self-test: corrupted r channel accepted\n");
            fails++;
        }
        if (!vvalid(ctx, z, 0, "self-test bound", bits, 0))
        {
            printf("FAIL self-test: bound violation accepted\n");
            fails++;
        }
        quiet = 0;
    }

    // worst case raw inputs h*k*N - 1
    mpz_mul_ui(t, N, (unsigned long)RNS_HEADROOM * fx->k);
    mpz_sub_ui(t, t, 1);
    rnsy_from_int(ctx, x, t);
    rnsy_vsqrmod_avx512(ctx, z, x);
    rnsy_mulmod_ref(ctx, rz, x, x);
    fails += vcmp(ctx, z, rz, fx->k, "worst-case sqr", bits, 0);
    rnsy_vmulmod_avx512(ctx, z, x, x);
    fails += vcmp(ctx, z, rz, fx->k, "worst-case mul", bits, 0);

    // powm against mpz_powm, full-size exponent
    mpz_urandomm(a, rs, N);
    mpz_urandomb(e, rs, bits);
    mpz_setbit(e, bits - 1);
    rnsy_to_monty(ctx, x, a);
    rnsy_vpowm_avx512(ctx, z, x, e);
    mpz_powm(t, a, e, N);
    fails += !check(ctx, z, t, "vec powm", bits, 0);

    printf("%5d bits: vector tests done\n", bits);

    rns_vfree(x); rns_vfree(y); rns_vfree(z); rns_vfree(w);
    rns_vfree(rx); rns_vfree(ry); rns_vfree(rz); rns_vfree(rw);
    rnsy_ctx_free(ctx);
    rnsy_fixed_free(fx);
    mpz_clear(N); mpz_clear(a); mpz_clear(e); mpz_clear(t);
    return fails;
}
#endif

#ifdef HAVE_VEC
// batch kernel: 8 numbers per call, n channels per base.  Each lane is
// checked against the reference applied to that lane's own inputs.
static int test_batch(int bits, int n, int iters, gmp_randstate_t rs)
{
    rnsy_fixed_t *fx = rnsy_fixed_init_n(bits, n);
    rnsy_ctx_t *ctx;
    mpz_t N, a, e, t;
    uint64_t *vx, *vy, *vz, *vw, *lx[RNSY_NCH], *ly[RNSY_NCH], *lz, *rz[RNSY_NCH], *rw[RNSY_NCH];
    int it, l, fails = 0;

    if (fx == NULL)
    {
        printf("%5d bits, n = %2d: no setup (needs s > 10, or constants exceed 52 bits), skipped\n", bits, n);
        return 0;
    }
    mpz_init(N); mpz_init(a); mpz_init(e); mpz_init(t);
    do
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        ctx = rnsy_ctx_init(fx, N);
    } while (ctx == NULL);

    vx = rnsyb_alloc(ctx); vy = rnsyb_alloc(ctx); vz = rnsyb_alloc(ctx); vw = rnsyb_alloc(ctx);
    lz = rnsy_alloc(ctx);
    for (l = 0; l < RNSY_NCH; l++)
    {
        lx[l] = rnsy_alloc(ctx); ly[l] = rnsy_alloc(ctx);
        rz[l] = rnsy_alloc(ctx); rw[l] = rnsy_alloc(ctx);
    }

#define BCHECK(v, ref, bm, what) \
    for (l = 0; l < RNSY_NCH; l++) { \
        rnsyb_vget_lane(ctx, lz, v, l); \
        fails += vcmp(ctx, lz, ref[l], bm, what, bits, it); }

    for (it = 0; it < iters; it++)
    {
        for (l = 0; l < RNSY_NCH; l++)
        {
            mpz_urandomm(a, rs, N);
            rnsy_to_monty(ctx, lx[l], a);
            mpz_urandomm(a, rs, N);
            rnsy_to_monty(ctx, ly[l], a);
            rnsyb_vset_lane(ctx, vx, l, lx[l]);
            rnsyb_vset_lane(ctx, vy, l, ly[l]);
        }

        rnsyb_vmulmod_avx512(ctx, vz, vx, vy);
        for (l = 0; l < RNSY_NCH; l++) rnsy_mulmod_ref(ctx, rz[l], lx[l], ly[l]);
        BCHECK(vz, rz, fx->k, "batch mul");

        rnsyb_vsqrmod_avx512(ctx, vw, vz);
        for (l = 0; l < RNSY_NCH; l++) rnsy_mulmod_ref(ctx, rw[l], rz[l], rz[l]);
        BCHECK(vw, rw, fx->k, "batch sqr");

        rnsyb_vaddmod_avx512(ctx, vx, vz, vw);
        for (l = 0; l < RNSY_NCH; l++) rnsy_addmod_ref(ctx, lx[l], rz[l], rw[l]);
        BCHECK(vx, lx, 2 * fx->k, "batch add");

        rnsyb_vsubmod_avx512(ctx, vy, vz, vw);
        for (l = 0; l < RNSY_NCH; l++) rnsy_submod_ref(ctx, ly[l], rz[l], rw[l]);
        BCHECK(vy, ly, 3 * fx->k, "batch sub");

        // in place, inputs up to 3kN; reference fed this kernel's own inputs
        for (l = 0; l < RNSY_NCH; l++)
        {
            rnsyb_vget_lane(ctx, lx[l], vx, l);
            rnsyb_vget_lane(ctx, ly[l], vy, l);
            rnsy_mulmod_ref(ctx, rz[l], lx[l], ly[l]);
            rnsy_mulmod_ref(ctx, rw[l], ly[l], ly[l]);
        }
        rnsyb_vmulmod_avx512(ctx, vx, vx, vy);
        rnsyb_vsqrmod_avx512(ctx, vy, vy);
        BCHECK(vx, rz, fx->k, "batch mul in place");
        BCHECK(vy, rw, fx->k, "batch sqr in place");
    }

    // worst case raw inputs h*k*N - 1 in every lane
    it = 0;
    mpz_mul_ui(t, N, (unsigned long)RNS_HEADROOM * fx->k);
    mpz_sub_ui(t, t, 1);
    rnsy_from_int(ctx, lx[0], t);
    for (l = 0; l < RNSY_NCH; l++)
    {
        rnsyb_vset_lane(ctx, vx, l, lx[0]);
        rnsy_mulmod_ref(ctx, rz[l], lx[0], lx[0]);
    }
    rnsyb_vsqrmod_avx512(ctx, vz, vx);
    BCHECK(vz, rz, fx->k, "batch worst-case");

    // powm against mpz_powm, a different base per lane, shared exponent
    mpz_urandomb(e, rs, bits);
    mpz_setbit(e, bits - 1);
    for (l = 0; l < RNSY_NCH; l++)
    {
        mpz_urandomm(a, rs, N);
        rnsy_to_monty(ctx, lx[l], a);
        rnsyb_vset_lane(ctx, vx, l, lx[l]);
    }
    rnsyb_vpowm_avx512(ctx, vz, vx, e);
    for (l = 0; l < RNSY_NCH; l++)
    {
        rnsy_from_monty(ctx, a, lx[l]);
        mpz_powm(t, a, e, N);
        rnsyb_vget_lane(ctx, lz, vz, l);
        fails += !check(ctx, lz, t, "batch powm", bits, l);
    }
#undef BCHECK

    printf("%5d bits, n = %2d: s = %2d, k = %2d, v = %d (stored %d), v' = %d (stored %d)%s: "
        "batch tests done\n", bits, n, fx->s, fx->k, fx->vB, fx->vBt, fx->vC, fx->vCt,
        fx->shifted ? ", shifted" : "");

    rns_vfree(vx); rns_vfree(vy); rns_vfree(vz); rns_vfree(vw); rns_vfree(lz);
    for (l = 0; l < RNSY_NCH; l++)
    {
        rns_vfree(lx[l]); rns_vfree(ly[l]); rns_vfree(rz[l]); rns_vfree(rw[l]);
    }
    rnsy_ctx_free(ctx);
    rnsy_fixed_free(fx);
    mpz_clear(N); mpz_clear(a); mpz_clear(e); mpz_clear(t);
    return fails;
}
#endif

int main(int argc, char **argv)
{
    int sizes[] = { 200, 512, 818, 1024, 1234, 1536, 2048, 2066, 3000, 4096, 4146 };
    int iters = (argc > 1) ? atoi(argv[1]) : 100;
    unsigned long seed = (argc > 2) ? strtoul(argv[2], NULL, 10) : (unsigned long)time(NULL);
    gmp_randstate_t rs;
    int s, fails = 0;

    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, seed);
    printf("seed %lu, %d iterations per size\n", seed, iters);
    for (s = 0; s < (int)(sizeof(sizes) / sizeof(sizes[0])); s++)
        fails += test_size(sizes[s], iters, rs);
#ifdef HAVE_VEC
    for (s = 0; s < (int)(sizeof(sizes) / sizeof(sizes[0])); s++)
        fails += test_vec(sizes[s], iters, rs);
    {
        // batch kernel: the cost-model n, plus a spread of n (some with
        // shifted constants), each also exercising the general-n reference
        int bsizes[] = { 200, 512, 1024, 1536, 2048, 3072, 4096 };
        int ns[] = { 2, 4, 6, 8, 9, 10, 11 };
        int bi, ni;
        for (bi = 0; bi < (int)(sizeof(bsizes) / sizeof(bsizes[0])); bi++)
        {
            int nc = rnsyb_choose_n(bsizes[bi]);
            printf("%5d bits: cost model picks n = %d\n", bsizes[bi], nc);
            fails += test_batch(bsizes[bi], nc, iters / 4 + 1, rs);
            for (ni = 0; ni < (int)(sizeof(ns) / sizeof(ns[0])); ni++)
                if (ns[ni] != nc)
                    fails += test_batch(bsizes[bi], ns[ni], iters / 8 + 1, rs);
        }
    }
#endif
#ifdef HAVE_VEC
    printf("vector kernels: %s\n", rnsy_vec_exact() ? "exact (RNSY_COLUMN_FOLD=0)" :
        "column fold (RNSY_COLUMN_FOLD=1), checked by congruence/consistency/bounds");
#endif
    printf("%s: %d failures\n", fails ? "FAILED" : "PASSED", fails);
    gmp_randclear(rs);
    return fails != 0;
}
