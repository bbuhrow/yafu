/*
 * Test the scalar RNS reference against GMP.
 *   rns_test [iterations] [seed]
 */

#include "rns.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int check(rns_ctx_t *ctx, uint64_t *x, mpz_t expect, const char *what,
    int bits, int iter)
{
    mpz_t a, raw, bound;
    int ok;

    mpz_init(a);
    mpz_init(raw);
    mpz_init(bound);
    rns_from_monty(ctx, a, x);
    ok = (mpz_cmp(a, expect) == 0);
    if (!ok)
        gmp_printf("FAIL %s: %d bits, iter %d\n  got %Zx\n  exp %Zx\n",
            what, bits, iter, a, expect);
    mpz_clear(a);
    mpz_clear(raw);
    mpz_clear(bound);
    return ok;
}

// the represented integer of a mulmod output must be < k*N
static int check_bound(rns_ctx_t *ctx, uint64_t *x, const char *what, int bits)
{
    mpz_t raw, bound;
    int ok;

    mpz_init(raw);
    mpz_init(bound);
    rns_to_int(ctx, raw, x);
    mpz_mul_ui(bound, ctx->N, ctx->fx->k);
    ok = (mpz_cmp(raw, bound) < 0);
    if (!ok)
        printf("FAIL bound %s: %d bits, output >= k*N\n", what, bits);
    mpz_clear(raw);
    mpz_clear(bound);
    return ok;
}

static int test_size(int bits, int iters, gmp_randstate_t rs)
{
    rns_fixed_t *fx = rns_fixed_init(bits);
    rns_ctx_t *ctx;
    mpz_t N, a, b, e, t;
    uint64_t *x, *y, *z, *w;
    int it, fails = 0;

    mpz_init(N); mpz_init(a); mpz_init(b); mpz_init(e); mpz_init(t);

    do
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        ctx = rns_ctx_init(fx, N);
    } while (ctx == NULL);

    printf("%5d bits: nB = %d, nC = %d, k = %d, d_max = %llu\n", bits,
        fx->nB, fx->nC, fx->k,
        (unsigned long long)(fx->db[fx->nB - 1] > fx->dc[fx->nC - 1] ?
            fx->db[fx->nB - 1] : fx->dc[fx->nC - 1]));

    x = rns_alloc(ctx); y = rns_alloc(ctx); z = rns_alloc(ctx); w = rns_alloc(ctx);

    for (it = 0; it < iters; it++)
    {
        // random mul / sqr / add / sub, compared against GMP
        mpz_urandomm(a, rs, N);
        mpz_urandomm(b, rs, N);
        rns_to_monty(ctx, x, a);
        rns_to_monty(ctx, y, b);

        rns_mulmod_ref(ctx, z, x, y);
        mpz_mul(e, a, b); mpz_mod(e, e, N);
        fails += !check(ctx, z, e, "mul", bits, it);
        fails += !check_bound(ctx, z, "mul", bits);

        rns_sqrmod_ref(ctx, w, z);
        mpz_mul(e, e, e); mpz_mod(e, e, N);
        fails += !check(ctx, w, e, "sqr", bits, it);

        // (z+w)*(z-w) with z, w < kN: sub input < 2kN, mul inputs < 4kN
        rns_addmod_ref(ctx, x, z, w);
        rns_submod_ref(ctx, y, z, w);
        rns_mulmod_ref(ctx, z, x, y);
        mpz_mul(t, a, b); mpz_mod(t, t, N);     // z
        mpz_mul(b, t, t); mpz_mod(b, b, N);     // w
        mpz_add(a, t, b);
        mpz_sub(t, t, b);
        mpz_mul(e, a, t); mpz_mod(e, e, N);
        fails += !check(ctx, z, e, "add/sub/mul", bits, it);
        fails += !check_bound(ctx, z, "add/sub/mul", bits);
    }

    // worst-case inputs: X = Y = h*k*N - 1, raw (not reduced mod N)
    mpz_mul_ui(t, N, (unsigned long)RNS_HEADROOM * fx->k);
    mpz_sub_ui(t, t, 1);
    rns_from_int(ctx, x, t);
    rns_mulmod_ref(ctx, z, x, x);
    mpz_mul(e, t, t);
    mpz_mul(e, e, ctx->MBinvmodN);  // raw*raw*M_B^{-1}; check via from_monty
    mpz_mul(e, e, ctx->MBinvmodN);
    mpz_mod(e, e, N);
    fails += !check(ctx, z, e, "worst-case", bits, 0);
    fails += !check_bound(ctx, z, "worst-case", bits);

    // long squaring chain: a^(2^iters)
    mpz_urandomm(a, rs, N);
    rns_to_monty(ctx, x, a);
    mpz_set(e, a);
    for (it = 0; it < iters * 10; it++)
    {
        rns_sqrmod_ref(ctx, x, x);
        mpz_mul(e, e, e);
        mpz_mod(e, e, N);
    }
    fails += !check(ctx, x, e, "sqr chain", bits, it);

    free(x); free(y); free(z); free(w);
    rns_ctx_free(ctx);
    rns_fixed_free(fx);
    mpz_clear(N); mpz_clear(a); mpz_clear(b); mpz_clear(e); mpz_clear(t);
    return fails;
}

#ifdef RNS_HAVE_AVX512
#define HAVE_VEC 1

// compare all 8 lanes of v against scalar residue vectors s[lane]
static int vcmp(rns_ctx_t *ctx, uint64_t *v, uint64_t **s, uint64_t *tmp,
    const char *what, int bits, int iter)
{
    int lane, ch, fails = 0;
    for (lane = 0; lane < RNS_VL; lane++)
    {
        rns_vget_lane(ctx, tmp, v, lane);
        for (ch = 0; ch < ctx->fx->nch; ch++)
        {
            if (tmp[ch] != s[lane][ch])
            {
                printf("FAIL vec %s: %d bits, iter %d, lane %d, ch %d: "
                    "%013llx != %013llx\n", what, bits, iter, lane, ch,
                    (unsigned long long)tmp[ch], (unsigned long long)s[lane][ch]);
                fails++;
                break;
            }
        }
    }
    return fails;
}

static int test_vec(int bits, int iters, gmp_randstate_t rs)
{
    rns_fixed_t *fx = rns_fixed_init(bits);
    rns_ctx_t *ctx;
    mpz_t N, a, e[RNS_VL];
    uint64_t *vx, *vy, *vz, *vw, *tmp;
    uint64_t *sx[RNS_VL], *sy[RNS_VL], *sz[RNS_VL], *sw[RNS_VL];
    int it, l, fails = 0;

    mpz_init(N); mpz_init(a);
    do
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        ctx = rns_ctx_init(fx, N);
    } while (ctx == NULL);

    vx = rns_valloc(ctx); vy = rns_valloc(ctx);
    vz = rns_valloc(ctx); vw = rns_valloc(ctx);
    tmp = rns_alloc(ctx);
    for (l = 0; l < RNS_VL; l++)
    {
        sx[l] = rns_alloc(ctx); sy[l] = rns_alloc(ctx);
        sz[l] = rns_alloc(ctx); sw[l] = rns_alloc(ctx);
        mpz_init(e[l]);
    }

    for (it = 0; it < iters; it++)
    {
        for (l = 0; l < RNS_VL; l++)
        {
            mpz_urandomm(a, rs, N);
            rns_to_monty(ctx, sx[l], a);
            mpz_urandomm(a, rs, N);
            rns_to_monty(ctx, sy[l], a);
            rns_vset_lane(ctx, vx, l, sx[l]);
            rns_vset_lane(ctx, vy, l, sy[l]);
        }

        rns_vmulmod_avx512(ctx, vz, vx, vy);
        for (l = 0; l < RNS_VL; l++) rns_mulmod_ref(ctx, sz[l], sx[l], sy[l]);
        fails += vcmp(ctx, vz, sz, tmp, "mul", bits, it);

        rns_vsqrmod_avx512(ctx, vw, vz);
        for (l = 0; l < RNS_VL; l++) rns_sqrmod_ref(ctx, sw[l], sz[l]);
        fails += vcmp(ctx, vw, sw, tmp, "sqr", bits, it);

        rns_vaddmod_avx512(ctx, vx, vz, vw);
        for (l = 0; l < RNS_VL; l++) rns_addmod_ref(ctx, sx[l], sz[l], sw[l]);
        fails += vcmp(ctx, vx, sx, tmp, "add", bits, it);

        rns_vsubmod_avx512(ctx, vy, vz, vw);
        for (l = 0; l < RNS_VL; l++) rns_submod_ref(ctx, sy[l], sz[l], sw[l]);
        fails += vcmp(ctx, vy, sy, tmp, "sub", bits, it);

        // in place: x = x*y, then y = y*y
        rns_vmulmod_avx512(ctx, vx, vx, vy);
        rns_vsqrmod_avx512(ctx, vy, vy);
        for (l = 0; l < RNS_VL; l++)
        {
            rns_mulmod_ref(ctx, sz[l], sx[l], sy[l]);
            rns_sqrmod_ref(ctx, sw[l], sy[l]);
        }
        fails += vcmp(ctx, vx, sz, tmp, "mul in place", bits, it);
        fails += vcmp(ctx, vy, sw, tmp, "sqr in place", bits, it);
    }

    // squaring chain checked against GMP in every lane
    for (l = 0; l < RNS_VL; l++)
    {
        mpz_urandomm(e[l], rs, N);
        rns_to_monty(ctx, sx[l], e[l]);
        rns_vset_lane(ctx, vx, l, sx[l]);
    }
    for (it = 0; it < iters * 10; it++)
        rns_vsqrmod_avx512(ctx, vx, vx);
    for (l = 0; l < RNS_VL; l++)
    {
        mpz_ui_pow_ui(a, 2, iters * 10);
        mpz_powm(e[l], e[l], a, N);
        rns_vget_lane(ctx, tmp, vx, l);
        fails += !check(ctx, tmp, e[l], "vec sqr chain", bits, l);
    }

    printf("%5d bits: vector tests done\n", bits);

    for (l = 0; l < RNS_VL; l++)
    {
        free(sx[l]); free(sy[l]); free(sz[l]); free(sw[l]);
        mpz_clear(e[l]);
    }
    rns_vfree(vx); rns_vfree(vy); rns_vfree(vz); rns_vfree(vw);
    free(tmp);
    rns_ctx_free(ctx);
    rns_fixed_free(fx);
    mpz_clear(N); mpz_clear(a);
    return fails;
}
#endif

static double now_sec(void)
{
    struct timespec ts;
#ifdef _WIN32
    timespec_get(&ts, TIME_UTC);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

// ns per squaring per number: scalar ref, vector (per lane), GMP mul+mod
static void bench(int bits, int reps, gmp_randstate_t rs)
{
    rns_fixed_t *fx = rns_fixed_init(bits);
    rns_ctx_t *ctx;
    mpz_t N, a, q;
    uint64_t *x;
    double t0, tref, tvec = 0, tgmp;
    int i;

    mpz_init(N); mpz_init(a); mpz_init(q);
    do
    {
        mpz_urandomb(N, rs, bits);
        mpz_setbit(N, bits - 1);
        mpz_setbit(N, 0);
        ctx = rns_ctx_init(fx, N);
    } while (ctx == NULL);

    mpz_urandomm(a, rs, N);
    x = rns_alloc(ctx);
    rns_to_monty(ctx, x, a);
    t0 = now_sec();
    for (i = 0; i < reps; i++)
        rns_sqrmod_ref(ctx, x, x);
    tref = (now_sec() - t0) / reps * 1e9;

#ifdef HAVE_VEC
    {
        uint64_t *vx = rns_valloc(ctx);
        int l;
        for (l = 0; l < RNS_VL; l++)
            rns_vset_lane(ctx, vx, l, x);
        t0 = now_sec();
        for (i = 0; i < reps; i++)
            rns_vsqrmod_avx512(ctx, vx, vx);
        tvec = (now_sec() - t0) / reps / RNS_VL * 1e9;
        rns_vfree(vx);
    }
#endif

    t0 = now_sec();
    for (i = 0; i < reps; i++)
    {
        mpz_mul(q, a, a);
        mpz_tdiv_r(a, q, N);
    }
    tgmp = (now_sec() - t0) / reps * 1e9;

    printf("%5d bits: sqrmod ns/number: ref %8.1f  avx512 %8.1f  gmp %8.1f\n",
        bits, tref, tvec, tgmp);

    free(x);
    rns_ctx_free(ctx);
    rns_fixed_free(fx);
    mpz_clear(N); mpz_clear(a); mpz_clear(q);
}

int main(int argc, char **argv)
{
    int sizes[] = { 64, 128, 256, 333, 512, 700, 1024, 1536, 2048, 3000, 4096 };
    int iters = (argc > 1) ? atoi(argv[1]) : 200;
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
        fails += test_vec(sizes[s], iters / 4 + 1, rs);
#endif

    printf("%s: %d failures\n", fails ? "FAILED" : "PASSED", fails);

    for (s = 0; s < (int)(sizeof(sizes) / sizeof(sizes[0])); s++)
        bench(sizes[s], 20000, rs);

    gmp_randclear(rs);
    return fails != 0;
}
