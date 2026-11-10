/*
 * Component timings of tinyecm_edw's scalar path (includes tinyecm_edw.c):
 * curve loads / builds, stage 1, stage 2 and the gcd check, per B1 level.
 *   tedw_parts [bits [reps]]
 */
#define _GNU_SOURCE
#include "../../factor/gmp-ecm/tinyecm_edw.c"
#include <stdio.h>
#include <time.h>

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main(int argc, char **argv)
{
    int bits = argc > 1 ? atoi(argv[1]) : 66;
    int reps = argc > 2 ? atoi(argv[2]) : 20000;
    gmp_randstate_t rs;
    mpz_t n, f;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 3);
    mpz_inits(n, f, NULL);
    mpz_urandomb(n, rs, (mp_bitcnt_t)bits);
    mpz_setbit(n, (mp_bitcnt_t)(bits - 1));
    mpz_setbit(n, 0);
    tedw_ctx ctx;
    tedw_ctx_init(&ctx, n);
    tedw_curve cv;
    tedw_wpt R;
    volatile uint64_t sink = 0;
    double t0, t1;

    t0 = now();
    for (int i = 0; i < reps; i++) { tedw_table_curve(i % 5, &cv, &ctx); sink += cv.P.X.lo; }
    t1 = now();
    printf("%d-bit n\ntable curve load   %8.0f ns\n", bits, (t1 - t0) * 1e9 / reps);
    t0 = now();
    for (int i = 0; i < reps; i++) { tedw_w_start(&R, 3 + (uint32_t)(i & 255), &ctx); sink += R.X.lo; }
    t1 = now();
    printf("GKL start [m0]Q    %8.0f ns (m0 in [3, 258])\n", (t1 - t0) * 1e9 / reps);
    tedw_w_start(&R, 100, &ctx);
    t0 = now();
    for (int i = 0; i < reps; i++) { tedw_gkl_curve(&R, &cv, &ctx); tedw_w_addq(&R, &ctx); sink += cv.P.X.lo; }
    t1 = now();
    printf("GKL curve + step   %8.0f ns\n", (t1 - t0) * 1e9 / reps);
    t0 = now();
    for (int i = 0; i < reps; i++) { tedw_check(cv.P.X, n, f); }
    t1 = now();
    printf("gcd check (mpz)    %8.0f ns\n", (t1 - t0) * 1e9 / reps);
    t0 = now();
    for (int i = 0; i < reps; i++) { R.X = tedw_mulredc(R.X, cv.P.Y, &ctx); }
    t1 = now();
    printf("mulredc (dependent)%8.1f ns\n", (t1 - t0) * 1e9 / reps);
    sink += R.X.lo;

    static const uint32_t lv[7] = { 27, 47, 70, 85, 125, 165, 205 };
    for (int k = 0; k < 7; k++) {
        edw_chain_t ch;
        const uint8_t *ba;
        int nb;
        tedw_pt P;
        edw_chain_for(lv[k], &ch);
        tedw_s2table(lv[k], &ba, &nb);
        int r2 = reps / 10;
        t0 = now();
        for (int i = 0; i < r2; i++) { tedw_stage1(&cv, &ch, &P, &ctx); sink += P.Z.lo; }
        t1 = now();
        double s1 = (t1 - t0) * 1e9 / r2;
        t0 = now();
        for (int i = 0; i < r2; i++) { tedw_u128 a = tedw_stage2(&P, &cv, lv[k], ba, nb, &ctx); sink += a.lo; }
        t1 = now();
        printf("B1 %3u: stage 1 %8.0f ns, stage 2 %8.0f ns\n", lv[k], s1, (t1 - t0) * 1e9 / r2);
    }
    if (sink == 42) printf(" ");
    return 0;
}
