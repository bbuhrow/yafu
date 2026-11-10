/*
 * microecm_edw: curves needed per input, and the effect of the small-constant
 * curve table (edw_tiny_curves.h) on speed.  Build microecm_edw.c with
 * -DUEDW_STATS and link microecm.c (for the getfactor_uecm reference) and GMP.
 *
 *   uedw_stats [num [seed]]
 * Per level (balanced semiprimes): the distribution of curves run per input
 * (the finding curve included), then us per input for getfactor_uecm and
 * getfactor_uecm_edw with 0 and with all table curves.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "gmp.h"
#include "microecm.h"
#include "microecm_edw.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static gmp_randstate_t grs;

static uint64_t rprime(int bits)
{
    mpz_t p;
    mpz_init(p);
    do {
        mpz_urandomb(p, grs, bits);
        mpz_setbit(p, bits - 1);
        mpz_nextprime(p, p);
    } while ((int)mpz_sizeinbase(p, 2) != bits);
    uint64_t r = mpz_get_ui(p);
    mpz_clear(p);
    return r;
}

int main(int argc, char **argv)
{
    int num = argc > 1 ? atoi(argv[1]) : 10000;
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, argc > 2 ? atoi(argv[2]) : 1);
    uint64_t *q = malloc(num * sizeof(uint64_t));
    uint64_t hist[65];
    static const int ub[7] = { 40, 44, 48, 52, 58, 62, 64 };
    static const int lv[7] = { 27, 47, 70, 85, 125, 165, 205 };
    uint64_t lcg = 7;

    printf("bits  B1 | curves per input: mean  P(1)  P(<=2) P(<=3) P(<=5) P(<=8) | us: uecm   edw(no table) edw(table) | speedup vs uecm, table on/off\n");
    for (int k = 0; k < 7; k++) {
        for (int i = 0; i < num; i++) {
            uint64_t n;
            do {
                n = rprime(ub[k] / 2) * rprime(ub[k] - ub[k] / 2);
            } while (n >> (ub[k] - 1) != 1);
            q[i] = n;
        }
        double tt[3];
        for (int mode = 0; mode < 3; mode++) {
            if (mode > 0) uecm_edw_set_table(mode == 1 ? 0 : 64);
            uecm_edw_hist(hist, 64);
            double t0 = now();
            for (int i = 0; i < num; i++) {
                uint64_t f = (mode == 0) ? getfactor_uecm(q[i], 0, &lcg)
                                         : getfactor_uecm_edw(q[i], 0, &lcg);
                if (!(f > 1 && f < q[i] && q[i] % f == 0) && mode == 2 && i < 3)
                    fprintf(stderr, "miss %llu\n", (unsigned long long)q[i]);
            }
            tt[mode] = (now() - t0) * 1e6 / num;
        }
        uecm_edw_hist(hist, 64);    // clear (counts from mode 2 only were wanted:
        // rerun mode 2 quickly to collect them cleanly)
        uecm_edw_set_table(64);
        for (int i = 0; i < num; i++) getfactor_uecm_edw(q[i], 0, &lcg);
        uecm_edw_hist(hist, 64);
        double tot = 0, mean = 0, cdf[9] = { 0 };
        for (int j = 0; j <= 64; j++) { tot += hist[j]; mean += j * (double)hist[j]; }
        for (int c = 1; c <= 8; c++) {
            double s = 0;
            for (int j = 0; j <= c; j++) s += hist[j];
            cdf[c] = s / tot;
        }
        printf("%4d %3d |  %5.2f  %.3f  %.3f  %.3f  %.3f  %.3f | %7.2f  %7.2f  %7.2f | %.3f %.3f\n",
            ub[k], lv[k], mean / tot, cdf[1], cdf[2], cdf[3], cdf[5], cdf[8],
            tt[0], tt[1], tt[2], tt[0] / tt[2], tt[0] / tt[1]);
    }
    return 0;
}
