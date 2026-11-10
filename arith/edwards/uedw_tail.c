/* Curves-per-input tail of one microecm_edw level on a semiprimes file:
   run up to `cap` curves (uecm_edw_level, no cap at the dispatch's 32/42)
   and print the mean, the fraction needing > k curves for several k, and
   the per-curve success rate fitted to the tail (k >= 10).
   usage: uedw_tail file B1 params [cap]      (run from the yafu root)
   gcc -O2 -march=native -DUEDW_NO_PM1 -I../../factor/gmp-ecm uedw_tail.c \
       ../../factor/gmp-ecm/microecm_edw.c -o uedw_tail */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "microecm_edw.h"

int main(int argc, char **argv)
{
    static uint64_t comp[200000];
    static uint32_t f1[200000], f2[200000];
    int B1 = atoi(argv[2]);
    int params = atoi(argv[3]);
    int cap = argc > 4 ? atoi(argv[4]) : 400;
    static int hist[1001];
    uint64_t lcg = 42;
    int i, num = 0, fail = 0;
    double sum = 0;
    FILE *in = fopen(argv[1], "r");
    while (num < 100000 && fscanf(in, "%lu,%u,%u", &comp[num], &f1[num], &f2[num]) == 3)
        num++;
    fclose(in);
    uecm_edw_set_params(params);
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (i = 0; i < num; i++) {
        int nc = 0;
        uint64_t g = uecm_edw_level(comp[i], B1, cap, &lcg, &nc);
        if (!(g == f1[i] || g == f2[i])) { fail++; nc = cap + 1; }
        hist[nc > 1000 ? 1000 : nc]++;
        sum += nc;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double us = 1e6 * ((t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec));
    printf("%s B1=%d params=%d: mean curves %.3f, %.2f us/input, %.3f us/curve, fails at cap %d: %d\n",
        argv[1], B1, params, sum / num, us / num, us / sum, cap, fail);
    int ks[] = { 10, 20, 32, 42, 50, 64, 80, 100 };
    for (int j = 0; j < 8; j++) {
        int k = ks[j], c = 0;
        for (int m = k + 1; m <= 1000; m++) c += hist[m];
        printf("  P(>%d) = %.2e (%d)\n", k, (double)c / num, c);
    }
    // geometric fit of the tail: rate = (inputs past 10) / (curves past 10)
    double past = 0, curv = 0;
    for (int m = 11; m <= 1000; m++) { past += hist[m]; curv += (double)hist[m] * (m - 10); }
    printf("  tail per-curve success ~ %.4f\n", past / curv);
    return 0;
}
