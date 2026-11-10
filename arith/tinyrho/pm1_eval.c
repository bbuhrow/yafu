/* micropm1 (getfactor_upm1) on the semiprimes files: success rate and
   time per input for each B1 (B2 = 25 B1), and the break-even against
   microecm_edw: a P-1 pre-pass pays if  t_pm1 < p_pm1 * t_ecm.
   usage: pm1_eval [num]       (run from the yafu root)
   build: arith/tinyrho/pm1_build.sh */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "microecm.h"
#include "microecm_edw.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main(int argc, char **argv)
{
    int num = argc > 1 ? atoi(argv[1]) : 20000;
    static uint64_t comp[200000];
    static uint32_t f1[200000], f2[200000];
    static const uint32_t b1s[4] = { 100, 333, 666, 1000 };
    int bits;
    uecm_edw_set_pm1(0);                // ECM alone
    printf("bits | ecm us | B1: pm1 us  found  (gain = found*ecm - pm1, us)\n");
    for (bits = 48; bits <= 64; bits += 4) {
        char fn[64];
        int n = 0, i, k;
        FILE *in;
        sprintf(fn, "semiprimes_%dbit.dat", bits);
        in = fopen(fn, "r");
        while (n < num && fscanf(in, "%lu,%u,%u", &comp[n], &f1[n], &f2[n]) == 3)
            n++;
        fclose(in);
        uint64_t lcg = 42;
        double t0 = now();
        for (i = 0; i < n; i++)
            (void)getfactor_uecm_edw(comp[i], 0, &lcg);
        double tecm = 1e6 * (now() - t0) / n;
        printf("%4d | %6.2f |", bits, tecm);
        for (k = 0; k < 4; k++) {
            int found = 0;
            t0 = now();
            for (i = 0; i < n; i++) {
                uint64_t f = getfactor_upm1(comp[i], b1s[k]);
                found += (f == f1[i] || f == f2[i]);
            }
            double tp = 1e6 * (now() - t0) / n;
            double pf = (double)found / n;
            printf(" %u: %5.2f %5.3f (%+.2f) |", b1s[k], tp, pf, pf * tecm - tp);
        }
        printf("\n");
    }
    return 0;
}
