/* A/B of the uecm_edw P-1 pre-pass: getfactor_uecm_edw with
   uecm_edw_set_pm1(0) / (1), alternating, best of `rounds`, on the
   semiprimes files 48..64 bits.  usage: pm1_ab [num [rounds]] (yafu root)
   build: as pm1_eval (pm1_build.sh), source pm1_ab.c */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "microecm_edw.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main(int argc, char **argv)
{
    int num = argc > 1 ? atoi(argv[1]) : 100000;
    int rounds = argc > 2 ? atoi(argv[2]) : 5;
    static uint64_t comp[200000];
    static uint32_t f1[200000], f2[200000];
    int bits;
    printf("bits | off us | on us | gain\n");
    for (bits = 48; bits <= 64; bits += 2) {
        char fn[64];
        int n = 0, i, r, on;
        double best[2] = { 1e30, 1e30 };
        FILE *in;
        sprintf(fn, "semiprimes_%dbit.dat", bits);
        in = fopen(fn, "r");
        while (n < num && fscanf(in, "%lu,%u,%u", &comp[n], &f1[n], &f2[n]) == 3)
            n++;
        fclose(in);
        for (r = 0; r < rounds; r++) {
            for (on = 0; on < 2; on++) {
                uint64_t lcg = 42;
                int good = 0;
                uecm_edw_set_pm1(on);
                double t0 = now();
                for (i = 0; i < n; i++) {
                    uint64_t f = getfactor_uecm_edw(comp[i], 0, &lcg);
                    good += (f == f1[i] || f == f2[i]);
                }
                double t = 1e6 * (now() - t0) / n;
                if (t < best[on]) best[on] = t;
                if (good != n) printf("  (%d bits pm1=%d: %d of %d)\n", bits, on, good, n);
            }
        }
        printf("%4d | %6.2f | %5.2f | %+.1f%%\n", bits, best[0], best[1],
            100.0 * (best[0] / best[1] - 1));
    }
    return 0;
}
