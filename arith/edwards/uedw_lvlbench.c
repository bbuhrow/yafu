/* Model calibration for the B1 / B2 study: one level of microecm_edw
   (uecm_edw_level, up to 48 curves) on n = p q with p of N bits and q of
   63 - N bits (too large to be found at these levels), per parameter set
   and B1: time per factor found, mean curves, fraction found.
   usage: uedw_lvlbench N num set:B1 [set:B1 ...]
   gcc -O2 -march=native -I../../factor/gmp-ecm -I$G/include uedw_lvlbench.c \
       ../../factor/gmp-ecm/microecm_edw.c $G/lib/libgmp.a */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "gmp.h"
#include "microecm_edw.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static gmp_randstate_t grs;

static void rprime(mpz_t p, int bits)
{
    do {
        mpz_urandomb(p, grs, bits);
        mpz_setbit(p, bits - 1);
        mpz_nextprime(p, p);
    } while ((int)mpz_sizeinbase(p, 2) != bits);
}

int main(int argc, char **argv)
{
    int N = atoi(argv[1]), num = atoi(argv[2]);
    uint64_t *q = malloc(num * sizeof(uint64_t)), *pp = malloc(num * sizeof(uint64_t));
    mpz_t p, r, n;
    int i, a, rep;
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, 7 + N);
    mpz_inits(p, r, n, NULL);
    for (i = 0; i < num; i++) {
        rprime(p, N);
        rprime(r, 63 - N);
        mpz_mul(n, p, r);
        q[i] = mpz_get_ui(n);
        pp[i] = mpz_get_ui(p);
    }
    for (a = 3; a < argc; a++) {
        int set = atoi(argv[a]);
        uint32_t B1 = (uint32_t)atoi(argv[a] + 2);
        double best = 1e30;
        long tc = 0, found = 0;
        uecm_edw_set_params(set);
        for (rep = 0; rep < 3; rep++) {
            uint64_t lcg = 12345 + rep;
            long c = 0, fd = 0;
            double t0 = now();
            for (i = 0; i < num; i++) {
                int nc;
                uint64_t f = uecm_edw_level(q[i], B1, 48, &lcg, &nc);
                c += nc;
                fd += (f == pp[i]);
            }
            double t1 = now();
            if (t1 - t0 < best) best = t1 - t0;
            tc = c;
            found = fd;
        }
        printf("N %d set %d B1 %3u: %8.0f ns per find  %6.0f ns per input  curves %.2f  found %.4f\n",
            N, set, B1, best * 1e9 / (found ? found : 1), best * 1e9 / num, (double)tc / num,
            (double)found / num);
        fflush(stdout);
    }
    return 0;
}
