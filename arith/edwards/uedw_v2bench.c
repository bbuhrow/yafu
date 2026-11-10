/* End-to-end check of the B1 / B2 study (plan doc Part I): getfactor_uecm_edw
   with the legacy parameter set (1) vs the re-tuned one (2), on balanced
   semiprimes n = p q of every even size 32..64 bits; time per input and the
   fraction factored.  The two sets run alternately (best of `rounds`).
   usage: uedw_v2bench [num [rounds [seed [unbalanced [bits0 bits1 [mode]]]]]]
   mode: 0 = getfactor_uecm_edw, 1 = the x8 list, 2 = the x4 list (link
   microecm_edw_x4.c for 1 / 2).
   unbalanced = 1: p of bits/2 - 4 bits (a smaller factor than the dispatch
   assumes).
   gcc -O2 -march=native -I../../factor/gmp-ecm -I$G/include uedw_v2bench.c \
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
    int num = argc > 1 ? atoi(argv[1]) : 4000;
    int rounds = argc > 2 ? atoi(argv[2]) : 3;
    int unbal = argc > 4 ? atoi(argv[4]) : 0;
    int mode = argc > 7 ? atoi(argv[7]) : 0;      // 0 scalar, 1 x8 list, 2 x4 list
    uint64_t *q = malloc(num * sizeof(uint64_t));
    uint64_t *fl = malloc(num * sizeof(uint64_t));
    mpz_t p, r, n;
    int bits, i, k, v;
    uint64_t lcg = 0x9e3779b97f4a7c15ULL;
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, argc > 3 ? atoi(argv[3]) : 1);
    mpz_inits(p, r, n, NULL);
    printf("bits | legacy us  found | v2 us  found | speedup\n");
    int b0 = argc > 5 ? atoi(argv[5]) : 32, b1 = argc > 6 ? atoi(argv[6]) : 64;
    for (bits = b0; bits <= b1; bits += 2) {
        int pb = unbal ? bits / 2 - 4 : bits / 2;
        for (i = 0; i < num; i++) {
            do {
                rprime(p, pb);
                rprime(r, bits - pb);
                mpz_mul(n, p, r);
            } while ((int)mpz_sizeinbase(n, 2) != bits);
            q[i] = mpz_get_ui(n);
        }
        double best[3] = { 1e30, 1e30, 1e30 };
        int ok[3] = { 0, 0, 0 };
        for (k = 0; k < rounds; k++) {
            for (v = 1; v <= 2; v++) {
                int okv = 0;
                uecm_edw_set_params(v);
                double t0 = now();
                if (mode == 0) {
                    for (i = 0; i < num; i++) {
                        uint64_t f = getfactor_uecm_edw(q[i], 0, &lcg);
                        okv += (f > 1 && f < q[i] && q[i] % f == 0);
                    }
                } else {
                    if (mode == 1)
                        getfactor_uecm_edw_x8_list(q, fl, num, &lcg);
                    else
                        getfactor_uecm_edw_x4_list(q, fl, num, &lcg);
                    for (i = 0; i < num; i++)
                        okv += (fl[i] > 1 && fl[i] < q[i] && q[i] % fl[i] == 0);
                }
                double t1 = now();
                if (t1 - t0 < best[v]) best[v] = t1 - t0;
                ok[v] = okv;
            }
        }
        printf("%4d | %8.2f %.4f | %8.2f %.4f | %.3f\n", bits, best[1] * 1e6 / num,
            (double)ok[1] / num, best[2] * 1e6 / num, (double)ok[2] / num, best[1] / best[2]);
        fflush(stdout);
    }
    return 0;
}
