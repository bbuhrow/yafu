/* curves run per input: getfactor_tecm (counted with a wrapper) vs getfactor_tecm_edw */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
static double now(void){struct timespec s;clock_gettime(CLOCK_MONOTONIC,&s);return s.tv_sec+1e-9*s.tv_nsec;}
#include "gmp.h"
#include "tinyecm_edw.h"
extern uint64_t tedw_stat_curves, tedw_stat_gkl, tedw_stat_starts;
int main(int argc, char **argv)
{
    int t = argc > 1 ? atoi(argv[1]) : 22, num = argc > 2 ? atoi(argv[2]) : 2000;
    gmp_randstate_t rs; mpz_t p, q, n, f; uint64_t lcg = 1;
    gmp_randinit_default(rs); gmp_randseed_ui(rs, 5); mpz_inits(p, q, n, f, NULL);
    double tt = 0;
    for (int i = 0; i < num; i++) {
        mpz_urandomb(p, rs, t); mpz_setbit(p, t - 1); mpz_nextprime(p, p);
        mpz_urandomb(q, rs, 2 * t); mpz_setbit(q, 2 * t - 1); mpz_nextprime(q, q);
        mpz_mul(n, p, q);
        double t0 = now();
        getfactor_tecm_edw(n, f, t, &lcg);
        tt += now() - t0;
    }
    printf("t %d: curves/input %.2f, gkl curves/input %.2f, starts/input %.3f, us/input %.2f, us/curve %.2f\n", t,
        (double)tedw_stat_curves / num, (double)tedw_stat_gkl / num, (double)tedw_stat_starts / num, tt * 1e6 / num, tt * 1e6 / (double)tedw_stat_curves);
    return 0;
}
