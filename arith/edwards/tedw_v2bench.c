/* B1 / B2 study, tinyecm_edw: the legacy parameter set (1) vs the re-tuned
   one (2), inputs as tiny_e2e.c (p of t bits, q of 2t bits, target t):
   scalar getfactor_tecm_edw, the x8 list and the x8 one-n version.
   usage: tedw_v2bench [num [rounds [seed]]]
   gcc -O2 -march=native -DUSE_AVX2 -DUSE_AVX512F [-DIFMA] -I... tedw_v2bench.c \
       ../../factor/gmp-ecm/tinyecm_edw.c $G/lib/libgmp.a */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "gmp.h"
#include "tinyecm_edw.h"

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
    int num = argc > 1 ? atoi(argv[1]) : 2000;
    int rounds = argc > 2 ? atoi(argv[2]) : 3;
    static const int tb[7] = { 20, 22, 24, 26, 29, 31, 32 };
    uint64_t *q = malloc(2 * num * sizeof(uint64_t)), *fo = malloc(2 * num * sizeof(uint64_t));
    mpz_t *nn = malloc(num * sizeof(mpz_t));
    mpz_t p, r, g, f;
    uint64_t lcg = 77;
    int k, i, rep, v, mode;
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, argc > 3 ? atoi(argv[3]) : 1);
    mpz_inits(p, r, g, f, NULL);
    for (i = 0; i < num; i++)
        mpz_init(nn[i]);
    printf("  t | mode    | legacy us  found | v2 us  found | speedup\n");
    for (k = 0; k < 7; k++) {
        int t = tb[k];
        for (i = 0; i < num; i++) {
            rprime(p, t);
            rprime(r, 2 * t);
            mpz_mul(nn[i], p, r);
            mpz_tdiv_r_2exp(g, nn[i], 52);
            q[2 * i] = mpz_get_ui(g);
            mpz_tdiv_q_2exp(g, nn[i], 52);
            q[2 * i + 1] = mpz_get_ui(g);
        }
        for (mode = 0; mode < 3; mode++) {
            static const char *mn[3] = { "scalar", "x8 list", "x8 one" };
            int n1 = (mode == 1) ? num : num / 2;
            double best[3] = { 1e30, 1e30, 1e30 };
            int ok[3] = { 0, 0, 0 };
            for (rep = 0; rep < rounds; rep++) {
                for (v = 1; v <= 2; v++) {
                    int okv = 0;
                    tecm_edw_set_params(v);
                    double t0 = now();
                    if (mode == 0) {
                        for (i = 0; i < n1; i++)
                            if (getfactor_tecm_edw(nn[i], f, t, &lcg) && mpz_cmp_ui(f, 1) > 0 &&
                                mpz_cmp(f, nn[i]) < 0 && mpz_divisible_p(nn[i], f))
                                okv++;
                    } else if (mode == 1) {
                        getfactor_tecm_edw_x8_list(q, fo, t, num, &lcg);
                    } else {
                        for (i = 0; i < n1; i++)
                            if (getfactor_tecm_edw_x8(nn[i], f, t, &lcg) && mpz_cmp_ui(f, 1) > 0 &&
                                mpz_cmp(f, nn[i]) < 0 && mpz_divisible_p(nn[i], f))
                                okv++;
                    }
                    double t1 = now();
                    if (mode == 1) {
                        for (i = 0; i < num; i++) {
                            mpz_set_ui(f, fo[2 * i + 1]);
                            mpz_mul_2exp(f, f, 64);
                            mpz_add_ui(f, f, fo[2 * i]);
                            if (mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, nn[i]) < 0 &&
                                mpz_divisible_p(nn[i], f))
                                okv++;
                        }
                    }
                    if (t1 - t0 < best[v]) best[v] = t1 - t0;
                    ok[v] = okv;
                }
            }
            printf("%3d | %-7s | %8.2f %.4f | %8.2f %.4f | %.3f\n", t, mn[mode],
                best[1] * 1e6 / n1, (double)ok[1] / n1, best[2] * 1e6 / n1, (double)ok[2] / n1,
                best[1] / best[2]);
            fflush(stdout);
        }
    }
    return 0;
}
