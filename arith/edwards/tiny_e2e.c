/*
 * End-to-end benchmark of microecm / tinyecm through their public entry
 * points, at every supported B1 level: random semiprimes, time per input
 * and the fraction factored.
 *
 *   uecm      getfactor_uecm         n = p q balanced, 40 .. 64 bits
 *   uecm8     getfactor_uecm_x8_list n = p q balanced, 40 .. 52 bits
 *   tecm      getfactor_tecm         p of t bits, q of 2t bits, target t
 *   tecm8     getfactor_tecm_x8_list same inputs, 2 x 52-bit limbs
 *
 * usage: tiny_e2e [num [seed]]
 * Build: link with microecm.c, tinyecm.c, arith/monty.c and GMP.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "gmp.h"
#include "microecm.h"
#include "tinyecm.h"
#ifdef WITH_EDW
#include "microecm_edw.h"
#endif

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
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, argc > 2 ? atoi(argv[2]) : 1);
    uint64_t *q = malloc(2 * num * sizeof(uint64_t));
    uint64_t *f = malloc(2 * num * sizeof(uint64_t));
    mpz_t *nn = malloc(num * sizeof(mpz_t));
    mpz_t p, r, g, gf;
    mpz_inits(p, r, g, gf, NULL);
    for (int i = 0; i < num; i++) mpz_init(nn[i]);
    uint64_t lcg = 0x123456789ULL;

    static const int ub[7] = { 40, 44, 48, 52, 58, 62, 64 };
    static const int ublv[7] = { 27, 47, 70, 85, 125, 165, 205 };
    printf("uecm (scalar 64), balanced semiprimes\n  bits  B1   us/input  found\n");
    for (int k = 0; k < 7; k++) {
        for (int i = 0; i < num; i++) {
            do {
                rprime(p, ub[k] / 2);
                rprime(r, ub[k] - ub[k] / 2);
                mpz_mul(nn[i], p, r);
            } while ((int)mpz_sizeinbase(nn[i], 2) != ub[k]);
            q[i] = mpz_get_ui(nn[i]);
        }
        int ok = 0;
        uint64_t h = 0;
        double t0 = now();
        for (int i = 0; i < num; i++) {
            uint64_t ff = getfactor_uecm(q[i], 0, &lcg);
            ok += (ff > 1 && ff < q[i] && q[i] % ff == 0);
            h = h * 1000003 + ff;
        }
        double t1 = now();
        printf("  %4d %4d  %8.2f  %.4f   hash %016llx\n", ub[k], ublv[k], (t1 - t0) * 1e6 / num,
            (double)ok / num, (unsigned long long)h);
#ifdef WITH_EDW
        {
            /* same inputs, Edwards microecm */
            int okE = 0;
            double e0 = now();
            for (int i = 0; i < num; i++) {
                uint64_t ff = getfactor_uecm_edw(q[i], 0, &lcg);
                okE += (ff > 1 && ff < q[i] && q[i] % ff == 0);
            }
            double e1 = now();
            printf("  %4d %4d  %8.2f  %.4f   (uecm_edw)   speedup %.3f\n", ub[k], ublv[k],
                (e1 - e0) * 1e6 / num, (double)okE / num, (t1 - t0) / (e1 - e0));
        }
#endif
        if (ub[k] <= 52) {
            /* same inputs through the x8 list (n < 2^52) */
            t0 = now();
            getfactor_uecm_x8_list(q, f, num, &lcg);
            t1 = now();
            ok = 0;
            h = 0;
            for (int i = 0; i < num; i++) {
                ok += (f[i] > 1 && f[i] < q[i] && q[i] % f[i] == 0);
                h = h * 1000003 + f[i];
            }
            printf("  %4d %4d  %8.2f  %.4f   hash %016llx (uecm x8 list)\n", ub[k], ublv[k],
                (t1 - t0) * 1e6 / num, (double)ok / num, (unsigned long long)h);
        }
    }

    static const int tb[7] = { 20, 22, 24, 26, 29, 31, 32 };
    printf("tecm: p of t bits, q of 2t bits\n  t  n-bits   us/input  found   | x8 list us/input  found\n");
    for (int k = 0; k < 7; k++) {
        int t = tb[k];
        for (int i = 0; i < num; i++) {
            rprime(p, t);
            rprime(r, 2 * t);
            mpz_mul(nn[i], p, r);
            mpz_tdiv_r_2exp(g, nn[i], 52);
            q[2 * i] = mpz_get_ui(g);
            mpz_tdiv_q_2exp(g, nn[i], 52);
            q[2 * i + 1] = mpz_get_ui(g);
        }
        int ok = 0;
        uint64_t hs = 0, h8 = 0;
        double t0 = now();
        for (int i = 0; i < num; i++) {
            if (getfactor_tecm(nn[i], gf, t, &lcg) && mpz_cmp_ui(gf, 1) > 0 &&
                mpz_cmp(gf, nn[i]) < 0 && mpz_divisible_p(nn[i], gf))
                ok++;
            hs = hs * 1000003 + mpz_get_ui(gf);
        }
        double t1 = now();
        double ts = (t1 - t0) * 1e6 / num;
        double fs = (double)ok / num;
        t0 = now();
        getfactor_tecm_x8_list(q, f, t, num, &lcg);
        t1 = now();
        ok = 0;
        for (int i = 0; i < num; i++) {
            mpz_set_ui(gf, f[2 * i + 1]);
            mpz_mul_2exp(gf, gf, 64);
            mpz_add_ui(gf, gf, f[2 * i]);
            if (mpz_cmp_ui(gf, 1) > 0 && mpz_cmp(gf, nn[i]) < 0 && mpz_divisible_p(nn[i], gf))
                ok++;
            h8 = h8 * 1000003 + f[2 * i];
        }
        printf("  %2d  %4d     %8.2f  %.4f  | %8.2f  %.4f   hash %016llx %016llx\n", t, 3 * t,
            ts, fs, (t1 - t0) * 1e6 / num, (double)ok / num,
            (unsigned long long)hs, (unsigned long long)h8);
    }
    return 0;
}
