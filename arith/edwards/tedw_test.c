/*
 * Test driver for factor/gmp-ecm/tinyecm_edw.c, used by uedw_oracle.py
 * (t* modes).  Reads commands from stdin, numbers in decimal:
 *   c n id B1 B2          -> stage f        (one scalar curve)
 *   s n id B1             -> X Z            (Montgomery image after stage 1)
 *   x B1 curves k n_1..   -> total f_1 ..   (x8 list, n < 2^104)
 *   o n B1 curves         -> found rounds f (x8 on one n)
 *   g n target_bits       -> found f        (getfactor_tecm_edw)
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "gmp.h"
#include "tinyecm_edw.h"

static char buf[1 << 16];

static int rd(mpz_t x)
{
    if (scanf("%65535s", buf) != 1) return 0;
    return mpz_set_str(x, buf, 10) == 0;
}

int main(void)
{
    char cmd;
    uint64_t lcg = 1;
    mpz_t n, f, X, Z, t;
    mpz_inits(n, f, X, Z, t, NULL);
    static uint64_t q[2 * 100000], fo[2 * 100000];
    while (scanf(" %c", &cmd) == 1) {
        if (cmd == 'c') {
            unsigned id, B1, B2;
            int stage;
            if (!rd(n) || scanf("%u %u %u", &id, &B1, &B2) != 3) break;
            tecm_edw_one_curve(n, id, B1, B2, f, &stage);
            gmp_printf("%d %Zd\n", stage, f);
        } else if (cmd == 's') {
            unsigned id, B1;
            if (!rd(n) || scanf("%u %u", &id, &B1) != 2) break;
            tecm_edw_stage1_point(n, id, B1, X, Z);
            gmp_printf("%Zd %Zd\n", X, Z);
        } else if (cmd == 'x') {
            unsigned B1, curves, k, i;
            if (scanf("%u %u %u", &B1, &curves, &k) != 3 || k > 100000) break;
            for (i = 0; i < k; i++) {
                if (!rd(n)) return 1;
                mpz_tdiv_r_2exp(t, n, 52);
                q[2 * i] = mpz_get_ui(t);
                mpz_tdiv_q_2exp(t, n, 52);
                q[2 * i + 1] = mpz_get_ui(t);
            }
            uint64_t tot = tecm_edw_x8_run(q, fo, k, B1, curves);
            printf("%llu", (unsigned long long)tot);
            for (i = 0; i < k; i++) {
                mpz_set_ui(t, fo[2 * i + 1]);
                mpz_mul_2exp(t, t, 64);
                mpz_add_ui(t, t, fo[2 * i]);
                gmp_printf(" %Zd", t);
            }
            printf("\n");
        } else if (cmd == 'o') {
            unsigned B1, curves;
            int rounds;
            if (!rd(n) || scanf("%u %u", &B1, &curves) != 2) break;
            int found = tecm_edw_x8_one(n, f, B1, curves, &rounds);
            gmp_printf("%d %d %Zd\n", found, rounds, f);
        } else if (cmd == 'g') {
            int tb;
            if (!rd(n) || scanf("%d", &tb) != 1) break;
            int found = getfactor_tecm_edw(n, f, tb, &lcg);
            gmp_printf("%d %Zd\n", found, f);
        }
        fflush(stdout);
    }
    return 0;
}
