/*
 * Test driver for factor/gmp-ecm/microecm_edw.c, used by uedw_oracle.py.
 * Reads commands from stdin, one per line:
 *   p n m          -> X Y Z T dn dd    (curve m mod n, out of Montgomery form)
 *   c n m B1 B2    -> stage factor     (one curve; stage 0 / 1 / 2)
 *   s n m B1       -> X Z              (Montgomery image after stage 1)
 *   f n            -> factor           (getfactor_uecm_edw, not arbitrary)
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "microecm_edw.h"

int main(void)
{
    char cmd;
    uint64_t lcg = 1;
    while (scanf(" %c", &cmd) == 1) {
        if (cmd == 'p') {
            uint64_t n, out[6];
            uint32_t m;
            if (scanf("%" SCNu64 " %u", &n, &m) != 2) break;
            uecm_edw_curve_params(n, m, out);
            printf("%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
                out[0], out[1], out[2], out[3], out[4], out[5]);
        } else if (cmd == 'c') {
            uint64_t n;
            uint32_t m, B1, B2;
            int stage;
            if (scanf("%" SCNu64 " %u %u %u", &n, &m, &B1, &B2) != 4) break;
            uint64_t f = uecm_edw_one_curve(n, m, B1, B2, &stage);
            printf("%d %" PRIu64 "\n", stage, f);
        } else if (cmd == 's') {
            uint64_t n, out[2];
            uint32_t m, B1;
            if (scanf("%" SCNu64 " %u %u", &n, &m, &B1) != 3) break;
            uecm_edw_stage1_point(n, m, B1, out);
            printf("%" PRIu64 " %" PRIu64 "\n", out[0], out[1]);
        } else if (cmd == 'x') {
            /* x B1 curves k n_1 .. n_k -> curves_run f_1 .. f_k (x8 list) */
            static uint64_t qs[100000], fs[100000];
            uint32_t B1, curves, k, i;
            if (scanf("%u %u %u", &B1, &curves, &k) != 3 || k > 100000) break;
            for (i = 0; i < k; i++)
                if (scanf("%" SCNu64, &qs[i]) != 1) return 1;
            uint64_t tot = uecm_edw_x8_run(qs, fs, k, B1, curves);
            printf("%" PRIu64, tot);
            for (i = 0; i < k; i++)
                printf(" %" PRIu64, fs[i]);
            printf("\n");
        } else if (cmd == 'f') {
            uint64_t n;
            if (scanf("%" SCNu64, &n) != 1) break;
            printf("%" PRIu64 "\n", getfactor_uecm_edw(n, 0, &lcg));
        }
        fflush(stdout);
    }
    return 0;
}
