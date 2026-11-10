/* Failures of getfactor_uecm_edw on the semiprimes_*bit.dat files, replaying
   test_dlp_composites() (top/test.c, microecm_marker) exactly: files 32..64
   in order, 100000 inputs each, one LCG (seed 42) across all files.
   usage: uedw_fails [params [extra]]   (run from the yafu root)
     params: 1 = legacy, 2 = re-tuned (default)
     extra:  if > 0, retry each failure with `extra` more curves at the same
             level (uecm_edw_level, fresh LCG) and report what finds it.
   Prints per file the fail count and each failing n with p, q, the largest
   prime factor of p-1, q-1 (P-1 reach) and the bits of p, q.
   gcc -O2 -march=native -I../../factor/gmp-ecm -I../../include -I../../ytools \
       -I$G/include uedw_fails.c ../../factor/gmp-ecm/microecm_edw.c \
       ../../factor/gmp-ecm/micropm1.c $G/lib/libgmp.a -o uedw_fails
   (or -DUEDW_NO_PM1 without micropm1.c / GMP: no P-1 pre-pass) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "microecm_edw.h"

static uint64_t lpf(uint64_t x)
{
    uint64_t p, big = 1;
    while ((x & 1) == 0) { x >>= 1; big = 2; }
    for (p = 3; p * p <= x; p += 2)
        while (x % p == 0) { x /= p; big = p; }
    return x > 1 ? (x > big ? x : big) : big;
}

static int nbits(uint64_t x) { int b = 0; while (x) { b++; x >>= 1; } return b; }

int main(int argc, char **argv)
{
    int params = argc > 1 ? atoi(argv[1]) : 2;
    int extra = argc > 2 ? atoi(argv[2]) : 0;
    static uint64_t comp[200000];
    static uint32_t f1[200000], f2[200000];
    uint64_t lcg = 42;
    int bits, tot = 0;
    uecm_edw_set_params(params);
    for (bits = 32; bits <= 64; bits += 2) {
        char fn[64];
        FILE *in;
        int i, nf = 0;
        sprintf(fn, "semiprimes_%dbit.dat", bits);
        in = fopen(fn, "r");
        if (!in) { printf("no %s\n", fn); return 1; }
        i = 0;
        while (!feof(in) && i < 200000) {
            if (fscanf(in, "%lu,%u,%u", &comp[i], &f1[i], &f2[i]) != 3) break;
            i++;
        }
        fclose(in);
        for (i = 0; i < 100000; i++) {
            uint64_t r = getfactor_uecm_edw(comp[i], 0, &lcg);
            if (r == f1[i] || r == f2[i])
                continue;
            nf++;
            if (bits >= 54) {
                uint64_t p = f1[i], q = f2[i];
                printf("  %2d %20lu = %10lu * %10lu  (%d,%d bits)  lpf(p-1)=%lu lpf(q-1)=%lu  ret=%lu",
                    bits, comp[i], p, q, nbits(p), nbits(q), lpf(p - 1), lpf(q - 1), r);
                if (extra > 0) {
                    uint64_t l2 = 12345 + i;
                    int nc;
                    uint32_t B1 = bits <= 56 ? 125 : bits <= 60 ? 165 : 205;
                    uint64_t g = uecm_edw_level(comp[i], B1, extra, &l2, &nc);
                    printf("  retry B1=%u: %s after %d", B1, g > 1 ? "found" : "none", nc);
                }
                printf("\n");
            }
        }
        tot += nf;
        printf("%d bits: %d fails\n", bits, nf);
    }
    printf("total %d\n", tot);
    return 0;
}
