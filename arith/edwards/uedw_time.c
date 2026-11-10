/* Per-curve time of microecm_edw's scalar code at each level: stage 1 only
   (B2 = B1) and stage 1 + 2, through uecm_edw_one_curve on table curve 0
   (cheap build), n = a 64-bit product of two 32-bit primes (finds rare).
   Output: B1, ns stage 1, ns stage 1+2, ns per call of a no-op level
   (context setup, for subtraction).  Calibrates b1b2_study.py's cost model.
   gcc -O2 -march=native -I../../factor/gmp-ecm uedw_time.c \
       ../../factor/gmp-ecm/microecm_edw.c */
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "microecm_edw.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main(void)
{
    static const uint32_t lv[7] = { 27, 47, 70, 85, 125, 165, 205 };
    static const uint32_t mult[7] = { 25, 25, 25, 25, 50, 50, 50 };
    /* 4294967291 * 4294967279 and a few others near 2^64 */
    static const uint64_t ns[4] = { 18446743979220271189ULL, 18446744030759878681ULL,
        14757395258967641293ULL, 12297829382473034411ULL };
    const int reps = 20000;
    int k, i, st;
    uint64_t sink = 0;
    for (k = 0; k < 7; k++) {
        double best1 = 1e9, best2 = 1e9;
        int trial;
        for (trial = 0; trial < 3; trial++) {
            double t0 = now();
            for (i = 0; i < reps; i++)
                sink += uecm_edw_one_curve(ns[i & 3], UECM_EDW_TABLE | 0, lv[k], lv[k], &st);
            double t1 = now();
            for (i = 0; i < reps; i++)
                sink += uecm_edw_one_curve(ns[i & 3], UECM_EDW_TABLE | 0, lv[k], mult[k] * lv[k], &st);
            double t2 = now();
            if ((t1 - t0) < best1) best1 = t1 - t0;
            if ((t2 - t1) < best2) best2 = t2 - t1;
        }
        printf("%u %.1f %.1f\n", lv[k], best1 * 1e9 / reps, best2 * 1e9 / reps);
    }
    fprintf(stderr, "%llu\n", (unsigned long long)sink);
    return 0;
}
