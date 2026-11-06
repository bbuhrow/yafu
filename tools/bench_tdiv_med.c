/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// Standalone benchmark of the medium-prime divisibility test in
// tdiv_medprimes_32k_new (factor/qs/tdiv_med_32k.c): the generic
// one-prime-at-a-time C loop against the 128-bit vector block, whose
// helpers (factor/qs/tdiv_med_sse.h) are SSE2 or SSE4.1 depending on
// USE_SSE41.  Both loops are copied from tdiv_medprimes_32k_new, except
// that a hit is recorded instead of trial dividing Q (DIVIDE_RESIEVED_PRIME_2),
// so only the test itself is timed.  The hit lists are checked for equality.
//
// Build (gcc flags as in the Makefile's gcc section), from the yafu root:
//   gcc -O2 -fno-tree-loop-vectorize -fno-tree-slp-vectorize -mno-sse4.1 -Ifactor/qs -o bench_sse2 tools/bench_tdiv_med.c
//   gcc -O2 -fno-tree-loop-vectorize -fno-tree-slp-vectorize -msse4.1 -DUSE_SSE41 -Ifactor/qs -o bench_sse41 tools/bench_tdiv_med.c
//
// Usage: bench [pmin [pmax [reps]]]
//   test primes in (pmin, pmax) (default 256..8192, i.e. from about
//   sieve_small_fb_start to fb_13bit_B); about half of the primes are kept,
//   like a real factor base.  4096 random block locations, reps passes.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>
#include "tdiv_med_sse.h"

#define NLOC 4096
#define MAXFB 2048

// the fb arrays, as yafu lays them out: uint16 primes/roots, uint32 consts
static uint16_t fb_prime[MAXFB + 8], fb_root1[MAXFB + 8], fb_root2[MAXFB + 8];
static uint32_t fb_minv[MAXFB + 8], fb_thresh[MAXFB + 8];
static uint32_t locs[NLOC];

// scalar: the generic #else loop of tdiv_medprimes_32k_new
__attribute__((noinline))
static uint32_t test_scalar(uint32_t block_loc, uint32_t start, uint32_t bound, uint16_t* hits)
{
    uint32_t i, nh = 0;
    uint32_t base = block_loc - 32768;

    for (i = start; i < bound; i++)
    {
        uint32_t prime = fb_prime[i];
        uint32_t root1 = fb_root1[i];
        uint32_t root2 = fb_root2[i];
        uint32_t minv = fb_minv[i];
        uint32_t thresh = fb_thresh[i];
        uint32_t off = base + ((uint32_t)prime << 15) + prime;

        if ((uint32_t)((off - root1) * minv) <= thresh ||
            (uint32_t)((off - root2) * minv) <= thresh)
        {
            hits[nh++] = (uint16_t)i;
        }
    }
    return nh;
}

// vector: the USE_SSE41 || D_HAS_SSE2 block of tdiv_medprimes_32k_new
__attribute__((noinline))
static uint32_t test_vec(uint32_t block_loc, uint32_t start, uint32_t bound, uint16_t* hits)
{
    uint32_t i = start, nh = 0;
    uint32_t base = block_loc - 32768;

    for (; i < bound; i++)
    {
        if ((i & 7) == 0) break;

        uint32_t prime = fb_prime[i];
        uint32_t root1 = fb_root1[i];
        uint32_t root2 = fb_root2[i];
        uint32_t minv = fb_minv[i];
        uint32_t thresh = fb_thresh[i];
        uint32_t off = base + ((uint32_t)prime << 15) + prime;

        if ((uint32_t)((off - root1) * minv) <= thresh ||
            (uint32_t)((off - root2) * minv) <= thresh)
        {
            hits[nh++] = (uint16_t)i;
        }
    }

    __m128i vbase = _mm_set1_epi32((int)(block_loc - 32768));

    for (; i + 4 < bound; i += 4)
    {
        __m128i vp = tdm_load4_u16(&fb_prime[i]);
        __m128i vr1 = tdm_load4_u16(&fb_root1[i]);
        __m128i vr2 = tdm_load4_u16(&fb_root2[i]);
        __m128i vminv = _mm_loadu_si128((const __m128i*) & fb_minv[i]);
        __m128i vthresh = _mm_loadu_si128((const __m128i*) & fb_thresh[i]);

        __m128i voff = _mm_add_epi32(vbase,
            _mm_add_epi32(_mm_slli_epi32(vp, 15), vp));

        __m128i x1 = tdm_mullo_epi32(_mm_sub_epi32(voff, vr1), vminv);
        __m128i x2 = tdm_mullo_epi32(_mm_sub_epi32(voff, vr2), vminv);

        __m128i h1 = tdm_cmple_epu32(x1, vthresh);
        __m128i h2 = tdm_cmple_epu32(x2, vthresh);
        unsigned mask = _mm_movemask_ps(_mm_castsi128_ps(_mm_or_si128(h1, h2)));

        while (mask) {
            int lane = __builtin_ctz(mask);
            mask &= mask - 1;
            hits[nh++] = (uint16_t)(i + lane);
        }
    }

    for (; i < bound; i++)
    {
        uint32_t prime = fb_prime[i];
        uint32_t root1 = fb_root1[i];
        uint32_t root2 = fb_root2[i];
        uint32_t minv = fb_minv[i];
        uint32_t thresh = fb_thresh[i];
        uint32_t off = base + ((uint32_t)prime << 15) + prime;

        if ((uint32_t)((off - root1) * minv) <= thresh ||
            (uint32_t)((off - root2) * minv) <= thresh)
        {
            hits[nh++] = (uint16_t)i;
        }
    }
    return nh;
}

typedef uint32_t (*test_fn)(uint32_t, uint32_t, uint32_t, uint16_t*);

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}

// best of 7 timed passes; returns ns per report, *cyc = TSC ticks per report
static double bench(test_fn fn, uint32_t start, uint32_t bound, int reps,
    uint64_t* sum, double* cyc)
{
    uint16_t hits[MAXFB];
    double best = 1e30, bestc = 1e30;
    int t, r, k;

    for (t = 0; t < 7; t++)
    {
        uint64_t s = 0;
        double t0 = now_ns();
        uint64_t c0 = __rdtsc();
        for (r = 0; r < reps; r++)
            for (k = 0; k < NLOC; k++)
            {
                uint32_t nh = fn(locs[k], start, bound, hits);
                s += nh + (nh ? hits[nh - 1] : 0);
            }
        uint64_t c1 = __rdtsc();
        double t1 = now_ns();
        double per = (t1 - t0) / ((double)reps * NLOC);
        double perc = (double)(c1 - c0) / ((double)reps * NLOC);
        if (per < best) best = per;
        if (perc < bestc) bestc = perc;
        *sum = s;
    }
    *cyc = bestc;
    return best;
}

static int isprime(uint32_t n)
{
    uint32_t d;
    if (n < 2) return 0;
    for (d = 2; d * d <= n; d++)
        if (n % d == 0) return 0;
    return 1;
}

int main(int argc, char** argv)
{
    uint32_t pmin = (argc > 1) ? (uint32_t)atoi(argv[1]) : 256;
    uint32_t pmax = (argc > 2) ? (uint32_t)atoi(argv[2]) : 8192;
    int reps = (argc > 3) ? atoi(argv[3]) : 20;
    uint32_t n = 0, p, k, start, bound;
    uint64_t s1, s2, nhits = 0;
    uint16_t h1[MAXFB], h2[MAXFB];
    double c1, c2;

    srand(1);

    // a few smaller primes first, so start is not a multiple of 8 and the
    // vector code's head loop runs as in yafu
    for (p = 3; p <= pmin && n < 5; p += 2)
        if (isprime(p)) { fb_prime[n++] = (uint16_t)p; }
    start = n;
    for (p = pmin + 1; (p < pmax) && (n < MAXFB); p++)
    {
        // keep about half: the primes where N is a QR
        if (!isprime(p) || (rand() & 1)) continue;
        fb_prime[n++] = (uint16_t)p;
    }
    bound = n;

    for (k = 0; k < n; k++)
    {
        uint32_t q = fb_prime[k], inv = q;
        inv *= 2 - q * inv; inv *= 2 - q * inv; inv *= 2 - q * inv; inv *= 2 - q * inv;
        fb_minv[k] = inv;
        fb_thresh[k] = 0xFFFFFFFFu / q;
        // roots as the sieve leaves them: in [0, p)
        fb_root1[k] = (uint16_t)(rand() % q);
        fb_root2[k] = (uint16_t)(rand() % q);
    }
    for (k = 0; k < NLOC; k++)
        locs[k] = (uint32_t)rand() & 32767;

    // correctness: identical hit lists for every location
    for (k = 0; k < NLOC; k++)
    {
        uint32_t a = test_scalar(locs[k], start, bound, h1);
        uint32_t b = test_vec(locs[k], start, bound, h2);
        if ((a != b) || memcmp(h1, h2, a * sizeof(uint16_t)))
        {
            printf("MISMATCH at loc %u: scalar %u hits, vector %u hits\n", locs[k], a, b);
            return 1;
        }
        nhits += a;
    }

    printf("vector helpers: %s\n",
#ifdef USE_SSE41
        "SSE4.1 (pmovzxwd / pmulld / pminud)"
#else
        "SSE2 (punpcklwd / 2x pmuludq + shuffles / biased pcmpgtd)"
#endif
    );
    printf("%u primes tested per report (%u..%u), %.2f hits per report, %d x %d reports\n",
        bound - start, fb_prime[start], fb_prime[bound - 1],
        (double)nhits / NLOC, reps, NLOC);

    double t1 = bench(test_scalar, start, bound, reps, &s1, &c1);
    double t2 = bench(test_vec, start, bound, reps, &s2, &c2);
    if (s1 != s2)
    {
        printf("checksum mismatch\n");
        return 1;
    }

    printf("scalar : %8.1f ns/report  %8.0f TSC/report  %6.2f TSC/prime\n",
        t1, c1, c1 / (bound - start));
    printf("vector : %8.1f ns/report  %8.0f TSC/report  %6.2f TSC/prime\n",
        t2, c2, c2 / (bound - start));
    printf("speedup: %.2fx\n", t1 / t2);
    return 0;
}
