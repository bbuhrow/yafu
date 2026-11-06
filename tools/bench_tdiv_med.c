/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// Standalone benchmark of the medium-prime divisibility tests of
// factor/qs/tdiv_med_32k.c, for x86 builds:
//   new scalar  - the generic loop of tdiv_medprimes_32k_new (32-bit Lemire test)
//   new vector  - its ISA blocks: 16 lanes (AVX512F), 8 (AVX2) or 4 (SSE4.1 /
//                 SSE2, helpers from factor/qs/tdiv_med_sse.h)
//   old scalar  - the 16-bit Barrett test of tdiv_medprimes_32k_old, one
//                 prime at a time
//   old vector  - tdiv_medprimes_32k_old's range cascade, 32 (AVX512BW) /
//                 16 (AVX2) / 8 (SSE2) lanes, from factor/qs/tdiv_med_old.h
// The new loops are copied from tdiv_medprimes_32k_new; the old ones and
// the SSE helpers are the headers yafu compiles.  A hit is recorded instead
// of trial dividing Q (DIVIDE_RESIEVED_PRIME_2), so only the test is timed,
// and every kernel's hit list is checked against the new scalar one.
//
// Build from the yafu root, one binary per ISA level (gcc flags as in the
// Makefile's gcc section; the -D chain as the Makefile sets it):
//   F="-O2 -fno-tree-loop-vectorize -fno-tree-slp-vectorize -Ifactor/qs"
//   gcc $F -mno-sse4.1 -o bench_sse2 tools/bench_tdiv_med.c
//   gcc $F -msse4.1 -DUSE_SSE41 -o bench_sse41 tools/bench_tdiv_med.c
//   gcc $F -mavx2 -DUSE_SSE41 -DUSE_AVX2 -o bench_avx2 tools/bench_tdiv_med.c
//   gcc $F -march=skylake-avx512 -DUSE_SSE41 -DUSE_AVX2 -DUSE_AVX512F -DUSE_AVX512BW -o bench_avx512 tools/bench_tdiv_med.c
// Add -DTDIV_OLD_LOC_BASE=32768 to time the original old test (it then
// misses the block-end hits, which the check reports).
//
// Usage: bench [pmin [reps]]
//   tests the fb primes from just above pmin (default 256, about
//   sieve_small_fb_start) to fb_13bit_B (the first index that is a multiple
//   of 8 with a prime > 8192); about half of all primes are kept, like a
//   real factor base.  4096 random block locations, reps passes (default 20).

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <x86intrin.h>

// what yafu's headers provide to the tdiv_med headers
#define SIQS_FORCEINLINE static inline __attribute__((always_inline))
#define _trail_zcnt __builtin_ctz
#define SIQS_HAVE_SINV16
#if defined(__SSE2__)
#define D_HAS_SSE2
#endif

#define MAXFB 2048
#define NLOC 4096

typedef struct
{
    uint16_t* prime;
    uint16_t* root1;
    uint16_t* root2;
} bench_fb_t;
#define TDM_FB_T bench_fb_t

#include "tdiv_med_sse.h"
#include "tdiv_med_old.h"

// the fb arrays as yafu lays them out (64-byte aligned, as xmalloc_align)
static _Alignas(64) uint16_t fb_prime[MAXFB + 64], fb_root1[MAXFB + 64], fb_root2[MAXFB + 64];
static _Alignas(64) uint32_t fb_minv[MAXFB + 64], fb_thresh[MAXFB + 64];
static _Alignas(64) uint16_t fb_sinv[MAXFB + 64], fb_corr[MAXFB + 64];
static uint32_t locs[NLOC];
static uint32_t g_start, g_bound, g_b[2];
static bench_fb_t g_fb = { fb_prime, fb_root1, fb_root2 };

// one prime of the new (Lemire) test
#define NEW_TEST(j) \
    do { \
        uint32_t _p = fb_prime[j]; \
        uint32_t _off = base + (_p << 15) + _p; \
        if ((uint32_t)((_off - fb_root1[j]) * fb_minv[j]) <= fb_thresh[j] || \
            (uint32_t)((_off - fb_root2[j]) * fb_minv[j]) <= fb_thresh[j]) \
            hits[nh++] = (uint16_t)(j); \
    } while (0)

// new scalar: the generic #else loop of tdiv_medprimes_32k_new
__attribute__((noinline))
static uint32_t new_scalar(uint32_t block_loc, uint16_t* hits)
{
    uint32_t i, nh = 0;
    uint32_t base = block_loc - 32768;

    for (i = g_start; i < g_bound; i++)
        NEW_TEST(i);
    return nh;
}

// new vector: the ISA blocks of tdiv_medprimes_32k_new
__attribute__((noinline))
static uint32_t new_vector(uint32_t block_loc, uint16_t* hits)
{
    uint32_t i = g_start, nh = 0;
    uint32_t base = block_loc - 32768;
    uint32_t bound = g_bound;

#if defined(USE_AVX512F)
    for (; i < bound; i++)
    {
        if ((i & 15) == 0) break;
        NEW_TEST(i);
    }

    __m512i vbase = _mm512_set1_epi32((int)(block_loc - 32768));

    for (; i + 16 < bound; i += 16)
    {
        __m512i vp = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) & fb_prime[i]));
        __m512i vr1 = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) & fb_root1[i]));
        __m512i vr2 = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) & fb_root2[i]));
        __m512i vminv = _mm512_load_si512(&fb_minv[i]);
        __m512i vthresh = _mm512_load_si512(&fb_thresh[i]);

        __m512i voff = _mm512_add_epi32(vbase,
            _mm512_add_epi32(_mm512_slli_epi32(vp, 15), vp));

        __m512i t1 = _mm512_sub_epi32(voff, vr1);
        __m512i t2 = _mm512_sub_epi32(voff, vr2);

        __mmask16 hit1 = _mm512_cmple_epu32_mask(_mm512_mullo_epi32(t1, vminv), vthresh);
        __mmask16 hit2 = _mm512_cmple_epu32_mask(_mm512_mullo_epi32(t2, vminv), vthresh);
        unsigned mask = (unsigned)(hit1 | hit2);

        while (mask)
        {
            int lane = _trail_zcnt(mask);
            mask &= mask - 1;
            hits[nh++] = (uint16_t)(i + lane);
        }
    }

#elif defined(USE_AVX2)
    __m256i vbase = _mm256_set1_epi32((int)(block_loc - 32768));

    for (; i < bound; i++)
    {
        if ((i & 7) == 0) break;
        NEW_TEST(i);
    }

    for (; i + 8 < bound; i += 8)
    {
        __m256i vp = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) & fb_prime[i]));
        __m256i vr1 = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) & fb_root1[i]));
        __m256i vr2 = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) & fb_root2[i]));
        __m256i vminv = _mm256_load_si256((const __m256i*) & fb_minv[i]);
        __m256i vthresh = _mm256_load_si256((const __m256i*) & fb_thresh[i]);

        __m256i voff = _mm256_add_epi32(vbase,
            _mm256_add_epi32(_mm256_slli_epi32(vp, 15), vp));

        __m256i x1 = _mm256_mullo_epi32(_mm256_sub_epi32(voff, vr1), vminv);
        __m256i x2 = _mm256_mullo_epi32(_mm256_sub_epi32(voff, vr2), vminv);

        __m256i h1 = _mm256_cmpeq_epi32(_mm256_min_epu32(x1, vthresh), x1);
        __m256i h2 = _mm256_cmpeq_epi32(_mm256_min_epu32(x2, vthresh), x2);
        unsigned mask = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_or_si256(h1, h2)));

        while (mask) {
            int lane = _trail_zcnt(mask);
            mask &= mask - 1;
            hits[nh++] = (uint16_t)(i + lane);
        }
    }

#elif defined(USE_SSE41) || defined(D_HAS_SSE2)
    for (; i < bound; i++)
    {
        if ((i & 7) == 0) break;
        NEW_TEST(i);
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
            int lane = _trail_zcnt(mask);
            mask &= mask - 1;
            hits[nh++] = (uint16_t)(i + lane);
        }
    }
#endif

    for (; i < bound; i++)
        NEW_TEST(i);
    return nh;
}

// old scalar: the 16-bit test of tdiv_medprimes_32k_old, one prime at a time
__attribute__((noinline))
static uint32_t old_scalar(uint32_t block_loc, uint16_t* hits)
{
    const TDM_FB_T* fbc = &g_fb;
    const uint16_t* sinv = fb_sinv;
    const uint16_t* corr = fb_corr;
    uint16_t* buf = hits;
    uint32_t i, nb = 0;

    for (i = g_start; i < g_b[0]; i++) TDM_OLD_SCALAR(i, 8);
    for (; i < g_b[1]; i++) TDM_OLD_SCALAR(i, 10);
    for (; i < g_bound; i++) TDM_OLD_SCALAR(i, 12);
    return nb;
}

// old vector: tdiv_medprimes_32k_old's three range calls
__attribute__((noinline))
static uint32_t old_vector(uint32_t block_loc, uint16_t* hits)
{
    uint32_t nb = 0;

    nb = tdm_old_range(&g_fb, fb_sinv, fb_corr, g_start, g_b[0], block_loc, 8, hits, nb);
    nb = tdm_old_range(&g_fb, fb_sinv, fb_corr, g_b[0], g_b[1], block_loc, 10, hits, nb);
    nb = tdm_old_range(&g_fb, fb_sinv, fb_corr, g_b[1], g_bound, block_loc, 12, hits, nb);
    return nb;
}

typedef uint32_t (*test_fn)(uint32_t, uint16_t*);

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}

// best of 7 timed passes; returns ns per report, *cyc = TSC ticks per report
static double bench(test_fn fn, int reps, uint64_t* sum, double* cyc)
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
                uint32_t nh = fn(locs[k], hits);
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
    int reps = (argc > 2) ? atoi(argv[2]) : 20;
    uint32_t n = 2, p, k, j, shift = 24;
    uint64_t nhits = 0;
    static const struct { const char* name; test_fn fn; } kern[] = {
        { "new scalar", new_scalar }, { "new vector", new_vector },
        { "old scalar", old_scalar }, { "old vector", old_vector } };
    const int nk = sizeof(kern) / sizeof(kern[0]);
    double tns[4], tcyc[4];
    uint64_t miss[4] = { 0, 0, 0, 0 }, extra[4] = { 0, 0, 0, 0 };

    srand(1);

    // fb as make_fb_siqs builds it: index 0, 1 unused, then about half of
    // the odd primes (those where N is a QR), up to fb_13bit_B
    fb_prime[0] = 1;
    fb_prime[1] = 2;
    g_start = 0;
    for (p = 3; n < MAXFB; p += 2)
    {
        if (!isprime(p) || ((p > 7) && (rand() & 1)))
            continue;
        if ((g_start == 0) && (p > pmin))
            g_start = n;
        // fb_13bit_B: the first multiple of 8 whose prime exceeds 8192
        if ((p > 8192) && ((n & 7) == 0))
            break;
        fb_prime[n++] = (uint16_t)p;
    }
    g_bound = n;

    for (j = 2; j < g_bound; j++)
    {
        uint32_t q = fb_prime[j], inv = q;
        inv *= 2 - q * inv; inv *= 2 - q * inv; inv *= 2 - q * inv; inv *= 2 - q * inv;
        fb_minv[j] = inv;
        fb_thresh[j] = 0xFFFFFFFFu / q;
        // roots as the sieve leaves them: in [0, p)
        fb_root1[j] = (uint16_t)(rand() % q);
        fb_root2[j] = (uint16_t)(rand() % q);

        // the USE_8X_MOD_ASM 16-bit tables (siqs_aux.c), for p >= 256
        if (q >= 256)
        {
            uint16_t sinv;
            if ((shift == 24) && (q > 1024) && (j % 8 == 0)) shift = 26;
            if ((shift == 26) && (q > 4096) && (j % 8 == 0)) shift = 28;
            sinv = (uint16_t)(((uint32_t)1 << shift) / q);
            if (floor((double)(1 << shift) / (double)q + 0.5) == (double)sinv)
                fb_corr[j] = 1;
            else
            {
                fb_corr[j] = 0;
                sinv++;
            }
            fb_sinv[j] = sinv;
        }
    }
    tdm_old_bounds(fb_prime, g_start, g_bound, g_b);

    for (k = 0; k < NLOC; k++)
        locs[k] = (uint32_t)rand() & 32767;
    // a few block locations where some progression ends exactly at the
    // block end (stored root 0), the case the original old test missed
    for (k = 0; k < 64; k++)
    {
        j = g_start + (uint32_t)rand() % (g_bound - g_start);
        fb_root1[j] = 0;
        locs[k * 61] = 32768 % fb_prime[j];
    }

    // correctness: every kernel against new scalar, every location
    for (k = 0; k < NLOC; k++)
    {
        uint16_t ref[MAXFB], h[MAXFB];
        uint32_t nr = new_scalar(locs[k], ref);
        int f;

        nhits += nr;
        for (f = 1; f < nk; f++)
        {
            uint32_t nh = kern[f].fn(locs[k], h), a = 0, b = 0;
            // both lists are ascending: count misses and extras
            while ((a < nr) || (b < nh))
            {
                if ((b >= nh) || ((a < nr) && (ref[a] < h[b]))) { miss[f]++; a++; }
                else if ((a >= nr) || (h[b] < ref[a])) { extra[f]++; b++; }
                else { a++; b++; }
            }
        }
    }

    printf("ISA: %s;  old vector lanes: %s\n",
#if defined(USE_AVX512F)
        "AVX-512 (new vector: 16 x 32-bit)",
#elif defined(USE_AVX2)
        "AVX2 (new vector: 8 x 32-bit)",
#elif defined(USE_SSE41)
        "SSE4.1 (new vector: 4 x 32-bit)",
#else
        "SSE2 (new vector: 4 x 32-bit, emulated pmulld/pminud)",
#endif
#if defined(USE_AVX512BW)
        "32/16/8 x 16-bit"
#elif defined(USE_AVX2)
        "16/8 x 16-bit"
#else
        "8 x 16-bit"
#endif
    );
    printf("fb primes [%u, %u) = %u..%u (%u primes), shift 26 from index %u, 28 from %u; "
        "old LOC_BASE %d\n", g_start, g_bound, fb_prime[g_start], fb_prime[g_bound - 1],
        g_bound - g_start, g_b[0], g_b[1], TDIV_OLD_LOC_BASE);
    printf("%d x %d reports, %.2f hits per report\n", reps, NLOC, (double)nhits / NLOC);

    for (k = 0; k < (uint32_t)nk; k++)
    {
        uint64_t s;
        tns[k] = bench(kern[k].fn, reps, &s, &tcyc[k]);
    }

    printf("%-11s %10s %10s %9s %10s  %s\n", "kernel", "ns/report", "TSC/report",
        "TSC/prime", "vs new sc.", "check");
    for (k = 0; k < (uint32_t)nk; k++)
    {
        char chk[64];
        if (k == 0) strcpy(chk, "reference");
        else if (miss[k] || extra[k])
            sprintf(chk, "%llu missed, %llu extra", (unsigned long long)miss[k],
                (unsigned long long)extra[k]);
        else strcpy(chk, "identical");
        printf("%-11s %10.1f %10.0f %9.2f %9.2fx  %s\n", kern[k].name, tns[k], tcyc[k],
            tcyc[k] / (g_bound - g_start), tns[0] / tns[k], chk);
    }
    printf("old vector vs new vector: %.2fx\n", tns[1] / tns[3]);

    for (k = 1; k < (uint32_t)nk; k++)
        if (miss[k] || extra[k])
            return (TDIV_OLD_LOC_BASE == 32768) ? 0 : 1;
    return 0;
}
