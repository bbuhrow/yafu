/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// Shared parts of the 16-bit Barrett ("old") tdiv_med kernels in
// tdiv_med_32k.c: the region bounds, the one-prime test and the x86
// 32/16/8-lane range cascade.  Also used by tools/bench_tdiv_med.c, so the
// benchmark times the same code.  Before including, define (or get from
// qs_impl.h): SIQS_FORCEINLINE, _trail_zcnt, SIQS_HAVE_SINV16 and the ISA
// macros (USE_AVX512BW, USE_AVX2, USE_SSE41, D_HAS_SSE2).  TDM_FB_T is the
// factor base type: anything with uint16_t* prime, root1, root2.

#ifndef TDIV_MED_OLD_H
#define TDIV_MED_OLD_H

#include <stdint.h>
#if !defined(__aarch64__)
#include <immintrin.h>
#endif

#ifndef TDM_FB_T
#define TDM_FB_T sieve_fb_compressed
#endif

// The old asm used LOC_BASE = 32768, which misses the hit whenever the
// progression lands exactly on the block end (stored root 0, test value p),
// i.e. about 1/p of the true hits per prime.  32767 makes q + 1 the exact
// ceiling: checked exhaustively for every fb prime < 2^13 with each shift
// make_fb_siqs can give it and every block_loc.  Define as 32768 to time
// or reproduce the original behavior.
#ifndef TDIV_OLD_LOC_BASE
#define TDIV_OLD_LOC_BASE 32767
#endif

// at most this many distinct fb primes > 256 can divide a Q of < 300 bits
#define TDIV_OLD_MAXHITS 64

#ifdef SIQS_HAVE_SINV16

// the shift regions of the 16-bit tables within [start, end): shift 24 in
// [start, b[0]), 26 in [b[0], b[1]), 28 in [b[1], end).  The shift steps
// at the first multiple of 8 whose prime exceeds 1024 (resp. 4096), as in
// make_fb_siqs.
static void tdm_old_bounds(const uint16_t* prime, uint32_t start, uint32_t end,
    uint32_t* b)
{
    uint32_t j = start & ~7u;

    while ((j < end) && (prime[j] <= 1024))
        j += 8;
    b[0] = (j < start) ? start : ((j > end) ? end : j);
    while ((j < end) && (prime[j] <= 4096))
        j += 8;
    b[1] = (j < start) ? start : ((j > end) ? end : j);
}

// one prime, the same 16-bit arithmetic as the vector lanes
#define TDM_OLD_SCALAR(j, xtra) \
    do { \
        uint16_t _t = (uint16_t)(TDIV_OLD_LOC_BASE - block_loc + corr[j]); \
        uint16_t _q = (uint16_t)((((uint32_t)_t * sinv[j]) >> 16) >> (xtra)); \
        uint16_t _y = (uint16_t)((uint32_t)(_q + 1) * fbc->prime[j] + block_loc - 32768); \
        if ((_y == fbc->root1[j]) || (_y == fbc->root2[j])) \
            buf[nb++] = (uint16_t)(j); \
    } while (0)

#if !defined(__aarch64__)

// test primes [i, end), all with shift 16 + xtra (a constant at each call),
// appending hit indices to buf; returns the new hit count
SIQS_FORCEINLINE uint32_t tdm_old_range(const TDM_FB_T* fbc,
    const uint16_t* sinv, const uint16_t* corr, uint32_t i, uint32_t end,
    uint32_t block_loc, const int xtra, uint16_t* buf, uint32_t nb)
{
#if defined(USE_AVX512BW)
    {
        __m512i vt = _mm512_set1_epi16((short)(TDIV_OLD_LOC_BASE - block_loc));
        __m512i vbl = _mm512_set1_epi16((short)(block_loc - 32768));
        __m512i vone = _mm512_set1_epi16(1);

        for (; i + 32 <= end; i += 32)
        {
            __m512i vp = _mm512_loadu_si512((const void*)(fbc->prime + i));
            __m512i t = _mm512_add_epi16(vt, _mm512_loadu_si512((const void*)(corr + i)));
            __m512i q = _mm512_srli_epi16(_mm512_mulhi_epu16(t,
                _mm512_loadu_si512((const void*)(sinv + i))), xtra);
            __m512i y = _mm512_add_epi16(_mm512_mullo_epi16(_mm512_add_epi16(q, vone), vp), vbl);
            uint32_t m = (uint32_t)(
                _mm512_cmpeq_epi16_mask(y, _mm512_loadu_si512((const void*)(fbc->root1 + i))) |
                _mm512_cmpeq_epi16_mask(y, _mm512_loadu_si512((const void*)(fbc->root2 + i))));

            while (m)
            {
                buf[nb++] = (uint16_t)(i + (uint32_t)_trail_zcnt(m));
                m &= m - 1;
            }
        }
    }
#endif

#if defined(USE_AVX2)
    {
        __m256i vt = _mm256_set1_epi16((short)(TDIV_OLD_LOC_BASE - block_loc));
        __m256i vbl = _mm256_set1_epi16((short)(block_loc - 32768));
        __m256i vone = _mm256_set1_epi16(1);

        for (; i + 16 <= end; i += 16)
        {
            __m256i vp = _mm256_loadu_si256((const __m256i*)(fbc->prime + i));
            __m256i t = _mm256_add_epi16(vt, _mm256_loadu_si256((const __m256i*)(corr + i)));
            __m256i q = _mm256_srli_epi16(_mm256_mulhi_epu16(t,
                _mm256_loadu_si256((const __m256i*)(sinv + i))), xtra);
            __m256i y = _mm256_add_epi16(_mm256_mullo_epi16(_mm256_add_epi16(q, vone), vp), vbl);
            __m256i h = _mm256_or_si256(
                _mm256_cmpeq_epi16(y, _mm256_loadu_si256((const __m256i*)(fbc->root1 + i))),
                _mm256_cmpeq_epi16(y, _mm256_loadu_si256((const __m256i*)(fbc->root2 + i))));
            // one movemask bit per byte: keep the high byte's bit of each lane
            uint32_t m = (uint32_t)_mm256_movemask_epi8(h) & 0xaaaaaaaau;

            while (m)
            {
                buf[nb++] = (uint16_t)(i + ((uint32_t)_trail_zcnt(m) >> 1));
                m &= m - 1;
            }
        }
    }
#endif

#if defined(USE_SSE41) || defined(D_HAS_SSE2)
    {
        // all SSE2 instructions
        __m128i vt = _mm_set1_epi16((short)(TDIV_OLD_LOC_BASE - block_loc));
        __m128i vbl = _mm_set1_epi16((short)(block_loc - 32768));
        __m128i vone = _mm_set1_epi16(1);

        for (; i + 8 <= end; i += 8)
        {
            __m128i vp = _mm_loadu_si128((const __m128i*)(fbc->prime + i));
            __m128i t = _mm_add_epi16(vt, _mm_loadu_si128((const __m128i*)(corr + i)));
            __m128i q = _mm_srli_epi16(_mm_mulhi_epu16(t,
                _mm_loadu_si128((const __m128i*)(sinv + i))), xtra);
            __m128i y = _mm_add_epi16(_mm_mullo_epi16(_mm_add_epi16(q, vone), vp), vbl);
            __m128i h = _mm_or_si128(
                _mm_cmpeq_epi16(y, _mm_loadu_si128((const __m128i*)(fbc->root1 + i))),
                _mm_cmpeq_epi16(y, _mm_loadu_si128((const __m128i*)(fbc->root2 + i))));
            uint32_t m = (uint32_t)_mm_movemask_epi8(h) & 0xaaaau;

            while (m)
            {
                buf[nb++] = (uint16_t)(i + ((uint32_t)_trail_zcnt(m) >> 1));
                m &= m - 1;
            }
        }
    }
#endif

    for (; i < end; i++)
        TDM_OLD_SCALAR(i, xtra);

    return nb;
}

#endif // !__aarch64__
#endif // SIQS_HAVE_SINV16

#endif
