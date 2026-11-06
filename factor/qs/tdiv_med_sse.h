/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// 128-bit helpers for the Lemire divisibility test in
// tdiv_medprimes_32k_new (tdiv_med_32k.c): SSE4.1 instructions when
// USE_SSE41 is defined, SSE2 equivalents otherwise.  Also used by
// tools/bench_tdiv_med.c, so the benchmark times the same code.

#ifndef TDIV_MED_SSE_H
#define TDIV_MED_SSE_H

#include <stdint.h>
#include <immintrin.h>

// 4 x uint16 at p, zero-extended to 4 x uint32 (reads only 8 bytes)
static inline __m128i tdm_load4_u16(const uint16_t* p)
{
    __m128i v = _mm_loadl_epi64((const __m128i*)p);
#ifdef USE_SSE41
    return _mm_cvtepu16_epi32(v);
#else
    return _mm_unpacklo_epi16(v, _mm_setzero_si128());
#endif
}

// low 32 bits of a * b, per lane
static inline __m128i tdm_mullo_epi32(__m128i a, __m128i b)
{
#ifdef USE_SSE41
    return _mm_mullo_epi32(a, b);
#else
    __m128i e = _mm_mul_epu32(a, b);                                     // lanes 0, 2
    __m128i o = _mm_mul_epu32(_mm_srli_epi64(a, 32), _mm_srli_epi64(b, 32)); // lanes 1, 3
    return _mm_unpacklo_epi32(_mm_shuffle_epi32(e, _MM_SHUFFLE(0, 0, 2, 0)),
        _mm_shuffle_epi32(o, _MM_SHUFFLE(0, 0, 2, 0)));
#endif
}

// all-ones in lanes where x <= t (unsigned)
static inline __m128i tdm_cmple_epu32(__m128i x, __m128i t)
{
#ifdef USE_SSE41
    return _mm_cmpeq_epi32(_mm_min_epu32(x, t), x);
#else
    // bias by 2^31 so a signed compare orders unsigned values
    const __m128i vsign = _mm_set1_epi32((int)0x80000000);
    __m128i gt = _mm_cmpgt_epi32(_mm_xor_si128(x, vsign), _mm_xor_si128(t, vsign));
    return _mm_xor_si128(gt, _mm_set1_epi32(-1));
#endif
}

#endif
