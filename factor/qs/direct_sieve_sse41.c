/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// SSE4.1 backend for the 4-lane large-prime kernels (direct_sieve_v4.h).
// Compiled into SSE4.1, AVX2 and AVX-512 builds.

#include <stdint.h>
#include <string.h>
#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "tdiv_macros_common.h"

#if defined(TRY_DIRECT_SIEVE) && defined(USE_SSE41)
#include <immintrin.h>

typedef __m128i v4t;

#define V4_LOAD(p)          _mm_loadu_si128((const __m128i*)(p))
#define V4_STORE(p, v)      _mm_storeu_si128((__m128i*)(p), (v))
#define V4_SET1(x)          _mm_set1_epi32((int)(x))
#define V4_SETR(a, b, c, d) _mm_setr_epi32(a, b, c, d)
#define V4_ADD(a, b)        _mm_add_epi32(a, b)
#define V4_SUB(a, b)        _mm_sub_epi32(a, b)
#define V4_AND(a, b)        _mm_and_si128(a, b)
#define V4_OR(a, b)         _mm_or_si128(a, b)
#define V4_ANDNOT(m, x)     _mm_andnot_si128(m, x)
#define V4_LT(a, b)         _mm_cmpgt_epi32(b, a)
#define V4_EQ(a, b)         _mm_cmpeq_epi32(a, b)
#define V4_MASK(m)          ((uint32_t)_mm_movemask_ps(_mm_castsi128_ps(m)))
#define V4_BLEND(a, b, m)   _mm_blendv_epi8(a, b, m)
#define V4_MIN(a, b)        _mm_min_epu32(a, b)
#define V4_MAX(a, b)        _mm_max_epu32(a, b)
#define V4_COMPRESS(v, m)   _mm_shuffle_epi8(v, _mm_loadu_si128((const __m128i*)V4_FN(v4_cmp_lut)[m]))
#define V4_FN(name)         name##_sse41

#include "direct_sieve_v4.h"

#endif
