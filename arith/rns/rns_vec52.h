/*
 * 52-bit vector products for the RNS kernels: LO52(acc, a, b) / HI52(acc,
 * a, b) add the low / high 52 bits of lo52(a) * lo52(b) to acc, exactly
 * like madd52lo/madd52hi (which ignore bits 52..63 of a and b).
 *
 * With IFMA they are the instructions.  Without (AVX-512F only) they are
 * emulated with FP64 FMAs, exact for any inputs: the low 52 bits of a and b
 * convert to double exactly, x*y < 2^104, and
 *   t  = fma_rz(x, y, 2^104)           bits(t)  = 2^104's bits + hi
 *   t' = fma(x, y, 2^104 + 2^52 - t)   bits(t') = 2^52's bits + lo
 * A paired LO52/HI52 on the same operands shares the conversions and the
 * first FMA after inlining (the compiler merges the identical pure
 * expressions).
 */

#ifndef RNS_VEC52_H
#define RNS_VEC52_H

#include "rns.h"

#ifdef RNS_HAVE_AVX512
#include <immintrin.h>

#ifdef RNS_HAVE_IFMA

#define LO52(acc, a, b) _mm512_madd52lo_epu64((acc), (a), (b))
#define HI52(acc, a, b) _mm512_madd52hi_epu64((acc), (a), (b))

#else

#define RNS_FP_RZ (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)
#define RNS_FP_BIAS_HI 0x4670000000000000ULL    // bits of 2^104
#define RNS_FP_BIAS_LO 0x4330000000000000ULL    // bits of 2^52

// low 52 bits of x, as an exact double
RNS_INLINE __m512d rns_cvt52(__m512i x)
{
    return _mm512_cvtepu64_pd(_mm512_and_si512(x,
        _mm512_set1_epi64((long long)0xfffffffffffffULL)));
}

// 2^104 + floor(x*y / 2^52) * 2^52  (x, y < 2^52)
RNS_INLINE __m512d rns_fp_hi(__m512d x, __m512d y)
{
    return _mm512_fmadd_round_pd(x, y,
        _mm512_castsi512_pd(_mm512_set1_epi64((long long)RNS_FP_BIAS_HI)), RNS_FP_RZ);
}

RNS_INLINE __m512i rns_madd52hi(__m512i acc, __m512i a, __m512i b)
{
    __m512d t = rns_fp_hi(rns_cvt52(a), rns_cvt52(b));
    return _mm512_add_epi64(acc, _mm512_sub_epi64(_mm512_castpd_si512(t),
        _mm512_set1_epi64((long long)RNS_FP_BIAS_HI)));
}

RNS_INLINE __m512i rns_madd52lo(__m512i acc, __m512i a, __m512i b)
{
    __m512d x = rns_cvt52(a), y = rns_cvt52(b);
    __m512d t = rns_fp_hi(x, y);
    t = _mm512_fmadd_round_pd(x, y, _mm512_sub_pd(_mm512_castsi512_pd(
        _mm512_set1_epi64((long long)(RNS_FP_BIAS_HI + 1))), t), RNS_FP_RZ);
    return _mm512_add_epi64(acc, _mm512_sub_epi64(_mm512_castpd_si512(t),
        _mm512_set1_epi64((long long)RNS_FP_BIAS_LO)));
}

#define LO52(acc, a, b) rns_madd52lo((acc), (a), (b))
#define HI52(acc, a, b) rns_madd52hi((acc), (a), (b))

/*
 * Fast path for the product loops: operands already converted (rns_cvt52),
 * both halves added as raw biased bit patterns:
 *   l += RNS_FP_BIAS_LO + lo,  h += RNS_FP_BIAS_HI + hi
 * The caller removes the biases once per accumulator with rns_fp_unbias
 * (counts of lo and hi products added to it; arithmetic mod 2^64).
 */
#define RNS_FP_MAC(l, h, x, y) do {                                         \
        __m512d _t = rns_fp_hi((x), (y));                                   \
        (h) = _mm512_add_epi64((h), _mm512_castpd_si512(_t));               \
        _t = _mm512_fmadd_round_pd((x), (y), _mm512_sub_pd(_mm512_castsi512_pd( \
            _mm512_set1_epi64((long long)(RNS_FP_BIAS_HI + 1))), _t), RNS_FP_RZ); \
        (l) = _mm512_add_epi64((l), _mm512_castpd_si512(_t));               \
    } while (0)

RNS_INLINE __m512i rns_fp_unbias(__m512i acc, int nlo, int nhi)
{
    uint64_t b = (uint64_t)(unsigned)nlo * RNS_FP_BIAS_LO +
        (uint64_t)(unsigned)nhi * RNS_FP_BIAS_HI;
    return _mm512_sub_epi64(acc, _mm512_set1_epi64((long long)b));
}

// broadcast lo52(c) as an exact double
RNS_INLINE __m512d rns_bcd52(uint64_t c)
{
    return _mm512_set1_pd((double)(c & 0xfffffffffffffULL));
}

#endif
#endif
#endif
