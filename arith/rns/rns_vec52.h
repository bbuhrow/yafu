/*
 * Vector layer for the RNS kernels: 64-bit-lane integer ops, and 52-bit
 * products.  A kernel file selects the width before including this:
 *   RNS_VW = 512 (default): AVX-512, 8 lanes (needs RNS_HAVE_AVX512)
 *   RNS_VW = 256: AVX2 + FMA, 4 lanes (needs RNS_HAVE_AVX2)
 * Types rns_vi / rns_vd, ops V_*.
 *
 * Products: LO52(acc, a, b) / HI52(acc, a, b) add the low / high 52 bits
 * of lo52(a) * lo52(b) to acc, exactly like madd52lo/madd52hi (which ignore
 * bits 52..63 of a and b).  With IFMA (512 only) they are the
 * instructions.  Otherwise (RNS_VEC_FP) they are emulated exactly with FP64
 * FMAs: the low 52 bits of a and b convert to double exactly, x*y < 2^104,
 *   t  = fma_rz(x, y, 2^104)           bits(t)  = 2^104's bits + hi
 *   t' = fma(x, y, 2^104 + 2^52 - t)   bits(t') = 2^52's bits + lo
 * 512 uses embedded round-to-zero; 256 has none, so its callers set MXCSR
 * to round toward zero around each kernel call (rns_fp_rz_enter/leave; the
 * kernel files are compiled with -frounding-math).  Every other FP op is
 * exact.  A paired LO52/HI52 on the same operands shares the conversions
 * and the first FMA after inlining.
 */

#ifndef RNS_VEC52_H
#define RNS_VEC52_H

#include "rns.h"

#ifndef RNS_VW
#define RNS_VW 512
#endif

#if (RNS_VW == 512 && defined(RNS_HAVE_AVX512)) || (RNS_VW == 256 && defined(RNS_HAVE_AVX2))
#define RNS_VEC_OK 1
#include <immintrin.h>

#define RNS_FP_BIAS_HI 0x4670000000000000ULL    // bits of 2^104
#define RNS_FP_BIAS_LO 0x4330000000000000ULL    // bits of 2^52

#if RNS_VW == 512

typedef __m512i rns_vi;
typedef __m512d rns_vd;
#define V_LANES         8
#define V_ADD(a, b)     _mm512_add_epi64((a), (b))
#define V_SUB(a, b)     _mm512_sub_epi64((a), (b))
#define V_AND(a, b)     _mm512_and_si512((a), (b))
#define V_OR(a, b)      _mm512_or_si512((a), (b))
#define V_SRLI(a, n)    _mm512_srli_epi64((a), (n))
#define V_SLLV(a, n)    _mm512_sllv_epi64((a), (n))
#define V_SRLV(a, n)    _mm512_srlv_epi64((a), (n))
#define V_SRAI52(a)     _mm512_srai_epi64((a), 52)
#define V_MULU32(a, b)  _mm512_mul_epu32((a), (b))
#define V_BC(x)         _mm512_set1_epi64((long long)(x))
#define V_ZERO          _mm512_setzero_si512()
#define V_LD(p)         _mm512_load_si512((const void *)(p))
#define V_ST(p, v)      _mm512_store_si512((void *)(p), (v))
#define V_LDD(p)        _mm512_load_pd((const void *)(p))
#define V_STD(p, v)     _mm512_store_pd((void *)(p), (v))
#define V_BCD(x)        _mm512_set1_pd(x)
#define V_SUBD(a, b)    _mm512_sub_pd((a), (b))
#define V_I2D(a)        _mm512_castsi512_pd(a)
#define V_D2I(a)        _mm512_castpd_si512(a)
// any lane >= 2^52
#define V_ANY_GE52(a)   (_mm512_cmpge_epu64_mask((a), V_BC(1ULL << 52)) != 0)
// per lane: tab[idx], tab = 16 words (aligned), idx < 16
#define V_LOOKUP16(tab, idx) _mm512_permutex2var_epi64(V_LD(tab), (idx), V_LD((tab) + 8))
#define V_FMA_RZ(x, y, z) _mm512_fmadd_round_pd((x), (y), (z), (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC))
#define V_CVT52(a)      _mm512_cvtepu64_pd(a)   // a < 2^52
#ifdef RNS_HAVE_IFMA
#define RNS_VEC_FP 0
#else
#define RNS_VEC_FP 1
#endif

#else   // 256

typedef __m256i rns_vi;
typedef __m256d rns_vd;
#define V_LANES         4
#define V_ADD(a, b)     _mm256_add_epi64((a), (b))
#define V_SUB(a, b)     _mm256_sub_epi64((a), (b))
#define V_AND(a, b)     _mm256_and_si256((a), (b))
#define V_OR(a, b)      _mm256_or_si256((a), (b))
#define V_SRLI(a, n)    _mm256_srli_epi64((a), (n))
#define V_SLLV(a, n)    _mm256_sllv_epi64((a), (n))
#define V_SRLV(a, n)    _mm256_srlv_epi64((a), (n))
// no vpsraq: floor(a / 2^52) for |a| < 2^62 via a biased logical shift
#define V_SRAI52(a)     V_SUB(V_SRLI(V_ADD((a), V_BC(1ULL << 62)), 52), V_BC(1ULL << 10))
#define V_MULU32(a, b)  _mm256_mul_epu32((a), (b))
#define V_BC(x)         _mm256_set1_epi64x((long long)(x))
#define V_ZERO          _mm256_setzero_si256()
#define V_LD(p)         _mm256_load_si256((const __m256i *)(const void *)(p))
#define V_ST(p, v)      _mm256_store_si256((__m256i *)(void *)(p), (v))
#define V_LDD(p)        _mm256_load_pd((const double *)(const void *)(p))
#define V_STD(p, v)     _mm256_store_pd((double *)(void *)(p), (v))
#define V_BCD(x)        _mm256_set1_pd(x)
#define V_SUBD(a, b)    _mm256_sub_pd((a), (b))
#define V_I2D(a)        _mm256_castsi256_pd(a)
#define V_D2I(a)        _mm256_castpd_si256(a)
RNS_INLINE int rns_any_ge52_256(__m256i a)
{
    __m256i h = _mm256_srli_epi64(a, 52);
    return !_mm256_testz_si256(h, h);
}
#define V_ANY_GE52(a)   rns_any_ge52_256(a)
#define V_LOOKUP16(tab, idx) _mm256_i64gather_epi64((const long long *)(const void *)(tab), (idx), 8)
// MXCSR is round-toward-zero inside the kernels (rns_fp_rz_enter)
#define V_FMA_RZ(x, y, z) _mm256_fmadd_pd((x), (y), (z))
// a < 2^52: bits(a | 2^52's bits) - 2^52, exact (no vcvtuqq2pd in AVX2)
#define V_CVT52(a)      V_SUBD(V_I2D(V_OR((a), V_BC(RNS_FP_BIAS_LO))), V_I2D(V_BC(RNS_FP_BIAS_LO)))
#define RNS_VEC_FP 1

RNS_INLINE unsigned int rns_fp_rz_enter(void)
{
    unsigned int csr = _mm_getcsr();
    _mm_setcsr((csr & ~0x6000u) | 0x6000u);
    return csr;
}
RNS_INLINE void rns_fp_rz_leave(unsigned int csr)
{
    _mm_setcsr(csr);
}

#endif

#if !RNS_VEC_FP

#define LO52(acc, a, b) _mm512_madd52lo_epu64((acc), (a), (b))
#define HI52(acc, a, b) _mm512_madd52hi_epu64((acc), (a), (b))

#else

// low 52 bits of x, as an exact double
RNS_INLINE rns_vd rns_cvt52(rns_vi x)
{
    return V_CVT52(V_AND(x, V_BC(0xfffffffffffffULL)));
}

// 2^104 + floor(x*y / 2^52) * 2^52  (x, y < 2^52)
RNS_INLINE rns_vd rns_fp_hi(rns_vd x, rns_vd y)
{
    return V_FMA_RZ(x, y, V_I2D(V_BC(RNS_FP_BIAS_HI)));
}

RNS_INLINE rns_vi rns_madd52hi(rns_vi acc, rns_vi a, rns_vi b)
{
    rns_vd t = rns_fp_hi(rns_cvt52(a), rns_cvt52(b));
    return V_ADD(acc, V_SUB(V_D2I(t), V_BC(RNS_FP_BIAS_HI)));
}

RNS_INLINE rns_vi rns_madd52lo(rns_vi acc, rns_vi a, rns_vi b)
{
    rns_vd x = rns_cvt52(a), y = rns_cvt52(b);
    rns_vd t = rns_fp_hi(x, y);
    t = V_FMA_RZ(x, y, V_SUBD(V_I2D(V_BC(RNS_FP_BIAS_HI + 1)), t));
    return V_ADD(acc, V_SUB(V_D2I(t), V_BC(RNS_FP_BIAS_LO)));
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
        rns_vd _t = rns_fp_hi((x), (y));                                    \
        (h) = V_ADD((h), V_D2I(_t));                                        \
        _t = V_FMA_RZ((x), (y), V_SUBD(V_I2D(V_BC(RNS_FP_BIAS_HI + 1)), _t)); \
        (l) = V_ADD((l), V_D2I(_t));                                        \
    } while (0)

RNS_INLINE rns_vi rns_fp_unbias(rns_vi acc, int nlo, int nhi)
{
    uint64_t b = (uint64_t)(unsigned)nlo * RNS_FP_BIAS_LO +
        (uint64_t)(unsigned)nhi * RNS_FP_BIAS_HI;
    return V_SUB(acc, V_BC(b));
}

// broadcast lo52(c) as an exact double
RNS_INLINE rns_vd rns_bcd52(uint64_t c)
{
    return V_BCD((double)(c & 0xfffffffffffffULL));
}

#endif
#endif
#endif
