/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// NEON (AArch64) backend for the 4-lane large-prime kernels
// (direct_sieve_v4.h).  Compiled on AArch64 builds with __ARM_NEON.
//
// For testing on x86, build with -DSIQS_NEON_EMU: the same backend is then
// compiled against arm_neon_emu.h, a scalar emulation of the handful of
// NEON intrinsics used here, and -siqsLPisa neon selects it, so the NEON
// backend can be checked bit-for-bit against the other kernels.

#include <stdint.h>
#include <string.h>
#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "tdiv_macros_common.h"

#if defined(TRY_DIRECT_SIEVE) && (defined(SIQS_NEON_EMU) || (defined(__aarch64__) && defined(__ARM_NEON)))

#if defined(SIQS_NEON_EMU)
#include "arm_neon_emu.h"
#else
#include <arm_neon.h>
#endif

typedef uint32x4_t v4t;

static inline uint32x4_t v4_setr_neon(uint32_t a, uint32_t b, uint32_t c, uint32_t d)
{
    uint32_t t[4] = { a, b, c, d };
    return vld1q_u32(t);
}

static inline uint32_t v4_mask_neon(uint32x4_t m)
{
    static const uint32_t bits[4] = { 1, 2, 4, 8 };
    return vaddvq_u32(vandq_u32(m, vld1q_u32(bits)));
}

#define V4_LOAD(p)          vld1q_u32((const uint32_t*)(p))
#define V4_STORE(p, v)      vst1q_u32((uint32_t*)(p), (v))
#define V4_SET1(x)          vdupq_n_u32((uint32_t)(x))
#define V4_SETR(a, b, c, d) v4_setr_neon(a, b, c, d)
#define V4_ADD(a, b)        vaddq_u32(a, b)
#define V4_SUB(a, b)        vsubq_u32(a, b)
#define V4_AND(a, b)        vandq_u32(a, b)
#define V4_OR(a, b)         vorrq_u32(a, b)
#define V4_ANDNOT(m, x)     vbicq_u32(x, m)
#define V4_LT(a, b)         vcltq_u32(a, b)
#define V4_EQ(a, b)         vceqq_u32(a, b)
#define V4_MASK(m)          v4_mask_neon(m)
#define V4_BLEND(a, b, m)   vbslq_u32(m, b, a)
#define V4_MIN(a, b)        vminq_u32(a, b)
#define V4_MAX(a, b)        vmaxq_u32(a, b)
#define V4_COMPRESS(v, m)   vreinterpretq_u32_u8(vqtbl1q_u8(vreinterpretq_u8_u32(v), \
                                vld1q_u8(V4_FN(v4_cmp_lut)[m])))
#define V4_FN(name)         name##_neon
#if !defined(SIQS_NEON_EMU)
// the NEON bucket fill (update_poly_roots_32k_neon.c), AArch64 only
#define V4_LP_BUCKETS       nextRoots_32k_lp_buckets_neon
#endif

#include "direct_sieve_v4.h"

#endif
