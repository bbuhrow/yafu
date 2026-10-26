/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// Shared helpers for the AArch64 NEON SIQS kernels: replacements for the
// x86 instructions NEON lacks (movemask, pmulhuw).  Include only under
// #if defined(__aarch64__) && defined(__ARM_NEON).

#ifndef SIQS_NEON_UTIL_H
#define SIQS_NEON_UTIL_H

#include <stdint.h>
#include <arm_neon.h>

// 64-bit movemask of four byte-compare results (lanes 0x00 or 0xff):
// bit i of the result is set iff byte i of m0:m1:m2:m3 is set, matching
// x86 pmovmskb on 64 consecutive bytes.  Each lane is masked to its bit
// position within its group of 8, and three rounds of pairwise adds fold
// 64 bytes down to 8.
static inline uint64_t neon_movemask64_u8(uint8x16_t m0, uint8x16_t m1,
    uint8x16_t m2, uint8x16_t m3)
{
    static const uint8_t bitpos[16] = {
        1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128 };
    const uint8x16_t b = vld1q_u8(bitpos);
    uint8x16_t s0 = vpaddq_u8(vandq_u8(m0, b), vandq_u8(m1, b));
    uint8x16_t s1 = vpaddq_u8(vandq_u8(m2, b), vandq_u8(m3, b));
    s0 = vpaddq_u8(s0, s1);
    s0 = vpaddq_u8(s0, s0);
    return vgetq_lane_u64(vreinterpretq_u64_u8(s0), 0);
}

// 16-byte compare result to a 64-bit "nibble mask": 4 bits per byte, all
// set iff the byte is set.  Cheaper than an exact movemask (one shrn);
// byte index = ctz(mask) / 4.
static inline uint64_t neon_nibblemask_u8(uint8x16_t m)
{
    return vget_lane_u64(vreinterpret_u64_u8(
        vshrn_n_u16(vreinterpretq_u16_u8(m), 4)), 0);
}

// unsigned 16x16 -> high 16 bits, 8 lanes (x86 pmulhuw)
static inline uint16x8_t neon_mulhi_u16(uint16x8_t a, uint16x8_t b)
{
    uint32x4_t lo = vmull_u16(vget_low_u16(a), vget_low_u16(b));
    uint32x4_t hi = vmull_high_u16(a, b);
    return vuzp2q_u16(vreinterpretq_u16_u32(lo), vreinterpretq_u16_u32(hi));
}

#endif
