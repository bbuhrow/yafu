// Scalar emulation of the AArch64 NEON intrinsics used by
// direct_sieve_neon.c, for testing that backend on x86 (-DSIQS_NEON_EMU).
// Semantics follow the Arm ACLE definitions; not a general-purpose header.

#ifndef ARM_NEON_EMU_H
#define ARM_NEON_EMU_H

#include <stdint.h>

typedef struct { uint32_t v[4]; } uint32x4_t;
typedef struct { uint8_t v[16]; } uint8x16_t;

static inline uint32x4_t vld1q_u32(const uint32_t* p)
{
    uint32x4_t r; int i; for (i = 0; i < 4; i++) r.v[i] = p[i]; return r;
}
static inline void vst1q_u32(uint32_t* p, uint32x4_t a)
{
    int i; for (i = 0; i < 4; i++) p[i] = a.v[i];
}
static inline uint8x16_t vld1q_u8(const uint8_t* p)
{
    uint8x16_t r; int i; for (i = 0; i < 16; i++) r.v[i] = p[i]; return r;
}
static inline uint32x4_t vdupq_n_u32(uint32_t x)
{
    uint32x4_t r; int i; for (i = 0; i < 4; i++) r.v[i] = x; return r;
}

#define NEON_EMU_BIN(name, expr) \
static inline uint32x4_t name(uint32x4_t a, uint32x4_t b) \
{ uint32x4_t r; int i; for (i = 0; i < 4; i++) r.v[i] = (expr); return r; }

NEON_EMU_BIN(vaddq_u32, a.v[i] + b.v[i])
NEON_EMU_BIN(vsubq_u32, a.v[i] - b.v[i])
NEON_EMU_BIN(vandq_u32, a.v[i] & b.v[i])
NEON_EMU_BIN(vorrq_u32, a.v[i] | b.v[i])
NEON_EMU_BIN(vbicq_u32, a.v[i] & ~b.v[i])
NEON_EMU_BIN(vcltq_u32, (a.v[i] < b.v[i]) ? 0xffffffffu : 0)
NEON_EMU_BIN(vceqq_u32, (a.v[i] == b.v[i]) ? 0xffffffffu : 0)
NEON_EMU_BIN(vminq_u32, (a.v[i] < b.v[i]) ? a.v[i] : b.v[i])
NEON_EMU_BIN(vmaxq_u32, (a.v[i] > b.v[i]) ? a.v[i] : b.v[i])
#undef NEON_EMU_BIN

// bitwise select: (m & b) | (~m & c)
static inline uint32x4_t vbslq_u32(uint32x4_t m, uint32x4_t b, uint32x4_t c)
{
    uint32x4_t r; int i;
    for (i = 0; i < 4; i++) r.v[i] = (m.v[i] & b.v[i]) | (~m.v[i] & c.v[i]);
    return r;
}

static inline uint32_t vaddvq_u32(uint32x4_t a)
{
    return a.v[0] + a.v[1] + a.v[2] + a.v[3];
}

// little-endian lane views
static inline uint8x16_t vreinterpretq_u8_u32(uint32x4_t a)
{
    uint8x16_t r; int i;
    for (i = 0; i < 16; i++) r.v[i] = (uint8_t)(a.v[i >> 2] >> (8 * (i & 3)));
    return r;
}
static inline uint32x4_t vreinterpretq_u32_u8(uint8x16_t a)
{
    uint32x4_t r; int i;
    for (i = 0; i < 4; i++)
        r.v[i] = (uint32_t)a.v[4 * i] | ((uint32_t)a.v[4 * i + 1] << 8) |
            ((uint32_t)a.v[4 * i + 2] << 16) | ((uint32_t)a.v[4 * i + 3] << 24);
    return r;
}

// table lookup: indices >= 16 give 0
static inline uint8x16_t vqtbl1q_u8(uint8x16_t t, uint8x16_t idx)
{
    uint8x16_t r; int i;
    for (i = 0; i < 16; i++) r.v[i] = (idx.v[i] < 16) ? t.v[idx.v[i]] : 0;
    return r;
}

#endif
