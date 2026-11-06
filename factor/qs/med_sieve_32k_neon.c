/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// AArch64 NEON med sieve (primes from start_prime to med_B, 32k blocks).
//
// med_sieveblock_32k_neon follows med_sieveblock_32k_sse41: scalar loops
// for p < 2^13, then 8 primes at a time (16-bit roots) in ranges where
// the number of hits per root is bounded:
//   8192  < p < 32768/3: 3 unconditional steps + 1 conditional
//   32768/3 < p < 2^14:  2 + 1
//   2^14  < p < 2^15:    1 + 1
//   2^15  < p < med_B:   1 conditional
// with scalar crossover loops between ranges.  Ranges need not start at a
// multiple of 8 (unaligned loads, scalar tails).
//
// a-primes above 8192 have prime = root = 0 (set_aprime_roots), so in the
// vector ranges they only hit sieve[0], which is never reported.
//
// logp: like the x86 vector sieves, one value per range (13/14/15).  Build
// with -DMED_NEON_EXACT_LOGP to use each prime's own fb->logp instead; the
// sieve is then bit-identical to med_sieveblock_32k's (except sieve[0]), so
// a whole run matches -noneon exactly (tools/neon_ab.sh).

#include "qs_impl.h"
#include "sieve_macros_32k.h"

#if defined(__aarch64__) && defined(__ARM_NEON)
#include "neon_util.h"

#ifdef MED_NEON_EXACT_LOGP
#define MED_NEON_SET_LOGP(c) \
    { uint32_t _k; for (_k = 0; _k < 8; _k++) lp[_k] = (uint8_t)fb->logp[i + _k]; }
#define MLP(k) lp[k]
#else
#define MED_NEON_SET_LOGP(c) logp = (c);
#define MLP(k) logp
#endif

// both roots of all 8 lanes hit; advance by p
#define MED_NEON_STEP \
    vst1q_u16(r1, vr1); \
    vst1q_u16(r2, vr2); \
    sieve[r1[0]] -= MLP(0); sieve[r2[0]] -= MLP(0); \
    sieve[r1[1]] -= MLP(1); sieve[r2[1]] -= MLP(1); \
    sieve[r1[2]] -= MLP(2); sieve[r2[2]] -= MLP(2); \
    sieve[r1[3]] -= MLP(3); sieve[r2[3]] -= MLP(3); \
    sieve[r1[4]] -= MLP(4); sieve[r2[4]] -= MLP(4); \
    sieve[r1[5]] -= MLP(5); sieve[r2[5]] -= MLP(5); \
    sieve[r1[6]] -= MLP(6); sieve[r2[6]] -= MLP(6); \
    sieve[r1[7]] -= MLP(7); sieve[r2[7]] -= MLP(7); \
    vr1 = vaddq_u16(vr1, vp); \
    vr2 = vaddq_u16(vr2, vp);

// lanes whose root is still in the block take a last hit (root1 <= root2,
// so a root2 hit implies a root1 hit); then advance the hit roots by p,
// subtract the blocksize, and re-sort so root1 <= root2 for the next block
#define MED_NEON_LANE(k) \
    if ((b1 >> (8 * k)) & 1) { \
        sieve[r1[k]] -= MLP(k); \
        if ((b2 >> (8 * k)) & 1) \
            sieve[r2[k]] -= MLP(k); \
    }

#define MED_NEON_FINAL_AND_STORE \
    m1 = vcleq_u16(vr1, vlast); \
    m2 = vcleq_u16(vr2, vlast); \
    vst1q_u16(r1, vr1); \
    vst1q_u16(r2, vr2); \
    b1 = vget_lane_u64(vreinterpret_u64_u8(vmovn_u16(m1)), 0); \
    b2 = vget_lane_u64(vreinterpret_u64_u8(vmovn_u16(m2)), 0); \
    MED_NEON_LANE(0) MED_NEON_LANE(1) MED_NEON_LANE(2) MED_NEON_LANE(3) \
    MED_NEON_LANE(4) MED_NEON_LANE(5) MED_NEON_LANE(6) MED_NEON_LANE(7) \
    vr1 = vsubq_u16(vaddq_u16(vr1, vandq_u16(m1, vp)), vblock); \
    vr2 = vsubq_u16(vaddq_u16(vr2, vandq_u16(m2, vp)), vblock); \
    vst1q_u16(fb->root1 + i, vminq_u16(vr1, vr2)); \
    vst1q_u16(fb->root2 + i, vmaxq_u16(vr1, vr2));

// one prime, any number of hits (crossover between vector ranges)
#define MED_NEON_SCALAR_1X \
    prime = fb->prime[i]; \
    root1 = fb->root1[i]; \
    root2 = fb->root2[i]; \
    logp = (uint8_t)fb->logp[i]; \
    if (prime == 0) \
        continue; \
    SIEVE_1X; \
    SIEVE_LAST; \
    UPDATE_ROOTS;

#define MED_NEON_LOAD \
    vp = vld1q_u16(fb->prime + i); \
    vr1 = vld1q_u16(fb->root1 + i); \
    vr2 = vld1q_u16(fb->root2 + i);

void med_sieveblock_32k_neon(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    uint32_t i;
    uint32_t med_B = full_fb->med_B;
    uint32_t prime, root1, root2, tmp, stop;
    uint8_t logp;
#ifdef MED_NEON_EXACT_LOGP
    uint8_t lp[8];
#endif
    uint16_t r1[8], r2[8];
    uint64_t b1, b2;
    uint16x8_t vp, vr1, vr2, m1, m2;
    const uint16x8_t vlast = vdupq_n_u16(32767);
    const uint16x8_t vblock = vdupq_n_u16(32768);
    // full fb primes: unaffected by the a-prime zeroing in fb->prime
    const uint32_t* plist = full_fb->list->prime;
    // range starts, clamped to med_B (small inputs can stop early).  Each
    // is the first multiple-of-8 index past its threshold, so the group
    // just below it can hold primes above the range's limit: the vector
    // loops check their limit on the group's last prime instead.
    uint32_t b13 = MIN(full_fb->fb_13bit_B, med_B);
    uint32_t b3 = MIN(full_fb->fb_32k_div3, med_B);
    uint32_t b14 = MIN(full_fb->fb_14bit_B, med_B);
    uint32_t b15 = MIN(full_fb->fb_15bit_B, med_B);

    // p < 2^13: many hits per root, scalar loops (as on x86)
    for (i = start_prime; i < b13; i++)
    {
        uint8_t* s2;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = (uint8_t)fb->logp[i];

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    // 8192 < p <= 32768/3: 3 steps stay in the block (r + 2p < 3p < 32768)
    for (; ((i + 8) <= b3) && (plist[i + 7] <= 10922); i += 8)
    {
        MED_NEON_SET_LOGP(13);
        MED_NEON_LOAD;
        MED_NEON_STEP;
        MED_NEON_STEP;
        MED_NEON_STEP;
        MED_NEON_FINAL_AND_STORE;
    }

    for (; i < b3; i++)
    {
        MED_NEON_SCALAR_1X;
    }

    // 32768/3 < p <= 2^14: 2 steps + 1 conditional
    for (; ((i + 8) <= b14) && (plist[i + 7] <= 16384); i += 8)
    {
        MED_NEON_SET_LOGP(14);
        MED_NEON_LOAD;
        MED_NEON_STEP;
        MED_NEON_STEP;
        MED_NEON_FINAL_AND_STORE;
    }

    for (; i < b14; i++)
    {
        MED_NEON_SCALAR_1X;
    }

    // 2^14 < p <= 2^15: 1 step + 1 conditional
    for (; ((i + 8) <= b15) && (plist[i + 7] <= 32768); i += 8)
    {
        MED_NEON_SET_LOGP(15);
        MED_NEON_LOAD;
        MED_NEON_STEP;
        MED_NEON_FINAL_AND_STORE;
    }

    for (; (i < med_B) && (plist[i] <= 32768); i++)
    {
        MED_NEON_SCALAR_1X;
    }

    // 2^15 < p < med_B: at most one hit per root
    for (; (i + 8) <= med_B; i += 8)
    {
        MED_NEON_SET_LOGP(15);
        MED_NEON_LOAD;
        MED_NEON_FINAL_AND_STORE;
    }

    for (; i < med_B; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = (uint8_t)fb->logp[i];

        SIEVE_BIG;
        UPDATE_ROOTS;
    }

    return;
}

// 8-bit lane mask (bit 8k set for lane k) of a 16-bit compare result
#define MED_NEON_LANEBITS(m) \
    (vget_lane_u64(vreinterpret_u64_u8(vmovn_u16(m)), 0) & 0x0101010101010101ULL)

// sink strip for the fixed passes: the last MED_NEON_SINK_BYTES of the
// previous block, which is finished before this one is sieved (the leading
// guard for the first block; see MED_SINK_STRIP in med_sieve_32k_avx2.c)
#define MED_NEON_SINK_BYTES 256

static const uint16_t med_neon_sink_lanes[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

// one fixed pass over an 8-prime group: lanes in mask mv whose root is in
// the block take a hit at sieve[root] and advance by p; the rest write the
// sink strip, so every pass does 16 RMWs and no branches
#define MED_NEON_FIXED_PASS(mv) \
    m1 = vandq_u16(mv, vcleq_u16(vr1, vlast)); \
    m2 = vandq_u16(mv, vcleq_u16(vr2, vlast)); \
    vst1q_u16(r1, vbslq_u16(m1, vaddq_u16(vr1, vsink_bias), vsink1)); \
    vst1q_u16(r2, vbslq_u16(m2, vaddq_u16(vr2, vsink_bias), vsink2)); \
    sbase[r1[0]] -= MLP(0); sbase[r2[0]] -= MLP(0); \
    sbase[r1[1]] -= MLP(1); sbase[r2[1]] -= MLP(1); \
    sbase[r1[2]] -= MLP(2); sbase[r2[2]] -= MLP(2); \
    sbase[r1[3]] -= MLP(3); sbase[r2[3]] -= MLP(3); \
    sbase[r1[4]] -= MLP(4); sbase[r2[4]] -= MLP(4); \
    sbase[r1[5]] -= MLP(5); sbase[r2[5]] -= MLP(5); \
    sbase[r1[6]] -= MLP(6); sbase[r2[6]] -= MLP(6); \
    sbase[r1[7]] -= MLP(7); sbase[r2[7]] -= MLP(7); \
    vsink1 = vandq_u16(vaddq_u16(vsink1, vsink_step), vsink_wrap); \
    vsink2 = vaddq_u16(vsink1, vsink_half); \
    vr1 = vaddq_u16(vr1, vandq_u16(vp, m1)); \
    vr2 = vaddq_u16(vr2, vandq_u16(vp, m2));

// one sink-strip pass over the selective lanes of the first part: lanes in
// mask m hit sieve[root] (sbase-relative, so biased), the rest write the
// strip; 16 RMWs and no branches
#define MED_NEON_SINK_RMW(m, va, vb) \
    vst1q_u16(r1, vbslq_u16(m, vaddq_u16(va, vsink_bias), vsink1)); \
    vst1q_u16(r2, vbslq_u16(m, vaddq_u16(vb, vsink_bias), vsink2)); \
    sbase[r1[0]] -= MLP(0); sbase[r2[0]] -= MLP(0); \
    sbase[r1[1]] -= MLP(1); sbase[r2[1]] -= MLP(1); \
    sbase[r1[2]] -= MLP(2); sbase[r2[2]] -= MLP(2); \
    sbase[r1[3]] -= MLP(3); sbase[r2[3]] -= MLP(3); \
    sbase[r1[4]] -= MLP(4); sbase[r2[4]] -= MLP(4); \
    sbase[r1[5]] -= MLP(5); sbase[r2[5]] -= MLP(5); \
    sbase[r1[6]] -= MLP(6); sbase[r2[6]] -= MLP(6); \
    sbase[r1[7]] -= MLP(7); sbase[r2[7]] -= MLP(7); \
    vsink1 = vandq_u16(vaddq_u16(vsink1, vsink_step), vsink_wrap); \
    vsink2 = vaddq_u16(vsink1, vsink_half);

// the same for one root per lane (the last root1 hits): 8 RMWs
#define MED_NEON_SINK_RMW8(m, va) \
    vst1q_u16(r1, vbslq_u16(m, vaddq_u16(va, vsink_bias), vsink1)); \
    sbase[r1[0]] -= MLP(0); sbase[r1[1]] -= MLP(1); \
    sbase[r1[2]] -= MLP(2); sbase[r1[3]] -= MLP(3); \
    sbase[r1[4]] -= MLP(4); sbase[r1[5]] -= MLP(5); \
    sbase[r1[6]] -= MLP(6); sbase[r1[7]] -= MLP(7); \
    vsink1 = vandq_u16(vaddq_u16(vsink1, vsink_step), vsink_wrap); \
    vsink2 = vaddq_u16(vsink1, vsink_half);

// NEON counterpart of med_sieveblock_32k_avx2_split:
//  - below fb_32k_div3, 8 primes at a time: dump both roots of every lane
//    while all valid lanes are in the block, then only the lanes still in
//    it, then a last root1 hit (logp of the group's first prime for all 8
//    lanes).  With small_scalar, primes below 2^13 use the scalar loops
//    instead (as med_sieveblock_32k_neon does) and the vector loop starts
//    at fb_13bit_B.  With sink_select, the "only the lanes still in it"
//    loops write every lane, missed lanes to the sink strip, instead of
//    looping over the set lanes (MED_SINK_STRIP in the AVX-512 version).
//  - above fb_32k_div3, groups of 8 with K fixed passes each, K from the
//    group's smallest prime: 32768/3 < p < 2^14: 3, 2^14 < p < 2^15: 2,
//    2^15 < p: 1 (lanes with p > 2^15 sit out the later passes, where a
//    wrapped root would look like another hit)
// Invalid lanes in the first part (a-primes: root 0xffff, or prime 0) are
// parked at 0xffff and dump to sieve[65535], as in the AVX2 version (guard
// space, or saved/restored by the caller with the direct layout).  The
// sink strip is sieve[-MED_NEON_SINK_BYTES, 0).
static inline __attribute__((always_inline)) void med_neon_split_body(
    uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, const int small_scalar, const int sink_select)
{
    uint32_t i;
    uint32_t med_B = full_fb->med_B;
    const uint32_t* plist = full_fb->list->prime;
    uint32_t prime, root1, root2, tmp, stop;
    uint32_t b13, b3, vbound;
    uint8_t logp;
#ifdef MED_NEON_EXACT_LOGP
    uint8_t lp[8];
#endif
    uint8_t* sbase = sieve - MED_NEON_SINK_BYTES;
    uint16_t r1[8], r2[8];
    uint64_t bm;
    uint16x8_t vp, vr1, vr2, m1, m2, vinit, vval1, vval2, nbig;
    const uint16x8_t vlast = vdupq_n_u16(32767);
    const uint16x8_t vblock = vdupq_n_u16(32768);
    const uint16x8_t vffff = vdupq_n_u16(0xffff);
    const uint16x8_t vsink_bias = vdupq_n_u16(MED_NEON_SINK_BYTES);
    const uint16x8_t vsink_half = vdupq_n_u16(8);
    const uint16x8_t vsink_step = vdupq_n_u16(16);
    const uint16x8_t vsink_wrap = vdupq_n_u16(MED_NEON_SINK_BYTES - 1);
    uint16x8_t vsink1 = vld1q_u16(med_neon_sink_lanes);
    uint16x8_t vsink2 = vaddq_u16(vsink1, vsink_half);

    b3 = MIN(full_fb->fb_32k_div3, med_B);
    b13 = MIN(full_fb->fb_13bit_B, b3);
    vbound = b3 & ~7u;

    i = start_prime;
    if (small_scalar)
    {
        // p < 2^13: scalar loops
        for (; i < b13; i++)
        {
            uint8_t* s2;

            prime = fb->prime[i];
            root1 = fb->root1[i];
            root2 = fb->root2[i];
            logp = (uint8_t)fb->logp[i];

            // invalid root (part of poly->a)
            if (prime == 0)
                continue;

            SIEVE_2X;
            SIEVE_1X;
            SIEVE_LAST;
            UPDATE_ROOTS;
        }
    }

    // up to the next multiple of 8
    for (; (i < b3) && (i & 7); i++)
    {
        uint8_t* s2;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = (uint8_t)fb->logp[i];

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    for (; i < vbound; i += 8)
    {
        MED_NEON_SET_LOGP((uint8_t)fb->logp[i]);
        MED_NEON_LOAD;

        // valid lanes: prime != 0 and root1 != 0xffff.  Zero the others'
        // primes and park their roots at 0xffff.
        vinit = vbicq_u16(vtstq_u16(vp, vp), vceqq_u16(vr1, vffff));
        vp = vandq_u16(vp, vinit);
        vr1 = vorrq_u16(vr1, vmvnq_u16(vinit));
        vr2 = vorrq_u16(vr2, vmvnq_u16(vinit));
        vval2 = vinit;

        // until a root2 leaves the block, dump every lane
        do
        {
            MED_NEON_STEP;
            vval2 = vbicq_u16(vval2, vcgtq_u16(vr2, vlast));
        } while (vmaxvq_u16(veorq_u16(vval2, vinit)) == 0);

        // then only the lanes whose root2 is still in the block
        vp = vandq_u16(vp, vval2);
        while (vmaxvq_u16(vval2) != 0)
        {
            if (sink_select)
            {
                MED_NEON_SINK_RMW(vval2, vr1, vr2);
            }
            else
            {
                vst1q_u16(r1, vr1);
                vst1q_u16(r2, vr2);
                bm = MED_NEON_LANEBITS(vval2);
                while (bm)
                {
                    uint32_t k = (uint32_t)__builtin_ctzll(bm) >> 3;
                    sieve[r2[k]] -= MLP(k);
                    sieve[r1[k]] -= MLP(k);
                    bm &= bm - 1;
                }
            }

            vr1 = vaddq_u16(vr1, vp);
            vr2 = vaddq_u16(vr2, vp);
            vval2 = vbicq_u16(vval2, vcgtq_u16(vr2, vlast));
            vp = vandq_u16(vp, vval2);
        }

        // last hit for root1s still in the block
        vval1 = vbicq_u16(vinit, vcgtq_u16(vr1, vlast));
        if (sink_select)
        {
            MED_NEON_SINK_RMW8(vval1, vr1);
        }
        else
        {
            vst1q_u16(r1, vr1);
            bm = MED_NEON_LANEBITS(vval1);
            while (bm)
            {
                uint32_t k = (uint32_t)__builtin_ctzll(bm) >> 3;
                sieve[r1[k]] -= MLP(k);
                bm &= bm - 1;
            }
        }

        // advance those root1s, reduce both roots, store sorted
        vp = vandq_u16(vld1q_u16(fb->prime + i), vval1);
        vr1 = vsubq_u16(vaddq_u16(vr1, vp), vblock);
        vr2 = vsubq_u16(vr2, vblock);
        vst1q_u16(fb->root1 + i, vminq_u16(vr1, vr2));
        vst1q_u16(fb->root2 + i, vmaxq_u16(vr1, vr2));
    }

    for (; i < b3; i++)
    {
        MED_NEON_SCALAR_1X;
    }

    // fixed passes for p > 32768 / 3
    for (; (i + 8) <= med_B; i += 8)
    {
        MED_NEON_SET_LOGP((uint8_t)fb->logp[i]);
        MED_NEON_LOAD;
        // p < 32768 (p = 0 included)
        nbig = vcleq_u16(vp, vlast);

        MED_NEON_FIXED_PASS(vffff);
        if (plist[i] < 32768)
        {
            MED_NEON_FIXED_PASS(nbig);
            if (plist[i] < 16384)
            {
                MED_NEON_FIXED_PASS(nbig);
            }
        }

        vr1 = vsubq_u16(vr1, vblock);
        vr2 = vsubq_u16(vr2, vblock);
        vst1q_u16(fb->root1 + i, vminq_u16(vr1, vr2));
        vst1q_u16(fb->root2 + i, vmaxq_u16(vr1, vr2));
    }

    for (; i < med_B; i++)
    {
        MED_NEON_SCALAR_1X;
    }

    return;
}

// all small primes vector, selective loops over set lanes
void med_sieveblock_32k_neon_split(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    med_neon_split_body(sieve, fb, full_fb, start_prime, 0, 0);
}

// all small primes vector, sink-strip selective loops
void med_sieveblock_32k_neon_splitsink(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    med_neon_split_body(sieve, fb, full_fb, start_prime, 0, 1);
}

// scalar below 2^13, then as split
void med_sieveblock_32k_neon_hybrid(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    med_neon_split_body(sieve, fb, full_fb, start_prime, 1, 0);
}

// scalar below 2^13, then as splitsink
void med_sieveblock_32k_neon_hybridsink(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    med_neon_split_body(sieve, fb, full_fb, start_prime, 1, 1);
}

#undef MED_NEON_SINK_RMW
#undef MED_NEON_SINK_RMW8
#undef MED_NEON_FIXED_PASS
#undef MED_NEON_LANEBITS
#undef MED_NEON_LOAD
#undef MED_NEON_SCALAR_1X
#undef MED_NEON_FINAL_AND_STORE
#undef MED_NEON_LANE
#undef MED_NEON_STEP

#endif
