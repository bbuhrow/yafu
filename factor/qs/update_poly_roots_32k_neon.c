/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// AArch64 NEON per-poly root updates.
//
// nextRoots_32k_neon_small: NEON version of nextRoots_32k_generic_small
// (primes below med_B).  The tiny primes below sieve_small_fb_start stay
// scalar; [sieve_small_fb_start, med_B) is done 8 primes at a time in
// 16-bit lanes, with results identical to the portable code.
//
// nextRoots_32k_lp_buckets_neon: NEON version of nextRoots_32k_lp_buckets
// (the bucket ranges above med_B).  Bit-identical buckets: see below.

#include <stdlib.h>
#include <string.h>
#include "qs_impl.h"
#include "poly_macros_32k.h"
#include "poly_macros_common.h"

#if defined(__aarch64__) && defined(__ARM_NEON)
#include "neon_util.h"

// 8 consecutive 32-bit values (all < 2^16) narrowed to 16-bit lanes
static inline uint16x8_t neon_narrow8_u32(const uint32_t* p)
{
    return vcombine_u16(vmovn_u32(vld1q_u32(p)), vmovn_u32(vld1q_u32(p + 4)));
}

void nextRoots_32k_neon_small(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    sieve_fb_compressed* fb_p = dconf->comp_sieve_p;
    sieve_fb_compressed* fb_n = dconf->comp_sieve_n;
    int* rootupdates = dconf->rootupdates;
    update_t update_data = dconf->update_data;
    uint32_t startprime = 2;
    uint32_t bound = sconf->factor_base->B;
    uint32_t med_B = sconf->factor_base->med_B;
    uint32_t sm_start = sconf->sieve_small_fb_start;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    int* ptr = &rootupdates[(uint32_t)(v - 1) * bound + startprime];
    uint32_t j, root1, root2, prime;

    // tiny primes: not sieved, so root order doesn't matter
    for (j = startprime; j < sm_start; j++, ptr++)
    {
        prime = update_data.prime[j];
        root1 = (uint32_t)update_data.firstroots1[j];
        root2 = (uint32_t)update_data.firstroots2[j];

        if (sign > 0)
        {
            COMPUTE_NEXT_ROOTS_P;
        }
        else
        {
            COMPUTE_NEXT_ROOTS_N;
        }

        update_data.firstroots1[j] = (int)root1;
        update_data.firstroots2[j] = (int)root2;

        fb_p->root1[j] = (uint16_t)root1;
        fb_p->root2[j] = (uint16_t)root2;
        fb_n->root1[j] = (uint16_t)(prime - root2);
        fb_n->root2[j] = (uint16_t)(prime - root1);
        if (fb_n->root1[j] == prime)
            fb_n->root1[j] = 0;
        if (fb_n->root2[j] == prime)
            fb_n->root2[j] = 0;
    }

    // [sieve_small_fb_start, med_B): 8 primes at a time.  sign > 0:
    // r = r - u (+p on borrow).  sign < 0: r = r + u mod p, computed as
    // r - (p - u) (+p on borrow) because r + u can exceed 16 bits.
    for (; (j + 8) <= med_B; j += 8, ptr += 8)
    {
        uint16x8_t vp = neon_narrow8_u32(update_data.prime + j);
        uint16x8_t vu = neon_narrow8_u32((const uint32_t*)ptr);
        uint16x8_t r1 = vld1q_u16(update_data.sm_firstroots1 + j);
        uint16x8_t r2 = vld1q_u16(update_data.sm_firstroots2 + j);
        uint16x8_t lo, hi;

        if (sign <= 0)
            vu = vsubq_u16(vp, vu);

        r1 = vaddq_u16(vsubq_u16(r1, vu), vandq_u16(vp, vcltq_u16(r1, vu)));
        r2 = vaddq_u16(vsubq_u16(r2, vu), vandq_u16(vp, vcltq_u16(r2, vu)));

        lo = vminq_u16(r1, r2);
        hi = vmaxq_u16(r1, r2);
        vst1q_u16(update_data.sm_firstroots1 + j, lo);
        vst1q_u16(update_data.sm_firstroots2 + j, hi);
        vst1q_u16(fb_p->root1 + j, lo);
        vst1q_u16(fb_p->root2 + j, hi);
        vst1q_u16(fb_n->root1 + j, vsubq_u16(vp, hi));
        vst1q_u16(fb_n->root2 + j, vsubq_u16(vp, lo));
    }

    for (; j < med_B; j++, ptr++)
    {
        prime = update_data.prime[j];
        root1 = (uint32_t)update_data.sm_firstroots1[j];
        root2 = (uint32_t)update_data.sm_firstroots2[j];

        if (sign > 0)
        {
            COMPUTE_NEXT_ROOTS_P;
        }
        else
        {
            COMPUTE_NEXT_ROOTS_N;
        }

        if (root2 < root1)
        {
            uint32_t t = root1;
            root1 = root2;
            root2 = t;
        }

        update_data.sm_firstroots1[j] = (uint16_t)root1;
        update_data.sm_firstroots2[j] = (uint16_t)root2;
        fb_p->root1[j] = (uint16_t)root1;
        fb_p->root2[j] = (uint16_t)root2;
        fb_n->root1[j] = (uint16_t)(prime - root2);
        fb_n->root2[j] = (uint16_t)(prime - root1);
    }

    return;
}

// Bucket fill for one root stream of prime j: every hit r, r+p, ... below
// the interval appends (j - bound_val) << 16 | (r & 32767) to bucket r>>15.
// "fixed": K passes with no branches; a miss writes the next free slot of
// bucket 0 and adds 0 to its count, so the buckets end up identical to the
// portable loop's.  K = interval / p_min + 1 bounds the hits of any root
// r < p of the group's primes.
#define NEON_FILL_STREAM_FIXED(r0, slice, num) \
    { \
        uint32_t _r = (r0), _t; \
        for (_t = 0; _t < K; _t++) \
        { \
            uint32_t _b = _r >> 15; \
            uint32_t _m = (_b < (uint32_t)numblocks); \
            _b = _m ? _b : 0; \
            slice[(_b << BUCKET_BITS) + num[_b]] = elem_hi | (_r & 32767); \
            num[_b] += _m; \
            _r += prime; \
        } \
    }

// the portable loop (FILL_ONE_PRIME_LOOP_P/N), one stream
#define NEON_FILL_STREAM_LOOP(r0, slice, num) \
    { \
        uint32_t _r = (r0); \
        uint32_t _b = _r >> 15; \
        while (_b < (uint32_t)numblocks) \
        { \
            slice[(_b << BUCKET_BITS) + num[_b]] = elem_hi | (_r & 32767); \
            num[_b]++; \
            _r += prime; \
            _b = _r >> 15; \
        } \
    }

// The root updates are done 4 primes at a time (32-bit lanes, the same
// formulas as the portable code; the bucket path keeps roots unsorted).
// The fills then go prime by prime in the portable order (root1, root2,
// -root1, -root2), each after its own CHECK_NEW_SLICE, so slice bounds,
// bucket contents and element order all match nextRoots_32k_lp_buckets.
// YAFU_NEON_BUCKETS=loop|fixed picks the fill (default fixed).
void nextRoots_32k_lp_buckets_neon(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    static int fill_mode = -1;   // 0 loop, 1 fixed
    fb_list* fb = sconf->factor_base;
    update_t update_data = dconf->update_data;
    lp_bucket* lp_bucket_p = dconf->buckets;
    uint32_t bound = fb->B;
    uint32_t sparse_B = fb->fb_blk8i;
    uint32_t large_B = fb->large_B;
    uint32_t stop = sconf->lp_direct ? large_B : bound;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    uint32_t* ptr = (uint32_t*)&dconf->rootupdates[(uint32_t)(v - 1) * bound + sparse_B];
    uint32_t* fr1 = (uint32_t*)update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)update_data.firstroots2;
    uint32_t j, q, root1, root2, prime, nroot1, nroot2, interval, K, elem_hi;
    int numblocks = (int)sconf->num_blocks;
    uint32_t* bptr;
    int bnum, room, k;
    int bound_index = 0;
    int check_bound;
    uint32_t bound_val = sparse_B;
    uint32_t* numptr_p, * numptr_n, * sliceptr_p, * sliceptr_n;
    uint8_t* slicelogp_ptr;
    uint32_t* slicebound_ptr;
    uint8_t logp = 0;
    uint32_t pa[4], a1[4], a2[4], n1[4], n2[4];

    if (fill_mode < 0)
    {
        const char* s = getenv("YAFU_NEON_BUCKETS");
        fill_mode = ((s != NULL) && (strcmp(s, "loop") == 0)) ? 0 : 1;
    }

    interval = (uint32_t)numblocks << 15;

    if ((lp_bucket_p->alloc_slices == 0) || (sparse_B >= stop))
    {
        if (lp_bucket_p->alloc_slices != 0)
            lp_bucket_p->num_slices = 0;
        return;
    }

    // bucket ranges [sparse_B, stop)
    lp_bucket_p->fb_bounds[0] = sparse_B;
    sliceptr_p = lp_bucket_p->list;
    sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);
    numptr_p = lp_bucket_p->num;
    numptr_n = lp_bucket_p->num + numblocks;
    for (j = 0; j < lp_bucket_p->list_size; j++)
        numptr_p[j] = 0;
    lp_bucket_p->num_slices = 0;
    slicelogp_ptr = lp_bucket_p->logp;
    slicebound_ptr = lp_bucket_p->fb_bounds;
    check_bound = (int)(sparse_B + BUCKET_ALLOC / 2);

    logp = update_data.logp[sparse_B - 1];
    for (j = sparse_B; (j + 4) <= stop; j += 4, ptr += 4)
    {
        uint32x4_t vp = vld1q_u32(update_data.prime + j);
        uint32x4_t r1 = vld1q_u32(fr1 + j);
        uint32x4_t r2 = vld1q_u32(fr2 + j);
        uint32x4_t u = vld1q_u32(ptr);

        if (sign > 0)
        {
            r1 = vaddq_u32(vsubq_u32(r1, u), vandq_u32(vp, vcltq_u32(r1, u)));
            r2 = vaddq_u32(vsubq_u32(r2, u), vandq_u32(vp, vcltq_u32(r2, u)));
        }
        else
        {
            r1 = vaddq_u32(r1, u);
            r2 = vaddq_u32(r2, u);
            r1 = vsubq_u32(r1, vandq_u32(vp, vcgeq_u32(r1, vp)));
            r2 = vsubq_u32(r2, vandq_u32(vp, vcgeq_u32(r2, vp)));
        }

        vst1q_u32(fr1 + j, r1);
        vst1q_u32(fr2 + j, r2);
        vst1q_u32(pa, vp);
        vst1q_u32(a1, r1);
        vst1q_u32(a2, r2);
        vst1q_u32(n1, vsubq_u32(vp, r1));
        vst1q_u32(n2, vsubq_u32(vp, r2));
        K = interval / pa[0] + 1;

        for (q = 0; q < 4; q++)
        {
            uint32_t jj = j + q;

            CHECK_NEW_SLICE(jj);

            prime = pa[q];
            if (jj < large_B)
            {
                // p < interval: possibly several hits per root
                elem_hi = (jj - bound_val) << 16;
                if (fill_mode)
                {
                    NEON_FILL_STREAM_FIXED(a1[q], sliceptr_p, numptr_p);
                    NEON_FILL_STREAM_FIXED(a2[q], sliceptr_p, numptr_p);
                    NEON_FILL_STREAM_FIXED(n1[q], sliceptr_n, numptr_n);
                    NEON_FILL_STREAM_FIXED(n2[q], sliceptr_n, numptr_n);
                }
                else
                {
                    NEON_FILL_STREAM_LOOP(a1[q], sliceptr_p, numptr_p);
                    NEON_FILL_STREAM_LOOP(a2[q], sliceptr_p, numptr_p);
                    NEON_FILL_STREAM_LOOP(n1[q], sliceptr_n, numptr_n);
                    NEON_FILL_STREAM_LOOP(n2[q], sliceptr_n, numptr_n);
                }
            }
            else
            {
                root1 = a1[q];
                root2 = a2[q];
                FILL_ONE_PRIME_P(jj);
                root1 = n1[q];
                root2 = n2[q];
                FILL_ONE_PRIME_N(jj);
            }
        }
    }

    // tail: the portable code
    for (; j < stop; j++, ptr++)
    {
        uint32_t uu = *ptr;

        CHECK_NEW_SLICE(j);

        prime = update_data.prime[j];
        root1 = fr1[j];
        root2 = fr2[j];

        if (sign > 0)
        {
            root1 = (uu > root1) ? root1 + prime - uu : root1 - uu;
            root2 = (uu > root2) ? root2 + prime - uu : root2 - uu;
        }
        else
        {
            root1 += uu;
            root2 += uu;
            if (root1 >= prime) root1 -= prime;
            if (root2 >= prime) root2 -= prime;
        }

        fr1[j] = root1;
        fr2[j] = root2;
        nroot1 = prime - root1;
        nroot2 = prime - root2;

        if (j < large_B)
        {
            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);
        }
        else
        {
            FILL_ONE_PRIME_P(j);
            root1 = nroot1;
            root2 = nroot2;
            FILL_ONE_PRIME_N(j);
        }
    }

    lp_bucket_p->num_slices = (uint32_t)(bound_index + 1);
    lp_bucket_p->logp[bound_index] = (uint8_t)(sconf->lp_direct ?
        fb->list->logprime[large_B - 1] : logp);

    return;
}

#undef NEON_FILL_STREAM_FIXED
#undef NEON_FILL_STREAM_LOOP

#endif
