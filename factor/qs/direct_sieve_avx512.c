/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// AVX-512 kernels for the large-prime strategies (see direct_sieve.c for the
// scalar versions and direct_sieve_avx2.c / direct_sieve_v4.h for the AVX2
// and 4-lane ones).  These are the reference kernels: every other ISA visits
// hits and divisions in the same order and does bit-identical work.
// Compiled into AVX-512BW builds.
//
// Each kernel is a SIQS_FORCEINLINE body taking a constant kmode, with thin
// per-kernel wrappers (_bitloop, _sink, _creg, _cmem) so the hot loops have
// no branches on the kernel choice.

#include <stdint.h>
#include <string.h>
#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "tdiv_macros_common.h"

#if defined(TRY_DIRECT_SIEVE) && defined(USE_AVX512BW)
#include <immintrin.h>

// ----------------------------------------------------------------------
// helpers: compress stores and the 16-prime scatter
// ----------------------------------------------------------------------

//#define ENABLE_AMDEYPC_FASTPATH 1

// append the lanes of v selected by m to buf[n...]; returns the new count.
// kmode (a compile-time constant at every call site) picks the form:
// SIQS_K_CREG compresses in a register and stores the whole vector (buf
// must have 16 entries of room past n), anything else uses compressstoreu.
SIQS_FORCEINLINE uint32_t direct_cstore16(uint32_t* buf, uint32_t n, __mmask16 m, __m512i v,
    const int kmode)
{
    if (kmode == SIQS_K_CREG)
        _mm512_storeu_si512((__m512i*)(buf + n), _mm512_maskz_compress_epi32(m, v));
    else
        _mm512_mask_compressstoreu_epi32(buf + n, m, v);
    return n + _mm_popcnt_u32(m);
}

// append the in-range lanes of 16 primes' 4 root vectors to buf[n...] as
// (logp << 24) | sieve position (pos side r, neg side interval + p - r;
// 2*interval < 2^24).  Returns the new count.
SIQS_FORCEINLINE uint32_t direct_compress16(uint32_t* buf, uint32_t n, uint32_t interval,
    __m512i vroot1, __m512i vroot2, __m512i vnroot1, __m512i vnroot2, uint8_t logp,
    const int kmode)
{
    __m512i vint = _mm512_set1_epi32(interval);
    __m512i vtag = _mm512_set1_epi32((uint32_t)logp << 24);
    __mmask16 m1 = _mm512_cmplt_epu32_mask(vroot1, vint);
    __mmask16 m2 = _mm512_cmplt_epu32_mask(vroot2, vint);
    __mmask16 m3 = _mm512_cmplt_epu32_mask(vnroot1, vint);
    __mmask16 m4 = _mm512_cmplt_epu32_mask(vnroot2, vint);
    __m512i v1 = _mm512_or_epi32(vroot1, vtag);
    __m512i v2 = _mm512_or_epi32(vroot2, vtag);
    __m512i v3 = _mm512_or_epi32(_mm512_add_epi32(vnroot1, vint), vtag);
    __m512i v4 = _mm512_or_epi32(_mm512_add_epi32(vnroot2, vint), vtag);

    n = direct_cstore16(buf, n, m1, v1, kmode);
    n = direct_cstore16(buf, n, m2, v2, kmode);
    n = direct_cstore16(buf, n, m3, v3, kmode);
    return direct_cstore16(buf, n, m4, v4, kmode);
}

// sink-stripe scatter of 16 primes' 4 root vectors into both intervals
// (pos side sieve[0, interval), neg side sieve[interval, 2*interval)).
// vroot = pos-side roots, vnroot = p - root (neg-side offsets).  Every lane
// is written: lanes outside the interval go to the stripe at
// sieve[2*interval + 32768 + *sink ...], which is allocated (second tail
// guard block) and never scanned.  *sink rotates through the stripe so
// consecutive miss RMWs don't chain through store forwarding.
SIQS_FORCEINLINE void direct_sink_scatter16(uint8_t* sieve, uint32_t interval,
    __m512i vroot1, __m512i vroot2, __m512i vnroot1, __m512i vnroot2,
    uint8_t logp, uint32_t* sink)
{
    ALIGNED_MEM uint32_t idx[64];
    const __m512i vlane = _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15);
    __m512i vint = _mm512_set1_epi32(interval);
    __m512i vs = _mm512_add_epi32(_mm512_set1_epi32(2 * interval + 32768 + *sink), vlane);
    __m512i v16 = _mm512_set1_epi32(16);
    __mmask16 m1 = _mm512_cmplt_epu32_mask(vroot1, vint);
    __mmask16 m2 = _mm512_cmplt_epu32_mask(vroot2, vint);
    __mmask16 m3 = _mm512_cmplt_epu32_mask(vnroot1, vint);
    __mmask16 m4 = _mm512_cmplt_epu32_mask(vnroot2, vint);

    _mm512_store_si512(idx, _mm512_mask_blend_epi32(m1, vs, vroot1));
    vs = _mm512_add_epi32(vs, v16);
    _mm512_store_si512(idx + 16, _mm512_mask_blend_epi32(m2, vs, vroot2));
    vs = _mm512_add_epi32(vs, v16);
    _mm512_store_si512(idx + 32, _mm512_mask_blend_epi32(m3, vs,
        _mm512_add_epi32(vnroot1, vint)));
    vs = _mm512_add_epi32(vs, v16);
    _mm512_store_si512(idx + 48, _mm512_mask_blend_epi32(m4, vs,
        _mm512_add_epi32(vnroot2, vint)));

#define DIRECT_SINK_RMW4(k) \
    sieve[idx[k]] -= logp; sieve[idx[k + 1]] -= logp; \
    sieve[idx[k + 2]] -= logp; sieve[idx[k + 3]] -= logp;
#define DIRECT_SINK_RMW16(k) \
    DIRECT_SINK_RMW4(k) DIRECT_SINK_RMW4(k + 4) \
    DIRECT_SINK_RMW4(k + 8) DIRECT_SINK_RMW4(k + 12)

    DIRECT_SINK_RMW16(0);
    DIRECT_SINK_RMW16(16);
    DIRECT_SINK_RMW16(32);
    DIRECT_SINK_RMW16(48);

#undef DIRECT_SINK_RMW16
#undef DIRECT_SINK_RMW4

    *sink += 64;
    if (*sink >= DIRECT_SINK_BYTES)
        *sink = 0;
}

// scatter one 16-prime group of [large_B, B) into both intervals with
// kernel kmode (a compile-time constant at each call site): pos hits at r,
// neg hits at p - r in sieve[interval, 2*interval).  BITLOOP and SINK write
// at once; CREG / CMEM append to cbuf (count *nc) and sieve it whenever it
// reaches DIRECT_CBUF_FLUSH, so callers must direct_compress_flush the rest
// after their loop.  *sink is the sink-stripe rotation (SINK only).
SIQS_FORCEINLINE void direct_scatter16(uint8_t* sieve, uint32_t interval,
    __m512i vprime, __m512i vroot1, __m512i vroot2, uint8_t logp, const int kmode,
    uint32_t* cbuf, uint32_t* nc, uint32_t* sink)
{
    __m512i vnroot1 = _mm512_sub_epi32(vprime, vroot1);
    __m512i vnroot2 = _mm512_sub_epi32(vprime, vroot2);

#ifdef ENABLE_AMDEYPC_FASTPATH
    __m512i vint = _mm512_set1_epi32(interval);
    //__m512i vtag = _mm512_set1_epi32((uint32_t)logp << 24);
    __mmask16 m1 = _mm512_cmplt_epu32_mask(vroot1, vint);
    __mmask16 m2 = _mm512_cmplt_epu32_mask(vroot2, vint);
    __mmask16 m3 = _mm512_cmplt_epu32_mask(vnroot1, vint);
    __mmask16 m4 = _mm512_cmplt_epu32_mask(vnroot2, vint);
    //__m512i v1 = _mm512_or_epi32(vroot1, vtag);
    //__m512i v2 = _mm512_or_epi32(vroot2, vtag);
    //__m512i v3 = _mm512_or_epi32(_mm512_add_epi32(vnroot1, vint), vtag);
    //__m512i v4 = _mm512_or_epi32(_mm512_add_epi32(vnroot2, vint), vtag);
    __m512i v1 = vroot1;
    __m512i v2 = vroot2;
    __m512i v3 = _mm512_add_epi32(vnroot1, vint);
    __m512i v4 = _mm512_add_epi32(vnroot2, vint);

    v1 = _mm512_maskz_compress_epi32(m1, v1);
    v2 = _mm512_maskz_compress_epi32(m2, v2);
    v3 = _mm512_maskz_compress_epi32(m3, v3);
    v4 = _mm512_maskz_compress_epi32(m4, v4);
    int idx1 = _mm_popcnt_u32(m1);
    int idx2 = _mm_popcnt_u32(m2);
    int idx3 = _mm_popcnt_u32(m3);
    int idx4 = _mm_popcnt_u32(m4);
    _mm512_storeu_si512((__m512i*)(cbuf + *nc), v1);
    _mm512_storeu_si512((__m512i*)(cbuf + *nc + idx1), v2);
    _mm512_storeu_si512((__m512i*)(cbuf + *nc + idx1 + idx2), v3);
    _mm512_storeu_si512((__m512i*)(cbuf + *nc + idx1 + idx2 + idx3), v4);
    *nc = *nc + idx1 + idx2 + idx3 + idx4;

    if (*nc >= DIRECT_CBUF_FLUSH)
    {
        direct_compress_flush_fast(sieve, cbuf, *nc, logp);
        *nc = 0;
    }

#else

    if ((kmode == SIQS_K_CREG) || (kmode == SIQS_K_CMEM))
    {
        *nc = direct_compress16(cbuf, *nc, interval, vroot1, vroot2,
            vnroot1, vnroot2, logp, kmode);
        if (*nc >= DIRECT_CBUF_FLUSH)
        {
            direct_compress_flush(sieve, cbuf, *nc);
            *nc = 0;
        }
    }
    else if (kmode == SIQS_K_SINK)
    {
        direct_sink_scatter16(sieve, interval, vroot1, vroot2, vnroot1, vnroot2,
            logp, sink);
    }
    else
    {
        ALIGNED_MEM uint32_t r_1[16];
        ALIGNED_MEM uint32_t r_2[16];
        ALIGNED_MEM uint32_t n_1[16];
        ALIGNED_MEM uint32_t n_2[16];
        uint8_t* sieven = sieve + interval;
        __m512i vint = _mm512_set1_epi32(interval);
        uint32_t m1 = _mm512_cmplt_epu32_mask(vroot1, vint);
        uint32_t m2 = _mm512_cmplt_epu32_mask(vroot2, vint);
        uint32_t m3 = _mm512_cmplt_epu32_mask(vnroot1, vint);
        uint32_t m4 = _mm512_cmplt_epu32_mask(vnroot2, vint);
        int pos;

        if ((m1 | m2 | m3 | m4) == 0)
            return;

        _mm512_store_si512(r_1, vroot1);
        _mm512_store_si512(r_2, vroot2);
        _mm512_store_si512(n_1, vnroot1);
        _mm512_store_si512(n_2, vnroot2);

        while (m1 > 0) {
            pos = _trail_zcnt(m1);
            sieve[r_1[pos]] -= logp;
            m1 = _reset_lsb(m1);
        }
        while (m2 > 0) {
            pos = _trail_zcnt(m2);
            sieve[r_2[pos]] -= logp;
            m2 = _reset_lsb(m2);
        }
        while (m3 > 0) {
            pos = _trail_zcnt(m3);
            sieven[n_1[pos]] -= logp;
            m3 = _reset_lsb(m3);
        }
        while (m4 > 0) {
            pos = _trail_zcnt(m4);
            sieven[n_2[pos]] -= logp;
            m4 = _reset_lsb(m4);
        }
    }

#endif
}

// ----------------------------------------------------------------------
// direct sieve pass 1: [large_B, B) into both whole intervals
// ----------------------------------------------------------------------
// pass 1 of the direct interval sieve (see direct_sieve.md): sieve every
// prime in [large_B, B) into both whole intervals, pos side at sieve[0,
// interval) and neg side at sieve[interval, 2*interval).  All these primes
// are > interval (large_B rounds up), so each root hits each side at most
// once.  Roots come straight from update_data.firstroots1/2, treated as an
// unordered pair: pos hit at r, neg hit at p - r (r == 0 gives p, which is
// out of range).  Roots are not modified here; pass 2
// (direct_resieve_update_body) re-derives the same hits and
// advances the roots.
SIQS_FORCEINLINE void direct_sieve_largeprimes_body(uint8_t* sieve, static_conf_t* sconf,
    dynamic_conf_t* dconf, const int kmode)
{
    fb_list* fb = sconf->factor_base;
    update_t update_data = dconf->update_data;
    uint32_t interval = sconf->num_blocks << 15;
    uint8_t* sieven = sieve + interval;
    uint32_t* primes = update_data.prime;
    uint32_t* fr1 = (uint32_t*)update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)update_data.firstroots2;
    uint32_t j;
    uint8_t logp;
    ALIGNED_MEM uint32_t cbuf[DIRECT_CBUF_SIZE];
    uint32_t nc = 0, sink = 0;

    for (j = fb->large_B; (j + 16) <= fb->B; j += 16)
    {
        // approximate the next 16 logp's as equal to the largest
        logp = update_data.logp[j + 15] + DIRECT_LOGP_BIAS;
        direct_scatter16(sieve, interval, _mm512_loadu_si512(primes + j),
            _mm512_loadu_si512(fr1 + j), _mm512_loadu_si512(fr2 + j), logp, kmode,
            cbuf, &nc, &sink);
    }
#ifdef ENABLE_AMDEYPC_FASTPATH
    direct_compress_flush_fast(sieve, cbuf, nc, update_data.logp[fb->B - 1]);
#else
    direct_compress_flush(sieve, cbuf, nc);
#endif
    // scalar tail (empty when B - large_B is a multiple of 16)
    for (; j < fb->B; j++)
    {
        uint32_t prime = primes[j];
        logp = update_data.logp[j] + DIRECT_LOGP_BIAS;

        if (fr1[j] < interval) sieve[fr1[j]] -= logp;
        if (fr2[j] < interval) sieve[fr2[j]] -= logp;
        if ((prime - fr1[j]) < interval) sieven[prime - fr1[j]] -= logp;
        if ((prime - fr2[j]) < interval) sieven[prime - fr2[j]] -= logp;
    }

    return;
}

void direct_sieve_largeprimes_avx512_bitloop(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_body(sieve, sconf, dconf, SIQS_K_BITLOOP);
}

void direct_sieve_largeprimes_avx512_sink(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_body(sieve, sconf, dconf, SIQS_K_SINK);
}

void direct_sieve_largeprimes_avx512_creg(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_body(sieve, sconf, dconf, SIQS_K_CREG);
}

void direct_sieve_largeprimes_avx512_cmem(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_body(sieve, sconf, dconf, SIQS_K_CMEM);
}

// ----------------------------------------------------------------------
// direct sieve fused pass: resieve, root update, next-poly scatter
// ----------------------------------------------------------------------
// (the survivor pool, direct_resieve_hit and direct_finalize_survivors are
// in direct_sieve.c)

// fused pass of the direct interval sieve: for every prime in [large_B, B),
// resieve both sides with the CURRENT poly's roots (the same values that
// were sieved, so exactly the sieve's hit set) against the survivor pool,
// then advance the roots to the next poly exactly as nextRoots would after
// nextB, and sieve the next poly's roots into 'sieve' (both whole
// intervals, already initialized by the caller).  Must run before nextB:
// nextB and nextRoots both read the gray code entry for dconf->numB, which
// is what we use here.  With sieve == NULL (last poly of this 'a') only the
// resieve is done; the roots are left as they are.
// kmode (constant in each wrapper below) is the scatter kernel for the
// next poly's hits (see direct_scatter16).
SIQS_FORCEINLINE void direct_resieve_update_body(static_conf_t* sconf, dynamic_conf_t* dconf,
    uint8_t* sieve, const int kmode)
{
    fb_list* fb = sconf->factor_base;
    direct_pool_t* pool = dconf->dpool;
    update_t update_data = dconf->update_data;
    uint32_t interval = pool->interval;
    uint32_t large_B = fb->large_B;
    uint32_t bound = fb->B;
    uint32_t* primes = update_data.prime;
    uint32_t* fr1 = (uint32_t*)update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)update_data.firstroots2;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    uint32_t* ptr = (uint32_t*)&dconf->rootupdates[(v - 1) * bound + large_B];
    int resieve = (pool->active > 0);
    uint8_t* sieven = (sieve == NULL) ? NULL : sieve + interval;
    uint32_t j;
    uint8_t logp;

    ALIGNED_MEM uint32_t r_1[16];
    ALIGNED_MEM uint32_t r_2[16];
    ALIGNED_MEM uint32_t n_1[16];
    ALIGNED_MEM uint32_t n_2[16];
    ALIGNED_MEM uint32_t p_16[16];

    __m512i vint = _mm512_set1_epi32(interval);
    __m512i vint2 = _mm512_set1_epi32(2 * interval);
    __m512i vzero = _mm512_setzero_si512();
    __m512i vone = _mm512_set1_epi32(1);
    __m512i v31 = _mm512_set1_epi32(31);
    // the bitmap as 32-bit words (little endian: bit pos is bit pos&31 of
    // word pos>>5, same as bit pos&63 of 64-bit word pos>>6)
    const int* bitmap32 = (const int*)pool->bitmap;
    // survivor filter: direct compares against the survivor positions when
    // there are few survivors (cost ~2 compares per survivor per 16 primes),
    // bitmap gathers otherwise (4 gathers per 16 primes regardless).
    int use_cmp = resieve && (pool->num <= DIRECT_CMP_MAX_SURVIVORS);
    uint32_t nsp = 0;   // pos-side survivors: pool->pos[0, nsp) < interval
    ALIGNED_MEM uint32_t cbuf[DIRECT_CBUF_SIZE];
    uint32_t nc = 0, sink = 0;

    if ((sieve == NULL) && (!resieve))
        return;
#if defined(DIRECT_EXIT_STATS)
    if (resieve)
        pool->ex_groups += (bound - large_B) / 16;
#endif

    if (use_cmp)
    {
        while ((nsp < pool->num) && (pool->pos[nsp] < interval))
            nsp++;
    }

#if DIRECT_PREFILTER_BITS > 0
    // bucket prefilter, one bit per 2^pf_shift positions of both intervals,
    // held in one (512 bits) or two (1024 bits) zmm registers
    ALIGNED_MEM uint32_t pf_words[32];
    __m512i vpf_lo = vzero, vpf_hi = vzero, vpf_shift;
    {
        uint32_t k, pf_shift = 0;
        while (((uint64_t)DIRECT_PREFILTER_BITS << pf_shift) < (2 * (uint64_t)interval))
            pf_shift++;
        memset(pf_words, 0, sizeof(pf_words));
        for (k = 0; k < pool->num; k++)
        {
            uint32_t b = pool->pos[k] >> pf_shift;
            pf_words[b >> 5] |= (1u << (b & 31));
        }
        vpf_lo = _mm512_load_si512(pf_words);
        vpf_hi = _mm512_load_si512(pf_words + 16);
        vpf_shift = _mm512_set1_epi32(pf_shift);
    }
#endif

    for (j = large_B; (j + 16) <= bound; j += 16, ptr += 16)
    {
        __m512i vprime = _mm512_loadu_si512(primes + j);
        __m512i vroot1 = _mm512_loadu_si512(fr1 + j);
        __m512i vroot2 = _mm512_loadu_si512(fr2 + j);
        __m512i vpval = _mm512_loadu_si512(ptr);
        __mmask16 mask1, mask2;

        if (resieve)
        {
            // interval positions of all 4 roots (neg side offset by interval)
            __m512i vnpos1 = _mm512_add_epi32(_mm512_sub_epi32(vprime, vroot1), vint);
            __m512i vnpos2 = _mm512_add_epi32(_mm512_sub_epi32(vprime, vroot2), vint);
            __mmask16 m1, m2, m3, m4;
            int pos;

            // most lanes miss the interval or land far from any survivor, so
            // find the lanes whose position is a survivor's before taking the
            // scalar path.  First the in-range masks...
            m1 = _mm512_cmplt_epu32_mask(vroot1, vint);
            m2 = _mm512_cmplt_epu32_mask(vroot2, vint);
            m3 = _mm512_cmplt_epu32_mask(vnpos1, vint2);
            m4 = _mm512_cmplt_epu32_mask(vnpos2, vint2);

#if DIRECT_PREFILTER_BITS > 0
            // ...then the register-resident bucket prefilter: keep only lanes
            // whose (2^pf_shift)-position bucket holds a survivor.  Lanes in
            // range always have bucket < DIRECT_PREFILTER_BITS.
#if DIRECT_PREFILTER_BITS == 1024
#define DIRECT_PF_WORD(vb) _mm512_permutex2var_epi32(vpf_lo, _mm512_srli_epi32(vb, 5), vpf_hi)
#else
#define DIRECT_PF_WORD(vb) _mm512_permutexvar_epi32(_mm512_srli_epi32(vb, 5), vpf_lo)
#endif
#define DIRECT_PREFILTER(m, vpos) \
            do { \
                __m512i vb = _mm512_srlv_epi32(vpos, vpf_shift); \
                m = _mm512_mask_test_epi32_mask(m, DIRECT_PF_WORD(vb), \
                    _mm512_sllv_epi32(vone, _mm512_and_epi32(vb, v31))); \
            } while (0)

            DIRECT_PREFILTER(m1, vroot1);
            DIRECT_PREFILTER(m2, vroot2);
            DIRECT_PREFILTER(m3, vnpos1);
            DIRECT_PREFILTER(m4, vnpos2);
#undef DIRECT_PREFILTER
#undef DIRECT_PF_WORD
#endif

#if 0
            if (use_cmp)
            {
                // few survivors: every prime here is > interval, so a root
                // hits survivor position s exactly when it equals s.  Compare
                // pos-side roots with pos-side survivors and neg-side
                // positions with neg-side survivors, skipping a side when no
                // lane is left by the masks above.
                uint32_t k;
                __mmask16 c1 = 0, c2 = 0, c3 = 0, c4 = 0;
                if (m1 | m2)
                {
                    for (k = 0; k < nsp; k++)
                    {
                        __m512i vs = _mm512_set1_epi32(pool->pos[k]);
                        c1 |= _mm512_mask_cmpeq_epu32_mask(m1, vroot1, vs);
                        c2 |= _mm512_mask_cmpeq_epu32_mask(m2, vroot2, vs);
                    }
                }
                if (m3 | m4)
                {
                    for (k = nsp; k < pool->num; k++)
                    {
                        __m512i vs = _mm512_set1_epi32(pool->pos[k]);
                        c3 |= _mm512_mask_cmpeq_epu32_mask(m3, vnpos1, vs);
                        c4 |= _mm512_mask_cmpeq_epu32_mask(m4, vnpos2, vs);
                    }
                }
                m1 = c1; m2 = c2; m3 = c3; m4 = c4;
            }
            else
#endif
            {
                // many survivors: test the survivor bitmap for the remaining
                // lanes (gather the 32-bit word holding each position's bit),
                // skipping the gather when no lane is left
#define DIRECT_BITMAP_TEST(m, vpos) \
                if (m) m = _mm512_mask_test_epi32_mask(m, \
                    _mm512_mask_i32gather_epi32(vzero, m, _mm512_srli_epi32(vpos, 5), bitmap32, 4), \
                    _mm512_sllv_epi32(vone, _mm512_and_epi32(vpos, v31)));

                DIRECT_BITMAP_TEST(m1, vroot1);
                DIRECT_BITMAP_TEST(m2, vroot2);
                DIRECT_BITMAP_TEST(m3, vnpos1);
                DIRECT_BITMAP_TEST(m4, vnpos2);
#undef DIRECT_BITMAP_TEST
            }

            if (m1 | m2 | m3 | m4)
            {
                _mm512_store_si512(r_1, vroot1);
                _mm512_store_si512(r_2, vroot2);
                _mm512_store_si512(n_1, vnpos1);
                _mm512_store_si512(n_2, vnpos2);
                _mm512_store_si512(p_16, vprime);

                while (m1 > 0) {
                    pos = _trail_zcnt(m1);
                    direct_resieve_hit(pool, r_1[pos], j + pos, p_16[pos]);
                    m1 = _reset_lsb(m1);
                }
                while (m2 > 0) {
                    pos = _trail_zcnt(m2);
                    direct_resieve_hit(pool, r_2[pos], j + pos, p_16[pos]);
                    m2 = _reset_lsb(m2);
                }
                while (m3 > 0) {
                    pos = _trail_zcnt(m3);
                    direct_resieve_hit(pool, n_1[pos], j + pos, p_16[pos]);
                    m3 = _reset_lsb(m3);
                }
                while (m4 > 0) {
                    pos = _trail_zcnt(m4);
                    direct_resieve_hit(pool, n_2[pos], j + pos, p_16[pos]);
                    m4 = _reset_lsb(m4);
                }
            }

            if (pool->active == 0)
            {
                //printf("resieve now inactive\n");
                resieve = 0;
#if defined(DIRECT_EXIT_STATS)
                // 16-groups of [large_B, B) left unresieved
                pool->ex_skip += (bound - j) / 16;
#endif
            }
        }

        if (sieve == NULL)
            continue;

        // advance to the next poly (same arithmetic as the L2 sections of
        // nextRoots_32k_avx2_small); the pair stays unordered
        if (sign > 0)
        {
            mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
            mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
            vroot1 = _mm512_sub_epi32(vroot1, vpval);
            vroot2 = _mm512_sub_epi32(vroot2, vpval);
            vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
        }
        else
        {
            vroot1 = _mm512_add_epi32(vroot1, vpval);
            vroot2 = _mm512_add_epi32(vroot2, vpval);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
        }
        _mm512_storeu_si512(fr1 + j, vroot1);
        _mm512_storeu_si512(fr2 + j, vroot2);

        // sieve the next poly's roots, both sides (as pass 1)
        // approximate the next 16 logp's as equal to the largest
        logp = update_data.logp[j + 15] + DIRECT_LOGP_BIAS;
        direct_scatter16(sieve, interval, vprime, vroot1, vroot2, logp, kmode,
            cbuf, &nc, &sink);
    }

#ifdef ENABLE_AMDEYPC_FASTPATH
    direct_compress_flush_fast(sieve, cbuf, nc, update_data.logp[bound - 1]);
#else
    // nc stays 0 when sieve == NULL
    direct_compress_flush(sieve, cbuf, nc);
#endif

    // scalar tail (empty when B - large_B is a multiple of 16)
    for (; j < bound; j++, ptr++)
    {
        uint32_t prime = primes[j];
        uint32_t root1 = fr1[j];
        uint32_t root2 = fr2[j];
        uint32_t u = *ptr;

        if (resieve)
        {
            if (root1 < interval) direct_resieve_hit(pool, root1, j, prime);
            if (root2 < interval) direct_resieve_hit(pool, root2, j, prime);
            if ((prime - root1) < interval)
                direct_resieve_hit(pool, interval + prime - root1, j, prime);
            if ((prime - root2) < interval)
                direct_resieve_hit(pool, interval + prime - root2, j, prime);
        }

        if (sieve == NULL)
            continue;

        if (sign > 0)
        {
            root1 = (u > root1) ? root1 + prime - u : root1 - u;
            root2 = (u > root2) ? root2 + prime - u : root2 - u;
        }
        else
        {
            root1 += u;
            root2 += u;
            if (root1 >= prime) root1 -= prime;
            if (root2 >= prime) root2 -= prime;
        }
        fr1[j] = root1;
        fr2[j] = root2;

        logp = update_data.logp[j] + DIRECT_LOGP_BIAS;
        if (root1 < interval) sieve[root1] -= logp;
        if (root2 < interval) sieve[root2] -= logp;
        if ((prime - root1) < interval) sieven[prime - root1] -= logp;
        if ((prime - root2) < interval) sieven[prime - root2] -= logp;
    }

    return;
}

void direct_resieve_update_avx512_bitloop(static_conf_t* sconf, dynamic_conf_t* dconf,
    uint8_t* sieve)
{
    direct_resieve_update_body(sconf, dconf, sieve, SIQS_K_BITLOOP);
}

void direct_resieve_update_avx512_sink(static_conf_t* sconf, dynamic_conf_t* dconf,
    uint8_t* sieve)
{
    direct_resieve_update_body(sconf, dconf, sieve, SIQS_K_SINK);
}

void direct_resieve_update_avx512_creg(static_conf_t* sconf, dynamic_conf_t* dconf,
    uint8_t* sieve)
{
    direct_resieve_update_body(sconf, dconf, sieve, SIQS_K_CREG);
}

void direct_resieve_update_avx512_cmem(static_conf_t* sconf, dynamic_conf_t* dconf,
    uint8_t* sieve)
{
    direct_resieve_update_body(sconf, dconf, sieve, SIQS_K_CMEM);
}

// ----------------------------------------------------------------------
// sparse sieve of one block: [med_B, fb_blk8i)
// ----------------------------------------------------------------------
// sparse sieve of one 32k block for fb [med_B, fb_blk8i): every prime is
// > 32768 and root1 <= root2, so each root hits the block at most once.
// kmode (constant in each wrapper below) picks the hit kernel:
// SIQS_K_BITLOOP loops over the in-range mask bits; SIQS_K_CREG / CMEM
// append the hits as (logp << 24) | block offset with compress stores and
// sieve them in batches.  All kernels do identical sieve updates.
SIQS_FORCEINLINE void large_sieveblock_body(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    const int kmode)
{
    uint32_t i;
    uint32_t prime, root1, root2;
    uint8_t logp;
    __m512i vblock = _mm512_set1_epi32(32768);
    uint32_t bound = full_fb->fb_blk8i;

    if (kmode == SIQS_K_BITLOOP)
    {
        for (i = full_fb->med_B; (i + 16) <= bound; i += 16)
        {
            __m512i vprime, vroot1, vroot2;
            __mmask16 valid_mask;
            uint32_t msk_2;
            int pos;

            vprime = _mm512_load_si512(fb->prime + i);
            vroot1 = _mm512_load_si512(fb->root1 + i);
            vroot2 = _mm512_load_si512(fb->root2 + i);
            // approximate the next 16 logp's as equal to the largest
            logp = fb->logp[i + 15];

            // roots will be increasingly less likely to land in a block
            // as their stride increases, but we should get several to
            // hit of the 16 root2's (root1's will follow - they are ordered).
            valid_mask = _mm512_cmplt_epu32_mask(vroot2, vblock);
            msk_2 = valid_mask;

            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[fb->root1[i + pos]] -= logp;
                sieve[fb->root2[i + pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }

            vroot1 = _mm512_mask_add_epi32(vroot1, valid_mask, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, valid_mask, vroot2, vprime);

            // now all larger roots are invalid.  Last iteration for
            // possibly still valid root1s.  If they are still valid,
            // record the sieve hit, advance them, and swap with the
            // other root
            valid_mask = _mm512_cmplt_epu32_mask(vroot1, vblock);
            msk_2 = valid_mask;

            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[fb->root1[i + pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }

            // reduce both roots and store back for the next block
            vroot1 = _mm512_mask_add_epi32(vroot1, valid_mask, vroot1, vprime);
            vroot1 = _mm512_sub_epi32(vroot1, vblock);
            vroot2 = _mm512_sub_epi32(vroot2, vblock);
            _mm512_store_si512(fb->root1 + i, _mm512_min_epu32(vroot1, vroot2));
            _mm512_store_si512(fb->root2 + i, _mm512_max_epu32(vroot1, vroot2));
        }
    }
    else
    {
        // root1 hits iff root1 < 32768, root2 likewise
        ALIGNED_MEM uint32_t cbuf[DIRECT_CBUF_SIZE];
        uint32_t nc = 0;

        for (i = full_fb->med_B; (i + 16) <= bound; i += 16)
        {
            __m512i vprime = _mm512_load_si512(fb->prime + i);
            __m512i vroot1 = _mm512_load_si512(fb->root1 + i);
            __m512i vroot2 = _mm512_load_si512(fb->root2 + i);
            __m512i vtag;
            __mmask16 m1 = _mm512_cmplt_epu32_mask(vroot1, vblock);
            __mmask16 m2 = _mm512_cmplt_epu32_mask(vroot2, vblock);

            // approximate the next 16 logp's as equal to the largest
            vtag = _mm512_set1_epi32((uint32_t)fb->logp[i + 15] << 24);
            nc = direct_cstore16(cbuf, nc, m1, _mm512_or_epi32(vroot1, vtag), kmode);
            nc = direct_cstore16(cbuf, nc, m2, _mm512_or_epi32(vroot2, vtag), kmode);

            // advance the roots that hit, reduce both for the next block
            vroot1 = _mm512_mask_add_epi32(vroot1, m1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, m2, vroot2, vprime);
            vroot1 = _mm512_sub_epi32(vroot1, vblock);
            vroot2 = _mm512_sub_epi32(vroot2, vblock);
            _mm512_store_si512(fb->root1 + i, _mm512_min_epu32(vroot1, vroot2));
            _mm512_store_si512(fb->root2 + i, _mm512_max_epu32(vroot1, vroot2));

            if (nc >= DIRECT_CBUF_FLUSH)
            {
                direct_compress_flush(sieve, cbuf, nc);
                nc = 0;
            }
        }
        direct_compress_flush(sieve, cbuf, nc);
    }

    // any primes left over after the vector loop (fewer than 16 of them;
    // zero when med_B and fb_blk8i are both 16-aligned) must still be sieved
    // and have their roots advanced.  If they are skipped, their roots fall
    // behind the block and the sparse tdiv matches them against the wrong
    // residue class, causing Q to be divided by non-divisors.  Same
    // per-block semantics as the vector loop: hits where the root is
    // < 32768, then both roots move to the next block in min/max order.
    for (; i < bound; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        while (root2 < (uint32_t)(32768))
        {
            sieve[root1] -= logp;
            sieve[root2] -= logp;
            root1 += prime;
            root2 += prime;
        }

        // root1 can still have one last hit inside the interval
        if (root1 < (uint32_t)(32768))
        {
            sieve[root1] -= logp;
            root1 += prime;
        }

        // wrap both roots into the next interval, store in min/max order
        root1 -= 32768;
        root2 -= 32768;
        fb->root1[i] = (root1 < root2) ? root1 : root2;
        fb->root2[i] = (root1 < root2) ? root2 : root1;
    }

    return;
}

void large_sieveblock_avx512_bitloop(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    large_sieveblock_body(sieve, fb, full_fb, SIQS_K_BITLOOP);
}

void large_sieveblock_avx512_creg(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    large_sieveblock_body(sieve, fb, full_fb, SIQS_K_CREG);
}

void large_sieveblock_avx512_cmem(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    large_sieveblock_body(sieve, fb, full_fb, SIQS_K_CMEM);
}

// ----------------------------------------------------------------------
// sparse-sieve tdiv: [med_B, fb_blk8i)
// ----------------------------------------------------------------------
// the sparse-sieve tdiv for [med_B, fb_blk8i).  kmode (constant in each
// wrapper below) picks how matching lanes are visited: SIQS_K_BITLOOP loops
// over the mask bits, CREG / CMEM compress the matching fb indices into a
// buffer first.  Matches are divided in ascending fb index order either
// way, so fb_offsets come out identical.
SIQS_FORCEINLINE void tdiv_large_body(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf, const int kmode)
{
    // we have flagged this sieve offset as likely to produce a relation
    // nothing left to do now but check and see.
    uint32_t i;
    uint32_t tmp, prime, root1, root2, report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb* fbc;
    uint32_t block_loc;
    uint32_t bound = sconf->factor_base->fb_blk8i;

    if (parity)
    {
        fbc = dconf->sieve_n;
    }
    else
    {
        fbc = dconf->sieve_p;
    }

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {
        if (!dconf->valid_Qs[report_num])
            continue;

        // pull the details of this report to get started.
        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        __m512i vblockloc = _mm512_set1_epi32(block_loc - 32768);

        if (kmode == SIQS_K_BITLOOP)
        {
            for (i = sconf->factor_base->med_B; (i + 16) <= bound; i += 16)
            {
                __m512i vprime, vroot1, vroot2;
                uint32_t msk_2;
                int pos;

                vprime = _mm512_load_si512(fbc->prime + i);
                vroot1 = _mm512_load_si512(fbc->root1 + i);
                vroot2 = _mm512_load_si512(fbc->root2 + i);

                vprime = _mm512_add_epi32(vprime, vblockloc);

                msk_2 = _mm512_cmpeq_epu32_mask(vprime, vroot1) |
                    _mm512_cmpeq_epu32_mask(vprime, vroot2);

                while (msk_2 > 0) {
                    pos = _trail_zcnt(msk_2);
                    DIVIDE_RESIEVED_PRIME(pos);
                    msk_2 = _reset_lsb(msk_2);
                }
            }
        }
        else
        {
            __m512i vidx = _mm512_add_epi32(_mm512_set1_epi32(sconf->factor_base->med_B),
                _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
            __m512i v16 = _mm512_set1_epi32(16);
            // fb indices of the matching primes, in ascending order;
            // divided out in batches of 64
            ALIGNED_MEM uint32_t mbuf[64 + 16];
            uint32_t nm = 0, k;

            for (i = sconf->factor_base->med_B; (i + 16) <= bound; i += 16)
            {
                __m512i vprime = _mm512_load_si512(fbc->prime + i);
                __m512i vroot1 = _mm512_load_si512(fbc->root1 + i);
                __m512i vroot2 = _mm512_load_si512(fbc->root2 + i);
                __mmask16 msk;

                vprime = _mm512_add_epi32(vprime, vblockloc);
                msk = _mm512_cmpeq_epu32_mask(vprime, vroot1) |
                    _mm512_cmpeq_epu32_mask(vprime, vroot2);
                nm = direct_cstore16(mbuf, nm, msk, vidx, kmode);
                vidx = _mm512_add_epi32(vidx, v16);

                if (nm >= 64)
                {
                    for (k = 0; k < nm; k++)
                    {
                        uint32_t j = mbuf[k];
                        DIVIDE_RESIEVED_PRIME_2(j);
                    }
                    nm = 0;
                }
            }

            for (k = 0; k < nm; k++)
            {
                uint32_t j = mbuf[k];
                DIVIDE_RESIEVED_PRIME_2(j);
            }
        }

        for ( ; i < bound; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            // these primes are all > blocksize, so we only need one step:
            // advance the offset to where it should be after the block, and
            // check whether either root is there now.  If so, this offset is
            // on the progression of the sieve for this prime.
            tmp = block_loc + prime - 32768;
            if (tmp == root1 || tmp == root2)
            {
                //it will divide Q(x).  do so as many times as we can.
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;
    }

    return;
}

void tdiv_large_avx512_bitloop(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    tdiv_large_body(parity, poly_id, bnum, sconf, dconf, SIQS_K_BITLOOP);
}

void tdiv_large_avx512_creg(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    tdiv_large_body(parity, poly_id, bnum, sconf, dconf, SIQS_K_CREG);
}

void tdiv_large_avx512_cmem(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    tdiv_large_body(parity, poly_id, bnum, sconf, dconf, SIQS_K_CMEM);
}

#endif
