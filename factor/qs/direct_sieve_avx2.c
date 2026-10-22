/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// AVX2 kernels for the large-prime strategies (see direct_sieve.c for the
// scalar versions and direct_sieve_avx512.c for the
// AVX-512 ones).  Each 16-prime group is handled as two 8-lane halves, so
// the logp approximation and fb boundaries match the other ISAs, and hits
// and divisions are visited in the same order as the AVX-512 kernels: all
// kernels do bit-identical work.  Compiled into AVX2 and AVX-512 builds.
//
// AVX2 has no compress instruction: SIQS_K_CREG compresses with a
// 256-entry permute table (vpermd), and SIQS_K_CMEM is treated as CREG.
// All values compared here (roots, primes, positions) are < 2^31, so the
// signed compares are exact.

#include <stdint.h>
#include <string.h>
#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "tdiv_macros_common.h"

#if defined(TRY_DIRECT_SIEVE) && defined(USE_AVX2)
#include <immintrin.h>

// lane indices of the set bits of each 8-bit mask, one byte each
static uint64_t direct_avx2_lut[256];
static int direct_avx2_lut_ready = 0;

void direct_avx2_init(void)
{
    uint32_t m, b, n;

    if (direct_avx2_lut_ready)
        return;

    for (m = 0; m < 256; m++)
    {
        uint64_t e = 0;
        for (b = 0, n = 0; b < 8; b++)
        {
            if ((m >> b) & 1)
            {
                e |= (uint64_t)b << (8 * n);
                n++;
            }
        }
        direct_avx2_lut[m] = e;
    }
    direct_avx2_lut_ready = 1;
}

#define MASK8(v) ((uint32_t)_mm256_movemask_ps(_mm256_castsi256_ps(v)))

// append the lanes of v selected by the 8-bit mask m to buf[n...] (buf
// needs 8 entries of room past n); returns the new count
SIQS_FORCEINLINE uint32_t direct_cstore8(uint32_t* buf, uint32_t n, uint32_t m, __m256i v)
{
    __m256i idx = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i*)&direct_avx2_lut[m]));
    _mm256_storeu_si256((__m256i*)(buf + n), _mm256_permutevar8x32_epi32(v, idx));
    return n + siqs_popcnt32(m);
}

// ----------------------------------------------------------------------
// direct scatter: one 16-prime group of [large_B, B) into both intervals
// (pos hits at r, neg hits at interval + p - r), as direct_scatter16
// ----------------------------------------------------------------------
SIQS_FORCEINLINE void direct_scatter16_avx2(uint8_t* sieve, uint32_t interval,
    const uint32_t* p, const uint32_t* r1, const uint32_t* r2, uint8_t logp,
    const int kmode, uint32_t* cbuf, uint32_t* nc, uint32_t* sink)
{
    __m256i vint = _mm256_set1_epi32(interval);
    __m256i vp0 = _mm256_loadu_si256((const __m256i*)p);
    __m256i vp1 = _mm256_loadu_si256((const __m256i*)(p + 8));
    __m256i a0 = _mm256_loadu_si256((const __m256i*)r1);
    __m256i a1 = _mm256_loadu_si256((const __m256i*)(r1 + 8));
    __m256i b0 = _mm256_loadu_si256((const __m256i*)r2);
    __m256i b1 = _mm256_loadu_si256((const __m256i*)(r2 + 8));
    __m256i na0 = _mm256_sub_epi32(vp0, a0);
    __m256i na1 = _mm256_sub_epi32(vp1, a1);
    __m256i nb0 = _mm256_sub_epi32(vp0, b0);
    __m256i nb1 = _mm256_sub_epi32(vp1, b1);
    __m256i ma0 = _mm256_cmpgt_epi32(vint, a0);
    __m256i ma1 = _mm256_cmpgt_epi32(vint, a1);
    __m256i mb0 = _mm256_cmpgt_epi32(vint, b0);
    __m256i mb1 = _mm256_cmpgt_epi32(vint, b1);
    __m256i mna0 = _mm256_cmpgt_epi32(vint, na0);
    __m256i mna1 = _mm256_cmpgt_epi32(vint, na1);
    __m256i mnb0 = _mm256_cmpgt_epi32(vint, nb0);
    __m256i mnb1 = _mm256_cmpgt_epi32(vint, nb1);

    // neg-side positions in the combined buffer
    na0 = _mm256_add_epi32(na0, vint);
    na1 = _mm256_add_epi32(na1, vint);
    nb0 = _mm256_add_epi32(nb0, vint);
    nb1 = _mm256_add_epi32(nb1, vint);

    if ((kmode == SIQS_K_CREG) || (kmode == SIQS_K_CMEM))
    {
        __m256i vtag = _mm256_set1_epi32((uint32_t)logp << 24);
        uint32_t n = *nc;
        n = direct_cstore8(cbuf, n, MASK8(ma0), _mm256_or_si256(a0, vtag));
        n = direct_cstore8(cbuf, n, MASK8(ma1), _mm256_or_si256(a1, vtag));
        n = direct_cstore8(cbuf, n, MASK8(mb0), _mm256_or_si256(b0, vtag));
        n = direct_cstore8(cbuf, n, MASK8(mb1), _mm256_or_si256(b1, vtag));
        n = direct_cstore8(cbuf, n, MASK8(mna0), _mm256_or_si256(na0, vtag));
        n = direct_cstore8(cbuf, n, MASK8(mna1), _mm256_or_si256(na1, vtag));
        n = direct_cstore8(cbuf, n, MASK8(mnb0), _mm256_or_si256(nb0, vtag));
        n = direct_cstore8(cbuf, n, MASK8(mnb1), _mm256_or_si256(nb1, vtag));
        if (n >= DIRECT_CBUF_FLUSH)
        {
            direct_compress_flush(sieve, cbuf, n);
            n = 0;
        }
        *nc = n;
    }
    else if (kmode == SIQS_K_SINK)
    {
        // misses go to a rotating stripe in the second tail guard block
        ALIGNED_MEM uint32_t idx[64];
        __m256i v8 = _mm256_set1_epi32(8);
        __m256i vs = _mm256_add_epi32(_mm256_set1_epi32(2 * interval + 32768 + *sink),
            _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
        int k;

        _mm256_store_si256((__m256i*)(idx + 0), _mm256_blendv_epi8(vs, a0, ma0));
        vs = _mm256_add_epi32(vs, v8);
        _mm256_store_si256((__m256i*)(idx + 8), _mm256_blendv_epi8(vs, a1, ma1));
        vs = _mm256_add_epi32(vs, v8);
        _mm256_store_si256((__m256i*)(idx + 16), _mm256_blendv_epi8(vs, b0, mb0));
        vs = _mm256_add_epi32(vs, v8);
        _mm256_store_si256((__m256i*)(idx + 24), _mm256_blendv_epi8(vs, b1, mb1));
        vs = _mm256_add_epi32(vs, v8);
        _mm256_store_si256((__m256i*)(idx + 32), _mm256_blendv_epi8(vs, na0, mna0));
        vs = _mm256_add_epi32(vs, v8);
        _mm256_store_si256((__m256i*)(idx + 40), _mm256_blendv_epi8(vs, na1, mna1));
        vs = _mm256_add_epi32(vs, v8);
        _mm256_store_si256((__m256i*)(idx + 48), _mm256_blendv_epi8(vs, nb0, mnb0));
        vs = _mm256_add_epi32(vs, v8);
        _mm256_store_si256((__m256i*)(idx + 56), _mm256_blendv_epi8(vs, nb1, mnb1));

        for (k = 0; k < 64; k += 8)
        {
            sieve[idx[k + 0]] -= logp; sieve[idx[k + 1]] -= logp;
            sieve[idx[k + 2]] -= logp; sieve[idx[k + 3]] -= logp;
            sieve[idx[k + 4]] -= logp; sieve[idx[k + 5]] -= logp;
            sieve[idx[k + 6]] -= logp; sieve[idx[k + 7]] -= logp;
        }

        *sink += 64;
        if (*sink >= DIRECT_SINK_BYTES)
            *sink = 0;
    }
    else
    {
        ALIGNED_MEM uint32_t t[8];
        uint32_t m;

#define DIRECT_BITLOOP8(vpos, vm) \
        m = MASK8(vm); \
        if (m) { \
            _mm256_store_si256((__m256i*)t, vpos); \
            while (m) { sieve[t[siqs_ctz32(m)]] -= logp; m &= m - 1; } \
        }

        DIRECT_BITLOOP8(a0, ma0);
        DIRECT_BITLOOP8(a1, ma1);
        DIRECT_BITLOOP8(b0, mb0);
        DIRECT_BITLOOP8(b1, mb1);
        DIRECT_BITLOOP8(na0, mna0);
        DIRECT_BITLOOP8(na1, mna1);
        DIRECT_BITLOOP8(nb0, mnb0);
        DIRECT_BITLOOP8(nb1, mnb1);
#undef DIRECT_BITLOOP8
    }
}

// direct pass 1
SIQS_FORCEINLINE void direct_sieve_largeprimes_avx2_body(uint8_t* sieve,
    static_conf_t* sconf, dynamic_conf_t* dconf, const int kmode)
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
        direct_scatter16_avx2(sieve, interval, primes + j, fr1 + j, fr2 + j,
            logp, kmode, cbuf, &nc, &sink);
    }
    direct_compress_flush(sieve, cbuf, nc);

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
}

void direct_sieve_largeprimes_avx2_bitloop(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_avx2_body(sieve, sconf, dconf, SIQS_K_BITLOOP);
}

void direct_sieve_largeprimes_avx2_sink(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_avx2_body(sieve, sconf, dconf, SIQS_K_SINK);
}

void direct_sieve_largeprimes_avx2_creg(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_avx2_body(sieve, sconf, dconf, SIQS_K_CREG);
}

// ----------------------------------------------------------------------
// fused resieve / root update / next-poly scatter
// ----------------------------------------------------------------------

// survivor test for 8 lanes at positions vpos with in-range lane mask vin:
// the 512-bit register prefilter (one bit per 2^shift positions), then a
// masked gather of the bitmap words for the lanes left.  Returns the
// 8-bit mask of lanes that sit on a survivor.
SIQS_FORCEINLINE uint32_t direct_rs8(__m256i vpos, __m256i vin, __m256i vpf_lo,
    __m256i vpf_hi, __m256i vshift, const int* bitmap32)
{
    __m256i vone = _mm256_set1_epi32(1);
    __m256i v31 = _mm256_set1_epi32(31);
    __m256i vb = _mm256_srlv_epi32(vpos, vshift);
    __m256i widx = _mm256_srli_epi32(vb, 5);
    __m256i w = _mm256_blendv_epi8(_mm256_permutevar8x32_epi32(vpf_lo, widx),
        _mm256_permutevar8x32_epi32(vpf_hi, widx),
        _mm256_cmpgt_epi32(widx, _mm256_set1_epi32(7)));
    __m256i m = _mm256_and_si256(vin, _mm256_cmpeq_epi32(vone,
        _mm256_and_si256(_mm256_srlv_epi32(w, _mm256_and_si256(vb, v31)), vone)));

    if (MASK8(m) == 0)
        return 0;

    w = _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), bitmap32,
        _mm256_srli_epi32(vpos, 5), m, 4);
    m = _mm256_and_si256(m, _mm256_cmpeq_epi32(vone,
        _mm256_and_si256(_mm256_srlv_epi32(w, _mm256_and_si256(vpos, v31)), vone)));
    return MASK8(m);
}

SIQS_FORCEINLINE void direct_resieve_update_avx2_body(static_conf_t* sconf,
    dynamic_conf_t* dconf, uint8_t* sieve, const int kmode)
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
    int resieve = (pool->num > 0);
    uint8_t* sieven = (sieve == NULL) ? NULL : sieve + interval;
    const int* bitmap32 = (const int*)pool->bitmap;
    uint32_t j;
    uint8_t logp;
    ALIGNED_MEM uint32_t cbuf[DIRECT_CBUF_SIZE];
    ALIGNED_MEM uint32_t pf_words[16];
    ALIGNED_MEM uint32_t pos16[16];
    uint32_t nc = 0, sink = 0;
    __m256i vint = _mm256_set1_epi32(interval);
    __m256i vone = _mm256_set1_epi32(1);
    __m256i vpf_lo, vpf_hi, vshift;

    if ((sieve == NULL) && (!resieve))
        return;

    // 512-bit prefilter over both intervals, built from this poly's pool
    {
        uint32_t k, pf_shift = 0;
        while (((uint64_t)512 << pf_shift) < (2 * (uint64_t)interval))
            pf_shift++;
        memset(pf_words, 0, sizeof(pf_words));
        for (k = 0; k < pool->num; k++)
        {
            uint32_t b = pool->pos[k] >> pf_shift;
            pf_words[b >> 5] |= (1u << (b & 31));
        }
        vpf_lo = _mm256_load_si256((const __m256i*)pf_words);
        vpf_hi = _mm256_load_si256((const __m256i*)(pf_words + 8));
        vshift = _mm256_set1_epi32(pf_shift);
    }

    for (j = large_B; (j + 16) <= bound; j += 16, ptr += 16)
    {
        __m256i vp0 = _mm256_loadu_si256((const __m256i*)(primes + j));
        __m256i vp1 = _mm256_loadu_si256((const __m256i*)(primes + j + 8));
        __m256i a0 = _mm256_loadu_si256((const __m256i*)(fr1 + j));
        __m256i a1 = _mm256_loadu_si256((const __m256i*)(fr1 + j + 8));
        __m256i b0 = _mm256_loadu_si256((const __m256i*)(fr2 + j));
        __m256i b1 = _mm256_loadu_si256((const __m256i*)(fr2 + j + 8));

        if (resieve)
        {
            // per root kind (pos root1, pos root2, neg root1, neg root2),
            // lanes 0..15 in order: the AVX-512 visiting order
            uint32_t h;

#define DIRECT_RS16(lo, hi, neg) \
            { \
                __m256i q0 = lo, q1 = hi, in0, in1; \
                if (neg) { q0 = _mm256_sub_epi32(vp0, q0); q1 = _mm256_sub_epi32(vp1, q1); } \
                in0 = _mm256_cmpgt_epi32(vint, q0); \
                in1 = _mm256_cmpgt_epi32(vint, q1); \
                if (neg) { q0 = _mm256_add_epi32(q0, vint); q1 = _mm256_add_epi32(q1, vint); } \
                h = direct_rs8(q0, in0, vpf_lo, vpf_hi, vshift, bitmap32) | \
                    (direct_rs8(q1, in1, vpf_lo, vpf_hi, vshift, bitmap32) << 8); \
                if (h) \
                { \
                    _mm256_store_si256((__m256i*)pos16, q0); \
                    _mm256_store_si256((__m256i*)(pos16 + 8), q1); \
                    while (h) \
                    { \
                        uint32_t k = siqs_ctz32(h); \
                        direct_resieve_hit(pool, pos16[k], j + k, primes[j + k]); \
                        h &= h - 1; \
                    } \
                } \
            }

            DIRECT_RS16(a0, a1, 0);
            DIRECT_RS16(b0, b1, 0);
            DIRECT_RS16(a0, a1, 1);
            DIRECT_RS16(b0, b1, 1);
#undef DIRECT_RS16
        }

        if (sieve == NULL)
            continue;

        // advance to the next poly; the pair stays unordered
        {
            __m256i u0 = _mm256_loadu_si256((const __m256i*)ptr);
            __m256i u1 = _mm256_loadu_si256((const __m256i*)(ptr + 8));

            if (sign > 0)
            {
                __m256i m;
                m = _mm256_cmpgt_epi32(u0, a0);
                a0 = _mm256_add_epi32(_mm256_sub_epi32(a0, u0), _mm256_and_si256(vp0, m));
                m = _mm256_cmpgt_epi32(u1, a1);
                a1 = _mm256_add_epi32(_mm256_sub_epi32(a1, u1), _mm256_and_si256(vp1, m));
                m = _mm256_cmpgt_epi32(u0, b0);
                b0 = _mm256_add_epi32(_mm256_sub_epi32(b0, u0), _mm256_and_si256(vp0, m));
                m = _mm256_cmpgt_epi32(u1, b1);
                b1 = _mm256_add_epi32(_mm256_sub_epi32(b1, u1), _mm256_and_si256(vp1, m));
            }
            else
            {
                __m256i pm0 = _mm256_sub_epi32(vp0, vone);
                __m256i pm1 = _mm256_sub_epi32(vp1, vone);
                a0 = _mm256_add_epi32(a0, u0);
                a1 = _mm256_add_epi32(a1, u1);
                b0 = _mm256_add_epi32(b0, u0);
                b1 = _mm256_add_epi32(b1, u1);
                a0 = _mm256_sub_epi32(a0, _mm256_and_si256(vp0, _mm256_cmpgt_epi32(a0, pm0)));
                a1 = _mm256_sub_epi32(a1, _mm256_and_si256(vp1, _mm256_cmpgt_epi32(a1, pm1)));
                b0 = _mm256_sub_epi32(b0, _mm256_and_si256(vp0, _mm256_cmpgt_epi32(b0, pm0)));
                b1 = _mm256_sub_epi32(b1, _mm256_and_si256(vp1, _mm256_cmpgt_epi32(b1, pm1)));
            }
            _mm256_storeu_si256((__m256i*)(fr1 + j), a0);
            _mm256_storeu_si256((__m256i*)(fr1 + j + 8), a1);
            _mm256_storeu_si256((__m256i*)(fr2 + j), b0);
            _mm256_storeu_si256((__m256i*)(fr2 + j + 8), b1);
        }

        // sieve the next poly's roots, both sides (as pass 1)
        logp = update_data.logp[j + 15] + DIRECT_LOGP_BIAS;
        direct_scatter16_avx2(sieve, interval, primes + j, fr1 + j, fr2 + j,
            logp, kmode, cbuf, &nc, &sink);
    }

    // nc stays 0 when sieve == NULL
    direct_compress_flush(sieve, cbuf, nc);

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
}

void direct_resieve_update_avx2_bitloop(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    direct_resieve_update_avx2_body(sconf, dconf, sieve, SIQS_K_BITLOOP);
}

void direct_resieve_update_avx2_sink(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    direct_resieve_update_avx2_body(sconf, dconf, sieve, SIQS_K_SINK);
}

void direct_resieve_update_avx2_creg(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    direct_resieve_update_avx2_body(sconf, dconf, sieve, SIQS_K_CREG);
}

// ----------------------------------------------------------------------
// sparse sieve and its tdiv, [med_B, fb_blk8i)
// ----------------------------------------------------------------------
SIQS_FORCEINLINE void large_sieveblock_avx2_body(uint8_t* sieve, sieve_fb* fb,
    fb_list* full_fb, const int kmode)
{
    uint32_t i;
    uint32_t prime, root1, root2;
    uint8_t logp;
    uint32_t bound = full_fb->fb_blk8i;
    __m256i vblock = _mm256_set1_epi32(32768);
    ALIGNED_MEM uint32_t cbuf[DIRECT_CBUF_SIZE];
    ALIGNED_MEM uint32_t t[8];
    uint32_t nc = 0;

    for (i = full_fb->med_B; (i + 16) <= bound; i += 16)
    {
        int h;
        // approximate the next 16 logp's as equal to the largest
        logp = fb->logp[i + 15];

        for (h = 0; h < 16; h += 8)
        {
            __m256i vp = _mm256_loadu_si256((const __m256i*)(fb->prime + i + h));
            __m256i v1 = _mm256_loadu_si256((const __m256i*)(fb->root1 + i + h));
            __m256i v2 = _mm256_loadu_si256((const __m256i*)(fb->root2 + i + h));
            __m256i m1 = _mm256_cmpgt_epi32(vblock, v1);
            __m256i m2 = _mm256_cmpgt_epi32(vblock, v2);

            if (kmode == SIQS_K_BITLOOP)
            {
                uint32_t m;
                m = MASK8(m1);
                if (m)
                {
                    _mm256_store_si256((__m256i*)t, v1);
                    while (m) { sieve[t[siqs_ctz32(m)]] -= logp; m &= m - 1; }
                }
                m = MASK8(m2);
                if (m)
                {
                    _mm256_store_si256((__m256i*)t, v2);
                    while (m) { sieve[t[siqs_ctz32(m)]] -= logp; m &= m - 1; }
                }
            }
            else
            {
                __m256i vtag = _mm256_set1_epi32((uint32_t)logp << 24);
                nc = direct_cstore8(cbuf, nc, MASK8(m1), _mm256_or_si256(v1, vtag));
                nc = direct_cstore8(cbuf, nc, MASK8(m2), _mm256_or_si256(v2, vtag));
            }

            // advance the roots that hit, reduce both for the next block
            v1 = _mm256_sub_epi32(_mm256_add_epi32(v1, _mm256_and_si256(vp, m1)), vblock);
            v2 = _mm256_sub_epi32(_mm256_add_epi32(v2, _mm256_and_si256(vp, m2)), vblock);
            _mm256_storeu_si256((__m256i*)(fb->root1 + i + h), _mm256_min_epu32(v1, v2));
            _mm256_storeu_si256((__m256i*)(fb->root2 + i + h), _mm256_max_epu32(v1, v2));
        }

        if (nc >= DIRECT_CBUF_FLUSH)
        {
            direct_compress_flush(sieve, cbuf, nc);
            nc = 0;
        }
    }
    direct_compress_flush(sieve, cbuf, nc);

    // primes left over after the 16-groups; same as the AVX-512 tail
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

        if (root1 < (uint32_t)(32768))
        {
            sieve[root1] -= logp;
            root1 += prime;
        }

        root1 -= 32768;
        root2 -= 32768;
        fb->root1[i] = (root1 < root2) ? root1 : root2;
        fb->root2[i] = (root1 < root2) ? root2 : root1;
    }
}

void large_sieveblock_avx2_bitloop(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    large_sieveblock_avx2_body(sieve, fb, full_fb, SIQS_K_BITLOOP);
}

void large_sieveblock_avx2_creg(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    large_sieveblock_avx2_body(sieve, fb, full_fb, SIQS_K_CREG);
}

SIQS_FORCEINLINE void tdiv_large_avx2_body(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf, const int kmode)
{
    uint32_t i;
    uint32_t tmp, prime, root1, root2, report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb* fbc = parity ? dconf->sieve_n : dconf->sieve_p;
    uint32_t block_loc;
    uint32_t bound = sconf->factor_base->fb_blk8i;

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {
        __m256i vbl;
        ALIGNED_MEM uint32_t mbuf[64 + 16];
        uint32_t nm = 0, k;

        if (!dconf->valid_Qs[report_num])
            continue;

        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];
        vbl = _mm256_set1_epi32(block_loc - 32768);

        for (i = sconf->factor_base->med_B; (i + 16) <= bound; i += 16)
        {
            uint32_t m[2];
            int h;

            for (h = 0; h < 2; h++)
            {
                __m256i t = _mm256_add_epi32(_mm256_loadu_si256(
                    (const __m256i*)(fbc->prime + i + 8 * h)), vbl);
                __m256i e = _mm256_or_si256(
                    _mm256_cmpeq_epi32(t, _mm256_loadu_si256((const __m256i*)(fbc->root1 + i + 8 * h))),
                    _mm256_cmpeq_epi32(t, _mm256_loadu_si256((const __m256i*)(fbc->root2 + i + 8 * h))));
                m[h] = MASK8(e);
            }

            if (kmode == SIQS_K_BITLOOP)
            {
                uint32_t msk = m[0] | (m[1] << 8);
                while (msk)
                {
                    uint32_t j = i + siqs_ctz32(msk);
                    DIVIDE_RESIEVED_PRIME_2(j);
                    msk &= msk - 1;
                }
            }
            else
            {
                __m256i vidx = _mm256_add_epi32(_mm256_set1_epi32(i),
                    _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
                nm = direct_cstore8(mbuf, nm, m[0], vidx);
                nm = direct_cstore8(mbuf, nm, m[1], _mm256_add_epi32(vidx, _mm256_set1_epi32(8)));
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
        }

        for (k = 0; k < nm; k++)
        {
            uint32_t j = mbuf[k];
            DIVIDE_RESIEVED_PRIME_2(j);
        }

        for (; i < bound; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            // these primes are all > blocksize, so we only need one step
            tmp = block_loc + prime - 32768;
            if (tmp == root1 || tmp == root2)
            {
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;
    }
}

void tdiv_large_avx2_bitloop(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    tdiv_large_avx2_body(parity, poly_id, bnum, sconf, dconf, SIQS_K_BITLOOP);
}

void tdiv_large_avx2_creg(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    tdiv_large_avx2_body(parity, poly_id, bnum, sconf, dconf, SIQS_K_CREG);
}

// ----------------------------------------------------------------------
// per-poly roots for [med_B, B): vector sparse update + portable buckets
// ----------------------------------------------------------------------
void nextRoots_32k_lp_sparse_avx2(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    fb_list* fb = sconf->factor_base;
    update_t update_data = dconf->update_data;
    uint32_t bound = fb->B;
    uint32_t med_B = fb->med_B;
    uint32_t sparse_B = fb->fb_blk8i;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    uint32_t* ptr = (uint32_t*)&dconf->rootupdates[(v - 1) * bound + med_B];
    uint32_t* fr1 = (uint32_t*)update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)update_data.firstroots2;
    uint32_t j, root1, root2, prime;
    __m256i vone = _mm256_set1_epi32(1);

    for (j = med_B; (j + 8) <= sparse_B; j += 8, ptr += 8)
    {
        __m256i vp = _mm256_loadu_si256((const __m256i*)(update_data.prime + j));
        __m256i r1 = _mm256_loadu_si256((const __m256i*)(fr1 + j));
        __m256i r2 = _mm256_loadu_si256((const __m256i*)(fr2 + j));
        __m256i u = _mm256_loadu_si256((const __m256i*)ptr);
        __m256i lo, hi;

        if (sign > 0)
        {
            __m256i m1 = _mm256_cmpgt_epi32(u, r1);
            __m256i m2 = _mm256_cmpgt_epi32(u, r2);
            r1 = _mm256_add_epi32(_mm256_sub_epi32(r1, u), _mm256_and_si256(vp, m1));
            r2 = _mm256_add_epi32(_mm256_sub_epi32(r2, u), _mm256_and_si256(vp, m2));
        }
        else
        {
            __m256i pm = _mm256_sub_epi32(vp, vone);
            r1 = _mm256_add_epi32(r1, u);
            r2 = _mm256_add_epi32(r2, u);
            r1 = _mm256_sub_epi32(r1, _mm256_and_si256(vp, _mm256_cmpgt_epi32(r1, pm)));
            r2 = _mm256_sub_epi32(r2, _mm256_and_si256(vp, _mm256_cmpgt_epi32(r2, pm)));
        }

        lo = _mm256_min_epu32(r1, r2);
        hi = _mm256_max_epu32(r1, r2);
        _mm256_storeu_si256((__m256i*)(fr1 + j), lo);
        _mm256_storeu_si256((__m256i*)(fr2 + j), hi);
        _mm256_storeu_si256((__m256i*)(dconf->sieve_p->root1 + j), lo);
        _mm256_storeu_si256((__m256i*)(dconf->sieve_p->root2 + j), hi);
        _mm256_storeu_si256((__m256i*)(dconf->sieve_n->root1 + j), _mm256_sub_epi32(vp, hi));
        _mm256_storeu_si256((__m256i*)(dconf->sieve_n->root2 + j), _mm256_sub_epi32(vp, lo));
    }

    for (; j < sparse_B; j++, ptr++)
    {
        uint32_t uu = *ptr;
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

        if (root2 < root1) { uint32_t t = root1; root1 = root2; root2 = t; }
        fr1[j] = root1;
        fr2[j] = root2;
        dconf->sieve_p->root1[j] = root1;
        dconf->sieve_p->root2[j] = root2;
        dconf->sieve_n->root1[j] = prime - root2;
        dconf->sieve_n->root2[j] = prime - root1;
    }
}

void nextRoots_32k_lp_avx2(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    nextRoots_32k_lp_sparse_avx2(sconf, dconf);
    nextRoots_32k_lp_buckets(sconf, dconf);
}

#endif
