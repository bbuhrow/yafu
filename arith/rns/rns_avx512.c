/*
 * RNS Montgomery modular multiplication, AVX-512 IFMA kernels.
 *
 * Eight independent numbers modulo the same N per vector (channel-major
 * layout, see rns.h), so every table entry is a broadcast and each
 * extension-matrix element is loaded once per 8 numbers.  The arithmetic
 * mirrors rns_ref.c step for step and produces bit-identical residues.
 */

#include "rns.h"
#include <stdlib.h>
#include <string.h>

#ifdef RNS_HAVE_AVX512IFMA
#include <immintrin.h>

#define VL RNS_VL

#define LO52(acc, a, b) _mm512_madd52lo_epu64((acc), (a), (b))
#define HI52(acc, a, b) _mm512_madd52hi_epu64((acc), (a), (b))
#define BC(p)           _mm512_set1_epi64((long long)(p))
#define LD(p)           _mm512_load_si512((const void *)(p))
#define ST(p, v)        _mm512_store_si512((void *)(p), (v))

// rows of the extension matrices processed together; 2*RB accumulators
#define RB 4

// reduce lo + hi*2^52 mod 2^52 - d to < 2^52 (see red116 in rns_ref.c)
static inline __m512i vred116(__m512i lo, __m512i hi, __m512i d)
{
    const __m512i m52 = BC(RNS_MASK52);
    const __m512i zero = _mm512_setzero_si512();
    __m512i a = _mm512_and_si512(lo, m52);
    __m512i bv = _mm512_add_epi64(hi, _mm512_srli_epi64(lo, 52));
    __m512i b0 = _mm512_and_si512(bv, m52);
    __m512i b1 = _mm512_srli_epi64(bv, 52);
    __m512i pl = LO52(zero, b0, d);
    __m512i ph = HI52(zero, b0, d);
    __m512i s;

    ph = LO52(ph, b1, d);                               // b1*d + ph < 2^52
    s = _mm512_add_epi64(_mm512_add_epi64(a, pl), LO52(zero, ph, d));
    s = LO52(_mm512_and_si512(s, m52), _mm512_srli_epi64(s, 52), d);
    s = LO52(_mm512_and_si512(s, m52), _mm512_srli_epi64(s, 52), d);
    return s;
}

// x*y mod 2^52 - d, x, y < 2^52 (hi < 2^52 so the b1 and lo>>52 terms vanish)
static inline __m512i vmulred(__m512i x, __m512i y, __m512i d)
{
    const __m512i m52 = BC(RNS_MASK52);
    const __m512i zero = _mm512_setzero_si512();
    __m512i lo = LO52(zero, x, y);
    __m512i hi = HI52(zero, x, y);
    __m512i pl = LO52(zero, hi, d);
    __m512i ph = HI52(zero, hi, d);
    __m512i s;

    s = _mm512_add_epi64(_mm512_add_epi64(lo, pl), LO52(zero, ph, d));
    s = LO52(_mm512_and_si512(s, m52), _mm512_srli_epi64(s, 52), d);
    s = LO52(_mm512_and_si512(s, m52), _mm512_srli_epi64(s, 52), d);
    return s;
}

static inline __m512i vredsmall(__m512i t, __m512i d)
{
    const __m512i m52 = BC(RNS_MASK52);
    t = LO52(_mm512_and_si512(t, m52), _mm512_srli_epi64(t, 52), d);
    t = LO52(_mm512_and_si512(t, m52), _mm512_srli_epi64(t, 52), d);
    return t;
}

void rns_vmulmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    rns_fixed_t *fx = ctx->fx;
    const int nB = fx->nB, nC = fx->nC;
    const __m512i zero = _mm512_setzero_si512();
    const __m512i m52 = BC(RNS_MASK52);
    uint64_t *xi = ctx->vscratch;
    uint64_t *xc = x + nB * VL, *yc = y + nB * VL;
    uint64_t *zc = z + nB * VL;
    __m512i ur, sr, alpha;
    int i, j, r;

    // 1-2: xi_i = |x_i y_i (-N^{-1} B_i^{-1})|_{b_i}
    for (i = 0; i < nB; i++)
    {
        __m512i d = BC(fx->db[i]);
        __m512i t = vmulred(LD(x + i * VL), LD(y + i * VL), d);
        ST(xi + i * VL, vmulred(t, BC(ctx->f1_b[i]), d));
    }

    // redundant channel: U_r = T_r M_B^{-1} + sum_i xi_i BtoCN_ri (mod 2^52).
    // Read x_r, y_r before any z store (z may alias x or y).
    {
        const uint64_t *row = ctx->BtoCN + nC * nB;
        __m512i tr = LO52(zero, LD(x + (nB + nC) * VL), LD(y + (nB + nC) * VL));
        ur = LO52(zero, tr, BC(fx->MBinv_r));
        for (i = 0; i < nB; i++)
            ur = LO52(ur, LD(xi + i * VL), BC(row[i]));
        ur = _mm512_and_si512(ur, m52);
    }

    // 3-4: first extension merged with the Montgomery division, in C.
    // Row block j..j+RB-1 reads xc/yc rows j..j+RB-1 before storing them.
    for (j = 0; j + RB <= nC; j += RB)
    {
        __m512i lo[RB], hi[RB];
        const uint64_t *row[RB];

        for (r = 0; r < RB; r++)
        {
            __m512i t = vmulred(LD(xc + (j + r) * VL), LD(yc + (j + r) * VL),
                BC(fx->dc[j + r]));
            lo[r] = LO52(zero, t, BC(fx->CjMBinv_c[j + r]));
            hi[r] = HI52(zero, t, BC(fx->CjMBinv_c[j + r]));
            row[r] = ctx->BtoCN + (j + r) * nB;
        }
        for (i = 0; i < nB; i++)
        {
            __m512i v = LD(xi + i * VL);
            for (r = 0; r < RB; r++)
            {
                __m512i e = BC(row[r][i]);
                lo[r] = LO52(lo[r], v, e);
                hi[r] = HI52(hi[r], v, e);
            }
        }
        for (r = 0; r < RB; r++)
            ST(zc + (j + r) * VL, vred116(lo[r], hi[r], BC(fx->dc[j + r])));
    }
    for (; j < nC; j++)
    {
        const uint64_t *row = ctx->BtoCN + j * nB;
        __m512i t = vmulred(LD(xc + j * VL), LD(yc + j * VL), BC(fx->dc[j]));
        __m512i lo = LO52(zero, t, BC(fx->CjMBinv_c[j]));
        __m512i hi = HI52(zero, t, BC(fx->CjMBinv_c[j]));
        for (i = 0; i < nB; i++)
        {
            __m512i v = LD(xi + i * VL);
            lo = LO52(lo, v, BC(row[i]));
            hi = HI52(hi, v, BC(row[i]));
        }
        ST(zc + j * VL, vred116(lo, hi, BC(fx->dc[j])));
    }

    // 5: exact second extension C -> B; alpha' from the redundant channel
    {
        const uint64_t *row = fx->CtoB + nB * nC;
        sr = zero;
        for (j = 0; j < nC; j++)
            sr = LO52(sr, LD(zc + j * VL), BC(row[j]));
        alpha = LO52(zero, _mm512_and_si512(_mm512_sub_epi64(sr, ur), m52),
            BC(fx->MCinv_r));
    }

    for (i = 0; i + RB <= nB; i += RB)
    {
        __m512i lo[RB], hi[RB];
        const uint64_t *row[RB];

        for (r = 0; r < RB; r++)
        {
            lo[r] = LO52(zero, alpha, BC(fx->negMC_b[i + r]));
            hi[r] = HI52(zero, alpha, BC(fx->negMC_b[i + r]));
            row[r] = fx->CtoB + (i + r) * nC;
        }
        for (j = 0; j < nC; j++)
        {
            __m512i v = LD(zc + j * VL);
            for (r = 0; r < RB; r++)
            {
                __m512i e = BC(row[r][j]);
                lo[r] = LO52(lo[r], v, e);
                hi[r] = HI52(hi[r], v, e);
            }
        }
        for (r = 0; r < RB; r++)
            ST(z + (i + r) * VL, vred116(lo[r], hi[r], BC(fx->db[i + r])));
    }
    for (; i < nB; i++)
    {
        const uint64_t *row = fx->CtoB + i * nC;
        __m512i lo = LO52(zero, alpha, BC(fx->negMC_b[i]));
        __m512i hi = HI52(zero, alpha, BC(fx->negMC_b[i]));
        for (j = 0; j < nC; j++)
        {
            __m512i v = LD(zc + j * VL);
            lo = LO52(lo, v, BC(row[j]));
            hi = HI52(hi, v, BC(row[j]));
        }
        ST(z + i * VL, vred116(lo, hi, BC(fx->db[i])));
    }

    ST(z + (nB + nC) * VL, ur);
}

void rns_vsqrmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x)
{
    rns_vmulmod_avx512(ctx, z, x, x);
}

void rns_vaddmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    rns_fixed_t *fx = ctx->fx;
    const int nB = fx->nB, nC = fx->nC;
    int ch;

    for (ch = 0; ch < nB + nC; ch++)
    {
        uint64_t d = (ch < nB) ? fx->db[ch] : fx->dc[ch - nB];
        __m512i t = _mm512_add_epi64(LD(x + ch * VL), LD(y + ch * VL));
        ST(z + ch * VL, vredsmall(t, BC(d)));
    }
    ST(z + ch * VL, _mm512_and_si512(
        _mm512_add_epi64(LD(x + ch * VL), LD(y + ch * VL)), BC(RNS_MASK52)));
}

void rns_vsubmod_avx512(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    // X - Y + SUBMULT*k*N; adding 2m keeps every channel sum positive
    rns_fixed_t *fx = ctx->fx;
    const int nB = fx->nB, nC = fx->nC;
    const uint64_t *s = ctx->subconst;
    int ch;

    for (ch = 0; ch < nB + nC; ch++)
    {
        uint64_t m = (ch < nB) ? fx->b[ch] : fx->c[ch - nB];
        __m512i t = _mm512_add_epi64(LD(x + ch * VL), BC(s[ch] + 2 * m));
        t = _mm512_sub_epi64(t, LD(y + ch * VL));
        ST(z + ch * VL, vredsmall(t, BC((1ULL << 52) - m)));
    }
    ST(z + ch * VL, _mm512_and_si512(_mm512_sub_epi64(
        _mm512_add_epi64(LD(x + ch * VL), BC(s[ch])), LD(y + ch * VL)),
        BC(RNS_MASK52)));
}

// window size minimizing squarings-free multiplies, as in AVX-ECM's get_winsize
static int rns_winsize(int bits)
{
    int size, best = 4, minmuls = 0x7fffffff;
    for (size = 2; size <= 8; size++)
    {
        int muls = (bits / size) + (1 << size);
        if (muls < minmuls)
        {
            minmuls = muls;
            best = size;
        }
    }
    return best;
}

void rns_vpowm_avx512(rns_ctx_t *ctx, uint64_t *d, uint64_t *b, mpz_t e)
{
    // left-to-right sliding window over odd powers b^1, b^3, ..., b^(2^k - 1)
    int bit = (int)mpz_sizeinbase(e, 2) - 1;
    int k = rns_winsize(bit);
    int ntab = 1 << (k - 1);
    uint64_t **g = (uint64_t **)calloc((size_t)ntab, sizeof(uint64_t *));
    uint64_t *b2 = rns_valloc(ctx);
    size_t bytes = (size_t)ctx->fx->nch * VL * sizeof(uint64_t);
    int i, j, lane;

    for (i = 0; i < ntab; i++)
        g[i] = rns_valloc(ctx);

    memcpy(g[0], b, bytes);
    rns_vsqrmod_avx512(ctx, b2, b);
    for (i = 1; i < ntab; i++)
        rns_vmulmod_avx512(ctx, g[i], g[i - 1], b2);

    for (lane = 0; lane < VL; lane++)
        rns_vset_lane(ctx, d, lane, ctx->one);

    while (bit >= 0)
    {
        if (mpz_tstbit(e, (mp_bitcnt_t)bit) == 0)
        {
            rns_vsqrmod_avx512(ctx, d, d);
            bit--;
        }
        else
        {
            // longest window ending in a 1 bit
            int a = bit - k + 1, thisk, bstr = 1;
            if (a < 0)
                a = 0;
            a = (int)mpz_scan1(e, (mp_bitcnt_t)a);
            thisk = bit - a + 1;
            for (j = 1; j < thisk; j++)
                bstr = (bstr << 1) | mpz_tstbit(e, (mp_bitcnt_t)(bit - j));
            for (j = 0; j < thisk; j++)
                rns_vsqrmod_avx512(ctx, d, d);
            rns_vmulmod_avx512(ctx, d, d, g[bstr >> 1]);
            bit -= thisk;
        }
    }

    for (i = 0; i < ntab; i++)
        rns_vfree(g[i]);
    free(g);
    rns_vfree(b2);
}

#endif
