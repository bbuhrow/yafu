/*
 * vecpoly kernels: residue conversion, NTTs, pointwise products, CRT and
 * coefficient add/sub, for 8-lane polynomials (see vecpoly.h).  One source
 * for every target through the rns_vec52.h vector layer:
 *   AVX-512 + IFMA   native 52-bit products
 *   AVX-512F         exact FP-emulated products (embedded round-to-zero)
 *   AVX2 + FMA       two 4-lane halves per 8-lane vector, FP-emulated
 *                    products with MXCSR set to round toward zero for the
 *                    duration of each call (build with -frounding-math)
 *
 * Residues are lazy, < 2^52 (p < 2^50): Shoup products by constants give
 * [0, 2p), Montgomery pointwise products [0, 2p), the forward NTT keeps
 * [0, 2p), the inverse NTT gives [0, 4p).
 */

#include "vecpoly.h"

#ifdef VP_HAVE_KERNELS

#if defined(RNS_HAVE_AVX2) && !defined(RNS_HAVE_AVX512)
#if defined(_MSC_VER) && !defined(__clang__)
#pragma fenv_access(on)
#elif defined(__clang__)
#pragma STDC FENV_ACCESS ON
#endif
#define RNS_VW 256
#else
#define RNS_VW 512
#endif
#include "../rns/rns_vec52.h"

#define VP_H (8 / V_LANES)      // vectors per 8-lane element
#define VP_MAXW 512             // max words per coefficient (stack arrays)

#if RNS_VW == 512
// x >= m ? x - m : x   (x, m < 2^63)
#define V_CSUB(x, m)        _mm512_min_epu64((x), _mm512_sub_epi64((x), (m)))
// 1 where t != 0, else 0
#define V_NZ1(t)            _mm512_maskz_set1_epi64(_mm512_test_epi64_mask((t), (t)), 1)
// per lane: sign(s) ? b : a
#define V_BLEND_SIGN(a, b, s) _mm512_mask_blend_epi64(_mm512_movepi64_mask(s), (a), (b))
#define V_CVTD(a)           _mm512_cvtepu64_pd(a)
#define V_FMAD(a, b, c)     _mm512_fmadd_pd((a), (b), (c))
#define V_ZEROD             _mm512_setzero_pd()
// floor(s + 1/8) for 0 <= s < 2^51
#define V_KAPPA(s)          _mm512_cvttpd_epu64(_mm512_add_pd((s), _mm512_set1_pd(0.125)))
#define VP_RZ_ENTER()
#define VP_RZ_LEAVE()
#else
#ifdef VP_INT_NTT
RNS_INLINE __m256i vp_csub256(__m256i x, __m256i m)
{
    __m256i t = _mm256_sub_epi64(x, m);
    return _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(t),
        _mm256_castsi256_pd(x), _mm256_castsi256_pd(t)));
}
#endif
#define V_CSUB(x, m)        vp_csub256((x), (m))
#define V_NZ1(t)            _mm256_add_epi64(_mm256_set1_epi64x(1), \
                                _mm256_cmpeq_epi64((t), _mm256_setzero_si256()))
#define V_BLEND_SIGN(a, b, s) _mm256_castpd_si256(_mm256_blendv_pd( \
                                _mm256_castsi256_pd(a), _mm256_castsi256_pd(b), _mm256_castsi256_pd(s)))
#define V_CVTD(a)           V_CVT52(a)
#define V_FMAD(a, b, c)     _mm256_fmadd_pd((a), (b), (c))
#define V_ZEROD             _mm256_setzero_pd()
// MXCSR is round-toward-zero: (s + 1/8) + 2^52 truncates to floor
#define V_KAPPA(s)          V_SUB(V_D2I(_mm256_add_pd(_mm256_add_pd((s), _mm256_set1_pd(0.125)), \
                                V_I2D(V_BC(RNS_FP_BIAS_LO)))), V_BC(RNS_FP_BIAS_LO))
#define VP_RZ_ENTER()       unsigned int vp_csr_ = rns_fp_rz_enter()
#define VP_RZ_LEAVE()       rns_fp_rz_leave(vp_csr_)
#endif

#define VM52 V_BC(VP_MASK52)

// products with pre-converted operands (FP: doubles, biased accumulation)
#if RNS_VEC_FP
typedef rns_vd vp_op;
#define VP_OP(x)            rns_cvt52(x)
#define VP_OPC(c)           rns_bcd52(c)
// constant table entry: tabd = the table as doubles (ctx->cvtd / crtd)
#define VP_OPT(tabi, tabd, k)  ((void)(tabi), V_BCD((tabd)[k]))
#define VP_MAC(l, h, x, y)  RNS_FP_MAC(l, h, x, y)
#define VP_UNBIAS_LO(l, cnt) ((l) = rns_fp_unbias((l), (cnt), 0))
#define VP_UNBIAS_HI(h, cnt) ((h) = rns_fp_unbias((h), 0, (cnt)))
#else
typedef rns_vi vp_op;
#define VP_OP(x)            (x)
#define VP_OPC(c)           V_BC(c)
#define VP_OPT(tabi, tabd, k)  ((void)(tabd), V_BC((tabi)[k]))
#define VP_MAC(l, h, x, y)  do { (l) = LO52((l), (x), (y)); (h) = HI52((h), (x), (y)); } while (0)
#define VP_UNBIAS_LO(l, cnt)
#define VP_UNBIAS_HI(h, cnt)
#endif

#if !RNS_VEC_FP || defined(VP_INT_NTT)       // integer NTT only
// Shoup: x * w mod p in [0, 2p), x < 2^52, wp = floor(w 2^52 / p)
RNS_INLINE rns_vi vp_shoup(rns_vi x, rns_vi w, rns_vi wp, rns_vi p)
{
    rns_vi q = HI52(V_ZERO, x, wp);
    rns_vi r = V_SUB(LO52(V_ZERO, x, w), LO52(V_ZERO, q, p));
    return V_AND(r, VM52);
}

// Montgomery: a * b * 2^-52 mod p in [0, 2p), a, b < 2p
RNS_INLINE rns_vi vp_mont(rns_vi a, rns_vi b, rns_vi p, rns_vi pneg)
{
    rns_vi lo = LO52(V_ZERO, a, b);
    rns_vi hi = HI52(V_ZERO, a, b);
    rns_vi m = LO52(V_ZERO, lo, pneg);
    return V_ADD(HI52(hi, m, p), V_NZ1(lo));
}
#endif

/*
 * Native floating-point NTT (non-IFMA builds; VP_INT_NTT forces the
 * integer one): NTT-domain values are doubles holding integers x with
 * |x| <= p (p < 2^50) after every kernel.  Products modulo p with FMA
 * (as FLINT's fft_small): h = a*b, l = a*b - h exactly, q = round(h/p),
 * r = (h - q p) + l = a*b - q p; every step is exact, |r| < p.
 * Twiddles and CRT scales are (w, w/p) pairs of doubles, |w| <= p/2.
 * The pointwise product is plain (no 2^-52), so the CRT scales drop 2^52.
 */
#if RNS_VEC_FP && !defined(VP_INT_NTT)
#define VP_FPNTT 1
#if RNS_VW == 512
#define D_ADD(a, b)         _mm512_add_pd((a), (b))
#define D_SUB(a, b)         _mm512_sub_pd((a), (b))
#define D_MUL(a, b)         _mm512_mul_pd((a), (b))
#define D_FMSUB(a, b, c)    _mm512_fmsub_pd((a), (b), (c))      // a b - c
#define D_FNMADD(a, b, c)   _mm512_fnmadd_pd((a), (b), (c))     // c - a b
#define D_RND(x)            _mm512_roundscale_pd((x), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
// y < 0 ? y + p : y
#define D_ADDNEG(y, p)      _mm512_mask_add_pd((y), _mm512_cmp_pd_mask((y), _mm512_setzero_pd(), _CMP_LT_OQ), (y), (p))
#else
#define D_ADD(a, b)         _mm256_add_pd((a), (b))
#define D_SUB(a, b)         _mm256_sub_pd((a), (b))
#define D_MUL(a, b)         _mm256_mul_pd((a), (b))
#define D_FMSUB(a, b, c)    _mm256_fmsub_pd((a), (b), (c))
#define D_FNMADD(a, b, c)   _mm256_fnmadd_pd((a), (b), (c))
#define D_RND(x)            _mm256_round_pd((x), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define D_ADDNEG(y, p)      _mm256_add_pd((y), _mm256_and_pd(_mm256_cmp_pd((y), _mm256_setzero_pd(), _CMP_LT_OQ), (p)))
#endif

// x - p round(x/p): |result| <= p/2 (+ rounding), x an integer < 2^53
RNS_INLINE rns_vd vpd_red(rns_vd x, rns_vd p, rns_vd pinv)
{
    return D_FNMADD(D_RND(D_MUL(x, pinv)), p, x);
}

// x w mod p for a constant w (|w| <= p/2, wp = w/p), |x| < 2^52:
// |result| <= 0.75 p for |x| <= 2p
RNS_INLINE rns_vd vpd_mulc(rns_vd x, rns_vd w, rns_vd wp, rns_vd p)
{
    rns_vd h = D_MUL(x, w);
    rns_vd l = D_FMSUB(x, w, h);
    rns_vd q = D_RND(D_MUL(x, wp));
    return D_ADD(D_FNMADD(q, p, h), l);
}

// a b mod p, |a|, |b| <= p: |result| <= 0.875 p
RNS_INLINE rns_vd vpd_mulm(rns_vd a, rns_vd b, rns_vd p, rns_vd pinv)
{
    rns_vd h = D_MUL(a, b);
    rns_vd l = D_FMSUB(a, b, h);
    rns_vd q = D_RND(D_MUL(h, pinv));
    return D_ADD(D_FNMADD(q, p, h), l);
}
#else
#define VP_FPNTT 0
#endif

void vpk_to_ntt_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L,
    const uint64_t *a, uint32_t len, int i0, int i1)
{
    int n = ctx->n, i, k, h;
    uint32_t t;
    vp_op d[VP_MAXW];
    VP_RZ_ENTER();

    for (t = 0; t < len; t++)
    {
        const uint64_t *c = a + (size_t)t * 8 * (size_t)n;
        for (h = 0; h < VP_H; h++)
        {
            int o = h * V_LANES;
            for (k = 0; k < n; k++)
                d[k] = VP_OP(V_LD(c + k * 8 + o));
            for (i = i0; i < i1; i++)
            {
                const uint64_t *cv = ctx->cvt + (size_t)i * (size_t)n;
                const double *cvd = ctx->cvtd + (size_t)i * (size_t)n;
                rns_vi p = V_BC(ctx->p[i]), pneg = V_BC(ctx->pneg[i]);
                rns_vi lo = V_ZERO, hi = V_ZERO, m, r, h0;
                for (k = 0; k < n; k++)
                    VP_MAC(lo, hi, d[k], VP_OPT(cv, cvd, k));
                VP_UNBIAS_LO(lo, n);
                VP_UNBIAS_HI(hi, n);
                // x 2^104 = hi 2^52 + lo; REDC by 2^104
                hi = V_ADD(hi, V_SRLI(lo, 52));
                lo = V_AND(lo, VM52);
                m = LO52(V_ZERO, lo, pneg);
                hi = V_ADD(HI52(hi, m, p), V_NZ1(lo));
                h0 = V_AND(hi, VM52);
                m = LO52(V_ZERO, h0, pneg);
                r = V_ADD(HI52(V_SRLI(hi, 52), m, p), V_NZ1(h0));
#if VP_FPNTT
                // r in [0, 2p) -> |r| <= p/2, as a double
                V_STD(VP_RES(v, L, i - i0, t) + o, vpd_red(V_CVT52(r),
                    V_BCD((double)ctx->p[i]), V_BCD(ctx->pinvd[i])));
#else
                V_ST(VP_RES(v, L, i - i0, t) + o, r);
#endif
            }
        }
    }
    for (i = i0; i < i1; i++)
        for (t = len; t < L; t++)
            for (h = 0; h < VP_H; h++)
                V_ST(VP_RES(v, L, i - i0, t) + h * V_LANES, V_ZERO);
    VP_RZ_LEAVE();
}

#if VP_FPNTT

/*
 * Transforms: radix-4 passes (two radix-2 levels per sweep over the data,
 * half the memory traffic), plus one radix-2 level when log2 L is odd
 * (first for DIF, last for DIT).  No MXCSR switch: every step is exact or
 * explicitly rounded.  Bounds, inputs |.| <= p:
 *   DIF, two levels: sums unreduced after the first (<= 2p), the second
 *     level reduces its sum (<= 4p -> p/2) and multiplies differences of
 *     at most 4p (-> 0.75p): out |.| <= p.
 *   DIT, two levels: first level unreduced (<= 1.5p), second level reduces
 *     both outputs (<= 2p -> p/2): out |.| <= p.
 */

// DIF butterfly: x' = x + y (reduced if red), y' = (x - y) w
#define VPD_FWD2(x, y, w, wp, red) do {                                 \
        rns_vd s_ = D_ADD((x), (y));                                    \
        (y) = vpd_mulc(D_SUB((x), (y)), (w), (wp), p);                  \
        (x) = (red) ? vpd_red(s_, p, pinv) : s_;                        \
    } while (0)

// DIT butterfly: t = y w, x' = x + t, y' = x - t (reduced if red)
#define VPD_INV2(x, y, w, wp, red) do {                                 \
        rns_vd t_ = vpd_mulc((y), (w), (wp), p);                        \
        rns_vd a_ = D_ADD((x), t_), b_ = D_SUB((x), t_);                \
        (x) = (red) ? vpd_red(a_, p, pinv) : a_;                        \
        (y) = (red) ? vpd_red(b_, p, pinv) : b_;                        \
    } while (0)

static int vp_log2(uint32_t L)
{
    int l = 0;
    while ((1u << l) < L) l++;
    return l;
}

void vpk_fwd_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L, int i0, int i1)
{
    int i, h, lg = vp_log2(L);
    uint32_t m, b, j, st;

    for (i = i0; i < i1; i++)
    {
        rns_vd p = V_BCD((double)ctx->p[i]), pinv = V_BCD(ctx->pinvd[i]);
        double *base = (double *)(v + (size_t)(i - i0) * L * 8);
        const double *tw = (const double *)ctx->tw + (size_t)i * ctx->Lmax;
        m = L >> 1;
        if (lg & 1)
        {
            // one radix-2 level (twiddles w_2m^j, table stride Lmax/2m)
            st = ctx->Lmax / (2 * m);
            for (b = 0; b < L; b += 2 * m)
                for (j = 0; j < m; j++)
                {
                    rns_vd w = V_BCD(tw[2 * j * st]), wp = V_BCD(tw[2 * j * st + 1]);
                    double *px = base + (size_t)(b + j) * 8, *py = px + (size_t)m * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vd x = V_LDD(px + o), y = V_LDD(py + o);
                        VPD_FWD2(x, y, w, wp, 1);
                        V_STD(px + o, x);
                        V_STD(py + o, y);
                    }
                }
            m >>= 1;
        }
        for (; m >= 2; m >>= 2)
        {
            // levels m and m/2 on a0..a3 = b + j + {0, m/2, m, 3m/2}
            uint32_t q = m / 2;
            st = ctx->Lmax / (2 * m);
            for (b = 0; b < L; b += 2 * m)
                for (j = 0; j < q; j++)
                {
                    rns_vd wa = V_BCD(tw[2 * j * st]), wap = V_BCD(tw[2 * j * st + 1]);
                    rns_vd wb = V_BCD(tw[2 * (j + q) * st]), wbp = V_BCD(tw[2 * (j + q) * st + 1]);
                    rns_vd wc = V_BCD(tw[4 * j * st]), wcp = V_BCD(tw[4 * j * st + 1]);
                    double *p0 = base + (size_t)(b + j) * 8, *p1 = p0 + (size_t)q * 8;
                    double *p2 = p0 + (size_t)m * 8, *p3 = p2 + (size_t)q * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vd a0 = V_LDD(p0 + o), a1 = V_LDD(p1 + o);
                        rns_vd a2 = V_LDD(p2 + o), a3 = V_LDD(p3 + o);
                        VPD_FWD2(a0, a2, wa, wap, 0);
                        VPD_FWD2(a1, a3, wb, wbp, 0);
                        VPD_FWD2(a0, a1, wc, wcp, 1);
                        VPD_FWD2(a2, a3, wc, wcp, 1);
                        V_STD(p0 + o, a0);
                        V_STD(p1 + o, a1);
                        V_STD(p2 + o, a2);
                        V_STD(p3 + o, a3);
                    }
                }
        }
    }
}

void vpk_inv_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L, int i0, int i1)
{
    int i, h;
    uint32_t m, b, j, st;

    for (i = i0; i < i1; i++)
    {
        rns_vd p = V_BCD((double)ctx->p[i]), pinv = V_BCD(ctx->pinvd[i]);
        double *base = (double *)(v + (size_t)(i - i0) * L * 8);
        const double *tw = (const double *)ctx->itw + (size_t)i * ctx->Lmax;
        for (m = 1; 2 * m < L; m <<= 2)
        {
            // levels m and 2m on a0..a3 = b + j + {0, m, 2m, 3m}
            st = ctx->Lmax / (4 * m);
            for (b = 0; b < L; b += 4 * m)
                for (j = 0; j < m; j++)
                {
                    rns_vd wa = V_BCD(tw[4 * j * st]), wap = V_BCD(tw[4 * j * st + 1]);
                    rns_vd wb = V_BCD(tw[2 * j * st]), wbp = V_BCD(tw[2 * j * st + 1]);
                    rns_vd wc = V_BCD(tw[2 * (j + m) * st]), wcp = V_BCD(tw[2 * (j + m) * st + 1]);
                    double *p0 = base + (size_t)(b + j) * 8, *p1 = p0 + (size_t)m * 8;
                    double *p2 = p1 + (size_t)m * 8, *p3 = p2 + (size_t)m * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vd a0 = V_LDD(p0 + o), a1 = V_LDD(p1 + o);
                        rns_vd a2 = V_LDD(p2 + o), a3 = V_LDD(p3 + o);
                        VPD_INV2(a0, a1, wa, wap, 0);
                        VPD_INV2(a2, a3, wa, wap, 0);
                        VPD_INV2(a0, a2, wb, wbp, 1);
                        VPD_INV2(a1, a3, wc, wcp, 1);
                        V_STD(p0 + o, a0);
                        V_STD(p1 + o, a1);
                        V_STD(p2 + o, a2);
                        V_STD(p3 + o, a3);
                    }
                }
        }
        if (m < L)
        {
            // last radix-2 level (m = L/2)
            st = ctx->Lmax / (2 * m);
            for (j = 0; j < m; j++)
            {
                rns_vd w = V_BCD(tw[2 * j * st]), wp = V_BCD(tw[2 * j * st + 1]);
                double *px = base + (size_t)j * 8, *py = px + (size_t)m * 8;
                for (h = 0; h < VP_H; h++)
                {
                    int o = h * V_LANES;
                    rns_vd x = V_LDD(px + o), y = V_LDD(py + o);
                    VPD_INV2(x, y, w, wp, 1);
                    V_STD(px + o, x);
                    V_STD(py + o, y);
                }
            }
        }
    }
}

// c = a b mod p (no Montgomery factor)
void vpk_pmul_r(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t L, int i0, int i1)
{
    int i, h;
    uint32_t t;

    for (i = i0; i < i1; i++)
    {
        rns_vd p = V_BCD((double)ctx->p[i]), pinv = V_BCD(ctx->pinvd[i]);
        size_t off = (size_t)(i - i0) * L * 8;
        for (t = 0; t < L; t++)
        {
            for (h = 0; h < VP_H; h++)
            {
                size_t o = off + (size_t)t * 8 + (size_t)(h * V_LANES);
                V_STD((double *)c + o, vpd_mulm(V_LDD(a + o), V_LDD(b + o), p, pinv));
            }
        }
    }
}

#else

/*
 * Integer transforms, radix-4 passes as above.  Ranges: DIF keeps [0, 2p)
 * (sum: conditional subtract of 2p; difference + 2p < 4p into the Shoup
 * product); DIT takes [0, 4p) to [0, 4p) (x reduced below 2p, t = y w in
 * [0, 2p)).  The same butterflies twice in registers keep the same ranges.
 */

#define VPI_FWD2(x, y, w, wp) do {                                      \
        rns_vi s_ = V_CSUB(V_ADD((x), (y)), p2);                        \
        (y) = vp_shoup(V_ADD(V_SUB((x), (y)), p2), (w), (wp), p);       \
        (x) = s_;                                                       \
    } while (0)

#define VPI_INV2(x, y, w, wp) do {                                      \
        rns_vi x_ = V_CSUB((x), p2);                                    \
        rns_vi t_ = vp_shoup((y), (w), (wp), p);                        \
        (x) = V_ADD(x_, t_);                                            \
        (y) = V_ADD(V_SUB(x_, t_), p2);                                 \
    } while (0)

static int vp_log2(uint32_t L)
{
    int l = 0;
    while ((1u << l) < L) l++;
    return l;
}

void vpk_fwd_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L, int i0, int i1)
{
    int i, h, lg = vp_log2(L);
    uint32_t m, b, j, st;
    VP_RZ_ENTER();

    for (i = i0; i < i1; i++)
    {
        rns_vi p = V_BC(ctx->p[i]);
        rns_vi p2 = V_BC(2 * ctx->p[i]);
        uint64_t *base = v + (size_t)(i - i0) * L * 8;
        const uint64_t *tw = ctx->tw + (size_t)i * ctx->Lmax;
        m = L >> 1;
        if (lg & 1)
        {
            st = ctx->Lmax / (2 * m);
            for (b = 0; b < L; b += 2 * m)
                for (j = 0; j < m; j++)
                {
                    rns_vi w = V_BC(tw[2 * j * st]), wp = V_BC(tw[2 * j * st + 1]);
                    uint64_t *px = base + (size_t)(b + j) * 8, *py = px + (size_t)m * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vi x = V_LD(px + o), y = V_LD(py + o);
                        VPI_FWD2(x, y, w, wp);
                        V_ST(px + o, x);
                        V_ST(py + o, y);
                    }
                }
            m >>= 1;
        }
        for (; m >= 2; m >>= 2)
        {
            uint32_t q = m / 2;
            st = ctx->Lmax / (2 * m);
            for (b = 0; b < L; b += 2 * m)
                for (j = 0; j < q; j++)
                {
                    rns_vi wa = V_BC(tw[2 * j * st]), wap = V_BC(tw[2 * j * st + 1]);
                    rns_vi wb = V_BC(tw[2 * (j + q) * st]), wbp = V_BC(tw[2 * (j + q) * st + 1]);
                    rns_vi wc = V_BC(tw[4 * j * st]), wcp = V_BC(tw[4 * j * st + 1]);
                    uint64_t *p0 = base + (size_t)(b + j) * 8, *p1 = p0 + (size_t)q * 8;
                    uint64_t *pp2 = p0 + (size_t)m * 8, *p3 = pp2 + (size_t)q * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vi a0 = V_LD(p0 + o), a1 = V_LD(p1 + o);
                        rns_vi a2 = V_LD(pp2 + o), a3 = V_LD(p3 + o);
                        VPI_FWD2(a0, a2, wa, wap);
                        VPI_FWD2(a1, a3, wb, wbp);
                        VPI_FWD2(a0, a1, wc, wcp);
                        VPI_FWD2(a2, a3, wc, wcp);
                        V_ST(p0 + o, a0);
                        V_ST(p1 + o, a1);
                        V_ST(pp2 + o, a2);
                        V_ST(p3 + o, a3);
                    }
                }
        }
    }
    VP_RZ_LEAVE();
}

void vpk_inv_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L, int i0, int i1)
{
    int i, h;
    uint32_t m, b, j, st;
    VP_RZ_ENTER();

    for (i = i0; i < i1; i++)
    {
        rns_vi p = V_BC(ctx->p[i]);
        rns_vi p2 = V_BC(2 * ctx->p[i]);
        uint64_t *base = v + (size_t)(i - i0) * L * 8;
        const uint64_t *tw = ctx->itw + (size_t)i * ctx->Lmax;
        for (m = 1; 2 * m < L; m <<= 2)
        {
            st = ctx->Lmax / (4 * m);
            for (b = 0; b < L; b += 4 * m)
                for (j = 0; j < m; j++)
                {
                    rns_vi wa = V_BC(tw[4 * j * st]), wap = V_BC(tw[4 * j * st + 1]);
                    rns_vi wb = V_BC(tw[2 * j * st]), wbp = V_BC(tw[2 * j * st + 1]);
                    rns_vi wc = V_BC(tw[2 * (j + m) * st]), wcp = V_BC(tw[2 * (j + m) * st + 1]);
                    uint64_t *p0 = base + (size_t)(b + j) * 8, *p1 = p0 + (size_t)m * 8;
                    uint64_t *pp2 = p1 + (size_t)m * 8, *p3 = pp2 + (size_t)m * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vi a0 = V_LD(p0 + o), a1 = V_LD(p1 + o);
                        rns_vi a2 = V_LD(pp2 + o), a3 = V_LD(p3 + o);
                        VPI_INV2(a0, a1, wa, wap);
                        VPI_INV2(a2, a3, wa, wap);
                        VPI_INV2(a0, a2, wb, wbp);
                        VPI_INV2(a1, a3, wc, wcp);
                        V_ST(p0 + o, a0);
                        V_ST(p1 + o, a1);
                        V_ST(pp2 + o, a2);
                        V_ST(p3 + o, a3);
                    }
                }
        }
        if (m < L)
        {
            st = ctx->Lmax / (2 * m);
            for (j = 0; j < m; j++)
            {
                rns_vi w = V_BC(tw[2 * j * st]), wp = V_BC(tw[2 * j * st + 1]);
                uint64_t *px = base + (size_t)j * 8, *py = px + (size_t)m * 8;
                for (h = 0; h < VP_H; h++)
                {
                    int o = h * V_LANES;
                    rns_vi x = V_LD(px + o), y = V_LD(py + o);
                    VPI_INV2(x, y, w, wp);
                    V_ST(px + o, x);
                    V_ST(py + o, y);
                }
            }
        }
    }
    VP_RZ_LEAVE();
}

void vpk_pmul_r(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t L, int i0, int i1)
{
    int i, h;
    uint32_t t;
    VP_RZ_ENTER();

    for (i = i0; i < i1; i++)
    {
        rns_vi p = V_BC(ctx->p[i]), pneg = V_BC(ctx->pneg[i]);
        size_t off = (size_t)(i - i0) * L * 8;
        for (t = 0; t < L; t++)
        {
            for (h = 0; h < VP_H; h++)
            {
                size_t o = off + (size_t)t * 8 + (size_t)(h * V_LANES);
                V_ST(c + o, vp_mont(V_LD(a + o), V_LD(b + o), p, pneg));
            }
        }
    }
    VP_RZ_LEAVE();
}
#endif  // !VP_FPNTT

void vpk_from_ntt(const vp_ctx_t *ctx, uint64_t *a, const uint64_t *v,
    uint32_t L, uint32_t start, uint32_t count)
{
    int n = ctx->n, K = ctx->K, i, k, h, l = 0;
    uint32_t t;
    rns_vi col[VP_MAXW + 2], hic[VP_MAXW + 2];
    const uint64_t *ys;
    VP_RZ_ENTER();

    while ((1u << l) < L) l++;
    ys = ctx->yscale + (size_t)l * (size_t)K * 2;

    for (t = 0; t < count; t++)
    {
        uint64_t *c = a + (size_t)t * 8 * (size_t)n;
        for (h = 0; h < VP_H; h++)
        {
            int o = h * V_LANES;
            rns_vd sum = V_ZEROD;
            rns_vi kap, m, carry;
            vp_op yo;

            for (k = 0; k <= n + 1; k++)
            {
                col[k] = V_ZERO;
                hic[k] = V_ZERO;
            }
            for (i = 0; i < K; i++)
            {
                const uint64_t *cr = ctx->crt + (size_t)i * (size_t)n;
                const double *crd = ctx->crtd + (size_t)i * (size_t)n;
#if VP_FPNTT
                // y = x * scale in [0, p), an exact double (= its VP_OP form)
                const double *ysd = (const double *)ys;
                rns_vd pd = V_BCD((double)ctx->p[i]);
                rns_vd y = D_ADDNEG(vpd_mulc(V_LDD(VP_RES(v, L, i, start + t) + o),
                    V_BCD(ysd[2 * i]), V_BCD(ysd[2 * i + 1]), pd), pd);
                sum = V_FMAD(y, V_BCD(ctx->pinvd[i]), sum);
                yo = y;
#else
                rns_vi y = vp_shoup(V_LD(VP_RES(v, L, i, start + t) + o),
                    V_BC(ys[2 * i]), V_BC(ys[2 * i + 1]), V_BC(ctx->p[i]));
                sum = V_FMAD(V_CVTD(y), V_BCD(ctx->pinvd[i]), sum);
                yo = VP_OP(y);
#endif
                for (k = 0; k < n; k++)
                    VP_MAC(col[k], hic[k + 1], yo, VP_OPT(cr, crd, k));
            }
            kap = V_KAPPA(sum);
            yo = VP_OP(kap);
            {
                const uint64_t *cr = ctx->crt + (size_t)K * (size_t)n;
                const double *crd = ctx->crtd + (size_t)K * (size_t)n;
                for (k = 0; k < n; k++)
                    VP_MAC(col[k], hic[k + 1], yo, VP_OPT(cr, crd, k));
            }
            for (k = 0; k < n; k++)
            {
                VP_UNBIAS_LO(col[k], K + 1);
                VP_UNBIAS_HI(hic[k + 1], K + 1);
            }
            for (k = 1; k <= n; k++)
                col[k] = V_ADD(col[k], hic[k]);

            // two Montgomery word reductions (divide by 2^104)
            m = LO52(V_ZERO, col[0], V_BC(ctx->n0));
            for (k = 0; k < n; k++)
            {
                rns_vi nw = V_BC(ctx->Nw[k]);
                col[k] = LO52(col[k], m, nw);
                col[k + 1] = HI52(col[k + 1], m, nw);
            }
            col[1] = V_ADD(col[1], V_SRLI(col[0], 52));
            m = LO52(V_ZERO, col[1], V_BC(ctx->n0));
            for (k = 0; k < n; k++)
            {
                rns_vi nw = V_BC(ctx->Nw[k]);
                col[k + 1] = LO52(col[k + 1], m, nw);
                col[k + 2] = HI52(col[k + 2], m, nw);
            }
            carry = V_SRLI(col[1], 52);
            for (k = 0; k < n; k++)
            {
                rns_vi x = V_ADD(col[k + 2], carry);
                V_ST(c + k * 8 + o, V_AND(x, VM52));
                carry = V_SRLI(x, 52);
            }
        }
    }
    VP_RZ_LEAVE();
}

void vpk_add(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t cnt)
{
    int n = ctx->n, k, h;
    uint32_t t;
    rns_vi s[VP_MAXW];

    for (t = 0; t < cnt; t++)
    {
        size_t base = (size_t)t * 8 * (size_t)n;
        for (h = 0; h < VP_H; h++)
        {
            int o = h * V_LANES;
            rns_vi cy = V_ZERO, bw = V_ZERO, sel;
            for (k = 0; k < n; k++)
            {
                rns_vi x = V_ADD(V_ADD(V_LD(a + base + k * 8 + o),
                    V_LD(b + base + k * 8 + o)), cy);
                cy = V_SRLI(x, 52);
                s[k] = V_AND(x, VM52);
            }
            for (k = 0; k < n; k++)
            {
                rns_vi x = V_SUB(V_SUB(s[k], V_BC(ctx->N2w[k])), bw);
                bw = V_SRLI(x, 63);
                V_ST(c + base + k * 8 + o, V_AND(x, VM52));
            }
            // borrow: sum < 2N, keep the sum
            sel = V_SUB(V_ZERO, bw);
            for (k = 0; k < n; k++)
            {
                rns_vi d = V_LD(c + base + k * 8 + o);
                V_ST(c + base + k * 8 + o, V_BLEND_SIGN(d, s[k], sel));
            }
        }
    }
}

void vpk_sub(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t cnt)
{
    int n = ctx->n, k, h;
    uint32_t t;
    rns_vi s[VP_MAXW];

    for (t = 0; t < cnt; t++)
    {
        size_t base = (size_t)t * 8 * (size_t)n;
        for (h = 0; h < VP_H; h++)
        {
            int o = h * V_LANES;
            rns_vi cy = V_ZERO, bw = V_ZERO, sel;
            for (k = 0; k < n; k++)
            {
                rns_vi x = V_SUB(V_SUB(V_LD(a + base + k * 8 + o),
                    V_LD(b + base + k * 8 + o)), bw);
                bw = V_SRLI(x, 63);
                s[k] = V_AND(x, VM52);
            }
            // borrow: a < b, add 2N
            sel = V_SUB(V_ZERO, bw);
            for (k = 0; k < n; k++)
            {
                rns_vi x = V_ADD(V_ADD(s[k], V_BC(ctx->N2w[k])), cy);
                cy = V_SRLI(x, 52);
                V_ST(c + base + k * 8 + o, V_BLEND_SIGN(s[k], V_AND(x, VM52), sel));
            }
        }
    }
}

void vpk_neg(const vp_ctx_t *ctx, uint64_t *c, const uint64_t *a, uint32_t cnt)
{
    // 0 - a: a scratch zero coefficient per call
    int n = ctx->n, k, h;
    uint32_t t;
    rns_vi s[VP_MAXW];

    for (t = 0; t < cnt; t++)
    {
        size_t base = (size_t)t * 8 * (size_t)n;
        for (h = 0; h < VP_H; h++)
        {
            int o = h * V_LANES;
            rns_vi cy = V_ZERO, bw = V_ZERO, sel;
            for (k = 0; k < n; k++)
            {
                rns_vi x = V_SUB(V_SUB(V_ZERO, V_LD(a + base + k * 8 + o)), bw);
                bw = V_SRLI(x, 63);
                s[k] = V_AND(x, VM52);
            }
            sel = V_SUB(V_ZERO, bw);
            for (k = 0; k < n; k++)
            {
                rns_vi x = V_ADD(V_ADD(s[k], V_BC(ctx->N2w[k])), cy);
                cy = V_SRLI(x, 52);
                V_ST(c + base + k * 8 + o, V_BLEND_SIGN(s[k], V_AND(x, VM52), sel));
            }
        }
    }
}

// 1: the NTT domain holds doubles (see VP_FPNTT); vp_ctx_init builds the
// twiddle and CRT-scale tables to match
int vpk_fp_ntt(void)
{
    return VP_FPNTT;
}

#else

int vpk_fp_ntt(void)
{
    return 0;
}

#endif
