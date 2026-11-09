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
RNS_INLINE __m256i vp_csub256(__m256i x, __m256i m)
{
    __m256i t = _mm256_sub_epi64(x, m);
    return _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(t),
        _mm256_castsi256_pd(x), _mm256_castsi256_pd(t)));
}
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
#define VP_MAC(l, h, x, y)  RNS_FP_MAC(l, h, x, y)
#define VP_UNBIAS_LO(l, cnt) ((l) = rns_fp_unbias((l), (cnt), 0))
#define VP_UNBIAS_HI(h, cnt) ((h) = rns_fp_unbias((h), 0, (cnt)))
#else
typedef rns_vi vp_op;
#define VP_OP(x)            (x)
#define VP_OPC(c)           V_BC(c)
#define VP_MAC(l, h, x, y)  do { (l) = LO52((l), (x), (y)); (h) = HI52((h), (x), (y)); } while (0)
#define VP_UNBIAS_LO(l, cnt)
#define VP_UNBIAS_HI(h, cnt)
#endif

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
                rns_vi p = V_BC(ctx->p[i]), pneg = V_BC(ctx->pneg[i]);
                rns_vi lo = V_ZERO, hi = V_ZERO, m, r, h0;
                for (k = 0; k < n; k++)
                    VP_MAC(lo, hi, d[k], VP_OPC(cv[k]));
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
                V_ST(VP_RES(v, L, i - i0, t) + o, r);
            }
        }
    }
    for (i = i0; i < i1; i++)
        for (t = len; t < L; t++)
            for (h = 0; h < VP_H; h++)
                V_ST(VP_RES(v, L, i - i0, t) + h * V_LANES, V_ZERO);
    VP_RZ_LEAVE();
}

void vpk_fwd_r(const vp_ctx_t *ctx, uint64_t *v, uint32_t L, int i0, int i1)
{
    int i, h;
    uint32_t m, b, j, st;
    VP_RZ_ENTER();

    for (i = i0; i < i1; i++)
    {
        rns_vi p = V_BC(ctx->p[i]);
        rns_vi p2 = V_BC(2 * ctx->p[i]);
        uint64_t *base = v + (size_t)(i - i0) * L * 8;
        for (m = L >> 1; m >= 1; m >>= 1)
        {
            // w_2m^j = w_Lmax^(j st), st = Lmax / 2m
            const uint64_t *tw = ctx->tw + (size_t)i * ctx->Lmax;
            st = ctx->Lmax / (2 * m);
            for (b = 0; b < L; b += 2 * m)
            {
                for (j = 0; j < m; j++)
                {
                    rns_vi w = V_BC(tw[2 * j * st]), wp = V_BC(tw[2 * j * st + 1]);
                    uint64_t *px = base + (size_t)(b + j) * 8;
                    uint64_t *py = px + (size_t)m * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vi x = V_LD(px + o), y = V_LD(py + o);
                        rns_vi s = V_CSUB(V_ADD(x, y), p2);
                        rns_vi d = V_ADD(V_SUB(x, y), p2);
                        V_ST(px + o, s);
                        V_ST(py + o, vp_shoup(d, w, wp, p));
                    }
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
        for (m = 1; m < L; m <<= 1)
        {
            const uint64_t *tw = ctx->itw + (size_t)i * ctx->Lmax;
            st = ctx->Lmax / (2 * m);
            for (b = 0; b < L; b += 2 * m)
            {
                for (j = 0; j < m; j++)
                {
                    rns_vi w = V_BC(tw[2 * j * st]), wp = V_BC(tw[2 * j * st + 1]);
                    uint64_t *px = base + (size_t)(b + j) * 8;
                    uint64_t *py = px + (size_t)m * 8;
                    for (h = 0; h < VP_H; h++)
                    {
                        int o = h * V_LANES;
                        rns_vi x = V_CSUB(V_LD(px + o), p2);
                        rns_vi tt = vp_shoup(V_LD(py + o), w, wp, p);
                        V_ST(px + o, V_ADD(x, tt));
                        V_ST(py + o, V_ADD(V_SUB(x, tt), p2));
                    }
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
                rns_vi y = vp_shoup(V_LD(VP_RES(v, L, i, start + t) + o),
                    V_BC(ys[2 * i]), V_BC(ys[2 * i + 1]), V_BC(ctx->p[i]));
                sum = V_FMAD(V_CVTD(y), V_BCD(ctx->pinvd[i]), sum);
                yo = VP_OP(y);
                for (k = 0; k < n; k++)
                    VP_MAC(col[k], hic[k + 1], yo, VP_OPC(cr[k]));
            }
            kap = V_KAPPA(sum);
            yo = VP_OP(kap);
            {
                const uint64_t *cr = ctx->crt + (size_t)K * (size_t)n;
                for (k = 0; k < n; k++)
                    VP_MAC(col[k], hic[k + 1], yo, VP_OPC(cr[k]));
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

#endif
