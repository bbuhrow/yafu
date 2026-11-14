/*
microecm_edw_x4.c: microecm_edw's vector list on AVX2, 4 inputs n < 2^52 at
a time (one per 64-bit lane).  Same curves, chains, stage 2 (B2 = 25 B1,
a lane's accumulator skips products that would be 0 mod n) and results as
the AVX-512 x8 list in microecm_edw.c; the curve code is edw_kern.h.

Arithmetic: residues are doubles holding integers in [0, n), R = 2^52.
The Montgomery product uses FMA with the default rounding (to nearest), so
no MXCSR change is needed:
  t = fma(x, y, 2^104),  hd = t - 2^104,  l = fma(x, y, -hd)
gives x y = hd + l exactly with hd a multiple of 2^52 and |l| <= 2^51;
the same with a bias of 3 2^103 splits signed products.  With
rho = -1/n mod 2^52, m = (l rho) mod 2^52 (signed, |m| <= 2^51) makes
x y + m n divisible by 2^52, and r = (x y + m n) / 2^52 lies in
[-n/2, 3n/2), which two sign-bit blends bring into [0, n).  Every other
FP operation is exact.  Needs FMA (built with -mfma in the Makefile).

Design and measurements: arith/micro_tiny_edwards_plan.md.

Dual-licensed, at your option, under the FreeBSD license or the Mozilla
Public License 2.0, as microecm.c.
*/

#include "microecm_edw.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifndef LIMB1_SCALAR_ONLY
#define LIMB1_SCALAR_ONLY
#endif
#include "limb1.h"          /* shared 64-bit kernels (arith/limb) */

#ifndef ALIGNED_MEM
#if defined(_MSC_VER)
#define ALIGNED_MEM __declspec(align(64))
#else
#define ALIGNED_MEM __attribute__((aligned(64)))
#endif
#endif

#if defined(_MSC_VER)
#define U4_INLINE __forceinline
#elif defined(__GNUC__)
#define U4_INLINE inline __attribute__((always_inline))
#else
#define U4_INLINE __inline
#endif

#ifdef USE_AVX2
#include <immintrin.h>
#include "edw_tiny_curves.h"
#include "edw_tiny_v2.h"

/* microecm_edw.c: the parameter set and its v2 levels */
const edw2_level_t *uedw_level2_pub(uint32_t B1);

#define U4_PD(bits) _mm256_castsi256_pd(_mm256_set1_epi64x((long long)(bits)))
#define U4_B104 0x4670000000000000ULL       // 2^104
#define U4_B3_103 0x4678000000000000ULL     // 3 2^103
#define U4_SC 0x3CB0000000000000ULL         // 2^-52

/* Montgomery product, inputs integers in [0, 2^52) with x y < 2^52 n.
   x y = hd + l1 (hd = H1 2^52, |l1| <= 2^51), m n = vd + l2 likewise, and
   l1 + l2 = 0 mod 2^52, so r = H1 + H2 + (l1 + l2) / 2^52.  This full form
   handles every case; uedw4_mulredc uses it only when |l1| = 2^51. */
#if defined(_MSC_VER)
static __declspec(noinline) __m256d uedw4_mulredc_full(__m256d x, __m256d y, __m256d n, __m256d rho)
#else
static __attribute__((noinline)) __m256d uedw4_mulredc_full(__m256d x, __m256d y, __m256d n, __m256d rho)
#endif
{
    const __m256d b1 = U4_PD(U4_B104);
    const __m256d b2 = U4_PD(U4_B3_103);
    __m256d t = _mm256_fmadd_pd(x, y, b1);
    __m256d hd = _mm256_sub_pd(t, b1);
    __m256d l1 = _mm256_fmsub_pd(x, y, hd);
    __m256d u = _mm256_fmadd_pd(l1, rho, b2);
    __m256d m = _mm256_fmsub_pd(l1, rho, _mm256_sub_pd(u, b2));
    __m256d v = _mm256_fmadd_pd(m, n, b2);
    __m256d vd = _mm256_sub_pd(v, b2);
    __m256d l2 = _mm256_fmsub_pd(m, n, vd);
    __m256d r = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(hd, vd), _mm256_add_pd(l1, l2)),
        U4_PD(U4_SC));
    __m256d rm = _mm256_sub_pd(r, n);
    r = _mm256_blendv_pd(rm, r, rm);
    return _mm256_blendv_pd(r, _mm256_add_pd(r, n), r);
}

/* The fast form: when |l1| < 2^51, l2 = -l1 exactly, so r = H1 + H2 comes
   straight from the high halves, r = t / 2^52 - 2^52 + v / 2^52 - 3 2^51,
   with r - n and r + n formed in parallel (latency ~ 9 FP ops instead of
   ~15).  A tie |l1| = 2^51 (probability ~ 2^-52 per lane) takes the full
   form. */
static U4_INLINE __m256d uedw4_mulredc(__m256d x, __m256d y, __m256d n, __m256d rho)
{
    const __m256d b1 = U4_PD(U4_B104);
    const __m256d b2 = U4_PD(U4_B3_103);
    const __m256d sc = U4_PD(U4_SC);
    __m256d t = _mm256_fmadd_pd(x, y, b1);
    __m256d hd = _mm256_sub_pd(t, b1);
    __m256d l1 = _mm256_fmsub_pd(x, y, hd);
    __m256d k1 = _mm256_fmadd_pd(t, sc, U4_PD(0xC344000000000000ULL));    // H1 - 3 2^51
    __m256d u = _mm256_fmadd_pd(l1, rho, b2);
    __m256d m = _mm256_fmsub_pd(l1, rho, _mm256_sub_pd(u, b2));
    __m256d v = _mm256_fmadd_pd(m, n, b2);
    __m256d r = _mm256_fmadd_pd(v, sc, k1);
    __m256d rm = _mm256_fmadd_pd(v, sc, _mm256_sub_pd(k1, n));
    __m256d rp = _mm256_fmadd_pd(v, sc, _mm256_add_pd(k1, n));
    __m256d tie = _mm256_cmp_pd(_mm256_and_pd(l1, U4_PD(0x7fffffffffffffffULL)),
        U4_PD(0x4320000000000000ULL), _CMP_EQ_OQ);
    if (_mm256_movemask_pd(tie))
        return uedw4_mulredc_full(x, y, n, rho);
    __m256d s = _mm256_blendv_pd(rm, r, rm);
    return _mm256_blendv_pd(s, rp, r);
}

static U4_INLINE __m256d uedw4_addmod(__m256d x, __m256d y, __m256d n)
{
    __m256d s = _mm256_add_pd(x, y);
    __m256d sm = _mm256_sub_pd(s, n);
    return _mm256_blendv_pd(sm, s, sm);
}

static U4_INLINE __m256d uedw4_submod(__m256d x, __m256d y, __m256d n)
{
    __m256d d = _mm256_sub_pd(x, y);
    return _mm256_blendv_pd(d, _mm256_add_pd(d, n), d);
}

/* A residue vector is UEDW4_W ymm registers (4 lanes each): with W = 2 the
   two halves are independent work for the out-of-order core, which hides
   the long latency of the FP product (the Edwards formulas have only 2-4
   independent products per step).  Measured: W = 1 is faster (W = 2
   spills), so it is the default. */
#ifndef UEDW4_W
#define UEDW4_W 1
#endif
#define U4_VL (4 * UEDW4_W)
#define U4_ALL ((1u << U4_VL) - 1)

typedef struct { __m256d h[UEDW4_W]; } u4v;

static U4_INLINE u4v u4_mul(u4v x, u4v y, const u4v *n, const u4v *rho)
{
    u4v r;
    int i;
    for (i = 0; i < UEDW4_W; i++)
        r.h[i] = uedw4_mulredc(x.h[i], y.h[i], n->h[i], rho->h[i]);
    return r;
}

static U4_INLINE u4v u4_add(u4v x, u4v y, const u4v *n)
{
    u4v r;
    int i;
    for (i = 0; i < UEDW4_W; i++)
        r.h[i] = uedw4_addmod(x.h[i], y.h[i], n->h[i]);
    return r;
}

static U4_INLINE u4v u4_sub(u4v x, u4v y, const u4v *n)
{
    u4v r;
    int i;
    for (i = 0; i < UEDW4_W; i++)
        r.h[i] = uedw4_submod(x.h[i], y.h[i], n->h[i]);
    return r;
}

static U4_INLINE u4v u4_zero(void)
{
    u4v r;
    int i;
    for (i = 0; i < UEDW4_W; i++)
        r.h[i] = _mm256_setzero_pd();
    return r;
}

/* acc * t, except in lanes where that is 0 mod n (microecm x8's rule) */
static U4_INLINE u4v u4_acc(u4v acc, u4v t, const u4v *n, const u4v *rho)
{
    u4v p = u4_mul(acc, t, n, rho);
    int i;
    for (i = 0; i < UEDW4_W; i++)
        acc.h[i] = _mm256_blendv_pd(p.h[i], acc.h[i],
            _mm256_cmp_pd(p.h[i], _mm256_setzero_pd(), _CMP_EQ_OQ));
    return acc;
}

/* all-ones in the lanes of the U4_VL-bit mask m */
static U4_INLINE u4v u4_lanes(unsigned m)
{
    const __m256i bits = _mm256_set_epi64x(8, 4, 2, 1);
    u4v r;
    int i;
    for (i = 0; i < UEDW4_W; i++) {
        __m256i v = _mm256_and_si256(_mm256_set1_epi64x((m >> (4 * i)) & 0xf), bits);
        r.h[i] = _mm256_castsi256_pd(_mm256_cmpeq_epi64(v, bits));
    }
    return r;
}

/* b in the lanes of mask, else a */
static U4_INLINE u4v u4_sel(u4v mask, u4v a, u4v b)
{
    int i;
    for (i = 0; i < UEDW4_W; i++)
        a.h[i] = _mm256_blendv_pd(a.h[i], b.h[i], mask.h[i]);
    return a;
}

static U4_INLINE u4v u4_load(const double *a)
{
    u4v r;
    int i;
    for (i = 0; i < UEDW4_W; i++)
        r.h[i] = _mm256_load_pd(a + 4 * i);
    return r;
}

static U4_INLINE void u4_store(double *a, u4v v)
{
    int i;
    for (i = 0; i < UEDW4_W; i++)
        _mm256_store_pd(a + 4 * i, v.h[i]);
}

typedef struct {
    u4v n, rho, one, Rsqr, R3;
    u4v c6, c10, c16, c32, c48;
} uedw_ctx4;

#define EK_T u4v
#define EK_CTX uedw_ctx4
#define EK_FN(x) uedw4_##x
#define EK_MUL(a, b) u4_mul((a), (b), &c->n, &c->rho)
#define EK_SQR(a) u4_mul((a), (a), &c->n, &c->rho)
#define EK_ADD(a, b) u4_add((a), (b), &c->n)
#define EK_SUB(a, b) u4_sub((a), (b), &c->n)
#define EK_ZERO u4_zero()
#define EK_ACC_DECL
#define EK_ACC(acc, t) { acc = u4_acc(acc, t, &c->n, &c->rho); }
#include "edw_kern.h"
#undef EK_T
#undef EK_CTX
#undef EK_FN
#undef EK_MUL
#undef EK_SQR
#undef EK_ADD
#undef EK_SUB
#undef EK_ZERO
#undef EK_ACC_DECL
#undef EK_ACC

#define M4(a, b) u4_mul((a), (b), &c->n, &c->rho)
#define A4(a, b) u4_add((a), (b), &c->n)
#define D4(a, b) u4_sub((a), (b), &c->n)

#define uedw4_inverse multiplicative_inverse   /* arith/limb/limb1.h */

static uint64_t uedw4_gcd(uint64_t u, uint64_t v)
{
    while (v) {
        uint64_t t = u % v;
        u = v;
        v = t;
    }
    return u;
}

/* the derived constants, from n, rho and one */
static void uedw4_ctx_finish(uedw_ctx4 *c)
{
    u4v c2 = A4(c->one, c->one);
    u4v c4 = A4(c2, c2);
    u4v c8 = A4(c4, c4);
    u4v r;
    int i;
    c->c16 = A4(c8, c8);
    c->c32 = A4(c->c16, c->c16);
    c->c48 = A4(c->c32, c->c16);
    c->c6 = A4(c4, c2);
    c->c10 = A4(c8, c2);
    r = c->c16;                                 // 2^4 R; R^2 = 2^52 R
    for (i = 4; i < 52; i++)
        r = A4(r, r);
    c->Rsqr = r;
    c->R3 = M4(r, r);
}

/* a signed constant per lane (|v| < 2^63) to Montgomery form */
static u4v uedw4_const(const int64_t *v, const uedw_ctx4 *c)
{
    ALIGNED_MEM double lo[U4_VL], hi[U4_VL];
    int i;
    unsigned neg = 0;
    for (i = 0; i < U4_VL; i++) {
        uint64_t a = (v[i] < 0) ? (uint64_t)0 - (uint64_t)v[i] : (uint64_t)v[i];
        lo[i] = (double)(a & 0x000fffffffffffffULL);
        hi[i] = (double)(a >> 52);
        if (v[i] < 0) neg |= 1u << i;
    }
    u4v r = A4(M4(u4_load(lo), c->Rsqr), M4(u4_load(hi), c->R3));
    return u4_sel(u4_lanes(neg), r, D4(u4_zero(), r));
}

static void uedw4_table_curve(const int *k, uedw4_curve *cv, const uedw_ctx4 *c)
{
    int64_t v[U4_VL];
    u4v w[6];
    int i, j;
    for (j = 0; j < 6; j++) {
        for (i = 0; i < U4_VL; i++)
            v[i] = edw_tiny_curves[k[i]][j];
        w[j] = uedw4_const(v, c);
    }
    uedw4_curve_from(w, cv, c);
}

static void uedw4_blend_curve(u4v m, uedw4_curve *cv, const uedw4_curve *s)
{
    cv->P.X = u4_sel(m, cv->P.X, s->P.X);
    cv->P.Y = u4_sel(m, cv->P.Y, s->P.Y);
    cv->P.Z = u4_sel(m, cv->P.Z, s->P.Z);
    cv->P.T = u4_sel(m, cv->P.T, s->P.T);
    cv->dn2 = u4_sel(m, cv->dn2, s->dn2);
    cv->dd = u4_sel(m, cv->dd, s->dd);
    cv->dd2 = u4_sel(m, cv->dd2, s->dd2);
    cv->a24n = u4_sel(m, cv->a24n, s->a24n);
    cv->a24d = u4_sel(m, cv->a24d, s->a24d);
}

static void uedw4_blend_w(u4v m, uedw4_wpt *R, const uedw4_wpt *S)
{
    R->X = u4_sel(m, R->X, S->X);
    R->Y = u4_sel(m, R->Y, S->Y);
    R->Z = u4_sel(m, R->Z, S->Z);
}

/* ECM on n[0..num-1] (all < 2^52), U4_VL at a time; f[i] = factor or 1.
   *ncurves (if not 0) gets the total number of curves run. */
static void uedw4_list(const uint64_t *n64, uint64_t *f, uint32_t num, uint32_t B1,
    uint32_t B2, uint32_t curves, uint64_t *ncurves)
{
    ALIGNED_MEM double nA[U4_VL], rhoA[U4_VL], oneA[U4_VL], st[U4_VL];
    uint64_t nU[U4_VL];
    int kcur[U4_VL], idx[U4_VL], kk[U4_VL];
    unsigned msk = 0;
    uint32_t j = 0, i;
    int lane;
    edw_chain_t ch;
    const uint8_t *barray;
    int numb;
    uedw_ctx4 ctx;
    const uedw_ctx4 *c = &ctx;
    uedw4_wpt R;
    uint64_t total = 0;
    const int ntab = EDW_TINY_NCURVES;

    for (i = 0; i < num; i++)
        f[i] = 1;
    const edw2_level_t *lv = 0;
    if (num == 0)
        return;
    if (uecm_edw_get_params() == 2) {
        if (!(lv = uedw_level2_pub(B1)))
            return;
        ch.ch = lv->ch;
        ch.len = lv->len;
        ch.m = lv->m;
        ch.top = lv->top;
        ch.tz = lv->tz;
    } else if (!edw_chain_for(B1, &ch) || !uedw_s2table25(B1, &barray, &numb)) {
        return;
    }
    // every lane starts with a valid modulus, so idle lanes compute harmlessly
    for (lane = 0; lane < U4_VL; lane++) {
        nU[lane] = n64[0];
        nA[lane] = (double)n64[0];
        rhoA[lane] = (double)((0 - uedw4_inverse(n64[0])) & 0x000fffffffffffffULL);
        oneA[lane] = (double)((1ULL << 52) % n64[0]);
        kcur[lane] = 0;
        idx[lane] = 0;
    }
    memset(&R, 0, sizeof(R));

    while (1) {
        unsigned lmsk = 0;
        for (lane = 0; lane < U4_VL; lane++) {
            if (!(msk & (1u << lane)) && j < num) {
                uint64_t n = n64[j];
                nU[lane] = n;
                nA[lane] = (double)n;
                rhoA[lane] = (double)((0 - uedw4_inverse(n)) & 0x000fffffffffffffULL);
                oneA[lane] = (double)((1ULL << 52) % n);
                kcur[lane] = 0;
                idx[lane] = (int)j;
                lmsk |= 1u << lane;
                j++;
            }
        }
        msk |= lmsk;
        if (msk == 0)
            break;
        if (lmsk) {
            ctx.n = u4_load(nA);
            ctx.rho = u4_load(rhoA);
            ctx.one = u4_load(oneA);
            uedw4_ctx_finish(&ctx);
        }

        // this round's curve in every lane: table curves, then GKL m = 3, 4, ...
        unsigned tabm = 0, newm = 0;
        for (lane = 0; lane < U4_VL; lane++) {
            kk[lane] = (kcur[lane] < ntab) ? kcur[lane] : 0;
            if (kcur[lane] < ntab) tabm |= 1u << lane;
            if (kcur[lane] == ntab) newm |= 1u << lane;
        }
        uedw4_curve cv, cvg;
        if (tabm)
            uedw4_table_curve(kk, &cv, c);
        if (tabm != U4_ALL) {
            u4v gm = u4_lanes(~tabm & U4_ALL);
            if (newm) {
                // [3]Q = (1205/49, -40664/343) = (413315 : -1992536 : 16807)
                int64_t q3[3][U4_VL];
                uedw4_wpt Q3;
                for (lane = 0; lane < U4_VL; lane++) {
                    q3[0][lane] = 413315;
                    q3[1][lane] = -1992536;
                    q3[2][lane] = 16807;
                }
                Q3.X = uedw4_const(q3[0], c);
                Q3.Y = uedw4_const(q3[1], c);
                Q3.Z = uedw4_const(q3[2], c);
                uedw4_blend_w(u4_lanes(newm), &R, &Q3);
            }
            uedw4_gkl_curve(&R, &cvg, c);
            uedw4_wpt R2 = R;
            uedw4_w_addq(&R2, c);
            uedw4_blend_w(gm, &R, &R2);
            if (!tabm)
                cv = cvg;
            else
                uedw4_blend_curve(gm, &cv, &cvg);
        }
        for (lane = 0; lane < U4_VL; lane++)
            kcur[lane]++;

        uedw4_pt P;
        uedw4_stage1(&cv, &ch, &P, c);
        u4_store(st, P.Z);
        for (lane = 0; lane < U4_VL; lane++) {
            if (msk & (1u << lane)) {
                total++;
                uint64_t g = uedw4_gcd(nU[lane], (uint64_t)st[lane]);
                if (g > 1 && g != nU[lane]) {
                    f[idx[lane]] = g;
                    msk &= ~(1u << lane);
                }
            }
        }
        if (B2 > B1 && msk) {
            u4v acc = lv ? uedw4_stage2_gen(&P, &cv, lv->s2, lv->s2len, c)
                : uedw4_stage2(&P, &cv, B1, barray, numb, c);
            u4_store(st, acc);
            for (lane = 0; lane < U4_VL; lane++) {
                if (msk & (1u << lane)) {
                    uint64_t g = uedw4_gcd(nU[lane], (uint64_t)st[lane]);
                    if (g > 1 && g != nU[lane]) {
                        f[idx[lane]] = g;
                        msk &= ~(1u << lane);
                    }
                }
            }
        }
        for (lane = 0; lane < U4_VL; lane++)
            if (kcur[lane] >= (int)curves)
                msk &= ~(1u << lane);
    }
    if (ncurves)
        *ncurves = total;
}

/* B1 / curves as microecm's uecm_dispatch_x8_list */
static void uedw4_params(int bits, uint32_t *B1, uint32_t *curves)
{
    if (uecm_edw_get_params() == 2 && bits <= 36) { *B1 = 20; *curves = 32; }
    else if (bits <= 40) { *B1 = 27; *curves = 32; }
    else if (bits <= 44) { *B1 = 47; *curves = 32; }
    else if (bits <= 48) { *B1 = 70; *curves = 32; }
    else if (bits <= 52) { *B1 = 85; *curves = 32; }
    else if (bits <= 58) { *B1 = 125; *curves = 32; }
    else if (bits <= 62) { *B1 = 165; *curves = 42; }
    else { *B1 = 205; *curves = 42; }
}
#endif /* USE_AVX2 */

void getfactor_uecm_edw_x4_list(uint64_t *q64, uint64_t *f64, uint32_t num_in, uint64_t *pran)
{
    uint32_t i;
    if (num_in == 0)
        return;
#ifdef USE_AVX2
    {
        uint64_t n = q64[0];
        int bits = 0;
        while (n) { bits++; n >>= 1; }          // all inputs about the same size
        if (bits <= 52) {
            uint32_t B1, curves;
            uedw4_params(bits, &B1, &curves);
            uedw4_list(q64, f64, num_in, B1, 25 * B1, curves, 0);
            return;
        }
    }
#endif
    for (i = 0; i < num_in; i++)
        f64[i] = getfactor_uecm_edw(q64[i], 0, pran);
}

uint64_t uecm_edw_x4_run(uint64_t *q64, uint64_t *f64, uint32_t num, uint32_t B1, uint32_t curves)
{
#ifdef USE_AVX2
    uint64_t total = 0;
    uedw4_list(q64, f64, num, B1, 25 * B1, curves, &total);
    return total;
#else
    (void)q64; (void)f64; (void)num; (void)B1; (void)curves;
    return 0;
#endif
}
