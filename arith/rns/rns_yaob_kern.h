/*
 * Yao-regime batch kernel body (8 numbers mod the same N, one per lane of
 * the 8-lane batch layout), written against the vector layer of
 * rns_vec52.h.  Included by rns_yaob_avx512.c (RNS_VW 512: one 8-lane
 * vector per digit) and rns_yaob_avx2.c (RNS_VW 256: the same layout as
 * two 4-lane halves, processed one after the other).
 *
 * Channels are the loop index (n per base, any n the constants allow), so
 * every channel modulus, rotated constant and extension constant is a
 * broadcast, the redundant-channel sums are ordinary lane-wise
 * accumulations, and alpha' is a vector (one value per number), applied
 * through a 16-entry table lookup (vpermi2q / a gather) instead of a
 * multiply.  Each channel's s digits are fully unrolled (s = 1..10,
 * compile-time); the channel loops are not.  With RNSY_COLUMN_FOLD = 0 the
 * kernel reproduces rns_yao_ref.c exactly.
 *
 * The including file defines RNSYB_NAME(f) (public function names) and
 * RNSYB_HALVES (1 or 2 halves per call) before including this.
 */

#include "rns_yao_vec.h"
#include <stdlib.h>
#include <string.h>

// acc[0..S+1] += sum_i a_i * tab_i, tab_i = f * 2^(52i) mod m, scalar digits
// tab[i*S + k] broadcast to all lanes
#if RNS_VEC_FP
INL void bcmul_acc(rns_vi *acc, const rns_vi *a, const uint64_t *tab, const int S)
{
    const double *dtab = (const double *)tab;   // FP builds store doubles
    rns_vd ad[MAXS];
    int i, k;
    UNROLL
    for (i = 0; i < S; i++)
        ad[i] = rns_cvt52(a[i]);
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (k = 0; k < S; k++)
            RNS_FP_MAC(acc[k], acc[k + 1], ad[i], V_BCD(dtab[i * S + k]));
    }
    // S lo halves into acc[0..S-1], S hi halves into acc[1..S]
    UNROLL
    for (k = 0; k <= S; k++)
        acc[k] = rns_fp_unbias(acc[k], (k < S) ? S : 0, (k > 0) ? S : 0);
}

// store lo52 of n digits as doubles (FP extension operands)
INL void vstore_d(uint64_t *p, const rns_vi *v, int n)
{
    int t;
    UNROLL
    for (t = 0; t < n; t++)
        V_STD((double *)(p + t * NCH), rns_cvt52(v[t]));
}

// acc[0..S] += sum_i x_i * c_i over n channels of S digits (an extension
// row); FP operands are the double copies, c_i broadcast as doubles
INL void bext_acc(rns_vi *acc, const uint64_t *xd, size_t chw,
    const uint64_t *c, size_t cstride, int n, const int S)
{
    int i, k;
    for (i = 0; i < n; i++)
    {
        rns_vd kc = rns_bcd52(c[(size_t)i * cstride]);
        UNROLL
        for (k = 0; k < S; k++)
        {
            rns_vd v = V_LDD((const double *)(xd + (size_t)i * chw + k * NCH));
            RNS_FP_MAC(acc[k], acc[k + 1], v, kc);
        }
    }
    UNROLL
    for (k = 0; k <= S; k++)
        acc[k] = rns_fp_unbias(acc[k], (k < S) ? n : 0, (k > 0) ? n : 0);
}
#else
INL void bcmul_acc(rns_vi *acc, const rns_vi *a, const uint64_t *tab, const int S)
{
    int i, k;
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (k = 0; k < S; k++)
        {
            rns_vi f = BC(tab[i * S + k]);
            acc[k] = LO52(acc[k], a[i], f);
            acc[k + 1] = HI52(acc[k + 1], a[i], f);
        }
    }
}

INL void bext_acc(rns_vi *acc, const uint64_t *x, size_t chw,
    const uint64_t *c, size_t cstride, int n, const int S)
{
    int i, k;
    for (i = 0; i < n; i++)
    {
        rns_vi kc = BC(c[(size_t)i * cstride]);
        UNROLL
        for (k = 0; k < S; k++)
        {
            rns_vi v = LD(x + (size_t)i * chw + k * NCH);
            acc[k] = LO52(acc[k], v, kc);
            acc[k + 1] = HI52(acc[k + 1], v, kc);
        }
    }
}
#define vstore_d vstore
#endif

// nd normalized digits shifted left by z bits (0 < z < 52); the caller
// guarantees nothing is shifted out of the top digit
INL void vshl(rns_vi *a, const int nd, int z)
{
    rns_vi vz = BC(z), vr = BC(52 - z);
    int t;
    UNROLL
    for (t = nd - 1; t > 0; t--)
        a[t] = V_OR(V_AND(V_SLLV(a[t], vz), M52),
            V_SRLV(a[t - 1], vr));
    a[0] = V_AND(V_SLLV(a[0], vz), M52);
}

// columns acc[0..S] of an extension sum -> digits [0..S+1], restore the
// trailing-zero shift, then fold to x < 2^w
INL void bext_finish(rns_vi *x, rns_vi *acc, int z, rns_vi d, const int S)
{
    int t;
    acc[S + 1] = vnorm(acc, S + 1);
    if (z)
        vshl(acc, S + 2, z);
    UNROLL
    for (t = 0; t < S; t++)
        x[t] = acc[t];
    vfold_top_x(x, acc[S], acc[S + 1], d, S, 1);
}

// x, y, z point at this call's first lane; lo = that lane (scratch offset)
INL void bmulmod_s(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y,
    const int S, const int SQR, const int lo)
{
    rnsy_fixed_t *fx = ctx->fx;
    const int n = fx->n;
    const size_t chw = (size_t)S * NCH;                 // words per channel
    uint64_t *xi = ctx->bscratch + lo;
    uint64_t *xc = x + (size_t)n * chw, *yc = y + (size_t)n * chw, *zc = z + (size_t)n * chw;
    // xi is in bscratch (FP: as doubles)
#if RNS_VEC_FP
    uint64_t *zcd = ctx->bscratch + (size_t)n * chw + lo;    // zc as doubles
#else
    uint64_t *zcd = zc;
#endif
    rns_vi a[MAXS], b[MAXS], t[MAXS], q[MAXS], acc[2 * MAXS + 2];
    rns_vi ur, sr, alpha;
    int i, j, k;

    // redundant channel, started before any z store (z may alias x or y)
    ur = LO52(ZERO, LO52(ZERO, LD(x + fx->rword * NCH), LD(y + fx->rword * NCH)),
        BC(fx->MBinv_r));

    // 1: B channels, xi_i = fold(fold(x_i y_i) f1_i)
    for (i = 0; i < n; i++)
    {
        rns_vi d = BC(fx->db[i]);
        vload(a, x + (size_t)i * chw, S);
        if (SQR)
            vsqr(acc, a, S);
        else
        {
            vload(b, y + (size_t)i * chw, S);
            vzero(acc, 2 * S + 1);
            vmul_acc(acc, a, b, S);
        }
        vfold_prod(t, acc, d, S);
        vzero(acc, S + 2);
        bcmul_acc(acc, t, ctx->bF1 + (size_t)(i * S * S), S);
        vfold_cprod(t, acc, d, S);
        vstore_d(xi + (size_t)i * chw, t, S);
        ur = LO52(ur, t[0], BC(ctx->BrN[i]));
    }
    ur = V_AND(ur, M52);

    // 2-3: per C channel, Q_j = fold(sum_i xi_i Bt_ij << zB_j), then
    // xi'_j = fold(fold(x'_j y'_j) g1_j + Q_j g2_j)
    sr = ZERO;
    for (j = 0; j < n; j++)
    {
        rns_vi d = BC(fx->dc[j]);

        vzero(acc, S + 1);
        bext_acc(acc, xi, chw, &fx->Bt[0][j], RNSY_MAXN, n, S);
        bext_finish(q, acc, fx->zB[j], d, S);

        vload(a, xc + (size_t)j * chw, S);
        if (SQR)
            vsqr(acc, a, S);
        else
        {
            vload(b, yc + (size_t)j * chw, S);
            vzero(acc, 2 * S + 1);
            vmul_acc(acc, a, b, S);
        }
        vfold_prod(t, acc, d, S);
        vzero(acc, S + 2);
        bcmul_acc(acc, t, ctx->bG1 + (size_t)(j * S * S), S);
        bcmul_acc(acc, q, ctx->bG2 + (size_t)(j * S * S), S);
        vfold_cprod(t, acc, d, S);
        vstore(zc + (size_t)j * chw, t, S);
#if RNS_VEC_FP
        vstore_d(zcd + (size_t)j * chw, t, S);
#endif
        sr = LO52(sr, t[0], BC(fx->Cr[j]));
    }

    // alpha' = (sum_j lo52(xi'_j) |C_j|_r - U_r) M_C^{-1} (mod 2^52), per lane
    alpha = LO52(ZERO, V_AND(V_SUB(sr, ur), M52), BC(fx->MCinv_r));

    // 5: U_i = fold(+-(sum_j xi'_j Ct_ji << zC_i) + alpha' |-M_C|_{b_i}), the
    // minus sign (n even) taken as K_i - sum
    for (i = 0; i < n; i++)
    {
        const uint64_t *K = ctx->bK + (size_t)(i * (S + 2));
        const uint64_t *tab = ctx->balpha + (size_t)(i * (S + 1) * 16);
        rns_vi borrow = ZERO;

        vzero(acc, S + 1);
        bext_acc(acc, zcd, chw, &fx->Ct[0][i], RNSY_MAXN, n, S);
        acc[S + 1] = vnorm(acc, S + 1);
        if (fx->zC[i])
            vshl(acc, S + 2, fx->zC[i]);
        if (fx->cneg)
        {
            // n even: |C_j|_{b_i} == -Ct, so take K_i - sum
            UNROLL
            for (k = 0; k < S + 2; k++)
            {
                rns_vi v = V_ADD(V_SUB(BC(K[k]), acc[k]), borrow);
                acc[k] = V_AND(v, M52);
                borrow = V_SRAI52(v);
            }
        }
        UNROLL
        for (k = 0; k <= S; k++)
            acc[k] = V_ADD(acc[k], V_LOOKUP16(tab + k * 16, alpha));
        vnorm(acc, S + 2);
        UNROLL
        for (k = 0; k < S; k++)
            t[k] = acc[k];
        vfold_top_x(t, acc[S], acc[S + 1], BC(fx->db[i]), S, 1);
        vstore(z + (size_t)i * chw, t, S);
    }

    ST(z + fx->rword * NCH, ur);
}

INL void baddsub_s(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y,
    const int S, const int SUB)
{
    rnsy_fixed_t *fx = ctx->fx;
    const int n = fx->n;
    const size_t chw = (size_t)S * NCH;
    int base, ch, k;

    for (base = 0; base < 2; base++)
    {
        for (ch = 0; ch < n; ch++)
        {
            size_t off = (size_t)(base * n + ch) * chw;
            rns_vi v[MAXS + 1], c = ZERO;
            rns_vi d = BC(base ? fx->dc[ch] : fx->db[ch]);
            const uint64_t *E = (base ? ctx->bEc : ctx->bEb) + (size_t)(ch * (S + 1));

            if (SUB)
            {
                // x + (subconst + 2m) - y, signed digit chain, top digit <= 3
                UNROLL
                for (k = 0; k <= S; k++)
                {
                    rns_vi e = V_ADD(BC(E[k]), c);
                    if (k < S)
                        e = V_SUB(V_ADD(e, LD(x + off + k * NCH)),
                            LD(y + off + k * NCH));
                    v[k] = V_AND(e, M52);
                    c = V_SRAI52(e);
                }
            }
            else
            {
                UNROLL
                for (k = 0; k < S; k++)
                    v[k] = V_ADD(LD(x + off + k * NCH), LD(y + off + k * NCH));
                v[S] = vnorm(v, S);
            }
            vfold_top(v, v[S], d, S);
            vstore(z + off, v, S);
        }
    }
    if (SUB)
        ST(z + fx->rword * NCH, V_AND(V_SUB(V_ADD(
            LD(x + fx->rword * NCH), BC(ctx->subconst[fx->rword])),
            LD(y + fx->rword * NCH)), M52));
    else
        ST(z + fx->rword * NCH, V_AND(V_ADD(
            LD(x + fx->rword * NCH), LD(y + fx->rword * NCH)), M52));
}

// per-s instances and dispatch: each instance covers the whole 8-lane
// batch vector (RNSYB_HALVES calls of the vector-width body)
#define HALF_STRIDE (NCH / RNSYB_HALVES)
#define INSTANCES(S) \
    static void bmul_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { int h; for (h = 0; h < RNSYB_HALVES; h++) \
            bmulmod_s(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, S, 0, h * HALF_STRIDE); } \
    static void bsqr_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { int h; for (h = 0; h < RNSYB_HALVES; h++) \
            bmulmod_s(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, S, 1, h * HALF_STRIDE); } \
    static void badd_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { int h; for (h = 0; h < RNSYB_HALVES; h++) \
            baddsub_s(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, S, 0); } \
    static void bsub_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { int h; for (h = 0; h < RNSYB_HALVES; h++) \
            baddsub_s(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, S, 1); }

INSTANCES(1) INSTANCES(2) INSTANCES(3) INSTANCES(4) INSTANCES(5)
INSTANCES(6) INSTANCES(7) INSTANCES(8) INSTANCES(9) INSTANCES(10)

typedef void (*opfn)(rnsy_ctx_t *, uint64_t *, uint64_t *, uint64_t *);
#define TABLE(op) static const opfn op##_tab[MAXS + 1] = { NULL, \
    op##_1, op##_2, op##_3, op##_4, op##_5, op##_6, op##_7, op##_8, op##_9, op##_10 }
TABLE(bmul); TABLE(bsqr); TABLE(badd); TABLE(bsub);

// ---------------------------------------------------------------------------
// Large s (> MAXS = RNSY_MAXS_UNROLL): loop forms with column-wise register
// accumulation (each output column summed in registers, stored once),
// scratch from ctx->bbig instead of stack arrays, and always the exact
// product fold (the column-fold bounds need column sums < 2^58, which fails
// from s ~ 30).  So the large-s kernel matches rns_yao_ref.c bit for bit in
// both RNSY_COLUMN_FOLD modes.  Channel products use karatsuba down to
// rnsyb_kara_min words (0: schoolbook only).
// ---------------------------------------------------------------------------

#if RNS_VEC_FP
typedef rns_vd big_op;
#define BOP(x)          rns_cvt52(x)
#define BMAC(l, h, x, y) RNS_FP_MAC(l, h, x, y)
#define BFIX(l, h, cnt) do { (l) = rns_fp_unbias((l), (cnt), 0); (h) = rns_fp_unbias((h), 0, (cnt)); } while (0)
#define BCK(c)          rns_bcd52(c)
#define BTAB(tab, idx)  V_BCD(((const double *)(tab))[idx])    // FP tables are doubles
#define BLDX(p)         V_LDD(p)                                // xi / zc double copies
#else
typedef rns_vi big_op;
#define BOP(x)          (x)
#define BMAC(l, h, x, y) do { (l) = LO52((l), (x), (y)); (h) = HI52((h), (x), (y)); } while (0)
#define BFIX(l, h, cnt) do { (void)(cnt); } while (0)
#define BCK(c)          BC(c)
#define BTAB(tab, idx)  BC((tab)[idx])
#define BLDX(p)         LD(p)
#endif

// acc[0..na+nb-1] += a * b as lo/hi column sums (operands converted)
static void big_mul_cols(rns_vi *acc, const big_op *a, int na, const big_op *b, int nb)
{
    int c;
    for (c = 0; c < na + nb - 1; c++)
    {
        int ilo = (c - (nb - 1) > 0) ? c - (nb - 1) : 0;
        int ihi = (c < na - 1) ? c : na - 1;
        rns_vi l0 = ZERO, h0 = ZERO, l1 = ZERO, h1 = ZERO;
        int i;
        for (i = ilo; i + 1 <= ihi; i += 2)
        {
            BMAC(l0, h0, a[i], b[c - i]);
            BMAC(l1, h1, a[i + 1], b[c - i - 1]);
        }
        if (i <= ihi)
            BMAC(l0, h0, a[i], b[c - i]);
        l0 = V_ADD(l0, l1);
        h0 = V_ADD(h0, h1);
        BFIX(l0, h0, ihi - ilo + 1);
        acc[c] = V_ADD(acc[c], l0);
        acc[c + 1] = V_ADD(acc[c + 1], h0);
    }
}

// acc[0..2n-1] += a^2: off-diagonal pairs doubled, plus the diagonal
static void big_sqr_cols(rns_vi *acc, const big_op *a, int n)
{
    int c;
    for (c = 0; c < 2 * n - 1; c++)
    {
        int ilo = (c - (n - 1) > 0) ? c - (n - 1) : 0;
        rns_vi l0 = ZERO, h0 = ZERO, l1 = ZERO, h1 = ZERO;
        int i, cnt = 0;
        for (i = ilo; 2 * (i + 1) < c; i += 2)
        {
            BMAC(l0, h0, a[i], a[c - i]);
            BMAC(l1, h1, a[i + 1], a[c - i - 1]);
            cnt += 2;
        }
        if (2 * i < c)
        {
            BMAC(l0, h0, a[i], a[c - i]);
            cnt++;
        }
        l0 = V_ADD(l0, l1);
        h0 = V_ADD(h0, h1);
        BFIX(l0, h0, cnt);
        l0 = V_ADD(l0, l0);
        h0 = V_ADD(h0, h0);
        if ((c & 1) == 0)
        {
            rns_vi l2 = ZERO, h2 = ZERO;
            BMAC(l2, h2, a[c / 2], a[c / 2]);
            BFIX(l2, h2, 1);
            l0 = V_ADD(l0, l2);
            h0 = V_ADD(h0, h2);
        }
        acc[c] = V_ADD(acc[c], l0);
        acc[c + 1] = V_ADD(acc[c + 1], h0);
    }
}

// r[0..2n] = a * b (a == b: square), normalized digits, r[2n] = 0;
// a, b normalized n-digit numbers; scr: scratch vectors (~10 n)
static void kprod(rns_vi *r, const rns_vi *a, const rns_vi *b, int n, int sqr, rns_vi *scr)
{
    int i;
    if (rnsyb_kara_min <= 0 || n < 2 * rnsyb_kara_min || n < 4)
    {
        big_op *ao = (big_op *)(void *)scr, *bo = ao + n;
        for (i = 0; i <= 2 * n; i++)
            r[i] = ZERO;
        for (i = 0; i < n; i++)
        {
            ao[i] = BOP(a[i]);
            if (!sqr)
                bo[i] = BOP(b[i]);
        }
        if (sqr)
            big_sqr_cols(r, ao, n);
        else
            big_mul_cols(r, ao, n, bo, n);
        r[2 * n] = vnorm(r, 2 * n);
        return;
    }
    {
        // a = a0 + a1 B^l: r = P0 + (P1 - P0 - P2) B^l + P2 B^2l,
        // P1 = (a0 + a1)(b0 + b1).  Subproducts are normalized, so the
        // middle term is an ordinary (nonnegative) multiprecision difference.
        const int l = (n + 1) / 2, h = n - l;
        rns_vi *sa = scr, *sb = sa + (l + 1), *p1 = sb + (l + 1);
        rns_vi *rest = p1 + (2 * l + 3);
        rns_vi c, bo;

        kprod(r, a, b, l, sqr, rest);                       // P0: r[0..2l]
        kprod(r + 2 * l, a + l, b + l, h, sqr, rest);       // P2: r[2l..2n]
        for (i = 0, c = ZERO; i < l; i++)
        {
            rns_vi v = V_ADD(a[i], c);
            if (i < h) v = V_ADD(v, a[l + i]);
            sa[i] = V_AND(v, M52);
            c = V_SRLI(v, 52);
        }
        sa[l] = c;
        if (!sqr)
        {
            for (i = 0, c = ZERO; i < l; i++)
            {
                rns_vi v = V_ADD(b[i], c);
                if (i < h) v = V_ADD(v, b[l + i]);
                sb[i] = V_AND(v, M52);
                c = V_SRLI(v, 52);
            }
            sb[l] = c;
        }
        kprod(p1, sa, sqr ? sa : sb, l + 1, sqr, rest);    // P1: p1[0..2l+2]
        // p1 -= P0 (r[0..2l-1]) and P2 (r[2l..2n-1]), signed borrow chain
        for (i = 0, bo = ZERO; i < 2 * l + 2; i++)
        {
            rns_vi v = V_ADD(p1[i], bo);
            if (i < 2 * l) v = V_SUB(v, r[i]);
            if (i < 2 * h) v = V_SUB(v, r[2 * l + i]);
            p1[i] = V_AND(v, M52);
            bo = V_SRAI52(v);
        }
        // r += mid * B^l (the sum is < B^(2n), so the carry dies out)
        for (i = 0, c = ZERO; l + i < 2 * n; i++)
        {
            rns_vi v = V_ADD(r[l + i], c);
            if (i < 2 * l + 2) v = V_ADD(v, p1[i]);
            r[l + i] = V_AND(v, M52);
            c = V_SRLI(v, 52);
        }
        r[2 * n] = ZERO;
    }
}

// columns acc[0..2S] (normalized digits) -> x < 2^w, the exact fold
static void big_fold_prod(rns_vi *x, rns_vi *acc, rns_vi d, int S, rns_vi *y)
{
    int t;
    acc[2 * S] = vnorm(acc, 2 * S);
    for (t = 0; t < S; t++)
        y[t] = acc[t];
    y[S] = ZERO;
    for (t = 0; t <= S; t++)
    {
        y[t] = LO52(y[t], acc[S + t], d);
        if (t < S)
            y[t + 1] = HI52(y[t + 1], acc[S + t], d);
    }
    vnorm(y, S + 1);
    for (t = 0; t < S; t++)
        x[t] = y[t];
    vfold_top(x, y[S], d, S);
}

// acc[0..S] += sum_i a_i * tab[i*S + k] (rotated constants), column-wise
static void big_cmul(rns_vi *acc, const rns_vi *a, const uint64_t *tab, int S, big_op *ao)
{
    int i, k;
    for (i = 0; i < S; i++)
        ao[i] = BOP(a[i]);
    for (k = 0; k < S; k++)
    {
        rns_vi l0 = ZERO, h0 = ZERO;
        for (i = 0; i < S; i++)
            BMAC(l0, h0, ao[i], BTAB(tab, (size_t)i * (size_t)S + (size_t)k));
        BFIX(l0, h0, S);
        acc[k] = V_ADD(acc[k], l0);
        acc[k + 1] = V_ADD(acc[k + 1], h0);
    }
}

// acc[0..S] += sum_i x_i * c_i over n channels (an extension row), column-wise
static void big_ext(rns_vi *acc, const uint64_t *xd, size_t chw,
    const uint64_t *c, size_t cstride, int n, int S)
{
    big_op kc[RNSY_MAXN];
    int i, k;
    for (i = 0; i < n; i++)
        kc[i] = BCK(c[(size_t)i * cstride]);
    for (k = 0; k < S; k++)
    {
        rns_vi l0 = ZERO, h0 = ZERO;
        for (i = 0; i < n; i++)
            BMAC(l0, h0, BLDX(xd + (size_t)i * chw + (size_t)k * NCH), kc[i]);
        BFIX(l0, h0, n);
        acc[k] = V_ADD(acc[k], l0);
        acc[k + 1] = V_ADD(acc[k + 1], h0);
    }
}

static void bmulmod_big(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y,
    const int SQR, const int lo)
{
    rnsy_fixed_t *fx = ctx->fx;
    const int n = fx->n, S = fx->s;
    const size_t chw = (size_t)S * NCH;
    uint64_t *xi = ctx->bscratch + lo;
    uint64_t *xc = x + (size_t)n * chw, *yc = y + (size_t)n * chw, *zc = z + (size_t)n * chw;
#if RNS_VEC_FP
    uint64_t *zcd = ctx->bscratch + (size_t)n * chw + lo;
#else
    uint64_t *zcd = zc;
#endif
    rns_vi *a = (rns_vi *)(void *)ctx->bbig;
    rns_vi *b = a + S, *t = b + S, *q = t + S, *acc = q + S, *yv = acc + (2 * S + 4);
    big_op *ao = (big_op *)(void *)(yv + (S + 2));
    rns_vi *r2 = (rns_vi *)(void *)(ao + S);
    rns_vi *kscr = r2 + (2 * S + 4);
    // karatsuba for all five products (also the constant ones: plain digits
    // of f1 / g1 / g2, a full product, then the exact fold) when it applies
    const int kc = (rnsyb_kara_min > 0) && (S >= 2 * rnsyb_kara_min) && (ctx->bkconst != NULL);
    const rns_vi *KF1 = (const rns_vi *)(const void *)ctx->bkconst;
    const rns_vi *KG1 = KF1 + (size_t)n * (size_t)S, *KG2 = KG1 + (size_t)n * (size_t)S;
    rns_vi ur, sr, alpha;
    int i, j, k;

    ur = LO52(ZERO, LO52(ZERO, LD(x + fx->rword * NCH), LD(y + fx->rword * NCH)),
        BC(fx->MBinv_r));

    for (i = 0; i < n; i++)
    {
        rns_vi d = BC(fx->db[i]);
        vload(a, x + (size_t)i * chw, S);
        if (!SQR)
            vload(b, y + (size_t)i * chw, S);
        kprod(acc, a, SQR ? a : b, S, SQR, kscr);
        big_fold_prod(t, acc, d, S, yv);
        if (kc)
        {
            kprod(acc, t, KF1 + (size_t)i * (size_t)S, S, 0, kscr);
            big_fold_prod(t, acc, d, S, yv);
        }
        else
        {
            vzero(acc, S + 2);
            big_cmul(acc, t, ctx->bF1 + (size_t)i * (size_t)S * (size_t)S, S, ao);
            vfold_cprod(t, acc, d, S);
        }
        vstore_d(xi + (size_t)i * chw, t, S);
        ur = LO52(ur, t[0], BC(ctx->BrN[i]));
    }
    ur = V_AND(ur, M52);

    sr = ZERO;
    for (j = 0; j < n; j++)
    {
        rns_vi d = BC(fx->dc[j]);
        vzero(acc, S + 1);
        big_ext(acc, xi, chw, &fx->Bt[0][j], RNSY_MAXN, n, S);
        bext_finish(q, acc, fx->zB[j], d, S);

        vload(a, xc + (size_t)j * chw, S);
        if (!SQR)
            vload(b, yc + (size_t)j * chw, S);
        kprod(acc, a, SQR ? a : b, S, SQR, kscr);
        big_fold_prod(t, acc, d, S, yv);
        if (kc)
        {
            // t*g1 + q*g2: two full products summed digit-wise (< 2^53
            // per digit), then one fold
            kprod(acc, t, KG1 + (size_t)j * (size_t)S, S, 0, kscr);
            kprod(r2, q, KG2 + (size_t)j * (size_t)S, S, 0, kscr);
            for (k = 0; k < 2 * S; k++)
                acc[k] = V_ADD(acc[k], r2[k]);
            big_fold_prod(t, acc, d, S, yv);
        }
        else
        {
            vzero(acc, S + 2);
            big_cmul(acc, t, ctx->bG1 + (size_t)j * (size_t)S * (size_t)S, S, ao);
            big_cmul(acc, q, ctx->bG2 + (size_t)j * (size_t)S * (size_t)S, S, ao);
            vfold_cprod(t, acc, d, S);
        }
        vstore(zc + (size_t)j * chw, t, S);
#if RNS_VEC_FP
        vstore_d(zcd + (size_t)j * chw, t, S);
#endif
        sr = LO52(sr, t[0], BC(fx->Cr[j]));
    }

    alpha = LO52(ZERO, V_AND(V_SUB(sr, ur), M52), BC(fx->MCinv_r));

    for (i = 0; i < n; i++)
    {
        const uint64_t *K = ctx->bK + (size_t)i * (size_t)(S + 2);
        const uint64_t *tab = ctx->balpha + (size_t)i * (size_t)(S + 1) * 16;
        rns_vi borrow = ZERO;

        vzero(acc, S + 1);
        big_ext(acc, zcd, chw, &fx->Ct[0][i], RNSY_MAXN, n, S);
        acc[S + 1] = vnorm(acc, S + 1);
        if (fx->zC[i])
            vshl(acc, S + 2, fx->zC[i]);
        if (fx->cneg)
        {
            for (k = 0; k < S + 2; k++)
            {
                rns_vi v = V_ADD(V_SUB(BC(K[k]), acc[k]), borrow);
                acc[k] = V_AND(v, M52);
                borrow = V_SRAI52(v);
            }
        }
        for (k = 0; k <= S; k++)
            acc[k] = V_ADD(acc[k], V_LOOKUP16(tab + k * 16, alpha));
        vnorm(acc, S + 2);
        for (k = 0; k < S; k++)
            t[k] = acc[k];
        vfold_top_x(t, acc[S], acc[S + 1], BC(fx->db[i]), S, 1);
        vstore(z + (size_t)i * chw, t, S);
    }

    ST(z + fx->rword * NCH, ur);
}

static void baddsub_big(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y, const int SUB)
{
    rnsy_fixed_t *fx = ctx->fx;
    const int n = fx->n, S = fx->s;
    const size_t chw = (size_t)S * NCH;
    rns_vi *v = (rns_vi *)(void *)ctx->bbig;
    int base, ch, k;

    for (base = 0; base < 2; base++)
    {
        for (ch = 0; ch < n; ch++)
        {
            size_t off = (size_t)(base * n + ch) * chw;
            rns_vi c = ZERO;
            rns_vi d = BC(base ? fx->dc[ch] : fx->db[ch]);
            const uint64_t *E = (base ? ctx->bEc : ctx->bEb) + (size_t)ch * (size_t)(S + 1);

            if (SUB)
            {
                for (k = 0; k <= S; k++)
                {
                    rns_vi e = V_ADD(BC(E[k]), c);
                    if (k < S)
                        e = V_SUB(V_ADD(e, LD(x + off + (size_t)k * NCH)), LD(y + off + (size_t)k * NCH));
                    v[k] = V_AND(e, M52);
                    c = V_SRAI52(e);
                }
            }
            else
            {
                for (k = 0; k < S; k++)
                    v[k] = V_ADD(LD(x + off + (size_t)k * NCH), LD(y + off + (size_t)k * NCH));
                v[S] = vnorm(v, S);
            }
            vfold_top(v, v[S], d, S);
            vstore(z + off, v, S);
        }
    }
    if (SUB)
        ST(z + fx->rword * NCH, V_AND(V_SUB(V_ADD(
            LD(x + fx->rword * NCH), BC(ctx->subconst[fx->rword])),
            LD(y + fx->rword * NCH)), M52));
    else
        ST(z + fx->rword * NCH, V_AND(V_ADD(
            LD(x + fx->rword * NCH), LD(y + fx->rword * NCH)), M52));
}

#define BIG_INSTANCE(name, call) \
    static void name(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { int h; for (h = 0; h < RNSYB_HALVES; h++) call; }
BIG_INSTANCE(bmul_big, bmulmod_big(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, 0, h * HALF_STRIDE))
BIG_INSTANCE(bsqr_big, bmulmod_big(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, 1, h * HALF_STRIDE))
BIG_INSTANCE(badd_big, baddsub_big(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, 0))
BIG_INSTANCE(bsub_big, baddsub_big(c, z + h * HALF_STRIDE, x + h * HALF_STRIDE, y + h * HALF_STRIDE, 1))


// the op for this context's s: unrolled instances up to MAXS, loop forms above
static int check_s(rnsy_ctx_t *ctx)
{
    if (ctx->fx->s < 1 || ctx->fx->s > RNSY_MAXS || ctx->fx->n + 2 > 16 ||
        (ctx->fx->s > MAXS && ctx->bbig == NULL))
    {
        fprintf(stderr, "rnsyb: no batch kernel for n = %d, s = %d\n",
            ctx->fx->n, ctx->fx->s);
        exit(1);
    }
    return ctx->fx->s;
}

#if RNS_VW == 256
// (the AVX-512 build defines this in rns_yao_avx512.c)
int rnsy_vec_exact(void)
{
    // (the large-s batch kernel with karatsuba, s >= 2 * rnsyb_kara_min,
    // is not exact: its constant products fold differently)
    return !RNSY_COLUMN_FOLD;
}
#endif

#define PICK(op) static opfn pick_##op(rnsy_ctx_t *ctx) \
    { int s_ = check_s(ctx); return (s_ <= MAXS) ? op##_tab[s_] : op##_big; }
PICK(bmul) PICK(bsqr) PICK(badd) PICK(bsub)

// 256: FP products (also in the add/sub folds) need MXCSR round-toward-zero
// for the duration of a call
#if RNS_VW == 256
#define RZ_ENTER unsigned int csr_ = rns_fp_rz_enter();
#define RZ_LEAVE rns_fp_rz_leave(csr_);
#else
#define RZ_ENTER
#define RZ_LEAVE
#endif

void RNSYB_NAME(rnsyb_vmulmod)(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    RZ_ENTER
    pick_bmul(ctx)(ctx, z, x, y);
    RZ_LEAVE
}

void RNSYB_NAME(rnsyb_vsqrmod)(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x)
{
    RZ_ENTER
    pick_bsqr(ctx)(ctx, z, x, x);
    RZ_LEAVE
}

void RNSYB_NAME(rnsyb_vaddmod)(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    RZ_ENTER
    pick_badd(ctx)(ctx, z, x, y);
    RZ_LEAVE
}

void RNSYB_NAME(rnsyb_vsubmod)(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    RZ_ENTER
    pick_bsub(ctx)(ctx, z, x, y);
    RZ_LEAVE
}

static int winsize(int bits)
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

void RNSYB_NAME(rnsyb_vpowm)(rnsy_ctx_t *ctx, uint64_t *d, uint64_t *b, mpz_t e)
{
    opfn mul = pick_bmul(ctx), sqr = pick_bsqr(ctx);
    int bit = (int)mpz_sizeinbase(e, 2) - 1;
    int k = winsize(bit), ntab = 1 << (k - 1), i, j, lane;
    size_t bytes = (size_t)ctx->fx->nwords * NCH * sizeof(uint64_t);
    uint64_t **g = (uint64_t **)calloc((size_t)ntab, sizeof(uint64_t *));
    uint64_t *b2 = rnsyb_alloc(ctx);

    for (i = 0; i < ntab; i++)
        g[i] = rnsyb_alloc(ctx);
    for (lane = 0; lane < NCH; lane++)
        rnsyb_vset_lane(ctx, d, lane, ctx->one);
    memcpy(g[0], b, bytes);
    {
        RZ_ENTER
        sqr(ctx, b2, b, b);
        for (i = 1; i < ntab; i++)
            mul(ctx, g[i], g[i - 1], b2);
        while (bit >= 0)
        {
            if (mpz_tstbit(e, (mp_bitcnt_t)bit) == 0)
            {
                sqr(ctx, d, d, d);
                bit--;
            }
            else
            {
                int a = bit - k + 1, thisk, bstr = 1;
                if (a < 0)
                    a = 0;
                a = (int)mpz_scan1(e, (mp_bitcnt_t)a);
                thisk = bit - a + 1;
                for (j = 1; j < thisk; j++)
                    bstr = (bstr << 1) | mpz_tstbit(e, (mp_bitcnt_t)(bit - j));
                for (j = 0; j < thisk; j++)
                    sqr(ctx, d, d, d);
                mul(ctx, d, d, g[bstr >> 1]);
                bit -= thisk;
            }
        }
        RZ_LEAVE
    }

    for (i = 0; i < ntab; i++)
        rns_vfree(g[i]);
    free(g);
    rns_vfree(b2);
}
