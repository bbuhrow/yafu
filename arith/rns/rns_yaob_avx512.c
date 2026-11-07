/*
 * RNS Montgomery modular multiplication, Yao regime, batch AVX-512 IFMA
 * kernels: 8 independent numbers modulo the same N per call, one per lane.
 *
 * Channels are the loop index (n per base, any n the constants allow), so
 * every channel modulus, rotated constant and extension constant is a
 * broadcast, the redundant-channel sums are ordinary lane-wise
 * accumulations, and alpha' is a vector (one value per number), applied
 * through a 16-entry table lookup (vpermi2q) instead of a multiply.
 * Each channel's s digits are fully unrolled (s = 1..10, compile-time);
 * the channel loops are not.  The per-channel arithmetic is the same as the
 * one-number kernel's, so with RNSY_COLUMN_FOLD = 0 both kernels reproduce
 * rns_yao_ref.c exactly.
 */

#include "rns_yao.h"
#include <immintrin.h>
#include <stdlib.h>
#include <string.h>

#if defined(__AVX512F__) && defined(__AVX512IFMA__)

#include "rns_yao_vec.h"

// acc[0..S+1] += sum_i a_i * tab_i, tab_i = f * 2^(52i) mod m, scalar digits
// tab[i*S + k] broadcast to all lanes
INL void bcmul_acc(__m512i *acc, const __m512i *a, const uint64_t *tab, const int S)
{
    int i, k;
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (k = 0; k < S; k++)
        {
            __m512i f = BC(tab[i * S + k]);
            acc[k] = LO52(acc[k], a[i], f);
            acc[k + 1] = HI52(acc[k + 1], a[i], f);
        }
    }
}

// nd normalized digits shifted left by z bits (0 < z < 52); the caller
// guarantees nothing is shifted out of the top digit
INL void vshl(__m512i *a, const int nd, int z)
{
    __m512i vz = BC(z), vr = BC(52 - z);
    int t;
    UNROLL
    for (t = nd - 1; t > 0; t--)
        a[t] = _mm512_or_si512(_mm512_and_si512(_mm512_sllv_epi64(a[t], vz), M52),
            _mm512_srlv_epi64(a[t - 1], vr));
    a[0] = _mm512_and_si512(_mm512_sllv_epi64(a[0], vz), M52);
}

// columns acc[0..S] of an extension sum -> digits [0..S+1], restore the
// trailing-zero shift, then fold to x < 2^w
INL void bext_finish(__m512i *x, __m512i *acc, int z, __m512i d, const int S)
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

INL void bmulmod_s(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y,
    const int S, const int SQR)
{
    rnsy_fixed_t *fx = ctx->fx;
    const int n = fx->n;
    const size_t chw = (size_t)S * NCH;                 // words per channel
    uint64_t *xi = ctx->bscratch;
    uint64_t *xc = x + n * chw, *yc = y + n * chw, *zc = z + n * chw;
    __m512i a[MAXS], b[MAXS], t[MAXS], q[MAXS], acc[2 * MAXS + 2];
    __m512i ur, sr, alpha;
    int i, j, k;

    // redundant channel, started before any z store (z may alias x or y)
    ur = LO52(ZERO, LO52(ZERO, LD(x + fx->rword * NCH), LD(y + fx->rword * NCH)),
        BC(fx->MBinv_r));

    // 1: B channels, xi_i = fold(fold(x_i y_i) f1_i)
    for (i = 0; i < n; i++)
    {
        __m512i d = BC(fx->db[i]);
        vload(a, x + i * chw, S);
        if (SQR)
            vsqr(acc, a, S);
        else
        {
            vload(b, y + i * chw, S);
            vzero(acc, 2 * S + 1);
            vmul_acc(acc, a, b, S);
        }
        vfold_prod(t, acc, d, S);
        vzero(acc, S + 2);
        bcmul_acc(acc, t, ctx->bF1 + (size_t)i * S * S, S);
        vfold_cprod(t, acc, d, S);
        vstore(xi + i * chw, t, S);
        ur = LO52(ur, t[0], BC(ctx->BrN[i]));
    }
    ur = _mm512_and_si512(ur, M52);

    // 2-3: per C channel, Q_j = fold(sum_i xi_i Bt_ij << zB_j), then
    // xi'_j = fold(fold(x'_j y'_j) g1_j + Q_j g2_j)
    sr = ZERO;
    for (j = 0; j < n; j++)
    {
        __m512i d = BC(fx->dc[j]);

        vzero(acc, S + 1);
        for (i = 0; i < n; i++)
        {
            __m512i kc = BC(fx->Bt[i][j]);
            UNROLL
            for (k = 0; k < S; k++)
            {
                __m512i v = LD(xi + i * chw + k * NCH);
                acc[k] = LO52(acc[k], v, kc);
                acc[k + 1] = HI52(acc[k + 1], v, kc);
            }
        }
        bext_finish(q, acc, fx->zB[j], d, S);

        vload(a, xc + j * chw, S);
        if (SQR)
            vsqr(acc, a, S);
        else
        {
            vload(b, yc + j * chw, S);
            vzero(acc, 2 * S + 1);
            vmul_acc(acc, a, b, S);
        }
        vfold_prod(t, acc, d, S);
        vzero(acc, S + 2);
        bcmul_acc(acc, t, ctx->bG1 + (size_t)j * S * S, S);
        bcmul_acc(acc, q, ctx->bG2 + (size_t)j * S * S, S);
        vfold_cprod(t, acc, d, S);
        vstore(zc + j * chw, t, S);
        sr = LO52(sr, t[0], BC(fx->Cr[j]));
    }

    // alpha' = (sum_j lo52(xi'_j) |C_j|_r - U_r) M_C^{-1} (mod 2^52), per lane
    alpha = LO52(ZERO, _mm512_and_si512(_mm512_sub_epi64(sr, ur), M52), BC(fx->MCinv_r));

    // 5: U_i = fold(+-(sum_j xi'_j Ct_ji << zC_i) + alpha' |-M_C|_{b_i}), the
    // minus sign (n even) taken as K_i - sum
    for (i = 0; i < n; i++)
    {
        const uint64_t *K = ctx->bK + (size_t)i * (S + 2);
        const uint64_t *tab = ctx->balpha + (size_t)i * (S + 1) * 16;
        __m512i borrow = ZERO;

        vzero(acc, S + 1);
        for (j = 0; j < n; j++)
        {
            __m512i kc = BC(fx->Ct[j][i]);
            UNROLL
            for (k = 0; k < S; k++)
            {
                __m512i v = LD(zc + j * chw + k * NCH);
                acc[k] = LO52(acc[k], v, kc);
                acc[k + 1] = HI52(acc[k + 1], v, kc);
            }
        }
        acc[S + 1] = vnorm(acc, S + 1);
        if (fx->zC[i])
            vshl(acc, S + 2, fx->zC[i]);
        if (fx->cneg)
        {
            // n even: |C_j|_{b_i} == -Ct, so take K_i - sum
            UNROLL
            for (k = 0; k < S + 2; k++)
            {
                __m512i v = _mm512_add_epi64(_mm512_sub_epi64(BC(K[k]), acc[k]), borrow);
                acc[k] = _mm512_and_si512(v, M52);
                borrow = _mm512_srai_epi64(v, 52);
            }
        }
        UNROLL
        for (k = 0; k <= S; k++)
            acc[k] = _mm512_add_epi64(acc[k], _mm512_permutex2var_epi64(
                LD(tab + k * 16), alpha, LD(tab + k * 16 + 8)));
        vnorm(acc, S + 2);
        UNROLL
        for (k = 0; k < S; k++)
            t[k] = acc[k];
        vfold_top_x(t, acc[S], acc[S + 1], BC(fx->db[i]), S, 1);
        vstore(z + i * chw, t, S);
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
            size_t off = (base * n + ch) * chw;
            __m512i v[MAXS + 1], c = ZERO;
            __m512i d = BC(base ? fx->dc[ch] : fx->db[ch]);
            const uint64_t *E = (base ? ctx->bEc : ctx->bEb) + (size_t)ch * (S + 1);

            if (SUB)
            {
                // x + (subconst + 2m) - y, signed digit chain, top digit <= 3
                UNROLL
                for (k = 0; k <= S; k++)
                {
                    __m512i e = _mm512_add_epi64(BC(E[k]), c);
                    if (k < S)
                        e = _mm512_sub_epi64(_mm512_add_epi64(e, LD(x + off + k * NCH)),
                            LD(y + off + k * NCH));
                    v[k] = _mm512_and_si512(e, M52);
                    c = _mm512_srai_epi64(e, 52);
                }
            }
            else
            {
                UNROLL
                for (k = 0; k < S; k++)
                    v[k] = _mm512_add_epi64(LD(x + off + k * NCH), LD(y + off + k * NCH));
                v[S] = vnorm(v, S);
            }
            vfold_top(v, v[S], d, S);
            vstore(z + off, v, S);
        }
    }
    if (SUB)
        ST(z + fx->rword * NCH, _mm512_and_si512(_mm512_sub_epi64(_mm512_add_epi64(
            LD(x + fx->rword * NCH), BC(ctx->subconst[fx->rword])),
            LD(y + fx->rword * NCH)), M52));
    else
        ST(z + fx->rword * NCH, _mm512_and_si512(_mm512_add_epi64(
            LD(x + fx->rword * NCH), LD(y + fx->rword * NCH)), M52));
}

// per-s instances and dispatch
#define INSTANCES(S) \
    static void bmul_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { bmulmod_s(c, z, x, y, S, 0); } \
    static void bsqr_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { bmulmod_s(c, z, x, y, S, 1); } \
    static void badd_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { baddsub_s(c, z, x, y, S, 0); } \
    static void bsub_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { baddsub_s(c, z, x, y, S, 1); }

INSTANCES(1) INSTANCES(2) INSTANCES(3) INSTANCES(4) INSTANCES(5)
INSTANCES(6) INSTANCES(7) INSTANCES(8) INSTANCES(9) INSTANCES(10)

typedef void (*opfn)(rnsy_ctx_t *, uint64_t *, uint64_t *, uint64_t *);
#define TABLE(op) static const opfn op##_tab[MAXS + 1] = { NULL, \
    op##_1, op##_2, op##_3, op##_4, op##_5, op##_6, op##_7, op##_8, op##_9, op##_10 }
TABLE(bmul); TABLE(bsqr); TABLE(badd); TABLE(bsub);

static int check_s(rnsy_ctx_t *ctx)
{
    if (ctx->fx->s < 1 || ctx->fx->s > MAXS || ctx->fx->n + 2 > 16)
    {
        fprintf(stderr, "rnsyb: no batch kernel for n = %d, s = %d\n",
            ctx->fx->n, ctx->fx->s);
        exit(1);
    }
    return ctx->fx->s;
}

void rnsyb_vmulmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    bmul_tab[check_s(ctx)](ctx, z, x, y);
}

void rnsyb_vsqrmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x)
{
    bsqr_tab[check_s(ctx)](ctx, z, x, x);
}

void rnsyb_vaddmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    badd_tab[check_s(ctx)](ctx, z, x, y);
}

void rnsyb_vsubmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    bsub_tab[check_s(ctx)](ctx, z, x, y);
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

void rnsyb_vpowm_avx512(rnsy_ctx_t *ctx, uint64_t *d, uint64_t *b, mpz_t e)
{
    int S = check_s(ctx);
    opfn mul = bmul_tab[S], sqr = bsqr_tab[S];
    int bit = (int)mpz_sizeinbase(e, 2) - 1;
    int k = winsize(bit), ntab = 1 << (k - 1), i, j, lane;
    size_t bytes = (size_t)ctx->fx->nwords * NCH * sizeof(uint64_t);
    uint64_t **g = (uint64_t **)malloc(ntab * sizeof(uint64_t *));
    uint64_t *b2 = rnsyb_alloc(ctx);

    for (i = 0; i < ntab; i++)
        g[i] = rnsyb_alloc(ctx);
    memcpy(g[0], b, bytes);
    sqr(ctx, b2, b, b);
    for (i = 1; i < ntab; i++)
        mul(ctx, g[i], g[i - 1], b2);

    for (lane = 0; lane < NCH; lane++)
        rnsyb_vset_lane(ctx, d, lane, ctx->one);
    while (bit >= 0)
    {
        if (mpz_tstbit(e, bit) == 0)
        {
            sqr(ctx, d, d, d);
            bit--;
        }
        else
        {
            int a = bit - k + 1, thisk, bstr = 1;
            if (a < 0)
                a = 0;
            a = (int)mpz_scan1(e, a);
            thisk = bit - a + 1;
            for (j = 1; j < thisk; j++)
                bstr = (bstr << 1) | mpz_tstbit(e, bit - j);
            for (j = 0; j < thisk; j++)
                sqr(ctx, d, d, d);
            mul(ctx, d, d, g[bstr >> 1]);
            bit -= thisk;
        }
    }

    for (i = 0; i < ntab; i++)
        rns_vfree(g[i]);
    free(g);
    rns_vfree(b2);
}

#endif
