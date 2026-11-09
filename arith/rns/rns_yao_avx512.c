/*
 * RNS Montgomery modular multiplication, Yao regime, AVX-512 IFMA kernels.
 *
 * One number per call.  Each base's 8 channels (2^w - d, w = 52*s) sit in
 * the 8 lanes, held as s digit vectors of 52 bits.  Channel products are
 * fully unrolled column sums (lo/hi IFMA into 64-bit columns, normalized
 * once), specialized per s = 1..10.  The extensions broadcast one channel's
 * digit from memory against a column of single-word Yao constants.
 * With RNSY_COLUMN_FOLD = 0 it mirrors rns_yao_ref.c exactly: every fold is
 * X -> (X mod 2^w) + (X>>w)*d.
 */

#include "rns_yao.h"
#include <stdlib.h>
#include <string.h>

#ifdef RNS_HAVE_AVX512
#include <immintrin.h>

#include "rns_yao_vec.h"

int rnsy_vec_exact(void)
{
    // (the large-s batch kernel with karatsuba, s >= 2 * rnsyb_kara_min,
    // is not exact: its constant products fold differently)
    return !RNSY_COLUMN_FOLD;
}

#ifdef RNSY_COUNT_SLOW
long rnsy_slow_count = 0;
#endif

// sum over source channels ch of (digits of channel ch, from memory src)
// times the constant column col[ch] -> columns acc[0..S]
INL void vext(__m512i *acc, const uint64_t *src, const uint64_t *col, const int S)
{
    int ch, t;
    vzero(acc, S + 1);
    UNROLL
    for (ch = 0; ch < NCH; ch++)
    {
        __m512i k = LD(col + ch * NCH);
        UNROLL
        for (t = 0; t < S; t++)
        {
            __m512i bx = BC(src[t * NCH + ch]);
            acc[t] = LO52(acc[t], bx, k);
            acc[t + 1] = HI52(acc[t + 1], bx, k);
        }
    }
}

INL void vmulmod_s(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y,
    const int S, const int SQR)
{
    rnsy_fixed_t *fx = ctx->fx;
    __m512i a[MAXS], b[MAXS], t[MAXS], q[MAXS], acc[2 * MAXS + 2];
    __m512i db = LD(ctx->vdb), dc = LD(ctx->vdc);
    uint64_t *xi = ctx->vscratch, *xr = ctx->vscratch + S * NCH;
    uint64_t rx = x[16 * S], ry = y[16 * S], ur, sr, alpha;
    int k;

    // 1: B, xi = fold(fold(x y) f1)
    vload(a, x, S);
    if (SQR)
        vsqr(acc, a, S);
    else
    {
        vload(b, y, S);
        vzero(acc, 2 * S + 1);
        vmul_acc(acc, a, b, S);
    }
    vfold_prod(t, acc, db, S);
    if (S >= 2)
    {
        vzero(acc, S + 2);
        vcmul_acc(acc, t, ctx->vF1, S);
        vfold_cprod(t, acc, db, S);
    }
    else
    {
        vload(b, ctx->vf1, S);
        vzero(acc, 2 * S + 1);
        vmul_acc(acc, t, b, S);
        vfold_prod(t, acc, db, S);
    }
    vstore(xi, t, S);

    // redundant channel: U_r = T_r M_B^{-1} + sum_i lo52(xi_i) BrN_i
    ur = rx * ry * fx->MBinv_r +
        (uint64_t)_mm512_reduce_add_epi64(LO52(ZERO, t[0], LD(ctx->vBrN)));
    ur &= RNS_MASK52;

    // 2: first extension Q_j = fold(sum_i xi_i Bt_ij)
    vext(acc, xi, ctx->vbt, S);
    vnorm(acc, S + 1);
    vfold_top(acc, acc[S], dc, S);
    for (k = 0; k < S; k++)
        q[k] = acc[k];

    // 3: C, xi' = fold(fold(x' y') g1 + Q g2)
    vload(a, x + S * NCH, S);
    if (SQR)
        vsqr(acc, a, S);
    else
    {
        vload(b, y + S * NCH, S);
        vzero(acc, 2 * S + 1);
        vmul_acc(acc, a, b, S);
    }
    vfold_prod(t, acc, dc, S);
    if (S >= 2)
    {
        vzero(acc, S + 2);
        vcmul_acc(acc, t, ctx->vG1, S);
        vcmul_acc(acc, q, ctx->vG2, S);
        vfold_cprod(t, acc, dc, S);
    }
    else
    {
        vzero(acc, 2 * S + 1);
        vload(b, ctx->vg1, S);
        vmul_acc(acc, t, b, S);
        vload(b, ctx->vg2, S);
        vmul_acc(acc, q, b, S);
        vfold_prod(t, acc, dc, S);
    }
    vstore(xr, t, S);

    // alpha' = (sum_j lo52(xi'_j) |C_j|_r - U_r) M_C^{-1}  (mod 2^52)
    sr = (uint64_t)_mm512_reduce_add_epi64(LO52(ZERO, t[0], LD(ctx->vCr)));
    alpha = ((sr - ur) * fx->MCinv_r) & RNS_MASK52;
    vstore(z + S * NCH, t, S);

    // 5: U = fold(K - sum_j xi'_j Ct_ji + alpha' |-M_C|_b)
    vext(acc, xr, ctx->vct, S);
    vnorm(acc, S + 1);
    {
        __m512i borrow = ZERO, va = BC(alpha);
        UNROLL
        for (k = 0; k <= S; k++)
        {
            __m512i v = _mm512_add_epi64(_mm512_sub_epi64(LD(ctx->vK + k * NCH), acc[k]), borrow);
            acc[k] = _mm512_and_si512(v, M52);
            borrow = _mm512_srai_epi64(v, 52);
        }
        UNROLL
        for (k = 0; k < S; k++)
        {
            __m512i m = LD(ctx->vnegMC + k * NCH);
            acc[k] = LO52(acc[k], va, m);
            acc[k + 1] = HI52(acc[k + 1], va, m);
        }
    }
    vnorm(acc, S + 1);
    vfold_top(acc, acc[S], db, S);
    vstore(z, acc, S);
    z[16 * S] = ur;
}

INL void vaddsub_s(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y,
    const int S, const int SUB)
{
    int base, k;
    for (base = 0; base < 2; base++)
    {
        __m512i v[MAXS + 1], c = ZERO;
        __m512i d = LD(base ? ctx->vdc : ctx->vdb);
        const uint64_t *E = base ? ctx->vEc : ctx->vEb;
        const uint64_t *xp = x + base * S * NCH, *yp = y + base * S * NCH;

        if (SUB)
        {
            // x + (subconst + 2m) - y, signed digit chain, top digit <= 3
            UNROLL
            for (k = 0; k <= S; k++)
            {
                __m512i e = _mm512_add_epi64(LD(E + k * NCH), c);
                if (k < S)
                    e = _mm512_sub_epi64(_mm512_add_epi64(e, LD(xp + k * NCH)), LD(yp + k * NCH));
                v[k] = _mm512_and_si512(e, M52);
                c = _mm512_srai_epi64(e, 52);
            }
        }
        else
        {
            UNROLL
            for (k = 0; k < S; k++)
                v[k] = _mm512_add_epi64(LD(xp + k * NCH), LD(yp + k * NCH));
            v[S] = vnorm(v, S);
        }
        vfold_top(v, v[S], d, S);
        vstore(z + base * S * NCH, v, S);
    }
    if (SUB)
        z[16 * S] = (x[16 * S] + ctx->subconst[16 * S] - y[16 * S]) & RNS_MASK52;
    else
        z[16 * S] = (x[16 * S] + y[16 * S]) & RNS_MASK52;
}

// per-s instances and dispatch
#define INSTANCES(S) \
    static void mul_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { vmulmod_s(c, z, x, y, S, 0); } \
    static void sqr_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { vmulmod_s(c, z, x, y, S, 1); } \
    static void add_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { vaddsub_s(c, z, x, y, S, 0); } \
    static void sub_##S(rnsy_ctx_t *c, uint64_t *z, uint64_t *x, uint64_t *y) \
        { vaddsub_s(c, z, x, y, S, 1); }

INSTANCES(1) INSTANCES(2) INSTANCES(3) INSTANCES(4) INSTANCES(5)
INSTANCES(6) INSTANCES(7) INSTANCES(8) INSTANCES(9) INSTANCES(10)

typedef void (*opfn)(rnsy_ctx_t *, uint64_t *, uint64_t *, uint64_t *);
#define TABLE(op) static const opfn op##_tab[MAXS + 1] = { NULL, \
    op##_1, op##_2, op##_3, op##_4, op##_5, op##_6, op##_7, op##_8, op##_9, op##_10 }
TABLE(mul); TABLE(sqr); TABLE(add); TABLE(sub);

static int check_s(rnsy_ctx_t *ctx)
{
    if (ctx->fx->s < 1 || ctx->fx->s > MAXS || ctx->vblock == NULL)
    {
        fprintf(stderr, "rnsy: no one-number vector kernel for n = %d, s = %d\n",
            ctx->fx->n, ctx->fx->s);
        exit(1);
    }
    return ctx->fx->s;
}

void rnsy_vmulmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    mul_tab[check_s(ctx)](ctx, z, x, y);
}

void rnsy_vsqrmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x)
{
    sqr_tab[check_s(ctx)](ctx, z, x, x);
}

void rnsy_vaddmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    add_tab[check_s(ctx)](ctx, z, x, y);
}

void rnsy_vsubmod_avx512(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    sub_tab[check_s(ctx)](ctx, z, x, y);
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

void rnsy_vpowm_avx512(rnsy_ctx_t *ctx, uint64_t *d, uint64_t *b, mpz_t e)
{
    int S = check_s(ctx);
    opfn mul = mul_tab[S], sqr = sqr_tab[S];
    int bit = (int)mpz_sizeinbase(e, 2) - 1;
    int k = winsize(bit), ntab = 1 << (k - 1), i, j;
    size_t bytes = (size_t)ctx->fx->nwords * sizeof(uint64_t);
    uint64_t **g = (uint64_t **)calloc((size_t)ntab, sizeof(uint64_t *));
    uint64_t *b2 = rnsy_alloc(ctx);

    for (i = 0; i < ntab; i++)
        g[i] = rnsy_alloc(ctx);
    memcpy(g[0], b, bytes);
    sqr(ctx, b2, b, b);
    for (i = 1; i < ntab; i++)
        mul(ctx, g[i], g[i - 1], b2);

    memcpy(d, ctx->one, bytes);
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

    for (i = 0; i < ntab; i++)
        rns_vfree(g[i]);
    free(g);
    rns_vfree(b2);
}

#endif
