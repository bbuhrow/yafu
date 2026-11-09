/*
 * Shared AVX-512 helpers for the Yao-regime kernels (rns_yao_avx512.c,
 * rns_yaob_avx512.c): carry normalization, the pseudo-Mersenne folds and
 * fully unrolled channel products.  Each vector holds one 52-bit digit of
 * 8 independent channel values (8 channels of one number, or one channel
 * of 8 numbers); S (words per channel) is a compile-time constant at every
 * call site.
 */

#ifndef RNS_YAO_VEC_H
#define RNS_YAO_VEC_H

#include "rns_yao.h"
#include "rns_vec52.h"     // vector layer (RNS_VW), LO52 / HI52

#define NCH RNSY_NCH
#define MAXS 10
#define INL RNS_INLINE
#define UNROLL RNS_UNROLL(32)

#define BC(p)           V_BC(p)
#define LD(p)           V_LD(p)
#define ST(p, v)        V_ST((p), (v))
#define ZERO            V_ZERO
#define M52             BC(RNS_MASK52)

// RNSY_COLUMN_FOLD = 1: fold variable products (x*y in B and C) at column
// level, splitting each unnormalized high column into lo52 and its small
// overflow, so one carry chain replaces two.  The residues stay valid
// (congruent, < 2^w, the same bounds) but no longer match rns_yao_ref.c bit
// for bit; rns_yao_test switches to congruence/consistency/bound checks.
// RNSY_COLUMN_FOLD = 0: the exact kernel.
#ifndef RNSY_COLUMN_FOLD
#define RNSY_COLUMN_FOLD 1
#endif


#ifdef RNSY_COUNT_SLOW
extern long rnsy_slow_count;
#endif

// carry-normalize n columns in place; returns the carry out
INL rns_vi vnorm(rns_vi *a, int n)
{
    rns_vi c = ZERO;
    int t;
    UNROLL
    for (t = 0; t < n; t++)
    {
        rns_vi v = V_ADD(a[t], c);
        a[t] = V_AND(v, M52);
        c = V_SRLI(v, 52);
    }
    return c;
}

// x[0..S-1] + top*2^w, digits normalized, top < 2^52  ->  < 2^w, by up to
// three folds (enough for every input this kernel produces).  top*d < 2^59,
// so its high half is tiny and fold 1 carries past digit 1 only when digit 1
// is within a few units of 2^52: the fast path checks that and stops, since
// the value is then already < 2^w (exactly where the reference stops too).
INL void vfold_top_x(rns_vi *x, rns_vi top, rns_vi top1, rns_vi d,
    const int S, const int TWO)
{
    rns_vi c;

    if (S == 1)
    {
        x[0] = LO52(x[0], top, d);
        c = V_ADD(vnorm(x, 1), HI52(ZERO, top, d));
        if (TWO)
            c = LO52(c, top1, d);           // + top1*d*2^52; c stays < 2^21
        x[0] = LO52(x[0], c, d);
        c = vnorm(x, 1);
        x[0] = LO52(x[0], c, d);
        vnorm(x, 1);
        return;
    }

    // fold 1: + top*d into digits 0 and 1 (+ top1*d*2^52 into digit 1 when
    // the top spans two digits; top1*d is tiny, its high half is 0)
    x[0] = LO52(x[0], top, d);
    x[1] = V_ADD(HI52(x[1], top, d), V_SRLI(x[0], 52));
    if (TWO)
        x[1] = LO52(x[1], top1, d);
    x[0] = V_AND(x[0], M52);
    if (RNS_LIKELY(!V_ANY_GE52(x[1])))
        return;

    // rare: the carry ripples on, possibly out of the top digit
#ifdef RNSY_COUNT_SLOW
    rnsy_slow_count++;
#endif
    c = vnorm(x + 1, S - 1);
    // fold 2: c in {0, 1}
    x[0] = LO52(x[0], c, d);
    c = vnorm(x, S);
    // fold 3: c in {0, 1}; after it the value is < 2^w
    x[0] = LO52(x[0], c, d);
    vnorm(x, S);
}

INL void vfold_top(rns_vi *x, rns_vi top, rns_vi d, const int S)
{
    vfold_top_x(x, top, ZERO, d, S, 0);
}

// acc[0..S+1] += sum_i a_i * tab_i, tab_i = f * 2^(52i) mod m (S digits each)
INL void vcmul_acc(rns_vi *acc, const rns_vi *a, const uint64_t *tab, const int S)
{
    int i, k;
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (k = 0; k < S; k++)
        {
            rns_vi f = LD(tab + (i * S + k) * NCH);
            acc[k] = LO52(acc[k], a[i], f);
            acc[k + 1] = HI52(acc[k + 1], a[i], f);
        }
    }
}

// columns acc[0..S] of a constant product (S >= 2) -> x < 2^w
INL void vfold_cprod(rns_vi *x, rns_vi *acc, rns_vi d, const int S)
{
    int t;
    acc[S + 1] = vnorm(acc, S + 1);         // value < 2S * 2^52 * 2^w
    UNROLL
    for (t = 0; t < S; t++)
        x[t] = acc[t];
    vfold_top_x(x, acc[S], acc[S + 1], d, S, 1);
}

// columns acc[0..2S-1] (a sum of up to two S x S products) -> x < 2^w
INL void vfold_prod(rns_vi *x, rns_vi *acc, rns_vi d, const int S)
{
    rns_vi y[MAXS + 2];
    rns_vi top;
    int t;

#if RNSY_COLUMN_FOLD
    // Products only fill columns 0..2S-1, each < 2^58.  Fold at the column
    // boundary instead of the value boundary: low columns pass through
    // unnormalized, and each high column c folds as lo52(c)*d (IFMA lo/hi)
    // plus (c >> 52)*d (< 2^13) one digit up.  Y = sum < 2^(w+13), so a
    // single carry chain leaves a top digit < 2^13.
    UNROLL
    for (t = 0; t < S; t++)
        y[t] = acc[t];
    y[S] = ZERO;
#if !RNS_VEC_FP
    UNROLL
    for (t = 0; t < S; t++)
    {
        rns_vi c = acc[S + t];
        rns_vi lo = V_AND(c, M52);
        y[t] = LO52(y[t], lo, d);
        y[t + 1] = HI52(y[t + 1], lo, d);
        y[t + 1] = LO52(y[t + 1], V_SRLI(c, 52), d);
    }
#else
    {
        // FP: one lo/hi pair per column (biases removed below); (c >> 52)*d
        // < 2^32 * 2^32 is exact with a 32x32 -> 64 integer multiply
        const rns_vd dd = rns_cvt52(d);
        UNROLL
        for (t = 0; t < S; t++)
        {
            rns_vi c = acc[S + t];
            RNS_FP_MAC(y[t], y[t + 1], rns_cvt52(c), dd);
            y[t + 1] = V_ADD(y[t + 1],
                V_MULU32(V_SRLI(c, 52), d));
        }
        UNROLL
        for (t = 0; t <= S; t++)
            y[t] = rns_fp_unbias(y[t], t < S, t > 0);
    }
#endif
    vnorm(y, S + 1);
#else
    acc[2 * S] = vnorm(acc, 2 * S);         // carry out <= 1
    UNROLL
    for (t = 0; t < S; t++)
        y[t] = acc[t];
    y[S] = ZERO;
    // fold 1: + (acc[S..2S]) * d
#if !RNS_VEC_FP
    UNROLL
    for (t = 0; t <= S; t++)
    {
        y[t] = LO52(y[t], acc[S + t], d);
        if (t < S)
            y[t + 1] = HI52(y[t + 1], acc[S + t], d);
    }
#else
    {
        const rns_vd dd = rns_cvt52(d);
        UNROLL
        for (t = 0; t < S; t++)
            RNS_FP_MAC(y[t], y[t + 1], rns_cvt52(acc[S + t]), dd);
        UNROLL
        for (t = 0; t <= S; t++)
            y[t] = rns_fp_unbias(y[t], t < S, t > 0);
        // acc[2S] <= 1: its product with d is exact in 64 bits
        y[S] = V_ADD(y[S], V_MULU32(acc[2 * S], d));
    }
#endif
    vnorm(y, S + 1);                        // value < 2^(w+52): no carry out
#endif
    top = y[S];
    UNROLL
    for (t = 0; t < S; t++)
        x[t] = y[t];
    vfold_top(x, top, d, S);
}

#if RNS_VEC_FP
// FP products: column c of an S x S product gets L(c) lo and L(c-1) hi
// halves, L(c) = number of i + j == c
INL int rns_ncol(int c, const int S)
{
    if (c < 0 || c > 2 * S - 2)
        return 0;
    return (c < S) ? c + 1 : 2 * S - 1 - c;
}

// acc[0..2S-1] += a * b, schoolbook column sums
INL void vmul_acc(rns_vi *acc, const rns_vi *a, const rns_vi *b, const int S)
{
    rns_vd ad[MAXS], bd[MAXS];
    int i, j;
    UNROLL
    for (i = 0; i < S; i++)
    {
        ad[i] = rns_cvt52(a[i]);
        bd[i] = rns_cvt52(b[i]);
    }
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (j = 0; j < S; j++)
            RNS_FP_MAC(acc[i + j], acc[i + j + 1], ad[i], bd[j]);
    }
    UNROLL
    for (i = 0; i < 2 * S; i++)
        acc[i] = rns_fp_unbias(acc[i], rns_ncol(i, S), rns_ncol(i - 1, S));
}

// acc[0..2S-1] = a^2: off-diagonal products once, doubled, plus diagonal
INL void vsqr(rns_vi *acc, const rns_vi *a, const int S)
{
    rns_vd ad[MAXS];
    int i, j;
    UNROLL
    for (i = 0; i <= 2 * S; i++)
        acc[i] = ZERO;
    UNROLL
    for (i = 0; i < S; i++)
        ad[i] = rns_cvt52(a[i]);
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (j = i + 1; j < S; j++)
            RNS_FP_MAC(acc[i + j], acc[i + j + 1], ad[i], ad[j]);
    }
    // off-diagonal counts: lo halves at c = i + j (i < j), hi at c + 1
    UNROLL
    for (i = 1; i < 2 * S; i++)
    {
        int nlo = 0, nhi = 0, k;
        for (k = 0; k < S; k++)
        {
            if (i - k > k && i - k < S) nlo++;
            if (i - 1 - k > k && i - 1 - k < S) nhi++;
        }
        acc[i] = rns_fp_unbias(acc[i], nlo, nhi);
        acc[i] = V_ADD(acc[i], acc[i]);
    }
    UNROLL
    for (i = 0; i < S; i++)
        RNS_FP_MAC(acc[2 * i], acc[2 * i + 1], ad[i], ad[i]);
    UNROLL
    for (i = 0; i < 2 * S; i++)
        acc[i] = rns_fp_unbias(acc[i], !(i & 1), i & 1);   // lo at 2i, hi at 2i+1
}
#else
// acc[0..2S-1] += a * b, schoolbook column sums
INL void vmul_acc(rns_vi *acc, const rns_vi *a, const rns_vi *b, const int S)
{
    int i, j;
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (j = 0; j < S; j++)
        {
            acc[i + j] = LO52(acc[i + j], a[i], b[j]);
            acc[i + j + 1] = HI52(acc[i + j + 1], a[i], b[j]);
        }
    }
}

// acc[0..2S-1] = a^2: off-diagonal products once, doubled, plus diagonal
INL void vsqr(rns_vi *acc, const rns_vi *a, const int S)
{
    int i, j;
    UNROLL
    for (i = 0; i <= 2 * S; i++)
        acc[i] = ZERO;
    UNROLL
    for (i = 0; i < S; i++)
    {
        UNROLL
        for (j = i + 1; j < S; j++)
        {
            acc[i + j] = LO52(acc[i + j], a[i], a[j]);
            acc[i + j + 1] = HI52(acc[i + j + 1], a[i], a[j]);
        }
    }
    UNROLL
    for (i = 1; i < 2 * S; i++)
        acc[i] = V_ADD(acc[i], acc[i]);
    UNROLL
    for (i = 0; i < S; i++)
    {
        acc[2 * i] = LO52(acc[2 * i], a[i], a[i]);
        acc[2 * i + 1] = HI52(acc[2 * i + 1], a[i], a[i]);
    }
}
#endif

INL void vload(rns_vi *v, const uint64_t *p, int n)
{
    int t;
    UNROLL
    for (t = 0; t < n; t++)
        v[t] = LD(p + t * NCH);
}

INL void vstore(uint64_t *p, const rns_vi *v, int n)
{
    int t;
    UNROLL
    for (t = 0; t < n; t++)
        ST(p + t * NCH, v[t]);
}

INL void vzero(rns_vi *v, int n)
{
    int t;
    UNROLL
    for (t = 0; t < n; t++)
        v[t] = ZERO;
}

#endif
