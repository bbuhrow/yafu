/*
microecm_edw.c: microecm with stage 1 on twisted Edwards curves.

Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla (the parts
taken from microecm.c: Montgomery arithmetic, gcd, the x-only stage 2 and
its pairing tables).  Dual-licensed, at your option, under the FreeBSD
license or the Mozilla Public License 2.0, exactly as microecm.c; see the
full license texts there.

Design and measurements: arith/micro_tiny_edwards_plan.md.  In short:
  - Curves: GKL family (iv), -x^2 + y^2 = 1 + d x^2 y^2, d = -e^4.  The
    sequence [m]Q on E': y^2 = x^3 - x^2 - 9x + 9, Q = (5, 8), is stepped
    in projective coordinates mod n (one mixed addition per curve), and
    t, e, (x1, y1) and d come out as homogeneous polynomials in t = tn/td:
    no inversion anywhere, d stays a fraction dn/dd.  Family (iv) finds
    9-14% more factors per curve than Suyama's at microecm's sizes.
  - Stage 1: the odd part of microecm's scalar s on a precomputed
    signed-digit chain (edw_tiny_chains.h, from arith/edwards/
    tiny_chains.py) in extended coordinates (a = -1: doubling 3M + 4S,
    cached addition 7M).  The doubling's independent products overlap,
    which beats the latency-bound x-only PRAC chain even though it costs
    more multiplications.  Then the hand-off to the Montgomery image,
    (X : Z) = (Z + Y : Z - Y), A24 = dd / (dd + dn), and the power of 2 of
    s as x-only doublings.
  - Stage 2: microecm's pairing stage 2 (D = 60, B2 = 25 or 50 B1) with a
    projective A24 in its few doublings.
*/

#include "microecm_edw.h"
#include "edw_tiny_chains.h"
#include "edw_tiny_curves.h"
#include "edw_tiny_v2.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

#ifndef ALIGNED_MEM
#if defined(_MSC_VER)
#define ALIGNED_MEM __declspec(align(64))
#else
#define ALIGNED_MEM __attribute__((aligned(64)))
#endif
#endif

#if defined(_MSC_VER)
#define UEDW_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define UEDW_INLINE inline __attribute__((always_inline))
#else
#define UEDW_INLINE __inline
#endif

/* ---------------------------------------------------------------- arith */
/* Montgomery multiplication with the positive inverse (Hurchalla), R = 2^64,
   inputs and outputs in [0, n); modular add / sub; inverse mod 2^64; binary
   gcd.  One shared copy of each in arith/limb/limb1.h -- this file used to
   carry its own.  Scalar layer only: the 52-bit vector code below is local. */
#ifndef LIMB1_SCALAR_ONLY
#define LIMB1_SCALAR_ONLY
#endif
#include "limb1.h"

#define uedw_mulredc                mulredc_pos
#define uedw_submod                 submod
#define uedw_addmod                 addmod
#define uedw_multiplicative_inverse multiplicative_inverse
#define uedw_bingcd64               bingcd64

/* per-modulus constants, Montgomery form */
typedef struct {
    uint64_t n, rho, one, Rsqr;
    uint64_t c6, c10, c16, c32, c48;
} uedw_ctx;

#define MUL(a, b) uedw_mulredc((a), (b), c->n, c->rho)
#define SQR(a) uedw_mulredc((a), (a), c->n, c->rho)
#define ADD(a, b) uedw_addmod((a), (b), c->n)
#define SUB(a, b) uedw_submod((a), (b), c->n)

static void uedw_ctx_init(uedw_ctx *c, uint64_t n)
{
    c->n = n;
    c->rho = uedw_multiplicative_inverse(n);
    c->one = ((uint64_t)0 - n) % n;            // R mod n
    uint64_t c2 = ADD(c->one, c->one);
    uint64_t c4 = ADD(c2, c2);
    uint64_t c8 = ADD(c4, c4);
    c->c16 = ADD(c8, c8);
    c->c32 = ADD(c->c16, c->c16);
    c->c48 = ADD(c->c32, c->c16);
    c->c6 = ADD(c4, c2);
    c->c10 = ADD(c8, c2);
    uint64_t t = SQR(c->c16);                   // R 2^8
    t = SQR(t);                                 // R 2^16
    t = SQR(t);                                 // R 2^32
    c->Rsqr = SQR(t);                           // R 2^64 = R^2
}

/* --------------------------------------------- GKL family (iv), mod n */
/* E': y^2 = x^3 - x^2 - 9x + 9, projective (X : Y : Z), Montgomery form */
typedef struct { uint64_t X, Y, Z; } uedw_wpt;

/* R + Q, Q = (5, 8) affine, R != +-Q.  lambda = u/v with u = 8Z - Y,
   v = 5Z - X; A = u^2 Z - v^2 (X + (a2 + 5) Z), a2 + 5 = 4;
   X3 = v A, Y3 = u (X v^2 - A) - Y v^3, Z3 = v^3 Z. */
static void uedw_w_addq(uedw_wpt *R, const uedw_ctx *c)
{
    uint64_t X = R->X, Y = R->Y, Z = R->Z;
    uint64_t z2 = ADD(Z, Z);
    uint64_t z4 = ADD(z2, z2);
    uint64_t u = SUB(ADD(z4, z4), Y);
    uint64_t v = SUB(ADD(z4, Z), X);
    uint64_t u2 = SQR(u);
    uint64_t v2 = SQR(v);
    uint64_t v3 = MUL(v2, v);
    uint64_t w1 = MUL(v2, X);
    uint64_t w2 = MUL(v2, Z);
    w2 = ADD(w2, w2);
    w2 = ADD(w2, w2);
    uint64_t A = SUB(SUB(MUL(u2, Z), w1), w2);
    R->X = MUL(v, A);
    R->Y = SUB(MUL(u, SUB(w1, A)), MUL(Y, v3));
    R->Z = MUL(v3, Z);
}

/* 2R.  lambda = w/s with w = 3X^2 - 2XZ - 9Z^2, s = 2YZ;
   A = w^2 Z - s^2 (2X - Z); X3 = s A, Y3 = w (X s^2 - A) - Y s^3,
   Z3 = s^3 Z. */
static void uedw_w_dbl(uedw_wpt *R, const uedw_ctx *c)
{
    uint64_t X = R->X, Y = R->Y, Z = R->Z;
    uint64_t X2 = SQR(X), Z2 = SQR(Z), XZ = MUL(X, Z);
    uint64_t w = ADD(ADD(X2, X2), X2);
    w = SUB(w, ADD(XZ, XZ));
    uint64_t z8 = ADD(Z2, Z2);
    z8 = ADD(z8, z8);
    z8 = ADD(z8, z8);
    w = SUB(w, ADD(z8, Z2));
    uint64_t s = MUL(Y, Z);
    s = ADD(s, s);
    uint64_t s2 = SQR(s);
    uint64_t s3 = MUL(s2, s);
    uint64_t A = SUB(MUL(SQR(w), Z), MUL(s2, SUB(ADD(X, X), Z)));
    uint64_t Xs2 = MUL(X, s2);
    R->X = MUL(s, A);
    R->Y = SUB(MUL(w, SUB(Xs2, A)), MUL(Y, s3));
    R->Z = MUL(s3, Z);
}

/* [m]Q, m >= 1 */
static void uedw_w_start(uedw_wpt *R, uint32_t m, const uedw_ctx *c)
{
    uint64_t c4 = ADD(ADD(c->one, c->one), ADD(c->one, c->one));
    R->X = ADD(c4, c->one);                     // 5
    R->Y = ADD(c4, c4);                         // 8
    R->Z = c->one;
    int top = 31;
    while (top > 0 && !((m >> top) & 1))
        top--;
    for (int i = top - 1; i >= 0; i--) {
        uedw_w_dbl(R, c);
        if ((m >> i) & 1)
            uedw_w_addq(R, c);
    }
}

typedef struct { uint64_t X, Y, Z, T; } uedw_ept;
typedef struct { uint64_t YmX, YpX, T2d, Z2; } uedw_cached;
typedef struct {
    uedw_ept P;                                 // start point, extended
    uint64_t dn2, dd, dd2;                      // 2 dn, dd, 2 dd (d = dn/dd)
    uint64_t a24n, a24d;                        // A24 = dd / (dd + dn)
} uedw_curve;

/* The family (iv) curve of the E' point R (GKL Thm 3.5, j = 4):
   t = (4x + 4)/(y - 4) = tn/td,
   e = t(t + 4)/(t^2 - 4),
   x1 = 2(t + 1)(t - 2)(t + 2) / (t (t + 4)(t^2 + 2t + 4)),
   y1 = (t^6 + 6t^5 + 10t^4 - 16t^3 - 48t^2 - 32t - 32) /
        (t (t + 4)(t^4 + 2t^3 + 2t^2 + 8t + 16)),
   all homogenized in (a, b) = (tn, td); d = -e^4; the start point in
   extended coordinates (x1n y1d : y1n x1d : x1d y1d : x1n y1n).
   Checked against arith/edwards/edw_ref.py curve_iv for m = 2..39. */
static void uedw_gkl_curve(const uedw_wpt *R, uedw_curve *cv, const uedw_ctx *c)
{
    uint64_t t = ADD(R->X, R->Z);
    t = ADD(t, t);
    uint64_t a = ADD(t, t);                     // 4 (X + Z)
    uint64_t z4 = ADD(R->Z, R->Z);
    z4 = ADD(z4, z4);
    uint64_t b = SUB(R->Y, z4);                 // Y - 4 Z
    uint64_t a2 = SQR(a), b2 = SQR(b), ab = MUL(a, b);

    uint64_t b2x = ADD(b, b);
    uint64_t xn = MUL(MUL(SUB(a, b2x), ADD(a, b2x)), ADD(a, b));
    xn = MUL(xn, b);
    xn = ADD(xn, xn);
    uint64_t aa4 = MUL(a, ADD(a, ADD(b2x, b2x)));     // a (a + 4b)
    uint64_t b2_4 = ADD(b2, b2);
    b2_4 = ADD(b2_4, b2_4);
    uint64_t xd = MUL(aa4, ADD(ADD(a2, ADD(ab, ab)), b2_4));

    uint64_t t1 = SQR(a2), t2 = MUL(a2, ab), t3 = MUL(a2, b2), t4 = MUL(ab, b2), t5 = SQR(b2);
    uint64_t q = ADD(t2, t3);
    q = ADD(t1, ADD(q, q));
    uint64_t r = ADD(t4, ADD(t5, t5));          // (t4 + 2 t5) * 8
    r = ADD(r, r);
    r = ADD(r, r);
    r = ADD(r, r);
    uint64_t yd = MUL(aa4, ADD(q, r));

    uint64_t b3 = MUL(b2, b), b4 = SQR(b2), b5 = MUL(b4, b), b6 = SQR(b3);
    uint64_t h = ADD(a, MUL(b, c->c6));
    h = ADD(MUL(h, a), MUL(b2, c->c10));
    h = SUB(MUL(h, a), MUL(b3, c->c16));
    h = SUB(MUL(h, a), MUL(b4, c->c48));
    h = SUB(MUL(h, a), MUL(b5, c->c32));
    uint64_t yn = SUB(MUL(h, a), MUL(b6, c->c32));

    uint64_t en2 = SQR(aa4), ed = SUB(a2, b2_4), ed2 = SQR(ed);
    uint64_t dn = SUB(0, SQR(en2));
    uint64_t dd = SQR(ed2);

    cv->P.X = MUL(xn, yd);
    cv->P.Y = MUL(yn, xd);
    cv->P.Z = MUL(xd, yd);
    cv->P.T = MUL(xn, yn);
    cv->dn2 = ADD(dn, dn);
    cv->dd = dd;
    cv->dd2 = ADD(dd, dd);
    cv->a24n = dd;
    cv->a24d = ADD(dd, dn);
}

/* ------------------------------------------------- Edwards arithmetic */
/* doubling, a = -1, signs folded: 3M + 4S (+1M for T) */
static UEDW_INLINE void uedw_dbl(uedw_ept *R, const uedw_ept *P, int want_t, const uedw_ctx *c)
{
    uint64_t A = SQR(P->X);
    uint64_t B = SQR(P->Y);
    uint64_t C = SQR(P->Z);
    uint64_t S = SQR(ADD(P->X, P->Y));
    C = ADD(C, C);
    uint64_t H = ADD(A, B);
    uint64_t E = SUB(H, S);
    uint64_t G = SUB(A, B);
    uint64_t F = ADD(C, G);
    R->X = MUL(E, F);
    R->Y = MUL(G, H);
    R->Z = MUL(F, G);
    if (want_t)
        R->T = MUL(E, H);
}

/* P +- Q with Q cached (unified add-2008-hwcd-3): 7M (+1M for T) */
static UEDW_INLINE void uedw_add(uedw_ept *R, const uedw_ept *P, const uedw_cached *Q,
    int neg, int want_t, const uedw_ctx *c)
{
    uint64_t A = MUL(SUB(P->Y, P->X), neg ? Q->YpX : Q->YmX);
    uint64_t B = MUL(ADD(P->Y, P->X), neg ? Q->YmX : Q->YpX);
    uint64_t C = MUL(P->T, Q->T2d);
    uint64_t D = MUL(P->Z, Q->Z2);
    uint64_t E = SUB(B, A);
    uint64_t H = ADD(B, A);
    uint64_t F = neg ? ADD(D, C) : SUB(D, C);
    uint64_t G = neg ? SUB(D, C) : ADD(D, C);
    R->X = MUL(E, F);
    R->Y = MUL(G, H);
    R->Z = MUL(F, G);
    if (want_t)
        R->T = MUL(E, H);
}

/* cached form, scaled by dd so that d = dn/dd needs no inversion: 4M */
static UEDW_INLINE void uedw_cache(uedw_cached *q, const uedw_ept *P, const uedw_curve *cv,
    const uedw_ctx *c)
{
    q->YmX = MUL(SUB(P->Y, P->X), cv->dd);
    q->YpX = MUL(ADD(P->Y, P->X), cv->dd);
    q->T2d = MUL(P->T, cv->dn2);
    q->Z2 = MUL(P->Z, cv->dd2);
}

/* ------------------------------------------- Montgomery x-only (stage 2) */
typedef struct { uint64_t X, Z; } uedw_pt;

/* 2 P with projective A24 = a24n / a24d: 2S + 4M */
static UEDW_INLINE void uedw_mdup(const uedw_pt *P, uedw_pt *R, const uedw_curve *cv,
    const uedw_ctx *c)
{
    uint64_t U = SQR(SUB(P->X, P->Z));
    uint64_t V = SQR(ADD(P->X, P->Z));
    uint64_t w = SUB(V, U);
    uint64_t Ub = MUL(U, cv->a24d);
    R->X = MUL(Ub, V);
    R->Z = MUL(w, ADD(Ub, MUL(w, cv->a24n)));
}

/* differential addition P1 + P2 given Pin = P1 - P2 (microecm's uecm_uadd) */
static UEDW_INLINE void uedw_uadd(const uedw_pt P1, const uedw_pt P2, const uedw_pt Pin,
    uedw_pt *Pout, const uedw_ctx *c)
{
    uint64_t diff1 = SUB(P1.X, P1.Z);
    uint64_t sum1 = ADD(P1.X, P1.Z);
    uint64_t diff2 = SUB(P2.X, P2.Z);
    uint64_t sum2 = ADD(P2.X, P2.Z);
    uint64_t tt1 = MUL(diff1, sum2);
    uint64_t tt2 = MUL(sum1, diff2);
    uint64_t tt3 = ADD(tt1, tt2);
    uint64_t tt4 = SUB(tt1, tt2);
    tt1 = SQR(tt3);
    tt2 = SQR(tt4);
    uint64_t tmpx = MUL(tt1, Pin.Z);
    uint64_t tmpz = MUL(tt2, Pin.X);
    Pout->X = tmpx;
    Pout->Z = tmpz;
}

/* ------------------------------------------------------------ stage 1 */
typedef struct {
    const int8_t *ch;
    int len, m, top, tz;
} uedw_chain;

static int uedw_chain_for(uint32_t B1, uedw_chain *c)
{
#define UEDW_CH(b) case b: c->ch = edw_chain_##b; c->len = EDW_CHAIN_##b##_LEN; \
    c->m = EDW_CHAIN_##b##_M; c->top = EDW_CHAIN_##b##_TOP; c->tz = EDW_CHAIN_##b##_TZ; return 1;
    switch (B1) {
        UEDW_CH(27) UEDW_CH(47) UEDW_CH(70) UEDW_CH(85) UEDW_CH(125) UEDW_CH(165) UEDW_CH(205)
    default:
        return 0;
    }
#undef UEDW_CH
}

/* [s]P for microecm's stage-1 scalar s at this B1; returns the Montgomery
   image (X : Z) of the result. */
static void uedw_stage1(const uedw_curve *cv, const uedw_chain *ch, uedw_pt *Pm,
    const uedw_ctx *c)
{
    uedw_cached tab[16];
    uedw_ept ext[16];
    int k, i;
    ext[0] = cv->P;
    uedw_cache(&tab[0], &ext[0], cv, c);
    if (ch->m > 1) {
        uedw_ept P2;
        uedw_cached c2;
        uedw_dbl(&P2, &cv->P, 1, c);
        uedw_cache(&c2, &P2, cv, c);
        for (k = 1; k <= (ch->m - 1) / 2; k++) {
            uedw_add(&ext[k], &ext[k - 1], &c2, 0, 1, c);
            uedw_cache(&tab[k], &ext[k], cv, c);
        }
    }
    uedw_ept R = ext[(ch->top - 1) / 2];
    for (i = 0; i < ch->len; i++) {
        int dg = ch->ch[i];
        if (dg == 0) {
            // T is needed only when an addition follows
            int t = (i + 1 < ch->len) && (ch->ch[i + 1] != 0);
            uedw_dbl(&R, &R, t, c);
        } else if (dg > 0) {
            uedw_add(&R, &R, &tab[(dg - 1) >> 1], 0, 0, c);
        } else {
            uedw_add(&R, &R, &tab[(-dg - 1) >> 1], 1, 0, c);
        }
    }
    // Montgomery image: u = (1 + y)/(1 - y), then the power of 2 in s
    Pm->X = ADD(R.Z, R.Y);
    Pm->Z = SUB(R.Z, R.Y);
    for (i = 0; i < ch->tz; i++)
        uedw_mdup(Pm, Pm, cv, c);
}

/* ------------------------------------------------------------ stage 2 */
static const int uedw_map[60] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

/* microecm.c's pre-paired stage-2 sequences (B2 = 25 B1; 50 B1 at B1 >= 125) */
static const int uedw_numb1_27 = 76;
static const uint8_t uedw_b1_27[76] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7 };
static const int uedw_numb1_47 = 121;
static const uint8_t uedw_b1_47[121] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7,1,19,37,47,
0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,
37,23,19,13,7,1,29,31,41,59,0,59,49,47,41,23,
19,17,13,7,1,43,53,0,59 };
static const int uedw_numb1_70 = 175;
static const uint8_t uedw_b1_70[175] = {
31,41,43,47,49,0,59,47,43,41,37,31,29,19,13,7,
1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,
0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,
41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,
31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,
29,31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,
53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,
53,49,47,43,31,29,23,11,17,0,47,43,41,37,31,23,
19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,
7,17,29,0,53,47,43,41,17,13,11,1,23,31,37 };
static const int uedw_numb1_85 = 225;
static const uint8_t uedw_b1_85[225] = {
    1,53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
    31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
    0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
    31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
    29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
    17,13,11,1,23,31,37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,
    43,41,31,17,7,1,11,13,19,29 };
static const int uedw_numb1_125_x50 = 645;
static const uint8_t uedw_b1_125_x50[645] = {
    23,19,13,11,1,7,17,29,31,41,43,47,49,0,59,47,
    43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,
    37,31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,
    23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,
    17,13,11,1,47,0,59,49,41,31,23,17,11,7,1,19,
    37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,
    49,43,37,23,19,13,7,1,29,31,41,59,0,59,49,47,
    41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,
    13,7,1,19,47,53,0,59,53,49,47,43,31,29,23,11,
    17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,
    59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
    17,13,11,1,23,31,37,49,0,53,47,43,41,29,19,7,
    1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,
    53,0,59,53,43,41,31,17,7,1,11,13,19,29,0,59,
    53,49,47,37,29,11,13,17,23,31,0,59,43,41,37,29,
    23,17,13,1,31,47,0,59,53,49,47,41,37,31,19,13,
    7,11,17,29,43,0,47,29,19,11,7,1,41,43,59,0,
    53,49,37,23,13,11,7,1,17,19,29,41,43,59,0,59,
    49,41,37,23,13,1,7,11,29,43,47,53,0,59,53,49,
    31,23,13,7,1,17,29,43,47,0,59,31,29,19,11,7,
    37,49,53,0,41,31,29,13,17,19,37,53,59,0,53,49,
    41,19,17,13,11,1,29,31,37,43,59,0,59,47,43,31,
    29,19,17,1,23,0,53,49,47,43,41,19,11,1,7,17,
    23,29,31,37,0,59,49,47,37,23,17,13,7,1,29,41,
    43,0,59,53,49,41,31,23,17,11,19,29,43,47,0,49,
    47,37,23,19,17,7,11,41,53,59,0,59,47,43,23,1,
    11,13,17,29,31,37,0,59,53,37,31,19,17,11,1,23,
    29,43,47,49,0,53,29,19,13,11,1,23,31,41,43,59,
    0,53,23,13,11,1,7,41,47,59,0,49,47,29,23,19,
    13,7,11,37,43,53,0,43,41,29,23,7,1,13,31,47,
    49,53,59,0,59,53,47,37,31,19,7,11,13,23,41,49,
    0,47,43,41,37,31,29,17,13,1,59,0,47,41,31,19,
    17,13,7,1,23,37,43,49,53,59,0,59,49,47,31,19,
    11,7,17,29,37,43,0,43,37,23,19,11,1,7,41,47,
    0,49,37,31,29,13,7,1,23,41,0,49,43,37,31,23,
    17,13,11,1,7,19,41,47,53,0,49,47,43,31,29,23,
    19,7,13,41,0,47,31,29,23,19,17,13,11,1,41,0,
    53,49,47,41,11,7,1,17,23,31,37,59,0,59,53,49,
    43,41,31,29,13,7,17,0,49,43,23,19,1,7,13,17,
    37,59,0,59,49,37,29,19,17,7,1,13,23,47,53,0,
    59,53,49,41,23 };
static const int uedw_numb1_165_x50 = 830;
static const uint8_t uedw_b1_165_x50[830] = {
    13,7,1,11,17,19,31,43,47,49,53,59,0,59,49,43,
    37,31,29,23,19,17,7,11,13,47,53,0,53,47,41,37,
    31,23,19,11,1,13,29,43,59,0,53,49,41,37,31,19,
    17,1,7,23,29,47,59,0,59,53,47,43,41,29,19,17,
    13,7,1,23,31,49,0,53,47,41,37,29,23,19,11,7,
    17,31,43,49,59,0,47,43,41,37,23,19,17,13,7,11,
    29,53,0,53,49,43,37,29,23,11,7,1,13,19,31,41,
    0,53,49,47,43,37,31,23,17,11,13,41,0,59,47,43,
    37,31,29,23,11,1,17,19,41,0,59,53,19,13,7,1,
    29,43,47,49,0,53,49,47,41,29,19,17,13,11,7,1,
    23,31,43,59,0,53,49,41,37,23,19,13,11,7,1,17,
    43,47,0,47,43,41,31,19,17,7,1,13,37,49,0,59,
    49,37,29,13,1,7,11,17,19,41,47,53,0,49,47,31,
    29,7,1,13,17,19,23,37,59,0,47,37,31,19,17,13,
    11,1,29,41,43,53,0,59,41,17,13,7,1,19,23,31,
    47,49,53,0,59,53,47,43,31,29,7,1,11,17,37,41,
    49,0,49,43,37,23,19,13,1,7,17,0,59,49,41,37,
    31,29,23,1,11,13,53,0,53,43,41,37,29,23,17,13,
    11,7,1,19,31,49,0,53,43,31,29,23,19,17,1,13,
    37,41,59,0,53,43,37,31,23,13,1,17,29,59,0,59,
    49,41,37,23,19,11,1,7,29,0,59,43,17,13,11,1,
    7,23,29,37,41,49,0,49,47,43,41,29,1,7,13,19,
    23,31,59,0,59,49,47,31,29,13,7,37,41,43,0,49,
    41,29,23,13,11,7,1,17,19,31,43,53,0,53,47,43,
    37,29,23,17,1,11,13,31,41,49,59,0,53,47,41,19,
    13,11,1,17,23,43,0,53,49,47,37,23,19,11,7,17,
    29,31,43,0,53,31,19,17,13,7,1,29,37,59,0,49,
    47,41,29,13,11,7,1,17,19,37,0,59,49,43,41,31,
    29,19,17,7,1,11,13,23,37,0,53,43,41,31,23,17,
    7,11,29,0,59,53,49,43,37,19,17,7,13,23,47,0,
    59,53,37,29,23,17,1,19,31,43,0,49,37,19,17,11,
    7,43,47,53,59,0,59,47,43,29,1,11,17,49,0,49,
    47,43,37,29,23,13,11,7,19,31,41,59,0,49,41,23,
    19,13,1,7,47,53,0,53,49,41,31,23,11,7,13,17,
    59,0,59,43,37,31,17,7,11,41,47,53,0,53,47,43,
    41,29,23,19,17,11,59,0,59,53,49,23,17,11,7,1,
    43,0,59,53,49,47,43,41,31,17,11,7,1,37,0,41,
    37,29,19,13,7,1,23,31,47,49,59,0,59,43,37,17,
    13,1,41,47,0,53,49,31,23,17,13,7,19,29,41,0,
    59,49,47,37,29,17,7,19,23,31,41,0,53,43,37,31,
    29,23,13,1,11,17,59,0,59,53,47,41,31,23,1,7,
    29,49,0,59,49,19,11,7,13,23,29,31,37,41,0,53,
    41,23,7,1,13,19,29,31,43,49,59,0,47,43,19,17,
    1,11,13,23,49,53,0,59,43,37,31,29,17,1,7,11,
    47,49,0,59,53,49,43,37,29,23,19,7,1,0,37,31,
    19,13,11,47,53,0,53,49,47,41,31,23,17,13,7,37,
    0,59,49,47,31,29,11,13,37,53,0,49,43,41,23,19,
    13,11,1,7,17,29,37,47,59,0,59,47,43,37,31,29,
    17,13,1,19,23,49,53,0,59,53,49,41,37,23,17,13,
    1,19,0,43,37,31,19,7,13,17,23,41,47,59,0,53,
    47,43,31,29,17,13,37,59,0,47,41,31,19,13,11,7,
    1,17,23,0,59,53,49,41,29,11,1,13,17,23 };
static const int uedw_numb1_205_x50 = 1013;
static const uint8_t uedw_b1_205_x50[1013] = {
    1,13,17,19,23,29,31,41,47,53,59,0,59,53,49,47,
    37,23,19,17,13,1,7,29,43,0,53,49,41,31,29,19,
    17,11,7,1,13,37,59,0,49,47,29,23,13,7,1,17,
    31,37,43,0,59,49,47,43,37,31,29,17,13,7,1,11,
    19,53,0,59,53,49,41,37,23,13,1,11,17,19,29,43,
    47,0,53,49,47,43,23,19,11,1,7,17,37,41,0,59,
    53,41,37,31,29,19,17,11,1,13,43,47,0,53,47,41,
    19,17,7,1,11,23,31,43,59,0,59,53,41,31,13,11,
    7,1,17,29,37,0,49,43,37,29,11,1,13,17,19,23,
    41,0,59,49,47,43,41,37,31,19,7,1,13,23,29,53,
    0,53,49,43,41,37,31,29,23,13,7,17,19,47,59,0,
    49,47,37,29,23,17,11,7,13,19,31,41,53,0,59,43,
    29,23,19,17,13,11,1,41,0,59,37,31,23,17,13,11,
    7,1,19,29,43,53,0,49,47,43,41,31,19,17,1,7,
    11,13,23,0,47,43,37,29,13,11,7,1,17,19,23,31,
    59,0,59,37,31,29,23,19,13,1,7,11,41,47,53,0,
    53,49,43,31,23,17,13,41,59,0,59,53,31,19,17,1,
    7,11,23,37,47,49,0,59,53,47,43,41,37,31,23,19,
    17,11,1,0,59,53,49,47,31,17,13,7,1,11,29,37,
    0,53,43,31,17,13,7,1,29,41,49,0,53,49,41,29,
    23,11,7,1,19,31,47,0,47,43,41,29,23,19,7,1,
    11,49,0,59,31,29,23,17,11,7,1,13,41,43,0,59,
    43,37,17,1,7,11,13,19,41,49,0,59,53,43,41,37,
    31,29,23,13,11,1,47,0,59,53,47,31,19,17,13,1,
    7,11,29,37,43,49,0,49,43,41,31,17,13,7,11,23,
    37,53,0,53,49,41,23,19,13,11,7,1,17,37,59,0,
    49,47,43,37,31,29,23,1,7,41,0,59,43,41,37,31,
    17,13,11,7,47,49,0,59,49,47,37,31,29,19,17,7,
    1,0,53,47,37,19,13,1,11,31,41,0,49,47,37,23,
    17,13,11,7,19,31,53,0,59,53,47,29,13,11,7,1,
    23,41,0,49,47,41,37,19,11,13,17,23,29,31,43,0,
    59,29,19,13,1,41,43,47,53,0,59,53,43,41,37,23,
    17,11,7,1,13,29,49,0,53,49,43,31,29,23,17,11,
    37,41,59,0,53,41,23,19,17,13,11,29,31,47,59,0,
    47,37,23,19,11,17,29,43,49,0,59,53,49,47,41,19,
    13,11,7,17,29,31,37,0,53,47,41,37,29,19,13,31,
    43,49,59,0,47,41,37,29,19,13,7,11,49,53,0,59,
    49,43,37,29,23,11,7,1,17,19,31,47,53,0,47,43,
    31,17,11,37,41,59,0,53,47,43,37,23,17,11,1,31,
    41,0,59,47,37,13,11,7,1,19,53,0,59,53,43,31,
    29,19,13,7,1,23,37,49,0,53,29,23,1,19,31,41,
    0,49,41,23,19,17,7,1,11,29,37,0,53,37,31,29,
    17,11,1,13,19,43,47,0,49,47,31,29,19,17,7,13,
    23,53,59,0,59,47,31,23,19,13,17,29,37,41,53,0,
    59,53,49,37,31,23,11,7,19,29,0,49,43,41,19,11,
    7,17,23,37,59,0,53,47,43,37,7,17,19,31,41,59,
    0,59,41,17,1,7,23,47,49,0,53,49,43,41,31,23,
    13,7,1,11,17,19,29,47,59,0,59,47,43,29,11,7,
    1,19,23,31,37,41,49,53,0,53,47,43,29,17,13,11,
    19,23,59,0,49,37,23,17,13,11,7,29,43,47,59,0,
    59,47,17,1,7,29,43,49,0,49,43,41,37,29,19,13,
    7,17,31,0,59,41,31,29,19,17,13,7,23,37,43,47,
    0,59,53,41,17,7,1,19,49,0,59,47,43,29,23,11,
    31,37,49,53,0,47,37,29,13,11,1,17,19,31,53,59,
    0,53,49,41,37,31,23,17,11,1,7,0,47,43,31,29,
    19,13,11,1,17,37,0,47,41,37,29,19,7,1,31,43,
    59,0,49,47,41,31,23,1,13,19,37,43,0,59,53,49,
    37,29,23,11,7,1,17,31,47,0,53,49,47,37,19,11,
    7,13,41,0,59,53,47,37,31,29,19,17,13,11,23,41,
    0,59,49,37,31,23,19,17,43,53,0,59,47,41,29,13,
    11,1,7,31,43,49,53,0,43,41,29,23,19,7,1,47,
    49,0,59,47,43,29,23,7,1,11,19,37,0,43,41,13,
    11,17,19,29,49,53,0,59,37,31,29,19,11,7,1,23,
    41,53,0,47,43 };

#if defined(USE_AVX512F) || defined(USE_AVX2)
/* the B2 = 25 B1 sequences at B1 >= 125, used by the vector stage 2 (as microecm's x8) */
static const int uedw_numb1_125 = 319;
static const uint8_t uedw_b1_125[319] = {
    23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
    29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,41,31,23,
    17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,31,41,59,0,59,49,
    47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,29,23,11,17,0,47,
    43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,17,13,11,1,23,31,
    37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,43,41,31,17,7,1,
    11,13,19,29,0,59,53,49,47,37,29,11,13,17,23,31,0,59,43,41,37,29,23,17,13,1,31,47,0,59,53,49,47,41,37,
    31,19,13,7,11,17,29,43,0,47,29,19,11,7,1,41,43,59,0,53,49,37,23,13,11,7,1,17,19,29,41,43,59,0,59,49,
    41,37,23,13,1,7,11,29,43,47,53,0,59,53,49,31,23,13,7,1,17,29,43,47,0,59,31,29,19,11,7,37,49,53 };
static const int uedw_numb1_165 = 435;
static const uint8_t uedw_b1_165[435] = {
    13,7,1,11,17,19,31,43,47,49,53,59,0,59,49,43,
    37,31,29,23,19,17,7,11,13,47,53,0,53,47,41,37,
    31,23,19,11,1,13,29,43,59,0,53,49,41,37,31,19,
    17,1,7,23,29,47,59,0,59,53,47,43,41,29,19,17,
    13,7,1,23,31,49,0,53,47,41,37,29,23,19,11,7,
    17,31,43,49,59,0,47,43,41,37,23,19,17,13,7,11,
    29,53,0,53,49,43,37,29,23,11,7,1,13,19,31,41,
    0,53,49,47,43,37,31,23,17,11,13,41,0,59,47,43,
    37,31,29,23,11,1,17,19,41,0,59,53,19,13,7,1,
    29,43,47,49,0,53,49,47,41,29,19,17,13,11,7,1,
    23,31,43,59,0,53,49,41,37,23,19,13,11,7,1,17,
    43,47,0,47,43,41,31,19,17,7,1,13,37,49,0,59,
    49,37,29,13,1,7,11,17,19,41,47,53,0,49,47,31,
    29,7,1,13,17,19,23,37,59,0,47,37,31,19,17,13,
    11,1,29,41,43,53,0,59,41,17,13,7,1,19,23,31,
    47,49,53,0,59,53,47,43,31,29,7,1,11,17,37,41,
    49,0,49,43,37,23,19,13,1,7,17,0,59,49,41,37,
    31,29,23,1,11,13,53,0,53,43,41,37,29,23,17,13,
    11,7,1,19,31,49,0,53,43,31,29,23,19,17,1,13,
    37,41,59,0,53,43,37,31,23,13,1,17,29,59,0,59,
    49,41,37,23,19,11,1,7,29,0,59,43,17,13,11,1,
    7,23,29,37,41,49,0,49,47,43,41,29,1,7,13,19,
    23,31,59,0,59,49,47,31,29,13,7,37,41,43,0,49,
    41,29,23,13,11,7,1,17,19,31,43,53,0,53,47,43,
    37,29,23,17,1,11,13,31,41,49,59,0,53,47,41,19,
    13,11,1,17,23,43,0,53,49,47,37,23,19,11,7,17,
    29,31,43,0,53,31,19,17,13,7,1,29,37,59,0,49,
    47,41,29 };
static const int uedw_numb1_205 = 511;
static const uint8_t uedw_b1_205[511] = {
    1,23,41,0,59,53,49,47,37,23,19,17,13,1,7,29,43,0,53,49,41,31,29,19,17,11,7,1,13,37,59,0,49,47,29,23,
    13,7,1,17,31,37,43,0,59,49,47,43,37,31,29,17,13,7,1,11,19,53,0,59,53,49,41,37,23,13,1,11,17,19,29,43,
    47,0,53,49,47,43,23,19,11,1,7,17,37,41,0,59,53,41,37,31,29,19,17,11,1,13,43,47,0,53,47,41,19,17,7,1,
    11,23,31,43,59,0,59,53,41,31,13,11,7,1,17,29,37,0,49,43,37,29,11,1,13,17,19,23,41,0,59,49,47,43,41,37,
    31,19,7,1,13,23,29,53,0,53,49,43,41,37,31,29,23,13,7,17,19,47,59,0,49,47,37,29,23,17,11,7,13,19,31,41,
    53,0,59,43,29,23,19,17,13,11,1,41,0,59,37,31,23,17,13,11,7,1,19,29,43,53,0,49,47,43,41,31,19,17,1,7,11,
    13,23,0,47,43,37,29,13,11,7,1,17,19,23,31,59,0,59,37,31,29,23,19,13,1,7,11,41,47,53,0,53,49,43,31,23,
    17,13,41,59,0,59,53,31,19,17,1,7,11,23,37,47,49,0,59,53,47,43,41,37,31,23,19,17,11,1,0,59,53,49,47,31,
    17,13,7,1,11,29,37,0,53,43,31,17,13,7,1,29,41,49,0,53,49,41,29,23,11,7,1,19,31,47,0,47,43,41,29,23,19,
    7,1,11,49,0,59,31,29,23,17,11,7,1,13,41,43,0,59,43,37,17,1,7,11,13,19,41,49,0,59,53,43,41,37,31,29,23,
    13,11,1,47,0,59,53,47,31,19,17,13,1,7,11,29,37,43,49,0,49,43,41,31,17,13,7,11,23,37,53,0,53,49,41,23,
    19,13,11,7,1,17,37,59,0,49,47,43,37,31,29,23,1,7,41,0,59,43,41,37,31,17,13,11,7,47,49,0,59,49,47,37,31,
    29,19,17,7,1,0,53,47,37,19,13,1,11,31,41,0,49,47,37,23,17,13,11,7,19,31,53,0,59,53,47,29,13,11,7,1,23,
    41,0,49,47,41,37,19,11,13,17,23,29,31,43,0,59,29,19,13,1,41,43,47,53,0,59,53,43,41,37,23,17,11,7,1,13,
    29,49 };
#endif

/* microecm's uecm_stage2 (D = 60) with projective A24 */
static uint64_t uedw_stage2(const uedw_pt *P, const uedw_curve *cv, uint32_t stg1_max,
    const uedw_ctx *c)
{
    int b;
    int i, j, k;
    uedw_pt Pa, Pd, Pad, Pb[18], pt4, pt5, pt6;
    uint64_t Pbprod[18];
    const uint8_t *barray = 0;
    int numb = 0;
    const int D = 60;

    Pb[1] = *P;
    Pbprod[1] = MUL(Pb[1].X, Pb[1].Z);
    uedw_mdup(P, &Pb[2], cv, c);                            // [2]Q
    uedw_uadd(Pb[1], Pb[2], Pb[1], &Pb[3], c);              // [3]Q, temporary
    uedw_mdup(&Pb[3], &pt6, cv, c);                         // [6]Q
    uedw_uadd(Pb[3], Pb[2], Pb[1], &pt5, c);                // [5]Q
    Pb[3] = pt5;
    uedw_uadd(pt6, pt5, Pb[1], &Pb[4], c);                  // [11]Q

    i = 3; k = 4; j = 5;
    while ((j + 12) < D) {
        uedw_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);   // [j+12]Q
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw_uadd(pt6, Pb[1], pt5, &Pb[3], c);                  // [7]Q
    i = 1; k = 3; j = 1;
    while ((j + 12) < D) {
        uedw_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw_uadd(Pb[9], Pb[10], Pb[2], &Pd, c);                // [60]Q
    uedw_mdup(&Pb[2], &pt4, cv, c);                         // [4]Q
    for (i = 3; i < 18; i++)
        Pbprod[i] = MUL(Pb[i].X, Pb[i].Z);

    uedw_uadd(Pb[uedw_map[17]], Pb[uedw_map[13]], pt4, &Pad, c);   // [30]Q
    uedw_uadd(Pd, Pad, Pad, &Pa, c);                        // [90]Q
    uedw_pt pt90 = Pa;
    uedw_pt pt60 = Pd;
    uedw_uadd(Pa, Pad, Pd, &Pd, c);                         // [120]Q
    uedw_uadd(Pd, Pad, Pa, &Pa, c);                         // [150]Q

    uint64_t acc = c->one;

    if (stg1_max < 70) {
        static const int steps27[16] = { 59,53,49,47,43,37,31,29,23,19,17,11,7,1,13,41 };
        static const int steps47[15] = { 43,37,31,29,23,19,17,11,7,1,13,41,47,49,59 };
        const int *steps = (stg1_max == 27) ? steps27 : steps47;
        int numsteps = (stg1_max == 27) ? 16 : 15;
        uint64_t pt90prod = MUL(pt90.X, pt90.Z);
        for (i = 0; i < numsteps; i++) {
            b = steps[i];
            uint64_t tt1 = SUB(pt90.X, Pb[uedw_map[b]].X);
            uint64_t tt2 = ADD(pt90.Z, Pb[uedw_map[b]].Z);
            uint64_t tt3 = MUL(tt1, tt2);
            tt1 = ADD(tt3, Pbprod[uedw_map[b]]);
            tt2 = SUB(tt1, pt90prod);
            uint64_t tmp = MUL(acc, tt2);
            if (tmp == 0)
                break;
            acc = tmp;
        }
    } else if (stg1_max == 70) {
        static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
        uint64_t pdprod = MUL(Pd.X, Pd.Z);
        for (i = 0; i < 15; i++) {
            b = steps[i];
            uint64_t tt1 = SUB(Pd.X, Pb[uedw_map[b]].X);
            uint64_t tt2 = ADD(Pd.Z, Pb[uedw_map[b]].Z);
            uint64_t tt3 = MUL(tt1, tt2);
            tt1 = ADD(tt3, Pbprod[uedw_map[b]]);
            tt2 = SUB(tt1, pdprod);
            uint64_t tmp = MUL(acc, tt2);
            if (tmp == 0)
                break;
            acc = tmp;
        }
    } else if (stg1_max == 165) {
        uedw_mdup(&pt90, &Pa, cv, c);                       // [180]Q
        Pad = pt60;
    } else if (stg1_max == 205) {
        uedw_uadd(Pd, pt90, Pad, &Pa, c);                   // [210]Q
        Pad = pt90;
    }

    uint64_t Paprod = MUL(Pa.X, Pa.Z);

    switch (stg1_max) {
    case 27: barray = uedw_b1_27; numb = uedw_numb1_27; break;
    case 47: barray = uedw_b1_47; numb = uedw_numb1_47; break;
    case 70: barray = uedw_b1_70; numb = uedw_numb1_70; break;
    case 85: barray = uedw_b1_85; numb = uedw_numb1_85; break;
    case 125: barray = uedw_b1_125_x50; numb = uedw_numb1_125_x50; break;
    case 165: barray = uedw_b1_165_x50; numb = uedw_numb1_165_x50; break;
    case 205: barray = uedw_b1_205_x50; numb = uedw_numb1_205_x50; break;
    default: return acc;
    }

    for (i = 0; i < numb; i++) {
        if (barray[i] == 0) {
            // giant step
            uedw_pt point = Pa;
            uedw_uadd(Pa, Pd, Pad, &Pa, c);
            Pad = point;
            Paprod = MUL(Pa.X, Pa.Z);
            i++;
        }
        // XrZd - XdZr = (Xr - Xd)(Zr + Zd) + XdZd - XrZr
        b = barray[i];
        uint64_t tt1 = SUB(Pa.X, Pb[uedw_map[b]].X);
        uint64_t tt2 = ADD(Pa.Z, Pb[uedw_map[b]].Z);
        uint64_t tt3 = MUL(tt1, tt2);
        tt1 = ADD(tt3, Pbprod[uedw_map[b]]);
        tt2 = SUB(tt1, Paprod);
        uint64_t tmp = MUL(acc, tt2);
        if (tmp == 0)
            break;
        acc = tmp;
    }
    return acc;
}

/* parameter set: 1 = the legacy levels (microecm's s, B2 = 25 / 50 B1 and
   dispatch), 2 = the re-tuned levels of edw_tiny_v2.h (plan doc Part I),
   the default */
static int uedw_params = 2;

static const edw2_level_t *uedw_level2(uint32_t B1)
{
    int i;
    for (i = 0; i < EDW2_NLEVELS; i++)
        if (edw2_levels[i].B1 == B1)
            return &edw2_levels[i];
    return 0;
}

/* stage 2 with a generic table (edw_tiny_v2.h): the same baby steps as
   uedw_stage2, then giant windows A = 30, 150, 270, ... from [30]Q (with
   [-90]Q as the difference); entry 0 = next window, else a cross term at
   A +- b.  Same accumulation rule (stop before a product 0 mod n). */
static uint64_t uedw_stage2_gen(const uedw_pt *P, const uedw_curve *cv, const uint8_t *s2,
    int s2len, const uedw_ctx *c)
{
    int i, j, k;
    uedw_pt Pa, Pd, Pad, Pb[18], pt4, pt5, pt6;
    uint64_t Pbprod[18];
    const int D = 60;

    Pb[1] = *P;
    Pbprod[1] = MUL(Pb[1].X, Pb[1].Z);
    uedw_mdup(P, &Pb[2], cv, c);
    uedw_uadd(Pb[1], Pb[2], Pb[1], &Pb[3], c);
    uedw_mdup(&Pb[3], &pt6, cv, c);
    uedw_uadd(Pb[3], Pb[2], Pb[1], &pt5, c);
    Pb[3] = pt5;
    uedw_uadd(pt6, pt5, Pb[1], &Pb[4], c);
    i = 3; k = 4; j = 5;
    while ((j + 12) < D) {
        uedw_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw_uadd(pt6, Pb[1], pt5, &Pb[3], c);
    i = 1; k = 3; j = 1;
    while ((j + 12) < D) {
        uedw_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw_uadd(Pb[9], Pb[10], Pb[2], &Pd, c);                // [60]Q
    uedw_mdup(&Pb[2], &pt4, cv, c);
    for (i = 3; i < 18; i++)
        Pbprod[i] = MUL(Pb[i].X, Pb[i].Z);
    uedw_uadd(Pb[uedw_map[17]], Pb[uedw_map[13]], pt4, &Pad, c);   // [30]Q
    uedw_uadd(Pd, Pad, Pad, &Pa, c);                        // [90]Q
    uedw_uadd(Pa, Pad, Pd, &Pd, c);                         // [120]Q

    uedw_pt cur = Pad, prev = Pa;                           // [30]Q, [90]Q = [-90]Q
    uint64_t curprod = MUL(cur.X, cur.Z);
    uint64_t acc = c->one;
    for (i = 0; i < s2len; i++) {
        int b = s2[i];
        if (b == 0) {
            uedw_pt nx;
            uedw_uadd(cur, Pd, prev, &nx, c);
            prev = cur;
            cur = nx;
            curprod = MUL(cur.X, cur.Z);
            continue;
        }
        uint64_t tt1 = SUB(cur.X, Pb[uedw_map[b]].X);
        uint64_t tt2 = ADD(cur.Z, Pb[uedw_map[b]].Z);
        uint64_t tt3 = MUL(tt1, tt2);
        tt1 = ADD(tt3, Pbprod[uedw_map[b]]);
        tt2 = SUB(tt1, curprod);
        uint64_t tmp = MUL(acc, tt2);
        if (tmp == 0)
            break;
        acc = tmp;
    }
    return acc;
}

/* the B2 = 25 B1 pair tables of the vector stage 2 (also for
   microecm_edw_x4.c); 0 if B1 is not a level */
int uedw_s2table25(uint32_t B1, const uint8_t **barray, int *numb)
{
#if !defined(USE_AVX512F) && !defined(USE_AVX2)
    (void)B1; (void)barray; (void)numb;
    return 0;
#else
    switch (B1) {
    case 27: *barray = uedw_b1_27; *numb = uedw_numb1_27; return 1;
    case 47: *barray = uedw_b1_47; *numb = uedw_numb1_47; return 1;
    case 70: *barray = uedw_b1_70; *numb = uedw_numb1_70; return 1;
    case 85: *barray = uedw_b1_85; *numb = uedw_numb1_85; return 1;
    case 125: *barray = uedw_b1_125; *numb = uedw_numb1_125; return 1;
    case 165: *barray = uedw_b1_165; *numb = uedw_numb1_165; return 1;
    case 205: *barray = uedw_b1_205; *numb = uedw_numb1_205; return 1;
    default: return 0;
    }
#endif
}

/* ------------------------------------------------------------- driver */
static int uedw_check(uint64_t x, uint64_t n, uint64_t *f)
{
    uint64_t g = uedw_bingcd64(n, x);
    if (g > 1 && g != n) {
        *f = g;
        return 1;
    }
    return 0;
}

/* The curve source of one input: first the curves of edw_tiny_curves.h
   (family (iv) curves whose constants fit in 63 bits: 6 multiplications to
   load), then GKL [m]Q for m = m0, m0 + 1, ... (~46 multiplications per
   curve, and a short scalar multiplication on E' for [m0]Q, done only if
   the table runs out: most small inputs need just a few curves). */
typedef struct {
    uedw_wpt R;             // [m]Q of the next GKL curve, once started
    uint32_t m0;
    int started;
    int tab;                // next table curve
    int ntab;               // table curves to use
} uedw_src;

static int uedw_ntab_default = EDW_TINY_NCURVES;

static void uedw_src_init(uedw_src *s, uint32_t m0)
{
    s->m0 = m0;
    s->started = 0;
    s->tab = 0;
    s->ntab = uedw_ntab_default;
}

static void uedw_table_curve(int i, uedw_curve *cv, const uedw_ctx *c)
{
    uint64_t v[6];
    int k;
    for (k = 0; k < 6; k++) {
        int64_t s = edw_tiny_curves[i][k];
        uint64_t a = (s < 0) ? (uint64_t)0 - (uint64_t)s : (uint64_t)s;
        uint64_t m = MUL(a, c->Rsqr);           // a < 2^63, so a R^2 < n R
        v[k] = (s < 0) ? SUB(0, m) : m;
    }
    cv->P.X = v[0];
    cv->P.Y = v[1];
    cv->P.Z = v[2];
    cv->P.T = v[3];
    cv->dn2 = ADD(v[4], v[4]);
    cv->dd = v[5];
    cv->dd2 = ADD(v[5], v[5]);
    cv->a24n = v[5];
    cv->a24d = ADD(v[5], v[4]);
}

static void uedw_next_curve(uedw_src *s, uedw_curve *cv, const uedw_ctx *c)
{
    if (s->tab < s->ntab) {
        uedw_table_curve(s->tab++, cv, c);
        return;
    }
    if (!s->started) {
        uedw_w_start(&s->R, s->m0, c);
        s->started = 1;
    }
    uedw_gkl_curve(&s->R, cv, c);
    uedw_w_addq(&s->R, c);                      // the next curve's point
}

#ifdef UEDW_STATS
#define UEDW_HIST_MAX 64
static uint64_t uedw_hist[UEDW_HIST_MAX + 1];
#endif

/* run curves from the source; returns the number of curves run; *f =
   factor or 1; *stage = stage of the find */
static int uedw_curves(uint64_t *f, uint32_t B1, uint32_t B2, uint32_t curves,
    const uedw_ctx *c, uedw_src *src, int *stage)
{
    uedw_chain ch;
    uint32_t curve;
    const edw2_level_t *lv = 0;
    *f = 1;
    *stage = 0;
    if (uedw_params == 2) {
        if (!(lv = uedw_level2(B1)))
            return 0;
        ch.ch = lv->ch;
        ch.len = lv->len;
        ch.m = lv->m;
        ch.top = lv->top;
        ch.tz = lv->tz;
    } else if (!uedw_chain_for(B1, &ch)) {
        return 0;
    }
    for (curve = 0; curve < curves; curve++) {
        uedw_curve cv;
        uedw_pt P;
        uedw_next_curve(src, &cv, c);
        uedw_stage1(&cv, &ch, &P, c);
        if (uedw_check(P.Z, c->n, f)) {
            *stage = 1;
            break;
        }
        if (B2 > B1) {
            uint64_t acc = lv ? uedw_stage2_gen(&P, &cv, lv->s2, lv->s2len, c)
                : uedw_stage2(&P, &cv, B1, c);
            if (uedw_check(acc, c->n, f)) {
                *stage = 2;
                break;
            }
        }
    }
    return (int)curve;
}

static uint64_t uedw_run(uint64_t *f, uint32_t B1, uint32_t curves, const uedw_ctx *c,
    uedw_src *src, int *ncurves)
{
    int stage;
    *ncurves += uedw_curves(f, B1, 25 * B1, curves, c, src, &stage) + (*f > 1);
    return *f;
}

static uint32_t uedw_lcg_m0(uint64_t *pran)
{
    // start of the GKL sequence: m0 in [3, 258] from the caller's LCG
    // (m = 1 and 2 are in the table)
    *pran = 6364136223846793005ULL * (*pran) + 1442695040888963407ULL;
    return 3 + (uint32_t)((*pran >> 40) & 0xff);
}

/* the P-1 pre-pass of uedw_dispatch (parameter set 2 only); compile with
   UEDW_NO_PM1 to build without micropm1.c */
static int uedw_pm1 = 1;
#ifndef UEDW_NO_PM1
uint64_t getfactor_upm1(uint64_t q64, uint32_t b1);     // micropm1.c
#endif

static uint64_t uedw_dispatch(uint64_t n, int targetBits, int arbitrary, uint64_t *pran)
{
    uedw_ctx ctx;
    uedw_src src;
    uint64_t f = 1;
    int nc = 0;
    uedw_ctx_init(&ctx, n);
    uedw_src_init(&src, uedw_lcg_m0(pran));

    if (arbitrary) {
        // fast attempts at possible small factors
        if (uedw_run(&f, 47, 1, &ctx, &src, &nc) > 1) goto done;
        if (uedw_run(&f, 70, 1, &ctx, &src, &nc) > 1) goto done;
        if (targetBits > 58)
            if (uedw_run(&f, 125, 1, &ctx, &src, &nc) > 1) goto done;
    }
    if (uedw_params == 2) {
#ifndef UEDW_NO_PM1
        // P-1 pre-pass (micropm1, B1 = 333, B2 = 25 B1): it finds 16-34% of
        // balanced n at 56-64 bits for ~3.7 us, a 2-3% saving per input
        // over ECM alone; break-even at 48-52 bits (arith/tinyrho/pm1_eval.c).
        // A larger B2 loses (stage 2 costs ~1 mulmod per prime).
        if (uedw_pm1 && targetBits > 54) {
            f = getfactor_upm1(n, 333);
            if (f > 1) goto done;
            f = 1;
        }
#endif
        // re-tuned levels (Part I): balanced n, factor ~ targetBits / 2.
        // Curve caps: the re-tuned levels trade per-curve success for
        // cheaper curves (64 bits: 0.19 vs legacy 0.24 per curve), so the
        // legacy caps of 32 / 42 left 1.5e-4 of 64-bit inputs unfactored.
        // The caps below give ~1e-7 at every size (uedw_tail.c); the extra
        // curves only run on inputs that would otherwise fail.
        if (targetBits <= 36) uedw_run(&f, 20, 32, &ctx, &src, &nc);
        else if (targetBits <= 40) uedw_run(&f, 27, 32, &ctx, &src, &nc);
        else if (targetBits <= 44) uedw_run(&f, 47, 32, &ctx, &src, &nc);
        else if (targetBits <= 48) uedw_run(&f, 70, 32, &ctx, &src, &nc);
        else if (targetBits <= 52) uedw_run(&f, 85, 32, &ctx, &src, &nc);
        else if (targetBits <= 56) uedw_run(&f, 125, 48, &ctx, &src, &nc);
        else if (targetBits <= 60) uedw_run(&f, 165, 80, &ctx, &src, &nc);
        else uedw_run(&f, 205, 80, &ctx, &src, &nc);
        goto done;
    }
    if (targetBits <= 40) uedw_run(&f, 27, 32, &ctx, &src, &nc);
    else if (targetBits <= 44) uedw_run(&f, 47, 32, &ctx, &src, &nc);
    else if (targetBits <= 48) uedw_run(&f, 70, 32, &ctx, &src, &nc);
    else if (targetBits <= 52) uedw_run(&f, 85, 32, &ctx, &src, &nc);
    else if (targetBits <= 58) uedw_run(&f, 125, 32, &ctx, &src, &nc);
    else if (targetBits <= 62) uedw_run(&f, 165, 42, &ctx, &src, &nc);
    else uedw_run(&f, 205, 42, &ctx, &src, &nc);
done:
#ifdef UEDW_STATS
    // curves run per input (the finding curve included), not thread-safe
    uedw_hist[(nc > UEDW_HIST_MAX) ? UEDW_HIST_MAX : nc]++;
#endif
    (void)nc;
    return f;
}

static int uedw_bits(uint64_t n)
{
    int b = 0;
    while (n) {
        b++;
        n >>= 1;
    }
    return b;
}

/* research hook: up to `curves` curves of one level (no pre-passes), the
   curve sequence from *pran as in getfactor_uecm_edw; *ncurves = curves run
   (the finding one included) */
uint64_t uecm_edw_level(uint64_t n, uint32_t B1, uint32_t curves, uint64_t *pran, int *ncurves)
{
    uedw_ctx ctx;
    uedw_src src;
    uint64_t f = 1;
    int nc = 0;
    uedw_ctx_init(&ctx, n);
    uedw_src_init(&src, uedw_lcg_m0(pran));
    uedw_run(&f, B1, curves, &ctx, &src, &nc);
    *ncurves = nc;
    return f;
}

uint64_t getfactor_uecm_edw(uint64_t q64, int is_arbitrary, uint64_t *pran)
{
    if (q64 % 2 == 0)
        return 2;
    return uedw_dispatch(q64, uedw_bits(q64), is_arbitrary, pran);
}

/* test hooks: a curve id >= 1 is GKL m = id; UECM_EDW_TABLE | i is curve i
   of edw_tiny_curves.h */
static void uedw_curve_by_id(uint32_t id, uedw_curve *cv, const uedw_ctx *c)
{
    if (id & UECM_EDW_TABLE) {
        uedw_table_curve((int)(id & ~UECM_EDW_TABLE) % EDW_TINY_NCURVES, cv, c);
    } else {
        uedw_wpt R;
        uedw_w_start(&R, id, c);
        uedw_gkl_curve(&R, cv, c);
    }
}

uint64_t uecm_edw_one_curve(uint64_t n, uint32_t id, uint32_t B1, uint32_t B2, int *stage)
{
    uedw_ctx ctx;
    uedw_src src;
    uint64_t f = 1;
    uedw_ctx_init(&ctx, n);
    if (id & UECM_EDW_TABLE) {
        uedw_src_init(&src, 3);
        src.tab = (int)(id & ~UECM_EDW_TABLE) % EDW_TINY_NCURVES;
        src.ntab = src.tab + 1;
    } else {
        uedw_src_init(&src, id);
        src.ntab = 0;
    }
    uedw_curves(&f, B1, B2, 1, &ctx, &src, stage);
    return f;
}

void uecm_edw_curve_params(uint64_t n, uint32_t id, uint64_t out[6])
{
    uedw_ctx ctx;
    uedw_curve cv;
    const uedw_ctx *c = &ctx;
    uedw_ctx_init(&ctx, n);
    uedw_curve_by_id(id, &cv, &ctx);
    // out of Montgomery form: x R / R
    out[0] = MUL(cv.P.X, 1);
    out[1] = MUL(cv.P.Y, 1);
    out[2] = MUL(cv.P.Z, 1);
    out[3] = MUL(cv.P.T, 1);
    out[4] = MUL(SUB(cv.a24d, cv.dd), 1);       // dn
    out[5] = MUL(cv.dd, 1);
}

void uecm_edw_stage1_point(uint64_t n, uint32_t id, uint32_t B1, uint64_t out[2])
{
    uedw_ctx ctx;
    uedw_curve cv;
    uedw_chain ch;
    uedw_pt P;
    const uedw_ctx *c = &ctx;
    uedw_ctx_init(&ctx, n);
    uedw_curve_by_id(id, &cv, &ctx);
    out[0] = out[1] = 0;
    if (uedw_params == 2) {
        const edw2_level_t *lv = uedw_level2(B1);
        if (!lv)
            return;
        ch.ch = lv->ch;
        ch.len = lv->len;
        ch.m = lv->m;
        ch.top = lv->top;
        ch.tz = lv->tz;
    } else if (!uedw_chain_for(B1, &ch)) {
        return;
    }
    uedw_stage1(&cv, &ch, &P, &ctx);
    out[0] = MUL(P.X, 1);
    out[1] = MUL(P.Z, 1);
}

void uecm_edw_set_params(int v)
{
    uedw_params = (v == 2) ? 2 : 1;
}

int uecm_edw_get_params(void)
{
    return uedw_params;
}

void uecm_edw_set_pm1(int on)
{
    uedw_pm1 = on ? 1 : 0;
}

const edw2_level_t *uedw_level2_pub(uint32_t B1)
{
    return uedw_level2(B1);
}

void uecm_edw_set_table(int ntab)
{
    uedw_ntab_default = (ntab < 0) ? 0 : (ntab > EDW_TINY_NCURVES) ? EDW_TINY_NCURVES : ntab;
}

int uecm_edw_hist(uint64_t *hist, int max)
{
#ifdef UEDW_STATS
    int i;
    for (i = 0; i <= max && i <= UEDW_HIST_MAX; i++) {
        hist[i] = uedw_hist[i];
        uedw_hist[i] = 0;
    }
    return UEDW_HIST_MAX;
#else
    (void)hist;
    (void)max;
    return -1;
#endif
}

/* ======================================================================
   x8: AVX-512, 8 inputs n < 2^52 at a time (one per lane), R = 2^52.
   The same curves, chains and stage 2 as the scalar code, vectorized:
   - per lane, the curve sequence is the table curves, then GKL [m]Q for
     m = 3, 4, ... (a fixed start: [3]Q is a small constant, so no lane
     needs a scalar multiplication on E');
   - table constants are loaded per lane as two 52-bit digits,
     a0 R^2 / R + a1 R^3 / R;
   - stage 2 is microecm's uecm_stage2_x8 (B2 = 25 B1 tables at every
     level; a lane's accumulator skips products that would be 0 mod n).
   ====================================================================== */
#ifdef USE_AVX512F
#include <immintrin.h>

typedef __m512i uedw_v8;

/* Montgomery multiplication (R = 2^52, positive inverse, inputs in [0, n)),
   modular add / sub and the 52 x 52 product halves, 8 lanes: shared, in
   arith/limb/limb1.h.  The non-IFMA FMA bias constants are local to those
   kernels, so this file no longer keeps (or initializes) its own. */
#define uedw8_mul52lo   mp_mul52lo_x8
#define uedw8_mul52hi   mp_mul52hi_x8
#define uedw8_mul52lohi mp_mul52lohi_x8
#define uedw8_mulredc   mulredc52_x8
#define uedw8_addmod    addmod52_x8
#define uedw8_submod    submod52_x8

typedef struct {
    uedw_v8 n, rho, one, Rsqr, R3;
    uedw_v8 c6, c10, c16, c32, c48;
} uedw_ctx8;

#define M8(a, b) uedw8_mulredc((a), (b), c->n, c->rho)
#define S8(a) uedw8_mulredc((a), (a), c->n, c->rho)
#define A8(a, b) uedw8_addmod((a), (b), c->n)
#define D8(a, b) uedw8_submod((a), (b), c->n)
#define Z8 _mm512_setzero_si512()

/* the derived constants, from n, rho, one and Rsqr */
static void uedw8_ctx_finish(uedw_ctx8 *c)
{
    __m512i c2 = A8(c->one, c->one);
    __m512i c4 = A8(c2, c2);
    __m512i c8 = A8(c4, c4);
    c->c16 = A8(c8, c8);
    c->c32 = A8(c->c16, c->c16);
    c->c48 = A8(c->c32, c->c16);
    c->c6 = A8(c4, c2);
    c->c10 = A8(c8, c2);
    c->R3 = M8(c->Rsqr, c->Rsqr);               // R^3 mod n
}

/* a signed constant < 2^63 per lane, to Montgomery form */
static __m512i uedw8_const(const int64_t *v, const uedw_ctx8 *c)
{
    ALIGNED_MEM uint64_t lo[8], hi[8];
    __mmask8 neg = 0;
    int i;
    for (i = 0; i < 8; i++) {
        uint64_t a = (v[i] < 0) ? (uint64_t)0 - (uint64_t)v[i] : (uint64_t)v[i];
        lo[i] = a & 0x000fffffffffffffULL;
        hi[i] = a >> 52;
        if (v[i] < 0) neg |= (__mmask8)(1 << i);
    }
    __m512i r = A8(M8(_mm512_load_si512(lo), c->Rsqr), M8(_mm512_load_si512(hi), c->R3));
    return _mm512_mask_mov_epi64(r, neg, D8(Z8, r));
}

typedef struct { uedw_v8 X, Y, Z; } uedw_wpt8;
typedef struct { uedw_v8 X, Y, Z, T; } uedw_ept8;
typedef struct { uedw_v8 YmX, YpX, T2d, Z2; } uedw_cached8;
typedef struct {
    uedw_ept8 P;
    uedw_v8 dn2, dd, dd2, a24n, a24d;
} uedw_curve8;
typedef struct { uedw_v8 X, Z; } uedw_pt8;

static void uedw8_w_addq(uedw_wpt8 *R, const uedw_ctx8 *c)
{
    __m512i X = R->X, Y = R->Y, Z = R->Z;
    __m512i z2 = A8(Z, Z);
    __m512i z4 = A8(z2, z2);
    __m512i u = D8(A8(z4, z4), Y);
    __m512i v = D8(A8(z4, Z), X);
    __m512i u2 = S8(u);
    __m512i v2 = S8(v);
    __m512i v3 = M8(v2, v);
    __m512i w1 = M8(v2, X);
    __m512i w2 = M8(v2, Z);
    w2 = A8(w2, w2);
    w2 = A8(w2, w2);
    __m512i A = D8(D8(M8(u2, Z), w1), w2);
    R->X = M8(v, A);
    R->Y = D8(M8(u, D8(w1, A)), M8(Y, v3));
    R->Z = M8(v3, Z);
}

static void uedw8_gkl_curve(const uedw_wpt8 *R, uedw_curve8 *cv, const uedw_ctx8 *c)
{
    __m512i t = A8(R->X, R->Z);
    t = A8(t, t);
    __m512i a = A8(t, t);
    __m512i z4 = A8(R->Z, R->Z);
    z4 = A8(z4, z4);
    __m512i b = D8(R->Y, z4);
    __m512i a2 = S8(a), b2 = S8(b), ab = M8(a, b);
    __m512i b2x = A8(b, b);
    __m512i xn = M8(M8(D8(a, b2x), A8(a, b2x)), A8(a, b));
    xn = M8(xn, b);
    xn = A8(xn, xn);
    __m512i aa4 = M8(a, A8(a, A8(b2x, b2x)));
    __m512i b2_4 = A8(b2, b2);
    b2_4 = A8(b2_4, b2_4);
    __m512i xd = M8(aa4, A8(A8(a2, A8(ab, ab)), b2_4));
    __m512i t1 = S8(a2), t2 = M8(a2, ab), t3 = M8(a2, b2), t4 = M8(ab, b2), t5 = S8(b2);
    __m512i q = A8(t2, t3);
    q = A8(t1, A8(q, q));
    __m512i r = A8(t4, A8(t5, t5));
    r = A8(r, r);
    r = A8(r, r);
    r = A8(r, r);
    __m512i yd = M8(aa4, A8(q, r));
    __m512i b3 = M8(b2, b), b4 = S8(b2), b5 = M8(b4, b), b6 = S8(b3);
    __m512i h = A8(a, M8(b, c->c6));
    h = A8(M8(h, a), M8(b2, c->c10));
    h = D8(M8(h, a), M8(b3, c->c16));
    h = D8(M8(h, a), M8(b4, c->c48));
    h = D8(M8(h, a), M8(b5, c->c32));
    __m512i yn = D8(M8(h, a), M8(b6, c->c32));
    __m512i en2 = S8(aa4), ed = D8(a2, b2_4), ed2 = S8(ed);
    __m512i dn = D8(Z8, S8(en2));
    __m512i dd = S8(ed2);
    cv->P.X = M8(xn, yd);
    cv->P.Y = M8(yn, xd);
    cv->P.Z = M8(xd, yd);
    cv->P.T = M8(xn, yn);
    cv->dn2 = A8(dn, dn);
    cv->dd = dd;
    cv->dd2 = A8(dd, dd);
    cv->a24n = dd;
    cv->a24d = A8(dd, dn);
}

/* table curve k[lane] in every lane */
static void uedw8_table_curve(const int *k, uedw_curve8 *cv, const uedw_ctx8 *c)
{
    int64_t v[8];
    __m512i w[6];
    int i, j;
    for (j = 0; j < 6; j++) {
        for (i = 0; i < 8; i++)
            v[i] = edw_tiny_curves[k[i]][j];
        w[j] = uedw8_const(v, c);
    }
    cv->P.X = w[0];
    cv->P.Y = w[1];
    cv->P.Z = w[2];
    cv->P.T = w[3];
    cv->dn2 = A8(w[4], w[4]);
    cv->dd = w[5];
    cv->dd2 = A8(w[5], w[5]);
    cv->a24n = w[5];
    cv->a24d = A8(w[5], w[4]);
}

static UEDW_INLINE void uedw8_dbl(uedw_ept8 *R, const uedw_ept8 *P, int want_t, const uedw_ctx8 *c)
{
    __m512i A = S8(P->X);
    __m512i B = S8(P->Y);
    __m512i C = S8(P->Z);
    __m512i S = S8(A8(P->X, P->Y));
    C = A8(C, C);
    __m512i H = A8(A, B);
    __m512i E = D8(H, S);
    __m512i G = D8(A, B);
    __m512i F = A8(C, G);
    R->X = M8(E, F);
    R->Y = M8(G, H);
    R->Z = M8(F, G);
    if (want_t)
        R->T = M8(E, H);
}

static UEDW_INLINE void uedw8_add(uedw_ept8 *R, const uedw_ept8 *P, const uedw_cached8 *Q,
    int neg, int want_t, const uedw_ctx8 *c)
{
    __m512i A = M8(D8(P->Y, P->X), neg ? Q->YpX : Q->YmX);
    __m512i B = M8(A8(P->Y, P->X), neg ? Q->YmX : Q->YpX);
    __m512i C = M8(P->T, Q->T2d);
    __m512i D = M8(P->Z, Q->Z2);
    __m512i E = D8(B, A);
    __m512i H = A8(B, A);
    __m512i F = neg ? A8(D, C) : D8(D, C);
    __m512i G = neg ? D8(D, C) : A8(D, C);
    R->X = M8(E, F);
    R->Y = M8(G, H);
    R->Z = M8(F, G);
    if (want_t)
        R->T = M8(E, H);
}

static UEDW_INLINE void uedw8_cache(uedw_cached8 *q, const uedw_ept8 *P, const uedw_curve8 *cv,
    const uedw_ctx8 *c)
{
    q->YmX = M8(D8(P->Y, P->X), cv->dd);
    q->YpX = M8(A8(P->Y, P->X), cv->dd);
    q->T2d = M8(P->T, cv->dn2);
    q->Z2 = M8(P->Z, cv->dd2);
}

static UEDW_INLINE void uedw8_mdup(const uedw_pt8 *P, uedw_pt8 *R, const uedw_curve8 *cv,
    const uedw_ctx8 *c)
{
    __m512i U = S8(D8(P->X, P->Z));
    __m512i V = S8(A8(P->X, P->Z));
    __m512i w = D8(V, U);
    __m512i Ub = M8(U, cv->a24d);
    R->X = M8(Ub, V);
    R->Z = M8(w, A8(Ub, M8(w, cv->a24n)));
}

static UEDW_INLINE void uedw8_uadd(const uedw_pt8 P1, const uedw_pt8 P2, const uedw_pt8 Pin,
    uedw_pt8 *Pout, const uedw_ctx8 *c)
{
    __m512i diff1 = D8(P1.X, P1.Z);
    __m512i sum1 = A8(P1.X, P1.Z);
    __m512i diff2 = D8(P2.X, P2.Z);
    __m512i sum2 = A8(P2.X, P2.Z);
    __m512i tt1 = M8(diff1, sum2);
    __m512i tt2 = M8(sum1, diff2);
    __m512i tt3 = A8(tt1, tt2);
    __m512i tt4 = D8(tt1, tt2);
    tt1 = S8(tt3);
    tt2 = S8(tt4);
    __m512i tmpx = M8(tt1, Pin.Z);
    __m512i tmpz = M8(tt2, Pin.X);
    Pout->X = tmpx;
    Pout->Z = tmpz;
}

static void uedw8_stage1(const uedw_curve8 *cv, const uedw_chain *ch, uedw_pt8 *Pm,
    const uedw_ctx8 *c)
{
    uedw_cached8 tab[16];
    uedw_ept8 ext[16];
    int k, i;
    ext[0] = cv->P;
    uedw8_cache(&tab[0], &ext[0], cv, c);
    if (ch->m > 1) {
        uedw_ept8 P2;
        uedw_cached8 c2;
        uedw8_dbl(&P2, &cv->P, 1, c);
        uedw8_cache(&c2, &P2, cv, c);
        for (k = 1; k <= (ch->m - 1) / 2; k++) {
            uedw8_add(&ext[k], &ext[k - 1], &c2, 0, 1, c);
            uedw8_cache(&tab[k], &ext[k], cv, c);
        }
    }
    uedw_ept8 R = ext[(ch->top - 1) / 2];
    for (i = 0; i < ch->len; i++) {
        int dg = ch->ch[i];
        if (dg == 0) {
            int t = (i + 1 < ch->len) && (ch->ch[i + 1] != 0);
            uedw8_dbl(&R, &R, t, c);
        } else if (dg > 0) {
            uedw8_add(&R, &R, &tab[(dg - 1) >> 1], 0, 0, c);
        } else {
            uedw8_add(&R, &R, &tab[(-dg - 1) >> 1], 1, 0, c);
        }
    }
    Pm->X = A8(R.Z, R.Y);
    Pm->Z = D8(R.Z, R.Y);
    for (i = 0; i < ch->tz; i++)
        uedw8_mdup(Pm, Pm, cv, c);
}

/* acc *= t, skipping lanes where the product would be 0 mod n */
#define UEDW8_ACC(t) do { __m512i tmp_ = M8(acc, (t)); \
    acc = _mm512_mask_mov_epi64(acc, _mm512_cmpneq_epu64_mask(tmp_, Z8), tmp_); } while (0)

static __m512i uedw8_cross(const uedw_pt8 *R, __m512i Rprod, const uedw_pt8 *Pb,
    __m512i Pbprod, const uedw_ctx8 *c)
{
    __m512i tt1 = D8(R->X, Pb->X);
    __m512i tt2 = A8(R->Z, Pb->Z);
    __m512i tt3 = M8(tt1, tt2);
    tt1 = A8(tt3, Pbprod);
    return D8(tt1, Rprod);
}

/* stage 2 with a generic table (edw_tiny_v2.h; as uedw_stage2_gen) */
static __m512i uedw8_stage2_gen(const uedw_pt8 *P, const uedw_curve8 *cv, const uint8_t *s2,
    int s2len, const uedw_ctx8 *c)
{
    int i, j, k;
    uedw_pt8 Pa, Pd, Pad, Pb[18], pt4, pt5, pt6;
    __m512i Pbprod[18];
    const int D = 60;

    Pb[1] = *P;
    Pbprod[1] = M8(Pb[1].X, Pb[1].Z);
    uedw8_mdup(P, &Pb[2], cv, c);
    uedw8_uadd(Pb[1], Pb[2], Pb[1], &Pb[3], c);
    uedw8_mdup(&Pb[3], &pt6, cv, c);
    uedw8_uadd(Pb[3], Pb[2], Pb[1], &pt5, c);
    Pb[3] = pt5;
    uedw8_uadd(pt6, pt5, Pb[1], &Pb[4], c);
    i = 3; k = 4; j = 5;
    while ((j + 12) < D) {
        uedw8_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw8_uadd(pt6, Pb[1], pt5, &Pb[3], c);
    i = 1; k = 3; j = 1;
    while ((j + 12) < D) {
        uedw8_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw8_uadd(Pb[9], Pb[10], Pb[2], &Pd, c);
    uedw8_mdup(&Pb[2], &pt4, cv, c);
    for (i = 3; i < 18; i++)
        Pbprod[i] = M8(Pb[i].X, Pb[i].Z);
    uedw8_uadd(Pb[uedw_map[17]], Pb[uedw_map[13]], pt4, &Pad, c);  // [30]Q
    uedw8_uadd(Pd, Pad, Pad, &Pa, c);                       // [90]Q
    uedw8_uadd(Pa, Pad, Pd, &Pd, c);                        // [120]Q

    uedw_pt8 cur = Pad, prev = Pa;
    __m512i curprod = M8(cur.X, cur.Z);
    __m512i acc = c->one;
    for (i = 0; i < s2len; i++) {
        int b = s2[i];
        if (b == 0) {
            uedw_pt8 nx;
            uedw8_uadd(cur, Pd, prev, &nx, c);
            prev = cur;
            cur = nx;
            curprod = M8(cur.X, cur.Z);
            continue;
        }
        UEDW8_ACC(uedw8_cross(&cur, curprod, &Pb[uedw_map[b]], Pbprod[uedw_map[b]], c));
    }
    return acc;
}

/* microecm's uecm_stage2_x8 (D = 60, B2 = 25 B1) with projective A24 */
static __m512i uedw8_stage2(const uedw_pt8 *P, const uedw_curve8 *cv, uint32_t stg1_max,
    const uedw_ctx8 *c)
{
    int b, i, j, k;
    uedw_pt8 Pa, Pd, Pad, Pb[18], pt4, pt5, pt6;
    __m512i Pbprod[18];
    const uint8_t *barray = 0;
    int numb = 0;
    const int D = 60;

    Pb[1] = *P;
    Pbprod[1] = M8(Pb[1].X, Pb[1].Z);
    uedw8_mdup(P, &Pb[2], cv, c);
    uedw8_uadd(Pb[1], Pb[2], Pb[1], &Pb[3], c);
    uedw8_mdup(&Pb[3], &pt6, cv, c);
    uedw8_uadd(Pb[3], Pb[2], Pb[1], &pt5, c);
    Pb[3] = pt5;
    uedw8_uadd(pt6, pt5, Pb[1], &Pb[4], c);
    i = 3; k = 4; j = 5;
    while ((j + 12) < D) {
        uedw8_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw8_uadd(pt6, Pb[1], pt5, &Pb[3], c);
    i = 1; k = 3; j = 1;
    while ((j + 12) < D) {
        uedw8_uadd(pt6, Pb[k], Pb[i], &Pb[uedw_map[j + 12]], c);
        i = k;
        k = uedw_map[j + 12];
        j += 6;
    }
    uedw8_uadd(Pb[9], Pb[10], Pb[2], &Pd, c);
    uedw8_mdup(&Pb[2], &pt4, cv, c);
    for (i = 3; i < 18; i++)
        Pbprod[i] = M8(Pb[i].X, Pb[i].Z);
    uedw8_uadd(Pb[uedw_map[17]], Pb[uedw_map[13]], pt4, &Pad, c);
    uedw8_uadd(Pd, Pad, Pad, &Pa, c);
    uedw_pt8 pt90 = Pa;
    uedw_pt8 pt60 = Pd;
    uedw8_uadd(Pa, Pad, Pd, &Pd, c);
    uedw8_uadd(Pd, Pad, Pa, &Pa, c);

    __m512i acc = c->one;
    if (stg1_max < 70) {
        static const int steps27[16] = { 59,53,49,47,43,37,31,29,23,19,17,11,7,1,13,41 };
        static const int steps47[15] = { 43,37,31,29,23,19,17,11,7,1,13,41,47,49,59 };
        const int *steps = (stg1_max == 27) ? steps27 : steps47;
        int numsteps = (stg1_max == 27) ? 16 : 15;
        __m512i pt90prod = M8(pt90.X, pt90.Z);
        for (i = 0; i < numsteps; i++) {
            b = steps[i];
            UEDW8_ACC(uedw8_cross(&pt90, pt90prod, &Pb[uedw_map[b]], Pbprod[uedw_map[b]], c));
        }
    } else if (stg1_max == 70) {
        static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
        __m512i pdprod = M8(Pd.X, Pd.Z);
        for (i = 0; i < 15; i++) {
            b = steps[i];
            UEDW8_ACC(uedw8_cross(&Pd, pdprod, &Pb[uedw_map[b]], Pbprod[uedw_map[b]], c));
        }
    } else if (stg1_max == 165) {
        uedw8_mdup(&pt90, &Pa, cv, c);
        Pad = pt60;
    } else if (stg1_max == 205) {
        uedw8_uadd(Pd, pt90, Pad, &Pa, c);
        Pad = pt90;
    }
    __m512i Paprod = M8(Pa.X, Pa.Z);
    switch (stg1_max) {
    case 27: barray = uedw_b1_27; numb = uedw_numb1_27; break;
    case 47: barray = uedw_b1_47; numb = uedw_numb1_47; break;
    case 70: barray = uedw_b1_70; numb = uedw_numb1_70; break;
    case 85: barray = uedw_b1_85; numb = uedw_numb1_85; break;
    case 125: barray = uedw_b1_125; numb = uedw_numb1_125; break;
    case 165: barray = uedw_b1_165; numb = uedw_numb1_165; break;
    case 205: barray = uedw_b1_205; numb = uedw_numb1_205; break;
    default: return acc;
    }
    for (i = 0; i < numb; i++) {
        if (barray[i] == 0) {
            uedw_pt8 point = Pa;
            uedw8_uadd(Pa, Pd, Pad, &Pa, c);
            Pad = point;
            Paprod = M8(Pa.X, Pa.Z);
            i++;
        }
        b = barray[i];
        UEDW8_ACC(uedw8_cross(&Pa, Paprod, &Pb[uedw_map[b]], Pbprod[uedw_map[b]], c));
    }
    return acc;
}

static uint64_t uedw8_mulmod_plain(uint64_t a, uint64_t b, uint64_t n)
{
#if defined(__SIZEOF_INT128__)
    return (uint64_t)(((unsigned __int128)a * b) % n);
#else
    // n < 2^52: double-and-add, sums stay below 2^53
    uint64_t r = 0, x = a % n;
    uint64_t y = b % n;
    while (y) {
        if (y & 1) r = (r + x) % n;
        x = (x + x) % n;
        y >>= 1;
    }
    return r;
#endif
}

static uint64_t uedw8_inverse52(uint64_t a)
{
    return uedw_multiplicative_inverse(a) & 0x000fffffffffffffULL;
}

/* ECM on n[0..num-1] (all < 2^52), 8 at a time; f[i] = factor or 1.
   *ncurves (if not 0) gets the total number of curves run. */
static void uedw8_list(const uint64_t *n64, uint64_t *f, uint32_t num, uint32_t B1,
    uint32_t B2, uint32_t curves, uint64_t *ncurves)
{
    ALIGNED_MEM uint64_t nA[8], rhoA[8], oneA[8], RsqA[8], st[8];
    int kcur[8], idx[8], k8[8];
    uint8_t msk = 0;
    uint32_t j = 0, i;
    int lane;
    uedw_chain ch;
    uedw_ctx8 ctx;
    const uedw_ctx8 *c = &ctx;
    uedw_wpt8 R;
    uint64_t total = 0;
    const int ntab = EDW_TINY_NCURVES;

    for (i = 0; i < num; i++)
        f[i] = 1;
    const edw2_level_t *lv = 0;
    if (num == 0)
        return;
    if (uedw_params == 2) {
        if (!(lv = uedw_level2(B1)))
            return;
        ch.ch = lv->ch;
        ch.len = lv->len;
        ch.m = lv->m;
        ch.top = lv->top;
        ch.tz = lv->tz;
    } else if (!uedw_chain_for(B1, &ch)) {
        return;
    }
    // every lane starts with a valid modulus, so idle lanes compute harmlessly
    for (lane = 0; lane < 8; lane++) {
        nA[lane] = n64[0];
        rhoA[lane] = uedw8_inverse52(n64[0]);
        oneA[lane] = (1ULL << 52) % n64[0];
        RsqA[lane] = uedw8_mulmod_plain(oneA[lane], oneA[lane], n64[0]);
        kcur[lane] = 0;
        idx[lane] = 0;
    }
    memset(&R, 0, sizeof(R));

    while (1) {
        uint8_t lmsk = 0;
        for (lane = 0; lane < 8; lane++) {
            if (!(msk & (1 << lane)) && j < num) {
                uint64_t n = n64[j];
                nA[lane] = n;
                rhoA[lane] = uedw8_inverse52(n);
                oneA[lane] = (1ULL << 52) % n;
                RsqA[lane] = uedw8_mulmod_plain(oneA[lane], oneA[lane], n);
                kcur[lane] = 0;
                idx[lane] = (int)j;
                lmsk |= (uint8_t)(1 << lane);
                j++;
            }
        }
        msk |= lmsk;
        if (msk == 0)
            break;
        if (lmsk) {
            ctx.n = _mm512_load_si512(nA);
            ctx.rho = _mm512_load_si512(rhoA);
            ctx.one = _mm512_load_si512(oneA);
            ctx.Rsqr = _mm512_load_si512(RsqA);
            uedw8_ctx_finish(&ctx);
        }

        // this round's curve in every lane: table curves, then GKL m = 3, 4, ...
        uint8_t tabm = 0, newm = 0;
        for (lane = 0; lane < 8; lane++) {
            k8[lane] = (kcur[lane] < ntab) ? kcur[lane] : 0;
            if (kcur[lane] < ntab) tabm |= (uint8_t)(1 << lane);
            if (kcur[lane] == ntab) newm |= (uint8_t)(1 << lane);
        }
        uedw_curve8 cv, cvg;
        if (tabm)
            uedw8_table_curve(k8, &cv, c);
        if ((uint8_t)~tabm) {
            if (newm) {
                // [3]Q = (1205/49, -40664/343) = (413315 : -1992536 : 16807)
                int64_t q3[3][8];
                for (lane = 0; lane < 8; lane++) {
                    q3[0][lane] = 413315;
                    q3[1][lane] = -1992536;
                    q3[2][lane] = 16807;
                }
                R.X = _mm512_mask_mov_epi64(R.X, newm, uedw8_const(q3[0], c));
                R.Y = _mm512_mask_mov_epi64(R.Y, newm, uedw8_const(q3[1], c));
                R.Z = _mm512_mask_mov_epi64(R.Z, newm, uedw8_const(q3[2], c));
            }
            uedw8_gkl_curve(&R, &cvg, c);
            uedw_wpt8 R2 = R;
            uedw8_w_addq(&R2, c);
            __mmask8 gm = (__mmask8)(uint8_t)~tabm;
            R.X = _mm512_mask_mov_epi64(R.X, gm, R2.X);
            R.Y = _mm512_mask_mov_epi64(R.Y, gm, R2.Y);
            R.Z = _mm512_mask_mov_epi64(R.Z, gm, R2.Z);
            if (!tabm) {
                cv = cvg;
            } else {
                __mmask8 m = gm;
                cv.P.X = _mm512_mask_mov_epi64(cv.P.X, m, cvg.P.X);
                cv.P.Y = _mm512_mask_mov_epi64(cv.P.Y, m, cvg.P.Y);
                cv.P.Z = _mm512_mask_mov_epi64(cv.P.Z, m, cvg.P.Z);
                cv.P.T = _mm512_mask_mov_epi64(cv.P.T, m, cvg.P.T);
                cv.dn2 = _mm512_mask_mov_epi64(cv.dn2, m, cvg.dn2);
                cv.dd = _mm512_mask_mov_epi64(cv.dd, m, cvg.dd);
                cv.dd2 = _mm512_mask_mov_epi64(cv.dd2, m, cvg.dd2);
                cv.a24n = _mm512_mask_mov_epi64(cv.a24n, m, cvg.a24n);
                cv.a24d = _mm512_mask_mov_epi64(cv.a24d, m, cvg.a24d);
            }
        }
        for (lane = 0; lane < 8; lane++)
            kcur[lane]++;

        uedw_pt8 P;
        uedw8_stage1(&cv, &ch, &P, c);
        _mm512_store_si512(st, P.Z);
        for (lane = 0; lane < 8; lane++) {
            if (msk & (1 << lane)) {
                total++;
                uint64_t g = uedw_bingcd64(nA[lane], st[lane]);
                if (g > 1 && g != nA[lane]) {
                    f[idx[lane]] = g;
                    msk &= (uint8_t)~(1 << lane);
                }
            }
        }
        if (B2 > B1 && msk) {
            __m512i acc = lv ? uedw8_stage2_gen(&P, &cv, lv->s2, lv->s2len, c)
                : uedw8_stage2(&P, &cv, B1, c);
            _mm512_store_si512(st, acc);
            for (lane = 0; lane < 8; lane++) {
                if (msk & (1 << lane)) {
                    uint64_t g = uedw_bingcd64(nA[lane], st[lane]);
                    if (g > 1 && g != nA[lane]) {
                        f[idx[lane]] = g;
                        msk &= (uint8_t)~(1 << lane);
                    }
                }
            }
        }
        for (lane = 0; lane < 8; lane++)
            if (kcur[lane] >= (int)curves)
                msk &= (uint8_t)~(1 << lane);
    }
    if (ncurves)
        *ncurves = total;
}

static void uedw8_init_consts(void)
{
}
#endif /* USE_AVX512F */

#ifdef USE_AVX512F
/* B1 / curves as microecm's uecm_dispatch_x8_list */
static void uedw_x8_params(int bits, uint32_t *B1, uint32_t *curves)
{
    if (uedw_params == 2 && bits <= 36) { *B1 = 20; *curves = 32; }
    else if (bits <= 40) { *B1 = 27; *curves = 32; }
    else if (bits <= 44) { *B1 = 47; *curves = 32; }
    else if (bits <= 48) { *B1 = 70; *curves = 32; }
    else if (bits <= 52) { *B1 = 85; *curves = 32; }
    else if (bits <= 58) { *B1 = 125; *curves = 32; }
    else if (bits <= 62) { *B1 = 165; *curves = 42; }
    else { *B1 = 205; *curves = 42; }
}
#endif

void getfactor_uecm_edw_x8_list(uint64_t *q64, uint64_t *f64, uint32_t num_in, uint64_t *pran)
{
    uint32_t i;
    if (num_in == 0)
        return;
#ifdef USE_AVX512F
    int bits = uedw_bits(q64[0]);               // all inputs about the same size
    if (bits <= 52) {
        uint32_t B1, curves;
        uedw8_init_consts();
        uedw_x8_params(bits, &B1, &curves);
        uedw8_list(q64, f64, num_in, B1, 25 * B1, curves, 0);
        return;
    }
#elif defined(USE_AVX2)
    // AVX2 builds: the x4 list (same curves and results)
    getfactor_uecm_edw_x4_list(q64, f64, num_in, pran);
    return;
#endif
    for (i = 0; i < num_in; i++)
        f64[i] = getfactor_uecm_edw(q64[i], 0, pran);
}

uint64_t uecm_edw_x8_run(uint64_t *q64, uint64_t *f64, uint32_t num, uint32_t B1, uint32_t curves)
{
#ifdef USE_AVX512F
    uint64_t total = 0;
    uedw8_init_consts();
    uedw8_list(q64, f64, num, B1, 25 * B1, curves, &total);
    return total;
#else
    (void)q64; (void)f64; (void)num; (void)B1; (void)curves;
    return 0;
#endif
}
