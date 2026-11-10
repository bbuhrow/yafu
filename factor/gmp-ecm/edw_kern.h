/*
edw_kern.h: the curve code of the Edwards micro/tiny ECM implementations,
written once against a small arithmetic interface and included once per
arithmetic type (see tinyecm_edw.c).  The design and the formulas are in
arith/micro_tiny_edwards_plan.md; microecm_edw.c has the same code for
64-bit residues written out directly.

This file is dual-licensed like microecm.c / tinyecm.c (FreeBSD or MPL 2.0).

Before including, define:
  EK_T                residue value type (Montgomery form, in [0, n))
  EK_CTX              context type with EK_T fields one, c6, c10, c16, c32, c48
  EK_FN(x)            name prefix for the functions and types defined here
  EK_MUL(a, b)        Montgomery product; may use the variable c (const EK_CTX *)
  EK_SQR(a)           Montgomery square
  EK_ADD(a, b), EK_SUB(a, b)    modular addition and subtraction
  EK_ZERO             the residue 0
  EK_ACC_DECL         declarations for EK_ACC (may be empty)
  EK_ACC(acc, t)      stage-2 accumulation of the cross product t into acc
                      (a statement; it may use break to leave the loop)
Defined here (with EK_FN names):
  types wpt (E' point), ept (extended Edwards point), cached, curve, pt
  (Montgomery x-only point); functions w_addq, w_dbl, w_addp, gkl_curve,
  curve_from, dbl, add, cache, mdup, uadd, stage1, stage2.
*/

#include <stdint.h>
#include "edw_tiny_chains.h"

#ifndef EDW_KERN_COMMON
#define EDW_KERN_COMMON

#if defined(__GNUC__) || defined(__clang__)
#define EK_UNUSED __attribute__((unused))
#else
#define EK_UNUSED
#endif

typedef struct {
    const int8_t *ch;
    int len, m, top, tz;
} edw_chain_t;

static EK_UNUSED int edw_chain_for(uint32_t B1, edw_chain_t *c)
{
#define EDW_CH(b) case b: c->ch = edw_chain_##b; c->len = EDW_CHAIN_##b##_LEN; \
    c->m = EDW_CHAIN_##b##_M; c->top = EDW_CHAIN_##b##_TOP; c->tz = EDW_CHAIN_##b##_TZ; return 1;
    switch (B1) {
        EDW_CH(27) EDW_CH(47) EDW_CH(70) EDW_CH(85) EDW_CH(125) EDW_CH(165) EDW_CH(205)
    default:
        return 0;
    }
#undef EDW_CH
}

static const int edw_s2map[60] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };
#endif

typedef struct { EK_T X, Y, Z; } EK_FN(wpt);
typedef struct { EK_T X, Y, Z, T; } EK_FN(ept);
typedef struct { EK_T YmX, YpX, T2d, Z2; } EK_FN(cached);
typedef struct {
    EK_FN(ept) P;                   // start point, extended coordinates
    EK_T dn2, dd, dd2;              // 2 dn, dd, 2 dd (d = dn / dd)
    EK_T a24n, a24d;                // A24 = dd / (dd + dn)
} EK_FN(curve);
typedef struct { EK_T X, Z; } EK_FN(pt);

#define MUL EK_MUL
#define SQR EK_SQR
#define ADD EK_ADD
#define SUB EK_SUB

/* E': y^2 = x^3 - x^2 - 9x + 9, projective.  R + Q, Q = (5, 8), R != +-Q */
static EK_UNUSED void EK_FN(w_addq)(EK_FN(wpt) *R, const EK_CTX *c)
{
    EK_T X = R->X, Y = R->Y, Z = R->Z;
    EK_T z2 = ADD(Z, Z);
    EK_T z4 = ADD(z2, z2);
    EK_T u = SUB(ADD(z4, z4), Y);
    EK_T v = SUB(ADD(z4, Z), X);
    EK_T u2 = SQR(u);
    EK_T v2 = SQR(v);
    EK_T v3 = MUL(v2, v);
    EK_T w1 = MUL(v2, X);
    EK_T w2 = MUL(v2, Z);
    w2 = ADD(w2, w2);
    w2 = ADD(w2, w2);
    EK_T A = SUB(SUB(MUL(u2, Z), w1), w2);
    R->X = MUL(v, A);
    R->Y = SUB(MUL(u, SUB(w1, A)), MUL(Y, v3));
    R->Z = MUL(v3, Z);
}

/* 2R */
static EK_UNUSED void EK_FN(w_dbl)(EK_FN(wpt) *R, const EK_CTX *c)
{
    EK_T X = R->X, Y = R->Y, Z = R->Z;
    EK_T X2 = SQR(X), Z2 = SQR(Z), XZ = MUL(X, Z);
    EK_T w = ADD(ADD(X2, X2), X2);
    w = SUB(w, ADD(XZ, XZ));
    EK_T z8 = ADD(Z2, Z2);
    z8 = ADD(z8, z8);
    z8 = ADD(z8, z8);
    w = SUB(w, ADD(z8, Z2));
    EK_T s = MUL(Y, Z);
    s = ADD(s, s);
    EK_T s2 = SQR(s);
    EK_T s3 = MUL(s2, s);
    EK_T A = SUB(MUL(SQR(w), Z), MUL(s2, SUB(ADD(X, X), Z)));
    EK_T Xs2 = MUL(X, s2);
    R->X = MUL(s, A);
    R->Y = SUB(MUL(w, SUB(Xs2, A)), MUL(Y, s3));
    R->Z = MUL(s3, Z);
}

/* R + S, both projective, R != +-S: lambda = u/v, u = Y2 Z1 - Y1 Z2,
   v = X2 Z1 - X1 Z2, W = Z1 Z2; A = u^2 W + v^2 W - v^2 (X1 Z2 + X2 Z1);
   X3 = v A, Y3 = u (X1 Z2 v^2 - A) - Y1 Z2 v^3, Z3 = v^3 W. */
static EK_UNUSED void EK_FN(w_addp)(EK_FN(wpt) *R, const EK_FN(wpt) *S, const EK_CTX *c)
{
    EK_T X1 = R->X, Y1 = R->Y, Z1 = R->Z;
    EK_T X1Z2 = MUL(X1, S->Z), Y1Z2 = MUL(Y1, S->Z);
    EK_T X2Z1 = MUL(S->X, Z1), Y2Z1 = MUL(S->Y, Z1);
    EK_T W = MUL(Z1, S->Z);
    EK_T u = SUB(Y2Z1, Y1Z2);
    EK_T v = SUB(X2Z1, X1Z2);
    EK_T u2 = SQR(u), v2 = SQR(v), v3 = MUL(v2, v);
    EK_T A = SUB(MUL(ADD(u2, v2), W), MUL(v2, ADD(X1Z2, X2Z1)));
    R->X = MUL(v, A);
    R->Y = SUB(MUL(u, SUB(MUL(X1Z2, v2), A)), MUL(Y1Z2, v3));
    R->Z = MUL(v3, W);
}

/* the GKL family (iv) curve of the E' point R (see microecm_edw.c) */
static EK_UNUSED void EK_FN(gkl_curve)(const EK_FN(wpt) *R, EK_FN(curve) *cv, const EK_CTX *c)
{
    EK_T t = ADD(R->X, R->Z);
    t = ADD(t, t);
    EK_T a = ADD(t, t);
    EK_T z4 = ADD(R->Z, R->Z);
    z4 = ADD(z4, z4);
    EK_T b = SUB(R->Y, z4);
    EK_T a2 = SQR(a), b2 = SQR(b), ab = MUL(a, b);
    EK_T b2x = ADD(b, b);
    EK_T xn = MUL(MUL(SUB(a, b2x), ADD(a, b2x)), ADD(a, b));
    xn = MUL(xn, b);
    xn = ADD(xn, xn);
    EK_T aa4 = MUL(a, ADD(a, ADD(b2x, b2x)));
    EK_T b2_4 = ADD(b2, b2);
    b2_4 = ADD(b2_4, b2_4);
    EK_T xd = MUL(aa4, ADD(ADD(a2, ADD(ab, ab)), b2_4));
    EK_T t1 = SQR(a2), t2 = MUL(a2, ab), t3 = MUL(a2, b2), t4 = MUL(ab, b2), t5 = SQR(b2);
    EK_T q = ADD(t2, t3);
    q = ADD(t1, ADD(q, q));
    EK_T r = ADD(t4, ADD(t5, t5));
    r = ADD(r, r);
    r = ADD(r, r);
    r = ADD(r, r);
    EK_T yd = MUL(aa4, ADD(q, r));
    EK_T b3 = MUL(b2, b), b4 = SQR(b2), b5 = MUL(b4, b), b6 = SQR(b3);
    EK_T h = ADD(a, MUL(b, c->c6));
    h = ADD(MUL(h, a), MUL(b2, c->c10));
    h = SUB(MUL(h, a), MUL(b3, c->c16));
    h = SUB(MUL(h, a), MUL(b4, c->c48));
    h = SUB(MUL(h, a), MUL(b5, c->c32));
    EK_T yn = SUB(MUL(h, a), MUL(b6, c->c32));
    EK_T en2 = SQR(aa4), ed = SUB(a2, b2_4), ed2 = SQR(ed);
    EK_T dn = SUB(EK_ZERO, SQR(en2));
    EK_T dd = SQR(ed2);
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

/* a curve from Montgomery-form constants X, Y, Z, T, dn, dd */
static EK_UNUSED void EK_FN(curve_from)(const EK_T *v, EK_FN(curve) *cv, const EK_CTX *c)
{
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

/* a = -1 doubling, signs folded: 3M + 4S (+1M for T) */
static EK_UNUSED void EK_FN(dbl)(EK_FN(ept) *R, const EK_FN(ept) *P, int want_t, const EK_CTX *c)
{
    EK_T A = SQR(P->X);
    EK_T B = SQR(P->Y);
    EK_T C = SQR(P->Z);
    EK_T S = SQR(ADD(P->X, P->Y));
    C = ADD(C, C);
    EK_T H = ADD(A, B);
    EK_T E = SUB(H, S);
    EK_T G = SUB(A, B);
    EK_T F = ADD(C, G);
    R->X = MUL(E, F);
    R->Y = MUL(G, H);
    R->Z = MUL(F, G);
    if (want_t)
        R->T = MUL(E, H);
}

/* P +- Q, Q cached: 7M (+1M for T) */
static EK_UNUSED void EK_FN(add)(EK_FN(ept) *R, const EK_FN(ept) *P, const EK_FN(cached) *Q,
    int neg, int want_t, const EK_CTX *c)
{
    EK_T A = MUL(SUB(P->Y, P->X), neg ? Q->YpX : Q->YmX);
    EK_T B = MUL(ADD(P->Y, P->X), neg ? Q->YmX : Q->YpX);
    EK_T C = MUL(P->T, Q->T2d);
    EK_T D = MUL(P->Z, Q->Z2);
    EK_T E = SUB(B, A);
    EK_T H = ADD(B, A);
    EK_T F = neg ? ADD(D, C) : SUB(D, C);
    EK_T G = neg ? SUB(D, C) : ADD(D, C);
    R->X = MUL(E, F);
    R->Y = MUL(G, H);
    R->Z = MUL(F, G);
    if (want_t)
        R->T = MUL(E, H);
}

static EK_UNUSED void EK_FN(cache)(EK_FN(cached) *q, const EK_FN(ept) *P, const EK_FN(curve) *cv,
    const EK_CTX *c)
{
    q->YmX = MUL(SUB(P->Y, P->X), cv->dd);
    q->YpX = MUL(ADD(P->Y, P->X), cv->dd);
    q->T2d = MUL(P->T, cv->dn2);
    q->Z2 = MUL(P->Z, cv->dd2);
}

/* Montgomery x-only doubling with projective A24: 2S + 4M */
static EK_UNUSED void EK_FN(mdup)(const EK_FN(pt) *P, EK_FN(pt) *R, const EK_FN(curve) *cv,
    const EK_CTX *c)
{
    EK_T U = SQR(SUB(P->X, P->Z));
    EK_T V = SQR(ADD(P->X, P->Z));
    EK_T w = SUB(V, U);
    EK_T Ub = MUL(U, cv->a24d);
    R->X = MUL(Ub, V);
    R->Z = MUL(w, ADD(Ub, MUL(w, cv->a24n)));
}

/* differential addition P1 + P2, Pin = P1 - P2 */
static EK_UNUSED void EK_FN(uadd)(const EK_FN(pt) P1, const EK_FN(pt) P2, const EK_FN(pt) Pin,
    EK_FN(pt) *Pout, const EK_CTX *c)
{
    EK_T diff1 = SUB(P1.X, P1.Z);
    EK_T sum1 = ADD(P1.X, P1.Z);
    EK_T diff2 = SUB(P2.X, P2.Z);
    EK_T sum2 = ADD(P2.X, P2.Z);
    EK_T tt1 = MUL(diff1, sum2);
    EK_T tt2 = MUL(sum1, diff2);
    EK_T tt3 = ADD(tt1, tt2);
    EK_T tt4 = SUB(tt1, tt2);
    tt1 = SQR(tt3);
    tt2 = SQR(tt4);
    EK_T tmpx = MUL(tt1, Pin.Z);
    EK_T tmpz = MUL(tt2, Pin.X);
    Pout->X = tmpx;
    Pout->Z = tmpz;
}

/* stage 1: the chain on the odd part of s, the hand-off to the Montgomery
   image, then the power of 2 */
static EK_UNUSED void EK_FN(stage1)(const EK_FN(curve) *cv, const edw_chain_t *ch, EK_FN(pt) *Pm,
    const EK_CTX *c)
{
    EK_FN(cached) tab[16];
    EK_FN(ept) ext[16];
    int k, i;
    ext[0] = cv->P;
    EK_FN(cache)(&tab[0], &ext[0], cv, c);
    if (ch->m > 1) {
        EK_FN(ept) P2;
        EK_FN(cached) c2;
        EK_FN(dbl)(&P2, &cv->P, 1, c);
        EK_FN(cache)(&c2, &P2, cv, c);
        for (k = 1; k <= (ch->m - 1) / 2; k++) {
            EK_FN(add)(&ext[k], &ext[k - 1], &c2, 0, 1, c);
            EK_FN(cache)(&tab[k], &ext[k], cv, c);
        }
    }
    EK_FN(ept) R = ext[(ch->top - 1) / 2];
    for (i = 0; i < ch->len; i++) {
        int dg = ch->ch[i];
        if (dg == 0) {
            int t = (i + 1 < ch->len) && (ch->ch[i + 1] != 0);
            EK_FN(dbl)(&R, &R, t, c);
        } else if (dg > 0) {
            EK_FN(add)(&R, &R, &tab[(dg - 1) >> 1], 0, 0, c);
        } else {
            EK_FN(add)(&R, &R, &tab[(-dg - 1) >> 1], 1, 0, c);
        }
    }
    Pm->X = ADD(R.Z, R.Y);
    Pm->Z = SUB(R.Z, R.Y);
    for (i = 0; i < ch->tz; i++)
        EK_FN(mdup)(Pm, Pm, cv, c);
}

static EK_UNUSED EK_T EK_FN(cross)(const EK_FN(pt) *R, EK_T Rprod, const EK_FN(pt) *Pb, EK_T Pbprod,
    const EK_CTX *c)
{
    EK_T tt1 = SUB(R->X, Pb->X);
    EK_T tt2 = ADD(R->Z, Pb->Z);
    EK_T tt3 = MUL(tt1, tt2);
    tt1 = ADD(tt3, Pbprod);
    return SUB(tt1, Rprod);
}

/* microecm / tinyecm's pairing stage 2 (D = 60), projective A24, with the
   pair table barray / numb of the caller and EK_ACC's accumulation rule */
static EK_UNUSED EK_T EK_FN(stage2)(const EK_FN(pt) *P, const EK_FN(curve) *cv, uint32_t stg1_max,
    const uint8_t *barray, int numb, const EK_CTX *c)
{
    int b, i, j, k;
    EK_FN(pt) Pa, Pd, Pad, Pb[18], pt4, pt5, pt6;
    EK_T Pbprod[18];
    const int D = 60;
    EK_ACC_DECL

    Pb[1] = *P;
    Pbprod[1] = MUL(Pb[1].X, Pb[1].Z);
    EK_FN(mdup)(P, &Pb[2], cv, c);
    EK_FN(uadd)(Pb[1], Pb[2], Pb[1], &Pb[3], c);
    EK_FN(mdup)(&Pb[3], &pt6, cv, c);
    EK_FN(uadd)(Pb[3], Pb[2], Pb[1], &pt5, c);
    Pb[3] = pt5;
    EK_FN(uadd)(pt6, pt5, Pb[1], &Pb[4], c);
    i = 3; k = 4; j = 5;
    while ((j + 12) < D) {
        EK_FN(uadd)(pt6, Pb[k], Pb[i], &Pb[edw_s2map[j + 12]], c);
        i = k;
        k = edw_s2map[j + 12];
        j += 6;
    }
    EK_FN(uadd)(pt6, Pb[1], pt5, &Pb[3], c);
    i = 1; k = 3; j = 1;
    while ((j + 12) < D) {
        EK_FN(uadd)(pt6, Pb[k], Pb[i], &Pb[edw_s2map[j + 12]], c);
        i = k;
        k = edw_s2map[j + 12];
        j += 6;
    }
    EK_FN(uadd)(Pb[9], Pb[10], Pb[2], &Pd, c);
    EK_FN(mdup)(&Pb[2], &pt4, cv, c);
    for (i = 3; i < 18; i++)
        Pbprod[i] = MUL(Pb[i].X, Pb[i].Z);
    EK_FN(uadd)(Pb[edw_s2map[17]], Pb[edw_s2map[13]], pt4, &Pad, c);
    EK_FN(uadd)(Pd, Pad, Pad, &Pa, c);
    EK_FN(pt) pt90 = Pa;
    EK_FN(pt) pt60 = Pd;
    EK_FN(uadd)(Pa, Pad, Pd, &Pd, c);
    EK_FN(uadd)(Pd, Pad, Pa, &Pa, c);

    EK_T acc = c->one;
    if (stg1_max < 70) {
        static const int steps27[16] = { 59,53,49,47,43,37,31,29,23,19,17,11,7,1,13,41 };
        static const int steps47[15] = { 43,37,31,29,23,19,17,11,7,1,13,41,47,49,59 };
        const int *steps = (stg1_max == 27) ? steps27 : steps47;
        int numsteps = (stg1_max == 27) ? 16 : 15;
        EK_T pt90prod = MUL(pt90.X, pt90.Z);
        for (i = 0; i < numsteps; i++) {
            b = steps[i];
            EK_T t_ = EK_FN(cross)(&pt90, pt90prod, &Pb[edw_s2map[b]], Pbprod[edw_s2map[b]], c);
            EK_ACC(acc, t_)
        }
    } else if (stg1_max == 70) {
        static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
        EK_T pdprod = MUL(Pd.X, Pd.Z);
        for (i = 0; i < 15; i++) {
            b = steps[i];
            EK_T t_ = EK_FN(cross)(&Pd, pdprod, &Pb[edw_s2map[b]], Pbprod[edw_s2map[b]], c);
            EK_ACC(acc, t_)
        }
    } else if (stg1_max == 165) {
        EK_FN(mdup)(&pt90, &Pa, cv, c);
        Pad = pt60;
    } else if (stg1_max == 205) {
        EK_FN(uadd)(Pd, pt90, Pad, &Pa, c);
        Pad = pt90;
    }
    EK_T Paprod = MUL(Pa.X, Pa.Z);
    for (i = 0; i < numb; i++) {
        if (barray[i] == 0) {
            EK_FN(pt) point = Pa;
            EK_FN(uadd)(Pa, Pd, Pad, &Pa, c);
            Pad = point;
            Paprod = MUL(Pa.X, Pa.Z);
            i++;
        }
        b = barray[i];
        EK_T t_ = EK_FN(cross)(&Pa, Paprod, &Pb[edw_s2map[b]], Pbprod[edw_s2map[b]], c);
        EK_ACC(acc, t_)
    }
    return acc;
}

#undef MUL
#undef SQR
#undef ADD
#undef SUB
