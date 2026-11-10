/*
 * Phase 0 microbenchmark (arith/micro_tiny_edwards_plan.md): microecm's
 * 64-bit scalar stage 1, PRAC on Suyama curves (the existing code) vs an
 * Edwards a = -1 extended-coordinate prototype on the chains from
 * tiny_chains.py (edw_tiny_chains.h).  Also 2-way interleaved versions of
 * both (two independent curves per call), which fill the multiplier's
 * latency the way a list / multi-curve driver could.
 *
 * Includes microecm.c, so it uses exactly that file's arithmetic and PRAC.
 *
 *   tiny_bench -t     correctness: the Edwards chain + hand-off must give
 *                     the same Montgomery x as PRAC on the Montgomery image
 *   tiny_bench [ncurves]   timings, ns per curve, per level
 */
#define _GNU_SOURCE
#include "../../factor/gmp-ecm/microecm.c"
#include "edw_tiny_chains.h"
#include <stdio.h>
#include <time.h>

static uint64_t xs = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd64(void)
{
    xs ^= xs << 13; xs ^= xs >> 7; xs ^= xs << 17;
    return xs;
}

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

typedef struct { uint64_t X, Y, Z, T; } ept;
typedef struct { uint64_t YmX, YpX, T2d, Z2; } ecached;

typedef struct {
    uint64_t n, rho, one, Rsqr;
} mctx;

static void mctx_init(mctx *c, uint64_t n)
{
    c->n = n;
    c->rho = uecm_multiplicative_inverse(n);
    c->one = ((uint64_t)0 - n) % n;
    uint64_t two = uecm_addmod(c->one, c->one, n);
    uint64_t four = uecm_addmod(two, two, n);
    uint64_t eight = uecm_addmod(four, four, n);
    uint64_t sixteen = uecm_addmod(eight, eight, n);
    uint64_t t = uecm_sqrredc(sixteen, n, c->rho);
    t = uecm_sqrredc(t, n, c->rho);
    t = uecm_sqrredc(t, n, c->rho);
    c->Rsqr = uecm_sqrredc(t, n, c->rho);
}

#define MUL(a, b) uecm_mulredc((a), (b), n, rho)
#define SQR(a) uecm_sqrredc((a), n, rho)
#define ADD(a, b) uecm_addmod((a), (b), n)
#define SUB(a, b) uecm_submod((a), (b), n)

/* doubling, a = -1, signs folded (avxecm_edwards_plan.md C0): 3M + 4S (+1M T) */
MICRO_ECM_FORCE_INLINE static void edw_dbl(ept *R, const ept *P, int want_t,
    uint64_t n, uint64_t rho)
{
    uint64_t A = SQR(P->X);
    uint64_t B = SQR(P->Y);
    uint64_t C = SQR(P->Z);
    uint64_t S = SQR(ADD(P->X, P->Y));
    C = ADD(C, C);
    uint64_t H = ADD(A, B);
    uint64_t E = SUB(H, S);
    uint64_t G = SUB(A, B);
    uint64_t Fv = ADD(C, G);
    R->X = MUL(E, Fv);
    R->Y = MUL(G, H);
    R->Z = MUL(Fv, G);
    if (want_t) R->T = MUL(E, H);
}

/* P +- Q, Q cached: 7M (+1M T) */
MICRO_ECM_FORCE_INLINE static void edw_add(ept *R, const ept *P, const ecached *Q,
    int neg, int want_t, uint64_t n, uint64_t rho)
{
    uint64_t A = MUL(SUB(P->Y, P->X), neg ? Q->YpX : Q->YmX);
    uint64_t B = MUL(ADD(P->Y, P->X), neg ? Q->YmX : Q->YpX);
    uint64_t C = MUL(P->T, Q->T2d);
    uint64_t D = MUL(P->Z, Q->Z2);
    uint64_t E = SUB(B, A);
    uint64_t H = ADD(B, A);
    uint64_t Fv = neg ? ADD(D, C) : SUB(D, C);
    uint64_t G = neg ? SUB(D, C) : ADD(D, C);
    R->X = MUL(E, Fv);
    R->Y = MUL(G, H);
    R->Z = MUL(Fv, G);
    if (want_t) R->T = MUL(E, H);
}

/* cached form scaled by dd (d = dn/dd, no inversion): 4M */
MICRO_ECM_FORCE_INLINE static void edw_cache(ecached *c, const ept *P, uint64_t dn2,
    uint64_t dd, uint64_t dd2, uint64_t n, uint64_t rho)
{
    c->YmX = MUL(SUB(P->Y, P->X), dd);
    c->YpX = MUL(ADD(P->Y, P->X), dd);
    c->T2d = MUL(P->T, dn2);
    c->Z2 = MUL(P->Z, dd2);
}

/* Montgomery x-only doubling with projective A24 = a24n / a24d: 2S + 4M */
MICRO_ECM_FORCE_INLINE static void mdup_proj(uint64_t *X, uint64_t *Z, uint64_t a24n,
    uint64_t a24d, uint64_t n, uint64_t rho)
{
    uint64_t U = SQR(SUB(*X, *Z));
    uint64_t V = SQR(ADD(*X, *Z));
    uint64_t w = SUB(V, U);
    uint64_t Ub = MUL(U, a24d);
    *X = MUL(Ub, V);
    *Z = MUL(w, ADD(Ub, MUL(w, a24n)));
}

typedef struct {
    const int8_t *ch;
    int len, m, top, tz;
} chain_t;

static chain_t chain_for(int B1)
{
    chain_t c;
#define CH(b) case b: c.ch = edw_chain_##b; c.len = EDW_CHAIN_##b##_LEN; \
    c.m = EDW_CHAIN_##b##_M; c.top = EDW_CHAIN_##b##_TOP; c.tz = EDW_CHAIN_##b##_TZ; break;
    switch (B1) {
        CH(27) CH(47) CH(70) CH(85) CH(125) CH(165) CH(205)
    default: c.ch = 0; c.len = 0; break;
    }
    return c;
}

typedef struct {
    ept P;          /* start point, extended, Montgomery form */
    uint64_t dn2, dd, dd2;     /* 2 dn, dd, 2 dd */
    uint64_t a24n, a24d;       /* A24 = dd / (dd + dn) */
} ecurve;

/* Edwards stage 1: chain on the odd part, hand-off, tz Montgomery doublings */
static void edw_stage1(const ecurve *cv, const chain_t *c, uint64_t *Xm, uint64_t *Zm,
    uint64_t n, uint64_t rho)
{
    ecached tab[16];
    ept ext[16];
    int k, i;
    ext[0] = cv->P;
    edw_cache(&tab[0], &ext[0], cv->dn2, cv->dd, cv->dd2, n, rho);
    if (c->m > 1) {
        ept P2;
        ecached c2;
        edw_dbl(&P2, &cv->P, 1, n, rho);
        edw_cache(&c2, &P2, cv->dn2, cv->dd, cv->dd2, n, rho);
        for (k = 1; k <= (c->m - 1) / 2; k++) {
            edw_add(&ext[k], &ext[k - 1], &c2, 0, 1, n, rho);
            edw_cache(&tab[k], &ext[k], cv->dn2, cv->dd, cv->dd2, n, rho);
        }
    }
    ept R = ext[(c->top - 1) / 2];
    for (i = 0; i < c->len; i++) {
        int dg = c->ch[i];
        if (dg == 0) {
            int t = (i + 1 < c->len) && (c->ch[i + 1] != 0);
            edw_dbl(&R, &R, t, n, rho);
        } else if (dg > 0) {
            edw_add(&R, &R, &tab[(dg - 1) >> 1], 0, 0, n, rho);
        } else {
            edw_add(&R, &R, &tab[(-dg - 1) >> 1], 1, 0, n, rho);
        }
    }
    uint64_t X = ADD(R.Z, R.Y), Z = SUB(R.Z, R.Y);
    for (i = 0; i < c->tz; i++)
        mdup_proj(&X, &Z, cv->a24n, cv->a24d, n, rho);
    *Xm = X;
    *Zm = Z;
}

/* ---- 2-way interleaved versions: two independent curves per call ---- */
typedef struct { uecm_pt p[2]; } pt2;

#define FOR2 for (int j_ = 0; j_ < 2; j_++)

MICRO_ECM_FORCE_INLINE static void uadd2(const uint64_t *rho, const uint64_t *n,
    const pt2 *P1, const pt2 *P2, const pt2 *Pin, pt2 *Pout)
{
    pt2 o;
    FOR2 uecm_uadd(rho[j_], n[j_], P1->p[j_], P2->p[j_], Pin->p[j_], &o.p[j_]);
    *Pout = o;
}

MICRO_ECM_FORCE_INLINE static void udup2(const uint64_t *s, const uint64_t *rho,
    const uint64_t *n, pt2 *P)
{
    FOR2 {
        uint64_t d1 = uecm_submod(P->p[j_].X, P->p[j_].Z, n[j_]);
        uint64_t s1 = uecm_addmod(P->p[j_].X, P->p[j_].Z, n[j_]);
        uecm_udup(s[j_], rho[j_], n[j_], s1, d1, &P->p[j_]);
    }
}

static void uprac2(const uint64_t *rho, const uint64_t *n, pt2 *P, uint64_t c, double v,
    const uint64_t *s)
{
    uint64_t d, e, r;
    int i;
    int shift = _trail_zcnt64(c);
    c = c >> shift;
    d = c;
    r = (uint64_t)((double)d * v + 0.5);
    d = c - r;
    e = 2 * r - c;
    pt2 pt1, pt2_, pt3, pt4, pt5, swp;
    pt1 = pt2_ = pt3 = *P;
    udup2(s, rho, n, &pt1);
    while (d != e) {
        if (d < e) { r = d; d = e; e = r; swp = pt1; pt1 = pt2_; pt2_ = swp; }
        if (d - e <= e / 4 && ((d + e) % 3) == 0) {
            d = (2 * d - e) / 3; e = (e - d) / 2;
            uadd2(rho, n, &pt1, &pt2_, &pt3, &pt4);
            uadd2(rho, n, &pt4, &pt1, &pt2_, &pt5);
            uadd2(rho, n, &pt2_, &pt4, &pt1, &pt2_);
            swp = pt1; pt1 = pt5; pt5 = swp;
        } else if (d - e <= e / 4 && (d - e) % 6 == 0) {
            d = (d - e) / 2;
            uadd2(rho, n, &pt1, &pt2_, &pt3, &pt2_);
            udup2(s, rho, n, &pt1);
        } else if ((d + 3) / 4 <= e) {
            d -= e;
            uadd2(rho, n, &pt2_, &pt1, &pt3, &pt4);
            swp = pt2_; pt2_ = pt4; pt4 = pt3; pt3 = swp;
        } else if ((d + e) % 2 == 0) {
            d = (d - e) / 2;
            uadd2(rho, n, &pt2_, &pt1, &pt3, &pt2_);
            udup2(s, rho, n, &pt1);
        } else if (d % 2 == 0) {
            d /= 2;
            uadd2(rho, n, &pt3, &pt1, &pt2_, &pt3);
            udup2(s, rho, n, &pt1);
        } else if (d % 3 == 0) {
            d = d / 3 - e;
            pt4 = pt1; udup2(s, rho, n, &pt4);
            uadd2(rho, n, &pt1, &pt2_, &pt3, &pt5);
            uadd2(rho, n, &pt4, &pt1, &pt1, &pt1);
            uadd2(rho, n, &pt4, &pt5, &pt3, &pt4);
            swp = pt3; pt3 = pt2_; pt2_ = pt4; pt4 = swp;
        } else if ((d + e) % 3 == 0) {
            d = (d - 2 * e) / 3;
            uadd2(rho, n, &pt1, &pt2_, &pt3, &pt4);
            uadd2(rho, n, &pt4, &pt1, &pt2_, &pt2_);
            pt4 = pt1; udup2(s, rho, n, &pt4);
            uadd2(rho, n, &pt1, &pt4, &pt1, &pt1);
        } else if ((d - e) % 3 == 0) {
            d = (d - e) / 3;
            uadd2(rho, n, &pt1, &pt2_, &pt3, &pt4);
            uadd2(rho, n, &pt3, &pt1, &pt2_, &pt3);
            swp = pt2_; pt2_ = pt4; pt4 = swp;
            pt4 = pt1; udup2(s, rho, n, &pt4);
            uadd2(rho, n, &pt1, &pt4, &pt1, &pt1);
        } else {
            e /= 2;
            uadd2(rho, n, &pt3, &pt2_, &pt1, &pt3);
            udup2(s, rho, n, &pt2_);
        }
    }
    uadd2(rho, n, &pt1, &pt2_, &pt3, P);
    for (i = 0; i < shift; i++) udup2(s, rho, n, P);
}

static void uprac_steps2(const uint64_t *rho, const uint64_t *n, pt2 *P, const uint8_t *steps,
    int nsteps, const uint64_t *s)
{
    pt2 pt1, pt2_, pt3, pt4;
    pt1 = pt2_ = pt3 = *P;
    for (int i = 0; i < nsteps; i++) {
        if (steps[i] == 0) {
            pt1 = pt2_ = pt3 = *P;
            udup2(s, rho, n, &pt1);
        } else if (steps[i] == 3) {
            uadd2(rho, n, &pt2_, &pt1, &pt3, &pt4);
            pt2 t = pt1; pt1 = pt4; pt4 = pt3; pt3 = pt2_; pt2_ = t;
        } else if (steps[i] == 4) {
            uadd2(rho, n, &pt2_, &pt1, &pt3, &pt4);
            pt2 t = pt2_; pt2_ = pt4; pt4 = pt3; pt3 = t;
        } else if (steps[i] == 5) {
            pt4 = pt1;
            uadd2(rho, n, &pt2_, &pt1, &pt3, &pt2_);
            pt1 = pt4; udup2(s, rho, n, &pt1);
        } else if (steps[i] == 6) {
            uadd2(rho, n, &pt1, &pt2_, &pt3, P);
        }
    }
}

static const uint8_t steps70[116] = {
    0,6,0,6,0,6,0,4,6,0,4,6,0,4,4,6,0,4,4,6,0,5,4,6,0,3,3,4,6,0,3,5,
    4,6,0,3,4,3,4,6,0,5,5,4,6,0,5,3,3,4,6,0,3,3,4,3,4,6,0,5,3,3,3,3,
    3,3,3,3,4,3,3,4,6,0,5,4,3,3,4,6,0,3,4,3,5,4,6,0,5,3,3,3,4,6,0,5,
    4,3,5,4,6,0,5,5,3,3,4,6,0,4,3,3,3,5,4,6 };
static const uint8_t steps85[146] = {
    0,6,0,6,0,6,0,6,0,4,6,0,4,6,0,4,4,6,0,4,4,6,0,5,4,6,0,3,3,4,
    6,0,3,5,4,6,0,3,4,3,4,6,0,5,5,4,6,0,5,3,3,4,6,0,3,3,4,3,4,6,
    0,4,3,4,3,5,3,3,3,3,3,3,3,3,4,6,0,3,3,3,3,3,3,3,3,3,4,3,4,3,
    4,6,0,3,4,3,5,4,6,0,5,3,3,3,4,6,0,5,4,3,5,4,6,0,4,3,3,3,5,4,
    6,0,4,3,5,3,3,4,6,0,3,3,3,3,5,4,6,0,3,3,3,4,3,3,4,6 };

#define G_ 0.618033988749894903
typedef struct { uint64_t c; double v; } pc;
static const pc l27[] = { {3,G_},{3,G_},{3,G_},{5,G_},{5,G_},{7,G_},{11,0.580178728295464130},
    {13,G_},{17,G_},{19,G_},{23,0.522786351415446049} };
static const pc l47[] = { {3,G_},{3,G_},{3,G_},{3,G_},{5,G_},{5,G_},{7,G_},{11,0.580178728295464130},
    {13,G_},{17,G_},{19,G_},{23,0.522786351415446049},{29,0.548409048446403258},{31,G_},
    {37,0.580178728295464130},{41,0.548409048446403258},{43,G_} };
static const pc t125[] = { {5,G_},{11,0.580178728295464130},{89,G_},{97,0.723606797749978936},
    {101,0.556250337855490828},{107,0.580178728295464130},{109,0.548409048446403258},{113,G_} };
static const pc t125b[] = { {61,0.522786351415446049},{103,0.632839806088706269} };
static const pc t165[] = { {7747,0.552188778811121},{131,G_},{14111,0.632839806088706},
    {20989,0.620181980807415},{157,0.640157392785047019},{163,0.551390822543526449} };
static const pc t165b[] = { {149,0.580178728295464130} };
static const pc t205[] = { {13,G_},{167,0.580178728295464130},{173,0.612429949509495031},{179,G_},
    {181,0.551390822543526449},{191,G_},{193,G_},{29353,0.580178728295464},{199,0.551390822543526449} };
static const pc t85[] = { {61,0.522786351415446049} };

#define RUNL(L) for (int q_ = 0; q_ < (int)(sizeof(L) / sizeof(L[0])); q_++) \
    uprac2(rho, n, P, L[q_].c, L[q_].v, s);

static void stage1_x2(const uint64_t *rho, const uint64_t *n, pt2 *P, int B1, const uint64_t *s)
{
    uint64_t q = 2;
    while (q < (uint64_t)B1 * 4) { udup2(s, rho, n, P); q *= 2; }
    if (B1 == 27) { RUNL(l27) }
    else if (B1 == 47) { RUNL(l47) }
    else if (B1 == 70) uprac_steps2(rho, n, P, steps70, 116, s);
    else {
        uprac_steps2(rho, n, P, steps85, 146, s);
        if (B1 == 85) { RUNL(t85) }
        else {
            RUNL(t125)
            if (B1 == 125) { RUNL(t125b) }
            else {
                RUNL(t165)
                if (B1 == 165) { RUNL(t165b) } else { RUNL(t205) }
            }
        }
    }
}

/* Edwards x2: the same chain on two curves, interleaved op by op */
static void edw_stage1_x2(const ecurve *cv, const chain_t *c, uint64_t *Xm, uint64_t *Zm,
    const uint64_t *n, const uint64_t *rho)
{
    ecached tab[2][16];
    ept ext[2][16];
    ept R[2];
    int k, i;
    FOR2 {
        ext[j_][0] = cv[j_].P;
        edw_cache(&tab[j_][0], &ext[j_][0], cv[j_].dn2, cv[j_].dd, cv[j_].dd2, n[j_], rho[j_]);
    }
    if (c->m > 1) {
        ept P2[2];
        ecached c2[2];
        FOR2 {
            edw_dbl(&P2[j_], &cv[j_].P, 1, n[j_], rho[j_]);
            edw_cache(&c2[j_], &P2[j_], cv[j_].dn2, cv[j_].dd, cv[j_].dd2, n[j_], rho[j_]);
        }
        for (k = 1; k <= (c->m - 1) / 2; k++) FOR2 {
            edw_add(&ext[j_][k], &ext[j_][k - 1], &c2[j_], 0, 1, n[j_], rho[j_]);
            edw_cache(&tab[j_][k], &ext[j_][k], cv[j_].dn2, cv[j_].dd, cv[j_].dd2, n[j_], rho[j_]);
        }
    }
    FOR2 R[j_] = ext[j_][(c->top - 1) / 2];
    for (i = 0; i < c->len; i++) {
        int dg = c->ch[i];
        if (dg == 0) {
            int t = (i + 1 < c->len) && (c->ch[i + 1] != 0);
            FOR2 edw_dbl(&R[j_], &R[j_], t, n[j_], rho[j_]);
        } else if (dg > 0) {
            FOR2 edw_add(&R[j_], &R[j_], &tab[j_][(dg - 1) >> 1], 0, 0, n[j_], rho[j_]);
        } else {
            FOR2 edw_add(&R[j_], &R[j_], &tab[j_][(-dg - 1) >> 1], 1, 0, n[j_], rho[j_]);
        }
    }
    FOR2 {
        uint64_t X = uecm_addmod(R[j_].Z, R[j_].Y, n[j_]), Z = uecm_submod(R[j_].Z, R[j_].Y, n[j_]);
        Xm[j_] = X; Zm[j_] = Z;
    }
    for (i = 0; i < c->tz; i++) FOR2
        mdup_proj(&Xm[j_], &Zm[j_], cv[j_].a24n, cv[j_].a24d, n[j_], rho[j_]);
}

/* ---- curve setup ---- */
static uint64_t mod_inv(uint64_t a, uint64_t n)
{
    uint64_t g;
    return uecm_modinv_64(a % n, n, &g);
}

static uint64_t to_m(const mctx *c, uint64_t a) { return uecm_mulredc(a % c->n, c->Rsqr, c->n, c->rho); }
static uint64_t from_m(const mctx *c, uint64_t a) { return uecm_mulredc(a, 1, c->n, c->rho); }

static uint64_t mulmod_plain(uint64_t a, uint64_t b, uint64_t n)
{
    return (uint64_t)(((unsigned __int128)a * b) % n);
}

/* random a = -1 Edwards curve through a random point: d = (y^2 - x^2 - 1)/(x^2 y^2),
   kept as dn / dd with dd random (projective), as the table curves would be */
static void random_ecurve(const mctx *c, ecurve *cv, uint64_t *dplain)
{
    uint64_t n = c->n;
    uint64_t x = rnd64() % n, y = rnd64() % n;
    uint64_t x2 = mulmod_plain(x, x, n), y2 = mulmod_plain(y, y, n);
    uint64_t num = uecm_submod(y2, x2, n);
    num = uecm_submod(num, 1, n);
    uint64_t den = mulmod_plain(x2, y2, n);
    uint64_t d = mulmod_plain(num, mod_inv(den, n), n);
    *dplain = d;
    uint64_t dd = rnd64() % n | 1;
    uint64_t dn = mulmod_plain(d, dd, n);
    uint64_t rho = c->rho;
    cv->P.X = to_m(c, x);
    cv->P.Y = to_m(c, y);
    cv->P.Z = c->one;
    cv->P.T = MUL(cv->P.X, cv->P.Y);
    cv->dd = to_m(c, dd);
    cv->dd2 = ADD(cv->dd, cv->dd);
    uint64_t dnm = to_m(c, dn);
    cv->dn2 = ADD(dnm, dnm);
    cv->a24n = cv->dd;
    cv->a24d = ADD(cv->dd, dnm);
}

static int is_prime64(uint64_t n)
{
    if (n < 2 || n % 2 == 0) return n == 2;
    uint64_t dd = n - 1; int r = 0;
    while (dd % 2 == 0) { dd /= 2; r++; }
    static const uint64_t bases[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    for (int i = 0; i < 12; i++) {
        uint64_t a = bases[i] % n, x = 1, b = dd, e = a;
        if (a == 0) continue;
        while (b) { if (b & 1) x = mulmod_plain(x, e, n); e = mulmod_plain(e, e, n); b >>= 1; }
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int j = 1; j < r; j++) { x = mulmod_plain(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}

static const int LV[7] = { 27, 47, 70, 85, 125, 165, 205 };

static int test(int ntests)
{
    int bad = 0, tot = 0;
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        for (int t = 0; t < ntests; t++) {
            uint64_t n;
            do { n = (rnd64() >> (t & 15)) | 1; } while (n < 1000 || !is_prime64(n));
            mctx c; mctx_init(&c, n);
            uint64_t rho = c.rho;
            ecurve cv; uint64_t d;
            random_ecurve(&c, &cv, &d);
            uint64_t Xe, Ze;
            edw_stage1(&cv, &ch, &Xe, &Ze, n, rho);
            /* Montgomery image + PRAC */
            uecm_pt P;
            P.X = ADD(cv.P.Z, cv.P.Y);
            P.Z = SUB(cv.P.Z, cv.P.Y);
            uint64_t a24 = mulmod_plain(1, mod_inv((1 + d) % n, n), n);   /* 1/(1+d) */
            uint64_t s = to_m(&c, a24);
            uecm_stage1(rho, n, &P, (uint64_t)B1, s);
            /* x2 versions */
            pt2 P2; uint64_t n2[2] = { n, n }, r2[2] = { rho, rho }, s2[2] = { s, s };
            P2.p[0].X = P2.p[1].X = ADD(cv.P.Z, cv.P.Y);
            P2.p[0].Z = P2.p[1].Z = SUB(cv.P.Z, cv.P.Y);
            stage1_x2(r2, n2, &P2, B1, s2);
            ecurve cv2[2] = { cv, cv };
            uint64_t Xe2[2], Ze2[2];
            edw_stage1_x2(cv2, &ch, Xe2, Ze2, n2, r2);
            int ok = MUL(Xe, P.Z) == MUL(P.X, Ze)
                && P2.p[1].X == P.X && P2.p[1].Z == P.Z
                && Xe2[1] == Xe && Ze2[1] == Ze;
            if (!ok && bad < 5)
                printf("mismatch B1 %d n %llu\n", B1, (unsigned long long)n);
            bad += !ok;
            tot++;
        }
    }
    printf("test: %d of %d mismatches\n", bad, tot);
    return bad;
}

#define NBATCH 1024

static void bench(int ncurves, int nbits)
{
    static uint64_t ns_[NBATCH], rh[NBATCH], ss[NBATCH];
    static uecm_pt Ps[NBATCH];
    static ecurve cvs[NBATCH];
    static mctx cs[NBATCH];
    volatile uint64_t sink = 0;
    for (int i = 0; i < NBATCH; i++) {
        uint64_t n = (rnd64() >> (64 - nbits)) | 1 | (1ULL << (nbits - 1));
        mctx_init(&cs[i], n);
        uint64_t d;
        random_ecurve(&cs[i], &cvs[i], &d);
        ns_[i] = n; rh[i] = cs[i].rho;
        Ps[i].X = to_m(&cs[i], rnd64()); Ps[i].Z = cs[i].one;
        ss[i] = to_m(&cs[i], rnd64());
    }
    printf("%d-bit n, %d curves per point; ns per curve (stage 1 only)\n", nbits, ncurves);
    printf("  B1 | PRAC x1  PRAC x2 | Edw x1   Edw x2 | x1 ratio  x2 ratio | Suyama build  Edw load\n");
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        double best[4] = { 1e30, 1e30, 1e30, 1e30 };
        for (int rep = 0; rep < 3; rep++) {
            double t0 = now();
            for (int k = 0; k < ncurves; k++) {
                int i = k & (NBATCH - 1);
                uecm_pt P = Ps[i];
                uecm_stage1(rh[i], ns_[i], &P, (uint64_t)B1, ss[i]);
                sink += P.Z;
            }
            double t1 = now();
            for (int k = 0; k < ncurves; k += 2) {
                int i = k & (NBATCH - 1);
                pt2 P; P.p[0] = Ps[i]; P.p[1] = Ps[i + 1];
                stage1_x2(&rh[i], &ns_[i], &P, B1, &ss[i]);
                sink += P.p[0].Z + P.p[1].Z;
            }
            double t2 = now();
            for (int k = 0; k < ncurves; k++) {
                int i = k & (NBATCH - 1);
                uint64_t X, Z;
                edw_stage1(&cvs[i], &ch, &X, &Z, ns_[i], rh[i]);
                sink += Z;
            }
            double t3 = now();
            for (int k = 0; k < ncurves; k += 2) {
                int i = k & (NBATCH - 1);
                uint64_t X[2], Z[2];
                edw_stage1_x2(&cvs[i], &ch, X, Z, &ns_[i], &rh[i]);
                sink += Z[0] + Z[1];
            }
            double t4 = now();
            double v[4] = { t1 - t0, t2 - t1, t3 - t2, t4 - t3 };
            for (int j = 0; j < 4; j++) if (v[j] < best[j]) best[j] = v[j];
        }
        for (int j = 0; j < 4; j++) best[j] *= 1e9 / ncurves;
        double ts2 = 1e30;
        for (int rep = 0; rep < 3; rep++) {
            double t0 = now();
            for (int k = 0; k < ncurves; k++) {
                int i = k & (NBATCH - 1);
                uecm_pt P = Ps[i];
                sink += uecm_stage2(&P, rh[i], ns_[i], (uint32_t)B1, ss[i], cs[i].one);
            }
            double t1 = now();
            if (t1 - t0 < ts2) ts2 = t1 - t0;
        }
        ts2 *= 1e9 / ncurves;
        printf("  [stage 2: %.0f ns]", ts2);
        /* curve build: Suyama (existing uecm_build: ~15M + modinv) vs loading
           a table curve (6 constants to Montgomery form, T = XY) */
        double tb = 0, tl = 0;
        if (li == 0) {
            uint64_t lcg = 1;
            int nb = ncurves * 4;
            double t0 = now();
            for (int k = 0; k < nb; k++) {
                int i = k & (NBATCH - 1);
                uecm_pt P; uint64_t s;
                uint64_t n = ns_[i];
                uint64_t five = ADD(ADD(cs[i].one, cs[i].one), ADD(ADD(cs[i].one, cs[i].one), cs[i].one));
                sink += uecm_build(&P, rh[i], ns_[i], &lcg, &s, five, cs[i].Rsqr);
                sink += s;
            }
            double t1 = now();
            for (int k = 0; k < nb; k++) {
                int i = k & (NBATCH - 1);
                uint64_t n = ns_[i], rho = rh[i];
                uint64_t X = MUL(rnd64() >> 2, cs[i].Rsqr), Y = MUL(k + 12345, cs[i].Rsqr);
                uint64_t Z = MUL(k + 777, cs[i].Rsqr), T = MUL(k + 99, cs[i].Rsqr);
                uint64_t dn = MUL(k + 5, cs[i].Rsqr), dd = MUL(k + 3, cs[i].Rsqr);
                sink += X ^ Y ^ Z ^ T ^ dn ^ dd;
            }
            double t2 = now();
            tb = (t1 - t0) * 1e9 / nb;
            tl = (t2 - t1) * 1e9 / nb;
            printf("%4d | %7.0f  %7.0f | %7.0f  %7.0f | %8.3f  %8.3f | %8.0f      %6.0f\n", B1,
                best[0], best[1], best[2], best[3], best[0] / best[2], best[1] / best[3], tb, tl);
        } else {
            printf("%4d | %7.0f  %7.0f | %7.0f  %7.0f | %8.3f  %8.3f |\n", B1,
                best[0], best[1], best[2], best[3], best[0] / best[2], best[1] / best[3]);
        }
    }
    if (sink == 42) printf(" ");
}

int main(int argc, char **argv)
{
    if (argc > 1 && strcmp(argv[1], "-t") == 0)
        return test(argc > 2 ? atoi(argv[2]) : 200) != 0;
    int nc = argc > 1 ? atoi(argv[1]) : 20000;
    int nb = argc > 2 ? atoi(argv[2]) : 64;
    bench(nc, nb);
    return 0;
}
