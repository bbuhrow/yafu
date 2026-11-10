/*
 * Phase 0 microbenchmark, tinyecm scalar (n <= 128 bits): the existing PRAC
 * stage 1 (ecm_stage1 in tinyecm.c) vs an Edwards a = -1 prototype on the
 * tiny_chains.py chains, both with the monty128 arithmetic of arith/monty.c
 * (mulmod128 / sqrmod128 / addmod128 / submod128, as tinyecm uses them).
 * Also the curve build (build_one_curve: Suyama + mpz_invert) vs a table load.
 *
 * With -DINLINE128, both PRAC (tinyecm.c's own code) and Edwards use the
 * inline __int128 CIOS arithmetic below instead of monty.c's out-of-line
 * functions (what a new scalar128 backend would do).
 *
 *   tiny_bench128 -t           Edwards chain + hand-off vs PRAC on the image
 *   tiny_bench128 [ncurves [bits]]
 */
#define _GNU_SOURCE
#ifdef INLINE128
#include "gmp.h"
#include "monty.h"
typedef unsigned __int128 u128x;

static inline __attribute__((always_inline)) void my_mulmod128(uint64_t *u, uint64_t *v,
    uint64_t *w, monty128_t *md)
{
    uint64_t n0 = md->n[0], n1 = md->n[1], rho = md->rho;
    u128x t;
    uint64_t t0, t1, t2, c, m;
    t = (u128x)u[0] * v[0]; t0 = (uint64_t)t; c = (uint64_t)(t >> 64);
    t = (u128x)u[0] * v[1] + c; t1 = (uint64_t)t; t2 = (uint64_t)(t >> 64);
    m = t0 * rho;
    t = (u128x)m * n0 + t0; c = (uint64_t)(t >> 64);
    t = (u128x)m * n1 + t1 + c; t0 = (uint64_t)t; c = (uint64_t)(t >> 64);
    t = (u128x)t2 + c; t1 = (uint64_t)t; t2 = (uint64_t)(t >> 64);
    t = (u128x)u[1] * v[0] + t0; t0 = (uint64_t)t; c = (uint64_t)(t >> 64);
    t = (u128x)u[1] * v[1] + t1 + c; t1 = (uint64_t)t; c = (uint64_t)(t >> 64);
    t2 += c;
    m = t0 * rho;
    t = (u128x)m * n0 + t0; c = (uint64_t)(t >> 64);
    t = (u128x)m * n1 + t1 + c; t0 = (uint64_t)t; c = (uint64_t)(t >> 64);
    t = (u128x)t2 + c; t1 = (uint64_t)t; t2 = (uint64_t)(t >> 64);
    u128x r = ((u128x)t1 << 64) | t0, nn = ((u128x)n1 << 64) | n0;
    u128x msk = (u128x)0 - (u128x)((t2 != 0) | (r >= nn));      /* branch-free */
    r -= nn & msk;
    w[0] = (uint64_t)r; w[1] = (uint64_t)(r >> 64);
}
static inline __attribute__((always_inline)) void my_sqrmod128(uint64_t *u, uint64_t *w,
    monty128_t *md)
{
    my_mulmod128(u, u, w, md);
}
static inline __attribute__((always_inline)) void my_addmod128(uint64_t *a, uint64_t *b,
    uint64_t *w, uint64_t *n)
{
    u128x x = ((u128x)a[1] << 64) | a[0], y = ((u128x)b[1] << 64) | b[0];
    u128x nn = ((u128x)n[1] << 64) | n[0];
    u128x s = x + y;
    u128x msk = (u128x)0 - (u128x)((s < x) | (s >= nn));
    s -= nn & msk;
    w[0] = (uint64_t)s; w[1] = (uint64_t)(s >> 64);
}
static inline __attribute__((always_inline)) void my_submod128(uint64_t *a, uint64_t *b,
    uint64_t *w, uint64_t *n)
{
    u128x x = ((u128x)a[1] << 64) | a[0], y = ((u128x)b[1] << 64) | b[0];
    u128x nn = ((u128x)n[1] << 64) | n[0];
    u128x s = x - y;
    u128x msk = (u128x)0 - (u128x)(x < y);
    s += nn & msk;
    w[0] = (uint64_t)s; w[1] = (uint64_t)(s >> 64);
}
#define mulmod128 my_mulmod128
#define sqrmod128 my_sqrmod128
#define addmod128 my_addmod128
#define submod128 my_submod128
#endif
#include "../../factor/gmp-ecm/tinyecm.c"
#include "edw_tiny_chains.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

typedef struct { uint64_t X[2], Y[2], Z[2], T[2]; } ept128;
typedef struct { uint64_t YmX[2], YpX[2], T2d[2], Z2[2]; } ecached128;

#define MUL(a, b, r) mulmod128((a), (b), (r), md)
#define SQR(a, r) sqrmod128((a), (r), md)
#define ADD(a, b, r) addmod128((a), (b), (r), md->n)
#define SUB(a, b, r) submod128((a), (b), (r), md->n)

static void edw_dbl128(ept128 *R, const ept128 *P, int want_t, monty128_t *md)
{
    uint64_t A[2], B[2], C[2], S[2], t[2], H[2], E[2], G[2], Fv[2];
    SQR(P->X, A);
    SQR(P->Y, B);
    SQR(P->Z, C);
    ADD(P->X, P->Y, t);
    SQR(t, S);
    ADD(C, C, C);
    ADD(A, B, H);
    SUB(H, S, E);
    SUB(A, B, G);
    ADD(C, G, Fv);
    MUL(E, Fv, R->X);
    MUL(G, H, R->Y);
    MUL(Fv, G, R->Z);
    if (want_t) MUL(E, H, R->T);
}

static void edw_add128(ept128 *R, const ept128 *P, const ecached128 *Q, int neg, int want_t,
    monty128_t *md)
{
    uint64_t A[2], B[2], C[2], Dz[2], t[2], E[2], H[2], Fv[2], G[2];
    SUB(P->Y, P->X, t);
    MUL(t, neg ? Q->YpX : Q->YmX, A);
    ADD(P->Y, P->X, t);
    MUL(t, neg ? Q->YmX : Q->YpX, B);
    MUL(P->T, Q->T2d, C);
    MUL(P->Z, Q->Z2, Dz);
    SUB(B, A, E);
    ADD(B, A, H);
    if (neg) { ADD(Dz, C, Fv); SUB(Dz, C, G); }
    else { SUB(Dz, C, Fv); ADD(Dz, C, G); }
    MUL(E, Fv, R->X);
    MUL(G, H, R->Y);
    MUL(Fv, G, R->Z);
    if (want_t) MUL(E, H, R->T);
}

static void edw_cache128(ecached128 *c, const ept128 *P, uint64_t *dn2, uint64_t *dd,
    uint64_t *dd2, monty128_t *md)
{
    uint64_t t[2];
    SUB(P->Y, P->X, t);
    MUL(t, dd, c->YmX);
    ADD(P->Y, P->X, t);
    MUL(t, dd, c->YpX);
    MUL(P->T, dn2, c->T2d);
    MUL(P->Z, dd2, c->Z2);
}

static void mdup_proj128(uint64_t *X, uint64_t *Z, uint64_t *a24n, uint64_t *a24d,
    monty128_t *md)
{
    uint64_t U[2], V[2], w[2], t[2], Ub[2];
    SUB(X, Z, t);
    SQR(t, U);
    ADD(X, Z, t);
    SQR(t, V);
    SUB(V, U, w);
    MUL(U, a24d, Ub);
    MUL(Ub, V, X);
    MUL(w, a24n, t);
    ADD(t, Ub, t);      /* addmod128 must not alias w with its 2nd operand */
    MUL(w, t, Z);
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
    ept128 P;
    uint64_t dn2[2], dd[2], dd2[2], a24n[2], a24d[2];
} ecurve128;

static void edw_stage1_128(ecurve128 *cv, const chain_t *c, uint64_t *Xm, uint64_t *Zm,
    monty128_t *md)
{
    ecached128 tab[16];
    ept128 ext[16];
    int k, i;
    ext[0] = cv->P;
    edw_cache128(&tab[0], &ext[0], cv->dn2, cv->dd, cv->dd2, md);
    if (c->m > 1) {
        ept128 P2;
        ecached128 c2;
        edw_dbl128(&P2, &cv->P, 1, md);
        edw_cache128(&c2, &P2, cv->dn2, cv->dd, cv->dd2, md);
        for (k = 1; k <= (c->m - 1) / 2; k++) {
            edw_add128(&ext[k], &ext[k - 1], &c2, 0, 1, md);
            edw_cache128(&tab[k], &ext[k], cv->dn2, cv->dd, cv->dd2, md);
        }
    }
    ept128 R = ext[(c->top - 1) / 2];
    for (i = 0; i < c->len; i++) {
        int dg = c->ch[i];
        if (dg == 0) {
            int t = (i + 1 < c->len) && (c->ch[i + 1] != 0);
            edw_dbl128(&R, &R, t, md);
        } else if (dg > 0) {
            edw_add128(&R, &R, &tab[(dg - 1) >> 1], 0, 0, md);
        } else {
            edw_add128(&R, &R, &tab[(-dg - 1) >> 1], 1, 0, md);
        }
    }
    ADD(R.Z, R.Y, Xm);
    SUB(R.Z, R.Y, Zm);
    for (i = 0; i < c->tz; i++)
        mdup_proj128(Xm, Zm, cv->a24n, cv->a24d, md);
}

static gmp_randstate_t grs;

static void to_m(monty128_t *md, mpz_t a, uint64_t *out)
{
    mpz_to_u128(a, out);
    to_monty128(md, out);
}

/* random a = -1 curve through a random point, d kept as dn/dd */
static void random_curve(monty128_t *md, mpz_t n, ecurve128 *cv, mpz_t d)
{
    mpz_t x, y, t, u, dd, dn;
    mpz_inits(x, y, t, u, dd, dn, NULL);
    mpz_urandomm(x, grs, n);
    mpz_urandomm(y, grs, n);
    mpz_mul(t, y, y); mpz_submul(t, x, x); mpz_sub_ui(t, t, 1);       /* y^2 - x^2 - 1 */
    mpz_mul(u, x, y); mpz_mul(u, u, u); mpz_mod(u, u, n);
    mpz_invert(u, u, n);
    mpz_mul(d, t, u); mpz_mod(d, d, n);
    mpz_urandomm(dd, grs, n);
    mpz_mul(dn, d, dd); mpz_mod(dn, dn, n);
    to_m(md, x, cv->P.X);
    to_m(md, y, cv->P.Y);
    mpz_set_ui(t, 1); to_m(md, t, cv->P.Z);
    mpz_mul(t, x, y); mpz_mod(t, t, n); to_m(md, t, cv->P.T);
    mpz_mul_2exp(t, dn, 1); mpz_mod(t, t, n); to_m(md, t, cv->dn2);
    to_m(md, dd, cv->dd);
    mpz_mul_2exp(t, dd, 1); mpz_mod(t, t, n); to_m(md, t, cv->dd2);
    to_m(md, dd, cv->a24n);
    mpz_add(t, dd, dn); mpz_mod(t, t, n); to_m(md, t, cv->a24d);
    mpz_clears(x, y, t, u, dd, dn, NULL);
}

static const int LV[7] = { 27, 47, 70, 85, 125, 165, 205 };

static int test(int ntests)
{
    int bad = 0, tot = 0;
    mpz_t n, d, t;
    mpz_inits(n, d, t, NULL);
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        for (int k = 0; k < ntests; k++) {
            mpz_urandomb(n, grs, 60 + (k % 68));
            mpz_setbit(n, 59 + (k % 68));
            mpz_nextprime(n, n);
            if (mpz_sizeinbase(n, 2) > 128) continue;
            uint64_t n128[2];
            monty128_t md;
            tinyecm_work work;
            mpz_to_u128(n, n128);
            monty128_init(&md, n128);
            tinyecm_work_init(&work);
            copy128(n128, work.n);
            work.stg1_max = B1;
            ecurve128 cv;
            random_curve(&md, n, &cv, d);
            uint64_t Xe[2], Ze[2];
            edw_stage1_128(&cv, &ch, Xe, Ze, &md);
            tinyecm_pt P;
            addmod128(cv.P.Z, cv.P.Y, P.X, md.n);
            submod128(cv.P.Z, cv.P.Y, P.Z, md.n);
            mpz_add_ui(t, d, 1); mpz_invert(t, t, n);
            to_m(&md, t, work.s);
            ecm_stage1(&md, &work, &P);
            uint64_t a[2], b[2];
            mulmod128(Xe, P.Z, a, &md);
            mulmod128(P.X, Ze, b, &md);
            int ok = a[0] == b[0] && a[1] == b[1];
            if (!ok && bad < 5) gmp_printf("mismatch B1 %d n %Zd\n", B1, n);
            bad += !ok;
            tot++;
        }
    }
    printf("test: %d of %d mismatches\n", bad, tot);
    mpz_clears(n, d, t, NULL);
    return bad;
}

#define NSET 256

static void bench(int ncurves, int nbits)
{
    static monty128_t md[NSET];
    static tinyecm_work wk[NSET];
    static ecurve128 cv[NSET];
    static tinyecm_pt P0[NSET];
    mpz_t n, d, t;
    mpz_inits(n, d, t, NULL);
    uint64_t sink = 0;
    for (int i = 0; i < NSET; i++) {
        mpz_urandomb(n, grs, nbits);
        mpz_setbit(n, nbits - 1);
        mpz_setbit(n, 0);
        uint64_t n128[2];
        mpz_to_u128(n, n128);
        monty128_init(&md[i], n128);
        tinyecm_work_init(&wk[i]);
        copy128(n128, wk[i].n);
        random_curve(&md[i], n, &cv[i], d);
        mpz_urandomm(t, grs, n); to_m(&md[i], t, wk[i].s);
        mpz_urandomm(t, grs, n); to_m(&md[i], t, P0[i].X);
        copy128(cv[i].P.Z, P0[i].Z);
    }
    printf("%d-bit n; ns per curve (stage 1)\n", nbits);
    printf("  B1 |    PRAC     Edw   ratio | build: Suyama  Edw load\n");
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        double best[2] = { 1e30, 1e30 };
        for (int i = 0; i < NSET; i++) wk[i].stg1_max = B1;
        for (int rep = 0; rep < 3; rep++) {
            double t0 = now();
            for (int k = 0; k < ncurves; k++) {
                int i = k & (NSET - 1);
                tinyecm_pt P = P0[i];
                ecm_stage1(&md[i], &wk[i], &P);
                sink += P.Z[0];
            }
            double t1 = now();
            for (int k = 0; k < ncurves; k++) {
                int i = k & (NSET - 1);
                uint64_t X[2], Z[2];
                edw_stage1_128(&cv[i], &ch, X, Z, &md[i]);
                sink += Z[0];
            }
            double t2 = now();
            if (t1 - t0 < best[0]) best[0] = t1 - t0;
            if (t2 - t1 < best[1]) best[1] = t2 - t1;
        }
        best[0] *= 1e9 / ncurves;
        best[1] *= 1e9 / ncurves;
        double ts2 = 1e30;
        for (int i = 0; i < NSET; i++) wk[i].stg2_max = 25 * B1;
        for (int rep = 0; rep < 3; rep++) {
            double t0 = now();
            for (int k = 0; k < ncurves; k++) {
                int i = k & (NSET - 1);
                tinyecm_pt P = P0[i];
                ecm_stage2(&P, &md[i], &wk[i]);
                sink += wk[i].stg2acc[0];
            }
            double t1 = now();
            if (t1 - t0 < ts2) ts2 = t1 - t0;
        }
        printf("  [stage 2: %.0f ns]", ts2 * 1e9 / ncurves);
        if (li == 0) {
            uint64_t lcg = 1;
            int nb = ncurves;
            double t0 = now();
            for (int k = 0; k < nb; k++) {
                int i = k & (NSET - 1);
                tinyecm_pt P;
                build_one_curve(&P, &md[i], &wk[i], 0, &lcg, 0);
                sink += wk[i].s[0];
            }
            double t1 = now();
            for (int k = 0; k < nb; k++) {
                int i = k & (NSET - 1);
                uint64_t c[6][2];
                for (int j = 0; j < 6; j++) {
                    c[j][0] = k * 7 + j + 3; c[j][1] = j;
                    to_monty128(&md[i], c[j]);
                }
                sink += c[0][0] ^ c[1][0] ^ c[2][0] ^ c[3][0] ^ c[4][0] ^ c[5][0];
            }
            double t2 = now();
            printf("%4d | %7.0f %7.0f %7.3f | %8.0f  %6.0f\n", B1, best[0], best[1],
                best[0] / best[1], (t1 - t0) * 1e9 / nb, (t2 - t1) * 1e9 / nb);
        } else {
            printf("%4d | %7.0f %7.0f %7.3f |\n", B1, best[0], best[1], best[0] / best[1]);
        }
    }
    if (sink == 42) printf(" ");
    mpz_clears(n, d, t, NULL);
}

int main(int argc, char **argv)
{
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, 12345);
    if (argc > 1 && strcmp(argv[1], "-t") == 0)
        return test(argc > 2 ? atoi(argv[2]) : 100) != 0;
    bench(argc > 1 ? atoi(argv[1]) : 20000, argc > 2 ? atoi(argv[2]) : 100);
    return 0;
}
