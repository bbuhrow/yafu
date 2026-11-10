/*
 * Phase 0 microbenchmark, tinyecm x8 (AVX-512, n <= 104 bits, 2 x 52-bit
 * digits): the existing PRAC stage 1 (ecm_stage1_x8) vs an Edwards a = -1
 * prototype on the tiny_chains.py chains, with tinyecm.c's own 104-bit
 * vector arithmetic.  Also stage 2 and the curve build (build_curves_104_x8:
 * Suyama + 8 mpz_invert).
 *
 * Built against a copy of tinyecm.c with -DIFMA104_SRC (see tiny_ifma104.py),
 * tecm_mulredc104_x8 does all its 52x52 products with madd52lo/hi instead
 * of the FP fmadd_round pairs (IFMA only); PRAC, stage 2 and Edwards all
 * use it then.
 *
 *   tiny_bench104 -t           Edwards chain + hand-off vs PRAC on the image
 *   tiny_bench104 [nvec [bits]]
 */
#define _GNU_SOURCE
#ifdef IFMA104_SRC
#include "tinyecm_ifma104.c"
#else
#include "../../factor/gmp-ecm/tinyecm.c"
#endif
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

typedef vec_u104_t W;
typedef struct { ALIGNED_MEM W X, Y, Z, T; } ept104;
typedef struct { ALIGNED_MEM W YmX, YpX, T2d, Z2; } ecached104;

#define MUL(a, b, r) tecm_mulredc104_x8((r), (a), (b), n, rho)
#define SQR(a, r) tecm_sqrredc104_x8((r), (a), n, rho)
#define ADDM(a, b, r) tecm_addmod104_x8((r), (a), (b), n)
#define SUBM(a, b, r) tecm_submod104_x8((r), (a), (b), n)

static void edw_dbl104(ept104 *R, ept104 *P, int want_t, W *n, W *rho)
{
    ALIGNED_MEM W A, B, C, S, t, H, E, G, Fv;
    SQR(&P->X, &A);
    SQR(&P->Y, &B);
    SQR(&P->Z, &C);
    ADDM(&P->X, &P->Y, &t);
    SQR(&t, &S);
    ADDM(&C, &C, &C);
    ADDM(&A, &B, &H);
    SUBM(&H, &S, &E);
    SUBM(&A, &B, &G);
    ADDM(&C, &G, &Fv);
    MUL(&E, &Fv, &R->X);
    MUL(&G, &H, &R->Y);
    MUL(&Fv, &G, &R->Z);
    if (want_t) MUL(&E, &H, &R->T);
}

static void edw_add104(ept104 *R, ept104 *P, ecached104 *Q, int neg, int want_t, W *n, W *rho)
{
    ALIGNED_MEM W A, B, C, Dz, t, E, H, Fv, G;
    SUBM(&P->Y, &P->X, &t);
    MUL(&t, neg ? &Q->YpX : &Q->YmX, &A);
    ADDM(&P->Y, &P->X, &t);
    MUL(&t, neg ? &Q->YmX : &Q->YpX, &B);
    MUL(&P->T, &Q->T2d, &C);
    MUL(&P->Z, &Q->Z2, &Dz);
    SUBM(&B, &A, &E);
    ADDM(&B, &A, &H);
    if (neg) { ADDM(&Dz, &C, &Fv); SUBM(&Dz, &C, &G); }
    else { SUBM(&Dz, &C, &Fv); ADDM(&Dz, &C, &G); }
    MUL(&E, &Fv, &R->X);
    MUL(&G, &H, &R->Y);
    MUL(&Fv, &G, &R->Z);
    if (want_t) MUL(&E, &H, &R->T);
}

static void edw_cache104(ecached104 *c, ept104 *P, W *dn2, W *dd, W *dd2, W *n, W *rho)
{
    ALIGNED_MEM W t;
    SUBM(&P->Y, &P->X, &t);
    MUL(&t, dd, &c->YmX);
    ADDM(&P->Y, &P->X, &t);
    MUL(&t, dd, &c->YpX);
    MUL(&P->T, dn2, &c->T2d);
    MUL(&P->Z, dd2, &c->Z2);
}

static void mdup_proj104(W *X, W *Z, W *a24n, W *a24d, W *n, W *rho)
{
    ALIGNED_MEM W U, V, w, t, Ub;
    SUBM(X, Z, &t);
    SQR(&t, &U);
    ADDM(X, Z, &t);
    SQR(&t, &V);
    SUBM(&V, &U, &w);
    MUL(&U, a24d, &Ub);
    MUL(&Ub, &V, X);
    MUL(&w, a24n, &t);
    ADDM(&t, &Ub, &t);
    MUL(&w, &t, Z);
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
    ALIGNED_MEM ept104 P;
    ALIGNED_MEM W dn2, dd, dd2, a24n, a24d;
} ecurve104;

static void edw_stage1_104(ecurve104 *cv, const chain_t *c, W *Xm, W *Zm, W *n, W *rho)
{
    ALIGNED_MEM ecached104 tab[16];
    ALIGNED_MEM ept104 ext[16];
    ALIGNED_MEM ept104 R;
    int k, i;
    ext[0] = cv->P;
    edw_cache104(&tab[0], &ext[0], &cv->dn2, &cv->dd, &cv->dd2, n, rho);
    if (c->m > 1) {
        ALIGNED_MEM ept104 P2;
        ALIGNED_MEM ecached104 c2;
        edw_dbl104(&P2, &cv->P, 1, n, rho);
        edw_cache104(&c2, &P2, &cv->dn2, &cv->dd, &cv->dd2, n, rho);
        for (k = 1; k <= (c->m - 1) / 2; k++) {
            edw_add104(&ext[k], &ext[k - 1], &c2, 0, 1, n, rho);
            edw_cache104(&tab[k], &ext[k], &cv->dn2, &cv->dd, &cv->dd2, n, rho);
        }
    }
    R = ext[(c->top - 1) / 2];
    for (i = 0; i < c->len; i++) {
        int dg = c->ch[i];
        if (dg == 0) {
            int t = (i + 1 < c->len) && (c->ch[i + 1] != 0);
            edw_dbl104(&R, &R, t, n, rho);
        } else if (dg > 0) {
            edw_add104(&R, &R, &tab[(dg - 1) >> 1], 0, 0, n, rho);
        } else {
            edw_add104(&R, &R, &tab[(-dg - 1) >> 1], 1, 0, n, rho);
        }
    }
    ADDM(&R.Z, &R.Y, Xm);
    SUBM(&R.Z, &R.Y, Zm);
    for (i = 0; i < c->tz; i++)
        mdup_proj104(Xm, Zm, &cv->a24n, &cv->a24d, n, rho);
}

static gmp_randstate_t grs;

/* lane i of w <- a * 2^104 mod n */
static void set_m(W *w, int i, mpz_t a, mpz_t n)
{
    mpz_t t;
    mpz_init(t);
    mpz_mul_2exp(t, a, 104);
    mpz_mod(t, t, n);
    mpz_to_u104(t, w, i);
    mpz_clear(t);
}

/* per lane: modulus, Montgomery constants, a random a = -1 curve through a
   random point (d = dn/dd kept projective), and 1/(1+d) for the PRAC check */
static void setup(monty104_x8_t *md, ecurve104 *cv, W *a24, W *rnd, int nbits, int prime)
{
    mpz_t n, R, x, y, t, u, d, dd, dn;
    mpz_inits(n, R, x, y, t, u, d, dd, dn, NULL);
    mpz_set_ui(R, 1);
    mpz_mul_2exp(R, R, 104);
    for (int i = 0; i < 8; i++) {
        do {
            mpz_urandomb(n, grs, nbits);
            mpz_setbit(n, nbits - 1);
            mpz_setbit(n, 0);
            if (prime) mpz_nextprime(n, n);
        } while (mpz_sizeinbase(n, 2) > 104 || mpz_sizeinbase(n, 2) < 20);
        mpz_to_u104(n, &md->n, i);
        mpz_mod(t, R, n);
        mpz_to_u104(t, &md->one, i);
        mpz_invert(t, n, R);
        mpz_to_u104(t, &md->rho, i);
        mpz_urandomm(x, grs, n);
        mpz_urandomm(y, grs, n);
        mpz_mul(t, y, y); mpz_submul(t, x, x); mpz_sub_ui(t, t, 1);
        mpz_mul(u, x, y); mpz_mul(u, u, u); mpz_mod(u, u, n);
        if (!mpz_invert(u, u, n)) mpz_set_ui(u, 1);
        mpz_mul(d, t, u); mpz_mod(d, d, n);
        mpz_urandomm(dd, grs, n);
        mpz_mul(dn, d, dd); mpz_mod(dn, dn, n);
        set_m(&cv->P.X, i, x, n);
        set_m(&cv->P.Y, i, y, n);
        mpz_set_ui(t, 1); set_m(&cv->P.Z, i, t, n);
        mpz_mul(t, x, y); mpz_mod(t, t, n); set_m(&cv->P.T, i, t, n);
        mpz_mul_2exp(t, dn, 1); mpz_mod(t, t, n); set_m(&cv->dn2, i, t, n);
        set_m(&cv->dd, i, dd, n);
        mpz_mul_2exp(t, dd, 1); mpz_mod(t, t, n); set_m(&cv->dd2, i, t, n);
        set_m(&cv->a24n, i, dd, n);
        mpz_add(t, dd, dn); mpz_mod(t, t, n); set_m(&cv->a24d, i, t, n);
        mpz_add_ui(t, d, 1);
        if (!mpz_invert(t, t, n)) mpz_set_ui(t, 1);
        set_m(a24, i, t, n);
        mpz_urandomm(t, grs, n);
        set_m(rnd, i, t, n);
    }
    mpz_clears(n, R, x, y, t, u, d, dd, dn, NULL);
}

static void init_consts(void)
{
    dbias = _mm512_castsi512_pd(tecm_set64(0x4670000000000000ULL));
    vbias1 = tecm_set64(0x4670000000000000ULL);
    vbias2 = tecm_set64(0x4670000000000001ULL);
    lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
}

static const int LV[7] = { 27, 47, 70, 85, 125, 165, 205 };

static int test(int ntests)
{
    int bad = 0, tot = 0;
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        for (int k = 0; k < ntests; k++) {
            ALIGNED_MEM monty104_x8_t md;
            ALIGNED_MEM ecurve104 cv;
            ALIGNED_MEM tinyecm_work_x8 work;
            ALIGNED_MEM W a24, rnd, Xe, Ze, a, b;
            ALIGNED_MEM tecm_pt_x8 P;
            setup(&md, &cv, &a24, &rnd, 30 + (k % 75), 1);
            edw_stage1_104(&cv, &ch, &Xe, &Ze, &md.n, &md.rho);
            tecm_addmod104_x8(&P.X, &cv.P.Z, &cv.P.Y, &md.n);
            tecm_submod104_x8(&P.Z, &cv.P.Z, &cv.P.Y, &md.n);
            work.s = a24;
            work.n = md.n;
            work.stg1_max = B1;
            ecm_stage1_x8(&md, &work, &P);
            tecm_mulredc104_x8(&a, &Xe, &P.Z, &md.n, &md.rho);
            tecm_mulredc104_x8(&b, &P.X, &Ze, &md.n, &md.rho);
            for (int i = 0; i < 8; i++) {
                int ok = a.data[0][i] == b.data[0][i] && a.data[1][i] == b.data[1][i];
                if (!ok && bad < 5) printf("mismatch B1 %d lane %d\n", B1, i);
                bad += !ok;
                tot++;
            }
        }
    }
    printf("test: %d of %d lanes mismatch\n", bad, tot);
    return bad;
}

#define NSET 32

static void bench(int nvec, int nbits)
{
    static ALIGNED_MEM monty104_x8_t md[NSET];
    static ALIGNED_MEM ecurve104 cv[NSET];
    static ALIGNED_MEM tinyecm_work_x8 wk[NSET];
    static ALIGNED_MEM tecm_pt_x8 P0[NSET];
    ALIGNED_MEM W a24, sinkv;
    uint64_t sink = 0;
    for (int k = 0; k < NSET; k++) {
        setup(&md[k], &cv[k], &a24, &P0[k].X, nbits, 0);
        wk[k].s = a24;
        wk[k].n = md[k].n;
        P0[k].Z = cv[k].P.Z;
    }
    printf("%d-bit n; ns per vector of 8 curves\n", nbits);
    printf("  B1 |    PRAC     Edw   ratio | stage 2 | build: Suyama x8\n");
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        double best[3] = { 1e30, 1e30, 1e30 };
        for (int k = 0; k < NSET; k++) { wk[k].stg1_max = B1; wk[k].stg2_max = 25 * B1; }
        for (int rep = 0; rep < 3; rep++) {
            double t0 = now();
            for (int k = 0; k < nvec; k++) {
                int i = k & (NSET - 1);
                ALIGNED_MEM tecm_pt_x8 P = P0[i];
                ecm_stage1_x8(&md[i], &wk[i], &P);
                sink += P.Z.data[0][0];
            }
            double t1 = now();
            for (int k = 0; k < nvec; k++) {
                int i = k & (NSET - 1);
                ALIGNED_MEM W X, Z;
                edw_stage1_104(&cv[i], &ch, &X, &Z, &md[i].n, &md[i].rho);
                sink += Z.data[0][0];
            }
            double t2 = now();
            for (int k = 0; k < nvec; k++) {
                int i = k & (NSET - 1);
                ALIGNED_MEM tecm_pt_x8 P = P0[i];
                ecm_stage2_x8(&P, &md[i], &wk[i]);
                sink += wk[i].stg2acc.data[0][0];
            }
            double t3 = now();
            double v[3] = { t1 - t0, t2 - t1, t3 - t2 };
            for (int j = 0; j < 3; j++) if (v[j] < best[j]) best[j] = v[j];
        }
        for (int j = 0; j < 3; j++) best[j] *= 1e9 / nvec;
        if (li == 0) {
            uint64_t lcg = 1;
            int nb = nvec;
            double t0 = now();
            for (int k = 0; k < nb; k++) {
                int i = k & (NSET - 1);
                ALIGNED_MEM tecm_pt_x8 P;
                build_curves_104_x8(&P, &md[i], &wk[i], &lcg, 0);
                sink += wk[i].s.data[0][0];
            }
            double t1 = now();
            for (int k = 0; k < NSET; k++) wk[k].s = cv[k].a24n;   /* restore a valid s */
            printf("%4d | %7.0f %7.0f %7.3f | %7.0f | %8.0f\n", B1, best[0], best[1],
                best[0] / best[1], best[2], (t1 - t0) * 1e9 / nb);
        } else {
            printf("%4d | %7.0f %7.0f %7.3f | %7.0f |\n", B1, best[0], best[1],
                best[0] / best[1], best[2]);
        }
    }
    if (sink == 42) printf(" ");
    (void)sinkv;
}

int main(int argc, char **argv)
{
    gmp_randinit_default(grs);
    gmp_randseed_ui(grs, 4321);
    init_consts();
    if (argc > 1 && strcmp(argv[1], "-t") == 0)
        return test(argc > 2 ? atoi(argv[2]) : 40) != 0;
    bench(argc > 1 ? atoi(argv[1]) : 5000, argc > 2 ? atoi(argv[2]) : 100);
    return 0;
}
