/*
 * AVX-ECM stage 1 on a = -1 twisted Edwards curves -x^2 + y^2 = 1 + d x^2 y^2
 * in extended coordinates (Hisil, Wong, Carter, Dawson 2008), with signed
 * windows over chunks of the stage-1 scalar.  Selected with -avxecm_curve 2;
 * the curves are Gelin-Kleinjung-Lenstra family (iv), built in avxecm.c
 * (vec_build_one_curve_edwards).  Design, bounds and results:
 * arith/avxecm_edwards_plan.md; reference: arith/edwards/edw_ref.py.
 *
 * Every field operation goes through the backend pointers (vecmulmod_ptr,
 * vecsqrmod_ptr, vecaddmod_ptr, vecsubmod_ptr), so this runs on the template
 * cores (IFMA, AVX-512F, AVX2) and on the RNS backend alike.  All 8 lanes
 * (curves) share the scalar, hence the digit string and the control flow.
 *
 * Operations (M = mul, S = sqr, a = add/sub):
 *   edw_dbl    3M + 4S + 6a, T only on request (+1M)
 *   edw_add    P + cached Q: 7M + 6a, T only on request (+1M)
 *   edw_cache  (Y-X, Y+X, 2d T, 2Z): 1M + 3a
 * A doubling computes T only when an addition follows; an addition is always
 * followed by a doubling and never computes T, except at the end of a chunk
 * (the next chunk's table needs T).
 *
 * At the end of stage 1 the point is handed to stage 2 as its Montgomery
 * image, X_M = Z + Y, Z_M = Z - Y (u = (1+y)/(1-y); work->s already holds
 * (A+2)/4 = 1/(1+d) from the curve build).  The stage-1 gcd uses the
 * Edwards X, which vanishes mod p at the identity (0 : 1 : 1) and at the
 * 2-torsion point (0 : -1 : 1).
 */

#include "avx_ecm.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EDW_W 6                     // signed-window width (16 table points)
#define EDW_CHUNK_BITS 16384        // scalar chunk: product of prime powers

static void edw_pt_init(edw_pt *P, uint32_t words)
{
    P->X = vecInit(words);
    P->Y = vecInit(words);
    P->Z = vecInit(words);
    P->T = vecInit(words);
}

static void edw_pt_free(edw_pt *P)
{
    vecFree(P->X);
    vecFree(P->Y);
    vecFree(P->Z);
    vecFree(P->T);
}

static void edw_cpt_init(edw_cpt *C, uint32_t words)
{
    C->ymx = vecInit(words);
    C->ypx = vecInit(words);
    C->t2d = vecInit(words);
    C->z2 = vecInit(words);
}

static void edw_cpt_free(edw_cpt *C)
{
    vecFree(C->ymx);
    vecFree(C->ypx);
    vecFree(C->t2d);
    vecFree(C->z2);
}

void avxecm_edw_alloc(ecm_work *work, uint32_t words)
{
    avxecm_edw_t *e = (avxecm_edw_t *)calloc(1, sizeof(avxecm_edw_t));
    int i;

    e->w = EDW_W;
    e->ntab = 1 << (EDW_W - 2);
    edw_pt_init(&e->R, words);
    edw_pt_init(&e->P2, words);
    edw_cpt_init(&e->c2, words);
    e->d2 = vecInit(words);
    e->tab = (edw_cpt *)malloc((size_t)e->ntab * sizeof(edw_cpt));
    for (i = 0; i < e->ntab; i++)
        edw_cpt_init(&e->tab[i], words);
    for (i = 0; i < AVXECM_EDW_NTMP; i++)
        e->t[i] = vecInit(words);
    e->diglen = EDW_CHUNK_BITS + 64;
    e->dig = (int8_t *)malloc((size_t)e->diglen);
    work->edw = e;
}

void avxecm_edw_free(ecm_work *work)
{
    avxecm_edw_t *e = work->edw;
    int i;

    if (e == NULL)
        return;
    edw_pt_free(&e->R);
    edw_pt_free(&e->P2);
    edw_cpt_free(&e->c2);
    vecFree(e->d2);
    for (i = 0; i < e->ntab; i++)
        edw_cpt_free(&e->tab[i]);
    free(e->tab);
    for (i = 0; i < AVXECM_EDW_NTMP; i++)
        vecFree(e->t[i]);
    free(e->dig);
    free(e);
    work->edw = NULL;
}

#define MUL(a, b, c) vecmulmod_ptr(a, b, c, work->n, work->tt5, mdata)
#define SQR(a, c) vecsqrmod_ptr(a, c, work->n, work->tt5, mdata)
#define ADD(a, b, c) vecaddmod_ptr(a, b, c, mdata)
#define SUB(a, b, c) vecsubmod_ptr(a, b, c, mdata)

// Q = 2P (Q may be P).  dbl-2008-hwcd with a = -1 and the signs folded:
// E' = -E, F' = -F, G' = -G, H' = -H leave X3 = EF, Y3 = GH, Z3 = FG,
// T3 = EH unchanged and need no negation:
//   A = X^2, B = Y^2, C = 2Z^2, H' = A + B, E' = H' - (X+Y)^2,
//   G' = A - B, F' = C + G'
// Every subtrahend is a product output (the radix add/sub bound, plan §D).
static void edw_dbl(vec_monty_t *mdata, ecm_work *work, edw_pt *P, edw_pt *Q, int want_t)
{
    avxecm_edw_t *e = work->edw;
    vec_bignum_t *A = e->t[0], *B = e->t[1], *C = e->t[2], *S = e->t[3];
    vec_bignum_t *E = e->t[4], *F = e->t[5], *G = e->t[6], *H = e->t[7];
    vec_bignum_t *u = e->t[8], *zz = e->t[9];

    SQR(P->X, A);
    SQR(P->Y, B);
    SQR(P->Z, zz);
    ADD(zz, zz, C);
    ADD(P->X, P->Y, u);
    SQR(u, S);
    ADD(A, B, H);
    SUB(H, S, E);
    SUB(A, B, G);
    ADD(C, G, F);
    MUL(E, F, Q->X);
    MUL(G, H, Q->Y);
    MUL(F, G, Q->Z);
    if (want_t)
        MUL(E, H, Q->T);
    e->ndbl++;
}

// R = P + s*C, C cached, s = +1 or -1 (R may be P).  add-2008-hwcd-3
// (unified, a = -1):
//   A = (Y1-X1)(Y2-X2), B = (Y1+X1)(Y2+X2), C = T1 2d T2, D = Z1 2 Z2,
//   E = B - A, F = D - C, G = D + C, H = B + A
// -C = (-X2, Y2, Z2, -T2): Y2-X2 and Y2+X2 swap, and 2d T2 changes sign,
// which swaps F and G.
static void edw_add(vec_monty_t *mdata, ecm_work *work, edw_pt *P, edw_cpt *Cq,
    int neg, edw_pt *R, int want_t)
{
    avxecm_edw_t *e = work->edw;
    vec_bignum_t *A = e->t[0], *B = e->t[1], *C = e->t[2], *D = e->t[3];
    vec_bignum_t *E = e->t[4], *F = e->t[5], *G = e->t[6], *H = e->t[7];
    vec_bignum_t *u = e->t[8], *v = e->t[9];

    SUB(P->Y, P->X, u);
    ADD(P->Y, P->X, v);
    MUL(u, neg ? Cq->ypx : Cq->ymx, A);
    MUL(v, neg ? Cq->ymx : Cq->ypx, B);
    MUL(P->T, Cq->t2d, C);
    MUL(P->Z, Cq->z2, D);
    SUB(B, A, E);
    ADD(B, A, H);
    if (neg)
    {
        ADD(D, C, F);
        SUB(D, C, G);
    }
    else
    {
        SUB(D, C, F);
        ADD(D, C, G);
    }
    MUL(E, F, R->X);
    MUL(G, H, R->Y);
    MUL(F, G, R->Z);
    if (want_t)
        MUL(E, H, R->T);
    e->nadd++;
}

// C = cached form of P (with T)
static void edw_cache(vec_monty_t *mdata, ecm_work *work, edw_pt *P, edw_cpt *C)
{
    avxecm_edw_t *e = work->edw;
    SUB(P->Y, P->X, C->ymx);
    ADD(P->Y, P->X, C->ypx);
    MUL(P->T, e->d2, C->t2d);
    ADD(P->Z, P->Z, C->z2);
}

// width-w NAF of k (as secp256k1_ecmult_wnaf): dig[i] in {0, +-1, +-3, ...,
// +-(2^(w-1) - 1)}, k = sum dig[i] 2^i.  Returns the index of the top
// (positive) digit.
static int edw_wnaf(int8_t *dig, int len, mpz_t k, int w)
{
    int bits = (int)mpz_sizeinbase(k, 2);
    int bit = 0, carry = 0, top = 0, n = bits + 1;

    if (n > len)
    {
        printf("avxecm_edw: scalar chunk too long (%d bits)\n", bits);
        exit(1);
    }
    memset(dig, 0, (size_t)n);
    while (bit < n)
    {
        int now, word, j;
        if ((int)mpz_tstbit(k, (mp_bitcnt_t)bit) == carry)
        {
            bit++;
            continue;
        }
        now = w;
        if (now > n - bit)
            now = n - bit;
        for (word = 0, j = now - 1; j >= 0; j--)
            word = (word << 1) | (int)mpz_tstbit(k, (mp_bitcnt_t)(bit + j));
        word += carry;
        carry = (word >> (w - 1)) & 1;
        word -= carry << w;
        dig[bit] = (int8_t)word;
        top = bit;
        bit += now;
    }
    return top;
}

// R = [k]R for one chunk k > 0.  R has T on entry and on exit.
static void edw_mul_chunk(vec_monty_t *mdata, ecm_work *work, mpz_t k)
{
    avxecm_edw_t *e = work->edw;
    edw_pt *R = &e->R;
    int i, top;

    top = edw_wnaf(e->dig, e->diglen, k, e->w);

    // table of cached odd multiples R, 3R, ..., (2 ntab - 1) R; R is the
    // running point while building it
    edw_cache(mdata, work, R, &e->tab[0]);
    edw_dbl(mdata, work, R, &e->P2, 1);
    edw_cache(mdata, work, &e->P2, &e->c2);
    for (i = 1; i < e->ntab; i++)
    {
        edw_add(mdata, work, R, &e->c2, 0, R, 1);
        edw_cache(mdata, work, R, &e->tab[i]);
    }

    // R = identity (0 : 1 : 1 : 0); the unified addition handles it, so the
    // top digit is an ordinary step
    vecClear(R->X);
    vecClear(R->T);
    vecCopy(mdata->one, R->Y);
    vecCopy(mdata->one, R->Z);
    for (i = top; i >= 0; i--)
    {
        int v = e->dig[i];
        if (i < top)
            edw_dbl(mdata, work, R, R, (v != 0) || (i == 0));
        if (v != 0)
            edw_add(mdata, work, R, &e->tab[(v < 0 ? -v : v) >> 1], v < 0, R, i == 0);
    }
}

// Stage 1 over the primes of one range, like vec_ecm_stage1 (avxecm.c):
// 2^k < b1, then q^j < b1 for primes[1..] < b1.  The prime powers are
// multiplied into chunks of about EDW_CHUNK_BITS bits, each a signed-window
// scalar multiplication of the Edwards point work->edw->R.  On return P holds
// the Montgomery image (Z + Y : Z - Y) for stage 2.
void avxecm_edw_stage1(vec_monty_t *mdata, ecm_work *work, ecm_pt *P,
    uint64_t b1, uint64_t *primes, uint64_t nump, int verbose)
{
    avxecm_edw_t *e = work->edw;
    mpz_t k;
    uint64_t q, i;

    mpz_init_set_ui(k, 1);
    e->ndbl = e->nadd = 0;

    q = 2;
    while (q < b1)
    {
        mpz_mul_2exp(k, k, 1);
        q *= 2;
    }

    for (i = 1; (i < nump) && (primes[i] < b1); i++)
    {
        uint64_t c = 1;
        q = primes[i];
        do {
            mpz_mul_ui(k, k, (unsigned long)q);
            c *= q;
        } while ((c * q) < b1);

        if (mpz_sizeinbase(k, 2) >= EDW_CHUNK_BITS - 64)
        {
            edw_mul_chunk(mdata, work, k);
            mpz_set_ui(k, 1);
        }

        if ((verbose > 1) && ((i & 8191) == 0))
        {
            printf("accumulating prime %" PRIu64 "\r", q);
            fflush(stdout);
        }
    }
    if (mpz_cmp_ui(k, 1) > 0)
        edw_mul_chunk(mdata, work, k);

    work->last_pid = (uint32_t)i;

    // Montgomery image for stage 2 (and the save_b1 / intermediate output)
    vecaddmod_ptr(e->R.Z, e->R.Y, P->X, mdata);
    vecsubmod_ptr(e->R.Z, e->R.Y, P->Z, mdata);

    if (verbose > 1)
    {
        printf("\nStage 1 (Edwards, w = %d) completed at prime %" PRIu64 " with %u point-adds "
            "and %u point-doubles\n", e->w, primes[i - 1], e->nadd, e->ndbl);
        fflush(stdout);
    }

    mpz_clear(k);
}
