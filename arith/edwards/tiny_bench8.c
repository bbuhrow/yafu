/*
 * Phase 0 microbenchmark, AVX-512 x8 microecm (n < 2^52): the existing PRAC
 * stage 1 (uecm_stage1_x8) vs an Edwards a = -1 prototype on the
 * tiny_chains.py chains, with microecm.c's own x8 arithmetic.  Also the
 * curve build (uecm_build_x8: Suyama + 8 scalar modinv) vs a table load.
 *
 * Build with -DUSE_AVX512F (and -DIFMA for the IFMA products).  Built
 * against a copy of microecm.c with -DM52LO_SRC (see tiny_bench8.sh), m is
 * computed with madd52lo instead of vpmullq (IFMA only).
 *
 *   tiny_bench8 -t      Edwards chain + hand-off vs PRAC on the image, per lane
 *   tiny_bench8 [nvec]  timings, ns per 8-curve vector
 */
#define _GNU_SOURCE
#ifdef M52LO_SRC
#include "microecm_m52lo.c"
#else
#include "../../factor/gmp-ecm/microecm.c"
#endif
#include "edw_tiny_chains.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

typedef __m512i V;

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

typedef struct { V X, Y, Z, T; } ept8;
typedef struct { V YmX, YpX, T2d, Z2; } ecached8;

#define MUL(a, b) uecm_mulredc_x8((a), (b), n, rho)
#define SQR(a) uecm_sqrredc_x8((a), n, rho)
#define ADD(a, b) uecm_addmod_x8((a), (b), n)
#define SUB(a, b) uecm_submod_x8((a), (b), n)

MICRO_ECM_FORCE_INLINE static void edw_dbl8(ept8 *R, const ept8 *P, int want_t, V n, V rho)
{
    V A = SQR(P->X);
    V B = SQR(P->Y);
    V C = SQR(P->Z);
    V S = SQR(ADD(P->X, P->Y));
    C = ADD(C, C);
    V H = ADD(A, B);
    V E = SUB(H, S);
    V G = SUB(A, B);
    V Fv = ADD(C, G);
    R->X = MUL(E, Fv);
    R->Y = MUL(G, H);
    R->Z = MUL(Fv, G);
    if (want_t) R->T = MUL(E, H);
}

MICRO_ECM_FORCE_INLINE static void edw_add8(ept8 *R, const ept8 *P, const ecached8 *Q,
    int neg, int want_t, V n, V rho)
{
    V A = MUL(SUB(P->Y, P->X), neg ? Q->YpX : Q->YmX);
    V B = MUL(ADD(P->Y, P->X), neg ? Q->YmX : Q->YpX);
    V C = MUL(P->T, Q->T2d);
    V D = MUL(P->Z, Q->Z2);
    V E = SUB(B, A);
    V H = ADD(B, A);
    V Fv = neg ? ADD(D, C) : SUB(D, C);
    V G = neg ? SUB(D, C) : ADD(D, C);
    R->X = MUL(E, Fv);
    R->Y = MUL(G, H);
    R->Z = MUL(Fv, G);
    if (want_t) R->T = MUL(E, H);
}

MICRO_ECM_FORCE_INLINE static void edw_cache8(ecached8 *c, const ept8 *P, V dn2, V dd, V dd2,
    V n, V rho)
{
    c->YmX = MUL(SUB(P->Y, P->X), dd);
    c->YpX = MUL(ADD(P->Y, P->X), dd);
    c->T2d = MUL(P->T, dn2);
    c->Z2 = MUL(P->Z, dd2);
}

MICRO_ECM_FORCE_INLINE static void mdup_proj8(V *X, V *Z, V a24n, V a24d, V n, V rho)
{
    V U = SQR(SUB(*X, *Z));
    V Vv = SQR(ADD(*X, *Z));
    V w = SUB(Vv, U);
    V Ub = MUL(U, a24d);
    *X = MUL(Ub, Vv);
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
    ept8 P;
    V dn2, dd, dd2, a24n, a24d;
} ecurve8;

static void edw_stage1_8(const ecurve8 *cv, const chain_t *c, V *Xm, V *Zm, V n, V rho)
{
    ecached8 tab[16];
    ept8 ext[16];
    int k, i;
    ext[0] = cv->P;
    edw_cache8(&tab[0], &ext[0], cv->dn2, cv->dd, cv->dd2, n, rho);
    if (c->m > 1) {
        ept8 P2;
        ecached8 c2;
        edw_dbl8(&P2, &cv->P, 1, n, rho);
        edw_cache8(&c2, &P2, cv->dn2, cv->dd, cv->dd2, n, rho);
        for (k = 1; k <= (c->m - 1) / 2; k++) {
            edw_add8(&ext[k], &ext[k - 1], &c2, 0, 1, n, rho);
            edw_cache8(&tab[k], &ext[k], cv->dn2, cv->dd, cv->dd2, n, rho);
        }
    }
    ept8 R = ext[(c->top - 1) / 2];
    for (i = 0; i < c->len; i++) {
        int dg = c->ch[i];
        if (dg == 0) {
            int t = (i + 1 < c->len) && (c->ch[i + 1] != 0);
            edw_dbl8(&R, &R, t, n, rho);
        } else if (dg > 0) {
            edw_add8(&R, &R, &tab[(dg - 1) >> 1], 0, 0, n, rho);
        } else {
            edw_add8(&R, &R, &tab[(-dg - 1) >> 1], 1, 0, n, rho);
        }
    }
    V X = ADD(R.Z, R.Y), Z = SUB(R.Z, R.Y);
    for (i = 0; i < c->tz; i++)
        mdup_proj8(&X, &Z, cv->a24n, cv->a24d, n, rho);
    *Xm = X;
    *Zm = Z;
}

/* ---- per-lane scalar helpers (R = 2^52) ---- */
static uint64_t mulmod_plain(uint64_t a, uint64_t b, uint64_t n)
{
    return (uint64_t)(((unsigned __int128)a * b) % n);
}
static uint64_t to_m52(uint64_t a, uint64_t n) { return mulmod_plain(a % n, (1ULL << 52) % n, n); }
static uint64_t inv_plain(uint64_t a, uint64_t n)
{
    uint64_t g;
    return uecm_modinv_64(a % n, n, &g);
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

typedef struct {
    uint64_t n[8], rho[8], one[8], Rsqr[8], five[8];
    uint64_t x[8], y[8], d[8], dd[8];
} lanes_t;

static void setup_lanes(lanes_t *L, int nbits, int prime, ecurve8 *cv, V *vn, V *vrho)
{
    ALIGNED_MEM uint64_t a[11][8];
    for (int i = 0; i < 8; i++) {
        uint64_t n;
        do {
            n = (rnd64() >> (64 - nbits)) | 1 | (1ULL << (nbits - 1));
        } while (prime && !is_prime64(n));
        L->n[i] = n;
        L->rho[i] = uecm_multiplicative_inverse52(n);
        L->one[i] = (1ULL << 52) % n;
        L->Rsqr[i] = mulmod_plain(L->one[i], L->one[i], n);
        L->five[i] = to_m52(5, n);
        uint64_t x = rnd64() % n, y = rnd64() % n;
        uint64_t x2 = mulmod_plain(x, x, n), y2 = mulmod_plain(y, y, n);
        uint64_t num = (y2 + 2 * n - x2 - 1) % n;
        uint64_t d = mulmod_plain(num, inv_plain(mulmod_plain(x2, y2, n), n), n);
        uint64_t dd = rnd64() % n | 1;
        uint64_t dn = mulmod_plain(d, dd, n);
        L->x[i] = x; L->y[i] = y; L->d[i] = d; L->dd[i] = dd;
        a[0][i] = to_m52(x, n);
        a[1][i] = to_m52(y, n);
        a[2][i] = L->one[i];
        a[3][i] = to_m52(mulmod_plain(x, y, n), n);
        a[4][i] = to_m52(2 * dn % n, n);
        a[5][i] = to_m52(dd, n);
        a[6][i] = to_m52(2 * dd % n, n);
        a[7][i] = to_m52(dd, n);
        a[8][i] = to_m52((dd + dn) % n, n);
    }
    cv->P.X = _mm512_loadu_si512(a[0]);
    cv->P.Y = _mm512_loadu_si512(a[1]);
    cv->P.Z = _mm512_loadu_si512(a[2]);
    cv->P.T = _mm512_loadu_si512(a[3]);
    cv->dn2 = _mm512_loadu_si512(a[4]);
    cv->dd = _mm512_loadu_si512(a[5]);
    cv->dd2 = _mm512_loadu_si512(a[6]);
    cv->a24n = _mm512_loadu_si512(a[7]);
    cv->a24d = _mm512_loadu_si512(a[8]);
    *vn = _mm512_loadu_si512(L->n);
    *vrho = _mm512_loadu_si512(L->rho);
}

static const int LV[7] = { 27, 47, 70, 85, 125, 165, 205 };

static void init_consts(void)
{
    dbias = _mm512_castsi512_pd(uecm_set64(0x4670000000000000ULL));
    vbias1 = uecm_set64(0x4670000000000000ULL);
    vbias2 = uecm_set64(0x4670000000000001ULL);
    lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
}

static int test(int ntests)
{
    int bad = 0, tot = 0;
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        for (int t = 0; t < ntests; t++) {
            lanes_t L; ecurve8 cv; V n, rho;
            setup_lanes(&L, 20 + (t % 32), 1, &cv, &n, &rho);
            V Xe, Ze;
            edw_stage1_8(&cv, &ch, &Xe, &Ze, n, rho);
            uecm_pt_x8 P;
            P.X = ADD(cv.P.Z, cv.P.Y);
            P.Z = SUB(cv.P.Z, cv.P.Y);
            ALIGNED_MEM uint64_t s[8];
            for (int i = 0; i < 8; i++)
                s[i] = to_m52(inv_plain((1 + L.d[i]) % L.n[i], L.n[i]), L.n[i]);
            uecm_stage1_x8(rho, n, &P, (uint64_t)B1, _mm512_loadu_si512(s));
            ALIGNED_MEM uint64_t a[8], b[8];
            _mm512_storeu_si512(a, MUL(Xe, P.Z));
            _mm512_storeu_si512(b, MUL(P.X, Ze));
            for (int i = 0; i < 8; i++) {
                int ok = a[i] == b[i];
                if (!ok && bad < 5)
                    printf("mismatch B1 %d lane %d n %llu\n", B1, i, (unsigned long long)L.n[i]);
                bad += !ok;
                tot++;
            }
        }
    }
    printf("test: %d of %d lanes mismatch\n", bad, tot);
    return bad;
}

#define NSET 64

static void bench(int nvec, int nbits)
{
    static lanes_t L[NSET];
    static ecurve8 cv[NSET];
    static V vn[NSET], vr[NSET], vs[NSET];
    static uecm_pt_x8 P0[NSET];
    V sink = _mm512_setzero_si512();
    for (int k = 0; k < NSET; k++) {
        setup_lanes(&L[k], nbits, 0, &cv[k], &vn[k], &vr[k]);
        ALIGNED_MEM uint64_t a[8], b[8];
        for (int i = 0; i < 8; i++) {
            a[i] = to_m52(rnd64(), L[k].n[i]);
            b[i] = to_m52(rnd64(), L[k].n[i]);
        }
        P0[k].X = _mm512_loadu_si512(a);
        P0[k].Z = cv[k].P.Z;
        vs[k] = _mm512_loadu_si512(b);
    }
    printf("%d-bit n; ns per vector of 8 curves (stage 1)\n", nbits);
    printf("  B1 |    PRAC     Edw   ratio | build: Suyama x8  Edw load x8\n");
    for (int li = 0; li < 7; li++) {
        int B1 = LV[li];
        chain_t ch = chain_for(B1);
        double best[2] = { 1e30, 1e30 };
        for (int rep = 0; rep < 3; rep++) {
            double t0 = now();
            for (int k = 0; k < nvec; k++) {
                int i = k & (NSET - 1);
                uecm_pt_x8 P = P0[i];
                uecm_stage1_x8(vr[i], vn[i], &P, (uint64_t)B1, vs[i]);
                sink = _mm512_xor_si512(sink, P.Z);
            }
            double t1 = now();
            for (int k = 0; k < nvec; k++) {
                int i = k & (NSET - 1);
                V X, Z;
                edw_stage1_8(&cv[i], &ch, &X, &Z, vn[i], vr[i]);
                sink = _mm512_xor_si512(sink, Z);
            }
            double t2 = now();
            if (t1 - t0 < best[0]) best[0] = t1 - t0;
            if (t2 - t1 < best[1]) best[1] = t2 - t1;
        }
        best[0] *= 1e9 / nvec;
        best[1] *= 1e9 / nvec;
        double ts2 = 1e30;
        for (int rep = 0; rep < 3; rep++) {
            double t0 = now();
            for (int k = 0; k < nvec; k++) {
                int i = k & (NSET - 1);
                uecm_pt_x8 P = P0[i];
                V a = uecm_stage2_x8(&P, vr[i], vn[i], (uint32_t)B1, vs[i], _mm512_loadu_si512(L[i].one));
                sink = _mm512_xor_si512(sink, a);
            }
            double t1 = now();
            if (t1 - t0 < ts2) ts2 = t1 - t0;
        }
        printf("  [stage 2: %.0f ns]", ts2 * 1e9 / nvec);
        if (li == 0) {
            uint64_t lcg = 1;
            int nb = nvec * 2;
            ALIGNED_MEM uint64_t lg[8];
            double t0 = now();
            for (int k = 0; k < nb; k++) {
                int i = k & (NSET - 1);
                uecm_pt_x8 P;
                V s = uecm_build_x8(&P, vr[i], vn[i], &lcg, lg,
                    _mm512_loadu_si512(L[i].five), _mm512_loadu_si512(L[i].Rsqr));
                sink = _mm512_xor_si512(sink, s);
            }
            double t1 = now();
            for (int k = 0; k < nb; k++) {
                int i = k & (NSET - 1);
                V n = vn[i], rho = vr[i], R2 = _mm512_loadu_si512(L[i].Rsqr);
                V c = _mm512_set1_epi64(k + 3);
                V X = MUL(c, R2), Y = MUL(_mm512_add_epi64(c, c), R2);
                V Z = MUL(_mm512_add_epi64(c, n), R2), T = MUL(X, Y);
                V dn = MUL(_mm512_sub_epi64(n, c), R2), dd = MUL(_mm512_srli_epi64(n, 3), R2);
                sink = _mm512_xor_si512(sink, _mm512_xor_si512(_mm512_xor_si512(X, Y),
                    _mm512_xor_si512(_mm512_xor_si512(Z, T), _mm512_xor_si512(dn, dd))));
            }
            double t2 = now();
            printf("%4d | %7.0f %7.0f %7.3f | %8.0f         %6.0f\n", B1, best[0], best[1],
                best[0] / best[1], (t1 - t0) * 1e9 / nb, (t2 - t1) * 1e9 / nb);
        } else {
            printf("%4d | %7.0f %7.0f %7.3f |\n", B1, best[0], best[1], best[0] / best[1]);
        }
    }
    ALIGNED_MEM uint64_t o[8];
    _mm512_storeu_si512(o, sink);
    if (o[0] == 42) printf(" ");
}

int main(int argc, char **argv)
{
    init_consts();
    if (argc > 1 && strcmp(argv[1], "-t") == 0)
        return test(argc > 2 ? atoi(argv[2]) : 100) != 0;
    bench(argc > 1 ? atoi(argv[1]) : 20000, argc > 2 ? atoi(argv[2]) : 50);
    return 0;
}
