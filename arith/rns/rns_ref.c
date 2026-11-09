/*
 * RNS Montgomery modular multiplication: setup and scalar reference.
 * See rns.h and arith/rns_design.md.
 */

#include "rns.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// 64x64 -> 128-bit products: unsigned __int128 where available, else the
// MSVC x64 intrinsics.
#if defined(__SIZEOF_INT128__)
typedef unsigned __int128 u128;
static inline void mul64(uint64_t x, uint64_t y, uint64_t *lo, uint64_t *hi)
{
    u128 p = (u128)x * y;
    *lo = (uint64_t)p;
    *hi = (uint64_t)(p >> 64);
}
// x*y mod m, for x, y < m
static inline uint64_t mulmod64(uint64_t x, uint64_t y, uint64_t m)
{
    return (uint64_t)(((u128)x * y) % m);
}
#elif defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
static inline void mul64(uint64_t x, uint64_t y, uint64_t *lo, uint64_t *hi)
{
    *lo = _umul128(x, y, hi);
}
static inline uint64_t mulmod64(uint64_t x, uint64_t y, uint64_t m)
{
    uint64_t hi, lo = _umul128(x, y, &hi), r;
    _udiv128(hi, lo, m, &r);    // hi < m since x, y < m
    return r;
}
#else
#error "rns_ref.c needs unsigned __int128 or the MSVC x64 intrinsics"
#endif

static void *rns_xmalloc(size_t sz)
{
    void *p = calloc(1, sz);
    if (p == NULL)
    {
        fprintf(stderr, "rns: out of memory\n");
        exit(1);
    }
    return p;
}

static uint64_t gcd64(uint64_t a, uint64_t b)
{
    while (b)
    {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static uint64_t mpz_mod_u52(mpz_t x, uint64_t m, mpz_t t)
{
    mpz_set_ui(t, m);
    mpz_mod(t, x, t);
    return mpz_get_ui(t);
}

static uint64_t mpz_low52(mpz_t x, mpz_t t)
{
    // x mod 2^52 for x of any sign
    mpz_fdiv_r_2exp(t, x, RNS_WBITS);
    return mpz_get_ui(t);
}

// ---------------------------------------------------------------------------
// channel reduction for m = 2^52 - d, d < 2^20.  Input lo + hi*2^52 with
// lo, hi < 2^62.  Output < 2^52, congruent mod m (partially reduced).
// This mirrors the intended IFMA kernel: a lazily accumulated lo/hi pair
// reduced once.
// ---------------------------------------------------------------------------
static inline uint64_t red116(uint64_t lo, uint64_t hi, uint64_t d)
{
    uint64_t a = lo & RNS_MASK52;
    uint64_t bv = hi + (lo >> 52);
    uint64_t b0 = bv & RNS_MASK52;
    uint64_t b1 = bv >> 52;                 // < 2^12
    uint64_t plo, phi;
    mul64(b0, d, &plo, &phi);
    uint64_t pl = plo & RNS_MASK52;
    uint64_t ph = (phi << 12) | (plo >> 52);    // < d
    // 2^52 == d, 2^104 == d^2 (mod m)
    uint64_t s = a + pl + (b1 * d + ph) * d;    // < 3*2^52
    s = (s & RNS_MASK52) + (s >> 52) * d;       // < 2^52 + 3d
    s = (s & RNS_MASK52) + (s >> 52) * d;       // < 2^52
    return s;
}

static inline uint64_t mulred(uint64_t x, uint64_t y, uint64_t d)
{
    uint64_t plo, phi;
    mul64(x, y, &plo, &phi);
    return red116(plo & RNS_MASK52, (phi << 12) | (plo >> 52), d);
}

// fold a value < 2^62 to < 2^52 (mod 2^52 - d)
static inline uint64_t redsmall(uint64_t t, uint64_t d)
{
    t = (t & RNS_MASK52) + (t >> 52) * d;
    t = (t & RNS_MASK52) + (t >> 52) * d;
    return t;
}

// lazy accumulator of 52x52 products as lo/hi 52-bit halves (madd52lo/hi)
#define ACC_MADD(lo, hi, x, y) do {                     \
        uint64_t _pl, _ph;                              \
        mul64((x), (y), &_pl, &_ph);                    \
        (lo) += _pl & RNS_MASK52;                       \
        (hi) += (_ph << 12) | (_pl >> 52);              \
    } while (0)

// ---------------------------------------------------------------------------
// fixed setup
// ---------------------------------------------------------------------------

// FCFS moduli 2^52 - d, d odd ascending, pairwise coprime
static int select_moduli(uint64_t *dpool, int count)
{
    uint64_t *m = (uint64_t *)rns_xmalloc((size_t)(count) * sizeof(uint64_t));
    int n = 0;
    uint64_t d;

    for (d = 1; n < count; d += 2)
    {
        uint64_t cand = (1ULL << 52) - d;
        int i, ok = 1;

        if (d >= (1ULL << 20))
        {
            free(m);
            return n;
        }
        for (i = 0; i < n; i++)
        {
            if (gcd64(cand, m[i]) != 1)
            {
                ok = 0;
                break;
            }
        }
        if (ok)
        {
            dpool[n] = d;
            m[n++] = cand;
        }
    }
    free(m);
    return n;
}

rns_fixed_t *rns_fixed_init(int maxbits)
{
    rns_fixed_t *fx = (rns_fixed_t *)rns_xmalloc(sizeof(rns_fixed_t));
    // every modulus exceeds 2^(52 - 2^-31); use 51.999 bits/channel
    const double wb = 51.999;
    int nB, nC, i, j, npool;
    uint64_t *dpool;
    mpz_t t, u, bound;

    fx->maxbits = maxbits;

    // M_B >= h^2 k^2 2^maxbits, k = nB + 2
    for (nB = 1; ; nB++)
    {
        double need = maxbits + 2.0 * log2((double)RNS_HEADROOM * (nB + 2)) + 1;
        if (wb * nB >= need)
            break;
    }
    // M_C > max U = k*N, with margin
    for (nC = 1; ; nC++)
    {
        double need = maxbits + log2(2.0 * (nB + 2)) + 1;
        if (wb * nC >= need)
            break;
    }

    fx->nB = nB;
    fx->nC = nC;
    fx->nch = nB + nC + 1;
    fx->k = nB + 2;

    npool = nB + nC;
    dpool = (uint64_t *)rns_xmalloc((size_t)(npool) * sizeof(uint64_t));
    if (select_moduli(dpool, npool) < npool)
    {
        fprintf(stderr, "rns: not enough moduli 2^52-d with d < 2^20\n");
        exit(1);
    }

    fx->b = (uint64_t *)rns_xmalloc((size_t)(nB) * sizeof(uint64_t));
    fx->db = (uint64_t *)rns_xmalloc((size_t)(nB) * sizeof(uint64_t));
    fx->c = (uint64_t *)rns_xmalloc((size_t)(nC) * sizeof(uint64_t));
    fx->dc = (uint64_t *)rns_xmalloc((size_t)(nC) * sizeof(uint64_t));

    // interleave the pool between B and C, the remainder to the larger base
    {
        int ib = 0, ic = 0, p = 0;
        while (ib < nB || ic < nC)
        {
            if (ib < nB && (ib <= ic || ic >= nC))
                fx->db[ib++] = dpool[p++];
            else
                fx->dc[ic++] = dpool[p++];
        }
    }
    for (i = 0; i < nB; i++)
        fx->b[i] = (1ULL << 52) - fx->db[i];
    for (j = 0; j < nC; j++)
        fx->c[j] = (1ULL << 52) - fx->dc[j];
    free(dpool);

    mpz_init(t);
    mpz_init(u);
    mpz_init(bound);
    mpz_init_set_ui(fx->MB, 1);
    mpz_init_set_ui(fx->MC, 1);
    for (i = 0; i < nB; i++)
        mpz_mul_ui(fx->MB, fx->MB, fx->b[i]);
    for (j = 0; j < nC; j++)
        mpz_mul_ui(fx->MC, fx->MC, fx->c[j]);

    // verify the bounds exactly
    mpz_set_ui(bound, (unsigned long)(RNS_HEADROOM * fx->k));
    mpz_mul(bound, bound, bound);
    mpz_mul_2exp(bound, bound, (mp_bitcnt_t)maxbits);
    if (mpz_cmp(fx->MB, bound) < 0)
    {
        fprintf(stderr, "rns: M_B bound check failed\n");
        exit(1);
    }
    mpz_set_ui(bound, (unsigned long)(2 * fx->k));
    mpz_mul_2exp(bound, bound, (mp_bitcnt_t)maxbits);
    if (mpz_cmp(fx->MC, bound) < 0)
    {
        fprintf(stderr, "rns: M_C bound check failed\n");
        exit(1);
    }

    fx->Bi = (mpz_t *)rns_xmalloc((size_t)(nB) * sizeof(mpz_t));
    fx->Biinv_b = (uint64_t *)rns_xmalloc((size_t)(nB) * sizeof(uint64_t));
    fx->Cjinv_c = (uint64_t *)rns_xmalloc((size_t)(nC) * sizeof(uint64_t));
    fx->BtoC = (uint64_t *)rns_xmalloc((size_t)((nC + 1) * nB) * sizeof(uint64_t));
    fx->CtoB = (uint64_t *)rns_xmalloc((size_t)((nB + 1) * nC) * sizeof(uint64_t));
    fx->negMC_b = (uint64_t *)rns_xmalloc((size_t)(nB) * sizeof(uint64_t));
    fx->CjMBinv_c = (uint64_t *)rns_xmalloc((size_t)(nC) * sizeof(uint64_t));

    for (i = 0; i < nB; i++)
    {
        mpz_init(fx->Bi[i]);
        mpz_divexact_ui(fx->Bi[i], fx->MB, fx->b[i]);
        mpz_set_ui(u, fx->b[i]);
        mpz_invert(t, fx->Bi[i], u);
        fx->Biinv_b[i] = mpz_get_ui(t);
        for (j = 0; j < nC; j++)
            fx->BtoC[j * nB + i] = mpz_mod_u52(fx->Bi[i], fx->c[j], t);
        fx->BtoC[nC * nB + i] = mpz_low52(fx->Bi[i], t);
    }

    for (j = 0; j < nC; j++)
    {
        mpz_t Cj;
        mpz_init(Cj);
        mpz_divexact_ui(Cj, fx->MC, fx->c[j]);
        mpz_set_ui(u, fx->c[j]);
        mpz_invert(t, Cj, u);
        fx->Cjinv_c[j] = mpz_get_ui(t);
        for (i = 0; i < nB; i++)
            fx->CtoB[i * nC + j] = mpz_mod_u52(Cj, fx->b[i], t);
        fx->CtoB[nB * nC + j] = mpz_low52(Cj, t);

        // |C_j * M_B^{-1}|_{c_j}
        mpz_invert(t, fx->MB, u);
        mpz_mul(t, t, Cj);
        fx->CjMBinv_c[j] = mpz_mod_u52(t, fx->c[j], u);
        mpz_clear(Cj);
    }

    for (i = 0; i < nB; i++)
    {
        mpz_neg(t, fx->MC);
        fx->negMC_b[i] = mpz_mod_u52(t, fx->b[i], u);
    }

    mpz_set_ui(u, 1);
    mpz_mul_2exp(u, u, RNS_WBITS);
    mpz_invert(t, fx->MC, u);
    fx->MCinv_r = mpz_get_ui(t);
    mpz_invert(t, fx->MB, u);
    fx->MBinv_r = mpz_get_ui(t);

    mpz_clear(t);
    mpz_clear(u);
    mpz_clear(bound);
    return fx;
}

void rns_fixed_free(rns_fixed_t *fx)
{
    int i;
    for (i = 0; i < fx->nB; i++)
        mpz_clear(fx->Bi[i]);
    free(fx->Bi);
    free(fx->b); free(fx->db);
    free(fx->c); free(fx->dc);
    free(fx->BtoC); free(fx->CtoB);
    free(fx->negMC_b); free(fx->CjMBinv_c);
    free(fx->Biinv_b); free(fx->Cjinv_c);
    mpz_clear(fx->MB);
    mpz_clear(fx->MC);
    free(fx);
}

// ---------------------------------------------------------------------------
// variable (N-dependent) setup
// ---------------------------------------------------------------------------

uint64_t *rns_alloc(rns_ctx_t *ctx)
{
    return (uint64_t *)rns_xmalloc((size_t)(ctx->fx->nch) * sizeof(uint64_t));
}

uint64_t *rns_aligned_alloc(size_t words)
{
    void *p;
#ifdef _WIN32
    p = _aligned_malloc(words * sizeof(uint64_t), 64);
#else
    if (posix_memalign(&p, 64, words * sizeof(uint64_t)) != 0)
        p = NULL;
#endif
    if (p == NULL)
    {
        fprintf(stderr, "rns: out of memory\n");
        exit(1);
    }
    memset(p, 0, words * sizeof(uint64_t));
    return (uint64_t *)p;
}

uint64_t *rns_valloc(rns_ctx_t *ctx)
{
    return rns_aligned_alloc((size_t)ctx->fx->nch * RNS_VL);
}

void rns_vfree(uint64_t *v)
{
#ifdef _WIN32
    _aligned_free(v);
#else
    free(v);
#endif
}

void rns_vset_lane(rns_ctx_t *ctx, uint64_t *v, int lane, uint64_t *x)
{
    int ch;
    for (ch = 0; ch < ctx->fx->nch; ch++)
        v[ch * RNS_VL + lane] = x[ch];
}

void rns_vget_lane(rns_ctx_t *ctx, uint64_t *x, uint64_t *v, int lane)
{
    int ch;
    for (ch = 0; ch < ctx->fx->nch; ch++)
        x[ch] = v[ch * RNS_VL + lane];
}

rns_ctx_t *rns_ctx_init(rns_fixed_t *fx, mpz_t N)
{
    rns_ctx_t *ctx;
    int nB = fx->nB, nC = fx->nC, i, j;
    mpz_t t, u, Ninv;

    if (mpz_even_p(N) || (int)mpz_sizeinbase(N, 2) > fx->maxbits)
    {
        fprintf(stderr, "rns: N must be odd and < 2^%d\n", fx->maxbits);
        return NULL;
    }

    mpz_init(t);
    mpz_init(u);
    mpz_init(Ninv);

    mpz_gcd(t, N, fx->MB);
    mpz_gcd(u, N, fx->MC);
    if (mpz_cmp_ui(t, 1) != 0 || mpz_cmp_ui(u, 1) != 0)
    {
        // the moduli are pairwise coprime but not prime; a shared factor
        // is a found factor of N, left to the caller (gcd with M_B*M_C)
        mpz_clear(t); mpz_clear(u); mpz_clear(Ninv);
        return NULL;
    }

    ctx = (rns_ctx_t *)rns_xmalloc(sizeof(rns_ctx_t));
    ctx->fx = fx;
    mpz_init_set(ctx->N, N);
    ctx->f1_b = (uint64_t *)rns_xmalloc((size_t)(nB) * sizeof(uint64_t));
    ctx->g2_c = (uint64_t *)rns_xmalloc((size_t)(nC) * sizeof(uint64_t));

    for (i = 0; i < nB; i++)
    {
        // |-N^{-1} * B_i^{-1}|_{b_i}
        mpz_set_ui(u, fx->b[i]);
        mpz_invert(Ninv, N, u);
        mpz_mul_ui(t, Ninv, fx->Biinv_b[i]);
        mpz_neg(t, t);
        ctx->f1_b[i] = mpz_mod_u52(t, fx->b[i], u);
    }

    for (j = 0; j < nC; j++)
    {
        // |N * M_B^{-1} * C_j^{-1}|_{c_j}
        mpz_set_ui(u, fx->c[j]);
        mpz_invert(t, fx->MB, u);
        mpz_mul(t, t, N);
        mpz_mul_ui(t, t, fx->Cjinv_c[j]);
        ctx->g2_c[j] = mpz_mod_u52(t, fx->c[j], u);
    }

    mpz_set_ui(u, 1);
    mpz_mul_2exp(u, u, RNS_WBITS);
    mpz_invert(t, fx->MB, u);
    mpz_mul(t, t, N);
    ctx->NMBinv_r = mpz_low52(t, u);

    // Gandino reorganization: fold |N M_B^{-1} C_j^{-1}|_{c_j} into the
    // first extension matrix, so the extension sum and the T term share one
    // lazy accumulator.  Row nC is the redundant channel (mod 2^52).
    ctx->BtoCN = (uint64_t *)rns_xmalloc((size_t)((nC + 1) * nB) * sizeof(uint64_t));
    for (j = 0; j < nC; j++)
        for (i = 0; i < nB; i++)
            ctx->BtoCN[j * nB + i] = mulmod64(fx->BtoC[j * nB + i],
                ctx->g2_c[j], fx->c[j]);
    for (i = 0; i < nB; i++)
        ctx->BtoCN[nC * nB + i] = (fx->BtoC[nC * nB + i] * ctx->NMBinv_r) & RNS_MASK52;

    mpz_init(ctx->MBmodN);
    mpz_init(ctx->MBinvmodN);
    mpz_mod(ctx->MBmodN, fx->MB, N);
    mpz_invert(ctx->MBinvmodN, fx->MB, N);

    ctx->subconst = rns_alloc(ctx);
    mpz_mul_ui(t, N, (unsigned long)(RNS_SUBMULT * fx->k));
    rns_from_int(ctx, ctx->subconst, t);

    ctx->one = rns_alloc(ctx);
    mpz_set_ui(t, 1);
    rns_to_monty(ctx, ctx->one, t);

    ctx->vscratch = rns_aligned_alloc((size_t)nB * RNS_VL);

    mpz_clear(t);
    mpz_clear(u);
    mpz_clear(Ninv);
    return ctx;
}

void rns_ctx_free(rns_ctx_t *ctx)
{
    mpz_clear(ctx->N);
    mpz_clear(ctx->MBmodN);
    mpz_clear(ctx->MBinvmodN);
    free(ctx->f1_b);
    free(ctx->g2_c);
    free(ctx->BtoCN);
    free(ctx->subconst);
    free(ctx->one);
    rns_vfree(ctx->vscratch);
    free(ctx);
}

// ---------------------------------------------------------------------------
// conversions (GMP, not performance critical)
// ---------------------------------------------------------------------------

void rns_from_int(rns_ctx_t *ctx, uint64_t *x, mpz_t X)
{
    rns_fixed_t *fx = ctx->fx;
    int nB = fx->nB, nC = fx->nC, i, j;
    mpz_t t;

    mpz_init(t);
    for (i = 0; i < nB; i++)
        x[i] = mpz_mod_u52(X, fx->b[i], t);
    for (j = 0; j < nC; j++)
        x[nB + j] = mulred(mpz_mod_u52(X, fx->c[j], t), fx->Cjinv_c[j], fx->dc[j]);
    x[nB + nC] = mpz_low52(X, t);
    mpz_clear(t);
}

void rns_to_int(rns_ctx_t *ctx, mpz_t X, uint64_t *x)
{
    // CRT from the B channels
    rns_fixed_t *fx = ctx->fx;
    int i;
    mpz_t t;

    mpz_init(t);
    mpz_set_ui(X, 0);
    for (i = 0; i < fx->nB; i++)
    {
        uint64_t v = mulred(x[i], fx->Biinv_b[i], fx->db[i]);
        mpz_mul_ui(t, fx->Bi[i], v);
        mpz_add(X, X, t);
    }
    mpz_mod(X, X, fx->MB);
    mpz_clear(t);
}

void rns_to_monty(rns_ctx_t *ctx, uint64_t *x, mpz_t a)
{
    mpz_t t;
    mpz_init(t);
    mpz_mul(t, a, ctx->MBmodN);
    mpz_mod(t, t, ctx->N);
    rns_from_int(ctx, x, t);
    mpz_clear(t);
}

void rns_from_monty(rns_ctx_t *ctx, mpz_t a, uint64_t *x)
{
    rns_to_int(ctx, a, x);
    mpz_mul(a, a, ctx->MBinvmodN);
    mpz_mod(a, a, ctx->N);
}

// ---------------------------------------------------------------------------
// scalar reference arithmetic
// ---------------------------------------------------------------------------

void rns_mulmod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    rns_fixed_t *fx = ctx->fx;
    int nB = fx->nB, nC = fx->nC, i, j;
    const uint64_t *xc = x + nB, *yc = y + nB;
    uint64_t xi[1024], xr[1024];    // xi_i, xi'_j
    uint64_t lo, hi, ur, sr, alpha;

    if (nB > 1024 || nC > 1024)
    {
        fprintf(stderr, "rns: too many channels for the reference\n");
        exit(1);
    }

    // 1-2: T = X*Y in B, xi_i = |T_i * (-N^{-1} B_i^{-1})|_{b_i}
    for (i = 0; i < nB; i++)
        xi[i] = mulred(mulred(x[i], y[i], fx->db[i]), ctx->f1_b[i], fx->db[i]);

    // 3-4: first extension merged with U = (T + Q^ N) / M_B, in C
    // (pre-scaled) and r:  xi'_j = T'_j |C_j M_B^{-1}| + sum_i xi_i BtoCN_ji,
    // one lazy accumulator per channel
    for (j = 0; j < nC; j++)
    {
        const uint64_t *row = ctx->BtoCN + j * nB;
        uint64_t t = mulred(xc[j], yc[j], fx->dc[j]);
        lo = hi = 0;
        ACC_MADD(lo, hi, t, fx->CjMBinv_c[j]);
        for (i = 0; i < nB; i++)
            ACC_MADD(lo, hi, xi[i], row[i]);
        xr[j] = red116(lo, hi, fx->dc[j]);
    }
    {
        const uint64_t *row = ctx->BtoCN + nC * nB;
        ur = (x[nB + nC] * y[nB + nC]) * fx->MBinv_r;
        for (i = 0; i < nB; i++)
            ur += xi[i] * row[i];
        ur &= RNS_MASK52;
    }

    // 5: exact second extension C -> B, alpha' from the redundant channel
    {
        const uint64_t *row = fx->CtoB + nB * nC;
        sr = 0;
        for (j = 0; j < nC; j++)
            sr += xr[j] * row[j];
        alpha = ((sr - ur) * fx->MCinv_r) & RNS_MASK52;
#ifdef RNS_DEBUG
        if (alpha > (uint64_t)nC + 1)
        {
            fprintf(stderr, "rns: alpha' = %llu out of range\n", (unsigned long long)alpha);
            exit(1);
        }
#endif
    }

    for (i = 0; i < nB; i++)
    {
        const uint64_t *row = fx->CtoB + i * nC;
        lo = hi = 0;
        for (j = 0; j < nC; j++)
            ACC_MADD(lo, hi, xr[j], row[j]);
        ACC_MADD(lo, hi, alpha, fx->negMC_b[i]);
        z[i] = red116(lo, hi, fx->db[i]);
    }
    for (j = 0; j < nC; j++)
        z[nB + j] = xr[j];
    z[nB + nC] = ur;
}

void rns_sqrmod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x)
{
    rns_mulmod_ref(ctx, z, x, x);
}

void rns_addmod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    rns_fixed_t *fx = ctx->fx;
    int nB = fx->nB, nC = fx->nC, i;

    for (i = 0; i < nB; i++)
        z[i] = redsmall(x[i] + y[i], fx->db[i]);
    for (i = 0; i < nC; i++)
        z[nB + i] = redsmall(x[nB + i] + y[nB + i], fx->dc[i]);
    z[nB + nC] = (x[nB + nC] + y[nB + nC]) & RNS_MASK52;
}

void rns_submod_ref(rns_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    // X - Y + SUBMULT*k*N.  2m > 2^52 > y keeps the channel sum positive.
    rns_fixed_t *fx = ctx->fx;
    int nB = fx->nB, nC = fx->nC, i;
    uint64_t *s = ctx->subconst;

    for (i = 0; i < nB; i++)
        z[i] = redsmall(x[i] + s[i] + 2 * fx->b[i] - y[i], fx->db[i]);
    for (i = 0; i < nC; i++)
        z[nB + i] = redsmall(x[nB + i] + s[nB + i] + 2 * fx->c[i] - y[nB + i], fx->dc[i]);
    z[nB + nC] = (x[nB + nC] + s[nB + nC] - y[nB + nC]) & RNS_MASK52;
}
