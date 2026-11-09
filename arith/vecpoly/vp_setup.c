/*
 * vecpoly setup: NTT primes, twiddles, conversion and CRT constants,
 * allocation, lane get/set.  See vecpoly.h.
 */

#include "vecpoly.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif

static uint64_t vp_mulmod(uint64_t a, uint64_t b, uint64_t p)
{
#if defined(__SIZEOF_INT128__)
    return (uint64_t)(((unsigned __int128)a * b) % p);
#else
    uint64_t hi, lo = _umul128(a, b, &hi), r;
    _udiv128(hi, lo, p, &r);
    return r;
#endif
}

// floor(w * 2^52 / p), w < p < 2^50
static uint64_t vp_shoup(uint64_t w, uint64_t p)
{
#if defined(__SIZEOF_INT128__)
    return (uint64_t)(((unsigned __int128)w << 52) / p);
#else
    uint64_t r;
    return _udiv128(w >> 12, w << 52, p, &r);
#endif
}

static uint64_t vp_powmod(uint64_t b, uint64_t e, uint64_t p)
{
    uint64_t r = 1;
    b %= p;
    while (e)
    {
        if (e & 1) r = vp_mulmod(r, b, p);
        b = vp_mulmod(b, b, p);
        e >>= 1;
    }
    return r;
}

static uint64_t vp_invmod(uint64_t a, uint64_t p)
{
    return vp_powmod(a, p - 2, p);
}

// -x^-1 mod 2^52, x odd
static uint64_t vp_neginv52(uint64_t x)
{
    uint64_t y = 1;
    int i;
    for (i = 0; i < 6; i++)
        y = y * (2 - x * y);        // Newton, doubles the correct bits
    return (0 - y) & VP_MASK52;
}

void *vp_aligned_alloc(size_t bytes)
{
    void *p = NULL;
    if (bytes == 0) bytes = 64;
#if defined(_MSC_VER) || defined(__MINGW32__)
    p = _aligned_malloc(bytes, 64);
#else
    if (posix_memalign(&p, 64, bytes) != 0) p = NULL;
#endif
    if (p == NULL)
    {
        fprintf(stderr, "vecpoly: out of memory (%zu bytes)\n", bytes);
        exit(1);
    }
    return p;
}

void vp_aligned_free(void *p)
{
    if (p == NULL) return;
#if defined(_MSC_VER) || defined(__MINGW32__)
    _aligned_free(p);
#else
    free(p);
#endif
}

uint64_t *vp_alloc(size_t ncoef, int n)
{
    size_t bytes = ncoef * 8 * (size_t)n * sizeof(uint64_t);
    uint64_t *a = (uint64_t *)vp_aligned_alloc(bytes);
    memset(a, 0, bytes);
    return a;
}

void vp_free(uint64_t *a)
{
    vp_aligned_free(a);
}

// portable for 32-bit unsigned long (MSVC + MPIR)
static void vp_mpz_set_u64(mpz_t x, uint64_t v)
{
    mpz_import(x, 1, -1, sizeof(uint64_t), 0, 0, &v);
}

static uint64_t vp_mpz_get_u64(mpz_t x)     // 0 <= x < 2^64
{
    uint64_t w = 0;
    size_t cnt = 0;
    mpz_export(&w, &cnt, -1, sizeof(uint64_t), 0, 0, x);
    return w;
}

static void vp_mpz_to_words(uint64_t *w, int n, mpz_t x)
{
    mpz_t t, u;
    int k;
    mpz_init_set(t, x);
    mpz_init(u);
    for (k = 0; k < n; k++)
    {
        mpz_tdiv_r_2exp(u, t, 52);
        w[k] = vp_mpz_get_u64(u);
        mpz_tdiv_q_2exp(t, t, 52);
    }
    mpz_clear(t);
    mpz_clear(u);
}

static uint64_t vp_mpz_mod_ui64(mpz_t x, uint64_t p)
{
    mpz_t t, u;
    uint64_t r;
    mpz_init(t);
    mpz_init(u);
    vp_mpz_set_u64(t, p);
    mpz_mod(u, x, t);
    r = vp_mpz_get_u64(u);
    mpz_clear(t);
    mpz_clear(u);
    return r;
}

void vp_get(vp_ctx_t *ctx, mpz_t x, const uint64_t *c, int lane)
{
    int k;
    mpz_t t;
    mpz_init(t);
    mpz_set_ui(x, 0);
    for (k = ctx->n - 1; k >= 0; k--)
    {
        mpz_mul_2exp(x, x, 52);
        vp_mpz_set_u64(t, c[k * 8 + lane]);
        mpz_add(x, x, t);
    }
    mpz_clear(t);
}

void vp_set(vp_ctx_t *ctx, uint64_t *c, mpz_t x, int lane)
{
    uint64_t *w = (uint64_t *)malloc((size_t)ctx->n * sizeof(uint64_t));
    int k;
    vp_mpz_to_words(w, ctx->n, x);
    for (k = 0; k < ctx->n; k++)
        c[k * 8 + lane] = w[k];
    free(w);
}

static int vp_is_prime(uint64_t p)
{
    mpz_t t;
    int r;
    mpz_init(t);
    vp_mpz_set_u64(t, p);
    r = mpz_probab_prime_p(t, 25);
    mpz_clear(t);
    return r;
}

vp_ctx_t *vp_ctx_init(mpz_t N, int n, int logLmax)
{
    vp_ctx_t *ctx;
    mpz_t P, R, t, u, bound;
    uint64_t c, p;
    int i, k, l, K;
    uint32_t Lmax;
    uint64_t *plist;

    if ((logLmax < 1) || (logLmax > VP_MAX_LOGL) || (n < 2) ||
        (mpz_sizeinbase(N, 2) + 3 > (size_t)(52 * n)) || mpz_even_p(N))
        return NULL;

    ctx = (vp_ctx_t *)calloc(1, sizeof(vp_ctx_t));
    ctx->n = n;
    ctx->logLmax = logLmax;
    ctx->Lmax = Lmax = 1u << logLmax;
    ctx->school_max = 8;
    mpz_init_set(ctx->N, N);

    mpz_init(P);
    mpz_init(R);
    mpz_init(t);
    mpz_init(u);
    mpz_init(bound);

    // bound: P > 4 * Lmax * R^2, R = 2^(52n)
    mpz_set_ui(R, 1);
    mpz_mul_2exp(R, R, (mp_bitcnt_t)(52 * n));
    mpz_mul(bound, R, R);
    mpz_mul_2exp(bound, bound, (mp_bitcnt_t)(logLmax + 2));

    // primes p = c * 2^20 + 1 < 2^50, largest first
    plist = (uint64_t *)malloc(4096 * sizeof(uint64_t));
    K = 0;
    mpz_set_ui(P, 1);
    for (c = (1ULL << (VP_PRIME_BITS - VP_MAX_LOGL)) - 1;
        c > (1ULL << (VP_PRIME_BITS - 1 - VP_MAX_LOGL)); c--)
    {
        p = (c << VP_MAX_LOGL) + 1;
        if (!vp_is_prime(p)) continue;
        plist[K++] = p;
        vp_mpz_set_u64(t, p);
        mpz_mul(P, P, t);
        if (mpz_cmp(P, bound) > 0) break;
        if (K == 4096) break;
    }
    if (mpz_cmp(P, bound) <= 0)
    {
        free(plist);
        free(ctx);
        return NULL;
    }
    ctx->K = K;

    ctx->p = (uint64_t *)malloc((size_t)K * sizeof(uint64_t));
    ctx->pneg = (uint64_t *)malloc((size_t)K * sizeof(uint64_t));
    ctx->pinvd = (double *)malloc((size_t)K * sizeof(double));
    ctx->cvt = (uint64_t *)malloc((size_t)K * (size_t)n * sizeof(uint64_t));
    ctx->crt = (uint64_t *)malloc((size_t)(K + 1) * (size_t)n * sizeof(uint64_t));
    ctx->yscale = (uint64_t *)malloc((size_t)(logLmax + 1) * (size_t)K * 2 * sizeof(uint64_t));
    ctx->tw = (uint64_t *)vp_aligned_alloc((size_t)K * Lmax * sizeof(uint64_t));
    ctx->itw = (uint64_t *)vp_aligned_alloc((size_t)K * Lmax * sizeof(uint64_t));
    ctx->Nw = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    ctx->N2w = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));

    vp_mpz_to_words(ctx->Nw, n, N);
    mpz_mul_2exp(t, N, 1);
    vp_mpz_to_words(ctx->N2w, n, t);
    ctx->n0 = vp_neginv52(ctx->Nw[0]);

    for (i = 0; i < K; i++)
    {
        uint64_t g, h, w, iw;
        uint32_t j;

        p = ctx->p[i] = plist[i];
        ctx->pneg[i] = vp_neginv52(p);
        ctx->pinvd[i] = 1.0 / (double)p;

        // conversion: 2^(52k + 104) mod p
        h = vp_powmod(2, 104, p);
        w = vp_powmod(2, 52, p);
        for (k = 0; k < n; k++)
        {
            ctx->cvt[i * n + k] = h;
            h = vp_mulmod(h, w, p);
        }

        // element of order 2^20
        for (g = 3; ; g++)
        {
            h = vp_powmod(g, (p - 1) >> VP_MAX_LOGL, p);
            if (vp_powmod(h, 1ULL << (VP_MAX_LOGL - 1), p) == p - 1)
                break;
        }

        // twiddles: powers of an Lmax-th root of unity, j < Lmax/2; a level
        // with 2m-th roots reads every (Lmax/2m)-th entry
        {
            uint64_t wm = vp_powmod(h, (1ULL << VP_MAX_LOGL) / Lmax, p);
            uint64_t iwm = vp_invmod(wm, p);
            uint64_t *tp = ctx->tw + (size_t)i * Lmax;
            uint64_t *itp = ctx->itw + (size_t)i * Lmax;
            w = 1;
            iw = 1;
            for (j = 0; j < Lmax / 2; j++)
            {
                tp[2 * j] = w;
                tp[2 * j + 1] = vp_shoup(w, p);
                itp[2 * j] = iw;
                itp[2 * j + 1] = vp_shoup(iw, p);
                w = vp_mulmod(w, wm, p);
                iw = vp_mulmod(iw, iwm, p);
            }
        }

        // CRT: (P/p_i)^-1 mod p_i, then per length L^-1 2^52 (P/p_i)^-1
        vp_mpz_set_u64(t, p);
        mpz_divexact(u, P, t);          // u = P / p_i
        h = vp_invmod(vp_mpz_mod_ui64(u, p), p);
        h = vp_mulmod(h, vp_powmod(2, 52, p), p);
        for (l = 0; l <= logLmax; l++)
        {
            uint64_t s = vp_mulmod(h, vp_invmod(vp_powmod(2, (uint64_t)l, p), p), p);
            ctx->yscale[((size_t)l * (size_t)K + (size_t)i) * 2] = s;
            ctx->yscale[((size_t)l * (size_t)K + (size_t)i) * 2 + 1] = vp_shoup(s, p);
        }

        // C_i = (P/p_i) R^-1 2^104 mod N
        mpz_invert(t, R, N);
        mpz_mul(u, u, t);
        mpz_mul_2exp(u, u, 104);
        mpz_mod(u, u, N);
        vp_mpz_to_words(ctx->crt + (size_t)i * (size_t)n, n, u);
    }

    // C_K = -P R^-1 2^104 mod N
    mpz_invert(t, R, N);
    mpz_mul(u, P, t);
    mpz_mul_2exp(u, u, 104);
    mpz_mod(u, u, N);
    mpz_sub(u, N, u);
    mpz_mod(u, u, N);
    vp_mpz_to_words(ctx->crt + (size_t)K * (size_t)n, n, u);

    free(plist);
    mpz_clear(P);
    mpz_clear(R);
    mpz_clear(t);
    mpz_clear(u);
    mpz_clear(bound);
    return ctx;
}

void vp_ctx_free(vp_ctx_t *ctx)
{
    if (ctx == NULL) return;
    mpz_clear(ctx->N);
    free(ctx->p);
    free(ctx->pneg);
    free(ctx->pinvd);
    free(ctx->cvt);
    free(ctx->crt);
    free(ctx->yscale);
    vp_aligned_free(ctx->tw);
    vp_aligned_free(ctx->itw);
    free(ctx->Nw);
    free(ctx->N2w);
    free(ctx);
}

vp_work_t *vp_work_init(vp_ctx_t *ctx)
{
    vp_work_t *w = (vp_work_t *)calloc(1, sizeof(vp_work_t));
    int n = ctx->n, i, k;
    size_t nttbytes = (size_t)ctx->K * ctx->Lmax * 8 * sizeof(uint64_t);
    size_t grpbytes = (size_t)VP_GRP * ctx->Lmax * 8 * sizeof(uint64_t);
    mpz_t t;

    w->ctx = ctx;
#if defined(USE_AVX512F) || defined(USE_AVX2)
    if (!vecmulmod52_tmpl_get(n, &w->mul, &w->sqr))
    {
        free(w);
        return NULL;
    }
#else
    free(w);
    return NULL;
#endif
    memset(&w->md, 0, sizeof(w->md));
    w->md.NWORDS = (uint32_t)n;
    w->md.MAXBITS = (uint32_t)(52 * n);
    w->md.vrho = (base_t *)vp_aligned_alloc(8 * sizeof(base_t));
    for (i = 0; i < 8; i++) w->md.vrho[i] = ctx->n0;
    w->md.fpwork = (base_t *)vp_aligned_alloc((size_t)3 * (size_t)n * 8 * sizeof(base_t) + 64);

    w->vn.data = (base_t *)vp_aligned_alloc((size_t)n * 8 * sizeof(base_t));
    w->vn.size = n;
    w->vn.WORDS_ALLOC = (uint32_t)n;
    for (k = 0; k < n; k++)
        for (i = 0; i < 8; i++)
            w->vn.data[k * 8 + i] = ctx->Nw[k];
    w->vs.data = (base_t *)vp_aligned_alloc((size_t)(2 * n + 4) * 8 * sizeof(base_t));
    w->vs.size = 2 * n + 4;
    w->vs.WORDS_ALLOC = (uint32_t)(2 * n + 4);

    w->ntt = (uint64_t *)vp_aligned_alloc(nttbytes);
    w->grp = (uint64_t *)vp_aligned_alloc(grpbytes);
    w->tmp[0] = vp_alloc(ctx->Lmax / 2, n);
    w->tmp[1] = vp_alloc(ctx->Lmax / 4 > 1 ? ctx->Lmax / 4 : 1, n);
    w->sm = vp_alloc(4 * VP_SCHOOL_CAP, n);

    // Montgomery one = R mod N
    w->one = vp_alloc(1, n);
    mpz_init(t);
    mpz_set_ui(t, 1);
    mpz_mul_2exp(t, t, (mp_bitcnt_t)(52 * n));
    mpz_mod(t, t, ctx->N);
    for (i = 0; i < 8; i++)
        vp_set(ctx, w->one, t, i);
    mpz_clear(t);
    return w;
}

void vp_work_free(vp_work_t *w)
{
    if (w == NULL) return;
    vp_aligned_free(w->md.vrho);
    vp_aligned_free(w->md.fpwork);
    vp_aligned_free(w->vn.data);
    vp_aligned_free(w->vs.data);
    vp_aligned_free(w->ntt);
    vp_aligned_free(w->grp);
    vp_free(w->tmp[0]);
    vp_free(w->tmp[1]);
    vp_free(w->sm);
    vp_free(w->one);
    free(w);
}

int vp_nprimes(int n, int logLmax)
{
    // P > 4 Lmax R^2 with every p within 2^-20 of 2^50
    return (104 * n + logLmax + 2) / VP_PRIME_BITS + 1;
}

size_t vp_ctx_bytes(int n, int K, int logLmax)
{
    size_t Lmax = (size_t)1 << logLmax;
    return 2 * (size_t)K * Lmax * sizeof(uint64_t)                 // tw, itw
        + ((size_t)(2 * K + 1) * (size_t)n + (size_t)(logLmax + 1) * (size_t)K * 2
            + (size_t)K * 3) * sizeof(uint64_t);
}

size_t vp_work_bytes(int n, int K, int logLmax)
{
    size_t Lmax = (size_t)1 << logLmax, cb = (size_t)64 * (size_t)n;
    return (size_t)(K + VP_GRP) * Lmax * 64                         // ntt, grp
        + (Lmax / 2 + Lmax / 4 + 4 * VP_SCHOOL_CAP + 1) * cb        // tmp, sm, one
        + (size_t)(6 * n + 16) * 64;
}
