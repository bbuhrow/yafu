/*
 * vecpoly polynomial algorithms (see vecpoly.h and
 * arith/avxecm_stage2_fft_design.md section 4).
 *
 * Memory: one NTT buffer of K * Lmax vectors (w->ntt) plus a small one for
 * the second operand of a full-length product, VP_GRP primes at a time
 * (w->grp).  Products of length L <= Lmax/2 keep both operands in w->ntt.
 *
 * Scratch use: w->tmp[0] (Lmax/2 coefficients) belongs to vp_mulmod_monic
 * and vp_recip; w->tmp[1] (Lmax/4) to vp_mul_monic and the evaluation
 * nodes; w->sm to the schoolbook paths.  None of these call each other
 * with live scratch.
 */

#include "vecpoly.h"
#include <string.h>
#include <stdlib.h>

#ifdef VP_HAVE_KERNELS

extern void *vp_aligned_alloc(size_t bytes);
extern void vp_aligned_free(void *p);

#define CO(a, j) VP_COEF(w->ctx, (a), (j))

static void vp_wrap(vec_bignum_t *v, const uint64_t *p, int n)
{
    v->data = (base_t *)p;
    v->size = n;
    v->WORDS_ALLOC = (uint32_t)n;
    v->signmask = 0;
}

static void vp_copy(vp_work_t *w, uint64_t *dst, const uint64_t *src, uint32_t cnt)
{
    if (dst != src)
        memmove(dst, src, (size_t)cnt * 8 * (size_t)w->ctx->n * sizeof(uint64_t));
}

static void vp_zero(vp_work_t *w, uint64_t *dst, uint32_t cnt)
{
    memset(dst, 0, (size_t)cnt * 8 * (size_t)w->ctx->n * sizeof(uint64_t));
}

// reverse a[0..cnt) in place
static void vp_reverse(vp_work_t *w, uint64_t *a, uint32_t cnt)
{
    uint64_t *t = w->sm;
    uint32_t i;
    for (i = 0; i < cnt / 2; i++)
    {
        vp_copy(w, t, CO(a, i), 1);
        vp_copy(w, CO(a, i), CO(a, cnt - 1 - i), 1);
        vp_copy(w, CO(a, cnt - 1 - i), t, 1);
    }
}

void vp_cmul(vp_work_t *w, uint64_t *c, const uint64_t *a, const uint64_t *b,
    uint32_t cnt)
{
    vec_bignum_t A, B, C;
    uint32_t t;
    int n = w->ctx->n;
    for (t = 0; t < cnt; t++)
    {
        vp_wrap(&A, CO(a, t), n);
        vp_wrap(&B, CO(b, t), n);
        vp_wrap(&C, CO(c, t), n);
        w->mul(&A, &B, &C, &w->vn, &w->vs, &w->md);
    }
}

// schoolbook: c[0..la+lb-1) = a * b.  c must not overlap a or b; uses the
// last coefficient of w->sm.
static void vp_school(vp_work_t *w, uint64_t *c, const uint64_t *a, uint32_t la,
    const uint64_t *b, uint32_t lb)
{
    uint64_t *t = CO(w->sm, 4 * VP_SCHOOL_CAP - 1);
    uint32_t i, j;
    vp_zero(w, c, la + lb - 1);
    for (i = 0; i < la; i++)
    {
        for (j = 0; j < lb; j++)
        {
            vp_cmul(w, t, CO(a, i), CO(b, j), 1);
            vpk_add(w->ctx, CO(c, i + j), CO(c, i + j), t, 1);
        }
    }
}

// schoolbook is used for min(la, lb) <= school_max, and when the product
// must go through w->sm, also la + lb <= 3 VP_SCHOOL_CAP
static int vp_use_school(vp_work_t *w, uint32_t la, uint32_t lb, int via_sm)
{
    uint32_t mn = la < lb ? la : lb;
    if ((int)mn > w->ctx->school_max)
        return 0;
    return !via_sm || (la + lb <= 3 * VP_SCHOOL_CAP);
}

static uint32_t vp_pow2(uint32_t x)
{
    uint32_t L = 1;
    while (L < x) L <<= 1;
    return L;
}

static void vp_check_len(vp_ctx_t *ctx, uint32_t L)
{
    if (L > ctx->Lmax)
    {
        fprintf(stderr, "vecpoly: transform length %u > Lmax %u\n", L, ctx->Lmax);
        exit(1);
    }
}

// w->ntt (length L, all primes, forward transformed) *= transform of
// b[0..lb), computed VP_GRP primes at a time in w->grp
static void vp_ntt_mul_by(vp_work_t *w, uint32_t L, const uint64_t *b, uint32_t lb)
{
    vp_ctx_t *ctx = w->ctx;
    int i0, i1;
    for (i0 = 0; i0 < ctx->K; i0 = i1)
    {
        uint64_t *v = w->ntt + (size_t)i0 * L * 8;
        i1 = i0 + VP_GRP < ctx->K ? i0 + VP_GRP : ctx->K;
        vpk_to_ntt_r(ctx, w->grp, L, b, lb, i0, i1);
        vpk_fwd_r(ctx, w->grp, L, i0, i1);
        vpk_pmul_r(ctx, v, v, w->grp, L, i0, i1);
    }
}

// w->ntt = transform of a * b, length L (no inverse transform)
static void vp_ntt_prod(vp_work_t *w, const uint64_t *a, uint32_t la,
    const uint64_t *b, uint32_t lb, uint32_t L)
{
    vp_ctx_t *ctx = w->ctx;
    vp_check_len(ctx, L);
    vpk_to_ntt(ctx, w->ntt, L, a, la);
    vpk_fwd(ctx, w->ntt, L);
    if (2 * L <= ctx->Lmax)
    {
        uint64_t *v2 = w->ntt + (size_t)ctx->K * L * 8;
        vpk_to_ntt(ctx, v2, L, b, lb);
        vpk_fwd(ctx, v2, L);
        vpk_pmul(ctx, w->ntt, w->ntt, v2, L);
    }
    else
        vp_ntt_mul_by(w, L, b, lb);
}

// cyclic convolution of length L of a[0..la) and b[0..lb), coefficients
// start..start+count-1 into c (c may overlap a or b)
static void vp_cyclic(vp_work_t *w, uint64_t *c, const uint64_t *a, uint32_t la,
    const uint64_t *b, uint32_t lb, uint32_t L, uint32_t start, uint32_t count)
{
    vp_ntt_prod(w, a, la, b, lb, L);
    vpk_inv(w->ctx, w->ntt, L);
    vpk_from_ntt(w->ctx, c, w->ntt, L, start, count);
}

// coefficients start..start+count-1 of a * b into c (c may overlap a or b),
// the product being exact in a cyclic convolution of length L
static void vp_mul_part(vp_work_t *w, uint64_t *c, const uint64_t *a, uint32_t la,
    const uint64_t *b, uint32_t lb, uint32_t L, uint32_t start, uint32_t count)
{
    if (vp_use_school(w, la, lb, 1))
    {
        vp_school(w, w->sm, a, la, b, lb);
        vp_copy(w, c, CO(w->sm, start), count);
    }
    else
        vp_cyclic(w, c, a, la, b, lb, L, start, count);
}

void vp_mul(vp_work_t *w, uint64_t *c, const uint64_t *a, uint32_t la,
    const uint64_t *b, uint32_t lb)
{
    if (vp_use_school(w, la, lb, 0))
        vp_school(w, c, a, la, b, lb);
    else
        vp_cyclic(w, c, a, la, b, lb, vp_pow2(la + lb - 1), 0, la + lb - 1);
}

void vp_mul_monic(vp_work_t *w, uint64_t *c, const uint64_t *a,
    const uint64_t *b, uint32_t d)
{
    // (X^d + a)(X^d + b) = X^2d + X^d (a + b) + a b
    uint64_t *s = w->tmp[1];
    vpk_add(w->ctx, s, a, b, d);
    vp_mul_part(w, c, a, d, b, d, 2 * d, 0, 2 * d - 1);
    vp_zero(w, CO(c, 2 * d - 1), 1);
    vpk_add(w->ctx, CO(c, d), CO(c, d), s, d);
}

void vp_from_roots(vp_work_t *w, uint64_t *out, uint64_t **tree,
    const uint64_t *roots, uint32_t m)
{
    uint32_t d, i;
    int s;

    if (m == 1)
    {
        vpk_neg(w->ctx, out, roots, 1);     // no tree levels
        return;
    }
    if (tree != NULL)
    {
        // level s lives in tree[s], or in out when tree[s] is NULL
        uint64_t *cur = tree[0];
        vpk_neg(w->ctx, cur, roots, m);
        for (d = 1, s = 0; 2 * d < m; d <<= 1, s++)
        {
            uint64_t *dst = tree[s + 1] != NULL ? tree[s + 1] : out;
            for (i = 0; i < m; i += 2 * d)
                vp_mul_monic(w, CO(dst, i), CO(cur, i), CO(cur, i + d), d);
            cur = dst;
        }
        vp_mul_monic(w, out, cur, CO(cur, d), d);   // in place when cur == out
        return;
    }

    vpk_neg(w->ctx, out, roots, m);
    for (d = 1; d < m; d <<= 1)
        for (i = 0; i < m; i += 2 * d)
            vp_mul_monic(w, CO(out, i), CO(out, i), CO(out, i + d), d);
}

void vp_recip(vp_work_t *w, uint64_t *inv, const uint64_t *f, uint32_t m,
    uint64_t *scratch)
{
    // S = 1 / rev(F) mod X^m by Newton, built in inv; then inv = rev(S)
    uint64_t *A = scratch, *S = inv;
    uint64_t *P = w->tmp[0], *T = CO(w->tmp[0], m / 2);
    uint32_t k, j;

    vp_copy(w, A, w->one, 1);
    for (j = 1; j < m; j++)
        vp_copy(w, CO(A, j), CO(f, m - j), 1);
    vp_copy(w, S, w->one, 1);

    for (k = 1; k < m; k <<= 1)
    {
        // e' = (A mod X^2k) * S_k, coefficients k..2k-1 (exact in a
        // length-2k cyclic product: the wrapped part lands below k)
        vp_mul_part(w, P, A, 2 * k, S, k, 2 * k, k, k);
        // T = S_k * e' mod X^k
        vp_mul_part(w, T, S, k, P, k, 2 * k, 0, k);
        vpk_neg(w->ctx, CO(S, k), T, k);
    }
    vp_reverse(w, inv, m);
}

size_t vp_modf_bytes(int n, int K, uint32_t m, int flags)
{
    size_t b = (size_t)m * 64 * (size_t)n + sizeof(vp_modf_t);
    if (flags & VP_MODF_KEEP_INV) b += (size_t)K * 2 * m * 64;
    if (flags & VP_MODF_KEEP_F) b += (size_t)K * m * 64;
    return b;
}

vp_modf_t *vp_modf_init(vp_work_t *w, uint64_t *f, uint32_t m, int flags,
    uint64_t *scratch)
{
    vp_ctx_t *ctx = w->ctx;
    vp_modf_t *mf = (vp_modf_t *)calloc(1, sizeof(vp_modf_t));
    uint64_t *sc = scratch != NULL ? scratch : vp_alloc(m, ctx->n);

    vp_check_len(ctx, 2 * m);
    mf->m = m;
    mf->f = f;
    mf->inv = vp_alloc(m, ctx->n);
    vp_recip(w, mf->inv, f, m, sc);
    if (scratch == NULL)
        vp_free(sc);
    if (flags & VP_MODF_KEEP_INV)
    {
        mf->ntt_inv = (uint64_t *)vp_aligned_alloc((size_t)ctx->K * 2 * m * 64);
        vpk_to_ntt(ctx, mf->ntt_inv, 2 * m, mf->inv, m);
        vpk_fwd(ctx, mf->ntt_inv, 2 * m);
    }
    if (flags & VP_MODF_KEEP_F)
    {
        mf->ntt_f = (uint64_t *)vp_aligned_alloc((size_t)ctx->K * m * 64);
        vpk_to_ntt(ctx, mf->ntt_f, m, f, m);
        vpk_fwd(ctx, mf->ntt_f, m);
    }
    return mf;
}

void vp_modf_free(vp_modf_t *mf)
{
    if (mf == NULL) return;
    vp_free(mf->inv);
    if (mf->ntt_inv != NULL) vp_aligned_free(mf->ntt_inv);
    if (mf->ntt_f != NULL) vp_aligned_free(mf->ntt_f);
    free(mf);
}

void vp_mod_monic(vp_work_t *w, vp_modf_t *mf, uint64_t *h, const uint64_t *g)
{
    vpk_sub(w->ctx, h, g, mf->f, mf->m);
}

// w->ntt (length 2m, forward) *= transform of inv
static void vp_ntt_mul_inv(vp_work_t *w, vp_modf_t *mf)
{
    if (mf->ntt_inv != NULL)
        vpk_pmul(w->ctx, w->ntt, w->ntt, mf->ntt_inv, 2 * mf->m);
    else
        vp_ntt_mul_by(w, 2 * mf->m, mf->inv, mf->m);
}

void vp_mulmod_monic(vp_work_t *w, vp_modf_t *mf, uint64_t *h, const uint64_t *g)
{
    vp_ctx_t *ctx = w->ctx;
    uint32_t m = mf->m, L = 2 * m;
    uint64_t *B = w->tmp[0];

    // P = h g + X^m h = Q F + R.  Q = (P_hi * inv)[m-1 .. 2m-2]; then, mod
    // X^m - 1: R = (P_lo + P_hi) - Q - (Q f mod X^m - 1), a length-m product
    vp_ntt_prod(w, h, m, g, m, L);
    vpk_inv(ctx, w->ntt, L);
    vpk_from_ntt(ctx, B, w->ntt, L, m, m - 1);
    vp_zero(w, CO(B, m - 1), 1);
    vpk_add(ctx, B, B, h, m);                           // B = P_hi
    vpk_from_ntt(ctx, h, w->ntt, L, 0, m);              // h = P_lo

    vpk_to_ntt(ctx, w->ntt, L, B, m);
    vpk_fwd(ctx, w->ntt, L);
    vpk_add(ctx, h, h, B, m);                           // h = P mod X^m - 1
    vp_ntt_mul_inv(w, mf);
    vpk_inv(ctx, w->ntt, L);
    vpk_from_ntt(ctx, B, w->ntt, L, m - 1, m);          // B = Q
    vpk_sub(ctx, h, h, B, m);

    vpk_to_ntt(ctx, w->ntt, m, B, m);
    vpk_fwd(ctx, w->ntt, m);
    if (mf->ntt_f != NULL)
        vpk_pmul(ctx, w->ntt, w->ntt, mf->ntt_f, m);
    else
        vp_ntt_mul_by(w, m, mf->f, m);
    vpk_inv(ctx, w->ntt, m);
    vpk_from_ntt(ctx, B, w->ntt, m, 0, m);              // B = Q f mod X^m - 1
    vpk_sub(ctx, h, h, B, m);
}

// one node of the transposed evaluation, in place: u (d = 2hh values)
// becomes uA (u[0..hh)) and uB (u[hh..2hh)) for the children a and b
// (monic, hh coefficients each):
//   uA[s] = u[s+hh] + sum_r b_r u[s+r];  uB likewise with a
static void vp_eval_node(vp_work_t *w, uint64_t *u, const uint64_t *a,
    const uint64_t *b, uint32_t hh)
{
    vp_ctx_t *ctx = w->ctx;
    uint32_t d = 2 * hh, i, r;

    if (((int)hh <= ctx->school_max) && (hh <= VP_SCHOOL_CAP))
    {
        uint64_t *uA = w->sm, *uB = CO(w->sm, hh), *t = CO(w->sm, 2 * hh);
        for (i = 0; i < hh; i++)
        {
            vp_copy(w, CO(uA, i), CO(u, i + hh), 1);
            vp_copy(w, CO(uB, i), CO(u, i + hh), 1);
            for (r = 0; r < hh; r++)
            {
                vp_cmul(w, t, CO(b, r), CO(u, i + r), 1);
                vpk_add(ctx, CO(uA, i), CO(uA, i), t, 1);
                vp_cmul(w, t, CO(a, r), CO(u, i + r), 1);
                vpk_add(ctx, CO(uB, i), CO(uB, i), t, 1);
            }
        }
        vp_copy(w, u, uA, d);
        return;
    }

    {
        // transform of u in the first half of w->ntt, products in the second
        uint64_t *rv = w->tmp[1];
        uint64_t *v1 = w->ntt, *v2 = w->ntt + (size_t)ctx->K * d * 8;
        vpk_to_ntt(ctx, v1, d, u, d);
        vpk_fwd(ctx, v1, d);

        for (i = 0; i < hh; i++)
            vp_copy(w, CO(rv, i), CO(b, hh - 1 - i), 1);
        vpk_to_ntt(ctx, v2, d, rv, hh);
        vpk_fwd(ctx, v2, d);
        vpk_pmul(ctx, v2, v2, v1, d);
        vpk_inv(ctx, v2, d);
        vpk_from_ntt(ctx, u, v2, d, hh - 1, hh);        // u[0..hh) is consumed
        vpk_add(ctx, u, u, CO(u, hh), hh);              // uA

        for (i = 0; i < hh; i++)
            vp_copy(w, CO(rv, i), CO(a, hh - 1 - i), 1);
        vpk_to_ntt(ctx, v2, d, rv, hh);
        vpk_fwd(ctx, v2, d);
        vpk_pmul(ctx, v2, v2, v1, d);
        vpk_inv(ctx, v2, d);
        vpk_from_ntt(ctx, rv, v2, d, hh - 1, hh);
        vpk_add(ctx, CO(u, hh), CO(u, hh), rv, hh);     // uB
    }
}

void vp_eval_tree(vp_work_t *w, vp_modf_t *mf, uint64_t **tree,
    uint64_t *vals, const uint64_t *h, uint64_t *scratch)
{
    vp_ctx_t *ctx = w->ctx;
    uint32_t m = mf->m, i, hh;
    int s;

    // root: U[j] = (h * inv)[2m-2-j], j = 0..m-1
    vp_copy(w, vals, h, m);
    if (m == 1)
    {
        vp_cmul(w, vals, vals, mf->inv, 1);
        return;
    }
    if (vp_use_school(w, m, m, 1))
    {
        vp_school(w, w->sm, vals, m, mf->inv, m);
        vp_copy(w, vals, CO(w->sm, m - 1), m);
    }
    else
    {
        vpk_to_ntt(ctx, w->ntt, 2 * m, vals, m);
        vpk_fwd(ctx, w->ntt, 2 * m);
        vp_ntt_mul_inv(w, mf);
        vpk_inv(ctx, w->ntt, 2 * m);
        vpk_from_ntt(ctx, vals, w->ntt, 2 * m, m - 1, m);
    }
    vp_reverse(w, vals, m);

    s = 0;
    while ((1u << s) < m) s++;
    for (s = s - 1; s >= 0; s--)
    {
        const uint64_t *lev = tree[s];
        hh = 1u << s;
        if (lev == NULL)
        {
            // rebuild level s from the nearest kept level below (tree[0]
            // is always kept), in place in scratch after the first step
            const uint64_t *lo;
            int b = s - 1, q;
            if (scratch == NULL)
            {
                fprintf(stderr, "vecpoly: vp_eval_tree needs scratch\n");
                exit(1);
            }
            while (tree[b] == NULL) b--;
            lo = tree[b];
            for (q = b + 1; q <= s; q++)
            {
                uint32_t dd = 1u << (q - 1);
                for (i = 0; i < m; i += 2 * dd)
                    vp_mul_monic(w, CO(scratch, i), CO(lo, i), CO(lo, i + dd), dd);
                lo = scratch;
            }
            lev = scratch;
        }
        for (i = 0; i < m; i += 2 * hh)
            vp_eval_node(w, CO(vals, i), CO(lev, i), CO(lev, i + hh), hh);
    }
}

void vp_prod(vp_work_t *w, uint64_t *prod, const uint64_t *a, uint32_t cnt)
{
    uint32_t t;
    vp_copy(w, prod, a, 1);
    for (t = 1; t < cnt; t++)
        vp_cmul(w, prod, prod, CO(a, t), 1);
}

#endif
