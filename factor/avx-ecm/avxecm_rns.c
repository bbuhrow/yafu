/*
Copyright (c) 2026, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * AVX-ECM arithmetic backends.
 *
 * AVX-ECM's curve code works on vec_bignum_t through the vec*mod_ptr
 * function pointers plus vecCopy / vecClear.  Two representations sit
 * behind them:
 *   radix 2^52 (mdata->rns == NULL): NWORDS digits per lane, Montgomery
 *     with R = 2^(52 NWORDS) (template / BFIPS / karatsuba cores);
 *   RNS (mdata->rns != NULL): the Yao batch layout of arith/rns
 *     (rnsyb_*: 2ns+1 words x 8 lanes, lane = curve), Montgomery with
 *     R = M_B.  The vectors are allocated with >= n*s words so vecCopy /
 *     vecClear cover the whole layout.
 *
 * Values cross between mpz and the vectors only through the helpers below
 * (curve setup, the stage-2 batch inversion, gcds), which do the right
 * thing for either representation.  RNS outputs are lazily reduced
 * (< k*N, k = n + 2); the curve code's add/sub chains stay within the RNS
 * bounds (mul inputs < 8kN, subtrahends < 2kN), see
 * arith/handoff_vector_arith.md.
 */

#include "avx_ecm.h"
#include "../../arith/rns/rns_yao.h"
#include <stdlib.h>
#include <string.h>

typedef struct
{
    rnsy_fixed_t *fx;       // size-class tables, shared by all threads
    int own_fx;             // this mdata frees fx
    rnsy_ctx_t *ctx;        // N-dependent constants + kernel scratch (per thread)
    uint64_t *x1;           // one number, single-number layout
    uint64_t *vt;           // one batch vector (add/sub temporary)
} avxecm_rns_t;

#define RNSB(md) ((avxecm_rns_t *)(md)->rns)

// ---------------------------------------------------------------------------
// backend-neutral helpers
// ---------------------------------------------------------------------------

void vec_lane_get(vec_monty_t *mdata, mpz_t dest, vec_bignum_t *v, int lane)
{
#ifdef RNSYB_HAVE_VEC
    if (mdata->rns != NULL)
    {
        avxecm_rns_t *r = RNSB(mdata);
        rnsyb_vget_lane(r->ctx, r->x1, v->data, lane);
        rnsy_to_int(r->ctx, dest, r->x1);
        return;
    }
#endif
    extract_bignum_from_vec_to_mpz(dest, v, lane, (int)mdata->NWORDS);
}

void vec_lane_set(vec_monty_t *mdata, vec_bignum_t *v, mpz_t src, int lane)
{
#ifdef RNSYB_HAVE_VEC
    if (mdata->rns != NULL)
    {
        avxecm_rns_t *r = RNSB(mdata);
        mpz_t t;
        mpz_init(t);
        mpz_mod(t, src, r->ctx->N);     // < N < M_B
        rnsy_from_int(r->ctx, r->x1, t);
        rnsyb_vset_lane(r->ctx, v->data, lane, r->x1);
        mpz_clear(t);
        return;
    }
#endif
    insert_mpz_to_vec(v, src, lane);
}

void vec_set_raw_ui(vec_monty_t *mdata, vec_bignum_t *v, uint64_t x)
{
    int j;
#ifdef RNSYB_HAVE_VEC
    if (mdata->rns != NULL)
    {
        avxecm_rns_t *r = RNSB(mdata);
        mpz_t t;
        mpz_init_set_ui(t, (unsigned long)x);
        rnsy_from_int(r->ctx, r->x1, t);
        for (j = 0; j < VECLEN_I; j++)
            rnsyb_vset_lane(r->ctx, v->data, j, r->x1);
        mpz_clear(t);
        return;
    }
#endif
    vecClear(v);
    for (j = 0; j < VECLEN_I; j++)
        v->data[j] = x;
    v->size = 1;
}

void vec_monty_mulR(vec_monty_t *mdata, mpz_t x, mpz_t n)
{
#ifdef RNSYB_HAVE_VEC
    if (mdata->rns != NULL)
    {
        mpz_mul(x, x, RNSB(mdata)->ctx->MBmodN);
        mpz_mod(x, x, n);
        return;
    }
#endif
    mpz_mul_2exp(x, x, DIGITBITS * mdata->NWORDS);
    mpz_tdiv_r(x, x, n);
}

void vec_set_monty_one(vec_monty_t *mdata, vec_bignum_t *v, mpz_t n)
{
    mpz_t t;
    mpz_init_set_ui(t, 1);
    if (mdata->isMersenne == 0)
        vec_monty_mulR(mdata, t, n);
#ifdef RNSYB_HAVE_VEC
    if (mdata->rns != NULL)
    {
        int j;
        for (j = 0; j < VECLEN_I; j++)
            vec_lane_set(mdata, v, t, j);
        mpz_clear(t);
        return;
    }
#endif
    broadcast_mpz_to_vec(v, t);
    mpz_clear(t);
}

#ifdef RNSYB_HAVE_VEC

// ---------------------------------------------------------------------------
// setup
// ---------------------------------------------------------------------------

int avxecm_rns_available(void)
{
    return 1;
}

void *avxecm_rns_fixed(mpz_t N, uint32_t *words)
{
    int bits = (int)mpz_sizeinbase(N, 2);
    int n = rnsyb_choose_n(bits);
    rnsy_fixed_t *fx;
    rnsy_ctx_t *ctx;

    if (n <= 0)
        return NULL;
    fx = rnsy_fixed_init_n(bits, n);
    if (fx == NULL)
        return NULL;
    // N must be coprime to the moduli (otherwise the gcd is a factor that
    // trial division should have found; just fall back)
    ctx = rnsy_ctx_init(fx, N);
    if (ctx == NULL)
    {
        rnsy_fixed_free(fx);
        return NULL;
    }
    rnsy_ctx_free(ctx);
    *words = (uint32_t)(fx->n * fx->s);
    return fx;
}

void avxecm_rns_fixed_info(void *fxv, int *n, int *s)
{
    rnsy_fixed_t *fx = (rnsy_fixed_t *)fxv;
    *n = fx->n;
    *s = fx->s;
}

int avxecm_rns_attach(vec_monty_t *mdata, void *fxv, int own_fx, mpz_t N)
{
    rnsy_fixed_t *fx = (rnsy_fixed_t *)fxv;
    avxecm_rns_t *r;
    rnsy_ctx_t *ctx = rnsy_ctx_init(fx, N);

    if (ctx == NULL)
        return 0;
    if ((uint32_t)fx->nwords > 2 * mdata->NWORDS + 1)
    {
        printf("avxecm_rns_attach: vectors too small for the RNS layout\n");
        exit(1);
    }
    r = (avxecm_rns_t *)malloc(sizeof(avxecm_rns_t));
    r->fx = fx;
    r->own_fx = own_fx;
    r->ctx = ctx;
    r->x1 = rnsy_alloc(ctx);
    r->vt = rnsyb_alloc(ctx);
    mdata->rns = r;
    return 1;
}

void avxecm_rns_detach(vec_monty_t *mdata)
{
    avxecm_rns_t *r = RNSB(mdata);
    if (r == NULL)
        return;
    rns_vfree(r->x1);
    rns_vfree(r->vt);
    rnsy_ctx_free(r->ctx);
    if (r->own_fx)
        rnsy_fixed_free(r->fx);
    free(r);
    mdata->rns = NULL;
}

void *avxecm_rns_get_fixed(vec_monty_t *mdata)
{
    return (mdata->rns != NULL) ? RNSB(mdata)->fx : NULL;
}

void avxecm_rns_free_fixed(void *fxv)
{
    rnsy_fixed_free((rnsy_fixed_t *)fxv);
}

// ---------------------------------------------------------------------------
// the vec*mod_ptr API on RNS batch vectors (n and s arguments are unused)
// ---------------------------------------------------------------------------

void avxecm_rns_mul(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    rnsyb_vmulmod(RNSB(mdata)->ctx, c->data, a->data, b->data);
}

void avxecm_rns_sqr(vec_bignum_t *a, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata)
{
    rnsyb_vsqrmod(RNSB(mdata)->ctx, c->data, a->data);
}

void avxecm_rns_add(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,
    vec_monty_t *mdata)
{
    rnsyb_vaddmod(RNSB(mdata)->ctx, c->data, a->data, b->data);
}

void avxecm_rns_sub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c,
    vec_monty_t *mdata)
{
    rnsyb_vsubmod(RNSB(mdata)->ctx, c->data, a->data, b->data);
}

// sum = a + b, diff = a - b; sum and diff may alias a and/or b
void avxecm_rns_addsub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *sum,
    vec_bignum_t *diff, vec_monty_t *mdata)
{
    avxecm_rns_t *r = RNSB(mdata);
    if ((sum != a) && (sum != b))
    {
        rnsyb_vaddmod(r->ctx, sum->data, a->data, b->data);
        rnsyb_vsubmod(r->ctx, diff->data, a->data, b->data);
    }
    else if ((diff != a) && (diff != b))
    {
        rnsyb_vsubmod(r->ctx, diff->data, a->data, b->data);
        rnsyb_vaddmod(r->ctx, sum->data, a->data, b->data);
    }
    else
    {
        rnsyb_vsubmod(r->ctx, r->vt, a->data, b->data);
        rnsyb_vaddmod(r->ctx, sum->data, a->data, b->data);
        memcpy(diff->data, r->vt, (size_t)r->fx->nwords * VECLEN * sizeof(base_t));
    }
}

#else

int avxecm_rns_available(void) { return 0; }
void *avxecm_rns_fixed(mpz_t N, uint32_t *words) { return NULL; }
void avxecm_rns_fixed_info(void *fxv, int *n, int *s) { *n = *s = 0; }
int avxecm_rns_attach(vec_monty_t *mdata, void *fxv, int own_fx, mpz_t N) { return 0; }
void avxecm_rns_detach(vec_monty_t *mdata) { }
void *avxecm_rns_get_fixed(vec_monty_t *mdata) { return NULL; }
void avxecm_rns_free_fixed(void *fxv) { }

#endif
