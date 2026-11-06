/*
kmpqs: Kleinjung's small-number MPQS from lasieve5 (gnfs4linux), ported
into yafu.

The algorithms are Copyright (C) 2001 Jens Franke, T. Kleinjung, and this
module is distributed under the terms of the GNU General Public Licence and
WITHOUT ANY WARRANTY.  See factor/lasieve5_64/COPYING.
*/

#include <stdlib.h>
#include <stdint.h>
#include "kmpqs.h"

/* implemented in kmpqs1.c / kmpqs3.c; ctx is mpqs_ctx_t / mpqs3_ctx_t */
void *kmpqs1_ctx_new(void);
void kmpqs1_ctx_free(void *ctx);
int kmpqs1_run(void *ctx, mpz_t N, size_t max_bits, mpz_t **factors);
void *kmpqs3_ctx_new(void);
void kmpqs3_ctx_free(void *ctx);
int kmpqs3_run(void *ctx, mpz_t N, size_t max_bits, mpz_t **factors);

#if defined(_MSC_VER)
#define KMPQS_TLS __declspec(thread)
#else
#define KMPQS_TLS __thread
#endif

int kmpqs_verbose = 0;

struct kmpqs_ctx {
	void *c1;	/* mpqs_ctx_t, created on first kmpqs1 use */
	void *c3;	/* mpqs3_ctx_t, created on first kmpqs3 use */
};

kmpqs_ctx_t *kmpqs_ctx_new(void)
{
	return (kmpqs_ctx_t *)calloc(1, sizeof(kmpqs_ctx_t));
}

void kmpqs_ctx_free(kmpqs_ctx_t *ctx)
{
	if (ctx == NULL)
		return;
	kmpqs1_ctx_free(ctx->c1);
	kmpqs3_ctx_free(ctx->c3);
	free(ctx);
}

int kmpqs1_factor_ctx(kmpqs_ctx_t *ctx, mpz_t N, size_t max_bits, mpz_t **factors)
{
	if (ctx->c1 == NULL && (ctx->c1 = kmpqs1_ctx_new()) == NULL)
		return KMPQS_ERR_INTERNAL;
	return kmpqs1_run(ctx->c1, N, max_bits, factors);
}

int kmpqs3_factor_ctx(kmpqs_ctx_t *ctx, mpz_t N, size_t max_bits, mpz_t **factors)
{
	if (ctx->c3 == NULL && (ctx->c3 = kmpqs3_ctx_new()) == NULL)
		return KMPQS_ERR_INTERNAL;
	return kmpqs3_run(ctx->c3, N, max_bits, factors);
}

int kmpqs_factor_ctx(kmpqs_ctx_t *ctx, mpz_t N, size_t max_bits, mpz_t **factors)
{
	/* lasieve's strategy: mpqs up to 96 bits, mpqs3 above */
	if (mpz_sizeinbase(N, 2) > KMPQS1_MAX_BITS)
		return kmpqs3_factor_ctx(ctx, N, max_bits, factors);
	return kmpqs1_factor_ctx(ctx, N, max_bits, factors);
}

/* ---- thread-local convenience API ---- */

static KMPQS_TLS kmpqs_ctx_t *kmpqs_tls_ctx = NULL;

static kmpqs_ctx_t *tls_ctx(void)
{
	if (kmpqs_tls_ctx == NULL)
		kmpqs_tls_ctx = kmpqs_ctx_new();
	return kmpqs_tls_ctx;
}

int kmpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
	kmpqs_ctx_t *ctx = tls_ctx();

	return ctx ? kmpqs_factor_ctx(ctx, N, max_bits, factors) : KMPQS_ERR_INTERNAL;
}

int kmpqs1_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
	kmpqs_ctx_t *ctx = tls_ctx();

	return ctx ? kmpqs1_factor_ctx(ctx, N, max_bits, factors) : KMPQS_ERR_INTERNAL;
}

int kmpqs3_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
	kmpqs_ctx_t *ctx = tls_ctx();

	return ctx ? kmpqs3_factor_ctx(ctx, N, max_bits, factors) : KMPQS_ERR_INTERNAL;
}

void kmpqs_thread_cleanup(void)
{
	kmpqs_ctx_free(kmpqs_tls_ctx);
	kmpqs_tls_ctx = NULL;
}

/* ---- parameter tables ---- */

int kmpqs1_param(int row, uint16_t *get, const uint16_t *set);
int kmpqs3_param(int row, uint16_t *get, const uint16_t *set);
int kmpqs1_row(void *ctx, mpz_t N);
int kmpqs3_row(void *ctx, mpz_t N);

int kmpqs_param(int variant, int row, unsigned short *get, const unsigned short *set)
{
	if (variant == 1) return kmpqs1_param(row, get, set);
	if (variant == 3) return kmpqs3_param(row, get, set);
	return -1;
}

int kmpqs_param_row(int variant, mpz_t N)
{
	kmpqs_ctx_t *ctx = tls_ctx();

	if (ctx == NULL)
		return -1;
	if (variant == 1) {
		if (ctx->c1 == NULL && (ctx->c1 = kmpqs1_ctx_new()) == NULL) return -1;
		return kmpqs1_row(ctx->c1, N);
	}
	if (variant == 3) {
		if (ctx->c3 == NULL && (ctx->c3 = kmpqs3_ctx_new()) == NULL) return -1;
		return kmpqs3_row(ctx->c3, N);
	}
	return -1;
}
