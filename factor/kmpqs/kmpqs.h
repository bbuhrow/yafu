/*
kmpqs: Kleinjung's small-number MPQS from lasieve5 (gnfs4linux), ported
into yafu.

The algorithms are Copyright (C) 2001 Jens Franke, T. Kleinjung, and this
module is distributed under the terms of the GNU General Public Licence and
WITHOUT ANY WARRANTY.  See factor/lasieve5_64/COPYING.
*/

#ifndef KMPQS_H
#define KMPQS_H

#include <stddef.h>
#include "gmp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Largest input (in bits) each variant accepts. */
#define KMPQS1_MAX_BITS 96
#define KMPQS3_MAX_BITS 148

/* Returned when an internal consistency check fails. */
#define KMPQS_ERR_INTERNAL (-100)

/*
 * Factor N completely (all factors prime, or prime powers).
 *
 * kmpqs1 handles N up to 96 bits; kmpqs3 handles roughly 80..148 bits;
 * kmpqs_factor picks one by size (kmpqs1 up to 96 bits, as lasieve does).
 * N must be odd, composite, not a perfect power, and have no factors below
 * a few thousand (as in lasieve, callers trial-divide first).
 *
 * Returns the number of factors found (>= 2) and points *factors at an
 * array inside the context; it stays valid until the next call with the
 * same context.  Returns 0 if a factor exceeds max_bits, and a negative
 * value on failure.
 */

/*
 * Explicit contexts.  A context holds all working state (about 0.3 MB for
 * kmpqs1 and 1.3 MB for kmpqs3, allocated on first use) and may be used by
 * one thread at a time.  Reusing one avoids the setup cost per call.
 */
typedef struct kmpqs_ctx kmpqs_ctx_t;

kmpqs_ctx_t *kmpqs_ctx_new(void);
void kmpqs_ctx_free(kmpqs_ctx_t *ctx);

int kmpqs_factor_ctx(kmpqs_ctx_t *ctx, mpz_t N, size_t max_bits, mpz_t **factors);
int kmpqs1_factor_ctx(kmpqs_ctx_t *ctx, mpz_t N, size_t max_bits, mpz_t **factors);
int kmpqs3_factor_ctx(kmpqs_ctx_t *ctx, mpz_t N, size_t max_bits, mpz_t **factors);

/*
 * Convenience versions that use a context private to the calling thread,
 * created on first use.  kmpqs_thread_cleanup() frees it; call it before a
 * thread exits if the thread used these, or the context leaks.
 */
int kmpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors);
int kmpqs1_factor(mpz_t N, size_t max_bits, mpz_t **factors);
int kmpqs3_factor(mpz_t N, size_t max_bits, mpz_t **factors);
void kmpqs_thread_cleanup(void);

/* Nonzero: print lasieve's warnings and internal errors to stderr. */
extern int kmpqs_verbose;

#ifdef __cplusplus
}
#endif

#endif
