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

/*
 * Factor N completely (all factors prime, or prime powers).
 *
 * kmpqs1_factor handles N up to 96 bits; kmpqs3_factor handles roughly
 * 80..148 bits.  N must be odd, composite, not a perfect power, and have no
 * factors below a few thousand (as in lasieve, callers trial-divide first).
 *
 * Returns the number of factors found (>= 2) and points *factors at an
 * internal array holding them; the array is overwritten by the next call
 * on the same thread.  Returns 0 if a factor exceeds max_bits, and a
 * negative value on failure.
 */
int kmpqs1_factor(mpz_t N, size_t max_bits, mpz_t **factors);
int kmpqs3_factor(mpz_t N, size_t max_bits, mpz_t **factors);

/* Picks kmpqs1 or kmpqs3 by the size of N. */
int kmpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors);

/* Nonzero: print lasieve's warnings and internal errors to stderr. */
extern int kmpqs_verbose;

#ifdef __cplusplus
}
#endif

#endif
