/*
kmpqs: Kleinjung's small-number MPQS from lasieve5 (gnfs4linux), ported
into yafu.

The algorithms are Copyright (C) 2001 Jens Franke, T. Kleinjung, and this
module is distributed under the terms of the GNU General Public Licence and
WITHOUT ANY WARRANTY.  See factor/lasieve5_64/COPYING.
*/

#include "kmpqs.h"

int kmpqs_verbose = 0;

int kmpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
	/* lasieve's strategy: mpqs up to 96 bits, mpqs3 above */
	if (mpz_sizeinbase(N, 2) > KMPQS1_MAX_BITS)
		return kmpqs3_factor(N, max_bits, factors);
	return kmpqs1_factor(N, max_bits, factors);
}
