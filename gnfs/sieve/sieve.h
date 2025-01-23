/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: sieve.h 638 2011-09-11 15:31:19Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _GNFS_SIEVE_SIEVE_H_
#define _GNFS_SIEVE_SIEVE_H_

#include <ms_batch_factor.h>
#include "gnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* factors smaller than the following are not printed */

#define MAX_SKIPPED_FACTOR 256

/* dump one NFS relation to the savefile */

void print_relation(savefile_t *savefile, int64 a, uint32 b, 
		uint32 *factors_r, uint32 num_factors_r, 
		uint32 large_prime_r[MAX_LARGE_PRIMES],
		uint32 *factors_a, uint32 num_factors_a, 
		uint32 large_prime_a[MAX_LARGE_PRIMES]);
	
/* convert k bits to another base */

uint32 fplog(uint32 k, double log_of_base);

/* find the size of f(a,b) both in bits and in base log_base */

int32 fplog_eval_poly(int64 a, uint32 b, mpz_t scratch,
			mpz_poly_t *f, double log_base,
			uint32 *bits);

/* compute a base of logarithms suitable for the current
   sieve line. The base chosen should be larger when the
   size of sieve values are smaller, in order to use up more
   of the dynamic range in one byte of the sieve array */

double get_log_base(mpz_poly_t *poly, 
			int64 a0, int64 a1, uint32 b);

uint32 read_last_line(msieve_obj *obj, mpz_t n);

void write_last_line(msieve_obj *obj, mpz_t n, uint32 b);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_SIEVE_SIEVE_H_ */
