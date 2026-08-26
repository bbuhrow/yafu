/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage1.h 1023 2018-08-19 00:30:42Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _STAGE1_H_
#define _STAGE1_H_

/* Interface for selecting GNFS polynomials whose top
   three coefficients are small. We use Kleinjung's improved
   algorithm presented at the 2008 CADO Factoring Workshop,
   with many modifications */

#include <poly_skew.h>
#include <cuda_xface.h>
#include <thread.h>


#ifdef __cplusplus
extern "C" {
#endif

#define MAX_POLYSELECT_DEGREE 6

#if MAX_POLY_DEGREE < MAX_POLYSELECT_DEGREE
#error "supported poly degree must be at least 6"
#endif

/* we try to limit the high-order algebraic poly coefficient 
   to the set with many small primes, to increase the odds 
   that the resulting polynomial will have unusually many 
   projective roots modulo small primes. */

#define HIGH_COEFF_PRIME_LIMIT 100

/* for practical reasons, we limit the size of product of
   small primes that we look for in the high-order algebraic
   coefficient. High coeffs larger than
   HIGH_COEFF_MULTIPLIER * HIGH_COEFF_SIEVE_LIMIT
   will have additional factors that may be large */

#define HIGH_COEFF_SIEVE_LIMIT 1e12

/*-----------------------------------------------------------------------*/

/* main structure used in stage 1 */

typedef struct {

	uint32 degree;

	mpz_t N;

	/* bound on the norm used in stage 1; this is the maximum
	   value of (poly coefficient i) * (optimal skew)^i across
	   all poly coefficients. The low-order poly coefficients
	   are bounded in size by m0, because the polynomial is
	   essentially N split into d+1 pieces. Our job is to find
	   a 'stage 1 hit' that obeys the norm bound even for the
	   high-order algebraic coefficients */

	double norm_max;

	/* the range on a_d, provided by calling code */

	mpz_t gmp_high_coeff_begin;
	mpz_t gmp_high_coeff_end;

	/* user override for leading coeff increment (0=auto) */

	uint32 high_coeff_multiplier;

	/* if nonzero, read leading coefficients from coeff_list.txt */

	uint32 use_coeff_list;

	/* if nonzero, stop the search once this many stage-1 polynomials
	   have been found (the num_polys= flag). poly_count is the running
	   total. On the GPU path each worker owns its own lock-free count in
	   device_thread_data_t.polys_found; when a leading coefficient
	   completes the worker publishes the sum here for the main thread to
	   read. That publish is a plain cross-thread store of a word-sized
	   advisory value (only ever advanced) rather than a synchronized
	   counter, which is fine because reaching the target only stops
	   starting new coefficients and is already approximate. The CPU path
	   is single-threaded and increments poly_count directly. poly_count
	   is only read, never used to drive obj->flags, so reaching the
	   target stops stage 1 without aborting the rest of the
	   factorization. */

	uint32 target_poly_count;
	uint32 poly_count;

	/* function to call when a collision is found */

	stage1_callback_t callback;
	void *callback_data;

	/* internal stuff */

	mpz_t tmp1;

} poly_search_t;


/* data for searching a single leading coefficient */

typedef struct {

	uint32 degree;

	mpz_t high_coeff; 

	/* (computed) bound on the third-highest algebraic
	   poly coefficient. Making this small is the only
	   thing stage 1 can do; the other coefficients can
	   only be optimized in stage 2 */

	double coeff_max;

	/* approximately (N/a_d) ^ (1/d) for input N, high
	   coefficient a_d and poly degree d */

	double m0;

	/* bound on the leading rational poly coefficient */

	double p_size_max;
	
	/* effective norm_max used for this coefficient (may be
	   dynamically computed to maximize special_q range) */

	double norm_max_effective;

	/* internal values used */

	mpz_t trans_N;
	mpz_t trans_m0;
	mpz_t m; 
	mpz_t p;
	mpz_t tmp1;
	mpz_t tmp2;
	mpz_t tmp3;
	
	/* number of stage-1 polynomials found for this leading coefficient;
	   reset when the coefficient starts and reported when it finishes */

	uint32 found_count;
} poly_coeff_t;

poly_coeff_t * poly_coeff_init(void);
void poly_coeff_free(poly_coeff_t *c);
void poly_coeff_copy(poly_coeff_t *dest, poly_coeff_t *src);

/*-----------------------------------------------------------------------*/

/* Kleinjung's algorithm essentially reduces to an
   all-against-all search between two large collections of
   arithmetic progressions, looking for pairs of progressions
   that satisfy a specific modular property. The following
   controls how much larger the largest element in a collection
   is, compared to the smallest */

#define P_SCALE 2.4

/* Rational leading coeffs of NFS polynomials are assumed 
   to be the product of 2 or 3 coprime groups of factors p; 
   each p is < 2^32 and the product of (powers of) one or more
   distinct primes. The arithmetic progressions mentioned 
   above are of the form r_i + k * p^2 for 'root' r_i. 
   p may be prime or composite.
   
   We get a collision when we can find two progressions 
   aprog1(k) = r_i1 + k*p1^2 and aprog2(k) = r_i2+k*p2^2, 
   and integers k1 and k2, such that

   	- p1 and p2 are coprime
	- aprog1(k1) = aprog2(k2)
	- the common value is "close to" m0
   
   If p has several prime factors p_i, the exact number of 
   roots that a given p has is the product of the number of 
   d_th roots of N modulo each p_i. For degree 4 polyomials, N 
   has either 0, 2, or 4 fourth roots mod p_i. For degree 5, 
   N has either 1 or 5 fifth roots mod p_i. For degree 6, 
   N has either 0, 2 or 6 sixth roots mod p_i. So especially
   when p is large and p_i are small, a composite p can 
   contribute a large number of progressions to the collision 
   search

   We will need to generate many p along with all their r_i
   fairly often, and need efficient methods to do so */

/* initialize and free the factory */

void * sieve_fb_alloc(void);

void sieve_fb_free(void *s_in);

/* initialize the factory for a new leading algebraic 
   coefficient. p is allowed to have small prime
   factors p_i between factor_min and factor_max, and N
   will have between fb_roots_min and fb_roots_max d_th
   roots modulo each of these p_i. Additionally, if fb_only
   is zero, a sieve is used to find any prime p which are
   larger than factor_max */

void sieve_fb_init(void *s_in, poly_coeff_t *coeff,
			uint32 factor_min, uint32 factor_max,
			uint32 fb_roots_min, uint32 fb_roots_max,
			uint32 fb_only);

/* set up for a run of p production. The generated p will 
   all be between p_min and p_max, and the number of roots
   for each p is between num_roots_min and num_roots_max 
   (bounded by MAX_ROOTS) */

#define MAX_ROOTS 128

void sieve_fb_reset(void *s_in, uint32 p_min, uint32 p_max,
			uint32 num_roots_min, uint32 num_roots_max);

/* count the (p, root) pairs the factory will produce for the
   given parameters, without computing any roots. Only cheap
   for factories restricted to smooth p (fb_only nonzero);
   returns 0 (unknown) if large prime p would also be produced.
   The factory is left reset and ready for a fresh run */

uint64 sieve_fb_count(void* s_in, uint32 p_min, uint32 p_max,
	uint32 num_roots_min, uint32 num_roots_max);

/* function that 'does something' when a single p 
   and all its roots is found */

typedef void (*root_callback)(uint32 p, uint32 num_roots, uint64 *roots, 
				void *extra);

/* find the next p and all of its roots. The code returns
   P_SEARCH_DONE if no more p exist, otherwise it calls 
   callback() and returns p.

   p which are products of small primes are found first, then
   large prime p (since prime p are slower and you may not want
   all of them). Prime p will be found in ascending order but 
   no order may be assumed for composite p returned by 
   consecutive calls */

#define P_SEARCH_DONE ((uint32)(-2))

uint32 sieve_fb_next(void *s_in, poly_coeff_t *c, 
			root_callback callback,
			void *extra);

/*-----------------------------------------------------------------------*/

/* what to do when the collision search finds a 'stage 1 hit' */

uint32
handle_collision(poly_coeff_t *c, uint64 p, uint32 special_q,
		uint64 special_q_root, int64 res);

/* main search routine */

#ifdef HAVE_CUDA_POLY

/* GPU search routine */

double sieve_lattice_gpu(msieve_obj *obj,
			poly_search_t *poly, 
			poly_coeff_t *c, 
			void *gpu_data,
			double deadline);

void * gpu_data_init(msieve_obj *obj, poly_search_t *poly);
void gpu_data_free(void *gpu_data);

#else

/* CPU search routine */

double sieve_lattice_cpu(msieve_obj *obj,
			poly_search_t *poly, 
			poly_coeff_t *c,
			double deadline);

#endif

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_H_ */
