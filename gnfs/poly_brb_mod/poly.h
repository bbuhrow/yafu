/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: poly.h 1025 2018-08-19 02:20:28Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _GNFS_POLY_POLY_H_
#define _GNFS_POLY_POLY_H_

#include <ms_common.h>
#include <dd.h>
#include <ddcomplex.h>
#include <integrate.h>
#include <polyroot.h>
#include "gnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

#if MAX_POLY_DEGREE < 6
#error "Polynomial generation assumes degree <= 6 allowed"
#endif

/* parameters */
typedef struct {
	double digits;
	double stage1_norm;
	double stage2_norm;
	double final_norm;
	uint32 deadline;
} poly_param_t;

/* used if polynomials will ever be generated in parallel */
#define POLY_HEAP_SIZE 1

/* when analyzing a polynomial's root properties, the
   bound on factor base primes that are checked */
#define PRIME_BOUND 2000

typedef struct {
	mpz_poly_t rpoly;
	mpz_poly_t apoly;
	double size_score;
	double root_score;
	double combined_score;
	double skewness;
	uint32 num_real_roots;
} poly_select_t;

/* main structure for poly selection */

typedef struct {
	poly_select_t *heap[POLY_HEAP_SIZE];
	uint32 heap_num_filled;

	integrate_t integ_aux;
	dickman_t dickman_aux;
} poly_config_t;

void poly_config_init(poly_config_t *config);
void poly_config_free(poly_config_t *config);

#define SIZE_EPS 1e-6

/* main routines */

void get_poly_params(msieve_obj *obj, mpz_t n,
			uint32 *degree_out, 
			poly_param_t *params_out);

void find_poly_core(msieve_obj *obj, mpz_t n,
			poly_param_t *params,
			poly_config_t *config,
			uint32 degree);

typedef struct {
	uint32 degree;
	double coeff[MAX_POLY_DEGREE + 1];
} dpoly_t;

typedef struct {
	uint32 degree;
	dd_t coeff[MAX_POLY_DEGREE + 1];
} ddpoly_t;

uint32 analyze_poly_size(integrate_t *integ_aux,
			ddpoly_t *rpoly, ddpoly_t *apoly, 
			double *result);

uint32 analyze_poly_murphy(integrate_t *integ_aux, dickman_t *dickman_aux,
			ddpoly_t *rpoly, double root_score_r,
			ddpoly_t *apoly, double root_score_a,
			double skewness, double *result,
			uint32 *num_real_roots);

uint32 analyze_poly_roots(mpz_poly_t *poly, uint32 prime_bound,
				double *result);

uint32 analyze_poly_roots_projective(mpz_poly_t *poly, 
				uint32 prime_bound,
				double *result);

void get_poly_combined_score(poly_select_t *poly);

void analyze_poly(poly_config_t *config, poly_select_t *poly);

void save_poly(poly_config_t *config, poly_select_t *poly);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_POLY_POLY_H_ */
