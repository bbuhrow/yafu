/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: stage2.h 732 2012-08-04 02:32:46Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _STAGE2_H_
#define _STAGE2_H_

#include <poly_skew.h>

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------*/
/* data used in the current polynomial */

typedef struct {
	mpz_t gmp_a[MAX_POLY_DEGREE + 1];
	mpz_t gmp_b[MAX_POLY_DEGREE + 1];
	mpz_t gmp_c[MAX_POLY_DEGREE + 1];
	mpz_t gmp_lina[2];
	mpz_t gmp_linb[2];
	mpz_t gmp_linc[2];
	mpz_t gmp_help1;
	mpz_t gmp_help2;
	mpz_t gmp_help3;
	mpz_t gmp_help4;
	mpz_t gmp_p;
	mpz_t gmp_d;
} curr_poly_t;

/*-----------------------------------------------------------------------*/
/* data for rating polynomial yield */

typedef struct {
	integrate_t integ_aux;
	dickman_t dickman_aux;
} assess_t;

void assess_init(assess_t *a);

void assess_free(assess_t *a);

uint32 stage2_root_score(uint32 deg1, mpz_t *coeff1, 
       			uint32 prime_bound, double *score,
       			uint32 projective_only);

/*-----------------------------------------------------------------------*/
/* routines for optimizing polynomials */

void optimize_initial(curr_poly_t *data, uint32 deg, double *pol_norm,
			uint32 skew_only);

double optimize_initial_deg6(double best[MAX_VARS], 
			curr_poly_t *c, uint32 degree);

void optimize_final(mpz_t x, mpz_t y, int64 z, poly_rootopt_t *data);

double optimize_basic(dpoly_t *apoly, double *best_skewness,
				double *best_translation);

/*-----------------------------------------------------------------------*/
/* data for the root sieve */

#define MAX_SIEVE_PRIME 100
#define ROOT_HEAP_SIZE 200
#define LOG_SCALE_FACTOR 1000
#define ROOT_SCORE_COARSE_MIN (-4.0)

typedef struct {
	uint16 resclass;
	uint16 start;
	uint16 step;
} sieve_root_t;

typedef struct {
	uint32 power;
	uint32 num_roots;
	double root_contrib;
	uint16 sieve_contrib;
	sieve_root_t *roots;
} sieve_power_t;

typedef struct {
	uint32 prime;
	uint32 num_powers;
	sieve_power_t *powers;

	uint32 contrib_array_size;
	uint32 contrib_array_offset;
	uint16 *contrib_array;
} sieve_prime_t;

typedef struct {
	mpz_t x;
	mpz_t y;
	int64 z;
	double score;
} rotation_t;

typedef struct {
	uint32 num_entries;
	uint32 max_entries;
	rotation_t *entries;
	void *extra;
} root_heap_t;

/* definitions for root sieve */

#define MAX_CRT_FACTORS 10

typedef struct {
	uint16 score;
	uint64 x;
	uint64 y;
	uint64 z;
} lattice_t;

#define MAX_ROOTS 64

typedef struct {
	uint8 power;
	uint8 num_roots;
	uint16 score[MAX_ROOTS];
	uint8 roots[MAX_ROOTS][3];
} hit_t;


void compute_lattices(hit_t *hitlist, uint32 num_lattice_primes,
			lattice_t *lattices, uint64 lattice_size,
			uint32 num_lattices, uint32 dim);
void compute_line_size(double max_norm, dpoly_t *apoly,
		  double dbl_p, double dbl_d, double direction[3],
		  double last_line_min_in, double last_line_max_in,
		  double *line_min, double *line_max);


/* sieve for selecting XYZ triplets */

typedef struct {
	uint64 lattice_size;
	sieve_prime_t lattice_primes[MAX_CRT_FACTORS];
	uint32 num_lattice_primes;

	uint32 num_lattices;
	lattice_t *lattices;

	int64 z_base;
	uint32 z_blocks;
	double *y_line_min;
	double *y_line_max;
} sieve_xyz_t;

void sieve_xyz_alloc(sieve_xyz_t *xyz);
void sieve_xyz_free(sieve_xyz_t *xyz);

/* root sieve for selecting XY planes */

typedef struct {
	uint64 lattice_size;
	sieve_prime_t lattice_primes[MAX_CRT_FACTORS];
	uint32 num_lattice_primes;

	mpz_t y_base;
	uint32 y_blocks;

	mpz_t mp_lattice_size;
	double dbl_lattice_size;

	/* degree 6 only */
	uint16 curr_score;
	dpoly_t apoly;
	mpz_t crt0;
	mpz_t crt1;
	mpz_t resclass_x;
	mpz_t resclass_y;
	mpz_t tmp1, tmp2, tmp3, tmp4;

	/* degree 4 and 5 only */
	uint32 num_lattices;
	lattice_t *lattices;
	double *x_line_min;
	double *x_line_max;

} sieve_xy_t;

void sieve_xy_alloc(sieve_xy_t *xy);
void sieve_xy_free(sieve_xy_t *xy);


/* root sieve for selecting X lines */

typedef struct {
	uint64 lattice_size;
	sieve_prime_t lattice_primes[MAX_CRT_FACTORS];
	uint32 num_lattice_primes;

	dpoly_t apoly;

	double last_line_min;
	double last_line_max;
	uint16 curr_score;

	mpz_t x_base;
	uint32 x_blocks;

	mpz_t mp_lattice_size;
	double dbl_lattice_size;
	mpz_t crt0;
	mpz_t crt1;
	mpz_t resclass;
	mpz_t tmp1;
} sieve_x_t;

void sieve_x_alloc(sieve_x_t *x);
void sieve_x_free(sieve_x_t *x);

#define UNROLL 4
#define DEFAULT_BLOCK_SIZE  8192

typedef struct {
	poly_rootopt_t *data;

	uint32 num_primes;
	sieve_prime_t *primes;

	double sieve_bias; 
	double random_root_score;

	uint16 *sieve_block;

	dpoly_t apoly;
	double dbl_p;
	double dbl_d;
	double max_norm;

	root_heap_t root_heap;

	sieve_xyz_t xyzdata;
	sieve_xy_t xydata;
	sieve_x_t xdata;

	mpz_t curr_x, curr_y;
	int64 curr_z;
} root_sieve_t;

void root_sieve_init(root_sieve_t *rs);
void root_sieve_free(root_sieve_t *rs);
void root_sieve_run(poly_rootopt_t *data, double curr_norm,
				double alpha_proj);

uint64 find_lattice_size_z(double line_length);
void sieve_xyz_run_deg6(root_sieve_t *rs, uint64 lattice_size,
			double line_min, double line_max);

uint64 find_lattice_size_y(double line_length);
void sieve_xy_run_deg5(root_sieve_t *rs, uint64 lattice_size,
			double line_min, double line_max);
void sieve_xy_run_deg6(root_sieve_t *rs);

uint64 find_lattice_size_x(mpz_t prev_lattice_size,
				double line_length);
void sieve_x_run_deg4(root_sieve_t *rs, uint64 lattice_size,
			double line_min, double line_max);
void sieve_x_run_deg5(root_sieve_t *rs);
void sieve_x_run_deg6(root_sieve_t *rs);

void root_sieve_line(root_sieve_t *rs);

void save_rotation(root_heap_t *heap, mpz_t x, mpz_t y,
		int64 z, float score);

/*-------------------------------------------------------------------------*/

/* data for optimizing a single (ad, p, d) triplet */

typedef struct {
	curr_poly_t curr_poly;
	root_sieve_t root_sieve;
	assess_t assess;
	double size_cutoff;
} stage2_curr_data_t;

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE2_H_ */
