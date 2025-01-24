/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: common.h 1025 2018-08-19 02:20:28Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _COMMON_H_
#define _COMMON_H_

#include <ms_msieve.h>
#include <ms_gmp_xface.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_MPI
	#define MAX_MPI_GRID_DIM 35

	#define MPI_TRY(x) \
	{								\
		int status;						\
		if ((status = (x)) != MPI_SUCCESS) {			\
			printf("MPI error at %s:%d\n", __FILE__, __LINE__);\
			MPI_Abort(MPI_COMM_WORLD, status);		\
		}							\
	}

	#ifndef MPI_ERR_ASSERT  /* MPI library may not have this */
	#define MPI_ERR_ASSERT 22
	#endif
#endif


/*--------------PRIME SIEVE RELATED DECLARATIONS ---------------------*/

/* many separate places in the code need a list
   of primes. Include such a list pregenerated, with
   the differences between primes stored */

#define PRECOMPUTED_PRIME_BOUND 100000
#define PRECOMPUTED_NUM_PRIMES 9592
extern const uint8 prime_delta[PRECOMPUTED_NUM_PRIMES];

typedef struct {
	uint32 p;
	uint32 r;
} prime_aux_t;

typedef struct {
	uint32 num_aux;
	prime_aux_t *aux;
	uint8 *sieve;
	uint32 curr_block;
	uint32 curr_off;
} prime_sieve_t;

void init_prime_sieve(prime_sieve_t *sieve, 
			uint32 min_prime,
			uint32 max_prime);

void free_prime_sieve(prime_sieve_t *sieve);

uint32 get_next_prime(prime_sieve_t *sieve);

typedef struct {
	uint32 *list;
	uint32 num_primes;
} prime_list_t;

void fill_prime_list(prime_list_t *prime_list, 
			uint32 max_size, 
			uint32 max_prime);

/*--------------DECLARATIONS FOR MANAGING FACTORS FOUND -----------------*/

typedef struct {
	mp_t factor;
	enum msieve_factor_type type;
} final_factor_t;

typedef struct {
	uint32 num_factors;
	final_factor_t *final_factors[256];
} factor_list_t;

void factor_list_init(factor_list_t *list);

uint32 factor_list_max_composite(factor_list_t *list);

void factor_list_free(mp_t *n, factor_list_t *list, msieve_obj *obj);

uint32 factor_list_add(msieve_obj *obj, 
			factor_list_t *list, 
			mp_t *new_factor);

/*----------------- DECLARATIONS FOR HASHTABLES -----------------*/

/* we use two separate hash functions for multi-word structures */

#define HASH1(word) ((uint32)(word) * (uint32)(2654435761UL))
#define HASH2(word) ((uint32)(word) * ((uint32)40499 * 65543))

/* structure used by hashtables. Each entry indexed by the hashtable
   has blob_words words, plus a 'next' pointer. Only the first 
   hash_words words are inserted into the table and used for lookups */

typedef struct {
	uint32 blob_words;
	uint32 hash_words;

	uint32 *hashtable;
	uint32 log2_hashtable_size;
	uint32 num_used;
	uint32 congestion_target;

	uint32 *match_array;
	uint32 match_array_size;
	uint32 match_array_alloc;
} hashtable_t;

/* hash_words = 0 means hash_words = blob_words */

void hashtable_init(hashtable_t *h, 
		    uint32 blob_words, uint32 hash_words);

void hashtable_close(hashtable_t *h);
void hashtable_free(hashtable_t *h);
size_t hashtable_sizeof(hashtable_t *h);

static INLINE void hashtable_reset(hashtable_t *h) {
	h->match_array_size = 1;
	h->num_used = 0;

	memset(h->hashtable, 0, 
		sizeof(uint32) << h->log2_hashtable_size);
}

static INLINE uint32 hashtable_get_num(hashtable_t *h) {
	return h->match_array_size - 1;
}

static INLINE void *hashtable_get_first(hashtable_t *h) {
	return h->match_array + h->blob_words + 1;
}

static INLINE void *hashtable_get_next(hashtable_t *h, void *prev) {
	return (uint32 *)prev + h->blob_words + 1;
}

static INLINE uint32 hash_function(uint32 *data, uint32 num_words) {

	uint32 hashval = 0;

	switch(num_words) {
	default:
	case 4:
		hashval ^= HASH2(data[3] ^ 0x0f0f0f0f);
	case 3:
		hashval ^= HASH1(data[2] ^ 0x0f0f0f0f);
	case 2:
		hashval ^= HASH2(data[1]);
	case 1:
		hashval ^= HASH1(data[0]);
	case 0:
		break;
	}

	return hashval;
}

/* return a pointer to the hashtable entry that matches
   blob[]. If there is no such entry, make *present zero 
   (if the pointer is non-NULL) and add blob[] to the
   hashtable. Return a pointer to the entry in the hashtable
   that matches blob[]. blob[] is assigned a counter value
   that starts from zero, and the counter for the entry
   matching blob[] is output in *ordinal_id (if non-NULL) */


void *hashtable_find(hashtable_t *h, void *blob, 
		uint32 *ordinal_id, uint32 *present);

/*--------------DECLARATIONS FOR FACTORING METHODS -----------------*/

/* perform trial factoring on an integer.
   Returns 1 if any factors were found, 0 if not */

uint32 trial_factor(msieve_obj *obj, mp_t *n, mp_t *reduced_n,
		    factor_list_t *factor_list);

/* Attempt to factor an integer using Pollard's Rho method.
   Returns 1 if any factors were found, 0 if not */

uint32 rho(msieve_obj *obj, mp_t *n, mp_t *reduced_n, 
		factor_list_t *factor_list);

/* Attempt to factor an integer using the P-1, P+1 and
   ECM implementations of the GMP_ECM library
   Returns 1 if any factors were found, 0 if not */

uint32 ecm_pp1_pm1(msieve_obj *obj, mp_t *n, mp_t *reduced_n, 
		   factor_list_t *factor_list);

/* number of bits in the input that causes a switchover
   to the rigorous tiny factoring methods below */

#define SMALL_COMPOSITE_CUTOFF_BITS 85

/* Factor a number up to 62 bits in size using the SQUFOF
   algorithm. Returns zero if the factorization failed for 
   whatever reason, otherwise returns one factor up to 31 bits.
   Note that the factor returned may be 1, indicating a
   trivial factorization you probably don't want */

uint32 squfof(mp_t *n);

/* Factor a number up to 85 bits in size using MPQS. 
   Returns 0 on failure and nonzero on success, with
   factor1 and factor2 filled in on success. NOT TO
   BE USED for high-throughput applications */

uint32 tinyqs(mp_t *n, mp_t *factor1, mp_t *factor2);

/* Factor a number using the full MPQS implementation. 
   Returns 1 if any factors were found and 0 if not */

uint32 factor_mpqs(msieve_obj *obj, mp_t *n, factor_list_t *factor_list);

/* Factor a number using GNFS. Returns
   1 if any factors were found and 0 if not */

uint32 factor_gnfs(msieve_obj *obj, mp_t *n, factor_list_t *factor_list);

#define MIN_NFS_BITS 264

/*--------------LINEAR ALGEBRA RELATED DECLARATIONS ---------------------*/

/* used whenever temporary arrays are needed to store
   lists of factors of relations */

#define TEMP_FACTOR_LIST_SIZE 100

/* Used to represent a list of relations */

typedef struct {
	uint32 num_relations;  /* number of relations in the cycle */
	uint32 *list;          /* list of offsets into an array of relations */
} la_cycle_t;

/* A column of the matrix */

#define MAX_COL_IDEALS 1000

typedef struct {
	uint32 *data;		/* The list of occupied rows in this column */
	uint32 weight;		/* Number of nonzero entries in this column */
	la_cycle_t cycle;       /* list of relations comprising this column */
} la_col_t;

/* merge src1[] and src2[] into merge_array[], assumed
   large enough to hold the merged result. Return the
   final number of elements in merge_array */

uint32 merge_relations(uint32 *merge_array,
		  uint32 *src1, uint32 n1,
		  uint32 *src2, uint32 n2);

uint64 * block_lanczos(msieve_obj *obj,
			uint32 nrows, uint32 max_nrows, uint32 start_row,
			uint32 num_dense_rows,
			uint32 ncols, uint32 max_ncols, uint32 start_col,
			la_col_t *cols, uint32 *deps_found);

uint64 count_matrix_nonzero(msieve_obj *obj,
			uint32 nrows, uint32 num_dense_rows,
			uint32 ncols, la_col_t *cols);

uint64 reduce_matrix(msieve_obj *obj, uint32 *nrows, 
		uint32 num_dense_rows, uint32 *ncols, 
		la_col_t *cols, uint32 num_excess);

#define MIN_REORDER_SIZE 200000

void reorder_matrix(msieve_obj *obj, uint32 **rowperm, uint32 **colperm);

void free_cycle_list(la_col_t *cycle_list, uint32 num_cycles);

void dump_cycles(msieve_obj *obj, la_col_t *cols, uint32 ncols);

void dump_matrix(msieve_obj *obj, 
		uint32 nrows, uint32 num_dense_rows,
		uint32 ncols, la_col_t *cols,
		uint64 num_nonzero);

void read_matrix(msieve_obj *obj, 
		uint32 *nrows, uint32 *max_nrows, uint32 *start_row,
		uint32 *num_dense_rows_out,
		uint32 *ncols, uint32 *max_ncols, uint32 *start_col,
		la_col_t **cols_out,
		uint32 *rowperm, uint32 *colperm);

void dump_dependencies(msieve_obj *obj, 
			uint64 *deps, uint32 ncols);

void read_cycles(msieve_obj *obj, 
		uint32 *num_cycles_out, 
		la_col_t **cycle_list_out, 
		uint32 dependency,
		uint32 *colperm);

/*-------------- MISCELLANEOUS STUFF ----------------------------------*/

#define POSITIVE 0
#define NEGATIVE 1

/* emit logging information */

void logprintf(msieve_obj *obj, char *fmt, ...);

/* convert an expression into an mp_t; returns 0 on
   success, negative value on failure */

int32 evaluate_expression(char *expr, mp_t *res);

/* remember a factor that was just found */

void add_next_factor(msieve_obj *obj, mp_t *n, 
		enum msieve_factor_type factor_type);

/* definitions for multidimensional minimization */

#define MAX_VARS 5

typedef double (*objective_func)(double v[MAX_VARS], void *extra);

double minimize(double p[MAX_VARS], uint32 ndim, 
			double ftol, uint32 max_iter,
			objective_func callback, void *extra);

double minimize_global(double p[MAX_VARS], uint32 ndim,
                        double limits[MAX_VARS][2],
			double tol, uint32 iter_limit,
			objective_func callback, 
			void *extra);

typedef double (*objective_func_grad)(double v[MAX_VARS], 
					double grad[MAX_VARS],
					void *extra);

double minimize_grad(double p[MAX_VARS], uint32 ndim, 
			double ftol, uint32 max_iter,
			objective_func callback, 
			objective_func_grad callback_grad,
			void *extra);

typedef double (*objective_func_hess)(double v[MAX_VARS], 
					double grad[MAX_VARS],
					double hess[MAX_VARS][MAX_VARS],
					void *extra);

double minimize_hess(double p[MAX_VARS], uint32 ndim, 
			double ftol, uint32 max_iter,
			objective_func callback, 
			objective_func_hess callback_hess,
			void *extra);

/* solve a linear system of size n; matrix and b are overwritten */

void solve_dmatrix(double matrix[MAX_VARS][MAX_VARS], 
			double x[MAX_VARS],
			double b[MAX_VARS],
			uint32 n);

/* Dickman's rho function */

#define DICKMAN_ACCURACY 1e-8

typedef struct {
	uint32 num_coeffs;
	uint32 coeff_offset;
} dickman_line_t;

typedef struct {
	dickman_line_t *lines;
	double *coeffs;
} dickman_t;

void dickman_init(dickman_t *aux);
void dickman_free(dickman_t *aux);
double dickman(dickman_t *aux, double arg);

#ifdef __cplusplus
}
#endif

#endif /* _COMMON_H_ */
