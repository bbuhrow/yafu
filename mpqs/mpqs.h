/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: mpqs.h 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _MPQS_MPQS_H_
#define _MPQS_MPQS_H_

/* An implementation of the Self-Initializing Multiple 
   Polynomial Quadratic Sieve algorithm for integer 
   factorization. */

#include <ms_common.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------- SIEVE RELATED DECLARATIONS ---------------------*/
/* There's a limit to how many factors can contribute to
   SIQS polynomials */

#define MAX_POLY_FACTORS 20

/* The size of a factor base prime is limited too */

#define MAX_FB_PRIME 0x3ffffff

/*  representing one prime in the factor base */

#define INVALID_ROOT ((uint32)(-1))

typedef struct {
	uint32 prime : 26;	/* the factor base prime */
	uint32 logprime : 5;	/* log2(prime)rounded to nearest */
	uint32 rcorrect : 1;    /* used to correct multiples by recip */
	uint32 root1;	        /* the two square roots of an MPQS */
	uint32 root2;	        /* polynomial (INVALID_ROOT means don't use) */
	uint32 recip;           /* integer reciprocal of prime */
} fb_t;

/* A packed representation of the above structure, used
   during the sieving phase only. Careful: the room allowed
   here for the factor base prime is less than is allowed
   in the ordinary factor base. That's okay, since this
   structure is only used for factor base primes less than
   the sieve block size */

typedef struct {
	uint32 prime : 24;	/* the factor base prime */
	uint8 logprime;		/* log2(prime)rounded to nearest */
	uint32 next_loc1;	/* the two next sieve locations for */
	uint32 next_loc2;	/*	the above two roots */
} packed_fb_t;

/* Offset 0 in the array of factor base primes is reserved
   for the "prime" -1. This allows the sieve to be over both
   positive and negative values. When actually dividing by
   factor base primes, start with the fb_t at this position in
   the factor base array: */

#define MIN_FB_OFFSET 1

/* structure representing one MPQS polynomial. Every
   relation found by the sieve must be associated with
   (i.e. must refer to) one of these. For sieve offset
   x, the MPQS polynomial is a * x + b */

typedef struct poly_t {
	uint32 a_idx;		/* offset into a list of 'a' values */
	signed_mp_t b;			/* the MPQS 'b' value */
} poly_t;

/* If a sieve value turns out to be smooth, the following is used
   to record information for use in later stages of the algorithm.
   One sieve value is packed into a relation_t, and an la_col_t
   contains a collection of sieve values */

#define POSITIVE 0
#define NEGATIVE 1

typedef struct relation_t {
	uint32 sieve_offset;	/* value of x in a*x+b */
	uint32 poly_idx;	/* pointer to this relation's MPQS poly */
	uint32 num_factors;	/* size of list of factors */
	uint32 large_prime[2];	/* for partial relations, the leftover 
					factor(s). Either factor may equal 1 */
	uint32 *fb_offsets;	/* The array of offsets into the factor base
					which contain primes dividing this
					sieve value. Duplicate offsets are
					possible. Sorted in ascending order */
} relation_t;

/* Once the sieving stage has collected enough relations,
   we'll need a few extra relations to make sure that the
   linear algebra stage finds linear dependencies */

#define NUM_EXTRA_RELATIONS 64

/* The sieving phase uses a hashtable to handle the large
   factor base primes. The hashtable is an array of type
   bucket_t, each entry of which contains an array 
   of type bucket_entry_t */

typedef struct {
	uint32 logprime : 8;       /* the logprime field from the factor base */
	uint32 prime_index : 24;   /* factor base array offset for entry */
	uint32 sieve_offset;       /* offset within one hash bin */
} bucket_entry_t;

typedef struct {
	uint32 num_alloc;     /* total list size this hashtable position */
	uint32 num_used;      /* number occupied at this hashtable position */
	bucket_entry_t *list; /* list of entries at this hashtable position */
} bucket_t;

/* Configuration of the sieving code requires passing the
   following structure */

typedef struct {
	uint32 bits;       /* size of integer this config info applies to */
	uint32 fb_size;    /* number of factor base primes */
	uint32 large_mult; /* the large prime multiplier */
	uint32 sieve_size; /* the size of the sieve (actual sieve is 2x this) */
} sieve_param_t;

/* The sieving code needs a special structure for determining
   the number of cycles in a collection of partial relations */

#define LOG2_CYCLE_HASH 22

typedef struct {
	uint32 next;
	uint32 prime;
	uint32 data;
	uint32 count;
} cycle_t;

/* To avoid huge parameter lists passed between sieving
   routines, all of the relevant data used in the sieving
   phase is packed into a single structure. Routines take
   the information they need out of this. */

typedef struct {
	msieve_obj *obj;     /* object controlling entire factorization */
	mp_t *n;             /* the number to factor (scaled by multiplier)*/
	uint32 multiplier;   /* small multiplier for n (may be composite) */

	/* bookkeeping information for sieving */

	uint32 sieve_block_size; /* bytes in one sieve block */
	uint8 *sieve_array;  /* scratch space used for one sieve block */
	fb_t *factor_base;       /* the factor base to use */
	packed_fb_t *packed_fb;  /* scratch space for a packed version of it */
	uint32 fb_size;          /* number of factor base primes (includes -1)*/

	uint32 num_sieve_blocks; /* number of sieve blocks in sieving interval*/
	uint32 poly_block;       /* number of FB primes in a block of work */
	uint32 fb_block;     /* number of polynomials simultaneously sieved */

	uint32 sieve_small_fb_start;  /* starting FB offset for sieving */
	uint32 sieve_large_fb_start;

	uint32 tf_small_recip1_cutoff;
	uint32 tf_small_recip2_cutoff;
	uint32 tf_med_recip1_cutoff;
	uint32 tf_med_recip2_cutoff;
	uint32 tf_large_cutoff;

	bucket_t *buckets;  /* hash bins for sieve values */
	uint32 cutoff1;          /* if log2(sieve value) exceeds this number,
				    a little trial division is performed */
	uint32 cutoff2;          /* if log2(sieve value) exceeds this number,
				    full trial division is performed */

	/* bookkeeping information for creating polynomials */

	uint32 *modsqrt_array;    /* a square root of n mod each FB prime */
	mp_t target_a;            /* optimal value of 'a' */
	uint32 a_bits;            /* required number of bits in the 'a' value 
	                                          of all MPQS polynomials */

	uint32 total_poly_a;    /* total number of polynomial 'a' values */
	mp_t *poly_a_list;        /* list of 'a' values for MPQS polys */
	poly_t *poly_list;   /* list of MPQS polynomials */

	mp_t curr_a;	  	  /* the current 'a' value */
	signed_mp_t *curr_b;    /* list of all the 'b' values for that 'a' */
	uint8 *next_poly_action;  /* see poly.c */

	uint32 num_poly_factors;  /* number of factors in poly 'a' value */
	uint32 factor_bounds[25]; /* bounds on FB offsets for primes that */
	                          /*        can be factors of MPQS polys */	
	uint32 poly_factors[MAX_POLY_FACTORS];  /* factorization of curr. 'a' */
	uint8 factor_bits[MAX_POLY_FACTORS]; /* size of each factor of 'a' */
	signed_mp_t poly_tmp_b[MAX_POLY_FACTORS];  /* temporary quantities */
	uint32 *poly_b_array;      /* precomputed values for all factor base
	                              primes, used to compute new polynomials */
	uint32 *poly_b_small[MAX_POLY_FACTORS];

	/* bookkeeping information for double large primes */

	uint32 large_prime_max;  /* the cutoff value for keeping a partial
				    relation; actual value, not a multiplier */
	mp_t max_fb2;          /* the square of the largest factor base prime */
	mp_t large_prime_max2; /* the cutoff value for factoring partials */

	relation_t *relation_list;     /* list of full/partial relations */
	uint32 num_relations;	/* number of relations in list */
	la_col_t *cycle_list;   /* cycles derived from relations */
	uint32 num_cycles;	/* number of cycles in list */

	cycle_t *cycle_table;      /* list of all the vertices in the graph */
	uint32 cycle_table_size;   /* number of vertices filled in the table */
	uint32 cycle_table_alloc;  /* number of cycle_t structures allocated */
	uint32 *cycle_hashtable;   /* hashtable to index into cycle_table */
	uint32 components;         /* connected components (see relation.c) */
	uint32 vertices;           /* vertices in graph (see relation.c) */

} sieve_conf_t;

/* attempt to trial factor one sieve value */

uint32 check_sieve_val(sieve_conf_t *conf,
		      int32 sieve_offset, 
		      uint32 bits,
		      mp_t *a,
		      signed_mp_t *b,
		      signed_mp_t *c,
		      uint32 poly_index,
		      bucket_t *hash_bucket);

/* the core of the sieving code */

typedef uint32 (*qs_core_sieve_fcn)(sieve_conf_t *conf,
				    uint32 poly_start,
				    uint32 num_poly);

#define DECLARE_SIEVE_FCN(name) 		\
	uint32 name(sieve_conf_t *conf, 	\
			uint32 poly_start,	\
			uint32 num_poly);

DECLARE_SIEVE_FCN(qs_core_sieve_generic_32k);
DECLARE_SIEVE_FCN(qs_core_sieve_generic_64k);

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	#define HAS_MSVC_SIEVE_CORE
	DECLARE_SIEVE_FCN(qs_core_sieve_vc8_32k);
	DECLARE_SIEVE_FCN(qs_core_sieve_vc8_64k);
#elif defined(GCC_ASM32X)
	DECLARE_SIEVE_FCN(qs_core_sieve_p2_64k);
	DECLARE_SIEVE_FCN(qs_core_sieve_p3_64k);
	DECLARE_SIEVE_FCN(qs_core_sieve_p4_64k);
	DECLARE_SIEVE_FCN(qs_core_sieve_pm_32k);
	DECLARE_SIEVE_FCN(qs_core_sieve_core_32k);
	DECLARE_SIEVE_FCN(qs_core_sieve_k7_64k);
	DECLARE_SIEVE_FCN(qs_core_sieve_k7xp_64k);
	DECLARE_SIEVE_FCN(qs_core_sieve_k8_64k);
#elif defined(GCC_ASM64X)
	DECLARE_SIEVE_FCN(qs_core_sieve_p4_64k);
	DECLARE_SIEVE_FCN(qs_core_sieve_core_32k);
	DECLARE_SIEVE_FCN(qs_core_sieve_k8_64k);
#endif

/* lightweight profiling of the sieve routines */
#ifdef SIEVE_TIMING
#define TIME1(var) { uint64 tmp_##var = read_clock();
#define TIME2(var) (var) += read_clock() - tmp_##var; }
extern uint64 next_poly_small_time;
extern uint64 next_poly_large_time;
extern uint64 plus_init_time;
extern uint64 bucket_time;
extern uint64 sieve_small_time;
extern uint64 sieve_large_time;
extern uint64 tf_plus_scan_time;
extern uint64 tf_total_time;
#else
#define TIME1(var) /* nothing */
#define TIME2(var) /* nothing */
#endif

/* pull out the large primes from a relation read from
   the savefile */

void read_large_primes(char *buf, uint32 *prime1, uint32 *prime2);

/* given the primes from a sieve relation, add
   that relation to the graph used for tracking
   cycles */

void add_to_cycles(sieve_conf_t *conf, 
			uint32 large_prime1, 
			uint32 large_prime2);

/* encapsulate all of the information conerning a sieve
   relation and dump it to the savefile */

void save_relation(sieve_conf_t *conf,
		  uint32 sieve_offset, 
		  uint32 *fb_offsets, 
		  uint32 num_factors, 
		  uint32 poly_index,
		  uint32 large_prime1,
		  uint32 large_prime2);

/* perform postprocessing on a list of relations */

void qs_filter_relations(sieve_conf_t *conf);

/* Initialize/free polynomial scratch data. Note that
   poly_free can only be called after the relation
   filtering phase */

void poly_init(sieve_conf_t *conf, uint32 sieve_size);
void poly_free(sieve_conf_t *conf);

/* compute a random polynomial 'a' value, and also
   compute all of the 'b' values, all of the precomputed
   quantities for the 'b' values, and all of the initial
   roots in the factor base */

void build_base_poly(sieve_conf_t *conf);

/* As above, except it's assumed that conf->poly_factors
   and conf->curr_a are already filled in. Factor base roots
   are not affected */

void build_derived_poly(sieve_conf_t *conf);

/* The main function to perform sieving.
	obj is the object controlling this factorization
	n is the number to be factored
	poly_list is a linked list of all the polynomials created
	factor_base contains the factor base (duh)
	modsqrt_array contains modular square roots of FB primes
	params is used to initialize and configure the sieve code
	multiplier is a small integer by which n is scaled
	relation_list is the list of relations from the sieving stage
	num_relations is the size of relation_list
	cycle_list is the list of cycles that the QS filtering code builds
	num_cycles is the number of cycles */

void do_sieving(msieve_obj *obj, mp_t *n, mp_t **poly_a_list, 
		poly_t **poly_list, fb_t *factor_base, 
		uint32 *modsqrt_array, sieve_param_t *params, 
		uint32 multiplier,
		relation_t **relation_list,
		uint32 *num_relations,
		la_col_t **cycle_list,
		uint32 *num_cycles);

void qs_free_relation_list(relation_t *list, uint32 num_relations);

/*--------------LINEAR ALGEBRA RELATED DECLARATIONS ---------------------*/

/* Find linear dependencies. The number of nontrivial dependencies
   found is returned
	obj is the object controlling this factorization
	fb_size is the size of the factor base
	bitfield is an array of num_cycles numbers. Bit i of word j tells
		whether cycle_list[j] is used in nullspace vector i. 
		Essentially, bitfield[] is a collection of 64 nullspace
		vectors packed together.
	relation_list is the list of relations from the sieving stage
	num_relations is the size of relation_list
	cycle_list is the list of cycles that the linear algebra
		code constructs (on input it is the list of cycles that
		the QS filtering code has found)
	num_cycles is the number of cycles. On input this is the size of
		cycle_list. The linear algebra code can change this number; 
		the only guarantee is that its final value is at least 
		fb_size + NUM_EXTRA_RELATIONS */

void solve_linear_system(msieve_obj *obj, uint32 fb_size, 
		    uint64 **bitfield, 
		    relation_t *relation_list, 
		    la_col_t *cycle_list,
		    uint32 *num_cycles);

/*-------------- MPQS SQUARE ROOT RELATED DECLARATIONS ---------------------*/

uint32 find_factors(msieve_obj *obj, mp_t *n, fb_t *factor_base, 
			uint32 fb_size, la_col_t *vectors, 
			uint32 vsize, relation_t *relation_list,
			uint64 *null_vectors, uint32 multiplier,
			mp_t *a_list, poly_t *poly_list, 
			factor_list_t *factor_list);
#ifdef __cplusplus
}
#endif

#endif /* _MPQS_MPQS_H_ */
