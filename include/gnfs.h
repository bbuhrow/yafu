/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#ifndef _GNFS_GNFS_H_
#define _GNFS_GNFS_H_

/* An implementation of the General Number Field
   Sieve algorithm for integer factorization. */

/* include basic stuff */

#include <common.h>
#include "msieve.h"
#include "yafu.h"
#include "factor.h"

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------- general stuff ---------------------------*/

#define MAX_POLY_DEGREE 8

/* representation of polynomials with multiple-
   precision coefficients. For polynomial p(x),
   element i of coeff[] gives the coefficient
   of x^i */

typedef struct {
	uint32 degree;
	mpz_t coeff[MAX_POLY_DEGREE + 1];

	/* scratch quantities for evaluating the homogeneous
	   form of poly */
	mpz_t tmp1, tmp2, tmp3;
} mpz_poly_t;

static INLINE void mpz_poly_init(mpz_poly_t * poly) {
	uint32 i;

	memset(poly, 0, sizeof(mpz_poly_t));

	mpz_init(poly->tmp1);
	mpz_init(poly->tmp2);
	mpz_init(poly->tmp3);
	for (i = 0; i <= MAX_POLY_DEGREE; i++)
		mpz_init_set_ui(poly->coeff[i], 0);
}

static INLINE void mpz_poly_free(mpz_poly_t * poly) {
	uint32 i;

	mpz_clear(poly->tmp1);
	mpz_clear(poly->tmp2);
	mpz_clear(poly->tmp3);
	for (i = 0; i <= MAX_POLY_DEGREE; i++)
		mpz_clear(poly->coeff[i]);
}

/* evaluate the homogeneous form of poly(x). If poly has
   degree d, then res = (b ^ d) * poly(a / b) */

void eval_poly(mpz_t res, int64 a, uint32 b, mpz_poly_t *poly);

typedef struct {
	int64 a;
	uint32 b;
} abpair_t;

/* Configuration for NFS parameters */

typedef struct {
	uint32 bits;       /* size of integer this config info applies to */
	uint32 rfb_limit;   /* largest rational factor base prime */
	uint32 afb_limit;   /* largest algebraic factor base prime */
	uint32 rfb_lp_size;   /* size of rational large primes */
	uint32 afb_lp_size;   /* size of algebraic large primes */
	uint64 sieve_size;  /* default length of sieve interval (actual 
			       interval is 2x this size) */
	int64 sieve_begin;  /* bounds of sieving interval; these default to */
	int64 sieve_end;    /* plus-or-minus sieve_size */
	double skewness;
} sieve_param_t;

/*---------------------- finite field poly stuff ---------------------*/

/* reduce the coefficients of _f modulo p, then compute
   all the x for which _f(x) = 0 mod p. The number of
   roots found, and the leading coefficient of _f mod p,
   is returned. If count_only is zero, the roots are 
   also returned in zeros[] */

uint32 poly_get_zeros(uint32 *zeros, 
			mpz_poly_t *_f, 
			uint32 p,
			uint32 *high_coeff,
			uint32 count_only);

/* like poly_get_zeros, except mult[i] is nonzero if
   zeros[i] is a multiple zero of _f mod p */

uint32 poly_get_zeros_and_mult(uint32 *zeros, 
			uint32 *mult,
			mpz_poly_t *_f, 
			uint32 p,
			uint32 *high_coeff);

/* return 1 if poly cannot be expressed as the product 
   of some other polynomials with coefficients modulo p,
   zero otherwise */

uint32 is_irreducible(mpz_poly_t *poly, uint32 p);

/* compute the inverse square root of the polynomial s_in,
   modulo the monic polynomial f_in, with all coefficients
   reduced modulo q. Returns 1 if the root is found and 
   zero otherwise */

uint32 inv_sqrt_mod_q(mpz_poly_t *res, mpz_poly_t *s_in, 
			mpz_poly_t *f_in, uint32 q, 
			uint32 *rand_seed1, uint32 *rand_seed2);

/*---------------------- factor base stuff ---------------------------*/

/* general entry in the factor base */

typedef struct {
	uint32 p;   /* prime for the entry */
	uint32 r;   /* the root of polynomial mod p for this
			entry, or p for projective roots */
} fb_entry_t;

/* rational and algebraic factor bases are treated
   the same as often as possible */

typedef struct {
	mpz_poly_t poly;        /* rational or algebraic polynomial */
	uint32 max_prime;       /* largest prime in the factor base */
	uint32 num_entries;     /* number of factor base entries */
	uint32 num_alloc;       /* amount allocated for FB entries */
	fb_entry_t *entries;    /* list of factor base entries */
} fb_side_t;

/* the NFS factor base */

typedef struct {
	fb_side_t rfb;    /* rational factor base */
	fb_side_t afb;    /* algebraic factor base */
} factor_base_t;

/* Given a factor base fb with the polynomials and the 
   maximum size primes filled in, fill in the rest of
   the entries in fb */

void create_factor_base(msieve_obj *obj, 
			factor_base_t *fb, 
			uint32 report_progress);

/* read / write / free a factor base */

int32 read_factor_base(msieve_obj *obj, mpz_t n,
		     sieve_param_t *params, factor_base_t *fb);

void write_factor_base(msieve_obj *obj, mpz_t n,
			sieve_param_t *params, factor_base_t *fb);

void free_factor_base(factor_base_t *fb);

/*---------------------- poly selection stuff ---------------------------*/

/* select NFS polynomials for factoring n, save
   in rat_poly and alg_poly */

int32 find_poly(msieve_obj *obj, mpz_t n);

/* attempt to read NFS polynomials from the factor 
   base file, save them and return 0 if successful.
   Skewness is ignored if NULL */

int32 read_poly(msieve_obj *obj, mpz_t n,
	       mpz_poly_t *rat_poly,
	       mpz_poly_t *alg_poly,
	       double *skewness);

/* unconditionally write the input NFS polynomials
   to a new factor base file. Skewness is ignored
   if < 0 */

void write_poly(msieve_obj *obj, mpz_t n,
	       mpz_poly_t *rat_poly,
	       mpz_poly_t *alg_poly,
	       double skewness);

/* determine the size and root properties of one polynomial */

void analyze_one_poly(msieve_obj *obj,
	       mpz_poly_t *rat_poly,
	       mpz_poly_t *alg_poly,
	       double skewness);

/*---------------------- sieving stuff ----------------------------------*/

/* external interface to perform sieving. The number of
   relations in the savefile at the time sieving completed
   is returned */

uint32 do_line_sieving(msieve_obj *obj, 
			sieve_param_t *params,
			mpz_t n, uint32 start_relations,
			uint32 max_relations);

/* the largest prime to be used in free relations */

#define FREE_RELATION_LIMIT (1 << 28)

/* add free relations to the savefile for this factorization;
   the candidates to add are bit entries in free_bits (which is
   compressed by a factor of 2). Returns the number of relations 
   that were added */

uint32 add_free_relations(msieve_obj *obj, factor_base_t *fb,
			  uint8 *free_bits);

/*---------------------- filtering stuff --------------------------------*/

/* external interface to filter relations. The return value is zero
   if filtering succeeded and the linear algebra can run, otherwise
   the estimated number of relations still needed before filtering
   could succeed is returned */

uint32 nfs_filter_relations(msieve_obj *obj, mpz_t n);

/*---------------------- linear algebra stuff ----------------------------*/

/* the minimum number of excess columns in the final
   matrix generated from relations. Note that the value 
   chosen contains a healthy fudge factor */

#define NUM_EXTRA_RELATIONS 200

#define TEMP_FACTOR_LIST_SIZE 100 
// bill: this is in Msieve's *current* common.h (but not the current gnfs.h), 
// which is completely different from yafu's common.h
// ditto on this struct, which is included by the next struct:
typedef struct {
	uint32 num_relations;  /* number of relations in the cycle */
	uint32 *list;          /* list of offsets into an array of relations */
} la_cycle_t;
// ditto on this struct:
typedef struct {
	uint32 *data;           /* The list of occupied rows in this column */
	uint32 weight;          /* Number of nonzero entries in this column */
	la_cycle_t cycle;       /* list of relations comprising this column */
} la_col_t;

/* external interface for NFS linear algebra */

void nfs_solve_linear_system(msieve_obj *obj, mpz_t n);

/* The largest prime ideal that is stored in compressed format
   when the matrix is built. Setting this to zero will cause
   all matrix rows to be stored in uncompressed format */

#define MAX_PACKED_PRIME 97

/*------------------------ square root stuff --------------------------*/

uint32 nfs_find_factors(msieve_obj *obj, mpz_t n, 
			factor_list_t *factor_list);

/*------------------- relation processing stuff --------------------------*/

#define RATIONAL_IDEAL 0
#define ALGEBRAIC_IDEAL 1

/* canonical representation of an ideal. NFS filtering
   and linear algebra will only work if the different
   ideals that occur in a factorization map to unique
   values of this structure. 
   
   Every ideal has a prime p and root r. To save space
   but still allow 48-bit p we store (p-1)/2 (thanks to
   Alex Kruppa for this trick) */

typedef struct {
	uint32 p_lo;  		/* (p - 1) / 2 (low 32 bits) */
	uint32 r_lo;            /* root for ideal (low 32 bits) */
	uint16 p_hi : 15;       /* (p - 1) / 2 (high 16 bits) */
	uint16 rat_or_alg : 1;  /* RATIONAL_IDEAL, ALGEBRAIC_IDEAL */
	uint16 r_hi;            /* root for ideal (high 16 bits) */
} ideal_t;

/* factors of relations are stored in a runlength-
   encoded format; 7 bits of each byte are used
   to store the factor, with the top bit set to
   1 if the corresponding byte is the most significant
   for the factor */

static INLINE uint32 compress_p(uint8 *array, 
				uint64 p, uint32 offset) {
	do {
		array[offset++] = p & 0x7f;
		p >>= 7;
	} while (p != 0);

	array[offset - 1] |= 0x80;
	return offset;
}

static INLINE uint64 decompress_p(uint8 *array, uint32 *offset_in) {

	uint32 offset = *offset_in;
	uint8 next_byte = array[offset++];
	uint64 p = next_byte & 0x7f;
	uint32 shift = 7;

	while (!(next_byte & 0x80)) {
		next_byte = array[offset++];
		p |= ((uint64)(next_byte & 0x7f)) << shift;
		shift += 7;
	} 
	
	*offset_in = offset;
	return p;
}

/* canonical representation of a relation, used in
   the NFS postprocessing phase */

#define COMPRESSED_P_MAX_SIZE 256

typedef struct relation_t {
	int64 a;               /* coordinates of relation; free relations */
	uint32 b;              /*   have b = 0 */
	uint32 rel_index;      /* line of savefile where relation occurs */
	uint8 num_factors_r;   /* number of rational factors */
	uint8 num_factors_a;   /* number of algebraic factors */
	uint8 *factors;       /* compressed list of factors; rational 
				 factors first, then algebraic */
} relation_t;

/* structure used to conveniently represent all of
   the large ideals that occur in a relation. The
   structure is far too large to be useful when 
   representing large groups of relations, so in
   these applications the data should be transferred
   to other containers once it is filled in here.
   Note that all of the rational ideals are listed
   first, then all of the algebraic ideals */

typedef struct {
	uint32 rel_index;         /* line of savefile where relation occurs */
	uint8 ideal_count;        /* count of large ideals */
	uint8 gf2_factors;        /* count of ideals not listed in ideal_list */
	ideal_t ideal_list[TEMP_FACTOR_LIST_SIZE];
} relation_lp_t;

/* convert a line of text into a relation_t, return 0
   if conversion succeeds. The array for factos pointed 
   to within 'r' should have at least COMPRESSED_P_MAX_SIZE 
   bytes. If 'compress' is nonzero then store only one 
   instance of any factors, and only if the factor occurs 
   in r an odd number of times */

int32 nfs_read_relation(char *buf, factor_base_t *fb, 
			relation_t *r, uint32 *array_size_out,
			uint32 compress, mpz_t scratch);

/* given a relation, find and list all of the rational
   ideals > filtmin_r and all of the algebraic ideals 
   whose prime exceeds filtmin_a. If these bounds are 
   zero then all ideals are listed */

uint32 find_large_ideals(relation_t *rel, relation_lp_t *out, 
			uint32 filtmin_r, uint32 filtmin_a);

/* Assuming a group of relations has been grouped together
   into a collection of cycles, read the collection of cycles
   from disk, read the relations they need, and link the two
   collections together. If 'compress' is nonzero then relations
   get only one instance of any factors, and only if the factor 
   occurs an odd number of times in the relation. If dependency
   is nonzero, only the cycles and relations required by that
   one dependency are read in. If fb is NULL, only the cycles
   (and not the relations they need) are read in */
   
void nfs_read_cycles(msieve_obj *obj, factor_base_t *fb, uint32 *ncols, 
			la_col_t **cols, uint32 *num_relations,
			relation_t **relation_list, uint32 compress,
			uint32 dependency);

void nfs_free_relation_list(relation_t *rlist, uint32 num_relations);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_GNFS_H_ */
