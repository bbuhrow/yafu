#ifndef _COFACTORIZE_H_
#define _COFACTORIZE_H_

#include <gmp.h>
#include "ytools.h"

typedef signed char s8;
typedef unsigned char u8;
typedef signed short s16;
typedef unsigned short u16;
typedef signed int s32;
typedef unsigned int u32;

#ifdef _MSC_VER
typedef signed __int64 s64;
typedef unsigned __int64 u64;
#ifndef INLINE
#define INLINE __inline
#endif
#else
typedef long long s64;
typedef unsigned long long u64;
#ifndef INLINE
#define INLINE __inline
#endif
#endif

/* maximum size pool of primes from which
   factor base is constructed */
#define NUM_PRIMES_TINY 1024

/* the number of dependencies the linear algebra
	will find */
#define NUM_EXTRA_RELATIONS_TINY 16

/* largest number of relations that can go into the
	linear algebra (includes relations combined from
	pairs of partial relations */
#define MAX_RELATIONS_TINY 512

		 /* the largest possible factor base */
#define MAX_FB_SIZE_TINY (MAX_RELATIONS_TINY - \
                          NUM_EXTRA_RELATIONS_TINY)

/* offset of the first valid factor base prime */
#define MIN_FB_OFFSET_TINY 1

/* offset of the first factor base prime
   actually contributing to the sieving */
#define MIN_FB_OFFSET_TO_SIEVE_TINY 7

   /* number of primes used when testing multipliers */
#define NUM_TEST_PRIMES_TINY 30

/* fudge factor to the target sieve value to account
   for not sieving with the smallest factor base primes.
   BRB: it is decidedly better to use smaller values
   for small problem sizes.  Above 116 bits needs yet
   another tier to increase accuracy. */
#define SMALL_PRIME_FUDGE_TINY 1
#define SMALL_PRIME_FUDGE 6

   /* maximum number of MPQS polynomials to be computed */
#define MAX_POLY_TINY 256

/* maximum number of FB primes that contribute to
   a single polynomial 'A' value */
#define MAX_POLY_FACTORS_TINY 5

/* the size of the sieve interval. Each polynomial will
	sieve over this many positive and negative values
	BRB: it is decidedly better to use smaller values
	for small problem sizes. Above 116 bits needs yet
	another tier to increase accuracy. */
#define SIEVE_SIZE_TINY 4096
#define SIEVE_SIZE 16384

/* value of the sieve root used when sieving is not
	to be performed for a given FB prime. Since this is
	larger than SIEVE_SIZE_TINY no special-case code
	is needed in the core sieve code */
#define DO_NOT_SIEVE_TINY 65535

/* maximum number of factors a relation can have (the
large prime is stored separately) */
#define MAX_FACTORS_TINY 40

/* partial relations are listed in the order
	in which they occur, and a hashtable matches
	up partial relations with the same large prime. */
#define LOG2_PARTIAL_TABLE_SIZE 10
#define LARGE_PRIME_HASH(x) (((u32)(x) * ((u32)40499 * 65543)) >> \
                                (32 - LOG2_PARTIAL_TABLE_SIZE))

/* number of collisions allowed in one hashtable entry */
#define LP_HASH_DEPTH_TINY 3

/* scale factor for all log values */
#define LOGPRIME_SCALE_TINY 2

/* maximum number of relations to be saved for
   resieving, used in place of trial factoring */
#define SIEVE_BATCH_SIZE_TINY 128

/* maximum size of the pool of FB primes that
	can appear in a polynomial 'A' value */
#define POLY_SELECT_BITS_TINY 12

#define POSITIVE 0
#define NEGATIVE 1




/* structure describing a single relation */

typedef struct {
	u32 large_prime;      /* the large prime (may be 1) */
	u16 fb_offsets[MAX_FACTORS_TINY]; /* offsets into FB of primes that
										 divide this relation */
	s16 sieve_offset;     /* the sieve offset of the relation */
	u8 poly_num;          /* ID of the poly that produce the relation */
	u8 num_factors;       /* number of factors from the factor base
							 (duplicates count) */
} tiny_relation;

/* structure describing one SIQS polynomial */

typedef struct {
	u16 a_fb_offsets[MAX_POLY_FACTORS_TINY];  /* factors of 'A' value */
	mpz_t b;                                  /* B value */
} tiny_poly;

/* main structure controlling the factorization */

typedef struct {

	/* basic stuff */

	mpz_t n;                          /* number to be factored */
	u32 multiplier;                   /* small multiplier of n */
	u16 multiplier_fb[2];             /* fb offsets of factors of multiplier */

	/* polynomial selection stuff */

	double target_a;                  /* the optimal size of poly A values */
	s32 poly_num;                     /* ID of current polynomial */
	s32 num_a_factors;                /* # of factors in poly 'A' values */
	s32 poly_select_idx;              /* ID of the combination of primes
										 that will make current A value */
	u16 poly_select_offsets[POLY_SELECT_BITS_TINY]; /* pool of primes for A */
	mpz_t poly_b_aux[MAX_POLY_FACTORS_TINY];      /* scratch values for com-
													 puting poly B values */
	tiny_poly poly_list[MAX_POLY_TINY];      /* list of SIQS polynomials */

	/* sieve stuff */

	double align_me;
	u8 sieve_block[2 * SIEVE_SIZE];  /* the sieve interval (8-byte aligned) */

	/* factor base stuff */

	s32 fb_size;                      /* number of FB primes */
	u16 prime_list[NUM_PRIMES_TINY];  /* complete list of primes from which
										 factor base is generated */
	float test_prime_contrib[NUM_TEST_PRIMES_TINY]; /* scratch space used in
													   multiplier selection */

													   /* relation stuff */

	s32 num_full_relations;   /* where next full relation will go */
	s32 partial_idx;          /* where next partial relation will go */
	s32 large_prime_max;      /* max value of a large prime */
	s32 error_bits;           /* value used for trial factoring cutoff */
	tiny_relation sieve_batch[SIEVE_BATCH_SIZE_TINY]; /* resieved relations */
	s32 sieve_hit_locs[SIEVE_BATCH_SIZE_TINY]; /* indices in the sieve block that need trial division */

	/* all relations that survive sieving are put in relation_list.
	   Full relations (and partial relations whose large prime has
	   occurred more than once) are stored in a list that grows up
	   from the beginning of the list, while partial relations that
	   have not been matched up yet are stored in a list growing down
	   from the end of relation_list. num_full_relations is the index
	   of the first free space for full relations, and partial_idx
	   does the same for unmatched partial relations. */

	tiny_relation relation_list[4 * MAX_RELATIONS_TINY];

	/* a hashtable is used to match up partial relations, using the
	   large prime as a hash key. The hashtable stores the index in
	   relation_list of the partial relation that connects up all the
	   other partial relations with the same large prime (those other
	   relations are treated as full relations) */

	u16 partial_hash[1 << LOG2_PARTIAL_TABLE_SIZE][LP_HASH_DEPTH_TINY];

	/* linear algebra stuff */

	u16 null_vectors[MAX_RELATIONS_TINY];
	u64 matrix[MAX_FB_SIZE_TINY][(MAX_RELATIONS_TINY + 63) / 64];

	/* data used during a factorization */

	ALIGNED_MEM  u8 glogprimes[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 gprimes[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 gmodsqrt[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 roots1[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 roots2[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 root_aux[MAX_POLY_FACTORS_TINY * MAX_FB_SIZE_TINY];      /* scratch value for initializing sieve roots */
	ALIGNED_MEM  u32 grecip[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 grecip16[MAX_FB_SIZE_TINY];

	mpz_t gmptmp1;
	mpz_t gmptmp2;
	mpz_t gmptmp3;
	mpz_t gmptmp4;
	mpz_t gmptmp5;

} tiny_qs_params;


// top level routines
tiny_qs_params * init_tinyqs(void);
u32 tinyqs(tiny_qs_params *g_params, mpz_t n, mpz_t factor1, mpz_t factor2);
tiny_qs_params *free_tinyqs(tiny_qs_params *g_params);

#endif /* !_COFACTORIZE_H_ */


#ifndef _COFACTORIZE_H_
#define _COFACTORIZE_H_

#include <gmp.h>

typedef signed char s8;
typedef unsigned char u8;
typedef signed short s16;
typedef unsigned short u16;
typedef signed int s32;
typedef unsigned int u32;

#ifdef _MSC_VER
typedef signed __int64 s64;
typedef unsigned __int64 u64;
#define INLINE _inline
#else
typedef long long s64;
typedef unsigned long long u64;
#define INLINE inline
#endif

static INLINE u64 gmp2u64(mpz_t src) 
{
  /* mpz_export is terribly slow */
  u64 ans = mpz_getlimbn(src, 0);
#if GMP_LIMB_BITS == 32
  if (mpz_size(src) >= 2)
    ans |= (u64)mpz_getlimbn(src, 1) << 32;
#endif
  return ans;
}

static INLINE void u64_2gmp(u64 src, mpz_t dest) 
{
#if GMP_LIMB_BITS == 64
  dest->_mp_d[0] = src;
  dest->_mp_size = (src ? 1 : 0);
#else
  /* mpz_import is terribly slow */
  mpz_set_ui(dest, (u32)(src >> 32));
  mpz_mul_2exp(dest, dest, 32);
  mpz_add_ui(dest, dest, (u32)src);
#endif
}

/* maximum size pool of primes from which
   factor base is constructed */
#define NUM_PRIMES_TINY 1024

   /* the number of dependencies the linear algebra
	  will find */
#define NUM_EXTRA_RELATIONS_TINY 16

	  /* largest number of relations that can go into the
		 linear algebra (includes relations combined from
		 pairs of partial relations */
#define MAX_RELATIONS_TINY 512

		 /* the largest possible factor base */
#define MAX_FB_SIZE_TINY (MAX_RELATIONS_TINY - \
                          NUM_EXTRA_RELATIONS_TINY)

/* offset of the first valid factor base prime */
#define MIN_FB_OFFSET_TINY 1

/* offset of the first factor base prime
   actually contributing to the sieving */
#define MIN_FB_OFFSET_TO_SIEVE_TINY 7

   /* number of primes used when testing multipliers */
#define NUM_TEST_PRIMES_TINY 30

/* fudge factor to the target sieve value to account
   for not sieving with the smallest factor base primes.
   BRB: it is decidedly better to use smaller values
   for small problem sizes.  Above 116 bits needs yet
   another tier to increase accuracy. */
#define SMALL_PRIME_FUDGE_TINY 1
#define SMALL_PRIME_FUDGE 6

   /* maximum number of MPQS polynomials to be computed */
#define MAX_POLY_TINY 256

/* maximum number of FB primes that contribute to
   a single polynomial 'A' value */
#define MAX_POLY_FACTORS_TINY 5

   /* the size of the sieve interval. Each polynomial will
	  sieve over this many positive and negative values
	  BRB: it is decidedly better to use smaller values
	  for small problem sizes. Above 116 bits needs yet
	  another tier to increase accuracy. */
#define SIEVE_SIZE_TINY 4096
#define SIEVE_SIZE 16384

	  /* value of the sieve root used when sieving is not
		 to be performed for a given FB prime. Since this is
		 larger than SIEVE_SIZE_TINY no special-case code
		 is needed in the core sieve code */
#define DO_NOT_SIEVE_TINY 65535

		 /* maximum number of factors a relation can have (the
			large prime is stored separately) */
#define MAX_FACTORS_TINY 40

			/* partial relations are listed in the order
			   in which they occur, and a hashtable matches
			   up partial relations with the same large prime. */
#define LOG2_PARTIAL_TABLE_SIZE 10
#define LARGE_PRIME_HASH(x) (((u32)(x) * ((u32)40499 * 65543)) >> \
                                (32 - LOG2_PARTIAL_TABLE_SIZE))

			   /* number of collisions allowed in one hashtable entry */
#define LP_HASH_DEPTH_TINY 3

/* scale factor for all log values */
#define LOGPRIME_SCALE_TINY 2

/* maximum number of relations to be saved for
   resieving, used in place of trial factoring */
#define SIEVE_BATCH_SIZE_TINY 128

   /* maximum size of the pool of FB primes that
	  can appear in a polynomial 'A' value */
#define POLY_SELECT_BITS_TINY 12

#define POSITIVE 0
#define NEGATIVE 1




/* structure describing a single relation */

typedef struct {
	u32 large_prime;      /* the large prime (may be 1) */
	u16 fb_offsets[MAX_FACTORS_TINY]; /* offsets into FB of primes that
										 divide this relation */
	s16 sieve_offset;     /* the sieve offset of the relation */
	u8 poly_num;          /* ID of the poly that produce the relation */
	u8 num_factors;       /* number of factors from the factor base
							 (duplicates count) */
} tiny_relation;

/* structure describing one SIQS polynomial */

typedef struct {
	u16 a_fb_offsets[MAX_POLY_FACTORS_TINY];  /* factors of 'A' value */
	mpz_t b;                                  /* B value */
} tiny_poly;

/* main structure controlling the factorization */

typedef struct {

	/* basic stuff */

	mpz_t n;                          /* number to be factored */
	u32 multiplier;                   /* small multiplier of n */
	u16 multiplier_fb[2];             /* fb offsets of factors of multiplier */

	/* polynomial selection stuff */

	double target_a;                  /* the optimal size of poly A values */
	s32 poly_num;                     /* ID of current polynomial */
	s32 num_a_factors;                /* # of factors in poly 'A' values */
	s32 poly_select_idx;              /* ID of the combination of primes
										 that will make current A value */
	u16 poly_select_offsets[POLY_SELECT_BITS_TINY]; /* pool of primes for A */
	mpz_t poly_b_aux[MAX_POLY_FACTORS_TINY];      /* scratch values for com-
													 puting poly B values */
	tiny_poly poly_list[MAX_POLY_TINY];      /* list of SIQS polynomials */

	/* sieve stuff */

	double align_me;
	u8 sieve_block[2 * SIEVE_SIZE];  /* the sieve interval (8-byte aligned) */

	/* factor base stuff */

	s32 fb_size;                      /* number of FB primes */
	u16 prime_list[NUM_PRIMES_TINY];  /* complete list of primes from which
										 factor base is generated */
	float test_prime_contrib[NUM_TEST_PRIMES_TINY]; /* scratch space used in
													   multiplier selection */

													   /* relation stuff */

	s32 num_full_relations;   /* where next full relation will go */
	s32 partial_idx;          /* where next partial relation will go */
	s32 large_prime_max;      /* max value of a large prime */
	s32 error_bits;           /* value used for trial factoring cutoff */
	tiny_relation sieve_batch[SIEVE_BATCH_SIZE_TINY]; /* resieved relations */
	s32 sieve_hit_locs[SIEVE_BATCH_SIZE_TINY]; /* indices in the sieve block that need trial division */

	/* all relations that survive sieving are put in relation_list.
	   Full relations (and partial relations whose large prime has
	   occurred more than once) are stored in a list that grows up
	   from the beginning of the list, while partial relations that
	   have not been matched up yet are stored in a list growing down
	   from the end of relation_list. num_full_relations is the index
	   of the first free space for full relations, and partial_idx
	   does the same for unmatched partial relations. */

	tiny_relation relation_list[4 * MAX_RELATIONS_TINY];

	/* a hashtable is used to match up partial relations, using the
	   large prime as a hash key. The hashtable stores the index in
	   relation_list of the partial relation that connects up all the
	   other partial relations with the same large prime (those other
	   relations are treated as full relations) */

	u16 partial_hash[1 << LOG2_PARTIAL_TABLE_SIZE][LP_HASH_DEPTH_TINY];

	/* linear algebra stuff */

	u16 null_vectors[MAX_RELATIONS_TINY];
	u64 matrix[MAX_FB_SIZE_TINY][(MAX_RELATIONS_TINY + 63) / 64];

	/* data used during a factorization */

	ALIGNED_MEM  u8 glogprimes[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 gprimes[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 gmodsqrt[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 roots1[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 roots2[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 root_aux[MAX_POLY_FACTORS_TINY * MAX_FB_SIZE_TINY];      /* scratch value for initializing sieve roots */
	ALIGNED_MEM  u32 grecip[MAX_FB_SIZE_TINY];
	ALIGNED_MEM  u16 grecip16[MAX_FB_SIZE_TINY];

	mpz_t gmptmp1;
	mpz_t gmptmp2;
	mpz_t gmptmp3;
	mpz_t gmptmp4;
	mpz_t gmptmp5;

} tiny_qs_params;


// top level routines
tiny_qs_params * init_tinyqs(void);
u32 tinyqs(tiny_qs_params *g_params, mpz_t n, mpz_t factor1, mpz_t factor2);

#endif /* !_COFACTORIZE_H_ */

