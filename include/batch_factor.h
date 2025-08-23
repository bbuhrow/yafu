/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: batch_factor.h 638 2011-09-11 15:31:19Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _BATCH_FACTOR_H_
#define _BATCH_FACTOR_H_

#include <stdint.h>
#include "ytools.h"
#include "common.h"
#include "cofactorize.h"
#include "gmp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* prototypes for the subsystem that batch-factors the
   portions of relations that contain large primes. The
   current system is designed for cofactors containing up to
   three large primes, though relations with any number of
   large primes will be batched and factored */

#define MAX_LARGE_PRIMES 4

//typedef void (*print_relation_t)(savefile_t *savefile, int64_t a, uint32_t b,
//			uint32_t *factors_r, uint32_t num_factors_r, 
//			uint32_t lp_r[MAX_LARGE_PRIMES],
//			uint32_t *factors_a, uint32_t num_factors_a, 
//			uint32_t lp_a[MAX_LARGE_PRIMES]);
//
//void print_relation_dummy(savefile_t *savefile, int64_t a, uint32_t b,
//    uint32_t *factors_r, uint32_t num_factors_r,
//    uint32_t lp_r[MAX_LARGE_PRIMES],
//    uint32_t *factors_a, uint32_t num_factors_a,
//    uint32_t lp_a[MAX_LARGE_PRIMES]);

typedef void (*print_relation_t)(void);

/* simplified representation of one relation. Note that
   any of the factors may be trivial */

typedef struct {
    uint64_t a;               // for siqs, index of a-poly in master list (was an int64)
	uint32_t b;               // for siqs, index of b-poly within the indicated a-poly
	uint32_t extra_f;
    int32_t signed_offset;    // for siqs, location of relation in the block (new to struct)
	uint8_t num_factors_r;     /* doesn't include large primes */
	uint8_t num_factors_a;     /* doesn't include large primes */
	uint8_t lp_r_num_words;     /* one number with this many words */
	uint8_t lp_a_num_words;     /* one number with this many words */
	uint32_t factor_list_word;  /* offset in an array of factors where
				            all of the above above appear, in order */
    uint8_t success;
    uint32_t lp_r[MAX_LARGE_PRIMES];
} cofactor_t;

/* main structure controlling batch factoring. The main goal
   of batch factoring is to compute gcd(each_relation, 
   product_of_many_primes) much faster than running a conventional
   factoring algorithm on each relation. We hope that the gcd
   splits a given relation into two approximately equal-size halves,
   so that each part is trivial (has one factor) or is easy (has two
   factors) to decompose; failing that, we hope that the gcd provides
   enough information to tell us whether going to the trouble of
   fully factoring the relation is unlikely to yield something useful.

   The algorithm depends critically on the fact that gcd(A,B) =
   gcd(A mod B, B), so that when B is massively smaller than A the
   former is easy to compute. We don't actually compute the gcd
   until each relation B has (A mod B) available, and it is the latter
   that Bernstein's algorithm computes */

typedef struct {
	mpz_t prime_product;  /* product of primes used in the gcd */

	tiny_qs_params *params;
	uint32_t num_uecm[4];		/* calls to uecm to:
									(0) split 2LP on the f1r side
									(1) split 2LP on the f2r side
									(2) split 2LP on the small side of TLP residues
									(3) split 2LP on the large side of TLP residues */
	uint32_t num_tecm;			/* calls to tecm to split 3LP on the f1r side */
	uint32_t num_tecm2;			/* calls to tecm to split 4LP on the f1r side */
	uint32_t num_qs;			/* tecm failures (in the future, handled by qs) */
	uint32_t num_attempt;		/* number of calls to check_relation (possibly involving ecm) */
	uint32_t num_success;       /* number of surviving relations */
	uint32_t num_abort[8];		/* stop because:
								0) f2r has one factor that's too big, 
								1) f2r has 2 factors, at least one of which is too big (residue > LPB^2),
								2) f2r contains a large prime: residue is less than (max-GCD-prime)^2 
								3) f2r has more than 2 factors; is > 64 bits
								4) prime double-word f2r
								5) non-useful f1r double-word split (factor larger than LPB)
								6) non-useful f2r double-word split (factor larger than LPB)
								7) non-useful f1r TLP split (factor larger than LPB), at any point in the process */
	uint32_t target_relations;  /* number of relations to batch up */
	uint64_t lp_cutoff_r;       /* maximum size of rational factors */
	mpz_t lp_cutoff_r2;        /* square of lp_cutoff_r */
	mpz_t lp_cutoff_r3;        /* cube of lp_cutoff_r */
	uint64_t lp_cutoff_a;       /* maximum size of algebraic factors */
    mpz_t lp_cutoff_a2         /* square of lp_cutoff_a */;
	mpz_t lp_cutoff_a3         /* cube of lp_cutoff_a */;
	mpz_t min_prime2;          /* the square of the smallest prime that occurs in prime_product */
    mpz_t max_prime2;          /* the square of the largest prime that occurs in prime_product */
	mpz_t max_prime3;          /* the cube of the largest prime that occurs in prime_product */

	uint32_t num_relations;     /* number relations currently batched */
	uint32_t num_relations_alloc;     /* relations allocated */
	cofactor_t *relations;    /* relations currently batched */

	uint32_t num_factors;       /* words for factors of batched relations */
	uint32_t num_factors_alloc; /* space for batched factors */
	uint32_t *factors;          /* factors of batched relations */

	//savefile_t *savefile;
	print_relation_t print_relation;

    //mpz_t small, large;
    mpz_t _large;
    mpz_t _small;
    mpz_t n, f1r, f2r, f1a, f2a, t0, t1;

    double conversion_ratio;
} relation_batch_t;

/* initialize the relation batch. Batch factoring uses all
   the primes from min_prime to max_prime. We assume min_prime
   is the smallest element possible in unfactored parts of
   relations, i.e. MIN(max factor base prime within rational
   or algebraic factor bases).
   
   In general we do not want max_prime as large as the 
   large prime cutoffs; making it smaller allows the batch 
   factoring to split most of the cofactors in relations that 
   contain large primes, or at least prove most relations to 
   be not worth the trouble to do so manually */

void relation_batch_init(FILE *logfile, relation_batch_t *rb,
    uint32_t min_prime, uint64_t max_prime,
    uint64_t lp_cutoff_r, uint64_t lp_cutoff_a,
    print_relation_t print_relation,
    int do_prime_product);

void relation_batch_free(relation_batch_t *rb);

/* add one relation to the batch. Note that the relation may
   contain any number of large primes, as long as the unfactored
   parts of the relation fit into an mp_t. Modifying the code to
   handle more than three large primes within unfactored_[ra]
   just involves modifying the base case of the recursion; maybe 
   that should be made into a callback */

void relation_batch_add(uint64_t a, uint32_t b, int32_t offset,
			uint32_t *factors_r, uint32_t num_factors_r, 
			mpz_t unfactored_r,
			uint32_t *factors_a, uint32_t num_factors_a, 
			mpz_t unfactored_a,
            mpz_t tmp_in,
			relation_batch_t *rb);

void check_batch_relation(relation_batch_t *rb,
    uint32_t index, uint64_t * lcg_state,
    mpz_t prime_product);
	
/* factor all the batched relations, saving all the ones whose
   largest {rational|algebraic} factors are all less than
   lp_cutoff_[ra] */

uint32_t relation_batch_run(relation_batch_t *rb, uint64_t *lcg_state);

#ifdef __cplusplus
}
#endif

#endif /* _BATCH_FACTOR_H_ */
