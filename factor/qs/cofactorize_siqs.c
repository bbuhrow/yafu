#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#ifndef _MSC_VER
#include <sys/time.h>
#endif
#include <time.h>
#include "common.h"
#include "types.h"
#include "cofactorize.h"
#include <immintrin.h>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif

// gcc -O3 -march=core-avx2 cofactorize_siqs.c -o tsiqs2 -I/sppdg/scratch/buhrow/projects/factoring/gmp/include -L /sppdg/scratch/buhrow/projects/factoring/gmp/lib/linux/x86_64 -lm -lgmp

//#define VERBOSE 1
//#define DO_RESIEVE
#define USE_8X_TD
// BRB: I didn't see much variation when all functions are static, or not.
// But, if they are static then it was much harder to read anything in
// gprof profiling output... so being able to quickly remove the static
// label on all functions is nice.
#define func_static static

#if (defined(GCC_ASM64X) || defined(__MINGW64__)) && !defined(FORCE_GENERIC) && !defined(TARGET_KNC)
#define USE_ASM
// todo: need to define alternate routines if this isn't defined...
#endif

/* High-throughput, low-overhead implementation of the self-
   initializing multiple polynomial quadratic sieve, optimized
   for small inputs (50-120 bits). Many of the ideas here are
   extensions of the remarkable MPQS code of F. Bahr, used 
   in the lattice sievers by Jens Franke */

/* seeds for random numbers */

static u32 rand_seed1 = 11111111;
static u32 rand_seed2 = 22222222;

#define RAND_MULT 2131995753

static u32 get_rand(u32 *seed1, u32 *seed2)
{
  /* A multiply-with-carry generator by George Marsaglia.
     The period is about 2^63. */

  u64 temp = (u64)(*seed1) * (u64)RAND_MULT + (u64)(*seed2);
  *seed1 = (u32)temp;
  *seed2 = (u32)(temp >> 32);
  return (u32)temp;
}

/* masks for picking out individual bits of 64-bit
   words, used for the linear algebra */

#define B(x) ((u64)(1) << (x))

static const u64 bitmask[] = {
  B( 0), B( 1), B( 2), B( 3), B( 4), B( 5), B( 6), B( 7),
  B( 8), B( 9), B(10), B(11), B(12), B(13), B(14), B(15),
  B(16), B(17), B(18), B(19), B(20), B(21), B(22), B(23),
  B(24), B(25), B(26), B(27), B(28), B(29), B(30), B(31),
  B(32), B(33), B(34), B(35), B(36), B(37), B(38), B(39),
  B(40), B(41), B(42), B(43), B(44), B(45), B(46), B(47),
  B(48), B(49), B(50), B(51), B(52), B(53), B(54), B(55),
  B(56), B(57), B(58), B(59), B(60), B(61), B(62), B(63),
};

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
  u8 sieve_block[2*SIEVE_SIZE];  /* the sieve interval (8-byte aligned) */

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
  u64 matrix[MAX_FB_SIZE_TINY][(MAX_RELATIONS_TINY+63) / 64];
} tiny_qs_params;

ALIGNED_MEM  u8 glogprimes[MAX_FB_SIZE_TINY];
ALIGNED_MEM  u16 gprimes[MAX_FB_SIZE_TINY];
ALIGNED_MEM  u16 gmodsqrt[MAX_FB_SIZE_TINY];
ALIGNED_MEM  u16 roots1[MAX_FB_SIZE_TINY];
ALIGNED_MEM  u16 roots2[MAX_FB_SIZE_TINY];
ALIGNED_MEM  u16 root_aux[MAX_POLY_FACTORS_TINY * MAX_FB_SIZE_TINY];      /* scratch value for initializing sieve roots */
ALIGNED_MEM  u32 grecip[MAX_FB_SIZE_TINY];
ALIGNED_MEM  u16 grecip16[MAX_FB_SIZE_TINY];

/* the following is reused across factorizations */

static tiny_qs_params *g_params = NULL;

/* The following utility routines are not really
   a performance bottleneck, but since they always
   deal with 16-bit data at most their input 
   datatypes should really by u16's. This will 
   make all the division and remainder operations 
   a lot faster */

/***********************************/
func_static s32 legendre_16(s32 a, s32 p)
/***********************************
Compute the Legendre symbol (a/p)
************************************/
{
  s32 tmp;
  s32 x = a;
  s32 y = p;
  s32 out = 1;

  while (x) {
    while ((x & 1) == 0) {
      x = x / 2;
      if ( (y & 7) == 3 || (y & 7) == 5 )
        out = -out;
    }

    tmp = x;
    x = y;
    y = tmp;

    if ( (x & 3) == 3 && (y & 3) == 3 )
      out = -out;

    x = x % y;
  }
  if (y == 1)
    return out;
  return 0;
}

/***********************************/
func_static s32 powm_16(s32 a, s32 b, s32 n) 
/***********************************
Compute a^b mod n
************************************/
{
  s32 res = 1;
  while (b) {
    if (b & 1)
      res = res * a % n;
    a = a * a % n;
    b = b >> 1;
  }
  return res;
}

/***********************************/
func_static s32 modinv_16(s32 a, s32 p)
/***********************************
High-speed modular inverse of 'a' mod 'p'
Thanks to the folks at www.mersenneforum.com
for coming up with this
************************************/
{
  s32 ps1, ps2, parity, dividend, divisor, rem, q, t;

  q = 1;
  rem = a;
  dividend = p;
  divisor = a;
  ps1 = 1;
  ps2 = 0;
  parity = 0;

  while (divisor > 1) {
    rem = dividend - divisor;
    t = rem - divisor;
    if (t >= 0) {
      q += ps1; rem = t; t -= divisor;
      if (t >= 0) {
        q += ps1; rem = t; t -= divisor;
        if (t >= 0) {
          q += ps1; rem = t; t -= divisor;
          if (t >= 0) {
            q += ps1; rem = t; t -= divisor;
            if (t >= 0) {
              q += ps1; rem = t; t -= divisor;
              if (t >= 0) {
                q += ps1; rem = t; t -= divisor;
                if (t >= 0) {
                  q += ps1; rem = t; t -= divisor;
                  if (t >= 0) {
                    q += ps1; rem = t;
                    if (rem >= divisor) {
                      q = dividend / divisor;
                      rem = dividend % divisor;
                      q *= ps1;
                    } } } } } } } }
    }
    q += ps2;
    parity = ~parity;
    dividend = divisor;
    divisor = rem;
    ps2 = ps1;
    ps1 = q;
  }
  
  if (parity == 0)
    return ps1;
  else
    return p - ps1;
}

/***********************************/
func_static s32 sqrtModP_16(s32 a, s32 p)
/***********************************
Compute the square root of 'a' mod 'p'
This is Algorithm 2.3.8 from Crandall & 
Pomerance, "Prime Numbers: A Computational
Perspective"
************************************/
{
  if ( (p & 7) == 3 || (p & 7) == 7 ) {
    return powm_16(a, (p+1)/4, p);
  }
  else if ( (p & 7) == 5 ) {
    s32 x, y;
    
    x = powm_16(a, (p+3)/8, p);
    if ((x * x % p) == a)
      return x;

    y = powm_16(2, (p-1)/4, p);
    return (s32)x * y % p;
  }
  else {
    s32 i, d0, d1, a1, s, t, m;

    d0 = get_rand(&rand_seed1, &rand_seed2) % p;
    while (legendre_16(d0, p) != -1)
      d0 = get_rand(&rand_seed1, &rand_seed2) % p;

    t = p - 1;
    s = 0;
    while (!(t & 1)) {
      s++;
      t = t / 2;
    }

    a1 = powm_16(a, t, p);
    d1 = powm_16(d0, t, p);

    for (i = 0, m = 0; i < s; i++) {
      s32 ad = powm_16(d1, m, p);
      ad = ad * a1 % p;
      ad = powm_16(ad, (u16)(1) << (s-1-i), p);
      if (ad == (p - 1))
        m += (1 << i);
    }

    a1 = powm_16(a, (t+1)/2, p);
    d1 = powm_16(d1, m/2, p);
    return a1 * d1 % p;
  }
}


/***********************************/
func_static void init_tinyqs(void)
/***********************************/
{
  s32 i, j, k, rem;
  tiny_qs_params *p;

  if (g_params)
    return;

  /* allocate the main structure */

  p = g_params = (tiny_qs_params *)malloc(sizeof(tiny_qs_params));
  mpz_init(p->n);

  /* fill in the pool of primes */

  p->prime_list[0] = 2;
  p->prime_list[1] = 3;
  for (i = 2, j = 5; i < NUM_PRIMES_TINY; j += 2) {
    for (k = 1, rem = 0; k < i; k++) {
      s32 prime = p->prime_list[k];
      rem = j % prime;
      if (prime * prime > j || rem == 0)
        break;
    }
    if (rem != 0)
      p->prime_list[i++] = j;
  }

  /* init the scratch values for polynomial 'B'
     value computations */

  for (i = 0; i < MAX_POLY_FACTORS_TINY; i++) {
    mpz_init(p->poly_b_aux[i]);
  }

  /* set up the list of sieve polynomials */

  for (i = 0; i < MAX_POLY_TINY; i++) {
    mpz_init(p->poly_list[i].b);
  }

  /* see the next routine for an explanation of what
     these quantities are */

  for (i = 1; i < NUM_TEST_PRIMES_TINY; i++) {
    p->test_prime_contrib[i] = 2 * log((double)p->prime_list[i]) / 
                               (p->prime_list[i] - 1) / M_LN2;
  }
}

/* Implementation of the modified Knuth-Schroeppel multiplier
   algorithm. This borrows ideas from at least four different
   sources, and seems to choose multipliers that are better on
   average than many of the other methods available.
   
   There are many misconceptions about what this algorithm is
   supposed to do. We want to multiply the input number n by a
   small odd squarefree constant k, chosen so that the factor base 
   for k * n contains as many small primes as possible. Since small primes
   occur more often than big ones, this makes sieve values smaller
   on average and so more likely to be smooth. We quantify this
   by measuring the average contribution of the first NUM_TEST_PRIMES_TINY
   primes to sieve values. There are two constraints: first, larger 
   multipliers mean a larger number to factor. Second, we can't spend 
   all day testing multipliers, so the set of multipliers to test should 
   be small. 

   The list of available multipliers depends on the value of n mod
   8, 3, and 5; each row of the table below gives the multipliers
   to try, pre-sorted by how well they approximately optimize sieving
   (the routine below computes a better approximation). Note that a
   multiplier of 1 (i.e. no multiplier) is always possible. Experiments
   show that 90% of the time the optimal multiplier is in one of the 
   first four columns of the table */

#define MAX_MULTIPLIERS 13                           /* for residue classes: */
static u8 mult_list[32][MAX_MULTIPLIERS] = {         /* mod 8  mod 3  mod 5 */
{ 1, 19, 61, 31, 21, 13,  7,  3, 73, 41,  5, 33, 37 }, /*  1      1      1 */
{ 1, 13,  7,  3, 73, 33, 37, 17, 57, 43,  5, 19, 15 }, /*  1      1      2 */
{ 1, 13,  7,  3, 73, 33, 37, 17, 57, 43,  5, 19, 15 }, /*  1      1      3 */
{ 1, 19, 61, 31, 21, 13,  7,  3, 73, 41,  5, 33, 37 }, /*  1      1      4 */
{ 1, 41,  5, 17, 11, 89, 29, 65, 21,  3, 59, 33, 35 }, /*  1      2      1 */
{ 1, 17,  5,  3, 33, 65, 57, 23, 41, 53, 47, 11, 89 }, /*  1      2      2 */
{ 1, 17,  5,  3, 33, 65, 57, 23, 41, 53, 47, 11, 89 }, /*  1      2      3 */
{ 1, 41,  5, 17, 11, 89, 29, 65, 21,  3, 59, 33, 35 }, /*  1      2      4 */
{ 1, 19,  3, 11, 31,  7, 51, 43, 15, 39, 61, 55, 21 }, /*  3      1      1 */
{ 1,  3,  7, 43, 19, 13, 37, 15, 55, 11, 73, 31, 35 }, /*  3      1      2 */
{ 1,  3,  7, 43, 19, 13, 37, 15, 55, 11, 73, 31, 35 }, /*  3      1      3 */
{ 1, 19,  3, 11, 31,  7, 51, 43, 15, 39, 61, 55, 21 }, /*  3      1      4 */
{ 1, 11,  3, 59, 35,  5, 51, 19, 29, 41, 15, 23, 39 }, /*  3      2      1 */
{ 1,  3, 11, 35,  5, 23, 17, 47,  7, 59, 43, 15, 53 }, /*  3      2      2 */
{ 1,  3, 11, 35,  5, 23, 17, 47,  7, 59, 43, 15, 53 }, /*  3      2      3 */
{ 1, 11,  3, 59, 35,  5, 51, 19, 29, 41, 15, 23, 39 }, /*  3      2      4 */
{ 1, 61, 21, 13,  5, 19, 37, 31, 29,  7,  3, 11, 15 }, /*  5      1      1 */
{ 1, 13, 37,  7,  3,  5, 73, 61, 21, 43, 33, 53, 17 }, /*  5      1      2 */
{ 1, 13, 37,  7,  3,  5, 73, 61, 21, 43, 33, 53, 17 }, /*  5      1      3 */
{ 1, 61, 21, 13,  5, 19, 37, 31, 29,  7,  3, 11, 15 }, /*  5      1      4 */
{ 1,  5, 29, 21, 11, 41, 53, 17, 89,  3, 59, 61, 65 }, /*  5      2      1 */
{ 1,  5, 53, 17,  3, 13, 29, 23, 21, 37, 47, 33, 11 }, /*  5      2      2 */
{ 1,  5, 53, 17,  3, 13, 29, 23, 21, 37, 47, 33, 11 }, /*  5      2      3 */
{ 1,  5, 29, 21, 11, 41, 53, 17, 89,  3, 59, 61, 65 }, /*  5      2      4 */
{ 1, 31,  7, 19, 15, 39, 55,  3, 11, 61, 21, 13, 51 }, /*  7      1      1 */
{ 1,  7,  3, 15, 13, 55, 31, 43, 23, 37, 19, 47, 73 }, /*  7      1      2 */
{ 1,  7,  3, 15, 13, 55, 31, 43, 23, 37, 19, 47, 73 }, /*  7      1      3 */
{ 1, 31,  7, 19, 15, 39, 55,  3, 11, 61, 21, 13, 51 }, /*  7      1      4 */
{ 1, 11,  5, 15, 23, 39,  3, 29, 47, 59, 31, 35,  7 }, /*  7      2      1 */
{ 1, 23,  3, 47,  7,  5, 15, 17, 11, 35, 53, 39, 33 }, /*  7      2      2 */
{ 1, 23,  3, 47,  7,  5, 15, 17, 11, 35, 53, 39, 33 }, /*  7      2      3 */
{ 1, 11,  5, 15, 23, 39,  3, 29, 47, 59, 31, 35,  7 }, /*  7      2      4 */
};

/***********************************/
func_static void find_multiplier_tiny(void)
/***********************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  u16 *prime_list = params->prime_list;
  u16 test_nmodp[NUM_TEST_PRIMES_TINY];
  s32 best_mult = 1;
  s32 nmod8 = mpz_get_ui(params->n) % 8;
  float best_score;
  u8 *mult_row;
  s32 num_tests;

  /* precompute information that will be needed 
     for all multipliers */

  for (i = 1; i < NUM_TEST_PRIMES_TINY; i++)
    test_nmodp[i] = mpz_tdiv_ui(params->n, prime_list[i]);

  /* find the row of the table that is approriate for this
     value of n */

  mult_row = mult_list[ test_nmodp[2] - 1 +
                        4*(test_nmodp[1] - 1) + 
		        8*(nmod8 / 2) ];

  /* test less than the whole row if n is small */

  num_tests = mpz_sizeinbase(params->n, 2) / 10;
  if (num_tests > MAX_MULTIPLIERS)
    num_tests = MAX_MULTIPLIERS;

  best_score = 1000.0;
  for (i = 0; i < num_tests; i++) {
    s32 curr_mult = mult_row[i];
    s32 knmod8 = (nmod8 * curr_mult) % 8;
    float score;

    /* measure the contribution of 2 as a factor of sieve
       values. The multiplier itself must also be taken into
       account in the score. 'score' is the correction that
       is implicitly applied to the size of sieve values; a
       negative score makes sieve values smaller, and so is 
       better. */

    if (knmod8 == 1)
      score = 0.5 * log((double)curr_mult) / M_LN2 - 2;
    else if (knmod8 == 5)
      score = 0.5 * log((double)curr_mult) / M_LN2 - 1;
    else
      score = 0.5 * log((double)curr_mult) / M_LN2 - 0.5;

    for (j = 1; j < NUM_TEST_PRIMES_TINY; j++) {
      s32 prime = prime_list[j];
      s32 knmodp = (s32)test_nmodp[j] * curr_mult % prime;

      /* if prime j is actually in the factor base 
         for k * n ... */

      if (legendre_16(knmodp, prime) != -1) {

        /* ...add its contribution. A prime p con-
           tributes log(p) to 1 in p sieve values, plus
           log(p) to 1 in p^2 sieve values, etc. The
           average contribution of all multiples of p 
           to a random sieve value is thus

           log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
           = (log(p) / p) * 1 / (1 - (1/p)) 
           = log(p) / (p-1)

           This contribution occurs once for each
           square root used for sieving. There are two
           roots for each factor base prime, unless
           the prime divides the multiplier. In that
           case there is only one root. The scores are
           premultiplied by 2.0, and logarithms are 
           in base 2 (though any base will do) */

        if (knmodp == 0)
          score -= 0.5 * params->test_prime_contrib[j];
        else
          score -= params->test_prime_contrib[j];
      }
    }
    if (score < best_score) {
      best_score = score;
      best_mult = curr_mult;
    }
  }

  /* from now on we will factor best_mult * n */

  params->multiplier = best_mult;
  mpz_mul_ui(params->n, params->n, best_mult);
}

/***********************************/
func_static s32 init_fb_tiny(s32 fb_size)
/***********************************/
{
  tiny_qs_params *params = g_params;
  u16 *prime_list = params->prime_list;
  s32 i, j, mult_idx;
  s32 shift = 20;

  i = MIN_FB_OFFSET_TINY;
  mult_idx = 0;

  gprimes[i] = 2;
  params->multiplier_fb[0] = 0;
  params->multiplier_fb[1] = 0;

  /* Keep setting up factor base primes until enough
     are found or the pool of primes runs out */

#ifdef USE_8X_TD
  for (i++, j = 1; i < MIN_FB_OFFSET_TO_SIEVE_TINY && j < NUM_PRIMES_TINY; j++) {
#else
  for (i++, j = 1; i < fb_size && j < NUM_PRIMES_TINY; j++) {
#endif
    u8 logp;
    
    s32 prime = prime_list[j];
    s32 nmodp = mpz_tdiv_ui(params->n, prime);

    if (legendre_16(nmodp, prime) != -1) {
      logp = (u8)(LOGPRIME_SCALE_TINY *
                             log((double)prime) / M_LN2 + 0.5);
      
      /* if the prime divides n, it is part of n's 
         multiplier and is treated separately */

      if (nmodp != 0) {
        gmodsqrt[i] = (u16)sqrtModP_16(nmodp, prime);
      }
      else {
        gmodsqrt[i] = DO_NOT_SIEVE_TINY;
        params->multiplier_fb[mult_idx++] = i;
      }
      
      gprimes[i] = prime;
      glogprimes[i] = logp;
      grecip[i] = (u32)(B(32) / (u64)prime);      

      i++;
    }
  }
  
#ifdef USE_8X_TD
  for ( ; i < fb_size && j < NUM_PRIMES_TINY; j++) {
    u8 logp;
    
    s32 prime = prime_list[j];
    s32 nmodp = mpz_tdiv_ui(params->n, prime);

    if (legendre_16(nmodp, prime) != -1) {
      logp = (u8)(LOGPRIME_SCALE_TINY *
                             log((double)prime) / M_LN2 + 0.5);
      
      /* if the prime divides n, it is part of n's 
         multiplier and is treated separately */

      if (nmodp != 0) {
        gmodsqrt[i] = (u16)sqrtModP_16(nmodp, prime);
      }
      else {
        gmodsqrt[i] = DO_NOT_SIEVE_TINY;
        params->multiplier_fb[mult_idx++] = i;
      }
      
      gprimes[i] = prime;
      glogprimes[i] = logp;  
      grecip16[i] = (u16)(((u32)1 << shift) / prime);

      i++;
    }
  }

#endif

  params->fb_size = i;
  return i;
}

#if defined( __INTEL_COMPILER) && defined(notdefined)
static int presieve_block(void)
{
  int i, j;
  tiny_qs_params *params = g_params;
  int fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  u16 prime;
  u16 offset1;
  u16 offset2;
  u64 *ptab;
  __m256i p1;
  __m256i p2;
  __m256i init;

  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
    prime = gprimes[i];

    if (prime > 23) 
      break;
    
    offset1 = (roots1[i] * 32) % prime;
    offset2 = (roots2[i] * 32) % prime;
    
    for (j = 0; j < SIEVE_SIZE_TINY / 32; j++, offset1++, offset2++)
    {
      if (offset1 == prime) offset1 = 0;
      if (offset2 == prime) offset2 = 0;
      ptab = smallp_tabs[prime];
    
      init = _mm256_loadu_si256((__m256i *)(&sieve_block[j * 32]));
      p1 = _mm256_load_si256((__m256i *)(&ptab[offset1 * 4]));
      p2 = _mm256_load_si256((__m256i *)(&ptab[offset2 * 4]));
      
      p1 = _mm256_sub_epi8(init, p1);
      p1 = _mm256_sub_epi8(p1, p2);
      
      _mm256_storeu_si256((__m256i *)(&sieve_block[j * 32]), p1);
    }
  }
  return i;  
}
#endif

/***********************************/
func_static void fill_sieve_block_tiny(int starti)
/***********************************
Core sieving routine
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  

  /* Note that since this code will only ever 
     factor small inputs, the sieve interval will
     always be ridiculously small and does not 
     need to be broken up into chunks. Further,
     the bottleneck with small inputs is the trial
     factoring of relations and not the sieving,
     so no crazy unrolling tricks are needed
     here either. */
  // BRB: however, ordering the roots and unrolling x2 does still help.
    
     
    for (i = starti; i < fb_size; i++) {
      s32 prime = gprimes[i];
      u8 logprime = glogprimes[i];
      s32 root1;
      s32 root2;
      
      root1 = roots1[i];
      root2 = roots2[i];
      
      while (root2 < SIEVE_SIZE_TINY) {
        sieve_block[root1] -= logprime;
        sieve_block[root2] -= logprime;
        root1 += prime;
        root2 += prime;
      }

      if (root1 < SIEVE_SIZE_TINY)
        sieve_block[root1] -= logprime;
      
    }
}

/***********************************/
func_static void fill_sieve_block(int starti)
/***********************************
Core sieving routine
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  

  /* Note that since this code will only ever 
     factor small inputs, the sieve interval will
     always be ridiculously small and does not 
     need to be broken up into chunks. Further,
     the bottleneck with small inputs is the trial
     factoring of relations and not the sieving,
     so no crazy unrolling tricks are needed
     here either. */
  // BRB: however, ordering the roots and unrolling x2 does still help.
    
     
    for (i = starti; i < fb_size; i++) {
      s32 prime = gprimes[i];
      u8 logprime = glogprimes[i];
      s32 root1;
      s32 root2;
      
      root1 = roots1[i];
      root2 = roots2[i];
      
      while (root2 < SIEVE_SIZE) {
        sieve_block[root1] -= logprime;
        sieve_block[root2] -= logprime;
        root1 += prime;
        root2 += prime;
      }

      if (root1 < SIEVE_SIZE)
        sieve_block[root1] -= logprime;
      
    }
}


#ifdef USE_ASM

#define PACKED_SIEVE_MASK ((u64)0x80808080 << 32 | 0x80808080)

#define SM_SIEVE_SCAN_64_VEC					\
		asm volatile (							\
			"vmovdqu (%1), %%ymm0   \n\t"		\
			"vpor 32(%1), %%ymm0, %%ymm0    \n\t"		\
			"vpmovmskb %%ymm0, %%r11   \n\t"	/* output results to 64 bit register */		\
			"testq %%r11, %%r11 \n\t"			/* AND, and set ZF */ \
			"jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"vmovdqu (%1), %%ymm0   \n\t"		/* else, we had hits, move sections of sieveblock back in */ \
			"vmovdqu 32(%1), %%ymm1   \n\t"		/* extract high bit masks from each byte */ \
			"vpmovmskb %%ymm1, %%r9d   \n\t"		/*  */		\
			"salq $32, %%r9		\n\t"			/*  */ \
			"vpmovmskb %%ymm0, %%r8d   \n\t"		/*  */		\
			"orq	%%r9,%%r8		\n\t"		/* r8 now holds 64 byte mask results, in order, from sieveblock */ \
			"xorq	%%r11,%%r11		\n\t"		/* initialize count of set bits */ \
			"xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
			"1:			\n\t"					/* top of bit scan loop */ \
			"bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
			"addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
			"movb	%%r10b, (%2, %%r11, 1) \n\t"		/* put the bit index into the output buffer */ \
			"shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
			"incq	%%r11		\n\t"			/* increment the count of set bits */ \
            "incq	%%r10		\n\t"			/* increment the index */ \
			"shrq	$1, %%r8 \n\t"				/* clear the bit */ \
			"testq	%%r8,%%r8	\n\t"			/* check if there are any more set bits */ \
			"jnz 1b		\n\t"					/* loop if so */ \
			"2:		\n\t"						/*  */ \
			"movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "=r"(result)						\
			: "r"(packed_sieve_block + i), "r"(buffer)	\
			: "xmm0", "xmm1", "r8", "r9", "r10", "r11", "rcx", "cc", "memory");
      
#else

#define SM_SIEVE_SCAN_64_VEC

#endif


/***********************************/
func_static s32 mark_sieve_block_tiny(void)
/***********************************
Walk through a filled-in sieve block and find 
the offsets correspodning to relations that
are probably useful
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, k;
  u8 *sieve_block = params->sieve_block;
  u64 *packed_sieve_block = (u64 *)params->sieve_block;
     

  /* standard technique for testing sieve locations
     in parallel: initialize each byte to the target
     sieve value, and subtract logs of the factor base
     primes instead of adding them. Sieve offsets that
     accumulate enough log values become negative, 
     and it's easy to simultaneously test for the top 
     bit in several bytes being set */

  for (i = j = 0; i < SIEVE_SIZE_TINY / 8; i += 8) {

    /* handle 64 bytes at a time */
    u32 result = 0;
    u8 buffer[SIEVE_BATCH_SIZE_TINY];
		
		SM_SIEVE_SCAN_64_VEC;

		if (result == 0)
			continue;

    /* at least one byte is a hit; go back and search
       the list one at a time. We treat the sieve interval
       as a hashtable, and associate entry j in the list
       of relations to be resieved (params->sieve_batch[])
       with a byte that is negative. The high-order bit of
       the byte is set to indicate that the low-order bits 
       mean something */

    for (k=0; k<result; k++) {    
      u32 thisloc = (i * 8) + (u32)buffer[k];
      u32 val = sieve_block[thisloc];
      if (val & 0x80) {
        if (j < SIEVE_BATCH_SIZE_TINY) {
          tiny_relation *r = params->sieve_batch + j;
          r->sieve_offset = thisloc;
          r->num_factors = 0;
#ifdef DO_RESIEVE
          sieve_block[thisloc] = j | 0x80;
#endif
          j++;
        }
        else
        {
          sieve_block[thisloc] = 0;
        }
      }
    }
  }
  
#ifdef VERBOSE
  printf("block sieve: %u survivors\n", j);
#endif

  return j;
}

/***********************************/
func_static s32 mark_sieve_block(void)
/***********************************
Walk through a filled-in sieve block and find 
the offsets correspodning to relations that
are probably useful
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, k;
  u8 *sieve_block = params->sieve_block;
  u64 *packed_sieve_block = (u64 *)params->sieve_block;
     

  /* standard technique for testing sieve locations
     in parallel: initialize each byte to the target
     sieve value, and subtract logs of the factor base
     primes instead of adding them. Sieve offsets that
     accumulate enough log values become negative, 
     and it's easy to simultaneously test for the top 
     bit in several bytes being set */

  for (i = j = 0; i < SIEVE_SIZE / 8; i += 8) {

    /* handle 64 bytes at a time */
    u32 result = 0;
    u8 buffer[SIEVE_BATCH_SIZE_TINY];
		
		SM_SIEVE_SCAN_64_VEC;

		if (result == 0)
			continue;

    /* at least one byte is a hit; go back and search
       the list one at a time. We treat the sieve interval
       as a hashtable, and associate entry j in the list
       of relations to be resieved (params->sieve_batch[])
       with a byte that is negative. The high-order bit of
       the byte is set to indicate that the low-order bits 
       mean something */

    for (k=0; k<result; k++) {    
      u32 thisloc = (i * 8) + (u32)buffer[k];
      u32 val = sieve_block[thisloc];
      if (val & 0x80) {
          tiny_relation *r = params->sieve_batch + j;
          r->sieve_offset = thisloc;
          r->num_factors = 0;
#ifdef DO_RESIEVE
          sieve_block[thisloc] = j | 0x80;
#endif
          j++;
      }
    }
  }
  
#ifdef VERBOSE
  printf("block sieve: %u survivors\n", j);
#endif

  return j;
}


/***********************************/
func_static void resieve_tiny()
/***********************************
Just like fill_sieve_block_tiny(), except
sieving is used to avoid trial division
on all the relations previously found
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;

  /* Note that even though this routine does only
     a little more work than fill_sieve_block_tiny(),
     it runs almost 3x slower */

  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
    s32 prime = gprimes[i];
    s32 root1 = roots1[i];
    s32 root2 = roots2[i];

    while (root2 < SIEVE_SIZE_TINY) {
      s32 val1 = sieve_block[root1];
      s32 val2 = sieve_block[root2];
      if (val1 & 0x80) {
        tiny_relation *r = params->sieve_batch + (val1 & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
      if (val2 & 0x80) {
        tiny_relation *r = params->sieve_batch + (val2 & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
      root1 += prime;
      root2 += prime;
    }

    if (root1 < SIEVE_SIZE_TINY) {
      s32 val1 = sieve_block[root1];
      if (val1 & 0x80) {
        tiny_relation *r = params->sieve_batch + (val1 & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
    }
  }
}

/***********************************/
func_static void resieve()
/***********************************
Just like fill_sieve_block_tiny(), except
sieving is used to avoid trial division
on all the relations previously found
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;

  /* Note that even though this routine does only
     a little more work than fill_sieve_block_tiny(),
     it runs almost 3x slower */

  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
    s32 prime = gprimes[i];
    s32 root1 = roots1[i];
    s32 root2 = roots2[i];

    while (root2 < SIEVE_SIZE) {
      s32 val1 = sieve_block[root1];
      s32 val2 = sieve_block[root2];
      if (val1 & 0x80) {
        tiny_relation *r = params->sieve_batch + (val1 & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
      if (val2 & 0x80) {
        tiny_relation *r = params->sieve_batch + (val2 & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
      root1 += prime;
      root2 += prime;
    }

    if (root1 < SIEVE_SIZE) {
      s32 val1 = sieve_block[root1];
      if (val1 & 0x80) {
        tiny_relation *r = params->sieve_batch + (val1 & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
    }
  }
}

#ifdef USE_ASM

#define MOD_INIT_8X												\
		__asm__ (														\
      "vmovdqu (%0), %%xmm1 \n\t"		/* move in block_loc */ \
			:														\
			: "r" (bl_locs)		\
			: "xmm1");

#define MOD_INIT_16X												\
		__asm__ (														\
      "vmovdqu (%0), %%ymm1 \n\t"		/* move in block_loc */ \
			:														\
			: "r" (bl_locs)		\
			: "xmm1");

#define MOD_CMP_8X_vec(xtra_bits)																		\
		__asm__ (																				\
			"vmovdqu (%1), %%xmm3 \n\t"		/* move in primes */							\
      "vmovdqa %%xmm1, %%xmm4 \n\t" /* copy sieve_offset */ \
			"vmovdqu (%3), %%xmm6 \n\t"		/* move in root1s */							\
			"vpmulhuw	(%2), %%xmm4, %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
			"vmovdqu (%4), %%xmm2 \n\t"		/* move in root2s */							\
			"vpsrlw	$" xtra_bits ", %%xmm4, %%xmm4 \n\t"		/* to get to total shift of 20 bits */			\
			"vpmullw	%%xmm3, %%xmm4, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
			"vpsubw	%%xmm4, %%xmm1, %%xmm4 \n\t"	/* substract from blockloc */						\
			"vpcmpeqw	%%xmm4, %%xmm6, %%xmm6 \n\t"	/* compare to root1s */						\
			"vpcmpeqw	%%xmm4, %%xmm2, %%xmm2 \n\t"	/* compare to root2s */						\
			"vpor	%%xmm6, %%xmm2, %%xmm2 \n\t"	/* combine compares */							\
			"vpmovmskb %%xmm2, %%r8 \n\t"		/* export to result */							\
      "andl   $0xaaaaaaaa,%%r8d   \n\t"   /* mask the bits we don't care about */ \
      "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
      "xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
      "1:			\n\t"					/* top of bit scan loop */ \
      "bsfl	%%r8d,%%ecx		\n\t"		/* put least significant set bit index into rcx */ \
      "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
      "addl	%%ecx,%%r10d	\n\t"			/* add in the offset of this index */ \
      "movl   %%r10d,%%r9d \n\t" \
      "shrl   $1,%%r9d \n\t"   \
      "addl   %6,%%r9d \n\t"   \
      "movw	%%r9w, (%5, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
      "shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
      "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
      "incq	%%r10		\n\t"			/* increment the index */ \
      "shrq	$1, %%r8 \n\t"				/* clear the bit */ \
      "jmp 1b		\n\t"					/* loop if so */ \
      "2:		\n\t"						/*  */ \
      "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "+r" (numdiv)																	\
			: "r" (gprimes + i), \
				"r" (grecip16 + i), "r" (roots1 + i), \
			  "r" (roots2 + i), "r"(buffer), "r"(i) \
			: "r9", "r8", "r10", "r11", "rcx", "xmm2", "xmm3", "xmm4", "xmm6", "memory", "cc");
      
#define MOD_CMP_16X_vec(xtra_bits)																		\
		__asm__ (																				\
			"vmovdqu (%1), %%ymm3 \n\t"		/* move in primes */							\
      "vmovdqa %%ymm1, %%ymm4 \n\t" /* copy sieve_offset */ \
        "vmovdqu (%3), %%ymm6 \n\t"		/* move in root1s */							\
        "vpmulhuw	(%2), %%ymm4, %%ymm4 \n\t"	/* (unsigned) multiply by inverses */		\
        "vmovdqu (%4), %%ymm2 \n\t"		/* move in root2s */							\
        "vpsrlw	$" xtra_bits ", %%ymm4, %%ymm4 \n\t"		/* to get to total shift of 20 bits */			\
        "vpmullw	%%ymm3, %%ymm4, %%ymm4 \n\t"	/* (signed) multiply by primes */				\
        "vpsubw	%%ymm4, %%ymm1, %%ymm4 \n\t"	/* substract from blockloc */						\
        "vpcmpeqw	%%ymm4, %%ymm6, %%ymm6 \n\t"	/* compare to root1s */						\
        "vpcmpeqw	%%ymm4, %%ymm2, %%ymm2 \n\t"	/* compare to root2s */						\
        "vpor	%%ymm6, %%ymm2, %%ymm2 \n\t"	/* combine compares */							\
        "vpmovmskb %%ymm2, %%r8 \n\t"		/* export to result */							\
      "andl   $0xaaaaaaaa,%%r8d   \n\t"   /* mask the bits we don't care about */ \
      "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
      "xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
      "1:			\n\t"					/* top of bit scan loop */ \
      "bsfl	%%r8d,%%ecx		\n\t"		/* put least significant set bit index into rcx */ \
      "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
      "addl	%%ecx,%%r10d	\n\t"			/* add in the offset of this index */ \
      "movl   %%r10d,%%r9d \n\t" \
      "shrl   $1,%%r9d \n\t"   \
      "addl   %6,%%r9d \n\t"   \
      "movw	%%r9w, (%5, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
      "shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
      "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
      "incq	%%r10		\n\t"			/* increment the index */ \
      "shrq	$1, %%r8 \n\t"				/* clear the bit */ \
      "jmp 1b		\n\t"					/* loop if so */ \
      "2:		\n\t"						/*  */ \
      "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "+r" (numdiv)																	\
			: "r" (gprimes + i), \
				"r" (grecip16 + i), "r" (roots1 + i), \
			  "r" (roots2 + i), "r"(buffer), "r"(i) \
			: "r9", "r8", "r10", "r11", "rcx", "ymm2", "ymm3", "ymm4", "ymm6", "memory", "cc");

#else

#define MOD_INIT_8X
#define MOD_INIT_16X
#define MOD_CMP_8X_vec(xtra_bits)
#define MOD_CMP_16X_vec(xtra_bits)

#endif



/***********************************/
func_static s32 check_sieve_val_tiny(mpz_t a, mpz_t b, mpz_t c, 
                                 tiny_relation *r,
				 s32 sign_of_index)
/***********************************
Trial factor a relation that survived sieving
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  s32 num_factors = 0;
  s32 sieve_offset = r->sieve_offset;
  tiny_relation *relation = params->relation_list +
                            params->num_full_relations;
  u16 *fb_offsets = relation->fb_offsets;
  static u8 initialized = 0;
  static mpz_t res, res2;
  
  
    
#ifndef DO_RESIEVE
  u8 logp = 0;
  
#ifdef USE_8X_TD

  ALIGNED_MEM u16 bl_locs[16];
  
  bl_locs[0] = sieve_offset;
  bl_locs[1] = sieve_offset;
  bl_locs[2] = sieve_offset;
  bl_locs[3] = sieve_offset;
  bl_locs[4] = sieve_offset;
  bl_locs[5] = sieve_offset;
  bl_locs[6] = sieve_offset;
  bl_locs[7] = sieve_offset;  
  
#ifdef USE_16X_TD
  bl_locs[8] = sieve_offset;
  bl_locs[9] = sieve_offset;
  bl_locs[10] = sieve_offset;
  bl_locs[11] = sieve_offset;
  bl_locs[12] = sieve_offset;
  bl_locs[13] = sieve_offset;
  bl_locs[14] = sieve_offset;
  bl_locs[15] = sieve_offset;  
#endif
  
#endif
#endif

  if (initialized == 0) {
    mpz_init(res);
    mpz_init(res2);
    initialized = 1;
  }

  /* form the polynomial value */

  mpz_mul_ui(res, a, sieve_offset);
  if (sign_of_index == POSITIVE)
    mpz_add(res, res, b);
  else
    mpz_sub(res, res, b);
  mpz_mul_ui(res, res, sieve_offset);
  mpz_add(res, res, c);
  if (mpz_sgn(res) < 0) {
    mpz_abs(res, res);
    fb_offsets[num_factors++] = 0;
  }

  /* extract powers of two */

  i = mpz_scan1(res, 0);
  if (i) {
    mpz_tdiv_q_2exp(res, res, i);
    do {
      if (num_factors >= MAX_FACTORS_TINY)
          return 0;
      fb_offsets[num_factors++] = MIN_FB_OFFSET_TINY;
#ifndef DO_RESIEVE
      logp++;
#endif
    } while (--i);
  }

  /* divide out the unsieved factor base primes */
  for (i = MIN_FB_OFFSET_TINY + 1; 
         i < MIN_FB_OFFSET_TO_SIEVE_TINY; i++) {
    
    s32 prime = gprimes[i];
    s32 root1 = roots1[i];
    s32 root2 = roots2[i];
    u32 recip = grecip[i];
    
    if (root1 == DO_NOT_SIEVE_TINY)
      continue;

    j = (s32)(((u64)sieve_offset * (u64)recip) >> 32);
    j = sieve_offset - j * prime;
    if (j >= prime)
      j -= prime;

    if (j == root1 || j == root2) {
      while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;
        fb_offsets[num_factors++] = i;
#ifndef DO_RESIEVE
        logp += glogprimes[i];
#endif
        mpz_swap(res, res2);
      }
    }
  }

  /* divide out the factors of the multiplier, 
     if any */
     
  for (i = 0; i < 2; i++) {
    if (params->multiplier_fb[i]) {
      s32 prime;
      j = params->multiplier_fb[i];
      prime = gprimes[j];

      while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;
        fb_offsets[num_factors++] = j;
#ifndef DO_RESIEVE
        logp += glogprimes[j];
#endif
        mpz_swap(res, res2);
      }
    }
  }

#ifdef DO_RESIEVE
  /* We should probably have been adding log values
     to the log of this relation in the previous loops,
     and testing that the complete log value now
     exceeds the trial factoring cutoff. However, 
     resieving has already found the remaining factors, 
     so we wouldn't save much time bailing out at 
     this point */

  for (i = 0; i < r->num_factors; i++) {
    s32 prime;
    j = r->fb_offsets[i];
    prime = gprimes[j];

    while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
      if (num_factors >= MAX_FACTORS_TINY)
          return 0;
      fb_offsets[num_factors++] = j;
      mpz_swap(res, res2);
    }
  }
#else
  
  // bail if we didn't find enough small factors
  if (logp < 6)
    return 0;
  
  #ifdef USE_8X_TD
  {
    s32 numdiv = 0;
    ALIGNED_MEM u16 buffer[32];   
    
    MOD_INIT_8X;
    
    // trial divide by the rest of the factor base 8 at a time using SSE41.
    // with blocksizes of 4096 or 16384 a constant shift value of 20 
    // seems to always work for primes > a few bits, which we should always
    // have with MIN_FB_OFFSET_TO_SIEVE_TINY == 7
    for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; 
           i < g_params->fb_size; i += 8) {
      MOD_CMP_8X_vec("4");
    }
    
    for (j = 0; j < numdiv; j++) {
      s32 prime;
      
      if (buffer[j] >= g_params->fb_size)
        break;
      
      if (roots1[buffer[j]] == DO_NOT_SIEVE_TINY)
        continue;
      
      prime = gprimes[buffer[j]];      
      while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;
        fb_offsets[num_factors++] = buffer[j];
        mpz_swap(res, res2);
      }
    }
  }

  #else
  
  // trial divide by the rest of the factor base
  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; 
         i < g_params->fb_size; i++) {
    
    s32 prime = gprimes[i];
    s32 root1 = roots1[i];
    s32 root2 = roots2[i];
    u32 recip = grecip[i];
    
    if (root1 == DO_NOT_SIEVE_TINY)
      continue;

    j = (s32)(((u64)sieve_offset * (u64)recip) >> 32);
    j = sieve_offset - j * prime;
    if (j >= prime)
      j -= prime;

    if (j == root1 || j == root2) {
      while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;
        fb_offsets[num_factors++] = i;
        mpz_swap(res, res2);
      }
    }
  }
  
  #endif
  
#endif

  /* start filling in the final relation */

  if (sign_of_index == NEGATIVE)
    sieve_offset = -sieve_offset;
  relation->sieve_offset = sieve_offset;
  relation->num_factors = num_factors;
  relation->poly_num = params->poly_num;

  if (mpz_cmp_ui(res, 1) == 0) {

    /* full relation; we're done */

    relation->large_prime = 1;
    params->num_full_relations++;
  }
  else if (mpz_cmp_ui(res, params->large_prime_max) < 0) {
    s32 lp = mpz_get_ui(res);
    u32 table_idx = LARGE_PRIME_HASH(lp);
    s32 partial_idx;

    /* partial relation; see if it has occurred already */

    relation->large_prime = lp;
    for (i = 0; i < LP_HASH_DEPTH_TINY; i++) {
      partial_idx = params->partial_hash[table_idx][i];
      if (partial_idx == 0xffff ||
          lp == params->relation_list[partial_idx].large_prime)
        break;
    }

    if (i == LP_HASH_DEPTH_TINY) {

      /* not found, and no room to store it */

      return 0;
    }
    else if (partial_idx == 0xffff) {

      /* not found, but the hashtable entry has
         room to keep it; transfer the relation to
         the partial list */

      params->relation_list[params->partial_idx] = *relation;
      params->partial_hash[table_idx][i] = params->partial_idx--;
    }
    else {

      /* large prime has matched, new relation can stay */

      params->num_full_relations++;
    }
  }

  /* make sure the 'heap' of full relations has not
     overflowed into the 'stack' of partial relations */

  if (params->num_full_relations >= params->partial_idx)
    return -1;
  return 0;
}

/***********************************/
func_static void init_siqs_tiny(void)
/***********************************
Initialize the subsystem for forming SIQS
sieve polynomials
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  s32 plus_idx, minus_idx;
  u32 fb_size = params->fb_size;
  u32 num_factors = params->num_a_factors;

  /* compute the optimal size of the factors of
     the polynomial 'A' value. We know how many
     primes it should have, and know the optimal
     A value that will minimize sieving time. 
     Assume further that all factors are the
     same size. First compute the factor size,
     then locate the factor base offset where
     it approximately occurs */

  j = (s32)(exp(log(params->target_a) / num_factors) + 0.5);
  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size - 1; i++) {
    if (gprimes[i] > j)
      break;
  }
  if (i == MIN_FB_OFFSET_TINY + 1)
    i++;

  /* polynomial A values are built by selecting from
     a pool of primes. There are POLY_SELECT_BITS_TINY
     primes in the pool, evenly distributed above and
     below the optimal factor base offset */

  memset(params->poly_select_offsets, 0,
         sizeof(params->poly_select_offsets));
  plus_idx = i;
  minus_idx = i-1;
  i = 0;
  
  while (1) {
    if (plus_idx < fb_size && 
        gmodsqrt[plus_idx] != DO_NOT_SIEVE_TINY) {
      params->poly_select_offsets[i] = plus_idx;
      if (++i == POLY_SELECT_BITS_TINY)
        break;
    }
    if (minus_idx > MIN_FB_OFFSET_TINY + 1 &&
        gmodsqrt[minus_idx] != DO_NOT_SIEVE_TINY) {
      params->poly_select_offsets[i] = minus_idx;
      if (++i == POLY_SELECT_BITS_TINY)
        break;
    }
    plus_idx++;
    minus_idx--;
  }
  
  /* polynomial selection will begin at offset
     zero of the tables below */

  params->poly_select_idx = 0;
}

/* A perpetual problem with SIQS is deciding which
   primes should make up the next polynomial A value.
   The selected set must multiply out to a value as
   close as possible to the optimal A value, but must
   be sufficiently different from previously selected
   sets that the odds of producing duplicate relations
   are low. And the set has to be computed quickly.

   Fortunately, we will only need a few polynomials
   so the sets to use can be precomputed. Each prime in
   the pool is assigned a bit in a bitfield. Consecutive
   bits in the bitfield refer to primes alternately 
   above and below the optimal factor size. The low-order
   bits correspond to factors near the optimal value, and
   the more significant bits march away from the optimal
   value. Hence, setting the bitfield to an integer 
   will select a unique set of primes, the number of which
   is the number of set bits in the integer. Small integer
   values of the bitfield will pick primes close to the
   optimal factor size, with later bitfield values selecting
   prime factors that march away from the optimal size.

   popcount[] gives the number of set bits for each value
   of the bitfield, and a_choice[] lists the bitfields
   themselves. A given factorization only uses one of the
   population count sizes from the table; bitfields are
   arranged so that low-order bits are set first, then
   higher-order bits are set. Every bitfield is different
   by at least two bits from all other bitfields with the
   same weight, and there are enough bitfields to generate
   256 polynomials, whether A values contain 3, 4, or 5 primes */
   
u8 popcount[] = {
       3,     4,     3,     5,     3,     4,
       3,     4,     3,     3,     4,     4,     
       3,     4,     5,     4,     5,     4,     
       4,     4,     4,     5,     5,     4,     
       4,     5,     5,     4,     5,     5,     
       5,     5,     3,     5,     5,     5,     
       5,     5,     3,     4,     3,     4,     
       4,     4,     3,     3,     4,     4,     
       4,     4,     3,     4,     4,     4,     
       4,     3,     4,     4,     3,     4,     
       4,     4,     4,     3,     4,     3,     
};

u16 a_choice[] = {
       0x007, 0x00f, 0x019, 0x01f, 0x02a, 0x033,
       0x034, 0x03c, 0x04c, 0x052, 0x055, 0x05a,
       0x061, 0x066, 0x067, 0x069, 0x079, 0x096,
       0x099, 0x0a5, 0x0aa, 0x0ab, 0x0b5, 0x0c3,
       0x0cc, 0x0cd, 0x0d3, 0x0f0, 0x12d, 0x133,
       0x14b, 0x155, 0x181, 0x187, 0x199, 0x1e1,
       0x22e, 0x256, 0x282, 0x303, 0x304, 0x30c,
       0x330, 0x3c0, 0x484, 0x502, 0x505, 0x50a,
       0x550, 0x5a0, 0x601, 0x606, 0x609, 0x660,
       0x690, 0x888, 0x906, 0x909, 0x910, 0x960,
       0x990, 0xa05, 0xa0a, 0xa20, 0xa50, 0xc40,
};


/***********************************/
func_static s32 find_poly_a(mpz_t a)
/***********************************
Compute the next polynomial A value
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, mask;
  s32 num_a_factors = params->num_a_factors;
  tiny_poly *poly = params->poly_list + params->poly_num;

  /* choose the next bitfield representing
     primes to use */

  for (i = params->poly_select_idx; i < sizeof(popcount); i++) {
    if (popcount[i] == num_a_factors)
      break;
  }
  if (i >= sizeof(popcount))
    return -1;
  mask = a_choice[i];
  params->poly_select_idx = i + 1;

  /* gather the chosen primes */

  for (i = j = 0; i < POLY_SELECT_BITS_TINY; i++) {
    if (!(mask & (1 << i)))
      continue;

    if (params->poly_select_offsets[i] == 0)
      return -2;

    poly->a_fb_offsets[j] = params->poly_select_offsets[i];
    if (++j == num_a_factors)
      break;
  }

  /* multiply them together */

  mpz_set_ui(a, 1);
  for (i = 0; i < num_a_factors; i++) {
    j = poly->a_fb_offsets[i];
    mpz_mul_ui(a, a, gprimes[j]);
  }

  return 0;
}

/***********************************/
func_static void find_first_poly_b(mpz_t a, mpz_t b, mpz_t c)
/***********************************
Compute the first of a list of polynomial
B values
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  s32 num_a_factors = params->num_a_factors;
  u32 fb_size = params->fb_size;
  tiny_poly *poly = params->poly_list + params->poly_num;

  mpz_set_ui(b, 0);

  /* fill in the auxiliary quantities needed to
     compute future B values */

  for (i = 0; i < num_a_factors; i++) {
    s32 g, prime = gprimes[poly->a_fb_offsets[i]]; 

    mpz_divexact_ui(params->poly_b_aux[i], a, prime);
    g = mpz_tdiv_ui(params->poly_b_aux[i], prime);
    g = modinv_16(g, prime);
    g = (s32)g * gmodsqrt[poly->a_fb_offsets[i]] % prime; 
    if (g > prime/2)
      g = prime - g;
    mpz_mul_ui(params->poly_b_aux[i], 
               params->poly_b_aux[i], g);
    mpz_add(b, b, params->poly_b_aux[i]);
    mpz_add(params->poly_b_aux[i],
            params->poly_b_aux[i],
            params->poly_b_aux[i]);
  }

  /* This first B is the sum of the auxiliary 
     quantities computed previously */

  mpz_set(poly->b, b);
  
  /* Form C, a helper for computing the value
     of a polynomial before trial factoring */

  mpz_mul(c, b, b);
  mpz_sub(c, c, params->n);
  mpz_divexact(c, c, a);

  /* initialize the factor base for sieving */

  // order these roots
  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size; i++) {    
    s32 prime = gprimes[i];
    s32 modsqrt = gmodsqrt[i];
    s32 amodp = mpz_tdiv_ui(a, prime);
    s32 bmodp = prime - mpz_tdiv_ui(b, prime);
    u16 root1;
    u16 root2;
    u16 tmp;
    
    if (modsqrt == DO_NOT_SIEVE_TINY) {

      /* factors of the multiplier never 
         contribute to sieving */
      root1 = DO_NOT_SIEVE_TINY;
      root2 = DO_NOT_SIEVE_TINY;

    }
    else if (amodp == 0) {

      /* factor base primes that divide the A value
         get one sieve root and not two */

      amodp = prime - mpz_tdiv_ui(c, prime);
      root1 = amodp * modinv_16(2 * bmodp % prime, prime) % prime;
      root2 = DO_NOT_SIEVE_TINY;
      
    }
    else {      
      
      /* handle all the other FB primes, including the
         initialization that allows the next 2^(num_a_factors-1)-1
         factor bases to initialize quickly */

      amodp = modinv_16(amodp, prime);
      root1 = amodp * (bmodp + modsqrt) % prime;
      root2 = amodp * (bmodp + prime - modsqrt) % prime;
            
      if (root1 > root2)
      {
        tmp = root1;
        root1 = root2;
        root2 = tmp;    
      }    
      
      for (j = 0; j < num_a_factors; j++) {
        bmodp = mpz_tdiv_ui(params->poly_b_aux[j], prime);
        root_aux[j * fb_size + i] = bmodp * amodp % prime;
      }
    }
    
    roots1[i] = root1;
    roots2[i] = root2;

  }
}
      

#ifdef USE_ASM
      
#define COMPUTE_16X_SMALL_PROOTS_AVX2	\
		asm (	\
			"movl   $2, %%eax \n\t"	/* eax = start */	\
			"vpxor	%%ymm6, %%ymm6, %%ymm6 \n\t" \
			"cmpl	%0, %%eax \n\t"	\
			"jge	1f \n\t"	\ 
			"0: \n\t"	\
			/* compute 16 new roots on the P side */	\
			"vmovdqu    (%4, %%rax, 2), %%ymm3 \n\t"			/* xmm3 = next 16 updates */	\
			"vmovdqu    (%1, %%rax, 2), %%ymm1 \n\t"			/* xmm1 = next 16 values of root1 */	\
			"vmovdqu    (%2, %%rax, 2), %%ymm2 \n\t"			/* xmm2 = next 16 values of root2 */	\
      "vmovdqu    (%3, %%rax, 2), %%ymm0 \n\t"			/* ymm0 = next 16 primes */	\
      "vpcmpgtw	%%ymm6, %%ymm1, %%ymm9 \n\t"				/* root1 > 0 ? 0xffff : 0 */ \
      "vpcmpgtw	%%ymm6, %%ymm2, %%ymm10 \n\t"				/* root2 > 0 ? 0xffff : 0 */ \
      "vpand	    %%ymm0,%%ymm9,  %%ymm9 \n\t"					/* copy prime to where roots > 0 */	\
			"vpand	    %%ymm0,%%ymm10,  %%ymm10 \n\t"					/* copy prime to where roots > 0 */	\
      "vpsubw	    %%ymm1, %%ymm9, %%ymm1 \n\t"					/* selective root1 = prime - root1 */	\
			"vpsubw	    %%ymm2, %%ymm10, %%ymm2 \n\t"					/* selective root2 = prime - root2 */	\
			"vmovdqa	%%ymm1, %%ymm4 \n\t"					/* copy r1 */ \
			"vmovdqa	%%ymm2, %%ymm5 \n\t"					/* copy r2 */ \
			"vmovdqa	%%ymm1, %%ymm7 \n\t"					/* copy r1 */ \
			"vmovdqa	%%ymm2, %%ymm8 \n\t"					/* copy r2 */ \
			"vpsubw	    %%ymm3, %%ymm1, %%ymm1 \n\t"					/* root1 -= ptr */	\
			"vpsubw	    %%ymm3, %%ymm2, %%ymm2 \n\t"					/* root2 -= ptr */	\
			"vpsubusw	%%ymm3, %%ymm4, %%ymm4 \n\t"					/* root1 -= ptr */	\
			"vpsubusw	%%ymm3, %%ymm5, %%ymm5 \n\t"					/* root2 -= ptr */	\
			"vpcmpeqw	%%ymm3, %%ymm7, %%ymm7 \n\t"				/* xmm4 := root1 == ptr ? 0xffff : 0 */ \
			"vpcmpeqw	%%ymm3, %%ymm8, %%ymm8 \n\t"				/* xmm5 := root2 == ptr ? 0xffff : 0 */ \
			"vpcmpeqw	%%ymm6, %%ymm4, %%ymm4 \n\t"				/* xmm4 := ptr >= root1 ? 0xffff : 0 */ \
			"vpcmpeqw	%%ymm6, %%ymm5, %%ymm5 \n\t"				/* xmm5 := ptr >= root2 ? 0xffff : 0 */ \
			"vpandn	    %%ymm4, %%ymm7, %%ymm7 \n\t"					/* ptr > root1 (greater than and not equal) */	\
			"vpandn	    %%ymm5, %%ymm8, %%ymm8 \n\t"					/* ptr > root1 (greater than and not equal) */	\
			"vpand	    %%ymm0, %%ymm7, %%ymm7 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpand	    %%ymm0, %%ymm8, %%ymm8 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpaddw	    %%ymm7, %%ymm1, %%ymm1 \n\t"					/* selectively add back prime (modular sub) */	\
			"vpaddw	    %%ymm8, %%ymm2, %%ymm2 \n\t"					/* selectively add back prime (modular sub) */	\
			"vmovdqa    %%ymm2, %%ymm5 \n\t"					/* xmm5 = root2 copy */	\
			"vpmaxuw	%%ymm1, %%ymm5, %%ymm5 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"vpminuw	%%ymm1, %%ymm2, %%ymm2 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
															/* now copy results to appropriate data structures */	\
															/* root1p always gets the smaller roots (LT) */	\
			"vmovdqu	%%ymm2, (%1, %%rax, 2) \n\t"			/* update root1p */	\
			"vmovdqu	%%ymm5, (%2, %%rax, 2) \n\t"			/* update root2p */	\
			"addl	$16, %%eax \n\t"	\
			"cmpl	%0, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "r"(fb_size), "r"(roots1), "r"(roots2), "r"(gprimes), "r"(row)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "rax", "cc", "memory");

  
#define COMPUTE_16X_SMALL_NROOTS_AVX2	\
		asm (	\
			"movl   $2, %%eax \n\t"	/* eax = start */	\
			"vpxor	%%ymm6, %%ymm6, %%ymm6 \n\t" \
			"cmpl	%0, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 16 new roots on the N side */	\
			"vmovdqu    (%4, %%rax, 2), %%ymm3 \n\t"			/* xmm3 = next 16 updates */	\
			"vmovdqu    (%1, %%rax, 2), %%ymm1 \n\t"			/* xmm1 = next 16 values of root1 */	\
			"vmovdqu    (%2, %%rax, 2), %%ymm2 \n\t"			/* xmm2 = next 16 values of root2 */	\
      "vmovdqu    (%3, %%rax, 2), %%ymm0 \n\t"			/* ymm0 = next 16 primes */	\
      "vpcmpgtw	%%ymm6, %%ymm1, %%ymm9 \n\t"				/* root1 > 0 ? 0xffff : 0 */ \
      "vpcmpgtw	%%ymm6, %%ymm2, %%ymm10 \n\t"				/* root2 > 0 ? 0xffff : 0 */ \
      "vpand	    %%ymm0,%%ymm9,  %%ymm9 \n\t"					/* copy prime to where roots > 0 */	\
			"vpand	    %%ymm0,%%ymm10,  %%ymm10 \n\t"					/* copy prime to where roots > 0 */	\
      "vpsubw	    %%ymm1, %%ymm9, %%ymm1 \n\t"					/* selective root1 = prime - root1 */	\
			"vpsubw	    %%ymm2, %%ymm10, %%ymm2 \n\t"					/* selective root2 = prime - root2 */	\
			"vmovdqa	  %%ymm1, %%ymm7 \n\t"					/* copy root1 */ \
			"vmovdqa	  %%ymm2, %%ymm8 \n\t"					/* copy root2 */ \
			"vpaddw	    %%ymm3, %%ymm1, %%ymm1 \n\t"					/* root1 += ptr */	\
			"vpaddw	    %%ymm3, %%ymm2, %%ymm2 \n\t"					/* root2 += ptr */	\
			"vpaddusw	  %%ymm3, %%ymm7, %%ymm7 \n\t"					/* root1 += ptr */	\
			"vpaddusw	  %%ymm3, %%ymm8, %%ymm8 \n\t"					/* root2 += ptr */	\
			"vmovdqa	  %%ymm0, %%ymm4 \n\t"					/* copy primes before comparison */ \
			"vmovdqa	  %%ymm0, %%ymm5 \n\t"					/* copy primes before comparison */ \
			"vpsubusw	  %%ymm7, %%ymm4, %%ymm4 \n\t"				/* ymm4 := prime - root1 (saturated subs to set masks for modadd) */ \
			"vpsubusw	  %%ymm8, %%ymm5, %%ymm5 \n\t"				/* ymm5 := prime - root2 (saturated subs to set masks for modadd) */ \
			"vpcmpeqw	  %%ymm6, %%ymm4, %%ymm4 \n\t"				/* ymm4 := root1 >= prime ? 1 : 0 */ \
			"vpcmpeqw	  %%ymm6, %%ymm5, %%ymm5 \n\t"				/* ymm5 := root2 >= prime ? 1 : 0 */ \
			"vpand	    %%ymm0, %%ymm4, %%ymm4 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpand	    %%ymm0, %%ymm5, %%ymm5 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpsubw	    %%ymm4, %%ymm1, %%ymm1 \n\t"					/* selectively sub back prime (modular add) */	\
			"vpsubw	    %%ymm5, %%ymm2, %%ymm2 \n\t"					/* selectively sub back prime (modular add) */	\
			"vmovdqa    %%ymm2, %%ymm5 \n\t"					/* ymm5 = root2 copy */	\
			"vpmaxuw	  %%ymm1, %%ymm5, %%ymm5 \n\t"					/* ymm5 = root2 > root1 ? root2 : root1 */	\
			"vpminuw	  %%ymm1, %%ymm2, %%ymm2 \n\t"					/* ymm2 = root2 < root1 ? root2 : root1 */	\
															/* now copy results to appropriate data structures */	\
															/* root1p always gets the smaller roots (LT) */	\
															/* root2p always gets the bigger roots (GT) */	\
      "vmovdqu	%%ymm2, (%1, %%rax, 2) \n\t"			/* update root1p */	\
			"vmovdqu	%%ymm5, (%2, %%rax, 2) \n\t"			/* update root2p */	\
			"addl	$16, %%eax \n\t"	\
			"cmpl	%0, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "r"(fb_size), "r"(roots1), "r"(roots2), "r"(gprimes), "r"(row)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "rax", "cc", "memory");
      
      

#define FLIP_ROOTS_AVX2	\
		asm (	\
			"movl   $2, %%eax \n\t"	/* eax = start */	\
			"vpxor	%%ymm6, %%ymm6, %%ymm6 \n\t" \
			"cmpl	%0, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* flip the roots 16 at a time */	\
			"vmovdqu    (%1, %%rax, 2), %%ymm1 \n\t"			/* xmm1 = next 16 values of root1 */	\
			"vmovdqu    (%2, %%rax, 2), %%ymm2 \n\t"			/* xmm2 = next 16 values of root2 */	\
      "vmovdqu    (%3, %%rax, 2), %%ymm0 \n\t"			/* ymm0 = next 16 primes */	\
      "vpcmpgtw	%%ymm6, %%ymm1, %%ymm9 \n\t"				/* root1 > 0 ? 0xffff : 0 */ \
      "vpcmpgtw	%%ymm6, %%ymm2, %%ymm10 \n\t"				/* root2 > 0 ? 0xffff : 0 */ \
      "vpcmpeqw	%%ymm9, %%ymm1, %%ymm7 \n\t"				/* root1 == 0xffff ? 0xffff : 0 */ \
      "vpcmpeqw	%%ymm10, %%ymm2, %%ymm8 \n\t"				/* root2 == 0xffff ? 0xffff : 0 */ \
      "vpandn	    %%ymm9, %%ymm7,  %%ymm9 \n\t"					/* root1 > 0 && roots != DO_NOT_SIEVE_TINY */	\
			"vpandn	    %%ymm10, %%ymm8,  %%ymm10 \n\t"					/* root1 > 0 && roots != DO_NOT_SIEVE_TINY */	\
      "vpand	    %%ymm0,%%ymm9,  %%ymm9 \n\t"					/* copy prime to where roots > 0 && roots != DO_NOT_SIEVE_TINY */	\
			"vpand	    %%ymm0,%%ymm10,  %%ymm10 \n\t"					/* copy prime to where roots > 0 && roots != DO_NOT_SIEVE_TINY */	\
      "vpsubw	    %%ymm1, %%ymm9, %%ymm1 \n\t"					/* selective root1 = prime - root1 */	\
			"vpsubw	    %%ymm2, %%ymm10, %%ymm2 \n\t"					/* selective root2 = prime - root2 */	\
      "vmovdqa    %%ymm2, %%ymm5 \n\t"					/* ymm5 = root2 copy */	\
			"vpmaxuw	  %%ymm1, %%ymm5, %%ymm5 \n\t"					/* ymm5 = root2 > root1 ? root2 : root1 */	\
			"vpminuw	  %%ymm1, %%ymm2, %%ymm2 \n\t"					/* ymm2 = root2 < root1 ? root2 : root1 */	\
															/* now copy results to appropriate data structures */	\
															/* root1p always gets the smaller roots (LT) */	\
															/* root2p always gets the bigger roots (GT) */	\
      "vmovdqu	%%ymm2, (%1, %%rax, 2) \n\t"			/* update root1p */	\
			"vmovdqu	%%ymm5, (%2, %%rax, 2) \n\t"			/* update root2p */	\
			"addl	$16, %%eax \n\t"	\
			"cmpl	%0, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "r"(fb_size), "r"(roots1), "r"(roots2), "r"(gprimes)	\
			: "xmm0", "xmm1", "xmm2", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10","rax", "cc", "memory");
      
#else
#define COMPUTE_16X_SMALL_PROOTS_AVX2
#define COMPUTE_16X_SMALL_NROOTS_AVX2
#define FLIP_ROOTS_AVX2
#endif
      
/***********************************/
func_static void find_next_poly_b(mpz_t a, mpz_t b, mpz_t c)
/***********************************
Initialize B values beyond the first
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  s32 num_a_factors = params->num_a_factors;
  s32 fb_size = params->fb_size;
  tiny_poly *poly = params->poly_list + params->poly_num;
  u32 mask = params->poly_num & ((1 << (num_a_factors-1)) - 1);
  u8 do_sub;
  u16 *row;
  /* current poly starts of with the previous poly */

  mpz_set(b, poly[-1].b);
  for (i = 0; i < num_a_factors; i++)
    poly[0].a_fb_offsets[i] = poly[-1].a_fb_offsets[i];

  /* determine the auxiliary B constant that comes next
     in Gray code order, and add to or subtract from
     the current B. This also determines which of the
     rows from the table of corrections are applied to
     the factor base */

  i = 0;
  while ((mask & (1 << i)) == 0)
    i++;

  row = root_aux + fb_size * i;

  do_sub = 0;
  if (mask & (1 << (i+1))) {
    // this corresponds to the the gray > 0, sign > 0, COMPUTE_NEXT_ROOTS_P case in yafu
    mpz_add(b, b, params->poly_b_aux[i]);
    do_sub = 1;
  }
  else {
    // this corresponds to the the gray < 0, sign < 0, COMPUTE_NEXT_ROOTS_N case in yafu
    mpz_sub(b, b, params->poly_b_aux[i]);
  }

  /* form the C helper value */

  mpz_mul(c, b, b);
  mpz_sub(c, c, params->n);
  mpz_divexact(c, c, a);

  /* set up the factor base for the next B */
  if (do_sub)
  {
    COMPUTE_16X_SMALL_PROOTS_AVX2;
  }
  else
  {
    COMPUTE_16X_SMALL_NROOTS_AVX2;
  }
  
  // now do the special cases
  if (params->multiplier_fb[0] > 0) {
    roots1[params->multiplier_fb[0]] = DO_NOT_SIEVE_TINY; 
    roots2[params->multiplier_fb[0]] = DO_NOT_SIEVE_TINY; 
  }
  if (params->multiplier_fb[1] > 0) {
    roots1[params->multiplier_fb[1]] = DO_NOT_SIEVE_TINY;
    roots2[params->multiplier_fb[1]] = DO_NOT_SIEVE_TINY;
  }
  
  for (i = 0; i < num_a_factors; i++)
  {
      s32 prime = gprimes[poly->a_fb_offsets[i]];

      /* sieving with root1 but not root 2 only
         happens if the prime divides 'A'. Compute
         the new sieve root manually */

      s32 cmodp = prime - mpz_tdiv_ui(c, prime);
      s32 bmodp = mpz_tdiv_ui(b, prime);
      if (mpz_sgn(b) > 0)
        bmodp = prime - bmodp;

      roots1[poly->a_fb_offsets[i]] = cmodp * modinv_16(2 * bmodp % prime, prime) % prime;
      roots2[poly->a_fb_offsets[i]] = DO_NOT_SIEVE_TINY;
      
    }    
    
  mpz_set(poly[0].b, b);
}

/***********************************/
func_static s32 sieve_next_poly_tiny(void)
/***********************************
Do all the sieving for one polynomial
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  s32 cutoff1, num_surviving;
  s32 poly_num = params->poly_num;
  s32 target_relations = params->fb_size + NUM_EXTRA_RELATIONS_TINY;
  static u8 initialized = 0;
  static mpz_t a, b, c;
  u16 tmp;
  
  if (initialized == 0) {
    mpz_init(a); mpz_init(b); mpz_init(c);
    initialized = 1;
  }

  /* generate the polynomial */

  if (!(poly_num & ((1 << (params->num_a_factors-1))-1))) {
    i = find_poly_a(a);
    if (i)
      return i;
    find_first_poly_b(a, b, c);
  }
  else {
    find_next_poly_b(a, b, c);
  }

  /* compute the cutoff beyond which trial factoring
     will be used on sieve values. */

#ifdef VERBOSE
  //gmp_printf("c = %Zu (%d bits)\n", c, mpz_sizeinbase(c, 2));
#endif
  // bitsize of c doesn't really change much poly to poly.
  // it may be beneficial to fix the size of cutoff1, to make
  // other things more predicable.
  cutoff1 = LOGPRIME_SCALE_TINY * (mpz_sizeinbase(c, 2) - 
                  params->error_bits - SMALL_PRIME_FUDGE_TINY - 1);

  /* the trial factoring code wants 2*B and not B */

  mpz_add(b, b, b);
  
  /* sieve over positive offsets, mark the most
     promising offsets, resieve to trial factor
     them all at once and then finish each in turn */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE_TINY);
  fill_sieve_block_tiny(MIN_FB_OFFSET_TO_SIEVE_TINY);
  // currently the presieve is not effective... maybe if it subbed 
  // more than one prime per 32B sieve block like the soe version...
  //i = presieve_block();
  //fill_sieve_block_tiny(i);
  
  num_surviving = mark_sieve_block_tiny();
  if (num_surviving) {
#ifdef DO_RESIEVE
    resieve_tiny();
#endif
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(a, b, c, params->sieve_batch + i, POSITIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        return 0;
    }
  }
  
  //FLIP_ROOTS_AVX2;
      
  /* flip the sieve roots from positive to
     negative values */
  for (i = MIN_FB_OFFSET_TINY + 1;  i < fb_size; i++) {
    
    s32 prime = gprimes[i];
    s32 root1;
    s32 root2;
    
    root1 = roots1[i];
    root2 = roots2[i];

    if (root1 != DO_NOT_SIEVE_TINY && root1)
      root1 = prime - root1;
    if (root2 != DO_NOT_SIEVE_TINY && root2)
      root2 = prime - root2;
    
    if (root1 > root2)
      {
        tmp = root1;
        root1 = root2;
        root2 = tmp;       
      }

    roots1[i] = root1;
    roots2[i] = root2;
  }
  
    
  /* repeat the sieve procedure for negative
     sieve offsets */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE_TINY);
  fill_sieve_block_tiny(MIN_FB_OFFSET_TO_SIEVE_TINY);
  num_surviving = mark_sieve_block_tiny();
  if (num_surviving) {
#ifdef DO_RESIEVE
    resieve_tiny();
#endif
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(a, b, c, params->sieve_batch + i, NEGATIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        return 0;
    }
  }
      
  return 0;
}

/***********************************/
func_static s32 sieve_next_poly(void)
/***********************************
Do all the sieving for one polynomial
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  s32 cutoff1, num_surviving;
  s32 poly_num = params->poly_num;
  s32 target_relations = params->fb_size + NUM_EXTRA_RELATIONS_TINY;
  static u8 initialized = 0;
  static mpz_t a, b, c;
  u16 tmp;
  
  if (initialized == 0) {
    mpz_init(a); mpz_init(b); mpz_init(c);
    initialized = 1;
  }

  /* generate the polynomial */

  if (!(poly_num & ((1 << (params->num_a_factors-1))-1))) {
    i = find_poly_a(a);
    if (i)
      return i;
    find_first_poly_b(a, b, c);
  }
  else {
    find_next_poly_b(a, b, c);
  }

  /* compute the cutoff beyond which trial factoring
     will be used on sieve values. */

#ifdef VERBOSE
  //gmp_printf("c = %Zu (%d bits)\n", c, mpz_sizeinbase(c, 2));
#endif
  // bitsize of c doesn't really change much poly to poly.
  // it may be beneficial to fix the size of cutoff1, to make
  // other things more predicable.
  cutoff1 = LOGPRIME_SCALE_TINY * (mpz_sizeinbase(c, 2) - 
                  params->error_bits - SMALL_PRIME_FUDGE - 1);

  /* the trial factoring code wants 2*B and not B */

  mpz_add(b, b, b);
  
  /* sieve over positive offsets, mark the most
     promising offsets, resieve to trial factor
     them all at once and then finish each in turn */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE);
  fill_sieve_block(MIN_FB_OFFSET_TO_SIEVE_TINY);
  num_surviving = mark_sieve_block();
  if (num_surviving) {
#ifdef DO_RESIEVE
    resieve();
#endif
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(a, b, c, params->sieve_batch + i, POSITIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        return 0;
    }
  }
  
  //FLIP_ROOTS_AVX2;
      
  /* flip the sieve roots from positive to
     negative values */
  for (i = MIN_FB_OFFSET_TINY + 1;  i < fb_size; i++) {
    
    s32 prime = gprimes[i];
    s32 root1;
    s32 root2;
    
    root1 = roots1[i];
    root2 = roots2[i];

    if (root1 != DO_NOT_SIEVE_TINY && root1)
      root1 = prime - root1;
    if (root2 != DO_NOT_SIEVE_TINY && root2)
      root2 = prime - root2;
    
    if (root1 > root2)
      {
        tmp = root1;
        root1 = root2;
        root2 = tmp;       
      }

    roots1[i] = root1;
    roots2[i] = root2;
  }
  
    
  /* repeat the sieve procedure for negative
     sieve offsets */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE);
  fill_sieve_block(MIN_FB_OFFSET_TO_SIEVE_TINY);
  num_surviving = mark_sieve_block();
  if (num_surviving) {
#ifdef DO_RESIEVE
    resieve();
#endif
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(a, b, c, params->sieve_batch + i, NEGATIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        return 0;
    }
  }
      
  return 0;
}

/***********************************/
func_static void solve_linear_system_tiny(void)
/***********************************
Find linear dependencies among a set of relations
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, k, start_row;
  s32 nrows = params->fb_size;
  s32 ncols = params->num_full_relations;
  s32 num_a_factors = params->num_a_factors;
  u16 rowperm[MAX_FB_SIZE_TINY];
  u16 pivot[MAX_FB_SIZE_TINY];
  s32 row = 0;

  memset(params->matrix, 0, sizeof(params->matrix));

  /* build the matrix; relations become columns, and
     pairs of matched partial relations fuse into 
     columns as well */

  for (i = 0; i < ncols; i++) {
    tiny_relation *r;
    tiny_poly *poly;
    for (j = 0; j < 2; j++) {

      r = params->relation_list + i;
      if (j == 1) {
        s32 hash_idx = LARGE_PRIME_HASH(r->large_prime);
        s32 partial_idx;
        for (k = 0; k < LP_HASH_DEPTH_TINY; k++) {
          partial_idx = params->partial_hash[hash_idx][k];
          if (params->relation_list[partial_idx].large_prime == r->large_prime)
            break;
        }
        r = params->relation_list + partial_idx;
      }
      poly = params->poly_list + r->poly_num;

      for (k = 0; k < r->num_factors; k++) {
        row = r->fb_offsets[k];
        params->matrix[row][i / 64] ^= bitmask[i % 64];
      }

      /* the factors in the polynomial A value
         figure into the matrix as well */

      for (k = 0; k < num_a_factors; k++) {
        row = poly->a_fb_offsets[k];
        params->matrix[row][i / 64] ^= bitmask[i % 64];
      }
      if (r->large_prime == 1)
        break;
    }
  }
  for (i = 0; i < nrows; i++)
    rowperm[i] = i;

  /* begin with a random vector of dependencies */

  for (i = 0; i < ncols; i++)
    params->null_vectors[i] = (u16)get_rand(
                      &rand_seed1, &rand_seed2);

  /* perform the elimination */

  for (i = start_row = 0; start_row < nrows && i < ncols; i++) {
    
    /* find the next pivot */

    for (j = start_row; j < nrows; j++) {
      row = rowperm[j];
      if (params->matrix[row][i / 64] & bitmask[i % 64])
        break;
    }
    if (j == nrows)
      continue;

    rowperm[j] = rowperm[start_row];
    rowperm[start_row] = row;
    pivot[start_row++] = i;

    /* eliminate it from the other rows */

    for (j++; j < nrows; j++) {
      s32 row2 = rowperm[j];
      if (params->matrix[row2][i / 64] & bitmask[i % 64]) {
        for (k = i / 64; k < (ncols + 63) / 64; k++) {
          params->matrix[row2][k] ^= params->matrix[row][k];
        }
      }
    }
  }

  /* perform back substitution */

  for (i = start_row - 1; i >= 0; i--) {
    u16 accum;
    row = rowperm[i];

    for (j = pivot[i] + 1, accum = 0; j < ncols; j++) {
      if (params->matrix[row][j / 64] & bitmask[j & 63])
        accum ^= params->null_vectors[j];
    }
    params->null_vectors[pivot[i]] = accum;
  }
}


/***********************************/
func_static u32 find_factors_tiny(mpz_t factor1, 
                             mpz_t factor2)
/***********************************
perform MPQS square root phase
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, k;
  u16 mask;
  u16 fb_counts[MAX_FB_SIZE_TINY];
  static mpz_t x, y, t0, t1;
  static u8 initialized = 0;

  if (initialized == 0) {
    mpz_init(x); mpz_init(y);
    mpz_init(t0); mpz_init(t1);
    initialized = 1;
  }

  /* for each dependency */

  for (mask = 1; mask; mask <<= 1) {

    memset(fb_counts, 0, sizeof(fb_counts));
    mpz_set_ui(x, 1);
    mpz_set_ui(y, 1);

    /* for each relation allowed in the dependency */
    for (i = 0; i < params->num_full_relations; i++) {

      if (!(params->null_vectors[i] & mask))
        continue;

      for (j = 0; j < 2; j++) {
        tiny_relation *r = params->relation_list + i;
        tiny_poly *poly;

        /* match up partials with the same large prime */

        if (j == 1) {
          s32 hash_idx = LARGE_PRIME_HASH(r->large_prime);
          s32 partial_idx;
          for (k = 0; k < LP_HASH_DEPTH_TINY; k++) {
            partial_idx = params->partial_hash[hash_idx][k];
            if (params->relation_list[partial_idx].large_prime == 
                r->large_prime)
              break;
          }
          r = params->relation_list + partial_idx;
          mpz_mul_ui(t0, y, r->large_prime);
          mpz_mod(y, t0, params->n);
        }
        poly = params->poly_list + r->poly_num;
  
        /* add the factors of this relation to the table
           of factors. Include the factors of A as well */

        for (k = 0; k < r->num_factors; k++)
          fb_counts[r->fb_offsets[k]]++;
  
        mpz_set_ui(t1, 1);
        for (k = 0; k < params->num_a_factors; k++) {
          s32 idx = poly->a_fb_offsets[k];
          fb_counts[idx]++;
          mpz_mul_ui(t1, t1, gprimes[idx]);
        }

        /* multiply A * sieve_offset + B into the left 
           side of the congruence */

        if (r->sieve_offset < 0) {
          mpz_mul_ui(t1, t1, -(r->sieve_offset));
          mpz_sub(t1, t1, poly->b);
        }
        else {
          mpz_mul_ui(t1, t1, r->sieve_offset);
          mpz_add(t1, t1, poly->b);
        }
        mpz_mul(t0, x, t1);
        mpz_mod(x, t0, params->n);

        if (r->large_prime == 1)
          break;
      }
    }

    /* Form the right side of the congruence; given its
       prime factorization, cut the exponent of each prime
       in half and perform a modular exponentiation */
    for (i = MIN_FB_OFFSET_TINY; i < params->fb_size; i++) {
      u16 mask2 = 0x8000;
      u16 exponent = fb_counts[i] / 2;
      u32 prime = gprimes[i];
      
      if (exponent == 0)
        continue;

      mpz_set_ui(t0, prime);
      while (!(exponent & mask2))
        mask2 >>= 1;

      for (mask2 >>= 1; mask2; mask2 >>= 1) {
        mpz_mul(t1, t0, t0);
        mpz_mod(t0, t1, params->n);
        if (exponent & mask2) {
          mpz_mul_ui(t1, t0, prime);
          mpz_mod(t0, t1, params->n);
        }
      }
      mpz_mul(t1, t0, y);
      mpz_mod(y, t1, params->n);
    }

    /* For x and y the halves of the congruence, 
       compute gcd(x+-y, n) */
    for (i = 0; i < 2; i++) {
      if (i == 0)
        mpz_add(t0, x, y);
      else
        mpz_sub(t0, x, y);

      mpz_gcd(t1, t0, params->n);
      if (mpz_cmp_ui(t1, 1) && mpz_cmp(t1, params->n)) {

        /* we've possibly found a nontrivial factor of n.
           Divide any factors of the multiplier out from
           both factors */

        u32 mult1 = 0;
        u32 mult2 = 0;

        if (params->multiplier_fb[0])
          mult1 = gprimes[params->multiplier_fb[0]];
        if (params->multiplier_fb[1])
          mult2 = gprimes[params->multiplier_fb[1]];

        
        mpz_divexact(t0, params->n, t1);
        if (mult1) {
          if (mpz_tdiv_ui(t0, mult1) == 0)
            mpz_divexact_ui(t0, t0, mult1);
          if (mpz_tdiv_ui(t1, mult1) == 0)
            mpz_divexact_ui(t1, t1, mult1);
        }
        if (mult2) {
          if (mpz_tdiv_ui(t0, mult2) == 0)
            mpz_divexact_ui(t0, t0, mult2);
          if (mpz_tdiv_ui(t1, mult2) == 0)
            mpz_divexact_ui(t1, t1, mult2);
        }

        /* If both remaining factors exceed unity, 
           we've factored n and can stop */
        if (mpz_cmp_ui(t0, 1) && mpz_cmp_ui(t1, 1)) {
          mpz_set(factor1, t0);
          mpz_set(factor2, t1);
          return 1;
        }
      }
    }

    /* otherwise try the next dependency */
  }

  return 0;
}


typedef struct {
  s32 fb_size;
  s32 num_poly_factors;
} tiny_qs_config;

/* factor base sizes and num_poly_factors for 
   50 to 116-bit factorizations.
   the row lookup is: ((bits - 50) / 4);  */
tiny_qs_config static_config[] = {
 { 40, 3 },
 { 50, 3 },
 { 66, 3 },
 { 66, 3 },
 { 82, 3 },
 { 98, 3 },
 { 114, 3 },
 { 130, 3 },
 { 146, 4 },
 { 146, 4 },
 { 162, 4 },
 { 178, 4 },
 { 226, 4 },
 { 274, 5 },
 { 354, 5 },
 { 418, 5 },
 { 482, 5 },
};

/***********************************/
u32 tinyqs(mpz_t n, mpz_t factor1, mpz_t factor2)
/***********************************
Main driver for MPQS factorization
Returns 1 and sets factor1 and factor2 if
successful, returns 0 otherwise
************************************/
{
  tiny_qs_params *params;
  s32 bits, status = 0;
  s32 fb_size;
  s32 bound;
  s32 large_prime_mult;
  tiny_qs_config *config;
  static mpz_t tmp;
  static u8 initialized = 0;

  if (initialized == 0) {
    mpz_init(tmp);
    initialized = 1;
  }

  /* make sure the input isn't a perfect square.
     We may also want to add a test for a perfect
     cube, but that's so unlikely it's probably
     not worth worrying about */

  if (mpz_root(tmp, n, 2) != 0) {
    mpz_set(factor1, tmp);
    mpz_set(factor2, tmp);
    return 1;
  }

  /* start the initialization */

  init_tinyqs();
  
  params = g_params;
  mpz_set(params->n, n);
  params->num_full_relations = 0;
  params->partial_idx = 4 * MAX_RELATIONS_TINY - 1;
  params->poly_num = 0;

  bits = mpz_sizeinbase(params->n, 2);
  find_multiplier_tiny();

  /* determine the factor base size and the
     number of primes in a polynomial A value */

  if (bits < 50)
    bits = 50;
  if (bits > 116)
    bits = 116;
  config = static_config + ((bits - 50) / 4);
  fb_size = config->fb_size;
  params->num_a_factors = config->num_poly_factors;

  /* build the factor base */

  fb_size = init_fb_tiny(fb_size);

  /* compute the optimal A value */

  mpz_sqrt(tmp, params->n);
  if (bits > 85)
    params->target_a = mpz_get_d(tmp) * M_SQRT2 / SIEVE_SIZE;
  else
    params->target_a = mpz_get_d(tmp) * M_SQRT2 / SIEVE_SIZE_TINY;
  init_siqs_tiny();

  /* compute the large prime cutoff and the
     size of the fudge factor needed to account
     for it in the sieving cutoff */

  if (bits == 116)
    large_prime_mult = 50;
  else
    large_prime_mult = 15;
  bound = gprimes[fb_size - 1];
  bound *= large_prime_mult;
  params->large_prime_max = bound;
  params->error_bits = (u32)(log(bound) / M_LN2 + 1);

  /* empty out the hashtable for partial relations */

  memset(params->partial_hash, 0xff, sizeof(params->partial_hash));
  
#ifdef VERBOSE
  gmp_printf("input = %Zu\n", n);
  printf("fbsize = %u, largest prime = %u\n", fb_size, gprimes[fb_size - 1]);
  printf("large_prime_mult = %u\n", large_prime_mult);
  printf("config->num_poly_factors = %u\n", config->num_poly_factors);
  printf("first sieved primes = %u, %u, %u\n", gprimes[7], gprimes[8], gprimes[9]);
#endif

  /* do the sieving! */

  if (bits > 85)
  {
    while (params->poly_num < MAX_POLY_TINY &&
           params->num_full_relations < fb_size + NUM_EXTRA_RELATIONS_TINY) {
      if (sieve_next_poly() != 0)
        return 0;
      params->poly_num++;
    }
  }
  else
  {
    while (params->poly_num < MAX_POLY_TINY &&
           params->num_full_relations < fb_size + NUM_EXTRA_RELATIONS_TINY) {
      if (sieve_next_poly_tiny() != 0)
        return 0;
      params->poly_num++;
    }
  }
  
  /* if enough relations were found, finish
     off the factorization */

  if (params->num_full_relations >= fb_size + NUM_EXTRA_RELATIONS_TINY) {
      solve_linear_system_tiny();
      status = find_factors_tiny(factor1, factor2);
  }

  return status;
}



