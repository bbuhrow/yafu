#ifndef _SIQS_IMPL_H_
#define _SIQS_IMPL_H_


#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/stat.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <signal.h>
#include "ytools.h"
#include "common.h"
#include "monty.h"
#include "cofactorize.h"
#include "gmp.h"
#include "factor.h"
#include "msieve_common.h"
#include "savefile.h"

// I've been unable to get this to run any faster, for
// either generic or knc codebases...
//#define USE_BATCHPOLY
//#define USE_BATCHPOLY_X2

// to use subset-sum searching, define this...
//#define USE_SS_SEARCH 

// and either this one
//#define USE_DIRECT_SIEVE_SS

// or this one, but not both
//#define USE_POLY_BUCKET_SS

#ifdef USE_SS_SEARCH
#define USE_POLY_BUCKET_PN_COMBINED_VARIATION

// an element of a poly-bucket is a 32-bit integer divided
// into two pieces of information.  The high X bits store an
// index into the factor base starting from the current slice
// fboffset.  Smaller slice bits means more slices are needed
// for a given number of primes treated by subset-sum.  
// The bottom 32-X bits are for storing the sieve hit.  Thus
// this field determines the maximum sieve region size.  
// Note, too many slices can be detrimental to performance, so
// it is a balancing act between allowing large sieve regions and
// the right amount of factor-base slices.  Further, different
// combinations of these may be best for differently sized inputs.
// 12 slice bits and 20 bits for sieve hits (max sieve interval
// of 16 blocks per side: 16 * 2 * 32768 = 2^20) seems to be
// a good general-purpose combination.
#define SS_SLICE_SIZE 4096      // these must agree
#define SS_SLICE_BITS 12        // these must agree
#define SS_SIGN_BIT (31 - SS_SLICE_BITS)
#define SS_MAX_ROOT (1 << (SS_SIGN_BIT))
#define SS_ROOT_MASK (SS_MAX_ROOT - 1)
#endif

// harebrained idea
#define NUM_ALP 3

// as part of analyzing 3lp parameterizations, we save off the 
// full residues in tdiv.  These will be fully factored later and sorted
// so that a synthetic data set resembling this real one can
// be constructed for deeper analysis of cycle generation for
// this parameterization.
//#define GATHER_RESIDUE_STATS

#ifdef _MSC_VER
// optionally define this or not depending on whether your hardware supports it.
// if defined, compile the sse41 functions into the fat binary.  the global
// flag HAS_SSE41 is set at runtime on compatible hardware to enable the functions
// to be used.  For gcc and mingw64 builds, USE_SSE41 is enabled in the makefile.
#ifdef USE_BATCHPOLY
#define FORCE_GENERIC 1
#else
#define USE_SSE41 1
#endif

//#define USE_AVX2 1
#endif

// assume we have SSE2 available
#define D_HAS_SSE2

#if (defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)) || \
    (defined(_MSC_VER) && defined(__clang__))
#define _MM_SCALE_1 1
#define _MM_SCALE_2 2
#endif

#if defined(USE_AVX2)
#ifdef _MSC_VER
#define CLEAN_AVX2 _mm256_zeroupper();
#else
#define CLEAN_AVX2 __asm__ volatile ("vzeroupper   \n\t");
#endif
#endif

#define USE_BATCH_FACTOR

#ifdef USE_BATCH_FACTOR
#include "batch_factor.h"
#endif

/************************* Common types and functions *****************/

#define BLOCKBITS 15
#define BLOCKSIZEm1 32767
#define BLOCKSIZE 32768

// store only largish smooth factors of relations; trial divide during filtering
#define SPARSE_STORE

#ifdef SPARSE_STORE
#define MAX_SMOOTH_PRIMES 200	//maximum number of factors for a smooth, including duplicates
#else
#define MAX_SMOOTH_PRIMES 200	//maximum number of factors for a smooth, including duplicates
#endif

#ifdef USE_DIRECT_SIEVE_SS
#define MAX_SIEVE_REPORTS 65536
#else
#define MAX_SIEVE_REPORTS 2048
#endif
#define MIN_FB_OFFSET 1
#define NUM_EXTRA_QS_RELATIONS 64
#define MAX_A_FACTORS 20
#define BUCKET_ALLOC 2048
#define BUCKET_BITS 11
#define BUCKET_ALLOCtxt "2048"
#define HALFBUCKET_ALLOCtxt "1024"
#define BUCKET_BITStxt "11"


// always use these optimizations using sse2
#if !defined (FORCE_GENERIC)
#define USE_8X_MOD_ASM 1
#define USE_RESIEVING
#define SSE2_RESIEVING 1
#endif

// multiplication by inverse constants
#define FOGSHIFT 24
#define FOGSHIFT_2 40

#if !defined (FORCE_GENERIC) && (defined (__MINGW64__) || (defined(__GNUC__) && defined(__x86_64__) && defined(GCC_ASM64X)))
    // requires 64 bit and gcc inline assembler syntax
#define USE_POLY_SSE2_ASM 1
#define USE_ASM_SMALL_PRIME_SIEVING
#endif

#if !defined (FORCE_GENERIC) && (defined(GCC_ASM64X) || defined(__MINGW64__) || (defined(_MSC_VER) && defined(USE_AVX2)))
    //assume we have sse2, set defines to use the sse2 code available
#define SIMD_SIEVE_SCAN 1
#define SIMD_SIEVE_SCAN_VEC 1

#elif !defined (FORCE_GENERIC) && (defined(GCC_ASM32X) || defined(__MINGW32__))
#define SIMD_SIEVE_SCAN 1

#elif !defined (FORCE_GENERIC) && defined(MSC_ASM32A)
#define SIMD_SIEVE_SCAN 1

#elif !defined (FORCE_GENERIC) && defined(_WIN64)
#define SIMD_SIEVE_SCAN 1
#endif




typedef struct
{
    uint32_t small_limit;
    uint32_t num_blocks;
    uint32_t large_mult;
    double fudge_factor;
    uint32_t num_extra_relations;
    uint32_t dlp_lower;
    uint32_t dlp_upper;
    int in_mem;
    int use_dlp;
} qs_params;

/*--------------LINEAR ALGEBRA RELATED DECLARATIONS ---------------------*/

/* Used to represent a list of relations */

typedef struct {
    uint32_t num_relations;  /* number of relations in the cycle */
    uint32_t* list;          /* list of offsets into an array of relations */
} qs_la_cycle_t;

/* A column of the matrix */

typedef struct {
    uint32_t* data;		/* The list of occupied rows in this column */
    uint32_t weight;		/* Number of nonzero entries in this column */
    qs_la_cycle_t cycle;       /* list of relations comprising this column */
} qs_la_col_t;


//holds only the info necessary during the sieve.  removing unnecessary info allows us to fit more
//factor base elements into cache during the sieve
typedef struct
{
    uint16_t* prime;
    uint16_t* root1;
    uint16_t* root2;
#ifdef LOGP_BITS
    uint8_t* logp;
#else
    uint16_t* logp;
#endif
} sieve_fb_compressed;

/************************* SIQS types and functions *****************/
#define MAXLP 4

typedef struct
{
    uint32_t large_prime[MAXLP];		//large prime in the pd.
    uint32_t sieve_offset;		//offset specifying Q (the quadratic polynomial)
    uint32_t apoly_idx;           // the index of the a-poly this relation uses
    uint32_t poly_idx;			//which b-poly this relation uses
    uint32_t parity;				//the sign of the offset (x) 0 is positive, 1 is negative
    uint32_t fb_offsets[MAX_SMOOTH_PRIMES];			//offsets of factor base primes dividing Q(offset).  
                                //note that other code limits the max # of fb primes to < 2^16
    uint32_t num_factors;			//number of factor base factors in the factorization of Q
} siqs_r;

typedef struct poly_t {
    uint32_t a_idx;				// offset into a list of 'a' values 
    mpz_t b;					// the MPQS 'b' value 
} poly_t;

typedef struct
{
    uint32_t* correction;
    uint32_t* small_inv;
    uint32_t* prime;
    uint32_t* logprime;
} tiny_fb_element_siqs;

typedef struct
{
#ifdef USE_8X_MOD_ASM
    uint16_t* correction;
    uint16_t* small_inv;
#else
    uint32_t* correction;
    uint32_t* small_inv;
#endif
    uint32_t* prime;
    uint32_t* logprime;
    uint64_t* binv;
} fb_element_siqs;

// we trade a few bit operations per prime 
// to be able to load both the fb_index and the sieve location in a single 
// 32bit load.  A bucket element is now stored using a uint32_t.  The top 16 bits 
// hold the fb_index and the bottom 16 bits hold the sieve hit index (loc)
#ifdef USEBUCKETSTRUCT
//used when bucket sieving
typedef struct
{
    uint16_t fb_index;	//a bucket only holds primes from a slice of the FB, must be < 2^16 primes
    uint16_t loc;			//blocksize must never exceed 2^16
} bucket_element;
#endif

typedef struct
{
    uint32_t* prime;
    int* firstroots1;
    int* firstroots2;
    uint16_t* sm_firstroots1;
    uint16_t* sm_firstroots2;
    uint8_t* logp;
} update_t;

typedef struct
{
    uint32_t* num;			    // array of counts in each bucket (1000 max buckets)
    uint32_t* fb_bounds;		// array of boundaries in the factorbase
    uint8_t* logp;			    // array of the logp values in each bucket
    uint32_t num_slices;		// the number of fb slices needed
    uint32_t alloc_slices;	    // the number of fb slices allocated
    uint32_t *num_slices_batch;	// the number of fb slices needed
    uint32_t list_size;		    // number of contiguous buckets allocated
    uint32_t* list;			    // contiguous space for all buckets
} lp_bucket;

typedef struct
{
    uint32_t* list;			// contiguous space for all xl-bucket elements
    uint32_t* slicenum;     // number of elements in each slice of the fb.
    uint32_t* sliceid;      // starting fb-index for the slice.
    uint8_t* slicelogp;		//  logp value for the slice
    uint32_t numslices;     // number of slices.   
    uint32_t alloc_slices;	// the number of fb slices allocated
} xlp_bucket_t;

typedef struct
{
    int* root;
    int* polynum;
    int size;
    int alloc;
} ss_set_t;

typedef struct
{
    uint32_t* elements;     // one huge list for all buckets
    uint32_t* size;         // number of filled elements in each bucket (a portion of the huge array).
    uint32_t alloc;         // the portion in the huge array that is allocated to each poly
    uint32_t numbuckets;    // number of buckets, need one for each b-poly
    uint32_t fboffset;      // starting fb offset for this slice
    uint8_t logp;           // logp for this bucket.
} ss_bucket_slice_t;

typedef struct
{
    mpz_t mpz_poly_a;
    mpz_t mpz_poly_b;
    mpz_t mpz_poly_c;
    int* qlisort;	//sorted list of indices into the factor base of factors of poly_a
    char* gray;		//2^s - 1 length array of signs for the graycode
    char* nu;		//2^s - 1 length array of Bl values to use in the graycode
    int s;
    int index;      // index within the master list of poly_a's
} siqs_poly;

typedef struct
{
    uint32_t B;					//number of primes in the entire factor base
    uint32_t small_B;				//index 1024
    uint32_t fb_10bit_B;			//index at which primes are bigger than 10 bits (and a multiple of 8)
    uint32_t fb_11bit_B;			//index at which primes are bigger than 11 bits (and a multiple of 8)
    uint32_t fb_12bit_B;			//index at which primes are bigger than 12 bits (and a multiple of 8)
    uint32_t fb_13bit_B;			//index at which primes are bigger than 13 bits (and a multiple of 8)
    uint32_t fb_32k_div3;
    uint32_t fb_14bit_B;			//index at which primes are bigger than 14 bits (and a multiple of 8)
    uint32_t fb_15bit_B;			//index at which primes are bigger than 15 bits (and a multiple of 8)
    uint32_t med_B;				//index at which primes are bigger than blocksize (and a multiple of 8)
    uint32_t large_B;				//index at which primes are bigger than entire sieve interval (and a multiple of 8)
    uint32_t x2_large_B;
    uint32_t ss_start_B;
    uint32_t num_ss_slices;
    uint32_t slice_size;
    tiny_fb_element_siqs* tinylist;
    fb_element_siqs* list;
} fb_list;

/* The sieving code needs a special structure for determining
   the number of cycles in a collection of partial relations */

#define QS_LOG2_CYCLE_HASH 22

typedef struct {
    uint32_t next;
    uint32_t prime;
    uint32_t data;
    uint32_t count;
} qs_cycle_t;

typedef struct {
    fact_obj_t* obj;			// passed in with info from 'outside'

    // the stuff in this structure is written once, then is read only
    // for the duration of the factorization, thus can be shared between
    // threads

    uint32_t digits_n;			// digits in n * multiplier
    uint32_t bits;				// bits in n * multiplier

    uint32_t* modsqrt_array;		// a square root of n mod each FB prime
    uint32_t multiplier;			// small multiplier for n (may be composite) 
    uint32_t knmod8;
    mpz_t n;					// the number to factor (scaled by multiplier)
    mpz_t sqrt_n;				// sqrt of n
    fb_list* factor_base;       // the factor base to use
    uint32_t* sieve_primes;            // sieve primes
    uint32_t num_sieve_blocks;	// number of sieve blocks in sieving interval
    uint32_t sieve_small_fb_start;// starting FB offset for sieving
    mpz_t target_a;				// optimal value of 'a' 

    double fudge_factor;		// used to compute the closnuf value
    uint32_t large_mult;			// used to compute the large prime cutoff
    uint32_t num_blocks;			// number of blocks to sieve on each side
    uint32_t num_extra_relations;	// number of extra relations to find
    uint32_t small_limit;			// upper limit of small prime variation
    uint32_t use_dlp;				// use double large primes? (0 or 1)
    uint32_t dlp_lower;			// lower bit range for dlp factorization attempts
    uint32_t dlp_upper;			// upper bit range for dlp factorization attempts

    uint32_t sieve_interval;		// one side of the sieve interval
    uint32_t qs_blocksize;		// blocksize of the sieve - only to be used
    uint32_t qs_blockbits;		// in non speed critical areas of the code	
    uint8_t blockinit;			// initial sieve value

    uint32_t pmax;				// largest prime in factor base
    int scan_unrolling;			// how many bytes to unroll the sieve scan

    uint32_t tf_small_cutoff;		// bit level to determine whether to bail early from tf
    uint32_t tf_closnuf;			// subject anything sieved beyond this to tf

    uint32_t tf_small_recip2_cutoff;
    uint32_t tf_med_recip1_cutoff;
    uint32_t tf_med_recip2_cutoff;
    uint32_t tf_large_cutoff;
    //uint32_t poly_batch_size;     // how many polys to update at one time?

    mpz_t prime_product;


    struct timeval totaltime_start;	// start time of this job

    uint32_t large_prime_max;		// the cutoff value for keeping a partial
                                    // relation; actual value, not a multiplier
    uint64_t max_fb2;				// the square of the largest factor base prime 
    uint64_t large_prime_max2;		// the cutoff value for factoring dlp-partials 
    double dlp_exp;

    double max_fb3;					// the cube of the largest factor base prime 
    double large_prime_max3;		// the cutoff value for factoring tlp-partials 
    double tlp_exp;                 // large_prime_max3 = large_prime_max ^ tlp_exp

    double max_fb4;					// largest factor base prime ^ 4
    double large_prime_max4;		// the cutoff value for factoring qlp-partials 
    double qlp_exp;                 // large_prime_max4 = large_prime_max ^ qlp_exp

    //master list of cycles
    qs_cycle_t* cycle_table;		/* list of all the vertices in the graph */
    uint32_t cycle_table_size;	/* number of vertices filled in the table */
    uint32_t cycle_table_alloc;	/* number of cycle_t structures allocated */
    uint32_t* cycle_hashtable;	/* hashtable to index into cycle_table */
    uint32_t components;			/* connected components (see relation.c) */
    uint32_t vertices;			/* vertices in graph (see relation.c) */
    uint32_t num_relations;		/* number of relations in list */
    uint32_t num_cycles;			/* number of cycles in list */
    uint32_t num_r;				// total relations found
    uint64_t num;					// sieve locations we've subjected to trial division
    uint32_t num_found;
    uint32_t num_slp;
    uint32_t num_full;
    uint32_t num_lp;

    //used to check on progress of the factorization	
    struct timeval update_start;	// time at which we last assessed the situation
    double t_update;			// ?
    double update_time;			// time at which we should next assess the situation
    uint32_t check_inc;			// for really short jobs, these kick in before the
    uint32_t check_total;			// timers can go off.  they use the number of
    uint32_t last_numfull;		// relations found since the last update as a guide
    uint32_t last_numpartial;		// to when to assess the situation
    uint32_t last_numcycles;		// used in computing rels/sec
    double last_fullrate;		// used in tlp to predict relations/batch
    double last_cycrate;		// used in tlp to predict relations/batch
    uint32_t num_needed;
    uint32_t num_expected;
    int charcount;				// characters on the screen
    uint32_t tot_poly;
    uint32_t failed_squfof;
    uint32_t attempted_squfof;
    uint32_t failed_cosiqs;
    uint32_t attempted_cosiqs;
    uint32_t dlp_outside_range;
    uint32_t dlp_prp;
    uint32_t dlp_useful;
    uint32_t tlp_outside_range;
    uint32_t tlp_prp;
    uint32_t tlp_useful;
    uint32_t qlp_outside_range;
    uint32_t qlp_prp;
    uint32_t qlp_useful;
    uint64_t total_reports;
    uint64_t total_surviving_reports;
    uint64_t total_blocks;
    uint32_t lp_scan_failures;
    uint64_t num_scatter_opp;
    uint64_t num_scatter;

    //master time record
    double t_time1;				// sieve time
    double t_time2;				// relation scanning and trial division
    double t_time3;				// polynomial calculations and large prime sieving
    double t_time4;				// extra?

    //stuff for sqrt
    factor_list_t factor_list;

    //these are used during linear algebra and sqrt root
    uint32_t total_poly_a;		// used during sieving: total number of polynomial 'a' values 
    mpz_t* poly_a_list;			// used during sieving: list of 'a' values for MPQS polys 
    uint32_t filt_total_poly_a;	// used during filtering: total number of polynomial 'a' values 
    mpz_t* filt_poly_a_list;	// used during filtering: list of 'a' values for MPQS polys 
    poly_t* poly_list;			// list of MPQS polynomials 
    uint32_t poly_list_alloc;
    uint32_t apoly_alloc;

    //these are used during filtering
    mpz_t curr_a;	  			// the current 'a' value in filtering
    mpz_t* curr_b;				// list of all the 'b' values for that 'a' 
    uint32_t bpoly_alloc;			// size of allocated curr_b
    siqs_r* relation_list;		// list of relations	
    qs_la_col_t* cycle_list;	// cycles derived from relations
    siqs_poly* curr_poly;		// current poly during filtering

    int is_restart;
    int is_tiny;
    int in_mem;
    int flag;

    //storage of relations found during in-mem sieving
    uint32_t buffered_rels;
    uint32_t buffered_rel_alloc;
    siqs_r* in_mem_relations;

    int do_periodic_tlp_filter;
    uint32_t est_raw_rels_needed;
    int do_batch;
    int batch_buffer_id;
    int batch_run_override;
#ifdef USE_BATCH_FACTOR
    // ping-pong relation batches so we can be processing
    // while loading data into one that's unused.
    relation_batch_t* rb;
    int num_alloc_rb;
    int num_active_rb;
    int max_active_rb;
#endif

#ifdef GATHER_RESIDUE_STATS
    FILE* residue_files[64];
#endif

} static_conf_t;

typedef struct {
    // the stuff in this structure is continuously being updated
    // during the course of the factorization, but we want it all
    // in one place so the data can be passed around easily
    int tid;

    //small prime sieving
    uint8_t* sieve;				// scratch space used for one sieve block 
    sieve_fb_compressed* comp_sieve_p;		// scratch space for a packed versions of fb
    sieve_fb_compressed* comp_sieve_n;		// for use during sieving smallish primes

    

    //large prime sieving
    update_t update_data;		// data for updating root values
    lp_bucket* buckets;			// bins holding sieve updates
    int* rootupdates;			// updates to apply to roots of primes
    uint16_t* sm_rootupdates;			// updates to apply to roots of primes
    uint32_t* mask2;

#ifdef USE_XLBUCKET
    xlp_bucket_t xl_pbucket;
    xlp_bucket_t xl_nbucket;
#endif

#if defined(USE_POLY_BUCKET_SS)
    ss_bucket_slice_t* ss_slices_p;
    ss_bucket_slice_t* ss_slices_n;
    int poly_buckets_allocated;
#elif defined( USE_DIRECT_SIEVE_SS )
    int ss_sieve_sz;
    uint16_t* report_ht_p;
    uint16_t* report_ht_n;
    int report_ht_size;
#endif

#ifdef USE_SS_SEARCH
    int* polymap;
    int* polynums;
    int* polyv;
    int* polysign;

    mpz_t polyb1;
    mpz_t polyb2;

    int using_ss_search;
    uint32_t num_ss_slices;

    int* firstroot1a;
    int* firstroot1b;
    int* firstroot2;
    uint8_t* ss_sieve_p;
    uint8_t* ss_sieve_n;
#endif

    //scratch
    mpz_t gmptmp1;
    mpz_t gmptmp2;
    mpz_t gmptmp3;

    // used in trial division
    uint16_t* mask;
    uint32_t* reports;			//sieve locations to submit to trial division
    uint32_t num_reports;
    mpz_t* Qvals;
    int* valid_Qs;				//which of the report are still worth persuing after SPV check
    uint32_t fb_offsets[MAX_SIEVE_REPORTS][MAX_SMOOTH_PRIMES];
    int* smooth_num;			//how many factors are there for each valid Q
    uint32_t failed_squfof;
    uint32_t attempted_squfof;
    uint32_t num_slp;
    uint32_t num_full;
    uint32_t dlp_outside_range;
    uint32_t dlp_prp;
    uint32_t dlp_useful;
    uint32_t failed_cosiqs;
    uint32_t attempted_cosiqs;
    uint32_t tlp_outside_range;
    uint32_t tlp_prp;
    uint32_t tlp_useful;
    uint32_t qlp_outside_range;
    uint32_t qlp_prp;
    uint32_t qlp_useful;
    uint64_t total_reports;
    uint64_t total_surviving_reports;
    uint64_t total_blocks;
    uint64_t num_scatter_opp;
    uint64_t num_scatter;
    uint32_t num_lp;

    // various things to support tlp co-factorization
    monty_t* mdata;		// monty brent attempt
    fact_obj_t* fobj2;	// smallmpqs attempt
    tiny_qs_params* cosiqs;

#ifdef USE_8X_MOD_ASM
    uint16_t* bl_sizes;
    uint16_t* bl_locs;
#endif

    //polynomial info during sieving
    siqs_poly* curr_poly;		// current poly during sieving	
    mpz_t* Bl;					// array of Bl values used to compute new B polys
    uint32_t tot_poly, numB, maxB;// polynomial counters
    int poly_batchsize;

    //storage of relations found during sieving
    uint32_t buffered_rels;
    uint32_t buffered_rel_alloc;
    siqs_r* relation_buf;

    uint64_t* unfactored_residue;
    uint64_t* residue_factors;
    uint32_t num_64bit_residue;

    uint32_t* polyscratch;
    uint16_t* corrections;

    //counters and timers
    uint32_t lp_scan_failures;
    uint64_t num;					// sieve locations we've subjected to trial division
    double rels_per_sec;

    char buf[512];
    uint32_t cutoff;

    uint32_t tf_small_cutoff;		// bit level to determine whether to bail early from tf
    int do_batch;
    int batch_run_override;
#ifdef USE_BATCH_FACTOR
    relation_batch_t rb;
#endif

    uint64_t lcg_state;

} dynamic_conf_t;

typedef struct {
    static_conf_t* sconf;
    dynamic_conf_t* dconf;
} thread_sievedata_t;

// sieving
void med_sieveblock_32k(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init);
void med_sieveblock_32k_sse41(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init);
void med_sieveblock_32k_avx2(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init);
void med_sieveblock_32k_avx512bw(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init);
extern void (*med_sieve_ptr)(uint8_t*, sieve_fb_compressed*, fb_list*, uint32_t, uint8_t);

void lp_sieveblock(uint8_t* sieve, uint32_t bnum, uint32_t numblocks,
    lp_bucket* lp, int side, dynamic_conf_t* dconf);
void lp_sieveblock_avx512f(uint8_t* sieve, uint32_t bnum, uint32_t numblocks,
    lp_bucket* lp, int side, dynamic_conf_t* dconf);
void lp_sieveblock_avx512bw(uint8_t* sieve, uint32_t bnum, uint32_t numblocks,
    lp_bucket* lp, int side, dynamic_conf_t* dconf);
extern void (*lp_sieveblock_ptr)(uint8_t* , uint32_t , uint32_t ,
    lp_bucket* , int , dynamic_conf_t* );

// sieving of primes using subset-sum
void lp_sieve_ss(uint8_t* sieve, int side, dynamic_conf_t* dconf);

// trial division
int check_relations_siqs_4_sse2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf);
int check_relations_siqs_8_sse2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf);
int check_relations_siqs_16_sse2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf);
int check_relations_siqs_4_avx2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf);
int check_relations_siqs_8_avx2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf);
int check_relations_siqs_16_avx2(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf);
int check_relations_siqs_16_avx512(uint32_t blocknum, uint8_t parity,
    static_conf_t* sconf, dynamic_conf_t* dconf);
extern int (*scan_ptr)(uint32_t, uint8_t, static_conf_t*, dynamic_conf_t*);

void filter_SPV(uint8_t parity, uint8_t* sieve, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);

void tdiv_LP_sse2(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void tdiv_LP_avx2(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void tdiv_LP_avx512(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
extern void (*tdiv_LP_ptr)(uint32_t, uint8_t, uint32_t,
    static_conf_t* , dynamic_conf_t* );

void tdiv_medprimes_32k(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void tdiv_medprimes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
// TBD
//void tdiv_medprimes_32k_knl(uint32_t* reports, uint32_t num_reports,
//    uint8_t parity, uint32_t poly_id, uint32_t bnum,
//    static_conf_t* sconf, dynamic_conf_t* dconf);
extern void (*tdiv_med_ptr)(uint8_t, uint32_t, uint32_t,
    static_conf_t*, dynamic_conf_t*);

void resieve_medprimes_32k(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void resieve_medprimes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void resieve_medprimes_32k_avx512bw(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
extern void (*resieve_med_ptr)(uint8_t, uint32_t, uint32_t,
    static_conf_t*, dynamic_conf_t*);

void trial_divide_Q_siqs(uint32_t report_num,
    uint8_t parity, uint32_t poly_id, uint32_t blocknum,
    static_conf_t* sconf, dynamic_conf_t* dconf);

void buffer_relation(uint32_t offset, uint32_t* large_prime, uint32_t num_factors,
    uint32_t* fb_offsets, uint32_t apoly_id, uint32_t poly_id, uint32_t parity,
    dynamic_conf_t* conf, uint32_t* polya_factors,
    uint32_t num_polya_factors, uint64_t unfactored_residue);

void save_relation_siqs(uint32_t offset, uint32_t* large_prime, uint32_t num_factors,
    uint32_t* fb_offsets, uint32_t poly_id, uint32_t apoly_id, uint32_t parity,
    static_conf_t* conf);

// trial division for primes sieved using subset-sum
void tdiv_SS(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);

// threading and misc
void siqs_sync(void* vptr);
void siqs_dispatch(void* vptr);
void siqs_work_fcn(void* vptr);
void siqs_thread_start(void* vptr);
void* process_poly(void* ptr);
int free_sieve(dynamic_conf_t* dconf);
int update_final(static_conf_t* sconf);
int update_check(static_conf_t* sconf);
int free_siqs(static_conf_t* sconf);
void free_filter_vars(static_conf_t* sconf);

//poly
void new_poly_a(static_conf_t* sconf, dynamic_conf_t* dconf);
void computeBl(static_conf_t* sconf, dynamic_conf_t* dconf, int needC);
void nextB(dynamic_conf_t* dconf, static_conf_t* sconf, int needC);

void firstRoots_32k(static_conf_t* sconf, dynamic_conf_t* dconf);
extern void (*firstRoots_ptr)(static_conf_t*, dynamic_conf_t*);

// functions implementing the subset-sum algorithm
void ss_search_setup(static_conf_t* sconf, dynamic_conf_t* dconf);
void ss_search_poly_buckets_sorted(static_conf_t* sconf, dynamic_conf_t* dconf);
void ss_search_poly_buckets_binned(static_conf_t* sconf, dynamic_conf_t* dconf);
void ss_search_clear(static_conf_t* sconf, dynamic_conf_t* dconf);
void ss_search_poly_buckets_2(static_conf_t* sconf, dynamic_conf_t* dconf, int set2_poly_id);
void ss_search_sort_set_1(static_conf_t* sconf, dynamic_conf_t* dconf);

void nextRoots_32k(static_conf_t* sconf, dynamic_conf_t* dconf);
void nextRoots_32k_sse41(static_conf_t* sconf, dynamic_conf_t* dconf);
void nextRoots_32k_avx2(static_conf_t* sconf, dynamic_conf_t* dconf);
void nextRoots_32k_avx2_small(static_conf_t* sconf, dynamic_conf_t* dconf);
void nextRoots_32k_avx2_intrin(static_conf_t* sconf, dynamic_conf_t* dconf);
extern void (*nextRoots_ptr)(static_conf_t*, dynamic_conf_t*);

void nextRoots_32k_generic_small(static_conf_t* sconf, dynamic_conf_t* dconf);
void nextRoots_32k_generic_polybatch(static_conf_t* sconf, dynamic_conf_t* dconf);

void nextRoots_32k_knl_bucket(static_conf_t* sconf, dynamic_conf_t* dconf);
void nextBigRoots_32k_knl_polybatch(static_conf_t* sconf, dynamic_conf_t* dconf);

void testfirstRoots_32k(static_conf_t* sconf, dynamic_conf_t* dconf);
extern void (*testRoots_ptr)(static_conf_t*, dynamic_conf_t*);

void batch_roots(int* rootupdates, int* firstroots1, int* firstroots2,
    siqs_poly* poly, uint32_t start_prime, fb_list* fb, uint32_t* primes);

//data I/O
uint32_t process_poly_a(static_conf_t* sconf);
int get_a_offsets(fb_list* fb, siqs_poly* poly, mpz_t tmp);
void generate_bpolys(static_conf_t* sconf, dynamic_conf_t* dconf, int maxB);
int process_rel(char* substr, fb_list* fb, mpz_t n,
    static_conf_t* sconf, fact_obj_t* obj, siqs_r* rel);
int restart_siqs(static_conf_t* sconf, dynamic_conf_t* dconf);
uint32_t qs_purge_singletons(fact_obj_t* obj, siqs_r* list,
    uint32_t num_relations,
    qs_cycle_t* table, uint32_t* hashtable);
uint32_t qs_purge_singletons3(fact_obj_t* obj, siqs_r* list,
    uint32_t num_relations,
    qs_cycle_t* table, uint32_t* hashtable);
uint32_t qs_purge_singletonsN(fact_obj_t* obj, siqs_r* list,
    uint32_t num_relations, uint32_t num_lp,
    qs_cycle_t* table, uint32_t* hashtable);
uint32_t qs_purge_duplicate_relations(fact_obj_t* obj,
    siqs_r* rlist,
    uint32_t num_relations);
uint32_t qs_purge_duplicate_relations3(fact_obj_t* obj,
    siqs_r* rlist,
    uint32_t num_relations);
void qs_enumerate_cycle(fact_obj_t* obj,
    qs_la_col_t* c,
    qs_cycle_t* table,
    qs_cycle_t* entry1, qs_cycle_t* entry2,
    uint32_t final_relation);
int yafu_sort_cycles(const void* x, const void* y);

//aux
uint8_t choose_multiplier_siqs(uint32_t B, mpz_t n);
int siqs_static_init(static_conf_t* sconf, int is_tiny);
int siqs_dynamic_init(dynamic_conf_t* dconf, static_conf_t* sconf);
int siqs_check_restart(dynamic_conf_t* dconf, static_conf_t* sconf);
uint32_t siqs_merge_data(dynamic_conf_t* dconf, static_conf_t* sconf);

void get_params(static_conf_t* sconf);
void get_gray_code(siqs_poly* poly);
void set_aprime_roots(static_conf_t* sconf, uint32_t val, int* qli, int s,
    sieve_fb_compressed* fb, int action);
void siqsexit(int sig);
int qcomp_siqs(const void* x, const void* y);
uint32_t make_fb_siqs(static_conf_t* sconf);
void get_dummy_params(int bits, uint32_t* B, uint32_t* M, uint32_t* NB);
void print_siqs_splash(dynamic_conf_t* dconf, static_conf_t* sconf);

// tiny variants of a few routines, that live in tinySIQS.c
int tiny_update_check(static_conf_t* sconf);
void* tiny_process_poly(void* ptr);

//test routines
int check_specialcase(FILE* sieve_log, fact_obj_t* fobj);
void test_polya(fb_list* fb, mpz_t n, mpz_t target_a);
int checkpoly_siqs(siqs_poly* poly, mpz_t n);
int checkBl(mpz_t n, uint32_t* qli, fb_list* fb, mpz_t* Bl, int s);
int check_relation(mpz_t a, mpz_t b, siqs_r* r, fb_list* fb, mpz_t n, int VFLAG, int knmod8);
int check_poly_at_loc(int loc, int parity, int prime, int polynum, int index,
    int* polyv, int* polysign, dynamic_conf_t* dconf, static_conf_t* sconf);


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

int qs_solve_linear_system(fact_obj_t* obj, uint32_t fb_size,
    uint64_t** bitfield,
    siqs_r* relation_list,
    qs_la_col_t* cycle_list,
    uint32_t* num_cycles);

/* merge src1[] and src2[] into merge_array[], assumed
   large enough to hold the merged result. Return the
   final number of elements in merge_array */

uint32_t qs_merge_relations(uint32_t* merge_array,
    uint32_t* src1, uint32_t n1,
    uint32_t* src2, uint32_t n2);

uint64_t* qs_block_lanczos(fact_obj_t* obj,
    uint32_t nrows,
    uint32_t num_dense_rows,
    uint32_t ncols,
    qs_la_col_t* cols,
    uint32_t* deps_found);

void count_qs_matrix_nonzero(fact_obj_t* obj,
    uint32_t nrows, uint32_t num_dense_rows,
    uint32_t ncols, qs_la_col_t* cols);

void reduce_qs_matrix(fact_obj_t* obj, uint32_t* nrows,
    uint32_t num_dense_rows, uint32_t* ncols,
    qs_la_col_t* cols, uint32_t num_excess);

void qs_free_cycle_list(qs_la_col_t* cycle_list, uint32_t num_cycles);

/*-------------- MPQS SQUARE ROOT RELATED DECLARATIONS ---------------------*/

uint32_t yafu_find_factors(fact_obj_t* obj, mpz_t n, fb_element_siqs* factor_base,
    uint32_t fb_size, qs_la_col_t* vectors,
    uint32_t vsize, siqs_r* relation_list,
    uint64_t* null_vectors, uint32_t multiplier,
    mpz_t* a_list, poly_t* poly_list,
    factor_list_t* factor_list);

uint32_t yafu_factor_list_add(fact_obj_t* obj,
    factor_list_t* list,
    mpz_t new_factor);

/*-------------- CYCLE FINDING RELATED DECLARATIONS ------------------------*/
/* pull out the large primes from a relation read from
   the savefile */
void rebuild_graph(static_conf_t* sconf, siqs_r* relation_list, int num_relations);
void rebuild_graph3(static_conf_t* sconf, siqs_r* relation_list, int num_relations);
void yafu_read_large_primes(char* buf, uint32_t* prime1, uint32_t* prime2);
void yafu_read_tlp(char* buf, uint32_t* primes);
int yafu_read_Nlp(char* buf, uint32_t* primes);

/* given the primes from a sieve relation, add
   that relation to the graph used for tracking
   cycles */

void yafu_add_to_cyclesN(static_conf_t* conf, uint32_t flags, uint32_t* primes);
void yafu_add_to_cycles3(static_conf_t* conf, uint32_t flags, uint32_t* primes);
void yafu_add_to_cycles(static_conf_t* conf, uint32_t flags, uint32_t prime1, uint32_t prime2);

qs_la_col_t* find_cycles(fact_obj_t* obj, uint32_t* hashtable, qs_cycle_t* table,
    siqs_r* relation_list, uint32_t num_relations, uint32_t* numcycles, uint32_t* numpasses);

qs_la_col_t* find_cycles3(fact_obj_t* obj, static_conf_t* sconf,
    siqs_r* relation_list, uint32_t num_relations, uint32_t* numcycles, uint32_t* numpasses);

/* perform postprocessing on a list of relations */
void yafu_qs_filter_relations(static_conf_t* sconf);

extern uint64_t* siqs_primes;
extern uint64_t siqs_nump;
extern uint64_t siqs_minp;
extern uint64_t siqs_maxp;




#endif
