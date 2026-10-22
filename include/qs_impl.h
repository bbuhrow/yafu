#ifndef _SIQS_IMPL_H_
#define _SIQS_IMPL_H_

// =====================================================================
// SIQS small-prime-variation (SPV) tuning experiment.
// See docs/siqs_spv_plan.md.  Commenting out this single define reverts
// the qs library to pre-experiment behavior with zero code deletion:
// every experimental edit in the qs library is wrapped as
//   #ifdef TRY_SMALLP_OPT  <new code>
//   #else                  <original code, preserved verbatim>
//   #endif
#define TRY_SMALLP_OPT 1


#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <time.h>
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

#ifdef HAVE_CUDA_BATCH_FACTOR
#include "gpu_cofactorization.h"
#endif

#ifdef HAVE_OCL_BATCH_FACTOR
#include "gpu_cofactorization_cl.h"
#endif

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

// harebrained idea, didn't work :(
#define NUM_ALP 0

// new largep bucket sorting idea
// #define TRY_COMPRESS_SORT_LARGEP 1

// The experimental sieves below (L2 / sparse / list / direct) have only
// AVX-512F/BW implementations so far.  Lower-ISA builds compile the stock
// bucket sieve instead.
#if defined(USE_AVX512F) && defined(USE_AVX512BW)
#define SIQS_EXPERIMENTAL_SIEVES
#endif

#if defined(SIQS_EXPERIMENTAL_SIEVES)
// L2 sieving
#define TRY_L2_SIEVE
// change the L2 sieve to a sparse sieve up to 8x blocksize primes
#define TRY_SPARSE_SIEVE

// experimental linked-list sieve for primes in [med_B, B): one list per
// sieve block; each node is one (prime, root) hit.  Sieving a block walks
// its list, and each visited hit advances its root by prime, re-appending
// to a later block's list if the new root is still inside the interval.
// Requires:
//  1. per-poly root updates extended to cover [med_B, B) (guarded edits in
//     update_poly_roots_32k_avx2.c and poly_roots_32k.c)
//  2. dconf->sieve_p/n arrays allocated with B entries instead of fb_blk8i
//     (guarded edit in SIQS.c)
//  3. one list_sieve_ctx_t per thread per side (list_sieve_alloc_32k),
//     wired into process_poly (sieve) and tdiv_scan.c (resieve)
#define TRY_LIST_SIEVE

// experimental direct interval sieve for the primes p > interval,
// fb indices [large_B, B) (see direct_sieve.md).  Both sides' whole
// intervals live in memory at once (pos then neg, 2*num_blocks blocks).
// Per poly:
//  1. the intervals arrive already holding [large_B, B)'s hits for this
//     poly (from the previous poly's fused pass, or for the first poly of
//     an 'a' from a standalone pass 1 straight from update_data.firstroots).
//  2. the block loops run med + sparse sieves and the scan per block as
//     before (scan sees complete sieve bytes); [med_B, large_B) keeps the
//     TRY_SPARSE_SIEVE machinery.  Survivors of the scan are moved into a
//     per-poly pool instead of being finalized.
//  3. the intervals are re-initialized and one fused pass over [large_B, B)
//     resieves the survivors (bitmap of their interval positions, guarded
//     divisions), updates the firstroots to the next poly, and sieves the
//     next poly's roots.
//  4. survivors are finalized (trial_divide_Q_siqs), then nextB.
// The bucket sieve is not used.  Requires TRY_L2_SIEVE + TRY_SPARSE_SIEVE,
// AVX512F, and check_relations_siqs_16_avx512 as the scan.  Not
// compatible with USE_SS_SEARCH / USE_BATCHPOLY / USE_XLBATCHPOLY.
// Takes precedence over TRY_LIST_SIEVE (which is undefined below).
#define TRY_DIRECT_SIEVE
#endif  // SIQS_EXPERIMENTAL_SIEVES

#if defined(TRY_DIRECT_SIEVE)
#undef TRY_LIST_SIEVE
#if !defined(TRY_L2_SIEVE) || !defined(TRY_SPARSE_SIEVE)
#error "TRY_DIRECT_SIEVE requires TRY_L2_SIEVE and TRY_SPARSE_SIEVE"
#endif
#if defined(USE_SS_SEARCH) || defined(USE_BATCHPOLY) || defined(USE_XLBATCHPOLY)
#error "TRY_DIRECT_SIEVE is not compatible with SS search or batch polys"
#endif
#if !defined(USE_AVX512F) || !defined(USE_AVX512BW)
#error "TRY_DIRECT_SIEVE requires USE_AVX512F and USE_AVX512BW"
#endif
// added to the logp the direct sieve uses for [large_B, B).  The bucket
// sieve uses one logp per slice, the logp of the slice's last prime, so it
// credits these primes slightly more than the per-16 approximation; a bias
// of 1 roughly matches it.
#ifndef DIRECT_LOGP_BIAS
#define DIRECT_LOGP_BIAS 0
#endif
// the large-prime resieve finds survivor hits by comparing roots against
// the survivor positions directly when a poly has at most this many
// survivors (both sides), and by gathering from the survivor bitmap
// otherwise.  0 = always gather.  Without the prefilter, compares won at
// ~10 survivors/poly and lost above ~17 (Rocket Lake, 77 digits); with the
// 1024-bit prefilter they gave no measurable gain over gathers, so the
// default is gathers only.  Kept selectable for other input sizes.
#ifndef DIRECT_CMP_MAX_SURVIVORS
#define DIRECT_CMP_MAX_SURVIVORS 0
#endif
// register-resident prefilter in front of the survivor compares/gathers:
// one bit per bucket of positions across both intervals says whether the
// bucket holds a survivor.  0 (off), 512 (one zmm) or 1024 (two zmm).
// 1024 measured best: -12% resieve time at NB 2/4, -9% at NB=8.
#ifndef DIRECT_PREFILTER_BITS
#define DIRECT_PREFILTER_BITS 1024
#endif
// 1: handle [med_B, large_B) with the stock bucket sieve
// (nextRoots_32k_knl_medbucket + lp_sieveblock + tdiv_LP) instead of the
// TRY_SPARSE_SIEVE per-block sparse sieve (large_sieveblock + tdiv_large).
// Implemented by setting fb_blk8i = med_B, which empties every sparse-range
// loop, while large_B stays the first 16-aligned index with prime > interval.
#ifndef DIRECT_MEDBUCKET
#define DIRECT_MEDBUCKET 0
#endif
// how the [large_B, B) roots are scattered into the intervals (standalone
// pass 1 and the fused pass's next-poly sieve), see Sink-Stripe Sieve
// Stores.md:
// 0: loop over the set bits of the in-range masks (trail_zcnt / reset_lsb).
// 1: "sink stripe": every one of the 64 root lanes per 16 primes does a
//    read-modify-write; lanes that miss the interval are blended to a
//    rotating address in an unscanned stripe past the neg interval (the
//    second tail guard block), so there are no data-dependent branches.
// 2: compress the in-range lanes (as logp << 24 | position) in registers
//    (maskz_compress) and append them to a buffer with unaligned stores;
//    sieve the buffer once it holds DIRECT_CBUF_FLUSH entries.
// 3: as 2, but append with mask_compressstoreu directly to memory.
#ifndef DIRECT_SCATTER_MODE
#define DIRECT_SCATTER_MODE 3
#endif
// size of the sink stripe in bytes: a multiple of 64, at most 32768
#ifndef DIRECT_SINK_BYTES
#define DIRECT_SINK_BYTES 1024
#endif
// modes 2/3: sieve the compressed hit buffer once it holds at least this
// many entries (1 = after every 16 primes).  Flat from 256 to 4096 on
// Rocket Lake (mode 3, 77-79 digits, NB 2/4/8); 128 is ~1-2% slower at NB=2.
#ifndef DIRECT_CBUF_FLUSH
#define DIRECT_CBUF_FLUSH 512
#endif
// 1: the sparse sieve for [med_B, large_B) (large_sieveblock_32k_avx512bw)
// appends its block hits to a buffer with compress stores (mode 2 style
// when DIRECT_SCATTER_MODE == 2, compressstoreu otherwise) and sieves the
// buffer in one loop, and tdiv_large_32k_avx2 collects each report's
// matching fb indices the same way.  0: trail_zcnt / reset_lsb loops.
// Only used when DIRECT_MEDBUCKET is 0.
#ifndef DIRECT_SPARSE_COMPRESS
#define DIRECT_SPARSE_COMPRESS 1
#endif
#endif

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
#endif

// assume we have SSE2 available on x86 architectures
#if defined(__x86_64__) && !defined(FORCE_GENERIC)
#define D_HAS_SSE2
#endif

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

typedef struct
{
    uint32_t* prime;
    uint32_t* root1;
    uint32_t* root2;
    uint32_t* logp;
} sieve_fb;

// context for the experimental linked-list sieve (TRY_LIST_SIEVE).
// One instance per thread per side (p/n).  Node ids encode a factor base
// index and which root of that prime the node represents:
//      node = 2 * (fb_index - start) + root_slot      (root_slot 0 or 1)
// so next[] has 2*(stop-start) entries.  The sieve and the resieve run the
// same procedure in lockstep, block by block, so each keeps its OWN set of
// per-block lists (both hold pending hits for future blocks at the same
// time and cannot share storage).  The sieve's lists are built from the
// current poly's fb->root1/2; the resieve's lists are built from the
// pristine roots saved in oroot1/oroot2 and advanced in rroot1/rroot2.
#define LISTSIEVE_NULL 0xffffffffu

// append node to the end of block b's list
#define LIST_SIEVE_APPEND(hd, tl, nx, b, node) \
    do { \
        if ((tl)[(b)] == LISTSIEVE_NULL) \
            (hd)[(b)] = (node); \
        else \
            (nx)[(tl)[(b)]] = (node); \
        (nx)[(node)] = LISTSIEVE_NULL; \
        (tl)[(b)] = (node); \
    } while (0)

typedef struct
{
    uint32_t start;         // first fb index covered by the list sieve
    uint32_t stop;          // one past the last fb index covered
    uint32_t numblocks;     // blocks per poly interval
    uint32_t* next;         // [2*(stop-start)] sieve list links
    uint32_t* head;         // [numblocks] sieve per-block list heads
    uint32_t* tail;         // [numblocks] sieve per-block list tails
    uint32_t* rnext;        // [2*(stop-start)] resieve list links
    uint32_t* rhead;        // [numblocks] resieve per-block list heads
    uint32_t* rtail;        // [numblocks] resieve per-block list tails
    uint32_t* oroot1;       // [stop-start] pristine roots saved at init
    uint32_t* oroot2;       // [stop-start]
    uint32_t* rroot1;       // [stop-start] resieve running roots
    uint32_t* rroot2;       // [stop-start]
    uint16_t* flags;        // [32768] resieve: report_num+1 marks, 0 = none
    int resieve_ready;      // resieve lists built for the current poly
    uint32_t rblock;        // next block index the resieve state expects
} list_sieve_ctx_t;

// per-poly pool of scan survivors for the direct interval sieve
// (TRY_DIRECT_SIEVE).  Survivors are appended in scan order (pos side
// blocks 0..nb-1, then neg side), so pos[] is sorted ascending.
//      pos = side * interval + block * 32768 + block_loc
typedef struct
{
    uint32_t num;           // survivors in the pool for the current poly
    uint32_t alloc;         // allocated entries
    uint32_t interval;      // num_blocks * 32768
    uint32_t* pos;          // [alloc] interval position, both sides
    mpz_t* Q;               // [alloc] partially factored Q
    int* smooth_num;        // [alloc] index of last stored fb offset
    uint32_t* fb_offsets;   // [alloc * MAX_SMOOTH_PRIMES]
    uint64_t* bitmap;       // [2 * interval / 64] marks survivor positions
} direct_pool_t;

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
#ifdef TRY_SMALLP_OPT
    // Lemire-style divisibility constants for the vectorized SPV check
    // in tdiv_small.c: invmod = p^-1 mod 2^32, thresh = (2^32-1)/p
    uint32_t* invmod;
    uint32_t* thresh;
#endif
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
    uint32_t* inv32;
    uint32_t* corr32;
    uint32_t* shift32;
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
#ifdef TRY_COMPRESS_SORT_LARGEP
    uint32_t lp_alloc_slices;	    // the number of large prime fb slices allocated
    uint32_t lp_num_slices;		    // the number of large prime fb slices currently used
    uint32_t* lp_id_p;			    // array of current id in each slice
    uint32_t* lp_id_n;			    // array of current id in each slice
    uint32_t* lp_fb_bounds;		    // array of fb offsets for each slice
    uint8_t* lp_logp;			    // array of logp values for each slice
    uint32_t* lp_list;			    // contiguous space for all lp entries
#endif
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
    uint32_t fb_blk8i;          // index at which primes are bigger than 8x a block (and a multiple of 8)
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

    uint64_t large_prime_max;		// the cutoff value for keeping a partial
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

    // master list of cycles
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

    // used to check on progress of the factorization	
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

    // master time record
    double t_time1;				// sieve time
    double t_time2;				// relation scanning and trial division
    double t_time3;				// polynomial calculations and large prime sieving
    double t_time4;				// extra?

    double t_mpsieve;
    double t_lpsieve;
    double t_tdiv;
    double t_polyupdate;

    // stuff for sqrt
    factor_list_t factor_list;

    // these are used during linear algebra and sqrt root
    uint32_t total_poly_a;		// used during sieving: total number of polynomial 'a' values 
    mpz_t* poly_a_list;			// used during sieving: list of 'a' values for MPQS polys 
    uint32_t filt_total_poly_a;	// used during filtering: total number of polynomial 'a' values 
    mpz_t* filt_poly_a_list;	// used during filtering: list of 'a' values for MPQS polys 
    poly_t* poly_list;			// list of MPQS polynomials 
    uint32_t poly_list_alloc;
    uint32_t apoly_alloc;

    // these are used during filtering
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

    // storage of relations found during in-mem sieving
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

#if defined(HAVE_CUDA_BATCH_FACTOR) || defined(HAVE_OCL_BATCH_FACTOR)

    // a gpu device context
    int which_gpu;
    device_ctx_t* gpu_dev_ctx;

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

    sieve_fb* sieve_p;		// scratch space for a packed versions of fb
    sieve_fb* sieve_n;		// for use during sieving medium primes
#if defined(TRY_LIST_SIEVE)
    list_sieve_ctx_t* lsctx_p;	// linked-list sieve context, positive side
    list_sieve_ctx_t* lsctx_n;	// linked-list sieve context, negative side
#endif
#if defined(TRY_DIRECT_SIEVE)
    direct_pool_t* dpool;		// scan survivors awaiting the large-prime resieve
#endif

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
    uint32_t num_alp;

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

    double t_mpsieve;
    double t_lpsieve;
    double t_tdiv;
    double t_polyupdate;

    uint64_t lcg_state;

#if defined(HAVE_CUDA_BATCH_FACTOR) || defined(HAVE_OCL_BATCH_FACTOR)
    device_thread_ctx_t* gpu_cofactor_ctx;
#endif

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
void L2_sieveblock_32k_avx512bw(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init);
void large_sieveblock_32k_avx512bw(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init);
list_sieve_ctx_t* list_sieve_alloc_32k(uint32_t start, uint32_t stop, uint32_t numblocks);
void list_sieve_free_32k(list_sieve_ctx_t* ctx);
void list_sieve_init_32k(sieve_fb* fb, fb_list* full_fb, list_sieve_ctx_t* ctx);
void list_sieveblock_32k_avx512bw(uint8_t* blk_sieve, uint32_t bnum,
    sieve_fb* fb, fb_list* full_fb, list_sieve_ctx_t* ctx);
void direct_sieve_largeprimes_32k(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf);

#if defined(TRY_DIRECT_SIEVE) && ((DIRECT_SCATTER_MODE >= 2) || DIRECT_SPARSE_COMPRESS)
#include <immintrin.h>
// buffer size for modes 2/3: a flush threshold plus room for one more
// 16-prime group (64 entries) and the last unaligned store's 16-lane tail
#define DIRECT_CBUF_SIZE (DIRECT_CBUF_FLUSH + 80)

// append the lanes of v selected by m to buf[n...]; returns the new count.
// buf must have 16 entries of room past n (the register-compress form
// stores a full vector).
static inline uint32_t direct_cstore16(uint32_t* buf, uint32_t n, __mmask16 m, __m512i v)
{
#if DIRECT_SCATTER_MODE == 2
    _mm512_storeu_si512((__m512i*)(buf + n), _mm512_maskz_compress_epi32(m, v));
#else
    _mm512_mask_compressstoreu_epi32(buf + n, m, v);
#endif
    return n + _mm_popcnt_u32(m);
}

// append the in-range lanes of 16 primes' 4 root vectors to buf[n...] as
// (logp << 24) | sieve position (pos side r, neg side interval + p - r;
// 2*interval < 2^24).  Returns the new count.
static inline uint32_t direct_compress16(uint32_t* buf, uint32_t n, uint32_t interval,
    __m512i vroot1, __m512i vroot2, __m512i vnroot1, __m512i vnroot2, uint8_t logp)
{
    __m512i vint = _mm512_set1_epi32(interval);
    __m512i vtag = _mm512_set1_epi32((uint32_t)logp << 24);
    __mmask16 m1 = _mm512_cmplt_epu32_mask(vroot1, vint);
    __mmask16 m2 = _mm512_cmplt_epu32_mask(vroot2, vint);
    __mmask16 m3 = _mm512_cmplt_epu32_mask(vnroot1, vint);
    __mmask16 m4 = _mm512_cmplt_epu32_mask(vnroot2, vint);
    __m512i v1 = _mm512_or_epi32(vroot1, vtag);
    __m512i v2 = _mm512_or_epi32(vroot2, vtag);
    __m512i v3 = _mm512_or_epi32(_mm512_add_epi32(vnroot1, vint), vtag);
    __m512i v4 = _mm512_or_epi32(_mm512_add_epi32(vnroot2, vint), vtag);

    n = direct_cstore16(buf, n, m1, v1);
    n = direct_cstore16(buf, n, m2, v2);
    n = direct_cstore16(buf, n, m3, v3);
    return direct_cstore16(buf, n, m4, v4);
}

// sieve the buffered hits
static inline void direct_compress_flush(uint8_t* sieve, const uint32_t* buf, uint32_t n)
{
    uint32_t k;
    for (k = 0; k < n; k++)
    {
        uint32_t e = buf[k];
        sieve[e & 0xffffff] -= (uint8_t)(e >> 24);
    }
}
#endif

#if defined(TRY_DIRECT_SIEVE) && (DIRECT_SCATTER_MODE == 1)
#include <immintrin.h>
// sink-stripe scatter of 16 primes' 4 root vectors into both intervals
// (pos side sieve[0, interval), neg side sieve[interval, 2*interval)).
// vroot = pos-side roots, vnroot = p - root (neg-side offsets).  Every lane
// is written: lanes outside the interval go to the stripe at
// sieve[2*interval + 32768 + *sink ...], which is allocated (second tail
// guard block) and never scanned.  *sink rotates through the stripe so
// consecutive miss RMWs don't chain through store forwarding.
static inline void direct_sink_scatter16(uint8_t* sieve, uint32_t interval,
    __m512i vroot1, __m512i vroot2, __m512i vnroot1, __m512i vnroot2,
    uint8_t logp, uint32_t* sink)
{
    ALIGNED_MEM uint32_t idx[64];
    const __m512i vlane = _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15);
    __m512i vint = _mm512_set1_epi32(interval);
    __m512i vs = _mm512_add_epi32(_mm512_set1_epi32(2 * interval + 32768 + *sink), vlane);
    __m512i v16 = _mm512_set1_epi32(16);
    __mmask16 m1 = _mm512_cmplt_epu32_mask(vroot1, vint);
    __mmask16 m2 = _mm512_cmplt_epu32_mask(vroot2, vint);
    __mmask16 m3 = _mm512_cmplt_epu32_mask(vnroot1, vint);
    __mmask16 m4 = _mm512_cmplt_epu32_mask(vnroot2, vint);

    _mm512_store_si512(idx, _mm512_mask_blend_epi32(m1, vs, vroot1));
    vs = _mm512_add_epi32(vs, v16);
    _mm512_store_si512(idx + 16, _mm512_mask_blend_epi32(m2, vs, vroot2));
    vs = _mm512_add_epi32(vs, v16);
    _mm512_store_si512(idx + 32, _mm512_mask_blend_epi32(m3, vs,
        _mm512_add_epi32(vnroot1, vint)));
    vs = _mm512_add_epi32(vs, v16);
    _mm512_store_si512(idx + 48, _mm512_mask_blend_epi32(m4, vs,
        _mm512_add_epi32(vnroot2, vint)));

#define DIRECT_SINK_RMW4(k) \
    sieve[idx[k]] -= logp; sieve[idx[k + 1]] -= logp; \
    sieve[idx[k + 2]] -= logp; sieve[idx[k + 3]] -= logp;
#define DIRECT_SINK_RMW16(k) \
    DIRECT_SINK_RMW4(k) DIRECT_SINK_RMW4(k + 4) \
    DIRECT_SINK_RMW4(k + 8) DIRECT_SINK_RMW4(k + 12)

    DIRECT_SINK_RMW16(0);
    DIRECT_SINK_RMW16(16);
    DIRECT_SINK_RMW16(32);
    DIRECT_SINK_RMW16(48);

#undef DIRECT_SINK_RMW16
#undef DIRECT_SINK_RMW4

    *sink += 64;
    if (*sink >= DIRECT_SINK_BYTES)
        *sink = 0;
}
#endif

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
void tdiv_LP_avx512_allreports(uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);

extern void (*tdiv_LP_ptr)(uint32_t, uint8_t, uint32_t,
    static_conf_t* , dynamic_conf_t* );

void tdiv_medprimes_32k(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void tdiv_medprimes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);

void tdiv_L2primes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void tdiv_large_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void resieve_large_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf, list_sieve_ctx_t* ctx);
direct_pool_t* direct_pool_alloc(uint32_t num_blocks);
void direct_pool_free(direct_pool_t* pool);
void direct_save_survivor(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf);
void direct_resieve_update_largeprimes_32k(static_conf_t* sconf, dynamic_conf_t* dconf,
    uint8_t* sieve);
void direct_finalize_survivors(static_conf_t* sconf, dynamic_conf_t* dconf);

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

void buffer_relation(uint32_t offset, uint64_t* large_prime, uint32_t num_factors,
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
void new_poly_a(static_conf_t* sconf, dynamic_conf_t* dconf, int gen_and_test_only);
void new_poly_a_test(static_conf_t* sconf, dynamic_conf_t* dconf, int gen_and_test_only);
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

void nextRoots_32k_knl_medbucket(static_conf_t* sconf, dynamic_conf_t* dconf,
    int* bi, int* cb, uint32_t* bv, int** uptr,
    uint32_t** np_p, uint32_t** np_n, uint32_t** sl_p, uint32_t** sl_n);
void nextRoots_32k_knl_medbucket_manyblocks(static_conf_t* sconf, dynamic_conf_t* dconf,
    int* bi, int* cb, uint32_t* bv, int** uptr,
    uint32_t** np_p, uint32_t** np_n, uint32_t** sl_p, uint32_t** sl_n);
void nextRoots_32k_knl_largebucket(static_conf_t* sconf, dynamic_conf_t* dconf,
    int bound_index, int check_bound, uint32_t bound_val, int* ptr,
    uint32_t* numptr_p, uint32_t* numptr_n, uint32_t* sliceptr_p, uint32_t* sliceptr_n);
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
void set_aprime_roots(static_conf_t* sconf, dynamic_conf_t* dconf, uint32_t val, 
    int* qli, int s, sieve_fb_compressed* fb, int action);
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
