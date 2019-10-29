/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 7/28/10
----------------------------------------------------------------------*/

#ifndef YAFU_SOE_H
#define YAFU_SOE_H

#include "yafu.h"
#include "arith.h"	//needed for spGCD and modinv_1
#include "util.h"
#include "immintrin.h"

#define USE_SOE_THREADPOOL


#ifdef _MSC_VER
// optionally define this or not depending on whether your hardware supports it.
// if defined, compile the sse41 functions into the fat binary.  the global
// flag HAS_SSE41 is set at runtime on compatible hardware to enable the functions
// to be used.  For gcc and mingw64 builds, USE_SSE41 is enabled in the makefile.
//#define USE_SSE41 1
//#define USE_AVX2 1
#endif

//#define BLOCKSIZE 32768
//#define FLAGSIZE 262144
//#define FLAGSIZEm1 262143
//#define FLAGBITS 18
//#define BUCKETSTARTI 33336

//#define BLOCKSIZE 1048576
//#define FLAGSIZE 8388608
//#define FLAGSIZEm1 8388607
//#define FLAGBITS 23
//#define BUCKETSTARTI 846248

//#define BLOCKSIZE 524288
//#define FLAGSIZE 4194304
//#define FLAGSIZEm1 4194303
//#define FLAGBITS 22
////#define BUCKETSTARTI 191308
//#define BUCKETSTARTI 443920
//#define BUCKETSTARTP 393216 

#define BITSINBYTE 8
#define MAXSIEVEPRIMECOUNT 100000000	//# primes less than ~2e9: limit of 2e9^2 = 4e18
//#define DO_SPECIAL_COUNT

enum soe_command {
	SOE_COMMAND_INIT,
	SOE_COMMAND_WAIT,
	SOE_COMMAND_SIEVE_AND_COUNT,
	SOE_COMMAND_SIEVE_AND_COMPUTE,
	SOE_COMPUTE_ROOTS,
	SOE_COMPUTE_PRIMES,
	SOE_COMPUTE_PRPS,
	SOE_COMMAND_END
};

#ifdef USE_AVX2

#ifdef USE_AVX512F
extern ALIGNED_MEM uint64 presieve_largemasks[16][173][8];
extern ALIGNED_MEM uint32 presieve_steps[32];
extern ALIGNED_MEM uint32 presieve_primes[32];
extern ALIGNED_MEM uint32 presieve_p1[32];

#else
// for storage of presieving lists from prime index 24 to 40 (97 to 173 inclusive)
extern ALIGNED_MEM uint64 presieve_largemasks[16][173][4];
extern ALIGNED_MEM uint32 presieve_steps[32];
extern ALIGNED_MEM uint32 presieve_primes[32];
extern ALIGNED_MEM uint32 presieve_p1[32];
#endif

// macros for Montgomery arithmetic - helpful for computing 
// division-less offsets once we have enough reuse (number of
// classes) to justify the setup costs.
#define _mm256_cmpge_epu32(a, b) \
        _mm256_cmpeq_epi32(_mm256_max_epu32(a, b), a)

#define _mm256_cmple_epu32(a, b) _mm256_cmpge_epu32(b, a)


static __inline __m256i CLEAR_HIGH_VEC(__m256i x)
{
    __m256i chi = _mm256_set1_epi64x(0x00000000ffffffff);
    return _mm256_and_si256(chi, x);
}

static __inline __m256i vec_redc(__m256i x64e, __m256i x64o, __m256i pinv, __m256i p)
{
    // uint32 m = (uint32)x * pinv;
    __m256i t1 = _mm256_shuffle_epi32(pinv, 0xB1);      // odd-index pinv in lo words
    __m256i even = _mm256_mul_epu32(x64e, pinv);
    __m256i odd = _mm256_mul_epu32(x64o, t1);
    __m256i t2;

    // x += (uint64)m * (uint64)p;
    t1 = _mm256_shuffle_epi32(p, 0xB1);      // odd-index p in lo words
    even = _mm256_add_epi64(x64e, _mm256_mul_epu32(even, p));
    odd = _mm256_add_epi64(x64o, _mm256_mul_epu32(odd, t1));

    // m = x >> 32;
    t1 = _mm256_blend_epi32(odd, _mm256_shuffle_epi32(even, 0xB1), 0x55);

    // if (m >= p) m -= p;
    t2 = _mm256_cmpge_epu32(t1, p); //_mm256_or_si256(_mm256_cmpgt_epi32(t1, p), _mm256_cmpeq_epi32(t1, p));
    t2 = _mm256_and_si256(p, t2);

    return _mm256_sub_epi32(t1, t2);
}

static __inline __m256i vec_to_monty(__m256i x, __m256i r2, __m256i pinv, __m256i p)
{
    //uint64 t = (uint64)x * (uint64)r2;
    __m256i t1 = _mm256_shuffle_epi32(x, 0xB1);
    __m256i t2 = _mm256_shuffle_epi32(r2, 0xB1);
    __m256i even = _mm256_mul_epu32(x, r2);
    __m256i odd = _mm256_mul_epu32(t1, t2);

    //return redc_loc(t, pinv, p);
    return vec_redc(even, odd, pinv, p);
}

#endif

typedef struct
{
	uint16 loc;
	uint8 mask;
	uint8 bnum;
} soe_bucket_t_old;

//typedef struct
//{
//	uint32 root;
//	uint32 prime;
//} soe_bucket_t;

typedef struct
{
	//uint32 prime;		// the prime, so that we don't have to also look in the
						// main prime array
	uint32 bitloc;		// bit location of the current hit
	uint32 next_pid;	// index of the next prime that hits in the current sieve
	uint32 p_div;		// prime / prodN
	uint8 p_mod;		// prime % prodN
	uint8 eacc;			// accumulated error
} soe_bitmap_p;

typedef struct
{
    int sync_count;
	uint32 *sieve_p;
	int *root;
	uint32 *lower_mod_prime;

    uint32 *pinv;       // montgomery inverse
    uint32 *r2modp;     // to go out of montgomery rep

	uint64 blk_r;
	uint64 blocks;
	uint64 partial_block_b;
	uint64 prodN;
	uint64 startprime;
	uint64 orig_hlimit;
	uint64 orig_llimit;
	uint64 pbound;
	uint64 pboundi;

	uint32 bucket_start_id;
	uint32 large_bucket_start_prime;
	uint32 num_bucket_primes;
    uint64 bitmap_lower_bound;
	uint32 bitmap_start_id;
	uint32 num_bitmap_primes;

	uint64 lowlimit;
	uint64 highlimit;
	uint64 numlinebytes;
	uint32 numclasses;
	uint32 *rclass;
	uint32 *special_count;
	uint32 num_special_bins;
	uint8 **lines;
	uint32 bucket_alloc;
	uint32 large_bucket_alloc;
	uint64 num_found;
#if defined(bitmap_BUCKET)
	soe_bitmap_p *bitmap_data;
	int **bitmap_ptrs;
#endif
	int only_count;
	mpz_t *offset;
	int sieve_range;
	uint64 min_sieved_val;

    // presieving stuff
    int presieve_max_id;

    // for small ranges we have only 2 residue classes which
    // is not enough calls to get_offsets() to justify the
    // Montgomery arithmetic setup costs (reduction for each prime
    // is only performed twice).
    int use_monty;

} soe_staticdata_t;

typedef struct
{
	uint64 *pbounds;
	uint32 *offsets;
	uint64 lblk_b;
	uint64 ublk_b;
	uint64 blk_b_sqrt;
	uint32 bucket_depth;
    uint32 blockstart;
    uint32 blockstop;

	uint32 bucket_alloc;
	uint32 *bucket_hits;
    uint64 **sieve_buckets;
	
	uint32 *special_count;
	uint32 num_special_bins;

	uint32 **large_sieve_buckets;
	uint32 *large_bucket_hits;
	uint32 bucket_alloc_large;
	
	uint64 *primes;
	uint32 largep_offset;
	uint64 min_sieved_val;

    // presieving stuff
    uint32 *presieve_scratch;

} soe_dynamicdata_t;

typedef struct {
	soe_dynamicdata_t ddata;
	soe_staticdata_t sdata;
	uint64 linecount;
	uint32 current_line;

    int tindex;
    int tstartup;

	// start and stop for computing roots
	uint32 startid, stopid;

	// stuff for computing PRPs
	mpz_t offset, lowlimit, highlimit, tmpz;

	/* fields for thread pool synchronization */
	volatile enum soe_command command;

#ifdef USE_SOE_THREADPOOL
    /* fields for thread pool synchronization */
    volatile int *thread_queue, *threads_waiting;

#if defined(WIN32) || defined(_WIN64)
    HANDLE thread_id;
    HANDLE run_event;

    HANDLE finish_event;
    HANDLE *queue_event;
    HANDLE *queue_lock;

#else
    pthread_t thread_id;
    pthread_mutex_t run_lock;
    pthread_cond_t run_cond;

    pthread_mutex_t *queue_lock;
    pthread_cond_t *queue_cond;
#endif

#else

#if defined(WIN32) || defined(_WIN64)
    HANDLE thread_id;
    HANDLE run_event;
    HANDLE finish_event;
#else
    pthread_t thread_id;
    pthread_mutex_t run_lock;
    pthread_cond_t run_cond;
#endif

#endif

} thread_soedata_t;

// for use with threadpool
typedef struct
{
    soe_staticdata_t *sdata;
    thread_soedata_t *ddata;
} soe_userdata_t;

static __inline uint32 redc_loc(uint64 x, uint32 pinv, uint32 p)
{
    uint32 m = (uint32)x * pinv;
    x += (uint64)m * (uint64)p;
    m = x >> 32;
    if (m >= p) m -= p;
    return m;
}

static __inline uint32 to_monty_loc(uint32 x, uint32 r2, uint32 pinv, uint32 p)
{
    uint64 t = (uint64)x * (uint64)r2;
    return redc_loc(t, pinv, p);
}


// top level sieving code
uint64 spSOE(uint32 *sieve_p, uint32 num_sp, mpz_t *offset, 
	uint64 lowlimit, uint64 *highlimit, int count, uint64 *primes);

// thread ready sieving functions
void sieve_line(thread_soedata_t *thread_data);
void sieve_line_avx2_32k(thread_soedata_t *thread_data);
void sieve_line_avx2_64k(thread_soedata_t *thread_data);
void sieve_line_avx2_128k(thread_soedata_t *thread_data);
void sieve_line_avx2_256k(thread_soedata_t *thread_data);
void sieve_line_avx512_32k(thread_soedata_t *thread_data);
void sieve_line_avx512_64k(thread_soedata_t *thread_data);
void sieve_line_avx512_128k(thread_soedata_t *thread_data);
void sieve_line_avx512_256k(thread_soedata_t *thread_data);
void sieve_line_avx512_512k(thread_soedata_t *thread_data);
void(*sieve_line_ptr)(thread_soedata_t *);


uint64 count_line(soe_staticdata_t *sdata, uint32 current_line);
void count_line_special(thread_soedata_t *thread_data);
uint32 compute_32_bytes(soe_staticdata_t *sdata,
    uint32 pcount, uint64 *primes, uint64 byte_offset);
uint64 primes_from_lineflags(soe_staticdata_t *sdata, thread_soedata_t *thread_data,
	uint32 start_count, uint64 *primes);
void get_offsets(thread_soedata_t *thread_data);
void getRoots(soe_staticdata_t *sdata, thread_soedata_t *thread_data);
void stop_soe_worker_thread(thread_soedata_t *t);
void start_soe_worker_thread(thread_soedata_t *t);
#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI soe_worker_thread_main(LPVOID thread_data);
#else
void *soe_worker_thread_main(void *thread_data);
#endif

// routines for finding small numbers of primes; seed primes for main SOE
uint32 tiny_soe(uint32 limit, uint32 *primes);

void test_soe(int upper);

// interface functions
uint64 *GetPRIMESRange(uint32 *sieve_p, uint32 num_sp, 
	mpz_t *offset, uint64 lowlimit, uint64 highlimit, uint64 *num_p);
uint64 *soe_wrapper(uint32 *sieve_p, uint32 num_sp, 
	uint64 lowlimit, uint64 highlimit, int count, uint64 *num_p);
uint64 *sieve_to_depth(uint32 *seed_p, uint32 num_sp, 
	mpz_t lowlimit, mpz_t highlimit, int count, int num_witnesses, uint64 *num_p);

// misc and helper functions
uint64 estimate_primes_in_range(uint64 lowlimit, uint64 highlimit);
void get_numclasses(uint64 highlimit, uint64 lowlimit, soe_staticdata_t *sdata);
int check_input(uint64 highlimit, uint64 lowlimit, uint32 num_sp, uint32 *sieve_p,
	soe_staticdata_t *sdata, mpz_t offset);
uint64 init_sieve(soe_staticdata_t *sdata);
void set_bucket_depth(soe_staticdata_t *sdata);
uint64 alloc_threaddata(soe_staticdata_t *sdata, thread_soedata_t *thread_data);
void do_soe_sieving(soe_staticdata_t *sdata, thread_soedata_t *thread_data, int count);
void finalize_sieve(soe_staticdata_t *sdata,
	thread_soedata_t *thread_data, int count, uint64 *primes);

void pre_sieve(soe_dynamicdata_t *ddata, soe_staticdata_t *sdata, uint8 *flagblock);
void pre_sieve_avx2(soe_dynamicdata_t *ddata, soe_staticdata_t *sdata, uint8 *flagblock);
void pre_sieve_avx512(soe_dynamicdata_t *ddata, soe_staticdata_t *sdata, uint8 *flagblock);
void (*pre_sieve_ptr)(soe_dynamicdata_t *, soe_staticdata_t *, uint8 *);

uint32 compute_8_bytes(soe_staticdata_t *sdata,
    uint32 pcount, uint64 *primes, uint64 byte_offset);
uint32 compute_8_bytes_bmi2(soe_staticdata_t *sdata,
    uint32 pcount, uint64 *primes, uint64 byte_offset);
uint32 (*compute_8_bytes_ptr)(soe_staticdata_t *, uint32, uint64 *, uint64);

// misc
void primesum(uint64 lower, uint64 upper);
void primesum_check12(uint64 lower, uint64 upper, uint64 startmod, z *squaresum, z *sum);
void primesum_check3(uint64 lower, uint64 upper, uint64 startmod, z *sum);

//masks for removing or reading single bits in a byte.  nmasks are simply
//the negation of these masks, and are filled in within the spSOE function.
static const uint8 masks[8] = {0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f};
static const uint8 nmasks[8] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };
uint32 max_bucket_usage;
uint64 GLOBAL_OFFSET;
int NO_STORE;

uint32 SOEBLOCKSIZE;
uint32 FLAGSIZE;
uint32 FLAGSIZEm1;
uint32 FLAGBITS;
uint32 BUCKETSTARTI;

#endif // YAFU_SOE_H
