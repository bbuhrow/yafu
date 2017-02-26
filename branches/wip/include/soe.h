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

#define USE_SOE_THREADPOOL

#ifdef _MSC_VER
// optionally define this or not depending on whether your hardware supports it.
// if defined, compile the sse41 functions into the fat binary.  the global
// flag HAS_SSE41 is set at runtime on compatible hardware to enable the functions
// to be used.  For gcc and mingw64 builds, USE_SSE41 is enabled in the makefile.
//#define USE_SSE41 1
//#define USE_AVX2 1
#endif

#define BLOCKSIZE 32768
#define FLAGSIZE 262144
#define FLAGSIZEm1 262143
#define FLAGBITS 18
#define BUCKETSTARTP 393216 
#define BUCKETSTARTI 33335
#define BITSINBYTE 8
#define MAXSIEVEPRIMECOUNT 100000000	//# primes less than ~2e9: limit of 2e9^2 = 4e18
//#define INPLACE_BUCKET 1
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
// for storage of presieving lists from prime index 24 to 40 (97 to 173 inclusive)
#ifdef __GNUC__
__attribute__((aligned(64))) uint64 presieve_largemasks[16][173][4];
__attribute__((aligned(64))) uint32 presieve_steps[32];
__attribute__((aligned(64))) uint32 presieve_primes[32];
__attribute__((aligned(64))) uint32 presieve_p1[32];
#else
__declspec(align(64)) uint64 presieve_largemasks[16][173][4];
__declspec(align(64)) uint32 presieve_steps[32];
__declspec(align(64)) uint32 presieve_primes[32];
__declspec(align(64)) uint32 presieve_p1[32];
#endif
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
} soe_inplace_p;

typedef struct
{
    int sync_count;
	uint32 *sieve_p;
	int *root;
	uint32 *lower_mod_prime;
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
	uint32 inplace_start_id;
	uint32 num_inplace_primes;

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
#if defined(INPLACE_BUCKET)
	soe_inplace_p *inplace_data;
	int **inplace_ptrs;
#endif
	int only_count;
	mpz_t *offset;
	int sieve_range;
	uint64 min_sieved_val;

    // presieving stuff
    int presieve_max_id;

} soe_staticdata_t;

typedef struct
{
	uint64 *pbounds;
	uint32 *offsets;
	uint64 lblk_b;
	uint64 ublk_b;
	uint64 blk_b_sqrt;
	uint32 bucket_depth;

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

// top level sieving code
uint64 spSOE(uint32 *sieve_p, uint32 num_sp, mpz_t *offset, 
	uint64 lowlimit, uint64 *highlimit, int count, uint64 *primes);

// thread ready sieving functions
void sieve_line(thread_soedata_t *thread_data);
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
void (*pre_sieve_ptr)(soe_dynamicdata_t *, soe_staticdata_t *, uint8 *);

uint32 compute_8_bytes(soe_staticdata_t *sdata,
    uint32 pcount, uint64 *primes, uint64 byte_offset);
uint32 compute_8_bytes_avx2(soe_staticdata_t *sdata,
    uint32 pcount, uint64 *primes, uint64 byte_offset);
uint32 (*compute_8_bytes_ptr)(soe_staticdata_t *, uint32, uint64 *, uint64);

// misc
void primesum(uint64 lower, uint64 upper);
void primesum_check12(uint64 lower, uint64 upper, uint64 startmod, z *squaresum, z *sum);
void primesum_check3(uint64 lower, uint64 upper, uint64 startmod, z *sum);

//masks for removing or reading single bits in a byte.  nmasks are simply
//the negation of these masks, and are filled in within the spSOE function.
static const uint8 masks[8] = {0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f};
uint8 nmasks[8];
uint32 max_bucket_usage;
uint64 GLOBAL_OFFSET;
int NO_STORE;

//progression of residue classes. each row
//is for a different prime mod prodN
//static int residue_pattern_mod30[8][8] = {
//	{1, 7, 11, 13, 17, 19, 23, 29},
//	{7, 19, 17, 1, 29, 13, 11, 23},
//	{11, 17, 1, 23, 7, 29, 13, 19},
//	{13, 1, 23, 19, 11, 7, 29, 17},
//	{17, 29, 7, 11, 19, 23, 1, 13},
//	{19, 13, 29, 7, 23, 1, 17, 11},
//	{23, 11, 13, 29, 1, 17, 19, 7},
//	{29, 23, 19, 17, 13, 11, 7, 1}};

//static int steps_mod30[8] = {6,4,2,4,2,4,6,2};

//same info differently organized.  row/col
//returns the next residue class when row
//is the current prime mod prodN and col is the
//current res class.
//static int residue_pattern_mod30[8][8] = {
//	{7, 11, 13, 17, 19, 23, 29, 1},
//	{29, 19, 23, 11, 1, 17, 7, 13},
//	{23, 29, 17, 19, 1, 11, 7, 13},
//	{23, 29, 7, 1, 13, 11, 19, 17},
//	{13, 11, 19, 17, 29, 23, 1, 7},
//	{17, 23, 19, 29, 11, 13, 1, 7},
//	{17, 23, 13, 29, 19, 7, 11, 1},
//	{29, 1, 7, 11, 13, 17, 19, 23}};
//
////steps offset
//static int steps_mod30[8][8] = {
//	{6, 4, 2, 4, 2, 4, 6, 2},
//	{4, 6, 6, 4, 2, 4, 2, 2},
//	{2, 2, 6, 6, 4, 2, 4, 4},
//	{4, 4, 2, 6, 2, 4, 2, 6},
//	{6, 2, 4, 2, 6, 2, 4, 4},
//	{4, 4, 2, 4, 6, 6, 2, 2},
//	{2, 2, 4, 2, 4, 6, 6, 4},
//	{2, 6, 4, 2, 4, 2, 4, 6}};
//
//// residue lookup
//static int resID_mod30[30] = {
//	-1, 0,-1,-1,-1,
//	-1,-1, 1,-1,-1,
//	-1, 2,-1, 3,-1,
//	-1,-1, 4,-1, 5,
//	-1,-1,-1, 6,-1,
//	-1,-1,-1,-1, 7};
//
//// next valid class lookup
//static int next_mod30[8][30] =	{
//	{257,1,1287,1031,775,519,263,7,779,523,267,11,269,13,785,529,273,17,275,19,791,535,279,23,1309,1053,797,541,285,29},
//	{263,1,791,529,267,531,269,7,797,535,273,11,275,13,1299,541,279,17,1303,19,779,1043,285,23,257,1047,785,523,787,29},
//	{267,1,269,1041,775,1043,273,7,275,513,781,11,279,13,785,519,787,17,285,19,257,525,1297,23,1299,529,263,531,769,29},
//	{269,1,779,541,273,513,275,7,785,1025,279,11,1297,13,791,523,285,17,257,19,797,529,769,23,263,1041,1281,535,267,29},
//	{273,1,275,519,1309,1037,279,7,797,525,769,11,285,13,257,531,775,17,1293,19,263,1053,781,23,267,541,269,513,787,29},
//	{275,1,797,523,279,525,1291,7,1293,529,285,11,257,13,779,535,781,17,263,19,785,541,267,23,269,1035,791,1037,273,29},
//	{279,1,779,531,781,1031,285,7,257,1035,787,11,1287,13,263,513,1291,17,267,19,269,519,769,23,273,523,275,525,775,29},
//	{285,1,257,513,769,1025,1281,7,263,519,775,11,267,13,269,525,781,17,273,19,275,531,787,23,279,535,791,1047,1303,29}};
//
//

#endif // YAFU_SOE_H
