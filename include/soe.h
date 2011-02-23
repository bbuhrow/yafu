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

#include "yafu.h"
#include "arith.h"	//needed for spGCD and modinv_1
#include "util.h"

#ifdef YAFU_64K
#define BLOCKSIZE 65536
#define FLAGSIZE 524288
#define FLAGSIZEm1 524287
#define FLAGBITS 19
#define BUCKETSTARTP 393216 //524288  //786432 
#define BUCKETSTARTI 33335 //43390  //62946 
#else
#define BLOCKSIZE 32768
#define FLAGSIZE 262144
#define FLAGSIZEm1 262143
#define FLAGBITS 18
#define BUCKETSTARTP 393216 
#define BUCKETSTARTI 33335
#endif

#define BUCKET_BUFFER 1024
#define BITSINBYTE 8
#define MAXSIEVEPRIMECOUNT 100000000	//# primes less than ~2e9: limit of 2e9^2 = 4e18
//#define INPLACE_BUCKET 1

//#define DO_SPECIAL_COUNT

enum soe_command {
	SOE_COMMAND_INIT,
	SOE_COMMAND_WAIT,
	SOE_COMMAND_SIEVE_AND_COUNT,
	SOE_COMMAND_SIEVE_AND_COMPUTE,
	SOE_COMMAND_END
};

typedef struct
{
	uint16 loc;
	uint8 mask;
	uint8 bnum;
} soe_bucket_t_old;

//create num_blocks * num_residues linked lists, one for each combination of residue and block.
//Then, for each residue, for each block, iterate through the appropriate list, AND its location, 
//compute the next block, residue, and location, and add it to the end of the appropriate linked list.
//******** this is only valid for primes with strides long enough to 
//******** make sure that they hit a new block each step in residue space

//update: nope doesn't work.  in the linesieve, we sieve all blocks of each residue class once.  but
//with the inplace sieve, the residue class changes as the prime is advanced through the real number
//line.  thus, we would need to sieve each block with each residue class, each of which uses a different
//segment of memory.  this would be hugely inefficient.
//example.  say we sieve all blocks in residue class 1, and then move on to residue class 7.  some primes
//in that class will move back to residue class 1 when advanced, but we've already sieved class 1, so they
//are missed.  thus class 7 would also need to sieve class 1, and would need the memory segments (blocks) from
//the class 1 line.  the last class would need to sieve all other classes and would need all other lines' 
//memory segments and would thus take forever.
typedef struct
{
	uint32 prime;		//the value of this prime
	uint32 id;			//the index of the prime in the bucket_prime array
	uint32 loc;				//location of the hit in the block
	uint32 next;			//relative offset to the next bucket_prime_t with the same residue and block
	uint32 steps;
	uint16 res;				//the residue class of this prime
} bucket_prime_t;			

typedef struct
{
	uint32 root;
	uint32 prime;
} soe_bucket_t;

typedef struct
{
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
	uint64 lowlimit;
	uint64 highlimit;
	uint64 numlinebytes;
	uint32 numclasses;
	uint32 *rclass;
	uint32 *special_count;
	uint32 num_special_bins;

#ifdef INPLACE_BUCKET
	bucket_prime_t **listptrs;		//array of pointers to bucket_prime_t's, one for each block and residue
	bucket_prime_t *bucket_primes;	
#endif
	uint32 inplace_startindex;
	uint32 inplace_startprime;
	uint32 *valid_residue;

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
	soe_bucket_t **sieve_buckets;
	
	uint32 *special_count;
	uint32 num_special_bins;

	uint32 **large_sieve_buckets;
	uint32 *large_bucket_hits;
	uint32 bucket_alloc_large;

	uint8 *line;
	uint64 *primes;
	uint32 largep_offset;

} soe_dynamicdata_t;

typedef struct {
	soe_dynamicdata_t ddata;
	soe_staticdata_t sdata;
	//uint8 **line;
	uint64 linecount;
	uint32 current_line;

	/* fields for thread pool synchronization */
	volatile enum soe_command command;

#if defined(WIN32) || defined(_WIN64)
	HANDLE thread_id;
	HANDLE run_event;
	HANDLE finish_event;
#else
	pthread_t thread_id;
	pthread_mutex_t run_lock;
	pthread_cond_t run_cond;
#endif

} thread_soedata_t;

//top level
uint64 spSOE(uint64 *primes, uint64 lowlimit, uint64 *highlimit, int count);

//thread ready sieving functions
void sieve_line(thread_soedata_t *thread_data);
void count_line(thread_soedata_t *thread_data);
void count_line_special(thread_soedata_t *thread_data);
void primes_from_lineflags(thread_soedata_t *thread_data);
void get_offsets(thread_soedata_t *thread_data);
void getRoots(soe_staticdata_t *sdata);
void stop_soe_worker_thread(thread_soedata_t *t, uint32 is_master_thread);
void start_soe_worker_thread(thread_soedata_t *t, uint32 is_master_thread);
#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI soe_worker_thread_main(LPVOID thread_data);
#else
void *soe_worker_thread_main(void *thread_data);
#endif

//routines for finding small numbers of primes; seed primes for main SOE
uint32 tiny_soe(uint32 limit, uint32 *primes);

//wrapper function
void GetPRIMESRange(uint64 lowlimit, uint64 highlimit);

//sieve an interval of numbers to a specified bit depth
int sieve_to_bitdepth(z *start, z *stop, int depth, uint32 *offsets);

//wrapper function
uint64 soe_wrapper(uint64 lowlimit, uint64 highlimit, int count);

//get the next prime after or before the one specified
void zNextPrime(z *n, z *p, int dir);

//misc
void primesum(uint64 lower, uint64 upper);
void primesum_check12(uint64 lower, uint64 upper, uint64 startmod, z *squaresum, z *sum);
void primesum_check3(uint64 lower, uint64 upper, uint64 startmod, z *sum);

//masks for removing or reading single bits in a byte.  nmasks are simply
//the negation of these masks, and are filled in within the spSOE function.
static const uint8 masks[8] = {0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f};
uint8 nmasks[8];
uint32 max_bucket_usage;

