/*
MIT License

Copyright (c) 2021 Ben Buhrow

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SOE_H
#define SOE_H

#include <stdlib.h>
#include <stdint.h>

#if defined(_WIN32)

//#include <windows.h>
//#include <process.h>

#endif

#include "gmp.h"
#include "ytools.h"

#define USE_SOE_THREADPOOL


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

typedef struct
{
	//uint32_t prime;		// the prime, so that we don't have to also look in the
						// main prime array
	uint32_t bitloc;		// bit location of the current hit
	uint32_t next_pid;	// index of the next prime that hits in the current sieve
	uint32_t p_div;		// prime / prodN
	uint8_t p_mod;		// prime % prodN
	uint8_t eacc;			// accumulated error
} soe_bitmap_p;

typedef struct
{
    int VFLAG;
    int THREADS;
    int sync_count;
	int witnesses;
	int userclasses;
	uint32_t *sieve_p;
    uint32_t num_sp;
	uint32_t alloc_sp;
	int *root;
	uint32_t *lower_mod_prime;

    uint32_t *pinv;       // montgomery inverse
    uint32_t *r2modp;     // to go out of montgomery rep

	uint64_t blk_r;
	uint64_t blocks;
	uint64_t partial_block_b;
	uint64_t prodN;
	uint64_t startprime;
	uint64_t orig_hlimit;
	uint64_t orig_llimit;
	uint64_t pbound;
	uint64_t pboundi;

	uint32_t bucket_start_id;
	uint64_t large_bucket_start_prime;
	uint32_t num_bucket_primes;
    uint64_t bitmap_lower_bound;
	uint32_t bitmap_start_id;
	uint32_t num_bitmap_primes;

	uint64_t lowlimit;
	uint64_t highlimit;
	uint64_t numlinebytes;
	uint32_t numclasses;
	uint32_t *rclass;
	uint32_t *special_count;
	uint32_t num_special_bins;
	uint8_t **lines;
	uint32_t bucket_alloc;
	uint32_t large_bucket_alloc;
	uint64_t num_found;
#if defined(bitmap_BUCKET)
	soe_bitmap_p *bitmap_data;
	int **bitmap_ptrs;
#endif
	int only_count;
	int analysis;		// 1 == compute all primes, 2 == twins, 3 == gaps
	int gapmin;			// minimum gap size to search for with analysis=3;
	int primeconst;		// prime constellation size to search for (default 2 == twins)
	mpz_t *offset;
	int sieve_range;
	uint64_t min_sieved_val;

    // presieving stuff
    int presieve_max_id;

    // for small ranges we have only 2 residue classes which
    // is not enough calls to get_offsets() to justify the
    // Montgomery arithmetic setup costs (reduction for each prime
    // is only performed twice).
    int use_monty;
	int is_main_sieve;

    // masks for removing or reading single bits in a byte.  nmasks are simply
    // the negation of these masks, and are filled in within the spSOE function.
    uint8_t masks[8];
    uint8_t nmasks[8];
    uint32_t masks32[32];
    uint32_t nmasks32[32];
    uint32_t max_bucket_usage;
    uint64_t GLOBAL_OFFSET;
    int NO_STORE;
    uint32_t SOEBLOCKSIZE;
    uint32_t FLAGSIZE;
    uint32_t FLAGSIZEm1;
    uint32_t FLAGBITS;
    uint32_t BUCKETSTARTI;
	int has_avx2;
	int has_bmi2;
	int has_bmi1;
	int has_avx512f;

} soe_staticdata_t;

typedef struct
{
	uint64_t *pbounds;
	uint32_t *offsets;
	uint64_t lblk_b;
	uint64_t ublk_b;
	uint64_t blk_b_sqrt;
	uint32_t bucket_depth;
    uint32_t blockstart;
    uint32_t blockstop;

	uint32_t bucket_alloc;
	uint32_t *bucket_hits;
    uint64_t **sieve_buckets;
	
	uint32_t *special_count;
	uint32_t num_special_bins;

	uint32_t **large_sieve_buckets;
	uint32_t *large_bucket_hits;
	uint32_t bucket_alloc_large;
	
	uint64_t *primes;
	uint32_t largep_offset;
	uint64_t min_sieved_val;

    // presieving stuff
    uint32_t *presieve_scratch;

	uint8_t* analysis_carry_data;		// data to bridge between analysis boundaries

} soe_dynamicdata_t;

typedef struct {
	soe_dynamicdata_t ddata;
	soe_staticdata_t sdata;
	uint64_t linecount;
	uint32_t current_line;

    int tindex;
    int tstartup;

	// start and stop for computing roots
	uint32_t startid, stopid;

	// stuff for computing PRPs
	mpz_t offset, lowlimit, highlimit, tmpz;

	/* fields for thread pool synchronization */
//	volatile enum soe_command command;
//
//#ifdef USE_SOE_THREADPOOL
//    /* fields for thread pool synchronization */
//    volatile int *thread_queue, *threads_waiting;
//
//#if (defined(WIN32) || defined(_WIN64))
//    HANDLE thread_id;
//    HANDLE run_event;
//
//    HANDLE finish_event;
//    HANDLE *queue_event;
//    HANDLE *queue_lock;
//
//#else
//    pthread_t thread_id;
//    pthread_mutex_t run_lock;
//    pthread_cond_t run_cond;
//
//    pthread_mutex_t *queue_lock;
//    pthread_cond_t *queue_cond;
//#endif
//
//#else
//
//#if defined(WIN32) || defined(_WIN64)
//    HANDLE thread_id;
//    HANDLE run_event;
//    HANDLE finish_event;
//#else
//    pthread_t thread_id;
//    pthread_mutex_t run_lock;
//    pthread_cond_t run_cond;
//#endif
//
//#endif

} thread_soedata_t;

// for use with threadpool
typedef struct
{
    soe_staticdata_t *sdata;
    thread_soedata_t *ddata;
} soe_userdata_t;


// interface functions
extern soe_staticdata_t* soe_init(int vflag, int threads, int blocksize);
extern void soe_finalize(soe_staticdata_t* sdata);
extern uint64_t* soe_wrapper(soe_staticdata_t* sdata, uint64_t lowlimit, uint64_t highlimit,
    int count, uint64_t* num_p, int PRIMES_TO_FILE, int PRIMES_TO_SCREEN);
extern uint64_t* sieve_to_depth(soe_staticdata_t* sdata,
    mpz_t lowlimit, mpz_t highlimit, int count, int num_witnesses, 
    uint64_t sieve_limit, uint64_t* num_p,
    int PRIMES_TO_FILE, int PRIMES_TO_SCREEN);


#endif // #ifndef SOE_H
