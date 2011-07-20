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
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/


/*
implements a very fast sieve of erathostenes.  Speed enhancements
include a variable mod30 or mod210 wheel, bit packing
cache blocking, various methods of loop unrolling, batch sieving
of small primes, and threading.  

cache blocking means that the sieve interval is split up into blocks,
and each block is sieved separately by all primes that are needed.
the size of the blocks and of the information about the sieving
primes and offsets into the next block are carefully chosen so
that everything fits into cache while sieving a block.

bit packing means that a number is represented as prime or not with
a single bit, rather than a byte or word.  this greatly reduces 
the storage requirements of the flags, which makes cache blocking 
more effective.

the mod30 or mod210 wheel means that depending
on the size of limit, either the primes 2,3, and 5 are presieved
or the primes 2,3,5 and 7 are presieved.  in other words, only the 
numbers not divisible by these small primes are 'written down' with flags.
this reduces the storage requirements for the prime flags, which in 
turn makes the cache blocking even more effective.  It also directly 
increases the speed because the slowest sieving steps have been removed.

*/

#include "yafu.h"
#include "soe.h"
#include "common.h"
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
#define MAXSIEVEPRIMECOUNT 189961812	//# primes less than 4e9: limit of 4e9^2 = 1.6e19

static int _SMALL_P[9] = {5,7,11,13,17,19,23,29,31};
static int _64_MOD_P[9] = {4,1,9,12,13,7,18,6,2};

static uint64 _5_MASKS[5] = {
	0xef7bdef7bdef7bde,
	0xdef7bdef7bdef7bd,
	0xbdef7bdef7bdef7b,
	0x7bdef7bdef7bdef7,
	0xf7bdef7bdef7bdef};

static uint64 _7_MASKS[7] = {
	0x7efdfbf7efdfbf7e,
	0xfdfbf7efdfbf7efd,
	0xfbf7efdfbf7efdfb,
	0xf7efdfbf7efdfbf7,
	0xefdfbf7efdfbf7ef,
	0xdfbf7efdfbf7efdf,
	0xbf7efdfbf7efdfbf};

static uint64 _11_MASKS[11] = {
	0xff7feffdffbff7fe,
	0xfeffdffbff7feffd,
	0xfdffbff7feffdffb,
	0xfbff7feffdffbff7,
	0xf7feffdffbff7fef,
	0xeffdffbff7feffdf,
	0xdffbff7feffdffbf,
	0xbff7feffdffbff7f,
	0x7feffdffbff7feff,
	0xffdffbff7feffdff,
	0xffbff7feffdffbff};

static uint64 _13_MASKS[13] = {
	0xffefff7ffbffdffe,
	0xffdffefff7ffbffd,
	0xffbffdffefff7ffb,
	0xff7ffbffdffefff7,
	0xfefff7ffbffdffef,
	0xfdffefff7ffbffdf,
	0xfbffdffefff7ffbf,
	0xf7ffbffdffefff7f,
	0xefff7ffbffdffeff,
	0xdffefff7ffbffdff,
	0xbffdffefff7ffbff,
	0x7ffbffdffefff7ff,
	0xfff7ffbffdffefff};

static uint64 _17_MASKS[17] = {
	0xfff7fffbfffdfffe,
	0xffeffff7fffbfffd,
	0xffdfffeffff7fffb,
	0xffbfffdfffeffff7,
	0xff7fffbfffdfffef,
	0xfeffff7fffbfffdf,
	0xfdfffeffff7fffbf,
	0xfbfffdfffeffff7f,
	0xf7fffbfffdfffeff,
	0xeffff7fffbfffdff,
	0xdfffeffff7fffbff,
	0xbfffdfffeffff7ff,
	0x7fffbfffdfffefff,
	0xffff7fffbfffdfff,
	0xfffeffff7fffbfff,
	0xfffdfffeffff7fff,
	0xfffbfffdfffeffff};

static uint64 _19_MASKS[19] = {
	0xfdffffbffff7fffe,
	0xfbffff7fffeffffd,
	0xf7fffeffffdffffb,
	0xeffffdffffbffff7,
	0xdffffbffff7fffef,
	0xbffff7fffeffffdf,
	0x7fffeffffdffffbf,
	0xffffdffffbffff7f,
	0xffffbffff7fffeff,
	0xffff7fffeffffdff,
	0xfffeffffdffffbff,
	0xfffdffffbffff7ff,
	0xfffbffff7fffefff,
	0xfff7fffeffffdfff,
	0xffeffffdffffbfff,
	0xffdffffbffff7fff,
	0xffbffff7fffeffff,
	0xff7fffeffffdffff,
	0xfeffffdffffbffff};

static uint64 _23_MASKS[23] = {
	0xffffbfffff7ffffe,
	0xffff7ffffefffffd,
	0xfffefffffdfffffb,
	0xfffdfffffbfffff7,
	0xfffbfffff7ffffef,
	0xfff7ffffefffffdf,
	0xffefffffdfffffbf,
	0xffdfffffbfffff7f,
	0xffbfffff7ffffeff,
	0xff7ffffefffffdff,
	0xfefffffdfffffbff,
	0xfdfffffbfffff7ff,
	0xfbfffff7ffffefff,
	0xf7ffffefffffdfff,
	0xefffffdfffffbfff,
	0xdfffffbfffff7fff,
	0xbfffff7ffffeffff,
	0x7ffffefffffdffff,
	0xfffffdfffffbffff,
	0xfffffbfffff7ffff,
	0xfffff7ffffefffff,
	0xffffefffffdfffff,
	0xffffdfffffbfffff};

static uint64 _29_MASKS[29] = {
	0xfbffffffdffffffe,
	0xf7ffffffbffffffd,
	0xefffffff7ffffffb,
	0xdffffffefffffff7,
	0xbffffffdffffffef,
	0x7ffffffbffffffdf,
	0xfffffff7ffffffbf,
	0xffffffefffffff7f,
	0xffffffdffffffeff,
	0xffffffbffffffdff,
	0xffffff7ffffffbff,
	0xfffffefffffff7ff,
	0xfffffdffffffefff,
	0xfffffbffffffdfff,
	0xfffff7ffffffbfff,
	0xffffefffffff7fff,
	0xffffdffffffeffff,
	0xffffbffffffdffff,
	0xffff7ffffffbffff,
	0xfffefffffff7ffff,
	0xfffdffffffefffff,
	0xfffbffffffdfffff,
	0xfff7ffffffbfffff,
	0xffefffffff7fffff,
	0xffdffffffeffffff,
	0xffbffffffdffffff,
	0xff7ffffffbffffff,
	0xfefffffff7ffffff,
	0xfdffffffefffffff};

//masks for removing or reading single bits in a byte.  nmasks are simply
//the negation of these masks, and are filled in within the spSOE function.
const uint8 masks[8] = {0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f};
uint8 nmasks[8];

uint32 max_bucket_usage;

#if defined(MSC_ASM64X) || defined(GCC_ASM64X)
	//enable loop unrolling only if we have extra registers
	#define SOE_UNROLL_LOOPS
#endif

void start_soe_worker_thread(thread_soedata_t *t, uint32 is_master_thread) 
{

	/* create a thread that will process a sieve line. The last line does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread) {
		return;
	}

	t->command = SOE_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
	t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
	t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	t->thread_id = CreateThread(NULL, 0, soe_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_create(&t->thread_id, NULL, soe_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != SOE_COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
}

void stop_soe_worker_thread(thread_soedata_t *t, uint32 is_master_thread)
{
	if (is_master_thread) {
		return;
	}

	t->command = SOE_COMMAND_END;
#if defined(WIN32) || defined(_WIN64)
	SetEvent(t->run_event);
	WaitForSingleObject(t->thread_id, INFINITE);
	CloseHandle(t->thread_id);
	CloseHandle(t->run_event);
	CloseHandle(t->finish_event);
#else
	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_join(t->thread_id, NULL);
	pthread_cond_destroy(&t->run_cond);
	pthread_mutex_destroy(&t->run_lock);
#endif
}

#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI soe_worker_thread_main(LPVOID thread_data) {
#else
void *soe_worker_thread_main(void *thread_data) {
#endif
	thread_soedata_t *t = (thread_soedata_t *)thread_data;

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == SOE_COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */

		if (t->command == SOE_COMMAND_SIEVE_AND_COUNT)
		{
			sieve_line(t);
#ifdef DO_SPECIAL_COUNT
			count_line_special(t);
#else
			count_line(t);
#endif
		}
		else if (t->command == SOE_COMMAND_SIEVE_AND_COMPUTE)
		{
			sieve_line(t);
			primes_from_lineflags(t);
		}
		else if (t->command == SOE_COMMAND_END)
			break;

		/* signal completion */

		t->command = SOE_COMMAND_WAIT;
#if defined(WIN32) || defined(_WIN64)
		SetEvent(t->finish_event);
#else
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}

uint64 spSOE(uint64 *primes, uint64 lowlimit, uint64 *highlimit, int count)
{
	/*
	finds primes up to 'limit' using the Sieve of Erathostenes
	using the mod30 or mod210 wheel, bit packing, and block sieving

	return the number of primes found less than 'limit'
	'limit' is modified to an integer slightly larger, due to 
	cache blocking and the wheel - it is easier to finish a block
	rather than check in mid loop to see if the limit has been
	exceeded.

	if count == 1, then the primes are simply counted, and not 
	explicitly calculated and saved in *primes.
	*/

	//*********************** VARIABLES ******************************//
	
	//variables used for the wheel
	uint64 numclasses,prodN,startprime;
	uint32 bucket_depth;

	//variables for blocking up the line structures
	uint64 numflags, numbytes, numlinebytes, bucket_alloc;

	//misc
	uint64 i,j,k,it=0,num_p=0;
	uint64 i1,i2,i3;
	uint32 sp;
	int pchar;
	uint64 allocated_bytes = 0;

	//structure of static info
	soe_staticdata_t sdata;

	//thread data holds all data needed during sieving
	thread_soedata_t *thread_data;		//an array of thread data objects
	uint64 *locprimes, *mergeprimes;

	//timing
	double t;
	struct timeval tstart, tstop;
	TIME_DIFF *	difference;

	//*********************** CODE ******************************//

	sdata.orig_hlimit = *highlimit;
	sdata.orig_llimit = lowlimit;

	if (*highlimit - lowlimit < 1000000)
		*highlimit = lowlimit + 1000000;

	if (*highlimit - lowlimit > 1000000000000)
	{
		printf("range too big\n");
		return 0;
	}

	//more efficient to sieve using mod210 when the range is big
	if ((*highlimit - lowlimit) > 400000000000)
	{
		numclasses=5760;
		prodN=30030;
		startprime=6;
	}	
	else if ((*highlimit - lowlimit) > 40000000000)
	{
		numclasses=480;
		prodN=2310;
		startprime=5;
	}	
	else if ((*highlimit - lowlimit) > 4000000000)
	{
		numclasses=48;
		prodN=210;
		startprime=4;		
	}
	else if ((*highlimit - lowlimit) > 100000000)
	{
		numclasses=8;
		prodN=30;
		startprime=3;
	}
	else
	{
		numclasses=2;
		prodN=6;
		startprime=2;
	}

	sdata.numclasses = numclasses;
	sdata.prodN = prodN;
	sdata.startprime = startprime;

	//create the selection masks
	for (i=0;i<BITSINBYTE;i++)
		nmasks[i] = ~masks[i];

	//store the primes used to sieve the rest of the flags
	//with the max sieve range set by the size of uint32, the number of primes
	//needed is fixed.
	//find the bound of primes we'll need to sieve with to get to limit
	sdata.pbound = (uint64)(sqrt((int64)(*highlimit)) + 1);

	//give these some initial storage
	locprimes = (uint64 *)calloc(1,sizeof(uint64));
	mergeprimes = (uint64 *)calloc(1,sizeof(uint64));

	//get primes to sieve with
	if (VFLAG > 2)
		gettimeofday(&tstart, NULL);
	
	if (sdata.pbound > 1000000)
	{
		// we need a lot of sieving primes... recurse using the fast routine
		// get temporary 64 bit prime storage
		// first estimate how many primes we think we'll find
		j = sdata.pbound;
		k = (uint64)((double)j/log((double)j)*1.2);
		locprimes = realloc(locprimes,k * sizeof(uint64));		

		// recursion
		sp = (uint32)spSOE(locprimes,0,&j,0);

		//check if we've found too many
		if (sp > MAXSIEVEPRIMECOUNT)
		{
			printf("input too high\n");
			free(locprimes);
			free(mergeprimes);
			free(sdata.sieve_p);
			return 0;
		}
		else
		{
			//allocate the sieving prime array
			sdata.sieve_p = (uint32 *)malloc(sp * sizeof(uint32));
			allocated_bytes += sp * sizeof(uint32);
			if (sdata.sieve_p == NULL)
			{
				printf("error allocating sieve_p\n");
				exit(-1);
			}
			else
			{
				if (VFLAG > 2)
					printf("allocated %u bytes for sieving primes\n",sp * sizeof(uint32));
			}
		}

		// copy into the 32 bit sieving prime array
		for (k=0; k<sp; k++)
			sdata.sieve_p[k] = (uint32)locprimes[k];
	}
	else	
	{
		//base case, maxP <= 1000000.  Need at most 78498 primes
		sdata.sieve_p = (uint32 *)malloc(78498 * sizeof(uint32));
		allocated_bytes += 78498 * sizeof(uint32);
		if (sdata.sieve_p == NULL)
		{
			printf("error allocating sieve_p\n");
			exit(-1);
		}
		else
		{
			if (VFLAG > 2)
				printf("allocated %u bytes for sieving primes\n",78498 * sizeof(uint32));
		}
		sp = tiny_soe(sdata.pbound,sdata.sieve_p);
	}

	if (count)
	{
		locprimes = (uint64 *)realloc(locprimes,1 * sizeof(uint64));
	}
	else
	{
		//allocate two arrays used for merging together primes found in different lines 
		//that will be used later.
		j = *highlimit;
		k = (uint64)((double)j/log((double)j)*1.2);
		if (VFLAG > 2)
		{
			printf("estimating storage for primes up to %lu\n",*highlimit);
			printf("allocating merge prime storage for %u primes\n",k);
		}
		locprimes = (uint64 *)realloc(locprimes,k * sizeof(uint64));
		mergeprimes = (uint64 *)realloc(mergeprimes,k * sizeof(uint64));
		if (locprimes == NULL)
		{
			printf("could not allocate storage for locprimes\n");
			exit(-1);
		}
		if (mergeprimes == NULL)
		{
			printf("could not allocate storage for mergeprimes\n");
			exit(-1);
		}
	
	}

	if (VFLAG > 2)
	{
		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);

		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);

		printf("elapsed time for seed primes = %6.4f\n",t);
	}
	
	sdata.pboundi = sp;

	//allocate the residue classes.  
	sdata.rclass = (uint32 *)malloc(numclasses * sizeof(uint32));
	allocated_bytes += numclasses * sizeof(uint32);

	//find the residue classes
	k=0;
	for (i=1;i<prodN;i++)
	{
		if (spGCD(i,(fp_digit)prodN) == 1)
		{
			sdata.rclass[k] = (uint32)i;
			//printf("%u ",i);
			k++;
		}
	}

	//temporarily set lowlimit to the first multiple of numclasses*prodN < lowlimit
	lowlimit = (lowlimit/(numclasses*prodN))*(numclasses*prodN);
	sdata.lowlimit = lowlimit;

	//reallocate flag structure for wheel and block sieving
	//starting at lowlimit, we need a flag for every 'numresidues' numbers out of 'prodN' up to 
	//limit.  round limit up to make this a whole number.
	numflags = (*highlimit - lowlimit)/prodN;
	numflags += ((numflags % prodN) != 0);
	numflags *= numclasses;

	//since we can pack 8 flags in a byte, we need numflags/8 bytes allocated.
	numbytes = numflags / BITSINBYTE + ((numflags % BITSINBYTE) != 0);

	//since there are 8 lines to sieve over, each line will contain (numflags/8)/8 bytes
	//so round numflags/8 up to the nearest multiple of 8
	numlinebytes = numbytes/numclasses + ((numbytes % numclasses) != 0);

	//we want an integer number of blocks, so round up to the nearest multiple of blocksize bytes
	i = 0;
	while (1)
	{
		i += BLOCKSIZE;
		if (i > numlinebytes)
			break;
	}
	numlinebytes = i;

	*highlimit = (uint64)((uint64)numlinebytes * (uint64)prodN * (uint64)BITSINBYTE + lowlimit);
	sdata.highlimit = *highlimit;
	sdata.numlinebytes = numlinebytes;

	//a block consists of BLOCKSIZE bytes of flags
	//which holds FLAGSIZE flags.
	sdata.blocks = numlinebytes/BLOCKSIZE;
	sdata.partial_block_b = (numlinebytes % BLOCKSIZE)*BITSINBYTE;

	sdata.blk_r = FLAGSIZE*prodN;
	it=0;

	//this could perhaps be multithreaded, but it probably isn't necessary.
	//find all xGCDs of prime with prodN.  These are used when finding offsets
	sdata.root = (int *)malloc(sdata.pboundi * sizeof(int));
	allocated_bytes += sdata.pboundi * sizeof(uint32);
	if (sdata.root == NULL)
	{
		printf("error allocating roots\n");
		exit(-1);
	}
	else
	{
		if (VFLAG > 2)
			printf("allocated %u bytes for roots\n",sdata.pboundi * sizeof(uint32));
	}

	sdata.lower_mod_prime = (uint32 *)malloc(sdata.pboundi * sizeof(uint32));
	allocated_bytes += sdata.pboundi * sizeof(uint32);
	if (sdata.lower_mod_prime == NULL)
	{
		printf("error allocating lower mod prime\n");
		exit(-1);
	}
	else
	{
		if (VFLAG > 2)
			printf("allocated %u bytes for lower mod prime\n",
			(uint32)sdata.pboundi * (uint32)sizeof(uint32));
	}

	if (sdata.pboundi > BUCKETSTARTI)
	{
		uint64 flagsperline = numlinebytes * 8;
		uint64 num_hits = 0;
		uint64 hits_per_bucket;

		//then we have primes bigger than BUCKETSTARTP
		bucket_depth = (uint32)(sdata.pboundi - BUCKETSTARTI);
		for (i=BUCKETSTARTI; i<sdata.pboundi; i++)
			num_hits += (flagsperline / sdata.sieve_p[i] + 1);

		//assume hits are evenly distributed among buckets.
		hits_per_bucket = num_hits / sdata.blocks;

		//add some margin for uneven bucket distribution
		hits_per_bucket = (uint64)((double)hits_per_bucket * 1.1);

		//set the bucket allocation amount
		bucket_alloc = hits_per_bucket;

		max_bucket_usage = 0;
	}
	else
	{
		bucket_depth = 0;
	}

	//uncomment this to disable bucket sieving
	//bucket_depth = 0;

	thread_data = (thread_soedata_t *)malloc(THREADS * sizeof(thread_soedata_t));
	allocated_bytes += THREADS * sizeof(thread_soedata_t);
	for (i=0; i<THREADS; i++)
	{
		thread_soedata_t *thread = thread_data + i;

		//allocate a bound for each block
		thread->ddata.pbounds = (uint64 *)malloc(
			(sdata.blocks + (sdata.partial_block_b > 0))*sizeof(uint64));
		allocated_bytes += (sdata.blocks + (sdata.partial_block_b > 0))*sizeof(uint64);

		thread->ddata.pbounds[0] = sdata.pboundi;

		//we'll need to store the offset into the next block for each prime.
		//actually only need those primes less than BUCKETSTARTP since bucket sieving
		//doesn't use the offset array.
		j = MIN(sp,BUCKETSTARTI);
		thread->ddata.offsets = (uint32 *)malloc(j * sizeof(uint32));
		allocated_bytes += j * sizeof(uint32);
		if (thread->ddata.offsets == NULL)
		{
			printf("error allocating offsets\n");
			exit(-1);
		}
		else
		{
			if (VFLAG > 2)
				printf("allocated %u bytes for offsets\n",j * sizeof(uint32));
		}

		//allocate the line for this thread
		thread->ddata.line = (uint8 *)malloc(numlinebytes * sizeof(uint8));
		allocated_bytes += numlinebytes * sizeof(uint8);
		if (thread->ddata.line == NULL)
		{
			printf("error allocating sieve line\n");
			exit(-1);
		}
		else
		{
			if (VFLAG > 2)
				printf("allocated %u bytes for sieve line\n",numlinebytes * sizeof(uint8));
		}

#ifdef DO_SPECIAL_COUNT
		//allocate an array so we can count intervals smaller than the sieve line
		j = (sdata.orig_hlimit - sdata.orig_llimit) / 1000000000;
		j += ((sdata.orig_hlimit - sdata.orig_llimit) % 1000000000 > 0);
		thread->ddata.num_special_bins = j;
		thread->ddata.special_count = (uint32 *)calloc(j, sizeof(uint32));
		sdata.num_special_bins = j;
		sdata.special_count = (uint32 *)calloc(j, sizeof(uint32));
		allocated_bytes += j * 2 * sizeof(uint32);
		if (thread->ddata.special_count == NULL)
		{
			printf("error allocating special line count\n");
			exit(-1);
		}
		else
		{
			if (VFLAG > 2)
				printf("allocated %u bytes for special line count\n",j * sizeof(uint32));
		}
#endif

		if (bucket_depth > BUCKET_BUFFER)
		{			
			//create a bucket for each block
			thread->ddata.sieve_buckets = (soe_bucket_t **)malloc(
				sdata.blocks * sizeof(soe_bucket_t *));
			allocated_bytes += sdata.blocks * sizeof(soe_bucket_t *);

			if (thread->ddata.sieve_buckets == NULL)
			{
				printf("error allocating buckets\n");
				exit(-1);
			}
			else
			{
				if (VFLAG > 2)
					printf("allocated %u bytes for bucket bases\n",sdata.blocks * sizeof(soe_bucket_t *));
			}

			//create a hit counter for each bucket
			thread->ddata.bucket_hits = (uint32 *)malloc(
				sdata.blocks * sizeof(uint32));
			allocated_bytes += sdata.blocks * sizeof(uint32);
			if (thread->ddata.bucket_hits == NULL)
			{
				printf("error allocating hit counters\n");
				exit(-1);
			}
			else
			{
				if (VFLAG > 2)
					printf("allocated %u bytes for hit counters\n",sdata.blocks * sizeof(uint32));
			}

			//each bucket must be able to hold a hit from every prime used above BUCKETSTARTP.
			//this is overkill, because every prime will not hit every bucket when
			//primes are greater than FLAGSIZE.  but we always write to the end of each
			//bucket so the depth probably doesn't matter much from a cache access standpoint,
			//just from a memory capacity standpoint, and we shouldn't be in any danger
			//of that as long as the input range is managed (should be <= 10B).
			thread->ddata.bucket_depth = bucket_depth;
			thread->ddata.bucket_alloc = bucket_alloc;

			for (j = 0; j < sdata.blocks; j++)
			{
				thread->ddata.sieve_buckets[j] = (soe_bucket_t *)malloc(
					bucket_alloc * sizeof(soe_bucket_t));
				allocated_bytes += bucket_alloc * sizeof(soe_bucket_t);

				if (thread->ddata.sieve_buckets[j] == NULL)
				{
					printf("error allocating buckets\n");
					exit(-1);
				}			
				
				thread->ddata.bucket_hits[j] = 0;
			}
			////thread->ddata.bucket_hits = 0;
			if (VFLAG > 2)
				printf("allocated %u bytes for buckets\n",sdata.blocks * bucket_alloc * sizeof(soe_bucket_t));
		}	
		else
			thread->ddata.bucket_depth = 0;

		//this threads' count of primes in its' line
		thread->linecount = 0;
		//share the common static data structure
		thread->sdata = sdata;
	}

	getRoots(&sdata);

	if (VFLAG > 2)
	{	
		
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		printf("sieving range %lu to %lu\n",lowlimit,*highlimit);
		printf("using %lu primes, max prime = %lu  \n",sdata.pboundi,sdata.pbound);
		printf("using %lu residue classes\n",numclasses);
		printf("lines have %lu bytes and %lu flags\n",numlinebytes,numlinebytes * 8);
		printf("lines broken into = %lu blocks of size %u\n",sdata.blocks,BLOCKSIZE);
		printf("blocks contain %u flags and cover %lu primes\n", FLAGSIZE, sdata.blk_r);
		if (bucket_depth > BUCKET_BUFFER)
		{
			printf("bucket sieving %u primes > %u\n",bucket_depth,BUCKETSTARTP);
			printf("allocating space for %lu hits per bucket\n",bucket_alloc);
		}
		printf("using %lu bytes for sieving storage\n",allocated_bytes);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		printf("sieving range %llu to %llu\n",lowlimit,*highlimit);
		printf("using %llu primes, max prime = %llu  \n",sdata.pboundi,sdata.pbound);
		printf("using %llu residue classes\n",numclasses);
		printf("lines have %llu bytes and %llu flags\n",numlinebytes,numlinebytes * 8);
		printf("lines broken into = %llu blocks of size %u\n",sdata.blocks,BLOCKSIZE);
		printf("blocks contain %u flags and cover %llu primes\n", FLAGSIZE, sdata.blk_r);
		if (bucket_depth > BUCKET_BUFFER)
		{
			printf("bucket sieving %u primes > %u\n",bucket_depth,BUCKETSTARTP);
			printf("allocating space for %llu hits per bucket\n",bucket_alloc);
		}
		printf("using %llu bytes for sieving storage\n",allocated_bytes);
#else
		printf("sieving range %I64u to %I64u\n",lowlimit,*highlimit);
		printf("using %I64u primes, max prime = %I64u  \n",sdata.pboundi,sdata.pbound);
		printf("using %I64u residue classes\n",numclasses);
		printf("lines have %I64u bytes and %I64u flags\n",numlinebytes,numlinebytes * 8);
		printf("lines broken into = %I64u blocks of size %u\n",sdata.blocks,BLOCKSIZE);
		printf("blocks contain %u flags and cover %I64u primes\n", FLAGSIZE, sdata.blk_r);
		if (bucket_depth > BUCKET_BUFFER)
		{
			printf("bucket sieving %u primes > %u\n",bucket_depth,BUCKETSTARTP);
			printf("allocating space for %I64u hits per bucket\n",bucket_alloc);
		}
		printf("using %I64u bytes for sieving storage\n",allocated_bytes);
#endif
			
	}

	/* activate the threads one at a time. The last is the
	   master thread (i.e. not a thread at all). */

	for (i = 0; i < THREADS - 1; i++)
		start_soe_worker_thread(thread_data + i, 0);

	start_soe_worker_thread(thread_data + i, 1);

	//main sieve, line by line
	k = 0;	//count total lines processed
	pchar = 0;
	num_p = 0;
	
	while (k < numclasses)
	{
		fflush(stdout);
		//assign lines to the threads
		j = 0;	//how many lines to process this pass
		for (i = 0; ((int)i < THREADS) && (k < numclasses); i++, k++)
		{
			thread_soedata_t *t = thread_data + i;

			t->current_line = (uint32)k;
			j++;
		}

		//process the lines
		for (i = 0; i < j; i++) 
		{
			thread_soedata_t *t = thread_data + i;

			if (i == j - 1) {
				
				if (count)
				{
					sieve_line(t);
#ifdef DO_SPECIAL_COUNT
					count_line_special(t);
#else
					count_line(t);
#endif
				}
				else
				{
					sieve_line(t);
					primes_from_lineflags(t);
				}
			}
			else {
				if (count)
					t->command = SOE_COMMAND_SIEVE_AND_COUNT;
				else
					t->command = SOE_COMMAND_SIEVE_AND_COMPUTE;
#if defined(WIN32) || defined(_WIN64)
				SetEvent(t->run_event);
#else
				pthread_cond_signal(&t->run_cond);
				pthread_mutex_unlock(&t->run_lock);
#endif
			}
		}

		//wait for each thread to finish
		for (i = 0; i < j; i++) {
			thread_soedata_t *t = thread_data + i;

			if (i < j - 1) {
#if defined(WIN32) || defined(_WIN64)
				WaitForSingleObject(t->finish_event, INFINITE);
#else
				pthread_mutex_lock(&t->run_lock);
				while (t->command != SOE_COMMAND_WAIT)
					pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
			}
		}

		//printf a progress report if counting
		if (count && VFLAG >= 0)
		{
			//don't print status if computing primes, because lots of routines within
			//yafu do this and they don't want this side effect
			for (i = 0; i<pchar; i++)
				printf("\b");
			pchar = printf("%d%%",(int)((double)k / (double)(numclasses) * 100.0));
			fflush(stdout);
		}

		
		if (count)
		{
			//accumulate the results from each line
#ifdef DO_SPECIAL_COUNT
			for (i = 0; i < j; i++)
			{
				int ix;
				num_p += thread_data[i].linecount;
				
				for (ix=0; ix < sdata.num_special_bins; ix++)
					sdata.special_count[ix] += thread_data[i].ddata.special_count[ix];
			}

#else
			for (i = 0; i < j; i++)
				num_p += thread_data[i].linecount;
#endif
		}
		else
		{
			//accumulate the primes from each line
			for (i=0; i< j; i++)
			{
				uint64 start_index = num_p;
				uint64 linecount = thread_data[i].linecount;
				uint64 index;
				
				
				if (num_p == 0)
				{
					//add the first primes to the mergelist
					for (index = 0; index < linecount; index++)
					{
						mergeprimes[index] = thread_data[i].ddata.primes[index];
						locprimes[index] = thread_data[i].ddata.primes[index];
					}
				}
				else
				{
					//merge this lines primes (linecount), with the previous merge (num_p)
					//this lines primes are ordered, and so is the previous merged list.

					i1 = i2 = i3 = 0;
					while (i1 < num_p && i2 < linecount) {
						if (locprimes[i1] < thread_data[i].ddata.primes[i2])
							mergeprimes[i3++] = locprimes[i1++];
						else
							mergeprimes[i3++] = thread_data[i].ddata.primes[i2++];
					}
					while (i1 < num_p)
						mergeprimes[i3++] = locprimes[i1++];
					while (i2 < linecount)
						mergeprimes[i3++] = thread_data[i].ddata.primes[i2++];

					for (i1=0; i1<(num_p + linecount); i1++)
						locprimes[i1] = mergeprimes[i1];

				}

				num_p += linecount;
				//then free the line's primes
				free(thread_data[i].ddata.primes);

			}
		}
	}

	//stop the worker threads and free stuff not needed anymore
	for (i=0; i<THREADS - 1; i++)
	{
		stop_soe_worker_thread(thread_data + i, 0);
		free(thread_data[i].ddata.offsets);
	}
	stop_soe_worker_thread(thread_data + i, 1);
	free(thread_data[i].ddata.offsets);

	if (count && VFLAG >= 0)
	{
		//don't print status if computing primes, because lots of routines within
		//yafu do this and they don't want this side effect
		for (i = 0; i<pchar; i++)
			printf("\b");
	}

	if (count)
	{
		//add in relevant sieving primes not captured in the flag arrays
		if (sdata.pbound > lowlimit)
		{
			i=0;
			while (i<thread_data[0].ddata.pbounds[0])
			{ 
				if (sdata.sieve_p[i] > lowlimit)
				{
					num_p++;
#ifdef DO_SPECIAL_COUNT
					sdata.special_count[0]++;
#endif
				}
				i++;
			}
		}
	}
	else
	{
		//the sieve primes are not in the line array, so they must be added
		//in if necessary
		uint64 start_index = num_p;
		it=0;
		
		//first count them
		if (sdata.pbound > lowlimit)
		{
			i=0;
			while (i<thread_data[0].ddata.pbounds[0])
			{ 
				if (sdata.sieve_p[i] > lowlimit)
					it++;
				i++;
			}
		}

		//merge the sieve primes (it), with the previous merge (num_p)
		//this lines primes are ordered, and so is the previous merged list.
		i1 = i2 = i3 = 0;
		while (i1 < num_p && i2 < it) {
			if (locprimes[i1] < sdata.sieve_p[i2])
				mergeprimes[i3++] = locprimes[i1++];
			else
				mergeprimes[i3++] = sdata.sieve_p[i2++];
		}
		while (i1 < num_p)
			mergeprimes[i3++] = locprimes[i1++];
		while (i2 < it)
			mergeprimes[i3++] = sdata.sieve_p[i2++];

		num_p += it;

		//copy the local mergeprimes into the output primes array
		memcpy(primes,mergeprimes,num_p * sizeof(uint64));

	}

	//printf("maximum bucket usage = %u\n",max_bucket_usage);
#ifdef DO_SPECIAL_COUNT
	//print the special counts
	for (i=0; i<sdata.num_special_bins; i++)
		printf("count in range %d = %u\n",i,sdata.special_count[i]);
#endif

	for (i=0; i<THREADS; i++)
	{
		thread_soedata_t *thread = thread_data + i;
		free(thread->ddata.pbounds);
		free(thread->ddata.line);
#ifdef DO_SPECIAL_COUNT
		free(thread->ddata.special_count);
#endif
	}

	if (bucket_depth > BUCKET_BUFFER)
	{
		for (i=0; i< THREADS; i++)
		{
			thread_soedata_t *thread = thread_data + i;

			free(thread->ddata.bucket_hits);
			for (j=0; j < thread->sdata.blocks; j++)
				free(thread->ddata.sieve_buckets[j]);
			free(thread->ddata.sieve_buckets);
		}
	}

	free(locprimes);
	free(mergeprimes);
	free(sdata.root);
	free(sdata.lower_mod_prime);
#ifdef DO_SPECIAL_COUNT
	free(sdata.special_count);
#endif
	free(thread_data);
	free(sdata.sieve_p);
	free(sdata.rclass);

	return num_p;
}

void primes_from_lineflags(thread_soedata_t *t)
{
	//the sieve primes are not in the line array, so they must be added
	//in if necessary
	soe_staticdata_t *sdata = &t->sdata;
	soe_dynamicdata_t *ddata = &t->ddata;
	uint8 *line = t->ddata.line;
	uint64 current_line = t->current_line;
	uint64 prime, num_alloc;
	uint64 i,j,it;

	it=0;

	count_line(t);
	num_alloc = t->linecount;
	ddata->primes = (uint64 *)malloc(num_alloc * sizeof(uint64));
	if (ddata->primes == NULL)
	{
		printf("failed to allocate primes array in primes_from_lineflags\n");
		exit(-1);
	}

	//this will find all the primes in the line.  when other lines are appended, the primes
	//will be out of order and will need to be sorted.
	for (i=0;i<sdata->numlinebytes;i++)
	{
		for (j=0;j<BITSINBYTE;j++)
		{
			if (line[i] & nmasks[j])
			{
				prime = sdata->prodN * (i*BITSINBYTE + j) + 
					sdata->rclass[current_line] + sdata->lowlimit;

				//only store the prime if it is within our requested bounds
				if ((prime >= sdata->orig_llimit) &&  (prime <= sdata->orig_hlimit))
				{
					ddata->primes[it] = prime;
					it++;
				}
			}
			
		}
	}

	if (it != t->linecount)
		printf("warning, counts do not match after computing primes\n");

	t->linecount = it;

	return;
}

void count_line(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint8 *line = thread_data->ddata.line;
	uint32 current_line = thread_data->current_line;
	uint64 numlinebytes = sdata->numlinebytes;
	uint64 lowlimit = sdata->lowlimit;
	uint64 prodN = sdata->prodN;
	uint64 *flagblock64 = (uint64 *)line;
	uint8 *flagblock;
	uint64 i, k, it = 0;
	int ix;

	//zero out any bits below the requested range
	for (i=lowlimit + sdata->rclass[current_line], ix=0; i < sdata->orig_llimit; i += prodN, ix++)
		line[ix >> 3] &= masks[ix & 7];
	
	//and any high bits above the requested range
	for (i=sdata->highlimit + sdata->rclass[current_line] - prodN, ix=0; i > sdata->orig_hlimit; i -= prodN, ix++)
		line[numlinebytes - 1 - (ix >> 3)] &= masks[7 - (ix & 7)];

	//then just count em
	for (i=0;i<(numlinebytes >> 3);i++)
	{
		/* Convert to 64-bit unsigned integer */    
		uint64 x = flagblock64[i];
	    
		/*  Employ bit population counter algorithm from Henry S. Warren's
		 *  "Hacker's Delight" book, chapter 5.   Added one more shift-n-add
		 *  to accomdate 64 bit values.
		 */
		x = x - ((x >> 1) & 0x5555555555555555);
		x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
		x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F;
		x = x + (x >> 8);
		x = x + (x >> 16);
		x = x + (x >> 32);

		it += (x & 0x000000000000003F);
	}

	//potentially misses the last few bytes
	//use the simpler baseline method to get these few
	flagblock = line;
	for (k=0; k<(numlinebytes & 7);k++)
	{
		it += (flagblock[(i<<3)+k] & nmasks[0]) >> 7;
		it += (flagblock[(i<<3)+k] & nmasks[1]) >> 6;
		it += (flagblock[(i<<3)+k] & nmasks[2]) >> 5;
		it += (flagblock[(i<<3)+k] & nmasks[3]) >> 4;
		it += (flagblock[(i<<3)+k] & nmasks[4]) >> 3;
		it += (flagblock[(i<<3)+k] & nmasks[5]) >> 2;
		it += (flagblock[(i<<3)+k] & nmasks[6]) >> 1;
		it += (flagblock[(i<<3)+k] & nmasks[7]);
	}

	thread_data->linecount = it;

	return;
}

void count_line_special(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint8 *line = thread_data->ddata.line;
	uint32 current_line = thread_data->current_line;
	uint64 numlinebytes = sdata->numlinebytes;
	uint64 lowlimit = sdata->lowlimit;
	uint64 prodN = sdata->prodN;
	uint64 *flagblock64 = (uint64 *)line;
	uint64 i, k, it, lower, upper;
	int ix;
	int64 start, stop;

	//zero out any bits below the requested range
	for (i=lowlimit + sdata->rclass[current_line], ix=0; i < sdata->orig_llimit; i += prodN, ix++)
		line[ix >> 3] &= masks[ix & 7];
	
	//and any high bits above the requested range
	for (i=sdata->highlimit + sdata->rclass[current_line] - prodN, ix=0; i > sdata->orig_hlimit; i -= prodN, ix++)
		line[numlinebytes - 1 - (ix >> 3)] &= masks[7 - (ix & 7)];

	//count each block of 1e9
	lower = sdata->orig_llimit;
	upper = lower;
	k = 0;
	thread_data->linecount = 0;
	while (upper != sdata->orig_hlimit)
	{
		//set the bounds for the next batch
		upper = upper + 1000000000; 
		if (upper > sdata->orig_hlimit)
			upper = sdata->orig_hlimit;

		//find the starting byte number.  first find the number of bits between the current lower
		//limit and the start of the line.
		start = (int64)((lower - lowlimit) / prodN);

		//we'll be counting in 64 bit chunks, so compute how many 64 bit chunks this is
		start /= 64;
		
		//start a little before the range, to account for rounding errors
		start -= 2;

		if (start < 0) start = 0;

		//find the stopping byte number: first find the number of bits between the current upper
		//limit and the start of the line.
		stop = (int64)((upper - lowlimit) / prodN);

		//we'll be counting in 64 bit chunks, so compute how many 64 bit chunks this is
		stop /= 64;
		
		//stop a little after the range, to account for rounding errors
		stop += 2;

		if (stop > (numlinebytes >> 3)) stop = (numlinebytes >> 3);

		//count these bytes
		it = 0;
		for (ix = start; ix < stop; ix++)
		{
			/* Convert to 64-bit unsigned integer */    
			uint64 x = flagblock64[ix];
		    
			/*  Employ bit population counter algorithm from Henry S. Warren's
			 *  "Hacker's Delight" book, chapter 5.   Added one more shift-n-add
			 *  to accomdate 64 bit values.
			 */
			x = x - ((x >> 1) & 0x5555555555555555);
			x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
			x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F;
			x = x + (x >> 8);
			x = x + (x >> 16);
			x = x + (x >> 32);

			it += (x & 0x000000000000003F);
		}

		//then correct the counts
		//zero out any bits below the requested range
		for (i= (start * 64) * prodN + sdata->rclass[current_line] + lowlimit, ix=0; i < lower; i += prodN, ix++)
		{
			if (line[(ix >> 3) + (start << 3)] & nmasks[ix & 7])
				it--;
		}
		
		//and any high bits above the requested range
		for (i=(stop * 64) * prodN + sdata->rclass[current_line] - prodN + lowlimit, ix=0; i > upper; i -= prodN, ix++)
		{
			if (line[(stop << 3) - 1 - (ix >> 3)] & nmasks[7 - (ix & 7)])
				it--;

		}
		
		//add the count to the special array
		thread_data->ddata.special_count[k] = it;
		thread_data->linecount += it;
		k++;
		lower = upper;
	}

}

void sieve_line(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint8 *line = thread_data->ddata.line;
	uint32 current_line = thread_data->current_line;
	
	//stuff for bucket sieving
	soe_bucket_t *bptr;
	soe_bucket_t **buckets;
	uint32 *nptr;
	uint32 linesize = FLAGSIZE * sdata->blocks, bnum;
	//uint32 *bptr;

	uint8 *flagblock;
	uint64 startprime = sdata->startprime;
	uint64 i,j,k;
	uint32 prime;
	uint32 maxP;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (uint32)(sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	//for the current line, find the offsets past the low limit
	get_offsets(thread_data);

#define SOE_UNROLL_LOOPS
#ifdef SOE_UNROLL_LOOPS

	flagblock = line;
	for (i=0;i<sdata->blocks;i++)
	{
		uint64 *flagblock64;
		int mask_step, mask_step2;
		int mask_num, mask_num2;

		//sieve the block with each effective prime
		//set all flags for this block, which also puts it into cache for the sieving
		//to follow
		memset(flagblock,255,BLOCKSIZE);

		flagblock64 = (uint64 *)flagblock;		
		
		//do the smallest primes in predetermined 64 bit batches
		if (startprime == 2)
		{
			for (k=0, mask_step = _64_MOD_P[0], mask_num = ddata->offsets[2],
				mask_step2 = _64_MOD_P[1], mask_num2 = ddata->offsets[3]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_5_MASKS[mask_num] & _7_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 5 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 7 + mask_num2;
			}
			ddata->offsets[2]= (uint32)mask_num;
			ddata->offsets[3]= (uint32)mask_num2;

			for (k=0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
				mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 11 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
			}
			ddata->offsets[4]= (uint32)mask_num;
			ddata->offsets[5]= (uint32)mask_num2;
		}
		else if (startprime == 3)
		{
			for (k=0, mask_step = _64_MOD_P[1], mask_num = ddata->offsets[3]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= _7_MASKS[mask_num];
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 7 + mask_num;
			}
			ddata->offsets[3]= (uint32)mask_num;

			for (k=0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
				mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 11 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
			}
			ddata->offsets[4]= (uint32)mask_num;
			ddata->offsets[5]= (uint32)mask_num2;
		}
		else if (startprime == 4)
		{

			for (k=0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
				mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 11 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
			}
			ddata->offsets[4]= (uint32)mask_num;
			ddata->offsets[5]= (uint32)mask_num2;
		}
		else if (startprime == 5)
		{
			for (k=0, mask_step = _64_MOD_P[3], mask_num = ddata->offsets[5];
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= _13_MASKS[mask_num];
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 13 + mask_num;
			}
			ddata->offsets[5]= (uint32)mask_num;
		}	

		for (k=0, mask_step = _64_MOD_P[4], mask_num = ddata->offsets[6],
			mask_step2 = _64_MOD_P[5], mask_num2 = ddata->offsets[7]; 
			k<FLAGSIZE >> 6; k++)
		{
			flagblock64[k] &= (_17_MASKS[mask_num] & _19_MASKS[mask_num2]);
			mask_num -= mask_step;
			if (mask_num < 0) mask_num = 17 + mask_num;
			mask_num2 -= mask_step2;
			if (mask_num2 < 0) mask_num2 = 19 + mask_num2;
		}
		ddata->offsets[6]= (uint32)mask_num;
		ddata->offsets[7]= (uint32)mask_num2;

		for (k=0, mask_step = _64_MOD_P[6], mask_num = ddata->offsets[8],
			mask_step2 = _64_MOD_P[7], mask_num2 = ddata->offsets[9]; 
			k<FLAGSIZE >> 6; k++)
		{
			flagblock64[k] &= (_23_MASKS[mask_num] & _29_MASKS[mask_num2]);
			mask_num -= mask_step;
			if (mask_num < 0) mask_num = 23 + mask_num;
			mask_num2 -= mask_step2;
			if (mask_num2 < 0) mask_num2 = 29 + mask_num2;
		}
		ddata->offsets[8]= (uint32)mask_num;
		ddata->offsets[9]= (uint32)mask_num2;	
		
		//one is not a prime
		if ((sdata->rclass[current_line] == 1) &&
			(sdata->lowlimit <= 1) && (i == 0))
			flagblock[0] &= 0xfe;

		//unroll the loop: all primes less than this max hit the interval at least 16 times
		maxP = FLAGSIZE >> 4;

		for (j=10;j<ddata->pbounds[i];j++)
		{
			uint32 tmpP;
			uint64 stop;
			uint64 p1,p2,p3;

			prime = sdata->sieve_p[j];
			if (prime > maxP)
				break;

			tmpP = prime << 4;
			stop = FLAGSIZE - tmpP + prime;
			k=ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
				flagblock[k>>3] &= masks[k&7];
				flagblock[(k+p1)>>3] &= masks[(k+p1)&7];
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];
				flagblock[(k+p3)>>3] &= masks[(k+p3)&7];
				k += (prime << 2);
			}

			for (;k<FLAGSIZE;k+=prime)
				flagblock[k>>3] &= masks[k&7];


			
			//if ((j >= 2) && (j <= 10))
			//{
			//	printf("actual = %lx; next offset = %u\n",flagblock64[0], k - FLAGSIZE);
			//}
			ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			
		}

		//unroll the loop: all primes less than this max hit the interval at least 8 times
		maxP = FLAGSIZE >> 3;

		for (;j<ddata->pbounds[i];j++)
		{
			uint32 tmpP;
			uint64 stop;
			uint64 p2, p4;

			prime = sdata->sieve_p[j];
			if (prime > maxP)
				break;

			tmpP = prime << 3;
			stop = FLAGSIZE - tmpP + prime;
			k=ddata->offsets[j];
			p2 = prime<<1;
			p4 = prime<<2;

			while (k < stop)
			{
				flagblock[k>>3] &= masks[k&7];								//0 * prime
				flagblock[(k+prime)>>3] &= masks[(k+prime)&7];				//1 * prime
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];					//2 * prime
				flagblock[(k+prime+p2)>>3] &= masks[(k+prime+p2)&7];		//3 * prime
				flagblock[(k+p4)>>3] &= masks[(k+p4)&7];					//4 * prime
				flagblock[(k+prime+p4)>>3] &= masks[(k+prime+p4)&7];		//5 * prime
				flagblock[(k+p2+p4)>>3] &= masks[(k+p2+p4)&7];				//6 * prime
				flagblock[(k+prime+p2+p4)>>3] &= masks[(k+prime+p2+p4)&7];	//7 * prime
				k += (prime << 3);											//advance
			}

			for (;k<FLAGSIZE;k+=prime)								//finish
				flagblock[k>>3] &= masks[k&7];

			ddata->offsets[j]= (uint32)(k - FLAGSIZE);
		}

		//unroll the loop: all primes less than this max hit the interval at least 4 times
		maxP = FLAGSIZE >> 2;

		for (;j<ddata->pbounds[i];j++)
		{
			uint32 tmpP;
			uint64 stop;
			uint64 p2;

			prime = sdata->sieve_p[j];
			if (prime > maxP)
				break;

			tmpP = prime << 2;
			stop = FLAGSIZE - tmpP + prime;
			k=ddata->offsets[j];
			p2 = prime<<1;
			while (k < stop)
			{
				flagblock[k>>3] &= masks[k&7];								//0 * prime
				flagblock[(k+prime)>>3] &= masks[(k+prime)&7];				//1 * prime
				flagblock[(k+p2)>>3] &= masks[(k+p2)&7];					//2 * prime
				flagblock[(k+prime+p2)>>3] &= masks[(k+prime+p2)&7];		//3 * prime
				k += (prime << 2);											//advance
			}

			for (;k<FLAGSIZE;k+=prime)								//finish
				flagblock[k>>3] &= masks[k&7];

			ddata->offsets[j]= (uint32)(k - FLAGSIZE);
		}

		
		if (ddata->bucket_depth > BUCKET_BUFFER)
		{
			for (;j<ddata->pbounds[i];j++)
			{
				prime = sdata->sieve_p[j];
				if (prime > BUCKETSTARTP)
					break;
				for (k=ddata->offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			}

			//finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];	
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			if (nptr[i] > max_bucket_usage)
				max_bucket_usage = nptr[i];

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j=0; j < (nptr[i] & (uint32)(~7)); j+=8)
			{				
				//unload 8 hits
				flagblock[(bptr[j + 0].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 0].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 1].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 1].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 2].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 2].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 3].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 3].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 4].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 4].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 5].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 5].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 6].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 6].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 7].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 7].root & FLAGSIZEm1) & 7];

				
				//then compute their next hit
				bptr[j + 0].root += bptr[j + 0].prime;	
				bptr[j + 1].root += bptr[j + 1].prime;		
				bptr[j + 2].root += bptr[j + 2].prime;	

				if (bptr[j + 0].root < linesize)			
				{	
					bnum = (bptr[j + 0].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 0].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 0].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 1].root < linesize)			
				{	
					bnum = (bptr[j + 1].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 1].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 1].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 2].root < linesize)			
				{	
					bnum = (bptr[j + 2].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 2].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 2].prime;
					nptr[bnum]++;						
				}	

				bptr[j + 3].root += bptr[j + 3].prime;		
				bptr[j + 4].root += bptr[j + 4].prime;	
				bptr[j + 5].root += bptr[j + 5].prime;		
				if (bptr[j + 3].root < linesize)			
				{	
					bnum = (bptr[j + 3].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 3].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 3].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 4].root < linesize)			
				{	
					bnum = (bptr[j + 4].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 4].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 4].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 5].root < linesize)			
				{	
					bnum = (bptr[j + 5].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 5].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 5].prime;
					nptr[bnum]++;						
				}	

				bptr[j + 6].root += bptr[j + 6].prime;		
				bptr[j + 7].root += bptr[j + 7].prime;		
				if (bptr[j + 6].root < linesize)			
				{	
					bnum = (bptr[j + 6].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 6].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 6].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 7].root < linesize)			
				{	
					bnum = (bptr[j + 7].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 7].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 7].prime;
					nptr[bnum]++;						
				}	
				
			}

			for (;j < nptr[i]; j++)
			{
				flagblock[(bptr[j].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j].root & FLAGSIZEm1) & 7];

				
				bptr[j].root += bptr[j].prime;		
				if (bptr[j].root < linesize)			
				{	
					bnum = (bptr[j].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j].prime;
					nptr[bnum]++;						
				}
				
			}

		}
		else
		{
			//finish with primes greater than (flagblocklimit >> 2)
			for (;j<ddata->pbounds[i];j++)
			{
				prime = sdata->sieve_p[j];
				for (k=ddata->offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			}
		}

		flagblock += BLOCKSIZE;
	}
	

#else

	flagblock = line;
	if (ddata->bucket_depth > BUCKET_BUFFER)
	{

		for (i=0;i<sdata->blocks;i++)
		{
			//set all flags for this block, which also puts it into cache for the sieving
			//to follow
			memset(flagblock,255,BLOCKSIZE);

			//one is not a prime
			if ((sdata->rclass[current_line] == 1) &&
				(sdata->lowlimit <= 1) && (i == 0))
				flagblock[0] = 0xfe;

		
			for (j = startprime; j<ddata->pbounds[i];j++)
			{
				prime = sdata->sieve_p[j];
				if (prime > BUCKETSTARTP)
					break;
				for (k=ddata->offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			}

			//finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			for (j=0; j < (ddata->bucket_hits[i] & (uint32)(~15)); j+=16)
			{				
				//unload 8 hits
				flagblock[(bptr[j + 0].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 0].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 1].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 1].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 2].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 2].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 3].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 3].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 4].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 4].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 5].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 5].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 6].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 6].root & FLAGSIZEm1) & 7];
				flagblock[(bptr[j + 7].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 7].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 8].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 8].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 9].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 9].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 10].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 10].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 11].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 11].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 12].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 12].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 13].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 13].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 14].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 14].root & FLAGSIZEm1) & 7];
				//flagblock[(bptr[j + 15].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 15].root & FLAGSIZEm1) & 7];

				//then compute their next hit
				bptr[j + 0].root += bptr[j + 0].prime;	
				bptr[j + 1].root += bptr[j + 1].prime;		
				bptr[j + 2].root += bptr[j + 2].prime;	

				if (bptr[j + 0].root < linesize)			
				{	
					bnum = (bptr[j + 0].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 0].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 0].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 1].root < linesize)			
				{	
					bnum = (bptr[j + 1].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 1].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 1].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 2].root < linesize)			
				{	
					bnum = (bptr[j + 2].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 2].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 2].prime;
					nptr[bnum]++;						
				}	

				bptr[j + 3].root += bptr[j + 3].prime;		
				bptr[j + 4].root += bptr[j + 4].prime;	
				bptr[j + 5].root += bptr[j + 5].prime;		
				if (bptr[j + 3].root < linesize)			
				{	
					bnum = (bptr[j + 3].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 3].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 3].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 4].root < linesize)			
				{	
					bnum = (bptr[j + 4].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 4].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 4].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 5].root < linesize)			
				{	
					bnum = (bptr[j + 5].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 5].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 5].prime;
					nptr[bnum]++;						
				}	

				bptr[j + 6].root += bptr[j + 6].prime;		
				bptr[j + 7].root += bptr[j + 7].prime;		
				if (bptr[j + 6].root < linesize)			
				{	
					bnum = (bptr[j + 6].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 6].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 6].prime;
					nptr[bnum]++;						
				}	

				if (bptr[j + 7].root < linesize)			
				{	
					bnum = (bptr[j + 7].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j + 7].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j + 7].prime;
					nptr[bnum]++;						
				}	

			}

			for (;j < ddata->bucket_hits[i]; j++)
			{
				flagblock[(bptr[j].root & FLAGSIZEm1) >> 3] &= masks[(bptr[j].root & FLAGSIZEm1) & 7];
				bptr[j].root += bptr[j].prime;		
				if (bptr[j].root < linesize)			
				{	
					bnum = (bptr[j].root >> FLAGBITS);
					buckets[bnum][nptr[bnum]].root = bptr[j].root;
					buckets[bnum][nptr[bnum]].prime = bptr[j].prime;
					nptr[bnum]++;						
				}	
			}

			flagblock += BLOCKSIZE;

		}
	}
	else
	{
		for (i=0;i<sdata->blocks;i++)
		{
			uint64 *flagblock64;
			int mask_step, mask_step2;
			int mask_num, mask_num2;

			//sieve the block with each effective prime
			//set all flags for this block, which also puts it into cache for the sieving
			//to follow
			memset(flagblock,255,BLOCKSIZE);

			flagblock64 = (uint64 *)flagblock;		
			
			//do the smallest primes in predetermined 64 bit batches
			if (startprime == 2)
			{
				for (k=0, mask_step = _64_MOD_P[0], mask_num = ddata->offsets[2]; 
					k<FLAGSIZE >> 6; k++)
				{
					flagblock64[k] &= _5_MASKS[mask_num];
					mask_num -= mask_step;
					if (mask_num < 0) mask_num = 5 + mask_num;
				}
				ddata->offsets[2]= (uint32)mask_num;

				for (k=0, mask_step = _64_MOD_P[1], mask_num = ddata->offsets[3]; 
					k<FLAGSIZE >> 6; k++)
				{
					flagblock64[k] &= _7_MASKS[mask_num];
					mask_num -= mask_step;
					if (mask_num < 0) mask_num = 7 + mask_num;
				}
				ddata->offsets[3]= (uint32)mask_num;

			}
			else if (startprime == 3)
			{
				for (k=0, mask_step = _64_MOD_P[1], mask_num = ddata->offsets[3]; 
					k<FLAGSIZE >> 6; k++)
				{
					flagblock64[k] &= _7_MASKS[mask_num];
					mask_num -= mask_step;
					if (mask_num < 0) mask_num = 7 + mask_num;
				}
				ddata->offsets[3]= (uint32)mask_num;
			}

			for (k=0, mask_step = _64_MOD_P[2], mask_num = ddata->offsets[4],
				mask_step2 = _64_MOD_P[3], mask_num2 = ddata->offsets[5]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_11_MASKS[mask_num] & _13_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 11 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 13 + mask_num2;
			}
			ddata->offsets[4]= (uint32)mask_num;
			ddata->offsets[5]= (uint32)mask_num2;

			for (k=0, mask_step = _64_MOD_P[4], mask_num = ddata->offsets[6],
				mask_step2 = _64_MOD_P[5], mask_num2 = ddata->offsets[7]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_17_MASKS[mask_num] & _19_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 17 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 19 + mask_num2;
			}
			ddata->offsets[6]= (uint32)mask_num;
			ddata->offsets[7]= (uint32)mask_num2;

			for (k=0, mask_step = _64_MOD_P[6], mask_num = ddata->offsets[8],
				mask_step2 = _64_MOD_P[7], mask_num2 = ddata->offsets[9]; 
				k<FLAGSIZE >> 6; k++)
			{
				flagblock64[k] &= (_23_MASKS[mask_num] & _29_MASKS[mask_num2]);
				mask_num -= mask_step;
				if (mask_num < 0) mask_num = 23 + mask_num;
				mask_num2 -= mask_step2;
				if (mask_num2 < 0) mask_num2 = 29 + mask_num2;
			}
			ddata->offsets[8]= (uint32)mask_num;
			ddata->offsets[9]= (uint32)mask_num2;	


			//set all flags for this block, which also puts it into cache for the sieving
			//to follow
			//memset(flagblock,255,BLOCKSIZE);

			//one is not a prime
			if ((sdata->rclass[current_line] == 1) &&
				(sdata->lowlimit <= 1) && (i == 0))
				flagblock[0] = 0xfe;

			for (j = 10;j<ddata->pbounds[i];j++)
			{
				prime = sdata->sieve_p[j];
				for (k=ddata->offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			}
			flagblock += BLOCKSIZE;
		}
	}


#endif

	return;
}

void getRoots(soe_staticdata_t *sdata)
{
	int prime, prodN;
	uint64 startprime;
	uint64 i;

	prodN = (int)sdata->prodN;
	startprime = sdata->startprime;

	for (i=startprime;i<sdata->pboundi;i++)
	{
		uint32 inv;
		prime = sdata->sieve_p[i];

		inv = modinv_1(prodN,prime);
		
		sdata->root[i] = prime - inv;
		sdata->lower_mod_prime[i] = (sdata->lowlimit + 1) % prime;
	}

	return;
}

void get_offsets(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;

	uint64 i,startprime = sdata->startprime, prodN = sdata->prodN, block=0;
	uint32 prime, root, bnum;
	uint64 tmp2;
	uint32 diff = sdata->rclass[thread_data->current_line] - 1;
	int s;

	//failsafe: set all blocks to sieve with all primes.  the loop below will overwrite
	//these with better limits according to the size of flags in the blocks.
	for (i=0; i<sdata->blocks; i++)
	{
		ddata->pbounds[i] = sdata->pboundi;
		//initialize bucket
		if (ddata->bucket_depth > BUCKET_BUFFER)
			ddata->bucket_hits[i] = 0;
	}

	for (i=startprime;i<sdata->pboundi;i++)
	{
		prime = sdata->sieve_p[i];
		if ((prime > BUCKETSTARTP) && (ddata->bucket_depth > BUCKET_BUFFER))
			break;

		//find the first multiple of the prime which is greater than 'block1' and equal
		//to the residue class mod 'prodN'.  

		//if the prime is greater than the limit at which it is necessary to sieve
		//a block, start that prime in the next block.
		if (sdata->sieve_p[i] > ddata->blk_b_sqrt)
		{
			ddata->pbounds[block] = i;
			block++;
			ddata->lblk_b = ddata->ublk_b + prodN;
			ddata->ublk_b += sdata->blk_r;
			ddata->blk_b_sqrt = (uint64)(sqrt((int64)(ddata->ublk_b + prodN))) + 1;
		}

		//solving the congruence: rclass[current_line] == kp mod prodN for k
		//eGCD gives r and s such that r*p + s*prodN = gcd(p,prodN).
		//then k = r*class/gcd(p,prodN) is a solution.
		//the gcd of p and prodN is always 1 by construction of prodN and choice of p.  
		//therefore k = r * class is a solution.  furthermore, since the gcd is 1, there
		//is only one solution.
		//xGCD_1((int)prime,(int)prodN,&r,&s,&tmp);
		s = sdata->root[i];

		//the lower block bound (lblk_b) times s can exceed 64 bits for large ranges,
		//so reduce mod p here as well as when finding the root.
		tmp2 =  (uint64)s * (ddata->lblk_b % (uint64)prime);

		ddata->offsets[i] = (uint32)(tmp2 % (uint64)prime);
	}

	if (ddata->bucket_depth > BUCKET_BUFFER)
	{
		soe_bucket_t **bptr;
		uint32 *nptr;
		uint32 linesize = FLAGSIZE * sdata->blocks;
		
		for (; i<sdata->pboundi-1; i+=2)
		{
			uint64 tmp3;
			uint32 p2, r2;
			int s2;
						
			prime = sdata->sieve_p[i];
			p2 = sdata->sieve_p[i+1];

			s = sdata->root[i];
			s2 = sdata->root[i+1];
			
			//we solved for lower_mod_prime while computing the modular inverse of
			//each prime, for the residue class 1.  add the difference between this
			//residue class and 1 before multiplying by the modular inverse.
			tmp2 = (uint64)s * (uint64)(sdata->lower_mod_prime[i] + diff);			
			tmp3 = (uint64)s2 * (uint64)(sdata->lower_mod_prime[i+1] + diff);

			root = (uint32)(tmp2 % (uint64)prime);
			r2 = (uint32)(tmp3 % (uint64)p2);

			nptr = ddata->bucket_hits;
			bptr = ddata->sieve_buckets;

			if (root < linesize)			
			{	
				bnum = (root >> FLAGBITS);
				bptr[bnum][nptr[bnum]].root = root;
				bptr[bnum][nptr[bnum]].prime = prime;
				nptr[bnum]++;	
			}	

			if (r2 < linesize)			
			{	
				bnum = (r2 >> FLAGBITS);
				bptr[bnum][nptr[bnum]].root = r2;
				bptr[bnum][nptr[bnum]].prime = p2;
				nptr[bnum]++;	
			}	
			
		}

		if (i<sdata->pboundi)
		{
			prime = sdata->sieve_p[i];

			s = sdata->root[i];
			
			tmp2 = (uint64)s * (uint64)(sdata->lower_mod_prime[i] + diff);

			root = (uint32)(tmp2 % (uint64)prime);

			nptr = ddata->bucket_hits;
			bptr = ddata->sieve_buckets;
			
			if (root < linesize)			
			{	
				bnum = (root >> FLAGBITS);
				bptr[bnum][nptr[bnum]].root = root;
				bptr[bnum][nptr[bnum]].prime = prime;
				nptr[bnum]++;	
			}	
			
		}

	}

	return;
}

uint32 tiny_soe(uint32 limit, uint32 *primes)
{
	//simple sieve of erathosthenes for small limits - not efficient
	//for large limits.
	uint8 *flags;
	uint32 prime;
	uint32 i,j;
	int it;

	//allocate flags
	flags = (uint8 *)malloc(limit/2 * sizeof(uint8));
	if (flags == NULL)
		printf("error allocating flags\n");
	memset(flags,1,limit/2);

	//find the sieving primes, don't bother with offsets, we'll need to find those
	//separately for each line in the main sieve.
	primes[0] = 2;
	it=1;
	
	//sieve using primes less than the sqrt of block1
	//flags are created only for odd numbers (mod2)
	for (i=1;i<(uint32)(sqrt(limit)/2+1);i++)
	{
		if (flags[i] > 0)
		{
			prime = (uint32)(2*i + 1);
			for (j=i+prime;j<limit/2;j+=prime)
				flags[j]=0;

			primes[it]=prime;
			it++;
		}
	}

	//now find all the prime flags and compute the sieving primes
	//the last few will exceed uint16, we can fix this later.
	for (;i<limit/2;i++)
	{
		if (flags[i] == 1)
		{
			primes[it] = (uint32)(2*i + 1);
			it++;
		}
	}

	free(flags);
	return it;
}

void GetPRIMESRange(uint64 lowlimit, uint64 highlimit)
{
	uint64 i;
	
	//reallocate array based on conservative estimate of the number of 
	//primes in the interval
	
	i = (uint64)(highlimit/log(highlimit));
	if (lowlimit != 0)
		i -= (uint64)(lowlimit/log(lowlimit));
	i += (highlimit - lowlimit) * 1.2;

	PRIMES = (uint64 *)realloc(PRIMES,(size_t) (i*sizeof(uint64)));
	if (PRIMES == NULL)
	{
		printf("unable to allocate %lu bytes for range %lu to %lu\n",
			(uint64)(i*sizeof(uint64)),lowlimit,highlimit);
		exit(1);
	}

	//reset the global constants
	P_MIN = lowlimit;
	P_MAX = highlimit; 

	//find the primes in the interval
	NUM_P = spSOE(PRIMES,lowlimit,&highlimit,0);

	return;
}

uint64 soe_wrapper(uint64 lowlimit, uint64 highlimit, int count)
{
	//public interface to the sieve.  necessary because in order to keep the 
	//sieve efficient it must sieve larger blocks of numbers than a user may want,
	//and because the program keeps a cache of primes on hand which may or may 
	//not contain the range of interest.  Manage this on-hand cache and any addition
	//sieving needed.
	//TODO: manage really large requested blocks by splitting the range up into
	//blocks of no more than 10B else memory requirements get outragous.
	uint64 retval, tmpl, tmph,i=0;
	uint64 maxrange = 1000000000000;

	if (highlimit < lowlimit)
	{
		printf("error: lowlimit must be less than highlimit\n");
		return 0;
	}

	if (count)
	{
		//this needs to be a range of at least 1e6
		if ((highlimit - lowlimit) < 1000000)
		{
			//maybe it is already in our list of cached primes
			if ((lowlimit >= P_MIN) && (highlimit <= P_MAX))
			{
				retval = 0;
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						retval++;
				}
			}
			else
			{
				//nope, go and get a new range.
				tmpl = lowlimit;
				tmph = tmpl + 1000000;

				//since this is a small range, we need to 
				//find a bigger range and count them.
				GetPRIMESRange(tmpl,tmph);
				retval = 0;
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						retval++;
				}
			}
		}
		else
		{
			//check for really big ranges
			if ((highlimit - lowlimit) > maxrange)
			{
				uint32 num_ranges = (uint32)((highlimit - lowlimit) / maxrange);
				uint64 remainder = (highlimit - lowlimit) % maxrange;
				uint32 j;
				
				retval = 0;
				tmpl = lowlimit;
				tmph = lowlimit + maxrange;
				for (j = 0; j < num_ranges; j++)
				{
					retval += spSOE(NULL,tmpl,&tmph,1);
					tmpl += maxrange;
					tmph = tmpl + maxrange;
				}
				
				tmph = tmpl + remainder;
				retval += spSOE(NULL,tmpl,&tmph,1);
			}
			else
			{
				//we're in a sweet spot already, just get the requested range
				retval = spSOE(NULL,lowlimit,&highlimit,1);
			}
		}

	}
	else
	{
		if (lowlimit < P_MIN || lowlimit > P_MAX || highlimit > P_MAX)
		{
			//requested range is not covered by the current range
			tmpl = lowlimit;
			tmph = highlimit;

			//this needs to be a range of at least 1e6
			if (tmph - tmpl < 1000000)
			{
				//there is slack built into the sieve limit, so go ahead and increase
				//the size of the interval to make it at least 1e6.
				tmph = tmpl + 1000000;

				//since this is a small range, we need to 
				//find a bigger range and count them.
				GetPRIMESRange(tmpl,tmph);
				retval = 0;
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						retval++;
				}
			}
			else
			{
				//we don't need to mess with the requested range,
				//so GetPRIMESRange will return the requested range directly
				//and the count will be in NUM_P
				GetPRIMESRange(lowlimit,highlimit);
				retval = NUM_P;
			}
		}
		else
		{
			// the requested range is covered by the current range
			// just count them
			retval = 0;
			for (i = 0; i < NUM_P; i++)
			{
				if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
					retval++;
			}
		}

		// now dump the requested range of primes to a file, or the
		// screen, both, or neither, depending on the state of a couple
		// global configuration variables
		if (PRIMES_TO_FILE)
		{
			FILE *out;
			out = fopen("primes.dat","w");
			if (out == NULL)
			{
				printf("can't open primes.dat for writing\n");
			}
			else
			{
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						fprintf(out,"%lu\n",PRIMES[i]);
				}
				fclose(out);
			}
		}

		if (PRIMES_TO_SCREEN)
		{
			for (i = 0; i < NUM_P; i++)
			{
				if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
					printf("%lu ",PRIMES[i]);
			}
			printf("\n");
		}

			
	}

	return retval;
}

void primesum_check12(uint64 lower, uint64 upper, uint64 startmod, z *squaresum, z *sum)
{
	z mp1;
	//various counters
	uint64 n64;
	uint32 j;
	uint64 pcount = 0;
	
	//stuff for timing
	double t;
	struct timeval tstart, tstop;
	TIME_DIFF *	difference;

	//mananging batches
	uint64 inc, tmpupper=0;
	
	//tracking the modulus to check
	uint64 squaremod, summod;

	//logfile
	FILE *out;

	//for fast divisibility checks
	uint64 powof2sqr, powof2m1sqr, powof2sum, powof2m1sum;

	zInit(&mp1);

	//initialize both the sum and square modulus
	n64 = squaremod = summod = startmod;

	//set batch size based on input range
	if (upper - lower > 1000000000)
		inc = 1000000000;
	else
		inc = upper - lower;

	//paranoia.
	if (startmod == 0)
		startmod = 10;

	//count the number of factors of 2 in the modulus
	powof2sqr = 0;
	while ((n64 & 1) == 0)
	{
		n64 >>= 1;
		powof2sqr++;
	}
	//store this power minus 1, for fast remaindering
	powof2m1sqr = (1 << powof2sqr) - 1;

	//track of the powers of two in the sum and sqr modulus separately.
	powof2m1sum = powof2m1sqr;
	powof2sum = powof2sqr;

	//lets not print all this stuff to the screen...
	PRIMES_TO_SCREEN = 0;
	PRIMES_TO_FILE = 0;
	
	tmpupper = lower;
	while (tmpupper != upper)
	{
		//set the bounds for the next batch
		tmpupper = lower + inc; 
		if (tmpupper > upper)
			tmpupper =  upper;

		gettimeofday(&tstart, NULL);		

		//get a new batch of primes.
		n64 = soe_wrapper(lower,tmpupper,0);

		//keep a running sum of how many we've found
		pcount += n64;

		//status update
		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("\nfound %lu primes in range %lu to %lu in %6.4f sec\n",NUM_P,lower,tmpupper,t);
				
		gettimeofday(&tstart, NULL);
		//now add up the squares.  NUM_P is set by the soe_wrapper to the number of primes
		//found the last time it was called.
		for (j=0; j<NUM_P; j++)
		{
			//soe wrapper should truncate the prime array so that we don't have to do this,
			//but doublecheck that we don't include primes outside the requested batch.
			if ((PRIMES[j] > tmpupper) || (PRIMES[j] < lower))
				break;

#if defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)
			
			ASM_G (
				"addq %%rax, (%%rcx) \n\t"
				"adcq $0, 8(%%rcx) \n\t"
				"mulq %1		\n\t"
				"addq %%rax, (%%rbx)		\n\t"
				"adcq %%rdx, 8(%%rbx)		\n\t"
				"adcq $0, 16(%%rbx)	\n\t"
				: 
				: "b"(squaresum->val), "a"(PRIMES[j]), "c"(sum->val)
				: "rdx", "memory", "cc");				

#else
			sp642z(PRIMES[j],&mp1);
			zSqr(&mp1,&mp1);
			zAdd(squaresum,&mp1,squaresum);
			zShortAdd(sum,PRIMES[j],sum);

#endif
						
			//check if the squaresum is 0 modulo the current power of 10
			if ((squaresum->val[0] & powof2m1sqr) == 0)
			{
				//we have the right number of twos.  do a full check for divisibility.
				//first adjust the bigint size, since the ASM doesn't do this and zShortMod
				//would like size to be correct.
				squaresum->size = 3;
				while (squaresum->val[squaresum->size - 1] == 0)
				{
					squaresum->size--;
					if (squaresum->size == 0)
						break;
				}

				if (squaresum->size == 0)
					squaresum->size = 1;
				
				while (zShortMod(squaresum,squaremod) == 0)
				{
					out = fopen("sum_of_squares.csv","a");
					printf("**** %lu divides prime square sum up to %lu ****\n",squaremod,PRIMES[j]);
					fprintf(out,
						"**** %lu divides prime square sum up to %lu, sum is %s ****\n",
						squaremod,PRIMES[j],z2decstr(squaresum,&gstr1));
					fclose(out);
					//start looking for the next power of 10
					squaremod *= 10;
					//recompute the fast divisibility check
					powof2sqr++;
					powof2m1sqr = (1 << powof2sqr) - 1;
				}
			}			

			//now check the sum.
			if ((sum->val[0] & powof2m1sum) == 0)
			{
				//we have the right number of twos.  do a full check for divisibility.
				//adjust the size, since the ASM doesn't do this and zShortMod
				//would like size to be correct.
				sum->size = 2;
				while (sum->val[sum->size - 1] == 0)
				{
					sum->size--;
					if (sum->size == 0)
						break;
				}

				if (sum->size == 0)
					sum->size = 1;

				while (zShortMod(sum,summod) == 0)
				{
					out = fopen("sum_of_squares.csv","a");
					printf("**** %lu divides prime sum up to %lu ****\n",summod,PRIMES[j]);
					fprintf(out,
						"**** %lu divides prime sum up to %lu, sum is %s ****\n",
						summod,PRIMES[j],z2decstr(sum,&gstr1));
					fclose(out);
					summod *= 10;
					powof2sum++;
					powof2m1sum = (1 << powof2sum) - 1;
				}
			}
		}

		//all done with this batch.  status update again.  first get the sum
		//sizes right, for the benefit of z2decstr
		squaresum->size = 3;
		while (squaresum->val[squaresum->size - 1] == 0)
		{
			squaresum->size--;
			if (squaresum->size == 0)
				break;
		}

		if (squaresum->size == 0)
			squaresum->size = 1;

		sum->size = 2;
		while (sum->val[sum->size - 1] == 0)
		{
			sum->size--;
			if (sum->size == 0)
				break;
		}

		if (sum->size == 0)
			sum->size = 1;

		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("sum complete in %6.4f sec, squaresum = %s, sum = %s\n",
			t,z2decstr(squaresum,&gstr1),z2decstr(sum,&gstr2));
		
		out = fopen("sum_of_squares.csv","a");
		fprintf(out,"%lu,%lu,%lu,%s,%s\n",
			upper,n64,pcount,z2decstr(sum,&gstr1),z2decstr(squaresum,&gstr2));
		fclose(out);
		
		//next range.
		lower = tmpupper;
	}

	zFree(&mp1);
	return;
}

void primesum_check3(uint64 lower, uint64 upper, uint64 startmod, z *sum)
{
	z mp1;
	uint32 i=0;
	uint64 n64;
	uint32 j;
	uint64 count, tmpupper=0;
	
	double t;
	struct timeval tstart, tstop;
	TIME_DIFF *	difference;
	uint64 inc, summod, pcount = 0;

	FILE *out;
	uint64 powof2sum, powof2m1sum;

	zInit(&mp1);
	zClear(&mp1);

	//initialize both the sum and square modulus
	n64 = summod = startmod;

	if (upper - lower > 1000000000)
	{
		inc = 1000000000;
		count = (upper - lower) / inc;
	}
	else
	{
		inc = upper - lower;
		count = 0;
	}

	if (startmod == 0)
		startmod = 10;

	powof2sum = 0;
	//count the number of factors of 2 in the modulus
	while ((n64 & 1) == 0)
	{
		n64 >>= 1;
		powof2sum++;
	}
	powof2m1sum = (1 << powof2sum) - 1;

	PRIMES_TO_SCREEN = 0;
	PRIMES_TO_FILE = 0;
	
	//for each chunk
	for (i = 0; i < count; i++)
	{
		tmpupper = lower + inc; 
		gettimeofday(&tstart, NULL);				
		n64 = soe_wrapper(lower,tmpupper,0);
		pcount += n64;
		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("\nfound %lu primes in range %lu to %lu in %6.4f sec\n",NUM_P,lower,tmpupper,t);
		
		
		gettimeofday(&tstart, NULL);
		//now add up the squares
		for (j=0; j<NUM_P; j++)
		{

			if ((PRIMES[j] > tmpupper) || (PRIMES[j] < lower))
				break;

#if defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)
			
			/* 
				fast method to cube a number and sum with a 3-word fixed precision
				bigint

				                p
				*               p
				  ---------------
				          d     a
				*               p
				  ---------------
				  dpd dpa+apd apa
				+  s2    s1    s0
				  ---------------
				   s2    s1    s0


				where d,a are the high,low words of p^2,
				apd,apa are the high,low words of p * a and,
				dpd,dpa are the high,low words of p * d.
				we then sum these three words with our fixed precision cumulative sum s2,s1,s0:
				s0 = s0 + apa
				s1 = s1 + dpa + apd + carry
				s2 = s2 + dpd + carry
			*/


			ASM_G (
				"movq %1, %%rcx \n\t"			/* store prime */
				"mulq %%rcx		\n\t"			/* square it */
				"movq %%rax, %%r8 \n\t"			/* save p^2 lo (a) */
				"movq %%rdx, %%r9 \n\t"			/* save p^2 hi (d) */
				"mulq %%rcx	\n\t"				/* p * a */
				"movq %%rax, %%r10 \n\t"		/* save p*a lo (apa) */
				"movq %%rdx, %%r11 \n\t"		/* save p*a hi (apd) */
				"movq %%r9, %%rax \n\t"			/* p * d */
				"mulq %%rcx \n\t"				/* lo part in rax (dpa), hi in rdx (dpd) */
				"addq %%r10, (%%rbx) \n\t"		/* sum0 = sum0 + apa */
				"adcq %%rax, 8(%%rbx) \n\t"		/* sum1 = sum1 + dpa + carry */
				"adcq %%rdx, 16(%%rbx) \n\t"	/* sum2 = sum2 + dpd + carry */
				"addq %%r11, 8(%%rbx) \n\t"		/* sum1 = sum1 + apd */
				"adcq $0, 16(%%rbx)	\n\t"		/* sum2 = sum2 + carry */
				: 
				: "b"(sum->val), "a"(PRIMES[j])
				: "rcx", "rdx", "r8", "r9", "r10", "r11", "memory", "cc");			

#else
			sp642z(PRIMES[j],&mp1);
			zSqr(&mp1,&mp1);
			zShortMul(&mp1,PRIMES[j],&mp1);
			zAdd(sum,&mp1,sum);

#endif
						
			//check if the sum is 0 modulo the current power of 10
			if ((sum->val[0] & powof2m1sum) == 0)
			{
				//adjust the size, since the ASM doesn't do this and zShortMod
				//would like size to be correct.
				sum->size = 3;
				while (sum->val[sum->size - 1] == 0)
				{
					sum->size--;
					if (sum->size == 0)
						break;
				}

				if (sum->size == 0)
					sum->size = 1;

				//then we have the right number of twos.  check for divisibility.
				while (zShortMod(sum,summod) == 0)
				{
					out = fopen("sum_of_cubes.csv","a");
					printf("**** %lu divides prime cube sum up to %lu, sum = %s ****\n",
						summod,PRIMES[j],z2decstr(sum,&gstr1));
					fprintf(out,
						"**** %lu divides prime cube sum up to %lu, sum is %s ****\n",
						summod,PRIMES[j],z2decstr(sum,&gstr1));
					fclose(out);
					summod *= 10;
					powof2sum++;
					powof2m1sum = (1 << powof2sum) - 1;
				}
			}
		}

		sum->size = 3;
		while (sum->val[sum->size - 1] == 0)
		{
			sum->size--;
			if (sum->size == 0)
				break;
		}

		if (sum->size == 0)
			sum->size = 1;

		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("sum complete in %6.4f sec, sum = %s\n",
			t,z2decstr(sum,&gstr2));
		
		out = fopen("sum_of_cubes.csv","a");
		fprintf(out,"%lu,%lu,%lu,%s\n",
			tmpupper,n64,pcount,z2decstr(sum,&gstr1));
		fclose(out);
		
		lower = tmpupper;
		tmpupper += inc;
	}

	//printf("done looping.  lower = %lu, upper = %lu, tmpupper = %lu\n",lower,upper,tmpupper);
	if (upper - lower > 0)
	{
		gettimeofday(&tstart, NULL);				
		n64 = soe_wrapper(lower,upper,0);
		pcount += n64;
		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("\nfound %lu primes in range %lu to %lu in %6.4f sec\n",NUM_P,lower,upper,t);		
		
		gettimeofday(&tstart, NULL);
		//now add up the squares
		for (j=0; j<NUM_P; j++)
		{

			if ((PRIMES[j] > upper) || (PRIMES[j] < lower))
				break;

#if defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

			ASM_G (
				"movq %1, %%rcx \n\t"
				"mulq %%rcx		\n\t"
				"movq %%rax, %%r8 \n\t" /* pp lo (a) */
				"movq %%rdx, %%r9 \n\t" /* pp hi (d) */
				"mulq %%rcx		\n\t"
				"movq %%rax, %%r10 \n\t" /* ap lo */
				"movq %%rdx, %%r11 \n\t" /* ap hi */
				"movq %%r9, %%rax \n\t" 
				"mulq %%rcx \n\t" /* dp lo (rax), dp hi (rdx) */
				"addq %%r10, (%%rbx)		\n\t"
				"adcq %%rax, 8(%%rbx)		\n\t"
				"adcq %%rdx, 16(%%rbx)	\n\t" /* sum(2) += f + carry */
				"addq %%r11, 8(%%rbx) \n\t" /* sum(1) += e */
				"adcq $0, 16(%%rbx)	\n\t" /* sum(2) += f + carry */
				: 
				: "b"(sum->val), "a"(PRIMES[j])
				: "rcx", "rdx", "r8", "r9", "r10", "r11", "memory", "cc");					
				
#else
			sp642z(PRIMES[j],&mp1);
			zSqr(&mp1,&mp1);
			zShortMul(&mp1,PRIMES[j],&mp1);
			zAdd(sum,&mp1,sum);

#endif

			//printf("%d:%lu:%s\n",j,PRIMES[j],z2decstr(sum,&gstr1));
			//check if the sum is 0 modulo the current power of 10
			if ((sum->val[0] & powof2m1sum) == 0)
			{
				//adjust the size, since the ASM doesn't do this and zShortMod
				//would like size to be correct.
				sum->size = 3;
				while (sum->val[sum->size - 1] == 0)
				{
					sum->size--;
					if (sum->size == 0)
						break;
				}

				if (sum->size == 0)
					sum->size = 1;

				//then we have the right number of twos.  check for divisibility.
				while (zShortMod(sum,summod) == 0)
				{
					out = fopen("sum_of_cubes.csv","a");
					printf("**** %lu divides prime cube sum up to %lu, sum = %s ****\n",
						summod,PRIMES[j],z2decstr(sum,&gstr1));
					fprintf(out,
						"**** %lu divides prime cube sum up to %lu, sum is %s ****\n",
						summod,PRIMES[j],z2decstr(sum,&gstr1));
					fclose(out);
					summod *= 10;
					powof2sum++;
					powof2m1sum = (1 << powof2sum) - 1;
				}
			}

		}

		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("sum complete in %6.4f sec, sum = %s\n",
			t,z2decstr(sum,&gstr2));
	}

	zFree(&mp1);
	return;
}

void primesum(uint64 lower, uint64 upper)
{
	z mp1,mp2,mp3,*squaresum,*sum;
	uint64 n64;
	uint32 j;
	uint64 tmpupper=0;
	
	double t;
	struct timeval tstart, tstop;
	TIME_DIFF *	difference;
	uint64 inc, pcount = 0;

	zInit(&mp1);
	zInit(&mp2);
	zInit(&mp3);
	zClear(&mp1);
	squaresum = &mp2;
	sum = &mp3;

	//set batch size based on input range
	if (upper - lower > 1000000000)
		inc = 1000000000;
	else
		inc = upper - lower;

	tmpupper = lower;
	while (tmpupper != upper)
	{
		//set the bounds for the next batch
		tmpupper = lower + inc; 
		if (tmpupper > upper)
			tmpupper =  upper;

		gettimeofday(&tstart, NULL);				
		n64 = soe_wrapper(lower,tmpupper,0);
		pcount += n64;
		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("\nfound %lu primes in range %lu to %lu in %6.4f sec\n",NUM_P,lower,tmpupper,t);
		
		
		gettimeofday(&tstart, NULL);
		//now add up the squares
		for (j=0; j<NUM_P; j++)
		{

			if ((PRIMES[j] > tmpupper) || (PRIMES[j] < lower))
				break;

#if defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)
			
			ASM_G (
				"addq %%rax, (%%rcx) \n\t"
				"adcq $0, 8(%%rcx) \n\t"
				"mulq %1		\n\t"
				"addq %%rax, (%%rbx)		\n\t"
				"adcq %%rdx, 8(%%rbx)		\n\t"
				"adcq $0, 16(%%rbx)	\n\t"
				: 
				: "b"(squaresum->val), "a"(PRIMES[j]), "c"(sum->val)
				: "rdx", "memory", "cc");				

#else
			sp642z(PRIMES[j],&mp1);
			zSqr(&mp1,&mp1);
			zAdd(squaresum,&mp1,squaresum);
			zShortAdd(sum,PRIMES[j],sum);

#endif
			
		}

		squaresum->size = 3;
		while (squaresum->val[squaresum->size - 1] == 0)
		{
			squaresum->size--;
			if (squaresum->size == 0)
				break;
		}

		if (squaresum->size == 0)
			squaresum->size = 1;

		sum->size = 2;
		while (sum->val[sum->size - 1] == 0)
		{
			sum->size--;
			if (sum->size == 0)
				break;
		}

		if (sum->size == 0)
			sum->size = 1;

		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("sum complete in %6.4f sec, sum = %s, squaresum = %s\n",
			t,z2decstr(sum,&gstr1),z2decstr(squaresum,&gstr2));

		lower = tmpupper;
		tmpupper += inc;
	}

	zFree(&mp1);
	zFree(&mp2);
	zFree(&mp3);
	return;
}


#ifdef NOTDEF
int sieve_to_bitdepth(z *startn, z *stopn, int depth, uint32 *offsets)
{
	//for now, just use a fixed depth of 2^16

	//sieve an interval of numbers from start to stop, using
	//all primes less than 2^depth
	
	//masks for removing single bits in a byte
	const uint8 masks[8] = {0x7F, 0xBF, 0xDF, 0xEF, 0xF7, 0xFB, 0xFD, 0xFE};
	//masks for selecting single bits in a byte
	uint8 nmasks[8];
	
	//variables used for the wheel
	uint32 numclasses,prodN,startprime,lcount;

	//variables for blocking up the line structures
	uint32 numflags, numbytes, numlinebytes; //,lineblocks,partial_limit;

	//misc
	uint32 i,j,k,it,sp;
	uint16 prime;
	uint8 *flagblock;
	z lowlimit;
	z highlimit, w1;
	uint32 correction,offsettemp;
	uint32 range;

	//sieving structures
	soe_sieve16_t sieve16;

	soe_t soe;

	//timing variables
	clock_t start, stop;
	double t;
	//clock_t start2, stop2;
	double t2;

	start = clock();

	zInit(&w1);
	zInit(&lowlimit);
	zInit(&highlimit);
	numclasses=8;
	prodN=30;
	startprime=3;

	//allocate the residue classes.  
	soe.rclass = (uint8 *)malloc(numclasses * sizeof(uint8));

	//create the selection masks
	for (i=0;i<BITSINBYTE;i++)
		nmasks[i] = ~masks[i];
	
	//we'll need to store the offset into the next block for each prime
	sieve16.offsets = (uint16 *)malloc(6542 * sizeof(uint16));
	if (sieve16.offsets == NULL)
		printf("error allocating offsets\n");

	//store the primes used to sieve the rest of the flags
	//with the max sieve range set by the size of uint32, the number of primes
	//needed is fixed.
	sieve16.sieve_p = (uint16 *)malloc(6542 * sizeof(uint16));
	if (sieve16.sieve_p == NULL)
		printf("error allocating sieve_p\n");

	//use these primes
	//GetPRIMESRange(0,10000000);
	
	//get primes to sieve with
	sp = tiny_soe(65536,sieve16.sieve_p);
	soe.pboundi = 6542;

	//find the residue classes
	k=0;
	for (i=1;i<prodN;i++)
	{
		if (spGCD(i,prodN) == 1)
		{
			soe.rclass[k] = (uint8)i;
			k++;
		}
	}
	
	//temporarily set lowlimit to the first multiple of numclasses*prodN < lowlimit
	zShortDiv(startn,numclasses*prodN,&lowlimit);
	zShortMul(&lowlimit,numclasses*prodN,&lowlimit);
	zSub(startn,&lowlimit,&w1);
	correction = w1.val[0];
	//lowlimit = (lowlimit/(numclasses*prodN))*(numclasses*prodN);

	//reallocate flag structure for wheel and block sieving
	//starting at lowlimit, we need a flag for every 'numresidues' numbers out of 'prodN' up to 
	//limit.  round limit up to make this a whole number.
	//numflags = (*highlimit - lowlimit)/prodN;
	zSub(stopn,&lowlimit,&w1);
	if (w1.size > 1)
	{
		printf("interval too large\n");
		zFree(&w1);
		zFree(&lowlimit);
		zFree(&highlimit);
		return 0;
	}
	range = w1.val[0];
	numflags = range/prodN;
	numflags += ((numflags % prodN) != 0);
	numflags *= numclasses;

	//since we can pack 8 flags in a byte, we need numflags/8 bytes allocated.
	numbytes = numflags / BITSINBYTE + ((numflags % BITSINBYTE) != 0);

	//since there are 8 lines to sieve over, each line will contain (numflags/8)/8 bytes
	//so round numflags/8 up to the nearest multiple of 8
	numlinebytes = numbytes/numclasses + ((numbytes % numclasses) != 0);
	k = (uint32)((uint64)numlinebytes * (uint64)prodN * (uint64)BITSINBYTE);
	zShortAdd(&lowlimit,k,&highlimit);

	start = clock();
	//we will sieve over this many bytes of flags, for each line, block by block.
	//allocate the lines
	soe.line = (uint8 **)malloc(numclasses * sizeof(uint8 *));
	for (i=0;i<numclasses;i++)
	{
		soe.line[i] = (uint8 *)malloc(numlinebytes * sizeof(uint8));
		memset(soe.line[i],255,numlinebytes);
	}

	//a block consists of 32768 bytes of flags
	//which holds 262144 flags.
	sieve16.blocks = numlinebytes/BLOCKSIZE;
	sieve16.partial_block_b = (numlinebytes % BLOCKSIZE)*BITSINBYTE;

	if (0)
	{	
		printf("found %d residue classes\n",numclasses);
		printf("numlinebytes = %d\n",numlinebytes);
		printf("lineblocks = %d\n",sieve16.blocks);
		printf("partial limit = %d\n",sieve16.partial_block_b);
	}
	
	sieve16.blk_r = FLAGSIZE*prodN;
	it=0;

	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	if (0)
		printf("elapsed time for init = %6.5f\n",t);
	
	t=t2=0;
	//main sieve, line by line
	start = clock();
	for (lcount=0;lcount<numclasses;lcount++)
	{
		//for the current line, find the offsets past the low limit
		for (i=startprime;i<soe.pboundi;i++)
		{
			prime = sieve16.sieve_p[i];
			//find the first multiple of the prime which is greater than 'block1' and equal
			//to the residue class mod 'prodN'.  

			//solving the congruence: rclass[lcount] == kp mod prodN for k
			//use the eGCD?

			//tmp = sieve16.lblk_b/prime;
			//tmp *= prime;
			zShortDiv(&lowlimit,prime,&w1);
			zShortMul(&w1,prime,&w1);
			do
			{
				//tmp += prime;
				zShortAdd(&w1,prime,&w1);
			} while (zShortMod(&w1,prodN) != soe.rclass[lcount]); 
			//while ((tmp % prodN) != soe.rclass[lcount]);

			//now find out how much bigger this is than 'lowlimit'
			//in steps of 'prodN'.  this is exactly the offset in flags.
			zSub(&w1,&lowlimit,&w1);
			zShortDiv(&w1,prodN,&w1);
			//sieve16.offsets[i] = (uint16)((tmp - sieve16.lblk_b)/prodN);
			sieve16.offsets[i] = (uint16)(w1.val[0]);
		}

		//now we have primes and offsets for every prime, for this line.
		//proceed to sieve the entire flag set, block by block, for this line	
		flagblock = soe.line[lcount];
		for (i=0;i<sieve16.blocks;i++)
		{
			//sieve the block with each prime
			for (j=startprime;j<soe.pboundi;j++)
			{
				prime = sieve16.sieve_p[j];
				for (k=sieve16.offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				sieve16.offsets[j]= (uint16)(k - FLAGSIZE);
			}

			flagblock += BLOCKSIZE;
		}

		//and the last partial block, don't worry about updating the offsets
		if (sieve16.partial_block_b > 0)
		{
			for (i=startprime;i<soe.pboundi;i++)
			{
				prime = sieve16.sieve_p[i];
				for (j=sieve16.offsets[i];j<sieve16.partial_block_b;j+=prime)
					flagblock[j>>3] &= masks[j&7];
			}
		}
	}
	stop = clock();
	t += (double)(stop - start)/(double)CLOCKS_PER_SEC;

	//now we have a block of flags, not all of which will be prime depending
	//on the bitdepth and magnitude of the range being sieved.
	
	//rearrange the sieve lines so they are in order with respect to the original
	//starting value requested, and return a single array of offsets from start.

	it=0;
	for (j=0;j<numclasses;j++)
	{
		flagblock = soe.line[j];
		//offsets will not be ordered, but qsort could fix this if necessary
		//for a nominal fee.
		for (k=0; k<numlinebytes;k++)
		{
			//find out how many 30 counts this is away from the beginning (bits)
			//add the residue class
			//subtract the correction
			//result is the offset from start
			if ((flagblock[k] & nmasks[0]) >> 7)
			{
				offsettemp = prodN*(k*BITSINBYTE) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
			if ((flagblock[k] & nmasks[1]) >> 6)
			{
				offsettemp = prodN*(k*BITSINBYTE + 1) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
			if ((flagblock[k] & nmasks[2]) >> 5)
			{
				offsettemp = prodN*(k*BITSINBYTE + 2) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
			if ((flagblock[k] & nmasks[3]) >> 4)
			{
				offsettemp = prodN*(k*BITSINBYTE + 3) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
			if ((flagblock[k] & nmasks[4]) >> 3)
			{
				offsettemp = prodN*(k*BITSINBYTE + 4) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
			if ((flagblock[k] & nmasks[5]) >> 2)
			{
				offsettemp = prodN*(k*BITSINBYTE + 5) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
			if ((flagblock[k] & nmasks[6]) >> 1)
			{
				offsettemp = prodN*(k*BITSINBYTE + 6) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
			if ((flagblock[k] & nmasks[7]))
			{
				offsettemp = prodN*(k*BITSINBYTE + 7) + soe.rclass[j];
				if (offsettemp > correction)
					offsets[it++] = offsettemp - correction;
			}
		}
	}

	//printf("found %d primes\n",it);
	//offsets = (uint32 *)realloc(offsets,it*sizeof(uint32));

	zFree(&w1);
	zFree(&lowlimit);
	zFree(&highlimit);
	return it;
}
#endif
