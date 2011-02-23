/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 7/1/10
----------------------------------------------------------------------*/

#include "soe.h"

uint64 spSOE(uint64 *primes, uint64 lowlimit, uint64 *highlimit, int count)
{
	/*
	if count == 1, then the primes are simply counted, and not 
	explicitly calculated and saved in *primes.
	*/

	//variables used for the wheel
	uint64 numclasses,prodN,startprime;
	uint32 bucket_depth;

	//variables for blocking up the line structures
	uint64 numflags, numbytes, numlinebytes, bucket_alloc, large_bucket_alloc;

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

	//*********************** BEGIN ******************************//

	sdata.orig_hlimit = *highlimit;
	sdata.orig_llimit = lowlimit;

	if (*highlimit - lowlimit < 1000000)
		*highlimit = lowlimit + 1000000;

	if (*highlimit - lowlimit > 1000000000000ULL)
	{
		printf("range too big\n");
		return 0;
	}

	//more efficient to sieve using mod210 when the range is big
	if ((*highlimit - lowlimit) > 400000000000ULL)
	{
		//numclasses=5760;
		//prodN=30030;
		//startprime=6;
		numclasses=480;
		prodN=2310;
		startprime=5;
	}	
	else if ((*highlimit - lowlimit) > 40000000000ULL)
	{
		numclasses=480;
		prodN=2310;
		startprime=5;
	}	
	else if ((*highlimit - lowlimit) > 4000000000ULL)
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

	//initialize
	sdata.inplace_startprime = FLAGSIZE * prodN;

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
	if (*highlimit > 4000000000000000000ULL)
	{
		printf("input too high\n");
		return 0;
	}

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
		locprimes = (uint64 *)realloc(locprimes,k * sizeof(uint64));	
		if (locprimes == NULL)
		{
			printf("error allocating locprimes\n");
			exit(-1);
		}

		// recursion
		sp = (uint32)spSOE(locprimes,0,&j,0);

		//check if we've found too many
		if (sp > MAXSIEVEPRIMECOUNT)
		{
			printf("input too high\n");
			free(locprimes);
			free(mergeprimes);
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
					printf("allocated %u bytes for sieving primes\n",sp * (uint32)sizeof(uint32));
			}
			
		}

		// copy into the 32 bit sieving prime array
		for (k=0; k<sp; k++)
		{
			if (locprimes[k] == 0)
				printf("found prime == 0 in locprimes at location %" PRIu64 "\n",k);
			sdata.sieve_p[k] = (uint32)locprimes[k];
		}

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
				printf("allocated %u bytes for sieving primes\n",78498 * (uint32)sizeof(uint32));
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
			printf("estimating storage for primes up to %" PRIu64 "\n",*highlimit);
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

	sdata.valid_residue = (uint32 *)malloc(prodN * sizeof(uint32));

	//find the residue classes
	k=0;
	sdata.valid_residue[0] = 0;
	for (i=1;i<prodN;i++)
	{
		if (spGCD(i,(fp_digit)prodN) == 1)
		{
			sdata.rclass[k] = (uint32)i;
			printf("%u ",i);
			k++;
			sdata.valid_residue[i] = i;
		}
		else
			sdata.valid_residue[i] = 0;
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

	//since there are N lines to sieve over, each line will contain (numflags/8)/N bytes
	//so round numflags/8 up to the nearest multiple of N
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

	//all this rounding has likely changed the desired high limit.  compute the new highlimit.
	//the orignial desired high limit is already recorded so the proper count will be returned.
	*highlimit = (uint64)((uint64)numlinebytes * (uint64)prodN * (uint64)BITSINBYTE + lowlimit);
	sdata.highlimit = *highlimit;
	sdata.numlinebytes = numlinebytes;

	//a block consists of BLOCKSIZE bytes of flags
	//which holds FLAGSIZE flags.
	sdata.blocks = numlinebytes/BLOCKSIZE;

	//each flag in a block is spaced prodN integers apart.  record the resulting size of the 
	//number line encoded in each block.
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
			printf("allocated %u bytes for roots\n",(uint32)(sdata.pboundi * sizeof(uint32)));
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
	
	//initialize.  defaults to this value if inplace sieving is not available
	sdata.inplace_startindex = sdata.pboundi;

	if (sdata.pboundi > BUCKETSTARTI)
	{
		//then we have primes bigger than BUCKETSTARTP - need to bucket sieve
		uint64 flagsperline = numlinebytes * 8;
		uint64 num_hits = 0;
		uint64 hits_per_bucket;
		
		for (i=BUCKETSTARTI; i<sdata.pboundi; i++)
		{

#ifdef INPLACE_BUCKET
			if (sdata.sieve_p[i] > sdata.inplace_startprime)
			{
				//at this point and above, we can use in-place bucket sieving
				sdata.inplace_startindex = i;
				break;
			}
#endif

			//condition to see if the current prime only hits the sieve interval once
			if ((sdata.sieve_p[i] * sdata.prodN) > (sdata.blk_r * sdata.blocks))
				break;
			num_hits += ((uint32)flagsperline / sdata.sieve_p[i] + 1);
		}

		//assume hits are evenly distributed among buckets.
		hits_per_bucket = num_hits / sdata.blocks;

		//add some margin
		hits_per_bucket = (uint64)((double)hits_per_bucket * 1.10);

		//set the bucket allocation amount
		bucket_alloc = hits_per_bucket;

		num_hits = 0;
		//printf("%u primes above large prime threshold\n",(uint32)(sdata.pboundi - i));
		for (; i<sdata.pboundi; i++)
		{
#ifdef INPLACE_BUCKET
			if (sdata.sieve_p[i] > sdata.inplace_startprime)
			{
				//at this point and above, we can use in-place bucket sieving
				sdata.inplace_startindex = i;
				break;
			}
#endif
			num_hits++;
		}

		//assume hits are evenly distributed among buckets.
		hits_per_bucket = num_hits / sdata.blocks;

		//add some margin
		hits_per_bucket = (uint64)((double)hits_per_bucket * 1.1);

		//set the bucket allocation amount, with a minimum of at least 25000
		//because small allocation amounts may violate the uniformity assumption
		//of hits per bucket
		if (num_hits > 0)
			large_bucket_alloc = MAX(hits_per_bucket,50000);
		else
			large_bucket_alloc = 0;

		bucket_depth = (uint32)(sdata.inplace_startindex - BUCKETSTARTI);
	}
	else
	{
		bucket_depth = 0;
	}

	//uncomment this to disable bucket sieving
	//bucket_depth = 0;

#ifdef INPLACE_BUCKET
	//force this to be single threaded for now
	THREADS = 1;
#endif

	thread_data = (thread_soedata_t *)malloc(THREADS * sizeof(thread_soedata_t));
	allocated_bytes += THREADS * sizeof(thread_soedata_t);
	for (i=0; i<THREADS; i++)
	{
		thread_soedata_t *thread = thread_data + i;

		//allocate a bound for each block
		thread->ddata.pbounds = (uint64 *)malloc(
			sdata.blocks * sizeof(uint64));
		allocated_bytes += sdata.blocks * sizeof(uint64);

		thread->ddata.pbounds[0] = sdata.pboundi;

		//we'll need to store the offset into the next block for each prime
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
				printf("allocated %u bytes for offsets\n",(uint32)(j * sizeof(uint32)));
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
				printf("allocated %u bytes for sieve line\n",
				(uint32)numlinebytes * (uint32)sizeof(uint8));
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
				printf("allocated %u bytes for special line count\n",(uint32)(j * sizeof(uint32)));
		}
#endif

		if (bucket_depth > BUCKET_BUFFER)
		{			

#ifdef INPLACE_BUCKET

			//allocate linked list pointers
			sdata.listptrs = (bucket_prime_t **)malloc(
				sdata.blocks * sdata.prodN * sizeof(bucket_prime_t *));

			for (j=0; j < sdata.blocks * sdata.prodN; j++)
				sdata.listptrs[j] = NULL;

			//allocate bucket_primes
			sdata.bucket_primes = (bucket_prime_t *)malloc(
				(sdata.pboundi - BUCKETSTARTI) * sizeof(bucket_prime_t));

			if (VFLAG > 2)
				printf("allocated %u bytes for bucket pointers\n",
				(uint32)sdata.blocks * sdata.prodN * sizeof(bucket_prime_t *));

			if (VFLAG > 2)
				printf("allocated %u bytes for bucket primes\n",
				(sdata.pboundi - BUCKETSTARTI) * sizeof(bucket_prime_t));

#endif
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
					printf("allocated %u bytes for bucket bases\n",
					(uint32)sdata.blocks * (uint32)sizeof(soe_bucket_t *));
			}

			if (large_bucket_alloc > 0)
			{
				thread->ddata.large_sieve_buckets = (uint32 **)malloc(
					sdata.blocks * sizeof(uint32 *));
				allocated_bytes += sdata.blocks * sizeof(uint32 *);

				if (thread->ddata.large_sieve_buckets == NULL)
				{
					printf("error allocating large buckets\n");
					exit(-1);
				}
				else
				{
					if (VFLAG > 2)
						printf("allocated %u bytes for large bucket bases\n",
						(uint32)sdata.blocks * (uint32)sizeof(uint32 *));
				}
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
					printf("allocated %u bytes for hit counters\n",
					(uint32)sdata.blocks * (uint32)sizeof(uint32));
			}

			if (large_bucket_alloc > 0)
			{
				thread->ddata.large_bucket_hits = (uint32 *)malloc(
					sdata.blocks * sizeof(uint32));
				allocated_bytes += sdata.blocks * sizeof(uint32);
				if (thread->ddata.large_bucket_hits == NULL)
				{
					printf("error allocating large hit counters\n");
					exit(-1);
				}
				else
				{
					if (VFLAG > 2)
						printf("allocated %u bytes for large hit counters\n",
						(uint32)sdata.blocks * (uint32)sizeof(uint32));
				}
			}

			//each bucket must be able to hold a hit from every prime used above BUCKETSTARTP.
			//this is overkill, because every prime will not hit every bucket when
			//primes are greater than FLAGSIZE.  but we always write to the end of each
			//bucket so the depth probably doesn't matter much from a cache access standpoint,
			//just from a memory capacity standpoint, and we shouldn't be in any danger
			//of that as long as the input range is managed (should be <= 10B).
			thread->ddata.bucket_depth = bucket_depth;
			thread->ddata.bucket_alloc = bucket_alloc;
			thread->ddata.bucket_alloc_large = large_bucket_alloc;

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

				if (large_bucket_alloc > 0)
				{
					thread->ddata.large_sieve_buckets[j] = (uint32 *)malloc(
						large_bucket_alloc * sizeof(uint32));
					allocated_bytes += large_bucket_alloc * sizeof(uint32);

					if (thread->ddata.large_sieve_buckets[j] == NULL)
					{
						printf("error allocating large buckets\n");
						exit(-1);
					}	

					thread->ddata.large_bucket_hits[j] = 0;
				}
							
			}

			if (VFLAG > 2)
				printf("allocated %u bytes for buckets\n",
				(uint32)sdata.blocks * (uint32)bucket_alloc * (uint32)sizeof(soe_bucket_t));

			if (VFLAG > 2)
				printf("allocated %u bytes for large buckets\n",
				(uint32)sdata.blocks * (uint32)large_bucket_alloc * (uint32)sizeof(uint32));

		}	
		else
			thread->ddata.bucket_depth = 0;

		//this threads' count of primes in its' line
		thread->linecount = 0;
		//share the common static data structure
		thread->sdata = sdata;
	}

	gettimeofday (&tstart, NULL);
	getRoots(&sdata);

	if (VFLAG > 2)
	{
		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);

		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);

		printf("elapsed time for computing roots = %6.4f\n",t);
	}

	if (VFLAG > 2)
	{	
		printf("sieving range %" PRIu64 " to %" PRIu64 "\n",lowlimit,*highlimit);
		printf("using %" PRIu64 " primes, max prime = %" PRIu64 "  \n",sdata.pboundi,sdata.pbound);
		printf("using %" PRIu64 " residue classes\n",numclasses);
		printf("lines have %" PRIu64 " bytes and %" PRIu64 " flags\n",numlinebytes,numlinebytes * 8);
		printf("lines broken into = %" PRIu64 " blocks of size %u\n",sdata.blocks,BLOCKSIZE);
		printf("blocks contain %u flags and cover %" PRIu64 " primes\n", FLAGSIZE, sdata.blk_r);
		if (bucket_depth > BUCKET_BUFFER)
		{
			printf("bucket sieving %u primes > %u\n",bucket_depth,BUCKETSTARTP);
			printf("allocating space for %" PRIu64 " hits per bucket\n",bucket_alloc);

#ifdef INPLACE_BUCKET
			printf("bucket sieving %u primes in-place\n",sdata.pboundi - sdata.inplace_startindex);
#endif

#ifdef DO_LARGE_BUCKETS
			printf("allocating space for %u hits per large bucket\n",large_bucket_alloc);
#endif
		}
		printf("using %" PRIu64 " bytes for sieving storage\n",allocated_bytes);
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
				//uint64 start_index = num_p;
				uint64 linecount = thread_data[i].linecount;
				uint64 index;
				
				if (linecount == 0)
				{
					printf("found no primes in line\n");
					continue;
				}
				
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

					//increase size of mergeprimes
					i1 = i2 = i3 = 0;
					while ((i1 < num_p) && (i2 < linecount)) {
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
					{
						locprimes[i1] = mergeprimes[i1];
					}

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
		//uint64 start_index = num_p;
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
		while ((i1 < num_p) && (i2 < it)) {
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
	{
		if (sdata.special_count[i] > 0)
			printf("count in range %d = %u\n",(int)i,sdata.special_count[i]);
	}
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

#ifdef INPLACE_BUCKET


#else
		for (i=0; i< THREADS; i++)
		{
			thread_soedata_t *thread = thread_data + i;

			free(thread->ddata.bucket_hits);
#ifdef DO_LARGE_BUCKETS
			free(thread->ddata.large_bucket_hits);
#endif
			for (j=0; j < thread->sdata.blocks; j++)
			{
				free(thread->ddata.sieve_buckets[j]);
#ifdef DO_LARGE_BUCKETS
				free(thread->ddata.large_sieve_buckets[j]);
#endif
			}
			free(thread->ddata.sieve_buckets);
#ifdef DO_LARGE_BUCKETS
			free(thread->ddata.large_sieve_buckets);
#endif
		}

#endif
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
		printf("\nfound %" PRIu64 " primes in range %" PRIu64 " to %" PRIu64 " in %6.4f sec\n",NUM_P,lower,tmpupper,t);
				
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
					printf("**** %" PRIu64 " divides prime square sum up to %" PRIu64 " ****\n",squaremod,PRIMES[j]);
					fprintf(out,
						"**** %" PRIu64 " divides prime square sum up to %" PRIu64 ", sum is %s ****\n",
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
					printf("**** %" PRIu64 " divides prime sum up to %" PRIu64 " ****\n",summod,PRIMES[j]);
					fprintf(out,
						"**** %" PRIu64 " divides prime sum up to %" PRIu64 ", sum is %s ****\n",
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
		fprintf(out,"%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%s,%s\n",
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
		printf("\nfound %" PRIu64 " primes in range %" PRIu64 " to %" PRIu64 " in %6.4f sec\n",NUM_P,lower,tmpupper,t);
		
		
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
					printf("**** %" PRIu64 " divides prime cube sum up to %" PRIu64 ", sum = %s ****\n",
						summod,PRIMES[j],z2decstr(sum,&gstr1));
					fprintf(out,
						"**** %" PRIu64 " divides prime cube sum up to %" PRIu64 ", sum is %s ****\n",
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
		fprintf(out,"%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%s\n",
			tmpupper,n64,pcount,z2decstr(sum,&gstr1));
		fclose(out);
		
		lower = tmpupper;
		tmpupper += inc;
	}

	//printf("done looping.  lower = %" PRIu64 ", upper = %" PRIu64 ", tmpupper = %" PRIu64 "\n",lower,upper,tmpupper);
	if (upper - lower > 0)
	{
		gettimeofday(&tstart, NULL);				
		n64 = soe_wrapper(lower,upper,0);
		pcount += n64;
		gettimeofday (&tstop, NULL);
		difference = my_difftime (&tstart, &tstop);
		t = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		printf("\nfound %" PRIu64 " primes in range %" PRIu64 " to %" PRIu64 " in %6.4f sec\n",NUM_P,lower,upper,t);		
		
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

			//printf("%d:%" PRIu64 ":%s\n",j,PRIMES[j],z2decstr(sum,&gstr1));
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
					printf("**** %" PRIu64 " divides prime cube sum up to %" PRIu64 ", sum = %s ****\n",
						summod,PRIMES[j],z2decstr(sum,&gstr1));
					fprintf(out,
						"**** %" PRIu64 " divides prime cube sum up to %" PRIu64 ", sum is %s ****\n",
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
		printf("\nfound %" PRIu64 " primes in range %" PRIu64 " to %" PRIu64 " in %6.4f sec\n",NUM_P,lower,tmpupper,t);
		
		
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

