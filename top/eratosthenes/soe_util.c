#include "soe.h"

uint64 estimate_primes_in_range(uint64 lowlimit, uint64 highlimit)
{
	uint64 hi_est, lo_est;

	hi_est = (uint64)(highlimit/log((double)highlimit));
	if (lowlimit > 1)
		lo_est = (uint64)(lowlimit/log((double)lowlimit));
	else
		lo_est = 0;

	return (uint64)((double)(hi_est - lo_est) * 1.2);
}

void get_numclasses(uint64 highlimit, uint64 lowlimit, soe_staticdata_t *sdata)
{
	uint64 numclasses, prodN, startprime;

	//more efficient to sieve using mod210 when the range is big
	if ((highlimit - lowlimit) > 40000000000ULL)
	{
		numclasses=480;
		prodN=2310;
		startprime=5;
	}	
	else if ((highlimit - lowlimit) > 4000000000ULL)
	{
		numclasses=48;
		prodN=210;
		startprime=4;		
	}
	else if ((highlimit - lowlimit) > 100000000)
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

	sdata->numclasses = numclasses;
	sdata->prodN = prodN;
	sdata->startprime = startprime;

	return;
}

int check_input(uint64 highlimit, uint64 lowlimit, uint32 num_sp, uint32 *sieve_p,
	soe_staticdata_t *sdata, mpz_t offset)
{
	int i;

	sdata->orig_hlimit = highlimit;
	sdata->orig_llimit = lowlimit;

	//the wrapper should handle this, but just in case we are called
	//directly and not via the wrapper...
	if (highlimit - lowlimit < 1000000)
		highlimit = lowlimit + 1000000;

	if ((highlimit - lowlimit) > 1000000000000ULL)
	{
		printf("range too big\n");
		return 1;
	}

	if (highlimit > 4000000000000000000ULL)
	{
		printf("input too high\n");
		return 1;
	}

	//set sieve primes in the local data structure to the ones that were passed in
	sdata->sieve_p = sieve_p;	
	
	if (offset == NULL)
	{
		//see if we were provided enough primes to do the job
		sdata->pbound = (uint64)(sqrt((int64)(highlimit)));

		if (sieve_p[num_sp - 1] < sdata->pbound)
		{
			printf("found %d primes, max = %u: not enough sieving primes\n", 
                num_sp, sieve_p[num_sp - 1]);

            //for (i = 0; i < num_sp; i++)
            //    printf("%u\n", sieve_p[i]);

			exit(1);
		}

		//find the highest index that we'll need.  Much of the rest of the code is 
		//sensitive to this.  Note that this could be slow for large numbers of
		//sieve primes... could replace with a binary search.
		for (i=0; i<num_sp; i++)
		{
			// stop when we have enough for this input
			if (sieve_p[i] > sdata->pbound)
				break;
		}
		sdata->pboundi = i;	
		
		sdata->offset = NULL;
		sdata->sieve_range = 0;
	}
	else
	{
		// for ranges with offsets, don't worry if we don't have enough
		// primes, but still check to see if we have too many.
		mpz_t tmpz;

		mpz_init(tmpz);
		mpz_add_ui(tmpz, offset, highlimit);
		mpz_sqrt(tmpz, tmpz);

		if (mpz_cmp_ui(tmpz, sieve_p[num_sp - 1]) < 0)
		{
			// then we were passed too many.  truncate the input list.
			sdata->pbound = mpz_get_64(tmpz);
			for (i=0; i<num_sp; i++)
			{
				// stop when we have enough for this input
				if (sieve_p[i] > sdata->pbound)
					break;
			}
			sdata->pboundi = i;	
		}
		else
		{
			// use all of 'em.
			sdata->pbound = sieve_p[num_sp - 1];
			sdata->pboundi = num_sp;
		}
		sdata->offset = offset;
		mpz_clear(tmpz);
		sdata->sieve_range = 1;
	}

	return 0;
}

uint64 init_sieve(soe_staticdata_t *sdata)
{
	int i,j,k;
	uint64 numclasses = sdata->numclasses;
	uint64 prodN = sdata->prodN;
	uint64 allocated_bytes = 0;
	uint64 lowlimit = sdata->orig_llimit; 
	uint64 highlimit = sdata->orig_hlimit;
	uint64 numflags, numbytes, numlinebytes;

	//create the selection masks
    for (i = 0; i < BITSINBYTE; i++)
    {
        nmasks[i] = ~masks[i];
    }

	//allocate the residue classes.  
	sdata->rclass = (uint32 *)malloc(numclasses * sizeof(uint32));
	allocated_bytes += numclasses * sizeof(uint32);

	//find the residue classes
	k=0;
	for (i=1;i<prodN;i++)
	{
		if (spGCD(i,(uint64)prodN) == 1)
		{
			sdata->rclass[k] = (uint32)i;
			k++;
		}
	}

	sdata->min_sieved_val = 1ULL << 63;
	
	//temporarily set lowlimit to the first multiple of numclasses*prodN < lowlimit
	if (sdata->sieve_range == 0)
	{
		lowlimit = (lowlimit/(numclasses*prodN))*(numclasses*prodN);
		sdata->lowlimit = lowlimit;
	}
	else
	{
		mpz_t tmpz, tmpz2;
		mpz_init(tmpz);
		mpz_init(tmpz2);

		//the start of the range of interest is controlled by offset, not lowlimit
		//figure out how it needs to change to accomodate sieving
		mpz_tdiv_q_ui(tmpz, *sdata->offset, numclasses * prodN);
		mpz_mul_ui(tmpz, tmpz, numclasses * prodN);		
		mpz_sub(tmpz2, *sdata->offset, tmpz);

		//raise the high limit by the amount the offset was lowered, so that
		//we allocate enough flags to cover the range of interest
		highlimit += mpz_get_ui(tmpz2);
		sdata->orig_hlimit += mpz_get_ui(tmpz2);

		//also raise the original lowlimit so that we don't include sieve primes
		//that we shouldn't when finalizing the process.
		sdata->orig_llimit += mpz_get_ui(tmpz2);

		//copy the new value to the pointer, which will get passed back to sieve_to_depth
		mpz_set(*sdata->offset, tmpz);
		mpz_clear(tmpz);
		mpz_clear(tmpz2);

		//set the lowlimit to 0; the real start of the range is controlled by offset
		sdata->lowlimit = 0;
	}

	//reallocate flag structure for wheel and block sieving
	//starting at lowlimit, we need a flag for every 'numresidues' numbers out of 'prodN' up to 
	//limit.  round limit up to make this a whole number.
	numflags = (highlimit - lowlimit)/prodN;
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
	//todo... did we round too much?  look into this.
	highlimit = (uint64)((uint64)numlinebytes * (uint64)prodN * (uint64)BITSINBYTE + lowlimit);
	sdata->highlimit = highlimit;
	sdata->numlinebytes = numlinebytes;

	//a block consists of BLOCKSIZE bytes of flags
	//which holds FLAGSIZE flags.
	sdata->blocks = numlinebytes/BLOCKSIZE;

	//each flag in a block is spaced prodN integers apart.  record the resulting size of the 
	//number line encoded in each block.
	sdata->blk_r = FLAGSIZE*prodN;

	//allocate space for the root of each sieve prime
	sdata->root = (int *)malloc(sdata->pboundi * sizeof(int));
	allocated_bytes += sdata->pboundi * sizeof(uint32);
	if (sdata->root == NULL)
	{
		printf("error allocating roots\n");
		exit(-1);
	}
	else
	{
		if (VFLAG > 2)
			printf("allocated %u bytes for roots\n",(uint32)(sdata->pboundi * sizeof(uint32)));
	}

	//compute the breakpoints at which we switch to other sieving methods	
	if (sdata->pboundi > BUCKETSTARTI)
	{
		sdata->bucket_start_id = BUCKETSTARTI;
		sdata->num_bucket_primes = sdata->pboundi - sdata->bucket_start_id;		
	}
	else
	{
		sdata->num_bucket_primes = 0;
		sdata->bucket_start_id = sdata->pboundi;
	}

	//any prime larger than this will only hit the interval once (in residue space)
	sdata->large_bucket_start_prime = sdata->blocks * FLAGSIZE;

	//block ranges for the various line sieving scenarios:
	//2 classes: block range = 1572864, approx min prime index = 119275
	//8 classes: block range = 7864320, min prime index = 531290
	//48 classes: block range = 55050240, min prime index = 3285304
	//480 classes: block range = 605552640, min prime index = 31599827
	sdata->inplace_start_id = sdata->pboundi;
	sdata->num_inplace_primes = 0;
#if defined(INPLACE_BUCKET)

	switch (sdata->numclasses)
	{
	case 2:
		if (sdata->pboundi > 119275)
		{
			sdata->inplace_start_id = 119275;
			sdata->num_inplace_primes = sdata->pboundi - sdata->inplace_start_id;		
			sdata->num_bucket_primes = sdata->inplace_start_id - sdata->bucket_start_id;
		}

		break;
	case 8:
		if (sdata->pboundi > 531290)
		{
			sdata->inplace_start_id = 531290;
			sdata->num_inplace_primes = sdata->pboundi - sdata->inplace_start_id;	
			sdata->num_bucket_primes = sdata->inplace_start_id - sdata->bucket_start_id;
		}

		break;
	case 48:
		if (sdata->pboundi > 3285304)
		{
			sdata->inplace_start_id = 3285304;
			sdata->num_inplace_primes = sdata->pboundi - sdata->inplace_start_id;
			sdata->num_bucket_primes = sdata->inplace_start_id - sdata->bucket_start_id;
		}

		break;
	case 480:
		if (sdata->pboundi > 31599827)
		{
			sdata->inplace_start_id = 31599827;
			sdata->num_inplace_primes = sdata->pboundi - sdata->inplace_start_id;
			sdata->num_bucket_primes = sdata->inplace_start_id - sdata->bucket_start_id;
		}

		break;
	default:
		printf("unknown number of classes\n");
		exit(1);

		break;
	}

	// allocate data structures for inplace sieving
	if (sdata->num_inplace_primes > 0)
	{	
		// set up a two dimensional array of pointers to elements of the sieving prime array
		// first dimension is block number
		// second dimension is class number
		// element is a pointer to a uint32 (a sieving prime)
		sdata->inplace_ptrs = (int **)malloc(sdata->blocks * sizeof(int *));
		for (i=0; i<sdata->blocks; i++)
		{
			int j;
			sdata->inplace_ptrs[i] = (int *)malloc(sdata->numclasses * sizeof(int));
			for (j=0; j<sdata->numclasses; j++)
				sdata->inplace_ptrs[i][j] = -1;
		}

		// allocate space for the inplace data we'll need
		sdata->inplace_data = (soe_inplace_p *)malloc(
			sdata->num_inplace_primes * sizeof(soe_inplace_p));
	}

#endif

	//these are only used by the bucket sieve
	//sdata->lower_mod_prime = (uint32 *)malloc(sdata->num_bucket_primes * sizeof(uint32));    
    //allocated_bytes += sdata->num_bucket_primes * sizeof(uint32);
    sdata->lower_mod_prime = (uint32 *)malloc(sdata->pboundi * sizeof(uint32));
    allocated_bytes += sdata->pboundi * sizeof(uint32);
	if (sdata->lower_mod_prime == NULL)
	{
		printf("error allocating lower mod prime\n");
		exit(-1);
	}
	else
	{
		if (VFLAG > 2)
			printf("allocated %u bytes for lower mod prime\n",
            (uint32)sdata->pboundi * (uint32)sizeof(uint32));
	}

	//allocate all of the lines if we are computing primes.  if we are
	//only counting them, just create the pointer array - lines will
	//be allocated as needed during sieving
	sdata->lines = (uint8 **)xmalloc_align(sdata->numclasses * sizeof(uint8 *));
	numbytes = 0;
	if (sdata->only_count)
	{
		//don't allocate anything now, but
		//provide an figure for the memory that will be allocated later
		numbytes = numlinebytes * sizeof(uint8) * THREADS;
	}
	else
	{
        //actually allocate all of the lines
        sdata->lines[0] = (uint8 *)xmalloc_align(numlinebytes * sdata->numclasses * sizeof(uint8));
        if (sdata->lines[0] == NULL)
        {
            printf("error allocated sieve lines\n");
            exit(-1);
        }
        numbytes += sdata->numclasses * numlinebytes * sizeof(uint8);

		for (i=0; i<sdata->numclasses; i++)
		{
            sdata->lines[i] = sdata->lines[0] + i * numlinebytes;
		}
	}

#ifdef USE_AVX2
    // during presieveing, storing precomputed lists will start to get unwieldy, so
    // generate the larger lists here.
    sdata->presieve_max_id = 40;

    for (j = 24; j < sdata->presieve_max_id; j++)
    {
        uint32 prime = sdata->sieve_p[j];

        // for each possible starting location
        for (i = 0; i < prime; i++)
        {
            uint64 interval[4] = { 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL };

            // sieve up to the bound, printing each 64-bit word as we fill it
            k = i;
            if (k >= 256) k -= 256;

            for (; k < 256; k += prime)
            {
                interval[k >> 6] &= ~(1ULL << k);
            }

            presieve_largemasks[j - 24][i][0] = interval[0];
            presieve_largemasks[j - 24][i][1] = interval[1];
            presieve_largemasks[j - 24][i][2] = interval[2];
            presieve_largemasks[j - 24][i][3] = interval[3];
        }        
    }

    for (j = 24; j < sdata->presieve_max_id; j++)
    {
        presieve_primes[j - 24] = sdata->sieve_p[j];
    }

    for (j = 24; j < sdata->presieve_max_id; j++)
    {
        presieve_p1[j - 24] = sdata->sieve_p[j] - 1;
    }

    for (j = 24; j < sdata->presieve_max_id; j++)
    {
        presieve_steps[j - 24] = 256 % sdata->sieve_p[j];
    }

#else
    sdata->presieve_max_id = 10;

#endif

	if (VFLAG > 2)
		printf("allocated %" PRIu64 " bytes for sieve lines\n",numbytes);
	allocated_bytes += numbytes;

	return allocated_bytes;
}

void set_bucket_depth(soe_staticdata_t *sdata)
{
	uint64 numlinebytes = sdata->numlinebytes;
	int i;

#if defined(INPLACE_BUCKET)
	//use normal bucket sieving until primes are larger than flagsize * prodN,
	//at which point we can use inplace sieving.
	//need to store prime div prodN and prime mod prodN, as well as maybe a roll-over
	//counter.  find the first hit as done below.  for each subsequent hit, find
	//the next residue class as done below.  find the next bit position by multiplying
	//steps_mod30 by prime div 30.  increment by the rollover, if necessary.  increment
	//rollover by prime mod 30.
	

		//uint32 **hitcount = (uint32 **)malloc(sdata->blocks * sizeof(uint32 *));
		//for (i=0; i<sdata->blocks; i++)
		//{
		//	int j;
		//	hitcount[i] = (uint32 *)malloc(sdata->numclasses * sizeof(uint32));
		//	for (j=0; j<sdata->numclasses; j++)
		//		hitcount[i][j] = 0;
		//}		

		//// now, for testing, rewind all of these pointers.  
		//// this is how we will sieve, eventually
		//for (i=0; i<sdata->blocks; i++)
		//{
		//	int j;
		//	for (j=0; j<sdata->numclasses; j++)
		//	{
		//		printf("block %d, class %d... ",i,j);

		//		if (sdata->inplace_ptrs[i][j] == -1)
		//			continue;
		//		while (sdata->next_p[sdata->inplace_ptrs[i][j]] != 0)
		//		{
		//			uint32 nextclass, pres;
		//			uint32 nextblock;
		//			uint64 nextbits, loc;
		//			int index = sdata->inplace_ptrs[i][j];
		//			int confirm;

		//			loc = sdata->prodN * ((i << FLAGBITS) + sdata->bitloc[index]) + 
		//				sdata->rclass[j] + sdata->lowlimit;

		//			if (loc % sdata->sieve_p[index] == 0)
		//			{
		//				confirm = 1;
		//				//printf("block %d, class %d, hit by prime %u, at flag %u; confirmed hit...\n",
		//				//	i,j,sdata->sieve_p[index],sdata->bitloc[index]);
		//			}
		//			else
		//			{
		//				confirm = 0;
		//				printf("block %d, class %d, hit by prime %u, at flag %u; error!!\n",
		//					i,j,sdata->sieve_p[index],sdata->bitloc[index]);	
		//				printf("computed loc was %" PRIu64 ", error is %u\n",
		//					loc,loc % sdata->sieve_p[index]);
		//				exit(1);
		//			}
		//			hitcount[i][j]++;
		//			
		//			//the next hit can be computed from the residue class progression
		//			//and the current prime's residue mod prodN
		//			//the correct progression is the row of residue_pattern_mod30 corresponding
		//			//to the current prime residue.  the current residue class gives the column.
		//			//the next column gives the next class hit, while the steps array gives
		//			//the next bit offset.
		//			if (sdata->sieve_p[index] == 393451)
		//			{
		//				int a;
		//				a = 2;
		//			}

		//			pres = resID_mod30[sdata->prime_res[index]];
		//			nextclass = residue_pattern_mod30[pres][j];
		//			nextbits = (uint64)sdata->bitloc[index] + (uint64)sdata->sieve_p[index] + 
		//				(uint64)steps_mod30[pres][j];
		//			nextblock = nextbits >> FLAGBITS;
		//			
		//			//follow the link to the next prime
		//			sdata->inplace_ptrs[i][j] = sdata->next_p[index];

		//			//if the next hit is within our interval
		//			if (nextblock < sdata->blocks)
		//			{
		//				//then reassign this prime to its next hit
		//				if (sdata->inplace_ptrs[nextblock][nextclass] == -1)
		//				{
		//					sdata->inplace_ptrs[nextblock][nextclass] = index;
		//					sdata->next_p[i] = 0;
		//				}
		//				else
		//				{
		//					sdata->next_p[i] = sdata->inplace_ptrs[nextblock][nextclass];
		//					sdata->inplace_ptrs[nextblock][nextclass] = index;						
		//				}
		//				sdata->bitloc[index] = nextbits & FLAGSIZEm1;
		//			}					
		//		}
		//		printf("processed %d hits\n",hitcount[i][j]);
		//	}
		//	
		//}

		//for (i=0; i<sdata->blocks; i++)
		//{
		//	int j;
		//	for (j=0; j<sdata->numclasses; j++)
		//		printf("hits in block %d, class %d = %u\n",i,j,hitcount[i][j]);
		//	free(hitcount[i]);
		//}
		//free(hitcount);

#endif

	if (sdata->num_bucket_primes > 0)
	{
		//then we have primes bigger than BUCKETSTARTP - need to bucket sieve
		uint64 flagsperline = numlinebytes * 8;
		uint64 num_hits = 0;
		uint64 hits_per_bucket;
		
		for (i = sdata->bucket_start_id; i < sdata->bucket_start_id + sdata->num_bucket_primes; i++)
		{
			//condition to see if the current prime only hits the sieve interval once
			if ((sdata->sieve_p[i] * sdata->prodN) > (sdata->blk_r * sdata->blocks))
				break;
			num_hits += ((uint32)flagsperline / sdata->sieve_p[i] + 1);
		}

		//assume hits are evenly distributed among buckets.
		hits_per_bucket = num_hits / sdata->blocks;

		//add some margin
		hits_per_bucket = (uint64)((double)hits_per_bucket * 1.10);

		//set the bucket allocation amount, with a minimum of at least 50000
		//because small allocation amounts may violate the uniformity assumption
		//of hits per bucket.  The idea is to set this right once, even if it is too big,
		//so that we don't have to keep checking for full buckets in the middle of
		//the bucket sieve (which would be slow)
		sdata->bucket_alloc = MAX(hits_per_bucket,50000);

		//now count primes that only hit the interval once
		num_hits = 0;
		for (; i < sdata->bucket_start_id + sdata->num_bucket_primes; i++)
			num_hits++;

		//assume hits are evenly distributed among buckets.
		hits_per_bucket = num_hits / sdata->blocks;

		//add some margin
		hits_per_bucket = (uint64)((double)hits_per_bucket * 1.1);

		if (num_hits > 0)
			sdata->large_bucket_alloc = MAX(hits_per_bucket,50000);
		else
			sdata->large_bucket_alloc = 0;

	}


	return;
}

uint64 alloc_threaddata(soe_staticdata_t *sdata, thread_soedata_t *thread_data)
{
	uint32 bucket_alloc = sdata->bucket_alloc;
	uint32 large_bucket_alloc = sdata->large_bucket_alloc;
	uint64 allocated_bytes = 0;
	uint32 bucket_depth = sdata->num_bucket_primes;
	int i,j;
	
	allocated_bytes += THREADS * sizeof(thread_soedata_t);
	for (i=0; i<THREADS; i++)
	{
		thread_soedata_t *thread = thread_data + i;

        // presieving scratch space
        thread->ddata.presieve_scratch = (uint32 *)xmalloc_align(8 * sizeof(uint32));

		//allocate a bound for each block
		thread->ddata.pbounds = (uint64 *)malloc(
			sdata->blocks * sizeof(uint64));
		allocated_bytes += sdata->blocks * sizeof(uint64);

		thread->ddata.pbounds[0] = sdata->pboundi;

		//we'll need to store the offset into the next block for each prime
		//actually only need those primes less than BUCKETSTARTP since bucket sieving
		//doesn't use the offset array.
		j = MIN(sdata->pboundi, BUCKETSTARTI);
        thread->ddata.offsets = (uint32 *)xmalloc_align(j * sizeof(uint32));
		allocated_bytes += j * sizeof(uint32);
		if (thread->ddata.offsets == NULL)
		{
			printf("error allocating offsets\n");
			exit(-1);
		}
		else
		{
			if (VFLAG > 2)
				printf("allocated %u bytes for offsets for %d sieving primes \n",
					(uint32)(j * sizeof(uint32)), j);
		}
		thread->ddata.bucket_depth = 0;

		if (bucket_depth > 0)
		{			
			//create a bucket for each block
			//thread->ddata.sieve_buckets = (soe_bucket_t **)malloc(
			//	sdata->blocks * sizeof(soe_bucket_t *));
			//allocated_bytes += sdata->blocks * sizeof(soe_bucket_t *);
            thread->ddata.sieve_buckets = (uint64 **)malloc(
                sdata->blocks * sizeof(uint64 *));
            allocated_bytes += sdata->blocks * sizeof(uint64 *);

			if (thread->ddata.sieve_buckets == NULL)
			{
				printf("error allocating buckets\n");
				exit(-1);
			}
			else
			{
                if (VFLAG > 2)
                {
                    printf("allocated %u bytes for bucket bases\n",
                        (uint32)sdata->blocks * (uint32)sizeof(uint64 *));
                }
			}

			if (large_bucket_alloc > 0)
			{
				thread->ddata.large_sieve_buckets = (uint32 **)malloc(
					sdata->blocks * sizeof(uint32 *));
				allocated_bytes += sdata->blocks * sizeof(uint32 *);

				if (thread->ddata.large_sieve_buckets == NULL)
				{
					printf("error allocating large buckets\n");
					exit(-1);
				}
				else
				{
					if (VFLAG > 2)
						printf("allocated %u bytes for large bucket bases\n",
							(uint32)sdata->blocks * (uint32)sizeof(uint32 *));
				}
			}
            else
            {
                thread->ddata.large_sieve_buckets = NULL;
            }

			//create a hit counter for each bucket
			thread->ddata.bucket_hits = (uint32 *)malloc(
				sdata->blocks * sizeof(uint32));
			allocated_bytes += sdata->blocks * sizeof(uint32);
			if (thread->ddata.bucket_hits == NULL)
			{
				printf("error allocating hit counters\n");
				exit(-1);
			}
			else
			{
                if (VFLAG > 2)
                {
                    printf("allocated %u bytes for hit counters\n",
                        (uint32)sdata->blocks * (uint32)sizeof(uint32));
                }
			}

			if (large_bucket_alloc > 0)
			{
				thread->ddata.large_bucket_hits = (uint32 *)malloc(
					sdata->blocks * sizeof(uint32));
				allocated_bytes += sdata->blocks * sizeof(uint32);
				if (thread->ddata.large_bucket_hits == NULL)
				{
					printf("error allocating large hit counters\n");
					exit(-1);
				}
				else
				{
                    if (VFLAG > 2)
                    {
                        printf("allocated %u bytes for large hit counters\n",
                            (uint32)sdata->blocks * (uint32)sizeof(uint32));
                    }
				}
			}
            else
            {
                thread->ddata.large_bucket_hits = NULL;
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

			for (j = 0; j < sdata->blocks; j++)
			{
				//thread->ddata.sieve_buckets[j] = (soe_bucket_t *)malloc(
				//	bucket_alloc * sizeof(soe_bucket_t));
				//allocated_bytes += bucket_alloc * sizeof(soe_bucket_t);

                thread->ddata.sieve_buckets[j] = (uint64 *)malloc(
                    bucket_alloc * sizeof(uint64));
                allocated_bytes += bucket_alloc * sizeof(uint64);

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
            {
                //printf("allocated %u bytes for buckets\n",
                  //  (uint32)sdata->blocks * (uint32)bucket_alloc * (uint32)sizeof(soe_bucket_t));
                printf("allocated %u bytes for buckets\n",
                    (uint32)sdata->blocks * (uint32)bucket_alloc * (uint32)sizeof(uint64));
            }

            if (VFLAG > 2)
            {
                printf("allocated %u bytes for large buckets\n",
                    (uint32)sdata->blocks * (uint32)large_bucket_alloc * (uint32)sizeof(uint32));
            }

		}	

		//this threads' count of primes in its' line
		thread->linecount = 0;
		//share the common static data structure
		thread->sdata = *sdata;
	}

	return allocated_bytes;
}

void build_table()
{
	// output the sequence of residues and sequence of steps for
	// each residue class
	int i,j,prodn,bytenum,resnum;
	uint64 resword[64];
	int res_pattern[8][8] = {
		{1, 7, 11, 13, 17, 19, 23, 29},
		{7, 19, 17, 1, 29, 13, 11, 23},
		{11, 17, 1, 23, 7, 29, 13, 19},
		{13, 1, 23, 19, 11, 7, 29, 17},
		{17, 29, 7, 11, 19, 23, 1, 13},
		{19, 13, 29, 7, 23, 1, 17, 11},
		{23, 11, 13, 29, 1, 17, 19, 7},
		{29, 23, 19, 17, 13, 11, 7, 1}};
	int steps[8] = {6,4,2,4,2,4,6,2};
	int res[8] = {1, 7, 11, 13, 17, 19, 23, 29};


	// build a table which, given the current class and a residue class with which
	// to increment by, returns the next residue class in a 8 class system as well
	// as the number of steps it will take to get there.  the next class is in the 
	// low order byte and the number of steps is in the high order byte of a uint16.
	// this is useful to initialize the inplace sieve.  see roots.c
	printf("\n{");
	for (i=0; i<8; i++)
	{
		// given increment residue class i (i.e., 1, 7, 11, 13, etc...)		
		for (j=0; j<30; j++)
		{
			// and starting residue mod 30 == j
			uint32 val = j;
			uint32 steps = 0;
			
			for (;;)
			{
				// print the next valid residue in (1,7,11,13, etc.)
				int k, breakout = 0;
				
				for (k=0; k<8; k++)
				{
					if (val == res[k])
					{
						breakout = 1;
						break;
					}
				}
				if (breakout)
					break;

				val = (val + res[i]) % 30;
				steps++;
			}
			//printf("(%u,%u),",val,steps);
			printf("%u,",val + (steps << 8));
		}
		printf("},\n");
	}
	printf("};\n\n");
	exit(0);


	printf("\n{");
	for (i=0; i<8; i++)
	{
		for (j=0; j<8; j++)
		{
			int k;
			for (k=0; k<8; k++)
			{
				if (res_pattern[i][k] == res[j])
				{
					if (k < 7)
						printf("%d, ",res_pattern[i][k+1]);
					else
						printf("%d, ",res_pattern[i][0]);
				}
			}
		}
		printf("}\n{");
	}
	printf("};\n\n");

	printf("\n{");
	for (i=0; i<8; i++)
	{
		for (j=0; j<8; j++)
		{
			int k;
			for (k=0; k<8; k++)
			{
				if (res_pattern[i][k] == res[j])
				{
					printf("%d, ",steps[k]);
				}
			}
		}
		printf("}\n{");
	}
	printf("};\n\n");

	exit(0);



	prodn = 30;
	bytenum = 0;
	resnum = 0;
	for (i=0; i<64; i++)
		resword[i] = 0;

	printf("\n");
	for (i=1; i<prodn; i++)
	{
		if (spGCD(i,prodn) == 1)
		{
			int steps = 0;
			printf("\n");
			
			for (j=i; j<=prodn*i; j += i)
			{
				int n = j % prodn;
				if (spGCD(n,prodn) == 1)
				{
					resword[resnum] |= (uint64)n << (bytenum * 8);
					bytenum++;
					printf("%02d,%d ",n,steps);
					steps = 0;
				}
				steps++;
			}
			resnum++;
		}
	}
	printf("\n");
	for (i=0; i<resnum; i++)
	{
		printf("%" PRIu64 " ", resword[i]);
		resword[i] = 0;
	}
	printf("\n");

	prodn = 210;
	bytenum = 0;
	resnum = 0;
	printf("\n");
	for (i=1; i<prodn; i++)
	{
		if (spGCD(i,prodn) == 1)
		{
			int steps = 0;
			printf("\n");
			
			for (j=i; j<=prodn*i; j += i)
			{
				int n = j % prodn;
				if (spGCD(n,prodn) == 1)
				{
					resword[resnum] |= (uint64)n << (bytenum * 8);
					bytenum++;
					printf("%03d,%d ",n,steps);
					steps = 0;
				}
				steps++;
			}
			resnum++;
		}
	}
	printf("\n");
	for (i=0; i<resnum; i++)
		printf("%" PRIu64 " ", resword[i]);
	printf("\n");

	return;
}


