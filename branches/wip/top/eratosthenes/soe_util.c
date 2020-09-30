#include "soe.h"
#include <immintrin.h>

#ifdef USE_AVX512F
ALIGNED_MEM uint64 presieve_largemasks[16][173][8];
ALIGNED_MEM uint32 presieve_steps[32];
ALIGNED_MEM uint32 presieve_primes[32];
ALIGNED_MEM uint32 presieve_p1[32];

#else
// for storage of presieving lists from prime index 24 to 40 (97 to 173 inclusive)
ALIGNED_MEM uint64 presieve_largemasks[16][173][4];
ALIGNED_MEM uint32 presieve_steps[32];
ALIGNED_MEM uint32 presieve_primes[32];
ALIGNED_MEM uint32 presieve_p1[32];
#endif
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

// row: numclasses 2 thru 480
// col: lowlimit, 10^15 thru 18
// rough tuning of where bitmap sieving is effective, as measured
// on a Intel Xeon CPU E5-2697 (AVX2 haswell) with one thread.
// numbers are the start index at which to start the bitmap sieve.
// Note: 480 classes seems to be too many to be efficient (all class
// lines must remain in memory during the bitmap sieve - this gets
// to be a lot for 480 classes!).  The same is true to a lesser
// extent for 48... but above 10^18 it is still a win.
// 
// the tuning is different on a 5122 Gold SkylakeX processor...
// likely it depends on cache size/speed.
#if 1
uint32 bitmap_bound_tab[4][5] = {
//	/*        10^14,     10^15,     10^16,     10^17,     10^18*/
//	/*2*/   { 200000,    700000,    500000,    700000,    200000  },
//	/*8*/   { 999999999, 700000,    700000,    700000,    700000 },
//	/*48*/  { 999999999, 999999999, 999999999, 999999999, 5000000   },
//	/*480*/ { 999999999, 999999999, 999999999, 999999999, 999999999 } };
#if 1
/*            10^14,     10^15,     10^16,     10^17,     10^18*/
/*2*/       { 200000,    700000,    500000,    700000,    500000  },
/*8*/       { 999999999, 999999999, 999999999, 999999999, 700000 },
/*48*/      { 999999999, 999999999, 999999999, 999999999, 999999999 },
/*480*/     { 999999999, 999999999, 999999999, 999999999, 999999999 } };
#else
/*            10^14,     10^15,     10^16,     10^17,     10^18*/
/*2*/   { 200000,    700000,    500000,    700000,    200000  },
/*8*/   { 999999999, 700000,    700000,    700000,    700000 },
/*48*/  { 999999999, 999999999, 999999999, 999999999, 5000000 },
/*480*/ { 999999999, 999999999, 999999999, 999999999, 999999999 } };
#endif
#else
uint32 bitmap_bound_tab[4][5] = {
	/*        10^14,     10^15,     10^16,     10^17,     10^18*/
	/*2*/   { 999999999, 999999999, 999999999, 999999999, 999999999 },
	/*8*/   { 999999999, 999999999, 999999999, 999999999, 999999999 },
	/*48*/  { 999999999, 999999999, 999999999, 999999999, 999999999 },
	/*480*/ { 999999999, 999999999, 999999999, 999999999, 999999999 } };
#endif



void get_numclasses(uint64 highlimit, uint64 lowlimit, soe_staticdata_t *sdata)
{
	uint64 numclasses, prodN, startprime;

    sdata->use_monty = 0;

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

	sieve_line_ptr = &sieve_line;
    FLAGBITS = 18;
    BUCKETSTARTI = 33336;

	//SOEBLOCKSIZE = 524288;
	FLAGSIZE = 8 * SOEBLOCKSIZE;
	FLAGSIZEm1 = FLAGSIZE - 1;

    switch (SOEBLOCKSIZE)
    {
    case 32768:
        FLAGBITS = 18;
        BUCKETSTARTI = 33336;

        // the avx2 version is faster on avx512 capable cpus...
#ifdef USE_AVX512Fa
        sieve_line_ptr = &sieve_line_avx512_32k;
#elif defined(USE_AVX2)
        sieve_line_ptr = &sieve_line_avx2_32k;
#endif
        break;
    case 65536:
        FLAGBITS = 19;
        BUCKETSTARTI = 43392;
        break;
    case 131072:
        FLAGBITS = 20;
        BUCKETSTARTI = 123040;
#ifdef USE_AVX512F
        sieve_line_ptr = &sieve_line_avx512_128k;
#elif defined(USE_AVX2)
        sieve_line_ptr = &sieve_line_avx2_128k;
#endif
        break;
    case 262144:
#ifdef USE_AVX512F
        sieve_line_ptr = &sieve_line_avx512_256k;
#elif defined(USE_AVX2)

#endif
        FLAGBITS = 21;
        BUCKETSTARTI = 233416;
        break;
    case 524288:
#ifdef USE_AVX512F
        sieve_line_ptr = &sieve_line_avx512_512k;
#elif defined(USE_AVX2)
        // the non-avx2 sieve is better, at least,
        // for huge offsets when you might be using this blocksize.
        //sieve_line_ptr = &sieve_line_avx2_512k;
#endif
        FLAGBITS = 22;
        BUCKETSTARTI = 443920;
        break;
    case 1048576:
        FLAGBITS = 23;
        BUCKETSTARTI = 846248;
        break;
    default:
        printf("Bad soe_block\n");
        exit(1);
    }

	//printf("Sieve Parameters:\nBLOCKSIZE = %u\nFLAGSIZE = %u\nFLAGBITS = %u\nBUCKETSTARTI = %u\n",
	//	SOEBLOCKSIZE, FLAGSIZE, FLAGBITS, BUCKETSTARTI);

	//more efficient to sieve using mod210 when the range is big
	if ((highlimit - lowlimit) > 40000000000ULL)
	{
        if (lowlimit < 100000000000000ULL)
        {
            numclasses = 480;
            prodN = 2310;
            startprime = 5;
        }
        else
        {
            numclasses = 48;
            prodN = 210;
            startprime = 4;
        }
#if defined(USE_AVX2)
        sdata->use_monty = 1;
#endif
	}	
	else if ((highlimit - lowlimit) > 4000000000ULL)
	{
        numclasses = 48;
        prodN = 210;
        startprime = 4;
#if defined(USE_AVX2)
        sdata->use_monty = 1;
#endif
	}
	else if ((highlimit - lowlimit) > 100000000)
	{
		numclasses=8;
		prodN=30;
		startprime=3;
#if defined(USE_AVX2)
        sdata->use_monty = 1;
#endif
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

	// the wrapper should handle this, but just in case we are called
	// directly and not via the wrapper...
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

	// set sieve primes in the local data structure to the ones that were passed in
	sdata->sieve_p = sieve_p;	
	
	if (offset == NULL)
	{
		//see if we were provided enough primes to do the job
		sdata->pbound = (uint64)(sqrt((int64)(highlimit)));

		if (sieve_p[num_sp - 1] < sdata->pbound)
		{
			printf("found %d primes, max = %u: not enough sieving primes\n", 
                num_sp, sieve_p[num_sp - 1]);
			exit(1);
		}

		// find the highest index that we'll need.  Much of the rest of the code is 
		// sensitive to this.  Note that this could be slow for large numbers of
		// sieve primes... could replace with a binary search.
		for (i=0; i<num_sp; i++)
		{
			// stop when we have enough for this input
            if (sieve_p[i] > sdata->pbound)
            {
                break;
            }
		}
		sdata->pboundi = i;	

#ifdef USE_AVX2
        // plus perhaps a few extra to get us to a convienient vector boundary
        while (sdata->pboundi & 7) sdata->pboundi++;
#endif
		
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
    int i, j, k;
    uint64 numclasses = sdata->numclasses;
    uint64 prodN = sdata->prodN;
    uint64 allocated_bytes = 0;
    uint64 lowlimit = sdata->orig_llimit;
    uint64 highlimit = sdata->orig_hlimit;
    uint64 numflags, numbytes, numlinebytes;


    // some groupings of residue classes and blocks per thread
    // 240 = 48 * 5
    // 240 = 8 * 30
    // 240 = 16 * 15
    // 256 = 8 * 32
    // 256 = 16 * 16
    // 272 = 16 * 17
    // 272 = 8 * 34
    

    // allocate the residue classes.  
    sdata->rclass = (uint32 *)malloc(numclasses * sizeof(uint32));
    allocated_bytes += numclasses * sizeof(uint32);

    // find the residue classes
    k = 0;
    for (i = 1; i < prodN; i++)
    {
        if (spGCD(i, (uint64)prodN) == 1)
        {
            sdata->rclass[k] = (uint32)i;
            k++;
        }
    }

    sdata->min_sieved_val = 1ULL << 63;

    // temporarily set lowlimit to the first multiple of numclasses*prodN < lowlimit
    if (sdata->sieve_range == 0)
    {
        lowlimit = (lowlimit / (numclasses*prodN))*(numclasses*prodN);
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

    // reallocate flag structure for wheel and block sieving.
    // we want the smallest number of blocks such that all integers in the requested
    // range are represented.  Each block contains 32768 * 8 bit-flags, each of
    // which represents integers spaced 'prodN' apart.
    sdata->blocks = (highlimit - lowlimit) / prodN / FLAGSIZE;
    if (((highlimit - lowlimit) / prodN) % FLAGSIZE != 0) sdata->blocks++;
    sdata->numlinebytes = sdata->blocks * SOEBLOCKSIZE;
    numlinebytes = sdata->numlinebytes;
    highlimit = (uint64)((uint64)sdata->numlinebytes * (uint64)prodN * (uint64)BITSINBYTE + lowlimit);
    sdata->highlimit = highlimit;

    // each flag in a block is spaced prodN integers apart.  record the resulting size of the 
    // number line encoded in each block.
    sdata->blk_r = FLAGSIZE*prodN;

    // compute the breakpoints at which we switch to other sieving methods	
    sdata->num_bitmap_primes = 0;
    sdata->bitmap_start_id = sdata->pboundi;
    if (sdata->pboundi > BUCKETSTARTI)
    {
        sdata->bucket_start_id = BUCKETSTARTI;

        // also see if a bitmap sieving step will be beneficial.
        sdata->bitmap_lower_bound = 99999999999999ULL;
        if (((int)log10(sdata->orig_llimit) - 14) >= 0)
        {
            sdata->bitmap_start_id = MIN(
                bitmap_bound_tab[sdata->startprime - 2][(int)log10(sdata->orig_llimit) - 14],
                sdata->pboundi);

            if ((sdata->bitmap_start_id < sdata->pboundi))
                sdata->num_bitmap_primes = sdata->pboundi - sdata->bitmap_start_id;
            else
                sdata->bitmap_start_id = sdata->pboundi;
        }

        // buckets operate up to the bitmap bound (which may be the 
        // equal to the sieve bound, if not benficial).
        sdata->num_bucket_primes = sdata->bitmap_start_id - sdata->bucket_start_id;
    }
    else
    {
        sdata->num_bucket_primes = 0;
        sdata->bucket_start_id = sdata->pboundi;
    }

    // any prime larger than this will only hit the interval once (in residue space)
    sdata->large_bucket_start_prime = sdata->blocks * FLAGSIZE;

    // allocate space for the root of each sieve prime (used by the bucket sieve)
    sdata->root = (int *)malloc(sdata->bitmap_start_id * sizeof(int));
    allocated_bytes += sdata->bitmap_start_id * sizeof(uint32);
    if (sdata->root == NULL)
    {
        printf("error allocating roots\n");
        exit(-1);
    }
    else
    {
        if (VFLAG > 2)
            printf("allocated %u bytes for roots\n", 
            (uint32)(sdata->bitmap_start_id * sizeof(uint32)));
    }

    sdata->r2modp = (uint32 *)xmalloc_align(sdata->bitmap_start_id * sizeof(uint32));
    allocated_bytes += sdata->bitmap_start_id * sizeof(uint32);
    if (sdata->r2modp == NULL)
    {
        printf("error allocating r2modp\n");
        exit(-1);
    }
    else
    {
        if (VFLAG > 2)
            printf("allocated %u bytes for r2modp\n", 
            (uint32)(sdata->bitmap_start_id * sizeof(uint32)));
    }

    // experimental montgomery arithmetic
    sdata->pinv = (uint32 *)xmalloc_align(sdata->bitmap_start_id * sizeof(uint32));
    allocated_bytes += sdata->bitmap_start_id * sizeof(uint32);
    if (sdata->pinv == NULL)
    {
        printf("error allocating pinv\n");
        exit(-1);
    }
    else
    {
        if (VFLAG > 2)
            printf("allocated %u bytes for pinv\n", 
            (uint32)(sdata->bitmap_start_id * sizeof(uint32)));
    }


    // these are used by the bucket sieve
    sdata->lower_mod_prime = (uint32 *)malloc(sdata->bitmap_start_id * sizeof(uint32));
    allocated_bytes += sdata->bitmap_start_id * sizeof(uint32);
    if (sdata->lower_mod_prime == NULL)
    {
        printf("error allocating lower mod prime\n");
        exit(-1);
    }
    else
    {
        if (VFLAG > 2)
            printf("allocated %u bytes for lower mod prime\n",
            (uint32)sdata->bitmap_start_id * (uint32)sizeof(uint32));
    }

    // allocate all of the lines if we are computing primes or if we
    // are using bitmap sieving.  otherwise lines will
    // be allocated as needed during sieving
    sdata->lines = (uint8 **)xmalloc_align(sdata->numclasses * sizeof(uint8 *));
    numbytes = 0;
    
    if ((sdata->only_count == 0) || (sdata->num_bitmap_primes > 0))
    {
        //actually allocate all of the lines as a continuous linear array of bytes
        sdata->lines[0] = (uint8 *)xmalloc_align(numlinebytes * sdata->numclasses * sizeof(uint8));
        if (sdata->lines[0] == NULL)
        {
            printf("error allocated sieve lines\n");
            exit(-1);
        }
        numbytes += sdata->numclasses * numlinebytes * sizeof(uint8);

        for (i = 0; i < sdata->numclasses; i++)
        {
            sdata->lines[i] = sdata->lines[0] + i * numlinebytes;
        }
    }
    else
    {
        //don't allocate anything now, but
        //provide an figure for the memory that will be allocated later
        numbytes = numlinebytes * sizeof(uint8) * THREADS;
    }


#ifdef USE_AVX2
    // during presieveing, storing precomputed lists will start to get unwieldy, so
    // generate the larger lists here.
#ifdef USE_AVX512Fa
#define DYNAMIC_BOUND 512
#else
#define DYNAMIC_BOUND 256
#endif

    sdata->presieve_max_id = 40;

    for (j = 24; j < sdata->presieve_max_id; j++)
    {
        uint32 prime = sdata->sieve_p[j];

        // for each possible starting location
        for (i = 0; i < prime; i++)
        {
            int x;
            uint64 interval[DYNAMIC_BOUND/64];

            for (x = 0; x < DYNAMIC_BOUND/64; x++)
                interval[x] = 0xffffffffffffffffULL;

            // sieve up to the bound, printing each 64-bit word as we fill it
            for (k = i; k < DYNAMIC_BOUND; k += prime)
            {
                interval[k >> 6] &= ~(1ULL << (k & 63));
            }

            //printf("largemask[%d][%d] = %016lx, %016lx, %016lx, %016lx", j-24, i,
            //    interval[0], interval[1], interval[2], interval[3]);
            presieve_largemasks[j - 24][i][0] = interval[0];
            presieve_largemasks[j - 24][i][1] = interval[1];
            presieve_largemasks[j - 24][i][2] = interval[2];
            presieve_largemasks[j - 24][i][3] = interval[3];
#ifdef USE_AVX512Fa
            //printf(", %016lx, %016lx, %016lx, %016lx", 
            //    interval[4], interval[5], interval[6], interval[7]);
            presieve_largemasks[j - 24][i][4] = interval[4];
            presieve_largemasks[j - 24][i][5] = interval[5];
            presieve_largemasks[j - 24][i][6] = interval[6];
            presieve_largemasks[j - 24][i][7] = interval[7];
#endif
            //printf("\n");

        }        
    }

    //printf("primes: ");
    for (j = 24; j < sdata->presieve_max_id; j++)
    {        
        presieve_primes[j - 24] = sdata->sieve_p[j];
        //printf("%u ", presieve_primes[j - 24]);
    }
    //printf("\n");

    for (j = 24; j < sdata->presieve_max_id; j++)
    {
        presieve_p1[j - 24] = sdata->sieve_p[j] - 1;
    }

    //printf("steps: ");
    for (j = 24; j < sdata->presieve_max_id; j++)
    {        
        presieve_steps[j - 24] = DYNAMIC_BOUND % sdata->sieve_p[j];
        //printf("%u ", presieve_steps[j - 24]);
    }
    //printf("\n");

#endif

#if defined(USE_BMI2) || defined(USE_AVX512F)
#ifdef __INTEL_COMPILER
    if (_may_i_use_cpu_feature(_FEATURE_BMI))
#elif defined(__GNUC__)
    if (__builtin_cpu_supports("bmi2"))
#else
    if (0)
#endif
    {
        compute_8_bytes_ptr = &compute_8_bytes_bmi2;
    }
    else
    {
        compute_8_bytes_ptr = &compute_8_bytes;
    }
#else
    compute_8_bytes_ptr = &compute_8_bytes;
#endif

#if defined(USE_AVX2)
#ifdef __INTEL_COMPILER
    if (_may_i_use_cpu_feature(_FEATURE_AVX2))
#elif defined(__GNUC__)
    if (__builtin_cpu_supports("avx2"))
#else
    if (1)
#endif
    {
        pre_sieve_ptr = &pre_sieve_avx2;
    }
    else
    {
        pre_sieve_ptr = &pre_sieve;
        sdata->presieve_max_id = 10;
    }

#ifdef USE_AVX512Fa
    pre_sieve_ptr = &pre_sieve_avx512;
#endif

#else
    // if we haven't built the code with AVX2 support, or if at runtime
    // we find that AVX2 isn't supported, use the portable version
    // of these routines.
    pre_sieve_ptr = &pre_sieve;
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

	if (sdata->num_bucket_primes > 0)
	{
		// then we have primes bigger than BUCKETSTARTP - need to bucket sieve
		uint64 flagsperline = numlinebytes * 8;
		uint64 num_hits = 0;
		uint64 hits_per_bucket;
		
		for (i = sdata->bucket_start_id; i < sdata->bucket_start_id + sdata->num_bucket_primes; i++)
		{
			// condition to see if the current prime only hits the sieve interval once
			if ((sdata->sieve_p[i] * sdata->prodN) > (sdata->blk_r * sdata->blocks))
				break;
			num_hits += ((uint32)flagsperline / sdata->sieve_p[i] + 1);
		}

		// assume hits are evenly distributed among buckets.
		hits_per_bucket = num_hits / sdata->blocks;

		// add some margin
		hits_per_bucket = (uint64)((double)hits_per_bucket * 1.10);

		// set the bucket allocation amount, with a minimum of at least 50000
		// because small allocation amounts may violate the uniformity assumption
		// of hits per bucket.  The idea is to set this right once, even if it is too big,
		// so that we don't have to keep checking for full buckets in the middle of
		// the bucket sieve (which would be slow)
		sdata->bucket_alloc = MAX(hits_per_bucket,50000);

		// now count primes that only hit the interval once
		num_hits = 0;
		for (; i < sdata->bucket_start_id + sdata->num_bucket_primes; i++)
			num_hits++;

		// assume hits are evenly distributed among buckets.
		hits_per_bucket = num_hits / sdata->blocks;

		// add some margin
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
        thread->ddata.presieve_scratch = (uint32 *)xmalloc_align(16 * sizeof(uint32));

		// allocate a bound for each block
		thread->ddata.pbounds = (uint64 *)malloc(
			sdata->blocks * sizeof(uint64));
		allocated_bytes += sdata->blocks * sizeof(uint64);

		thread->ddata.pbounds[0] = sdata->pboundi;

		// we'll need to store the offset into the next block for each prime.
		// actually only need those primes less than BUCKETSTARTP since bucket sieving
		// doesn't use the offset array.
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

