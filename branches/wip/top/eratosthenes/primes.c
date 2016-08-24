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
#include "threadpool.h"
#include <immintrin.h>

//for testing one of 8 bits in a byte in one of 8 lines.
//bit num picks the row, lines num picks the col.	
const uint64 nmasks64[8][8] = {
    { 1ULL, 256ULL, 65536ULL, 16777216ULL, 4294967296ULL, 1099511627776ULL, 281474976710656ULL, 72057594037927936ULL },
    { 2ULL, 512ULL, 131072ULL, 33554432ULL, 8589934592ULL, 2199023255552ULL, 562949953421312ULL, 144115188075855872ULL },
    { 4ULL, 1024ULL, 262144ULL, 67108864ULL, 17179869184ULL, 4398046511104ULL, 1125899906842624ULL, 288230376151711744ULL },
    { 8ULL, 2048ULL, 524288ULL, 134217728ULL, 34359738368ULL, 8796093022208ULL, 2251799813685248ULL, 576460752303423488ULL },
    { 16ULL, 4096ULL, 1048576ULL, 268435456ULL, 68719476736ULL, 17592186044416ULL, 4503599627370496ULL, 1152921504606846976ULL },
    { 32ULL, 8192ULL, 2097152ULL, 536870912ULL, 137438953472ULL, 35184372088832ULL, 9007199254740992ULL, 2305843009213693952ULL },
    { 64ULL, 16384ULL, 4194304ULL, 1073741824ULL, 274877906944ULL, 70368744177664ULL, 18014398509481984ULL, 4611686018427387904ULL },
    { 128ULL, 32768ULL, 8388608ULL, 2147483648ULL, 549755813888ULL, 140737488355328ULL, 36028797018963968ULL, 9223372036854775808ULL } };


void compute_primes_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *t = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = t->sdata;

    // launch one range of computation for each thread.  don't really
    // need a threadpool for this, but the infrastructure is there...
    if (sdata->sync_count < THREADS)
    {
        tdata->work_fcn_id = 0;
        sdata->sync_count++;
    }
    else
    {
        tdata->work_fcn_id = tdata->num_work_fcn;
    }

    return;
}

void compute_primes_work_fcn(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *udata = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];
    int i;

    if (THREADS > 1)
    {
        t->linecount = 0;
    }
    
#ifdef USE_AVX2
    if ((sdata->numclasses == 2) && (sdata->lowlimit == 0))
    {
        //printf("using avx2 nc2 code from byte offset %d to %d in steps of 32 bytes\n",
          //  t->startid, t->stopid);
        for (i = t->startid; i < t->stopid; i += 32)
        {
            t->linecount = compute_32_bytes(sdata, t->linecount, t->ddata.primes, i);
        }
    }
    else
    {
        for (i = t->startid; i < t->stopid; i += 8)
        {
            t->linecount = compute_8_bytes(sdata, t->linecount, t->ddata.primes, i);
        }
    }
#else
    for (i = t->startid; i < t->stopid; i += 8)
    {
        t->linecount = compute_8_bytes(sdata, t->linecount, t->ddata.primes, i);
    }
#endif

    return;
}


uint64 primes_from_lineflags(soe_staticdata_t *sdata, thread_soedata_t *thread_data,
	uint32 start_count, uint64 *primes)
{
	//compute primes using all of the sieved lines we have stored
	uint32 pcount = start_count;	
	uint64 i;
	int j;
	uint32 range, lastid;

    //timing
    double t;
    struct timeval tstart, tstop;
    TIME_DIFF *	difference;

    // threading structures
    tpool_t *tpool_data;
    soe_userdata_t udata;

    if (VFLAG > 1)
    {
        gettimeofday(&tstart, NULL);
    }

	// each thread needs to work on a number of bytes that is divisible by 32
	range = sdata->numlinebytes / THREADS;
	range -= (range % 32);
	lastid = 0;

    // divvy up the line bytes
    for (i = 0; i < THREADS; i++)
    {
        thread_soedata_t *t = thread_data + i;

        t->sdata = *sdata;
        t->startid = lastid;
        t->stopid = t->startid + range;
        lastid = t->stopid;

        if (VFLAG > 2)
        {
            printf("thread %d finding primes from byte offset %u to %u\n",
                (int)i, t->startid, t->stopid);
        }
    }

    // allocate a temporary array for each thread's primes
    if (THREADS > 1)
    {
        uint64 memchunk;

        if (sdata->sieve_range)
        {
            // then just split the overall range into equal parts
            memchunk = (sdata->orig_hlimit - sdata->orig_llimit) / THREADS + THREADS;

            for (i = 0; i < THREADS; i++)
            {
                thread_soedata_t *t = thread_data + i;

                t->ddata.primes = (uint64 *)malloc(memchunk * sizeof(uint64));
            }
        }
        else
        {
            // then estimate the number of primes we'll find in each chunk.
            // it's important to do this chunk by chunk, because in some cases
            // the number of primes changes rapidly as a function of offset
            // from the start of the range (i.e., when start is 0)
            uint64 hi_est, lo_est;
            uint64 tmplo = sdata->orig_llimit;
            uint64 tmphi;
            uint64 chunk = 8 * sdata->numlinebytes / THREADS;

            chunk *= sdata->prodN;
            tmphi = tmplo + chunk;
            for (i = 0; i < THREADS; i++)
            {
                thread_soedata_t *t = thread_data + i;

                hi_est = (uint64)(tmphi / log((double)tmphi));
                if (tmplo > 1)
                    lo_est = (uint64)(tmplo / log((double)tmplo));
                else
                    lo_est = 0;

                memchunk = (uint64)((double)(hi_est - lo_est) * 1.25);

                if (VFLAG > 2)
                {
                    printf("allocating temporary space for %" PRIu64 " primes between %" PRIu64 " and %" PRIu64 "\n",
                        memchunk, tmplo, tmphi);
                }

                t->ddata.primes = (uint64 *)malloc(memchunk * sizeof(uint64));

                tmplo += chunk;
                tmphi += chunk;
            }
        }
    }
    else
    {
        // with just one thread, don't bother with creating a temporary array
        thread_data[0].ddata.primes = primes;
    }

    udata.sdata = sdata;
    udata.ddata = thread_data;
    tpool_data = tpool_setup(THREADS, NULL, NULL, NULL,
        &compute_primes_dispatch, &udata);

    if (THREADS == 1)
    {
        thread_data->linecount = pcount;
        compute_primes_work_fcn(tpool_data);
    }
    else
    {
        sdata->sync_count = 0;
        tpool_add_work_fcn(tpool_data, &compute_primes_work_fcn);
        tpool_go(tpool_data);
    }
    free(tpool_data);    

	// now combine all of the temporary arrays, if necessary
	if (THREADS > 1)
	{
		pcount = start_count;
		for (j = 0; j < THREADS; j++)
		{
			thread_soedata_t *t = thread_data + j;

			if (t->linecount == 0)
			{
				free(t->ddata.primes);
				continue;
			}

			if (VFLAG > 2)
				printf("adding %" PRIu64 " primes found in thread %d\n", t->linecount, j);

			memcpy(primes + GLOBAL_OFFSET + pcount, t->ddata.primes, t->linecount * sizeof(uint64));

			pcount += t->linecount;
			free(t->ddata.primes);
		}
	}
    else
    {
        pcount = thread_data[0].linecount;
    }

	// and finally, get primes from any residual portion of the line arrays
	// using a direct method
	if (lastid != sdata->numlinebytes)
	{
		if (VFLAG > 2)
			printf("adding primes from byte offset %u to %u\n", 
				lastid, (uint32)sdata->numlinebytes);

		for (i = lastid; i < sdata->numlinebytes; i+=8)
		{
			pcount = compute_8_bytes(sdata, pcount, primes, i);		
		}
	}

    if (VFLAG > 1)
    {
        gettimeofday(&tstop, NULL);

        difference = my_difftime(&tstart, &tstop);
        t = ((double)difference->secs + (double)difference->usecs / 1000000);
        free(difference);

        if (VFLAG > 2)
        {
            printf("time to compute primes = %1.2f\n", t);
        }
    }

	return pcount;
}

uint32 compute_32_bytes(soe_staticdata_t *sdata,
    uint32 pcount, uint64 *primes, uint64 byte_offset)
{
    int b;
    uint64 prime;
    uint64 prodN = sdata->prodN;
    uint8 **lines = sdata->lines;
    uint64 ohigh = sdata->orig_hlimit;

    if ((byte_offset & 32767) == 0)
    {
        if (VFLAG > 1)
        {
            printf("computing: %d%%\r", (int)
                ((double)byte_offset / (double)(sdata->numlinebytes) * 100.0));
            fflush(stdout);
        }
    }

#ifdef USE_AVX2

    // here is the 2 line version
    if (1)
    {
        // here is the process for 2 lines.
        // load the first 8 dwords from line 1 and the first 8 dwords from line 2.
        // vreg1 = {l1,0, l1,1, l1,2, ... l1,7}
        // vreg2 = {l2,0, l2,1, l2,2, ... l2,7}
        // where li,j is the ith line, jth dword.
        // the bits are again ordered columnwise, so we do parallel bitcount
        // to determine where each dword computation should write the prime array.
        // i.e., if l1,0 has 9 set bits (of 32) and l2,0 has 5 set bits then
        // primes emerging from l1,1 are written to the primes array offset by 9+5 locations.
        // other than the larger number of offsets to keep track of compared to
        // the 8-line case, everything proceeds similarly.  However, a benefit is we
        // don't need to use vgather... just two vector loads.
        // Also we compute 8*4 = 32 bytes of each line instead of 8, so we'll
        // need fewer iterations of this function.
        // need a vector of the indices of the memory elements to load,
        // relative to the base address.  we set the base address to the
        // first line's byte offset, so the indices are all just multiples
        // of the length of a line.
        __m256i vtmp1;
        __m256i vlinenums;
        __m256i vlinenums2;
        __m256i vmasks = _mm256_set1_epi32(0x1);
        __m256i vprodN = _mm256_set1_epi32((uint32)prodN);
        __m256i vbitoffset;
        __m256i vbits, vbits2;
        __m256i vtmp2, vtmp3, vtmp4;
#if defined(__GNUC__)
        __attribute__((aligned(64))) uint32 tmpstore[8];
        __attribute__((aligned(64))) uint32 tmpstore2[8];
#else
        __declspec(align(64)) uint32 tmpstore[8];
        __declspec(align(64)) uint32 tmpstore2[8];
#endif
        int i;
        int poffset1[8];
        //int poffset2[8];
        int pcounts1[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        //int pcounts2[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        uint32 bmask, bmask2;

        // all primes from the first vector are in class 1 and
        // all primes from the second vector are in class 5.
        vlinenums = _mm256_set1_epi32(1);
        vlinenums2 = _mm256_set1_epi32(5);

        // load 256 bits of each line
        vbits = _mm256_load_si256((__m256i *)(&lines[0][byte_offset]));
        vtmp4 = _mm256_load_si256((__m256i *)(&lines[0][byte_offset + sdata->numlinebytes]));

        // count the primes in each vector
        {
            __m256i v5 = _mm256_set1_epi32(0x55555555);
            __m256i v3 = _mm256_set1_epi32(0x33333333);
            __m256i v0f = _mm256_set1_epi32(0x0F0F0F0F);
            __m256i v3f = _mm256_set1_epi32(0x0000003F);
            vbits2 = vbits;
            vtmp1 = _mm256_srli_epi32(vbits2, 1);
            vtmp1 = _mm256_and_si256(vtmp1, v5);
            vbits2 = _mm256_sub_epi32(vbits2, vtmp1);
            vtmp1 = _mm256_and_si256(vbits2, v3);
            vtmp2 = _mm256_srli_epi32(vbits2, 2);
            vtmp2 = _mm256_and_si256(vtmp2, v3);
            vbits2 = _mm256_add_epi32(vtmp2, vtmp1);
            vtmp1 = _mm256_srli_epi32(vbits2, 4);
            vbits2 = _mm256_add_epi32(vbits2, vtmp1);
            vbits2 = _mm256_and_si256(vbits2, v0f);
            vtmp1 = _mm256_srli_epi32(vbits2, 8);
            vbits2 = _mm256_add_epi32(vbits2, vtmp1);
            vtmp1 = _mm256_srli_epi32(vbits2, 16);
            vbits2 = _mm256_add_epi32(vbits2, vtmp1);
            vbits2 = _mm256_and_si256(vbits2, v3f);
            _mm256_store_si256((__m256i *)tmpstore, vbits2);
        }

        // count the primes in each vector
        {
            __m256i v5 = _mm256_set1_epi32(0x55555555);
            __m256i v3 = _mm256_set1_epi32(0x33333333);
            __m256i v0f = _mm256_set1_epi32(0x0F0F0F0F);
            __m256i v3f = _mm256_set1_epi32(0x0000003F);
            vbits2 = vtmp4;
            vtmp1 = _mm256_srli_epi32(vbits2, 1);
            vtmp1 = _mm256_and_si256(vtmp1, v5);
            vbits2 = _mm256_sub_epi32(vbits2, vtmp1);
            vtmp1 = _mm256_and_si256(vbits2, v3);
            vtmp2 = _mm256_srli_epi32(vbits2, 2);
            vtmp2 = _mm256_and_si256(vtmp2, v3);
            vbits2 = _mm256_add_epi32(vtmp2, vtmp1);
            vtmp1 = _mm256_srli_epi32(vbits2, 4);
            vbits2 = _mm256_add_epi32(vbits2, vtmp1);
            vbits2 = _mm256_and_si256(vbits2, v0f);
            vtmp1 = _mm256_srli_epi32(vbits2, 8);
            vbits2 = _mm256_add_epi32(vbits2, vtmp1);
            vtmp1 = _mm256_srli_epi32(vbits2, 16);
            vbits2 = _mm256_add_epi32(vbits2, vtmp1);
            vbits2 = _mm256_and_si256(vbits2, v3f);
            _mm256_store_si256((__m256i *)tmpstore2, vbits2);
        }

        // map each dword to the correct location in the primes array
        // using the counts we just computed.
        poffset1[0] = pcount;
        poffset1[1] = poffset1[0] + tmpstore[0] + tmpstore2[0];
        poffset1[2] = poffset1[1] + tmpstore[1] + tmpstore2[1];
        poffset1[3] = poffset1[2] + tmpstore[2] + tmpstore2[2];
        poffset1[4] = poffset1[3] + tmpstore[3] + tmpstore2[3];
        poffset1[5] = poffset1[4] + tmpstore[4] + tmpstore2[4];
        poffset1[6] = poffset1[5] + tmpstore[5] + tmpstore2[5];
        poffset1[7] = poffset1[6] + tmpstore[6] + tmpstore2[6];

#ifdef PRINT_DEBUG
        {
            printf("at byte offset %d, here are the dword prime counts with initial count %u\n", 
                byte_offset, pcount);
            printf("class 1:\n");
            for (i = 0; i < 32; i++)
            {
                printf("%02x", lines[0][byte_offset + i]);
                if (i % 4 == 3) printf(" ");
            }
            printf("\nclass 5:\n");
            for (i = 0; i < 32; i++)
            {
                printf("%02x", lines[0][byte_offset + sdata->numlinebytes + i]);
                if (i % 4 == 3) printf(" ");
            }
            printf("\ncounts1:\n");
            for (i = 0; i < 8; i++)
            {
                printf("%d ", tmpstore[i]);
            }
            printf("\ncounts2:\n");
            for (i = 0; i < 8; i++)
            {
                printf("%d ", tmpstore2[i]);
            }
            printf("\n");
            printf("\ncumulative offsets1:\n");
            for (i = 0; i < 8; i++)
            {
                printf("%d ", poffset1[i]);
            }
            printf("\ncumulative offsets2:\n");
            for (i = 0; i < 8; i++)
            {
                printf("%d ", poffset2[i]);
            }
            printf("\n");            
        }
#endif


        // each dword has a different bit offset.  they are the same
        // for each line (which have different class numbers).
        b = (uint32)byte_offset << 3;
        vbitoffset = _mm256_setr_epi32(b, b + 32, b + 64, b + 96,
            b + 128, b + 160, b + 192, b + 224);

        vbits2 = vtmp4;

        // for each of the 32 bits, compute a prime if the bit is set.
        for (i = 0; i < 32; i++)
        {
            //int id;
            int j;

            // isolate the current bit and compare if equal to 1.            
            vtmp2 = _mm256_cmpeq_epi32(_mm256_and_si256(vmasks, vbits), vmasks);
            vtmp4 = _mm256_cmpeq_epi32(_mm256_and_si256(vmasks, vbits2), vmasks);

            // extract topmost bit of every byte of the comparison mask
            bmask = _mm256_movemask_epi8(vtmp2);
            bmask2 = _mm256_movemask_epi8(vtmp4);

            // next bitmap location
            vbits = _mm256_srli_epi32(vbits, 1);
            vbits2 = _mm256_srli_epi32(vbits2, 1);

            // skip if all non-prime
            if ((bmask == 0) && (bmask2 == 0))
                continue;

            // compute the primes.  The only part of the calculation 
            // that involves more than 32 bits is the addition of
            // lowlimit.  This is because when we compute primes the
            // range is not allowed to exceed 32-bits in size.
            vtmp1 = _mm256_set1_epi32(i);
            vtmp1 = _mm256_add_epi32(vtmp1, vbitoffset);
            vtmp1 = _mm256_mullo_epi32(vtmp1, vprodN);
            vtmp3 = _mm256_add_epi32(vtmp1, vlinenums2);
            vtmp1 = _mm256_add_epi32(vtmp1, vlinenums);
            

            // Use the result to conditionally store the primes we computed.
            vtmp2 = _mm256_and_si256(vtmp2, vtmp1);
            vtmp4 = _mm256_and_si256(vtmp4, vtmp3);

            // one vector store.  (extracting directly from the vector is not profitable...)
            _mm256_store_si256((__m256i *)tmpstore, vtmp2);
            _mm256_store_si256((__m256i *)tmpstore2, vtmp4);

            // search the masked primes for non-zero values and copy
            // them to the output array.
            // the two vectors are contiguous along columns of dwords, so we
            // look at the dwords one at a time.  could probably also do a
            // similar tzcnt/blsr thing as the 8-class case if we first
            // OR together the bmasks...
            for (j = 0; j < 8; j++)
            {
                if (bmask & 0x1)
                {
                    prime = tmpstore[j];
                    if (prime <= ohigh)
                    {
                        primes[GLOBAL_OFFSET + poffset1[j] + pcounts1[j]] = prime;
                        pcounts1[j]++;
                    }
                }

                if (bmask2 & 0x1)
                {
                    prime = tmpstore2[j];
                    if (prime <= ohigh)
                    {
                        primes[GLOBAL_OFFSET + poffset1[j] + pcounts1[j]] = prime;
                        pcounts1[j]++;
                    }
                }

                bmask >>= 4;
                bmask2 >>= 4;
            }

        }

        for (i = 0; i < 8; i++)
        {
            pcount += pcounts1[i];
        }

#ifdef PRINT_DEBUG
        printf("pcount is now %d\n", pcount);
#endif
    }

#ifdef PRINT_DEBUG
    exit(1);
#endif
#endif

    return pcount;
}

uint32 compute_8_bytes(soe_staticdata_t *sdata, 
	uint32 pcount, uint64 *primes, uint64 byte_offset)
{
	int b;	
	uint32 current_line;
	//8 bytes from each of up to 48 sieve lines are packed into these words
	uint64 cache_word[64];	
	uint64 prime;
	uint32 nc = sdata->numclasses;
	uint32 *rclass = sdata->rclass;
	uint64 lowlimit = sdata->lowlimit;
	uint64 prodN = sdata->prodN;
	uint8 **lines = sdata->lines;
	uint64 olow = sdata->orig_llimit;
	uint64 ohigh = sdata->orig_hlimit;
		
	if ((byte_offset & 32767) == 0)
	{
		if (VFLAG > 1)
		{
			printf("computing: %d%%\r",(int)
				((double)byte_offset / (double)(sdata->numlinebytes) * 100.0));
			fflush(stdout);
		}
	}

    // AVX2 version:
    // vgatherdd from the i'th column of 8 lines into one vector.
    // in a loop of 32, AND each element with a mask of (1 << i).
    // compute 8 primes (need aux vectors of 3, prodN, and rclass 
    // for each loaded line.
    // continue until all lines are finished, then repeat for column
    // i + 1 of 32-bit entries.

#ifdef USE_AVX2

    // here is the 8 line version
    if ((nc == 8) && (lowlimit == 0))
    {
        // need a vector of the indices of the memory elements to load,
        // relative to the base address.  we set the base address to the
        // first line's byte offset, so the indices are all just multiples
        // of the length of a line.
        __m256i vtmp1 = _mm256_set1_epi32((uint32)sdata->numlinebytes);
        __m256i vlinenums = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
        __m256i voffsets = _mm256_mullo_epi32(vtmp1, vlinenums);
        __m256i vmasks = _mm256_set1_epi32(0x1);
        __m256i vprodN = _mm256_set1_epi32((uint32)prodN);
        __m256i vbitoffset = _mm256_set1_epi32((uint32)byte_offset << 3);
        //__m256i vlowlimit = _mm256_set1_epi64x(lowlimit);
        __m256i vbits, vbits2;
        __m256i vtmp2, vtmp3, vtmp4;
        __m256i vbitoffset2;
#if defined(__GNUC__)
        __attribute__((aligned(64))) uint32 tmpstore[8];
#else
        __declspec(align(64)) uint32 tmpstore[8];
#endif
        int *startaddr = (int *)(&lines[0][byte_offset]);
        int i;
        int poffset = pcount;
        uint32 pcount2 = 0;
        uint32 bmask, bmask2;

        // reuse linenums now that we've computed the offsets for each line
        vlinenums = _mm256_setr_epi32(1, 7, 11, 13, 17, 19, 23, 29);

        // gather the first 32-bits of each of our 8 lines in one vector
        vbits = _mm256_i32gather_epi32(startaddr, voffsets, 1);

        // count the primes in this vector so we know where to start adding
        // primes from the next vector
        {
            __m256i v5 = _mm256_set1_epi32(0x55555555);
            __m256i v3 = _mm256_set1_epi32(0x33333333);
            __m256i v0f = _mm256_set1_epi32(0x0F0F0F0F);
            __m256i v3f = _mm256_set1_epi32(0x0000003F);
            vbits2 = vbits;
            vtmp1 = _mm256_srli_epi64(vbits2, 1);
            vtmp1 = _mm256_and_si256(vtmp1, v5);
            vbits2 = _mm256_sub_epi64(vbits2, vtmp1);
            vtmp1 = _mm256_and_si256(vbits2, v3);
            vtmp2 = _mm256_srli_epi64(vbits2, 2);
            vtmp2 = _mm256_and_si256(vtmp2, v3);
            vbits2 = _mm256_add_epi64(vtmp2, vtmp1);
            vtmp1 = _mm256_srli_epi64(vbits2, 4);
            vbits2 = _mm256_add_epi64(vbits2, vtmp1);
            vbits2 = _mm256_and_si256(vbits2, v0f);
            vtmp1 = _mm256_srli_epi64(vbits2, 8);
            vbits2 = _mm256_add_epi64(vbits2, vtmp1);
            vtmp1 = _mm256_srli_epi64(vbits2, 16);
            vbits2 = _mm256_add_epi64(vbits2, vtmp1);
            vtmp1 = _mm256_srli_epi64(vbits2, 32);
            vbits2 = _mm256_add_epi64(vbits2, vtmp1);
            vbits2 = _mm256_and_si256(vbits2, v3f);
            _mm256_store_si256((__m256i *)tmpstore, vbits2);
            poffset += tmpstore[0] + tmpstore[2] + tmpstore[4] + tmpstore[6];
        }

        startaddr = (int *)(&lines[0][byte_offset + 4]);
        vbits2 = _mm256_i32gather_epi32(startaddr, voffsets, 1);
        vbitoffset2 = _mm256_set1_epi32((byte_offset + 4) << 3);     

        // for each of the 32 bits, compute a prime if the bit is set.
        // we operate on each line's 32-bit block simultaneously and via
        // the gather operation, the output primes will be contiguous.
        for (i = 0; i < 32; i++)
        {            
            int id;

            // isolate the current bit and compare if equal to 1.            
            vtmp2 = _mm256_cmpeq_epi32(_mm256_and_si256(vmasks, vbits), vmasks);
            vtmp4 = _mm256_cmpeq_epi32(_mm256_and_si256(vmasks, vbits2), vmasks);

            // extract topmost bit of every byte of the comparison mask
            bmask = _mm256_movemask_epi8(vtmp2);
            bmask2 = _mm256_movemask_epi8(vtmp4);

            // next bitmap location
            vbits = _mm256_srli_epi32(vbits, 1);
            vbits2 = _mm256_srli_epi32(vbits2, 1);

            // skip if all non-prime
            if ((bmask == 0) && (bmask2 == 0))
                continue;

            // compute the primes.  The only part of the calculation 
            // that involves more than 32 bits is the addition of
            // lowlimit.  This is because when we compute primes the
            // range is not allowed to exceed 32-bits in size.
            vtmp1 = _mm256_set1_epi32(i);
            vtmp3 = _mm256_add_epi32(vtmp1, vbitoffset2);
            vtmp1 = _mm256_add_epi32(vtmp1, vbitoffset);            
            vtmp1 = _mm256_mullo_epi32(vtmp1, vprodN);
            vtmp3 = _mm256_mullo_epi32(vtmp3, vprodN);
            vtmp1 = _mm256_add_epi32(vtmp1, vlinenums);
            vtmp3 = _mm256_add_epi32(vtmp3, vlinenums);

            // now we have to split the 8x vector into two 4x vectors
            // and add in the lowlimit.  for now assume lowlimit is zero...
            // should create more than one case to handle this.



            // Use the result to conditionally store the primes we computed.
            vtmp2 = _mm256_and_si256(vtmp2, vtmp1);
            vtmp4 = _mm256_and_si256(vtmp4, vtmp3);

            // one vector store.  (extracting directly from the vector is not profitable...)
            _mm256_store_si256((__m256i *)tmpstore, vtmp2);

            // search the masked primes for non-zero values and copy
            // them to the output array.
            bmask = bmask & 0x88888888;
            bmask2 = bmask2 & 0x88888888;

            while ((id = _tzcnt_u32(bmask)) < 32)
            {
                prime = tmpstore[id / 4];

                // we only use this routine when we need primes for sieving
                // larger ranges (make sure of this), so we can simplify some
                // of the checking...
                if (prime <= ohigh)
                {
                    primes[GLOBAL_OFFSET + pcount++] = prime;
                }

                bmask = _blsr_u32(bmask);
            }

            // one vector store.  (extracting directly from the vector is not profitable...)
            _mm256_store_si256((__m256i *)tmpstore, vtmp4);

            while ((id = _tzcnt_u32(bmask2)) < 32)
            {
                prime = tmpstore[id / 4];

                // we only use this routine when we need primes for sieving
                // larger ranges (make sure of this), so we can simplify some
                // of the checking...
                if (prime <= ohigh)
                {
                    primes[GLOBAL_OFFSET + poffset + pcount2] = prime;
                    pcount2++;
                }

                bmask2 = _blsr_u32(bmask2);
            }
        }

        pcount += pcount2;

    }
    else
    {
        // todo: make a 2 line AVX2 version
        //get 8 bytes from each residue class and pack into a series of 64 bit words.
        //then we can raster across those 64 bits much more efficiently
        memset(cache_word, 0, 64 * sizeof(uint64));
        for (current_line = 0; current_line < nc; current_line++)
        {
            //put 8 bytes from the current line into each of 8 different 64 bit words.
            //shift the byte left according to the current line mod 8 so that
            //each 64 bit word will eventually hold bytes from up to 8 different lines.
            //the bytes from the current line are spaced 8 words apart, so that there is
            //room to store up to 64 lines (capacity enough for the 48 line case mod 210).
            uint32 line_div8 = current_line >> 3;
            uint32 line_mod8 = current_line & 7;
            cache_word[line_div8] |= ((uint64)lines[current_line][byte_offset] << (line_mod8 << 3));
            cache_word[8 + line_div8] |= ((uint64)lines[current_line][byte_offset + 1] << (line_mod8 << 3));
            cache_word[16 + line_div8] |= ((uint64)lines[current_line][byte_offset + 2] << (line_mod8 << 3));
            cache_word[24 + line_div8] |= ((uint64)lines[current_line][byte_offset + 3] << (line_mod8 << 3));
            cache_word[32 + line_div8] |= ((uint64)lines[current_line][byte_offset + 4] << (line_mod8 << 3));
            cache_word[40 + line_div8] |= ((uint64)lines[current_line][byte_offset + 5] << (line_mod8 << 3));
            cache_word[48 + line_div8] |= ((uint64)lines[current_line][byte_offset + 6] << (line_mod8 << 3));
            cache_word[56 + line_div8] |= ((uint64)lines[current_line][byte_offset + 7] << (line_mod8 << 3));
        }

        //for each bit
        for (b = 0; b < 64; b++)
        {
            for (current_line = 0; current_line < nc; current_line++)
            {
                //compute the prime at this location if it is flagged and 
                //within our original boundaries.  
                //if (lines[current_line][i] & nmasks[b])
                //if (cache_word & nmasks64[b])
                //if (cache_word & ((uint64)nmasks[b] << (current_line << 3)))
                //if (cache_word[current_line >> 3] & nmasks64[b][current_line & 7])
                //select the appropriate word according to the bit and line.
                //then 'and' it with the appropriate mask according to the bit and line.
                //all these bit operations are cheaper than continually fetching new
                //bytes from many different lines (cache optimization)
                if (cache_word[((b >> 3) << 3) + (current_line >> 3)] &
                    nmasks64[b & 7][current_line & 7])
                {
                    prime = prodN * ((byte_offset << 3) + b) + rclass[current_line] + lowlimit;

                    if ((prime >= olow) && (prime <= ohigh))
                    {
                        if (NO_STORE)
                            pcount++;
                        else
                            primes[GLOBAL_OFFSET + pcount++] = prime;
                    }
                }
            }
        }

    }
    
#else

	//get 8 bytes from each residue class and pack into a series of 64 bit words.
	//then we can raster across those 64 bits much more efficiently
	memset(cache_word, 0, 64 * sizeof(uint64));
	for (current_line = 0; current_line < nc; current_line++)
	{
		//put 8 bytes from the current line into each of 8 different 64 bit words.
		//shift the byte left according to the current line mod 8 so that
		//each 64 bit word will eventually hold bytes from up to 8 different lines.
		//the bytes from the current line are spaced 8 words apart, so that there is
		//room to store up to 64 lines (capacity enough for the 48 line case mod 210).
		uint32 line_div8 = current_line >> 3;
		uint32 line_mod8 = current_line & 7;
		cache_word[line_div8] |= ((uint64)lines[current_line][byte_offset] << (line_mod8 << 3));
		cache_word[8 + line_div8] |= ((uint64)lines[current_line][byte_offset+1] << (line_mod8 << 3));
		cache_word[16 + line_div8] |= ((uint64)lines[current_line][byte_offset+2] << (line_mod8 << 3));
		cache_word[24 + line_div8] |= ((uint64)lines[current_line][byte_offset+3] << (line_mod8 << 3));
		cache_word[32 + line_div8] |= ((uint64)lines[current_line][byte_offset+4] << (line_mod8 << 3));
		cache_word[40 + line_div8] |= ((uint64)lines[current_line][byte_offset+5] << (line_mod8 << 3));
		cache_word[48 + line_div8] |= ((uint64)lines[current_line][byte_offset+6] << (line_mod8 << 3));
		cache_word[56 + line_div8] |= ((uint64)lines[current_line][byte_offset+7] << (line_mod8 << 3));
	}

	//for each bit
	for (b = 0; b < 64; b++)
	{
		for (current_line = 0; current_line < nc; current_line++)
		{
			//compute the prime at this location if it is flagged and 
			//within our original boundaries.  
			//if (lines[current_line][i] & nmasks[b])
			//if (cache_word & nmasks64[b])
			//if (cache_word & ((uint64)nmasks[b] << (current_line << 3)))
			//if (cache_word[current_line >> 3] & nmasks64[b][current_line & 7])
			//select the appropriate word according to the bit and line.
			//then 'and' it with the appropriate mask according to the bit and line.
			//all these bit operations are cheaper than continually fetching new
			//bytes from many different lines (cache optimization)
			if (cache_word[((b >> 3) << 3) + (current_line >> 3)] & 
				nmasks64[b & 7][current_line & 7])
			{
				prime = prodN * ((byte_offset << 3) + b) + rclass[current_line] + lowlimit;

				if ((prime >= olow) && (prime <= ohigh))
				{
					if (NO_STORE)
						pcount++;
					else
						primes[GLOBAL_OFFSET + pcount++] = prime;
				}
			}
		}
	}

#endif


    // alternate (slightly slower) AVX2 code
#ifdef NOTDEF
    // for each of the 32 bits, compute a prime if the bit is set.
    // we operate on each line's 32-bit block simultaneously and via
    // the gather operation, the output primes will be contiguous.
    for (i = 0; i < 32; i++)
    {
        int id;

        // isolate the current bit and compare if equal to 1.            
        vtmp2 = _mm256_cmpeq_epi32(_mm256_and_si256(vmasks, vbits), _mm256_set1_epi32(1));

        // extract topmost bit of every byte of the comparison mask
        bmask = _mm256_movemask_epi8(vtmp2);

        // next bitmap location
        vbits = _mm256_srli_epi32(vbits, 1);

        // skip if all non-prime
        if (bmask == 0)
            continue;

        // compute the primes.  The only part of the calculation 
        // that involves more than 32 bits is the addition of
        // lowlimit.  This is because when we compute primes the
        // range is not allowed to exceed 32-bits in size.
        vtmp1 = _mm256_set1_epi32(i);
        vtmp1 = _mm256_add_epi32(vtmp1, vbitoffset);
        vtmp1 = _mm256_mullo_epi32(vtmp1, vprodN);
        vtmp1 = _mm256_add_epi32(vtmp1, vlinenums);

        // now we have to split the 8x vector into two 4x vectors
        // and add in the lowlimit.  for now assume lowlimit is zero...
        // should create more than one case to handle this.



        // Use the result to conditionally store the primes we computed.
        vtmp2 = _mm256_and_si256(vtmp2, vtmp1);

        // one vector store.  (extracting directly from the vector is not profitable...)
        _mm256_store_si256((__m256i *)tmpstore, vtmp2);

        // search the masked primes for non-zero values and copy
        // them to the output array.
        bmask = bmask & 0x88888888;

        while ((id = _tzcnt_u32(bmask)) < 32)
        {
            prime = tmpstore[id / 4];

            // we only use this routine when we need primes for sieving
            // larger ranges (make sure of this), so we can simplify some
            // of the checking...
            if (prime <= ohigh)
            {
                primes[GLOBAL_OFFSET + pcount++] = prime;
            }

            bmask = _blsr_u32(bmask);
        }
    }

    // gather the second 32-bits of each of our 8 lines in one vector
    startaddr = (int *)(&lines[0][byte_offset + 4]);
    vbits = _mm256_i32gather_epi32(startaddr, voffsets, 1);
    vbitoffset = _mm256_set1_epi32((byte_offset + 4) << 3);

    // for each of the 32 bits, compute a prime if the bit is set.
    // we operate on each line's 32-bit block simultaneously and via
    // the gather operation, the output primes will be contiguous.
    for (i = 0; i < 32; i++)
    {
        int id;

        // isolate the current bit and compare if equal to 1.            
        vtmp2 = _mm256_cmpeq_epi32(_mm256_and_si256(vmasks, vbits), _mm256_set1_epi32(1));

        // extract topmost bit of every byte of the comparison mask
        bmask = _mm256_movemask_epi8(vtmp2);

        // next bitmap location
        vbits = _mm256_srli_epi32(vbits, 1);

        // skip if all non-prime
        if (bmask == 0)
            continue;

        // compute the primes.  The only part of the calculation 
        // that involves more than 32 bits is the addition of
        // lowlimit.  This is because when we compute primes the
        // range is not allowed to exceed 32-bits in size.
        // prime = prodN * ((byte_offset << 3) + b) + rclass[current_line] + lowlimit;
        vtmp1 = _mm256_set1_epi32(i);
        vtmp1 = _mm256_add_epi32(vtmp1, vbitoffset);
        vtmp1 = _mm256_mullo_epi32(vtmp1, vprodN);
        vtmp1 = _mm256_add_epi32(vtmp1, vlinenums);

        // now we have to split the 8x vector into two 4x vectors
        // and add in the lowlimit.  for now assume lowlimit is zero...
        // should create more than one case to handle this.     


        // Use the result to conditionally store the primes we computed.
        vtmp2 = _mm256_and_si256(vtmp2, vtmp1);

        // one vector store.  (extracting directly from the vector is not profitable...)
        _mm256_store_si256((__m256i *)tmpstore, vtmp2);

        // search the masked primes for non-zero values and copy
        // them to the output array.
        bmask = bmask & 0x88888888;

        while ((id = _tzcnt_u32(bmask)) < 32)
        {
            prime = tmpstore[id / 4];

            // we only use this routine when we need primes for sieving
            // larger ranges (make sure of this), so we can simplify some
            // of the checking...
            if (prime <= ohigh)
            {
                primes[GLOBAL_OFFSET + pcount++] = prime;
            }

            bmask = _blsr_u32(bmask);
        }

    }
#endif


	return pcount;
}

