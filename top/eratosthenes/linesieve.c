/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 7/28/10
----------------------------------------------------------------------*/

#include "soe.h"
#include <immintrin.h>



void sieve_line(thread_soedata_t *thread_data)
{
	//extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;	
	uint32 current_line = thread_data->current_line;
	uint8 *line = thread_data->sdata.lines[current_line];
	
	//stuff for bucket sieving
    uint64 *bptr;
    uint64 **buckets;

	uint32 *nptr;
	uint32 linesize = FLAGSIZE * sdata->blocks, bnum;

	uint8 *flagblock;
	uint64 i,j,k;
	uint32 prime;
	uint32 maxP;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (uint32)(sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	//for the current line, find the offsets of each small/med prime past the low limit
	//and bucket sieve large primes
	get_offsets(thread_data);

	flagblock = line;
	for (i=0;i<sdata->blocks;i++)
	{			
#ifdef USE_AVX2
        __m256i vflagsizem1 = _mm256_set1_epi32(FLAGSIZEm1);
        __m256i v31 = _mm256_set1_epi32(31);
        __m256i vfull = _mm256_set1_epi32(0xffffffff);
        __m256i vone = _mm256_set1_epi32(1);

#ifdef __INTEL_COMPILER
        __declspec(align(64)) uint32 t[8];
        __declspec(align(64)) uint32 t2[8];
#else
        uint32 *t;
        uint32 *t2;
        t = (uint32 *)xmalloc_align(8 * sizeof(uint32));
        t2 = (uint32 *)xmalloc_align(8 * sizeof(uint32));
#endif

        uint32 *flagblock32 = (uint32 *)flagblock;
#endif

		//set all flags for this block, which also puts it into cache for the sieving
		//to follow
		memset(flagblock,255,BLOCKSIZE);			
		
		//smallest primes use special methods
		pre_sieve(ddata, sdata, flagblock);
		
		//one is not a prime
		if (sdata->sieve_range == 0)
		{
			if ((sdata->rclass[current_line] == 1) &&
				(sdata->lowlimit <= 1) && (i == 0))
				flagblock[0] &= 0xfe;
		}
		else
		{
			if ((sdata->rclass[current_line] == 1) &&
				(mpz_cmp_ui(*sdata->offset, 1) <= 0) && (i == 0))
				flagblock[0] &= 0xfe;
		}		

#ifdef USE_AVX2
        // AVX2 presieves more, so we start with the 24th prime (97)
        j = 24;
#else
        j = 10;
#endif

        //unroll the loop: all primes less than this max hit the interval at least 16 times
        maxP = FLAGSIZE >> 4;

		for (;j<ddata->pbounds[i];j++)
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

		
		if (ddata->bucket_depth > 0)
		{
			//finish the primes greater than (flagblocklimit >> 2) that we
			//didn't unroll prior to proceeding with the bucket sieving.
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

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j=0; j < (nptr[i] & (uint32)(~7)); j+=8)
			{				
				//unload 8 hits

                flagblock[(bptr[j + 0] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 0] & FLAGSIZEm1) & 7];
                flagblock[(bptr[j + 1] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 1] & FLAGSIZEm1) & 7];
                flagblock[(bptr[j + 2] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 2] & FLAGSIZEm1) & 7];
                flagblock[(bptr[j + 3] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 3] & FLAGSIZEm1) & 7];
                flagblock[(bptr[j + 4] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 4] & FLAGSIZEm1) & 7];
                flagblock[(bptr[j + 5] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 5] & FLAGSIZEm1) & 7];
                flagblock[(bptr[j + 6] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 6] & FLAGSIZEm1) & 7];
                flagblock[(bptr[j + 7] & FLAGSIZEm1) >> 3] &= masks[(bptr[j + 7] & FLAGSIZEm1) & 7];
				
                bptr[j + 0] += (bptr[j + 0] >> 32);	
                bptr[j + 1] += (bptr[j + 1] >> 32);		
                bptr[j + 2] += (bptr[j + 2] >> 32);	
                bptr[j + 3] += (bptr[j + 3] >> 32);
                bptr[j + 4] += (bptr[j + 4] >> 32);
                bptr[j + 5] += (bptr[j + 5] >> 32);
                bptr[j + 6] += (bptr[j + 6] >> 32);
                bptr[j + 7] += (bptr[j + 7] >> 32);

				//then compute their next hit and update the roots while they are
				//still fresh in the cache
                
                if ((uint32)bptr[j + 0] < linesize)
                {
                    bnum = ((uint32)bptr[j + 0] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 0];
                    nptr[bnum]++;
                }

                if ((uint32)bptr[j + 1] < linesize)
                {
                    bnum = ((uint32)bptr[j + 1] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 1];
                    nptr[bnum]++;
                }

                if ((uint32)bptr[j + 2] < linesize)
                {
                    bnum = ((uint32)bptr[j + 2] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 2];
                    nptr[bnum]++;
                }

                if ((uint32)bptr[j + 3] < linesize)
                {
                    bnum = ((uint32)bptr[j + 3] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 3];
                    nptr[bnum]++;
                }

                if ((uint32)bptr[j + 4] < linesize)
                {
                    bnum = ((uint32)bptr[j + 4] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 4];
                    nptr[bnum]++;
                }

                if ((uint32)bptr[j + 5] < linesize)
                {
                    bnum = ((uint32)bptr[j + 5] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 5];
                    nptr[bnum]++;
                }

                if ((uint32)bptr[j + 6] < linesize)
                {
                    bnum = ((uint32)bptr[j + 6] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 6];
                    nptr[bnum]++;
                }

                if ((uint32)bptr[j + 7] < linesize)
                {
                    bnum = ((uint32)bptr[j + 7] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j + 7];
                    nptr[bnum]++;
                }
                

			}

			//finish up those that didn't fit into a group of 8 hits
			for (;j < nptr[i]; j++)
			{
                flagblock[((uint32)bptr[j] & FLAGSIZEm1) >> 3] &= masks[((uint32)bptr[j] & FLAGSIZEm1) & 7];
                
                bptr[j] += (bptr[j] >> 32);
                if ((uint32)bptr[j] < linesize)
                {
                    bnum = ((uint32)bptr[j] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j];
                    nptr[bnum]++;
                }
                
			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
            {                
				uint32 *large_bptr = ddata->large_sieve_buckets[i];	
				uint32 *large_nptr = ddata->large_bucket_hits;

                for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{				
					//unload 8 hits
#ifdef NOTUSE_AVX2
                    __m256i vbuck, vbuck2;
                    __m256i vt1;
                    __m256i vt2;                
                    __m128i vext;

                    // The AVX2 is almost exactly the same speed as the non-AVX2 code...
                    // the bottleneck is not in computing the indices.
                    // keep it here for future reference: scatter might help eventually.
                    vbuck = _mm256_load_si256((__m256i *)(large_bptr + j));
                    vt1 = _mm256_and_si256(vbuck, vflagsizem1);
                    vt2 = _mm256_and_si256(vt1, v31);
                    vt2 = _mm256_sllv_epi32(vone, vt2);        // bit location
                    vt2 = _mm256_andnot_si256(vt2, vfull);      // sieve &= not(bit location)
                    vt1 = _mm256_srai_epi32(vt1, 5);
                    _mm256_store_si256((__m256i *)t, vt1);
                    //_mm256_store_si256((__m256i *)t2, vt2);
                    vext = _mm256_extracti128_si256(vt2, 0);

                    flagblock32[t[0]] &= _mm_extract_epi32(vext, 0); //t2[0];
                    flagblock32[t[1]] &= _mm_extract_epi32(vext, 1); //t2[1];
                    flagblock32[t[2]] &= _mm_extract_epi32(vext, 2); //t2[2];
                    flagblock32[t[3]] &= _mm_extract_epi32(vext, 3); //t2[3];
                    vext = _mm256_extracti128_si256(vt2, 1);
                    flagblock32[t[4]] &= _mm_extract_epi32(vext, 0); //t2[4];
                    flagblock32[t[5]] &= _mm_extract_epi32(vext, 1); //t2[5];
                    flagblock32[t[6]] &= _mm_extract_epi32(vext, 2); //t2[6];
                    flagblock32[t[7]] &= _mm_extract_epi32(vext, 3); //t2[7];

                    vbuck = _mm256_load_si256((__m256i *)(large_bptr + j + 8));
                    vt1 = _mm256_and_si256(vbuck, vflagsizem1);
                    vt2 = _mm256_and_si256(vt1, v31);
                    vt2 = _mm256_sllv_epi32(vone, vt2);        // bit location
                    vt2 = _mm256_andnot_si256(vt2, vfull);      // sieve &= not(bit location)
                    vt1 = _mm256_srai_epi32(vt1, 5);
                    
                    _mm256_store_si256((__m256i *)t, vt1);
                    //_mm256_store_si256((__m256i *)t2, vt2);
                    vext = _mm256_extracti128_si256(vt2, 0);

                    flagblock32[t[0]] &= _mm_extract_epi32(vext, 0); //t2[0];
                    flagblock32[t[1]] &= _mm_extract_epi32(vext, 1); //t2[1];
                    flagblock32[t[2]] &= _mm_extract_epi32(vext, 2); //t2[2];
                    flagblock32[t[3]] &= _mm_extract_epi32(vext, 3); //t2[3];
                    vext = _mm256_extracti128_si256(vt2, 1);
                    flagblock32[t[4]] &= _mm_extract_epi32(vext, 0); //t2[4];
                    flagblock32[t[5]] &= _mm_extract_epi32(vext, 1); //t2[5];
                    flagblock32[t[6]] &= _mm_extract_epi32(vext, 2); //t2[6];
                    flagblock32[t[7]] &= _mm_extract_epi32(vext, 3); //t2[7];

#else
					flagblock[(large_bptr[j + 0] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 0] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 1] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 1] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 2] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 2] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 3] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 3] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 4] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 4] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 5] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 5] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 6] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 6] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 7] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 7] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 8] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 8] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 9] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 9] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 10] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 10] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 11] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 11] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 12] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 12] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 13] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 13] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 14] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 14] & FLAGSIZEm1) & 7];
					flagblock[(large_bptr[j + 15] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j + 15] & FLAGSIZEm1) & 7];
#endif
				}

				for (;j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & FLAGSIZEm1) >> 3] &= masks[(large_bptr[j] & FLAGSIZEm1) & 7];

			}

		}
		else
		{
			//didn't need to use a bucket sieve
			//finish with primes greater than (flagblocklimit >> 2) that we
			//didn't unroll.
			for (;j<ddata->pbounds[i];j++)
			{
				prime = sdata->sieve_p[j];
				for (k=ddata->offsets[j];k<FLAGSIZE;k+=prime)
					flagblock[k>>3] &= masks[k&7];

				ddata->offsets[j]= (uint32)(k - FLAGSIZE);
			}
		}

		if (i == 0)
		{
			for (j = 0; j < FLAGSIZE; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					break;
				}
			}
		}

		flagblock += BLOCKSIZE;

#ifdef USE_AVX2
#ifndef __INTEL_COMPILER
        align_free(t);
        align_free(t2);
#endif
#endif

	}

	return;
}



