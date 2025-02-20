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

#include "soe.h"
#include "soe_impl.h"
#include "ytools.h"
#if defined(_MSC_VER) && defined(__clang__)
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif
#include <math.h>
#include <string.h>

#define BUCKET_UPDATE2(i,bits) \
    bnum = ((uint32_t)bptr[j + i] >> (bits)); \
    buckets[bnum][nptr[bnum]] = bptr[j + i]; \
    nptr[bnum]++;

// addressing words vs. bytes and methods for setting bits
// have different sweet spots for different compilers...
// these were tested with 10^11 on AVX512 5122 gold cpu
// using the avx2_32k linesieve and avx2 presieving.
#if defined(_INTEL_COMPILER) || defined(__clang__)
#define BITLOGIC32  /* 15.68 */
//#define BITMASKS8  
#elif defined(_MSC_VER)
#define BITMASKS8  /* 17.56 */
#else // gcc, mingw64-gcc
//#define BITMASKS8  /* 17.91 */
#define BITMASKS32 /* 16.69 */
//#define BITLOGIC8  /* 19.02 */
//#define BITLOGIC32   /* 17.20 */
#endif

// sieve all blocks of a line, i.e., a row of the sieve area.
void sieve_line(thread_soedata_t *thread_data)
{
	// extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;	
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = thread_data->sdata.lines[current_line];
    int FLAGSIZE = sdata->FLAGSIZE;
    int FLAGSIZEm1 = sdata->FLAGSIZEm1;
    int FLAGBITS = sdata->FLAGBITS;
    int SOEBLOCKSIZE = sdata->SOEBLOCKSIZE;
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;
	
	// stuff for bucket sieving
    uint64_t *bptr;
    uint64_t **buckets;

	uint32_t *nptr;
	uint32_t linesize = FLAGSIZE * sdata->blocks, bnum;

	uint8_t *flagblock;
	uint64_t i,j,k;
	uint32_t prime;
    int stopid;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	// for the current line, find the offsets of each small/med prime past the low limit
	// and bucket sieve large primes
	get_offsets(thread_data);

    flagblock = line;
	for (i = 0; i < sdata->blocks; i++)
	{
        if (sdata->num_bitmap_primes == 0)
        {
            // set all flags for this block, which also puts it into cache for the sieving
            // to follow.  only do this if we are not bitmap sieving, because
            // that happens before the linesieve and we don't want to overwrite it.
            memset(flagblock, 255, SOEBLOCKSIZE);
        }

        if (sdata->sieve_range == 0)
        {
			// one is not a prime
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
		
		// smallest primes use special methods
		pre_sieve_ptr(ddata, sdata, flagblock);

        // start where presieving left off, which is different for various cpus.
        j = sdata->presieve_max_id;

        // unroll the loop: all primes less than this max hit the interval at least 16 times
        // at first glance it looks like AVX2 operations to compute the indices
        // might be helpful, but since we can't address memory with SIMD registers
        // it actually isn't.  Things might be different with a scatter operation.
		stopid = MIN(ddata->pbounds[i], 1901); //22998); // 1901);
		for (;j<stopid;j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1,p2,p3;

            // we store byte masks to speed things up.  The byte masks
            // are addressed by index % 8.  However, (k+p) % 8 is
            // the same as (k+8p) % 8 so it suffices to compute and
            // store in registers (k+np) % 8 for n = 0:7.  Thanks to 
            // tverniquet for discovering this optimization.
            // The precomputed mask trick works well here because many
            // of these primes actually hit the interval many more
            // than 16 times, thus we get a lot of reuse out of
            // the masks.  In subsequent loops this isn't true.
            uint8_t m0;
            uint8_t m1;
            uint8_t m2;
            uint8_t m3;
            uint8_t m4;
            uint8_t m5;
            uint8_t m6;
            uint8_t m7;

			prime = sdata->sieve_p[j];

			tmpP = prime << 4;
			stop = FLAGSIZE - tmpP + prime;
			k=ddata->offsets[j];            

            p1 = prime;
            p2 = p1 + prime;
            p3 = p2 + prime;

            m0 = masks[k & 7];
            m1 = masks[(k + p1) & 7];
            m2 = masks[(k + p2) & 7];
            m3 = masks[(k + p3) & 7];
            m4 = masks[(k + 4 * prime) & 7];
            m5 = masks[(k + 5 * prime) & 7];
            m6 = masks[(k + 6 * prime) & 7];
            m7 = masks[(k + 7 * prime) & 7];

			while (k < stop)
			{
                flagblock[k >> 3] &= m0;
                flagblock[(k + p1) >> 3] &= m1;
                flagblock[(k + p2) >> 3] &= m2;
                flagblock[(k + p3) >> 3] &= m3;
                k += (prime << 2);
                flagblock[k >> 3] &= m4;
                flagblock[(k + p1) >> 3] &= m5;
                flagblock[(k + p2) >> 3] &= m6;
                flagblock[(k + p3) >> 3] &= m7;
                k += (prime << 2);
                flagblock[k >> 3] &= m0;
                flagblock[(k + p1) >> 3] &= m1;
                flagblock[(k + p2) >> 3] &= m2;
                flagblock[(k + p3) >> 3] &= m3;
                k += (prime << 2);
                flagblock[k >> 3] &= m4;
                flagblock[(k + p1) >> 3] &= m5;
                flagblock[(k + p2) >> 3] &= m6;
                flagblock[(k + p3) >> 3] &= m7;
                k += (prime << 2);
			}

			for (;k<FLAGSIZE;k+=prime)
				flagblock[k>>3] &= masks[k&7];

			ddata->offsets[j]= k - FLAGSIZE;
			
		}

		// unroll the loop: all primes less than this max hit the interval at least 8 times
		stopid = MIN(ddata->pbounds[i], 3513); //43388); // 3513);
        for (; j<stopid; j++)
        {
            uint32_t tmpP;
            uint64_t stop;
            uint64_t p1, p2, p3;

            prime = sdata->sieve_p[j];

            tmpP = prime << 3;
            stop = FLAGSIZE - tmpP + prime;
            k = ddata->offsets[j];
            p1 = prime;
            p2 = p1 + prime;
            p3 = p2 + prime;

            while (k < stop)
            {
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
            }

            for (; k<FLAGSIZE; k += prime)
                flagblock[k >> 3] &= masks[k & 7];

            ddata->offsets[j] = (uint32_t)(k - FLAGSIZE);

        }

        // unroll the loop: all primes less than this max hit the interval at least 4 times
		stopid = MIN(ddata->pbounds[i], 6543); //82023); // 
        for (; j<stopid; j++)
        {
            uint32_t tmpP;
            uint64_t stop;
            uint64_t p1, p2, p3;

            prime = sdata->sieve_p[j];

            tmpP = prime << 2;
            stop = FLAGSIZE - tmpP + prime;
            k = ddata->offsets[j];
            p1 = prime;
            p2 = p1 + prime;
            p3 = p2 + prime;
            while (k < stop)
            {
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
            }

            for (; k<FLAGSIZE; k += prime)
                flagblock[k >> 3] &= masks[k & 7];

            ddata->offsets[j] = (uint32_t)(k - FLAGSIZE);

        }

        // primes are getting fairly big now, unrolling is less useful.
        // keep going up to the large prime bound.
        //stopid = MIN(ddata->pbounds[i], 23002);
		stopid = ddata->pbounds[i];
        for (; j<stopid; j++)
        {
            prime = sdata->sieve_p[j];

			if (prime > FLAGSIZE)
				break;

#ifdef _INTEL_COMPILER
#pragma nounroll
#endif
            for (k = ddata->offsets[j]; k < FLAGSIZE; k += prime)
            {
                flagblock[k >> 3] &= masks[k & 7];
            }

            ddata->offsets[j] = (uint32_t)(k - FLAGSIZE);
        }

        for (; j<ddata->pbounds[i]; j++)
        {
            k = ddata->offsets[j];
            if (ddata->offsets[j] < FLAGSIZE)
            {
                flagblock[k >> 3] &= masks[k & 7];
                k += sdata->sieve_p[j];
            }
            ddata->offsets[j] = (uint32_t)(k - FLAGSIZE);
        }

		// finally, fill any primes in this block's bucket
		if (ddata->bucket_depth > 0)
		{		
			bptr = ddata->sieve_buckets[i];	
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;            

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j=0; j < (nptr[i] & (uint32_t)(~7)); j+=8)
			{		
				// unload 8 hits

                flagblock[(bptr[j + 0] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 0] & 7];
                flagblock[(bptr[j + 1] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 1] & 7];
                flagblock[(bptr[j + 2] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 2] & 7];
                flagblock[(bptr[j + 3] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 3] & 7];
                flagblock[(bptr[j + 4] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 4] & 7];
                flagblock[(bptr[j + 5] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 5] & 7];
                flagblock[(bptr[j + 6] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 6] & 7];
                flagblock[(bptr[j + 7] & FLAGSIZEm1) >> 3] &= masks[bptr[j + 7] & 7];
				
                // then compute their next hit and update the roots while they are
                // still fresh in the cache

#define BUCKET_UPDATE(i) \
    bnum = ((uint32_t)bptr[j + i] >> FLAGBITS); \
    buckets[bnum][nptr[bnum]] = bptr[j + i]; \
    nptr[bnum]++;

                bptr[j + 0] += (bptr[j + 0] >> 32);	
                bptr[j + 1] += (bptr[j + 1] >> 32);		
                bptr[j + 2] += (bptr[j + 2] >> 32);	
                bptr[j + 3] += (bptr[j + 3] >> 32);
                bptr[j + 4] += (bptr[j + 4] >> 32);
                bptr[j + 5] += (bptr[j + 5] >> 32);
                bptr[j + 6] += (bptr[j + 6] >> 32);
                bptr[j + 7] += (bptr[j + 7] >> 32);

                if ((uint32_t)bptr[j + 0] < linesize) { BUCKET_UPDATE(0) }
                if ((uint32_t)bptr[j + 1] < linesize) { BUCKET_UPDATE(1) }
                if ((uint32_t)bptr[j + 2] < linesize) { BUCKET_UPDATE(2) }
                if ((uint32_t)bptr[j + 3] < linesize) { BUCKET_UPDATE(3) }
                if ((uint32_t)bptr[j + 4] < linesize) { BUCKET_UPDATE(4) }
                if ((uint32_t)bptr[j + 5] < linesize) { BUCKET_UPDATE(5) }
                if ((uint32_t)bptr[j + 6] < linesize) { BUCKET_UPDATE(6) }
                if ((uint32_t)bptr[j + 7] < linesize) { BUCKET_UPDATE(7) }
				
			}

			//finish up those that didn't fit into a group of 8 hits
			for (;j < nptr[i]; j++)
			{
                flagblock[(bptr[j] & FLAGSIZEm1) >> 3] &= masks[bptr[j] & 7];
                
                bptr[j] += (bptr[j] >> 32);
                if ((uint32_t)bptr[j] < linesize)
                {
                    bnum = ((uint32_t)bptr[j] >> FLAGBITS);
                    buckets[bnum][nptr[bnum]] = bptr[j];
                    nptr[bnum]++;
                }
                
			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
            {    
				uint32_t *large_bptr = ddata->large_sieve_buckets[i];	
				uint32_t *large_nptr = ddata->large_bucket_hits;

                for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{				
					// unload 16 hits
					flagblock[(large_bptr[j + 0] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 0] & 7];
					flagblock[(large_bptr[j + 1] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 1] & 7];
					flagblock[(large_bptr[j + 2] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 2] & 7];
					flagblock[(large_bptr[j + 3] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 3] & 7];
					flagblock[(large_bptr[j + 4] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 4] & 7];
					flagblock[(large_bptr[j + 5] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 5] & 7];
					flagblock[(large_bptr[j + 6] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 6] & 7];
					flagblock[(large_bptr[j + 7] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 7] & 7];
					flagblock[(large_bptr[j + 8] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 8] & 7];
					flagblock[(large_bptr[j + 9] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 9] & 7];
					flagblock[(large_bptr[j + 10] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 10] & 7];
					flagblock[(large_bptr[j + 11] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 11] & 7];
					flagblock[(large_bptr[j + 12] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 12] & 7];
					flagblock[(large_bptr[j + 13] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 13] & 7];
					flagblock[(large_bptr[j + 14] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 14] & 7];
					flagblock[(large_bptr[j + 15] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j + 15] & 7];
				}

				for (;j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & FLAGSIZEm1) >> 3] &= masks[large_bptr[j] & 7];

			}

		}

		if (i == 0)
		{
			// find the smallest value indicated as prime by the main sieve.
			for (j = 0; j < FLAGSIZE; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					//printf("bit position %d is set in class %d, p = %u\n",
					//	j, sdata->rclass[current_line], ddata->min_sieved_val);
					break;
				}
			}
		}

		flagblock += SOEBLOCKSIZE;
	}

	return;
}

#ifdef USE_AVX512F
void sieve_line_avx512_32k(thread_soedata_t *thread_data)
{
	// extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = thread_data->sdata.lines[current_line];
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

	// stuff for bucket sieving
	uint64_t *bptr;
	uint64_t **buckets;

	uint32_t *nptr;
	uint32_t linesize = 262144 * sdata->blocks, bnum;

	uint8_t *flagblock;
	uint64_t i, j, k;
	uint32_t prime;
	//uint32_t maxP;
	int stopid;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	// for the current line, find the offsets of each small/med prime past the low limit
	// and bucket sieve large primes
	get_offsets(thread_data);

    flagblock = line;
	for (i = 0; i < sdata->blocks; i++)
	{
        uint32_t* flagblock32 = (uint32_t*)flagblock;

		ALIGNED_MEM uint32_t t[16];
		ALIGNED_MEM uint32_t t2[16];

		if (sdata->num_bitmap_primes == 0)
		{
			// set all flags for this block, which also puts it into cache for the sieving
			// to follow.  only do this if we are not bitmap sieving, because
			// that happens before the linesieve and we don't want to overwrite it.
			memset(flagblock, 255, 32768);
		}

        // one is not a prime
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

		// smallest primes use special methods
		pre_sieve_ptr(ddata, sdata, flagblock);

		// start where presieving left off, which is different for various cpus.
		j = sdata->presieve_max_id;

		// unroll the loop: all primes less than this max hit the interval at least 16 times
		//maxP = 262144 >> 4;

		// at first glance it looks like AVX2 operations to compute the indices
		// might be helpful, but since we can't address memory with SIMD registers
		// it actually isn't.  Things might be different with a scatter operation.
		stopid = MIN(ddata->pbounds[i], 1901); //22998); // 1901);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			// we store byte masks to speed things up.  The byte masks
			// are addressed by index % 8.  However, (k+p) % 8 is
			// the same as (k+8p) % 8 so it suffices to compute and
			// store in registers (k+np) % 8 for n = 0:7.  Thanks to 
			// tverniquet for discovering this optimization.
			// The precomputed mask trick works well here because many
			// of these primes actually hit the interval many more
			// than 16 times, thus we get a lot of reuse out of
			// the masks.  In subsequent loops this isn't true.
			uint8_t m0;
			uint8_t m1;
			uint8_t m2;
			uint8_t m3;
			uint8_t m4;
			uint8_t m5;
			uint8_t m6;
			uint8_t m7;

			prime = sdata->sieve_p[j];

			tmpP = prime << 4;
			stop = 262144 - tmpP + prime;
			k = ddata->offsets[j];

			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			m0 = masks[k & 7];
			m1 = masks[(k + p1) & 7];
			m2 = masks[(k + p2) & 7];
			m3 = masks[(k + p3) & 7];
			m4 = masks[(k + 4 * prime) & 7];
			m5 = masks[(k + 5 * prime) & 7];
			m6 = masks[(k + 6 * prime) & 7];
			m7 = masks[(k + 7 * prime) & 7];

			while (k < stop)
			{
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
			}

			for (; k < 262144; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = k - 262144;

		}

		// unroll the loop: all primes less than this max hit the interval at least 8 times
		//maxP = 262144 >> 3;

		stopid = MIN(ddata->pbounds[i], 3513); //43388); // 3513);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 3;
			stop = 262144 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 262144; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 262144);

		}

		// unroll the loop: all primes less than this max hit the interval at least 4 times
		//maxP = 262144 >> 2;

		stopid = MIN(ddata->pbounds[i], 6543); //82023); // 
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 2;
			stop = 262144 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 262144; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 262144);

		}

		// primes are getting fairly big now, unrolling is less useful.
		// keep going up to the large prime bound.
		//stopid = MIN(ddata->pbounds[i], 23002);
		stopid = ddata->pbounds[i];
		for (; j < stopid; j++)
		{
			prime = sdata->sieve_p[j];

			if (prime > 262144)
				break;

#ifdef _INTEL_COMPILER
#pragma nounroll
#endif
			for (k = ddata->offsets[j]; k < 262144; k += prime)
			{
				flagblock[k >> 3] &= masks[k & 7];
			}

			ddata->offsets[j] = (uint32_t)(k - 262144);
		}

		for (; j < ddata->pbounds[i]; j++)
		{
			k = ddata->offsets[j];
			if (ddata->offsets[j] < 262144)
			{
				flagblock[k >> 3] &= masks[k & 7];
				k += sdata->sieve_p[j];
			}
			ddata->offsets[j] = (uint32_t)(k - 262144);
		}

		if (ddata->bucket_depth > 0)
		{
			__m512i vt1, vt3;      // temp vectors
			__m512i vlinesize = _mm512_set1_epi32(linesize);
			__mmask16 cmp1;

			// finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j = 0; j < (nptr[i] & (uint32_t)(~7)); j += 8)
			{
				// unload 8 hits

				flagblock[(bptr[j + 0] & 262143) >> 3] &= masks[bptr[j + 0] & 7];
				flagblock[(bptr[j + 1] & 262143) >> 3] &= masks[bptr[j + 1] & 7];
				flagblock[(bptr[j + 2] & 262143) >> 3] &= masks[bptr[j + 2] & 7];
				flagblock[(bptr[j + 3] & 262143) >> 3] &= masks[bptr[j + 3] & 7];
				flagblock[(bptr[j + 4] & 262143) >> 3] &= masks[bptr[j + 4] & 7];
				flagblock[(bptr[j + 5] & 262143) >> 3] &= masks[bptr[j + 5] & 7];
				flagblock[(bptr[j + 6] & 262143) >> 3] &= masks[bptr[j + 6] & 7];
				flagblock[(bptr[j + 7] & 262143) >> 3] &= masks[bptr[j + 7] & 7];

				// then compute their next hit and update the roots while they are
				// still fresh in the cache
				vt1 = _mm512_loadu_si512((__m512i *)(&bptr[j]));
				vt3 = _mm512_srli_epi64(vt1, 32);
				vt1 = _mm512_add_epi64(vt1, vt3);
				_mm512_storeu_si512((__m512i *)(&bptr[j]), vt1);
				cmp1 = _mm512_cmpgt_epu32_mask(vlinesize, vt1);

				if (cmp1 & 0x1) { BUCKET_UPDATE2(0, 18) }
				if (cmp1 & 0x4) { BUCKET_UPDATE2(1, 18) }
				if (cmp1 & 0x10) { BUCKET_UPDATE2(2, 18) }
				if (cmp1 & 0x40) { BUCKET_UPDATE2(3, 18) }
				if (cmp1 & 0x100) { BUCKET_UPDATE2(4, 18) }
				if (cmp1 & 0x400) { BUCKET_UPDATE2(5, 18) }
				if (cmp1 & 0x1000) { BUCKET_UPDATE2(6, 18) }
				if (cmp1 & 0x4000) { BUCKET_UPDATE2(7, 18) }

			}

			//finish up those that didn't fit into a group of 8 hits
			for (; j < nptr[i]; j++)
			{
				flagblock[(bptr[j] & 262143) >> 3] &= masks[bptr[j] & 7];

				bptr[j] += (bptr[j] >> 32);
				if ((uint32_t)bptr[j] < linesize)
				{
					bnum = ((uint32_t)bptr[j] >> 18);
					buckets[bnum][nptr[bnum]] = bptr[j];
					nptr[bnum]++;
				}

			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
			{
				__m512i vbuck;
				__m512i vt1;
				__m512i vt2;
				__m512i v262143 = _mm512_set1_epi32(262143);
				__m512i v31 = _mm512_set1_epi32(31);
				__m512i vfull = _mm512_set1_epi32(0xffffffff);
				__m512i vone = _mm512_set1_epi32(1);

				uint32_t *large_bptr = ddata->large_sieve_buckets[i];
				uint32_t *large_nptr = ddata->large_bucket_hits;

				for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{
					// unload 16 hits
					// The AVX2 is almost exactly the same speed as the non-AVX2 code...
					// the bottleneck is not in computing the indices.
					// keep it here for future reference: scatter might help eventually.
					vbuck = _mm512_loadu_si512((__m256i *)(large_bptr + j));
					vt1 = _mm512_and_epi32(vbuck, v262143);
					vt2 = _mm512_and_epi32(vt1, v31);
					vt2 = _mm512_sllv_epi32(vone, vt2);        // bit location
					vt2 = _mm512_andnot_epi32(vt2, vfull);      // sieve &= not(bit location)
					vt1 = _mm512_srli_epi32(vt1, 5);

					// this is not collision free... and it's not
					// much faster, if at all, from the other way.
					//vbuck = _mm512_i32gather_epi32(vt1, flagblock32, _MM_SCALE_4);
					//vbuck = _mm512_and_epi32(vext, vt2);
					//_mm512_i32scatter_epi32(flagblock32, vt1, vbuck, _MM_SCALE_4);

					_mm512_storeu_si512((__m256i *)t, vt1);
					_mm512_storeu_si512((__m256i *)t2, vt2);

					flagblock32[t[0]] &= t2[0];
					flagblock32[t[1]] &= t2[1];
					flagblock32[t[2]] &= t2[2];
					flagblock32[t[3]] &= t2[3];
					flagblock32[t[4]] &= t2[4];
					flagblock32[t[5]] &= t2[5];
					flagblock32[t[6]] &= t2[6];
					flagblock32[t[7]] &= t2[7];
					flagblock32[t[8]] &= t2[8];
					flagblock32[t[9]] &= t2[9];
					flagblock32[t[10]] &= t2[10];
					flagblock32[t[11]] &= t2[11];
					flagblock32[t[12]] &= t2[12];
					flagblock32[t[13]] &= t2[13];
					flagblock32[t[14]] &= t2[14];
					flagblock32[t[15]] &= t2[15];

				}

				for (; j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & 262143) >> 3] &= masks[large_bptr[j] & 7];

			}

		}

		if (i == 0)
		{
			for (j = 0; j < 262144; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					break;
				}
			}
		}

		flagblock += 32768;
	}

	return;
}

void sieve_line_avx512_128k(thread_soedata_t *thread_data)
{
	// extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = thread_data->sdata.lines[current_line];
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

	// stuff for bucket sieving
	uint64_t *bptr;
	uint64_t **buckets;

	uint32_t *nptr;
	uint32_t linesize = 1048576 * sdata->blocks, bnum;

	uint8_t *flagblock;
	uint64_t i, j, k;
	uint32_t prime;
	//uint32_t maxP;
	int stopid;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	// for the current line, find the offsets of each small/med prime past the low limit
	// and bucket sieve large primes
	get_offsets(thread_data);

	flagblock = line;
	for (i = 0; i < sdata->blocks; i++)
	{
        uint32_t* flagblock32 = (uint32_t*)flagblock;

		ALIGNED_MEM uint32_t t[16];
		ALIGNED_MEM uint32_t t2[16];

		if (sdata->num_bitmap_primes == 0)
		{
			// set all flags for this block, which also puts it into cache for the sieving
			// to follow.  only do this if we are not bitmap sieving, because
			// that happens before the linesieve and we don't want to overwrite it.
			memset(flagblock, 255, 131072);
		}

		// smallest primes use special methods
		pre_sieve_ptr(ddata, sdata, flagblock);

		// one is not a prime
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

		// start where presieving left off, which is different for various cpus.
		j = sdata->presieve_max_id;

		// unroll the loop: all primes less than this max hit the interval at least 16 times
		//maxP = 1048576 >> 4;

		// at first glance it looks like AVX2 operations to compute the indices
		// might be helpful, but since we can't address memory with SIMD registers
		// it actually isn't.  Things might be different with a scatter operation.
		stopid = MIN(ddata->pbounds[i], 6540); // 1901);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			// we store byte masks to speed things up.  The byte masks
			// are addressed by index % 8.  However, (k+p) % 8 is
			// the same as (k+8p) % 8 so it suffices to compute and
			// store in registers (k+np) % 8 for n = 0:7.  Thanks to 
			// tverniquet for discovering this optimization.
			// The precomputed mask trick works well here because many
			// of these primes actually hit the interval many more
			// than 16 times, thus we get a lot of reuse out of
			// the masks.  In subsequent loops this isn't true.
			uint8_t m0;
			uint8_t m1;
			uint8_t m2;
			uint8_t m3;
			uint8_t m4;
			uint8_t m5;
			uint8_t m6;
			uint8_t m7;

			prime = sdata->sieve_p[j];

			tmpP = prime << 4;
			stop = 1048576 - tmpP + prime;
			k = ddata->offsets[j];

			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			m0 = masks[k & 7];
			m1 = masks[(k + p1) & 7];
			m2 = masks[(k + p2) & 7];
			m3 = masks[(k + p3) & 7];
			m4 = masks[(k + 4 * prime) & 7];
			m5 = masks[(k + 5 * prime) & 7];
			m6 = masks[(k + 6 * prime) & 7];
			m7 = masks[(k + 7 * prime) & 7];

			while (k < stop)
			{
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
			}

			for (; k < 1048576; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = k - 1048576;

		}

		// unroll the loop: all primes less than this max hit the interval at least 8 times
		//maxP = 1048576 >> 3;

		stopid = MIN(ddata->pbounds[i], 12249); // 3513);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 3;
			stop = 1048576 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 1048576; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 1048576);

		}

		// unroll the loop: all primes less than this max hit the interval at least 4 times
		//maxP = 1048576 >> 2;

		stopid = MIN(ddata->pbounds[i], 22998); // 
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 2;
			stop = 1048576 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 1048576; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 1048576);

		}

		// primes are getting fairly big now, unrolling is less useful.
		// keep going up to the large prime bound.
		//stopid = MIN(ddata->pbounds[i], 23002);
		stopid = ddata->pbounds[i];
		for (; j < stopid; j++)
		{
			prime = sdata->sieve_p[j];

			if (prime > 1048576)
				break;

#ifdef _INTEL_COMPILER
#pragma nounroll
#endif
			for (k = ddata->offsets[j]; k < 1048576; k += prime)
			{
				flagblock[k >> 3] &= masks[k & 7];
			}

			ddata->offsets[j] = (uint32_t)(k - 1048576);
		}

		for (; j < ddata->pbounds[i]; j++)
		{
			k = ddata->offsets[j];
			if (ddata->offsets[j] < 1048576)
			{
				flagblock[k >> 3] &= masks[k & 7];
				k += sdata->sieve_p[j];
			}
			ddata->offsets[j] = (uint32_t)(k - 1048576);
		}

		if (ddata->bucket_depth > 0)
		{
			__m512i vt1, vt3;      // temp vectors
			__m512i vlinesize = _mm512_set1_epi32(linesize);
			__mmask16 cmp1;

			// finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j = 0; j < (nptr[i] & (uint32_t)(~7)); j += 8)
			{
				// unload 8 hits

				flagblock[(bptr[j + 0] & 1048575) >> 3] &= masks[bptr[j + 0] & 7];
				flagblock[(bptr[j + 1] & 1048575) >> 3] &= masks[bptr[j + 1] & 7];
				flagblock[(bptr[j + 2] & 1048575) >> 3] &= masks[bptr[j + 2] & 7];
				flagblock[(bptr[j + 3] & 1048575) >> 3] &= masks[bptr[j + 3] & 7];
				flagblock[(bptr[j + 4] & 1048575) >> 3] &= masks[bptr[j + 4] & 7];
				flagblock[(bptr[j + 5] & 1048575) >> 3] &= masks[bptr[j + 5] & 7];
				flagblock[(bptr[j + 6] & 1048575) >> 3] &= masks[bptr[j + 6] & 7];
				flagblock[(bptr[j + 7] & 1048575) >> 3] &= masks[bptr[j + 7] & 7];

				// then compute their next hit and update the roots while they are
				// still fresh in the cache
				vt1 = _mm512_loadu_si512((__m512i *)(&bptr[j]));
				vt3 = _mm512_srli_epi64(vt1, 32);
				vt1 = _mm512_add_epi64(vt1, vt3);
				_mm512_storeu_si512((__m512i *)(&bptr[j]), vt1);
				cmp1 = _mm512_cmpgt_epu32_mask(vlinesize, vt1);

				if (cmp1 & 0x1) { BUCKET_UPDATE2(0, 20) }
				if (cmp1 & 0x4) { BUCKET_UPDATE2(1, 20) }
				if (cmp1 & 0x10) { BUCKET_UPDATE2(2, 20) }
				if (cmp1 & 0x40) { BUCKET_UPDATE2(3, 20) }
				if (cmp1 & 0x100) { BUCKET_UPDATE2(4, 20) }
				if (cmp1 & 0x400) { BUCKET_UPDATE2(5, 20) }
				if (cmp1 & 0x1000) { BUCKET_UPDATE2(6, 20) }
				if (cmp1 & 0x4000) { BUCKET_UPDATE2(7, 20) }

			}

			//finish up those that didn't fit into a group of 8 hits
			for (; j < nptr[i]; j++)
			{
				flagblock[(bptr[j] & 1048575) >> 3] &= masks[bptr[j] & 7];

				bptr[j] += (bptr[j] >> 32);
				if ((uint32_t)bptr[j] < linesize)
				{
					bnum = ((uint32_t)bptr[j] >> 20);
					buckets[bnum][nptr[bnum]] = bptr[j];
					nptr[bnum]++;
				}

			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
			{
				__m512i vbuck;
				__m512i vt1;
				__m512i vt2;
				__m512i v1048575 = _mm512_set1_epi32(1048575);
				__m512i v31 = _mm512_set1_epi32(31);
				__m512i vfull = _mm512_set1_epi32(0xffffffff);
				__m512i vone = _mm512_set1_epi32(1);

				uint32_t *large_bptr = ddata->large_sieve_buckets[i];
				uint32_t *large_nptr = ddata->large_bucket_hits;

				for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{
					// unload 16 hits
					// The AVX2 is almost exactly the same speed as the non-AVX2 code...
					// the bottleneck is not in computing the indices.
					// keep it here for future reference: scatter might help eventually.
					vbuck = _mm512_loadu_si512((__m256i *)(large_bptr + j));
					vt1 = _mm512_and_epi32(vbuck, v1048575);
					vt2 = _mm512_and_epi32(vt1, v31);
					vt2 = _mm512_sllv_epi32(vone, vt2);        // bit location
					vt2 = _mm512_andnot_epi32(vt2, vfull);      // sieve &= not(bit location)
					vt1 = _mm512_srli_epi32(vt1, 5);

					// this is not collision free... and it's not
					// much faster, if at all, from the other way.
					//vbuck = _mm512_i32gather_epi32(vt1, flagblock32, _MM_SCALE_4);
					//vbuck = _mm512_and_epi32(vext, vt2);
					//_mm512_i32scatter_epi32(flagblock32, vt1, vbuck, _MM_SCALE_4);

					_mm512_storeu_si512((__m256i *)t, vt1);
					_mm512_storeu_si512((__m256i *)t2, vt2);

					flagblock32[t[0]] &= t2[0];
					flagblock32[t[1]] &= t2[1];
					flagblock32[t[2]] &= t2[2];
					flagblock32[t[3]] &= t2[3];
					flagblock32[t[4]] &= t2[4];
					flagblock32[t[5]] &= t2[5];
					flagblock32[t[6]] &= t2[6];
					flagblock32[t[7]] &= t2[7];
					flagblock32[t[8]] &= t2[8];
					flagblock32[t[9]] &= t2[9];
					flagblock32[t[10]] &= t2[10];
					flagblock32[t[11]] &= t2[11];
					flagblock32[t[12]] &= t2[12];
					flagblock32[t[13]] &= t2[13];
					flagblock32[t[14]] &= t2[14];
					flagblock32[t[15]] &= t2[15];

				}

				for (; j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & 1048575) >> 3] &= masks[large_bptr[j] & 7];

			}

		}

		if (i == 0)
		{
			for (j = 0; j < 1048576; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					break;
				}
			}
		}

		flagblock += 131072;
	}

	return;
}

void sieve_line_avx512_256k(thread_soedata_t *thread_data)
{
	// extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = thread_data->sdata.lines[current_line];
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

	// stuff for bucket sieving
	uint64_t *bptr;
	uint64_t **buckets;

	uint32_t *nptr;
	uint32_t linesize = 2097152 * sdata->blocks, bnum;

	uint8_t *flagblock;
	uint64_t i, j, k;
	uint32_t prime;
	//uint32_t maxP;
	int stopid;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	// for the current line, find the offsets of each small/med prime past the low limit
	// and bucket sieve large primes
	get_offsets(thread_data);

	flagblock = line;
	for (i = 0; i < sdata->blocks; i++)
	{
        uint32_t* flagblock32 = (uint32_t*)flagblock;

		ALIGNED_MEM uint32_t t[16];
		ALIGNED_MEM uint32_t t2[16];

		if (sdata->num_bitmap_primes == 0)
		{
			// set all flags for this block, which also puts it into cache for the sieving
			// to follow.  only do this if we are not bitmap sieving, because
			// that happens before the linesieve and we don't want to overwrite it.
			memset(flagblock, 255, 262144);
		}

		// smallest primes use special methods
		pre_sieve_ptr(ddata, sdata, flagblock);

		// one is not a prime
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

		// start where presieving left off, which is different for various cpus.
		j = sdata->presieve_max_id;

		// unroll the loop: all primes less than this max hit the interval at least 16 times
		//maxP = 2097152 >> 4;

		// at first glance it looks like AVX2 operations to compute the indices
		// might be helpful, but since we can't address memory with SIMD registers
		// it actually isn't.  Things might be different with a scatter operation.
		stopid = MIN(ddata->pbounds[i], 12249); // 1901);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			// we store byte masks to speed things up.  The byte masks
			// are addressed by index % 8.  However, (k+p) % 8 is
			// the same as (k+8p) % 8 so it suffices to compute and
			// store in registers (k+np) % 8 for n = 0:7.  Thanks to 
			// tverniquet for discovering this optimization.
			// The precomputed mask trick works well here because many
			// of these primes actually hit the interval many more
			// than 16 times, thus we get a lot of reuse out of
			// the masks.  In subsequent loops this isn't true.
			uint8_t m0;
			uint8_t m1;
			uint8_t m2;
			uint8_t m3;
			uint8_t m4;
			uint8_t m5;
			uint8_t m6;
			uint8_t m7;

			prime = sdata->sieve_p[j];

			tmpP = prime << 4;
			stop = 2097152 - tmpP + prime;
			k = ddata->offsets[j];

			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			m0 = masks[k & 7];
			m1 = masks[(k + p1) & 7];
			m2 = masks[(k + p2) & 7];
			m3 = masks[(k + p3) & 7];
			m4 = masks[(k + 4 * prime) & 7];
			m5 = masks[(k + 5 * prime) & 7];
			m6 = masks[(k + 6 * prime) & 7];
			m7 = masks[(k + 7 * prime) & 7];

			while (k < stop)
			{
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
			}

			for (; k < 2097152; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = k - 2097152;

		}

		// unroll the loop: all primes less than this max hit the interval at least 8 times
		//maxP = 2097152 >> 3;

		stopid = MIN(ddata->pbounds[i], 22998); // 3513);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 3;
			stop = 2097152 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 2097152; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 2097152);

		}

		// unroll the loop: all primes less than this max hit the interval at least 4 times
		//maxP = 2097152 >> 2;

		stopid = MIN(ddata->pbounds[i], 43388); // 
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 2;
			stop = 2097152 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 2097152; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 2097152);

		}

		// primes are getting fairly big now, unrolling is less useful.
		// keep going up to the large prime bound.
		//stopid = MIN(ddata->pbounds[i], 23002);
		stopid = ddata->pbounds[i];
		for (; j < stopid; j++)
		{
			prime = sdata->sieve_p[j];

			if (prime > 2097152)
				break;

#ifdef _INTEL_COMPILER
#pragma nounroll
#endif
			for (k = ddata->offsets[j]; k < 2097152; k += prime)
			{
				flagblock[k >> 3] &= masks[k & 7];
			}

			ddata->offsets[j] = (uint32_t)(k - 2097152);
		}

		for (; j < ddata->pbounds[i]; j++)
		{
			k = ddata->offsets[j];
			if (ddata->offsets[j] < 2097152)
			{
				flagblock[k >> 3] &= masks[k & 7];
				k += sdata->sieve_p[j];
			}
			ddata->offsets[j] = (uint32_t)(k - 2097152);
		}

		if (ddata->bucket_depth > 0)
		{
			__m512i vt1, vt3;      // temp vectors
			__m512i vlinesize = _mm512_set1_epi32(linesize);
			__mmask16 cmp1;

			// finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j = 0; j < (nptr[i] & (uint32_t)(~7)); j += 8)
			{
				// unload 8 hits

				flagblock[(bptr[j + 0] & 2097151) >> 3] &= masks[bptr[j + 0] & 7];
				flagblock[(bptr[j + 1] & 2097151) >> 3] &= masks[bptr[j + 1] & 7];
				flagblock[(bptr[j + 2] & 2097151) >> 3] &= masks[bptr[j + 2] & 7];
				flagblock[(bptr[j + 3] & 2097151) >> 3] &= masks[bptr[j + 3] & 7];
				flagblock[(bptr[j + 4] & 2097151) >> 3] &= masks[bptr[j + 4] & 7];
				flagblock[(bptr[j + 5] & 2097151) >> 3] &= masks[bptr[j + 5] & 7];
				flagblock[(bptr[j + 6] & 2097151) >> 3] &= masks[bptr[j + 6] & 7];
				flagblock[(bptr[j + 7] & 2097151) >> 3] &= masks[bptr[j + 7] & 7];

				// then compute their next hit and update the roots while they are
				// still fresh in the cache
				vt1 = _mm512_loadu_si512((__m512i *)(&bptr[j]));
				vt3 = _mm512_srli_epi64(vt1, 32);
				vt1 = _mm512_add_epi64(vt1, vt3);
				_mm512_storeu_si512((__m512i *)(&bptr[j]), vt1);
				cmp1 = _mm512_cmpgt_epu32_mask(vlinesize, vt1);

				if (cmp1 & 0x1) { BUCKET_UPDATE2(0, 21) }
				if (cmp1 & 0x4) { BUCKET_UPDATE2(1, 21) }
				if (cmp1 & 0x10) { BUCKET_UPDATE2(2, 21) }
				if (cmp1 & 0x40) { BUCKET_UPDATE2(3, 21) }
				if (cmp1 & 0x100) { BUCKET_UPDATE2(4, 21) }
				if (cmp1 & 0x400) { BUCKET_UPDATE2(5, 21) }
				if (cmp1 & 0x1000) { BUCKET_UPDATE2(6, 21) }
				if (cmp1 & 0x4000) { BUCKET_UPDATE2(7, 21) }

			}

			//finish up those that didn't fit into a group of 8 hits
			for (; j < nptr[i]; j++)
			{
				flagblock[(bptr[j] & 2097151) >> 3] &= masks[bptr[j] & 7];

				bptr[j] += (bptr[j] >> 32);
				if ((uint32_t)bptr[j] < linesize)
				{
					bnum = ((uint32_t)bptr[j] >> 21);
					buckets[bnum][nptr[bnum]] = bptr[j];
					nptr[bnum]++;
				}

			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
			{
				__m512i vbuck;
				__m512i vt1;
				__m512i vt2;
				__m512i v2097151 = _mm512_set1_epi32(2097151);
				__m512i v31 = _mm512_set1_epi32(31);
				__m512i vfull = _mm512_set1_epi32(0xffffffff);
				__m512i vone = _mm512_set1_epi32(1);

				uint32_t *large_bptr = ddata->large_sieve_buckets[i];
				uint32_t *large_nptr = ddata->large_bucket_hits;

				for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{
					// unload 16 hits
					// The AVX2 is almost exactly the same speed as the non-AVX2 code...
					// the bottleneck is not in computing the indices.
					// keep it here for future reference: scatter might help eventually.
					vbuck = _mm512_loadu_si512((__m256i *)(large_bptr + j));
					vt1 = _mm512_and_epi32(vbuck, v2097151);
					vt2 = _mm512_and_epi32(vt1, v31);
					vt2 = _mm512_sllv_epi32(vone, vt2);        // bit location
					vt2 = _mm512_andnot_epi32(vt2, vfull);      // sieve &= not(bit location)
					vt1 = _mm512_srli_epi32(vt1, 5);

					// this is not collision free... and it's not
					// much faster, if at all, from the other way.
					//vbuck = _mm512_i32gather_epi32(vt1, flagblock32, _MM_SCALE_4);
					//vbuck = _mm512_and_epi32(vext, vt2);
					//_mm512_i32scatter_epi32(flagblock32, vt1, vbuck, _MM_SCALE_4);

					_mm512_storeu_si512((__m256i *)t, vt1);
					_mm512_storeu_si512((__m256i *)t2, vt2);

					flagblock32[t[0]] &= t2[0];
					flagblock32[t[1]] &= t2[1];
					flagblock32[t[2]] &= t2[2];
					flagblock32[t[3]] &= t2[3];
					flagblock32[t[4]] &= t2[4];
					flagblock32[t[5]] &= t2[5];
					flagblock32[t[6]] &= t2[6];
					flagblock32[t[7]] &= t2[7];
					flagblock32[t[8]] &= t2[8];
					flagblock32[t[9]] &= t2[9];
					flagblock32[t[10]] &= t2[10];
					flagblock32[t[11]] &= t2[11];
					flagblock32[t[12]] &= t2[12];
					flagblock32[t[13]] &= t2[13];
					flagblock32[t[14]] &= t2[14];
					flagblock32[t[15]] &= t2[15];

				}

				for (; j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & 2097151) >> 3] &= masks[large_bptr[j] & 7];

			}

		}

		if (i == 0)
		{
			for (j = 0; j < 2097152; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					break;
				}
			}
		}

		flagblock += 262144;
	}

	return;
}

void sieve_line_avx512_512k(thread_soedata_t *thread_data)
{
	// extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = thread_data->sdata.lines[current_line];
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

	// stuff for bucket sieving
	uint64_t *bptr;
	uint64_t **buckets;

	uint32_t *nptr;
	uint32_t linesize = 4194304 * sdata->blocks, bnum;

	uint8_t *flagblock;
	uint64_t i, j, k;
	uint32_t prime;
	//uint32_t maxP;
	int stopid;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	// for the current line, find the offsets of each small/med prime past the low limit
	// and bucket sieve large primes
	get_offsets(thread_data);

	flagblock = line;
	for (i = 0; i < sdata->blocks; i++)
	{
        uint32_t* flagblock32 = (uint32_t*)flagblock;

		ALIGNED_MEM uint32_t t[16];
		ALIGNED_MEM uint32_t t2[16];

		if (sdata->num_bitmap_primes == 0)
		{
			// set all flags for this block, which also puts it into cache for the sieving
			// to follow.  only do this if we are not bitmap sieving, because
			// that happens before the linesieve and we don't want to overwrite it.
			memset(flagblock, 255, 524288);
		}

		// smallest primes use special methods
		pre_sieve_ptr(ddata, sdata, flagblock);

		// one is not a prime
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

		// start where presieving left off, which is different for various cpus.
		j = sdata->presieve_max_id;

		// unroll the loop: all primes less than this max hit the interval at least 16 times
		//maxP = 4194304 >> 4;

		// at first glance it looks like AVX2 operations to compute the indices
		// might be helpful, but since we can't address memory with SIMD registers
		// it actually isn't.  Things might be different with a scatter operation.
		stopid = MIN(ddata->pbounds[i], 22998); // 1901);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			// we store byte masks to speed things up.  The byte masks
			// are addressed by index % 8.  However, (k+p) % 8 is
			// the same as (k+8p) % 8 so it suffices to compute and
			// store in registers (k+np) % 8 for n = 0:7.  Thanks to 
			// tverniquet for discovering this optimization.
			// The precomputed mask trick works well here because many
			// of these primes actually hit the interval many more
			// than 16 times, thus we get a lot of reuse out of
			// the masks.  In subsequent loops this isn't true.
			uint8_t m0;
			uint8_t m1;
			uint8_t m2;
			uint8_t m3;
			uint8_t m4;
			uint8_t m5;
			uint8_t m6;
			uint8_t m7;

			prime = sdata->sieve_p[j];

			tmpP = prime << 4;
			stop = 4194304 - tmpP + prime;
			k = ddata->offsets[j];

			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			m0 = masks[k & 7];
			m1 = masks[(k + p1) & 7];
			m2 = masks[(k + p2) & 7];
			m3 = masks[(k + p3) & 7];
			m4 = masks[(k + 4 * prime) & 7];
			m5 = masks[(k + 5 * prime) & 7];
			m6 = masks[(k + 6 * prime) & 7];
			m7 = masks[(k + 7 * prime) & 7];

			while (k < stop)
			{
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
			}

			for (; k < 4194304; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = k - 4194304;

		}

		// unroll the loop: all primes less than this max hit the interval at least 8 times
		//maxP = 4194304 >> 3;

		stopid = MIN(ddata->pbounds[i], 43388); // 3513);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 3;
			stop = 4194304 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 4194304; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 4194304);

		}

		// unroll the loop: all primes less than this max hit the interval at least 4 times
		//maxP = 4194304 >> 2;

		stopid = MIN(ddata->pbounds[i], 82023); // 
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 2;
			stop = 4194304 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 4194304; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 4194304);

		}

		// primes are getting fairly big now, unrolling is less useful.
		// keep going up to the large prime bound.
		//stopid = MIN(ddata->pbounds[i], 23002);
		stopid = ddata->pbounds[i];
		for (; j < stopid; j++)
		{
			prime = sdata->sieve_p[j];

			if (prime > 4194304)
				break;

#ifdef _INTEL_COMPILER
#pragma nounroll
#endif
			for (k = ddata->offsets[j]; k < 4194304; k += prime)
			{
				flagblock[k >> 3] &= masks[k & 7];
			}

			ddata->offsets[j] = (uint32_t)(k - 4194304);
		}

		for (; j < ddata->pbounds[i]; j++)
		{
			k = ddata->offsets[j];
			if (ddata->offsets[j] < 4194304)
			{
				flagblock[k >> 3] &= masks[k & 7];
				k += sdata->sieve_p[j];
			}
			ddata->offsets[j] = (uint32_t)(k - 4194304);
		}

		if (ddata->bucket_depth > 0)
		{
			__m512i vt1, vt3;      // temp vectors
			__m512i vlinesize = _mm512_set1_epi32(linesize);
			__mmask16 cmp1;

			// finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j = 0; j < (nptr[i] & (uint32_t)(~7)); j += 8)
			{
				// unload 8 hits

				flagblock[(bptr[j + 0] & 4194303) >> 3] &= masks[bptr[j + 0] & 7];
				flagblock[(bptr[j + 1] & 4194303) >> 3] &= masks[bptr[j + 1] & 7];
				flagblock[(bptr[j + 2] & 4194303) >> 3] &= masks[bptr[j + 2] & 7];
				flagblock[(bptr[j + 3] & 4194303) >> 3] &= masks[bptr[j + 3] & 7];
				flagblock[(bptr[j + 4] & 4194303) >> 3] &= masks[bptr[j + 4] & 7];
				flagblock[(bptr[j + 5] & 4194303) >> 3] &= masks[bptr[j + 5] & 7];
				flagblock[(bptr[j + 6] & 4194303) >> 3] &= masks[bptr[j + 6] & 7];
				flagblock[(bptr[j + 7] & 4194303) >> 3] &= masks[bptr[j + 7] & 7];

				// then compute their next hit and update the roots while they are
				// still fresh in the cache
				vt1 = _mm512_loadu_si512((__m512i *)(&bptr[j]));
				vt3 = _mm512_srli_epi64(vt1, 32);
				vt1 = _mm512_add_epi64(vt1, vt3);
				_mm512_storeu_si512((__m512i *)(&bptr[j]), vt1);
				cmp1 = _mm512_cmpgt_epu32_mask(vlinesize, vt1);

				if (cmp1 & 0x1) { BUCKET_UPDATE2(0, 22) }
				if (cmp1 & 0x4) { BUCKET_UPDATE2(1, 22) }
				if (cmp1 & 0x10) { BUCKET_UPDATE2(2, 22) }
				if (cmp1 & 0x40) { BUCKET_UPDATE2(3, 22) }
				if (cmp1 & 0x100) { BUCKET_UPDATE2(4, 22) }
				if (cmp1 & 0x400) { BUCKET_UPDATE2(5, 22) }
				if (cmp1 & 0x1000) { BUCKET_UPDATE2(6, 22) }
				if (cmp1 & 0x4000) { BUCKET_UPDATE2(7, 22) }

			}

			//finish up those that didn't fit into a group of 8 hits
			for (; j < nptr[i]; j++)
			{
				flagblock[(bptr[j] & 4194303) >> 3] &= masks[bptr[j] & 7];

				bptr[j] += (bptr[j] >> 32);
				if ((uint32_t)bptr[j] < linesize)
				{
					bnum = ((uint32_t)bptr[j] >> 22);
					buckets[bnum][nptr[bnum]] = bptr[j];
					nptr[bnum]++;
				}

			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
			{
				__m512i vbuck;
				__m512i vt1;
				__m512i vt2;
				__m512i v4194303 = _mm512_set1_epi32(4194303);
				__m512i v31 = _mm512_set1_epi32(31);
				__m512i vfull = _mm512_set1_epi32(0xffffffff);
				__m512i vone = _mm512_set1_epi32(1);

				uint32_t *large_bptr = ddata->large_sieve_buckets[i];
				uint32_t *large_nptr = ddata->large_bucket_hits;

				for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{
					// unload 16 hits
					// The AVX2 is almost exactly the same speed as the non-AVX2 code...
					// the bottleneck is not in computing the indices.
					// keep it here for future reference: scatter might help eventually.
					vbuck = _mm512_loadu_si512((__m256i *)(large_bptr + j));
					vt1 = _mm512_and_epi32(vbuck, v4194303);
					vt2 = _mm512_and_epi32(vt1, v31);
					vt2 = _mm512_sllv_epi32(vone, vt2);        // bit location
					vt2 = _mm512_andnot_epi32(vt2, vfull);      // sieve &= not(bit location)
					vt1 = _mm512_srli_epi32(vt1, 5);

					// this is not collision free... and it's not
					// much faster, if at all, from the other way.
					//vbuck = _mm512_i32gather_epi32(vt1, flagblock32, _MM_SCALE_4);
					//vbuck = _mm512_and_epi32(vext, vt2);
					//_mm512_i32scatter_epi32(flagblock32, vt1, vbuck, _MM_SCALE_4);

					_mm512_storeu_si512((__m256i *)t, vt1);
					_mm512_storeu_si512((__m256i *)t2, vt2);

					flagblock32[t[0]] &= t2[0];
					flagblock32[t[1]] &= t2[1];
					flagblock32[t[2]] &= t2[2];
					flagblock32[t[3]] &= t2[3];
					flagblock32[t[4]] &= t2[4];
					flagblock32[t[5]] &= t2[5];
					flagblock32[t[6]] &= t2[6];
					flagblock32[t[7]] &= t2[7];
					flagblock32[t[8]] &= t2[8];
					flagblock32[t[9]] &= t2[9];
					flagblock32[t[10]] &= t2[10];
					flagblock32[t[11]] &= t2[11];
					flagblock32[t[12]] &= t2[12];
					flagblock32[t[13]] &= t2[13];
					flagblock32[t[14]] &= t2[14];
					flagblock32[t[15]] &= t2[15];

				}

				for (; j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & 4194303) >> 3] &= masks[large_bptr[j] & 7];

			}

		}

		if (i == 0)
		{
			for (j = 0; j < 4194304; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					break;
				}
			}
		}

		flagblock += 524288;
	}

	return;
}

#endif

#if defined(USE_AVX2)
void sieve_line_avx2_32k(thread_soedata_t *thread_data)
{
	// extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = thread_data->sdata.lines[current_line];
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;
#ifdef BITMASKS32
    uint32_t* masks32 = sdata->masks32;
#endif

	// stuff for bucket sieving
	uint64_t *bptr;
	uint64_t **buckets;

	uint32_t *nptr;
	uint32_t linesize = 262144 * sdata->blocks, bnum;

	uint8_t *flagblock;
	uint64_t i, j, k;
	uint32_t prime;
	//uint32_t maxP;
	int stopid;

	//printf("lowlimit = %lu, range = %lu\n", sdata->lowlimit, sdata->blk_r);

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;

	mpz_t gmp_sqrt;
	mpz_init(gmp_sqrt);

	//ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	mpz_set_ui(gmp_sqrt, ddata->ublk_b + sdata->prodN);
	mpz_sqrt(gmp_sqrt, gmp_sqrt);
	ddata->blk_b_sqrt = mpz_get_ui(gmp_sqrt) + 1;

	mpz_clear(gmp_sqrt);

	// for the current line, find the offsets of each small/med prime past the low limit
	// and bucket sieve large primes
	get_offsets(thread_data);

	flagblock = line;
	for (i = 0; i < sdata->blocks; i++)
	{
		uint32_t *flagblock32 = (uint32_t *)flagblock;

		ALIGNED_MEM uint32_t t[8];

		if (sdata->num_bitmap_primes == 0)
		{
			// set all flags for this block, which also puts it into cache for the sieving
			// to follow.  only do this if we are not bitmap sieving, because
			// that happens before the linesieve and we don't want to overwrite it.
			memset(flagblock, 255, 32768);
		}

		// smallest primes use special methods
		pre_sieve_ptr(ddata, sdata, flagblock);

		// one is not a prime
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

		// start where presieving left off, which is different for various cpus.
        j = sdata->presieve_max_id;

		// at first glance it looks like AVX2 operations to compute the indices
		// might be helpful, but since we can't address memory with SIMD registers
		// it actually isn't.  Things might be different with a scatter operation.
		stopid = MIN(ddata->pbounds[i], 1901); //22998); // 1901);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			// we store byte masks to speed things up.  The byte masks
			// are addressed by index % 8.  However, (k+p) % 8 is
			// the same as (k+8p) % 8 so it suffices to compute and
			// store in registers (k+np) % 8 for n = 0:7.  Thanks to 
			// tverniquet for discovering this optimization.
			// The precomputed mask trick works well here because many
			// of these primes actually hit the interval many more
			// than 16 times, thus we get a lot of reuse out of
			// the masks.  In subsequent loops this isn't true.
			uint8_t m0;
			uint8_t m1;
			uint8_t m2;
			uint8_t m3;
			uint8_t m4;
			uint8_t m5;
			uint8_t m6;
			uint8_t m7;

			prime = sdata->sieve_p[j];

			tmpP = prime << 4;
			stop = 262144 - tmpP + prime;
			k = ddata->offsets[j];

			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			m0 = masks[k & 7];
			m1 = masks[(k + p1) & 7];
			m2 = masks[(k + p2) & 7];
			m3 = masks[(k + p3) & 7];
			m4 = masks[(k + 4 * prime) & 7];
			m5 = masks[(k + 5 * prime) & 7];
			m6 = masks[(k + 6 * prime) & 7];
			m7 = masks[(k + 7 * prime) & 7];

			while (k < stop)
			{
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
			}

            for (; k < 262144; k += prime)
            {
#ifdef BITMASKS8
                flagblock[k >> 3] &= masks[k & 7];
#elif defined(BITMASKS32)
                flagblock32[k >> 5] &= masks32[k & 31];
#elif defined(BITLOGIC8)
                flagblock[k >> 3] &= ~(1 << (k & 7));
#else
                flagblock32[k >> 5] &= ~(1 << (k & 31));
#endif
            }

			ddata->offsets[j] = k - 262144;

		}

		// unroll the loop: all primes less than this max hit the interval at least 8 times
		//maxP = 262144 >> 3;

		stopid = MIN(ddata->pbounds[i], 3513); //43388); // 3513);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 3;
			stop = 262144 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			while (k < stop)
			{
#ifdef BITMASKS8
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
#elif defined(BITMASKS32)
                flagblock32[k >> 5] &= masks32[k & 31];
                flagblock32[(k + p1) >> 5] &= masks32[(k + p1) & 31];
                flagblock32[(k + p2) >> 5] &= masks32[(k + p2) & 31];
                flagblock32[(k + p3) >> 5] &= masks32[(k + p3) & 31];
                k += (prime << 2);
                flagblock32[k >> 5] &= masks32[k & 31];
                flagblock32[(k + p1) >> 5] &= masks32[(k + p1) & 31];
                flagblock32[(k + p2) >> 5] &= masks32[(k + p2) & 31];
                flagblock32[(k + p3) >> 5] &= masks32[(k + p3) & 31];
                k += (prime << 2);
#elif defined(BITLOGIC8)
                flagblock[k >> 3] &= ~(1 << (k & 7));
                flagblock[(k + p1) >> 3] &= ~(1 << ((k + p1) & 7));
                flagblock[(k + p2) >> 3] &= ~(1 << ((k + p2) & 7));
                flagblock[(k + p3) >> 3] &= ~(1 << ((k + p3) & 7));
                k += (prime << 2);
                flagblock[k >> 3] &= ~(1 << (k & 7));
                flagblock[(k + p1) >> 3] &= ~(1 << ((k + p1) & 7));
                flagblock[(k + p2) >> 3] &= ~(1 << ((k + p2) & 7));
                flagblock[(k + p3) >> 3] &= ~(1 << ((k + p3) & 7));
                k += (prime << 2);
#else
                flagblock32[k >> 5] &= ~(1 << (k & 31));
                flagblock32[(k + p1) >> 5] &= ~(1 << ((k + p1) & 31));
                flagblock32[(k + p2) >> 5] &= ~(1 << ((k + p2) & 31));
                flagblock32[(k + p3) >> 5] &= ~(1 << ((k + p3) & 31));
                k += (prime << 2);
                flagblock32[k >> 5] &= ~(1 << (k & 31));
                flagblock32[(k + p1) >> 5] &= ~(1 << ((k + p1) & 31));
                flagblock32[(k + p2) >> 5] &= ~(1 << ((k + p2) & 31));
                flagblock32[(k + p3) >> 5] &= ~(1 << ((k + p3) & 31));
                k += (prime << 2);
#endif
			}

            for (; k < 262144; k += prime)
            {
#ifdef BITMASKS8
                flagblock[k >> 3] &= masks[k & 7];
#elif defined(BITMASKS32)
                flagblock32[k >> 5] &= masks32[k & 31];
#elif defined(BITLOGIC8)
                flagblock[k >> 3] &= ~(1 << (k & 7));
#else
                flagblock32[k >> 5] &= ~(1 << (k & 31));
#endif
            }

			ddata->offsets[j] = (uint32_t)(k - 262144);

		}

		// unroll the loop: all primes less than this max hit the interval at least 4 times
		//maxP = 262144 >> 2;

		stopid = MIN(ddata->pbounds[i], 6543); //82023); // 
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 2;
			stop = 262144 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{

#ifdef BITMASKS8
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
#elif defined(BITMASKS32)
                flagblock32[k >> 5] &= masks32[k & 31];
                flagblock32[(k + p1) >> 5] &= masks32[(k + p1) & 31];
                flagblock32[(k + p2) >> 5] &= masks32[(k + p2) & 31];
                flagblock32[(k + p3) >> 5] &= masks32[(k + p3) & 31];
                k += (prime << 2);
#elif defined(BITLOGIC8)
                flagblock[k >> 3] &= ~(1 << (k & 7));
                flagblock[(k + p1) >> 3] &= ~(1 << ((k + p1) & 7));
                flagblock[(k + p2) >> 3] &= ~(1 << ((k + p2) & 7));
                flagblock[(k + p3) >> 3] &= ~(1 << ((k + p3) & 7));
                k += (prime << 2);
#else
                flagblock32[k >> 5] &= ~(1 << (k & 31));
                flagblock32[(k + p1) >> 5] &= ~(1 << ((k + p1) & 31));
                flagblock32[(k + p2) >> 5] &= ~(1 << ((k + p2) & 31));
                flagblock32[(k + p3) >> 5] &= ~(1 << ((k + p3) & 31));
                k += (prime << 2);
#endif
			}

            for (; k < 262144; k += prime)
            {
#ifdef BITMASKS8
                flagblock[k >> 3] &= masks[k & 7];
#elif defined(BITMASKS32)
                flagblock32[k >> 5] &= masks32[k & 31];
#elif defined(BITLOGIC8)
                flagblock[k >> 3] &= ~(1 << (k & 7));
#else
                flagblock32[k >> 5] &= ~(1 << (k & 31));
#endif
            }

			ddata->offsets[j] = (uint32_t)(k - 262144);

		}

		// primes are getting fairly big now, unrolling is less useful.
		// keep going up to the large prime bound.
		stopid = ddata->pbounds[i];
        for (; j < stopid; j++)
        {
            prime = sdata->sieve_p[j];

            if (prime > 262144)
                break;

            for (k = ddata->offsets[j]; k < 262144; k += prime)
            {
#ifdef BITMASKS8
                flagblock[k >> 3] &= masks[k & 7];
#elif defined(BITMASKS32)
                flagblock32[k >> 5] &= masks32[k & 31];
#elif defined(BITLOGIC8)
                flagblock[k >> 3] &= ~(1 << (k & 7));
#else
                flagblock32[k >> 5] &= ~(1 << (k & 31));
#endif
            }

            ddata->offsets[j] = (uint32_t)(k - 262144);
        }

		for (; j < ddata->pbounds[i]; j++)
		{
			k = ddata->offsets[j];
			if (ddata->offsets[j] < 262144)
			{
#ifdef BITMASKS8
                flagblock[k >> 3] &= masks[k & 7];
#elif defined(BITMASKS32)
                flagblock32[k >> 5] &= masks32[k & 31];
#elif defined(BITLOGIC8)
                flagblock[k >> 3] &= ~(1 << (k & 7));
#else
                flagblock32[k >> 5] &= ~(1 << (k & 31));
#endif
				k += sdata->sieve_p[j];
			}
			ddata->offsets[j] = (uint32_t)(k - 262144);
		}

		// Now that we have sieved tiny, small, and medium-sized primes
		// it is time to dump in the buckets.  Larger primes were sorted
		// into buckets (one bucket per block per residue class) in the
		// get_offsets function in the file offsets.c.
		// This is where AVX2 provides some benefit.
		if (ddata->bucket_depth > 0)
		{
			__m256i vt1, vt2, vt3, vt4;      // temp vectors
			__m256i vlinesize = _mm256_set1_epi32(linesize);
			int cmp1, cmp2;

			// finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j = 0; j < (nptr[i] & (uint32_t)(~7)); j += 8)
			{
				// unload 8 hits

				flagblock[(bptr[j + 0] & 262143) >> 3] &= masks[bptr[j + 0] & 7];
				flagblock[(bptr[j + 1] & 262143) >> 3] &= masks[bptr[j + 1] & 7];
				flagblock[(bptr[j + 2] & 262143) >> 3] &= masks[bptr[j + 2] & 7];
				flagblock[(bptr[j + 3] & 262143) >> 3] &= masks[bptr[j + 3] & 7];
				flagblock[(bptr[j + 4] & 262143) >> 3] &= masks[bptr[j + 4] & 7];
				flagblock[(bptr[j + 5] & 262143) >> 3] &= masks[bptr[j + 5] & 7];
				flagblock[(bptr[j + 6] & 262143) >> 3] &= masks[bptr[j + 6] & 7];
				flagblock[(bptr[j + 7] & 262143) >> 3] &= masks[bptr[j + 7] & 7];

				// then compute their next hit and update the roots while they are
				// still fresh in the cache
				vt1 = _mm256_loadu_si256((__m256i *)(&bptr[j]));
				vt2 = _mm256_loadu_si256((__m256i *)(&bptr[j + 4]));
				vt3 = _mm256_srli_epi64(vt1, 32);
				vt4 = _mm256_srli_epi64(vt2, 32);
				vt1 = _mm256_add_epi64(vt1, vt3);
				vt2 = _mm256_add_epi64(vt2, vt4);
				_mm256_storeu_si256((__m256i *)(&bptr[j]), vt1);
				_mm256_storeu_si256((__m256i *)(&bptr[j + 4]), vt2);
				vt1 = _mm256_cmpgt_epi32(vlinesize, vt1);
				vt2 = _mm256_cmpgt_epi32(vlinesize, vt2);
				cmp1 = _mm256_movemask_epi8(vt1);
				cmp2 = _mm256_movemask_epi8(vt2);

				if (cmp1 & 0x1) { BUCKET_UPDATE2(0, 18) }
				if (cmp1 & 0x0100) { BUCKET_UPDATE2(1, 18) }
				if (cmp1 & 0x010000) { BUCKET_UPDATE2(2, 18) }
				if (cmp1 & 0x01000000) { BUCKET_UPDATE2(3, 18) }
				if (cmp2 & 0x01) { BUCKET_UPDATE2(4, 18) }
				if (cmp2 & 0x0100) { BUCKET_UPDATE2(5, 18) }
				if (cmp2 & 0x010000) { BUCKET_UPDATE2(6, 18) }
				if (cmp2 & 0x01000000) { BUCKET_UPDATE2(7, 18) }
			}

			//finish up those that didn't fit into a group of 8 hits
			for (; j < nptr[i]; j++)
			{
				flagblock[(bptr[j] & 262143) >> 3] &= masks[bptr[j] & 7];

				bptr[j] += (bptr[j] >> 32);
				if ((uint32_t)bptr[j] < linesize)
				{
					bnum = ((uint32_t)bptr[j] >> 18);
					buckets[bnum][nptr[bnum]] = bptr[j];
					nptr[bnum]++;
				}

			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			// Todo: see if it is any faster to not sort into block buckets and just
			// dump all primes into the whole line at once.  Maybe we can even sort them
			// first.
			if (ddata->largep_offset > 0)
			{
				__m256i vbuck;
				__m256i vt1;
				__m256i vt2;
				__m128i vext;
				__m256i v262143 = _mm256_set1_epi32(262143);
				__m256i v31 = _mm256_set1_epi32(31);
				__m256i vfull = _mm256_set1_epi32(0xffffffff);
				__m256i vone = _mm256_set1_epi32(1);
				uint32_t *large_bptr = ddata->large_sieve_buckets[i];
				uint32_t *large_nptr = ddata->large_bucket_hits;

				for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{
					// unload 16 hits
					// The AVX2 is almost exactly the same speed as the non-AVX2 code...
					// the bottleneck is not in computing the indices.
					// keep it here for future reference: scatter might help eventually.
					vbuck = _mm256_loadu_si256((__m256i *)(large_bptr + j));
					vt1 = _mm256_and_si256(vbuck, v262143);
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

					vbuck = _mm256_loadu_si256((__m256i *)(large_bptr + j + 8));
					vt1 = _mm256_and_si256(vbuck, v262143);
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

				}

				for (; j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & 262143) >> 3] &= masks[large_bptr[j] & 7];

			}

		}

		if (i == 0)
		{
			for (j = 0; j < 262144; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					break;
				}
			}
		}

		flagblock += 32768;
	}


	return;
}

void sieve_line_avx2_128k(thread_soedata_t *thread_data)
{
	// extract stuff from the thread data structure
	soe_dynamicdata_t *ddata = &thread_data->ddata;
	soe_staticdata_t *sdata = &thread_data->sdata;
	uint32_t current_line = thread_data->current_line;
	uint8_t *line = thread_data->sdata.lines[current_line];
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

	// stuff for bucket sieving
	uint64_t *bptr;
	uint64_t **buckets;

	uint32_t *nptr;
	uint32_t linesize = 1048576 * sdata->blocks, bnum;

	uint8_t *flagblock;
	uint64_t i, j, k;
	uint32_t prime;
	//uint32_t maxP;
	int stopid;

	ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
	ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
	ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

	// for the current line, find the offsets of each small/med prime past the low limit
	// and bucket sieve large primes
	get_offsets(thread_data);

	flagblock = line;
	for (i = 0; i < sdata->blocks; i++)
	{
		uint32_t *flagblock32 = (uint32_t *)flagblock;

		ALIGNED_MEM uint32_t t[8];

		if (sdata->num_bitmap_primes == 0)
		{
			// set all flags for this block, which also puts it into cache for the sieving
			// to follow.  only do this if we are not bitmap sieving, because
			// that happens before the linesieve and we don't want to overwrite it.
			memset(flagblock, 255, 131072);
		}

		// smallest primes use special methods
		pre_sieve_ptr(ddata, sdata, flagblock);

		// one is not a prime
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

		// start where presieving left off, which is different for various cpus.
		j = sdata->presieve_max_id;

		// unroll the loop: all primes less than this max hit the interval at least 16 times
		//maxP = 1048576 >> 4;

		// at first glance it looks like AVX2 operations to compute the indices
		// might be helpful, but since we can't address memory with SIMD registers
		// it actually isn't.  Things might be different with a scatter operation.
		stopid = MIN(ddata->pbounds[i], 6542); //22998); // 1901);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			// we store byte masks to speed things up.  The byte masks
			// are addressed by index % 8.  However, (k+p) % 8 is
			// the same as (k+8p) % 8 so it suffices to compute and
			// store in registers (k+np) % 8 for n = 0:7.  Thanks to 
			// tverniquet for discovering this optimization.
			// The precomputed mask trick works well here because many
			// of these primes actually hit the interval many more
			// than 16 times, thus we get a lot of reuse out of
			// the masks.  In subsequent loops this isn't true.
			uint8_t m0;
			uint8_t m1;
			uint8_t m2;
			uint8_t m3;
			uint8_t m4;
			uint8_t m5;
			uint8_t m6;
			uint8_t m7;

			prime = sdata->sieve_p[j];

			tmpP = prime << 4;
			stop = 1048576 - tmpP + prime;
			k = ddata->offsets[j];

			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			m0 = masks[k & 7];
			m1 = masks[(k + p1) & 7];
			m2 = masks[(k + p2) & 7];
			m3 = masks[(k + p3) & 7];
			m4 = masks[(k + 4 * prime) & 7];
			m5 = masks[(k + 5 * prime) & 7];
			m6 = masks[(k + 6 * prime) & 7];
			m7 = masks[(k + 7 * prime) & 7];

			while (k < stop)
			{
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
				flagblock[k >> 3] &= m0;
				flagblock[(k + p1) >> 3] &= m1;
				flagblock[(k + p2) >> 3] &= m2;
				flagblock[(k + p3) >> 3] &= m3;
				k += (prime << 2);
				flagblock[k >> 3] &= m4;
				flagblock[(k + p1) >> 3] &= m5;
				flagblock[(k + p2) >> 3] &= m6;
				flagblock[(k + p3) >> 3] &= m7;
				k += (prime << 2);
			}

			for (; k < 1048576; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = k - 1048576;

		}

		// unroll the loop: all primes less than this max hit the interval at least 8 times
		//maxP = 1048576 >> 3;

		stopid = MIN(ddata->pbounds[i], 12251); //43388); // 3513);
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 3;
			stop = 1048576 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;

			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 1048576; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 1048576);

		}

		// unroll the loop: all primes less than this max hit the interval at least 4 times
		//maxP = 1048576 >> 2;

		stopid = MIN(ddata->pbounds[i], 23000); //82023); // 6543
		for (; j < stopid; j++)
		{
			uint32_t tmpP;
			uint64_t stop;
			uint64_t p1, p2, p3;

			prime = sdata->sieve_p[j];

			tmpP = prime << 2;
			stop = 1048576 - tmpP + prime;
			k = ddata->offsets[j];
			p1 = prime;
			p2 = p1 + prime;
			p3 = p2 + prime;
			while (k < stop)
			{
				flagblock[k >> 3] &= masks[k & 7];
				flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
				flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
				flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
				k += (prime << 2);
			}

			for (; k < 1048576; k += prime)
				flagblock[k >> 3] &= masks[k & 7];

			ddata->offsets[j] = (uint32_t)(k - 1048576);

		}

		// primes are getting fairly big now, unrolling is less useful.
		// keep going up to the large prime bound.
		//stopid = MIN(ddata->pbounds[i], 23002);
		stopid = ddata->pbounds[i];
		for (; j < stopid; j++)
		{
			prime = sdata->sieve_p[j];

			if (prime > 1048576)
				break;

#ifdef _INTEL_COMPILER
#pragma nounroll
#endif
			for (k = ddata->offsets[j]; k < 1048576; k += prime)
			{
				flagblock[k >> 3] &= masks[k & 7];
			}

			ddata->offsets[j] = (uint32_t)(k - 1048576);
		}

		for (; j < ddata->pbounds[i]; j++)
		{
			k = ddata->offsets[j];
			if (ddata->offsets[j] < 1048576)
			{
				flagblock[k >> 3] &= masks[k & 7];
				k += sdata->sieve_p[j];
			}
			ddata->offsets[j] = (uint32_t)(k - 1048576);
		}

		if (ddata->bucket_depth > 0)
		{
			__m256i vt1, vt2, vt3, vt4;      // temp vectors
			__m256i vlinesize = _mm256_set1_epi32(linesize);
			int cmp1, cmp2;

			// finally, fill any primes in this block's bucket
			bptr = ddata->sieve_buckets[i];
			buckets = ddata->sieve_buckets;
			nptr = ddata->bucket_hits;

			//printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
			for (j = 0; j < (nptr[i] & (uint32_t)(~7)); j += 8)
			{
				// unload 8 hits

				flagblock[(bptr[j + 0] & 1048575) >> 3] &= masks[bptr[j + 0] & 7];
				flagblock[(bptr[j + 1] & 1048575) >> 3] &= masks[bptr[j + 1] & 7];
				flagblock[(bptr[j + 2] & 1048575) >> 3] &= masks[bptr[j + 2] & 7];
				flagblock[(bptr[j + 3] & 1048575) >> 3] &= masks[bptr[j + 3] & 7];
				flagblock[(bptr[j + 4] & 1048575) >> 3] &= masks[bptr[j + 4] & 7];
				flagblock[(bptr[j + 5] & 1048575) >> 3] &= masks[bptr[j + 5] & 7];
				flagblock[(bptr[j + 6] & 1048575) >> 3] &= masks[bptr[j + 6] & 7];
				flagblock[(bptr[j + 7] & 1048575) >> 3] &= masks[bptr[j + 7] & 7];

				// then compute their next hit and update the roots while they are
				// still fresh in the cache
				vt1 = _mm256_loadu_si256((__m256i *)(&bptr[j]));
				vt2 = _mm256_loadu_si256((__m256i *)(&bptr[j + 4]));
				vt3 = _mm256_srli_epi64(vt1, 32);
				vt4 = _mm256_srli_epi64(vt2, 32);
				vt1 = _mm256_add_epi64(vt1, vt3);
				vt2 = _mm256_add_epi64(vt2, vt4);
				_mm256_storeu_si256((__m256i *)(&bptr[j]), vt1);
				_mm256_storeu_si256((__m256i *)(&bptr[j + 4]), vt2);
				vt1 = _mm256_cmpgt_epi32(vlinesize, vt1);
				vt2 = _mm256_cmpgt_epi32(vlinesize, vt2);
				cmp1 = _mm256_movemask_epi8(vt1);
				cmp2 = _mm256_movemask_epi8(vt2);

				if (cmp1 & 0x1) { BUCKET_UPDATE2(0, 20) }
				if (cmp1 & 0x0100) { BUCKET_UPDATE2(1, 20) }
				if (cmp1 & 0x010000) { BUCKET_UPDATE2(2, 20) }
				if (cmp1 & 0x01000000) { BUCKET_UPDATE2(3, 20) }
				if (cmp2 & 0x01) { BUCKET_UPDATE2(4, 20) }
				if (cmp2 & 0x0100) { BUCKET_UPDATE2(5, 20) }
				if (cmp2 & 0x010000) { BUCKET_UPDATE2(6, 20) }
				if (cmp2 & 0x01000000) { BUCKET_UPDATE2(7, 20) }
			}

			//finish up those that didn't fit into a group of 8 hits
			for (; j < nptr[i]; j++)
			{
				flagblock[(bptr[j] & 1048575) >> 3] &= masks[bptr[j] & 7];

				bptr[j] += (bptr[j] >> 32);
				if ((uint32_t)bptr[j] < linesize)
				{
					bnum = ((uint32_t)bptr[j] >> 20);
					buckets[bnum][nptr[bnum]] = bptr[j];
					nptr[bnum]++;
				}

			}

			// repeat the dumping of bucket primes, this time with very large primes
			// that only hit the interval once.  thus, we don't need to update the root
			// with the next hit, and we can do more at once because each bucket hit is smaller
			if (ddata->largep_offset > 0)
			{
				__m256i vbuck;
				__m256i vt1;
				__m256i vt2;
				__m128i vext;
				__m256i v1048575 = _mm256_set1_epi32(1048575);
				__m256i v31 = _mm256_set1_epi32(31);
				__m256i vfull = _mm256_set1_epi32(0xffffffff);
				__m256i vone = _mm256_set1_epi32(1);
				uint32_t *large_bptr = ddata->large_sieve_buckets[i];
				uint32_t *large_nptr = ddata->large_bucket_hits;

				for (j = 0; j < (large_nptr[i] - 16); j += 16)
				{
					// unload 16 hits
					// The AVX2 is almost exactly the same speed as the non-AVX2 code...
					// the bottleneck is not in computing the indices.
					// keep it here for future reference: scatter might help eventually.
					vbuck = _mm256_loadu_si256((__m256i *)(large_bptr + j));
					vt1 = _mm256_and_si256(vbuck, v1048575);
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

					vbuck = _mm256_loadu_si256((__m256i *)(large_bptr + j + 8));
					vt1 = _mm256_and_si256(vbuck, v1048575);
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

				}

				for (; j < large_nptr[i]; j++)
					flagblock[(large_bptr[j] & 1048575) >> 3] &= masks[large_bptr[j] & 7];

			}

		}

		if (i == 0)
		{
			for (j = 0; j < 1048576; j++)
			{
				if (flagblock[j >> 3] & nmasks[j & 7])
				{
					ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
					break;
				}
			}
		}

		flagblock += 131072;
	}

	return;
}

void sieve_line_avx2_512k(thread_soedata_t* thread_data)
{
    // extract stuff from the thread data structure
    soe_dynamicdata_t* ddata = &thread_data->ddata;
    soe_staticdata_t* sdata = &thread_data->sdata;
    uint32_t current_line = thread_data->current_line;
    uint8_t* line = thread_data->sdata.lines[current_line];
    uint8_t* masks = sdata->masks;
    uint8_t* nmasks = sdata->nmasks;

    // stuff for bucket sieving
    uint64_t* bptr;
    uint64_t** buckets;

    uint32_t* nptr;
    uint32_t linesize = 4194304 * sdata->blocks, bnum;

    uint8_t* flagblock;
    uint64_t i, j, k;
    uint32_t prime;
    //uint32_t maxP;
    int stopid;

    ddata->lblk_b = sdata->lowlimit + sdata->rclass[current_line];
    ddata->ublk_b = sdata->blk_r + ddata->lblk_b - sdata->prodN;
    ddata->blk_b_sqrt = (sqrt(ddata->ublk_b + sdata->prodN)) + 1;

    // for the current line, find the offsets of each small/med prime past the low limit
    // and bucket sieve large primes
    get_offsets(thread_data);

    flagblock = line;
    for (i = 0; i < sdata->blocks; i++)
    {
        uint32_t* flagblock32 = (uint32_t*)flagblock;

        ALIGNED_MEM uint32_t t[8];

        if (sdata->num_bitmap_primes == 0)
        {
            // set all flags for this block, which also puts it into cache for the sieving
            // to follow.  only do this if we are not bitmap sieving, because
            // that happens before the linesieve and we don't want to overwrite it.
            memset(flagblock, 255, 524288);
        }

        // smallest primes use special methods
        pre_sieve_ptr(ddata, sdata, flagblock);

        // one is not a prime
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

        // start where presieving left off, which is different for various cpus.
        j = sdata->presieve_max_id;

        // unroll the loop: all primes less than this max hit the interval at least 16 times
        //maxP = 4194304 >> 4;

        // at first glance it looks like AVX2 operations to compute the indices
        // might be helpful, but since we can't address memory with SIMD registers
        // it actually isn't.  Things might be different with a scatter operation.
        stopid = MIN(ddata->pbounds[i], 23000); //22998); // 1901);
        for (; j < stopid; j++)
        {
            uint32_t tmpP;
            uint64_t stop;
            uint64_t p1, p2, p3;

            // we store byte masks to speed things up.  The byte masks
            // are addressed by index % 8.  However, (k+p) % 8 is
            // the same as (k+8p) % 8 so it suffices to compute and
            // store in registers (k+np) % 8 for n = 0:7.  Thanks to 
            // tverniquet for discovering this optimization.
            // The precomputed mask trick works well here because many
            // of these primes actually hit the interval many more
            // than 16 times, thus we get a lot of reuse out of
            // the masks.  In subsequent loops this isn't true.
            uint8_t m0;
            uint8_t m1;
            uint8_t m2;
            uint8_t m3;
            uint8_t m4;
            uint8_t m5;
            uint8_t m6;
            uint8_t m7;

            prime = sdata->sieve_p[j];

            tmpP = prime << 4;
            stop = 4194304 - tmpP + prime;
            k = ddata->offsets[j];

            p1 = prime;
            p2 = p1 + prime;
            p3 = p2 + prime;

            m0 = masks[k & 7];
            m1 = masks[(k + p1) & 7];
            m2 = masks[(k + p2) & 7];
            m3 = masks[(k + p3) & 7];
            m4 = masks[(k + 4 * prime) & 7];
            m5 = masks[(k + 5 * prime) & 7];
            m6 = masks[(k + 6 * prime) & 7];
            m7 = masks[(k + 7 * prime) & 7];

            while (k < stop)
            {
                flagblock[k >> 3] &= m0;
                flagblock[(k + p1) >> 3] &= m1;
                flagblock[(k + p2) >> 3] &= m2;
                flagblock[(k + p3) >> 3] &= m3;
                k += (prime << 2);
                flagblock[k >> 3] &= m4;
                flagblock[(k + p1) >> 3] &= m5;
                flagblock[(k + p2) >> 3] &= m6;
                flagblock[(k + p3) >> 3] &= m7;
                k += (prime << 2);
                flagblock[k >> 3] &= m0;
                flagblock[(k + p1) >> 3] &= m1;
                flagblock[(k + p2) >> 3] &= m2;
                flagblock[(k + p3) >> 3] &= m3;
                k += (prime << 2);
                flagblock[k >> 3] &= m4;
                flagblock[(k + p1) >> 3] &= m5;
                flagblock[(k + p2) >> 3] &= m6;
                flagblock[(k + p3) >> 3] &= m7;
                k += (prime << 2);
            }

            for (; k < 4194304; k += prime)
                flagblock[k >> 3] &= masks[k & 7];

            ddata->offsets[j] = k - 4194304;

        }

        // unroll the loop: all primes less than this max hit the interval at least 8 times
        //maxP = 4194304 >> 3;

        stopid = MIN(ddata->pbounds[i], 43390); //43388); // 3513);
        for (; j < stopid; j++)
        {
            uint32_t tmpP;
            uint64_t stop;
            uint64_t p1, p2, p3;

            prime = sdata->sieve_p[j];

            tmpP = prime << 3;
            stop = 4194304 - tmpP + prime;
            k = ddata->offsets[j];
            p1 = prime;
            p2 = p1 + prime;
            p3 = p2 + prime;

            while (k < stop)
            {
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
            }

            for (; k < 4194304; k += prime)
                flagblock[k >> 3] &= masks[k & 7];

            ddata->offsets[j] = (uint32_t)(k - 4194304);

        }

        // unroll the loop: all primes less than this max hit the interval at least 4 times
        //maxP = 4194304 >> 2;

        stopid = MIN(ddata->pbounds[i], 82025); //82023); // 6543
        for (; j < stopid; j++)
        {
            uint32_t tmpP;
            uint64_t stop;
            uint64_t p1, p2, p3;

            prime = sdata->sieve_p[j];

            tmpP = prime << 2;
            stop = 4194304 - tmpP + prime;
            k = ddata->offsets[j];
            p1 = prime;
            p2 = p1 + prime;
            p3 = p2 + prime;
            while (k < stop)
            {
                flagblock[k >> 3] &= masks[k & 7];
                flagblock[(k + p1) >> 3] &= masks[(k + p1) & 7];
                flagblock[(k + p2) >> 3] &= masks[(k + p2) & 7];
                flagblock[(k + p3) >> 3] &= masks[(k + p3) & 7];
                k += (prime << 2);
            }

            for (; k < 4194304; k += prime)
                flagblock[k >> 3] &= masks[k & 7];

            ddata->offsets[j] = (uint32_t)(k - 4194304);

        }

        // primes are getting fairly big now, unrolling is less useful.
        // keep going up to the large prime bound.
        //stopid = MIN(ddata->pbounds[i], 23002);
        stopid = ddata->pbounds[i];
        for (; j < stopid; j++)
        {
            prime = sdata->sieve_p[j];

            if (prime > 4194304)
                break;

#ifdef _INTEL_COMPILER
#pragma nounroll
#endif
            for (k = ddata->offsets[j]; k < 4194304; k += prime)
            {
                flagblock[k >> 3] &= masks[k & 7];
            }

            ddata->offsets[j] = (uint32_t)(k - 4194304);
        }

        for (; j < ddata->pbounds[i]; j++)
        {
            k = ddata->offsets[j];
            if (ddata->offsets[j] < 4194304)
            {
                flagblock[k >> 3] &= masks[k & 7];
                k += sdata->sieve_p[j];
            }
            ddata->offsets[j] = (uint32_t)(k - 4194304);
        }

        if (ddata->bucket_depth > 0)
        {
            __m256i vt1, vt2, vt3, vt4;      // temp vectors
            __m256i vlinesize = _mm256_set1_epi32(linesize);
            int cmp1, cmp2;

            // finally, fill any primes in this block's bucket
            bptr = ddata->sieve_buckets[i];
            buckets = ddata->sieve_buckets;
            nptr = ddata->bucket_hits;

            //printf("unloading %d hits in block %d of line %d\n",nptr[i],i,thread_data->current_line);
            for (j = 0; j < (nptr[i] & (uint32_t)(~7)); j += 8)
            {
                // unload 8 hits

                flagblock[(bptr[j + 0] & 4194303) >> 3] &= masks[bptr[j + 0] & 7];
                flagblock[(bptr[j + 1] & 4194303) >> 3] &= masks[bptr[j + 1] & 7];
                flagblock[(bptr[j + 2] & 4194303) >> 3] &= masks[bptr[j + 2] & 7];
                flagblock[(bptr[j + 3] & 4194303) >> 3] &= masks[bptr[j + 3] & 7];
                flagblock[(bptr[j + 4] & 4194303) >> 3] &= masks[bptr[j + 4] & 7];
                flagblock[(bptr[j + 5] & 4194303) >> 3] &= masks[bptr[j + 5] & 7];
                flagblock[(bptr[j + 6] & 4194303) >> 3] &= masks[bptr[j + 6] & 7];
                flagblock[(bptr[j + 7] & 4194303) >> 3] &= masks[bptr[j + 7] & 7];

                // then compute their next hit and update the roots while they are
                // still fresh in the cache
                vt1 = _mm256_loadu_si256((__m256i*)(&bptr[j]));
                vt2 = _mm256_loadu_si256((__m256i*)(&bptr[j + 4]));
                vt3 = _mm256_srli_epi64(vt1, 32);
                vt4 = _mm256_srli_epi64(vt2, 32);
                vt1 = _mm256_add_epi64(vt1, vt3);
                vt2 = _mm256_add_epi64(vt2, vt4);
                _mm256_storeu_si256((__m256i*)(&bptr[j]), vt1);
                _mm256_storeu_si256((__m256i*)(&bptr[j + 4]), vt2);
                vt1 = _mm256_cmpgt_epi32(vlinesize, vt1);
                vt2 = _mm256_cmpgt_epi32(vlinesize, vt2);
                cmp1 = _mm256_movemask_epi8(vt1);
                cmp2 = _mm256_movemask_epi8(vt2);

                if (cmp1 & 0x1) { BUCKET_UPDATE2(0, 22) }
                if (cmp1 & 0x0100) { BUCKET_UPDATE2(1, 22) }
                if (cmp1 & 0x010000) { BUCKET_UPDATE2(2, 22) }
                if (cmp1 & 0x01000000) { BUCKET_UPDATE2(3, 22) }
                if (cmp2 & 0x01) { BUCKET_UPDATE2(4, 22) }
                if (cmp2 & 0x0100) { BUCKET_UPDATE2(5, 22) }
                if (cmp2 & 0x010000) { BUCKET_UPDATE2(6, 22) }
                if (cmp2 & 0x01000000) { BUCKET_UPDATE2(7, 22) }
            }

            //finish up those that didn't fit into a group of 8 hits
            for (; j < nptr[i]; j++)
            {
                flagblock[(bptr[j] & 4194303) >> 3] &= masks[bptr[j] & 7];

                bptr[j] += (bptr[j] >> 32);
                if ((uint32_t)bptr[j] < linesize)
                {
                    bnum = ((uint32_t)bptr[j] >> 22);
                    buckets[bnum][nptr[bnum]] = bptr[j];
                    nptr[bnum]++;
                }

            }

            // repeat the dumping of bucket primes, this time with very large primes
            // that only hit the interval once.  thus, we don't need to update the root
            // with the next hit, and we can do more at once because each bucket hit is smaller
            if (ddata->largep_offset > 0)
            {
                __m256i vbuck;
                __m256i vt1;
                __m256i vt2;
                __m128i vext;
                __m256i v4194303 = _mm256_set1_epi32(4194303);
                __m256i v31 = _mm256_set1_epi32(31);
                __m256i vfull = _mm256_set1_epi32(0xffffffff);
                __m256i vone = _mm256_set1_epi32(1);
                uint32_t* large_bptr = ddata->large_sieve_buckets[i];
                uint32_t* large_nptr = ddata->large_bucket_hits;

                for (j = 0; j < (large_nptr[i] - 16); j += 16)
                {
                    // unload 16 hits
                    // The AVX2 is almost exactly the same speed as the non-AVX2 code...
                    // the bottleneck is not in computing the indices.
                    // keep it here for future reference: scatter might help eventually.
                    vbuck = _mm256_loadu_si256((__m256i*)(large_bptr + j));
                    vt1 = _mm256_and_si256(vbuck, v4194303);
                    vt2 = _mm256_and_si256(vt1, v31);
                    vt2 = _mm256_sllv_epi32(vone, vt2);        // bit location
                    vt2 = _mm256_andnot_si256(vt2, vfull);      // sieve &= not(bit location)
                    vt1 = _mm256_srai_epi32(vt1, 5);
                    _mm256_store_si256((__m256i*)t, vt1);
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

                    vbuck = _mm256_loadu_si256((__m256i*)(large_bptr + j + 8));
                    vt1 = _mm256_and_si256(vbuck, v4194303);
                    vt2 = _mm256_and_si256(vt1, v31);
                    vt2 = _mm256_sllv_epi32(vone, vt2);        // bit location
                    vt2 = _mm256_andnot_si256(vt2, vfull);      // sieve &= not(bit location)
                    vt1 = _mm256_srai_epi32(vt1, 5);

                    _mm256_store_si256((__m256i*)t, vt1);
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

                }

                for (; j < large_nptr[i]; j++)
                    flagblock[(large_bptr[j] & 4194303) >> 3] &= masks[large_bptr[j] & 7];

            }

        }

        if (i == 0)
        {
            for (j = 0; j < 4194304; j++)
            {
                if (flagblock[j >> 3] & nmasks[j & 7])
                {
                    ddata->min_sieved_val = sdata->prodN * j + sdata->rclass[current_line] + sdata->lowlimit;
                    break;
                }
            }
        }

        flagblock += 524288;
    }

    return;
}
#endif

