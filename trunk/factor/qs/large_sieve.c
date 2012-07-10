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

#include "yafu.h"
#include "qs.h"

#if defined(_MSC_VER)
	#include <mmintrin.h>
#endif

void lp_sieveblock(uint8 *sieve, uint32 bnum, uint32 numblocks,
		lp_bucket *lp, int side)
{

	uint32 i,j,lpnum,basebucket;
	uint8 logp;
	uint32 *bptr;

	//finally, dump the buckets into the now cached 
	//sieve block in prefetched batches
	//the starting address for this slice and bucket can be computed as an offset from
	//lp->list, the list of bucket hits.  each block is a bucket, and each bucket gets
	//2^BUCKET_BITS bytes of storage.  negative buckets are arranged immediately after
	//positive buckets, thus we have numblocks positive buckets followed by numblocks
	//negative buckets in every slice.  each slice is a complete set of such buckets.  slices
	//are contiguous in memory, so the whole thing is physically one giant linear list of
	//bucket hits which we have subdivided.  
	bptr = lp->list + (bnum << BUCKET_BITS);
	if (side)
	{
		bptr += (numblocks << BUCKET_BITS);
		basebucket = numblocks;
	}
	else
		basebucket = 0;

	//use x8 when cache line has 32 bytes
	//use x16 when chache line has 64 bytes
#if defined(CACHE_LINE_64)
	//printf("cache_line_64\n");
	for (j=0;j<lp->num_slices;j++)
	{
		lpnum = *(lp->num + bnum + basebucket);
		//printf("dumping %d primes from slice %d, bucket %d\n",lpnum, j, bnum);
		logp = *(lp->logp + j);
		for (i = 0; i < (lpnum & (uint32)(~15)); i += 16)
		{
#if defined(MANUAL_PREFETCH)
#if defined(GCC_ASM64X) || defined (__MINGW32__)
			//_mm_prefetch(bptr + i + 16,_MM_HINT_NTA);
			__asm ("prefetchnta %0	\n\t"
				: 
				:"m"(bptr + i + 16));
#else
			_mm_prefetch(bptr + i + 16,_MM_HINT_NTA);
#endif
#endif

			sieve[bptr[i  ] & 0x0000ffff] -= logp;
			sieve[bptr[i+1] & 0x0000ffff] -= logp;
			sieve[bptr[i+2] & 0x0000ffff] -= logp;
			sieve[bptr[i+3] & 0x0000ffff] -= logp;
			sieve[bptr[i+4] & 0x0000ffff] -= logp;
			sieve[bptr[i+5] & 0x0000ffff] -= logp;
			sieve[bptr[i+6] & 0x0000ffff] -= logp;
			sieve[bptr[i+7] & 0x0000ffff] -= logp;
			sieve[bptr[i+8] & 0x0000ffff] -= logp;
			sieve[bptr[i+9] & 0x0000ffff] -= logp;
			sieve[bptr[i+10] & 0x0000ffff] -= logp;
			sieve[bptr[i+11] & 0x0000ffff] -= logp;
			sieve[bptr[i+12] & 0x0000ffff] -= logp;
			sieve[bptr[i+13] & 0x0000ffff] -= logp;
			sieve[bptr[i+14] & 0x0000ffff] -= logp;
			sieve[bptr[i+15] & 0x0000ffff] -= logp;
		}

		for (;i<lpnum;i++)
			sieve[bptr[i] & 0x0000ffff] -= logp;
		
		//point to the next slice of primes
		bptr += (numblocks << (BUCKET_BITS + 1));
		basebucket += (numblocks << 1);
	}
#else

	for (j=0;j<lp->num_slices;j++)
	{
		lpnum = *(lp->num + bnum + basebucket);
		//printf("dumping %d primes from slice %d, bucket %d\n",lpnum, j, bnum);
		logp = *(lp->logp + j);
		for (i = 0; i < (lpnum & (uint32)(~7)); i += 8)
		//the slices can be considered stacks; the highest indices are put in last.
		//so it makes sense to take these out first, as they are more likely to 
		//still be in cache.
		//for (i = lpnum; i > (lpnum & 15); i -=16)
		{
			
#if defined(MANUAL_PREFETCH)
#if defined(GCC_ASM64X) || defined (__MINGW32__)
			//_mm_prefetch(bptr + i + 16,_MM_HINT_NTA);
			
#else
			//_mm_prefetch(bptr + i + 8,_MM_HINT_NTA);
#endif
#endif

			sieve[bptr[i  ] & 0x0000ffff] -= logp;
			sieve[bptr[i+1] & 0x0000ffff] -= logp;
			sieve[bptr[i+2] & 0x0000ffff] -= logp;
			sieve[bptr[i+3] & 0x0000ffff] -= logp;
			sieve[bptr[i+4] & 0x0000ffff] -= logp;
			sieve[bptr[i+5] & 0x0000ffff] -= logp;
			sieve[bptr[i+6] & 0x0000ffff] -= logp;
			sieve[bptr[i+7] & 0x0000ffff] -= logp;
		}

		for (;i<lpnum;i++)
			sieve[bptr[i] & 0x0000ffff] -= logp;
		
		//point to the next slice of primes
		bptr += (numblocks << (BUCKET_BITS + 1));
		basebucket += (numblocks << 1);
	}

#endif

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	SIEVE_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif

	return;
}

