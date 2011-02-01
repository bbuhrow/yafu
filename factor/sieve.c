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

#if defined(_MSC_VER)
	#define ADDRESS_SIEVE(j) (sieve[j])
	#define ADDRESS_LOC(j) (bptr[j].loc)
#else
	#define ADDRESS_SIEVE(j) (*(sieve + (j)))
	#define ADDRESS_LOC(j) ((bptr + j)->loc)
#endif

void lp_sieveblock(uint8 *sieve, sieve_fb_compressed *fb, uint32 med_B, uint32 bnum, uint32 numblocks,
							 lp_bucket *lp, uint32 start_prime, uint8 s_init, int side)
{
	uint32 prime, root1, root2, tmp, stop, p16;
	uint32 i,j,lpnum,basebucket;
	uint8 logp;
	sieve_fb_compressed *fbptr;
	uint32 *bptr;
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	//initialize the block
	memset(sieve,s_init,BLOCKSIZE);

	p16 = BLOCKSIZE >> 1;
	for (i=start_prime;i<med_B;i++)
	{	
		uint8 *s2;
		fbptr = fb + i;
		prime = fbptr->prime_and_logp & 0xFFFF;
		root1 = fbptr->roots & 0xFFFF;
		root2 = fbptr->roots >> 16;
		logp = fbptr->prime_and_logp >> 16;

		//if we are past the blocksize, bail out because there are faster methods
		if (prime > p16)
			break;

		stop = BLOCKSIZE - prime;
		s2 = sieve + prime;

		while (root2 < stop)
		{
			sieve[root1] -= logp;
			sieve[root2] -= logp;
			s2[root1] -= logp;
			s2[root2] -= logp;
			root1 += (prime << 1);
			root2 += (prime << 1);
		}

		while (root2 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			ADDRESS_SIEVE(root2) -= logp;
			root1 += prime;
			root2 += prime;
		}

		//don't forget the last proot1[i], and compute the roots for the next block
		if (root1 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			root1 += prime;
			//root1 will be bigger on the next iteration, switch them now
			tmp = root2;
			root2 = root1;
			root1 = tmp;
		}
			
		fbptr->roots = (uint32)(((root2 - BLOCKSIZE) << 16) | (root1 - BLOCKSIZE));
	}

	for (;i<med_B;i++)
	{	
		fbptr = fb + i;
		prime = fbptr->prime_and_logp & 0xFFFF;
		root1 = fbptr->roots & 0xFFFF;
		root2 = fbptr->roots >> 16;
		logp = fbptr->prime_and_logp >> 16;

		//if we are past the blocksize, bail out because there are faster methods
		if (prime > BLOCKSIZE)
			break;

		while (root2 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			ADDRESS_SIEVE(root2) -= logp;
			root1 += prime;
			root2 += prime;
		}

		//don't forget the last proot1[i], and compute the roots for the next block
		if (root1 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			root1 += prime;
			//root1 will be bigger on the next iteration, switch them now
			tmp = root2;
			root2 = root1;
			root1 = tmp;
		}
			
		fbptr->roots = (uint32)(((root2 - BLOCKSIZE) << 16) | (root1 - BLOCKSIZE));
	}

	//if there are primes left bigger than the blocksize, this will take
	//care of them.  if not, it doesn't run at all.
	for (;i<med_B;i++)
	{	
		fbptr = fb + i;
		prime = fbptr->prime_and_logp & 0xFFFF;
		root1 = fbptr->roots & 0xFFFF;
		root2 = fbptr->roots >> 16;
		logp = fbptr->prime_and_logp >> 16;

		if (root1 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			root1 += prime;

			if (root2 < BLOCKSIZE)
			{
				ADDRESS_SIEVE(root2) -= logp;
				root2 += prime;
			}
			else
			{
				tmp=root2;
				root2=root1;
				root1=tmp;
			}
		}

		fbptr->roots = (uint32)(((root2 - BLOCKSIZE) << 16) | (root1 - BLOCKSIZE));
	}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		SIEVE_STG1 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

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

}

void test_block_siqs(uint8 *sieve, sieve_fb *fb, uint32 start_prime)
{
	//sieve the small primes over a block, in order to estimate their average
	//contribution to a sieve location.  there are probably analytical methods
	//to do this, but this is fast, easy, and gives decent results.
	uint32 prime, root1, root2;
	uint32 i,j;
	uint8 *inner_sieve;
	uint8 logp;
	sieve_fb *fbptr;
	
	//initialize block
	memset(sieve,0,BLOCKSIZE);

	inner_sieve = sieve;

	//break block into sections to keep in L1 cache
	//blocksize should be a multiple of inner_blocksize
	for (i=0;i<BLOCKSIZE;i+=INNER_BLOCKSIZE)
	{
		//for some small number of primes (those that fit in L1 cache/2
		for (j=2; j<start_prime; j++)
		{
			fbptr = fb + j;
			prime = fbptr->prime;
			root1 = fbptr->root1;
			root2 = fbptr->root2;
			logp = fbptr->logprime;

			while (root2 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				inner_sieve[root2] += logp;
				root1 += prime;
				root2 += prime;
			}

			//don't forget the last proot1[i], and compute the roots for the next block
			if (root1 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				root1 += prime;
				//don't update the root pointers... this is just a test.
			}
		}
		//move inner_sieve to the right to paint a new stripe of the small fb primes in sieve
		inner_sieve += INNER_BLOCKSIZE;
	}
	return;
}
