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

#include "qs_impl.h"

#if defined(_MSC_VER)
	#include <mmintrin.h>
#endif


#ifdef USE_AVX512F
#include <immintrin.h>
#define NUM_LANES 16
#endif


void lp_sieveblock(uint8_t *sieve, uint32_t bnum, uint32_t numblocks,
		lp_bucket *lp, int side, dynamic_conf_t * dconf)
{

	uint32_t i,j,lpnum,basebucket;
	uint8_t logp;
	uint32_t *bptr;

#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    int poly_offset = (dconf->numB % dconf->poly_batchsize) - 2;
    int pnum;

    if (dconf->numB == 1)
    {
        poly_offset = 0;
    }
    else if (poly_offset < 0)
    {
        poly_offset += dconf->poly_batchsize;
    }
    pnum = poly_offset;
    poly_offset = poly_offset * 2 * numblocks * dconf->buckets->alloc_slices;

    //printf("begin large_sieve on side %d with poly %d\n", side, pnum);
    
#else 
    int pnum = 0;
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

#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    bptr = lp->list + (bnum << BUCKET_BITS) + poly_offset * BUCKET_ALLOC;
#else
	bptr = lp->list + (bnum << BUCKET_BITS);
#endif

    if (side)
    {
        bptr += (numblocks << BUCKET_BITS);
        basebucket = numblocks;
    }
    else
    {
        basebucket = 0;
    }

	//use x8 when cache line has 32 bytes
	//use x16 when chache line has 64 bytes
#if defined(USE_BATCHPOLY_X2)
    for (j = 0; j < dconf->buckets->num_slices_batch[pnum]; j++)
#else
	for (j = 0; j < lp->num_slices; j++)
#endif
	{
#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
        lpnum = *(lp->num + bnum + basebucket + poly_offset);
#else
		lpnum = *(lp->num + bnum + basebucket);
        //lpnum = bptr[0];
#endif
		logp = *(lp->logp + j);

        for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
        {

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

	return;
}

#ifdef USE_AVX512F
void lp_sieveblock_avx512f(uint8_t* sieve, uint32_t bnum, uint32_t numblocks,
    lp_bucket* lp, int side, dynamic_conf_t* dconf)
{

    uint32_t i, j, lpnum, basebucket;
    uint8_t logp;
    uint32_t* bptr;

    __m512i vmask = _mm512_set1_epi32(0x0000ffff);
    __m512i vlomask = _mm512_set1_epi32(0x000000ff);
    __m512i vhimask = _mm512_set1_epi32(0xffffff00);
    __m512i vbuckets, vnextbuckets, vlosieve, vhisieve;
    __m512i vblock = _mm512_set1_epi32(BLOCKSIZE);
    __mmask16 mask16;


#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    int poly_offset = (dconf->numB % dconf->poly_batchsize) - 2;
    int pnum;

    if (dconf->numB == 1)
    {
        poly_offset = 0;
    }
    else if (poly_offset < 0)
    {
        poly_offset += dconf->poly_batchsize;
    }
    pnum = poly_offset;
    poly_offset = poly_offset * 2 * numblocks * dconf->buckets->alloc_slices;

    //printf("begin large_sieve on side %d with poly %d\n", side, pnum);

#else 
    int pnum = 0;
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

#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    bptr = lp->list + (bnum << BUCKET_BITS) + poly_offset * BUCKET_ALLOC;
#else
    bptr = lp->list + (bnum << BUCKET_BITS);
#endif

    if (side)
    {
        bptr += (numblocks << BUCKET_BITS);
        basebucket = numblocks;
    }
    else
    {
        basebucket = 0;
    }

    //use x8 when cache line has 32 bytes
    //use x16 when chache line has 64 bytes
#if defined(USE_BATCHPOLY_X2)
    for (j = 0; j < dconf->buckets->num_slices_batch[pnum]; j++)
#else
    for (j = 0; j < lp->num_slices; j++)
#endif
    {
#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
        lpnum = *(lp->num + bnum + basebucket + poly_offset);
#else
        lpnum = *(lp->num + bnum + basebucket);
        //lpnum = bptr[0];
#endif
        logp = *(lp->logp + j);

#ifdef DO_VLP_OPT
        if (lp->fb_bounds[j] & 0x80000000)
        {
            //lp->fb_bounds[j] ^= 0x80000000;
            //printf("found flag in lpsieve, commencing vlp sieve\n");
            break;
        }
#endif


        //#if defined(USE_BATCHPOLY_X2)
        //        printf("lp sieve: dumping %d primes from slice %d of %d for poly %d, bucket/block %d\n", 
        //            lpnum, j, dconf->buckets->num_slices_batch[pnum], pnum, bnum);
        //#else
        //        printf("lp sieve: dumping %d primes from slice %d of %d for poly %d, bucket/block %d\n",
        //            lpnum, j, dconf->buckets->num_slices, pnum, bnum);
        //#endif

        __m512i vlogp = _mm512_set1_epi32(logp);
        vnextbuckets = _mm512_and_epi32(vmask, _mm512_load_epi32((__m512i*)(&bptr[0])));

        for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
        {
            vbuckets = vnextbuckets;
            if ((i + 16) < lpnum)
            {
                vnextbuckets = _mm512_and_epi32(vmask, _mm512_load_epi32((__m512i*)(&bptr[i + 16])));
#ifdef USE_AVX512PF
                _mm512_prefetch_i32scatter_ps(sieve, vnextbuckets, _MM_SCALE_1, _MM_HINT_T0);
#endif
            }

            // ignore conflicts...
            vhisieve = _mm512_i32gather_epi32(vbuckets, sieve, _MM_SCALE_1);
            vlosieve = _mm512_and_epi32(vhisieve, vlomask);
            vhisieve = _mm512_and_epi32(vhisieve, vhimask);
            vlosieve = _mm512_sub_epi32(vlosieve, vlogp);
            vlosieve = _mm512_or_epi32(vhisieve, _mm512_and_epi32(vlosieve, vlomask));
            _mm512_i32scatter_epi32(sieve, vbuckets, vlosieve, _MM_SCALE_1);

        }

        for (; i < lpnum; i++)
            sieve[bptr[i] & 0x0000ffff] -= logp;

        //point to the next slice of primes
        bptr += (numblocks << (BUCKET_BITS + 1));
        basebucket += (numblocks << 1);
    }

    return;
}
#endif

#ifdef USE_AVX512BW
void lp_sieveblock_avx512bw(uint8_t* sieve, uint32_t bnum, uint32_t numblocks,
    lp_bucket* lp, int side, dynamic_conf_t* dconf)
{

    uint32_t i, j, lpnum, basebucket;
    uint8_t logp;
    uint32_t* bptr;

    __m512i vmask = _mm512_set1_epi32(0x0000ffff);
    __m512i vlomask = _mm512_set1_epi32(0x000000ff);
    __m512i vhimask = _mm512_set1_epi32(0xffffff00);
    __m512i vbuckets, vnextbuckets, vlosieve, vhisieve;
    __m512i vblock = _mm512_set1_epi32(BLOCKSIZE);
    __mmask16 mask16;

#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    int poly_offset = (dconf->numB % dconf->poly_batchsize) - 2;
    int pnum;

    if (dconf->numB == 1)
    {
        poly_offset = 0;
    }
    else if (poly_offset < 0)
    {
        poly_offset += dconf->poly_batchsize;
    }
    pnum = poly_offset;
    poly_offset = poly_offset * 2 * numblocks * dconf->buckets->alloc_slices;

    //printf("begin large_sieve on side %d with poly %d\n", side, pnum);

#else 
    int pnum = 0;
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

#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
    bptr = lp->list + (bnum << BUCKET_BITS) + poly_offset * BUCKET_ALLOC;
#else
    bptr = lp->list + (bnum << BUCKET_BITS);
#endif

    if (side)
    {
        bptr += (numblocks << BUCKET_BITS);
        basebucket = numblocks;
    }
    else
    {
        basebucket = 0;
    }

    //use x8 when cache line has 32 bytes
    //use x16 when chache line has 64 bytes
#if defined(USE_BATCHPOLY_X2)
    for (j = 0; j < dconf->buckets->num_slices_batch[pnum]; j++)
#else
    for (j = 0; j < lp->num_slices; j++)
#endif
    {
#if defined( USE_BATCHPOLY) || defined(USE_BATCHPOLY_X2)
        lpnum = *(lp->num + bnum + basebucket + poly_offset);
#else
        lpnum = *(lp->num + bnum + basebucket);
        //lpnum = bptr[0];
#endif
        logp = *(lp->logp + j);

#ifdef DO_VLP_OPT
        if (lp->fb_bounds[j] & 0x80000000)
        {
            //lp->fb_bounds[j] ^= 0x80000000;
            //printf("found flag in lpsieve, commencing vlp sieve\n");
            break;
        }
#endif


        //#if defined(USE_BATCHPOLY_X2)
        //        printf("lp sieve: dumping %d primes from slice %d of %d for poly %d, bucket/block %d\n", 
        //            lpnum, j, dconf->buckets->num_slices_batch[pnum], pnum, bnum);
        //#else
        //        printf("lp sieve: dumping %d primes from slice %d of %d for poly %d, bucket/block %d\n",
        //            lpnum, j, dconf->buckets->num_slices, pnum, bnum);
        //#endif

        __m512i vlogp = _mm512_set1_epi32(logp);
        vnextbuckets = _mm512_and_epi32(vmask, _mm512_load_epi32((__m512i*)(&bptr[0])));

        for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
        {
            vbuckets = vnextbuckets;
            if ((i + 16) < lpnum)
            {
                vnextbuckets = _mm512_and_epi32(vmask, _mm512_load_epi32((__m512i*)(&bptr[i + 16])));
            }

            // ignore conflicts...
            vhisieve = _mm512_i32gather_epi32(vbuckets, sieve, _MM_SCALE_1);
            vlosieve = _mm512_sub_epi8(vhisieve, vlogp);
            _mm512_i32scatter_epi32(sieve, vbuckets, vlosieve, _MM_SCALE_1);
        }

        for (; i < lpnum; i++)
            sieve[bptr[i] & 0x0000ffff] -= logp;

        //point to the next slice of primes
        bptr += (numblocks << (BUCKET_BITS + 1));
        basebucket += (numblocks << 1);
    }

    return;
}
#endif