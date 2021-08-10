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

#ifdef USE_XLBUCKET

    //ALIGNED_MEM uint32_t b1[2048];
    uint32_t sid;
    uint32_t shifted_block = bnum << 15;
    __m512i vmask15 = _mm512_set1_epi32(0x00007fff);

    // now dump in the xl-primes.  we have to check that
    // each hits the block first.  hoping this is less-bad
    // than sorting into blocks during poly update.
    // update: seems to be working, but this is slower than 
    // sorting into blocks during poly update.
    if (side == 0)
    {
        bptr = dconf->xl_pbucket.list; // +dconf->xl_pbucket.slicenum[bnum];

        //for (sid = 0; sid < dconf->xl_pbucket.numslices; sid++)
        {
            //printf("starting lpsieve on pslice %d, xl-prime id %u\n", sid, dconf->xl_pbucket.sliceid[sid]);
            
            //uint32_t numhit = 0;
            //lpnum = dconf->xl_pbucket.slicenum[sid];
            //logp = dconf->xl_pbucket.slicelogp[sid];

            //printf("compressing %u elements with logp = %u\n", lpnum, logp);

            //for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            //{
            //    __m512i v1 = _mm512_load_epi32(bptr + i);
            //    __mmask16 mask1 = _mm512_cmp_epu32_mask(
            //        _mm512_and_epi32(v1, _mm512_set1_epi32(0x001F8000)),
            //        _mm512_set1_epi32(shifted_block), _MM_CMPINT_EQ);
            //
            //    ALIGNED_MEM uint32_t b1[16];
            //
            //    //_mm512_mask_compressstoreu_epi32(b1, mask1, v1);
            //    _mm512_store_epi32(b1, v1);
            //
            //    uint32_t m1 = mask1;
            //
            //    while (m1 > 0)
            //    {
            //        int idx = _trail_zcnt(m1);
            //        sieve[bptr[i + idx] & 0x00007fff] -= logp;
            //        m1 = _reset_lsb(m1);
            //    }
            //
            //for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            //{
            //    __m512i v1 = _mm512_load_epi32(bptr + i);
            //    __mmask16 mask1 = _mm512_cmp_epu32_mask(
            //        _mm512_and_epi32(v1, _mm512_set1_epi32(0x001F8000)),
            //        _mm512_set1_epi32(shifted_block), _MM_CMPINT_EQ);
            //
            //    _mm512_mask_compressstoreu_epi32(b1 + numhit, mask1, _mm512_and_epi32(v1, vmask15));
            //    numhit += _mm_popcnt_u32(mask1);
            //}
            //
            //for (; i < lpnum; i++)
            //    if ((bptr[i] & 0x001F8000) == shifted_block) b1[numhit++] = bptr[i] & 0x00007fff;
            //
            //lpnum = numhit;
            //
            //printf("processing %u elements in block %d\n", lpnum, bnum);
            //for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            //{
            //    sieve[b1[i + 0] ]-= logp;
            //    sieve[b1[i + 1] ]-= logp;
            //    sieve[b1[i + 2] ]-= logp;
            //    sieve[b1[i + 3] ]-= logp;
            //    sieve[b1[i + 4] ]-= logp;
            //    sieve[b1[i + 5] ]-= logp;
            //    sieve[b1[i + 6] ]-= logp;
            //    sieve[b1[i + 7] ]-= logp;
            //    sieve[b1[i + 8] ]-= logp;
            //    sieve[b1[i + 9] ]-= logp;
            //    sieve[b1[i + 10]] -= logp;
            //    sieve[b1[i + 11]] -= logp;
            //    sieve[b1[i + 12]] -= logp;
            //    sieve[b1[i + 13]] -= logp;
            //    sieve[b1[i + 14]] -= logp;
            //    sieve[b1[i + 15]] -= logp;
                //if ((bptr[i + 0 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 0] & 0x00007fff] -= logp;
                //if ((bptr[i + 1 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 1] & 0x00007fff] -= logp;
                //if ((bptr[i + 2 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 2] & 0x00007fff] -= logp;
                //if ((bptr[i + 3 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 3] & 0x00007fff] -= logp;
                //if ((bptr[i + 4 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 4] & 0x00007fff] -= logp;
                //if ((bptr[i + 5 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 5] & 0x00007fff] -= logp;
                //if ((bptr[i + 6 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 6] & 0x00007fff] -= logp;
                //if ((bptr[i + 7 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 7] & 0x00007fff] -= logp;
                //if ((bptr[i + 8 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 8] & 0x00007fff] -= logp;
                //if ((bptr[i + 9 ] & 0x001F8000) == shifted_block) sieve[bptr[i + 9] & 0x00007fff] -= logp;
                //if ((bptr[i + 10] & 0x001F8000) == shifted_block) sieve[bptr[i + 10] & 0x00007fff] -= logp;
                //if ((bptr[i + 11] & 0x001F8000) == shifted_block) sieve[bptr[i + 11] & 0x00007fff] -= logp;
                //if ((bptr[i + 12] & 0x001F8000) == shifted_block) sieve[bptr[i + 12] & 0x00007fff] -= logp;
                //if ((bptr[i + 13] & 0x001F8000) == shifted_block) sieve[bptr[i + 13] & 0x00007fff] -= logp;
                //if ((bptr[i + 14] & 0x001F8000) == shifted_block) sieve[bptr[i + 14] & 0x00007fff] -= logp;
                //if ((bptr[i + 15] & 0x001F8000) == shifted_block) sieve[bptr[i + 15] & 0x00007fff] -= logp;
            //}

            //for (; i < lpnum; i++)
            //    sieve[b1[i]] -= logp;

            logp = 21;
            lpnum = dconf->xl_pbucket.sliceid[0]; // [bnum] ;
            //printf("dumping %u xlp entries in block %d side %d\n", lpnum, bnum, side);

            for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            {
                //__m512i v1 = _mm512_loadu_epi32(bptr + i);
                //ALIGNED_MEM uint32_t locs[16];
                //
                //_mm512_storeu_epi32(locs, _mm512_and_epi32(v1, vmask15));
                //
                //sieve[locs[0]] -= 21;
                //sieve[locs[1]] -= 21;
                //sieve[locs[2]] -= 21;
                //sieve[locs[3]] -= 21;
                //sieve[locs[4]] -= 21;
                //sieve[locs[5]] -= 21;
                //sieve[locs[6]] -= 21;
                //sieve[locs[7]] -= 21;
                //sieve[locs[8]] -= 21;
                //sieve[locs[9]] -= 21;
                //sieve[locs[10]] -= 21;
                //sieve[locs[11]] -= 21;
                //sieve[locs[12]] -= 21;
                //sieve[locs[13]] -= 21;
                //sieve[locs[14]] -= 21;
                //sieve[locs[15]] -= 21;

                __m512i v1 = _mm512_load_epi32(bptr + i);
                __mmask16 mask1 = _mm512_cmp_epu32_mask(
                    _mm512_and_epi32(v1, _mm512_set1_epi32(0x00038000)),
                    _mm512_set1_epi32(shifted_block), _MM_CMPINT_EQ);
                
                ALIGNED_MEM uint32_t b1[16];
                
                //_mm512_mask_compressstoreu_epi32(b1, mask1, v1);
                _mm512_store_epi32(b1, v1);
                
                uint32_t m1 = mask1;
                
                while (m1 > 0)
                {
                    int idx = _trail_zcnt(m1);
                    sieve[bptr[i + idx] & 0x00007fff] -= logp;
                    m1 = _reset_lsb(m1);
                }

            }

            //for (; i < lpnum; i++)
            //    sieve[bptr[i] & 0x00007fff] -= 21;

            for (; i < lpnum; i++)
                if ((bptr[i] & 0x0038000) == shifted_block) sieve[bptr[i] & 0x00007fff] -= logp;
        }
    }
    else
    {
        bptr = dconf->xl_nbucket.list; // +dconf->xl_nbucket.slicenum[bnum];

        //for (sid = 0; sid < dconf->xl_nbucket.numslices; sid++)
        {
            //printf("starting lpsieve on nslice %d, xl-prime id %u\n", sid, dconf->xl_nbucket.sliceid[sid]);
            //uint32_t numhit = 0;
            //lpnum = dconf->xl_nbucket.slicenum[sid];
            //logp = dconf->xl_nbucket.slicelogp[sid];

            //printf("compressing %u elements with logp = %u\n", lpnum, logp);

            
            //for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            //{
            //    __m512i v1 = _mm512_load_epi32(bptr + i);
            //    __mmask16 mask1 = _mm512_cmp_epu32_mask(
            //        _mm512_and_epi32(v1, _mm512_set1_epi32(0x001F8000)),
            //        _mm512_set1_epi32(shifted_block), _MM_CMPINT_EQ);
            //
            //    ALIGNED_MEM uint32_t b1[16];
            //
            //    //_mm512_mask_compressstoreu_epi32(b1, mask1, v1);
            //    _mm512_store_epi32(b1, v1);
            //
            //    uint32_t m1 = mask1;
            //
            //    while (m1 > 0)
            //    {
            //        int idx = _trail_zcnt(m1);
            //        sieve[bptr[i + idx] & 0x00007fff] -= logp;
            //        m1 = _reset_lsb(m1);
            //    }
            //
            //for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            //{
            //    __m512i v1 = _mm512_load_epi32(bptr + i);
            //    __mmask16 mask1 = _mm512_cmp_epu32_mask(
            //        _mm512_and_epi32(v1, _mm512_set1_epi32(0x001F8000)),
            //        _mm512_set1_epi32(shifted_block), _MM_CMPINT_EQ);
            //
            //    _mm512_mask_compressstoreu_epi32(b1 + numhit, mask1, _mm512_and_epi32(v1, vmask15));
            //    numhit += _mm_popcnt_u32(mask1);
            //}
            //
            //for (; i < lpnum; i++)
            //    if ((bptr[i] & 0x001F8000) == shifted_block) b1[numhit++] = bptr[i] & 0x00007fff;
            //
            //lpnum = numhit;
            //
            ////printf("processing %u elements with logp = %u\n", lpnum, logp);
            //printf("processing %u elements in block %u\n", lpnum, bnum);
            //for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            //{
            //    sieve[b1[i + 0]] -= logp;
            //    sieve[b1[i + 1]] -= logp;
            //    sieve[b1[i + 2]] -= logp;
            //    sieve[b1[i + 3]] -= logp;
            //    sieve[b1[i + 4]] -= logp;
            //    sieve[b1[i + 5]] -= logp;
            //    sieve[b1[i + 6]] -= logp;
            //    sieve[b1[i + 7]] -= logp;
            //    sieve[b1[i + 8]] -= logp;
            //    sieve[b1[i + 9]] -= logp;
            //    sieve[b1[i + 10]] -= logp;
            //    sieve[b1[i + 11]] -= logp;
            //    sieve[b1[i + 12]] -= logp;
            //    sieve[b1[i + 13]] -= logp;
            //    sieve[b1[i + 14]] -= logp;
            //    sieve[b1[i + 15]] -= logp;

                //if ((bptr[i + 0] & 0x001F8000) == shifted_block) sieve[bptr[i + 0] & 0x00007fff] -= logp;
                //if ((bptr[i + 1] & 0x001F8000) == shifted_block) sieve[bptr[i + 1] & 0x00007fff] -= logp;
                //if ((bptr[i + 2] & 0x001F8000) == shifted_block) sieve[bptr[i + 2] & 0x00007fff] -= logp;
                //if ((bptr[i + 3] & 0x001F8000) == shifted_block) sieve[bptr[i + 3] & 0x00007fff] -= logp;
                //if ((bptr[i + 4] & 0x001F8000) == shifted_block) sieve[bptr[i + 4] & 0x00007fff] -= logp;
                //if ((bptr[i + 5] & 0x001F8000) == shifted_block) sieve[bptr[i + 5] & 0x00007fff] -= logp;
                //if ((bptr[i + 6] & 0x001F8000) == shifted_block) sieve[bptr[i + 6] & 0x00007fff] -= logp;
                //if ((bptr[i + 7] & 0x001F8000) == shifted_block) sieve[bptr[i + 7] & 0x00007fff] -= logp;
                //if ((bptr[i + 8] & 0x001F8000) == shifted_block) sieve[bptr[i + 8] & 0x00007fff] -= logp;
                //if ((bptr[i + 9] & 0x001F8000) == shifted_block) sieve[bptr[i + 9] & 0x00007fff] -= logp;
                //if ((bptr[i + 10] & 0x001F8000) == shifted_block) sieve[bptr[i + 10] & 0x00007fff] -= logp;
                //if ((bptr[i + 11] & 0x001F8000) == shifted_block) sieve[bptr[i + 11] & 0x00007fff] -= logp;
                //if ((bptr[i + 12] & 0x001F8000) == shifted_block) sieve[bptr[i + 12] & 0x00007fff] -= logp;
                //if ((bptr[i + 13] & 0x001F8000) == shifted_block) sieve[bptr[i + 13] & 0x00007fff] -= logp;
                //if ((bptr[i + 14] & 0x001F8000) == shifted_block) sieve[bptr[i + 14] & 0x00007fff] -= logp;
                //if ((bptr[i + 15] & 0x001F8000) == shifted_block) sieve[bptr[i + 15] & 0x00007fff] -= logp;
            //}

            //for (; i < lpnum; i++)
            //    sieve[b1[i]] -= logp;

            logp = 21;
            lpnum = dconf->xl_nbucket.sliceid[0]; // [bnum] ;
            //printf("dumping %u xlp entries in block %d side %d\n", lpnum, bnum, side);

            for (i = 0; (uint32_t)i < (lpnum & (uint32_t)(~15)); i += 16)
            {
                //__m512i v1 = _mm512_loadu_epi32(bptr + i);
                //ALIGNED_MEM uint32_t locs[16];
                //
                //_mm512_storeu_epi32(locs, _mm512_and_epi32(v1, vmask15));
                //
                //sieve[locs[0]] -= 21;
                //sieve[locs[1]] -= 21;
                //sieve[locs[2]] -= 21;
                //sieve[locs[3]] -= 21;
                //sieve[locs[4]] -= 21;
                //sieve[locs[5]] -= 21;
                //sieve[locs[6]] -= 21;
                //sieve[locs[7]] -= 21;
                //sieve[locs[8]] -= 21;
                //sieve[locs[9]] -= 21;
                //sieve[locs[10]] -= 21;
                //sieve[locs[11]] -= 21;
                //sieve[locs[12]] -= 21;
                //sieve[locs[13]] -= 21;
                //sieve[locs[14]] -= 21;
                //sieve[locs[15]] -= 21;

                __m512i v1 = _mm512_load_epi32(bptr + i);
                __mmask16 mask1 = _mm512_cmp_epu32_mask(
                    _mm512_and_epi32(v1, _mm512_set1_epi32(0x00038000)),
                    _mm512_set1_epi32(shifted_block), _MM_CMPINT_EQ);

                ALIGNED_MEM uint32_t b1[16];

                //_mm512_mask_compressstoreu_epi32(b1, mask1, v1);
                _mm512_store_epi32(b1, v1);

                uint32_t m1 = mask1;

                while (m1 > 0)
                {
                    int idx = _trail_zcnt(m1);
                    sieve[bptr[i + idx] & 0x00007fff] -= logp;
                    m1 = _reset_lsb(m1);
                }
            }

            //for (; i < lpnum; i++)
            //    sieve[bptr[i] & 0x00007fff] -= 21;

            for (; i < lpnum; i++)
                if ((bptr[i] & 0x0038000) == shifted_block) sieve[bptr[i] & 0x00007fff] -= logp;
        }
    }

    

#endif

    return;
}
#endif