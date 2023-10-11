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

#ifdef DEBUGPRINT_BATCHPOLY
    printf("begin large_sieve on side %d with poly %d\n", side, pnum);
#endif

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

#ifdef DEBUGPRINT_BATCHPOLY
    printf("begin large_sieve_avx512f on side %d with poly %d\n", side, pnum);
#endif

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

#ifdef DEBUGPRINT_BATCHPOLY
    printf("begin large_sieve_avx512bw on side %d with poly %d... ", side, pnum);
#endif

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

#ifdef DEBUGPRINT_BATCHPOLY
    printf("complete.\n"); fflush(stdout);
#endif
    return;
}
#endif

#ifdef USE_LINKED_LIST_SS
void lp_sieve_interval_ss(uint8_t* sieve, int side, dynamic_conf_t* dconf)
{
    // with linked links by poly
    int i;
    int j;
    //uint32_t numentries = 0;

    if (side == 0)
    {
        // find the first hit for this poly.  should probably instead store
        // a list of first locations for each poly.
        uint32_t eidx = 0;
        //while ((dconf->ss_slices_p[0].element[eidx] & 0xffff) != dconf->numB)
        //{
        //    eidx++;
        //}
        eidx = dconf->poly_ll_first_p[dconf->numB];

        //printf("first element index for poly %u side %d = %u\n", dconf->numB, side, eidx);
        uint32_t root;

        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint64_t* elements = dconf->ss_slices_p[i].element;
            uint8_t logp = dconf->ss_slices_p[i].logp;
            uint32_t nextentry;
            uint32_t slicesize = dconf->ss_slices_p[i].size;
            //__m512i vbuckets, vnextbuckets, vlosieve, vhisieve;
            //__m512i vlogp = _mm512_set1_epi32(logp);
            //ALIGNED_MEM uint32_t roots[16];
            //int j;

            do
            {
                //if ((dconf->ss_slices_p[i].element[eidx] & 0xffff) != dconf->numB)
                //{
                //    printf("error: next element index %u is not the same poly index (now %u)!\n",
                //        eidx, (dconf->ss_slices_p[i].element[eidx] & 0xffff));
                //    printf("faulty entry = %lx\n", dconf->ss_slices_p[i].element[eidx - nextentry]);
                //    printf("current entry = %lx\n", dconf->ss_slices_p[i].element[eidx]);
                //    printf("correctly processed %u entries for poly %d\n",
                //        numentries, dconf->numB);
                //    exit(1);
                //}

                //__mmask16 m = 0;
                //for (j = 0; j < 16; j++)
                //{
                //    nextentry = elements[eidx] >> 48;
                //    _mm_prefetch(&elements[eidx + nextentry], _MM_HINT_T1);
                //    roots[j] = (elements[eidx] & 0xffff0000) >> 16;
                //    m |= (1 << j);
                //    eidx += nextentry;
                //    if ((eidx >= slicesize) || (nextentry == 0xffff))
                //    {
                //        break;
                //    }
                //}
                //
                //vbuckets = _mm512_load_epi32(roots);
                //vhisieve = _mm512_mask_i32gather_epi32(vbuckets, m, vbuckets, sieve, _MM_SCALE_1);
                //vlosieve = _mm512_sub_epi8(vhisieve, vlogp);
                //_mm512_mask_i32scatter_epi32(sieve, m, vbuckets, vlosieve, _MM_SCALE_1);

                nextentry = elements[eidx] >> 40;
                _mm_prefetch(&elements[eidx + nextentry], _MM_HINT_T1);

                root = (elements[eidx] & 0xffffff);
                sieve[root] += logp;
                eidx += nextentry;

                //numentries++;
            } while ((eidx < slicesize) && (nextentry != 0xffffff));




            if (nextentry == 0xffffff)
            {
                // last entry for this poly
                break;
            }
            else
            {
                // more entries in the next slice.  adjust the element index
                // for the remainder of this slice's size.
                //printf("finished processing %d entries from slice %d, adjusting "
                //    "next entry index %d by size of slice %d (%d) to %d\n",
                //    numentries, i, eidx, i, dconf->ss_slices_p[i].size,
                //    eidx - dconf->ss_slices_p[i].size);

                eidx -= dconf->ss_slices_p[i].size;
            }

            //printf("dumped %d sieve hits with logp %d for poly %d on side %d bucket %d, bucket index now %d\n",
            //    dconf->ss_slices_p[i].curr_poly_num, logp, dconf->numB, side, i, curr_poly_idx);
        }
    }
    else
    {
        // find the first hit for this poly.  should probably instead store
        // a list of first locations for each poly.
        uint32_t eidx = 0;
        //while ((dconf->ss_slices_n[0].element[eidx] & 0xffff) != dconf->numB)
        //{
        //    eidx++;
        //}
        eidx = dconf->poly_ll_first_n[dconf->numB];
        //printf("first element index for poly %u side %d = %u\n", dconf->numB, side, eidx);
        uint32_t root;

        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint64_t* elements = dconf->ss_slices_n[i].element;
            uint8_t logp = dconf->ss_slices_n[i].logp;
            uint32_t nextentry;
            uint32_t slicesize = dconf->ss_slices_n[i].size;

            do
            {
                //if ((dconf->ss_slices_n[i].element[eidx] & 0xffff) != dconf->numB)
                //{
                //    printf("error: next element index %u is not the same poly index (now %u)!\n",
                //        eidx, (dconf->ss_slices_n[i].element[eidx] & 0xffff));
                //    printf("faulty entry = %lx\n", dconf->ss_slices_n[i].element[eidx - nextentry]);
                //    printf("current entry = %lx\n", dconf->ss_slices_n[i].element[eidx]);
                //    printf("correctly processed %u entries for poly %d\n",
                //        numentries, dconf->numB);
                //    exit(1);
                //}
                nextentry = elements[eidx] >> 40;
                _mm_prefetch(&elements[eidx + nextentry], _MM_HINT_T1);

                root = (elements[eidx] & 0xffffff);
                sieve[root] += logp;
                eidx += nextentry;
            } while ((eidx < slicesize) && (nextentry != 0xffffff));

            if (nextentry == 0xffffff)
            {
                // last entry for this poly
                break;
            }
            else
            {
                // more entries in the next slice.  adjust the element index
                // for the remainder of this slice's size.
                //printf("finished processing %d entries from slice %d, adjusting "
                //    "next entry index %d by size of slice %d (%d) to %d\n",
                //    numentries, i, eidx, i, dconf->ss_slices_n[i].size,
                //    eidx - dconf->ss_slices_n[i].size);

                eidx -= dconf->ss_slices_n[i].size;
            }

            //printf("dumped %d sieve hits with logp %d for poly %d on side %d bucket %d, bucket index now %d\n",
            //    dconf->ss_slices_p[i].curr_poly_num, logp, dconf->numB, side, i, curr_poly_idx);
        }
    }

    return;
}

void lp_sieve_ss_linked_links(uint8_t* sieve, int side, dynamic_conf_t* dconf)
{
    // with linked links by poly
    int i;
    int j;
    //uint32_t numentries = 0;

    if (side == 0)
    {
        // find the first hit for this poly.  should probably instead store
        // a list of first locations for each poly.
        uint32_t eidx = 0;
        //while ((dconf->ss_slices_p[0].element[eidx] & 0xffff) != dconf->numB)
        //{
        //    eidx++;
        //}
        eidx = dconf->poly_ll_first_p[dconf->numB];

        //printf("first element index for poly %u side %d = %u\n", dconf->numB, side, eidx);
        uint32_t root;

        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint64_t* elements = dconf->ss_slices_p[i].element;           
            uint8_t logp = dconf->ss_slices_p[i].logp;
            uint32_t nextentry;
            uint32_t slicesize = dconf->ss_slices_p[i].size;
            //__m512i vbuckets, vnextbuckets, vlosieve, vhisieve;
            //__m512i vlogp = _mm512_set1_epi32(logp);
            //ALIGNED_MEM uint32_t roots[16];
            //int j;

            do 
            {
                //if ((dconf->ss_slices_p[i].element[eidx] & 0xffff) != dconf->numB)
                //{
                //    printf("error: next element index %u is not the same poly index (now %u)!\n",
                //        eidx, (dconf->ss_slices_p[i].element[eidx] & 0xffff));
                //    printf("faulty entry = %lx\n", dconf->ss_slices_p[i].element[eidx - nextentry]);
                //    printf("current entry = %lx\n", dconf->ss_slices_p[i].element[eidx]);
                //    printf("correctly processed %u entries for poly %d\n",
                //        numentries, dconf->numB);
                //    exit(1);
                //}
                
                //__mmask16 m = 0;
                //for (j = 0; j < 16; j++)
                //{
                //    nextentry = elements[eidx] >> 48;
                //    _mm_prefetch(&elements[eidx + nextentry], _MM_HINT_T1);
                //    roots[j] = (elements[eidx] & 0xffff0000) >> 16;
                //    m |= (1 << j);
                //    eidx += nextentry;
                //    if ((eidx >= slicesize) || (nextentry == 0xffff))
                //    {
                //        break;
                //    }
                //}
                //
                //vbuckets = _mm512_load_epi32(roots);
                //vhisieve = _mm512_mask_i32gather_epi32(vbuckets, m, vbuckets, sieve, _MM_SCALE_1);
                //vlosieve = _mm512_sub_epi8(vhisieve, vlogp);
                //_mm512_mask_i32scatter_epi32(sieve, m, vbuckets, vlosieve, _MM_SCALE_1);

                nextentry = elements[eidx] >> 40;
                _mm_prefetch(&elements[eidx + nextentry], _MM_HINT_T1);

                root = (elements[eidx] & 0xffffff);
                sieve[root] -= logp;
                eidx += nextentry;
                
                //numentries++;
            } while ((eidx < slicesize) && (nextentry != 0xffffff));


            

            if (nextentry == 0xffffff)
            {
                // last entry for this poly
                break;
            }
            else
            {
                // more entries in the next slice.  adjust the element index
                // for the remainder of this slice's size.
                //printf("finished processing %d entries from slice %d, adjusting "
                //    "next entry index %d by size of slice %d (%d) to %d\n",
                //    numentries, i, eidx, i, dconf->ss_slices_p[i].size,
                //    eidx - dconf->ss_slices_p[i].size);

                eidx -= dconf->ss_slices_p[i].size;
            }

            //printf("dumped %d sieve hits with logp %d for poly %d on side %d bucket %d, bucket index now %d\n",
            //    dconf->ss_slices_p[i].curr_poly_num, logp, dconf->numB, side, i, curr_poly_idx);
        }
    }
    else
    {
        // find the first hit for this poly.  should probably instead store
        // a list of first locations for each poly.
        uint32_t eidx = 0;
        //while ((dconf->ss_slices_n[0].element[eidx] & 0xffff) != dconf->numB)
        //{
        //    eidx++;
        //}
        eidx = dconf->poly_ll_first_n[dconf->numB];
        //printf("first element index for poly %u side %d = %u\n", dconf->numB, side, eidx);
        uint32_t root;

        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint64_t* elements = dconf->ss_slices_n[i].element;
            uint8_t logp = dconf->ss_slices_n[i].logp;
            uint32_t nextentry;
            uint32_t slicesize = dconf->ss_slices_n[i].size;

            do
            {
                //if ((dconf->ss_slices_n[i].element[eidx] & 0xffff) != dconf->numB)
                //{
                //    printf("error: next element index %u is not the same poly index (now %u)!\n",
                //        eidx, (dconf->ss_slices_n[i].element[eidx] & 0xffff));
                //    printf("faulty entry = %lx\n", dconf->ss_slices_n[i].element[eidx - nextentry]);
                //    printf("current entry = %lx\n", dconf->ss_slices_n[i].element[eidx]);
                //    printf("correctly processed %u entries for poly %d\n",
                //        numentries, dconf->numB);
                //    exit(1);
                //}
                nextentry = elements[eidx] >> 40;
                _mm_prefetch(&elements[eidx + nextentry], _MM_HINT_T1);

                root = (elements[eidx] & 0xffffff);
                sieve[root] -= logp;
                eidx += nextentry;
            } while ((eidx < slicesize) && (nextentry != 0xffffff));

            if (nextentry == 0xffffff)
            {
                // last entry for this poly
                break;
            }
            else
            {
                // more entries in the next slice.  adjust the element index
                // for the remainder of this slice's size.
                //printf("finished processing %d entries from slice %d, adjusting "
                //    "next entry index %d by size of slice %d (%d) to %d\n",
                //    numentries, i, eidx, i, dconf->ss_slices_n[i].size,
                //    eidx - dconf->ss_slices_n[i].size);

                eidx -= dconf->ss_slices_n[i].size;
            }

            //printf("dumped %d sieve hits with logp %d for poly %d on side %d bucket %d, bucket index now %d\n",
            //    dconf->ss_slices_p[i].curr_poly_num, logp, dconf->numB, side, i, curr_poly_idx);
        }
    }

    return;
}
#endif

#ifdef USE_SORTED_LIST_SS
void lp_sieve_ss(uint8_t* sieve, int side, dynamic_conf_t* dconf)
{
    // _sorted_buckets
    int i;
    int j;

    if (side == 0)
    {    
        //printf("dumping %d sieve hits on side %d for poly %d\n",
        //    dconf->ss_size_p[dconf->numB], side, dconf->numB);
        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint64_t* elements = dconf->ss_slices_p[i].element;
            uint32_t curr_poly_idx = dconf->ss_slices_p[i].curr_poly_idx;
            uint32_t root;
            uint8_t logp = dconf->ss_slices_p[i].logp;

            while (((elements[curr_poly_idx] >> 40) == dconf->numB) &&
                (curr_poly_idx < dconf->ss_slices_p[i].size))
            {
                root = (elements[curr_poly_idx] & 0xffff);
                sieve[root] -= logp;
                curr_poly_idx++;
            }
            
            dconf->ss_slices_p[i].curr_poly_num = curr_poly_idx - dconf->ss_slices_p[i].curr_poly_idx;
            dconf->ss_slices_p[i].curr_poly_idx = curr_poly_idx;
            //printf("dumped %d sieve hits with logp %d for poly %d on side %d bucket %d, bucket index now %d\n",
            //    dconf->ss_slices_p[i].curr_poly_num, logp, dconf->numB, side, i, curr_poly_idx);
        }
    }
    else
    {
        //printf("dumping %d sieve hits on side %d for poly %d\n",
        //    dconf->ss_size_n[dconf->numB], side, dconf->numB);
        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint64_t* elements = dconf->ss_slices_n[i].element;
            uint32_t curr_poly_idx = dconf->ss_slices_n[i].curr_poly_idx;
            uint32_t root;
            uint8_t logp = dconf->ss_slices_n[i].logp;

            while (((elements[curr_poly_idx] >> 40) == dconf->numB) &&
                (curr_poly_idx < dconf->ss_slices_p[i].size))
            {
                root = (elements[curr_poly_idx] & 0xffff);
                sieve[root] -= logp;
                curr_poly_idx++;
            }

            dconf->ss_slices_n[i].curr_poly_num = curr_poly_idx - dconf->ss_slices_n[i].curr_poly_idx;
            dconf->ss_slices_n[i].curr_poly_idx = curr_poly_idx;
            //printf("dumped %d sieve hits with logp %d for poly %d on side %d bucket %d, bucket index now %d\n",
            //    dconf->ss_slices_n[i].curr_poly_num, logp, dconf->numB, side, i, curr_poly_idx);
        }
    }

    return;
}
#endif

#ifdef USE_POLY_BUCKET_SS

#define TRY_PN_BUCKET_COMBINE_SIEVE

void lp_sieve_ss(uint8_t* sieve, int side, dynamic_conf_t* dconf)
{
    int i;
    int j;

    //int pidx = dconf->numB; // dconf->polymap[dconf->numB];
    int pidx = dconf->polymap[dconf->numB];
    int bucketalloc = dconf->ss_slices_p[0].alloc;

    // if the mapped binary-encoded poly isn't in this block of 
    // poly buckets then just skip large prime sieving.
    //if ((pidx < dconf->ss_slices_p[0].curr_poly_idx) ||
    //    (pidx >= (dconf->ss_slices_p[0].curr_poly_idx + (1 << dconf->ss_set1.size))))
    //    return;

    //printf("commencing large_sieve_ss on side %d, pidx %u (gray-index %d, set2 instance %d), sizes %d,%d\n",
    //    side, pidx, dconf->numB, polymask / (1 << dconf->ss_set2.size),
    //    (1 << dconf->ss_set1.size), (1 << dconf->ss_set2.size));

    //printf("mapping b-index %d to bucket %d\n", dconf->numB, pidx);
    __m512i vrootmask = _mm512_set1_epi32(0x1ffff);
    __m512i vposmask = _mm512_set1_epi32(1 << 17);

    if ((side & 1) == 0)
    {
        for (i = 0; i < dconf->num_ss_slices; i++)
        {
            uint32_t* bucketelements = dconf->ss_slices_p[i].elements + pidx * bucketalloc;
            uint32_t root;
            uint8_t logp = dconf->ss_slices_p[i].logp;

#ifdef TRY_PN_BUCKET_COMBINE_SIEVE

            for (j = 0; j < dconf->ss_slices_p[i].size[pidx]; j += 16)
            {
                __m512i vr = _mm512_loadu_epi32(bucketelements + j);
                __mmask16 mpos = ~_mm512_test_epi32_mask(vr, vposmask);
             
                vr = _mm512_and_epi32(vr, vrootmask);
                ALIGNED_MEM uint32_t roots[16];
            
                //_mm512_mask_compressstoreu_epi32(roots, mpos, vr);
                //int num = _mm_popcnt_u32(mpos) - 1;
                //while (num >= 0)
                //{
                //    sieve[roots[num--]] -= logp;
                //}
            
                _mm512_store_epi32(roots, vr);
                
                while (mpos > 0)
                {
                    int idx = _trail_zcnt(mpos);
                    sieve[roots[idx]] -= logp;
                    mpos = _reset_lsb(mpos);
                
                }
            }

            for ( ; j < dconf->ss_slices_p[i].size[pidx]; j++)
            {
                if ((bucketelements[j] & 0x20000) == 0)
                {
                    sieve[bucketelements[j] & 0x1ffff] -= logp;
                }
            }

#else
            for (j = 0; j < dconf->ss_slices_p[i].size[pidx]; j++)
            {
                root = (bucketelements[j] & 0x3ffff);
                sieve[root] -= logp;
            }
#endif
        }

    }
    else
    {
        for (i = 0; i < dconf->num_ss_slices; i++)
        {
#ifdef TRY_PN_BUCKET_COMBINE_SIEVE
            uint32_t* bucketelements = dconf->ss_slices_p[i].elements + pidx * bucketalloc;
            uint32_t root;
            uint8_t logp = dconf->ss_slices_p[i].logp;

            for (j = 0; j < dconf->ss_slices_p[i].size[pidx]; j += 16)
            {
                __m512i vr = _mm512_loadu_epi32(bucketelements + j);
                __mmask16 mpos = _mm512_test_epi32_mask(vr, vposmask);
                vr = _mm512_and_epi32(vr, vrootmask);
                ALIGNED_MEM uint32_t roots[16];
                
                //_mm512_mask_compressstoreu_epi32(roots, mpos, vr);
                //int num = _mm_popcnt_u32(mpos) - 1;
                //while (num >= 0)
                //{
                //    sieve[roots[num--]] -= logp;
                //}
            
                _mm512_store_epi32(roots, vr);
                
                while (mpos > 0)
                {
                    int idx = _trail_zcnt(mpos);
                    sieve[roots[idx]] -= logp;
                    mpos = _reset_lsb(mpos);
                }
            }

            for ( ; j < dconf->ss_slices_p[i].size[pidx]; j++)
            {
                if ((bucketelements[j] & 0x20000))
                {
                    sieve[bucketelements[j] & 0x1ffff] -= logp;
                }
            }
           
#else
            uint32_t* bucketelements = dconf->ss_slices_n[i].elements + pidx * 16384;
            uint32_t root;
            uint8_t logp = dconf->ss_slices_n[i].logp + 1;

            for (j = 0; j < dconf->ss_slices_n[i].size[pidx]; j++)
            {
                root = (bucketelements[j] & 0x3ffff);
                sieve[root] -= logp;
            }
#endif
        }

    }

    return;
}
#endif