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
#include "ytools.h"
#include "common.h"
#include "poly_macros_32k.h"
#include "poly_macros_common.h"
#include <stdlib.h>

#ifdef USE_AVX512F

#include <immintrin.h>

// prefetch instructions available for KNL    
#ifdef USE_AVX512PF
//#define KNL_SCATTER_PREFETCH _MM_HINT_T1
//#define KNL_SCATTER_PREFETCH_NUM _MM_HINT_T1
#endif

// when not enough primes in a set of 16 hit an interval
// it is better to just loop through the ones that hit
// using _trail_zcnt and update buckets one at a time.
//#define KNL_XLARGE_METHOD_1_CUTOFF 7
//#define XLARGE_BUCKET_DEBUG

#define SCATTER_COMPRESSED_VECTOR_P \
for (k = 0; k < idx; k++) { \
    bptr = sliceptr_p + (b1[k] << BUCKET_BITS) + numptr_p[b1[k]]; \
    *bptr = e1[k]; \
    numptr_p[b1[k]]++; \
}

#define SCATTER_COMPRESSED_VECTOR_N \
for (k = 0; k < idx; k++) { \
    bptr = sliceptr_n + (b1[k] << BUCKET_BITS) + numptr_n[b1[k]]; \
    *bptr = e1[k]; \
    numptr_n[b1[k]]++; \
}

int qcomp_xlbucket(const void* x, const void* y)
{
    uint32_t* xx = (uint32_t *)x;
    uint32_t* yy = (uint32_t *)y;
    uint32_t a = (xx[0] >> 15) & 0x7;
    uint32_t b = (yy[0] >> 15) & 0x7;

    if (a > b)
        return 1;
    else if (a == b)
        return 0;
    else
        return -1;
}

void nextRoots_32k_knl_bucket(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //update the roots 
    int *rootupdates = dconf->rootupdates;

    update_t update_data = dconf->update_data;

    uint32_t startprime = 2;
    uint32_t bound = sconf->factor_base->B;

    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];

#ifdef USE_BATCHPOLY_X2
    char* nu = dconf->curr_poly->nu;
    char* gray = dconf->curr_poly->gray;
    int numB = dconf->numB;
    uint32_t poly_offset = 2 * sconf->num_blocks * dconf->buckets->alloc_slices;
    int p, pnum;
    uint32_t gmask[64];
#endif

#ifdef USE_XLBUCKET
    xlp_bucket_t *xlp_bucket_p = &dconf->xl_pbucket;
    xlp_bucket_t* xlp_bucket_n = &dconf->xl_nbucket;
#endif

    int* ptr;
    lp_bucket *lp_bucket_p = dconf->buckets;
    uint32_t med_B = sconf->factor_base->med_B;
    uint32_t large_B = sconf->factor_base->large_B;
    uint32_t xlarge_B = sconf->factor_base->x2_large_B;

    uint32_t j, interval;
    int k,numblocks,idx;
    //uint32_t root1, prime;

    int bound_index=0;
    int check_bound = BUCKET_ALLOC/2 - 1;
    uint32_t bound_val = med_B;
    uint32_t *numptr_p, *numptr_n, *sliceptr_p,*sliceptr_n;

    uint8_t* slicelogp_ptr = NULL;
    uint32_t* slicebound_ptr = NULL;
    uint32_t *bptr;
    int room;
    // int bnum, idx_n;

    uint8_t logp=0;

    __m512i vshifted_index = _mm512_setr_epi32(
        0, 65536, 131072, 196608, 
        262144,327680, 393216, 458752, 
        524288,589824, 655360, 720896, 
        786432,851968, 917504, 983040);

    __m512i vshifted21_index = _mm512_setr_epi32(
        0  << 21, 1  << 21, 2  << 21, 3  << 21,
        4  << 21, 5  << 21, 6  << 21, 7  << 21,
        8  << 21, 9  << 21, 10 << 21, 11 << 21,
        12 << 21, 13 << 21, 14 << 21, 15 << 21);

    __m512i vone = _mm512_set1_epi32(1);

    __m512i vnroot1, vnroot2;
    //__m512i vone = _mm512_set1_epi32(1), vzero = _mm512_setzero_epi32();
    __m512i vprime, vroot1, vroot2, vpval, vbnum1, vbnum2, vinterval;
    __m512i velement1, velement2, vblockm1 = _mm512_set1_epi32(32767);
    __mmask16 mask1, mask2;

    // for scatter/gather:
    //__mmask16 mNotProc, mconflictfree, munique;;
    // __m512i vaddr, vidx, vnum, vbsize;
    
    ALIGNED_MEM uint32_t e1[32]; //e1[65536];
    ALIGNED_MEM uint32_t b1[32]; //b1[65536];
    ALIGNED_MEM uint32_t b2[32]; //b1[65536];

    numblocks = sconf->num_blocks;
    interval = numblocks << 15;
    vinterval = _mm512_set1_epi32(interval);

    if (lp_bucket_p->alloc_slices != 0) // != NULL)
    {

#ifdef USE_BATCHPOLY_X2
        // every N iterations we do the bucket sieve for extra large primes.
        // So we only reset the buckets every N iterations.  Otherwise, pick
        // up where we left off.
        if ((dconf->numB % dconf->poly_batchsize) == 1)
        {
            //printf("reset all %u bucket counts for batch\n", dconf->buckets->list_size);

            // reset bucket counts for all polys
            numptr_p = lp_bucket_p->num;
            for (j = 0; j < lp_bucket_p->list_size; j++)
            {
                numptr_p[j] = 0;
            }

            //printf("reset slice counts for batch\n");

            // and the number of slices used per poly.
            for (j = 0; j < dconf->poly_batchsize; j++)
            {
                lp_bucket_p->num_slices_batch[j] = 0;
            }

        }

        // the id of the zero-indexed polynomial in this batch
        pnum = (dconf->numB % dconf->poly_batchsize) - 1;
        if (pnum < 0)
            pnum += dconf->poly_batchsize;

        sliceptr_p = lp_bucket_p->list + pnum * poly_offset * BUCKET_ALLOC;
        sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS) + pnum * poly_offset * BUCKET_ALLOC;

        numptr_p = lp_bucket_p->num + pnum * poly_offset;
        numptr_n = lp_bucket_p->num + numblocks + pnum * poly_offset;

        // if this isn't the first poly in a batch, then we resume
        // where the previous xlarge bucket sieve stopped.
        bound_index = lp_bucket_p->num_slices_batch[pnum];
        slicelogp_ptr = lp_bucket_p->logp + pnum * lp_bucket_p->alloc_slices;
        slicebound_ptr = lp_bucket_p->fb_bounds + pnum * lp_bucket_p->alloc_slices;

        sliceptr_p += bound_index * (numblocks << (BUCKET_BITS + 1));
        sliceptr_n += bound_index * (numblocks << (BUCKET_BITS + 1));
        numptr_p += bound_index * (numblocks << 1);
        numptr_n += bound_index * (numblocks << 1);

        slicebound_ptr[bound_index] = med_B;

        //printf("begin bucket sieve for poly %d with bi = %u, bv = %u, cb = %u\n",
        //    pnum, bound_index, bound_val, med_B + BUCKET_ALLOC / 2);

#else

        lp_bucket_p->fb_bounds[0] = med_B;

        sliceptr_p = lp_bucket_p->list;
        sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

        numptr_p = lp_bucket_p->num;
        numptr_n = lp_bucket_p->num + numblocks;

        slicelogp_ptr = lp_bucket_p->logp;
        slicebound_ptr = lp_bucket_p->fb_bounds;

        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
        {
            numptr_p[j] = 0;
        }
        lp_bucket_p->num_slices = 0;

#endif

    }
    else
    {
        sliceptr_p = NULL;
        sliceptr_n = NULL;
        numptr_p = NULL;
        numptr_n = NULL;
    }

    k=0;
    ptr = &rootupdates[(v-1) * bound + med_B];	

#ifdef USE_BATCHPOLY_X2
    for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
    {
        if (gray[numB + p] > 0)
        {
            gmask[p] = 0xffff;
        }
        else
        {
            gmask[p] = 0;
        }
    }
#endif

    if (sign > 0)
    {        

        //bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC/2;

        logp = update_data.logp[med_B];
        for (j=med_B;j<large_B;j+=16,ptr+=16)				
        {			
            CHECK_NEW_SLICE(j);

            vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i *)ptr);
            mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
            mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
            vroot1 = _mm512_sub_epi32(vroot1, vpval);
            vroot2 = _mm512_sub_epi32(vroot2, vpval);
            vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i *)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i *)(&update_data.firstroots2[j]), vroot2);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);

            // we loop over roots below, so have to compute these before we
            // clobber them.
            vnroot1 = _mm512_sub_epi32(vprime, vroot1);
            vnroot2 = _mm512_sub_epi32(vprime, vroot2);

#ifdef _MSC_VER
            // msvc currently has a bug in the new _mm512_mask_compressstoreu_epi32
            // intrinsic, so we fall back on this almost-as-fast code.

            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m = _blsr_u32(m);
                }

                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            }

            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            }

            mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            }

            mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot2, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            }

#else




            // idea: don't store the prime indices, just the block locations.
            // this will make bucket sorting and lp_sieve way easier, but lp_tdiv
            // harder.  Could do lp_tdiv with a resieve maybe?
            // I think this actually has promise, but is quite a bit of coding.
            // lp_tdiv is trivial right now while bucket sorting is by far the
            // most time consuming process.  Anything we can do to offset that 
            // could be a nice win.
            // looking at what lp_tdiv would involve... update_data.firstroots
            // holds the starting position of the prime in this sieve interval.
            // so we could iterate that by prime and check for equality with
            // the candidate block position, like in the resieve code.  In the
            // case of very large primes, would only need to check, no iteration
            // at all.  But there are many many such factor base primes, so the
            // scanning process will take much longer.  
            // larger numbers typically have < 1 report per block, so at least
            // the scanning won't happen too often.
            // maybe could even keep the current lp_bucket structure for primes
            // less than the interval size, and only switch to the new structure
            // for primes larger than the interval size.  make a vlp_bucket type
            // and vlp_tdiv.
            // ok, another thought.  for vlp's, don't even do a block sieve, do
            // an interval sieve.  Then we don't have to scatter into appropriate
            // blocks, just compress-write the roots of anything thats hits the
            // interval.  During the block sieve, walk through the list and 
            // increment logp of anything in the current block.  That will be
            // a little slower, but bucket sieving should be massively faster.
            // could also keep the prime index info so we can keep mostly the
            // same lp_tdiv.

            __m512i vebase = velement1;
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot1, vblockm1));
            velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot2, vblockm1));

            do
            {
                //SCATTER_COMPRESSED_VECTOR_P;

                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_p + (k << BUCKET_BITS) + numptr_p[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_p[k] += (idx + _mm_popcnt_u32(mask2));
                }

                vroot1 = _mm512_add_epi32(vroot1, vprime);
                vroot2 = _mm512_add_epi32(vroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

                if ((mask1 | mask2) == 0)
                {
                    break;
                }

                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot1, vblockm1));
                velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot2, vblockm1));

                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot1, vblockm1));
                velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot2, vblockm1));

            } while (1); // (idx > 0);

            //mask1 = 0xffff;
            //mask2 = 0xffff;
            //idx = 32;
            //
            //_mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
            //_mm512_store_epi32((__m512i*)e1,
            //    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));
            //
            //_mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vnroot2, 15));
            //_mm512_store_epi32((__m512i*)(&(e1[16])),
            //    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

            vbnum1 = _mm512_srli_epi32(vnroot1, 15);
            vbnum2 = _mm512_srli_epi32(vnroot2, 15);
            velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot1, vblockm1));
            velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot2, vblockm1));

            do
            {
                //SCATTER_COMPRESSED_VECTOR_N;
                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_n + (k << BUCKET_BITS) + numptr_n[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_n[k] += (idx + _mm_popcnt_u32(mask2));
                }

                vnroot1 = _mm512_add_epi32(vnroot1, vprime);
                vnroot2 = _mm512_add_epi32(vnroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT); 
                
                if ((mask1 | mask2) == 0)
                {
                    break;
                }

                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot1, vblockm1));
                velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot2, vblockm1));

            } while (1);
                        
#endif

        }

        logp = update_data.logp[j - 1];

#if defined( USE_BATCHPOLY_X2 ) || defined(USE_XLBUCKET)
        for (j = large_B; j < xlarge_B; j += 16, ptr += 16)
#else
        for (j = large_B; j < xlarge_B; j += 16, ptr += 16)
#endif
        {
            //int i;
            int k;

            CHECK_NEW_SLICE(j);

            vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i*)ptr);
            mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
            mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
            vroot1 = _mm512_sub_epi32(vroot1, vpval);
            vroot2 = _mm512_sub_epi32(vroot2, vpval);
            vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

            dconf->num_scatter_opp++;

            idx = _mm_popcnt_u32(mask1);
            if ((idx >= KNL_XLARGE_METHOD_1_CUTOFF))
            {
                //vbnum1 = _mm512_load_epi32((__m512i*)b1);
                //velement1 = _mm512_load_epi32((__m512i*)e1);
                //mask1 = (1 << idx) - 1;

                dconf->num_scatter++;
                // get the conflicted indices.  Note that the mask is a *writemask*... all
                // lanes will participate in the read and conflict instruction.  
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);

                // identify the location of the last index of each unique block.  This location
                // will eventually hold the largest increment for each unique block.
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);

                // gather the num_p values for the blocks that are hit
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

                // try scatter prefetching?
#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                // start out having "processed" those indices that don't hit the interval.
                mNotProc = mask1;

                while (mNotProc != 0)
                {

                    // set a mask of non-conflicted indices that are also interval hits
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    // write the bucket elements.
                     // instead of scatter in the loop, compress-write to a register and do scatter once at the end?
                     // yes, this is better.  in fact, we just need to update vnum in a conflict-free manner.
                     //_mm512_mask_i32scatter_epi32(sliceptr_p, mconflictfree & ~mprocessed, vaddr, velement1, _MM_SCALE_4);                           

                     // adjust the conflict masks (clear the conflicted bits that we've just identified)
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);

                    // remove the lanes we wrote this time
                    mNotProc &= (~mconflictfree);

                    // increment the counters that are still conflicted
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }

                // compute the bucket indices
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, velement1, _MM_SCALE_4);

                // update the num_p values (increment)
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
            }
            else
            {
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);

                SCATTER_COMPRESSED_VECTOR_P;
            }

            idx = _mm_popcnt_u32(mask2);
            if ((idx >= KNL_XLARGE_METHOD_1_CUTOFF))
            {
                vbnum1 = vbnum2;
                velement1 = velement2;
                mask1 = mask2;

                //vbnum1 = _mm512_load_epi32((__m512i*)b1);
                //velement1 = _mm512_load_epi32((__m512i*)e1);
                //mask1 = (1 << idx) - 1;

                dconf->num_scatter++;
                // get the conflicted indices.  Note that the mask is a *writemask*... all
                // lanes will participate in the read and conflict instruction.  
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);

                // identify the location of the last index of each unique block.  This location
                // will eventually hold the largest increment for each unique block.
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);

                // gather the num_p values for the blocks that are hit
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

                // try scatter prefetching?
#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                // start out having "processed" those indices that don't hit the interval.
                mNotProc = mask1;

                while (mNotProc != 0)
                {

                    // set a mask of non-conflicted indices that are also interval hits
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    // write the bucket elements.
                     // instead of scatter in the loop, compress-write to a register and do scatter once at the end?
                     // yes, this is better.  in fact, we just need to update vnum in a conflict-free manner.
                     //_mm512_mask_i32scatter_epi32(sliceptr_p, mconflictfree & ~mprocessed, vaddr, velement1, _MM_SCALE_4);                           

                     // adjust the conflict masks (clear the conflicted bits that we've just identified)
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);

                    // remove the lanes we wrote this time
                    mNotProc &= (~mconflictfree);

                    // increment the counters that are still conflicted
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }

                // compute the bucket indices
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, velement1, _MM_SCALE_4);

                // update the num_p values (increment)
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
            }
            else
            {
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask2, vbnum2);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask2, velement2);

                SCATTER_COMPRESSED_VECTOR_P;
            }

#elif defined(_MSC_VER)

            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

            // msvc currently has a bug in the new _mm512_mask_compressstoreu_epi32
            // intrinsic, so we fall back on this almost-as-fast code.
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask1 ^= (1 << idx);
            }

#else

            //if (1)
            //{
            //    mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            //    mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            //    idx = _mm_popcnt_u32(mask1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
            //    idx += _mm_popcnt_u32(mask2);
            //
            //    SCATTER_COMPRESSED_VECTOR_P;
            //}
            //else
            {
                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_p + (k << BUCKET_BITS) + numptr_p[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_p[k] += (idx + _mm_popcnt_u32(mask2));
                }
            }
#endif

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            dconf->num_scatter_opp++;
            idx = _mm_popcnt_u32(mask1);
            if ((idx >= KNL_XLARGE_METHOD_1_CUTOFF))
            {
                //vbnum1 = _mm512_load_epi32((__m512i*)b1);
                //velement1 = _mm512_load_epi32((__m512i*)e1);
                //mask1 = (1 << idx) - 1;

            // msvc currently has a bug in the new _mm512_mask_compressstoreu_epi32
            // intrinsic, so we fall back on this almost-as-fast code.
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
            }

#else
            //if (1)
            //{
            //    mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            //    mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            //    idx = _mm_popcnt_u32(mask1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
            //    idx += _mm_popcnt_u32(mask2);
            //
            //    SCATTER_COMPRESSED_VECTOR_N;
            //}
            //else
            {
                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_n + (k << BUCKET_BITS) + numptr_n[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_n[k] += (idx + _mm_popcnt_u32(mask2));
                }
            }
#endif


        }
#endif

        // above this bound primes are much bigger and therefore
        // less likely to hit the interval.  For the ones that do,
        // there will be more blocks involved (because primes only
        // get this big for bigger inputs) so less likely to have
        // multiple hits per block.  For these reasons it is more
        // efficient to use the individual scatter rather than 
        // the block-based compress-store.
        
#ifndef USE_XLBUCKET
        for (j = xlarge_B; j < bound; j += 16, ptr += 16)
        {
            //int i;
            int k;

            CHECK_NEW_SLICE(j);

            vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i*)ptr);
            mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
            mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
            vroot1 = _mm512_sub_epi32(vroot1, vpval);
            vroot2 = _mm512_sub_epi32(vroot2, vpval);
            vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));

#if defined(_MSC_VER)

            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

            // msvc currently has a bug in the new _mm512_mask_compressstoreu_epi32
            // intrinsic, so we fall back on this almost-as-fast code.
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask1 ^= (1 << idx);
            }

#else

            //if (1)
            {
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                idx = _mm_popcnt_u32(mask1);
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                idx += _mm_popcnt_u32(mask2);
            
                SCATTER_COMPRESSED_VECTOR_P;
            }
            
#endif

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));

#if defined(_MSC_VER)

            // msvc currently has a bug in the new _mm512_mask_compressstoreu_epi32
            // intrinsic, so we fall back on this almost-as-fast code.
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
            }

#else
            //if (1)
            {
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                idx = _mm_popcnt_u32(mask1);
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                idx += _mm_popcnt_u32(mask2);
            
                SCATTER_COMPRESSED_VECTOR_N;
            }
            
            
#endif


        }
#endif
    }
    else
    {
        // all of this else case is identical to the if case except we
        // calculate the roots slightly differently (modular addition
        // versus modular subtraction).

        // note: the reason for the top-level if-else is one of the main
        // reasons why this is fast - it gets all of the polynomial
        // switching related branching out of the inner loops of the bucket 
        // sort.
        //bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC / 2;


        logp = update_data.logp[med_B];
        for (j = med_B; j<large_B; j += 16, ptr += 16)
        {
            CHECK_NEW_SLICE(j);

            // load the current roots, associated primes, and polynomial
            // update info (ptr).  Do the polynomial update and then
            // mask off values that are larger than the sieving interval.
            vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i *)ptr);
            vroot1 = _mm512_add_epi32(vroot1, vpval);
            vroot2 = _mm512_add_epi32(vroot2, vpval);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i *)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i *)(&update_data.firstroots2[j]), vroot2);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);

            // we loop over roots below, so have to compute these before we
            // clobber them.
            vnroot1 = _mm512_sub_epi32(vprime, vroot1);
            vnroot2 = _mm512_sub_epi32(vprime, vroot2);

#ifdef _MSC_VER

            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m = _reset_lsb(m);
                }

                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            }

            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            }

            mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            }

            mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot2, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                while (m > 0)
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            }

#else

            //mask1 = 0xffff;
            //mask2 = 0xffff;
            //idx = 32;
            //
            //_mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot1, 15));
            //_mm512_store_epi32((__m512i*)e1,
            //    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));
            //
            //_mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vroot2, 15));
            //_mm512_store_epi32((__m512i*)(&(e1[16])),
            //    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

            __m512i vebase = velement1;
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot1, vblockm1));
            velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot2, vblockm1));

            do
            {
                //SCATTER_COMPRESSED_VECTOR_P;

                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_p + (k << BUCKET_BITS) + numptr_p[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_p[k] += (idx + _mm_popcnt_u32(mask2));
                }

                vroot1 = _mm512_add_epi32(vroot1, vprime);
                vroot2 = _mm512_add_epi32(vroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

                if ((mask1 | mask2) == 0)
                {
                    break;
                }

                //_mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vroot1, 15));
                //_mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                //    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));
                //
                //idx = _mm_popcnt_u32(mask1);
                //_mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vroot2, 15));
                //_mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                //    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));
                //
                //idx += _mm_popcnt_u32(mask2);
                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot1, vblockm1));
                velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot2, vblockm1));

            } while (1); // (idx > 0);

            //mask1 = 0xffff;
            //mask2 = 0xffff;
            //idx = 32;
            //
            //_mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
            //_mm512_store_epi32((__m512i*)e1,
            //    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));
            //
            //_mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vnroot2, 15));
            //_mm512_store_epi32((__m512i*)(&(e1[16])),
            //    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

            vbnum1 = _mm512_srli_epi32(vnroot1, 15);
            vbnum2 = _mm512_srli_epi32(vnroot2, 15);
            velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot1, vblockm1));
            velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot2, vblockm1));

            do
            {
                //SCATTER_COMPRESSED_VECTOR_N;

                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_n + (k << BUCKET_BITS) + numptr_n[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_n[k] += (idx + _mm_popcnt_u32(mask2));
                }

                vnroot1 = _mm512_add_epi32(vnroot1, vprime);
                vnroot2 = _mm512_add_epi32(vnroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

                if ((mask1 | mask2) == 0)
                {
                    break;
                }

                //_mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vnroot1, 15));
                //_mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                //    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));
                //
                //idx = _mm_popcnt_u32(mask1);
                //_mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vnroot2, 15));
                //_mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                //    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));
                //
                //idx += _mm_popcnt_u32(mask2);
                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot1, vblockm1));
                velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot2, vblockm1));

            } while (1); //(idx > 0);

#endif


        }

        logp = update_data.logp[j - 1];


#ifdef DO_VLP_OPT

        // force a new slice
        if ((lp_bucket_p->list != NULL) && (bound > large_B))
        {
            j = large_B;
            logp = update_data.logp[j];
            lp_bucket_p->logp[bound_index] = logp;
            bound_index++;
            // signal to large_sieve that we are now doing VLP_OPT.
            // we are not storing prime indices anyway, so this value 
            // is useless for lp_tdiv.
            lp_bucket_p->fb_bounds[bound_index] = 0xffffffff;
            bound_val = j;
            sliceptr_p += (numblocks << (BUCKET_BITS + 1));
            sliceptr_n += (numblocks << (BUCKET_BITS + 1));
            numptr_p += (numblocks << 1);
            numptr_n += (numblocks << 1);
            check_bound = bound_val + ((BUCKET_ALLOC * numblocks) >> 1);
        }

        for (j = large_B; j < bound; j += 16, ptr += 16)
        {
            int i;
            int k;

            if (j >= check_bound)
            {
                room = BUCKET_ALLOC * numblocks - MAX(numptr_p[0], numptr_n[0]);
                /* if it is filled close to the allocation, start recording in a new set of buckets */
                if (room < 64)
                {
                    //printf("nextroots-: bucket full, added %u,%u elements, starting new slice %d\n",
                    //    numptr_p[0], numptr_n[0], bound_index + 1);
                    logp = update_data.logp[j];
                    lp_bucket_p->logp[bound_index] = logp;
                    bound_index++;
                    lp_bucket_p->fb_bounds[bound_index] = j;
                    bound_val = j;
                    sliceptr_p += (numblocks << (BUCKET_BITS + 1));
                    sliceptr_n += (numblocks << (BUCKET_BITS + 1));
                    numptr_p += (numblocks << 1);
                    numptr_n += (numblocks << 1);
                    check_bound += ((BUCKET_ALLOC * numblocks) >> 1);
                }
                else
                {
                    check_bound += room >> 1;
                }
            }
            else if ((j - bound_val) >= 65536)
            {
                //printf("nextroots-: prime slice limit, added %u,%u elements, starting new slice %d\n",
                //    numptr_p[0], numptr_n[0], bound_index + 1);
                lp_bucket_p->logp[bound_index] = logp;
                bound_index++;
                lp_bucket_p->fb_bounds[bound_index] = j;
                bound_val = j;
                sliceptr_p += (numblocks << (BUCKET_BITS + 1));
                sliceptr_n += (numblocks << (BUCKET_BITS + 1));
                numptr_p += (numblocks << 1);
                numptr_n += (numblocks << 1);
                check_bound += ((BUCKET_ALLOC * numblocks) >> 1);
            }

            vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i*)ptr);
            vroot1 = _mm512_add_epi32(vroot1, vpval);
            vroot2 = _mm512_add_epi32(vroot2, vpval);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

            idx = numptr_p[0];
            _mm512_mask_compressstoreu_epi32((__m512i*)(&(sliceptr_p[idx])), mask1, vroot1);
            idx += _mm_popcnt_u32(mask1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(&(sliceptr_p[idx])), mask2, vroot2);
            idx += _mm_popcnt_u32(mask2);
            numptr_p[0] = idx;

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

            idx = numptr_n[0];
            _mm512_mask_compressstoreu_epi32((__m512i*)(&(sliceptr_n[idx])), mask1, vroot1);
            idx += _mm_popcnt_u32(mask1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(&(sliceptr_n[idx])), mask2, vroot2);
            idx += _mm_popcnt_u32(mask2);
            numptr_n[0] = idx;
        }

#else

#if defined( USE_BATCHPOLY_X2 ) || defined(USE_XLBUCKET)
        for (j = large_B; j < xlarge_B; j += 16, ptr += 16)
#else
        for (j = large_B; j < xlarge_B; j += 16, ptr += 16)
#endif
        {
            //int i;

            CHECK_NEW_SLICE(j);

            vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i *)ptr);
            vroot1 = _mm512_add_epi32(vroot1, vpval);
            vroot2 = _mm512_add_epi32(vroot2, vpval);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i *)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i *)(&update_data.firstroots2[j]), vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            dconf->num_scatter_opp++;

            idx = _mm_popcnt_u32(mask1);
            if ((idx >= KNL_XLARGE_METHOD_1_CUTOFF))
            {
                //vbnum1 = _mm512_load_epi32((__m512i*)b1);
                //velement1 = _mm512_load_epi32((__m512i*)e1);
                //mask1 = (1 << idx) - 1;

                dconf->num_scatter++;
                // get the conflicted indices.  Note that the mask is a *writemask*... all
                // lanes will participate in the read and conflict instruction.  
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);

                // identify the location of the last index of each unique block.  This location
                // will eventually hold the largest increment for each unique block.
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);

                // gather the num_p values for the blocks that are hit
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

                // try scatter prefetching?
#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                // start out having "processed" those indices that don't hit the interval.
                mNotProc = mask1;

                while (mNotProc != 0)
                {

                    // set a mask of non-conflicted indices that are also interval hits
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    // write the bucket elements.
                     // instead of scatter in the loop, compress-write to a register and do scatter once at the end?
                     // yes, this is better.  in fact, we just need to update vnum in a conflict-free manner.
                     //_mm512_mask_i32scatter_epi32(sliceptr_p, mconflictfree & ~mprocessed, vaddr, velement1, _MM_SCALE_4);                           

                     // adjust the conflict masks (clear the conflicted bits that we've just identified)
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);

                    // remove the lanes we wrote this time
                    mNotProc &= (~mconflictfree);

                    // increment the counters that are still conflicted
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }

                // compute the bucket indices
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, velement1, _MM_SCALE_4);

                // update the num_p values (increment)
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
                }
            else
            {
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);

                SCATTER_COMPRESSED_VECTOR_P;
            }

            idx = _mm_popcnt_u32(mask2);
            if ((idx >= KNL_XLARGE_METHOD_1_CUTOFF))
            {
                vbnum1 = vbnum2;
                velement1 = velement2;
                mask1 = mask2;

                //vbnum1 = _mm512_load_epi32((__m512i*)b1);
                //velement1 = _mm512_load_epi32((__m512i*)e1);
                //mask1 = (1 << idx) - 1;

                dconf->num_scatter++;
                // get the conflicted indices.  Note that the mask is a *writemask*... all
                // lanes will participate in the read and conflict instruction.  
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);

                // identify the location of the last index of each unique block.  This location
                // will eventually hold the largest increment for each unique block.
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);

                // gather the num_p values for the blocks that are hit
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

                // try scatter prefetching?
#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                // start out having "processed" those indices that don't hit the interval.
                mNotProc = mask1;

                while (mNotProc != 0)
                {

                    // set a mask of non-conflicted indices that are also interval hits
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    // write the bucket elements.
                     // instead of scatter in the loop, compress-write to a register and do scatter once at the end?
                     // yes, this is better.  in fact, we just need to update vnum in a conflict-free manner.
                     //_mm512_mask_i32scatter_epi32(sliceptr_p, mconflictfree & ~mprocessed, vaddr, velement1, _MM_SCALE_4);                           

                     // adjust the conflict masks (clear the conflicted bits that we've just identified)
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);

                    // remove the lanes we wrote this time
                    mNotProc &= (~mconflictfree);

                    // increment the counters that are still conflicted
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }

                // compute the bucket indices
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, velement1, _MM_SCALE_4);

                // update the num_p values (increment)
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
            }
            else
            {
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask2, vbnum2);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask2, velement2);

                SCATTER_COMPRESSED_VECTOR_P;
            }
#elif defined(_MSC_VER)

            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask1 ^= (1 << idx);
            }

#else
            //{
            //    mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            //    mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            //    idx = _mm_popcnt_u32(mask1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
            //    idx += _mm_popcnt_u32(mask2);
            //
            //    SCATTER_COMPRESSED_VECTOR_P;
            //}
            //else
            {
                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_p + (k << BUCKET_BITS) + numptr_p[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_p[k] += (idx + _mm_popcnt_u32(mask2));
                }
            }
#endif

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            dconf->num_scatter_opp++;
            idx = _mm_popcnt_u32(mask1);
            if ((idx >= KNL_XLARGE_METHOD_1_CUTOFF))
            {
                //vbnum1 = _mm512_load_epi32((__m512i*)b1);
                //velement1 = _mm512_load_epi32((__m512i*)e1);
                //mask1 = (1 << idx) - 1;

                dconf->num_scatter++;
                // get the conflicted indices.  Note that the mask is a *writemask*... all
                // lanes will participate in the read and conflict instruction.  
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);

                // identify the location of the last index of each unique block.  This location
                // will eventually hold the largest increment for each unique block.
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);

                // gather the num_p values for the blocks that are hit
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

                // try scatter prefetching?
#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                // start out having "processed" those indices that don't hit the interval.
                mNotProc = mask1;

                while (mNotProc != 0)
                {

                    // set a mask of non-conflicted indices that are also interval hits
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    // write the bucket elements.
                     // instead of scatter in the loop, compress-write to a register and do scatter once at the end?
                     // yes, this is better.  in fact, we just need to update vnum in a conflict-free manner.
                     //_mm512_mask_i32scatter_epi32(sliceptr_p, mconflictfree & ~mprocessed, vaddr, velement1, _MM_SCALE_4);                           

                     // adjust the conflict masks (clear the conflicted bits that we've just identified)
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);

                    // remove the lanes we wrote this time
                    mNotProc &= (~mconflictfree);

                    // increment the counters that are still conflicted
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }

                // compute the bucket indices
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_n, mask1, vaddr, velement1, _MM_SCALE_4);

                // update the num_p values (increment)
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_n, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
            }
            else
            {
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);

                SCATTER_COMPRESSED_VECTOR_N;
            }

            idx = _mm_popcnt_u32(mask2);
            if ((idx >= KNL_XLARGE_METHOD_1_CUTOFF))
            {
                vbnum1 = vbnum2;
                velement1 = velement2;
                mask1 = mask2;

                //vbnum1 = _mm512_load_epi32((__m512i*)b1);
                //velement1 = _mm512_load_epi32((__m512i*)e1);
                //mask1 = (1 << idx) - 1;

                dconf->num_scatter++;
                // get the conflicted indices.  Note that the mask is a *writemask*... all
                // lanes will participate in the read and conflict instruction.  
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);

                // identify the location of the last index of each unique block.  This location
                // will eventually hold the largest increment for each unique block.
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);

                // gather the num_p values for the blocks that are hit
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

                // try scatter prefetching?
#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                // start out having "processed" those indices that don't hit the interval.
                mNotProc = mask1;

                while (mNotProc != 0)
                {

                    // set a mask of non-conflicted indices that are also interval hits
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    // write the bucket elements.
                     // instead of scatter in the loop, compress-write to a register and do scatter once at the end?
                     // yes, this is better.  in fact, we just need to update vnum in a conflict-free manner.
                     //_mm512_mask_i32scatter_epi32(sliceptr_p, mconflictfree & ~mprocessed, vaddr, velement1, _MM_SCALE_4);                           

                     // adjust the conflict masks (clear the conflicted bits that we've just identified)
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);

                    // remove the lanes we wrote this time
                    mNotProc &= (~mconflictfree);

                    // increment the counters that are still conflicted
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }

                // compute the bucket indices
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_n, mask1, vaddr, velement1, _MM_SCALE_4);

                // update the num_p values (increment)
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_n, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
            }
            else
            {
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask2, vbnum2);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask2, velement2);

                SCATTER_COMPRESSED_VECTOR_N;
            }
#elif defined(_MSC_VER)
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
            }

#else

            //{
            //    mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            //    mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            //    idx = _mm_popcnt_u32(mask1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
            //    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
            //    idx += _mm_popcnt_u32(mask2);
            //
            //    SCATTER_COMPRESSED_VECTOR_N;
            //}
            //else
            {
                __m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                for (k = 0; k < numblocks; k++)
                {
                    mask1 = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    mask2 = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_n + (k << BUCKET_BITS) + numptr_n[k];
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), mask2, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_n[k] += (idx + _mm_popcnt_u32(mask2));
                }
            }
#endif


        }
#endif

#ifndef USE_XLBUCKET
        for (j = xlarge_B; j < bound; j += 16, ptr += 16)
        {
            //int i;

            CHECK_NEW_SLICE(j);

            vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i*)ptr);
            vroot1 = _mm512_add_epi32(vroot1, vpval);
            vroot2 = _mm512_add_epi32(vroot2, vpval);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));

#if defined(_MSC_VER)

            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask1 ^= (1 << idx);
            }

#else
            {
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                idx = _mm_popcnt_u32(mask1);
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                idx += _mm_popcnt_u32(mask2);
            
                SCATTER_COMPRESSED_VECTOR_P;
            }
            
#endif

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));


#if defined(_MSC_VER)
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            _mm512_store_epi32((__m512i*)b1, vbnum1);
            _mm512_store_epi32((__m512i*)e1, velement1);

            while (mask1 > 0)
            {
                idx = _trail_zcnt(mask1);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
            }

#else

            {
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                idx = _mm_popcnt_u32(mask1);
                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                idx += _mm_popcnt_u32(mask2);
            
                SCATTER_COMPRESSED_VECTOR_N;
            }
            
#endif


        }
#endif
    }

#ifdef USE_BATCHPOLY_X2



    // need to slightly change how buckets operate.  the buckets associated with the 
    // current poly will be filled with primes up to fb->large_B while the buckets
    // associated with the rest of the polys in this batch will be new.  So bound_index,
    // (which tracks the current slice), bound_val (which tracks the prime index that
    // starts each slice), and check_bound (which control when we look at the buckets
    // to potentially start a new slice), should be different for each poly.  


    // every N iterations we do the bucket sieve for extra large primes
    if ((dconf->numB % dconf->poly_batchsize) == 1)
    {
        uint32_t bi[16];    // bound_index, per poly (16 max)
        uint32_t cb[16];    // check_bound, per poly (16 max)
        uint32_t bv[16];    // bound_val, per poly (16 max)

        for (k = 1; k < dconf->poly_batchsize; k++)
        {
            bi[k] = 0;
            cb[k] = j + BUCKET_ALLOC / 2;
            bv[k] = slicebound_ptr[k * lp_bucket_p->alloc_slices] = j;
        }
        
        // check how much room is left for the current poly
        room = 0;
        for (k = 0; k < numblocks; k++)
        {
            if (*(numptr_p + k) > room)
                room = *(numptr_p + k);
            if (*(numptr_n + k) > room)
                room = *(numptr_n + k);
        }

        // fill in the check bound and bound index for the current poly
        cb[pnum] = j + (BUCKET_ALLOC - room) / 2;
        bi[pnum] = bound_index;
        bv[pnum] = bound_val;

        //printf("begin xlarge bucket sieve at j = %u\n", j);
        ////for (k = 0; k < dconf->poly_batchsize; k++)
        //{
        //    int jj;
        //
        //    for (jj = 0; jj < lp_bucket_p->list_size; jj++)
        //    {
        //        int np = jj / (2 * numblocks * lp_bucket_p->alloc_slices);
        //        if ((jj % (2 * numblocks * lp_bucket_p->alloc_slices)) == 0)
        //        {
        //            printf("poly %d: num slices = %u\n", np, bi[np]);
        //        }
        //        if ((jj % (2 * numblocks)) == 0)
        //        {
        //            int ns = (jj - (np * 2 * numblocks * lp_bucket_p->alloc_slices)) / (2 * numblocks);
        //            printf("\tslice %d: fb_bounds = %u, logp = %u\n", ns,
        //                lp_bucket_p->fb_bounds[np * lp_bucket_p->alloc_slices + ns],
        //                lp_bucket_p->logp[np * lp_bucket_p->alloc_slices + ns]);
        //        }
        //        printf("\t\tblock %d: num = %u\n", jj - (np * 2 * numblocks * lp_bucket_p->alloc_slices),
        //            lp_bucket_p->num[jj]);
        //
        //    }
        //}

        slicelogp_ptr = lp_bucket_p->logp; // +pnum * lp_bucket_p->alloc_slices;
        slicebound_ptr = lp_bucket_p->fb_bounds; // +pnum * lp_bucket_p->alloc_slices;

        // bucket sieve over primes and polys
        logp = update_data.logp[j - 1];
        for (; j < bound; j += 16)
        {
            int p;

            vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));

            // todo:
            // compute all roots for N polynomials.
            // compress-store those less than interval to temporary location.
            // adjust vbnum by the poly offset so scatter/gather grabs from
            // blocks across the various poly storage locations.
            // do 1 scatter/gather pass across the temp storage location to
            // distribute roots for all polynomials that hit the interval.

            for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
            {

                __mmask16 gm = gmask[p];

                //slicelogp_ptr = lp_bucket_p->logp + p * lp_bucket_p->alloc_slices;
                //slicebound_ptr = lp_bucket_p->fb_bounds + p * lp_bucket_p->alloc_slices;

                sliceptr_p = lp_bucket_p->list + p * poly_offset * BUCKET_ALLOC;
                sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS) + p * poly_offset * BUCKET_ALLOC;

                numptr_p = lp_bucket_p->num + p * poly_offset;
                numptr_n = lp_bucket_p->num + numblocks + p * poly_offset;

                sliceptr_p += bi[p] * (numblocks << (BUCKET_BITS + 1));
                sliceptr_n += bi[p] * (numblocks << (BUCKET_BITS + 1));
                numptr_p += bi[p] * (numblocks << 1);
                numptr_n += bi[p] * (numblocks << 1);


                CHECK_NEW_SLICE_BATCH_2(j);

                vpval = _mm512_load_epi32((__m512i*)(&rootupdates[(nu[numB + p] - 1) * bound + j]));

                //if (gray[numB + p] > 0)
                mask1 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_mask_sub_epi32(vroot1, gm, vroot1, vpval);
                vroot2 = _mm512_mask_sub_epi32(vroot2, gm, vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, gm & mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, gm & mask2, vroot2, vprime);

                //else
                vroot1 = _mm512_mask_add_epi32(vroot1, ~gm, vroot1, vpval);
                vroot2 = _mm512_mask_add_epi32(vroot2, ~gm, vroot2, vpval);
                mask1 |= _mm512_mask_cmp_epu32_mask(~gm, vroot1, vprime, _MM_CMPINT_GE);
                mask2 |= _mm512_mask_cmp_epu32_mask(~gm, vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, ~gm & mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, ~gm & mask2, vroot2, vprime);


                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bv[p]) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_P;
                }

                // and the -side roots; same story.
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bv[p]) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_N;
                }

                // advance pointers
                sliceptr_p += poly_offset * BUCKET_ALLOC;
                sliceptr_n += poly_offset * BUCKET_ALLOC;
                numptr_p += poly_offset;
                numptr_n += poly_offset;

            }

            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);

            // reset pointers
            sliceptr_p -= p * poly_offset * BUCKET_ALLOC;
            sliceptr_n -= p * poly_offset * BUCKET_ALLOC;
            numptr_p -= p * poly_offset;
            numptr_n -= p * poly_offset;
        }

        // need to track the N bound_index values per poly while sieving and
        // record the final values here.
        if (lp_bucket_p->list != NULL)
        {

            // record the final number of slices and logp value per poly.
            for (k = 0; k < dconf->poly_batchsize; k++)
            {
                lp_bucket_p->num_slices_batch[k] = bi[k] + 1;
                lp_bucket_p->logp[k * lp_bucket_p->alloc_slices + bi[k]] = logp;
            }

            //int jj;
            //printf("end xlarge bucket sieve at j = %u\n", j);
            //for (jj = 0; jj < lp_bucket_p->list_size; jj++)
            //{
            //    int np = jj / (2 * numblocks * lp_bucket_p->alloc_slices);
            //    if ((jj % (2 * numblocks * lp_bucket_p->alloc_slices)) == 0)
            //    {
            //        printf("poly %d: num slices = %u\n", np, bi[np]);
            //    }
            //    if ((jj % (2 * numblocks)) == 0)
            //    {
            //        int ns = (jj - (np * 2 * numblocks * lp_bucket_p->alloc_slices)) / (2 * numblocks);
            //        printf("\tslice %d: fb_bounds = %u, logp = %u\n", ns,
            //            lp_bucket_p->fb_bounds[np * lp_bucket_p->alloc_slices + ns],
            //            lp_bucket_p->logp[np * lp_bucket_p->alloc_slices + ns]);
            //    }
            //    printf("\t\tblock %d: num = %u\n", jj - (np * 2 * numblocks * lp_bucket_p->alloc_slices),
            //        lp_bucket_p->num[jj]);
            //}
        }
    }
    else
    {
        if (lp_bucket_p->list != NULL)
        {
            lp_bucket_p->num_slices_batch[pnum] = bound_index + 1;
            lp_bucket_p->logp[pnum * lp_bucket_p->alloc_slices + bound_index] = logp;

            //printf("end bucket sieve at j = %u and bound_index = %u\n", j, bound_index + 1);
            //for (j = 0; j < lp_bucket_p->alloc_slices; j++)
            //{
            //    printf("\tslice %d: fb_bounds = %u, logp = %u\n", j,
            //        lp_bucket_p->fb_bounds[pnum * lp_bucket_p->alloc_slices + j],
            //        lp_bucket_p->logp[pnum * lp_bucket_p->alloc_slices + j]);
            //}
            
        }
    }

#elif defined(USE_XLBUCKET)


    // test this idea:
    // when primes start to get much larger than the interval, 
    // define a single "bucket" equal in size to the interval.
    // then bucket "sorting" becomes much easier at the expense
    // of more work in lp_sieve and lp_tdiv.  Since the bucket
    // sort routine is currently the longest pole in the tent,
    // maybe this re-balancing will be a net win.

    // also change the bucket structure.
    // the bucket structure is now:
    // 15 bits for the residue
    // 6 bits for the block (allowing up to 64 blocks per side)
    // 11 bits for the fb diff index.
    // this last removes the need for slices.  We store the difference
    // in index since the last prime that hit this block.  That means
    // at least 1 of every 2048 primes needs to hit a block or the diff
    // will go in the weeds.

    uint32_t lastidp, lastidn;
    uint32_t nump = 0;
    uint32_t numn = 0;

    xlp_bucket_p->numslices = 1;
    xlp_bucket_n->numslices = 1;
    lastidp = j;
    lastidn = j;
    for (k = 0; k < numblocks; k++)
    {
        xlp_bucket_p->sliceid[k] = 0;
        xlp_bucket_n->sliceid[k] = 0;
        xlp_bucket_p->slicenum[k] = 0;
        xlp_bucket_n->slicenum[k] = 0;
    }

//#define DEBUG_XLBUCKET
//#define DEBUG_XLBUCKET2

#ifdef DEBUG_XLBUCKET
    printf("startid=%u, logp=%u\n", j, update_data.logp[j]);
#endif

    if (sign > 0)
    {
        for (; j < bound; j += 16, ptr += 16)
        {
            // new slice?
            //if (xlp_bucket_p->slicenum[xlp_bucket_p->numslices - 1] > 4294967295) // 2015
            //{
            //    xlp_bucket_p->sliceid[xlp_bucket_p->numslices] = lastidp = j;
            //    xlp_bucket_p->slicenum[xlp_bucket_p->numslices] = 0;
            //    xlp_bucket_p->slicelogp[xlp_bucket_p->numslices] = update_data.logp[j];
            //    xlp_bucket_p->numslices++;
            //}


            //ALIGNED_MEM uint64_t b64[8];
            //printf("pid: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08d ", j + k);
            //}
            //printf("\n");
            //printf("p  : ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.prime[j+k]);
            //}
            //printf("\n");
            


            vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i*)ptr);
            mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
            mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
            vroot1 = _mm512_sub_epi32(vroot1, vpval);
            vroot2 = _mm512_sub_epi32(vroot2, vpval);
            vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);
            //vbnum1 = _mm512_srli_epi32(vroot1, 15);
            //vbnum2 = _mm512_srli_epi32(vroot2, 15);

            //printf("r1+: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.firstroots1[j + k]);
            //}
            //printf("\n");
            //printf("r2+: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.firstroots2[j + k]);
            //}
            //printf("\n");
            //printf("r1-: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.prime[j + k] - update_data.firstroots1[j + k]);
            //}
            //printf("\n");
            //printf("r2-: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.prime[j + k] - update_data.firstroots2[j + k]);
            //}
            //printf("\n");

            // index since the last prime to hit the interval
            //velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidp) << 21), vshifted21_index);
            // residue in the block
            //velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            //velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            // block in the interval, must not be greater than 64 or will 
            // overwrite some of the index bits.
            //velement2 = _mm512_or_epi32(velement2, _mm512_slli_epi32(vbnum2, 15));
            //velement1 = _mm512_or_epi32(velement1, _mm512_slli_epi32(vbnum1, 15));

            velement1 = velement2 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidp) << 18), vshifted18_index);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            velement1 = _mm512_or_epi32(velement1, vroot1);
            velement2 = _mm512_or_epi32(velement2, vroot2);

#ifdef DEBUG_XLBUCKET
            
            _mm512_store_epi64(b64, velement1);
            printf("e1: ");
            for (k = 0; k < 8; k++)
            {
                printf("%16lx ", b64[k]);
            }
            printf("\n");
            _mm512_store_epi64(b64, velement2);
            printf("e2: ");
            for (k = 0; k < 8; k++)
            {
                printf("%16lx ", b64[k]);
            }
            printf("\n");
            _mm512_store_epi64(b64, velement3);
            printf("e3: ");
            for (k = 0; k < 8; k++)
            {
                printf("%16lx ", b64[k]);
            }
            printf("\n");
            _mm512_store_epi64(b64, velement4);
            printf("e4: ");
            for (k = 0; k < 8; k++)
            {
                printf("%16lx ", b64[k]);
            }
            printf("\n");
            
#endif

            
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump), mask1, velement1);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump + idx), mask2, velement2);
            idx = _mm_popcnt_u32(mask1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump), mask1, velement1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump + idx), mask2, velement2);
            idx += _mm_popcnt_u32(mask2);

            //printf("%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n",
            //    mask8_1, _mm_popcnt_u32(mask8_1), mask8_2, _mm_popcnt_u32(mask8_2),
            //    mask8_3, _mm_popcnt_u32(mask8_3), mask8_4, _mm_popcnt_u32(mask8_4));
            //
            //printf("xlp_bucket_p->list: ");
            //for (k = 0; k < idx; k++)
            //{
            //    printf("%16lx ", xlp_bucket_p->list[nump + k]);
            //}
            //printf("\n");
            //
            //exit(1);

#ifdef DEBUG_XLBUCKET2
            printf("list: ");
            for (k = 0; k < idx; k++)
            {
                printf("%08x ", xlp_bucket_p->list[nump + k]);
            }
            printf("\n");
#endif

            nump += idx;
            //xlp_bucket_p->slicenum[xlp_bucket_p->numslices - 1] += idx;

#ifdef DEBUG_XLBUCKET
            printf("now nump=%u, logp=%u, numslices=%u, slicenum=%u, sliceid=%u\n",
                nump, xlp_bucket_p->slicelogp[xlp_bucket_p->numslices - 1],
                xlp_bucket_p->numslices, xlp_bucket_p->slicenum[xlp_bucket_p->numslices - 1],
                xlp_bucket_p->sliceid[xlp_bucket_p->numslices - 1]);
#endif

            // new slice?
            //if (xlp_bucket_n->slicenum[xlp_bucket_n->numslices - 1] > 4294967295) // 2015
            //{
            //    xlp_bucket_n->sliceid[xlp_bucket_n->numslices] = lastidn = j;
            //    xlp_bucket_n->slicenum[xlp_bucket_n->numslices] = 0;
            //    xlp_bucket_n->slicelogp[xlp_bucket_n->numslices] = update_data.logp[j];
            //    xlp_bucket_n->numslices++;
            //}

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            //vbnum1 = _mm512_srli_epi32(vroot1, 15);
            //vbnum2 = _mm512_srli_epi32(vroot2, 15);
            //
            //// index since the last prime to hit the interval
            //velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidn) << 21), vshifted21_index);
            //// residue in the block
            //velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            //velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            //// block in the interval, must not be greater than 64 or will 
            //// overwrite some of the index bits.
            //velement2 = _mm512_or_epi32(velement2, _mm512_slli_epi32(vbnum2, 15));
            //velement1 = _mm512_or_epi32(velement1, _mm512_slli_epi32(vbnum1, 15));
            //mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            //mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

            //vbase1 = _mm512_add_epi64(_mm512_set1_epi64((uint64_t)((j - lastidn)) << 21), v64shifted21_index_e);
            //vbase1 = _mm512_or_epi64(_mm512_set1_epi64((uint64_t)update_data.logp[j] << 53), vbase1);
            //vbase2 = _mm512_add_epi64(_mm512_set1_epi64((uint64_t)((j - lastidn)) << 21), v64shifted21_index_o);
            //vbase2 = _mm512_or_epi64(_mm512_set1_epi64((uint64_t)update_data.logp[j] << 53), vbase2);
            //velement1 = _mm512_maskz_shuffle_epi32(0x5555, vroot1, 0xe4);
            //velement3 = _mm512_maskz_shuffle_epi32(0x5555, vroot1, 0xb1);
            //velement2 = _mm512_maskz_shuffle_epi32(0x5555, vroot2, 0xe4);
            //velement4 = _mm512_maskz_shuffle_epi32(0x5555, vroot2, 0xb1);
            //mask8_1 = _mm512_cmp_epu64_mask(velement1, vinterval64, _MM_CMPINT_LT);
            //mask8_2 = _mm512_cmp_epu64_mask(velement2, vinterval64, _MM_CMPINT_LT);
            //mask8_3 = _mm512_cmp_epu64_mask(velement3, vinterval64, _MM_CMPINT_LT);
            //mask8_4 = _mm512_cmp_epu64_mask(velement4, vinterval64, _MM_CMPINT_LT);
            //velement1 = _mm512_or_epi64(velement1, vbase1);
            //velement3 = _mm512_or_epi64(velement3, vbase2);
            //velement2 = _mm512_or_epi64(velement2, vbase1);
            //velement4 = _mm512_or_epi64(velement4, vbase2);

            velement1 = velement2 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidn) << 18), vshifted18_index);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            velement1 = _mm512_or_epi32(velement1, vroot1);
            velement2 = _mm512_or_epi32(velement2, vroot2);

#ifdef DEBUG_XLBUCKET
            printf("%08x (%u hits)\n%08x (%u hits)\n",
                mask1, _mm_popcnt_u32(mask1), mask2, _mm_popcnt_u32(mask2));
            //_mm512_store_epi32(b1, velement1);
            //printf("e1: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", b1[k]);
            //}
            //printf("\n");
            //_mm512_store_epi32(b1, velement2);
            //printf("e2: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", b1[k]);
            //}
            //printf("\n");
#endif

            //idx = _mm_popcnt_u32(mask1);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn), mask1, velement1);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn + idx), mask2, velement2);

            //idx = _mm_popcnt_u32(mask8_1);
            //_mm512_mask_compressstoreu_epi64((__m512i*)(xlp_bucket_n->list + numn), mask8_1, velement1);
            //_mm512_mask_compressstoreu_epi64((__m512i*)(xlp_bucket_n->list + numn + idx), mask8_2, velement2);
            //idx += _mm_popcnt_u32(mask8_2);
            //_mm512_mask_compressstoreu_epi64((__m512i*)(xlp_bucket_n->list + numn + idx), mask8_3, velement3);
            //idx += _mm_popcnt_u32(mask8_3);
            //_mm512_mask_compressstoreu_epi64((__m512i*)(xlp_bucket_n->list + numn + idx), mask8_4, velement4);
            //idx += _mm_popcnt_u32(mask8_4);

            //printf("%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n",
            //    mask8_1, _mm_popcnt_u32(mask8_1), mask8_2, _mm_popcnt_u32(mask8_2),
            //    mask8_3, _mm_popcnt_u32(mask8_3), mask8_4, _mm_popcnt_u32(mask8_4));
            //
            //printf("xlp_bucket_n->list: ");
            //for (k = 0; k < idx; k++)
            //{
            //    printf("%16lx ", xlp_bucket_n->list[numn + k]);
            //}
            //printf("\n");

            //idx += _mm_popcnt_u32(mask2);

            idx = _mm_popcnt_u32(mask1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn), mask1, velement1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn + idx), mask2, velement2);
            idx += _mm_popcnt_u32(mask2);

#ifdef DEBUG_XLBUCKET2
            printf("list: ");
            for (k = 0; k < idx; k++)
            {
                printf("%08x ", xlp_bucket_n->list[numn + k]);
            }
            printf("\n");
#endif

            numn += idx;
            //xlp_bucket_n->slicenum[xlp_bucket_n->numslices - 1] += idx;

#ifdef DEBUG_XLBUCKET
            printf("now numn=%u, logp=%u, numslices=%u, slicenum=%u, sliceid=%u\n",
                numn, xlp_bucket_n->slicelogp[xlp_bucket_n->numslices - 1],
                xlp_bucket_n->numslices, xlp_bucket_n->slicenum[xlp_bucket_n->numslices - 1],
                xlp_bucket_n->sliceid[xlp_bucket_n->numslices - 1]);
#endif
        }
    }
    else
    {
        for (; j < bound; j += 16, ptr += 16)
        {
            // new slice?
            //if (xlp_bucket_p->slicenum[xlp_bucket_p->numslices - 1] > 4294967295) // 2015
            //{
            //    xlp_bucket_p->sliceid[xlp_bucket_p->numslices] = lastidp = j;
            //    xlp_bucket_p->slicenum[xlp_bucket_p->numslices] = 0;
            //    xlp_bucket_p->slicelogp[xlp_bucket_p->numslices] = update_data.logp[j];
            //    xlp_bucket_p->numslices++;
            //}

            //ALIGNED_MEM uint64_t b64[8];
            //printf("pid: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08d ", j + k);
            //}
            //printf("\n");
            //printf("p  : ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.prime[j + k]);
            //}
            //printf("\n");

            vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
            vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
            vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));
            vpval = _mm512_load_epi32((__m512i*)ptr);
            vroot1 = _mm512_add_epi32(vroot1, vpval);
            vroot2 = _mm512_add_epi32(vroot2, vpval);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);

            //printf("r1+: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.firstroots1[j + k]);
            //}
            //printf("\n");
            //printf("r2+: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.firstroots2[j + k]);
            //}
            //printf("\n");
            //printf("r1-: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.prime[j + k] - update_data.firstroots1[j + k]);
            //}
            //printf("\n");
            //printf("r2-: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", update_data.prime[j + k] - update_data.firstroots2[j + k]);
            //}
            //printf("\n");

            //vbnum1 = _mm512_srli_epi32(vroot1, 15);
            //vbnum2 = _mm512_srli_epi32(vroot2, 15);
            //// index since the last prime to hit the interval
            //velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidp) << 21), vshifted21_index);
            //// residue in the block
            //velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            //velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            //// block in the interval, must not be greater than 64 or will 
            //// overwrite some of the index bits.
            //velement2 = _mm512_or_epi32(velement2, _mm512_slli_epi32(vbnum2, 15));
            //velement1 = _mm512_or_epi32(velement1, _mm512_slli_epi32(vbnum1, 15));
            //mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            //mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            
            velement1 = velement2 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidp) << 18), vshifted18_index);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            velement1 = _mm512_or_epi32(velement1, vroot1);
            velement2 = _mm512_or_epi32(velement2, vroot2);


#ifdef DEBUG_XLBUCKET
            printf("%08x (%u hits)\n%08x (%u hits)\n",
                mask1, _mm_popcnt_u32(mask1), mask2, _mm_popcnt_u32(mask2));
            //_mm512_store_epi32(b1, velement1);
            //printf("e1: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", b1[k]);
            //}
            //printf("\n");
            //_mm512_store_epi32(b1, velement2);
            //printf("e2: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", b1[k]);
            //}
            //printf("\n");
#endif

            //idx = _mm_popcnt_u32(mask1);
            ////_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump), mask1, velement1);
            ////_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump + idx), mask2, velement2);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump), mask1, vroot1);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump + idx), mask2, vroot2);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->fbid + nump), mask1,
            //    _mm512_add_epi32(_mm512_set1_epi32(j - lastidp), vindex));
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->fbid + nump + idx), mask2,
            //    _mm512_add_epi32(_mm512_set1_epi32(j - lastidp), vindex));
            //idx += _mm_popcnt_u32(mask2);

            idx = _mm_popcnt_u32(mask1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump), mask1, velement1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_p->list + nump + idx), mask2, velement2);
            idx += _mm_popcnt_u32(mask2);

            //printf("%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n",
            //    mask8_1, _mm_popcnt_u32(mask8_1), mask8_2, _mm_popcnt_u32(mask8_2),
            //    mask8_3, _mm_popcnt_u32(mask8_3), mask8_4, _mm_popcnt_u32(mask8_4));
            //
            //printf("xlp_bucket_p->list: ");
            //for (k = 0; k < idx; k++)
            //{
            //    printf("%16lx ", xlp_bucket_p->list[nump + k]);
            //}
            //printf("\n");

#ifdef DEBUG_XLBUCKET2
            printf("list: ");
            for (k = 0; k < idx; k++)
            {
                printf("%08x ", xlp_bucket_p->list[nump + k]);
            }
            printf("\n");
#endif

            nump += idx;
            //xlp_bucket_p->slicenum[xlp_bucket_p->numslices - 1] += idx;

#ifdef DEBUG_XLBUCKET
            printf("now nump=%u, logp=%u, numslices=%u, slicenum=%u, sliceid=%u\n",
                nump, xlp_bucket_p->slicelogp[xlp_bucket_p->numslices - 1],
                xlp_bucket_p->numslices, xlp_bucket_p->slicenum[xlp_bucket_p->numslices - 1],
                xlp_bucket_p->sliceid[xlp_bucket_p->numslices - 1]);
#endif

            // new slice?
            //if (xlp_bucket_n->slicenum[xlp_bucket_n->numslices - 1] > 4294967295) // 2015
            //{
            //    xlp_bucket_n->sliceid[xlp_bucket_n->numslices] = lastidn = j;
            //    xlp_bucket_n->slicenum[xlp_bucket_n->numslices] = 0;
            //    xlp_bucket_n->slicelogp[xlp_bucket_n->numslices] = update_data.logp[j];
            //    xlp_bucket_n->numslices++;
            //}

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            //vbnum1 = _mm512_srli_epi32(vroot1, 15);
            //vbnum2 = _mm512_srli_epi32(vroot2, 15);
            //// index since the last prime to hit the interval
            //velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidn) << 21), vshifted21_index);
            //// residue in the block
            //velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            //velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            //// block in the interval, must not be greater than 64 or will 
            //// overwrite some of the index bits.
            //velement2 = _mm512_or_epi32(velement2, _mm512_slli_epi32(vbnum2, 15));
            //velement1 = _mm512_or_epi32(velement1, _mm512_slli_epi32(vbnum1, 15));
            //mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            //mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

            velement1 = velement2 = _mm512_add_epi32(_mm512_set1_epi32((j - lastidn) << 18), vshifted18_index);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            velement1 = _mm512_or_epi32(velement1, vroot1);
            velement2 = _mm512_or_epi32(velement2, vroot2);

#ifdef DEBUG_XLBUCKET
            printf("%08x (%u hits)\n%08x (%u hits)\n",
                mask1, _mm_popcnt_u32(mask1), mask2, _mm_popcnt_u32(mask2));
            //_mm512_store_epi32(b1, velement1);
            //printf("e1: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", b1[k]);
            //}
            //printf("\n");
            //_mm512_store_epi32(b1, velement2);
            //printf("e2: ");
            //for (k = 0; k < 16; k++)
            //{
            //    printf("%08x ", b1[k]);
            //}
            //printf("\n");
#endif

            //idx = _mm_popcnt_u32(mask1);
            ////_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn), mask1, velement1);
            ////_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn + idx), mask2, velement2);
            //
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn), mask1, vroot1);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn + idx), mask2, vroot2);
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->fbid + numn), mask1,
            //    _mm512_add_epi32(_mm512_set1_epi32(j - lastidn), vindex));
            //_mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->fbid + numn + idx), mask2,
            //    _mm512_add_epi32(_mm512_set1_epi32(j - lastidn), vindex));

            //idx += _mm_popcnt_u32(mask2);

            idx = _mm_popcnt_u32(mask1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn), mask1, velement1);
            _mm512_mask_compressstoreu_epi32((__m512i*)(xlp_bucket_n->list + numn + idx), mask2, velement2);
            idx += _mm_popcnt_u32(mask2);

            //printf("%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n%08x (%u hits)\n",
            //    mask8_1, _mm_popcnt_u32(mask8_1), mask8_2, _mm_popcnt_u32(mask8_2),
            //    mask8_3, _mm_popcnt_u32(mask8_3), mask8_4, _mm_popcnt_u32(mask8_4));
            //
            //printf("xlp_bucket_n->list: ");
            //for (k = 0; k < idx; k++)
            //{
            //    printf("%16lx ", xlp_bucket_n->list[numn + k]);
            //}
            //printf("\n");

#ifdef DEBUG_XLBUCKET2
            printf("list: ");
            for (k = 0; k < idx; k++)
            {
                printf("%08x ", xlp_bucket_n->list[numn + k]);
            }
            printf("\n");
#endif

            numn += idx;
            //xlp_bucket_n->slicenum[xlp_bucket_n->numslices - 1] += idx;

#ifdef DEBUG_XLBUCKET
            printf("now numn=%u, logp=%u, numslices=%u, slicenum=%u, sliceid=%u\n",
                numn, xlp_bucket_n->slicelogp[xlp_bucket_n->numslices - 1], 
                xlp_bucket_n->numslices, xlp_bucket_n->slicenum[xlp_bucket_n->numslices - 1],
                xlp_bucket_n->sliceid[xlp_bucket_n->numslices - 1]);
#endif
        }
    }

    xlp_bucket_p->slicenum[0] = 0;
    xlp_bucket_n->slicenum[0] = 0;
    xlp_bucket_p->sliceid[0] = nump;
    xlp_bucket_n->sliceid[0] = numn;
    // sort the lists by block then count the elements in each block
    
    //qsort(xlp_bucket_p->list, nump, sizeof(uint32_t), qcomp_xlbucket);
    //qsort(xlp_bucket_n->list, numn, sizeof(uint32_t), qcomp_xlbucket);
    ////printf("sorting %u entries in xlp_bucket_p\n", nump);
    ////printf("sorting %u entries in xlp_bucket_n\n", numn);
    //
    //uint32_t curr_blk = 0;
    //xlp_bucket_p->sliceid[curr_blk] = 0;
    //xlp_bucket_p->slicenum[curr_blk] = 0;
    //for (k = 0; k < nump; k++)
    //{
    //    if (((xlp_bucket_p->list[k] >> 15) & 0x7) > curr_blk)
    //    {
    //        //printf("block %u has %u entries in xlp_bucket_p starting at offset %u\n", 
    //        //    curr_blk, xlp_bucket_p->sliceid[curr_blk], xlp_bucket_p->slicenum[curr_blk]);
    //        curr_blk++;
    //        xlp_bucket_p->sliceid[curr_blk] = 0;
    //        xlp_bucket_p->slicenum[curr_blk] = k;
    //    }
    //    xlp_bucket_p->sliceid[curr_blk]++;
    //}
    ////printf("block %u has %u entries in xlp_bucket_p starting at offset %u\n",
    ////    curr_blk, xlp_bucket_p->sliceid[curr_blk], xlp_bucket_p->slicenum[curr_blk]);
    //curr_blk++;
    //xlp_bucket_p->sliceid[curr_blk] = 0;
    //xlp_bucket_p->slicenum[curr_blk] = k;
    //
    //
    //curr_blk = 0;
    //xlp_bucket_n->sliceid[curr_blk] = 0;
    //xlp_bucket_n->slicenum[curr_blk] = 0;
    //for (k = 0; k < numn; k++)
    //{
    //    if (((xlp_bucket_n->list[k] >> 15) & 0x7) > curr_blk)
    //    {
    //        //printf("block %u has %u entries in xlp_bucket_n starting at offset %u\n",
    //        //    curr_blk, xlp_bucket_n->sliceid[curr_blk], xlp_bucket_n->slicenum[curr_blk]);
    //        curr_blk++;
    //        xlp_bucket_n->sliceid[curr_blk] = 0;
    //        xlp_bucket_n->slicenum[curr_blk] = k;
    //    }
    //    xlp_bucket_n->sliceid[curr_blk]++;
    //}
    ////printf("block %u has %u entries in xlp_bucket_p starting at offset %u\n",
    ////    curr_blk, xlp_bucket_n->sliceid[curr_blk], xlp_bucket_n->slicenum[curr_blk]);
    //curr_blk++;
    //xlp_bucket_n->sliceid[curr_blk] = 0;
    //xlp_bucket_n->slicenum[curr_blk] = k;
    

#endif


    
    if (lp_bucket_p->list != NULL)
    {
        // all done bucket sieving, record final logp and number of slices.
#if defined (USE_BATCHPOLY_X2)
        
#else
        lp_bucket_p->num_slices = bound_index + 1;
        slicelogp_ptr[bound_index] = logp;

        
        //for (k = 0; k < dconf->poly_batchsize; k++)
        //{
        //    int jj;
        //    printf("poly %d: num slices = %u\n", k, lp_bucket_p->num_slices);
        //    //for (jj = 0; jj < lp_bucket_p->alloc_slices; jj++)
        //    //{
        //    //    printf("\tslice %d: fb_bounds = %u, logp = %u\n", jj,
        //    //        lp_bucket_p->fb_bounds[k * lp_bucket_p->alloc_slices + jj],
        //    //        lp_bucket_p->logp[k * lp_bucket_p->alloc_slices + jj]);
        //    //    int i;
        //    //    for (i = 0; i < numblocks; i++)
        //    //    {
        //    //        printf("\t\tblock %d: pnum = %u\n", i,
        //    //            lp_bucket_p->num[jj * lp_bucket_p->alloc_slices + i]);
        //    //    }
        //    //    for (; i < 2 * numblocks; i++)
        //    //    {
        //    //        printf("\t\tblock %d: nnum = %u\n", i - numblocks,
        //    //            lp_bucket_p->num[jj * lp_bucket_p->alloc_slices + i]);
        //    //    }
        //    //}
        //
        //    for (jj = 0; jj < lp_bucket_p->list_size; jj++)
        //    {
        //        if ((jj % (2 * numblocks)) == 0)
        //        {
        //            int ns = jj / (2 * numblocks);
        //            printf("\tslice %d: fb_bounds = %u, logp = %u\n", ns,
        //                lp_bucket_p->fb_bounds[ns], lp_bucket_p->logp[ns]);
        //        }
        //        printf("\t\tblock %d: num = %u\n", jj,
        //            lp_bucket_p->num[jj]);
        //
        //    }
        //}
#endif
        
    }

    return;
}

#define COMPUTE_NEXT_ROOTS_BATCH_P(i) \
        root1 = (int)root1 - rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root2 = (int)root2 - rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root1 += ((root1 >> 31) * prime); \
        root2 += ((root2 >> 31) * prime);

#define COMPUTE_NEXT_ROOTS_BATCH_N(i) \
        root1 = (int)root1 + rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root2 = (int)root2 + rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root1 -= ((root1 >= prime) * prime); \
        root2 -= ((root2 >= prime) * prime); \


#ifdef USE_BATCHPOLY
void nextRoots_32k_knl_polybatch(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    int *rootupdates = dconf->rootupdates;
    update_t update_data = dconf->update_data;

    uint32_t startprime = 2;
    uint32_t bound = sconf->factor_base->B;
    char *nu = dconf->curr_poly->nu;
    char *gray = dconf->curr_poly->gray;
    int numB = dconf->numB;
    uint32_t poly_offset = 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

    lp_bucket *lp_bucket_p = dconf->buckets;
    uint32_t med_B = sconf->factor_base->med_B;
    uint32_t large_B = sconf->factor_base->large_B;

    uint32_t j, interval;
    int k, numblocks, idx;
    uint32_t root1, root2, nroot1, nroot2, prime;

    int bound_index = 0;
    int check_bound = BUCKET_ALLOC / 2 - 1;
    uint32_t bound_val = med_B;
    uint32_t *numptr_p, *numptr_n, *sliceptr_p, *sliceptr_n;

    uint32_t *bptr;
    int bnum, room;
    uint8_t logp = 0;
    int p;

    __m512i vshifted_index = _mm512_setr_epi32(
        0, 65536, 131072, 196608, 
        262144,327680, 393216, 458752, 
        524288,589824, 655360, 720896, 
        786432,851968, 917504, 983040);

    __m512i vnroot1, vnroot2, vzero = _mm512_setzero_epi32(), vone = _mm512_set1_epi32(1);
    __m512i vprime, vroot1, vroot2, vpval, vbnum1, vbnum2, vinterval, vindex, vnum, vaddr;
    __m512i vtmproot1, vtmproot2, vtmp;
    __m512i velement1, velement2, vbsize, vidx, vblockm1 = _mm512_set1_epi32(32767);
    __mmask16 mask1, mask2, mconflictfree, munique;
    __mmask16 mNotProc;

    uint32_t gmask[64];
    
    ALIGNED_MEM uint32_t e1[32];
    ALIGNED_MEM uint32_t e2[32];
    ALIGNED_MEM uint32_t b1[32];
    ALIGNED_MEM uint32_t b2[32];


#if 0
    // good for up to 16 polynomial batches...
    ALIGNED_MEM uint32_t prootstore1[256];
    ALIGNED_MEM uint32_t prootstore2[256];
    ALIGNED_MEM uint32_t nrootstore1[256];
    ALIGNED_MEM uint32_t nrootstore2[256];
    int num_prootstore1;
    int num_prootstore2;
    int num_nrootstore1;
    int num_nrootstore2;
#endif

    numblocks = sconf->num_blocks;
    interval = numblocks << 15;
    vinterval = _mm512_set1_epi32(interval);

    if (lp_bucket_p->alloc_slices != 0) // != NULL)
    {
        lp_bucket_p->fb_bounds[0] = med_B;

        sliceptr_p = lp_bucket_p->list;
        sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

        numptr_p = lp_bucket_p->num;
        numptr_n = lp_bucket_p->num + numblocks;

        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
            numptr_p[j] = 0;

        lp_bucket_p->num_slices = 0;

    }
    else
    {
        sliceptr_p = NULL;
        sliceptr_n = NULL;
        numptr_p = NULL;
        numptr_n = NULL;
    }

    for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
    {
        if (gray[numB + p] > 0)
        {
            gmask[p] = 0xffff;
        }
        else
        {
            gmask[p] = 0;
        }
    }
   
    bound_index = 0;
    bound_val = med_B;
    check_bound = med_B + BUCKET_ALLOC / 2;

    logp = update_data.logp[med_B];
    for (j = med_B; j < large_B; j += 16)
    {
        
        CHECK_NEW_SLICE_BATCH(j);        
        
        vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
        vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
        vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));
        
#if 1
        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {

            __mmask16 gm = gmask[p];

            vpval = _mm512_load_epi32((__m512i*)(&rootupdates[(nu[numB + p] - 1) * bound + j]));

            //vpval = _mm512_mask_add_epi32(vpval, gm, vpval, _mm512_set1_epi32(0x80000000));
            //vtmp = _mm512_mask_add_epi32(vprime, gm, vprime, _mm512_set1_epi32(0x80000000));
            //
            //vroot1 = _mm512_add_epi32(vroot1, vpval);
            //vroot2 = _mm512_add_epi32(vroot2, vpval);
            //mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            //mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            //vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vtmp);
            //vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vtmp);


            //if (gray[numB + p] > 0)
            mask1 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot1, _MM_CMPINT_GT);
            mask2 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot2, _MM_CMPINT_GT);
            vroot1 = _mm512_mask_sub_epi32(vroot1, gm, vroot1, vpval);
            vroot2 = _mm512_mask_sub_epi32(vroot2, gm, vroot2, vpval);
            vroot1 = _mm512_mask_add_epi32(vroot1, gm & mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, gm & mask2, vroot2, vprime);
            
            //else
            vroot1 = _mm512_mask_add_epi32(vroot1, ~gm, vroot1, vpval);
            vroot2 = _mm512_mask_add_epi32(vroot2, ~gm, vroot2, vpval);
            mask1 |= _mm512_mask_cmp_epu32_mask(~gm, vroot1, vprime, _MM_CMPINT_GE);
            mask2 |= _mm512_mask_cmp_epu32_mask(~gm, vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, ~gm & mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, ~gm & mask2, vroot2, vprime);


            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);

            // we loop over roots below, so have to compute these before we
            // clobber them.
            vnroot1 = _mm512_sub_epi32(vprime, vroot1);
            vnroot2 = _mm512_sub_epi32(vprime, vroot2);


            mask1 = 0xffff;
            mask2 = 0xffff;
            idx = 32;

            _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot1, 15));
            _mm512_store_epi32((__m512i*)e1,
                _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

            _mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vroot2, 15));
            _mm512_store_epi32((__m512i*)(&(e1[16])),
                _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

            vtmproot1 = vroot1;
            vtmproot2 = vroot2;

            do
            {
                SCATTER_COMPRESSED_VECTOR_P;

                vtmproot1 = _mm512_mask_add_epi32(vtmproot1, mask1, vtmproot1, vprime);
                vtmproot2 = _mm512_mask_add_epi32(vtmproot2, mask2, vtmproot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vtmproot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vtmproot2, vinterval, _MM_CMPINT_LT);

                if ((mask1 | mask2) == 0)
                {
                    break;
                }

                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vtmproot1, 15));
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot1, vblockm1)));

                idx = _mm_popcnt_u32(mask1);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vtmproot2, 15));
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot2, vblockm1)));

                idx += _mm_popcnt_u32(mask2);
            } while (idx > 0);

            mask1 = 0xffff;
            mask2 = 0xffff;
            idx = 32;

            _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
            _mm512_store_epi32((__m512i*)e1,
                _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

            _mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vnroot2, 15));
            _mm512_store_epi32((__m512i*)(&(e1[16])),
                _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

            do
            {
                SCATTER_COMPRESSED_VECTOR_N;

                vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

                if ((mask1 | mask2) == 0)
                {
                    break;
                }

                _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vnroot1, 15));
                _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                idx = _mm_popcnt_u32(mask1);
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vnroot2, 15));
                _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                idx += _mm_popcnt_u32(mask2);
            } while (idx > 0);


            // advance pointers
            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;
        }

#else

        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            if (gray[numB + p] > 0)
            {
                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));
                mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_sub_epi32(vroot1, vpval);
                vroot2 = _mm512_sub_epi32(vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);                
                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
                
                // we loop over roots below, so have to compute these before we
                // clobber them.
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                
#if 0
                // do the first 16
                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif
                
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask1;
                while (mNotProc != 0)   
                {     
					mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
					vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
					mNotProc &= (~mconflictfree);
					vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
                  
                // then the rest as they start to drop out of the interval
                vtmproot1 = _mm512_add_epi32(vroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vtmproot1, vinterval, _MM_CMPINT_LT);
                while (mask1 > 0)
                {
                    __mmask16 m = mask1;

                    _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vtmproot1, 15));
                    _mm512_store_epi32((__m512i *)e1, 
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot1, vblockm1)));

                    while (m > 0)  
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                        *bptr = e1[idx];
                        numptr_p[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vtmproot1 = _mm512_mask_add_epi32(vtmproot1, mask1, vtmproot1, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vtmproot1, vinterval, _MM_CMPINT_LT);
                }
                
                // do the first 16 of mask2
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));              
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

                vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask2;
                while (mNotProc != 0)   
                {     
					mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
					vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
					mNotProc &= (~mconflictfree);
					vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask2, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique & mask2, vbnum2, vnum, _MM_SCALE_4);

                // then the rest as they start to drop out of the interval
                vtmproot2 = _mm512_add_epi32(vroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vtmproot2, vinterval, _MM_CMPINT_LT);
                while (mask2 > 0)
                {
                    __mmask16 m = mask2;

                    _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vtmproot2, 15));
                    _mm512_store_epi32((__m512i *)e1, 
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot2, vblockm1)));

                    while (m > 0)  
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                        *bptr = e1[idx];
                        numptr_p[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vtmproot2 = _mm512_mask_add_epi32(vtmproot2, mask2, vtmproot2, vprime);
                    mask2 = _mm512_cmp_epu32_mask(vtmproot2, vinterval, _MM_CMPINT_LT);
                }
                
                // do the first 16 of nmask1
                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif
                
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask1;
                while (mNotProc != 0)   
                {     
					mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
					vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
					mNotProc &= (~mconflictfree);
					vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_n, mask1, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_n, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
                  
                // then the rest as they start to drop out of the interval
                vnroot1 = _mm512_add_epi32(vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                while (mask1 > 0)
                {
                    __mmask16 m = mask1;

                    _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot1, 15));
                    _mm512_store_epi32((__m512i *)e1, 
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                    while (m > 0)  
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                        *bptr = e1[idx];
                        numptr_n[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                }
                
                // do the first 16 of nmask2
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1));              
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

                vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask2;
                while (mNotProc != 0)   
                {     
					mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
					vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
					mNotProc &= (~mconflictfree);
					vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_n, mask2, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_n, munique & mask2, vbnum2, vnum, _MM_SCALE_4);

                // then the rest as they start to drop out of the interval
                vnroot2 = _mm512_add_epi32(vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                while (mask2 > 0)
                {
                    __mmask16 m = mask2;

                    _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot2, 15));
                    _mm512_store_epi32((__m512i *)e1, 
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                    while (m > 0)  
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                        *bptr = e1[idx];
                        numptr_n[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                    mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                }

#else

                mask1 = 0xffff;
                mask2 = 0xffff;
                idx = 32;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

                _mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i*)(&(e1[16])),
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

                vtmproot1 = vroot1;
                vtmproot2 = vroot2;

                do
                {
                    SCATTER_COMPRESSED_VECTOR_P;

                    vtmproot1 = _mm512_mask_add_epi32(vtmproot1, mask1, vtmproot1, vprime);
                    vtmproot2 = _mm512_mask_add_epi32(vtmproot2, mask2, vtmproot2, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vtmproot1, vinterval, _MM_CMPINT_LT);
                    mask2 = _mm512_cmp_epu32_mask(vtmproot2, vinterval, _MM_CMPINT_LT);

                    if ((mask1 | mask2) == 0)
                    {
                        break;
                    }

                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vtmproot1, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot1, vblockm1)));

                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vtmproot2, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot2, vblockm1)));

                    idx += _mm_popcnt_u32(mask2);
                } while (idx > 0);

                mask1 = 0xffff;
                mask2 = 0xffff;
                idx = 32;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                _mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vnroot2, 15));
                _mm512_store_epi32((__m512i*)(&(e1[16])),
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                do
                {
                    SCATTER_COMPRESSED_VECTOR_N;

                    vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                    vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                    mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

                    if ((mask1 | mask2) == 0)
                    {
                        break;
                    }

                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vnroot1, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vnroot2, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                    idx += _mm_popcnt_u32(mask2);
                } while (idx > 0);

#endif

            }
            else
            {
                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));
                vroot1 = _mm512_add_epi32(vroot1, vpval);
                vroot2 = _mm512_add_epi32(vroot2, vpval);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);               
                
                // we loop over roots below, so have to compute these before we
                // clobber them.
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);

#if 0
                // do the first 16
                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask1;
                while (mNotProc != 0)
                {
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                    mNotProc &= (~mconflictfree);
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique& mask1, vbnum1, vnum, _MM_SCALE_4);

                // then the rest as they start to drop out of the interval
                vtmproot1 = _mm512_add_epi32(vroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vtmproot1, vinterval, _MM_CMPINT_LT);
                while (mask1 > 0)
                {
                    __mmask16 m = mask1;

                    _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vtmproot1, 15));
                    _mm512_store_epi32((__m512i*)e1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot1, vblockm1)));

                    while (m > 0)
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                        *bptr = e1[idx];
                        numptr_p[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vtmproot1 = _mm512_mask_add_epi32(vtmproot1, mask1, vtmproot1, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vtmproot1, vinterval, _MM_CMPINT_LT);
                }

                // do the first 16
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

                vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask2;
                while (mNotProc != 0)
                {
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                    mNotProc &= (~mconflictfree);
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_p, mask2, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_p, munique& mask2, vbnum2, vnum, _MM_SCALE_4);

                // then the rest as they start to drop out of the interval
                vtmproot2 = _mm512_add_epi32(vroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vtmproot2, vinterval, _MM_CMPINT_LT);
                while (mask2 > 0)
                {
                    __mmask16 m = mask2;

                    _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vtmproot2, 15));
                    _mm512_store_epi32((__m512i*)e1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot2, vblockm1)));

                    while (m > 0)
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                        *bptr = e1[idx];
                        numptr_p[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vtmproot2 = _mm512_mask_add_epi32(vtmproot2, mask2, vtmproot2, vprime);
                    mask2 = _mm512_cmp_epu32_mask(vtmproot2, vinterval, _MM_CMPINT_LT);
                }

                // do the first 16
                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask1;
                while (mNotProc != 0)
                {
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                    mNotProc &= (~mconflictfree);
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_n, mask1, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_n, munique& mask1, vbnum1, vnum, _MM_SCALE_4);

                // then the rest as they start to drop out of the interval
                vnroot1 = _mm512_add_epi32(vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                while (mask1 > 0)
                {
                    __mmask16 m = mask1;

                    _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
                    _mm512_store_epi32((__m512i*)e1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                    while (m > 0)
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                        *bptr = e1[idx];
                        numptr_n[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                }

                // do the first 16
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                vidx = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1));
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

                vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask2;
                while (mNotProc != 0)
                {
                    mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                    vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                    mNotProc &= (~mconflictfree);
                    vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                }
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_i32scatter_epi32(sliceptr_n, mask2, vaddr, vidx, _MM_SCALE_4);
                vnum = _mm512_add_epi32(vnum, vone);
                _mm512_mask_i32scatter_epi32(numptr_n, munique& mask2, vbnum2, vnum, _MM_SCALE_4);

                // then the rest as they start to drop out of the interval
                vnroot2 = _mm512_add_epi32(vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                while (mask2 > 0)
                {
                    __mmask16 m = mask2;

                    _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot2, 15));
                    _mm512_store_epi32((__m512i*)e1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                    while (m > 0)
                    {
                        idx = _trail_zcnt(m);
                        bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                        *bptr = e1[idx];
                        numptr_n[b1[idx]]++;
                        m = _reset_lsb(m);
                    }

                    vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                    mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                }

#else
                mask1 = 0xffff;
                mask2 = 0xffff;
                idx = 32;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

                _mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i*)(&(e1[16])),
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

                vtmproot1 = vroot1;
                vtmproot2 = vroot2;

                do
                {
                    SCATTER_COMPRESSED_VECTOR_P;

                    vtmproot1 = _mm512_mask_add_epi32(vtmproot1, mask1, vtmproot1, vprime);
                    vtmproot2 = _mm512_mask_add_epi32(vtmproot2, mask2, vtmproot2, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vtmproot1, vinterval, _MM_CMPINT_LT);
                    mask2 = _mm512_cmp_epu32_mask(vtmproot2, vinterval, _MM_CMPINT_LT);

                    if ((mask1 | mask2) == 0)
                    {
                        break;
                    }

                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vtmproot1, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot1, vblockm1)));

                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vtmproot2, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vtmproot2, vblockm1)));

                    idx += _mm_popcnt_u32(mask2);
                } while (idx > 0);

                mask1 = 0xffff;
                mask2 = 0xffff;
                idx = 32;

                _mm512_store_epi32((__m512i*)b1, _mm512_srli_epi32(vnroot1, 15));
                _mm512_store_epi32((__m512i*)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                _mm512_store_epi32((__m512i*)(&(b1[16])), _mm512_srli_epi32(vnroot2, 15));
                _mm512_store_epi32((__m512i*)(&(e1[16])),
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                do
                {
                    SCATTER_COMPRESSED_VECTOR_N;

                    vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                    vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                    mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                    mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

                    if ((mask1 | mask2) == 0)
                    {
                        break;
                    }

                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, _mm512_srli_epi32(vnroot1, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, _mm512_srli_epi32(vnroot2, 15));
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2,
                        _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                    idx += _mm_popcnt_u32(mask2);
                } while (idx > 0);

#endif
                

                
            }



            // advance pointers
            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;
        }



#endif
        
        _mm512_store_epi32((__m512i *)(&update_data.firstroots1[j]), vroot1);
        _mm512_store_epi32((__m512i *)(&update_data.firstroots2[j]), vroot2);

        // reset pointers
        sliceptr_p -= p * poly_offset * BUCKET_ALLOC;
        sliceptr_n -= p * poly_offset * BUCKET_ALLOC;
        numptr_p -= p * poly_offset;
        numptr_n -= p * poly_offset;

    }

    logp = update_data.logp[j - 1];
    for (j = large_B; j < bound; j += 16)
    {
        int p;

        CHECK_NEW_SLICE_BATCH(j);

        vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
        vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
        vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));
                
        // todo:
        // compute all roots for N polynomials.
        // compress-store those less than interval to temporary location.
        // adjust vbnum by the poly offset so scatter/gather grabs from
        // blocks across the various poly storage locations.
        // do 1 scatter/gather pass across the temp storage location to
        // distribute roots for all polynomials that hit the interval.
        /*
        
        // the elements that we store and the location of the bucket (via block num)
        // both depend on vroot, so we need to store the full roots.
        // the location of the bucket also depends on the polynomial number; we
        // add (poly_offset * BUCKET_ALLOC) to the sliceptr and poly_offset to the numptr
        // for each poly increment.  So we'd need to store the poly_offset as well, 
        // for each hit.  Or, we could store bnum = (root >> 15) + poly_offset 
        // and element = root & blocksizeM1?
        // BUCKET_ALLOC = 2048, so we'd store bnum = (root >> 15) + (poly_offset * BUCKET_ALLOC)?
        // look at how scatter work to see if we need to divide that last part by sizeof(uint32_t)
        // or something.
        // and that doesn't address the separate offset that numptr would need...
        
        num_prootstore1 = 0;
        num_prootstore2 = 0;
        num_nrootstore1 = 0;
        num_nrootstore2 = 0;
        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            if (gray[numB + p] > 0)
            {                
                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));
                mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_sub_epi32(vroot1, vpval);
                vroot2 = _mm512_sub_epi32(vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);                
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&prootstore1[num_prootstore1], mask1, vroot1);
                _mm512_mask_compressstoreu_epi32(&prootstore2[num_prootstore2], mask2, vroot2);
                num_prootstore1 += _mm_popcnt_u32(mask1);
                num_prootstore2 += _mm_popcnt_u32(mask2);
                
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&nrootstore1[num_nrootstore1], mask1, vnroot1);
                _mm512_mask_compressstoreu_epi32(&nrootstore2[num_nrootstore2], mask2, vnroot2);
                num_nrootstore1 += _mm_popcnt_u32(mask1);
                num_nrootstore2 += _mm_popcnt_u32(mask2);
            }
            else
            {
                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));
                vroot1 = _mm512_add_epi32(vroot1, vpval);
                vroot2 = _mm512_add_epi32(vroot2, vpval);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&prootstore1[num_prootstore1], mask1, vroot1);
                _mm512_mask_compressstoreu_epi32(&prootstore2[num_prootstore2], mask2, vroot2);
                num_prootstore1 += _mm_popcnt_u32(mask1);
                num_prootstore2 += _mm_popcnt_u32(mask2);
                
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&nrootstore1[num_nrootstore1], mask1, vnroot1);
                _mm512_mask_compressstoreu_epi32(&nrootstore2[num_nrootstore2], mask2, vnroot2);
                num_nrootstore1 += _mm_popcnt_u32(mask1);
                num_nrootstore2 += _mm_popcnt_u32(mask2);
            }
            
            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;
        }
        
        // now a single scatter/gather pass over the accumulated roots.  
        
        
        */
        
        
        
        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            //if (gray[numB + p] > 0)
            {         
                __mmask16 gm = gmask[p];

                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));


                //mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                //mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                //vroot1 = _mm512_sub_epi32(vroot1, vpval);
                //vroot2 = _mm512_sub_epi32(vroot2, vpval);
                //vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                //vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);     

                //if (gray[numB + p] > 0)
                mask1 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_mask_sub_epi32(vroot1, gm, vroot1, vpval);
                vroot2 = _mm512_mask_sub_epi32(vroot2, gm, vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, gm & mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, gm & mask2, vroot2, vprime);

                //else
                vroot1 = _mm512_mask_add_epi32(vroot1, ~gm, vroot1, vpval);
                vroot2 = _mm512_mask_add_epi32(vroot2, ~gm, vroot2, vpval);
                mask1 |= _mm512_mask_cmp_epu32_mask(~gm, vroot1, vprime, _MM_CMPINT_GE);
                mask2 |= _mm512_mask_cmp_epu32_mask(~gm, vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, ~gm & mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, ~gm & mask2, vroot2, vprime);


                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

#if 0
#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                    munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask1;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, velement1, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    _mm512_mask_i32scatter_epi32(numptr_p, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b1, vbnum1);
                    _mm512_store_epi32((__m512i*)e1, velement1);

                    while (mask1 > 0)
                    {
                        idx = _trail_zcnt(mask1);
                        bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                        *bptr = e1[idx];
                        numptr_p[b1[idx]]++;
                        mask1 = _reset_lsb(mask1);
                    }
                }

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                    munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask2;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_i32scatter_epi32(sliceptr_p, mask2, vaddr, velement2, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    _mm512_mask_i32scatter_epi32(numptr_p, munique & mask2, vbnum2, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b2, vbnum2);
                    _mm512_store_epi32((__m512i*)e2, velement2);

                    while (mask2 > 0)
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                        *bptr = e2[idx];
                        numptr_p[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }
#else
                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_P;
                }
#endif



                // and the -side roots; same story.
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);


#if 0
#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    // latency 26
                    vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                    // latency ??
                    munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                    // latency 9
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask1;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    // latency 17
                    _mm512_mask_i32scatter_epi32(sliceptr_n, mask1, vaddr, velement1, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    // latency 17
                    _mm512_mask_i32scatter_epi32(numptr_n, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b1, vbnum1);
                    _mm512_store_epi32((__m512i*)e1, velement1);

                    while (mask1 > 0)
                    {
                        idx = _trail_zcnt(mask1);
                        bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                        *bptr = e1[idx];
                        numptr_n[b1[idx]]++;
                        mask1 = _reset_lsb(mask1);
                    }
                }

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                    munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask2;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_i32scatter_epi32(sliceptr_n, mask2, vaddr, velement2, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    _mm512_mask_i32scatter_epi32(numptr_n, munique & mask2, vbnum2, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b2, vbnum2);
                    _mm512_store_epi32((__m512i*)e2, velement2);

                    while (mask2 > 0)
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                        *bptr = e2[idx];
                        numptr_n[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }
#else
                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_N;
                }
#endif


                
                
            }

#if 0
            else
            {
                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));
                vroot1 = _mm512_add_epi32(vroot1, vpval);
                vroot2 = _mm512_add_epi32(vroot2, vpval);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
         
#if 0
#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                    munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask1;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    _mm512_mask_i32scatter_epi32(sliceptr_p, mask1, vaddr, velement1, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    _mm512_mask_i32scatter_epi32(numptr_p, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b1, vbnum1);
                    _mm512_store_epi32((__m512i*)e1, velement1);

                    while (mask1 > 0)
                    {
                        idx = _trail_zcnt(mask1);
                        bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                        *bptr = e1[idx];
                        numptr_p[b1[idx]]++;
                        mask1 = _reset_lsb(mask1);
                    }
                }

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                    munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask2;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_i32scatter_epi32(sliceptr_p, mask2, vaddr, velement2, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    _mm512_mask_i32scatter_epi32(numptr_p, munique & mask2, vbnum2, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b2, vbnum2);
                    _mm512_store_epi32((__m512i*)e2, velement2);

                    while (mask2 > 0)
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                        *bptr = e2[idx];
                        numptr_p[b2[idx]]++;
                        mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
                    }
                }
#else
                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_P;
                }
#endif



                // and the -side roots; same story.
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);

#if 0
#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                    munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask1;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                    _mm512_mask_i32scatter_epi32(sliceptr_n, mask1, vaddr, velement1, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    _mm512_mask_i32scatter_epi32(numptr_n, munique & mask1, vbnum1, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b1, vbnum1);
                    _mm512_store_epi32((__m512i*)e1, velement1);

                    while (mask1 > 0)
                    {
                        idx = _trail_zcnt(mask1);
                        bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                        *bptr = e1[idx];
                        numptr_n[b1[idx]]++;
                        mask1 = _reset_lsb(mask1);
                    }
                }

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

#ifdef KNL_XLARGE_METHOD_1_CUTOFF
                if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
                {
                    vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                    munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                    vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
#endif

                    mNotProc = mask2;
                    while (mNotProc != 0)
                    {
                        mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                        vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                        mNotProc &= (~mconflictfree);
                        vnum = _mm512_mask_add_epi32(vnum, mNotProc, vnum, vone);
                    }
                    vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                    _mm512_mask_i32scatter_epi32(sliceptr_n, mask2, vaddr, velement2, _MM_SCALE_4);
                    vnum = _mm512_add_epi32(vnum, vone);
                    _mm512_mask_i32scatter_epi32(numptr_n, munique & mask2, vbnum2, vnum, _MM_SCALE_4);
                }
                else
#endif
                {
                    _mm512_store_epi32((__m512i*)b2, vbnum2);
                    _mm512_store_epi32((__m512i*)e2, velement2);

                    while (mask2 > 0)
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                        *bptr = e2[idx];
                        numptr_n[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }

#else
                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_N;
                }
#endif


                
            }
#endif

            // advance pointers
            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;

        }   

        _mm512_store_epi32((__m512i *)(&update_data.firstroots1[j]), vroot1);
        _mm512_store_epi32((__m512i *)(&update_data.firstroots2[j]), vroot2);

        // reset pointers
        sliceptr_p -= p * poly_offset * BUCKET_ALLOC;
        sliceptr_n -= p * poly_offset * BUCKET_ALLOC;
        numptr_p -= p * poly_offset;
        numptr_n -= p * poly_offset;

    }

    if (lp_bucket_p->list != NULL)
    {
        lp_bucket_p->num_slices = bound_index + 1;
        lp_bucket_p->logp[bound_index] = logp;
    }

    return;
}
#endif

//#endif

#ifdef USE_XLBATCHPOLY
#include <immintrin.h>

void nextBigRoots_32k_knl_polybatch(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    int* rootupdates = dconf->rootupdates;
    update_t update_data = dconf->update_data;

    uint32_t startprime = 2;
    uint32_t bound = sconf->factor_base->B;
    char* nu = dconf->curr_poly->nu;
    char* gray = dconf->curr_poly->gray;
    int numB = dconf->numB;
    uint32_t poly_offset = 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

    lp_bucket* lp_bucket_p = dconf->buckets;
    uint32_t med_B = sconf->factor_base->med_B;
    uint32_t large_B = sconf->factor_base->large_B;
    uint32_t xlB = sconf->factor_base->x2_large_B;

    uint32_t j, interval;
    int k, numblocks, idx;
    uint32_t root1, root2, nroot1, nroot2, prime;

    int bound_index = 0;
    int check_bound = BUCKET_ALLOC / 2 - 1;
    uint32_t bound_val = med_B;
    uint32_t* numptr_p, * numptr_n, * sliceptr_p, * sliceptr_n;

    uint32_t* bptr;
    int bnum, room;
    uint8_t logp = 0;
    int p;

    __m512i vshifted_index = _mm512_setr_epi32(
        0, 65536, 131072, 196608,
        262144, 327680, 393216, 458752,
        524288, 589824, 655360, 720896,
        786432, 851968, 917504, 983040);

    __m512i vnroot1, vnroot2, vzero = _mm512_setzero_epi32(), vone = _mm512_set1_epi32(1);
    __m512i vprime, vroot1, vroot2, vpval, vbnum1, vbnum2, vinterval, vindex, vnum, vaddr;
    __m512i vtmproot1, vtmproot2, vtmp;
    __m512i velement1, velement2, vbsize, vidx, vblockm1 = _mm512_set1_epi32(32767);
    __mmask16 mask1, mask2, mconflictfree, munique;
    __mmask16 mNotProc;

    uint32_t gmask[64];

    ALIGNED_MEM uint32_t e1[32];
    ALIGNED_MEM uint32_t e2[32];
    ALIGNED_MEM uint32_t b1[32];
    ALIGNED_MEM uint32_t b2[32];


#if 0
    // good for up to 16 polynomial batches...
    ALIGNED_MEM uint32_t prootstore1[256];
    ALIGNED_MEM uint32_t prootstore2[256];
    ALIGNED_MEM uint32_t nrootstore1[256];
    ALIGNED_MEM uint32_t nrootstore2[256];
    int num_prootstore1;
    int num_prootstore2;
    int num_nrootstore1;
    int num_nrootstore2;
#endif

    numblocks = sconf->num_blocks;
    interval = numblocks << 15;
    vinterval = _mm512_set1_epi32(interval);

    if (lp_bucket_p->alloc_slices != 0) // != NULL)
    {
        lp_bucket_p->fb_bounds[0] = med_B;

        sliceptr_p = lp_bucket_p->list;
        sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

        numptr_p = lp_bucket_p->num;
        numptr_n = lp_bucket_p->num + numblocks;

        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
            numptr_p[j] = 0;

        lp_bucket_p->num_slices = 0;

    }
    else
    {
        sliceptr_p = NULL;
        sliceptr_n = NULL;
        numptr_p = NULL;
        numptr_n = NULL;
    }

    for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
    {
        if (gray[numB + p] > 0)
        {
            gmask[p] = 0xffff;
        }
        else
        {
            gmask[p] = 0;
        }
    }

    bound_index = 0;
    bound_val = med_B;
    check_bound = med_B + BUCKET_ALLOC / 2;

    logp = update_data.logp[j - 1];
    for (j = xlB; j < bound; j += 16)
    {
        int p;

        CHECK_NEW_SLICE_BATCH(j);

        vprime = _mm512_load_epi32((__m512i*)(&update_data.prime[j]));
        vroot1 = _mm512_load_epi32((__m512i*)(&update_data.firstroots1[j]));
        vroot2 = _mm512_load_epi32((__m512i*)(&update_data.firstroots2[j]));

        // todo:
        // compute all roots for N polynomials.
        // compress-store those less than interval to temporary location.
        // adjust vbnum by the poly offset so scatter/gather grabs from
        // blocks across the various poly storage locations.
        // do 1 scatter/gather pass across the temp storage location to
        // distribute roots for all polynomials that hit the interval.
        /*

        // the elements that we store and the location of the bucket (via block num)
        // both depend on vroot, so we need to store the full roots.
        // the location of the bucket also depends on the polynomial number; we
        // add (poly_offset * BUCKET_ALLOC) to the sliceptr and poly_offset to the numptr
        // for each poly increment.  So we'd need to store the poly_offset as well,
        // for each hit.  Or, we could store bnum = (root >> 15) + poly_offset
        // and element = root & blocksizeM1?
        // BUCKET_ALLOC = 2048, so we'd store bnum = (root >> 15) + (poly_offset * BUCKET_ALLOC)?
        // look at how scatter work to see if we need to divide that last part by sizeof(uint32_t)
        // or something.
        // and that doesn't address the separate offset that numptr would need...

        num_prootstore1 = 0;
        num_prootstore2 = 0;
        num_nrootstore1 = 0;
        num_nrootstore2 = 0;
        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            if (gray[numB + p] > 0)
            {
                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));
                mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_sub_epi32(vroot1, vpval);
                vroot2 = _mm512_sub_epi32(vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&prootstore1[num_prootstore1], mask1, vroot1);
                _mm512_mask_compressstoreu_epi32(&prootstore2[num_prootstore2], mask2, vroot2);
                num_prootstore1 += _mm_popcnt_u32(mask1);
                num_prootstore2 += _mm_popcnt_u32(mask2);

                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&nrootstore1[num_nrootstore1], mask1, vnroot1);
                _mm512_mask_compressstoreu_epi32(&nrootstore2[num_nrootstore2], mask2, vnroot2);
                num_nrootstore1 += _mm_popcnt_u32(mask1);
                num_nrootstore2 += _mm_popcnt_u32(mask2);
            }
            else
            {
                vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));
                vroot1 = _mm512_add_epi32(vroot1, vpval);
                vroot2 = _mm512_add_epi32(vroot2, vpval);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&prootstore1[num_prootstore1], mask1, vroot1);
                _mm512_mask_compressstoreu_epi32(&prootstore2[num_prootstore2], mask2, vroot2);
                num_prootstore1 += _mm_popcnt_u32(mask1);
                num_prootstore2 += _mm_popcnt_u32(mask2);

                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
                _mm512_mask_compressstoreu_epi32(&nrootstore1[num_nrootstore1], mask1, vnroot1);
                _mm512_mask_compressstoreu_epi32(&nrootstore2[num_nrootstore2], mask2, vnroot2);
                num_nrootstore1 += _mm_popcnt_u32(mask1);
                num_nrootstore2 += _mm_popcnt_u32(mask2);
            }

            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;
        }

        // now a single scatter/gather pass over the accumulated roots.


        */

        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            //if (gray[numB + p] > 0)
            {
                __mmask16 gm = gmask[p];

                vpval = _mm512_load_epi32((__m512i*)(&rootupdates[(nu[numB + p] - 1) * bound + j]));


                //mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                //mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                //vroot1 = _mm512_sub_epi32(vroot1, vpval);
                //vroot2 = _mm512_sub_epi32(vroot2, vpval);
                //vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                //vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);     

                //if (gray[numB + p] > 0)
                mask1 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_mask_sub_epi32(vroot1, gm, vroot1, vpval);
                vroot2 = _mm512_mask_sub_epi32(vroot2, gm, vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, gm & mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, gm & mask2, vroot2, vprime);

                //else
                vroot1 = _mm512_mask_add_epi32(vroot1, ~gm, vroot1, vpval);
                vroot2 = _mm512_mask_add_epi32(vroot2, ~gm, vroot2, vpval);
                mask1 |= _mm512_mask_cmp_epu32_mask(~gm, vroot1, vprime, _MM_CMPINT_GE);
                mask2 |= _mm512_mask_cmp_epu32_mask(~gm, vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, ~gm & mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, ~gm & mask2, vroot2, vprime);


                vbnum1 = _mm512_srli_epi32(vroot1, 15);
                vbnum2 = _mm512_srli_epi32(vroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

#if 0

                {
                    _mm512_store_epi32((__m512i*)b1, vbnum1);
                    _mm512_store_epi32((__m512i*)e1, velement1);

                    while (mask1 > 0)
                    {
                        idx = _trail_zcnt(mask1);
                        bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                        *bptr = e1[idx];
                        numptr_p[b1[idx]]++;
                        mask1 = _reset_lsb(mask1);
                    }
                }

                {
                    _mm512_store_epi32((__m512i*)b2, vbnum2);
                    _mm512_store_epi32((__m512i*)e2, velement2);

                    while (mask2 > 0)
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                        *bptr = e2[idx];
                        numptr_p[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }
#else
                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_P;
                }
#endif



                // and the -side roots; same story.
                vnroot1 = _mm512_sub_epi32(vprime, vroot1);
                vnroot2 = _mm512_sub_epi32(vprime, vroot2);
                vbnum1 = _mm512_srli_epi32(vnroot1, 15);
                vbnum2 = _mm512_srli_epi32(vnroot2, 15);
                velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);


#if 0

                {
                    _mm512_store_epi32((__m512i*)b1, vbnum1);
                    _mm512_store_epi32((__m512i*)e1, velement1);

                    while (mask1 > 0)
                    {
                        idx = _trail_zcnt(mask1);
                        bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                        *bptr = e1[idx];
                        numptr_n[b1[idx]]++;
                        mask1 = _reset_lsb(mask1);
                    }
                }

                {
                    _mm512_store_epi32((__m512i*)b2, vbnum2);
                    _mm512_store_epi32((__m512i*)e2, velement2);

                    while (mask2 > 0)
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                        *bptr = e2[idx];
                        numptr_n[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }
#else
                {
                    idx = _mm_popcnt_u32(mask1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)b1, mask1, vbnum1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)e1, mask1, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(b1[idx])), mask2, vbnum2);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(&(e1[idx])), mask2, velement2);
                    idx += _mm_popcnt_u32(mask2);

                    SCATTER_COMPRESSED_VECTOR_N;
                }
#endif




            }

            // advance pointers
            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;

        }

        _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
        _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);

        // reset pointers
        sliceptr_p -= p * poly_offset * BUCKET_ALLOC;
        sliceptr_n -= p * poly_offset * BUCKET_ALLOC;
        numptr_p -= p * poly_offset;
        numptr_n -= p * poly_offset;

    }

    if (lp_bucket_p->list != NULL)
    {
        lp_bucket_p->num_slices = bound_index + 1;
        lp_bucket_p->logp[bound_index] = logp;
    }

    return;
}

#endif