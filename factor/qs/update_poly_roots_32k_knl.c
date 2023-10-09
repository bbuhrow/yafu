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
    uint8_t logp=0;

    __m512i vshifted_index = _mm512_setr_epi32(
        0, 65536, 131072, 196608,
        262144, 327680, 393216, 458752,
        524288, 589824, 655360, 720896,
        786432, 851968, 917504, 983040);

    __m512i vone = _mm512_set1_epi32(1);
    __m512i vtwo = _mm512_set1_epi32(2);
    __m512i vnroot1, vnroot2;
    __m512i vprime, vroot1, vroot2, vpval, vbnum1, vbnum2, vinterval;
    __m512i velement1, velement2, vblockm1 = _mm512_set1_epi32(32767);
    __mmask16 mask1, mask2;
    
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

#ifdef DEBUGPRINT_BATCHPOLY
        printf("begin bucket sieve for poly %d with bi = %u, bv = %u, cb = %u... ",
            pnum, bound_index, bound_val, med_B + BUCKET_ALLOC / 2);
#endif

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

            } while (1);

            vbnum1 = _mm512_srli_epi32(vnroot1, 15);
            vbnum2 = _mm512_srli_epi32(vnroot2, 15);
            velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot1, vblockm1));
            velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot2, vblockm1));

            do
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
        //logp = update_data.logp[large_B];

#if defined( USE_BATCHPOLY_X2 )
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

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.

            {
                //mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                //mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

                //_mm512_store_epi32((__m512i*)b1, vbnum1);
                //while (mask1 > 0)
                //{
                //    idx = _trail_zcnt(mask1);
                //    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                //    __mmask16 maska = _mm512_cmpeq_epu32_mask(vbnum1, _mm512_set1_epi32(b1[idx]));
                //    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, maska, velement1);
                //    numptr_p[b1[idx]] += _mm_popcnt_u32(maska);
                //    mask1 &= (~maska);
                //}
                //
                //_mm512_store_epi32((__m512i*)b1, vbnum2);
                //while (mask2 > 0)
                //{
                //    idx = _trail_zcnt(mask2);
                //    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                //    __mmask16 maska = _mm512_cmpeq_epu32_mask(vbnum2, _mm512_set1_epi32(b1[idx]));
                //    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, maska, velement2);
                //    numptr_p[b1[idx]] += _mm_popcnt_u32(maska);
                //    mask2 &= (~maska);
                //}

                //_mm512_store_epi32((__m512i*)b1, vbnum1);
                //_mm512_store_epi32((__m512i*)e1, velement1);
                //
                //while (mask1 > 0)
                //{
                //    idx = _trail_zcnt(mask1);
                //    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                //    *bptr = e1[idx];
                //    numptr_p[b1[idx]]++;
                //    mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
                //}
                //
                //_mm512_store_epi32((__m512i*)b1, vbnum2);
                //_mm512_store_epi32((__m512i*)e1, velement2);
                //
                //while (mask2 > 0)
                //{
                //    idx = _trail_zcnt(mask2);
                //    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                //    *bptr = e1[idx];
                //    numptr_p[b1[idx]]++;
                //    mask2 = _reset_lsb(mask2); //mask1 ^= (1 << idx);
                //}

                //__m512i vbmask = _mm512_xor_epi32(vpval, vpval);
                //for (k = 0; k < numblocks; k++)
                //{
                //    __mmask16 maska = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                //    __mmask16 maskb = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                //    bptr = sliceptr_p + (k << BUCKET_BITS) + numptr_p[k];
                //    idx = _mm_popcnt_u32(maska);
                //    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, maska, velement1);
                //    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), maskb, velement2);
                //    vbmask = _mm512_add_epi32(vbmask, vone);
                //    numptr_p[k] += (idx + _mm_popcnt_u32(maskb));
                //}

                __m512i vbmask = _mm512_xor_epi32(vbnum1, vbnum1);
                for (k = 0; k < numblocks; k++)
                {
                    __mmask16 maska = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    __mmask16 maskb = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_p + (k << BUCKET_BITS) + numptr_p[k];
                    idx = _mm_popcnt_u32(maska);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, maska, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), maskb, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_p[k] += (idx + _mm_popcnt_u32(maskb));
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
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
            }

            _mm512_store_epi32((__m512i*)b1, vbnum2);
            _mm512_store_epi32((__m512i*)e1, velement2);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while (mask2 > 0)
            {
                idx = _trail_zcnt(mask2);
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask2 = _reset_lsb(mask2); //mask1 ^= (1 << idx);
            }

#else

            {

                __m512i vbmask = _mm512_xor_epi32(vbnum1, vbnum1);
                for (k = 0; k < numblocks; k++)
                {
                    __mmask16 maska = _mm512_cmp_epu32_mask(vbnum1, vbmask, _MM_CMPINT_EQ);
                    __mmask16 maskb = _mm512_cmp_epu32_mask(vbnum2, vbmask, _MM_CMPINT_EQ);
                    bptr = sliceptr_n + (k << BUCKET_BITS) + numptr_n[k];
                    idx = _mm_popcnt_u32(maska);
                    _mm512_mask_compressstoreu_epi32((__m512i*)bptr, maska, velement1);
                    _mm512_mask_compressstoreu_epi32((__m512i*)(bptr + idx), maskb, velement2);
                    vbmask = _mm512_add_epi32(vbmask, vone);
                    numptr_n[k] += (idx + _mm_popcnt_u32(maskb));
                }
            }
#endif


        }


        // above this bound primes are much bigger and therefore
        // less likely to hit the interval.  For the ones that do,
        // there will be more blocks involved (because primes only
        // get this big for bigger inputs) so less likely to have
        // multiple hits per block.  For these reasons it is more
        // efficient to use the individual scatter rather than 
        // the block-based compress-store.
        
#if !defined(USE_BATCHPOLY_X2)

#if defined(USE_SS_SEARCH)

        for (j = xlarge_B; j < sconf->factor_base->ss_start_B; j += 16, ptr += 16)
#else
        for (j = xlarge_B; j < bound; j += 16, ptr += 16)
#endif
        {
            //int i;
            int k;

            CHECK_NEW_SLICE(j);

            //if (j == xlarge_B)
            //{
            //    printf("v = %d, sign = %d, update = %d\n", v, sign, *ptr);
            //}

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

            //if (j == xlarge_B)
            //{
            //    char sym1, sym2;
            //
            //    sym1 = ' ';
            //    sym2 = ' ';
            //
            //    if (update_data.firstroots1[j] < interval)
            //        sym1 = '*';
            //    if (update_data.firstroots2[j] < interval)
            //        sym2 = '*';
            //
            //    printf("new roots = %d%c,%d%c\n",
            //        update_data.firstroots1[j], sym1, update_data.firstroots2[j], sym2);
            //}

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

            __m512i vebase = velement1;
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot1, vblockm1));
            velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vroot2, vblockm1));

            do
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

            } while (1);

            vbnum1 = _mm512_srli_epi32(vnroot1, 15);
            vbnum2 = _mm512_srli_epi32(vnroot2, 15);
            velement1 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot1, vblockm1));
            velement2 = _mm512_or_epi32(vebase, _mm512_and_epi32(vnroot2, vblockm1));

            do
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
        //logp = update_data.logp[large_B];

#if defined( USE_BATCHPOLY_X2 )
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
                __m512i vbmask = _mm512_xor_epi32(vbnum1, vbnum1);
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
                __m512i vbmask = _mm512_xor_epi32(vbnum1, vbnum1);
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

#if !defined(USE_BATCHPOLY_X2)

#if defined(USE_SS_SEARCH)
        for (j = xlarge_B; j < sconf->factor_base->ss_start_B; j += 16, ptr += 16)
#else
        for (j = xlarge_B; j < bound; j += 16, ptr += 16)
#endif
        {
            //int i;

            CHECK_NEW_SLICE(j);

            //if (j == xlarge_B)
            //{
            //    printf("v = %d, sign = %d, update = %d\n", v, sign, *ptr);
            //}

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

            //if (j == xlarge_B)
            //{
            //    char sym1, sym2;
            //
            //    sym1 = ' ';
            //    sym2 = ' ';
            //
            //    if (update_data.firstroots1[j] < interval)
            //        sym1 = '*';
            //    if (update_data.firstroots2[j] < interval)
            //        sym2 = '*';
            //
            //    printf("new roots = %d%c,%d%c\n",
            //        update_data.firstroots1[j], sym1, update_data.firstroots2[j], sym2);
            //}

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


    // every N iterations we do the bucket sieve for extra large primes, which means we
    // need to slightly change how buckets operate.  
    // When this routine runs, the buckets associated with the 
    // current poly will be filled with primes up to fb->large_B while the buckets
    // associated with the rest of the polys in this batch will be new (empty).  So bound_index,
    // (which tracks the current slice), bound_val (which tracks the prime index that
    // starts each slice), and check_bound (which control when we look at the buckets
    // to potentially start a new slice), will need to be different for each poly.  

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
                // a mask for the sign of the Gray code (controls whether we
                // need to do a modadd or modsub for each of the polys).
                // (computing both using this mask is ~20% faster than using
                // a branch.)
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
                //{
                    mask1 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot1, _MM_CMPINT_GT);
                    mask2 = _mm512_mask_cmp_epu32_mask(gm, vpval, vroot2, _MM_CMPINT_GT);
                    vroot1 = _mm512_mask_sub_epi32(vroot1, gm, vroot1, vpval);
                    vroot2 = _mm512_mask_sub_epi32(vroot2, gm, vroot2, vpval);
                    vroot1 = _mm512_mask_add_epi32(vroot1, gm & mask1, vroot1, vprime);
                    vroot2 = _mm512_mask_add_epi32(vroot2, gm & mask2, vroot2, vprime);
                //}
                //else
                //{
                    vroot1 = _mm512_mask_add_epi32(vroot1, ~gm, vroot1, vpval);
                    vroot2 = _mm512_mask_add_epi32(vroot2, ~gm, vroot2, vpval);
                    mask1 |= _mm512_mask_cmp_epu32_mask(~gm, vroot1, vprime, _MM_CMPINT_GE);
                    mask2 |= _mm512_mask_cmp_epu32_mask(~gm, vroot2, vprime, _MM_CMPINT_GE);
                    vroot1 = _mm512_mask_sub_epi32(vroot1, ~gm & mask1, vroot1, vprime);
                    vroot2 = _mm512_mask_sub_epi32(vroot2, ~gm & mask2, vroot2, vprime);
                //}

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
            }

            // now we can store the updated roots for these primes
            _mm512_store_epi32((__m512i*)(&update_data.firstroots1[j]), vroot1);
            _mm512_store_epi32((__m512i*)(&update_data.firstroots2[j]), vroot2);
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

#ifdef DEBUGPRINT_BATCHPOLY
    printf("complete.\n"); fflush(stdout);
#endif
    return;
}



#endif

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


#if defined(USE_BATCHPOLY) && defined(USE_AVX512F)
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
