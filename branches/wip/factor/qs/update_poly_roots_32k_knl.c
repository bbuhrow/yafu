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
#include "util.h"
#include "common.h"
#include "poly_macros_32k.h"
#include "poly_macros_common.h"

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
#define KNL_XLARGE_METHOD_1_CUTOFF 7
//#define XLARGE_BUCKET_DEBUG

void nextRoots_32k_knl_bucket(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //update the roots 
    int *rootupdates = dconf->rootupdates;

    update_t update_data = dconf->update_data;

    uint32 startprime = 2;
    uint32 bound = sconf->factor_base->B;

    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    int *ptr;

    lp_bucket *lp_bucket_p = dconf->buckets;
    uint32 med_B = sconf->factor_base->med_B;
    uint32 large_B = sconf->factor_base->large_B;

    uint32 j, interval;
    int k,numblocks,idx;
    uint32 root1, root2, prime;

    int bound_index=0;
    int check_bound = BUCKET_ALLOC/2 - 1;
    uint32 bound_val = med_B;
    uint32 *numptr_p, *numptr_n, *sliceptr_p,*sliceptr_n;

    uint32 *bptr;
    int bnum, room;

    uint8 logp=0;
    polysieve_t helperstruct;

    __m512i vshifted_index = _mm512_setr_epi32(
        0, 65536, 131072, 196608, 
        262144,327680, 393216, 458752, 
        524288,589824, 655360, 720896, 
        786432,851968, 917504, 983040);

    __m512i vnroot1, vnroot2, vzero = _mm512_setzero_epi32(), vone = _mm512_set1_epi32(1);
    __m512i vprime, vroot1, vroot2, vpval, vbnum1, vbnum2, vinterval, vindex, vnum, vaddr;
    __m512i velement1, velement2, vbsize, vidx, vblockm1 = _mm512_set1_epi32(32767);
    __mmask16 mask1, mask2, mconflictfree, munique;
    __mmask16 mNotProc;
    
    ALIGNED_MEM uint32 e1[16];
    ALIGNED_MEM uint32 e2[16];
    ALIGNED_MEM uint32 b1[16];
    ALIGNED_MEM uint32 b2[16];

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

    k=0;
    ptr = &rootupdates[(v-1) * bound + med_B];	

    if (sign > 0)
    {        

        bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC/2;

		if (VFLAG > 2)
		{
			printf("commencing 1\n");
		}

        logp = update_data.logp[med_B];
        for (j=med_B;j<large_B;j+=16,ptr+=16)				
        {			
            int i;

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
            

            // do the first 16
#ifdef noUSE_AVX512PF
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
            vroot1 = _mm512_add_epi32(vroot1, vprime);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
#else
            mask1 = 0xffff;
#endif
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i *)e1, 
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
            
            // do the first 16
#ifdef noUSE_AVX512PF
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
            vroot2 = _mm512_add_epi32(vroot2, vprime);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
#else
            mask2 = 0xffff;
#endif
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i *)e1, 
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
            
            // do the first 16
#ifdef noUSE_AVX512PF
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
#else
            mask1 = 0xffff;
#endif
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
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            }
            
            // do the first 16
#ifdef noUSE_AVX512PF
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
#else
            mask2 = 0xffff;
#endif
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
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            }
                        
        }

		if (VFLAG > 2)
		{
			printf("commencing 2\n");
		}

        logp = update_data.logp[j - 1];
        for (j = large_B; j<bound; j += 16, ptr += 16)
        {
            int i;
            int k;

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
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);



            /*
            // https://gcc.gnu.org/wiki/cauldron2014?action=AttachFile&do=get&target=Cauldron14_AVX-512_Vector_ISA_Kirill_Yukhin_20140711.pdf
            // conflict detect logic.  
            // we do a comparison on the result of vconflictd to make a mask.  also, the
            // index must be updated to remove processed lanes every iteration so that
            // they don't participate in subsequent conflict instructions.  otherwise
            // the loop becomes infinite (will always find the same set of conflicts).
            index = vload &B[i] // Load 16 B[i]
            pending_elem = 0xFFFF; // all still remaining
            do {
              curr_elem = get_conflict_free_subset(index, pending_elem)
              old_val = vgather {curr_elem} A, index // Grab A[B[i]]
              new_val = vadd old_val, +1.0 // Compute new values
              vscatter A {curr_elem}, index, new_val // Update A[B[i]]
              pending_elem = pending_elem ^ curr_elem // remove done idx
            } while (pending_elem)
              
            // here is how to do it faster, using the vector info returned by the vconflictd instruction.
            // each lane will hold a bitmask of lanes ahead of it that are conflicted.  lanes
            // with 0 are unconflicted.  
            // process lanes with a 0 (cmp to make a mask) that haven't already been processed.
            // remove the first bit of all non-zero lanes as these have now been processed.
            // repeat
            // will only have to do vconflict once and vgather/vscatter of block nums once.
            // 
            */

#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
            {                
				
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
                #ifdef XLARGE_BUCKET_DEBUG
                    printf("begin CFGS loop1: unique mask = %08x\n", munique);
                    k = 0;
                #endif

                while (mNotProc != 0)   
                {     
                  
                     // set a mask of non-conflicted indices that are also interval hits
                     mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                     
                     #ifdef XLARGE_BUCKET_DEBUG
              _mm512_store_epi32((__m512i *)b1, vbnum1);
              printf("CF mask = %08x for processed %08x and blocknums: ", mconflictfree, mprocessed);
              for (i = 15; i >= 0; i--)
                printf("%u ", b1[i]);
              printf("\n");
              printf("vconflict: \n");
              _mm512_store_epi32((__m512i *)e1, vindex);
              for (i = 15; i >= 0; i--)
                printf("%u ", e1[i]);
              printf("\n");
              
              printf("vnum: \n");
              _mm512_store_epi32((__m512i *)e1, vnum);
              for (i = 15; i >= 0; i--)
                printf("%u ", e1[i]);
              printf("\n");
              
              k++;
              if (k > 16)
                break;
                     #endif
			         
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
                
                #ifdef XLARGE_BUCKET_DEBUG
                printf("final vnum: \n");
                _mm512_store_epi32((__m512i *)e1, vnum);
                for (i = 15; i >= 0; i--)
                    printf("%u ", e1[i]);
                printf("\n");
                exit(1);
                #endif
            
            }
            else
            {
              
                 _mm512_store_epi32((__m512i *)b1, vbnum1);
                 _mm512_store_epi32((__m512i *)e1, velement1);
			     
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

            }



#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
            {

                vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_p, _MM_SCALE_4);
			    
#ifdef KNL_SCATTER_PREFETCH
				vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
				_mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
				//_mm512_mask_prefetch_i32scatter_ps(numptr_p, munique & mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
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
			{
				_mm512_store_epi32((__m512i *)b2, vbnum2);
				_mm512_store_epi32((__m512i *)e2, velement2);

				while (mask2 > 0)  
				{
				    idx = _trail_zcnt(mask2);
				    bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
				    *bptr = e2[idx];
				    numptr_p[b2[idx]]++;
				    mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
				}
			}
            
            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);



#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
            {
				vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
				munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
				vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
            vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
            _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
            //_mm512_mask_prefetch_i32scatter_ps(numptr_n, munique & mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
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
            {
				_mm512_store_epi32((__m512i *)b1, vbnum1);
				_mm512_store_epi32((__m512i *)e1, velement1);
				
				while (mask1 > 0)  
				{
				    idx = _trail_zcnt(mask1);
				    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
				    *bptr = e1[idx];
				    numptr_n[b1[idx]]++;
				    mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
				}
            }
            
#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
            {
				vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
				munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
				vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
            vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
            _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
            //_mm512_mask_prefetch_i32scatter_ps(numptr_n, munique & mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
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
            {
				_mm512_store_epi32((__m512i *)b2, vbnum2);
				_mm512_store_epi32((__m512i *)e2, velement2);
				
				while (mask2 > 0)  
				{
				    idx = _trail_zcnt(mask2);
				    bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
				    *bptr = e2[idx];
				    numptr_n[b2[idx]]++;
				    mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
				}
            }

        }
        
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
        bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC / 2;

		if (VFLAG > 2)
		{
			printf("commencing 1n\n");
		}

        logp = update_data.logp[med_B];
        for (j = med_B; j<large_B; j += 16, ptr += 16)
        {
            int i;

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


            // do the first 16
#ifdef noUSE_AVX512PF
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
            //_mm512_mask_prefetch_i32scatter_ps(numptr_p, munique & mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
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
            vroot1 = _mm512_add_epi32(vroot1, vprime);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
#else
            mask1 = 0xffff;
#endif
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i *)e1,
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

                while (m > 0)  
                {
                    idx = _trail_zcnt(m);
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            }

            // do the first 16
#ifdef noUSE_AVX512PF
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
            //_mm512_mask_prefetch_i32scatter_ps(numptr_p, munique & mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
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
            vroot2 = _mm512_add_epi32(vroot2, vprime);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
#else
            mask2 = 0xffff;
#endif
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i *)e1,
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

            // do the first 16
#ifdef noUSE_AVX512PF
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
            //_mm512_mask_prefetch_i32scatter_ps(numptr_n, munique & mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
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
#else
            mask1 = 0xffff;
#endif
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
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            }

            // do the first 16
#ifdef noUSE_AVX512PF
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
            // _mm512_mask_prefetch_i32scatter_ps(numptr_n, munique & mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
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
#else
            mask2 = 0xffff;
#endif
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
                    m = _reset_lsb(m); //m ^= (1 << idx);
                }

                vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            }
            
 
        }

		if (VFLAG > 2)
		{
			printf("commencing 2n\n");
		}

        logp = update_data.logp[j - 1];
        for (j = large_B; j<bound; j += 16, ptr += 16)
        {
            int i;

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
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);


    
#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
            {
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
                //_mm512_mask_prefetch_i32scatter_ps(numptr_p, munique & mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

                mNotProc = mask1;
                while (mNotProc != 0)   
                {     
                  mconflictfree = _mm512_cmp_epu32_mask(vindex, vzero, _MM_CMPINT_EQ);
                  vindex = _mm512_mask_and_epi32(vindex, mNotProc, _mm512_set1_epi32(~mconflictfree), vindex);
                  mNotProc &= (~mconflictfree);
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

                _mm512_store_epi32((__m512i *)b1, vbnum1);
                _mm512_store_epi32((__m512i *)e1, velement1);

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

            }


#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
            {
                vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_p, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_p, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
                //_mm512_mask_prefetch_i32scatter_ps(numptr_p, munique & mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
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
            {
                _mm512_store_epi32((__m512i *)b2, vbnum2);
                _mm512_store_epi32((__m512i *)e2, velement2);

                while (mask2 > 0)  
                {
                    idx = _trail_zcnt(mask2);
                    bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                    *bptr = e2[idx];
                    numptr_p[b2[idx]]++;
                    mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
                }
            }

            // and the -side roots; same story.
            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);



#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask1) > KNL_XLARGE_METHOD_1_CUTOFF)
            {
                vindex = _mm512_mask_conflict_epi32(vzero, mask1, vbnum1);
                munique = ~_mm512_mask_reduce_or_epi32(mask1, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask1, vbnum1, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum1, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask1, vaddr, 4, KNL_SCATTER_PREFETCH);
                //_mm512_mask_prefetch_i32scatter_ps(numptr_n, munique & mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
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
            {
                _mm512_store_epi32((__m512i *)b1, vbnum1);
                _mm512_store_epi32((__m512i *)e1, velement1);

                while (mask1 > 0)  
                {
                    idx = _trail_zcnt(mask1);
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    mask1 = _reset_lsb(mask1); //mask1 ^= (1 << idx);
                }
            }

#ifdef KNL_SCATTER_PREFETCH_NUM
            _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
#endif

            if (_mm_popcnt_u32(mask2) > KNL_XLARGE_METHOD_1_CUTOFF)
            {
                vindex = _mm512_mask_conflict_epi32(vzero, mask2, vbnum2);
                munique = ~_mm512_mask_reduce_or_epi32(mask2, vindex);
                vnum = _mm512_mask_i32gather_epi32(vzero, mask2, vbnum2, numptr_n, _MM_SCALE_4);

#ifdef KNL_SCATTER_PREFETCH
                vaddr = _mm512_add_epi32(_mm512_slli_epi32(vbnum2, BUCKET_BITS), vnum);
                _mm512_mask_prefetch_i32scatter_ps(sliceptr_n, mask2, vaddr, 4, KNL_SCATTER_PREFETCH);
                //_mm512_mask_prefetch_i32scatter_ps(numptr_n, munique & mask2, vbnum2, 4, KNL_SCATTER_PREFETCH);
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
            {
                _mm512_store_epi32((__m512i *)b2, vbnum2);
                _mm512_store_epi32((__m512i *)e2, velement2);

                while (mask2 > 0)  
                {
                    idx = _trail_zcnt(mask2);
                    bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                    *bptr = e2[idx];
                    numptr_n[b2[idx]]++;
                    mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
                }
            }

        }

    }

    if (lp_bucket_p->list != NULL)
    {
        lp_bucket_p->num_slices = bound_index + 1;
        lp_bucket_p->logp[bound_index] = logp;
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



void nextRoots_32k_knl_polybatch(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    int *rootupdates = dconf->rootupdates;
    update_t update_data = dconf->update_data;

    uint32 startprime = 2;
    uint32 bound = sconf->factor_base->B;
    char *nu = dconf->curr_poly->nu;
    char *gray = dconf->curr_poly->gray;
    int numB = dconf->numB;
    uint32 poly_offset = 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

    lp_bucket *lp_bucket_p = dconf->buckets;
    uint32 med_B = sconf->factor_base->med_B;
    uint32 large_B = sconf->factor_base->large_B;

    uint32 j, interval;
    int k, numblocks, idx;
    uint32 root1, root2, nroot1, nroot2, prime;

    int bound_index = 0;
    int check_bound = BUCKET_ALLOC / 2 - 1;
    uint32 bound_val = med_B;
    uint32 *numptr_p, *numptr_n, *sliceptr_p, *sliceptr_n;

    uint32 *bptr;
    int bnum, room;
    uint8 logp = 0;

    __m512i vshifted_index = _mm512_setr_epi32(
        0, 65536, 131072, 196608, 
        262144,327680, 393216, 458752, 
        524288,589824, 655360, 720896, 
        786432,851968, 917504, 983040);

    __m512i vnroot1, vnroot2, vzero = _mm512_setzero_epi32(), vone = _mm512_set1_epi32(1);
    __m512i vprime, vroot1, vroot2, vpval, vbnum1, vbnum2, vinterval, vindex, vnum, vaddr;
    __m512i vtmproot1, vtmproot2;
    __m512i velement1, velement2, vbsize, vidx, vblockm1 = _mm512_set1_epi32(32767);
    __mmask16 mask1, mask2, mconflictfree, munique;
    __mmask16 mNotProc;
    
    ALIGNED_MEM uint32 e1[16];
    ALIGNED_MEM uint32 e2[16];
    ALIGNED_MEM uint32 b1[16];
    ALIGNED_MEM uint32 b2[16];
    
    // good for up to 16 polynomial batches...
    ALIGNED_MEM uint32 prootstore1[256];
    ALIGNED_MEM uint32 prootstore2[256];
    ALIGNED_MEM uint32 nrootstore1[256];
    ALIGNED_MEM uint32 nrootstore2[256];
    int num_prootstore1;
    int num_prootstore2;
    int num_nrootstore1;
    int num_nrootstore2;

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
   
    bound_index = 0;
    bound_val = med_B;
    check_bound = med_B + BUCKET_ALLOC / 2;

    logp = update_data.logp[med_B];
    for (j = med_B; j < large_B; j += 16)
    {
        int p;
        
        CHECK_NEW_SLICE_BATCH(j);        
        
        vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
        vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
        vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));
        
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
            
                
                
            }

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
        // look at how scatter work to see if we need to divide that last part by sizeof(uint32)
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
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

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
                {
                    _mm512_store_epi32((__m512i *)b1, vbnum1);
                    _mm512_store_epi32((__m512i *)e1, velement1);

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
                {
                    _mm512_store_epi32((__m512i *)b2, vbnum2);
                    _mm512_store_epi32((__m512i *)e2, velement2);

                    while (mask2 > 0)  
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                        *bptr = e2[idx];
                        numptr_p[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }

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

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

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
                {
                    _mm512_store_epi32((__m512i *)b1, vbnum1);
                    _mm512_store_epi32((__m512i *)e1, velement1);

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
                {
                    _mm512_store_epi32((__m512i *)b2, vbnum2);
                    _mm512_store_epi32((__m512i *)e2, velement2);

                    while (mask2 > 0)  
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                        *bptr = e2[idx];
                        numptr_n[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }
                
                
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
                velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
                velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
         
                
#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_p, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

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
                {
                    _mm512_store_epi32((__m512i *)b1, vbnum1);
                    _mm512_store_epi32((__m512i *)e1, velement1);

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
                {
                    _mm512_store_epi32((__m512i *)b2, vbnum2);
                    _mm512_store_epi32((__m512i *)e2, velement2);

                    while (mask2 > 0)  
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                        *bptr = e2[idx];
                        numptr_p[b2[idx]]++;
                        mask2 = _reset_lsb(mask2); //mask2 ^= (1 << idx);
                    }
                }

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

#ifdef KNL_SCATTER_PREFETCH_NUM
                _mm512_mask_prefetch_i32scatter_ps(numptr_n, mask1, vbnum1, 4, KNL_SCATTER_PREFETCH);
#endif

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
                {
                    _mm512_store_epi32((__m512i *)b1, vbnum1);
                    _mm512_store_epi32((__m512i *)e1, velement1);

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
                {
                    _mm512_store_epi32((__m512i *)b2, vbnum2);
                    _mm512_store_epi32((__m512i *)e2, velement2);

                    while (mask2 > 0)  
                    {
                        idx = _trail_zcnt(mask2);
                        bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                        *bptr = e2[idx];
                        numptr_n[b2[idx]]++;
                        mask2 = _reset_lsb(mask2);
                    }
                }
                
                
            }

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
