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

#ifdef TARGET_KNC

#include <immintrin.h>

void nextRoots_32k_knc_small(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //update the roots 
    sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
    sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
    int *rootupdates = dconf->rootupdates;

    update_t update_data = dconf->update_data;

    uint32 startprime = 2;
    uint32 bound = sconf->factor_base->B;

    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    int *ptr;
    uint32 med_B = sconf->factor_base->med_B;

    uint32 j;
    int k;
    uint32 root1, root2, prime;
    uint8 logp = 0;

    k = 0;
    ptr = &rootupdates[(v - 1) * bound + startprime];

    if (sign > 0)
    {
        for (j = startprime; j<sconf->sieve_small_fb_start; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            //we don't sieve these, so ordering doesn't matter
            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            fb_p->root1[j] = (uint16)root1;
            fb_p->root2[j] = (uint16)root2;
            fb_n->root1[j] = (uint16)(prime - root2);
            fb_n->root2[j] = (uint16)(prime - root1);
            if (fb_n->root1[j] == prime)
                fb_n->root1[j] = 0;
            if (fb_n->root2[j] == prime)
                fb_n->root2[j] = 0;

        }

        // do one at a time up to the 10bit boundary, where
        // we can start doing things 8 at a time and be
        // sure we can use aligned moves (static_data_init).		
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint32)update_data.sm_firstroots1[j];
            root2 = (uint32)update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16)root2;
                update_data.sm_firstroots2[j] = (uint16)root1;

                fb_p->root1[j] = (uint16)root2;
                fb_p->root2[j] = (uint16)root1;
                fb_n->root1[j] = (uint16)(prime - root1);
                fb_n->root2[j] = (uint16)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16)root1;
                update_data.sm_firstroots2[j] = (uint16)root2;

                fb_p->root1[j] = (uint16)root1;
                fb_p->root2[j] = (uint16)root2;
                fb_n->root1[j] = (uint16)(prime - root2);
                fb_n->root2[j] = (uint16)(prime - root1);
            }
        }

        ptr = &dconf->rootupdates[(v - 1) * bound + sconf->factor_base->fb_10bit_B];
        for (j = sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
        {
            if ((j & 15) == 0)
                break;

            prime = update_data.prime[j];
            root1 = update_data.sm_firstroots1[j];
            root2 = update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16)root2;
                update_data.sm_firstroots2[j] = (uint16)root1;

                fb_p->root1[j] = (uint16)root2;
                fb_p->root2[j] = (uint16)root1;
                fb_n->root1[j] = (uint16)(prime - root1);
                fb_n->root2[j] = (uint16)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16)root1;
                update_data.sm_firstroots2[j] = (uint16)root2;

                fb_p->root1[j] = (uint16)root1;
                fb_p->root2[j] = (uint16)root2;
                fb_n->root1[j] = (uint16)(prime - root2);
                fb_n->root2[j] = (uint16)(prime - root1);
            }
        }

        for (; j < sconf->factor_base->med_B; j += 16, ptr += 16)
        {
            __m512i vprime, vroot1, vroot2, vpval;
            __m512i vmax, vmin;
            __mmask16 mask1, mask2;

            vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
            vroot1 = _mm512_extload_epi32((__m512i *)(&update_data.sm_firstroots1[j]),
                _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            vroot2 = _mm512_extload_epi32((__m512i *)(&update_data.sm_firstroots2[j]),
                _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            vpval = _mm512_load_epi32((__m512i *)ptr);
            mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
            mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
            vroot1 = _mm512_sub_epi32(vroot1, vpval);
            vroot2 = _mm512_sub_epi32(vroot2, vpval);
            vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
            vmax = _mm512_max_epu32(vroot1, vroot2);
            vmin = _mm512_min_epu32(vroot1, vroot2);

            _mm512_extstore_epi32((__m512i *)(&update_data.sm_firstroots1[j]), vmin,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
            _mm512_extstore_epi32((__m512i *)(&update_data.sm_firstroots2[j]), vmax,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

            _mm512_extstore_epi32((__m512i *)(&fb_p->root1[j]), vmin,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
            _mm512_extstore_epi32((__m512i *)(&fb_p->root2[j]), vmax,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

            _mm512_extstore_epi32((__m512i *)(&fb_n->root1[j]), _mm512_sub_epi32(vprime, vmax),
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
            _mm512_extstore_epi32((__m512i *)(&fb_n->root2[j]), _mm512_sub_epi32(vprime, vmin),
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

        }

    }
    else
    {

        for (j = startprime; j<sconf->sieve_small_fb_start; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            //we don't sieve these, so ordering doesn't matter
            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            fb_p->root1[j] = (uint16)root1;
            fb_p->root2[j] = (uint16)root2;
            fb_n->root1[j] = (uint16)(prime - root2);
            fb_n->root2[j] = (uint16)(prime - root1);
            if (fb_n->root1[j] == prime)
                fb_n->root1[j] = 0;
            if (fb_n->root2[j] == prime)
                fb_n->root2[j] = 0;

        }

        // do one at a time up to the 10bit boundary, where
        // we can start doing things 8 at a time and be
        // sure we can use aligned moves (static_data_init).	
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint32)update_data.sm_firstroots1[j];
            root2 = (uint32)update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16)root2;
                update_data.sm_firstroots2[j] = (uint16)root1;

                fb_p->root1[j] = (uint16)root2;
                fb_p->root2[j] = (uint16)root1;
                fb_n->root1[j] = (uint16)(prime - root1);
                fb_n->root2[j] = (uint16)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16)root1;
                update_data.sm_firstroots2[j] = (uint16)root2;

                fb_p->root1[j] = (uint16)root1;
                fb_p->root2[j] = (uint16)root2;
                fb_n->root1[j] = (uint16)(prime - root2);
                fb_n->root2[j] = (uint16)(prime - root1);
            }
        }


        ptr = &dconf->rootupdates[(v - 1) * bound + sconf->factor_base->fb_10bit_B];
        for (j = sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
        {
            if ((j % 16) == 0)
                break;

            prime = update_data.prime[j];
            root1 = update_data.sm_firstroots1[j];
            root2 = update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16)root2;
                update_data.sm_firstroots2[j] = (uint16)root1;

                fb_p->root1[j] = (uint16)root2;
                fb_p->root2[j] = (uint16)root1;
                fb_n->root1[j] = (uint16)(prime - root1);
                fb_n->root2[j] = (uint16)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16)root1;
                update_data.sm_firstroots2[j] = (uint16)root2;

                fb_p->root1[j] = (uint16)root1;
                fb_p->root2[j] = (uint16)root2;
                fb_n->root1[j] = (uint16)(prime - root2);
                fb_n->root2[j] = (uint16)(prime - root1);
            }
        }

        for (; j < sconf->factor_base->med_B; j += 16, ptr += 16)
        {
            int i;
            __m512i vprime, vroot1, vroot2, vpval;
            __m512i vmax, vmin;
            __mmask16 mask1, mask2;

            vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
            vroot1 = _mm512_extload_epi32((__m512i *)(&update_data.sm_firstroots1[j]),
                _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            vroot2 = _mm512_extload_epi32((__m512i *)(&update_data.sm_firstroots2[j]),
                _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);
            vpval = _mm512_load_epi32((__m512i *)ptr);
            vroot1 = _mm512_add_epi32(vroot1, vpval);
            vroot2 = _mm512_add_epi32(vroot2, vpval);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
            vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            vmax = _mm512_max_epu32(vroot1, vroot2);
            vmin = _mm512_min_epu32(vroot1, vroot2);

            _mm512_extstore_epi32((__m512i *)(&update_data.sm_firstroots1[j]), vmin,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
            _mm512_extstore_epi32((__m512i *)(&update_data.sm_firstroots2[j]), vmax,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

            _mm512_extstore_epi32((__m512i *)(&fb_p->root1[j]), vmin,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
            _mm512_extstore_epi32((__m512i *)(&fb_p->root2[j]), vmax,
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

            _mm512_extstore_epi32((__m512i *)(&fb_n->root1[j]), _mm512_sub_epi32(vprime, vmax),
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);
            _mm512_extstore_epi32((__m512i *)(&fb_n->root2[j]), _mm512_sub_epi32(vprime, vmin),
                _MM_DOWNCONV_EPI32_UINT16, _MM_HINT_NONE);

        }

    }


    return;
}


void nextRoots_32k_knc_bucket(static_conf_t *sconf, dynamic_conf_t *dconf)
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

    __m512i vnroot1, vnroot2;
    __m512i vprime, vroot1, vroot2, vpval, vbnum1, vbnum2, vinterval;
    __m512i velement1, velement2, vbsize, vidx, vblockm1 = _mm512_set1_epi32(32767);
    __mmask16 mask1, mask2;
    __declspec(aligned(64)) uint32 e1[16];
    __declspec(aligned(64)) uint32 e2[16];
    __declspec(aligned(64)) uint32 b1[16];
    __declspec(aligned(64)) uint32 b2[16];

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

            // the first set of roots is guaranteed to have at least one
            // hit in the interval.
            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot1, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
            }

            // after that, roots eventually start dropping out of the interval
            // as they are iterated, in no particular order (because the
            // beginning roots are essentially random).
            vroot1 = _mm512_add_epi32(vroot1, vprime);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            }

            // ditto for root2.
            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot2, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
            }

            vroot2 = _mm512_add_epi32(vroot2, vprime);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            }
            
            // and now the negative side roots that we computed earlier.
            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot1, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
            }

            vnroot1 = _mm512_add_epi32(vnroot1, vprime);
            mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot1, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            }

            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot2, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
            }

            vnroot2 = _mm512_add_epi32(vnroot2, vprime);
            mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot2, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            }
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


            _mm512_store_epi32((__m512i *)b1, vbnum1);
            _mm512_store_epi32((__m512i *)b2, vbnum2);
            _mm512_store_epi32((__m512i *)e1, velement1);
            _mm512_store_epi32((__m512i *)e2, velement2);

            // extra big roots are much easier because they hit at most once
            // in the entire +side interval.  no need to iterate.
            while ((idx = _mm_tzcnt_32(mask1)) < 16)
            {
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask1 ^= (1 << idx);
            }

            while ((idx = _mm_tzcnt_32(mask2)) < 16)
            {
                bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                *bptr = e2[idx];
                numptr_p[b2[idx]]++;
                mask2 ^= (1 << idx);
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

            _mm512_store_epi32((__m512i *)b1, vbnum1);
            _mm512_store_epi32((__m512i *)b2, vbnum2);
            _mm512_store_epi32((__m512i *)e1, velement1);
            _mm512_store_epi32((__m512i *)e2, velement2);

            while ((idx = _mm_tzcnt_32(mask1)) < 16)
            {
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask1 ^= (1 << idx);
            }

            while ((idx = _mm_tzcnt_32(mask2)) < 16)
            {
                bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                *bptr = e2[idx];
                numptr_n[b2[idx]]++;
                mask2 ^= (1 << idx);
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

        logp = update_data.logp[med_B];
        for (j = med_B; j<large_B; j += 16, ptr += 16)
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
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);


            
            vnroot1 = _mm512_sub_epi32(vprime, vroot1);
            vnroot2 = _mm512_sub_epi32(vprime, vroot2);

            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot1, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
            }

            vroot1 = _mm512_add_epi32(vroot1, vprime);
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot1, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            }

            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot2, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
            }

            vroot2 = _mm512_add_epi32(vroot2, vprime);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vroot2, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                    *bptr = e1[idx];
                    numptr_p[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);
            }
            
            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot1, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
            }

            vnroot1 = _mm512_add_epi32(vnroot1, vprime);
            mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            while (mask1 > 0)
            {
                __mmask16 m = mask1;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot1, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vnroot1 = _mm512_mask_add_epi32(vnroot1, mask1, vnroot1, vprime);
                mask1 = _mm512_cmp_epu32_mask(vnroot1, vinterval, _MM_CMPINT_LT);
            }

            _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot2, 15));
            _mm512_store_epi32((__m512i *)e1, 
                _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

            for (idx = 0; idx < 16; idx++)
            {
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
            }

            vnroot2 = _mm512_add_epi32(vnroot2, vprime);
            mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            while (mask2 > 0)
            {
                __mmask16 m = mask2;

                _mm512_store_epi32((__m512i *)b1, _mm512_srli_epi32(vnroot2, 15));
                _mm512_store_epi32((__m512i *)e1, 
                    _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, vblockm1)));

                while ((idx = _mm_tzcnt_32(m)) < 16)
                {
                    bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                    *bptr = e1[idx];
                    numptr_n[b1[idx]]++;
                    m ^= (1 << idx);
                }

                vnroot2 = _mm512_mask_add_epi32(vnroot2, mask2, vnroot2, vprime);
                mask2 = _mm512_cmp_epu32_mask(vnroot2, vinterval, _MM_CMPINT_LT);
            }

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

            // if i do a gather/scatter i get fewer rels/sec.  This is
            // because gather/scatter does not handle collisions... i.e., if
            // more than one index is equal.  if this happens the index is
            // only written once and we effectively lose the relation that did
            // not get a bucket element written.  therefore we still have to do
            // this in a loop.  AVX512 provides instructions that will help.
            _mm512_store_epi32((__m512i *)b1, vbnum1);
            _mm512_store_epi32((__m512i *)b2, vbnum2);
            _mm512_store_epi32((__m512i *)e1, velement1);
            _mm512_store_epi32((__m512i *)e2, velement2);

            while ((idx = _mm_tzcnt_32(mask1)) < 16)
            {
                bptr = sliceptr_p + (b1[idx] << BUCKET_BITS) + numptr_p[b1[idx]];
                *bptr = e1[idx];
                numptr_p[b1[idx]]++;
                mask1 ^= (1 << idx);
            }

            while ((idx = _mm_tzcnt_32(mask2)) < 16)
            {
                bptr = sliceptr_p + (b2[idx] << BUCKET_BITS) + numptr_p[b2[idx]];
                *bptr = e2[idx];
                numptr_p[b2[idx]]++;
                mask2 ^= (1 << idx);
            }

            vroot1 = _mm512_sub_epi32(vprime, vroot1);
            vroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, vblockm1));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, vblockm1));
            mask1 = _mm512_cmp_epu32_mask(vroot1, vinterval, _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, vinterval, _MM_CMPINT_LT);


            _mm512_store_epi32((__m512i *)b1, vbnum1);
            _mm512_store_epi32((__m512i *)b2, vbnum2);
            _mm512_store_epi32((__m512i *)e1, velement1);
            _mm512_store_epi32((__m512i *)e2, velement2);

            while ((idx = _mm_tzcnt_32(mask1)) < 16)
            {
                bptr = sliceptr_n + (b1[idx] << BUCKET_BITS) + numptr_n[b1[idx]];
                *bptr = e1[idx];
                numptr_n[b1[idx]]++;
                mask1 ^= (1 << idx);
            }

            while ((idx = _mm_tzcnt_32(mask2)) < 16)
            {
                bptr = sliceptr_n + (b2[idx] << BUCKET_BITS) + numptr_n[b2[idx]];
                *bptr = e2[idx];
                numptr_n[b2[idx]]++;
                mask2 ^= (1 << idx);
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

// experimental code for processing multiple polynomials simultaneously
// (e.g., as is done in msieve).  I have never been able to make experiments
// like this faster than my specialized bucket sort, and this is no
// exception.  
// left here for "someday".
void nextRoots_32k_knc_polybatch(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //update the roots 
    int *rootupdates = dconf->rootupdates;

    update_t update_data = dconf->update_data;

    uint32 bound = sconf->factor_base->B;

    char *nu = dconf->curr_poly->nu;
    char *gray = dconf->curr_poly->gray;
    int numB = dconf->numB;
    uint32 poly_offset = 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

    lp_bucket *lp_bucket_p = dconf->buckets;
    uint32 med_B = sconf->factor_base->med_B;
    uint32 large_B = sconf->factor_base->large_B;

    uint32 j, interval;
    int k, numblocks;
    uint32 root1, root2, prime;

    int bound_index = 0;
    int check_bound = BUCKET_ALLOC / 2 - 1;
    uint32 bound_val = med_B;
    uint32 *numptr_p, *numptr_n, *sliceptr_p, *sliceptr_n;

    uint32 *bptr;
    int bnum, room;

    uint8 logp = 0;
    
    __declspec(aligned(64)) uint32 e1[16];
    __declspec(aligned(64)) uint32 e2[16];
    __declspec(aligned(64)) uint32 b1[16];
    __declspec(aligned(64)) uint32 b2[16];


    __m512i vshifted_index = _mm512_setr_epi32(
        0, 65536, 131072, 196608,
        262144, 327680, 393216, 458752,
        524288, 589824, 655360, 720896,
        786432, 851968, 917504, 983040);    

    numblocks = sconf->num_blocks;
    interval = numblocks << 15;

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
    for (j = med_B; j<large_B; j += 16)
    {
        int i;
        int p;

        __m512i vprime, vroot1, vroot2, vnroot1, vnroot2, vpval, vbnum1, vbnum2;
        __m512i velement1, velement2, vbsize, vidx;
        __mmask16 mask1, mask2;

        CHECK_NEW_SLICE_BATCH(j);

        vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
        vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
        vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));

        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {            
            vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));

            if (gray[numB + p] > 0)
            {
                mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_sub_epi32(vroot1, vpval);
                vroot2 = _mm512_sub_epi32(vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
            }
            else
            {
                vroot1 = _mm512_add_epi32(vroot1, vpval);
                vroot2 = _mm512_add_epi32(vroot2, vpval);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            }
            
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);

            _mm512_store_epi32((__m512i *)b1, vroot1);
            _mm512_store_epi32((__m512i *)b2, vroot2);
            _mm512_store_epi32((__m512i *)e1, velement1);

            for (i = 0; i < 16; i++)
            {
                while (b1[i] < interval)
                {
                    bnum = b1[i] >> 15;
                    bptr = sliceptr_p + (bnum << BUCKET_BITS) + numptr_p[bnum];
                    *bptr = e1[i] | (b1[i] & 32767);
                    numptr_p[bnum]++;
                    b1[i] += update_data.prime[j + i];
                }

                while (b2[i] < interval)
                {
                    bnum = b2[i] >> 15;
                    bptr = sliceptr_p + (bnum << BUCKET_BITS) + numptr_p[bnum];
                    *bptr = e1[i] | (b2[i] & 32767);
                    numptr_p[bnum]++;
                    b2[i] += update_data.prime[j + i];
                }
            }


            vnroot1 = _mm512_sub_epi32(vprime, vroot1);
            vnroot2 = _mm512_sub_epi32(vprime, vroot2);

            _mm512_store_epi32((__m512i *)b1, vnroot1);
            _mm512_store_epi32((__m512i *)b2, vnroot2);

            for (i = 0; i < 16; i++)
            {
                while (b1[i] < interval)
                {
                    bnum = b1[i] >> 15;
                    bptr = sliceptr_n + (bnum << BUCKET_BITS) + numptr_n[bnum];
                    *bptr = e1[i] | (b1[i] & 32767);
                    numptr_n[bnum]++;
                    b1[i] += update_data.prime[j + i];
                }

                while (b2[i] < interval)
                {
                    bnum = b2[i] >> 15;
                    bptr = sliceptr_n + (bnum << BUCKET_BITS) + numptr_n[bnum];
                    *bptr = e1[i] | (b2[i] & 32767);
                    numptr_n[bnum]++;
                    b2[i] += update_data.prime[j + i];
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


    logp = update_data.logp[large_B];
    for (j = large_B; j<bound; j += 16)
    {
        int i;
        int p;

        __m512i vprime, vroot1, vroot2, vnroot1, vnroot2, vpval, vbnum1, vbnum2;
        __m512i velement1, velement2, vbsize, vidx;
        __mmask16 mask1, mask2;

        CHECK_NEW_SLICE_BATCH(j);

        vprime = _mm512_load_epi32((__m512i *)(&update_data.prime[j]));
        vroot1 = _mm512_load_epi32((__m512i *)(&update_data.firstroots1[j]));
        vroot2 = _mm512_load_epi32((__m512i *)(&update_data.firstroots2[j]));

        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            vpval = _mm512_load_epi32((__m512i *)(&rootupdates[(nu[numB + p] - 1) * bound + j]));

            if (gray[numB + p] > 0)
            {
                mask1 = _mm512_cmp_epu32_mask(vpval, vroot1, _MM_CMPINT_GT);
                mask2 = _mm512_cmp_epu32_mask(vpval, vroot2, _MM_CMPINT_GT);
                vroot1 = _mm512_sub_epi32(vroot1, vpval);
                vroot2 = _mm512_sub_epi32(vroot2, vpval);
                vroot1 = _mm512_mask_add_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_add_epi32(vroot2, mask2, vroot2, vprime);
            }
            else
            {
                vroot1 = _mm512_add_epi32(vroot1, vpval);
                vroot2 = _mm512_add_epi32(vroot2, vpval);
                mask1 = _mm512_cmp_epu32_mask(vroot1, vprime, _MM_CMPINT_GE);
                mask2 = _mm512_cmp_epu32_mask(vroot2, vprime, _MM_CMPINT_GE);
                vroot1 = _mm512_mask_sub_epi32(vroot1, mask1, vroot1, vprime);
                vroot2 = _mm512_mask_sub_epi32(vroot2, mask2, vroot2, vprime);
            }
            
            vbnum1 = _mm512_srli_epi32(vroot1, 15);
            vbnum2 = _mm512_srli_epi32(vroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot2, _mm512_set1_epi32(32767)));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vroot1, _mm512_set1_epi32(32767)));
            mask1 = _mm512_cmp_epu32_mask(vroot1, _mm512_set1_epi32(interval), _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vroot2, _mm512_set1_epi32(interval), _MM_CMPINT_LT);

            _mm512_store_epi32((__m512i *)b1, vbnum1);
            _mm512_store_epi32((__m512i *)b2, vbnum2);
            _mm512_store_epi32((__m512i *)e1, velement1);
            _mm512_store_epi32((__m512i *)e2, velement2);

            for (i = 0; i < 16; i++)
            {
                if (mask1 & (1 << i))
                {
                    bptr = sliceptr_p + (b1[i] << BUCKET_BITS) + numptr_p[b1[i]];
                    *bptr = e1[i];
                    numptr_p[b1[i]]++;
                }

                if (mask2 & (1 << i))
                {
                    bptr = sliceptr_p + (b2[i] << BUCKET_BITS) + numptr_p[b2[i]];
                    *bptr = e2[i];
                    numptr_p[b2[i]]++;
                }
            }

            vnroot1 = _mm512_sub_epi32(vprime, vroot1);
            vnroot2 = _mm512_sub_epi32(vprime, vroot2);
            vbnum1 = _mm512_srli_epi32(vnroot1, 15);
            vbnum2 = _mm512_srli_epi32(vnroot2, 15);
            velement1 = _mm512_add_epi32(_mm512_set1_epi32((j - bound_val) << 16), vshifted_index);
            velement2 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot2, _mm512_set1_epi32(32767)));
            velement1 = _mm512_or_epi32(velement1, _mm512_and_epi32(vnroot1, _mm512_set1_epi32(32767)));
            mask1 = _mm512_cmp_epu32_mask(vnroot1, _mm512_set1_epi32(interval), _MM_CMPINT_LT);
            mask2 = _mm512_cmp_epu32_mask(vnroot2, _mm512_set1_epi32(interval), _MM_CMPINT_LT);

            _mm512_store_epi32((__m512i *)b1, vbnum1);
            _mm512_store_epi32((__m512i *)b2, vbnum2);
            _mm512_store_epi32((__m512i *)e1, velement1);
            _mm512_store_epi32((__m512i *)e2, velement2);

            for (i = 0; i < 16; i++)
            {
                if (mask1 & (1 << i))
                {
                    bptr = sliceptr_n + (b1[i] << BUCKET_BITS) + numptr_n[b1[i]];
                    *bptr = e1[i];
                    numptr_n[b1[i]]++;
                }

                if (mask2 & (1 << i))
                {
                    bptr = sliceptr_n + (b2[i] << BUCKET_BITS) + numptr_n[b2[i]];
                    *bptr = e2[i];
                    numptr_n[b2[i]]++;
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
