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
#include "factor.h"
#include "util.h"
#include "common.h"
#include "tdiv_macros_common.h"

#ifdef TARGET_KNC

#include <immintrin.h>

#define INIT_RESIEVE \
        vcorr = _mm512_set1_epi32(32768 - block_loc); \
        vroot1 = _mm512_extload_epi32((__m512i *)(&fbc->root1[i]),             \
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);      \
        vroot2 = _mm512_extload_epi32((__m512i *)(&fbc->root2[i]),             \
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);      \
        vprimes = _mm512_extload_epi32((__m512i *)(&fbc->prime[i]),             \
            _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NONE);      \
        vroot1 = _mm512_add_epi32(vroot1, vcorr); \
        vroot2 = _mm512_add_epi32(vroot2, vcorr);

#define RESIEVE_8X_14BIT_MAX \
        vroot1 = _mm512_sub_epi32(vroot1, vprimes); \
        vroot2 = _mm512_sub_epi32(vroot2, vprimes); \
        resmask = _mm512_cmp_epu32_mask(vzero, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot2, _MM_CMPINT_EQ); \
        vroot1 = _mm512_sub_epi32(vroot1, vprimes); \
        vroot2 = _mm512_sub_epi32(vroot2, vprimes); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot2, _MM_CMPINT_EQ); \
        vroot1 = _mm512_sub_epi32(vroot1, vprimes); \
        vroot2 = _mm512_sub_epi32(vroot2, vprimes); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot2, _MM_CMPINT_EQ); \
        vroot1 = _mm512_sub_epi32(vroot1, vprimes); \
        vroot2 = _mm512_sub_epi32(vroot2, vprimes); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot2, _MM_CMPINT_EQ);

#define RESIEVE_8X_15BIT_MAX \
        vroot1 = _mm512_sub_epi32(vroot1, vprimes); \
        vroot2 = _mm512_sub_epi32(vroot2, vprimes); \
        resmask = _mm512_cmp_epu32_mask(vzero, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot2, _MM_CMPINT_EQ); \
        vroot1 = _mm512_sub_epi32(vroot1, vprimes); \
        vroot2 = _mm512_sub_epi32(vroot2, vprimes); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot2, _MM_CMPINT_EQ);

#define RESIEVE_8X_16BIT_MAX \
        vroot1 = _mm512_sub_epi32(vroot1, vprimes); \
        vroot2 = _mm512_sub_epi32(vroot2, vprimes); \
        resmask = _mm512_cmp_epu32_mask(vzero, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vzero, vroot2, _MM_CMPINT_EQ);


//#define SIQSDEBUG 1

/*
We are given an array of bytes that has been sieved.  The basic trial 
division strategy is as follows:

1) Scan through the array and 'mark' locations that meet criteria 
indicating they may factor completely over the factor base.  

2) 'Filter' the marked locations by trial dividing by small primes
that we did not sieve.  These primes are all less than 256.  If after
removing small primes the location does not meet another set of criteria,
remove it from the 'marked' list (do not subject it to further trial
division).

3) Divide out primes from the factor base between 256 and 2^13 or 2^14, 
depending on the version (2^13 for 32k version, 2^14 for 64k).  

4) Resieve primes between 2^{13|14} and 2^16, max.  

5) Primes larger than 2^16 will have been bucket sieved.  Remove these
by scanning the buckets for sieve hits equal to the current block location.

6) If applicable/appropriate, factor a remaining composite with squfof

this file contains code implementing 4)


*/

void resieve_medprimes_32k_knc(uint8 parity, uint32 poly_id, uint32 bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	int i;
	uint32 bound, report_num;
	int smooth_num;
	uint32 *fb_offsets;
	sieve_fb_compressed *fbc;
	fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
	uint32 block_loc;
	uint16 *corrections = dconf->corrections;


    __m512i vzero = _mm512_setzero_epi32();


	fullfb_ptr = fullfb;
	if (parity)
	{
		fbc = dconf->comp_sieve_n;
	}
	else
	{
		fbc = dconf->comp_sieve_p;
	}		

#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	for (report_num = 0; report_num < dconf->num_reports; report_num++)
	{

		if (!dconf->valid_Qs[report_num])
			continue;

		// pull the details of this report to get started.
		fb_offsets = &dconf->fb_offsets[report_num][0];
		smooth_num = dconf->smooth_num[report_num];
		block_loc = dconf->reports[report_num];
		
		// where tdiv_medprimes left off
		i = sconf->factor_base->fb_13bit_B;

        bound = sconf->factor_base->fb_14bit_B;        

        while ((uint32)i < bound)
        {
            int p = (int)fbc->prime[i];
            int r1 = (int)fbc->root1[i] + 32768 - block_loc;
            int r2 = (int)fbc->root2[i] + 32768 - block_loc;

            if ((i & 15) == 0)
                break;

            r1 -= p;
            r2 -= p;
            if (r1 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            if (r2 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            r1 -= p;
            r2 -= p;
            if (r1 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            if (r2 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            r1 -= p;
            r2 -= p;
            if (r1 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            if (r2 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            r1 -= p;
            r2 -= p;
            if (r1 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            if (r2 == 0){ DIVIDE_RESIEVED_PRIME(0); i++; continue;}
            i++;
        }

        /*
        if (((bound - i) & 15) != 0)
        {            
            bound -= 8;
        }
        */

        while ((uint32)i < bound)
        {
            __mmask16 resmask;
            __m512i vroot1, vroot2, vcorr, vprimes;
            int idx;

            INIT_RESIEVE;
            RESIEVE_8X_14BIT_MAX;

            //if (resmask == 0)
            //{
            //    i += 16;
            //    continue;
            //}

            while ((idx = _mm_tzcnt_32(resmask)) < 16) {
                DIVIDE_RESIEVED_PRIME(idx);
                resmask ^= (1 << idx);
            }

            i += 16;
        }

        bound = sconf->factor_base->fb_15bit_B;        

        while ((uint32)i < bound)
        {
            __mmask16 resmask;
            __m512i vroot1, vroot2, vcorr, vprimes;
            int idx;

            INIT_RESIEVE;
            RESIEVE_8X_15BIT_MAX;

            while ((idx = _mm_tzcnt_32(resmask)) < 16) {
                DIVIDE_RESIEVED_PRIME(idx);
                resmask ^= (1 << idx);
            }

            i += 16;
        }

        bound = sconf->factor_base->med_B;        

        while ((uint32)i < bound)
        {
            __mmask16 resmask;
            __m512i vroot1, vroot2, vcorr, vprimes;
            int idx;

            INIT_RESIEVE;
            RESIEVE_8X_16BIT_MAX;
            
            while ((idx = _mm_tzcnt_32(resmask)) < 16) {
                DIVIDE_RESIEVED_PRIME(idx);
                resmask ^= (1 << idx);
            }

            i += 16;
        }

		// after either resieving or standard trial division, record
		// how many factors we've found so far.
		dconf->smooth_num[report_num] = smooth_num;	

	}
			
#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
    TF_STG4 += my_difftime (&qs_timing_start, &qs_timing_stop);
#endif

	return;
}

#endif
