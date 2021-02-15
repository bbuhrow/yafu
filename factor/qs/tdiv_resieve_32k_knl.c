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

#include "common.h"

#if defined( USE_AVX512F )

#include "qs_impl.h"
#include <immintrin.h>

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


void resieve_medprimes_32k_knl(uint32_t *reports, uint32_t num_reports, 
    uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //we have flagged this sieve offset as likely to produce a relation
    //nothing left to do now but check and see.
    int i, j, k;
    int idx, id;
    uint32_t bound, report_num;
    int smooth_num;
    uint32_t *fb_offsets;
    sieve_fb_compressed *fbc;
    fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
    uint32_t block_loc;
    uint16_t *corrections = dconf->corrections;
    uint16_t buffer[16];
    uint32_t result = 0;
    __m512i vlomask = _mm512_set1_epi32(0x0000ffff);
    __m512i vdoubleblksz = _mm512_set1_epi32(0x80008000);
    __m512i zero = _mm512_setzero_epi32();
    __mmask16 meq1, meq2, lanemask1, lanemask2;

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

    // printf("running knl resieve targetting blocks: ");
    // for (j = 0; j < num_reports; j++)
      // printf("%u ", dconf->reports[reports[j]]);
    // printf("\n");
    
        
    for (j = 0; j < num_reports; j++)
    {
        __m512i vblock = _mm512_set1_epi32(dconf->reports[reports[j]]);  

        meq1 = meq2 = 0;        
        lanemask1 = 0xffff;
        lanemask2 = 0xffff;
            
        // where tdiv_medprimes left off
        for (i = sconf->factor_base->fb_13bit_B; i < sconf->factor_base->med_B; i += 32)
        {
            // load 32 16-bit roots/primes.
            __m512i vprime = _mm512_load_epi32(fbc->prime + i);
            __m512i vroot1 = _mm512_load_epi32(fbc->root1 + i);
            __m512i vroot2 = _mm512_load_epi32(fbc->root2 + i);               
            
            // add blocksize and subtract primes.
            __m512i vresieve16a, vresieve16b, vprime16;
            __m512i vresieve1 = _mm512_add_epi32(vroot1, vdoubleblksz);
            __m512i vresieve2 = _mm512_add_epi32(vroot2, vdoubleblksz);
        
            if (fbc->prime[i] > sconf->factor_base->fb_14bit_B)
                break;
          
            vresieve1 = _mm512_sub_epi32(vresieve1, vprime);
            vresieve2 = _mm512_sub_epi32(vresieve2, vprime);
            
            // process the evens
            // check if equal to all reported block locations
            vprime16 = _mm512_and_epi32(vprime, vlomask);
        
        
            vresieve16a = _mm512_and_epi32(vresieve1, vlomask);
            vresieve16b = _mm512_and_epi32(vresieve2, vlomask);

            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);

            // scan to find all set bits and trial divide with the appropriate prime
            while ((idx = _tzcnt_u32(meq1)) < 16)
            {
                id = i+idx*2;
                while (mpz_tdiv_ui(dconf->Qvals[reports[j]], fbc->prime[id]) == 0)
                {
                    dconf->fb_offsets[reports[j]][++dconf->smooth_num[reports[j]]] = id;
                    mpz_tdiv_q_ui(dconf->Qvals[reports[j]], dconf->Qvals[reports[j]], fbc->prime[id]);
                }
                meq1 ^= (1 << idx);
            }                           
        
            // process the odds
            vprime16 = _mm512_srli_epi32(vprime, 16);

            vresieve16a = _mm512_srli_epi32(vresieve1, 16);
            vresieve16b = _mm512_srli_epi32(vresieve2, 16);            

            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);

        }

        for (i = sconf->factor_base->fb_13bit_B; i < sconf->factor_base->med_B; i += 32)
        {
            // load 32 16-bit roots/primes.
            __m512i vprime = _mm512_load_epi32(fbc->prime + i);
            __m512i vroot1 = _mm512_load_epi32(fbc->root1 + i);
            __m512i vroot2 = _mm512_load_epi32(fbc->root2 + i);               
            
            // add blocksize and subtract primes.
            __m512i vresieve16a, vresieve16b, vprime16;
            __m512i vresieve1 = _mm512_add_epi32(vroot1, vdoubleblksz);
            __m512i vresieve2 = _mm512_add_epi32(vroot2, vdoubleblksz);
        
            if (fbc->prime[i] > sconf->factor_base->fb_14bit_B)
                break;
          
            vresieve1 = _mm512_sub_epi32(vresieve1, vprime);
            vresieve2 = _mm512_sub_epi32(vresieve2, vprime);
            
            // process the evens
            // check if equal to all reported block locations
            vprime16 = _mm512_and_epi32(vprime, vlomask);
        
        
            vresieve16a = _mm512_and_epi32(vresieve1, vlomask);
            vresieve16b = _mm512_and_epi32(vresieve2, vlomask);

            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
                    
            // process the odds
            vprime16 = _mm512_srli_epi32(vprime, 16);

            vresieve16a = _mm512_srli_epi32(vresieve1, 16);
            vresieve16b = _mm512_srli_epi32(vresieve2, 16);            

            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);
            
            // step backward by prime
            vresieve16a = _mm512_mask_sub_epi32(vresieve16a, lanemask1, vresieve16a, vprime16);
            vresieve16b = _mm512_mask_sub_epi32(vresieve16b, lanemask2, vresieve16b, vprime16);
            
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);

        }
        
        for (i = sconf->factor_base->fb_13bit_B; i < sconf->factor_base->med_B; i += 32)
        {
            // load 32 16-bit roots/primes.
            __m512i vprime = _mm512_load_epi32(fbc->prime + i);
            __m512i vroot1 = _mm512_load_epi32(fbc->root1 + i);
            __m512i vroot2 = _mm512_load_epi32(fbc->root2 + i);               
            
            // add blocksize and subtract primes.
            __m512i vresieve16a, vresieve16b, vprime16;
            __m512i vresieve1 = _mm512_add_epi32(vroot1, vdoubleblksz);
            __m512i vresieve2 = _mm512_add_epi32(vroot2, vdoubleblksz);
        
            if (fbc->prime[i] > sconf->factor_base->fb_14bit_B)
                break;
          
            vresieve1 = _mm512_sub_epi32(vresieve1, vprime);
            vresieve2 = _mm512_sub_epi32(vresieve2, vprime);
            
            // process the evens
            // check if equal to all reported block locations
            vprime16 = _mm512_and_epi32(vprime, vlomask);
        
        
            vresieve16a = _mm512_and_epi32(vresieve1, vlomask);
            vresieve16b = _mm512_and_epi32(vresieve2, vlomask);

            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq1 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);

        
            // process the odds
            vprime16 = _mm512_srli_epi32(vprime, 16);

            vresieve16a = _mm512_srli_epi32(vresieve1, 16);
            vresieve16b = _mm512_srli_epi32(vresieve2, 16);            

            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask1, vblock, vresieve16a, _MM_CMPINT_EQ);
            meq2 |= _mm512_mask_cmp_epi32_mask(lanemask2, vblock, vresieve16b, _MM_CMPINT_EQ);

        }
        
        
    }
        
        
        
#ifdef QS_TIMING
    gettimeofday(&qs_timing_stop, NULL);
    TF_STG4 += ytools_difftime(&qs_timing_start, &qs_timing_stop);
#endif

    return;
}

#endif // USE_AVX2
