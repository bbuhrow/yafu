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
#include "ytools.h"
#include "common.h"
#include "tdiv_macros_common.h"
#include "tdiv_macros_32k.h"

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

this file contains code implementing 3)


*/

#ifdef USE_AVX512F

#include <immintrin.h>

#define MOD_INIT_16X	\
		vblksz = _mm512_set1_epi32(32768); \
    vblkloc = _mm512_set1_epi32(block_loc); \
    vindex = _mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);

#define MOD_CMP_16X(i,xtra_bits)	\
        if (1) { \
        int idx; __mmask16 resmask; \
        vprimes = _mm512_i32gather_epi32(vindex, (__m512i *)(&fbc->prime[i]), _MM_SCALE_2);      \
        vcorr = _mm512_i32gather_epi32(vindex, (__m512i *)(&fullfb_ptr->correction[i]), _MM_SCALE_2);      \
        vinv = _mm512_i32gather_epi32(vindex, (__m512i *)(&fullfb_ptr->small_inv[i]), _MM_SCALE_2);      \
        vroot1 = _mm512_i32gather_epi32(vindex, (__m512i *)(&fbc->root1[i]), _MM_SCALE_2);      \
        vroot2 = _mm512_i32gather_epi32(vindex, (__m512i *)(&fbc->root2[i]), _MM_SCALE_2);      \
        vcorr = _mm512_add_epi32(vcorr, _mm512_sub_epi32(vblksz, vblkloc)); \
        vcorr = _mm512_mullo_epi32(vcorr, vinv);    \
        vcorr = _mm512_srli_epi32(vcorr, xtra_bits);    \
        vcorr = _mm512_mullo_epi32(vcorr, vprimes); \
        vcorr = _mm512_add_epi32(vcorr, _mm512_add_epi32(vblkloc, vprimes));    \
        vcorr = _mm512_sub_epi32(vcorr, vblksz);    \
        resmask = _mm512_cmp_epu32_mask(vcorr, vroot1, _MM_CMPINT_EQ); \
        resmask |= _mm512_cmp_epu32_mask(vcorr, vroot2, _MM_CMPINT_EQ); \
        while (_BitScanForward(&idx, resmask)) {    \
            buffer32[tmp3++] = (i) + idx; \
            resmask ^= (1 << idx); \
                        } }


void tdiv_medprimes_32k_knl(uint8 parity, uint32 poly_id, uint32 bnum,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //we have flagged this sieve offset as likely to produce a relation
    //nothing left to do now but check and see.
    int i;
    uint32 bound, tmp, prime, root1, root2, report_num;
    int smooth_num;
    uint32 *fb_offsets;
    sieve_fb_compressed *fbc;
    fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
    uint32 block_loc;
    uint16 buffer[32];
    uint32 tmp3 = 0;
    int r;
    int starti;


    __m512i vblksz, vblkloc, vprimes, vcorr, vinv, vroot1, vroot2, vindex;
    uint32 buffer32[32];
    int k;


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

        // for each report, we trial divide, and then either trial divide or
        // resieve.  for the first trial division step, we have either
        // unrolled C routines or SIMD assembly routines to choose from.  The
        // second trial division step only has a C routine - the more optimized
        // path is resieving.
        //
        // the basic idea of trial division is to test if the block location
        // in question lies on the arithmetic progression of a prime.  By the time
        // we get to this routine, the arithmetic progression has been reset for
        // the next sieve block, so we have to do a few manipulations to revert
        // to the "real" progression (add the blocksize back to the roots).  
        // there are various methods for doing the test depending on the size of
        // the primes (and therefore how many times it could have possibly landed
        // in the sieve block).  The most straightforward, and the method the SIMD
        // assembly uses, is to see if the roots (adjusted for the "real" progression)
        // minus the block location in question, divided by the prime, is zero.  if
        // so, this block location is on the arithmetic progression of that prime,
        // and we can proceed to trial divide the prime into Q(x) for this sieve hit.
        // 
        // the basic idea of resieving is to start from the end of the "real"
        // arithmetic progression, repeatedly subtract each prime, and test after
        // each subtraction if we've hit the sieve location in question.  if so,
        // we know this location is on the prime's arithmetic progression and can 
        // proceed to trial divide.  For "large enough" primes, this is very efficient
        // because we only need to do a few subtractions and tests instead of a 
        // division.  Since we are not really doing division tests, and instead are
        // doing multiplication by inverses, and futhermore since we might be doing those
        // multiplications 8 at a time using SIMD, resieving is only a win for the 
        // very largest primes less than 16 bits in size.  
        //
        // for OS/architecture/compilers where resieving isn't implemented, there are
        // further trial division steps instead.  These are more efficient than
        // the "check for exact division of a difference" method described above,
        // but are only implemented in portable C.  See code below for more detail.

        // pull the details of this report to get started.
        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        i = sconf->sieve_small_fb_start;

        // 8x trial division
        bound = sconf->factor_base->fb_10bit_B;

        // single-up test until i is a multiple of 16
        while (((uint32)i < bound) && ((i & 15) != 0))
        {
            prime = fbc->prime[i];

            //tmp = distance from this sieve block offset to the end of the block
            tmp = 32768 - block_loc;

            //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
            //block, which is the state of the sieve now.
            tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                * (uint64)fullfb_ptr->small_inv[i]) >> 24);
            tmp = block_loc + tmp*prime;
            tmp = tmp - 32768;

            //tmp = advance the offset to where it should be after the interval, and
            //check to see if that's where either of the roots are now.  if so, then
            //this offset is on the progression of the sieve for this prime
            if ((tmp == fbc->root1[i]) || (tmp == fbc->root2[i])) {
                buffer32[tmp3++] = i;
            }
            i++;
        }


        starti = i;

        bound = sconf->factor_base->fb_10bit_B;
        if (((bound - i) & 15) != 0)
        {            
            bound -= 8;
        }

        MOD_INIT_16X;

        while (i < bound)
        {
            MOD_CMP_16X(i,24);            
            i += 16;
        }
        
        bound = sconf->factor_base->fb_12bit_B;

        // advance past the 10-bit boundary until we are aligned again
        if (i != sconf->factor_base->fb_10bit_B)
        {
            for (k = 0; k < 8; k++, i++)
            {
                uint32 b = block_loc;
                prime = fbc->prime[i];

                tmp = 32768 - block_loc;
                tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                    * (uint64)fullfb_ptr->small_inv[i]) >> 24);
                tmp = block_loc + tmp*prime;
                tmp = tmp - 32768;

                if ((tmp == fbc->root1[i]) || (tmp == fbc->root2[i])) {
                    buffer32[tmp3++] = i;
                }
            }
            for (k = 0; k < 8; k++, i++)
            {
                uint32 b = block_loc;
                prime = fbc->prime[i];

                tmp = 32768 - block_loc;
                tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                    * (uint64)fullfb_ptr->small_inv[i]) >> 26);
                tmp = block_loc + tmp*prime;
                tmp = tmp - 32768;

                if ((tmp == fbc->root1[i]) || (tmp == fbc->root2[i])) {
                    buffer32[tmp3++] = i;
                }
            }
        }

        if (((bound - i) & 15) != 0)
        {
            bound -= 8;
        }

        while (i < bound)
        {
            MOD_CMP_16X(i,26);
            i += 16;
        }

        bound = sconf->factor_base->fb_13bit_B;

        // advance past the 12-bit boundary until we are aligned again
        if (i != sconf->factor_base->fb_12bit_B)
        {
            for (k = 0; k < 8; k++, i++)
            {
                uint32 b = block_loc;
                prime = fbc->prime[i];

                tmp = 32768 - block_loc;
                tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                    * (uint64)fullfb_ptr->small_inv[i]) >> 26);
                tmp = block_loc + tmp*prime;
                tmp = tmp - 32768;

                if ((tmp == fbc->root1[i]) || (tmp == fbc->root2[i])) {
                    buffer32[tmp3++] = i;
                }
            }
            for (k = 0; k < 8; k++, i++)
            {
                uint32 b = block_loc;
                prime = fbc->prime[i];

                tmp = 32768 - block_loc;
                tmp = 1 + (uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
                    * (uint64)fullfb_ptr->small_inv[i]) >> 28);
                tmp = block_loc + tmp*prime;
                tmp = tmp - 32768;

                if ((tmp == fbc->root1[i]) || (tmp == fbc->root2[i])) {
                    buffer32[tmp3++] = i;
                }
            }
        }

        if (((bound - i) & 15) != 0)
        {
            bound -= 8;
        }

        while (i < bound)
        {
            MOD_CMP_16X(i,28);
            i += 16;
        }

        bound = sconf->factor_base->fb_13bit_B;
        while ((uint32)i < bound)
        {
            uint32 b = block_loc;
            prime = fbc->prime[i];

            while (b < 32768) {b += prime;} tmp = b - 32768;
            if ((tmp == fbc->root1[i]) || (tmp == fbc->root2[i])) {
                buffer32[tmp3++] = i;
            }
            i++;
        }

        for (r = 0; r < tmp3; r++)
        {
            while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[buffer32[r]]) == 0)
            {
                fb_offsets[++smooth_num] = buffer32[r];
                mpz_tdiv_q_ui(dconf->Qvals[report_num],
                    dconf->Qvals[report_num], fbc->prime[buffer32[r]]);
            }
        }

        // either after 8x SSE2 ASM, or standard trial division, record
        // how many factors we've found so far
        dconf->smooth_num[report_num] = smooth_num;

    }

#ifdef QS_TIMING
    gettimeofday(&qs_timing_stop, NULL);
    TF_STG2 += ytools_difftime(&qs_timing_start, &qs_timing_stop);
#endif

    return;
}


#endif

