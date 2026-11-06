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
#include "tdiv_macros_common.h"
#include "tdiv_macros_32k.h"


#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#elif defined( USE_AVX2 ) || defined( USE_AVX512F ) || defined( USE_SSE41 ) || defined( D_HAS_SSE2 )
#include <immintrin.h>
#endif

#if !defined(__aarch64__) && (defined(USE_SSE41) || defined(D_HAS_SSE2))
#include "tdiv_med_sse.h"
#endif

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

void tdiv_medprimes_32k(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //we have flagged this sieve offset as likely to produce a relation
    //nothing left to do now but check and see.
    int i;
    uint32_t bound, tmp, prime, root1, root2, report_num;
    int smooth_num;
    uint32_t *fb_offsets;
    sieve_fb_compressed *fbc;
    fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
    uint32_t block_loc;
    uint16_t buffer[32];
    uint32_t tmp3 = 0;
    int r;

#ifdef USE_8X_MOD_ASM
    uint16_t *bl_sizes = dconf->bl_sizes;
    uint16_t *bl_locs = dconf->bl_locs;
#endif

    fullfb_ptr = fullfb;
    if (parity)
    {
        fbc = dconf->comp_sieve_n;
    }
    else
    {
        fbc = dconf->comp_sieve_p;
    }

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

        // single-up test until i is a multiple of 8
        while ((uint32_t)i < bound && ((i & 7) != 0))
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            //tmp = distance from this sieve block offset to the end of the block
            tmp = 32768 - block_loc;

            //tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
            //block, which is the state of the sieve now.
            tmp = 1 + (uint32_t)(((uint64_t)(tmp + fullfb_ptr->correction[i])
                * (uint64_t)fullfb_ptr->small_inv[i]) >> FOGSHIFT);
            tmp = block_loc + tmp*prime;
            tmp = tmp - 32768;

            //tmp = advance the offset to where it should be after the interval, and
            //check to see if that's where either of the roots are now.  if so, then
            //this offset is on the progression of the sieve for this prime
            if (tmp == root1 || tmp == root2)
            {
                //it will divide Q(x).  do so as many times as we can.
                DIVIDE_ONE_PRIME;
            }
            i++;
        }

#if defined(USE_8X_MOD_ASM)

        bl_sizes[0] = 32768;
        bl_sizes[1] = 32768;
        bl_sizes[2] = 32768;
        bl_sizes[3] = 32768;
        bl_sizes[4] = 32768;
        bl_sizes[5] = 32768;
        bl_sizes[6] = 32768;
        bl_sizes[7] = 32768;

        bl_locs[0] = block_loc;
        bl_locs[1] = block_loc;
        bl_locs[2] = block_loc;
        bl_locs[3] = block_loc;
        bl_locs[4] = block_loc;
        bl_locs[5] = block_loc;
        bl_locs[6] = block_loc;
        bl_locs[7] = block_loc;

        tmp3 = 0;

        MOD_INIT_8X;

        while ((uint32_t)i < bound)
        {


#if defined(D_HAS_SSE2) && defined(GCC_ASM64X)

            MOD_CMP_8X_vec("8");

            i += 8;

#else

            uint32_t tmp3 = 0;
            int r;

#ifdef _MSC_VER
            MOD_CMP_8X(8);
#else
            MOD_CMP_8X("8");
#endif


            if (tmp3 == 0)
            {
                i += 8;
                continue;
            }

            if (tmp3 & 0x2)
            {
                DIVIDE_RESIEVED_PRIME(0);
            }

            if (tmp3 & 0x8)
            {
                DIVIDE_RESIEVED_PRIME(1);
            }

            if (tmp3 & 0x20)
            {
                DIVIDE_RESIEVED_PRIME(2);
            }

            if (tmp3 & 0x80)
            {
                DIVIDE_RESIEVED_PRIME(3);
            }

            if (tmp3 & 0x200)
            {
                DIVIDE_RESIEVED_PRIME(4);
            }

            if (tmp3 & 0x800)
            {
                DIVIDE_RESIEVED_PRIME(5);
            }

            if (tmp3 & 0x2000)
            {
                DIVIDE_RESIEVED_PRIME(6);
            }

            if (tmp3 & 0x8000)
            {
                DIVIDE_RESIEVED_PRIME(7);
            }

            i += 8;

#endif


        }

        bound = sconf->factor_base->fb_12bit_B;

        while ((uint32_t)i < bound)
        {

#if defined(D_HAS_SSE2) && defined(GCC_ASM64X)

            MOD_CMP_8X_vec("10");

            i += 8;

#else


            uint32_t tmp3 = 0;
            int r;

#ifdef _MSC_VER
            MOD_CMP_8X(10);
#else
            MOD_CMP_8X("10");
#endif


            if (tmp3 == 0)
            {
                i += 8;
                continue;
            }

            if (tmp3 & 0x2)
            {
                DIVIDE_RESIEVED_PRIME(0);
            }

            if (tmp3 & 0x8)
            {
                DIVIDE_RESIEVED_PRIME(1);
            }

            if (tmp3 & 0x20)
            {
                DIVIDE_RESIEVED_PRIME(2);
            }

            if (tmp3 & 0x80)
            {
                DIVIDE_RESIEVED_PRIME(3);
            }

            if (tmp3 & 0x200)
            {
                DIVIDE_RESIEVED_PRIME(4);
            }

            if (tmp3 & 0x800)
            {
                DIVIDE_RESIEVED_PRIME(5);
            }

            if (tmp3 & 0x2000)
            {
                DIVIDE_RESIEVED_PRIME(6);
            }

            if (tmp3 & 0x8000)
            {
                DIVIDE_RESIEVED_PRIME(7);
            }

            i += 8;

#endif

        }

        bound = sconf->factor_base->fb_13bit_B;

        while ((uint32_t)i < bound)
        {

#if defined(D_HAS_SSE2) && defined(GCC_ASM64X)

            MOD_CMP_8X_vec("12");

            i += 8;

#else


            uint32_t tmp3 = 0;
            int r;

#ifdef _MSC_VER
            MOD_CMP_8X(12);
#else
            MOD_CMP_8X("12");
#endif


            if (tmp3 == 0)
            {
                i += 8;
                continue;
            }

            if (tmp3 & 0x2)
            {
                DIVIDE_RESIEVED_PRIME(0);
            }

            if (tmp3 & 0x8)
            {
                DIVIDE_RESIEVED_PRIME(1);
            }

            if (tmp3 & 0x20)
            {
                DIVIDE_RESIEVED_PRIME(2);
            }

            if (tmp3 & 0x80)
            {
                DIVIDE_RESIEVED_PRIME(3);
            }

            if (tmp3 & 0x200)
            {
                DIVIDE_RESIEVED_PRIME(4);
            }

            if (tmp3 & 0x800)
            {
                DIVIDE_RESIEVED_PRIME(5);
            }

            if (tmp3 & 0x2000)
            {
                DIVIDE_RESIEVED_PRIME(6);
            }

            if (tmp3 & 0x8000)
            {
                DIVIDE_RESIEVED_PRIME(7);
            }

            i += 8;

#endif

        }

        TDIV_MED_CLEAN;

#if defined(D_HAS_SSE2) && defined(GCC_ASM64X)
        for (r = 0; r < tmp3; r++)
        {
            DIVIDE_RESIEVED_PRIME_2((buffer[r]));
        }
#endif

#else

        // the bound before resieving takes over
        bound = sconf->factor_base->fb_13bit_B;

        // generic C version of trial division (by multiplication) up through 13 bits.
        for (; i<bound; i++)
        {
            uint32_t tmp1 = 32768 - block_loc;
            uint64_t q64;

            tmp1 = tmp1 + fullfb_ptr->correction[i];
            q64 = (uint64_t)tmp1 * (uint64_t)fullfb_ptr->small_inv[i];
            tmp1 = q64 >> FOGSHIFT;
            tmp1 = tmp1 + 1;
            tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];

            tmp1 = tmp1 - 32768;
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            // to vectorize, probably need to mark this prime in some kind
            // of byte array, then go back and scan/divide afterwards.
            // the vectorizer won't like this iterative function call.
            if (tmp1 == root1 || tmp1 == root2)
            {
                //it will divide Q(x).  do so as many times as we can.
                DIVIDE_ONE_PRIME;
            }
        }

#endif

        // either after 8x SSE2 ASM, or standard trial division, record
        // how many factors we've found so far
        dconf->smooth_num[report_num] = smooth_num;

    }

    return;
}


void tdiv_medprimes_32k_new(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    // we have flagged this sieve offset as likely to produce a relation
    // nothing left to do now but check and see.
    uint32_t i;
    uint32_t prime, root1, root2, report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb_compressed* fbc;
    fb_element_siqs* fullfb_ptr, * fullfb = sconf->factor_base->list;
    uint32_t block_loc;

    fullfb_ptr = fullfb;
    if (parity)
    {
        fbc = dconf->comp_sieve_n;
    }
    else
    {
        fbc = dconf->comp_sieve_p;
    }

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

#if defined(USE_AVX512F)

        uint32_t base = block_loc - 32768;   // wraps mod 2^32; the offset below keeps t >= 0

        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            if ((i & 15) == 0) break;

            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];

            // off = block_loc - 32768 + 32769*p, a multiple-of-p shift that keeps t >= 0
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh
            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                // it will divide Q(x).  do so as many times as we can.
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        __m512i vbase = _mm512_set1_epi32((int)(block_loc - 32768));

        for (; i + 16 < sconf->factor_base->fb_13bit_B; i += 16)
        {
            __m512i vp = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) & fbc->prime[i]));
            __m512i vr1 = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) & fbc->root1[i]));
            __m512i vr2 = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) & fbc->root2[i]));
            __m512i vminv = _mm512_load_si512(&fullfb_ptr->minv32[i]);
            __m512i vthresh = _mm512_load_si512(&fullfb_ptr->thresh32[i]);

            // voff = block_loc - 32768 + 32769*p   (32769*p = (p << 15) + p)
            __m512i voff = _mm512_add_epi32(vbase,
                _mm512_add_epi32(_mm512_slli_epi32(vp, 15), vp));

            __m512i t1 = _mm512_sub_epi32(voff, vr1);
            __m512i t2 = _mm512_sub_epi32(voff, vr2);

            __mmask16 hit1 = _mm512_cmple_epu32_mask(_mm512_mullo_epi32(t1, vminv), vthresh);
            __mmask16 hit2 = _mm512_cmple_epu32_mask(_mm512_mullo_epi32(t2, vminv), vthresh);
            unsigned mask = (unsigned)(hit1 | hit2);

            while (mask)
            {
                int lane = _trail_zcnt(mask);
                mask &= mask - 1;

                DIVIDE_RESIEVED_PRIME_2(i + lane);
            }
        }

        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];

            // off = block_loc - 32768 + 32769*p, a multiple-of-p shift that keeps t >= 0
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh
            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                // it will divide Q(x).  do so as many times as we can.
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        dconf->smooth_num[report_num] = smooth_num;


#elif defined(USE_AVX2)

        __m256i vbase = _mm256_set1_epi32((int)(block_loc - 32768));

        uint32_t base = block_loc - 32768;   // wraps mod 2^32; the offset below keeps t >= 0

        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            if ((i & 7) == 0) break;

            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];

            // off = block_loc - 32768 + 32769*p, a multiple-of-p shift that keeps t >= 0
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh
            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                // it will divide Q(x).  do so as many times as we can.
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        for (; i + 8 < sconf->factor_base->fb_13bit_B; i += 8)
        {
            __m256i vp = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) & fbc->prime[i]));
            __m256i vr1 = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) & fbc->root1[i]));
            __m256i vr2 = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) & fbc->root2[i]));
            __m256i vminv = _mm256_load_si256((const __m256i*) & fullfb_ptr->minv32[i]);
            __m256i vthresh = _mm256_load_si256((const __m256i*) & fullfb_ptr->thresh32[i]);

            __m256i voff = _mm256_add_epi32(vbase,
                _mm256_add_epi32(_mm256_slli_epi32(vp, 15), vp));

            __m256i x1 = _mm256_mullo_epi32(_mm256_sub_epi32(voff, vr1), vminv);
            __m256i x2 = _mm256_mullo_epi32(_mm256_sub_epi32(voff, vr2), vminv);

            __m256i h1 = _mm256_cmpeq_epi32(_mm256_min_epu32(x1, vthresh), x1);
            __m256i h2 = _mm256_cmpeq_epi32(_mm256_min_epu32(x2, vthresh), x2);
            unsigned mask = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_or_si256(h1, h2)));

            while (mask) {
                int lane = _trail_zcnt(mask);
                mask &= mask - 1;
                DIVIDE_RESIEVED_PRIME_2(i + lane);
            }
        }

        // scalar tail loop for the remaining < 8 primes
        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];

            // off = block_loc - 32768 + 32769*p, a multiple-of-p shift that keeps t >= 0
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh
            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                // it will divide Q(x).  do so as many times as we can.
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        dconf->smooth_num[report_num] = smooth_num;

#elif defined(USE_SSE41) || defined(D_HAS_SSE2)

        uint32_t base = block_loc - 32768;   // wraps mod 2^32; the offset below keeps t >= 0

        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            if ((i & 7) == 0) break;

            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];

            // off = block_loc - 32768 + 32769*p, a multiple-of-p shift that keeps t >= 0
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh
            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                // it will divide Q(x).  do so as many times as we can.
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        __m128i vbase = _mm_set1_epi32((int)(block_loc - 32768));

        for (; i + 4 < sconf->factor_base->fb_13bit_B; i += 4)
        {
            // 4 x uint16 = 64 bits, zero-extended to 4 x uint32
            __m128i vp = tdm_load4_u16(&fbc->prime[i]);
            __m128i vr1 = tdm_load4_u16(&fbc->root1[i]);
            __m128i vr2 = tdm_load4_u16(&fbc->root2[i]);
            __m128i vminv = _mm_loadu_si128((const __m128i*) & fullfb_ptr->minv32[i]);
            __m128i vthresh = _mm_loadu_si128((const __m128i*) & fullfb_ptr->thresh32[i]);

            // voff = block_loc - 32768 + 32769*p
            __m128i voff = _mm_add_epi32(vbase,
                _mm_add_epi32(_mm_slli_epi32(vp, 15), vp));

            __m128i x1 = tdm_mullo_epi32(_mm_sub_epi32(voff, vr1), vminv);
            __m128i x2 = tdm_mullo_epi32(_mm_sub_epi32(voff, vr2), vminv);

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh (unsigned)
            __m128i h1 = tdm_cmple_epu32(x1, vthresh);
            __m128i h2 = tdm_cmple_epu32(x2, vthresh);
            unsigned mask = _mm_movemask_ps(_mm_castsi128_ps(_mm_or_si128(h1, h2)));

            while (mask) {
                int lane = _trail_zcnt(mask);
                mask &= mask - 1;
                DIVIDE_RESIEVED_PRIME_2(i + lane);
            }
        }

        // scalar tail loop for the remaining < 4 primes
        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];

            // off = block_loc - 32768 + 32769*p, a multiple-of-p shift that keeps t >= 0
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh
            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                // it will divide Q(x).  do so as many times as we can.
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        dconf->smooth_num[report_num] = smooth_num;

#else

        uint32_t base = block_loc - 32768;   // wraps mod 2^32; the offset below keeps t >= 0

        for (; i < sconf->factor_base->fb_13bit_B; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];

            // off = block_loc - 32768 + 32769*p, a multiple-of-p shift that keeps t >= 0
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            // x divisible by p  <=>  x * minv (mod 2^32) <= thresh
            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                // it will divide Q(x).  do so as many times as we can.
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        dconf->smooth_num[report_num] = smooth_num;


#endif

    }

    return;
}


#if defined(__aarch64__) && defined(__ARM_NEON)

// AArch64 NEON version of tdiv_medprimes_32k_new: the same Lemire test,
// x divisible by p  <=>  x * minv (mod 2^32) <= thresh, on 8 primes per
// iteration (two uint32x4 halves, so the two vmulq chains overlap).  NEON
// has unsigned compares and no alignment rules, so no head loop is needed.
// Primes are visited in ascending order, as in the x86 kernels, so
// fb_offsets come out identical.  Selected under HAS_NEON (-noneon reverts
// to tdiv_medprimes_32k_new's portable loop).
void tdiv_medprimes_32k_neon(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    uint32_t i;
    uint32_t prime, root1, root2, report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb_compressed* fbc;
    fb_element_siqs* fullfb_ptr = sconf->factor_base->list;
    uint32_t block_loc;
    uint32_t bound = sconf->factor_base->fb_13bit_B;
    static const uint8_t lanebit[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };
    const uint8x8_t vlanebit = vld1_u8(lanebit);

    if (parity)
    {
        fbc = dconf->comp_sieve_n;
    }
    else
    {
        fbc = dconf->comp_sieve_p;
    }

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {
        if (!dconf->valid_Qs[report_num])
            continue;

        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        // wraps mod 2^32; off = base + 32769*p keeps t >= 0 (see _new)
        uint32_t base = block_loc - 32768;
        uint32x4_t vbase = vdupq_n_u32(base);

        i = sconf->sieve_small_fb_start;

        for (; i + 8 <= bound; i += 8)
        {
            uint16x8_t p16 = vld1q_u16(&fbc->prime[i]);
            uint16x8_t r1_16 = vld1q_u16(&fbc->root1[i]);
            uint16x8_t r2_16 = vld1q_u16(&fbc->root2[i]);

            uint32x4_t plo = vmovl_u16(vget_low_u16(p16));
            uint32x4_t phi = vmovl_high_u16(p16);

            // voff = block_loc - 32768 + 32769*p   (32769*p = (p << 15) + p)
            uint32x4_t offlo = vaddq_u32(vbase, vsliq_n_u32(plo, plo, 15));
            uint32x4_t offhi = vaddq_u32(vbase, vsliq_n_u32(phi, phi, 15));

            uint32x4_t minvlo = vld1q_u32(&fullfb_ptr->minv32[i]);
            uint32x4_t minvhi = vld1q_u32(&fullfb_ptr->minv32[i + 4]);
            uint32x4_t thrlo = vld1q_u32(&fullfb_ptr->thresh32[i]);
            uint32x4_t thrhi = vld1q_u32(&fullfb_ptr->thresh32[i + 4]);

            uint32x4_t x1lo = vmulq_u32(vsubq_u32(offlo, vmovl_u16(vget_low_u16(r1_16))), minvlo);
            uint32x4_t x1hi = vmulq_u32(vsubq_u32(offhi, vmovl_high_u16(r1_16)), minvhi);
            uint32x4_t x2lo = vmulq_u32(vsubq_u32(offlo, vmovl_u16(vget_low_u16(r2_16))), minvlo);
            uint32x4_t x2hi = vmulq_u32(vsubq_u32(offhi, vmovl_high_u16(r2_16)), minvhi);

            uint32x4_t hlo = vorrq_u32(vcleq_u32(x1lo, thrlo), vcleq_u32(x2lo, thrlo));
            uint32x4_t hhi = vorrq_u32(vcleq_u32(x1hi, thrhi), vcleq_u32(x2hi, thrhi));

            // narrow the 8 lane masks to 8 bytes (0x00 / 0xff); hits are
            // rare, so test the 64-bit word before building a lane mask
            uint8x8_t h8 = vmovn_u16(vcombine_u16(vmovn_u32(hlo), vmovn_u32(hhi)));
            if (vget_lane_u64(vreinterpret_u64_u8(h8), 0) == 0)
                continue;

            uint32_t mask = vaddv_u8(vand_u8(h8, vlanebit));
            while (mask)
            {
                uint32_t lane = (uint32_t)__builtin_ctz(mask);
                mask &= mask - 1;
                DIVIDE_RESIEVED_PRIME_2(i + lane);
            }
        }

        // scalar tail loop for the remaining < 8 primes
        for (; i < bound; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            uint32_t minv = fullfb_ptr->minv32[i];
            uint32_t thresh = fullfb_ptr->thresh32[i];
            uint32_t off = base + ((uint32_t)prime << 15) + prime;

            if ((uint32_t)((off - root1) * minv) <= thresh ||
                (uint32_t)((off - root2) * minv) <= thresh)
            {
                DIVIDE_RESIEVED_PRIME_2(i);
            }
        }

        dconf->smooth_num[report_num] = smooth_num;
    }

    return;
}

#endif
