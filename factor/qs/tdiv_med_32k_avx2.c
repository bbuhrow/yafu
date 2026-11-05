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


#if defined( USE_AVX2 )

#include "tdiv_macros_common.h"
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

this file contains code implementing 3)


*/


void tdiv_medprimes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    // we have flagged this sieve offset as likely to produce a relation
    // nothing left to do now but check and see.
    uint32_t i;
    uint32_t prime, root1, root2, report_num;
    int smooth_num;
    uint32_t *fb_offsets;
    sieve_fb_compressed *fbc;
    fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
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

        CLEAN_AVX2;

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
                DIVIDE_ONE_PRIME;
            }
        }

        __m512i vbase = _mm512_set1_epi32((int)(block_loc - 32768));

        for (; i + 16 <= sconf->factor_base->fb_13bit_B; i += 16)
        {
            __m512i vp = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) &fbc->prime[i]));
            __m512i vr1 = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) &fbc->root1[i]));
            __m512i vr2 = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i*) &fbc->root2[i]));
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
                DIVIDE_ONE_PRIME;
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
                DIVIDE_ONE_PRIME;
            }
        }

        for (; i + 8 <= sconf->factor_base->fb_13bit_B; i += 8)
        {
            __m256i vp = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) &fbc->prime[i]));
            __m256i vr1 = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) &fbc->root1[i]));
            __m256i vr2 = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*) &fbc->root2[i]));
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
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;

#elif defined(USE_SSE41)

        __m128i vbase = _mm_set1_epi32((int)(block_loc - 32768));

        for (; i + 4 <= sconf->factor_base->fb_13bit_B; i += 4)
        {
            // 4 x uint16 = 64 bits, zero-extended to 4 x uint32
            __m128i vp = _mm_cvtepu16_epi32(_mm_loadu_si128((const __m128i*) &fbc->prime[i]));
            __m128i vr1 = _mm_cvtepu16_epi32(_mm_loadu_si128((const __m128i*) &fbc->root1[i]));
            __m128i vr2 = _mm_cvtepu16_epi32(_mm_loadu_si128((const __m128i*) &fbc->root2[i]));
            __m128i vminv = _mm_loadu_si128((const __m128i*) &fullfb_ptr->minv32[i]);
            __m128i vthresh = _mm_loadu_si128((const __m128i*) &fullfb_ptr->thresh32[i]);

            // voff = block_loc - 32768 + 32769*p
            __m128i voff = _mm_add_epi32(vbase,
                _mm_add_epi32(_mm_slli_epi32(vp, 15), vp));

            __m128i x1 = _mm_mullo_epi32(_mm_sub_epi32(voff, vr1), vminv);
            __m128i x2 = _mm_mullo_epi32(_mm_sub_epi32(voff, vr2), vminv);

            // unsigned x <= thresh  <=>  min_epu32(x, thresh) == x
            __m128i h1 = _mm_cmpeq_epi32(_mm_min_epu32(x1, vthresh), x1);
            __m128i h2 = _mm_cmpeq_epi32(_mm_min_epu32(x2, vthresh), x2);
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
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;

#else

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
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;


#endif

    }

    return;
}

#endif // USE_AVX2
