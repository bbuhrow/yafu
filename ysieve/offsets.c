/*
MIT License

Copyright (c) 2021 Ben Buhrow

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "soe.h"
#if defined(_MSC_VER) && defined(__clang__)
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif
#include "ytools.h"
#include <math.h>


#define NO_64U_REM
//#define U64_REM_ONLY

#ifdef USE_AVX2
// macros for Montgomery arithmetic - helpful for computing 
// division-less offsets once we have enough reuse (number of
// classes) to justify the setup costs.
#define _mm256_cmpge_epu32(a, b) \
        _mm256_cmpeq_epi32(_mm256_max_epu32(a, b), a)

#define _mm256_cmple_epu32(a, b) _mm256_cmpge_epu32(b, a)


static __inline __m256i CLEAR_HIGH_VEC(__m256i x)
{
    __m256i chi = _mm256_set1_epi64x(0x00000000ffffffff);
    return _mm256_and_si256(chi, x);
}

static __inline __m256i vec_redc(__m256i x64e, __m256i x64o, __m256i pinv, __m256i p)
{
    // uint32_t m = (uint32_t)x * pinv;
    __m256i t1 = _mm256_shuffle_epi32(pinv, 0xB1);      // odd-index pinv in lo words
    __m256i even = _mm256_mul_epu32(x64e, pinv);
    __m256i odd = _mm256_mul_epu32(x64o, t1);
    __m256i t2;

    // x += (uint64_t)m * (uint64_t)p;
    t1 = _mm256_shuffle_epi32(p, 0xB1);      // odd-index p in lo words
    even = _mm256_add_epi64(x64e, _mm256_mul_epu32(even, p));
    odd = _mm256_add_epi64(x64o, _mm256_mul_epu32(odd, t1));

    // m = x >> 32;
    t1 = _mm256_blend_epi32(odd, _mm256_shuffle_epi32(even, 0xB1), 0x55);

    // if (m >= p) m -= p;
    t2 = _mm256_cmpge_epu32(t1, p); //_mm256_or_si256(_mm256_cmpgt_epi32(t1, p), _mm256_cmpeq_epi32(t1, p));
    t2 = _mm256_and_si256(p, t2);

    return _mm256_sub_epi32(t1, t2);
}

static __inline __m256i vec_to_monty(__m256i x, __m256i r2, __m256i pinv, __m256i p)
{
    //uint64_t t = (uint64_t)x * (uint64_t)r2;
    __m256i t1 = _mm256_shuffle_epi32(x, 0xB1);
    __m256i t2 = _mm256_shuffle_epi32(r2, 0xB1);
    __m256i even = _mm256_mul_epu32(x, r2);
    __m256i odd = _mm256_mul_epu32(t1, t2);

    return vec_redc(even, odd, pinv, p);
}

#endif

void get_offsets(thread_soedata_t *thread_data)
{
    //extract stuff from the thread data structure
    soe_dynamicdata_t *ddata = &thread_data->ddata;
    soe_staticdata_t *sdata = &thread_data->sdata;

    uint64_t i, startprime = sdata->startprime, prodN = sdata->prodN, block = 0;
    uint32_t prime, root, bnum;
    // implement a special case for lowlimit = 0?  
    // Then lmp[i] + diff = rclass b/c lmp = 1.
    uint32_t diff = sdata->rclass[thread_data->current_line] - 1;
    uint64_t tmp2;
    int s;
    int FLAGSIZE = sdata->FLAGSIZE;
    int FLAGBITS = sdata->FLAGBITS;

    // failsafe: set all blocks to sieve with all primes.  the loop below will overwrite
    // these with better limits according to the size of flags in the blocks.
    ddata->largep_offset = 0;

    for (i = 0; i < sdata->blocks; i++)
    {
        ddata->pbounds[i] = sdata->bucket_start_id;

        //initialize bucket
        if (ddata->bucket_depth > 0)
            ddata->bucket_hits[i] = 0;
    }

    if (sdata->sieve_range == 0)
    {
        //uint32_t *lmp = sdata->lower_mod_prime;

        mpz_t gmp_sqrt;
        mpz_init(gmp_sqrt);

        for (i = startprime; i < sdata->bucket_start_id; i++)
        {
            prime = sdata->sieve_p[i];

            // find the first multiple of the prime which is greater than the first sieve location 
            // and also equal to the residue class mod 'prodN'.  
            // we need to solve the congruence: rclass[current_line] == kp mod prodN for k
            // xGCD gives r and s such that r*p + s*prodN = gcd(p,prodN).
            // then k = r*class/gcd(p,prodN) is a solution.
            // the gcd of p and prodN is always 1 by construction of prodN and choice of p.  
            // therefore k = r * class is a solution.  furthermore, since the gcd is 1, there
            // is only one solution.  
            // xGCD_1((int)prime,(int)prodN,&r,&s,&tmp);

            // To speed things up we solve and store modinv(prodN, prime) for every prime (only
            // needs to be done once, in roots.c).  Then to get the offset for the current block
            // we just need to multiply the stored root with the starting sieve location (mod p).	

            // if the prime is greater than the limit at which it is necessary to sieve
            // a block, start that prime in the next block.
            if (sdata->sieve_p[i] > ddata->blk_b_sqrt)
            {
                //printf("lblk_b = %lu, blk_b_sqrt = %lu, pbounds block %lu = %lu (%u)\n", 
                //    ddata->lblk_b, ddata->blk_b_sqrt, block, i, sdata->sieve_p[i]);
                ddata->pbounds[block] = i;

                // if (block < (sdata->blocks - 0)) // <-- sometimes crashes, but correct ranges
                // if (block < (sdata->blocks - 1)) // <-- no crashes, but incorrect ranges
                if (block < (sdata->blocks - 1))
                    block++;
                ddata->lblk_b = ddata->ublk_b + prodN;
                ddata->ublk_b += sdata->blk_r;
                //ddata->blk_b_sqrt = (uint64_t)(sqrt((int64_t)(ddata->ublk_b + prodN))) + 1;
                mpz_set_ui(gmp_sqrt, ddata->ublk_b + prodN);
                mpz_sqrt(gmp_sqrt, gmp_sqrt);
                ddata->blk_b_sqrt = mpz_get_ui(gmp_sqrt) + 1;
            }

            s = sdata->root[i];

            // the lower block bound (lblk_b) times s can exceed 64 bits for large ranges,
            // so reduce mod p here as well.
            tmp2 = (uint64_t)s * (ddata->lblk_b % (uint64_t)prime);

            // tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);
            ddata->offsets[i] = (uint32_t)(tmp2 % (uint64_t)prime);
        }

        mpz_clear(gmp_sqrt);
    }
    else
    {
        uint32_t modp;
        mpz_t lowz, sqrtz;
        mpz_init(lowz);
        mpz_init(sqrtz);
        mpz_set(lowz, sdata->offset);
        mpz_add_ui(lowz, lowz, ddata->lblk_b);

        mpz_set(sqrtz, lowz);
        mpz_add_ui(sqrtz, sqrtz, sdata->blk_r);
        mpz_sqrt(sqrtz, sqrtz);
        mpz_add_ui(sqrtz, sqrtz, 1);

        // if we're sieving with an offset, use all of the primes for each block
        // and just find the offset into the first block
        for (i = startprime; i < sdata->bucket_start_id; i++)
        {
            prime = sdata->sieve_p[i];
            s = sdata->root[i];

            if (mpz_cmp_ui(sqrtz, sdata->sieve_p[i]) <= 0)
            {
                ddata->pbounds[block] = i;
                block++;
                mpz_add_ui(lowz, lowz, sdata->blk_r);
                mpz_set(sqrtz, lowz);
                mpz_add_ui(sqrtz, sqrtz, sdata->blk_r);
                mpz_sqrt(sqrtz, sqrtz);
                mpz_add_ui(sqrtz, sqrtz, 1);
            }

            modp = mpz_tdiv_ui(lowz, prime);
            tmp2 = (uint64_t)s * (uint64_t)modp;
            ddata->offsets[i] = (uint32_t)(tmp2 % (uint64_t)prime);
            //gmp_printf("p = %u, o = %u, r = %d, lblk_b = %Zd, modp = %u\n", 
            //	prime, ddata->offsets[i], s, tmpz, modp);
        }

        mpz_clear(lowz);
        mpz_clear(sqrtz);
    }

    if (ddata->bucket_depth > 0)
    {
        uint64_t** bptr;

        uint32_t* nptr;
        uint32_t linesize = FLAGSIZE * sdata->blocks;
        uint32_t* lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;

        nptr = ddata->bucket_hits;
        bptr = ddata->sieve_buckets;

#if defined(USE_AVX2)

        if (sdata->use_monty)
        {
            // AVX2 (in soe_util.c) has arranged things so that the end of the
            // sieve prime array is at an index divisible by 8.
            for (; i < sdata->bitmap_start_id; i += 8)
            {
                __m256i vp = _mm256_loadu_si256((__m256i*)(&sdata->sieve_p[i]));
                __m256i vpinv = _mm256_loadu_si256((__m256i*)(&sdata->pinv[i]));
                __m256i vr2 = _mm256_loadu_si256((__m256i*)(&sdata->r2modp[i]));
                __m256i vr = _mm256_loadu_si256((__m256i*)(&sdata->root[i]));
                __m256i vlmp = _mm256_loadu_si256((__m256i*)(&lmp[i]));
                __m256i vdiff = _mm256_set1_epi32(diff);
                __m256i t1, t2, t3;
                __m256i even, odd;
                ALIGNED_MEM uint64_t tmp[8];

                // condition to see if the current prime only hits the sieve interval once
                if (sdata->sieve_p[i] > sdata->large_bucket_start_prime)
                {
                    ddata->largep_offset = i;
                    break;
                }

                // lmp[i:i+7] + diff
                vlmp = _mm256_add_epi32(vlmp, vdiff);

                // lmp to monty     
                t1 = vec_to_monty(vlmp, vr2, vpinv, vp);

                //tmp2 = (uint64_t)s * (uint64_t)t;
                t3 = _mm256_shuffle_epi32(t1, 0xB1);
                t2 = _mm256_shuffle_epi32(vr, 0xB1);

                even = _mm256_mul_epu32(vr, t1);
                odd = _mm256_mul_epu32(t3, t2);

                // reduce
                vr = vec_redc(even, odd, vpinv, vp);

                // take out of monty rep
                vr = vec_redc(CLEAR_HIGH_VEC(vr), 
                    CLEAR_HIGH_VEC(_mm256_shuffle_epi32(vr, 0xB1)), vpinv, vp);

                //t1 = _mm256_set1_epi32(linesize);
                //t1 = _mm256_or_si256(_mm256_cmpgt_epi32(t1, vr), _mm256_cmpeq_epi32(t1, vr));
                //mask = ~_mm256_movemask_epi8(t1);
                //t1 = _mm256_srai_epi32(vr, FLAGBITS);
                //t2 = _mm256_cmpgt_epi32(t1, _mm256_set1_epi32(sdata->blocks-1));
                //mask = ~_mm256_movemask_epi8(t2);
                //_mm256_store_si256((__m256i *)tmp32, t1);

                // recombine and store for distribution to the buckets
                t2 = _mm256_blend_epi32(_mm256_shuffle_epi32(vp, 0xB1), vr, 0x55);
                t3 = _mm256_blend_epi32(vp, _mm256_shuffle_epi32(vr, 0xB1), 0x55);

                _mm256_store_si256((__m256i*)tmp, t2);         // even (prime | root)
                _mm256_store_si256((__m256i*)(&tmp[4]), t3);   // odd  (prime | root)

                // bucket sort the roots
                root = (uint32_t)tmp[0];
                if (root < linesize)
                    //if (mask & 0x1)
                {
                    bnum = ((uint32_t)tmp[0] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[0];
                    nptr[bnum]++;
                }

                root = (uint32_t)tmp[1];
                if (root < linesize)
                    //if (mask & 0x10)
                {
                    bnum = ((uint32_t)tmp[1] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[1];
                    nptr[bnum]++;
                }

                root = (uint32_t)tmp[2];
                if (root < linesize)
                    //if (mask & 0x100)
                {
                    bnum = ((uint32_t)tmp[2] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[2];
                    nptr[bnum]++;
                }

                root = (uint32_t)tmp[3];
                if (root < linesize)
                    //if (mask & 0x1000)
                {
                    bnum = ((uint32_t)tmp[3] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[3];
                    nptr[bnum]++;
                }

                root = (uint32_t)tmp[4];
                if (root < linesize)
                    //if (mask & 0x10000)
                {
                    bnum = ((uint32_t)tmp[4] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[4];
                    nptr[bnum]++;
                }

                root = (uint32_t)tmp[5];
                if (root < linesize)
                    //if (mask & 0x100000)
                {
                    bnum = ((uint32_t)tmp[5] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[5];
                    nptr[bnum]++;
                }

                root = (uint32_t)tmp[6];
                if (root < linesize)
                    //if (mask & 0x1000000)
                {
                    bnum = ((uint32_t)tmp[6] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[6];
                    nptr[bnum]++;
                }

                root = (uint32_t)tmp[7];
                if (root < linesize)
                    //if (mask & 0x10000000)
                {
                    bnum = ((uint32_t)tmp[7] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[7];
                    nptr[bnum]++;
                }
            }
        }
        else
        {
            // AVX2 (in soe_util.c) has arranged things so that the end of the
            // sieve prime array is at an index divisible by 8.
            for (; i < sdata->bitmap_start_id; i += 2)
            {
                uint64_t tmp3;
                uint32_t p2, r2;
                int s2;

                prime = sdata->sieve_p[i];
                p2 = sdata->sieve_p[i + 1];

                // condition to see if the current prime only hits the sieve interval once
                if (prime > sdata->large_bucket_start_prime)
                {
                    ddata->largep_offset = i;
                    break;
                }

                s = sdata->root[i];
                s2 = sdata->root[i + 1];

                // we solved for lower_mod_prime while computing the modular inverse of
                // each prime, for the residue class 1.  add the difference between this
                // residue class and 1 before multiplying by the modular inverse to find the offset.

                tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);
                tmp3 = (uint64_t)s2 * (uint64_t)(lmp[i + 1] + diff);

                root = (uint32_t)(tmp2 % (uint64_t)prime);
                r2 = (uint32_t)(tmp3 % (uint64_t)p2);

                // It is faster to update during
                // linesieve than doing it all here in a loop.
                // measured 6/2016
                if (root < linesize)
                {
                    bnum = (root >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = ((uint64_t)prime << 32) | (uint64_t)root;
                    nptr[bnum]++;
                }

                if (r2 < linesize)
                {
                    bnum = (r2 >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = ((uint64_t)p2 << 32) | (uint64_t)r2;
                    nptr[bnum]++;
                }

            }

            for (; i < sdata->bitmap_start_id; i++)
            {
                prime = sdata->sieve_p[i];

                // condition to see if the current prime only hits the sieve interval once
                if (prime > sdata->large_bucket_start_prime)
                {
                    ddata->largep_offset = i;
                    break;
                }

                s = sdata->root[i];

                // we solved for lower_mod_prime while computing the modular inverse of
                // each prime, for the residue class 1.  add the difference between this
                // residue class and 1 before multiplying by the modular inverse to find the offset.

                tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);
                root = (uint32_t)(tmp2 % (uint64_t)prime);

                // It is faster to update during
                // linesieve than doing it all here in a loop.
                // measured 6/2016
                if (root < linesize)
                {
                    bnum = (root >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = ((uint64_t)prime << 32) | (uint64_t)root;
                    nptr[bnum]++;
                }
            }
        }

#else

        for (; i < sdata->bitmap_start_id - 1; i += 2)
        {
            uint64_t tmp3;
            uint32_t p2, r2;
            int s2;

            prime = sdata->sieve_p[i];
            p2 = sdata->sieve_p[i + 1];

            // condition to see if the current prime only hits the sieve interval once
            if (prime > sdata->large_bucket_start_prime)
            {
                ddata->largep_offset = i;
                break;
            }

            s = sdata->root[i];
            s2 = sdata->root[i + 1];

            // we solved for lower_mod_prime while computing the modular inverse of
            // each prime, for the residue class 1.  add the difference between this
            // residue class and 1 before multiplying by the modular inverse to find the offset.
            tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);
            tmp3 = (uint64_t)s2 * (uint64_t)(lmp[i + 1] + diff);

            root = (uint32_t)(tmp2 % (uint64_t)prime);
            r2 = (uint32_t)(tmp3 % (uint64_t)p2);

            // It is faster to update during
            // linesieve than doing it all here in a loop.
            // measured 6/2016
            if (root < linesize)
            {
                bnum = (root >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64_t)prime << 32) | (uint64_t)root;
                nptr[bnum]++;
            }

            if (r2 < linesize)
            {
                bnum = (r2 >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64_t)p2 << 32) | (uint64_t)r2;
                nptr[bnum]++;
            }

        }

        if ((i < sdata->bitmap_start_id) && (ddata->largep_offset == 0))
        {
            uint32_t* lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;
            prime = sdata->sieve_p[i];

            s = sdata->root[i];

            tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);
            root = (uint32_t)(tmp2 % (uint64_t)prime);

            nptr = ddata->bucket_hits;
            bptr = ddata->sieve_buckets;

            if (root < linesize)
            {
                bnum = (root >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64_t)prime << 32) | (uint64_t)root;
                nptr[bnum]++;
            }

        }
#endif



        if (ddata->largep_offset > 0)
        {
            // primes greater than the entire sieve interval, thus they
            // at most hit one block and we don't need to save the prime
            // itself since it doesn't need to be advanced.
            uint32_t** large_bptr;
            uint32_t* large_nptr;
            uint32_t* lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;

            large_nptr = ddata->large_bucket_hits;
            large_bptr = ddata->large_sieve_buckets;

            for (i = 0; i < sdata->blocks; i++)
            {
                //initialize bucket
                large_nptr[i] = 0;
            }

#if defined(USE_AVX2)

            if ((sdata->has_avx2) && (sdata->use_monty))
            {
                // AVX2 (in soe_util.c) has arranged things so that the end of the
                // sieve prime array is at an index divisible by 8.
                for (i = ddata->largep_offset; i < sdata->bitmap_start_id; i += 8)
                {
                    __m256i vp = _mm256_loadu_si256((__m256i*)(&sdata->sieve_p[i]));
                    __m256i vpinv = _mm256_loadu_si256((__m256i*)(&sdata->pinv[i]));
                    __m256i vr2 = _mm256_loadu_si256((__m256i*)(&sdata->r2modp[i]));
                    __m256i vr = _mm256_loadu_si256((__m256i*)(&sdata->root[i]));
                    __m256i vlmp = _mm256_loadu_si256((__m256i*)(&lmp[i]));
                    __m256i vdiff = _mm256_set1_epi32(diff);
                    __m256i t1, t2, t3;
                    __m256i even, odd;
                    ALIGNED_MEM uint32_t tmp[8];
                    uint32_t mask;

                    // lmp[i:i+7] + diff
                    vlmp = _mm256_add_epi32(vlmp, vdiff);

                    // lmp to monty     
                    t1 = vec_to_monty(vlmp, vr2, vpinv, vp);

                    //tmp2 = (uint64_t)s * (uint64_t)t;
                    t3 = _mm256_shuffle_epi32(t1, 0xB1);
                    t2 = _mm256_shuffle_epi32(vr, 0xB1);

                    even = _mm256_mul_epu32(vr, t1);
                    odd = _mm256_mul_epu32(t3, t2);

                    // reduce
                    vr = vec_redc(even, odd, vpinv, vp);

                    // take out of monty rep
                    vr = vec_redc(CLEAR_HIGH_VEC(vr), CLEAR_HIGH_VEC(_mm256_shuffle_epi32(vr, 0xB1)), vpinv, vp);

                    t1 = _mm256_set1_epi32(linesize);
                    t1 = _mm256_cmpge_epu32(vr, t1);
                    mask = ~(_mm256_movemask_epi8(t1));

                    _mm256_storeu_si256((__m256i*)tmp, vr);

                    if (mask & 0x1)
                    {
                        bnum = (tmp[0] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[0];
                        large_nptr[bnum]++;
                    }

                    if (mask & 0x10)
                    {
                        bnum = (tmp[1] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[1];
                        large_nptr[bnum]++;
                    }

                    if (mask & 0x100)
                    {
                        bnum = (tmp[2] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[2];
                        large_nptr[bnum]++;
                    }

                    if (mask & 0x1000)
                    {
                        bnum = (tmp[3] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[3];
                        large_nptr[bnum]++;
                    }

                    if (mask & 0x10000)
                    {
                        bnum = (tmp[4] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[4];
                        large_nptr[bnum]++;
                    }

                    if (mask & 0x100000)
                    {
                        bnum = (tmp[5] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[5];
                        large_nptr[bnum]++;
                    }

                    if (mask & 0x1000000)
                    {
                        bnum = (tmp[6] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[6];
                        large_nptr[bnum]++;
                    }

                    if (mask & 0x10000000)
                    {
                        bnum = (tmp[7] >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = tmp[7];
                        large_nptr[bnum]++;
                    }
                }
            }
            else
            {
                // AVX2 (in soe_util.c) has arranged things so that the end of the
                // sieve prime array is at an index divisible by 8.
                for (i = ddata->largep_offset; i < sdata->bitmap_start_id; i++)
                {
                    prime = sdata->sieve_p[i];
                    s = sdata->root[i];

                    // we solved for lower_mod_prime while computing the modular inverse of
                    // each prime, for the residue class 1.  add the difference between this
                    // residue class and 1 before multiplying by the modular inverse.
                    // could use (_mm256_mul_epu32 --> VPMULUDQ)
                    tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);

                    // would need custom solution
                    root = (uint32_t)(tmp2 % (uint64_t)prime);

                    // gather may help, but writes would need to be done 1 by 1.
                    if (root < linesize)
                    {
                        bnum = (root >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = root;
                        large_nptr[bnum]++;
                    }
                }
            }


#else

            for (i = ddata->largep_offset; i < sdata->bitmap_start_id - 1; i += 2)
            {
                uint64_t tmp3;
                uint32_t p2, r2;
                int s2;

                prime = sdata->sieve_p[i];
                p2 = sdata->sieve_p[i + 1];

                s = sdata->root[i];
                s2 = sdata->root[i + 1];

                // we solved for lower_mod_prime while computing the modular inverse of
                // each prime, for the residue class 1.  add the difference between this
                // residue class and 1 before multiplying by the modular inverse.
                // could use (_mm256_mul_epu32 --> VPMULUDQ)
                tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);
                tmp3 = (uint64_t)s2 * (uint64_t)(lmp[i + 1] + diff);

                // would need custom solution
                root = (uint32_t)(tmp2 % (uint64_t)prime);
                r2 = (uint32_t)(tmp3 % (uint64_t)p2);

                // gather may help, but writes would need to be done 1 by 1.
                if (root < linesize)
                {
                    bnum = (root >> FLAGBITS);
                    large_bptr[bnum][large_nptr[bnum]] = root;
                    large_nptr[bnum]++;
                }

                if (r2 < linesize)
                {
                    bnum = (r2 >> FLAGBITS);
                    large_bptr[bnum][large_nptr[bnum]] = r2;
                    large_nptr[bnum]++;
                }

            }

            if (i < sdata->bitmap_start_id)
            {
                prime = sdata->sieve_p[i];
                s = sdata->root[i];

                tmp2 = (uint64_t)s * (uint64_t)(lmp[i] + diff);
                root = (uint32_t)(tmp2 % (uint64_t)prime);

                if (root < linesize)
                {
                    bnum = (root >> FLAGBITS);
                    large_bptr[bnum][large_nptr[bnum]] = root;
                    large_nptr[bnum]++;
                }

            }

#endif

        }

    }

	return;
}
