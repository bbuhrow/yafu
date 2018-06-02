/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 7/1/10
----------------------------------------------------------------------*/

#include "soe.h"
#include "immintrin.h"


#define NO_64U_REM
//#define U64_REM_ONLY

void get_offsets(thread_soedata_t *thread_data)
{
    //extract stuff from the thread data structure
    soe_dynamicdata_t *ddata = &thread_data->ddata;
    soe_staticdata_t *sdata = &thread_data->sdata;

    uint64 i, startprime = sdata->startprime, prodN = sdata->prodN, block = 0;
    uint32 prime, root, bnum;
    uint32 diff = sdata->rclass[thread_data->current_line] - 1;
    uint64 tmp2;
    int s;

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
        uint32 *lmp = sdata->lower_mod_prime;

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
                ddata->pbounds[block] = i;
                block++;
                ddata->lblk_b = ddata->ublk_b + prodN;
                ddata->ublk_b += sdata->blk_r;
                ddata->blk_b_sqrt = (uint64)(sqrt((int64)(ddata->ublk_b + prodN))) + 1;
            }

            s = sdata->root[i];

            // the lower block bound (lblk_b) times s can exceed 64 bits for large ranges,
            // so reduce mod p here as well.
            tmp2 = (uint64)s * (ddata->lblk_b % (uint64)prime);

            // tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
            ddata->offsets[i] = (uint32)(tmp2 % (uint64)prime);
        }
    }
    else
    {
        uint32 modp;
        mpz_t lowz, sqrtz;
        mpz_init(lowz);
        mpz_init(sqrtz);
        mpz_set(lowz, *sdata->offset);
        mpz_add_ui(lowz, lowz, ddata->lblk_b);

        mpz_set(sqrtz, lowz);
        mpz_add_ui(sqrtz, sqrtz, sdata->blk_r);
        mpz_sqrt(sqrtz, sqrtz);
        mpz_add_ui(sqrtz, sqrtz, 1);
        //mpz_set_ui(tmpz, ddata->lblk_b);

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
            tmp2 = (uint64)s * (uint64)modp;
            ddata->offsets[i] = (uint32)(tmp2 % (uint64)prime);
            //gmp_printf("p = %u, o = %u, r = %d, lblk_b = %Zd, modp = %u\n", 
            //	prime, ddata->offsets[i], s, tmpz, modp);
        }

        mpz_clear(lowz);
        mpz_clear(sqrtz);
    }

    if (ddata->bucket_depth > 0)
    {
        uint64 **bptr;

        uint32 *nptr;
        uint32 linesize = FLAGSIZE * sdata->blocks;
        uint32 *lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;

        nptr = ddata->bucket_hits;
        bptr = ddata->sieve_buckets;

#if defined(USE_AVX2)

        if (sdata->use_monty)
        {
            // AVX2 (in soe_util.c) has arranged things so that the end of the
            // sieve prime array is at an index divisible by 8.
            for (; i < sdata->bitmap_start_id; i += 8)
            {
                __m256i vp = _mm256_loadu_si256((__m256i *)(&sdata->sieve_p[i]));
                __m256i vpinv = _mm256_loadu_si256((__m256i *)(&sdata->pinv[i]));
                __m256i vr2 = _mm256_loadu_si256((__m256i *)(&sdata->r2modp[i]));
                __m256i vr = _mm256_loadu_si256((__m256i *)(&sdata->root[i]));
                __m256i vlmp = _mm256_loadu_si256((__m256i *)(&lmp[i]));
                __m256i vdiff = _mm256_set1_epi32(diff);
                __m256i t1, t2, t3;
                __m256i even, odd;
                ALIGNED_MEM uint64 tmp[8];
                ALIGNED_MEM uint32 tmp32[8];
                uint32 mask;

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

                //tmp2 = (uint64)s * (uint64)t;
                t3 = _mm256_shuffle_epi32(t1, 0xB1);
                t2 = _mm256_shuffle_epi32(vr, 0xB1);

                even = _mm256_mul_epu32(vr, t1);
                odd = _mm256_mul_epu32(t3, t2);

                // reduce
                vr = vec_redc(even, odd, vpinv, vp);

                // take out of monty rep
                vr = vec_redc(CLEAR_HIGH_VEC(vr), CLEAR_HIGH_VEC(_mm256_shuffle_epi32(vr, 0xB1)), vpinv, vp);

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

                _mm256_store_si256((__m256i *)tmp, t2);         // even (prime | root)
                _mm256_store_si256((__m256i *)(&tmp[4]), t3);   // odd  (prime | root)

                // bucket sort the roots
                root = (uint32)tmp[0];
                if (root < linesize)
                    //if (mask & 0x1)
                {
                    bnum = ((uint32)tmp[0] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[0];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[1];
                if (root < linesize)
                    //if (mask & 0x10)
                {
                    bnum = ((uint32)tmp[1] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[1];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[2];
                if (root < linesize)
                    //if (mask & 0x100)
                {
                    bnum = ((uint32)tmp[2] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[2];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[3];
                if (root < linesize)
                    //if (mask & 0x1000)
                {
                    bnum = ((uint32)tmp[3] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[3];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[4];
                if (root < linesize)
                    //if (mask & 0x10000)
                {
                    bnum = ((uint32)tmp[4] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[4];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[5];
                if (root < linesize)
                    //if (mask & 0x100000)
                {
                    bnum = ((uint32)tmp[5] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[5];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[6];
                if (root < linesize)
                    //if (mask & 0x1000000)
                {
                    bnum = ((uint32)tmp[6] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[6];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[7];
                if (root < linesize)
                    //if (mask & 0x10000000)
                {
                    bnum = ((uint32)tmp[7] >> FLAGBITS);
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
                uint64 tmp3;
                uint32 p2, r2;
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

                tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
                tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff);

                root = (uint32)(tmp2 % (uint64)prime);
                r2 = (uint32)(tmp3 % (uint64)p2);

                // It is faster to update during
                // linesieve than doing it all here in a loop.
                // measured 6/2016
                if (root < linesize)
                {
                    bnum = (root >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = ((uint64)prime << 32) | (uint64)root;
                    nptr[bnum]++;
                }

                if (r2 < linesize)
                {
                    bnum = (r2 >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = ((uint64)p2 << 32) | (uint64)r2;
                    nptr[bnum]++;
                }

            }
        }

#else

        for (; i < sdata->bitmap_start_id - 1; i += 2)
        {
            uint64 tmp3;
            uint32 p2, r2;
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

            // Maybe it can still be a win in some cases to use Monty-arith even without AVX2?
            // Until we can benchmark that, just skip it unless AVX2 is enabled.
#ifdef NOT_USE_MONTY

            {
                uint32 m, d, t;

                // inputs to monty     
                t = to_monty_loc(lmp[i] + diff, sdata->r2modp[i], sdata->pinv[i], prime);
                //s = to_monty_loc(s, sdata->r2modp[i], sdata->pinv[i], prime);
                tmp2 = (uint64)s * (uint64)t;

                // reduce
                root = redc_loc(tmp2, sdata->pinv[i], prime);

                // take out of monty rep
                root = redc_loc(root, sdata->pinv[i], prime);

                t = to_monty_loc(lmp[i+1] + diff, sdata->r2modp[i+1], sdata->pinv[i+1], p2);
                //s2 = to_monty_loc(s2, sdata->r2modp[i+1], sdata->pinv[i+1], p2);
                tmp3 = (uint64)s2 * (uint64)t;

                // reduce
                r2 = redc_loc(tmp3, sdata->pinv[i+1], p2);

                // take out of monty rep
                r2 = redc_loc(r2, sdata->pinv[i+1], p2);
            }

#else
            tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
            tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff);

            root = (uint32)(tmp2 % (uint64)prime);
            r2 = (uint32)(tmp3 % (uint64)p2);
#endif

            // It is faster to update during
            // linesieve than doing it all here in a loop.
            // measured 6/2016
            if (root < linesize)
            {
                bnum = (root >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64)prime << 32) | (uint64)root;
                nptr[bnum]++;
            }

            if (r2 < linesize)
            {
                bnum = (r2 >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64)p2 << 32) | (uint64)r2;
                nptr[bnum]++;
            }

        }

        if ((i < sdata->bitmap_start_id) && (ddata->largep_offset == 0))
        {
            uint32 *lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;
            prime = sdata->sieve_p[i];

            s = sdata->root[i];

            tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
            root = (uint32)(tmp2 % (uint64)prime);

            nptr = ddata->bucket_hits;
            bptr = ddata->sieve_buckets;

            if (root < linesize)
            {
                bnum = (root >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64)prime << 32) | (uint64)root;
                nptr[bnum]++;
            }

        }
#endif



        if (ddata->largep_offset > 0)
        {
            // primes greater than the entire sieve interval, thus they
            // at most hit one block and we don't need to save the prime
            // itself since it doesn't need to be advanced.
            uint32 **large_bptr;
            uint32 *large_nptr;
            uint32 *lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;

            large_nptr = ddata->large_bucket_hits;
            large_bptr = ddata->large_sieve_buckets;

            for (i = 0; i < sdata->blocks; i++)
            {
                //initialize bucket
                large_nptr[i] = 0;
            }

#if defined(USE_AVX2)

            if (sdata->use_monty)
            {
                // AVX2 (in soe_util.c) has arranged things so that the end of the
                // sieve prime array is at an index divisible by 8.
                for (i = ddata->largep_offset; i<sdata->bitmap_start_id; i+=8)
                {
                    __m256i vp = _mm256_loadu_si256((__m256i *)(&sdata->sieve_p[i]));
                    __m256i vpinv = _mm256_loadu_si256((__m256i *)(&sdata->pinv[i]));
                    __m256i vr2 = _mm256_loadu_si256((__m256i *)(&sdata->r2modp[i]));
                    __m256i vr = _mm256_loadu_si256((__m256i *)(&sdata->root[i]));
                    __m256i vlmp = _mm256_loadu_si256((__m256i *)(&lmp[i]));
                    __m256i vdiff = _mm256_set1_epi32(diff);
                    __m256i t1, t2, t3;
                    __m256i even, odd;
                    ALIGNED_MEM uint32 tmp[8];
                    uint32 mask;
                    int j;

                    // lmp[i:i+7] + diff
                    vlmp = _mm256_add_epi32(vlmp, vdiff);

                    // lmp to monty     
                    t1 = vec_to_monty(vlmp, vr2, vpinv, vp);

                    //tmp2 = (uint64)s * (uint64)t;
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

                    _mm256_storeu_si256((__m256i *)tmp, vr);

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
                for (i = ddata->largep_offset; i<sdata->bitmap_start_id; i+=2)
                {
                    uint64 tmp3;
                    uint32 p2, r2;
                    int s2;

                    prime = sdata->sieve_p[i];
                    p2 = sdata->sieve_p[i + 1];

                    s = sdata->root[i];
                    s2 = sdata->root[i + 1];

                    // we solved for lower_mod_prime while computing the modular inverse of
                    // each prime, for the residue class 1.  add the difference between this
                    // residue class and 1 before multiplying by the modular inverse.
                    // could use (_mm256_mul_epu32 --> VPMULUDQ)
                    tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
                    tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff);

                    // would need custom solution
                    root = (uint32)(tmp2 % (uint64)prime);
                    r2 = (uint32)(tmp3 % (uint64)p2);

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
            }


#else

			for (i = ddata->largep_offset; i<sdata->bitmap_start_id-1; i+=2)
			{
				uint64 tmp3;
				uint32 p2, r2;
				int s2;
							
				prime = sdata->sieve_p[i];
				p2 = sdata->sieve_p[i+1];

				s = sdata->root[i];
				s2 = sdata->root[i+1];
				
				// we solved for lower_mod_prime while computing the modular inverse of
				// each prime, for the residue class 1.  add the difference between this
				// residue class and 1 before multiplying by the modular inverse.
                // could use (_mm256_mul_epu32 --> VPMULUDQ)
				tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
				tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff);

                // would need custom solution
				root = (uint32)(tmp2 % (uint64)prime);
				r2 = (uint32)(tmp3 % (uint64)p2);

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

                tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
                root = (uint32)(tmp2 % (uint64)prime);

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

void get_offsets2(thread_soedata_t *thread_data)
{
    //extract stuff from the thread data structure
    soe_dynamicdata_t *ddata = &thread_data->ddata;
    soe_staticdata_t *sdata = &thread_data->sdata;

    uint64 i, startprime = sdata->startprime, prodN = sdata->prodN;
    uint32 block = ddata->blockstart;
    uint32 prime, root, bnum;
    uint32 diff = sdata->rclass[thread_data->current_line] - 1;
    uint64 tmp2;
    int s;

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
        uint32 *lmp = sdata->lower_mod_prime;

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
                ddata->pbounds[block] = i;
                block++;
                ddata->lblk_b = ddata->ublk_b + prodN;
                ddata->ublk_b += sdata->blk_r;
                ddata->blk_b_sqrt = (uint64)(sqrt((int64)(ddata->ublk_b + prodN))) + 1;
            }

            s = sdata->root[i];

            // the lower block bound (lblk_b) times s can exceed 64 bits for large ranges,
            // so reduce mod p here as well.
            tmp2 = (uint64)s * (ddata->lblk_b % (uint64)prime);

            // tmp2 = (uint64)s * (uint64)(lmp[i] + diff);
            ddata->offsets[i] = (uint32)(tmp2 % (uint64)prime);
        }
    }
    else
    {
        uint32 modp;
        mpz_t lowz, sqrtz;
        mpz_init(lowz);
        mpz_init(sqrtz);
        mpz_set(lowz, *sdata->offset);
        mpz_add_ui(lowz, lowz, ddata->lblk_b);

        mpz_set(sqrtz, lowz);
        mpz_add_ui(sqrtz, sqrtz, sdata->blk_r);
        mpz_sqrt(sqrtz, sqrtz);
        mpz_add_ui(sqrtz, sqrtz, 1);
        //mpz_set_ui(tmpz, ddata->lblk_b);

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
            tmp2 = (uint64)s * (uint64)modp;
            ddata->offsets[i] = (uint32)(tmp2 % (uint64)prime);
            //gmp_printf("p = %u, o = %u, r = %d, lblk_b = %Zd, modp = %u\n", 
            //	prime, ddata->offsets[i], s, tmpz, modp);
        }

        mpz_clear(lowz);
        mpz_clear(sqrtz);
    }

    if (ddata->bucket_depth > 0)
    {
        uint64 **bptr;

        uint32 *nptr;
        uint32 rangesize = FLAGSIZE * (ddata->blockstop - ddata->blockstart);
        uint32 *lmp = sdata->lower_mod_prime;

        nptr = ddata->bucket_hits;
        bptr = ddata->sieve_buckets;

#if defined(USE_AVX2)

        if (sdata->use_monty)
        {
            // AVX2 (in soe_util.c) has arranged things so that the end of the
            // sieve prime array is at an index divisible by 8.
            for (; i < sdata->bitmap_start_id; i += 8)
            {
                __m256i vp = _mm256_loadu_si256((__m256i *)(&sdata->sieve_p[i]));
                __m256i vpinv = _mm256_loadu_si256((__m256i *)(&sdata->pinv[i]));
                __m256i vr2 = _mm256_loadu_si256((__m256i *)(&sdata->r2modp[i]));
                __m256i vr = _mm256_loadu_si256((__m256i *)(&sdata->root[i]));
                __m256i vlmp = _mm256_loadu_si256((__m256i *)(&lmp[i]));
                __m256i vdiff = _mm256_set1_epi32(diff + FLAGSIZE * ddata->blockstart);
                __m256i t1, t2, t3;
                __m256i even, odd;
                ALIGNED_MEM uint64 tmp[8];
                ALIGNED_MEM uint32 tmp32[8];
                uint32 mask;

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

                //tmp2 = (uint64)s * (uint64)t;
                t3 = _mm256_shuffle_epi32(t1, 0xB1);
                t2 = _mm256_shuffle_epi32(vr, 0xB1);

                even = _mm256_mul_epu32(vr, t1);
                odd = _mm256_mul_epu32(t3, t2);

                // reduce
                vr = vec_redc(even, odd, vpinv, vp);

                // take out of monty rep
                vr = vec_redc(CLEAR_HIGH_VEC(vr), CLEAR_HIGH_VEC(_mm256_shuffle_epi32(vr, 0xB1)), vpinv, vp);

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

                _mm256_store_si256((__m256i *)tmp, t2);         // even (prime | root)
                _mm256_store_si256((__m256i *)(&tmp[4]), t3);   // odd  (prime | root)

                // bucket sort the roots
                root = (uint32)tmp[0];
                if (root < rangesize)
                    //if (mask & 0x1)
                {
                    bnum = ((uint32)tmp[0] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[0];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[1];
                if (root < rangesize)
                    //if (mask & 0x10)
                {
                    bnum = ((uint32)tmp[1] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[1];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[2];
                if (root < rangesize)
                    //if (mask & 0x100)
                {
                    bnum = ((uint32)tmp[2] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[2];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[3];
                if (root < rangesize)
                    //if (mask & 0x1000)
                {
                    bnum = ((uint32)tmp[3] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[3];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[4];
                if (root < rangesize)
                    //if (mask & 0x10000)
                {
                    bnum = ((uint32)tmp[4] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[4];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[5];
                if (root < rangesize)
                    //if (mask & 0x100000)
                {
                    bnum = ((uint32)tmp[5] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[5];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[6];
                if (root < rangesize)
                    //if (mask & 0x1000000)
                {
                    bnum = ((uint32)tmp[6] >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = tmp[6];
                    nptr[bnum]++;
                }

                root = (uint32)tmp[7];
                if (root < rangesize)
                    //if (mask & 0x10000000)
                {
                    bnum = ((uint32)tmp[7] >> FLAGBITS);
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
                uint64 tmp3;
                uint32 p2, r2;
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

                tmp2 = (uint64)s * (uint64)(lmp[i] + diff + FLAGSIZE * ddata->blockstart);
                tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff + FLAGSIZE * ddata->blockstart);

                root = (uint32)(tmp2 % (uint64)prime);
                r2 = (uint32)(tmp3 % (uint64)p2);

                // It is faster to update during
                // linesieve than doing it all here in a loop.
                // measured 6/2016
                if (root < rangesize)
                {
                    bnum = (root >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = ((uint64)prime << 32) | (uint64)root;
                    nptr[bnum]++;
                }

                if (r2 < rangesize)
                {
                    bnum = (r2 >> FLAGBITS);
                    bptr[bnum][nptr[bnum]] = ((uint64)p2 << 32) | (uint64)r2;
                    nptr[bnum]++;
                }

            }
        }

#else

        for (; i < sdata->bitmap_start_id - 1; i += 2)
        {
            uint64 tmp3;
            uint32 p2, r2;
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

            // Maybe it can still be a win in some cases to use Monty-arith even without AVX2?
            // Until we can benchmark that, just skip it unless AVX2 is enabled.
#ifdef NOT_USE_MONTY

            {
                uint32 m, d, t;

                // inputs to monty     
                t = to_monty_loc(lmp[i] + diff, sdata->r2modp[i], sdata->pinv[i], prime);
                //s = to_monty_loc(s, sdata->r2modp[i], sdata->pinv[i], prime);
                tmp2 = (uint64)s * (uint64)t;

                // reduce
                root = redc_loc(tmp2, sdata->pinv[i], prime);

                // take out of monty rep
                root = redc_loc(root, sdata->pinv[i], prime);

                t = to_monty_loc(lmp[i + 1] + diff, sdata->r2modp[i + 1], sdata->pinv[i + 1], p2);
                //s2 = to_monty_loc(s2, sdata->r2modp[i+1], sdata->pinv[i+1], p2);
                tmp3 = (uint64)s2 * (uint64)t;

                // reduce
                r2 = redc_loc(tmp3, sdata->pinv[i + 1], p2);

                // take out of monty rep
                r2 = redc_loc(r2, sdata->pinv[i + 1], p2);
            }

#else
            tmp2 = (uint64)s * (uint64)(lmp[i] + diff + FLAGSIZE * ddata->blockstart);
            tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff + FLAGSIZE * ddata->blockstart);

            root = (uint32)(tmp2 % (uint64)prime);
            r2 = (uint32)(tmp3 % (uint64)p2);
#endif

            // It is faster to update during
            // linesieve than doing it all here in a loop.
            // measured 6/2016
            if (root < rangesize)
            {
                bnum = (root >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64)prime << 32) | (uint64)root;
                nptr[bnum]++;
            }

            if (r2 < rangesize)
            {
                bnum = (r2 >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64)p2 << 32) | (uint64)r2;
                nptr[bnum]++;
            }

        }

        if ((i < sdata->bitmap_start_id) && (ddata->largep_offset == 0))
        {
            uint32 *lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;
            prime = sdata->sieve_p[i];

            s = sdata->root[i];

            tmp2 = (uint64)s * (uint64)(lmp[i] + diff + FLAGSIZE * ddata->blockstart);
            root = (uint32)(tmp2 % (uint64)prime);

            nptr = ddata->bucket_hits;
            bptr = ddata->sieve_buckets;

            if (root < rangesize)
            {
                bnum = (root >> FLAGBITS);
                bptr[bnum][nptr[bnum]] = ((uint64)prime << 32) | (uint64)root;
                nptr[bnum]++;
            }

        }
#endif



        if (ddata->largep_offset > 0)
        {
            // primes greater than the entire sieve interval, thus they
            // at most hit one block and we don't need to save the prime
            // itself since it doesn't need to be advanced.
            uint32 **large_bptr;
            uint32 *large_nptr;
            uint32 *lmp = sdata->lower_mod_prime;// -sdata->bucket_start_id;

            large_nptr = ddata->large_bucket_hits;
            large_bptr = ddata->large_sieve_buckets;

            for (i = 0; i < sdata->blocks; i++)
            {
                //initialize bucket
                large_nptr[i] = 0;
            }

#if defined(USE_AVX2)

            if (sdata->use_monty)
            {
                // AVX2 (in soe_util.c) has arranged things so that the end of the
                // sieve prime array is at an index divisible by 8.
                for (i = ddata->largep_offset; i<sdata->bitmap_start_id; i += 8)
                {
                    __m256i vp = _mm256_loadu_si256((__m256i *)(&sdata->sieve_p[i]));
                    __m256i vpinv = _mm256_loadu_si256((__m256i *)(&sdata->pinv[i]));
                    __m256i vr2 = _mm256_loadu_si256((__m256i *)(&sdata->r2modp[i]));
                    __m256i vr = _mm256_loadu_si256((__m256i *)(&sdata->root[i]));
                    __m256i vlmp = _mm256_loadu_si256((__m256i *)(&lmp[i]));
                    __m256i vdiff = _mm256_set1_epi32(diff + FLAGSIZE * ddata->blockstart);
                    __m256i t1, t2, t3;
                    __m256i even, odd;
                    ALIGNED_MEM uint32 tmp[8];
                    uint32 mask;
                    int j;

                    // lmp[i:i+7] + diff
                    vlmp = _mm256_add_epi32(vlmp, vdiff);

                    // lmp to monty     
                    t1 = vec_to_monty(vlmp, vr2, vpinv, vp);

                    //tmp2 = (uint64)s * (uint64)t;
                    t3 = _mm256_shuffle_epi32(t1, 0xB1);
                    t2 = _mm256_shuffle_epi32(vr, 0xB1);

                    even = _mm256_mul_epu32(vr, t1);
                    odd = _mm256_mul_epu32(t3, t2);

                    // reduce
                    vr = vec_redc(even, odd, vpinv, vp);

                    // take out of monty rep
                    vr = vec_redc(CLEAR_HIGH_VEC(vr), CLEAR_HIGH_VEC(_mm256_shuffle_epi32(vr, 0xB1)), vpinv, vp);

                    t1 = _mm256_set1_epi32(rangesize);
                    t1 = _mm256_cmpge_epu32(vr, t1);
                    mask = ~(_mm256_movemask_epi8(t1));

                    _mm256_storeu_si256((__m256i *)tmp, vr);

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
                for (i = ddata->largep_offset; i<sdata->bitmap_start_id; i += 2)
                {
                    uint64 tmp3;
                    uint32 p2, r2;
                    int s2;

                    prime = sdata->sieve_p[i];
                    p2 = sdata->sieve_p[i + 1];

                    s = sdata->root[i];
                    s2 = sdata->root[i + 1];

                    // we solved for lower_mod_prime while computing the modular inverse of
                    // each prime, for the residue class 1.  add the difference between this
                    // residue class and 1 before multiplying by the modular inverse.
                    // could use (_mm256_mul_epu32 --> VPMULUDQ)
                    tmp2 = (uint64)s * (uint64)(lmp[i] + diff + FLAGSIZE * ddata->blockstart);
                    tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff + FLAGSIZE * ddata->blockstart);

                    // would need custom solution
                    root = (uint32)(tmp2 % (uint64)prime);
                    r2 = (uint32)(tmp3 % (uint64)p2);

                    // gather may help, but writes would need to be done 1 by 1.
                    if (root < rangesize)
                    {
                        bnum = (root >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = root;
                        large_nptr[bnum]++;
                    }

                    if (r2 < rangesize)
                    {
                        bnum = (r2 >> FLAGBITS);
                        large_bptr[bnum][large_nptr[bnum]] = r2;
                        large_nptr[bnum]++;
                    }
                }
            }


#else

            for (i = ddata->largep_offset; i<sdata->bitmap_start_id - 1; i += 2)
            {
                uint64 tmp3;
                uint32 p2, r2;
                int s2;

                prime = sdata->sieve_p[i];
                p2 = sdata->sieve_p[i + 1];

                s = sdata->root[i];
                s2 = sdata->root[i + 1];

                // we solved for lower_mod_prime while computing the modular inverse of
                // each prime, for the residue class 1.  add the difference between this
                // residue class and 1 before multiplying by the modular inverse.
                // could use (_mm256_mul_epu32 --> VPMULUDQ)
                tmp2 = (uint64)s * (uint64)(lmp[i] + diff + FLAGSIZE * ddata->blockstart);
                tmp3 = (uint64)s2 * (uint64)(lmp[i + 1] + diff + FLAGSIZE * ddata->blockstart);

                // would need custom solution
                root = (uint32)(tmp2 % (uint64)prime);
                r2 = (uint32)(tmp3 % (uint64)p2);

                // gather may help, but writes would need to be done 1 by 1.
                if (root < rangesize)
                {
                    bnum = (root >> FLAGBITS);
                    large_bptr[bnum][large_nptr[bnum]] = root;
                    large_nptr[bnum]++;
                }

                if (r2 < rangesize)
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

                tmp2 = (uint64)s * (uint64)(lmp[i] + diff + FLAGSIZE * ddata->blockstart);
                root = (uint32)(tmp2 % (uint64)prime);

                if (root < rangesize)
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