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
#include "soe_impl.h"
#include "gmp.h"
#include "ytools.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "threadpool.h"



// known issues:
// sieve_to_depth crashes with offsets too big (somewhere above 2^80)
// sieve_to_depth with analysis=2, numclasses=480 and computing PRPs is *way* too slow,
//      something is wrong.  compared to numclasses=48, which seems normal.  same
//		number of candidates in both cases... weird.  Could be a threading thing; 
//		computing PRPs with one thread seems about the same speed for each numclasses.
//		"2^70" "2^70+10^11" -v -v -v -t 8 -p 1000000000 -w 1 -c 480 -b 131072 -a 2
//		"2^70" "2^70+10^11" -v -v -v -t 8 -p 1000000000 -w 1 -c 48 -b 131072 -a 2





void compute_prps_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *t = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = t->sdata;

    if (sdata->sync_count < sdata->THREADS)
    {
        tdata->work_fcn_id = 0;
        sdata->sync_count++;
    }
    else
    {
        tdata->work_fcn_id = tdata->num_work_fcn;
    }

    return;
}

void compute_prps_work_fcn(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    soe_userdata_t *udata = (soe_userdata_t *)tdata->user_data;
    soe_staticdata_t *sdata = udata->sdata;
    thread_soedata_t *t = &udata->ddata[tdata->tindex];
	int witnesses = sdata->witnesses;
    int i;



    
#if 0

#ifndef IFMA
	dbias = _mm512_castsi512_pd(set64(0x4670000000000000ULL));
	vbias1 = set64(0x4670000000000000ULL);
	vbias2 = set64(0x4670000000000001ULL);
	vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
#endif

    t->linecount = 0;
	for (i = t->startid; i < t->stopid - 8; i += 8)
	{
		if ((((i - t->startid) & 8191) == 0) && (sdata->VFLAG > 0))
		{
			printf("thread %d progress: %d%%\r", tdata->tindex,
				(int)((double)(i - t->startid) / (double)(t->stopid - t->startid) * 100.0));
			fflush(stdout);
		}

		ALIGNED_MEM uint64_t n8[16];
		int j;
		uint8_t valid_msk = 0;

		for (j = 0; j < 8; j++)
		{
			mpz_add_ui(t->tmpz, t->offset, t->ddata.primes[i + j - t->startid]);

			if ((mpz_cmp(t->tmpz, t->lowlimit) >= 0) && (mpz_cmp(t->highlimit, t->tmpz) >= 0))
			{
				n8[j] = mpz_get_ui(t->tmpz) & 0xfffffffffffffull;
				mpz_tdiv_q_2exp(t->tmpz, t->tmpz, 52);
				n8[j + 8] = mpz_get_ui(t->tmpz) & 0xfffffffffffffull;
				valid_msk |= (1 << j);
			}
		}

		uint8_t prpmask = valid_msk & MR_2sprp_104x8(n8);
		t->linecount += _mm_popcnt_u32(prpmask);

		//for (j = 0; j < 8; j++)
		//{
		//	if (prpmask & (1 << j))
		//	{
		//		t->ddata.primes[t->linecount++] = t->ddata.primes[i + j - t->startid];
		//	}
		//}

	}

	for ( ; i < t->stopid; i++)
#else
	for (i = t->startid; i < t->stopid; i++)
#endif
    {
        if (((i & 8191) == 0) && (sdata->VFLAG > 0))
        {
            printf("thread %d progress: %d%%\r", tdata->tindex, 
                (int)((double)(i - t->startid) / (double)(t->stopid - t->startid) * 100.0));
            fflush(stdout);
        }

        mpz_add_ui(t->tmpz, t->offset, t->ddata.primes[i - t->startid]);
        if ((mpz_cmp(t->tmpz, t->lowlimit) >= 0) && (mpz_cmp(t->highlimit, t->tmpz) >= 0))
        {
			//gmp_printf("candidate %Zd... ", t->tmpz);
            //if (mpz_extrastrongbpsw_prp(t->tmpz))
            if (mpz_probab_prime_p(t->tmpz, witnesses))
            {
				if (sdata->analysis == 2)
				{
					// also need to check the twin
					//gmp_printf("and twin %Zd is...", t->tmpz);
					mpz_add_ui(t->tmpz, t->tmpz, 2);
					if (mpz_probab_prime_p(t->tmpz, sdata->witnesses))
					{
						t->ddata.primes[t->linecount++] = t->ddata.primes[i - t->startid];
						//printf("prime!\n");
					}
				}
				else
				{
					//t->ddata.primes[t->linecount++] = t->ddata.primes[i - t->startid];
					t->linecount++;
				}
				//printf("prime!\n");
            }
			else
			{
				//printf("not prime\n");
			}
        }
    }

    return;
}

soe_staticdata_t* soe_init(int vflag, int threads, int blocksize)
{
    soe_staticdata_t* sdata;

    sdata = (soe_staticdata_t*)malloc(sizeof(soe_staticdata_t));

	//int i;
	//int pn = 2310 * 13;
	//int nc = 0;
	//for (i = 1; i < pn; i++)
	//{
	//	if (gcd_1(i, pn) == 1)
	//	{
	//		printf("class: %d, inv: %d\n", i,
	//			pn - modinv1(i, pn));
	//		nc++;
	//	}		
	//}
	//printf("found %d classes\n", nc);
	//exit(1);

    // bootstrap the sieve
    sdata->sieve_p = (uint32_t*)xmalloc(65536 * sizeof(uint32_t));
    sdata->num_sp = tiny_soe(66000, sdata->sieve_p);

	// default settings: compute primes with default number of classes
	sdata->witnesses = 0;
	sdata->userclasses = 0;
	sdata->analysis = 1;
	sdata->gapmin = 0;

    sdata->VFLAG = vflag;
    sdata->THREADS = threads;
    if (blocksize > 1024)
        sdata->SOEBLOCKSIZE = blocksize;
    else
        sdata->SOEBLOCKSIZE = blocksize << 10;

	mpz_init(sdata->offset);

    return sdata;
}

void soe_finalize(soe_staticdata_t* sdata)
{
    free(sdata->sieve_p);
	mpz_clear(sdata->offset);
	free(sdata);
    return;
}

uint64_t *GetPRIMESRange(soe_staticdata_t* sdata, 
	mpz_t offset, uint64_t lowlimit, uint64_t highlimit, uint64_t *num_p)
{
	uint64_t i;
	uint64_t hi_est, lo_est;
	uint64_t maxrange = 10000000000ULL;
	uint64_t *primes = NULL;
	
	//reallocate output array based on conservative estimate of the number of 
	//primes in the interval
	if (mpz_cmp_ui(offset, 0) > 0)
	{
		mpz_t a, b;
		mpz_init(a);
		mpz_init(b);
		mpz_add_ui(a, offset, lowlimit);
		mpz_add_ui(b, offset, highlimit);
		i = mpz_estimate_primes_in_range(a, b);
		mpz_clear(a);
		mpz_clear(b);
		primes = (uint64_t *)realloc(primes, (size_t) (i * sizeof(uint64_t)));

		printf("allocating space for an estimated %lu primes in requested range\n", i);

		if (primes == NULL)
		{
            if (mpz_cmp_ui(offset, 0) > 0)
            {
                printf("unable to allocate %" PRIu64 " bytes for range %" PRIu64 " to %" PRIu64 "\n",
                    (uint64_t)(i * sizeof(uint64_t)), lowlimit, highlimit);
            }
            else
            {
                printf("unable to allocate %" PRIu64 " bytes \n",
                    (uint64_t)(i * sizeof(uint64_t)));
            }
			exit(1);
		}
	}
	else
	{
		i = estimate_primes_in_range(lowlimit, highlimit);
		primes = (uint64_t *)xrealloc(primes, (size_t) (i * sizeof(uint64_t)));
	}

	//check for really big ranges ('big' is different here than when we are counting
	//primes because there are higher memory demands when computing primes)
	if (0) //((highlimit - lowlimit) > maxrange)
	{
		uint64_t tmpl, tmph, tmpcount = 0;
		uint32_t num_ranges = (uint32_t)((highlimit - lowlimit) / maxrange);
		uint64_t remainder = (highlimit - lowlimit) % maxrange;
		uint32_t j;
				
		sdata->GLOBAL_OFFSET = 0;
		tmpl = lowlimit;
        // maxrange - 1, so that we don't count the upper
        // limit twice (again on the next iteration's lower bound).
		tmph = lowlimit + maxrange - 1;
		for (j = 0; j < num_ranges; j++)
		{
			tmpcount += spSOE(sdata, offset, tmpl, &tmph, 0, primes);
			tmpl += maxrange;
			tmph = tmpl + maxrange - 1;
            sdata->GLOBAL_OFFSET = tmpcount;
		}
				
		tmph = tmpl + remainder;
		tmpcount += spSOE(sdata, offset, tmpl, &tmph, 0, primes);
		*num_p = tmpcount;
	}
	else
	{
		//find the primes in the interval
        sdata->GLOBAL_OFFSET = 0;
        if (sdata->VFLAG > 1)
        {
            printf("generating primes in range %" PRIu64 " : %" PRIu64 "\n", 
                lowlimit, highlimit);
        }
		*num_p = spSOE(sdata, offset, lowlimit, &highlimit, 0, primes);
	}

	return primes;
}

uint64_t *soe_wrapper(soe_staticdata_t* sdata, uint64_t lowlimit, uint64_t highlimit, 
    int count, uint64_t* num_p, int PRIMES_TO_FILE, int PRIMES_TO_SCREEN)
{
	//public interface to the sieve.  
	uint64_t retval, tmpl, tmph, i;
	uint32_t max_p;	
	mpz_t offset;
	uint64_t *primes = NULL;

    sdata->only_count = count;
	//if (count)
	//{
	//	sdata->analysis = 0;
	//}

    if (highlimit < lowlimit)
    {
        printf("error: lowlimit must be less than highlimit\n");
        *num_p = 0;
        return primes;
    }

	mpz_t a, b;
	mpz_init(a);
	mpz_init(b);
	mpz_init(offset);
	mpz_set_ui(offset, 0);
	mpz_set_ui(a, highlimit);
	mpz_set_ui(b, sdata->sieve_p[sdata->num_sp - 1]);
	mpz_mul_ui(b, b, sdata->sieve_p[sdata->num_sp - 1]);
	retval = (mpz_cmp(a, b) > 0);
	

	//printf("highlimit = %lu, num_sp = %u, spmax = %u\n",
	//	highlimit, sdata->num_sp, sdata->sieve_p[sdata->num_sp - 1]);

	//if (highlimit > (sdata->sieve_p[sdata->num_sp-1] * sdata->sieve_p[sdata->num_sp-1]))
	if (retval)
	{
		//then we need to generate more sieving primes
		uint32_t range_est;

		//allocate array based on conservative estimate of the number of 
		//primes in the interval	
		mpz_sqrt(b, a);
		max_p = (uint32_t)mpz_get_ui(b) + 65536;
		range_est = (uint32_t)estimate_primes_in_range(0, (uint64_t)max_p);

        if (sdata->VFLAG > 1)
        {
            printf("generating more sieving primes in range 0 : %u \n", max_p);
            printf("allocating %u bytes \n", range_est);
        }

        sdata->sieve_p = (uint32_t *)xrealloc(sdata->sieve_p, 
            (size_t) (range_est * sizeof(uint32_t)));

		//find the sieving primes using the seed primes
        sdata->NO_STORE = 0;
		sdata->is_main_sieve = 0;
		primes = GetPRIMESRange(sdata, offset, 0, max_p, &retval);

        if (sdata->VFLAG > 1)
        {
            printf("found %u sieving primes\n", (uint32_t)retval);
        }

        for (i = 0; i < retval; i++)
        {
            sdata->sieve_p[i] = (uint32_t)primes[i];
        }

        sdata->num_sp = (uint32_t)retval;
		free(primes);
		primes = NULL;
	}

	mpz_clear(a);
	mpz_clear(b);

	if (count)
	{
		sdata->is_main_sieve = 1;
		*num_p = spSOE(sdata, offset, lowlimit, &highlimit, count, NULL);
	}
	else
	{
		sdata->is_main_sieve = 1;
		primes = GetPRIMESRange(sdata, offset, lowlimit, highlimit, num_p);

		// now dump the requested range of primes to a file, or the
		// screen, both, or neither, depending on the state of a couple
		// global configuration variables
		if (PRIMES_TO_FILE)
		{
			FILE *out;
			out = fopen("primes.dat","w");
			if (out == NULL)
			{
				printf("fopen error: %s\n", strerror(errno));
				printf("can't open primes.dat for writing\n");
			}
			else
			{
				for (i = 0; i < *num_p; i++)
				{
					if (sdata->analysis == 1)
					{
						if ((primes[i] >= lowlimit) && (primes[i] <= highlimit))
						{
							fprintf(out, "%" PRIu64 "\n", primes[i]);
						}
					}
					else if (sdata->analysis == 2)
					{
						if ((primes[i] >= lowlimit) && (primes[i] <= highlimit))
						{
							fprintf(out, "(%" PRIu64 ", %" PRIu64 ")\n", 
								primes[i], primes[i] + 2);
						}
						else
						{
							printf("out of bounds prime %lu found in list\n", primes[i]);
						}
					}
				}
				fclose(out);
			}
		}

		if (PRIMES_TO_SCREEN)
		{
			for (i = 0; i < *num_p; i++)
			{
                if ((primes[i] >= lowlimit) && (primes[i] <= highlimit))
                {
                    printf("%" PRIu64 " ", primes[i]);
                }
			}
			printf("\n");
		}			
	}

	return primes;
}

uint64_t *sieve_to_depth(soe_staticdata_t* sdata,
	mpz_t lowlimit, mpz_t highlimit, int count, int num_witnesses, 
    uint64_t sieve_limit, uint64_t *num_p,
    int PRIMES_TO_FILE, int PRIMES_TO_SCREEN)
{
	// public interface to a routine which will sieve a range of integers
	// with the supplied primes and either count or compute the values
	// that survive.  Basically, it is just the sieve, but with no
	// guareentees that what survives the sieving is prime.  The idea is to 
	// remove cheap composites.
	uint64_t retval, i, range, tmpl, tmph, num_sp_needed;
	uint64_t *values = NULL;
	mpz_t tmpz;
	mpz_t offset;

	if (mpz_cmp(highlimit, lowlimit) <= 0)
	{
		printf("error: lowlimit must be less than highlimit\n");
		*num_p = 0;
		return values;
	}	

	mpz_init(tmpz);
	mpz_init(offset);
	mpz_set(offset, lowlimit);
	mpz_sub(tmpz, highlimit, lowlimit);
	range = mpz_get_ui(tmpz);

	// sieve with the range requested, down to a minimum of 1000.
	// (so that we don't run into issues with the presieve.)
	if ((sieve_limit + 10000) < sdata->sieve_p[sdata->num_sp - 1])
	{
		// if the seed primes we were provided are enough then
		// just use them
	}
	else
	{
		// then we need to generate more sieving primes
		uint32_t range_est;

		// allocate array based on conservative estimate of the number of 
		// primes in the interval	
		range_est = (uint32_t)estimate_primes_in_range(0, sieve_limit + 10000);

		if (sdata->VFLAG > 1)
		{
			printf("generating more sieving primes in range 0 : %u \n", sieve_limit);
			printf("allocating %u bytes \n", range_est);
		}

		sdata->sieve_p = (uint32_t*)xrealloc(sdata->sieve_p,
			(size_t)(range_est * sizeof(uint32_t)));

		// find the sieving primes using the seed primes.
		// we find slightly more than we have been requested to use,
		// because the vector routines in the sieve may want to slightly 
		// increase the number that we sieve with in order to work
		// with full vectors.
		mpz_t tmp_offset;
		mpz_init(tmp_offset);
		mpz_set_ui(tmp_offset, 0);

		sdata->NO_STORE = 0;
		sdata->is_main_sieve = 0;
		uint64_t *primes = GetPRIMESRange(sdata, tmp_offset, 0, sieve_limit + 10000, &retval);

		// from the oversieved list of primes, find how many will satisfy the 
		// requested sieve_primes limit.
		num_sp_needed = retval - 1;
		while (primes[num_sp_needed] > sieve_limit)
		{
			num_sp_needed--;
		}
#ifdef USE_AVX2
		num_sp_needed += 8;
#else
		num_sp_needed++;
#endif

		if (sdata->VFLAG > 1)
		{
			printf("found %u sieving primes, max prime = %u\n", 
				(uint32_t)num_sp_needed, primes[num_sp_needed - 1]);
		}

		// copy over all the primes
		sdata->sieve_p = (uint32_t*)xrealloc(sdata->sieve_p, retval * sizeof(uint32_t));
		for (i = 0; i < retval; i++)
		{
			sdata->sieve_p[i] = (uint32_t)primes[i];
		}

		// but list the number needed, which should be slightly smaller than those found.
		sdata->num_sp = (uint32_t)num_sp_needed;
		sdata->alloc_sp = retval;
		free(primes);
		mpz_clear(tmp_offset);
		primes = NULL;
	}

	if (count)
	{
		gmp_printf("commencing sieve over interval %Zd + (%lu:%lu) with %u sieve primes\n",
			offset, 0, range, sdata->num_sp);

		if (num_witnesses > 0)
		{
			sdata->witnesses = num_witnesses;
			switch (num_witnesses)
			{
			case 1: printf("verifying candidates with base-2 Fermat PRP check\n"); break;
			default: printf("verifying candidates with base-2 strong PRP (MR) check\n"); break;
			}
			
		}
	}
	else
	{
		gmp_printf("generating primes in interval %Zd + (%lu:%lu) with %u sieve primes and num_witness = %d\n",
			offset, 0, range, sdata->num_sp, num_witnesses);
	}

	if (count)
	{
		sdata->is_main_sieve = 1;
		*num_p = spSOE(sdata, offset, 0, &range, 1, NULL);
	}
	else
	{
		sdata->is_main_sieve = 1;
		values = GetPRIMESRange(sdata, offset, 0, range, num_p);

		if (num_witnesses > 0)
		{
			thread_soedata_t *thread_data;		//an array of thread data objects
			uint32_t lastid;
			int j;

            // threading structures
            tpool_t *tpool_data;
            soe_userdata_t udata;

			//allocate thread data structure
			thread_data = (thread_soedata_t *)malloc(sdata->THREADS * sizeof(thread_soedata_t));
			
			// conduct PRP tests on all surviving values
            if (sdata->VFLAG > 0)
            {
                gmp_printf("starting PRP tests with %d witnesses on "
                    "%" PRIu64 " surviving candidates using %d threads\n",
                    num_witnesses, *num_p, sdata->THREADS);
				fflush(stdout);
            }

			range = *num_p / sdata->THREADS;
			lastid = 0;

			// divvy up the range
			for (j = 0; j < sdata->THREADS; j++)
			{
				thread_soedata_t *t = thread_data + j;
				
				t->startid = lastid;
				t->stopid = t->startid + range;
				lastid = t->stopid;

                if (sdata->VFLAG > 2)
                {
                    printf("thread %d computing PRPs from %u to %u\n",
                        (int)j, t->startid, t->stopid);
					fflush(stdout);
                }
			}

			// the last one gets any leftover
            if (thread_data[sdata->THREADS - 1].stopid != (uint32_t)*num_p)
            {
                thread_data[sdata->THREADS - 1].stopid = (uint32_t)*num_p;
            }

			if (sdata->THREADS == 1)
			{
				mpz_t tmpz;
				mpz_init(tmpz);

				retval = 0;

#if 0
				for (i = 0; i < range - 8; i += 8)
				{
					if (((i & 8191) == 0) && (sdata->VFLAG > 0))
					{
						printf("progress: %d%%\r",
							(int)((double)(i) / (double)(range) * 100.0));
						fflush(stdout);
					}

					ALIGNED_MEM uint64_t n8[16];
					uint8_t loc_msk = 0;
					int j;

					for (j = 0; j < 8; j++)
					{
						mpz_add_ui(tmpz, offset, values[i + j]);

						n8[j] = mpz_get_ui(tmpz) & 0xfffffffffffffull;
						mpz_tdiv_q_2exp(tmpz, tmpz, 52);
						n8[j + 8] = mpz_get_ui(tmpz) & 0xfffffffffffffull;

					}

					uint8_t prpmask = MR_2sprp_104x8(n8);
					retval += _mm_popcnt_u32(prpmask);
				}

				for ( ; i < range; i++)
#else
				for (i = 0; i < range; i += 8)
#endif
				{
					if (((i & 8191) == 0) && (sdata->VFLAG > 0))
					{
					    printf("progress: %d%%\r", 
					        (int)((double)(i) / (double)(range) * 100.0));
					    fflush(stdout);
					}

					mpz_add_ui(tmpz, offset, values[i]);
					if ((mpz_cmp(tmpz, lowlimit) >= 0) && (mpz_cmp(highlimit, tmpz) >= 0))
					{
						//gmp_printf("candidate %Zd is...", tmpz);
						//if (mpz_extrastrongbpsw_prp(t->tmpz))
						if (mpz_probab_prime_p(tmpz, num_witnesses))
						{
							if (sdata->analysis == 2)
							{
								// also need to check the twin
								mpz_add_ui(tmpz, tmpz, 2);
								if (mpz_probab_prime_p(tmpz, num_witnesses))
								{
									values[retval++] = values[i];
								}
							}
							else
							{
								values[retval++] = values[i];
							}
							//printf("prime!\n");
						}
						else
						{
							//printf("not prime\n");
						}
					}
				}
				mpz_clear(tmpz);
				*num_p = retval;
			}
			else
			{
				// allocate space for stuff in the threads
				for (j = 0; j < sdata->THREADS; j++)
				{
					thread_soedata_t* t = thread_data + j;

					mpz_init(t->tmpz);
					mpz_init(t->offset);
					mpz_init(t->lowlimit);
					mpz_init(t->highlimit);
					mpz_set(t->offset, offset);
					mpz_set(t->lowlimit, lowlimit);
					mpz_set(t->highlimit, highlimit);
					t->current_line = (uint64_t)num_witnesses;

					t->ddata.primes = (uint64_t*)malloc((t->stopid - t->startid) * sizeof(uint64_t));
					for (i = t->startid; i < t->stopid; i++)
					{
						t->ddata.primes[i - t->startid] = values[i];
					}
				}

				// now run with the threads.  don't really need the 
				// threadpool since we are statically dividing up the range
				// to test, but it is easy so we use it.
				udata.sdata = &thread_data->sdata;
				udata.ddata = thread_data;
				tpool_data = tpool_setup(sdata->THREADS, NULL, NULL, NULL,
					&compute_prps_dispatch, &udata);

				thread_data->sdata.sync_count = 0;
				tpool_add_work_fcn(tpool_data, &compute_prps_work_fcn);
				tpool_go(tpool_data);

				free(tpool_data);

				// combine results and free stuff
				retval = 0;
				for (i = 0; i < sdata->THREADS; i++)
				{
					thread_soedata_t* t = thread_data + i;

					for (j = 0; j < t->linecount; j++)
					{
						values[retval++] = t->ddata.primes[j];
					}

					free(t->ddata.primes);
					mpz_clear(t->tmpz);
					mpz_clear(t->offset);
					mpz_clear(t->lowlimit);
					mpz_clear(t->highlimit);
				}

				free(thread_data);

				*num_p = retval;
			}
			
            if (sdata->VFLAG > 0)
            {
                printf("found %" PRIu64 " PRPs\n", *num_p);
            }
			
		}

		if (mpz_cmp(offset, lowlimit) != 0)
		{
			// sieving needed to change the lower sieve limit.  adjust the returned
			// values accordingly.
			uint64_t a;
			mpz_sub(tmpz, lowlimit, offset);
			a = mpz_get_ui(tmpz);

            for (i = 0; i < *num_p; i++)
            {
                values[i] -= a;
            }
		}

		// now dump the requested range of primes to a file, or the
		// screen, both, or neither, depending on the state of a couple
		// global configuration variables
		if (PRIMES_TO_FILE)
		{
			FILE *out;
			if (num_witnesses > 0)
				out = fopen("prp_values.dat", "w");
			else
				out = fopen("sieved_values.dat","w");

			if (out == NULL)
			{
				printf("fopen error: %s\n", strerror(errno));
				printf("can't open file for writing\n");
			}
			else
			{
				for (i = 0; i < *num_p; i++)
				{
					//mpz_add_ui(tmpz, *offset, values[i]);
					mpz_add_ui(tmpz, lowlimit, values[i]);
                    if ((mpz_cmp(tmpz, lowlimit) >= 0) && (mpz_cmp(highlimit, tmpz) >= 0))
                    {
                        char* buf = mpz_get_str(NULL, 10, tmpz);
                        fprintf(out, "%s\n", buf);
                        free(buf);
                    }
				}
				fclose(out);
			}
		}

		if (PRIMES_TO_SCREEN)
		{
			for (i = 0; i < *num_p; i++)
			{
				//mpz_add_ui(tmpz, *offset, values[i]);
				mpz_add_ui(tmpz, lowlimit, values[i]);
                if ((mpz_cmp(tmpz, lowlimit) >= 0) && (mpz_cmp(highlimit, tmpz) >= 0))
                {
                    gmp_printf("%Zd\n", tmpz);
                }
			}
			printf("\n");
		}			
	}

	mpz_clear(tmpz);
	mpz_clear(offset);
	//free(offset);

	return values;
}


