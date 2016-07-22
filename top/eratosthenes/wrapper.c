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

uint64 *GetPRIMESRange(uint32 *sieve_p, uint32 num_sp, 
	mpz_t *offset, uint64 lowlimit, uint64 highlimit, uint64 *num_p)
{
	uint64 i;
	uint64 hi_est, lo_est;
	uint64 maxrange = 10000000000ULL;
	uint64 *primes = NULL;
	
	//reallocate output array based on conservative estimate of the number of 
	//primes in the interval
	if (offset != NULL)
	{
		i = (highlimit - lowlimit);
		primes = (uint64 *)realloc(primes, (size_t) (i * sizeof(uint64)));
		if (primes == NULL)
		{
			if (offset == NULL)
				printf("unable to allocate %" PRIu64 " bytes for range %" PRIu64 " to %" PRIu64 "\n",
					(uint64)(i * sizeof(uint64)),lowlimit,highlimit);
			else
				printf("unable to allocate %" PRIu64 " bytes \n",
					(uint64)(i * sizeof(uint64)));
			exit(1);
		}
	}
	else
	{
		hi_est = (uint64)(highlimit/log((double)highlimit));
		if (lowlimit > 1)
			lo_est = (uint64)(lowlimit/log((double)lowlimit));
		else
			lo_est = 0;

		i = (uint64)((double)(hi_est - lo_est) * 1.25);

		if (1) //(!NO_STORE)
		{
			primes = (uint64 *)realloc(primes,(size_t) (i * sizeof(uint64)));
			if (primes == NULL)
			{
				printf("unable to allocate %" PRIu64 " bytes for range %" PRIu64 " to %" PRIu64 "\n",
					(uint64)(i * sizeof(uint64)),lowlimit,highlimit);
				exit(1);
			}
		}
	}

	//check for really big ranges ('big' is different here than when we are counting
	//primes because there are higher memory demands when computing primes)
	if ((highlimit - lowlimit) > maxrange)
	{
		uint64 tmpl, tmph, tmpcount = 0;
		uint32 num_ranges = (uint32)((highlimit - lowlimit) / maxrange);
		uint64 remainder = (highlimit - lowlimit) % maxrange;
		uint32 j;
				
		GLOBAL_OFFSET = 0;
		tmpl = lowlimit;
		tmph = lowlimit + maxrange;
		for (j = 0; j < num_ranges; j++)
		{
			tmpcount += spSOE(sieve_p, num_sp, offset, tmpl, &tmph, 0, primes);
			tmpl += maxrange;
			tmph = tmpl + maxrange;
			GLOBAL_OFFSET = tmpcount;
		}
				
		tmph = tmpl + remainder;
		tmpcount += spSOE(sieve_p, num_sp, offset, tmpl, &tmph, 0, primes);
		*num_p = tmpcount;
	}
	else
	{
		//find the primes in the interval
		GLOBAL_OFFSET = 0;
		*num_p = spSOE(sieve_p, num_sp, offset, lowlimit, &highlimit, 0, primes);
	}

	return primes;
}

uint64 *soe_wrapper(uint32 *seed_p, uint32 num_sp, 
	uint64 lowlimit, uint64 highlimit, int count, uint64 *num_p)
{
	//public interface to the sieve.  
	uint64 retval, tmpl, tmph, i;
	uint32 max_p;	
	uint32 *sieve_p;
	uint64 *primes = NULL;

	if (highlimit < lowlimit)
	{
		printf("error: lowlimit must be less than highlimit\n");
		*num_p = 0;
		return primes;
	}	

	if (highlimit > (seed_p[num_sp-1] * seed_p[num_sp-1]))
	{
		//then we need to generate more sieving primes
		uint32 range_est;
	
		//allocate array based on conservative estimate of the number of 
		//primes in the interval	
		max_p = (uint32)sqrt((int64)(highlimit)) + 65536;
		range_est = (uint32)estimate_primes_in_range(0, (uint64)max_p);
		sieve_p = (uint32 *)xmalloc_align((size_t) (range_est * sizeof(uint32)));

		if (sieve_p == NULL)
		{
			printf("unable to allocate %u bytes for %u sieving primes\n",
				range_est * (uint32)sizeof(uint32), range_est);
			exit(1);
		}

		//find the sieving primes using the seed primes
		NO_STORE = 0;
		primes = GetPRIMESRange(seed_p, num_sp, NULL, 0, max_p, &retval);
		for (i=0; i<retval; i++)
			sieve_p[i] = (uint32)primes[i];
		printf("found %u sieving primes\n",(uint32)retval);
		num_sp = (uint32)retval;
		free(primes);
		primes = NULL;
		//NO_STORE = 1;
	}
	else
	{
		//seed primes are enough
        sieve_p = (uint32 *)xmalloc_align((size_t)(num_sp * sizeof(uint32)));
		//NO_STORE = 1;

		if (sieve_p == NULL)
		{
			printf("unable to allocate %u bytes for %u sieving primes\n",
				num_sp * (uint32)sizeof(uint32), num_sp);
			exit(1);
		}

		for (i=0; i<num_sp; i++)
			sieve_p[i] = seed_p[i];
	}

	if (count)
	{
		//this needs to be a range of at least 1e6
		if ((highlimit - lowlimit) < 1000000)
		{
			//go and get a new range.
			tmpl = lowlimit;
			tmph = tmpl + 1000000;

			//since this is a small range, we need to 
			//find a bigger range and count them.
			primes = GetPRIMESRange(sieve_p, num_sp, NULL, tmpl, tmph, &retval);

			*num_p = 0;
			//count how many are in the original range of interest
			for (i = 0; i < retval; i++)
			{
				if (primes[i] >= lowlimit && primes[i] <= highlimit)
					(*num_p)++;
			}
			free(primes);
			primes = NULL;
		}
		else
		{
			//check for really big ranges
			uint64 maxrange = 100000000000ULL;

			if ((highlimit - lowlimit) > maxrange)
			{
				uint32 num_ranges = (uint32)((highlimit - lowlimit) / maxrange);
				uint64 remainder = (highlimit - lowlimit) % maxrange;
				uint32 j;
				//to get time per range
				double t_time;
				struct timeval start, stop;
				TIME_DIFF *	difference;
				
				*num_p = 0;
				tmpl = lowlimit;
				tmph = lowlimit + maxrange;
				gettimeofday (&start, NULL);

				for (j = 0; j < num_ranges; j++)
				{
					*num_p += spSOE(sieve_p, num_sp, NULL, tmpl, &tmph, 1, NULL);

					gettimeofday (&stop, NULL);
					difference = my_difftime (&start, &stop);

					t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
					free(difference);

					if (VFLAG > 1)
						printf("so far, found %" PRIu64 " primes in %1.1f seconds\n",*num_p, t_time);
					tmpl += maxrange;
					tmph = tmpl + maxrange;
				}
				
				if (remainder > 0)
				{
					tmph = tmpl + remainder;
					*num_p += spSOE(sieve_p, num_sp, NULL, tmpl, &tmph, 1, NULL);
				}
				if (VFLAG > 1)
					printf("so far, found %" PRIu64 " primes\n",*num_p);
			}
			else
			{
				//we're in a sweet spot already, just get the requested range
				*num_p = spSOE(sieve_p, num_sp, NULL, lowlimit, &highlimit, 1, NULL);
			}
		}

	}
	else
	{
		tmpl = lowlimit;
		tmph = highlimit;

		//this needs to be a range of at least 1e6
		if ((tmph - tmpl) < 1000000)
		{
			//there is slack built into the sieve limit, so go ahead and increase
			//the size of the interval to make it at least 1e6.
			tmph = tmpl + 1000000;

			//since this is a small range, we need to 
			//find a bigger range and count them.
			primes = GetPRIMESRange(sieve_p, num_sp, NULL, tmpl, tmph, &retval);
			*num_p = 0;
			for (i = 0; i < retval; i++)
			{
				if (primes[i] >= lowlimit && primes[i] <= highlimit)
					(*num_p)++;
			}

		}
		else
		{
			//we don't need to mess with the requested range,
			//so GetPRIMESRange will return the requested range directly
			//and the count will be in NUM_P
			primes = GetPRIMESRange(sieve_p, num_sp, NULL, lowlimit, highlimit, num_p);
		}

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
					if (primes[i] >= lowlimit && primes[i] <= highlimit)
						fprintf(out,"%" PRIu64 "\n",primes[i]);
				}
				fclose(out);
			}
		}

		if (PRIMES_TO_SCREEN)
		{
			for (i = 0; i < *num_p; i++)
			{
				if (primes[i] >= lowlimit && primes[i] <= highlimit)
					printf("%" PRIu64 " ",primes[i]);
			}
			printf("\n");
		}			
	}

	align_free(sieve_p);
	return primes;
}

uint64 *sieve_to_depth(uint32 *seed_p, uint32 num_sp, 
	mpz_t lowlimit, mpz_t highlimit, int count, int num_witnesses, uint64 *num_p)
{
	//public interface to a routine which will sieve a range of integers
	//with the supplied primes and either count or compute the values
	//that survive.  Basically, it is just the sieve, but with no
	//guareentees that what survives the sieving is prime.  The idea is to 
	//remove cheap composites.
	uint64 retval, i, range, tmpl, tmph;
	uint64 *values = NULL;
	mpz_t tmpz;
	mpz_t *offset;

	if (mpz_cmp(highlimit, lowlimit) <= 0)
	{
		printf("error: lowlimit must be less than highlimit\n");
		*num_p = 0;
		return values;
	}	

	offset = (mpz_t *)malloc(sizeof(mpz_t));
	mpz_init(tmpz);
	mpz_init(*offset);
	mpz_set(*offset, lowlimit);
	mpz_sub(tmpz, highlimit, lowlimit);
	range = mpz_get_64(tmpz);

	if (count)
	{
		//this needs to be a range of at least 1e6
		if (range < 1000000)
		{
			//go and get a new range.
			tmpl = 0;
			tmph = 1000000;

			//since this is a small range, we need to 
			//find a bigger range and count them.
			values = GetPRIMESRange(seed_p, num_sp, offset, tmpl, tmph, &retval);

			*num_p = 0;
			//count how many are in the original range of interest
			for (i = 0; i < retval; i++)
			{
				mpz_add_ui(tmpz, *offset, values[i]);
				if ((mpz_cmp(tmpz, lowlimit) >= 0) && (mpz_cmp(highlimit, tmpz) >= 0))
					(*num_p)++;
			}
			free(values);
			values = NULL;
		}
		else
		{
			//check for really big ranges
			uint64 maxrange = 100000000000ULL;

			if (range > maxrange)
			{
				uint32 num_ranges = (uint32)(range / maxrange);
				uint64 remainder = range % maxrange;
				uint32 j;
				
				*num_p = 0;
				tmpl = 0;
				tmph = tmpl + maxrange;
				for (j = 0; j < num_ranges; j++)
				{
					*num_p += spSOE(seed_p, num_sp, offset, tmpl, &tmph, 1, NULL);
					if (VFLAG > 1)
						printf("so far, found %" PRIu64 " primes\n",*num_p);
					tmpl += maxrange;
					tmph = tmpl + maxrange;
				}
				
				if (remainder > 0)
				{
					tmph = tmpl + remainder;
					*num_p += spSOE(seed_p, num_sp, offset, tmpl, &tmph, 1, NULL);
				}
				if (VFLAG > 1)
					printf("so far, found %" PRIu64 " primes\n",*num_p);
			}
			else
			{
				//we're in a sweet spot already, just get the requested range
				*num_p = spSOE(seed_p, num_sp, offset, 0, &range, 1, NULL);
			}
		}

	}
	else
	{
		//this needs to be a range of at least 1e6
		if (range < 1000000)
		{
			//there is slack built into the sieve limit, so go ahead and increase
			//the size of the interval to make it at least 1e6.
			tmpl = 0;
			tmph = tmpl + 1000000;

			//since this is a small range, we need to 
			//find a bigger range and count them.
			values = GetPRIMESRange(seed_p, num_sp, offset, tmpl, tmph, &retval);
			*num_p = 0;
			for (i = 0; i < retval; i++)
			{
				mpz_add_ui(tmpz, *offset, values[i]);
				if ((mpz_cmp(tmpz, lowlimit) >= 0) && (mpz_cmp(highlimit, tmpz) >= 0))
					(*num_p)++;
			}

		}
		else
		{
			//we don't need to mess with the requested range,
			//so GetPRIMESRange will return the requested range directly
			//and the count will be in NUM_P
			values = GetPRIMESRange(seed_p, num_sp, offset, 0, range, num_p);

		}

		if (num_witnesses > 0)
		{
			int pchar = 0;
			thread_soedata_t *thread_data;		//an array of thread data objects
			uint32 lastid;
			int j;

			//allocate thread data structure
			thread_data = (thread_soedata_t *)malloc(THREADS * sizeof(thread_soedata_t));
			
			// conduct PRP tests on all surviving values
			if (VFLAG > 0)
				printf("starting PRP tests with %d witnesses on %" PRIu64 " surviving candidates\n", 
					num_witnesses, *num_p);

			// start the threads
			for (i = 0; i < THREADS - 1; i++)
				start_soe_worker_thread(thread_data + i);

			//start_soe_worker_thread(thread_data + i, 1);

			range = *num_p / THREADS;
			lastid = 0;

			// divvy up the range
			for (j = 0; j < THREADS; j++)
			{
				thread_soedata_t *t = thread_data + j;
				
				t->startid = lastid;
				t->stopid = t->startid + range;
				lastid = t->stopid;

				if (VFLAG > 2)
					printf("thread %d computing PRPs from %u to %u\n", 
						(int)i, t->startid, t->stopid);
			}

			// the last one gets any leftover
			if (thread_data[THREADS-1].stopid != (uint32)*num_p)
				thread_data[THREADS-1].stopid = (uint32)*num_p;

			// allocate space for stuff in the threads
			if (THREADS == 1)
			{
				thread_data[0].ddata.primes = values;
			}
			else
			{
				for (j = 0; j < THREADS; j++)
				{
					thread_soedata_t *t = thread_data + j;

					mpz_init(t->tmpz);
					mpz_init(t->offset);
					mpz_init(t->lowlimit);
					mpz_init(t->highlimit);
					mpz_set(t->offset, *offset);
					mpz_set(t->lowlimit, lowlimit);
					mpz_set(t->highlimit, highlimit);
					t->current_line = (uint64)num_witnesses;

					t->ddata.primes = (uint64 *)malloc((t->stopid - t->startid) * sizeof(uint64));
					for (i = t->startid; i < t->stopid; i++)
						t->ddata.primes[i - t->startid] = values[i];
				}
			}

			// now run with the threads
			for (j = 0; j < THREADS; j++)
			{
				thread_soedata_t *t = thread_data + j;

				if (j == (THREADS - 1)) 
				{	
					t->linecount = 0;
					for (i = t->startid; i < t->stopid; i++)
					{
						if (((i & 128) == 0) && (VFLAG > 0))
						{
							int k;
							for (k = 0; k<pchar; k++)
								printf("\b");
							pchar = printf("progress: %d%%",(int)((double)i / (double)(*num_p) * 100.0));
							fflush(stdout);
						}

						mpz_add_ui(tmpz, *offset, t->ddata.primes[i - t->startid]);
						if ((mpz_cmp(tmpz, lowlimit) >= 0) && (mpz_cmp(highlimit, tmpz) >= 0))
						{
							if (is_mpz_prp(tmpz))
								t->ddata.primes[t->linecount++] = t->ddata.primes[i - t->startid];
						}
					}
				}
				else
				{
					t->command = SOE_COMPUTE_PRPS;

#if defined(WIN32) || defined(_WIN64)
					SetEvent(t->run_event);
#else
					pthread_cond_signal(&t->run_cond);
					pthread_mutex_unlock(&t->run_lock);
#endif

				}
			}

			//wait for each thread to finish
			for (i = 0; i < THREADS; i++) 
			{
				thread_soedata_t *t = thread_data + i;

				if (i < (THREADS - 1)) 
				{
#if defined(WIN32) || defined(_WIN64)
					WaitForSingleObject(t->finish_event, INFINITE);
#else
					pthread_mutex_lock(&t->run_lock);
					while (t->command != SOE_COMMAND_WAIT)
						pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
				}
			}

			//stop the worker threads
			for (i=0; i<THREADS - 1; i++)
				stop_soe_worker_thread(thread_data + i);

			// combine results and free stuff
			if (THREADS == 1)
			{
				retval = thread_data[0].linecount;
			}
			else
			{
				retval = 0;
				for (i=0; i<THREADS; i++)
				{
					thread_soedata_t *t = thread_data + i;

					for (j=0; j < t->linecount; j++)
						values[retval++] = t->ddata.primes[j];

					free(t->ddata.primes);
					mpz_clear(t->tmpz);
					mpz_clear(t->offset);
					mpz_clear(t->lowlimit);
					mpz_clear(t->highlimit);
				}
			}

			free(thread_data);

			if (VFLAG > 0)
			{
				int k;
				for (k = 0; k<pchar; k++)
					printf("\b");
			}

			*num_p = retval;
			if (VFLAG > 0)
				printf("found %" PRIu64 " PRPs\n", *num_p);
			
		}

		if (mpz_cmp(*offset, lowlimit) != 0)
		{
			// sieving needed to change the lower sieve limit.  adjust the returned
			// values accordingly.
			uint64 a;
			mpz_sub(tmpz, lowlimit, *offset);
			a = mpz_get_64(tmpz);

			for (i=0; i < *num_p; i++)
				values[i] -= a;
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
						gmp_fprintf(out,"%Zd\n",tmpz);
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
					gmp_printf("%Zd\n",tmpz);
			}
			printf("\n");
		}			
	}

	mpz_clear(tmpz);
	mpz_clear(*offset);
	free(offset);

	return values;
}


