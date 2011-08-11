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

void GetPRIMESRange(uint64 lowlimit, uint64 highlimit)
{
	uint64 i;
	
	//reallocate array based on conservative estimate of the number of 
	//primes in the interval
	
	i = (uint64)(highlimit/log(highlimit));
	if (lowlimit > 1)
		i -= (uint64)(lowlimit/log(lowlimit));
	i += (highlimit - lowlimit) * 1.2;

	PRIMES = (uint64 *)realloc(PRIMES,(size_t) (i*sizeof(uint64)));
	if (PRIMES == NULL)
	{
		printf("unable to allocate %" PRIu64 " bytes for range %" PRIu64 " to %" PRIu64 "\n",
			(uint64)(i*sizeof(uint64)),lowlimit,highlimit);
		exit(1);
	}

	//reset the global constants
	P_MIN = lowlimit;
	P_MAX = highlimit; 

	//find the primes in the interval
	NUM_P = spSOE(PRIMES,lowlimit,&highlimit,0);

	return;
}

uint64 soe_wrapper(uint64 lowlimit, uint64 highlimit, int count)
{
	//public interface to the sieve.  necessary because in order to keep the 
	//sieve efficient it must sieve larger blocks of numbers than a user may want,
	//and because the program keeps a cache of primes on hand which may or may 
	//not contain the range of interest.  Manage this on-hand cache and any addition
	//sieving needed.
	//TODO: manage really large requested blocks by splitting the range up into
	//blocks of no more than 10B else memory requirements get outragous.
	uint64 retval, tmpl, tmph,i=0;
	uint64 maxrange = 1000000000000ULL;

	if (highlimit < lowlimit)
	{
		printf("error: lowlimit must be less than highlimit\n");
		return 0;
	}

	if (count)
	{
		//this needs to be a range of at least 1e6
		if ((highlimit - lowlimit) < 1000000)
		{
			//maybe it is already in our list of cached primes
			if ((lowlimit >= P_MIN) && (highlimit <= P_MAX))
			{
				retval = 0;
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						retval++;
				}
			}
			else
			{
				//nope, go and get a new range.
				tmpl = lowlimit;
				tmph = tmpl + 1000000;

				//since this is a small range, we need to 
				//find a bigger range and count them.
				GetPRIMESRange(tmpl,tmph);
				retval = 0;
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						retval++;
				}
			}
		}
		else
		{
			//check for really big ranges
			if ((highlimit - lowlimit) > maxrange)
			{
				uint32 num_ranges = (uint32)((highlimit - lowlimit) / maxrange);
				uint64 remainder = (highlimit - lowlimit) % maxrange;
				uint32 j;
				
				retval = 0;
				tmpl = lowlimit;
				tmph = lowlimit + maxrange;
				for (j = 0; j < num_ranges; j++)
				{
					retval += spSOE(NULL,tmpl,&tmph,1);
					tmpl += maxrange;
					tmph = tmpl + maxrange;
				}
				
				tmph = tmpl + remainder;
				retval += spSOE(NULL,tmpl,&tmph,1);
			}
			else
			{
				//we're in a sweet spot already, just get the requested range
				retval = spSOE(NULL,lowlimit,&highlimit,1);
			}
		}

	}
	else
	{
		if (lowlimit < P_MIN || lowlimit > P_MAX || highlimit > P_MAX)
		{
			//requested range is not covered by the current range
			tmpl = lowlimit;
			tmph = highlimit;

			//this needs to be a range of at least 1e6
			if (tmph - tmpl < 1000000)
			{
				//there is slack built into the sieve limit, so go ahead and increase
				//the size of the interval to make it at least 1e6.
				tmph = tmpl + 1000000;

				//since this is a small range, we need to 
				//find a bigger range and count them.
				GetPRIMESRange(tmpl,tmph);
				retval = 0;
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						retval++;
				}
			}
			else
			{
				//we don't need to mess with the requested range,
				//so GetPRIMESRange will return the requested range directly
				//and the count will be in NUM_P
				GetPRIMESRange(lowlimit,highlimit);
				retval = NUM_P;
			}
		}
		else
		{
			// the requested range is covered by the current range
			// just count them
			retval = 0;
			for (i = 0; i < NUM_P; i++)
			{
				if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
					retval++;
			}
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
				printf("can't open primes.dat for writing\n");
			}
			else
			{
				for (i = 0; i < NUM_P; i++)
				{
					if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
						fprintf(out,"%" PRIu64 "\n",PRIMES[i]);
				}
				fclose(out);
			}
		}

		if (PRIMES_TO_SCREEN)
		{
			for (i = 0; i < NUM_P; i++)
			{
				if (PRIMES[i] >= lowlimit && PRIMES[i] <= highlimit)
					printf("%" PRIu64 " ",PRIMES[i]);
			}
			printf("\n");
		}

			
	}

	return retval;
}