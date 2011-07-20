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

void primes_from_lineflags(thread_soedata_t *t)
{
	//the sieve primes are not in the line array, so they must be added
	//in if necessary
	soe_staticdata_t *sdata = &t->sdata;
	soe_dynamicdata_t *ddata = &t->ddata;
	uint8 *line = t->ddata.line;
	uint64 current_line = t->current_line;
	uint64 prime, num_alloc;
	uint64 i,j,it;

	it=0;

	count_line(t);
	num_alloc = t->linecount;
	ddata->primes = (uint64 *)malloc(num_alloc * sizeof(uint64));
	if (ddata->primes == NULL)
	{
		printf("failed to allocate primes array in primes_from_lineflags\n");
		exit(-1);
	}

	//this will find all the primes in the line.  when other lines are appended, the primes
	//will be out of order and will need to be sorted.
	for (i=0;i<sdata->numlinebytes;i++)
	{
		for (j=0;j<BITSINBYTE;j++)
		{
			if (line[i] & nmasks[j])
			{
				prime = sdata->prodN * (i*BITSINBYTE + j) + 
					sdata->rclass[current_line] + sdata->lowlimit;

				//only store the prime if it is within our requested bounds
				if ((prime >= sdata->orig_llimit) &&  (prime <= sdata->orig_hlimit))
				{
					ddata->primes[it] = prime;
					it++;
				}
				//else
					//printf("prime %" PRIu64 " is outside of original limits of %" PRIu64 " and %" PRIu64"\n",
					//prime, sdata->orig_llimit, sdata->orig_hlimit);
			}
			
		}
	}

	if (it != t->linecount)
		printf("warning, counts do not match after computing primes\n");

	t->linecount = it;

	return;
}


