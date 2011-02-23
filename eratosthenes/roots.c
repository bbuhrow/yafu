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

void getRoots(soe_staticdata_t *sdata)
{
	int prime, prodN, j;
	uint64 startprime;
	uint64 i;

	prodN = (int)sdata->prodN;
	startprime = sdata->startprime;

	for (i=startprime;i<sdata->pboundi;i++)
	{
		uint32 inv;
		prime = sdata->sieve_p[i];

#ifdef INPLACE_BUCKET
		if (i == sdata->inplace_startindex)
			break;
#endif

		inv = modinv_1(prodN,prime);
		
		sdata->root[i] = prime - inv;
		sdata->lower_mod_prime[i] = (sdata->lowlimit + 1) % prime;
	}

#ifdef INPLACE_BUCKET
	//for each bucket prime, compute its starting block,residue, and sieve location
	//and add it to the appropriate linked list
	j = 0;
	for (i=sdata->inplace_startindex; i<sdata->pboundi; i++, j++)
	{
		uint64 root;
		uint32 p = sdata->sieve_p[i];
		uint32 block, loc, residue, steps;

		sdata->bucket_primes[j].prime = p;
		sdata->bucket_primes[j].id = j;

		//compute the location of the first hit on the number line 
		//that is past the beginning of our interval
		root = (sdata->lowlimit / (uint64)p + 1) * (uint64)p;				

		//translate this to a hit on one of the residue lines within the interval
		//begin by finding the residue class of the hit
		residue = root % sdata->prodN;

		//then find the number of steps on the residue line -- this is the bit offset
		//of the hit in this residue line
		steps = (root - sdata->lowlimit - residue) / sdata->prodN;									

		//to prevent having to do the divisions every time we step this prime on the
		//number line, record the number of steps and the residual step that this
		//prime takes in residue space.  we need to know the residue class of the prime 
		//itself anyway.
		sdata->bucket_primes[j].steps = p / prodN;
		sdata->bucket_primes[j].res = p % prodN;

		//now that we have all this info, iterate the first hit until it is on
		//a residue line that we are sieving
		while (!sdata->valid_residue[residue])
		{
			steps += sdata->bucket_primes[j].steps;
			residue += sdata->bucket_primes[j].res;
			if (residue >= sdata->prodN)
			{
				residue -= sdata->prodN;
				steps++;
			}
		}					

		//translate the final hit to block,offset notation
		block = steps >> FLAGBITS;
		sdata->bucket_primes[j].loc = steps & FLAGSIZEm1;

		//debug: check for divisibility
		root = sdata->lowlimit + residue + (uint64)((block << FLAGBITS) + sdata->bucket_primes[j].loc) * (uint64)prodN;
		////root = sdata->lowlimit + residue + (uint64)steps * (uint64)prodN;
		if (root % p != 0)
			printf("prime %u does not divide this location\n",p);

		if (block >= sdata->blocks)
			continue;

		if (sdata->listptrs[residue * sdata->blocks + block] == NULL)
		{
			sdata->listptrs[residue * sdata->blocks + block] = &sdata->bucket_primes[j];
			sdata->bucket_primes[j].next = (uint32)-1;
		}
		else
		{			
			bucket_prime_t *lastp = sdata->listptrs[residue * sdata->blocks + block];
			sdata->bucket_primes[j].next = lastp->id;
			sdata->listptrs[residue * sdata->blocks + block] = &sdata->bucket_primes[j];
		}
	}


#endif

	return;
}
