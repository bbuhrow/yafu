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
		if (i == BUCKETSTARTI)
			break;
#else

		inv = modinv_1(prodN,prime);
		
		sdata->root[i] = prime - inv;
		sdata->lower_mod_prime[i] = (sdata->lowlimit + 1) % prime;
#endif
	}

#ifdef INPLACE_BUCKET
	//for each bucket prime, compute its starting block,residue, and sieve location
	//and add it to the appropriate linked list
	j = 0;
	for (i=BUCKETSTARTI; i<sdata->pboundi; i++, j++)
	{
		uint64 root;
		uint32 p = sdata->sieve_p[i];
		uint32 block, loc, residue;

		sdata->bucket_primes[j].primeid = i;
		root = (sdata->lowlimit / (uint64)p + 1) * (uint64)p;
		block = (root - sdata->lowlimit) >> FLAGBITS;
		loc = (root - sdata->lowlimit) & FLAGSIZEm1;
		residue = root % sdata->prodN;
		sdata->bucket_primes[j].loc = loc;
		sdata->bucket_primes[j].res = p % prodN;

		if (sdata->listptrs[block * sdata->prodN + residue] == NULL)
		{
			sdata->listptrs[block * sdata->prodN + residue] = &sdata->bucket_primes[j];
			sdata->bucket_primes[j].next = (uint16)-1;
		}
		else
		{
			bucket_prime_t *lastp = sdata->listptrs[block * sdata->prodN + residue];
			sdata->bucket_primes[j].next = i - lastp->primeid;
			sdata->listptrs[block * sdata->prodN + residue] = &sdata->bucket_primes[j];
		}
	}

#endif

	return;
}
