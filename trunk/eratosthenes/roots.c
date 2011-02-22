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
	int prime, prodN;
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

#endif

	return;
}
