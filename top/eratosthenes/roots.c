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

		//sieving requires that we find the offset of each sieve prime in each block 
		//that we sieve.  We are more restricted in choice of offset because we
		//sieve residue classes.  A good way to find the offset is the extended 
		//euclidean algorithm, which reads ax + by = gcd(a,b),
		//where a = prime, b = prodN, and therefore gcd = 1.  
		//since a and b are coprime, y is the multiplicative inverse of prodN modulo prime.
		//This value is a constant, so compute it here in order to facilitate 
		//finding offsets later.

		//solve prodN ^ -1 % p 
		inv = modinv_1(prodN,prime);
		sdata->root[i] = prime - inv;

		//we can also speed things up by computing and storing the residue
		//mod p of the first sieve location in the first residue class.  This provides
		//a speedup by pulling this constant (involving a division) out of a critical loop
		//when finding offsets of bucket sieved primes.
		sdata->lower_mod_prime[i] = (sdata->lowlimit + 1) % prime;
	}

	return;
}
