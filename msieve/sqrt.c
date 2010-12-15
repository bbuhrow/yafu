/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 6/27/07

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.
--------------------------------------------------------------------*/

#include "qs.h"

/*--------------------------------------------------------------------*/
uint32 find_factors(fact_obj_t *obj, z *n, 
		fb_element_siqs *factor_base, uint32 fb_size,
		la_col_t *vectors, uint32 vsize, 
		siqs_r *relation_list,
		uint64 *null_vectors, uint32 multiplier,
		z *poly_a_list, poly_t *poly_list,
		factor_list_t *factor_list) {

	/* Perform the square root phase of MPQS. null_vectors
	   contains 64 linear dependencies of the relations in
	   vectors[], which has vsize elements. The code constructs
	   two numbers X and Y (mod n) that are each the product
	   of a huge number of factors derived from vectors[]. It
	   then computes gcd(X+-Y, n) for each linear dependency.
	   gcd values that are not 1 or n are a nontrivial factor
	   of n (or a product of nontrivial factors). For n the
	   product of two prime factors this will happen 2/3 of the 
	   time per dependency, on average.
	   
	   More specifically, vectors[] is a list of relations,
	   each of the form

   	   (a[i] * x[i] + b[i])^2 = prod(factors[i][]) mod n

	   X is a product of the (a[i] * x[i] + b[i]) and Y is
	   the product of sqrt(prod(factors[i][])), both mod n. 
	   The code never needs to calculate explicit square roots, 
	   because it knows the complete factorization of Y. Whether 
	   relation i contributes to dependency j is determined by 
	   the value of the j_th bit of null_vectors[i].

	   Because this implementation uses the double large prime
	   variation of MPQS, a single relation can actually be composed
	   of the product of several entities of the above form. In that
	   case, all entities are multiplied into X and Y (or all
	   are skipped).

	   Note that the code doesn't stop with one nontrivial
	   factor; it prints them all. If you go to so much work
	   and the other dependencies are there for free, why not
	   use them? */

	z factor, x, y, tmp, tmp2, tmp3, sum, tmpn;
	uint32 i, j, k, m;
	uint64 mask;
	uint32 *fb_counts;
	uint32 large_primes[200], num_large_primes;
	uint32 num_relations, prime;
	siqs_r *relation;
	uint32 factor_found = 0;
	int bits;

	zInit(&factor);
	zInit(&x);
	zInit(&y);
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);
	zInit(&sum);
	zInit(&tmpn);
	
	fb_counts = (uint32 *)malloc(fb_size * sizeof(uint32));
	zClear(&factor);
	factor.size = 1;
	zCopy(n, &tmpn);

	bits = 0;
	/* For each dependency */
	for (mask = 1; mask; mask <<= 1) {
		memset(fb_counts, 0, fb_size * sizeof(uint32));
		zClear(&x);
		zClear(&y);
		x.size = y.size = x.val[0] = y.val[0] = 1;
		/* For each sieve relation */

		
		//printf("dependency %d contains relations:\n",dnum);
		for (i = 0; i < vsize; i++) {

			/* If the relation is not scheduled to
			   contribute to x and y, skip it */

			if (!(null_vectors[i] & mask))
				continue;
			
			/* compute the number of sieve_values */

			num_large_primes = 0;
			num_relations = vectors[i].cycle.num_relations;

			/* for all sieve values */

			for (j = 0; j < num_relations; j++) {
				z *a, *b;
				poly_t *poly;
				uint32 sieve_offset;
				uint32 sign_of_index;

				relation = &relation_list[vectors[i].cycle.list[j]];
				
				/* reconstruct a[i], b[i], x[i] and
				   the sign of x[i]. Drop the subscript
				   from here on. */

				poly = poly_list + relation->poly_idx;
				b = &poly->b;
				a = poly_a_list + poly->a_idx;
				//sieve_offset = relation->sieve_offset & 
				//				0x7fffffff;
				//sign_of_index = relation->sieve_offset >> 31;
				sieve_offset = relation->sieve_offset;
				sign_of_index = relation->parity;

				/* Form (a * sieve_offset + b). Note that 
				   sieve_offset can be negative; in that
				   case the minus sign is implicit. We don't
				   have to normalize mod n because there
				   are an even number of negative values
				   to multiply together */
	
				zShortMul(a,sieve_offset,&sum);

				if (sign_of_index == POSITIVE)
				{
					zAdd(&sum,b,&tmp);
					zCopy(&tmp,&sum);
				}
				else
				{
					zSub(&sum,b,&tmp);
					zCopy(&tmp,&sum);
				}
	
				/* multiply the sum into x */
	
				zModMul(&x,&sum,n,&x);
	
				/* do not multiply the factors associated 
				   with this relation into y; instead, just 
				   update the count for each factor base 
				   prime. Unlike ordinary MPQS, the list
				   of factors is for the complete factor-
				   ization of a*(a*x^2+b*x+c), so the 'a' 
				   in front need not be treated separately */
	
				for (k = 0; k < relation->num_factors; k++)
					fb_counts[relation->fb_offsets[k]]++;
					
				/* if the sieve value contains one or more
				   large primes, accumulate them in a 
				   dedicated table. Do not multiply them
				   into y until all of the sieve values
				   for this relation have been processed */

				for (k = 0; k < 2; k++) {
					prime = relation->large_prime[k];
					if (prime == 1)
						continue;

					for (m = 0; m < num_large_primes; m++) {
						if (prime == large_primes[2*m]){
							large_primes[2*m+1]++;
							break;
						}
					}
					if (m == num_large_primes) {
						large_primes[2*m] = prime;
						large_primes[2*m+1] = 1;
						num_large_primes++;
					}
				}
			}

			for (j = 0; j < num_large_primes; j++) {
				for (k = 0; k < large_primes[2*j+1]/2; k++) {
					factor.val[0] = large_primes[2*j];
					zModMul(&y,&factor,n,&y);  
				}
			}
		}

		/* For each factor base prime p, compute 
			p ^ ((number of times p occurs in y) / 2) mod n
		   then multiply it into y. This is enormously
		   more efficient than multiplying by one p at a time */

		for (i = MIN_FB_OFFSET; i < fb_size; i++) {
			uint32 mask2 = 0x80000000;
			uint32 exponent = fb_counts[i] / 2;
			uint32 prime = factor_base->prime[i];

			
			if (fb_counts[i] &0x1)
				printf("odd exponent found\n");
				

			if (exponent == 0)
				continue;

			zClear(&tmp);
			tmp.size = 1;
			tmp.val[0] = factor_base->prime[i]; 
			factor.val[0] = prime;

			while (!(exponent & mask2))
				mask2 >>= 1;
			for (mask2 >>= 1; mask2; mask2 >>= 1) {
				zModMul(&tmp,&tmp,n,&tmp2);
				zCopy(&tmp2,&tmp);
				
				if (exponent & mask2) {
					zModMul(&tmp,&factor,n,&tmp); 
				}
			}
			zModMul(&tmp,&y,n,&y);  
		}

		/* compute gcd(x+y, n). If it's not 1 or n, save it 
		   (and stop processing dependencies if the product 
		   of all the probable prime factors found so far equals 
		   n). See the comments in Pari's MPQS code for a proof 
		   that it isn't necessary to also check gcd(x-y, n) */

		zAdd(&x,&y,&tmp);
		zLEGCD(&tmp,n,&tmp2);
		zCopy(&tmp2,&tmp);
		if (zCompare(&tmp, n) != 0 && !isOne(&tmp)) {

			/* remove any factors of the multiplier 
			   before saving tmp, and don't save at all
			   if tmp contains *only* multiplier factors */
			if (multiplier > 1) {
				uint32 ignore_me = spGCD(multiplier,
						zShortMod(&tmp, multiplier));
				if (ignore_me > 1) {
					zShortDiv(&tmp, ignore_me, &tmp2);
					zCopy(&tmp2,&tmp);
					if (isOne(&tmp))
						continue;
				}
			}

			//ignore composite factors for now...
			if (!isPrime(&tmp))
			{
				continue;
				printf("prp%d = %s\n",ndigits(&tmp),z2decstr(&tmp,&gstr1));
			}

			//add the factor to our global list
			bits = factor_list_add(obj, factor_list, &tmp);

			//check if only the multiplier remains
			if (abs(bits) < 8)
				break;

			//divide the factor out of our number
			zDiv(&tmpn, &tmp, &tmp2, &tmp3);

			//check if the remaining number is prime
			if (isPrime(&tmp2))
			{
				//add it to our global factor list
				//printf("remaining cofactor is prime\n");
				bits = factor_list_add(obj, factor_list, &tmp2);

				//then bail
				break;
			}

			//divide out the multiplier from the remaining number
			if (multiplier > 1) {
				uint32 ignore_me = spGCD(multiplier,
						zShortMod(&tmp2, multiplier));
				if (ignore_me > 1) {
					zShortDiv(&tmp2, ignore_me, &tmp);

					//check again if the remaining number is prime
					if (isPrime(&tmp))
					{
						//add it to our global factor list
						//printf("remaining cofactor is prime after removing multiplier\n");
						bits = factor_list_add(obj, factor_list, &tmp);

						//then bail
						break;
					}
				}
			}

			

			//paranoia
			if (!isZero(&tmp3))
				printf("factor doesn't divide n!\n");


		}
	}

	free(fb_counts);
	zFree(&factor);
	zFree(&x);
	zFree(&y);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&sum);
	return factor_found;
}

