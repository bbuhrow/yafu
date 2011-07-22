/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"
#include "util.h"
#include "common.h"
#include "gmp_xface.h"

//#define POLYA_DEBUG

void new_poly_a(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	/*the goal of this routine is to generate a new poly_a value from elements of the factor base
	subject to a few constraints.  first, the number of fb elements used should always be greater than
	3, and should grow based on the size of a.  second, the elements each should be greater than 2000, 
	to prevent relation redundancy and to prevent decrease the probability of finding smooth relations
	(because we don't sieve with the primes making up 'a').  third, the elements should be as small
	as possible, subject to condition 2.  fourth, the elements making up each 'a' should be different by 
	at least 2 element from every other 'a', to prevent relation redundancy.
		
	to start, determine approximately how many elements can be used given the target 'a', say this 
	number is s.  pick s-1 elements from the factor base such that the last element will need to be larger
	than any of the others.  then choose the best last value so the the actual 'a' is as close as
	possible to the target 'a'.

	when picking the elements, randomly pick the first s-1 elements from a pool of fb elements, then tailor 
	the last one.  when done, compare to all previous 'a' elements, and change individual values of the 
	new one so that it is sufficiently different from all others.  i assume the permutations of element
	choices will essentially never run dry.  this seems reasonable.
	*/
	
	/*
	try it like this:
	set the pool of elements to be the primes between 500 and 1500 (average 1000).  on average there are
	about 70 such primes.  we will pick the first s-1 primes from this pool, and the last one will be higher than
	1500.  

	here are estimates for the number of elements given digits in n
	ndigits		adigits		elements
	40			14			5
	50			19			7
	60			24			9
	70			29			11
	80			34			13
	etc.

	pick the first s-1 digits randomly from the pool.  with n=50 digits, s-1 = 6.  there are about 130e6 
	combinations in picking 6 elements out of 70, so the list should not run dry.  pick the last element out
	of the general factor base, but be sure that it has an index less than small_B.  compare to previous 
	element choices for 'a' and redo if too similar.  with that many combinations, I'm betting re-doing won't 
	happen often.
	*/

	//unpack stuff from the job data structure
	siqs_poly *poly = dconf->curr_poly;
	z *target_a = &sconf->target_a;
	fb_list *fb = sconf->factor_base;

	z tmp, tmp2, tmp3, *poly_a = &poly->poly_a;
	int j, *qli = poly->qlisort, *s = &poly->s;
	uint32 i,randindex = 0, mindiff,a1,poly_low_found=0,target_bits;
	uint32 potential_a_factor = 0, found_a_factor;
	uint32 afact[20];
	double target_mul = 0.9;
	int too_close, min_ratio;
	FILE *sieve_log = sconf->obj->logfile;
	uint32 upper_polypool_index, lower_polypool_index;

	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);

	//determine polypool indexes.  
	//this really should be done once after generating the factor base
	//these will be set to more appropriate values below
	lower_polypool_index = 2;
	upper_polypool_index = fb->small_B - 1;

	if (sconf->bits < 130)
	{
		// don't worry so much about generating a poly close to the target,
		// just make sure the factors of poly_a are all relatively large and
		// disperse to keep duplicates low.
		uint32 hi = fb->B-1;
		uint32 lo = (fb->B-1) / 4;
		uint32 id;
		int doover;

		id = lo + (uint32)((double)(hi - lo) * (double)rand() / (double)RAND_MAX);
		sp2z(fb->list->prime[id],poly_a);
		qli[0] = id;

		j=1;
		while (zCompare(poly_a,target_a) < 0)
		{
			int k = 0;
			id = lo + (uint32)((double)(hi - lo) * (double)rand() / (double)RAND_MAX);

			doover = 0;
			for (k=0; k < j; k++)
			{
				if (id == qli[k])
				{
					doover = 1;
					break;
				}
			}

			if (doover)
				continue;

			zShortMul(poly_a,fb->list->prime[id],poly_a);
			qli[j++] = id;
		}

		*s = j;

#ifdef POLYA_DEBUG
			printf("A id/factors = %d:%u, %d:%u, %d:%u, A = %s\n", 
				qli[0],fb->list->prime[qli[0]], 
				qli[1],fb->list->prime[qli[1]], 
				qli[2],fb->list->prime[qli[2]], 
				z2decstr(poly_a,&gstr1));
#endif

		goto done;

	}
	else
	{
		for (i=0;i<fb->small_B;i++)
		{
			if ((fb->list->prime[i] > 1000) && (poly_low_found == 0))
			{
				lower_polypool_index = i;
				poly_low_found=1;
			}

			if (fb->list->prime[i] > 4000)
			{
				upper_polypool_index = i-1;
				break;
			}
		}
		upper_polypool_index = fb->small_B - 1;

		//brute force the poly to be somewhat close to the target
		target_bits = (uint32)((double)zBits(target_a) * target_mul);
		too_close = 10;
		min_ratio = 1000;
	}


	while (1)
	{
		//generate poly_a's until the residue is 'small enough'
#ifdef POLYA_DEBUG
		printf("*******new trial a********\n");
#endif

		//sp2z(1,poly_a);
		zCopy(&zOne,poly_a);
		*s=0;
		for (;;)
		{
			//randomly pick a new unique factor
			found_a_factor = 0;
			while (!found_a_factor)
			{
				randindex = (uint32)spRand((fp_digit)lower_polypool_index,
					(fp_digit)upper_polypool_index);
				//randindex = lower_polypool_index + 
				//	(uint32)((upper_polypool_index-lower_polypool_index) * (double)rand() / (double)RAND_MAX);
				potential_a_factor = fb->list->prime[randindex];
				//make sure we haven't already randomly picked this one
				found_a_factor = 1;
				for (j=0;j<*s;j++)
				{
					if (afact[j] == potential_a_factor)
					{
						found_a_factor = 0;
						break;
					}
				}
			}
			
			//build up poly_a
			zShortMul(poly_a,potential_a_factor,poly_a);
#ifdef POLYA_DEBUG
			printf("afactor %d = %u\n",*s,potential_a_factor);
#endif
			afact[*s]=potential_a_factor;
			qli[*s] = randindex;
			*s = *s + 1;
			//compute how close we are to target_a
			j = zBits(target_a) - zBits(poly_a);
			if (j < too_close)
			{
				//too close, we want the last factor to be between 15 and 10 bits
#ifdef POLYA_DEBUG
				printf("target_a too close for last factor\n");
#endif
				zCopy(&zOne,poly_a);
				*s=0;
				continue;
			}
			else if (j < (too_close + 5))
			{
				//close enough to pick a last factor
#ifdef POLYA_DEBUG
				printf("picking last factor\n");
#endif
				break;
			}
		}

		//at this point, poly_a is too small by one factor, find the closest factor
		zCopy(target_a,&tmp);
		zDiv(&tmp,poly_a,&tmp2,&tmp3);

		mindiff = 0xffffffff;
		a1 = tmp2.val[0];
		if (a1 < min_ratio)
		{
#ifdef POLYA_DEBUG
			printf("ratio = %u, starting over\n",a1);
#endif
			continue;
		}

		randindex = 0;
		for (i=0;i<fb->small_B;i++)
		{
			if ((uint32)abs(a1 - fb->list->prime[i]) < mindiff)
			{
				mindiff = abs(a1 - fb->list->prime[i]);
				randindex = i;
			}
		}
		//randindex should be the index of the best prime
		//check to make sure it's unique
		found_a_factor = 0;
		do
		{
			potential_a_factor = fb->list->prime[randindex];
			//make sure we haven't already randomly picked this one
			found_a_factor = 1;
			for (j=0;j<*s;j++)
			{
				if (afact[j] == potential_a_factor)
				{
					found_a_factor = 0;
					break;
				}
			}
			if (!found_a_factor)
			{
				//this one is taken.  for now, just try the next bigger one
				randindex++;
			}
		} while (!found_a_factor);

		if (randindex > fb->small_B)
		{
#ifdef POLYA_DEBUG
			printf("last prime in poly_a > small_B\n");
#endif
			continue;
		}

		zShortMul(poly_a,fb->list->prime[randindex],poly_a);
#ifdef POLYA_DEBUG
		printf("afactor %d = %u\n",*s,fb->list[randindex].prime);
#endif
		afact[*s] = fb->list->prime[randindex];
		qli[*s] = randindex;
		*s = *s + 1;

		//check if 'close enough'
		zSub(target_a,poly_a,&tmp);

		if ((uint32)zBits(&tmp) < target_bits)
		{ 
			// if not a duplicate
			found_a_factor = 0;
			for (j=0; j< (int)sconf->total_poly_a; j++)
			{
				if (zCompare(poly_a,&sconf->poly_a_list[j]) == 0)
				{
					found_a_factor = 1;
					break;
				}
			}

			if (found_a_factor)
			{
				//increase the target bound, so it is easier to find a factor.
				//very rarely, inputs seem to generate many duplicates, and
				//in that case we make it easier to find a non-duplicate
				if (target_bits > 1000)
				{
					printf("running away.  POLYPOOL bounds were: %u to %u (%d primes)\nkilling... \n",
						fb->list->prime[lower_polypool_index],fb->list->prime[upper_polypool_index],
						upper_polypool_index - lower_polypool_index);
					exit(-1);
				}

				target_bits++;
				printf("poly %s is a duplicate of #%d = %s\n",
					z2decstr(poly_a,&gstr1),j,z2decstr(&sconf->poly_a_list[j],&gstr2));
				printf("rejecting duplicate poly_a, new target = %d\n",target_bits);
				printf("primes in a: ");
				for (i=0;i<*s;i++)
					printf("%u, ",fb->list->prime[qli[i]]);
				printf("\n");
				logprint(sieve_log,"rejecting duplicate poly_a, new target = %d\n",target_bits);
				continue;
			}
			else break;
			
			//check that this poly has at least 2 factors different from all
			//previous polys.  this requires all previous polys to be factored, since
			//we don't store the factors, just the polya coefficient, but trial
			//division is fast.


		}
	}

done:

	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);

	//record this a in the list
	sconf->poly_a_list = (z *)realloc(sconf->poly_a_list,
		(sconf->total_poly_a + 1) * sizeof(z));
	zInit(&sconf->poly_a_list[sconf->total_poly_a]);
	zCopy(poly_a,&sconf->poly_a_list[sconf->total_poly_a]);
	mp2gmp(poly_a, poly->mpz_poly_a);

	//sort the indices of factors of 'a'
	qsort(poly->qlisort,poly->s,sizeof(int),&qcomp_int);
	memset(&poly->qlisort[poly->s], 255, (MAX_A_FACTORS - poly->s) * sizeof(int));	

#ifdef POLYA_DEBUG
		printf("done generating poly_a\n");
#endif

	return;
}

void computeBl(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//ql = array of factors of a
	//Bl = array of generated Bl values
	//notation of polynomials, here and elsewhere, generally follows
	//contini's notation

	uint32 root1, root2, prime, gamma;
	uint32 amodql;	//(a/ql)^-1 mod ql = inv(a/ql mod ql) mod ql
	siqs_poly *poly = dconf->curr_poly;
	uint32 *modsqrt = sconf->modsqrt_array;
	fb_list *fb = sconf->factor_base;
	z *n = &sconf->n;
	int i, s = poly->s, *qli = poly->qlisort;
	z *Bl = dconf->Bl;

	//initialize b
	zCopy(&zZero,&poly->poly_b);

	for (i=0;i<s;i++)
	{
		prime = fb->list->prime[qli[i]];
		root1 = modsqrt[qli[i]];
		root2 = prime - root1; 
		
		zShortDiv(&poly->poly_a,(fp_digit)prime,&dconf->qstmp1);
		amodql = (uint32)zShortMod(&dconf->qstmp1,(fp_digit)prime);
		amodql = modinv_1(amodql,prime);

		//the primes will all be < 65536, so we can multiply safely
		gamma = (root1 * amodql) % prime;

		//check if the other root makes gamma smaller
		if (gamma > (prime>>1))
			gamma = prime-gamma;
		
		//qstmp1 holds a/prime
		zShortMul(&dconf->qstmp1,(fp_digit)gamma,&Bl[i]);

		//build up b
		zAdd(&poly->poly_b,&Bl[i],&dconf->qstmp1);
		//double Bl (the rest of the code wants it that way)
		zShiftLeft(&Bl[i],&Bl[i],1);
		zCopy(&dconf->qstmp1,&poly->poly_b);
	}

	//now that we have b, compute c = (b*b - n)/a
	zSqr(&poly->poly_b,&dconf->qstmp1);
	zSub(&dconf->qstmp1,n,&dconf->qstmp2);
	dconf->qstmp2.size *= -1;
	zDiv(&dconf->qstmp2,&poly->poly_a,&poly->poly_c,&dconf->qstmp1);
	poly->poly_c.size *= -1;
	mp2gmp(&poly->poly_b, poly->mpz_poly_b);
	mp2gmp(&poly->poly_c, poly->mpz_poly_c);

	/*
	// check the mpz representation
	// b^2 - c * a = n
	mpz_mul(dconf->gmptmp1, poly->mpz_poly_c, poly->mpz_poly_a);
	mpz_mul(dconf->gmptmp2, poly->mpz_poly_b, poly->mpz_poly_b);
	mpz_sub(dconf->gmptmp1, dconf->gmptmp2, dconf->gmptmp1);
	gmp2mp(dconf->gmptmp1, &dconf->qstmp1);
	if (zCompare(&dconf->qstmp1, &sconf->n) != 0)
	{
		printf("mpz poly not correct!\n");
		printf("%s, %s, %s\n", z2decstr(&poly->poly_a, &gstr1), z2decstr(&poly->poly_b, &gstr2), z2decstr(&poly->poly_c, &gstr3));
		gmp_printf("%Zd, %Zd, %Zd\n", poly->mpz_poly_a, poly->mpz_poly_b, poly->mpz_poly_c);
	}
	*/

	return;
}

void nextB(dynamic_conf_t *dconf, static_conf_t *sconf)
{
	//compute the ith b value for this polya
	//using a Gray code
	//b_i+1 = bi + 2*(-1)^ceil(i/2^v)*Bv
	//where 2^v is the highest power of 2 that divides 2*i
	//notation of polynomials, here and elsewhere, generally follows
	//contini's notation
	z *tmp, *tmp2;
	uint32 Bnum = dconf->numB;
	z *Bl = dconf->Bl;
	siqs_poly *poly = dconf->curr_poly;
	z *n = &sconf->n;

	tmp = &dconf->qstmp1;
	tmp2 = &dconf->qstmp2;

	//compute the next b
	if (poly->gray[Bnum] < 0)
	{
		zSub(&poly->poly_b,&Bl[poly->nu[Bnum]-1],tmp2);
		zCopy(tmp2,&poly->poly_b);
	}
	else
	{
		zAdd(&poly->poly_b,&Bl[poly->nu[Bnum]-1],tmp2);
		zCopy(tmp2,&poly->poly_b);
	}
	
	//now that we have b, compute c = (b*b - n)/a
	zSqr(&poly->poly_b,tmp);
	zSub(tmp,n,tmp2);
	tmp2->size *= -1;
	zDiv(tmp2,&poly->poly_a,&poly->poly_c,tmp);
	poly->poly_c.size *= -1;
	mp2gmp(&poly->poly_b, poly->mpz_poly_b);
	mp2gmp(&poly->poly_c, poly->mpz_poly_c);	
	
	// check the mpz representation
	// b^2 - c * a = n
	/*
	mpz_mul(dconf->gmptmp1, poly->mpz_poly_c, poly->mpz_poly_a);
	mpz_mul(dconf->gmptmp2, poly->mpz_poly_b, poly->mpz_poly_b);
	mpz_sub(dconf->gmptmp1, dconf->gmptmp2, dconf->gmptmp1);
	gmp2mp(dconf->gmptmp1, &dconf->qstmp1);
	if (zCompare(&dconf->qstmp1, &sconf->n) != 0)
	{
		printf("mpz poly not correct!\n");
		printf("%s, %s, %s\n", z2decstr(&poly->poly_a, &gstr1), z2decstr(&poly->poly_b, &gstr2), z2decstr(&poly->poly_c, &gstr3));
		gmp_printf("%Zd, %Zd, %Zd\n", poly->mpz_poly_a, poly->mpz_poly_b, poly->mpz_poly_c);
	}
	*/

	return;
}
