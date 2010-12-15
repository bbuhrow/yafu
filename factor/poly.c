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
//
//typedef struct 
//{
//	//read/write data inputs
//	uint32 *numptr_n;
//	uint32 *numptr_p;
//	bucket_element *sliceptr_n;
//	bucket_element *sliceptr_p;
//	update_t *update_data;
//	lp_bucket *lp_bucket_p;
//	int *ptr;
//
//	//read only inputs:
//	uint64 large_B;
//	uint64 B;
//	uint64 interval;
//	int64 numblocks;
//
//	//read/write words
//	uint64 bound_val;
//	int64 bound_index;
//	int64 check_bound;
//	uint64 logp;
//
//} polysieve_t;

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
	uint32 i,randindex, mindiff,a1,poly_low_found=0,target_bits;
	uint32 potential_a_factor, found_a_factor;
	uint32 afact[20];
	double target_mul = 0.9;
	FILE *sieve_log = sconf->obj->logfile;

	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);

	//determine polypool indexes.  
	//this really should be done once after generating the factor base
	UPPER_POLYPOOL_INDEX = fb->small_B - 1;

	for (i=0;i<fb->small_B;i++)
	{
		if ((fb->list->prime[i] > 1000) && (poly_low_found == 0))
		{
			LOWER_POLYPOOL_INDEX = i;
			poly_low_found=1;
		}

		if (fb->list->prime[i] > 4000)
		{
			UPPER_POLYPOOL_INDEX = i-1;
			break;
		}
	}
	UPPER_POLYPOOL_INDEX = fb->small_B - 1;

	//brute force the poly to be somewhat close to the target
	target_bits = (uint32)((double)zBits(target_a) * target_mul);

	while (1)
	{
		//generate poly_a's until the residue is 'small enough'
		//printf("*******new trial a********\n");

		//sp2z(1,poly_a);
		zCopy(&zOne,poly_a);
		*s=0;
		for (;;)
		{
			//randomly pick a new unique factor
			found_a_factor = 0;
			while (!found_a_factor)
			{
				randindex = (uint32)spRand((fp_digit)LOWER_POLYPOOL_INDEX,
					(fp_digit)UPPER_POLYPOOL_INDEX);
				//randindex = LOWER_POLYPOOL_INDEX + 
				//	(uint32)((UPPER_POLYPOOL_INDEX-LOWER_POLYPOOL_INDEX) * (double)rand() / (double)RAND_MAX);
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
			//printf("afactor %d = %u\n",*s,potential_a_factor);
			afact[*s]=potential_a_factor;
			qli[*s] = randindex;
			*s = *s + 1;
			//compute how close we are to target_a
			j = zBits(target_a) - zBits(poly_a);
			if (j < 10)
			{
				//too close, we want the last factor to be between 15 and 10 bits
				zCopy(&zOne,poly_a);
				*s=0;
				continue;
			}
			else if (j < 15)
			{
				//close enough to pick a last factor
				break;
			}
		}

		//at this point, poly_a is too small by one factor, find the closest factor
		zCopy(target_a,&tmp);
		zDiv(&tmp,poly_a,&tmp2,&tmp3);

		mindiff = 0xffffffff;
		a1 = tmp2.val[0];
		if (a1 < 1000)
			continue;

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
			//printf("last prime in poly_a > small_B\n");
			continue;
		}

		zShortMul(poly_a,fb->list->prime[randindex],poly_a);
		//printf("afactor %d = %u\n",*s,fb->list[randindex].prime);
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
						fb->list->prime[LOWER_POLYPOOL_INDEX],fb->list->prime[UPPER_POLYPOOL_INDEX],
						UPPER_POLYPOOL_INDEX,LOWER_POLYPOOL_INDEX);
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

	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);

	//record this a in the list
	sconf->poly_a_list = (z *)realloc(sconf->poly_a_list,
		(sconf->total_poly_a + 1) * sizeof(z));
	zInit(&sconf->poly_a_list[sconf->total_poly_a]);
	zCopy(poly_a,&sconf->poly_a_list[sconf->total_poly_a]);

	//sort the indices of factors of 'a'
	qsort(poly->qlisort,poly->s,sizeof(int),&qcomp_int);
	memset(&poly->qlisort[poly->s], 255, (MAX_A_FACTORS - poly->s) * sizeof(int));


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

	return;
}

void testfirstRoots(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//the roots are computed using a and b as follows:
	//(+/-t - b)(a)^-1 mod p
	//where the t values are the roots to t^2 = N mod p, found by shanks_tonelli
	//when constructing the factor base.
	//assume b > t

	//compute the roots as if we were actually going to use this, but don't save
	//anything.  We are just trying to determine the size needed for each large 
	//prime bucket by sieving over just the first bucket

	uint32 i,logp;
	int root1, root2, prime, amodp, bmodp, inv, bnum,numblocks;
	int lpnum,last_bound;

	//unpack stuff from the job data
	siqs_poly *poly = dconf->curr_poly;
	fb_list *fb = sconf->factor_base;
	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 *modsqrt = sconf->modsqrt_array;

	numblocks = sconf->num_blocks;

	lpnum = 0;
	dconf->buckets->alloc_slices = 1;

	//extreme estimate for number of slices
	i = (sconf->factor_base->B - sconf->factor_base->med_B) / 512;

	last_bound = fb->med_B;
	for (i=fb->med_B;i<fb->B;i++)
	{
		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 
		logp = fb->list->logprime[i];

		amodp = (int)zShortMod(&poly->poly_a,prime);
		bmodp = (int)zShortMod(&poly->poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		root1 = (int)root1 - bmodp;
		if (root1 < 0) root1 += prime;

		root2 = (int)root2 - bmodp;
		if (root2 < 0) root2 += prime;
	
		root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

		//just need to do this once, because the next step of prime will be 
		//into a different bucket
		bnum = root1 >> BLOCKBITS;
		if (bnum == 0)
			lpnum++;

		//repeat for the other root
		bnum = root2 >> BLOCKBITS;
		if (bnum == 0)
			lpnum++;

		if ((uint32)lpnum > (double)BUCKET_ALLOC * 0.75)
		{
			//we want to allocate more slices than we will probably need
			//assume alloc/2 is a safe amount of slack
			lp_bucket_p->alloc_slices++;
			lpnum = 0;
		}

		if (i - last_bound == 65536)
		{
			//when prime are really big, we may cross this boundary
			//before the buckets fill up
			lp_bucket_p->alloc_slices++;
			lpnum = 0;
			last_bound = i;
		}
	}

	return;
}

#define FILL_ONE_PRIME_P(i)					\
	if (root1 < interval)					\
	{										\
		bnum = root1 >> BLOCKBITS;			\
		bptr = sliceptr_p +					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum];					\
		/*bptr->fb_index = i - bound_val;	*/	\
		/*bptr->loc = root1 & BLOCKSIZEm1;	*/ \
		*bptr = ((i - bound_val) << 16) | (root1 & BLOCKSIZEm1); \
		numptr_p[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = sliceptr_p +					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum];					\
		/* bptr->fb_index = i - bound_val;*/		\
		/* bptr->loc = root2 & BLOCKSIZEm1;*/	\
		*bptr = ((i - bound_val) << 16) | (root2 & BLOCKSIZEm1); \
		numptr_p[bnum]++;					\
	}

#define FILL_ONE_PRIME_P2(i)					\
	if (root1 < interval)					\
	{										\
		bnum = root1 >> BLOCKBITS;			\
		bptr = next_loc_p[bnum];				\
		bptr->fb_index = i - bound_val;		\
		bptr->loc = root1 & BLOCKSIZEm1;	\
		next_loc_p[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = next_loc_p[bnum];				\
		bptr->fb_index = i - bound_val;		\
		bptr->loc = root2 & BLOCKSIZEm1;	\
		next_loc_p[bnum]++;					\
	}

#define FILL_ONE_PRIME_N(i)						\
	if (root1 < interval)						\
	{											\
		bnum = root1 >> BLOCKBITS;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		/*bptr->fb_index = i - bound_val;*/		\
		/*bptr->loc = root1 & BLOCKSIZEm1;*/	\
		*bptr = ((i - bound_val) << 16) | (root1 & BLOCKSIZEm1); \
		numptr_n[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		/*bptr->fb_index = i - bound_val;*/		\
		/*bptr->loc = root2 & BLOCKSIZEm1;*/	\
		*bptr = ((i - bound_val) << 16) | (root2 & BLOCKSIZEm1); \
		numptr_n[bnum]++;					\
	}			

#define FILL_ONE_PRIME_N3(i)						\
	if (root1 < interval)						\
	{											\
		bnum = root1 >> BLOCKBITS;			\
		bptr = sliceptr_p +	MAX_NUM_BLOCKS_PROD_BUCKET_BITS + \
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum+MAX_NUM_BLOCKS];					\
		bptr->fb_index = i - bound_val;		\
		bptr->loc = root1 & BLOCKSIZEm1;	\
		numptr_p[bnum+MAX_NUM_BLOCKS]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = sliceptr_p +	MAX_NUM_BLOCKS_PROD_BUCKET_BITS + 					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum+MAX_NUM_BLOCKS];					\
		bptr->fb_index = i - bound_val;		\
		bptr->loc = root2 & BLOCKSIZEm1;	\
		numptr_p[bnum+MAX_NUM_BLOCKS]++;					\
	}			

#define FILL_ONE_PRIME_N2(i)						\
	if (root1 < interval)						\
	{											\
		bnum = root1 >> BLOCKBITS;			\
		bptr = next_loc_n[bnum];				\
		bptr->fb_index = i - bound_val;		\
		bptr->loc = root1 & BLOCKSIZEm1;	\
		next_loc_n[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = next_loc_n[bnum];				\
		bptr->fb_index = i - bound_val;		\
		bptr->loc = root2 & BLOCKSIZEm1;	\
		next_loc_n[bnum]++;					\
	}			

#define FILL_ONE_PRIME_LOOP_P(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		/*bptr->fb_index = i - bound_val;*/			\
		/*bptr->loc = root1 & BLOCKSIZEm1;*/		\
		*bptr = ((i - bound_val) << 16) | (root1 & BLOCKSIZEm1); \
		numptr_p[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> BLOCKBITS;				\
	}											\
	bnum = root2 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		/*bptr->fb_index = i - bound_val;*/			\
		/*bptr->loc = root2 & BLOCKSIZEm1;*/		\
		*bptr = ((i - bound_val) << 16) | (root2 & BLOCKSIZEm1); \
		numptr_p[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> BLOCKBITS;				\
	} 

#define FILL_ONE_PRIME_LOOP_P2(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = next_loc_p[bnum];				\
		bptr->fb_index = i - bound_val;			\
		bptr->loc = root1 & BLOCKSIZEm1;		\
		next_loc_p[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> BLOCKBITS;				\
	}											\
	bnum = root2 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = next_loc_p[bnum];				\
		bptr->fb_index = i - bound_val;			\
		bptr->loc = root2 & BLOCKSIZEm1;		\
		next_loc_p[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> BLOCKBITS;				\
	} 

#define FILL_ONE_PRIME_LOOP_N(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_n +						\
			(bnum << BUCKET_BITS) +				\
			numptr_n[bnum];					\
		/*bptr->fb_index = i - bound_val;*/			\
		/*bptr->loc = root1 & BLOCKSIZEm1;*/		\
		*bptr = ((i - bound_val) << 16) | (root1 & BLOCKSIZEm1); \
		numptr_n[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> BLOCKBITS;				\
	}											\
	bnum = root2 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_n +						\
			(bnum << BUCKET_BITS) +				\
			numptr_n[bnum];					\
		/*bptr->fb_index = i - bound_val;*/			\
		/*bptr->loc = root2 & BLOCKSIZEm1;*/		\
		*bptr = ((i - bound_val) << 16) | (root2 & BLOCKSIZEm1); \
		numptr_n[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> BLOCKBITS;				\
	} 

#define FILL_ONE_PRIME_LOOP_N3(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +	MAX_NUM_BLOCKS_PROD_BUCKET_BITS +			\
			numptr_p[bnum+MAX_NUM_BLOCKS];					\
		bptr->fb_index = i - bound_val;			\
		bptr->loc = root1 & BLOCKSIZEm1;		\
		numptr_p[bnum+MAX_NUM_BLOCKS]++;					\
		root1 += prime;							\
		bnum = root1 >> BLOCKBITS;				\
	}											\
	bnum = root2 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +	MAX_NUM_BLOCKS_PROD_BUCKET_BITS +					\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum+MAX_NUM_BLOCKS];					\
		bptr->fb_index = i - bound_val;			\
		bptr->loc = root2 & BLOCKSIZEm1;		\
		numptr_p[bnum+MAX_NUM_BLOCKS]++;					\
		root2 += prime;							\
		bnum = root2 >> BLOCKBITS;				\
	} 

#define FILL_ONE_PRIME_LOOP_N2(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = next_loc_n[bnum];				\
		bptr->fb_index = i - bound_val;			\
		bptr->loc = root1 & BLOCKSIZEm1;		\
		next_loc_n[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> BLOCKBITS;				\
	}											\
	bnum = root2 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = next_loc_n[bnum];				\
		bptr->fb_index = i - bound_val;			\
		bptr->loc = root2 & BLOCKSIZEm1;		\
		next_loc_n[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> BLOCKBITS;				\
	} 

#define CHECK_NEW_SLICE(j)									\
	if (j >= check_bound)							\
	{														\
		room = 0;											\
		for (k=0;k<numblocks;k++)							\
		{													\
			if (*(numptr_p + k) > room)						\
				room = *(numptr_p + k);						\
			if (*(numptr_n + k) > room)						\
				room = *(numptr_n + k);						\
		}													\
		room = BUCKET_ALLOC - room;							\
		if (room < 32)										\
		{													\
			lp_bucket_p->logp[bound_index] = logp;			\
			bound_index++;									\
			lp_bucket_p->fb_bounds[bound_index] = j;		\
			bound_val = j;									\
			sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
			sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
			numptr_p += (numblocks << 1);							\
			numptr_n += (numblocks << 1);							\
			check_bound += BUCKET_ALLOC >> 1;					\
		}													\
		else												\
			check_bound += room >> 1;						\
	}										\
	else if ((j - bound_val) >= 65536)		\
	{										\
		lp_bucket_p->logp[bound_index] = logp;			\
		bound_index++;									\
		lp_bucket_p->fb_bounds[bound_index] = j;		\
		bound_val = j;									\
		sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
		sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
		numptr_p += (numblocks << 1);							\
		numptr_n += (numblocks << 1);							\
		check_bound += BUCKET_ALLOC >> 1;					\
	}

#define CHECK_NEW_SLICE3(j)									\
	if (j >= check_bound)							\
	{														\
		room = 0;											\
		for (k=0;k<numblocks;k++)							\
		{													\
			if (*(numptr_p + k) > room)						\
				room = *(numptr_p + k);						\
			if (*(numptr_p + k + MAX_NUM_BLOCKS) > room)						\
				room = *(numptr_p + k + MAX_NUM_BLOCKS);						\
		}													\
		room = BUCKET_ALLOC - room;							\
		if (room < 32)										\
		{													\
			lp_bucket_p->logp[bound_index] = logp;			\
			bound_index++;									\
			lp_bucket_p->fb_bounds[bound_index] = j;		\
			bound_val = j;									\
			sliceptr_p += (MAX_NUM_BLOCKS << (BUCKET_BITS + 1));		\
			numptr_p += (MAX_NUM_BLOCKS << 1);							\
			check_bound += BUCKET_ALLOC >> 1;					\
		}													\
		else												\
			check_bound += room >> 1;						\
	}										\
	else if ((j - bound_val) >= 65536)		\
	{										\
		lp_bucket_p->logp[bound_index] = logp;			\
		bound_index++;									\
		lp_bucket_p->fb_bounds[bound_index] = j;		\
		bound_val = j;									\
		sliceptr_p += (MAX_NUM_BLOCKS << (BUCKET_BITS + 1));		\
		numptr_p += (MAX_NUM_BLOCKS << 1);							\
		check_bound += BUCKET_ALLOC >> 1;					\
	}

#if defined(_MSC_VER)

#define ADDRESS_ROOT_1 (firstroots1[j])
#define ADDRESS_ROOT_2 (firstroots2[j])

#else

#define ADDRESS_ROOT_1 (*(firstroots1 + j))
#define ADDRESS_ROOT_2 (*(firstroots2 + j))

#endif

#if defined(MSC_ASM32A)
	#define COMPUTE_NEXT_ROOTS_P	\
	do {	\
		uint32 update = *ptr;	\
		ASM_M {					\
			ASM_M xor ebx, ebx	\
			ASM_M xor ecx, ecx	\
			ASM_M mov eax, root1	\
			ASM_M mov edx, root2	\
			ASM_M sub eax, update	\
			ASM_M cmovc ebx, prime	\
			ASM_M sub edx, update	\
			ASM_M cmovc ecx, prime	\
			ASM_M add eax, ebx	\
			ASM_M add edx, ecx	\
			ASM_M mov root1, eax	\
			ASM_M mov root2, edx}	\
		} while (0);		

	#define COMPUTE_NEXT_ROOTS_N	\
	do {	\
		uint32 update = *ptr;	\
		ASM_M {					\
			ASM_M mov eax, root1	\
			ASM_M mov edx, root2	\
			ASM_M mov ebx, eax		\
			ASM_M add ebx, update	\
			ASM_M mov ecx, edx		\
			ASM_M add ecx, update	\
			ASM_M sub eax, prime	\
			ASM_M sub edx, prime	\
			ASM_M add eax, update	\
			ASM_M cmovae eax, ebx	\
			ASM_M add edx, update	\
			ASM_M cmovae edx, ecx	\
			ASM_M mov root1, eax	\
			ASM_M mov root2, edx}	\
		} while (0);
	#define COMPUTE_FIRST_ROOTS	\
	do {	\
		ASM_M {					\
			ASM_M xor ebx, ebx	\
			ASM_M xor ecx, ecx	\
			ASM_M mov eax, root1	\
			ASM_M mov edx, root2	\
			ASM_M sub eax, bmodp	\
			ASM_M cmovc ebx, prime	\
			ASM_M sub edx, bmodp	\
			ASM_M cmovc ecx, prime	\
			ASM_M add eax, ebx	\
			ASM_M add edx, ecx	\
			ASM_M mov root1, eax	\
			ASM_M mov root2, edx}	\
		} while (0);		



#elif defined(GCC_ASM64X)

	#define COMPUTE_NEXT_ROOTS_P						\
		ASM_G (											\
			"xorl %%r8d, %%r8d		\n\t"	/*r8d = 0*/	\
			"xorl %%r9d, %%r9d		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - ptr*/	\
			"cmovc %3, %%r8d		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - ptr*/	\
			"cmovc %3, %%r9d		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%r8d, %%eax		\n\t"		\
			"addl %%r9d, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "g"(*ptr), "g"(prime)		\
			: "r8", "r9", "cc");	

	#define COMPUTE_NEXT_ROOTS_N		\
		ASM_G (							\
			"movl %%eax, %%r8d		\n\t"	\
			"addl %2, %%r8d			\n\t"	/*r8d = root1 + ptr*/		\
			"movl %%edx, %%r9d		\n\t"								\
			"addl %2, %%r9d			\n\t"	/*r9d = root2 + ptr*/		\
			"subl %3, %%eax			\n\t"	/*root1 = root1 - prime*/	\
			"subl %3, %%edx			\n\t"	/*root2 = root2 - prime*/	\
			"addl %2, %%eax			\n\t"	/*root1 + ptr*/				\
			"cmovae %%r8d, %%eax	\n\t"	/*other caluclation if no overflow*/	\
			"addl %2, %%edx			\n\t"	/*root2 + ptr*/							\
			"cmovae %%r9d, %%edx	\n\t"	/*other caluclation if no overflow*/	\
			: "+a"(root1), "+d"(root2)		\
			: "g"(*ptr), "g"(prime)			\
			: "r8", "r9", "cc");

	#define COMPUTE_FIRST_ROOTS			\
		ASM_G (											\
			"xorl %%r8d, %%r8d		\n\t"	/*r8d = 0*/	\
			"xorl %%r9d, %%r9d		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - bmodp*/	\
			"cmovc %3, %%r8d		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - bmodp*/	\
			"cmovc %3, %%r9d		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%r8d, %%eax		\n\t"		\
			"addl %%r9d, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "r"(bmodp), "g"(prime)		\
			: "r8", "r9", "cc");	

	#define COMPUTE_NEXT_ROOTS_X2_P						\
		ASM_G (											\
			"xorl %%r8d, %%r8d		\n\t"	/*r8d = 0*/	\
			"xorl %%r9d, %%r9d		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - ptr*/	\
			"cmovc %3, %%r8d		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - ptr*/	\
			"cmovc %3, %%r9d		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%r8d, %%eax		\n\t"		\
			"addl %%r9d, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "g"(*ptr2), "g"(prime)		\
			: "r8", "r9", "cc");	

	#define COMPUTE_NEXT_ROOTS_X2_N		\
		ASM_G (							\
			"movl %%eax, %%r8d		\n\t"	\
			"addl %2, %%r8d			\n\t"	/*r8d = root1 + ptr*/		\
			"movl %%edx, %%r9d		\n\t"								\
			"addl %2, %%r9d			\n\t"	/*r9d = root2 + ptr*/		\
			"subl %3, %%eax			\n\t"	/*root1 = root1 - prime*/	\
			"subl %3, %%edx			\n\t"	/*root2 = root2 - prime*/	\
			"addl %2, %%eax			\n\t"	/*root1 + ptr*/				\
			"cmovae %%r8d, %%eax	\n\t"	/*other caluclation if no overflow*/	\
			"addl %2, %%edx			\n\t"	/*root2 + ptr*/							\
			"cmovae %%r9d, %%edx	\n\t"	/*other caluclation if no overflow*/	\
			: "+a"(root1), "+d"(root2)		\
			: "g"(*ptr2), "g"(prime)			\
			: "r8", "r9", "cc");


#else
	#define COMPUTE_NEXT_ROOTS_P		\
		root1 = (int)root1 - *ptr;		\
		root2 = (int)root2 - *ptr;		\
		root1 += ((root1 >> 31) * prime);			\
		root2 += ((root2 >> 31) * prime);	

	#define COMPUTE_NEXT_ROOTS_N		\
		root1 = (int)root1 + *ptr;		\
		root2 = (int)root2 + *ptr;		\
		root1 -= ((root1 >= prime) * prime);	\
		root2 -= ((root2 >= prime) * prime);	

	#define COMPUTE_FIRST_ROOTS		\
		root1 = (int)root1 - bmodp;		\
		root2 = (int)root2 - bmodp;		\
		if (root1 < 0) root1 += prime;			\
		if (root2 < 0) root2 += prime;

#endif

/*
firstRoots and nextRoots do two things: prepare all primes
for sieving for a new polynomial according to the siqs
algorithm, and pre-sieve large primes.  The preparation of
primes amounts to finding the roots for the current (a,b)
polynomial pair and generally follows contini's outline using
gray codes.  These first roots are stored in a data structure
along with the prime and logp value of size 16 bytes.  
Along with these, the vector rootupdates is computed which 
holds an increment to be applied to each prime's roots for
each subsesquest b poly.  Then for each new polynomial, the value
stored in firstroots1 or firstroots2 is incremented by
the value in rootupdates and re-stored.  

When sorting roots for each new polynomial the initial root
info is first loaded, and since all this info is adjacent to each
other in a small structure, it is all cached into the same
cache line.  this reduces memory access a bit beyond previous
versions.  side note has idea for further reducing memory
access and cache pollution.

== side note
perhaps some asm code and/or
loop unrolling will allow more than one structure to be loaded
at once from the 64 byte cache line and saved in registers
before the cache line is updated when the roots are pushed
out to buckets.  this will reduce memory access by the need 
to reload the root update cache line for every prime, and
instead to go that same data saved in registers.  probably only
useful for 64 bit systems with extended register set... unless
we can store some in multimedia registers...
== end side note


after getting the starting root for the poly as outlined, 
we proceed one of three ways depending on the size of the
prime.  For those primes less than a predetermined 'large'
size, we load up the roots into the sieving factor base
structure.  This structure is used during the sieve routine
to do the actual sieving, as well as in the trial division
routine to know where to look for divisible entries of the 
sieve array for a particular prime.  For primes larger than
this bound, but less than the size of the entire sieve 
interval, we proceed with large prime sieving.  Large prime
sieving is based on the notion of bucket sorting.  This
means we create buckets, one for each block of the sieve
interval, and sort the roots into those buckets.  A fast
right shift by the (power of 2 sized) blocksize will tell
us what block a given prime's root will fall in.  We dump 
the location in that block (given by a fast AND operation)
into the bucket, along with the index of that prime in the
factor base.  We proceed until the primes become larger than
the size of the entire sieve interval.  This is true for a
great many primes, because since this method of changing 
polynomials is so fast we can keep the sieve interval very
small.  Once this boundary is crossed, the idea is the same
but instead of looping through the buckets we just need
to compute the one single bucket into which that prime will
fall for this poly.  once we are done sorting all the primes
into buckets, sieving of a block amounts to just dumping 
the contents of the appropriate bucket into the block.  This
is very fast because essentially every access of the sieve
block is a L1 cache hit as we linearly walk through the 
contents of the bucket.  The final note is related to 
the concept of factor base slices.  This is done in order 
to keep the size of the buckets small, and close to L1.  
The size of each slice depends on the 
number of buckets in the interval as well as the frequency
with which a prime hits a bucket.  This is best explained
by example.  Say that we have determined that the sieve 
interval needs to be 10 blocks of size 32k.  When sorting, 
ideally all buckets fit in the 32k L1, thus each
bucket gets 3.2kB.  Each element of a bucket is 4 bytes
(consisting of a 2 byte location and a 2 byte factor base 
offset).  Thus there is room for about 819 elements in each
bucket.  But additions to a bucket are essentially random,
and the frequency of additions to a given bucket depends on
the size of the primes being considered.  Primes only slightly
larger than the blocksize, for instance, might hit roughly
every bucket in the interval while a much larger prime will
only hit one of them.  So as the primes get bigger the buckets
fill more slowly.  We create a new slice of the factor base
when the buckets are nearly full.  The new slice is just 
another set of 10 buckets immediately adjacent to the previous
10 in memory.  We remember the boundaries of the slices as
we go, and all factor base offsets are given relative to the
slice boundary, rather than to the start of the factor base.
that way, we can get away with only using 16 bits to hold
the offset.  16 bits to hold the block location limits the
size of a blocksize to 2^16 or less.  and 16 blocks to hold
the FB offset limits the size of the FB slice, although this
only applies for the largest factorizations.  for each slice
we also need to record the logp, which for small slices is
generally constant for all primes in the slice.  thus
we don't need to add logp to the bucket element for each
prime, further reducing the cache footprint of a bucket
element.  In practice, it seems to be faster to not restrict
the buckets to fit entirely in L1, but to be a little larger.
Perhaps this is because if the buckets are smaller, then we
need more slices, and looping over the slices during the sieve
has overhead which negates the benefits.  probably it is most
important to not exceed a page size.  So instead, the bucket 
size is fixed using a compile time constant, but the 
filling procedure stays the same.  We also kept lists of 
buckets for pos and neg sieve blocks adjacent to each other
in linear memory, rather than allocating a set of buckets for
pos side and another for neg, because the sorting loop bounces
back and forth between dumping in positive roots and negative
ones, and its beneficial to keep the two sets of buckets close
to each other within a slice.  This is the origin of the fixed
offset of (numblocks << BUCKET_BITS) between pos and neg
bucket pointers.

This awful mess makes for one fast sieve.
*/

void firstRoots(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//the roots are computed using a and b as follows:
	//(+/-t - b)(a)^-1 mod p
	//where the t values are the roots to t^2 = N mod p, found by shanks_tonelli
	//when constructing the factor base.
	//assume b > t

	//unpack stuff from the job data structures
	siqs_poly *poly = dconf->curr_poly;
	fb_list *fb = sconf->factor_base;
	z *Bl = dconf->Bl;
	uint32 start_prime = 2;
	int *rootupdates = dconf->rootupdates;
#ifdef UPDATEDATA_AOS
	update_t *update_data = dconf->update_data;
#else
	update_t update_data = dconf->update_data;
#endif
	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 *modsqrt = sconf->modsqrt_array;

	//locals
	uint32 i, interval;
	uint8 logp;
	int root1, root2, prime, amodp, bmodp, inv, x, bnum,j,numblocks;
	int s = poly->s;
	int bound_index = 0, k;
	uint32 bound_val = fb->med_B;
#ifdef USEBUCKETSTRUCT
	bucket_element *bptr,*sliceptr_p,*sliceptr_n;
#else
	uint32 *bptr, *sliceptr_p, *sliceptr_n;
#endif
	uint32 *numptr_p, *numptr_n;
	int check_bound = BUCKET_ALLOC/2 - 1, room;

	numblocks = sconf->num_blocks;
	interval = numblocks << BLOCKBITS;

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->fb_bounds[0] = fb->med_B;

#ifdef USEMAXBLOCKS
		sliceptr_p = lp_bucket_p->list;
		//sliceptr_n = lp_bucket_p->list + MAX_NUM_BLOCKS_PROD_BUCKET_BITS;

		numptr_p = lp_bucket_p->num;
		//numptr_n = lp_bucket_p->num + MAX_NUM_BLOCKS;
		//reset lp_buckets
		for (i=0;i< (2*MAX_NUM_BLOCKS*lp_bucket_p->alloc_slices) ;i++)
			numptr_p[i] = 0;
#else
		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		//reset lp_buckets
		for (i=0;i< (2*numblocks*lp_bucket_p->alloc_slices) ;i++)
			numptr_p[i] = 0;
#endif

		lp_bucket_p->num_slices = 0;
	}
	else
	{
		sliceptr_p = NULL;
		sliceptr_n = NULL;
		numptr_p = NULL;
		numptr_n = NULL;
	}

	for (i=start_prime;i<fb->med_B;i++)
	{
		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)zShortMod(&poly->poly_a,prime);
		bmodp = (int)zShortMod(&poly->poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

		if (root2 < root1)
		{
#ifdef UPDATEDATA_AOS
			update_data[i].firstroots1 = root2;
			update_data[i].firstroots2 = root1;
#else
			update_data.firstroots1[i] = root2;
			update_data.firstroots2[i] = root1;
#endif
			//update_data[i].firstroots1 = (int)((root1 << 16) | root2);
			fb_p[i].roots = (uint32)((root1 << 16) | root2);
			fb_n[i].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
		}
		else
		{
#ifdef UPDATEDATA_AOS
			update_data[i].firstroots1 = root1;
			update_data[i].firstroots2 = root2;
#else
			update_data.firstroots1[i] = root1;
			update_data.firstroots2[i] = root2;
#endif
			//update_data[i].firstroots1 = (int)((root2 << 16) | root1);
			fb_p[i].roots = (uint32)((root2 << 16) | root1);
			fb_n[i].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
		}

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)zShortMod(&Bl[j],prime);
			x = (int)((int64)x * (int64)inv % (int64)prime);
			rootupdates[(j)*fb->B+i] = x;
		}
	}

	check_bound = fb->med_B + BUCKET_ALLOC/2;
	for (i=fb->med_B;i<fb->large_B;i++)
	{
		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1; 
		logp = fb->list->logprime[i];

		amodp = (int)zShortMod(&poly->poly_a,prime);
		bmodp = (int)zShortMod(&poly->poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

#ifdef USEMANYSLICES
#ifdef USEMAXBLOCKS
		CHECK_NEW_SLICE3(i);
#else
		CHECK_NEW_SLICE(i);
#endif
#endif

#ifdef UPDATEDATA_AOS
		update_data[i].firstroots1 = root1;
		update_data[i].firstroots2 = root2;
#else
		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;
#endif

		FILL_ONE_PRIME_LOOP_P(i);

#ifdef UPDATEDATA_AOS
		root1 = (prime - update_data[i].firstroots1);
		root2 = (prime - update_data[i].firstroots2);
#else
		root1 = (prime - update_data.firstroots1[i]);
		root2 = (prime - update_data.firstroots2[i]);
#endif
	
#ifdef USEMAXBLOCKS
		FILL_ONE_PRIME_LOOP_N3(i);
#else
		FILL_ONE_PRIME_LOOP_N(i);
#endif

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)zShortMod(&Bl[j],prime);
			x = (int)((int64)x * (int64)inv % (int64)prime);
			rootupdates[(j)*fb->B+i] = x;
		}
		//we could do the following:
		//loop over very large primes
		//	compute the s root updates for this prime
		//	loop over (a few of) the 2^(s-1) polys
		//		compute the roots for this poly
		//		check the 4 candidates against interval size
		//		only put those that survive in a matrix where:
		//			rows are bucket elements for 1 poly
		//			cols are different polys
		//	end
		//end
		//
		//we then sort a row using a cache friendly quicksort and the 
		//elements are ready to be dumped into a sieve.  this eliminates
		//the vlp bucket sorting in nextroots.
		//matrix rows are contiguous in memory
		//matrix cols are not - but we always write to the head of each column
		//when looping over polys, so they should still be in L1.
		//we benefit by not needing to store the rootupdates array and
		//by not needing to update the firstroots fields in update data
		//as often (only at the end of the loop over polys).  I suppose
		//we could either do a linear dump of bucket elements followed by a 
		//sort or we could directly bucket sort... hmmm.


	}

	for (i=fb->large_B;i<fb->B;i++)
	{
		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1; 
		logp = fb->list->logprime[i];

		amodp = (int)zShortMod(&poly->poly_a,prime);
		bmodp = (int)zShortMod(&poly->poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

#ifdef USEMANYSLICES
#ifdef USEMAXBLOCKS
		CHECK_NEW_SLICE3(i);
#else
		CHECK_NEW_SLICE(i);
#endif
#endif

#ifdef UPDATEDATA_AOS
		update_data[i].firstroots1 = root1;
		update_data[i].firstroots2 = root2;
#else
		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;
#endif

		FILL_ONE_PRIME_P(i);

		root1 = (prime - root1);
		root2 = (prime - root2);

#ifdef USEMAXBLOCKS
		FILL_ONE_PRIME_N3(i);
#else
		FILL_ONE_PRIME_N(i);
#endif

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		//s is the number of primes in 'a'
		for (j=0;j<s;j++)
		{
			x = (int)zShortMod(&Bl[j],prime);
			x = (int)((int64)x * (int64)inv % (int64)prime);
			rootupdates[(j)*fb->B+i] = x;
		}
	}

	if (lp_bucket_p->list != NULL)
		lp_bucket_p->num_slices = bound_index + 1;
	

	return;
}

//this is in the poly library, even though the bulk of the time is spent
//bucketizing large primes, because it's where the roots of a poly are updated
void nextRoots(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//update the roots 

	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	int *rootupdates = dconf->rootupdates;

#ifdef UPDATEDATA_AOS
	update_t *update_data = dconf->update_data;
#else
	update_t update_data = dconf->update_data;
#endif

	uint32 startprime = 2;
	uint32 bound = sconf->factor_base->B;

	char v = dconf->curr_poly->nu[dconf->numB];
	char sign = dconf->curr_poly->gray[dconf->numB];
	int *ptr;

	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 med_B = sconf->factor_base->med_B;
	uint32 large_B = sconf->factor_base->large_B;

	uint32 j, interval;
	int k,bnum,numblocks,room;
	uint32 root1, root2, prime;
	//uint32 *scratch;
	//uint32 ptroffset;
	
	int bound_index=0;
	int check_bound = BUCKET_ALLOC/2 - 1;
	uint32 bound_val = med_B;
#ifdef USEBUCKETSTRUCT
	bucket_element *bptr,*sliceptr_p,*sliceptr_n;
#else
	uint32 *bptr, *sliceptr_p,*sliceptr_n;
#endif
	uint32 *numptr_p, *numptr_n;
	uint8 logp=0;
	//polysieve_t helperstruct;

	//scratch = (uint32 *)memalign(64,4 * sizeof(uint32));
	//scratch[0] = 1;
	//scratch[1] = 1;
	//scratch[2] = 1;
	//scratch[3] = 1;

	numblocks = sconf->num_blocks;
	interval = numblocks << BLOCKBITS;
	
	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->fb_bounds[0] = med_B;

#ifdef USEMAXBLOCKS
		sliceptr_p = lp_bucket_p->list;

		numptr_p = lp_bucket_p->num;

		//reuse this for a sec...
		prime = 2*MAX_NUM_BLOCKS*lp_bucket_p->alloc_slices;

		//reset lp_buckets
		for (j=0;j<prime;j++)
			numptr_p[j] = 0;
#else
		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		
		//reuse this for a sec...
		prime = 2*numblocks*lp_bucket_p->alloc_slices;

		//reset lp_buckets
		for (j=0;j<prime;j++)
			numptr_p[j] = 0;
#endif
	
		lp_bucket_p->num_slices = 0;

	}
	else
	{
		sliceptr_p = NULL;
		sliceptr_n = NULL;
		numptr_p = NULL;
		numptr_n = NULL;
	}

	k=0;
	ptr = &rootupdates[(v-1) * bound + startprime];
	//ptroffset = (v-1) * bound + startprime;

	if (sign > 0)
	{

		// update all roots, then sort into buckets or the compressed sieve fb.
		// first get to a spot where we can use movdqa on rootupdates and the update data fields
		// this means that the index into all of these arrays is even (a multiple of 64 bytes)
		// (v-1)   bound   startprime   ptroffset   initial j   possible?
		//  odd     odd       odd         even      odd            n
		//  odd     odd       even        odd       even           n
		//  odd     even      odd         odd       odd            y
		//  odd     even      even        even      even           y
		//  even    odd       odd         odd       odd            y
		//  even    odd       even        even      even           y
		//  even    even      odd         odd       odd            y
		//  even    even      even        even      even           y
		// the above truth table suggests that bound should be forced to be even
		//
		// argh.  the condition moves, etc, just seem to make this unsuitable for
		// multimedia operations.  at least without pmulld.

		//if (startprime & 0x1)
		//{
		//	j = startprime;
		//	prime = update_data.prime[j];
		//	root1 = update_data.firstroots1[j];
		//	root2 = update_data.firstroots2[j];
		//	COMPUTE_NEXT_ROOTS_P;
		//	update_data.firstroots1[j] = root1;
		//	update_data.firstroots2[j] = root2;
		//	j++;
		//	ptr++;
		//}
		//else
		//	j = startprime;		


		//ASM_G (											\
		//	"xorl %%r8d, %%r8d		\n\t"	/*r8d = 0*/	\
		//	"xorl %%r9d, %%r9d		\n\t"	/*r9d = 0*/	\
		//	"subl %2, %%eax			\n\t"	/*root1 - ptr*/	\
		//	"cmovc %3, %%r8d		\n\t"	/*prime into r8 if overflow*/	\
		//	"subl %2, %%edx			\n\t"	/*root2 - ptr*/	\
		//	"cmovc %3, %%r9d		\n\t"	/*prime into r9 if overflow*/	\
		//	"addl %%r8d, %%eax		\n\t"		\
		//	"addl %%r9d, %%edx		\n\t"		\
		//	: "+a"(root1), "+d"(root2)			\
		//	: "g"(*ptr), "g"(prime)		\
		//	: "r8", "r9", "cc");	


		//asm volatile (
		//	"movdqa (%1), %%xmm1 \n\t" /* mov roots1 into sse2 reg */
		//	"movdqa (%2), %%xmm2 \n\t" /* mov roots2 into sse2 reg */
		//	"movdqa (%0), %%xmm3 \n\t" /* mov update values into sse2 reg */
		//	"movdqa (%3), %%xmm4 \n\t" /* mov primes into sse2 reg */
		//	"psubd %%xmm3, %%xmm1 \n\t" /* compute new root1s */
		//	"psubd %%xmm3, %%xmm2 \n\t" /* compute new root2s */
		//	"movdqa %%xmm1, %%xmm6 \n\t" /* copy new roots */
		//	"movdqa %%xmm2, %%xmm7 \n\t" /* copy new roots */
		//	"psrad 31, %%xmm1 \n\t" /* shift right all doublewords 31 bits */
		//	"psrad 31, %%xmm2 \n\t" /* shift right all doublewords 31 bits */
		//	/*"pand (%4), %%xmm1 \n\t" /* clear possible sign bits from shifting */
		//	/*"pand (%4), %%xmm2 \n\t" /* clear possible sign bits from shifting */
		//	/*"pmulld %%xmm4, %%xmm1 \n\t" /* multiply by primes */
		//	/*"pmulld %%xmm4, %%xmm2 \n\t" /* multiply by primes */
		//	"pmovmskb %%xmm1, %%eax \n\t" /* extract sign bits */
		//	"pmovmskb %%xmm2, %%ebx \n\t" /* extract sign bits */
		//	"paddd %%xmm6, %%xmm1 \n\t" /* add correction to new roots */
		//	"paddd %%xmm7, %%xmm2 \n\t" /* add correction to new roots */
		//	"movdqa %%xmm1, (%1) \n\t" /* store new root1s back to memory */
		//	"movdqa %%xmm2, (%2) \n\t" /* store new root2s back to memory */
		//	: 
		//	: "r"(ptr), "r"(update_data.firstroots1 + j), "r"(update_data.firstroots2 + j), "r"(update_data.prime + j), "r"(scratch), "g"(j), "g"(bound)
		//	: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");



#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=startprime;j<med_B;j++,ptr++)
		{
			
#ifdef UPDATEDATA_AOS
			prime = update_data[j].prime;
			root1 = update_data[j].firstroots1;
			root2 = update_data[j].firstroots2;
#else
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
#endif

			COMPUTE_NEXT_ROOTS_P;
			
			if (root2 < root1)
			{
#ifdef UPDATEDATA_AOS
				update_data[j].firstroots1 = root2;
				update_data[j].firstroots2 = root1;
#else
				update_data.firstroots1[j] = root2;
				update_data.firstroots2[j] = root1;
#endif
				fb_p[j].roots = (uint32)((root1 << 16) | root2);
				fb_n[j].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
			}
			else
			{
#ifdef UPDATEDATA_AOS
				update_data[j].firstroots1 = root1;
				update_data[j].firstroots2 = root2;
#else
				update_data.firstroots1[j] = root1;
				update_data.firstroots2[j] = root2;
#endif
				fb_p[j].roots = (uint32)((root2 << 16) | root1);
				fb_n[j].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
			}
		}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		
		for (j=med_B;j<large_B;j++,ptr++)
		{
#ifdef UPDATEDATA_AOS
			prime = update_data[j].prime;
			logp = update_data[j].logp;

			root1 = update_data[j].firstroots1;
			root2 = update_data[j].firstroots2;
#else
			prime = update_data.prime[j];
			logp = update_data.logp[j];

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
#endif

			COMPUTE_NEXT_ROOTS_P;

#ifdef UPDATEDATA_AOS
			update_data[j].firstroots1 = root1;
			update_data[j].firstroots2 = root2;
#else
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;
#endif

#ifdef USEMANYSLICES
#ifdef USEMAXBLOCKS
			CHECK_NEW_SLICE3(j);
#else
			CHECK_NEW_SLICE(j);
#endif
#endif

			FILL_ONE_PRIME_LOOP_P(j);

#ifdef UPDATEDATA_AOS
			root1 = (prime - update_data[j].firstroots1);
			root2 = (prime - update_data[j].firstroots2);
#else
			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);
#endif

#ifdef USEMAXBLOCKS
			FILL_ONE_PRIME_LOOP_N3(j);
#else
			FILL_ONE_PRIME_LOOP_N(j);
#endif
		}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif
			

		for (j=large_B;j<bound;j++,ptr++)				
		{			
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_P;

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);		
			root2 = (prime - root2);	

			FILL_ONE_PRIME_N(j);

		}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

	}
	else
	{

#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=startprime;j<med_B;j++,ptr++)
		{
			
#ifdef UPDATEDATA_AOS
			root1 = update_data[j].firstroots1;
			root2 = update_data[j].firstroots2;
			prime = update_data[j].prime;
#else
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];
#endif

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
#ifdef UPDATEDATA_AOS
				update_data[j].firstroots1 = root2;
				update_data[j].firstroots2 = root1;
#else
				update_data.firstroots1[j] = root2;
				update_data.firstroots2[j] = root1;
#endif
				//update_data[j].firstroots1 = (int)((root1 << 16) | root2);
				fb_p[j].roots = (uint32)((root1 << 16) | root2);
				fb_n[j].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
			}
			else
			{
#ifdef UPDATEDATA_AOS
				update_data[j].firstroots1 = root1;
				update_data[j].firstroots2 = root2;
#else
				update_data.firstroots1[j] = root1;
				update_data.firstroots2[j] = root2;
#endif
				//update_data[j].firstroots1 = (int)((root2 << 16) | root1);
				fb_p[j].roots = (uint32)((root2 << 16) | root1);
				fb_n[j].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
			}
		}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		
		for (j=med_B;j<large_B;j++,ptr++)
		{
#ifdef UPDATEDATA_AOS
			prime = update_data[j].prime;
			logp = update_data[j].logp;

			root1 = update_data[j].firstroots1;
			root2 = update_data[j].firstroots2;
#else
			prime = update_data.prime[j];
			logp = update_data.logp[j];

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
#endif

			COMPUTE_NEXT_ROOTS_N;

#ifdef UPDATEDATA_AOS
			update_data[j].firstroots1 = root1;
			update_data[j].firstroots2 = root2;
#else
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;
#endif

#ifdef USEMANYSLICES
#ifdef USEMAXBLOCKS
			CHECK_NEW_SLICE3(j);
#else
			CHECK_NEW_SLICE(j);
#endif
#endif

			FILL_ONE_PRIME_LOOP_P(j);

#ifdef UPDATEDATA_AOS
			root1 = (prime - update_data[j].firstroots1);
			root2 = (prime - update_data[j].firstroots2);
#else
			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);
#endif

#ifdef USEMAXBLOCKS
			FILL_ONE_PRIME_LOOP_N3(j);
#else
			FILL_ONE_PRIME_LOOP_N(j);
#endif

		}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=large_B;j<bound;j++,ptr++)				
		{				
#ifdef USEMANYSLICES
#ifdef USEMAXBLOCKS
			CHECK_NEW_SLICE3(j);
#else
			CHECK_NEW_SLICE(j);
#endif
#endif

#ifdef UPDATEDATA_AOS
			prime = update_data[j].prime;			
			logp = update_data[j].logp;				
			root1 = update_data[j].firstroots1;	
			root2 = update_data[j].firstroots2;	
#else
			prime = update_data.prime[j];			
			logp = update_data.logp[j];				
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	
#endif

			COMPUTE_NEXT_ROOTS_N;				

#ifdef UPDATEDATA_AOS
			update_data[j].firstroots1 = root1;	
			update_data[j].firstroots2 = root2;	
#else
			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	
#endif

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
#ifdef USEMAXBLOCKS
			FILL_ONE_PRIME_N3(j);
#else
			FILL_ONE_PRIME_N(j);
#endif

		}

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

	}

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->num_slices = bound_index + 1;
		lp_bucket_p->logp[bound_index] = logp;
	}

	//free(scratch);
	return;
}
