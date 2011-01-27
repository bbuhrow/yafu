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

typedef struct 
{
	//read/write data inputs
	uint32 *numptr_n;
	uint32 *numptr_p;
	uint32 *sliceptr_n;
	uint32 *sliceptr_p;
	uint32 *update_data_prime;
	int *update_data_root1;
	int *update_data_root2;
	uint8 *update_data_logp;
	lp_bucket *lp_bucket_p;
	int *ptr;

	//read only inputs:
	uint32 large_B;
	uint32 B;
	uint32 interval;
	int numblocks;

	//read/write words
	uint32 bound_val;
	int bound_index;
	int check_bound;
	uint8 logp;

} polysieve_t;

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
		/* *bptr = ((i - bound_val) << 16) | (root1 & BLOCKSIZEm1); */ \
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
		numptr_p[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = sliceptr_p +					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_p[bnum]++;					\
	}

#define FILL_ONE_PRIME_N(i)						\
	if (root1 < interval)						\
	{											\
		bnum = root1 >> BLOCKBITS;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
		numptr_n[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_n[bnum]++;					\
	}			

#define FILL_ONE_PRIME_LOOP_P(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
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
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_p[bnum]++;					\
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
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
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
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_n[bnum]++;					\
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
			logp = update_data.logp[j];						\
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


#if defined(_MSC_VER)

	#define ADDRESS_ROOT_1 (firstroots1[j])
	#define ADDRESS_ROOT_2 (firstroots2[j])

	#define COMPUTE_4_PROOTS(j)								\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				ptrs = _mm_load_si128((__m128i *)(&rootupdates[(v-1) * bound + j])); \
				root1s = _mm_load_si128((__m128i *)(update_data.firstroots1 + j)); \
				root1s = _mm_sub_epi32(root1s, ptrs); 	 					/* root1 -= ptr */ \
				root2s = _mm_load_si128((__m128i *)(update_data.firstroots2 + j)); \
				root2s = _mm_sub_epi32(root2s, ptrs); 	 					/* root2 -= ptr */ \
				tmp1 = _mm_xor_si128(tmp1, tmp1); 							/* zero xmm4 */ \
				tmp2 = _mm_xor_si128(tmp2, tmp2);							/* zero xmm5 */ \
				primes = _mm_load_si128((__m128i *)(update_data.prime + j)); \
				tmp1 = _mm_cmpgt_epi32(tmp1,root1s); 						/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
				tmp2 = _mm_cmpgt_epi32(tmp2,root2s); 						/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
				tmp1 = _mm_and_si128(tmp1, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				tmp2 = _mm_and_si128(tmp2, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				root1s = _mm_add_epi32(root1s, tmp1); 						/* selectively add back prime (modular subtract) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots1 + j),root1s);		/* save new root1 values */ \
				root2s = _mm_add_epi32(root2s, tmp2); 						/* selectively add back prime (modular subtract) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots2 + j),root2s); 		/* save new root2 values */ \
			} while (0);

	#define COMPUTE_4_NROOTS(j)								\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				ptrs = _mm_load_si128((__m128i *)(&rootupdates[(v-1) * bound + j])); \
				root1s = _mm_load_si128((__m128i *)(update_data.firstroots1 + j)); \
				root1s = _mm_add_epi32(root1s, ptrs); 	 					/* root1 += ptr */ \
				root2s = _mm_load_si128((__m128i *)(update_data.firstroots2 + j)); \
				root2s = _mm_add_epi32(root2s, ptrs); 	 					/* root2 += ptr */ \
				tmp1 = _mm_shuffle_epi32(root1s, 0xe4); 					/* copy root1 to xmm4 */ \
				tmp2 = _mm_shuffle_epi32(root2s, 0xe4);						/* copy root2 to xmm5 */ \
				primes = _mm_load_si128((__m128i *)(update_data.prime + j)); \
				tmp1 = _mm_cmpgt_epi32(tmp1,primes); 						/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
				tmp2 = _mm_cmpgt_epi32(tmp2,primes); 						/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
				tmp1 = _mm_and_si128(tmp1, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				tmp2 = _mm_and_si128(tmp2, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				root1s = _mm_sub_epi32(root1s, tmp1); 						/* selectively sub back prime (modular addition) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots1 + j),root1s);		/* save new root1 values */ \
				root2s = _mm_sub_epi32(root2s, tmp2); 						/* selectively sub back prime (modular addition) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots2 + j),root2s); 		/* save new root2 values */ \
			} while (0);


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

	#define COMPUTE_4_PROOTS(j)								\
		ASM_G (											\
			"movdqa (%%rax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%rcx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%rdx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%rbx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%rcx) \n\t"			/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%rdx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	

	#define COMPUTE_4_NROOTS(j)								\
		ASM_G (											\
			"movdqa (%%rax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%rcx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%rdx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%rbx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%rcx) \n\t"			/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%rdx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	

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

#elif defined(GCC_ASM32X)

	#define COMPUTE_NEXT_ROOTS_P						\
		ASM_G (											\
			"xorl %%ecx, %%ecx		\n\t"	/*r8d = 0*/	\
			"xorl %%edi, %%edi		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - ptr*/	\
			"cmovc %3, %%ecx		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - ptr*/	\
			"cmovc %3, %%edi		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%ecx, %%eax		\n\t"		\
			"addl %%edi, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "g"(*ptr), "g"(prime)		\
			: "ecx", "edi", "cc");	

	#define COMPUTE_4_PROOTS(j)								\
		ASM_G (											\
			"movdqa (%%eax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%ecx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%edx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%ebx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%ecx) \n\t"			/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%edx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	

	#define COMPUTE_4_NROOTS(j)								\
		ASM_G (											\
			"movdqa (%%eax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%ecx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%edx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%ebx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%ecx) \n\t"			/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%edx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	

	#define COMPUTE_NEXT_ROOTS_N		\
		ASM_G (							\
			"movl %%eax, %%ecx		\n\t"	\
			"addl %2, %%ecx			\n\t"	/*r8d = root1 + ptr*/		\
			"movl %%edx, %%edi		\n\t"								\
			"addl %2, %%edi			\n\t"	/*r9d = root2 + ptr*/		\
			"subl %3, %%eax			\n\t"	/*root1 = root1 - prime*/	\
			"subl %3, %%edx			\n\t"	/*root2 = root2 - prime*/	\
			"addl %2, %%eax			\n\t"	/*root1 + ptr*/				\
			"cmovae %%ecx, %%eax	\n\t"	/*other caluclation if no overflow*/	\
			"addl %2, %%edx			\n\t"	/*root2 + ptr*/							\
			"cmovae %%edi, %%edx	\n\t"	/*other caluclation if no overflow*/	\
			: "+a"(root1), "+d"(root2)		\
			: "g"(*ptr), "g"(prime)			\
			: "ecx", "edi", "cc");

	#define COMPUTE_FIRST_ROOTS			\
		ASM_G (											\
			"xorl %%ecx, %%ecx		\n\t"	/*r8d = 0*/	\
			"xorl %%edi, %%edi		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - bmodp*/	\
			"cmovc %3, %%ecx		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - bmodp*/	\
			"cmovc %3, %%edi		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%ecx, %%eax		\n\t"		\
			"addl %%edi, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "r"(bmodp), "g"(prime)		\
			: "ecx", "edi", "cc");	

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
	update_t update_data = dconf->update_data;
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
	uint32 *bptr, *sliceptr_p, *sliceptr_n;
	uint32 *numptr_p, *numptr_n;
	int check_bound = BUCKET_ALLOC/2 - 1, room;

	numblocks = sconf->num_blocks;
	interval = numblocks << BLOCKBITS;

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->fb_bounds[0] = fb->med_B;

		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		//reset lp_buckets
		for (i=0;i< (2*numblocks*lp_bucket_p->alloc_slices) ;i++)
			numptr_p[i] = 0;

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
			update_data.firstroots1[i] = root2;
			update_data.firstroots2[i] = root1;
			fb_p[i].roots = (uint32)((root1 << 16) | root2);
			fb_n[i].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
		}
		else
		{
			update_data.firstroots1[i] = root1;
			update_data.firstroots2[i] = root2;
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
		CHECK_NEW_SLICE(i);

		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1; 
		//logp = fb->list->logprime[i];

		amodp = (int)zShortMod(&poly->poly_a,prime);
		bmodp = (int)zShortMod(&poly->poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		//fb_offset = (i - bound_val) << 16;
		root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);
		
		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;

		FILL_ONE_PRIME_LOOP_P(i);

		root1 = (prime - update_data.firstroots1[i]);
		root2 = (prime - update_data.firstroots2[i]);

		FILL_ONE_PRIME_LOOP_N(i);

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

	for (i=fb->large_B;i<fb->B;i++)
	{
		CHECK_NEW_SLICE(i);

		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1; 
		//logp = fb->list->logprime[i];

		amodp = (int)zShortMod(&poly->poly_a,prime);
		bmodp = (int)zShortMod(&poly->poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		root1 = (uint32)((uint64)inv * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)inv * (uint64)root2 % (uint64)prime);

		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;

		FILL_ONE_PRIME_P(i);

		root1 = (prime - root1);
		root2 = (prime - root2);

		FILL_ONE_PRIME_N(i);

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

	update_t update_data = dconf->update_data;

	uint32 startprime = 2;
	uint32 bound = sconf->factor_base->B;

	char v = dconf->curr_poly->nu[dconf->numB];
	char sign = dconf->curr_poly->gray[dconf->numB];
	int *ptr;

	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 med_B = sconf->factor_base->med_B;
	uint32 large_B = sconf->factor_base->large_B;

	uint32 j, interval; //, fb_offset;
	int k,bnum,numblocks,room;
	uint32 root1, root2, prime;

	int bound_index=0;
	int check_bound = BUCKET_ALLOC/2 - 1;
	uint32 bound_val = med_B;
	uint32 *bptr, *sliceptr_p,*sliceptr_n;
	uint32 *numptr_p, *numptr_n;
	uint8 logp=0;
	polysieve_t helperstruct;

	numblocks = sconf->num_blocks;
	interval = numblocks << BLOCKBITS;
	
	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->fb_bounds[0] = med_B;

		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		
		//reuse this for a sec...
		prime = 2*numblocks*lp_bucket_p->alloc_slices;

		//reset lp_buckets
		for (j=0;j<prime;j++)
			numptr_p[j] = 0;
	
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

	if (sign > 0)
	{
#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=startprime;j<med_B;j++,ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;
			
			if (root2 < root1)
			{
				update_data.firstroots1[j] = root2;
				update_data.firstroots2[j] = root1;
				fb_p[j].roots = (uint32)((root1 << 16) | root2);
				fb_n[j].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
			}
			else
			{
				update_data.firstroots1[j] = root1;
				update_data.firstroots2[j] = root2;
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
		

		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[med_B-1];

		if (med_B % 16 != 0)
		{
			printf("med_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((large_B - med_B) % 16 != 0)
		{
			printf("med range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];						//72
		helperstruct.large_B = med_B;				//80
		helperstruct.B = large_B;					//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
			"5:		\n\t"							\
				/* now we are in the room loop */ \
				/* room is in register ecx */ \
			"movl   (%%r10,%%rbx,4),%%eax \n\t"		/* value at numptr_p + k */ \
			"movl   (%%r11,%%rbx,4),%%edx \n\t"		/* value at numptr_n + k */ \
			"cmpl   %%ecx,%%eax \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%eax,%%ecx \n\t"				/* new value of room if so */ \
			"cmpl   %%ecx,%%edx \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%edx,%%ecx \n\t"				/* new value of room if so */ \
			"incq   %%rbx \n\t"						/* increment counter */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"	/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
				/* ================================================ */	\
				/* =========== ITERATION 1               ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 2               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 3               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 4               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory", "cc");

		// refresh local pointers and constants before entering the next loop
		numptr_n = helperstruct.numptr_n;
		numptr_p = helperstruct.numptr_p;
		sliceptr_n = helperstruct.sliceptr_n;
		sliceptr_p = helperstruct.sliceptr_p;
		ptr = helperstruct.ptr;	
		bound_val = helperstruct.bound_val;	
		check_bound = helperstruct.check_bound;
		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;


#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_PROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;
		}


#else
		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;
			//fb_offset = (j - bound_val) << 16;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif
			
		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[large_B-1];

		if (large_B % 16 != 0)
		{
			printf("large_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((bound - large_B) % 16 != 0)
		{
			printf("large range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = large_B;				//80
		helperstruct.B = bound;						//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
			"5:		\n\t"							\
				/* now we are in the room loop */ \
				/* room is in register ecx */ \
			"movl   (%%r10,%%rbx,4),%%eax \n\t"		/* value at numptr_p + k */ \
			"movl   (%%r11,%%rbx,4),%%edx \n\t"		/* value at numptr_n + k */ \
			"cmpl   %%ecx,%%eax \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%eax,%%ecx \n\t"				/* new value of room if so */ \
			"cmpl   %%ecx,%%edx \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%edx,%%ecx \n\t"				/* new value of room if so */ \
			"incq   %%rbx \n\t"						/* increment counter */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
				/* ================================================ */	\
				/* =========== ITERATION 1               ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 2               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 3               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 4               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory", "cc");

		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;


#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_PROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);
			
			j++;
		}


#else
		logp = update_data.logp[j-1];
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

#endif

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
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.firstroots1[j] = root2;
				update_data.firstroots2[j] = root1;
				fb_p[j].roots = (uint32)((root1 << 16) | root2);
				fb_n[j].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
			}
			else
			{
				update_data.firstroots1[j] = root1;
				update_data.firstroots2[j] = root2;
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
		
		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[med_B-1];

		if (med_B % 16 != 0)
		{
			printf("med_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((large_B - med_B) % 16 != 0)
		{
			printf("med range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = med_B;				//80
		helperstruct.B = large_B;					//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
			"5:		\n\t"							\
				/* now we are in the room loop */ \
				/* room is in register ecx */ \
			"movl   (%%r10,%%rbx,4),%%eax \n\t"		/* value at numptr_p + k */ \
			"movl   (%%r11,%%rbx,4),%%edx \n\t"		/* value at numptr_n + k */ \
			"cmpl   %%ecx,%%eax \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%eax,%%ecx \n\t"				/* new value of room if so */ \
			"cmpl   %%ecx,%%edx \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%edx,%%ecx \n\t"				/* new value of room if so */ \
			"incq   %%rbx \n\t"						/* increment counter */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
				/* ================================================ */	\
				/* =========== ITERATION 1               ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 2               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 3               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 4               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"	 						/* beginning of loop */ \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d,%%r8d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r8d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%r12d, %%r9d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r9d \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movl	%%r12d,%%edx \n\t"				/* copy prime to edx */ \
			"movl	%%r12d,%%eax \n\t"				/* copy prime to eax */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r12d,%%r13d \n\t"				/* copy root to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%r12d \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%r12d \n\t"		/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root2 > interval? */ \
			"jae	6f \n\t"						/* jump out if so */ \
			"1:		\n\t"							/* beginning of do loop */ \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%edx,%%r13d \n\t"				/* copy root1 to r13 */ \
			"andl   $" BLOCKSIZEm1txt ",%%r13d \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%r13d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"addl	%%eax, %%edx \n\t"				/* increment root by prime */ \
			"cmpl   88(%%rsi,1),%%edx \n\t"			/* root > interval? */ \
			"jb		1b \n\t"						/* repeat if necessary */ \
			"6:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory", "cc");

		// refresh local pointers and constants before entering the next loop
		numptr_n = helperstruct.numptr_n;
		numptr_p = helperstruct.numptr_p;
		sliceptr_n = helperstruct.sliceptr_n;
		sliceptr_p = helperstruct.sliceptr_p;
		ptr = helperstruct.ptr;	
		bound_val = helperstruct.bound_val;	
		check_bound = helperstruct.check_bound;
		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_NROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;
		}


#else

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[large_B-1];
		
		if (large_B % 16 != 0)
		{
			printf("large_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((bound - large_B) % 16 != 0)
		{
			printf("large range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = large_B;				//80
		helperstruct.B = bound;						//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
			"5:		\n\t"							\
				/* now we are in the room loop */ \
				/* room is in register ecx */ \
			"movl   (%%r10,%%rbx,4),%%eax \n\t"		/* value at numptr_p + k */ \
			"movl   (%%r11,%%rbx,4),%%edx \n\t"		/* value at numptr_n + k */ \
			"cmpl   %%ecx,%%eax \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%eax,%%ecx \n\t"				/* new value of room if so */ \
			"cmpl   %%ecx,%%edx \n\t"				/* *(numptr_p + k) > room ? */ \
			"cmova  %%edx,%%ecx \n\t"				/* new value of room if so */ \
			"incq   %%rbx \n\t"						/* increment counter */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t" /* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
				/* ================================================ */	\
				/* =========== ITERATION 1               ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 2               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"		/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 3               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== ITERATION 4               ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm0,%%r12d \n\t"				/* extract prime from xmm0 */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"movd	%%xmm1,%%r8d \n\t"				/* extract root1,1 from xmm1 */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root2 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* extract root2,1 from xmm2 */ \
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
				/* all code paths at this point have:	*/ \
				/* root1		-> r8d					*/ \
				/* root2		-> r9d					*/ \
				/* r12			-> prime				*/ \
				/* r15			 -> j					*/ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    1f \n\t" \
			"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    3f \n\t" \
			"movl   %%r9d,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"movl   %%r9d,%%eax \n\t"				/* mov root2 to eax */ \
			"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== INVERT ROOTS FOR NEG SIDE ========== */	\
				/* ================================================ */	\
			"movl	%%r12d,%%edx \n\t"				/* copy prime to ebx */ \
			"subl   %%r8d,%%r12d \n\t"				/* root1 (r12d) = (prime - root1);	*/	\
			"subl   %%r9d,%%edx \n\t"				/* root2 (edx) = (prime - root2);	*/	\
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%r12d \n\t"				/* root1 > interval? */ \
			"jae    5f \n\t" \
			"movl   %%r12d,%%ebx \n\t"				/* copy root1 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%r12d \n\t"	/* root1 & BLOCKSIZEm1 */ \
			"orl	%%r12d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2 ========== */	\
				/* ================================================ */	\
			"cmpl   %%r13d,%%edx \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" \
			"movl   %%edx,%%ebx \n\t"				/* copy root2 to ebx */ \
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
			"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
			"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
			"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_n[bnum] */ \
			"subl   96(%%rsi,1),%%edi \n\t"			/* j - bound_val */ \
			"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
			"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
			"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_n[bnum] */ \
			"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
			"andl   $" BLOCKSIZEm1txt ",%%edx \n\t"	/* root2 & BLOCKSIZEm1 */ \
			"orl	%%edx,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
			"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
			"7:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $1,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory", "cc");

		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_NROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);
			
			j++;
		}


#else

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound;j++,ptr++)				
		{				
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_N;		

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
			FILL_ONE_PRIME_N(j);
		}

#endif

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

	return;
}
