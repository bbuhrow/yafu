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

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

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


#if defined(MSC_ASM32A)

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

	#define COMPUTE_FIRST_ROOTS		\
		root1 = (int)root1 - bmodp;		\
		root2 = (int)root2 - bmodp;		\
		if (root1 < 0) root1 += prime;			\
		if (root2 < 0) root2 += prime;

#endif

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

	for (i=start_prime;i<sconf->sieve_small_fb_start;i++)
	{
		uint64 q64, tmp, t2;

		prime = fb->tinylist->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		// reuse integer inverse of prime that we've calculated for use
		// in trial division stage
		// inv * root1 % prime
		t2 = (uint64)inv * (uint64)root1;
		tmp = t2 + (uint64)fb->tinylist->correction[i];
		q64 = tmp * (uint64)fb->tinylist->small_inv[i];
		tmp = q64 >> 32; 
		root1 = t2 - tmp * prime;

		// inv * root2 % prime
		t2 = (uint64)inv * (uint64)root2;
		tmp = t2 + (uint64)fb->tinylist->correction[i];
		q64 = tmp * (uint64)fb->tinylist->small_inv[i];
		tmp = q64 >> 32; 
		root2 = t2 - tmp * prime;
	
		//we don't sieve these primes, so ordering doesn't matter
		update_data.firstroots1[i] = root1;
		update_data.firstroots2[i] = root2;
#ifdef USE_COMPRESSED_FB
		fb_p[i].roots = (uint32)((root1 << 16) | root2);
		if (root2 == 0)
			fb_n[i].roots = 0;
		else
			fb_n[i].roots = (uint32)((prime - root2) << 16);
		if (root1 == 0)
			fb_n[i].roots |= 0;
		else
			fb_n[i].roots |= (uint32)(prime - root1);
#else
		fb_p->root1[i] = (uint16)root1;
		fb_p->root2[i] = (uint16)root2;
		fb_n->root1[i] = (uint16)(prime - root2);
		fb_n->root2[i] = (uint16)(prime - root1);
		//if we were sieving, this would double count the location on the 
		//positive side.  but since we're not, its easier to check for inclusion
		//on the progression if we reset the negative root to zero if it is == prime
		if (fb_n->root1[i] == prime)
			fb_n->root1[i] = 0;
		if (fb_n->root2[i] == prime)
			fb_n->root2[i] = 0;
#endif

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j],prime);
			
			// x * inv % prime
			t2 = (uint64)inv * (uint64)x;
			tmp = t2 + (uint64)fb->tinylist->correction[i];
			q64 = tmp * (uint64)fb->tinylist->small_inv[i];
			tmp = q64 >> 32; 
			x = t2 - tmp * prime;

			rootupdates[(j)*fb->B+i] = x;
		}
	}

	for (i=sconf->sieve_small_fb_start;i<fb->med_B;i++)
	{
		uint64 small_inv, correction;
		uint64 q64, tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		// compute integer inverse of prime for use in mod operations in this
		// function.
		small_inv = ((uint64)1 << 48) / (uint64)prime;
		if (floor((double)((uint64)1 << 48) / (double)prime + 0.5) ==
						(double)small_inv) {
			correction = 1;
		}
		else {
			correction = 0;
			small_inv++;
		}

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS

		// inv * root1 % prime
		t2 = (uint64)inv * (uint64)root1;
		tmp = t2 + correction;
		q64 = tmp * small_inv;
		tmp = q64 >> 48; 
		root1 = t2 - tmp * prime;

		// inv * root2 % prime
		t2 = (uint64)inv * (uint64)root2;
		tmp = t2 + correction;
		q64 = tmp * small_inv;
		tmp = q64 >> 48; 
		root2 = t2 - tmp * prime;

		if (root2 < root1)
		{
			update_data.firstroots1[i] = root2;
			update_data.firstroots2[i] = root1;
#ifdef USE_COMPRESSED_FB
			fb_p[i].roots = (uint32)((root1 << 16) | root2);
			fb_n[i].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
#else
			fb_p->root1[i] = (uint16)root2;
			fb_p->root2[i] = (uint16)root1;
			fb_n->root1[i] = (uint16)(prime - root1);
			fb_n->root2[i] = (uint16)(prime - root2);
#endif
		}
		else
		{
			update_data.firstroots1[i] = root1;
			update_data.firstroots2[i] = root2;
#ifdef USE_COMPRESSED_FB
			fb_p[i].roots = (uint32)((root2 << 16) | root1);
			fb_n[i].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
#else
			fb_p->root1[i] = (uint16)root1;
			fb_p->root2[i] = (uint16)root2;
			fb_n->root1[i] = (uint16)(prime - root2);
			fb_n->root2[i] = (uint16)(prime - root1);
#endif
		}

		//for this factor base prime, compute the rootupdate value for all s
		//Bl values.  amodp holds a^-1 mod p
		//the rootupdate value is given by 2*Bj*amodp
		//Bl[j] now holds 2*Bl
		for (j=0;j<s;j++)
		{
			x = (int)mpz_tdiv_ui(dconf->Bl[j],prime);

			// x * inv % prime
			t2 = (uint64)inv * (uint64)x;
			tmp = t2 + correction;
			q64 = tmp * small_inv;
			tmp = q64 >> 48; 
			x = t2 - tmp * prime;

			rootupdates[(j)*fb->B+i] = x;
		}
	}

	check_bound = fb->med_B + BUCKET_ALLOC/2;
	logp = fb->list->logprime[fb->med_B-1];
	for (i=fb->med_B;i<fb->large_B;i++)
	{
		//uint64 small_inv, correction;
		//uint64 q64, tmp, t2;

		CHECK_NEW_SLICE(i);

		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1; 

		//small_inv = ((uint64)1 << 40) / (uint64)prime;
		//if (floor((double)((uint64)1 << 40) / (double)prime + 0.5) ==
		//				(double)small_inv) {
		//	correction = 1;
		//}
		//else {
		//	correction = 0;
		//	small_inv++;
		//}

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

		//find a^-1 mod p = inv(a mod p) mod p
		inv = modinv_1(amodp,prime);

		COMPUTE_FIRST_ROOTS
	
		//t2 = (uint64)inv * (uint64)root1;
		//tmp = t2 + correction;
		//q64 = tmp * small_inv;
		//tmp = q64 >> 40; 
		//root1 = t2 - tmp * prime;

		//if ((t2 % (uint64)prime) != root1)
		//{
		//	printf("failed for prime = %u, amodpinv = %u\n",
		//		prime, inv);
		//	exit(1);
		//}

		//t2 = (uint64)inv * (uint64)root2;
		//tmp = t2 + correction;
		//q64 = tmp * small_inv;
		//tmp = q64 >> 40; 
		//root2 = t2 - tmp * prime;

		//if ((t2 % (uint64)prime) != root2)
		//{
		//	printf("failed for prime = %u, amodpinv = %u\n",
		//		prime, inv);
		//	exit(1);
		//}

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
			x = (int)mpz_tdiv_ui(dconf->Bl[j], prime);
			x = (int)((int64)x * (int64)inv % (int64)prime);

			//t2 = (uint64)inv * (uint64)x;
			//tmp = t2 + correction;
			//q64 = tmp * small_inv;
			//tmp = q64 >> 48; 
			//x = t2 - tmp * prime;

			rootupdates[(j)*fb->B+i] = x;
		}

	}

	logp = fb->list->logprime[fb->large_B-1];
	for (i=fb->large_B;i<fb->B;i++)
	{
		CHECK_NEW_SLICE(i);

		prime = fb->list->prime[i];
		root1 = modsqrt[i];
		root2 = prime - root1; 

		amodp = (int)mpz_tdiv_ui(poly->mpz_poly_a,prime);
		bmodp = (int)mpz_tdiv_ui(poly->mpz_poly_b,prime);

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
			x = (int)mpz_tdiv_ui(dconf->Bl[j], prime);
			x = (int)((int64)x * (int64)inv % (int64)prime);
			rootupdates[(j)*fb->B+i] = x;
		}
	}

	if (lp_bucket_p->list != NULL)
		lp_bucket_p->num_slices = bound_index + 1;
	

	return;
}

