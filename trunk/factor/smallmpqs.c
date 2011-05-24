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
#include "factor.h"
#include "soe.h"
#include "util.h"
#include "common.h"

typedef struct
{
	uint32 prime_and_logp;
	uint32 roots;					//root1 is stored in the lower 16 bits, root2 in the upper 16
} smpqs_sieve_fb;

static void smpqs_sieve_block(uint8 *sieve, smpqs_sieve_fb *fb, uint32 start_prime, 
	uint8 s_init, fb_list *fullfb);

static int smpqs_check_relations(uint32 sieve_interval, uint32 blocknum, uint8 sieve[] ,z *n, 
								mpqs_poly *poly, uint8 closnuf,
								smpqs_sieve_fb *fb,fb_list *fullfb, mpqs_rlist *full, 
								mpqs_rlist *partial,uint32 cutoff,
								uint8 small_cutoff, uint32 start_prime, uint32 parity, 
								uint32 *num, int numpoly);

static void smpqs_trial_divide_Q(z *Q,smpqs_sieve_fb *fb,mpqs_rlist *full,mpqs_rlist *partial,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, fb_list *fullfb, 
						  uint32 cutoff, uint8 small_cutoff, uint32 start_prime, int numpoly, 
						  uint32 parity, uint8 closnuf, uint8 bits);

static void smpqs_save_relation(mpqs_rlist *list, uint32 offset, uint32 largeprime, uint32 num_factors, 
						  uint32 rnum, uint16 *fboffset, int numpoly, uint32 parity);

void smpqs_make_fb_mpqs(fb_list *fb, uint32 *modsqrt, z *n);
void smpqs_computeB(mpqs_poly *poly, uint32 polyd, z *n);
void smpqs_nextD(uint32 *polyd, uint32 *polyd_id, z *n);
void smpqs_computeRoots(mpqs_poly *poly, fb_list *fb, uint32 *modsqrt, 
	smpqs_sieve_fb *fbp, smpqs_sieve_fb *fbn, uint32 start_prime);
int qcomp_smpqs(const void *x, const void *y);
uint8 smpqs_choose_multiplier(z *n, uint32 fb_size);
int smpqs_BlockGauss(mpqs_rlist *full, mpqs_rlist *partial, z *apoly, z *bpoly,
			fb_list *fb, z *n, int mul, 
			z *factors,uint32 *num_factor);

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	//these compilers support SIMD 
	#define SIMD_SIEVE_SCAN 1
	#define SCAN_CLEAN asm volatile("emms");	
	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			asm volatile (		\
				"movdqa (%1), %%xmm0   \n\t"		\
				"orpd 16(%1), %%xmm0    \n\t"		\
				"pmovmskb %%xmm0, %0   \n\t"		\
				: "=r"(result)						\
				: "r"(sieveblock + j), "0"(result)	\
				: "%xmm0");

	#else
		#define SCAN_16X	\
			result = 1;	/*dont know what compiler this is. force the normal method*/
		#undef SIMD_SIEVE_SCAN
	#endif

#elif defined(MSC_ASM32A)
	#define SIMD_SIEVE_SCAN 1
	#define SCAN_CLEAN ASM_M {emms};
	#define SSE2_RESIEVING 1

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
					ASM_M pmovmskb ecx, xmm0			\
					ASM_M mov result, ecx};			\
			} while (0);

	#else
		#undef SIMD_SIEVE_SCAN
	#endif

#elif defined(_WIN64)

	#define SIMD_SIEVE_SCAN 1
	#define SSE2_RESIEVING 1
	#define SCAN_CLEAN /*nothing*/

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			do	{						\
				__m128i local_block;	\
				__m128i local_block2;	\
				local_block = _mm_load_si128(sieveblock + j); \
				local_block2 = _mm_load_si128(sieveblock + j + 2); \
				local_block = _mm_or_si128(local_block, local_block2); \
				result = _mm_movemask_epi8(local_block); \
			} while (0);


	#else
		#undef SIMD_SIEVE_SCAN
	#endif

#endif

#define SCAN_MASK 0x8080808080808080ULL


void smpqs_make_fb_mpqs(fb_list *fb, uint32 *modsqrt, z *n)
{
	//finds the factor base primes, and computes the solutions to the congruence x^2 = N mod p
	//for the QS, these are the starting positions of the sieve relative to the sqrt of N.
	//for the MPQS, additional work using the polynomial coefficents and these congruences 
	//needs to be done to compute the starting positions of the sieve.

	int i;
	uint32 b,j,r,k;
	uint32 prime, root1, root2;
	fp_digit f;
	uint8 logp;

	//the 0th element in the fb is always  2, so start searching with 3
	j=2; i=1;
	while (j<fb->B)
	{
		r = (uint32)zShortMod(n,(fp_digit)spSOEprimes[i]);
		if (r == 0)
		{
			//p divides n, which means it divides the multiplier.
			//we can still use it, but it only has one solution to x^2 == n mod p instead
			//of two.  just divide its logprime in half.
			//we also can't find the root using shanks-tonelli, but it will be very small
			//because the multiplier is very small, so just use brute force.
			prime = (uint32)spSOEprimes[i];
			b = (uint32)zShortMod(n,(fp_digit)prime);
			k=0;
			while (1)
			{
				if (((k*k) % prime) == b)
					break;
				k++;
			}
			root1 = k;
			root2 = prime - k;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5)/2;

			//fill in factor base
			fb->list->prime[j] = prime;
			modsqrt[j] = root1;
			fb->list->logprime[j] = logp;

			//store a couple things so we can replace single precision
			//mods with shifts and muls during trial division
			if (prime < 256)
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 32) / (uint64)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			else
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 40) / (uint64)prime);
				if (floor(256 * MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			
			j++;
			i++;
			continue;
		}

		b = jacobi_1((fp_digit)r,(fp_digit)spSOEprimes[i]);
		if (b==1)
		{
			//this prime works
			prime = (uint32)spSOEprimes[i];
			ShanksTonelli_1((fp_digit)r,(fp_digit)prime,&f);
			root1 = (uint32)f;
			root2 = prime - root1;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5);

			//fill in factor base
			fb->list->prime[j] = prime;
			modsqrt[j] = root1;
			fb->list->logprime[j] = logp;

			//store a couple things so we can replace single precision
			//mods with shifts and muls during trial division
			if (prime < 256)
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 32) / (uint64)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			else
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 40) / (uint64)prime);
				if (floor(256 * MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			
			j++;
		}
		i++;
	}

	return;
}

void smallmpqs(fact_obj_t *fobj)
{
	z *n = &fobj->qs_obj.n;
	mpqs_rlist *full, *partial;
	fb_list *fb;
	smpqs_sieve_fb *fb_sieve_p,*fb_sieve_n;
	mpqs_poly *poly;
	uint32 *modsqrt;
	
	z *factors;
	z tmp, tmp2, tmp3, sqrt_n;
	z *apoly, *bpoly;

	double t_time;
	struct timeval tstart, tend;
	TIME_DIFF *	difference;
	fp_digit fpt;
	uint32 numpoly, polyalloc, polyd, polyd_id;
	uint32 mul,i,j, j2;
	uint32 pmax;							//largest prime in factor base
	uint32 cutoff;
	uint32 sieve_interval;
	uint32 start_prime;
	uint32 num, max_f;
	int digits_n,charcount;
	uint32 num_factors;
	uint8 *sieve;							//sieve values
	uint8 s_init;							//initial sieve value
	uint8 closnuf, small_bits, max_bits;

	if (n->val[0] == 1 && n->size == 1)
		return;

	if ((n->val[0] % 2) == 0)
	{
		printf("n must be odd\n");
		return;
	}

	//if (isPrime(n))
	//{
	//	n->type = PRP;
	//	zCopy(&zOne,n);
	//	return;
	//}

	zInit(&tmp);
	
	// factor base bound
	i = zBits(n);
	//crude tuning of sieve interval based on digits in n
	if (i < 50)
	{
		j = sp_shanks_loop(n,fobj);	
		fobj->qs_obj.num_factors += 2;
		fobj->qs_obj.factors = (z *)realloc(fobj->qs_obj.factors, 
			fobj->qs_obj.num_factors * sizeof(z));
		zInit(&fobj->qs_obj.factors[fobj->qs_obj.num_factors - 2]);
		sp2z(j, &fobj->qs_obj.factors[fobj->qs_obj.num_factors - 2]);

		zShortDiv(n,j,&tmp);
		zInit(&fobj->qs_obj.factors[fobj->qs_obj.num_factors - 1]);
		zCopy(&tmp, &fobj->qs_obj.factors[fobj->qs_obj.num_factors - 1]);

		zFree(&tmp);
		return;
	}
	else if (i < 60)
	{
		j = 40;
		sieve_params.num_blocks = 1;
	}
	else if (i < 70)
	{
		j = 55;
		sieve_params.num_blocks = 1;
	}
	else if (i < 80)
	{
		j = 75;
		sieve_params.num_blocks = 1;
	}
	else if (i < 90)
	{
		j = 95;
		sieve_params.num_blocks = 2;
	}
	else if (i < 100)
	{
		j = 160;
		sieve_params.num_blocks = 2;
	}
	else if (i < 120)
	{
		j = 250;
		sieve_params.num_blocks = 3;
	}
	else
	{
		printf("input too big\n");
		zFree(&tmp);
		return;
	}

	gettimeofday(&tstart, NULL);
	zInit(&tmp2);
	zInit(&tmp3);
	
	//default mpqs parameters
	sieve_params.fudge_factor = 1.3;
	sieve_params.large_mult = 50;
	sieve_params.num_extra_relations = 16;

	//allocate the space for the factor base
	fb = (fb_list *)malloc(sizeof(fb_list));
	
	//set fb size from above
	fb->B = j;

	//compute the number of digits in n 
	digits_n = ndigits(n);

	//allocate storage for relations based on the factor base size
	max_f = fb->B + 3*sieve_params.num_extra_relations;	
	full = (mpqs_rlist *)malloc((size_t)(sizeof(mpqs_rlist)));
	full->allocated = max_f;
	full->num_r = 0;
	full->act_r = 0;
	full->list = (mpqs_r **)malloc((size_t) (max_f * sizeof(mpqs_r *)));

	//we will typically also generate max_f/2 * 10 partials (empirically determined)
	partial = (mpqs_rlist *)malloc((size_t)(sizeof(mpqs_rlist)));
	partial->allocated = 10*fb->B;
	partial->num_r = 0;
	partial->act_r = 0;
	partial->list = (mpqs_r **)malloc((size_t) (10*fb->B* sizeof(mpqs_r *)));

	//set the sieve interval.  this depends on the size of n, but for now, just fix it.  as more data
	//is gathered, use some sort of table lookup.
	sieve_interval = BLOCKSIZE*sieve_params.num_blocks;

	//allocate the space for the factor base
	modsqrt = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list = (fb_element_siqs *)malloc((size_t)(sizeof(fb_element_siqs)));
	fb->list->correction = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list->prime = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list->small_inv = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list->logprime = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb_sieve_p = (smpqs_sieve_fb *)malloc((size_t)(fb->B * sizeof(smpqs_sieve_fb)));
	fb_sieve_n = (smpqs_sieve_fb *)malloc((size_t)(fb->B * sizeof(smpqs_sieve_fb)));
	
	//allocate the sieve
	sieve = (uint8 *)malloc((size_t) (BLOCKSIZE * sizeof(uint8)));

	//allocate the current polynomial
	poly = (mpqs_poly *)malloc(sizeof(mpqs_poly));
	zInit(&poly->poly_a);
	zInit(&poly->poly_b);
	zInit(&poly->poly_c);
	zInit(&poly->poly_d);

	//allocate the polynomial lists
	polyalloc = 32;
	apoly = (z *)malloc(polyalloc * sizeof(z));
	bpoly = (z *)malloc(polyalloc * sizeof(z));
	for (i=0;i<polyalloc;i++)
	{
		zInit(&apoly[i]);
		zInit(&bpoly[i]);
	}

	//find multiplier
	mul = (uint32)smpqs_choose_multiplier(n,fb->B);
	zShortMul(n,mul,&tmp);
	zCopy(&tmp,n);

	//find new sqrt_n
	zInit(&sqrt_n);
	zNroot(n,&sqrt_n,2);

	//find upper bound of Q values
	zShiftRight(&tmp,n,1);
	zNroot(&tmp,&tmp2,2);
	zShortMul(&tmp2,sieve_interval,&tmp);
	max_bits = (uint8)zBits(&tmp);
	
	//compute the first polynominal 'a' value.  we'll need it before creating the factor base in order
	//to find the first roots
	//'a' values should be as close as possible to sqrt(2n)/M, they should be a quadratic residue mod N (d/N) = 1,
	//and be a prime congruent to 3 mod 4.  this last requirement is so that b values can be computed without using the 
	//shanks-tonelli algorithm, and instead use faster methods.
	//since a = d^2, find a d value near to sqrt(sqrt(2n)/M)
	zShiftLeft(&tmp,n,1);
	zNroot(&tmp,&tmp2,2);
	zShortDiv(&tmp2,sieve_interval,&tmp2);

	//make sure the first primes up to 10M are computed
	if (P_MAX < 10000000)
		GetPRIMESRange(0,10001000);

	zNroot(&tmp2,&poly->poly_d,2);
	if (!(poly->poly_d.val[0] & 1))
		zAdd(&poly->poly_d,&zOne,&poly->poly_d);

	polyd = poly->poly_d.val[0];
	zNextPrime_1(polyd, &fpt, &tmp, 1);
	polyd = fpt;

	for (i=0; i<NUM_P; i++)
		if (polyd == PRIMES[i])
			break;
	polyd_id = i;

	//alternating poly_a values to either side of the ideal poly_a makes the average poly_a
	//closer to ideal, and hence more effective at producing relations.  I've verified that the poly_a's produced
	//are about half the size of those going in only one direction.
	//in practice though, i don't see much of a difference in timing.
	//this is because the percent difference to the ideal poly_a is tiny in either case.  
	//printf("computing next D from %s\n",z2decstr(&poly->poly_d,&gstr1));
	smpqs_nextD(&polyd,&polyd_id,n);
	//printf("computeB\n");
	smpqs_computeB(poly,polyd,n);

	fb->list->prime[0] = 1;
	fb->list->prime[1] = 2;

	//construct the factor base, and copy to the sieve factor base
	smpqs_make_fb_mpqs(fb,modsqrt,n);

	for (i=2;i<fb->B;i++)
	{
		fb_sieve_p[i].prime_and_logp = (fb->list->prime[i] << 16) | (fb->list->logprime[i]);
		fb_sieve_n[i].prime_and_logp = (fb->list->prime[i] << 16) | (fb->list->logprime[i]);
	}

	//find the root locations of the factor base primes for this poly
	smpqs_computeRoots(poly,fb,modsqrt,fb_sieve_p,fb_sieve_n,2);

	pmax = fb->list->prime[fb->B-1];
	cutoff = pmax * sieve_params.large_mult;

	//compute the number of bits in M/2*sqrt(N/2), the approximate value
	//of residues in the sieve interval
	//sieve locations greater than this are worthy of trial dividing
	closnuf = (uint8)(double)((zBits(n) - 1)/2);
	closnuf += (uint8)(log((double)sieve_interval/2)/log(2.0));
	closnuf -= (uint8)(sieve_params.fudge_factor * log(cutoff) / log(2.0));
	
	closnuf += 6;

	//small prime variation -- hand tuned small_bits correction (likely could be better)
	small_bits = 7;
	closnuf -= small_bits;
	start_prime = 7;
	
	s_init = closnuf;

	//print some info to the screen and the log file
	if (VFLAG > 0)
	{
		printf("n = %d digits, %d bits\n",digits_n,zBits(n));
		printf("==== sieve params ====\n");
		printf("factor base: %d primes (max prime = %u)\n",fb->B,pmax);
		printf("large prime cutoff: %u (%d * pmax)\n",cutoff,sieve_params.large_mult);
		printf("sieve interval: %d blocks of size %d\n",sieve_interval/BLOCKSIZE,BLOCKSIZE);
		printf("multiplier is %u\n",mul);
		printf("trial factoring cutoff at %d bits\n",closnuf);
		printf("==== sieving in progress ====\n");
	}

	numpoly = 0;
	num = 0;
	charcount=0;
	while (1)
	{
		//copy current poly into the poly lists
		if (numpoly < polyalloc)
		{
			zCopy(&poly->poly_a,&apoly[numpoly]);
			zCopy(&poly->poly_b,&bpoly[numpoly]);
		}
		else
		{
			// get more space for the polys, if needed
			polyalloc *= 2;
			apoly = (z *)realloc(apoly,polyalloc * sizeof(z));
			bpoly = (z *)realloc(bpoly,polyalloc * sizeof(z));
			for (i=numpoly; i<polyalloc; i++)
			{
				zInit(&apoly[i]);
				zInit(&bpoly[i]);
			}
			zCopy(&poly->poly_a,&apoly[numpoly]);
			zCopy(&poly->poly_b,&bpoly[numpoly]);
		}

		for (j2=0; j2 < sieve_params.num_blocks; j2++)
		{
			//the minimum value for the current poly_a and poly_b occur at offset (-b + sqrt(N))/a
			//compute minimum bits for Q
			//zSub(&sqrt_n,&poly->poly_b,&tmp);
			//zDiv(&tmp,&poly->poly_a,&tmp2,&tmp3);
			smpqs_sieve_block(sieve,fb_sieve_p,start_prime,s_init,fb);

			i = smpqs_check_relations(sieve_interval,j2,sieve,n,poly,
				s_init,fb_sieve_p,fb,full,partial,cutoff,small_bits,
				start_prime,0,&num,numpoly);

			smpqs_sieve_block(sieve,fb_sieve_n,start_prime,s_init,fb);
		
			i = smpqs_check_relations(sieve_interval,j2,sieve,n,poly,
				s_init,fb_sieve_n,fb,full,partial,cutoff,small_bits,
				start_prime,1,&num,numpoly);
		
			if (partial->num_r > 0)
			{
				//check the partials for full relations
				qsort(partial->list,partial->num_r,sizeof(mpqs_r *),&qcomp_smpqs);
				j=0;
				for (i=0;i<partial->num_r-1;i++)
				{
					if (partial->list[i]->largeprime == partial->list[i+1]->largeprime)
						j++;
				}
				partial->act_r = j;
				
				if (j+(full->num_r) >= fb->B + sieve_params.num_extra_relations) 
				{
					//we've got enough total relations to stop
					goto done;
				}
			}
		}

		//next polynomial
		smpqs_nextD(&polyd,&polyd_id,n);
		smpqs_computeB(poly,polyd,n);
		smpqs_computeRoots(poly,fb,modsqrt,fb_sieve_p,fb_sieve_n,2);

		numpoly++;
	}

done:

	if (VFLAG > 0)
		printf("%d relations found: %d full + %d from %d partial, using %d polys\n",
			partial->act_r+full->num_r,full->num_r,partial->act_r,partial->num_r,numpoly);

	gettimeofday (&tend, NULL);
	difference = my_difftime (&tstart, &tend);

	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (VFLAG > 0)
		printf("QS elapsed time = %6.4f seconds.\n",t_time);

	//can free sieving structures now
	free(sieve);
	free(fb_sieve_p);
	free(fb_sieve_n);
	zFree(&poly->poly_a);
	zFree(&poly->poly_b);
	zFree(&poly->poly_c);
	zFree(&poly->poly_d);
	free(poly);
	free(modsqrt);

	num_factors=0;
	factors = (z *)malloc(MAX_FACTORS * sizeof(z));
	for (i=0;i<MAX_FACTORS;i++)
		zInit(&factors[i]);

	gettimeofday(&tstart,NULL);
	i = smpqs_BlockGauss(full,partial,apoly,bpoly,fb,n,mul,
		factors,&num_factors);

	gettimeofday (&tend, NULL);
	difference = my_difftime (&tstart, &tend);

	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);
	
	if (VFLAG > 0)
		printf("Gauss elapsed time = %6.4f seconds.\n",t_time);

	fobj->qs_obj.num_factors += num_factors;
	fobj->qs_obj.factors = (z *)realloc(fobj->qs_obj.factors, 
		fobj->qs_obj.num_factors * sizeof(z));

	for(i=0;i<num_factors;i++)
	{
		zInit(&fobj->qs_obj.factors[fobj->qs_obj.num_factors - i - 1]);
		zCopy(&factors[i], &fobj->qs_obj.factors[fobj->qs_obj.num_factors - i - 1]);
	}

	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&sqrt_n);

	for (i=0;i<MAX_FACTORS;i++)
		zFree(&factors[i]);
	free(factors);

	for (i=0;i<full->num_r;i++)
	{
		free(full->list[i]->fboffset);
		free(full->list[i]);
	}
	free(full->list);
	free(full);

	for (i=0;i<partial->num_r;i++)
	{
		free(partial->list[i]->fboffset);
		free(partial->list[i]);
	}
	free(partial->list);
	free(partial);

	for (i=0;i<polyalloc;i++)
	{
		zFree(&apoly[i]);
		zFree(&bpoly[i]);
	}
	free(apoly);
	free(bpoly);

	free(fb->list->correction);
	free(fb->list->prime);
	free(fb->list->logprime);
	free(fb->list->small_inv);
	free(fb->list);
	free(fb);

	return;
}

static uint8 smpqs_mult_list[] =
	{1, 2, 3, 5, 7, 10, 11, 13, 15, 17, 19, 
	 23, 26, 29, 30, 31, 43, 59, 67, 73};

uint8 smpqs_choose_multiplier(z *n, uint32 fb_size) 
{
	uint32 i, j;
	uint32 num_primes = MIN(2 * fb_size, 30);
	double best_score;
	uint8 best_mult;
	double scores[20];
	uint32 num_multipliers;
	double log2n = zlog(n);

	/* measure the contribution of 2 as a factor of sieve
	   values. The multiplier itself must also be taken into
	   account in the score. scores[i] is the correction that
	   is implicitly applied to the size of sieve values for
	   multiplier i; a negative score makes sieve values 
	   smaller, and so is better */

	for (i = 0; i < 20; i++) {
		uint8 curr_mult = smpqs_mult_list[i];
		uint8 knmod8 = (uint8)((curr_mult * n->val[0]) % 8);
		double logmult = log((double)curr_mult);

		/* only consider multipliers k such than
		   k*n will not overflow an mp_t */

		if (log2n + logmult > (32 * MAX_DIGITS - 2) * LN2)
			break;

		scores[i] = 0.5 * logmult;
		switch (knmod8) {
		case 1:
			scores[i] -= 2 * LN2;
			break;
		case 5:
			scores[i] -= LN2;
			break;
		case 3:
		case 7:
			scores[i] -= 0.5 * LN2;
			break;
		/* even multipliers start with a handicap */
		}
	}
	num_multipliers = i;

	/* for the rest of the small factor base primes */
	for (i = 1; i < num_primes; i++) {
		uint32 prime = (uint32)spSOEprimes[i];
		double contrib = log((double)prime) / (prime - 1);
		uint32 modp = (uint32)zShortMod(n, prime);

		for (j = 0; j < num_multipliers; j++) {
			uint8 curr_mult = smpqs_mult_list[j];
			uint32 knmodp = (modp * curr_mult) % prime;

			/* if prime i is actually in the factor base
			   for k * n ... */

			if (knmodp == 0 || jacobi_1(knmodp, prime) == 1) {

				/* ...add its contribution. A prime p con-
				   tributes log(p) to 1 in p sieve values, plus
				   log(p) to 1 in p^2 sieve values, etc. The
				   average contribution of all multiples of p 
				   to a random sieve value is thus

				   log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
				   = (log(p) / p) * 1 / (1 - (1/p)) 
				   = log(p) / (p-1)

				   This contribution occurs once for each
				   square root used for sieving. There are two
				   roots for each factor base prime, unless
				   the prime divides k*n. In that case there 
				   is only one root */

				//printf("scores[%d] = %f\n",j,contrib);
				if (knmodp == 0)
					scores[j] -= contrib;
				else
					scores[j] -= 2 * contrib;
			}
		}

	}

	/* use the multiplier that generates the best score */
	best_score = 1000.0;
	best_mult = 1;
	for (i = 0; i < num_multipliers; i++) {
		
		double score = scores[i];
		if (score < best_score) {
			best_score = score;
			best_mult = smpqs_mult_list[i];
		}
	}
	return best_mult;
}

static int smpqs_check_relations(uint32 sieve_interval, uint32 blocknum, uint8 sieve[], z *n, mpqs_poly *poly, uint8 closnuf,
						smpqs_sieve_fb *fb,fb_list *fullfb, mpqs_rlist *full, mpqs_rlist *partial, 
						uint32 cutoff, uint8 small_cutoff, uint32 start_prime, uint32 parity, uint32 *num, int numpoly)
{
	z Q,t1,t2,t3;
	uint32 offset,i,j,k;
	uint32 neg;
	uint64 *sieveblock;
	uint32 limit = BLOCKSIZE >> 3;
	
	sieveblock = (uint64 *)sieve;
	zInit(&Q);
	zInit(&t1);
	zInit(&t2);
	zInit(&t3);

	//check for relations
	for (j=0;j<limit;j+=4)	
	{

#ifdef SIMD_SIEVE_SCAN

		uint32 result = 0;

		SIEVE_SCAN_32;

		if (result == 0)
			continue;

#else
		uint64 mask = SCAN_MASK;

		if (((sieveblock[j] | sieveblock[j+1] | 
			sieveblock[j+2] | sieveblock[j+3]) & 
		      mask) == (uint64)(0))
			continue;
#endif
	
#if defined(SIMD_SIEVE_SCAN)
		// make it safe to perform floating point
		SCAN_CLEAN;
#endif

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<4; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64)(0))
				continue;

			//at least one passed the check, find which one(s) and pass to 
			//trial division stage
			for (k=0;k<8;k++)
			{
				uint32 thisloc = ((j+i)<<3) + k;
				if ((sieve[thisloc] & 0x80) == 0)
					continue;

				(*num)++;

				offset = (blocknum<<BLOCKBITS) + thisloc;
				zShiftLeft_1(&t2,&poly->poly_b);

				zShortMul(&poly->poly_a,offset,&t1);
				if (parity)
					zSub(&t1,&t2,&t3);
				else
					zAdd(&t1,&t2,&t3);

				zShortMul(&t3,offset,&t1);
				zAdd(&t1,&poly->poly_c,&Q);
				if (Q.size < 0)
					neg = 1;
				else
					neg = 0;
				Q.size = abs(Q.size);

				smpqs_trial_divide_Q(&Q,fb,full,partial,sieve,offset,
					thisloc,neg,fullfb,cutoff,small_cutoff,start_prime,
					numpoly,parity,closnuf,sieve[thisloc]);
			}
		}
	}

#if defined(SIMD_SIEVE_SCAN)
	// make it safe to perform floating point
	SCAN_CLEAN;
#endif

	zFree(&Q);
	zFree(&t1);
	zFree(&t2);
	zFree(&t3);
	return full->num_r;
}

#define DIVIDE_ONE_PRIME \
	do \
	{						\
		fboffset[++smooth_num] = i;	\
		zShortDiv32(&Q32,prime,&Q32);			\
	} while (zShortMod32(&Q32,prime) == 0);

static void smpqs_trial_divide_Q(z *Q, smpqs_sieve_fb *fb, mpqs_rlist *full, mpqs_rlist *partial,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, fb_list *fullfb, uint32 cutoff,
						  uint8 small_cutoff, uint32 start_prime, int numpoly, uint32 parity,uint8 closnuf,uint8 bits)
{
	smpqs_sieve_fb *fbptr;
	uint32 i,num_f,num_p;
	uint32 root1,root2,prime;
	uint32 r;
	int smooth_num;
	uint16 fboffset[MAX_SMOOTH_PRIMES];
	uint8 logp;
	z32 Q32;

	num_f = full->num_r;
	num_p = partial->num_r;
	
	//we have two signs to worry about.  the sign of the offset tells us how to calculate ax + b, while
	//the sign of Q(x) tells us how to factor Q(x) (with or without a factor of -1)
	//the square root phase will need to know both.  fboffset holds the sign of Q(x).  the sign of the 
	//offset is stored standalone in the relation structure.
	if (sign)
		fboffset[0] = 1;
	else
		fboffset[0] = 0;

	smooth_num=0;
	bits = (255-sieve[j]) + closnuf + 1;

	//take care of powers of two
	while (!(Q->val[0] & 1))
	{
		uint64 mask1;

		//right shift Q
		mask1 = Q->val[1] & 0x1;
		Q->val[1] >>= 1;
		Q->val[0] = (Q->val[0] >> 1) | (mask1 << 63);
		if (Q->val[1] == 0)
			Q->size = 1;

		fboffset[++smooth_num] = 1;
	}

	zInit32(&Q32);
#if BITS_PER_DIGIT == 32
		for (i=0; i<abs(Q->size); i++)
			Q32.val[i] = Q->val[i];
		Q32.size = Q->size;
		Q32.type = Q->type;

#else
		z64_to_z32(Q,&Q32);
#endif

	//completely unrolled trial division by primes that we have not 
	//sieved via the small prime variation
	{
		uint64 q64;
		uint32 tmp1;

		fbptr = fb + 2;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[2];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[2];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		//if ((tmp1 == root1 || tmp1 == root2) || 
		//	(root1 == prime && tmp1 == 0) || (root2 == prime && tmp1 == 0))
		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 2;
				zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 3;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[3];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[3];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 3;
				zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 4;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[4];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[4];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 4;
				zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 5;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[5];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[5];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 5;
				zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 6;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[6];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[6];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 6;
				zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (zShortMod32(&Q32,prime) == 0);
		}
	}

	if (bits < (closnuf + small_cutoff))
		return;

	i=start_prime;
	while (i < fullfb->B)
	{
		uint32 tmp;
		uint64 q64;

		fbptr = fb + i;
		root1 = (fbptr->roots >> 16) + BLOCKSIZE - j;	
		root2 = (fbptr->roots & 0xffff) + BLOCKSIZE - j;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		if (prime > 256)
			break;

		if (Q->size == 1)
		{
			if (Q->val[0] < prime)
				break;
		}

		if (root2 >= prime)
		{
			tmp = root2 + fullfb->list->correction[i];
			q64 = (uint64)tmp * (uint64)fullfb->list->small_inv[i];
			tmp = q64 >> 32; 
			//at this point tmp1 is offset / prime
			tmp = root2 - tmp * prime;

			if (tmp == 0)
			{
				DIVIDE_ONE_PRIME;			
			}
			else if (root1 >= prime)
			{
				tmp = root1 + fullfb->list->correction[i];
				q64 = (uint64)tmp * (uint64)fullfb->list->small_inv[i];
				tmp = q64 >> 32; 
				//at this point tmp1 is offset / prime
				tmp = root1 - tmp * prime;

				if (tmp == 0)
					DIVIDE_ONE_PRIME;			
			}

		}

		i++;
	}

	while (i < fullfb->B)
	{
		uint32 tmp;
		uint64 q64;

		fbptr = fb + i;
		root1 = (fbptr->roots >> 16) + BLOCKSIZE - j;		
		root2 = (fbptr->roots & 0xffff) + BLOCKSIZE - j;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		if (Q->size == 1)
		{
			if (Q->val[0] < prime)
				break;
		}

		if (root2 >= prime)
		{
			tmp = root2 + fullfb->list->correction[i];
			q64 = (uint64)tmp * (uint64)fullfb->list->small_inv[i];
			tmp = q64 >> 40; 
			//at this point tmp1 is offset / prime
			tmp = root2 - tmp * prime;

			if (tmp == 0)
			{
				DIVIDE_ONE_PRIME;			
			}
			else if (root1 >= prime)
			{
				tmp = root1 + fullfb->list->correction[i];
				q64 = (uint64)tmp * (uint64)fullfb->list->small_inv[i];
				tmp = q64 >> 40; 
				//at this point tmp1 is offset / prime
				tmp = root1 - tmp * prime;

				if (tmp == 0)
				{
					DIVIDE_ONE_PRIME;			
				}
			}

		}
		i++;
	}

	//check if it completely factored by looking at the unfactored portion in tmp
	if ((Q32.size == 1) && (Q32.val[0] == 1))
	{
		if (full->num_r == full->allocated) 
		{
			//printf("\nreallocating fulls\n");
			r = full->allocated;
			full->allocated *= 2;
			full->list = (mpqs_r **)realloc(full->list, 
					full->allocated * sizeof(mpqs_r *));
		}
		smpqs_save_relation(full,offset,1,smooth_num+1,num_f,fboffset,numpoly,parity);
	}
	else if ((Q32.size == 1)  && (Q32.val[0] < cutoff))
	{
		smpqs_save_relation(partial,offset,Q32.val[0],smooth_num+1,num_p,fboffset,numpoly,parity);

		if (partial->num_r == partial->allocated) 
		{
			//printf("\nreallocating partials\n");
			r = partial->allocated;
			partial->allocated *= 2;
			partial->list = (mpqs_r **)realloc(partial->list, 
					partial->allocated * sizeof(mpqs_r *));
		}
	}

	zFree32(&Q32);
	return;
}

static void smpqs_save_relation(mpqs_rlist *list, uint32 offset, uint32 largeprime, uint32 num_factors, 
						  uint32 rnum, uint16 *fboffset, int numpoly, uint32 parity)
{
	uint32 i;
	list->list[rnum] = (mpqs_r *)malloc(sizeof(mpqs_r));
	list->list[rnum]->fboffset = (uint16 *)malloc(num_factors*sizeof(uint16));
	for (i=0;i<num_factors;i++)
		list->list[rnum]->fboffset[i] = fboffset[i];
	
	list->list[rnum]->offset = offset;
	list->list[rnum]->largeprime = largeprime;
	list->list[rnum]->parity = parity;
	list->list[rnum]->num_factors = (uint8)(num_factors);
	list->list[rnum]->polynum = numpoly;
	list->num_r++;
	return;
}

static void smpqs_sieve_block(uint8 *sieve, smpqs_sieve_fb *fb, uint32 start_prime, 
	uint8 s_init, fb_list *fullfb)
{
	uint32 prime, root1, root2, tmp, stop;
	uint32 B=fullfb->B;
	uint32 i;
	uint8 logp, *s2;
	smpqs_sieve_fb *fbptr;

	//initialize block
	memset(sieve,s_init,BLOCKSIZE);
	
	//we've now filled the entire block with the small fb primes, proceed with the rest
	for (i=start_prime;i<B;i++)
	{	
		fbptr = fb + i;
		prime = fbptr->prime_and_logp >> 16;
		root1 = fbptr->roots >> 16;
		root2 = fbptr->roots & 0xffff;
		logp = fbptr->prime_and_logp & 0xff;

		stop = BLOCKSIZE - prime;
		s2 = sieve + prime;

		while (root2 < stop)
		{
			sieve[root1] -= logp;
			sieve[root2] -= logp;
			s2[root1] -= logp;
			s2[root2] -= logp;
			root1 += (prime << 1);
			root2 += (prime << 1);
		}

		while (root2 < BLOCKSIZE)
		{
			sieve[root1] -= logp;
			sieve[root2] -= logp;
			root1 += prime;
			root2 += prime;
		}

		//don't forget the last proot1[i], and compute the roots for the next block
		if (root1 < BLOCKSIZE)
		{
			sieve[root1] -= logp;
			root1 += prime;
			//root1 will be bigger on the next iteration, switch them now
			tmp = root2;
			root2 = root1;
			root1 = tmp;
		}
			
		fbptr->roots = ((root1 - BLOCKSIZE) << 16) | (root2 - BLOCKSIZE);
	}

	return;
}

void smpqs_nextD(uint32 *polyd, uint32 *polyd_id, z *n)
{
	uint32 r;
	do 
	{
		(*polyd_id)++;
		*polyd = PRIMES[*polyd_id];
		r = zShortMod(n,*polyd);
	} while ((jacobi_1(r,*polyd) != 1) || ((*polyd & 3) != 3));

	return;
}

void smpqs_computeB(mpqs_poly *poly, uint32 polyd, z *n)
{
	//using poly_d, compute poly_b and poly_a = poly_d^2
	//int i;
	z t1, t3, t4, t5, t6;
	fp_digit ut1, pa, pb;
	fp_digit ult1;

	zInit(&t1);
	zInit(&t3);
	zInit(&t4);
	zInit(&t5);
	zInit(&t6);
	//poly_a = d^2.  we just found d.  also compute b using Hegel theorem and lifting.

	//t0 = n^(d-3)/4 mod d
	zModExp_1(n,(polyd-3) >> 2,polyd,&ut1);	

	//t1 = h1 = n*t0 mod d
	zShortMul(n,ut1,&t3);		
	ut1 = zShortMod(&t3,polyd);
	sp2z(ut1,&t1);

	//t3 = n - h1^2
	zShortSub(n,ut1 * ut1,&t3);		

	//t4 = (n - h1^2)/d
	zShortDiv(&t3,polyd,&t4);	

	//t4 = (n - h1^2)/d mod d
	ult1 = zShortMod(&t4,polyd);	
	sp2z(ult1,&t4);

	//t5 = 2*h1;
	//zShiftLeft(&t5,&t1,1);		
	sp2z(ut1 * 2, &t5);

	//compute t6 = (2*h1)^-1 mod d = (2*h1)^(d-2) mod d
	zModExp_1(&t5,polyd-2,polyd,&ut1);
	sp2z(ut1,&t6);

	//compute t3 = h2 = ((2*h1)^-1 * (n - h1^2)/d) mod d
	zMul(&t6,&t4,&t5);		
	ut1 = zShortMod(&t5,polyd);
	sp2z(ut1,&t3);

	//compute t5 = h1 + h2*D
	zShortMul(&t3,polyd,&t4);	
	zAdd(&t4,&t1,&t5);

	//we're now done with d, so compute a = d^2
	pa = (uint64)polyd * (uint64)polyd;
	sp2z(pa,&poly->poly_a);

	//compute b = h1 + h2*d mod a
	pb = zShortMod(&t5,pa);	

	//make sure b < a/2
	if (pb > (pa >> 1))
		sp2z(pa - pb,&poly->poly_b);
	else
		sp2z(pb,&poly->poly_b);

	//now that we have b, compute c = (b*b - n)/a
	zSqr(&poly->poly_b,&t1);
	zSub(&t1,n,&t3);
	zShortDiv(&t3,pa,&poly->poly_c);

	zFree(&t1);
	zFree(&t3);
	zFree(&t4);
	zFree(&t5);
	zFree(&t6);
	return;
}

void smpqs_computeRoots(mpqs_poly *poly, fb_list *fb, uint32 *modsqrt, 
	smpqs_sieve_fb *fbp, smpqs_sieve_fb *fbn, uint32 start_prime)
{
	//the roots are computed using a and b as follows:
	//(+/-t - b)(a)^-1 mod p
	//assume b > t
	//uint32 root1, root2, prime, amodp;
	int root1, root2, prime, amodp, bmodp, x;
	uint32 i;
	
	for (i=2;i<7;i++)
	{
		uint64 q64, tmp, t2;

		prime = fb->list->prime[i]; 
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		bmodp = zShortMod(&poly->poly_b,prime);
		x = (int)root1 - bmodp;
		if (x < 0) x += prime; root1 = x;
		x = (int)root2 - bmodp;
		if (x < 0) x += prime; root2 = x;

		amodp = zShortMod(&poly->poly_a,prime);
		amodp = modinv_1(amodp,prime);

		//root1 = (uint32)((uint64)amodp * (uint64)root1 % (uint64)prime);
		t2 = (uint64)amodp * (uint64)root1;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root1 = t2 - tmp * prime;

		//root2 = (uint32)((uint64)amodp * (uint64)root2 % (uint64)prime);	
		t2 = (uint64)amodp * (uint64)root2;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root2 = t2 - tmp * prime;

		// we don't sieve these primes, so ordering doesn't matter
		fbp[i].roots = (uint32)((root1 << 16) | root2);
		if (root2 == 0) fbn[i].roots = 0;
		else fbn[i].roots = (uint32)((prime - root2) << 16);
		if (root1 == 0) fbn[i].roots |= 0;
		else fbn[i].roots |= (uint32)(prime - root1);
	}

	for (i=7;i<fb->B;i++)
	{
		//take the fast method of computing the inverse from lenstra...
		//root1 = fb->list[i].c1;
		//root2 = fb->list[i].c2;
		//prime = fb->list[i].prime;
		uint64 q64, tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		if (prime > 256)
			break;

		bmodp = zShortMod(&poly->poly_b,prime);
		x = (int)root1 - bmodp;
		if (x < 0) x += prime;
		root1 = x;
		x = (int)root2 - bmodp;
		if (x < 0) x += prime;
		root2 = x;
	
		//now (t - b) mod p is in root1 and (-t - b) mod p is in root2
		//find a^-1 mod p = inv(a mod p) mod p
		amodp = zShortMod(&poly->poly_a,prime);
		amodp = modinv_1(amodp,prime);

		//root1 = (uint32)((uint64)amodp * (uint64)root1 % (uint64)prime);
		t2 = (uint64)amodp * (uint64)root1;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root1 = t2 - tmp * prime;

		//root2 = (uint32)((uint64)amodp * (uint64)root2 % (uint64)prime);	
		t2 = (uint64)amodp * (uint64)root2;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root2 = t2 - tmp * prime;

		if (root2 < root1)
		{
			fbp[i].roots = (uint32)((root2 << 16) | root1);
			fbn[i].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
		}
		else
		{
			fbp[i].roots = (uint32)((root1 << 16) | root2);
			fbn[i].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
		}
	}

	for ( ;i<fb->B;i++)
	{
		//take the fast method of computing the inverse from lenstra...
		//root1 = fb->list[i].c1;
		//root2 = fb->list[i].c2;
		//prime = fb->list[i].prime;
		uint64 q64, tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		bmodp = zShortMod(&poly->poly_b,prime);
		x = (int)root1 - bmodp;
		if (x < 0) x += prime;
		root1 = x;
		x = (int)root2 - bmodp;
		if (x < 0) x += prime;
		root2 = x;
	
		//now (t - b) mod p is in root1 and (-t - b) mod p is in root2
		//find a^-1 mod p = inv(a mod p) mod p
		amodp = zShortMod(&poly->poly_a,prime);
		amodp = modinv_1(amodp,prime);

		//root1 = (uint32)((uint64)amodp * (uint64)root1 % (uint64)prime);
		t2 = (uint64)amodp * (uint64)root1;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 40; 
		root1 = t2 - tmp * prime;

		//root2 = (uint32)((uint64)amodp * (uint64)root2 % (uint64)prime);	
		t2 = (uint64)amodp * (uint64)root2;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 40; 
		root2 = t2 - tmp * prime;

		if (root2 < root1)
		{
			fbp[i].roots = (uint32)((root2 << 16) | root1);
			fbn[i].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
		}
		else
		{
			fbp[i].roots = (uint32)((root1 << 16) | root2);
			fbn[i].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
		}
	}
	return;
}

int qcomp_smpqs(const void *x, const void *y)
{
	mpqs_r **xx = (mpqs_r **)x;
	mpqs_r **yy = (mpqs_r **)y;
	
	if (xx[0]->largeprime > yy[0]->largeprime)
		return 1;
	else if (xx[0]->largeprime == yy[0]->largeprime)
		return 0;
	else
		return -1;
}

void fastDiv(z *n, z *d, z *q, z *r)
{
	//compute q = n / d and r = n % d
	//given that we know that n is at most 3 digits
	//and d is at most 2 digits

	fp_digit nn[2], d1, d2, q1, q2, r1, r2;

	



	return;
}

static uint64 smpqs_bitValRead64(uint64 **m, int row, int col);

int smpqs_BlockGauss(mpqs_rlist *full, mpqs_rlist *partial, z *apoly, z *bpoly,
			fb_list *fb, z *n, int mul, 
			z *factors,uint32 *num_factor)
{
	int i,j,k,l,a,q,polynum;
	int *bl;
	uint8 **m;		//matrix of the powers of the prime decompositions of the relations over the factor base
	uint64 **m2_64;	//m mod 2, packed into 32 bit words
	uint64 **aug_64;	//matrix to store the permutations of the rows of m2, packed into 32 bit words
	uint32 largep, bool_val, B = fb->B;
	uint32 *partial_index;
	int num_f,num_p;
	int num_r,num_col,num_col_aug,set_continue;
	const int blocksz = 64;

	uint32 *pd;
	uint32 r;
	z zx, zy, tmp, tmp2, tmp3, tmp4, nn,tmp_a,input;
	
	zInit(&zx);
	zInit(&zy);
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);
	zInit(&tmp4);
	zInit(&nn);
	zInit(&input);
	zInit(&tmp_a);

	num_f = full->num_r;
	num_p = partial->act_r;

	num_r = full->num_r + partial->act_r;
	num_col = (uint32)((B/blocksz)+1);
	num_col_aug = (uint32)(num_r/blocksz+1);

	//allocate storage based on total number of relations.
	pd = (uint32 *)malloc(B * sizeof(uint32));
	partial_index = (uint32 *)malloc(num_p * sizeof(uint32));

	aug_64 = (uint64 **)malloc(num_r * sizeof(uint64 *));
	for (i=0; i<num_r; i++)
		aug_64[i] = (uint64 *)malloc(num_col_aug * sizeof(uint64));

	m2_64 = (uint64 **)malloc(num_r * sizeof(uint64 *));
	for (i=0; i<num_r; i++)
		m2_64[i] = (uint64 *)malloc(num_col * sizeof(uint64));

	m = (uint8 **)malloc(num_r * sizeof(uint8 *));
	for (i=0; i<num_r; i++)
		m[i] = (uint8 *)malloc(B * sizeof(uint8));

	bl = (int *)malloc(num_r * sizeof(int));

	//write fulls to m
	for (i=0;i<num_f;i++)
	{
		//Initialize
		for (j=0;j<(int)B;j++)
			m[i][j] = 0;

		//copy the pd's of the fboffsets to the correct location in m
		//offset 0 is special - indicates the parity of the offset
		m[i][0] = (uint8)full->list[i]->fboffset[0];
		j=1;
		while (j<full->list[i]->num_factors)
		{
			m[i][full->list[i]->fboffset[j]]++;
			j++;
		}
	}

	//write fulls from partials to m, probably also redundant to do it this way?
	largep = partial->list[0]->largeprime;
	j=num_f;
	for (i=1;i<(int)partial->num_r;i++)
	{
		if (partial->list[i]->largeprime == largep)
		{
			//this partial's largep is the same as the one before, add the pd's and copy to m
			for (k=0;k<(int)B;k++)
				m[j][k] = 0;

			//do the factor of -1
			m[j][0] = (uint8)partial->list[i-1]->fboffset[0];
			//then the rest
			k=1;
			while (k<partial->list[i-1]->num_factors) 
			{
				m[j][partial->list[i-1]->fboffset[k]]++;
				k++;
			}
			//factor of -1
			m[j][0] += partial->list[i]->fboffset[0];
			//the rest
			k=1;
			while (k<partial->list[i]->num_factors)
			{
				m[j][partial->list[i]->fboffset[k]]++;
				k++;
			}
			//remember the index of the partial that made this full relation, we'll need it later
			partial_index[j-num_f]=i;
			//increment the relation counter
			j++;
		}
		largep = partial->list[i]->largeprime;
	}

	//construct the bit matrix
	for (i=0;i<num_r;i++)
	{
		for (j=0;j<num_col;j++)
		{
			m2_64[i][j] = 0;
			for (k=0;k<blocksz;k++)
			{
				if ((blocksz*j+k) < (int)B)
					m2_64[i][j] |= ((uint64)((uint64)m[i][blocksz*j+k]%2) << k);
			}
		}
	}

	//construct augmented matrix
	for (i=0;i<num_r;i++)
	{
		for (j=0;j<num_col_aug;j++)
		{
			aug_64[i][j] = 0;
			for (k=0;k<blocksz;k++)
			{
				if ((blocksz*j+k)==i)
					aug_64[i][j] = ((uint64)(1) << (uint64)(k));
			}
		}
	}

	*num_factor=0;

	// remove the multiplier from the input
	zShortDiv(n,mul,&input);

	//initialize blacklist
	for (i=0;i<num_r;i++) bl[i] = 0;
	//search over all columns, right to left (more sparse on the right side)
	for (i=B-1;i>=0;i--)
	{
		//and all rows
		for (j=0;j<num_r;j++)
		{
			//if the j'th row, i'th bit is 1 and not blacklisted, continue
			bool_val = (smpqs_bitValRead64(m2_64,j,i) != 0) && (bl[j] == 0);
			if (bool_val)
			{
				//add the j'th row mod 2 to all rows after it with a 1 in the ith column
				for (k=j+1;k<num_r;k++)
				{
					bool_val = (smpqs_bitValRead64(m2_64,k,i) != 0) && (bl[k] == 0);
					if (bool_val)
					{
						//found one in the k'th row.  add to the j'th row starting at column i.
						//record the addition in the augmented matrix
						for (l=(uint32)(i/blocksz);l>=0;l--)
							m2_64[k][l] = m2_64[j][l] ^ m2_64[k][l];
						for (l=0;l<num_col_aug;l++)
							aug_64[k][l] = aug_64[k][l] ^ aug_64[j][l];
						
						//then check if the row is all zeros
						a=0;
						for (l=(uint32)(i/blocksz);l>=0;l--)
							a = a || m2_64[k][l];

						if (a==0)
						{
							//initialize solution vector
							for (l=0;l<(int)B;l++) pd[l] = 0;

							//found a potential solution. check it.
							for (l=0;l<num_r;l++)
							{
								bool_val = smpqs_bitValRead64(aug_64,k,l) != 0;
								if (bool_val)
								{
									//then the l'th row of m was involved
									for (q=0;q<(int)B;q++)
										pd[q] += m[l][q];
								}
							}

							//compute x mod n
							zCopy(&zOne,&zy);
							zCopy(&zOne,&zx);
							for (l=0;l<num_r;l++)
							{
								bool_val = smpqs_bitValRead64(aug_64,k,l) != 0;
								if (bool_val)
								{
									//printf("accumulating relation %d\n",l);
									//then the l'th relation is involved in the product of relations
									if (l >= num_f)
									{
										uint64 pa,pb,d1,d2;
										//if l >= num_f, then this row refers to a relation generated from two partials.
										//we'll need to go back to the two partial locations to find the two offsets to
										//multiply together
										//luckily, we've remembered the index in the complete list of partials that 
										//created this full relation
											
										//our relation is of the form (ax + b)^2 == a(ax^2 + 2bx + c) mod n
										//(ax^2 + 2bx + c) is what we trial divided, and we remembered
										//a and b, so we can form the left hand side easily
											
										// poly_a and b will fit in a uint64 in smallmpqs

										//recreate poly_b from poly_a (store instead??)
										polynum = partial->list[partial_index[l-num_f]]->polynum;
										pb = bpoly[polynum].val[0];
											
										//compute Q1(x)
										pa = apoly[polynum].val[0];
										d1 = sqrt(pa);
										
										zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]]->offset,&tmp);
										if (partial->list[partial_index[l-num_f]]->parity)
											zShortSub(&tmp,pb,&tmp);
										else
											zShortAdd(&tmp,pb,&tmp);

										//include 'a'
										zShortMul(&zy,d1,&tmp2);
										zDiv(&tmp2,n,&tmp3,&zy);

										//compute Q2(x)
										polynum = partial->list[partial_index[l-num_f]-1]->polynum;
										pb = bpoly[polynum].val[0];
											
										//compute Q1(x)
										pa = apoly[polynum].val[0];
										d2 = sqrt(pa);
										
										zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]-1]->offset,&tmp3);
										if (partial->list[partial_index[l-num_f]-1]->parity)
											zShortSub(&tmp3,pb,&tmp2);
										else
											zShortAdd(&tmp3,pb,&tmp2);

										//compute Q(x1)*Q(x2)
										zMul(&tmp,&tmp2,&tmp4);	
										zMul(&zx,&tmp4,&tmp);		//accumulate with previous terms
										zDiv(&tmp,n,&tmp2,&zx);		//mod n

										//include the large prime in mp_y
										zShortMul(&zy,partial->list[partial_index[l-num_f]]->largeprime,&tmp2);
										zDiv(&tmp2,n,&tmp3,&zy);

										//include 'a'
										zShortMul(&zy,d2,&tmp2);	
										zDiv(&tmp2,n,&tmp3,&zy);
									}
									else
									{
										uint64 pa,pb,d1;

										//recreate poly_b from poly_a (store instead??)
										polynum = full->list[l]->polynum;
										pb = bpoly[polynum].val[0];
											
										//compute Q1(x)
										zShortMul(&apoly[polynum],full->list[l]->offset,&tmp);
										pa = apoly[polynum].val[0];
										d1 = sqrt(pa);

										if (full->list[l]->parity)
											zShortSub(&tmp,pb,&nn);
										else
											zShortAdd(&tmp,pb,&nn);

										zMul(&zx,&nn,&tmp);			//accumulate with previous terms
										zDiv(&tmp,n,&tmp2,&zx);		//mod n

										zShortMul(&zy,d1,&tmp2);	//sqrt(a) = d is part of mp_y
										zDiv(&tmp2,n,&tmp3,&zy);
									}
								}
							}

							//compute y mod n
							//ignore the factor of -1 in this operation
							for (l=1;l<(int)B;l++)
							{
								if (pd[l] > 0)
								{
									sp2z(fb->list->prime[l],&tmp);
									//pd tracks the exponents of the smooth factors.  we know they are all even
									//at this point.  we don't want to compute pd^2, so divide by 2.
									//computing the explicit exponentiation and then reducing is
									//slightly faster than doing modexp in smallmpqs.
									zExp(pd[l]/2,&tmp,&tmp2);
									zDiv(&tmp2,n,&tmp4,&tmp3);
									zMul(&tmp3,&zy,&tmp4);
									zDiv(&tmp4,n,&tmp2,&zy);
								}
							}

							//split this off into a subroutine... also look for all non-trivial factors if one is composite
							//compute gcd(x-y,n)
							zSub(&zx,&zy,&tmp);
							zLEGCD(&tmp,n,&nn);
							
							if ((r = (uint32)zShortDiv(&nn,mul,&tmp)) == 0)
							{
								//mul divides this factor
								zCopy(&tmp,&nn);
							}

							if (!(nn.val[0] & 0x1))
							{
								//if it's not odd (factors of 2 creep in there, for some reason
								//remove the 2's
								while (!(nn.val[0] & 0x1))
									zShiftRight_1(&nn,&nn);
							}							

							if ((zCompare(&nn,&zOne) > 0) && (zCompare(&nn,&input) < 0))
							{

								//check that we havent' already found this one
								set_continue = 0;
								for (l=0;l<(int)*num_factor;l++)
								{
									if (zCompare(&nn,&factors[l]) == 0)
										set_continue = 1;
								}
								if (set_continue)
									continue;

								zCopy(&nn,&factors[*num_factor]);

								(*num_factor)++;
								if (*num_factor > MAX_FACTORS)
								{
									printf("max number of factors found in block gauss\n");
									goto free;
								}
								//check if we're done by accumulating all factors and comparing to n
								zCopy(&factors[0],&nn);
								for (l=1;l<(int)*num_factor;l++)
								{
									zCopy(&factors[l],&tmp);
									zMul(&tmp,&nn,&tmp2);
									zCopy(&tmp2,&nn);
								}
								if (zBits(&nn) + 10 >= zBits(n))
								{
									//+ 10 accounts for the multiplier in n
									//found all factors, done
									goto free;
								}

								//check the other factor
								zCopy(&input,&tmp);
								zDiv(&tmp,&nn,&tmp2,&tmp3);

								zCopy(&tmp2,&tmp);
	
								//check that we havent' already found this one
								set_continue = 0;
								for (l=0;l<(int)*num_factor;l++)
								{
									if (zCompare(&tmp,&factors[l]) == 0)
										set_continue = 1;
								}
								if (set_continue)
									continue;

								zCopy(&tmp,&factors[*num_factor]);

								(*num_factor)++;
								if (*num_factor > MAX_FACTORS)
								{
									printf("max number of factors found in block gauss\n");
									goto free;
								}
								//check if we're done by accumulating all factors and comparing to n
								zCopy(&factors[0],&nn);
								for (l=1;l<(int)*num_factor;l++)
								{
									zCopy(&factors[l],&tmp);
									zMul(&tmp,&nn,&tmp2);
									zCopy(&tmp2,&nn);
								}
								if (zBits(&nn) + 10 >= zBits(n))
								{
									//+ 10 accounts for the multiplier in n
									//found all factors, done
									goto free;
								}
							} //if non-trivial factor
						} //if a == 0
					} //if found in k'th row
				} //add jth row mod 2 to all appropriate rows after it
				//blacklist the j'th row
				bl[j] = 1;
			} //if not blacklisted
		} //for all rows
	} //for all columns

	printf("matrix exhausted\n");
	r = (uint32)zShortDiv(n,mul,&tmp);
	for (i=0;(uint32)i<*num_factor;i++)
	{
		zCopy(&tmp,&nn);
		zDiv(&nn,&factors[i],&tmp,&tmp2);
	}

free:
	free(pd);
	free(partial_index);
	for (i=0; i<num_r; i++)
		free(aug_64[i]);
	free(aug_64);
	for (i=0; i<num_r; i++)
		free(m2_64[i]);
	free(m2_64);
	for (i=0; i<num_r; i++)
		free(m[i]);
	free(m);

	free(bl);
	zFree(&zx);
	zFree(&zy);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&tmp4);
	zFree(&nn);
	zFree(&tmp_a);
	return 0;
}

static uint64 smpqs_masks64[64] = {0x1,0x2,0x4,0x8,
							0x10,0x20,0x40,0x80,
							0x100,0x200,0x400,0x800,
							0x1000,0x2000,0x4000,0x8000,
							0x10000,0x20000,0x40000,0x80000,
							0x100000,0x200000,0x400000,0x800000,
							0x1000000,0x2000000,0x4000000,0x8000000,
							0x10000000,0x20000000,0x40000000,0x80000000ULL,
							0x100000000ULL,0x200000000ULL,0x400000000ULL,0x800000000ULL,
							0x1000000000ULL,0x2000000000ULL,0x4000000000ULL,0x8000000000ULL,
							0x10000000000ULL,0x20000000000ULL,0x40000000000ULL,0x80000000000ULL,
							0x100000000000ULL,0x200000000000ULL,0x400000000000ULL,0x800000000000ULL,
							0x1000000000000ULL,0x2000000000000ULL,0x4000000000000ULL,0x8000000000000ULL,
							0x10000000000000ULL,0x20000000000000ULL,0x40000000000000ULL,0x80000000000000ULL,
							0x100000000000000ULL,0x200000000000000ULL,0x400000000000000ULL,0x800000000000000ULL,
							0x1000000000000000ULL,0x2000000000000000ULL,0x4000000000000000ULL,0x8000000000000000ULL};

static uint64 smpqs_bitValRead64(uint64 **m, int row, int col)
{
	//col is the column in 0 to B-1 representation
	//read the bit in the packed 64 bit representation of the appropriate row
	//don't bother to check the bounds of m w.r.t row and col, assume caller knows what it's doing
	//return 0 if bit not set, 1 << bit offset otherwize
	int offset, mcol;
	mcol = col >> 6;
	offset = col & 63;
	return (m[row][mcol] & smpqs_masks64[offset]);
}





