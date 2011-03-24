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

static void smpqs_sieve_block(uint8 *sieve, sieve_fb *fb, uint32 start_prime, uint8 s_init, fb_list_mpqs *fullfb);
static void smpqs_test_block(uint8 *sieve, fb_element_mpqs *fb, uint32 start_prime);
static int smpqs_check_relations(uint32 sieve_interval, uint32 blocknum, uint8 sieve[] ,z *n, 
								mpqs_poly *poly, uint8 closnuf,
								sieve_fb *fb,fb_list_mpqs *fullfb,mpqs_rlist *full, mpqs_rlist *partial,uint32 cutoff,
								uint8 small_cutoff, uint32 start_prime, uint32 parity, uint32 *num, int numpoly);
static void smpqs_trial_divide_Q(z *Q,sieve_fb *fb,mpqs_rlist *full,mpqs_rlist *partial,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, fb_list_mpqs *fullfb, uint32 cutoff,
						  uint8 small_cutoff, uint32 start_prime, int numpoly, uint32 parity,
						  uint8 closnuf, uint8 bits);
static void smpqs_save_relation(mpqs_rlist *list, uint32 offset, uint32 largeprime, uint32 num_factors, 
						  uint32 rnum, uint16 *fboffset, int numpoly, uint32 parity);
void smpqs_make_fb_mpqs(fb_list_mpqs *fb, z *n);
void smpqs_computeB(mpqs_poly *poly, z *n);
void smpqs_nextD(z *poly_d, z *n);
void smpqs_computeRoots(mpqs_poly *poly, fb_list_mpqs *fb, uint32 start_prime);
int qcomp_smpqs(const void *x, const void *y);

void smpqs_make_fb_mpqs(fb_list_mpqs *fb, z *n)
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
		//printf("%s mod %u = %u\n",z2decstr(n,&gstr1),spSOEprimes[i],zShortMod(n,spSOEprimes[i]));
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
			fb->list[j].prime = prime;

			if (root2 > root1)
			{
				fb->list[j].c1 = root1;
				fb->list[j].c2 = root2;
			}
			else
			{
				fb->list[j].c1 = root2;
				fb->list[j].c2 = root1;
			}
			fb->list[j].logprime = logp;
			fb->list[j].inv = modinv_1(prime,1);
			
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
			fb->list[j].prime = prime;
			if (root2 > root1)
			{
				fb->list[j].c1 = root1;
				fb->list[j].c2 = root2;
			}
			else
			{
				fb->list[j].c1 = root2;
				fb->list[j].c2 = root1;
			}
			fb->list[j].logprime = logp;
			fb->list[j].inv = modinv_1(prime,1);
			
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
	fb_list_mpqs *fb;
	sieve_fb *fb_sieve_p,*fb_sieve_n;
	mpqs_poly *poly;
	
	z *factors;
	z tmp, tmp2, tmp3, sqrt_n;
	z *apoly, *bpoly;

	uint32 numpoly, polyalloc;

	double ln_n,sum,avg,sd;
	double t_time;
	struct timeval tstart, tend;
	TIME_DIFF *	difference;
	uint32 mul,i,j,j2;
	uint32 pmax;							//largest prime in factor base
	uint32 cutoff, last_numfull, last_numpartial, num_needed, num_expected;
	uint32 sieve_interval;
	uint32 start_prime, offset1, min_block;
	uint32 check_total, check_inc, num, max_f;
	int digits_n,charcount;
	uint32 num_factors;
	uint8 *sieve;							//sieve values
	uint8 s_init, blockinit;				//initial sieve value
	uint8 closnuf, small_bits, max_bits, blockbits, smallprime = 0;
	FILE *sieve_log;

	if (n->val[0] == 1 && n->size == 1)
		return;

	if ((n->val[0] % 2) == 0)
	{
		printf("n must be odd\n");
		return;
	}

	if (isPrime(n))
	{
		n->type = PRP;
		//add_to_factor_list(fobj, n);
		zCopy(&zOne,n);
		return;
	}

	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);

	//default mpqs parameters
	sieve_params.fudge_factor = 1.3;
	sieve_params.large_mult = 30;
	sieve_params.num_blocks = 40;
	sieve_params.num_extra_relations = 64;
	sieve_params.small_limit = 100;

	zNroot(n,&tmp,2);
	zSqr(&tmp,&tmp2);
	if (zCompare(&tmp2,n) == 0)
	{
		if (isPrime(&tmp))
		{
			tmp.type = PRP;
			//add_to_factor_list(fobj, &tmp);
			//add_to_factor_list(fobj, &tmp);
		}
		else
		{
			tmp.type = COMPOSITE;
			//add_to_factor_list(fobj, &tmp);
			//add_to_factor_list(fobj, &tmp);
		}
		zCopy(&zOne,n);
		zFree(&tmp);
		zFree(&tmp2);
		zFree(&tmp3);
		return;
	}

	//allocate the space for the factor base
	fb = (fb_list_mpqs *)malloc(sizeof(fb_list_mpqs));

	//calculate approximate factor base bound.  Should do some experiments to determine if this really is optimal.
	//I divide by two because only approximately half of the primes up to the bound will be quadratic residues
	//mod n, and thus be eligible for the factor base.
	i = zBits(n) - 1;
	ln_n = (double)i * log(2.0);
	//calculate the optimal largest prime in the factor base.  for mpqs, this is too high.
	fb->B = (uint32)(exp(.5 * sqrt(ln_n * log(ln_n))));
	
	//compute the number of digits in n 
	digits_n = ndigits(n);

	//swag at better limit
	if (digits_n < 40)
		fb->B = (uint32)((double)fb->B / 1.5);
	else if (digits_n < 55)
		fb->B = (uint32)((double)fb->B / 2.8);
	else if (digits_n < 68)
		fb->B = (uint32)((double)fb->B / 2.8);
	else
		fb->B = (uint32)((double)fb->B / 3.5);	//2.8

	//determine the number of primes we'll need such that the largest prime is close to optimal
	for (i=0;i<fb->B;i++)
	{
		if (spSOEprimes[i] > fb->B)
			break;
	}
	//because only about half the primes will be suitable
	fb->B = i/2;

	//crude tuning of sieve interval based on digits in n
	sieve_params.num_blocks = 1;
	sieve_params.num_extra_relations = 16;

	//compute how often to check our list of partial relations
	//and update the gui.
	check_inc = fb->B/20;
	check_total = check_inc;

	//allocate storage for relations based on the factor base size
	max_f = fb->B + sieve_params.num_extra_relations;	//num_f should never get this high, should find at least a few fulls from partials.
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
	fb->list = (fb_element_mpqs *)malloc((size_t)(fb->B * sizeof(fb_element_mpqs)));
	fb_sieve_p = (sieve_fb *)malloc((size_t)(fb->B * sizeof(sieve_fb)));
	fb_sieve_n = (sieve_fb *)malloc((size_t)(fb->B * sizeof(sieve_fb)));
	
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
	mul = (uint32)choose_multiplier(n,fb->B);
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

	zNroot(&tmp2,&poly->poly_d,2);
	if (!(poly->poly_d.val[0] & 1))
		zAdd(&poly->poly_d,&zOne,&poly->poly_d);

	//alternating poly_a values to either side of the ideal poly_a makes the average poly_a
	//closer to ideal, and hence more effective at producing relations.  I've verified that the poly_a's produced
	//are about half the size of those going in only one direction.
	//in practice though, i don't see much of a difference in timing.
	//this is because the percent difference to the ideal poly_a is tiny in either case.  
	smpqs_nextD(&poly->poly_d,n);
	smpqs_computeB(poly,n);

	//printf("a: %s\nb: %s\nc: %s\n",z2decstr(&poly->poly_a,&gstr1),z2decstr(&poly->poly_b,&gstr2),z2decstr(&poly->poly_c,&gstr3));

	fb->list[0].prime = 1;
	fb->list[1].prime = 2;

	//construct the factor base, and copy to the sieve factor base
	smpqs_make_fb_mpqs(fb,n);
	for (i=2;i<fb->B;i++)
	{
		fb_sieve_p[i].prime = fb->list[i].prime;
		fb_sieve_p[i].logprime = fb->list[i].logprime;
		fb_sieve_n[i].prime = fb->list[i].prime;
		fb_sieve_n[i].logprime = fb->list[i].logprime;
	}

	//find the root locations of the factor base primes for this poly
	smpqs_computeRoots(poly,fb,2);

	/*
	printf("fb dump\nprimes = \n");
	for (i=0;i<fb->B;i++)
		printf("%u ",fb->list[i].prime);
	printf("\n\n");
	*/

	//compute sieving limits
	fb->small_B = MIN(fb->B,((INNER_BLOCKSIZE)/(sizeof(sieve_fb))));
	for (i=fb->small_B;i<fb->B;i++)
	{
		if (fb->list[i].prime > BLOCKSIZE)
			break;
	}
	fb->med_B = i;

	pmax = fb->list[fb->B-1].prime;
	cutoff = pmax * sieve_params.large_mult;

	//how do I do a multiple precision natural logarithm?
	//convert from a logarithm that is easy to do, i.e. log2
	//log_e(x) = log_2(x) * log_e(2)
	//log_e(x) = (bits(x) - 1) * log_e(2)
	//ln_n = (double)(mpBits(n) - 1) * log(2.0);
	
	//small prime variation.  don't turn this on unless n is sufficiently large, because if it is too
	//small, the benefit gained from not sieving the small primes is generally outweighed by the extra
	//relation checking that occurs.  for now, if the closnuf cutoff is > 90, then use the small prime variation
	//initialize the sieve array to the average value of the logs of all the primes we are not sieving.  not
	//perfect as an estimate to the actual contribution of small primes, but a good working start.
	
	//compute the number of bits in M/2*sqrt(N/2), the approximate value
	//of residues in the sieve interval
	//sieve locations greater than this are worthy of trial dividing
	closnuf = (uint8)(double)((zBits(n) - 1)/2);
	closnuf += (uint8)(log((double)sieve_interval/2)/log(2.0));
	closnuf -= (uint8)(sieve_params.fudge_factor * log(cutoff) / log(2.0));
	
	if (closnuf >= 70)
	{
		smallprime = 1;
		//test the contribution of the small primes to the sieve.  
		for (i=2;i<fb->B;i++)
		{
			if (fb->list[i].prime > sieve_params.small_limit)
				break;
		}

		start_prime = i;
		smpqs_test_block(sieve,fb->list,start_prime);
		sum=0;
		for (i=0;i<BLOCKSIZE;i++)
			sum += sieve[i];
		avg = sum/BLOCKSIZE;
		sum = 0;
		for (i=0;i<BLOCKSIZE;i++)
			sum += ((sieve[i] - avg) * (sieve[i] - avg));
		sd = sum/BLOCKSIZE;
		sd = sqrt(sd);
		small_bits = (uint8)(avg + sd);
		closnuf -= small_bits;
		s_init=0;
	}
	else
	{
		s_init = 0;
		start_prime = 2;
		small_bits = 0;
	}

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
		if (smallprime)
			printf("small prime variation in effect, with initial value %d\n",small_bits);
		printf("trial factoring cutoff at %d bits\n",closnuf);
		printf("==== sieving in progress ====\n");
	}

	numpoly = 0;
	num = 0;
	last_numfull = 0;
	last_numpartial = 0;
	charcount=0;
	while (1)
	{
		//checkpoly(poly,n);
		//copy current poly into the poly lists
		if (numpoly < polyalloc)
		{
			zCopy(&poly->poly_a,&apoly[numpoly]);
			zCopy(&poly->poly_b,&bpoly[numpoly]);
		}
		else
		{
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

		//copy root info into the fb_sieve
		for (i=2;i<fb->B;i++)
		{
			fb_sieve_n[i].root1 = fb->list[i].nroot1;
			fb_sieve_n[i].root2 = fb->list[i].nroot2;
		}
		for (i=2;i<fb->B;i++)
		{
			fb_sieve_p[i].root1 = fb->list[i].proot1;
			fb_sieve_p[i].root2 = fb->list[i].proot2;
		}

		//the minimum value for the current poly_a and poly_b occur at offset (-b + sqrt(N))/a
		//compute minimum bits for Q
		//poly->poly_b.size *= -1;
		zSub(&sqrt_n,&poly->poly_b,&tmp);
		//poly->poly_b.size *= -1;
		zDiv(&tmp,&poly->poly_a,&tmp2,&tmp3);
		min_block = tmp2.val[0]/BLOCKSIZE;
		
		//compute an approximate bit value for Q based on the average offset for the minimum block
		offset1 = (min_block<<BLOCKBITS) + BLOCKSIZE/2;
		zShiftLeft(&tmp2,&poly->poly_b,1);
		zShortMul(&poly->poly_a,offset1,&tmp);
		zAdd(&tmp,&tmp2,&tmp3);
		zShortMul(&tmp3,offset1,&tmp);
		zAdd(&tmp,&poly->poly_c,&tmp2);
		blockbits = zBits(&tmp2);
		
		//sieve each block with the current polynomial.  
		for (j2=0;j2<(sieve_interval/BLOCKSIZE);j2++)
		{	
			//adjust the frequency of searching for factorable Q's if we are near the 
			//minimum Q for this polynomial
			if (abs(j2 - min_block) < 1)
				blockinit = s_init - ((max_bits - blockbits));
			else if (abs(j2 - min_block) < 2)
				blockinit = s_init - ((max_bits - blockbits)/2);
			else if (abs(j2 - min_block) < 5)
				blockinit = s_init - ((max_bits - blockbits)/4);
			else 
				blockinit = s_init;
				
			smpqs_sieve_block(sieve,fb_sieve_p,start_prime,blockinit,fb);

			i = smpqs_check_relations(sieve_interval,j2,sieve,n,poly,blockinit,fb_sieve_p,fb,full,partial,cutoff,small_bits,start_prime,0,&num,numpoly);

			smpqs_sieve_block(sieve,fb_sieve_n,start_prime,blockinit,fb);

			i = smpqs_check_relations(sieve_interval,j2,sieve,n,poly,blockinit,fb_sieve_n,fb,full,partial,cutoff,small_bits,start_prime,1,&num,numpoly);

			if (i >= check_total)
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
				else
				{
					//need to keep sieving.  since the last time we checked, we've found
					//(full->num_r + partial->act_r) - (last_numfull + last_numpartial)
					//relations.  assume we'll find this many next time (actually should be slightly more), and 
					//scale how much longer we need to sieve to hit the target.
					//we only count the full relations to determine when to check.

					//if the number expected to be found next time puts us over the needed amount, scale the 
					//check total appropriately.  otherwise, just increment check_total by check_inc
					if (((full->num_r + partial->act_r) + (last_numfull + last_numpartial)) > (fb->B + sieve_params.num_extra_relations))
					{
						num_needed = (fb->B + sieve_params.num_extra_relations) - (partial->act_r + full->num_r);
						num_expected = (last_numfull + last_numpartial);
						check_total += (uint32)((double)check_inc * (double)num_needed / (double)num_expected);
						check_total += sieve_params.num_extra_relations;
					}
					else
						check_total += check_inc;
				}
				last_numfull = full->num_r;
				last_numpartial = partial->act_r;
			}
		}
		
		//next polynomial
		smpqs_nextD(&poly->poly_d,n);
		smpqs_computeB(poly,n);
		smpqs_computeRoots(poly,fb,2);

		numpoly++;
	}

done:

	//can free sieving structures now
	free(sieve);
	free(fb_sieve_p);
	free(fb_sieve_n);
	zFree(&poly->poly_a);
	zFree(&poly->poly_b);
	zFree(&poly->poly_c);
	zFree(&poly->poly_d);
	free(poly);

	num_factors=0;
	factors = (z *)malloc(MAX_FACTORS * sizeof(z));
	for (i=0;i<MAX_FACTORS;i++)
		zInit(&factors[i]);

	i = BlockGauss_MPQS(full,partial,apoly,bpoly,fb,n,mul,sieve_log,factors,&num_factors,0);

	for(i=0;i<num_factors;i++)
	{
		zDiv(n,&factors[i],&tmp,&tmp2);
		zCopy(&tmp,n);
		//non-trivial factor found
		if (isPrime(&factors[i]))
		{
			factors[i].type = PRP;
			//add_to_factor_list(fobj, &factors[i]);
		}
		else
		{
			factors[i].type = COMPOSITE;
			//add_to_factor_list(fobj, &factors[i]);
		}
	}
	zShortDiv(n,mul,&tmp);
	zCopy(&tmp,n);

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

	free(fb->list);
	free(fb);

	return;
}

static int smpqs_check_relations(uint32 sieve_interval, uint32 blocknum, uint8 sieve[], z *n, mpqs_poly *poly, uint8 closnuf,
						sieve_fb *fb,fb_list_mpqs *fullfb, mpqs_rlist *full, mpqs_rlist *partial, 
						uint32 cutoff, uint8 small_cutoff, uint32 start_prime, uint32 parity, uint32 *num, int numpoly)
{
	z Q,t1,t2,t3;
	uint32 offset,j,k;
	uint32 neg;
	uint64 *sieveblock;
	uint64 mask = ((uint64)0x80808080 << 32) | (uint64)0x80808080;
	uint32 limit = BLOCKSIZE >> 3;
	
	sieveblock = (uint64 *)sieve;
	zInit(&Q);
	zInit(&t1);
	zInit(&t2);
	zInit(&t3);
	//check for relations
	for (j=0;j<limit;j++)
	{
		if ((sieveblock[j] & mask) == (uint64)(0))
			continue;

		//else figure out which one's need to be checked
		for (k=0;k<8;k++)
		{
			if ((sieve[8*j + k] & 0x80) == 0)
				continue;

			(*num)++;

			//this one is close enough, compute 
			//Q(x) = (ax + b)^2 - N, where x is the sieve index
			//(a*x +/- 2b)*x + c;
			//offset = blocknum*BLOCKSIZE + 8*j+k;
			//printf("computing Q for block %d, offset %d\n",blocknum,8*j+k);
			offset = (blocknum<<BLOCKBITS) + (j<<3) + k;
			zShiftLeft(&t2,&poly->poly_b,1);

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

			smpqs_trial_divide_Q(&Q,fb,full,partial,sieve,offset,8*j+k,neg,fullfb,cutoff,small_cutoff,start_prime,numpoly,parity,closnuf,sieve[8*j+k]);
		}
	}
	zFree(&Q);
	zFree(&t1);
	zFree(&t2);
	zFree(&t3);
	return full->num_r;
}

static void smpqs_trial_divide_Q(z *Q, sieve_fb *fb, mpqs_rlist *full, mpqs_rlist *partial,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, fb_list_mpqs *fullfb, uint32 cutoff,
						  uint8 small_cutoff, uint32 start_prime, int numpoly, uint32 parity,uint8 closnuf,uint8 bits)
{
	sieve_fb *fbptr;
	uint32 i,num_f,num_p;
	uint32 root1,root2,prime,B=fullfb->B, med_B=fullfb->med_B;
	uint32 r;
	int smooth_num;
	uint16 fboffset[MAX_SMOOTH_PRIMES];
	uint8 logp;
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

	//printf("trial dividing Q\n");
	bits = (255 - bits) + closnuf + 1;
	smooth_num=0;
	if (start_prime > 2)
	{	
		//these are the primes we ignored in the small_prime variation.
		//trial divide by all of them, and if they sufficiently divide Q, proceed with
		//the rest of the trial division
		i=2;
		while (i < start_prime)
		{
			fbptr = fb + i;
			prime = fbptr->prime;
			root1 = fbptr->root1;
			root2 = fbptr->root2;
			logp = fbptr->logprime;
				
			if (((offset - root1) % prime) == 0)
			{
				do
				{
					fboffset[++smooth_num] = (uint16)i;
					zShortDiv(Q,prime,Q);
					bits += logp;
					r = zShortMod(Q,prime);
				} while (r == 0);
			}
			else if (((offset - root2) % prime) == 0)
			{
				do
				{
					fboffset[++smooth_num] = (uint16)i;
					zShortDiv(Q,prime,Q);
					bits += logp;
					r = zShortMod(Q,prime);
				} while (r == 0);
			}
			i++;
		}
	}

	if (bits < (closnuf + small_cutoff))
		return;

	//take care of powers of two
	while (!(Q->val[0] & 1))
	{
		zShiftRight(Q,Q,1);
		fboffset[++smooth_num] = 1;
	}

	i=start_prime;
	//do the primes less than the blocksize.  primes bigger than the blocksize can be handled
	//even more efficiently.
	//a couple of observations from jasonp:
	//if a prime divides Q(x), then this index (j) and either
	//root1 or root2 are on the same arithmetic progression.  this we can
	//test with a single precision mod operation
	//also, the way we have written this, if prime is larger than either root, it can't divide Q(x)
	//so we can avoid even the single precision mod much of the time
	//4/12 elusive bug fixed.   trial division by primes greater than Q is a no-no.
	while (i < med_B)
	{
		fbptr = fb + i;
		root1 = fbptr->root1 + BLOCKSIZE - j;		//after sieving a block, the root is updated for the start of the next block
		root2 = fbptr->root2 + BLOCKSIZE - j;		//get it back on the current block's progression

		prime = fbptr->prime;
		logp = fbptr->logprime;

		if (Q->size == 1)
		{
			if (Q->val[0] < prime)
				break;
		}

		if (root2 >= prime)
		{
			//r2 is bigger than prime, it could be on the progression, check it.
			if (!(root2 % prime))
			{
				//it is, so it will divide Q(x).  do so as many times as we can.
				do
				{
					fboffset[++smooth_num] = (uint16)i;
					zShortDiv(Q,prime,Q);
					sieve[j] += logp;
					r = zShortMod(Q,prime);
				} while (r == 0);
				if (sieve[j] == closnuf)
					goto done;
			}
			else if ((root1 >= prime) && (!(root1 % prime)))
			{
				//r2 was a bust, but root1 met the criteria.  divide Q(x).		
				do
				{
					fboffset[++smooth_num] = (uint16)i;
					zShortDiv(Q,prime,Q);
					sieve[j] += logp;
					r = zShortMod(Q,prime);
				} while (r == 0);
				if (sieve[j] == closnuf)
					goto done;
			}
		}
		i++;
	}

	//for the primes bigger than the blocksize, the progression for
	//each root can only touch the block once.  so we can do a simple comparison and
	//not divide at all most of the time
	while (i < B)
	{
		fbptr = fb + i;
		root1 = fbptr->root1 + BLOCKSIZE - j;		//after sieveing a block, the root is updated for the start of the next block
		root2 = fbptr->root2 + BLOCKSIZE - j;		//get it back on the current block's progression
		prime = fbptr->prime;
		logp = fbptr->logprime;

		if (Q->size == 1)
		{
			if (Q->val[0] < prime)
				break;
		}

		if ((root1 == prime) || (root2 == prime))
		{
			do
			{
				fboffset[++smooth_num] = (uint16)i;
				zShortDiv(Q,prime,Q);
				sieve[j] += logp;
				r = zShortMod(Q,prime);
			} while (r == 0);
			if (sieve[j] == closnuf)
				goto done;
		}
		i++;
	}

done:

	//check if it completely factored by looking at the unfactored portion in tmp
	if ((Q->size == 1) && (Q->val[0] == 1))
	{
		//printf("saving full...\n");
		//fflush(stdin);
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
	else if ((Q->size == 1)  && (Q->val[0] < cutoff))
	{
		//printf("saving partial...\n");
		smpqs_save_relation(partial,offset,Q->val[0],smooth_num+1,num_p,fboffset,numpoly,parity);
		if (partial->num_r == partial->allocated) 
		{
			//printf("\nreallocating partials\n");
			r = partial->allocated;
			partial->allocated *= 2;
			partial->list = (mpqs_r **)realloc(partial->list, 
					partial->allocated * sizeof(mpqs_r *));
		}
	}
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

static void smpqs_sieve_block(uint8 *sieve, sieve_fb *fb, uint32 start_prime, uint8 s_init, fb_list_mpqs *fullfb)
{
	uint32 prime, root1, root2, tmp;
	uint32 B=fullfb->B, small_B=fullfb->small_B, med_B=fullfb->med_B;
	uint32 i,j;
	uint8 *inner_sieve;
	uint8 logp;
	sieve_fb *fbptr;

	//initialize block
	memset(sieve,s_init,BLOCKSIZE);
	
	//get pointers pointed in the right direction
	inner_sieve = sieve;

	//break block into sections to keep in L1 cache
	//blocksize should be a multiple of inner_blocksize
	for (i=0;i<BLOCKSIZE;i+=INNER_BLOCKSIZE)
	{
		//for some small number of primes (those that fit in L1 cache/2
		for (j=start_prime; j<small_B; j++)
		{
			fbptr = fb + j;
			prime = fbptr->prime;
			root1 = fbptr->root1;
			root2 = fbptr->root2;
			logp = fbptr->logprime;

			while (root2 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] -= logp;
				inner_sieve[root2] -= logp;
				root1 += prime;
				root2 += prime;
			}

			//don't forget the last proot1[i], and compute the roots for the next block
			if (root1 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] -= logp;
				root1 += prime;
				//root1 will be bigger on the next iteration, switch them now
				tmp = root2;
				root2 = root1;
				root1 = tmp;
			}

			fbptr->root1 = root1 - INNER_BLOCKSIZE;
			fbptr->root2 = root2 - INNER_BLOCKSIZE;
		}
		//move inner_sieve to the right to paint a new stripe of the small fb primes in sieve
		inner_sieve += INNER_BLOCKSIZE;
	}
	
	//we've now filled the entire block with the small fb primes, proceed with the rest
	for (i=small_B;i<med_B;i++)
	{	
		fbptr = fb + i;
		prime = fbptr->prime;
		root1 = fbptr->root1;
		root2 = fbptr->root2;
		logp = fbptr->logprime;

		//sieve over both roots
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
			
		fbptr->root1 = root1 - BLOCKSIZE;
		fbptr->root2 = root2 - BLOCKSIZE;
	}

	for (i=med_B;i<B;i++)
	{	
		fbptr = fb + i;
		prime = fbptr->prime;
		root1 = fbptr->root1;
		root2 = fbptr->root2;
		logp = fbptr->logprime;

		if (root1 < BLOCKSIZE)
		{
			sieve[root1] -= logp;
			root1 += prime;

			if (root2 < BLOCKSIZE)
			{
				sieve[root2] -= logp;
				root2 += prime;
			}
			else
			{
				tmp=root2;
				root2=root1;
				root1=tmp;
			}
		}

		fbptr->root1 = root1 - BLOCKSIZE;
		fbptr->root2 = root2 - BLOCKSIZE;
	}

	return;
}

static void smpqs_test_block(uint8 *sieve, fb_element_mpqs *fb, uint32 start_prime)
{
	uint32 prime, root1, root2;
	uint32 i,j;
	uint8 *inner_sieve;
	uint8 logp;
	fb_element_mpqs *fbptr;
	
	//initialize block
	memset(sieve,0,BLOCKSIZE);

	//get pointers pointed in the right direction
	inner_sieve = sieve;

	//break block into sections to keep in L1 cache
	//blocksize should be a multiple of inner_blocksize
	for (i=0;i<BLOCKSIZE;i+=INNER_BLOCKSIZE)
	{
		//for some small number of primes (those that fit in L1 cache/2
		for (j=2; j<start_prime; j++)
		{
			fbptr = fb + j;
			prime = fbptr->prime;
			root1 = fbptr->proot1;
			root2 = fbptr->proot2;
			logp = fbptr->logprime;

			while (root2 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				inner_sieve[root2] += logp;
				root1 += prime;
				root2 += prime;
			}

			//don't forget the last proot1[i], and compute the roots for the next block
			if (root1 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				root1 += prime;
				//don't update the root pointers... this is just a test.
			}
		}
		//move inner_sieve to the right to paint a new stripe of the small fb primes in sieve
		inner_sieve += INNER_BLOCKSIZE;
	}
	return;
}

void smpqs_nextD(z *poly_d, z *n)
{
	//the current poly_d is passed in.
	//compare it to opt_d and find the next poly_d, return in poly_d.
	//I'd like this to alternate between d's above and below opt_d in magnitude, as then they'd stay closer
	//to opt_d for longer.  for now, just go in the postive direction
	z t1;

	zInit(&t1);
	//the current poly_d is prime, so check and increment to find the next one that meets all the criteria
	do 
	{
		/*
		if (poly_d->val[0] <= 0xFFFFFFFD)
			poly_d->val[0] +=2;	//to find the next odd number
		else
			zShortAddp(poly_d,2,poly_d);
			*/
		zNextPrime(poly_d,&t1,1);	//return the next prime number, which may be the current number
		zCopy(&t1,poly_d);
	} while ((zJacobi(n,poly_d) != 1) || (zShortMod(poly_d,4) != 3));

	zFree(&t1);
	return;
}

void smpqs_computeB(mpqs_poly *poly, z *n)
{
	//using poly_d, compute poly_b and poly_a = poly_d^2
	//int i;
	z tmp, t0, t1, t2, t3, t4, t5, t6, t7;

	zInit(&tmp);
	zInit(&t0);
	zInit(&t1);
	zInit(&t2);
	zInit(&t3);
	zInit(&t4);
	zInit(&t5);
	zInit(&t6);
	zInit(&t7);
	//poly_a = d^2.  we just found d.  also compute b using Hegel theorem and lifting.

	//first compute (d-3)/4
	zCopy(&poly->poly_d,&t1);
	if (t1.val[0] >= 3)
		t1.val[0] -= 3;
	else
	{
		zShortSub(&t1,3,&t3);
		zCopy(&t3,&t1);
	}
	zShortDiv(&t1,4,&t1);

	//then t0 = n^(d-3)/4 mod d
	zModExp(n,&t1,&poly->poly_d,&t0);	

	//t1 = h1 = n*t0 mod d
	zMul(&t0,n,&t3);		
	zDiv(&t3,&poly->poly_d,&t2,&t1);

	//t2 = t1^2
	zSqr(&t1,&t2);		

	//t3 = n - h1^2
	zSub(n,&t2,&t3);		

	//t4 = (n - h1^2)/d
	zDiv(&t3,&poly->poly_d,&t4,&t5);	

	//t4 = (n - h1^2)/d mod d
	zDiv(&t4,&poly->poly_d,&t3,&tmp);	
	zCopy(&tmp,&t4);

	//t5 = 2*h1;
	zShiftLeft(&t5,&t1,1);		

	//compute d-2 for the inversion
	zCopy(&poly->poly_d,&tmp);
	if (tmp.val[0] >= 2)			
		tmp.val[0] -= 2;
	else
	{
		zShortSub(&tmp,2,&t7);
		zCopy(&t7,&tmp);
	}

	//compute t6 = (2*h1)^-1 mod d = (2*h1)^(d-2) mod d
	zModExp(&t5,&tmp,&poly->poly_d,&t6);	

	//compute t3 = h2 = ((2*h1)^-1 * (n - h1^2)/d) mod d
	zMul(&t6,&t4,&t5);		
	zDiv(&t5,&poly->poly_d,&t6,&t3);

	//compute t5 = h1 + h2*D
	zMul(&t3,&poly->poly_d,&t4);	
	zAdd(&t4,&t1,&t5);

	//we're now done with d, so compute a = d^2
	zSqr(&poly->poly_d,&poly->poly_a);		

	//compute b = h1 + h2*d mod a
	zDiv(&t5,&poly->poly_a,&t6,&poly->poly_b);	

	//make sure b < a/2
	zShiftRight(&t1,&poly->poly_a,1);
	if (zCompare(&poly->poly_b,&t1) > 0)
	{
		zSub(&poly->poly_a,&poly->poly_b,&t1);
		zCopy(&t1,&poly->poly_b);
	}

	//now that we have b, compute c = (b*b - n)/a
	zSqr(&poly->poly_b,&t1);
	zSub(&t1,n,&t2);
	zDiv(&t2,&poly->poly_a,&poly->poly_c,&t1);

	zFree(&tmp);
	zFree(&t0);
	zFree(&t1);
	zFree(&t2);
	zFree(&t3);
	zFree(&t4);
	zFree(&t5);
	zFree(&t6);
	zFree(&t7);
	return;
}

void smpqs_computeRoots(mpqs_poly *poly, fb_list_mpqs *fb, uint32 start_prime)
{
	//the roots are computed using a and b as follows:
	//(+/-t - b)(a)^-1 mod p
	//assume b > t
	//uint32 root1, root2, prime, amodp;
	int root1, root2, prime, amodp, bmodp, x;
	uint32 i;

	for (i=start_prime;i<fb->B;i++)
	{
		//take the fast method of computing the inverse from lenstra...
		root1 = fb->list[i].c1;
		root2 = fb->list[i].c2;
		prime = fb->list[i].prime;
		/*
		zShortSub(&poly->poly_b,root1,&w2);
		zShortSub(&poly->poly_b,root2,&w3);

		root1 = zShortMod(&w2,prime);
		root2 = zShortMod(&w3,prime);

		//but b is bigger than t, so we need to take the negative of these solutions.
		root1 = prime - root1;
		root2 = prime - root2;
		*/

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

		root1 = (uint32)((uint64)amodp * (uint64)root1 % (uint64)prime);
		root2 = (uint32)((uint64)amodp * (uint64)root2 % (uint64)prime);		

		if (root2 < root1)
		{
			fb->list[i].proot1 = root2;
			fb->list[i].proot2 = root1;
			fb->list[i].nroot1 = prime - root1;
			fb->list[i].nroot2 = prime - root2;
		}
		else
		{
			fb->list[i].proot1 = root1;
			fb->list[i].proot2 = root2;
			fb->list[i].nroot1 = prime - root2;
			fb->list[i].nroot2 = prime - root1;
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

