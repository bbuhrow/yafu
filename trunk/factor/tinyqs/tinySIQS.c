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
#include "gmp_xface.h"

int tiny_static_init(static_conf_t *sconf, mpz_t n);
int tiny_dynamic_init(dynamic_conf_t *dconf, static_conf_t *sconf);
uint32 make_tiny_fb(static_conf_t *sconf);
uint32 tiny_merge_data(dynamic_conf_t *dconf, static_conf_t *sconf);
void tiny_poly_a(static_conf_t *sconf, dynamic_conf_t *dconf);
void tiny_get_params(static_conf_t *sconf);
int qcomp_tiny(const void *x, const void *y);
int free_tiny(static_conf_t *sconf, dynamic_conf_t *dconf);
int tiny_BlockGauss(siqs_r *rlist, uint32 num_r,
			fb_list *fb, mpz_t n, int mul, 
			mpz_t *factors, uint32 *num_factor);

#define POLYA_DEBUG 1

// a slimmed down version of siqs, meant to run very quickly on tiny
// inputs (less than 115 bits, say).
void tinySIQS(mpz_t n, mpz_t *factors, uint32 *num_factors)
{
	// factor n into its prime factors, returning the factors and
	// how many of them we found

	//master control structure
	static_conf_t *static_conf;
	dynamic_conf_t *dconf;

	//some locals
	uint32 num_needed, num_found;
	uint32 j;
	int i;

	// we assume the input is not even, prime, a perfect power, or divisible
	// by small primes.  Just check the input size to see if we can
	// use squfof instead.
	if (mpz_sizeinbase(n, 2) > 115)
	{
		printf("input too big\n");
		return;
	}

	if (mpz_sizeinbase(n, 2) < 60)
	{
		mpz_t ztmp;
		mpz_init(ztmp);

		j = sp_shanks_loop(n, NULL);	
		if (j > 1)
		{
			mpz_set_64(ztmp, j);
			//add_to_factor_list(fobj, ztmp);

			mpz_tdiv_q_ui(n, n, j);
			//add_to_factor_list(fobj, fobj->qs_obj.gmp_n);

			mpz_set_ui(n, 1);			
		}
		
		return;
	}	

	//initialize the data objects
	static_conf = (static_conf_t *)malloc(sizeof(static_conf_t));
	dconf = (dynamic_conf_t *)malloc(sizeof(dynamic_conf_t));
	static_conf->obj = NULL;

	//fill in the factorization object
	static_conf->bits = mpz_sizeinbase(n, 2);
	static_conf->digits_n = gmp_base10(n);

	if (VFLAG > 0)
		printf("\nstarting tinySIQS on c%d: %s\n",static_conf->digits_n, 
		mpz_conv2str(&gstr1.s, 10, n));

	//get best parameters, multiplier, and factor base for the job
	//initialize and fill out the static part of the job data structure
	tiny_static_init(static_conf, n);
	static_conf->in_mem = 1;
	static_conf->is_tiny = 1;

	//allocate structures for use in sieving with threads
	tiny_dynamic_init(dconf, static_conf);

	//start the process
	num_needed = static_conf->factor_base->B + static_conf->num_extra_relations;
	num_found = 0;
	static_conf->total_poly_a = -1;

    while (num_found < num_needed) 
	{
		thread_sievedata_t t;		

		// generate a new poly A value for the thread we pulled out of the queue
		// using its dconf.  this is done by the master thread because it also 
		// stores the coefficients in a master list
		static_conf->total_poly_a++;
		tiny_poly_a(static_conf,dconf);

		//do some work
		t.dconf = dconf;
		t.sconf = static_conf;
		tiny_process_poly(&t);

		// combine partials and count total relations
		qsort(dconf->relation_buf, dconf->buffered_rels, sizeof(siqs_r), 
			&qcomp_tiny);

		static_conf->num_relations = 0;
		static_conf->num_cycles = 0;
		for (i=0; i < dconf->buffered_rels; i++)
		{
			if (dconf->relation_buf[i].large_prime[0] == 1)
				static_conf->num_relations++;
			else if ((i > 0) && 
				(dconf->relation_buf[i].large_prime[0] == 
				dconf->relation_buf[i-1].large_prime[0]))
				static_conf->num_cycles++;
		}
		num_found = static_conf->num_cycles + static_conf->num_relations;

		if (VFLAG > 0)
			printf("%d rels found: %d full + "
				"%d from %d partial (%d examined)\n",
				static_conf->num_cycles + static_conf->num_relations,
				static_conf->num_relations, static_conf->num_cycles, 
				dconf->buffered_rels - static_conf->num_relations,
				dconf->num);
	}	

	free_sieve(dconf);

	if (VFLAG > 0)
		printf("\n==== post processing stage (block gauss) ====\n");

	static_conf->factor_list.num_factors = 0;
	//tiny_BlockGauss(dconf->relation_buf, dconf->buffered_rels,
	//	static_conf->factor_base, n, static_conf->multiplier, 
	//	factors, num_factors);

	//free everything else
	free_tiny(static_conf, dconf);
	free(dconf);
	free(static_conf);

	return;
}

int free_tiny(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	uint32 i;

	//current poly info used during filtering
	free(sconf->curr_poly->gray);
	free(sconf->curr_poly->nu);
	free(sconf->curr_poly->qlisort);
	mpz_clear(sconf->curr_poly->mpz_poly_a);
	mpz_clear(sconf->curr_poly->mpz_poly_b);
	mpz_clear(sconf->curr_poly->mpz_poly_c);
	free(sconf->curr_poly);
	mpz_clear(sconf->curr_a);	
	free(sconf->modsqrt_array);
	align_free(sconf->factor_base->list->prime);
	align_free(sconf->factor_base->list->small_inv);
	align_free(sconf->factor_base->list->correction);
	align_free(sconf->factor_base->list->logprime);
	align_free(sconf->factor_base->tinylist->prime);
	align_free(sconf->factor_base->tinylist->small_inv);
	align_free(sconf->factor_base->tinylist->correction);
	align_free(sconf->factor_base->tinylist->logprime);
	align_free(sconf->factor_base->list);
	align_free(sconf->factor_base->tinylist);
	free(sconf->factor_base);

	//while freeing the list of factors, divide them out of the input
	for (i=0;i<sconf->factor_list.num_factors;i++)
	{
		mpz_t tmp;
		mpz_init(tmp);

		//convert the factor
		mp_t2gmp(&sconf->factor_list.final_factors[i]->factor,tmp);

		//divide it out
		mpz_tdiv_q(sconf->obj->qs_obj.gmp_n, sconf->obj->qs_obj.gmp_n, tmp);

		//log it
		add_to_factor_list(sconf->obj, tmp);
		
		mpz_clear(tmp);
		free(sconf->factor_list.final_factors[i]);
	}

	mpz_clear(sconf->sqrt_n);
	mpz_clear(sconf->n);
	mpz_clear(sconf->target_a);
	
	for (i=0; i<sconf->total_poly_a; i++)
		mpz_clear(sconf->poly_a_list[i]);
	free(sconf->poly_a_list);

	for (i=0; i<dconf->buffered_rels; i++)
		free(dconf->relation_buf[i].fb_offsets);
	free(dconf->relation_buf);

	return 0;
}

int qcomp_tiny(const void *x, const void *y)
{
	siqs_r *xx = (siqs_r *)x;
	siqs_r *yy = (siqs_r *)y;
	
	if (xx->large_prime[0] > yy->large_prime[0])
		return 1;
	else if (xx->large_prime[0] == yy->large_prime[0])
		return 0;
	else
		return -1;
}

void tiny_poly_a(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unpack stuff from the job data structure
	siqs_poly *poly = dconf->curr_poly;
	mpz_ptr target_a = sconf->target_a;
	fb_list *fb = sconf->factor_base;

	mpz_t tmp, tmp2, tmp3;
	mpz_ptr poly_a = poly->mpz_poly_a;
	int j, *qli = poly->qlisort, *s = &poly->s;
	uint32 i;
	uint32 upper_polypool_index, lower_polypool_index;

	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init(tmp3);

	//determine polypool indexes.  
	//this really should be done once after generating the factor base
	//these will be set to more appropriate values below
	lower_polypool_index = 2;
	upper_polypool_index = fb->small_B - 1;

	if (sconf->bits < 115)
	{
		uint32 hi = fb->B-1;
		uint32 lo = sconf->sieve_small_fb_start; //(fb->B-1) / 4;
		uint32 id;
		uint64 p;
		int doover, k;

		// we have to be careful not to pick primes that fall into the
		// SPV region.  This is because polya roots are not computed correctly
		// but we rely on special signals to the trial division engine to
		// not attempt the fast trial division methods on these primes, or
		// at least to guareentee that they will fail.  These special signals
		// only work on primes outside the SPV region, however, so we have
		// to limit our prime pool to those larger than the SPV cutoff.
#ifdef POLYA_DEBUG
		printf("lo pool index = %d, prime = %u\n", lo, fb->list->prime[lo]);
		printf("hi pool index = %d, prime = %u\n", hi, fb->list->prime[hi]);
#endif

		*s = 3;
		do {		
			// pick two factors randomly
#ifdef POLYA_DEBUG
			printf("building new candidate poly...\n");
			fflush(stdout);
#endif
			j = 0;
			mpz_set_ui(poly->mpz_poly_a, 1);
			for (i=0; i<2; i++)
			{
				id = lo + (uint32)((double)(hi - lo) * (double)rand() / (double)RAND_MAX);

				// make sure they are unique
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
					break;

				mpz_mul_ui(poly->mpz_poly_a, poly->mpz_poly_a, fb->list->prime[id]);
				qli[j++] = id;
			}

			if (doover)
				continue;

			// then a 3rd to try to get close to target_a
			mpz_tdiv_q(tmp, target_a, poly->mpz_poly_a);
			p = mpz_get_ui(tmp);
			if (p > fb->list->prime[fb->B - 1])
			{
				id = 0 + (uint32)((double)(10 - 0) * (double)rand() / (double)RAND_MAX);
				id = fb->B - 1 - id;
				mpz_mul_ui(poly->mpz_poly_a, poly->mpz_poly_a, 
					fb->list->prime[id]);
#ifdef POLYA_DEBUG
				printf("final factor is above fb bound; chose id = %d, prime = %u\n",
					id, fb->list->prime[id]);
#endif
			}
			else if (p < fb->list->prime[lo])
			{
				id = lo + (uint32)((double)(10 - 0) * (double)rand() / (double)RAND_MAX);
				mpz_mul_ui(poly->mpz_poly_a, poly->mpz_poly_a, 
					fb->list->prime[id]);	
#ifdef POLYA_DEBUG
				printf("final factor is below SPV bound; chose id = %d, prime = %u\n",
					id, fb->list->prime[id]);
#endif
			}
			else
			{
				// quarter the rest of the fb and pick one randomly from the 
				// best region
				uint32 q = (hi - lo) / 4;
				if ((p > fb->list->prime[lo]) && 
					(p <= fb->list->prime[lo + q-1]))
				{
					id = lo + 
						(uint32)((double)(q) * (double)rand() / (double)RAND_MAX);
					mpz_mul_ui(poly->mpz_poly_a, poly->mpz_poly_a, 
						fb->list->prime[id]);
#ifdef POLYA_DEBUG
				printf("final factor from region 1; chose id = %d, prime = %u\n",
					id, fb->list->prime[id]);
#endif
				}
				else if ((p > fb->list->prime[lo + q]) && 
					(p <= fb->list->prime[lo + q+q-1]))
				{
					id = lo + q + 
						(uint32)((double)(q) * (double)rand() / (double)RAND_MAX);
					mpz_mul_ui(poly->mpz_poly_a, poly->mpz_poly_a, 
						fb->list->prime[id]);
#ifdef POLYA_DEBUG
				printf("final factor from region 2; chose id = %d, prime = %u\n",
					id, fb->list->prime[id]);
#endif
				}
				else if ((p > fb->list->prime[lo + q+q]) && 
					(p <= fb->list->prime[lo + q+q+q-1]))
				{
					id = lo + (q+q) + 
						(uint32)((double)(q) * (double)rand() / (double)RAND_MAX);
					mpz_mul_ui(poly->mpz_poly_a, poly->mpz_poly_a, 
						fb->list->prime[id]);
#ifdef POLYA_DEBUG
				printf("final factor from region 3; chose id = %d, prime = %u\n",
					id, fb->list->prime[id]);
#endif
				}
				else
				{
					id = lo + (q+q+q) + 
						(uint32)((double)(q) * (double)rand() / (double)RAND_MAX);
					mpz_mul_ui(poly->mpz_poly_a, poly->mpz_poly_a, 
						fb->list->prime[id]);
#ifdef POLYA_DEBUG
				printf("final factor from region 4; chose id = %d, prime = %u\n",
					id, fb->list->prime[id]);
#endif
				}
			}

			// make sure it is unique
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

			qli[j++] = id;

			mpz_sub(tmp, target_a, poly->mpz_poly_a);
			if (mpz_sgn(tmp) < 0)
				mpz_neg(tmp, tmp);

		} while ( mpz_cmp(tmp, target_a) > 0);

#ifdef POLYA_DEBUG
		gmp_printf("Target A: %Zd\n A indices: %u %u %u\n A factors = %u, %u, %u\n A = %Zd\n", 
			target_a,
			poly->qlisort[0],
			poly->qlisort[1],
			poly->qlisort[2],
			fb->list->prime[poly->qlisort[0]], 
			fb->list->prime[poly->qlisort[1]], 
			fb->list->prime[poly->qlisort[2]], 
			poly->mpz_poly_a);
#endif

	}

	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(tmp3);

	//record this a in the list
	sconf->poly_a_list = (mpz_t *)realloc(sconf->poly_a_list,
		(sconf->total_poly_a + 1) * sizeof(mpz_t));
	mpz_init(sconf->poly_a_list[sconf->total_poly_a]);
	mpz_set(sconf->poly_a_list[sconf->total_poly_a], poly_a);

	//sort the indices of factors of 'a'
	qsort(poly->qlisort,poly->s,sizeof(int),&qcomp_int);
	memset(&poly->qlisort[poly->s], 255, (MAX_A_FACTORS - poly->s) * sizeof(int));	

#ifdef POLYA_DEBUG
		printf("done generating poly_a\n");
#endif

	return;
}

uint32 make_tiny_fb(static_conf_t *sconf)
{
	//finds the factor base primes, and computes the solutions to the congruence x^2 = N mod p
	//for the QS, these are the starting positions of the sieve relative to the sqrt of N.
	//for the MPQS, additional work using the polynomial coefficents and these congruences 
	//needs to be done to compute the starting positions of the sieve.

	//locals
	int i;
	uint32 b,j,r,k;
	uint32 prime, root1, root2;
	uint8 logp;
	fp_digit f;
	uint32 shift = 24;

	//unpack stuff from static data structure
	fb_list *fb = sconf->factor_base;
	mpz_ptr n = sconf->n;
	uint32 mul = sconf->multiplier;
	uint32 *modsqrt = sconf->modsqrt_array;

	//the 0th and 1st elements in the fb are always 1 and 2, so start searching with 3
	j=2; i=1;
	while (j<fb->B)
	{
		prime = (uint32)spSOEprimes[i];
		r = mpz_tdiv_ui(n, prime);
		if (r == 0)
		{
			if (mul % prime != 0)
			{
				//prime doesn't divide the multiplier, so
				//this prime divides the input, divide it out and bail
				mpz_tdiv_q_ui(n, n, prime);
				return prime;
			}

			//p divides n, which means it divides the multiplier.
			//we can still use it, but it only has one solution to x^2 == n mod p instead
			//of two.  just divide its logprime in half.
			//we also can't find the root using shanks-tonelli, but it will be very small
			//because the multiplier is very small, so just use brute force.
			b = mpz_tdiv_ui(n, prime);
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
			//the shift determines how big the dividends can be - it should be 
			//greater than the max possible dividend.  The shift should not be too
			//big though.  Max word size (inverse * dividend) >> shift should leave
			//more bits than the max divisor.
			if (prime < 256)
			{
				fb->tinylist->prime[j] = prime;
				fb->tinylist->logprime[j] = logp;

				fb->tinylist->small_inv[j] = (uint32)(((uint64)1 << 32) / (uint64)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->tinylist->small_inv[j]) {
					fb->tinylist->correction[j] = 1;
				}
				else {
					fb->tinylist->correction[j] = 0;
					fb->tinylist->small_inv[j]++;
				}
			}
			else
			{
				if ((shift == 24) && 
					(prime > 1024) && 
					(j % 8 == 0))
					shift = 26;

				if ((shift == 26) && 
					(prime > 4096) && 
					(j % 8 == 0))
					shift = 28;

				fb->list->small_inv[j] = (uint16)(((uint32)1 << shift) / prime);
				if (floor((double)(1 << shift) / (double)prime + 0.5) ==
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

		b = jacobi_1((fp_digit)r,(fp_digit)prime);
		if (b==1)
		{
			//this prime works
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
			//this is very fragile... need better range checking.
			if (prime < 256)
			{
				fb->tinylist->prime[j] = prime;
				fb->tinylist->logprime[j] = logp;

				fb->tinylist->small_inv[j] = (uint32)(((uint64)1 << 32) / (uint64)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->tinylist->small_inv[j]) {
					fb->tinylist->correction[j] = 1;
				}
				else {
					fb->tinylist->correction[j] = 0;
					fb->tinylist->small_inv[j]++;
				}
			}
			else
			{
				if ((shift == 24) && 
					(prime > 1024) && 
					(j % 8 == 0))
					shift = 26;

				if ((shift == 26) && 
					(prime > 4096) && 
					(j % 8 == 0))
					shift = 28;

				fb->list->small_inv[j] = (uint16)(((uint32)1 << shift) / prime);
				if (floor((double)(1 << shift) / (double)prime + 0.5) ==
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

	return 0;
}

#define TINY_PARAM_ROWS 8
void tiny_get_params(static_conf_t *sconf)
{
	int bits,i;
	double scale;
	fb_list *fb = sconf->factor_base;

	//parameter table
	//bits, fb primes, lp mulitplier, 64k blocks
	//adjustment in v1.27 - more primes and less blocks for numbers > ~80 digits
	//also different scaling for numbers bigger than 100 digits (constant increase
	//of 20% per line)
	int param_table[TINY_PARAM_ROWS][4] = {
		{50,	30,	30,	1},
		{60,	36,	40,	1},
		{70,	50,	40,	1},
		{80,	80,	40,	1},
		{90,	120,	40,	1},
		{100,	175,	50,	1},
		{110,	275,	50,	1},	
		{120,	375,	50,	1},
	};

	//linear interpolation according to bit size to determine
	//factor base bound.  use the closest parameter for lp multiplier
	//and number of blocks.
	bits = sconf->bits;

	fb->B = 0;
	if (bits <= param_table[0][0])
	{
		scale = (double)bits / (double)param_table[0][0];
		fb->B = (uint32)(scale * (double)(param_table[0][1]));		
		sconf->large_mult = 40;
		sconf->num_blocks = 1;
	}
	else
	{
		for (i=0;i<TINY_PARAM_ROWS;i++)
		{
			if (bits > param_table[i][0] && bits <= param_table[i+1][0])
			{
				scale = (double)(param_table[i+1][0] - bits) /
					(double)(param_table[i+1][0] - param_table[i][0]);
				fb->B = param_table[i+1][1] - 
					(uint32)(scale * (double)(param_table[i+1][1] - param_table[i][1]));
				
				sconf->large_mult = (uint32)((param_table[i+1][2] + param_table[i][2])/2.0 + 0.5);
				sconf->num_blocks = (uint32)((param_table[i+1][3] + param_table[i][3])/2.0 + 0.5);
			}
		}
	}

	if (fb->B == 0)
	{
		//off the end of the table, extrapolate based on the slope of 
		//the last two

		scale = (double)(param_table[TINY_PARAM_ROWS-1][1] - param_table[TINY_PARAM_ROWS-2][1]) /
			(double)(param_table[TINY_PARAM_ROWS-1][0] - param_table[TINY_PARAM_ROWS-2][0]);
		fb->B = (uint32)(((double)bits - param_table[TINY_PARAM_ROWS-1][0]) * 
			scale + param_table[TINY_PARAM_ROWS-1][1]);
		sconf->large_mult = param_table[TINY_PARAM_ROWS-1][2];	//reuse last one

		scale = (double)(param_table[TINY_PARAM_ROWS-1][3] - param_table[TINY_PARAM_ROWS-2][3]) /
			(double)(param_table[TINY_PARAM_ROWS-1][0] - param_table[TINY_PARAM_ROWS-2][0]);
		sconf->num_blocks = (uint32)(((double)bits - param_table[TINY_PARAM_ROWS-1][0]) * 
			scale + param_table[TINY_PARAM_ROWS-1][3]);

	}

	// minimum factor base - for use with really small inputs.
	// not efficient, but needed for decent poly selection
	//if (fb->B < 250)
//		fb->B = 250;

	return;
}

int tiny_static_init(static_conf_t *sconf, mpz_t n)
{
	//find the best parameters, multiplier, and factor base
	//for the input.  This is an iterative job because we may find
	//a factor of the input while finding the factor base, in which case
	//we log that fact, and start over with new parameters, multiplier etc.
	//also allocate space for things that don't change during the sieving
	//process.  this scratch space is shared among all threads.

	uint32 i;
	uint32 closnuf;
	double sum, avg, sd;

	//default parameters
	sconf->fudge_factor = 1.3;
	sconf->large_mult = 0; //30;
	sconf->num_blocks = 40;
	sconf->num_extra_relations = 32;
	sconf->small_limit = 256;
	sconf->use_dlp = 0;

	// sieve core functions are fixed
	firstRoots_ptr = &firstRoots_32k;
	nextRoots_ptr = &nextRoots_32k;
	testRoots_ptr = &testfirstRoots_32k;
	med_sieve_ptr = &med_sieveblock_32k;
	tdiv_med_ptr = &tdiv_medprimes_32k;
	resieve_med_ptr = &resieve_medprimes_32k;
	sconf->qs_blocksize = 32768;
	sconf->qs_blockbits = 15;

	//allocate the space for the factor base structure
	sconf->factor_base = (fb_list *)malloc(sizeof(fb_list));

	//allocate space for a copy of input number in the job structure
	mpz_init(sconf->n);
	mpz_set(sconf->n, n);

	//initialize some constants
	mpz_init(sconf->sqrt_n);
	mpz_init(sconf->target_a);

	//initialize the bookkeeping for tracking partial relations
	sconf->components = 0;
	sconf->vertices = 0;
	sconf->num_cycles = 0;
	sconf->num_relations = 0;

	//look up some parameters - tuned for 64k L1 cache
	tiny_get_params(sconf);

	//the number of primes in the factor base is rounded up to the 
	//next multiple of 16, so that we can use aligned moves to speed
	//computations of root updates
	sconf->factor_base->B += (16 - (sconf->factor_base->B % 16));

	//allocate the space for the factor base elements
	sconf->factor_base->list = (fb_element_siqs *)xmalloc_align(
		(size_t)(sizeof(fb_element_siqs)));
	sconf->factor_base->tinylist = (tiny_fb_element_siqs *)xmalloc_align(
		(size_t)(sizeof(tiny_fb_element_siqs)));

	sconf->modsqrt_array = (uint32 *)malloc(
		sconf->factor_base->B * sizeof(uint32));
	sconf->factor_base->list->prime = (uint32 *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(uint32)));

	sconf->factor_base->tinylist->prime = (uint32 *)xmalloc_align(
		(size_t)(256 * sizeof(uint32)));
	sconf->factor_base->tinylist->small_inv = (uint32 *)xmalloc_align(
		(size_t)(256 * sizeof(uint32)));
	sconf->factor_base->tinylist->correction = (uint32 *)xmalloc_align(
		(size_t)(256 * sizeof(uint32)));
	sconf->factor_base->tinylist->logprime = (uint32 *)xmalloc_align(
		(size_t)(256 * sizeof(uint32)));

	sconf->factor_base->list->small_inv = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(uint16)));
	sconf->factor_base->list->correction = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(uint16)));
		
	sconf->factor_base->list->logprime = (uint32 *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(uint32)));

	//find multiplier
	sconf->multiplier = (uint32)choose_multiplier_siqs(sconf->factor_base->B, sconf->n);
	mpz_mul_ui(sconf->n, sconf->n, sconf->multiplier);

	//sconf holds n*mul, so update its digit count and number of bits
	sconf->digits_n = gmp_base10(sconf->n);
	sconf->bits = mpz_sizeinbase(sconf->n, 2);

	//find sqrt_n
	mpz_sqrt(sconf->sqrt_n, sconf->n);

	//construct the factor base - input will not have any small factors
	make_tiny_fb(sconf);

	//the first two fb primes are always the same
	sconf->factor_base->list->prime[0] = 1;		//actually represents -1
	sconf->factor_base->list->prime[1] = 2;
	sconf->num_blocks = 1;

	//set the sieve interval.  this many blocks on each side of 0
	sconf->sieve_interval = sconf->qs_blocksize;

	//compute sieving limits
	sconf->factor_base->small_B = MIN(
		sconf->factor_base->B,1024); //((INNER_BLOCKSIZE)/(sizeof(sieve_fb))));

	//test the contribution of the small primes to the sieve.  
	for (i = 2; i < sconf->factor_base->B; i++)
	{
		if (sconf->factor_base->list->prime[i] > sconf->small_limit)
			break;
	}
	sconf->sieve_small_fb_start = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 4 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 1024)  &&
			(i % 8 == 0)) break;
	}
	sconf->factor_base->fb_10bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 4 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 2048)  &&
			(i % 8 == 0)) break;
	}
	sconf->factor_base->fb_11bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 4 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 4096)  &&
			(i % 8 == 0)) break;
	}
	sconf->factor_base->fb_12bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 4 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 8192)  &&
			(i % 8 == 0)) break;
	}
	sconf->factor_base->fb_13bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 14 bits.  
		//wait until the index is a multiple of 4 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > 16384)  &&
			(i % 8 == 0)) 
		{
			i -= 8;
			break;
		}
	}	
	sconf->factor_base->fb_14bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 15 bits.  
		//wait until the index is a multiple of 4 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > 32768)  &&
			(i % 8 == 0)) break;
	}
	sconf->factor_base->fb_15bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceed the blocksize.  
		//wait until the index is a multiple of 16 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > (uint32)(1.5 * (double)sconf->qs_blocksize))  &&
			(i % 16 == 0))
			break;

		//or 2^16, whichever is smaller
		if (sconf->factor_base->list->prime[i] > 65536)
		{
			i -= i%16;
			break;
		}

		//or, of course, the entire factor base (loop condition)
	}
	sconf->factor_base->med_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceed the size of the sieve 
		//interval.  wait until the index is a multiple of 16 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > sconf->sieve_interval) &&
			(i % 16 == 0))
		{
			i -= 16;
			break;
		}
	}
	sconf->factor_base->large_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		if (sconf->factor_base->list->prime[i] > 2*sconf->sieve_interval)
			break;
	}
	sconf->factor_base->x2_large_B = i;

	//a couple limits
	sconf->pmax = sconf->factor_base->list->prime[sconf->factor_base->B-1];
	sconf->large_prime_max = sconf->pmax * sconf->large_mult;

	//based on the size of the input, determine how to proceed.
	scan_ptr = &check_relations_siqs_1;
	sconf->scan_unrolling = 8;
	sconf->use_dlp = 0;

	//'a' values should be as close as possible to sqrt(2n)/M in order to make
	//values of g_{a,b}(x) as uniform as possible
	mpz_mul_2exp(sconf->target_a, sconf->n, 1);
	mpz_sqrt(sconf->target_a, sconf->target_a);
	mpz_tdiv_q_ui(sconf->target_a, sconf->target_a, sconf->sieve_interval); 

	//compute the number of bits in M/2*sqrt(N/2), the approximate value
	//of residues in the sieve interval.  Then subtract some slack.
	//sieve locations greater than this are worthy of trial dividing
	closnuf = (uint8)(double)((sconf->bits - 1)/2);
	closnuf += (uint8)(log((double)sconf->sieve_interval/2)/log(2.0));
	closnuf -= (uint8)(sconf->fudge_factor * log(sconf->large_prime_max) / log(2.0));

	//contribution of all small primes we're skipping to a block's
	//worth of sieving... compute the average per sieve location
	sum = 0;
	for (i = 2; i < sconf->sieve_small_fb_start; i++)
	{
		uint32 prime = sconf->factor_base->list->prime[i];

		sum += (double)sconf->factor_base->list->logprime[i] * 
			(double)sconf->qs_blocksize / (double)prime;
	}
	avg = 2*sum/sconf->qs_blocksize;
	//this was observed to be the typical standard dev. for one block of 
	//test sieving... wrap this magic number along with several others into
	//one empirically determined fudge factor...
	sd = sqrt(28);

	//this appears to work fairly well... paper mentioned doing it this
	//way... find out and reference here.
	sconf->tf_small_cutoff = (uint8)(avg + 2.5*sd);
	closnuf -= sconf->tf_small_cutoff;	//correction to the previous estimate

	sconf->blockinit = closnuf;
	sconf->tf_closnuf = closnuf;

	//needed during filtering
	mpz_init(sconf->curr_a);
	sconf->curr_poly = (siqs_poly *)malloc(sizeof(siqs_poly));
	mpz_init(sconf->curr_poly->mpz_poly_a); //, sconf->bits);
	mpz_init(sconf->curr_poly->mpz_poly_b); //, sconf->bits);
	mpz_init(sconf->curr_poly->mpz_poly_c); //, sconf->bits);
	sconf->curr_poly->qlisort = (int *)malloc(MAX_A_FACTORS*sizeof(int));
	sconf->curr_poly->gray = (char *) malloc( 65536 * sizeof(char));
	sconf->curr_poly->nu = (char *) malloc( 65536 * sizeof(char));

	//initialize a list of all poly_a values used 
	sconf->poly_a_list = (mpz_t *)malloc(sizeof(mpz_t));

	sconf->total_poly_a = 0;	//track number of A polys used
	sconf->num_r = 0;			//total relations found
	
	sconf->tot_poly = 0;		//track total number of polys
	sconf->num = 0;				//sieve locations subjected to trial division

	//no factors so far...
	sconf->factor_list.num_factors = 0;

	return 0;
}

int tiny_dynamic_init(dynamic_conf_t *dconf, static_conf_t *sconf)
{
	//allocate the dynamic structure which hold scratch space for 
	//various things used during sieving.
	uint32 i;

	//workspace bigints
	mpz_init(dconf->gmptmp1); 
	mpz_init(dconf->gmptmp2);
	mpz_init(dconf->gmptmp3);

	//this stuff changes with every new poly
	//allocate a polynomial structure which will hold the current
	//set of polynomial coefficients (a,b,c) and other info
	dconf->curr_poly = (siqs_poly *)malloc(sizeof(siqs_poly));
	mpz_init(dconf->curr_poly->mpz_poly_a);
	mpz_init(dconf->curr_poly->mpz_poly_b);
	mpz_init(dconf->curr_poly->mpz_poly_c);
	dconf->curr_poly->qlisort = (int *)malloc(MAX_A_FACTORS*sizeof(int));
	dconf->curr_poly->gray = (char *) malloc( 65536 * sizeof(char));
	dconf->curr_poly->nu = (char *) malloc( 65536 * sizeof(char));

	// point the dynamic relation storage to the static storage - we'll
	// be storing them all in memory...
	dconf->relation_buf = (siqs_r *)malloc(32768 * sizeof(siqs_r));
	dconf->buffered_rel_alloc = 32768;
	dconf->buffered_rels = 0;

	//allocate the sieving factor bases
	dconf->comp_sieve_p = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));
	dconf->comp_sieve_n = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));

	dconf->comp_sieve_p->prime = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_p->root1 = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_p->root2 = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_p->logp = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));

	dconf->comp_sieve_n->prime = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_n->root1 = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_n->root2 = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_n->logp = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));

	dconf->fb_sieve_p = (sieve_fb *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(sieve_fb)));
	dconf->fb_sieve_n = (sieve_fb *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(sieve_fb)));
	
	dconf->update_data.sm_firstroots1 = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->update_data.sm_firstroots2 = (uint16 *)xmalloc_align(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->update_data.firstroots1 = (int *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(int)));
	dconf->update_data.firstroots2 = (int *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(int)));
	dconf->update_data.prime = (uint32 *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(uint32)));
	dconf->update_data.logp = (uint8 *)xmalloc_align(
		(size_t)(sconf->factor_base->B * sizeof(uint8)));
	dconf->rootupdates = (int *)xmalloc_align(
		(size_t)(MAX_A_FACTORS * sconf->factor_base->B * sizeof(int)));
	dconf->sm_rootupdates = (uint16 *)xmalloc_align(
		(size_t)(MAX_A_FACTORS * sconf->factor_base->B * sizeof(uint16)));
	
	//allocate the sieve
	dconf->sieve = (uint8 *)xmalloc_align(
		(size_t) (sconf->qs_blocksize * sizeof(uint8)));

	//allocate the Bl array, space for MAX_Bl bigint numbers
	dconf->Bl = (mpz_t *)malloc(MAX_A_FACTORS * sizeof(mpz_t));
	for (i=0;i<MAX_A_FACTORS;i++)
		mpz_init(dconf->Bl[i]);

	//copy the unchanging part to the sieving factor bases
	for (i = 2; i < sconf->factor_base->med_B; i++)
	{
		uint32 p = sconf->factor_base->list->prime[i];
		uint32 lp = sconf->factor_base->list->logprime[i];

		dconf->comp_sieve_p->logp[i] = (uint8)lp;
		dconf->comp_sieve_p->prime[i] = (uint16)p;
		dconf->comp_sieve_n->logp[i] = (uint8)lp;
		dconf->comp_sieve_n->prime[i] = (uint16)p;

		dconf->fb_sieve_p[i].prime = p;
		dconf->fb_sieve_p[i].logprime = lp;
		dconf->fb_sieve_n[i].prime = p;
		dconf->fb_sieve_n[i].logprime = lp;
		dconf->update_data.prime[i] = p;
		dconf->update_data.logp[i] = lp;
	}

	for (; i < sconf->factor_base->B; i++)
	{
		dconf->fb_sieve_p[i].prime = sconf->factor_base->list->prime[i];
		dconf->fb_sieve_p[i].logprime = sconf->factor_base->list->logprime[i];
		dconf->fb_sieve_n[i].prime = sconf->factor_base->list->prime[i];
		dconf->fb_sieve_n[i].logprime = sconf->factor_base->list->logprime[i];
		dconf->update_data.prime[i] = sconf->factor_base->list->prime[i];
		dconf->update_data.logp[i] = sconf->factor_base->list->logprime[i];
	}

	// we will not be using bucket sieving
	dconf->buckets = (lp_bucket *)malloc(sizeof(lp_bucket));
	dconf->buckets->list = NULL;
	dconf->buckets->alloc_slices = 0;
	dconf->buckets->num_slices = 0;

	//used in trial division to mask out the fb_index portion of bucket entries, so that
	//multiple block locations can be searched for in parallel using SSE2 instructions
	dconf->mask = (uint16 *)xmalloc_align(8 * sizeof(uint16));

	dconf->mask[1] = 0xFFFF;
	dconf->mask[3] = 0xFFFF;
	dconf->mask[5] = 0xFFFF;
	dconf->mask[7] = 0xFFFF;

#ifdef SSE2_RESIEVING
	dconf->corrections = (uint16 *)xmalloc_align(8 * sizeof(uint16));
#endif

	// array of sieve locations scanned from the sieve block that we
	// will submit to trial division.  make it the size of a sieve block 
	// in the pathological case that every sieve location is a report
	dconf->reports = (uint32 *)malloc(MAX_SIEVE_REPORTS * sizeof(uint32));
	dconf->num_reports = 0;
		
#ifdef USE_YAFU_TDIV
	dconf->Qvals32 = (z32 *)malloc(MAX_SIEVE_REPORTS * sizeof(z32));
	for (i=0; i<MAX_SIEVE_REPORTS; i++)
		zInit32(&dconf->Qvals32[i]);
#endif
	dconf->Qvals = (mpz_t *)malloc(MAX_SIEVE_REPORTS * sizeof(mpz_t));
	for (i=0; i<MAX_SIEVE_REPORTS; i++)
		mpz_init(dconf->Qvals[i]);

	dconf->valid_Qs = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));
	dconf->smooth_num = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));

	//initialize some counters
	dconf->tot_poly = 0;		//track total number of polys
	dconf->num = 0;				//sieve locations subjected to trial division

	return 0;
}

void *tiny_process_poly(void *ptr)
{
	//top level sieving function which performs all work for a single
	//new a coefficient.  has pthread calling conventions, meant to be
	//used in a multi-threaded environment
	thread_sievedata_t *thread_data = (thread_sievedata_t *)ptr;
	static_conf_t *sconf = thread_data->sconf;
	dynamic_conf_t *dconf = thread_data->dconf;

	//unpack stuff from the job data structure
	sieve_fb_compressed *fb_sieve_p = dconf->comp_sieve_p;	
	sieve_fb_compressed *fb_sieve_n = dconf->comp_sieve_n;
	siqs_poly *poly = dconf->curr_poly;
	uint8 *sieve = dconf->sieve;
	fb_list *fb = sconf->factor_base;
	lp_bucket *buckets = dconf->buckets;
	uint32 start_prime = sconf->sieve_small_fb_start;
	uint32 num_blocks = sconf->num_blocks;	
	uint8 blockinit = sconf->blockinit;

	//locals
	uint32 i;

	//this routine is handed a dconf structure which already has a
	//new poly a coefficient (and some supporting data).  continue from
	//there, first initializing the gray code...
	//update the gray code
	get_gray_code(dconf->curr_poly);

	//update roots, etc.
	dconf->maxB = 1<<(dconf->curr_poly->s-1);
	dconf->numB = 1;
	computeBl(sconf,dconf);

	firstRoots_ptr(sconf,dconf);

	//loop over each possible b value, for the current a value
	for ( ; dconf->numB < dconf->maxB; dconf->numB++, dconf->tot_poly++)
	{
		//setting these to be invalid means the last entry of every block will be ignored
		//so we're throwing away 1/blocksize relations, potentially, but gaining 
		//a faster sieve routine.
		uint32 invalid_root_marker = 0xFFFFFFFF;

		for (i=0; i < num_blocks; i++)
		{
			//set the roots for the factors of a such that
			//they will not be sieved.  we haven't found roots for them
			set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p, 1);
			med_sieve_ptr(sieve, fb_sieve_p, fb, start_prime, blockinit);
			lp_sieveblock(sieve, i, num_blocks, buckets, 0);

			//set the roots for the factors of a to force the following routine
			//to explicitly trial divide since we haven't found roots for them
			set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p, 0);
			scan_ptr(i,0,sconf,dconf);

			//set the roots for the factors of a such that
			//they will not be sieved.  we haven't found roots for them
			set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n, 1);
			med_sieve_ptr(sieve, fb_sieve_n, fb, start_prime, blockinit);
			lp_sieveblock(sieve, i, num_blocks, buckets, 1);

			//set the roots for the factors of a to force the following routine
			//to explicitly trial divide since we haven't found roots for them
			set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n, 0);
			scan_ptr(i,1,sconf,dconf);			

		}

		//next polynomial
		//use the stored Bl's and the gray code to find the next b
		nextB(dconf,sconf);
		//and update the roots
		nextRoots_ptr(sconf, dconf);

	}

	return 0;
}

static uint64 tiny_bitValRead64(uint64 **m, int row, int col);

int tiny_BlockGauss(siqs_r *rlist, uint32 num_r,
			fb_list *fb, mpz_t n, int mul, 
			mpz_t *factors, uint32 *num_factor)
{
	int i,j,k,l,a,q,polynum;
	int *bl;
	uint8 **m;		//matrix of the powers of the prime decompositions of the relations over the factor base
	uint64 **m2_64;	//m mod 2, packed into 32 bit words
	uint64 **aug_64;	//matrix to store the permutations of the rows of m2, packed into 32 bit words
	uint32 largep, bool_val, B = fb->B;
	uint32 *partial_index;
	int num_f,num_p;
	int num_col,num_col_aug,set_continue;
	const int blocksz = 64;
	uint64 *apoly, *bpoly;

	uint32 *pd;
	uint32 r;
	mpz_t zx, zy, tmp, tmp2, tmp3, tmp4, nn,tmp_a,input,zmul;
	
	mpz_init(zx);
	mpz_init(zy);
	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init(tmp3);
	mpz_init(tmp4);
	mpz_init(nn);
	mpz_init(input);
	mpz_init(tmp_a);
	mpz_init(zmul);

	// for every stored a_poly, derive the resulting b_poly's and 
	// store both a_poly and the b_poly's in a list

	num_col = (uint32)((B/blocksz)+1);
	num_col_aug = (uint32)(num_r/blocksz+1);

	//allocate storage based on total number of relations.
	pd = (uint32 *)malloc(B * sizeof(uint32));
	partial_index = (uint32 *)malloc(num_r * sizeof(uint32));

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
	for (i=0;i<num_r;i++)
	{
		if (rlist[i].large_prime[0] != 1)
			break;

		//Initialize
		for (j=0;j<(int)B;j++)
			m[i][j] = 0;

		//copy the pd's of the fboffsets to the correct location in m
		//offset 0 is special - indicates the parity of the offset
		m[i][0] = (uint8)rlist[i].fb_offsets[0];
		j=1;
		while (j<rlist[i].num_factors)
		{
			m[i][rlist[i].fb_offsets[j]]++;
			j++;
		}
	}
	num_f = i;

	//write fulls from partials to m, probably also redundant to do it this way?
	largep = rlist[i++].large_prime[0];
	j=num_f;
	for (; i < num_r; i++)
	{
		if (rlist[i].large_prime[0] == largep)
		{
			//this partial's largep is the same as the one before, add the pd's and copy to m
			for (k=0;k<(int)B;k++)
				m[j][k] = 0;

			//do the factor of -1
			m[j][0] = (uint8)rlist[i-1].fb_offsets[0];
			//then the rest
			k=1;
			while (k<rlist[i-1].num_factors) 
			{
				m[j][rlist[i-1].fb_offsets[k]]++;
				k++;
			}

			//factor of -1
			m[j][0] += rlist[i].fb_offsets[0];
			//the rest
			k=1;
			while (k<rlist[i].num_factors)
			{
				m[j][rlist[i].fb_offsets[k]]++;
				k++;
			}

			//remember the index of the partial that made this full relation, we'll need it later
			partial_index[j-num_f]=i;
			//increment the relation counter
			j++;
		}
		largep = rlist[i].large_prime[0];
	}
	num_p = j - num_f;

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
	mpz_tdiv_q_ui(input, n, mul); //zShortDiv(n,mul,&input);
	mpz_set_ui(zmul, mul); //sp2z(mul,&zmul);

	//initialize blacklist
	for (i=0;i<num_r;i++) bl[i] = 0;
	//search over all columns, right to left (more sparse on the right side)
	for (i=B-1;i>=0;i--)
	{
		//and all rows
		for (j=0;j<num_r;j++)
		{
			//if the j'th row, i'th bit is 1 and not blacklisted, continue
			bool_val = (tiny_bitValRead64(m2_64,j,i) != 0) && (bl[j] == 0);
			//bool_val = (((m2_64[(j)][(i >> 6)]) & (1ULL << ((uint64)i & 63ULL))) && (bl[j] == 0));
			if (bool_val)
			{
				//add the j'th row mod 2 to all rows after it with a 1 in the ith column
				for (k=j+1;k<num_r;k++)
				{
					bool_val = (tiny_bitValRead64(m2_64,k,i) != 0) && (bl[k] == 0);
					//bool_val = (((m2_64[(k)][(i >> 6)]) & (1ULL << ((uint64)i & 63ULL))) && (bl[k] == 0));
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
								bool_val = tiny_bitValRead64(aug_64,k,l) != 0;
								//bool_val = (((aug_64[(k)][(l >> 6)]) & (1ULL << ((uint64)l & 63ULL)))) != 0;
								if (bool_val)
								{
									//then the l'th row of m was involved
									for (q=0;q<(int)B;q++)
										pd[q] += m[l][q];
								}
							}

							//compute x mod n
							mpz_set_ui(zy, 1); //sm_zcopy(&zOne,&zy);
							mpz_set_ui(zx, 1); //sm_zcopy(&zOne,&zx);
							for (l=0;l<num_r;l++)
							{
								bool_val = tiny_bitValRead64(aug_64,k,l) != 0;
								//bool_val = (((aug_64[(k)][(l >> 6)]) & (1ULL << ((uint64)l & 63ULL)))) != 0;
								if (bool_val)
								{
									//printf("accumulating relation %d\n",l);
									//then the l'th relation is involved in the product of relations
									if (l >= num_f)
									{
										uint64 pa;
										uint32 d1,d2;
										//if l >= num_f, then this row refers to a relation generated from two partials.
										//we'll need to go back to the two partial locations to find the two offsets to
										//multiply together
										//luckily, we've remembered the index in the complete list of partials that 
										//created this full relation
											
										//our relation is of the form (ax + b)^2 == a(ax^2 + 2bx + c) mod n
										//(ax^2 + 2bx + c) is what we trial divided, and we remembered
										//a and b, so we can form the left hand side easily

										// poly_a index
										polynum = rlist[partial_index[l-num_f]].poly_idx >> 16;
											
										//compute Q1(x)
										pa = apoly[polynum];
										d1 = (uint32)sqrt((int64)pa);
										
										mpz_set_64(tmp, apoly[polynum]);
										mpz_set_64(tmp3, bpoly[polynum]);
										mpz_mul_ui(tmp, tmp, rlist[partial_index[l-num_f]].sieve_offset);
										if (rlist[partial_index[l-num_f]].parity)
											mpz_sub(tmp, tmp, tmp3); 
										else
											mpz_add(tmp, tmp, tmp3);

										//include 'a'
										mpz_mul_ui(tmp2, zy, d1);
										mpz_tdiv_r(zy, tmp2, n);

										//compute Q2(x)
										polynum = rlist[partial_index[l-num_f]-1].poly_idx >> 16;
											
										//compute Q1(x)
										pa = apoly[polynum]; 
										d2 = (uint32)sqrt((int64)pa);
										
										mpz_set_64(tmp3, apoly[polynum]);
										mpz_set_64(tmp4, bpoly[polynum]);
										mpz_mul_ui(tmp3, tmp3, rlist[partial_index[l-num_f]-1].sieve_offset); 
										if (rlist[partial_index[l-num_f]-1].parity)
											mpz_sub(tmp2, tmp3, tmp4); 
										else
											mpz_add(tmp2, tmp3, tmp4); 

										//compute Q(x1)*Q(x2)
										mpz_mul(tmp4, tmp, tmp2);
										mpz_mul(tmp, zx, tmp4);	//accumulate with previous terms
										mpz_tdiv_r(zx, tmp, n);	//mod n

										//include the large prime in mp_y
										mpz_mul_ui(tmp2, zy, rlist[partial_index[l-num_f]].large_prime[0]); 
										mpz_tdiv_r(zy, tmp2, n); 

										//include 'a'
										mpz_mul_ui(tmp2, zy, d2);	
										mpz_tdiv_r(zy, tmp2, n);
									}
									else
									{
										uint64 pa,d1;

										//recreate poly_b from poly_a (store instead??)
										polynum = rlist[l].poly_idx >> 16;
											
										//compute Q1(x)
										mpz_set_64(tmp, apoly[polynum]);
										mpz_set_64(tmp4, bpoly[polynum]);
										mpz_mul_ui(tmp, tmp, rlist[l].sieve_offset); //zShortMul(&apoly[polynum],full->list[l]->offset,&tmp);
										pa = apoly[polynum]; //apoly[polynum].val[0];
										d1 = (uint32)sqrt((int64)pa);

										if (rlist[l].parity)
											mpz_sub(nn, tmp, tmp4); //zShortSub(&tmp,pb,&nn);
										else
											mpz_add(nn, tmp, tmp4); //zShortAdd(&tmp,pb,&nn);

										mpz_mul(tmp, zx, nn); //zMul(&zx,&nn,&tmp);			//accumulate with previous terms
										mpz_tdiv_r(zx, tmp, n); //zDiv(&tmp,n,&tmp2,&zx);		//mod n

										mpz_mul_ui(tmp2, zy, d1); //zShortMul(&zy,d1,&tmp2);	//sqrt(a) = d is part of mp_y
										mpz_tdiv_r(zy, tmp2, n); //zDiv(&tmp2,n,&tmp3,&zy);
									}
								}
							}

							//compute y mod n
							//ignore the factor of -1 in this operation
							for (l=1;l<(int)B;l++)
							{
								if (pd[l] > 0)
								{
									mpz_set_ui(tmp, fb->list->prime[l]); //sp2z(fb->list->prime[l],&tmp);
									//pd tracks the exponents of the smooth factors.  we know they are all even
									//at this point.  we don't want to compute pd^2, so divide by 2.
									//computing the explicit exponentiation and then reducing is
									//slightly faster than doing modexp in smallmpqs.
									//zExp(pd[l]/2,&tmp,&tmp2);
									//zDiv(&tmp2,n,&tmp4,&tmp3);
									mpz_powm_ui(tmp3, tmp, pd[l] / 2, n);
									mpz_mul(tmp4, tmp3, zy); //zMul(&tmp3,&zy,&tmp4);
									mpz_tdiv_r(zy, tmp4, n); //zDiv(&tmp4,n,&tmp2,&zy);
								}
							}

							//split this off into a subroutine... also look for all non-trivial factors if one is composite
							//compute gcd(x-y,n)
							mpz_sub(tmp, zx, zy); //zSub(&zx,&zy,&tmp);
							mpz_gcd(nn, tmp, n); //zLEGCD(&tmp,n,&nn);

							//gmp_printf("gcd is %Zd\n",nn);
							
							/* remove any factors of the multiplier 
							   before saving tmp, and don't save at all
							   if tmp contains *only* multiplier factors */
							if (mul > 1) {
								uint32 ignore_me = spGCD(mul,
										mpz_tdiv_ui(nn, mul)); //zShortMod(&nn, mul));
								if (ignore_me > 1) {
									mpz_tdiv_q_ui(nn, nn, ignore_me); //zShortDiv(&nn, ignore_me, &tmp2);
									if (mpz_cmp_ui(nn, 1) == 0)
										continue;
								}								
							}								

							
							if ((mpz_cmp_ui(nn, 1) > 0) && (mpz_cmp(nn,input) < 0))
							{

								//gmp_printf("checking factor %Zd\n",nn);
								if (mpz_probab_prime_p(nn,5))
								{

									// sometime we find small primes that don't divide the input.
									// ignore these
									if (mpz_sizeinbase(nn,2) < 32)
									{
										if (mpz_tdiv_ui(input, mpz_get_ui(nn)) != 0)
											continue;

										//if (zShortMod(&input,nn.val[0]) != 0)
											//continue;
									}

									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (mpz_cmp(nn,factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									mpz_set(factors[*num_factor],nn);

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}

									//check if we're done by accumulating all factors and comparing to n
									mpz_set(nn,factors[0]);
									for (l=1;l<(int)*num_factor;l++)
										mpz_mul(nn,factors[l],nn); //,&nn);
									if (mpz_cmp(nn,input) == 0)
									{
										//found all factors, done
										goto free;
									}
								}

								//check the other factor
								//sm_zcopy(&input,&tmp);
								//zDiv(&tmp,&nn,&tmp2,&tmp3);
								//sm_zcopy(&tmp2,&tmp);
								
								mpz_tdiv_q(tmp, input, nn);								
	
								//gmp_printf("checking factor %Zd\n",tmp);
								if (mpz_probab_prime_p(tmp,5))
								{
									// sometime we find small primes that don't divide the input.
									// ignore these
									if (mpz_sizeinbase(nn,2) < 32)
									{
										if (mpz_tdiv_ui(input, mpz_get_ui(tmp)) != 0)
											continue;

										//if (zShortMod(&input,nn.val[0]) != 0)
											//continue;
									}

									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (mpz_cmp(tmp,factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									mpz_set(factors[*num_factor],tmp);

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}

									//check if we're done by accumulating all factors and comparing to n
									mpz_set(tmp,factors[0]);
									for (l=1;l<(int)*num_factor;l++)
										mpz_mul(tmp,factors[l],tmp); //,&nn);
									if (mpz_cmp(tmp,input) == 0)
									{
										//found all factors, done
										goto free;
									}
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
	r = mpz_tdiv_q_ui(tmp, n, mul); //r = (uint32)zShortDiv(n,mul,&tmp);
	for (i=0;(uint32)i<*num_factor;i++)
	{
		//sm_zcopy(&tmp,&nn);
		//zDiv(&nn,&factors[i],&tmp,&tmp2);
		mpz_tdiv_q(tmp, tmp, factors[i]);
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
	mpz_clear(zx);
	mpz_clear(zy);
	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
	mpz_clear(tmp4);
	mpz_clear(nn);
	mpz_clear(tmp_a);
	mpz_clear(input);
	mpz_clear(zmul);
	return 0;
}

static uint64 tiny_masks64[64] = {0x1,0x2,0x4,0x8,
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



static uint64 tiny_bitValRead64(uint64 **m, int row, int col)
{
	//col is the column in 0 to B-1 representation
	//read the bit in the packed 64 bit representation of the appropriate row
	//don't bother to check the bounds of m w.r.t row and col, assume caller knows what it's doing
	//return 0 if bit not set, 1 << bit offset otherwize
	int offset, mcol;
	mcol = col >> 6;
	offset = col & 63;
	return (m[row][mcol] & tiny_masks64[offset]);
}



