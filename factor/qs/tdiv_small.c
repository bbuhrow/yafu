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
#include "util.h"
#include "common.h"
#include "gmp_xface.h"

/*
We are given an array of bytes that has been sieved.  The basic trial 
division strategy is as follows:

1) Scan through the array and 'mark' locations that meet criteria 
indicating they may factor completely over the factor base.  

2) 'Filter' the marked locations by trial dividing by small primes
that we did not sieve.  These primes are all less than 256.  If after
removing small primes the location does not meet another set of criteria,
remove it from the 'marked' list (do not subject it to further trial
division).

3) Divide out primes from the factor base between 256 and 2^13 or 2^14, 
depending on the version (2^13 for 32k version, 2^14 for 64k).  

4) Resieve primes between 2^{13|14} and 2^16, max.  

5) Primes larger than 2^16 will have been bucket sieved.  Remove these
by scanning the buckets for sieve hits equal to the current block location.

6) If applicable/appropriate, factor a remaining composite with squfof

this file contains code implementing 2)

*/


#ifdef SPARSE_STORE
#define DIVIDE_ONE_PRIME(x) \
    do	\
    {	\
        mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);	\
        bits += logp;	\
    } while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0);
#else
#define DIVIDE_ONE_PRIME(x) \
	do	\
    	{	\
		dconf->fb_offsets[report_num][++smooth_num] = (x);	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);	\
		bits += logp;	\
	} while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0);
#endif

//#define DO_4X_SPV 1

void filter_SPV(uint8 parity, uint8 *sieve, uint32 poly_id, uint32 bnum, 
				static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	int i;
	uint32 bound, tmp, prime, root1, root2;
	int smooth_num;
	sieve_fb_compressed *fbc;
	tiny_fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->tinylist;
	uint8 logp, bits;
	uint32 tmp1, tmp2, tmp3, tmp4, offset, report_num;

	fullfb_ptr = fullfb;
	if (parity)
	{
		fbc = dconf->comp_sieve_n;
	}
	else
	{
		fbc = dconf->comp_sieve_p;
	}

	//the minimum value for the current poly_a and poly_b occur at offset (-b + sqrt(N))/a
	//make it slightly easier for a number to go through full trial division for
	//nearby blocks, since these offsets are more likely to factor over the FB.
	mpz_sub(dconf->gmptmp1, sconf->sqrt_n, dconf->curr_poly->mpz_poly_b);
	mpz_tdiv_q(dconf->gmptmp1, dconf->gmptmp1, dconf->curr_poly->mpz_poly_a);
	if (mpz_sgn(dconf->gmptmp1) < 0)
		mpz_neg(dconf->gmptmp1, dconf->gmptmp1);

	mpz_tdiv_q_2exp(dconf->gmptmp1, dconf->gmptmp1, sconf->qs_blockbits);
	if (abs(bnum - mpz_get_ui(dconf->gmptmp1)) == 0)
		dconf->tf_small_cutoff = sconf->tf_small_cutoff - 5;
	else if (abs(bnum - mpz_get_ui(dconf->gmptmp1)) == 1)
		dconf->tf_small_cutoff = sconf->tf_small_cutoff - 3;
	else 
		dconf->tf_small_cutoff = sconf->tf_small_cutoff;

#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	for (report_num = 0; report_num < dconf->num_reports; report_num++)
	{
		uint64 q64;

		//this one qualifies to check further, log that fact.
		dconf->num++;

		smooth_num = -1;

		//this one is close enough, compute 
		//Q(x)/a = (ax + b)^2 - N, where x is the sieve index
		//Q(x)/a = (ax + 2b)x + c;	
		offset = (bnum << sconf->qs_blockbits) + dconf->reports[report_num];

		//multiple precision arithmetic.  all the qstmp's are a global hack
		//but I don't want to Init/Free millions of times if I don't have to.

		mpz_mul_2exp(dconf->gmptmp2, dconf->curr_poly->mpz_poly_b, 1); 
		mpz_mul_ui(dconf->gmptmp1, dconf->curr_poly->mpz_poly_a, offset);

		if (parity)
			mpz_sub(dconf->Qvals[report_num], dconf->gmptmp1, dconf->gmptmp2); 
		else
			mpz_add(dconf->Qvals[report_num], dconf->gmptmp1, dconf->gmptmp2); 

		mpz_mul_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], offset); 
		mpz_add(dconf->Qvals[report_num], dconf->Qvals[report_num], dconf->curr_poly->mpz_poly_c); 

		if (mpz_sgn(dconf->Qvals[report_num]) < 0)
		{
			mpz_neg(dconf->Qvals[report_num], dconf->Qvals[report_num]);
		}

		//we have two signs to worry about.  the sign of the offset tells us how to calculate ax + b, while
		//the sign of Q(x) tells us how to factor Q(x) (with or without a factor of -1)
		//the square root phase will need to know both.  fboffset holds the sign of Q(x).  the sign of the 
		//offset is stored standalone in the relation structure.
	
		//compute the bound for small primes.  if we can't find enough small
		//primes, then abort the trial division early because it is likely to fail to
		//produce even a partial relation.
		bits = sieve[dconf->reports[report_num]];
		bits = (255 - bits) + sconf->tf_closnuf + 1;


		//take care of powers of two
		while (mpz_even_p(dconf->Qvals[report_num]))
		{
			//zShiftRight32_x(Q,Q,1);
			mpz_tdiv_q_2exp(dconf->Qvals[report_num], dconf->Qvals[report_num], 1);

#ifndef SPARSE_STORE
			dconf->fb_offsets[report_num][++smooth_num] = 1;
#endif
			bits++;
		}

		i=2;
		//explicitly trial divide by small primes which we have not
		//been sieving.  because we haven't been sieving, their progressions
		//have not been updated and thus we can't use the faster methods
		//seen below.  fortunately, there shouldn't be many of these to test.
		//to speed things up, use multiplication by inverse rather than 
		//division, and do things in batches of 4 so we can use
		//the whole cache line at once (16 byte structure)


		//do the small primes in optimized batches of 4
		bound = (sconf->sieve_small_fb_start - 4);
		
		while ((uint32)i < bound)
		{
			tmp1 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> 32; 
			//at this point tmp1 is offset / prime
			tmp1 = offset - tmp1 * fullfb_ptr->prime[i];
			//now tmp1 is offset % prime

			tmp2 = offset + fullfb_ptr->correction[i+1];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> 32; 
			tmp2 = offset - tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = offset + fullfb_ptr->correction[i+2];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> 32;
			tmp3 = offset - tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = offset + fullfb_ptr->correction[i+3];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >> 32; 
			tmp4 = offset - tmp4 * fullfb_ptr->prime[i+3];
			
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			
			if (tmp1 == root1 || tmp1 == root2)
			{
				prime = fbc->prime[i];
				logp = fbc->logp[i];
                DIVIDE_ONE_PRIME(i);
			}

			root1 = fbc->root1[i+1];
			root2 = fbc->root2[i+1];
			
			if (tmp2 == root1 || tmp2 == root2)
			{
				prime = fbc->prime[i+1];
				logp = fbc->logp[i+1];
                DIVIDE_ONE_PRIME(i+1);
			}

			root1 = fbc->root1[i+2];
			root2 = fbc->root2[i+2];
			
			if (tmp3 == root1 || tmp3 == root2)
			{
				prime = fbc->prime[i+2];
				logp = fbc->logp[i+2];
                DIVIDE_ONE_PRIME(i+2);
			}

			root1 = fbc->root1[i+3];
			root2 = fbc->root2[i+3];
			
			if (tmp4 == root1 || tmp4 == root2)
			{
				prime = fbc->prime[i+3];
				logp = fbc->logp[i+3];
                DIVIDE_ONE_PRIME(i+3);
			}

            i += 4;
		}

		//finish up the rest of the small primes
		while ((uint32)i < sconf->sieve_small_fb_start)
		{
			uint64 q64;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
			
			//this is just offset % prime (but divisionless!)
			tmp = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp * (uint64)fullfb_ptr->small_inv[i];
			tmp = q64 >>  32; 
			tmp = offset - tmp * prime;

			//if offset % prime == either root, it's on the progression.  also
			//need to check for the case if root1 or root2 == prime at the same
			//time as offset mod prime = 0.  for small primes, this happens fairly
			//often.  the simple offset % prime check will miss these cases.
			if (tmp == root1 || tmp == root2)
			{
                DIVIDE_ONE_PRIME(i);
			}
			i++;
		}

		if (bits < (sconf->tf_closnuf + dconf->tf_small_cutoff))
			dconf->valid_Qs[report_num] = 0;
		else
			dconf->valid_Qs[report_num] = 1;

		dconf->smooth_num[report_num] = smooth_num;
	}

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
    TF_STG1 += my_difftime (&qs_timing_start, &qs_timing_stop);
	gettimeofday(&qs_timing_start, NULL);
#endif

	return;
}


