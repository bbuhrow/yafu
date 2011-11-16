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

#ifdef USE_YAFU_TDIV
#define DIVIDE_ONE_PRIME \
	do	\
	{	\
		dconf->fb_offsets[report_num][++smooth_num] = i;	\
		zShortDiv32(tmp32, prime, tmp32);	\
		bits += logp;	\
	} while (zShortMod32(tmp32, prime) == 0);
#else
#define DIVIDE_ONE_PRIME \
	do	\
	{	\
		dconf->fb_offsets[report_num][++smooth_num] = i;	\
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
	sieve_fb *fb;
	sieve_fb_compressed *fbc;
	tiny_fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->tinylist;
	uint8 logp, bits;
	uint32 tmp1, tmp2, tmp3, tmp4, offset, report_num;
#ifdef USE_COMPRESSED_FB
	sieve_fb_compressed *fbptr;
#endif

#ifdef DO_4X_SPV
	uint32 *offsetarray, *mask1, *mask2;
#endif

	fullfb_ptr = fullfb;
	if (parity)
	{
		fb = dconf->fb_sieve_n;
		fbc = dconf->comp_sieve_n;
	}
	else
	{
		fb = dconf->fb_sieve_p;
		fbc = dconf->comp_sieve_p;
	}

	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	for (report_num = 0; report_num < dconf->num_reports; report_num++)
	{
		uint64 q64;
#ifdef USE_YAFU_TDIV
		z32 *tmp32 = &dconf->Qvals32[report_num];
#endif

		//this one qualifies to check further, log that fact.
		dconf->num++;

		smooth_num = -1;

		//this one is close enough, compute 
		//Q(x)/a = (ax + b)^2 - N, where x is the sieve index
		//Q(x)/a = (ax + 2b)x + c;	
		offset = (bnum << BLOCKBITS) + dconf->reports[report_num];

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
			dconf->fb_offsets[report_num][++smooth_num] = 0;
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

#ifdef USE_YAFU_TDIV
		mpz_to_z32(dconf->Qvals[report_num], tmp32);

		//take care of powers of two
		while ((tmp32->val[0] & 0x1) == 0)
		{
			zShiftRight32_x(tmp32, tmp32, 1);
			dconf->fb_offsets[report_num][++smooth_num] = 1;
			bits++;
		}

#else


		//take care of powers of two
		while (mpz_even_p(dconf->Qvals[report_num]))
		{
			//zShiftRight32_x(Q,Q,1);
			mpz_tdiv_q_2exp(dconf->Qvals[report_num], dconf->Qvals[report_num], 1);
			dconf->fb_offsets[report_num][++smooth_num] = 1;
			bits++;
		}
#endif

#ifdef DO_4X_SPV
		// this won't work because to do the division via multiplication by inverse, the
		// precomputed inverses need to be 32 bits in order to accomodate values of offset greater than
		// 16 bits.  thus root1, root2, and prime are uint16's, but small_inv and correction are 
		// uint32's and the vectors don't line up.  plus, we have to do a bunch of shifting and 
		// masking in order to get around the lack of a pmulhuD instruction which limits the
		// speedup of the SIMD approach.
		//
		// potential solution: if at the end of every sieve block we advance the values of the roots
		// of these lowest primes into the next block then we can use the block offset rather than
		// the global sieve array offset.  the block offset is always less than 16 bits, so this would
		// perhaps allow the precomputed inverses to be crammed into 16 bit values.  then we can use
		// pmulhuw and do 8x small primes at once.  this is also potentially useful for medprimes up to
		// the resieving bound.
		//
		// is it worth it?  for a test c80, we spend 3.2 seconds in this stage, which is 2% of the total
		// sieving time.  if doing 8x at once makes things (conservatively) 4x faster, then we can
		// save 2.4 sec, for approx a 1.5% speedup overall.  meh.

		offsetarray = xmalloc_align(4 * sizeof(uint32));
		mask1 = xmalloc_align(4 * sizeof(uint32));
		mask2 = xmalloc_align(4 * sizeof(uint32));
		offsetarray[0] = offset; offsetarray[1] = offset; 
		offsetarray[2] = offset; offsetarray[3] = offset;
		mask1[0] = 0; mask1[1] = 0xffffffff; mask1[2] = 0; mask1[3] = 0xffffffff;
		mask2[0] = 0xffffffff; mask2[1] = 0; mask2[2] = 0xffffffff; mask2[3] = 0;


		// do i=2 and i=3
		tmp1 = offset + fullfb_ptr->correction[2];
		q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[2];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * fullfb_ptr->prime[2];

		prime = fbc->prime[2];
		root1 = fbc->root1[2];
		root2 = fbc->root2[2];
		logp = fbc->logp[2];

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				dconf->fb_offsets[report_num][++smooth_num] = 2;
				zShortDiv32(Q,prime,Q);
				bits += logp;
			} while (zShortMod32(Q,prime) == 0);
		}

		tmp1 = offset + fullfb_ptr->correction[3];
		q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[3];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * fullfb_ptr->prime[3];

		prime = fbc->prime[3];
		root1 = fbc->root1[3];
		root2 = fbc->root2[3];
		logp = fbc->logp[3];

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				dconf->fb_offsets[report_num][++smooth_num] = 3;
				zShortDiv32(Q,prime,Q);
				bits += logp;
			} while (zShortMod32(Q,prime) == 0);
		}

		// start i=4 in batches of 4 (aligned memory boundary)
		i=4;
		bound = (sconf->sieve_small_fb_start - 4);		
		while ((uint32)i < bound)
		{
			uint32 result;

			ASM_G (
				"movdqa	(%1), %%xmm0 \n\t" /* offset into xmm0 */
				"movdqa (%2), %%xmm1 \n\t" /* correction into xmm1 */
				"movdqa (%3), %%xmm2 \n\t" /* inv into xmm2 */
				"movdqa (%4), %%xmm3 \n\t" /* primes into xmm3 */
				"movdqa (%5), %%xmm4 \n\t" /* root1 into xmm4 */
				"movdqa (%6), %%xmm5 \n\t" /* root2 into xmm5 */
				"movdqa (%7), %%xmm8 \n\t" /* mask1 into xmm8 */
				"movdqa (%8), %%xmm9 \n\t" /* mask2 into xmm9 */
				"paddd %%xmm0, %%xmm1 \n\t" /* add offset with correction, overwrite correction */
				"movdqa %%xmm2, %%xmm6 \n\t" /* copy of inv */
				"psrldq $4, %%xmm6 \n\t" /* shift right 4 bytes */
				"pmuludq %%xmm1, %%xmm2 \n\t" /* tmp * inv3,1 */
				"pmuludq %%xmm1, %%xmm6 \n\t" /* tmp * inv4,2 */
				"psrldq $4, %%xmm2 \n\t" /* shift right 4 bytes */
				"psrldq $4, %%xmm6 \n\t" /* shift right 4 bytes */
				"movdqa %%xmm3, %%xmm7 \n\t" /* copy of primes */
				"psrldq $4, %%xmm7 \n\t" /* shift right 4 bytes */
				"pmuludq %%xmm3, %%xmm2 \n\t" /* primes3,1 * tmp3,1 */
				"pmuludq %%xmm7, %%xmm6 \n\t" /* primes4,2 * tmp4,2 */
				"psrldq $4, %%xmm2 \n\t" /* shift right 4 bytes */
				"pand %%xmm8, %%xmm6 \n\t" /* zero dwords 3 and 1 of xmm6 */
				"pand %%xmm9, %%xmm2 \n\t" /* zero dwords 4 and 2 of xmm2 */
				"por %%xmm2, %%xmm6 \n\t" /* recombine high dwords of products */
				"psubd %%xmm6, %%xmm0 \n\t" /* final subtract to accomplish the mod operation */
				"pcmpeqd %%xmm0, %%xmm4 \n\t" /* test equality to root1 */
				"pcmpeqd %%xmm0, %%xmm5 \n\t" /* test equality to root2 */
				"por %%xmm4, %%xmm5 \n\t" /* combine tests */
				"pmovmskb %%xmm5, %0 \n\t" /* result to output variable */
				: "=r"(result)
				: "r"(offsetarray), "r"(fullfb_ptr->correction + i), "r"(fullfb_ptr->small_inv + i), "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(mask1), "r"(mask2)
				: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9"); 

			if (result & 0x8)
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i;
					zShortDiv32(Q,fbc->prime[i],Q);
					bits += fbc->logp[i];
				} while (zShortMod32(Q,fbc->prime[i]) == 0);
			}
			if (result & 0x80)
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i+1;
					zShortDiv32(Q,fbc->prime[i+1],Q);
					bits += fbc->logp[i+1];
				} while (zShortMod32(Q,fbc->prime[i+1]) == 0);
			}
			if (result & 0x800)
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i+2;
					zShortDiv32(Q,fbc->prime[i+2],Q);
					bits += fbc->logp[i+2];
				} while (zShortMod32(Q,fbc->prime[i+2]) == 0);
			}
			if (result & 0x8000)
			{
				do
				{
					dconf->fb_offsets[report_num][++smooth_num] = i+3;
					zShortDiv32(Q,fbc->prime[i+3],Q);
					bits += fbc->logp[i+3];
				} while (zShortMod32(Q,fbc->prime[i+3]) == 0);
			}

			i += 4;

		}

		free(offsetarray);
		free(mask1);
		free(mask2);

#else

		i=2;
		//explicitly trial divide by small primes which we have not
		//been sieving.  because we haven't been sieving, their progressions
		//have not been updated and thus we can't use the faster methods
		//seen below.  fortunately, there shouldn't be many of these to test
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
			i++;

			tmp2 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i];
			tmp2 = q64 >> 32; 
			tmp2 = offset - tmp2 * fullfb_ptr->prime[i];
			i++;

			tmp3 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i];
			tmp3 = q64 >> 32;
			tmp3 = offset - tmp3 * fullfb_ptr->prime[i];
			i++;

			tmp4 = offset + fullfb_ptr->correction[i];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i];
			tmp4 = q64 >> 32; 
			tmp4 = offset - tmp4 * fullfb_ptr->prime[i];
			
			i -= 3;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if (tmp1 == root1 || tmp1 == root2)
			{
				DIVIDE_ONE_PRIME
			}

			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if (tmp2 == root1 || tmp2 == root2)
			{
				DIVIDE_ONE_PRIME
			}

			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if (tmp3 == root1 || tmp3 == root2)
			{
				DIVIDE_ONE_PRIME
			}

			i++;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif

			if (tmp4 == root1 || tmp4 == root2)
			{
				DIVIDE_ONE_PRIME
			}
			i++;
		}
		
#endif		

		//finish up the rest of the small primes
		while ((uint32)i < sconf->sieve_small_fb_start)
		{
			uint64 q64;

#ifdef USE_COMPRESSED_FB
			fbptr = fbc + i;
			prime = fbptr->prime_and_logp & 0xFFFF;
			root1 = fbptr->roots & 0xFFFF;
			root2 = fbptr->roots >> 16;
			logp = fbptr->prime_and_logp >> 16;
#else
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];
			logp = fbc->logp[i];
#endif
			
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
				DIVIDE_ONE_PRIME
			}
			i++;
		}

		if (bits < (sconf->tf_closnuf + sconf->tf_small_cutoff))
			dconf->valid_Qs[report_num] = 0;
		else
			dconf->valid_Qs[report_num] = 1;

		dconf->smooth_num[report_num] = smooth_num;
	}

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG1 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);

	gettimeofday(&qs_timing_start, NULL);
#endif

	return;
}


