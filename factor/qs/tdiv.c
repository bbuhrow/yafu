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

//#define SIQSDEBUG 1

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

this file contains code implementing 6) as well as other auxiliary routines


*/



void trial_divide_Q_siqs(uint32 report_num,  uint8 parity, 
						 uint32 poly_id, uint32 bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	uint64 q64, f64;
	int j,it;
	uint32 prime;
	int smooth_num;
	uint32 *fb_offsets;
	uint32 polya_factors[20];
	sieve_fb *fb;
	uint32 offset, block_loc;

#if !defined(TDIV_GMP)
	z32 *Q = &dconf->Qvals[report_num];
#endif
	fb_offsets = &dconf->fb_offsets[report_num][0];
	smooth_num = dconf->smooth_num[report_num];
	block_loc = dconf->reports[report_num];
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	offset = (bnum << BLOCKBITS) + block_loc;

	if (parity)
		fb = dconf->fb_sieve_n;
	else
		fb = dconf->fb_sieve_p;

	//check for additional factors of the a-poly factors
	//make a separate list then merge it with fb_offsets
	it=0;	//max 20 factors allocated for - should be overkill
	for (j = 0; (j < dconf->curr_poly->s) && (it < 20); j++)
	{
		//fbptr = fb + dconf->curr_poly->qlisort[j];
		//prime = fbptr->prime;
		prime = fb[dconf->curr_poly->qlisort[j]].prime;

#if defined(TDIV_GMP)
		while ((mpz_tdiv_ui(dconf->Qvals[report_num],prime) == 0) && (it < 20))
		{
			mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);
			polya_factors[it++] = dconf->curr_poly->qlisort[j];
		}
#else
		while ((zShortMod32(Q, prime) == 0) && (it < 20))
		{
			zShortDiv32(Q, prime, Q);
			polya_factors[it++] = dconf->curr_poly->qlisort[j];
		}

#endif
	}

	//check if it completely factored by looking at the unfactored portion in tmp
#if defined(TDIV_GMP)
	if ((mpz_size(dconf->Qvals[report_num]) == 1) && 
		((uint32)mpz_get_ui(dconf->Qvals[report_num]) < sconf->large_prime_max))
#else
	if ((Q->size == 1) && 
		(Q->val[0] < sconf->large_prime_max))
#endif
	{
		uint32 large_prime[2];
		
#if defined(TDIV_GMP)
		large_prime[0] = (uint32)mpz_get_ui(dconf->Qvals[report_num]); //Q->val[0];
#else
		large_prime[0] = Q->val[0];
#endif
		large_prime[1] = 1;

		//add this one
		buffer_relation(offset,large_prime,smooth_num+1,
			fb_offsets,poly_id,parity,dconf,polya_factors,it);

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG6 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif

		return;
	}

	if (sconf->use_dlp == 0)
		return;

#if defined(TDIV_GMP)
	//quick check if Q is way too big for DLP (more than 64 bits)
	
	if (mpz_sizeinbase(dconf->Qvals[report_num], 2) >= 64)
		return;

	q64 = mpz_get_64(dconf->Qvals[report_num]);
#else
	if (Q->size >= 3)
		return;

	q64 = (uint64)Q->val[0] + ((uint64)Q->val[1] << 32);
#endif


	if ((q64 > sconf->max_fb2) && (q64 < sconf->large_prime_max2))
	{	
		//quick prime check: compute 2^(residue-1) mod residue.  
		uint64 res;

#if BITS_PER_DIGIT == 32
		sp642z(q64,&dconf->qstmp1);
		zShortSub(&dconf->qstmp1,1,&dconf->qstmp2);
		zModExp(&zTwo,&dconf->qstmp2,&dconf->qstmp1,&dconf->qstmp3);
		res = z264(&dconf->qstmp3);
#else
		spModExp(2, q64 - 1, q64, &res);
#endif

		//if equal to 1, assume it is prime.  this may be wrong sometimes, but we don't care.
		//more important to quickly weed out probable primes than to spend more time to be
		//more sure.
		if (res == 1)
		{
#ifdef QS_TIMING
			gettimeofday (&qs_timing_stop, NULL);
			qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

			TF_STG6 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
			free(qs_timing_diff);
#endif
			dconf->dlp_prp++;
			return;
		}
		
		//try to find a double large prime
		dconf->attempted_squfof++;
		sp642z(q64, &dconf->qstmp1);
		f64 = sp_shanks_loop(&dconf->qstmp1,sconf->obj);
		if (f64 > 1 && f64 != q64)
		{
			uint32 large_prime[2];

			large_prime[0] = (uint32)f64;
			large_prime[1] = (uint32)(q64 / f64);

			if (large_prime[0] < sconf->large_prime_max 
				&& large_prime[1] < sconf->large_prime_max)
			{
				//add this one
				dconf->dlp_useful++;
				buffer_relation(offset,large_prime,smooth_num+1,
					fb_offsets,poly_id,parity,dconf,polya_factors,it);
			}
				
		}
		else
			dconf->failed_squfof++;
	}
	else
		dconf->dlp_outside_range++;

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG6 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif
	
	return;
}

void buffer_relation(uint32 offset, uint32 *large_prime, uint32 num_factors, 
						  uint32 *fb_offsets, uint32 poly_id, uint32 parity,
						  dynamic_conf_t *conf, uint32 *polya_factors, 
						  uint32 num_polya_factors)
{
	//put this relations's info into a temporary buffer which
	//will get merged with other such buffers (if multi-threaded) and
	//dumped to file once the threads are joined.
	siqs_r *rel;
	uint32 i, j, k;

	//first check that this relation won't overflow the buffer
	if (conf->buffered_rels >= conf->buffered_rel_alloc)
	{
		printf("reallocating relation buffer\n");
		conf->relation_buf = (siqs_r *)realloc(conf->relation_buf, 
			conf->buffered_rel_alloc * 2 * sizeof(siqs_r));
		if (conf->relation_buf == NULL)
		{
			printf("error re-allocating temporary storage of relations\n");
			exit(-1);
		}
		conf->buffered_rel_alloc *= 2;
	}

	//then stick all the info in the buffer
	rel = conf->relation_buf + conf->buffered_rels;
	
	rel->sieve_offset = offset;
	rel->parity = parity;
	rel->poly_idx = poly_id;

	rel->fb_offsets = (uint32 *)malloc(
		(num_polya_factors + num_factors) * sizeof(uint32));

	//merge in extra factors of the apoly factors
	i = j = k = 0;
	while (k < num_factors && j < num_polya_factors) {
		if (fb_offsets[k] < polya_factors[j]) {
			rel->fb_offsets[i++] = fb_offsets[k++];
		}
		else if (fb_offsets[k] > polya_factors[j]) {
			rel->fb_offsets[i++] = polya_factors[j++];
		}
		else {
			rel->fb_offsets[i++] = fb_offsets[k++];
			rel->fb_offsets[i++] = polya_factors[j++];
		}
	}
	while (k < num_factors)
		rel->fb_offsets[i++] = fb_offsets[k++];
	while (j < num_polya_factors)
		rel->fb_offsets[i++] = polya_factors[j++];
		
	rel->num_factors = num_factors + num_polya_factors;
	rel->large_prime[0] = large_prime[0];
	rel->large_prime[1] = large_prime[1];

	conf->buffered_rels++;
	return;
}

void save_relation_siqs(uint32 offset, uint32 *large_prime, uint32 num_factors, 
						  uint32 *fb_offsets, uint32 poly_id, uint32 parity,
						  static_conf_t *conf)
{
	char buf[LINE_BUF_SIZE];
	fact_obj_t *obj = conf->obj;
	uint32 i, k;

	//store to file
	i = sprintf(buf, "R ");

	if (parity)
		i += sprintf(buf + i, "-%x ", offset);
	else
		i += sprintf(buf + i, "%x ", offset);
	
	i += sprintf(buf + i, "%x ", poly_id);
	
	k = 0;
	while (k < num_factors)
		i += sprintf(buf + i, "%x ", fb_offsets[k++]);

	if (large_prime[0] < large_prime[1])
		i += sprintf(buf + i, "L %x %x\n", large_prime[0], large_prime[1]);
	else
		i += sprintf(buf + i, "L %x %x\n", large_prime[1], large_prime[0]);

	qs_savefile_write_line(&obj->qs_obj.savefile, buf);

	/* for partial relations, also update the bookeeping for
	   tracking the number of fundamental cycles */

	if (large_prime[0] != large_prime[1]) {
		yafu_add_to_cycles(conf, obj->flags, large_prime[0], large_prime[1]);
		conf->num_cycles++;
	}
	else {
		conf->num_relations++;
	}

	return;
}

