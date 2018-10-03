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
    sieve_fb_compressed *fbc = dconf->comp_sieve_p;
	uint32 offset, block_loc;
	fb_offsets = &dconf->fb_offsets[report_num][0];
	smooth_num = dconf->smooth_num[report_num];
	block_loc = dconf->reports[report_num];
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	offset = (bnum << sconf->qs_blockbits) + block_loc;

#ifdef USE_YAFU_TDIV
	z32_to_mpz(&dconf->Qvals32[report_num], dconf->Qvals[report_num]);
#endif

	//check for additional factors of the a-poly factors
	//make a separate list then merge it with fb_offsets
	it=0;	//max 20 factors allocated for - should be overkill
	for (j = 0; (j < dconf->curr_poly->s) && (it < 20); j++)
	{
        prime = fbc->prime[dconf->curr_poly->qlisort[j]]; // .prime;

		while ((mpz_tdiv_ui(dconf->Qvals[report_num],prime) == 0) && (it < 20))
		{
			mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);

#ifndef SPARSE_STORE
			polya_factors[it++] = dconf->curr_poly->qlisort[j];
#endif
		}
	}

	// check if it completely factored by looking at the unfactored portion in tmp
	if ((mpz_size(dconf->Qvals[report_num]) == 1) && 
		(mpz_cmp_ui(dconf->Qvals[report_num], sconf->large_prime_max) < 0))
	{
		uint32 large_prime[2];
		
		large_prime[0] = (uint32)mpz_get_ui(dconf->Qvals[report_num]); //Q->val[0];
		large_prime[1] = 1;

		//add this one
		if (sconf->is_tiny)
		{	
			// we need to encode both the a_poly and b_poly index
			// in poly_id
			poly_id |= (sconf->total_poly_a << 16);
			buffer_relation(offset,large_prime,smooth_num+1,
				fb_offsets,poly_id,parity,dconf,polya_factors,it, 1);
		}
        else
        {
            buffer_relation(offset, large_prime, smooth_num + 1,
                fb_offsets, poly_id, parity, dconf, polya_factors, it, 1);
        }

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
        TF_STG6 +=  my_difftime (&qs_timing_start, &qs_timing_stop);
#endif

		return;
	}

	if (sconf->use_dlp == 0)
		return;

	//quick check if Q is way too big for DLP (more than 64 bits)	
    if (mpz_sizeinbase(dconf->Qvals[report_num], 2) >= 64)
    {
        dconf->dlp_outside_range++;
        return;
    }

	q64 = mpz_get_64(dconf->Qvals[report_num]);

	if ((q64 > sconf->max_fb2) && (q64 < sconf->large_prime_max2))
	{	
		//quick prime check: compute 2^(residue-1) mod residue.  
		uint64 res;

#if BITS_PER_DIGIT == 32
		mpz_set_64(dconf->gmptmp1, q64);
		mpz_set_64(dconf->gmptmp2, 2);
		mpz_set_64(dconf->gmptmp3, q64-1);

		mpz_powm(dconf->gmptmp1, dconf->gmptmp2, dconf->gmptmp3, dconf->gmptmp1);
		res = mpz_get_64(dconf->gmptmp1);
#elif defined (FORCE_GENERIC) && !defined(TARGET_KNC)
        mpz_set_64(dconf->gmptmp1, q64);
        mpz_set_64(dconf->gmptmp2, 2);
        mpz_set_64(dconf->gmptmp3, q64-1);

        mpz_powm(dconf->gmptmp1, dconf->gmptmp2, dconf->gmptmp3, dconf->gmptmp1);
        res = mpz_get_64(dconf->gmptmp1);
#else
		// meh, not really any faster, but fun to write...
		res = spPRP2(q64);
#endif

		//if equal to 1, assume it is prime.  this may be wrong sometimes, but we don't care.
		//more important to quickly weed out probable primes than to spend more time to be
		//more sure.
		if (res == 1)
		{
#ifdef QS_TIMING
			gettimeofday (&qs_timing_stop, NULL);
            TF_STG6 += my_difftime(&qs_timing_start, &qs_timing_stop);
#endif
			dconf->dlp_prp++;
			return;
		}
		
		//try to find a double large prime
#ifdef USE_VEC_SQUFOF

        buffer_relation(offset, NULL, smooth_num + 1,
            fb_offsets, poly_id, parity, dconf, polya_factors, it, q64);

#else
		dconf->attempted_squfof++;
		
#ifdef _MSC_VER
        mpz_set_64(dconf->gmptmp1, q64);
		f64 = sp_shanks_loop(dconf->gmptmp1, sconf->obj);
#else
        f64 = spbrent(q64, 1, 1024);
#endif
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
                buffer_relation(offset, large_prime, smooth_num + 1,
                    fb_offsets, poly_id, parity, dconf, polya_factors, it, 1);
            }
				
		}
		else
		{
			dconf->failed_squfof++;
		}

#endif


	}
    else
    {
        //gmp_printf("residue %Zd (%" PRIu64 ") is outside range\n", dconf->Qvals[report_num], q64);
        dconf->dlp_outside_range++;
    }

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
    TF_STG6 += my_difftime(&qs_timing_start, &qs_timing_stop);
#endif
	
	return;
}

void buffer_relation(uint32 offset, uint32 *large_prime, uint32 num_factors,
    uint32 *fb_offsets, uint32 poly_id, uint32 parity,
    dynamic_conf_t *conf, uint32 *polya_factors,
    uint32 num_polya_factors, uint64 unfactored_residue)
{
    // put this relations's info into a temporary buffer which
    // will get merged with other such buffers (if multi-threaded) and
    // dumped to file once the threads are joined.
    siqs_r *rel;
    uint32 i, j, k;

    // first check that this relation won't overflow the buffer
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

    // then stick all the info in the buffer
    rel = conf->relation_buf + conf->buffered_rels;

    rel->sieve_offset = offset;
    rel->parity = parity;
    rel->poly_idx = poly_id;

    if ((num_polya_factors + num_factors) > MAX_SMOOTH_PRIMES)
    {
        printf("error: too many smooth primes!\n");
        exit(234);
    }

#ifdef SPARSE_STORE
    // extra factors of polya are not added to the list of factors. they will be added
    // during filtering on only the relations that survive singleton removal.
    memcpy(rel->fb_offsets, fb_offsets, sizeof(uint32) * num_factors);

#else
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
#endif

    rel->num_factors = num_factors + num_polya_factors;

#ifdef USE_VEC_SQUFOF

    if (unfactored_residue > 1)
    {
        if (conf->num_64bit_residue >= 4096)
        {
            printf("uh-oh, too many residues\n"); fflush(stdout);
        }

        //printf("adding %lu to unfactored residue list in position %d, relation %d\n", 
        //    unfactored_residue, conf->num_64bit_residue, conf->buffered_rels);

        conf->unfactored_residue[conf->num_64bit_residue] = unfactored_residue;        

        // use this to signify we need to add factors later.
        rel->large_prime[0] = 0xffffffff;
        conf->num_64bit_residue++;
    }
    else
    {
        rel->large_prime[0] = large_prime[0];
        rel->large_prime[1] = large_prime[1];
    }

#else
    
    rel->large_prime[0] = large_prime[0];
	rel->large_prime[1] = large_prime[1];

#endif


	conf->buffered_rels++;
	return;
}

void save_relation_siqs(uint32 offset, uint32 *large_prime, uint32 num_factors, 
						  uint32 *fb_offsets, uint32 poly_id, uint32 parity,
						  static_conf_t *conf)
{
	char buf[1024];
	fact_obj_t *obj = conf->obj;
	uint32 i, k;

	if (conf->in_mem)
	{
		// copy info to sconf relation structure
		siqs_r *r;

		//first check that this relation won't overflow the buffer
		if (conf->buffered_rels >= conf->buffered_rel_alloc)
		{
			printf("reallocating in-mem relation storage...\n");
			conf->in_mem_relations = (siqs_r *)realloc(conf->in_mem_relations, 
				conf->buffered_rel_alloc * 2 * sizeof(siqs_r));
			if (conf->in_mem_relations == NULL)
			{
				printf("error re-allocating in-memory storage of relations\n");
				exit(-1);
			}
			conf->buffered_rel_alloc *= 2;
		}

		r = conf->in_mem_relations + conf->buffered_rels++;

		r->large_prime[0] = large_prime[0];
		r->large_prime[1] = large_prime[1];
		r->num_factors = num_factors;
		r->poly_idx = poly_id;
		r->parity = parity;
		r->sieve_offset = offset;
		for (i=0; i<num_factors; i++)
			r->fb_offsets[i] = fb_offsets[i];

	}
	else
	{

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
	}

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

