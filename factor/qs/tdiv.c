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

#include "qs.h"
#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "cofactorize.h"
#ifdef USE_BATCH_FACTOR
#include "batch_factor.h"
#endif
#include "microecm.h"
#include "tinyecm.h"

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


int split_3lp_tdiv(mpz_t candidate3lp, mpz_t _1, mpz_t _2, mpz_t _3,
    uint32_t max_factor, uint64_t* lcg)
{

    if (getfactor_tecm(candidate3lp, _1, 0, lcg) > 0)
    {
        if ((mpz_sizeinbase(_1, 2) <= 32) && (mpz_get_ui(_1) < max_factor))
        {
            mpz_tdiv_q(_2, candidate3lp, _1);

            // if the remaining residue is obviously too big, we're done.
            if (mpz_sizeinbase(_2, 2) > 64) //((max_primebits[s] * 2)))
            {
                return 0;
            }

            // check if the residue is prime.  could again use
            // a cheaper method.
            if (mpz_probab_prime_p(_2, 1) > 0)
            {
                if ((mpz_sizeinbase(_2, 2) <= 32) && (mpz_get_ui(_2) < max_factor))
                    return 2;
                else
                    return 0;
            }

            // ok, so we have extracted one suitable factor, and the 
            // cofactor is not prime and a suitable size.  Do more work to 
            // split the cofactor.
            // todo: target this better based on expected factor size.
            uint64_t q64;
            uint64_t f64;

            q64 = mpz_get_ui(_2);
            f64 = getfactor_uecm(q64, 0, lcg);
            mpz_set_ui(_3, f64);

            if (f64 > 1)
            {
                mpz_tdiv_q_ui(_2, _2, f64);

                if (mpz_get_ui(_2) > max_factor) {
                    return 0;
                }
                if (mpz_get_ui(_3) > max_factor) {
                    return 0;
                }
                if (mpz_probab_prime_p(_1, 1) == 0)
                {
                    return 0;
                }
                if (mpz_probab_prime_p(_2, 1) == 0)
                {
                    return 0;
                }
                if (mpz_probab_prime_p(_3, 1) == 0)
                {
                    return 0;
                }

                return 3;
            }
        }
        else
        {
            // check if the factor is prime.  could again use
            // a cheaper method.
            if (mpz_probab_prime_p(_1, 1) > 0)
            {
                // if the factor is obviously too big, give up.  This isn't a
                // failure since we haven't expended much effort yet.
                return 0;
            }
            else
            {
                // tecm found a composite first factor.
                // if it is obviously too big, we're done.
                if (mpz_sizeinbase(_1, 2) > 64) //((max_primebits[s] * 2)))
                {
                    return 0;
                }

                // isolate the 2nd smaller factor, and check its size.
                mpz_tdiv_q(_2, candidate3lp, _1);

                if (mpz_sizeinbase(_2, 2) > 32) //(max_primebits[s]))
                {
                    return 0;
                }

                // todo: target this better based on expected factor size.
                uint64_t q64;
                uint64_t f64;

                q64 = mpz_get_ui(_1);
                f64 = getfactor_uecm(q64, 0, lcg);
                mpz_set_ui(_3, f64);


                if (f64 > 1)
                {
                    mpz_tdiv_q_ui(_1, _1, f64);

                    if (mpz_get_ui(_1) > max_factor) {
                        return 0;
                    }
                    if (mpz_get_ui(_3) > max_factor) {
                        return 0;
                    }
                    if (mpz_probab_prime_p(_1, 1) == 0)
                    {
                        return 0;
                    }
                    if (mpz_probab_prime_p(_2, 1) == 0)
                    {
                        return 0;
                    }
                    if (mpz_probab_prime_p(_3, 1) == 0)
                    {
                        return 0;
                    }
                    return 3;
                }
                else
                {
                    return 0;
                }
            }
        }
    }
    else
    {
        // if ecm can't find a factor, give up.  
        // unless this is a DLP with lpbr/a > 32... i.e., if the
        // large factor size is greater than 64 bits but less than
        // lpbr/a * 2.  In that case run mpqs... or tecm with
        // greater effort.

        return 0;
    }

    return 0;
}


void trial_divide_Q_siqs(uint32_t report_num,  uint8_t parity, 
						 uint32_t poly_id, uint32_t bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	uint64_t q64, f64;
	int j, i, it;
	uint32_t prime;
	int dlp_done = 0;
	int smooth_num;
	uint32_t *fb_offsets;
	uint32_t polya_factors[20];
    sieve_fb_compressed *fbc = dconf->comp_sieve_p;
	uint32_t offset, block_loc;
	fb_offsets = &dconf->fb_offsets[report_num][0];
	smooth_num = dconf->smooth_num[report_num];
	block_loc = dconf->reports[report_num];
    int VFLAG = sconf->obj->VFLAG;
    int THREADS = sconf->obj->THREADS;

	offset = (bnum << sconf->qs_blockbits) + block_loc;

    //gmp_printf("commencing tdiv on %Zd\n", dconf->Qvals[report_num]);

	// check for additional factors of the a-poly factors
	// make a separate list then merge it with fb_offsets
	it=0;	// max 20 factors allocated for - should be overkill
	for (j = 0; (j < dconf->curr_poly->s - dconf->num_alp) && (it < 20); j++)
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

#if 0 //def SPARSE_STORE
    // edit: removed because med_B seems to be a little too high
    // of a bound for non-stored factors.  Puts quite a bit of 
    // extra work in the filtering code.  Current bound is the
    // small prime bound.

    // only actually store indices of primes larger than med_B.
    // If the relation ultimately proves useful during filtering
    // then we will trial divide it there to recover its divisors,
    // including any factors of the A poly.
    for (j = 0, i = 0; j <= smooth_num; j++)
    {
        if (fb_offsets[j] >= sconf->factor_base->med_B)
        {
            fb_offsets[i++] = fb_offsets[j];
        }
    }
    smooth_num = i - 1;
#endif
    

	// check if it completely factored by looking at the unfactored portion in tmp
	if ((mpz_size(dconf->Qvals[report_num]) == 1) && 
		(mpz_cmp_ui(dconf->Qvals[report_num], sconf->large_prime_max) < 0))
	{
        // save this slp (single large prime)
		uint32_t large_prime[4];
		
		large_prime[0] = (uint32_t)mpz_get_ui(dconf->Qvals[report_num]); //Q->val[0];
		//if (dconf->num_alp == 1)
        //    large_prime[1] = dconf->curr_poly->qlisort[dconf->curr_poly->s - 1];
        //else
            large_prime[1] = 1;
		large_prime[2] = 1;
        large_prime[3] = 1;

#ifdef GATHER_RESIDUE_STATS

        // as part of analyzing 3lp parameterizations, save off this residue.  
        // These will be fully factored later and sorted
        // so that a synthetic data set resembling this real one can
        // be constructed for deeper analysis of cycle generation for
        // this parameterization.
        fprintf(sconf->residue_files[dconf->tid],
            "%u\n", large_prime[0]);

#endif

        if (dconf->num_alp == 1)
        {
            if (large_prime[0] == 1)
                dconf->num_slp++;
            else
                dconf->dlp_useful++;
        }
        else
        {
            if (large_prime[0] == 1)
                dconf->num_full++;
            else
                dconf->num_slp++;
        }

		// add this one
        buffer_relation(offset, large_prime, smooth_num + 1,
            fb_offsets, dconf->curr_poly->index, poly_id, 
            parity, dconf, polya_factors, it, 1);

		return;
	}

    // if we are not considering dlps then we are done.
	if (sconf->use_dlp == 0)
		return;

	// quick check if Q is way too big for DLP (more than 64 bits)	
	if (mpz_sizeinbase(dconf->Qvals[report_num], 2) < 64)
	{
		uint64_t res;

        // not obviously too big.  So now check if we are actually
        // within the defined bounds for DLPs.
		q64 = mpz_get_64(dconf->Qvals[report_num]);

		if ((q64 > sconf->max_fb2) && (q64 < sconf->large_prime_max2))
		{
            // processing DLP's only with the batch GCD is slower...
            if ((sconf->use_dlp >= 2) && (dconf->do_batch) && 0)
            {
                int32_t soffset = offset;

                if (parity)
                {
                    soffset *= -1;
                }

                mpz_set_ui(dconf->gmptmp1, 0);

                // use this field to record how many we've batched.
                dconf->attempted_cosiqs++;
                relation_batch_add(dconf->curr_poly->index, poly_id, soffset, fb_offsets, smooth_num + 1,
                    dconf->Qvals[report_num], NULL, 0, dconf->gmptmp1, NULL, &dconf->rb);

                // start processing the relations in the indicated buffer if running
                // single threaded.  otherwise, let the threading dispatcher assign 
                // the processing to a free thread.
                if ((dconf->batch_run_override > 0) ||
                    ((THREADS == 1) && (sconf->rb[0].num_relations > sconf->rb[0].target_relations)))
                {
                    struct timeval start;
                    struct timeval stop;
                    double ttime;
                    relation_batch_t* rb;
                    int i;

                    if (THREADS == 1)
                        rb = &sconf->rb[0];
                    else
                        rb = &sconf->rb[dconf->batch_run_override - 1];
                    
                    if (VFLAG > 1)
                    {
                        printf("now processing %u relations in batch %d in thread %d\n",
                            rb->num_relations, dconf->batch_run_override, dconf->tid);
                    }

                    gettimeofday(&start, NULL);
                    relation_batch_run(rb, &dconf->lcg_state);
                    gettimeofday(&stop, NULL);

                    ttime = ytools_difftime(&start, &stop);
                    if (VFLAG > 1)
                    {
                        printf("relation_batch_run took %1.4f sec producing %u tlp's\n",
                            ttime, rb->num_success);
                    }

                    rb->conversion_ratio =
                        (double)rb->num_success / (double)rb->num_relations;

                    // take our new tlp relations and buffer them to be
                    // saved out to the data file.
                    for (i = 0; i < rb->num_relations; i++)
                    {
                        cofactor_t* c = rb->relations + i;
                        uint32_t* f = rb->factors + c->factor_list_word;

                        if (c->success)
                        {
                            uint8_t parity = c->signed_offset < 0 ? 1 : 0;

                            if (c->success == 4)
                                dconf->qlp_useful++;
                            else if (c->success == 3)
                                dconf->tlp_useful++;
                            else if (c->success == 2)
                                dconf->dlp_useful++;
                            else if (c->success == 1)
                                dconf->num_slp++;

                            buffer_relation(abs(c->signed_offset), c->lp_r, c->num_factors_r,
                                f, c->a, c->b, parity, dconf, NULL, 0, 1);
                        }
                    }

                    if (VFLAG > 1)
                    {
                        printf("done processing batch %d in thread %d, found %u relations\n",
                            dconf->batch_run_override, dconf->tid, rb->num_success);
                    }

                    // clear the relation batch now that we are done processing it.
                    rb->num_relations = 0;
                    rb->num_success = 0;
                    rb->num_factors = 0;

                    // signal we are done processing the batch.
                    if (THREADS > 1)
                        dconf->batch_run_override = -1;
                }

                // if batch factoring, we're done now.
                return;
            }


#ifdef GATHER_RESIDUE_STATS

            // as part of analyzing 3lp parameterizations, save off this residue.  
            // These will be fully factored later and sorted
            // so that a synthetic data set resembling this real one can
            // be constructed for deeper analysis of cycle generation for
            // this parameterization.
            fprintf(sconf->residue_files[dconf->tid], "%lu\n", q64);

#endif

			//quick prime check: compute 2^(residue-1) mod residue.  

#if defined(_MSC_VER) || (BITS_PER_DIGIT == 32)
			mpz_set_64(dconf->gmptmp1, q64);
			mpz_set_64(dconf->gmptmp2, 2);
			//mpz_set_64(dconf->gmptmp3, q64 - 1);

			mpz_powm_ui(dconf->gmptmp1, dconf->gmptmp2, q64 - 1, dconf->gmptmp1);
			res = mpz_get_64(dconf->gmptmp1);
#elif defined (FORCE_GENERIC)
			mpz_set_64(dconf->gmptmp1, q64);
			mpz_set_64(dconf->gmptmp2, 2);
			mpz_set_64(dconf->gmptmp3, q64 - 1);

			mpz_powm(dconf->gmptmp1, dconf->gmptmp2, dconf->gmptmp3, dconf->gmptmp1);
			res = mpz_get_64(dconf->gmptmp1);
#else
			// meh, not really any faster, but fun to write...
			res = spPRP2(q64);
            //res = pow2m(q64 - 1, q64);
#endif

			// if equal to 1, assume it is prime.  this may be wrong sometimes, but we don't care.
			// more important to quickly weed out probable primes than to spend more time to be
			// more sure.
			if (res == 1)
			{
				dconf->dlp_prp++;
				return;
			}

			// try to find a double large prime.
            // now with superfast ecm, squfof, rho, etc are obsolete.
            // but we still for now use "attempted_squfof" to indicate
            // that we are attempting to factor the potential dlp 
            // residue.
			dconf->attempted_squfof++;
            f64 = getfactor_uecm(q64, 0, &dconf->lcg_state); // do_uecm(q64);
			
			if ((f64 > 1) && (f64 != q64) && ((q64 % f64) == 0))
			{
				uint32_t large_prime[4];

                if ((f64 < (uint64_t)sconf->large_prime_max) && 
                    ((q64 / f64) < (uint64_t)sconf->large_prime_max))
				{
					//add this one
                    large_prime[0] = (uint32_t)f64;
                    large_prime[1] = (uint32_t)(q64 / f64);
                    if (dconf->num_alp == 1)
                    {
                        large_prime[2] = dconf->curr_poly->qlisort[dconf->curr_poly->s - 1];
                        dconf->tlp_useful++;
                    }
                    else
                    {
                        large_prime[2] = 1;
                        dconf->dlp_useful++;
                    }
                    large_prime[3] = 1;

					
					buffer_relation(offset, large_prime, smooth_num + 1,
						fb_offsets, dconf->curr_poly->index, poly_id, parity,
                        dconf, polya_factors, it, 1);
				}
			}
			else
			{
				dconf->failed_squfof++;
			}

			// whether we found a DLP or not, we are done checking.
			return;
		}
		else
		{
			dconf->dlp_outside_range++;
			
			// too big for DLP, but too small for TLP (if active).
			return;
		}
	}
	else
    {
        dconf->dlp_outside_range++;
    }

    // if we are not considering TLP, then we are done.
	if (sconf->use_dlp < 2)
		return;
    
	// quick check if Q is obviously too big.
	if (mpz_sizeinbase(dconf->Qvals[report_num], 2) < 96)
	{
		double qfloat = mpz_get_d(dconf->Qvals[report_num]);

//#define OUTPUT_TLP_ATTEMPT_DETAILS

        // not obviously too big, see if it is actually within
        // the defined tlp bounds.
		if ((qfloat > sconf->max_fb3) && (qfloat < sconf->large_prime_max3))
		{

#ifdef OUTPUT_TLP_ATTEMPT_DETAILS
			FILE *fid;
			char fname[20];
			sprintf(fname, "tlp_attempts.dat");
#endif

#ifdef GATHER_RESIDUE_STATS

            // as part of analyzing 3lp parameterizations, save off this residue.  
            // These will be fully factored later and sorted
            // so that a synthetic data set resembling this real one can
            // be constructed for deeper analysis of cycle generation for
            // this parameterization.
            gmp_fprintf(sconf->residue_files[dconf->tid],
                "%Zd\n", dconf->Qvals[report_num]);

#endif

            if (dconf->do_batch)
            {
                int32_t soffset = offset;
                
                if (parity)
                {
                    soffset *= -1;
                }

                if (dconf->num_alp == 1)
                {
                    //mpz_set_ui(dconf->gmptmp1, dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                    if (dconf->curr_poly->qlisort[dconf->curr_poly->s - 1] < sconf->pmax)
                    {
                        printf("last qli entry too small: %u\n", dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                    }
                    //mpz_mul_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                    //    dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                    mpz_set_ui(dconf->gmptmp1, dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                }
                else
                {
                    mpz_set_ui(dconf->gmptmp1, 0);
                }

                // use this field to record how many we've batched.
                dconf->attempted_cosiqs++;

                relation_batch_add(dconf->curr_poly->index, poly_id, soffset, fb_offsets, smooth_num + 1,
                    dconf->Qvals[report_num], NULL, 0, dconf->gmptmp1, NULL, &dconf->rb);

                // the relation batch persists across polynomials (we save enough
                // info to know which polynomial the relation belongs to).  When we've 
                // reached our target watermark the following batch processing will
                // save valid relations off to a buffer.  That buffer in turn
                // will get merged into the top-level relation structure and
                // eventual TLP filtering.  Reset the batch structure when finished.

                // rough in a flag that we can eventually toggle from calling code,
                // forcing the batch to be processed (for example on abort or after
                // determination that we have enough relations for final filtering).
                // right now this flag isn't used and any relations in the batch
                // buffer are lost on abort or when post-processing starts.

                // one other thing we could do is switch to non-batch processing
                // when post-processing is getting close... or incrementally 
                // lower the target_relations watermark.

                // start processing the relations in the indicated buffer if running
                // single threaded.  otherwise, let the threading dispatcher assign 
                // the processing to a free thread.
                if ((dconf->batch_run_override > 0) ||
                    ((THREADS == 1) && (sconf->rb[0].num_relations > sconf->rb[0].target_relations)))
                {
                    struct timeval start;
                    struct timeval stop;
                    double ttime;
                    relation_batch_t *rb;
                    int i;

                    // if single threaded and we've reached the number of target relations.
                    // or if multi-threaded and the dispatcher has released a buffer
                    // for processing by the next available thread (us).
                    if (THREADS == 1)
                        rb = &sconf->rb[0];
                    else
                        rb = &sconf->rb[dconf->batch_run_override - 1];

                    if (VFLAG > 0)
                    {
                        printf("\nnow processing %u relations in batch %d in thread %d\n",
                            rb->num_relations, dconf->batch_run_override, dconf->tid);

                        if ((dconf->batch_run_override - 1) > 0)
                        {
                            rb->num_uecm[0] = rb->num_uecm[1] = rb->num_uecm[2] = rb->num_uecm[3] = 0;
                            rb->num_tecm = rb->num_tecm2 = 0;
                            rb->num_qs = 0;
                            for (i = 0; i < 8; i++)
                            {
                                rb->num_abort[i] = 0;
                            }
                        }
                    }

                    gettimeofday(&start, NULL);
                    relation_batch_run(rb, &dconf->lcg_state);
                    gettimeofday(&stop, NULL);

                    ttime = ytools_difftime(&start, &stop);
                    if (VFLAG > 0)
                    {
                        printf("\nrelation_batch_run took %1.4f sec producing %u tlp's\n", 
                            ttime, rb->num_success);

                        sconf->rb[0].num_attempt += rb->num_relations;
                        if ((dconf->batch_run_override - 1) > 0)
                        {
                            sconf->rb[0].num_success += rb->num_success;
                            sconf->rb[0].num_uecm[0] += rb->num_uecm[0];
                            sconf->rb[0].num_uecm[1] += rb->num_uecm[1];
                            sconf->rb[0].num_uecm[2] += rb->num_uecm[2];
                            sconf->rb[0].num_uecm[3] += rb->num_uecm[3];
                            sconf->rb[0].num_tecm += rb->num_tecm;
                            sconf->rb[0].num_tecm2 += rb->num_tecm2;
                            sconf->rb[0].num_qs += rb->num_qs;
                            for (i = 0; i < 8; i++)
                            {
                                sconf->rb[0].num_abort[i] += rb->num_abort[i];
                            }
                        }
                        printf("\tattempt: %u\n", sconf->rb[0].num_attempt);
                        printf("\tuecm1/uecm2/uecm3/uecm4/tecm1/tecm2/qs: %u,%u,%u,%u,%u,%u,%u\n", 
                            sconf->rb[0].num_uecm[0], sconf->rb[0].num_uecm[1], sconf->rb[0].num_uecm[2],
                            sconf->rb[0].num_uecm[3], sconf->rb[0].num_tecm, sconf->rb[0].num_tecm2, 
                            sconf->rb[0].num_qs);
                        printf("\taborts: %u,%u,%u,%u,%u,%u,%u,%u\n",
                            sconf->rb[0].num_abort[0], sconf->rb[0].num_abort[1], sconf->rb[0].num_abort[2],
                            sconf->rb[0].num_abort[3], sconf->rb[0].num_abort[4], sconf->rb[0].num_abort[5],
                            sconf->rb[0].num_abort[6], sconf->rb[0].num_abort[7]);
                    }

                    rb->conversion_ratio =
                          (double)rb->num_success / (double)rb->num_relations;

                    // take our new tlp relations and buffer them to be
                    // saved out to the data file.
                    for (i = 0; i < rb->num_relations; i++)
                    {
                        cofactor_t *c = rb->relations + i;
                        uint32_t *f = rb->factors + c->factor_list_word;

                        if (c->success)
                        {
                            uint8_t parity = c->signed_offset < 0 ? 1 : 0;

                            c->extra_f = c->a >> 32;
                            c->a &= 0xffffffffull;

                            if (c->extra_f > 0)
                            {
                                //printf("buffering TLP + extra factor %u\n", c->extra_f);
                                c->lp_r[c->success++] = c->extra_f;
                            }

                            if (c->success == 4)
                            {
                                dconf->qlp_useful++;
                            }
                            else if (c->success == 3)
                                dconf->tlp_useful++;
                            else if (c->success == 2)
                                dconf->dlp_useful++;
                            else if (c->success == 1)
                                dconf->num_slp++;

                            buffer_relation(abs(c->signed_offset), c->lp_r, c->num_factors_r,
                                f, c->a, c->b, parity, dconf, NULL, 0, 1);
                        }
                    }

                    if (VFLAG > 1)
                    {
                        printf("done processing batch %d in thread %d, found %u relations\n",
                            dconf->batch_run_override, dconf->tid, rb->num_success);  
                    }

                    // clear the relation batch now that we are done processing it.
                    rb->num_relations = 0;
                    rb->num_success = 0;
                    rb->num_factors = 0;

                    // signal we are done processing the batch.
                    if (THREADS > 1)
                        dconf->batch_run_override = -1;
                }
                
                // if batch factoring, we're done now.
                return;
            }
		}
		else
		{
			dconf->tlp_outside_range++;
		}
	}
	else
	{
		dconf->tlp_outside_range++;
	}

    // if we are not considering QLP, then we are done.
    if (sconf->use_dlp < 3)
        return;
	
    // quick check if Q is obviously too big.
    if ((mpz_sizeinbase(dconf->Qvals[report_num], 2) < 128) && (dconf->num_alp == 0))
    {
        double qfloat = mpz_get_d(dconf->Qvals[report_num]);

        //#define OUTPUT_TLP_ATTEMPT_DETAILS

        // not obviously too big, see if it is actually within
        // the defined qlp bounds.
        if ((qfloat > sconf->max_fb4) && (qfloat < sconf->large_prime_max4))
        {

#ifdef OUTPUT_TLP_ATTEMPT_DETAILS
            FILE* fid;
            char fname[20];
            sprintf(fname, "tlp_attempts.dat");
#endif

#ifdef GATHER_RESIDUE_STATS

            // as part of analyzing 3lp parameterizations, save off this residue.  
            // These will be fully factored later and sorted
            // so that a synthetic data set resembling this real one can
            // be constructed for deeper analysis of cycle generation for
            // this parameterization.
            gmp_fprintf(sconf->residue_files[dconf->tid],
                "%Zd\n", dconf->Qvals[report_num]);

#endif

            if (dconf->do_batch)
            {
                int32_t soffset = offset;

                if (parity)
                {
                    soffset *= -1;
                }

                if (dconf->num_alp == 1)
                {
                    printf("****** warning 4+1 LPs not supported yet\n");

                    //mpz_set_ui(dconf->gmptmp1, dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                    if (dconf->curr_poly->qlisort[dconf->curr_poly->s - 1] < sconf->pmax)
                    {
                        printf("last qli entry too small: %u\n", dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                    }
                    //mpz_mul_ui(dconf->Qvals[report_num], dconf->Qvals[report_num],
                    //    dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                    mpz_set_ui(dconf->gmptmp1, dconf->curr_poly->qlisort[dconf->curr_poly->s - 1]);
                }
                else
                {
                    mpz_set_ui(dconf->gmptmp1, 0);
                }

                // use this field to record how many we've batched.
                dconf->attempted_cosiqs++;
                relation_batch_add(dconf->curr_poly->index, poly_id, soffset, fb_offsets, smooth_num + 1,
                    dconf->Qvals[report_num], NULL, 0, dconf->gmptmp1, NULL, &dconf->rb);

                // the relation batch persists across polynomials (we save enough
                // info to know which polynomial the relation belongs to).  When we've 
                // reached our target watermark the following batch processing will
                // save valid relations off to a buffer.  That buffer in turn
                // will get merged into the top-level relation structure and
                // eventual TLP filtering.  Reset the batch structure when finished.

                // rough in a flag that we can eventually toggle from calling code,
                // forcing the batch to be processed (for example on abort or after
                // determination that we have enough relations for final filtering).
                // right now this flag isn't used and any relations in the batch
                // buffer are lost on abort or when post-processing starts.

                // one other thing we could do is switch to non-batch processing
                // when post-processing is getting close... or incrementally 
                // lower the target_relations watermark.

                // start processing the relations in the indicated buffer if running
                // single threaded.  otherwise, let the threading dispatcher assign 
                // the processing to a free thread.
                if ((dconf->batch_run_override > 0) ||
                    ((THREADS == 1) && (sconf->rb[0].num_relations > sconf->rb[0].target_relations)))
                {
                    struct timeval start;
                    struct timeval stop;
                    double ttime;
                    relation_batch_t* rb;
                    int i;

                    // if single threaded and we've reached the number of target relations.
                    // or if multi-threaded and the dispatcher has released a buffer
                    // for processing by the next available thread (us).
                    if (THREADS == 1)
                        rb = &sconf->rb[0];
                    else
                        rb = &sconf->rb[dconf->batch_run_override - 1];

                    if (VFLAG > 0)
                    {
                        printf("\nnow processing %u relations in batch %d in thread %d\n",
                            rb->num_relations, dconf->batch_run_override, dconf->tid);

                        if ((dconf->batch_run_override - 1) > 0)
                        {
                            rb->num_uecm[0] = rb->num_uecm[1] = rb->num_uecm[2] = rb->num_uecm[3] = 0;
                            rb->num_tecm = rb->num_tecm2 = 0;
                            rb->num_qs = 0;
                            for (i = 0; i < 8; i++)
                            {
                                rb->num_abort[i] = 0;
                            }
                        }
                    }

                    gettimeofday(&start, NULL);
                    relation_batch_run(rb, &dconf->lcg_state);
                    gettimeofday(&stop, NULL);

                    ttime = ytools_difftime(&start, &stop);
                    if (VFLAG > 0)
                    {
                        printf("\nrelation_batch_run took %1.4f sec producing %u tlp's\n",
                            ttime, rb->num_success);

                        sconf->rb[0].num_attempt += rb->num_relations;
                        if ((dconf->batch_run_override - 1) > 0)
                        {
                            sconf->rb[0].num_success += rb->num_success;
                            sconf->rb[0].num_uecm[0] += rb->num_uecm[0];
                            sconf->rb[0].num_uecm[1] += rb->num_uecm[1];
                            sconf->rb[0].num_uecm[2] += rb->num_uecm[2];
                            sconf->rb[0].num_uecm[3] += rb->num_uecm[3];
                            sconf->rb[0].num_tecm += rb->num_tecm;
                            sconf->rb[0].num_tecm2 += rb->num_tecm2;
                            sconf->rb[0].num_qs += rb->num_qs;
                            for (i = 0; i < 8; i++)
                            {
                                sconf->rb[0].num_abort[i] += rb->num_abort[i];
                            }
                        }
                        printf("\tattempt: %u\n", sconf->rb[0].num_attempt);
                        printf("\tuecm1/uecm2/uecm3/uecm4/tecm1/tecm2/qs: %u,%u,%u,%u,%u,%u,%u\n",
                            sconf->rb[0].num_uecm[0], sconf->rb[0].num_uecm[1], sconf->rb[0].num_uecm[2],
                            sconf->rb[0].num_uecm[3], sconf->rb[0].num_tecm, sconf->rb[0].num_tecm2,
                            sconf->rb[0].num_qs);
                        printf("\taborts: %u,%u,%u,%u,%u,%u,%u,%u\n",
                            sconf->rb[0].num_abort[0], sconf->rb[0].num_abort[1], sconf->rb[0].num_abort[2],
                            sconf->rb[0].num_abort[3], sconf->rb[0].num_abort[4], sconf->rb[0].num_abort[5],
                            sconf->rb[0].num_abort[6], sconf->rb[0].num_abort[7]);
                    }

                    rb->conversion_ratio =
                        (double)rb->num_success / (double)rb->num_relations;

                    // take our new tlp relations and buffer them to be
                    // saved out to the data file.
                    for (i = 0; i < rb->num_relations; i++)
                    {
                        cofactor_t* c = rb->relations + i;
                        uint32_t* f = rb->factors + c->factor_list_word;

                        if (c->success)
                        {
                            uint8_t parity = c->signed_offset < 0 ? 1 : 0;

                            c->extra_f = c->a >> 32;
                            c->a &= 0xffffffffull;

                            if (c->extra_f > 0)
                            {
                                //printf("buffering TLP + extra factor %u\n", c->extra_f);
                                c->lp_r[c->success++] = c->extra_f;
                            }

                            if (c->success == 4)
                            {
                                dconf->qlp_useful++;
                            }
                            else if (c->success == 3)
                                dconf->tlp_useful++;
                            else if (c->success == 2)
                                dconf->dlp_useful++;
                            else if (c->success == 1)
                                dconf->num_slp++;

                            buffer_relation(abs(c->signed_offset), c->lp_r, c->num_factors_r,
                                f, c->a, c->b, parity, dconf, NULL, 0, 1);
                        }
                    }

                    if (VFLAG > 1)
                    {
                        printf("done processing batch %d in thread %d, found %u relations\n",
                            dconf->batch_run_override, dconf->tid, rb->num_success);
                    }

                    // clear the relation batch now that we are done processing it.
                    rb->num_relations = 0;
                    rb->num_success = 0;
                    rb->num_factors = 0;

                    // signal we are done processing the batch.
                    if (THREADS > 1)
                        dconf->batch_run_override = -1;
                }

                // if batch factoring, we're done now.
                return;
            }

        }
        else
        {
            dconf->qlp_outside_range++;
        }
    }
    else
    {
        dconf->qlp_outside_range++;
    }

	return;
}

void buffer_relation(uint32_t offset, uint32_t *large_prime, uint32_t num_factors,
    uint32_t *fb_offsets, uint32_t apoly_id, uint32_t bpoly_id, uint32_t parity,
    dynamic_conf_t *conf, uint32_t *polya_factors,
    uint32_t num_polya_factors, uint64_t unfactored_residue)
{
    // put this relations's info into a temporary buffer which
    // will get merged with other such buffers (if multi-threaded) and
    // dumped to file once the threads are joined.
    siqs_r *rel;
    uint32_t i, j, k;

    // first check that this relation won't overflow the buffer
    if (conf->buffered_rels >= conf->buffered_rel_alloc)
    {
        //printf("reallocating relation buffer\n");
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
    rel->apoly_idx = apoly_id;
    rel->poly_idx = bpoly_id;

    if ((num_polya_factors + num_factors) > MAX_SMOOTH_PRIMES)
    {
        printf("error: too many smooth primes!\n");
        exit(234);
    }

#ifdef SPARSE_STORE
    // extra factors of polya are not added to the list of factors. they will be added
    // during filtering on only the relations that survive singleton removal.
    memcpy(rel->fb_offsets, fb_offsets, sizeof(uint32_t) * num_factors);
    rel->num_factors = num_factors;

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

    rel->num_factors = num_factors + num_polya_factors;
#endif

    rel->large_prime[0] = large_prime[0];
	rel->large_prime[1] = large_prime[1];
	rel->large_prime[2] = large_prime[2];
    rel->large_prime[3] = large_prime[3];

	conf->buffered_rels++;
	return;
}

void save_relation_siqs(uint32_t offset, uint32_t *large_prime, uint32_t num_factors, 
    uint32_t *fb_offsets, uint32_t poly_id, uint32_t apoly_id, uint32_t parity,
	static_conf_t *conf)
{
	char buf[1024];
	fact_obj_t *obj = conf->obj;
	uint32_t i, k, buf_offset;
	uint32_t lp[MAXLP];
    uint32_t numlp = conf->num_lp + NUM_ALP;

	if (conf->in_mem)
	{
		// copy info to sconf relation structure
		siqs_r *r;

		//first check that this relation won't overflow the buffer
		if (conf->buffered_rels >= conf->buffered_rel_alloc)
		{
			//printf("reallocating in-mem relation storage...\n");
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

        // insertion sort 3lp or more
        if (conf->num_lp > 2)
        {
            i = 1;
            while (i < conf->num_lp)
            {
                uint32_t j = i;
                while ((j > 0) && (large_prime[j - 1] > large_prime[j]))
                {
                    uint32_t tmp = large_prime[j - 1];
                    large_prime[j - 1] = large_prime[j];
                    large_prime[j] = tmp;
                    j = j - 1;
                }
                i = i + 1;
            }

            for (i = 0; i < conf->num_lp; i++)
            {
                r->large_prime[i] = large_prime[i];
            }

        }
        else if (conf->num_lp == 2)
        {
            for (i = 0; i < 4; i++)
            {
                r->large_prime[i] = 1;
            }

            if (large_prime[0] < large_prime[1])
            {
                r->large_prime[0] = large_prime[0];
                r->large_prime[1] = large_prime[1];
            }
            else
            {
                r->large_prime[1] = large_prime[0];
                r->large_prime[0] = large_prime[1];
            }
        }
        else
        {
            for (i = 1; i < 4; i++)
            {
                r->large_prime[i] = 1;
            }

            r->large_prime[0] = large_prime[0];
        }

		r->num_factors = num_factors;
        r->apoly_idx = apoly_id;
		r->poly_idx = poly_id;
		r->parity = parity;
		r->sieve_offset = offset;
        if (num_factors > MAX_SMOOTH_PRIMES)
        {
            printf("error: num_factors (%u) exceeds maximum (%u)!\n",
                num_factors, MAX_SMOOTH_PRIMES);
        }

		for (i=0; i < MIN(num_factors, MAX_SMOOTH_PRIMES); i++)
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

        buf_offset = i;

		// insertion sort 3lp or more
        if (numlp > 2)
        {
            i = 1;
            while (i < numlp)
            {
                uint32_t j = i;
                while ((j > 0) && (large_prime[j - 1] > large_prime[j]))
                {
                    uint32_t tmp = large_prime[j - 1];
                    large_prime[j - 1] = large_prime[j];
                    large_prime[j] = tmp;
                    j = j - 1;
                }
                i = i + 1;
            }

            if (numlp == 3)
            {
                sprintf(buf + buf_offset, "L %x %x %x\n", large_prime[0], large_prime[1], large_prime[2]);
            }
            else if (numlp == 4)
            {
                sprintf(buf + buf_offset, "L %x %x %x %x\n", large_prime[0], large_prime[1],
                    large_prime[2], large_prime[3]);
            }
        }
        else
        {
            if (large_prime[0] < large_prime[1])
                sprintf(buf + buf_offset, "L %x %x\n", large_prime[0], large_prime[1]);
            else
                sprintf(buf + buf_offset, "L %x %x\n", large_prime[1], large_prime[0]);
        }

		savefile_write_line(&obj->qs_obj.savefile, buf);		
	}

	/* for partial relations, also update the bookeeping for
		   tracking the number of fundamental cycles */

	if (conf->use_dlp >= 2)
	{
        // when using TLP, cycle counting is a lot more complicated.
        // hand it off here.
		yafu_add_to_cyclesN(conf, 0, lp);
	}
	else
	{
		if (large_prime[0] != large_prime[1]) {
			yafu_add_to_cycles(conf, obj->flags, large_prime[0], large_prime[1]);
			conf->num_cycles++;
		}
		else {
			conf->num_relations++;
		}
	}

	return;
}

