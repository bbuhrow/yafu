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
#include "msieve_common.h"

/*
this file contains code to implement filtering of relation sets prior
to matrix construction, and restarting of siqs from a data file
*/

uint32_t process_poly_a(static_conf_t *sconf)
{
	//given a poly a value, and some aux info about the factorization
	//generate all poly b values associated with that a and
	//store both a and all the b's in conf
	//also check that the generated polys are valid.

	//we will be reusing some routines here that are normally used
	//during sieving, and expect a dynamic_conf structure as input which
	//we don't have here.  So we need to create one for use in this 
	//routine only.
	dynamic_conf_t *dconf;
	int maxB,j,i;

	dconf = (dynamic_conf_t *)malloc(sizeof(dynamic_conf_t));

	//just initialize the stuff needed for processing poly's, namely, 
	//the curr_poly structure, the Bl array, and a few temp bigints.
	mpz_init(dconf->gmptmp1);
	mpz_init(dconf->gmptmp2);

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

	//allocate the Bl array, space for MAX_Bl bigint numbers
	dconf->Bl = (mpz_t *)malloc(MAX_A_FACTORS * sizeof(mpz_t));
	for (i=0;i<MAX_A_FACTORS;i++)
		mpz_init(dconf->Bl[i]);

	mpz_set(dconf->curr_poly->mpz_poly_a, sconf->curr_a);
	
	//then compute all the 'b' poly's for this 'a'
	//and add them to the b-list

	//first, need to factorize this 'a'
	j = get_a_offsets(sconf->factor_base, dconf->curr_poly, dconf->gmptmp1);
	if (j)
	{
		//then this poly a might be corrupted - we couldn't factor it over the factor base
		free(dconf->curr_poly->gray);
		free(dconf->curr_poly->nu);
		free(dconf->curr_poly->qlisort);
		mpz_clear(dconf->curr_poly->mpz_poly_a);
		mpz_clear(dconf->curr_poly->mpz_poly_b);
		mpz_clear(dconf->curr_poly->mpz_poly_c);
		free(dconf->curr_poly);

		for (i=0;i<MAX_A_FACTORS;i++)
			mpz_clear(dconf->Bl[i]);
		free(dconf->Bl);

		//workspace bigints
		mpz_clear(dconf->gmptmp1);
		mpz_clear(dconf->gmptmp2);
		free(dconf);
		return 0;
	}

	//then initialize the gray code
	get_gray_code(dconf->curr_poly);

	//then compute all the Bl's
	//the first 'b' poly comes with computeBl
	computeBl(sconf, dconf);

	//compute how many 'b' values we can get from this 'a'
	maxB = 1 << (dconf->curr_poly->s - 1);

	//now we copy all the b coefficients over to sconf where they are
	//needed by the rest of the filtering routine.
	//make sure there is enough room for them.  curr_b could
	//currently be allocated for a number of poly smaller or bigger than
	//maxB

	//free any we won't be needing
	for (j = 0; (uint32_t)j < sconf->bpoly_alloc; j++)
		mpz_clear(sconf->curr_b[j]);

	//reallocate the size of the array
	sconf->curr_b = (mpz_t *)realloc(sconf->curr_b, maxB * sizeof(mpz_t));

	//allocate any additional we need
	for (j = 0; j < maxB; j++)
		mpz_init(sconf->curr_b[j]);
	sconf->bpoly_alloc = maxB;

	//generate all the b polys
	generate_bpolys(sconf, dconf, maxB);

	//we'll need to remember some things about the current poly,
	//so copy those over to sconf first...
	for (j=0; j<dconf->curr_poly->s; j++)
		sconf->curr_poly->qlisort[j] = dconf->curr_poly->qlisort[j];
	sconf->curr_poly->s = dconf->curr_poly->s;

	//then free the temp dynamic struct
	free(dconf->curr_poly->gray);
	free(dconf->curr_poly->nu);
	free(dconf->curr_poly->qlisort);
	mpz_clear(dconf->curr_poly->mpz_poly_a);
	mpz_clear(dconf->curr_poly->mpz_poly_b);
	mpz_clear(dconf->curr_poly->mpz_poly_c);
	free(dconf->curr_poly);

	for (i=0;i<MAX_A_FACTORS;i++)
		mpz_clear(dconf->Bl[i]);
	free(dconf->Bl);

	//workspace bigints
	mpz_clear(dconf->gmptmp1);
	mpz_clear(dconf->gmptmp2);
	free(dconf);

	return maxB;
}

int get_a_offsets(fb_list *fb, siqs_poly *poly, mpz_t tmp)
{
	int j,k;
	uint64_t r;

	mpz_set(tmp, poly->mpz_poly_a);
	k=0;
	poly->s = 0;
	while ((mpz_cmp_ui(tmp, 1) != 0) && (siqs_primes[k] < 65536))
	{
		r = mpz_tdiv_ui(tmp,(uint64_t)siqs_primes[k]);
		
		if (r != 0)
			k++;
		else
		{
			for (j=1;(uint32_t)j < fb->B;j++)
				if (fb->list->prime[j] == siqs_primes[k])
					break;
			if ((uint32_t)j >= fb->B)
			{
				//then we didn't find this factor in the fb, thus the 
				//read in A is probably bad.
				printf("bad poly A encountered in savefile\n");
				return 1;
			}
			poly->qlisort[poly->s] = j;
			poly->s++;
			mpz_tdiv_q_ui(tmp, tmp, (uint64_t)siqs_primes[k]);
		}
	}

	return 0;
}

void generate_bpolys(static_conf_t *sconf, dynamic_conf_t *dconf, int maxB)
{
	//given poly, which contains the first value of poly_b, and
	//info needed to generate the rest of the poly b's (namely,
	//the Bl vector), generate the rest of the Bl's and put
	//them in an output vector, which has been initialized 
	//by the callee.

	int numB=1;

	//iterate through all b's
	for ( ; numB < maxB; numB++)
	{
		//zCopy(&dconf->curr_poly->poly_b,&sconf->curr_b[numB - 1]);
		mpz_set(sconf->curr_b[numB - 1], dconf->curr_poly->mpz_poly_b);
		dconf->numB = numB;
		nextB(dconf, sconf);
	}

	return;
}

void td_and_merge_relation(fb_list *fb, mpz_t n,
    static_conf_t *sconf, fact_obj_t *obj, 
    siqs_r *r_out, siqs_r *rel)
{
    int i, j, k, err_code = 0;
    mpz_t Q;
    mpz_init(Q);

    // trial divide to find divisible primes less than med_B.
    // this will include small primes and primes dividing poly_a.

    //Q(x)/a = (ax + b)^2 - N, where x is the sieve index
    mpz_mul_ui(Q, sconf->curr_a, rel->sieve_offset);
    if (rel->parity)
        mpz_sub(Q, Q, sconf->curr_b[rel->poly_idx]);
    else
        mpz_add(Q, Q, sconf->curr_b[rel->poly_idx]);
    mpz_mul(Q, Q, Q);
    mpz_sub(Q, Q, n);

    if (mpz_sgn(Q) < 0)
        mpz_neg(Q, Q);

    for (i = 1, k = 0; i < sconf->sieve_small_fb_start; i++)
    {
        uint32_t prime = sconf->factor_base->list->prime[i];
        while (mpz_tdiv_ui(Q, prime) == 0)
        {
            r_out->fb_offsets[k++] = i;
            mpz_tdiv_q_ui(Q, Q, prime);
        }
    }

    mpz_tdiv_q(Q, Q, sconf->curr_a);

    // merge in the list of polya coefficients to the list of
    // relation factors that we read in.  both lists are sorted,
    // so this mergesort works quickly.  We didn't add extra factors
    // of polya during sieving, so search for and add them here as
    // they are merged into the complete list of factors for this relation.
    i = j = 0;
    while ((i < (int)rel->num_factors) && (j < sconf->curr_poly->s)) {
        uint32_t prime;

        if (rel->fb_offsets[i] < sconf->curr_poly->qlisort[j]) {
            r_out->fb_offsets[k++] = rel->fb_offsets[i++];
        }
        else if (rel->fb_offsets[i] > sconf->curr_poly->qlisort[j]) {
            // add in the one factor that will always be there because a | Q
            r_out->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
            prime = sconf->factor_base->list->prime[sconf->curr_poly->qlisort[j]];

            // then test and add more if we can
            while (mpz_tdiv_ui(Q, prime) == 0)
            {
                r_out->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
                mpz_tdiv_q_ui(Q, Q, prime);
            }
            j++;
        }
        else {
            r_out->fb_offsets[k] = rel->fb_offsets[i++];
            r_out->fb_offsets[k + 1] = sconf->curr_poly->qlisort[j];
            prime = sconf->factor_base->list->prime[sconf->curr_poly->qlisort[j]];
            k += 2;
            while (mpz_tdiv_ui(Q, prime) == 0)
            {
                r_out->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
                mpz_tdiv_q_ui(Q, Q, prime);
            }
            j++;
        }
    }

    while (i < (int)rel->num_factors)
    {
        r_out->fb_offsets[k++] = rel->fb_offsets[i++];
    }

    while (j < sconf->curr_poly->s)
    {
        r_out->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
        uint32_t prime = sconf->factor_base->list->prime[sconf->curr_poly->qlisort[j]];

        // then test and add more if we can
        while (mpz_tdiv_ui(Q, prime) == 0)
        {
            r_out->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
            mpz_tdiv_q_ui(Q, Q, prime);
        }
        j++;
    }

    r_out->sieve_offset = rel->sieve_offset;
    r_out->large_prime[0] = rel->large_prime[0];
    r_out->large_prime[1] = rel->large_prime[1];
    r_out->large_prime[2] = rel->large_prime[2];
    r_out->parity = rel->parity;
    r_out->num_factors = k;
    r_out->poly_idx = rel->poly_idx;
    mpz_clear(Q);

    if (sconf->use_dlp == 2)
    {
        if (!check_relation(sconf->curr_a,
            sconf->curr_b[r_out->poly_idx], r_out, fb, n, obj->VFLAG))
        {
            yafu_add_to_cycles3(sconf, 0, r_out->large_prime);
        }
        else
        {
            //printf("relation string: %s\n",instr);
            err_code = 1;
        }
    }
    else
    {
        if (!check_relation(sconf->curr_a,
            sconf->curr_b[rel->poly_idx], r_out, fb, n, obj->VFLAG))
        {
            if (rel->large_prime[0] != r_out->large_prime[1]) {
                yafu_add_to_cycles(sconf, obj->flags, 
                    r_out->large_prime[0], r_out->large_prime[1]);
                sconf->num_cycles++;
            }
            else {
                sconf->num_relations++;
            }
        }
        else
        {
            //printf("relation string: %s\n",instr);
            err_code = 1;
        }
    }
    return;
}

int process_rel(char *substr, fb_list *fb, mpz_t n,
				 static_conf_t *sconf, fact_obj_t*obj, siqs_r *rel)
{
	char *nextstr;
	uint32_t lp[3];
	uint32_t this_offset, this_id, this_num_factors, this_parity, this_val;
	uint32_t fb_offsets[MAX_SMOOTH_PRIMES];
	int i,j,k, err_code = 0;

    // assumes SPARSE_STORE
    mpz_t Q;
    mpz_init(Q);

    //read the offset
    this_offset = strtoul(substr, &nextstr, 16);	//convert
    substr = nextstr;

    this_id = strtoul(substr, &nextstr, 16);
    substr = nextstr;

    if (this_offset & 0x80000000)
    {
        this_parity = 1;
        this_offset = (~this_offset) + 1;
    }
    else
    {
        this_parity = 0;
    }

    j = 0;
    do
    {
        this_val = strtoul(substr, &nextstr, 16);
        if (this_val == (uint32_t)(-1))
        {
            printf("error parsing relation: strtoul returned error code\n");
            continue;
        }
        substr = nextstr;
        fb_offsets[j] = this_val;
        j++;
    } while (substr[1] != 'L');
    this_num_factors = j;

    substr += 2;
    this_val = strtoul(substr, &nextstr, 16);
    substr = nextstr;
    lp[0] = this_val;

    this_val = strtoul(substr, &nextstr, 16);
    substr = nextstr;
    lp[1] = this_val;

    if (sconf->use_dlp == 2)
    {
        this_val = strtoul(substr, &nextstr, 16);
        substr = nextstr;
        lp[2] = this_val;
    }
    else
    {
        lp[2] = 1;
    }

    // combine the factors of the sieve value with
    // the factors of the polynomial 'a' value; the 
    // linear algebra code has to know about both.
    // Because both lists are sorted, this is just
    // a merge operation

    // trial divide to find divisible primes less than med_B.
    // this will include small primes and primes dividing poly_a.

    //Q(x)/a = (ax + b)^2 - N, where x is the sieve index
    mpz_mul_ui(Q, sconf->curr_a, this_offset);
    if (this_parity)
        mpz_sub(Q, Q, sconf->curr_b[this_id]);
    else
        mpz_add(Q, Q, sconf->curr_b[this_id]);
    mpz_mul(Q, Q, Q);
    mpz_sub(Q, Q, n);

    if (mpz_sgn(Q) < 0)
        mpz_neg(Q, Q);


    for (i = 1, k = 0; i < sconf->sieve_small_fb_start; i++)
    {
        uint32_t prime = sconf->factor_base->list->prime[i];
        while (mpz_tdiv_ui(Q, prime) == 0)
        {
            rel->fb_offsets[k++] = i;
            mpz_tdiv_q_ui(Q, Q, prime);
        }
    }

    mpz_tdiv_q(Q, Q, sconf->curr_a);

    // merge in the list of polya coefficients to the list of
    // relation factors that we read in.  both lists are sorted,
    // so this mergesort works quickly.  We didn't add extra factors
    // of polya during sieving, so search for and add them here as
    // they are merged into the complete list of factors for this relation.
    i = j = 0;
    while ((i < (int)this_num_factors) && (j < sconf->curr_poly->s)) {
        uint32_t prime;

        if (fb_offsets[i] < sconf->curr_poly->qlisort[j]) {
            rel->fb_offsets[k++] = fb_offsets[i++];
        }
        else if (fb_offsets[i] > sconf->curr_poly->qlisort[j]) {
            // add in the one factor that will always be there because a | Q
            rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
            prime = sconf->factor_base->list->prime[sconf->curr_poly->qlisort[j]];

            // then test and add more if we can
            while (mpz_tdiv_ui(Q, prime) == 0)
            {
                rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
                mpz_tdiv_q_ui(Q, Q, prime);
            }
            j++;
        }
        else {
            rel->fb_offsets[k] = fb_offsets[i++];
            rel->fb_offsets[k + 1] = sconf->curr_poly->qlisort[j];
            prime = sconf->factor_base->list->prime[sconf->curr_poly->qlisort[j]];
            k += 2;
            while (mpz_tdiv_ui(Q, prime) == 0)
            {
                rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
                mpz_tdiv_q_ui(Q, Q, prime);
            }
            j++;
        }
    }

    while (i < (int)this_num_factors)
        rel->fb_offsets[k++] = fb_offsets[i++];

    while (j < sconf->curr_poly->s)
    {
        rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
        uint32_t prime = sconf->factor_base->list->prime[sconf->curr_poly->qlisort[j]];

        // then test and add more if we can
        while (mpz_tdiv_ui(Q, prime) == 0)
        {
            rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j];
            mpz_tdiv_q_ui(Q, Q, prime);
        }
        j++;
    }

    this_num_factors = k;
    mpz_clear(Q);

	rel->sieve_offset = this_offset;
	rel->large_prime[0] = lp[0];
	rel->large_prime[1] = lp[1];
	rel->large_prime[2] = lp[2];
	rel->parity = this_parity;
    rel->num_factors = this_num_factors; 
	rel->poly_idx = this_id;

	if (sconf->use_dlp == 2)
	{
		if (!check_relation(sconf->curr_a,
			sconf->curr_b[rel->poly_idx], rel, fb, n, obj->VFLAG))
		{
			yafu_add_to_cycles3(sconf, 0, lp);
		}
		else
		{
			//printf("relation string: %s\n",instr);
			err_code = 1;
		}
	}
	else
	{
		if (!check_relation(sconf->curr_a,
			sconf->curr_b[rel->poly_idx], rel, fb, n, obj->VFLAG))
		{
			if (lp[0] != lp[1]) {
				yafu_add_to_cycles(sconf, obj->flags, lp[0], lp[1]);
				sconf->num_cycles++;
			}
			else {
				sconf->num_relations++;
			}
		}
		else
		{
			//printf("relation string: %s\n",instr);
			err_code = 1;
		}
	}


	return err_code;	//error code, if there is one.
}

int restart_siqs(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	int i,j;
	char *str, *substr;
	FILE *data;
	uint32_t lp[2],pmax = sconf->large_prime_max / sconf->large_mult;
    fact_obj_t*fobj = sconf->obj;

	str = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	data = fopen(sconf->obj->qs_obj.siqs_savefile,"r");
	i=0;
	j=0;
	
	if (data != NULL)
	{	
		fgets(str,1024,data);
		substr = str + 2;
		mpz_set_str(dconf->gmptmp1, substr, 0); //str2hexz(substr,&dconf->qstmp1);
		// check against the input to SIQS, i.e., does not have a 
		// multiplier applied (the file saved N does not include the multiplier).
		if (mpz_cmp(dconf->gmptmp1, sconf->obj->qs_obj.gmp_n) == 0)
		{
			if (fobj->VFLAG > 1)
				printf("restarting siqs from saved data set\n");

			fflush(stdout);
			fflush(stderr);

			if (sconf->use_dlp == 2)
			{
				siqs_r *relation_list;
				qs_la_col_t *cycle_list;
				uint32_t num_cycles;
				uint32_t *hashtable = sconf->cycle_hashtable;
				qs_cycle_t *table = sconf->cycle_table;
				uint32_t num_relations;
				uint32_t i = 0, passes;
				uint32_t curr_a_idx, curr_poly_idx, curr_rel;
				uint32_t curr_expected, curr_saved, curr_cycle;
				uint32_t all_relations;
				uint32_t total_poly_a = 0;
				uint32_t poly_saved;
				uint32_t *plist0;
				uint32_t *plist1;
				uint32_t *plist2;
				int j;

				printf("reading relations\n");

				// we don't know beforehand how many rels to expect, so start
				// with some amount and allow it to increase as we read them
				relation_list = (siqs_r *)xmalloc(10000 * sizeof(siqs_r));
				plist0 = (uint32_t *)xmalloc(10000 * sizeof(uint32_t));
				plist1 = (uint32_t *)xmalloc(10000 * sizeof(uint32_t));
				plist2 = (uint32_t *)xmalloc(10000 * sizeof(uint32_t));
				curr_rel = 10000;
				while (!feof(data)) {
					char *start;

					fgets(str, 1024, data);
					substr = str + 2;

					switch (str[0]) {
					case 'A':
						total_poly_a++;
						break;

					case 'R':
						start = strchr(str, 'L');
						if (start != NULL) {
							uint32_t primes[3];
							yafu_read_tlp(start, primes);
							if (i == curr_rel) {
								curr_rel = 3 * curr_rel / 2;
								relation_list = (siqs_r *)xrealloc(
									relation_list,
									curr_rel *
									sizeof(siqs_r));

								plist0 = (uint32_t *)xrealloc(plist0, curr_rel * sizeof(uint32_t));
								plist1 = (uint32_t *)xrealloc(plist1, curr_rel * sizeof(uint32_t));
								plist2 = (uint32_t *)xrealloc(plist2, curr_rel * sizeof(uint32_t));
							}

							//printf("found primes %u,%u,%u on line %u, current allocation: %u\n",
							//	primes[0], primes[1], primes[2], i, curr_rel);

							relation_list[i].poly_idx = i;
							relation_list[i].num_factors = lcg_rand_32_range(0, 1000000000, &dconf->lcg_state);
							relation_list[i].large_prime[0] = primes[0];
							relation_list[i].large_prime[1] = primes[1];
							relation_list[i].large_prime[2] = primes[2];

							plist0[0] = primes[0];
							plist1[1] = primes[1];
							plist2[2] = primes[2];
							i++;
						}
						break;
					case 'N':
						break;
					}
				}

				all_relations = i;
				num_relations = i;

				printf("read %u relations\n", i);
				printf("building graph\n");

				rebuild_graph3(sconf, relation_list, num_relations);

				hashtable = sconf->cycle_hashtable;
				table = sconf->cycle_table;

				printf("commencing singleton removal\n");

				printf("cycle table size = %u, vertices = %u, components = %u\n",
					sconf->cycle_table_size, sconf->vertices, sconf->components);

				num_relations = qs_purge_singletons3(sconf->obj, relation_list,
					num_relations, table, hashtable);

				printf("%u relations survived singleton removal\n", num_relations);

				if (num_relations > 0)
				{
					printf("commencing duplicate removal\n");
					
					num_relations = qs_purge_duplicate_relations3(sconf->obj,
						relation_list, num_relations);
					
					printf("%u relations survived duplicate removal\n", num_relations);

					printf("building reduced graph\n");

					rebuild_graph3(sconf, relation_list, num_relations);

					printf("cycle table size = %u, edges = %u, vertices = %u, components = %u\n",
						sconf->cycle_table_size, sconf->num_cycles,
						sconf->vertices, sconf->components);

					printf("commencing cycle find\n");

					num_cycles = sconf->vertices;
					cycle_list = find_cycles3(sconf->obj, sconf, relation_list,
						num_relations, &num_cycles, &passes);

					sconf->last_numfull = sconf->num_relations;
					sconf->last_numpartial = sconf->num_r - sconf->num_relations;
					sconf->last_numcycles = num_relations;

                    double percent_complete =
                        ((double)sconf->last_numfull + (double)sconf->last_numpartial) /
                        (double)sconf->factor_base->B;

                    if (sconf->digits_n < 120)
                    {
                        if (percent_complete < 0.2)
                        {
                            sconf->check_inc = 0.05 * sconf->factor_base->B;
                        }
                        else if (percent_complete < 0.33)
                        {
                            sconf->check_inc = 0.03 * sconf->factor_base->B;
                        }
                        else if (percent_complete < 0.5)
                        {
                            sconf->check_inc = 0.02 * sconf->factor_base->B;
                        }
                        else
                        {
                            sconf->check_inc = 0.01 * sconf->factor_base->B;
                        }
                    }

                    if (sconf->digits_n > 120)
                    {
                        if (percent_complete < 0.2)
                        {
                            sconf->check_inc = 0.03 * sconf->factor_base->B;
                        }
                        else if (percent_complete < 0.33)
                        {
                            sconf->check_inc = 0.02 * sconf->factor_base->B;
                        }
                        else if (percent_complete < 0.5)
                        {
                            sconf->check_inc = 0.01 * sconf->factor_base->B;
                        }
                        else
                        {
                            sconf->check_inc = 0.005 * sconf->factor_base->B;
                        }
                    }

					for (j = 0; j < num_cycles; j++)
					{
						free(cycle_list[j].cycle.list);
					}
					free(cycle_list);

					printf("found %u cycles in %u passes\n", num_cycles, passes);
					printf("restoring full %u relation set\n", all_relations);

					for (j = 0; j < all_relations; j++)
					{
						relation_list[j].large_prime[0] = plist0[0];
						relation_list[j].large_prime[1] = plist1[1];
						relation_list[j].large_prime[2] = plist2[2];
					}
					num_relations = all_relations;

					printf("rebuilding full graph\n");

					rebuild_graph3(sconf, relation_list, num_relations);
				}
				else
				{
					printf("restoring full %u relation set\n", all_relations);

					for (j = 0; j < all_relations; j++)
					{
						relation_list[j].large_prime[0] = plist0[0];
						relation_list[j].large_prime[1] = plist1[1];
						relation_list[j].large_prime[2] = plist2[2];
					}
					num_relations = all_relations;

					printf("rebuilding full graph\n");

					rebuild_graph3(sconf, relation_list, num_relations);
					sconf->num_cycles = num_relations;
				}

				free(plist0);
				free(plist1);
				free(plist2);
				free(relation_list);
			}
			else
			{
				
				while (1)
				{
					//read a line
					if (feof(data))
						break;
					fgets(str, 1024, data);
					substr = str + 2;

					if (str[0] == 'R')
					{
						//process a relation
						//just trying to figure out how many relations we have
						//so read in the large primes and add to cycles
						substr = strchr(substr, 'L');
						yafu_read_large_primes(substr, lp, lp + 1);
						if (sconf->use_dlp)
						{
							//if ((lp[0] > 1) && (lp[0] < pmax))
							//{
							//	j++;
							//	continue;
							//}
							//if ((lp[1] > 1) && (lp[1] < pmax))
							//{
							//	j++;
							//	continue;
							//}
						}
						if (lp[0] != lp[1])
						{
							yafu_add_to_cycles(sconf, sconf->obj->flags, lp[0], lp[1]);
							sconf->num_cycles++;
						}
						else {
							sconf->num_relations++;
						}
					}
					else if (str[0] == 'A')
					{
						i++;
					}
				}

				sconf->num_r = sconf->num_relations +
					sconf->num_cycles +
					sconf->components - sconf->vertices;

				if (fobj->VFLAG > 0)
				{
					printf("%d relations found: %d full + "
						"%d from %d partial\n",
						sconf->num_r, sconf->num_relations,
						sconf->num_cycles +
						sconf->components - sconf->vertices,
						sconf->num_cycles);
					printf("threw away %d relations with large primes too small\n", j);
					fflush(stdout);
					sconf->last_numfull = sconf->num_relations;
					sconf->last_numcycles = sconf->num_cycles;
				}
			}
		}	
		fclose(data);
	}
	free(str);

	return 0;
}

#define QS_HASH_MULT ((uint32_t)(2654435761UL))
#define QS_HASH_ADD ((uint32_t)(18932479UL))
#define QS_HASH(a) (((a) * QS_HASH_MULT) >> (32 - QS_LOG2_CYCLE_HASH))
#define PBR_HASH(a) ((((a) + QS_HASH_ADD) * QS_HASH_MULT) >> (32 - QS_LOG2_CYCLE_HASH))

static qs_cycle_t *get_table_entry(qs_cycle_t *table, uint32_t *hashtable,
    uint32_t prime, uint32_t new_entry_offset) {

    /* return a pointer to a unique qs_cycle_t specific
    to 'prime'. The value of 'prime' is hashed and
    the result used to index 'table'. If prime does
    not appear in table, specify 'new_entry_offset'
    as the table entry for prime.

    Values of 'prime' which hash to the same offset
    in hashtable are connected by a linked list of
    offsets in 'table'. */

    uint32_t offset, first_offset;
    qs_cycle_t *entry = NULL;

    first_offset = QS_HASH(prime);
    offset = hashtable[first_offset];

    /* follow the list of entries corresponding to
    primes that hash to 'offset', stopping if either
    the list runs out or we find an entry that's
    dedicated to 'prime' already in the table */

    while (offset != 0) {
        entry = table + offset;
        if (entry->prime == prime)
            break;
        offset = entry->next;
    }

    /* if an entry was not found, initialize a new one */

    if (offset == 0) {
        entry = table + new_entry_offset;
        entry->next = hashtable[first_offset];
        entry->prime = prime;
        entry->data = 0;
        entry->count = 0;
        hashtable[first_offset] = new_entry_offset;
    }

    return entry;
}


static uint32_t add_to_hashtable3(qs_cycle_t *table, uint32_t *hashtable,
    uint32_t *primes,
    uint32_t default_table_entry,
    uint32_t *components, uint32_t *vertices) {

    /* update the list of cycles to reflect the presence
    of a partial relation with large primes 0..2

    There are three quantities to track, the number of
    edges, components and vertices. The number of cycles
    in the graph is then e + c - v. There is an edge for
	each partial relation, and one vertex for each prime
    that appears in the graph (these are easy to count).

    A connected component is a group of primes that are
    all reachable from each other via edges in the graph.
    All of the vertices that are in a cycle must belong
    to the same connected component. Think of a component
    as a tree of vertices; one prime is chosen arbitrarily
    as the 'root' of that tree.

    The number of new primes added to the graph (0, 1, 2 or 3)
    is returned 
	
	To handle the TLP case we create a connected component for
	the two smallest primes, and an edge between that
	and the last prime.  Adding a unconnected TLP is then
	cycle-neutral as we get +1 edge, +2 components, 
	and -3 vertices.

	a key realization is that not all cycles need to have 
	primes with even exponents... the matrix step will sort
	all of that out. the cycles just do some pre-combining
	so that the matrix is not too huge.
	*/

    uint32_t root[3];
    uint32_t root1, root2, root3;
    uint32_t i;
    uint32_t num_new_entries = 0;

    /* for each prime */
    for (i = 0; i < 3; i++) {
        uint32_t prime = primes[i];
        uint32_t offset;
        qs_cycle_t *entry;

        /* retrieve the qs_cycle_t corresponding to that
        prime from the graph (or allocate a new one) */

        entry = get_table_entry(table, hashtable, prime, default_table_entry);
        entry->count++;

        if (entry->data == 0) {

            /* if this prime has not occurred in the graph
            before, increment the number of vertices and
            the number of components, then make the table
            entry point to itself. */
            //printf("new vertex: %u\n",prime);
            num_new_entries++;
            default_table_entry++;

            offset = entry - table;
            entry->data = offset;
            (*components)++;
            (*vertices)++;
        }
        else {
            /* the prime is already in the table, which
            means we can follow a linked list of pointers
            to other primes until we reach a qs_cycle_t
            that points to itself. This last qs_cycle_t is
            the 'root' of the connected component that
            contains 'prime'. Save its value */
            //printf("data found: %u\n",prime);

            qs_cycle_t *first_entry, *next_entry;

            first_entry = entry;
            next_entry = table + entry->data;
            while (entry != next_entry) {
                entry = next_entry;
                next_entry = table + next_entry->data;
            }

            /* Also perform path compression: now that we
            know the value of the root for this prime,
            make all of the pointers in the primes we
            visited along the way point back to this root.
            This will speed up future root lookups */

            offset = entry->data;
            entry = first_entry;
            next_entry = table + entry->data;
            while (entry != next_entry) {
                entry->data = offset;
                entry = next_entry;
                next_entry = table + next_entry->data;
            }
        }

        root[i] = offset;
    }

    /* If the roots for the primes are different,
    then they lie within separate connected components.
    We're about to connect this edge (relation) to one of these
    components, and the presence of the other primes
    means that the other components are about to be
    merged together. Hence the total number of components
    in the graph goes down by one. 
	
	
	Here are the possible tlp cases:
	1) three roots
		make a connected component with the smallest prime as root.
		subtract two components for a total of 1.  This along with the
		2 added edges and 3 new vertices gives 0 net additional cycles
	2) two roots
		join the two new vertices to the existing root (component)
		and subtract 2 components for a total of +0 components.
		Along with the 2 added edges and 2 new vertices gives 0
		net additional cycles.
	3) one root
		just do path compression (set all roots to minimum root),
		yielding 2 cycles.
	
	
	
	
	*/

	

    root1 = root[0];
    root2 = root[1];
    root3 = root[2];
    //printf("discovered roots: %u, %u, %u\n", root1, root2, root3);

	if ((root1 != root2) && (root2 != root3) && (root1 != root3))
	{
		// all different roots.  we will combine them all next, so
		// subtract off two of the components
		(*components) -= 2;
	}
	else if ((root1 != root2) || (root2 != root3) || (root1 != root3))
	{
		// two different roots, after combining there will be only one.
		(*components)--;
	}
	else
	{
		// if neither of the above cases are true then we must have
		// that root1 == root2 == root3. 
		// do nothing.
	}

    //printf("num components now %u\n", (*components));

    /* We have to merge any components that are now
    connected by this edge (relation). The surviving component is
    the component whose representative prime is smallest;
    since small primes are more common, this will give
    the smaller root more edges, and will potentially
    increase the number of cycles the graph contains */

    if (table[root1].prime < table[root2].prime)
    {
        if (table[root1].prime < table[root3].prime)
        {
            // r1 < r2 && r1 < r3
            table[root2].data = root1;
            table[root3].data = root1;
        }
        else
        {
            // r1 < r2 && r3 <= r1
            table[root2].data = root3;
            table[root1].data = root3;
        }
    }        
    else
    {
        if (table[root2].prime < table[root3].prime)
        {
            // r2 <= r1 && r2 < r3
            table[root1].data = root2;
            table[root3].data = root2;
        }
        else
        {
            // r2 <= r1 && r3 <= r2
            table[root2].data = root3;
            table[root1].data = root3;
        }
    }


    return num_new_entries;
}

void yafu_add_to_cycles3(static_conf_t *conf, uint32_t flags, uint32_t *primes) {

    /* Top level routine for updating the graph of partial
    relations */

    uint32_t table_size = conf->cycle_table_size;
    uint32_t table_alloc = conf->cycle_table_alloc;
    qs_cycle_t *table = conf->cycle_table;
    uint32_t *hashtable = conf->cycle_hashtable;

    /* if we don't actually want to count cycles,
    just increment the number of vertices. This is
    equivalent to saying that the cycle-counting
    code will never detect any cycles */

    if (flags & MSIEVE_FLAG_SKIP_QS_CYCLES) {
        conf->vertices++;
        return;
    }

    /* make sure there's room for new primes */

    if (table_size + 3 >= table_alloc) {
        table_alloc = conf->cycle_table_alloc = 2 * table_alloc;
        conf->cycle_table = (qs_cycle_t *)xrealloc(conf->cycle_table,
            table_alloc * sizeof(qs_cycle_t));
        table = conf->cycle_table;
    }

    conf->cycle_table_size += add_to_hashtable3(table, hashtable,
        primes,
        table_size,
        &conf->components,
        &conf->vertices);

	if ((primes[0] != primes[1]) && (primes[0] != primes[2]) && (primes[1] != primes[2]))
		conf->num_cycles += 3;
	else if ((primes[0] != primes[1]) || (primes[0] != primes[2]) || (primes[1] != primes[2]))
		conf->num_cycles += 1;

}


/**********************************************************
These 3 routines are used to add a relation to the cycle
tree every time one is found after sieving
(e.g., add_to_cycles is called in save_relation)
**********************************************************/

/*--------------------------------------------------------------------*/
static uint32_t add_to_hashtable(qs_cycle_t *table, uint32_t *hashtable, 
			uint32_t prime1, uint32_t prime2, 
			uint32_t default_table_entry, 
			uint32_t *components, uint32_t *vertices) {

	/* update the list of cycles to reflect the presence
	   of a partial relation with large primes 'prime1'
	   and 'prime2'.

	   There are three quantities to track, the number of
	   edges, components and vertices. The number of cycles 
	   in the graph is then e + c - v. There is one edge for
	   each partial relation, and one vertex for each prime
	   that appears in the graph (these are easy to count).

	   A connected component is a group of primes that are
	   all reachable from each other via edges in the graph.
	   All of the vertices that are in a cycle must belong
	   to the same connected component. Think of a component
	   as a tree of vertices; one prime is chosen arbitrarily
	   as the 'root' of that tree.
	   
	   The number of new primes added to the graph (0, 1, or 2)
	   is returned */

	uint32_t root[2];
	uint32_t root1, root2;
	uint32_t i;
	uint32_t num_new_entries = 0;

	/* for each prime */
	

	for (i = 0; i < 2; i++) {
		uint32_t prime = ((i == 0) ? prime1 : prime2);
		uint32_t offset; 
		qs_cycle_t *entry;

		/* retrieve the qs_cycle_t corresponding to that
		   prime from the graph (or allocate a new one) */

		entry = get_table_entry(table, hashtable,
					prime, default_table_entry);
		entry->count++;

		if (entry->data == 0) {

			/* if this prime has not occurred in the graph
			   before, increment the number of vertices and
			   the number of components, then make the table
			   entry point to itself. */
			//fprintf(globaltest,"new vertex: %u\n",prime);
			num_new_entries++;
			default_table_entry++;

			offset = entry - table;
			entry->data = offset;
			(*components)++;
			(*vertices)++;
		}
		else {
			/* the prime is already in the table, which
			   means we can follow a linked list of pointers
			   to other primes until we reach a qs_cycle_t
			   that points to itself. This last qs_cycle_t is
			   the 'root' of the connected component that
			   contains 'prime'. Save its value */
			//fprintf(globaltest,"data found: %u\n",prime);

			qs_cycle_t *first_entry, *next_entry;

			first_entry = entry;
			next_entry = table + entry->data;
			while (entry != next_entry) {
				entry = next_entry;
				next_entry = table + next_entry->data;
			}
				
			/* Also perform path compression: now that we
			   know the value of the root for this prime,
			   make all of the pointers in the primes we
			   visited along the way point back to this root.
			   This will speed up future root lookups */

			offset = entry->data;
			entry = first_entry;
			next_entry = table + entry->data;
			while (entry != next_entry) {
				entry->data = offset;
				entry = next_entry;
				next_entry = table + next_entry->data;
			}
		}

		root[i] = offset;
	}
				
	/* If the roots for prime1 and prime2 are different,
	   then they lie within separate connected components.
	   We're about to connect this edge to one of these
	   components, and the presence of the other prime
	   means that these two components are about to be
	   merged together. Hence the total number of components
	   in the graph goes down by one. */

	root1 = root[0];
	root2 = root[1];
	if (root1 != root2)
		(*components)--;
	
	/* This partial relation represents an edge in the
	   graph; we have to attach this edge to one or the
	   other of the connected components. Attach it to
	   the component whose representative prime is smallest;
	   since small primes are more common, this will give
	   the smaller root more edges, and will potentially
	   increase the number of cycles the graph contains */

	if (table[root1].prime < table[root2].prime)
		table[root2].data = root1;
	else
		table[root1].data = root2;
	
	return num_new_entries;
}

/*--------------------------------------------------------------------*/
void yafu_add_to_cycles(static_conf_t *conf, uint32_t flags, uint32_t prime1, uint32_t prime2) {

	/* Top level routine for updating the graph of partial
	   relations */

	uint32_t table_size = conf->cycle_table_size;
	uint32_t table_alloc = conf->cycle_table_alloc;
	qs_cycle_t *table = conf->cycle_table;
	uint32_t *hashtable = conf->cycle_hashtable;

	/* if we don't actually want to count cycles,
	   just increment the number of vertices. This is
	   equivalent to saying that the cycle-counting
	   code will never detect any cycles */
	
	if (flags & MSIEVE_FLAG_SKIP_QS_CYCLES) {
		conf->vertices++;
		return;
	}

	/* make sure there's room for new primes */

	if (table_size + 2 >= table_alloc) {
		table_alloc = conf->cycle_table_alloc = 2 * table_alloc;
		conf->cycle_table = (qs_cycle_t *)xrealloc(conf->cycle_table,
						table_alloc * sizeof(qs_cycle_t));
		table = conf->cycle_table;
	}

	conf->cycle_table_size += add_to_hashtable(table, hashtable, 
						prime1, prime2, 
						table_size, 
						&conf->components, 
						&conf->vertices);
}


qs_la_col_t * find_cycles(fact_obj_t*obj, uint32_t *hashtable, qs_cycle_t *table,
	siqs_r *relation_list, uint32_t num_relations, uint32_t *numcycles, uint32_t *numpasses)
{
	qs_la_col_t *cycle_list;
	uint32_t i, start, curr_cycle, passes;
	uint32_t num_cycles = *numcycles;

	/* The idea behind the cycle-finding code is this: the
	   graph is composed of a bunch of connected components,
	   and each component contains one or more cycles. To
	   find the cycles, you build the 'spanning tree' for
	   each component.

	   Think of the spanning tree as a binary tree; there are
	   no cycles in it because leaves are only connected to a
	   common root and not to each other. Any time you connect
	   together two leaves of the tree, though, a cycle is formed.
	   So, for a spanning tree like this:

			 1
			 o
		    / \
		2  o   o  3
		  / \   \
		 o   o   o
		 4   5   6

	   if you connect leaves 4 and 5 you get a cycle (4-2-5). If
	   you connect leaves 4 and 6 you get another cycle (4-2-1-3-6)
	   that will reuse two of the nodes in the first cycle. It's
	   this reuse that makes double large primes so powerful.

	   For our purposes, every edge in the tree above represents
	   a partial relation. Every edge that would create a cycle
	   comes from another partial relation. So to find all the cycles,
	   you begin with the roots of all of the connected components,
	   and then iterate through the list of partial relations until
	   all have been 'processed'. A partial relation is considered
	   processed when one or both of its primes is in the tree. If
	   one prime is present then the relation gets added to the tree;
	   if both primes are present then the relation creates one cycle
	   but is *not* added to the tree.

	   It's really great to see such simple ideas do something so
	   complicated as finding cycles (and doing it very quickly) */

	   /* First traverse the entire graph and remove any vertices
		  that are not the roots of connected components (i.e.
		  remove any primes whose cycle_t entry does not point
		  to itself */

	for (i = 0; i < (1 << QS_LOG2_CYCLE_HASH); i++) {
		uint32_t offset = hashtable[i];

		while (offset != 0) {
			qs_cycle_t *entry = table + offset;

			if (offset != entry->data)
				entry->data = 0;
			offset = entry->next;
		}
	}

	cycle_list = (qs_la_col_t *)xmalloc(num_cycles * sizeof(qs_la_col_t));

	/* keep going until either all cycles are found, all
	   relations are processed, or cycles stop arriving.
	   Normally these conditions all occur at the same time */

	for (start = passes = curr_cycle = 0; start < num_relations &&
		curr_cycle < num_cycles; passes++) {

		/* The list of relations up to index 'start' is con-
		   sidered processed. For all relations past that... */

		uint32_t start_cycles = curr_cycle;

		for (i = start; i < num_relations &&
			curr_cycle < num_cycles; i++) {

			qs_cycle_t *entry1, *entry2;
			siqs_r rtmp = relation_list[i];

			if (rtmp.large_prime[0] == rtmp.large_prime[1]) {

				/* this is a full relation, and forms a
				   cycle just by itself. Move it to position
				   'start' of the relation list and increment
				   'start'. The relation is now frozen at
				   that position */

				qs_la_col_t *c = cycle_list + curr_cycle++;

				if (i != start)
				{
					relation_list[i] = relation_list[start];
				}
				relation_list[start] = rtmp;

				/* build a trivial cycle for the relation */

				c->cycle.num_relations = 1;
				c->cycle.list = (uint32_t *)
					xmalloc(sizeof(uint32_t));
				c->cycle.list[0] = start++;
				continue;
			}

			/* retrieve the cycle_t entries associated
			   with the large primes in relation r. */

			entry1 = get_table_entry(table, hashtable,
				rtmp.large_prime[0], 0);
			entry2 = get_table_entry(table, hashtable,
				rtmp.large_prime[1], 0);

			/* if both vertices do not point to other
			   vertices, then neither prime has been added
			   to the graph yet, and r must remain unprocessed */

			if (entry1->data == 0 && entry2->data == 0)
				continue;

			/* if one or the other prime is part of the
			   graph, add r to the graph. The vertex not in
			   the graph points to the vertex that is, and
			   this entry also points to the relation that
			   is associated with rtmp.

			   If both primes are in the graph, recover the
			   cycle this generates */

			if (entry1->data == 0) {
				entry1->data = entry2 - table;
				entry1->count = start;
			}
			else if (entry2->data == 0) {
				entry2->data = entry1 - table;
				entry2->count = start;
			}
			else {
				int j;

				qs_la_col_t *c = cycle_list + curr_cycle;
				c->cycle.list = NULL;
				qs_enumerate_cycle(obj, c, table, entry1,
					entry2, start);

				//if (c->cycle.num_relations >= 4)
				//	printf("large cycle\n");
				//
				//printf("relations in this cycle:\n");
				//for (j = 0; j < c->cycle.num_relations; j++)
				//{
				//	siqs_r *thisr = &relation_list[c->cycle.list[j]];
				//	printf("poly_idx = %u, offset = %u, lp1 = %u, lp2 = %u\n", 
				//		thisr->poly_idx, thisr->sieve_offset, 
				//		thisr->large_prime[0], thisr->large_prime[1]);
				//}

				if (c->cycle.list)
					curr_cycle++;
			}

			/* whatever happened above, the relation is
			   processed now; move it to position 'start'
			   of the relation list and increment 'start'.
			   The relation is now frozen at that position */

			if (i != start)
			{
				relation_list[i] = relation_list[start];
			}
			relation_list[start++] = rtmp;
		}

		/* If this pass did not find any new cycles, then
		   we've reached steady state and are finished */

		if (curr_cycle == start_cycles)
			break;
	}

	*numcycles = curr_cycle;
	*numpasses = passes;
	return cycle_list;
}

// data structures for hash tables used to build cycles in TLP variation
typedef struct
{
	// relation-by-prime record.
	uint32_t prime;

	// the hashtable is smallish and imperfect, so include a pointer
	// to the next main-table location associated with this hash
	uint32_t next;

	// this record stores the relation indices that contain this prime
	uint32_t num_rids;
	uint32_t alloc_rids;
	uint32_t *rids;
} rbp_t;

typedef struct
{
	// primes-by-relation record
	uint32_t rid;

	// the hashtable is smallish and imperfect, so include a pointer
	// to the next main-table location associated with this hash
	uint32_t next;

	// hold primes appearing in this relation
	uint32_t num_primes;
	uint32_t alloc_primes;
	uint32_t *primes;

	// and relations in this chain
	uint32_t chain_sz;
	uint32_t chain_alloc;
	uint32_t *chain;
} pbr_t;

rbp_t *new_rbp_entry(rbp_t *table, uint32_t *hashtable,
	uint32_t prime, uint32_t *new_entry_offset) {

	/* return a pointer to a unique qs_cycle_t specific
	to 'prime'. The value of 'prime' is hashed and
	the result used to index 'table'. If prime does
	not appear in table, specify 'new_entry_offset'
	as the table entry for prime.

	Values of 'prime' which hash to the same offset
	in hashtable are connected by a linked list of
	offsets in 'table'. */

	uint32_t offset, first_offset;
	rbp_t *entry = NULL;

	first_offset = QS_HASH(prime);
	offset = hashtable[first_offset];

	/* follow the list of entries corresponding to
	primes that hash to 'offset', stopping if either
	the list runs out or we find an entry that's
	dedicated to 'prime' already in the table */

	while (offset != 0) {
		entry = table + offset;
		if (entry->prime == prime)
			break;
		offset = entry->next;
	}

	/* if an entry was not found, initialize a new one */

	if (offset == 0) {
		entry = table + *new_entry_offset;
		entry->next = hashtable[first_offset];
		entry->prime = prime;
		entry->rids = (uint32_t *)xcalloc(4, sizeof(uint32_t));
		entry->alloc_rids = 4;
		entry->num_rids = 0;
		hashtable[first_offset] = *new_entry_offset;
		(*new_entry_offset)++;
	}

	return entry;
}

pbr_t *new_pbr_entry(pbr_t *table, uint32_t *hashtable,
	uint32_t rid, uint32_t *new_entry_offset) {

	/* return a pointer to a unique qs_cycle_t specific
	to 'prime'. The value of 'prime' is hashed and
	the result used to index 'table'. If prime does
	not appear in table, specify 'new_entry_offset'
	as the table entry for prime.

	Values of 'prime' which hash to the same offset
	in hashtable are connected by a linked list of
	offsets in 'table'. */

	uint32_t offset, first_offset;
	pbr_t *entry = NULL;

	first_offset = PBR_HASH(rid);
	offset = hashtable[first_offset];

	/* follow the list of entries corresponding to
	primes that hash to 'offset', stopping if either
	the list runs out or we find an entry that's
	dedicated to 'prime' already in the table */

	while (offset != 0) {
		entry = table + offset;
		if (entry->rid == rid)
			break;
		offset = entry->next;
	}

	/* if an entry was not found, initialize a new one */

	if (offset == 0) {
		entry = table + *new_entry_offset;
		entry->next = hashtable[first_offset];
		entry->rid = rid;
		entry->alloc_primes = 3;
		entry->primes = (uint32_t *)xcalloc(entry->alloc_primes, sizeof(uint32_t));
		entry->num_primes = 0;
		entry->chain = (uint32_t *)xcalloc(2, sizeof(uint32_t));
		entry->chain_alloc = 2;
		entry->chain_sz = 0;
		entry->chain[0] = rid;
		hashtable[first_offset] = *new_entry_offset;
		(*new_entry_offset)++;
	}

	return entry;
}

rbp_t *get_rbp_entry(rbp_t *table, uint32_t *hashtable, uint32_t prime) {

	/* return a pointer to a unique qs_cycle_t specific
	to 'prime'. The value of 'prime' is hashed and
	the result used to index 'table'. If prime does
	not appear in table, specify 'new_entry_offset'
	as the table entry for prime.

	Values of 'prime' which hash to the same offset
	in hashtable are connected by a linked list of
	offsets in 'table'. */

	uint32_t offset, first_offset;
	rbp_t *entry = NULL;

	first_offset = QS_HASH(prime);
	offset = hashtable[first_offset];

	/* follow the list of entries corresponding to
	primes that hash to 'offset', stopping if either
	the list runs out or we find an entry that's
	dedicated to 'prime' already in the table */

	while (offset != 0) {
		entry = table + offset;
		if (entry->prime == prime)
			break;
		offset = entry->next;
	}

	/* if an entry was not found, initialize a new one */

	if (offset == 0) {
		return NULL;
	}

	return entry;
}

pbr_t *get_pbr_entry(pbr_t *table, uint32_t *hashtable, uint32_t rid) {

	/* return a pointer to a unique qs_cycle_t specific
	to 'prime'. The value of 'prime' is hashed and
	the result used to index 'table'. If prime does
	not appear in table, specify 'new_entry_offset'
	as the table entry for prime.

	Values of 'prime' which hash to the same offset
	in hashtable are connected by a linked list of
	offsets in 'table'. */

	uint32_t offset, first_offset;
	pbr_t *entry = NULL;

	first_offset = PBR_HASH(rid);
	offset = hashtable[first_offset];

	/* follow the list of entries corresponding to
	primes that hash to 'offset', stopping if either
	the list runs out or we find an entry that's
	dedicated to 'prime' already in the table */

	while (offset != 0) {
		entry = table + offset;
		if (entry->rid == rid)
			break;
		offset = entry->next;
	}

	/* if an entry was not found, initialize a new one */

	if (offset == 0) {
		return NULL;
	}

	return entry;
}

qs_la_col_t * find_cycles3(fact_obj_t*fobj, static_conf_t *sconf,
	siqs_r *relation_list, uint32_t num_relations, uint32_t *numcycles, uint32_t *numpasses)
{
	// assume that we've kept a backup of all relation data, so feel free to modify it.
	qs_la_col_t *cycle_list;
	uint32_t i, j, k, start, curr_cycle;
	siqs_r *rtmp;
	uint32_t *lp;
	uint32_t *pbr_hashtable;
	uint32_t pbr_table_size;
	uint32_t pbr_table_alloc;
	pbr_t *pbr_table;
	uint32_t *rbp_hashtable;
	uint32_t rbp_table_size;
	uint32_t rbp_table_alloc;
	rbp_t *rbp_table;
	int done;
	uint32_t max_length = 0;
	uint32_t numfull = 0;
	uint32_t cycle_alloc = (num_relations - *numcycles) * 2;
	uint32_t num3lp = 0;

	/*
		Each relation is read in turn and two hash tables built.
		The first, called rbp for "relations-by-primes", contains one
		record per prime, each record keyed by a prime being a linked
		list of relations which contain that prime.
		The other, pbr, contains two records per relation that are
		keyed by the line number of the relation: a linked list of 
		primes appearing in that relation, and a linked list called
		the "chain", initially empty, to hold relations which are in
		the process of being formed into a cycle.
	*/

	pbr_hashtable = (uint32_t *)xcalloc(
		(size_t)(1 << QS_LOG2_CYCLE_HASH),
		sizeof(uint32_t));
	pbr_table_size = 1;
	pbr_table_alloc = num_relations + 100;
	pbr_table = (pbr_t *)xmalloc(
		pbr_table_alloc * sizeof(pbr_t));

	rbp_hashtable = (uint32_t *)xcalloc(
		(size_t)(1 << QS_LOG2_CYCLE_HASH),
		sizeof(uint32_t));
	rbp_table_size = 1;
	rbp_table_alloc = *numcycles + 100;
	rbp_table = (rbp_t *)xmalloc(
		rbp_table_alloc * sizeof(rbp_t));

	// we won't have more than (num_relations - num_primes) cycles
	cycle_list = (qs_la_col_t *)xmalloc(
		cycle_alloc * sizeof(qs_la_col_t));
	curr_cycle = 0;

	printf("commencing rbp/pbr table initialization\n");
	for (i = 1; i < num_relations; i++)
	{
		pbr_t *pbr_entry;
		rtmp = &relation_list[i];
		lp = rtmp->large_prime;

		pbr_entry = new_pbr_entry(pbr_table, pbr_hashtable, i, &pbr_table_size);

		if (pbr_table_size == pbr_table_alloc)
		{
			printf("========== had to reallocate pbr table\n");
			pbr_table_alloc *= 2;
			pbr_table = (pbr_t *)xrealloc(pbr_table, pbr_table_alloc * sizeof(pbr_t));
		}

		// special cases:
		if (lp[1] == lp[2])
		{
			if (lp[0] > 1)
			{
				// list as a single large prime relation
				printf("special case single large prime\n");
				rbp_t *rbp_entry = new_rbp_entry(rbp_table, rbp_hashtable, lp[0], &rbp_table_size);

				rbp_entry->rids[rbp_entry->num_rids] = i;
				rbp_entry->num_rids++;

				pbr_entry->primes[pbr_entry->num_primes] = lp[0];
				pbr_entry->num_primes++;
				continue;
			}
			else
			{
				// list as a full relation
				//printf("special case full\n");

				// build a cycle and don't add the relation to the table.
				qs_la_col_t *c = cycle_list + curr_cycle;

				// we have a cycle
				curr_cycle++;

				if (curr_cycle > cycle_alloc)
					printf("====== cycle alloc failure\n");

				c->cycle.num_relations = 1;
				c->cycle.list = (uint32_t *)xmalloc(c->cycle.num_relations *
					sizeof(uint32_t));

				c->cycle.list[0] = i;
				numfull++;

				continue;
			}
		}

		for (j = 0; j < 3; j++)
		{
			if (lp[j] > 1)
			{
				rbp_t *rbp_entry = new_rbp_entry(rbp_table, rbp_hashtable, lp[j], &rbp_table_size);

				if (rbp_table_size == rbp_table_alloc)
				{
					printf("=========== had to reallocate rbp table\n");
					rbp_table_alloc *= 2;
					rbp_table = (rbp_t *)xrealloc(rbp_table, rbp_table_alloc * sizeof(pbr_t));
				}

				// add this relation to the rbp list
				if (rbp_entry->num_rids == rbp_entry->alloc_rids)
				{
					rbp_entry->alloc_rids *= 2;
					rbp_entry->rids = (uint32_t *)xrealloc(rbp_entry->rids,
						rbp_entry->alloc_rids * sizeof(uint32_t));
				}
				rbp_entry->rids[rbp_entry->num_rids] = i;
				rbp_entry->num_rids++;

				// add this prime to the pbr list
				if (pbr_entry->num_primes == pbr_entry->alloc_primes)
				{
					pbr_entry->alloc_primes *= 2;
					pbr_entry->primes = (uint32_t *)xrealloc(pbr_entry->primes,
						pbr_entry->alloc_primes * sizeof(uint32_t));
				}
				pbr_entry->primes[pbr_entry->num_primes] = lp[j];
				pbr_entry->num_primes++;
			}
		}
	}

	printf("found %u full relations\n", numfull);
	printf("pbr table size = %u\n", pbr_table_size);
	printf("rbp table size = %u\n", rbp_table_size);

	/*
	Once all relations have been read and the hash tables built, we begin
	growing chains of relations until a cycle is formed, which is then emitted as
	a line containing the line numbers of the relations concerned.  Repeated linear
	sweeps through the pbr table are made: if the referenced relation r0 is a par, i.e.,
	its list of primes consists of a single element p, the list of relations containing
	p is retrieved from the rbp table. Each relation ri in the list, other than r0, is
	dealt with in turn.  If the relation ri is a par, the pair form a cycle; r0, ri and
	their respective chains (if non-empty) are emitted.  Otherwise, the chain of r0
	and r0 itself is appended to the chain of ri, and the prime p is deleted from the 
	prime-list in ri. When all the list of relations containing p has been processed 
	in this manner, the entry keyed by r0 is deleted from pbr and the entry for r0 
	keyed by p is deleted from rbp.

	The above procedure is repeated until no further changes to the hash tables are made.
	*/
	
	done = 0;
	*numpasses = 0;
	
	while (!done)
	{
		printf("commencing cycle formation pass %d\n", *numpasses + 1);
		done = 1;
		for (i = 1; i < pbr_table_size; i++)
		{
			pbr_t *pbr_entry = pbr_table + i;
			rbp_t *rbp_entry;
			pbr_t *pbr_entry2;

			if (pbr_entry->num_primes == 1)
			{
				if ((rbp_entry = get_rbp_entry(rbp_table, rbp_hashtable,
					pbr_entry->primes[0])) == NULL)
				{
					printf("couldn't find rbp entry associated with prime %u\n",
						pbr_entry->primes[0]);
					exit(1);
				}

				//printf("rbp entry for prime %u of relation %d lists %u relations\n",
				//	pbr_entry->primes[0], i, rbp_entry->num_rids);

				for (j = 0; j < rbp_entry->num_rids; j++)
				{
					if (rbp_entry->rids[j] == i)
						continue;

					if ((pbr_entry2 = get_pbr_entry(pbr_table, pbr_hashtable,
						rbp_entry->rids[j])) == NULL)
					{
						printf("couldn't find pbr entry associated with rid %u\n",
							rbp_entry->rids[j]);
						printf("rid %d lp's are: %u,%u,%u\n", rbp_entry->rids[j],
							relation_list[rbp_entry->rids[j]].large_prime[0],
							relation_list[rbp_entry->rids[j]].large_prime[1],
							relation_list[rbp_entry->rids[j]].large_prime[2]);

						exit(1);

						printf("pbr %d p's are: ", rbp_entry->rids[j]);
						for (k = 0; k < pbr_table[rbp_entry->rids[j]].num_primes; k++)
						{
							printf("%u ", pbr_table[rbp_entry->rids[j]].primes[k]);
						}
						printf("\n");
						printf("HASH(%u) = %u\n", rbp_entry->rids[j], PBR_HASH(rbp_entry->rids[j]));
						printf("pbr_hashtable[%u] = %u\n", PBR_HASH(rbp_entry->rids[j]),
							pbr_hashtable[PBR_HASH(rbp_entry->rids[j])]);
						exit(1);
					}

					if (pbr_entry2->num_primes == 1)
					{
						uint32_t rid;
						uint32_t length = pbr_entry->chain_sz +
							pbr_entry2->chain_sz + 2;
						qs_la_col_t *c = cycle_list + curr_cycle;
						uint32_t m = 0;
						int has3lp = 0;

						if (length > max_length)
							max_length = length;

						// we have a cycle
						curr_cycle++;

						if (curr_cycle > cycle_alloc)
							printf("====== cycle alloc failure\n");
						
						/* Now that we know how many relations are in the
						cycle, allocate space to remember them */

						c->cycle.num_relations = length;
						c->cycle.list = (uint32_t *)xmalloc(c->cycle.num_relations *
							sizeof(uint32_t));

						rid = pbr_entry2->rid;
						c->cycle.list[m++] = i;
						c->cycle.list[m++] = rid;

						if ((relation_list[i].large_prime[0] > 1) &&
							(relation_list[i].large_prime[1] > 1) &&
							(relation_list[i].large_prime[2] > 1) &&
							(has3lp == 0))
						{
							num3lp++;
							has3lp = 1;
						}

						if ((relation_list[rid].large_prime[0] > 1) &&
							(relation_list[rid].large_prime[1] > 1) &&
							(relation_list[rid].large_prime[2] > 1) &&
							(has3lp == 0))
						{
							num3lp++;
							has3lp = 1;
						}

						if (fobj->VFLAG > 2)
						{
							printf("found cycle of length %d\n", length);

							printf("relation %08d: %u,%u,%u\n", i,
								relation_list[i].large_prime[0],
								relation_list[i].large_prime[1],
								relation_list[i].large_prime[2]);
							
							printf("relation %08d: %u,%u,%u\n", rid,
								relation_list[rid].large_prime[0],
								relation_list[rid].large_prime[1],
								relation_list[rid].large_prime[2]);
						}

						for (k = 0; k < pbr_entry->chain_sz; k++)
						{
							rid = pbr_entry->chain[k];
							c->cycle.list[m++] = rid;
							if ((relation_list[rid].large_prime[0] > 1) &&
								(relation_list[rid].large_prime[1] > 1) &&
								(relation_list[rid].large_prime[2] > 1) &&
								(has3lp == 0))
							{
								num3lp++;
								has3lp = 1;
							}

							if (fobj->VFLAG > 2)
							{
								printf("relation %08d: %u,%u,%u\n", rid,
									relation_list[rid].large_prime[0],
									relation_list[rid].large_prime[1],
									relation_list[rid].large_prime[2]);
							}
						}

						for (k = 0; k < pbr_entry2->chain_sz; k++)
						{
							rid = pbr_entry2->chain[k];
							c->cycle.list[m++] = rid;
							if ((relation_list[rid].large_prime[0] > 1) &&
								(relation_list[rid].large_prime[1] > 1) &&
								(relation_list[rid].large_prime[2] > 1) && 
								(has3lp == 0))
							{
								num3lp++;
								has3lp = 1;
							}

							if (fobj->VFLAG > 2)
							{
								printf("relation %08d: %u,%u,%u\n", rid,
									relation_list[rid].large_prime[0],
									relation_list[rid].large_prime[1],
									relation_list[rid].large_prime[2]);
							}
						}


					}
					else
					{
						// Otherwise, the chain of r0
						// and r0 itself is appended to the chain of ri, 
						// and the prime p is deleted from the prime-list in ri.
						uint32_t offset = pbr_entry2->chain_sz;

						if ((offset + pbr_entry->chain_sz + 1) >= pbr_entry2->chain_alloc)
						{
							pbr_entry2->chain_alloc = (pbr_entry2->chain_alloc + pbr_entry->chain_sz + 2);
							pbr_entry2->chain = (uint32_t *)xrealloc(pbr_entry2->chain,
								pbr_entry2->chain_alloc * sizeof(uint32_t));
						}

						for (k = 0; k < pbr_entry->chain_sz; k++)
						{
							pbr_entry2->chain[offset + k] = pbr_entry->chain[k];
						}
						pbr_entry2->chain[offset + k] = i;
						pbr_entry2->chain_sz += (pbr_entry->chain_sz + 1);

						for (k = 0; k < pbr_entry2->num_primes; k++)
						{
							if (pbr_entry2->primes[k] == pbr_entry->primes[0])
								break;
						}

						if (k == pbr_entry2->num_primes)
						{
							printf("error: didn't find prime %u in pbr_entry2 prime list\n",
								pbr_entry->primes[0]);
							exit(1);
						}

						for (; k < pbr_entry2->num_primes - 1; k++)
						{
							pbr_entry2->primes[k] = pbr_entry2->primes[k + 1];
						}
						pbr_entry2->num_primes--;
					}
				}

				// delete pbr[r0] and rbp[p], or ensure we never visit them again.
				// if linear traversal speed becomes an issue, can probably make
				// some sort of list pointer so we can skip them.
				pbr_entry->num_primes = 0;
				rbp_entry->num_rids = 0;
				done = 0;
			}
		}

		(*numpasses)++;
	}

	sconf->num_r = curr_cycle;
	sconf->num_relations = numfull;
	sconf->num_cycles = num_relations;

	printf("expected %u cycles, found %u cycles\n", 
		pbr_table_size - rbp_table_size, curr_cycle);
	printf("found %u cycles from partial relations\n", curr_cycle - numfull);
	printf("maximum cycle length = %u\n", max_length);
	printf("%1.1f%% of cycles from partials involve a tlp\n", (double)num3lp / (double)(curr_cycle - numfull) * 100.0);

	if (fobj->logfile != NULL)
	{
		logprint(fobj->logfile, "expected %u cycles, found %u cycles\n",
			pbr_table_size - rbp_table_size, curr_cycle);
		logprint(fobj->logfile, "found %u cycles from partial relations\n", curr_cycle - numfull);
		logprint(fobj->logfile, "maximum cycle length = %u\n", max_length);
		logprint(fobj->logfile, "%1.1f%% of cycles from partials involve a tlp\n", 
			(double)num3lp / (double)(curr_cycle - numfull) * 100.0);

	}

	for (i = 1; i < pbr_table_size; i++)
	{
		free(pbr_table[i].chain);
		free(pbr_table[i].primes);
	}

	for (i = 1; i < rbp_table_size; i++)
	{
		free(rbp_table[i].rids);
	}

	free(rbp_hashtable);
	free(rbp_table);
	free(pbr_hashtable);
	free(pbr_table);

	*numcycles = curr_cycle;
	return cycle_list;
}

/*******************************************************************************
These functions are used after sieving is complete to read in all
relations and find/optimize all the cycles
*******************************************************************************/
#define NUM_CYCLE_BINS 8

void yafu_qs_filter_relations(static_conf_t *sconf) {

	/* Perform all of the postprocessing on the list
	   of relations from the sieving phase. There are
	   two main jobs, reading in all the relations that
	   will be used and then determining the list of 
	   cycles in which partial relations appear. Care
	   should be taken to avoid wasting huge amounts of
	   memory */

	fact_obj_t *fobj = sconf->obj;
	uint32_t *hashtable = sconf->cycle_hashtable;
	qs_cycle_t *table = sconf->cycle_table;
	uint32_t num_derived_poly;
	uint32_t *final_poly_index;
	uint32_t num_relations, num_cycles, num_poly;
	qs_la_col_t *cycle_list;
	siqs_r *relation_list;

	uint32_t i, passes, start;
	uint32_t curr_a_idx, curr_poly_idx, curr_rel; 
	uint32_t curr_expected, curr_saved, curr_cycle; 
	uint32_t total_poly_a;
	uint32_t poly_saved;
	uint32_t cycle_bins[NUM_CYCLE_BINS+1] = {0};
	char buf[LINE_BUF_SIZE];
	char *subbuf;
	int first, last_poly;
	uint32_t this_rel = 0;

 	/* Rather than reading all the relations in and 
	   then removing singletons, read only the large 
	   primes of each relation into an initial list,
	   remove the singletons, and then only read in
	   the relations that survive. This avoids reading
	   in useless relations (and usually the polynomials 
	   they would need) */

	i = 0;
	total_poly_a = 0;

	if (!sconf->in_mem)
	{
		/* skip over the first line */
		qs_savefile_open(&fobj->qs_obj.savefile, SAVEFILE_READ);
		qs_savefile_read_line(buf, sizeof(buf), &fobj->qs_obj.savefile);

		//we don't know beforehand how many rels to expect, so start
		//with some amount and allow it to increase as we read them
		relation_list = (siqs_r *)xmalloc(10000 * sizeof(siqs_r));
		curr_rel = 10000;
		while (!qs_savefile_eof(&fobj->qs_obj.savefile)) {
			char *start;

			switch (buf[0]) {
			case 'A':
				total_poly_a++;
				break;

			case 'R':
				start = strchr(buf, 'L');
				if (start != NULL) {

					if (i == curr_rel) {
						curr_rel = 3 * curr_rel / 2;
						relation_list = (siqs_r *)xrealloc(
							relation_list,
							curr_rel *
							sizeof(siqs_r));
					}

					if (sconf->use_dlp == 2)
					{ 
						uint32_t primes[3];
						yafu_read_tlp(start, primes);

						relation_list[i].poly_idx = i;
						relation_list[i].large_prime[0] = primes[0];
						relation_list[i].large_prime[1] = primes[1];
						relation_list[i].large_prime[2] = primes[2];
					}
					else
					{
						uint32_t prime1, prime2;
						yafu_read_large_primes(start, &prime1, &prime2);
						
						relation_list[i].poly_idx = i;
						relation_list[i].large_prime[0] = prime1;
						relation_list[i].large_prime[1] = prime2;
					}
                    relation_list[i].apoly_idx = total_poly_a - 1;
					i++;
				}
				break;
			}

			qs_savefile_read_line(buf, sizeof(buf), &fobj->qs_obj.savefile);
		}
		num_relations = i;
	}
	else
	{
		relation_list = (siqs_r *)xmalloc(sconf->buffered_rels * sizeof(siqs_r));
		for (i=0; i<sconf->buffered_rels; i++)
		{
            relation_list[i].poly_idx = i;
            relation_list[i].apoly_idx = sconf->in_mem_relations[i].apoly_idx;
			relation_list[i].large_prime[0] = sconf->in_mem_relations[i].large_prime[0];
			relation_list[i].large_prime[1] = sconf->in_mem_relations[i].large_prime[1];
            relation_list[i].large_prime[2] = sconf->in_mem_relations[i].large_prime[2];
		}
		total_poly_a = sconf->total_poly_a;
		num_relations = sconf->buffered_rels;
	}

    if (fobj->VFLAG > 0)
    {
        printf("read %d relations\n", num_relations);
    }

	// re-filtering
    if (sconf->charcount == 42)
    {
        rebuild_graph(sconf, relation_list, num_relations);
    }

	if (sconf->use_dlp == 2)
	{
		// tlp variation may not have built graph as it progressed... do it now.
		rebuild_graph3(sconf, relation_list, num_relations);

		num_relations = qs_purge_singletons3(fobj, relation_list, num_relations,
			table, hashtable);

		// the loop below will rebuild the graph again with relations that
		// survived singleton removal.  Clear the tables.
		memset(sconf->cycle_hashtable, 0, sizeof(uint32_t) * (1 << QS_LOG2_CYCLE_HASH));
		sconf->num_cycles = 0;
		sconf->vertices = 0;
		sconf->components = 0;
		sconf->cycle_table_size = 1;
		memset(sconf->cycle_table, 0, sconf->cycle_table_alloc * sizeof(qs_cycle_t));
	}
	else
	{
		num_relations = qs_purge_singletons(fobj, relation_list, num_relations,
			table, hashtable);
	}

    //printf("relations surviving singleton removal:\n");
    //for (i = 0; i < num_relations; i++)
    //{
    //    printf("pos: %u, polya = %u, lp = %x,%x\n", relation_list[i].poly_idx,
    //        relation_list[i].apoly_idx, relation_list[i].large_prime[0],
    //        relation_list[i].large_prime[1]);
    //
    //}

	relation_list = (siqs_r *)xrealloc(relation_list, num_relations * 
							sizeof(siqs_r));

	/* Now we know how many relations to expect. Also
	   initialize the lists of polynomial 'a' and 'b' values */

	num_poly = 10000;

	sconf->poly_list = (poly_t *)xmalloc(num_poly * sizeof(poly_t));

	// if in-mem, don't need a new a_list
	if (!sconf->in_mem)
	{		
		sconf->total_poly_a = total_poly_a;
		sconf->poly_a_list = (mpz_t *)xmalloc(total_poly_a * sizeof(mpz_t));
		for (i=0; i<total_poly_a; i++)
			mpz_init(sconf->poly_a_list[i]);
	}

	final_poly_index = (uint32_t *)xmalloc(1024 * 
						sizeof(uint32_t));
	
	/* initialize the running counts of relations and
	   polynomials */

	i = 0;
	last_poly = -1;
	curr_expected = 0;
	curr_saved = 0;
	curr_rel = (uint32_t)(-1);
	curr_poly_idx = (uint32_t)(-1);
	curr_a_idx = (uint32_t)(-1);
	poly_saved = 0;
	sconf->poly_list_alloc = 0;
	if (fobj->VFLAG > 0)
		printf("attempting to read %u relations\n", num_relations);

	/* Read in the relations and the polynomials they use
	   at the same time. */

	if (!sconf->in_mem)
		qs_savefile_rewind(&fobj->qs_obj.savefile);
	else
		this_rel = 0;

	first = 1;
	while (curr_expected < num_relations) {
		char *tmp;
		uint32_t bad_A_val = 0;
		siqs_r *r;
		siqs_r *rel = NULL;

		/* read in the next entity */
		if (!sconf->in_mem)
		{
			if (qs_savefile_eof(&fobj->qs_obj.savefile))
				break;
			qs_savefile_read_line(buf, sizeof(buf), &fobj->qs_obj.savefile);
		}
		else
		{
            if (this_rel == sconf->buffered_rels)
            {
                printf("end of in-memory relations\n");
                break;
            }

			rel = sconf->in_mem_relations + this_rel++;

            //printf("this rel: %u of %u, curr_rel: %u, curr_expected: %u, "
            //    "curr_saved: %u, num_relations: %u, poly a,b = %u,%u\n",
            //    this_rel - 1, sconf->buffered_rels, curr_rel,
            //    curr_expected, curr_saved, num_relations, rel->poly_idx, rel->apoly_idx);

            if (rel->apoly_idx != last_poly)
                buf[0] = 'A';
            else
                buf[0] = 'R';
		}

		switch (buf[0]) {
		case 'A':
			/* Read in a new 'a' value */
			/* build all of the 'b' values associated with it */
			if (!sconf->in_mem)
			{
				subbuf = buf + 2;	//skip the A and a space
				mpz_set_str(sconf->curr_a, subbuf, 0);
                curr_a_idx++;
                //printf("New poly_a = %s, index %u\n", subbuf, curr_a_idx + 1);
			}
			else
			{
                //gmp_printf("New poly_a = %Zx, index %u of %u\n", 
                //    sconf->poly_a_list[rel->apoly_idx],
                //    rel->apoly_idx, sconf->total_poly_a);
                last_poly = rel->apoly_idx;
				this_rel--;
				mpz_set(sconf->curr_a, sconf->poly_a_list[rel->apoly_idx]);
                curr_a_idx = rel->apoly_idx;
			}

            
			num_derived_poly = process_poly_a(sconf);

			if (num_derived_poly == 0)
			{
				//this is an error indicating a bad poly a.  skip all relations
				//until we see the next A
				bad_A_val = 1;
				continue;
			}
			else
				bad_A_val = 0;

            if (!sconf->in_mem)
			    mpz_set(sconf->poly_a_list[curr_a_idx], sconf->curr_a); 

			/* all 'b' values start off unused */
			final_poly_index = (uint32_t *)xrealloc(final_poly_index,
				num_derived_poly * sizeof(uint32_t));
			memset(final_poly_index, -1, num_derived_poly *
							sizeof(uint32_t));

			break;

		case 'R':
				/* handle a new relation. First find the 
			   large primes; these will determine
	     		   if a relation is full or partial */
			if (bad_A_val)
				continue;

			// corrupted rel?
			if (!sconf->in_mem)
			{
				tmp = strchr(buf, 'L');
				if (tmp == NULL)
					break;
			}
			else
			{
				if (rel->large_prime[0] == 0)
					break;
			}

			/* Check if this relation is needed. If it
			   survived singleton removal then its 
			   ordinal ID will be in the next entry 
			   of relation_list. 
		   
			   First move up the index of relation_list 
			   until the relation index to check is >= 
			   the one we have (it may have gotten behind 
			   because relations were corrupted) */
			
			curr_rel++;

			while (curr_expected < num_relations &&
				relation_list[curr_expected].poly_idx <
						curr_rel) {
				curr_expected++;
			}

			/* now check if the relation should be saved */

			if (curr_expected >= num_relations ||
			    relation_list[curr_expected].poly_idx != curr_rel)
				break;

			curr_expected++;

			if (!sconf->in_mem)
			{
				/* convert the ASCII text of the relation to a
				relation_t, verifying correctness in the process */
				r = relation_list + curr_saved;
				subbuf = buf + 2;	//skip over the R and a space

                //printf("saving buffered rel @ pos %u: %s",
                //    curr_saved, subbuf);

				if (process_rel(subbuf, sconf->factor_base,
					sconf->n, sconf, sconf->obj, r)) {

						if (fobj->logfile != NULL)
							logprint(fobj->logfile, "failed to read relation %d\n", 
								curr_expected - 1);
						if (fobj->VFLAG > 1)
							printf("failed to read relation %d\n", 
								curr_expected - 1);
					break;
				}
				// process_rel calls add_to_hashtable, which can grow the size of the cycle table.
				// update our local pointer to the table every time in case this happens.
				table = sconf->cycle_table;
			}
			else
			{
                r = relation_list + curr_saved;

                //char c = ' ';
                //
                //if (rel->parity)
                //    c = '-';
                //printf("saving in-mem rel %u @ pos %u, %c%x %u,%u L %x %x\n", 
                //    this_rel-1, curr_saved, c, rel->sieve_offset, rel->apoly_idx,
                //    rel->poly_idx, rel->large_prime[0], rel->large_prime[1]);

                td_and_merge_relation(sconf->factor_base,
                    sconf->n, sconf, sconf->obj, r, rel);

                table = sconf->cycle_table;
			}

			curr_saved++;

			/* if necessary, save the b value corresponding 
			   to this relation */

			if (final_poly_index[r->poly_idx] == (uint32_t)(-1)) {

				if (i == num_poly) {
					num_poly *= 2;
					sconf->poly_list = (poly_t *)xrealloc(
							sconf->poly_list,
							num_poly *
							sizeof(poly_t));
				}

				sconf->poly_list[i].a_idx = curr_a_idx;
				mpz_init(sconf->poly_list[i].b);
				mpz_set(sconf->poly_list[i].b, sconf->curr_b[r->poly_idx]);
				sconf->poly_list_alloc++;
				final_poly_index[r->poly_idx] = i;
				r->poly_idx = i++;
			}
			else {
				r->poly_idx = final_poly_index[r->poly_idx];
			}

			break;  /* done with this relation */

		}
	}

	/* update the structures with the counts of relations
	   and polynomials actually recovered */

	num_relations = curr_saved;
	if (fobj->logfile != NULL)
	{
		logprint(fobj->logfile, "recovered %u relations\n", num_relations);
		logprint(fobj->logfile, "recovered %u polynomials\n", i);
	}
	if (fobj->VFLAG > 0)
	{
		printf("recovered %u relations\n", num_relations);
		printf("recovered %u polynomials\n", i);
	}

	if (!sconf->in_mem)
		qs_savefile_close(&fobj->qs_obj.savefile);

	free(final_poly_index);
	sconf->poly_list = (poly_t *)xrealloc(sconf->poly_list,
					   i * sizeof(poly_t));

	/* begin the cycle generation process by purging
	   duplicate relations. For the sake of consistency, 
	   always rebuild the graph afterwards */

	if (sconf->use_dlp == 2)
	{
		num_relations = qs_purge_duplicate_relations3(fobj,
			relation_list, num_relations);

		num_cycles = sconf->vertices;

		memset(sconf->cycle_hashtable, 0, sizeof(uint32_t) * (1 << QS_LOG2_CYCLE_HASH));
		sconf->vertices = 0;
		sconf->components = 0;
		sconf->cycle_table_size = 1;
		memset(sconf->cycle_table, 0, sconf->cycle_table_alloc * sizeof(qs_cycle_t));

		rebuild_graph3(sconf, relation_list, num_relations);

		cycle_list = find_cycles3(fobj, sconf, relation_list,
			num_relations, &num_cycles, &passes);
	}
	else
	{
		num_relations = qs_purge_duplicate_relations(fobj,
			relation_list, num_relations);

		memset(hashtable, 0, sizeof(uint32_t) * (1 << QS_LOG2_CYCLE_HASH));
		sconf->vertices = 0;
		sconf->components = 0;
		sconf->cycle_table_size = 1;

		rebuild_graph(sconf, relation_list, num_relations);

		/* compute the number of cycles to expect. Note that
		   this number includes cycles from both full and partial
		   relations (the cycle for a full relation is trivial) */

		num_cycles = num_relations + sconf->components - sconf->vertices;

		if (fobj->logfile != NULL)
			logprint(fobj->logfile, "attempting to build %u cycles\n", num_cycles);
		if (fobj->VFLAG > 0)
		{
			printf("attempting to build %u cycles\n", num_cycles);
			fflush(stdout);
		}

		cycle_list = find_cycles(fobj, hashtable, table, relation_list,
			num_relations, &num_cycles, &passes);
	}

	if (fobj->logfile != NULL)
		logprint(fobj->logfile, "found %u cycles from %u relations in %u passes\n", 
			num_cycles, num_relations, passes);
	if (fobj->VFLAG > 0)
		printf("found %u cycles from %u relations in %u passes\n", 
			num_cycles, num_relations, passes);
	
	/* sort the list of cycles so that the cycles with
	   the largest number of relations will come last. 
	   If the linear algebra code skips any cycles it
	   can easily skip the most dense cycles */

	qsort(cycle_list, (size_t)num_cycles, sizeof(qs_la_col_t), yafu_sort_cycles);

	sconf->relation_list = relation_list;
	sconf->num_relations = num_relations;
	sconf->cycle_list = cycle_list;
	sconf->num_cycles = num_cycles;

	/* print out a histogram of cycle lengths for infor-
	   mational purposes */

	for (i = 0; i < num_cycles; i++) {
		num_relations = cycle_list[i].cycle.num_relations;

		if (num_relations >= NUM_CYCLE_BINS)
			cycle_bins[NUM_CYCLE_BINS]++;
		else
			cycle_bins[num_relations - 1]++;
	}

	if (fobj->logfile != NULL)
		logprint(fobj->logfile, "distribution of cycle lengths:\n");

	if (fobj->VFLAG > 0)
		printf("distribution of cycle lengths:\n");
	for (i = 0; i < NUM_CYCLE_BINS; i++) {
		if (cycle_bins[i]) {
			if (fobj->logfile != NULL)
				logprint(fobj->logfile, "   length %d : %d\n", 
					i + 1, cycle_bins[i]);
			if (fobj->VFLAG > 0)
				printf("   length %d : %d\n", 
					i + 1, cycle_bins[i]);
		}
	}
	if (cycle_bins[i])
	{
		if (fobj->logfile != NULL)
			logprint(fobj->logfile, "   length %u+: %u\n", i + 1, cycle_bins[i]);

		if (fobj->VFLAG > 0)
			printf("   length %u+: %u\n", i + 1, cycle_bins[i]);
	}
	
	if (fobj->logfile != NULL)
		logprint(fobj->logfile, "largest cycle: %u relations\n",
			cycle_list[num_cycles-1].cycle.num_relations);

	if (fobj->VFLAG > 0)
		printf("largest cycle: %u relations\n",
			cycle_list[num_cycles-1].cycle.num_relations);
}

static int compare_relations(const void *x, const void *y) {

	/* Callback used to sort a list of sieve relations.
	   Sorting is by size of large primes, then by number
	   of factors, then by factor values. Only the first
	   rule is needed for ordinary MPQS, but with self-
	   initialization we have to detect duplicate relations,
	   and this is easier if they are sorted as described */

	siqs_r *xx = (siqs_r *)x;
	siqs_r *yy = (siqs_r *)y;
	uint32_t i;

	if (xx->large_prime[1] > yy->large_prime[1])
		return 1;
	if (xx->large_prime[1] < yy->large_prime[1])
		return -1;

	if (xx->large_prime[0] > yy->large_prime[0])
		return 1;
	if (xx->large_prime[0] < yy->large_prime[0])
		return -1;

	if (xx->num_factors > yy->num_factors)
		return 1;
	if (xx->num_factors < yy->num_factors)
		return -1;

	for (i = 0; i < xx->num_factors; i++) {
		if (xx->fb_offsets[i] > yy->fb_offsets[i])
			return 1;
		if (xx->fb_offsets[i] < yy->fb_offsets[i])
			return -1;
	}
	return 0;
}

static int compare_relations3(const void *x, const void *y) {

	/* Callback used to sort a list of sieve relations.
	   Sorting is by size of large primes, then by number
	   of factors, then by factor values. Only the first
	   rule is needed for ordinary MPQS, but with self-
	   initialization we have to detect duplicate relations,
	   and this is easier if they are sorted as described */

	siqs_r *xx = (siqs_r *)x;
	siqs_r *yy = (siqs_r *)y;
	uint32_t i;

	if (xx->large_prime[2] > yy->large_prime[2])
		return 1;
	if (xx->large_prime[2] < yy->large_prime[2])
		return -1;

	if (xx->large_prime[1] > yy->large_prime[1])
		return 1;
	if (xx->large_prime[1] < yy->large_prime[1])
		return -1;

	if (xx->large_prime[0] > yy->large_prime[0])
		return 1;
	if (xx->large_prime[0] < yy->large_prime[0])
		return -1;

	if (xx->num_factors > yy->num_factors)
		return 1;
	if (xx->num_factors < yy->num_factors)
		return -1;

	for (i = 0; i < xx->num_factors; i++) {
		if (xx->fb_offsets[i] > yy->fb_offsets[i])
			return 1;
		if (xx->fb_offsets[i] < yy->fb_offsets[i])
			return -1;
	}
	return 0;
}

/*--------------------------------------------------------------------*/
uint32_t qs_purge_duplicate_relations(fact_obj_t*fobj,
				siqs_r *rlist, 
				uint32_t num_relations) {

	uint32_t i, j;
	
	/* remove duplicates from rlist */

	if (num_relations < 2)
		return num_relations;

	qsort(rlist, (size_t)num_relations,
		sizeof(siqs_r), compare_relations);

	for (i = 1, j = 0; i < num_relations; i++) {
		if (compare_relations(rlist + j, rlist + i) == 0)
		{
			printf("relations %d and %d with polyidx %d,%d and polyidx %d,%d are the same\n", 
				j, i, rlist[j].apoly_idx, rlist[j].poly_idx, 
                rlist[i].apoly_idx, rlist[i].poly_idx);
        }
        else
        {
            j++;
            if (j != i)
            {
                rlist[j] = rlist[i];
            }
        }
	}

	j++;
	if (j != num_relations)
	{
		if (fobj->logfile != NULL)
			logprint(fobj->logfile, "freed %d duplicate relations\n", 
					num_relations - j);
		if (fobj->VFLAG > 0)
			printf("freed %d duplicate relations\n", 
					num_relations - j);
	}

	return j;
}

uint32_t qs_purge_duplicate_relations3(fact_obj_t*obj,
	siqs_r *rlist,
	uint32_t num_relations) {

	uint32_t i, j;

	/* remove duplicates from rlist */

	if (num_relations < 2)
		return num_relations;

	qsort(rlist, (size_t)num_relations,
		sizeof(siqs_r), compare_relations);

	for (i = 1, j = 0; i < num_relations; i++) {
		if (compare_relations3(rlist + j, rlist + i) == 0)
		{
			printf("relations with polyidx %d and polyidx %d are the same\n",
				rlist[j].poly_idx, rlist[i].poly_idx);
		}
		else
		{
			j++;
			if (j != i)
			{
				rlist[j] = rlist[i];
			}
		}
	}

	j++;
	if (j != num_relations)
	{
		if (obj->logfile != NULL)
			logprint(obj->logfile, "freed %d duplicate relations\n",
				num_relations - j);
		if (obj->VFLAG > 0)
			printf("freed %d duplicate relations\n",
				num_relations - j);
	}

	return j;
}

void yafu_read_large_primes(char *buf, uint32_t *prime1, uint32_t *prime2) {

	char *next_field;
	uint32_t p1, p2;

	*prime1 = p1 = 1;
	*prime2 = p2 = 2;
	if (*buf != 'L')
		return;

	buf++;
	while (isspace(*buf))
		buf++;
	if (isxdigit(*buf)) {
		p1 = strtoul(buf, &next_field, 16);
		buf = next_field;
	}
	else {
		return;
	}

	while (isspace(*buf))
		buf++;
	if (isxdigit(*buf))
		p2 = strtoul(buf, &next_field, 16);
	
	if (p1 < p2) {
		*prime1 = p1;
		*prime2 = p2;
	}
	else {
		*prime1 = p2;
		*prime2 = p1;
	}
}

void yafu_read_tlp(char *buf, uint32_t *primes) {

	char *next_field;
	uint32_t p1, p2, p3;

	primes[0] = p1 = 1;
	primes[1] = p2 = 2;
	primes[2] = p3 = 3;

	if (*buf != 'L')
	{
		printf("input buf doesn't start with 'L'\n");
		return;
	}

	buf++;
	while (isspace(*buf))
		buf++;
	if (isxdigit(*buf)) {
		p1 = strtoul(buf, &next_field, 16);
		buf = next_field;
	}
	else {
		printf("field %c not an xdigit\n", *buf);
		return;
	}

	while (isspace(*buf))
		buf++;
	if (isxdigit(*buf)) {
		p2 = strtoul(buf, &next_field, 16);
		buf = next_field;
	}
	else {
		printf("field %c not an xdigit\n", *buf);
		return;
	}

	while (isspace(*buf))
		buf++;
	if (isxdigit(*buf))
		p3 = strtoul(buf, &next_field, 16);

	primes[0] = p1;
	primes[1] = p2;
	primes[2] = p3;

	return;
}

uint32_t qs_purge_singletons(fact_obj_t*fobj, siqs_r *list,
				uint32_t num_relations,
				qs_cycle_t *table, uint32_t *hashtable) {
	
	/* given a list of relations and the graph from the
	   sieving stage, remove any relation that contains
	   a prime that only occurs once in the graph. Because
	   removing a relation removes a second prime as well,
	   this process must be iterated until no more relations
	   are removed */

	uint32_t num_left;
	uint32_t i, j, k;
	uint32_t passes = 0;

	if (fobj->VFLAG > 0)
		printf("begin singleton removal with %u relations\n", num_relations);
	if (fobj->logfile != NULL)
		logprint(fobj->logfile, "begin singleton removal with %u relations\n", num_relations);

	do {
		num_left = num_relations;

		/* for each relation */
		//printf("now at pass %d\n", passes);

		for (i = j = 0; i < num_relations; i++) {
			siqs_r *r = list + i;
			uint32_t prime;
			qs_cycle_t *entry;

			/* full relations always survive */

			if (r->large_prime[0] == r->large_prime[1]) {
                if (j != i)
                {
                    list[j] = list[i];
                }
                j++;
				continue;
			}

			/* for each prime in that relation */

			for (k = 0; k < 2; k++) {
				prime = r->large_prime[k];
				entry = get_table_entry(table, hashtable,
							prime, 0);

				/* if the relation is due to be removed,
				   decrement the count of its other
				   primes in the graph. The following is
				   specialized for two primes */

				if (entry->count < 2) {
					prime = r->large_prime[k ^ 1];
					entry = get_table_entry(table, 
								hashtable, 
								prime, 0);
					entry->count--;
					break;
				}
			}

            if (k == 2)
            {
                if (j != i)
                {
                    list[j] = list[i];
                }
                j++;
            }
		}
		num_relations = j;
		passes++;

	} while (num_left != num_relations);
			
	if (fobj->logfile != NULL)
		logprint(fobj->logfile, "reduce to %u relations in %u passes\n", 
				num_left, passes);
	if (fobj->VFLAG > 0)
		printf("reduce to %u relations in %u passes\n", 
				num_left, passes);
	return num_left;
}

uint32_t qs_purge_singletons3(fact_obj_t*fobj, siqs_r *list,
	uint32_t num_relations,
	qs_cycle_t *table, uint32_t *hashtable) {

	/* given a list of relations and the graph from the
	   sieving stage, remove any relation that contains
	   a prime that only occurs once in the graph. Because
	   removing a relation removes a second prime as well,
	   this process must be iterated until no more relations
	   are removed */

	uint32_t num_left;
	uint32_t i, j, k;
	uint32_t passes = 0;

	if (fobj->VFLAG > 0)
		printf("begin singleton removal with %u relations\n", 
			num_relations);
	if (fobj->logfile != NULL)
		logprint(fobj->logfile, "begin singleton removal with %u relations\n", num_relations);

	// start with all unique possible primes in table

	do {
		num_left = num_relations;

		/* for each relation */
		printf("now at pass %d: %u relations\n", passes, num_relations);

		for (i = j = 0; i < num_relations; i++) {
			siqs_r *r = list + i;
			uint32_t prime;
			qs_cycle_t *entry;

			/* full relations always survive */

			if ((r->large_prime[0] == 1) && (r->large_prime[1] == r->large_prime[2])) {
				if (j != i)
				{
					list[j] = list[i];
				}
				j++;
				continue;
			}

			/* for each prime in that relation */

			for (k = 0; k < 3; k++) {
				prime = r->large_prime[k];
				entry = get_table_entry(table, hashtable,
					prime, 0);

				/* if the relation is due to be removed,
				   decrement the count of its other
				   primes in the graph. */

				if (entry->count < 2) {
					int m;
					for (m = 0; m < 3; m++)
					{
						if (m == k) continue;
						prime = r->large_prime[m];
						entry = get_table_entry(table,
							hashtable,
							prime, 0);
						entry->count--;
					}
					break;
				}
			}

			if (k == 3)
			{
				// relation survived
				if (j != i)
				{
					list[j] = list[i];
				}
				j++;
			}
		}
		num_relations = j;
		passes++;

	} while (num_left != num_relations);

	if (fobj->logfile != NULL)
		logprint(fobj->logfile, "reduce to %u relations in %u passes\n",
			num_left, passes);
	if (fobj->VFLAG > 0)
		printf("reduce to %u relations in %u passes\n",
			num_left, passes);
	return num_left;
}


/*--------------------------------------------------------------------*/
void qs_enumerate_cycle(fact_obj_t*obj,
			    qs_la_col_t *c, 
			    qs_cycle_t *table,
			    qs_cycle_t *entry1, qs_cycle_t *entry2,
			    uint32_t final_relation) {

	/* given two entries out of the hashtable, corresponding
	   to two distinct primes, generate the list of relations
	   that participate in the cycle that these two primes
	   have just created. final_relation is the relation
	   to which the two primes belong, and the completed cycle
	   is packed into 'c' */

	uint32_t traceback1[100];
	uint32_t traceback2[100];
	uint32_t num1, num2;
	uint32_t i, j;

	/* Follow each cycle_t back up the graph until
	   the root component for this cycle is reached.
	   For each prime encountered along the way, save
	   the offset of the relation containing that prime */

	num1 = 0;
	while (entry1 != table + entry1->data) {
		
		//printf("entry1 step %u: entry %u points to prime %u, next step is %u\n",
		//	num1,entry1->data,table[entry1->data].prime,table[entry1->data].data);
		
		if (num1 >= 100) {
			if (obj->logfile != NULL)
				logprint(obj->logfile, "warning: cycle too long, "
					"skipping it\n");
			printf("warning: cycle too long, "
					"skipping it\n");
			return;
		}
		traceback1[num1++] = entry1->count;
		entry1 = table + entry1->data;
	}

	num2 = 0;
	while (entry2 != table + entry2->data) {
		//printf("entry2 step %u: entry %u points to prime %u, next step is %u\n",
		//	num2,entry2->data,table[entry2->data].prime,table[entry2->data].data);
		if (num2 >= 100) {
			if (obj->logfile != NULL)
				logprint(obj->logfile, "warning: cycle too long, "
					"skipping it\n");
			printf("warning: cycle too long, "
					"skipping it\n");
			return;
		}
		traceback2[num2++] = entry2->count;
		entry2 = table + entry2->data;
	}

	/* Now walk backwards through the lists, until
	   either one list runs out or a relation is
	   encountered that does not appear in both lists */

	while (num1 > 0 && num2 > 0) {
		if (traceback1[num1 - 1] != traceback2[num2 - 1])
			break;
		num1--; 
		num2--;
	}

	/* Now that we know how many relations are in the
	   cycle, allocate space to remember them */

	c->cycle.num_relations = num1 + num2 + 1;
	c->cycle.list = (uint32_t *)xmalloc(c->cycle.num_relations * 
					sizeof(uint32_t));
	
	/* Combine the two lists of relations */
	for (i = 0; i < num1; i++)
		c->cycle.list[i] = traceback1[i];

	for (j = 0; j < num2; j++, i++)
		c->cycle.list[i] = traceback2[j];

	/* Add the relation that created the cycle in the
	   first place */

	c->cycle.list[i] = final_relation;
}

void qs_enumerate_cycle3(fact_obj_t*obj,
	qs_la_col_t *c,
	qs_cycle_t *table,
	qs_cycle_t *entry1, qs_cycle_t *entry2, qs_cycle_t *entry3,
	uint32_t final_relation) {

	/* given two entries out of the hashtable, corresponding
	   to two distinct primes, generate the list of relations
	   that participate in the cycle that these two primes
	   have just created. final_relation is the relation
	   to which the two primes belong, and the completed cycle
	   is packed into 'c' */

	uint32_t traceback1[100];
	uint32_t traceback2[100];
	uint32_t num1, num2;
	uint32_t i, j;

	/* Follow each cycle_t back up the graph until
	   the root component for this cycle is reached.
	   For each prime encountered along the way, save
	   the offset of the relation containing that prime */

	num1 = 0;
	while (entry1 != table + entry1->data) {
		printf("entry1 step %u: entry %u points to prime %u, next step is %u\n",
			num1, entry1->data, table[entry1->data].prime, table[entry1->data].data);
		if (num1 >= 100) {
			if (obj->logfile != NULL)
				logprint(obj->logfile, "warning: cycle too long, "
					"skipping it\n");
			printf("warning: cycle too long, "
				"skipping it\n");
			return;
		}
		traceback1[num1++] = entry1->count;
		entry1 = table + entry1->data;
	}

	num2 = 0;
	while (entry2 != table + entry2->data) {
		printf("entry2 step %u: entry %u points to prime %u, next step is %u\n",
			num2, entry2->data, table[entry2->data].prime, table[entry2->data].data);
		if (num2 >= 100) {
			if (obj->logfile != NULL)
				logprint(obj->logfile, "warning: cycle too long, "
					"skipping it\n");
			printf("warning: cycle too long, "
				"skipping it\n");
			return;
		}
		traceback2[num2++] = entry2->count;
		entry2 = table + entry2->data;
	}

	/* Now walk backwards through the lists, until
	   either one list runs out or a relation is
	   encountered that does not appear in both lists */

	while (num1 > 0 && num2 > 0) {
		if (traceback1[num1 - 1] != traceback2[num2 - 1])
			break;
		num1--;
		num2--;
	}

	/* Now that we know how many relations are in the
	   cycle, allocate space to remember them */

	c->cycle.num_relations = num1 + num2 + 1;
	c->cycle.list = (uint32_t *)xmalloc(c->cycle.num_relations *
		sizeof(uint32_t));

	/* Combine the two lists of relations */
	for (i = 0; i < num1; i++)
		c->cycle.list[i] = traceback1[i];

	for (j = 0; j < num2; j++, i++)
		c->cycle.list[i] = traceback2[j];

	/* Add the relation that created the cycle in the
	   first place */

	c->cycle.list[i] = final_relation;
}

/*--------------------------------------------------------------------*/
int yafu_sort_cycles(const void *x, const void *y) {
	qs_la_col_t *xx = (qs_la_col_t *)x;
	qs_la_col_t *yy = (qs_la_col_t *)y;

	/* Callback for sorting a list of cycles by the
	   number of relations each contains */

	return xx->cycle.num_relations - yy->cycle.num_relations;
}
