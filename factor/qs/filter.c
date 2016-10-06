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


/*
this file contains code to implement filtering of relation sets prior
to matrix construction, and restarting of siqs from a data file
*/

uint32 process_poly_a(static_conf_t *sconf)
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
	for (j = 0; (uint32)j < sconf->bpoly_alloc; j++)
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
	fp_digit r;

	mpz_set(tmp, poly->mpz_poly_a);
	k=0;
	poly->s = 0;
	while ((mpz_cmp_ui(tmp, 1) != 0) && (spSOEprimes[k] < 65536))
	{
		r = mpz_tdiv_ui(tmp,(fp_digit)spSOEprimes[k]);
		
		if (r != 0)
			k++;
		else
		{
			for (j=1;(uint32)j < fb->B;j++)
				if (fb->list->prime[j] == spSOEprimes[k])
					break;
			if ((uint32)j >= fb->B)
			{
				//then we didn't find this factor in the fb, thus the 
				//read in A is probably bad.
				printf("bad poly A encountered in savefile\n");
				return 1;
			}
			poly->qlisort[poly->s] = j;
			poly->s++;
			mpz_tdiv_q_ui(tmp, tmp, (fp_digit)spSOEprimes[k]);
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



int process_rel(char *substr, fb_list *fb, mpz_t n,
				 static_conf_t *sconf, fact_obj_t *obj, siqs_r *rel)
{
	char *nextstr;
	uint32 lp[2];
	uint32 this_offset, this_id, this_num_factors, this_parity, this_val;
	uint32 fb_offsets[MAX_SMOOTH_PRIMES];
	int i,j,k, err_code = 0;

#ifdef SPARSE_STORE
    mpz_t Q;
    mpz_init(Q);
#endif

	//read the offset
	this_offset = strtoul(substr,&nextstr,HEX);	//convert
	substr = nextstr;

	this_id = strtoul(substr,&nextstr,HEX);
	substr = nextstr;

	if (this_offset & 0x80000000)
	{
		this_parity = 1;
		this_offset = (~this_offset) + 1;
	}
	else
		this_parity = 0;

	j=0;
	do
	{
		this_val = strtoul(substr,&nextstr,HEX);
		if (this_val == (uint32)(-1))
			printf("error parsing relation: strtoul returned error code\n");
		substr = nextstr;
		fb_offsets[j] = this_val;
		j++;
	} while (substr[1] != 'L');
	this_num_factors = j;

	substr += 2;
	this_val = strtoul(substr,&nextstr,HEX);
	substr = nextstr;
	lp[0] = this_val;

	this_val = strtoul(substr,&nextstr,HEX);
	substr = nextstr;
    lp[1] = this_val;

    // combine the factors of the sieve value with
    // the factors of the polynomial 'a' value; the 
    // linear algebra code has to know about both.
    // Because both lists are sorted, this is just
    // a merge operation

#ifdef SPARSE_STORE
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
        uint32 prime = sconf->factor_base->list->prime[i];
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
        uint32 prime;

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
        uint32 prime = sconf->factor_base->list->prime[sconf->curr_poly->qlisort[j]];

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

#else

    // merge in the list of polya coefficients to the list of
    // relation factors that we read in.  both lists are sorted,
    // so this mergesort works quickly.
	i = j = k = 0;
	while (i < (int)this_num_factors && j < sconf->curr_poly->s) {
		if (fb_offsets[i] < sconf->curr_poly->qlisort[j]) {
			rel->fb_offsets[k++] = fb_offsets[i++];
		}
		else if (fb_offsets[i] > sconf->curr_poly->qlisort[j]) {
			rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j++];
		}
		else {
			rel->fb_offsets[k] = fb_offsets[i++];
			rel->fb_offsets[k+1] = sconf->curr_poly->qlisort[j++];
			k += 2;
		}
	}
	while (i < (int)this_num_factors)
		rel->fb_offsets[k++] = fb_offsets[i++];
	while (j < sconf->curr_poly->s)
		rel->fb_offsets[k++] = sconf->curr_poly->qlisort[j++];

    this_num_factors += sconf->curr_poly->s

#endif
	
	rel->sieve_offset = this_offset;
	rel->large_prime[0] = lp[0];
	rel->large_prime[1] = lp[1];
	rel->parity = this_parity;
    rel->num_factors = this_num_factors; 
	rel->poly_idx = this_id;

	if (!check_relation(sconf->curr_a,
			sconf->curr_b[rel->poly_idx], rel, fb, n))
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


	return err_code;	//error code, if there is one.
}

int restart_siqs(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	int i,j;
	char *str, *substr;
	FILE *data;
	uint32 lp[2],pmax = sconf->large_prime_max / sconf->large_mult;
	//fact_obj_t *obj = sconf->obj;

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
			if (VFLAG > 1)
				printf("restarting siqs from saved data set\n");
			fflush(stdout);
			fflush(stderr);
			while (1)
			{
				//read a line
				if (feof(data))
					break;
				fgets(str,1024,data);
				substr = str + 2;

				if (str[0] == 'R')
				{	
					//process a relation
					//just trying to figure out how many relations we have
					//so read in the large primes and add to cycles
					substr = strchr(substr,'L');
					yafu_read_large_primes(substr,lp,lp+1);
					if (sconf->use_dlp)
					{
						if ((lp[0] > 1) && (lp[0] < pmax))
						{
							j++;
							continue;
						}
						if ((lp[1] > 1) && (lp[1] < pmax))
						{
							j++;
							continue;
						}
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

			if (VFLAG > 0)
			{
				printf("%d relations found: %d full + "
					"%d from %d partial\n",
					sconf->num_r,sconf->num_relations,
					sconf->num_cycles +
					sconf->components - sconf->vertices,
					sconf->num_cycles);
				printf("threw away %d relations with large primes too small\n",j);
				fflush(stdout);
				sconf->last_numfull = sconf->num_relations;
				sconf->last_numcycles = sconf->num_cycles;
			}

		}	
		fclose(data);
	}
	free(str);

	return 0;
}

#define QS_HASH_MULT ((uint32)(2654435761UL))
#define QS_HASH(a) (((a) * QS_HASH_MULT) >> (32 - QS_LOG2_CYCLE_HASH))

/**********************************************************
These 3 routines are used to add a relation to the cycle
tree every time one is found after sieving
(e.g., add_to_cycles is called in save_relation)
**********************************************************/
static qs_cycle_t *get_table_entry(qs_cycle_t *table, uint32 *hashtable,
				uint32 prime, uint32 new_entry_offset) {

	/* return a pointer to a unique qs_cycle_t specific
	   to 'prime'. The value of 'prime' is hashed and
	   the result used to index 'table'. If prime does
	   not appear in table, specify 'new_entry_offset'
	   as the table entry for prime.

	   Values of 'prime' which hash to the same offset
	   in hashtable are connected by a linked list of
	   offsets in 'table'. */

	uint32 offset, first_offset;
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

/*--------------------------------------------------------------------*/
static uint32 add_to_hashtable(qs_cycle_t *table, uint32 *hashtable, 
			uint32 prime1, uint32 prime2, 
			uint32 default_table_entry, 
			uint32 *components, uint32 *vertices) {

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

	uint32 root[2];
	uint32 root1, root2;
	uint32 i;
	uint32 num_new_entries = 0;

	/* for each prime */
	

	for (i = 0; i < 2; i++) {
		uint32 prime = ((i == 0) ? prime1 : prime2);
		uint32 offset; 
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
void yafu_add_to_cycles(static_conf_t *conf, uint32 flags, uint32 prime1, uint32 prime2) {

	/* Top level routine for updating the graph of partial
	   relations */

	uint32 table_size = conf->cycle_table_size;
	uint32 table_alloc = conf->cycle_table_alloc;
	qs_cycle_t *table = conf->cycle_table;
	uint32 *hashtable = conf->cycle_hashtable;

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

	fact_obj_t *obj = sconf->obj;
	uint32 *hashtable = sconf->cycle_hashtable;
	qs_cycle_t *table = sconf->cycle_table;
	uint32 num_derived_poly;
	uint32 *final_poly_index;
	uint32 num_relations, num_cycles, num_poly;
	qs_la_col_t *cycle_list;
	siqs_r *relation_list;

	uint32 i, passes, start;
	uint32 curr_a_idx, curr_poly_idx, curr_rel; 
	uint32 curr_expected, curr_saved, curr_cycle; 
	uint32 total_poly_a;
	uint32 poly_saved;
	uint32 cycle_bins[NUM_CYCLE_BINS+1] = {0};
	char buf[LINE_BUF_SIZE];
	char *subbuf;
	uint32 last_id;
	int first, last_poly;
	uint32 this_rel = 0;

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
		qs_savefile_open(&obj->qs_obj.savefile, SAVEFILE_READ);
		qs_savefile_read_line(buf, sizeof(buf), &obj->qs_obj.savefile);

		//we don't know beforehand how many rels to expect, so start
		//with some amount and allow it to increase as we read them
		relation_list = (siqs_r *)xmalloc(10000 * sizeof(siqs_r));
		curr_rel = 10000;
		while (!qs_savefile_eof(&obj->qs_obj.savefile)) {
			char *start;

			switch (buf[0]) {
			case 'A':
				total_poly_a++;
				break;

			case 'R':
				start = strchr(buf, 'L');
				if (start != NULL) {
					uint32 prime1, prime2;
					yafu_read_large_primes(start, &prime1, &prime2);
					if (i == curr_rel) {
						curr_rel = 3 * curr_rel / 2;
						relation_list = (siqs_r *)xrealloc(
								relation_list,
								curr_rel *
								sizeof(siqs_r));
					}
					relation_list[i].poly_idx = i;
					relation_list[i].large_prime[0] = prime1;
					relation_list[i].large_prime[1] = prime2;
					i++;
				}
				break;
			}

			qs_savefile_read_line(buf, sizeof(buf), &obj->qs_obj.savefile);
		}
		num_relations = i;
	}
	else
	{
		relation_list = (siqs_r *)xmalloc(sconf->buffered_rels * sizeof(siqs_r));
		for (i=0; i<sconf->buffered_rels; i++)
		{
			relation_list[i].poly_idx = i;
			relation_list[i].large_prime[0] = sconf->in_mem_relations[i].large_prime[0];
			relation_list[i].large_prime[1] = sconf->in_mem_relations[i].large_prime[1];
		}
		total_poly_a = sconf->total_poly_a;
		num_relations = sconf->buffered_rels;
	}
		
	num_relations = qs_purge_singletons(obj, relation_list, num_relations,
					table, hashtable);

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

	final_poly_index = (uint32 *)xmalloc(1024 * 
						sizeof(uint32));
	
	/* initialize the running counts of relations and
	   polynomials */

	i = 0;
	last_id = 999;
	last_poly = -1;
	curr_expected = 0;
	curr_saved = 0;
	curr_rel = (uint32)(-1);
	curr_poly_idx = (uint32)(-1);
	curr_a_idx = (uint32)(-1);
	poly_saved = 0;
	sconf->poly_list_alloc = 0;
	if (VFLAG > 0)
		printf("attempting to read %u relations\n", num_relations);

	/* Read in the relations and the polynomials they use
	   at the same time. */

	if (!sconf->in_mem)
		qs_savefile_rewind(&obj->qs_obj.savefile);
	else
		this_rel = 0;

	first = 1;
	while (curr_expected < num_relations) {
		char *tmp;
		uint32 bad_A_val = 0;
		siqs_r *r;
		siqs_r *rel;

		/* read in the next entity */
		if (!sconf->in_mem)
		{
			if (qs_savefile_eof(&obj->qs_obj.savefile))
				break;
			qs_savefile_read_line(buf, sizeof(buf), &obj->qs_obj.savefile);
		}
		else
		{
			//printf("this rel: %u, curr_rel: %u, curr_expected: %u, "
			//	"curr_saved: %u, num_relations: %u, buffered_rels: %u, i: %u\n",
			//	this_rel, curr_rel, curr_expected, curr_saved, num_relations, 
			//	sconf->buffered_rels, i);

			if (this_rel == sconf->buffered_rels)
				break;

			rel = sconf->in_mem_relations + this_rel++;
			if (rel->poly_idx < last_id)
				buf[0] = 'A';
			else
				buf[0] = 'R';

			last_id = rel->poly_idx;
		}

		switch (buf[0]) {
		case 'A':
			/* Read in a new 'a' value */
			/* build all of the 'b' values associated with it */
			if (!sconf->in_mem)
			{
				subbuf = buf + 2;	//skip the A and a space
				mpz_set_str(sconf->curr_a, subbuf, 0);
			}
			else
			{
				last_poly++;
				this_rel--;
				mpz_set(sconf->curr_a, sconf->poly_a_list[last_poly]);
			}

			curr_a_idx++;			
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

			mpz_set(sconf->poly_a_list[curr_a_idx], sconf->curr_a); 

			/* all 'b' values start off unused */
			final_poly_index = (uint32 *)xrealloc(final_poly_index,
				num_derived_poly * sizeof(uint32));
			memset(final_poly_index, -1, num_derived_poly *
							sizeof(uint32));

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
				if (process_rel(subbuf, sconf->factor_base,
					sconf->n, sconf, sconf->obj, r)) {

						if (obj->logfile != NULL)
							logprint(obj->logfile, "failed to read relation %d\n", 
								curr_expected - 1);
						if (VFLAG > 0)
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
				//combine the factors of the sieve value with
				//  the factors of the polynomial 'a' value; the 
				//  linear algebra code has to know about both.
				//  Because both lists are sorted, this is just
				//  a merge operation 
				int ii,jj,kk;				
				r = relation_list + curr_saved;

				//r->fb_offsets = (uint32 *)xmalloc(
				//	(rel->num_factors + sconf->curr_poly->s) * sizeof(uint32));

				ii = jj = kk = 0;
				while (ii < (int)rel->num_factors && jj < sconf->curr_poly->s) {
					if (rel->fb_offsets[ii] < sconf->curr_poly->qlisort[jj]) {
						r->fb_offsets[kk++] = rel->fb_offsets[ii++];
					}
					else if (rel->fb_offsets[ii] > sconf->curr_poly->qlisort[jj]) {
						r->fb_offsets[kk++] = sconf->curr_poly->qlisort[jj++];
					}
					else {
						r->fb_offsets[kk] = rel->fb_offsets[ii++];
						r->fb_offsets[kk+1] = sconf->curr_poly->qlisort[jj++];
						kk += 2;
					}
				}
				while (ii < (int)rel->num_factors)
					r->fb_offsets[kk++] = rel->fb_offsets[ii++];
				while (jj < sconf->curr_poly->s)
					r->fb_offsets[kk++] = sconf->curr_poly->qlisort[jj++];
	
				r->num_factors = rel->num_factors + sconf->curr_poly->s;
				r->large_prime[0] = rel->large_prime[0];
				r->large_prime[1] = rel->large_prime[1];
				r->parity = rel->parity;
				r->sieve_offset = rel->sieve_offset;
				r->poly_idx = rel->poly_idx;

				if (!check_relation(sconf->curr_a,
						sconf->curr_b[r->poly_idx], r, sconf->factor_base, sconf->n))
				{
					if (r->large_prime[0] != r->large_prime[1]) {
						yafu_add_to_cycles(sconf, obj->flags, 
							r->large_prime[0], r->large_prime[1]);
						sconf->num_cycles++;
					}
					else {
						sconf->num_relations++;
					}
				}

			}

			curr_saved++;

			/* if necessary, save the b value corresponding 
			   to this relation */

			if (final_poly_index[r->poly_idx] == (uint32)(-1)) {

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
	if (obj->logfile != NULL)
	{
		logprint(obj->logfile, "recovered %u relations\n", num_relations);
		logprint(obj->logfile, "recovered %u polynomials\n", i);
	}
	if (VFLAG > 0)
	{
		printf("recovered %u relations\n", num_relations);
		printf("recovered %u polynomials\n", i);
	}

	if (!sconf->in_mem)
		qs_savefile_close(&obj->qs_obj.savefile);

	free(final_poly_index);
	sconf->poly_list = (poly_t *)xrealloc(sconf->poly_list,
					   i * sizeof(poly_t));

	/* begin the cycle generation process by purging
	   duplicate relations. For the sake of consistency, 
	   always rebuild the graph afterwards */

	num_relations = qs_purge_duplicate_relations(obj, 
				relation_list, num_relations);

	memset(hashtable, 0, sizeof(uint32) * (1 << QS_LOG2_CYCLE_HASH));
	sconf->vertices = 0;
	sconf->components = 0;
	sconf->cycle_table_size = 1;

	for (i = 0; i < num_relations; i++) {
		siqs_r *r = relation_list + i;
		if (r->large_prime[0] != r->large_prime[1]) {		
			yafu_add_to_cycles(sconf, sconf->obj->flags, r->large_prime[0], 
					r->large_prime[1]);
		}
	}
	
	/* compute the number of cycles to expect. Note that
	   this number includes cycles from both full and partial
	   relations (the cycle for a full relation is trivial) */

	num_cycles = num_relations + sconf->components - sconf->vertices;

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
		uint32 offset = hashtable[i];		

		while (offset != 0) {
			qs_cycle_t *entry = table + offset;

			if (offset != entry->data)
				entry->data = 0;
			offset = entry->next;
		}
	}

	if (obj->logfile != NULL)
		logprint(obj->logfile, "attempting to build %u cycles\n", num_cycles);
	if (VFLAG > 0)
	{
		printf("attempting to build %u cycles\n", num_cycles);
		fflush(stdout);
	}
	cycle_list = (qs_la_col_t *)xmalloc(num_cycles * sizeof(qs_la_col_t));

	/* keep going until either all cycles are found, all
	   relations are processed, or cycles stop arriving. 
	   Normally these conditions all occur at the same time */

	for (start = passes = curr_cycle = 0; start < num_relations && 
			curr_cycle < num_cycles; passes++) {

		/* The list of relations up to index 'start' is con-
		   sidered processed. For all relations past that... */

		uint32 start_cycles = curr_cycle;

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
				c->cycle.list = (uint32 *)
						xmalloc(sizeof(uint32));
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
				qs_la_col_t *c = cycle_list + curr_cycle;
				c->cycle.list = NULL;
				qs_enumerate_cycle(obj, c, table, entry1,
						entry2, start);
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
	num_cycles = curr_cycle;

	if (obj->logfile != NULL)
		logprint(obj->logfile, "found %u cycles in %u passes\n", num_cycles, passes);
	if (VFLAG > 0)
		printf("found %u cycles in %u passes\n", num_cycles, passes);
	
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

	if (obj->logfile != NULL)
		logprint(obj->logfile, "distribution of cycle lengths:\n");

	if (VFLAG > 0)
		printf("distribution of cycle lengths:\n");
	for (i = 0; i < NUM_CYCLE_BINS; i++) {
		if (cycle_bins[i]) {
			if (obj->logfile != NULL)
				logprint(obj->logfile, "   length %d : %d\n", 
					i + 1, cycle_bins[i]);
			if (VFLAG > 0)
				printf("   length %d : %d\n", 
					i + 1, cycle_bins[i]);
		}
	}
	if (cycle_bins[i])
	{
		if (obj->logfile != NULL)
			logprint(obj->logfile, "   length %u+: %u\n", i + 1, cycle_bins[i]);

		if (VFLAG > 0)
			printf("   length %u+: %u\n", i + 1, cycle_bins[i]);
	}
	
	if (obj->logfile != NULL)
		logprint(obj->logfile, "largest cycle: %u relations\n",
			cycle_list[num_cycles-1].cycle.num_relations);

	if (VFLAG > 0)
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
	uint32 i;

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
uint32 qs_purge_duplicate_relations(fact_obj_t *obj,
				siqs_r *rlist, 
				uint32 num_relations) {

	uint32 i, j;
	
	/* remove duplicates from rlist */

	if (num_relations < 2)
		return num_relations;

	qsort(rlist, (size_t)num_relations,
		sizeof(siqs_r), compare_relations);

	for (i = 1, j = 0; i < num_relations; i++) {
		if (compare_relations(rlist + j, rlist + i) == 0)
		{

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
		if (VFLAG > 0)
			printf("freed %d duplicate relations\n", 
					num_relations - j);
	}

	return j;
}

void yafu_read_large_primes(char *buf, uint32 *prime1, uint32 *prime2) {

	char *next_field;
	uint32 p1, p2;

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

uint32 qs_purge_singletons(fact_obj_t *obj, siqs_r *list, 
				uint32 num_relations,
				qs_cycle_t *table, uint32 *hashtable) {
	
	/* given a list of relations and the graph from the
	   sieving stage, remove any relation that contains
	   a prime that only occurs once in the graph. Because
	   removing a relation removes a second prime as well,
	   this process must be iterated until no more relations
	   are removed */

	uint32 num_left;
	uint32 i, j, k;
	uint32 passes = 0;

	if (VFLAG > 0)
		printf("begin with %u relations\n", num_relations);
	if (obj->logfile != NULL)
		logprint(obj->logfile, "begin with %u relations\n", num_relations);

	do {
		num_left = num_relations;

		/* for each relation */

		for (i = j = 0; i < num_relations; i++) {
			siqs_r *r = list + i;
			uint32 prime;
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
			
	if (obj->logfile != NULL)
		logprint(obj->logfile, "reduce to %u relations in %u passes\n", 
				num_left, passes);
	if (VFLAG > 0)
		printf("reduce to %u relations in %u passes\n", 
				num_left, passes);
	return num_left;
}

/*--------------------------------------------------------------------*/
void qs_enumerate_cycle(fact_obj_t *obj, 
			    qs_la_col_t *c, 
			    qs_cycle_t *table,
			    qs_cycle_t *entry1, qs_cycle_t *entry2,
			    uint32 final_relation) {

	/* given two entries out of the hashtable, corresponding
	   to two distinct primes, generate the list of relations
	   that participate in the cycle that these two primes
	   have just created. final_relation is the relation
	   to which the two primes belong, and the completed cycle
	   is packed into 'c' */

	uint32 traceback1[100];
	uint32 traceback2[100];
	uint32 num1, num2;
	uint32 i, j;

	/* Follow each cycle_t back up the graph until
	   the root component for this cycle is reached.
	   For each prime encountered along the way, save
	   the offset of the relation containing that prime */

	num1 = 0;
	while (entry1 != table + entry1->data) {
		//printf("entry1 step %u: entry %u points to prime %u, next step is %u\n",
			//num1,entry1->data,table[entry1->data].prime,table[entry1->data].data);
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
			//num2,entry2->data,table[entry2->data].prime,table[entry2->data].data);
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
	c->cycle.list = (uint32 *)xmalloc(c->cycle.num_relations * 
					sizeof(uint32));
	
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
