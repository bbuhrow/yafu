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

// a slimmed down version of siqs, meant to run very quickly on tiny
// inputs (less than 115 bits, say).
void tinySIQS(fact_obj_t *fobj)
{
	//the input fobj->N and this 'n' are pointers to memory which holds
	//the input number to this function.  a copy in different memory
	//will also be made, and modified by a multiplier.  don't confuse
	//these two.
	//input expected in fobj->qs_obj->gmp_n	

	//master control structure
	static_conf_t *static_conf;
	dynamic_conf_t *dconf;

	//stuff for lanczos
	qs_la_col_t *cycle_list;
	uint32 num_cycles = 0;
	uint64 *bitfield = NULL;
	siqs_r *relation_list;

	//some locals
	uint32 num_needed, num_found;
	uint32 alldone;
	uint32 j;
	int i;
	clock_t start, stop;
	double t_time;
	struct timeval myTVend, optstart;
	TIME_DIFF *	difference;
	int updatecode = 0;

	int tmpTHREADS = THREADS;

	// checking savefile
	FILE *data;
	char tmpstr[GSTR_MAXSIZE];

	fobj->logfile = NULL;

	// we assume the input is not even, prime, a perfect power, or divisible
	// by small primes.  Just check the input size to see if we can
	// use squfof instead.
	if (mpz_sizeinbase(fobj->qs_obj.gmp_n, 2) > 115)
	{
		printf("input too big\n");
		return;
	}

	if (mpz_sizeinbase(fobj->qs_obj.gmp_n, 2) < 60)
	{
		mpz_t ztmp;
		mpz_init(ztmp);

		j = sp_shanks_loop(fobj->qs_obj.gmp_n, fobj);	
		if (j > 1)
		{
			mpz_set_64(ztmp, j);
			add_to_factor_list(fobj, ztmp);

			mpz_tdiv_q_ui(fobj->qs_obj.gmp_n, fobj->qs_obj.gmp_n, j);
			add_to_factor_list(fobj, fobj->qs_obj.gmp_n);

			mpz_set_ui(fobj->qs_obj.gmp_n, 1);			
		}
		
		return;
	}

	THREADS = 1;

	//At this point, we are committed to doing qs on the input
	//we need to:
	//1.) notify the screen and log
	//2.) start watching for an abort signal
	//3.) initialize data objects
	//4.) get ready to find the factor base

	//fill in the factorization object
	fobj->bits = mpz_sizeinbase(fobj->qs_obj.gmp_n, 2);
	fobj->digits = gmp_base10(fobj->qs_obj.gmp_n);
	fobj->qs_obj.savefile.name = (char *)malloc(80 * sizeof(char));
	strcpy(fobj->savefile_name,fobj->qs_obj.siqs_savefile);

	//initialize the data objects
	static_conf = (static_conf_t *)malloc(sizeof(static_conf_t));
	dconf = (dynamic_conf_t *)malloc(sizeof(dynamic_conf_t));
	static_conf->obj = fobj;

	//initialize offset for savefile buffer
	savefile_buf_off = 0;	
	
	//function pointer to the sieve array scanner
	scan_ptr = NULL;

	//start a counter for the whole job
	gettimeofday(&static_conf->totaltime_start, NULL);

	if (VFLAG >= 0)
		printf("\nstarting SIQS on c%d: %s\n",fobj->digits, 
		mpz_conv2str(&gstr1.s, 10, fobj->qs_obj.gmp_n));

	//get best parameters, multiplier, and factor base for the job
	//initialize and fill out the static part of the job data structure
	siqs_static_init(static_conf, 1);

	static_conf->in_mem = 1;
	static_conf->in_mem_relations = (siqs_r *)malloc(32768 * sizeof(siqs_r));
	static_conf->buffered_rel_alloc = 32768;
	static_conf->buffered_rels = 0;

	//allocate structures for use in sieving with threads
	siqs_dynamic_init(dconf, static_conf);

	//check if a savefile exists for this number, and if so load the data
	//into the master data structure
	alldone = siqs_check_restart(dconf, static_conf);
	print_siqs_splash(dconf, static_conf);

	//start the process
	num_needed = static_conf->factor_base->B + static_conf->num_extra_relations;
	num_found = static_conf->num_r;
	static_conf->total_poly_a = -1;

    while (1)
	{
		thread_sievedata_t t;

		// Check whether the thread has any results to collect. This should only be false at the 
		// very beginning, when the thread hasn't actually done anything yet.
		if (dconf->buffered_rels)
		{				
			num_found = siqs_merge_data(dconf,static_conf);

			// free sieving structure
			for (j=0; j<dconf->buffered_rels; j++)
				free(dconf->relation_buf[j].fb_offsets);
			dconf->num = 0;
			dconf->tot_poly = 0;
			dconf->buffered_rels = 0;

			//check whether to continue or not, and update the screen
			tiny_update_check(static_conf);
		}

		// if we have enough relations, or if there was a break signal, stop dispatching
		// any more threads
		if (num_found < num_needed) 
		{
			// generate a new poly A value for the thread we pulled out of the queue
			// using its dconf.  this is done by the master thread because it also 
			// stores the coefficients in a master list
			static_conf->total_poly_a++;
			new_poly_a(static_conf,dconf);
		}
		else
			break;


		//do some work
		t.dconf = dconf;
		t.sconf = static_conf;
		tiny_process_poly(&t);

	}	

	free_sieve(dconf);
	free(dconf->relation_buf);
	
	//finialize savefile
	qs_savefile_flush(&static_conf->obj->qs_obj.savefile);
	qs_savefile_close(&static_conf->obj->qs_obj.savefile);		
	
	tiny_update_final(static_conf);

	if (!static_conf->in_mem)
	{
		//we don't need the poly_a_list anymore... free it so the other routines
		//can use it
		for (i=0;i<static_conf->total_poly_a + 1;i++)
			mpz_clear(static_conf->poly_a_list[i]);
		free(static_conf->poly_a_list);
	}
	
	gettimeofday (&myTVend, NULL);
	difference = my_difftime (&static_conf->totaltime_start, &myTVend);

	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (VFLAG > 0)
	{
		printf("QS elapsed time = %6.4f seconds.\n",t_time);
		//printf("Predicted MAX_DIFF = %u, Actual MAX_DIFF = %u\n",MAX_DIFF,MAX_DIFF2);
		printf("\n==== post processing stage (msieve-1.38) ====\n");
	}

	fobj->qs_obj.qs_time = t_time;	
	
	start = clock();

	//return;

	//filter the relation set and get ready for linear algebra
	//all the polys and relations are on disk.
	//read them in and build the cycle list to send to LA.
	
	//initialize the b list for the current a.  qs_filter_relations
	//will change this as needed.
	static_conf->curr_b = (mpz_t *)malloc(2 * sizeof(mpz_t));
	for (i = 0; i < 2; i++)
		mpz_init(static_conf->curr_b[i]);
	static_conf->bpoly_alloc = 2;

	//load and filter relations and polys.
	yafu_qs_filter_relations(static_conf);

	cycle_list = static_conf->cycle_list;
	num_cycles = static_conf->num_cycles;
	relation_list = static_conf->relation_list;

	//solve the system of equations
	qs_solve_linear_system(static_conf->obj, static_conf->factor_base->B, 
		&bitfield, relation_list, cycle_list, &num_cycles);

	stop = clock();
	static_conf->t_time1 = (double)(stop - start)/(double)CLOCKS_PER_SEC;

	start = clock();
	//sqrt stage
	if (bitfield != NULL && num_cycles > 0) 
	{
	
		yafu_find_factors(static_conf->obj, static_conf->n, static_conf->factor_base->list, 
			static_conf->factor_base->B, cycle_list, num_cycles, 
			relation_list, bitfield, static_conf->multiplier, 
			static_conf->poly_a_list, static_conf->poly_list, 
			&static_conf->factor_list);
					
	}

	stop = clock();
	static_conf->t_time2= (double)(stop - start)/(double)CLOCKS_PER_SEC;

	gettimeofday (&myTVend, NULL);
	difference = my_difftime (&static_conf->totaltime_start, &myTVend);

	static_conf->t_time3 = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);
	
	if (VFLAG > 0)
	{
		printf("Lanczos elapsed time = %6.4f seconds.\n",static_conf->t_time1);
		printf("Sqrt elapsed time = %6.4f seconds.\n",static_conf->t_time2);
		
	}

	if (VFLAG >= 0)
		printf("SIQS elapsed time = %6.4f seconds.\n",static_conf->t_time3);

	fobj->qs_obj.total_time = static_conf->t_time3;

	static_conf->cycle_list = cycle_list;
	static_conf->num_cycles = num_cycles;

	//free stuff used during filtering
	free_filter_vars(static_conf);

done:

	//free everything else
	free_siqs(static_conf);
	free(dconf);
	free(static_conf);

	tmpTHREADS = THREADS;

	return;
}


void *tiny_process_poly(void *ptr)
//void process_hypercube(static_conf_t *sconf,dynamic_conf_t *dconf)
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
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p);
			med_sieve_ptr(sieve, fb_sieve_p, fb, start_prime, blockinit);
			lp_sieveblock(sieve, i, num_blocks, buckets, 0);

			//set the roots for the factors of a to force the following routine
			//to explicitly trial divide since we haven't found roots for them
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p);
			scan_ptr(i,0,sconf,dconf);

			//set the roots for the factors of a such that
			//they will not be sieved.  we haven't found roots for them
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n);
			med_sieve_ptr(sieve, fb_sieve_n, fb, start_prime, blockinit);
			lp_sieveblock(sieve, i, num_blocks, buckets, 1);

			//set the roots for the factors of a to force the following routine
			//to explicitly trial divide since we haven't found roots for them
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n);
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

int tiny_update_check(static_conf_t *sconf)
{
	//check to see if we should combine partial relations
	//found and check to see if we're done.  also update the screen.
	//this happens after every polynomial is sieved, since tiny numbers
	//go so fast.

	uint32 num_full = sconf->num_relations;
	uint32 check_total = sconf->check_total;
	uint32 check_inc = sconf->check_inc;
	int i;
	fb_list *fb = sconf->factor_base;
	int retcode = 0;

	//update status on screen
	sconf->num_r = sconf->num_relations + 
	sconf->num_cycles +
	sconf->components - sconf->vertices;

	//also change rel sum to update_rels below...
	if (VFLAG >= 0)
	{
		uint32 update_rels;

		//in order to keep rels/sec from going mad when relations
		//are reloaded on a restart, just use the number of
		//relations we've found since the last update.  don't forget
		//to initialize last_numfull and partial when loading
		update_rels = sconf->num_relations + sconf->num_cycles - 
			sconf->last_numfull - sconf->last_numcycles;
		sconf->charcount = printf("%d rels found: %d full + "
			"%d from %d partial\n",
			sconf->num_r,sconf->num_relations,
			sconf->num_cycles +
			sconf->components - sconf->vertices,
			sconf->num_cycles);

		fflush(stdout);
	}
		
	if (sconf->num_r >= fb->B + sconf->num_extra_relations) 
	{
		//we've got enough total relations to stop
		return 1;
	}	

	return 0;
}

int tiny_update_final(static_conf_t *sconf)
{
	mpz_t tmp1;
	struct timeval myTVend;
	TIME_DIFF *	difference;

	mpz_init(tmp1);

	if (VFLAG >= 0)
	{
		mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
		mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
		mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
		mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block	

		if (VFLAG > 0)
		{
			printf("\n\nsieving required %d total polynomials\n",
				sconf->tot_poly);
			gmp_printf("trial division touched %d sieve locations out of %Zd\n",
				sconf->num, tmp1);
		}
		else
			printf("\n\n");

		fflush(stdout);
		fflush(stderr);
	}
	
	gettimeofday (&myTVend, NULL);
	difference = my_difftime (&sconf->totaltime_start, &myTVend);

	sconf->obj->qs_obj.rels_per_sec = (double)(sconf->num_relations + sconf->num_cycles) /
		((double)difference->secs + (double)difference->usecs / 1000000);

	mpz_clear(tmp1);
	free(difference);

	return 0;
}


