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
#include "threadpool.h"
#include "cofactorize.h"

//#define VDEBUG
// opt debug will print out some info relevant to the process used to 
// optimize the tf_small_cutoff value.  This value is varied slightly during the
// first few polynomials and the value which maximizes relation discover rate
// is chosen.
//#define OPT_DEBUG

typedef struct
{
    thread_sievedata_t *thread_data;

    // adaptive tf_small_cutoff variables needed across syncronizations
    int orig_value;
    int num_meas;
    double results[3];
    int averaged_polys;
    double rels_per_sec_avg;
#ifdef OPT_DEBUG
    FILE *optfile;
#endif

    // other stuff that we need a global copy of
    uint32 num_needed;
    uint32 num_found;
    int updatecode;
    uint32 avg_rels_per_acoeff;
} siqs_userdata_t;

void siqs_start(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    siqs_userdata_t *udata = tdata->user_data;
    static_conf_t *static_conf = udata->thread_data[0].sconf;
    fact_obj_t *fobj = static_conf->obj;
    int i;

#ifdef OPT_DEBUG
    udata->optfile = fopen("optfile.csv", "a");
    fprintf(udata->optfile, "Optimization Debug File\n");
    fprintf(udata->optfile, "Detected cpu %d, with L1 = %d bytes, L2 = %d bytes\n",
        yafu_get_cpu_type(), L1CACHE, L2CACHE);
    gmp_fprintf(udata->optfile, "Starting SIQS on c%d: %s\n\n",
        fobj->digits, mpz_conv2str(&gstr1.s, 10, fobj->qs_obj.gmp_n));
    fprintf(udata->optfile, "Meas #,Poly A #, Avg Rels/Poly/Sec, small_tf_cutoff\n");
#endif

    udata->results[0] = 0.0;
    udata->results[1] = 0.0;
    udata->results[2] = 0.0;

    udata->averaged_polys = 0;
    udata->rels_per_sec_avg = 0.0;
    udata->updatecode = 0;

    // get best parameters, multiplier, and factor base for the job
    // initialize and fill out the static part of the job data structure
    siqs_static_init(udata->thread_data[0].sconf, 0);

    // allocate structures for use in sieving with threads
    for (i = 0; i < THREADS; i++)
    {
        siqs_dynamic_init(udata->thread_data[i].dconf, static_conf);
    }

    // check if a savefile exists for this number, and if so load the data
    // into the master data structure
    siqs_check_restart(udata->thread_data[0].dconf, static_conf);
    print_siqs_splash(udata->thread_data[0].dconf, static_conf);

    // start the process
    udata->num_needed = static_conf->factor_base->B + static_conf->num_extra_relations;
    static_conf->num_needed = udata->num_needed;
    udata->num_found = static_conf->num_r;
    static_conf->total_poly_a = -1;
    udata->num_meas = 0;
    udata->orig_value = static_conf->tf_small_cutoff;

    if (fobj->qs_obj.gbl_override_small_cutoff_flag)
    {
        printf("overriding small TF cutoff at %d\n", fobj->qs_obj.gbl_override_small_cutoff);
        fobj->qs_obj.no_small_cutoff_opt = 1;
        static_conf->tf_small_cutoff = fobj->qs_obj.gbl_override_small_cutoff;
    }

    return;
}

void siqs_sync(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    siqs_userdata_t *udata = tdata->user_data;
    thread_sievedata_t *t = udata->thread_data;
    static_conf_t *static_conf = udata->thread_data[0].sconf;
    fact_obj_t *fobj = static_conf->obj;
    
    int tid = tdata->tindex;
    
    // adaptive tf_small_cutoff variables    
    int num_avg = 10;    
    int orig_value = udata->orig_value;
    int j;
    
    // Check whether the thread has any results to collect and merge them if so.
    if (t[tid].dconf->buffered_rels)
    {
        udata->num_found = siqs_merge_data(t[tid].dconf, static_conf);        
        static_conf->num_found = udata->num_found;

        if (fobj->qs_obj.no_small_cutoff_opt == 0)
        {
            int	poly_start_num = 0;

            if (udata->num_meas < 3)
            {
                if (udata->averaged_polys >= num_avg)
                {
                    poly_start_num = static_conf->total_poly_a;
                    udata->results[udata->num_meas] = udata->rels_per_sec_avg /
                        (double)udata->averaged_polys;

#ifdef OPT_DEBUG
                    fprintf(udata->optfile, "%d,%d,%f,%d\n", udata->num_meas, 
                        static_conf->total_poly_a,
                        udata->results[udata->num_meas], 
                        static_conf->tf_small_cutoff);
#endif

                    udata->rels_per_sec_avg = 0.0;
                    udata->averaged_polys = 0;

                    if (udata->num_meas == 0)
                    {
                        static_conf->tf_small_cutoff = orig_value + 5;
                    }
                    else if (udata->num_meas == 1)
                    {
                        static_conf->tf_small_cutoff = orig_value - 5;
                    }
                    else
                    {
                        //we've got our three measurements, make a decision.
                        //the experimental results need to be convincingly better
                        //in order to switch to a different value.  2% for numbers
                        //with one LP, 5% for DLP.
                        if (static_conf->use_dlp)
                            udata->results[0] *= 1.05;
                        else
                            udata->results[0] *= 1.02;

                        if (udata->results[0] > udata->results[1])
                        {
                            if (udata->results[0] > udata->results[2])
                                static_conf->tf_small_cutoff = orig_value;
                            else
                                static_conf->tf_small_cutoff = orig_value - 5;
                        }
                        else
                        {
                            if (udata->results[1] > udata->results[2])
                                static_conf->tf_small_cutoff = orig_value + 5;
                            else
                                static_conf->tf_small_cutoff = orig_value - 5;
                        }

                        //printf("\nfinal value for tf_small_cutoff = %u\n",
                        //    static_conf->tf_small_cutoff);
#ifdef OPT_DEBUG
                        fprintf(udata->optfile, "final value = %d\n", static_conf->tf_small_cutoff);
                        fprintf(udata->optfile, "\n\n");
                        fclose(udata->optfile);
#endif
                    }

                    udata->num_meas++;
                }
                else
                {
                    udata->rels_per_sec_avg += t[tid].dconf->rels_per_sec;
                    udata->averaged_polys++;
                }
            }
        }

        // free sieving structure
        for (j = 0; j < t[tid].dconf->buffered_rels; j++)
            t[tid].dconf->relation_buf[j].num_factors = 0;

        t[tid].dconf->num = 0;
        t[tid].dconf->tot_poly = 0;
        t[tid].dconf->buffered_rels = 0;
        t[tid].dconf->attempted_squfof = 0;
        t[tid].dconf->failed_squfof = 0;
        t[tid].dconf->dlp_outside_range = 0;
        t[tid].dconf->dlp_prp = 0;
        t[tid].dconf->dlp_useful = 0;
        t[tid].dconf->total_blocks = 0;
        t[tid].dconf->total_reports = 0;
        t[tid].dconf->total_surviving_reports = 0;
        t[tid].dconf->lp_scan_failures = 0;
        t[tid].dconf->num_64bit_residue = 0;
    }

    return;
}

void siqs_dispatch(void *vptr)
{
    tpool_t *tdata = (tpool_t *)vptr;
    siqs_userdata_t *udata = tdata->user_data;
    thread_sievedata_t *t = udata->thread_data;
    static_conf_t *static_conf = t[0].sconf;
    fact_obj_t *fobj = static_conf->obj;
    int tid = tdata->tindex;    
    double avg_rels_per_acoeff;
    double in_flight_rels;

    //check whether to continue or not, and update the screen
    udata->updatecode = update_check(static_conf);
    if ((udata->num_found > 0) && ((int)static_conf->total_poly_a > 0))
    {
        avg_rels_per_acoeff = (double)udata->num_found / (double)(static_conf->total_poly_a + 1);

        if (static_conf->is_restart == 1)
        {
            // we'd have to separately track the relations found since the
            // restart for this to work right.  that impacts many things,
            // so for now just don't worry about overshoot due to multi-threading.
            in_flight_rels = 0;
        }
        else
        {
            in_flight_rels = (tdata->num_threads - 1) * 0.9 * avg_rels_per_acoeff;
        }
        //printf("found %u rels so far in %d 'A' polys, avg rels_per_poly = %f, rels in flight = %f\n",
        //    udata->num_found, static_conf->total_poly_a + 1, avg_rels_per_acoeff, in_flight_rels);
    }
    else
    {
        in_flight_rels = 0;
    }

    // if we have enough relations, or if there was a break signal, stop dispatching
    // any more threads
    if ((udata->updatecode == 0) && 
        ((udata->num_found + in_flight_rels) < udata->num_needed))
    {
        // generate a new poly A value for the thread we pulled out of the queue
        // using its dconf.  this is done by the master thread because it also 
        // stores the coefficients in a master list
        static_conf->total_poly_a++;
        new_poly_a(static_conf, t[tid].dconf);
        tdata->work_fcn_id = 0;
    }
    else
    {
        static_conf->flag = 1;
        tdata->work_fcn_id = tdata->num_work_fcn;        
    }


    return;
}

void SIQS(fact_obj_t *fobj)
{
	// the input fobj->N and this 'n' are pointers to memory which holds
	// the input number to this function.  a copy in different memory
	// will also be made, and modified by a multiplier.  don't confuse
	// these two.
	// input expected in fobj->qs_obj->gmp_n	

	// thread data holds all data needed during sieving
    tpool_t *tpool_data;
    siqs_userdata_t udata;
	thread_sievedata_t *thread_data;		//an array of thread data objects

	// master control structure
	static_conf_t *static_conf;

	// stuff for lanczos
	qs_la_col_t *cycle_list;
	uint32 num_cycles = 0;
	uint64 *bitfield = NULL;
	siqs_r *relation_list;

	// log file
	FILE *sieve_log;

	// some locals	
	int i;
	clock_t start, stop;
	double t_time;
	struct timeval myTVend;

	// checking savefile
	FILE *data;
	char tmpstr[GSTR_MAXSIZE];

	// logfile for this factorization
	// must ensure it is only written to by main thread
	if (fobj->qs_obj.flags != 12345)
	{
		fobj->logfile = fopen(fobj->flogname,"a");
		sieve_log = fobj->logfile;
	}
	else
	{
		fobj->logfile = NULL;
		sieve_log = fobj->logfile;
	}

	// check for special cases and bail if there is one
	if ((i = check_specialcase(fobj->logfile,fobj)) > 0)
	{
		if (i == 1)
		{
			if (sieve_log != NULL)
				fclose(sieve_log);
		}
		return;
	}	

	// check to see if a siqs savefile exists for this input	
	data = fopen(fobj->qs_obj.siqs_savefile,"r");

	if (data != NULL)
	{	
		char *substr;
		mpz_t tmpz;
		mpz_t g;

		//read in the number from the savefile
		mpz_init(tmpz);
		mpz_init(g);

		fgets(tmpstr, GSTR_MAXSIZE, data);
		substr = tmpstr + 2;
		mpz_set_str(tmpz, substr, 0);	//auto detect the base

		if (resume_check_input_match(tmpz, fobj->qs_obj.gmp_n, g))
		{
			// remove any common factor so the input exactly matches
			// the file.  
            // TODO: we should also factor the common factor 'g' 
            // somehow, since the user may be expecting the complete
            // factorization of the input instead of just the 
            // qs-resumed bit.
            if (mpz_cmp_ui(g, 1) > 0)
            {
                add_to_factor_list(fobj, g);
            }
			mpz_tdiv_q(fobj->qs_obj.gmp_n, fobj->qs_obj.gmp_n, g);
			mpz_set(fobj->N, fobj->qs_obj.gmp_n);
		}
		mpz_clear(tmpz);
		mpz_clear(g);
		fclose(data);
	}

	// then do a small amount of trial division
	// which will add anything found to the global factor list
	// this is only called with the main thread
	mpz_set(fobj->div_obj.gmp_n, fobj->qs_obj.gmp_n);
	fobj->div_obj.print = 0;
	fobj->div_obj.limit = 10000;
	zTrial(fobj);
	mpz_set(fobj->qs_obj.gmp_n, fobj->div_obj.gmp_n);

	// At this point, we are committed to doing qs on the input
	// we need to:
	// 1.) notify the screen and log
	// 2.) start watching for an abort signal
	// 3.) initialize data objects
	// 4.) get ready to find the factor base

	// fill in the factorization object
	fobj->bits = mpz_sizeinbase(fobj->qs_obj.gmp_n, 2);
	fobj->digits = gmp_base10(fobj->qs_obj.gmp_n);
	fobj->qs_obj.savefile.name = (char *)malloc(80 * sizeof(char));
	strcpy(fobj->savefile_name,fobj->qs_obj.siqs_savefile);

	// initialize the data objects both shared (static) and 
    // per-thread (dynamic)
	static_conf = (static_conf_t *)malloc(sizeof(static_conf_t));
	static_conf->obj = fobj;

	thread_data = (thread_sievedata_t *)malloc(THREADS * sizeof(thread_sievedata_t));
	for (i=0; i<THREADS; i++)
	{
		thread_data[i].dconf = (dynamic_conf_t *)malloc(sizeof(dynamic_conf_t));
		thread_data[i].sconf = static_conf;
	}

	// initialize the flag to watch for interrupts, and set the
	// pointer to the function to call if we see a user interrupt
	SIQS_ABORT = 0;
	signal(SIGINT,siqsexit);

	// initialize global offset for savefile buffer
	savefile_buf_off = 0;	
		
	// start a counter for the whole job
	gettimeofday(&static_conf->totaltime_start, NULL);

    if (VFLAG >= 0)
    {
        printf("\nstarting SIQS on c%d: %s\n", fobj->digits,
            mpz_conv2str(&gstr1.s, 10, fobj->qs_obj.gmp_n));
    }

	if (sieve_log != NULL)
	{
		logprint(sieve_log,"starting SIQS on c%d: %s\n",fobj->digits,
			mpz_conv2str(&gstr1.s, 10, fobj->qs_obj.gmp_n));
		logprint(sieve_log,"random seeds: %u, %u\n",g_rand.hi, g_rand.low);
		fflush(sieve_log);
	}

    udata.thread_data = thread_data;

    tpool_data = tpool_setup(THREADS, NULL, NULL,
        &siqs_sync, &siqs_dispatch, &udata);
    tpool_add_work_fcn(tpool_data, &process_poly);

    // this function is not run by tpool. It puts things into the user
    // data portion of the tpool structure and runs a few initialization 
    // routines in an attempt to clean up this toplevel function.
    siqs_start(tpool_data);

    if (THREADS == 1)
    {
        // it is noticably faster to remove the tpool overhead
        // if we just have one thread.  This is basically what
        // tpool_go() does without all of the threading overhead.
        // todo: maybe this *should* be what tpool_go() does when
        // num_threads == 1...
        while (1)
        {
            siqs_sync(tpool_data);
            siqs_dispatch(tpool_data);
            if (tpool_data->work_fcn_id == 0)
            {
                process_poly(tpool_data);
            }
            else
            {
                break;
            }
        }
    }
    else
    {
        tpool_go(tpool_data);
    }

    // this just holds pointers to other stuff that still exists.
    // we can get rid of the wrapper now...
    free(tpool_data);

	//stop worker threads
    for (i = 0; i < THREADS; i++)
    {
        free_sieve(thread_data[i].dconf);
        free(thread_data[i].dconf->relation_buf);
    }
	
	//finialize savefile
	qs_savefile_flush(&static_conf->obj->qs_obj.savefile);
	qs_savefile_close(&static_conf->obj->qs_obj.savefile);		
	
	update_final(static_conf);

    if (udata.updatecode == 2)
    {
        goto done;
    }

	//we don't need the poly_a_list anymore... free it so the other routines
	//can use it (unless we are doing in-mem)
	if (!static_conf->in_mem)
	{
		for (i=0;i<static_conf->total_poly_a + 1;i++)
			mpz_clear(static_conf->poly_a_list[i]);
		free(static_conf->poly_a_list);
	}
	
	gettimeofday (&myTVend, NULL);
    t_time = my_difftime(&static_conf->totaltime_start, &myTVend);

	if (VFLAG > 0)
	{
		printf("QS elapsed time = %6.4f seconds.\n",t_time);
		printf("\n==== post processing stage (msieve-1.38) ====\n");
	}

#if defined (TARGET_KNC) || defined(TARGET_KNL)
    // for now, just do timing on the sieving portion.  LA is really slow.
    exit(1);
#endif

	fobj->qs_obj.qs_time = t_time;	
	
	start = clock();

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
	yafu_qs_filter_relations(thread_data[0].sconf);

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
					
        free(bitfield);
	}

	stop = clock();
	static_conf->t_time2= (double)(stop - start)/(double)CLOCKS_PER_SEC;

	gettimeofday (&myTVend, NULL);
    static_conf->t_time3 = my_difftime(&static_conf->totaltime_start, &myTVend);
	
	if (VFLAG > 0)
	{
		printf("Lanczos elapsed time = %6.4f seconds.\n",static_conf->t_time1);
		printf("Sqrt elapsed time = %6.4f seconds.\n",static_conf->t_time2);
		
	}

	if (VFLAG >= 0)
		printf("SIQS elapsed time = %6.4f seconds.\n",static_conf->t_time3);

	fobj->qs_obj.total_time = static_conf->t_time3;

	if (sieve_log != NULL)
	{
		logprint(sieve_log,"Lanczos elapsed time = %6.4f seconds.\n",static_conf->t_time1);
		logprint(sieve_log,"Sqrt elapsed time = %6.4f seconds.\n",static_conf->t_time2);
		logprint(sieve_log,"SIQS elapsed time = %6.4f seconds.\n",static_conf->t_time3);
		logprint(sieve_log,"\n");
		logprint(sieve_log,"\n");
	}

	static_conf->cycle_list = cycle_list;
	static_conf->num_cycles = num_cycles;

	//free stuff used during filtering
	free_filter_vars(static_conf);

done:

	//free everything else
	free_siqs(thread_data[0].sconf);

    if (sieve_log != NULL)
    {
        fclose(sieve_log);
    }

	for (i=0; i<THREADS; i++)
	{
		free(thread_data[i].dconf);
	}
	free(static_conf);
	free(thread_data);

	//reset signal handler to default (no handler).
	signal(SIGINT,NULL);

	return;
}

void *process_poly(void *vptr)
{
    // top-level sieving function which performs all work for a single
    // new a coefficient.  has pthread calling conventions, meant to be
    // used in a multi-threaded environment
    tpool_t *tdata = (tpool_t *)vptr;
    siqs_userdata_t *udata = tdata->user_data;
    thread_sievedata_t *thread_data = &udata->thread_data[tdata->tindex];
    static_conf_t *sconf = thread_data->sconf;
    dynamic_conf_t *dconf = thread_data->dconf;

    // unpack stuff from the job data structure
    sieve_fb_compressed *fb_sieve_p = dconf->comp_sieve_p;
    sieve_fb_compressed *fb_sieve_n = dconf->comp_sieve_n;
    siqs_poly *poly = dconf->curr_poly;
    uint8 *sieve = dconf->sieve;
    fb_list *fb = sconf->factor_base;
    lp_bucket *buckets = dconf->buckets;
    uint32 start_prime = sconf->sieve_small_fb_start;
    uint32 num_blocks = sconf->num_blocks;
    uint8 blockinit = sconf->blockinit;

    // locals
    uint32 i;

    // to get relations per second
    double t_time;
    struct timeval start, stop, st;

    // this routine is handed a dconf structure which already has a
    // new poly a coefficient (and some supporting data).  continue from
    // there, first initializing the gray code...
    gettimeofday(&start, NULL);

    // used to print a little more status info for huge jobs.
    if (sconf->digits_n > 110)
        gettimeofday(&st, NULL);

    // update the gray code
    get_gray_code(dconf->curr_poly);

    // update roots, etc.
    dconf->maxB = 1 << (dconf->curr_poly->s - 1);
    dconf->numB = 1;
    computeBl(sconf, dconf);

    //printf("\n");
    //for (i = 1; i < dconf->maxB; i++)
    //{
    //    printf("poly = %d, nu = %d, sign = %d\n", i,
    //        dconf->curr_poly->nu[i], dconf->curr_poly->gray[i]);
    //}
    //printf("\n");


    firstRoots_ptr(sconf, dconf);

    // loop over each possible b value, for the current a value
    for (; dconf->numB < dconf->maxB; dconf->numB++, dconf->tot_poly++)
    {
        uint32 invalid_root_marker = 0xFFFFFFFF;

        for (i = 0; i < num_blocks; i++)
        {
            // set the roots for the factors of a such that
            // they will not be sieved.  we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p, 1);
            med_sieve_ptr(sieve, fb_sieve_p, fb, start_prime, blockinit);
            lp_sieveblock(sieve, i, num_blocks, buckets, 0, dconf);

            // set the roots for the factors of a to force the following routine
            // to explicitly trial divide since we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p, 0);
            scan_ptr(i, 0, sconf, dconf);

            // set the roots for the factors of a such that
            // they will not be sieved.  we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n, 1);
            med_sieve_ptr(sieve, fb_sieve_n, fb, start_prime, blockinit);
            lp_sieveblock(sieve, i, num_blocks, buckets, 1, dconf);

            // set the roots for the factors of a to force the following routine
            // to explicitly trial divide since we haven't found roots for them
            set_aprime_roots(sconf, invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n, 0);
            scan_ptr(i, 1, sconf, dconf);
        }

        // print a little more status info for huge jobs.
        if (THREADS < 32)
        {
            if (sconf->digits_n > 110)
            {
                gettimeofday(&stop, NULL);
                t_time = my_difftime(&st, &stop);

                if (t_time > 5)
                {
                    // print some status
                    gettimeofday(&stop, NULL);
                    t_time = my_difftime(&start, &stop);

                    printf("Bpoly %u of %u: buffered %u rels, checked %u (%1.2f rels/sec)\n",
                        dconf->numB, dconf->maxB, dconf->buffered_rels, dconf->num,
                        (double)dconf->buffered_rels / t_time);

                    // reset the timer
                    gettimeofday(&st, NULL);
                }
            }
        }

        // next polynomial
        // use the stored Bl's and the gray code to find the next b
        nextB(dconf, sconf);

#ifdef USE_BATCHPOLY

#ifdef TARGET_KNC
        nextRoots_32k_knc_small(sconf, dconf);

        // every N iterations we do the bucket sieve
        if ((dconf->numB % dconf->poly_batchsize) == 1)
        {
            nextRoots_32k_knc_polybatch(sconf, dconf);
        }
#elif defined(USE_AVX512F)
        // every iteration we update the small-med prime's roots
        nextRoots_32k_avx2_small(sconf, dconf);

        // every N iterations we do the bucket sieve
        if ((dconf->numB % dconf->poly_batchsize) == 1)
        {
            nextRoots_32k_knl_polybatch(sconf, dconf);
        }

#else
        // every iteration we update the small-med prime's roots
        nextRoots_32k_generic_small(sconf, dconf);

        // every N iterations we do the bucket sieve
        if ((dconf->numB % dconf->poly_batchsize) == 1)
        {
            nextRoots_32k_generic_polybatch(sconf, dconf);
        }
#endif

#else

#ifdef TARGET_KNC

        nextRoots_32k_knc_small(sconf, dconf);
        nextRoots_32k_knc_bucket(sconf, dconf);

#elif USE_AVX512F

        nextRoots_32k_avx2_small(sconf, dconf);
        nextRoots_32k_knl_bucket(sconf, dconf);

#else
        // and update the roots
        nextRoots_ptr(sconf, dconf);
#endif

#endif

        if (THREADS >= 32)
        {
            // other threads may have got us past the threshold while we
            // are still working on this 'A' poly... check if we can stop.
            if ((sconf->num_found > sconf->num_needed) || (sconf->flag == 1))
            {
                break;
            }
        }
    }



#ifdef USE_VEC_SQUFOF
    // vector SQUFOF if necessary
    if ((sconf->use_dlp) && (dconf->num_64bit_residue > 0))
    {
        uint64 *f = dconf->residue_factors;
        uint32 f32;
        int j = 0;
        siqs_r *rel;

        //printf("attempting vector squfof on %d residues... ", dconf->num_64bit_residue);
#if defined(__INTEL_COMPILER)
        f32 = par_shanks_loop(dconf->unfactored_residue, f, 
            dconf->num_64bit_residue);
#else
        for (i = 0; i < dconf->num_64bit_residue; i++)
        {
            f[i] = spbrent(dconf->unfactored_residue[i], 1, 1024);
        }
#endif

        //printf("vector squfof reported %d successes\n", f32);

        dconf->attempted_squfof += dconf->num_64bit_residue;
        for (i=0; i < dconf->buffered_rels; i++)
        {
            rel = dconf->relation_buf + i;

            // buffered_rels contains non-dlp relations as well, so check
            // if this one was a dlp.  Only dlp relations will have
            // been sent to squfof, and they appear in the factor list
            // in the same order as in the relation buffer.
            if (rel->large_prime[0] == 0xffffffff)
            {                
                // get the next factorization
                f32 = f[j];

                //printf("relation %d found factor %u of input %lu in position %d\n", 
                //    i, f32, dconf->unfactored_residue[j], j);

                if (f32 > 1)
                {
                    rel->large_prime[0] = f32;
                    rel->large_prime[1] = dconf->unfactored_residue[j] / f32;

                    if ((rel->large_prime[0] < sconf->large_prime_max) &&
                        (rel->large_prime[1] < sconf->large_prime_max))
                    {
                        //add this one
                        dconf->dlp_useful++;
                    }
                    else
                    {
                        // mark it as failed so we don't write it to the savefile
                        rel->large_prime[0] = 0xffffffff;
                    }

                }
                else
                {
                    dconf->failed_squfof++;

                    // mark it as failed so we don't write it to the savefile
                    rel->large_prime[0] = 0xffffffff;
                }

                j++;
            }
        }
    }
#endif



	gettimeofday (&stop, NULL);
    t_time = my_difftime(&start, &stop);

	dconf->rels_per_sec = (double)dconf->buffered_rels / t_time;

	//printf("average utilization of buckets in slices\n");
	//for (i=0; i<20; i++)
	//	printf("%d: %1.1f ",i,average_primes_per_slice[i]);
	//printf("\n");

	//unlock_thread_from_core();
	return 0;
}

uint32 siqs_merge_data(dynamic_conf_t *dconf, static_conf_t *sconf)
{
	//the sconf structure holds the master list of relations and cycles.
	//merge everything we found in the last round of sieving into
	//this table and save relations out to disk
	uint32 i;
	siqs_r *rel;
	char buf[1024];

	// save the A value.
	if (!sconf->in_mem)
	{
		gmp_sprintf(buf,"A 0x%Zx\n", dconf->curr_poly->mpz_poly_a);
		qs_savefile_write_line(&sconf->obj->qs_obj.savefile,buf);
	}

#ifdef TARGET_KNC
    // spam the screen
    //printf("saving %d buffered relations\n", dconf->buffered_rels);
#endif

	//save the data and merge into master cycle structure
	for (i=0; i<dconf->buffered_rels; i++)
	{
		rel = dconf->relation_buf + i;

#ifdef USE_VEC_SQUFOF
        // rarely, squfof will fail to factor a dlp.  
        // just skip these.
        if (rel->large_prime[0] == 0xffffffff)
        {
            continue;
        }
#endif
		save_relation_siqs(rel->sieve_offset,rel->large_prime,
			rel->num_factors, rel->fb_offsets, rel->poly_idx, 
			rel->parity, sconf);
	}

	//update some progress indicators
    sconf->total_blocks += dconf->total_blocks;
    sconf->total_reports += dconf->total_reports;
    sconf->total_surviving_reports += dconf->total_surviving_reports;
	sconf->num += dconf->num;
	sconf->tot_poly += dconf->tot_poly;
	sconf->failed_squfof += dconf->failed_squfof;
	sconf->attempted_squfof += dconf->attempted_squfof;
	sconf->dlp_outside_range += dconf->dlp_outside_range;
	sconf->dlp_prp += dconf->dlp_prp;
	sconf->dlp_useful += dconf->dlp_useful;
    sconf->lp_scan_failures += dconf->lp_scan_failures;

	//compute total relations found so far
	sconf->num_r = sconf->num_relations + 
		sconf->num_cycles +
		sconf->components - sconf->vertices;

	return sconf->num_r;
}

int siqs_check_restart(dynamic_conf_t *dconf, static_conf_t *sconf)
{
	fact_obj_t *obj = sconf->obj;
	char buf[1024];
	int state = 0;
    sconf->is_restart = 0;

	// if we want to do an in-memory factorization, then 
	// ignore the current state of the savefile and don't
	// prepare it for use.
	if (sconf->in_mem)
		return 0;

	//we're now almost ready to start, but first
	//check if this number has had work done
	restart_siqs(sconf,dconf);
	
	if ((uint32)sconf->num_r >= (sconf->factor_base->B + sconf->num_extra_relations)) 
	{
		//we've got enough total relations to stop		
		qs_savefile_open(&obj->qs_obj.savefile,SAVEFILE_APPEND);	
		dconf->buckets->list = NULL;
		//signal that we should proceed to post-processing
		state = 1;
	}
	else if (sconf->num_r > 0)
	{
		//we've got some relations, but not enough to finish.
		//whether or not this is a big job, it needed to be resumed
		//once so treat it as if it will need to be again.  use the savefile.
		qs_savefile_open(&obj->qs_obj.savefile,SAVEFILE_APPEND);

        // don't try to do optimization of the small cutoff... it is not 
        // designed to cope with starting with a bunch of relations and poly_a's
        obj->qs_obj.no_small_cutoff_opt = 1;
        sconf->is_restart = 1;
	}
	else
	{
		//no relations found, get ready for new factorization
		//we'll be writing to the savefile as we go, so get it ready
		qs_savefile_open(&obj->qs_obj.savefile,SAVEFILE_WRITE);
		gmp_sprintf(buf,"N 0x%Zx\n", sconf->obj->qs_obj.gmp_n);
		qs_savefile_write_line(&obj->qs_obj.savefile,buf);
		qs_savefile_flush(&obj->qs_obj.savefile);
		qs_savefile_close(&obj->qs_obj.savefile);
		//and get ready for collecting relations
		qs_savefile_open(&obj->qs_obj.savefile,SAVEFILE_APPEND);
	}

	return state;
}

void print_siqs_splash(dynamic_conf_t *dconf, static_conf_t *sconf)
{
    //print some info to the screen and the log file
    char inst_set[16];

#ifdef TARGET_KNC

    strcpy(inst_set, "KNC");

#elif defined(USE_AVX2)
    if (HAS_AVX2)
    {
#if defined (USE_AVX512F)
        strcpy(inst_set, "AVX512");
#else
        strcpy(inst_set, "AVX2");
#endif
    }
    else if (HAS_SSE41)
    {
        strcpy(inst_set, "SSE4.1");
    }
    else if (HAS_SSE2)
    {
        strcpy(inst_set, "SSE2");
    }
    else
    {
        strcpy(inst_set, "generic C");
    }
#elif defined(USE_SSE41)
    if (HAS_SSE41)
    {
        strcpy(inst_set, "SSE4.1");
    }
    else if (HAS_SSE2)
    {
        strcpy(inst_set, "SSE2");
    }
    else
    {
        strcpy(inst_set, "generic C");
    }
#elif defined(HAS_SSE2)
    if (HAS_SSE2)
    {
        strcpy(inst_set, "SSE2");
    }
    else
    {
        strcpy(inst_set, "generic C");
    }
#else
    strcpy(inst_set, "generic C");
#endif

    if (VFLAG > 0)
    {
        printf("\n==== sieve params ====\n");
        printf("n = %d digits, %d bits\n", sconf->digits_n, sconf->bits);
        printf("factor base: %d primes (max prime = %u)\n", sconf->factor_base->B, sconf->pmax);
        printf("single large prime cutoff: %u (%d * pmax)\n",
            sconf->large_prime_max, sconf->large_mult);
        if (sconf->use_dlp)
        {
            printf("double large prime range from %d to %d bits\n",
                sconf->dlp_lower, sconf->dlp_upper);
            printf("double large prime range from %" PRIu64 " to %" PRIu64 "\n",
                sconf->max_fb2, sconf->large_prime_max2);
        }
        if (dconf->buckets->list != NULL)
        {
            printf("allocating %d large prime slices of factor base\n",
                dconf->buckets->alloc_slices);
            printf("buckets hold %d elements\n", BUCKET_ALLOC);
#ifdef USE_BATCHPOLY
            printf("processing polynomials in batches of %d\n", dconf->poly_batchsize);
#endif
            printf("large prime hashtables have %d bytes\n", 
                dconf->buckets->list_size * BUCKET_ALLOC * sizeof(uint32));
        }
        printf("using %s enabled 32k sieve core\n", inst_set);
        printf("sieve interval: %d blocks of size %d\n",
            sconf->num_blocks, sconf->qs_blocksize);
        printf("polynomial A has ~ %d factors\n",
            (int)mpz_sizeinbase(sconf->target_a, 2) / 11);
        printf("using multiplier of %u\n", sconf->multiplier);
        printf("using SPV correction of %d bits, starting at offset %d\n",
            sconf->tf_small_cutoff, sconf->sieve_small_fb_start);
		printf("trial factoring cutoff at %d bits\n",sconf->tf_closnuf);
	}

	if (VFLAG >= 0)
	{
		if (THREADS == 1)
		{
			printf("\n==== sieving in progress (1 thread): %7u relations needed ====\n",
				sconf->factor_base->B + sconf->num_extra_relations);
			printf(  "====           Press ctrl-c to abort and save state           ====\n");
		}
		else
		{
			printf("\n==== sieving in progress (%3d threads): %7u relations needed ====\n",
				THREADS,sconf->factor_base->B + sconf->num_extra_relations);
			printf(  "====             Press ctrl-c to abort and save state            ====\n");
		}
	}

	if (sconf->obj->logfile != NULL)
	{
		logprint(sconf->obj->logfile,"==== sieve params ====\n");
		logprint(sconf->obj->logfile,"n = %d digits, %d bits\n",
			sconf->digits_n,sconf->bits);
		logprint(sconf->obj->logfile,"factor base: %d primes (max prime = %u)\n",
			sconf->factor_base->B,sconf->pmax);
		logprint(sconf->obj->logfile,"single large prime cutoff: %u (%d * pmax)\n",
			sconf->large_prime_max,sconf->large_mult);
		if (sconf->use_dlp)
		{
			logprint(sconf->obj->logfile,"double large prime range from %d to %d bits\n",
				sconf->dlp_lower,sconf->dlp_upper);
			logprint(sconf->obj->logfile,"double large prime cutoff: %" PRIu64 "\n",
				sconf->large_prime_max2);
		}
		if (dconf->buckets->list != NULL)
		{
			logprint(sconf->obj->logfile,"allocating %d large prime slices of factor base\n",
				dconf->buckets->alloc_slices);
			logprint(sconf->obj->logfile,"buckets hold %d elements\n",BUCKET_ALLOC);
#ifdef USE_BATCHPOLY
            logprint(sconf->obj->logfile, "processing polynomials in batches of %d\n", 
                dconf->poly_batchsize);
#endif
            logprint(sconf->obj->logfile, "large prime hashtables have %d bytes\n", 
                dconf->buckets->list_size * BUCKET_ALLOC * sizeof(uint32));
		}
        logprint(sconf->obj->logfile,"using %s enabled 32k sieve core\n", inst_set);
		logprint(sconf->obj->logfile,"sieve interval: %d blocks of size %d\n",
			sconf->num_blocks,sconf->qs_blocksize);
		logprint(sconf->obj->logfile,"polynomial A has ~ %d factors\n", 
			mpz_sizeinbase(sconf->target_a, 2) / 11);
		logprint(sconf->obj->logfile,"using multiplier of %u\n",sconf->multiplier);
		logprint(sconf->obj->logfile,"using SPV correction of %d bits, starting at offset %d\n",
			sconf->tf_small_cutoff,sconf->sieve_small_fb_start);
		logprint(sconf->obj->logfile,"trial factoring cutoff at %d bits\n",
			sconf->tf_closnuf);
		if (THREADS == 1)
			logprint(sconf->obj->logfile,"==== sieving started (1 thread) ====\n");
		else
			logprint(sconf->obj->logfile,"==== sieving started (%2d threads) ====\n",THREADS);
	}

	return;
}

int siqs_dynamic_init(dynamic_conf_t *dconf, static_conf_t *sconf)
{
	//allocate the dynamic structure which hold scratch space for 
	//various things used during sieving.
	uint32 i, memsize;

	if (VFLAG > 2)
	{
		printf("memory usage during sieving:\n");
	}

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
	if (VFLAG > 2)
	{
		memsize = MAX_A_FACTORS * sizeof(int);
		memsize += 65536 * sizeof(char);
		memsize += 65536 * sizeof(char);
		printf("\tcurr_poly structure: %d bytes\n",memsize);
	}

	//initialize temporary storage of relations
	dconf->relation_buf = (siqs_r *)malloc(32768 * sizeof(siqs_r));
	dconf->buffered_rel_alloc = 32768;
	dconf->buffered_rels = 0;

	if (VFLAG > 2)
	{
		memsize = 32768 * sizeof(siqs_r);
		printf("\trelation buffer: %d bytes\n",memsize);
	}

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
        (size_t)(MAX_A_FACTORS * sconf->factor_base->med_B * sizeof(uint16)));


	if (VFLAG > 2)
	{
		memsize = sconf->factor_base->med_B * sizeof(sieve_fb_compressed) * 2;
		printf("\tfactor bases: %d bytes\n",memsize);
	}

	if (VFLAG > 2)
	{
		memsize = sconf->factor_base->B * sizeof(int)  * 2;
		memsize += sconf->factor_base->B * sizeof(uint32);
		memsize += sconf->factor_base->B * sizeof(uint8);
		printf("\tupdate data: %d bytes\n",memsize);
	}
	
	// allocate the sieve with a guard region.  this lets
    // us be a little more sloppy with vector-based sieving
    // and resieving.
	dconf->sieve = (uint8 *)xmalloc_align(
		(size_t) (sconf->qs_blocksize * 4 * sizeof(uint8)));
    dconf->sieve = &dconf->sieve[2 * sconf->qs_blocksize];

	if (VFLAG > 2)
	{
		memsize = sconf->qs_blocksize * sizeof(uint8);
		printf("\tsieve: %d bytes\n",memsize);
	}

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

		dconf->update_data.prime[i] = p;
		dconf->update_data.logp[i] = lp;
	}

	for (; i < sconf->factor_base->B; i++)
	{
		dconf->update_data.prime[i] = sconf->factor_base->list->prime[i];
		dconf->update_data.logp[i] = sconf->factor_base->list->logprime[i];
	}

	//check if we should use bucket sieving, and allocate structures if so
	if (sconf->factor_base->B > sconf->factor_base->med_B)
	{
		dconf->buckets = (lp_bucket *)malloc(sizeof(lp_bucket));

		//test to see how many slices we'll need.
		testRoots_ptr(sconf,dconf);

#ifdef USE_BATCHPOLY
        // this should be a function of the L2 size.
#ifdef TARGET_KNC
        dconf->poly_batchsize = MAX(2, 30000000 / 4 /
            (2 * sconf->num_blocks * dconf->buckets->alloc_slices *
            BUCKET_ALLOC * sizeof(uint32)));
#else
        dconf->poly_batchsize = MAX(2, L2CACHE / 2 / THREADS / 
            (2 * sconf->num_blocks * dconf->buckets->alloc_slices *
            BUCKET_ALLOC * sizeof(uint32)));
#endif

        dconf->poly_batchsize = 8;

#else
        dconf->poly_batchsize = 1;
#endif

		//initialize the bucket lists and auxilary info.
        dconf->buckets->num = (uint32 *)xmalloc_align(dconf->poly_batchsize * 
			2 * sconf->num_blocks * dconf->buckets->alloc_slices * sizeof(uint32));
		dconf->buckets->fb_bounds = (uint32 *)malloc(
			dconf->buckets->alloc_slices * sizeof(uint32));
		dconf->buckets->logp = (uint8 *)calloc(
			dconf->buckets->alloc_slices, sizeof(uint8));
        dconf->buckets->list_size = dconf->poly_batchsize * 2 * 
            sconf->num_blocks * dconf->buckets->alloc_slices;
		
		//now allocate the buckets
        dconf->buckets->list = (uint32 *)xmalloc_align(
            dconf->buckets->list_size *
			BUCKET_ALLOC * sizeof(uint32));
	}
	else
	{
        dconf->poly_batchsize = 1;
		dconf->buckets = (lp_bucket *)malloc(sizeof(lp_bucket));
		dconf->buckets->list = NULL;
		dconf->buckets->alloc_slices = 0;
		dconf->buckets->num_slices = 0;
	}

	if (VFLAG > 2)
	{
        memsize = dconf->buckets->list_size * sizeof(uint32);
		memsize += dconf->buckets->alloc_slices * sizeof(uint32);
		memsize += dconf->buckets->alloc_slices * sizeof(uint8);
        memsize += dconf->buckets->list_size * BUCKET_ALLOC * sizeof(uint32);
		printf("\tbucket data: %d bytes\n",memsize);
	}

	//used in trial division to mask out the fb_index portion of bucket entries, so that
	//multiple block locations can be searched for in parallel using SSE2 instructions
	dconf->mask = (uint16 *)xmalloc_align(16 * sizeof(uint16));
    dconf->mask2 = (uint32 *)xmalloc_align(8 * sizeof(uint32));

	dconf->mask[1] = 0xFFFF;
	dconf->mask[3] = 0xFFFF;
	dconf->mask[5] = 0xFFFF;
	dconf->mask[7] = 0xFFFF;
	dconf->mask[9] = 0xFFFF;
	dconf->mask[11] = 0xFFFF;
	dconf->mask[13] = 0xFFFF;
	dconf->mask[15] = 0xFFFF;

    dconf->mask2[0] = 0x0000ffff;
    dconf->mask2[1] = 0x0000ffff;
    dconf->mask2[2] = 0x0000ffff;
    dconf->mask2[3] = 0x0000ffff;
    dconf->mask2[4] = 0x0000ffff;
    dconf->mask2[5] = 0x0000ffff;
    dconf->mask2[6] = 0x0000ffff;
    dconf->mask2[7] = 0x0000ffff;


	// used in SIMD optimized resiever
	dconf->corrections = (uint16 *)xmalloc_align(16 * sizeof(uint16));

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
		mpz_init(dconf->Qvals[i]); //, 2*sconf->bits);

	dconf->valid_Qs = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));
	dconf->smooth_num = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));
	dconf->failed_squfof = 0;
	dconf->attempted_squfof = 0;
	dconf->dlp_outside_range = 0;
	dconf->dlp_prp = 0;
	dconf->dlp_useful = 0;

#ifdef USE_VEC_SQUFOF
    dconf->unfactored_residue = (uint64 *)malloc(4096 * sizeof(uint64));
    dconf->residue_factors = (uint64 *)malloc(4096 * sizeof(uint64));
    dconf->num_64bit_residue = 0;
#endif

	// used in SIMD optimized versions of tdiv_med
#ifdef USE_8X_MOD_ASM
	dconf->bl_sizes = (uint16 *)xmalloc_align(16 * sizeof(uint16));
	dconf->bl_locs = (uint16 *)xmalloc_align(16 * sizeof(uint16));
#endif

    dconf->polyscratch = (uint32 *)xmalloc_align(16 * sizeof(uint32));

	//initialize some counters
	dconf->tot_poly = 0;		//track total number of polys
	dconf->num = 0;				//sieve locations subjected to trial division
    dconf->lp_scan_failures = 0;
    dconf->total_reports = 0;
    dconf->total_surviving_reports = 0;
    dconf->total_blocks = 0;

	return 0;
}

int siqs_static_init(static_conf_t *sconf, int is_tiny)
{
	//find the best parameters, multiplier, and factor base
	//for the input.  This is an iterative job because we may find
	//a factor of the input while finding the factor base, in which case
	//we log that fact, and start over with new parameters, multiplier etc.
	//also allocate space for things that don't change during the sieving
	//process.  this scratch space is shared among all threads.

	fact_obj_t *obj = sconf->obj;
	uint32 i, memsize;
	uint32 closnuf;
	double sum, avg, sd;
    uint32 dlp_cutoff;

    // this pretty much has to stay "8".  the reason is that many of the specialized routines
    // have picky requirements about how large or small the primes can be for them
    // to work correctly.  It seems we have trouble finding indices within the factor
    // base that simultaneously provide the correct size boundaries and are also
    // aligned to a multiple of 16.
    // When using AVX2, we'll maybe have to perform a round or two of SSE41 in a couple
    // places to cross over size boundaries within the factor base, if those boundaries
    // are not divisible by 16.
	int nump = 8;

	if (VFLAG > 2)
	{
		printf("static memory usage:\n");
	}

    


	// some things work different if the input is tiny
	sconf->is_tiny = is_tiny;

	//default parameters
	sconf->fudge_factor = 1.3;
	sconf->large_mult = 30;
	sconf->num_blocks = 40;
	sconf->num_extra_relations = 64;
	sconf->small_limit = 256;
	sconf->use_dlp = 0;

    // function pointer to the sieve array scanner
    scan_ptr = NULL;


    // relevant:
    // https://software.intel.com/en-us/node/523363
    // https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/X86-Built-in-Functions.html

	// sieve core functions
    firstRoots_ptr = &firstRoots_32k;
    nextRoots_ptr = &nextRoots_32k;

    // if the yafu library was both compiled with SSE41 code (USE_SSE41), and the user's 
    // machine has SSE41 instructions (HAS_SSE41), then proceed with 4.1.

#if defined(TARGET_KNC)
    nextRoots_ptr = &nextRoots_32k_knc_small;

#elif defined(USE_AVX2)
    if (HAS_AVX2)
    {
        nextRoots_ptr = &nextRoots_32k_avx2;
    }
    else if (HAS_SSE41)
    {
        nextRoots_ptr = &nextRoots_32k_sse41;
    }

#elif defined(USE_SSE41)
	if (HAS_SSE41)
	{
		nextRoots_ptr = &nextRoots_32k_sse41;
	}

#endif
		
	testRoots_ptr = &testfirstRoots_32k;

    med_sieve_ptr = &med_sieveblock_32k;
	// if the yafu library was both compiled with AVX2 code (USE_AVX2), and the user's 
	// machine has AVX2 instructions (HAS_AVX2), then proceed with AVX2.
#if defined(USE_AVX2)
	if (HAS_AVX2)
	{
		med_sieve_ptr = &med_sieveblock_32k_avx2;
	}
	else if (HAS_SSE41)
	{
		med_sieve_ptr = &med_sieveblock_32k_sse41;
	}

#elif defined(USE_SSE41)
	if (HAS_SSE41)
	{
		med_sieve_ptr = &med_sieveblock_32k_sse41;
	}

#endif

    tdiv_med_ptr = &tdiv_medprimes_32k;
    resieve_med_ptr = &resieve_medprimes_32k;

#if defined(USE_AVX2)
	if (HAS_AVX2)
	{
		tdiv_med_ptr = &tdiv_medprimes_32k_avx2;
		resieve_med_ptr = &resieve_medprimes_32k_avx2;
	}
#endif


#if defined(TARGET_KNC)
    tdiv_med_ptr = &tdiv_medprimes_32k_knc;
    resieve_med_ptr = &resieve_medprimes_32k_knc;
#endif

#if defined(USE_AVX512F)
    //tdiv_med_ptr = &tdiv_medprimes_32k_knl;
    resieve_med_ptr = &resieve_medprimes_32k_avx2;
#endif
		
	sconf->qs_blocksize = 32768;
	sconf->qs_blockbits = 15;

	//allocate the space for the factor base structure
	sconf->factor_base = (fb_list *)malloc(sizeof(fb_list));

	//allocate space for a copy of input number in the job structure
	mpz_init(sconf->n); //zInit(&sconf->n);
	mpz_set(sconf->n, obj->qs_obj.gmp_n);

	//initialize some constants
	mpz_init(sconf->sqrt_n);
	mpz_init(sconf->target_a);

	//initialize the bookkeeping for tracking partial relations
	sconf->components = 0;
	sconf->vertices = 0;
	sconf->num_cycles = 0;
	sconf->num_relations = 0;
	//force this to happen for now, eventually should implement this flag
	if (1 || !(sconf->obj->flags & MSIEVE_FLAG_SKIP_QS_CYCLES)) {
		sconf->cycle_hashtable = (uint32 *)xcalloc(
					(size_t)(1 << QS_LOG2_CYCLE_HASH),
					sizeof(uint32));
		sconf->cycle_table_size = 1;
		sconf->cycle_table_alloc = 10000;
		sconf->cycle_table = (qs_cycle_t *)xmalloc(
			sconf->cycle_table_alloc * sizeof(qs_cycle_t));
	}

	if (VFLAG > 2)
	{
		memsize = (1 << QS_LOG2_CYCLE_HASH) * sizeof(uint32);
		printf("\tinitial cycle hashtable: %d bytes\n",memsize);
		memsize = sconf->cycle_table_alloc * sizeof(qs_cycle_t);
		printf("\tinitial cycle table: %d bytes\n",memsize);
	}

	while (1)
	{
		//look up some parameters - tuned for 64k L1 cache
		get_params(sconf);

		//the number of primes in the factor base is rounded up to the 
		//next multiple of 16, so that we can use aligned moves to speed
		//computations of root updates
		sconf->factor_base->B += (16 - (sconf->factor_base->B % 16));

		//allocate the space for the factor base elements
		sconf->factor_base->list = (fb_element_siqs *)xmalloc_align(
			(size_t)(sizeof(fb_element_siqs)));
		sconf->factor_base->tinylist = (tiny_fb_element_siqs *)xmalloc_align(
			(size_t)(sizeof(tiny_fb_element_siqs)));

		sconf->modsqrt_array = (uint32 *)xmalloc_align(
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

#ifdef USE_8X_MOD_ASM
		sconf->factor_base->list->small_inv = (uint16 *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint16)));
		sconf->factor_base->list->correction = (uint16 *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint16)));
#else
		sconf->factor_base->list->small_inv = (uint32 *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint32)));
		sconf->factor_base->list->correction = (uint32 *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint32)));
#endif
		
		sconf->factor_base->list->logprime = (uint32 *)xmalloc_align(
			(size_t)(sconf->factor_base->B * sizeof(uint32)));


		if (VFLAG > 2)
		{
			memsize = sconf->factor_base->B * sizeof(uint32) * 5;
			printf("\tfactor base: %d bytes\n",memsize);
		}

		//find multiplier
		sconf->multiplier = (uint32)choose_multiplier_siqs(sconf->factor_base->B, sconf->n);
		mpz_mul_ui(sconf->n, sconf->n, sconf->multiplier);

		//sconf holds n*mul, so update its digit count and number of bits
		sconf->digits_n = gmp_base10(sconf->n);
		sconf->bits = mpz_sizeinbase(sconf->n, 2);

		//find sqrt_n
		mpz_sqrt(sconf->sqrt_n, sconf->n);

		//construct the factor base - divide out any primes which 
		//factor n other than mul
		//and go back to get_params if any are found
		if ((i = make_fb_siqs(sconf)) == 0)
			break;
		else
		{
			mpz_t tmpz;
			mpz_init(tmpz);

			//i is already divided out of n.  record the factor we found
			uint64_2gmp(i, tmpz);
			add_to_factor_list(sconf->obj, tmpz);
			mpz_clear(tmpz);

			//and remove the multiplier we may have added, so that
			//we can try again to build a factor base.
			mpz_tdiv_q_ui(sconf->n, sconf->n, sconf->multiplier);
			free(sconf->modsqrt_array);
			align_free(sconf->factor_base->list->prime);
			align_free(sconf->factor_base->list->small_inv);
			align_free(sconf->factor_base->list->correction);
			align_free(sconf->factor_base->list->logprime);
			align_free(sconf->factor_base->list);	
		}
	}

#ifdef VDEBUG
		printf("found factor base of %u elements; max = %u\n",
			sconf->factor_base->B,sconf->factor_base->list->prime[sconf->factor_base->B-1]);
#endif

	//the first two fb primes are always the same
	sconf->factor_base->list->prime[0] = 1;		//actually represents -1
	sconf->factor_base->list->prime[1] = 2;

    sconf->sieve_primes = (uint32 *)xmalloc_align(
        (size_t)(sconf->factor_base->B * sizeof(uint32)));

    for (i = 2; i < sconf->factor_base->B; i++)
    {
        sconf->sieve_primes[i] = sconf->factor_base->list->prime[i];
    }

	//adjust for various architectures
	if (sconf->qs_blocksize == 32768)
		sconf->num_blocks *= 2;
	else if (sconf->qs_blocksize == 65536)
		sconf->num_blocks *= 1;
	else
	{
		printf("unknown block size!\n");
		exit(1);
	}

	if (sconf->num_blocks < 1)
		sconf->num_blocks = 1;

	//and adjust for really small jobs
	if (sconf->bits <= 140)
		sconf->num_blocks = 1;

	//set the sieve interval.  this many blocks on each side of 0
	sconf->sieve_interval = sconf->qs_blocksize * sconf->num_blocks;

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

    // ======================================================================
    // set thresholds at various bit levels.  Note: many places in the code
    // are extremely sensitive to these thresholds; the program can hang or
    // crash unexpectedly if set incorrectly.  Adjust carefully!
    // ======================================================================

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 10 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 1024)  &&
			(i % nump == 0)) break;
	}
	sconf->factor_base->fb_10bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 11 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 2048)  &&
			(i % nump == 0)) break;
	}
	sconf->factor_base->fb_11bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 12 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
        if ((sconf->factor_base->list->prime[i] > 4096) &&
            (i % nump == 0)) break;
	}
	sconf->factor_base->fb_12bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
        if ((sconf->factor_base->list->prime[i] > 8192) &&
            (i % nump == 0)) break;
	}
	sconf->factor_base->fb_13bit_B = i;

	// the prime at which we want to start testing the asm code below is the point
	// at which we have at most 3 increments of the root.  this point starts when
	// prime > 32768 / 3 = 10922.
	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 13 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		//don't let med_B grow larger than 1.5 * the blocksize
		if ((sconf->factor_base->list->prime[i] > 10922)  &&
			(i % nump == 0)) break;
	}
	// put the upper bound just before prime exceeds blocksize/3, required by SSE2 sieving
	sconf->factor_base->fb_32k_div3 = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 14 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > 16384)  &&
            (i % nump == 0)) break;
	}	
	// put the upper bound just before prime exceeds 14 bits, required by SSE2 sieving
	sconf->factor_base->fb_14bit_B = i;

	for (; i < sconf->factor_base->B; i++)
	{
		//find the point at which factor base primes exceeds 15 bits.  
		//wait until the index is a multiple of 8 so that we can enter
		//this region of primes aligned on a 16 byte boundary and thus be able to use
		//movdqa
		if ((sconf->factor_base->list->prime[i] > 32768)  &&
            (i % nump == 0)) break;
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

	if (VFLAG > 1)
	{
		printf("fb bounds\n\tsmall: %u\n\tSPV: %u\n\t10bit: %u\n\t11bit: %u\n\t12bit: %u\n\t"
			"13bit: %u\n\t32k div 3: %u\n\t14bit: %u\n\t15bit: %u\n\tmed: %u\n\tlarge: %u\n"
            "\tlarge_x2: %u\n\tall: %u\n",
			sconf->factor_base->small_B,
			sconf->sieve_small_fb_start,
			sconf->factor_base->fb_10bit_B,
			sconf->factor_base->fb_11bit_B,
			sconf->factor_base->fb_12bit_B,
			sconf->factor_base->fb_13bit_B,
			sconf->factor_base->fb_32k_div3,
			sconf->factor_base->fb_14bit_B,
			sconf->factor_base->fb_15bit_B,
			sconf->factor_base->med_B,
			sconf->factor_base->large_B,
            sconf->factor_base->x2_large_B,
			sconf->factor_base->B);

		printf("start primes\n\tSPV: %u\n\t10bit: %u\n\t11bit: %u\n\t12bit: %u\n\t"
			"13bit: %u\n\t32k div 3: %u\n\t14bit: %u\n\t15bit: %u\n\tmed: %u\n\tlarge: %u\n"
            "\tlarge_x2: %u\n\tall: %u\n",
			sconf->factor_base->list->prime[sconf->sieve_small_fb_start - 1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_10bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_11bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_12bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_13bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_32k_div3-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_14bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_15bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->med_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->large_B-1],
            sconf->factor_base->list->prime[sconf->factor_base->x2_large_B - 1],
            sconf->factor_base->list->prime[sconf->factor_base->B]);
	}

	//a couple limits
	sconf->pmax = sconf->factor_base->list->prime[sconf->factor_base->B-1];

	if ((4294967295ULL / sconf->large_mult) < sconf->pmax)
	{
		// job is so big that pmax * default large_mult won't fit in 32 bits
		// reduce large_mult accordingly
		sconf->large_mult = 4294967295ULL / sconf->pmax;
		sconf->large_prime_max = sconf->pmax * sconf->large_mult;
	}
	else
		sconf->large_prime_max = sconf->pmax * sconf->large_mult;

	// based on the size of the input, determine how to proceed.
    // this should maybe be tuned based on machine type and/or other
    // factors as well, not just instruction set used during compile.
#if defined(USE_AVX2) || defined(USE_SSE41)
    dlp_cutoff = 70;
#else
    dlp_cutoff = 77;
#endif

    // could maybe someday change this w.r.t input size... for now
    // just test it out...
    //sconf->poly_batch_size = 1;

    if ((sconf->digits_n >= dlp_cutoff) || sconf->obj->qs_obj.gbl_force_DLP)
	{
		sconf->use_dlp = 1;
		scan_ptr = &check_relations_siqs_16;
		sconf->scan_unrolling = 128;
	}
	else
	{
		if (sconf->digits_n < 30)
		{
			scan_ptr = &check_relations_siqs_4;
			sconf->scan_unrolling = 32;
		}
		else if (sconf->digits_n < 60)
		{
			scan_ptr = &check_relations_siqs_4;
			sconf->scan_unrolling = 32;
		}
		else
		{
			scan_ptr = &check_relations_siqs_8;
			sconf->scan_unrolling = 64;
		}
		sconf->use_dlp = 0;
	}

#if defined(TARGET_KNC)
    // so far have only implemented this one
    scan_ptr = &check_relations_siqs_16;
    sconf->scan_unrolling = 128;
#endif

	qs_savefile_init(&obj->qs_obj.savefile, sconf->obj->qs_obj.siqs_savefile);

	//if we're using dlp, compute the range of residues which will
	//be subjected to factorization beyond trial division
	if (sconf->use_dlp)
	{
		sconf->max_fb2 = (uint64)sconf->pmax * (uint64)sconf->pmax;
		sconf->dlp_lower = spBits(sconf->max_fb2); 
		sconf->large_prime_max2 = (uint64)pow((double)sconf->large_prime_max,1.8);
		sconf->dlp_upper = spBits(sconf->large_prime_max2);
	}

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
	
    //if (sconf->digits_n > 60)
    //{
    //    closnuf -= 2;
    //}
	//it pays to trial divide a little more for big jobs...
	//"optimized" is used rather loosely here, but these corrections
	//were observed to make things faster.  
    //if ((sconf->digits_n >= 70) && (sconf->digits_n < 75))
    //{
    //    closnuf += 1;	//optimized for 75 digit num
    //}
    //else if ((sconf->digits_n >= 75) && (sconf->digits_n < 79))
    //{
    //    closnuf -= 2;	//optimized for 75 digit num
    //}
    //else if ((sconf->digits_n >= 78) && (sconf->digits_n < 81))
    //{
    //    closnuf -= 2;	//optimized for 80 digit num
    //}
    //else if ((sconf->digits_n >= 79) && (sconf->digits_n < 85))
    //{
    //    closnuf -= 3;	//optimized for 82 digit num	
    //}

    // inputs larger than this use a different method (below)


#ifdef QS_TIMING
	printf("%d primes not sieved in SPV\n",sconf->sieve_small_fb_start);
	printf("%d primes in small_B range\n",
		sconf->factor_base->med_B - sconf->sieve_small_fb_start);
	printf("%d primes in med_B range\n",
		sconf->factor_base->large_B - sconf->factor_base->med_B);
	printf("%d primes in large_B range\n",
		sconf->factor_base->B - sconf->factor_base->large_B);
	printf("detailing QS timing profiling enabled\n");

	TF_STG1 = 0;
	TF_STG2 = 0;
	TF_STG3 = 0;
	TF_STG4 = 0;
	TF_STG5 = 0;
	TF_STG6 = 0;
	TF_SPECIAL = 0;
	SIEVE_STG1 = 0;
	SIEVE_STG2 = 0;
	POLY_STG0 = 0;
	POLY_STG1 = 0;
	POLY_STG2 = 0;
	POLY_STG3 = 0;
	POLY_STG4 = 0;
	COUNT = 0;

#endif

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

    if ((sconf->use_dlp == 1) || sconf->obj->qs_obj.gbl_force_DLP)
    {
        //empirically, these were observed to work fairly well.
        if(sconf->digits_n < 82)
            closnuf = sconf->digits_n + 7;
        else if (sconf->digits_n < 90)
            closnuf = sconf->digits_n + 5;
        else if (sconf->digits_n < 95)
            closnuf = sconf->digits_n + 3;
        else if (sconf->digits_n < 100)
            closnuf = sconf->digits_n + 1;
        else
            closnuf = sconf->digits_n;
    }
    else
    {
        closnuf -= sconf->tf_small_cutoff;	//correction to the previous estimate
    }

	if (sconf->obj->qs_obj.gbl_override_tf_flag)
	{
		closnuf = sconf->obj->qs_obj.gbl_override_tf;
		printf("overriding with new closnuf = %d\n",closnuf);
	}

	// need the highest bit to be clear in order to scan the sieve array efficiently
	// this means we may check more sieve reports than we otherwise would have, but
	// only for the very largest of jobs which might actually benefit from doing so anyway.
	if (closnuf >= 128)
		closnuf = 127;

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

	//compute how often to check our list of partial relations and update the gui.
	sconf->check_inc = sconf->factor_base->B/10;
	sconf->check_total = sconf->check_inc;
	sconf->update_time = 5;
    sconf->flag = 0;

	//get ready for the factorization and screen updating
	sconf->t_update=0;
	sconf->last_numfull = 0;
	sconf->last_numpartial = 0;
	sconf->last_numcycles = 0;
	gettimeofday(&sconf->update_start, NULL);
	//sconf->update_start = clock();

	sconf->failed_squfof = 0;
	sconf->attempted_squfof = 0;
	sconf->dlp_outside_range = 0;
	sconf->dlp_prp = 0;
	sconf->dlp_useful = 0;
	sconf->total_poly_a = 0;	//track number of A polys used
	sconf->num_r = 0;			//total relations found
	sconf->charcount = 0;		//characters on the screen

	sconf->t_time1 = 0;			//sieve time
	sconf->t_time2 = 0;			//relation scanning and trial division
	sconf->t_time3 = 0;			//polynomial calculations and large prime sieving
	sconf->t_time4 = 0;			//extra?
	
	sconf->tot_poly = 0;		//track total number of polys
	sconf->num = 0;				//sieve locations subjected to trial division
    sconf->total_reports = 0;
    sconf->total_surviving_reports = 0;
    sconf->total_blocks = 0;
    sconf->lp_scan_failures = 0;

	//no factors so far...
	sconf->factor_list.num_factors = 0;

    // test: use in-memory relation storage below a certain digit level?
    // no.  maybe for tinysiqs someday.
    sconf->in_mem = 0;

	return 0;
}

int update_check(static_conf_t *sconf)
{
	//check to see if we should combine partial relations
	//found and check to see if we're done.  also update the screen.
	//this happens one of two ways, either we have found more than a 
	//certain amount of relations since the last time we checked, or if
	//a certain amount of time has elapsed.
	mpz_t tmp1;
	struct timeval update_stop;
	uint32 num_full = sconf->num_relations;
	uint32 check_total = sconf->check_total;
	uint32 check_inc = sconf->check_inc;
    double update_time = sconf->update_time;
	double t_update, t_time;	
	//int i;
	fb_list *fb = sconf->factor_base;
	int retcode = 0;

    // if we are working on a large problem or have a lot of
    // threads, update as fast as we can.
    if ((sconf->digits_n >= 90) || (THREADS >= 32))
        update_time = 0;

	mpz_init(tmp1);

	gettimeofday(&update_stop, NULL);
    t_update = my_difftime(&sconf->update_start, &update_stop);

	if ((num_full >= check_total) || (t_update > update_time))
	{
		//watch for an abort
		if (SIQS_ABORT)
		{
			//for fun, compute the total number of locations sieved over
			mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
			mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
			mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
			mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block	

			printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				sconf->t_time1,sconf->t_time2,sconf->t_time3);
			gmp_printf("trial division touched %d sieve locations out of %Zd\n",
				sconf->num, tmp1);
			fflush(stdout);
			fflush(stderr);

			print_factors(sconf->obj);
			exit(1);
		}

		if 	(sconf->obj->qs_obj.gbl_override_rel_flag && 
			((num_full + sconf->num_cycles) > sconf->obj->qs_obj.gbl_override_rel))
		{
			printf("\nMax specified relations found\n");
			mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
			mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
			mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
			mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block	

			printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				sconf->t_time1,sconf->t_time2,sconf->t_time3);
			gmp_printf("trial division touched %d sieve locations out of %Zd\n",
				sconf->num, tmp1);
			fflush(stdout);
			fflush(stderr);
			
			sconf->obj->qs_obj.gbl_override_rel = num_full + sconf->num_cycles;

			return 2;
		}

        t_time = my_difftime(&sconf->totaltime_start, &update_stop);
		if (sconf->obj->qs_obj.gbl_override_time_flag &&
			(t_time > sconf->obj->qs_obj.gbl_override_time))
		{
			printf("\nMax specified time limit reached\n");

			mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
			mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
			mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
			mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block				

			printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				sconf->t_time1,sconf->t_time2,sconf->t_time3);
			gmp_printf("trial division touched %d sieve locations out of %Zd\n",
				sconf->num, tmp1);
			fflush(stdout);
			fflush(stderr);
			
			return 2;
		}

		// update status on screen
		sconf->num_r = sconf->num_relations + 
		sconf->num_cycles +
		sconf->components - sconf->vertices;

		// also change rel sum to update_rels below...
		if (VFLAG >= 0)
		{
			printf("%d rels found: %d full + "
				"%d from %d partial, (%6.2f rels/sec)\r",
				sconf->num_r,sconf->num_relations,
				sconf->num_cycles +
				sconf->components - sconf->vertices,
				sconf->num_cycles,
				(double)(sconf->num_relations + sconf->num_cycles) / t_time);

			fflush(stdout);
		}

        sconf->update_start.tv_sec = update_stop.tv_sec;
        sconf->update_start.tv_usec = update_stop.tv_usec;
		sconf->t_update = 0;
		
		if (sconf->num_r >= (fb->B + sconf->num_extra_relations))
		{
			//we've got enough total relations to stop
			mpz_clear(tmp1);
			return 1;
		}
		else
		{
			// need to keep sieving.  since the last time we checked, we've found
			// (full->num_r + partial->act_r) - (last_numfull + last_numpartial)
			// relations.  assume we'll find this many next time, and 
			// scale how much longer we need to sieve to hit the target.

			sconf->num_expected = sconf->num_r - sconf->last_numfull - sconf->last_numpartial;
			// if the number expected to be found next time puts us over the needed amount, scale the 
			// check total appropriately.  otherwise, just increment check_total by check_inc
			if ((sconf->num_r + sconf->num_expected) > (fb->B + sconf->num_extra_relations))
			{
				sconf->num_needed = (fb->B + sconf->num_extra_relations) - sconf->num_r;
				sconf->check_total += 
					(uint32)((double)check_inc * (double)sconf->num_needed / (double)sconf->num_expected);
				sconf->check_total += sconf->num_extra_relations;
				sconf->update_time *= (double)sconf->num_needed / (double)sconf->num_expected;
				
				// always go at least one more second.
				if (sconf->update_time < 1)
					sconf->update_time = 1;
			}
			else
				sconf->check_total += sconf->check_inc;
		}
		sconf->last_numfull = num_full;
		sconf->last_numcycles = sconf->num_cycles;
		sconf->last_numpartial = sconf->num_cycles + sconf->components - sconf->vertices;		
	}

	mpz_clear(tmp1);
	return retcode;
}

int update_final(static_conf_t *sconf)
{
	FILE *sieve_log = sconf->obj->logfile;
	mpz_t tmp1;
	struct timeval myTVend;
    double t_time;

	mpz_init(tmp1);

	if (VFLAG >= 0)
	{
		mpz_set_ui(tmp1, sconf->tot_poly);					//total number of polys
		mpz_mul_ui(tmp1, tmp1, sconf->num_blocks);	//number of blocks
		mpz_mul_2exp(tmp1, tmp1, 1);			//pos and neg sides
		mpz_mul_2exp(tmp1, tmp1, sconf->qs_blockbits);	//sieve locations per block	

		if (VFLAG > 0)
		{
			printf("\n\nsieving required %d total polynomials (%d 'A' polynomials)\n",
				sconf->tot_poly, sconf->total_poly_a);
			gmp_printf("trial division touched %d sieve locations out of %Zd\n",
				sconf->num, tmp1);

			if (sconf->use_dlp)
				printf("squfof: %u failures, %u attempts, %u outside range, %u prp, %u useful\n", 
					sconf->failed_squfof, sconf->attempted_squfof, 
					sconf->dlp_outside_range, sconf->dlp_prp, sconf->dlp_useful);

            printf("total reports = %u, total surviving reports = %u\ntotal blocks sieved = %u,"
                "avg surviving reports per block = %1.2f\n", sconf->total_reports,
                sconf->total_surviving_reports, sconf->total_blocks,
                (float)sconf->total_surviving_reports / (float)sconf->total_blocks);

#ifdef USE_AVX2
            printf("large prime scan failures = %u\n", sconf->lp_scan_failures);
#endif

		}
		else
			printf("\n\n");

		if (sieve_log != NULL)
			logprint(sieve_log,"trial division touched %d sieve locations out of %s\n",
				sconf->num, mpz_conv2str(&gstr1.s, 10, tmp1));
		if (sconf->use_dlp)
				logprint(sieve_log, "squfof: %u failures, %u attempts, %u outside range, %u prp, %u useful\n", 
					sconf->failed_squfof, sconf->attempted_squfof, 
					sconf->dlp_outside_range, sconf->dlp_prp, sconf->dlp_useful);

#ifdef QS_TIMING

		printf("sieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
			SIEVE_STG1+SIEVE_STG2,
			TF_STG1+TF_STG2+TF_STG3+TF_STG4+TF_STG5+TF_STG6,
			POLY_STG0+POLY_STG1+POLY_STG2+POLY_STG3+POLY_STG4);

		if (sieve_log != NULL)
			logprint(sieve_log,"sieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				SIEVE_STG1+SIEVE_STG2,
				TF_STG1+TF_STG2+TF_STG3+TF_STG4+TF_STG5+TF_STG6,
				POLY_STG0+POLY_STG1+POLY_STG2+POLY_STG3+POLY_STG4);

		printf("timing for SPV check = %1.3f\n",TF_STG1);
		printf("timing for small prime trial division = %1.3f\n",TF_STG2);
		printf("timing for medium prime trial division = %1.3f\n",TF_STG3+TF_STG4);
		printf("timing for medium prime resieving test = %1.3f\n",TF_SPECIAL);
		printf("timing for large prime trial division = %1.3f\n",TF_STG5);
		printf("timing for LP splitting + buffering = %1.3f\n",TF_STG6);
		printf("timing for poly a generation = %1.3f\n",POLY_STG0);
		printf("timing for poly roots init = %1.3f\n",POLY_STG1);
		printf("timing for poly update small primes = %1.3f\n",POLY_STG2);
		printf("timing for poly sieve medium primes = %1.3f\n",POLY_STG3);
		printf("timing for poly sieve large primes = %1.3f\n",POLY_STG4);
		printf("timing for sieving small/medium primes = %1.3f\n",SIEVE_STG1);
		printf("timing for sieving large primes = %1.3f\n",SIEVE_STG2);

		if (sieve_log != NULL)
		{
			logprint(sieve_log,"timing for SPV check = %1.3f\n",TF_STG1);
			logprint(sieve_log,"timing for small prime trial division = %1.3f\n",TF_STG2);
			logprint(sieve_log,"timing for medium prime trial division = %1.3f\n",TF_STG3+TF_STG4);
			logprint(sieve_log,"timing for large prime trial division = %1.3f\n",TF_STG5);
			logprint(sieve_log,"timing for LP splitting + buffering = %1.3f\n",TF_STG6);
			logprint(sieve_log,"timing for poly a generation = %1.3f\n",POLY_STG0);
			logprint(sieve_log,"timing for poly roots init = %1.3f\n",POLY_STG1);
			logprint(sieve_log,"timing for poly update small primes = %1.3f\n",POLY_STG2);
			logprint(sieve_log,"timing for poly sieve medium primes = %1.3f\n",POLY_STG3);
			logprint(sieve_log,"timing for poly sieve large primes = %1.3f\n",POLY_STG4);
			logprint(sieve_log,"timing for sieving small/medium primes = %1.3f\n",SIEVE_STG1);
			logprint(sieve_log,"timing for sieving large primes = %1.3f\n",SIEVE_STG2);
		}
		
		
#endif
		fflush(stdout);
		fflush(stderr);
	}

	if (sieve_log != NULL)
	{
		logprint(sieve_log,"%d relations found: %d full + "
			"%d from %d partial, using %d polys (%d A polys)\n",
			sconf->num_r,sconf->num_relations,
			sconf->num_cycles +
			sconf->components - sconf->vertices,
			sconf->num_cycles,sconf->tot_poly,sconf->total_poly_a);
	}
	
	gettimeofday (&myTVend, NULL);
    t_time = my_difftime(&sconf->totaltime_start, &myTVend);

	if (sieve_log != NULL)
	{
		logprint(sieve_log,"on average, sieving found %1.2f rels/poly and %1.2f rels/sec\n",
			(double)(sconf->num_relations + sconf->num_cycles)/(double)sconf->tot_poly,
			(double)(sconf->num_relations + sconf->num_cycles) / t_time);
		logprint(sieve_log,"trial division touched %d sieve locations out of %s\n",
				sconf->num, mpz_conv2str(&gstr1.s, 10, tmp1));
		logprint(sieve_log,"==== post processing stage (msieve-1.38) ====\n");
	}

	sconf->obj->qs_obj.rels_per_sec = 
        (double)(sconf->num_relations + sconf->num_cycles) / t_time;

	if (sieve_log != NULL)	
		fflush(sieve_log);

	//sieve_log = fopen(flogname,"a");
	mpz_clear(tmp1);

	return 0;
}

int free_sieve(dynamic_conf_t *dconf)
{
	uint32 i;

	//can free sieving structures now
    dconf->sieve = dconf->sieve - 2 * 32768;
	align_free(dconf->sieve);
	align_free(dconf->comp_sieve_p->prime);
	align_free(dconf->comp_sieve_p->root1);
	align_free(dconf->comp_sieve_p->root2);
	align_free(dconf->comp_sieve_p->logp);
	align_free(dconf->comp_sieve_n->prime);
	align_free(dconf->comp_sieve_n->root1);
	align_free(dconf->comp_sieve_n->root2);
	align_free(dconf->comp_sieve_n->logp);
	free(dconf->comp_sieve_p);
	free(dconf->comp_sieve_n);

	align_free(dconf->rootupdates);
	align_free(dconf->sm_rootupdates);

	align_free(dconf->update_data.sm_firstroots1);
	align_free(dconf->update_data.sm_firstroots2);
	align_free(dconf->update_data.firstroots1);
	align_free(dconf->update_data.firstroots2);
	align_free(dconf->update_data.prime);
	align_free(dconf->update_data.logp);

	if (dconf->buckets->list != NULL)
	{
		align_free(dconf->buckets->list);
		free(dconf->buckets->fb_bounds);
		free(dconf->buckets->logp);
		align_free(dconf->buckets->num);
		free(dconf->buckets);
	}
	else
		free(dconf->buckets);

	//support data on the poly currently being sieved
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
	mpz_clear(dconf->gmptmp3);
	
	align_free(dconf->mask);
    align_free(dconf->mask2);

	//free sieve scan report stuff
	free(dconf->reports);
#ifdef USE_YAFU_TDIV
	for (i=0; i<MAX_SIEVE_REPORTS; i++)
		zFree32(&dconf->Qvals32[i]);
	free(dconf->Qvals32);
#endif
	
	for (i=0; i<MAX_SIEVE_REPORTS; i++)
	{
		mpz_clear(dconf->Qvals[i]);
	}
	free(dconf->Qvals);	
	free(dconf->valid_Qs);
	free(dconf->smooth_num);

#ifdef USE_8X_MOD_ASM
	align_free(dconf->bl_locs);
	align_free(dconf->bl_sizes);
#endif
	align_free(dconf->corrections);
    align_free(dconf->polyscratch);

#ifdef USE_VEC_SQUFOF
    free(dconf->unfactored_residue);
    free(dconf->residue_factors);
#endif

	return 0;
}

void free_filter_vars(static_conf_t *sconf)
{
	int i;

	//list of a values used first to track all a coefficients
	//generated during sieving for duplication, then used again during
	//filtering
	for (i=0; (uint32)i < sconf->total_poly_a; i++)
		mpz_clear(sconf->poly_a_list[i]);
	free(sconf->poly_a_list);

	//cycle table stuff created at the beginning of the factorization
	free(sconf->cycle_hashtable);
	free(sconf->cycle_table);

	//list of polys used in filtering and sqrt
	if (sconf->curr_b != NULL)
	{
		for (i=0;(uint32)i < sconf->bpoly_alloc;i++)
			mpz_clear(sconf->curr_b[i]);
		free(sconf->curr_b);
	}

	if (sconf->poly_list != NULL)
	{
		for (i=0;(uint32)i < sconf->poly_list_alloc;i++)
			mpz_clear(sconf->poly_list[i].b);
		free(sconf->poly_list);
	}

	if (sconf->relation_list != NULL)
	{
		//free post-processed relations
		//for (i=0; (uint32)i < sconf->num_relations; i++)
		//	free(sconf->relation_list[i].fb_offsets);
		free(sconf->relation_list);
	}

	for (i=0;i<sconf->num_cycles;i++)
	{
		if (sconf->cycle_list[i].cycle.list != NULL)
			free(sconf->cycle_list[i].cycle.list);
		if (sconf->cycle_list[i].data != NULL)
			free(sconf->cycle_list[i].data);
	}
	if (sconf->cycle_list != NULL)
		free(sconf->cycle_list);

	return;
}

int free_siqs(static_conf_t *sconf)
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
    align_free(sconf->modsqrt_array);
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
    align_free(sconf->sieve_primes);

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

	if (sconf->in_mem)
		free(sconf->in_mem_relations);

	mpz_clear(sconf->sqrt_n);
	mpz_clear(sconf->n);
	mpz_clear(sconf->target_a);

	//free(sconf->obj->qs_obj.savefile.name);
	qs_savefile_free(&sconf->obj->qs_obj.savefile);
    

	return 0;
}

uint8 choose_multiplier_siqs(uint32 B, mpz_t n) 
{
	uint32 i, j;
	uint32 num_primes = MIN(2 * B, NUM_TEST_PRIMES);
	double best_score;
	uint8 best_mult;
	double scores[NUM_MULTIPLIERS];
	uint32 num_multipliers;
	double log2n = zlog(n);

	/* measure the contribution of 2 as a factor of sieve
	   values. The multiplier itself must also be taken into
	   account in the score. scores[i] is the correction that
	   is implicitly applied to the size of sieve values for
	   multiplier i; a negative score makes sieve values 
	   smaller, and so is better */

	for (i = 0; i < NUM_MULTIPLIERS; i++) {
		uint8 curr_mult = mult_list[i];
		uint8 knmod8 = (uint8)((curr_mult * mpz_get_ui(n)) % 8);
		double logmult = log((double)curr_mult);

		/* only consider multipliers k such than
		   k*n will not overflow an mp_t */

		if (log2n + logmult > (32 * MAX_DIGITS - 2) * LN2)
			break;

		scores[i] = 0.5 * logmult;
		switch (knmod8) {
		case 1:
			scores[i] -= 2 * LN2;
			break;
		case 5:
			scores[i] -= LN2;
			break;
		case 3:
		case 7:
			scores[i] -= 0.5 * LN2;
			break;
		/* even multipliers start with a handicap */
		}
	}
	num_multipliers = i;

	/* for the rest of the small factor base primes */

	for (i = 1; i < num_primes; i++) {
		uint32 prime = (uint32)spSOEprimes[i];
		double contrib = log((double)prime) / (prime - 1);
		uint32 modp = (uint32)mpz_tdiv_ui(n, prime);

		for (j = 0; j < num_multipliers; j++) {
			uint8 curr_mult = mult_list[j];
			//uint32 knmodp = mp_modmul_1(modp, curr_mult, prime);
			uint32 knmodp = (modp * curr_mult) % prime;

			/* if prime i is actually in the factor base
			   for k * n ... */

			if (knmodp == 0 || jacobi_1(knmodp, prime) == 1) {

				/* ...add its contribution. A prime p con-
				   tributes log(p) to 1 in p sieve values, plus
				   log(p) to 1 in p^2 sieve values, etc. The
				   average contribution of all multiples of p 
				   to a random sieve value is thus

				   log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
				   = (log(p) / p) * 1 / (1 - (1/p)) 
				   = log(p) / (p-1)

				   This contribution occurs once for each
				   square root used for sieving. There are two
				   roots for each factor base prime, unless
				   the prime divides k*n. In that case there 
				   is only one root */

				if (knmodp == 0)
					scores[j] -= contrib;
				else
					scores[j] -= 2 * contrib;
			}
		}

	}

	/* use the multiplier that generates the best score */

	best_score = 1000.0;
	best_mult = 1;
	for (i = 0; i < num_multipliers; i++) {
		double score = scores[i];
		if (score < best_score) {
			best_score = score;
			best_mult = mult_list[i];
		}
	}
	return best_mult;
}

