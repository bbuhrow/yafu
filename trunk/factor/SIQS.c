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

//#define VDEBUG
// opt debug will print out some info relevant to the process used to 
// optimize the tf_small_cutoff value.  This value is varied slightly during the
// first few polynomials and the value which maximizes relation discover rate
// is chosen.
//#define OPT_DEBUG

void SIQS(fact_obj_t *fobj)
{
	//the input fobj->N and this 'n' are pointers to memory which holds
	//the input number to this function.  a copy in different memory
	//will also be made, and modified by a multiplier.  don't confuse
	//these two.
	z *n = &fobj->qs_obj.n;
	z tmp1, tmp2;

	//thread data holds all data needed during sieving
	thread_sievedata_t *thread_data;		//an array of thread data objects

	//master control structure
	static_conf_t *static_conf;

	//stuff for lanczos
	qs_la_col_t *cycle_list;
	uint32 num_cycles = 0;
	uint64 *bitfield = NULL;
	siqs_r *relation_list;

	//log file
	FILE *sieve_log;

#ifdef OPT_DEBUG
	FILE *optfile;
#endif

	//some locals
	uint32 num_needed, num_found;
	uint32 alldone;
	uint32 j;
	int i;
	clock_t start, stop;
	double t_time, opttime;
	struct timeval myTVend, optstart, optstop;
	TIME_DIFF *	difference;
	int updatecode = 0;

	//adaptive tf_small_cutoff variables
	int avg_rels_per_A = 0;
	int results[3] = {0,0,0};
	int num_avg = 10;
	int num_start = 0;
	int num_meas = 1;
	int orig_value;
	int	poly_start_num = 0;

	//logfile for this factorization
	//must ensure it is only written to by main thread
	fobj->logfile = fopen(flogname,"a");
	sieve_log = fobj->logfile;

	//first a small amount of trial division
	//which will add anything found to the global factor list
	//this is only called with the main thread
	zCopy(n,&fobj->div_obj.n);
	zTrial(10000,0,fobj);
	zCopy(&fobj->div_obj.n,n);

	//check for special cases and bail if there is one
	if ((i = check_specialcase(n,fobj->logfile,fobj)) > 0)
	{
		if (i == 1)
			fclose(sieve_log);
		return;
	}

	//At this point, we are committed to doing qs on the input
	//we need to:
	//1.) notify the screen and log
	//2.) start watching for an abort signal
	//3.) initialize data objects
	//4.) get ready to find the factor base

	//experimental
	//MAX_DIFF = 0;
	//MAX_DIFF2 = 0;
	/*for (i=0; i<20; i++)
	{
		total_primes_per_slice[i] = 0;
		count_polys_using_slice[i] = 0;
		average_primes_per_slice[i] = 0;
	}*/

	//init locals
	zInit(&tmp1);
	zInit(&tmp2);

	//fill in the factorization object
	fobj->bits = zBits(n);
	fobj->digits = ndigits(n);
	fobj->qs_obj.savefile.name = (char *)malloc(80 * sizeof(char));
	strcpy(fobj->savefile_name,siqs_savefile);

	//initialize the data objects
	static_conf = (static_conf_t *)malloc(sizeof(static_conf_t));
	static_conf->obj = fobj;

	thread_data = (thread_sievedata_t *)malloc(THREADS * sizeof(thread_sievedata_t));
	for (i=0; i<THREADS; i++)
	{
		thread_data[i].dconf = (dynamic_conf_t *)malloc(sizeof(dynamic_conf_t));
		thread_data[i].sconf = static_conf;
		//thread_data[i].tindex = i;
	}

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	SIQS_ABORT = 0;
	signal(SIGINT,siqsexit);

	//initialize offset for savefile buffer
	savefile_buf_off = 0;	
	
	//function pointer to the sieve array scanner
	scan_ptr = NULL;

	//start a counter for the whole job
	gettimeofday(&static_conf->totaltime_start, NULL);

	if (VFLAG >= 0)
		printf("\nstarting SIQS on c%d: %s\n",fobj->digits,z2decstr(n,&gstr1));
	logprint(sieve_log,"starting SIQS on c%d: %s\n",fobj->digits,z2decstr(n,&gstr1));
	logprint(sieve_log,"random seeds: %u, %u\n",g_rand.hi, g_rand.low);
	fflush(sieve_log);

	//get best parameters, multiplier, and factor base for the job
	//initialize and fill out the static part of the job data structure
	siqs_static_init(static_conf);

	//allocate structures for use in sieving with threads
	for (i=0; i<THREADS; i++)
		siqs_dynamic_init(thread_data[i].dconf, static_conf);

	//check if a savefile exists for this number, and if so load the data
	//into the master data structure
	alldone = siqs_check_restart(thread_data[0].dconf, static_conf);

	//start the process
	num_needed = static_conf->factor_base->B + static_conf->num_extra_relations;
	num_found = static_conf->num_r;
	static_conf->total_poly_a = -1;

#ifdef OPT_DEBUG
	optfile = fopen("optfile.csv","a");
	fprintf(optfile,"Optimization Debug File\n");
	fprintf(optfile,"Detected cpu %d, with L1 = %d bytes, L2 = %d bytes\n",
		yafu_get_cpu_type(),L1CACHE,L2CACHE);
	fprintf(optfile,"Starting SIQS on c%d: %s\n\n",fobj->digits,z2decstr(n,&gstr1));
	fprintf(optfile,"Meas #,Poly A #, Avg Rels/Poly/Sec, small_tf_cutoff\n");
#endif
	gettimeofday(&optstart, NULL);

	num_meas = 0;
	orig_value = static_conf->tf_small_cutoff;

	/* activate the threads one at a time. The last is the
	   master thread (i.e. not a thread at all). */

	for (i = 0; i < THREADS - 1; i++)
		start_worker_thread(thread_data + i, 0);

	start_worker_thread(thread_data + i, 1);

	//main loop
	while (num_found < num_needed)
	{
		//sequentially generate N threads worth of new poly A values,
		//using each thread's dconf.  this routine also stores the coefficients
		//in a master list (which is why this part is single threaded, for now).

		for (i=0; i<THREADS; i++)
		{
#ifdef QS_TIMING
			gettimeofday (&qs_timing_start, NULL);
#endif

			static_conf->total_poly_a++;
			//get a new 'a'
			new_poly_a(static_conf,thread_data[i].dconf);

#ifdef QS_TIMING
			gettimeofday (&qs_timing_stop, NULL);
			qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

			POLY_STG0 += 
				((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
			free(qs_timing_diff);
#endif
		}

		//process the polynomials
		for (i = 0; i < THREADS; i++) 
		{
			thread_sievedata_t *t = thread_data + i;

			if (i == THREADS - 1) {
				process_poly(t);
			}
			else {
				t->command = COMMAND_RUN;
#if defined(WIN32) || defined(_WIN64)
				SetEvent(t->run_event);
#else
				pthread_cond_signal(&t->run_cond);
				pthread_mutex_unlock(&t->run_lock);
#endif
			}
		}

		/* wait for each thread to finish */

		for (i = 0; i < THREADS; i++) {
			thread_sievedata_t *t = thread_data + i;

			if (i < THREADS - 1) {
#if defined(WIN32) || defined(_WIN64)
				WaitForSingleObject(t->finish_event, INFINITE);
#else
				pthread_mutex_lock(&t->run_lock);
				while (t->command != COMMAND_WAIT)
					pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
			}
		}

		//merge the results of each thread's work into master structure	
		for (i=0; i<THREADS; i++)
		{
			num_found = siqs_merge_data(thread_data[i].dconf,static_conf);
		}

		if (NO_SIQS_OPT == 0) 
		{			
			if (static_conf->total_poly_a > num_avg && num_meas < 3)
			{
			
				if (static_conf->total_poly_a - poly_start_num >= num_avg)
				{
					poly_start_num = static_conf->total_poly_a;

					gettimeofday (&optstop, NULL);
					difference = my_difftime (&optstart, &optstop);

					opttime = 
						((double)difference->secs + (double)difference->usecs / 1000000);
					free(difference);

					//new avg
					avg_rels_per_A = 
						(static_conf->num_relations + static_conf->num_cycles - num_start)/num_avg/opttime;

					results[num_meas] = avg_rels_per_A;
#ifdef OPT_DEBUG
					fprintf(optfile,"%d,%d,%d,%d\n",num_meas,static_conf->total_poly_a,
						avg_rels_per_A, static_conf->tf_small_cutoff);
#endif
						
					if (num_meas == 0)
					{
						static_conf->tf_small_cutoff = orig_value + 5;
					}
					else if (num_meas == 1)
					{
						static_conf->tf_small_cutoff = orig_value - 5;
					}
					else
					{
						//we've got our three measurements, make a decision.
						//the experimental results need to be several bigger than
						//our original guess in order to override it.
						//make it much harder to change if using DLP
						if (static_conf->use_dlp)
							results[0] += 16;
						else
							results[0] += 4;
						if (results[0] > results[1])
						{
							if (results[0] > results[2])
							{
								static_conf->tf_small_cutoff = orig_value;
							}
							else
							{
								static_conf->tf_small_cutoff = orig_value - 5;
							}
						}
						else
						{
							if (results[1] > results[2])
							{
								static_conf->tf_small_cutoff = orig_value + 5;
							}
							else
							{
								static_conf->tf_small_cutoff = orig_value - 5;
							}
						}	

#ifdef OPT_DEBUG
						fprintf(optfile,"final value = %d\n",static_conf->tf_small_cutoff);
#endif
					}

					num_meas++;

					num_start = static_conf->num_relations + static_conf->num_cycles;	

					gettimeofday(&optstart, NULL);

				}
			}
		}

		//free sieving structure
		for (i=0; i<THREADS; i++)
		{
			//free the relations we merged, and reset some counters/timers
			for (j=0; j<thread_data[i].dconf->buffered_rels; j++)
				free(thread_data[i].dconf->relation_buf[j].fb_offsets);
			thread_data[i].dconf->num = 0;
			thread_data[i].dconf->tot_poly = 0;
		}

		//check whether to continue or not, and update the screen
		updatecode = update_check(static_conf);
		if (updatecode)
			break;

		//if we need to continue, allocate the sieving structure again
		for (i=0; i<THREADS; i++)
		{
			//get ready to collect some more relations
			thread_data[i].dconf->buffered_rels = 0;
		}
	}

#ifdef OPT_DEBUG
	fprintf(optfile,"\n\n");
	fclose(optfile);
#endif		

	//stop worker threads
	for (i=0; i<THREADS - 1; i++)
	{
		//static_conf->tot_poly += thread_data[i].dconf->tot_poly;
		stop_worker_thread(thread_data + i, 0);
		free_sieve(thread_data[i].dconf);
		free(thread_data[i].dconf->relation_buf);
	}
	stop_worker_thread(thread_data + i, 1);
	free_sieve(thread_data[i].dconf);
	free(thread_data[i].dconf->relation_buf);			
	
	//finialize savefile
	qs_savefile_flush(&static_conf->obj->qs_obj.savefile);
	qs_savefile_close(&static_conf->obj->qs_obj.savefile);	

	//for fun, compute the total number of locations sieved over
	sp2z(static_conf->tot_poly,&tmp1);					//total number of polys
	zShortMul(&tmp1,static_conf->num_blocks,&tmp1);	//number of blocks
	zShiftLeft(&tmp1,&tmp1,1);			//pos and neg sides
	zShiftLeft(&tmp1,&tmp1,BLOCKBITS);	//sieve locations per block
	
	update_final(static_conf);

	if (updatecode == 2)
		goto done;

	//we don't need the poly_a_list anymore... free it so the other routines
	//can use it
	for (i=0;i<static_conf->total_poly_a + 1;i++)
		zFree(&static_conf->poly_a_list[i]);
	free(static_conf->poly_a_list);
	
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
	static_conf->curr_b = (z *)malloc(2 * sizeof(z));
	for (i = 0; i < 2; i++)
		zInit(&static_conf->curr_b[i]);
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
	
		yafu_find_factors(static_conf->obj, &static_conf->n, static_conf->factor_base->list, 
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

	logprint(sieve_log,"Lanczos elapsed time = %6.4f seconds.\n",static_conf->t_time1);
	logprint(sieve_log,"Sqrt elapsed time = %6.4f seconds.\n",static_conf->t_time2);
	logprint(sieve_log,"SIQS elapsed time = %6.4f seconds.\n",static_conf->t_time3);
	logprint(sieve_log,"\n");
	logprint(sieve_log,"\n");

	static_conf->cycle_list = cycle_list;
	static_conf->num_cycles = num_cycles;

	//free stuff used during filtering
	free_filter_vars(static_conf);

done:

	//free everything else
	free_siqs(thread_data[0].sconf);

	fclose(sieve_log);
	for (i=0; i<THREADS; i++)
	{
		free(thread_data[i].dconf);
	}
	free(static_conf);
	free(thread_data);
	zFree(&tmp1);
	zFree(&tmp2);
	//free(threads);

	//reset signal handler to default (no handler).
	signal(SIGINT,NULL);

	return;
}

void start_worker_thread(thread_sievedata_t *t, 
				uint32 is_master_thread) {

	/* create a thread that will process a polynomial The last poly does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread) {
		//matrix_thread_init(t);
		return;
	}

	t->command = COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
	t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
	t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	t->thread_id = CreateThread(NULL, 0, worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_create(&t->thread_id, NULL, worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
}

void stop_worker_thread(thread_sievedata_t *t,
				uint32 is_master_thread)
{
	if (is_master_thread) {
		//matrix_thread_free(t);
		return;
	}

	t->command = COMMAND_END;
#if defined(WIN32) || defined(_WIN64)
	SetEvent(t->run_event);
	WaitForSingleObject(t->thread_id, INFINITE);
	CloseHandle(t->thread_id);
	CloseHandle(t->run_event);
	CloseHandle(t->finish_event);
#else
	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_join(t->thread_id, NULL);
	pthread_cond_destroy(&t->run_cond);
	pthread_mutex_destroy(&t->run_lock);
#endif
}

#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI worker_thread_main(LPVOID thread_data) {
#else
void *worker_thread_main(void *thread_data) {
#endif
	thread_sievedata_t *t = (thread_sievedata_t *)thread_data;

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);		
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */

		if (t->command == COMMAND_RUN)
			process_poly(t);
		else if (t->command == COMMAND_END)
			break;

		/* signal completion */

		t->command = COMMAND_WAIT;
#if defined(WIN32) || defined(_WIN64)
		SetEvent(t->finish_event);		
#else
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

	//matrix_thread_free(t);

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}

void *process_poly(void *ptr)
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

	//lock_thread_to_core();

#ifdef QS_TIMING
	gettimeofday (&qs_timing_start, NULL);
#endif

	//update the gray code
	get_gray_code(dconf->curr_poly);

	//update roots, etc.
	dconf->maxB = 1<<(dconf->curr_poly->s-1);
	dconf->numB = 1;
	computeBl(sconf,dconf);
	firstRoots(sconf,dconf);

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);
	POLY_STG1 += 
		((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif

	//loop over each possible b value, for the current a value
	for ( ; dconf->numB < dconf->maxB; dconf->numB++, dconf->tot_poly++)
	{
		int skip_x2_poly = dconf->numB & 0x1;
		//setting these to be invalid means the last entry of every block will be ignored
		//so we're throwing away 1/blocksize relations, potentially, but gaining 
		//a faster sieve routine.
		uint32 invalid_root_marker = 0xFFFFFFFF; //(BLOCKSIZEm1 << 16) | BLOCKSIZEm1;

		for (i=0; i < num_blocks; i++)
		{
			//set the roots for the factors of a such that
			//they will not be sieved.  we haven't found roots for them
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p);
			lp_sieveblock(sieve,fb_sieve_p,fb->med_B, i, num_blocks, buckets,
				start_prime,blockinit,0);

			//set the roots for the factors of a to force the following routine
			//to explicitly trial divide since we haven't found roots for them
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_p);
			scan_ptr(i,0,sconf,dconf);

			//set the roots for the factors of a such that
			//they will not be sieved.  we haven't found roots for them
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n);
			lp_sieveblock(sieve,fb_sieve_n,fb->med_B, i, num_blocks, buckets,
				start_prime,blockinit,1);

			//set the roots for the factors of a to force the following routine
			//to explicitly trial divide since we haven't found roots for them
			set_aprime_roots(invalid_root_marker, poly->qlisort, poly->s, fb_sieve_n);
			scan_ptr(i,1,sconf,dconf);			

		}

		//exit(1);

		//next polynomial
		//use the stored Bl's and the gray code to find the next b
		nextB(dconf,sconf);
		//and update the roots
		nextRoots(sconf, dconf);

	}
	
	//printf("average utilizaiton of buckets in slices\n");
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

	//save the A value
	sprintf(buf,"A %s\n",z2hexstr(&dconf->curr_poly->poly_a,&gstr1));
	qs_savefile_write_line(&sconf->obj->qs_obj.savefile,buf);

	//save the data and merge into master cycle structure
	for (i=0; i<dconf->buffered_rels; i++)
	{
		rel = dconf->relation_buf + i;
		save_relation_siqs(rel->sieve_offset,rel->large_prime,
			rel->num_factors, rel->fb_offsets, rel->poly_idx, 
			rel->parity, sconf);
	}

	//update some progress indicators
	sconf->num += dconf->num;
	sconf->tot_poly += dconf->tot_poly;

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

	//we're now almost ready to start, but first
	//check if this number has had work done
	restart_siqs(sconf,dconf);
	
	if ((uint32)sconf->num_r >= sconf->factor_base->B + sconf->num_extra_relations) 
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
	}
	else
	{
		//no relations found, get ready for new factorization
		//we'll be writing to the savefile as we go, so get it ready
		qs_savefile_open(&obj->qs_obj.savefile,SAVEFILE_WRITE);
		sprintf(buf,"N %s\n",z2hexstr(&sconf->n,&gstr1));
		qs_savefile_write_line(&obj->qs_obj.savefile,buf);
		qs_savefile_flush(&obj->qs_obj.savefile);
		qs_savefile_close(&obj->qs_obj.savefile);
		//and get ready for collecting relations
		qs_savefile_open(&obj->qs_obj.savefile,SAVEFILE_APPEND);
	}

	//print some info to the screen and the log file
	if (VFLAG > 0)
	{	
		printf("\n==== sieve params ====\n");
		printf("n = %d digits, %d bits\n",sconf->digits_n,sconf->bits);
		printf("factor base: %d primes (max prime = %u)\n",sconf->factor_base->B,sconf->pmax);
		printf("single large prime cutoff: %u (%d * pmax)\n",
			sconf->large_prime_max,sconf->large_mult);
		if (sconf->use_dlp)
		{
			printf("double large prime range from %d to %d bits\n",
				sconf->dlp_lower,sconf->dlp_upper);
			printf("double large prime cutoff: %s\n",
				z2decstr(&sconf->large_prime_max2,&gstr1));
		}
		if (dconf->buckets->list != NULL)
		{
			printf("allocating %d large prime slices of factor base\n",
				dconf->buckets->alloc_slices);
			printf("buckets hold %d elements\n",BUCKET_ALLOC);
		}
		printf("sieve interval: %d blocks of size %d\n",
			sconf->num_blocks,BLOCKSIZE);
		printf("polynomial A has ~ %d factors\n",zBits(&sconf->target_a)/11);
		printf("using multiplier of %u\n",sconf->multiplier);
		printf("using small prime variation correction of %d bits\n",
			sconf->tf_small_cutoff);
//#if defined(__MINGW32__)
//	#if defined(HAS_SSE2)
//		printf("using SSE2 for trial division and x%d sieve scanning\n",
//			sconf->scan_unrolling);
//	#elif defined(HAS_MMX)
//		printf("using MMX for trial division and x%d sieve scanning\n",
//			sconf->scan_unrolling);
//	#else
//		printf("using generic trial division and x%d sieve scanning\n",
//			sconf->scan_unrolling);
//	#endif
#if defined(__GNUC__)
	#if defined(HAS_SSE2)
		printf("using SSE2 for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#elif defined(HAS_MMX)
		printf("using MMX for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#else
		printf("using generic trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#endif
	#if defined(SSE2_RESIEVEING)
		#ifdef YAFU_64K
			printf("using SSE2 for resieving 14-16 bit primes\n");
		#else
			printf("using SSE2 for resieving 13-16 bit primes\n");
		#endif
	#endif
#elif defined(_MSC_VER)
	#if defined(HAS_SSE2)
		printf("using SSE2 for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#elif defined(HAS_MMX)
		printf("using MMX for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#else
		printf("using generic trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#endif

	#if defined(HAS_SSE2)
		#ifdef YAFU_64K
			printf("using SSE2 for resieving 14-16 bit primes\n");
		#else
			printf("using SSE2 for resieving 13-16 bit primes\n");
		#endif
	#endif
#else	/* compiler not recognized*/
	
#endif
#if defined(CACHE_LINE_64) && defined(MANUAL_PREFETCH)
		printf("using 64 byte cache line prefetching\n");
#elif defined(MANUAL_PREFETCH)
		printf("using 32 byte cache line prefetching\n");
#endif
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
			printf("\n==== sieving in progress (%2d threads): %7u relations needed ====\n",
				THREADS,sconf->factor_base->B + sconf->num_extra_relations);
			printf(  "====            Press ctrl-c to abort and save state            ====\n");
		}
	}

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
		logprint(sconf->obj->logfile,"double large prime cutoff: %s\n",
			z2decstr(&sconf->large_prime_max2,&gstr1));
	}
	if (dconf->buckets->list != NULL)
	{
		logprint(sconf->obj->logfile,"allocating %d large prime slices of factor base\n",
			dconf->buckets->alloc_slices);
		logprint(sconf->obj->logfile,"buckets hold %d elements\n",BUCKET_ALLOC);
	}
	logprint(sconf->obj->logfile,"sieve interval: %d blocks of size %d\n",
		sconf->num_blocks,BLOCKSIZE);
	logprint(sconf->obj->logfile,"polynomial A has ~ %d factors\n",zBits(&sconf->target_a)/11);
	logprint(sconf->obj->logfile,"using multiplier of %u\n",sconf->multiplier);
	logprint(sconf->obj->logfile,"using small prime variation correction of %d bits\n",
		sconf->tf_small_cutoff);
//#if defined(__MINGW32__)
//	#if defined(HAS_SSE2)
//		logprint(sconf->obj->logfile,"using SSE2 for trial division and x%d sieve scanning\n",
//			sconf->scan_unrolling);
//	#elif defined(HAS_MMX)
//		logprint(sconf->obj->logfile,"using MMX for trial division and x%d sieve scanning\n",
//			sconf->scan_unrolling);
//	#else
//		logprint(sconf->obj->logfile,"using generic trial division and x%d sieve scanning\n",
//			sconf->scan_unrolling);
//	#endif
#if defined(__GNUC__)
	#if defined(HAS_SSE2)
		logprint(sconf->obj->logfile,"using SSE2 for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#elif defined(HAS_MMX)
		logprint(sconf->obj->logfile,"using MMX for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#else
		logprint(sconf->obj->logfile,"using generic trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#endif
#elif defined(_MSC_VER)
	#if defined(HAS_SSE2)
		logprint(sconf->obj->logfile,"using SSE2 for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#elif defined(HAS_MMX)
		logprint(sconf->obj->logfile,"using MMX for trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#else
		logprint(sconf->obj->logfile,"using generic trial division and x%d sieve scanning\n",
			sconf->scan_unrolling);
	#endif
#else	/* compiler not recognized*/
	
#endif	
#if defined(CACHE_LINE_64) && defined(MANUAL_PREFETCH)
	logprint(sconf->obj->logfile,"using 64 byte cache line prefetching\n");
#elif defined(MANUAL_PREFETCH)
	logprint(sconf->obj->logfile,"using 32 byte cache line prefetching\n");
#endif
	logprint(sconf->obj->logfile,"trial factoring cutoff at %d bits\n",
		sconf->tf_closnuf);
	if (THREADS == 1)
		logprint(sconf->obj->logfile,"==== sieving started (1 thread) ====\n");
	else
		logprint(sconf->obj->logfile,"==== sieving started (%2d threads) ====\n",THREADS);


	return state;
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
	zInit(&dconf->qstmp1);
	zInit(&dconf->qstmp2);
	zInit(&dconf->qstmp3);
	zInit(&dconf->qstmp4);
	zInit32(&dconf->qstmp32);

	//this stuff changes with every new poly
	//allocate a polynomial structure which will hold the current
	//set of polynomial coefficients (a,b,c) and other info
	dconf->curr_poly = (siqs_poly *)malloc(sizeof(siqs_poly));
	zInit(&dconf->curr_poly->poly_a);
	zInit(&dconf->curr_poly->poly_b);
	zInit(&dconf->curr_poly->poly_c);
	dconf->curr_poly->qlisort = (int *)malloc(MAX_A_FACTORS*sizeof(int));
	dconf->curr_poly->gray = (char *) malloc( 65536 * sizeof(char));
	dconf->curr_poly->nu = (char *) malloc( 65536 * sizeof(char));
	if (VFLAG > 2)
	{
		memsize = dconf->curr_poly->poly_a.alloc * sizeof(fp_digit) * 3;
		memsize += MAX_A_FACTORS*sizeof(int);
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
#if defined (_MSC_VER) || defined(__MINGW32__)

#ifdef USE_COMPRESSED_FB
	dconf->comp_sieve_p = (sieve_fb_compressed *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(sieve_fb_compressed)),64);
	dconf->comp_sieve_n = (sieve_fb_compressed *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(sieve_fb_compressed)),64);
#else
	dconf->comp_sieve_p = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));
	dconf->comp_sieve_n = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));
	dconf->comp_sieve_p->prime = (uint16 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)),64);
	dconf->comp_sieve_p->root1 = (uint16 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)),64);
	dconf->comp_sieve_p->root2 = (uint16 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)),64);
	dconf->comp_sieve_p->logp = (uint8 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint8)),64);

	dconf->comp_sieve_n->prime = (uint16 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)),64);
	dconf->comp_sieve_n->root1 = (uint16 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)),64);
	dconf->comp_sieve_n->root2 = (uint16 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)),64);
	dconf->comp_sieve_n->logp = (uint8 *)_aligned_malloc(
		(size_t)(sconf->factor_base->med_B * sizeof(uint8)),64);
#endif
	dconf->fb_sieve_p = (sieve_fb *)_aligned_malloc(
		(size_t)(sconf->factor_base->B * sizeof(sieve_fb)),64);
	dconf->fb_sieve_n = (sieve_fb *)_aligned_malloc(
		(size_t)(sconf->factor_base->B * sizeof(sieve_fb)),64);
	//allocate storage for the update data needed when changing polys
	dconf->update_data.firstroots1 = (int *)_aligned_malloc(
		sconf->factor_base->B * sizeof(int),64);
	dconf->update_data.firstroots2 = (int *)_aligned_malloc(
		sconf->factor_base->B * sizeof(int),64);
	dconf->update_data.prime = (uint32 *)_aligned_malloc(
		sconf->factor_base->B * sizeof(uint32),64);
	dconf->update_data.logp = (uint8 *)_aligned_malloc(
		sconf->factor_base->B * sizeof(uint8),64);
	dconf->rootupdates = (int *)_aligned_malloc(
		MAX_A_FACTORS * sconf->factor_base->B * sizeof(int),64);

#else

#ifdef USE_COMPRESSED_FB
	dconf->comp_sieve_p = (sieve_fb_compressed *)memalign(64,
			(size_t)(sconf->factor_base->med_B * sizeof(sieve_fb_compressed)));
	dconf->comp_sieve_n = (sieve_fb_compressed *)memalign(64,
			(size_t)(sconf->factor_base->med_B * sizeof(sieve_fb_compressed)));
#else
	dconf->comp_sieve_p = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));
	dconf->comp_sieve_n = (sieve_fb_compressed *)malloc(sizeof(sieve_fb_compressed));
	dconf->comp_sieve_p->prime = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_p->root1 = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_p->root2 = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
#ifdef LOGP_BITS
	//printf("logp values have 8 bits\n");
	dconf->comp_sieve_p->logp = (uint8 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint8)));
#else
	//printf("logp values have 16 bits\n");
	dconf->comp_sieve_p->logp = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
#endif

	dconf->comp_sieve_n->prime = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_n->root1 = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
	dconf->comp_sieve_n->root2 = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
#ifdef LOGP_BITS
	//printf("logp values have 8 bits\n");
	dconf->comp_sieve_n->logp = (uint8 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint8)));
#else
	//printf("logp values have 16 bits\n");
	dconf->comp_sieve_n->logp = (uint16 *)memalign(64,
		(size_t)(sconf->factor_base->med_B * sizeof(uint16)));
#endif

#endif
	dconf->fb_sieve_p = (sieve_fb *)memalign(64,
			(size_t)(sconf->factor_base->B * sizeof(sieve_fb)));
	dconf->fb_sieve_n = (sieve_fb *)memalign(64,
			(size_t)(sconf->factor_base->B * sizeof(sieve_fb)));
	//allocate storage for the update data needed when changing polys
	dconf->update_data.firstroots1 = (int *)memalign(64,
		sconf->factor_base->B * sizeof(int));
	dconf->update_data.firstroots2 = (int *)memalign(64,
		sconf->factor_base->B * sizeof(int));
	dconf->update_data.prime = (uint32 *)memalign(64,
		sconf->factor_base->B * sizeof(uint32));
	dconf->update_data.logp = (uint8 *)memalign(64,
		sconf->factor_base->B * sizeof(uint8));
	dconf->rootupdates = (int *)memalign(64,
		MAX_A_FACTORS * sconf->factor_base->B * sizeof(int));

#endif

	if (VFLAG > 2)
	{
		memsize = sconf->factor_base->med_B * sizeof(sieve_fb_compressed) * 2;
		memsize += sconf->factor_base->B * sizeof(sieve_fb) * 2;
		printf("\tfactor bases: %d bytes\n",memsize);
	}

	if (VFLAG > 2)
	{
		memsize = sconf->factor_base->B * sizeof(int)  * 2;
		memsize += sconf->factor_base->B * sizeof(uint32);
		memsize += sconf->factor_base->B * sizeof(uint8);
		printf("\tupdate data: %d bytes\n",memsize);
	}
	
	//allocate the sieve
#if defined (_MSC_VER) || defined (__MINGW32__)
	dconf->sieve = (uint8 *)_aligned_malloc(
		(size_t) (BLOCKSIZE * sizeof(uint8)),64);
#else
	dconf->sieve = (uint8 *)memalign(64,
		(size_t) (BLOCKSIZE * sizeof(uint8)));
#endif

	if (VFLAG > 2)
	{
		memsize = BLOCKSIZE * sizeof(uint8);
		printf("\tsieve: %d bytes\n",memsize);
	}

	//allocate the Bl array, space for MAX_Bl bigint numbers
	dconf->Bl = (z *)malloc(MAX_A_FACTORS * sizeof(z));
	for (i=0;i<MAX_A_FACTORS;i++)
		zInit(&dconf->Bl[i]);

	if (VFLAG > 2)
	{
		memsize = MAX_A_FACTORS * dconf->Bl[0].alloc * sizeof(fp_digit);
		printf("\tBl data: %d bytes\n",memsize);
	}

	//copy the unchanging part to the sieving factor bases
	for (i = 2; i < sconf->factor_base->med_B; i++)
	{
		uint32 p = sconf->factor_base->list->prime[i];
		uint32 lp = sconf->factor_base->list->logprime[i];

#ifdef USE_COMPRESSED_FB
		dconf->comp_sieve_p[i].prime_and_logp = lp << 16;
		dconf->comp_sieve_p[i].prime_and_logp |= (uint16)p;
		dconf->comp_sieve_n[i].prime_and_logp = lp << 16;
		dconf->comp_sieve_n[i].prime_and_logp |= (uint16)p;
#else
		dconf->comp_sieve_p->logp[i] = (uint8)lp;
		dconf->comp_sieve_p->prime[i] = (uint16)p;
		dconf->comp_sieve_n->logp[i] = (uint8)lp;
		dconf->comp_sieve_n->prime[i] = (uint16)p;
#endif

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

	//check if we should use bucket sieving, and allocate structures if so
	if (sconf->factor_base->B > sconf->factor_base->med_B)
	{
		dconf->buckets = (lp_bucket *)malloc(sizeof(lp_bucket));

		//test to see how many slices we'll need.
		testfirstRoots(sconf,dconf);

		//initialize the bucket lists and auxilary info.
		dconf->buckets->num = (uint32 *)malloc(
			2 * sconf->num_blocks * dconf->buckets->alloc_slices * sizeof(uint32));
		dconf->buckets->fb_bounds = (uint32 *)malloc(
			dconf->buckets->alloc_slices * sizeof(uint32));
		dconf->buckets->logp = (uint8 *)malloc(
			dconf->buckets->alloc_slices * sizeof(uint8));
		dconf->buckets->list_size = 2 * sconf->num_blocks * dconf->buckets->alloc_slices;
		
		//now allocate the buckets
#if defined(_MSC_VER) || defined (__MINGW32__)
		dconf->buckets->list = (uint32 *)_aligned_malloc(
			2 * sconf->num_blocks * dconf->buckets->alloc_slices * 
			BUCKET_ALLOC * sizeof(uint32),64);
#else
		dconf->buckets->list = (uint32 *)memalign(64,
			2 * sconf->num_blocks * dconf->buckets->alloc_slices * 
			BUCKET_ALLOC * sizeof(uint32));
#endif

	}
	else
	{
		dconf->buckets = (lp_bucket *)malloc(sizeof(lp_bucket));
		dconf->buckets->list = NULL;
		dconf->buckets->alloc_slices = 0;
		dconf->buckets->num_slices = 0;
	}

	if (VFLAG > 2)
	{
		memsize = 2 * sconf->num_blocks * dconf->buckets->alloc_slices * sizeof(uint32);
		memsize += dconf->buckets->alloc_slices * sizeof(uint32);
		memsize += dconf->buckets->alloc_slices * sizeof(uint8);
		memsize += 2 * sconf->num_blocks * dconf->buckets->alloc_slices * BUCKET_ALLOC * sizeof(uint32);
		printf("\tbucket data: %d bytes\n",memsize);
	}

	//used in trial division to mask out the fb_index portion of bucket entries, so that
	//multiple block locations can be searched for in parallel using SSE2 instructions
#if defined (_MSC_VER) || defined (__MINGW32__)
	dconf->mask = (uint16 *)_aligned_malloc(8 * sizeof(uint16),16);
#elif defined (__GNUC__)
	dconf->mask = (uint16 *)memalign(16,8 * sizeof(uint16));
//#else
//	dconf->mask = (uint16 *)malloc(8 * sizeof(uint16));
#endif
	dconf->mask[1] = 0xFFFF;
	dconf->mask[3] = 0xFFFF;
	dconf->mask[5] = 0xFFFF;
	dconf->mask[7] = 0xFFFF;

	// array of sieve locations scanned from the sieve block that we
	// will submit to trial division.  make it the size of a sieve block 
	// in the pathological case that every sieve location is a report
	dconf->reports = (uint32 *)malloc(MAX_SIEVE_REPORTS * sizeof(uint32));
	dconf->num_reports = 0;
	dconf->Qvals = (z32 *)malloc(MAX_SIEVE_REPORTS * sizeof(z32));
	for (i=0; i<MAX_SIEVE_REPORTS; i++)
		zInit32(&dconf->Qvals[i]);
	dconf->valid_Qs = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));
	dconf->smooth_num = (int *)malloc(MAX_SIEVE_REPORTS * sizeof(int));

	//initialize some counters
	dconf->tot_poly = 0;		//track total number of polys
	dconf->num = 0;				//sieve locations subjected to trial division

	return 0;
}

int siqs_static_init(static_conf_t *sconf)
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
	z tmp1, tmp2;

	if (VFLAG > 2)
	{
		printf("static memory usage:\n");
	}

	//local bigints
	zInit(&tmp1);
	zInit(&tmp2);

	//default parameters
	sconf->fudge_factor = 1.3;
	sconf->large_mult = 30;
	sconf->num_blocks = 40;
	sconf->num_extra_relations = 64;
	sconf->small_limit = 256;
	sconf->use_dlp = 0;

	//allocate the space for the factor base structure
	sconf->factor_base = (fb_list *)malloc(sizeof(fb_list));

	//allocate space for a copy of input number in the job structure
	zInit(&sconf->n);
	zCopy(&obj->qs_obj.n,&sconf->n);

	//initialize some constants
	zInit(&sconf->sqrt_n);
	zInit(&sconf->target_a);
	zInit(&sconf->max_fb2);
	zInit(&sconf->large_prime_max2);

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

		//if (sconf->factor_base->B & 0x1)
		//	sconf->factor_base->B++;	//force this to be even for fast root updating

		//allocate the space for the factor base elements
#if defined(_MSC_VER) || defined (__MINGW32__)
		sconf->factor_base->list = (fb_element_siqs *)_aligned_malloc(
			(size_t)(sizeof(fb_element_siqs)),64);
		sconf->modsqrt_array = (uint32 *)malloc(
			sconf->factor_base->B * sizeof(uint32));
		sconf->factor_base->list->prime = (uint32 *)_aligned_malloc(
			(size_t)(sconf->factor_base->B * sizeof(uint32)),64);
		sconf->factor_base->list->small_inv = (uint32 *)_aligned_malloc(
			(size_t)(sconf->factor_base->B * sizeof(uint32)),64);
		sconf->factor_base->list->correction = (uint32 *)_aligned_malloc(
			(size_t)(sconf->factor_base->B * sizeof(uint32)),64);
		sconf->factor_base->list->logprime = (uint32 *)_aligned_malloc(
			(size_t)(sconf->factor_base->B * sizeof(uint32)),64);
//#elif defined(__MINGW32__)
//		sconf->factor_base->list = (fb_element_siqs *)malloc(
//			(size_t)(sizeof(fb_element_siqs)));
//		sconf->modsqrt_array = (uint32 *)malloc(
//			sconf->factor_base->B * sizeof(uint32));
//		sconf->factor_base->list->prime = (uint32 *)malloc(
//			(size_t)(sconf->factor_base->B * sizeof(uint32)));
//		sconf->factor_base->list->small_inv = (uint32 *)malloc(
//			(size_t)(sconf->factor_base->B * sizeof(uint32)));
//		sconf->factor_base->list->correction = (uint32 *)malloc(
//			(size_t)(sconf->factor_base->B * sizeof(uint32)));
//		sconf->factor_base->list->logprime = (uint32 *)malloc(
//			(size_t)(sconf->factor_base->B * sizeof(uint32)));
#else
		sconf->factor_base->list = (fb_element_siqs *)memalign(64,
			(size_t)(sizeof(fb_element_siqs)));
		sconf->modsqrt_array = (uint32 *)malloc(
			sconf->factor_base->B * sizeof(uint32));
		sconf->factor_base->list->prime = (uint32 *)memalign(64,
			(size_t)(sconf->factor_base->B * sizeof(uint32)));
		sconf->factor_base->list->small_inv = (uint32 *)memalign(64,
			(size_t)(sconf->factor_base->B * sizeof(uint32)));
		sconf->factor_base->list->correction = (uint32 *)memalign(64,
			(size_t)(sconf->factor_base->B * sizeof(uint32)));
		sconf->factor_base->list->logprime = (uint32 *)memalign(64,
			(size_t)(sconf->factor_base->B * sizeof(uint32)));
#endif

		if (VFLAG > 2)
		{
			memsize = sconf->factor_base->B * sizeof(uint32) * 5;
			printf("\tfactor base: %d bytes\n",memsize);
		}

		//find multiplier
		sconf->multiplier = (uint32)choose_multiplier_siqs(sconf->factor_base->B, &sconf->n);

		zShortMul(&sconf->n,sconf->multiplier,&tmp1);
		zCopy(&tmp1,&sconf->n);

		//sconf holds n*mul, so update its digit count and number of bits
		sconf->digits_n = ndigits(&sconf->n);
		sconf->bits = zBits(&sconf->n);

		//find sqrt_n
		zNroot(&sconf->n,&sconf->sqrt_n,2);

		//construct the factor base - divide out any primes which 
		//factor n other than mul
		//and go back to get_params if any are found
		if ((i = make_fb_siqs(sconf)) == 0)
			break;
		else
		{
			//i is already divided out of n.  record the factor we found
			sp2z(i,&tmp1);
			tmp1.type = PRIME;
			add_to_factor_list(sconf->obj, &tmp1);
			tmp1.type = UNKNOWN;
			//and remove the multiplier we may have added, so that
			//we can try again to build a factor base.
			zShortDiv(&sconf->n,sconf->multiplier,&sconf->n);
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

	//adjust for various architectures
	if (BLOCKSIZE == 32768)
		sconf->num_blocks *= 2;
	if (BLOCKSIZE == 16384)
		sconf->num_blocks *= 4;

	if (sconf->num_blocks < 1)
		sconf->num_blocks = 1;

	//and adjust for really small jobs
	if (sconf->bits <= 140)
		sconf->num_blocks = 1;

	//set the sieve interval.  this many blocks on each side of 0
	sconf->sieve_interval = BLOCKSIZE * sconf->num_blocks;

	//compute sieving limits
	sconf->factor_base->small_B = MIN(
		sconf->factor_base->B,1024); //((INNER_BLOCKSIZE)/(sizeof(sieve_fb))));

	for (i = 2; i < sconf->factor_base->B; i++)
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
			(i % 8 == 0)) break;
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
		if ((sconf->factor_base->list->prime[i] > (uint32)(1.5 * (double)BLOCKSIZE))  &&
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
		printf("fb bounds\n\tsmall: %u\n\t13bit: %u\n\t14bit: %u\n\t15bit: %u\n\tmed: %u\n\tlarge: %u\n\tall: %u\n",
			sconf->factor_base->small_B,
			sconf->factor_base->fb_13bit_B,
			sconf->factor_base->fb_14bit_B,
			sconf->factor_base->fb_15bit_B,
			sconf->factor_base->med_B,
			sconf->factor_base->large_B,
			sconf->factor_base->B);

		printf("start primes\n\t13bit: %u\n\t14bit: %u\n\t15bit: %u\n\tmed: %u\n\tlarge: %u\n",
			sconf->factor_base->list->prime[sconf->factor_base->fb_13bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_14bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->fb_15bit_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->med_B-1],
			sconf->factor_base->list->prime[sconf->factor_base->large_B-1]);
	}
	//exit(1);

	//a couple limits
	sconf->pmax = sconf->factor_base->list->prime[sconf->factor_base->B-1];
	sconf->large_prime_max = sconf->pmax * sconf->large_mult;

	//based on the size of the input, determine how to proceed.
	if (sconf->digits_n > 81 || gbl_force_DLP)
	{
		sconf->use_dlp = 1;
		scan_ptr = &check_relations_siqs_16;
		sconf->scan_unrolling = 128;
	}
	else
	{
		if (sconf->digits_n < 60)
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
	qs_savefile_init(&obj->qs_obj.savefile, siqs_savefile);

	//if we're using dlp, compute the range of residues which will
	//be subjected to factorization beyond trial division
	if (sconf->use_dlp)
	{
		sp2z(sconf->factor_base->list->prime[sconf->factor_base->B - 1],&tmp1);
		zSqr(&tmp1,&tmp1);
		zCopy(&tmp1,&sconf->max_fb2);
		sconf->dlp_lower = zBits(&tmp1) + 1;
		sum = pow((double)sconf->large_prime_max,1.8);
		dbl2z(sum,&tmp1);
		zCopy(&tmp1,&sconf->large_prime_max2);
		sconf->dlp_upper = zBits(&tmp1);
	}

	//'a' values should be as close as possible to sqrt(2n)/M, 
	zShiftLeft(&tmp1,&sconf->n,1);
	zNroot(&tmp1,&tmp2,2);
	zShortDiv(&tmp2,sconf->sieve_interval,&sconf->target_a);

	//compute the number of bits in M/2*sqrt(N/2), the approximate value
	//of residues in the sieve interval.  Then subtract some slack.
	//sieve locations greater than this are worthy of trial dividing
	closnuf = (uint8)(double)((zBits(&sconf->n) - 1)/2);
	closnuf += (uint8)(log((double)sconf->sieve_interval/2)/log(2.0));
	closnuf -= (uint8)(sconf->fudge_factor * log(sconf->large_prime_max) / log(2.0));
	
	//it pays to trial divide a little more for big jobs...
	//"optimized" is used rather loosely here, but these corrections
	//were observed to make things faster.  
	if (sconf->digits_n >=70 && sconf->digits_n < 75)
		closnuf += 1;	//optimized for 75 digit num
	else if (sconf->digits_n >=75 && sconf->digits_n < 78)
		closnuf -= 0;	//optimized for 75 digit num
	else if (sconf->digits_n >=78 && sconf->digits_n < 81)
		closnuf -= 2;	//optimized for 80 digit num
	else if (sconf->digits_n >=81 && sconf->digits_n < 84)
		closnuf -= 3;	//optimized for 82 digit num
	else if (sconf->digits_n >=84 && sconf->digits_n < 87)
		closnuf -= 4;	//optimized for 85 digit num
	else if (sconf->digits_n >=87 && sconf->digits_n < 90)
		closnuf -= 4;
	else if (sconf->digits_n >= 90)
	{
		//this is guesswork, but we definately want to expend more effort
		//to find relations for a given polynomial.  
		closnuf -= 5;
		for (i = sconf->digits_n; i > 90; i--)
			closnuf--;
	}

	//now add in any extra for double large prime variation, if active.
	//also dependant on input size... haven't figured out a good relationship
	//yet, so these are "handtuned" for now.
	if (sconf->use_dlp)
	{
		if (sconf->digits_n >=82 && sconf->digits_n < 85)
			closnuf -= 4;	
		else if (sconf->digits_n >=85 && sconf->digits_n < 90)
			closnuf -= 7;	
		else if (sconf->digits_n >=90 && sconf->digits_n < 95)
			closnuf -= 10;	
		else if (sconf->digits_n >=95 && sconf->digits_n < 100)
			closnuf -= 13;	
		else if (sconf->digits_n >=100)
			closnuf -= 17;	
	}

	//test the contribution of the small primes to the sieve.  
	for (i = 2; i < sconf->factor_base->B; i++)
	{
		if (sconf->factor_base->list->prime[i] > sconf->small_limit)
			break;
	}
	sconf->sieve_small_fb_start = i;

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
			(double)BLOCKSIZE / (double)prime;
	}
	avg = 2*sum/BLOCKSIZE;
	//this was observed to be the typical standard dev. for one block of 
	//test sieving... wrap this magic number along with several others into
	//one empirically determined fudge factor...
	sd = sqrt(28);

	//this appears to work fairly well... paper mentioned doing it this
	//way... find out and reference here.
	sconf->tf_small_cutoff = (uint8)(avg + 2.5*sd);
	closnuf -= sconf->tf_small_cutoff;

	if (gbl_override_tf_flag)
	{
		closnuf = gbl_override_tf;
		printf("overriding with new closnuf = %d\n",closnuf);
	}

	sconf->blockinit = closnuf;
	sconf->tf_closnuf = closnuf;

	//needed during filtering
	zInit(&sconf->curr_a);
	sconf->curr_poly = (siqs_poly *)malloc(sizeof(siqs_poly));
	zInit(&sconf->curr_poly->poly_a);
	zInit(&sconf->curr_poly->poly_b);
	zInit(&sconf->curr_poly->poly_c);
	sconf->curr_poly->qlisort = (int *)malloc(MAX_A_FACTORS*sizeof(int));
	sconf->curr_poly->gray = (char *) malloc( 65536 * sizeof(char));
	sconf->curr_poly->nu = (char *) malloc( 65536 * sizeof(char));

	//initialize a list of all poly_a values used 
	sconf->poly_a_list = (z *)malloc(sizeof(z));

	//compute how often to check our list of partial relations and update the gui.
	sconf->check_inc = sconf->factor_base->B/10;
	sconf->check_total = sconf->check_inc;
	sconf->update_time = 5;

	//get ready for the factorization and screen updating
	sconf->t_update=0;
	sconf->last_numfull = 0;
	sconf->last_numpartial = 0;
	sconf->last_numcycles = 0;
	gettimeofday(&sconf->update_start, NULL);
	//sconf->update_start = clock();

	sconf->total_poly_a = 0;	//track number of A polys used
	sconf->num_r = 0;			//total relations found
	sconf->charcount = 0;		//characters on the screen

	sconf->t_time1 = 0;			//sieve time
	sconf->t_time2 = 0;			//relation scanning and trial division
	sconf->t_time3 = 0;			//polynomial calculations and large prime sieving
	sconf->t_time4 = 0;			//extra?
	
	sconf->tot_poly = 0;		//track total number of polys
	sconf->num = 0;				//sieve locations subjected to trial division

	//no factors so far...
	sconf->factor_list.num_factors = 0;

	zFree(&tmp1);
	zFree(&tmp2);
	return 0;
}

int update_check(static_conf_t *sconf)
{
	//check to see if we should combine partial relations
	//found and check to see if we're done.  also update the screen.
	//this happens one of two ways, either we have found more than a 
	//certain amount of relations since the last time we checked, or if
	//a certain amount of time has elapsed.
	z tmp1, *qstmp1;
	struct timeval update_stop;
	TIME_DIFF *	difference;
	uint32 num_full = sconf->num_relations;
	uint32 check_total = sconf->check_total;
	uint32 check_inc = sconf->check_inc;
	double update_time = sconf->update_time;
	double t_update;
	int tot_poly = sconf->tot_poly, i;
	fb_list *fb = sconf->factor_base;
	int retcode = 0;

	zInit(&tmp1);
	qstmp1 = &tmp1;

	gettimeofday(&update_stop, NULL);
	difference = my_difftime (&sconf->update_start, &update_stop);
	t_update = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (num_full >= check_total || t_update > update_time)
	{
		//watch for an abort
		if (SIQS_ABORT)
		{
			sp2z(tot_poly,qstmp1);									
			zShortMul(qstmp1,sconf->num_blocks,qstmp1);		
			zShiftLeft(qstmp1,qstmp1,1);							
			zShiftLeft(qstmp1,qstmp1,BLOCKBITS);				

			printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				sconf->t_time1,sconf->t_time2,sconf->t_time3);
			printf("trial division touched %d sieve locations out of %s\n",
				sconf->num,z2decstr(qstmp1,&gstr1));
			fflush(stdout);
			fflush(stderr);

			print_factors(sconf->obj);
			exit(1);
		}

		if 	(gbl_override_rel_flag && 
			((num_full + sconf->num_cycles) > gbl_override_rel))
		{
			printf("\nMax specified relations found\n");
			sp2z(tot_poly,qstmp1);									
			zShortMul(qstmp1,sconf->num_blocks,qstmp1);		
			zShiftLeft(qstmp1,qstmp1,1);							
			zShiftLeft(qstmp1,qstmp1,BLOCKBITS);	

			printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				sconf->t_time1,sconf->t_time2,sconf->t_time3);
			printf("trial division touched %d sieve locations out of %s\n",
				sconf->num,z2decstr(qstmp1,&gstr1));
			fflush(stdout);
			fflush(stderr);
			
			gbl_override_rel = num_full + sconf->num_cycles;

			return 2;
		}

		difference = my_difftime (&sconf->totaltime_start, &update_stop);
		if (gbl_override_time_flag &&
			(((double)difference->secs + (double)difference->usecs / 1000000) > 
			gbl_override_time))
		{
			printf("\nMax specified time limit reached\n");

			sp2z(tot_poly,qstmp1);									
			zShortMul(qstmp1,sconf->num_blocks,qstmp1);		
			zShiftLeft(qstmp1,qstmp1,1);							
			zShiftLeft(qstmp1,qstmp1,BLOCKBITS);				

			printf("\nsieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
				sconf->t_time1,sconf->t_time2,sconf->t_time3);
			printf("trial division touched %d sieve locations out of %s\n",
				sconf->num,z2decstr(qstmp1,&gstr1));
			fflush(stdout);
			fflush(stderr);
			
			return 2;
		}
		free(difference);

		//update status on screen
		sconf->num_r = sconf->num_relations + 
		sconf->num_cycles +
		sconf->components - sconf->vertices;

		//difference = my_difftime (&sconf->update_start, &update_stop);
		//also change rel sum to update_rels below...
		difference = my_difftime (&sconf->totaltime_start, &update_stop);
		if (VFLAG >= 0)
		{
			uint32 update_rels;
			for (i=0; i < sconf->charcount; i++)
				printf("\b");
			//in order to keep rels/sec from going mad when relations
			//are reloaded on a restart, just use the number of
			//relations we've found since the last update.  don't forget
			//to initialize last_numfull and partial when loading
			update_rels = sconf->num_relations + sconf->num_cycles - 
				sconf->last_numfull - sconf->last_numcycles;
			sconf->charcount = printf("%d rels found: %d full + "
				"%d from %d partial, (%6.2f rels/sec)",
				sconf->num_r,sconf->num_relations,
				sconf->num_cycles +
				sconf->components - sconf->vertices,
				sconf->num_cycles,
				(double)(sconf->num_relations + sconf->num_cycles) /
				((double)difference->secs + (double)difference->usecs / 1000000));
				//(double)(num_full + sconf->num_cycles) /
				//((double)(clock() - sconf->totaltime_start)/(double)CLOCKS_PER_SEC));
			fflush(stdout);
		}
		free(difference);

		gettimeofday(&sconf->update_start, NULL);
		sconf->t_update = 0;
		
		if (sconf->num_r >= fb->B + sconf->num_extra_relations) 
		{
			//we've got enough total relations to stop
			zFree(&tmp1);
			return 1;
		}
		else
		{
			//need to keep sieving.  since the last time we checked, we've found
			//(full->num_r + partial->act_r) - (last_numfull + last_numpartial)
			//relations.  assume we'll find this many next time, and 
			//scale how much longer we need to sieve to hit the target.

			sconf->num_expected = sconf->num_r - sconf->last_numfull - sconf->last_numpartial;
			//if the number expected to be found next time puts us over the needed amount, scale the 
			//check total appropriately.  otherwise, just increment check_total by check_inc
			if ((sconf->num_r + sconf->num_expected) > (fb->B + sconf->num_extra_relations))
			{
				sconf->num_needed = (fb->B + sconf->num_extra_relations) - sconf->num_r;
				sconf->check_total += 
					(uint32)((double)check_inc * (double)sconf->num_needed / (double)sconf->num_expected);
				sconf->check_total += sconf->num_extra_relations;
				sconf->update_time *= (double)sconf->num_needed / (double)sconf->num_expected;
				
				//always go at least one more second.
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

	zFree(&tmp1);
	return retcode;
}

int update_final(static_conf_t *sconf)
{
	FILE *sieve_log = sconf->obj->logfile;
	z qstmp1;
	struct timeval myTVend;
	TIME_DIFF *	difference;

	zInit(&qstmp1);

	if (VFLAG >= 0)
	{
		sp2z(sconf->tot_poly,&qstmp1);									
		zShortMul(&qstmp1,sconf->num_blocks,&qstmp1);		
		zShiftLeft(&qstmp1,&qstmp1,1);							
		zShiftLeft(&qstmp1,&qstmp1,BLOCKBITS);

		if (VFLAG > 0)
			printf("\n\nsieving required %d total polynomials\ntrial division touched %d sieve locations out of %s\n",
				sconf->tot_poly, sconf->num,z2decstr(&qstmp1,&gstr1));
		else
			printf("\n\n");

		logprint(sieve_log,"trial division touched %d sieve locations out of %s\n",
			sconf->num,z2decstr(&qstmp1,&gstr1));

#ifdef QS_TIMING

		printf("sieve time = %6.4f, relation time = %6.4f, poly_time = %6.4f\n",
			SIEVE_STG1+SIEVE_STG2,
			TF_STG1+TF_STG2+TF_STG3+TF_STG4+TF_STG5+TF_STG6,
			POLY_STG0+POLY_STG1+POLY_STG2+POLY_STG3+POLY_STG4);
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
		
		
#endif
		fflush(stdout);
		fflush(stderr);
	}

	logprint(sieve_log,"%d relations found: %d full + "
		"%d from %d partial, using %d polys (%d A polys)\n",
		sconf->num_r,sconf->num_relations,
		sconf->num_cycles +
		sconf->components - sconf->vertices,
		sconf->num_cycles,sconf->tot_poly,sconf->total_poly_a);
	
	gettimeofday (&myTVend, NULL);
	difference = my_difftime (&sconf->totaltime_start, &myTVend);

	logprint(sieve_log,"on average, sieving found %1.2f rels/poly and %1.2f rels/sec\n",
		(double)(sconf->num_relations + sconf->num_cycles)/(double)sconf->tot_poly,
		(double)(sconf->num_relations + sconf->num_cycles) /
		((double)difference->secs + (double)difference->usecs / 1000000));
	logprint(sieve_log,"trial division touched %d sieve locations out of %s\n",
			sconf->num,z2decstr(&qstmp1,&gstr1));
	logprint(sieve_log,"==== post processing stage (msieve-1.38) ====\n");

	sconf->obj->qs_obj.rels_per_sec = (double)(sconf->num_relations + sconf->num_cycles) /
		((double)difference->secs + (double)difference->usecs / 1000000);

	fflush(sieve_log);
	//sieve_log = fopen(flogname,"a");
	zFree(&qstmp1);
	free(difference);

	return 0;
}

int free_sieve(dynamic_conf_t *dconf)
{
	uint32 i;

	//can free sieving structures now
	align_free(dconf->sieve);
	align_free(dconf->fb_sieve_p);
	align_free(dconf->fb_sieve_n);
#ifdef USE_COMPRESSED_FB
	align_free(dconf->comp_sieve_p);
	align_free(dconf->comp_sieve_n);
#else
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
#endif
	align_free(dconf->rootupdates);

	if (dconf->buckets->list != NULL)
	{
		align_free(dconf->buckets->list);
		free(dconf->buckets->fb_bounds);
		free(dconf->buckets->logp);
		free(dconf->buckets->num);
		free(dconf->buckets);
	}

	//support data on the poly currently being sieved
	free(dconf->curr_poly->gray);
	free(dconf->curr_poly->nu);
	free(dconf->curr_poly->qlisort);
	zFree(&dconf->curr_poly->poly_a);
	zFree(&dconf->curr_poly->poly_b);
	zFree(&dconf->curr_poly->poly_c);
	free(dconf->curr_poly);

	for (i=0;i<MAX_A_FACTORS;i++)
		zFree(&dconf->Bl[i]);
	free(dconf->Bl);

	//workspace bigints
	zFree(&dconf->qstmp1);
	zFree(&dconf->qstmp2);
	zFree(&dconf->qstmp3);
	zFree(&dconf->qstmp4);
	zFree32(&dconf->qstmp32);
	align_free(dconf->mask);

	//free sieve scan report stuff
	free(dconf->reports);
	for (i=0; i<100; i++)
		zFree32(&dconf->Qvals[i]);
	free(dconf->Qvals);
	free(dconf->valid_Qs);
	free(dconf->smooth_num);

	//free post-processed relations
	//for (i=0; (uint32)i < dconf->buffered_rels; i++)
	//	free(dconf->relation_buf[i].fb_offsets);
	//free(dconf->relation_buf);

	return 0;
}

void free_filter_vars(static_conf_t *sconf)
{
	int i;

	//list of a values used first to track all a coefficients
	//generated during sieving for duplication, then used again during
	//filtering
	for (i=0; (uint32)i < sconf->total_poly_a; i++)
		zFree(&sconf->poly_a_list[i]);
	free(sconf->poly_a_list);

	//cycle table stuff created at the beginning of the factorization
	free(sconf->cycle_hashtable);
	free(sconf->cycle_table);

	//list of polys used in filtering and sqrt
	if (sconf->curr_b != NULL)
	{
		for (i=0;(uint32)i < sconf->bpoly_alloc;i++)
			zFree(&sconf->curr_b[i]);
		free(sconf->curr_b);
	}

	if (sconf->poly_list != NULL)
	{
		for (i=0;(uint32)i < sconf->poly_list_alloc;i++)
			zFree(&sconf->poly_list[i].b);
		free(sconf->poly_list);
	}

	if (sconf->relation_list != NULL)
	{
		//free post-processed relations
		for (i=0; (uint32)i < sconf->num_relations; i++)
			free(sconf->relation_list[i].fb_offsets);
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
	z tmp1, tmp2;

	zInit(&tmp1);
	zInit(&tmp2);

	//current poly info used during filtering
	free(sconf->curr_poly->gray);
	free(sconf->curr_poly->nu);
	free(sconf->curr_poly->qlisort);
	zFree(&sconf->curr_poly->poly_a);
	zFree(&sconf->curr_poly->poly_b);
	zFree(&sconf->curr_poly->poly_c);
	free(sconf->curr_poly);
	zFree(&sconf->curr_a);	
	free(sconf->modsqrt_array);
	align_free(sconf->factor_base->list->prime);
	align_free(sconf->factor_base->list->small_inv);
	align_free(sconf->factor_base->list->correction);
	align_free(sconf->factor_base->list->logprime);
	align_free(sconf->factor_base->list);
	free(sconf->factor_base);

	//while freeing the list of factors, divide them out of the input
	zShortDiv(&sconf->n,sconf->multiplier,&sconf->n);
	zCopy(&sconf->n,&tmp1);
	for (i=0;i<sconf->factor_list.num_factors;i++)
	{
		z tmpz;
		char buf[32 * MAX_MP_WORDS+1];
		zInit(&tmpz);		

		mp_t2z(&sconf->factor_list.final_factors[i]->factor,&tmpz);
//		printf("extracting z = %s\n",z2decstr(&tmpz,&gstr1));

		zDiv(&tmp1,&tmpz,&sconf->n,&tmp2);
		if (isPrime(&tmpz))
		{
			tmpz.type = PRP;
			add_to_factor_list(sconf->obj, &tmpz);
		}
		else
		{
			tmpz.type = COMPOSITE;
			add_to_factor_list(sconf->obj, &tmpz);
		}
		zFree(&tmpz);
		free(sconf->factor_list.final_factors[i]);
		zCopy(&sconf->n,&tmp1);
	}

	zCopy(&sconf->n,&sconf->obj->qs_obj.n);
	zFree(&sconf->sqrt_n);
	zFree(&sconf->n);
	zFree(&sconf->max_fb2);
	zFree(&sconf->large_prime_max2);
	zFree(&sconf->target_a);
	zFree(&tmp1);
	zFree(&tmp2);

	//free(sconf->obj->qs_obj.savefile.name);
	qs_savefile_free(&sconf->obj->qs_obj.savefile);

	return 0;
}

uint8 choose_multiplier_siqs(uint32 B, z *n) 
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
		uint8 knmod8 = (uint8)((curr_mult * n->val[0]) % 8);
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
		uint32 modp = (uint32)zShortMod(n, prime);

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

