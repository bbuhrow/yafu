/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

       				   --bbuhrow@gmail.com 12/6/2012
----------------------------------------------------------------------*/

#include "nfs_impl.h"
#include "threadpool.h"

#ifdef USE_NFS

static const poly_deadline_t time_limits[] = {
    //  bits, seconds
        {248, 1 * 60},		// 74 digits
        {264, 2 * 60},		// 80 digits
        {304, 6 * 60},		// 92 digits
        {320, 15 * 60},		// 97 digits
        {348, 30 * 60},		// 105 digits
        {365, 1 * 3600},	// 110 digits
        {383, 2 * 3600},	// 116 digits
        {399, 4 * 3600},	// 120 digits
        {416, 8 * 3600},	// 126 digits
        {433, 16 * 3600},	// 131 digits
        {449, 32 * 3600},	// 135 digits
        {466, 64 * 3600},	// 140 digits
        {482, 100 * 3600},	// 146 digits
        {498, 200 * 3600},	// 150 digits
        {514, 300 * 3600},	// 155 digits
};

#define NUM_TIME_LIMITS sizeof(time_limits)/sizeof(time_limits[0])

snfs_t * snfs_find_form(fact_obj_t *fobj)
{
	snfs_t* poly;

	// allow larger one to come in, because some polynomials forms are significantly reduced
	// (i.e., x^(2n))
	if (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 2) > 2*MAX_SNFS_BITS)
	{
		if (fobj->VFLAG >= 0)
			printf("nfs: n is too large for snfs, skipping snfs poly select\n");
		return NULL;
	}	

	poly = (snfs_t*)malloc(sizeof(snfs_t));
	if( !poly )
	{
		printf("gen: initial malloc failed\n");
		exit(-1);
	}
	snfs_init(poly);

	if (poly->form_type == SNFS_NONE)
	{
		if (fobj->VFLAG >= 0) printf("nfs: searching for brent special forms...\n");
		// if this is a factor() run, restore the original input number so that we 
		// can detect these forms
		if (fobj->autofact_obj.autofact_active)
		{
			mpz_set(fobj->nfs_obj.snfs_cofactor, fobj->nfs_obj.gmp_n);
			mpz_set(fobj->nfs_obj.gmp_n, fobj->N);
		}

		find_brent_form(fobj, poly);

		if (fobj->autofact_obj.autofact_active)
		{
			mpz_set(fobj->nfs_obj.gmp_n, fobj->nfs_obj.snfs_cofactor);
		}
	}

	if (poly->form_type == SNFS_NONE)
	{
		if (fobj->VFLAG >= 0) printf("nfs: searching for homogeneous cunningham special forms...\n");
		find_hcunn_form(fobj, poly);	
	}

	if (poly->form_type == SNFS_NONE)
	{
		if (fobj->VFLAG >= 0) printf("nfs: searching for XYYXF special forms...\n");
		find_xyyxf_form(fobj, poly);
	}

    if (poly->form_type == SNFS_NONE)
    {
        if (fobj->VFLAG >= 0) printf("nfs: searching for direct special forms...\n");
        find_direct_form(fobj, poly);
    }

	if (poly->form_type == SNFS_NONE)
	{
		printf("nfs: couldn't find special form\n");
		
		snfs_clear(poly);
		free(poly);
		return NULL;
	}

	return poly;
}

int snfs_choose_poly(fact_obj_t* fobj, nfs_job_t* job)
{
	snfs_t* poly, * polys = NULL;
	int i, npoly;
	FILE *f;
	nfs_job_t *jobs, *best;
	int retcode = 0;

	if ((poly = snfs_find_form(fobj)) == NULL)
	{
		// form detection failed
		job->snfs = NULL;
		return 1;
	}

	// once form detection is done, do some quick trial division
	mpz_set(fobj->div_obj.gmp_n, poly->n);
	zTrial(fobj);
	mpz_set(poly->n, fobj->div_obj.gmp_n);

	// with the form detected, create a good polynomial
	if (poly->form_type == SNFS_XYYXF)
		polys = gen_xyyxf_poly(fobj, poly, &npoly);
	else
		polys = gen_brent_poly(fobj, poly, &npoly);

	if (npoly == 0)
	{
		printf("nfs: no snfs polynomial with small coefficients found\n");
		job->snfs = NULL;
		snfs_clear(poly);
		free(poly);
		return 1;
	}
    else if (npoly == 1)
    {
        printf("nfs: found %d polynomial\n", npoly);
    }
    else if (npoly > 1)
    {
        printf("nfs: found %d polynomials, selecting best\n", npoly);
    }

	// we've now measured the difficulty for poly's of all common degrees possibly formed
	// in several different ways.  now we have a decision to make based largely on difficulty and 
	// degree.  We want to pick low difficulty, but only if the degree allows the norms on 
	// both sides to be approximately equal.  Sometimes multiple degrees satisfy this requirement
	// approximately equally in which case only test-sieving can really resolve the difference.
	snfs_scale_difficulty(polys, npoly, fobj->VFLAG);
    if (npoly > 1)
    {
        npoly = snfs_rank_polys(fobj, polys, npoly);
    }

	if (fobj->VFLAG > 0 && npoly > 1)
	{		
		int np = MIN(npoly, NUM_SNFS_POLYS);

        if (fobj->LOGFLAG)
        {
            f = fopen(fobj->flogname, "a");
        }
        else
        {
            f = NULL;
        }
		
		printf("\ngen: ========================================================\n"
			"gen: best %d polynomials:\n"
			"gen: ========================================================\n\n", np);
		
		if (f != NULL) logprint(f, "gen: best %d polynomials:\n", np);

		for (i=0; i<np; i++)
		{
			print_snfs(&polys[i], stdout);
			if (f != NULL) print_snfs(&polys[i], f);
		}

		if (f != NULL) fclose(f);
	}

	// initialize job structures to pass into test sieving.  test sieving may modify parameters
	// of the job to optimize rels/q.
	jobs = (nfs_job_t*)malloc(npoly * sizeof(nfs_job_t));
	if( !jobs )
	{
		printf("out of memory\n");
		exit(-1);
	}
	memset(jobs, 0, npoly * sizeof(nfs_job_t));
	
	// assign initial parameters
	for (i=0; i<npoly; i++)
	{
		jobs[i].poly = polys[i].poly;
		jobs[i].snfs = &polys[i];
		get_ggnfs_params(fobj, &jobs[i]);
		skew_snfs_params(fobj, &jobs[i]);
	}

    //printf("gnfs size = %d, size n + 3 = %d, snfs = %d\n", est_gnfs_size(&jobs[0]),
    //    (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) + 3), fobj->nfs_obj.snfs);

	if ((est_gnfs_size(&jobs[0]) > (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) + 3)))
	{
        if (fobj->nfs_obj.snfs)
        {
            // the input is probably faster by gnfs, and the user hasn't specifically chosen
            // snfs, so don't trial sieve and return the gnfs preference
            printf("gen: WARNING: input probably faster by gnfs but proceeding with snfs-specified job\n");
        }
        else
        {
            // faster by gnfs despite found snfs form, skip to gnfs processing.
            retcode = 1;
            copy_job(&jobs[0], job);
            goto cleanup;
        }
	}

	// return best job, which contains the best poly
	best = snfs_test_sieve(fobj, polys, MIN(NUM_SNFS_POLYS,npoly), jobs);

	if (fobj->VFLAG > 0)
	{
        if (fobj->LOGFLAG)
        {
            f = fopen(fobj->flogname, "a");
        }
        else
        {
            f = NULL;
        }

		printf("\ngen: ========================================================\n");
		printf("gen: selected polynomial:\n");
		printf("gen: ========================================================\n\n");

		if (f != NULL) logprint(f, "gen: selected polynomial:\n");

		print_snfs(best->snfs, stdout);
		if (f != NULL) print_snfs(best->snfs, f);
		if (f != NULL) fclose(f);
	}

	// copy the best job to the object that will be returned from this function
	copy_job(best, job);

	// re-compute the min_rels, in case that changed
	nfs_set_min_rels(job);

	// make the file and fill it
	snfs_make_job_file(fobj, job);
	fill_job_file(fobj, job, PARAM_FLAG_ALL);

	// also create a .fb file
	ggnfs_to_msieve(fobj, job);

	// set the starting q
	job->startq = job->snfs->poly->side == RATIONAL_SPQ ? job->rlim/2 : job->alim/2;

cleanup:
    for (i = 0; i < npoly; i++)
    {
        snfs_clear(&polys[i]);
    }
	snfs_clear(poly);
	free(poly);
	free(polys);
	free(jobs);

	return retcode;
}

void do_msieve_polyselect(fact_obj_t *fobj, msieve_obj *obj, nfs_job_t *job, 
	mp_t *mpN, factor_list_t *factor_list)
{
	FILE *logfile;
	uint32_t flags;
    double bestscore = 0.;
    double quality_mult = 1.;
    char quality[8];

	//an array of thread data objects
	nfs_threaddata_t *thread_data;

	// thread work-queue controls
	int threads_working = 0;
	int *thread_queue, *threads_waiting;
#if defined(WIN32) || defined(_WIN64)
	HANDLE queue_lock = NULL;
	HANDLE *queue_events = NULL;
#else
	pthread_mutex_t queue_lock;
	pthread_cond_t queue_cond;
#endif

	int i,j,is_startup;
	char syscmd[GSTR_MAXSIZE + 4];
	char master_polyfile[GSTR_MAXSIZE + 2];
	uint64_t start = 0, range = 0;
	uint32_t deadline, estimated_range_time, this_range_time, total_time, num_ranges;
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	double t_time;
    double e0;
    int sysreturn;

	//file into which we will combine all of the thread results
	snprintf(master_polyfile, GSTR_MAXSIZE + 2, "%s.p",fobj->nfs_obj.outputfile);

	if (job->last_leading_coeff == 0)
	{
		//make sure we are starting from scratch
		remove(master_polyfile);
	}

	// figure out how long poly selection should take
	i = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 2);

    if (i < 363) 		/* <= 110 digits */
    {
        fobj->nfs_obj.pref_degree = 4;
    }
    else if (i < 726) 		/* 110-220 digits */
    {
        fobj->nfs_obj.pref_degree = 5;
    }
    else				/* 220+ digits */
    {
        fobj->nfs_obj.pref_degree = 6;
    }

	for (j = 0; j < NUM_TIME_LIMITS; j++) 
    {
        if (i < time_limits[j].bits)
        {
            break;
        }
	}
	
    if (j == NUM_TIME_LIMITS) 
    {
		deadline = time_limits[j-1].seconds;
	}
	else 
    {
		const poly_deadline_t *low = &time_limits[j-1];
		const poly_deadline_t *high = &time_limits[j];
		uint32_t dist = high->bits - low->bits;
		deadline = (uint32_t)(
			 ((double)low->seconds * (high->bits - i) +
			  (double)high->seconds * (i - low->bits)) / dist);
	}    

	// initialize the variable tracking the total time spent (over all threads)
	// so far in the polynomial selection phase
	total_time = 0;

	// if we're resuming a job in this phase, this will be non-zero
	if (job->poly_time > 0)
	{
		if (job->poly_time >= deadline)
		{
			// we've already spent as much time as we need.  
			// parse the master .p file and find the best poly	
			find_best_msieve_poly(fobj, job, 1);
			//also create a .fb file
			ggnfs_to_msieve(fobj, job);

			return;
		}
		else
		{
			if (fobj->VFLAG > 0)
				printf("nfs: resuming poly search: reducing deadline by %u seconds\n",
					job->poly_time);

            if (fobj->LOGFLAG)
            {
                logfile = fopen(fobj->flogname, "a");
                if (logfile == NULL)
                {
                    printf("could not open yafu logfile for appending\n");
                    printf("fopen error: %s\n", strerror(errno));
                }
                else
                {
                    logprint(logfile, "nfs: resuming poly search, reducing deadline by %u seconds\n",
                        job->poly_time);
                    logprint(logfile, "nfs: resuming poly search, starting at last coefficient %u\n",
                        job->last_leading_coeff);
                    fclose(logfile);
                }
            }

			// pick up where we left off.
			total_time = job->poly_time;

			// also need to pro-rate the deadline, since that is what is used
			// to determine when to stop
			deadline -= job->poly_time;
		}
	}		

	if (fobj->nfs_obj.poly_option == 0)
	{
		// 'fast' search.  scale by number of threads
		deadline /= fobj->THREADS;
	}

    // if a command line option 'min', 'avg', or 'good' is
    // specified, then allow early abort of poly search when
    // a polynomial of corresponding quality is found, according
    // to this heuristic:
    // http://www.mersenneforum.org/showthread.php?t=16994
    { /* SB: tried L[1/3,c] fit; it is no better than this */
        int digits = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);
        int degree = fobj->nfs_obj.pref_degree;

        e0 = (digits >= 121) ? (0.0607 * digits + 2.25) :
            (0.0526 * digits + 3.23);

        if (degree == 4) e0 = 0.0625 * digits + 1.69;
        e0 = exp(-log(10) * e0);
#ifdef HAVE_CUDA
        e0 *= 1.15;
#endif
        /* seen exceptional polys with +40% but that's */
        /* rare. The fit is good for 88..232 digits */
    }

    if (fobj->VFLAG > 0)
    {
        printf("nfs: setting deadline of %u seconds\n", deadline);
        printf("nfs: expecting degree %d poly E from %.2le to > %.2le\n",
            fobj->nfs_obj.pref_degree, e0, 1.15 * e0);
    }

    if (fobj->nfs_obj.poly_option < 3)
    {
        // if we want to use the full time, set the 
        // quality mulitplier really high.  If we find
        // one above this value surely it is ok to stop,
        // even if e.g., 'deep' was specified.
        quality_mult = 1.4;
        strcpy(quality, "awesome");
    }
    if (fobj->nfs_obj.poly_option == 3)
    {
        quality_mult = 1.0;
        strcpy(quality, "min");
    }
    else if (fobj->nfs_obj.poly_option == 5)
    {
        quality_mult = 1.072;
        strcpy(quality, "good");        
    }
    else
    {
        quality_mult = 1.036;
        strcpy(quality, "avg");
    }

	//start a counter for the poly selection
	gettimeofday(&startt, NULL);

	//init each thread data structure with info needed to do poly search on a range
	thread_data = (nfs_threaddata_t *)malloc(fobj->THREADS * sizeof(nfs_threaddata_t));

	// allocate the queue of threads waiting for work
	thread_queue = (int *)malloc(fobj->THREADS * sizeof(int));
	threads_waiting = (int *)malloc(sizeof(int));

	if (fobj->THREADS > 1)
	{
#if defined(WIN32) || defined(_WIN64)
		queue_lock = CreateMutex( 
			NULL,              // default security attributes
			FALSE,             // initially not owned
			NULL);             // unnamed mutex
		queue_events = (HANDLE *)malloc(fobj->THREADS * sizeof(HANDLE));
#else
		pthread_mutex_init(&queue_lock, NULL);
		pthread_cond_init(&queue_cond, NULL);
#endif
	}

	//set flags to do poly selection
	flags = 0;
	if (fobj->VFLAG > 0)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_POLY1;
	flags = flags | MSIEVE_FLAG_NFS_POLYSIZE;
	flags = flags | MSIEVE_FLAG_NFS_POLYROOT;
	
	for (i = 0; i < fobj->THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;		
		t->fobj = fobj;

		// create thread data with dummy range for now
		init_poly_threaddata(t, obj, mpN, factor_list, i, flags,
			(uint64_t)1, (uint64_t)1001);

		//give this thread a unique index
		t->tindex = i;
		t->is_poly_select = 1;

		t->thread_queue = thread_queue;
        t->threads_waiting = threads_waiting;

		if (fobj->THREADS > 1)
		{
#if defined(WIN32) || defined(_WIN64)
			// assign a pointer to the mutex
			t->queue_lock = &queue_lock;
			t->queue_event = &queue_events[i];
#else
			t->queue_lock = &queue_lock;
			t->queue_cond = &queue_cond;
#endif
		}

	}

	// log that we are starting
    if (fobj->LOGFLAG)
    {
        logfile = fopen(fobj->flogname, "a");
        if (logfile == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open yafu logfile for appending\n");
        }
        else
        {
            logprint(logfile, "nfs: commencing poly selection with %d threads\n", fobj->THREADS);
            logprint(logfile, "nfs: setting deadline of %u seconds\n", deadline);
            logprint(logfile, "nfs: expecting degree %d poly E from %.2le to > %.2le\n",
                fobj->nfs_obj.pref_degree, e0, 1.15 * e0);
            if (fobj->nfs_obj.poly_option > 2)
            {
                logprint(logfile, "nfs: searching for %s quality poly E > %.2le\n",
                    quality, e0 * quality_mult);
            }
            fclose(logfile);
        }
    }

	// determine the start and range values
	get_polysearch_params(fobj, &start, &range);	

	if (fobj->THREADS > 1)
	{
		// Activate the worker threads one at a time. 
		// Initialize the work queue to say all threads are waiting for work
		for (i = 0; i < fobj->THREADS; i++)
		{
			nfs_start_worker_thread(thread_data + i, 2);
			thread_queue[i] = i;
		}
	}
	*threads_waiting = fobj->THREADS;

	if (fobj->THREADS > 1)
	{
#if defined(WIN32) || defined(_WIN64)
		// nothing
#else
		pthread_mutex_lock(&queue_lock);
#endif
	}	

    //printf("============= starting %d thread poly search with option %d  =================\n", 
    //    fobj->THREADS, fobj->nfs_obj.poly_option);

	estimated_range_time = 0;
	is_startup = 1;	
	num_ranges = 0;
	while (1)
	{

		// Process threads until there are no more waiting for their results to be collected
		while (*threads_waiting > 0)
		{
			// one or more threads have just finished 
			// (or, on first loop, nothing has started yet)
			int tid;
			nfs_threaddata_t *t;	

			// this thread is done, so decrement the count of working threads
			if (!is_startup)
				threads_working--;

			if (fobj->THREADS > 1)
			{
				// Pop a waiting thread off the queue (OK, it's stack not a queue)
#if defined(WIN32) || defined(_WIN64)

				WaitForSingleObject( 
					queue_lock,    // handle to mutex
					INFINITE);  // no time-out interval
#endif
  
				tid = thread_queue[--(*threads_waiting)];

#if defined(WIN32) || defined(_WIN64)
				ReleaseMutex(queue_lock);
#endif
			}
			else
				tid = 0;


            //printf("syncing poly from thread %d, %d threads waiting, %d threads working\n", 
            //    tid, *threads_waiting, threads_working);

			// pointer to this thread's data
			t = thread_data + tid;

			//check the total time spent so far
			gettimeofday(&stopt, NULL);
            t_time = ytools_difftime(&startt, &stopt);

			//update the estimated range time			
            this_range_time = ytools_difftime(&t->thread_start_time, &stopt);

			if (!is_startup)
			{
				total_time += this_range_time;
				num_ranges++;
				estimated_range_time = this_range_time; //(uint32)(total_time / (double)num_ranges);
			}			

			if (!is_startup)
			{
				FILE *fid;

				// if the main poly file doesn't exist yet, put the input number
				// at the top.  this is used to help restart jobs in the polyfind phase
				fid = fopen(master_polyfile, "r");
				if (fid == NULL)

				{
					fid = fopen(master_polyfile, "w");
					gmp_fprintf(fid, "n: %Zd\n", fobj->nfs_obj.gmp_n);
					fclose(fid);
				}
				else
					fclose(fid);

				// combine thread's output with main file
#if defined(WIN32)
				{
					int a;

					// test for cat
					sprintf(syscmd,"cat %s.p >> %s 2> nul",t->polyfilename,master_polyfile);
					a = system(syscmd);
	
					if (a)		
					{
						char tmp[80];
						sprintf(tmp, "%s.p",t->polyfilename);
						win_file_concat(tmp, master_polyfile);
					}
				}

#else
				sprintf(syscmd,"cat %s.p >> %s",t->polyfilename,master_polyfile);
				sysreturn = system(syscmd);
#endif
				// then stick on the current total elasped time
				// this is used to help restart jobs in the polyfind phase
				fid = fopen(master_polyfile, "a");
				fprintf(fid, "time: %u\n", total_time);
				fclose(fid);

			}

			// remove each thread's .p file after it's copied
			snprintf(syscmd, GSTR_MAXSIZE + 4, "%s.p",t->polyfilename);
			remove(syscmd);	

			// also remove the temporary log file
			snprintf(syscmd, GSTR_MAXSIZE + 4, "%s.%d",fobj->nfs_obj.logfile,tid);
			remove(syscmd);

			// and the temporary fb file
			snprintf(syscmd, GSTR_MAXSIZE + 4, "%s.%d",fobj->nfs_obj.fbfile,tid);
			remove(syscmd);

			// free data
			free(t->logfilename);
			free(t->polyfilename);
			msieve_obj_free(t->obj);

			// exit the loop without restarting if abort is asserted.
			if (NFS_ABORT)
				break;

            // if user has specified "good enough" option then check if
            // we've found one and stop if so
            if (!is_startup)
            {
                if (fobj->nfs_obj.poly_option > 2)
                {
                    bestscore = find_best_msieve_poly(fobj, job, 0);
                }
                else
                {
                    bestscore = 0.;
                }
            }

            if (fobj->VFLAG > 0)
            {
                printf("nfs: best score is currently %1.3le\n", bestscore);
            }

            if (bestscore < (e0 * quality_mult))
            {
                // if we can re-start the thread such that it is likely to finish before the
                // deadline, go ahead and do so.
                if ((uint32_t)t_time + estimated_range_time <= deadline)
                {

                    // unless the user has specified a custom range search, in which case
                    // this thread is done
                    if ((fobj->nfs_obj.polyrange > 0) && !is_startup)
                    {
                        printf("thread %d finished custom range search\n", tid);
                    }
                    else
                    {
                        init_poly_threaddata(t, obj, mpN, factor_list, tid, flags,
                            start, start + range);

                        if (fobj->nfs_obj.poly_option == 2)
                        {
                            // 'deep' searching assigns all threads to the same range
                            start = start;
                        }
                        else
                        {
                            // 'fast' and 'wide' searches aim to cover more ground with extra threads.
                            // assign next range
                            start += range;
                        }

                        // signal the job to start
                        if (fobj->THREADS > 1)
                        {
                            // send the thread a signal to start processing the poly we just generated for it
#if defined(WIN32) || defined(_WIN64)
                            thread_data[tid].command = NFS_COMMAND_RUN_POLY;
                            SetEvent(thread_data[tid].run_event);
#else
                            pthread_mutex_lock(&thread_data[tid].run_lock);
                            thread_data[tid].command = NFS_COMMAND_RUN_POLY;
                            pthread_cond_signal(&thread_data[tid].run_cond);
                            pthread_mutex_unlock(&thread_data[tid].run_lock);
#endif
                        }

                        // this thread is now busy, so increment the count of working threads
                        threads_working++;
                    }
                }
            }
            else
            {
                // announce we are finishing and don't restart the thread.
                if (fobj->VFLAG > 0)
                {
                    printf("nfs: found poly better than %s quality\n", quality);
                }
            }

			if (fobj->THREADS == 1)
				*threads_waiting = 0;

		}

		// after starting all ranges for the first time, reset this flag
        if (is_startup)
        {
            is_startup = 0;
        }

		// if all threads are done, break out
		if (threads_working == 0)
			break;

		if (fobj->THREADS > 1)
		{
			// wait for a thread to finish and put itself in the waiting queue
#if defined(WIN32) || defined(_WIN64)
			j = WaitForMultipleObjects(
                fobj->THREADS,
				queue_events,
				FALSE,
				INFINITE);
#else
			pthread_cond_wait(&queue_cond, &queue_lock);
#endif
		}
		else
		{
			//do some work
			nfs_threaddata_t *t = thread_data + 0;

			polyfind_launcher(t);
			*threads_waiting = 1;
		}

	}	

	gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);

	if (fobj->nfs_obj.polyrange > 0)
	{		
		printf("custom range search complete in %6.4f seconds\n",t_time);
	}
    else if (fobj->VFLAG >= 0)
    {
        printf("elapsed time: %6.4f seconds (%u second deadline); poly select done\n",
            t_time, deadline);
    }
	
    if (fobj->LOGFLAG)
    {
        logfile = fopen(fobj->flogname, "a");
        if (logfile == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open yafu logfile for appending\n");
        }
        else
        {
            logprint(logfile, "nfs: completed %u ranges of size %" PRIu64 " in %6.4f seconds\n",
                num_ranges, range, t_time);
            fclose(logfile);
        }
    }

	//stop worker threads
	for (i = 0; i < fobj->THREADS; i++)
	{
		//static_conf->tot_poly += thread_data[i].dconf->tot_poly;
		if (fobj->THREADS > 1)
			nfs_stop_worker_thread(thread_data + i, 2);
	}

	//free the thread structure
	free(thread_data);		
	free(thread_queue);
	free(threads_waiting);

#if defined(WIN32) || defined(_WIN64)
	if (fobj->THREADS > 1)
		free(queue_events);
#endif

	if (!NFS_ABORT)
	{
		// parse the master .p file and find the best poly	
		find_best_msieve_poly(fobj, job, 1);
		// also create a .fb file
		ggnfs_to_msieve(fobj, job);
	}

	return;
}

void init_poly_threaddata(nfs_threaddata_t *t, msieve_obj *obj, 
	mp_t *mpN, factor_list_t *factor_list, int tid, uint32_t flags,
	uint64_t start, uint64_t stop)
{
	fact_obj_t *fobj = t->fobj;
	char *nfs_args = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	t->logfilename = (char *)malloc(80 * sizeof(char));
	t->polyfilename = (char *)malloc(80 * sizeof(char));
	t->fbfilename = (char *)malloc(80 * sizeof(char));
	sprintf(nfs_args, "min_coeff=%" PRIu64 " max_coeff=%" PRIu64, start, stop);
	sprintf(t->polyfilename,"%s.%d",fobj->nfs_obj.outputfile,tid);
	sprintf(t->logfilename,"%s.%d",fobj->nfs_obj.logfile,tid);
	sprintf(t->fbfilename,"%s.%d",fobj->nfs_obj.fbfile,tid);

	//make sure there isn't a preexisting fb file
	remove(t->fbfilename);

	//create an msieve_obj.  for poly select, the intermediate output file should be specified in
	//the savefile field
	t->obj = msieve_obj_new(obj->input, flags, t->polyfilename, t->logfilename, t->fbfilename, 
		fobj->seed1, fobj->seed2, (uint32_t)0,
		9, (uint32_t)fobj->L1CACHE, (uint32_t)fobj->L2CACHE,
        (uint32_t)fobj->THREADS, (uint32_t)0, nfs_args);

	//pointers to things that are static during poly select
	t->mpN = mpN;
	t->factor_list = factor_list;
	gettimeofday(&t->thread_start_time, NULL);

	return;
}

void get_polysearch_params(fact_obj_t *fobj, uint64_t*start, uint64_t*range)
{

	//search smallish chunks of the space in parallel until we've hit our deadline
	if (fobj->nfs_obj.polystart > 0)
		*start = fobj->nfs_obj.polystart;
	else
		*start = 8192ULL;		// default leading coefficient
	
	if (fobj->nfs_obj.polyrange > 0)
	{
		// search custom range of leading coefficient.  This effectively ignores the deadline
		// because the deadline is only checked after each range is done.  if we assign the
		// entire range up front, the deadline is never checked before polysearch finishes.
		*range = fobj->nfs_obj.polyrange / fobj->THREADS;

		// sanity check
		if (*range == 0) *range = 1;
	}
	else
	{
		// loop on a default range until the time deadline is reached
		*range = (uint64_t)fobj->nfs_obj.polybatch;
	}

	return;
}

void *polyfind_launcher(void *ptr)
{
	// top level function which performs poly search on a range of leading A coefficient values.
	// has pthread calling conventions, meant to be used in a multi-threaded environment
	nfs_threaddata_t *t = (nfs_threaddata_t *)ptr;

	// remove any temporary relation files
	remove(t->polyfilename);
    MySleep(100);

	if (t->fobj->VFLAG >= 0)
	{
		uint64_t start, stop; // one instance where the new msieve api is rather a pain
		sscanf(t->obj->nfs_args, "min_coeff=%" PRIu64 " max_coeff=%" PRIu64, &start, &stop);
		printf("nfs: thread %d commencing polynomial search over range: %" PRIu64 " - %" PRIu64"\n",
			t->tindex, start, stop);
		fflush(stdout);
	}

	// start polyfind
	factor_gnfs(t->obj, t->mpN, t->factor_list);

    printf("**** finished poly work in thread %d\n", t->tindex);

	return 0;
}

#endif
