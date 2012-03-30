#include "nfs.h"

void do_msieve_polyselect(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job, 
	mp_t *mpN, factor_list_t *factor_list)
{
	FILE *logfile;
	uint32 flags;

	//an array of thread data objects
	nfs_threaddata_t *thread_data;		

	// thread work-queue controls
    int threads_working = 0;
    int *thread_queue, *threads_waiting;
#if defined(WIN32) || defined(_WIN64)
	HANDLE queue_lock;
	HANDLE *queue_events = NULL;
#else
    pthread_mutex_t queue_lock;
    pthread_cond_t queue_cond;
#endif

	int i,j,is_startup;
	char syscmd[1024];
	char master_polyfile[80];
	uint64 start = 0, range = 0;
	uint32 deadline, estimated_range_time, this_range_time, total_time, num_ranges;
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	TIME_DIFF *	difference;
	double t_time;

	//file into which we will combine all of the thread results
	sprintf(master_polyfile,"%s.p",fobj->nfs_obj.outputfile);

	if (job->last_leading_coeff == 0)
	{
		//make sure we are starting from scratch
		remove(master_polyfile);
	}

	// figure out how long poly selection should take
	i = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 2);
	for (j = 0; j < NUM_TIME_LIMITS; j++) {
		if (i < time_limits[j].bits)
			break;
	}
	if (j == NUM_TIME_LIMITS) {
		deadline = time_limits[j-1].seconds;
	}
	else {
		const poly_deadline_t *low = &time_limits[j-1];
		const poly_deadline_t *high = &time_limits[j];
		uint32 dist = high->bits - low->bits;
		deadline = (uint32)(
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
			if (VFLAG > 0)
				printf("resuming poly search: reducing deadline by %u seconds\n",
					job->poly_time);

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
		deadline /= THREADS;
	}

	if (VFLAG > 0)
		printf("setting deadline of %u seconds\n",deadline);

	//start a counter for the poly selection
	gettimeofday(&startt, NULL);

	//init each thread data structure with info needed to do poly search on a range
	thread_data = (nfs_threaddata_t *)malloc(THREADS * sizeof(nfs_threaddata_t));

	// allocate the queue of threads waiting for work
    thread_queue = (int *)malloc(THREADS * sizeof(int));
    threads_waiting = (int *)malloc(sizeof(int));

	if (THREADS > 1)
	{
#if defined(WIN32) || defined(_WIN64)
		queue_lock = CreateMutex( 
			NULL,              // default security attributes
			FALSE,             // initially not owned
			NULL);             // unnamed mutex
		queue_events = (HANDLE *)malloc(THREADS * sizeof(HANDLE));
#else
		pthread_mutex_init(&queue_lock, NULL);
		pthread_cond_init(&queue_cond, NULL);
#endif
	}

	//set flags to do poly selection
	flags = 0;
	if (VFLAG > 0)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_POLY1;
	flags = flags | MSIEVE_FLAG_NFS_POLY2;		
	
	for (i=0; i<THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;		
		t->fobj = fobj;

		// create thread data with dummy range for now
		init_poly_threaddata(t, obj, mpN, factor_list, i, flags,
			(uint64)1, (uint64)1001);

		//give this thread a unique index
		t->tindex = i;
		t->is_poly_select = 1;

		t->thread_queue = thread_queue;
        t->threads_waiting = threads_waiting;

		if (THREADS > 1)
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
	logfile = fopen(fobj->flogname, "a");
	if (logfile == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open yafu logfile for appending\n");
	}
	else
	{
		logprint(logfile, "nfs: commencing poly selection with %d threads\n",THREADS);
		logprint(logfile, "nfs: setting deadline of %u seconds\n",deadline);
		fclose(logfile);
	}

	// determine the start and range values
	get_polysearch_params(fobj, &start, &range);	

	if (THREADS > 1)
	{
		// Activate the worker threads one at a time. 
		// Initialize the work queue to say all threads are waiting for work
		for (i = 0; i < THREADS; i++) 
		{
			nfs_start_worker_thread(thread_data + i, 2);
			thread_queue[i] = i;
		}
	}
    *threads_waiting = THREADS;

	if (THREADS > 1)
	{
#if defined(WIN32) || defined(_WIN64)
		// nothing
#else
		pthread_mutex_lock(&queue_lock);
#endif
	}	

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
						
			if (THREADS > 1)
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

			// pointer to this thread's data
			t = thread_data + tid;

			//check the total time spent so far
			gettimeofday(&stopt, NULL);
			difference = my_difftime (&startt, &stopt);

			t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
			free(difference);	

			//update the estimated range time			
			difference = my_difftime (&t->thread_start_time, &stopt);

			this_range_time = ((double)difference->secs + (double)difference->usecs / 1000000);			
			free(difference);	

			if (!is_startup)
			{
				total_time += this_range_time;
				num_ranges++;
				estimated_range_time = this_range_time; //(uint32)(total_time / (double)num_ranges);
			}			

			if (!is_startup)
			{
				FILE *fid;

				//if the main poly file doesn't exist yet, put the input number
				//at the top.  this is used to help restart jobs in the polyfind phase
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
					sprintf(syscmd,"cat %s.p >> %s",t->polyfilename,master_polyfile);
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
				system(syscmd);
#endif
				//then stick on the current total elasped time
				//this is used to help restart jobs in the polyfind phase
				fid = fopen(master_polyfile, "a");
				fprintf(fid, "time: %u\n", total_time);
				fclose(fid);

			}

			// remove each thread's .p file after it's copied
			sprintf(syscmd, "%s.p",t->polyfilename); 
			remove(syscmd);	

			// also remove the temporary log file
			sprintf(syscmd,"%s.%d",fobj->nfs_obj.logfile,tid);
			remove(syscmd);

			// and the temporary fb file
			sprintf(syscmd,"%s.%d",fobj->nfs_obj.fbfile,tid);
			remove(syscmd);

			//free data
			free(t->logfilename);
			free(t->polyfilename);
			msieve_obj_free(t->obj);			
						
			// exit the loop without restarting if abort is asserted.
			if (NFS_ABORT)
				break;

			// if we can re-start the thread such that it is likely to finish before the
			// deadline, go ahead and do so.
			if ((uint32)t_time + estimated_range_time <= deadline)
			{

				// unless the user has specified a custom range search, in which case
				// this thread is done
				if ((fobj->nfs_obj.polyrange > 0) && !is_startup)
				{
					printf("thread %d finished custom range search\n",tid);
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
					if (THREADS > 1)
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

			if (THREADS == 1)
				*threads_waiting = 0;

		}

		// after starting all ranges for the first time, reset this flag
		if (is_startup)
			is_startup = 0;

		// if all threads are done, break out
		if (threads_working == 0)
			break;

		if (THREADS > 1)
		{
			// wait for a thread to finish and put itself in the waiting queue
#if defined(WIN32) || defined(_WIN64)
			j = WaitForMultipleObjects(
				THREADS,
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
	difference = my_difftime (&startt, &stopt);

	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);	

	if (fobj->nfs_obj.polyrange > 0)
	{		
		printf("custom range search complete in %6.4f seconds\n",t_time);
	}
	else if (VFLAG >= 0)
		printf("elapsed time of %6.4f seconds exceeds %u second deadline; poly select done\n",
			t_time,deadline);
	
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

	//stop worker threads
	for (i=0; i<THREADS; i++)
	{
		//static_conf->tot_poly += thread_data[i].dconf->tot_poly;
		if (THREADS > 1)
			nfs_stop_worker_thread(thread_data + i, 2);
	}

	//free the thread structure
	free(thread_data);		
	free(thread_queue);
    free(threads_waiting);

#if defined(WIN32) || defined(_WIN64)
	if (THREADS > 1)
		free(queue_events);
#endif

	if (!NFS_ABORT)
	{
		// parse the master .p file and find the best poly	
		find_best_msieve_poly(fobj, job, 1);
		//also create a .fb file
		ggnfs_to_msieve(fobj, job);
	}

	return;
}

void init_poly_threaddata(nfs_threaddata_t *t, msieve_obj *obj, 
	mp_t *mpN, factor_list_t *factor_list, int tid, uint32 flags, 
	uint64 start, uint64 stop)
{
	fact_obj_t *fobj = t->fobj;
	t->logfilename = (char *)malloc(80 * sizeof(char));
	t->polyfilename = (char *)malloc(80 * sizeof(char));
	t->fbfilename = (char *)malloc(80 * sizeof(char));
	sprintf(t->polyfilename,"%s.%d",fobj->nfs_obj.outputfile,tid);
	sprintf(t->logfilename,"%s.%d",fobj->nfs_obj.logfile,tid);
	sprintf(t->fbfilename,"%s.%d",fobj->nfs_obj.fbfile,tid);

	//make sure there isn't a preexisting fb file
	remove(t->fbfilename);

	//create an msieve_obj.  for poly select, the intermediate output file should be specified in
	//the savefile field		
	t->obj = msieve_obj_new(obj->input, flags, t->polyfilename, t->logfilename, t->fbfilename, 
		g_rand.low, g_rand.hi, (uint32)0, start, stop, 
		obj->cpu, (uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);

	//pointers to things that are static during poly select
	t->mpN = mpN;
	t->factor_list = factor_list;
	gettimeofday(&t->thread_start_time, NULL);

	return;
}

void get_polysearch_params(fact_obj_t *fobj, uint64 *start, uint64 *range)
{

	//search smallish chunks of the space in parallel until we've hit our deadline
	if (fobj->nfs_obj.polystart > 0)
		*start = fobj->nfs_obj.polystart;
	else
		*start = 1ULL;
	
	if (fobj->nfs_obj.polyrange > 0)
	{
		// search custom range of leading coefficient.  This effectively ignores the deadline
		// because the deadline is only checked after each range is done.  if we assign the
		// entire range up front, the deadline is never checked before polysearch finishes.
		*range = fobj->nfs_obj.polyrange / THREADS;

		// sanity check
		if (*range == 0) *range = 1;
	}
	else
	{
		// loop on a default range until the time deadline is reached
		*range = (uint64)fobj->nfs_obj.polybatch;
	}

	return;
}

void *polyfind_launcher(void *ptr)
{
	//top level function which performs poly search on a range of leading A coefficient values.
	//has pthread calling conventions, meant to be used in a multi-threaded environment
	nfs_threaddata_t *t = (nfs_threaddata_t *)ptr;

	//remove any temporary relation files
	remove(t->polyfilename);
		
	if (VFLAG >= 0)
	{
		printf("nfs: commencing polynomial search over range: %" PRIu64 " - %" PRIu64"\n",
			t->obj->nfs_lower, t->obj->nfs_upper);
		fflush(stdout);
	}

	//start polyfind
	factor_gnfs(t->obj, t->mpN, t->factor_list);	

	return 0;
}
