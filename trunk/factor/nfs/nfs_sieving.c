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

#include "nfs.h"

#ifdef USE_NFS

char* time_from_secs(unsigned long time)
{ // perhaps should go in a different file
	char *str, buf[80];
	unsigned long d;
	str = (char*)malloc(80); // can't return a char[] because it's local
	// (but how to free it :P)
	if( !str )
	{
		printf("malloc failed\n");
		exit(-1);
	}
	
	if( time > 3600*24 )
	{
		d = time / (3600*24);
		time %= (3600*24);
		sprintf(str, "%lu day%s ", d, d>1?"s":"");
	}
	if( time > 3600 )
	{
		d = time / 3600;
		time %= 3600;
		sprintf(buf, "%luh ", d);
		strcat(str, buf);
	}
	if( time > 60 )
	{
		d = time / 60;
		time %= 60;
		sprintf(buf, "%lum ", d);
		strcat(str, buf);
	}
	sprintf(buf, "%lus", time);
	strcat(str, buf);
	return str;
}

int test_sieve(fact_obj_t* fobj, void* args, int njobs, int are_files)
/* if(are_files), then treat args as a char** list of (external) polys
 * else args is a nfs_job_t* array of job structs
 * the latter is preferred */
// @ben: the idea is a new yafu function "testsieve(n, ...)" where args are
// an arbitrary list of files of external polys
{
	uint32 count;
	int i, minscore_id = 0;
	double* score = (double*)malloc(njobs * sizeof(double));
	double t_time, min_score = 999999999.;
	char orig_name[GSTR_MAXSIZE]; // don't clobber fobj->nfs_obj.job_infile
	char* time;
	
	char** filenames; // args
	nfs_job_t* jobs; // args
	
	struct timeval stop, stop2;	// stop time of this job
	struct timeval start, start2;	// start time of this job
	TIME_DIFF *	difference;

	if( score == NULL )
	{
		printf("Couldn't alloc memory!\n");
		exit(-1);
	}

	// check to make sure we can find ggnfs sievers
	if (check_for_sievers(fobj, 0) == 1)
	{
		printf("test: can't find ggnfs lattice sievers - aborting test sieving\n");
		return -1;
	}

	gettimeofday(&start2, NULL);

	strcpy(orig_name, fobj->nfs_obj.job_infile);
	// necessary because parse/fill_job_file() get filename from fobj
	if( are_files )
	{ // read files into job structs (get fblim)

		filenames = (char**) args;
		jobs = (nfs_job_t*)malloc(njobs * sizeof(nfs_job_t));
		if( jobs == NULL )
		{
			printf("Couldn't alloc memory!\n");
			exit(-1);
		}
		memset(jobs, 0, njobs*sizeof(nfs_job_t));

		for(i = 0; i < njobs; i++)
		{
			uint32 missing_params;
			strcpy(fobj->nfs_obj.job_infile, filenames[i]);

			missing_params = parse_job_file(fobj, jobs+i); // get fblim
			get_ggnfs_params(fobj, jobs+i); // get siever
			if( missing_params )
			{
				if( VFLAG >= 0 )
					printf("test: warning: \"%s\" is missing some paramters (%#X). filling them.\n",
						filenames[i], missing_params);
				fill_job_file(fobj, jobs+i, missing_params);
			}
			/*if( !jobs[i].poly )
			{ // not snfs (detected by parse_job_file())
				printf("branch taken\n");
				jobs[i].poly = (mpz_polys_t*)malloc(sizeof(mpz_polys_t));
				if( !jobs[i].poly )
				{
					printf("derp: out of memory!\n");
					exit(-1);
				}
				mpz_polys_init(jobs[i].poly);
				// if user doesn't specify, go with algebraic *shrug*
				jobs[i].poly->side = fobj->nfs_obj.sq_side < 0 ? RATIONAL_SPQ : ALGEBRAIC_SPQ;
			}*/
		}
	}
	else
	{ // create poly files
		jobs = (nfs_job_t *) args;
		filenames = (char**)malloc(njobs*sizeof(char*));
		if( !filenames )
		{
			printf("malloc derped it up!\n");
			exit(-1);
		}
		for(i = 0; i < njobs; i++)
		{
			FILE* out;
			filenames[i] = (char*)malloc(GSTR_MAXSIZE);
			if( !filenames[i] )
			{
				printf("malloc failed\n");
				exit(-1);
			}
			sprintf(filenames[i], "test-sieve-%d.poly", i);
			out = fopen(filenames[i], "w");
			if( !out )
			{
				printf("test: couldn't open %s for writing, aborting test sieve\n", filenames[i]);
				return -1;
			}
			if( jobs[i].snfs )
				print_snfs(jobs[i].snfs, out);
			else
			{
				gmp_fprintf(out, "n: %Zd\n", fobj->nfs_obj.gmp_n);
				print_poly(jobs[i].poly, out);
			}
			fclose(out);
			strcpy(fobj->nfs_obj.job_infile, filenames[i]);
			fill_job_file(fobj, jobs+i, PARAM_FLAG_ALL);
		} 
		// that seems like a lot more code than it should be
	}
	strcpy(fobj->nfs_obj.job_infile, orig_name);

	// now we can get to the actual testing
	for(i = 0; i < njobs; i++)
	{
		char syscmd[GSTR_MAXSIZE], tmpbuf[GSTR_MAXSIZE], side[32];
		FILE* in;
		double est;

		if( jobs[i].poly->side == RATIONAL_SPQ)
		{
			sprintf(side, "rational");
			jobs[i].startq = jobs[i].rlim; // no reason to test sieve *inside* the fb
		}
		else
		{
			sprintf(side, "algebraic");
			jobs[i].startq = jobs[i].alim; // ditto
		}

		//create the afb/rfb - we don't want the time it takes to do this to
		//pollute the sieve timings
		sprintf(syscmd, "%s -b %s -k -c 0 -F", jobs[i].sievername, filenames[i]);
		if (VFLAG > 0) printf("\ntest: generating factor bases\n");
		gettimeofday(&start, NULL);
		system(syscmd);
		gettimeofday(&stop, NULL);
		difference = my_difftime (&start, &stop);
		t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);
		if (VFLAG > 0) printf("test: fb generation took %6.4f seconds\n", t_time);
		MySleep(.1);
		sprintf(tmpbuf, "%s.afb.0", filenames[i]);
		remove(tmpbuf);
		// are we really supposed to be removing the factorbase files before the test?

		//start the test
		sprintf(syscmd,"%s%s -%c %s -f %u -c %u -o %s.out",
			jobs[i].sievername, VFLAG>0?" -v":"", side[0], filenames[i], jobs[i].startq, 5000, filenames[i]);
		if (VFLAG > 0) printf("test: commencing test sieving of polynomial %d on the %s side over range %u-%u\n", i, 
			side, jobs[i].startq, jobs[i].startq + 5000);
		gettimeofday(&start, NULL);
		system(syscmd);
		gettimeofday(&stop, NULL);
		difference = my_difftime (&start, &stop);
		t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);

		//count relations
		sprintf(tmpbuf, "%s.out", filenames[i]);
		in = fopen(tmpbuf, "r");
		if( !in )
			score[i] = 999999999.;
		else
		{
			count = 0;
			while( fgets(syscmd, GSTR_MAXSIZE, in) )
				count++;
			fclose(in);

			score[i] = t_time / count;

			est = (score[i] * jobs[i].min_rels * 1.25) / THREADS; 
			// be conservative about estimates
		}

		if( score[i] < min_score )
		{
			minscore_id = i;
			min_score = score[i];
			if (VFLAG > 0) printf("test: new best score of %1.6f sec/rel, estimated total sieving time = %s (with %d threads)\n", 
				score[i], time=time_from_secs((unsigned long)est), THREADS); 

     		}
		else
			if (VFLAG > 0) printf("test: score was %1.6f sec/rel, estimated total sieving time = %s (with %d threads)\n", 
				score[i], time=time_from_secs((unsigned long)est), THREADS);


		remove(tmpbuf); // clean up after ourselves
		sprintf(tmpbuf, "%s", filenames[i]);
		remove(tmpbuf);
		free(time);
	}

	// clean up memory allocated
	if( are_files )
	{
		for(i = 0; i < njobs; i++)
		{
			mpz_polys_free(jobs[i].poly);
			free(jobs[i].poly);
		}
		free(jobs);
	}
	else
	{
		for(i = 0; i < njobs; i++)
			free(filenames[i]);
		free(filenames);
	}

	gettimeofday(&stop2, NULL);
	difference = my_difftime (&start2, &stop2);
	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);			
	if (VFLAG > 0) printf("test: test sieving took %1.2f seconds\n", t_time);

	return minscore_id;
}

void do_sieving(fact_obj_t *fobj, nfs_job_t *job)
{
	nfs_threaddata_t *thread_data;		//an array of thread data objects
	int i;
	FILE *logfile;

	thread_data = (nfs_threaddata_t *)malloc(THREADS * sizeof(nfs_threaddata_t));
	for (i=0; i<THREADS; i++)
	{
		sprintf(thread_data[i].outfilename, "rels%d.dat", i);
		thread_data[i].job.poly = job->poly; // no sense copying the whole struct
		thread_data[i].job.rlim = job->rlim;
		thread_data[i].job.alim = job->alim;
		thread_data[i].job.rlambda = job->rlambda;
		thread_data[i].job.alambda = job->alambda;
		thread_data[i].job.lpbr = job->lpbr;
		thread_data[i].job.lpba = job->lpba;
		thread_data[i].job.mfbr = job->mfbr;
		thread_data[i].job.mfba = job->mfba;
		if (fobj->nfs_obj.rangeq > 0)
			thread_data[i].job.qrange = ceil((double)fobj->nfs_obj.rangeq / (double)THREADS);
		else
			thread_data[i].job.qrange = ceil((double)job->qrange / (double)THREADS);
		thread_data[i].job.min_rels = job->min_rels;
		thread_data[i].job.current_rels = job->current_rels;
		thread_data[i].siever = fobj->nfs_obj.siever;
		thread_data[i].job.startq = job->startq;
		strcpy(thread_data[i].job.sievername, job->sievername);

		thread_data[i].tindex = i;
		thread_data[i].is_poly_select = 0;
		thread_data[i].fobj = fobj;
	}

	/* activate the threads one at a time. The last is the
			master thread (i.e. not a thread at all). */		
	for (i = 0; i < THREADS - 1; i++)
		nfs_start_worker_thread(thread_data + i, 0);

	nfs_start_worker_thread(thread_data + i, 1);			

	// load threads with work
	for (i = 0; i < THREADS; i++) 
	{
		nfs_threaddata_t *t = thread_data + i;
		t->job.startq = job->startq;
		job->startq += t->job.qrange;
	}

	logfile = fopen(fobj->flogname, "a");
	if (logfile == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open yafu logfile for appending\n");
	}
	else
	{
		logprint(logfile, "nfs: commencing lattice sieving with %d threads\n",THREADS);
		fclose(logfile);
	}

	// create a new lasieve process in each thread and watch it
	for (i = 0; i < THREADS; i++) 
	{
		nfs_threaddata_t *t = thread_data + i;

		if (i == THREADS - 1) {				
			lasieve_launcher(t);
		}
		else {
			t->command = NFS_COMMAND_RUN;
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
		nfs_threaddata_t *t = thread_data + i;

		if (i < THREADS - 1) {
#if defined(WIN32) || defined(_WIN64)
			WaitForSingleObject(t->finish_event, INFINITE);
#else
			pthread_mutex_lock(&t->run_lock);
			while (t->command != NFS_COMMAND_WAIT)
				pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
		}
	}
		
	//combine output
	for (i = 0; i < THREADS; i++) 
	{
		nfs_threaddata_t *t = thread_data + i;
		savefile_concat(t->outfilename,fobj->nfs_obj.outputfile,fobj->nfs_obj.mobj);

		// accumulate relation counts
		job->current_rels += thread_data[i].job.current_rels;

		remove(thread_data[i].outfilename);			
	}

	//stop worker threads
	for (i=0; i<THREADS - 1; i++)
	{
		nfs_stop_worker_thread(thread_data + i, 0);
	}
	nfs_stop_worker_thread(thread_data + i, 1);

	//free the thread structure
	free(thread_data);

	return;
}

void *lasieve_launcher(void *ptr)
{
	//top level sieving function which performs all work for a single
	//new a coefficient.  has pthread calling conventions, meant to be
	//used in a multi-threaded environment
	nfs_threaddata_t *thread_data = (nfs_threaddata_t *)ptr;
	fact_obj_t *fobj = thread_data->fobj;
	char syscmd[GSTR_MAXSIZE], tmpstr[GSTR_MAXSIZE], side[GSTR_MAXSIZE];	
	FILE *fid;
	int cmdret;

	sprintf(side, (thread_data->job.poly->side == ALGEBRAIC_SPQ) ? 
				"algebraic" : "rational"); // gotta love ?:

	//remove any temporary relation files
	remove(thread_data->outfilename);
		
	//start ggnfs binary
	sprintf(syscmd,"%s%s -%c %s -f %u -c %u -o %s -n %d",
			thread_data->job.sievername, VFLAG>0?" -v":"", *side, // hehe (*side == side[0])
			fobj->nfs_obj.job_infile, thread_data->job.startq, 
			thread_data->job.qrange, thread_data->outfilename, thread_data->tindex);

	if (VFLAG >= 0)
	{
		printf("nfs: commencing %s side lattice sieving over range: %u - %u\n",
			side, thread_data->job.startq, thread_data->job.startq + thread_data->job.qrange);
	}
	if (VFLAG > 1) printf("syscmd: %s\n", syscmd);
	if (VFLAG > 1) fflush(stdout);
	cmdret = system(syscmd);

	// a ctrl-c abort signal is caught by the system command, and nfsexit never gets called.
	// so check for abnormal exit from the system command.
	// -1073741819 is apparently what ggnfs returns when it crashes, which
	// we don't want to interpret as an abort.
	if (!((cmdret == 0) || cmdret == -1073741819))
		if( NFS_ABORT < 1 )
			NFS_ABORT = 1;

	// count the relations produced
	MySleep(100);
	fid = fopen(thread_data->outfilename,"r");
	thread_data->job.current_rels = 0;
	while (fgets(tmpstr, GSTR_MAXSIZE, fid) != NULL)
		thread_data->job.current_rels++;
	fclose(fid);

	return 0;
}

#endif
