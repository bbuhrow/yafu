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

// potential enhancements:

// 2)
// run a series of test "refinements" based on difficulty.  The idea
// being twofold: a) for lower difficulty even a 2k range of spq might
// be overkill and b) for some polys the differences in sieve speed
// is obvious after only a few hundred spq.  So run a very small range
// first, and if we have a clear winner go with it.  if not, run more
// tests depending on both the "closeness" of the testing and the
// difficulty.
// 3)
// in windows, we count the number of ctrl-c's and force quit after 2.
// this defeats ctrl-c'ing out of more than one test sieve.  while test
// sieving, on windows, we need to allow more than two ctrl-c's

int test_sieve(fact_obj_t* fobj, void* args, int njobs, int are_files)
/* if(are_files), then treat args as a char** list of (external) polys
 * else args is a nfs_job_t* array of job structs
 * the latter is preferred */
// @ben: the idea is a new yafu function "testsieve(n, ...)" where args are
// an arbitrary list of files of external polys
{
	uint32_t count;
	int i, minscore_id = 0;
	double* score = (double*)malloc(njobs * sizeof(double));
	double t_time, min_score = 999999999.;
	char orig_name[GSTR_MAXSIZE]; // don't clobber fobj->nfs_obj.job_infile
	char time[80];
	uint32_t spq_range = 1000, actual_range;
	FILE *flog = NULL;
    int sysreturn;

	char** filenames; // args
	nfs_job_t* jobs; // args
	
	struct timeval stop, stop2;	// stop time of this job
	struct timeval start, start2;	// start time of this job

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
			uint32_t missing_params;
			strcpy(fobj->nfs_obj.job_infile, filenames[i]);

			missing_params = parse_job_file(fobj, jobs+i); // get fblim
			get_ggnfs_params(fobj, jobs+i); // get siever
			if( missing_params )
			{
				if( fobj->VFLAG >= 0 )
					printf("test: warning: \"%s\" is missing some paramters (%#X). filling them.\n",
						filenames[i], missing_params);
				fill_job_file(fobj, jobs+i, missing_params);
			}
			// adjust a/rlim, lpbr/a, and mfbr/a if advantageous
			skew_snfs_params(fobj, jobs+i);
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

        if (fobj->VFLAG > 0)
            printf("\ntest: trial sieving %s\n", filenames[i]);

		// should probably scale the range of special-q to test based
		// on input difficulty, but not sure how to do that easily...

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

        if (fobj->LOGFLAG)
        {
            flog = fopen(fobj->flogname, "a");
            if (flog == NULL)
            {
                printf("could not open %s to append!\n", fobj->flogname);
                printf("disabling logging...\n");
                fobj->LOGFLAG = 0;
            }
        }

		//create the afb/rfb - we don't want the time it takes to do this to
		//pollute the sieve timings		
		sprintf(syscmd, "%s -b %s -k -c 0 -F", jobs[i].sievername, filenames[i]);
		if (fobj->VFLAG > 0) printf("\ntest: generating factor bases\n");
		gettimeofday(&start, NULL);
		sysreturn = system(syscmd);
		gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);

		if (fobj->VFLAG > 0) printf("test: fb generation took %6.4f seconds\n", t_time);
		logprint(flog, "test: fb generation took %6.4f seconds\n", t_time);
		MySleep(.1);

		//start the test
		sprintf(syscmd,"%s%s -%c %s -f %u -c %u -o %s.out",
			jobs[i].sievername, fobj->VFLAG>0?" -v":"", side[0], filenames[i], jobs[i].startq, spq_range, filenames[i]);

		if (fobj->VFLAG > 0) printf("test: commencing test sieving of polynomial %d on the %s side over range %u-%u\n", i, 
			side, jobs[i].startq, jobs[i].startq + spq_range);
		logprint(flog, "test: commencing test sieving of polynomial %d on the %s side over range %u-%u\n", i, 
			side, jobs[i].startq, jobs[i].startq + spq_range);
		
        if (fobj->LOGFLAG)
        {
            print_job(&jobs[i], flog);
            fclose(flog);
        }

		gettimeofday(&start, NULL);
        sysreturn = system(syscmd);
		gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);
		
		//count relations
		sprintf(tmpbuf, "%s.out", filenames[i]);
		in = fopen(tmpbuf, "r");
		actual_range = 0;
		count = 0;
		if( !in )
		{
			score[i] = 999999999.;
			
			//est = 7*365*24*3600; // 7 years seems like a nice round number
		}
		else
		{
			// scan the data file and
			// 1) count the relations
			// 2) save the last four relations to a buffer in order to extract the last processed
			//		special-q.
			// we need both 1) and 2) to compute yield correctly.

			char **lines, *ptr, tmp[GSTR_MAXSIZE];
			int line;
			int j;

			lines = (char **)malloc(4 * sizeof(char *));
			for (j=0; j < 4; j++)
				lines[j] = (char *)malloc(GSTR_MAXSIZE * sizeof(char));

			line = 0;
			count = 0;
			while (1)
			{
				// read a line into the next position of the circular buffer
				ptr = fgets(tmp, GSTR_MAXSIZE, in);
				if (ptr == NULL) 
					break;

				// quick check that it might be a valid line
				if (strlen(tmp) > 30)
				{
					// wrap
					if (++line > 3) line = 0;
					// then copy
					strcpy(lines[line], tmp);
				}

				count++;
			}
			fclose(in);

			line = get_spq(lines, line, fobj);
			actual_range = line - jobs[i].startq;
			if (fobj->VFLAG > 0)
				printf("test: found %u relations in a range of %u special-q\n", 
				count, actual_range);

			if (actual_range > spq_range) 
				actual_range = spq_range;

			for (j=0; j < 4; j++)
				free(lines[j]);
			free(lines);

			score[i] = t_time / count;

			// use estimated sieving time to rank, not sec/rel, since the latter
			// is a function of parameterization and therefore not directly comparable
			// to each other.
			score[i] = (score[i] * jobs[i].min_rels * 1.25) / fobj->THREADS;
			// be conservative about estimates
		}

        if (fobj->LOGFLAG)
        {
            flog = fopen(fobj->flogname, "a");
        }

		if( score[i] < min_score )
		{
			minscore_id = i;
			min_score = score[i];
			if (fobj->VFLAG > 0) printf("test: new best estimated total sieving time = %s (with %d threads)\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);
			logprint(flog, "test: new best estimated total sieving time = %s (with %d threads)\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);

			// edit lbpr/a depending on test results.  we target something around 2 rels/Q.
			// could also change siever version in more extreme cases.
			if (count > 4*actual_range)
			{
				if (fobj->VFLAG > 0)
					printf("test: yield greater than 4x/spq, reducing lpbr/lpba\n");
				jobs[i].lpba--;
				jobs[i].lpbr--;
				jobs[i].mfba -= 2;
				jobs[i].mfbr -= 2;
			}

			if (count > 8*actual_range)
			{
				char *pos;
				int siever;

				pos = strstr(jobs[i].sievername, "gnfs-lasieve4I");
				siever = (pos[14] - 48) * 10 + (pos[15] - 48);

				if (fobj->VFLAG > 0)
					printf("test: yield greater than 8x/spq, reducing siever version\n");

				switch (siever)
				{
				case 11:
					if (fobj->VFLAG > 0) printf("test: siever version cannot be decreased further\n");
					jobs[i].snfs->siever = 11;
					break;

				case 12:
					pos[15] = '1';
					jobs[i].snfs->siever = 11;
					break;

				case 13:
					pos[15] = '2';
					jobs[i].snfs->siever = 12;
					break;

				case 14:
					pos[15] = '3';
					jobs[i].snfs->siever = 13;
					break;

				case 15:
					pos[15] = '4';
					jobs[i].snfs->siever = 14;
					break;

				case 16:
					pos[15] = '5';
					jobs[i].snfs->siever = 15;
					break;
				}
			}

			if (count < actual_range)
			{
				if (fobj->VFLAG > 0)
					printf("test: yield less than 1x/spq, increasing lpbr/lpba\n");
				
				jobs[i].lpba++;
				jobs[i].lpbr++;
				jobs[i].mfba += 2;
				jobs[i].mfbr += 2;
			}

			if (count < (actual_range/2))
			{
				char *pos;
				int siever;

				pos = strstr(jobs[i].sievername, "gnfs-lasieve4I");
				siever = (pos[14] - 48) * 10 + (pos[15] - 48);

				if (fobj->VFLAG > 0)
					printf("test: yield less than 1x/2*spq, increasing siever version\n");

				switch (siever)
				{
				case 16:
					if (fobj->VFLAG > 0) printf("test: siever version cannot be increased further\n");
					jobs[i].snfs->siever = 16;
					break;

				case 15:
					pos[15] = '6';
					jobs[i].snfs->siever = 16;
					break;

				case 14:
					pos[15] = '5';
					jobs[i].snfs->siever = 15;
					break;

				case 13:
					pos[15] = '4';
					jobs[i].snfs->siever = 14;
					break;

				case 12:
					pos[15] = '3';
					jobs[i].snfs->siever = 13;
					break;

				case 11:
					pos[15] = '2';
					jobs[i].snfs->siever = 12;
					break;
				}
			}
     	}
		else
		{
			if (fobj->VFLAG > 0) printf("test: estimated total sieving time = %s (with %d threads)\n\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);
			logprint(flog, "test: estimated total sieving time = %s (with %d threads)\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);
		}

        if (fobj->LOGFLAG)
        {
            fclose(flog);
        }

		remove(tmpbuf); // clean up after ourselves

        if (!are_files)
        {
            sprintf(tmpbuf, "%s", filenames[i]);
            remove(tmpbuf);
        }

		sprintf(tmpbuf, "%s.afb.0", filenames[i]);
		remove(tmpbuf);
	}

	// clean up memory allocated
	if( are_files )
	{
        // TODO: need to write parameter adjustments to file
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

    if (fobj->LOGFLAG)
    {
        flog = fopen(fobj->flogname, "a");
        gettimeofday(&stop2, NULL);
        t_time = ytools_difftime(&start2, &stop2);

        if (fobj->VFLAG > 0) printf("test: test sieving took %1.2f seconds\n", t_time);
        logprint(flog, "test: test sieving took %1.2f seconds\n", t_time);
    }

	return minscore_id;
}

void do_sieving(fact_obj_t *fobj, nfs_job_t *job)
{
	nfs_threaddata_t *thread_data;		//an array of thread data objects
	int i;
	FILE *fid;
	FILE *logfile;

	thread_data = (nfs_threaddata_t *)malloc(fobj->THREADS * sizeof(nfs_threaddata_t));
	for (i = 0; i < fobj->THREADS; i++)
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
        // if the user specified a range of q to sieve, then divvy up that
        // range over the threads.  otherwise divvy up the qrange that was
        // determined to be good for this input.
        // todo: in the user-specified case the range may be large which
        // can cause problems with abort-restarts because only a small fraction
        // of each range actually gets done.  better would be to use the qrange
        // that was determined to be good for this input and run multiple 
        // such ranges until the user-range is met.
		if (fobj->nfs_obj.rangeq > 0)
			thread_data[i].job.qrange = ceil((double)fobj->nfs_obj.rangeq / (double)fobj->THREADS);
		else
			thread_data[i].job.qrange = ceil((double)job->qrange / (double)fobj->THREADS);
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
	for (i = 0; i < fobj->THREADS - 1; i++)
		nfs_start_worker_thread(thread_data + i, 0);

	nfs_start_worker_thread(thread_data + i, 1);			

	// load threads with work
	for (i = 0; i < fobj->THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;
		t->job.startq = job->startq;
		job->startq += t->job.qrange;
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
            logprint(logfile, "nfs: commencing lattice sieving with %d threads\n", fobj->THREADS);
            fclose(logfile);
        }
    }

	// create a new lasieve process in each thread and watch it
	for (i = 0; i < fobj->THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;

		if (i == fobj->THREADS - 1) {
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

	for (i = 0; i < fobj->THREADS; i++) {
		nfs_threaddata_t *t = thread_data + i;

		if (i < fobj->THREADS - 1) {
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
	for (i = 0; i < fobj->THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;
		savefile_concat(t->outfilename,fobj->nfs_obj.outputfile,fobj->nfs_obj.mobj);

		// accumulate relation counts
		job->current_rels += thread_data[i].job.current_rels;

		remove(thread_data[i].outfilename);
	}

	if ((fid = fopen("rels.add", "r")) != NULL)
	{
		char tmpstr[1024];
		uint32_t count = 0;

		while (fgets(tmpstr, GSTR_MAXSIZE, fid) != NULL)
			count++;
		fclose(fid);

		if (fobj->VFLAG > 0) printf("nfs: adding %u rels from rels.add\n",count);

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
                logprint(logfile, "nfs: adding %u rels from rels.add\n", count);
                fclose(logfile);
            }
        }

		savefile_concat("rels.add",fobj->nfs_obj.outputfile,fobj->nfs_obj.mobj);
		remove("rels.add");
	}

	//stop worker threads
	for (i = 0; i < fobj->THREADS - 1; i++)
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
		
	//start ggnfs binary - new win64 ASM enabled binaries current have a problem with this:
	//sprintf(syscmd,"%s%s -%c %s -f %u -c %u -o %s -n %d",
	//		thread_data->job.sievername, fobj->VFLAG>0?" -v":"", *side,
	//		fobj->nfs_obj.job_infile, thread_data->job.startq, 
	//		thread_data->job.qrange, thread_data->outfilename, thread_data->tindex);

    // todo: add command line input of arbitrary argument string to append to this command
	// but not this:
	snprintf(syscmd, GSTR_MAXSIZE, "%s%s -f %u -c %u -o %s -n %d -%c %s ",
			thread_data->job.sievername, fobj->VFLAG>0?" -v":"", thread_data->job.startq, 
			thread_data->job.qrange, thread_data->outfilename, thread_data->tindex,
			*side, fobj->nfs_obj.job_infile);

	if (fobj->VFLAG >= 0)
	{
		printf("nfs: commencing %s side lattice sieving over range: %u - %u\n",
			side, thread_data->job.startq, thread_data->job.startq + thread_data->job.qrange);
	}
	if (fobj->VFLAG > 1) printf("syscmd: %s\n", syscmd);
	if (fobj->VFLAG > 1) fflush(stdout);
	cmdret = system(syscmd);

	// a ctrl-c abort signal is caught by the system command, and nfsexit never gets called.
	// so check for abnormal exit from the system command.
	// -1073741819 is apparently what ggnfs returns when it crashes, which
	// we don't want to interpret as an abort.
    if (!((cmdret == 0) || cmdret == -1073741819))
    {
        if (NFS_ABORT < 1)
        {
            NFS_ABORT = 1;
        }
    }

	// count the relations produced
	MySleep(100);
	fid = fopen(thread_data->outfilename,"r");
	if (fid != NULL)
	{
		thread_data->job.current_rels = 0;
		while (fgets(tmpstr, GSTR_MAXSIZE, fid) != NULL)
			thread_data->job.current_rels++;
		fclose(fid);
	}
	else
	{
		printf("nfs: could not open output file, possibly bad path to siever\n");
	}

	return 0;
}

#endif
