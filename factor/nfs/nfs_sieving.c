#include "nfs.h"

void do_sieving(fact_obj_t *fobj, ggnfs_job_t *job)
{
	nfs_threaddata_t *thread_data;		//an array of thread data objects
	int i;
	FILE *logfile;

	thread_data = (nfs_threaddata_t *)malloc(THREADS * sizeof(nfs_threaddata_t));
	for (i=0; i<THREADS; i++)
	{
		sprintf(thread_data[i].outfilename, "rels%d.dat", i);
		thread_data[i].job.fblim = job->fblim;
		thread_data[i].job.lambda = job->lambda;
		thread_data[i].job.lpb = job->lpb;
		thread_data[i].job.mfb = job->mfb;
		if (fobj->nfs_obj.rangeq > 0)
			thread_data[i].job.qrange = ceil((double)fobj->nfs_obj.rangeq / (double)THREADS);
		else
			thread_data[i].job.qrange = job->qrange;
		thread_data[i].job.min_rels = job->min_rels;
		thread_data[i].job.current_rels = job->current_rels;
		thread_data[i].job.siever = job->siever;
		thread_data[i].job.startq = job->startq;
		strcpy(thread_data[i].job.sievername, job->sievername);

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
		printf("could not open yafu logfile for appending\n");
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
//void process_hypercube(static_conf_t *sconf,dynamic_conf_t *dconf)
{
	//top level sieving function which performs all work for a single
	//new a coefficient.  has pthread calling conventions, meant to be
	//used in a multi-threaded environment
	nfs_threaddata_t *thread_data = (nfs_threaddata_t *)ptr;
	fact_obj_t *fobj = thread_data->fobj;
	char syscmd[GSTR_MAXSIZE], tmpstr[GSTR_MAXSIZE];	
	FILE *fid;
	int cmdret;

	//remove any temporary relation files
	remove(thread_data->outfilename);
		
	//start ggnfs binary
	if (fobj->nfs_obj.sq_side)
		sprintf(syscmd,"%s -a %s -f %u -c %u -o %s",
			thread_data->job.sievername, fobj->nfs_obj.job_infile, thread_data->job.startq, 
			thread_data->job.qrange, thread_data->outfilename);
	else
		sprintf(syscmd,"%s -r %s -f %u -c %u -o %s",
			thread_data->job.sievername, fobj->nfs_obj.job_infile, thread_data->job.startq, 
			thread_data->job.qrange, thread_data->outfilename);

	if (VFLAG >= 0)
	{
		if (fobj->nfs_obj.sq_side)
			printf("nfs: commencing algebraic side lattice sieving over range: %u - %u\n",
				thread_data->job.startq, thread_data->job.startq + thread_data->job.qrange);
		else
			printf("nfs: commencing rational side lattice sieving over range: %u - %u\n",
				thread_data->job.startq, thread_data->job.startq + thread_data->job.qrange);
	}
	cmdret = system(syscmd);

	// count the relations produced
	MySleep(100);
	fid = fopen(thread_data->outfilename,"r");
	thread_data->job.current_rels = 0;
	while (fgets(tmpstr, GSTR_MAXSIZE, fid) != NULL)
		thread_data->job.current_rels++;
	fclose(fid);

	return 0;
}

