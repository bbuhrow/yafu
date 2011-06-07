#include "gnfs.h"
#include "yafu_string.h"
#include "arith.h"
#include "factor.h"
#include "qs.h"

//----------------------- LOCAL DATA TYPES -----------------------------------//

/* used to place a deadline on how long polynomial 
   selection will run. Note that the time budget is
   independent of CPU speed; faster CPUs will simply
   search more of the polynomial space */

typedef struct {
	uint32 bits;
	uint32 seconds;
} poly_deadline_t;

static const poly_deadline_t time_limits[] = {
	{264, 4 * 60},
	{304, 8 * 60},
	{320, 15 * 60},
	{348, 30 * 60},
	{365, 1 * 3600},
	{383, 2 * 3600},
	{399, 4 * 3600},
	{416, 8 * 3600},
	{433, 16 * 3600},
	{449, 32 * 3600},
	{466, 64 * 3600},
	{482, 100 * 3600},
	{498, 200 * 3600},
	{514, 300 * 3600},
};

#define NUM_TIME_LIMITS sizeof(time_limits)/sizeof(time_limits[0])


enum nfs_thread_command {
	NFS_COMMAND_INIT,
	NFS_COMMAND_WAIT,
	NFS_COMMAND_RUN,
	NFS_COMMAND_RUN_POLY,
	NFS_COMMAND_END
};

typedef struct
{
	uint32 fblim;
	uint32 lpb;
	uint32 mfb;
	float lambda;
	uint32 siever; 
	uint32 qrange;
	char sievername[1024];
	uint32 startq;
	uint32 min_rels;
	uint32 current_rels;
} ggnfs_job_t;

typedef struct {
	// stuff for parallel ggnfs sieving
	char outfilename[80];
	ggnfs_job_t job;

	// stuff for parallel msieve poly select
	char *polyfilename, *logfilename;
	uint64 poly_lower;
	uint64 poly_upper;
	msieve_obj *obj;
	mp_t *mpN;
	factor_list_t *factor_list;

	/* fields for thread pool synchronization */
	volatile enum nfs_thread_command command;

#if defined(WIN32) || defined(_WIN64)
	HANDLE thread_id;
	HANDLE run_event;
	HANDLE finish_event;
#else
	pthread_t thread_id;
	pthread_mutex_t run_lock;
	pthread_cond_t run_cond;
#endif

} nfs_threaddata_t;

//----------------------- LOCAL FUNCTIONS -------------------------------------//
void *lasieve_launcher(void *ptr);
void *polyfind_launcher(void *ptr);
void find_best_msieve_poly(z *N, ggnfs_job_t *job);
void msieve_to_ggnfs(ggnfs_job_t *job);
void ggnfs_to_msieve(ggnfs_job_t *job);
void get_ggnfs_params(z *N, ggnfs_job_t *job);
int check_existing_files(z *N, uint32 *last_spq, ggnfs_job_t *job);
void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj);
uint32 get_spq(char *line);
uint32 do_msieve_filtering(msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN);
void do_msieve_polyselect(msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN, factor_list_t *factor_list);
void get_polysearch_params(uint64 *start, uint64 *range);
void do_sieving(ggnfs_job_t *job);
void win_file_concat(char *filein, char *fileout);
void nfs_stop_worker_thread(nfs_threaddata_t *t,
				uint32 is_master_thread);
void nfs_start_worker_thread(nfs_threaddata_t *t, 
				uint32 is_master_thread);
void nfsexit(int sig);
#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI nfs_worker_thread_main(LPVOID thread_data);
#else
void *nfs_worker_thread_main(void *thread_data);
#endif

int NFS_ABORT;

#define USE_NFS
#ifdef USE_NFS

void nfsexit(int sig)
{
	printf("\nAborting... \n");
	NFS_ABORT = 1;
	return;
}


//----------------------- NFS ENTRY POINT ------------------------------------//
void test_msieve_gnfs(fact_obj_t *fobj)
{
	z *N = &fobj->qs_obj.n;
	char *input;
	str_t input_str;
	msieve_obj *obj = NULL;
	uint32 max_relations = 0;
	uint32 seed1 = g_rand.low;
	uint32 seed2 = g_rand.hi;
	uint64 nfs_lower = 0;
	uint64 nfs_upper = 0;
	enum cpu_type cpu = yafu_get_cpu_type();
	uint32 mem_mb = 0;
	uint32 which_gpu = 0;
	mp_t mpN;
	factor_list_t factor_list;
	uint32 flags = 0;
	ggnfs_job_t job;
	uint32 relations_needed = 1;	
	uint32 startq, qrange;
	int is_continuation;
	uint32 last_specialq = 0;
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	TIME_DIFF *	difference;
	double t_time;
	FILE *logfile;
	int statenum;
	char tmpstr[GSTR_MAXSIZE];	
	int process_done;

	//below a certain amount, revert to SIQS
	if (ndigits(N) < MIN_NFS_DIGITS)
	{
		SIQS(fobj);
		return;
	}
		
	//first a small amount of trial division
	//which will add anything found to the global factor list
	//this is only called with the main thread
	if (VFLAG > 0)
		printf("nfs: commencing trial factoring\n");
	zCopy(N,&fobj->div_obj.n);	
	zTrial(10000,0,fobj);
	zCopy(&fobj->div_obj.n,N);

	if (isPrime(N))
	{
		add_to_factor_list(fobj, N);
		zCopy(&zOne,N);

		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			if (VFLAG >= 0)
				printf("PRP%d = %s\n",ndigits(N),z2decstr(N,&gstr1));
			logprint(logfile, "PRP%d = %s\n",ndigits(N),z2decstr(N,&gstr1));
			fclose(logfile);
		}		
		return;
	}

	//check for inconsistent options
	if (((GGNFS_STARTQ > 0) && !(GGNFS_RANGEQ > 0)) ||
		(!(GGNFS_STARTQ > 0) && (GGNFS_RANGEQ > 0)))
	{
		printf("-f and -c must be specified together\n");
		exit(1);
	}

	if (((GGNFS_STARTQ > 0) || (GGNFS_RANGEQ > 0)) && GGNFS_POLY_ONLY)
	{
		printf("bad options: -np with -f or -c\n");
		exit(1);
	}

	if (((GGNFS_STARTQ > 0) || (GGNFS_RANGEQ > 0)) && GGNFS_POST_ONLY)
	{
		printf("bad options: -nc with -f or -c\n");
		exit(1);
	}

	if (GGNFS_SIEVE_ONLY && GGNFS_POST_ONLY)
	{
		printf("bad options: -nc with -ns\n");
		exit(1);
	}

	if (GGNFS_SIEVE_ONLY && GGNFS_POLY_ONLY)
	{
		printf("bad options: -np with -ns\n");
		exit(1);
	}

	if (GGNFS_POLY_ONLY && GGNFS_POST_ONLY)
	{
		printf("bad options: -nc with -np\n");
		exit(1);
	}

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	NFS_ABORT = 0;
	signal(SIGINT,nfsexit);

	//start a counter for the whole job
	gettimeofday(&start, NULL);
	
	//nfs state machine:
	statenum = 0;
	process_done = 0;
	while (!process_done)
	{
		switch (statenum)
		{
		case 0: //"init":

			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				if (VFLAG >= 0)
					printf("nfs: commencing gnfs on c%d: %s\n",ndigits(N), z2decstr(N,&gstr1));
				logprint(logfile, "nfs: commencing gnfs on c%d: %s\n",ndigits(N), z2decstr(N,&gstr1));
				fclose(logfile);
			}

			//write the input bigint as a string
			sInit(&input_str);
			input = z2decstr(N, &input_str);	

			//create an msieve_obj
			obj = msieve_obj_new(input, flags, GGNFS_OUTPUTFILE, GGNFS_LOGFILE, GGNFS_FBFILE, seed1, seed2,
				max_relations, nfs_lower, nfs_upper, cpu, L1CACHE, L2CACHE, THREADS, mem_mb, which_gpu, 0.0);

			//convert input to msieve bigint notation and initialize a list of factors
			z2mp_t(N,&mpN);
			factor_list_init(&factor_list);

			//check if this is a continuation of a factorization (defined when
			//a .fb file already exists for this N, or both a .fb and msieve.dat
			//file exist for this N, or if a ggnfs.job file exists for this N.
			//if any of the .job, .fb and/or .dat exist but for conflicting N's,
			//complain and stop.  else, start new job.
			job.current_rels = 0;
			is_continuation = check_existing_files(N, &last_specialq, &job);

			//find best job parameters
			get_ggnfs_params(N,&job);						

			//determine sieving start value and range.
			if (GGNFS_RANGEQ > 0)
				qrange = ceil((double)GGNFS_RANGEQ / (double)THREADS);
			else
				qrange = job.qrange;

			if (is_continuation > 0)
			{		
				if (last_specialq == 0)
				{
					if (VFLAG >= 0)
						printf("nfs: continuing job - could not determine last special q; using default startq\n");

					logfile = fopen(flogname, "a");
					if (logfile == NULL)
						printf("could not open yafu logfile for appending\n");
					else
					{
						logprint(logfile, "nfs: continuing job - could not determine last special q; using default startq\n");
						fclose(logfile);
					}

					if (GGNFS_STARTQ > 0)
						startq = GGNFS_STARTQ;
					else
						startq = job.fblim / 2;

					// next step is sieving
					statenum = 2;
				}
				else
				{
					if (VFLAG >= 0)
						printf("nfs: found ~ %u relations, continuing job at specialq = %u\n",
						job.current_rels,last_specialq);

					logfile = fopen(flogname, "a");
					if (logfile == NULL)
						printf("could not open yafu logfile for appending\n");
					else
					{
						logprint(logfile, "nfs: found ~ %u relations, continuing job at specialq = %u\n",
							job.current_rels,last_specialq);
						fclose(logfile);
					}

					if (GGNFS_STARTQ > 0)
					{
						startq = GGNFS_STARTQ;
						statenum = 2;
					}
					else if (GGNFS_SIEVE_ONLY)
					{
						startq = last_specialq;
						statenum = 2;
					}
					else
					{
						startq = last_specialq;

						// since we apparently found some relations, see if it is enough
						statenum = 3;
					}
				}
			}
			else if (is_continuation == 0)
			{
				// new factorization				
				if ((GGNFS_STARTQ > 0) || GGNFS_SIEVE_ONLY)
				{
					printf("no job file exists for this input\n");
					exit(1);
				}
				startq = job.fblim / 2;
				statenum = 1;								
			}
			else
			{
				// failed restart
				process_done = 1;
				break;
			}

			//load the job structure with the starting Q.
			job.startq = startq;
			job.qrange = qrange;

			//override with custom commands, if applicable
			if (GGNFS_POLY_ONLY)
				statenum = 1;
			else if (GGNFS_SIEVE_ONLY)
				statenum = 2;
			else if (GGNFS_POST_ONLY)
				statenum = 3;

			break;

		case 1: //"polysearch":
			
			do_msieve_polyselect(obj, &job, &mpN, &factor_list);
			
			if (GGNFS_POLY_ONLY)
				process_done = 1;
			statenum = 2;
			break;

		case 2: //"sieve":

			do_sieving(&job);

			//see how we're doing
			if (GGNFS_SIEVE_ONLY || (GGNFS_RANGEQ > 0))
				process_done = 1;
			else
			{
				if (job.current_rels >= job.min_rels)
				{
					if (VFLAG > 0)
						printf("found %u relations, need at least %u, proceeding with filtering ...\n",
						job.current_rels, job.min_rels);
				
					statenum = 3;
				}
				else
				{
					if (VFLAG > 0)
						printf("found %u relations, need at least %u, continuing with sieving ...\n",
						job.current_rels, job.min_rels);

					statenum = 2;
				}
			}

			break;

		case 3: //"filter":

			relations_needed = do_msieve_filtering(obj, &job, &mpN);
			if (relations_needed == 0)
				statenum = 4;	//proceed to linear algebra
			else
			{
				if (GGNFS_POST_ONLY)
					process_done = 1;
				else
					statenum = 2;	//more sieving
			}
			break;

		case 4: //case "linalg":

			//msieve: build matrix
			flags = 0;
			flags = flags | MSIEVE_FLAG_USE_LOGFILE;
			if (VFLAG > 0)
				flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
			flags = flags | MSIEVE_FLAG_NFS_LA;
			obj->flags = flags;

			if (VFLAG >= 0)
				printf("nfs: commencing msieve linear algebra\n");
	
			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "nfs: commencing msieve linear algebra\n");
				fclose(logfile);
			}

			nfs_solve_linear_system(obj, &mpN);			
			statenum = 5;
			break;

		case 5: //case "sqrt":

			//msieve: find factors
			flags = 0;
			flags = flags | MSIEVE_FLAG_USE_LOGFILE;
			if (VFLAG > 0)
				flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
			flags = flags | MSIEVE_FLAG_NFS_SQRT;
			obj->flags = flags;

			if (VFLAG >= 0)
				printf("nfs: commencing msieve sqrt\n");

			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "nfs: commencing msieve sqrt\n");
				fclose(logfile);
			}
			nfs_find_factors(obj, &mpN, &factor_list);
			
			extract_factors(&factor_list,fobj);
			if (zCompare(N,&zOne) == 0)
				statenum = 6;		//completely factored, clean up everything
			else
				statenum = 7;		//not factored completely, keep files and stop
			break;

		case 6: //case "cleanup":
			
			remove(GGNFS_OUTPUTFILE);
			remove(GGNFS_FBFILE);
			sprintf(tmpstr, "%s.p",GGNFS_OUTPUTFILE);	remove(tmpstr);			
			sprintf(tmpstr, "%s.br",GGNFS_OUTPUTFILE);	remove(tmpstr);
			sprintf(tmpstr, "%s.cyc",GGNFS_OUTPUTFILE);	remove(tmpstr);
			sprintf(tmpstr, "%s.dep",GGNFS_OUTPUTFILE);	remove(tmpstr);
			sprintf(tmpstr, "%s.hc",GGNFS_OUTPUTFILE);	remove(tmpstr);
			sprintf(tmpstr, "%s.mat",GGNFS_OUTPUTFILE);	remove(tmpstr);	
			sprintf(tmpstr, "%s.lp",GGNFS_OUTPUTFILE);	remove(tmpstr);
			sprintf(tmpstr, "%s.d",GGNFS_OUTPUTFILE);	remove(tmpstr);

			gettimeofday(&stop, NULL);

			difference = my_difftime (&start, &stop);

			t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
			free(difference);	

			if (VFLAG >= 0)
				printf("NFS elapsed time = %6.4f seconds.\n",t_time);

			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "NFS elapsed time = %6.4f seconds.\n",t_time);
				fclose(logfile);
			}

			// free stuff
			msieve_obj_free(obj);
			sFree(&input_str);
			statenum = 7;
			break;

		case 7:
			process_done = 1;
			break;

		default:
			printf("unknown state, bailing\n");
			exit(1);

		}

		//after every state, check elasped time against a specified timeout value
		gettimeofday(&stop, NULL);
		difference = my_difftime (&start, &stop);

		t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);	
		if ((GGNFS_TIMEOUT > 0) && (t_time > (double)GGNFS_TIMEOUT))
		{
			if (VFLAG >= 0)
				printf("NFS timeout after %6.4f seconds.\n",t_time);

			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "NFS timeout after %6.4f seconds.\n",t_time);
				fclose(logfile);
			}
			process_done = 1;
		}

		if (NFS_ABORT)
			process_done = 1;
	}
	
	return;
}

void do_sieving(ggnfs_job_t *job)
{
	nfs_threaddata_t *thread_data;		//an array of thread data objects
	int i;
	//uint32 startq = job.startq;
	FILE *logfile;
	char syscmd[1024];

	thread_data = (nfs_threaddata_t *)malloc(THREADS * sizeof(nfs_threaddata_t));
	for (i=0; i<THREADS; i++)
	{
		sprintf(thread_data[i].outfilename, "rels%d.dat", i);
		thread_data[i].job.fblim = job->fblim;
		thread_data[i].job.lambda = job->lambda;
		thread_data[i].job.lpb = job->lpb;
		thread_data[i].job.mfb = job->mfb;
		if (GGNFS_RANGEQ > 0)
			thread_data[i].job.qrange = ceil((double)GGNFS_RANGEQ / (double)THREADS);
		else
			thread_data[i].job.qrange = job->qrange;
		thread_data[i].job.min_rels = job->min_rels;
		thread_data[i].job.current_rels = job->current_rels;
		thread_data[i].job.siever = job->siever;
		thread_data[i].job.startq = job->startq;
		strcpy(thread_data[i].job.sievername, job->sievername);
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

	logfile = fopen(flogname, "a");
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
		
#if defined(WIN32)
		int a;

		// test for cat
		sprintf(syscmd,"cat %s >> %s",t->outfilename,GGNFS_OUTPUTFILE);
		a = system(syscmd);
	
		if (a)		
			win_file_concat(t->outfilename,GGNFS_OUTPUTFILE);

#else
		sprintf(syscmd,"cat %s >> %s",t->outfilename,GGNFS_OUTPUTFILE);
		system(syscmd);
#endif
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

void win_file_concat(char *filein, char *fileout)
{
	FILE *in, *out;
	//printf("for optimal performance, consider installing unix utilities for windows:\n");
	//printf("http://unxutils.sourceforge.net/ \n");

	in = fopen(filein,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading\n",filein);
		exit(-1);
	}

	out = fopen(fileout,"a");
	if (out == NULL)
	{
		printf("could not open %s for appending\n",fileout);
		exit(-1);
	}

	while (!feof(in))
	{
		char tmpline[GSTR_MAXSIZE], *tmpptr;
		tmpptr = fgets(tmpline, GSTR_MAXSIZE, in);					
		if (tmpptr == NULL)
			break;
		else
			fputs(tmpline, out);
	}
	fclose(in);
	fclose(out);

	return;
}

void do_msieve_polyselect(msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN, factor_list_t *factor_list)
{
	FILE *logfile;
	uint32 flags;
	z N;

	nfs_threaddata_t *thread_data;		//an array of thread data objects
	int i,j;
	char syscmd[1024];
	char master_polyfile[80];
	uint64 start = 0, range = 0;
	uint32 deadline;
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	TIME_DIFF *	difference;
	double t_time;

	//file into which we will combine all of the thread results
	sprintf(master_polyfile,"%s.p",GGNFS_OUTPUTFILE);

	//make sure we are starting from scratch
	remove(master_polyfile);

	/* figure out how long poly selection should take */
	zInit(&N);
	mp_t2z(mpN,&N);
	i = zBits(&N);
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

	if (GGNFS_POLY_OPTION == 0)
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

	//set flags to do poly selection
	flags = 0;
	if (VFLAG > 0)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_POLY1;
	flags = flags | MSIEVE_FLAG_NFS_POLY2;		
	
	for (i=0; i<THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;		

		t->logfilename = (char *)malloc(80 * sizeof(char));
		t->polyfilename = (char *)malloc(80 * sizeof(char));
		sprintf(t->polyfilename,"%s.%d",GGNFS_OUTPUTFILE,i);
		sprintf(t->logfilename,"%s.%d",GGNFS_LOGFILE,i);

		//create an msieve_obj.  for poly select, the intermediate output file should be specified in
		//the savefile field		
		t->obj = msieve_obj_new(obj->input, flags, t->polyfilename, t->logfilename, GGNFS_FBFILE, 
			g_rand.low, g_rand.hi, (uint32)0, (uint64)1, (uint64)1001, 
			obj->cpu, (uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);

		//pointers to things that are static during poly select
		t->mpN = mpN;
		t->factor_list = factor_list;

	}

	// log that we are starting
	logfile = fopen(flogname, "a");
	if (logfile == NULL)
		printf("could not open yafu logfile for appending\n");
	else
	{
		logprint(logfile, "nfs: commencing poly selection with %d threads\n",THREADS);
		logprint(logfile, "nfs: setting deadline of %u seconds\n",deadline);
		fclose(logfile);
	}

	/* activate the threads one at a time. The last is the
				master thread (i.e. not a thread at all). */		
	for (i = 0; i < THREADS - 1; i++)
		nfs_start_worker_thread(thread_data + i, 0);

	nfs_start_worker_thread(thread_data + i, 1);			

	// determine the start and range values
	get_polysearch_params(&start, &range);	

	while (1)
	{
		//assign some work
		for (i=0; i<THREADS; i++)
		{
			nfs_threaddata_t *t = thread_data + i;		

			//create temporary storage files
			t->logfilename = (char *)malloc(80 * sizeof(char));
			t->polyfilename = (char *)malloc(80 * sizeof(char));
			sprintf(t->polyfilename,"%s.%d",GGNFS_OUTPUTFILE,i);
			sprintf(t->logfilename,"%s.%d",GGNFS_LOGFILE,i);

			//create an msieve_obj.  for poly select, the intermediate output file should be specified in
			//the savefile field.  it may seem inefficient to continually allocate msieve_obj and then
			//destroy them within the timed loop, but this works.  the interface is pretty finicky - most
			//other attempts failed.
			t->obj = msieve_obj_new(obj->input, flags, t->polyfilename, t->logfilename, GGNFS_FBFILE, 
				g_rand.low, g_rand.hi, (uint32)0, (uint64)start, (uint64)(start + range), 
				obj->cpu, (uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);

			//pointers to things that are static during poly select
			t->mpN = mpN;
			t->factor_list = factor_list;
			
			if (GGNFS_POLY_OPTION == 2)
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
		}

		// launch a new poly selection process in each thread
		for (i = 0; i < THREADS; i++) 
		{
			nfs_threaddata_t *t = thread_data + i;

			if (i == THREADS - 1) {		
				polyfind_launcher(t);
			}
			else {
				t->command = NFS_COMMAND_RUN_POLY;
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
		
			// concatenate all of the .p files produced into a master .p file
#if defined(WIN32)
			int a;

			// test for cat
			sprintf(syscmd,"cat %s.p >> %s",t->polyfilename,master_polyfile);
			//printf("%s\n",syscmd);
			a = system(syscmd);
	
			if (a)		
			{
				char tmp[80];
				sprintf(tmp, "%s.p",t->polyfilename);
				win_file_concat(tmp, master_polyfile);
			}

#else
			sprintf(syscmd,"cat %s.p >> %s",t->polyfilename,master_polyfile);
			system(syscmd);
#endif

			// remove each thread's .p file after copied
			//printf("removing %s\n",thread_data[i].polyfilename);
			sprintf(syscmd, "%s.p",thread_data[i].polyfilename); 
			remove(syscmd);	

			// also remove the temporary log file
			sprintf(syscmd,"%s.%d",GGNFS_LOGFILE,i);
			remove(syscmd);
		}

		//free data
		for (i=0; i<THREADS - 1; i++)
		{
			nfs_threaddata_t *t = thread_data + i;

			free(t->logfilename);
			free(t->polyfilename);
			msieve_obj_free(t->obj);
		}
		msieve_obj_free((thread_data + i)->obj);

		remove(GGNFS_FBFILE);

		//check the time
		gettimeofday(&stopt, NULL);
		difference = my_difftime (&startt, &stopt);

		t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);	

		if ((uint32)t_time >= deadline)
		{
			printf("deadline of %u seconds met; poly select done\n",deadline);
			break;
		}

		if (GGNFS_POLYRANGE > 0)
		{
			printf("custom range search complete in %6.4f seconds\n",t_time);
			break;
		}

		if (GGNFS_POLY_OPTION == 2)
		{
			// 'deep' searching assigns all threads to the same range
			start += range;
		}

		if (NFS_ABORT)
			break;
	}			

	//stop worker threads
	for (i=0; i<THREADS - 1; i++)
	{
		nfs_threaddata_t *t = thread_data + i;
		nfs_stop_worker_thread(t, 0);		
	}
	nfs_stop_worker_thread(thread_data + i, 1);

	//free the thread structure
	free(thread_data);		

	// parse the master .p file and find the best poly	
	find_best_msieve_poly(&N, job);
	//also create a .fb file
	ggnfs_to_msieve(job);
	zFree(&N);

	return;
}

void get_polysearch_params(uint64 *start, uint64 *range)
{

	//search smallish chunks of the space in parallel until we've hit our deadline
	if (GGNFS_POLYSTART > 0)
		*start = GGNFS_POLYSTART;
	else
		*start = 1ULL;
	
	if (GGNFS_POLYRANGE > 0)
	{
		// search custom range of leading coefficient.  This effectively ignores the deadline
		// because the deadline is only checked after each range is done.  if we assign the
		// entire range up front, the deadline is never checked before polysearch finishes.
		*range = GGNFS_POLYRANGE / THREADS;

		// sanity check
		if (*range == 0) *range = 1;
	}
	else
	{
		// loop on a default range until the time deadline is reached
		*range = (uint64)GGNFS_POLYBATCH;
	}

	return;
}

uint32 do_msieve_filtering(msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN)
{
	FILE *tmp, *logfile;
	uint32 relations_needed;
	uint32 flags = 0;

	flags = flags | MSIEVE_FLAG_USE_LOGFILE;
	if (VFLAG > 0)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_FILTER;
	obj->flags = flags;

	//msieve: filter
	if (VFLAG >= 0)
		printf("nfs: commencing msieve filtering\n");

	logfile = fopen(flogname, "a");
	if (logfile == NULL)
		printf("could not open yafu logfile for appending\n");
	else
	{
		logprint(logfile, "nfs: commencing msieve filtering\n");
		fclose(logfile);
	}

	tmp = fopen(GGNFS_FBFILE,"r");
	if (tmp == NULL)
		ggnfs_to_msieve(job);
	else
		fclose(tmp);

	printf("%s\n",mp_print(mpN, 10, NULL, gstr1.s));
	relations_needed = nfs_filter_relations(obj, mpN);

	return relations_needed;
}

void nfs_start_worker_thread(nfs_threaddata_t *t, 
				uint32 is_master_thread) {

	/* create a thread that will process a polynomial The last poly does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread) {
		//matrix_thread_init(t);
		return;
	}

	t->command = NFS_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
	t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
	t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	t->thread_id = CreateThread(NULL, 0, nfs_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_create(&t->thread_id, NULL, nfs_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != NFS_COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
}

void nfs_stop_worker_thread(nfs_threaddata_t *t,
				uint32 is_master_thread)
{
	if (is_master_thread) {
		//matrix_thread_free(t);
		return;
	}

	t->command = NFS_COMMAND_END;
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
DWORD WINAPI nfs_worker_thread_main(LPVOID thread_data) {
#else
void *nfs_worker_thread_main(void *thread_data) {
#endif
	nfs_threaddata_t *t = (nfs_threaddata_t *)thread_data;

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);		
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == NFS_COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */

		if (t->command == NFS_COMMAND_RUN)
			lasieve_launcher(t);
		else if (t->command == NFS_COMMAND_RUN_POLY)
			polyfind_launcher(t);
		else if (t->command == NFS_COMMAND_END)
			break;

		/* signal completion */

		t->command = NFS_COMMAND_WAIT;
#if defined(WIN32) || defined(_WIN64)
		SetEvent(t->finish_event);		
#else
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

#if defined(WIN32) || defined(_WIN64)
	return 0;
#else
	return NULL;
#endif
}

void *lasieve_launcher(void *ptr)
//void process_hypercube(static_conf_t *sconf,dynamic_conf_t *dconf)
{
	//top level sieving function which performs all work for a single
	//new a coefficient.  has pthread calling conventions, meant to be
	//used in a multi-threaded environment
	nfs_threaddata_t *thread_data = (nfs_threaddata_t *)ptr;
	char syscmd[GSTR_MAXSIZE], tmpstr[GSTR_MAXSIZE];	
	FILE *fid;
	int cmdret;

	//remove any temporary relation files
	remove(thread_data->outfilename);
		
	//start ggnfs binary
	if (GGNFS_SQ_SIDE)
		sprintf(syscmd,"%s -a %s -f %u -c %u -o %s",
			thread_data->job.sievername, GGNFS_JOB_INFILE, thread_data->job.startq, 
			thread_data->job.qrange, thread_data->outfilename);
	else
		sprintf(syscmd,"%s -r %s -f %u -c %u -o %s",
			thread_data->job.sievername, GGNFS_JOB_INFILE, thread_data->job.startq, 
			thread_data->job.qrange, thread_data->outfilename);

	if (VFLAG >= 0)
	{
		if (GGNFS_SQ_SIDE)
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
			//thread_data->poly_lower, thread_data->poly_upper);
			t->obj->nfs_lower, t->obj->nfs_upper);
		fflush(stdout);
	}

	//start polyfind
	//factor_integer(t->obj->input, t->obj->flags, GGNFS_OUTPUTFILE, GGNFS_LOGFILE, GGNFS_FBFILE, 
	//	&g_rand.low, &g_rand.hi, 0, t->obj->nfs_lower, t->obj->nfs_upper, t->obj->cpu, 
	//	L1CACHE, L2CACHE, THREADS, 0, 0, 0.0);
	factor_gnfs(t->obj, t->mpN, t->factor_list);

	return 0;
}


void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj)
{
	z tmp1, tmp2, *N = &fobj->qs_obj.n;
	int i;
	FILE *logfile;
	char c[3];

	// extract the factors
	zInit(&tmp1);
	zInit(&tmp2);
	zCopy(N,&tmp1);
	for (i=0;i<factor_list->num_factors;i++)
	{
		z tmpz;
		zInit(&tmpz);				
		
		mp_t2z(&factor_list->final_factors[i]->factor,&tmpz);

		zDiv(&tmp1,&tmpz,N,&tmp2);
		if (isPrime(&tmpz))
		{
			tmpz.type = PRP;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"prp");
		}
		else
		{
			tmpz.type = COMPOSITE;
			add_to_factor_list(fobj, &tmpz);
			strcpy(c,"C");
		}

		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "%s%d = %s\n",c,ndigits(&tmpz),z2decstr(&tmpz,&gstr1));
			fclose(logfile);
		}		

		zFree(&tmpz);		
		zCopy(N,&tmp1);
	}

	if (zCompare(&tmp1,&zOne) > 0)
	{
		char c[3];
		if (isPrime(&tmp1))
		{
			tmp1.type = PRP;
			add_to_factor_list(fobj, &tmp1);
			strcpy(c,"prp");
			zCopy(&zOne, N);
		}
		else
		{
			tmp1.type = COMPOSITE;
			add_to_factor_list(fobj, &tmp1);
			strcpy(c,"C");
			zCopy(&zOne, N);
		}

		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "%s%d = %s\n",c,ndigits(&tmp1),z2decstr(&tmp1,&gstr1));
			fclose(logfile);
		}		
	}

	zFree(&tmp1);
	zFree(&tmp2);

	return;
}


int check_existing_files(z *N, uint32 *last_spq, ggnfs_job_t *job)
{
	//check to see if there is a .job for this number
	FILE *in, *logfile;
	char line[GSTR_MAXSIZE], *ptr;
	int ans;
	z tmpz;


	in = fopen(GGNFS_JOB_INFILE,"r");
	if (in == NULL)
		return 0;	// no input job file.  not a restart.

	ptr = fgets(line,GSTR_MAXSIZE,in);
	if (ptr == NULL)
		return 0;	// job file is empty.  not a restart.

	fclose(in);

	if (line[0] != 'n')
		return 0;	// malformed job file.  not a restart.

	zInit(&tmpz);
	z2decstr(N,&gstr1);
	str2hexz(line+3,&tmpz);
	ans = zCompare(&tmpz,N);
	if (ans == 0)
	{
		// ok, we have a job file for the current input.  this is a restart.
		ans = 1;

		if (VFLAG > 0)
			printf("nfs: commencing NFS restart\n");

		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "nfs: commencing NFS restart\n");
			fclose(logfile);
		}

		// attempt to open data file
		in = fopen(GGNFS_OUTPUTFILE,"r");
		if (in == NULL)
		{
			// data file doesn't exist, return flag to start sieving from the beginning
			*last_spq = 0;
		}
		else if (GGNFS_RESTART_FLAG)
		{
			// found it.  read the last specialq
			fclose(in);

			if (VFLAG > 0)
				printf("nfs: previous data file found - commencing search for last special-q\n");

			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "nfs: previous data file found - commencing search for last special-q\n");
				fclose(logfile);
			}

			//tail isn't good enough, because prior filtering steps could have inserted
			//free relations, which don't have a special q to read.
			//our task is to find the last valid line of the file, and then read
			//the last numeric entry from that line.
			if (1)
			{
				char line2[GSTR_MAXSIZE];

				in = fopen(GGNFS_OUTPUTFILE,"r");
				ptr = fgets(line, GSTR_MAXSIZE, in);
				if (ptr == NULL)
				{
					fclose(in);
					*last_spq = 0;
					return ans;
				}

				// crawl through the entire data file to find the next to last line
				while (!feof(in))
				{
					//throw away the next three lines.  if we encounter an end of file
					//in doing so, examine the last valid line
					ptr = fgets(line2, GSTR_MAXSIZE, in);
					if (ptr == NULL)
					{
						*last_spq = get_spq(line);
						return ans;
					}					

					ptr = fgets(line2, GSTR_MAXSIZE, in);
					if (ptr == NULL)
					{
						*last_spq = get_spq(line);
						return ans;
					}

					ptr = fgets(line2, GSTR_MAXSIZE, in);
					if (ptr == NULL)
					{
						*last_spq = get_spq(line);
						return ans;
					}

					//get another line, and check to make sure its valid.  if not,
					//keep going
					ptr = fgets(line2, GSTR_MAXSIZE, in);
					if (ptr == NULL)
					{
						*last_spq = get_spq(line);
						return ans;
					}
					else
					{
						if (strlen(line2) > 30)
							strcpy(line, line2);
					}
					job->current_rels += 4;

				}
			}
		}
		else
		{
			printf("must specify -R to restart when a savefile already exists\n");	

			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "nfs: refusing to restart without -R option\n");
				fclose(logfile);
			}

			*last_spq = 0;
			ans = -1;
		}

	}
	else
	{
		// job file is for a different input.  not a restart.
		ans = 0;
		*last_spq = 0;
	}

	zFree(&tmpz);
	return ans;

}

uint32 get_spq(char *line)
{	
	//next to the last line is in line2
	int i;
	uint32 ans;
	FILE *logfile;

	if (VFLAG > 0)
		printf("nfs: parsing special-q from line: %s\n",line);

	logfile = fopen(flogname, "a");
	if (logfile == NULL)
		printf("could not open yafu logfile for appending\n");
	else
	{
		logprint(logfile, "nfs: parsing special-q from line: %s\n",line);
		fclose(logfile);
	}

	ans = 0;
	//note: this assumes algebraic side sieving
	for (i=strlen(line); i>=0; i--)
	{
		if (line[i] == ',')
			break;
	}
	sscanf(line+i+1,"%x",&ans);
	return ans;
}

void find_best_msieve_poly(z *N, ggnfs_job_t *job)
{
	// parse a msieve.dat.p file to find the best polynomial (based on e score)
	// output this as a ggnfs polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE], *ptr;
	double score, bestscore = 0;
	int count, bestline = 0, i;

	sprintf(line, "%s.p",GGNFS_OUTPUTFILE);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",line);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for reading!\n",line);
			fclose(logfile);
		}
		exit(1);
	}

	// read and count polys of the file
	if (VFLAG > 0)
		printf("searching for best poly...\n");

	count = 0;
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;		

		if (line[0] == '#')
		{
			count++;
			if (VFLAG > 0)
				printf("found poly: %s",line);
			// the comment line for a new polynomial.  read characters until we 
			// find the e-score designator
			i=1;
			while (i < (strlen(line) - 1))
			{
				if (line[i] == 'e' && line[i+1] == ' ')
					break;
				if (line[i] == 'e' && line[i+1] == 9)
					break;
				i++;
			}

			if (line[i] != 'e')
			{
				//didn't find it.  ignore this poly
				continue;
			}
			else
			{
				//found it, read the score
				sscanf(line + i + 2, "%le", &score);
				//printf("reading: %s",line + i + 2);
				if (score > bestscore)
				{
					bestscore = score;
					bestline = count;
					if (VFLAG > 0)
						printf("new best score = %e, new best line = %d\n",bestscore, bestline);
				}
			}
		}
		
	}
	fclose(in);

	//open it again
	sprintf(line, "%s.p",GGNFS_OUTPUTFILE);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",GGNFS_OUTPUTFILE);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for reading!\n",GGNFS_OUTPUTFILE);
			fclose(logfile);
		}
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(GGNFS_JOB_INFILE,"w");
	if (out == NULL)
	{
		printf("could not open %s for writing!\n",GGNFS_JOB_INFILE);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for writing!\n",GGNFS_JOB_INFILE);
			fclose(logfile);
		}
		exit(1);
	}

	//go to the bestline
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;		

		if (line[0] == '#')
			bestline--;

		if (bestline == 0)
		{
			if (VFLAG > 0)
				printf("best poly: \n%s",line);

			logfile = fopen(flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "best poly: %s",line);
				fclose(logfile);
			}

			break;
		}
	}

	//copy n into the job file
	z2decstr(N,&gstr1);
	fprintf(out, "n: %s\n",gstr1.s);

	if (VFLAG > 0)
		printf("n: %s\n",gstr1.s);

	// copy out the poly
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;

		if (line[0] == '#')
			break;		
		
		fputs(line,out);

		if (VFLAG > 0)
			printf("%s",line);
	}

	// and copy in the job parameters
	fprintf(out,"rlim: %u\n",job->fblim);
	fprintf(out,"alim: %u\n",job->fblim);
	fprintf(out,"lpbr: %u\n",job->lpb);
	fprintf(out,"lpba: %u\n",job->lpb);
	fprintf(out,"mfbr: %u\n",job->mfb);
	fprintf(out,"mfba: %u\n",job->mfb);
	fprintf(out,"rlambda: %f\n",job->lambda);
	fprintf(out,"alambda: %f\n",job->lambda);

	fclose(in);
	fclose(out);

	return;
}

void msieve_to_ggnfs(ggnfs_job_t *job)
{
	// convert a msieve.fb polynomial into a ggnfs polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;

	in = fopen(GGNFS_FBFILE,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",GGNFS_FBFILE);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for reading!\n",GGNFS_FBFILE);
			fclose(logfile);
		}
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(GGNFS_JOB_INFILE,"w");
	if (out == NULL)
	{
		printf("could not open %s for writing!\n",GGNFS_JOB_INFILE);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for writing!\n",GGNFS_JOB_INFILE);
			fclose(logfile);
		}
		exit(1);
	}

	//translate the polynomial info
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;

		if (line[0] == 'N')
			sprintf(outline, "n: %s",line + 2);
		else if(line[0] == 'S')
			sprintf(outline, "skew: %s",line + 5);
		else if (line[0] == 'R')
			sprintf(outline, "Y%c: %s",line[1], line + 3);
		else if (line[0] == 'A')
			sprintf(outline, "c%c: %s",line[1], line + 3);
		else
			strcpy(outline, line);	//copy the line (probably white space)

		fputs(outline,out);
	}

	// and copy in the job parameters
	fprintf(out,"rlim: %u\n",job->fblim);
	fprintf(out,"alim: %u\n",job->fblim);
	fprintf(out,"lpbr: %u\n",job->lpb);
	fprintf(out,"lpba: %u\n",job->lpb);
	fprintf(out,"mfbr: %u\n",job->mfb);
	fprintf(out,"mfba: %u\n",job->mfb);
	fprintf(out,"rlambda: %f\n",job->lambda);
	fprintf(out,"alambda: %f\n",job->lambda);

	fclose(in);
	fclose(out);

	return;
}

void ggnfs_to_msieve(ggnfs_job_t *job)
{
	// convert a ggnfs.job file into a msieve.fb polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;

	in = fopen(GGNFS_JOB_INFILE,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",GGNFS_JOB_INFILE);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for reading!\n",GGNFS_JOB_INFILE);
			fclose(logfile);
		}
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(GGNFS_FBFILE,"w");
	if (out == NULL)
	{
		printf("could not open %s for writing!\n",GGNFS_FBFILE);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for writing!\n",GGNFS_FBFILE);
			fclose(logfile);
		}
		exit(1);
	}

	//translate the polynomial info
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;

		if (line[0] == 'n')
			sprintf(outline, "N %s",line + 3);
		else if(line[0] == 's')
			sprintf(outline, "SKEW %s",line + 5);
		else if (line[0] == 'Y')
			sprintf(outline, "R%c %s",line[1], line + 4);
		else if (line[0] == 'c')
			sprintf(outline, "A%c %s",line[1], line + 4);
		else
			strcpy(outline, "");

		fputs(outline,out);
	}

	fclose(in);
	fclose(out);

	return;
}

//digits, r/alim, lpbr/a, mfbr/a, r/alambda, siever, rels
//entries based on statistics gathered from many factorizations done
//over the years by myself and others, and from here:
//http://www.mersenneforum.org/showthread.php?t=12365
#define GGNFS_TABLE_ROWS 15
static double ggnfs_table[GGNFS_TABLE_ROWS][8] = {
	{85,  900000,   24, 48, 2.1, 11, 0, 10000},
	{90,  1200000,  25, 50, 2.3, 11, 0, 20000},
	{95,  1500000,  25, 50, 2.5, 12, 0, 20000},
	{100, 1800000,  26, 52, 2.5, 12, 0, 40000},
	{105, 2500000,  26, 52, 2.5, 12, 0, 40000},
	{110, 3200000,  26, 52, 2.5, 13, 0, 40000},
	{115, 4500000,  27, 54, 2.5, 13, 0, 80000},
	{120, 5000000,  27, 54, 2.5, 13, 0, 80000},
	{125, 5500000,  27, 54, 2.5, 13, 0, 80000},
	{130, 6000000,  27, 54, 2.5, 13, 0, 80000},
	{135, 8000000,  27, 54, 2.5, 14, 0, 80000},
	{140, 12000000, 28, 56, 2.5, 14, 0, 160000},
	{145, 15000000, 28, 56, 2.5, 14, 0, 160000},
	{150, 20000000, 29, 58, 2.5, 14, 0, 320000},
	{155, 30000000, 29, 58, 2.5, 15, 0, 320000}
};

void get_ggnfs_params(z *N, ggnfs_job_t *job)
{
	// based on the size/difficulty of the input number, determine "good" parameters
	// for the following: factor base limits, factor base large prime bound, trial
	// division fudge factors, trial division cutoffs, ggnfs lattice siever area, 
	// and expected number of relations needed.  linearly interpolate between table
	// entries.  keep last valid entry off the ends of the table.  This will produce
	// increasingly poor choices as one goes farther off the table, but you should be
	// doing things by hand by then anyway.
	int i, d = ndigits(N);
	double scale;
	FILE *test, *logfile;
	
	job->min_rels = 0;
	for (i=0; i<GGNFS_TABLE_ROWS - 1; i++)
	{
		if (d > ggnfs_table[i][0] && d <= ggnfs_table[i+1][0])
		{
			scale = (double)(ggnfs_table[i+1][0] - d) /
				(double)(ggnfs_table[i+1][0] - ggnfs_table[i][0]);

			//interp
			job->fblim = ggnfs_table[i+1][1] - 
				(uint32)(scale * (double)(ggnfs_table[i+1][1] - ggnfs_table[i][1]));

			//pick closest entry
			if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
				job->lpb = ggnfs_table[i][2];
			else
				job->lpb = ggnfs_table[i+1][2];
			
			//pick closest entry
			if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
				job->mfb = ggnfs_table[i][3];
			else
				job->mfb = ggnfs_table[i+1][3];

			//pick closest entry
			if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
				job->lambda = ggnfs_table[i][4];
			else
				job->lambda = ggnfs_table[i+1][4];

			//pick closest entry
			if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
				job->siever = ggnfs_table[i][5];
			else
				job->siever = ggnfs_table[i+1][5];

			job->min_rels = pow(2.0,(double)job->lpb) / log(pow(2.0,(double)job->lpb));

			//interp
			job->qrange = ggnfs_table[i+1][7] - 
				(uint32)(scale * (double)(ggnfs_table[i+1][7] - ggnfs_table[i][7]));
		}
	}

	if (job->min_rels == 0)
	{
		//couldn't find a table entry
		if (d <= ggnfs_table[0][0])
		{
			job->fblim = ggnfs_table[0][1];
			job->lpb = ggnfs_table[0][2];
			job->mfb = ggnfs_table[0][3];
			job->lambda = ggnfs_table[0][4];
			job->siever = ggnfs_table[0][5];
			job->min_rels = pow(2.0,(double)job->lpb) / log(pow(2.0,(double)job->lpb));
			job->qrange = ggnfs_table[0][7];
		}
		else
		{
			job->fblim = ggnfs_table[GGNFS_TABLE_ROWS-1][1];
			job->lpb = ggnfs_table[GGNFS_TABLE_ROWS-1][2];
			job->mfb = ggnfs_table[GGNFS_TABLE_ROWS-1][3];
			job->lambda = ggnfs_table[GGNFS_TABLE_ROWS-1][4];
			job->siever = ggnfs_table[GGNFS_TABLE_ROWS-1][5];
			job->min_rels = pow(2.0,(double)job->lpb) / log(pow(2.0,(double)job->lpb));
			job->qrange = ggnfs_table[GGNFS_TABLE_ROWS-1][7];
		}
	}

	switch (job->siever)
	{
	case 11:
		sprintf(job->sievername, "%sgnfs-lasieve4I11e", ggnfs_dir);
		break;
	case 12:
		sprintf(job->sievername, "%sgnfs-lasieve4I12e", ggnfs_dir);
		break;
	case 13:
		sprintf(job->sievername, "%sgnfs-lasieve4I13e", ggnfs_dir);
		break;
	case 14:
		sprintf(job->sievername, "%sgnfs-lasieve4I14e", ggnfs_dir);
		break;
	case 15:
		sprintf(job->sievername, "%sgnfs-lasieve4I15e", ggnfs_dir);
		break;
	}

#if defined(WIN32)
	sprintf(job->sievername, "%s.exe", job->sievername);
#endif

	// test for existence of the siever
	test = fopen(job->sievername, "rb");
	if (test == NULL)
	{
		printf("could not find %s, bailing\n",job->sievername);
		logfile = fopen(flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not find %s, bailing\n",job->sievername);
			fclose(logfile);
		}
		exit(-1);
	}

	return;
}
	
#else

void test_msieve_gnfs(fact_obj_t *fobj)
{
	printf("gnfs has been disabled\n");
	return;
}

#endif

