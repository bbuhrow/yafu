#include "gnfs.h"
#include "yafu_string.h"
#include "arith.h"
#include "factor.h"
#include "qs.h"


#define USE_NFS

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
	char *polyfilename, *logfilename, *fbfilename;
	uint64 poly_lower;
	uint64 poly_upper;
	msieve_obj *obj;
	mp_t *mpN;
	factor_list_t *factor_list;
	struct timeval thread_start_time;
	fact_obj_t *fobj;

	int tindex;
	int is_poly_select;

	/* fields for thread pool synchronization */
	volatile enum nfs_thread_command command;
	volatile int *thread_queue, *threads_waiting;

#if defined(WIN32) || defined(_WIN64)
	HANDLE thread_id;
	HANDLE run_event;
	HANDLE finish_event;

	HANDLE *queue_event;
	HANDLE *queue_lock;
#else
	pthread_t thread_id;
	pthread_mutex_t run_lock;
	pthread_cond_t run_cond;

	pthread_mutex_t *queue_lock;
	pthread_cond_t *queue_cond;
#endif

} nfs_threaddata_t;

//----------------------- LOCAL FUNCTIONS -------------------------------------//
void *lasieve_launcher(void *ptr);
void *polyfind_launcher(void *ptr);
void find_best_msieve_poly(fact_obj_t *fobj, z *N, ggnfs_job_t *job);
void msieve_to_ggnfs(fact_obj_t *fobj, ggnfs_job_t *job);
void ggnfs_to_msieve(fact_obj_t *fobj, ggnfs_job_t *job);
void get_ggnfs_params(fact_obj_t *fobj, z *N, ggnfs_job_t *job);
int check_existing_files(fact_obj_t *fobj, z *N, uint32 *last_spq, ggnfs_job_t *job);
void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj);
uint32 get_spq(char **lines, int last_line, fact_obj_t *fobj);
uint32 do_msieve_filtering(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN);
void do_msieve_polyselect(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN, factor_list_t *factor_list);
void get_polysearch_params(fact_obj_t *fobj, uint64 *start, uint64 *range);
void init_poly_threaddata(nfs_threaddata_t *t, msieve_obj *obj, 
	mp_t *mpN, factor_list_t *factor_list, int tid, uint32 flags, uint64 start, uint64 stop);
void do_sieving(fact_obj_t *fobj, ggnfs_job_t *job);
void savefile_concat(char *filein, char *fileout, msieve_obj *mobj);
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


#ifdef USE_NFS

void nfsexit(int sig)
{
	printf("\nAborting... \n");
	NFS_ABORT = 1;
	return;
}


//----------------------- NFS ENTRY POINT ------------------------------------//
void nfs(fact_obj_t *fobj)
{
	z *N = &fobj->nfs_obj.n;
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
	mpz_t gmpz;
	int i;

	//below a certain amount, revert to SIQS
	if (ndigits(N) < fobj->nfs_obj.min_digits)
	{
		zCopy(N, &fobj->qs_obj.n);
		SIQS(fobj);
		return;
	}	
		
	//first a small amount of trial division
	//which will add anything found to the global factor list
	//this is only called with the main thread
	if (VFLAG > 0)
		printf("nfs: commencing trial factoring\n");
	zCopy(N,&fobj->div_obj.n);	
	fobj->div_obj.print = 0;
	fobj->div_obj.limit = 10000;
	zTrial(fobj);
	zCopy(&fobj->div_obj.n,N);

	if (isPrime(N))
	{
		add_to_factor_list(fobj, N);
		zCopy(&zOne,N);

		logfile = fopen(fobj->flogname, "a");
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

	if (isSquare(N))
	{
		zNroot(N,N,2);
		N->type = PRP;
		add_to_factor_list(fobj, N);
		logfile = fopen(fobj->flogname, "a");
		if (logfile != NULL)
			logprint(logfile,"prp%d = %s\n",ndigits(N),z2decstr(N,&gstr1));
		add_to_factor_list(fobj, N);
		if (logfile != NULL)
			logprint(logfile,"prp%d = %s\n",ndigits(N),z2decstr(N,&gstr1));
		zCopy(&zOne,N);
		fclose(logfile);
		return;
	}

	mpz_init(gmpz);
	mpz_import(gmpz, abs(N->size), -1, sizeof(fp_digit), 
		0, (size_t)0, N->val);
	if (N->size < 0)
		mpz_neg(gmpz, gmpz);
	i = mpz_perfect_power_p(gmpz);
	mpz_clear(gmpz);

	if (i)
	{
		printf("input is a perfect power\n");
		add_to_factor_list(fobj, N);
		logfile = fopen(fobj->flogname, "a");
		if (logfile != NULL)
		{
			logprint(logfile,"input is a perfect power\n");
			logprint(logfile,"c%d = %s\n",ndigits(N),z2decstr(N,&gstr1));
		}
		fclose(logfile);
		return;
	}

	//check for inconsistent options
	if (((fobj->nfs_obj.startq > 0) && !(fobj->nfs_obj.rangeq > 0)) ||
		(!(fobj->nfs_obj.startq > 0) && (fobj->nfs_obj.rangeq > 0)))
	{
		printf("-f and -c must be specified together\n");
		exit(1);
	}

	if (((fobj->nfs_obj.startq > 0) || (fobj->nfs_obj.rangeq > 0)) && fobj->nfs_obj.poly_only)
	{
		printf("bad options: -np with -f or -c\n");
		exit(1);
	}

	if (((fobj->nfs_obj.startq > 0) || (fobj->nfs_obj.rangeq > 0)) && fobj->nfs_obj.post_only)
	{
		printf("bad options: -nc with -f or -c\n");
		exit(1);
	}

	if (fobj->nfs_obj.sieve_only && fobj->nfs_obj.post_only)
	{
		printf("bad options: -nc with -ns\n");
		exit(1);
	}

	if (fobj->nfs_obj.sieve_only && fobj->nfs_obj.poly_only)
	{
		printf("bad options: -np with -ns\n");
		exit(1);
	}

	if (fobj->nfs_obj.poly_only && fobj->nfs_obj.post_only)
	{
		printf("bad options: -nc with -np\n");
		exit(1);
	}

	//find best job parameters
	get_ggnfs_params(fobj, N,&job);			

	//if we are going to be doing sieving, check for the sievers
	if (!(fobj->nfs_obj.poly_only || fobj->nfs_obj.post_only))
	{
		FILE *test;
		// test for existence of the siever
		test = fopen(job.sievername, "rb");
		if (test == NULL)
		{
			printf("WARNING: could not find %s, reverting to siqs!\n",job.sievername);
			logfile = fopen(fobj->flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "WARNING: could not find %s, reverting to siqs!\n",job.sievername);
				fclose(logfile);
			}
			zCopy(N, &fobj->qs_obj.n);
			SIQS(fobj);
			return;
		}
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

			logfile = fopen(fobj->flogname, "a");
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
			//this will initialize the savefile to the outputfile name provided
			obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
				fobj->nfs_obj.fbfile, seed1, seed2, max_relations, nfs_lower, nfs_upper, cpu, 
				L1CACHE, L2CACHE, THREADS, mem_mb, which_gpu, 0.0);

			fobj->nfs_obj.mobj = obj;

			//convert input to msieve bigint notation and initialize a list of factors
			z2mp_t(N,&mpN);
			factor_list_init(&factor_list);

			//see if we can resume a factorization based on the combination of input number,
			//.job file, .fb file, .dat file, and/or .p file.  else, start new job.
			job.current_rels = 0;
			is_continuation = check_existing_files(fobj, N, &last_specialq, &job);						

			//determine sieving start value and range.
			if (fobj->nfs_obj.rangeq > 0)
				qrange = ceil((double)fobj->nfs_obj.rangeq / (double)THREADS);
			else
				qrange = job.qrange;

			if (is_continuation > 0)
			{		
				if (last_specialq == 0)
				{
					if (VFLAG >= 0)
						printf("nfs: continuing job - could not determine last special q; using default startq\n");

					logfile = fopen(fobj->flogname, "a");
					if (logfile == NULL)
						printf("could not open yafu logfile for appending\n");
					else
					{
						logprint(logfile, "nfs: continuing job - could not determine last special q; using default startq\n");
						fclose(logfile);
					}

					if (fobj->nfs_obj.startq > 0)
						startq = fobj->nfs_obj.startq;
					else
						startq = job.fblim / 2;

					// next step is sieving
					statenum = 2;
				}
				else
				{
					if (VFLAG >= 0)
						printf("nfs: found %u relations, continuing job at specialq = %u\n",
						job.current_rels,last_specialq);

					logfile = fopen(fobj->flogname, "a");
					if (logfile == NULL)
						printf("could not open yafu logfile for appending\n");
					else
					{
						logprint(logfile, "nfs: found %u relations, continuing job at specialq = %u\n",
							job.current_rels,last_specialq);
						fclose(logfile);
					}

					if (fobj->nfs_obj.startq > 0)
					{
						startq = fobj->nfs_obj.startq;
						statenum = 2;
					}
					else if (fobj->nfs_obj.sieve_only)
					{
						startq = last_specialq;
						statenum = 2;
					}
					else
					{
						startq = last_specialq;

						// since we apparently found some relations, see if it is enough
						statenum = 8;
					}
				}
			}
			else if (is_continuation == 0)
			{
				// new factorization				
				if ((fobj->nfs_obj.startq > 0) || fobj->nfs_obj.sieve_only)
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
			if (fobj->nfs_obj.poly_only)
				statenum = 1;
			else if (fobj->nfs_obj.sieve_only)
				statenum = 2;
			else if (fobj->nfs_obj.post_only)
				statenum = 3;

			break;

		case 1: //"polysearch":
			
			do_msieve_polyselect(fobj, obj, &job, &mpN, &factor_list);
			
			if (fobj->nfs_obj.poly_only)
				process_done = 1;
			statenum = 2;
			break;

		case 2: //"sieve":

			do_sieving(fobj, &job);

			//see how we're doing
			if (fobj->nfs_obj.sieve_only || (fobj->nfs_obj.rangeq > 0))
				process_done = 1;
			else
			{
				statenum = 8;
			}

			break;

		case 3: //"filter":

			relations_needed = do_msieve_filtering(fobj, obj, &job, &mpN);
			if (relations_needed == 0)
				statenum = 4;	//proceed to linear algebra
			else
			{
				if (fobj->nfs_obj.post_only)
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
	
			logfile = fopen(fobj->flogname, "a");
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

			logfile = fopen(fobj->flogname, "a");
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
			
			remove(fobj->nfs_obj.outputfile);
			remove(fobj->nfs_obj.fbfile);
			sprintf(tmpstr, "%s.p",fobj->nfs_obj.outputfile);	remove(tmpstr);			
			sprintf(tmpstr, "%s.br",fobj->nfs_obj.outputfile);	remove(tmpstr);
			sprintf(tmpstr, "%s.cyc",fobj->nfs_obj.outputfile);	remove(tmpstr);
			sprintf(tmpstr, "%s.dep",fobj->nfs_obj.outputfile);	remove(tmpstr);
			sprintf(tmpstr, "%s.hc",fobj->nfs_obj.outputfile);	remove(tmpstr);
			sprintf(tmpstr, "%s.mat",fobj->nfs_obj.outputfile);	remove(tmpstr);	
			sprintf(tmpstr, "%s.lp",fobj->nfs_obj.outputfile);	remove(tmpstr);
			sprintf(tmpstr, "%s.d",fobj->nfs_obj.outputfile);	remove(tmpstr);

			gettimeofday(&stop, NULL);

			difference = my_difftime (&start, &stop);

			t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
			free(difference);	

			if (VFLAG >= 0)
				printf("NFS elapsed time = %6.4f seconds.\n",t_time);

			logfile = fopen(fobj->flogname, "a");
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

		case 8:		// check to see if we should try filtering
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
		if ((fobj->nfs_obj.timeout > 0) && (t_time > (double)fobj->nfs_obj.timeout))
		{
			if (VFLAG >= 0)
				printf("NFS timeout after %6.4f seconds.\n",t_time);

			logfile = fopen(fobj->flogname, "a");
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

void savefile_concat(char *filein, char *fileout, msieve_obj *mobj)
{
	FILE *in;

	in = fopen(filein,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading\n",filein);
		exit(-1);
	}

	savefile_open(&mobj->savefile, SAVEFILE_APPEND);

	while (!feof(in))
	{
		char tmpline[GSTR_MAXSIZE], *tmpptr;
		tmpptr = fgets(tmpline, GSTR_MAXSIZE, in);					
		if (tmpptr == NULL)
			break;
		else
			savefile_write_line(&mobj->savefile, tmpline);
	}
	fclose(in);

	savefile_flush(&mobj->savefile);
	savefile_close(&mobj->savefile);

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

void do_msieve_polyselect(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN, factor_list_t *factor_list)
{
	FILE *logfile;
	uint32 flags;
	z N;

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
		printf("could not open yafu logfile for appending\n");
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
	total_time = 0;
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

			if (!is_startup)
			{
				// combine thread's output with main file
#if defined(WIN32)
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

#else
				sprintf(syscmd,"cat %s.p >> %s",t->polyfilename,master_polyfile);
				system(syscmd);
#endif

			}

			// remove each thread's .p file after copied
			sprintf(syscmd, "%s.p",t->polyfilename); 
			remove(syscmd);	

			// also remove the temporary log file
			sprintf(syscmd,"%s.%d",fobj->nfs_obj.logfile,tid);
			remove(syscmd);

			//free data
			free(t->logfilename);
			free(t->polyfilename);
			msieve_obj_free(t->obj);			

			//check the time
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
				estimated_range_time = (uint32)(total_time / (double)num_ranges);
			}
						
			// check for an abort
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
		printf("could not open yafu logfile for appending\n");
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

	// parse the master .p file and find the best poly	
	find_best_msieve_poly(fobj, &N, job);
	//also create a .fb file
	ggnfs_to_msieve(fobj, job);
	zFree(&N);

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

uint32 do_msieve_filtering(fact_obj_t *fobj, msieve_obj *obj, ggnfs_job_t *job, mp_t *mpN)
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

	logfile = fopen(fobj->flogname, "a");
	if (logfile == NULL)
		printf("could not open yafu logfile for appending\n");
	else
	{
		logprint(logfile, "nfs: commencing msieve filtering\n");
		fclose(logfile);
	}

	tmp = fopen(fobj->nfs_obj.fbfile,"r");
	if (tmp == NULL)
		ggnfs_to_msieve(fobj, job);
	else
		fclose(tmp);

	printf("%s\n",mp_print(mpN, 10, NULL, gstr1.s));
	relations_needed = nfs_filter_relations(obj, mpN);

	return relations_needed;
}

void nfs_start_worker_thread(nfs_threaddata_t *t, 
				uint32 is_master_thread) {

	/* create a thread that will process a polynomial. The last poly does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread == 1) {
		return;
	}

	t->command = NFS_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
		
	// specific to different structure of poly selection threading
	if (is_master_thread == 2)
	{
		t->run_event = CreateEvent(NULL, FALSE, FALSE, NULL);
		t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
		*t->queue_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	}
	else
	{
		t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
		t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	}

	t->thread_id = CreateThread(NULL, 0, nfs_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	if (is_master_thread == 0)
	{
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
	}

	pthread_create(&t->thread_id, NULL, nfs_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != NFS_COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);

	if (is_master_thread == 2)
		pthread_mutex_unlock(&t->run_lock);

#endif

}

void nfs_stop_worker_thread(nfs_threaddata_t *t,
				uint32 is_master_thread)
{
	if (is_master_thread == 1) {
		return;
	}

	t->command = NFS_COMMAND_END;
#if defined(WIN32) || defined(_WIN64)
	SetEvent(t->run_event);
	WaitForSingleObject(t->thread_id, INFINITE);
	CloseHandle(t->thread_id);
	CloseHandle(t->run_event);
	CloseHandle(t->finish_event);

	// specific to different structure of poly selection threading
	if (is_master_thread == 2)
		CloseHandle(*t->queue_event);
#else
	if (is_master_thread == 2)
		pthread_mutex_lock(&t->run_lock);

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

	/*
    * Respond to the master thread that we're ready for work. If we had any thread-
    * specific initialization which needed to be done, it would go before this signal.
    */

	// specific to different structure of poly selection threading
	if (t->is_poly_select)
	{
#if defined(WIN32) || defined(_WIN64)
		t->command = NFS_COMMAND_WAIT;
		SetEvent(t->finish_event);
#else
		pthread_mutex_lock(&t->run_lock);
		t->command = NFS_COMMAND_WAIT;
		pthread_cond_signal(&t->run_cond);
		pthread_mutex_unlock(&t->run_lock);
#endif
	}

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

		if (t->is_poly_select)
		{
#if defined(WIN32) || defined(_WIN64)

			WaitForSingleObject( 
				*t->queue_lock,    // handle to mutex
				INFINITE);  // no time-out interval
 			
			t->thread_queue[(*(t->threads_waiting))++] = t->tindex;

			SetEvent(*t->queue_event);			

			ReleaseMutex(*t->queue_lock);
		
#else
			pthread_mutex_unlock(&t->run_lock);

			// lock the work queue and insert my thread ID into it
			// this tells the master that my results should be collected
			// and I should be dispatched another polynomial
			pthread_mutex_lock(t->queue_lock);
			t->thread_queue[(*(t->threads_waiting))++] = t->tindex;
			pthread_cond_signal(t->queue_cond);
			pthread_mutex_unlock(t->queue_lock);
#endif
		}
		else
		{

#if defined(WIN32) || defined(_WIN64)
			SetEvent(t->finish_event);		
#else
			pthread_cond_signal(&t->run_cond);
			pthread_mutex_unlock(&t->run_lock);
#endif

		}
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
	//factor_integer(t->obj->input, t->obj->flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, fobj->nfs_obj.fbfile, 
	//	&g_rand.low, &g_rand.hi, 0, t->obj->nfs_lower, t->obj->nfs_upper, t->obj->cpu, 
	//	L1CACHE, L2CACHE, THREADS, 0, 0, 0.0);
	factor_gnfs(t->obj, t->mpN, t->factor_list);	

	return 0;
}


void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj)
{
	z tmp1, tmp2, *N = &fobj->nfs_obj.n;
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

		logfile = fopen(fobj->flogname, "a");
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

		logfile = fopen(fobj->flogname, "a");
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


int check_existing_files(fact_obj_t *fobj, z *N, uint32 *last_spq, ggnfs_job_t *job)
{
	// see if we can resume a factorization based on the combination of input number,
	// .job file, .fb file, .dat file, and/or .p file.  else, start new job.
	// be very cautious about overwriting .dat files, as building these up can 
	// represent a *large* amount of work.
	// the only time .dat files should be touched is if the user forces us to (with -R)
	// or if an NFS process completes normally (in which case, we delete the .dat files).

	FILE *in, *logfile;
	char line[GSTR_MAXSIZE], *ptr;
	int ans;
	z tmpz;
	msieve_obj *mobj = fobj->nfs_obj.mobj;

	in = fopen(fobj->nfs_obj.job_infile,"r");
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
		// ok, we have a job file for the current input.  this is a restart of sieving
		ans = 1;

		if (VFLAG > 0)
			printf("nfs: commencing NFS restart\n");

		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "nfs: commencing NFS restart\n");
			fclose(logfile);
		}

		// attempt to open data file
		if (!savefile_exists(&mobj->savefile))
		{
			// data file doesn't exist, return flag to start sieving from the beginning
			*last_spq = 0;
		}
		else if (fobj->nfs_obj.restart_flag)
		{
			// found it.  read the last specialq
			if (VFLAG > 0)
				printf("nfs: previous data file found - commencing search for last special-q\n");

			logfile = fopen(fobj->flogname, "a");
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
			//furthermore, the last valid line could be incomplete if a job is interrupted.  so
			//we must be able to check the next to last valid line as well.
			//finally, we must be able to parse the special q from the line, which could have
			//been from a rational side or algebraic side sieve.  unfortunately, there is
			//no way to tell what previous jobs may have done 
			//the procedure is as follows:  parse a line and grab both the alg and rat hex values
			//if the number we picked is actually a special-q, then it will have numbers of similar
			//size (or the same) in identical positions in other lines.  it should also be above
			//a certain bound, and be prime.  if it meets all of these criteria, then we probably
			//used the right interpretation.
			if (1)
			{
				char **lines, tmp[GSTR_MAXSIZE];
				int line;
				int i;

				lines = (char **)malloc(4 * sizeof(char *));
				for (i=0; i < 4; i++)
					lines[i] = (char *)malloc(GSTR_MAXSIZE * sizeof(char));

				savefile_open(&mobj->savefile, SAVEFILE_READ);
				savefile_read_line(lines[0], GSTR_MAXSIZE, &mobj->savefile);
				if (savefile_eof(&mobj->savefile))
				{
					savefile_close(&mobj->savefile);
					*last_spq = 0;
					for (i=0; i < 4; i++)
						free(lines[i]);
					free(lines);
					return ans;
				}

				// crawl through the entire data file to find the next to last line
				// TODO: can't we start from the end of the file somehow?
				line = 0;
				while (!feof(in))
				{
					// read a line into the next position of the circular buffer
					savefile_read_line(tmp, GSTR_MAXSIZE, &mobj->savefile);
					if (savefile_eof(&mobj->savefile))
						break;					

					// quick check that it might be a valid line
					if (strlen(tmp) > 30)
					{
						// wrap
						if (++line > 3) line = 0;
						// then copy
						strcpy(lines[line], tmp);					
					}

					// while we are at it, count the lines
					job->current_rels++;

				}
				savefile_close(&mobj->savefile);

				// now we are done and we have a buffer with the last 4 valid lines
				// throw away the last one, which may be malformed, and extract
				// the special q from the other 3.
				for (i=0; i<4; i++)
					printf("line %d = %s\n",i, lines[i]);
				*last_spq = get_spq(lines, line, fobj);

				for (i=0; i < 4; i++)
					free(lines[i]);
				free(lines);
			}
		}
		else
		{
			printf("must specify -R to restart when a savefile already exists\n");	

			logfile = fopen(fobj->flogname, "a");
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

uint32 get_spq(char **lines, int last_line, fact_obj_t *fobj)
{	
	//the last 4 valid lines are passed in
	int i, line, count;
	double var[2], avg[2];
	uint32 ans;
	FILE *logfile;
	uint32 rat[3], alg[3];
	char *ptr;

	if (VFLAG > 0)
		printf("nfs: parsing special-q\n");

	logfile = fopen(fobj->flogname, "a");
	if (logfile == NULL)
		printf("could not open yafu logfile for appending\n");
	else
	{
		logprint(logfile, "nfs: parsing special-q\n");
		fclose(logfile);
	}

	ans = 0;
	// grab the entry in both the rational side and algebraic side
	// special-q locations from 3 different lines
	line = last_line;
	for (count=0; count < 3; count++)
	{
		if (++line > 3) line = 0;

		// find a,b delimiter
		ptr = strstr(lines[line], ":");
		if (ptr == NULL)
		{
			rat[count] = alg[count] = 0;
			continue;
		}

		// find rat side delimiter
		ptr = strstr(ptr + 1, ":");
		if (ptr == NULL)
		{
			rat[count] = alg[count] = 0;
			continue;
		}

		// grab rat side entry
		for (i = (ptr - lines[line]) - 1; i >= 0; i--)
			if (lines[line][i] == ',')
				break;

		printf("parsing rat side spq from %s\n",lines[line]);
		sscanf(lines[line] + i + 1, "%x", &rat[count]);
		printf("found %x\n", rat[count]);

		// grab alg side entry
		for (i= strlen(lines[line]) - 1; i >= 0; i--)
			if (lines[line][i] == ',')
				break;

		printf("parsing alg side spq from %s\n",lines[line]);
		sscanf(lines[line] + i + 1, "%x", &alg[count]);
		printf("found %x\n", alg[count]);
	}

	// now we gotta make a decision.
	// if any two of the entries are less than 10000, eliminate
	// that side

	if (((rat[0] < 10000) && (rat[1] < 10000)) ||
		((rat[0] < 10000) && (rat[2] < 10000)) ||
		((rat[1] < 10000) && (rat[1] < 10000)))
	{
		// guess that it is an algebraic line
		return alg[0];
	}
	else if (((alg[0] < 10000) && (alg[1] < 10000)) ||
		((alg[0] < 10000) && (alg[2] < 10000)) ||
		((alg[1] < 10000) && (alg[1] < 10000)))
	{
		// guess that it is a rational line
		return rat[0];
	}
	
	// if all 3 are the same for either, then pick that one
	if ((rat[0] == rat[1]) && (rat[1] == rat[2]))
	{
		return rat[0];
	}
	else if ((alg[0] == alg[1]) && (alg[1] == alg[2]))
	{
		return alg[0];
	}

	// else, pick the one with lowest variance
	avg[0] = ((double)rat[0] + (double)rat[1] + (double)rat[2]) / 3;
	avg[1] = ((double)alg[0] + (double)alg[1] + (double)alg[2]) / 3;

	var[0] = (pow((double)rat[0] - avg[0],2) + 
		pow((double)rat[1] - avg[0],2) + 
		pow((double)rat[2] - avg[0],2)) / 3;
	var[1] = (pow((double)alg[0] - avg[1],2) + 
		pow((double)alg[1] - avg[1],2) + 
		pow((double)alg[2] - avg[1],2)) / 3;

	if (var[0] < var[1])
		return rat[0];
	else
		return alg[0];

}

void find_best_msieve_poly(fact_obj_t *fobj, z *N, ggnfs_job_t *job)
{
	// parse a msieve.dat.p file to find the best polynomial (based on e score)
	// output this as a ggnfs polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE], *ptr;
	double score, bestscore = 0;
	int count, bestline = 0, i;

	sprintf(line, "%s.p",fobj->nfs_obj.outputfile);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",line);
		logfile = fopen(fobj->flogname, "a");
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
	sprintf(line, "%s.p",fobj->nfs_obj.outputfile);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",fobj->nfs_obj.outputfile);
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for reading!\n",fobj->nfs_obj.outputfile);
			fclose(logfile);
		}
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(fobj->nfs_obj.job_infile,"w");
	if (out == NULL)
	{
		printf("could not open %s for writing!\n",fobj->nfs_obj.job_infile);
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for writing!\n",fobj->nfs_obj.job_infile);
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

			logfile = fopen(fobj->flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "nfs: best poly = %s",line);
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

void msieve_to_ggnfs(fact_obj_t *fobj, ggnfs_job_t *job)
{
	// convert a msieve.fb polynomial into a ggnfs polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;

	in = fopen(fobj->nfs_obj.fbfile,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",fobj->nfs_obj.fbfile);
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for reading!\n",fobj->nfs_obj.fbfile);
			fclose(logfile);
		}
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(fobj->nfs_obj.job_infile,"w");
	if (out == NULL)
	{
		printf("could not open %s for writing!\n",fobj->nfs_obj.job_infile);
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for writing!\n",fobj->nfs_obj.job_infile);
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

void ggnfs_to_msieve(fact_obj_t *fobj, ggnfs_job_t *job)
{
	// convert a ggnfs.job file into a msieve.fb polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;

	in = fopen(fobj->nfs_obj.job_infile,"r");
	if (in == NULL)
	{
		printf("could not open %s for reading!\n",fobj->nfs_obj.job_infile);
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for reading!\n",fobj->nfs_obj.job_infile);
			fclose(logfile);
		}
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(fobj->nfs_obj.fbfile,"w");
	if (out == NULL)
	{
		printf("could not open %s for writing!\n",fobj->nfs_obj.fbfile);
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			logprint(logfile, "could not open %s for writing!\n",fobj->nfs_obj.fbfile);
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
		//else if(line[0] == 's')
		//	sprintf(outline, "SKEW %s",line + 5);
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

void get_ggnfs_params(fact_obj_t *fobj, z *N, ggnfs_job_t *job)
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

	if (fobj->nfs_obj.siever > 0)
	{
		switch (fobj->nfs_obj.siever)
		{
		case 11:
			sprintf(job->sievername, "%sgnfs-lasieve4I11e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 12:
			sprintf(job->sievername, "%sgnfs-lasieve4I12e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 13:
			sprintf(job->sievername, "%sgnfs-lasieve4I13e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 14:
			sprintf(job->sievername, "%sgnfs-lasieve4I14e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 15:
			sprintf(job->sievername, "%sgnfs-lasieve4I15e", fobj->nfs_obj.ggnfs_dir);
			break;
		}
	}
	else
	{
		switch (job->siever)
		{
		case 11:
			sprintf(job->sievername, "%sgnfs-lasieve4I11e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 12:
			sprintf(job->sievername, "%sgnfs-lasieve4I12e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 13:
			sprintf(job->sievername, "%sgnfs-lasieve4I13e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 14:
			sprintf(job->sievername, "%sgnfs-lasieve4I14e", fobj->nfs_obj.ggnfs_dir);
			break;
		case 15:
			sprintf(job->sievername, "%sgnfs-lasieve4I15e", fobj->nfs_obj.ggnfs_dir);
			break;
		}
	}

#if defined(WIN32)
	sprintf(job->sievername, "%s.exe", job->sievername);
#endif

	return;
}
	
#else

void nfs(fact_obj_t *fobj)
{
	printf("gnfs has not been enabled\n");

	SIQS(fobj);

	return;
}

#endif

