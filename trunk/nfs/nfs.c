#include "gnfs.h"
#include "yafu_string.h"
#include "arith.h"
#include "factor.h"

enum nfs_thread_command {
	NFS_COMMAND_INIT,
	NFS_COMMAND_WAIT,
	NFS_COMMAND_RUN,
	NFS_COMMAND_RUN_TRANS,
	NFS_COMMAND_END
};

typedef struct
{
	uint32 fblim;
	uint32 lpb;
	uint32 mfb;
	float lambda;
	uint32 siever; 
	uint32 rels;
	uint32 qrange;
	char sievername[1024];
	uint32 startq;
} ggnfs_job_t;

typedef struct {
	char outfilename[80];
	ggnfs_job_t job;

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

void nfs_stop_worker_thread(nfs_threaddata_t *t,
				uint32 is_master_thread);
void nfs_start_worker_thread(nfs_threaddata_t *t, 
				uint32 is_master_thread);
#if defined(WIN32) || defined(_WIN64)
DWORD WINAPI nfs_worker_thread_main(LPVOID thread_data);
#else
void *nfs_worker_thread_main(void *thread_data);
#endif

void *lasieve_launcher(void *ptr);

void find_best_msieve_poly(z *N, ggnfs_job_t *job);
void msieve_to_ggnfs(ggnfs_job_t *job);
void get_ggnfs_params(z *N, ggnfs_job_t *job);
int check_existing_files(z *N, uint32 *last_spq);
void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj);

void test_msieve_gnfs(fact_obj_t *fobj)
{
	z *N = &fobj->qs_obj.n;
	char *input;
	str_t input_str;
	msieve_obj *obj;
	uint32 max_relations = 0;
	uint32 seed1 = g_rand.low;
	uint32 seed2 = g_rand.hi;
	uint64 nfs_lower = 0;
	uint64 nfs_upper = 0;
	enum cpu_type cpu = yafu_get_cpu_type();
	uint32 mem_mb = 0;
	uint32 which_gpu = 0;
	mp_t mpN;
	int status;
	factor_list_t factor_list;
	uint32 flags;
	FILE *fid;
	int i;
	ggnfs_job_t job;
	uint32 relations_needed = 1;	
	uint32 startq, qrange;
	int is_continuation;
	uint32 last_specialq = 0;
	nfs_threaddata_t *thread_data;		//an array of thread data objects
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	TIME_DIFF *	difference;
	double t_time;

	//start a counter for the whole job
	gettimeofday(&start, NULL);

	//write the input bigint as a string
	sInit(&input_str);
	input = z2decstr(N, &input_str);
	
	//set flags to do poly selection
	flags = 0;
	flags = flags | MSIEVE_FLAG_USE_LOGFILE;
	if (VFLAG > 1)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_POLY1;
	flags = flags | MSIEVE_FLAG_NFS_POLY2;

	//create an msieve_obj
	obj = msieve_obj_new(input, flags, NULL, NULL, NULL, seed1, seed2,
		max_relations, nfs_lower, nfs_upper, cpu, L1CACHE, L2CACHE, THREADS, mem_mb, which_gpu);

	//convert input to msieve bigint notation and initialize a list of factors
	z2mp_t(N,&mpN);
	factor_list_init(&factor_list);

	//check if this is a continuation of a factorization (defined when
	//a .fb file already exists for this N, or both a .fb and msieve.dat
	//file exist for this N, or if a ggnfs.job file exists for this N.
	//if any of the .job, .fb and/or .dat exist but for conflicting N's,
	//complain and stop.  else, start new job.
	is_continuation = check_existing_files(N, &last_specialq);

	//find best job parameters
	get_ggnfs_params(N,&job);

	qrange = job.qrange;
	if (is_continuation)
	{		
		if (last_specialq == 0)
		{
			printf("nfs: continuing job - could not determine last special q; using default startq\n");
			startq = job.fblim / 2;
		}
		else
		{
			printf("nfs: continuing job at specialq = %u\n",last_specialq);
			startq = last_specialq;
		}
	}
	else
	{
		//remove any previous .p file
		remove("msieve.dat.p");

		//msieve: polyselect
		if (VFLAG >= 0)
			printf("nfs: msieve poly select\n");
		status = factor_gnfs(obj, &mpN, &factor_list);		

		//convert best poly to ggnfs job file
		//check if polyselect made a .fb file
		fid = fopen("msieve.fb","r");
		if (fid == NULL)
		{
			// no, it didn't.  parse the .p file instead.
			find_best_msieve_poly(N, &job);
		}
		else
		{
			// yes it did, close the file and call the conversion utility
			fclose(fid);
			msieve_to_ggnfs(&job);
		}
		startq = job.fblim / 2;
	}	

	thread_data = (nfs_threaddata_t *)malloc(THREADS * sizeof(nfs_threaddata_t));
	for (i=0; i<THREADS; i++)
	{
		sprintf(thread_data[i].outfilename, "rels%d.dat", i);
		thread_data[i].job.fblim = job.fblim;
		thread_data[i].job.lambda = job.lambda;
		thread_data[i].job.lpb = job.lpb;
		thread_data[i].job.mfb = job.mfb;
		thread_data[i].job.qrange = job.qrange;
		thread_data[i].job.rels = job.rels;
		thread_data[i].job.siever = job.siever;
		thread_data[i].job.startq = job.startq;
		strcpy(thread_data[i].job.sievername, job.sievername);
	}

	/* activate the threads one at a time. The last is the
		  master thread (i.e. not a thread at all). */		
	for (i = 0; i < THREADS - 1; i++)
		nfs_start_worker_thread(thread_data + i, 0);

	nfs_start_worker_thread(thread_data + i, 1);

	//loop until filtering produces a matrix
	while (relations_needed != 0)
	{							
		char syscmd[1024];

		// load threads with work
		for (i = 0; i < THREADS; i++) 
		{
			nfs_threaddata_t *t = thread_data + i;
			t->job.startq = startq;
			startq += qrange;
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
		
#ifdef WIN32
			sprintf(syscmd,"type %s >> msieve.dat",t->outfilename);
#else
			sprintf(syscmd,"cat %s >> msieve.dat",t->outfilename);
#endif
			system(syscmd);
		}

		//set flags to do filtering
		flags = 0;
		flags = flags | MSIEVE_FLAG_USE_LOGFILE;
		if (VFLAG > 1)
			flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
		flags = flags | MSIEVE_FLAG_NFS_FILTER;
		obj->flags = flags;

		//msieve: filter
		if (VFLAG >= 0)
			printf("nfs: msieve filtering\n");
		relations_needed = nfs_filter_relations(obj, &mpN);
	}

	//stop worker threads
	for (i=0; i<THREADS - 1; i++)
	{
		nfs_stop_worker_thread(thread_data + i, 0);
	}
	nfs_stop_worker_thread(thread_data + i, 1);

	//don't need the threads anymore
	free(thread_data);

	//msieve: build matrix
	flags = 0;
	flags = flags | MSIEVE_FLAG_USE_LOGFILE;
	if (VFLAG > 1)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_LA;
	obj->flags = flags;

	if (VFLAG >= 0)
		printf("nfs: msieve linear algebra\n");
	nfs_solve_linear_system(obj, &mpN);			

	//msieve: find factors
	flags = 0;
	flags = flags | MSIEVE_FLAG_USE_LOGFILE;
	if (VFLAG > 1)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
	flags = flags | MSIEVE_FLAG_NFS_SQRT;
	obj->flags = flags;

	if (VFLAG >= 0)
		printf("nfs: msieve sqrt\n");
	nfs_find_factors(obj, &mpN, &factor_list);
	
	extract_factors(&factor_list,fobj);
	
	//start a counter for the whole job
	gettimeofday(&stop, NULL);

	difference = my_difftime (&start, &stop);

	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);	

	if (VFLAG >= 0)
		printf("NFS elapsed time = %6.4f seconds.\n",t_time);

	// free stuff
	msieve_obj_free(obj);
	sFree(&input_str);

	return;
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

	//matrix_thread_free(t);

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
	char syscmd[1024];	

	//remove any temporary relation files
	remove(thread_data->outfilename);
		
	//start ggnfs binary
	sprintf(syscmd,"%s -a ggnfs.job -f %u -c %u -o %s",
		thread_data->job.sievername, thread_data->job.startq, 
		thread_data->job.qrange, thread_data->outfilename);

	if (VFLAG >= 0)
	{
		printf("nfs: sieving with command: %s\n",syscmd);
	}
	system(syscmd);

	return 0;
}

void extract_factors(factor_list_t *factor_list, fact_obj_t *fobj)
{
	z tmp1, tmp2, *N = &fobj->qs_obj.n;
	int i;

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
		}
		else
		{
			tmpz.type = COMPOSITE;
			add_to_factor_list(fobj, &tmpz);
		}

		zFree(&tmpz);		
		//free(factor_list.final_factors[i]);
		zCopy(N,&tmp1);
	}

	if (zCompare(&tmp1,&zOne) > 0)
	{
		if (isPrime(&tmp1))
		{
			tmp1.type = PRP;
			add_to_factor_list(fobj, &tmp1);
			zCopy(&zOne, N);
		}
		else
		{
			tmp1.type = COMPOSITE;
			add_to_factor_list(fobj, &tmp1);
			zCopy(&zOne, N);
		}
	}

	zFree(&tmp1);
	zFree(&tmp2);

	return;
}


int check_existing_files(z *N, uint32 *last_spq)
{
	//check to see if there is a ggnfs.job for this number
	FILE *in;
	char line[GSTR_MAXSIZE], *ptr;
	int ans;
	z tmpz;


	in = fopen("ggnfs.job","r");
	if (in == NULL)
		return 0;

	ptr = fgets(line,GSTR_MAXSIZE,in);
	if (ptr == NULL)
		return 0;

	fclose(in);

	if (line[0] != 'n')
		return 0;

	zInit(&tmpz);
	z2decstr(N,&gstr1);
	str2hexz(line+3,&tmpz);
	ans = zCompare(&tmpz,N);
	if (ans == 0)
	{
		ans = 1;

		// attempt to open data file
		in = fopen("msieve.dat","r");
		if (in == NULL)
		{
			// didn't exist, return flag to start sieving from the beginning
			*last_spq = 0;
		}
		else
		{
			// found it.  read the last specialq
			char syscmd[1024];
			FILE *data;
			int a;

			fclose(in);
#ifdef WIN32
			a = system("tail msieve.dat > rels.dat");
			// no good tail command native to windows.  roll something similar here which will
			// likely be much slower
			// first, test for the existence of unxutils						
			if (a)
			{
				printf("for optimal performance, consider installing unix utilities for windows:\n");
				printf("http://unxutils.sourceforge.net/");

				in = fopen("msieve.dat","r");
				//set the file position indicated several kilobytes back from the end of the file
				//if this returns non-zero, then this was forbidden, likely because the file wasn't 
				//that big.  If that's the case, give up finding the last special q.
				a = fseek(in, -4096, SEEK_END);
				if (a)
				{
					fclose(in);
					*last_spq = 0;
					return ans;
				}
				else
				{
					// we were able to set the file position indicator back a ways.  we don't know
					// where it is exactly, so get and throw away a line
					ptr = fgets(line, GSTR_MAXSIZE, in);
					if (ptr == NULL)
					{
						fclose(in);
						*last_spq = 0;
						return ans;
					}

					// then get a line which we'll try to extract last_spq from
					ptr = fgets(line, GSTR_MAXSIZE, in);
					if (ptr == NULL)
					{
						fclose(in);
						*last_spq = 0;
						return ans;
					}
					else
					{
						int i;
						*last_spq = 0;
						//note: this assumes algebraic side sieving
						for (i=strlen(line); i>=0; i--)
						{
							if (line[i] == ',')
								break;
						}
						sscanf(line+i+1,"%x",last_spq);
						//printf("read last specialq value of %s\n",line+i+1);				
						fclose(in);
						return ans;
					}
				}
			}
			
			// if we got to here, tail worked: unxutils or something similar must be present.
			// share the linux code below.
#else
			system("tail msieve.dat > rels.dat");
#endif
			
			data = fopen("rels.dat","r");
			ptr = fgets(line,GSTR_MAXSIZE,data);
			if (ptr == NULL)
				*last_spq = 0;
			else
			{
				int i;
				*last_spq = 0;
				//note: this assumes algebraic side sieving
				for (i=strlen(line); i>=0; i--)
				{
					if (line[i] == ',')
						break;
				}
				sscanf(line+i+1,"%x",last_spq);
				//printf("read last specialq value of %s\n",line+i+1);								
			}
			fclose(data);
		}

	}
	else
	{
		ans = 0;
		*last_spq = 0;
	}

	zFree(&tmpz);
	return ans;

}

void find_best_msieve_poly(z *N, ggnfs_job_t *job)
{
	// parse a msieve.dat.p file to find the best polynomial (based on e score)
	// output this as a ggnfs polynomial file
	FILE *in, *out;
	char line[GSTR_MAXSIZE], *ptr;
	double score, bestscore = 0;
	int count, bestline, i;

	in = fopen("msieve.dat.p","r");
	if (in == NULL)
	{
		printf("could not open msieve.dat.p for reading!\n");
		exit(1);
	}

	// read and count polys of the file
	count = 0;
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;		

		if (line[0] == '#')
		{
			count++;
			//printf("found poly: %s",line);
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
					//printf("new best score = %e, new best line = %d\n",bestscore, bestline);
				}
			}
		}
		
	}
	fclose(in);

	//open it again
	in = fopen("msieve.dat.p","r");
	if (in == NULL)
	{
		printf("could not open msieve.dat.p for reading!\n");
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen("ggnfs.job","w");
	if (out == NULL)
	{
		printf("could not open ggnfs.job for writing!\n");
		exit(1);
	}

	//go to the bestline
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;		

		if (line[0] == '#')
			count--;

		if (count == 0)
			break;
	}

	//copy n into the job file
	z2decstr(N,&gstr1);
	fprintf(out, "n: %s\n",gstr1.s);

	// copy out the poly
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;

		if (line[0] == '#')
			break;		
		
		fputs(line,out);
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
	FILE *in, *out;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;

	in = fopen("msieve.fb","r");
	if (in == NULL)
	{
		printf("could not open msieve.fb for reading!\n");
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen("ggnfs.job","w");
	if (out == NULL)
	{
		printf("could not open ggnfs.job for writing!\n");
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

//digits, r/alim, lpbr/a, mfbr/a, r/alambda, siever, rels
//entries based on statistics gathered from many factorizations done
//over the years by myself and others, and from here:
//http://www.mersenneforum.org/showthread.php?t=12365
#define GGNFS_TABLE_ROWS 15
static double ggnfs_table[GGNFS_TABLE_ROWS][8] = {
	{85,  900000,   24, 48, 2.5, 11, 0, 10000},
	{90,  1200000,  25, 50, 2.5, 11, 0, 20000},
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
	FILE *test;
	
	job->rels = 0;
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
			job->lambda = 2.5;

			//pick closest entry
			if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
				job->siever = ggnfs_table[i][5];
			else
				job->siever = ggnfs_table[i+1][5];

			job->rels = 0.8 * 2 * pow(2.0,(double)job->lpb) / log(pow(2.0,(double)job->lpb));

			//interp
			job->qrange = ggnfs_table[i+1][7] - 
				(uint32)(scale * (double)(ggnfs_table[i+1][7] - ggnfs_table[i][7]));
		}
	}

	if (job->rels == 0)
	{
		//couldn't find a table entry
		if (d <= ggnfs_table[0][0])
		{
			job->fblim = ggnfs_table[0][1];
			job->lpb = ggnfs_table[0][2];
			job->mfb = ggnfs_table[0][3];
			job->lambda = ggnfs_table[0][4];
			job->siever = ggnfs_table[0][5];
			job->rels = 0.8 * 2 * pow(2.0,(double)job->lpb) / log(pow(2.0,(double)job->lpb));
			job->qrange = ggnfs_table[0][7];
		}
		else
		{
			job->fblim = ggnfs_table[GGNFS_TABLE_ROWS-1][1];
			job->lpb = ggnfs_table[GGNFS_TABLE_ROWS-1][2];
			job->mfb = ggnfs_table[GGNFS_TABLE_ROWS-1][3];
			job->lambda = ggnfs_table[GGNFS_TABLE_ROWS-1][4];
			job->siever = ggnfs_table[GGNFS_TABLE_ROWS-1][5];
			job->rels = 0.8 * 2 * pow(2.0,(double)job->lpb) / log(pow(2.0,(double)job->lpb));
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

#ifdef WIN32
	sprintf(job->sievername, "%s.exe", job->sievername);
#endif

	// test for existence of the siever
	test = fopen(job->sievername, "rb");
	if (test == NULL)
	{
		printf("could not find %s, bailing\n",job->sievername);
		exit(-1);
	}

	return;
}
	



