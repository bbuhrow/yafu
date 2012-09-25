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
       				   --bbuhrow@gmail.com 11/2/10
----------------------------------------------------------------------*/

#include "yafu_ecm.h"
#include "factor.h"
#include "yafu.h"
#include "calc.h"
#include "yafu_string.h"

int ecm_loop(fact_obj_t *fobj)
{
	//expects the input in ecm_obj->gmp_n
	ecm_thread_data_t *thread_data;		//an array of thread data objects
	mpz_t d,t;
	FILE *flog;
	int i,j;
	double est_time;
	double batch_time;
	double avg_batch_time;
	int num_batches;

	struct timeval stop;	
	struct timeval start;	
	TIME_DIFF *	difference;
	double t_time;

	//maybe make this an input option: whether or not to stop after
	//finding a factor in the middle of running a requested batch of curves
	int total_curves_run;
	int bail_on_factor = 1;
	int bail = 0;
	int input_digits = gmp_base10(fobj->ecm_obj.gmp_n);

	if (ecm_check_input(fobj) == 0)
		return 0;

	//ok, having gotten this far we are now ready to run the requested
	//curves.  initialize the needed data structures, then split
	//the curves up over N threads.  The main thread will farm out different
	//sigmas to the worker threads.

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	ECM_ABORT = 0;
	signal(SIGINT,ecmexit);

	//init ecm process
	ecm_process_init(fobj);

	//init local big ints
	mpz_init(d);
	mpz_init(t);

	thread_data = (ecm_thread_data_t *)malloc(THREADS * sizeof(ecm_thread_data_t));
	for (i=0; i<THREADS; i++)
	{
		thread_data[i].fobj = fobj;
		thread_data[i].thread_num = i;
		ecm_thread_init(&thread_data[i]);
	}

	//round numcurves up so that each thread has something to do every iteration.
	//this prevents a single threaded "cleanup" round after the multi-threaded rounds
	//to finish the leftover requested curves.
	fobj->ecm_obj.num_curves += ((fobj->ecm_obj.num_curves % THREADS) ? 
		(THREADS - (fobj->ecm_obj.num_curves % THREADS)) : 0);

	if (VFLAG >= 0)
	{
		for (i=0, total_curves_run=0; i<THREADS; i++)
            total_curves_run += thread_data[i].curves_run;

		printf("ecm: %d/%d curves on C%d, ",
			total_curves_run, fobj->ecm_obj.num_curves, 
			(int)gmp_base10(fobj->ecm_obj.gmp_n));
		print_B1B2(fobj, NULL);
		printf("\r");
		fflush(stdout);
	}

	/* activate the threads one at a time. The last is the
	   master thread (i.e. not a thread at all). */

	for (i = 0; i < THREADS - 1; i++)
		ecm_start_worker_thread(thread_data + i, 0);

	ecm_start_worker_thread(thread_data + i, 1);

	//split the requested curves up among the specified number of threads. 
	num_batches = 0;
	t_time = 0.0;
	for (j=0; j < fobj->ecm_obj.num_curves / THREADS; j++)
	{
		//watch for an abort
		if (ECM_ABORT)
		{
			print_factors(fobj);
			exit(1);
		}

		// start a counter for this batch of curves
		// for larger B1s
		if (fobj->ecm_obj.B1 > 48000)
			gettimeofday(&start, NULL);

		//do work on different sigmas
		for (i=0; i<THREADS; i++)
		{
			ecm_get_sigma(&thread_data[i]);

			if (i == THREADS - 1) {
				ecm_do_one_curve(&thread_data[i]);
			}
			else {
				thread_data[i].command = ECM_COMMAND_RUN;
#if defined(WIN32) || defined(_WIN64)
				SetEvent(thread_data[i].run_event);
#else
				pthread_cond_signal(&thread_data[i].run_cond);
				pthread_mutex_unlock(&thread_data[i].run_lock);
#endif
			}
		}

		//wait for threads to finish
		for (i=0; i<THREADS; i++)
		{
			if (i < THREADS - 1) {
#if defined(WIN32) || defined(_WIN64)
				WaitForSingleObject(thread_data[i].finish_event, INFINITE);
#else
				pthread_mutex_lock(&thread_data[i].run_lock);
				while (thread_data[i].command != ECM_COMMAND_WAIT)
					pthread_cond_wait(&thread_data[i].run_cond, &thread_data[i].run_lock);
#endif
			}
		}

		for (i=0; i<THREADS; i++)
		{
			//look at the result of each curve and see if we're done
			if ((mpz_cmp_ui(thread_data[i].gmp_factor, 1) > 0)
				&& (mpz_cmp(thread_data[i].gmp_factor, fobj->ecm_obj.gmp_n) < 0))
			{						
				//non-trivial factor found
				//since we could be doing many curves in parallel,
				//it's possible more than one curve found this factor at the
				//same time.  We divide out factors as we find them, so
				//just check if this factor still divides n
				mpz_tdiv_qr(t, d, fobj->ecm_obj.gmp_n, thread_data[i].gmp_factor);
				if (mpz_cmp_ui(d, 0) == 0)
				{
					//yes, it does... proceed to record the factor
					mpz_set(fobj->ecm_obj.gmp_n, t);
					ecm_deal_with_factor(&thread_data[i]);

					//we found a factor and might want to stop,
					//but should check all threads output first
					if (bail_on_factor)
						bail = 1;
					else if (mpz_probab_prime_p(fobj->ecm_obj.gmp_n, NUM_WITNESSES))
						bail = 1;
					else
					{
						//found a factor and the cofactor is composite.
						//the user has specified to keep going with ECM until the 
						//curve counts are finished thus:
						//we need to re-initialize with a different modulus.  this is
						//independant of the thread data initialization
						ecm_process_free(fobj);
						ecm_process_init(fobj);
					}
				}
			}

			thread_data[i].curves_run++;
		}

		if (bail)
			goto done;

		if (VFLAG >= 0)
		{
			for (i=0, total_curves_run=0; i<THREADS; i++)
				total_curves_run += thread_data[i].curves_run;			

			printf("ecm: %d/%d curves on C%d, ",
				total_curves_run, fobj->ecm_obj.num_curves, 
				(int)gmp_base10(fobj->ecm_obj.gmp_n));

			print_B1B2(fobj, NULL);
			
			// stop counter for this batch of curves
			// for larger B1s
			if (fobj->ecm_obj.B1 > 48000)
			{
				gettimeofday(&stop, NULL);
				difference = my_difftime (&start, &stop);
				batch_time = ((double)difference->secs + (double)difference->usecs / 1000000);
				free(difference);
				num_batches++;
				t_time += batch_time;
				avg_batch_time = t_time / (double)num_batches;
				est_time = (double)(fobj->ecm_obj.num_curves / THREADS - j) * avg_batch_time;

				printf(", ETA: %1.0f sec ", est_time);
			}
			printf("\r");
			fflush(stdout);
		}

	}

done:
	if (VFLAG >= 0)
		printf("\n");

	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for appending\n",fobj->flogname);
		return 0;
	}

	for (i=0, total_curves_run=0; i<THREADS; i++)
		total_curves_run += thread_data[i].curves_run;

	logprint(flog,"Finished %d curves using Lenstra ECM method on C%d input, ",
		total_curves_run,input_digits);

	print_B1B2(fobj, flog);
	fprintf(flog, "\n");

	fclose(flog);

	//stop worker threads
	for (i=0; i<THREADS - 1; i++)
	{
		ecm_stop_worker_thread(thread_data + i, 0);
	}
	ecm_stop_worker_thread(thread_data + i, 1);

	for (i=0; i<THREADS; i++)
		ecm_thread_free(&thread_data[i]);
	free(thread_data);

	mpz_clear(d);
	mpz_clear(t);
	signal(SIGINT,NULL);
	
	ecm_process_free(fobj);

	return total_curves_run;
}

int ecm_deal_with_factor(ecm_thread_data_t *thread_data)
{
	FILE *flog;
	fact_obj_t *fobj = thread_data->fobj;
	int curves_run = thread_data->curves_run;
	int thread_num = thread_data->thread_num;

	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for appending\n",fobj->flogname);
		return 0;
	}

	if (mpz_probab_prime_p(thread_data->gmp_factor, NUM_WITNESSES))
	{
		add_to_factor_list(fobj, thread_data->gmp_factor);

		if (VFLAG > 0)
			gmp_printf("\necm: found prp%d factor = %Zd\n", 
			gmp_base10(thread_data->gmp_factor),
			thread_data->gmp_factor);

		logprint(flog,"prp%d = %s (curve %d stg%d B1=%u sigma=%u thread=%d)\n",
			gmp_base10(thread_data->gmp_factor),
			mpz_conv2str(&gstr1.s, 10, thread_data->gmp_factor),
			curves_run+1, thread_data->stagefound,
			fobj->ecm_obj.B1, thread_data->sigma, thread_num);
	}
	else
	{
		add_to_factor_list(fobj, thread_data->gmp_factor);
		
		if (VFLAG > 0)
			gmp_printf("\necm: found c%d factor = %Zd\n", 
			gmp_base10(thread_data->gmp_factor),
			thread_data->gmp_factor);

		logprint(flog,"c%d = %s (curve %d stg%d B1=%u sigma=%u thread=%d)\n",
			gmp_base10(thread_data->gmp_factor),
			mpz_conv2str(&gstr1.s, 10, thread_data->gmp_factor),
			curves_run+1, thread_data->stagefound,
			fobj->ecm_obj.B1, thread_data->sigma, thread_num);
	}

	fclose(flog);

	return 1;
}

int ecm_get_sigma(ecm_thread_data_t *thread_data)
{
	fact_obj_t *fobj = thread_data->fobj;
	z tmp;

	zInit(&tmp);
	if (fobj->ecm_obj.sigma != 0)
	{
		if (thread_data->curves_run == 0 &&  fobj->ecm_obj.num_curves > 1)
			printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
		thread_data->sigma = fobj->ecm_obj.sigma;
	}
	else if (get_uvar("sigma",&tmp))
	{
		thread_data->sigma = spRand(6,MAX_DIGIT);
	}
	else
	{
		if (thread_data->curves_run == 0 &&  fobj->ecm_obj.num_curves > 1)
			printf("WARNING: work will be duplicated with sigma fixed and numcurves > 1\n");
		thread_data->sigma = (uint32)tmp.val[0];
	}

	zFree(&tmp);
	return 0;
}

int ecm_check_input(fact_obj_t *fobj)
{
	FILE *flog;

	//open the log file
	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for appending\n",fobj->flogname);
		return 0;
	}

	//check for trivial cases
	if (mpz_cmp_ui(fobj->ecm_obj.gmp_n, 0) == 0)
	{
		logprint(flog,"Trivial input == 0 in ECM\n");
		fclose(flog);
		return 0;
	}

	if (mpz_cmp_ui(fobj->ecm_obj.gmp_n, 1) == 0)
	{
		logprint(flog,"Trivial input == 1 in ECM\n");
		fclose(flog);
		return 0;
	}

	if (mpz_tdiv_ui(fobj->ecm_obj.gmp_n, 3) == 0)
	{
		mpz_t tmp;
		mpz_init(tmp);
		mpz_set_ui(tmp, 3);
		mpz_tdiv_q_ui(fobj->ecm_obj.gmp_n, fobj->ecm_obj.gmp_n, 3);
		add_to_factor_list(fobj, tmp);
		logprint(flog,"Trivial factor of 3 found in ECM\n");
		fclose(flog);
		mpz_clear(tmp);
		return 0;
	}

	if (mpz_tdiv_ui(fobj->ecm_obj.gmp_n, 2) == 0)
	{
		mpz_t tmp;
		mpz_init(tmp);
		mpz_set_ui(tmp, 2);
		mpz_tdiv_q_ui(fobj->ecm_obj.gmp_n, fobj->ecm_obj.gmp_n, 2);
		add_to_factor_list(fobj, tmp);
		logprint(flog,"Trivial factor of 2 found in ECM\n");
		fclose(flog);
		mpz_clear(tmp);
		return 0;
	}

	if (mpz_probab_prime_p(fobj->ecm_obj.gmp_n, NUM_WITNESSES))
	{
		//maybe have an input flag to optionally not perform
		//PRP testing (useful for really big inputs)
		add_to_factor_list(fobj, fobj->ecm_obj.gmp_n);
		logprint(flog,"prp%d = %s\n", gmp_base10(fobj->ecm_obj.gmp_n), 
			mpz_conv2str(&gstr1.s, 10, fobj->ecm_obj.gmp_n));		
		mpz_set_ui(fobj->ecm_obj.gmp_n, 1);
		fclose(flog);
		return 0;
	}

	//close the log file for until we have something further to report
	fclose(flog);

	return 1;
}

void ecm_start_worker_thread(ecm_thread_data_t *t, uint32 is_master_thread) {

	/* create a thread that will process a polynomial The last poly does 
	   not get its own thread (the current thread handles it) */

	if (is_master_thread) {
		return;
	}

	t->command = ECM_COMMAND_INIT;
#if defined(WIN32) || defined(_WIN64)
	t->run_event = CreateEvent(NULL, FALSE, TRUE, NULL);
	t->finish_event = CreateEvent(NULL, FALSE, FALSE, NULL);
	t->thread_id = CreateThread(NULL, 0, ecm_worker_thread_main, t, 0, NULL);

	WaitForSingleObject(t->finish_event, INFINITE); /* wait for ready */
#else
	pthread_mutex_init(&t->run_lock, NULL);
	pthread_cond_init(&t->run_cond, NULL);

	pthread_cond_signal(&t->run_cond);
	pthread_mutex_unlock(&t->run_lock);
	pthread_create(&t->thread_id, NULL, ecm_worker_thread_main, t);

	pthread_mutex_lock(&t->run_lock); /* wait for ready */
	while (t->command != ECM_COMMAND_WAIT)
		pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
}

void ecm_stop_worker_thread(ecm_thread_data_t *t, uint32 is_master_thread)
{
	if (is_master_thread) {
		return;
	}

	t->command = ECM_COMMAND_END;
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
DWORD WINAPI ecm_worker_thread_main(LPVOID thread_data) {
#else
void *ecm_worker_thread_main(void *thread_data) {
#endif
	ecm_thread_data_t *t = (ecm_thread_data_t *)thread_data;

	while(1) {

		/* wait forever for work to do */
#if defined(WIN32) || defined(_WIN64)
		WaitForSingleObject(t->run_event, INFINITE);		
#else
		pthread_mutex_lock(&t->run_lock);
		while (t->command == ECM_COMMAND_WAIT) {
			pthread_cond_wait(&t->run_cond, &t->run_lock);
		}
#endif
		/* do work */

		if (t->command == ECM_COMMAND_RUN)
			ecm_do_one_curve(t);
		else if (t->command == ECM_COMMAND_END)
			break;

		/* signal completion */

		t->command = ECM_COMMAND_WAIT;
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


void ecm_process_init(fact_obj_t *fobj)
{
	//initialize things which all threads will need when using
	//GMP-ECM
	TMP_THREADS = THREADS;
	TMP_STG2_MAX = fobj->ecm_obj.B2;

	if (strcmp(fobj->ecm_obj.ecm_path, "") != 0)
		fobj->ecm_obj.use_external = 1;
	
	if (THREADS > 1)
	{
		if (!fobj->ecm_obj.use_external)
		{
			if (VFLAG >= 2)
				printf("GMP-ECM does not support multiple threads... running single threaded\n");
			THREADS = 1;
		}
		else
		{
			if (fobj->ecm_obj.B1 < 40000 || fobj->ecm_obj.num_curves == 1)
			{
				THREADS = 1;
				fobj->ecm_obj.use_external = 0;
			}
		}
	}
	else
	{
		if (fobj->ecm_obj.B1 < 40000)
			fobj->ecm_obj.use_external = 0;
	}

	return;
}

void ecm_thread_init(ecm_thread_data_t *tdata)
{
	//initialize everything for all threads using GMP-ECM
	mpz_init(tdata->gmp_n);
	mpz_init(tdata->gmp_factor);
	ecm_init(tdata->params);
	gmp_randseed_ui(tdata->params->rng, get_rand(&g_rand.low, &g_rand.hi));
	mpz_set(tdata->gmp_n, tdata->fobj->ecm_obj.gmp_n);
	tdata->params->method = ECM_ECM;
	tdata->curves_run = 0;
		
	return;
}

void ecm_thread_free(ecm_thread_data_t *tdata)
{
	ecm_clear(tdata->params);
	mpz_clear(tdata->gmp_n);
	mpz_clear(tdata->gmp_factor);

	if (tdata->fobj->ecm_obj.use_external)
	{
		// remove temp file specific to this thread num
		sprintf(tdata->tmp_output, "_yafu_ecm_tmp%d.out", tdata->thread_num);
		remove(tdata->tmp_output);
	}

	return;
}

void ecm_process_free(fact_obj_t *fobj)
{
	THREADS = TMP_THREADS;
	fobj->ecm_obj.B2 = TMP_STG2_MAX;
	return;
}

void *ecm_do_one_curve(void *ptr)
{
	//unpack the data structure and stuff inside it
	ecm_thread_data_t *thread_data = (ecm_thread_data_t *)ptr;
	fact_obj_t *fobj = thread_data->fobj;

	if (!fobj->ecm_obj.use_external)
	{
		int status;

		thread_data->params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
		if (VFLAG >= 3)
			thread_data->params->verbose = VFLAG - 2;		
		mpz_set_ui(thread_data->params->x, (unsigned long)0);
		mpz_set_ui(thread_data->params->sigma, thread_data->sigma);

		// set the thread local copy of n
		mpz_set(thread_data->gmp_n, fobj->ecm_obj.gmp_n);

		if (fobj->ecm_obj.stg2_is_default == 0)
		{
			//not default, tell gmp-ecm to use the requested B2
			//printf("using requested B2 value\n");
			uint64_2gmp(fobj->ecm_obj.B2, thread_data->params->B2);
		}

		status = ecm_factor(thread_data->gmp_factor, thread_data->gmp_n,
				fobj->ecm_obj.B1, thread_data->params);		

		//printf ("used B2: ");
		//mpz_out_str (stdout, 10, thread_data->params->B2);
		//printf ("\n");

		//NOTE: this required a modification to the GMP-ECM source code in ecm.c
		//in order to get the automatically computed B2 value out of the
		//library
		//gmp2mp(thread_data->params->B2,&thread_data->factor);
		//ECM_STG2_MAX = z264(&thread_data->factor);

		//the return value is the stage the factor was found in, if no error
		thread_data->stagefound = status;
	}
	else
	{
		char *cmd;
		FILE *fid;
		char line[1024];
		char *ptr;
		char *tmpstr;
		int retcode;

		// let mpz figure out and allocate the string
		tmpstr = mpz_get_str(tmpstr, 10, thread_data->gmp_n);

		// allocate the appropriately sized command string
		cmd = (char *)malloc((strlen(tmpstr) + strlen(fobj->ecm_obj.ecm_path) + 100) 
			* sizeof(char));

		// external executable was specified
		sprintf(thread_data->tmp_output, "_yafu_ecm_tmp%d.out", thread_data->thread_num);

		// build system command
		//"echo \042 %s \042 | %s %u >> %s\n",
		sprintf(cmd, "echo %s | %s -sigma %u %u > %s\n", 
			tmpstr, fobj->ecm_obj.ecm_path, thread_data->sigma, fobj->ecm_obj.B1, 
			thread_data->tmp_output);

		// run system command
		retcode = system(cmd);

		free(tmpstr);
		free(cmd);

		// this is what I observed ecm returning on ctrl-c.  hopefully it is portable.
		if (retcode == 33280)
			ECM_ABORT = 1;
		
		// parse output file
		fid = fopen(thread_data->tmp_output, "r");
		while ((fid != NULL) && (!feof(fid)))
		{
			char fact[1024];

			fgets(line, 1024, fid);
			if (line == NULL)
				break;

			ptr = strstr(line, "**********");
			if (ptr == NULL)
				continue;

			// found a factor.  search for the :
			ptr = strstr(line, ":");
			if (ptr == NULL)
				continue;

			// the character prior to this is the stage, and the rest of the line
			// after it is the factor
			sscanf(ptr-2,"%d",&thread_data->stagefound);

			strcpy(fact, ptr+1);

			mpz_set_str(thread_data->gmp_factor, fact, 10);
			//str2hexz(fact, &thread_data->factor);

			break;
		}
		fclose(fid);

	}

	return 0;
}

// function definitions
void ecmexit(int sig)
{
	printf("\nAborting...\n");
	ECM_ABORT = 1;
	return;
}

int print_B1B2(fact_obj_t *fobj, FILE *fid)
{
	int i;
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (fobj->ecm_obj.B1 % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000000000, suffix);
	}
	else if (fobj->ecm_obj.B1 % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000000, suffix);
	}
	else if (fobj->ecm_obj.B1 % 1000 == 0)
	{
		suffix = 'K';
		sprintf(stg1str,"%u%c",fobj->ecm_obj.B1 / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",fobj->ecm_obj.B1);
	}

	if (fobj->ecm_obj.stg2_is_default == 0)
	{
		if (fobj->ecm_obj.B2 % 1000000000 == 0)
		{
			suffix = 'B';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000000000, suffix);
		}
		else if (fobj->ecm_obj.B2 % 1000000 == 0)
		{
			suffix = 'M';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000000, suffix);
		}
		else if (fobj->ecm_obj.B2 % 1000 == 0)
		{
			suffix = 'K';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->ecm_obj.B2 / 1000, suffix);
		}
		else
		{
			sprintf(stg2str,"%" PRIu64 "",fobj->ecm_obj.B2);
		}
	}
	else
		sprintf(stg2str, "gmp-ecm default");

	if (fid == NULL)
		i = printf("B1=%s, B2=%s",stg1str,stg2str);
	else
		i = fprintf(fid,"B1=%s, B2=%s",stg1str,stg2str);

	return i;
}



