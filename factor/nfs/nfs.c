#include "nfs.h"
#include "gmp_xface.h"

#ifdef USE_NFS

msieve_obj *obj_ptr;

void nfsexit(int sig)
{
	printf("\nReceived signal %d... please wait\n",sig);
	NFS_ABORT = 1;

	if (obj_ptr != NULL)
	{
		printf("setting flag\n");
		obj_ptr->flags |= MSIEVE_FLAG_STOP_SIEVING;
	}

	return;
}

//----------------------- NFS ENTRY POINT ------------------------------------//
void nfs(fact_obj_t *fobj)
{
	//expect the input in fobj->nfs_obj.gmp_n
	char *input;
	msieve_obj *obj = NULL;
	uint64 nfs_lower = 0;
	uint64 nfs_upper = 0;
	enum cpu_type cpu = yafu_get_cpu_type();
	mp_t mpN;
	factor_list_t factor_list;
	uint32 flags = 0;
	ggnfs_job_t job;
	uint32 relations_needed = 1;	
	uint32 last_specialq = 0;
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	struct timeval bstop;	// stop time of sieving batch
	struct timeval bstart;	// start time of sieving batch
	TIME_DIFF *	difference;
	double t_time;
	uint32 pre_batch_rels = 0;
	char tmpstr[GSTR_MAXSIZE];	
	int process_done;
	enum nfs_state_e nfs_state;

	obj_ptr = NULL;
	//below a certain amount, revert to SIQS
	if (gmp_base10(fobj->nfs_obj.gmp_n) < fobj->nfs_obj.min_digits)
	{
		mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
		SIQS(fobj);
		mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
		return;
	}	
		
	if (mpz_probab_prime_p(fobj->nfs_obj.gmp_n, NUM_WITNESSES))
	{
		add_to_factor_list(fobj, fobj->nfs_obj.gmp_n);
		
		if (VFLAG >= 0)
			gmp_printf("PRP%d = %Zd\n", gmp_base10(fobj->nfs_obj.gmp_n),
				fobj->nfs_obj.gmp_n);
		
		logprint_oc(fobj->flogname, "a", "PRP%d = %s\n",
			gmp_base10(fobj->nfs_obj.gmp_n),
			mpz_conv2str(&gstr1.s, 10, fobj->nfs_obj.gmp_n));	

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
		return;
	}

	if (mpz_perfect_square_p(fobj->nfs_obj.gmp_n))
	{
		mpz_sqrt(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n);

		add_to_factor_list(fobj, fobj->nfs_obj.gmp_n);
		logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
			gmp_base10(fobj->nfs_obj.gmp_n),
			mpz_conv2str(&gstr1.s, 10, fobj->nfs_obj.gmp_n));

		add_to_factor_list(fobj, fobj->nfs_obj.gmp_n);
		logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
			gmp_base10(fobj->nfs_obj.gmp_n),
			mpz_conv2str(&gstr1.s, 10, fobj->nfs_obj.gmp_n));

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
		return;
	}

	if (mpz_perfect_power_p(fobj->nfs_obj.gmp_n))
	{
		FILE *flog;
		uint32 j;

		if (VFLAG > 0)
			printf("input is a perfect power\n");
		
		factor_perfect_power(fobj, fobj->nfs_obj.gmp_n);

		flog = fopen(fobj->flogname, "a");
		if (flog == NULL)
		{
			printf("fopen error: %s\n", strerror(errno));
			printf("could not open %s for appending\n", fobj->flogname);
			exit(1);
		}
		
		logprint(flog,"input is a perfect power\n");

		for (j=0; j<fobj->num_factors; j++)
		{
			uint32 k;
			for (k=0; k<fobj->fobj_factors[j].count; k++)
			{
				logprint(flog,"prp%d = %s\n",gmp_base10(fobj->fobj_factors[j].factor), 
					mpz_conv2str(&gstr1.s, 10, fobj->fobj_factors[j].factor));
			}
		}

		fclose(flog);
		return;
	}

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	NFS_ABORT = 0;
	signal(SIGINT,nfsexit);

	//start a counter for the whole job
	gettimeofday(&start, NULL);
	
	//nfs state machine:
	input = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	nfs_state = NFS_STATE_INIT;
	process_done = 0;
	while (!process_done)
	{
		switch (nfs_state)
		{
		case NFS_STATE_INIT:
			// initialize some job parameters
			job.current_rels = 0;
			job.lpb = 0;
			job.mfb = 0;
			job.lambda = 0;
			job.type = 0;		// 0==GNFS, 1==SNFS
			job.size = 0;		// snfs difficulty
			job.fblim = 0;

			// write the input bigint as a string			
			input = mpz_conv2str(&input, 10, fobj->nfs_obj.gmp_n);

			// create an msieve_obj
			// this will initialize the savefile to the outputfile name provided
			obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
				fobj->nfs_obj.fbfile, g_rand.low, g_rand.hi, (uint32)0, nfs_lower, nfs_upper, cpu, 
				(uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);
			fobj->nfs_obj.mobj = obj;

			// initialize these before checking existing files.  If poly
			// select is resumed they will be changed by check_existing_files.
			job.last_leading_coeff = 0;
			job.poly_time = 0;

			// determine what to do next based on the state of various files.
			// this will set job.current_rels if it finds any
			nfs_state = check_existing_files(fobj, &last_specialq, &job);

			// get best parameters for the job - call this after checking the input
			// against the savefile because the input could be changed (i.e., reduced
			// by some factor).
			// get_ggnfs_params(fobj, &job);
			// instead call this only as necessary -- when done with poly select,
			// or when starting/resuming sieving

			// override the table-driven parameters with the supplied parameters,
			// if the user has supplied their own .job file
			//if (fobj->nfs_obj.user_job)
			//	parse_job_file(fobj, &job);

			// before we get started, check to make sure we can find ggnfs sievers
			// if we are going to be doing sieving
			if (check_for_sievers(fobj) == 1)
				nfs_state = NFS_STATE_DONE;

			if (VFLAG >= 0 && nfs_state != NFS_STATE_DONE)
				gmp_printf("nfs: commencing nfs on c%d: %Zd\n",
					gmp_base10(fobj->nfs_obj.gmp_n), fobj->nfs_obj.gmp_n);

			if (nfs_state != NFS_STATE_DONE)
				logprint_oc(fobj->flogname, "a", "nfs: commencing gnfs on c%d: %s\n",
					gmp_base10(fobj->nfs_obj.gmp_n),
					mpz_conv2str(&gstr1.s, 10, fobj->nfs_obj.gmp_n));

			// convert input to msieve bigint notation and initialize a list of factors
			gmp2mp_t(fobj->nfs_obj.gmp_n,&mpN);
			factor_list_init(&factor_list);
			
			if (fobj->nfs_obj.rangeq > 0)
				job.qrange = ceil((double)fobj->nfs_obj.rangeq / (double)THREADS);		

			break;

		case NFS_STATE_POLY:
			
			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_POLY))
				do_msieve_polyselect(fobj, obj, &job, &mpN, &factor_list);

			nfs_state = NFS_STATE_SIEVE;
			break;

		case NFS_STATE_SIEVE:

			pre_batch_rels = job.current_rels;
			gettimeofday(&bstart, NULL);

			if (((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE)) &&
				!(fobj->nfs_obj.nfs_phases & NFS_DONE_SIEVING))
				do_sieving(fobj, &job);

			// if this has been previously marked, then go ahead and exit.
			if (fobj->nfs_obj.nfs_phases & NFS_DONE_SIEVING)
				process_done = 1;
			
			// if user specified -ns with a fixed start and range,
			// then mark that we're done sieving.  
			if (fobj->nfs_obj.rangeq > 0)
			{
				// we're done sieving the requested range, but there may be
				// more phases to check, so don't exit yet
				fobj->nfs_obj.nfs_phases |= NFS_DONE_SIEVING;
			}
				
			// then move on to the next phase
			nfs_state = NFS_STATE_FILTCHECK;

			break;

		case NFS_STATE_FILTER:

			// if we've flagged not to do filtering, then assume we have
			// enough relations and move on to linear algebra
			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_FILTER))
				relations_needed = do_msieve_filtering(fobj, obj, &job);
			else
				relations_needed = 0;

			if (relations_needed == 0)
				nfs_state = NFS_STATE_LINALG;
			else
				nfs_state = NFS_STATE_SIEVE;

			break;

		case NFS_STATE_LINALG:

			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA_RESUME))
			{
				// msieve: build matrix
				flags = 0;
				flags = flags | MSIEVE_FLAG_USE_LOGFILE;
				if (VFLAG > 0)
					flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
				flags = flags | MSIEVE_FLAG_NFS_LA;

				// add restart flag if requested
				if (fobj->nfs_obj.nfs_phases & NFS_PHASE_LA_RESUME)
					flags |= MSIEVE_FLAG_NFS_LA_RESTART;

				obj->flags = flags;

				if (VFLAG >= 0)
					printf("nfs: commencing msieve linear algebra\n");

				logprint_oc(fobj->flogname, "a", "nfs: commencing msieve linear algebra\n");

				// use a different number of threads for the LA, if requested
				if (LATHREADS > 0)
				{
					msieve_obj_free(obj);
					obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
						fobj->nfs_obj.fbfile, g_rand.low, g_rand.hi, (uint32)0, nfs_lower, nfs_upper, cpu, 
						(uint32)L1CACHE, (uint32)L2CACHE, (uint32)LATHREADS, (uint32)0, (uint32)0, 0.0);
				}

				// try this hack - store a pointer to the msieve obj so that
				// we can change a flag on abort in order to interrupt the LA.
				obj_ptr = obj;
				nfs_solve_linear_system(obj, fobj->nfs_obj.gmp_n);
				if (obj_ptr->flags & MSIEVE_FLAG_STOP_SIEVING)
					nfs_state = NFS_STATE_DONE;
				else
					nfs_state = NFS_STATE_SQRT;

				// set the msieve threads back to where it was if we used
				// a different amount for linalg
				if (LATHREADS > 0)
				{
					msieve_obj_free(obj);
					obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
						fobj->nfs_obj.fbfile, g_rand.low, g_rand.hi, (uint32)0, nfs_lower, nfs_upper, cpu, 
						(uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);
				}
			
				obj_ptr = NULL;
			}
			else
				nfs_state = NFS_STATE_SQRT;

			break;

		case NFS_STATE_SQRT:

			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_SQRT))
			{
				// msieve: find factors
				flags = 0;
				flags = flags | MSIEVE_FLAG_USE_LOGFILE;
				if (VFLAG > 0)
					flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
				flags = flags | MSIEVE_FLAG_NFS_SQRT;
				obj->flags = flags;

				if (VFLAG >= 0)
					printf("nfs: commencing msieve sqrt\n");

				logprint_oc(fobj->flogname, "a", "nfs: commencing msieve sqrt\n");

				// try this hack - store a pointer to the msieve obj so that
				// we can change a flag on abort in order to interrupt the LA.
				obj_ptr = obj;

				nfs_find_factors(obj, fobj->nfs_obj.gmp_n, &factor_list);

				obj_ptr = NULL;
			
				extract_factors(&factor_list,fobj);

				if (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) == 0)
					nfs_state = NFS_STATE_CLEANUP;		//completely factored, clean up everything
				else
					nfs_state = NFS_STATE_DONE;		//not factored completely, keep files and stop
			}
			else
				nfs_state = NFS_STATE_DONE;		//not factored completely, keep files and stop

			break;

		case NFS_STATE_CLEANUP:
			
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
			sprintf(tmpstr, "%s.mat.chk",fobj->nfs_obj.outputfile);	remove(tmpstr);

			gettimeofday(&stop, NULL);

			difference = my_difftime (&start, &stop);

			t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
			free(difference);	

			if (VFLAG >= 0)
				printf("NFS elapsed time = %6.4f seconds.\n",t_time);

			logprint_oc(fobj->flogname, "a", "NFS elapsed time = %6.4f seconds.\n",t_time);
			logprint_oc(fobj->flogname, "a", "\n");
			logprint_oc(fobj->flogname, "a", "\n");

			// free stuff			
			nfs_state = NFS_STATE_DONE;
			break;

		case NFS_STATE_DONE:
			process_done = 1;
			break;

		case NFS_STATE_FILTCHECK:
			if (job.current_rels >= job.min_rels)
			{
				if (VFLAG > 0)
					printf("found %u relations, need at least %u, proceeding with filtering ...\n",
					job.current_rels, job.min_rels);
				
				nfs_state = NFS_STATE_FILTER;
			}
			else
			{
				// compute eta by dividing how many rels we have left to find
				// by the average time per relation.  we have average time
				// per relation because we've saved the time it took to do 
				// the last batch of sieving and we know how many relations we
				// found in that batch.
				uint32 est_time;

				gettimeofday(&bstop, NULL);
				difference = my_difftime (&bstart, &bstop);
				t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
				free(difference);

				est_time = (uint32)((job.min_rels - job.current_rels) * 
					(t_time / (job.current_rels - pre_batch_rels)));				

				// if the user doesn't want to sieve, then we can't make progress.
				if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
					(fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE))
				{
					if (VFLAG > 0)
						printf("found %u relations, need at least %u "
							"(filtering ETA: %uh %um), continuing with sieving ...\n",
							job.current_rels, job.min_rels, est_time / 3600, 
							(est_time % 3600) / 60);

					nfs_state = NFS_STATE_SIEVE;
				}
				else
				{
					if (VFLAG > 0)
						printf("found %u relations, need at least %u "
							"(filtering ETA: %uh %um), sieving not selected, finishing ...\n",
							job.current_rels, job.min_rels, est_time / 3600, 
							(est_time % 3600) / 60);

					nfs_state = NFS_STATE_DONE;
				}
			}
			break;

		case NFS_STATE_STARTNEW:

			//job.startq = job.fblim / 2;
			nfs_state = NFS_STATE_POLY;		

			// create a new directory for this job 
//#ifdef _WIN32
//			sprintf(tmpstr, "%s\%s", fobj->nfs_obj.ggnfs_dir, 
//				mpz_conv2str(&gstr1.s, 10, fobj->nfs_obj.gmp_n));
//			mkdir(tmpstr);
//#else
//			sprintf(tmpstr, "%s/%s", fobj->nfs_obj.ggnfs_dir, 
//				mpz_conv2str(&gstr1.s, 10, fobj->nfs_obj.gmp_n));
//			mkdir(tmpstr, S_IRWXU);
//#endif
			
			// point msieve and ggnfs to the new directory
//#ifdef _WIN32
//			sprintf(fobj->nfs_obj.outputfile, "%s\%s", 
//				tmpstr, fobj->nfs_obj.outputfile);
//			sprintf(fobj->nfs_obj.logfile, "%s\%s", 
//				tmpstr, fobj->nfs_obj.logfile);
//			sprintf(fobj->nfs_obj.fbfile, "%s\%s", 
//				tmpstr, fobj->nfs_obj.fbfile);
//#else
//			sprintf(fobj->nfs_obj.outputfile, "%s%s", 
//				tmpstr, fobj->nfs_obj.outputfile);
//			sprintf(fobj->nfs_obj.logfile, "%s%s", 
//				tmpstr, fobj->nfs_obj.logfile);
//			sprintf(fobj->nfs_obj.fbfile, "%s%s", 
//				tmpstr, fobj->nfs_obj.fbfile);
//
//#endif
//
//			msieve_obj_free(fobj->nfs_obj.mobj);
//			obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
//				fobj->nfs_obj.fbfile, g_rand.low, g_rand.hi, (uint32)0, nfs_lower, nfs_upper, cpu, 
//				(uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);
//			fobj->nfs_obj.mobj = obj;
//
//			printf("output: %s\n", fobj->nfs_obj.mobj->savefile.name);
//			printf("log: %s\n", fobj->nfs_obj.mobj->logfile_name);
//			printf("fb: %s\n", fobj->nfs_obj.mobj->nfs_fbfile_name);

			break;

		case NFS_STATE_RESUMESIEVE:
			if (last_specialq == 0) 
 			{				
				// no data file, so this is the first time using this job file
				uint32 missing_params = PARAM_FLAG_NONE;
				
				parse_job_file(fobj, &job, &missing_params);
				
				// set min_rels.  
				get_ggnfs_params(fobj, &job);
				
				fill_job_file(fobj, &job, missing_params);
				// if any ggnfs params are missing, fill them
				// this means the user can provide an SNFS poly or external GNFS poly, 
				// but let YAFU choose the other params
				// this won't override any params in the file.
			
				if (fobj->nfs_obj.startq > 0)
				{
					job.startq = fobj->nfs_obj.startq;
					if (VFLAG >= 0)
						printf("nfs: continuing with sieving at user specified special-q %u\n",
							job.startq);

					logprint_oc(fobj->flogname, "a", "nfs: continuing with sieving at user "
						"specified special-q %u\n", job.startq);
				}
				else
				{
					job.startq = job.fblim / 2;
					if (VFLAG >= 0)
						printf("nfs: continuing with sieving - could not determine "
							"last special q; using default startq\n");

					logprint_oc(fobj->flogname, "a", "nfs: continuing with sieving "
						"- could not determine last special q; using default startq\n");
				}

				// next step is sieving
				nfs_state = NFS_STATE_SIEVE;
			}
			else
			{				
				parse_job_file(fobj, &job, NULL);
				// set min_rels.  
				get_ggnfs_params(fobj, &job);

				if (fobj->nfs_obj.startq > 0)
				{
					job.startq = fobj->nfs_obj.startq;
					nfs_state = NFS_STATE_SIEVE;
				}
				else
				{
					job.startq = last_specialq;

					// since we apparently found some relations, and we done
					// have any special plan in place from the user,
					// see if we have enough to finish
					nfs_state = NFS_STATE_FILTCHECK;
				}

				if (VFLAG >= 0)
					printf("nfs: found %u relations, need at least %u, "
					"continuing job at specialq = %u\n",
					job.current_rels, job.min_rels, job.startq);

				logprint_oc(fobj->flogname, "a", "nfs: found %u relations, "
					"need at least %u, continuing job at specialq = %u\n",
					job.current_rels, job.min_rels, job.startq);

			}

			// if there is a job file and the user has specified -np, print
			// this warning.
			if (fobj->nfs_obj.nfs_phases & NFS_PHASE_POLY)
			{
				printf("WARNING: .job file exists!  If you really want to redo poly selection,"
					" delete the .job file.\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}

			break;

		case NFS_STATE_RESUMEPOLY:
			if (VFLAG > 1) printf("nfs: resuming poly select\n");
			fobj->nfs_obj.polystart = job.last_leading_coeff;

			//load the job structure with the starting Q.
			//job.startq = job.fblim / 2;

			nfs_state = NFS_STATE_POLY;

			break;

		default:
			printf("unknown state, bailing\n");
			// non ideal solution to infinite loop in factor() if we return without factors
			// (should return error code instead)
			NFS_ABORT = 1;
			break;

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

			logprint_oc(fobj->flogname, "a", "NFS timeout after %6.4f seconds.\n",t_time);
			process_done = 1;
		}

		if (NFS_ABORT)
		{
			print_factors(fobj);
			exit(1);
		}
	}

	//reset signal handler to default (no handler).
	signal(SIGINT,NULL);

	if (obj != NULL)
		msieve_obj_free(obj);
	free(input);

	return;
}

#else

void nfs(fact_obj_t *fobj)
{
	printf("gnfs has not been enabled\n");

	mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
	SIQS(fobj);
	mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);

	return;
}

#endif

int check_for_sievers(fact_obj_t *fobj)
{
	// if we are going to be doing sieving, check for the sievers
	if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
		(fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE))
	{
		FILE *test;
		char name[1024];
		int found, i;

		for (i=11; i<16; i++)
		{
			sprintf(name, "%sggnfs-lasieve4I%de", fobj->nfs_obj.ggnfs_dir, i);
			// test for existence of the siever
			test = fopen(name, "rb");
			if (test != NULL)
			{
				found = 1;
				fclose(test);
				break;
			}
		}

		if (!found)
		{
			printf("WARNING: could not find ggnfs sievers, reverting to siqs!\n");
			logprint_oc(fobj->flogname, "a", "WARNING: could not find ggnfs sievers, "
				"reverting to siqs!\n");

			mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
			SIQS(fobj);
			mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
			return 1;
		}
		else 
			return 0;
	}

	return 0;
}

//entries based on statistics gathered from many factorizations done
//over the years by myself and others, and from here:
//http://www.mersenneforum.org/showthread.php?t=12365
#define GGNFS_TABLE_ROWS 15
static double ggnfs_table[GGNFS_TABLE_ROWS][8] = {
/* note: min_rels column is no longer used - it is equation based and	*/
/* is filled in by get_ggnfs_params										*/
/* columns:																*/
/* digits, r/alim, lpbr/a, mfbr/a, r/alambda, siever, min-rels, q-range */
	{85,  900000,   24, 48, 2.1, 11, 0, 10000},
	{90,  1200000,  25, 50, 2.3, 11, 0, 20000},
	{95,  1500000,  25, 50, 2.5, 12, 0, 40000},
	{100, 1800000,  26, 52, 2.5, 12, 0, 40000},
	{105, 2500000,  26, 52, 2.5, 12, 0, 80000},
	{110, 3200000,  26, 52, 2.5, 13, 0, 80000},
	{115, 4500000,  27, 54, 2.5, 13, 0, 160000},
	{120, 5000000,  27, 54, 2.5, 13, 0, 160000},
	{125, 5500000,  27, 54, 2.5, 13, 0, 160000},
	{130, 6000000,  27, 54, 2.5, 13, 0, 320000},
	{135, 8000000,  27, 54, 2.5, 14, 0, 320000},
	{140, 12000000, 28, 56, 2.5, 14, 0, 320000},
	{145, 15000000, 28, 56, 2.5, 14, 0, 640000},
	{150, 20000000, 29, 58, 2.5, 14, 0, 640000},
	{155, 30000000, 29, 58, 2.5, 15, 0, 640000}
};


void get_ggnfs_params(fact_obj_t *fobj, ggnfs_job_t *job)
{
	// based on the size/difficulty of the input number, determine "good" parameters
	// for the following: factor base limits, factor base large prime bound, trial
	// division fudge factors, trial division cutoffs, ggnfs lattice siever area, 
	// and expected number of relations needed.  linearly interpolate between table
	// entries.  keep last valid entry off the ends of the table.  This will produce
	// increasingly poor choices as one goes farther off the table, but you should be
	// doing things by hand by then anyway.
	uint32 i, d = gmp_base10(fobj->nfs_obj.gmp_n);
	double scale;
	double fudge; // sundae :)
	int found = 0;
	uint32 lpb;

	//job->min_rels = 0;
	if (job->type == 0 && job->size != 0 && job->size != d && VFLAG > 0)
		printf("nfs: warning: size param in job file does not match size of "
			"number, ignoring param\n");
	
	if (job->type > 0) // SNFS
	{
		if (job->size == 0)
		{	
			if (VFLAG > 1)
				printf("nfs: detected snfs job but no snfs difficulty; "
					"assuming size of number is the snfs difficulty\n");
		}
		else
			d = job->size;

		// http://www.mersenneforum.org/showpost.php?p=312701&postcount=2
		i = 0.56*d + 30; 
		if (VFLAG > 0)
			printf("nfs: guessing snfs difficulty %d is roughly equal to "
				"gnfs difficulty %d\n", d, i);
		d = i;
		if (fobj->nfs_obj.sq_side == 0) // not overridden by the user
			fobj->nfs_obj.sq_side = -1; // rational side
	}

	for (i=0; i<GGNFS_TABLE_ROWS - 1; i++)
	{
		if (d > ggnfs_table[i][0] && d <= ggnfs_table[i+1][0])
		{
			scale = (double)(ggnfs_table[i+1][0] - d) /
				(double)(ggnfs_table[i+1][0] - ggnfs_table[i][0]);

			//interp
			//job->fblim = ggnfs_table[i+1][1] - 
			//	(uint32)(scale * (double)(ggnfs_table[i+1][1] - ggnfs_table[i][1]));
			if (job->fblim == 0)
				job->fblim = ggnfs_table[i+1][1] - 
					(uint32)(scale * (double)(ggnfs_table[i+1][1] - ggnfs_table[i][1]));

			//pick closest entry
			//if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
			//	job->lpb = ggnfs_table[i][2];
			//else
			//	job->lpb = ggnfs_table[i+1][2];
			if (job->lpb == 0)
			{
				if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
					job->lpb = ggnfs_table[i][2];
				else
					job->lpb = ggnfs_table[i+1][2];
			}
			
			//pick closest entry
			//if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
			//	job->mfb = ggnfs_table[i][3];
			//else
			//	job->mfb = ggnfs_table[i+1][3];
			if (job->mfb == 0)
			{
				//pick closest entry
				if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
					job->mfb = ggnfs_table[i][3];
				else
					job->mfb = ggnfs_table[i+1][3];
			}

			//pick closest entry
			//if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
			//	job->lambda = ggnfs_table[i][4];
			//else
			//	job->lambda = ggnfs_table[i+1][4];
			if (job->lambda == 0)
			{
				//pick closest entry
				if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
					job->lambda = ggnfs_table[i][4];
				else
					job->lambda = ggnfs_table[i+1][4];
			}

			//pick closest entry
			if ((d - ggnfs_table[i][0]) < (ggnfs_table[i+1][0] - d))
				job->siever = ggnfs_table[i][5];
			else
				job->siever = ggnfs_table[i+1][5];

			//interp
			job->qrange = ggnfs_table[i+1][7] - 
				(uint32)(scale * (double)(ggnfs_table[i+1][7] - ggnfs_table[i][7]));

			found = 1;
		}
	}

	if (found == 0)
	{
		//couldn't find a table entry
		if (d <= ggnfs_table[0][0])
		{
			//job->fblim = ggnfs_table[0][1];
			//job->lpb = ggnfs_table[0][2];
			//job->mfb = ggnfs_table[0][3];
			//job->lambda = ggnfs_table[0][4];
			if (job->fblim == 0) job->fblim = ggnfs_table[0][1];
			if (job->lpb == 0) job->lpb = ggnfs_table[0][2];
			if (job->mfb == 0) job->mfb = ggnfs_table[0][3];
			if (job->lambda == 0) job->lambda = ggnfs_table[0][4];
			job->siever = ggnfs_table[0][5];
			job->qrange = ggnfs_table[0][7];
		}
		else
		{
			//job->fblim = ggnfs_table[GGNFS_TABLE_ROWS-1][1];
			//job->lpb = ggnfs_table[GGNFS_TABLE_ROWS-1][2];
			//job->mfb = ggnfs_table[GGNFS_TABLE_ROWS-1][3];
			//job->lambda = ggnfs_table[GGNFS_TABLE_ROWS-1][4];
			if (job->fblim == 0) job->fblim = ggnfs_table[GGNFS_TABLE_ROWS-1][1];
			if (job->lpb == 0) job->lpb = ggnfs_table[GGNFS_TABLE_ROWS-1][2];
			if (job->mfb == 0) job->mfb = ggnfs_table[GGNFS_TABLE_ROWS-1][3];
			if (job->lambda == 0) job->lambda = ggnfs_table[GGNFS_TABLE_ROWS-1][4];
			job->siever = ggnfs_table[GGNFS_TABLE_ROWS-1][5];
			job->qrange = ggnfs_table[GGNFS_TABLE_ROWS-1][7];
		}
	}

	job->min_rels = 0;
	for (i = 0; i < 2; i++)
	{
		// these are always the same now... but someday maybe they won't be.  this
		// routine will handle that when we get there.
		if (i == 0)
			lpb = job->lpb;
		else
			lpb = job->lpb;

		// appoximate min_rels:
		// http://ggnfs.svn.sourceforge.net/viewvc/ggnfs/trunk/tests/factMsieve.pl?r1=374&r2=416
		// http://www.mersenneforum.org/showpost.php?p=294055&postcount=25
		switch (lpb)
		{
		case 24:
			fudge = 0.40;
			break;
		case 25:
			fudge = 0.50;
			break;
		case 26:
			fudge = 0.60;
			break;
		case 27:
			fudge = 0.70;
			break;
		case 28:
			fudge = 0.76;
			break;
		case 29:
			fudge = 0.84;
			break;
		case 30:
			fudge = 0.89;
			break;
		case 31:
			fudge = 0.91;
			break;
		case 32:
			fudge = 0.95;
			break;
		case 33:
			fudge = 0.98;
			break;
		default:
			fudge = 0.40;
			break;
		}

		job->min_rels += (uint32)(fudge * (
			pow(2.0,(double)lpb) / log(pow(2.0,(double)lpb))));
	}


	// if the user specified a siever, use that one, else use the one
	// determined from the lookup table
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

void parse_job_file(fact_obj_t *fobj, ggnfs_job_t *job, uint32* missing_params)
{
	FILE *in;
	//uint32 lpbr = 0, lpba = 0;
	uint32 lpbr = 0, lpba = 0, mfbr = 0, mfba = 0, alim = 0, rlim = 0;
	char line[1024];
	float alambda = 0, rlambda = 0;

	in = fopen(fobj->nfs_obj.job_infile, "r");
	if (in == NULL)
	{
		printf("nfs: couldn't open job file, using default min_rels\n");
		return;
	}

	while (!feof(in))
	{
		char *substr, *ptr;
		
		ptr = fgets(line, 1024, in);

		// bail if we couldn't read anything
		if (ptr == NULL)
			break;

		substr = strstr(line, "lpbr:");

		if (substr != NULL)
		{
			lpbr = strtoul(substr + 5, NULL, 10);
			continue;
		}

		substr = strstr(line, "lpba:");

		if (substr != NULL)
		{
			lpba = strtoul(substr + 5, NULL, 10);
			continue;
		}

		substr = strstr(line, "type:");		
		if (substr != NULL)
 		{
			job->type = (NULL != strstr(substr + 5, "snfs")); // case sensitive
			if (VFLAG > 0 && job->type > 0)
				printf("nfs: found type: snfs\n");
			continue;
		}
		
		substr = strstr(line, "size:");
		if (substr != NULL)
		{
			job->size = strtoul(substr + 5, NULL, 10);
			if (VFLAG > 0)
				printf("nfs: found size: %u\n", job->size);
			continue;
		}
	
		substr = strstr(line, "mfbr:");
		if (substr != NULL)
		{
			mfbr = strtoul(substr + 5, NULL, 10);
			continue;
		}
		
		substr = strstr(line, "mfba:");
		if (substr != NULL)
		{
			mfba = strtoul(substr + 5, NULL, 10);
			continue;
		}
		
		substr = strstr(line, "rlim:");		
		if (substr != NULL)
		{
			rlim = strtoul(substr + 5, NULL, 10);
			continue;
		}
		
		substr = strstr(line, "alim:");		
		if (substr != NULL)
		{
			alim = strtoul(substr + 5, NULL, 10);
			continue;
		}
		
		substr = strstr(line, "rlambda:");		
		if (substr != NULL)
		{
			sscanf(substr + 8, "%f", &rlambda); //strtof(substr + 8, NULL);
			continue;
		}
		
		substr = strstr(line, "alambda:");		
		if (substr != NULL)
		{
			sscanf(substr + 8, "%f", &alambda); //strtof(substr + 8, NULL);
			continue;
		}
	}

	if ((lpba > 0) && (lpbr > 0))
	{
		if (job->type) // assume -a for gnfs, -r for snfs
			job->lpb = lpbr;
		else
			job->lpb = lpba;

		if (VFLAG > 0)
			printf("nfs: parsed lpbr = %u, lpba = %u, type = %s, setting lpb = %u\n",
				lpbr, lpba, job->type ? "snfs" : "gnfs", job->lpb);
	}
	else if ((lpba == 0) && (lpbr == 0))
	{
		if (missing_params) *missing_params |= PARAM_FLAG_LPB;
	}
	
	if ((mfba > 0) && (mfbr > 0))
	{
		if (job->type) // assume -a for gnfs, -r for snfs
			job->mfb = mfbr;
		else
			job->mfb = mfba;
		
		if (VFLAG > 2)
			printf("nfs: parsed mfbr = %u, mfba = %u, setting mfb = %u\n",
				mfbr, mfba, job->mfb);
	}
	else if ((mfba == 0) && (mfbr == 0))
	{
		if (missing_params) *missing_params |= PARAM_FLAG_MFB;
	}
	
	if ((alim > 0) && (rlim > 0))
	{
		if (job->type) // assume -a for gnfs, -r for snfs
			job->fblim = rlim;
		else
			job->fblim = alim;
		
		if (VFLAG > 2)
			printf("nfs: parsed rlim = %u, alim = %u, setting fblim = %u\n",
				rlim, alim, job->fblim);
	}
	else if ((alim == 0) && (rlim == 0))
	{
		if (missing_params) *missing_params |= PARAM_FLAG_FBLIM;
	}
	
	if ((alambda > 0) && (rlambda > 0))
	{
		if (job->type) // assume -a for gnfs, -r for snfs
			job->lambda = rlambda;
		else
			job->lambda = alambda;
		
		if (VFLAG > 2)
			printf("nfs: parsed rlim = %.1f, alim = %.1f, setting lambda = %.1f\n",
				rlambda, alambda, job->lambda);
	}
	else if ((alambda == 0) && (rlambda == 0))
	{
		if (missing_params) *missing_params |= PARAM_FLAG_LAMBDA;
	}


	/*
		if ((lpba > 0) && (lpbr > 0))
		{
			uint32 lpb;
			int i;
			job->min_rels = 0;

			for (i = 0; i < 2; i++)
			{
				if (i == 0)
					lpb = lpbr;
				else
					lpb = lpba;

				// appoximate min_rels:
				// http://ggnfs.svn.sourceforge.net/viewvc/ggnfs/trunk/tests/factMsieve.pl?r1=374&r2=416
				// http://www.mersenneforum.org/showpost.php?p=294055&postcount=25
				switch (lpb)
				{
				case 24:
					fudge = 0.35;
					break;
				case 25:
					fudge = 0.50;
					break;
				case 26:
					fudge = 0.60;
					break;
				case 27:
					fudge = 0.70;
					break;
				case 28:
					fudge = 0.76;
					break;
				case 29:
					fudge = 0.84;
					break;
				case 30:
					fudge = 0.89;
					break;
				case 31:
					fudge = 0.91;
					break;
				case 32:
					fudge = 0.95;
					break;
				case 33:
					fudge = 0.98;
					break;
				default:
					fudge = 0.2;
					break;
				}

				job->min_rels += (uint32)(fudge * (
					pow(2.0,(double)lpb) / log(pow(2.0,(double)lpb))));
			}

			printf("nfs: parsed lpbr = %u, lpba = %u; setting min_rels = %u\n",
				lpbr, lpba, job->min_rels);
			break;
		}			
	}
	*/


	fclose(in);

	return;
}
	
void fill_job_file(fact_obj_t *fobj, ggnfs_job_t *job, uint32 missing_params)
{
	if (missing_params != PARAM_FLAG_NONE)
	{
		//printf("Missing params: %d\n", missing_params);
		FILE* out = fopen(fobj->nfs_obj.job_infile, "a");
		if (out == NULL)
		{
			printf("nfs: couldn't fill job file, will try sieving anyway\n");
			return;
		}
		else if (VFLAG > 0)
			printf("nfs: job file is missing params, filling them\n");

		// make sure we start on a new line
		fprintf(out, "\n");

		if (missing_params & PARAM_FLAG_FBLIM)
		{
			fprintf(out,"rlim: %u\n",job->fblim);
			fprintf(out,"alim: %u\n",job->fblim);
		}
		
		if (missing_params & PARAM_FLAG_LPB)
		{
			fprintf(out,"lpbr: %u\n",job->lpb);
			fprintf(out,"lpba: %u\n",job->lpb);
		}
		
		if (missing_params & PARAM_FLAG_MFB)
		{
			fprintf(out,"mfbr: %u\n",job->mfb);
			fprintf(out,"mfba: %u\n",job->mfb);
		}
		
		if (missing_params & PARAM_FLAG_LAMBDA)
		{
			fprintf(out,"rlambda: %.1f\n",job->lambda);
			fprintf(out,"alambda: %.1f\n",job->lambda);
		}
		
		fclose(out);
	}
 	
	return;
}
