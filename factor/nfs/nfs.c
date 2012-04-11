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
	uint32 startq, qrange = 1;
	int is_continuation;
	uint32 last_specialq = 0;
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	TIME_DIFF *	difference;
	double t_time;
	//FILE *logfile;
	int statenum;
	char tmpstr[GSTR_MAXSIZE];	
	int process_done;

	obj_ptr = NULL;
	//below a certain amount, revert to SIQS
	if (gmp_base10(fobj->nfs_obj.gmp_n) < fobj->nfs_obj.min_digits)
	{
		mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
		SIQS(fobj);
		mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
		return;
	}	
		
	// remove this - it messes with detecting whether or not 
	// this is a restart of a job.
	////first a small amount of trial division
	////which will add anything found to the global factor list
	////this is only called with the main thread
	//if (VFLAG > 0)
	//	printf("nfs: commencing trial factoring\n");
	//mpz_set(fobj->div_obj.gmp_n, fobj->nfs_obj.gmp_n);
	//fobj->div_obj.print = 0;
	//fobj->div_obj.limit = 10000;
	//zTrial(fobj);
	//mpz_set(fobj->nfs_obj.gmp_n, fobj->div_obj.gmp_n);

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

	//check for inconsistent options
	if (((fobj->nfs_obj.startq > 0) && !(fobj->nfs_obj.rangeq > 0)) ||
		(!(fobj->nfs_obj.startq > 0) && (fobj->nfs_obj.rangeq > 0)))
	{
		printf("-f and -c must be specified together\n");
		print_factors(fobj);
		exit(1);
	}

	if (((fobj->nfs_obj.startq > 0) || (fobj->nfs_obj.rangeq > 0)) && fobj->nfs_obj.poly_only)
	{
		printf("bad options: -np with -f or -c\n");
		print_factors(fobj);
		exit(1);
	}

	if (((fobj->nfs_obj.startq > 0) || (fobj->nfs_obj.rangeq > 0)) && fobj->nfs_obj.post_only)
	{
		printf("bad options: -nc with -f or -c\n");
		print_factors(fobj);
		exit(1);
	}

	if (fobj->nfs_obj.sieve_only && fobj->nfs_obj.post_only)
	{
		printf("bad options: -nc with -ns\n");
		print_factors(fobj);
		exit(1);
	}

	if (fobj->nfs_obj.sieve_only && fobj->nfs_obj.poly_only)
	{
		printf("bad options: -np with -ns\n");
		print_factors(fobj);
		exit(1);
	}

	if (fobj->nfs_obj.poly_only && fobj->nfs_obj.post_only)
	{
		printf("bad options: -nc with -np\n");
		print_factors(fobj);
		exit(1);
	}

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	NFS_ABORT = 0;
	signal(SIGINT,nfsexit);

	//start a counter for the whole job
	gettimeofday(&start, NULL);
	
	//nfs state machine:
	input = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	statenum = 0;
	process_done = 0;
	while (!process_done)
	{
		switch (statenum)
		{
		case 0: //"init":
			//initialize some job parameters
			job.current_rels = 0;

			//write the input bigint as a string			
			input = mpz_conv2str(&input, 10, fobj->nfs_obj.gmp_n);

			//create an msieve_obj
			//this will initialize the savefile to the outputfile name provided
			obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
				fobj->nfs_obj.fbfile, g_rand.low, g_rand.hi, (uint32)0, nfs_lower, nfs_upper, cpu, 
				(uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);
			fobj->nfs_obj.mobj = obj;

			// initialize these before checking existing files.  If poly
			// select is resumed they will be changed by check_existing_files.
			job.last_leading_coeff = 0;
			job.poly_time = 0;

			//determine what to do next based on the state of various files
			//this will set job.current_rels if it finds any
			is_continuation = check_existing_files(fobj, &last_specialq, &job);

			//get best parameters for the job - call this after checking the input
			//against the savefile because the input could be changed (i.e., reduced
			// by some factor).
			get_ggnfs_params(fobj, &job);			

			//if we are going to be doing sieving, check for the sievers
			if (!(fobj->nfs_obj.poly_only || fobj->nfs_obj.post_only))
			{
				FILE *test;
				// test for existence of the siever
				test = fopen(job.sievername, "rb");
				if (test == NULL)
				{
					printf("WARNING: could not find %s, reverting to siqs!\n",job.sievername);
					logprint_oc(fobj->flogname, "a", "WARNING: could not find %s, reverting to siqs!\n",
						job.sievername);

					mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
					SIQS(fobj);
					mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
					return;
				}
			}

			if (VFLAG >= 0)
				gmp_printf("nfs: commencing gnfs on c%d: %Zd\n",
					gmp_base10(fobj->nfs_obj.gmp_n), fobj->nfs_obj.gmp_n);

			logprint_oc(fobj->flogname, "a", "nfs: commencing gnfs on c%d: %s\n",
				gmp_base10(fobj->nfs_obj.gmp_n),
				mpz_conv2str(&gstr1.s, 10, fobj->nfs_obj.gmp_n));

			//convert input to msieve bigint notation and initialize a list of factors
			gmp2mp_t(fobj->nfs_obj.gmp_n,&mpN);
			factor_list_init(&factor_list);
			
			if (fobj->nfs_obj.rangeq > 0)
				qrange = ceil((double)fobj->nfs_obj.rangeq / (double)THREADS);
			else
				qrange = job.qrange;

			if (is_continuation == 0)
			{
				// didn't find a .job file or a .p file.  This is a new job
				statenum = 9;
			}
			else if (is_continuation == 1)
			{
				// found a .job file, so we are resuming a job having already
				// found a polynomial
				statenum = 10;
			}
			else if (is_continuation > 1)
			{
				// didn't find a .job file, but did find and parse a .p file.
				// resume poly selection
				statenum = 11;
			}
			else
			{
				// something went wrong.  bail.
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;				
			}						

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

			relations_needed = do_msieve_filtering(fobj, obj, &job);
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

			// add restart flag if requested
			if (fobj->nfs_obj.la_restart)
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

			//try this hack - store a pointer to the msieve obj so that
			//we can change a flag on abort in order to interrupt the LA.
			obj_ptr = obj;
			nfs_solve_linear_system(obj, fobj->nfs_obj.gmp_n);
			if (obj_ptr->flags & MSIEVE_FLAG_STOP_SIEVING)
				statenum = 7;
			else
				statenum = 5;

			// reset it back to where it was
			if (LATHREADS > 0)
			{
				msieve_obj_free(obj);
				obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
					fobj->nfs_obj.fbfile, g_rand.low, g_rand.hi, (uint32)0, nfs_lower, nfs_upper, cpu, 
					(uint32)L1CACHE, (uint32)L2CACHE, (uint32)THREADS, (uint32)0, (uint32)0, 0.0);
			}
			
			obj_ptr = NULL;
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

			logprint_oc(fobj->flogname, "a", "nfs: commencing msieve sqrt\n");

			//try this hack - store a pointer to the msieve obj so that
			//we can change a flag on abort in order to interrupt the LA.
			obj_ptr = obj;

			nfs_find_factors(obj, fobj->nfs_obj.gmp_n, &factor_list);

			obj_ptr = NULL;
			
			extract_factors(&factor_list,fobj);

			if (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) == 0)
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

		case 9:		// new factorization			

			startq = job.fblim / 2;
			statenum = 1;								

			//load the job structure with the starting Q.
			job.startq = startq;
			job.qrange = qrange;

			//deal with custom commands
			if (fobj->nfs_obj.poly_only)
				statenum = 1;
			else if (fobj->nfs_obj.sieve_only || (fobj->nfs_obj.startq > 0))
			{
				printf("poly selection must be performed before sieving is possible\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}
			else if (fobj->nfs_obj.post_only)
			{
				printf("poly selection must be performed before post processing is possible\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}
			else if (fobj->nfs_obj.la_restart)
			{
				printf("poly selection must be performed before post processing is possible\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}

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

		case 10:	// resume .job
			if (last_specialq == 0)
			{				
				if (fobj->nfs_obj.startq > 0)
				{
					startq = fobj->nfs_obj.startq;
					if (VFLAG >= 0)
						printf("nfs: continuing with sieving at user specified special-q %u\n",
							startq);

					logprint_oc(fobj->flogname, "a", "nfs: continuing with sieving at user "
						"specified special-q %u\n", startq);
				}
				else
				{
					startq = job.fblim / 2;
					if (VFLAG >= 0)
						printf("nfs: continuing with sieving - could not determine "
							"last special q; using default startq\n");

					logprint_oc(fobj->flogname, "a", "nfs: continuing with sieving "
						"- could not determine last special q; using default startq\n");
				}

				// next step is sieving
				statenum = 2;
			}
			else
			{				
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

				if (VFLAG >= 0)
					printf("nfs: found %u relations, continuing job at specialq = %u\n",
					job.current_rels, startq);

				logprint_oc(fobj->flogname, "a", "nfs: found %u relations, "
					"continuing job at specialq = %u\n",
					job.current_rels, startq);

			}

			//load the job structure with the starting Q.
			job.startq = startq;
			job.qrange = qrange;

			//override with custom commands, if applicable
			if (fobj->nfs_obj.poly_only)
			{
				printf("WARNING: .job file exists!  If you really want to redo poly selection,"
					" delete the .job file.\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}
			else if (fobj->nfs_obj.sieve_only)
				statenum = 2;
			else if (fobj->nfs_obj.post_only == 1)
				statenum = 3;
			else if (fobj->nfs_obj.post_only == 2)
				statenum = 4;
			else if (fobj->nfs_obj.post_only == 3)
				statenum = 5;
			else if (fobj->nfs_obj.la_restart)
				statenum = 4;

			break;

		case 11:		// resume poly select			
			if (VFLAG > 1) printf("nfs: resuming poly select\n");
			fobj->nfs_obj.polystart = job.last_leading_coeff;

			//load the job structure with the starting Q.
			startq = job.fblim / 2;
			job.startq = startq;
			job.qrange = qrange;

			//override with custom commands, if applicable
			if (fobj->nfs_obj.poly_only)
				statenum = 1;
			else if (fobj->nfs_obj.sieve_only)			
			{
				printf("poly selection must be performed before sieving is possible\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}
			else if (fobj->nfs_obj.startq > 0)
			{
				startq = fobj->nfs_obj.startq;
				printf("sieving will start at special-q %u when poly selection finishes\n",
					startq);
			}
			else if (fobj->nfs_obj.post_only)
			{
				printf("poly selection must be performed before post processing is possible\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}
			else if (fobj->nfs_obj.la_restart)
			{
				printf("poly selection must be performed before post processing is possible\n");
				// non ideal solution to infinite loop in factor() if we return without factors
				// (should return error code instead)
				NFS_ABORT = 1;
				process_done = 1;
			}

			statenum = 1;

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


//digits, r/alim, lpbr/a, mfbr/a, r/alambda, siever, rels
//entries based on statistics gathered from many factorizations done
//over the years by myself and others, and from here:
//http://www.mersenneforum.org/showthread.php?t=12365
#define GGNFS_TABLE_ROWS 15
static double ggnfs_table[GGNFS_TABLE_ROWS][8] = {
	{85,  900000,   24, 48, 2.1, 11, 0, 10000},
	{90,  1200000,  25, 50, 2.3, 11, 0, 10000},
	{95,  1500000,  25, 50, 2.5, 12, 0, 10000},
	{100, 1800000,  26, 52, 2.5, 12, 0, 20000},
	{105, 2500000,  26, 52, 2.5, 12, 0, 20000},
	{110, 3200000,  26, 52, 2.5, 13, 0, 20000},
	{115, 4500000,  27, 54, 2.5, 13, 0, 20000},
	{120, 5000000,  27, 54, 2.5, 13, 0, 20000},
	{125, 5500000,  27, 54, 2.5, 13, 0, 40000},
	{130, 6000000,  27, 54, 2.5, 13, 0, 40000},
	{135, 8000000,  27, 54, 2.5, 14, 0, 40000},
	{140, 12000000, 28, 56, 2.5, 14, 0, 40000},
	{145, 15000000, 28, 56, 2.5, 14, 0, 40000},
	{150, 20000000, 29, 58, 2.5, 14, 0, 80000},
	{155, 30000000, 29, 58, 2.5, 15, 0, 80000}
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
	int i, d = gmp_base10(fobj->nfs_obj.gmp_n);
	double scale;
	double fudge;
	int found = 0;

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
			job->fblim = ggnfs_table[0][1];
			job->lpb = ggnfs_table[0][2];
			job->mfb = ggnfs_table[0][3];
			job->lambda = ggnfs_table[0][4];
			job->siever = ggnfs_table[0][5];
			job->qrange = ggnfs_table[0][7];
		}
		else
		{
			job->fblim = ggnfs_table[GGNFS_TABLE_ROWS-1][1];
			job->lpb = ggnfs_table[GGNFS_TABLE_ROWS-1][2];
			job->mfb = ggnfs_table[GGNFS_TABLE_ROWS-1][3];
			job->lambda = ggnfs_table[GGNFS_TABLE_ROWS-1][4];
			job->siever = ggnfs_table[GGNFS_TABLE_ROWS-1][5];
			job->qrange = ggnfs_table[GGNFS_TABLE_ROWS-1][7];
		}
	}

	// appoximate min_rels.  factMsieve.pl resource, except we always
	// use LPBR == LPBA (for now... need to change to allow independent)
	// http://ggnfs.svn.sourceforge.net/viewvc/ggnfs/trunk/tests/factMsieve.pl?r1=374&r2=416
	// http://www.mersenneforum.org/showpost.php?p=294055&postcount=25
	switch (job->lpb)
	{
	case 24:
		fudge = 0.2;
		break;
	case 25:
		fudge = 0.35;
		break;
	case 26:
		fudge = 0.54;
		break;
	case 27:
		fudge = 0.63;
		break;
	case 28:
		fudge = 0.69;
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
	default:
		fudge = 0.2;
		break;
	}

	job->min_rels = (uint32)(2 * fudge * (
		pow(2.0,(double)job->lpb) / log(pow(2.0,(double)job->lpb))));

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
	
