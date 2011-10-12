#include "nfs.h"
#include "gmp_xface.h"

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
	//expect the input in fobj->nfs_obj.gmp_n
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
	if (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) < fobj->nfs_obj.min_digits)
	{
		mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
		SIQS(fobj);
		mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
		return;
	}	
		
	//first a small amount of trial division
	//which will add anything found to the global factor list
	//this is only called with the main thread
	if (VFLAG > 0)
		printf("nfs: commencing trial factoring\n");
	mpz_set(fobj->div_obj.gmp_n, fobj->nfs_obj.gmp_n);
	fobj->div_obj.print = 0;
	fobj->div_obj.limit = 10000;
	zTrial(fobj);
	mpz_set(fobj->nfs_obj.gmp_n, fobj->div_obj.gmp_n);

	if (mpz_probab_prime_p(fobj->nfs_obj.gmp_n, NUM_WITNESSES))
	{
		z tmpz;
		zInit(&tmpz);
		gmp2mp(fobj->nfs_obj.gmp_n, &tmpz);
		add_to_factor_list(fobj, &tmpz);
		
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
			printf("could not open yafu logfile for appending\n");
		else
		{
			if (VFLAG >= 0)
				gmp_printf("PRP%d = %Zd\n",mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10),
				fobj->nfs_obj.gmp_n);
			logprint(logfile, "PRP%d = %s\n",
				mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->nfs_obj.gmp_n));
			fclose(logfile);
		}		

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
		zFree(&tmpz);
		return;
	}

	if (mpz_perfect_square_p(fobj->nfs_obj.gmp_n))
	{
		z tmpz;
		zInit(&tmpz);

		mpz_sqrt(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n);
		gmp2mp(fobj->nfs_obj.gmp_n, &tmpz);

		add_to_factor_list(fobj, &tmpz);
		logfile = fopen(fobj->flogname, "a");
		if (logfile != NULL)
			logprint(logfile,"prp%d = %s\n",
				mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->nfs_obj.gmp_n));
		add_to_factor_list(fobj, &tmpz);
		if (logfile != NULL)
			logprint(logfile,"prp%d = %s\n",
				mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->nfs_obj.gmp_n));

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
		zFree(&tmpz);
		fclose(logfile);
		return;
	}

	if (mpz_perfect_power_p(fobj->nfs_obj.gmp_n))
	{
		printf("input is a perfect power\n");
		logfile = fopen(fobj->flogname, "a");
		if (logfile != NULL)
		{
			logprint(logfile,"input is a perfect power\n");
			logprint(logfile,"c%d = %s\n",
				mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10),
				mpz_get_str(gstr1.s, 10, fobj->nfs_obj.gmp_n));
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
			logfile = fopen(fobj->flogname, "a");
			if (logfile == NULL)
				printf("could not open yafu logfile for appending\n");
			else
			{
				logprint(logfile, "WARNING: could not find %s, reverting to siqs!\n",job.sievername);
				fclose(logfile);
			}
			mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
			SIQS(fobj);
			mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
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
					gmp_printf("nfs: commencing gnfs on c%d: %Zd\n",
						mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10), fobj->nfs_obj.gmp_n);
				logprint(logfile, "nfs: commencing gnfs on c%d: %s\n",
					mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10),
					mpz_get_str(gstr1.s, 10, fobj->nfs_obj.gmp_n));
				fclose(logfile);
			}

			//write the input bigint as a string
			sInit(&input_str);
			input = mpz_get_str(input, 10, fobj->nfs_obj.gmp_n);

			//create an msieve_obj
			//this will initialize the savefile to the outputfile name provided
			obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile, 
				fobj->nfs_obj.fbfile, seed1, seed2, max_relations, nfs_lower, nfs_upper, cpu, 
				L1CACHE, L2CACHE, THREADS, mem_mb, which_gpu, 0.0);

			fobj->nfs_obj.mobj = obj;

			//convert input to msieve bigint notation and initialize a list of factors
			gmp2mp_t(fobj->nfs_obj.gmp_n,&mpN);
			factor_list_init(&factor_list);

			//see if we can resume a factorization based on the combination of input number,
			//.job file, .fb file, .dat file, and/or .p file.  else, start new job.
			job.current_rels = 0;
			is_continuation = check_existing_files(fobj, &last_specialq, &job);						

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

void get_ggnfs_params(fact_obj_t *fobj, ggnfs_job_t *job)
{
	// based on the size/difficulty of the input number, determine "good" parameters
	// for the following: factor base limits, factor base large prime bound, trial
	// division fudge factors, trial division cutoffs, ggnfs lattice siever area, 
	// and expected number of relations needed.  linearly interpolate between table
	// entries.  keep last valid entry off the ends of the table.  This will produce
	// increasingly poor choices as one goes farther off the table, but you should be
	// doing things by hand by then anyway.
	int i, d = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);
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
	
