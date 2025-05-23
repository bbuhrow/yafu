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

#include <stdio.h>
#include <gmp.h>
#include "nfs_impl.h"
#include "threadpool.h"
#include <math.h>



#ifdef USE_NFS

// until we can integrate the codebases its easier to just 
// copy this from poly_params.c
const double params_deg4[11][5] = {

	{ 80, 3.00E+013, 2.00E+013, 1.00E-007,  1 * 60},
	{ 85, 3.00E+014, 4.00E+013, 6.50E-008,  1 * 60},
	{ 90, 2.50E+015, 5.00E+014, 3.80E-008,  3 * 60},
	{ 92, 9.96E+015, 4.87E+013, 2.80E-008,  4 * 60},
	{ 95, 2.92E+016, 1.75E+014, 1.85E-008,  6 * 60},
	{100, 1.52E+017, 1.13E+015, 8.75E-009, 12 * 60},
	{105, 8.97E+017, 7.44E+015, 4.40E-009, 23 * 60},
	{110, 6.60E+018, 3.00E+016, 1.40E-009, 45 * 60},
	{115, 1.00E+019, 1.00E+017, 4.09E-010, 90 * 60},
	{120, 3.00E+020, 5.00E+017, 2.14E-010, 180 * 60},
	{125, 1.00E+022, 1.00E+099, 1.12E-010, 8 * 3600},
};

const double params_deg5[24][5] = {

	{100, 9.89E+015, 4.00E+012, 2.95E-009,   12 * 60},
	{105, 4.60E+016, 1.56E+013, 1.60E-009,   23 * 60},
	{110, 2.63E+017, 9.80E+013, 8.76E-010,   45 * 60},
	{115, 1.50E+018, 5.37E+014, 4.87E-010,   90 * 60},
	{120, 8.75E+018, 2.96E+015, 2.61E-010,  180 * 60},
	{125, 4.99E+019, 1.57E+016, 1.32E-010,  360 * 60},
	{130, 2.62E+020, 8.24E+016, 6.57E-011,  720 * 60},
	{135, 1.04E+021, 3.90E+017, 3.59E-011, 1320 * 60},
	{140, 4.43E+021, 2.02E+018, 1.75E-011, 2520 * 60},
	{145, 1.77E+022, 1.03E+019, 8.70E-012, 4900 * 60},
	{150, 7.09E+022, 5.25E+019, 4.35E-012, 9000 * 60},
	{159, 2.00E+024, 2.00E+022, 1.00E-012, 300 * 3600},
	{165, 8.00E+024, 2.00E+023, 5.00E-013, 300 * 3600},
	{170, 5.00E+025, 1.58E+024, 1.50E-013, 300 * 3600},
	{175, 3.00E+026, 1.00E+025, 1.00E-013, 300 * 3600},
	{180, 1.80E+027, 5.36E+025, 7.00E-014, 300 * 3600},
	{185, 1.00E+028, 3.12E+026, 2.00E-014, 300 * 3600},
	{190, 6.00E+028, 1.82E+027, 4.00E-015, 300 * 3600},
	{197, 1.00E+030, 1.00E+029, 2.00E-015, 300 * 3600},
	{200, 3.10E+030, 1.10E+029, 1.50E-015, 300 * 3600},
	{205, 2.00E+031, 5.70E+029, 5.50E-016, 300 * 3600},
	{210, 1.00E+032, 3.00E+030, 1.90E-016, 300 * 3600},
	{215, 6.00E+032, 1.50E+031, 7.00E-017, 300 * 3600},
	{220, 2.40E+033, 7.70E+031, 3.00E-017, 300 * 3600},
};


static const poly_deadline_t time_limits[] = {
    //  bits, seconds
        {248, 15},		    // 74 digits 1
        {264, 30},		    // 80 digits 2
        {304, 60},		    // 92 digits 6
        {320, 3 * 60},		// 97 digits 15
        {348, 9 * 60},		// 105 digits 30
        {365, 1 * 3600},	// 110 digits
        {383, 2 * 3600},	// 116 digits
        {399, 4 * 3600},	// 120 digits
        {416, 8 * 3600},	// 126 digits
        {433, 16 * 3600},	// 131 digits
        {449, 32 * 3600},	// 135 digits
        {466, 64 * 3600},	// 140 digits
        {482, 100 * 3600},	// 146 digits
        {498, 200 * 3600},	// 150 digits
        {514, 300 * 3600},	// 155 digits
};

#define NUM_TIME_LIMITS sizeof(time_limits)/sizeof(time_limits[0])

snfs_t * snfs_find_form(fact_obj_t *fobj)
{
	snfs_t* poly;
    struct timeval stopt;	// stop time of this job
    struct timeval startt;	// start time of this job
    double t_time;

    gettimeofday(&startt, NULL);

	// allow larger one to come in, because some polynomials forms are significantly reduced
	// (i.e., x^(2n))
	if (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 2) > 2*MAX_SNFS_BITS)
	{
		if (fobj->VFLAG >= 0)
			printf("nfs: n is too large for snfs, skipping snfs poly select\n");
		return NULL;
	}	

	poly = (snfs_t*)malloc(sizeof(snfs_t));
	if( !poly )
	{
		printf("gen: initial malloc failed\n");
		exit(-1);
	}
	snfs_init(poly);

	if (poly->form_type == SNFS_NONE)
	{
		if (fobj->VFLAG >= 0) printf("nfs: searching for brent special forms...\n");
		// if this is a factor() run, restore the original input number so that we 
		// can detect these forms on the original input
		if (fobj->autofact_obj.autofact_active)
		{
			mpz_set(fobj->nfs_obj.snfs_cofactor, fobj->nfs_obj.gmp_n);
			mpz_set(fobj->nfs_obj.gmp_n, fobj->N);
		}

		find_brent_form(fobj, poly);

		if (fobj->autofact_obj.autofact_active)
		{
			mpz_set(fobj->nfs_obj.gmp_n, fobj->nfs_obj.snfs_cofactor);
		}
	}

	if (poly->form_type == SNFS_NONE)
	{
		if (fobj->VFLAG >= 0) printf("nfs: searching for homogeneous cunningham special forms...\n");
		find_hcunn_form(fobj, poly);	
	}

	if (poly->form_type == SNFS_NONE)
	{
		if (fobj->VFLAG >= 0) printf("nfs: searching for XYYXF special forms...\n");
		find_xyyxf_form(fobj, poly);
	}

    if (poly->form_type == SNFS_NONE)
    {
        if (fobj->VFLAG >= 0) printf("nfs: searching for direct special forms...\n");
        // if this is a factor() run, restore the original input number so that we 
        // can detect these forms
        if (fobj->autofact_obj.autofact_active)
        {
            mpz_set(fobj->nfs_obj.snfs_cofactor, fobj->nfs_obj.gmp_n);
            mpz_set(fobj->nfs_obj.gmp_n, fobj->N);
        }

        find_direct_form(fobj, poly);

        if (fobj->autofact_obj.autofact_active)
        {
            mpz_set(fobj->nfs_obj.gmp_n, fobj->nfs_obj.snfs_cofactor);
        }
    }

	if (poly->form_type == SNFS_NONE)
	{
		//if (fobj->VFLAG >= 0) printf("nfs: searching for lucas special forms...\n");
		//find_lucas_form(fobj, poly);
	}

    gettimeofday(&stopt, NULL);
    t_time = ytools_difftime(&startt, &stopt);
	if (fobj->VFLAG >= 0) printf("nfs: snfs form detection took %lf seconds\n", t_time);

	if (poly->form_type == SNFS_NONE)
	{
		if (fobj->VFLAG >= 0) printf("nfs: couldn't find special form\n");
		
		snfs_clear(poly);
		free(poly);
		return NULL;
	}

	return poly;
}

int snfs_choose_poly(fact_obj_t* fobj, nfs_job_t* job)
{
	snfs_t* poly, * polys = NULL;
	int i, npoly;
	FILE *f;
	nfs_job_t *jobs, *best;
	int retcode = 0;

	if ((poly = snfs_find_form(fobj)) == NULL)
	{
		// form detection failed
		job->snfs = NULL;
		return 0;
	}

	// once form detection is done, do some quick trial division
	mpz_set(fobj->div_obj.gmp_n, poly->n);
	zTrial(fobj);
	mpz_set(poly->n, fobj->div_obj.gmp_n);

	// with the form detected, create a good polynomial
	if (poly->form_type == SNFS_XYYXF)
	{
		polys = gen_xyyxf_poly(fobj, poly, &npoly);
	}
	else if (poly->form_type == SNFS_LUCAS)
	{
		polys = gen_lucas_poly(fobj, poly, &npoly);
	}
	else
	{
		polys = gen_brent_poly(fobj, poly, &npoly);
	}

	if (npoly == 0)
	{
		job->snfs = NULL;
		snfs_clear(poly);
		free(poly);

		if (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) == 0)
		{
			if (fobj->VFLAG >= 0)
			{
				printf("nfs: finishing with prp primitive poly\n");
			}
			return 2;
		}
		else
		{
			printf("nfs: no snfs polynomial with small coefficients found\n");
		}
		
		// better by gnfs - unable to construct a snfs poly.
		return 0;
	}
    else if (npoly == 1)
    {
		if (fobj->VFLAG >= 0)
		{
			printf("nfs: found %d polynomial\n", npoly);
		}
    }
    else if (npoly > 1)
    {
		if (fobj->VFLAG >= 0)
		{
			printf("nfs: found %d polynomials, selecting best\n", npoly);
		}
    }

	// we've now measured the difficulty for poly's of all common degrees possibly formed
	// in several different ways.  now we have a decision to make based largely on difficulty and 
	// degree.  We want to pick low difficulty, but only if the degree allows the norms on 
	// both sides to be approximately equal.  Sometimes multiple degrees satisfy this requirement
	// approximately equally in which case only test-sieving can really resolve the difference.
	snfs_scale_difficulty(polys, npoly, fobj->VFLAG);
    if (npoly > 1)
    {
        npoly = snfs_rank_polys(fobj, polys, npoly);
    }

	if (fobj->VFLAG > 0 && npoly > 1)
	{		
		int np = MIN(npoly, NUM_SNFS_POLYS);

        if (fobj->LOGFLAG)
        {
            f = fopen(fobj->flogname, "a");
        }
        else
        {
            f = NULL;
        }
		
		printf("\ngen: ========================================================\n"
			"gen: best %d polynomials:\n"
			"gen: ========================================================\n\n", np);
		
		if (f != NULL) logprint(f, "gen: best %d polynomials:\n", np);

		for (i=0; i<np; i++)
		{
			print_snfs(&polys[i], stdout);
			if (f != NULL) print_snfs(&polys[i], f);
		}

		if (f != NULL) fclose(f);
	}

	// initialize job structures to pass into test sieving.  test sieving may modify parameters
	// of the job to optimize rels/q.
	jobs = (nfs_job_t*)malloc(npoly * sizeof(nfs_job_t));
	if( !jobs )
	{
		printf("out of memory\n");
		exit(-1);
	}
	memset(jobs, 0, npoly * sizeof(nfs_job_t));
	
	// assign initial parameters
	for (i=0; i<npoly; i++)
	{
		jobs[i].poly = polys[i].poly;
		jobs[i].snfs = &polys[i];
		get_ggnfs_params(fobj, &jobs[i]);
		skew_snfs_params(fobj, &jobs[i]);
	}

    //printf("gnfs size = %d, size n + 3 = %d, snfs = %d\n", est_gnfs_size(&jobs[0]),
    //    (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) + 3), fobj->nfs_obj.snfs);

	if ((est_gnfs_size(&jobs[0]) > (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) + 3)))
	{
        if (fobj->nfs_obj.snfs)
        {
            // the input is probably faster by gnfs, and the user hasn't specifically chosen
            // snfs, so don't trial sieve and return the gnfs preference
            printf("gen: WARNING: input probably faster by gnfs but proceeding with snfs-specified job\n");
        }
        else
        {
            // faster by gnfs despite found snfs form, skip to gnfs processing.
            retcode = 1;
            copy_job(&jobs[0], job);
            goto cleanup;
        }
	}

	// return best job, which contains the best poly
	best = snfs_test_sieve(fobj, polys, MIN(NUM_SNFS_POLYS,npoly), jobs, 0);

	if (fobj->VFLAG > 0)
	{
        if (fobj->LOGFLAG)
        {
            f = fopen(fobj->flogname, "a");
        }
        else
        {
            f = NULL;
        }

		printf("\ngen: ========================================================\n");
		printf("gen: selected polynomial:\n");
		printf("gen: ========================================================\n\n");

		if (f != NULL) logprint(f, "gen: selected polynomial:\n");

		print_snfs(best->snfs, stdout);
		if (f != NULL) print_snfs(best->snfs, f);
		if (f != NULL) fclose(f);
	}

	// if the norms are close, also check both sides of the selected poly
	if (fabs(log10(best->snfs->rnorm) - log10(best->snfs->anorm)) < 5.0)
	{
		enum special_q_e side = best->poly->side;
		double score = best->test_score;

		best->poly->side = best->poly->side == RATIONAL_SPQ ? ALGEBRAIC_SPQ : RATIONAL_SPQ;

		printf("gen: selected polynomial has close norms (rat=%1.2e, alg=%1.2e)\ngen: testing opposite side\n",
			best->snfs->rnorm, best->snfs->anorm);

		best = snfs_test_sieve(fobj, best->snfs, 1, best, 1);

		if (best->test_score >= score)
		{
			// wasn't better, keep original side.
			best->poly->side = side;
		}
	}

    int do_skew_opt = 0;
    if (do_skew_opt)
    {
        // if requested, dither the skew to attempt to find 
        // a higher score.
        double bestmurph = best->poly->murphy;
        double bestskew = best->poly->skew;
        double origskew = best->poly->skew;

        for (i = 0; i < 100; i++)
        {
            best->poly->skew = origskew * (1 + (0.4 * (rand() / RAND_MAX) - 0.2));
            //printf("on iteration %d trying skew %lf: ", i, best->poly->skew);
            analyze_one_poly_xface(best->snfs);
            if (best->poly->murphy > bestmurph)
            {
                bestskew = best->poly->skew;
                bestmurph = best->poly->murphy;
                printf("on iteration %d found better skew: %lf with murphy score %le\n", 
                    i, bestskew, bestmurph);
            }
        }

        best->poly->skew = bestskew;
    }

	// copy the best job to the object that will be returned from this function
	copy_job(best, job);

	// re-compute the min_rels, in case LPB's changed
	nfs_set_min_rels(job);

	// make the file and fill it
	snfs_make_job_file(fobj, job);
	fill_job_file(fobj, job, PARAM_FLAG_ALL);

	// also create a .fb file
	ggnfs_to_msieve(fobj, job);

	// set the starting q: no longer needed, set by get_ggnfs_params
	// if (fobj->nfs_obj.startq > 0)
	// 	job->startq = fobj->nfs_obj.startq;
	// else
	// 	job->startq = job->snfs->poly->side == RATIONAL_SPQ ? job->rlim/2 : job->alim/2;

cleanup:
    for (i = 0; i < npoly; i++)
    {
        snfs_clear(&polys[i]);
    }
	snfs_clear(poly);
	free(poly);
	free(polys);
	free(jobs);

	return retcode;
}

void split_file(int nthreads, char* base_filename, const char* file_extension)
{
	// split a file into N files, with each file getting 1/Nth the lines
	int i, j;
	char line[8192];
	int count = 0;
	int lines_per_file = 0;
	char fname[80];
	char* strptr;

	sprintf(fname, "%s.%s", base_filename, file_extension);
	FILE* fid = fopen(fname, "r");
	if (fid != NULL)
	{
		// count the lines in the file
		while (!feof(fid))
		{
			strptr = fgets(line, 8192, fid);
			if (strptr == NULL)
			{
				break;
			}
			count++;
		}
		fclose(fid);
	}

	printf("nfs: Found %d lines in %s\n", count, fname);
	if (nthreads == 0)
	{
		printf("invalid number of threads, must be > 0\n");
		exit(1);
	}
	lines_per_file = count / nthreads;

	fid = fopen(fname, "r");
	if (fid != NULL)
	{
		for (i = 0; i < nthreads - 1; i++)
		{
			// note, dealing the lines out like cards, one at
			// a time to N files, will balance the load better.
			// but it requires having N files open at once.  Maybe 
			// it's not a big deal but that makes me nervous somehow.
			// so we do it one file at a time with a good-enough
			// load balancing.
			sprintf(fname, "%s.%d.%s", base_filename, i, file_extension);
			FILE* fid_out = fopen(fname, "w");
			if (fid_out != NULL)
			{
				for (j = 0; j < lines_per_file; j++)
				{
					// need the ability to parse an arbitrary length line
					char line[8192];
					strptr = fgets(line, 8192, fid);
					if (strptr == NULL)
					{
						break;
					}
					fputs(line, fid_out);
					if (feof(fid))
					{
						break;
					}
				}
				fclose(fid_out);
			}
			else
			{
				printf("could not open %s to write\n", fname);
			}

		}
		sprintf(fname, "%s.%d.%s", base_filename, i, file_extension);
		FILE* fid_out = fopen(fname, "w");
		if (fid_out != NULL)
		{
			for (j = 0; j < lines_per_file; j++)
			{
				// need the ability to parse an arbitrary length line
				char line[8192];
				strptr = fgets(line, 8192, fid);
				if (strptr == NULL)
				{
					break;
				}
				fputs(line, fid_out);
				if (feof(fid))
				{
					break;
				}
			}
			fclose(fid_out);
		}
		else
		{
			printf("could not open %s to write\n", fname);
		}
	}


	return;
}

void do_msieve_polyselect(fact_obj_t *fobj, msieve_obj *obj, nfs_job_t *job, 
	mp_t *mpN, factor_list_t *factor_list)
{
	FILE *logfile;
	uint32_t flags;
	double oldbest = 0.;
    double bestscore = 0.;
    double quality_mult = 1.;
    char quality[8];
	int have_new_best = 0;
	int poly_time_exceeded = 0;

	//an array of thread data objects
	nfs_threaddata_t *thread_data;

	// thread work-queue controls
	int threads_working = 0;
	int *thread_queue, *threads_waiting;
#if defined(WIN32) || defined(_WIN64)
	HANDLE queue_lock = NULL;
	HANDLE *queue_events = NULL;
#else
	pthread_mutex_t queue_lock;
	pthread_cond_t queue_cond;
#endif

	int special_polyfind = 0;		// if user has specified np1, nps, or npr
	int i,j,is_startup;
	char syscmd[GSTR_MAXSIZE + 4];
	char master_polyfile[GSTR_MAXSIZE + 2];
	char polyfile_extension[8];
	uint64_t start = 0, range = 0;
	uint32_t deadline, estimated_range_time, this_range_time, total_time, num_ranges;
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	double t_time;
    double e0;
    int sysreturn;

	//file into which we will combine all of the thread results
	strcpy(polyfile_extension, "p");
	snprintf(master_polyfile, GSTR_MAXSIZE + 2, "%s.p",fobj->nfs_obj.outputfile);

	if (job->last_leading_coeff == 0)
	{
		//make sure we are starting from scratch
		remove(master_polyfile);
	}

	// figure out how long poly selection should take
	i = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 2);
	// compute digits the same way msieve does.
	double digits = log(mpz_get_d(fobj->nfs_obj.gmp_n)) / log(10.0);

    if (digits < 108.0) 		/* <= 110 digits */
    {
        fobj->nfs_obj.pref_degree = 4;
    }
    else if (digits < 220.0) 		/* 110-220 digits */
    {
        fobj->nfs_obj.pref_degree = 5;
    }
    else				/* 220+ digits */
    {
        fobj->nfs_obj.pref_degree = 6;
    }

	for (j = 0; j < NUM_TIME_LIMITS; j++) 
    {
        if (i < time_limits[j].bits)
        {
            break;
        }
	}
	
    if (j == NUM_TIME_LIMITS) 
    {
		deadline = time_limits[j-1].seconds;
	}
	else 
    {
		const poly_deadline_t *low = &time_limits[j-1];
		const poly_deadline_t *high = &time_limits[j];
		uint32_t dist = high->bits - low->bits;
		deadline = (uint32_t)(
			 ((double)low->seconds * (high->bits - i) +
			  (double)high->seconds * (i - low->bits)) / dist);
	}    

	// initialize the variable tracking the total time spent (over all threads)
	// so far in the polynomial selection phase
	total_time = 0;

	// if we're resuming a job in this phase, this will be non-zero
	if (job->poly_time > 0)
	{
		if (job->poly_time >= deadline)
		{
			// we've already spent as much time as we need.  
			// parse the master .p file and find the best poly	
			find_best_msieve_poly(fobj, job, fobj->nfs_obj.job_infile, 1);
			//also create a .fb file
			ggnfs_to_msieve(fobj, job);

			return;
		}
		else
		{
			if (fobj->VFLAG > 0)
				printf("nfs: resuming poly search: reducing deadline by %u seconds\n",
					job->poly_time);

            if (fobj->LOGFLAG)
            {
                logfile = fopen(fobj->flogname, "a");
                if (logfile == NULL)
                {
                    printf("could not open yafu logfile for appending\n");
                    printf("fopen error: %s\n", strerror(errno));
                }
                else
                {
                    logprint(logfile, "nfs: resuming poly search, reducing deadline by %u seconds\n",
                        job->poly_time);
                    logprint(logfile, "nfs: resuming poly search, starting at last coefficient %u\n",
                        job->last_leading_coeff);
                    fclose(logfile);
                }
            }

			// pick up where we left off.
			total_time = job->poly_time;

			// also need to pro-rate the deadline, since that is what is used
			// to determine when to stop
			deadline -= job->poly_time;
		}
	}		

    // now we always do "fast", i.e., divide deadline by number of threads,
    // and search for an "avg" poly score by default.
	if (1) //fobj->nfs_obj.poly_option == 0)
	{
		// 'fast' search.  scale by number of threads
		deadline /= fobj->THREADS;
	}

	if ((fobj->nfs_obj.timeout < deadline) && (fobj->nfs_obj.timeout > 1.0))
	{
		deadline = fobj->nfs_obj.timeout;
		if (fobj->VFLAG > 0)
		{
			printf("nfs: reducing poly deadline to specified gnfs timeout of %u seconds\n",
				deadline);
		}
	}

    // if a command line option 'min', 'avg', or 'good' is
    // specified, then allow early abort of poly search when
    // a polynomial of corresponding quality is found, according
    // to this heuristic:
	
#if 0
	// from:
    // http://www.mersenneforum.org/showthread.php?t=16994
    { /* SB: tried L[1/3,c] fit; it is no better than this */
        int digits = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);
        int degree = fobj->nfs_obj.pref_degree;

        e0 = (digits >= 121) ? (0.0607 * digits + 2.25) :
            (0.0526 * digits + 3.23);

        if (degree == 4) e0 = 0.0625 * digits + 1.69;
        e0 = exp(-log(10) * e0);
#ifdef HAVE_CUDA
        e0 *= 1.15;
#endif

#else
	{ /* SB: tried L[1/3,c] fit; it is no better than this */
		// from msieve source...
		int digits = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);
		int degree = fobj->nfs_obj.pref_degree;

		e0 = 0.0625 * digits + 1.69;
		if (degree > 4)
			e0 = (digits >= 121) ?
			(0.0635 * digits + 1.608) :
			(0.0526 * digits + 3.23);
		e0 = exp(-log(10) * e0);
	#ifdef HAVE_CUDA
		e0 *= 1.15;
	#endif

#endif
        /* seen exceptional polys with +40% but that's */
        /* rare. The fit is good for 88..232 digits */

		// best msieve poly scores on mersenneforum
		// https://www.mersenneforum.org/node/15489?p=875709#post875709
		// Once again Batalov found some nice trendlines :
		// 
		// Score_deg5 = 10 ^ (-0.06248 * d - 1.75205) R² = 0.99911
		// 
		// Score_deg6 = 10 ^ (-0.05828 * d - 2.54967) R² = 0.99604


    }

    if (fobj->VFLAG > 0)
    {
        printf("nfs: setting deadline of %u seconds\n", deadline);
        printf("nfs: expecting degree %d poly E from %.2le to > %.2le\n",
            fobj->nfs_obj.pref_degree, e0, 1.15 * e0);
    }

    if (fobj->nfs_obj.poly_option < 3)
    {
        // if we want to use the full time, set the 
        // quality mulitplier really high.  If we find
        // one above this value surely it is ok to stop,
        // even if e.g., 'deep' was specified.
        quality_mult = 1.4;
        strcpy(quality, "awesome");

        fobj->nfs_obj.poly_option == 4;
        printf("psearch options deep, fast, and wide are deprecated, using psearch=avg\n");
    }

    if (fobj->nfs_obj.poly_option == 3)
    {
        quality_mult = 1.0;
        strcpy(quality, "min");
    }
    else if (fobj->nfs_obj.poly_option == 5)
    {
        quality_mult = 1.072;
        strcpy(quality, "good");        
    }
    else
    {
        quality_mult = 1.036;
        strcpy(quality, "avg");
    }

	//start a counter for the poly selection
	gettimeofday(&startt, NULL);

	//init each thread data structure with info needed to do poly search on a range
	thread_data = (nfs_threaddata_t *)malloc(fobj->THREADS * sizeof(nfs_threaddata_t));

	// allocate the queue of threads waiting for work
	thread_queue = (int *)malloc(fobj->THREADS * sizeof(int));
	threads_waiting = (int *)malloc(sizeof(int));

	if (fobj->THREADS > 1)
	{
#if defined(WIN32) || defined(_WIN64)
		queue_lock = CreateMutex( 
			NULL,              // default security attributes
			FALSE,             // initially not owned
			NULL);             // unnamed mutex
		queue_events = (HANDLE *)malloc(fobj->THREADS * sizeof(HANDLE));
#else
		pthread_mutex_init(&queue_lock, NULL);
		pthread_cond_init(&queue_cond, NULL);
#endif
	}


	//set flags to do poly selection
	flags = 0;
	if (fobj->VFLAG > 1)
		flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;

	if (fobj->nfs_obj.np1)
	{
		flags = flags | MSIEVE_FLAG_NFS_POLY1;
		strcpy(polyfile_extension, "m");
		snprintf(master_polyfile, GSTR_MAXSIZE + 2, "%s.m", fobj->nfs_obj.outputfile);
		special_polyfind = 1;
	}
	
	if (fobj->nfs_obj.nps)
	{
		flags = flags | MSIEVE_FLAG_NFS_POLYSIZE;
		strcpy(polyfile_extension, "ms");
		snprintf(master_polyfile, GSTR_MAXSIZE + 2, "%s.ms", fobj->nfs_obj.outputfile);
		special_polyfind = 1;
		// split the main .m file into N parts
		if (fobj->VFLAG > 0 ) 
			printf("nfs: splitting file: %s.m into %d piece(s)\n",
				fobj->nfs_obj.outputfile, fobj->THREADS);
		split_file(fobj->THREADS, fobj->nfs_obj.outputfile, "m");
	}

	if (fobj->nfs_obj.npr)
	{
		flags = flags | MSIEVE_FLAG_NFS_POLYROOT;
		strcpy(polyfile_extension, "p");
		snprintf(master_polyfile, GSTR_MAXSIZE + 2, "%s.p", fobj->nfs_obj.outputfile);
		special_polyfind = 1;
		// split the main .ms file into N parts
		split_file(fobj->THREADS, fobj->nfs_obj.outputfile, "ms");
	}

	if ((fobj->nfs_obj.np1 == 0) && (fobj->nfs_obj.nps == 0) && (fobj->nfs_obj.npr == 0))
	{
		// no special stages selected for polyselect
		flags = flags | MSIEVE_FLAG_NFS_POLY1;
		flags = flags | MSIEVE_FLAG_NFS_POLYSIZE;
		flags = flags | MSIEVE_FLAG_NFS_POLYROOT;
		strcpy(polyfile_extension, "p");
		snprintf(master_polyfile, GSTR_MAXSIZE + 2, "%s.p", fobj->nfs_obj.outputfile);
		special_polyfind = 0;
	}

	if (job->last_leading_coeff == 0)
	{
		//make sure we are starting from scratch
		remove(master_polyfile);
	}


	for (i = 0; i < fobj->THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;		
		t->fobj = fobj;
		t->task = TASK_POLY;

#ifdef HAVE_CUDA
		if (special_polyfind == 0)
		{
			if (fobj->THREADS > 1)
			{
				// if this is a normal multithreaded polyfind but we are using a gpu for stage 1,
				// then split the threads into different tasks.
				// one thread runs the gpu and stores results into files by leading coefficient.
				// the other threads run nps and npr on those files as they become available.


			}
		}
#endif

		// create thread data with dummy range for now
		init_poly_threaddata(t, obj, mpN, factor_list, i, flags,
			deadline, (uint64_t)1, (uint64_t)1001);

		//give this thread a unique index
		t->tindex = i;
		t->is_poly_select = 1;

		t->thread_queue = thread_queue;
        t->threads_waiting = threads_waiting;

		if (fobj->THREADS > 1)
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
            logprint(logfile, "nfs: commencing poly selection with %d threads\n", fobj->THREADS);
            logprint(logfile, "nfs: setting deadline of %u seconds\n", deadline);
            logprint(logfile, "nfs: expecting degree %d poly E from %.2le to > %.2le\n",
                fobj->nfs_obj.pref_degree, e0, 1.15 * e0);
            if (fobj->nfs_obj.poly_option > 2)
            {
                logprint(logfile, "nfs: searching for %s quality poly E > %.2le\n",
                    quality, e0 * quality_mult);
            }
            fclose(logfile);
        }
    }

	// determine the start and range values
	get_polysearch_params(fobj, &start, &range);	

	if (fobj->THREADS > 1)
	{
		// Activate the worker threads one at a time. 
		// Initialize the work queue to say all threads are waiting for work
		for (i = 0; i < fobj->THREADS; i++)
		{
			nfs_start_worker_thread(thread_data + i, 2);
			thread_queue[i] = i;
		}
	}
	*threads_waiting = fobj->THREADS;

	if (fobj->THREADS > 1)
	{
#if defined(WIN32) || defined(_WIN64)
		// nothing
#else
		pthread_mutex_lock(&queue_lock);
#endif
	}	

    //printf("============= starting %d thread poly search with option %d  =================\n", 
    //    fobj->THREADS, fobj->nfs_obj.poly_option);

	char tmpoutfile[80];
	double total_est_sec = 1e20;
	estimated_range_time = 0;
	is_startup = 1;	
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

			if (fobj->THREADS > 1)
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
			{
				tid = 0;
			}

            //printf("syncing poly from thread %d, %d threads waiting, %d threads working\n", 
            //    tid, *threads_waiting, threads_working);

			// pointer to this thread's data
			t = thread_data + tid;

			// regardless of the task we finished, reset whether or
			// not we found a new best poly. we will check the
			// status below if applicable for this thread's finished task.
			have_new_best = 0;

			// syncronize based on the task we just finished
			switch (t->task)
			{
			case TASK_POLY:
			{
				// check the total time spent so far
				gettimeofday(&stopt, NULL);
				t_time = ytools_difftime(&startt, &stopt);

				// update the estimated range time			
				this_range_time = ytools_difftime(&t->thread_start_time, &stopt);

				if (!is_startup)
				{
					total_time += this_range_time;
					num_ranges++;
					estimated_range_time = this_range_time;
				}

				if (!is_startup)
				{
					FILE* fid;

					// if the main poly file doesn't exist yet, put the input number
					// at the top.  this is used to help restart jobs in the polyfind phase
					fid = fopen(master_polyfile, "r");
					if (fid == NULL)

					{
						fid = fopen(master_polyfile, "w");
						gmp_fprintf(fid, "n: %Zd\n", fobj->nfs_obj.gmp_n);
						fclose(fid);
					}
					else
						fclose(fid);

					// combine thread's output with main file
#if defined(WIN32)
					{
						int a;

						// test for cat
						sprintf(syscmd, "cat %s.%s >> %s 2> nul",
							t->polyfilename, polyfile_extension, master_polyfile);
						a = system(syscmd);

						if (a)
						{
							char tmp[80];
							sprintf(tmp, "%s.%s", t->polyfilename, polyfile_extension);
							win_file_concat(tmp, master_polyfile);
						}
					}

#else
					sprintf(syscmd, "cat %s.%s >> %s",
						t->polyfilename, polyfile_extension, master_polyfile);
					sysreturn = system(syscmd);
#endif
					// then stick on the current total elasped time
					// this is used to help restart jobs in the polyfind phase
					if (!special_polyfind)
					{
						fid = fopen(master_polyfile, "a");
						fprintf(fid, "time: %u\n", total_time);
						fclose(fid);
					}

					if (fobj->nfs_obj.nps)
					{
						// remove this threads' portion of the .m file
						snprintf(syscmd, GSTR_MAXSIZE + 4, "%s.%s", t->polyfilename, "m");
						remove(syscmd);
					}
				}

				// remove each thread's .p file after it's copied
				snprintf(syscmd, GSTR_MAXSIZE + 4, "%s.%s", t->polyfilename, polyfile_extension);
				remove(syscmd);

				// also remove the temporary log file
				snprintf(syscmd, GSTR_MAXSIZE + 4, "%s.%d", fobj->nfs_obj.logfile, tid);
				remove(syscmd);

				// and the temporary fb file
				snprintf(syscmd, GSTR_MAXSIZE + 4, "%s.%d", fobj->nfs_obj.fbfile, tid);
				remove(syscmd);

				// free data
				free(t->logfilename);
				free(t->polyfilename);
				msieve_obj_free(t->obj);
				t->obj = NULL;

				// if user has specified "good enough" option then check if
				// we've found one and stop if so
				if ((!is_startup) && (!special_polyfind))
				{
					bestscore = find_best_msieve_poly(fobj, job, fobj->nfs_obj.job_infile, 0);
					if ((bestscore > oldbest) && (oldbest > 1e-30))
					{
						printf("=== new best score %1.4e > old best score %1.4e\n",
							bestscore, oldbest);
						oldbest = bestscore;
						have_new_best = 1;
					}
					else if (oldbest < 1e-30)
					{
						// do at least one round of searching before test
						// sieving, if applicable for this input.
						oldbest = bestscore;
					}
				}
			}
			break;

			case TASK_POLY_TEST_SIEVE:
			{
				// the estimated (single-threaded) total time:
				total_est_sec = ((double)t->job.min_rels / t->test_time) / fobj->THREADS;

				logprint_oc(fobj->flogname, "a", "thread %d test sieve estimating %1.2f seconds to find "
					" %u rels with %d threads at avg rate %1.2f/sec\n",
					t->tindex, total_est_sec, t->job.min_rels, fobj->THREADS, t->test_time);

				remove(t->outfilename);
				remove(t->job_infile_name);
				strncpy(t->outfilename, tmpoutfile, 80);
				sprintf(t->job_infile_name, "%s", fobj->nfs_obj.job_infile);
			}
			break;

			}

			// exit the loop without restarting the thread if abort is asserted.
			if (NFS_ABORT)
			{
				break;
			}

			// now decide what to do next.  If we haven't found a good enough
			// poly yet then go in here.
            if (bestscore < (e0 * quality_mult))
            {
                // if we can re-start the thread such that it is likely to finish before the
                // deadline, go ahead and do so.  Also make sure we at least have
                // one poly before quitting.
                if (((uint32_t)t_time + estimated_range_time <= deadline) ||
					(bestscore < 1e-30)) 
                {
					// we have time to run a new range or we haven't found 
					// a poly yet.  First check a couple special conditions before
					// restarting on a new range.
                    if ((fobj->nfs_obj.polyrange > 0) && !is_startup)
                    {
						// unless the user has specified a custom range search, in which case
						// this thread is done
                        printf("nfs: thread %d finished custom range search\n", tid);
						t->task = TASK_DONE;
                    }
					else if ((special_polyfind) && !is_startup)
					{
						// or if we are doing -np1, -nps, or -npr
						printf("nfs: thread %d finished special polyfind\n", tid);
						t->task = TASK_DONE;
					}
                    else
                    {
						// almost ready to restart this thread on a new range.
						// 
						// one last thing to check:
						// if we just found a new best score and the size of this 
						// job is such that substantial poly-select is involved,
						// do a test sieve of the new best polynomial.
						struct timeval teststart;
						struct timeval teststop;

						if (have_new_best && (!poly_time_exceeded) && 
							(mpz_sizeinbase(fobj->nfs_obj.gmp_n, 2) > fobj->nfs_obj.poly_testsieve))
						{
							// generate a job file for the new best scoring poly
							int tmplogflag = fobj->LOGFLAG;
							fobj->LOGFLAG = 0;
							sprintf(t->job_infile_name, "poly_test_sieve.%d.job", tid);
							bestscore = find_best_msieve_poly(fobj, job, t->job_infile_name, 1);
							fobj->LOGFLAG = tmplogflag;

							if (fobj->VFLAG >= 0)
							{
								printf("nfs: thread %d running test sieve of new best "
									"poly with score %1.4e\n",
									tid, bestscore);
							}
							logprint_oc(fobj->flogname, "a", 
								"nfs: thread %d running test sieve of new best poly "
								"with score %1.4e\n",
								tid, bestscore);

							// restart this thread running a test sieve to a temp data file
							strncpy(tmpoutfile, t->outfilename, 79);
							sprintf(t->outfilename, "poly_test_sieve.%d.dat", tid);

							t->job.poly = NULL;			// NULL poly will get initialized by copy_job
							copy_job(job, &t->job);
							t->task = TASK_POLY_TEST_SIEVE;
						}
						else
						{
							// prepare for a standard restart of this thread with a 
							// new set of leading coefficients.
							gettimeofday(&teststop, NULL);
							t->task = TASK_POLY;
						}

						// if we've used less than X% of the estimated sieve time, continue
						// with polyselect.
						double curr_poly_time = ytools_difftime(&startt, &teststop);
						double est_percent_poly = (curr_poly_time / total_est_sec) * 100;

						if (est_percent_poly < fobj->nfs_obj.poly_percent_max)
						{
							if (t->task == TASK_POLY)
							{
								// initialize the thread for poly select on
								// a new range of coefficients.
								init_poly_threaddata(t, obj, mpN, factor_list, tid, flags,
									deadline, start, start + range);

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
							}

							// signal the job to start
							if (fobj->THREADS > 1)
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
						else
						{
							poly_time_exceeded = 1;
							deadline = 0;

							if (fobj->VFLAG >= 0)
							{
								printf("nfs: current poly select time of %1.2f sec is %1.2f%% of estimated "
									"sieve time (> %d%%), stopping polyselect\n",
									curr_poly_time, est_percent_poly, 
									fobj->nfs_obj.poly_percent_max);
							}
							logprint_oc(fobj->flogname, "a",
								"nfs: current poly select time of %1.2f sec is %1.2f%% of estimated "
								"sieve time (> %d%%), stopping polyselect\n",
								curr_poly_time, est_percent_poly,
								fobj->nfs_obj.poly_percent_max);

							// send a stop signal to all active threads doing standard polyselect.
							for (i = 0; i < fobj->THREADS; i++)
							{
								if ((thread_data[i].obj != NULL) && (!special_polyfind))
								{
									thread_data[i].obj->flags |= MSIEVE_FLAG_STOP_SIEVING;
								}
							}
						}
                    }
                }
				else
				{
					if ((fobj->VFLAG > 0) && (!special_polyfind))
					{
						printf("nfs: range will not finish before deadline, "
							"thread stopping with bestscore = %1.4e\n", bestscore);
					}
					t->task = TASK_DONE;
				}
            }
            else
            {
				// otherwise this thread is done and we send a stop signal 
				// to all active threads doing standard polyselect.
				for (i = 0; i < fobj->THREADS; i++)
				{
					if ((thread_data[i].obj != NULL) && (!special_polyfind))
					{
						thread_data[i].obj->flags |= MSIEVE_FLAG_STOP_SIEVING;
					}
				}
				
                // announce we are finishing and don't restart the thread.
				if ((fobj->VFLAG > 0) && (!special_polyfind))
                {
                    printf("nfs: found poly better than %s quality (e = %1.4e > %1.4e)\n", 
						quality, bestscore, e0 * quality_mult);
                }
            }

			if (fobj->THREADS == 1)
				*threads_waiting = 0;

		}

		// after starting all ranges for the first time, reset this flag
        if (is_startup)
        {
            is_startup = 0;
        }

		// if all threads are done, break out
		if (threads_working == 0)
			break;

		if (fobj->THREADS > 1)
		{
			// wait for a thread to finish and put itself in the waiting queue
#if defined(WIN32) || defined(_WIN64)
			j = WaitForMultipleObjects(
                fobj->THREADS,
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
    t_time = ytools_difftime(&startt, &stopt);

	if (fobj->nfs_obj.polyrange > 0)
	{		
		printf("custom range search complete in %6.4f seconds\n",t_time);
	}
    else if (fobj->VFLAG >= 0)
    {
        printf("elapsed time: %6.4f seconds (%u second deadline); poly select done\n",
            t_time, deadline);
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
            logprint(logfile, "nfs: completed %u ranges of size %" PRIu64 " in %6.4f seconds\n",
                num_ranges, range, t_time);
            fclose(logfile);
        }
    }

	//stop worker threads
	for (i = 0; i < fobj->THREADS; i++)
	{
		//static_conf->tot_poly += thread_data[i].dconf->tot_poly;
		if (fobj->THREADS > 1)
			nfs_stop_worker_thread(thread_data + i, 2);
	}

	//free the thread structure
	free(thread_data);		
	free(thread_queue);
	free(threads_waiting);

#if defined(WIN32) || defined(_WIN64)
	if (fobj->THREADS > 1)
		free(queue_events);
#endif

	if (!NFS_ABORT)
	{
		// parse the master .p file and find the best poly	
		find_best_msieve_poly(fobj, job, fobj->nfs_obj.job_infile, 1);
		// also create a .fb file
		ggnfs_to_msieve(fobj, job);
	}

	return;
}

void get_default_poly4_norms(double digits, double *norm1, double *norm2, double *min_e) {

    uint32_t i;
    uint32_t num_default_entries = 11;
    const double* low, * high;
    double j, k, dist;
    double max_digits;

    /* if the input is too small (large), give up with an answer that will
     be increasingly wrong the further from the ends of the table it gets */

    if (digits < params_deg4[0][0])
    {
        *norm1 = params_deg4[0][1];
        *norm2 = params_deg4[0][2];
        *min_e = params_deg4[0][3];
        return;
    }

    max_digits = params_deg4[num_default_entries - 1][0];
    if (digits >= max_digits)
    {
        *norm1 = params_deg4[num_default_entries - 1][1];
        *norm2 = params_deg4[num_default_entries - 1][2];
        *min_e = params_deg4[num_default_entries - 1][3];
        return;
    }

    /* Otherwise the parameters to use are a weighted average
       of the two table entries the input falls between */

    for (i = 0; i < num_default_entries - 1; i++) {
        if (digits < params_deg4[i + 1][0])
            break;
    }

    low = params_deg4[i];
    high = params_deg4[i + 1];
    dist = high[0] - low[0];
    j = digits - low[0];
    k = high[0] - digits;

    /* use exponential interpolation */
    *norm1 = exp((log(low[1]) * k +
        log(high[1]) * j) / dist);
    *norm2 = exp((log(low[2]) * k +
        log(high[2]) * j) / dist);
    *min_e = exp((log(low[3]) * k +
        log(high[3]) * j) / dist);
}

void get_default_poly5_norms(double digits, double* norm1, double* norm2, double* min_e) {

	uint32_t i;
	uint32_t num_default_entries = 24;
	const double* low, * high;
	double j, k, dist;
	double max_digits;

	/* if the input is too small (large), give up with an answer that will
	 be increasingly wrong the further from the ends of the table it gets */

	if (digits < params_deg5[0][0])
	{
		*norm1 = params_deg5[0][1];
		*norm2 = params_deg5[0][2];
		*min_e = params_deg5[0][3];
		return;
	}

	max_digits = params_deg5[num_default_entries - 1][0];
	if (digits >= max_digits)
	{
		*norm1 = params_deg5[num_default_entries - 1][1];
		*norm2 = params_deg5[num_default_entries - 1][2];
		*min_e = params_deg5[num_default_entries - 1][3];
		return;
	}

	/* Otherwise the parameters to use are a weighted average
	   of the two table entries the input falls between */

	for (i = 0; i < num_default_entries - 1; i++) {
		if (digits < params_deg5[i + 1][0])
			break;
	}

	low = params_deg5[i];
	high = params_deg5[i + 1];
	dist = high[0] - low[0];
	j = digits - low[0];
	k = high[0] - digits;

	//printf("poly-parmas lo entry: %1.4f\n", low[0]);
	//printf("poly-parmas hi entry: %1.4f\n", high[0]);
	//printf("dist, j-weight, k-weight: %1.4f, %1.4f, %1.4f\n", dist, j, k);

	/* use exponential interpolation */
	*norm1 = exp((log(low[1]) * k +
		log(high[1]) * j) / dist);
	*norm2 = exp((log(low[2]) * k +
		log(high[2]) * j) / dist);
	*min_e = exp((log(low[3]) * k +
		log(high[3]) * j) / dist);

	//printf("norm1 lo,hi,interp = %1.4e, %1.4e, %1.4e\n", low[1], high[1], *norm1);
	//printf("norm2 lo,hi,interp = %1.4e, %1.4e, %1.4e\n", low[2], high[2], *norm2);
	//printf("min_e lo,hi,interp = %1.4e, %1.4e, %1.4e\n", low[3], high[3], *min_e);

}

void init_poly_threaddata(nfs_threaddata_t *t, msieve_obj *obj, 
	mp_t *mpN, factor_list_t *factor_list, int tid, uint32_t flags,
	uint32_t deadline, uint64_t start, uint64_t stop)
{
	fact_obj_t *fobj = t->fobj;
	char *nfs_args = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	// compute digits the same way msieve does.
    double digits = log(mpz_get_d(t->fobj->nfs_obj.gmp_n)) / log(10.0); //gmp_base10(t->fobj->nfs_obj.gmp_n);
    int deadline_per_coeff;
	int degree;

    // this is the old deadline table used in msieve prior to version 1023.
    // if we just use the deadline per thread then we spend all our time
    // on a single range of coefficients per thread.  
    // It's a question of searching really deep in a small range of coefficients
    // or scanning lightly through a wider range of coefficients.  The latter
    // is what yafu used to do prior to 1023 so that's what this emulates.
    // Testing with the new approach seems to show that we often find a perfectly
    // acceptable poly fairly quickly, then spend a long time finishing
    // the search for at best an incremental improvement in score.  With
    // larger inputs that might be ok, but for most use below say c130 it
    // seems wasteful.
    //if (digits <= 100.0)
    //    deadline_per_coeff = 5;
    //else if (digits <= 105.0)
    //    deadline_per_coeff = 20;
    //else if (digits <= 110.0)
    //    deadline_per_coeff = 30;
    //else if (digits <= 120.0)
    //    deadline_per_coeff = 50;
    //else if (digits <= 130.0)
    //    deadline_per_coeff = 100;
    //else if (digits <= 140.0)
    //    deadline_per_coeff = 200;
    //else if (digits <= 150.0)
    //    deadline_per_coeff = 400;
    //else if (digits <= 175.0)
    //    deadline_per_coeff = 800;
    //else if (digits <= 200.0)
    //    deadline_per_coeff = 1600;
    //else
    //    deadline_per_coeff = 3200;

	deadline_per_coeff = deadline;

	t->logfilename = (char *)malloc(80 * sizeof(char));
	t->polyfilename = (char *)malloc(80 * sizeof(char));
	t->fbfilename = (char *)malloc(80 * sizeof(char));

    
	double norm1, norm2, min_e;
	if (digits < 108.0)
    {
        get_default_poly4_norms(digits, &norm1, &norm2, &min_e);
		degree = 4;
    }
    else
    {
		get_default_poly5_norms(digits, &norm1, &norm2, &min_e);
		degree = 5;
    }

	// we want to make sure we actually find some polynomials when running on 
	// really small inputs.  The default msieve values don't seem to allow
	// enough polys to be found... here we tweak them a little bit.
	if (digits < 115.0)
	{
		norm1 *= 0.8;
		min_e *= 0.9;

		sprintf(nfs_args, "min_coeff=%" PRIu64 " max_coeff=%" PRIu64 " poly_deadline=%d "
			"stage1_norm=%1.4e stage2_norm=%1.4e min_evalue=%1.4e",
			start, stop, deadline_per_coeff, norm1, norm2, min_e);
	}
	else
	{
		sprintf(nfs_args, "min_coeff=%" PRIu64 " max_coeff=%" PRIu64 " poly_deadline=%d",
			start, stop, deadline_per_coeff);
	}

	if ((t->fobj->VFLAG > 0) && (tid == 0))
	{
		printf("nfs: flags = %08x\n", flags);
		printf("nfs: stage 1 norm = %0.4le\n", norm1);
		printf("nfs: stage 2 norm = %0.4le\n", norm2);
		printf("nfs: min E score  = %0.4le\n", min_e);
		printf("nfs: degree = %d\n", degree);
	}

	sprintf(t->polyfilename,"%s.%d",fobj->nfs_obj.outputfile,tid);
	sprintf(t->logfilename,"%s.%d",fobj->nfs_obj.logfile,tid);
	sprintf(t->fbfilename,"%s.%d",fobj->nfs_obj.fbfile,tid);

	//make sure there isn't a preexisting fb file
	remove(t->fbfilename);

	//create an msieve_obj.  for poly select, the intermediate output file should be specified in
	//the savefile field
	t->obj = msieve_obj_new(obj->input, flags, t->polyfilename, t->logfilename, t->fbfilename, 
		fobj->seed1, fobj->seed2, (uint32_t)0,
		9, (uint32_t)fobj->L1CACHE, (uint32_t)fobj->L2CACHE,
        (uint32_t)fobj->THREADS, (uint32_t)0, nfs_args);

	//pointers to things that are static during poly select
	t->mpN = mpN;
	t->factor_list = factor_list;
	gettimeofday(&t->thread_start_time, NULL);

	return;
}

void get_polysearch_params(fact_obj_t *fobj, uint64_t*start, uint64_t*range)
{
	//search smallish chunks of the space in parallel until we've hit our deadline
	if (fobj->nfs_obj.polystart > 0)
		*start = fobj->nfs_obj.polystart;
	else if (gmp_base10(fobj->nfs_obj.gmp_n) <= 120)
		*start = 120ULL;		// default leading coefficient
	else
		*start = 2048ULL;		// default leading coefficient
	
	if (fobj->nfs_obj.polyrange > 0)
	{
		// search custom range of leading coefficient.  This effectively ignores the deadline
		// because the deadline is only checked after each range is done.  if we assign the
		// entire range up front, the deadline is never checked before polysearch finishes.
		*range = fobj->nfs_obj.polyrange / fobj->THREADS;

		// sanity check
		if (*range == 0) *range = 1;
	}
	else
	{
		// loop on a default range until the time deadline is reached
		*range = (uint64_t)fobj->nfs_obj.polybatch;
	}

	return;
}

void *polyfind_launcher(void *ptr)
{
	// top level function which performs poly search on a range of leading A coefficient values.
	// has pthread calling conventions, meant to be used in a multi-threaded environment
	nfs_threaddata_t *t = (nfs_threaddata_t *)ptr;
	fact_obj_t* fobj = t->fobj;

	switch (t->task)
	{
	case TASK_POLY:
		// remove any temporary relation files
		remove(t->polyfilename);
		MySleep(100);

		if (t->fobj->VFLAG >= 0)
		{
			uint64_t start, stop; // one instance where the new msieve api is rather a pain
			sscanf(t->obj->nfs_args, "min_coeff=%" PRIu64 " max_coeff=%" PRIu64, &start, &stop);
			printf("nfs: thread %d commencing polynomial search over range: %" PRIu64 " - %" PRIu64"\n",
				t->tindex, start, stop);
			fflush(stdout);
		}

		// start polyfind
		factor_gnfs(t->obj, t->mpN, t->factor_list);
		break;

	case TASK_POLY_TEST_SIEVE:
		
		{
			struct timeval stopt;	// stop time of this job
			struct timeval startt;	// start time of this job
			double rels_per_sec[3];
			double this_time;
			uint32_t num_ranges;

			// to get here, poly select found a new best poly, which means
			// that find_best_msieve_poly() was run, which means that 
			// get_ggnfs_params() was called, which filled the job with
			// default/specified parameters, including startq and min_rels.
			// the test sieve runs 3 intervals of 1000 special-q (scale
			// as a function of size?) and finds the area under the
			// rels/sec graph formed by those data points to estimate
			// the total time needed to sieve with this polynomial.
			t->test_time = 0;
			t->job.qrange = 1000;

			// run with default startq
			gettimeofday(&startt, NULL);
			lasieve_launcher((void*)t);
			gettimeofday(&stopt, NULL);

			this_time = ytools_difftime(&startt, &stopt);
			rels_per_sec[0] = t->job.current_rels / this_time;

			// how many 1k ranges would be needed assuming 
			// production is constant at this rate?
			num_ranges = (t->job.min_rels / t->job.current_rels) + 1;

			// note: don't write the logfile here because we're in a thread.
			if (fobj->VFLAG > 0)
			{
				printf("nfs: thread %d test sieve found %1.2f rels/sec at q = %u\n",
					t->tindex, rels_per_sec[0], t->job.startq);
				printf("nfs: estimated total q-range needed is: %u-%u\n",
					t->job.startq, num_ranges * 1000 + t->job.startq);
			}

			// now test at num_ranges/2 and num_ranges to estimate 
			// the changing production rate.
			t->job.startq += num_ranges / 2 * 1000;
			t->job.current_rels = 0;
			gettimeofday(&startt, NULL);
			lasieve_launcher((void*)t);
			gettimeofday(&stopt, NULL);

			this_time = ytools_difftime(&startt, &stopt);
			rels_per_sec[1] = t->job.current_rels / this_time;
			if (fobj->VFLAG > 0)
			{
				printf("nfs: thread %d test sieve found %1.2f rels/sec at q = %u\n",
					t->tindex, rels_per_sec[1], t->job.startq);
			}

			t->job.startq += num_ranges / 2 * 1000;
			t->job.current_rels = 0;
			gettimeofday(&startt, NULL);
			lasieve_launcher((void*)t);
			gettimeofday(&stopt, NULL);

			this_time = ytools_difftime(&startt, &stopt);
			rels_per_sec[2] = t->job.current_rels / this_time;
			if (fobj->VFLAG > 0)
			{
				printf("nfs: thread %d test sieve found %1.2f rels/sec at q = %u\n",
					t->tindex, rels_per_sec[2], t->job.startq);
			}

			// could find the slopes between these points but there
			// is likely so much noise between 1k ranges that we'd
			// probably be fooling ourselves that the slopes are real.
			// find the average rels/sec
			rels_per_sec[0] = (rels_per_sec[0] + rels_per_sec[1] + rels_per_sec[2]) / 3.0;

			// and the estimated (single-threaded) total time:
			t->test_time = (double)t->job.min_rels / rels_per_sec[0];

			if (fobj->VFLAG > 0)
			{
				printf("nfs: thread %d test sieve estimating %1.2f seconds to find "
					" %u rels at avg rate %1.2f/sec\n",
					t->tindex, t->test_time, t->job.min_rels, rels_per_sec[0]);
			}

			// store the average rate, so the main thread can print the above
			// information to the logfile.
			t->test_time = rels_per_sec[0];
		}

		break;

	}

    //printf("**** finished poly work in thread %d\n", t->tindex);

	return 0;
}

#endif
