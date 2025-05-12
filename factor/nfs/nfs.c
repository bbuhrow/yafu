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
#include "nfs.h"
#include "nfs_impl.h"
#include "gmp_xface.h"
#include <signal.h>
#include <stdlib.h>
#include "ytools.h"
#include <math.h>

#ifndef _WIN32
#include <unistd.h>
#endif

#ifdef _MSC_VER
#include <direct.h>		// _getcwd
#endif

#ifdef USE_NFS

int NFS_ABORT;
int IGNORE_NFS_ABORT;
msieve_obj *obj_ptr;

void set_ggnfs_tables(fact_obj_t *fobj);

void nfsexit(int sig)
{

// idea and following comment "borrowed" from Oliver Weihe's mfaktc...
#ifdef _WIN32
	/* Windows resets the signal handler to the default action once it is
	invoked so we just register it again. */
  	signal(sig, nfsexit);
#endif

	if (IGNORE_NFS_ABORT)
		return;

	if( NFS_ABORT < 1 )
	{
		printf("\nReceived signal %d... please wait\n",sig);
		printf("If you quit again, YAFU will exit immediately but you will LOSE ALL UNSAVED PROGRESS.\n");
	}
	else
		printf("\tStopping immediately!\n");

	NFS_ABORT++;

	if (obj_ptr != NULL)
	{
		printf("setting flag\n");
		obj_ptr->flags |= MSIEVE_FLAG_STOP_SIEVING;
	}
	
	if (NFS_ABORT > 1) exit(1);

	return;
}

int nfs_check_special_case(fact_obj_t *fobj)
{
    int size = gmp_base10(fobj->nfs_obj.gmp_n);
	// below a certain amount, revert to SIQS

	//gmp_printf("checking special cases for nfs input (C%d): %Zd\n", size, fobj->nfs_obj.gmp_n);

	if ((size < fobj->autofact_obj.qs_snfs_xover) && (fobj->autofact_obj.has_snfs_form > 0))
	{
		if (fobj->VFLAG >= 0)
		{
			printf("nfs: snfs size %u < qs_snfs_xover %1.2f (form %u), using SIQS\n",
				size, fobj->autofact_obj.qs_snfs_xover, fobj->autofact_obj.has_snfs_form);
		}
		mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
		SIQS(fobj);
		mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
		return 1;
	}

	if ((size < fobj->autofact_obj.qs_gnfs_xover) && (fobj->autofact_obj.has_snfs_form < 1))
	{
		if (fobj->VFLAG >= 0)
		{
			printf("nfs: gnfs size %u < qs_gnfs_xover %1.2f, using SIQS\n",
				size, fobj->autofact_obj.qs_gnfs_xover);
		}
		mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
		SIQS(fobj);
		mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
		return 1;
	}

	if (is_mpz_prp(fobj->nfs_obj.gmp_n, fobj->NUM_WITNESSES))
	{
		add_to_factor_list(fobj->factors, fobj->nfs_obj.gmp_n,
            fobj->VFLAG, fobj->NUM_WITNESSES);
		
		if (fobj->VFLAG >= 0)
			gmp_printf("PRP%d = %Zd\n", gmp_base10(fobj->nfs_obj.gmp_n),
				fobj->nfs_obj.gmp_n);
		
        char* s = mpz_get_str(NULL, 10, fobj->nfs_obj.gmp_n);
		logprint_oc(fobj->flogname, "a", "PRP%d = %s\n",
			gmp_base10(fobj->nfs_obj.gmp_n), s);	
        free(s);

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
		return 1;
	}

	if (mpz_perfect_square_p(fobj->nfs_obj.gmp_n))
	{
		mpz_sqrt(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n);

		add_to_factor_list(fobj->factors, fobj->nfs_obj.gmp_n, 
            fobj->VFLAG, fobj->NUM_WITNESSES);
        char* s = mpz_get_str(NULL, 10, fobj->nfs_obj.gmp_n);
		logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
			gmp_base10(fobj->nfs_obj.gmp_n), s);

		add_to_factor_list(fobj->factors, fobj->nfs_obj.gmp_n, 
            fobj->VFLAG, fobj->NUM_WITNESSES);
		logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
			gmp_base10(fobj->nfs_obj.gmp_n), s);
        free(s);

		mpz_set_ui(fobj->nfs_obj.gmp_n, 1);
		return 1;
	}

	if (mpz_perfect_power_p(fobj->nfs_obj.gmp_n))
	{
		FILE *flog;
		uint32_t j;

		if (fobj->VFLAG > 0)
			printf("input is a perfect power\n");
		
		factor_perfect_power(fobj, fobj->nfs_obj.gmp_n);

		logprint_oc(fobj->flogname, "a", "input is a perfect power\n");

		for (j=0; j<fobj->factors->num_factors; j++)
		{
			uint32_t k;
            char* s = mpz_get_str(NULL, 10, fobj->factors->factors[j].factor);
			for (k=0; k<fobj->factors->factors[j].count; k++)
			{
				logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
                    gmp_base10(fobj->factors->factors[j].factor), s);
			}
            free(s);
		}

		return 1;
	}

	return 0;
}


// other stuff to look into:
// poly spinning: https://mersenneforum.org/showthread.php?p=565598
// skew optimization (iteration around theoretical skew to find local max of E-score)
// gnfs poly select should probably always be "fast", and min/avg/good/ should apply after that.
// more options to control test sieving (number of points, q-range for each, number to test)

void checkFp(FILE* fp, char* name) {
	if (fp == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not find %s, bailing\n", name);
		exit(-1);
	}
}

void checkFilePresence(char* name) {
	FILE *fp = fopen(name, "rb");
	checkFp(fp, name);
	fclose(fp);
}

//----------------------- NFS ENTRY POINT ------------------------------------//
void nfs(fact_obj_t *fobj)
{
	// expect the input in fobj->nfs_obj.gmp_n
	char *input;
	msieve_obj *obj = NULL;
	char *nfs_args = NULL; // unused as yet
	enum cpu_type cpu = ytools_get_cpu_type();
	mp_t mpN;
	factor_list_t factor_list;
	uint32_t flags = 0;
	nfs_job_t job;
	uint32_t relations_needed = 1;
	uint32_t last_specialq = 0;
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	struct timeval bstop;	// stop time of sieving batch
	struct timeval bstart;	// start time of sieving batch
	struct timeval ustop;	// utility stop time 
	struct timeval ustart;	// utility start time
	double t_time;
	uint32_t pre_batch_rels = 0;
	char tmpstr[GSTR_MAXSIZE];
	int process_done;
	enum nfs_state_e nfs_state;

	// initialize some job parameters
	memset(&job, 0, sizeof(nfs_job_t));

	obj_ptr = NULL;

	if (nfs_check_special_case(fobj))
		return;

	set_ggnfs_tables(fobj);

	if (fobj->nfs_obj.filearg[0] != '\0')
	{
		if (fobj->VFLAG > 0) printf("test: starting trial sieving\n");
		trial_sieve(fobj);
		return;
	}

	// initialize the flag to watch for interrupts, and set the
	// pointer to the function to call if we see a user interrupt
	NFS_ABORT = 0;
	signal(SIGINT,nfsexit);

	// start a counter for the whole job
	gettimeofday(&start, NULL);

	// nfs state machine:
	input = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	nfs_state = NFS_STATE_INIT;
	process_done = 0;

	// Used to predict CADO work file names
	// https://gitlab.inria.fr/cado-nfs/cado-nfs/-/blob/master/scripts/cadofactor/toplevel.py?ref_type=heads#L169
	int sizeOfN = gmp_base10(fobj->nfs_obj.gmp_n);
	int cadoPower;
	if (sizeOfN < 200) {
		cadoPower = ((sizeOfN + 2) / 5) * 5;
	} else {
		cadoPower = ((sizeOfN + 5) / 10) * 10;
	}

	while (!process_done)
	{
        char* s;

		switch (nfs_state)
		{
		case NFS_STATE_INIT:
			// write the input bigint as a string
			//input = mpz_conv2str(&input, 10, fobj->nfs_obj.gmp_n);

            s = mpz_get_str(NULL, 10, fobj->nfs_obj.gmp_n);
            strcpy(input, s);
            free(s);

			// create an msieve_obj
			// this will initialize the savefile to the outputfile name provided
			obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile,
				fobj->nfs_obj.fbfile, fobj->seed1, fobj->seed2, (uint32_t)0, 9,
				(uint32_t)fobj->L1CACHE, (uint32_t)fobj->L2CACHE, 
                (uint32_t)fobj->THREADS, (uint32_t)0, nfs_args);
			fobj->nfs_obj.mobj = obj;

			// initialize these before checking existing files.  If poly
			// select is resumed they will be changed by check_existing_files.
			job.last_leading_coeff = 0;
			job.poly_time = 0;
			job.use_max_rels = 0;
			job.snfs = NULL;

			if (fobj->nfs_obj.cadoMsieve) {
				// CADO-NFS takes care of resuming progress
				nfs_state = NFS_STATE_POLY;
			} else {
				// determine what to do next based on the state of various files.
				// this will set job.current_rels if it finds any
				nfs_state = check_existing_files(fobj, &last_specialq, &job);
			}

			// before we get started, check to make sure we can find ggnfs sievers
			// if we are going to be doing sieving
            if (check_for_sievers(fobj, 1) == 1)
            {
                nfs_state = NFS_STATE_DONE;
            }

            if ((fobj->VFLAG >= 0) && (nfs_state != NFS_STATE_DONE))
            {
                gmp_printf("nfs: commencing nfs on c%d: %Zd\n",
                    gmp_base10(fobj->nfs_obj.gmp_n), fobj->nfs_obj.gmp_n);
            }

            if (nfs_state != NFS_STATE_DONE)
            {
                char* s = mpz_get_str(NULL, 10, fobj->nfs_obj.gmp_n);
                logprint_oc(fobj->flogname, "a", "nfs: commencing nfs on c%d: %s\n",
                    gmp_base10(fobj->nfs_obj.gmp_n), s);
                free(s);
            }

			// convert input to msieve bigint notation and initialize a list of factors
			gmp2mp_t(fobj->nfs_obj.gmp_n, &mpN);
			factor_list_init(&factor_list);
            factor_list_add(obj, &factor_list, &mpN);

            if (fobj->nfs_obj.rangeq > 0)
            {
				// user specified a starting and stopping Q, so we sieve
				// that Q-range only (rest of job is default).
                job.qrange = ceil((double)fobj->nfs_obj.rangeq / (double)fobj->THREADS);
            }

			break;
		case NFS_STATE_POLY:

			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_POLY))
			{
				mpz_t orig_n;
				int better_by_gnfs = 0;
				int check_gnfs = 0;

				mpz_init(orig_n);
				mpz_set(orig_n, fobj->nfs_obj.gmp_n);

				// always check snfs forms (it is fast)
				better_by_gnfs = snfs_choose_poly(fobj, &job);

				if (mpz_cmp(orig_n, fobj->nfs_obj.gmp_n) != 0)
				{
					// the number changed during snfs poly detection,
					// i.e., it was algebraically factored.  Need
					// to reevaluate the input against GNFS/SIQS xovers.
					check_gnfs = 1;
				}

				mpz_clear(orig_n);

				if (better_by_gnfs == 2)
				{
					// primitive factor detection reduced the input
					// to a (probable) prime number.  We are done.
					nfs_state = NFS_STATE_CLEANUP;
					break;
				}

				if (fobj->VFLAG > 1)
				{
					printf("nfs: snfs poly find flags are\n\tsnfs found = %d\n"
						"\tgnfs specified = %d\n\tsnfs found, better by gnfs = %d\n"
						"\tsnfs specified = %d\n\tsnfs algebraic factor removed, recheck xovers = %d\n",
						(job.snfs == NULL), fobj->nfs_obj.gnfs, better_by_gnfs,
						fobj->nfs_obj.snfs, check_gnfs);
				}

				if ((job.snfs == NULL) || fobj->nfs_obj.gnfs ||
					(better_by_gnfs && !fobj->nfs_obj.snfs) ||
					check_gnfs)
				{ 
					// either we never were doing snfs, or the user selected gnfs,
					// or snfs form detect failed.
					if (fobj->nfs_obj.snfs)
					{
						// if the latter then bail with an error because the user 
						// explicitly wants to run snfs...
						printf("nfs: failed to find snfs polynomial!\n");
						printf("nfs: removing the -snfs option could allow the number to be "
							"completed by GNFS or SIQS\n");
						exit(-1);
					}

					// check if this is really a siqs job
					if ((!(job.snfs == NULL)) &&
						(mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) < fobj->autofact_obj.qs_snfs_xover))
					{
						if (fobj->VFLAG >= 0)
						{
							printf("nfs: input snfs form is better done by siqs: "
								"difficulty = %1.2f, size = %d, actual size = %d\n",
								job.snfs->difficulty, est_gnfs_size(&job),
								(int)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10));
							fflush(stdout);
						}

						logprint_oc(fobj->flogname, "a", "nfs: input snfs form is better done by siqs"
							": difficulty = %1.2f, size = %d, actual size = %d\n",
							job.snfs->difficulty, est_gnfs_size(&job), (int)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10));

						mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
						SIQS(fobj);
						mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);

						nfs_state = NFS_STATE_CLEANUP;
						break;
					}
					else if (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) < fobj->autofact_obj.qs_gnfs_xover)
					{

						if (fobj->VFLAG >= 0)
						{
							printf("nfs: non-snfs input of size %d is better done by siqs\n",
								(int)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10));
							fflush(stdout);
						}

						logprint_oc(fobj->flogname, "a", 
							"nfs: non-snfs input of size %d is better done by siqs\n",
							(int)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10));

						mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
						SIQS(fobj);
						mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);

						nfs_state = NFS_STATE_CLEANUP;
						break;
					}
					else if (better_by_gnfs)
					{
						if (fobj->VFLAG >= 0)
						{
							printf("nfs: input snfs form is better done by gnfs: "
								"difficulty = %1.2f, size = %d, actual size = %d\n",
								job.snfs->difficulty, est_gnfs_size(&job), (int)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10));
						}
						logprint_oc(fobj->flogname, "a", "nfs: input snfs form is better done by gnfs"
							": difficulty = %1.2f, size = %d, actual size = %d\n",
							job.snfs->difficulty, est_gnfs_size(&job), (int)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10));

						// get rid of the snfs poly since we'll be proceeding with gnfs.
						snfs_clear(job.snfs);
						free(job.snfs);
						job.snfs = NULL;
					}

					// prepare for gnfs polyselect
					job.alambda = 0.0;
					job.rlambda = 0.0;
					job.lpba = 0;
					job.lpbr = 0;
					job.mfba = 0;
					job.mfbr = 0;
					job.alim = 0;
					job.rlim = 0;
					fobj->nfs_obj.siever = 0;

					if (fobj->nfs_obj.cadoMsieve) {
						// Let CADO-NFS find the poly
						nfs_state = NFS_STATE_CADO;
						break;
					}

					// init job.poly for gnfs
					job.poly = (mpz_polys_t*)malloc(sizeof(mpz_polys_t));
					if (job.poly == NULL)
					{
						printf("nfs: couldn't allocate memory!\n");
						exit(-1);
					}
					mpz_polys_init(job.poly);
					job.poly->rat.degree = 1; // maybe way off in the future this isn't true
					// assume gnfs for now
					job.poly->side = ALGEBRAIC_SPQ;
                    job.poly->size = (double)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);

					gettimeofday(&ustart, NULL);
					do_msieve_polyselect(fobj, obj, &job, &mpN, &factor_list);
					gettimeofday(&ustop, NULL);
					t_time = ytools_difftime(&ustart, &ustop);
					fobj->nfs_obj.poly_time += t_time;
				}
				else
				{
					fobj->nfs_obj.snfs = 1;
					mpz_set(fobj->nfs_obj.gmp_n, job.snfs->n);
				}
			}

			nfs_state = NFS_STATE_SIEVE;
			break;

		case NFS_STATE_SIEVE:
			if (fobj->nfs_obj.cadoMsieve) {
				nfs_state = NFS_STATE_CADO;
				break;
			}

			pre_batch_rels = job.current_rels;
			gettimeofday(&bstart, NULL);

			// sieve if the user has requested to (or by default).  else,
			// set the done sieving flag.  this will prevent some infinite loops,
			// for instance if we only want to post-process, but filtering 
			// doesn't produce a matrix.  if we don't want to sieve in that case,
			// then we're done.
            if (((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
                (fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE)) &&
                !(fobj->nfs_obj.nfs_phases & NFS_DONE_SIEVING))
            {
                // this is not a threadpool, so there can be imbalances between
                // threads on multi-threaded runs where we end up waiting
                // for the last one to finish.
                // todo: make do_sieving a threadpool that incorporates the 
                // logic below (including timeout!) and the logic of NFS_STATE_FILTCHECK 
                // so that we can keep sieving until it is actually time to
                // filter.
				gettimeofday(&ustart, NULL);
                do_sieving_nfs(fobj, &job);
				gettimeofday(&ustop, NULL);
				t_time = ytools_difftime(&ustart, &ustop);
				fobj->nfs_obj.sieve_time += t_time;
            }
            else
            {
                fobj->nfs_obj.nfs_phases |= NFS_DONE_SIEVING;
            }

			// if this has been previously marked, then go ahead and exit.
            if (fobj->nfs_obj.nfs_phases & NFS_DONE_SIEVING)
            {
                process_done = 1;
            }
			
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

		case NFS_STATE_CADO: {
#if defined(WIN32) || defined(_WIN64)
			printf("cadoMsieve is not available on Windows! Bailing\n");
			exit(-1);
#endif

			FILE *fp;
			char buffer[1024];

			if (job.snfs != NULL) {
				// Select param file based on converted difficulty
				cadoPower = est_gnfs_size(&job);
				cadoPower -= cadoPower % 5;
			}

			FILE *dat = fopen("nfs.dat", "a+");

			// First run
			if (job.current_rels == 0) {
				if (fgetc(dat) == 'N') {
					// https://www.geeksforgeeks.org/c-program-count-number-lines-file/
					for (char c = fgetc(dat); !feof(dat); c = fgetc(dat)) {
						if (c == '\n') {
							job.current_rels++;
						}
					}

					// Don't count the header
					job.current_rels--;

					if (job.current_rels > job.min_rels) {
						// Immediately try to filter again
						job.min_rels = job.current_rels;
					}

					printf("nfs: found %d relations in nfs.dat\n", job.current_rels);
				} else {
					// Write the N header
					sprintf(buffer, "N %s\n", input);
					fwrite(buffer, sizeof(char), strlen(buffer), dat);
					printf("nfs: created nfs.dat\n");
				}
			}

			// Check for cado-nfs.py presence
			sprintf(buffer, "%scado-nfs.py", fobj->nfs_obj.cado_dir);
			checkFilePresence(buffer);

			printf("nfs: calling cado-nfs\n");
			char syscmd[8192];

			// Specify path, param file and N
			sprintf(syscmd, "%s %sparameters/factor/params.c%d N=%s ", buffer, fobj->nfs_obj.cado_dir, cadoPower, input);

			// Roundup of fobj->num_threads / 2
			int nrclients = (fobj->num_threads / 2) + (obj->num_threads % 2);
			// Spawn clients based on fobj->num_threads
			// Use sprintf to append to end of syscmd
			sprintf(syscmd + strlen(syscmd), "slaves.nrclients=%d slaves.hostnames=localhost ", nrclients);

			// Stop after sieving enough relations
			sprintf(syscmd + strlen(syscmd), "tasks.filter.run=false ");

			// Reduce time before retry from 10 seconds to 1 second
			sprintf(syscmd + strlen(syscmd), "slaves.downloadretry=1 ");

			if (job.min_rels != 0) {
				// Specify relations wanted
				sprintf(syscmd + strlen(syscmd), "tasks.sieve.rels_wanted=%d ", job.min_rels);
			}

			// SNFS, hell yeah!
			// Information source: https://www.mersenneforum.org/showthread.php?t=24842
			if (job.snfs != NULL) {
				// Create poly file
				FILE *poly = fopen("nfs.poly", "w");

				sprintf(buffer, "n: %s\n", input);
				sprintf(buffer + strlen(buffer), "skew: %f\n", job.snfs->poly->skew);

				for (int deg = MAX_POLY_DEGREE; deg >= 0; deg--) {
					s = mpz_get_str(NULL, 10, job.snfs->c[deg]);
					sprintf(buffer + strlen(buffer), "c%d: %s\n", deg, s);
					free(s);
				}

				for (int deg = 1; deg >= 0; deg--) {
					s = mpz_get_str(NULL, 10, job.snfs->poly->rat.coeff[deg]);
					sprintf(buffer + strlen(buffer), "Y%d: %s\n", deg, s);
					free(s);
				}

				fwrite(buffer, sizeof(char), strlen(buffer), poly);
				fclose(poly);

				// Make CADO use our poly
				sprintf(syscmd + strlen(syscmd), "tasks.polyselect.import=nfs.poly ");

				int sqside;
				if (job.snfs->poly->side == RATIONAL_SPQ) {
					sqside = 0;
				} else {
					sqside = 1;
				}

				// Specify whether to sieve rational or algebraic side
				sprintf(syscmd + strlen(syscmd), "tasks.sieve.sqside=%d ", sqside);

				// CADO seems to know the right sieving parameters, so we ignore the params from nfs.job
				/*
				// Specify sieving parameters
				// Apparently 0=r, 1=a
				sprintf(syscmd + strlen(syscmd), "tasks.lim0=%d tasks.lim1=%d ", job.rlim, job.alim);
				sprintf(syscmd + strlen(syscmd), "tasks.lpb0=%d tasks.lpb1=%d ", job.lpbr, job.lpba);
				sprintf(syscmd + strlen(syscmd), "tasks.sieve.mfb0=%d tasks.sieve.mfb1=%d ", job.mfbr, job.mfba);
				sprintf(syscmd + strlen(syscmd), "tasks.sieve.lambda0=%f tasks.sieve.lambda1=%f ", job.rlambda, job.alambda);
				*/
			}

			// Specify work directory as an absolute path
			char cwdBuf[4097];
#ifdef _MSC_VER
			_getcwd(cwdBuf, 4097);
#else
			getcwd(cwdBuf, 4097);
#endif
			sprintf(syscmd + strlen(syscmd), "-w %s/cadoWorkdir", cwdBuf);

			printf("nfs: cmdline: %s\n", syscmd);
			system(syscmd);

			// Check for convert_poly presence
			checkFilePresence(fobj->nfs_obj.convert_poly_path);

			// Ensure CADO did find the poly
			sprintf(buffer, "./cadoWorkdir/c%d.poly", cadoPower);
			checkFilePresence(buffer);

			printf("nfs: calling convert_poly to create nfs.fb from c*.poly\n");
			// TODO: Support Windows lol
			sprintf(syscmd, "%s -of msieve < ./cadoWorkdir/c%d.poly > nfs.fb", fobj->nfs_obj.convert_poly_path, cadoPower);
			system(syscmd);

			printf("nfs: appending CADO relations into nfs.dat\n");

			// By default, there are no cross-platform way of listing all files under a directory, so we have to find filenames through the log file
			sprintf(buffer, "cadoWorkdir/c%d.log", cadoPower);
			FILE *logFile = fopen(buffer, "r");
			checkFp(logFile, buffer);

			char logLine[3072];
			while (fgets(logLine, 3072, logFile)) {
				// Check if this logLine has relation filename
				char *match = " relations in '";
				char *filename = strstr(logLine, match);
				if (filename == NULL) continue;

				// Move beyond " relations in '" text
				filename += strlen(match);
				// Read ptr until next ' character
				strtok(filename, "'");

				// filename is now something like "/home/nyancat/Tools/yafu-combined/cadoWorkdir/c85.upload/c85.146453-147000.s_yzjb7c.gz"
				// printf("Extracting %s\n", filename);

				// Extract the relations into nfs.cado
				// Don't keep the file, as that would cause us to read it again the next run
				// Suppress stderr with 2>/dev/null
				sprintf(syscmd, "gunzip -c %s 1>nfs.cado 2>/dev/null && rm %s", filename, filename);
				system(syscmd);
				checkFilePresence("nfs.cado");

				// Read the relations
				FILE *relatFile = fopen("nfs.cado", "r");
				checkFp(relatFile, "nfs.cado");

				char relatLine[1024];
				while (fgets(relatLine, 1024, relatFile)) {
					// Is this line a comment?
					if (relatLine[0] == '#') continue;

					// fgets include \n already
					fwrite(relatLine, sizeof(char), strlen(relatLine), dat);
					job.current_rels++;
				}
				fclose(relatFile);
			}
			fclose(logFile);
			fclose(dat);

			// min_rels is not set by YAFU during GNFS
			if (job.min_rels == 0) {
				job.min_rels = job.current_rels;
			}

			printf("nfs: now have %d relations\n", job.current_rels);
			nfs_state = NFS_STATE_FILTER;
			break;
		}

		case NFS_STATE_FILTER:

			// if we've flagged not to do filtering, then assume we have
			// enough relations and move on to linear algebra
			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_FILTER))
			{
				gettimeofday(&ustart, NULL);
				relations_needed = do_msieve_filtering(fobj, obj, &job);
				gettimeofday(&ustop, NULL);
				t_time = ytools_difftime(&ustart, &ustop);
				fobj->nfs_obj.filter_time += t_time;
			}
			else
				relations_needed = 0;

			if (relations_needed == 0)
				nfs_state = NFS_STATE_LINALG;
			else
			{
				// if we filtered, but didn't produce a matrix, raise the target
				// min rels by a few percent.  this will prevent too frequent
				// filtering attempts while allowing the q_batch size to remain small.
				if (job.current_rels > job.min_rels)
					job.min_rels = job.current_rels * fobj->nfs_obj.filter_min_rels_nudge;
				else
					job.min_rels *= fobj->nfs_obj.filter_min_rels_nudge;

				if (fobj->VFLAG > 0)
					printf("nfs: raising min_rels by %1.2f percent to %u\n", 
					100*(fobj->nfs_obj.filter_min_rels_nudge-1), job.min_rels);

				logprint_oc(fobj->flogname, "a", "nfs: raising min_rels by %1.2f percent to %u\n", 
					100*(fobj->nfs_obj.filter_min_rels_nudge-1), job.min_rels);

				nfs_state = NFS_STATE_SIEVE;
			}

            if (fobj->VFLAG > 0)
            {
                gettimeofday(&stop, NULL);
                t_time = ytools_difftime(&start, &stop);
                printf("Elapsed time is now %6.4f seconds.\n", t_time);
            }

			break;

		case NFS_STATE_LINALG:

			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA_RESUME))
			{
				// msieve: build matrix
				flags = 0;
				flags = flags | MSIEVE_FLAG_USE_LOGFILE;
				if (fobj->VFLAG > 0)
					flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
				flags = flags | MSIEVE_FLAG_NFS_LA;

				// add restart flag if requested
				if (fobj->nfs_obj.nfs_phases & NFS_PHASE_LA_RESUME)
					flags |= MSIEVE_FLAG_NFS_LA_RESTART;

				obj->flags = flags;

				if (fobj->VFLAG >= 0)
					printf("nfs: commencing msieve linear algebra\n");

				logprint_oc(fobj->flogname, "a", "nfs: commencing msieve linear algebra\n");

				// use a different number of threads for the LA, if requested
				if (fobj->LATHREADS > 0)
				{
					msieve_obj_free(obj);
					obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile,
						fobj->nfs_obj.fbfile, fobj->seed1, fobj->seed2, (uint32_t)0, 9,
						(uint32_t)fobj->L1CACHE, (uint32_t)fobj->L2CACHE,
                        (uint32_t)fobj->LATHREADS, (uint32_t)0, nfs_args);
				}

				// try this hack - store a pointer to the msieve obj so that
				// we can change a flag on abort in order to interrupt the LA.
				obj_ptr = obj;

				gettimeofday(&ustart, NULL);
				nfs_solve_linear_system(obj, fobj->nfs_obj.gmp_n);
				gettimeofday(&ustop, NULL);
				t_time = ytools_difftime(&ustart, &ustop);
				fobj->nfs_obj.la_time += t_time;
				

				if (obj_ptr->flags & MSIEVE_FLAG_STOP_SIEVING)
					nfs_state = NFS_STATE_DONE;
				else
				{
					// check for a .dat.deps file.  if we don't have one, assume
					// its because the job was way oversieved and only trivial
					// dependencies were found.  try again from filtering with
					// 20% less relations.
					FILE *t;
					sprintf(tmpstr, "%s.dep", fobj->nfs_obj.outputfile);
					if ((t = fopen(tmpstr, "r")) == NULL)
					{
						if (job.use_max_rels > 0)
						{
							// we've already tried again with an attempted fix to the trivial
							// dependencies problem, so either that wasn't the problem or
							// it didn't work.  either way, give up.
							printf("nfs: no dependency file retry failed\n");
							fobj->flags |= FACTOR_INTERRUPT;
							nfs_state = NFS_STATE_DONE;
						}
						else
						{
							// this should be sufficient to produce a matrix, but not too much
							// to trigger the assumed failure mode.							
							if (job.min_rels == 0)
							{
								// if min_rels is not set, then we need to parse the .job file to compute it.
								parse_job_file(fobj, &job);
								nfs_set_min_rels(&job);
							}
							job.use_max_rels = job.min_rels * 1.5;
							printf("nfs: no dependency file found - trying again with %u relations\n",
								job.use_max_rels);
							nfs_state = NFS_STATE_FILTER;
						}
					}
					else
					{
						fclose(t);
						nfs_state = NFS_STATE_SQRT;
					}
				}

				// set the msieve threads back to where it was if we used
				// a different amount for linalg
				if (fobj->LATHREADS > 0)
				{
					msieve_obj_free(obj);
					obj = msieve_obj_new(input, flags, fobj->nfs_obj.outputfile, fobj->nfs_obj.logfile,
						fobj->nfs_obj.fbfile, fobj->seed1, fobj->seed2, (uint32_t)0, 9,
						(uint32_t)fobj->L1CACHE, (uint32_t)fobj->L2CACHE, 
                        (uint32_t)fobj->THREADS, (uint32_t)0, nfs_args);
				}

				obj_ptr = NULL;
			}
			else // not doing linalg
				nfs_state = NFS_STATE_SQRT;

            if (fobj->VFLAG > 0)
            {
                gettimeofday(&stop, NULL);
                t_time = ytools_difftime(&start, &stop);
                printf("Elapsed time is now %6.4f seconds.\n", t_time);
            }
			break;

		case NFS_STATE_SQRT:

            // TODO:
            // when the input has more than 2 factors we can continue 
            // with dependencies until the input is completely factored.
            // as it is now, we stop after the first one, potentially leaving
            // a significant number behind (and deleting all data so that
            // it can't be found manually either).

			if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_SQRT))
			{
				uint32_t retcode;

				// msieve: find factors
				flags = 0;
				flags = flags | MSIEVE_FLAG_USE_LOGFILE;
				if (fobj->VFLAG > 0)
					flags = flags | MSIEVE_FLAG_LOG_TO_STDOUT;
				flags = flags | MSIEVE_FLAG_NFS_SQRT;
				obj->flags = flags;

				if (fobj->VFLAG >= 0)
					printf("nfs: commencing msieve sqrt\n");

				logprint_oc(fobj->flogname, "a", "nfs: commencing msieve sqrt\n");

				// try this hack - store a pointer to the msieve obj so that
				// we can change a flag on abort in order to interrupt the sqrt.
				gmp_sprintf(obj->input, "%Zd", fobj->nfs_obj.gmp_n);
				obj_ptr = obj;

				gettimeofday(&ustart, NULL);
				retcode = nfs_find_factors(obj, fobj->nfs_obj.gmp_n, &factor_list);
				gettimeofday(&ustop, NULL);
				t_time = ytools_difftime(&ustart, &ustop);
				fobj->nfs_obj.sqrt_time += t_time;

				obj_ptr = NULL;

				if (retcode)
				{
					extract_factors(&factor_list,fobj);

					if (mpz_cmp_ui(fobj->nfs_obj.gmp_n, 1) == 0)
						nfs_state = NFS_STATE_CLEANUP;		//completely factored, clean up everything
					else
						nfs_state = NFS_STATE_DONE;		//not factored completely, keep files and stop
				}
				else
				{
					if (fobj->VFLAG >= 0)
					{
						printf("nfs: failed to find factors... possibly no dependencies found\n");
						printf("nfs: continuing with sieving\n");
					}

					logprint_oc(fobj->flogname, "a", "nfs: failed to find factors... "
						"possibly no dependencies found\n"
						"nfs: continuing with sieving\n");

					nfs_state = NFS_STATE_SIEVE;
				}				
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
			sprintf(tmpstr, "%s.ranges", fobj->nfs_obj.outputfile);	remove(tmpstr);
			sprintf(tmpstr, "%s.mat.chk",fobj->nfs_obj.outputfile);	remove(tmpstr);

			if (fobj->nfs_obj.cadoMsieve) {
				remove("nfs.poly");
				remove("nfs.cado");
				sprintf(tmpstr, "c%d.db", cadoPower);    	remove(tmpstr);
				sprintf(tmpstr, "c%d.db-shm", cadoPower);	remove(tmpstr);
				sprintf(tmpstr, "c%d.db-wal", cadoPower);	remove(tmpstr);
				system("rm -r cadoWorkdir");
			}

			gettimeofday(&stop, NULL);

            t_time = ytools_difftime(&start, &stop);

			if (fobj->VFLAG >= 0)
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

		case NFS_STATE_EXIT:
			// state DONE means we are done with nfs, which could return control
			// to factor().  But some states of finishing NFS are ill-suited to 
			// continuing with factor(), like if we detected an existing poly file
			// or data file and are refusing to continue.  In those cases user
			// intervention is required.  So we need to exit.
			exit(0);
			break;

		case NFS_STATE_FILTCHECK:
			if (job.current_rels >= job.min_rels)
			{
				if (fobj->VFLAG >= 0)
				{
					printf("nfs: found %u relations, need at least %u, proceeding with filtering ...\n",
						job.current_rels, job.min_rels);
				}
				
				nfs_state = NFS_STATE_FILTER;
			}
			else
			{
				// compute eta by dividing how many rels we have left to find
				// by the average time per relation.  we have average time
				// per relation because we've saved the time it took to do 
				// the last batch of sieving and we know how many relations we
				// found in that batch.
				uint32_t est_time;

				gettimeofday(&bstop, NULL);
                t_time = ytools_difftime(&bstart, &bstop);

				est_time = (uint32_t)((job.min_rels - job.current_rels) * 
					(t_time / (job.current_rels - pre_batch_rels)));				

				// if the user doesn't want to sieve, then we can't make progress.
				if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
					(fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE))
				{
					if (fobj->VFLAG >= 0)
					{
						printf("nfs: found %u relations, need at least %u "
							"(filtering ETA: %uh %um), continuing with sieving ...\n", // uh... um... hmm... idk *shrug*
							job.current_rels, job.min_rels, est_time / 3600,
							(est_time % 3600) / 60);
					}

					logprint_oc(fobj->flogname, "a", "nfs: found %u relations, need at least %u "
						"(filtering ETA: %uh %um), continuing with sieving ...\n", // uh... um... hmm... idk *shrug*
						job.current_rels, job.min_rels, est_time / 3600,
						(est_time % 3600) / 60);

					nfs_state = NFS_STATE_SIEVE;
				}
				else
				{
					if (fobj->VFLAG >= 0)
					{
						printf("nfs: found %u relations, need at least %u "
							"(filtering ETA: %uh %um), sieving not selected, finishing ...\n",
							job.current_rels, job.min_rels, est_time / 3600,
							(est_time % 3600) / 60);
					}

					nfs_state = NFS_STATE_DONE;
				}
			}

            if (fobj->VFLAG > 0)
            {
                gettimeofday(&stop, NULL);
                t_time = ytools_difftime(&start, &stop);
                printf("Elapsed time is now %6.4f seconds.\n", t_time);
            }

			break;

		case NFS_STATE_STARTNEW:

			nfs_state = NFS_STATE_POLY;		

			char fname[1024];
			sprintf(fname, "%s.ranges", fobj->nfs_obj.outputfile);
			FILE *fid = fopen(fname, "w");
			if (fid != NULL)
			{
				gmp_fprintf(fid, "%Zd\n", fobj->nfs_obj.gmp_n);
				fclose(fid);
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

			// should really be "resume_job", since we do more than just resume sieving...
		case NFS_STATE_RESUMESIEVE:

			// last_specialq == 0 if:
			// 1) user specifies -R and -ns with params
			// 2) user specifies post processing steps only
			// 3) user wants to resume sieving (either with a solo -ns or no arguments)
			//		but no data file or special-q was found
			// 4) -R was not specified (but then we won't be in this state, we'll be in DONE)
			// last_specialq > 1 if:
			// 5) user wants to resume sieving (either with a solo -ns or no arguments)
			//		and a data file and special-q was found
			// 6) it contains poly->time info (in which case we'll be in NFS_STATE_RESUMEPOLY)
            printf("nfs: resumesieve; last_spq = %u, nfs_phases = %u\n",
                last_specialq, fobj->nfs_obj.nfs_phases);

			strcpy(tmpstr, "");
			if ((last_specialq == 0) &&
				((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE)))
 			{				
                uint32_t missing_params;
                int betterskew;
				// this if-block catches cases 1 and 3 from above

                // init the poly: the job file parser will try to read
                // it if possible
                job.poly = (mpz_polys_t*)malloc(sizeof(mpz_polys_t));
                if (job.poly == NULL)
                {
                    printf("nfs: couldn't allocate memory!\n");
                    exit(-1);
                }
                mpz_polys_init(job.poly);
                job.poly->rat.degree = 1; // maybe way off in the future this isn't true
                // assume gnfs for now
                job.poly->side = ALGEBRAIC_SPQ;
                job.poly->size = (double)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);

				missing_params = parse_job_file(fobj, &job);

                // if we read a snfs job, also put any poly info
                // into the snfs poly.
                if (job.snfs != NULL)
                {
                    int i;
                    copy_mpz_polys_t(job.poly, job.snfs->poly);
                    mpz_set(job.snfs->n, fobj->nfs_obj.gmp_n);
                    for (i = job.poly->alg.degree; i >= 0; i--)
                        mpz_set(job.snfs->c[i], job.snfs->poly->alg.coeff[i]);
                }
				
				// set min_rels and other parameters.
                // this can also effect what side we sieve on and the
                // siever version.
                betterskew = get_ggnfs_params(fobj, &job);
                if ((job.snfs != NULL) && betterskew)
                {
                    missing_params |= PARAM_FLAG_IMPROVED_SKEW;
                }

                if ((job.snfs != NULL) && (missing_params > 0))
                {
                    // if the poly needs it, adjust any parameters that we filled.
                    skew_snfs_params(fobj, &job);
                }
				
				fill_job_file(fobj, &job, missing_params);
				// if any ggnfs params are missing, fill them
				// this means the user can provide an SNFS poly or external GNFS poly, 
				// but let YAFU choose the other params
				// this won't override any params in the file.

				char side = (job.poly->side == ALGEBRAIC_SPQ) ? 'a' : 'r';
				qrange_data_t* qrange_data = sort_completed_ranges(fobj, &job);
				qrange_t* next_range = get_next_range(qrange_data, side);

				if (job.current_rels >= job.min_rels)
				{
					nfs_state = NFS_STATE_FILTCHECK;
					continue;
				}
				
				if (next_range->qrange_start > 0)
				{
					sprintf(tmpstr, "nfs: resuming with sieving using ranges file.  "
						"Next range starts at special-q %u\n", next_range->qrange_start);
				}
				else if (fobj->nfs_obj.startq > 0)
				{
					job.startq = fobj->nfs_obj.startq;

					sprintf(tmpstr, "nfs: resuming with sieving at user specified special-q %u\n",
						job.startq);
				}
				else
				{
					// this is a guess, may be completely wrong
					//job.startq = (job.poly->side == RATIONAL_SPQ ? job.rlim : job.alim) / 2;

					sprintf(tmpstr, "nfs: continuing with sieving - could not determine "
						"last special q; using default startq = %u\n", job.startq);
				}

				// next step is sieving
				nfs_state = NFS_STATE_SIEVE;
			}
			else if ((last_specialq == 0) &&
				((fobj->nfs_obj.nfs_phases & NFS_PHASE_FILTER) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA_RESUME) ||
				(fobj->nfs_obj.nfs_phases & NFS_PHASE_SQRT)))
			{
				// this if-block catches case 2 from above
				// with these options we don't check for the last special-q, so this isn't
				// really a new factorization
				if ((fobj->nfs_obj.nfs_phases & NFS_PHASE_FILTER))
				{
					nfs_state = NFS_STATE_FILTCHECK;
					sprintf(tmpstr, "nfs: resuming with filtering\n");
				}
				else if ((fobj->nfs_obj.nfs_phases & NFS_PHASE_LA) ||
					(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA_RESUME))
				{
					nfs_state = NFS_STATE_LINALG;
					sprintf(tmpstr, "nfs: resuming with linear algebra\n");
				}
				else if (fobj->nfs_obj.nfs_phases & NFS_PHASE_SQRT)
				{
					nfs_state = NFS_STATE_SQRT;
					sprintf(tmpstr, "nfs: resuming with sqrt\n");
				}

			}
			else // data file already exists
			{				
				// this if-block catches case 5 from above
                // init the poly: the job file parser will try to read
                // it if possible
                job.poly = (mpz_polys_t*)malloc(sizeof(mpz_polys_t));
                if (job.poly == NULL)
                {
                    printf("nfs: couldn't allocate memory!\n");
                    exit(-1);
                }
                mpz_polys_init(job.poly);
                job.poly->rat.degree = 1; // maybe way off in the future this isn't true
                // assume gnfs for now
                job.poly->side = ALGEBRAIC_SPQ;
                job.poly->size = (double)mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);

				(void) parse_job_file(fobj, &job);
				
                // if we read a snfs job, also put any poly info
                // into the snfs poly.
                if (job.snfs != NULL)
                {
                    int i;
                    copy_mpz_polys_t(job.poly, job.snfs->poly);
                    mpz_set(job.snfs->n, fobj->nfs_obj.gmp_n);
                    for (i = job.poly->alg.degree; i >= 0; i--)
                        mpz_set(job.snfs->c[i], job.snfs->poly->alg.coeff[i]);
                }

                // set min_rels and other parameters.
                // this can also effect what side we sieve on and the
                // siever version.
                get_ggnfs_params(fobj, &job);

				if (fobj->nfs_obj.startq > 0)
				{
					// user wants to resume sieving at a specified point,
					// regardless of what the last found special-q may be.
					job.startq = fobj->nfs_obj.startq;
					nfs_state = NFS_STATE_SIEVE;
					sprintf(tmpstr, "nfs: resuming with sieving at special-q = %u\n",
						job.startq);
				}
				else
				{
					// resume at the last found special-q
					job.startq = last_specialq;

					// we found some relations - find the appropriate state
					// given user input
					if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
						(fobj->nfs_obj.nfs_phases & NFS_PHASE_FILTER))
					{
						nfs_state = NFS_STATE_FILTCHECK;
						sprintf(tmpstr, "nfs: resuming with filtering\n");
					}
					else if (fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE)
					{
						nfs_state = NFS_STATE_SIEVE;
						sprintf(tmpstr, "nfs: resuming with sieving at special-q = %u\n",
							last_specialq);
					}
					else if ((fobj->nfs_obj.nfs_phases & NFS_PHASE_LA) ||
						(fobj->nfs_obj.nfs_phases & NFS_PHASE_LA_RESUME))
					{
						nfs_state = NFS_STATE_LINALG;
						sprintf(tmpstr, "nfs: resuming with linear algebra\n");
					}
					else if (fobj->nfs_obj.nfs_phases & NFS_PHASE_SQRT)
					{
						nfs_state = NFS_STATE_SQRT;
						sprintf(tmpstr, "nfs: resuming with sqrt\n");
					}

				}				
			}

			if (fobj->VFLAG >= 0)
				printf("%s", tmpstr);

			logprint_oc(fobj->flogname, "a", "%s", tmpstr);

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
			if (fobj->VFLAG > 1) printf("nfs: resuming poly select\n");
			fobj->nfs_obj.polystart = job.last_leading_coeff;

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
        t_time = ytools_difftime(&start, &stop);

		if ((fobj->nfs_obj.timeout > 0) && (t_time > (double)fobj->nfs_obj.timeout))
		{
			if (fobj->VFLAG >= 0)
				printf("NFS timeout after %6.4f seconds.\n",t_time);

			logprint_oc(fobj->flogname, "a", "NFS timeout after %6.4f seconds.\n",t_time);
			process_done = 1;
		}

		if (NFS_ABORT)
		{
			if (obj != NULL)
				msieve_obj_free(obj);
			free(input);

			if (job.snfs)
			{
				snfs_clear(job.snfs);
				free(job.snfs);
			}
			else if (job.poly)
			{
				mpz_polys_free(job.poly);
				free(job.poly);
			}

			print_factors(fobj,fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES, fobj->OBASE);
			exit(1);
		}

		if (process_done && (nfs_state != NFS_STATE_DONE))
		{
			// if some state other than STATE_DONE set process_done, then
			// it had a good reason to.  Don't allow nfs() or factor() to continue, just leave.
			if (obj != NULL)
				msieve_obj_free(obj);
			free(input);

			if (job.snfs)
			{
				snfs_clear(job.snfs);
				free(job.snfs);
			}
			else if (job.poly)
			{
				mpz_polys_free(job.poly);
				free(job.poly);
			}

			print_factors(fobj, fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES, fobj->OBASE);
			exit(0);
		}
	}

	//reset signal handler to default (no handler).
	signal(SIGINT,NULL);

	gettimeofday(&stop, NULL);
	t_time = ytools_difftime(&start, &stop);
	fobj->nfs_obj.ttime = t_time;

	if (obj != NULL)
		msieve_obj_free(obj);
	free(input);
	
	if( job.snfs )
	{
		snfs_clear(job.snfs);
		free(job.snfs);
	} 
	else if( job.poly )
	{
		mpz_polys_free(job.poly);
		free(job.poly);
	}

	return;
}

int check_for_sievers(fact_obj_t *fobj, int revert_to_siqs)
{
	// if we are going to be doing sieving, check for the sievers
	if ((fobj->nfs_obj.nfs_phases == NFS_DEFAULT_PHASES) ||
		(fobj->nfs_obj.nfs_phases & NFS_PHASE_SIEVE))
	{
		FILE *test;
		char name[1024];
		int found, i;

		for (i=11; i<=16; i++)
		{
			sprintf(name, "%sggnfs-lasieve4I%de", fobj->nfs_obj.ggnfs_dir, i);
#if defined(WIN32)
			sprintf(name, "%s.exe", name);
#endif
			// test for existence of the siever
			test = fopen(name, "rb");
			if (test != NULL)
			{
				found = 1;
				fclose(test);
				break;
			}

            sprintf(name, "%sgnfs-lasieve4I%de", fobj->nfs_obj.ggnfs_dir, i);
#if defined(WIN32)
            sprintf(name, "%s.exe", name);
#endif
            // test for existence of the siever
            test = fopen(name, "rb");
            if (test != NULL)
            {
                found = 1;
                fclose(test);
                break;
            }
		}

		if (!found && revert_to_siqs)
		{
			printf("WARNING: could not find ggnfs sievers, reverting to siqs!\n");
			logprint_oc(fobj->flogname, "a", "WARNING: could not find ggnfs sievers, "
				"reverting to siqs!\n");

			mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
			SIQS(fobj);
			mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);
			return 1;
		}
		else if (!found)
			return 1;
		else
			return 0;
	}

	return 0;
}

int est_gnfs_size(nfs_job_t *job)
{
	if (job->snfs == NULL)
		return 999999;
	else if ((job->snfs->difficulty < 0) || (job->snfs->difficulty > 512))
		return 999999;
    else
    {
        if (job->snfs->sdifficulty > job->snfs->difficulty)
            return  0.56 * job->snfs->sdifficulty + 30;
        else
            return  0.56 * job->snfs->difficulty + 30;
    }
		
}

int est_gnfs_size_via_poly(snfs_t *job)
{
	if ((job->difficulty < 0) || (job->difficulty > 512))
		return 999999;
	else
		return  0.56*job->difficulty + 30;
}

double** gnfs_table;
double** snfs_table;
int gnfs_table_rows;
int snfs_table_rows;
int gnfs_table_rdy;
int snfs_table_rdy;


//entries based on statistics gathered from many factorizations done
//over the years by myself and others, and from here:
//http://www.mersenneforum.org/showthread.php?t=12365
// Note: this original table will become the new SNFS table.
// GNFS jobs will use the data contributed by Gimarel (below).
double ggnfs_table_orig[GGNFS_TABLE_ROWS][GGNFS_TABLE_COLS] = {
/* note: min_rels column is no longer used - it is equation based and	*/
/* is filled in by get_ggnfs_params					*/
/* columns:								*/
/* digits, rlim, alim, lpbr, lpba, mfbr, mfba, rlambda, alambda, siever, start-q, q-range, minrels */
    {75,  300000,   300000,   23, 23,  46, 46, 2.0, 2.0, 11, 150000, 1000   , 0},
    {80,  600000,   600000,   24, 24,  48, 48, 2.1, 2.1, 11, 300000, 1000   , 0},
    {85,  900000,   900000,   24, 24,  48, 48, 2.1, 2.1, 11, 450000, 1000   , 0},
	{90,  1200000,  1200000,  25, 25,  50, 50, 2.3, 2.3, 11, 600000, 1000   , 0},
	{95,  1500000,  1500000,  25, 25,  50, 50, 2.5, 2.5, 12, 750000, 2000   , 0},
	{100, 1800000,  1800000,  26, 26,  52, 52, 2.5, 2.5, 12, 900000, 2000   , 0},
	{105, 2500000,  2500000,  26, 26,  52, 52, 2.5, 2.5, 12, 1250000, 2000  , 0},
	{110, 3200000,  3200000,  26, 26,  52, 52, 2.5, 2.5, 13, 1600000, 4000  , 0},
	{115, 4500000,  4500000,  27, 27,  54, 54, 2.5, 2.5, 13, 2250000, 4000  , 0},
	{120, 5500000,  5500000,  27, 27,  54, 54, 2.5, 2.5, 13, 2750000, 4000  , 0},
	{125, 7000000,  7000000,  27, 27,  54, 54, 2.5, 2.5, 13, 3500000, 4000  , 0},
	{130, 9000000,  9000000,  28, 28,  56, 56, 2.5, 2.5, 13, 4500000, 8000  , 0},
	{135, 11500000, 11500000, 28, 28,  56, 56, 2.6, 2.6, 14, 5750000, 8000  , 0},
	{140, 14000000, 14000000, 28, 28,  56, 56, 2.6, 2.6, 14, 7000000, 8000  , 0},
	{145, 19000000, 19000000, 28, 28,  56, 56, 2.6, 2.6, 14, 9500000, 8000  , 0},
	{150, 25000000, 25000000, 29, 29,  58, 58, 2.6, 2.6, 14, 12500000, 16000, 0},
	{155, 32000000, 32000000, 29, 29,  58, 58, 2.6, 2.6, 14, 16000000, 16000, 0},	
	{160, 40000000, 40000000, 30, 30,  60, 60, 2.6, 2.6, 14, 20000000, 16000, 0},	// snfs 232
	{165, 49000000, 49000000, 30, 30,  60, 60, 2.6, 2.6, 14, 24500000, 16000, 0},	// 241
	{170, 59000000, 59000000, 31, 31,  62, 62, 2.6, 2.6, 14, 29500000, 32000, 0},	// 250
	{175, 70000000, 70000000, 31, 31,  62, 62, 2.6, 2.6, 15, 35000000, 32000, 0},	// 259
	{180, 82000000, 82000000, 31, 31,  62, 62, 2.6, 2.6, 15, 41000000, 32000, 0},	// 267
	{185, 100000000,100000000, 32, 32, 64, 64, 2.6, 2.6, 16, 50000000, 32000, 0}
};

// contributed by Gimarel
/* columns:								*/
/* digits, rlim, alim, lpbr, lpba, mfbr, mfba, rlambda, alambda, siever, start-q, q-range, minrels */
double ggnfs_table_Gimarel[GGNFS_TABLE_ROWS_NEW][GGNFS_TABLE_COLS_NEW] = {
	{91, 350000, 840000, 25, 25, 50, 50, 2.400, 2.400, 12, 210000, 2000, 1460000		},
	{92, 375000, 900000, 25, 25, 50, 50, 2.400, 2.400, 12, 225000, 2000, 1510000		},
	{93, 400000, 960000, 25, 25, 50, 50, 2.400, 2.400, 12, 240000, 2000, 1560000		},
	{94, 430000, 1080000, 25, 25, 50, 50, 2.400, 2.400, 12, 270000, 2000, 1670000		},
	{95, 470000, 1200000, 25, 25, 50, 50, 2.400, 2.400, 12, 300000, 2000, 1780000		},
	{96, 500000, 1360000, 26, 26, 52, 52, 2.423, 2.423, 12, 340000, 2000, 1890000		},
	{97, 530000, 1440000, 26, 26, 52, 52, 2.423, 2.423, 12, 360000, 2000, 2190000		},
	{98, 560000, 1520000, 26, 26, 52, 52, 2.423, 2.423, 12, 380000, 2000, 2490000		},
	{99, 590000, 1620000, 26, 26, 52, 52, 2.423, 2.423, 12, 405000, 2000, 2790000		},
	{100, 620000, 1720000, 26, 26, 52, 52, 2.423, 2.423, 12, 430000, 2000, 3090000		},
	{101, 650000, 1680000, 27, 27, 54, 54, 2.407, 2.407, 12, 420000, 2000, 4400000		},
	{102, 650000, 1780000, 27, 27, 54, 54, 2.407, 2.407, 12, 445000, 2000, 4600000		},
	{103, 650000, 1880000, 27, 27, 54, 54, 2.407, 2.407, 12, 470000, 2000, 4800000		},
	{104, 650000, 1980000, 27, 27, 54, 54, 2.407, 2.407, 12, 495000, 2000, 5000000		},
	{105, 650000, 2080000, 27, 27, 54, 54, 2.407, 2.407, 12, 520000, 2000, 5250000		},
	{106, 700000, 2000000, 29, 29, 58, 58, 2.500, 2.500, 12, 500000, 5000, 9200000		},
	{107, 700000, 2000000, 29, 29, 58, 58, 2.500, 2.500, 12, 500000, 5000, 9600000		},
	{108, 700000, 2000000, 29, 29, 58, 58, 2.500, 2.500, 12, 500000, 5000, 10000000		},
	{109, 700000, 2000000, 29, 29, 58, 58, 2.500, 2.500, 12, 500000, 5000, 10400000		},
	{110, 700000, 2000000, 29, 29, 58, 58, 2.500, 2.500, 12, 500000, 5000, 10800000		},
	{111, 700000, 2000000, 29, 29, 58, 58, 2.500, 2.500, 12, 500000, 5000, 11200000		},
	{112, 740000, 2240000, 29, 29, 58, 58, 2.500, 2.500, 12, 560000, 5000, 12000000		},
	{113, 780000, 2480000, 29, 29, 58, 58, 2.500, 2.500, 12, 620000, 5000, 12800000		},
	{114, 820000, 2720000, 29, 29, 58, 58, 2.500, 2.500, 12, 680000, 5000, 13600000		},
	{115, 860000, 2960000, 29, 29, 58, 58, 2.500, 2.500, 12, 740000, 5000, 14400000		},
	{116, 900000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 15200000		},
	{117, 1000000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 16200000	},
	{118, 1100000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 17200000	},
	{119, 1200000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 18200000	},
	{120, 1300000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 19200000	},
	{121, 1400000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 20200000	},
	{122, 1500000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 21200000	},
	{123, 1600000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 22200000	},
	{124, 1700000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 23200000	},
	{125, 1800000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 24200000	},
	{126, 1900000, 3200000, 29, 29, 58, 58, 2.500, 2.500, 12, 800000, 5000, 25200000	},
	{127, 2350000, 7680000, 29, 29, 58, 58, 2.500, 2.500, 13, 1920000, 5000, 24600000	},
	{128, 2500000, 7200000, 30, 30, 60, 60, 2.500, 2.500, 13, 1800000, 5000, 32000000	},
	{129, 2670000, 7880000, 30, 30, 60, 60, 2.500, 2.500, 13, 1970000, 5000, 34000000	},
	{130, 2840000, 8560000, 30, 30, 60, 60, 2.500, 2.500, 13, 2140000, 5000, 36000000	},
	{131, 2800000, 8000000, 30, 30, 60, 60, 2.533, 2.533, 13, 2000000, 5000, 40000000	},
	{132, 3080000, 8560000, 30, 30, 60, 60, 2.533, 2.533, 13, 2140000, 5000, 42000000	},
	{133, 3360000, 9120000, 30, 30, 60, 60, 2.533, 2.533, 13, 2280000, 5000, 44000000	},
	{134, 3640000, 9680000, 30, 30, 60, 60, 2.533, 2.533, 13, 2420000, 5000, 46000000	},
	{135, 3920000, 10240000, 30, 30, 60, 60, 2.533, 2.533, 13, 2560000, 5000, 48000000	},
	{136, 4200000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 50000000	},
	{137, 4460000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 52000000	},
	{138, 4720000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 54000000	},
	{139, 4980000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 56000000	},
	{140, 5240000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 58000000	},
	{141, 5500000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 60000000	},
	{142, 5760000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 62000000	},
	{143, 6020000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 64000000	},
	{144, 6280000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 66000000	},
	{145, 6540000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 68000000	},
	{146, 6800000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 70000000	},
	{147, 7060000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 72000000	},
	{148, 7320000, 10800000, 30, 30, 60, 60, 2.533, 2.533, 13, 2700000, 5000, 74000000	},

	{150, 25000000, 25000000, 29, 29,  58, 58, 2.6, 2.6, 14, 12500000, 16000, 0},
	{155, 32000000, 32000000, 29, 29,  58, 58, 2.6, 2.6, 14, 16000000, 16000, 0},
	{160, 40000000, 40000000, 30, 30,  60, 60, 2.6, 2.6, 14, 20000000, 16000, 0},	// snfs 232
	{165, 49000000, 49000000, 30, 30,  60, 60, 2.6, 2.6, 14, 24500000, 16000, 0},	// 241
	{170, 59000000, 59000000, 31, 31,  62, 62, 2.6, 2.6, 14, 29500000, 32000, 0},	// 250
	{175, 70000000, 70000000, 31, 31,  62, 62, 2.6, 2.6, 15, 35000000, 32000, 0},	// 259
	{180, 82000000, 82000000, 31, 31,  62, 62, 2.6, 2.6, 15, 41000000, 32000, 0},	// 267
	{185, 100000000,100000000, 32, 32, 64, 64, 2.6, 2.6, 16, 50000000, 32000, 0}
};

double* parse_params_file(char* filename, int *numrows)
{
	FILE* fid = fopen(filename, "r");
	double* table = NULL;

	if (fid == NULL)
	{
		printf("could not open %s to read\n", filename);
		return 0;
	}
	else
	{
		int status = 0;
		char line[1024];

		*numrows = 0;
		while (!feof(fid))
		{
			char* ptr = fgets(line, 1024, fid);

			if (ptr == NULL)
				break;

			if (strlen(line) < 10)
				continue;

			if ((line[0] == '/') && (line[1] == '/'))
				continue;

			// read 12 doubles from this line.  If sscanf indicates
			// it didn't read that many, abort and use the default table
			// after warning the user that something is wrong.
			double rl, al;
			uint32_t d, rlim, alim, lpbr, lpba, mfbr, mfba, sa, qr, sq, mr;
			int n = sscanf(line, "%u,%u,%u,%u,%u,%u,%u,%lf,%lf,%u,%u,%u,%u",
				&d, &rlim, &alim, &lpbr, &lpba, &mfbr, &mfba, &rl, &al, &sa, &sq, &qr, &mr);

			if (n == (GGNFS_TABLE_COLS - 1))
			{
				if ((*numrows) == 0)
					table = (double*)xmalloc(GGNFS_TABLE_COLS * sizeof(double));
				else
					table = (double*)xrealloc(table, ((*numrows) + 1) * 
						GGNFS_TABLE_COLS * sizeof(double));

				table[(*numrows) * GGNFS_TABLE_COLS + 0] = d;
				table[(*numrows) * GGNFS_TABLE_COLS + 1] = rlim;
				table[(*numrows) * GGNFS_TABLE_COLS + 2] = alim;
				table[(*numrows) * GGNFS_TABLE_COLS + 3] = lpbr;
				table[(*numrows) * GGNFS_TABLE_COLS + 4] = lpba;
				table[(*numrows) * GGNFS_TABLE_COLS + 5] = mfbr;
				table[(*numrows) * GGNFS_TABLE_COLS + 6] = mfba;
				table[(*numrows) * GGNFS_TABLE_COLS + 7] = rl;
				table[(*numrows) * GGNFS_TABLE_COLS + 8] = al;
				table[(*numrows) * GGNFS_TABLE_COLS + 9] = sa;
				table[(*numrows) * GGNFS_TABLE_COLS + 10] = sq;
				table[(*numrows) * GGNFS_TABLE_COLS + 11] = qr;
				table[(*numrows) * GGNFS_TABLE_COLS + 12] = 0;

				(*numrows)++;
			}
			else if (n == GGNFS_TABLE_COLS)
			{
				if ((*numrows) == 0)
					table = (double*)xmalloc(GGNFS_TABLE_COLS * sizeof(double));
				else
					table = (double*)xrealloc(table, ((*numrows) + 1) *
						GGNFS_TABLE_COLS * sizeof(double));

				table[(*numrows) * GGNFS_TABLE_COLS + 0] = d;
				table[(*numrows) * GGNFS_TABLE_COLS + 1] = rlim;
				table[(*numrows) * GGNFS_TABLE_COLS + 2] = alim;
				table[(*numrows) * GGNFS_TABLE_COLS + 3] = lpbr;
				table[(*numrows) * GGNFS_TABLE_COLS + 4] = lpba;
				table[(*numrows) * GGNFS_TABLE_COLS + 5] = mfbr;
				table[(*numrows) * GGNFS_TABLE_COLS + 6] = mfba;
				table[(*numrows) * GGNFS_TABLE_COLS + 7] = rl;
				table[(*numrows) * GGNFS_TABLE_COLS + 8] = al;
				table[(*numrows) * GGNFS_TABLE_COLS + 9] = sa;
				table[(*numrows) * GGNFS_TABLE_COLS + 10] = sq;
				table[(*numrows) * GGNFS_TABLE_COLS + 11] = qr;
				table[(*numrows) * GGNFS_TABLE_COLS + 12] = mr;

				(*numrows)++;
			}
			else
			{
				printf("Format error in params_file: %s, using built-in params table\n",
					filename);
				if ((*numrows) > 0)
				{
					free(table);
				}
				return NULL;
			}
		}

		fclose(fid);
		return table;
	}

	return NULL;
}

void set_ggnfs_tables(fact_obj_t* fobj)
{
	double* parsed_table = NULL;
	int parsed_rows;
	int i, j;

	// if a parameters file was provided, attempt to parse it.
	if (strlen(fobj->nfs_obj.params_file) > 0)
	{
		parsed_table = parse_params_file(fobj->nfs_obj.params_file, &parsed_rows);

		if (parsed_table != NULL)
		{
			gnfs_table = (double**)xmalloc(parsed_rows * sizeof(double*));
			for (i = 0; i < parsed_rows; i++)
			{
				gnfs_table[i] = (double*)xmalloc(GGNFS_TABLE_COLS * sizeof(double));
				for (j = 0; j < GGNFS_TABLE_COLS; j++)
				{
					if ((j == 7) || (j == 8))
						gnfs_table[i][j] = parsed_table[i * GGNFS_TABLE_COLS + j];
					else
						gnfs_table[i][j] = floor(parsed_table[i * GGNFS_TABLE_COLS + j] + 0.5);
				}
			}
			gnfs_table_rows = parsed_rows;

			if (fobj->VFLAG >= 0)
			{
				printf("nfs: successfully parsed %d nfs parameter entries from %s\n",
					parsed_rows, fobj->nfs_obj.params_file);
			}
		}
	}

	if (parsed_table == NULL)
	{
		// no user supplied table, use Gimarel's table (extended with orig table data > 148 digits).
		gnfs_table = (double**)xmalloc(GGNFS_TABLE_ROWS_NEW * sizeof(double*));
		for (i = 0; i < GGNFS_TABLE_ROWS_NEW; i++)
		{
			gnfs_table[i] = (double*)xmalloc(GGNFS_TABLE_COLS * sizeof(double));
			for (j = 0; j < GGNFS_TABLE_COLS; j++)
			{
				if ((j == 7) || (j == 8))
					gnfs_table[i][j] = ggnfs_table_Gimarel[i][j];
				else
					gnfs_table[i][j] = floor(ggnfs_table_Gimarel[i][j] + 0.5);
			}
		}

		gnfs_table_rows = GGNFS_TABLE_ROWS_NEW;
	}

	// snfs tables uses the original default table.
	snfs_table = (double**)xmalloc(GGNFS_TABLE_ROWS * sizeof(double*));
	for (i = 0; i < GGNFS_TABLE_ROWS; i++)
	{
		snfs_table[i] = (double*)xmalloc(GGNFS_TABLE_COLS * sizeof(double));
		for (j = 0; j < GGNFS_TABLE_COLS; j++)
		{
			if ((j == 7) || (j == 8))
				snfs_table[i][j] = ggnfs_table_orig[i][j];
			else
				snfs_table[i][j] = floor(ggnfs_table_orig[i][j] + 0.5);
		}
	}

	snfs_table_rows = GGNFS_TABLE_ROWS;

	gnfs_table_rdy = 1;
	snfs_table_rdy = 1;

	if (fobj->VFLAG > 2)
	{
		printf("gnfs table:\n");
		for (i = 0; i < gnfs_table_rows; i++)
		{
			for (j = 0; j < GGNFS_TABLE_COLS; j++)
			{
				if ((j == 7) || (j == 8))
					printf("%1.4f,\t", gnfs_table[i][j]);
				else
					printf("%d,\t", (int)floor(gnfs_table[i][j] + 0.5));
			}
			printf("\n");
		}
		printf("snfs table:\n");
		for (i = 0; i < snfs_table_rows; i++)
		{
			for (j = 0; j < GGNFS_TABLE_COLS; j++)
			{
				if ((j == 7) || (j == 8))
					printf("%1.4f,\t", snfs_table[i][j]);
				else
					printf("%d,\t", (int)floor(snfs_table[i][j] + 0.5));
			}
			printf("\n");
		}
	}

	if (parsed_table != NULL)
	{
		free(parsed_table);
	}

	if (fobj->VFLAG >= 0)
	{
		printf("nfs: gnfs parameters table has %d rows spanning %d-%d digits\n", 
			gnfs_table_rows, (int)gnfs_table[0][0], (int)gnfs_table[gnfs_table_rows - 1][0]);
		printf("nfs: snfs parameters table has %d rows spanning %d-%d digits\n", 
			snfs_table_rows, (int)snfs_table[0][0], (int)snfs_table[snfs_table_rows - 1][0]);
	}

	return;
}

int get_ggnfs_params(fact_obj_t *fobj, nfs_job_t *job)
{
	// based on the size/difficulty of the input number, determine "good" parameters
	// for the following: factor base limits, factor base large prime bound, trial
	// division fudge factors, trial division cutoffs, ggnfs lattice siever area, 
	// and expected number of relations needed.  linearly interpolate between table
	// entries.  keep last valid entry off the ends of the table.  This will produce
	// increasingly poor choices as one goes farther off the table, but you should be
	// doing things by hand by then anyway.
    // if we find a better skew for snfs polys then return 1 else return 0;
	uint32_t i, d = gmp_base10(fobj->nfs_obj.gmp_n);
	double scale;
	int found = 0;
	uint32_t lpb = 0, mfb = 0, fblim = 0, siever = 0;
	double lambda;
    int betterskew = 0;

	/*
	if (job->snfs == N && job->size != 0 && job->size != d && fobj->VFLAG > 0)
		printf("nfs: warning: size param in job file does not match size of "
			"number, ignoring param\n");	*/

	if (job->snfs != NULL)
	{
        if ((job->snfs->sdifficulty == 0) && (job->snfs->difficulty == 0) &&
            (fabs(mpz_get_d(job->snfs->poly->m)) > 0) && 
            (job->snfs->poly->alg.degree > 0))
        {
            if (fobj->VFLAG > 0)
            {
                //printf("nfs: detected snfs job but no snfs difficulty; "
                //    "assuming size of number is the snfs difficulty\n");

                printf("nfs: detected snfs job but no snfs difficulty\n");
                printf("nfs: using m and poly coefficients to compute difficulty\n");
            }

            job->snfs->difficulty = log10(fabs(mpz_get_d(job->snfs->poly->m))) +
                log10(mpz_get_d(job->snfs->poly->alg.coeff[job->snfs->poly->alg.degree]));

            if (fobj->VFLAG > 0)
            {
                printf("nfs: snfs difficulty: %lf\n", job->snfs->difficulty);
            }
        }
        else
        {
            printf("nfs: using provided snfs difficulty %lf\n", job->snfs->difficulty);
            //printf("nfs: sdifficulty %lf\n", job->snfs->sdifficulty);
            //gmp_printf("nfs: m = %Zd\n", job->snfs->poly->m);
            //printf("nfs: alg degree %u\n", job->snfs->poly->alg.degree);
            d = job->snfs->difficulty;
        }
            
        check_poly(job->snfs, fobj->VFLAG);

        if (fobj->VFLAG > 0)
        {
            printf("nfs: analyzing poly: \n");
        }

        approx_norms(job->snfs);
            
        if (fobj->VFLAG > 0)
        {
            printf("nfs: anorm = %le, rnorm = %le, size: %le, alpha = %le, murphyE = %le\n",
                job->snfs->anorm, job->snfs->rnorm, job->snfs->poly->size, 
                job->snfs->poly->alpha, job->snfs->poly->murphy);
        }

        // this can effect what side we sieve on and the
        // siever version.
        snfs_scale_difficulty(job->snfs, 1, fobj->VFLAG);
            
        if (fobj->VFLAG > 0)
        {
            printf("nfs: snfs skewed difficulty: %lf\n", job->snfs->sdifficulty);
			printf("nfs: snfs special-q side: %s\n", 
				job->poly->side == RATIONAL_SPQ ? "rational" : "algebraic");
        }

        d = job->snfs->sdifficulty;

        int do_skew_opt = 0;
        
        if (do_skew_opt && (job->snfs->poly->skew > 0))
        {
            // if requested, dither the skew to attempt to find 
            // a higher score.
            double bestmurph = job->snfs->poly->murphy;
            double bestskew = job->snfs->poly->skew;
            double origskew = job->snfs->poly->skew;
            double skew1percent = origskew * 0.01;
            int sign = 1;

            if (fobj->VFLAG > 0)
                printf("nfs: optimizing skew\n");

            i = 0;
            do
            {                
                job->snfs->poly->skew += skew1percent;
                
                //printf("on iteration %d trying skew %lf: ", i++, job->snfs->poly->skew);
                analyze_one_poly_xface(job->snfs);
                //printf("murphy = %le\n", job->snfs->poly->murphy);
                if (job->snfs->poly->murphy > bestmurph)
                {
                    bestskew = job->snfs->poly->skew;
                    bestmurph = job->snfs->poly->murphy;
                    if (fobj->VFLAG > 0)
                    {
                        printf("nfs: found better skew: %lf with murphy score %le\n",
                            bestskew, bestmurph);
                    }
                }
                else if ((job->snfs->poly->murphy < bestmurph) && (sign > 0))
                {
                    sign = -1;
                    job->snfs->poly->skew = origskew;
                    skew1percent *= -1.0;
                }
                else if ((job->snfs->poly->murphy < bestmurph) && (sign < 0))
                {
                    sign = 0;
                }
                i++;
                if (i >= 100)
                {
                    if (fobj->VFLAG > 0)
                    {
                        printf("nfs: giving up skew optimization after %u iterations\n", i);
                    }
                    break;
                }
            } while (sign != 0);

            if (bestskew > origskew)
            {
                betterskew = 1;
                job->snfs->poly->skew = bestskew;
            }
        }
		

		// http://www.mersenneforum.org/showpost.php?p=312701&postcount=2
		i = est_gnfs_size(job);
		if (fobj->VFLAG > 0)
			printf("nfs: guessing snfs difficulty %d is roughly equal to "
				"gnfs difficulty %d\n", d, i);
		d = i;
	}
    else if (fobj->nfs_obj.snfs)
    {
        printf("nfs: user passed snfs switch, but the job file does not specify snfs\n"
            "nfs: will continue as gnfs\n");
    }

	if (job->poly == NULL)
	{ // always be sure we can choose which side to sieve
		if (job->snfs != NULL)
		{
			job->poly = job->snfs->poly;
		}
		else
		{
			job->poly = (mpz_polys_t*)malloc(sizeof(mpz_polys_t));
			if (job->poly == NULL)
			{
				printf("nfs: couldn't allocate memory!\n");
				exit(-1);
			}
			mpz_polys_init(job->poly);
			job->poly->rat.degree = 1;
			// if we got in here without a poly, then it's probably a standard
			// gnfs resume (any snfs job would have already been detected and poly
			// properly set before calling this func)
			job->poly->side = ALGEBRAIC_SPQ;
		}
	}
	else if ((job->snfs != NULL) && (job->poly != job->snfs->poly))
	{ // *shrug* it could happen I suppose
		mpz_polys_free(job->poly);
		free(job->poly);
		job->poly = job->snfs->poly;
	}

    if (fobj->nfs_obj.sq_side != 0)
    {
        // user override
        job->poly->side = fobj->nfs_obj.sq_side > 0 ? ALGEBRAIC_SPQ : RATIONAL_SPQ;
		printf("nfs: user override of special-q side: to %s\n",
			job->poly->side == RATIONAL_SPQ ? "rational" : "algebraic");
    }

    if (fobj->VFLAG > 0)
    {
        printf("nfs: creating ggnfs job parameters for input of size %u\n", d);
        fflush(stdout);
    }

	double** table;
	int table_rows;

	if (job->snfs != NULL)
	{
		table = snfs_table;
		table_rows = snfs_table_rows;
	}
	else
	{
		table = gnfs_table;
		table_rows = gnfs_table_rows;
	}

	for (i=0; i< table_rows - 1; i++)
	{
		if (d > table[i][0] && d <= table[i+1][0])
		{
			scale = (double)(table[i+1][0] - d) /
				(double)(table[i+1][0] - table[i][0]);

			//interp
			fblim = table[i+1][1] - 
				(uint32_t)(scale * (double)(table[i+1][1] - table[i][1]));
			if (job->rlim == 0) job->rlim = fblim;

			fblim = table[i + 1][2] -
				(uint32_t)(scale * (double)(table[i + 1][2] - table[i][2]));
			if (job->alim == 0) job->alim = fblim;

			//pick closest entry
			if ((d - table[i][0]) < (table[i+1][0] - d))
				lpb = table[i][3];
			else
				lpb = table[i+1][3];
			if (job->lpbr == 0) job->lpbr = lpb;

			if ((d - table[i][0]) < (table[i + 1][0] - d))
				lpb = table[i][4];
			else
				lpb = table[i + 1][4];
			if (job->lpba == 0) job->lpba = lpb;

			//pick closest entry
			if ((d - table[i][0]) < (table[i+1][0] - d))
				mfb = table[i][5];
			else
				mfb = table[i+1][5];
			if (job->mfbr == 0) job->mfbr = mfb;

			if ((d - table[i][0]) < (table[i + 1][0] - d))
				mfb = table[i][6];
			else
				mfb = table[i + 1][6];
			if (job->mfba == 0) job->mfba = mfb;

			//pick closest entry
			if ((d - table[i][0]) < (table[i+1][0] - d))
				lambda = table[i][7];
			else
				lambda = table[i+1][7];
			if (job->rlambda < 1e-6) job->rlambda = lambda;

			if ((d - table[i][0]) < (table[i + 1][0] - d))
				lambda = table[i][8];
			else
				lambda = table[i + 1][8];
			if (job->alambda < 1e-6) job->alambda = lambda;

			//pick closest entry
			if ((d - table[i][0]) < (table[i+1][0] - d))
				siever = table[i][9];
			else
				siever = table[i+1][9];
			// if no user specified siever - pick one from the table
			if (fobj->nfs_obj.siever == 0) fobj->nfs_obj.siever = siever;

			//interp
			job->startq = table[i + 1][10] -
				(uint32_t)(scale * (double)(table[i + 1][10] - table[i][10]));

			//interp
			job->qrange = table[i+1][11] - 
				(uint32_t)(scale * (double)(table[i+1][11] - table[i][11]));

			//interp
			job->min_rels = table[i + 1][12] -
				(uint32_t)(scale * (double)(table[i + 1][12] - table[i][12]));

			found = 1;
		}
	}

	if (found == 0)
	{
		//couldn't find a table entry
		if (d <= table[0][0])
		{
			fblim = table[0][1];
			if (job->rlim == 0) job->rlim = fblim;

			fblim = table[0][2];
			if (job->alim == 0) job->alim = fblim;

			lpb = table[0][3];
			if (job->lpbr == 0) job->lpbr = lpb;

			lpb = table[0][4];
			if (job->lpba == 0) job->lpba = lpb;

			mfb = table[0][5];
			if (job->mfbr == 0) job->mfbr = mfb;

			mfb = table[0][6];
			if (job->mfba == 0) job->mfba = mfb;

			lambda = table[0][7];
			if (job->rlambda < 1e-6) job->rlambda = lambda;
			
			lambda = table[0][8];
			if (job->alambda < 1e-6) job->alambda = lambda;

			siever = table[0][9];
			if (fobj->nfs_obj.siever == 0) fobj->nfs_obj.siever = siever;

			job->startq = table[0][10];
			job->qrange = table[0][11];
			job->min_rels = table[0][12];
		}
		else
		{
			fblim = table[table_rows -1][1];
			if (job->rlim == 0) job->rlim = fblim;

			fblim = table[table_rows - 1][2];
			if (job->alim == 0) job->alim = fblim;

			lpb = table[table_rows -1][3];
			if (job->lpbr == 0) job->lpbr = lpb;

			lpb = table[table_rows - 1][4];
			if (job->lpba == 0) job->lpba = lpb;

			mfb = table[table_rows -1][5];
			if (job->mfbr == 0) job->mfbr = mfb;

			mfb = table[table_rows - 1][6];
			if (job->mfba == 0) job->mfba = mfb;

			lambda = table[table_rows -1][7];
			if (job->rlambda < 1e-6) job->rlambda = lambda;

			lambda = table[table_rows - 1][8];
			if (job->alambda < 1e-6) job->alambda = lambda;

			siever = table[table_rows -1][9];
			if (fobj->nfs_obj.siever == 0) fobj->nfs_obj.siever = siever;

			job->startq = table[table_rows - 1][10];
			job->qrange = table[table_rows - 1][11];
			job->min_rels = table[table_rows - 1][12];
		}
	}

	// swap parameters if side is not algebraic.  This assumes the param table data
	// was assembled for algebraic-side sieving.  
	if (0) //job->poly->side == RATIONAL_SPQ)
	{
		uint32_t tmp;
		double dtmp;

		tmp = job->rlim;
		job->rlim = job->alim;
		job->alim = tmp;

		tmp = job->lpbr;
		job->lpbr = job->lpba;
		job->lpba = tmp;

		tmp = job->mfbr;
		job->mfbr = job->mfba;
		job->mfba = tmp;

		dtmp = job->rlambda;
		job->rlambda = job->alambda;
		job->alambda = dtmp;
	}

	job->test_score = 9999999.0;		// haven't tested it yet.
	
	// if there is no min_rels in the table, use the equation
	nfs_set_min_rels(job);

	// command line switch override
	if (fobj->nfs_obj.minrels > 0)
	{
		if (fobj->VFLAG > 0)
		{
			logprint_oc(fobj->flogname, "a",
				"nfs: overriding default min_rels = %u with user supplied min_rels = %u\n",
				job->min_rels, fobj->nfs_obj.minrels);
			printf("nfs: overriding default min_rels = %u with user supplied min_rels = %u\n",
				job->min_rels, fobj->nfs_obj.minrels);
		}
		job->min_rels = fobj->nfs_obj.minrels;
	}

	if (fobj->nfs_obj.startq > 0)
	{
		// user specified startq, either by itself or as part of
		// a custom range.  override the table start-q
		if (fobj->VFLAG > 0)
		{
			logprint_oc(fobj->flogname, "a",
				"nfs: overriding default startq = %u with user supplied startq = %u\n",
				job->startq, fobj->nfs_obj.startq);
			printf("nfs: overriding default startq = %u with user supplied startq = %u\n",
				job->startq, fobj->nfs_obj.startq);
		}
		job->startq = fobj->nfs_obj.startq;
	}

	sprintf(job->sievername, "%sgnfs-lasieve4I%de", fobj->nfs_obj.ggnfs_dir, fobj->nfs_obj.siever);
#if defined(WIN32)
	sprintf(job->sievername, "%s.exe", job->sievername);
#endif

	return betterskew;
}

void trial_sieve(fact_obj_t* fobj)
{
	char** filenames = (char**)malloc(100*sizeof(char*));
	char* ptr, * arg = fobj->nfs_obj.filearg;
	int i = 0, me;

	if (fobj->VFLAG < 0) fobj->VFLAG = 0;

	while((ptr = strchr(arg, ','))) // this sort of thing is what's absolutely brilliant about C
	{
		filenames[i] = (char*)malloc(GSTR_MAXSIZE*sizeof(char));
        memset(filenames[i], 0, GSTR_MAXSIZE);
		//printf("pointer: %p\n", filenames[i]);
		strncpy(filenames[i++], arg, ptr-arg);
		arg = ptr + 1;
	}
	filenames[i] = (char*)malloc(GSTR_MAXSIZE*sizeof(char));
    memset(filenames[i], 0, GSTR_MAXSIZE);
	strcpy(filenames[i++], arg);

	me = test_sieve(fobj, filenames, i, 1);

	printf("test: \"%s\" is the fastest poly\n", filenames[me]);
	for(me = 0; me < i; me++)
		free(filenames[me]);
	free(filenames);
}

void nfs_set_min_rels(nfs_job_t *job)
{
	int i;
	double fudge; // sundae :)
	uint32_t lpb;

	if (job->min_rels == 0)
	{
		job->min_rels = 0;
		for (i = 0; i < 2; i++)
		{
			// these are always the same now... but someday maybe they won't be.  this
			// routine will handle that when we get there.
			// edit: we're there!!
			if (i == 0)
				lpb = job->lpbr;
			else
				lpb = job->lpba;

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

			job->min_rels += (uint32_t)(fudge * (
				pow(2.0, (double)lpb) / log(pow(2.0, (double)lpb))));

			//printf("nfs: lpb = %u, fudge = %lf, min_rels is now %u\n", lpb, fudge, job->min_rels);
		}
	}
	
}

void copy_mpz_polys_t(mpz_polys_t *src, mpz_polys_t *dest)
{
	int i;

	dest->skew = src->skew;
	dest->side = src->side;
	dest->murphy = src->murphy;
	mpz_set(dest->m, src->m);
	for (i=0; i<MAX_POLY_DEGREE+1; i++)
	{
		mpz_set(dest->rat.coeff[i], src->rat.coeff[i]);
		mpz_set(dest->alg.coeff[i], src->alg.coeff[i]);
	}
    dest->alg.degree = src->alg.degree;
    dest->rat.degree = src->rat.degree;

	return;
}

void copy_job(nfs_job_t *src, nfs_job_t *dest)
{
	if (dest->poly == NULL)
	{
		dest->poly = (mpz_polys_t *)malloc(sizeof(mpz_polys_t));
		mpz_polys_init(dest->poly);
	}

	copy_mpz_polys_t(src->poly, dest->poly);
	dest->rlim = src->rlim;
	dest->alim = src->alim;
	dest->lpbr = src->lpbr;
	dest->lpba = src->lpba;
	dest->mfbr = src->mfbr;
	dest->mfba = src->mfba;
	dest->rlambda = src->rlambda;
	dest->alambda = src->alambda;
	dest->qrange = src->qrange;
	strcpy(dest->sievername, src->sievername);
	dest->startq = src->startq;
	dest->min_rels = src->min_rels;
	dest->current_rels = src->current_rels;
	dest->poly_time = src->poly_time;
	dest->last_leading_coeff = src->last_leading_coeff;
	dest->use_max_rels = src->use_max_rels;
	dest->test_score = src->test_score;

	if (src->snfs != NULL)
	{
		if (dest->snfs == NULL)
		{
			dest->snfs = (snfs_t *)malloc(sizeof(snfs_t));
			snfs_init(dest->snfs);
		}
		snfs_copy_poly(src->snfs, dest->snfs);
	}
	return;
}

void mpz_polys_init(mpz_polys_t* poly) {
    mpz_poly_init(&poly->rat);
    mpz_poly_init(&poly->alg);
    mpz_init(poly->m);
    poly->skew = 0;
    poly->murphy = 0.;
    poly->size = 0.;
    poly->rroots = 0;
    poly->side = ALGEBRAIC_SPQ; // snfs routines will override if necessary
}

void mpz_polys_free(mpz_polys_t* poly) {
    mpz_poly_free(&poly->rat);
    mpz_poly_free(&poly->alg);
    mpz_clear(poly->m);
}

#else

void nfs(fact_obj_t *fobj)
{
	printf("nfs has not been enabled\n");

	mpz_set(fobj->qs_obj.gmp_n, fobj->nfs_obj.gmp_n);
	SIQS(fobj);
	mpz_set(fobj->nfs_obj.gmp_n, fobj->qs_obj.gmp_n);

	return;
}

#endif

