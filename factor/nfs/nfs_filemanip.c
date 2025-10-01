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
#include <stdlib.h>
#include "nfs_impl.h"
#include "gmp_xface.h"
#include <math.h>

#ifdef USE_NFS

void savefile_concat(char *filein, char *fileout, msieve_obj *mobj)
{
	FILE *in;

	in = fopen(filein,"r");
	if (in == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading\n",filein);
		exit(-1);
	}

	savefile_open(&mobj->savefile, SAVEFILE_APPEND);

	uint32_t count = 0;
	while (!feof(in))
	{
		char tmpline[GSTR_MAXSIZE], *tmpptr;
		tmpptr = fgets(tmpline, GSTR_MAXSIZE, in);
		if (tmpptr == NULL)
			break;
		else
		{
			savefile_write_line(&mobj->savefile, tmpline);
			count++;
		}
	}
	//printf("copied %u lines from %s to %s\n", count, filein, mobj->savefile.name);
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
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading\n",filein);
		exit(-1);
	}

	out = fopen(fileout,"a");
	if (out == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
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

enum nfs_state_e check_existing_files(fact_obj_t *fobj, uint32_t*last_spq, nfs_job_t *job)
{
	// see if we can resume a factorization based on the combination of input number,
	// .job file, .fb file, .dat file, and/or .p file.  else, start new job.
	// be very cautious about overwriting .dat files, as building these up can 
	// represent a *large* amount of work.
	// the only time .dat files should be touched is if the user forces us to (with -R)
	// or if an NFS process completes normally (in which case, we delete the .dat files).

	/*
	logic:

	1)
	.dat exists and no -R -> bail: insist user get rid of .dat or specify -R
	(this is first in order to prevent the case where there is a stale .dat file
	that will cause lots of -11 relation errors when a new job finishes, and also
	to prevent valid .dat files from being overwritten)

	2)
	job file missing, empty, malformed, or not matching -> check for poly in progress (3)
	else
	job file present and matching -> check for data file in progress (4)

	3)
	poly check:
	.p file missing, empty, malformed, or not matching -> start poly search at beginning
	else
	.p file present and matching -> get last leading coefficient and time
		return leading coefficient, put time in *last_spq.  calling routine
		will know what to do when it sees a large value returned

	4)
	data check:
	no data file -> start sieving at beginning
	data file present -> get last special q and count relations

	*/

	FILE *in, *logfile;
	char line[GSTR_MAXSIZE], *ptr;
	int do_poly_check, do_data_check;
	enum nfs_state_e ans;
	msieve_obj *mobj = fobj->nfs_obj.mobj;

	// 1) check for .dat file and resume flag
	if (savefile_exists(&mobj->savefile))
	{
		if (!fobj->nfs_obj.restart_flag)
		{
			printf("A data file (.dat) exists in the current directory.  It must either be "
				"removed or -R specified to resume nfs\n");

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
                    logprint(logfile, "nfs: refusing to resume without -R option\n");
                    fclose(logfile);
                }
            }

			// if we are inside factor, don't try to continue past this error
			fobj->flags |= FACTOR_INTERRUPT;
			*last_spq = 0;
			return NFS_STATE_DONE;
		}
	}

	do_data_check = 0;
	do_poly_check = 0;
	ans = NFS_STATE_STARTNEW;

	if (fobj->VFLAG > 0) printf("nfs: checking for job file - ");
	in = fopen(fobj->nfs_obj.job_infile,"r");
	if (in == NULL)
	{
		do_poly_check = 1;			// no job file.
		if (fobj->VFLAG > 0) printf("no job file found\n");
	}
	else
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
		{
			do_poly_check = 1;		// job file empty.
			if (fobj->VFLAG > 0) printf("job file empty\n");
			fclose(in);
		}
		else
		{
			if (fobj->VFLAG > 0) printf("job file found, testing for matching input\n");
			while (ptr != NULL)
			{
				if (line[0] == '#')
					ptr = fgets(line,GSTR_MAXSIZE,in);
				else if (line[0] != 'n')
				{
					do_poly_check = 1;	// malformed job file.
					if (fobj->VFLAG > 0) printf("nfs: malformed job file, first non-comment "
						"line should contain n: \n");
					break;
				}
				else
				{
					// check to see if the n in the file matches
					mpz_t tmp, g;
					mpz_init(tmp);
					mpz_init(g);

					mpz_set_str(tmp, line + 3, 0);
					if (resume_check_input_match(tmp, fobj->nfs_obj.gmp_n, g, fobj->VFLAG))
					{
						// divide out any common factor and copy the result to
						// other data structures
						if (mpz_cmp_ui(g, 1) > 0)
							add_to_factor_list(fobj->factors, g, fobj->VFLAG, fobj->NUM_WITNESSES);
						mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, g);
						mpz_set(fobj->N, fobj->nfs_obj.gmp_n);
                        char* s = mpz_get_str(NULL, 10, fobj->nfs_obj.gmp_n);
                        strcpy(fobj->nfs_obj.mobj->input, s);
                        free(s);

						do_data_check = 1;		// job file matches: check for data file
						if (fobj->VFLAG > 0) 
							printf("nfs: number in job file matches input\n");
					}
					else
					{
						if (fobj->VFLAG > 0)
							printf("nfs: number in job file does not match input\n");
						do_poly_check = 1;		// no match: check for poly file
					}
					mpz_clear(tmp);
					mpz_clear(g);
					break;
				}
			}
			if (ptr == NULL)
				do_poly_check = 1;		// malformed job file (empty)

			fclose(in);
		}
	}

	if (do_poly_check)
	{
		char master_polyfile[GSTR_MAXSIZE + 2];
		int do_poly_parse = 0;

		if (fobj->nfs_obj.nfs_phases != NFS_DEFAULT_PHASES && !(fobj->nfs_obj.nfs_phases & NFS_PHASE_POLY))
		{
			printf("nfs: error, no job file and polyselect not specified\n");
			return NFS_STATE_DONE;
		}

		if (fobj->VFLAG > 0)
			printf("nfs: checking for poly file - ");
		snprintf(master_polyfile, GSTR_MAXSIZE + 2, "%s.p",fobj->nfs_obj.outputfile);

		in = fopen(master_polyfile,"r");
		if (in == NULL)
		{
			if (fobj->VFLAG > 0) 
				printf("no poly file found\n");
			return NFS_STATE_STARTNEW;			// no .p file.
		}
		else
		{
			ptr = fgets(line,GSTR_MAXSIZE,in);
			if (ptr == NULL)
			{
				if (fobj->VFLAG > 0) 
					printf("poly file empty\n");
				return NFS_STATE_STARTNEW;		// .p file empty.
			}
			else
			{
				if (fobj->VFLAG > 0) printf("poly file found, testing for matching input\n");
				fclose(in);

				if (line[0] != 'n')
				{
					if (fobj->VFLAG > 0) 
						printf("nfs: malformed poly file, first line should contain n: \n");
					return NFS_STATE_STARTNEW;	// malformed .p file.
				}
				else
				{
					mpz_t tmp, g;
					mpz_init(tmp);
					mpz_init(g);
					mpz_set_str(tmp, line + 3, 0);
					if (resume_check_input_match(tmp, fobj->nfs_obj.gmp_n, g, fobj->VFLAG))
					{
						// divide out any common factor and copy the result to
						// other data structures
						if (mpz_cmp_ui(g, 1) > 0)
							add_to_factor_list(fobj->factors, g, fobj->VFLAG, fobj->NUM_WITNESSES);
						mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, g);
						mpz_set(fobj->N, fobj->nfs_obj.gmp_n);
                        char* s = mpz_get_str(NULL, 10, fobj->nfs_obj.gmp_n);
						strcpy(fobj->nfs_obj.mobj->input, s);
                        free(s);

						if (fobj->VFLAG > 0) 
							printf("nfs: poly file matches input\n");
						do_poly_parse = 1;		// .p file matches: parse .p file
					}
					mpz_clear(tmp);
					mpz_clear(g);
				}
			}
		}

		if (do_poly_parse)
		{
			if (fobj->nfs_obj.restart_flag)
			{
				if (fobj->VFLAG > 0) printf("nfs: finding best poly in poly file\n");
				find_best_msieve_poly(fobj, job, fobj->nfs_obj.job_infile, 0);
				*last_spq = job->poly_time;
				if (fobj->VFLAG > 0) printf("nfs: last leading coefficient was %u\n", 
					job->last_leading_coeff);
				return NFS_STATE_RESUMEPOLY;
			}
			else
			{
				printf("nfs: must specify -R to resume when a polyfile already exists\n");	

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
                        logprint(logfile, "nfs: refusing to resume poly select without -R option\n");
                        fclose(logfile);
                    }
                }

				*last_spq = 0;

				if (fobj->autofact_obj.autofact_active)
					return NFS_STATE_EXIT;
				else
					return NFS_STATE_DONE;
			}
		}
		else
			return NFS_STATE_STARTNEW;
	}
	else if (do_data_check)
	{
		// ok, we have a job file for the current input.  this is a restart of sieving
		ans = NFS_STATE_RESUMESIEVE;

		printf("nfs: checking for data file\n");
		// attempt to open data file
		if (!savefile_exists(&mobj->savefile))
		{
			// data file doesn't exist, return flag to start sieving from the beginning
            *last_spq = 0;
            if (fobj->VFLAG > 0)
            {
                printf("nfs: no data file found\n");
            }
		}
		else if (fobj->nfs_obj.restart_flag)
		{
			printf("nfs: previous data file found\n");
			// either restart from the end of the data file or from the specified 
			// sieve range
			//if (fobj->nfs_obj.sieve_only && (fobj->nfs_obj.startq > 0))
			if (fobj->nfs_obj.nfs_phases == NFS_PHASE_SIEVE && (fobj->nfs_obj.startq > 0))
			{
				if (fobj->VFLAG > 0)
					printf("nfs: user specified special-q range of %u-%u, "
					"skipping search for last special-q\n", 
					fobj->nfs_obj.startq, fobj->nfs_obj.startq + fobj->nfs_obj.rangeq);

				*last_spq = fobj->nfs_obj.startq;
				return NFS_STATE_RESUMESIEVE;
			}

			//if (fobj->nfs_obj.la_restart || fobj->nfs_obj.post_only)
			// if ( ! (fobj->nfs_obj.nfs_phases & ( ~(NFS_PHASE_FILTER | NFS_PHASE_LA | 
			//		NFS_PHASE_SQRT | NFS_PHASE_LA_RESUME) )) )
			if( fobj->nfs_obj.nfs_phases != NFS_DEFAULT_PHASES &&
				!(fobj->nfs_obj.nfs_phases & (NFS_PHASE_POLY | NFS_PHASE_SIEVE)))
			// if (not default) and not (poly or sieve)
			{
				if (fobj->VFLAG > 0)
					printf("nfs: user specified post processing only, "
					"skipping search for last special-q\n");

				*last_spq = 0;
				return NFS_STATE_RESUMESIEVE;
			}

			//if (fobj->VFLAG > 0)
			//	printf("nfs: commencing search for last special-q\n");

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
                    logprint(logfile, "nfs: previous data file found - "
                        "commencing search for last special-q\n");
                    fclose(logfile);
                }
            }

			// *last_spq = 0;
			// uint32_t totalrels;
			// qrange_data_t* qrange_data = sort_completed_ranges(fobj, &totalrels);
		}
		else
		{
			printf("nfs: must specify -R to resume when a savefile already exists\n");

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
                    logprint(logfile, "nfs: refusing to resume without -R option\n");
                    fclose(logfile);
                }
            }

			fobj->flags |= FACTOR_INTERRUPT;
			*last_spq = 0;
			

			if (fobj->autofact_obj.autofact_active)
				ans = NFS_STATE_EXIT;
			else
				ans = NFS_STATE_DONE;
		}

	}
	else
	{
		// job file is for a different input.  not a restart.
		*last_spq = 0;
		return NFS_STATE_STARTNEW;

	}

	return ans;

}

uint32_t get_spq(char **lines, int last_line, fact_obj_t *fobj)
{	
	//the last 4 valid lines are passed in
	int i, line, count;
	double var[2], avg[2];
	FILE *logfile;
	uint32_t rat[3], alg[3];
	char *ptr;

	if (fobj->VFLAG > 0)
		printf("nfs: parsing special-q from .dat file\n");

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
            logprint(logfile, "nfs: parsing special-q from .dat file\n");
            fclose(logfile);
        }
    }

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

		if (fobj->VFLAG > 1)
			printf("parsing rat side spq from %s\n",lines[line]);
		sscanf(lines[line] + i + 1, "%x", &rat[count]);
		
		if (fobj->VFLAG > 1)
			printf("found %x\n", rat[count]);

		// grab alg side entry
		for (i= strlen(lines[line]) - 1; i >= 0; i--)
			if (lines[line][i] == ',')
				break;

		if (fobj->VFLAG > 1)
			printf("parsing alg side spq from %s\n",lines[line]);
		sscanf(lines[line] + i + 1, "%x", &alg[count]);
		if (fobj->VFLAG > 1)
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

double find_best_msieve_poly(fact_obj_t *fobj, nfs_job_t *job, 
	char *jobfile_name,	int write_jobfile)
{
	// parse a msieve.dat.p file to find the best polynomial (based on e score)
	// and optionally output this as a ggnfs polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE + 2], *ptr;
	double score, bestscore = 0;
	int count, bestline = 0, i;
	uint32_t highest_c4 = 0, highest_c5 = 0;

	snprintf(line, GSTR_MAXSIZE + 2, "%s.p",fobj->nfs_obj.outputfile);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading!\n",line);
		return 0.0;
	}

	// read and count polys of the file
	// if (fobj->VFLAG > 0)
	// 	printf("searching for best poly...\n");

	count = 0;
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;		

		if (line[0] == '#')
		{
			count++;
			if (fobj->VFLAG > 2)
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
				}
			}
		}
		else if ((line[0] == 'c') && (line[1] == '4'))
		{
			uint32_t coeff;

			// get the c4 coefficient
			if (strlen(line) < 5)
				continue;	// not long enough to hold a line of format "c4: x"

			sscanf(line + 3, "%u", &coeff);
			if (coeff > highest_c4)
				highest_c4 = coeff;
		}
		else if ((line[0] == 'c') && (line[1] == '5'))
		{
			uint32_t coeff;

			// get the c5 coefficient
			if (strlen(line) < 5)
				continue;	// not long enough to hold a line of format "c5: x"

			sscanf(line + 3, "%u", &coeff);
			if (coeff > highest_c5)
				highest_c5 = coeff;
		}
		else if (line[0] == 't')
		{
			// get the time in seconds
			if (strlen(line) < 7)
				continue;	// not long enough to hold a line of format "time: x"

			//printf("found time record: %s",line);
			sscanf(line + 5, "%u", &job->poly_time);
		}		
	}
	fclose(in);

	if (fobj->VFLAG > 0)
		printf("nfs: best score = %0.4e, best poly = %d of %d\n", bestscore, bestline, count);

	if (highest_c5 > 0)
		job->last_leading_coeff = highest_c5;
	else 
		job->last_leading_coeff = highest_c4;

	if (!write_jobfile)
		return bestscore;

	// open it again
	snprintf(line, GSTR_MAXSIZE + 2, "%s.p",fobj->nfs_obj.outputfile);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading!\n",fobj->nfs_obj.outputfile);
		exit(1);
	}

	get_ggnfs_params(fobj, job);
	//job->startq = fobj->nfs_obj.sq_side < 0 ? job->rlim/2 : job->alim/2;
	// use alim if side not specified

	// always overwrites previous job files!
	//out = fopen(fobj->nfs_obj.job_infile,"w");
	out = fopen(jobfile_name, "w");
	if (out == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for writing!\n", jobfile_name);
		exit(1);
	}

	// go to the bestline
	char polyscore[256];
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;		

		if (line[0] == '#')
			bestline--;

		if (bestline == 0)
		{
			if (fobj->VFLAG > 0)
				printf("best poly: \n%s",line);

			strncpy(polyscore, line, 256);

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
                    logprint(logfile, "nfs: best poly = %s", line);
                    fclose(logfile);
                }
            }

			break;
		}
	}

	// copy n into the job file
	gmp_fprintf(out, "n: %Zd\n",fobj->nfs_obj.gmp_n);

	// and the score info
	fprintf(out, "%s", polyscore);

	if (fobj->VFLAG > 0)
		gmp_printf("n: %Zd\n",fobj->nfs_obj.gmp_n);

	logfile = NULL;
	if (fobj->LOGFLAG)
	{
		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
		{
			printf("error could not open logfile %s to append, using stdout\n", fobj->flogname);
			logfile = stdout;
		}
	}

	// read and record the poly so we can do some analysis on it.
	int max_deg = 0;
	snfs_t spoly;
	snfs_init(&spoly);

	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;

		// stop when we get to the next poly
		if (line[0] == '#')
			break;

		// prevent a "time" line from being copied into the .job file
		// ggnfs will ignore it, but it can be prevented here.
		if (strlen(line) > 4)
		{
			if ((line[0] == 't') && (line[1] == 'i') &&
				(line[2] == 'm') && (line[3] == 'e'))
				continue;
		}

		switch (line[0])
		{
		case 's':
			sscanf(line + 6, "%lf", &spoly.poly->skew);
			break;
		case 'c':
			switch (line[1])
			{
			case '0': mpz_set_str(spoly.poly->alg.coeff[0], line + 3, 10); max_deg = 0; break;
			case '1': mpz_set_str(spoly.poly->alg.coeff[1], line + 3, 10); max_deg = 1; break;
			case '2': mpz_set_str(spoly.poly->alg.coeff[2], line + 3, 10); max_deg = 2; break;
			case '3': mpz_set_str(spoly.poly->alg.coeff[3], line + 3, 10); max_deg = 3; break;
			case '4': mpz_set_str(spoly.poly->alg.coeff[4], line + 3, 10); max_deg = 4; break;
			case '5': mpz_set_str(spoly.poly->alg.coeff[5], line + 3, 10); max_deg = 5; break;
			case '6': mpz_set_str(spoly.poly->alg.coeff[6], line + 3, 10); max_deg = 6; break;
			default: break;
			}
			break;
		case 'Y':
			switch (line[1])
			{
			case '0': mpz_set_str(spoly.poly->rat.coeff[0], line + 3, 10); break;
			case '1': mpz_set_str(spoly.poly->rat.coeff[1], line + 3, 10); break;
			default: break;
			}
			break;
		default:
			break;
		}

		// don't write the poly here... we'll do that after
		// the whole thing is parsed into the poly structure and analyzed.
		//fputs(line, out);

		if (fobj->VFLAG > 0)
			printf("%s",line);
		
		if ((fobj->LOGFLAG) && (logfile != NULL))
		{
			logprint(logfile, "%s", line);
		}
	}

	// now that we have the poly in a structure, 
	// get norm estimates and print those too.
	mpz_set(spoly.n, fobj->nfs_obj.gmp_n);
	spoly.poly->alg.degree = max_deg;
	spoly.poly->rat.degree = 1;
	spoly.difficulty = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);
	// reverse calculation to get equivalent snfs difficulty for approx_norms...
	spoly.difficulty = (spoly.difficulty - 30) / 0.56;		
	spoly.valid = 1;

	approx_norms(&spoly);
	fprintf(out, "# anorm: %1.2e, rnorm: %1.2e\n", spoly.anorm, spoly.rnorm);

	// and print the poly to the file
	fprintf(out, "skew: %1.2f\n", spoly.poly->skew);
	for (i = 0; i <= max_deg; i++)
	{
		gmp_fprintf(out, "c%d: %Zd\n", i, spoly.poly->alg.coeff[i]);
	}
	gmp_fprintf(out, "Y0: %Zd\n", spoly.poly->rat.coeff[0]);
	gmp_fprintf(out, "Y1: %Zd\n", spoly.poly->rat.coeff[1]);

	// gnfs jobs don't need the poly info, it is all 
	// contained in the .job file.
	snfs_clear(&spoly);

	// and copy in the job parameters
	fprintf(out,"rlim: %u\n",job->rlim);
	fprintf(out,"alim: %u\n",job->alim);
	fprintf(out,"lpbr: %u\n",job->lpbr);
	fprintf(out,"lpba: %u\n",job->lpba);
	fprintf(out,"mfbr: %u\n",job->mfbr);
	fprintf(out,"mfba: %u\n",job->mfba);
	fprintf(out,"rlambda: %.4f\n",job->rlambda);
	fprintf(out,"alambda: %.4f\n",job->alambda);

	if ((fobj->LOGFLAG) && (logfile != NULL))
	{
		logprint(logfile, "rlim: %u\n", job->rlim);
		logprint(logfile, "alim: %u\n", job->alim);
		logprint(logfile, "lpbr: %u\n", job->lpbr);
		logprint(logfile, "lpba: %u\n", job->lpba);
		logprint(logfile, "mfbr: %u\n", job->mfbr);
		logprint(logfile, "mfba: %u\n", job->mfba);
		logprint(logfile, "rlambda: %.4f\n", job->rlambda);
		logprint(logfile, "alambda: %.4f\n", job->alambda);
	}

	fclose(in);
	fclose(out);

	if ((fobj->LOGFLAG) && (logfile != NULL))
	{
		fclose(logfile);
	}

	return bestscore;
}

void msieve_to_ggnfs(fact_obj_t *fobj, nfs_job_t *job)
{
	// convert a msieve.fb polynomial into a ggnfs polynomial file
	FILE *in, *out;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;

	in = fopen(fobj->nfs_obj.fbfile,"r");
	if (in == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading!\n",fobj->nfs_obj.fbfile);
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(fobj->nfs_obj.job_infile,"w");
	if (out == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for writing!\n",fobj->nfs_obj.job_infile);
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
	fprintf(out,"rlim: %u\n",job->rlim);
	fprintf(out,"alim: %u\n",job->alim);
	fprintf(out,"lpbr: %u\n",job->lpbr);
	fprintf(out,"lpba: %u\n",job->lpba);
	fprintf(out,"mfbr: %u\n",job->mfbr);
	fprintf(out,"mfba: %u\n",job->mfba);
	fprintf(out,"rlambda: %.1f\n",job->rlambda);
	fprintf(out,"alambda: %.1f\n",job->alambda);

	fclose(in);
	fclose(out);

	return;
}

void ggnfs_to_msieve(fact_obj_t *fobj, nfs_job_t *job)
{
	// convert a ggnfs.job file into a msieve.fb polynomial file
	FILE *in, *out;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;
	// sometimes reading rat.degree doesn't work... does it get set?
	char rats_printed = 2; //job->poly->rat.degree+1; 

	in = fopen(fobj->nfs_obj.job_infile,"r");
	if (in == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading!\n",fobj->nfs_obj.job_infile);
		exit(1);
	}

	//always overwrites previous job files!
	out = fopen(fobj->nfs_obj.fbfile,"w");
	if (out == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for writing!\n",fobj->nfs_obj.fbfile);
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
		else if (line[0] == 'Y' && rats_printed > 0) 
		{
			rats_printed--; 
			sprintf(outline, "R%c %s",line[1], line + 4);
		}
		else if (line[0] == 'c')
			sprintf(outline, "A%c %s",line[1], line + 4);
		else if (line[0] == 'm' && line[1] == ':' && rats_printed > 0) 
		{
			rats_printed=0;
			// need to do something different here for non-linear rational poly...
			sprintf(outline, "R1 1\nR0 -%s", line + 3);
		}
		else
			strcpy(outline, "");

		fputs(outline,out);
	}

	fclose(in);
	fclose(out);

	return;
}

int select_side_from_norms(double a, double r, FILE *flog, int VFLAG)
{
	int side;

	logprint(flog, "anorm = %1.2e, rnorm = %1.2e\n", a, r);

	if ((fabs(r / a) > 1e3)) // && (fobj->nfs_obj.sq_side == 0))
	{
		if (VFLAG > 0)
		{
			printf("nfs: found rnorm = %1.2e, anorm = %1.2e\n"
				"nfs: setting rational side sieve based on abs norm ratio r/a = %1.2e\n",
				r, a, fabs(r / a));
		}
		side = RATIONAL_SPQ;
		logprint(flog, "choosing rational side based on abs norm ratio r/a = %1.2e\n", fabs(r / a));
	}
	else // if (fobj->nfs_obj.sq_side == 0)
	{
		if (VFLAG > 0)
		{
			printf("nfs: found rnorm = %1.2e, anorm = %1.2e\n"
				"nfs: setting algebraic side sieve based on abs norm ratio a/r = %1.2e\n",
				r, a, fabs(a / r));
		}
		logprint(flog, "choosing algebraic side based on abs norm ratio a/r = %1.2e\n", fabs(a / r));
		side = ALGEBRAIC_SPQ;
	}

	return side;
}

uint32_t parse_job_file(fact_obj_t *fobj, nfs_job_t *job)
{
	FILE *in;
	uint32_t missing_params = 0;
	uint32_t lpbr = 0, lpba = 0, mfbr = 0, mfba = 0, alim = 0, rlim = 0, size = 0, siever = 0;
    int y0 = 0, y1 = 0, has_m = 0;
	char line[1024];
	float alambda = 0.0, rlambda = 0.0, skew = 0.0;
    int info1 = 0, info2 = 0, info3 = 0;
    int is_snfs = 0;
	enum special_q_e side = NEITHER_SPQ;
	
	
	in = fopen(fobj->nfs_obj.job_infile, "r");
	if (in == NULL)
	{
		printf("nfs: couldn't open job file, using default min_rels\n");
		return 0;
	}

	FILE* flog = fopen(fobj->flogname, "a");

	if (flog != NULL)
	{
		logprint(flog, "parsed parameters from job file %s\n", fobj->nfs_obj.job_infile);
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
			logprint(flog, "lpbr: %u\n", lpbr);
			continue;
		}

		substr = strstr(line, "lpba:");

		if (substr != NULL)
		{
			lpba = strtoul(substr + 5, NULL, 10);
			logprint(flog, "lpba: %u\n", lpba);
			continue;
		}

		substr = strstr(line, "type:");		
		if (substr != NULL)
 		{
 			if (strstr(substr + 5, "snfs")) // case sensitive
 			{
 				job->snfs = malloc(sizeof(snfs_t));
                is_snfs = 1;
 				if (job->snfs == NULL)
 				{
 					printf("nfs: couldn't allocate memory!\n");
 					exit(-1);
 				} 
                else if (fobj->VFLAG > 0)
                {
                    printf("nfs: found type: snfs\n");
                    printf("nfs: initializing snfs poly\n");
                }
				logprint(flog, "type: snfs\n");
 				snfs_init(job->snfs);
 			}
			continue;
		}

        {
            int coeff;
            int found = 0;
            for (coeff = 0; coeff < MAX_POLY_DEGREE; coeff++)
            {
                char cstr[8];
                sprintf(cstr, "c%d:", coeff);
                substr = strstr(line, cstr);
                if (substr != NULL)
                {
                    gmp_sscanf(line + 3, "%Zd", job->poly->alg.coeff[coeff]);
                    if (fobj->VFLAG > 0)
                    {
                        gmp_printf("nfs: found c[%d]: %Zd\n", coeff, job->poly->alg.coeff[coeff]);
                    }
                    found = 1;
                    if ((coeff > job->poly->alg.degree) &&
                        (mpz_cmp_ui(job->poly->alg.coeff[coeff], 0) > 0))
                    {
                        job->poly->alg.degree = coeff;
                    }
					logprint(flog, "c%d: %s", coeff, line + 3);
                    break;
                }
            }
            if (found)
                continue;
        }

        substr = strstr(line, "Y0:");
        if (substr != NULL)
        {
            gmp_sscanf(line + 3, "%Zd", job->poly->rat.coeff[0]);
            if (fobj->VFLAG > 0)
            {
                gmp_printf("nfs: found Y0: %Zd\n", job->poly->rat.coeff[0]);
            }
			logprint(flog, "Y0: %s", line + 3);
            y0 = 1;
            continue;
        }

        substr = strstr(line, "Y1:");
        if (substr != NULL)
        {
            gmp_sscanf(line + 3, "%Zd", job->poly->rat.coeff[1]);
            if (fobj->VFLAG > 0)
            {
                gmp_printf("nfs: found Y1: %Zd\n", job->poly->rat.coeff[1]);
            }
			logprint(flog, "Y1: %s", line + 3);
            job->poly->rat.degree = 1;
            y1 = 1;
            continue;
        }

        
        if ((line[0] == 'm') && (line[1] == ':'))
        {
            gmp_sscanf(line + 2, "%Zd", job->poly->m);
            if (fobj->VFLAG > 0)
            {
                gmp_printf("nfs: found m: %Zd\n", job->poly->m);
            }
			logprint(flog, "m: %s", line + 2);
            has_m = 1;
            continue;
        }

        substr = strstr(line, "skew:");
        if (substr != NULL)
        {
            sscanf(line + 5, "%lf", &job->poly->skew);
            skew = job->poly->skew;
            if (fobj->VFLAG > 0)
            {
                printf("nfs: found skew: %lf\n", job->poly->skew);
            }
			logprint(flog, "skew: %1.2f\n", skew);
        }
		
		substr = strstr(line, "size:");
		if (substr != NULL)
		{
			uint32_t difficulty = size = strtoul(substr + 5, NULL, 10);

            if (is_snfs)
            {
                job->snfs->difficulty = (double)difficulty;
            }
            else
            {
                printf("nfs: warning: found size line for non-snfs job.\n"
                    "nfs: if this is a SNFS job, type=snfs should appear prior to the size line\n");
            }

            if (fobj->VFLAG > 0)
            {
                printf("nfs: found size: %u\n", difficulty);
            }
			logprint(flog, "size: %u\n", difficulty);

			continue;
		}
	
		substr = strstr(line, "mfbr:");
		if (substr != NULL)
		{
			mfbr = strtoul(substr + 5, NULL, 10);
			logprint(flog, "mfbr: %u\n", mfbr);
			continue;
		}
		
		substr = strstr(line, "mfba:");
		if (substr != NULL)
		{
			mfba = strtoul(substr + 5, NULL, 10);
			logprint(flog, "mfba: %u\n", mfba);
			continue;
		}
		
		substr = strstr(line, "rlim:");		
		if (substr != NULL)
		{
			rlim = strtoul(substr + 5, NULL, 10);
			logprint(flog, "rlim: %u\n", rlim);
			continue;
		}
		
		substr = strstr(line, "alim:");		
		if (substr != NULL)
		{
			alim = strtoul(substr + 5, NULL, 10);
			logprint(flog, "alim: %u\n", alim);
			continue;
		}
		
		substr = strstr(line, "rlambda:");		
		if (substr != NULL)
		{
			sscanf(substr + 8, "%f", &rlambda); //strtof(substr + 8, NULL);
			logprint(flog, "rlambda: %1.4f\n", rlambda);
			continue;
		}
		
		substr = strstr(line, "alambda:");		
		if (substr != NULL)
		{
			sscanf(substr + 8, "%f", &alambda); //strtof(substr + 8, NULL);
			logprint(flog, "alambda: %1.4f\n", alambda);
			continue;
		}

		substr = strstr(line, "# siever:");		
		if (substr != NULL)
		{
			sscanf(substr + 9, "%u", &siever);
			if (fobj->VFLAG > 0)
				printf("nfs: found siever gnfs-lasieve4I%ue\n", siever);
			logprint(flog, "siever: %u\n", siever);
			continue;
		}
		
		substr = strstr(line, "algebraic");
		if (substr != NULL)
		{
			side = ALGEBRAIC_SPQ;
            info2 = 1;
			if (fobj->VFLAG > 0)
				printf("nfs: found side: algebraic\n");
			logprint(flog, "side: algebraic\n");
			continue;
		}
		
		substr = strstr(line, "rational");
		if (substr != NULL)
		{
			side = RATIONAL_SPQ;
            info2 = 1;
			if (fobj->VFLAG > 0)
				printf("nfs: found side: rational\n");
			logprint(flog, "side: rational\n");
			continue;
		}

        substr = strstr(line, "combined");
        if (substr != NULL)
        {
			double e;
			int nume = sscanf(substr + 9, "%lf", &e);
            info3 = 1;

            if (fobj->VFLAG > 0)
                printf("nfs: found murphyE score\n");

			if (nume == 1)
				logprint(flog, "combined score (MurphyE): %1.4e\n", e);
            continue;
        }

        substr = strstr(line, "anorm");
        if (substr != NULL)
        {
			double a, r;
			int numa, numr;
			numa = sscanf(substr + 6, "%lf", &a);
			
			substr = strstr(line, "rnorm");
			numr = sscanf(substr + 6, "%lf", &r);

			if ((numa < 0) || (numr < 0))
			{
				if (fobj->VFLAG > 0)
				{
					printf("nfs: found norm line but could not parse norm info\n"
						"nfs: using default algebraic side sieve\n");
				}
				side = ALGEBRAIC_SPQ;
			}
			else
			{
				side = select_side_from_norms(a, r, flog, fobj->VFLAG);
			}
            info1 = 1;
            
            continue;
        }
	}

	if (flog != NULL)
	{
		fclose(flog);
	}

    if (has_m)
    {

    }
    else if (y0 && y1)
    {
        // generate an m from the rational poly if we have the coefficients.
        // This does not take into account algebraic factors, if any.
        if (mpz_sgn(job->poly->rat.coeff[1]) < 0)
        {
            mpz_mul_si(job->poly->rat.coeff[1], job->poly->rat.coeff[1], -1);
            mpz_invert(job->poly->m, job->poly->rat.coeff[1], fobj->nfs_obj.gmp_n);
            mpz_mul(job->poly->m, job->poly->m, job->poly->rat.coeff[0]);
            //mpz_mod(job->poly->m, job->poly->m, fobj->nfs_obj.gmp_n);
            mpz_mul_si(job->poly->rat.coeff[1], job->poly->rat.coeff[1], -1);
        }
        else
        {
            mpz_mul_si(job->poly->rat.coeff[0], job->poly->rat.coeff[0], -1);
            mpz_invert(job->poly->m, job->poly->rat.coeff[1], fobj->nfs_obj.gmp_n);
            mpz_mul(job->poly->m, job->poly->m, job->poly->rat.coeff[0]);
            //mpz_mod(job->poly->m, job->poly->m, fobj->nfs_obj.gmp_n);
            mpz_mul_si(job->poly->rat.coeff[0], job->poly->rat.coeff[0], -1);
        }
    }

	if (siever > 0)
	{
		// if no command line siever specified
		if (fobj->nfs_obj.siever == 0)
		{
			// so get_ggnfs_params doesn't overwrite it
			fobj->nfs_obj.siever = siever;
			sprintf(job->sievername, "gnfs-lasieve4I%ue", siever);
#ifdef WIN32
			sprintf(job->sievername, "%s.exe", job->sievername);
#endif
		}

	}

	if (lpbr > 0)
	{
		if (job->lpbr == 0)
			job->lpbr = lpbr;
	}
    else
    {
        missing_params |= PARAM_FLAG_LPBR;
    }

	if (lpba > 0)
	{
		if (job->lpba == 0)
			job->lpba = lpba;
	}
    else
    {
        missing_params |= PARAM_FLAG_LPBA;
    }

    if (fobj->VFLAG > 0)
    {
        printf("nfs: parsed lpbr = %u, lpba = %u\n", lpbr, lpba);
    }

	if (mfbr > 0)
	{
		if (job->mfbr == 0)
			job->mfbr = mfbr;
	}
    else
    {
        missing_params |= PARAM_FLAG_MFBR;
    }

	if (mfba > 0)
	{
		if (job->mfba == 0)
			job->mfba = mfba;
	}
    else
    {
        missing_params |= PARAM_FLAG_MFBA;
    }


	if (rlim > 0)
	{
		if (job->rlim == 0)
			job->rlim = rlim;
	}
    else
    {
        missing_params |= PARAM_FLAG_RLIM;
    }

	if (alim > 0)
	{
		if (job->alim == 0)
			job->alim = alim;
	}
    else
    {
        missing_params |= PARAM_FLAG_ALIM;
    }


	if (rlambda > 0)
	{
		if (job->rlambda == 0)
			job->rlambda = rlambda;
	}
    else
    {
        missing_params |= PARAM_FLAG_RLAMBDA;
    }

	if (alambda > 0)
	{
		if (job->alambda == 0)
			job->alambda = alambda;
	}
    else
    {
        missing_params |= PARAM_FLAG_ALAMBDA;
    }

    if (info1 == 0)
    {
        missing_params |= PARAM_FLAG_INFO1;
    }

    if (info2 == 0)
    {
        missing_params |= PARAM_FLAG_INFO2;
    }

    if (info3 == 0)
    {
        missing_params |= PARAM_FLAG_INFO3;
    }

	if (size > 0)
	{
        if (is_snfs)
        {
            job->snfs->sdifficulty = job->snfs->difficulty = size;
        }
        else
        {
            printf("nfs: found a size parameter but not snfs type: gnfs job will ignore\n");
        }
	}
    else if (is_snfs)
    {
        // default for snfs jobs with no supplied size info.  
        // we will try to guess it later (in get_gnfs_params)
        job->snfs->sdifficulty = job->snfs->difficulty = 0;
    }
	
	if (side != NEITHER_SPQ)
	{
		if (job->snfs != NULL)
		{
			job->snfs->poly->side = side;
		}
		else
		{
			if (job->poly == NULL)
			{ 
				// always be sure we can choose which side to sieve
				job->poly = (mpz_polys_t*)malloc(sizeof(mpz_polys_t));
				if (job->poly == NULL)
				{
					printf("nfs: couldn't allocate memory!\n");
					exit(-1);
				}
				mpz_polys_init(job->poly);
				job->poly->rat.degree = 1;
			}
			
			job->poly->side = side;
		}		
	}

	fclose(in);

	return missing_params;
}
	
void fill_job_file(fact_obj_t *fobj, nfs_job_t *job, uint32_t missing_params)
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
		
        if (fobj->VFLAG > 0)
			printf("nfs: job file is missing params, filling them\n");

		// make sure we start on a new line if we are filling anything
		fprintf(out, "\n");

		if (missing_params & PARAM_FLAG_RLIM)
			fprintf(out,"rlim: %u\n",job->rlim);
		if (missing_params & PARAM_FLAG_ALIM)
			fprintf(out,"alim: %u\n",job->alim);
		
		if (missing_params & PARAM_FLAG_LPBR)
			fprintf(out,"lpbr: %u\n",job->lpbr);
		if (missing_params & PARAM_FLAG_LPBA)
			fprintf(out,"lpba: %u\n",job->lpba);
		
		if (missing_params & PARAM_FLAG_MFBR)
			fprintf(out,"mfbr: %u\n",job->mfbr);
		if (missing_params & PARAM_FLAG_MFBA)
			fprintf(out,"mfba: %u\n",job->mfba);
		
		if (missing_params & PARAM_FLAG_RLAMBDA)
			fprintf(out,"rlambda: %.4lf\n",job->rlambda);
		if (missing_params & PARAM_FLAG_ALAMBDA)
			fprintf(out,"alambda: %.4lf\n",job->alambda);

		if (missing_params & PARAM_FLAG_INFO1)
		{
			if (job->snfs == NULL)
			{
				snfs_t spoly;

				printf("filling info1: norms and side for gnfs deg %d,%d poly\n",
					job->poly->alg.degree, job->poly->rat.degree); fflush(stdout);

				snfs_init(&spoly);
				mpz_set(spoly.n, fobj->nfs_obj.gmp_n);
				copy_mpz_polys_t(job->poly, spoly.poly);
				spoly.difficulty = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);
				// reverse calculation to get equivalent snfs difficulty for approx_norms...
				spoly.difficulty = (spoly.difficulty - 30) / 0.56;
				spoly.valid = 1;

				approx_norms(&spoly);
				fprintf(out, "# anorm: %1.2e, rnorm: %1.2e\n", spoly.anorm, spoly.rnorm);

				FILE* flog = fopen(fobj->flogname, "a");
				if (flog == NULL)
				{
					flog = stdout;
				}

				job->poly->side = select_side_from_norms(spoly.anorm, spoly.rnorm, flog, fobj->VFLAG);

				if (flog != stdout)
				{
					fclose(flog);
				}

				snfs_clear(&spoly);
			}
			else
			{
				printf("filling info1: norms and side for snfs deg %d,%d poly\n",
					job->snfs->poly->alg.degree, job->snfs->poly->rat.degree); fflush(stdout);

				approx_norms(job->snfs);

				fprintf(out, "# anorm: %1.2e, rnorm: %1.2e\n", job->snfs->anorm, job->snfs->rnorm);

				FILE* flog = fopen(fobj->flogname, "a");
				if (flog == NULL)
				{
					flog = stdout;
				}
				job->poly->side = select_side_from_norms(job->snfs->anorm, job->snfs->rnorm, flog, fobj->VFLAG);
				if (flog != stdout)
				{
					fclose(flog);
				}
			}

		}
		
		fclose(out);
	}
}

void print_poly(mpz_polys_t* poly, FILE *out)
{
	// print the poly to stdout
	int i;
	fprintf(out, "skew: %1.4f\n", poly->skew);
	for (i=MAX_POLY_DEGREE; i>=0; i--)
		if (mpz_cmp_si(poly->alg.coeff[i], 0) != 0) 
			gmp_fprintf(out, "c%d: %Zd\n", i, poly->alg.coeff[i]);
	gmp_fprintf(out, "Y1: %Zd\n", poly->rat.coeff[1]);
	gmp_fprintf(out, "Y0: %Zd\n", poly->rat.coeff[0]);
    // if( mpz_cmp_si(poly->m, 0) != 0 ) gmp_fprintf(out, "m: %Zd\n", poly->m);
}

void print_job(nfs_job_t *job, FILE *out)
{
	mpz_polys_t* poly = job->poly;

	// print the poly to the supplied file stream
	int i;
	fprintf(out, "skew: %1.4f\n", poly->skew);
	for (i=MAX_POLY_DEGREE; i>=0; i--)
		if (mpz_cmp_si(poly->alg.coeff[i], 0) != 0) 
			gmp_fprintf(out, "c%d: %Zd\n", i, poly->alg.coeff[i]);
	gmp_fprintf(out, "Y1: %Zd\n", poly->rat.coeff[1]);
	gmp_fprintf(out, "Y0: %Zd\n", poly->rat.coeff[0]);
	// if( mpz_cmp_si(poly->m, 0) != 0 ) gmp_fprintf(out, "m: %Zd\n", poly->m);
	fprintf(out, "rlim: %u\nalim: %u\n", job->rlim, job->alim);
	fprintf(out, "mfbr: %u\nmfba: %u\n", job->mfbr, job->mfba);
	fprintf(out, "lpbr: %u\nlpba: %u\n", job->lpbr, job->lpba);
	fprintf(out, "rlambda: %1.2f\nalambda: %1.2f\n", job->rlambda, job->alambda);
}

#endif


