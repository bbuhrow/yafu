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

#include "nfs.h"
#include "gmp_xface.h"

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


enum nfs_state_e check_existing_files(fact_obj_t *fobj, uint32 *last_spq, ggnfs_job_t *job)
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

			*last_spq = 0;
			return NFS_STATE_DONE;
		}
	}

	do_data_check = 0;
	do_poly_check = 0;
	ans = NFS_STATE_STARTNEW;

	if (VFLAG > 0) printf("nfs: checking for job file - ");
	in = fopen(fobj->nfs_obj.job_infile,"r");
	if (in == NULL)
	{
		do_poly_check = 1;			// no job file.
		if (VFLAG > 0) printf("no job file found\n");
	}
	else
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
		{
			do_poly_check = 1;		// job file empty.
			if (VFLAG > 0) printf("job file empty\n");
			fclose(in);
		}
		else
		{
			while (ptr != NULL)
			{
				if (line[0] == '#')
					ptr = fgets(line,GSTR_MAXSIZE,in);
				else if (line[0] != 'n')
				{
					do_poly_check = 1;	// malformed job file.
					if (VFLAG > 0) printf("malformed job file, first non-comment "
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
					if (resume_check_input_match(tmp, fobj->nfs_obj.gmp_n, g))
					{
						// divide out any common factor and copy the result to
						// other data structures
						if (mpz_cmp_ui(g, 1) > 0)
							add_to_factor_list(fobj, g);
						mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, g);
						gmp2mp(fobj->nfs_obj.gmp_n, &fobj->N);	
						mpz_conv2str(&fobj->nfs_obj.mobj->input, 10, fobj->nfs_obj.gmp_n);

						do_data_check = 1;		// job file matches: check for data file
						if (VFLAG > 0) 
							printf("nfs: number in job file matches input\n");
					}
					else
					{
						if (VFLAG > 0)
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
		char master_polyfile[80];
		int do_poly_parse = 0;

		if (VFLAG > 0) 
			printf("nfs: checking for poly file - ");
		sprintf(master_polyfile,"%s.p",fobj->nfs_obj.outputfile);

		in = fopen(master_polyfile,"r");
		if (in == NULL)
		{
			if (VFLAG > 0) 
				printf("nfs: no poly file found\n");
			return NFS_STATE_STARTNEW;			// no .p file.
		}
		else
		{
			ptr = fgets(line,GSTR_MAXSIZE,in);
			if (ptr == NULL)
			{
				if (VFLAG > 0) 
					printf("nfs: poly file empty\n");
				return NFS_STATE_STARTNEW;		// .p file empty.
			}
			else
			{
				fclose(in);

				if (line[0] != 'n')
				{
					if (VFLAG > 0) 
						printf("nfs: malformed poly file, first line should contain n: \n");
					return NFS_STATE_STARTNEW;	// malformed .p file.
				}
				else
				{
					mpz_t tmp, g;
					mpz_init(tmp);
					mpz_init(g);
					mpz_set_str(tmp, line + 3, 0);
					if (resume_check_input_match(tmp, fobj->nfs_obj.gmp_n, g))
					{
						// divide out any common factor and copy the result to
						// other data structures
						if (mpz_cmp_ui(g, 1) > 0)
							add_to_factor_list(fobj, g);
						mpz_tdiv_q(fobj->nfs_obj.gmp_n, fobj->nfs_obj.gmp_n, g);
						gmp2mp(fobj->nfs_obj.gmp_n, &fobj->N);	
						mpz_conv2str(&fobj->nfs_obj.mobj->input, 10, fobj->nfs_obj.gmp_n);

						if (VFLAG > 0) 
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
				if (VFLAG > 0) printf("nfs: finding best poly in poly file\n");
				find_best_msieve_poly(fobj, job, 0);
				*last_spq = job->poly_time;
				if (VFLAG > 0) printf("nfs: last leading coefficient was %u\n", 
					job->last_leading_coeff);
				return NFS_STATE_RESUMEPOLY;
			}
			else
			{
				printf("nfs: must specify -R to resume when a polyfile already exists\n");	

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

				*last_spq = 0;
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
		
		//if (VFLAG > 0)
		//	printf("nfs: commencing NFS restart\n");

		logfile = fopen(fobj->flogname, "a");
		if (logfile == NULL)
		{
			printf("fopen error: %s\n", strerror(errno));
			printf("could not open yafu logfile for appending\n");
		}
		else
		{
			logprint(logfile, "nfs: commencing NFS restart\n");
			fclose(logfile);
		}

		printf("nfs: checking for data file\n");
		// attempt to open data file
		if (!savefile_exists(&mobj->savefile))
		{
			// data file doesn't exist, return flag to start sieving from the beginning
			if (VFLAG > 0)
				printf("nfs: no data file found\n");
			*last_spq = 0;
		}
		else if (fobj->nfs_obj.restart_flag)
		{
			// either restart from the end of the data file or from the specified 
			// sieve range
			if (fobj->nfs_obj.sieve_only && (fobj->nfs_obj.startq > 0))
			{
				if (VFLAG > 0)
					printf("nfs: user specified special-q range of %u-%u, "
					"skipping search for last special-q\n", 
					fobj->nfs_obj.startq, fobj->nfs_obj.rangeq);
				
				*last_spq = 0;
				return ans;
			}

			if (fobj->nfs_obj.la_restart || fobj->nfs_obj.post_only)
			{
				if (VFLAG > 0)
					printf("nfs: user specified post processing only, "
					"skipping search for last special-q\n");
				
				*last_spq = 0;
				return ans;
			}

			if (VFLAG > 0)
				printf("nfs: previous data file found - "
				"commencing search for last special-q\n");

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

			//tail isn't good enough, because prior filtering steps could have inserted
			//free relations, which don't have a special q to read.
			//our task is to find the last valid line of the file and parse it.
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
				//while (!feof(in))
				while (1)
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
			printf("must specify -R to resume when a savefile already exists\n");	

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

			*last_spq = 0;
			ans = NFS_STATE_DONE;
		}

	}
	else
	{
		// job file is for a different input.  not a restart.
		ans = NFS_STATE_STARTNEW;
		*last_spq = 0;
	}

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
		printf("nfs: parsing special-q from .dat file\n");

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

		if (VFLAG > 1)
			printf("parsing rat side spq from %s\n",lines[line]);
		sscanf(lines[line] + i + 1, "%x", &rat[count]);
		
		if (VFLAG > 1)
			printf("found %x\n", rat[count]);

		// grab alg side entry
		for (i= strlen(lines[line]) - 1; i >= 0; i--)
			if (lines[line][i] == ',')
				break;

		if (VFLAG > 1)
			printf("parsing alg side spq from %s\n",lines[line]);
		sscanf(lines[line] + i + 1, "%x", &alg[count]);
		if (VFLAG > 1)
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

void find_best_msieve_poly(fact_obj_t *fobj, ggnfs_job_t *job, int write_jobfile)
{
	// parse a msieve.dat.p file to find the best polynomial (based on e score)
	// output this as a ggnfs polynomial file
	FILE *in, *out, *logfile;
	char line[GSTR_MAXSIZE], *ptr;
	double score, bestscore = 0;
	int count, bestline = 0, i;
	uint32 highest_c4 = 0, highest_c5 = 0;

	sprintf(line, "%s.p",fobj->nfs_obj.outputfile);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading!\n",line);
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
			if (VFLAG > 2)
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
		else if ((line[0] == 'c') && (line[1] == '4'))
		{
			uint32 coeff;

			// get the c4 coefficient
			if (strlen(line) < 5)
				continue;	// not long enough to hold a line of format "c4: x"

			sscanf(line + 3, "%u", &coeff);
			if (coeff > highest_c4)
				highest_c4 = coeff;
		}
		else if ((line[0] == 'c') && (line[1] == '5'))
		{
			uint32 coeff;

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

	if (highest_c5 > 0)
		job->last_leading_coeff = highest_c5;
	else 
		job->last_leading_coeff = highest_c4;

	if (!write_jobfile)
		return;

	//open it again
	sprintf(line, "%s.p",fobj->nfs_obj.outputfile);
	in = fopen(line,"r");
	if (in == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for reading!\n",fobj->nfs_obj.outputfile);
		exit(1);
	}

	get_ggnfs_params(fobj, job);
	job->startq = job->fblim / 2;

	//always overwrites previous job files!
	out = fopen(fobj->nfs_obj.job_infile,"w");
	if (out == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open %s for writing!\n",fobj->nfs_obj.job_infile);
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
			{
				printf("fopen error: %s\n", strerror(errno));
				printf("could not open yafu logfile for appending\n");
			}
			else
			{
				logprint(logfile, "nfs: best poly = %s",line);
				fclose(logfile);
			}

			break;
		}
	}

	//copy n into the job file
	gmp_fprintf(out, "n: %Zd\n",fobj->nfs_obj.gmp_n);

	if (VFLAG > 0)
		gmp_printf("n: %Zd\n",fobj->nfs_obj.gmp_n);

	// copy out the poly
	while (!feof(in))
	{
		ptr = fgets(line,GSTR_MAXSIZE,in);
		if (ptr == NULL)
			break;

		if (line[0] == '#')
			break;		
		
		//prevent a "time" line from being copied into the .job file
		//ggnfs will ignore it, but it can be prevented here.
		if (strlen(line) > 4)
		{
			if ((line[0] == 't') && (line[1] == 'i') &&
				(line[2] == 'm') && (line[3] == 'e'))
				continue;
		}

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
	fprintf(out,"rlambda: %.1f\n",job->lambda);
	fprintf(out,"alambda: %.1f\n",job->lambda);

	fclose(in);
	fclose(out);

	return;
}

void msieve_to_ggnfs(fact_obj_t *fobj, ggnfs_job_t *job)
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
	fprintf(out,"rlim: %u\n",job->fblim);
	fprintf(out,"alim: %u\n",job->fblim);
	fprintf(out,"lpbr: %u\n",job->lpb);
	fprintf(out,"lpba: %u\n",job->lpb);
	fprintf(out,"mfbr: %u\n",job->mfb);
	fprintf(out,"mfba: %u\n",job->mfb);
	fprintf(out,"rlambda: %.1f\n",job->lambda);
	fprintf(out,"alambda: %.1f\n",job->lambda);

	fclose(in);
	fclose(out);

	return;
}

void ggnfs_to_msieve(fact_obj_t *fobj, ggnfs_job_t *job)
{
	// convert a ggnfs.job file into a msieve.fb polynomial file
	FILE *in, *out;
	char line[GSTR_MAXSIZE], outline[GSTR_MAXSIZE], *ptr;

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
		//else if(line[0] == 's')
		//	sprintf(outline, "SKEW %s",line + 5);
		else if (line[0] == 'Y')
			sprintf(outline, "R%c %s",line[1], line + 4);
		else if (line[0] == 'c')
			sprintf(outline, "A%c %s",line[1], line + 4);
		else if (line[0] == 'm' && line[1] == ':')
			sprintf(outline, "R1 1\nR0 -%s", line + 3);
		else
			strcpy(outline, "");

		fputs(outline,out);
	}

	fclose(in);
	fclose(out);

	return;
}

#endif
