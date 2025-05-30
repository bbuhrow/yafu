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
#include "nfs_impl.h"
#ifdef __INTEL_LLVM_COMPILER
#include <pthread.h>
#endif
#include <math.h>

#ifdef USE_NFS

// potential enhancements:

// 2)
// run a series of test "refinements" based on difficulty.  The idea
// being twofold: a) for lower difficulty even a 2k range of spq might
// be overkill and b) for some polys the differences in sieve speed
// is obvious after only a few hundred spq.  So run a very small range
// first, and if we have a clear winner go with it.  if not, run more
// tests depending on both the "closeness" of the testing and the
// difficulty.
// 3)
// in windows, we count the number of ctrl-c's and force quit after 2.
// this defeats ctrl-c'ing out of more than one test sieve.  while test
// sieving, on windows, we need to allow more than two ctrl-c's


int qcomp_qrange(const void* x, const void* y)
{
	qrange_t* xx = (qrange_t*)x;
	qrange_t* yy = (qrange_t*)y;

	if ((*xx).qrange_start > (*yy).qrange_start)
		return 1;
	else if ((*xx).qrange_start == (*yy).qrange_start)
		return 0;
	else
		return -1;
}

void print_ranges(qrange_data_t* qrange_data)
{
	int i;
	printf("rational side completed ranges:\n");
	
	if (qrange_data->num_r > 0)
	{
		printf("\t%u -> ", qrange_data->qranges_r[0].qrange_start);
	}

	for (i = 1; i < qrange_data->num_r; i++)
	{
		if (qrange_data->qranges_r[i - 1].qrange_end !=
			qrange_data->qranges_r[i].qrange_start)
		{
			printf("%u\n\t%u -> \n", qrange_data->qranges_r[i - 1].qrange_end,
				qrange_data->qranges_r[i].qrange_start);
		}
	}

	if (qrange_data->num_r > 0)
	{
		printf("%u\n", qrange_data->qranges_r[i - 1].qrange_end);
	}

	printf("algebraic side completed ranges:\n");

	if (qrange_data->num_a > 0)
	{
		printf("\t%u -> ", qrange_data->qranges_a[0].qrange_start);
	}
	for (i = 1; i < qrange_data->num_a; i++)
	{
		if (qrange_data->qranges_a[i - 1].qrange_end !=
			qrange_data->qranges_a[i].qrange_start)
		{
			printf("%u\n\t%u -> \n", qrange_data->qranges_a[i - 1].qrange_end,
				qrange_data->qranges_a[i].qrange_start);
		}
	}

	if (qrange_data->num_a > 0)
	{
		printf("%u\n", qrange_data->qranges_a[i - 1].qrange_end);
	}
	return;
}

qrange_data_t* sort_completed_ranges(fact_obj_t* fobj, nfs_job_t *job)
{
	qrange_data_t* qrange_data;
	char buf[1024];
	FILE* fid;
	uint32_t numranges = 0;
	uint32_t totalrels = 0;

	qrange_data = (qrange_data_t*)xmalloc(sizeof(qrange_data_t));
	qrange_data->qranges_r = (qrange_t*)xmalloc(16 * sizeof(qrange_t));
	qrange_data->qranges_a = (qrange_t*)xmalloc(16 * sizeof(qrange_t));
	qrange_data->alloc_r = 16;
	qrange_data->alloc_a = 16;
	qrange_data->num_r = 0;
	qrange_data->num_a = 0;

	// make a sorted list of completed q-ranges
	sprintf(buf, "%s.ranges", fobj->nfs_obj.outputfile);
	fid = fopen(buf, "r");
	if (fid != NULL)
	{
		char side;
		uint32_t startq;
		uint32_t rangeq;
		uint32_t rels;
		int start = 1;

		if (fobj->VFLAG > 0)
		{
			printf("nfs: parsing %s.ranges for previously completed special-q\n", 
				fobj->nfs_obj.outputfile);
		}


		while (~feof(fid))
		{
			fgets(buf, 1024, fid);

			if (feof(fid))
			{
				break;
			}

			if (strlen(buf) < 10)
				continue;

			if (start)
			{
				mpz_t gmpn, gmpd;
				mpz_init(gmpn);
				mpz_init(gmpd);
				gmp_sscanf(buf, "%Zd", gmpn);
				
				mpz_tdiv_r(gmpd, gmpn, fobj->nfs_obj.gmp_n);
				if (mpz_cmp_ui(gmpd, 0) != 0)
				{
					printf("nfs: number in ranges file is not a divisor of the input\n");
					gmp_printf("nfs: read:  %Zd\n", gmpn);
					gmp_printf("nfs: input: %Zd\n", fobj->nfs_obj.gmp_n);
					printf("nfs: resetting ranges file with current input\n");
					fclose(fid);

					char newbuf[1024];
					sprintf(newbuf, "%s.ranges.bkup", fobj->nfs_obj.outputfile);
					rename(buf, newbuf);

					fid = fopen(buf, "w");
					gmp_fprintf(fid, "%Zd\n", fobj->nfs_obj.gmp_n);
					fclose(fid);

					fid = fopen(buf, "r");
					mpz_clear(gmpn);
					continue;
				}

				mpz_clear(gmpn);
				mpz_clear(gmpd);
				start = 0;
				continue;
			}
			else
			{
				sscanf(buf, "%c,%u,%u,%u", &side, &startq, &rangeq, &rels);
			}

			totalrels += rels;
			numranges++;
			//printf("parsed: %c,%u,%u,%u\n", side, startq, rangeq, rels);

			if (side == 'r')
			{
				if (qrange_data->num_r == qrange_data->alloc_r)
				{
					qrange_data->alloc_r *= 2;
					qrange_data->qranges_r = (qrange_t*)xrealloc(qrange_data->qranges_r,
						qrange_data->alloc_r * sizeof(qrange_t));
				}

				qrange_data->qranges_r[qrange_data->num_r].qrange_start = startq;
				qrange_data->qranges_r[qrange_data->num_r].qrange_end = startq + rangeq;
				qrange_data->num_r++;
			}
			else if (side == 'a')
			{
				if (qrange_data->num_a == qrange_data->alloc_a)
				{
					qrange_data->alloc_a *= 2;
					qrange_data->qranges_a = (qrange_t*)xrealloc(qrange_data->qranges_a,
						qrange_data->alloc_a * sizeof(qrange_t));
				}

				qrange_data->qranges_a[qrange_data->num_a].qrange_start = startq;
				qrange_data->qranges_a[qrange_data->num_a].qrange_end = startq + rangeq;
				qrange_data->num_a++;
			}
			else
			{
				printf("unrecognized side '%c' in completed range data\n", side);
			}

			
		}
		fclose(fid);
	}

	qsort(qrange_data->qranges_a, qrange_data->num_a, sizeof(qrange_t), &qcomp_qrange);
	qsort(qrange_data->qranges_r, qrange_data->num_r, sizeof(qrange_t), &qcomp_qrange);

	job->current_rels = totalrels;
	if (fobj->nfs_obj.rangeq > 0)
	{
		// user specified range: split threads over the entire range
		// and put bounds on new assignments.
		printf("nfs: configuring custom q-range %u-%u, splitting over %u threads\n",
			job->startq, job->startq + fobj->nfs_obj.rangeq, fobj->THREADS);
		qrange_data->thread_qrange =
			ceil((double)fobj->nfs_obj.rangeq / (double)fobj->THREADS);
		qrange_data->minq = job->startq;
		qrange_data->maxq = job->startq + fobj->nfs_obj.rangeq;
	}
	else
	{
		qrange_data->thread_qrange = job->qrange;
		qrange_data->maxq = 0xffffffff;

		if (fobj->nfs_obj.startq > 0)
			qrange_data->minq = fobj->nfs_obj.startq;
		else
			qrange_data->minq = 0;
	}

	if (fobj->VFLAG > 0)
	{
		print_ranges(qrange_data);
		printf("nfs: found %d previously completed ranges\n", numranges);
		printf("nfs: ranges file indicated %u previously found relations\n", totalrels);
	}

	return qrange_data;
}

void insert_range(qrange_data_t* qrange_data, char side, uint32_t start, uint32_t range)
{
	int i;
	if (side == 'a')
	{
		if (qrange_data->num_a == qrange_data->alloc_a)
		{
			qrange_data->alloc_a *= 2;
			qrange_data->qranges_a = (qrange_t*)xrealloc(qrange_data->qranges_a,
				qrange_data->alloc_a * sizeof(qrange_t));
		}

		qrange_data->qranges_a[qrange_data->num_a].qrange_start = start;
		qrange_data->qranges_a[qrange_data->num_a].qrange_end = start + range;
		qrange_data->num_a++;

		qsort(qrange_data->qranges_a, qrange_data->num_a, sizeof(qrange_t), &qcomp_qrange);
	}
	else
	{
		if (qrange_data->num_r == qrange_data->alloc_r)
		{
			qrange_data->alloc_r *= 2;
			qrange_data->qranges_r = (qrange_t*)xrealloc(qrange_data->qranges_r,
				qrange_data->alloc_r * sizeof(qrange_t));
		}

		qrange_data->qranges_r[qrange_data->num_r].qrange_start = start;
		qrange_data->qranges_r[qrange_data->num_r].qrange_end = start + range;
		qrange_data->num_r++;

		qsort(qrange_data->qranges_r, qrange_data->num_r, sizeof(qrange_t), &qcomp_qrange);
	}
	return;
}

qrange_t* get_next_range(qrange_data_t* qrange_data, char side)
{
	qrange_t* qrange = (qrange_t *)xmalloc(sizeof(qrange_t));
	int i;

	qrange->qrange_start = 0;
	qrange->qrange_end = 0;

	if (side == 'a')
	{
		for (i = 0; i < qrange_data->num_a - 1; i++)
		{
			if (qrange_data->qranges_a[i + 1].qrange_start >
				(qrange_data->qranges_a[i].qrange_end + 10))
			{
				// an unfinished range... work towards finishing it.
				qrange->qrange_start = qrange_data->qranges_a[i].qrange_end;
				qrange->qrange_end = qrange->qrange_start + qrange_data->thread_qrange;
				if (qrange->qrange_end > qrange_data->qranges_a[i + 1].qrange_start)
					qrange->qrange_end = qrange_data->qranges_a[i + 1].qrange_start;
				return qrange;
			}
		}

		if (qrange_data->num_a > 0)
		{
			qrange->qrange_start = qrange_data->qranges_a[i].qrange_end;
			qrange->qrange_end = qrange_data->qranges_a[i].qrange_end +
				qrange_data->thread_qrange;
		}
	}
	else
	{
		for (i = 0; i < qrange_data->num_r - 1; i++)
		{
			if (qrange_data->qranges_r[i + 1].qrange_start >
				(qrange_data->qranges_r[i].qrange_end + 10))
			{
				// an unfinished range... work towards finishing it.
				qrange->qrange_start = qrange_data->qranges_r[i].qrange_end;
				qrange->qrange_end = qrange->qrange_start + qrange_data->thread_qrange;
				if (qrange->qrange_end > qrange_data->qranges_r[i + 1].qrange_start)
					qrange->qrange_end = qrange_data->qranges_r[i + 1].qrange_start;
				return qrange;
			}
		}

		if (qrange_data->num_r > 0)
		{
			qrange->qrange_start = qrange_data->qranges_r[i].qrange_end;
			qrange->qrange_end = qrange_data->qranges_r[i].qrange_end +
				qrange_data->thread_qrange;
		}
	}

	if (qrange->qrange_start < qrange_data->minq)
	{
		qrange->qrange_start = qrange_data->minq;
		if (qrange->qrange_end < qrange->qrange_start)
			qrange->qrange_end = qrange->qrange_start;
	}

	if (qrange->qrange_end > qrange_data->maxq)
		qrange->qrange_end = qrange_data->maxq;

	return qrange;
}

int test_sieve(fact_obj_t* fobj, void* args, int njobs, int are_files)
/* if(are_files), then treat args as a char** list of (external) polys
 * else args is a nfs_job_t* array of job structs
 * the latter is preferred */
// @ben: the idea is a new yafu function "testsieve(n, ...)" where args are
// an arbitrary list of files of external polys
{
	uint32_t count;
	int i, minscore_id = 0;
	double* score = (double*)malloc(njobs * sizeof(double));
	double t_time, min_score = 999999999.;
	char orig_name[GSTR_MAXSIZE]; // don't clobber fobj->nfs_obj.job_infile
	char time[80];
	uint32_t spq_range = 1000, actual_range;
	FILE *flog = NULL;
    int sysreturn;

	char** filenames; // args
	nfs_job_t* jobs; // args
	
	struct timeval stop, stop2;	// stop time of this job
	struct timeval start, start2;	// start time of this job

	if( score == NULL )
	{
		printf("Couldn't alloc memory!\n");
		exit(-1);
	}

	// check to make sure we can find ggnfs sievers
	if (check_for_sievers(fobj, 0) == 1)
	{
		printf("test: can't find ggnfs lattice sievers - aborting test sieving\n");
		return -1;
	}

	gettimeofday(&start2, NULL);

	strcpy(orig_name, fobj->nfs_obj.job_infile);
	// necessary because parse/fill_job_file() get filename from fobj
	if( are_files )
	{ // read files into job structs (get fblim)

		filenames = (char**) args;
		jobs = (nfs_job_t*)malloc(njobs * sizeof(nfs_job_t));
		if( jobs == NULL )
		{
			printf("Couldn't alloc memory!\n");
			exit(-1);
		}
		memset(jobs, 0, njobs*sizeof(nfs_job_t));

		for(i = 0; i < njobs; i++)
		{
			uint32_t missing_params;
			strcpy(fobj->nfs_obj.job_infile, filenames[i]);

			missing_params = parse_job_file(fobj, jobs+i); // get fblim
			get_ggnfs_params(fobj, jobs+i); // get siever
			if( missing_params )
			{
				if( fobj->VFLAG >= 0 )
					printf("test: warning: \"%s\" is missing some paramters (%#X). filling them.\n",
						filenames[i], missing_params);
				fill_job_file(fobj, jobs+i, missing_params);
			}
			// adjust a/rlim, lpbr/a, and mfbr/a if advantageous
			skew_snfs_params(fobj, jobs+i);
		}
	}
	else
	{ // create poly files
		jobs = (nfs_job_t *) args;
		filenames = (char**)malloc(njobs*sizeof(char*));
		if( !filenames )
		{
			printf("malloc derped it up!\n");
			exit(-1);
		}
		for(i = 0; i < njobs; i++)
		{
			FILE* out;
			filenames[i] = (char*)malloc(GSTR_MAXSIZE);
			if( !filenames[i] )
			{
				printf("malloc failed\n");
				exit(-1);
			}
            sprintf(filenames[i], "test-sieve-%d.poly", i);
			out = fopen(filenames[i], "w");
			if( !out )
			{
				printf("test: couldn't open %s for writing, aborting test sieve\n", filenames[i]);
				return -1;
			}
			if( jobs[i].snfs )
				print_snfs(jobs[i].snfs, out);
			else
			{
				gmp_fprintf(out, "n: %Zd\n", fobj->nfs_obj.gmp_n);
				print_poly(jobs[i].poly, out);
			}
			fclose(out);
			strcpy(fobj->nfs_obj.job_infile, filenames[i]);
			fill_job_file(fobj, jobs+i, PARAM_FLAG_ALL);
		} 
		// that seems like a lot more code than it should be
	}
	strcpy(fobj->nfs_obj.job_infile, orig_name);

	// now we can get to the actual testing
	for(i = 0; i < njobs; i++)
	{
		char syscmd[GSTR_MAXSIZE], tmpbuf[GSTR_MAXSIZE], side[32];
		FILE* in;

        if (fobj->VFLAG > 0)
            printf("\ntest: trial sieving %s\n", filenames[i]);

		// should probably scale the range of special-q to test based
		// on input difficulty, but not sure how to do that easily...

		if( jobs[i].poly->side == RATIONAL_SPQ)
		{
			sprintf(side, "rational");
			jobs[i].startq = jobs[i].rlim / 2; // no reason to test sieve *inside* the fb
		}
		else
		{
			sprintf(side, "algebraic");
			//jobs[i].startq = jobs[i].alim; // ditto
		}

        if (fobj->LOGFLAG)
        {
            flog = fopen(fobj->flogname, "a");
            if (flog == NULL)
            {
                printf("could not open %s to append!\n", fobj->flogname);
                printf("disabling logging...\n");
                fobj->LOGFLAG = 0;
            }
        }

		//create the afb/rfb - we don't want the time it takes to do this to
		//pollute the sieve timings		
		sprintf(syscmd, "%s -b %s -k -c 0 -F", jobs[i].sievername, filenames[i]);
		if (fobj->VFLAG > 0) printf("\ntest: generating factor bases\n");
		gettimeofday(&start, NULL);
		sysreturn = system(syscmd);
		gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);

		if (fobj->VFLAG > 0) printf("test: fb generation took %6.4f seconds\n", t_time);
		logprint(flog, "test: fb generation took %6.4f seconds\n", t_time);
		MySleep(.1);

		//start the test
		sprintf(syscmd,"%s%s -%c %s -f %u -c %u -o %s.out",
			jobs[i].sievername, fobj->VFLAG>0?" -v":"", side[0], filenames[i], jobs[i].startq, spq_range, filenames[i]);

		if (fobj->VFLAG > 0) printf("test: commencing test sieving of polynomial %d on the %s side over range %u-%u\n", i, 
			side, jobs[i].startq, jobs[i].startq + spq_range);
		logprint(flog, "test: commencing test sieving of polynomial %d on the %s side over range %u-%u\n", i, 
			side, jobs[i].startq, jobs[i].startq + spq_range);
		
        if (fobj->LOGFLAG)
        {
            print_job(&jobs[i], flog);
            if (flog != NULL) fclose(flog);
        }

		gettimeofday(&start, NULL);
        sysreturn = system(syscmd);
		gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);
		
		//count relations
		sprintf(tmpbuf, "%s.out", filenames[i]);
		in = fopen(tmpbuf, "r");
		actual_range = 0;
		count = 0;
		if( !in )
		{
			score[i] = 999999999.;
			
			//est = 7*365*24*3600; // 7 years seems like a nice round number
		}
		else
		{
			// scan the data file and
			// 1) count the relations
			// 2) save the last four relations to a buffer in order to extract the last processed
			//		special-q.
			// we need both 1) and 2) to compute yield correctly.

			char **lines, *ptr, tmp[GSTR_MAXSIZE];
			int line;
			int j;

			lines = (char **)malloc(4 * sizeof(char *));
			for (j=0; j < 4; j++)
				lines[j] = (char *)malloc(GSTR_MAXSIZE * sizeof(char));

			line = 0;
			count = 0;
			while (1)
			{
				// read a line into the next position of the circular buffer
				ptr = fgets(tmp, GSTR_MAXSIZE, in);
				if (ptr == NULL) 
					break;

				// quick check that it might be a valid line
				if (strlen(tmp) > 30)
				{
					// wrap
					if (++line > 3) line = 0;
					// then copy
					strcpy(lines[line], tmp);
				}

				count++;
			}
			fclose(in);

			line = get_spq(lines, line, fobj);
			actual_range = line - jobs[i].startq;
			if (fobj->VFLAG > 0)
				printf("test: found %u relations in a range of %u special-q\n", 
				count, spq_range);

			if (actual_range > spq_range) 
				actual_range = spq_range;

			for (j=0; j < 4; j++)
				free(lines[j]);
			free(lines);

			score[i] = t_time / count;

			// use estimated sieving time to rank, not sec/rel, since the latter
			// is a function of parameterization and therefore not directly comparable
			// to each other.
			score[i] = (score[i] * jobs[i].min_rels * 1.25) / fobj->THREADS;
			// be conservative about estimates
		}

		jobs[i].test_score = score[i];

        if (fobj->LOGFLAG)
        {
            flog = fopen(fobj->flogname, "a");
        }

		if( score[i] < min_score )
		{
			minscore_id = i;
			min_score = score[i];
			if (fobj->VFLAG > 0) printf("test: new best estimated total sieving time = %s (with %d threads)\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);
			logprint(flog, "test: new best estimated total sieving time = %s (with %d threads)\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);

			if (0)
			{
				// don't edit by default.  Maybe allow with a user option.
				// edit lbpr/a depending on test results.  we target something around 2 rels/Q.
				// could also change siever version in more extreme cases.

				if (count > 4 * actual_range)
				{
					if (fobj->VFLAG > 0)
						printf("test: yield greater than 4x/spq, reducing lpbr/lpba\n");
					jobs[i].lpba--;
					jobs[i].lpbr--;
					jobs[i].mfba -= 2;
					jobs[i].mfbr -= 2;
				}

				if (count > 8 * actual_range)
				{
					char* pos;
					int siever;

					pos = strstr(jobs[i].sievername, "gnfs-lasieve4I");
					siever = (pos[14] - 48) * 10 + (pos[15] - 48);

					if (fobj->VFLAG > 0)
						printf("test: yield greater than 8x/spq, reducing siever version\n");

					switch (siever)
					{
					case 11:
						if (fobj->VFLAG > 0) printf("test: siever version cannot be decreased further\n");
						jobs[i].snfs->siever = 11;
						break;

					case 12:
						pos[15] = '1';
						jobs[i].snfs->siever = 11;
						break;

					case 13:
						pos[15] = '2';
						jobs[i].snfs->siever = 12;
						break;

					case 14:
						pos[15] = '3';
						jobs[i].snfs->siever = 13;
						break;

					case 15:
						pos[15] = '4';
						jobs[i].snfs->siever = 14;
						break;

					case 16:
						pos[15] = '5';
						jobs[i].snfs->siever = 15;
						break;
					}
				}

				if (count < actual_range)
				{
					if (fobj->VFLAG > 0)
						printf("test: yield less than 1x/spq, increasing lpbr/lpba\n");

					jobs[i].lpba++;
					jobs[i].lpbr++;
					jobs[i].mfba += 2;
					jobs[i].mfbr += 2;
				}

				if (count < (actual_range / 2))
				{
					char* pos;
					int siever;

					pos = strstr(jobs[i].sievername, "gnfs-lasieve4I");
					siever = (pos[14] - 48) * 10 + (pos[15] - 48);

					if (fobj->VFLAG > 0)
						printf("test: yield less than 1x/2*spq, increasing siever version\n");

					switch (siever)
					{
					case 16:
						if (fobj->VFLAG > 0) printf("test: siever version cannot be increased further\n");
						jobs[i].snfs->siever = 16;
						break;

					case 15:
						pos[15] = '6';
						jobs[i].snfs->siever = 16;
						break;

					case 14:
						pos[15] = '5';
						jobs[i].snfs->siever = 15;
						break;

					case 13:
						pos[15] = '4';
						jobs[i].snfs->siever = 14;
						break;

					case 12:
						pos[15] = '3';
						jobs[i].snfs->siever = 13;
						break;

					case 11:
						pos[15] = '2';
						jobs[i].snfs->siever = 12;
						break;
					}
				}

			}
     	}
		else
		{
			if (fobj->VFLAG > 0) printf("test: estimated total sieving time = %s (with %d threads)\n\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);
			logprint(flog, "test: estimated total sieving time = %s (with %d threads)\n", 
				time_from_secs(time, (unsigned long)score[i]), fobj->THREADS);
		}

        if (fobj->LOGFLAG)
        {
            if (flog != NULL) fclose(flog);
        }

		remove(tmpbuf); // clean up after ourselves

        if (!are_files)
        {
            sprintf(tmpbuf, "%s", filenames[i]);
            remove(tmpbuf);
        }

		sprintf(tmpbuf, "%s.afb.0", filenames[i]);
		remove(tmpbuf);
	}

	// clean up memory allocated
	if( are_files )
	{
        // TODO: need to write parameter adjustments to file
		for(i = 0; i < njobs; i++)
		{
			mpz_polys_free(jobs[i].poly);
			free(jobs[i].poly);
		}
		free(jobs);
	}
	else
	{
		for(i = 0; i < njobs; i++)
			free(filenames[i]);
		free(filenames);
	}

    if (fobj->LOGFLAG)
    {
        flog = fopen(fobj->flogname, "a");
        gettimeofday(&stop2, NULL);
        t_time = ytools_difftime(&start2, &stop2);

        if (fobj->VFLAG > 0) printf("test: test sieving took %1.2f seconds\n", t_time);
        logprint(flog, "test: test sieving took %1.2f seconds\n", t_time);
    }

	return minscore_id;
}

void do_sieving_nfs(fact_obj_t *fobj, nfs_job_t *job)
{
	nfs_threaddata_t *thread_data;
	int i;
	FILE *fid;
	FILE *logfile;
	char side;
	uint32_t totalrels;
	uint32_t custom_qstart = 0;
	uint32_t custom_qrange = 0;
	int make_afb = 0;

	side = (job->poly->side == ALGEBRAIC_SPQ) ? 'a' : 'r';

	thread_data = (nfs_threaddata_t *)malloc(fobj->THREADS * sizeof(nfs_threaddata_t));

	qrange_data_t* qrange_data = sort_completed_ranges(fobj, job);	

	if (fobj->nfs_obj.rangeq > 0)
	{
		custom_qstart = fobj->nfs_obj.startq;
		custom_qrange = qrange_data->thread_qrange;
	}

	for (i = 0; i < fobj->THREADS; i++)
	{
		sprintf(thread_data[i].outfilename, "rels%d.dat", i);
		sprintf(thread_data[i].job_infile_name, "%s", fobj->nfs_obj.job_infile);
		thread_data[i].job.poly = job->poly;
		thread_data[i].job.rlim = job->rlim;
		thread_data[i].job.alim = job->alim;
		thread_data[i].job.rlambda = job->rlambda;
		thread_data[i].job.alambda = job->alambda;
		thread_data[i].job.lpbr = job->lpbr;
		thread_data[i].job.lpba = job->lpba;
		thread_data[i].job.mfbr = job->mfbr;
		thread_data[i].job.mfba = job->mfba;
		thread_data[i].inflight = 0;

		qrange_t* qrange = get_next_range(qrange_data, side);

		if (custom_qstart > 0)
		{
			thread_data[i].job.startq = custom_qstart;
			thread_data[i].job.qrange = custom_qrange;
			custom_qstart += custom_qrange;
		}
		else
		{
			if ((qrange_data->num_a == 0) && (qrange_data->num_r == 0))
			{
				thread_data[i].job.startq = job->startq;
				thread_data[i].job.qrange = qrange_data->thread_qrange;
			}
			else
			{
				thread_data[i].job.startq = qrange->qrange_start;
				thread_data[i].job.qrange = qrange->qrange_end - qrange->qrange_start;
			}
		}

		// make this range unavailable to other threads
		insert_range(qrange_data, side, thread_data[i].job.startq, thread_data[i].job.qrange);
		free(qrange);

		thread_data[i].job.min_rels = job->min_rels;
		thread_data[i].job.current_rels = job->current_rels;
		thread_data[i].siever = fobj->nfs_obj.siever;
		
		strcpy(thread_data[i].job.sievername, job->sievername);

		thread_data[i].tindex = i;
		thread_data[i].is_poly_select = 0;
		thread_data[i].fobj = fobj;

		// have all data assigned, thread is now in-flight
		if (thread_data[i].job.qrange > 0)
		{
			if (i == (fobj->THREADS - 1))
			{
				nfs_start_worker_thread(thread_data + i, 1);
			}
			else
			{
				nfs_start_worker_thread(thread_data + i, 0);
			}
		}
		else
		{
			printf("could not get a valid assignment for sieving, "
				"thread %d is inactive\n", i);
		}
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
            logprint(logfile, "nfs: commencing lattice sieving with %d threads\n", fobj->THREADS);
            fclose(logfile);
        }
    }

	if (0)
	{
		// it would be faster to write the .afb once rather than
		// recompute it for every qrange that is run.  But if it
		// doesn't get deleted when the job completes that could
		// be problematic.  Need a way to query if existing .afb's
		// are valid for the current run.  Not sure how to do that yet.
		char syscmd[1024];
		
		sprintf(syscmd, "%s -b %s -k -c 0 -F", job->sievername, fobj->nfs_obj.job_infile);

		if (fobj->VFLAG > 1) printf("syscmd: %s\n", syscmd);
		if (fobj->VFLAG > 1) fflush(stdout);
		system(syscmd);
	}

	// create a new lasieve process in each thread and watch it
	for (i = 0; i < fobj->THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;

		if (!t->inflight)
			continue;

		if (i == fobj->THREADS - 1) {
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
	for (i = 0; i < fobj->THREADS; i++) {
		nfs_threaddata_t *t = thread_data + i;

		if (!t->inflight)
			continue;

		if (i < fobj->THREADS - 1) {
#if defined(WIN32) || defined(_WIN64)
			WaitForSingleObject(t->finish_event, INFINITE);
#else
			pthread_mutex_lock(&t->run_lock);
			while (t->command != NFS_COMMAND_WAIT)
				pthread_cond_wait(&t->run_cond, &t->run_lock);
#endif
		}
	}
		
	// combine output and log the range
	for (i = 0; i < fobj->THREADS; i++)
	{
		nfs_threaddata_t *t = thread_data + i;

		if (!t->inflight)
			continue;

		savefile_concat(t->outfilename,fobj->nfs_obj.outputfile,fobj->nfs_obj.mobj);

		uint32_t last_spq = t->job.qrange;
		if (NFS_ABORT > 0)
		{
			// try reading the lasieve5 ".last_spqX" file.  The 
			// following is copied from gnfs-lasieve4e to duplicate
			// the file naming convention.
			char* ofn = xmalloc(256);
			FILE* of;
			char *hn = xmalloc(128);
			int ret;

#if defined(WIN32)

			int sysname_sz = 128;
			GetComputerName((LPWSTR)hn, (LPDWORD)&sysname_sz);
			ret = 0;

#else

			ret = gethostname(hn, 127);

#endif

			if (ret == 0) sprintf(ofn, "%s.%s.last_spq%d", fobj->nfs_obj.job_infile, hn, i);
			else sprintf(ofn, "%s.unknown_host.last_spq%d", fobj->nfs_obj.job_infile, i);
			free(hn);

			if ((of = fopen(ofn, "r")) != 0) {
				fscanf(of, "%u", &last_spq);

				if (fobj->VFLAG > 0)
				{
					printf("nfs: parsed last_spq = %u in thread %d\n", last_spq, i);
				}
				if (last_spq < t->job.startq)
				{
					if (fobj->VFLAG > 0)
					{
						printf("nfs: last_spq is too small for range %u - %u, discarding.  Restarts may be incorrect.\n",
							t->job.startq, t->job.startq + t->job.qrange);
					}
					last_spq = t->job.qrange;
				}
				else if (last_spq > (t->job.startq + t->job.qrange))
				{
					if (fobj->VFLAG > 0)
					{
						printf("nfs: last_spq is too big for range %u - %u.  Assuming range was completed. Restarts may be incorrect.\n",
							t->job.startq, t->job.startq + t->job.qrange);
					}
					last_spq = t->job.qrange;
				}
				else
				{
					last_spq = last_spq - t->job.startq;
				}
				fclose(of);
			}
			else
			{
				if (fobj->VFLAG > 0)
				{
					printf("nfs: could not find file %s to parse last_spq for thread %d\n", ofn, i);
					printf("nfs: commencing search of data file\n");
				}
				if (1)
				{
					// file didn't exist.  Try parsing the special-q from
					// the data file, the lasieve4 way.
					char** lines = (char**)malloc(4 * sizeof(char*));
					char tmp[GSTR_MAXSIZE];
					int line = 0;

					for (line = 0; line < 4; line++)
						lines[line] = (char*)malloc(GSTR_MAXSIZE * sizeof(char));

					sprintf(ofn, "rels%d.dat", i);

					if ((of = fopen(ofn, "r")) != 0) {
						
						while (1)
						{
							// read a line into the next position of the circular buffer
							fgets(tmp, GSTR_MAXSIZE, of);	

							// quick check that it might be a valid line
							if (strlen(tmp) > 30)
							{
								// wrap
								if (++line > 3) line = 0;
								// then copy
								strcpy(lines[line], tmp);
							}

							if (feof(of))
								break;
						}
						fclose(of);

						// now we are done and we have a buffer with the last 4 valid lines
						// throw away the last one, which may be malformed, and extract
						// the special q from the other 3.
						for (line = 0; line < 4; line++)
						{
							if (fobj->VFLAG > 0)
							{
								printf("nfs: parsed line %d = %s", line, lines[line]);
							}
						}

						last_spq = get_spq(lines, line, fobj);

						if (fobj->VFLAG > 0)
						{
							printf("nfs: parsed last_spq = %u in thread %d\n", last_spq, i);
						}
						if (last_spq < t->job.startq)
						{
							if (fobj->VFLAG > 0)
							{
								printf("nfs: last_spq is too small for range %u - %u, discarding.  Restarts may be incorrect.\n",
									t->job.startq, t->job.startq + t->job.qrange);
							}
							last_spq = t->job.qrange;
						}
						else if (last_spq > (t->job.startq + t->job.qrange))
						{
							if (fobj->VFLAG > 0)
							{
								printf("nfs: last_spq is too big for range %u - %u.  Assuming range was completed. Restarts may be incorrect.\n",
									t->job.startq, t->job.startq + t->job.qrange);
							}
							last_spq = t->job.qrange;
						}
						else
						{
							last_spq = last_spq - t->job.startq;
						}
					}

					for (line = 0; line < 4; line++)
						free(lines[line]);
					free(lines);

				}
			}
			
			free(ofn);

		}

		char fname[1024];
		sprintf(fname, "%s.ranges", fobj->nfs_obj.outputfile);
		fid = fopen(fname, "a");
		if (fid != NULL)
		{
			fprintf(fid, "%c,%u,%u,%u\n", side, t->job.startq, 
				last_spq, t->job.current_rels);
			fclose(fid);
		}
		else
		{
			printf("nfs: could not open %s for modifying, "
				"progress may not be tracked correctly\n", fname);
		}

		// accumulate relation counts
		job->current_rels += thread_data[i].job.current_rels;

		remove(thread_data[i].outfilename);
	}

	if ((fid = fopen("rels.add", "r")) != NULL)
	{
		char tmpstr[1024];
		uint32_t count = 0;

		while (fgets(tmpstr, GSTR_MAXSIZE, fid) != NULL)
			count++;
		fclose(fid);

		if (fobj->VFLAG > 0) printf("nfs: adding %u rels from rels.add\n",count);

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
                logprint(logfile, "nfs: adding %u rels from rels.add\n", count);
                fclose(logfile);
            }
        }

		savefile_concat("rels.add",fobj->nfs_obj.outputfile,fobj->nfs_obj.mobj);
		remove("rels.add");

        job->current_rels += count;
	}

	//stop worker threads
	for (i = 0; i < fobj->THREADS - 1; i++)
	{
		if (!thread_data[i].inflight)
			continue;

		nfs_stop_worker_thread(thread_data + i, 0);
	}
	nfs_stop_worker_thread(thread_data + i, 1);

	//free the thread structure
	free(thread_data);
	free(qrange_data->qranges_a);
	free(qrange_data->qranges_r);
	free(qrange_data);

	return;
}


void *lasieve_launcher(void *ptr)
{
	// launch a gnfs-lasieve job
	nfs_threaddata_t *thread_data = (nfs_threaddata_t *)ptr;
	fact_obj_t *fobj = thread_data->fobj;
	char syscmd[GSTR_MAXSIZE], tmpstr[GSTR_MAXSIZE], side[GSTR_MAXSIZE];	
	FILE *fid;
	int cmdret;

	sprintf(side, (thread_data->job.poly->side == ALGEBRAIC_SPQ) ? 
				"algebraic" : "rational"); // gotta love ?:

	//remove any temporary relation files
	remove(thread_data->outfilename);
		
	//start ggnfs binary - new win64 ASM enabled binaries current have a problem with this:
	//sprintf(syscmd,"%s%s -%c %s -f %u -c %u -o %s -n %d",
	//		thread_data->job.sievername, fobj->VFLAG>0?" -v":"", *side,
	//		fobj->nfs_obj.job_infile, thread_data->job.startq, 
	//		thread_data->job.qrange, thread_data->outfilename, thread_data->tindex);

    // todo: add command line input of arbitrary argument string to append to this command
	// but not this:
	snprintf(syscmd, GSTR_MAXSIZE, "%s%s -f %u -c %u -o %s -n %d -%c %s ",
			thread_data->job.sievername, fobj->VFLAG>0?" -v":"", thread_data->job.startq, 
			thread_data->job.qrange, thread_data->outfilename, thread_data->tindex,
			*side, thread_data->job_infile_name);

	if (fobj->VFLAG >= 0)
	{
		printf("nfs: commencing %s side lattice sieving over range: %u - %u\n",
			side, thread_data->job.startq, thread_data->job.startq + thread_data->job.qrange);
	}
	if (fobj->VFLAG > 1) printf("syscmd: %s\n", syscmd);
	if (fobj->VFLAG > 1) fflush(stdout);
	cmdret = system(syscmd);

	// a ctrl-c abort signal is caught by the system command, and nfsexit never gets called.
	// so check for abnormal exit from the system command.
	// -1073741819 is apparently what ggnfs returns when it crashes, which
	// we don't want to interpret as an abort.
    if (!((cmdret == 0) || cmdret == -1073741819))
    {
        if (NFS_ABORT < 1)
        {
            NFS_ABORT = 1;
        }
    }

	// count the relations produced
	MySleep(100);
	fid = fopen(thread_data->outfilename,"r");
	if (fid != NULL)
	{
		thread_data->job.current_rels = 0;
		while (fgets(tmpstr, GSTR_MAXSIZE, fid) != NULL)
			thread_data->job.current_rels++;
		fclose(fid);
	}
	else
	{
		printf("nfs: could not open output file %s, possibly bad path to siever\n",
			thread_data->outfilename);
	}

	return 0;
}

#endif
