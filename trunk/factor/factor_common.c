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
       				   --bbuhrow@gmail.com 3/26/10
----------------------------------------------------------------------*/

#include "yafu.h"
#include "soe.h"
#include "factor.h"
#include "qs.h"
#include "yafu_ecm.h"
#include "util.h"
#include "yafu_string.h"

int isUnique(z *in, z *list, int sz);
uint32 nCk(uint32 n, uint32 k);
void cycle_permute(z *list, int sz);

void init_factobj(fact_obj_t *fobj)
{
	// initialize global stuff in fobj
	fobj->seed1 = g_rand.low;
	fobj->seed2 = g_rand.hi;
	get_cache_sizes(&fobj->cache_size1,&fobj->cache_size2);
	fobj->flags = 0;
	strcpy(fobj->logname,flogname);
	fobj->num_threads = 1;		//read from input arguments
	zInit(&fobj->N);
	sInit(&fobj->str_N);

	// initialize stuff for rho
	fobj->rho_obj.num_factors = 0;
	fobj->rho_obj.factors = (z *)malloc(sizeof(z));
	fobj->rho_obj.num_poly = 3;
	fobj->rho_obj.polynomials = (uint32 *)malloc(fobj->rho_obj.num_poly * sizeof(uint32));
	zInit(&fobj->rho_obj.n);
	sInit(&fobj->rho_obj.in);
	sInit(&fobj->rho_obj.out);

	// initialize stuff for pm1
	fobj->pm1_obj.num_factors = 0;
	fobj->pm1_obj.factors = (z *)malloc(sizeof(z));
	zInit(&fobj->pm1_obj.n);
	sInit(&fobj->pm1_obj.in);
	sInit(&fobj->pm1_obj.out);

	// initialize stuff for pp1
	fobj->pp1_obj.num_factors = 0;
	fobj->pp1_obj.factors = (z *)malloc(sizeof(z));
	zInit(&fobj->pp1_obj.n);
	sInit(&fobj->pp1_obj.in);
	sInit(&fobj->pp1_obj.out);

	// initialize stuff for ecm
	fobj->ecm_obj.num_factors = 0;
	fobj->ecm_obj.factors = (z *)malloc(sizeof(z));
	zInit(&fobj->ecm_obj.n);
	sInit(&fobj->ecm_obj.in);
	sInit(&fobj->ecm_obj.out);

	// initialize stuff for squfof
	fobj->squfof_obj.num_factors = 0;
	fobj->squfof_obj.factors = (z *)malloc(sizeof(z));
	zInit(&fobj->squfof_obj.n);
	sInit(&fobj->squfof_obj.in);
	sInit(&fobj->squfof_obj.out);

	// initialize stuff for qs
	fobj->qs_obj.num_factors = 0;
	fobj->qs_obj.factors = (z *)malloc(sizeof(z));
	zInit(&fobj->qs_obj.n);
	sInit(&fobj->qs_obj.in);
	sInit(&fobj->qs_obj.out);

	// initialize stuff for trial division
	fobj->div_obj.num_factors = 0;
	fobj->div_obj.factors = (z *)malloc(sizeof(z));
	zInit(&fobj->div_obj.n);
	sInit(&fobj->div_obj.in);
	sInit(&fobj->div_obj.out);

	
	return;
}

void free_factobj(fact_obj_t *fobj)
{
	uint32 i;

	// free any factors found in rho
	for (i=0; i<fobj->rho_obj.num_factors; i++)
		zFree(&fobj->rho_obj.factors[i]);
	free(fobj->rho_obj.factors);

	// then free other stuff in rho
	free(fobj->rho_obj.polynomials);
	sFree(&fobj->rho_obj.in);
	sFree(&fobj->rho_obj.out);
	zFree(&fobj->rho_obj.n);

	// free any factors found in pm1
	for (i=0; i<fobj->pm1_obj.num_factors; i++)
		zFree(&fobj->pm1_obj.factors[i]);
	free(fobj->pm1_obj.factors);

	// then free other stuff in pm1
	sFree(&fobj->pm1_obj.in);
	sFree(&fobj->pm1_obj.out);
	zFree(&fobj->pm1_obj.n);

	// free any factors found in pp1
	for (i=0; i<fobj->pp1_obj.num_factors; i++)
		zFree(&fobj->pp1_obj.factors[i]);
	free(fobj->pp1_obj.factors);

	// then free other stuff in pp1
	sFree(&fobj->pp1_obj.in);
	sFree(&fobj->pp1_obj.out);
	zFree(&fobj->pp1_obj.n);

	// free any factors found in ecm
	for (i=0; i<fobj->ecm_obj.num_factors; i++)
		zFree(&fobj->ecm_obj.factors[i]);
	free(fobj->ecm_obj.factors);

	// then free other stuff in ecm
	sFree(&fobj->ecm_obj.in);
	sFree(&fobj->ecm_obj.out);
	zFree(&fobj->ecm_obj.n);

	// free any factors found in squfof
	for (i=0; i<fobj->squfof_obj.num_factors; i++)
		zFree(&fobj->squfof_obj.factors[i]);
	free(fobj->squfof_obj.factors);

	// then free other stuff in squfof
	sFree(&fobj->squfof_obj.in);
	sFree(&fobj->squfof_obj.out);
	zFree(&fobj->squfof_obj.n);

	// free any factors found in qs
	for (i=0; i<fobj->qs_obj.num_factors; i++)
		zFree(&fobj->qs_obj.factors[i]);
	free(fobj->qs_obj.factors);

	// then free other stuff in qs
	sFree(&fobj->qs_obj.in);
	sFree(&fobj->qs_obj.out);
	zFree(&fobj->qs_obj.n);

	// free any factors found in div
	for (i=0; i<fobj->div_obj.num_factors; i++)
		zFree(&fobj->div_obj.factors[i]);
	free(fobj->div_obj.factors);

	// then free other stuff in div
	sFree(&fobj->div_obj.in);
	sFree(&fobj->div_obj.out);
	zFree(&fobj->div_obj.n);

	//free general fobj stuff
	zFree(&fobj->N);
	sFree(&fobj->str_N);


	return;
}

void record_new_factor(fact_obj_t *fobj, char *method, z *n)
{
	//switch on 'method'
	if (strcmp(method,"rho"))
	{
		fobj->rho_obj.num_factors++;
		fobj->rho_obj.factors = (z *)realloc(fobj->rho_obj.factors,
			fobj->rho_obj.num_factors * sizeof(z));
		zInit(&fobj->rho_obj.factors[fobj->rho_obj.num_factors - 1]);
		zCopy(n,&fobj->rho_obj.factors[fobj->rho_obj.num_factors - 1]);
	}
	else if(strcmp(method,"pm1"))
	{
		fobj->pm1_obj.num_factors++;
		fobj->pm1_obj.factors = (z *)realloc(fobj->pm1_obj.factors,
			fobj->pm1_obj.num_factors * sizeof(z));
		zInit(&fobj->pm1_obj.factors[fobj->pm1_obj.num_factors - 1]);
		zCopy(n,&fobj->pm1_obj.factors[fobj->pm1_obj.num_factors - 1]);
	}
	else if(strcmp(method,"pp1"))
	{
		fobj->pp1_obj.num_factors++;
		fobj->pp1_obj.factors = (z *)realloc(fobj->pp1_obj.factors,
			fobj->pp1_obj.num_factors * sizeof(z));
		zInit(&fobj->pp1_obj.factors[fobj->pp1_obj.num_factors - 1]);
		zCopy(n,&fobj->pp1_obj.factors[fobj->pp1_obj.num_factors - 1]);
	}
	else if(strcmp(method,"ecm"))
	{
		fobj->ecm_obj.num_factors++;
		fobj->ecm_obj.factors = (z *)realloc(fobj->ecm_obj.factors,
			fobj->ecm_obj.num_factors * sizeof(z));
		zInit(&fobj->ecm_obj.factors[fobj->ecm_obj.num_factors - 1]);
		zCopy(n,&fobj->ecm_obj.factors[fobj->ecm_obj.num_factors - 1]);
	}
	else if(strcmp(method,"div"))
	{
		fobj->div_obj.num_factors++;
		fobj->div_obj.factors = (z *)realloc(fobj->div_obj.factors,
			fobj->div_obj.num_factors * sizeof(z));
		zInit(&fobj->div_obj.factors[fobj->div_obj.num_factors - 1]);
		zCopy(n,&fobj->div_obj.factors[fobj->div_obj.num_factors - 1]);
	}
	else if(strcmp(method,"squfof"))
	{
		fobj->squfof_obj.num_factors++;
		fobj->squfof_obj.factors = (z *)realloc(fobj->squfof_obj.factors,
			fobj->squfof_obj.num_factors * sizeof(z));
		zInit(&fobj->squfof_obj.factors[fobj->squfof_obj.num_factors - 1]);
		zCopy(n,&fobj->squfof_obj.factors[fobj->squfof_obj.num_factors - 1]);
	}
	else if(strcmp(method,"qs"))
	{
		fobj->qs_obj.num_factors++;
		fobj->qs_obj.factors = (z *)realloc(fobj->qs_obj.factors,
			fobj->qs_obj.num_factors * sizeof(z));
		zInit(&fobj->qs_obj.factors[fobj->qs_obj.num_factors - 1]);
		zCopy(n,&fobj->qs_obj.factors[fobj->qs_obj.num_factors - 1]);
	}
	else
	{
		printf("error: unknown method specified in record_new_factor\n");
		return;
	}

	return;
}


void save_state(int stage, z *n, z *m, z *mm, z *cc, int i, int method)
{
	FILE *outstate;
	uint32 limit;

	//save the result so far, so that future runs can resume using it.
	if (stage == 1)
	{
		if (method == 0)
		{
			outstate = fopen("pm1_stg1_save.dat","w");
			limit = POLLARD_STG1_MAX;
		}
		else
		{
			outstate = fopen("pp1_stg1_save.dat","w");
			limit = WILL_STG1_MAX;
		}

		//print n, the limit, and the residue
		fprintf(outstate,"%s\n",z2hexstr(n,&gstr1));
		fprintf(outstate,"%u\n",limit);
		fprintf(outstate,"%s\n",z2hexstr(m,&gstr1));
		fclose(outstate);
	}
	else
	{
		if (method == 0)
		{
			outstate = fopen("pm1_stg2_save.dat","w");
			limit = POLLARD_STG1_MAX;
		}
		else
		{
			outstate = fopen("pp1_stg2_save.dat","w");
			limit = WILL_STG1_MAX;
		}

		//print n, the limit, the resume point, and the residues
		fprintf(outstate,"%s\n",z2hexstr(n,&gstr1));
		fprintf(outstate,"%u,%u\n",limit,PRIMES[i]);
		fprintf(outstate,"%s\n",z2hexstr(cc,&gstr1));
		fprintf(outstate,"%s\n",z2hexstr(m,&gstr1));
		fprintf(outstate,"%s\n",z2hexstr(mm,&gstr1));

		fclose(outstate);
	}
	return;
}

void recover_stg1(z *n, z *m, int method)
{
	FILE *instate;
	int j;
	char *str;

	if (method == 0)
		instate = fopen("pm1_stg1_save.dat","r");
	else
		instate = fopen("pp1_stg1_save.dat","r");

	if (instate == NULL)
	{
		printf("no data file found\n");
		return;
	}

	str = (char *)malloc(1024*sizeof(char));

	//read in n
	fgets(str,1024,instate);
	str2hexz(str,n);

	fgets(str,1024,instate);
	j = atoi(str);

	//read in the residue
	fgets(str,1024,instate);
	str2hexz(str,m);

	fclose(instate);
	//make sure primes starting with where we left off are 
	//computed and ready to use
	GetPRIMESRange(j-1,j+10001000);
	free(str);
	return;
}

void recover_stg2(z *n, z *cc, z *m, z *mm, int method)
{
	FILE *instate;
	char c;
	int i,j;
	char *str;

	if (method == 0)
		instate = fopen("pm1_stg2_save.dat","r");
	else
		instate = fopen("pp1_stg2_save.dat","r");

	if (instate == NULL)
	{
		printf("no data file found\n");
		return;
	}

	str = (char *)malloc(1024*sizeof(char));

	//read in n
	fgets(str,1024,instate);
	str2hexz(str,n);

	//read in stage 1,2 max
	fgets(str,1024,instate);
	sscanf(str,"%u%c%u",&i,&c,&j);

	//read in the residues
	fgets(str,1024,instate);
	str2hexz(str,cc);

	fgets(str,1024,instate);
	str2hexz(str,m);

	fgets(str,1024,instate);
	str2hexz(str,mm);

	fclose(instate);

	//monty_init(n);
	GetPRIMESRange(j-1,j+10001000);
	free(str);
	return;
}

void add_to_factor_list(z *n)
{
	//stick the number n into the global factor list
	uint32 i;
	int found = 0;

	if (NUM_GLOBAL_FACTORS >= MAX_GLOBAL_FACTORS)
	{
		MAX_GLOBAL_FACTORS *= 2;
		global_factors = (factor_t *)realloc(global_factors,
			MAX_GLOBAL_FACTORS * sizeof(factor_t));
	}

	//look to see if this factor is already in the list
	for (i=0;i<NUM_GLOBAL_FACTORS && !found; i++)
	{
		if (zCompare(n,&global_factors[i].factor) == 0)
		{
			found = 1;
			global_factors[i].count++;
			return;
		}
	}

	//else, put it in the list
	zInit(&global_factors[NUM_GLOBAL_FACTORS].factor);
	zCopy(n,&global_factors[NUM_GLOBAL_FACTORS].factor);
	global_factors[NUM_GLOBAL_FACTORS].count = 1;
	NUM_GLOBAL_FACTORS++;

	return;
}

void free_factor_list(void)
{
	uint32 i;

	//clear this info
	for (i=0; i<NUM_GLOBAL_FACTORS; i++)
		zFree(&global_factors[i].factor);
	NUM_GLOBAL_FACTORS = 0;

	return;
}

void print_factors(void)
{
	uint32 i;
	int j;

	//always print factors unless complete silence is requested
	if (VFLAG >= 0)
	{
		printf("\n\n***factors found***\n\n");

		for (i=0; i<NUM_GLOBAL_FACTORS; i++)
		{
			if (global_factors[i].factor.type == COMPOSITE)
			{
				for (j=0;j<global_factors[i].count;j++)
					printf("C%d = %s\n",ndigits(&global_factors[i].factor),
						z2decstr(&global_factors[i].factor,&gstr1));
			}
			else if (global_factors[i].factor.type == PRP)
			{
				for (j=0;j<global_factors[i].count;j++)
					printf("PRP%d = %s\n",ndigits(&global_factors[i].factor),
						z2decstr(&global_factors[i].factor,&gstr1));
			}
			else if (global_factors[i].factor.type == PRIME)
			{
				for (j=0;j<global_factors[i].count;j++)
					printf("P%d = %s\n",ndigits(&global_factors[i].factor),
						z2decstr(&global_factors[i].factor,&gstr1));
			}
		}
	}

	return;
}

double get_qs_time_estimate(double freq, int digits)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think siqs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;

	cpu = get_cpu_type();
	switch (cpu)
	{
		case 0:
		case 1:
		case 2:
		case 3:
			if (VFLAG >= 2)
				printf("***** using 32bit windows p3 data for QS time estimation\n");
			//use p3 windows 32 bit @ 3.4 GHz estimate
			estimate = 0.0000967708 * exp(0.1981730566 * digits);
			//scale with frequency
			estimate = estimate * 3400.0 / freq; 

			break;
		case 4:
			
		case 5:
			if (VFLAG >= 2)
				printf("***** using 32bit windows p4 data for QS time estimation\n");
			//use p4 windows 32 bit @ 3.8 GHz estimate
			estimate = 0.0000836146 * exp(0.1984945754 * digits);
			//scale with frequency
			estimate = estimate * 3800.0 / freq; 
			 
			break;
		case 6:
			//subdivide into 32bit or 64bit architecture
			if (BITS_PER_DIGIT == 32)
			{
				if (VFLAG >= 2)
					printf("***** using 64bit windows core2 data for QS time estimation\n");
				estimate = 0.0000338875 * exp(0.2004412754 * digits);
			}
			else
			{
				if (VFLAG >= 2)
					printf("***** using 64bit linux core2 data for QS time estimation\n");
				estimate = 0.0000245199 * exp(0.2030407751 * digits);
			}

			//scale with frequency
			estimate = estimate * 3000.0 / freq; 

			break;
		case 7:		
		case 8:
		case 9:
			//subdivide into 32bit or 64bit architecture
			if (VFLAG >= 2)
				printf("***** using 64bit linux opteron data for QS time estimation\n");
			if (BITS_PER_DIGIT == 32)
				estimate = 0.0000213524 * exp(0.2066333088 * digits);
			else
				estimate = 0.0000213524 * exp(0.2066333088 * digits);

			//scale with frequency
			estimate = estimate * 2800.0 / freq; 

			break;
		case 10:
			//nehalem
			if (VFLAG >= 2)
				printf("***** using 64bit linux nehalem data for QS time estimation\n");

			estimate = 0.00002307333 * exp(0.2011087 * digits);

			//scale with frequency
			estimate = estimate * 2930.0 / freq; 

			break;

		default:
			if (VFLAG >= 2)
				printf("***** cpu type not found, ecm runtime will not be optimized\n");
			estimate =  0;
			break;
	}

	//adjust for multi-threaded qs
	//if we assume threading is perfect, we'll get a smaller estimate for
	//qs than we can really achieve, resulting in less ECM, so fudge it a bit
	if (THREADS > 1)
	{
		switch (cpu)
		{
		case 0:
		case 1:
		case 2:	
		case 3:
		case 4:
		case 5:			
		case 6:
		case 7:
		case 8:
			estimate = estimate / ((double)THREADS * 0.75);
			break;
		case 9:
		case 10:
			estimate = estimate / ((double)THREADS * 0.90);
			break;

		default:
			estimate = estimate / ((double)THREADS * 0.75);
			break;
		}
	}

	return estimate;
}

void factor(fact_obj_t *fobj)
{
	//run a varity of factoring algorithms on b
	//return any composite number left over
	//the factoring routines will build up a list of factors

	z *b = &fobj->N;
	int i, no_ecm_opt = 0;
	int curves = 1;
	int next = 0;
	int done = 0;
	uint32 auto_increasing_B1 = 10000000;
	uint64 auto_increasing_B2 = 1000000000;
	uint32 tmp1;
	uint64 tmp2;
	uint64 startticks, endticks;
	double freq, ecm_time = 0, qs_est_time = 0;
	double base_time_per_curve, time_per_curve;
	FILE *flog;
	struct timeval start, stop, start2, stop2;
	double t_time;
	TIME_DIFF *	difference;

	/*
	next = 0 -> williams
	next = 1 -> pollard
	next = 2 -> ecm 15 digit
	next = 3 -> ecm 20 digit
	next = 4 -> siqs
	next = 5 -> ecm 25 digit
	*/

	freq = MEAS_CPU_FREQUENCY;

	if (zCompare(b,&zOne) <= 0)
		return;

	gettimeofday(&start, NULL);

	flog = fopen(flogname,"a");
	logprint(flog,"\n");
	logprint(flog,"****************************\n");
	logprint(flog,"Starting factorization of %s\n",z2decstr(b,&gstr1));
	logprint(flog,"****************************\n");

	//if (VFLAG > 0)
	//	printf("***** cpu looks to be about %f MHz\n",freq);
	//logprint(flog,"cpu looks to be about %f MHz\n",freq);

	fclose(flog);

	AUTO_FACTOR=1;
	
	if (VFLAG >= 0)
	{
		printf("factoring %s\n\n",z2decstr(b,&gstr2));
		printf("div: primes less than %d\n",10000);
	}
	PRIME_THRESHOLD = 100000000;
	zCopy(b,&fobj->div_obj.n);
	zTrial(10000,0,fobj);
	zCopy(&fobj->div_obj.n,b);
	
	if (!(zCompare(b,&zOne) == 0))
	{
		if (VFLAG >= 0)
			printf("fmt: %d iterations\n",FMTMAX);
		zCopy(b,&fobj->N);
		zFermat(FMTMAX,fobj);
		zCopy(&fobj->N,b);
	}

	if (!(zCompare(b,&zOne) == 0))
	{
		zCopy(b,&fobj->rho_obj.n);
		brent_loop(fobj);
		zCopy(&fobj->rho_obj.n,b);

		if (!(b->size == 1 && b->val[0] == 1))
		{
			while (!done)
			{
				i = zBits(b);
				if (i < 100)
				{
					zCopy(b,&fobj->N);
					pQS(fobj);
					zCopy(&fobj->N,b);
					done = 1;
				}
				else if (i <= 135)
				{
					zCopy(b,&fobj->N);
					MPQS(fobj);
					zCopy(&fobj->N,b);
					done = 1;
				}
				else
				{
					switch(next)
					{
					case 0:
						tmp1 = WILL_STG1_MAX;
						tmp2 = WILL_STG2_MAX;
						WILL_STG1_MAX=20000;
						WILL_STG2_MAX=1000000;
						zCopy(b,&fobj->pp1_obj.n);
						williams_loop(3,fobj);
						zCopy(&fobj->pp1_obj.n,b);
						WILL_STG1_MAX=tmp1;
						WILL_STG2_MAX=tmp2;
						if (zBits(b) > 160)
							next++;
						else
							next = 4;
						break;
					case 1:
						tmp1 = POLLARD_STG1_MAX;
						tmp2 = POLLARD_STG2_MAX;
						POLLARD_STG1_MAX=100000;
						POLLARD_STG2_MAX=5000000;
						zCopy(b,&fobj->pm1_obj.n);
						pollard_loop(fobj);
						zCopy(&fobj->pm1_obj.n,b);
						POLLARD_STG1_MAX=tmp1;
						POLLARD_STG2_MAX=tmp2;
						if (zBits(b) > 180)
						{
							//the next step is ECM, so get an estimate of how
							//long the job would take to finish with siqs first
							next++;
							qs_est_time = get_qs_time_estimate(freq,ndigits(b));
							if (VFLAG >= 2)
								printf("***** qs time estimate = %f seconds\n",qs_est_time);

							//if qs_est_time is 0, set a flag to not do ECM
							//curve count optimization
							if (qs_est_time == 0)
								no_ecm_opt = 1;

						}
						else
							next = 4;
						break;
					case 2:
						//ecm to 15 digits
						//use these curves to estimate ECM speed on this
						//machine.  a better way would be to time each individual
						//curve and stop when time_elapsed > .25 * qs_est_time, but
						//right now don't have a good way to run single curves 
						//efficiently
						startticks = read_clock();
						tmp1 = ECM_STG1_MAX;
						tmp2 = ECM_STG2_MAX;
						ECM_STG1_MAX=2000;
						ECM_STG2_MAX=200000;
						curves = ecm_loop(b,25,fobj);
						ECM_STG1_MAX=tmp1;
						ECM_STG2_MAX=tmp2;
						endticks = read_clock();

						//estimate time per curve for this ecm level
						base_time_per_curve = (double)(endticks - startticks) 
							/ (freq * 1e6) / (double)curves;

						//record ecm time spent so far
						ecm_time = base_time_per_curve * curves;

						if (VFLAG >= 2)
							printf("***** 15 digit level took %f seconds.\n",
								base_time_per_curve * curves);

						//continue with next level of ecm.  decide how
						//many curves to do based on expected ECM time vs. 
						//expected QS time and how much time we've been 
						//ECMing so far
						
						//first scale the base time by the next level ECM
						//curve size
						time_per_curve = base_time_per_curve * 11000 / 2000;

						//estimate how many curves we can do in less than
						//0.25 * qs_est_time seconds have elapsed, taking into account
						//the time we've spent in ecm till now, as well.
						qs_est_time = get_qs_time_estimate(freq,ndigits(b));
						if (VFLAG >= 2)
							printf("***** qs time estimate = %f seconds\n",qs_est_time);

						curves = (int)((0.25 * qs_est_time - ecm_time)
							/ time_per_curve);

						//set the next number of curves at 11000 level
						next = 3;
						if (no_ecm_opt)
						{
							if (zBits(b) > 220)
								curves = 90;
							else
								next = 4;
						}
						else if (curves <= 0)
						{
							if (ndigits(b) < 130)
								next = 4;		//done ecm'ing, number ready for siqs
							else
								curves = 90;	//curves may have overflowed, set to max
						}
						else
						{
							if (curves > 90)
								curves = 90;	//cap to this digit limit
						}
						
						if (VFLAG >= 2)
							printf("***** estimating %d more curves can "
								"be run at 20 digit level\n",curves);

						break;
					case 3:
						//ecm to 20 digits
						startticks = read_clock();
						tmp1 = ECM_STG1_MAX;
						tmp2 = ECM_STG2_MAX;
						ECM_STG1_MAX=11000;
						ECM_STG2_MAX=1100000;
						curves = ecm_loop(b,curves,fobj);
						ECM_STG1_MAX=tmp1;
						ECM_STG2_MAX=tmp2;
						endticks = read_clock();

						//estimate time per curve for this ecm level
						base_time_per_curve = (double)(endticks - startticks) 
							/ (freq * 1e6) / (double)curves;

						//record ecm time spent so far
						ecm_time += (base_time_per_curve * curves);

						if (VFLAG >= 2)
							printf("***** 20 digit level took %f seconds.\n",
								base_time_per_curve * curves);

						
						//continue with next level of ecm.  decide how
						//many curves to do based on expected ECM time vs. 
						//expected QS time and how much time we've been 
						//ECMing so far
						
						//first scale the base time by the next level ECM
						//curve size
						time_per_curve = base_time_per_curve * 50000 / 11000;

						//estimate how many curves we can do in less than
						//0.25 * qs_est_time seconds have elapsed, taking into account
						//the time we've spent in ecm till now, as well.
						qs_est_time = get_qs_time_estimate(freq,ndigits(b));
						if (VFLAG >= 2)
							printf("***** qs time estimate = %f seconds\n",qs_est_time);

						curves = (int)((0.25 * qs_est_time - ecm_time)
							/ time_per_curve);

						//set the next number of curves at 50000 level
						next = 5;
						if (no_ecm_opt)
						{
							if (zBits(b) > 260)
								curves = 200;
							else
								next = 4;
						}
						else if (curves <= 0)
						{
							if (ndigits(b) < 130)
								next = 4;		//done ecm'ing, number ready for siqs
							else
								curves = 200;	//curves may have overflowed, set to max
						}
						else
						{
							if (curves > 200)
								curves = 200;	//cap to this digit limit
						}

						if (VFLAG >= 2)
							printf("***** estimating %d more curves can "
								"be run at 25 digit level\n",curves);

						break;
					case 4:
						gettimeofday(&start2,NULL);
						zCopy(b,&fobj->N);
						SIQS(fobj);
						zCopy(&fobj->N,b);
						gettimeofday(&stop2,NULL);
						difference = my_difftime (&start2, &stop2);
						t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
						free(difference);
						if (VFLAG > 0)
							printf("ECM/SIQS ratio was = %f\n",ecm_time/t_time);
						done = 1;
						break;
					case 5:
						//ecm to 25 digits
						startticks = read_clock();
						tmp1 = ECM_STG1_MAX;
						tmp2 = ECM_STG2_MAX;
						ECM_STG1_MAX=50000;
						ECM_STG2_MAX=5000000;
						curves = ecm_loop(b,curves,fobj);
						ECM_STG1_MAX=tmp1;
						ECM_STG2_MAX=tmp2;
						endticks = read_clock();

						//estimate time per curve for this ecm level
						base_time_per_curve = (double)(endticks - startticks) 
							/ (freq * 1e6) / (double)curves;

						//record ecm time spent so far
						ecm_time += (base_time_per_curve * curves);

						if (VFLAG >= 2)
							printf("***** 25 digit level took %f seconds.\n",
								base_time_per_curve * curves);

						//continue with next level of ecm.  decide how
						//many curves to do based on expected ECM time vs. 
						//expected QS time and how much time we've been 
						//ECMing so far
						
						//first scale the base time by the next level ECM
						//curve size
						time_per_curve = base_time_per_curve * 250000 / 50000;

						//estimate how many curves we can do in less than
						//0.25 * qs_est_time seconds have elapsed, taking into account
						//the time we've spent in ecm till now, as well.
						qs_est_time = get_qs_time_estimate(freq,ndigits(b));
						if (VFLAG >= 2)
							printf("***** qs time estimate = %f seconds\n",qs_est_time);

						curves = (int)((0.25 * qs_est_time - ecm_time)
							/ time_per_curve);

						//set the next number of curves at 250000 level
						next = 6;
						if (no_ecm_opt)
						{
							if (zBits(b) > 300)
								curves = 400;
							else
								next = 4;
						}
						else if (curves <= 0)
						{
							if (ndigits(b) < 130)
								next = 4;		//done ecm'ing, number ready for siqs
							else
								curves = 400;	//curves may have overflowed, set to max
						}
						else
						{
							if (curves > 400)
								curves = 400;	//cap to this digit limit
						}

						if (VFLAG >= 2)
							printf("***** estimating %d more curves can "
								"be run at 30 digit level\n",curves);

						break;
					case 6:
						//ecm to 30 digits
						startticks = read_clock();
						tmp1 = ECM_STG1_MAX;
						tmp2 = ECM_STG2_MAX;
						ECM_STG1_MAX=250000;
						ECM_STG2_MAX=25000000;
						curves = ecm_loop(b,curves,fobj);
						ECM_STG1_MAX=tmp1;
						ECM_STG2_MAX=tmp2;
						endticks = read_clock();

						//estimate time per curve for this ecm level
						base_time_per_curve = (double)(endticks - startticks) 
							/ (freq * 1e6) / (double)curves;

						//record ecm time spent so far
						ecm_time += (base_time_per_curve * curves);

						if (VFLAG >= 2)
							printf("***** 30 digit level took %f seconds.\n",
								base_time_per_curve * curves);

						//continue with next level of ecm.  decide how
						//many curves to do based on expected ECM time vs. 
						//expected QS time and how much time we've been 
						//ECMing so far
						
						//first scale the base time by the next level ECM
						//curve size
						time_per_curve = base_time_per_curve * 1000000 / 250000;

						//estimate how many curves we can do in less than
						//0.25 * qs_est_time seconds have elapsed, taking into account
						//the time we've spent in ecm till now, as well.
						qs_est_time = get_qs_time_estimate(freq,ndigits(b));
						if (VFLAG >= 2)
							printf("***** qs time estimate = %f seconds\n",qs_est_time);

						curves = (int)((0.25 * qs_est_time - ecm_time)
							/ time_per_curve);

						//set the next number of curves at 1000000 level
						next = 7;
						if (no_ecm_opt)
						{
							if (zBits(b) >= 330)
								curves = 1000;
							else
								next = 4;
						}
						else if (curves <= 0)
						{
							if (ndigits(b) < 130)
								next = 4;		//done ecm'ing, number ready for siqs
							else
								curves = 1000;	//curves may have overflowed, set to max
						}
						else
						{
							if (curves > 1000)
								curves = 1000;	//cap to this digit limit
						}

						if (VFLAG >= 2)
							printf("***** estimating %d more curves can "
								"be run at 35 digit level\n",curves);

						break;
					case 7:
						//ecm to 35 digits
						tmp1 = ECM_STG1_MAX;
						tmp2 = ECM_STG2_MAX;
						ECM_STG1_MAX=1000000;
						ECM_STG2_MAX=100000000;
						curves = ecm_loop(b,curves,fobj);
						ECM_STG1_MAX=tmp1;
						ECM_STG2_MAX=tmp2;

						if (ndigits(b) > 130)
						{
							if (VFLAG >= 2)
								printf("***** input too big for SIQS, continuing ECM\n");
							auto_increasing_B1 = 10000000;
							auto_increasing_B2 = 1000000000;
							curves = 3000;
							next = 8;
						}
						else
							next = 4;
						break;

					case 8:
						tmp1 = ECM_STG1_MAX;
						tmp2 = ECM_STG2_MAX;
						ECM_STG1_MAX=auto_increasing_B1;
						ECM_STG2_MAX=auto_increasing_B2;
						curves = ecm_loop(b,curves,fobj);
						ECM_STG1_MAX=tmp1;
						ECM_STG2_MAX=tmp2;

						if (ndigits(b) > 130)
						{
							if (VFLAG >= 2)
								printf("***** input too big for SIQS, continuing ECM\n");
							auto_increasing_B1 *= 10;
							auto_increasing_B2 *= 10;
							curves *= 3;
							next = 8;
						}
						else
							next = 4;
						
						break;
					}
				}
				if (isPrime(b))
				{
					b->type = PRP;
					add_to_factor_list(b);
					flog = fopen(flogname,"a");
					logprint(flog,"prp%d cofactor = %s\n",ndigits(b),z2decstr(b,&gstr1));
					fclose(flog);
					zCopy(&zOne,b);
				}

				if (isOne(b))
					done = 1;
			}
		}
	}
	if (!isOne(b))
	{
		b->type = COMPOSITE;
		flog = fopen(flogname,"a");
		logprint(flog,"c%d cofactor = %s\n",ndigits(b),z2decstr(b,&gstr1));
		fclose(flog);
		add_to_factor_list(b);
	}
	zCopy(b,&fobj->N);

	gettimeofday (&stop, NULL);
	difference = my_difftime (&start, &stop);
	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (VFLAG >= 0)
		printf("Total factoring time = %6.4f seconds\n",t_time);

	flog = fopen(flogname,"a");
	if (flog == NULL)
		printf("Could not open %s for appending\n",flogname);
	else
	{
		logprint(flog,"Total factoring time = %6.4f seconds\n",t_time);
		fclose(flog);
	}

	AUTO_FACTOR=0;
	return;
}

/*
Thanks to Jason Papadopoulos
*/

/* Implementation of the modified Knuth-Schroeppel multiplier
   algorithm. This borrows ideas from at least four different
   sources, and seems to choose multipliers that are better on
   average than many of the other methods available.
   
   There seem to be many misconceptions about what this algorithm
   is supposed to do. We want to multiply the input number n by a
   small squarefree constant k, chosen so that the factor base 
   for k * n contains as many small primes as possible. Since small primes
   occur more often than big ones, this makes sieve values smaller
   on average and so more likely to be smooth. We quantify this
   by measuring the average contribution of the first NUM_TEST_PRIMES
   primes to sieve values. There are two constraints: first, larger 
   multipliers mean a larger number to factor. Second, we can't spend 
   all day testing multipliers, so the set of multipliers to test should 
   be small. */

uint8 choose_multiplier(z *n, uint32 fb_size) 
{
	uint32 i, j;
	uint32 num_primes = MIN(2 * fb_size, NUM_TEST_PRIMES);
	double best_score;
	uint8 best_mult;
	double scores[NUM_MULTIPLIERS];
	uint32 num_multipliers;
	double log2n = zlog(n);

	/* measure the contribution of 2 as a factor of sieve
	   values. The multiplier itself must also be taken into
	   account in the score. scores[i] is the correction that
	   is implicitly applied to the size of sieve values for
	   multiplier i; a negative score makes sieve values 
	   smaller, and so is better */

	for (i = 0; i < NUM_MULTIPLIERS; i++) {
		uint8 curr_mult = mult_list[i];
		uint8 knmod8 = (uint8)((curr_mult * n->val[0]) % 8);
		double logmult = log((double)curr_mult);

		/* only consider multipliers k such than
		   k*n will not overflow an mp_t */

		if (log2n + logmult > (32 * MAX_DIGITS - 2) * LN2)
			break;

		scores[i] = 0.5 * logmult;
		switch (knmod8) {
		case 1:
			scores[i] -= 2 * LN2;
			break;
		case 5:
			scores[i] -= LN2;
			break;
		case 3:
		case 7:
			scores[i] -= 0.5 * LN2;
			break;
		/* even multipliers start with a handicap */
		}
	}
	num_multipliers = i;

	/* for the rest of the small factor base primes */
	for (i = 1; i < num_primes; i++) {
		uint32 prime = (uint32)spSOEprimes[i];
		double contrib = log((double)prime) / (prime - 1);
		uint32 modp = (uint32)zShortMod(n, prime);

		for (j = 0; j < num_multipliers; j++) {
			uint8 curr_mult = mult_list[j];
			//uint32 knmodp = mp_modmul_1(modp, curr_mult, prime);
			uint32 knmodp = (modp * curr_mult) % prime;

			/* if prime i is actually in the factor base
			   for k * n ... */

			//printf("prime %u, mult %u, knmodp = %u; trying jacobi_1\n",prime,curr_mult,knmodp);
			if (knmodp == 0 || jacobi_1(knmodp, prime) == 1) {

				/* ...add its contribution. A prime p con-
				   tributes log(p) to 1 in p sieve values, plus
				   log(p) to 1 in p^2 sieve values, etc. The
				   average contribution of all multiples of p 
				   to a random sieve value is thus

				   log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
				   = (log(p) / p) * 1 / (1 - (1/p)) 
				   = log(p) / (p-1)

				   This contribution occurs once for each
				   square root used for sieving. There are two
				   roots for each factor base prime, unless
				   the prime divides k*n. In that case there 
				   is only one root */

				//printf("scores[%d] = %f\n",j,contrib);
				if (knmodp == 0)
					scores[j] -= contrib;
				else
					scores[j] -= 2 * contrib;
			}
		}

	}

	/* use the multiplier that generates the best score */
	best_score = 1000.0;
	best_mult = 1;
	for (i = 0; i < num_multipliers; i++) {
		
		double score = scores[i];
		if (score < best_score) {
			best_score = score;
			best_mult = mult_list[i];
		}
	}
	return best_mult;
}

void aliquot(z *input, fact_obj_t *fobj)
{
	//compute the aliquot sequence of input
	//n_i+1 = sigma(n_i) - n_i
	//sum of proper divisors
	char bkup[80];
	FILE *log;
	z tmp, next, *proper, *factors, *seq, w1, w2, w3;
	uint32 i;
	int j,k,num_prop,num_seq;
	int proper_alloc,seq_alloc,num_fact,it;

	//use a different log file built from the input
	strcpy(bkup,flogname);
	sprintf(flogname,"aliquot_%u.log",(uint32)input->val[0]);
	log = fopen(flogname,"w");

	zInit(&tmp);
	zInit(&next);
	zInit(&w1);
	zInit(&w2);
	zInit(&w3);
	proper = (z *)malloc(8 * sizeof(z));
	proper_alloc = 8;
	factors = (z *)malloc(sizeof(z));
	seq = (z *)malloc(8 * sizeof(z));
	zInit(&seq[0]);
	zCopy(input,&seq[0]);
	seq_alloc = 8;
	num_seq = 1;
	num_prop = 0;
	VFLAG = 0;

	zCopy(input,&next);
	it = 0;
	while (!isOne(&next))
	{
		printf("iteration %d: %s\n",it++,z2decstr(&next,&gstr1));
		if (isPrime(&next))
		{
			printf("input %s is prime!\n",z2decstr(&next,&gstr1));
			break;
		}
		zCopy(&next,input);
		//get the divisors of the number
		zCopy(&next,&fobj->N);
		factor(fobj);
		zCopy(&fobj->N,&next);

		//put them in a simple list, and take care of any composite factors
		for (i=0; i<NUM_GLOBAL_FACTORS; i++)
		{
			if (global_factors[i].factor.type == COMPOSITE)
			{
				printf("composite factor detected, refactoring\n");
				zCopy(&global_factors[i].factor,&tmp);
				zCopy(&tmp,&fobj->N);
				factor(fobj);
				zCopy(&fobj->N,&tmp);
				if (zCompare(&tmp,&global_factors[i].factor) == 0)
				{
					printf("refactoring failed\n");
					exit(-1);
				}
			}
		}

		num_fact = 0;
		for (i=0; i<NUM_GLOBAL_FACTORS; i++)
		{
			if (global_factors[i].factor.type == COMPOSITE)
				continue;

			num_fact += global_factors[i].count;
		}

		zCopy(input,&w1);
		factors = (z *)realloc(factors,num_fact * sizeof(z));
		k=0;
		for (i=0; i<NUM_GLOBAL_FACTORS; i++)
		{
			if (global_factors[i].factor.type == COMPOSITE)
				continue;

			for (j=0;j<global_factors[i].count;j++)
			{
				zInit(&factors[k]);
				zCopy(&global_factors[i].factor,&factors[k]);
				//divide factor out of input
				zDiv(&w1,&factors[k],&w2,&w3);
				zCopy(&w2,&w1);
				//make sure it is a real divisor
				if (!isZero(&w3))
				{
					printf("improper divisor found\n");
					exit(-1);
				}
				k++;
			}
		}

		//make sure we found all the factors
		if (!isOne(&w1))
		{
			printf("factorization not completed\n");
			exit(-1);
		}

		num_prop = 0;
		for (i=0; i < num_fact; i++)
		{
			//multiply this factor by all others in the
			//proper divisor list, and add the result if unique
			k=num_prop;
			for (j=0; j<k; j++)
			{
				zMul(&factors[i],&proper[j],&tmp);
				if (isUnique(&tmp,proper,num_prop))
				{
					zInit(&proper[num_prop]);
					zCopy(&tmp,&proper[num_prop++]);
					if (num_prop >= proper_alloc)
					{
						proper_alloc *= 2;
						proper = realloc(proper,proper_alloc * sizeof(z));
					}
				}
			}

			//now add the factor itself, if unique
			if (isUnique(&factors[i],proper,num_prop))
			{
				zInit(&proper[num_prop]);
				zCopy(&factors[i],&proper[num_prop++]);
				if (num_prop >= proper_alloc)
				{
					proper_alloc *= 2;
					proper = realloc(proper,proper_alloc * sizeof(z));
				}
			}
		}
		zCopy(&zOne,&next);

		//sum them up
		for (i = 0; i<num_prop; i++)
			zAdd(&next,&proper[i],&next);

		zSub(&next,input,&next);

		for (i=0; i<num_fact; i++)
			zFree(&factors[i]);
		for (i=0; i<num_prop; i++)
			zFree(&proper[i]);

		free_factor_list();

		if (!isUnique(&next,seq,num_seq))
		{
			printf("cycle detected\n");
			break;
		}
		else
		{
			zInit(&seq[num_seq]);
			zCopy(&next,&seq[num_seq++]);
			if (num_seq >= seq_alloc)
			{
				seq_alloc *= 2;
				seq = realloc(seq,seq_alloc * sizeof(z));
			}
		}
	}

	printf("sequence terminated\n");

	VFLAG = 1;
	free(proper);
	free(factors);
	for (i=0; i<num_seq; i++)
		zFree(&seq[i]);
	free(seq);
	zFree(&tmp);
	zFree(&next);
	zFree(&w1);
	zFree(&w2);
	zFree(&w3);
	fclose(log);
	strcpy(flogname,bkup);
	return;
}

uint32 nCk(uint32 n, uint32 k)
{
	//nCk = n!/(n-k)!k!
	//works up to k = 20
	uint32 i;

	uint64 nn = n;
	uint64 kk;

	nn = 1;
	for (i = (n-k+1); i <= n; i++)
		nn *= (uint64)i;

	kk = 1;
	for (i = 2; i <= k; i++)
		kk *= (uint64)i;

	return (uint32)(nn/kk);
}


int isUnique(z *in, z *list, int sz)
{
	//return 1 if in is different from anything in list
	//return 0 otherwise
	int i;

	for (i=0; i<sz; i++)
	{
		if (zCompare(in,&list[i]) == 0)
			return 0;
	}

	return 1;
}

void cycle_permute(z *list, int sz)
{
	//replace list = {z1, z2, z3... zn}
	//with {zn, z1, z2... zn-1}
	z tmp;
	int i;

	zInit(&tmp);
	zCopy(&list[sz-1],&tmp);

	for (i=sz-2; i>=0; i--)
		zCopy(&list[i],&list[i+1]);

	zCopy(&tmp,&list[0]);

	zFree(&tmp);
	return;
}

