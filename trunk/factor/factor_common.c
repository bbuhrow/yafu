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

// what would be neat: a file format definition and library of functions 
// for recording the amount of work done on a number.  the library would
// be able to parse the file to determine the next optimal factorization step
// based on the recorded work history.  the library would also be able to
// merge work files and display them nicely.  this would allow people to 
// pass the work they have done to another person and have that work 
// immediately taken into account by the library to determine the next optimal 
// factorization step.  a centralized repository (factordb) could then merge
// work files to facilitate large scale contributions for difficult numbers.
typedef struct
{
	double trialdiv_time;
	double fermat_time;
	double rho_time;
	double pp1_lvl1_time_per_curve;
	double pp1_lvl2_time_per_curve;
	double pp1_lvl3_time_per_curve;
	double pm1_lvl1_time_per_curve;
	double pm1_lvl2_time_per_curve;
	double pm1_lvl3_time_per_curve;
	double ecm_15digit_time_per_curve;
	double ecm_20digit_time_per_curve;
	double ecm_25digit_time_per_curve;
	double ecm_30digit_time_per_curve;
	double ecm_35digit_time_per_curve;
	double ecm_autoinc_time_per_curve;
} method_timing_t;

enum work_method {
	trialdiv_work,
	fermat_work,
	rho_work,
	ecm_curve,
	pp1_curve,
	pm1_curve,
	qs_work,
	nfs_work
};

enum factorization_state {
	state_idle,
	state_trialdiv,
	state_fermat,
	state_rho,
	state_pp1_lvl1,
	state_pp1_lvl2,
	state_pp1_lvl3,
	state_pm1_lvl1,
	state_pm1_lvl2,
	state_pm1_lvl3,
	state_ecm_15digit,
	state_ecm_20digit,
	state_ecm_25digit,
	state_ecm_30digit,
	state_ecm_35digit,
	state_ecm_auto_increasing,
	state_qs,
	state_nfs
};

static int refactor_depth = 0;
double total_time = 0;
uint32 auto_increasing_curves = 3000;
uint32 auto_increasing_B1 = 10000000;
uint64 auto_increasing_B2 = 1000000000;

int isUnique(z *in, z *list, int sz);
uint32 nCk(uint32 n, uint32 k);
void cycle_permute(z *list, int sz);

// local function to do requested curve based factorization
double do_work(enum work_method method, uint32 B1, uint64 B2, int *work, 
	z *b, fact_obj_t *fobj);
int check_if_done(fact_obj_t *fobj, z *N);
enum factorization_state scale_requested_work(method_timing_t *method_times, 
	enum factorization_state fact_state, int *next_work, double time_available, z *N);
int switch_to_qs(z *N, double *time_available);

void init_factobj(fact_obj_t *fobj)
{
	// initialize global stuff in fobj
	fobj->seed1 = g_rand.low;
	fobj->seed2 = g_rand.hi;
	yafu_get_cache_sizes(&fobj->cache_size1,&fobj->cache_size2);
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

	//global list of factors
	fobj->allocated_factors = 256;
	fobj->fobj_factors = (factor_t *)malloc(256 * sizeof(factor_t));
	fobj->num_factors = 0;

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

	clear_factor_list(fobj);
	free(fobj->fobj_factors);

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

void add_to_factor_list(fact_obj_t *fobj, z *n)
{
	//stick the number n into the global factor list
	uint32 i;
	int found = 0;

	if (fobj->num_factors >= fobj->allocated_factors)
	{
		fobj->allocated_factors *= 2;
		fobj->fobj_factors = (factor_t *)realloc(fobj->fobj_factors,
			fobj->allocated_factors * sizeof(factor_t));
	}

	//look to see if this factor is already in the list
	for (i=0;i<fobj->num_factors && !found; i++)
	{
		if (zCompare(n,&fobj->fobj_factors[i].factor) == 0)
		{
			found = 1;
			fobj->fobj_factors[i].count++;
			return;
		}
	}

	//else, put it in the list
	zInit(&fobj->fobj_factors[fobj->num_factors].factor);
	zCopy(n,&fobj->fobj_factors[fobj->num_factors].factor);
	fobj->fobj_factors[fobj->num_factors].count = 1;
	fobj->num_factors++;

	return;
}

void delete_from_factor_list(fact_obj_t *fobj, z *n)
{
	//remove the number n from the global factor list
	uint32 i;

	//find the factor
	for (i=0;i<fobj->num_factors; i++)
	{
		if (zCompare(n,&fobj->fobj_factors[i].factor) == 0)
		{
			int j;
			// copy everything above this in the list back one position
			for (j=i; j<fobj->num_factors-1; j++)
			{
				zCopy(&fobj->fobj_factors[j+1].factor, &fobj->fobj_factors[j].factor);
				fobj->fobj_factors[j].count = fobj->fobj_factors[j+1].count;
			}
			// remove the last one in the list
			fobj->fobj_factors[j].count = 0;
			zFree(&fobj->fobj_factors[j].factor);

			fobj->num_factors--;
			break;
		}
	}

	return;
}

void clear_factor_list(fact_obj_t *fobj)
{
	uint32 i;

	//clear this info
	for (i=0; i<fobj->num_factors; i++)
	{
		fobj->fobj_factors[i].count = 0;
		zFree(&fobj->fobj_factors[i].factor);
	}
	fobj->num_factors = 0;

	return;
}

void print_factors(fact_obj_t *fobj)
{
	uint32 i;
	int j;
	z tmp, tmp2, tmp3, tmp4;

	//always print factors unless complete silence is requested
	if (VFLAG >= 0)
	{
		printf("\n\n***factors found***\n\n");
		zInit(&tmp);
		zCopy(&zOne,&tmp);

		for (i=0; i<fobj->num_factors; i++)
		{
			if (fobj->fobj_factors[i].factor.type == COMPOSITE)
			{
				for (j=0;j<fobj->fobj_factors[i].count;j++)
				{
					zMul(&tmp,&fobj->fobj_factors[i].factor,&tmp);
					printf("C%d = %s\n",ndigits(&fobj->fobj_factors[i].factor),
						z2decstr(&fobj->fobj_factors[i].factor,&gstr1));
				}
			}
			else if (fobj->fobj_factors[i].factor.type == PRP)
			{
				for (j=0;j<fobj->fobj_factors[i].count;j++)
				{
					printf("PRP%d = %s\n",ndigits(&fobj->fobj_factors[i].factor),
						z2decstr(&fobj->fobj_factors[i].factor,&gstr1));
					zMul(&tmp,&fobj->fobj_factors[i].factor,&tmp);
				}
			}
			else if (fobj->fobj_factors[i].factor.type == PRIME)
			{
				for (j=0;j<fobj->fobj_factors[i].count;j++)
				{
					printf("P%d = %s\n",ndigits(&fobj->fobj_factors[i].factor),
						z2decstr(&fobj->fobj_factors[i].factor,&gstr1));
					zMul(&tmp,&fobj->fobj_factors[i].factor,&tmp);
				}
			}
		}

		if (zCompare(&fobj->N, &zOne) > 0)
		{
			// non-trivial N remaining... compute and display the known co-factor
			zInit(&tmp2);
			zInit(&tmp3);
			zInit(&tmp4);
			zCopy(&fobj->N, &tmp2);
			zDiv(&tmp2, &tmp, &tmp3, &tmp4);
			if (zCompare(&tmp3,&zOne) != 0)
			{
//				printf("original N = %s\naccumulated value = %s\n",
	//				z2decstr(&fobj->N, &gstr1),z2decstr(&tmp, &gstr2));
				if (isPrime(&tmp3))
				{
					printf("\n***co-factor***\nPRP%d = %s\n",
						ndigits(&tmp3), z2decstr(&tmp3, &gstr1));
				}
				else
				{
					printf("\n***co-factor***\nC%d = %s\n",
						ndigits(&tmp3), z2decstr(&tmp3, &gstr1));
				}
			}			
			zFree(&tmp2);
			zFree(&tmp3);
			zFree(&tmp4);
		}
		zFree(&tmp);

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

	cpu = yafu_get_cpu_type();
	//if we have tuning info, use that instead
	if (QS_MULTIPLIER != 0 && QS_EXPONENT != 0 && QS_TUNE_FREQ != 0)
	{
		if (VFLAG >= 2)
			printf("***** using tuning data for QS time estimation\n");
		estimate = QS_MULTIPLIER * exp(QS_EXPONENT * digits);

		//scale with frequency
		estimate = estimate * QS_TUNE_FREQ / freq; 
	}
	else
	{		
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
				estimate = DBL_MAX;
				break;
		}
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

double get_gnfs_time_estimate(double freq, int digits)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think gnfs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;

	cpu = yafu_get_cpu_type();
	//if we have tuning info, use that instead
	if (GNFS_MULTIPLIER != 0 && GNFS_EXPONENT != 0 && GNFS_TUNE_FREQ != 0)
	{
		if (VFLAG >= 2)
			printf("***** using tuning data for GNFS time estimation\n");
		estimate = GNFS_MULTIPLIER * exp(GNFS_EXPONENT * digits);

		//scale with frequency
		estimate = estimate * GNFS_TUNE_FREQ / freq; 
	}
	else
		estimate = DBL_MAX;

	//adjust for multi-threaded nfs
	//if we assume threading is perfect, we'll get a smaller estimate for
	//nfs than we can really achieve, resulting in less ECM, so fudge it a bit
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

double do_work(enum work_method method, uint32 B1, uint64 B2, int *work, 
	z *b, fact_obj_t *fobj)
{
	uint32 tmp1;
	uint64 tmp2;
	uint64 startticks, endticks;
	double time_per_unit_work;	
	//struct timeval start2, stop2;
	//double t_time;
	//TIME_DIFF *	difference;

	startticks = yafu_read_clock();	

	switch (method)
	{
	case trialdiv_work:		
		if (VFLAG >= 0)
			printf("div: primes less than %d\n",10000);
		PRIME_THRESHOLD = 100000000;
		zCopy(b,&fobj->div_obj.n);
		zTrial(B1,0,fobj);
		zCopy(&fobj->div_obj.n,b);		
		break;

	case rho_work:
		zCopy(b,&fobj->rho_obj.n);
		brent_loop(fobj);
		zCopy(&fobj->rho_obj.n,b);
		break;

	case fermat_work:
		if (VFLAG >= 0)
			printf("fmt: %d iterations\n",B1);
		zCopy(b,&fobj->div_obj.n);
		zFermat(B1,fobj);
		zCopy(&fobj->div_obj.n,b);		
		break;

	case ecm_curve:
		tmp1 = ECM_STG1_MAX;
		tmp2 = ECM_STG2_MAX;
		ECM_STG1_MAX=B1;
		ECM_STG2_MAX=B2;
		*work = ecm_loop(b,*work,fobj);
		ECM_STG1_MAX=tmp1;
		ECM_STG2_MAX=tmp2;
		break;

	case pp1_curve:
		tmp1 = WILL_STG1_MAX;
		tmp2 = WILL_STG2_MAX;
		WILL_STG1_MAX=B1;
		WILL_STG2_MAX=B2;
		zCopy(b,&fobj->pp1_obj.n);
		williams_loop(*work,fobj);
		zCopy(&fobj->pp1_obj.n,b);
		WILL_STG1_MAX=tmp1;
		WILL_STG2_MAX=tmp2;
		break;

	case pm1_curve:
		tmp1 = POLLARD_STG1_MAX;
		tmp2 = POLLARD_STG2_MAX;
		POLLARD_STG1_MAX=B1;
		POLLARD_STG2_MAX=B2;
		zCopy(b,&fobj->pm1_obj.n);
		pollard_loop(fobj);
		zCopy(&fobj->pm1_obj.n,b);
		POLLARD_STG1_MAX=tmp1;
		POLLARD_STG2_MAX=tmp2;
		break;

	case qs_work:
		//gettimeofday(&start2,NULL);
		zCopy(b,&fobj->qs_obj.n);
		SIQS(fobj);
		zCopy(&fobj->qs_obj.n,b);
		break;
		//gettimeofday(&stop2,NULL);
		//difference = my_difftime (&start2, &stop2);
		//t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		//free(difference);

	case nfs_work:
		zCopy(b,&fobj->qs_obj.n);
		test_msieve_gnfs(fobj);
		zCopy(&fobj->qs_obj.n,b);
		break;

	default:
		printf("nothing to do for method %d\n",method);
		break;
	}

	endticks = yafu_read_clock();

	//estimate time per curve for these completed curves
	time_per_unit_work = (double)(endticks - startticks) 
		/ (MEAS_CPU_FREQUENCY * 1e6) / (double)(*work);

	return time_per_unit_work;
}

int check_if_done(fact_obj_t *fobj, z *N)
{
	int i, done = 0;
	z tmp;

	zInit(&tmp);
	zCopy(&zOne,&tmp);

	//printf("checking if prod of factors = %s\n",z2decstr(N,&gstr1));

	// check if the number is completely factorized
	for (i=0; i<fobj->num_factors; i++)
	{		
		int j;
		for (j=0; j<fobj->fobj_factors[i].count; j++)
		{
			zMul(&tmp,&fobj->fobj_factors[i].factor,&tmp);
			//printf("accumulating factor %s = %s\n",z2decstr(&global_factors[i].factor,&gstr1), 
				//z2decstr(&tmp,&gstr2));
		}		
	}

	if (zCompare(N,&tmp) == 0)
	{		
		// yes, they are equal.  make sure everything is prp or prime.
		done = 1;
		for (i=0; i<fobj->num_factors; i++)
		{
			if (fobj->fobj_factors[i].factor.type == COMPOSITE)
			{
				refactor_depth++;
				if (refactor_depth > 3)
				{
					printf("too many refactorization attempts, aborting\n");
				}
				else
				{
					//printf("ignoring refactorization of composite factor\n");
					fact_obj_t *fobj_refactor;
					int j;

					printf("\nComposite result found, starting re-factorization\n");

					// load the new fobj with this number
					fobj_refactor = (fact_obj_t *)malloc(sizeof(fact_obj_t));
					init_factobj(fobj_refactor);
					zCopy(&fobj->fobj_factors[i].factor,&fobj_refactor->N);

					// recurse on factor
					factor(fobj_refactor);

					// remove the factor from the original list
					delete_from_factor_list(fobj, &fobj->fobj_factors[i].factor);

					// add all factors found during the refactorization
					for (j=0; j< fobj_refactor->num_factors; j++)
					{
						int k;
						for (k=0; k < fobj_refactor->fobj_factors[j].count; k++)
							add_to_factor_list(fobj, &fobj_refactor->fobj_factors[j].factor);
					}

					// factorization completed if we got to here.  reset the recursion limit.
					refactor_depth = 0;

					// free temps
					free_factobj(fobj_refactor);
					free(fobj_refactor);
				}
			}
		}
	}

	zFree(&tmp);
	return done;
}

int switch_to_qs(z *N, double *time_available)
{
	// compare the total time spent so far with the estimate of how long it 
	// would take to finish using qs and decide whether or not to switch over
	// to qs.
	int decision, sizeN;
	double qs_est_time, nfs_est_time;
	
	// if the size of N is small enough, always switch to qs
	sizeN = zBits(N);
	if (sizeN < 135)
	{
		*time_available = 0;
		decision = 1;
	}
	else
	{
		qs_est_time = get_qs_time_estimate(MEAS_CPU_FREQUENCY,ndigits(N));
		if (VFLAG >= 2)
			printf("***** qs time estimate = %f seconds\n",qs_est_time);

		nfs_est_time = get_gnfs_time_estimate(MEAS_CPU_FREQUENCY,ndigits(N));
		if (VFLAG >= 2)
			printf("***** gnfs time estimate = %f seconds\n",nfs_est_time);

		//proceed with whichever estimate is smaller
		if (qs_est_time <= nfs_est_time)
		{		
			// if qs_est_time is very large, then we don't have a good estimate.  flag the caller 
			// of this fact
			if (qs_est_time > 1e9)
			{
				decision = 0;
				*time_available = -1;
			}		
			else if (total_time > TARGET_ECM_QS_RATIO * qs_est_time)
			{
				// if the total time we've spent so far is greater than a fraction of the time
				// we estimate it would take QS to finish, switch to qs.  
				decision = 1;
				*time_available = 0;
			}
			else
			{
				// otherwise, return the amount of time we have left before the switchover.
				decision = 0;
				*time_available = (TARGET_ECM_QS_RATIO * qs_est_time) - total_time;
			}
		}
		else
		{
			if (nfs_est_time > 1e9)
			{
				decision = 0;
				*time_available = -1;
			}		
			else if (total_time > TARGET_ECM_GNFS_RATIO * nfs_est_time)
			{
				// if the total time we've spent so far is greater than a fraction of the time
				// we estimate it would take gnfs to finish, switch to qs.  
				decision = 2;
				*time_available = 0;
			}
			else
			{
				// otherwise, return the amount of time we have left before the switchover.
				decision = 0;
				*time_available = (TARGET_ECM_GNFS_RATIO * nfs_est_time) - total_time;
			}
		}

	}

	return decision;
}

enum factorization_state scale_requested_work(method_timing_t *method_times, 
	enum factorization_state fact_state, int *next_work, double time_available, z *N)
{
	enum factorization_state new_state = fact_state;
	double base_time_per_curve, time_per_curve;
	int default_curves_15digit = 25;
	int default_curves_20digit = 90;
	int default_curves_25digit = 200;
	int default_curves_30digit = 400;
	int default_curves_35digit = 1000;

	if (time_available < 0)
	{
		// this is a flag from the switch_to_qs function indicating that we
		// can't make a good estimate about how long qs will take.  thus we also
		// don't know how much time is left.  we have to make a decision about how
		// much and what kind of work to do based on the size of the input and our
		// current factorization state.
		int sizeN = zBits(N);

		switch (fact_state)
		{
		case state_trialdiv:			
		case state_fermat:
		case state_rho:
		case state_pp1_lvl1:
		default:
			// any of these states requested: just do it
			break;
		
		case state_pp1_lvl2:
			// equivalent to ecm @ 30 digits
			if (sizeN < 300)
				new_state = state_qs;
			break;

		case state_pp1_lvl3:
			// equivalent to ecm @ 35 digits
			if (sizeN < 330)
				new_state = state_qs;
			break;

		case state_pm1_lvl1:
			if (sizeN < 160)
				new_state = state_qs;
			break;

		case state_pm1_lvl2:
			// equivalent to ecm @ 30 digits
			if (sizeN < 300)
				new_state = state_qs;
			break;

		case state_pm1_lvl3:
			// equivalent to ecm @ 35 digits
			if (sizeN < 330)
				new_state = state_qs;
			break;

		case state_ecm_15digit:
			if (sizeN < 180)
				new_state = state_qs;
			else
				*next_work = default_curves_15digit;
			break;

		case state_ecm_20digit:
			if (sizeN < 220)
				new_state = state_qs;
			else
				*next_work = default_curves_20digit;
			break;

		case state_ecm_25digit:
			if (sizeN < 260)
				new_state = state_qs;
			else
				*next_work = default_curves_25digit;
			break;

		case state_ecm_30digit:
			if (sizeN < 300)
				new_state = state_qs;
			else
				*next_work = default_curves_30digit;
			break;

		case state_ecm_35digit:
			if (sizeN < 330)
				new_state = state_qs;
			else
				*next_work = default_curves_35digit;
			break;

		case state_ecm_auto_increasing:
			if (ndigits(N) < 120)
				new_state = state_qs;
			else
				*next_work = auto_increasing_curves;
			break;

		}
	}
	else
	{
		// use the timing information recorded for previous states to decide how much
		// work we can do in the requested state		

		switch (fact_state)
		{
		case state_trialdiv:			
		case state_fermat:
		case state_rho:
		case state_pp1_lvl1:		
		case state_pm1_lvl1:
		default:
			// any of these states requested: just do it
			break;

		case state_pp1_lvl2:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->pp1_lvl1_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 10;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > 3)
				*next_work = 3;

			if (VFLAG >= 2)
				printf("***** estimating %d more curves can "
					"be run at p+1 level 2\n",*next_work);
			break;

		case state_pp1_lvl3:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->pp1_lvl2_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 10;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > 3)
				*next_work = 3;

			if (VFLAG >= 2)
				printf("***** estimating %d more curves can "
					"be run at p+1 level 3\n",*next_work);
			break;

		case state_pm1_lvl2:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->pm1_lvl1_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 10;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > 1)
				*next_work = 1;

			if (VFLAG >= 2)
				printf("***** estimating %d more curve can "
					"be run at p-1 level 2\n",*next_work);
			break;

		case state_pm1_lvl3:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->pm1_lvl2_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 10;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > 1)
				*next_work = 1;

			if (VFLAG >= 2)
				printf("***** estimating %d more curve can "
					"be run at p-1 level 3\n",*next_work);
			break;

		case state_ecm_15digit:
			*next_work = default_curves_15digit;
			break;

		case state_ecm_20digit:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->ecm_15digit_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 11000 / 2000;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > default_curves_20digit)
				*next_work = default_curves_20digit;

			if (VFLAG >= 2)
				printf("***** estimating %d more curves can "
					"be run at 20 digit level\n",*next_work);

			break;

		case state_ecm_25digit:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->ecm_20digit_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 50000 / 11000;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > default_curves_25digit)
				*next_work = default_curves_25digit;

			if (VFLAG >= 2)
				printf("***** estimating %d more curves can "
					"be run at 25 digit level\n",*next_work);

			break;

		case state_ecm_30digit:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->ecm_25digit_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 250000 / 50000;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > default_curves_30digit)
				*next_work = default_curves_30digit;

			if (VFLAG >= 2)
				printf("***** estimating %d more curves can "
					"be run at 30 digit level\n",*next_work);

			break;

		case state_ecm_35digit:
			//estimate time per curve for this ecm level
			base_time_per_curve = method_times->ecm_30digit_time_per_curve;
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 1000000 / 250000;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > default_curves_35digit)
				*next_work = default_curves_35digit;

			if (VFLAG >= 2)
				printf("***** estimating %d more curves can "
					"be run at 35 digit level\n",*next_work);

			break;

		case state_ecm_auto_increasing:
			//estimate time per curve for this ecm level
			if (method_times->ecm_autoinc_time_per_curve == 0)
			{
				base_time_per_curve = method_times->ecm_35digit_time_per_curve;
				auto_increasing_curves = default_curves_35digit * 3;
				auto_increasing_B1 = 10000000;
				auto_increasing_B2 = 1000000000;
			}
			else
			{
				base_time_per_curve = method_times->ecm_autoinc_time_per_curve;
				auto_increasing_curves *= 3;
				auto_increasing_B1 *= 10;
				auto_increasing_B2 *= 10;
			}
					
			//first scale the base time by the next level ECM
			//curve size
			time_per_curve = base_time_per_curve * 10;

			// estimate curves we can do
			*next_work = (int)(time_available / time_per_curve);

			// sanity check
			if (*next_work <= 0)
				*next_work = 1;
			else if (*next_work > auto_increasing_curves)
				*next_work = auto_increasing_curves;

			if (VFLAG >= 2)
				printf("***** estimating %d more curves can "
					"be run at auto increasing level\n",*next_work);

			break;

		}

	}

	return new_state;
}

void factor(fact_obj_t *fobj)
{
	//run a varity of factoring algorithms on b
	//return any composite number left over
	//the factoring routines will build up a list of factors

	z *b;
	z origN;
	z copyN;
	enum factorization_state fact_state = state_idle;
	method_timing_t method_times;
	int curves = 1;
	int decision = 0;
	int done = 0;	
	FILE *flog;
	struct timeval start, stop;
	double t_time;
	TIME_DIFF *	difference;

	zInit(&origN);
	zInit(&copyN);
	zCopy(&fobj->N,&origN);
	zCopy(&fobj->N,&copyN);
	b = &copyN;

	if (zCompare(b,&zOne) <= 0)
	{
		zFree(&copyN);
		zFree(&origN);
		return;
	}
	
	gettimeofday(&start, NULL);

	flog = fopen(flogname,"a");
	logprint(flog,"\n");
	logprint(flog,"****************************\n");
	logprint(flog,"Starting factorization of %s\n",z2decstr(b,&gstr1));
	logprint(flog,"****************************\n");
	fclose(flog);

	AUTO_FACTOR=1;

	if (VFLAG >= 0)
		printf("factoring %s\n\n",z2decstr(b,&gstr2));

	// initialize time per curve
	method_times.ecm_autoinc_time_per_curve = 0;
	
	while (!done)
	{
		switch (fact_state)
		{
		case state_idle:
			// if this is the top level call, reset the total factorization time
			if (refactor_depth == 0)
				total_time = 0;

			// start things up
			fact_state = state_trialdiv;
			break;

		case state_trialdiv:
			curves = 1;
			t_time = do_work(trialdiv_work, 10000, -1, &curves, b, fobj);
			method_times.trialdiv_time = t_time;
			total_time += t_time * curves;
			fact_state = state_fermat;
			break;

		case state_fermat:
			curves = 1;
			t_time = do_work(fermat_work, FMTMAX, -1, &curves, b, fobj);
			method_times.fermat_time = t_time;
			total_time += t_time * curves;
			fact_state = state_rho;
			break;

		case state_rho:
			curves = 1;
			t_time = do_work(rho_work, -1, -1, &curves, b, fobj);
			method_times.rho_time = t_time;
			total_time += t_time * curves;
			fact_state = state_pp1_lvl1;			
			break;

		case state_pp1_lvl1:
			curves = 3;
			t_time = do_work(pp1_curve, 20000, 2000000, &curves, b, fobj);
			method_times.pp1_lvl1_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_pm1_lvl1;
			break;

		case state_pp1_lvl2:
			t_time = do_work(pp1_curve, 200000, 20000000, &curves, b, fobj);
			method_times.pp1_lvl2_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_pm1_lvl2;
			break;

		case state_pp1_lvl3:
			t_time = do_work(pp1_curve, 2000000, 200000000, &curves, b, fobj);
			method_times.pp1_lvl3_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_pm1_lvl3;
			break;

		case state_pm1_lvl1:
			curves = 1;
			t_time = do_work(pm1_curve, 100000, 10000000, &curves, b, fobj);
			method_times.pm1_lvl1_time_per_curve = t_time;
			total_time += t_time * curves;
			if (NO_ECM)
			{
				if (QS_GNFS_XOVER > 0)				
					decision = QS_GNFS_XOVER;
				else
					decision = 95;

				if (ndigits(b) > decision)
					fact_state = state_nfs;
				else
					fact_state = state_qs;
			}
			else
				fact_state = state_ecm_15digit;

			break;

		case state_pm1_lvl2:
			t_time = do_work(pm1_curve, 1000000, 100000000, &curves, b, fobj);
			method_times.pm1_lvl2_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_30digit;
			break;

		case state_pm1_lvl3:
			t_time = do_work(pm1_curve, 10000000, 1000000000, &curves, b, fobj);
			method_times.pm1_lvl3_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_35digit;
			break;

		case state_ecm_15digit:
			t_time = do_work(ecm_curve, 2000, 200000, &curves, b, fobj);
			method_times.ecm_15digit_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_20digit;			
			break;

		case state_ecm_20digit:
			t_time = do_work(ecm_curve, 11000, 1100000, &curves, b, fobj);
			method_times.ecm_20digit_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_25digit;
			break;

		case state_ecm_25digit:
			t_time = do_work(ecm_curve, 50000, 5000000, &curves, b, fobj);
			method_times.ecm_25digit_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_pp1_lvl2;
			break;

		case state_ecm_30digit:
			t_time = do_work(ecm_curve, 250000, 25000000, &curves, b, fobj);
			method_times.ecm_30digit_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_pp1_lvl3;
			break;

		case state_ecm_35digit:
			t_time = do_work(ecm_curve, 1000000, 100000000, &curves, b, fobj);
			method_times.ecm_35digit_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_auto_increasing;
			break;

		case state_ecm_auto_increasing:
			t_time = do_work(ecm_curve, auto_increasing_B1, auto_increasing_B2, &curves, b, fobj);
			method_times.ecm_autoinc_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_auto_increasing;
			break;

		case state_qs:
			curves = 1;
			t_time = do_work(qs_work, -1, -1, &curves, b, fobj);
			if (VFLAG > 0)
				printf("ECM/SIQS ratio was = %f\n",total_time/t_time);
			total_time += t_time * curves;
			break;

		case state_nfs:			
			curves = 1;
			t_time = do_work(nfs_work, -1, -1, &curves, b, fobj);
			if (VFLAG > 0)
				printf("ECM/NFS ratio was = %f\n",total_time/t_time);
			total_time += t_time * curves;
			break;

		}

		// first, check if we're done
		done = check_if_done(fobj, &origN);		

		if (!done)
		{
			// if we're not done, decide what to do next: either the default next state or
			// switch to qs.
			decision = switch_to_qs(b, &t_time);
			//printf("total time = %f\n",total_time);
			//printf("time available = %f\n",t_time);
			if (decision)
			{
				// either the number is small enough to finish off, or it would be better, 
				// timewise, to switch
				if (decision == 2)
					fact_state = state_nfs;
				else
					fact_state = state_qs;
			}
			else
			{
				// continue with the next default state.
				// use the amount of time available to determine how much work to schedule
				// for the next state
				fact_state = scale_requested_work(&method_times, fact_state, &curves, t_time, b);

			}
		}
	}
		
	if (!isOne(b))
	{
		b->type = COMPOSITE;
		flog = fopen(flogname,"a");
		logprint(flog,"c%d cofactor = %s\n",ndigits(b),z2decstr(b,&gstr1));
		fclose(flog);
		add_to_factor_list(fobj,b);
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
	zFree(&origN);
	zFree(&copyN);
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
		for (i=0; i<fobj->num_factors; i++)
		{
			if (fobj->fobj_factors[i].factor.type == COMPOSITE)
			{
				printf("composite factor detected, refactoring\n");
				zCopy(&fobj->fobj_factors[i].factor,&tmp);
				zCopy(&tmp,&fobj->N);
				factor(fobj);
				zCopy(&fobj->N,&tmp);
				if (zCompare(&tmp,&fobj->fobj_factors[i].factor) == 0)
				{
					printf("refactoring failed\n");
					exit(-1);
				}
			}
		}

		num_fact = 0;
		for (i=0; i<fobj->num_factors; i++)
		{
			if (fobj->fobj_factors[i].factor.type == COMPOSITE)
				continue;

			num_fact += fobj->fobj_factors[i].count;
		}

		zCopy(input,&w1);
		factors = (z *)realloc(factors,num_fact * sizeof(z));
		k=0;
		for (i=0; i<fobj->num_factors; i++)
		{
			if (fobj->fobj_factors[i].factor.type == COMPOSITE)
				continue;

			for (j=0;j<fobj->fobj_factors[i].count;j++)
			{
				zInit(&factors[k]);
				zCopy(&fobj->fobj_factors[i].factor,&factors[k]);
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
						proper = (z *)realloc(proper,proper_alloc * sizeof(z));
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
					proper = (z *)realloc(proper,proper_alloc * sizeof(z));
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

		//free_factor_list(fobj);

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
				seq = (z *)realloc(seq,seq_alloc * sizeof(z));
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

