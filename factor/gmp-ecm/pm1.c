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
       				   --bbuhrow@gmail.com 11/2/10
----------------------------------------------------------------------*/

#include "yafu.h"
#include "yafu_ecm.h"
#include "soe.h"
#include "factor.h"
#include "util.h"
#include <gmp_xface.h>
#include <ecm.h>


// these are used by the top level function, so both YAFU and GMP-ECM
// paths must use these prototypes
void pm1_init(fact_obj_t *fobj);
void pm1_finalize(fact_obj_t *fobj);
void pm1exit(int sig);
int pm1_wrapper(fact_obj_t *fobj);
void pm1_print_B1_B2(fact_obj_t *fobj, FILE *flog);


uint64 TMP_STG2_MAX;

typedef struct
{
	mpz_t gmp_n, gmp_factor;
	ecm_params params;
	uint32 sigma;
	int stagefound;

} ecm_pm1_data_t;

ecm_pm1_data_t pm1_data;

void pm1_init(fact_obj_t *fobj)
{
	mpz_init(pm1_data.gmp_n);
	mpz_init(pm1_data.gmp_factor);
	ecm_init(pm1_data.params);
	//gmp_randseed_ui(tdata->params->rng, 
	//	get_rand(&obj->seed1, &obj->seed2));

	pm1_data.params->method = ECM_PM1;
	//pm1_data.params->verbose = 1;

	TMP_STG2_MAX = fobj->pm1_obj.B2;

	return;
}

void pm1_finalize(fact_obj_t *fobj)
{
	ecm_clear(pm1_data.params);
	mpz_clear(pm1_data.gmp_n);
	mpz_clear(pm1_data.gmp_factor);

	fobj->pm1_obj.B2 = TMP_STG2_MAX;
	
	return;
}

int pm1_wrapper(fact_obj_t *fobj)
{
	z *n, *f;	
	int status;
#if defined(_WIN64) && BITS_PER_DIGIT == 32
	size_t count;
#endif

	f = &fobj->pp1_obj.factors[0];
	n = &fobj->pp1_obj.n;

	pm1_data.params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
	if (VFLAG >= 3)
		pm1_data.params->verbose = VFLAG - 2;		

#if defined(_WIN64) && BITS_PER_DIGIT == 32
	mpz_import(pm1_data.gmp_n, (size_t)(abs(n->size)), -1, sizeof(uint32), 
		0, (size_t)0, n->val);
#else
	//wrapper for YAFU bigints and call to gmp-ecm
	mp2gmp(n, pm1_data.gmp_n);
#endif

	if (fobj->pm1_obj.stg2_is_default == 0)
	{
		//not default, tell gmp-ecm to use the requested B2
		//printf("using requested B2 value\n");
		sp642z(fobj->pm1_obj.B2,f);
		mp2gmp(f,pm1_data.params->B2);
		zClear(f);
	}

	status = ecm_factor(pm1_data.gmp_factor, pm1_data.gmp_n,
			fobj->pm1_obj.B1, pm1_data.params);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
	zClear(n);
	mpz_export(n->val, &count, -1, sizeof(uint32),
			0, (size_t)0, pm1_data.gmp_n);
	n->size = count;
#else
	//update n: not sure if gmp-ecm modifies it
	gmp2mp(pm1_data.gmp_n, n);
#endif

	//NOTE: this required a modification to the GMP-ECM source code in pm1.c
	//in order to get the automatically computed B2 value out of the
	//library
	//gmp2mp(pm1_data.params->B2,f);
	//POLLARD_STG2_MAX = z264(f);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
	zClear(f);
	mpz_export(f->val, &count, -1, sizeof(uint32),
			0, (size_t)0, pm1_data.gmp_factor);
	f->size = count;
#else
	//pull out any factor found
	gmp2mp(pm1_data.gmp_factor,f);
#endif

	//the return value is the stage the factor was found in, if no error
	pm1_data.stagefound = status;

	return status;
}


// top level routine: the only one visible to the rest of the program
void pollard_loop(fact_obj_t *fobj)
{
	//run pollard's p-1 algorithm once on the input, using a 
	//32 bit random base
	z *n = &fobj->pm1_obj.n;
	z d,*f,t;   
	//int i;
	uint32 base;
	FILE *flog;
	clock_t start, stop;
	double tt;

	//check for trivial cases
	if (isOne(n) || isZero(n))
	{
		n->type = COMPOSITE;
		return;
	}
	if (zCompare(n,&zTwo) == 0)
	{
		n->type = PRIME;
		return;
	}	

	//open the log file
	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->flogname);
		return;
	}
		
	start = clock();

	if (isPrime(n))
	{
		n->type = PRP;
		add_to_factor_list(fobj, n);
		logprint(flog,"prp%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
		stop = clock();
		tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		zCopy(&zOne,n);
		return;
	}

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	PM1_ABORT = 0;
	signal(SIGINT,pm1exit);

	zInit(&d);
	fobj->pm1_obj.num_factors = 1;
	zInit(&fobj->pm1_obj.factors[0]);
	f = &fobj->pm1_obj.factors[0];
	zInit(&t);

	pm1_init(fobj);
		
	//base = spRand(3,0xFFFFFFFF);

	pm1_print_B1_B2(fobj,flog);
	pm1_wrapper(fobj);
		
	if (zCompare(f,&zOne) > 0 && zCompare(f,n) < 0)
	{
		//non-trivial factor found
		stop = clock();
		tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		if (isPrime(f))
		{
			f->type = PRP;
			add_to_factor_list(fobj, f);
			//log result
			if (VFLAG > 0)
				printf("pm1: found prp%d factor = %s\n",ndigits(f),z2decstr(f,&gstr1));
			logprint(flog,"prp%d = %s\n",
				ndigits(f),z2decstr(f,&gstr2));
		}
		else
		{
			f->type = COMPOSITE;
			add_to_factor_list(fobj, f);
			//log result
			if (VFLAG > 0)
					printf("pm1: found c%d factor = %s\n",ndigits(f),z2decstr(f,&gstr1));
			logprint(flog,"c%d = %s\n",
				ndigits(f),z2decstr(f,&gstr2));
		}
		start = clock();

		//reduce input
		zDiv(n,f,&t,&d);
		zCopy(&t,n);
	}

	fclose(flog);

	pm1_finalize(fobj);

	//watch for an abort
	if (PM1_ABORT)
	{
		print_factors(fobj);
		exit(1);
	}

	signal(SIGINT,NULL);
	zFree(&d);
	zFree(&t);
	return;
}

void pm1_print_B1_B2(fact_obj_t *fobj, FILE *flog)
{
	z *n = &fobj->pm1_obj.n;
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (fobj->pm1_obj.B1 % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",fobj->pm1_obj.B1 / 1000000000, suffix);
	}
	else if (fobj->pm1_obj.B1 % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",fobj->pm1_obj.B1 / 1000000, suffix);
	}
	else if (fobj->pm1_obj.B1 % 1000 == 0)
	{
		suffix = 'K';
		sprintf(stg1str,"%u%c",fobj->pm1_obj.B1 / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",fobj->pm1_obj.B1);
	}

	if (fobj->pm1_obj.stg2_is_default == 0)
	{
		if (fobj->pm1_obj.B2 % 1000000000 == 0)
		{
			suffix = 'B';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->pm1_obj.B2 / 1000000000, suffix);
		}
		else if (fobj->pm1_obj.B2 % 1000000 == 0)
		{
			suffix = 'M';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->pm1_obj.B2 / 1000000, suffix);
		}
		else if (fobj->pm1_obj.B2 % 1000 == 0)
		{
			suffix = 'K';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->pm1_obj.B2 / 1000, suffix);
		}
		else
		{
			sprintf(stg2str,"%" PRIu64 "",fobj->pm1_obj.B2);
		}
	}
	else
		sprintf(stg2str, "gmp-ecm default");

	if (VFLAG >= 0)
	{
		printf("pm1: starting B1 = %s, B2 = %s on C%d",
			stg1str,stg2str,ndigits(n));
		fflush(stdout);	
	}
	logprint(flog,"pm1: starting B1 = %s, B2 = %s on C%d\n",
		stg1str,stg2str,ndigits(n));

#if defined(HAVE_GMP) && defined(HAVE_GMP_ECM)
		//need a new line to make screen output look right, when
		//using GMP-ECM, because the "processed" status is not printed
	if (VFLAG >= 0)
		printf("\n");
#endif
	
	return;
}

void pm1exit(int sig)
{
	printf("\nAborting...\n");
	PM1_ABORT = 1;
	return;
}

