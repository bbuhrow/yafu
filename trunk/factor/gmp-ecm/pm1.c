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
	int status;

	mpz_set(pm1_data.gmp_n, fobj->pm1_obj.mpz_n);

	pm1_data.params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
	if (VFLAG >= 3)
		pm1_data.params->verbose = VFLAG - 2;		

	if (fobj->pm1_obj.stg2_is_default == 0)
	{
		//not default, tell gmp-ecm to use the requested B2
		//printf("using requested B2 value\n");
		uint64_2gmp(fobj->pm1_obj.B2, pm1_data.params->B2);
	}

	status = ecm_factor(pm1_data.gmp_factor, pm1_data.gmp_n,
			fobj->pm1_obj.B1, pm1_data.params);

	mpz_set(fobj->pm1_obj.mpz_n, pm1_data.gmp_n);

	//NOTE: this required a modification to the GMP-ECM source code in pp1.c
	//in order to get the automatically computed B2 value out of the
	//library
	//gmp2mp(pp1_data.params->B2,f);
	//WILL_STG2_MAX = z264(f);

	mpz_set(fobj->pm1_obj.mpz_f, pm1_data.gmp_factor);

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
	mpz_t d,t;
	FILE *flog;
	clock_t start, stop;
	double tt;
	z f;
		
	mp2gmp(n, fobj->pm1_obj.mpz_n);

	//check for trivial cases
	if ((mpz_cmp_ui(fobj->pm1_obj.mpz_n, 1) == 0) || (mpz_cmp_ui(fobj->pm1_obj.mpz_n, 0) == 0))
	{
		n->type = COMPOSITE;
		return;
	}

	if (mpz_cmp_ui(fobj->pm1_obj.mpz_n, 2) == 0)
	{
		n->type = PRIME;
		return;
	}	
		
	start = clock();

	//open the log file
	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->flogname);
		return;
	}

	if (mpz_probab_prime_p(fobj->pm1_obj.mpz_n, NUM_WITNESSES))
	{
		zInit(&f);

		logprint(flog,"prp%d = %s\n",mpz_sizeinbase(fobj->pm1_obj.mpz_n,10),
			mpz_get_str(gstr1.s, 10, fobj->pm1_obj.mpz_n));

		gmp2mp(fobj->pm1_obj.mpz_n, &f);
		add_to_factor_list(fobj, &f);

		stop = clock();
		tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;			
		mpz_set_ui(fobj->pm1_obj.mpz_n, 1);
		gmp2mp(fobj->pm1_obj.mpz_n, n);
		
		zFree(&f);
		return;
	}

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	PM1_ABORT = 0;
	signal(SIGINT,pm1exit);

	//initialize some local args
	mpz_init(d);
	mpz_init(t);	
	zInit(&f);

	pm1_init(fobj);
		
	pm1_print_B1_B2(fobj,flog);
	pm1_wrapper(fobj);
		
	//check to see if 'f' is non-trivial
	if ((mpz_cmp_ui(fobj->pm1_obj.mpz_f, 1) > 0)
		&& (mpz_cmp(fobj->pm1_obj.mpz_f, fobj->pm1_obj.mpz_n) < 0))
	{				
		gmp2mp(fobj->pm1_obj.mpz_f, &f);

		//non-trivial factor found
		stop = clock();
		tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		//check if the factor is prime
		if (mpz_probab_prime_p(fobj->pm1_obj.mpz_f, NUM_WITNESSES))
		{
			add_to_factor_list(fobj, &f);

			if (VFLAG > 0)
				gmp_printf("pm1: found prp%d factor = %Zd\n",
				mpz_sizeinbase(fobj->pm1_obj.mpz_f, 10),fobj->pm1_obj.mpz_f);

			logprint(flog,"prp%d = %s\n",
				mpz_sizeinbase(fobj->pm1_obj.mpz_f, 10),
				mpz_get_str(gstr1.s, 10, fobj->pm1_obj.mpz_f));
		}
		else
		{
			add_to_factor_list(fobj, &f);
			if (VFLAG > 0)
				gmp_printf("pm1: found c%d factor = %Zd\n",
				mpz_sizeinbase(fobj->pm1_obj.mpz_f, 10),fobj->pm1_obj.mpz_f);

			logprint(flog,"c%d = %s\n",
				mpz_sizeinbase(fobj->pm1_obj.mpz_f, 10),
				mpz_get_str(gstr1.s, 10, fobj->pm1_obj.mpz_f));
		}
		start = clock();

		//reduce input
		mpz_tdiv_q(fobj->pm1_obj.mpz_n, fobj->pm1_obj.mpz_n, fobj->pm1_obj.mpz_f);
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
	gmp2mp(fobj->pm1_obj.mpz_n, n);
	mpz_clear(d);
	mpz_clear(t);
	zFree(&f);

	return;
}

void pm1_print_B1_B2(fact_obj_t *fobj, FILE *flog)
{
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
			stg1str,stg2str,mpz_sizeinbase(fobj->pm1_obj.mpz_n, 10));
		fflush(stdout);	
	}
	logprint(flog,"pm1: starting B1 = %s, B2 = %s on C%d\n",
		stg1str,stg2str,mpz_sizeinbase(fobj->pm1_obj.mpz_n, 10));

		//need a new line to make screen output look right, when
		//using GMP-ECM, because the "processed" status is not printed
	if (VFLAG >= 0)
		printf("\n");

	return;
}

void pm1exit(int sig)
{
	printf("\nAborting...\n");
	PM1_ABORT = 1;
	return;
}

