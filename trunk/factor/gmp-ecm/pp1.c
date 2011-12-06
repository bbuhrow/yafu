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
void pp1_init(fact_obj_t *fobj);
void pp1_finalize(fact_obj_t *fobj);
void pp1_print_B1_B2(fact_obj_t *fobj, FILE *flog);
int pp1_wrapper(fact_obj_t *fobj);
void pp1exit(int sig);
uint64 TMP_STG2_MAX;

typedef struct
{
	mpz_t gmp_n, gmp_factor;
	ecm_params params;
	uint32 sigma;
	int stagefound;

} ecm_pp1_data_t;

ecm_pp1_data_t pp1_data;

void pp1_init(fact_obj_t *fobj)
{
	mpz_init(pp1_data.gmp_n);
	mpz_init(pp1_data.gmp_factor);
	ecm_init(pp1_data.params);
	//gmp_randseed_ui(tdata->params->rng, 
	//	get_rand(&obj->seed1, &obj->seed2));

	pp1_data.params->method = ECM_PP1;
	//pp1_data.params->verbose = 1;

	TMP_STG2_MAX = fobj->pp1_obj.B2; //WILL_STG2_MAX;

	PP1_ABORT = 0;
	signal(SIGINT,pp1exit);

	return;
}

void pp1_finalize(fact_obj_t *fobj)
{
	ecm_clear(pp1_data.params);
	mpz_clear(pp1_data.gmp_n);
	mpz_clear(pp1_data.gmp_factor);

	//WILL_STG2_MAX = TMP_STG2_MAX;
	fobj->pp1_obj.B2 = TMP_STG2_MAX;
	signal(SIGINT,NULL);

	return;
}

int pp1_wrapper(fact_obj_t *fobj)
{
	int status;

	mpz_set(pp1_data.gmp_n, fobj->pp1_obj.gmp_n);

	pp1_data.params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
	if (VFLAG >= 3)
		pp1_data.params->verbose = VFLAG - 2;		

	if (fobj->pp1_obj.stg2_is_default == 0)
	{
		//not default, tell gmp-ecm to use the requested B2
		//printf("using requested B2 value\n");
		uint64_2gmp(fobj->pp1_obj.B2, pp1_data.params->B2);
	}

	status = ecm_factor(pp1_data.gmp_factor, pp1_data.gmp_n,
			fobj->pp1_obj.B1, pp1_data.params);

	mpz_set(fobj->pp1_obj.gmp_n, pp1_data.gmp_n);

	//NOTE: this required a modification to the GMP-ECM source code in pp1.c
	//in order to get the automatically computed B2 value out of the
	//library
	//gmp2mp(pp1_data.params->B2,f);
	//WILL_STG2_MAX = z264(f);

	mpz_set(fobj->pp1_obj.gmp_f, pp1_data.gmp_factor);

	//the return value is the stage the factor was found in, if no error
	pp1_data.stagefound = status;

	return status;
}

void williams_loop(fact_obj_t *fobj)
{
	//use william's p+1 algorithm 'trials' times on n.
	//we allow multiple trials since we may not always get a p+1 group
	//with our random choice of base
	//expects the input in pp1_obj->gmp_n
	mpz_t d,t;
	int i,it,trials = fobj->pp1_obj.numbases;
	FILE *flog;
	clock_t start, stop;
	double tt;
		
	//check for trivial cases
	if ((mpz_cmp_ui(fobj->pp1_obj.gmp_n, 1) == 0) || (mpz_cmp_ui(fobj->pp1_obj.gmp_n, 0) == 0))
		return;

	if (mpz_cmp_ui(fobj->pp1_obj.gmp_n, 2) == 0)
		return;

	//open the log file
	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->flogname);
		return;
	}

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	PP1_ABORT = 0;
	signal(SIGINT,pp1exit);

	//initialize some local args
	mpz_init(d);
	mpz_init(t);	

	pp1_init(fobj);

	i=0;
	while (i < trials)
	{
		//watch for an abort
		if (PP1_ABORT)
		{
			print_factors(fobj);
			exit(1);
		}

		start = clock();
		if (mpz_probab_prime_p(fobj->pp1_obj.gmp_n, NUM_WITNESSES))
		{
			logprint(flog,"prp%d = %s\n", gmp_base10(fobj->pp1_obj.gmp_n),
				mpz_get_str(gstr1.s, 10, fobj->pp1_obj.gmp_n));

			add_to_factor_list(fobj, fobj->pp1_obj.gmp_n);

			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;			
			mpz_set_ui(fobj->pp1_obj.gmp_n, 1);
			break;
		}
		
		fobj->pp1_obj.base = spRand(3,MAX_DIGIT);

		pp1_print_B1_B2(fobj,flog);

		it = pp1_wrapper(fobj);
		
		//check to see if 'f' is non-trivial
		if ((mpz_cmp_ui(fobj->pp1_obj.gmp_f, 1) > 0)
			&& (mpz_cmp(fobj->pp1_obj.gmp_f, fobj->pp1_obj.gmp_n) < 0))
		{				
			//non-trivial factor found
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;

			//check if the factor is prime
			if (mpz_probab_prime_p(fobj->pp1_obj.gmp_f, NUM_WITNESSES))
			{
				add_to_factor_list(fobj, fobj->pp1_obj.gmp_f);

				if (VFLAG > 0)
					gmp_printf("pp1: found prp%d factor = %Zd\n",
					gmp_base10(fobj->pp1_obj.gmp_f),fobj->pp1_obj.gmp_f);

				logprint(flog,"prp%d = %s\n",
					gmp_base10(fobj->pp1_obj.gmp_f),
					mpz_get_str(gstr1.s, 10, fobj->pp1_obj.gmp_f));
			}
			else
			{
				add_to_factor_list(fobj, fobj->pp1_obj.gmp_f);

				if (VFLAG > 0)
					gmp_printf("pp1: found c%d factor = %Zd\n",
					gmp_base10(fobj->pp1_obj.gmp_f),fobj->pp1_obj.gmp_f);

				logprint(flog,"c%d = %s\n",
					gmp_base10(fobj->pp1_obj.gmp_f),
					mpz_get_str(gstr1.s, 10, fobj->pp1_obj.gmp_f));
			}
			start = clock();

			//reduce input
			mpz_tdiv_q(fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_f);

			i++;
			break;
		}

		i++;
	}

	fclose(flog);

	pp1_finalize(fobj);
	signal(SIGINT,NULL);
	mpz_clear(d);
	mpz_clear(t);

	return;
}

void pp1_print_B1_B2(fact_obj_t *fobj, FILE *flog)
{
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (fobj->pp1_obj.B1 % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",fobj->pp1_obj.B1 / 1000000000, suffix);
	}
	else if (fobj->pp1_obj.B1 % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",fobj->pp1_obj.B1 / 1000000, suffix);
	}
	else if (fobj->pp1_obj.B1 % 1000 == 0)
	{
		suffix = 'K';
		sprintf(stg1str,"%u%c",fobj->pp1_obj.B1 / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",fobj->pp1_obj.B1);
	}

	if (fobj->pp1_obj.stg2_is_default == 0)
	{
		if (fobj->pp1_obj.B2 % 1000000000 == 0)
		{
			suffix = 'B';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->pp1_obj.B2 / 1000000000, suffix);
		}
		else if (fobj->pp1_obj.B2 % 1000000 == 0)
		{
			suffix = 'M';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->pp1_obj.B2 / 1000000, suffix);
		}
		else if (fobj->pp1_obj.B2 % 1000 == 0)
		{
			suffix = 'K';
			sprintf(stg2str,"%" PRIu64 "%c",fobj->pp1_obj.B2 / 1000, suffix);
		}
		else
		{
			sprintf(stg2str,"%" PRIu64 "",fobj->pp1_obj.B2);
		}
	}
	else
		sprintf(stg2str, "gmp-ecm default");

	if (VFLAG >= 0)
	{
		printf("pp1: starting B1 = %s, B2 = %s on C%d"
			,stg1str,stg2str, (int)gmp_base10(fobj->pp1_obj.gmp_n));
		fflush(stdout);	
	}
	logprint(flog, "pp1: starting B1 = %s, B2 = %s on C%d\n"
		,stg1str,stg2str, (int)gmp_base10(fobj->pp1_obj.gmp_n));

		//need a new line to make screen output look right, when
		//using GMP-ECM, because the "processed" status is not printed
	if (VFLAG >= 0)
		printf("\n");

	return;
}

void pp1exit(int sig)
{
	printf("\nAborting...\n");
	PP1_ABORT = 1;
	return;
}

