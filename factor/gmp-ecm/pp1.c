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

#include "yafu_ecm.h"
#include "soe.h"
#include "factor.h"
#include "ytools.h"
#include <gmp_xface.h>
#include <ecm.h>
#include <signal.h>
#include <math.h>

// these are used by the top level function, so both YAFU and GMP-ECM
// paths must use these prototypes
void pp1_init(fact_obj_t *fobj);
void pp1_finalize(fact_obj_t *fobj);
void pp1_print_B1_B2(fact_obj_t *fobj);
int pp1_wrapper(fact_obj_t *fobj);
void pp1exit(int sig);

uint64_t PP1_TMP_STG2_MAX;
int PP1_ABORT;

typedef struct
{
	mpz_t gmp_n, gmp_factor;
	ecm_params params;
	uint32_t sigma;
	int stagefound;

} ecm_pp1_data_t;

ecm_pp1_data_t pp1_data;

void pp1_init(fact_obj_t *fobj)
{
	mpz_init(pp1_data.gmp_n);
	mpz_init(pp1_data.gmp_factor);
	ecm_init(pp1_data.params);
	pp1_data.params->method = ECM_PP1;
	PP1_TMP_STG2_MAX = fobj->pp1_obj.B2;
	PP1_ABORT = 0;
	signal(SIGINT,pp1exit);

	return;
}

void pp1_finalize(fact_obj_t *fobj)
{
	ecm_clear(pp1_data.params);
	mpz_clear(pp1_data.gmp_n);
	mpz_clear(pp1_data.gmp_factor);
	fobj->pp1_obj.B2 = PP1_TMP_STG2_MAX;
	signal(SIGINT,NULL);

	return;
}

int pp1_wrapper(fact_obj_t *fobj)
{
	int status;

	// timing variables
	struct timeval stopt;	// stop time of this job
	struct timeval startt;	// start time of this job
	double t_time;

	gettimeofday(&startt, NULL);

	mpz_set(pp1_data.gmp_n, fobj->pp1_obj.gmp_n);

	pp1_data.params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
	if (fobj->VFLAG >= 3)
		pp1_data.params->verbose = fobj->VFLAG - 2;		

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

	gettimeofday(&stopt, NULL);
	t_time = ytools_difftime(&startt, &stopt);
	fobj->pp1_obj.ttime = t_time;

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
    clock_t start = clock();
    clock_t stop = clock();
	double tt;
		
	//check for trivial cases
	if ((mpz_cmp_ui(fobj->pp1_obj.gmp_n, 1) == 0) || (mpz_cmp_ui(fobj->pp1_obj.gmp_n, 0) == 0))
		return;

	if (mpz_cmp_ui(fobj->pp1_obj.gmp_n, 2) == 0)
		return;

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	//PP1_ABORT = 0;
	//signal(SIGINT,pp1exit);

	//initialize some local args
	mpz_init(d);
	mpz_init(t);	


#if !defined(USE_AVX512F) || defined(__MINGW32__)
    if (1)
#else
    if ((fobj->ecm_obj.prefer_gmpecm) ||
        (!fobj->HAS_AVX512F))
#endif
    {
        // prepare for gmp-ecm p+1
        pp1_init(fobj);

        i = 0;
        while (i < trials)
        {
            //watch for an abort
            if (PP1_ABORT)
            {
                print_factors(fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES);
                exit(1);
            }

            start = clock();
            if (is_mpz_prp(fobj->pp1_obj.gmp_n, fobj->NUM_WITNESSES))
            {
                char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_n);
                logprint_oc(fobj->flogname, "a", "prp%d = %s\n", gmp_base10(fobj->pp1_obj.gmp_n), s);
                free(s);
                add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_n,
                    fobj->VFLAG, fobj->NUM_WITNESSES);

                stop = clock();
                tt = (double)(stop - start) / (double)CLOCKS_PER_SEC;
                mpz_set_ui(fobj->pp1_obj.gmp_n, 1);
                break;
            }

            fobj->pp1_obj.base = lcg_rand_32_range(3, 0xffffffff, &fobj->pp1_obj.lcg_state);

            pp1_print_B1_B2(fobj);

            it = pp1_wrapper(fobj);

            //check to see if 'f' is non-trivial
            if ((mpz_cmp_ui(fobj->pp1_obj.gmp_f, 1) > 0)
                && (mpz_cmp(fobj->pp1_obj.gmp_f, fobj->pp1_obj.gmp_n) < 0))
            {
                //non-trivial factor found
                stop = clock();
                tt = (double)(stop - start) / (double)CLOCKS_PER_SEC;

                //check if the factor is prime
                if (is_mpz_prp(fobj->pp1_obj.gmp_f, fobj->NUM_WITNESSES))
                {
                    add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_f,
                        fobj->VFLAG, fobj->NUM_WITNESSES);

                    if (fobj->VFLAG > 0)
                        gmp_printf("pp1: found prp%d factor = %Zd\n",
                            gmp_base10(fobj->pp1_obj.gmp_f), fobj->pp1_obj.gmp_f);

                    char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_f);
                    logprint_oc(fobj->flogname, "a", "prp%d = %s\n",
                        gmp_base10(fobj->pp1_obj.gmp_f), s);
                    free(s);
                }
                else
                {
                    add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_f,
                        fobj->VFLAG, fobj->NUM_WITNESSES);

                    if (fobj->VFLAG > 0)
                        gmp_printf("pp1: found c%d factor = %Zd\n",
                            gmp_base10(fobj->pp1_obj.gmp_f), fobj->pp1_obj.gmp_f);

                    char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_f);
                    logprint_oc(fobj->flogname, "a", "c%d = %s\n",
                        gmp_base10(fobj->pp1_obj.gmp_f), s);
                    free(s);
                }
                start = clock();

                //reduce input
                mpz_tdiv_q(fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_n, fobj->pp1_obj.gmp_f);

                i++;
                break;
            }

            i++;
        }

        // clean up after gmp-ecm P+1
        pp1_finalize(fobj);
    }
    else
    {
        // run vecPPM1
        if (is_mpz_prp(fobj->pp1_obj.gmp_n, fobj->NUM_WITNESSES))
        {
            char* s = mpz_get_str(NULL, 10, fobj->pp1_obj.gmp_n);
            logprint_oc(fobj->flogname, "a", "prp%d = %s\n", gmp_base10(fobj->pp1_obj.gmp_n), s);
            free(s);
            add_to_factor_list(fobj->factors, fobj->pp1_obj.gmp_n,
                fobj->VFLAG, fobj->NUM_WITNESSES);

            stop = clock();
            tt = (double)(stop - start) / (double)CLOCKS_PER_SEC;
            mpz_set_ui(fobj->pp1_obj.gmp_n, 1);
        }
        else
        {
            pp1_print_B1_B2(fobj);

            vecPP1(fobj);

        }
    }

	//signal(SIGINT,NULL);
	mpz_clear(d);
	mpz_clear(t);

	return;
}

void pp1_print_B1_B2(fact_obj_t *fobj)
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

	if (fobj->VFLAG >= 0)
	{
		printf("pp1: starting B1 = %s, B2 = %s on C%d"
			,stg1str,stg2str, (int)gmp_base10(fobj->pp1_obj.gmp_n));
		fflush(stdout);	
	}

	logprint_oc(fobj->flogname, "a", "pp1: starting B1 = %s, B2 = %s on C%d\n"
		,stg1str,stg2str, (int)gmp_base10(fobj->pp1_obj.gmp_n));

		//need a new line to make screen output look right, when
		//using GMP-ECM, because the "processed" status is not printed
	if (fobj->VFLAG >= 0)
		printf("\n");

	return;
}

void pp1exit(int sig)
{
	printf("\nAborting...\n");
	PP1_ABORT = 1;
	return;
}

