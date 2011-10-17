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

// local function to do requested curve based factorization
double get_qs_time_estimate(fact_obj_t *fobj, double freq, int bits);
double get_gnfs_time_estimate(fact_obj_t *fobj, double freq, int digits);

double do_work(enum work_method method, uint32 B1, uint64 B2, int *work, 
	mpz_t b, fact_obj_t *fobj);
int check_if_done(fact_obj_t *fobj, mpz_t N);
enum factorization_state scale_requested_work(method_timing_t *method_times, 
	enum factorization_state fact_state, int *next_work, double time_available, mpz_t N);
int switch_to_qs(fact_obj_t *fobj, mpz_t N, double *time_available, int force_switch);

void init_factobj(fact_obj_t *fobj)
{
	// get space for everything
	alloc_factobj(fobj);

	// initialize global stuff in fobj
	fobj->seed1 = g_rand.low;
	fobj->seed2 = g_rand.hi;
	yafu_get_cache_sizes(&fobj->cache_size1,&fobj->cache_size2);
	fobj->flags = 0;
	fobj->num_threads = 1;		//read from input arguments	
	strcpy(fobj->flogname,"factor.log");	

	// initialize stuff for rho	
	fobj->rho_obj.iterations = 1000;		
	fobj->rho_obj.curr_poly = 0;	

	// initialize stuff for pm1	
	fobj->pm1_obj.B1 = 100000;
	fobj->pm1_obj.B2 = 10000000;
	fobj->pm1_obj.stg2_is_default = 1;	

	// initialize stuff for pp1	
	fobj->pp1_obj.B1 = 20000;
	fobj->pp1_obj.B2 = 1000000;
	fobj->pp1_obj.stg2_is_default = 1;	

	// initialize stuff for ecm	
	fobj->ecm_obj.B1 = 11000;
	fobj->ecm_obj.B2 = 1100000;
	fobj->ecm_obj.stg2_is_default = 1;
	fobj->ecm_obj.sigma = 0;
	fobj->ecm_obj.num_curves = 90;
#ifdef FORK_ECM
	fobj->ecm_obj.curves_run = NULL;
#else
	fobj->ecm_obj.curves_run = 0;
#endif	
	// unlike ggnfs, ecm does not *require* external binaries.  
	// an empty string indicates the use of the built-in GMP-ECM hooks, while
	// a non-empty string (filled in by the user) will indicate the use of
	// an external binary
	strcpy(fobj->ecm_obj.ecm_path,"");
	fobj->ecm_obj.use_external = 0;

	// initialize stuff for squfof
	fobj->squfof_obj.num_factors = 0;	

	// initialize stuff for qs	
	fobj->qs_obj.gbl_override_B_flag = 0;
	fobj->qs_obj.gbl_override_B = 0;
	fobj->qs_obj.gbl_override_blocks_flag = 0;
	fobj->qs_obj.gbl_override_blocks = 0 ;
	fobj->qs_obj.gbl_override_lpmult_flag = 0;
	fobj->qs_obj.gbl_override_lpmult = 0;
	fobj->qs_obj.gbl_override_rel_flag = 0;
	fobj->qs_obj.gbl_override_rel = 0;
	fobj->qs_obj.gbl_override_tf_flag = 0;
	fobj->qs_obj.gbl_override_tf = 0;
	fobj->qs_obj.gbl_override_time_flag = 0;
	fobj->qs_obj.gbl_override_time = 0;
	fobj->qs_obj.flags = 0;
	fobj->qs_obj.gbl_force_DLP = 0;
	fobj->qs_obj.qs_exponent = 0;
	fobj->qs_obj.qs_multiplier = 0;
	fobj->qs_obj.qs_tune_freq = 0;
	fobj->qs_obj.no_small_cutoff_opt = 0;
	strcpy(fobj->qs_obj.siqs_savefile,"siqs.dat");	

	// initialize stuff for trial division	
	fobj->div_obj.print = 0;
	fobj->div_obj.limit = 10000;
	fobj->div_obj.fmtlimit = 1000000;
	
	//initialize stuff for nfs	
	fobj->nfs_obj.gnfs_exponent = 0;
	fobj->nfs_obj.gnfs_multiplier = 0;
	fobj->nfs_obj.gnfs_tune_freq = 0;
	fobj->nfs_obj.min_digits = 85;
	fobj->nfs_obj.siever = 0;							//default, use automatic selection
	fobj->nfs_obj.startq = 0;							//default, not used
	fobj->nfs_obj.rangeq = 0;							//default, not used
	fobj->nfs_obj.polystart = 0;						//default, not used
	fobj->nfs_obj.polyrange = 0;						//default, not used
	strcpy(fobj->nfs_obj.outputfile,"nfs.dat");			//default
	strcpy(fobj->nfs_obj.logfile,"nfs.log");			//default
	strcpy(fobj->nfs_obj.fbfile,"nfs.fb");				//default
	fobj->nfs_obj.sq_side = 1;							//default = algebraic
	fobj->nfs_obj.timeout = 0;							//default, not used
	strcpy(fobj->nfs_obj.job_infile,"nfs.job");			//default
	fobj->nfs_obj.sieve_only = 0;						//default = no
	fobj->nfs_obj.poly_only = 0;						//default = no
	fobj->nfs_obj.post_only = 0;						//default = no
	fobj->nfs_obj.la_restart = 0;						//default = no
	fobj->nfs_obj.poly_option = 0;						//default = fast search
															//1 = wide
															//2 = deep
	fobj->nfs_obj.restart_flag = 0;						//default = not a restart
	fobj->nfs_obj.polybatch = 250;						//default	
#if defined(_WIN64)
	strcpy(fobj->nfs_obj.ggnfs_dir,".\\");
#elif defined(WIN32)
	strcpy(fobj->nfs_obj.ggnfs_dir,".\\");
#else
	strcpy(fobj->nfs_obj.ggnfs_dir,"./");
#endif

	//initialize autofactor object
	//whether we want to output certain info to their own files...
	fobj->autofact_obj.want_output_primes = 0;
	fobj->autofact_obj.want_output_factors = 0;
	fobj->autofact_obj.want_output_unfactored = 0;
	fobj->autofact_obj.want_output_expressions = 1;
	fobj->autofact_obj.qs_gnfs_xover = 0;
	fobj->autofact_obj.want_only_1_factor = 0;
	fobj->autofact_obj.no_ecm = 0;
	fobj->autofact_obj.target_ecm_qs_ratio = 0.25;
	fobj->autofact_obj.target_ecm_gnfs_ratio = 0.25;
	fobj->autofact_obj.target_ecm_snfs_ratio = 0.20;

	//pretesting plan used by factor()
	fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
	strcpy(fobj->autofact_obj.plan_str,"normal");
	fobj->autofact_obj.only_pretest = 0;
	fobj->autofact_obj.autofact_active = 0;

	return;
}

void free_factobj(fact_obj_t *fobj)
{
	uint32 i;

	// free stuff in rho
	free(fobj->rho_obj.polynomials);
	mpz_clear(fobj->rho_obj.gmp_n);
	mpz_clear(fobj->rho_obj.gmp_f);

	// free stuff in pm1
	mpz_clear(fobj->pm1_obj.gmp_n);
	mpz_clear(fobj->pm1_obj.gmp_f);

	// free stuff in pp1
	mpz_clear(fobj->pp1_obj.gmp_n);
	mpz_clear(fobj->pp1_obj.gmp_f);

	// free any factors found in ecm
	for (i=0; i<fobj->ecm_obj.num_factors; i++)
		zFree(&fobj->ecm_obj.factors[i]);
	free(fobj->ecm_obj.factors);

	// then free other stuff in ecm
	mpz_clear(fobj->ecm_obj.gmp_n);
	mpz_clear(fobj->ecm_obj.gmp_f);

	// free any factors found in squfof
	for (i=0; i<fobj->squfof_obj.num_factors; i++)
		zFree(&fobj->squfof_obj.factors[i]);
	free(fobj->squfof_obj.factors);

	// then free other stuff in squfof
	mpz_clear(fobj->squfof_obj.gmp_n);
	mpz_clear(fobj->squfof_obj.gmp_f);

	// free any factors found in qs
	for (i=0; i<fobj->qs_obj.num_factors; i++)
		zFree(&fobj->qs_obj.factors[i]);
	free(fobj->qs_obj.factors);

	// then free other stuff in qs
	mpz_clear(fobj->qs_obj.gmp_n);

	// free any factors found in div
	for (i=0; i<fobj->div_obj.num_factors; i++)
		zFree(&fobj->div_obj.factors[i]);
	free(fobj->div_obj.factors);

	// then free other stuff in div
	mpz_clear(fobj->div_obj.gmp_n);
	mpz_clear(fobj->div_obj.gmp_f);

	// free any factors found in nfs
	for (i=0; i<fobj->nfs_obj.num_factors; i++)
		zFree(&fobj->nfs_obj.factors[i]);
	free(fobj->nfs_obj.factors);

	// then free other stuff in nfs
	mpz_clear(fobj->nfs_obj.gmp_n);

	//free general fobj stuff
	zFree(&fobj->N);
	sFree(&fobj->str_N);

	clear_factor_list(fobj);
	free(fobj->fobj_factors);

	return;
}

void alloc_factobj(fact_obj_t *fobj)
{
	zInit(&fobj->N);
	sInit(&fobj->str_N);	

	fobj->rho_obj.num_poly = 3;
	fobj->rho_obj.polynomials = (uint32 *)malloc(fobj->rho_obj.num_poly * sizeof(uint32));
	fobj->rho_obj.polynomials[0] = 3;
	fobj->rho_obj.polynomials[1] = 2;
	fobj->rho_obj.polynomials[2] = 1;
	mpz_init(fobj->rho_obj.gmp_n);
	mpz_init(fobj->rho_obj.gmp_f);

	mpz_init(fobj->pm1_obj.gmp_n);
	mpz_init(fobj->pm1_obj.gmp_f);
	
	mpz_init(fobj->pp1_obj.gmp_n);
	mpz_init(fobj->pp1_obj.gmp_f);

	fobj->ecm_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->ecm_obj.gmp_n);
	mpz_init(fobj->ecm_obj.gmp_f);

	fobj->squfof_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->squfof_obj.gmp_n);
	mpz_init(fobj->squfof_obj.gmp_f);

	fobj->qs_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->qs_obj.gmp_n);

	fobj->div_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->div_obj.gmp_n);
	mpz_init(fobj->div_obj.gmp_f);

	fobj->nfs_obj.factors = (z *)malloc(sizeof(z));
	mpz_init(fobj->nfs_obj.gmp_n);

	fobj->allocated_factors = 256;
	fobj->fobj_factors = (factor_t *)malloc(256 * sizeof(factor_t));

	fobj->ecm_obj.num_factors = 0;	
	fobj->qs_obj.num_factors = 0;	
	fobj->div_obj.num_factors = 0;	
	fobj->nfs_obj.num_factors = 0;	
	fobj->num_factors = 0;

	return;
}

void reset_factobj(fact_obj_t *fobj)
{
	// keep all of the settings in fobj, but do an init/free cycle on all
	// allocated structures
	free_factobj(fobj);
	alloc_factobj(fobj);

	return;
}

void add_to_factor_list(fact_obj_t *fobj, mpz_t n)
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
		if (mpz_cmp(n, fobj->fobj_factors[i].factor) == 0)
		{
			found = 1;
			fobj->fobj_factors[i].count++;
			return;
		}
	}

	//else, put it in the list
	mpz_init(fobj->fobj_factors[fobj->num_factors].factor);
	mpz_set(fobj->fobj_factors[fobj->num_factors].factor, n);
	fobj->fobj_factors[fobj->num_factors].count = 1;
	if (mpz_probab_prime_p(n, NUM_WITNESSES))
	{
		if (mpz_cmp_ui(n, 100000000) <= 0)
			fobj->fobj_factors[fobj->num_factors].type = PRIME;
		else
			fobj->fobj_factors[fobj->num_factors].type = PRP;
	}
	else
		fobj->fobj_factors[fobj->num_factors].type = COMPOSITE;

	fobj->num_factors++;

	return;
}

void delete_from_factor_list(fact_obj_t *fobj, mpz_t n)
{
	//remove the number n from the global factor list
	uint32 i;

	//find the factor
	for (i=0;i<fobj->num_factors; i++)
	{
		if (mpz_cmp(n,fobj->fobj_factors[i].factor) == 0)
		{
			int j;
			// copy everything above this in the list back one position
			for (j=i; j<fobj->num_factors-1; j++)
			{
				mpz_set(fobj->fobj_factors[j].factor, fobj->fobj_factors[j+1].factor);
				fobj->fobj_factors[j].count = fobj->fobj_factors[j+1].count;
			}
			// remove the last one in the list
			fobj->fobj_factors[j].count = 0;
			mpz_clear(fobj->fobj_factors[j].factor);

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
		mpz_clear(fobj->fobj_factors[i].factor);
	}
	fobj->num_factors = 0;

	return;
}

void print_factors(fact_obj_t *fobj)
{
	uint32 i;
	int j;
	mpz_t tmp, tmp2;

	//always print factors unless complete silence is requested
	if (VFLAG >= 0)
	{
		printf("\n\n***factors found***\n\n");
		mpz_init(tmp);
		mpz_set_ui(tmp, 1);

		for (i=0; i<fobj->num_factors; i++)
		{
			if (fobj->fobj_factors[i].type == COMPOSITE)
			{
				for (j=0;j<fobj->fobj_factors[i].count;j++)
				{
					mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
					gmp_printf("C%d = %Zd\n",mpz_sizeinbase(fobj->fobj_factors[i].factor, 10),
						fobj->fobj_factors[i].factor);
				}
			}
			else if (fobj->fobj_factors[i].type == PRP)
			{
				for (j=0;j<fobj->fobj_factors[i].count;j++)
				{
					mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
					gmp_printf("PRP%d = %Zd\n",mpz_sizeinbase(fobj->fobj_factors[i].factor, 10),
						fobj->fobj_factors[i].factor);
				}
			}
			else if (fobj->fobj_factors[i].type == PRIME)
			{
				for (j=0;j<fobj->fobj_factors[i].count;j++)
				{
					mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
					gmp_printf("P%d = %Zd\n",mpz_sizeinbase(fobj->fobj_factors[i].factor, 10),
						fobj->fobj_factors[i].factor);
				}
			}
			else
			{
				//type not set, determine it now
				if (mpz_probab_prime_p(fobj->fobj_factors[i].factor, NUM_WITNESSES))
				{
					for (j=0;j<fobj->fobj_factors[i].count;j++)
					{
						mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
						gmp_printf("PRP%d = %Zd\n",mpz_sizeinbase(fobj->fobj_factors[i].factor, 10),
							fobj->fobj_factors[i].factor);
					}
				}
				else
				{
					for (j=0;j<fobj->fobj_factors[i].count;j++)
					{
						mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
						gmp_printf("C%d = %Zd\n",mpz_sizeinbase(fobj->fobj_factors[i].factor, 10),
							fobj->fobj_factors[i].factor);
					}
				}
			}
		}

		mpz_init(tmp2);
		mp2gmp(&fobj->N, tmp2);
		if (mpz_cmp_ui(tmp2, 1) > 0)
		{
			// non-trivial N remaining... compute and display the known co-factor
			mpz_tdiv_q(tmp2, tmp2, tmp);
			if (mpz_cmp_ui(tmp2, 1) > 0)
			{
				if (mpz_probab_prime_p(tmp2, NUM_WITNESSES))
				{
					gmp_printf("\n***co-factor***\nPRP%d = %Zd\n",
						mpz_sizeinbase(tmp2, 10), tmp2);
				}
				else
				{
					gmp_printf("\n***co-factor***\nC%d = %Zd\n",
						mpz_sizeinbase(tmp2, 10), tmp2);
				}
			}			
		}
		mpz_clear(tmp);
		mpz_clear(tmp2);
	}

	return;
}

double get_qs_time_estimate(fact_obj_t *fobj, double freq, int digits)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think siqs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;

	cpu = yafu_get_cpu_type();
	//if we have tuning info, use that instead
	if (fobj->qs_obj.qs_multiplier != 0 && fobj->qs_obj.qs_exponent != 0 && fobj->qs_obj.qs_tune_freq != 0)
	{
		if (VFLAG >= 2)
			printf("***** using tuning data for QS time estimation\n");
		estimate = fobj->qs_obj.qs_multiplier * exp(fobj->qs_obj.qs_exponent * digits);

		//scale with frequency
		estimate = estimate * fobj->qs_obj.qs_tune_freq / freq; 
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

double get_gnfs_time_estimate(fact_obj_t *fobj, double freq, int digits)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think gnfs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;

	cpu = yafu_get_cpu_type();
	//if we have tuning info, use that instead
	if (fobj->nfs_obj.gnfs_multiplier != 0 && fobj->nfs_obj.gnfs_exponent != 0 && fobj->nfs_obj.gnfs_tune_freq != 0)
	{
		if (VFLAG >= 2)
			printf("***** using tuning data for GNFS time estimation\n");
		estimate = fobj->nfs_obj.gnfs_multiplier * exp(fobj->nfs_obj.gnfs_exponent * digits);

		//scale with frequency
		estimate = estimate * fobj->nfs_obj.gnfs_tune_freq / freq; 
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
	mpz_t b, fact_obj_t *fobj)
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
			printf("div: primes less than %d\n",B1);
		fobj->prime_threshold = B1 * B1;
		mpz_set(fobj->div_obj.gmp_n, b);
		fobj->div_obj.print = 1;
		fobj->div_obj.limit = B1;
		zTrial(fobj);
		mpz_set(b, fobj->div_obj.gmp_n);
		break;

	case rho_work:
		mpz_set(fobj->rho_obj.gmp_n,b);
		brent_loop(fobj);
		mpz_set(b,fobj->rho_obj.gmp_n);
		break;

	case fermat_work:
		if (VFLAG >= 0)
			printf("fmt: %d iterations\n",B1);
		mpz_set(fobj->div_obj.gmp_n,b);
		zFermat(B1,fobj);
		mpz_set(b,fobj->div_obj.gmp_n);
		break;

	case ecm_curve:
		tmp1 = fobj->ecm_obj.B1;
		tmp2 = fobj->ecm_obj.B2;
		fobj->ecm_obj.B1 = B1;
		fobj->ecm_obj.B2 = B2;
		fobj->ecm_obj.num_curves = *work;
		mpz_set(fobj->ecm_obj.gmp_n, b);
		*work = ecm_loop(fobj);
		mpz_set(b, fobj->ecm_obj.gmp_n);
		fobj->ecm_obj.B1 = tmp1;
		fobj->ecm_obj.B2 = tmp2;
		break;

	case pp1_curve:
		tmp1 = fobj->pp1_obj.B1;
		tmp2 = fobj->pp1_obj.B2;
		fobj->pp1_obj.B1 = B1;
		fobj->pp1_obj.B2 = B2;
		mpz_set(fobj->pp1_obj.gmp_n,b);
		fobj->pp1_obj.numbases = *work;
		williams_loop(fobj);
		mpz_set(b,fobj->pp1_obj.gmp_n);
		fobj->pp1_obj.B1 = tmp1;
		fobj->pp1_obj.B2 = tmp2;
		break;

	case pm1_curve:
		tmp1 = fobj->pm1_obj.B1;
		tmp2 = fobj->pm1_obj.B2;
		fobj->pm1_obj.B1 = B1;
		fobj->pm1_obj.B2 = B2;
		mpz_set(fobj->pm1_obj.gmp_n,b);
		pollard_loop(fobj);
		mpz_set(b,fobj->pm1_obj.gmp_n);
		fobj->pm1_obj.B1 = tmp1;
		fobj->pm1_obj.B2 = tmp2;
		break;

	case qs_work:
		//gettimeofday(&start2,NULL);
		mpz_set(fobj->qs_obj.gmp_n,b);
		SIQS(fobj);
		mpz_set(b,fobj->qs_obj.gmp_n);
		break;
		//gettimeofday(&stop2,NULL);
		//difference = my_difftime (&start2, &stop2);
		//t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		//free(difference);

	case nfs_work:
		mpz_set(fobj->nfs_obj.gmp_n,b);
		nfs(fobj);
		mpz_set(b,fobj->nfs_obj.gmp_n);
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

int check_if_done(fact_obj_t *fobj, mpz_t N)
{
	int i, done = 0;
	mpz_t tmp;

	mpz_init(tmp);
	mpz_set_ui(tmp, 1);

	/* if the user only wants to find one factor, check for that here... */
	if (fobj->autofact_obj.want_only_1_factor && (fobj->num_factors >= 1))
	{
		done = 1;
		mpz_clear(tmp);
		return done;
	}

	// check if the number is completely factorized
	for (i=0; i<fobj->num_factors; i++)
	{		
		int j;
		for (j=0; j<fobj->fobj_factors[i].count; j++)
			mpz_mul(tmp, tmp, fobj->fobj_factors[i].factor);
	}

	if (mpz_cmp(N,tmp) == 0)
	{		
		// yes, they are equal.  make sure everything is prp or prime.
		done = 1;
		for (i=0; i<fobj->num_factors; i++)
		{
			if (!mpz_probab_prime_p(fobj->fobj_factors[i].factor, NUM_WITNESSES))
			{
				refactor_depth++;
				if (refactor_depth > 3)
				{
					printf("too many refactorization attempts, aborting\n");
					break;
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
					gmp2mp(fobj->fobj_factors[i].factor, &fobj_refactor->N);

					// recurse on factor
					factor(fobj_refactor);

					// remove the factor from the original list
					delete_from_factor_list(fobj, fobj->fobj_factors[i].factor);

					// add all factors found during the refactorization
					for (j=0; j< fobj_refactor->num_factors; j++)
					{
						int k;
						for (k=0; k < fobj_refactor->fobj_factors[j].count; k++)
							add_to_factor_list(fobj, fobj_refactor->fobj_factors[j].factor);
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

	mpz_clear(tmp);
	return done;
}

int switch_to_qs(fact_obj_t *fobj, mpz_t N, double *time_available, int force_switch)
{
	// compare the total time spent so far with the estimate of how long it 
	// would take to finish using qs and decide whether or not to switch over
	// to qs.
	int decision, sizeN, digitsN = mpz_sizeinbase(N, 10);
	double qs_est_time, nfs_est_time;
	
	// if the size of N is small enough, always switch to qs
	sizeN = mpz_sizeinbase(N, 2);
	if (sizeN < 135)
	{
		*time_available = 0;
		decision = 1;
	}
	else
	{
		qs_est_time = get_qs_time_estimate(fobj, MEAS_CPU_FREQUENCY, digitsN);
		if (VFLAG >= 2)
			printf("***** qs time estimate = %lg seconds\n",qs_est_time);

		nfs_est_time = get_gnfs_time_estimate(fobj, MEAS_CPU_FREQUENCY, digitsN);
		if (VFLAG >= 2)
			printf("***** gnfs time estimate = %lg seconds\n",nfs_est_time);

		//proceed with whichever estimate is smaller
		if (qs_est_time <= nfs_est_time)
		{		
			if (force_switch)
			{
				//calling code is forcing a decision to switch
				if ((fobj->nfs_obj.gnfs_exponent == 0) && (digitsN > 95))
				{
					//est time decision was to use qs, but the size is high enough and we're using
					//qs time only because nfs hasn't been tuned, so use nfs instead
					decision = 2;
					*time_available = 0;
				}
				else
				{
					decision = 1;
					*time_available = 0;
				}
			}
			else if (qs_est_time > 1000000000.0)
			{
				// if qs_est_time is very large, then we don't have a good estimate.  flag the caller 
				// of this fact
				decision = 0;
				*time_available = -1;
			}		
			else if (total_time > (fobj->autofact_obj.target_ecm_qs_ratio * qs_est_time))
			{
				// if the total time we've spent so far is greater than a fraction of the time
				// we estimate it would take QS to finish, switch to qs.  
				if ((fobj->nfs_obj.gnfs_exponent == 0) && (digitsN > 95))
				{
					//est time decision was to use qs, but the size is high enough and we're using
					//qs time only because nfs hasn't been tuned, so use nfs instead
					decision = 2;
					*time_available = 0;
				}
				else
				{
					decision = 1;
					*time_available = 0;
				}
			}
			else
			{
				// otherwise, return the amount of time we have left before the switchover.
				decision = 0;
				*time_available = (fobj->autofact_obj.target_ecm_qs_ratio * qs_est_time) - total_time;
			}
		}
		else
		{
			if (force_switch)
			{
				decision = 2;
				*time_available = 0;				
			}
			else if (nfs_est_time > 1e9)
			{
				decision = 0;
				*time_available = -1;
			}		
			else if (total_time > fobj->autofact_obj.target_ecm_gnfs_ratio * nfs_est_time)
			{
				// if the total time we've spent so far is greater than a fraction of the time
				// we estimate it would take gnfs to finish, switch to gnfs.  
				decision = 2;
				*time_available = 0;
			}
			else
			{
				// otherwise, return the amount of time we have left before the switchover.
				decision = 0;
				*time_available = (fobj->autofact_obj.target_ecm_gnfs_ratio * nfs_est_time) - total_time;
			}
		}

	}

	return decision;
}

enum factorization_state scale_requested_work(method_timing_t *method_times, 
	enum factorization_state fact_state, int *next_work, double time_available, mpz_t N)
{
	enum factorization_state new_state = fact_state;
	double base_time_per_curve, time_per_curve;
	int default_curves_15digit = 25;
	int default_curves_20digit = 90;
	int default_curves_25digit = 200;
	int default_curves_30digit = 400;
	int default_curves_35digit = 1000;
	char state_str[100];

	if (time_available < 0)
	{
		// this is a flag from the switch_to_qs function indicating that we
		// can't make a good estimate about how long qs will take.  thus we also
		// don't know how much time is left.  we have to make a decision about how
		// much and what kind of work to do based on the size of the input and our
		// current factorization state.
		int sizeN = mpz_sizeinbase(N, 2);

		switch (fact_state)
		{
		case state_trialdiv:			
			strcpy(state_str,"trial division");
			break;

		case state_fermat:
			strcpy(state_str,"Fermat");
			break;

		case state_rho:
			strcpy(state_str,"Rho");
			break;

		case state_pp1_lvl1:
			strcpy(state_str,"P+1");
			break;
		
		case state_pp1_lvl2:
			// equivalent to ecm @ 30 digits			
			if (sizeN < 300)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = 3;
				strcpy(state_str,"P+1");			
			}

			break;

		case state_pp1_lvl3:
			// equivalent to ecm @ 35 digits
			if (sizeN < 330)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = 3;
				strcpy(state_str,"P+1");		
			}

			break;

		case state_pm1_lvl1:
			if (sizeN < 160)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = 1;
				strcpy(state_str,"P-1");		
			}

			break;

		case state_pm1_lvl2:
			// equivalent to ecm @ 30 digits
			if (sizeN < 300)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = 1;
				strcpy(state_str,"P-1");		
			}

			break;

		case state_pm1_lvl3:
			// equivalent to ecm @ 35 digits
			if (sizeN < 330)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = 1;
				strcpy(state_str,"P-1");		
			}

			break;

		case state_ecm_15digit:
			if (sizeN < 180)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = default_curves_15digit;
				strcpy(state_str,"ECM");		
			}

			break;

		case state_ecm_20digit:
			if (sizeN < 220)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = default_curves_20digit;
				strcpy(state_str,"ECM");		
			}

			break;

		case state_ecm_25digit:
			if (sizeN < 260)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = default_curves_25digit;
				strcpy(state_str,"ECM");		
			}

			break;

		case state_ecm_30digit:
			if (sizeN < 300)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = default_curves_30digit;
				strcpy(state_str,"ECM");		
			}

			break;

		case state_ecm_35digit:
			if (sizeN < 330)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = default_curves_35digit;
				strcpy(state_str,"ECM");		
			}

			break;

		case state_ecm_auto_increasing:
			if (mpz_sizeinbase(N, 10) < 120)
			{
				new_state = state_qs;
				strcpy(state_str,"SIQS");
			}
			else
			{
				*next_work = auto_increasing_curves;
				strcpy(state_str,"ECM");		
			}

			break;

		default:
			printf("unknown factorization state\n");
			exit(-1);

		}

		if (VFLAG >= 2)
			printf("***** qs/nfs time estimate too high or not available, proceeding to %s\n",state_str);
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

	mpz_t b, origN, copyN;
	enum factorization_state fact_state = state_idle;
	method_timing_t method_times;
	int curves = 1;
	int decision = 0;
	int min_pretest_done = 0;
	int done = 0;	
	FILE *flog;
	struct timeval start, stop;
	double t_time;
	TIME_DIFF *	difference;
	int user_defined_ecm_b2 = fobj->ecm_obj.stg2_is_default;
	int user_defined_pp1_b2 = fobj->pp1_obj.stg2_is_default;
	int user_defined_pm1_b2 = fobj->pm1_obj.stg2_is_default;
	int force_switch = 0;
	FILE *data;
	char tmpstr[GSTR_MAXSIZE];

	//factor() always ignores user specified B2 values
	fobj->ecm_obj.stg2_is_default = 1;
	fobj->pp1_obj.stg2_is_default = 1;
	fobj->pm1_obj.stg2_is_default = 1;

	mpz_init(origN);
	mpz_init(copyN);
	mpz_init(b);

	mp2gmp(&fobj->N, origN);
	mpz_set(copyN, origN);
	mpz_set(b, origN);

	if (mpz_cmp_ui(b,1) <= 0)
	{
		mpz_clear(copyN);
		mpz_clear(origN);
		mpz_clear(b);
		return;
	}
	
	gettimeofday(&start, NULL);

	flog = fopen(fobj->flogname,"a");
	logprint(flog,"\n");
	logprint(flog,"****************************\n");
	logprint(flog,"Starting factorization of %s\n",mpz_get_str(gstr1.s, 10, b));
	logprint(flog,"****************************\n");
	fclose(flog);

	fobj->autofact_obj.autofact_active = 1;

	if (VFLAG >= 0)
	{
		gmp_printf("factoring %Zd\n",b);
		printf("using pretesting plan: %s\n\n",fobj->autofact_obj.plan_str);
	}

	// initialize time per curve
	method_times.ecm_autoinc_time_per_curve = 0;
	total_time = 0;

	//default choice
	fact_state = state_trialdiv;

	//check to see if a siqs savefile exists for this input	
	data = fopen(fobj->qs_obj.siqs_savefile,"r");

	if (data != NULL)
	{	
		char *substr;
		mpz_t tmpz;

		//read in the number from the savefile
		mpz_init(tmpz);
		fgets(tmpstr,1024,data);
		substr = tmpstr + 2;
		mpz_set_str(tmpz, substr, 0);	//auto detect the base

		// gcd required because the savefile may have a multiplier applied
		mpz_gcd(tmpz, tmpz, b);

		if (mpz_cmp(tmpz, b) == 0)
		{
			if (VFLAG > 1)
				printf("found siqs savefile, resuming siqs\n");

			//override default choice
			fact_state = state_qs;
		}
		mpz_clear(tmpz);
	}

	while (!done)
	{
		switch (fact_state)
		{
		case state_trialdiv:
			curves = 1;
			t_time = do_work(trialdiv_work, 10000, -1, &curves, b, fobj);
			method_times.trialdiv_time = t_time;
			total_time += t_time * curves;
			fact_state = state_fermat;
			break;

		case state_fermat:
			curves = 1;
			t_time = do_work(fermat_work, fobj->div_obj.fmtlimit, -1, &curves, b, fobj);
			method_times.fermat_time = t_time;
			total_time += t_time * curves;
			fact_state = state_rho;
			break;

		case state_rho:
			curves = 1;
			t_time = do_work(rho_work, -1, -1, &curves, b, fobj);
			method_times.rho_time = t_time;
			total_time += t_time * curves;
			
			//after trial division, fermat, and rho, we're ready to 
			//consider qs methods.  all pretest plans do at least this much work.
			min_pretest_done = 1;

			//where to go from here depends on the pretest plan in place
			if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_NONE)
				force_switch = 1;
			else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_DEEP)
				fact_state = state_ecm_25digit;
			else
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
			t_time = do_work(pp1_curve, 1250000, 125000000, &curves, b, fobj);
			method_times.pp1_lvl2_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_pm1_lvl2;
			break;

		case state_pp1_lvl3:
			t_time = do_work(pp1_curve, 5000000, 500000000, &curves, b, fobj);
			method_times.pp1_lvl3_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_pm1_lvl3;
			break;

		case state_pm1_lvl1:
			curves = 1;
			t_time = do_work(pm1_curve, 100000, 10000000, &curves, b, fobj);
			method_times.pm1_lvl1_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_15digit;

			//if we are not doing any ecm, force a switch to a sieve method
			if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_NOECM)
				force_switch = 1;
				
			break;

		case state_pm1_lvl2:
			t_time = do_work(pm1_curve, 2500000, 250000000, &curves, b, fobj);
			method_times.pm1_lvl2_time_per_curve = t_time;
			total_time += t_time * curves;
			fact_state = state_ecm_30digit;

			//if we are doing light pretesting, force a switch to a sieve method
			if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_LIGHT)
				force_switch = 1;

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

			//for deep pretesting, skip the pp1/pm1 and go right to 30 digit
			//ecm, after doing this curve at 25 digits to get timing info.
			if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_DEEP)
				fact_state = state_ecm_30digit;
			else
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

		default:
			printf("unknown factorization state\n");
			exit(-1);

		}		

		// first, check if we're done
		done = check_if_done(fobj, origN);	

		// paranoia
		if (mpz_cmp_ui(b, 0) == 0)
		{
			printf("b = 0, exiting\n");
			done = 1;
		}

		if ((!done && min_pretest_done) || force_switch)
		{
			// if we're not done, decide what to do next: either the default next state or
			// switch to qs.
			decision = switch_to_qs(fobj, b, &t_time, force_switch);
			//printf("total time = %f\n",total_time);
			//printf("time available = %f\n",t_time);
			if (decision)
			{
				// either the number is small enough to finish off, or it would be better, 
				// timewise, to switch
				if (fobj->autofact_obj.only_pretest)
				{
					// we're ready to go to a sieve method, but we only want to 
					// pretest.  so we bail.
					done = 1;
				}
				else
				{
					if (decision == 2)
						fact_state = state_nfs;
					else
						fact_state = state_qs;
				}
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

	if (fobj->num_factors >= 1) 
	{
		//If the only factor in our array == N, then N is prime or prp...
		if (fobj->autofact_obj.want_output_primes && (mpz_cmp(fobj->fobj_factors[0].factor,origN) == 0))
		{
			if ((fobj->autofact_obj.op_file = fopen(fobj->autofact_obj.op_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.op_str);
			else
			{
				if (fobj->autofact_obj.want_output_expressions)
					fprintf(fobj->autofact_obj.op_file, "%s\n", fobj->str_N.s);
				else
					gmp_fprintf(fobj->autofact_obj.op_file, "%Zd\n", origN);
				if (fclose(fobj->autofact_obj.op_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.op_str);
			}
		}

		//If the first factor in the array != N, then is composite and we have factors...
		if (fobj->autofact_obj.want_output_factors && (mpz_cmp(fobj->fobj_factors[0].factor,origN) != 0))
		{
			if ((fobj->autofact_obj.of_file = fopen(fobj->autofact_obj.of_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.of_str);
			else
			{
				int i;
				//fprintf(fobj->autofact_obj.of_file, "%s\n", z2decstr(&origN,&gstr1));
				if (fobj->autofact_obj.want_output_expressions)
					fprintf(fobj->autofact_obj.of_file, "(%s)", fobj->str_N.s);
				else
					gmp_fprintf(fobj->autofact_obj.of_file, "%Zd", origN);
				for (i=0; i<fobj->num_factors; i++)
				{
					gmp_fprintf(fobj->autofact_obj.of_file, "/%Zd", fobj->fobj_factors[i].factor);
					if (fobj->fobj_factors[i].count > 1)
						fprintf(fobj->autofact_obj.of_file, "^%d", fobj->fobj_factors[i].count);
					//fprintf(fobj->autofact_obj.of_file, "\n");
				}
				fprintf(fobj->autofact_obj.of_file,"\n");
				if (fclose(fobj->autofact_obj.of_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.of_str);
			}
		}
	}
	else //assume: composite with no known factors... (need to clarify)
	{
		if (fobj->autofact_obj.want_output_unfactored)
		{
			if ((fobj->autofact_obj.ou_file = fopen(fobj->autofact_obj.ou_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.ou_str);
			else
			{
				if (fobj->autofact_obj.want_output_expressions)
					fprintf(fobj->autofact_obj.ou_file, "%s\n", fobj->str_N.s);
				else
					gmp_fprintf(fobj->autofact_obj.ou_file, "%s\n", origN);
				if (fclose(fobj->autofact_obj.ou_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.ou_str);
			}
		}
	}

		
	if (mpz_cmp_ui(b, 1) != 0)
	{
		if (mpz_probab_prime_p(b, NUM_WITNESSES))
		{
			flog = fopen(fobj->flogname,"a");
			logprint(flog,"prp%d cofactor = %s\n",mpz_sizeinbase(b, 10),
				mpz_get_str(gstr1.s, 10, b));
			fclose(flog);
			add_to_factor_list(fobj,b);
		}
		else
		{
			flog = fopen(fobj->flogname,"a");
			logprint(flog,"c%d cofactor = %s\n",mpz_sizeinbase(b, 10),
				mpz_get_str(gstr1.s, 10, b));
			fclose(flog);
			add_to_factor_list(fobj,b);
		}
	}
	gmp2mp(b,&fobj->N);

	gettimeofday (&stop, NULL);
	difference = my_difftime (&start, &stop);
	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (VFLAG >= 0)
		printf("Total factoring time = %6.4f seconds\n",t_time);

	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
		printf("Could not open %s for appending\n",fobj->flogname);
	else
	{
		logprint(flog,"Total factoring time = %6.4f seconds\n",t_time);
		fclose(flog);
	}

	fobj->autofact_obj.autofact_active=0;

	//restore flags
	fobj->ecm_obj.stg2_is_default = user_defined_ecm_b2;
	fobj->pp1_obj.stg2_is_default = user_defined_pp1_b2;
	fobj->pm1_obj.stg2_is_default = user_defined_pm1_b2;

	mpz_clear(origN);
	mpz_clear(copyN);
	mpz_clear(b);
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


