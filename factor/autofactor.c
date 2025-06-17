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

#include <stdio.h>
#include "yafu.h"
#include "soe.h"
#include "factor.h"
#include "qs.h"
#include "nfs.h"
#include "nfs_impl.h"
#include "yafu_ecm.h"
#include "ytools.h"
#include "mpz_aprcl.h"
#include "cmdOptions.h"
#include <stdint.h>
#include "autofactor.h"
#include "gmp.h"
#include "ecm.h"
#include "microecm.h"

#include <math.h>

/* produced using ecm -v -v -v for the various B1 bounds (default B2).
/	Thanks A. Schindel !
/
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */
#define NUM_ECM_LEVELS 12
static int ecm_levels[NUM_ECM_LEVELS] = {
	15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 };
double ecm_data[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
	/*t15, 2000,	*/	{30,		12,			7,			5,			3,			2,			2,			2,			2,		1,		1},
	/*t20, 11000,	*/	{844,		74,			21,			8,			5,			3,			2,			2,			2,		2,		1},
	/*t25, 50000,	*/	{58129,		1539,		214,		50,			20,			11,			7,			5,			4,		3,		3},
	/*t30, 250000,	*/	{6711967,	49962,		3288,		430,		118,		54,			26,			14,			10,		8,		6},
	/*t35, 1E+06,	*/	{1.20E+09,	2292278,	68422,		4914,		904,		322,		122,		54,			34,		23,		15},
	/*t40, 3E+06,	*/	{2.90E+12,	1.40E+08,	1849287,	70293,		8613,		2350,		681,		242,		135,	82,		47},
	/*t45, 11E+06,	*/	{9.00E+99,	1.10E+10,	6.10E+07,	1214949,	97057,		20265,		4480,		1263,		613,	333,	168},
	/*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	2.50E+09,	2.50E+07,	1270662,	199745,		33652,		7404,		3133,	1512,	661},
	/*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	1.30E+11,	5.90E+08,	1.90E+07,	2246256,	283939,		48714,		17769,	7643,	2865},
	/*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	5.80E+16,	1.60E+10,	3.20E+08,	2.80E+07,	2655154,	350439,		111196,	42017,	13611},
	/*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	8.20E+21,	2.70E+13,	6.10E+09,	4.00E+08,	2.70E+07,	2768535,	751771,	250214,	69408} };

/* produced using ecm -v -v -v for the various B1 bounds using param 0 (default B2).
/ version 7.0.5
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */
double ecm_data_param0[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
	/*t15, 2000,	*/	{34,		7,			3,			2,			2,			1,			1,			1,			1,		1,		1},
	/*t20, 11000,	*/	{1282,		86,		    21,			8,			5,			3,			2,			2,			2,		2,		1},
	/*t25, 50000,	*/	{95034,	    1864,		214,		50,			20,			11,			7,			5,			4,		3,		3},
	/*t30, 250000,	*/	{1.2e+07,	62295,		3283,		430,		118,		54,			26,			14,			10,		8,		6},
	/*t35, 1E+06,	*/	{2.20E+09,	2924742,	69076,		4911,		910,		324,		122,		55,			34,		23,		15},
	/*t40, 3E+06,	*/	{9.00E+99,	1.80E+08,	1847472,	70940,		8615,		2351,		686,		246,		135,	82,		47},
	/*t45, 11E+06,	*/	{9.00E+99,	1.50E+10,	6.20E+07,	1226976,	97096,		20272,		4482,		1286,		614,	335,	168},
	/*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	2.50E+09,	2.50E+07,	1281819,	201449,		33676,		7557,		3135,	1521,	661},
	/*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	1.30E+11,	5.80E+08,	1.90E+07,	2247436,	284176,		49831,		17884,	7650,	2867},
	/*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	5.90E+16,	1.60E+10,	3.10E+08,	2.80E+07,	2657998,	361851,		111314,	42057,	13623},
	/*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	8.40E+21,	2.70E+13,	6.20E+09,	3.90E+08,	2.70E+07,	2844041,	752662,	250476,	69471} };



/* produced using ecm -v -v -v for the various B1 bounds using param 1 (default B2).
/ version 7.0.5
/					2k			11k			50k			250k		1M			3M			11M			43M			110M	260M	850M */
double ecm_data_param1[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
	/*t15, 2000,	*/	{43,		8,			4,			2,			2,			1,			1,			1,			1,		1,		1},
	/*t20, 11000,	*/	{1743,		107,		24,			9,			5,			4,			3,			2,			2,		2,		1},
	/*t25, 50000,	*/	{134769,	2402,		261,		58,			22,			13,			8,			5,			4,		3,		3},
	/*t30, 250000,	*/	{1.7e+07,	82576,		4108,		513,		137,		61,			29,			16,			11,		8,		6},
	/*t35, 1E+06,	*/	{3.30E+09,	3967858,	88265,		6022,		1071,		374,		138,		61,			38,		25,		17},
	/*t40, 3E+06,	*/	{9.00E+99,	2.50E+08,	2402639,	87544,		10283,		2753,		788,		278,		151,	91,		52},
	/*t45, 11E+06,	*/	{9.00E+99,	2.0E+10,	8.10E+07,	1534319,	118226,		24017,		5208,		1459,		692,	373,	185},
	/*t50, 44E+06,	*/	{9.00E+99,	9.00E+99,	3.30E+09,	3.10E+07,	1565171,	241048,		39497,		8704,		3583,	1709,	737},
	/*t55, 110E+06,	*/	{9.00E+99,	9.00E+99,	2.30E+11,	7.50E+08,	2.40E+07,	2713723,	336066,		57844,		20479,	8656,	3220},
	/*t60, 260E+06,	*/	{9.00E+99,	9.00E+99,	1.50E+17,	2.00E+10,	3.90E+08,	3.40E+07,	3167410,	419970,		128305,	47888,	15391},
	/*t65, 850E+06,	*/	{9.00E+99,	9.00E+99,	2.10E+22,	6.70E+13,	7.70E+09,	4.80E+08,	3.20E+07,	3346252,	872747,	288516,	78923} };



/* produced using ecm -v -power 1 -param 0 for the various B1 bounds (B2=B1*100).
/  Thanks tbusby!
/
/                          2k       11k      50k      250k     1M       3M       11M      43M      110M     260M     850M     2900M */
double avx_ecm_data[NUM_ECM_LEVELS][NUM_ECM_LEVELS] = {
	/*t15, 2000,    */    {32,      8,       4,       3,       2,       2,       2,       1,       1,       1,       1,       1     },
	/*t20, 11000,   */    {1203,    98,      26,      11,      7,       5,       4,       3,       3,       2,       2,       2     },
	/*t25, 50000,   */    {88933,   2174,    281,     74,      32,      20,      12,      8,       7,       5,       4,       4     },
	/*t30, 250000,  */    {1.1E+07, 73653,   4455,    671,     206,     100,     50,      28,      20,      15,      11,      8     },
	/*t35, 1E+06,   */    {2.0E+09, 3495192, 95330,   7965,    1676,    635,     250,     115,     73,      51,      33,      23    },
	/*t40, 3E+06,   */    {9.0E+99, 2.2E+08, 2610023, 117616,  16589,   4858,    1492,    550,     308,     194,     111,     68    },
	/*t45, 11E+06,  */    {9.0E+99, 1.8E+10, 8.8e+07, 2088438, 194098,  43492,   10273,   3014,    1481,    835,     419,     228   },
	/*t50, 44E+06,  */    {9.0E+99, 9.0E+99, 3.7E+09, 4.4E+07, 2626303, 445887,  80191,   18579,   7942,    3995,    1748,    838   },
	/*t55, 110E+06, */    {9.0E+99, 9.0E+99, 2.7E+11, 1.0E+09, 4.1E+07, 5150239, 699132,  126904,  46946,   20992,   7961,    3352  },
	/*t60, 260E+06, */    {9.0E+99, 9.0E+99, 9.0E+99, 2.9E+10, 6.7E+08, 6.6e+07, 6727121, 949988,  302925,  119976,  39223,   14455 },
	/*t65, 850E+06, */    {9.0E+99, 9.0E+99, 9.0E+99, 9.0E+99, 1.4E+10, 9.7e+08, 7.1E+07, 7721832, 2115419, 739454,  207648,  66688 } };
/*t70, 260E+06,     {9.0E+99, 9.0E+99, 9.0E+99, 9.0E+99, 9.0E+99, 1.5e+10, 7.6E+08, 6.7E+07, 1.6E+07, 4880638, 1173260, 327240} */



typedef struct
{	
	// total effort so far
	double total_time;
	double qs_time;
	double nfs_time;
	double trialdiv_time;
	double fermat_time;
	double rho_time;
	double pp1_time;
	double pm1_time;
	double ecm_time;
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
	double ecm_40digit_time_per_curve;
	double ecm_45digit_time_per_curve;
	double ecm_50digit_time_per_curve;
	double ecm_55digit_time_per_curve;
	double ecm_60digit_time_per_curve;
	double ecm_65digit_time_per_curve;
    double initial_work;

	// amount of work we've done in various areas
	uint32_t  tdiv_limit;
	uint32_t  fermat_iterations;
	uint32_t  rho_iterations;
	uint32_t  rho_bases;
	uint32_t  pp1_lvl1_curves;
	uint32_t  pm1_lvl1_curves;
	uint32_t  pp1_lvl2_curves;
	uint32_t  pm1_lvl2_curves;
	uint32_t  pp1_lvl3_curves;
	uint32_t  pm1_lvl3_curves;
	uint32_t  ecm_15digit_curves;
	uint32_t  ecm_20digit_curves;
	uint32_t  ecm_25digit_curves;
	uint32_t  ecm_30digit_curves;
	uint32_t  ecm_35digit_curves;
	uint32_t  ecm_40digit_curves;
	uint32_t  ecm_45digit_curves;
	uint32_t  ecm_50digit_curves;
	uint32_t  ecm_55digit_curves;
	uint32_t  ecm_60digit_curves;
	uint32_t  ecm_65digit_curves;
	int min_pretest_done;

	double tlevels[NUM_ECM_LEVELS];

	// max amount of work we'll allow in various areas.
	// to be filled in during init, or overriden by user
	uint32_t  tdiv_max_limit;
	uint32_t  fermat_max_iterations;
	uint32_t  rho_max_iterations;
	uint32_t  rho_max_bases;
	uint32_t  pp1_max_lvl1_curves;
	uint32_t  pm1_max_lvl1_curves;
	uint32_t  pp1_max_lvl2_curves;
	uint32_t  pm1_max_lvl2_curves;
	uint32_t  pp1_max_lvl3_curves;
	uint32_t  pm1_max_lvl3_curves;
	uint32_t  ecm_max_15digit_curves;
	uint32_t  ecm_max_20digit_curves;
	uint32_t  ecm_max_25digit_curves;
	uint32_t  ecm_max_30digit_curves;
	uint32_t  ecm_max_35digit_curves;
	uint32_t  ecm_max_40digit_curves;
	uint32_t  ecm_max_45digit_curves;
	uint32_t  ecm_max_50digit_curves;	
	uint32_t  ecm_max_55digit_curves;	
	uint32_t  ecm_max_60digit_curves;	
	uint32_t  ecm_max_65digit_curves;	

	// current parameters
	uint32_t  B1;
	uint64_t  B2;
	uint32_t  curves;

} factor_work_t;

enum factorization_state {
	state_idle,
	state_trialdiv,
	state_fermat,
	state_rho,
	state_pp1_lvl1,
	state_pm1_lvl1,
	state_pp1_lvl2,
	state_pm1_lvl2,
	state_pp1_lvl3,
	state_pm1_lvl3,
	state_ecm_15digit,
	state_ecm_20digit,
	state_ecm_25digit,
	state_ecm_30digit,
	state_ecm_35digit,
	state_ecm_40digit,
	state_ecm_45digit,
	state_ecm_50digit,
	state_ecm_55digit,
	state_ecm_60digit,
	state_ecm_65digit,
	state_qs,
	state_nfs,
	state_done
};

// local functions to do state based factorization
double get_qs_time_estimate(fact_obj_t *fobj, mpz_t b);
double get_gnfs_time_estimate(fact_obj_t *fobj, mpz_t b);
void do_work(enum factorization_state method, factor_work_t *fwork, mpz_t b, fact_obj_t *fobj);
enum factorization_state schedule_work(factor_work_t *fwork, mpz_t b, fact_obj_t *fobj);
int check_if_done(fact_obj_t *fobj, factor_work_t* fwork, mpz_t N);
uint32_t  get_ecm_curves_done(factor_work_t *fwork, enum factorization_state state);
uint32_t  set_ecm_curves_done(factor_work_t *fwork, enum factorization_state state, uint32_t  curves_done);
uint32_t  get_max_ecm_curves(factor_work_t *fwork, enum factorization_state state);
void set_work_params(factor_work_t *fwork, enum factorization_state state);
int check_tune_params(fact_obj_t *fobj);
enum factorization_state get_next_state(factor_work_t *fwork, fact_obj_t *fobj);
double compute_ecm_work_done(factor_work_t *fwork, int disp, FILE *log, int VFLAG, int LOGFLAG);
void init_factor_work(factor_work_t *fwork, fact_obj_t *fobj);
void interp_and_set_curves(factor_work_t *fwork, fact_obj_t *fobj, 
	enum factorization_state state, double work_done,
	double target_digits, int log_results);
void write_factor_json(fact_obj_t* fobj, factor_work_t* fwork,
	struct timeval* start, struct timeval* stop);

double get_qs_time_estimate(fact_obj_t *fobj, mpz_t b)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think siqs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;
	double freq = fobj->MEAS_CPU_FREQUENCY;
	int digits = gmp_base10(b);

	cpu = ytools_get_cpu_type();

	estimate = fobj->qs_obj.qs_multiplier * exp(fobj->qs_obj.qs_exponent * digits);
	estimate = estimate * fobj->qs_obj.qs_tune_freq / freq; 	

	//adjust for multi-threaded qs
	//if we assume threading is perfect, we'll get a smaller estimate for
	//qs than we can really achieve, resulting in less ECM, so fudge it a bit
	if (fobj->THREADS > 1)
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
			estimate = estimate / ((double)fobj->THREADS * 0.75);
			break;
		case 9:
		case 10:
			estimate = estimate / ((double)fobj->THREADS * 0.90);
			break;

		default:
			estimate = estimate / ((double)fobj->THREADS * 0.75);
			break;
		}
	}

	if (fobj->VFLAG >= 2)
		printf("fac: QS time estimation from tune data = %1.2f sec\n", estimate);

	return estimate;
}

double get_gnfs_time_estimate(fact_obj_t *fobj, mpz_t b)
{
	//using rough empirical scaling equations, number size, information
	//on cpu type, architecture, speed, and compilation options, 
	//compute how long we think gnfs would take to finish a factorization
	enum cpu_type cpu;
	double estimate;
	double freq = fobj->MEAS_CPU_FREQUENCY;
	int digits = gmp_base10(b);

	cpu = ytools_get_cpu_type();
	
	estimate = fobj->nfs_obj.gnfs_multiplier * exp(fobj->nfs_obj.gnfs_exponent * digits);
	estimate = estimate * fobj->nfs_obj.gnfs_tune_freq / freq; 

	//adjust for multi-threaded nfs
	//if we assume threading is perfect, we'll get a smaller estimate for
	//nfs than we can really achieve, resulting in less ECM, so fudge it a bit
	if (fobj->THREADS > 1)
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
			estimate = estimate / ((double)fobj->THREADS * 0.75);
			break;
		case 9:
		case 10:
			estimate = estimate / ((double)fobj->THREADS * 0.90);
			break;

		default:
			estimate = estimate / ((double)fobj->THREADS * 0.75);
			break;
		}
	}

	if (fobj->VFLAG >= 2)
		printf("fac: GNFS time estimation from tune data = %1.2f sec\n", estimate);

	return estimate;
}

void do_work(enum factorization_state method, factor_work_t *fwork, 
	mpz_t b, fact_obj_t *fobj)
{
	uint32_t  tmp1;
	uint64_t  tmp2;	
	struct timeval tstart, tstop;
	double t_time;
	uint32_t  curves_done;

	if (0) //mpz_cmp_ui(b, 1) <= 0)
	{
		gmp_printf("asked to do work on input = %Zd\n", b);
		printf("here are the factors I know about:\n");
		print_factors(fobj, fobj->factors, fobj->N, fobj->VFLAG, fobj->NUM_WITNESSES, fobj->OBASE);
		gmp_printf("here was the original input: %Zd\n", fobj->N);
		printf("please report this bug\n");
		exit(1);
	}

	gettimeofday(&tstart, NULL);

	switch (method)
	{
	case state_trialdiv:

        // if larger than a small bound, do a perfect power check
        fobj->prime_threshold = fwork->tdiv_max_limit * fwork->tdiv_max_limit;
        
        if ((mpz_cmp_ui(b, fobj->prime_threshold) > 1) && 
            mpz_perfect_power_p(b))
        {
            if (fobj->VFLAG > 0)
            {
                printf("fac: input is a perfect power\n");
            }

            logprint_oc(fobj->flogname, "a", "input is a perfect power\n");
            factor_perfect_power(fobj, b);

            mpz_set(fobj->N, b);
            break;
        }

        // then do all of the tdiv work requested
        if (fobj->VFLAG >= 0)
            printf("div: primes less than %d\n", fwork->tdiv_max_limit);
        
        mpz_set(fobj->div_obj.gmp_n, b);
        fobj->div_obj.print = 1;
        fobj->div_obj.limit = fwork->tdiv_max_limit;
        zTrial(fobj);
        mpz_set(b, fobj->div_obj.gmp_n);

        // record the work done
        fwork->tdiv_limit = fwork->tdiv_max_limit;

        // measure time for this completed work
        gettimeofday(&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

        fwork->trialdiv_time = t_time;
        fwork->total_time += t_time;

		break;

	case state_rho:

		if (mpz_perfect_square_p(b))
		{
			if (fobj->VFLAG > 0)
				printf("fac: input is a perfect square\n");

			mpz_sqrt(b, b);

			add_to_factor_list(fobj->factors, b,
				fobj->VFLAG, fobj->NUM_WITNESSES);

			add_to_factor_list(fobj->factors, b,
				fobj->VFLAG, fobj->NUM_WITNESSES);

			mpz_set_ui(b, 1);

			// measure time for this completed work
			gettimeofday(&tstop, NULL);
			t_time = ytools_difftime(&tstart, &tstop);

			fwork->rho_time = t_time;
			fwork->total_time += t_time;
			break;
		}

		if (mpz_perfect_power_p(b))
		{
			if (fobj->VFLAG > 0)
				printf("fac: input is a perfect power\n");

			factor_perfect_power(fobj, b);

			// measure time for this completed work
			gettimeofday(&tstop, NULL);
			t_time = ytools_difftime(&tstart, &tstop);

			fwork->rho_time = t_time;
			fwork->total_time += t_time;
			break;
		}

		// do all of the rho work requested
		mpz_set(fobj->rho_obj.gmp_n,b);
		brent_loop(fobj);
		mpz_set(b,fobj->rho_obj.gmp_n);

		// record the work done
		fwork->rho_bases = fwork->rho_max_bases;
		fwork->rho_iterations = fwork->rho_max_iterations;

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

		fwork->rho_time = t_time;
		fwork->total_time += t_time;
		break;

	case state_fermat:
		// do all of the fermat work requested
		if (fobj->VFLAG >= 0)
			printf("fmt: %d iterations\n", fwork->fermat_max_iterations);
		mpz_set(fobj->div_obj.gmp_n,b);
		zFermat(fwork->fermat_max_iterations, 1, fobj);
		mpz_set(b,fobj->div_obj.gmp_n);

		// record the work done
		fwork->fermat_iterations = fwork->fermat_max_iterations;

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

		fwork->fermat_time = t_time;
		fwork->total_time += t_time;

		break;

	case state_ecm_15digit:
	case state_ecm_20digit:
	case state_ecm_25digit:
	case state_ecm_30digit:
	case state_ecm_35digit:
	case state_ecm_40digit:
	case state_ecm_45digit:
	case state_ecm_50digit:
	case state_ecm_55digit:
	case state_ecm_60digit:
	case state_ecm_65digit:
		tmp1 = fobj->ecm_obj.B1;
		tmp2 = fobj->ecm_obj.B2;
		fobj->ecm_obj.B1 = fwork->B1;
		fobj->ecm_obj.B2 = fwork->B2;
		fobj->ecm_obj.num_curves = fwork->curves;
		mpz_set(fobj->ecm_obj.gmp_n, b);
		curves_done = ecm_loop(fobj);
		mpz_set(b, fobj->ecm_obj.gmp_n);
		fobj->ecm_obj.B1 = tmp1;
		fobj->ecm_obj.B2 = tmp2;

		// record the work done
		set_ecm_curves_done(fwork, method, 
			get_ecm_curves_done(fwork, method) + curves_done);
		
		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

		fwork->ecm_time += t_time;
		fwork->total_time += t_time;
		break;

	case state_pp1_lvl1:
	case state_pp1_lvl2:
	case state_pp1_lvl3:
		tmp1 = fobj->pp1_obj.B1;
		tmp2 = fobj->pp1_obj.B2;
		fobj->pp1_obj.B1 = fwork->B1;
		fobj->pp1_obj.B2 = fwork->B2;
		mpz_set(fobj->pp1_obj.gmp_n,b);
		fobj->pp1_obj.numbases = fwork->curves;
		williams_loop(fobj);
		mpz_set(b,fobj->pp1_obj.gmp_n);
		fobj->pp1_obj.B1 = tmp1;
		fobj->pp1_obj.B2 = tmp2;

		// record the work done
		if (method == state_pp1_lvl1)
			fwork->pp1_lvl1_curves = fwork->curves;
		else if (method == state_pp1_lvl2)
			fwork->pp1_lvl2_curves = fwork->curves;
		else if (method == state_pp1_lvl3)
			fwork->pp1_lvl3_curves = fwork->curves;

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

		fwork->pp1_time += t_time;
		fwork->total_time += t_time;
		break;

	case state_pm1_lvl1:
	case state_pm1_lvl2:
	case state_pm1_lvl3:
		tmp1 = fobj->pm1_obj.B1;
		tmp2 = fobj->pm1_obj.B2;
		fobj->pm1_obj.B1 = fwork->B1;
		fobj->pm1_obj.B2 = fwork->B2;
		mpz_set(fobj->pm1_obj.gmp_n,b);
		pollard_loop(fobj);
		mpz_set(b,fobj->pm1_obj.gmp_n);
		fobj->pm1_obj.B1 = tmp1;
		fobj->pm1_obj.B2 = tmp2;

		// record the work done
		if (method == state_pm1_lvl1)
			fwork->pm1_lvl1_curves = fwork->curves;
		else if (method == state_pm1_lvl2)
			fwork->pm1_lvl2_curves = fwork->curves;
		else if (method == state_pm1_lvl3)
			fwork->pm1_lvl3_curves = fwork->curves;
		
		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

		fwork->pm1_time += t_time;
		fwork->total_time += t_time;
		break;

	case state_qs:
		mpz_set(fobj->qs_obj.gmp_n,b);
		SIQS(fobj);
		mpz_set(b,fobj->qs_obj.gmp_n);

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

		fwork->qs_time = t_time;
		if (fobj->VFLAG > 0)
			printf("pretesting / qs ratio was %1.2f\n", 
				fwork->total_time / t_time); 
		break;

	case state_nfs:
		mpz_set(fobj->nfs_obj.gmp_n,b);
		nfs(fobj);
		mpz_set(b,fobj->nfs_obj.gmp_n);

		// measure time for this completed work
		gettimeofday (&tstop, NULL);
        t_time = ytools_difftime(&tstart, &tstop);

		fwork->nfs_time = t_time;
		if (fobj->VFLAG > 0)
			printf("pretesting / nfs ratio was %1.2f\n", 
				fwork->total_time / t_time); 

		break;

	default:
		printf("nothing to do for method %d\n", method);
		break;
	}

	return;
}

int check_if_done(fact_obj_t *fobj, factor_work_t* fwork, mpz_t N)
{
	int i, done = 0;
	mpz_t tmp;

	mpz_init(tmp);
	mpz_set_ui(tmp, 1);

	/* if the user only wants to find one factor, check for that here... */
	if (fobj->autofact_obj.want_only_1_factor && (fobj->factors->num_factors >= 1))
	{
		done = 1;
		mpz_clear(tmp);
		return done;
	}

	// more generally, stop after finding k factors
	if ((fobj->autofact_obj.stopk > 0) &&
		(fobj->factors->total_factors >= fobj->autofact_obj.stopk))
	{
		done = 1;
		mpz_clear(tmp);
		return done;
	}

	// check if the number is completely factorized
	for (i=0; i<fobj->factors->num_factors; i++)
	{		
		int j;
		for (j=0; j<fobj->factors->factors[i].count; j++)
			mpz_mul(tmp, tmp, fobj->factors->factors[i].factor);
	}

	if (mpz_cmp(N,tmp) == 0)
	{		
		// yes, they are equal.  make sure everything is prp or prime.
		done = 0;
		while (!done)
		{
			done = 1;
			for (i=0; i<fobj->factors->num_factors; i++)
			{
				if (is_mpz_prp(fobj->factors->factors[i].factor, fobj->NUM_WITNESSES) == 0)
				{	
					// can still do pretesting on composite factors, for instance
					// if primitive factor detection found a large composite factor.
					// 
					if (fobj->autofact_obj.only_pretest > 1)
					{
						if (fobj->autofact_obj.ecm_total_work_performed >= fobj->autofact_obj.only_pretest)
						{
							printf("fac: completed work %1.2f > %d, composite refactorization skipped\n",
								fobj->autofact_obj.ecm_total_work_performed, fobj->autofact_obj.only_pretest);
							// if we are pretesting and have already done all of
							// the ecm work specified then we are done.
							done = 1;
							break;
						}
						else
						{
							printf("fac: completed work %1.2f of %d, attempting composite refactorization\n",
								fobj->autofact_obj.ecm_total_work_performed, fobj->autofact_obj.only_pretest);
						}
					}
					
					if (fobj->refactor_depth > 3)
					{
						printf("too many refactorization attempts, aborting\n");
						done = 1;
						break;
					}
					else
					{
						fact_obj_t *fobj_refactor;
						int j;

						if (fobj->VFLAG > 0)
							printf("\nComposite result found, starting re-factorization\n");

						//gmp_printf("current factorization input: N = %Zd\n", fobj->N);
						//printf("here are the current factors of N I know about: \n");
						//print_factors(fobj, fobj->factors, fobj->N, 1, 1, 10);

						// load the new fobj with this number
						fobj_refactor = (fact_obj_t *)malloc(sizeof(fact_obj_t));
						init_factobj(fobj_refactor);
                        copy_factobj(fobj_refactor, fobj);

						mpz_set(fobj_refactor->N, fobj->factors->factors[i].factor);
                        fobj_refactor->refactor_depth = fobj->refactor_depth + 1;
						fobj_refactor->autofact_obj.initial_work = 
							fobj->autofact_obj.ecm_total_work_performed;
						
						// these will get recombined in the original fobj; 
						// no need to output them twice.
						fobj_refactor->autofact_obj.want_output_factors = 0;
						fobj_refactor->autofact_obj.want_output_primes = 0;
						fobj_refactor->autofact_obj.want_output_unfactored = 0;
						factor(fobj_refactor);

						// remember the ecm work we performed
						fobj->autofact_obj.ecm_total_work_performed =
							fobj->autofact_obj.initial_work =
							fobj_refactor->autofact_obj.ecm_total_work_performed;

						// original count: if > 1, need to add multiples of
						// each factor found.
						int ocount = fobj->factors->factors[i].count;

						// remove the factor from the original list
						delete_from_factor_list(fobj->factors, fobj->factors->factors[i].factor);

						// add all factors found during the refactorization
						for (j=0; j< fobj_refactor->factors->num_factors; j++)
						{
							int k;
							for (k=0; k < fobj_refactor->factors->factors[j].count * ocount; k++)
								add_to_factor_list(fobj->factors, 
                                    fobj_refactor->factors->factors[j].factor,
                                    fobj->VFLAG, fobj->NUM_WITNESSES);
						}

						// free temps
						free_factobj(fobj_refactor);
						free(fobj_refactor);

						// check again, since this factorization could have added new
						// composite factors
						done = 0;
					}
				}
			}
		}
	}

	mpz_clear(tmp);
	return done;
}

enum factorization_state get_next_state(factor_work_t *fwork, fact_obj_t *fobj)
{
	enum factorization_state next_state;

	// check each state's completed work against the maximum.
	// return the first one not complete.
	if (fwork->tdiv_limit < fwork->tdiv_max_limit)
		return state_trialdiv;	// always do this state if not done
	else if (fwork->fermat_iterations < fwork->fermat_max_iterations)
		return state_fermat;	// always do this state if not done
	else if (fwork->rho_bases < fwork->rho_max_bases)
		return state_rho;		// always do this state if not done
	else if (fwork->pp1_lvl1_curves < fwork->pp1_max_lvl1_curves)
		next_state = state_pp1_lvl1;
	else if (fwork->pm1_lvl1_curves < fwork->pm1_max_lvl1_curves)
		next_state = state_pm1_lvl1;
	else if (fwork->ecm_15digit_curves < fwork->ecm_max_15digit_curves)
		next_state = state_ecm_15digit;
	else if (fwork->ecm_20digit_curves < fwork->ecm_max_20digit_curves)
		next_state = state_ecm_20digit;
	else if (fwork->ecm_25digit_curves < fwork->ecm_max_25digit_curves)
		next_state = state_ecm_25digit;
	else if (fwork->pp1_lvl2_curves < fwork->pp1_max_lvl2_curves)
		next_state = state_pp1_lvl2;
	else if (fwork->pm1_lvl2_curves < fwork->pm1_max_lvl2_curves)
		next_state = state_pm1_lvl2;
	else if (fwork->ecm_30digit_curves < fwork->ecm_max_30digit_curves)
		next_state = state_ecm_30digit;
	else if (fwork->pp1_lvl3_curves < fwork->pp1_max_lvl3_curves)
		next_state = state_pp1_lvl3;
	else if (fwork->pm1_lvl3_curves < fwork->pm1_max_lvl3_curves)
		next_state = state_pm1_lvl3;
	else if (fwork->ecm_35digit_curves < fwork->ecm_max_35digit_curves)
		next_state = state_ecm_35digit;
	else if (fwork->ecm_40digit_curves < fwork->ecm_max_40digit_curves)
		next_state = state_ecm_40digit;
	else if (fwork->ecm_45digit_curves < fwork->ecm_max_45digit_curves)
		next_state = state_ecm_45digit;
	else if (fwork->ecm_50digit_curves < fwork->ecm_max_50digit_curves)
		next_state = state_ecm_50digit;
    else if (fwork->ecm_55digit_curves < fwork->ecm_max_55digit_curves)
        next_state = state_ecm_55digit;
    else if (fwork->ecm_60digit_curves < fwork->ecm_max_60digit_curves)
        next_state = state_ecm_60digit;
    else if (fwork->ecm_65digit_curves < fwork->ecm_max_65digit_curves)
        next_state = state_ecm_65digit;
	else
		next_state = state_nfs;

	// modify according to user preferences if necessary
	switch (next_state)
	{
		case state_ecm_15digit:
		case state_ecm_20digit:
		case state_ecm_25digit:
		case state_ecm_30digit:
		case state_ecm_35digit:
		case state_ecm_40digit:
		case state_ecm_45digit:
		case state_ecm_50digit:
		case state_ecm_55digit:
		case state_ecm_60digit:
		case state_ecm_65digit:
			if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_NOECM)
				next_state = state_nfs;
			break;
			
		default:

			break;
	}

	if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_NONE)
		next_state = state_nfs;

	return next_state;
}

int check_tune_params(fact_obj_t *fobj)
{
	if (fobj->qs_obj.qs_multiplier == 0 || 
		fobj->qs_obj.qs_exponent == 0 || 
		fobj->qs_obj.qs_tune_freq == 0 ||
		fobj->nfs_obj.gnfs_multiplier == 0 || 
		fobj->nfs_obj.gnfs_exponent == 0 || 
		fobj->nfs_obj.gnfs_tune_freq == 0)
	{
        if (fobj->VFLAG > 0)
        {
            printf("fac: check tune params contained invalid parameter(s), ignoring tune info.\n");
        }

        if (fobj->VFLAG > 2)
        {
            printf("\tqs_mult = %e\n", fobj->qs_obj.qs_multiplier);
            printf("\tqs_exp = %e\n", fobj->qs_obj.qs_exponent);
            printf("\tqs_freq = %e\n", fobj->qs_obj.qs_tune_freq);
            printf("\tnfs_mult = %e\n", fobj->nfs_obj.gnfs_multiplier);
            printf("\tnfs_exp = %e\n", fobj->nfs_obj.gnfs_exponent);
            printf("\tnfs_freq = %e\n", fobj->nfs_obj.gnfs_tune_freq);
        }
		return 0;
	}

	return 1;
}

void set_work_params(factor_work_t *fwork, enum factorization_state state)
{
    switch (state)
	{
	case state_pp1_lvl1:
        fwork->B1 = 25000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 0;
		break;

	case state_pp1_lvl2:
        fwork->B1 = 750000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 0;
		break;

	case state_pp1_lvl3:
        fwork->B1 = 2500000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 0;
		break;

	case state_pm1_lvl1:
		fwork->B1 = 150000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 1;
		break;

	case state_pm1_lvl2:
		fwork->B1 = 3750000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 1;
		break;

	case state_pm1_lvl3:
		fwork->B1 = 15000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = 1;
		break;

	case state_ecm_15digit:
		fwork->B1 = 2000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_15digit_curves;
		break;

	case state_ecm_20digit:
		fwork->B1 = 11000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_20digit_curves;
		break;

	case state_ecm_25digit:
		fwork->B1 = 50000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_25digit_curves;
		break;

	case state_ecm_30digit:
		fwork->B1 = 250000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_30digit_curves;
		break;

	case state_ecm_35digit:
		fwork->B1 = 1000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_35digit_curves;
		break;

	case state_ecm_40digit:
		fwork->B1 = 3000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_40digit_curves;
		break;

	case state_ecm_45digit:
		fwork->B1 = 11000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_45digit_curves;
		break;

	case state_ecm_50digit:
		fwork->B1 = 43000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_50digit_curves;
		break;

	case state_ecm_55digit:
		fwork->B1 = 110000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_55digit_curves;
		break;

	case state_ecm_60digit:
		fwork->B1 = 260000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_60digit_curves;
		break;

	case state_ecm_65digit:
		fwork->B1 = 850000000;
		fwork->B2 = 0;	//gmp-ecm default
		fwork->curves = fwork->ecm_max_65digit_curves;
		break;

	default:
		fwork->B1 = 0;	//error condition
		fwork->B2 = 0;
		fwork->curves = 0;
		break;

	}

	return;
}

uint32_t  get_ecm_curves_done(factor_work_t *fwork, enum factorization_state state)
{
	uint32_t  curves_done;

	switch (state)
	{
	case state_ecm_15digit:
		curves_done = fwork->ecm_15digit_curves;
		break;
	case state_ecm_20digit:
		curves_done = fwork->ecm_20digit_curves;
		break;
	case state_ecm_25digit:
		curves_done = fwork->ecm_25digit_curves;
		break;
	case state_ecm_30digit:
		curves_done = fwork->ecm_30digit_curves;
		break;
	case state_ecm_35digit:
		curves_done = fwork->ecm_35digit_curves;
		break;
	case state_ecm_40digit:
		curves_done = fwork->ecm_40digit_curves;
		break;
	case state_ecm_45digit:
		curves_done = fwork->ecm_45digit_curves;
		break;
	case state_ecm_50digit:
		curves_done = fwork->ecm_50digit_curves;
		break;
	case state_ecm_55digit:
		curves_done = fwork->ecm_55digit_curves;
		break;
	case state_ecm_60digit:
		curves_done = fwork->ecm_60digit_curves;
		break;
	case state_ecm_65digit:
		curves_done = fwork->ecm_65digit_curves;
		break;
	default:
		curves_done = 0;
		break;
	}

	return curves_done;
}

uint32_t  set_ecm_curves_done(factor_work_t *fwork, enum factorization_state state, uint32_t  curves_done)
{
	switch (state)
	{
	case state_ecm_15digit:
		fwork->ecm_15digit_curves = curves_done;
		break;
	case state_ecm_20digit:
		fwork->ecm_20digit_curves = curves_done;
		break;
	case state_ecm_25digit:
		fwork->ecm_25digit_curves = curves_done;
		break;
	case state_ecm_30digit:
		fwork->ecm_30digit_curves = curves_done;
		break;
	case state_ecm_35digit:
		fwork->ecm_35digit_curves = curves_done;
		break;
	case state_ecm_40digit:
		fwork->ecm_40digit_curves = curves_done;
		break;
	case state_ecm_45digit:
		fwork->ecm_45digit_curves = curves_done;
		break;
	case state_ecm_50digit:
		fwork->ecm_50digit_curves = curves_done;
		break;
	case state_ecm_55digit:
		fwork->ecm_55digit_curves = curves_done;
		break;
	case state_ecm_60digit:
		fwork->ecm_60digit_curves = curves_done;
		break;
	case state_ecm_65digit:
		fwork->ecm_65digit_curves = curves_done;
		break;
	default:
		printf("don't know how to set curves for state %d\n", state);
		exit(1);
		break;
	}

	return curves_done;
}

uint32_t  get_max_ecm_curves(factor_work_t *fwork, enum factorization_state state)
{
	uint32_t  max_curves;

	switch (state)
	{
	case state_ecm_15digit:
		max_curves = fwork->ecm_max_15digit_curves;
		break;
	case state_ecm_20digit:
		max_curves = fwork->ecm_max_20digit_curves;
		break;
	case state_ecm_25digit:
		max_curves = fwork->ecm_max_25digit_curves;
		break;
	case state_ecm_30digit:
		max_curves = fwork->ecm_max_30digit_curves;
		break;
	case state_ecm_35digit:
		max_curves = fwork->ecm_max_35digit_curves;
		break;
	case state_ecm_40digit:
		max_curves = fwork->ecm_max_40digit_curves;
		break;
	case state_ecm_45digit:
		max_curves = fwork->ecm_max_45digit_curves;
		break;
	case state_ecm_50digit:
		max_curves = fwork->ecm_max_50digit_curves;
		break;
	case state_ecm_55digit:
		max_curves = fwork->ecm_max_55digit_curves;
		break;
	case state_ecm_60digit:
		max_curves = fwork->ecm_max_60digit_curves;
		break;
	case state_ecm_65digit:
		max_curves = fwork->ecm_max_65digit_curves;
		break;
	default:
		max_curves = 0;
		break;
	}

	return max_curves;
}

double compute_ecm_work_done(factor_work_t *fwork, int disp_levels, FILE *log, 
    int VFLAG, int LOGFLAG)
{
	// there is probably a more elegant way to do this involving dickman's function
	// or something, but we can get a reasonable estimate using empirical data
	// for our fixed set of B1/B2 values.
	double *tlevels = fwork->tlevels;
	uint32_t  curves_done;
	int i, j;

    if (LOGFLAG && (log != NULL))
    {
        logprint(log, "ecm work completed:\n");
    }

	// compute the %done of each tlevel
	for (i=0; i < NUM_ECM_LEVELS; i++)
	{
		enum factorization_state k;
		tlevels[i] = 0;
		
		for (k=state_ecm_15digit, j=0; k <= state_ecm_65digit; k++, j++)
		{            
			curves_done = get_ecm_curves_done(fwork, k);			
			tlevels[i] += (double)curves_done / ecm_data[i][j];
		}

        if ((VFLAG >= 1) && disp_levels && (tlevels[i] > 0.01))
        {
            printf("\tt%d: %1.2f\n", ecm_levels[i], tlevels[i]);
        }

        if (LOGFLAG && (log != NULL) && (tlevels[i] > 0.01))
        {
            logprint(log, "\tt%d: %1.2f\n", ecm_levels[i], tlevels[i]);
        }
	}

	// find the first one less than 1
    for (i = 0; i < NUM_ECM_LEVELS; i++)
    {
        if (tlevels[i] < 1)
        {
            break;
        }
    }

	// estimate the t level done by extrapolating between this and the previous one
	// assuming they are all spaced 5 digits apart.
    if (i == 0)
    {
        return 0;
    }
    else
    {
        return ecm_levels[i - 1] + 5 * tlevels[i];
    }

}

enum factorization_state schedule_work(factor_work_t *fwork, mpz_t b, fact_obj_t *fobj)
{
	int have_tune;
    int i;
	enum factorization_state next_state;
	int numdigits = gmp_base10(b);
	double target_digits;
	double work_done;
	FILE *flog;
	snfs_t *poly;

	// get the next factorization state that hasn't been completed
	next_state = get_next_state(fwork, fobj);

	// make sure a minimum amount of work is done
	if (next_state == state_trialdiv ||
		next_state == state_fermat ||
		next_state == state_rho)
	{
		return next_state;
	}

	// check to see if 'tune' has been run or not
	have_tune = check_tune_params(fobj);		

	// set target pretesting depth, depending on user selection and whether or not
	// the input is both big enough and snfsable.  First we need to figure out
	// if the input is SNFSable, if we haven't already (has_snfs_form < 0), and
	// if large enough.
    if ((numdigits >= fobj->autofact_obj.qs_snfs_xover) && (fobj->autofact_obj.has_snfs_form < 0))
    {
        mpz_set(fobj->nfs_obj.gmp_n, b);
#ifdef USE_NFS

		if (fobj->nfs_obj.skip_snfs_check)
		{
			poly = NULL;
		}
		else
		{
			poly = snfs_find_form(fobj);
		}

        if (poly != NULL)
        {
			printf("fac: found form %d\n", (int)poly->form_type);
            fobj->autofact_obj.has_snfs_form = (int)poly->form_type;

			if (poly->form_type == SNFS_BRENT)
			{
				// are there algebraic factors we can extract?
				find_primitive_factor(fobj, poly, fobj->primes, fobj->num_p, fobj->VFLAG);

				// if we found any, set the autofactor target to whatever is left over.
				mpz_set(b, fobj->nfs_obj.gmp_n);
				gmp_printf("fac: continuing autofactor on residue %Zd\n", b);
			}
			else if (poly->form_type == SNFS_LUCAS)
			{
				// are there algebraic factors we can extract?
				//find_primitive_factor_lucas(fobj, poly, fobj->primes, fobj->num_p, fobj->VFLAG);

				// if we found any, set the autofactor target to whatever is left over.
				mpz_set(b, fobj->nfs_obj.gmp_n);
				gmp_printf("fac: continuing autofactor on residue %Zd\n", b);
			}

            // The actual poly is not needed now, so just get rid of it.
            snfs_clear(poly);
            free(poly);
        }
        else
        {
            fobj->autofact_obj.has_snfs_form = 0;
        }
#else
		fobj->autofact_obj.has_snfs_form = 0;
#endif
	}
	
	// set an ECM target depth based on any user-supplied plan
    if (fobj->autofact_obj.only_pretest > 1)
    {
        target_digits = fobj->autofact_obj.only_pretest;
    }
    else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_DEEP)
    {
        target_digits = 1. * (double)numdigits / 3.;
    }
    else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_LIGHT)
    {
        target_digits = 2. * (double)numdigits / 9.;
    }
    else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_CUSTOM)
    {
        target_digits = (double)numdigits * fobj->autofact_obj.target_pretest_ratio;
    }
    else
    {
        target_digits = 4. * (double)numdigits / 13.;
    }

#ifdef USE_NFS
	if ((fobj->autofact_obj.has_snfs_form > 0) && (fobj->nfs_obj.gnfs == 0) && 
        (strcmp(fobj->autofact_obj.plan_str,"custom") != 0)) // &&
        //(fobj->autofact_obj.only_pretest <= 1))
	{
		// 1) we found a snfs polynomial for the input.
		// 2) user has not specifically chosen gnfs.
		// 3) user has not already adjusted the ECM plan.
        // 4) user has not specified an only-pretest bound.
		// So it's looking like this will be an SNFS job and we should scale back the ECM effort .
		// however we also need to check if easier by gnfs... this involves a little more work.
		// This work will be duplicated if we actually get to SNFS (i.e., we don't find an
		// ECM factor), but this goes fast.
		// temporarily set verbosity silent so we don't spam the screen.  
		int tmpV = fobj->VFLAG;
		snfs_t *polys = NULL;
		int npoly;
		int gnfs_size;

		if (fobj->VFLAG > 0)
		{
			printf("fac: generating an SNFS polynomial to assess ECM effort\n");
		}

        //fobj->VFLAG = -1;

		mpz_set(fobj->nfs_obj.gmp_n, b);
		poly = snfs_find_form(fobj);

		// with the form detected, create a good polynomial
        if (poly->form_type == SNFS_XYYXF)
        {
            polys = gen_xyyxf_poly(fobj, poly, &npoly);
        }
		else if (poly->form_type == SNFS_LUCAS)
		{
			polys = gen_lucas_poly(fobj, poly, &npoly);
		}
		else
        {
            polys = gen_brent_poly(fobj, poly, &npoly);
        }

		// it's possible that poly generation found a primitive factor
		// so we need to update our input.  If not this should
		// have no effect.
		mpz_set(b, fobj->nfs_obj.gmp_n);

		// it's also possible that a detected primitive factor is prime,
		// and completes the factorization.  If so we are done.
		if (mpz_cmp_ui(b, 1) == 0)
		{
			snfs_clear(poly);
			free(poly);
			for (i = 0; i < npoly; i++)
			{
				snfs_clear(&polys[i]);
			}
			free(polys);
			return next_state;
		}

		//gmp_printf("after polygen, b = %Zd, nfs->gmp_n = %Zd, primitive = %Zd\n",
		//	b, fobj->nfs_obj.gmp_n, poly->primitive);

		// then scale and rank them
		snfs_scale_difficulty(polys, npoly, fobj->VFLAG);
		npoly = snfs_rank_polys(fobj, polys, npoly);

        fobj->VFLAG = tmpV;

		// and test the best one, compared to gnfs or qs, depending on 
        // which one will run.  
		if (mpz_cmp_ui(poly->primitive, 0) > 0)
		{
			// if gen_brent_poly found a primitive factor, it is put into
			// the factor list.  if what's left over is below the QS threshold
			// then we can stop with the SNFS detection
			gnfs_size = mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10);
			if (gnfs_size < fobj->autofact_obj.qs_snfs_xover)
			{
				// if the primitive factor size is small enough
				// we don't need to bother with NFS at all (S or G).
				// don't consider the qs/snfs cutoff any more
				fobj->autofact_obj.has_snfs_form = 0;
			}
		}
		else
		{
			gnfs_size = 999999;
		}

		if (npoly > 0)
		{
			// pick the smaller of the gnfs-equivalent SNFS size of the best
			// polynomial or the size of the primitive factor, if one was found.
			gnfs_size = MIN(gnfs_size, est_gnfs_size_via_poly(&polys[0]));
		}
		else if (fobj->VFLAG >= 0)
		{
			printf("fac: found no SNFS polynomials\n");
		}
		
		if (gnfs_size < fobj->autofact_obj.qs_snfs_xover)
		{
			// this is the case where we're keeping a small primitive 
			// factor to continue working on.  reassess target_digits for ECM.
			numdigits = gnfs_size;
			
			if (fobj->autofact_obj.only_pretest > 1)
			{

			}
			else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_DEEP)
			{
				target_digits = 1. * (double)numdigits / 3.;
			}
			else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_LIGHT)
			{
				target_digits = 2. * (double)numdigits / 9.;
			}
			else if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_CUSTOM)
			{
				target_digits = (double)numdigits * fobj->autofact_obj.target_pretest_ratio;
			}
			else
			{
				target_digits = 4. * (double)numdigits / 13.;
			}

			if (fobj->VFLAG >= 0)
			{
				if (fobj->autofact_obj.only_pretest > 1)
				{
					printf("fac: ecm effort is %1.2f for primitive factor of size %d\n",
						target_digits, numdigits);
					printf("fac: ecm effort maintained at %1.2f due to pretest condition\n",
						target_digits);
				}
				else
				{
					printf("fac: ecm effort reset to %1.2f for primitive factor of size %d\n",
						target_digits, numdigits);
				}
			}

		}
		else if ((gnfs_size <= (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) + 3)) &&
			(npoly > 0))
		{
			// Finally - to the best of our knowledge this will be a SNFS job.
			// Since we are in factor(), we'll proceed with any ecm required, but adjust 
			// the plan ratio in accord with the snfs job.
			if (fobj->VFLAG >= 0)
			{
				if (fobj->autofact_obj.only_pretest > 1)
				{
					printf("fac: ecm effort is %1.2f: input has snfs form\n",
						target_digits / 1.2857);
					printf("fac: ecm effort maintained at %1.2f due to pretest condition\n",
						target_digits);
				}
				else
				{
					printf("fac: ecm effort reduced from %1.2f to %1.2f: input has snfs form\n",
						target_digits, target_digits / 1.2857);
				}
			}

			if (fobj->autofact_obj.only_pretest > 1)
			{

			}
			else
			{
				target_digits /= 1.2857;
			}
		}
        else
        {
            if (mpz_sizeinbase(fobj->nfs_obj.gmp_n, 10) < fobj->autofact_obj.qs_gnfs_xover)
            {
                // don't consider the qs/snfs cutoff any more
                fobj->autofact_obj.has_snfs_form = 0;

                if (fobj->VFLAG >= 0)
                {
                    printf("fac: ecm effort maintained at %1.2f: input better by qs\n",
                        target_digits);
                }
            }
            else
            {
                if (fobj->VFLAG >= 0)
                {
                    printf("fac: ecm effort maintained at %1.2f: input better by gnfs\n",
                        target_digits);
                }
            }
        }

        // don't need the poly anymore
        snfs_clear(poly);
        free(poly);

        // or the list of best polys
        for (i = 0; i<npoly; i++)
        {
            snfs_clear(&polys[i]);
        }
        free(polys);
	}
	else if (fobj->nfs_obj.gnfs == 1) 
	{
		if (fobj->VFLAG > 0) printf("fac: skipping SNFS polygen for ECM effort detection - GNFS specified\n");
	}
	else if (strcmp(fobj->autofact_obj.plan_str, "custom") == 0)
	{
		if (fobj->VFLAG > 0) printf("fac: skipping SNFS polygen for ECM effort detection - custom pretest plan specified\n");
	}
	//else if (fobj->autofact_obj.only_pretest > 1)
	//{
	//	if (fobj->VFLAG > 0) printf("fac: skipping SNFS polygen for ECM effort detection - pretest bound specified\n");
	//}

#endif

	// get the current amount of work done - only print status prior to 
	// ecm steps
	switch (next_state)
	{
		case state_ecm_15digit:
		case state_ecm_20digit:
		case state_ecm_25digit:
		case state_ecm_30digit:
		case state_ecm_35digit:
		case state_ecm_40digit:
		case state_ecm_45digit:
		case state_ecm_50digit:
		case state_ecm_55digit:
		case state_ecm_60digit:
		case state_ecm_65digit:
			if (fobj->VFLAG >= 1)
				printf("fac: setting target pretesting digits to %1.6f\n", target_digits);
			
			work_done = compute_ecm_work_done(fwork, 1, NULL, fobj->VFLAG, fobj->LOGFLAG);
			fobj->autofact_obj.ecm_total_work_performed = work_done;
			
			if (fobj->VFLAG >= 1)
				printf("fac: estimated sum of completed work is t%1.6f\n", work_done);

			break;

		default:
			work_done = compute_ecm_work_done(fwork, 0, NULL, fobj->VFLAG, fobj->LOGFLAG);
			fobj->autofact_obj.ecm_total_work_performed = work_done;
			break;
	}

	// if there is a trivial amount of ecm to do, skip directly to a sieve method
	int is_trivial = 0;
	if ((target_digits < 15) && (numdigits <= 45))
	{
		if (fobj->VFLAG > 0)
			printf("fac: trivial ECM work to do... skipping to sieve method\n");
		next_state = state_nfs;
		is_trivial = 1;
	}

	// handle the case where the next state is a sieve method
	if ((next_state == state_nfs) || (work_done > target_digits) ||
		((work_done > fobj->autofact_obj.only_pretest) && 
		(fobj->autofact_obj.only_pretest > 1)) ||
		is_trivial)
	{
		logprint_oc(fobj->flogname, "a", "final ECM pretested depth: %1.6f\n", work_done);

		// if the user specified -pretest, with or without arguments,
		// we should stop factoring now that ecm is done.  this covers the
		// case where the user specified a pretest work amount that was
		// too large as determined by factor
		if ((fobj->autofact_obj.only_pretest) && !is_trivial)
		{
			logprint_oc(fobj->flogname, "a", "scheduler: pretesting active, now finishing\n");
			return state_done;
		}

		logprint_oc(fobj->flogname, "a", "scheduler: switching to sieve method\n");

		if (!have_tune || fobj->autofact_obj.prefer_xover)
		{
			// use a hard cutoff - within reason
            if ((((numdigits > fobj->autofact_obj.qs_snfs_xover) && 
                (fobj->autofact_obj.has_snfs_form)) ||
                (numdigits > fobj->autofact_obj.qs_gnfs_xover)) &&
				(numdigits >= 75))
				next_state = state_nfs;
			else
				next_state = state_qs;
		}
		else
		{
			double qs_time_est, gnfs_time_est;
		
			// compute the time to factor using estimates derived during 'tune'.
			qs_time_est = get_qs_time_estimate(fobj, b);
			gnfs_time_est = get_gnfs_time_estimate(fobj, b);

            if (fobj->VFLAG > 0)
            {
                printf("fac: tune params predict %1.2f sec for SIQS and %1.2f sec for NFS\n",
                    qs_time_est, gnfs_time_est);
            }

            if (qs_time_est < gnfs_time_est)
            {
                if (fobj->VFLAG > 0)
                {
                    printf("fac: tune params scheduling SIQS work\n");
                }
				next_state = state_qs;
            }
            else
            {
                if (fobj->VFLAG > 0)
                {
                    printf("fac: tune params scheduling NFS work\n");
                }
				next_state = state_nfs;
            }
		}
	}

	if (next_state == state_nfs) 
	{
		if ((fobj->autofact_obj.max_nfs > 0) && (gmp_base10(b) > fobj->autofact_obj.max_nfs))
		{
			if (fobj->VFLAG > 0)
			{
				printf("fac: NFS job size larger than specified maximum, finishing\n");
			}
			return state_done;
		}
		else
		{
			return next_state;
		}
	}

	if (next_state == state_qs)
	{
		if ((fobj->autofact_obj.max_siqs > 0) && (gmp_base10(b) > fobj->autofact_obj.max_siqs))
		{
			if (fobj->VFLAG > 0)
			{
				printf("fac: SIQS job size larger than specified maximum, finishing\n");
			}
			return state_done;
		}
		else
		{
			return next_state;
		}
	}

	// set the work parameters for the current state
	set_work_params(fwork, next_state);

	switch (next_state)
	{
		case state_ecm_15digit:
		case state_ecm_20digit:
		case state_ecm_25digit:
		case state_ecm_30digit:
		case state_ecm_35digit:
		case state_ecm_40digit:
		case state_ecm_45digit:
		case state_ecm_50digit:
		case state_ecm_55digit:
		case state_ecm_60digit:
		case state_ecm_65digit:
			// figure out how many curves at this level need to be done 
			// to get to the target level
			interp_and_set_curves(fwork, fobj, next_state, work_done,
				target_digits, fobj->LOGFLAG);

			break;

		default:
			// non-ecm curves are set with set_work_params above
			break;
	}

	return next_state;
}

void interp_and_set_curves(factor_work_t *fwork, fact_obj_t *fobj, 
	enum factorization_state state, double work_done,
	double target_digits, int log_results)
{
	// do a binary search on the target state's amount of work.
	// probably there is a more elegant way to compute this, but this seems
	// to work.
	double work_low, work_high, work;
	uint32_t  tmp_curves;

	// if there is a user specified pretest value, use it, regardless if it
    // means over or under ecm'ing something.
    if (fobj->autofact_obj.only_pretest > 1)
    {
        //target_digits = MIN(target_digits, fobj->autofact_obj.only_pretest);
        target_digits = fobj->autofact_obj.only_pretest;
    }


	work_low = get_ecm_curves_done(fwork, state);
	work_high = get_max_ecm_curves(fwork, state);		
	work = (work_low + work_high) / 2;

    if (fobj->VFLAG >= 1)
    {
        printf("fac: work done at B1=%u: %1.0f curves, max work = %1.0f curves\n",
            fwork->B1, work_low, work_high);
    }

	tmp_curves = work_low;		
	while ((work_high - work_low) > 1)
	{
        double compute;

		set_ecm_curves_done(fwork, state, (uint32_t )work);       
        compute = compute_ecm_work_done(fwork, 0, NULL, fobj->VFLAG, fobj->LOGFLAG);

		if (compute > target_digits)
		{
			work_high = work;
			work = (work_high + work_low) / 2;							
		}
		else					
		{
			work_low = work;
			work = (work_high + work_low) / 2;							
		}
	}

	
    set_ecm_curves_done(fwork, state, tmp_curves);    
	fwork->curves = (uint32_t )ceil(work);

    if ((tmp_curves + fwork->curves) > get_max_ecm_curves(fwork, state))
    {
        fwork->curves = get_max_ecm_curves(fwork, state) - tmp_curves;
    }

    if ((fobj->VFLAG >= 1) && fobj->LOGFLAG)
    {
        printf("fac: %u more curves at B1=%u needed to get to t%1.6f\n",
            fwork->curves, fwork->B1, target_digits);
    }

    if (fobj->LOGFLAG)
	{
		FILE *flog;
		flog = fopen(fobj->flogname,"a");
		if (flog != NULL)
		{
			logprint(flog, "current ECM pretesting depth: %1.6f\n", work_done);
			logprint(flog, "scheduled %u curves at B1=%u toward target "
				"pretesting depth of %1.6f\n", fwork->curves, fwork->B1, target_digits);
			fclose(flog);
		}
	}

	return;
}

void init_factor_work(factor_work_t *fwork, fact_obj_t *fobj)
{
	enum factorization_state interp_state = state_idle;
	// initialize max allowed work fields (note: maybe this structure should
	// be visible to the top level driver so that the user can edit values in it).
	// default values taken from gmp-ecm README, version 6.3
	fwork->ecm_max_15digit_curves = 30;		//2k
	fwork->ecm_max_20digit_curves = 74;		//11k
	fwork->ecm_max_25digit_curves = 214;	//50k
	fwork->ecm_max_30digit_curves = 430;	//250k
	fwork->ecm_max_35digit_curves = 904;	//1M
	fwork->ecm_max_40digit_curves = 2350;	//3M
	fwork->ecm_max_45digit_curves = 4480;	//11M
	fwork->ecm_max_50digit_curves = 7553;	//43M
	fwork->ecm_max_55digit_curves = 17769;	//110M
	fwork->ecm_max_60digit_curves = 42017;	//260M
	fwork->ecm_max_65digit_curves = 69408;	//850M
	fwork->tdiv_max_limit = fobj->div_obj.limit;
    fwork->fermat_iterations = 0;
	fwork->fermat_max_iterations = fobj->div_obj.fmtlimit;
	fwork->rho_max_bases = 3;
	fwork->rho_max_iterations = fobj->rho_obj.iterations;
	fwork->pp1_max_lvl1_curves = 0;
	fwork->pp1_max_lvl2_curves = 0;
	fwork->pp1_max_lvl3_curves = 0;
	fwork->pm1_max_lvl1_curves = 1;
	fwork->pm1_max_lvl2_curves = 1;
	fwork->pm1_max_lvl3_curves = 1;	
	fwork->total_time = 0;	
	fwork->trialdiv_time = 0;
	fwork->rho_time = 0;
	fwork->pp1_time = 0;
	fwork->pm1_time = 0;
	fwork->ecm_time = 0;
	fwork->qs_time = 0;
	fwork->nfs_time = 0;

	fwork->rho_bases = 0;

	fwork->pp1_lvl1_curves = 0;
	fwork->pp1_lvl2_curves = 0;
	fwork->pp1_lvl3_curves = 0;

	fwork->pm1_lvl1_curves = 0;
	fwork->pm1_lvl2_curves = 0;
	fwork->pm1_lvl3_curves = 0;

	fwork->ecm_15digit_curves = 0;
	fwork->ecm_20digit_curves = 0;
	fwork->ecm_25digit_curves = 0;
	fwork->ecm_30digit_curves = 0;
	fwork->ecm_35digit_curves = 0;
	fwork->ecm_40digit_curves = 0;
	fwork->ecm_45digit_curves = 0;
	fwork->ecm_50digit_curves = 0;
	fwork->ecm_55digit_curves = 0;
	fwork->ecm_60digit_curves = 0;
	fwork->ecm_65digit_curves = 0;

	int i;
	for (i = 0; i < NUM_ECM_LEVELS; i++)
	{
		fwork->tlevels[i] = 0.0;
	}
    
	// preload work structure with curves appropriate to the amount
	// of specified initial work
	if (fwork->initial_work >= 60.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_max_45digit_curves = 0;
		fwork->ecm_max_50digit_curves = 0;
		fwork->ecm_max_55digit_curves = 0;
		fwork->ecm_60digit_curves = get_max_ecm_curves(fwork, state_ecm_60digit);				
		interp_state = state_ecm_65digit;
	}
    else if (fwork->initial_work >= 55.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_max_45digit_curves = 0;
		fwork->ecm_max_50digit_curves = 0;
		fwork->ecm_55digit_curves = get_max_ecm_curves(fwork, state_ecm_55digit);				
		interp_state = state_ecm_60digit;
	}
    else if (fwork->initial_work >= 50.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_max_45digit_curves = 0;
		fwork->ecm_50digit_curves = get_max_ecm_curves(fwork, state_ecm_50digit);			
		interp_state = state_ecm_55digit;
	}
    else if (fwork->initial_work >= 45.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_max_40digit_curves = 0;
		fwork->ecm_45digit_curves = get_max_ecm_curves(fwork, state_ecm_45digit);				
		interp_state = state_ecm_50digit;
	}
    else if (fwork->initial_work >= 40.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_max_35digit_curves = 0;
		fwork->ecm_40digit_curves = get_max_ecm_curves(fwork, state_ecm_40digit);		
		interp_state = state_ecm_45digit;
	}
    else if (fwork->initial_work >= 35.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pp1_max_lvl3_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl3_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;		
		fwork->ecm_max_30digit_curves = 0;
		fwork->ecm_35digit_curves = get_max_ecm_curves(fwork, state_ecm_35digit);
		interp_state = state_ecm_40digit;
	}
    else if (fwork->initial_work >= 30.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pp1_max_lvl2_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl2_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_max_25digit_curves = 0;
		fwork->ecm_30digit_curves = get_max_ecm_curves(fwork, state_ecm_30digit);
		interp_state = state_ecm_35digit;
	}
    else if (fwork->initial_work >= 25.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_max_20digit_curves = 0;
		fwork->ecm_25digit_curves = get_max_ecm_curves(fwork, state_ecm_25digit);
		interp_state = state_ecm_30digit;
	}
    else if (fwork->initial_work >= 20.0)
	{
		fwork->pp1_max_lvl1_curves = 0;
		fwork->pm1_max_lvl1_curves = 0;
		fwork->ecm_max_15digit_curves = 0;
		fwork->ecm_20digit_curves = get_max_ecm_curves(fwork, state_ecm_20digit);
		interp_state = state_ecm_25digit;
	}
    else
    {
        interp_state = state_idle;
    }

    if (interp_state != state_idle)
    {

        // initializing with an indicated amount of work.  we are using this function
        // to try to figure out how many curves at the current level this is
        // equivalent too.  But the function is normally used to compute the number
        // of curves left toward the target, which could be either a digit level or
        // a pretest value.  Conclusion: temporarily set any pretest value to zero,
        // so that this function finds the equivalent curves to the indicated work.
        double tmp_pretest = fobj->autofact_obj.only_pretest;

        fobj->autofact_obj.only_pretest = 0;
        interp_and_set_curves(fwork, fobj, interp_state, fwork->initial_work,
            fwork->initial_work, 0);
            
        // restore any pretest value
        fobj->autofact_obj.only_pretest = tmp_pretest;
            
        // then fill in the equivalent curves to the indicated work amount.
        switch (interp_state)
        {
        case state_ecm_25digit: fwork->ecm_25digit_curves = fwork->curves; break;
        case state_ecm_30digit: fwork->ecm_30digit_curves = fwork->curves; break;
        case state_ecm_35digit: fwork->ecm_35digit_curves = fwork->curves; break;
        case state_ecm_40digit: fwork->ecm_40digit_curves = fwork->curves; break;
        case state_ecm_45digit: fwork->ecm_45digit_curves = fwork->curves; break;
        case state_ecm_50digit: fwork->ecm_50digit_curves = fwork->curves; break;
        case state_ecm_55digit: fwork->ecm_55digit_curves = fwork->curves; break;
        case state_ecm_60digit: fwork->ecm_60digit_curves = fwork->curves; break;
        case state_ecm_65digit: fwork->ecm_65digit_curves = fwork->curves; break;
        }

    }
	
	return;
}

int check_for_exit_on_factor(fact_obj_t* fobj)
{
	// request: https://www.mersenneforum.org/showpost.php?p=624156&postcount=262
	// -stoplt n : Stop after finding a factor with Less than n digits
	// -stople n : Stop after finding a factor with Less than or Equal to n digits
	// -stopeq n : Stop after finding a factor with n digits
	// -stopge n : Stop after finding a factor with Greater than or Equal to n digits
	// -stopgt n : Stop after finding a factor with Greater than n digits
	// -stopbase b : Base to use for stopXY options(default 10, range: 2 <= b <= 62)
	// -stopprime  : add constraint that number is also prime
	// ie : the bases supported by "mpz_get_str"

	yfactor_list_t* factors = fobj->factors;
	int base = fobj->autofact_obj.stopbase;
	int i;

	if (fobj->autofact_obj.check_stop_conditions == 0)
		return 0;

	for (i = 0; i < factors->num_factors; i++)
	{
		int sz = mpz_sizeinbase(factors->factors[i].factor, base);

		if (fobj->autofact_obj.stopprime && (factors->factors[i].type != (PRP | PRIME)))
		{
			// We require the factor to be prime and it's not, so don't need
			// to check the other conditions.
			continue;
		}

		if ((sz == fobj->autofact_obj.stopeq) && (fobj->autofact_obj.stopeq > 0))
		{
			if (fobj->VFLAG > 0)
			{
				printf("fac: found factor == %d digits in base %d, stopping.\n", 
					fobj->autofact_obj.stopeq, base);
			}
			logprint_oc(fobj->flogname, "a", "found factor == %d digits in base %d, stopping.",
				sz, base);
			return 1;
		}
		
		if ((sz <= fobj->autofact_obj.stople) && (fobj->autofact_obj.stople > 0))
		{
			if (fobj->VFLAG > 0)
			{
				printf("fac: found factor <= %d digits in base %d, stopping.\n", 
					fobj->autofact_obj.stople, base);
			}
			logprint_oc(fobj->flogname, "a", "found factor <= %d digits in base %d, stopping.",
				sz, base);
			return 1;
		}

		if ((sz >= fobj->autofact_obj.stopge) && (fobj->autofact_obj.stopge > 0))
		{
			if (fobj->VFLAG > 0)
			{
				printf("fac: found factor >= %d digits in base %d, stopping.\n", 
					fobj->autofact_obj.stopge, base);
			}
			logprint_oc(fobj->flogname, "a", "found factor >= %d digits in base %d, stopping.",
				sz, base);
			return 1;
		}

		if ((sz < fobj->autofact_obj.stoplt) && (fobj->autofact_obj.stoplt > 0))
		{
			if (fobj->VFLAG > 0)
			{
				printf("fac: found factor < %d digits in base %d, stopping.\n", 
					fobj->autofact_obj.stoplt, base);
			}
			logprint_oc(fobj->flogname, "a", "found factor < %d digits in base %d, stopping.",
				sz, base);
			return 1;
		}

		if ((sz > fobj->autofact_obj.stopgt) && (fobj->autofact_obj.stopgt > 0))
		{
			if (fobj->VFLAG > 0)
			{
				printf("fac: found factor > %d digits in base %d, stopping.\n", 
					fobj->autofact_obj.stopgt, base);
			}
			logprint_oc(fobj->flogname, "a", "found factor > %d digits in base %d, stopping.",
				sz, base);
			return 1;
		}
	}
	
	return 0;
}

void factor(fact_obj_t *fobj)
{
	//run a varity of factoring algorithms on b.
	//return any composite number left over.
	//the factoring routines will build up a list of factors.

	mpz_t b, origN, copyN;
	enum factorization_state fact_state;
	factor_work_t fwork;
	FILE *flog;
	struct timeval start, stop;
	double t_time;
	int user_defined_ecm_b2 = fobj->ecm_obj.stg2_is_default;
	int user_defined_pp1_b2 = fobj->pp1_obj.stg2_is_default;
	int user_defined_pm1_b2 = fobj->pm1_obj.stg2_is_default;
	FILE *data;
	char tmpstr[GSTR_MAXSIZE];
	int quit_after_sieve_method = 0;
    double initial_work = fobj->autofact_obj.initial_work;

	//factor() always ignores user specified B2 values
	fobj->ecm_obj.stg2_is_default = 1;
	fobj->pp1_obj.stg2_is_default = 1;
	fobj->pm1_obj.stg2_is_default = 1;

	mpz_init(origN);
	mpz_init(copyN);
	mpz_init(b);

	mpz_set(origN, fobj->N);
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

    if (fobj->LOGFLAG)
    {
        flog = fopen(fobj->flogname, "a");
        if (flog == NULL)
        {
            printf("could not open %s to append\n", fobj->flogname);
            flog = NULL;
        }
        else
        {
            logprint(flog, "\n");
            logprint(flog, "****************************\n");
        }

		char *s;
		s = mpz_get_str(NULL, 10, b);
		// use ... when we have very big numbers?
		logprint(flog,"Starting factorization of %s\n", s);
		free(s);
	}
    else
    {
        // calls to logprint won't use this because they
        // are also protected by LOGFLAG
        flog = NULL;
    }

	logprint(flog,"using pretesting plan: %s\n",fobj->autofact_obj.plan_str);
    if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_CUSTOM)
    {
        logprint(flog, "custom pretest ratio is: %1.4f\n",
            fobj->autofact_obj.target_pretest_ratio);
    }
    if (fobj->autofact_obj.only_pretest > 1)
    {
        logprint(flog, "custom pretesting limit is: %d\n",
            fobj->autofact_obj.only_pretest);
    }

	if (check_tune_params(fobj))
	{
        if (fobj->autofact_obj.prefer_xover)
        {
            logprint(flog, "using specified qs/gnfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_gnfs_xover);
            logprint(flog, "using specified qs/snfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_snfs_xover);
        }
        else
        {
            logprint(flog, "using tune info for qs/gnfs crossover\n");
        }
	}
    else
    {
        logprint(flog, "no tune info: using qs/gnfs crossover of %1.0f digits\n",
            fobj->autofact_obj.qs_gnfs_xover);
        logprint(flog, "no tune info: using qs/snfs crossover of %1.0f digits\n",
            fobj->autofact_obj.qs_snfs_xover);
    }

    // if the user input a scaling factor rather than a digit level
    // then compute the effective digit number for this input.
    if (fobj->autofact_obj.initial_work < 1.0)
    {
        initial_work = fobj->autofact_obj.initial_work * mpz_sizeinbase(fobj->N, 10);
    }

    // put the initial work done into the work structure.  It's important
    // to not modify the autofact_obj.initial_work element because if it was
    // input as a scaling factor then modifying it would destroy the scaling
    // factor for other inputs on this run (potentially a batch job).
    fwork.initial_work = initial_work;

    if (initial_work > 0.0)
    {       
        logprint(flog, "input indicated to have been pretested to t%1.2f\n",
            initial_work);
    }

	logprint(flog,"****************************\n");
	if (flog != NULL) fclose(flog);

	fobj->autofact_obj.autofact_active = 1;

	if (fobj->VFLAG >= 0)
	{
		gmp_printf("fac: factoring %Zd\n",b);
		printf("fac: using pretesting plan: %s\n",fobj->autofact_obj.plan_str);
		if (fobj->autofact_obj.yafu_pretest_plan == PRETEST_CUSTOM)
			printf("fac: custom pretest ratio is: %1.4f\n",fobj->autofact_obj.target_pretest_ratio);
		if (fobj->autofact_obj.only_pretest > 1)
			printf("fac: custom pretesting limit is: %d\n",fobj->autofact_obj.only_pretest);
		if (check_tune_params(fobj))
		{
            if (fobj->autofact_obj.prefer_xover)
            {
                printf("fac: using specified qs/gnfs crossover of %1.0f digits\n",
                    fobj->autofact_obj.qs_gnfs_xover);
                printf("fac: using specified qs/snfs crossover of %1.0f digits\n",
                    fobj->autofact_obj.qs_snfs_xover);
            }
            else
            {
                printf("fac: using tune info for qs/gnfs crossover\n");
            }
		}
        else
        {
            printf("fac: no tune info: using qs/gnfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_gnfs_xover);
            printf("fac: no tune info: using qs/snfs crossover of %1.0f digits\n",
                fobj->autofact_obj.qs_snfs_xover);
        }

        if (initial_work > 0.0)
        {
            printf("fac: input indicated to have been pretested to t%1.2f\n",
                initial_work);
        }

	}	

	init_factor_work(&fwork, fobj);

	//starting point of factorization effort
	fact_state = state_trialdiv;

	//check to see if a siqs savefile exists for this input	
	data = fopen(fobj->qs_obj.siqs_savefile,"r");

	if (data != NULL)
	{	
		char *substr;
		mpz_t tmpz;
		mpz_t g;

		//read in the number from the savefile
		mpz_init(tmpz);
		mpz_init(g);

		fgets(tmpstr,1024,data);
		substr = tmpstr + 2;
		mpz_set_str(tmpz, substr, 0);	//auto detect the base

		if (resume_check_input_match(tmpz, b, g, fobj->VFLAG))
		{
			if (fobj->VFLAG > 0)
				printf("fac: found siqs savefile, resuming siqs\n");

            // if the inputs don't match exactly, resume siqs on the exact
            // number in the savefile and put the cofactor (prime or composite)
            // into the factor list.  If composite it will get refactored.
            add_to_factor_list(fobj->factors, g, fobj->VFLAG, fobj->NUM_WITNESSES);

            mpz_set(b, tmpz);

			//override default choice
			fact_state = state_qs;

			// if for some reason qs doesn't find factors (such as
			// a user specified time out), don't continue ecm-ing, etc.
			quit_after_sieve_method = 1;
		}
		mpz_clear(tmpz);
		mpz_clear(g);
		fclose(data);
	}

	//check to see if a nfs job file exists for this input	
	data = fopen(fobj->nfs_obj.job_infile,"r");

	if (data != NULL)
	{	
		char *substr;
		mpz_t tmpz;
		mpz_t g;

		//read in the number from the job file
		mpz_init(tmpz);
		mpz_init(g);

		// may not be sufficient, in extreme cases...
		fgets(tmpstr,1024,data);
		substr = tmpstr + 2;
		mpz_set_str(tmpz, substr, 0);	//auto detect the base

		if (resume_check_input_match(tmpz, b, g, fobj->VFLAG))
		{
            // check if this is a snfsable number.  If the input is
            // really small and we don't check this, the resume
            // may default back to siqs.
            if ((mpz_sizeinbase(b,10) >= fobj->autofact_obj.qs_snfs_xover) && 
                (fobj->autofact_obj.has_snfs_form < 0))
            {
                snfs_t *poly;
                mpz_set(fobj->nfs_obj.gmp_n, b);
#ifdef USE_NFS
                poly = snfs_find_form(fobj);

                if (poly != NULL)
                {
                    fobj->autofact_obj.has_snfs_form = (int)poly->form_type;
                    // The actual poly is not needed now, so just get rid of it.
                    snfs_clear(poly);
                    free(poly);

                    if (fobj->VFLAG > 0)
                        printf("fac: found nfs job file and snfs form, resuming snfs\n");
                }
                else
                {
                    fobj->autofact_obj.has_snfs_form = 0;
                    if (fobj->VFLAG > 0)
                        printf("fac: found nfs job file, resuming nfs\n");
                }
#else
                fobj->autofact_obj.has_snfs_form = 0;
#endif
            }
            else
            {
                if (fobj->VFLAG > 0)
                    printf("fac: found nfs job file, resuming nfs\n");
            }		

			// remove any common factor so the input exactly matches
			// the file
			mpz_tdiv_q(b, b, g);
			mpz_set(fobj->N, b);
			mpz_set(origN, b);
			mpz_set(copyN, b);

			//override default choice
			fact_state = state_nfs;

			// if for some reason nfs doesn't find factors (such as
			// a user specified time out or -ns, -nc, etc.), 
			// don't continue ecm-ing, etc.
			quit_after_sieve_method = 1;
		}
		mpz_clear(tmpz);
		mpz_clear(g);
		fclose(data);
	}

	// state machine to factor the number using a variety of methods
	while (fact_state != state_done)
	{	
        // do the next item of work
		do_work(fact_state, &fwork, b, fobj);

        // check if we are done:
        // * number is completely factored
        // * sieve method was performed and either finished or was interrupted.
		// * one of the exit-on-factor-found conditions is met
        if (check_if_done(fobj, &fwork, origN) || check_for_exit_on_factor(fobj) ||
            (quit_after_sieve_method &&
            ((fact_state == state_qs) ||
            (fact_state == state_nfs))) ||
            ((fact_state == state_nfs) &&
            (fobj->flags == FACTOR_INTERRUPT)))
        {
            fact_state = state_done;
        }
        else if ((fact_state >= state_ecm_15digit) && (fact_state <= state_ecm_65digit))
        {
            // if we ran ecm, check the ecm exit code and
            // handle appropriately.
            if (fobj->ecm_obj.exit_cond == ECM_EXIT_ABORT)
            {
                FILE *flog;
                double work_done;

                if (fobj->LOGFLAG)
                {
                    flog = fopen(fobj->flogname, "a");
                }
                work_done = compute_ecm_work_done(&fwork, 1, flog, fobj->VFLAG, fobj->LOGFLAG);
				fobj->autofact_obj.ecm_total_work_performed = work_done;
                if (fobj->LOGFLAG)
                {
                    logprint(flog, "\testimated sum of completed work is t%1.2f\n", work_done);
                    if (flog != NULL) fclose(flog);
                }
                fact_state = state_done;
            }
        }

        // if not done, figure out the next item of work.
        if (fact_state != state_done)
        {
            fact_state = schedule_work(&fwork, b, fobj);
        }
	}

	// optionally record output in one or more file formats
	if (fobj->factors->num_factors >= 1)
	{
		// If the only factor in our array == N, then N is prime or prp...
		if (fobj->autofact_obj.want_output_primes && 
            (mpz_cmp(fobj->factors->factors[0].factor,origN) == 0))
		{
			if ((fobj->autofact_obj.op_file = fopen(fobj->autofact_obj.op_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.op_str);
			else
			{
				if (fobj->autofact_obj.want_output_expressions)
					gmp_fprintf(fobj->autofact_obj.op_file, "%Zd\n", fobj->N);
				else
					gmp_fprintf(fobj->autofact_obj.op_file, "%Zd\n", origN);
				if (fclose(fobj->autofact_obj.op_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.op_str);
			}
		}

		// If the first factor in the array != N, then is composite and we have factors...
		if (fobj->autofact_obj.want_output_factors &&
            (mpz_cmp(fobj->factors->factors[0].factor,origN) != 0))
		{
			if ((fobj->autofact_obj.of_file = fopen(fobj->autofact_obj.of_str, "a")) == NULL)
				printf(" ***Error: unable to open %s\n", fobj->autofact_obj.of_str);
			else
			{
				int i;
				//fprintf(fobj->autofact_obj.of_file, "%s\n", z2decstr(&origN,&gstr1));
				if (fobj->autofact_obj.want_output_expressions)
					gmp_fprintf(fobj->autofact_obj.of_file, "(%Zd)", fobj->N);
				else
					gmp_fprintf(fobj->autofact_obj.of_file, "%Zd", origN);
				for (i=0; i<fobj->factors->num_factors; i++)
				{
					gmp_fprintf(fobj->autofact_obj.of_file, "/%Zd", 
                        fobj->factors->factors[i].factor);
					if (fobj->factors->factors[i].count > 1)
						fprintf(fobj->autofact_obj.of_file, "^%d", 
                            fobj->factors->factors[i].count);
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
					gmp_fprintf(fobj->autofact_obj.ou_file, "%Zd\n", fobj->N);
				else
					gmp_fprintf(fobj->autofact_obj.ou_file, "%s\n", origN);
				if (fclose(fobj->autofact_obj.ou_file) != 0)
					printf(" ***Error: problem closing file %s\n", fobj->autofact_obj.ou_str);
			}
		}
	}

		
	if (mpz_cmp_ui(b, 1) != 0)
	{
		add_to_factor_list(fobj->factors, b, fobj->VFLAG, fobj->NUM_WITNESSES);

		int fid = find_in_factor_list(fobj->factors, b);

		if (fid >= 0)
		{
			char* s = mpz_get_str(NULL, 10, b);
			char prefix[10];

			switch (fobj->factors->factors[fid].type)
			{
			case PRIME:
				sprintf(prefix, "P");
				break;
			case PRP:
				sprintf(prefix, "prp");
				break;
			case COMPOSITE:
				sprintf(prefix, "c");
				break;
			case UNKNOWN:
				sprintf(prefix, "U");
				break;
			default:
				sprintf(prefix, "U");
				break;
			}

			logprint_oc(fobj->flogname, "a", "%s%d cofactor = %s\n", prefix, gmp_base10(b), s);
			free(s);
		}
		else
		{
			printf("failed to find cofactor in factor list!\n");
		}
	}
    
	mpz_set(fobj->N, b);

	gettimeofday (&stop, NULL);
    t_time = ytools_difftime(&start, &stop);
	fobj->autofact_obj.ttime = t_time;

	if (fobj->VFLAG >= 0)
		printf("Total factoring time = %6.4f seconds\n",t_time);

	logprint_oc(fobj->flogname, "a", "Total factoring time = %6.4f seconds\n",t_time);

	fobj->autofact_obj.autofact_active=0;

	if (fobj->refactor_depth == 0)
	{
		write_factor_json(fobj, &fwork, &start, &stop);
	}

	//restore flags
	fobj->ecm_obj.stg2_is_default = user_defined_ecm_b2;
	fobj->pp1_obj.stg2_is_default = user_defined_pp1_b2;
	fobj->pm1_obj.stg2_is_default = user_defined_pm1_b2;    

	mpz_clear(origN);
	mpz_clear(copyN);
	mpz_clear(b);
	return;
}

#if 1
int factor_tiny(mpz_t in, mpz_t* out,
	uint64_t* primes, uint64_t nump, uint64_t* prng)
{
	// factor input 'in', which is assumed to be <= 128 bits in size.
	// avoid the overhead associated with the main factor() routine.
	// utilize an input list of primes for trial division.
	// also accept a 64-bit PRNG seed for LCG-RNG
	mpz_t gmpf;
	mpz_init(gmpf);

	// first a bit of trial division.
	int k = 0;
	int numout = 0;
	while ((mpz_cmp_ui(in, 1) > 0) && (primes[k] < 10000) && (k < nump))
	{
		uint64_t q = primes[k];
		uint64_t r = mpz_tdiv_ui(in, q);

		if (r != 0)
		{
			k++;
		}
		else
		{
			mpz_tdiv_q_ui(in, in, q);
			mpz_init(out[numout]);
			mpz_set_64(out[numout], q);
			numout++;
		}
	}

	// survived TD, proceed to ECM
	// this is the lasieve5 cofactorization strategy, modified
	// to handle an arbitrary number of small factors.
	while (mpz_cmp_ui(in, 1) > 0)
	{
		if (mpz_probab_prime_p(in, 1) > 0)
		{
			// prime residue
			mpz_init(out[numout]);
			mpz_set(out[numout], in);
			numout++;
			break;
		}

		if (mpz_sizeinbase(in, 2) <= 64) {
			uint64_t n64 = mpz_get_ui(in);
			uint64_t f = getfactor_uecm(n64, 1, prng);
			if (f > 1)
			{
				if (prp_uecm(f) == 0)
				{
					// found a composite factor.  try P-1 and rho on the factor.
					uint64_t f1 = getfactor_upm1(f, 33);
					if (f1 > 1) {
						if (prp_uecm(f1) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f1);
							numout++;
							mpz_tdiv_q_ui(in, in, f1);
							if (prp_uecm(f / f1) == 1)
							{
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f / f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f / f1);
							}
							continue;
						}
					}
					f1 = getfactor_upm1(f, 100);
					if (f1 > 1) {
						if (prp_uecm(f1) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f1);
							numout++;
							mpz_tdiv_q_ui(in, in, f1);
							if (prp_uecm(f / f1) == 1)
							{
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f / f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f / f1);
							}
							continue;
						}
					}
					f1 = getfactor_upm1(f, 333);
					if (f1 > 1) {
						if (prp_uecm(f1) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f1);
							numout++;
							mpz_tdiv_q_ui(in, in, f1);
							if (prp_uecm(f / f1) == 1)
							{
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f / f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f / f1);
							}
							continue;
						}
					}
					int imax = 64;
					int found = 0;
					for (; imax < 8192; imax *= 2)
					{
						f1 = spbrent64(f, imax);
						if (f1 > 1) {
							if (prp_uecm(f1) == 1)
							{
								found = 1;
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f1);
								if (prp_uecm(f / f1) == 1)
								{
									mpz_init(out[numout]);
									mpz_set_ui(out[numout], f / f1);
									numout++;
									mpz_tdiv_q_ui(in, in, f / f1);
								}
								break;
							}
						}
					}
					if (!found)
					{
						printf("failed to split composite factor %lu of input %lu\n", f, n64);
						break;
					}
				}
				else
				{
					mpz_init(out[numout]);
					mpz_set_ui(out[numout], f);
					numout++;
					mpz_tdiv_q_ui(in, in, f);
				}
			}
			else
			{
				// uecm failed. try PM1, rho, then MPQS.
				f = getfactor_upm1(n64, 33);
				if (f > 1) {
					if (prp_uecm(f) == 1)
					{
						mpz_init(out[numout]);
						mpz_set_ui(out[numout], f);
						numout++;
						mpz_tdiv_q_ui(in, in, f);
						continue;
					}
				}
				f = getfactor_upm1(n64, 100);
				if (f > 1) {
					if (prp_uecm(f) == 1)
					{
						mpz_init(out[numout]);
						mpz_set_ui(out[numout], f);
						numout++;
						mpz_tdiv_q_ui(in, in, f);
						continue;
					}
				}
				f = getfactor_upm1(n64, 333);
				if (f > 1) {
					if (prp_uecm(f) == 1)
					{
						mpz_init(out[numout]);
						mpz_set_ui(out[numout], f);
						numout++;
						mpz_tdiv_q_ui(in, in, f);
						continue;
					}
				}
				int imax = 64;
				for (; imax < 8192; imax *= 2)
				{
					f = spbrent64(n64, imax);
					if (f > 1) {
						if (prp_uecm(f) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f);
							numout++;
							mpz_tdiv_q_ui(in, in, f);
							break;
						}
					}
				}
				printf("failed to find factor of %lu\n", n64);
				break;
			}
		}
		else
		{
#if 0
			if (getfactor_tecm(in, gmpf,
				mpz_sizeinbase(in, 2) / 3 - 2, &prng) > 0)
			{
				if (mpz_sizeinbase(gmpf, 2) <= max_primebits[s])
				{
					mpz_tdiv_q(fac[1], large_factors[s], fac[0]);

					// if the remaining residue is obviously too big, we're done.
					if (mpz_sizeinbase(fac[1], 2) > ((max_primebits[s] * 2)))
					{
						nf = 0;
						goto done;
					}

					// check if the residue is prime.  could again use
					// a cheaper method.
					if (mpz_probab_prime_p(fac[1], 1) > 0)
					{
						if (mpz_sizeinbase(fac[1], 2) <= max_primebits[s])
						{
							// we just completed a DLP factorization involving
							// 2 primes whos product was > 64 bits.
							nf = 2;
							goto done;
						}
						nf = 0;
						goto done;
					}

					// ok, so we have extracted one suitable factor, and the 
					// cofactor is not prime and a suitable size.  Do more work to 
					// split the cofactor.
					// todo: target this better based on expected factor size.
					uint64_t q64;
					uint64_t f64;
					if (mpz_sizeinbase(fac[1], 2) <= 64)
					{
						q64 = mpz_get_ui(fac[1]);
						f64 = getfactor_uecm(q64, 0, &pran);
						mpz_set_ui(fac[2], f64);
					}
					else
					{
						// we have a composite residue > 64 bits.  
						// use ecm first with high effort.
						getfactor_tecm(fac[1], fac[2], 32, &pran);
					}
					f64 = mpz_get_ui(fac[2]);

					if (f64 > 1)
					{
						mpz_tdiv_q_ui(fac[1], fac[1], f64);
						nf = 3;

						if (mpz_sizeinbase(fac[1], 2) > max_primebits[s]) {
							nf = 0;
						}
						if (mpz_sizeinbase(fac[2], 2) > max_primebits[s]) {
							nf = 0;
						}
						if (mpz_probab_prime_p(fac[0], 1) == 0)
						{
							nf = 0;
						}
						if (mpz_probab_prime_p(fac[1], 1) == 0)
						{
							nf = 0;
						}
						if (mpz_probab_prime_p(fac[2], 1) == 0)
						{
							nf = 0;
						}
					}
					else
					{
						// uecm/tecm failed, which does sometimes happen
						nf = mpqs_factor(fac[1], max_primebits[s], &fac);
						if (nf == 2)
						{
							// fac is now set to mpqs's statically allocated
							// set of mpz_t's.  copy in the one we found by ecm.
							nf = 3;
							mpz_set(fac[2], uecm_factors[0]);
						}
						else
						{
							nf = 0;
						}
					}
				}
				else
				{
					// check if the factor is prime.  could again use
					// a cheaper method.
					if (mpz_probab_prime_p(fac[0], 1) > 0)
					{
						// if the factor is obviously too big, give up.  This isn't a
						// failure since we haven't expended much effort yet.
						nf = 0;
					}
					else
					{
						// tecm found a composite first factor.
						// if it is obviously too big, we're done.
						if (mpz_sizeinbase(fac[0], 2) > ((max_primebits[s] * 2)))
						{
							nf = 0;
							goto done;
						}

						// isolate the 2nd smaller factor, and check its size.
						mpz_tdiv_q(fac[1], large_factors[s], fac[0]);

						if (mpz_sizeinbase(fac[1], 2) > (max_primebits[s]))
						{
							nf = 0;
							goto done;
						}

						// todo: target this better based on expected factor size.
						uint64_t q64;
						uint64_t f64;
						if (mpz_sizeinbase(fac[0], 2) <= 64)
						{
							q64 = mpz_get_ui(fac[0]);
							f64 = getfactor_uecm(q64, 0, &pran);
							mpz_set_ui(fac[2], f64);
						}
						else
						{
							// we have a composite residue > 64 bits.  
							// use ecm with high effort first
							getfactor_tecm(fac[0], fac[2], 32, &pran);
						}
						f64 = mpz_get_ui(fac[2]);

						if (f64 > 1)
						{
							mpz_tdiv_q_ui(fac[0], fac[0], f64);
							nf = 3;

							if (mpz_sizeinbase(fac[0], 2) > max_primebits[s]) {
								nf = 0;
							}
							if (mpz_sizeinbase(fac[2], 2) > max_primebits[s]) {
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[0], 1) == 0)
							{
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[1], 1) == 0)
							{
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[2], 1) == 0)
							{
								nf = 0;
							}

						}
						else
						{
							// uecm/tecm failed, which does sometimes happen
							nf = mpqs_factor(fac[0], max_primebits[s], &fac);
							if (nf == 2)
							{
								// fac is now set to mpqs's statically allocated
								// set of mpz_t's.  copy in the one we found by ecm.
								nf = 3;
								mpz_set(fac[2], uecm_factors[1]);
							}
							else
							{
								nf = 0;
							}
						}
					}
				}
			}
			else
			{
				// if ecm can't find a factor, give up.  
				// unless this is a DLP with lpbr/a > 32... i.e., if the
				// large factor size is greater than 64 bits but less than
				// lpbr/a * 2.  In that case run mpqs... or tecm with
				// greater effort.


#if 0
				if (mpz_sizeinbase(large_factors[s1], 2) <= (max_primebits[s1] * 2))
				{
					if (getfactor_tecm(large_factors[s1], factor1, 33, &pran) > 0)
					{
						if (mpz_sizeinbase(factor1, 2) <= max_primebits[s1])
						{
							mpz_tdiv_q(factor2, large_factors[s1], factor1);

							// check if the residue is prime.  could again use
							// a cheaper method.
							if (mpz_probab_prime_p(factor2, 1) > 0)
							{
								if (mpz_sizeinbase(factor2, 2) <= max_primebits[s1])
								{
									// we just completed a DLP factorization involving
									// 2 primes whos product was > 64 bits.
									mpz_set(large_primes[s1][0], factor1);
									mpz_set(large_primes[s1][1], factor2);
									nlp[s1] = 2;
								}
								else
									break;
							}
							else
								break;
						}
						else
							break;
					}
					else
						break;
				}
				else
					break;
#else

				if (mpz_sizeinbase(large_factors[s], 2) <= (max_primebits[s] * 2))
				{
					nf = mpqs_factor(large_factors[s], max_primebits[s], &fac);
				}
				else
				{
#if 0
					// try for a lucky p-1 hit on the 3LP before we go?
					// testing on an input with LPB=33 and 3LP enabled
					// saw that p-1 finds lots of factors but the residues
					// are all (99.9%) large primes.  I.e., exactly the
					// kind of inputs we want to not waste time on.
					if (getfactor_tpm1(large_factors[s], fac[0], 333))
					{
						mpz_tdiv_q(fac[1], large_factors[s], fac[0]);
						if (mpz_sizeinbase(fac[1], 2) <= max_primebits[s])
						{
							gmp_printf("P-1 Success! %Zd = %Zd * %Zd\n",
								large_factors[s], fac[0], fac[1]);
						}
						else if (mpz_probab_prime_p(fac[1], 1) == 0)
						{
							gmp_printf("Residue %Zd with %d bits is composite\n",
								fac[1], mpz_sizeinbase(fac[1], 2));
							gmp_printf("3LP = ");

							mpz_set(fac[2], fac[0]);
							nf = 1 + mpqs_factor(fac[2], max_primebits[s], &fac);

							for (i = 0; i < nf; i++)
								gmp_printf("%Zd ", fac[i]);
							printf("\n");
						}
					}
#else
					nf = 0;
#endif
				}
#endif
			}
#endif
		}
	}


	mpz_clear(gmpf);
	return numout;
}

#endif

void write_factor_json(fact_obj_t* fobj, factor_work_t *fwork,
	struct timeval *start, struct timeval *stop)
{
	FILE* fid = fopen(fobj->factor_json_name, "a");
	yfactor_list_t* flist = fobj->factors;
	int i;
	int j;
	int nump = 0;
	int numprp = 0;
	int numc = 0;
	int printedp = 0;
	int printedprp = 0;
	int printedc = 0;
	char lf = '\n';

	if (fobj->autofact_obj.json_pretty)
		lf = '\n';
	else
		lf = ' ';

	if (fid != NULL)
	{
		fprintf(fid, "{%c", lf);
		fprintf(fid, "\t\"input-expression\":\"%s\",%c", fobj->input_str, lf);
		gmp_fprintf(fid, "\t\"input-decimal\":\"%Zd\",%c", fobj->input_N, lf);

		if (fobj->argc > 1)
		{
			fprintf(fid, "\t\"input-argument-string\":\"");
			for (i = 1; i < fobj->argc; i++)
			{
				fprintf(fid, "%s ", fobj->argv[i]);
			}
			fprintf(fid, "\",%c", lf);
		}

		for (i = 0; i < flist->num_factors; i++)
		{
			if (flist->factors[i].type == PRIME)
				nump++;
		}
		for (i = 0; i < flist->num_factors; i++)
		{
			if (flist->factors[i].type == PRP)
				numprp++;
		}
		for (i = 0; i < flist->num_factors; i++)
		{
			if (flist->factors[i].type == COMPOSITE)
				numc++;
		}

		if (nump > 0)
		{
			int* printed = (int*)xmalloc(flist->num_factors * sizeof(int));
			for (i = 0; i < flist->num_factors; i++)
			{
				printed[i] = 0;
			}

			fprintf(fid, "\t\"factors-prime\":[");
			for (i = 0; i < flist->num_factors; i++)
			{
				// if already printed, move on
				if (printed[i])
					continue;

				int k;
				if (flist->factors[i].type == PRIME)
				{
					// if there is a smaller number than this one, print it first (sorts prime factors)
					k = i;
					for (j = i+1; j < flist->num_factors; j++)
					{
						if ((mpz_cmp(flist->factors[j].factor, flist->factors[i].factor) < 0) &&
							(printed[j] == 0))
						{
							k = j;		// set next index to print
							i = i - 1;	// consider this number again on the next iteration
							break;
						}
					}

					// don't redo APR-CL calculations already performed by add_to_factor_list
					for (j = 0; j < flist->factors[k].count - 1; j++)
					{
						gmp_fprintf(fid, "\"%Zd\",", flist->factors[k].factor);
					}

					if (printedp == (nump - 1))
					{
						gmp_fprintf(fid, "\"%Zd\"],%c", flist->factors[k].factor, lf);
					}
					else
					{
						gmp_fprintf(fid, "\"%Zd\",", flist->factors[k].factor);
					}
					printedp++;
					printed[k] = 1;
				}
			}
			free(printed);
		}

		if (numprp > 0)
		{
			fprintf(fid, "\t\"factors-prp\":[");
			for (i = 0; i < flist->num_factors; i++)
			{
				if (flist->factors[i].type == PRP)
				{
					for (j = 0; j < flist->factors[i].count - 1; j++)
					{
						gmp_fprintf(fid, "\"%Zd\",", flist->factors[i].factor);
					}
					if (printedprp == (numprp - 1))
					{
						gmp_fprintf(fid, "\"%Zd\"],%c", flist->factors[i].factor, lf);
					}
					else
					{
						gmp_fprintf(fid, "\"%Zd\",", flist->factors[i].factor);
					}
					printedprp++;
				}
			}
		}

		if (numc > 0)
		{
			fprintf(fid, "\t\"factors-composite\":[");
			for (i = 0; i < flist->num_factors; i++)
			{
				if (flist->factors[i].type == COMPOSITE)
				{
					for (j = 0; j < flist->factors[i].count - 1; j++)
					{
						gmp_fprintf(fid, "\"%Zd\",", flist->factors[i].factor);
					}
					if (printedc == (numc - 1))
					{
						gmp_fprintf(fid, "\"%Zd\"],%c", flist->factors[i].factor, lf);
					}
					else
					{
						gmp_fprintf(fid, "\"%Zd\",", flist->factors[i].factor);
					}
					printedc++;
				}
			}
		}

		if (fwork->pm1_lvl1_curves > 0)
		{
			fprintf(fid, "\t\"pm1-curves\" : {\"150000\":%d", fwork->pm1_lvl1_curves);
		}
		if (fwork->pm1_lvl2_curves > 0)
		{
			fprintf(fid, ",\"3750000\":%d", fwork->pm1_lvl2_curves);
		}
		if (fwork->pm1_lvl3_curves > 0)
		{
			fprintf(fid, ",\"15000000\":%d", fwork->pm1_lvl3_curves);
		}
		if (fwork->pm1_lvl1_curves > 0) fprintf(fid, "},%c", lf);

		if (fwork->pp1_lvl1_curves > 0)
		{
			fprintf(fid, "\t\"pp1-curves\" : {\"25000\":%d", fwork->pp1_lvl1_curves);
		}
		if (fwork->pp1_lvl2_curves > 0)
		{
			fprintf(fid, ",\"750000\":%d", fwork->pp1_lvl2_curves);
		}
		if (fwork->pp1_lvl3_curves > 0)
		{
			fprintf(fid, ",\"2500000\":%d", fwork->pp1_lvl3_curves);
		}
		if (fwork->pp1_lvl1_curves > 0) fprintf(fid, "},%c", lf);


		if (fwork->tlevels[0] > 0.01)
		{
			fprintf(fid, "\t\"ecm-curves\" : {");
			if (fwork->ecm_15digit_curves > 0) fprintf(fid, "\"2000\":%d", fwork->ecm_15digit_curves);
			if (fwork->ecm_20digit_curves > 0) fprintf(fid, ",\"11000\":%d", fwork->ecm_20digit_curves);
			if (fwork->ecm_25digit_curves > 0) fprintf(fid, ",\"50000\":%d", fwork->ecm_25digit_curves);
			if (fwork->ecm_30digit_curves > 0) fprintf(fid, ",\"250000\":%d", fwork->ecm_30digit_curves);
			if (fwork->ecm_35digit_curves > 0) fprintf(fid, ",\"1000000\":%d", fwork->ecm_35digit_curves);
			if (fwork->ecm_40digit_curves > 0) fprintf(fid, ",\"3000000\":%d", fwork->ecm_40digit_curves);
			if (fwork->ecm_45digit_curves > 0) fprintf(fid, ",\"11000000\":%d", fwork->ecm_45digit_curves);
			if (fwork->ecm_50digit_curves > 0) fprintf(fid, ",\"43000000\":%d", fwork->ecm_50digit_curves);
			if (fwork->ecm_55digit_curves > 0) fprintf(fid, ",\"110000000\":%d", fwork->ecm_55digit_curves);
			if (fwork->ecm_60digit_curves > 0) fprintf(fid, ",\"260000000\":%d", fwork->ecm_60digit_curves);
			if (fwork->ecm_65digit_curves > 0) fprintf(fid, ",\"850000000\":%d", fwork->ecm_65digit_curves);
			fprintf(fid, "},%c", lf);

			fprintf(fid, "\t\"ecm-levels\" : {");
			int level = 15;
			for (i = 0; i < NUM_ECM_LEVELS - 1; i++)
			{
				if ((fwork->tlevels[i] > 0.01) && (fwork->tlevels[i + 1] > 0.01))
				{
					fprintf(fid, "\"t%d\":%1.2f,", level, fwork->tlevels[i]);
				}
				else if (fwork->tlevels[i] > 0.01)
				{
					fprintf(fid, "\"t%d\":%1.2f},%c", level, fwork->tlevels[i], lf);
				}
				level += 5;
			}
			if (fwork->tlevels[i] > 0.01)
			{
				fprintf(fid, "\"t%d\":%1.2f},%c", level, fwork->tlevels[i], lf);
			}

			double work_done = compute_ecm_work_done(fwork, 0, NULL, -1, 0);

			if (work_done > 0.0)
			{
				fprintf(fid, "\"ecm-sum\":%1.2f,%c", work_done, lf);
			}
		}

		fprintf(fid, "\t\"runtime\" : {\"total\":%1.4f", fobj->autofact_obj.ttime);
		if (fobj->ecm_obj.ttime > 0.00001) fprintf(fid, ", \"ecm\":%1.4f", fobj->ecm_obj.ttime);
		if (fobj->pm1_obj.ttime > 0.00001) fprintf(fid, ", \"pm1\":%1.4f", fobj->pm1_obj.ttime);
		if (fobj->pp1_obj.ttime > 0.00001) fprintf(fid, ", \"pp1\":%1.4f", fobj->pp1_obj.ttime);
		if (fobj->qs_obj.total_time > 0.00001) fprintf(fid, ", \"siqs\":%1.4f", fobj->qs_obj.total_time);
		if (fobj->nfs_obj.ttime > 0.00001) fprintf(fid, ", \"nfs-total\":%1.4f", fobj->nfs_obj.ttime);
		if (fobj->nfs_obj.poly_time > 0.00001) fprintf(fid, ", \"nfs-poly\":%1.4f", fobj->nfs_obj.poly_time);
		if (fobj->nfs_obj.sieve_time > 0.00001) fprintf(fid, ", \"nfs-sieve\":%1.4f", fobj->nfs_obj.sieve_time);
		if (fobj->nfs_obj.filter_time > 0.00001) fprintf(fid, ", \"nfs-filter\":%1.4f", fobj->nfs_obj.filter_time);
		if (fobj->nfs_obj.la_time > 0.00001) fprintf(fid, ", \"nfs-la\":%1.4f", fobj->nfs_obj.la_time);
		if (fobj->nfs_obj.sqrt_time > 0.00001) fprintf(fid, ", \"nfs-sqrt\":%1.4f", fobj->nfs_obj.sqrt_time);
		//if (fobj->nfs_obj.poly_time > 0.0) fprintf(fid, ", \"nfs-poly-score\":%1.4f", fobj->nfs_obj.sqrt_time);
		fprintf(fid, "},%c", lf);

		char buffer[30];
		time_t curtime;

		curtime = start->tv_sec;
		strftime(buffer, 30, "%Y-%m-%d  %T", localtime(&curtime));
		fprintf(fid, "\t\"time-start\" : \"%s\",%c", buffer, lf);
		curtime = stop->tv_sec;
		strftime(buffer, 30, "%Y-%m-%d  %T", localtime(&curtime));
		fprintf(fid, "\t\"time-end\" : \"%s\",%c", buffer, lf);


		fprintf(fid, "\t\"info\":{");

#if defined(_MSC_VER) && defined(__clang_version__)
		fprintf(fid, "\"compiler\":\"MSVC %d, %s\",", _MSC_VER, __clang_version__);
#elif defined(_MSC_VER)
		fprintf(fid, "\"compiler\":\"MSVC %d\",", _MSC_VER);
#elif defined (__INTEL_COMPILER)
		fprintf(fid, "\"compiler\":\"INTEL %d\",", __INTEL_COMPILER);
#elif defined(__clang_version__)
		fprintf(fid, "\"compiler\":\"%s\",", __clang_version__);
#elif defined (__GNUC__)
		fprintf(fid, "\"compiler\":\"GNUC %d\",", __GNUC__);
#endif

#ifdef _MSC_MPIR_VERSION
#ifdef ECM_VERSION
		fprintf(fid, "\"ECM-version\":\"%s\",\"MPIR-version\":\"%s\",", ECM_VERSION,
			_MSC_MPIR_VERSION);
#elif defined(VERSION)

		fprintf(fid, "\"ECM-version\":\"%s\",\"MPIR-version\":\"%s\",", VERSION,
			_MSC_MPIR_VERSION);
#endif
#else
#ifdef ECM_VERSION
		fprintf(fid, "\"ECM-version\":\"%s\",\"GMP-version\":\"%d.%d.%d\",", ECM_VERSION,
			__GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
#endif

#endif
		fprintf(fid, "\"yafu-version\":\"%s\"}%c}\n", YAFU_VERSION_STRING, lf);

		fclose(fid);
	}
	else
	{
		printf("could not open %s to append\n", fobj->factor_json_name);
		exit(1);
	}

	return;
}



