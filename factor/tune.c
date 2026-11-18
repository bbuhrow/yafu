#include "arith.h"
#include <inttypes.h>
#include "factor.h"
#include "qs.h"
#include "gnfs-yafu.h"
#include "gmp_xface.h"
#include "yafu_ecm.h"
#include <stdint.h>
#include <math.h>
#include <time.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif

#ifdef __MINGW32__
#include <sys/time.h>
#endif

//----------------------- LOCAL DECLARATIONS ----------------------------------//
#define NUM_SIQS_PTS 10
#define NUM_GNFS_PTS 6
#define NUM_ECM_PTS 5
#define BASE_e 2.718281828459045

//----------------------- LOCAL FUNCTIONS -------------------------------------//
double best_fit(double *x, double *y, int numpts, 
	double *slope, double *intercept, int type);
void update_INI(double mult, double exponent, double mult2, 
    double exponent2, double xover,
	double b1slope, double b1intercept, double size_exp, 
	double cpu_freq, char *cpu_str);
void make_job_file(char *sname, uint32_t *startq, uint32_t *qrange, char *inputstr, int inputnum, fact_obj_t *fobj);
void check_siever(fact_obj_t* fobj, char* sname, int siever);

//----------------------- TUNE ENTRY POINT ------------------------------------//
void factor_tune(fact_obj_t *inobj)
{
	// perform trial runs on a set of inputs with known solutions using siqs and nfs.
	// compute crossover points and time estimation coefficients.
	// write all this info to the .ini file.
	// should produce roughly accurate predictions about how long qs and nfs
	// jobs will take for a given system in similar operating conditions as when this
	// routine was run (i.e., this should be run under the same nominal load from
	// other programs as you normally would expect to exist when running YAFU).
	// if any of (msieve, ggnfs, or yafu)'s routines are improved significantly, this 
	// test should be repeated (and possibly, the baseline "actuals" below should be
	// recomputed).  Likewise if the system specifications (other than cpu frequency, which
	// is assumed to scale linearly) are modified (new memory, MOBO, CPU), this test
	// should be repeated.
	char siqslist[NUM_SIQS_PTS][200];
	char nfslist[NUM_GNFS_PTS][200];
    char sievername[1024];
	int i, tmpT;
    mpz_t n;
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	double t_time;

	//uint32_t siqs_actualrels[NUM_SIQS_PTS] = {17136, 32337,
		//63709, 143984, 240663, 568071, 793434, 1205061, 1595268};
	UNUSED_VAR uint32_t siqs_actualrels[NUM_SIQS_PTS] = {17136, 32337,
		63709, 143984, 242825, 589192, 1178529, 1716902, 2457335, 3511771 }; // 2397232, 3293925
	UNUSED_VAR uint32_t siqs_actualrels_direct_avx512[NUM_SIQS_PTS] = { 27532, 44122,
		87583, 137216, 288934, 459212, 678979, 1133328, 1859293, 2675540 }; // 2397232, 3293925
	UNUSED_VAR uint32_t siqs_actualrels_direct[NUM_SIQS_PTS] = { 33129, 59729,
		72805, 130026, 211588, 281348, 651696, 1263068, 1923778, 2675540 }; // 2397232, 3293925
    uint32_t siqs_tf_small_cutoff[NUM_SIQS_PTS] = { 15, 20,
        15, 13, 16, 20, 18, 19, 20, 21 }; // 20, 20
	uint32_t siqs_testrels[NUM_SIQS_PTS] = { 10000, 10000, 10000 , 10000 , 10000,
		10000, 10000 , 20000 , 20000 , 20000 };

	double siqs_extraptime[NUM_SIQS_PTS];
	double siqs_sizes[NUM_SIQS_PTS] = {60, 65, 70, 75, 80, 85, 90, 95, 100, 105};

	uint32_t gnfs_actualrels[NUM_GNFS_PTS] = {1929527, 4074867, 4783410, 
		4969315, 5522597, 5783845};
	double gnfs_extraptime[NUM_GNFS_PTS];
	double gnfs_sizes[NUM_GNFS_PTS] = {85, 90, 95, 100, 105, 110};

	double ecm_extraptime1[3][NUM_ECM_PTS];
	double ecm_sizes[NUM_ECM_PTS];

	double a, b, a2, b2, fit, xover;
	double a3 = 0.0, b3 = 0.0;	// stay 0 when ECM tuning is skipped
	double avg_exp = 0.0;
	uint32_t count;
	char tmpbuf[GSTR_MAXSIZE];

    mpz_init(n);
	tmpT = inobj->THREADS;

	if (inobj->THREADS != 1)
		printf("Setting THREADS = 1 for tuning\n");
    inobj->THREADS = 1;

    if (inobj->VFLAG >= 0)
        printf("checking for NFS sievers... ");

    check_siever(inobj, sievername, 11);
    check_siever(inobj, sievername, 12);
    check_siever(inobj, sievername, 13);

    if (inobj->VFLAG >= 0)
        printf("done.\n");

#if defined(USE_AVX2) || defined(USE_SSE41)
    // as of version 1.35, c70 and above use DLP in SIQS for these CPUs:
    // need to modify the actual relations needed.
    siqs_actualrels[2] = 97812;
    siqs_actualrels[3] = 243854;
    siqs_actualrels[4] = 418664;
#endif

	// siqs: start with c60, increment by 5 digits, up to a c100
	// this will allow determination of NFS/QS crossover as well as provide enough
	// info to generate an equation for QS time estimation
	strcpy(siqslist[0], "349594255864176572614071853194924838158088864370890996447417");
	strcpy(siqslist[1], "34053408309992030649212497354061832056920539397279047809781589871");
	strcpy(siqslist[2], "6470287906463336878241474855987746904297564226439499503918586590778209");
	strcpy(siqslist[3], "281396163585532137380297959872159569353696836686080935550459706878100362721");
	strcpy(siqslist[4], "33727095233132290409342213138708322681737322487170896778164145844669592994743377");
	strcpy(siqslist[5], "1877138824359859508015524119652506869600959721781289179190693027302028679377371001561");
	strcpy(siqslist[6], "427351849650748515507228344120452096326780093349980867041485502247153375067354165128307841");
	strcpy(siqslist[7], "48404068520546498995797968938385187958997290617596242601254422967869040251141325866025672337021");
	strcpy(siqslist[8], "1802716097522165018257858828415111497060066282677325501816640492782221110851604465066510547671104729");
	strcpy(siqslist[9], "466734409955806375058988820327650664396976790744285564594552020197119774272189758795312820988691316775181");

	/*
	c95 dlp (default, 8t, inmem): 224/238
	c100 dlp (default, 8t, inmem): 520/543 (-siqsNB 6 502/525)
	c105 dlp (default, 8t): 2042/2090
	
	c95 tlp (tbd, inmem):
	c100 tlp: 495/509 (35.9M batched, lpb-27 mfbt-2.6 B-60k TF-95 TFsm-19 BT-250k, NB-8 Bdiv 3)
	c100 tlp: 472/487 (36.5M batched, lpb-27 mfbt-2.65 B-60k TF-95 TFsm-19 BT-250k NB-6 Bdiv 3)
	c100 tlp: 476/490 (36.5M batched, lpb-27 mfbt-2.65 B-55k TF-95 TFsm-19 BT-250k NB-6 Bdiv 3)
	c100 tlp: 487/499 (30.8M batched, lpb-26 mfbt-2.7 B-60k TF-95 TFsm-19 BT-250k NB-6 Bdiv 3)
 	
	c105 tlp: 1584/1635 (85M batched, lpb-27 mfbt-2.7 B-100k TF-100 TFsm-19 BT-500k NB-6 Bdiv 3)
	*/


	// nfs: start with c85, increment by 5 digits, up to C110
	// this will allow determination of NFS/QS crossover
	// to do NFS time estimation, probably need to go much higher - say c155.
	strcpy(nfslist[0], "1877138824359859508015524119652506869600959721781289179190693027302028679377371001561");
	strcpy(nfslist[1], "427351849650748515507228344120452096326780093349980867041485502247153375067354165128307841");
	strcpy(nfslist[2], "48404068520546498995797968938385187958997290617596242601254422967869040251141325866025672337021");
	strcpy(nfslist[3], "1802716097522165018257858828415111497060066282677325501816640492782221110851604465066510547671104729");
	strcpy(nfslist[4], "466734409955806375058988820327650664396976790744285564594552020197119774272189758795312820988691316775181");
	strcpy(nfslist[5], "48178889479314834847826896738914354061668125063983964035428538278448985505047157633738779051249185304620494013");

	// c135 from TMPQS
	// 116915434112329921568236283928181979297762987646390347857868153872054154807376462439621333455331738807075404918922573575454310187518221

	// RSA-120
	// 227010481295437363334259960947493668895875336466084780038173258247009162675779735389791151574049166747880487470296548479

	// RSA-129 (The magic words are squeamish ossifrage)
	// 114381625757888867669235779976146612010218296721242362562561842935706935245733897830597123563958705058989075147599290026879543541

	fact_obj_t *fobj = inobj;

	//init_factobj(&fobj);
	//copy_factobj(&fobj, inobj, 1);

	//goto tune_ecm;

	// for each of the siqs inputs
	for (i = 0; i < NUM_SIQS_PTS; i++)
	{
		reset_factobj(fobj);

		// measure how long it takes to gather a fixed number of relations 		
		mpz_set_str(n, siqslist[i], 10);

		fobj->qs_obj.gbl_override_rel_flag = 1;
		fobj->qs_obj.gbl_override_rel = siqs_testrels[i];

		// also set the tf_small_cutoff to its known best value so we
		// don't pollute the measurement with the optimization process.
		fobj->qs_obj.gbl_override_small_cutoff_flag = 1;
		fobj->qs_obj.gbl_override_small_cutoff = siqs_tf_small_cutoff[i];

		gettimeofday(&start, NULL);
		mpz_set(fobj->qs_obj.gmp_n, n);
		SIQS(fobj);
		gettimeofday(&stop, NULL);
		t_time = ytools_difftime(&start, &stop);

#if defined( USE_AVX512F ) //&& defined()
		siqs_extraptime[i] = t_time * siqs_actualrels_direct_avx512[i] / fobj->qs_obj.gbl_override_rel;
#else
		// the number of relations actually gathered is stored in gbl_override_rel
		siqs_extraptime[i] = t_time * siqs_actualrels_direct[i] / fobj->qs_obj.gbl_override_rel;
		//siqs_extraptime[i] = t_time * siqs_actualrels[i] / fobj->qs_obj.gbl_override_rel;
#endif

		// add a guess at the linalg + sqrt time: something reasonable seems to be about
		// 2% of the total sieve time
		siqs_extraptime[i] += 0.02 * siqs_extraptime[i];

		printf("elapsed time for ~%dk relations of c%zu = %6.4f seconds.\n",
			siqs_testrels[i] / 1000, mpz_sizeinbase(n, 10), t_time);
		printf("extrapolated time for complete factorization = %6.4f seconds\n", siqs_extraptime[i]);
	}

	fit = best_fit(siqs_sizes, siqs_extraptime, NUM_SIQS_PTS, &a, &b, 2);
	printf("best linear fit is ln(y) = %g * x + %g\nR^2 = %g\n",a,b,fit);
	printf("best exponential fit is y = %g * exp(%g * x)\n",pow(BASE_e,b),a);

	// for each of the gnfs inputs
	for (i=0; i < NUM_GNFS_PTS; i++)
	{
		char syscmd[1024];
		FILE *in;
		uint32_t startq, qrange;		
		double t_time2, d;

		// remove previous tests
#ifndef _MSC_VER
		//printf("checking if tunerels.out exists\n");
		if (access("tunerels.out", F_OK) == 0)
#endif
		{
			//printf("tunerels exists, checking for write access\n");
#ifndef _MSC_VER
			if (access("tunerels.out", W_OK) == 0)
#endif
			{
				//printf("write access confirmed, attempting to remove it\n");
				remove("tunerels.out");
			}
		}

#ifndef _MSC_VER
		//printf("checking if tune.job.afb.0 exists\n");
		if (access("tune.job.afb.0", F_OK) == 0)
#endif
		{
			//printf("tune.job.afb.0 exists, checking for write access\n");
#ifndef _MSC_VER
			if (access("tune.job.afb.0", W_OK) == 0)
#endif
			{
				//printf("write access confirmed, attempting to remove it\n");
				remove("tune.job.afb.0");
			}
		}

		MySleep(100);

		// make a job file for each input
		make_job_file(sievername, &startq, &qrange, nfslist[i], i, inobj);

		// measure how long it takes to generate the afb... for fun.		
		gettimeofday(&start, NULL);

		// create the afb - we don't want the time it takes to do this to
		// pollute the sieve timings
		sprintf(syscmd,"%s -b tune.job -k -c 0 -F", sievername);

		printf("nfs: commencing construction of afb\n");

		if (system(syscmd) != 0) { /* failure shows up when the output is checked */ }
		gettimeofday(&stop, NULL);
        t_time2 = ytools_difftime(&start, &stop);
		
		printf("afb generation took %6.4f seconds.\n",t_time2);
		remove("tune.job.afb.0");
		MySleep(100);

		// measure how long it takes to sieve a fixed range of special-q	
		gettimeofday(&start, NULL);

		// start the test
		sprintf(syscmd,"%s -f %u -c %u -o tunerels.out -a tune.job",
			sievername, startq, qrange);
		printf("nfs: commencing lattice sieving over range: %u - %u\n",
			startq, startq + qrange);
		if (system(syscmd) != 0) { /* failure shows up when the output is checked */ }
		gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);

		// count relations
		in = fopen("tunerels.out","r");
		if (in != NULL)
		{
			count = 0;
			while (fgets(tmpbuf, GSTR_MAXSIZE, in) != NULL)
				count++;
			fclose(in);
		}
		else
			count = 1;	//no divide by zero

		printf("nfs: elapsed time for %u relations of c%d = %6.4f seconds.\n",
			count,(uint32_t)gnfs_sizes[i],t_time - t_time2);

		// extrapolate
		gnfs_extraptime[i] = (double)gnfs_actualrels[i] * (t_time - t_time2) / (double)count;
		printf("nfs: extrapolated time for sieving = %6.4f seconds\n",gnfs_extraptime[i]);

		// add estimated linalg and sqrt time.  
		// reasonable estimate is 10% of the total runtime?
		d = gnfs_extraptime[i] * 0.1;
		gnfs_extraptime[i] += d;
		printf("adding %g seconds for linalg and sqrt\n",d);	

		// add avg poly time 
		// estimated at 5% of sieving time
		d = gnfs_extraptime[i] * 0.05; // gnfs_max_poly_time[i] * 3600 / 2;
		gnfs_extraptime[i] += d;
		printf("adding %g seconds for average polyselect time\n",d);

		printf("nfs: estimated time for complete factorization = %6.4f seconds\n",gnfs_extraptime[i]);
		
	}

	// set THREADS back the way it was
    inobj->THREADS = tmpT;

	fit = best_fit(gnfs_sizes, gnfs_extraptime, NUM_GNFS_PTS, &a2, &b2, 2);
	printf("best linear fit is ln(y) = %g * x + %g\nR^2 = %g\n",a2,b2,fit);
	printf("best exponential fit is y = %g * exp(%g * x)\n",pow(BASE_e,b2),a2);

	xover = (b2 - b) / (a - a2);
	printf("QS/NFS crossover occurs at %2.1f digits\n",xover);

	goto done;

tune_ecm:

	{
		// for each of the ecm inputs
		gmp_randstate_t gmp_randstate;
		gmp_randinit_default(gmp_randstate);
		mpz_t base_100b;
		int ecm_size = 200;
		int curves_run;

		mpz_init(base_100b);
		mpz_urandomb(base_100b, gmp_randstate, 100);
		mpz_nextprime(base_100b, base_100b);
		
		for (i = 0; i < NUM_ECM_PTS; i++, ecm_size += 208)
		{
			reset_factobj(fobj);

			mpz_urandomb(n, gmp_randstate, (mp_bitcnt_t)(ecm_size - 100));
			mpz_nextprime(n, n);
			mpz_mul(n, n, base_100b);
			ecm_sizes[i] = (double)mpz_sizeinbase(n, 2);

			// measure how long the total process takes for B1=10000, 100000, 1000000
			fobj->ecm_obj.B1 = 100000;
			fobj->ecm_obj.num_curves = 1;
			
			gettimeofday(&start, NULL);
			mpz_set(fobj->ecm_obj.gmp_n, n);
			mpz_set(fobj->N, n);
			curves_run = ecm_loop(fobj);
			gettimeofday(&stop, NULL);
			t_time = ytools_difftime(&start, &stop);
			
			printf("elapsed time per curve for B1=%" PRIu64 "k on c%zu = %6.4f seconds.\n",
				fobj->ecm_obj.B1 / 1000, mpz_sizeinbase(n, 10), t_time / curves_run);
			
			ecm_extraptime1[0][i] = t_time / curves_run;
			
			fobj->ecm_obj.B1 = 1000000;
			fobj->ecm_obj.num_curves = 1;
			
			gettimeofday(&start, NULL);
			mpz_set(fobj->ecm_obj.gmp_n, n);
			mpz_set(fobj->N, n);
			curves_run = ecm_loop(fobj);
			gettimeofday(&stop, NULL);
			t_time = ytools_difftime(&start, &stop);
			
			printf("elapsed time per curve for B1=%" PRIu64 "k on c%zu = %6.4f seconds.\n",
				fobj->ecm_obj.B1 / 1000, mpz_sizeinbase(n, 10), t_time / curves_run);
			
			ecm_extraptime1[1][i] = t_time / curves_run;

			fobj->ecm_obj.B1 = 10000000;
			fobj->ecm_obj.num_curves = 1;

			gettimeofday(&start, NULL);
			mpz_set(fobj->ecm_obj.gmp_n, n);
			mpz_set(fobj->N, n);
			curves_run = ecm_loop(fobj);
			gettimeofday(&stop, NULL);
			t_time = ytools_difftime(&start, &stop);

			printf("elapsed time per curve for B1=%" PRIu64 "k on c%zu = %6.4f seconds.\n",
				fobj->ecm_obj.B1 / 1000, mpz_sizeinbase(n, 10), t_time / curves_run);

			ecm_extraptime1[2][i] = t_time / curves_run;
		}

		mpz_clear(base_100b);
		gmp_randclear(gmp_randstate);
		int j;

		double exp_consts[3];
		double exp_exps[3];
		for (j = 0; j < 3; j++)
		{
			double etime[NUM_ECM_PTS];
			for (i = 0; i < NUM_ECM_PTS; i++)
			{
				printf("%1.0f, %1.4f\n", ecm_sizes[i], ecm_extraptime1[j][i]);
				etime[i] = ecm_extraptime1[j][i];
			}

			fit = best_fit(ecm_sizes, etime, NUM_ECM_PTS, &a3, &b3, 2);
			exp_consts[j] = pow(BASE_e, b3);
			exp_exps[j] = a3;
			avg_exp += a3;
			printf("best linear fit over size is ln(y) = %g * x + %g\nR^2 = %g\n", a3, b3, fit);
			printf("best exponential fit over size is y = %g * exp(%g * x)\n", 
				exp_consts[j], exp_exps[j]);
		}

		// now do an average of the exponent and a linear fit to the constant and that's
		// the extrapolation as a function of size and B1.
		double b1vals[3] = {1e+5, 1e+6, 1e+7};
		fit = best_fit(b1vals, exp_consts, 3, &a3, &b3, 1);
		avg_exp /= 3.0;
		printf("best linear fit over B1 is y = %g * x + %g\nR^2 = %g\n", a3, b3, fit);
		printf("average size exponent over B1 is %g\n", avg_exp);

	}

	// AVX-ECM: radix-2^52 / RNS arithmetic crossover (also ecmtune())
	ecm_rns_tune(fobj);

done:

	//write the coefficients to the .ini file
	update_INI(pow(BASE_e, b), a, pow(BASE_e, b2), a2, xover,
		a3, b3, avg_exp, inobj->MEAS_CPU_FREQUENCY, inobj->CPU_ID_STR);


	//free_factobj(&fobj);
	//clear_factor_list(fobj);
	mpz_clear(n);
	return;
}


void make_job_file(char *sname, uint32_t *startq, uint32_t *qrange, char *inputstr, int inputnum, fact_obj_t *fobj)
{
	FILE *out;
	int siever;

	out = fopen("tune.job","w");
	if (out == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open tune.job for writing!\n");
		exit(1);
	}

	printf("nfs: commencing job file construction for c%d\n",(int)strlen(inputstr));

	fprintf(out, "n: %s\n",inputstr);
	switch (inputnum)
	{
	case 0:
		//85 digit info:
		fprintf(out, "type: gnfs\n");
		fprintf(out, "skew: 117904.22\n");
		fprintf(out, "Y0: -154288192969328825805\n");
		fprintf(out, "Y1: 86313418387\n");
		fprintf(out, "c0: -19278127997978752776944\n");
		fprintf(out, "c1: -10879961095353705327\n");
		fprintf(out, "c2: 123502072695578\n");
		fprintf(out, "c3: 1026425568\n");
		fprintf(out, "c4: 3312\n");
		fprintf(out, "rlim: 900000\n");
		fprintf(out, "alim: 900000\n");
		fprintf(out, "lpbr: 24\n");
		fprintf(out, "lpba: 24\n");
		fprintf(out, "mfbr: 48\n");
		fprintf(out, "mfba: 48\n");
		fprintf(out, "rlambda: 2.5\n");
		fprintf(out, "alambda: 2.5\n");
		siever = 11;
		*startq = 450000;
		*qrange = 2000;
		break;
	case 1:
		//90 digit info:
		fprintf(out, "type: gnfs\n");
		fprintf(out, "skew: 568902.33\n");
		fprintf(out, "Y0: -3993643540969129401980\n");
		fprintf(out, "Y1: 679422774953\n");
		fprintf(out, "c0: -55713309271074415620388179\n");
		fprintf(out, "c1: -274978535439659352141\n");
		fprintf(out, "c2: -10068765363614\n");
		fprintf(out, "c3: -20194136\n");
		fprintf(out, "c4: 1680\n");
		fprintf(out, "rlim: 1200000\n");
		fprintf(out, "alim: 1200000\n");
		fprintf(out, "lpbr: 25\n");
		fprintf(out, "lpba: 25\n");
		fprintf(out, "mfbr: 50\n");
		fprintf(out, "mfba: 50\n");
		fprintf(out, "rlambda: 2.5\n");
		fprintf(out, "alambda: 2.5\n");
		siever = 11;
		*startq = 600000;
		*qrange = 2000;
		break;
	case 2:
		//95 digit info:
		fprintf(out, "type: gnfs\n");
		fprintf(out, "skew: 1924450.27\n");
		fprintf(out, "Y0: -109555701149202546782193\n");
		fprintf(out, "Y1: 1831799378543\n");
		fprintf(out, "c0: -270645055551195008101213600\n");
		fprintf(out, "c1: 2685721124630647750266\n");
		fprintf(out, "c2: -1414126465881023\n");
		fprintf(out, "c3: 123061242\n");
		fprintf(out, "c4: 336\n");
		fprintf(out, "rlim: 1500000\n");
		fprintf(out, "alim: 1500000\n");
		fprintf(out, "lpbr: 25\n");
		fprintf(out, "lpba: 25\n");
		fprintf(out, "mfbr: 50\n");
		fprintf(out, "mfba: 50\n");
		fprintf(out, "rlambda: 2.5\n");
		fprintf(out, "alambda: 2.5\n");
		siever = 12;
		*startq = 750000;
		*qrange = 2000;
		break;
	case 3:
		//100 digit info:
		fprintf(out, "type: gnfs\n");
		fprintf(out, "skew: 3626061.93\n");
		fprintf(out, "Y0: -1635701959972760374016205\n");
		fprintf(out, "Y1: 81224154764317\n");
		fprintf(out, "c0: 17481065715537621465301512424\n");
		fprintf(out, "c1: 18163435376729925884930\n");
		fprintf(out, "c2: -3760568072527089\n");
		fprintf(out, "c3: -3386217900\n");
		fprintf(out, "c4: 252\n");
		fprintf(out, "rlim: 1800000\n");
		fprintf(out, "alim: 1800000\n");
		fprintf(out, "lpbr: 26\n");
		fprintf(out, "lpba: 26\n");
		fprintf(out, "mfbr: 52\n");
		fprintf(out, "mfba: 52\n");
		fprintf(out, "rlambda: 2.5\n");
		fprintf(out, "alambda: 2.5\n");
		siever = 12;
		*startq = 900000;
		*qrange = 5000;
		break;
	case 4:
		//105 digit info:
		fprintf(out, "type: gnfs\n");
		fprintf(out, "skew: 18508.57\n");
		fprintf(out, "Y0: -164285488596487844605\n");
		fprintf(out, "Y1: 96982771931\n");
		fprintf(out, "c0: 3745552009373958707561976\n");
		fprintf(out, "c1: -171516434496742996334\n");
		fprintf(out, "c2: -98669563931139619\n");
		fprintf(out, "c3: 246374333386\n");
		fprintf(out, "c4: 139871450\n");
		fprintf(out, "c5: 3900\n");
		fprintf(out, "rlim: 2500000\n");
		fprintf(out, "alim: 2500000\n");
		fprintf(out, "lpbr: 26\n");
		fprintf(out, "lpba: 26\n");
		fprintf(out, "mfbr: 52\n");
		fprintf(out, "mfba: 52\n");
		fprintf(out, "rlambda: 2.5\n");
		fprintf(out, "alambda: 2.5\n");
		siever = 12;
		*startq = 1250000;
		*qrange = 5000;
		break;
	case 5:
		//110 digit info:
		fprintf(out, "type: gnfs\n");
		fprintf(out, "skew: 32397.25\n");
		fprintf(out, "Y0: -1237038659848269396424\n");
		fprintf(out, "Y1: 551747168957\n");
		fprintf(out, "c0: -135941441485299623821377975\n");
		fprintf(out, "c1: 18442770087320538690649\n");
		fprintf(out, "c2: 1144416907608563521\n");
		fprintf(out, "c3: -42817211718353\n");
		fprintf(out, "c4: -347973258\n");
		fprintf(out, "c5: 16632\n");
		fprintf(out, "rlim: 3200000\n");
		fprintf(out, "alim: 3200000\n");
		fprintf(out, "lpbr: 26\n");
		fprintf(out, "lpba: 26\n");
		fprintf(out, "mfbr: 52\n");
		fprintf(out, "mfba: 52\n");
		fprintf(out, "rlambda: 2.5\n");
		fprintf(out, "alambda: 2.5\n");
		siever = 13;
		*startq = 1600000;
		*qrange = 5000;
		break;
	default:
		printf("unknown input in gnfs tuning\n");
		exit(-1);

	}

    check_siever(fobj, sname, siever);
	
	fclose(out);

	return;
}

void check_siever(fact_obj_t *fobj, char* sname, int siever)
{
    FILE* test;

    switch (siever)
    {
    case 11:
        sprintf(sname, "%sgnfs-lasieve4I11e", fobj->nfs_obj.ggnfs_dir);
        break;
    case 12:
        sprintf(sname, "%sgnfs-lasieve4I12e", fobj->nfs_obj.ggnfs_dir);
        break;
    case 13:
        sprintf(sname, "%sgnfs-lasieve4I13e", fobj->nfs_obj.ggnfs_dir);
        break;
    case 14:
        sprintf(sname, "%sgnfs-lasieve4I14e", fobj->nfs_obj.ggnfs_dir);
        break;
    case 15:
        sprintf(sname, "%sgnfs-lasieve4I15e", fobj->nfs_obj.ggnfs_dir);
        break;
    }

#if defined(WIN32)
    sprintf(sname, "%s.exe", sname);
#endif

    // test for existence of the siever
    test = fopen(sname, "rb");
    if (test == NULL)
    {
        printf("fopen error: %s\n", strerror(errno));
        printf("could not find %s, bailing\n", sname);
        exit(-1);
    }
	fclose(test);

    return;
}

void update_INI(double mult, double exponent, double mult2, 
    double exponent2, double xover, 
	double b1slope, double b1intercept, double size_exp,
	double cpu_freq, char* cpu_str)
{
	FILE *in, *out;
	int i,j;

	char str[GSTR_MAXSIZE];
	char str2[GSTR_MAXSIZE];
	char newline[GSTR_MAXSIZE];
	char *key, *ptr, cpustr[80], osstr[80];
	int len, found_entry = 0;

	in = fopen("yafu.ini","r");
	if (in == NULL)
	{
		printf("could not open yafu.ini for reading!\n");
		printf("creating yafu.ini...\n");
		in = fopen("yafu.ini", "w");
		if (in == NULL)
		{
			printf("could not open yafu.ini for writing either - bailing\n");
			exit(1);
		}
		in = freopen("yafu.ini", "r", in);
	}

	out = fopen("_tmp.ini","w");
	if (out == NULL)
	{
		printf("could not open _tmp.ini for writing!");
		exit(-1);
	}

	while (fgets(str,GSTR_MAXSIZE,in) != NULL)
	{
		//if first character is a % sign, just copy it
		if (str[0] == '%')
		{
			fputs(str, out);
			continue;
		}
		else if (strlen(str) <= 1)
		{
			fputs(str, out);
			continue;
		}
		else if (found_entry)
		{
			//already found the right entry, just copy the rest of the file
			fputs(str, out);
			continue;
		}

		//remember original line
		strcpy(str2, str);

		//if last character of line is newline, remove it
		do 
		{
			len = (int)strlen(str);
			if (str[len - 1] == 10)
				str[len - 1] = '\0';
			else if (str[len - 1] == 13)
				str[len - 1] = '\0';
			else
				break;
		} while (len > 0);

		//read keyword by looking for an equal sign
		key = strtok(str,"=");

		if (key == NULL)
		{
			//printf("Invalid line in yafu.ini, use Keyword=Value pairs"
			//	"See docfile.txt for valid keywords");
			// just write it out and keep going
			//fputs(str, out);
			fprintf(out, "%s\n", str);
			continue;
		}
		else if (strcmp(key,"tune_info") == 0)
		{
			//this should return the rest of the string after the first =
			ptr = strtok((char *)0,"=");

			//read up to the first comma - this is the cpu id string
			j=0;
			for (i=0; i<strlen(ptr); i++)
			{
				if (ptr[i] == 10) break;
				if (ptr[i] == 13) break;
				if (ptr[i] == ',') break;
				cpustr[j++] = ptr[i];
			}
			cpustr[j] = '\0';
			i++;

			//read up to the next comma - this is the OS string
			j=0;
			for ( ; i<strlen(ptr); i++)
			{
				if (ptr[i] == 10) break;
				if (ptr[i] == 13) break;
				if (ptr[i] == ',') break;
				osstr[j++] = ptr[i];
			}
			osstr[j] = '\0';

			printf("found OS = %s and CPU = %s in tune_info field\n",osstr, cpustr);


#if defined(_WIN64)
			if ((strcmp(cpustr, cpu_str) == 0) && (strcmp(osstr, "WIN64") == 0))
			{
				printf("Replacing tune_info entry for %s - %s\n",osstr,cpustr);
				found_entry = 1;
				sprintf(newline, "tune_info=%s,WIN64,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
                    cpu_str,mult,exponent,mult2,exponent2,xover,
					b1slope, b1intercept, size_exp, cpu_freq);
				fputs(newline, out);
			}
#elif defined(WIN32)
			if ((strcmp(cpustr, cpu_str) == 0) && (strcmp(osstr, "WIN32") == 0))
			{
				printf("Replacing tune_info entry for %s - %s\n",osstr,cpustr);
				found_entry = 1;
				sprintf(newline, "tune_info=%s,WIN32,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
					cpu_str, mult, exponent, mult2, exponent2, xover,
					b1slope, b1intercept, size_exp, cpu_freq);
				fputs(newline, out);
			}
#else 
			if ((strcmp(cpustr, cpu_str) == 0) && (strcmp(osstr, "LINUX64") == 0))
			{
				printf("Replacing tune_info entry for %s - %s\n",osstr,cpustr);
				found_entry = 1;
				sprintf(newline, "tune_info=%s,LINUX64,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
					cpu_str, mult, exponent, mult2, exponent2, xover,
					b1slope, b1intercept, size_exp, cpu_freq);
				fputs(newline, out);
			}
#endif
			else
			{
				//just write the original line
				fputs(str2, out);
			}
		}
		else
		{
			//just write the original line
			fputs(str2, out);
		}
	}

	if (!found_entry)
	{
#if defined(_WIN64)
		printf("Adding tune_info entry for WIN64 - %s\n", cpu_str);
		sprintf(newline, "\ntune_info=%s,WIN64,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
			cpu_str, mult, exponent, mult2, exponent2, xover,
			b1slope, b1intercept, size_exp, cpu_freq);
		fputs(newline, out);
#elif defined(WIN32)
		printf("Adding tune_info entry for WIN32 - %s\n", cpu_str);
		sprintf(newline, "\ntune_info=%s,WIN32,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
			cpu_str, mult, exponent, mult2, exponent2, xover,
			b1slope, b1intercept, size_exp, cpu_freq);
		fputs(newline, out);
#else 
		printf("Adding tune_info entry for LINUX64 - %s\n", cpu_str);
		sprintf(newline, "\ntune_info=%s,LINUX64,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
			cpu_str, mult, exponent, mult2, exponent2, xover,
			b1slope, b1intercept, size_exp, cpu_freq);
		fputs(newline, out);
#endif
	}

	fclose(in);
	fclose(out);

	// swap old with new
	remove("yafu.ini");
	rename("_tmp.ini", "yafu.ini");

	return;
}

double best_fit(double *x, double *y, int numpts, 
	double *slope, double *intercept, int type)
{
	// given vectors of x and y values, compute the best fit
	// to the data and return the slope and y-intercept of the line
	// return the R^2 value.
	// if type == 1, does a linear fit.
	// if type == 2, does an exponential fit (by taking the ln of the y-data)
	// follows: http://mathworld.wolfram.com/LeastSquaresFitting.html
	double avgX, avgY, varX, varY, cov, *yy;
	int i;

	yy = (double *)malloc((unsigned long)numpts * sizeof(double));

	for (i = 0; i < numpts; i++)
	{
		if (type == 1)
		{
			yy[i] = y[i];
		}
		else
		{
			yy[i] = log(y[i]);
		}
	}
		

	avgX = 0;
	for (i=0; i<numpts; i++)
		avgX += x[i];
	avgX = avgX / (double)numpts;

	avgY = 0;
	for (i=0; i<numpts; i++)
		avgY += yy[i];
	avgY = avgY / (double)numpts;

	varX = 0;
	for (i=0; i<numpts; i++)
		varX += (x[i] * x[i]);
	varX = varX - (double)numpts * avgX * avgX;
	varX /= (double)numpts;

	varY = 0;
	for (i=0; i<numpts; i++)
		varY += (yy[i] * yy[i]);
	varY = varY - (double)numpts * avgY * avgY;
	varY /= (double)numpts;

	cov = 0;
	for (i=0; i<numpts; i++)
		cov += (x[i] * yy[i]);
	cov = cov - (double)numpts * avgX * avgY;
	cov /= (double)numpts;

	//printf("average x = %g, average y = %g, varX = %g, varY = %g, cov = %g\n",
	//	avgX, avgY, varX, varY, cov);
	
	*slope = cov / varX;
	*intercept = avgY - *slope * avgX;

	free(yy);
	return cov * cov / (varX * varY);
}

//----------------------- SIQS PARAMETER TUNE ---------------------------------//
// siqstune(): pick this machine's SIQS parameters from a table of competing
// parameter sets (siqs_candidates.txt, exported by siqs_tune.py: several sets
// per digit size, each with the number of full and partial relations it
// needed to finish, measured on the test number n in the row).  Rank 0 is yafu
// as it runs without flags and is run that way here; every other row (rank 1:
// the built-in parameters on the bucket pipeline they were tuned for, with the
// formula TFSm; 2+: tuned sets) is forced exactly as listed, so the partials
// it needs reproduce.  At each digit size every set sieves each of those
// numbers for -siqsT seconds (default 5); the rate of partial relations after
// a short warmup extrapolates to a full sieving time (on the same n the
// partials needed reproduce to ~0.5%), plus the table's post-processing time
// scaled by the sieving speed seen here.  Tables without test numbers fall
// back to a number made here and the fulls rate.  Inputs small
// enough to finish within the sample time are simply timed.  The fastest set
// per size is written to the siqs params file (siqs_params.txt), which SIQS
// then uses for its defaults; without that file the built-in table is used.

#define SIQS_CANDIDATES_FILE "siqs_candidates.txt"
#define SIQSTUNE_DEFAULT_SECONDS 5
#define SIQSTUNE_MAX_FULL_RUNS 3
#define SIQSTUNE_MAX_NUMBERS 16

typedef struct
{
    int ok;
    int completed;          // finished inside the sample time: timed directly
    int runs;
    double est;             // estimated (or measured) total time
    double t_warm, t_end;   // the sample: warmup snapshot and stop
    uint32_t full_warm, full_end, partial_warm, partial_end;
    double fulls_need, partials_need;   // table values scaled to this input
    double est_fulls;       // sieving time extrapolated from fulls
    double est_partials;    // ... and from partials
    int by_partials;        // est used the partials (same input as the table)
} siqstune_meas_t;

// extrapolate a count that grows linearly after the warmup to 'need'
static double siqstune_extrap(double t_warm, double t_end, double c_warm,
    double c_end, double need)
{
    double rate;
    if ((t_end <= t_warm) || (c_end <= c_warm) || (need <= 0))
        return -1.0;
    rate = (c_end - c_warm) / (t_end - t_warm);
    return t_warm + (need - c_warm) / rate;
}

// a balanced semiprime with exactly 'digits' digits
static void siqstune_semiprime(mpz_t n, int digits, gmp_randstate_t st)
{
    mpz_t lo, t, p, q;

    mpz_init(lo);
    mpz_init(t);
    mpz_init(p);
    mpz_init(q);
    mpz_ui_pow_ui(lo, 10, (unsigned long)(digits - 1));

    do
    {
        // target t in [10^(d-1), 10^d); p a little below sqrt(t), q ~ t/p
        mpz_mul_ui(t, lo, 9);
        mpz_urandomm(t, st, t);
        mpz_add(t, t, lo);
        mpz_sqrt(p, t);
        mpz_tdiv_q_2exp(q, p, 2);
        mpz_urandomm(q, st, q);
        mpz_sub(p, p, q);
        mpz_nextprime(p, p);
        mpz_tdiv_q(q, t, p);
        mpz_nextprime(q, q);
        mpz_mul(n, p, q);
    } while ((gmp_base10(n) != digits) || (mpz_cmp(p, q) == 0));

    mpz_clear(lo);
    mpz_clear(t);
    mpz_clear(p);
    mpz_clear(q);
}

// force exactly the row's parameters (none for rank 0: yafu without flags)
static void siqstune_set_params(qs_obj_t* q, siqs_param_row_t* r)
{
    uint32_t has = (r->rank == 0) ? 0 : r->has;

    q->gbl_override_B_flag = (has & SIQSP_B) ? 1 : 0;
    q->gbl_override_B = r->B;
    q->gbl_override_lpmult_flag = (has & SIQSP_M) ? 1 : 0;
    q->gbl_override_lpmult = r->M;
    q->gbl_override_blocks_flag = (has & SIQSP_NB) ? 1 : 0;
    q->gbl_override_blocks = r->NB;
    q->gbl_override_tf_flag = (has & SIQSP_TF) ? 1 : 0;
    q->gbl_override_tf = r->TF;
    q->gbl_override_small_cutoff_flag = (has & SIQSP_TFSM) ? 1 : 0;
    q->gbl_override_small_cutoff = r->TFSm;
    q->no_small_cutoff_opt = 0;
    q->gbl_override_mfbd = (has & SIQSP_MFBD) ? r->MFBD : 0.0;
    q->gbl_override_mfbt = (has & SIQSP_MFBT) ? r->MFBT : 0.0;
    q->gbl_override_spv_limit_flag = (has & SIQSP_SPV) ? 1 : 0;
    q->gbl_override_spv_limit = r->SPV;
    q->gbl_siqs_mid = (has & SIQSP_MID) ? r->mid : SIQS_STRAT_AUTO;
    q->gbl_siqs_large = (has & SIQSP_LARGE) ? r->large : SIQS_STRAT_AUTO;
    q->gbl_siqs_kernel = (has & SIQSP_KERNEL) ? r->kernel : SIQS_STRAT_AUTO;
    q->gbl_override_lpb = 0;
    q->gbl_override_rel_flag = 0;
}

static int siqstune_measure(fact_obj_t* fobj, mpz_t n, siqs_param_row_t* r,
    int same_n, double seconds, siqstune_meas_t* m)
{
    qs_obj_t* q = &fobj->qs_obj;
    double times[SIQSTUNE_MAX_FULL_RUNS];
    double spent = 0;
    int i;

    memset(m, 0, sizeof(siqstune_meas_t));

    for (i = 0; i < SIQSTUNE_MAX_FULL_RUNS; i++)
    {
        double scale, est_sieve, post;

        siqstune_set_params(q, r);
        q->gbl_override_time_flag = 1;
        q->gbl_override_time = (uint32_t)ceil(seconds);
        q->tune_warmup = (seconds / 4 > 0.5) ? seconds / 4 : 0.5;
        q->tune_t_warm = q->tune_t_end = 0;
        q->tune_full_warm = q->tune_partial_warm = 0;
        q->tune_full_end = q->tune_partial_end = 0;
        q->tune_needed = 0;
        q->total_time = 0;

        reset_factobj(fobj);
        mpz_set(q->gmp_n, n);
        SIQS(fobj);
        m->runs++;

        if (q->total_time > 0)
        {
            // finished inside the sample: that's the measurement.  Repeat
            // while it's cheap, and report the median.
            times[i] = q->total_time;
            spent += q->total_time;
            m->completed = 1;
            if (spent + q->total_time > seconds)
            {
                i++;
                break;
            }
            continue;
        }

        // stopped at the time limit: extrapolate the rates measured after the
        // warmup to the counts this set needs to finish.  When the row was
        // measured on this same input (its n column), its partials count
        // reproduces to ~0.5% and, being ~20x the fulls count, gives the less
        // noisy rate: partials are the measure.  On any other input the
        // partials needed vary ~10% between inputs of a size but the fulls
        // only ~3%, so fall back to fulls, scaled by the relations needed
        // (B + extra) here vs. there, since the defaults' factor base size
        // varies within a digit size.
        if (q->tune_t_end <= 0)
            return 0;

        scale = ((r->needed > 0) && (q->tune_needed > 0)) ?
            (double)q->tune_needed / r->needed : 1.0;
        m->fulls_need = (r->fulls > 0) ? r->fulls * scale : -1.0;
        m->partials_need = (r->partials > 0) ? r->partials * scale : -1.0;
        m->t_warm = q->tune_t_warm;
        m->t_end = q->tune_t_end;
        m->full_warm = q->tune_full_warm;
        m->full_end = q->tune_full_end;
        m->partial_warm = q->tune_partial_warm;
        m->partial_end = q->tune_partial_end;
        m->est_fulls = siqstune_extrap(m->t_warm, m->t_end, m->full_warm,
            m->full_end, m->fulls_need);
        m->est_partials = siqstune_extrap(m->t_warm, m->t_end, m->partial_warm,
            m->partial_end, m->partials_need);

        m->by_partials = same_n && (m->est_partials > 0);
        if (m->by_partials)
            est_sieve = m->est_partials;
        else
            est_sieve = (m->est_fulls > 0) ? m->est_fulls : m->est_partials;
        if (est_sieve <= 0)
            return 0;
        if (est_sieve < q->tune_t_end)
            est_sieve = q->tune_t_end;

        // post-processing: the table's time, scaled by how fast this machine
        // sieves compared to the one that made the table
        if ((r->t_post >= 0) && (r->t_sieve > 0))
            post = r->t_post * est_sieve / r->t_sieve;
        else
            post = 0.03 * est_sieve;

        m->ok = 1;
        m->est = est_sieve + post;
        return 1;
    }

    // median of the completed runs
    {
        int j, k;
        for (j = 1; j < i; j++)
        {
            double t = times[j];
            for (k = j; (k > 0) && (times[k - 1] > t); k--)
                times[k] = times[k - 1];
            times[k] = t;
        }
        m->ok = 1;
        m->est = times[i / 2];
    }
    return 1;
}

void siqs_param_tune(fact_obj_t* fobj)
{
    qs_obj_t* q = &fobj->qs_obj;
    qs_obj_t saved;
    siqs_param_row_t* cand;
    siqs_param_row_t* out;
    siqstune_meas_t* meas;
    double* rank_sum;
    int* rank_cnt;
    int* rank_full;
    int* rank_row;
    int ncand, nout = 0, i, j, maxrank;
    int saved_vflag = fobj->VFLAG;
    double seconds;
    gmp_randstate_t st;
    mpz_t n;
    char buf[256];
    char comment[1024];

    ncand = siqs_read_param_table(SIQS_CANDIDATES_FILE, &cand);
    if (ncand == 0)
    {
        printf("siqstune: no parameter sets found in %s (export them with "
            "'siqs_tune.py export')\n", SIQS_CANDIDATES_FILE);
        return;
    }

    seconds = q->gbl_override_time_flag ? (double)q->gbl_override_time :
        (double)SIQSTUNE_DEFAULT_SECONDS;

    // the measurements need full control of the SIQS options; the user's
    // are put back at the end
    memcpy(&saved, q, sizeof(qs_obj_t));
    q->no_tuned_params = 1;
    q->inmem_cutoff = 1000;         // never touch (or resume from) a savefile
    if (fobj->VFLAG <= 0)
        fobj->VFLAG = -1;

    printf("siqstune: %d parameter sets from %s, %1.0f second samples, "
        "%d thread%s\n", ncand, SIQS_CANDIDATES_FILE, seconds, fobj->THREADS,
        fobj->THREADS > 1 ? "s" : "");

    out = (siqs_param_row_t*)malloc((unsigned long)ncand * sizeof(siqs_param_row_t));
    meas = (siqstune_meas_t*)malloc((unsigned long)ncand * sizeof(siqstune_meas_t));
    mpz_init(n);
    gmp_randinit_default(st);

    for (i = 0, maxrank = 0; i < ncand; i++)
        if (cand[i].rank > maxrank)
            maxrank = cand[i].rank;
    rank_sum = (double*)malloc((unsigned long)(maxrank + 1) * sizeof(double));
    rank_cnt = (int*)malloc((unsigned long)(maxrank + 1) * sizeof(int));
    rank_full = (int*)malloc((unsigned long)(maxrank + 1) * sizeof(int));
    rank_row = (int*)malloc((unsigned long)(maxrank + 1) * sizeof(int));

    for (i = 0; i < ncand; i++)
    {
        int d = cand[i].digits;
        int best = -1, ref, k, rk, ntests;
        const char* nums[SIQSTUNE_MAX_NUMBERS];
        int nnum = 0;

        // each digit size once, in file order
        for (j = 0; j < i; j++)
            if (cand[j].digits == d)
                break;
        if (j < i)
            continue;

        // the test numbers: the ones the table's rows were measured on.  A
        // table without them gets one made here, which depends only on its
        // size, so a size measures the same input whatever else is in the file.
        for (j = i; j < ncand; j++)
        {
            if ((cand[j].digits != d) || (cand[j].n[0] == '\0'))
                continue;
            for (k = 0; k < nnum; k++)
                if (strcmp(nums[k], cand[j].n) == 0)
                    break;
            if ((k == nnum) && (nnum < SIQSTUNE_MAX_NUMBERS))
                nums[nnum++] = cand[j].n;
        }
        ntests = (nnum > 0) ? nnum : 1;

        for (rk = 0; rk <= maxrank; rk++)
        {
            rank_sum[rk] = 0;
            rank_cnt[rk] = 0;
            rank_full[rk] = 0;
            rank_row[rk] = -1;
        }

        for (k = 0; k < ntests; k++)
        {
            const char* nstr = (nnum > 0) ? nums[k] : "";

            if (nnum > 0)
            {
                mpz_set_str(n, nstr, 10);
                gmp_printf("\n==== c%d: %Zd ====\n", d, n);
            }
            else
            {
                gmp_randseed_ui(st, (unsigned long)(20261022 + d));
                siqstune_semiprime(n, d, st);
                gmp_printf("\n==== c%d: %Zd ====\n(the table has no test "
                    "number for c%d: estimating from fulls)\n", d, n, d);
            }

            for (j = i; j < ncand; j++)
            {
                siqstune_meas_t* m = &meas[j];

                if ((cand[j].digits != d) || (strcmp(cand[j].n, nstr) != 0))
                    continue;

                siqs_format_params(&cand[j], buf, 256);
                if (cand[j].rank == 0)
                    printf("rank 0: (yafu without flags)\n");
                else if (cand[j].rank == 1)
                    printf("rank 1: (built-in params, bucket pipeline) %s\n", buf);
                else
                    printf("rank %d: %s\n", cand[j].rank, buf);
                fflush(stdout);

                if (!siqstune_measure(fobj, n, &cand[j], nnum > 0, seconds, m))
                {
                    printf("    no estimate (missing counts in table, or no "
                        "progress in the sample)\n");
                    continue;
                }

                if (m->completed)
                    printf("    finished in the sample: %1.3f s (median of %d)\n",
                        m->est, m->runs);
                else
                {
                    printf("    sample %1.2f-%1.2f s: partials %u -> %u of %1.0f, "
                        "fulls %u -> %u of %1.0f\n", m->t_warm, m->t_end,
                        m->partial_warm, m->partial_end, m->partials_need,
                        m->full_warm, m->full_end, m->fulls_need);
                    printf("    sieving est. %1.2f s from partials%s, %1.2f s from "
                        "fulls%s -> estimated %1.2f s total\n",
                        m->est_partials, m->by_partials ? " (used)" : "",
                        m->est_fulls, m->by_partials ? "" : " (used)", m->est);
                }
                fflush(stdout);

                rk = cand[j].rank;
                rank_sum[rk] += m->est;
                rank_cnt[rk]++;
                if (m->completed)
                    rank_full[rk]++;
                if (rank_row[rk] < 0)
                    rank_row[rk] = j;
            }
        }

        // the fastest set on average over the test numbers, among those
        // measured on all of them.  Rank 0 (no flags) runs yafu's adaptive
        // TFSm tuning, which tries TFSm +-5 during the first polynomials --
        // a large part of a short sample -- and settles differently from run
        // to run, so its extrapolated estimate is unreliable (17% low in a
        // c70 test).  It only competes where it was timed outright, or when
        // nothing else was measured; rank 1, the same built-in parameters
        // forced, is its reproducible stand-in.
        for (rk = 1; rk <= maxrank; rk++)
        {
            if ((rank_cnt[rk] == ntests) && ((best < 0) ||
                (rank_sum[rk] < rank_sum[best])))
                best = rk;
        }
        if ((rank_cnt[0] == ntests) && ((best < 0) ||
            ((rank_full[0] == ntests) && (rank_sum[0] < rank_sum[best]))))
            best = 0;
        if ((rank_cnt[0] == ntests) && (rank_full[0] < ntests))
            printf("(rank 0's estimate is for reference only: adaptive TFSm "
                "makes it unreliable)\n");

        if (best < 0)
        {
            printf("c%d: no usable measurements\n", d);
            continue;
        }

        // ratios are relative to the built-in parameters: rank 1 where
        // measured, else rank 0
        ref = (rank_cnt[1] == ntests) ? 1 : ((rank_cnt[0] == ntests) ? 0 : -1);
        out[nout] = cand[rank_row[best]];
        out[nout].est_time = rank_sum[best] / ntests;
        out[nout].ratio = (ref >= 0) ? rank_sum[best] / rank_sum[ref] : -1.0;
        siqs_format_params(&out[nout], buf, 256);
        printf("c%d best: rank %d, %1.2f s", d, best, out[nout].est_time);
        if (ntests > 1)
            printf(" (mean of %d numbers)", ntests);
        if (ref >= 0)
            printf(" (%1.3f x rank %d)", out[nout].ratio, ref);
        printf(": %s%s\n", (best == 0) ? "(yafu without flags)" :
            (best == 1) ? "(built-in params, bucket pipeline) " : "",
            (best == 0) ? "" : buf);
        nout++;
    }

    free(rank_sum);
    free(rank_cnt);
    free(rank_full);
    free(rank_row);

    fobj->VFLAG = saved_vflag;

    if (nout > 0)
    {
        snprintf(comment, sizeof(comment),
            "# siqstune: %s, %d thread%s, %1.0f s samples, from %s\n"
            "# one row per digit size: the fastest parameter set here.  rank 0 "
            "rows mean yafu without\n# flags won; SIQS then uses its built-in "
            "table at that size.  rank 1 = the built-in\n# params on the "
            "bucket pipeline, 2+ = tuned sets.  est_time = estimated total "
            "seconds,\n# ratio = est_time / rank 1's (else rank 0's) est_time.\n",
            fobj->CPU_ID_STR, fobj->THREADS, fobj->THREADS > 1 ? "s" : "",
            seconds, SIQS_CANDIDATES_FILE);

        if (siqs_write_param_table(saved.siqs_params_file, out, nout, comment))
            printf("\nsiqstune: wrote %d rows to %s; SIQS will use them for its "
                "default parameters\n", nout, saved.siqs_params_file);
        else
            printf("\nsiqstune: could not write %s\n", saved.siqs_params_file);

        // the same winners as rows for the built-in tables (siqs_params.c),
        // for the table matching this cpu's kernel ISA and the rows'
        // large-prime strategy/kernel.  0 = not tabulated: the formula/default.
        printf("\nsiqstune: as built-in table rows (siqs_params.c; "
            "bits, B, M, NB, TF, TFSm, MFBD, MFBT), %s:\n", fobj->CPU_ID_STR);
        for (i = 0; i < nout; i++)
        {
            siqs_format_builtin_row(&out[i], buf, 256);
            printf("    %s\n", buf);
        }
    }

    // put the user's options back
    q->gbl_override_B_flag = saved.gbl_override_B_flag;
    q->gbl_override_B = saved.gbl_override_B;
    q->gbl_override_lpmult_flag = saved.gbl_override_lpmult_flag;
    q->gbl_override_lpmult = saved.gbl_override_lpmult;
    q->gbl_override_blocks_flag = saved.gbl_override_blocks_flag;
    q->gbl_override_blocks = saved.gbl_override_blocks;
    q->gbl_override_tf_flag = saved.gbl_override_tf_flag;
    q->gbl_override_tf = saved.gbl_override_tf;
    q->gbl_override_small_cutoff_flag = saved.gbl_override_small_cutoff_flag;
    q->gbl_override_small_cutoff = saved.gbl_override_small_cutoff;
    q->no_small_cutoff_opt = saved.no_small_cutoff_opt;
    q->gbl_override_mfbd = saved.gbl_override_mfbd;
    q->gbl_override_mfbt = saved.gbl_override_mfbt;
    q->gbl_override_spv_limit_flag = saved.gbl_override_spv_limit_flag;
    q->gbl_override_spv_limit = saved.gbl_override_spv_limit;
    q->gbl_siqs_mid = saved.gbl_siqs_mid;
    q->gbl_siqs_large = saved.gbl_siqs_large;
    q->gbl_siqs_kernel = saved.gbl_siqs_kernel;
    q->gbl_override_lpb = saved.gbl_override_lpb;
    q->gbl_override_rel_flag = saved.gbl_override_rel_flag;
    q->gbl_override_time_flag = saved.gbl_override_time_flag;
    q->gbl_override_time = saved.gbl_override_time;
    q->inmem_cutoff = saved.inmem_cutoff;
    q->no_tuned_params = saved.no_tuned_params;
    q->tune_warmup = 0;

    gmp_randclear(st);
    mpz_clear(n);
    free(meas);
    free(out);
    free(cand);
    return;
}

// ---------------------------------------------------------------------------
// ecmtune(): AVX-ECM radix-2^52 / RNS crossover for this machine
// ---------------------------------------------------------------------------

// replace (or append) a "key=value" line in yafu.ini
UNUSED_FUNC static void update_INI_key(const char *key, const char *value)
{
    FILE *in, *out;
    char str[GSTR_MAXSIZE];
    size_t klen = strlen(key);
    int found = 0;

    in = fopen("yafu.ini", "r");
    out = fopen("_tmp.ini", "w");
    if (out == NULL)
    {
        printf("could not open _tmp.ini for writing!\n");
        if (in != NULL)
            fclose(in);
        return;
    }
    if (in != NULL)
    {
        while (fgets(str, GSTR_MAXSIZE, in) != NULL)
        {
            if ((strncmp(str, key, klen) == 0) && (str[klen] == '='))
            {
                if (!found)
                    fprintf(out, "%s=%s\n", key, value);
                found = 1;
            }
            else
            {
                fputs(str, out);
            }
        }
        fclose(in);
    }
    if (!found)
        fprintf(out, "%s=%s\n", key, value);
    fclose(out);

    remove("yafu.ini");
    rename("_tmp.ini", "yafu.ini");
    printf("wrote %s=%s to yafu.ini\n", key, value);
}

#if defined(USE_AVX512F) || defined(USE_AVX2)

// seconds for one 8-curve AVX-ECM batch on n (1 thread), rns_bits as given
static double ecmtune_time(fact_obj_t *fobj, mpz_t n, int rns_bits,
    uint64_t B1, uint64_t sigma)
{
    struct timeval start, stop;
    uint32_t curves_run = 0;
    int numfactors = 0;
    double best = 1e30;
    int r;

    for (r = 0; r < 2; r++)
    {
        mpz_set(fobj->ecm_obj.gmp_n, n);
        fobj->ecm_obj.avx_rns_bits = rns_bits;
        fobj->ecm_obj.sigma = sigma;
        gettimeofday(&start, NULL);
        vec_ecm_main(fobj, 8, B1, 0, 1, &numfactors, -1, 0, &curves_run);
        gettimeofday(&stop, NULL);
        double t = ytools_difftime(&start, &stop);
        if (t < best)
            best = t;
    }
    return best;
}

void ecm_rns_tune(fact_obj_t *fobj)
{
    // inputs: products of 256-bit primes (no curve can find a factor at
    // this B1), sizes 768..4096 bits; time radix-2^52 vs RNS per size
    const uint64_t B1 = 20000;
    const uint64_t sigma = 1234567;
    int sizes[32], ns = 0, bits, i;
    double tr[32], tn[32];
    uint64_t saved_sigma = fobj->ecm_obj.sigma;
    int saved_vflag = fobj->VFLAG;
    int saved_threads = fobj->THREADS;
    double xover = 0;
    char val[64];
    gmp_randstate_t rs;
    mpz_t n, p;

    printf("ecmtune: AVX-ECM radix-2^52 vs RNS arithmetic, 8 curves, B1 = %" PRIu64
        ", 1 thread\n", B1);
    printf("  bits |   radix (s)     RNS (s) |  RNS/radix\n");
    if (!avxecm_rns_available())
    {
        printf("RNS kernels are not in this build: disabling RNS\n");
        update_INI_key("avxecm_rns_bits", "0");
        return;
    }

    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 42);
    mpz_init(n);
    mpz_init(p);
    fobj->VFLAG = -1;
    fobj->THREADS = 1;

    for (bits = 768; bits <= 4096; bits += 256)
        sizes[ns++] = bits;

    for (i = 0; i < ns; i++)
    {
        // n = product of 256-bit primes and a last prime of 257..512 bits,
        // exactly sizes[i] bits: no small factor (it would end the curves
        // early, and RNS needs N coprime to its composite moduli)
        do
        {
            mpz_set_ui(n, 1);
            while ((int)mpz_sizeinbase(n, 2) + 512 <= sizes[i])
            {
                mpz_urandomb(p, rs, 256);
                mpz_setbit(p, 255);
                mpz_nextprime(p, p);
                mpz_mul(n, n, p);
            }
            mpz_urandomb(p, rs, (mp_bitcnt_t)(unsigned)(sizes[i] - (int)mpz_sizeinbase(n, 2) + 1));
            mpz_setbit(p, (mp_bitcnt_t)(unsigned)(sizes[i] - (int)mpz_sizeinbase(n, 2)));
            mpz_nextprime(p, p);
            mpz_mul(n, n, p);
        } while ((int)mpz_sizeinbase(n, 2) != sizes[i]);

        tr[i] = ecmtune_time(fobj, n, 0, B1, sigma);
        {
            // is RNS possible for this n at all? (else a forced run is radix)
            uint32_t w;
            void *fx = avxecm_rns_fixed(n, &w);
            if (fx != NULL)
            {
                avxecm_rns_free_fixed(fx);
                tn[i] = ecmtune_time(fobj, n, 1, B1, sigma);
                printf("  %4d | %10.4f  %10.4f |  %8.3f\n", sizes[i], tr[i], tn[i], tn[i] / tr[i]);
            }
            else
            {
                tn[i] = 1e30;
                printf("  %4d | %10.4f         n/a |       n/a\n", sizes[i], tr[i]);
            }
        }
        fflush(stdout);
    }

    // crossover: the RNS/radix ratio falls below 1 and stays there; linear
    // interpolation of the ratio between the bracketing sizes
    for (i = ns - 1; i >= 0; i--)
    {
        if (tn[i] >= tr[i])
            break;
    }
    if (i == ns - 1)
    {
        // RNS is not faster at the largest measured size (RNS stops near
        // 4100 bits anyway): never use it
        xover = 0;
        printf("RNS is not faster at any size up to %d bits: disabling it\n", sizes[ns - 1]);
    }
    else if (i < 0)
    {
        xover = sizes[0];
        printf("RNS is faster at every size measured: crossover at %d bits\n", sizes[0]);
    }
    else
    {
        double r0 = tn[i] / tr[i], r1 = tn[i + 1] / tr[i + 1];
        xover = sizes[i] + (sizes[i + 1] - sizes[i]) * (r0 - 1.0) / (r0 - r1);
        printf("RNS is faster from about %1.0f bits\n", xover);
    }

    sprintf(val, "%d", (int)(xover + 0.5));
    update_INI_key("avxecm_rns_bits", val);

    fobj->ecm_obj.avx_rns_bits = (int)(xover + 0.5);
    fobj->ecm_obj.sigma = saved_sigma;
    fobj->VFLAG = saved_vflag;
    fobj->THREADS = saved_threads;
    mpz_set_ui(fobj->ecm_obj.gmp_n, 1);
    mpz_clear(n);
    mpz_clear(p);
    gmp_randclear(rs);
}

#else

void ecm_rns_tune(fact_obj_t *fobj)
{
    printf("ecmtune: AVX-ECM needs an AVX-512 or AVX2 build; nothing to tune\n");
}

#endif

// ---------------------------------------------------------------------------
// bayestune(): extended tune for -plan bayesian (docs/ecm_bayes.md), kept
// separate from tune() so that users who don't want the Bayesian pretest
// plan can skip it.  It writes bayes_tune.txt:
//   ECM: wall-clock seconds per curve with all threads, through ecm_loop()
//   (the same engine, curve family, stage 2 and threading that factor() will
//   use), for input sizes C60..C310 and B1 = 50k..3M, plus the 1-thread vs
//   N-thread efficiency at one point.
//   SIQS: whole-job wall time on random balanced semiprimes (two from C85),
//   C60..C95, in memory; written as a log-linear fit from C70.
//   NFS: test sieving of reference polys C120-C170 (bt_nfs_refs): three
//   blocks of special-q across the range each job needs, integrated to a
//   sieve time and scaled to a whole job.  bayestune(what) runs parts only
//   (1 ECM, 2 SIQS, 4 NFS, 8 NFS plus NFS@Home-parameter comparison) and
//   keeps the other parts' lines; "nfs" lines can also come from
//   docs/nfs_job_data/fit_sieve_times.py.
// bayes_tune_load() / bayes_ecm_sec_per_curve() / bayes_sieve_sec() read it.
// ---------------------------------------------------------------------------

#define BAYES_TUNE_FILE "bayes_tune.txt"
#define BT_MAXSIZES 64
#define BT_MAXB1 8

static struct
{
    int loaded;
    int threads;
    char cpu[80];
    int nsizes;
    double digits[BT_MAXSIZES];
    int nb1[BT_MAXSIZES];
    double b1[BT_MAXSIZES][BT_MAXB1];
    double spc[BT_MAXSIZES][BT_MAXB1];
    // whole-job sieve times: siqs (0), nfs (1), stored as thread-seconds
    // (wall s x the thread count the line was measured with)
    int nsieve[2];
    double sdig[2][BT_MAXSIZES];
    double ssec[2][BT_MAXSIZES];
    int sthr[2][BT_MAXSIZES];
} bt_ecm;

// n = p * q with p, q near digits/2 digits, n exactly 'digits' digits
static void bt_semiprime(mpz_t n, int digits, gmp_randstate_t rs)
{
    mpz_t p, q, t;
    mpz_init(p);
    mpz_init(q);
    mpz_init(t);
    do
    {
        int pbits = (int)((double)digits * 3.321928 / 2.0);
        mpz_urandomb(p, rs, (mp_bitcnt_t)pbits);
        mpz_setbit(p, (mp_bitcnt_t)(pbits - 1));
        mpz_nextprime(p, p);

        // t uniform in [1.5, 9] * 10^(digits-1)
        mpz_ui_pow_ui(t, 10, (unsigned long)(digits - 1));
        mpz_mul_ui(t, t, 15 + (unsigned long)gmp_urandomm_ui(rs, 75));
        mpz_tdiv_q_ui(t, t, 10);
        mpz_tdiv_q(q, t, p);
        mpz_nextprime(q, q);
        mpz_mul(n, p, q);
    } while (gmp_base10(n) != digits);
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(t);
}

// wall seconds per curve for one ecm_loop() call of ncurves curves at B1;
// *found = 1 if a factor turned up (the time is then not usable)
static double bt_ecm_run(fact_obj_t *fobj, mpz_t n, uint64_t B1,
    uint32_t ncurves, int *curves, double *B2, int *found)
{
    struct timeval start, stop;
    int c;
    double t;

    reset_factobj(fobj);
    mpz_set(fobj->N, n);
    mpz_set(fobj->ecm_obj.gmp_n, n);
    fobj->ecm_obj.B1 = B1;
    fobj->ecm_obj.num_curves = ncurves;
    fobj->ecm_obj.last_cp.B2 = 0;

    gettimeofday(&start, NULL);
    c = ecm_loop(fobj);
    gettimeofday(&stop, NULL);
    t = ytools_difftime(&start, &stop);

    *found = (mpz_cmp(fobj->ecm_obj.gmp_n, n) != 0);
    *curves = c;
    *B2 = fobj->ecm_obj.last_cp.B2;
    return (c > 0) ? t / (double)c : 0.0;
}

// time one (size, B1) point, drawing a new input if a curve finds a factor
static double bt_ecm_point(fact_obj_t *fobj, mpz_t n, int digits,
    uint64_t B1, uint32_t ncurves, gmp_randstate_t rs, int *curves,
    double *B2)
{
    int tries, found = 0;
    double spc = 0.0;
    for (tries = 0; tries < 5; tries++)
    {
        spc = bt_ecm_run(fobj, n, B1, ncurves, curves, B2, &found);
        if (!found)
            return spc;
        bt_semiprime(n, digits, rs);
    }
    return -1.0;
}

// what bayestune(what) measures (bit mask; bayestune() = ECM | SIQS | NFS)
#define BT_DO_ECM 1
#define BT_DO_SIQS 2
#define BT_DO_NFS 4
#define BT_DO_NFS_COMPARE 8

// whole GNFS job time per second of sieving: poly select ~7% of the job,
// filtering + LA + sqrt ~33% of the sieve time (complete C120/C130/C140 yafu
// runs on an i5-11600K, 6 threads: 1.48, 1.36, 1.46)
#define BT_NFS_TOTAL_PER_SIEVE 1.43
// target wall time of each thread's test-sieve block
#define BT_NFS_BLOCK_SEC 120.0

// Reference GNFS jobs for the NFS test sieving: a typical polynomial per
// size (C120-C140: yafu's own poly select on the local calibration runs;
// C150+: the median-Murphy-E NFS@Home poly of docs/nfs_job_data/polys), plus
// for C150+ the parameters NFS@Home sieved with and the median raw relation
// count of NFS@Home jobs with those parameters at that size (jobs.csv), for
// the comparison mode.
typedef struct
{
    int digits;
    const char *name;
    const char *poly;
    const char *own;        // NFS@Home parameter lines (NULL: none)
    int own_siever;
    uint32_t own_q0;
    double own_rels;
} bt_nfs_ref_t;

static const bt_nfs_ref_t bt_nfs_refs[] = {
    { 120, "local C120 run",
      "n: 326026769054469373585656403973021038065594659042330457601126530142220255457270935289475287783176236229586981677615616841\n"
      "skew: 395839.98\n"
      "c0: 213424040346523771560021371285\n"
      "c1: 3382706438461988268171541\n"
      "c2: 52509510896041622473\n"
      "c3: -46114004981635\n"
      "c4: -504885636\n"
      "c5: 156\n"
      "Y0: -291310524262297515725154\n"
      "Y1: 341740905217679\n", NULL, 0, 0, 0 },
    { 130, "local C130 run",
      "n: 1793745623039204287496874646909253275177942377082382889603744117550679912273343992708689975487443314894743928421441817425415799161\n"
      "skew: 1794046.19\n"
      "c0: -214170324545115532473123989784192\n"
      "c1: 2362969886060313593941759804\n"
      "c2: -1413744886029776133688\n"
      "c3: -1410085693472259\n"
      "c4: 338170380\n"
      "c5: 180\n"
      "Y0: -25099388576195464055193145\n"
      "Y1: 5316472426412777\n", NULL, 0, 0, 0 },
    { 140, "local C140 run",
      "n: 11667371413438610779651027524046401627437336389434906334441488004602840071783461929613705043678146111557846512173371540646602509533990141837\n"
      "skew: 2591213.05\n"
      "c0: 13659687464312813382518506327484240\n"
      "c1: 68389010700344902433972729920\n"
      "c2: 115351739183662847273118\n"
      "c3: 38226763455403805\n"
      "c4: -17582733488\n"
      "c5: 1200\n"
      "Y0: -1576107201356115944555987331\n"
      "Y1: 33926256534213953\n", NULL, 0, 0, 0 },
    { 150, "NFS@Home C150_829332_3736",
      "n: 177732960376900606849058002878068392783762464885373262455240888929003248313962586024387471460430146330150149604688127330913894732763747367564849904673\n"
      "skew: 6935452.02\n"
      "c0: -23341990681687169106151987572257851485\n"
      "c1: 34681493022024207656054533609801\n"
      "c2: -804632905643326178855684\n"
      "c3: -3624274055195821611\n"
      "c4: -26294732501\n"
      "c5: 10200\n"
      "Y0: -111746405762145756395883594856\n"
      "Y1: 1343118053814467\n",
      "rlim: 67000000\nalim: 67000000\nlpbr: 29\nlpba: 29\nmfbr: 58\nmfba: 58\n"
      "rlambda: 2.6\nalambda: 2.6\n", 14, 22000000, 64.9e6 },
    { 160, "NFS@Home C160_79_97_30",
      "n: 7147899195279076005126628539754280064346946453412379189519824669359603847180966162236261056367430639886117774702656241878940665411756142451546604262530551006287\n"
      "skew: 97406052.37\n"
      "c0: 104777786762760223299872473905044429299760\n"
      "c1: 16875649516643952940668249495488548\n"
      "c2: -24127128300196636843374592\n"
      "c3: -15496623762132245911\n"
      "c4: 8270171990\n"
      "c5: 120\n"
      "Y0: -35892145829078601856632218857329\n"
      "Y1: 66002327228495723\n",
      "rlim: 67000000\nalim: 67000000\nlpbr: 30\nlpba: 30\nmfbr: 60\nmfba: 60\n"
      "rlambda: 2.6\nalambda: 2.6\n", 14, 18000000, 130.5e6 },
    { 170, "NFS@Home C170_2381080_10275",
      "n: 14025988976053049132159502332472688247751009029636890956825855204662666432255313357443754674215769618001836605564383554679441797813037930016520085017530420942230022367911\n"
      "skew: 67841913.93\n"
      "c0: 12180003646158940602804670732376704254985385\n"
      "c1: 637176744251647249541970539059383766\n"
      "c2: -13644590466155451804249707265\n"
      "c3: -625708083069497050582\n"
      "c4: 4662579178504\n"
      "c5: 18480\n"
      "Y0: -946338041270337441832323675036756\n"
      "Y1: 297290394217753127\n",
      "rlim: 67000000\nalim: 67000000\nlpbr: 30\nlpba: 30\nmfbr: 60\nmfba: 60\n"
      "rlambda: 2.5\nalambda: 2.5\n", 15, 20000000, 126.5e6 },
};

// Relations yafu ends up gathering, as a multiple of its min_rels, where
// min_rels comes from nfs_set_min_rels()'s formula (table rows from C150 have
// none).  From the NFS@Home filtering brackets (docs/nfs_job_data/summary.md
// section 1): unique need / formula = ~0.95 at lpb 29, ~0.84 at lpb 30, ~0.77
// at lpb 31-32; raw = unique / (1 - ~17% duplicates); yafu filters at min_rels
// and adds 5% after each failed attempt.
static double bt_nfs_need_corr(int digits, int lpb)
{
    if (digits < 149)
        return 1.0;     // table min_rels, matched by the local runs
    if (lpb <= 29)
        return 1.10;
    if (lpb == 30)
        return 1.05;
    return 1.0;
}

// Sieve time from q0 until 'need' relations: yield y (rels per q) and wall
// time s (seconds per q, all threads) linear in q between the measured
// points, constant outside them.
static double bt_nfs_integrate(const double *qs, const double *y, const double *s,
    int np, double q0, double need, double *qend)
{
    double q = q0, rels = 0, t = 0, dq;
    int it = 0;

    dq = need / y[0] / 4000.0;
    if (dq < 1.0)
        dq = 1.0;
    while ((rels < need) && (it++ < 1000000))
    {
        double qm = q + dq / 2, yy, ss;
        int k;
        if (qm <= qs[0])
        {
            yy = y[0];
            ss = s[0];
        }
        else if (qm >= qs[np - 1])
        {
            yy = y[np - 1];
            ss = s[np - 1];
        }
        else
        {
            for (k = 1; k < np - 1; k++)
                if (qm < qs[k])
                    break;
            yy = y[k - 1] + (y[k] - y[k - 1]) * (qm - qs[k - 1]) / (qs[k] - qs[k - 1]);
            ss = s[k - 1] + (s[k] - s[k - 1]) * (qm - qs[k - 1]) / (qs[k] - qs[k - 1]);
        }
        if (yy <= 0)
            break;
        rels += yy * dq;
        t += ss * dq;
        q += dq;
    }
    *qend = q;
    return t;
}

// Test-sieve one reference job (own = 0: yafu's parameters; 1: the NFS@Home
// job's).  All threads sieve at once, each its own block at its own special-q
// spread over the range the job needs, so that each siever's time and yield
// is an independent sample at its q under a full machine load (and no thread
// waits for another at the end of a block).  Writes comment lines; returns
// the whole-job time estimate in seconds, or -1.
static double bt_nfs_ref(fact_obj_t *fobj, const bt_nfs_ref_t *r, int own,
    int threads, FILE *out, struct timeval *tstart)
{
    char job[4096];
    char desc[256];
    nfs_bayes_info_t info;
    double need, corr = 1.0, q0, startup, secq, y0, R, S, sieve, qend, total;
    double *qs, *y, *s, *secs;
    uint32_t *qv, *rels;
    int i, k, lpb, npts, rounds, T = (threads > 0) ? threads : 1;
    struct timeval tnow;

    snprintf(job, sizeof(job), "%s%s", r->poly, own ? r->own : "");
    if (nfs_bayes_setup(fobj, job, own ? r->own_siever : 0, &info) != 0)
        return -1.0;

    lpb = (info.lpbr > info.lpba) ? info.lpbr : info.lpba;
    if (own)
    {
        need = r->own_rels;
        q0 = r->own_q0;
        snprintf(desc, sizeof(desc), "need %1.1fM raw (median of NFS@Home jobs "
            "with these parameters)", need / 1e6);
    }
    else
    {
        corr = bt_nfs_need_corr(r->digits, lpb);
        need = info.min_rels * corr;
        q0 = info.startq;
        snprintf(desc, sizeof(desc), "need %1.1fM raw (min_rels %1.1fM x %1.2f)",
            need / 1e6, info.min_rels / 1e6, corr);
    }

    // at least 3 sample points; with fewer threads, several rounds
    rounds = (T >= 3) ? 1 : (3 + T - 1) / T;
    npts = rounds * T;
    qs = (double *)malloc((size_t)npts * sizeof(double));
    y = (double *)malloc((size_t)npts * sizeof(double));
    s = (double *)malloc((size_t)npts * sizeof(double));
    secs = (double *)malloc((size_t)T * sizeof(double));
    qv = (uint32_t *)malloc((size_t)T * sizeof(uint32_t));
    rels = (uint32_t *)malloc((size_t)T * sizeof(uint32_t));

    // siever start-up (factor base) time: one special-q per thread
    for (i = 0; i < T; i++)
        qv[i] = (uint32_t)q0 + (uint32_t)i;
    nfs_bayes_spread(fobj, qv, 1, secs, rels);
    for (startup = 0, i = 0; i < T; i++)
        startup += secs[i] / T;

    // pilot near the start: 100 special-q per thread, to size the blocks
    // and estimate the span of special-q the job needs
    for (i = 0; i < T; i++)
        qv[i] = (uint32_t)q0 + 100 * (uint32_t)i;
    nfs_bayes_spread(fobj, qv, 100, secs, rels);
    for (secq = 0, y0 = 0, i = 0; i < T; i++)
    {
        secq += ((secs[i] - startup > 0.2 * secs[i]) ? secs[i] - startup : 0.2 * secs[i]) /
            (100.0 * T);
        y0 += rels[i] / (100.0 * T);
    }
    if (y0 <= 0)
    {
        nfs_bayes_done(fobj);
        free(qs); free(y); free(s); free(secs); free(qv); free(rels);
        return -1.0;
    }
    R = 10.0 * floor(BT_NFS_BLOCK_SEC / secq / 10.0 + 0.5);
    if (R < 50) R = 50;
    if (R > 200000) R = 200000;
    // yields fall with q, so the start's yield underestimates the span:
    // sample to 1.2x it (the integration holds the last point beyond)
    S = 1.2 * need / y0;

    printf("bayestune: C%d %s, %s params: lim %1.1fM/%1.1fM lpb %d/%d mfb %d/%d %de, %s\n",
        r->digits, r->name, own ? "NFS@Home" : "yafu", info.rlim / 1e6, info.alim / 1e6,
        info.lpbr, info.lpba, info.mfbr, info.mfba, info.siever, desc);
    fprintf(out, "# nfs%s C%d %s, %s params: lim %1.1fM/%1.1fM lpb %d/%d mfb %d/%d "
        "lambda %1.2f/%1.2f %de side %c, %s; siever start-up %1.1f s; %d blocks of "
        "%1.0f q, one per thread\n",
        own ? " compare" : "", r->digits, r->name, own ? "NFS@Home" : "yafu",
        info.rlim / 1e6, info.alim / 1e6, info.lpbr, info.lpba, info.mfbr, info.mfba,
        info.rlambda, info.alambda, info.siever, info.side, desc, startup, npts, R);

    for (k = 0; k < rounds; k++)
    {
        for (i = 0; i < T; i++)
        {
            int p = k * T + i;
            qs[p] = q0 + ((npts > 1) ? S * p / (npts - 1) : 0.0);
            qv[i] = (uint32_t)qs[p];
        }
        nfs_bayes_spread(fobj, qv, (uint32_t)R, secs, rels);
        for (i = 0; i < T; i++)
        {
            int p = k * T + i;
            // that siever's own rate; s = wall seconds per q with all T
            // threads sieving at that rate
            double busy = (secs[i] - startup > 0.5 * secs[i]) ? secs[i] - startup : 0.5 * secs[i];
            y[p] = rels[i] / R;
            s[p] = busy / R / T;
            fprintf(out, "#   q %1.3fM: %1.3f rels/q, %1.1f rels/s per thread "
                "(%u rels in %1.1f s)\n", qs[p] / 1e6, y[p], rels[i] / busy, rels[i], secs[i]);
        }
        gettimeofday(&tnow, NULL);
        for (i = 0; i < T; i++)
        {
            int p = k * T + i;
            printf("   q %8.3fM: %5.2f rels/q, %6.1f rels/s per thread\n", qs[p] / 1e6, y[p],
                y[p] / s[p] / T);
        }
        printf("   (%1.0f s)\n", ytools_difftime(tstart, &tnow));
        fflush(stdout);
        fflush(out);
    }
    nfs_bayes_done(fobj);

    sieve = bt_nfs_integrate(qs, y, s, npts, q0, need, &qend);
    total = sieve * BT_NFS_TOTAL_PER_SIEVE;
    printf("   sieve q %1.2fM-%1.2fM: %1.0f s at %d threads; whole job ~%1.0f s (%1.1f days)\n",
        q0 / 1e6, qend / 1e6, sieve, T, total, total / 86400);
    fprintf(out, "#   sieve q %1.2fM-%1.2fM: %1.0f s; whole job x%1.2f = %1.0f s%s\n",
        q0 / 1e6, qend / 1e6, sieve, BT_NFS_TOTAL_PER_SIEVE, total,
        own ? " (comparison only, not used)" : "");
    fflush(out);
    free(qs); free(y); free(s); free(secs); free(qv); free(rels);
    return total;
}

// returns the number of reference jobs measured
static int bt_nfs_section(fact_obj_t *fobj, FILE *out, int threads, int compare,
    struct timeval *tstart)
{
    const int nref = (int)(sizeof(bt_nfs_refs) / sizeof(bt_nfs_refs[0]));
    double T[16];
    int i, ok = 0;

    printf("bayestune: NFS test sieving, %d threads: %d reference jobs, about %1.0f "
        "minutes each%s\n", threads, nref, (BT_NFS_BLOCK_SEC + 30.0) / 60.0,
        compare ? " (twice from C150: yafu and NFS@Home parameters)" : "");
    fprintf(out, "# nfs times: bayestune test sieving of reference polys; whole job = "
        "sieve x %1.2f\n", BT_NFS_TOTAL_PER_SIEVE);
    for (i = 0; i < nref; i++)
    {
        T[i] = bt_nfs_ref(fobj, &bt_nfs_refs[i], 0, threads, out, tstart);
        if (T[i] > 0)
            ok++;
        if (compare && (bt_nfs_refs[i].own != NULL))
            bt_nfs_ref(fobj, &bt_nfs_refs[i], 1, threads, out, tstart);
    }
    if (ok == 0)
    {
        printf("bayestune: NFS test sieving failed (sievers in ggnfs_dir?)\n");
        fprintf(out, "# nfs: test sieving failed\n");
        return 0;
    }
    fprintf(out, "# nfs <digits> <wall s, whole job> <threads>\n");
    for (i = 0; i < nref; i++)
        if (T[i] > 0)
            fprintf(out, "nfs %d %.6g %d\n", bt_nfs_refs[i].digits, T[i], threads);
    fflush(out);
    return ok;
}

static void bt_append(char **buf, size_t *len, const char *line)
{
    size_t l = strlen(line);
    *buf = (char *)realloc(*buf, *len + l + 2);
    memcpy(*buf + *len, line, l);
    *len += l;
    if ((*buf)[*len - 1] != '\n')
        (*buf)[(*len)++] = '\n';
    (*buf)[*len] = '\0';
}

void bayes_tune(fact_obj_t *fobj, int what)
{
    static const int sizes[] = { 60, 70, 80, 90, 100, 110, 120, 130, 140,
        150, 160, 170, 180, 190, 200, 220, 250, 280, 310 };
    static const uint64_t b1s[] = { 50000, 250000, 1000000, 3000000 };
    const int nsz = (int)(sizeof(sizes) / sizeof(sizes[0]));
    const int nb1 = (int)(sizeof(b1s) / sizeof(b1s[0]));
    int saved_vflag = fobj->VFLAG;
    int saved_logflag = fobj->LOGFLAG;
    int threads = fobj->THREADS;
    uint32_t ncurves;
    int i, j;
    FILE *out;
    time_t now = time(NULL);
    gmp_randstate_t rs;
    mpz_t n;
    struct timeval tstart, tnow;
    // the previous file, by section: parts not measured this time are kept
    char *old[4] = { NULL, NULL, NULL, NULL };      // header, ecm, siqs, nfs
    size_t oldlen[4] = { 0, 0, 0, 0 };

    if (what <= 0)
        what = BT_DO_ECM | BT_DO_SIQS | BT_DO_NFS;
    if (what & BT_DO_NFS_COMPARE)
        what |= BT_DO_NFS;

    {
        FILE *f = fopen(BAYES_TUNE_FILE, "r");
        char line[1024];
        if (f != NULL)
        {
            while (fgets(line, sizeof(line), f) != NULL)
            {
                int s = 0;
                if ((strncmp(line, "ecm ", 4) == 0) || (strncmp(line, "ecm_eff ", 8) == 0) ||
                    (strncmp(line, "# ecm", 5) == 0))
                    s = 1;
                else if ((strncmp(line, "siqs ", 5) == 0) || (strncmp(line, "# siqs", 6) == 0))
                    s = 2;
                else if ((strncmp(line, "nfs ", 4) == 0) || (strncmp(line, "# nfs", 5) == 0) ||
                    (strncmp(line, "#   ", 4) == 0) || (strncmp(line, "# sieve", 7) == 0))
                    s = 3;
                bt_append(&old[s], &oldlen[s], line);
            }
            fclose(f);
        }
    }

    out = fopen(BAYES_TUNE_FILE, "w");
    if (out == NULL)
    {
        printf("bayestune: could not open %s for writing\n", BAYES_TUNE_FILE);
        for (i = 0; i < 4; i++)
            free(old[i]);
        return;
    }

    if (threads < 1)
        threads = 1;
    // two 8-curve batches per thread (AVX-ECM runs curves in vector batches):
    // one batch leaves per-call setup (threads, stage 2 tables) too visible
    ncurves = 16 * (uint32_t)threads;

    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 31415);
    mpz_init(n);
    fobj->VFLAG = -1;
    fobj->LOGFLAG = 0;
    gettimeofday(&tstart, NULL);

    // header: new when ECM is measured (the ECM lines are per the header's
    // thread count), else the previous one
    if ((what & BT_DO_ECM) || (old[0] == NULL))
    {
        fprintf(out, "# yafu bayestune: timing data for -plan bayesian (docs/ecm_bayes.md)\n");
        fprintf(out, "# re-run bayestune() after changing threads, CPU, or ECM options\n");
        fprintf(out, "version: 1\n");
        fprintf(out, "cpu: %s\n", fobj->CPU_ID_STR);
        fprintf(out, "threads: %d\n", threads);
        fprintf(out, "date: %s", ctime(&now));
        fprintf(out, "ecm_options: avxecm_curve=%d avx_stg2_fft=%d ext_xover=%" PRIu64
            " prefer_gmpecm=%d\n", fobj->ecm_obj.avx_curve, fobj->ecm_obj.avx_stg2_fft,
            fobj->ecm_obj.ecm_ext_xover, fobj->ecm_obj.prefer_gmpecm);
    }
    else
        fputs(old[0], out);

    if (!(what & BT_DO_ECM))
    {
        if (old[1] != NULL)
            fputs(old[1], out);
    }
    else
    {
        fprintf(out, "# ecm <digits> <bits> <B1> <curves> <wall s per curve, all threads> <B2 (0: unknown)>\n");

        printf("bayestune: ECM wall-clock time per curve, %d threads, %u curves per point\n",
            threads, ncurves);
        printf("   C  |        B1 |  s/curve |       B2 | elapsed\n");

        for (i = 0; i < nsz; i++)
        {
            bt_semiprime(n, sizes[i], rs);
            for (j = 0; j < nb1; j++)
            {
                int curves;
                double B2, spc;

                // keep the factors out of reach (C60-C80 have 30-40 digit
                // factors) and the run time down at the largest sizes
                if ((b1s[j] > 1000000) && ((sizes[i] < 90) || (sizes[i] > 200)))
                    continue;

                spc = bt_ecm_point(fobj, n, sizes[i], b1s[j], ncurves, rs, &curves, &B2);
                gettimeofday(&tnow, NULL);
                if (spc <= 0)
                {
                    printf(" %4d | %9" PRIu64 " | factor found 5 times, skipped\n",
                        sizes[i], b1s[j]);
                    continue;
                }
                printf(" %4d | %9" PRIu64 " | %8.4f | %8.3g | %6.0f s\n", sizes[i], b1s[j],
                    spc, B2, ytools_difftime(&tstart, &tnow));
                fflush(stdout);
                fprintf(out, "ecm %d %d %" PRIu64 " %d %.6g %.6g\n", sizes[i],
                    (int)mpz_sizeinbase(n, 2), b1s[j], curves, spc, B2);
                fflush(out);
            }
        }

        // thread efficiency at one point: 1 thread (one batch) vs all threads
        if (threads > 1)
        {
            int curves1, curvesN;
            double B2, s1, sN;
            bt_semiprime(n, 150, rs);
            fobj->THREADS = 1;
            s1 = bt_ecm_point(fobj, n, 150, 250000, 16, rs, &curves1, &B2);
            fobj->THREADS = threads;
            sN = bt_ecm_point(fobj, n, 150, 250000, ncurves, rs, &curvesN, &B2);
            if ((s1 > 0) && (sN > 0))
            {
                printf("bayestune: ECM thread efficiency at C150, B1=250k: %1.3f "
                    "(%1.4f s/curve on 1 thread, %1.4f on %d)\n",
                    s1 / (sN * threads), s1, sN, threads);
                fprintf(out, "ecm_eff 150 250000 %d %.6g %.6g %.4f\n", threads, s1, sN,
                    s1 / (sN * threads));
            }
        }
    }

    // SIQS: whole jobs on random semiprimes, C60..C95.  In memory (no
    // savefile, so a siqs.dat in this directory is left alone).  Larger sizes
    // are extrapolated from the last two lines; NFS takes over near C100.
    if (!(what & BT_DO_SIQS))
    {
        if (old[2] != NULL)
            fputs(old[2], out);
    }
    else
    {
        static const int qsizes[] = { 60, 65, 70, 75, 80, 85, 90, 95 };
        const int nq = (int)(sizeof(qsizes) / sizeof(qsizes[0]));
        uint32_t saved_inmem = fobj->qs_obj.inmem_cutoff;
        double sx = 0, sy = 0, sxx = 0, sxy = 0, a = 0, b = 0;
        double lnt[8];
        int nfit = 0;

        fprintf(out, "# siqs times: bayestune, random semiprimes (raw times in comments)\n");
        printf("bayestune: SIQS whole-job time, %d threads\n", threads);
        fobj->qs_obj.inmem_cutoff = 200;
        for (i = 0; i < nq; i++)
        {
            // single inputs vary by about +-30% (multiplier, factor base):
            // two inputs from C85, where the fit below matters most
            int reps = (qsizes[i] >= 85) ? 2 : 1;
            lnt[i] = 0;
            for (j = 0; j < reps; j++)
            {
                struct timeval s0, s1;
                double tq;

                bt_semiprime(n, qsizes[i], rs);
                reset_factobj(fobj);
                mpz_set(fobj->N, n);
                mpz_set(fobj->qs_obj.gmp_n, n);
                gettimeofday(&s0, NULL);
                SIQS(fobj);
                gettimeofday(&s1, NULL);
                tq = ytools_difftime(&s0, &s1);
                gettimeofday(&tnow, NULL);
                printf("   C%d: %8.2f s | %6.0f s\n", qsizes[i], tq,
                    ytools_difftime(&tstart, &tnow));
                fflush(stdout);
                fprintf(out, "# siqs raw %d %.6g\n", qsizes[i], tq);
                fflush(out);
                lnt[i] += log(tq > 1e-6 ? tq : 1e-6) / reps;
            }
            if (qsizes[i] >= 70)
            {
                sx += qsizes[i];
                sy += lnt[i];
                sxx += (double)qsizes[i] * qsizes[i];
                sxy += qsizes[i] * lnt[i];
                nfit++;
            }
        }
        fobj->qs_obj.inmem_cutoff = saved_inmem;

        // write a least-squares fit of ln t from C70 (below that, fixed
        // overheads bend the curve) so that one slow or fast input doesn't
        // set the extrapolation slope toward the NFS crossover
        if (nfit > 1)
        {
            b = (nfit * sxy - sx * sy) / (nfit * sxx - sx * sx);
            a = (sy - b * sx) / nfit;
        }
        fprintf(out, "# siqs <digits> <wall s, whole job> <threads>: measured below C70, "
            "fit ln t = %.5g + %.5g D (x2 per %.1f digits) from C70\n", a, b,
            (b > 0) ? log(2.0) / b : 0.0);
        for (i = 0; i < nq; i++)
        {
            double v = ((qsizes[i] >= 70) && (nfit > 1)) ? exp(a + b * qsizes[i]) : exp(lnt[i]);
            fprintf(out, "siqs %d %.6g %d\n", qsizes[i], v, threads);
        }
        if ((nfit > 1) && (b > 0))
            printf("bayestune: SIQS fit from C70: x2 per %.1f digits, C100 ~ %1.0f s\n",
                log(2.0) / b, exp(a + b * 100));
    }

    // NFS: test sieving of reference jobs (bt_nfs_refs)
    if (!(what & BT_DO_NFS))
    {
        if (old[3] != NULL)
            fputs(old[3], out);
    }
    else if ((bt_nfs_section(fobj, out, threads, (what & BT_DO_NFS_COMPARE) != 0,
        &tstart) == 0) && (old[3] != NULL))
    {
        // nothing measured: keep the previous nfs lines
        fprintf(out, "# nfs: previous lines kept\n");
        fputs(old[3], out);
    }

    fclose(out);
    fobj->VFLAG = saved_vflag;
    fobj->LOGFLAG = saved_logflag;
    fobj->THREADS = threads;
    reset_factobj(fobj);
    mpz_clear(n);
    gmp_randclear(rs);
    for (i = 0; i < 4; i++)
        free(old[i]);

    gettimeofday(&tnow, NULL);
    printf("bayestune: wrote %s (%1.0f s)\n", BAYES_TUNE_FILE,
        ytools_difftime(&tstart, &tnow));

    bt_ecm.loaded = 0;
    bayes_tune_load(fobj);
}

int bayes_tune_load(fact_obj_t *fobj)
{
    FILE *in;
    char line[1024];

    memset(&bt_ecm, 0, sizeof(bt_ecm));
    in = fopen(BAYES_TUNE_FILE, "r");
    if (in == NULL)
        return 0;

    while (fgets(line, sizeof(line), in) != NULL)
    {
        int d, bits, curves;
        double B1, spc, B2;

        if (strncmp(line, "threads:", 8) == 0)
            bt_ecm.threads = atoi(line + 8);
        else if (strncmp(line, "cpu:", 4) == 0)
        {
            char *s = line + 4;
            while (*s == ' ') s++;
            strncpy(bt_ecm.cpu, s, sizeof(bt_ecm.cpu) - 1);
            bt_ecm.cpu[strcspn(bt_ecm.cpu, "\r\n")] = '\0';
        }
        else if ((strncmp(line, "siqs ", 5) == 0) || (strncmp(line, "nfs ", 4) == 0))
        {
            // siqs|nfs <digits> <wall s> [threads; default: the threads: header]
            int m = (line[0] == 's') ? 0 : 1;
            int thr = 0;
            double dd, sec;
            if ((sscanf(line + (m ? 4 : 5), "%lf %lf %d", &dd, &sec, &thr) >= 2) &&
                (sec > 0) && (bt_ecm.nsieve[m] < BT_MAXSIZES))
            {
                // keep sorted by digits
                int k = bt_ecm.nsieve[m]++;
                while ((k > 0) && (bt_ecm.sdig[m][k - 1] > dd))
                {
                    bt_ecm.sdig[m][k] = bt_ecm.sdig[m][k - 1];
                    bt_ecm.ssec[m][k] = bt_ecm.ssec[m][k - 1];
                    bt_ecm.sthr[m][k] = bt_ecm.sthr[m][k - 1];
                    k--;
                }
                bt_ecm.sdig[m][k] = dd;
                bt_ecm.ssec[m][k] = sec;
                bt_ecm.sthr[m][k] = (thr > 0) ? thr : 0;
            }
        }
        else if ((strncmp(line, "ecm ", 4) == 0) &&
            (sscanf(line + 4, "%d %d %lf %d %lf %lf", &d, &bits, &B1, &curves,
                &spc, &B2) == 6) && (spc > 0))
        {
            int k = bt_ecm.nsizes - 1;
            if ((k < 0) || (bt_ecm.digits[k] != (double)d))
            {
                if (bt_ecm.nsizes == BT_MAXSIZES)
                    continue;
                k = bt_ecm.nsizes++;
                bt_ecm.digits[k] = d;
            }
            if (bt_ecm.nb1[k] < BT_MAXB1)
            {
                bt_ecm.b1[k][bt_ecm.nb1[k]] = B1;
                bt_ecm.spc[k][bt_ecm.nb1[k]] = spc;
                bt_ecm.nb1[k]++;
            }
        }
    }
    fclose(in);

    // sieve lines to thread-seconds (lines without a thread count use the
    // header's)
    {
        int m, k;
        for (m = 0; m < 2; m++)
            for (k = 0; k < bt_ecm.nsieve[m]; k++)
            {
                if (bt_ecm.sthr[m][k] == 0)
                    bt_ecm.sthr[m][k] = (bt_ecm.threads > 0) ? bt_ecm.threads : 1;
                bt_ecm.ssec[m][k] *= (double)bt_ecm.sthr[m][k];
            }
    }

    if ((bt_ecm.nsizes == 0) && (bt_ecm.nsieve[0] < 2) && (bt_ecm.nsieve[1] < 2))
        return 0;
    bt_ecm.loaded = 1;

    if ((fobj != NULL) && (fobj->VFLAG >= 0))
    {
        if (bt_ecm.threads != fobj->THREADS)
            printf("bayestune: %s was measured with %d threads, running with %d: "
                "ECM times scaled linearly; re-run bayestune() for accuracy\n",
                BAYES_TUNE_FILE, bt_ecm.threads, fobj->THREADS);
        if (strcmp(bt_ecm.cpu, fobj->CPU_ID_STR) != 0)
            printf("bayestune: %s was measured on '%s', this is '%s'\n",
                BAYES_TUNE_FILE, bt_ecm.cpu, fobj->CPU_ID_STR);
    }
    return 1;
}

// log-log interpolation in B1 at size index k (extrapolation uses the end
// segment's slope, clamped to [0.8, 1.5])
static double bt_spc_at(int k, double B1)
{
    int m = bt_ecm.nb1[k], j;
    double x = log(B1), x0, x1, y0, y1, s;

    if (m == 1)
        return bt_ecm.spc[k][0] * B1 / bt_ecm.b1[k][0];
    for (j = 1; j < m - 1; j++)
        if (B1 < bt_ecm.b1[k][j])
            break;
    x0 = log(bt_ecm.b1[k][j - 1]);
    x1 = log(bt_ecm.b1[k][j]);
    y0 = log(bt_ecm.spc[k][j - 1]);
    y1 = log(bt_ecm.spc[k][j]);
    s = (y1 - y0) / (x1 - x0);
    if ((B1 < bt_ecm.b1[k][0]) || (B1 > bt_ecm.b1[k][m - 1]))
    {
        if (s < 0.8) s = 0.8;
        if (s > 1.5) s = 1.5;
    }
    return exp(y0 + s * (x - x0));
}

double bayes_ecm_sec_per_curve(fact_obj_t *fobj, double digits, double B1)
{
    int k;
    double t;

    if (!bt_ecm.loaded && !bayes_tune_load(NULL))
        return -1.0;
    if (bt_ecm.nsizes == 0)
        return -1.0;

    if (bt_ecm.nsizes == 1)
    {
        double r = digits / bt_ecm.digits[0];
        t = bt_spc_at(0, B1) * r * r;
    }
    else
    {
        for (k = 1; k < bt_ecm.nsizes - 1; k++)
            if (digits < bt_ecm.digits[k])
                break;
        {
            double d0 = bt_ecm.digits[k - 1], d1 = bt_ecm.digits[k];
            double y0 = log(bt_spc_at(k - 1, B1)), y1 = log(bt_spc_at(k, B1));
            if (digits < d0)
                t = exp(y0) * (digits / d0) * (digits / d0);        // ~ size^2
            else if (digits > d1)
                t = exp(y1) * (digits / d1) * (digits / d1);
            else
                t = exp(y0 + (y1 - y0) * (digits - d0) / (d1 - d0));
        }
    }

    if ((fobj != NULL) && (bt_ecm.threads > 0) && (fobj->THREADS > 0) &&
        (fobj->THREADS != bt_ecm.threads))
        t *= (double)bt_ecm.threads / (double)fobj->THREADS;
    return t;
}

// Whole-job sieve time (wall s, fobj->THREADS threads) for a composite of
// 'digits' digits by SIQS (method 0) or GNFS (method 1), from the "siqs" /
// "nfs" lines of bayes_tune.txt: log-linear interpolation, extrapolation with
// the slope of the nearest two points, scaled linearly from each line's thread
// count.  -1 if there are fewer than two points.
double bayes_sieve_sec(fact_obj_t *fobj, int method, double digits)
{
    int n, k;
    double d0, d1, y0, y1, t;

    if ((method < 0) || (method > 1))
        return -1.0;
    if (!bt_ecm.loaded && !bayes_tune_load(NULL))
        return -1.0;
    n = bt_ecm.nsieve[method];
    if (n < 2)
        return -1.0;
    // no extrapolation more than 5 digits below the first line (the caller
    // then uses the tune() fit): small-size NFS behaves differently (lpb
    // steps at C101/C106, poly and LA overheads), and the lines there are
    // an extrapolation of larger runs anyway
    if (digits < bt_ecm.sdig[method][0] - 5.0)
        return -1.0;
    for (k = 1; k < n - 1; k++)
        if (digits < bt_ecm.sdig[method][k])
            break;
    d0 = bt_ecm.sdig[method][k - 1];
    d1 = bt_ecm.sdig[method][k];
    y0 = log(bt_ecm.ssec[method][k - 1]);
    y1 = log(bt_ecm.ssec[method][k]);
    t = exp(y0 + (y1 - y0) * (digits - d0) / (d1 - d0));

    // thread-seconds -> wall seconds at this run's thread count
    if ((fobj != NULL) && (fobj->THREADS > 0))
        t /= (double)fobj->THREADS;
    else if (bt_ecm.threads > 0)
        t /= (double)bt_ecm.threads;
    return t;
}
