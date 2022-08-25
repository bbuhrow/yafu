#include "arith.h"
#include "factor.h"
#include "qs.h"
#include "gnfs.h"
#include "gmp_xface.h"
#include <stdint.h>

//----------------------- LOCAL DECLARATIONS ----------------------------------//
#define NUM_SIQS_PTS 9
#define NUM_GNFS_PTS 6
#define BASE_e 2.718281828459045

//----------------------- LOCAL FUNCTIONS -------------------------------------//
double best_linear_fit(double *x, double *y, int numpts, 
	double *slope, double *intercept);
void update_INI(double mult, double exponent, double mult2, 
    double exponent2, double xover, double cpu_freq, char *cpu_str);
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
	char siqslist[9][200];
	char nfslist[6][200];
    char sievername[1024];
	int i, tmpT;
    mpz_t n;
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	double t_time;

	//uint32_t siqs_actualrels[NUM_SIQS_PTS] = {17136, 32337,
		//63709, 143984, 240663, 568071, 793434, 1205061, 1595268};
	uint32_t siqs_actualrels[NUM_SIQS_PTS] = {17136, 32337,
		63709, 143984, 242825, 589192, 847299, 1272852, 1709598}; // 2397232, 3293925
    uint32_t siqs_tf_small_cutoff[NUM_SIQS_PTS] = { 15, 20,
        15, 13, 16, 20, 18, 19, 20 }; // 20, 20

	double siqs_extraptime[NUM_SIQS_PTS];
	double siqs_sizes[NUM_SIQS_PTS] = {60, 65, 70, 75, 80, 85, 90, 95, 100};

	uint32_t gnfs_actualrels[NUM_GNFS_PTS] = {1929527, 4074867, 4783410, 
		4969315, 5522597, 5783845};
	double gnfs_extraptime[NUM_GNFS_PTS];
	double gnfs_sizes[NUM_GNFS_PTS] = {85, 90, 95, 100, 105, 110};
	double gnfs_max_poly_time[NUM_GNFS_PTS] = {0.09, 0.12, 0.21, 0.34, 0.5, 1.0};	//HRS

	double a, b, a2, b2, fit, xover;
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

	//siqs: start with c60, increment by 5 digits, up to a c100
	//this will allow determination of NFS/QS crossover as well as provide enough
	//info to generate an equation for QS time estimation
	strcpy(siqslist[0],"349594255864176572614071853194924838158088864370890996447417");
	strcpy(siqslist[1],"34053408309992030649212497354061832056920539397279047809781589871");
	strcpy(siqslist[2],"6470287906463336878241474855987746904297564226439499503918586590778209");
	strcpy(siqslist[3],"281396163585532137380297959872159569353696836686080935550459706878100362721");
	strcpy(siqslist[4],"33727095233132290409342213138708322681737322487170896778164145844669592994743377");
	strcpy(siqslist[5],"1877138824359859508015524119652506869600959721781289179190693027302028679377371001561");
	strcpy(siqslist[6],"427351849650748515507228344120452096326780093349980867041485502247153375067354165128307841");
	strcpy(siqslist[7],"48404068520546498995797968938385187958997290617596242601254422967869040251141325866025672337021");
	strcpy(siqslist[8],"1802716097522165018257858828415111497060066282677325501816640492782221110851604465066510547671104729");

	//nfs: start with c85, increment by 5 digits, up to C110
	//this will allow determination of NFS/QS crossover
	//to do NFS time estimation, probably need to go much higher - say c155.
	strcpy(nfslist[0],"1877138824359859508015524119652506869600959721781289179190693027302028679377371001561");
	strcpy(nfslist[1],"427351849650748515507228344120452096326780093349980867041485502247153375067354165128307841");
	strcpy(nfslist[2],"48404068520546498995797968938385187958997290617596242601254422967869040251141325866025672337021");
	strcpy(nfslist[3],"1802716097522165018257858828415111497060066282677325501816640492782221110851604465066510547671104729");
	strcpy(nfslist[4],"466734409955806375058988820327650664396976790744285564594552020197119774272189758795312820988691316775181");
	strcpy(nfslist[5],"48178889479314834847826896738914354061668125063983964035428538278448985505047157633738779051249185304620494013");

	// c135 from TMPQS
	// 116915434112329921568236283928181979297762987646390347857868153872054154807376462439621333455331738807075404918922573575454310187518221

	// RSA-120
	// 227010481295437363334259960947493668895875336466084780038173258247009162675779735389791151574049166747880487470296548479

	// RSA-129 (The magic words are squeamish ossifrage)
	// 114381625757888867669235779976146612010218296721242362562561842935706935245733897830597123563958705058989075147599290026879543541

	// for each of the siqs inputs
	for (i=0; i<NUM_SIQS_PTS; i++)
	{
		fact_obj_t* fobj = inobj; // (fact_obj_t*)malloc(sizeof(fact_obj_t));
		//init_factobj(fobj);

		// measure how long it takes to gather a fixed number of relations 		
        mpz_set_str(n, siqslist[i], 10);
		fobj->qs_obj.gbl_override_rel_flag = 1;
		fobj->qs_obj.gbl_override_rel = 10000;	

        // also set the tf_small_cutoff to its known best value so we
        // don't pollute the measurement with the optimization process.
        fobj->qs_obj.gbl_override_small_cutoff_flag = 1;
        fobj->qs_obj.gbl_override_small_cutoff = siqs_tf_small_cutoff[i];

		gettimeofday(&start, NULL);
        mpz_set(fobj->qs_obj.gmp_n, n);
		SIQS(fobj);
		gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);

		// the number of relations actually gathered is stored in gbl_override_rel
		siqs_extraptime[i] = t_time * siqs_actualrels[i] / fobj->qs_obj.gbl_override_rel;

		// add a guess at the linalg + sqrt time: something reasonable seems to be about
		// 2% of the total sieve time
		siqs_extraptime[i] += 0.02 * siqs_extraptime[i];

		printf("elapsed time for ~10k relations of c%d = %6.4f seconds.\n",
            mpz_sizeinbase(n, 10), t_time);
		printf("extrapolated time for complete factorization = %6.4f seconds\n",siqs_extraptime[i]);

		clear_factor_list(fobj->factors);

		//free_factobj(fobj);
		//free(fobj);
	}

	fit = best_linear_fit(siqs_sizes, siqs_extraptime, NUM_SIQS_PTS, &a, &b);
	printf("best linear fit is ln(y) = %g * x + %g\nR^2 = %g\n",a,b,fit);
	printf("best exponential fit is y = %g * exp(%g * x)\n",pow(BASE_e,b),a);

	// for each of the gnfs inputs
	for (i=0; i<NUM_GNFS_PTS; i++)
	{
		char syscmd[1024];
		FILE *in;
		uint32_t startq, qrange;		
		double t_time2, d;

		//remove previous tests
		remove("tunerels.out");
		remove("tune.job.afb.0");
		MySleep(.1);

		//make a job file for each input
		make_job_file(sievername, &startq, &qrange, nfslist[i], i, inobj);

		//measure how long it takes to generate the afb... for fun.		
		gettimeofday(&start, NULL);

		//create the afb - we don't want the time it takes to do this to
		//pollute the sieve timings
		sprintf(syscmd,"%s -b tune.job -k -c 0 -F", sievername);

		printf("nfs: commencing construction of afb\n");

		system(syscmd);
		gettimeofday(&stop, NULL);
        t_time2 = ytools_difftime(&start, &stop);
		
		printf("afb generation took %6.4f seconds.\n",t_time2);
		remove("tune.job.afb.0");
		MySleep(.1);

		//measure how long it takes to sieve a fixed range of special-q	
		gettimeofday(&start, NULL);

		//start the test
		sprintf(syscmd,"%s -a tune.job -f %u -c %u -o tunerels.out",
			sievername, startq, qrange);
		printf("nfs: commencing lattice sieving over range: %u - %u\n",
			startq, startq + qrange);
		system(syscmd);
		gettimeofday(&stop, NULL);
        t_time = ytools_difftime(&start, &stop);

		//count relations
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

		//extrapolate
		gnfs_extraptime[i] = (double)gnfs_actualrels[i] * (t_time - t_time2) / (double)count;
		printf("nfs: extrapolated time for sieving = %6.4f seconds\n",gnfs_extraptime[i]);

		//add estimated linalg and sqrt time.  
		//reasonable estimate is 10% of the total runtime?
		d = gnfs_extraptime[i] * 0.1;
		gnfs_extraptime[i] += d;
		printf("adding %g seconds for linalg and sqrt\n",d);	

		//add avg poly time 
		// estimated at 5% of sieving time
		d = gnfs_extraptime[i] * 0.05; // gnfs_max_poly_time[i] * 3600 / 2;
		gnfs_extraptime[i] += d;
		printf("adding %g seconds for average polyselect time\n",d);

		printf("nfs: estimated time for complete factorization = %6.4f seconds\n",gnfs_extraptime[i]);
		
	}

	// set THREADS back the way it was
    inobj->THREADS = tmpT;

	fit = best_linear_fit(gnfs_sizes, gnfs_extraptime, NUM_GNFS_PTS, &a2, &b2);
	printf("best linear fit is ln(y) = %g * x + %g\nR^2 = %g\n",a2,b2,fit);
	printf("best exponential fit is y = %g * exp(%g * x)\n",pow(BASE_e,b2),a2);

	xover = (b2 - b) / (a - a2);
	printf("QS/NFS crossover occurs at %2.1f digits\n",xover);

	//write the coefficients to the .ini file
	update_INI(pow(BASE_e,b),a,pow(BASE_e,b2),a2,xover,
        inobj->MEAS_CPU_FREQUENCY, inobj->CPU_ID_STR);

	mpz_clear(n);
	return;
}

/*
//digits, r/alim, lpbr/a, mfbr/a, r/alambda, siever, rels
//entries based on statistics gathered from many factorizations done
//over the years by myself and others, and from here:
//http://www.mersenneforum.org/showthread.php?t=12365
#define GGNFS_TABLE_ROWS 15
static double ggnfs_table[GGNFS_TABLE_ROWS][8] = {
	{85,  900000,   24, 48, 2.5, 11, 0, 10000},
	{90,  1200000,  25, 50, 2.5, 11, 0, 20000},
	{95,  1500000,  25, 50, 2.5, 12, 0, 20000},
	{100, 1800000,  26, 52, 2.5, 12, 0, 40000},
	{105, 2500000,  26, 52, 2.5, 12, 0, 40000},
	{110, 3200000,  26, 52, 2.5, 13, 0, 40000},
	{115, 4500000,  27, 54, 2.5, 13, 0, 80000},
	{120, 5000000,  27, 54, 2.5, 13, 0, 80000},
	{125, 5500000,  27, 54, 2.5, 13, 0, 80000},
	{130, 6000000,  27, 54, 2.5, 13, 0, 80000},
	{135, 8000000,  27, 54, 2.5, 14, 0, 80000},
	{140, 12000000, 28, 56, 2.5, 14, 0, 160000},
	{145, 15000000, 28, 56, 2.5, 14, 0, 160000},
	{150, 20000000, 29, 58, 2.5, 14, 0, 320000},
	{155, 30000000, 29, 58, 2.5, 15, 0, 320000}
};
*/

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

	//fast test
	//*qrange = 1000;

	//more accurate test
	// *qrange = 5000;

	*qrange = 2000;

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
    double exponent2, double xover, double cpu_freq, char* cpu_str)
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
			len = strlen(str);
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
			fputs(str, out);
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
				sprintf(newline, "tune_info=%s,WIN64,%lg,%lg,%lg,%lg,%lg,%lg\n",
                    cpu_str,mult,exponent,mult2,exponent2,xover, cpu_freq);
				fputs(newline, out);
			}
#elif defined(WIN32)
			if ((strcmp(cpustr, cpu_str) == 0) && (strcmp(osstr, "WIN32") == 0))
			{
				printf("Replacing tune_info entry for %s - %s\n",osstr,cpustr);
				found_entry = 1;
				sprintf(newline, "tune_info=%s,WIN32,%lg,%lg,%lg,%lg,%lg,%lg\n",
                    cpu_str,mult,exponent,mult2,exponent2,xover, cpu_freq);
				fputs(newline, out);
			}
#else 
			if ((strcmp(cpustr, cpu_str) == 0) && (strcmp(osstr, "LINUX64") == 0))
			{
				printf("Replacing tune_info entry for %s - %s\n",osstr,cpustr);
				found_entry = 1;
				sprintf(newline, "tune_info=%s,LINUX64,%lg,%lg,%lg,%lg,%lg,%lg\n",
                    cpu_str,mult,exponent,mult2,exponent2,xover, cpu_freq);
				fputs(newline, out);
			}
#endif
			else
			{
				//just write the line
				fputs(str2, out);
			}
		}
		else
		{
			//just write the line
			fputs(str2, out);
		}
	}

	if (!found_entry)
	{
#if defined(_WIN64)
		printf("Adding tune_info entry for WIN64 - %s\n", cpu_str);
		sprintf(newline, "tune_info=%s,WIN64,%lg,%lg,%lg,%lg,%lg,%lg\n",
            cpu_str,mult,exponent,mult2,exponent2,xover, cpu_freq);
		fputs(newline, out);
#elif defined(WIN32)
		printf("Adding tune_info entry for WIN32 - %s\n", cpu_str);
		sprintf(newline, "tune_info=%s,WIN32,%lg,%lg,%lg,%lg,%lg,%lg\n",
            cpu_str,mult,exponent,mult2,exponent2,xover, cpu_freq);
		fputs(newline, out);
#else 
		printf("Adding tune_info entry for LINUX64 - %s\n", cpu_str);
		sprintf(newline, "tune_info=%s,LINUX64,%lg,%lg,%lg,%lg,%lg,%lg\n",
            cpu_str,mult,exponent,mult2,exponent2,xover, cpu_freq);
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

double best_linear_fit(double *x, double *y, int numpts, 
	double *slope, double *intercept)
{
	// given vectors of x and y values, compute the best linear fit
	// to the data and return the slope and y-intercept of the line
	// return the R^2 value.
	// follows: http://mathworld.wolfram.com/LeastSquaresFitting.html
	double avgX, avgY, varX, varY, cov, *yy;
	int i;

	yy = (double *)malloc(numpts * sizeof(double));

	// linearize the y-axis
	for (i=0; i<numpts; i++)
		yy[i] = log(y[i]);

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


