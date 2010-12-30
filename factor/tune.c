#include "yafu_string.h"
#include "arith.h"
#include "factor.h"
#include "qs.h"
#include "gnfs.h"

#define NUM_SIQS_PTS 5
double best_linear_fit(double *x, double *y, int numpts, 
	double *slope, double *intercept);

void factor_tune(void)
{
	// perform trial runs on a set of inputs with known solutions using siqs and nfs
	// compute crossover points and time estimation coefficients
	// write all this info to the .ini file
	char siqslist[9][200];
	char nfslist[6][200];
	z n;
	int i;
	uint32 save_gbl_override_rel;
	int save_gbl_override_rel_flag;
	struct timeval stop;	// stop time of this job
	struct timeval start;	// start time of this job
	TIME_DIFF *	difference;
	double t_time;
	uint32 siqs_actualrels[NUM_SIQS_PTS] = {17136, 32337,
		63709, 143984, 240663}; //, 568071, 793434, 1205061};
	double siqs_extraptime[NUM_SIQS_PTS];
	double siqs_sizes[NUM_SIQS_PTS] = {60, 65, 70, 75, 80}; //, 85, 90, 95, 100};
	double a, b, fit;

	zInit(&n);

	//siqs: start with c60, increment by 5 digits, up to a c100
	//this will allow determination of NFS/QS crossover as well as provide enough
	//info to generate a linear equation for QS time estimation
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
	//c85: 1877138824359859508015524119652506869600959721781289179190693027302028679377371001561
	//requires 1929527 relations
	//c90: 427351849650748515507228344120452096326780093349980867041485502247153375067354165128307841
	//requires 4074867 relations
	//c95: 48404068520546498995797968938385187958997290617596242601254422967869040251141325866025672337021
	//requires 4783410 relations
	//c100: 1802716097522165018257858828415111497060066282677325501816640492782221110851604465066510547671104729
	//requires 4969315 relations
	//c105: 466734409955806375058988820327650664396976790744285564594552020197119774272189758795312820988691316775181
	//requires 5522597 relations
	//c110: 48178889479314834847826896738914354061668125063983964035428538278448985505047157633738779051249185304620494013
	//requires 5783845 relations

	save_gbl_override_rel = gbl_override_rel;
	save_gbl_override_rel_flag = gbl_override_rel_flag;
		
	gbl_override_rel_flag = 1;

	// for each of the siqs inputs
	for (i=0; i<NUM_SIQS_PTS; i++)
	{
		fact_obj_t *fobj = (fact_obj_t *)malloc(sizeof(fact_obj_t));
		init_factobj(fobj);
		zCopy(&fobj->fobj_factors[i].factor,&fobj->N);

		//measure how long it takes to gather a fixed number of relations 		
		str2hexz(siqslist[i],&n);
		gbl_override_rel = 10000;	
		gettimeofday(&start, NULL);
		zCopy(&n,&fobj->qs_obj.n);
		SIQS(fobj);
		gettimeofday(&stop, NULL);
		difference = my_difftime (&start, &stop);
		t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
		free(difference);			
		// the number of relations actually gathered is stored in gbl_override_rel
		siqs_extraptime[i] = t_time * siqs_actualrels[i] / gbl_override_rel;
		printf("elapsed time for ~10k relations of c%d = %6.4f seconds.\n",ndigits(&n),t_time);
		printf("extrapolated time for complete factorization = %6.4f seconds\n",siqs_extraptime[i]);

		clear_factor_list(fobj);

		free_factobj(fobj);
		free(fobj);
	}

	gbl_override_rel = save_gbl_override_rel;		
	gbl_override_rel_flag = save_gbl_override_rel_flag;

	fit = best_linear_fit(siqs_sizes, siqs_extraptime, NUM_SIQS_PTS, &a, &b);
	printf("best linear fit is y = %g * x + %g\nR^2 = %g\n",pow(2.718281828459045,a),
		pow(2.718281828459045,b),fit);

	zFree(&n);
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

	printf("average x = %g, average y = %g, varX = %g, varY = %g, cov = %g\n",
		avgX, avgY, varX, varY, cov);
	
	*slope = cov / varX;
	*intercept = avgY - *slope * avgX;

	free(yy);
	return cov * cov / (varX * varY);
}


