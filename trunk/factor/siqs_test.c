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
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"
#include "factor.h"
#include "tfm.h"
#include "util.h"

int check_relation(z *a, z *b, siqs_r *r, fb_list *fb, z *n)
{
	int offset, lp[2], parity, num_factors;
	int j,retval;
	z Q, RHS,t1,t2;

	zInit(&Q);
	zInit(&RHS);
	zInit(&t1);
	zInit(&t2);

	offset = r->sieve_offset;
	lp[0] = r->large_prime[0];
	lp[1] = r->large_prime[1];
	parity = r->parity;
	num_factors = r->num_factors;

	RHS.val[0] = lp[0];
	zShortMul(&RHS,lp[1],&RHS);
	for (j=0; j<num_factors; j++)
	{
		zShortMul(&RHS,fb->list->prime[r->fb_offsets[j]],&t1);
		zCopy(&t1,&RHS);
	}

	//Q(x)/a = (ax + b)^2 - N, where x is the sieve index
	zShortMul(a,offset,&t1);
	if (parity)
		zSub(&t1,b,&t2);
	else
		zAdd(&t1,b,&t2);
	zMul(&t2,&t2,&t1);
	zSub(&t1,n,&Q);

	retval = 0;
	Q.size = abs(Q.size);
	if (zCompare(&Q,&RHS) != 0)
	{
		//printf("failure to equate relation\n");
		//printf("%s %s\n",z2decstr(&Q,&gstr1),z2decstr(&RHS,&gstr2));
		retval = 1;
	}

	zFree(&Q);
	zFree(&RHS);
	zFree(&t1);
	zFree(&t2);
	return retval;
}

int check_specialcase(z *n, FILE *sieve_log, fact_obj_t *fobj)
{
	z w1,w2;

	zInit(&w1);
	zInit(&w2);
	//check for some special cases of input number
	//sieve_log is passed in already open, and should return open
	if ((n->val[0] % 2) == 0)
	{
		printf("input must be odd\n");
		zFree(&w1);
		zFree(&w2);
		return 1;
	}

	if (isPrime(n))
	{
		n->type = PRP;
		add_to_factor_list(fobj, n);
		if (sieve_log != NULL)
			logprint(sieve_log,"prp%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
		zCopy(&zOne,n);
		zFree(&w1);
		zFree(&w2);
		return 1;
	}

	if (zBits(n) < 115)
	{
		//run MPQS, as SIQS doesn't work for smaller inputs
		//MPQS will take over the log file, so close it now.
		int i;

		if (sieve_log != NULL)
			fclose(sieve_log);

		// we've verified that the input is not odd or prime.  also
		// do some very quick trial division before calling smallmpqs, which
		// does none of these things.
		for (i=1; i<25; i++)
		{
			if (zShortMod(n, spSOEprimes[i]) == 0)
				zShortDiv(n, spSOEprimes[i], n);
		}

		zCopy(n,&fobj->qs_obj.n);
		smallmpqs(fobj);
		for (i=0; i<fobj->qs_obj.num_factors; i++)
		{
			if (isPrime(&fobj->qs_obj.factors[i]))
			{
				fobj->qs_obj.factors[i].type = PRP;
				add_to_factor_list(fobj, &fobj->qs_obj.factors[i]);
			}
			else
			{
				fobj->qs_obj.factors[i].type = COMPOSITE;
				add_to_factor_list(fobj, &fobj->qs_obj.factors[i]);
			}
			zCopy(&fobj->qs_obj.n,&w1);
			zDiv(&w1,&fobj->qs_obj.factors[i],&fobj->qs_obj.n,&w2);
			zFree(&fobj->qs_obj.factors[i]);
		}
		fobj->qs_obj.num_factors = 0;
		zCopy(&fobj->qs_obj.n,n);
		zFree(&w1);
		zFree(&w2);
		return 2;	//tells SIQS to not try to close the logfile
	}

	if (ndigits(n) > 150)
	{
		printf("input too big for SIQS\n");
		zFree(&w1);
		zFree(&w2);
		return 1;
	}

	zNroot(n,&w1,2);
	zMul(&w1,&w1,&w2);
	if (zCompare(&w2,n) == 0)
	{
		printf("input is a perfect square\n");
		if (isPrime(&w1))
		{
			if (sieve_log != NULL)
			{
				logprint(sieve_log,"prp%d = %s\n",ndigits(&w1),z2decstr(&w1,&gstr1));
				logprint(sieve_log,"prp%d = %s\n",ndigits(&w1),z2decstr(&w1,&gstr1));
			}
			w1.type = PRP;
			add_to_factor_list(fobj, &w1);
			add_to_factor_list(fobj, &w1);
		}
		else
		{
			if (sieve_log != NULL)
			{
				logprint(sieve_log,"c%d = %s\n",ndigits(&w1),z2decstr(&w1,&gstr1));
				logprint(sieve_log,"c%d = %s\n",ndigits(&w1),z2decstr(&w1,&gstr1));
			}
			w1.type = COMPOSITE;
			add_to_factor_list(fobj, &w1);
			add_to_factor_list(fobj, &w1);
		}
		zCopy(&zOne,n);
		zFree(&w1);
		zFree(&w2);
		return 1;
	}

	zFree(&w1);
	zFree(&w2);
	return 0;
}

int checkpoly_siqs(siqs_poly *poly, z *n)
{
	//check that b^2 == N mod a
	//and that c == (b*b - n)/a
	z t1,t2,t3,t4;

	zInit(&t1);
	zInit(&t2);
	zInit(&t3);
	zInit(&t4);

	zCopy(n,&t1);
	zDiv(&t1,&poly->poly_a,&t2,&t3);

	zMul(&poly->poly_b,&poly->poly_b,&t2);
	zDiv(&t2,&poly->poly_a,&t1,&t4);

	if (zCompare(&t3,&t4) != 0)
	{
		printf("\nError in checkpoly: %s^2 !== N mod %s\n",z2decstr(&poly->poly_b,&gstr1),
			z2decstr(&poly->poly_a,&gstr2));
		if (poly->poly_b.size < 0)
			printf("b is negative\n");
	}

	if (zJacobi(n,&poly->poly_a) != 1)
		printf("\nError in checkpoly: (a|N) != 1\n");

	zSqr(&poly->poly_b,&t1);
	zSub(&t1,n,&t2);
	zDiv(&t2,&poly->poly_a,&t1,&t3);

	if (zCompare(&t1,&poly->poly_c) != 0)
		printf("\nError in checkpoly: c != (b^2 - n)/a\n");
	
	zFree(&t1);
	zFree(&t2);
	zFree(&t3);
	zFree(&t4);
	return 0;
}

int checkBl(z *n, uint32 *qli, fb_list *fb, z *Bl, int s)
{
	//check that Bl^2 == N mod ql and Bl == 0 mod qj for 1<=j<=s, j != l
	int i,j,p,q;
	z t1,t2;

	zInit(&t1);
	zInit(&t2);

	for (i=0;i<s;i++)
	{
		p = fb->list->prime[qli[i]];

		zShiftRight(&t2,&Bl[i],1);
		zSqr(&t2,&t1);
		if (zShortMod(&t1,p)%zShortMod(n,p) != 0)
			printf("%s^2 mod %u != N mod %u\n",z2decstr(&Bl[i],&gstr1),p,p);

		for (j=0;j<s;j++)
		{
			if (j==i)
				continue;

			q = fb->list->prime[qli[j]];
			zShiftRight(&t2,&Bl[i],1);
			if (zShortMod(&t2,q) != 0)
				printf("%s mod %u != 0\n",z2decstr(&Bl[i],&gstr1),q);
		}
	}

	zFree(&t1);
	zFree(&t2);
	return 0;
}

void test_polya(fb_list *fb, z *n, z *target_a)
{
	uint32 i, cumulativebits=0;
	z tmp, tmp2, tmp3, poly_a;
	clock_t start, stop;
	double t;
	//int *qli, s;

	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);
	zInit(&poly_a);

	start = clock();
	printf("SIQS poly gen:\n");
	printf("target_a bits = %d\n",zBits(target_a));
	for (i=0;i<100;i++)
	{
		//new_poly_a(&poly_a,target_a,n,fb,qli,&s);

		//print target_a, poly_a, and their difference
		//printf("target_a = ");
		//zHex2Dec(target_a,&tmp);
		//zPrintDec(&tmp);
		//printf("\n");
		//printf("poly_a(%d) = ",zBits(&poly_a));
		//zHex2Dec(&poly_a,&tmp);
		//zPrintDec(&tmp);
		//printf("\n");
		
		zSub(target_a,&poly_a,&tmp);
		printf("error = %d bits = ",zBits(&tmp));
		cumulativebits += zBits(&tmp);
		//zHex2Dec(&tmp,&tmp2);
		//zPrintDec(&tmp2);
		//printf("\n");
	}
	printf("average bit difference: %4.2f\n",(double)cumulativebits/100.0);
	stop = clock();
	t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
	printf("Total elapsed time = %6.4f seconds.\n",t);

	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&poly_a);
	return;
}

void siqsbench(fact_obj_t *fobj)
{
	//run through the list of benchmark siqs factorizations

	char list[10][200];
	char oldflogname[80];
	enum cpu_type cpu;
	z n;
	FILE *log;
	int i;

	zInit(&n);
	strcpy(oldflogname,flogname);

	strcpy(list[0],"405461849292216354219321922871108605045931309");
	strcpy(list[1],"29660734457033883936073030405220515257819037444591");
	strcpy(list[2],"1592739380299790264815370870751381488173344970710983383");
	strcpy(list[3],"349594255864176572614071853194924838158088864370890996447417");
	strcpy(list[4],"34053408309992030649212497354061832056920539397279047809781589871");
	strcpy(list[5],"6470287906463336878241474855987746904297564226439499503918586590778209");
	strcpy(list[6],"281396163585532137380297959872159569353696836686080935550459706878100362721");
	strcpy(list[7],"33727095233132290409342213138708322681737322487170896778164145844669592994743377");
	strcpy(list[8],"7456482836301983072751757080079609980368702375378513429852397523678294751191007081");
	strcpy(list[9],"1877138824359859508015524119652506869600959721781289179190693027302028679377371001561");

	strcpy(flogname,"bench.log");
	log = fopen(flogname,"a");
	if (log == NULL)
	{
		printf("couldn't open %s for writing\n",flogname);
		exit(1);
	}

	cpu = yafu_get_cpu_type();
	fprintf(log,"detected cpu %d, with L1 = %d bytes, L2 = %d bytes\n",cpu,L1CACHE,L2CACHE);
#if defined(TFM_X86) || defined(TFM_X86_MSVC)
	fprintf(log,"Initialized with Tom's Fast Math (x86-32 asm)\n\n");
#elif defined(TFM_X86_64)
	fprintf(log,"Initialized with Tom's Fast Math (x86-64 asm)\n\n");
#else
	fprintf(log,"Initialized as (x86-32 generic)...\n\n");
#endif

	fclose(log);

	for (i=0; i<10; i++)
	{
		str2hexz(list[i],&n);
		zCopy(&n,&fobj->qs_obj.n);
		SIQS(fobj);
		clear_factor_list(fobj);
	}

	strcpy(flogname,oldflogname);

	zFree(&n);
	return;
}

void get_dummy_params(int bits, uint32 *B, uint32 *M, uint32 *NB)
{
	int i;
	double scale;

	//parameter table
	//bits, fb primes, lp mulitplier, 64k blocks
	int param_table[22][4] = {
		{140,	600,	40,	1},
		{149,	875,	40,	1},
		{165,	1228,	50,	1},
		{181,	2247,	50,	1},
		{198,	3485,	60,	2},
		{215,	6357,	60,	2},	
		{232,	12132,	70,	3},
		{248,	26379,	80,	4},
		{265,	42871,	90,	6},
		{281,	55137,	100,	8},
		{298,	65244,	120,	10},
		{310,	78247,	120,	12},
		{320,	90678,	140,	14},
		{330,	105000, 150,    18},
		{340,	125000, 150,    21},
		{350,	155000, 150,    25},
		{360,	195000, 150,    29},
		{370,	250000, 150,    34},
		{380,	310000, 150,    40},
		{390,	380000, 150,    47},
		{400,	460000, 150,    55},
		{410,	550000, 150,    64},
	};

	//linear interpolation according to bit size to determine
	//factor base bound.  use the closest parameter for lp multiplier
	//and number of blocks.

	if (bits <= param_table[0][0])
	{
		scale = (double)bits / (double)param_table[0][0];
		*B = (uint32)(scale * (double)(param_table[0][1]));		
		*M = 40;
		*NB = 1;
	}
	else
	{
		for (i=0;i<22;i++)
		{
			if (bits > param_table[i][0] && bits <= param_table[i+1][0])
			{
				scale = (double)(param_table[i+1][0] - bits) /
					(double)(param_table[i+1][0] - param_table[i][0]);
				*B = param_table[i+1][1] - 
					(uint32)(scale * (double)(param_table[i+1][1] - param_table[i][1]));
				
				*M = (uint32)((param_table[i+1][2] + param_table[i][2])/2.0 + 0.5);
				*NB = (uint32)((param_table[i+1][3] + param_table[i][3])/2.0 + 0.5);
			}
		}
	}

	if (*B == 0)
	{
		//off the end of the table, extrapolate based on the slope of 
		//the last two

		scale = (double)(param_table[21][1] - param_table[20][1]) /
			(double)(param_table[21][0] - param_table[20][0]);
		*B = (uint32)(((double)bits - param_table[21][0]) * 
			scale + param_table[21][1]);
		*M = param_table[21][2];	//reuse last one

		scale = (double)(param_table[21][3] - param_table[20][3]) /
			(double)(param_table[21][0] - param_table[20][0]);
		*NB = (uint32)(((double)bits - param_table[21][0]) * 
			scale + param_table[21][3]);
	}

	return;
}

void siqstune(int bits)
{
	FILE *res;
	char fstr[80];
	int i,x,y,w,n;
	z input;
	uint32 B, M, NB;
	//dimension 1 = number of primes
	//dimension 2 = number of blocks
	//dimension 3 = large prime multiplier
	uint32 Bvec[5], NBvec[3];
	fact_obj_t *fobj;
	double rels_per_sec, qs_time, tot_time;
	double rels_per_sec_norm[45], qs_time_norm[45], tot_time_norm[45];

	zInit(&input);

	for (n=150; n <= 275; n+=25)
	{
		for (i=0; i<5; i++)
		{
			//pick a random N bits number
			build_RSA(n, &input);

			//get the default parameters for this input
			get_dummy_params(n, &B, &M, &NB);

			//adjust for various architectures
			if (BLOCKSIZE == 32768)
				NB *= 2;
			if (BLOCKSIZE == 16384)
				NB *= 4;

			//and adjust for really small jobs
			if (bits <= 140)
				NB = 1;
			
			//build up the test vector.  
			//vary # of blocks +/- 20%
			//vary # of primes +/- 10% and +/- 5%
			//vary lp mult +/- 20%
			//Bvec[0] = (uint32)((double)B * 0.9);
			Bvec[0] = (uint32)((double)B * 0.95); 
			Bvec[1] = (uint32)((double)B * 1);
			Bvec[2] = (uint32)((double)B * 1.05); 
			//Bvec[4] = (uint32)((double)B * 1.1);

			//Mvec[0] = (uint32)(floor((double)M * 0.8));
			//Mvec[1] = (uint32)((double)M * 1); 
			//Mvec[2] = (uint32)(ceil((double)M * 1.2));

			NBvec[0] = (uint32)MAX(floor(((double)NB * 0.8)),1);
			NBvec[1] = (uint32)MAX(((double)NB * 1),1); 
			NBvec[2] = (uint32)MAX(ceil(((double)NB * 1.2)),1);

			sprintf(fstr,"tune_results_%d.csv",n);
			res = fopen(fstr,"a");
			fprintf(res,"N: %s\n",z2decstr(&input,&gstr1));
			fclose(res);

			//loop over the test vectors
			for (x=0; x<3; x++)
			{
				for (y=0; y<3; y++)
				{
					//for (w=0; w<3; w++)
					//{

						//new factorization object
						fobj = (fact_obj_t *)malloc(sizeof(fact_obj_t));
						init_factobj(fobj);

						//force these parameters into SIQS
						gbl_override_B_flag = 1;
						gbl_override_B = Bvec[x];

						gbl_override_blocks_flag = 1;
						gbl_override_blocks = NBvec[y];

						//gbl_override_lpmult_flag = 1;
						//gbl_override_lpmult = Mvec[w];

						//notification of progress
						//printf("\n\n================ RUN %d of %d =================\n\n",
						//	i*45 + x*9 + y*3 + w + 1,3*5*3*3);
						printf("\n\n================ RUN %d of %d =================\n\n",
							i*9 + x*3 + y + 1,5*3*3);

						//run siqs
						zCopy(&input,&fobj->qs_obj.n);
						SIQS(fobj);

						//get timing from fact_obj
						rels_per_sec = fobj->qs_obj.rels_per_sec;
						qs_time = fobj->qs_obj.qs_time;
						tot_time = fobj->qs_obj.total_time;
						rels_per_sec_norm[i*9 + x*3 + y] = fobj->qs_obj.rels_per_sec;
						qs_time_norm[i*9 + x*3 + y] = fobj->qs_obj.qs_time;
						tot_time_norm[i*9 + x*3 + y] = fobj->qs_obj.total_time;					

						//clean up
						system("rm siqs.dat");
						free_factobj(fobj);
						free(fobj);

					//}
				}
			}

			sprintf(fstr,"tune_results_%d.csv",n);
			res = fopen(fstr,"a");
			fprintf(res,"# primes,# blocks,rels/sec,qs time,tot time,norm rels/sec,"
				"norm qs time,norm tot time\n");
			for (y=0; y<3; y++)
			{
				for (w=0; w<3; w++)
				{
					/*
					fprintf(res,"%d,%d,%d,%f,%f,%f,%f,%f,%f\n",
						Bvec[x],
						NBvec[y],
						Mvec[w],
						rels_per_sec_norm[x*9 + y*3 + w],
						qs_time_norm[x*9 + y*3 + w],
						tot_time_norm[x*9 + y*3 + w],
						rels_per_sec_norm[x*9 + y*3 + w] / rels_per_sec_norm[2*9 + 1*3 + 1],
						qs_time_norm[x*9 + y*3 + w] / qs_time_norm[2*9 + 1*3 + 1],
						tot_time_norm[x*9 + y*3 + w] / tot_time_norm[2*9 + 1*3 + 1]);
						*/
					fprintf(res,"%d,%d,%f,%f,%f,%f,%f,%f\n",
						Bvec[y],
						NBvec[w],
						rels_per_sec_norm[i*9 + y*3 + w],
						qs_time_norm[i*9 + y*3 + w],
						tot_time_norm[i*9 + y*3 + w],
						rels_per_sec_norm[i*9 + y*3 + w] / rels_per_sec_norm[i*9 + 1*3 + 1],
						qs_time_norm[i*9 + y*3 + w] / qs_time_norm[i*9 + 1*3 + 1],
						tot_time_norm[i*9 + y*3 + w] / tot_time_norm[i*9 + 1*3 + 1]);
				}
			}
			fclose(res);
		}

		//write to results file
		sprintf(fstr,"tune_results_%d.csv",n);
		res = fopen(fstr,"a");
		fprintf(res,"AVERAGES\n");
		fprintf(res,"# primes,# blocks,avg rels/sec norm,avg qs time norm,avg tot time norm\n");
		for (x=0; x<3; x++)
		{
			for (y=0; y<3; y++)
			{
				rels_per_sec = 0;
				qs_time = 0;
				tot_time = 0;
				for (w=0; w<5; w++)
				{
					rels_per_sec += rels_per_sec_norm[w*9 + x*3 + y] / rels_per_sec_norm[w*9 + 1*3 + 1];
					qs_time += qs_time_norm[w*9 + x*3 + y] / qs_time_norm[w*9 + 1*3 + 1];
					tot_time += tot_time_norm[w*9 + x*3 + y] / tot_time_norm[w*9 + 1*3 + 1];
				}
				fprintf(res,"%d,%d,%f,%f,%f\n",
					Bvec[x],
					NBvec[y],
					rels_per_sec / 5,
					qs_time / 5,
					tot_time / 5);
			}
		}
		fclose(res);
	}

	zFree(&input);
	return;
}

