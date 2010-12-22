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
#include "monty.h"
#include "factor.h"
#include "util.h"

/*
implements pollard's p-1 algorihm.  inspired by Bressoud's book, and
various papers by montgomery and kruppa
*/

// these are used by the top level function, so both YAFU and GMP-ECM
// paths must use these prototypes
void pm1_init();
void pm1_finalize();
void pm1exit(int sig);
int mpollard(z *n, uint32 c, z *f);
void pm1_print_B1_B2(z *n, FILE *flog);

// declarations/definitions when using YAFU's ECM
#if !defined(HAVE_GMP) || !defined(HAVE_GMP_ECM)

	typedef struct
	{
		z G2;
		z G1;
		z G0;
		int first_use;
	} brent_suyama_t;

	brent_suyama_t bs_table;

	void table_init(z *b, uint32 a, uint32 D, z *n);
	void next_fa(z *fa, z *n);
	

	void pm1_init()
	{
		if (POLLARD_STG2_MAX > 0xEE6B2800)
		{
			fprintf(stderr,"primes greater than 4e9 not supported in P-1, reducing.\n");
			POLLARD_STG2_MAX = 0xEE6B2800;
		}

		//initialize global brent-suyama table
		zInit(&bs_table.G0);
		zInit(&bs_table.G1);
		zInit(&bs_table.G2);

		return;
	}

	void pm1_finalize()
	{
		zFree(&bs_table.G0);
		zFree(&bs_table.G1);
		zFree(&bs_table.G2);

		return;
	}

	int mpollard(z *n, uint32 c, z *f)
	{
		//run pollard's p-1 algorithm on n, returning any factor found in f, or else return
		//0 in f.
		//use base c

		z m,d,mm,t1,cc,s;
		z *t;
		int D, a, b;
		int brent_suyama=1;
		int pairing=1;
		int i,j,blocknum;
		int q,charcount,paired=0;
		int exitcode = 0;
		int resume1 = 0, resume2 = 0, save = 0;
		int pm1_blocksz = 1024;

		clock_t start, stop;
		double tt1=0;
		clock_t start2, stop2;
		double tt2=0;
		clock_t update_start, update_stop;
		double t_update,update_time=5.0;

		uint8 marks[2310];
		uint8 nmarks[2310];
		z *powarray;
		FILE *instate;
		char *str;

		if (POLLARD_STG1_MAX <= 2310)
		{
			if (POLLARD_STG1_MAX <= 210)
			{
				printf("POLLARD_STG1_MAX too low\n");
				return 1;
			}
			D = 210;
		}
		else
			D = 2310;

		zInit(&m);
		zInit(&s);
		zInit(&d);
		zInit(&t1);
		zInit(&mm);
		zInit(&cc);
		
		//build an array to hold values of f(b)
		t = (z *)malloc(D * sizeof(z));
		for (j=0;j<D;j++)
			zInit(&t[j]);

		//build an array to hold products of powers
		powarray = (z *)malloc(8192*sizeof(z));
		for (j=0;j<8192;j++)
			zInit(&powarray[j]);

		start = clock();

		//
		str = (char *)malloc(1024*sizeof(char));

		//check to see if a stage 2 resume file exists for this input
		instate = fopen("pm1_stg2_save.dat","r");
		if (instate != NULL)
		{
			if (!feof(instate))
			{
				//read in n.  this might not read in everything... needs
				//to check if there are more than 1024 characters in n
				fgets(str,1024,instate);
				str2hexz(str,&m);
				if (zCompare(&m,n) == 0)
				{
					printf("restarting stage 2 from saved data set\n");
					resume2 = 1;
				}
			}
			fclose(instate);
		}

		//check to see if a stage 1 resume file exists for this input, 
		//only if there is not a stage 2
		instate = fopen("pm1_stg1_save.dat","r");
		if (instate != NULL && !resume2)
		{
			if (!feof(instate))
			{
				//read in n.  this might not read in everything... needs
				//to check if there are more than 1024 characters in n
				fgets(str,1024,instate);
				str2hexz(str,&m);
				if (zCompare(&m,n) == 0)
				{
					printf("restarting stage 1 from saved data set\n");
					resume1 = 1;
				}
			}
			fclose(instate);
		}

		if (resume1)
		{
			recover_stg1(n,&m,POLLARD_METHOD);
			monty_init(n);
		}
		else if (resume2)
		{
			charcount=0;
			i=0;
			recover_stg2(n,&cc,&m,&mm,POLLARD_METHOD);
			monty_init(n);
			goto resume_stg2;
		}
		else
		{
			monty_init(n);
			m.val[0] = c;
			to_monty(&m,n);

			//make sure the first primes up to 10M are computed
			GetPRIMESRange(0,10001000);
		}

		if (resume1)
			printf("resuming at prime %u\n",PRIMES[0]);

		i=0;
		charcount=0;
		blocknum=0;
		update_start = clock();
		while (PRIMES[i] < POLLARD_STG1_MAX)
		{
			/*
			//not simple, doesn't work yet
			//for each block of 8192 primes, fill the array and multiply
			limit = 8192*(blocknum+1);
			for (j=0;(i<limit) && (PRIMES[i]<POLLARD_STG1_MAX) && (j<8192);i++,j++)
			{
				if ((uint32)i >= NUM_P)
				{
					GetPRIMESRange(PRIMES[i-1]-1,PRIMES[i-1]+10001000);
					i=0;
				}
				q = PRIMES[i];
				powarray[j].val[0] = (uint32)pow(q,floor(log(POLLARD_STG1_MAX)/log(q)));
				powarray[j].size = 1;
			}
			
			max_index = j-1;
			limit = (max_index)/2;

			while (1)
			{
				//multiply neighbors
				for (j=0;j<limit;j++)
				{
					//if there is no second number to multiply, copy it down instead
					//handles non power of 2 block sizes.
					if (2*j+1 > max_index)
						zCopy(&powarray[2*j],&powarray[j]);
					else
					{
						zMul(&powarray[2*j],&powarray[2*j+1],&t1);
						zCopy(&t1,&powarray[j]);
					}
				}
				if (powarray[0].size > MAX_DIGITS/2 - 1 || limit == 1)
					break;

				max_index = limit-1;
				limit >>= 1;
			}

			//exponentiate with each of these
			for (j=0;j<limit;j++)
			{
				//sliding window based binary power ladder on a large exponent
				zmModExpw(&m,&powarray[j],&t1,9);
				zCopy(&t1,&m);
			}

			*/

			//slightly less simple - 6% improvement or so.
			//for some block of primes, multiply all powers of those primes which
			//are less than the stg1 bound together into one big number.  then
			//use a sliding window exponentiation method, which is noticably faster
			//when there are a lot of bits to work with and the window can 
			//be made large.
			zCopy(&zOne,&mm);
			for (j=0;(PRIMES[i]<POLLARD_STG1_MAX) && (j<pm1_blocksz);i++,j++)
			{
				//compute more primes if we need to.  it is more efficient to 
				//get a large chunk, even if we won't use them all.
				if ((uint32)i >= NUM_P)
				{
					GetPRIMESRange(PRIMES[i-1]-1,PRIMES[i-1]+10001000);
					i=0;
				}

				q = PRIMES[i];
				while ((POLLARD_STG1_MAX/PRIMES[i]) > (uint32)q) q *= PRIMES[i];

				zShortMul(&mm,(fp_digit)q,&mm);
			}

			//sliding window based binary power ladder on a large exponent.
			//the window size has been loosely hand optimized for the prime power
			//block size of 1024.
			zmModExpw(&m,&mm,&t1,n,9);
			zCopy(&t1,&m);

			//check for a factor and update the screen (under verbose) every 
			//'update_time' seconds.
			update_stop = clock();
			t_update = (double)(update_stop - update_start)/(double)CLOCKS_PER_SEC;
			if (t_update > update_time)
			{
				//verbose: update the screen with the last prime processed.
				if (VFLAG >= 0)
				{
					for (j=0;j<charcount;j++)
						printf("\b");
					charcount = printf(", processed < %u",PRIMES[i]);
					fflush(stdout);
				}

				update_start = clock();
				t_update = 0;

				//check for an abort
				if (PM1_ABORT)
				{
					save_state(STG1,n,&m,NULL,NULL,0,POLLARD_METHOD);
					printf("state saved\n");
					goto free;
				}

				//check gcd for a factor
				monty_sub(&m,&montyconst.one,&d,n);

				zLEGCD(&d,n,f);
				if (zCompare(f,&zOne) > 0)
				{
					if (zCompare(f,n) == 0)
						zCopy(&zZero,f);
					exitcode=1;
					//record stage 1 time
					stop = clock();
					tt1 = (double)(stop-start)/(double)CLOCKS_PER_SEC;
					goto free;
				}
			}
		
			//next block
			blocknum++;
		}

		//all done with stage 1.  we may not have hit the 5 sec update point
		//recently, so do a final factor check here.
		if (VFLAG >= 0)
		{
			for (j=0;j<charcount;j++)
				printf("\b");
			charcount = printf(", processed < %u",PRIMES[i]);
			fflush(stdout);
		}

		monty_sub(&m,&montyconst.one,&d,n);

		zLEGCD(&d,n,f);
		if (zCompare(f,&zOne) > 0)
		{
			if (zCompare(f,n) == 0)
				zCopy(&zZero,f);
			exitcode=1;
			stop = clock();
			tt1 = (double)(stop-start)/(double)CLOCKS_PER_SEC;
			goto free;
		}

		stop = clock();
		tt1 = (double)(stop-start)/(double)CLOCKS_PER_SEC;

		//optionally, save the stage1 state.
		if (save)
			save_state(STG1,n,&m,NULL,NULL,0,POLLARD_METHOD);

		//start here if resuming from stage2 savefile
	resume_stg2:
		if (resume2)
			printf("resuming at prime %u",PRIMES[0]);
		
		start2 = clock();
		//initialize for stage 2 using m, the residue from stage 1.
		//i indicates the first prime bigger than POLLARD_STG1_MAX
		if (brent_suyama)
		{
			//we use a function of b, f(b), rather than b itself.
			//the function is f(b) = b^s, where s is even.  for now pick s=2 
			//other choices include other even s, and dickman's rho function, parameter s.
			d.size=1;
			for (j=1;j<D;j++)
			{
				if (spGCD(j,D) == 1)
				{
					d.val[0] = j;
					zExp(2,&d,&t1);
					zmModExp(&m,&t1,&t[j],n);	//find and store all x^f(b)
				}
			}
		}
		else
		{
			d.size=1;
			for (j=1;j<D;j++)
			{
				if (spGCD(j,D) == 1)
				{
					d.val[0] = j;
					zmModExp(&m,&d,&t[j],n);	//find and store all x^b
				}
			}
		}

		if (brent_suyama)
		{
			//we use a function of a, f(a), rather than a itself.
			//the function is f(a) = a^s, where s is even.  for now pick s=2
			//other choices include other even s, and dickman's rho function, parameter s.
			a = PRIMES[i]/D + (PRIMES[i]%D != 0);
			a *= D;
			d.size=1;
			d.val[0] = a;
			zExp(2,&d,&t1);
			zmModExp(&m,&t1,&mm,n);	//initialize to x^f(a)
			//init table to for rapidly stepping through polynomial evaluations
			//of our brent suyama function
			table_init(&m,a,D,n);
		}
		else
		{
			a = PRIMES[i]/D + (PRIMES[i]%D != 0);
			a *= D;
			d.val[0] = a;
			zmModExp(&m,&d,&mm,n);	//initialize to x^a
			d.val[0] = D;
			zmModExp(&m,&d,&s,n);	//s is the giant step = x^D
		}

		zCopy(&montyconst.one,&cc);
		
		memset(&marks,0,2310*sizeof(uint8));
		memset(&nmarks,0,2310*sizeof(uint8));
		paired=0;
		update_start = clock();
		t_update=0;
		//stage 2
		while (PRIMES[i]<POLLARD_STG2_MAX)
		{
			b = a - PRIMES[i];
			if (pairing)
			{
				if (!marks[b])
				{
					//not marked, so doesn't have a match on the other side of the previous a
					//accumulate it, and mark the next range of a
					monty_sub(&mm,&t[b],&d,n);
					monty_mul(&d,&cc,&t1,n);
					zCopy(&t1,&cc);
					nmarks[D-b]=1;
				}
				else
				{
					//it's marked, so don't need to accumulate it - it has a pair that is already
					//accumulated.
					paired++;
				}
			}
			else
			{
				monty_sub(&mm,&t[b],&d,n);
				monty_mul(&d,&cc,&t1,n);
				zCopy(&t1,&cc);
			}
			i++;
			
			if ((uint32)i >= NUM_P)
			{
				GetPRIMESRange(PRIMES[i-1]+1,PRIMES[i-1]+10001000);
				i=0;
			}

			if (PRIMES[i] > (uint32)a)
			{
				//set marks = nextmarks, then clear nextmarks
				if (brent_suyama)
					next_fa(&mm,n);
				else
				{
					//giant step
					monty_mul(&mm,&s,&t1,n);
					zCopy(&t1,&mm);
				}

				if (pairing)
				{
					memcpy(&marks,&nmarks,2310*sizeof(uint8));
					memset(&nmarks,0,2310*sizeof(uint8));
				}

				//next range
				a += D;
			}

			update_stop = clock();
			t_update = (double)(update_stop - update_start)/(double)CLOCKS_PER_SEC;
			if (t_update > update_time)
			{
				if (PM1_ABORT)
				{
					save_state(STG2,n,&m,&mm,&cc,i,POLLARD_METHOD);
					printf("state saved\n");
					goto free;
				}

				update_start = clock();
				t_update = 0;

				if (VFLAG >= 0)
				{
					for (j=0;j<charcount;j++)
						printf("\b");
					charcount = printf(", processed < %u",PRIMES[i]);
					fflush(stdout);
				}

				zLEGCD(&cc,n,f);
				if (zCompare(f,&zOne) > 0)
				{
					if (zCompare(f,n) == 0)
					{
						zCopy(&zZero,f);
					}
					exitcode=1;
					stop2 = clock();
					tt2 = (double)(stop2-start2)/(double)CLOCKS_PER_SEC;
					goto free;
				}
			}
		}
		if (VFLAG >= 0)
		{
			for (j=0;j<charcount;j++)
				printf("\b");
			charcount = printf(", processed < %u",PRIMES[i]);
			fflush(stdout);
		}

		zLEGCD(&cc,n,f);
		if (zCompare(f,&zOne) > 0)
		{
			if (zCompare(f,n) == 0)
			{
				zCopy(&zZero,f);
			}
			exitcode=1;
			stop2 = clock();
			tt2 = (double)(stop2-start2)/(double)CLOCKS_PER_SEC;
			goto free;
		}

		stop2 = clock();
		tt2 = (double)(stop2-start2)/(double)CLOCKS_PER_SEC;

		if (save)
			save_state(STG2,n,&m,&mm,&cc,i,POLLARD_METHOD);

	free:

		if (VFLAG >= 0)
			printf("\n");

		for (j=0;j<8192;j++)
			zFree(&powarray[j]);
		free(powarray);
		zFree(&m);
		zFree(&d);
		zFree(&t1);
		zFree(&mm);
		zFree(&cc);
		zFree(&s);
		for (j=0;j<D;j++)
			zFree(&t[j]);
		free(t);
		free(str);
		monty_free();

		if (exitcode == 0)
			zCopy(&zZero,f);

		return i;
	}

	void table_init(z *b, uint32 a, uint32 D, z *n)
	{
		/*
		build a table to facilitate computing sucessive values of f(a+jh)
		where a and h are fixed

		per montgomery, define
		Gd(n) = n^d
		Gi(n) = Gi+1(n+h) - Gi+1(n), i = d-1, d-2, ... 0

		for example for G(n) = n^2, n0 = 10, h = 3 we have
		G2(n) = 100
		G1(n) = G2(n+h) - G2(n) = 13^2 - 10^2 = 69
		G0(n) = G1(n+h) - G1(n) = G2(n+2h) - 2G2(n+h) + G2(n) = 16^2 - 2*13^2 + 10^2 = 18

		G0 is always constant

		note that we get the next value(s) in the sequence for free during the build of
		the table, in the example, they are G2(n+h) = 169 and G2(n+2h) = 256

		sucessive values can be computed as follows
		for now, stick with the case n^2.  some day make more general...
		G(n+kh) = G2 + 2*G1 + G0
		G2 = G2 + G1
		G1 = G1 + G0
		this is accomplished in the function next_fa  
		*/

		z t,t2,t3;
		zInit(&t);
		zInit(&t2);
		zInit(&t3);

		//compute b^G2(a)
		t.size = 1;
		t.val[0] = a;
		zExp(2,&t,&t2);
		zmModExp(b,&t2,&bs_table.G2,n);

		//compute b^G1(a)
		t.size = 1;
		t.val[0] = a+D;
		zExp(2,&t,&t3);		
		zSub(&t3,&t2,&t);		//(a+D)^2 - a^2	
		zmModExp(b,&t,&bs_table.G1,n);

		//compute b^G0(a)
		zAdd(&t3,&t,&t2);		//2*(a+D)^2 - a^2	
		t.size = 1;
		t.val[0] = a+2*D;
		zExp(2,&t,&t3);
		zSub(&t3,&t2,&t);		//(a+2D)^2 - (2*(a+D)^2 - a^2)
		zmModExp(b,&t,&bs_table.G0,n);

		bs_table.first_use = 1;

		zFree(&t);
		zFree(&t2);
		zFree(&t3);
		return;
	}

	void next_fa(z *fa, z *n)
	{
		//given a and h, return f(a+h), using the prebuilt table
		//update the table

		//there are n+1 table entries, where n is the order of the brent-suyama function

		//for n=2 compute f(a+h) = G2*(G1^2)*G0 mod N

		z t,t2;
		zInit(&t);
		zInit(&t2);

		if (bs_table.first_use)
		{
			//first call to next fa, simplified result
			monty_mul(&bs_table.G2,&bs_table.G1,&t,n);		//G2*G1 = next G2

			zCopy(&t,fa);									//G2*G1 = next fa
			bs_table.first_use = 0;
		}
		else
		{
			monty_mul(&bs_table.G2,&bs_table.G1,&t,n);		//G2*G1 = next G2
			zCopy(&t,&bs_table.G2);

			monty_mul(&bs_table.G1,&bs_table.G0,&t2,n);		//G1*G0 = next G1
			zCopy(&t2,&bs_table.G1);

			monty_mul(&t,&t2,fa,n);							//G2*G1*G1*G0 = next fa
		}

		zFree(&t);
		zFree(&t2);
		return;
	}	

#else

	// declarations and definitions used when GMP-ECM and GMP are 
	// available
	#include <gmp_xface.h>
	#include <ecm.h>

	uint64 TMP_STG2_MAX;

	typedef struct
	{
		mpz_t gmp_n, gmp_factor;
		ecm_params params;
		uint32 sigma;
		int stagefound;

	} ecm_pm1_data_t;

	ecm_pm1_data_t pm1_data;

	void pm1_init()
	{
		mpz_init(pm1_data.gmp_n);
		mpz_init(pm1_data.gmp_factor);
		ecm_init(pm1_data.params);
		//gmp_randseed_ui(tdata->params->rng, 
		//	get_rand(&obj->seed1, &obj->seed2));

		pm1_data.params->method = ECM_PM1;
		//pm1_data.params->verbose = 1;

		TMP_STG2_MAX = POLLARD_STG2_MAX;

		return;
	}

	void pm1_finalize()
	{
		ecm_clear(pm1_data.params);
		mpz_clear(pm1_data.gmp_n);
		mpz_clear(pm1_data.gmp_factor);

		POLLARD_STG2_MAX = TMP_STG2_MAX;

		//needs an extra return to make the screen output look right
		if (VFLAG >= 0)
			printf("\n");
		
		return;
	}

	int mpollard(z *n, uint32 c, z *f)
	{
		int status;
		size_t count;

		pm1_data.params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
		//pm1_data.params->verbose = 2;

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		mpz_import(pm1_data.gmp_n, (size_t)(abs(n->size)), -1, sizeof(uint32), 
			0, (size_t)0, n->val);
#else
		//wrapper for YAFU bigints and call to gmp-ecm
		mp2gmp(n, pm1_data.gmp_n);
#endif

		if (PM1_STG2_ISDEFAULT == 0)
		{
			//not default, tell gmp-ecm to use the requested B2
			//printf("using requested B2 value\n");
			sp642z(POLLARD_STG2_MAX,f);
			mp2gmp(f,pm1_data.params->B2);
			zClear(f);
		}

		status = ecm_factor(pm1_data.gmp_factor, pm1_data.gmp_n,
				POLLARD_STG1_MAX, pm1_data.params);

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


#endif

// top level routine: the only one visible to the rest of the program
void pollard_loop(fact_obj_t *fobj)
{
	//run pollard's p-1 algorithm once on the input, using a 
	//32 bit random base
	z *n = &fobj->pm1_obj.n;
	z d,f,t;   
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
	flog = fopen(fobj->logname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->logname);
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
	zInit(&f);
	zInit(&t);

	pm1_init();
		
	base = spRand(3,0xFFFFFFFF);

	pm1_print_B1_B2(n,flog);
	mpollard(n,base,&f);
		
	if (zCompare(&f,&zOne) > 0 && zCompare(&f,n) < 0)
	{
		//non-trivial factor found
		stop = clock();
		tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;
		if (isPrime(&f))
		{
			f.type = PRP;
			add_to_factor_list(fobj, &f);
			//log result
			logprint(flog,"prp%d = %s\n",
				ndigits(&f),z2decstr(&f,&gstr2));
		}
		else
		{
			f.type = COMPOSITE;
			add_to_factor_list(fobj, &f);
			//log result
			logprint(flog,"c%d = %s\n",
				ndigits(&f),z2decstr(&f,&gstr2));
		}
		start = clock();

		//reduce input
		zDiv(n,&f,&t,&d);
		zCopy(&t,n);
	}

	fclose(flog);

	pm1_finalize();

	//watch for an abort
	if (PM1_ABORT)
	{
		print_factors(fobj);
		exit(1);
	}

	signal(SIGINT,NULL);
	zFree(&d);
	zFree(&f);
	zFree(&t);
	return;
}

void pm1_print_B1_B2(z *n, FILE *flog)
{
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (POLLARD_STG1_MAX % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",POLLARD_STG1_MAX / 1000000000, suffix);
	}
	else if (POLLARD_STG1_MAX % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",POLLARD_STG1_MAX / 1000000, suffix);
	}
	else if (POLLARD_STG1_MAX % 1000 == 0)
	{
		suffix = 'K';
		sprintf(stg1str,"%u%c",POLLARD_STG1_MAX / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",POLLARD_STG1_MAX);
	}

	if (POLLARD_STG2_MAX % 1000000000 == 0)
	{
		suffix = 'B';
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu%c",POLLARD_STG2_MAX / 1000000000, suffix);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu%c",POLLARD_STG2_MAX / 1000000000, suffix);
#else
		sprintf(stg2str,"%I64u%c",POLLARD_STG2_MAX / 1000000000, suffix);
#endif
	}
	else if (POLLARD_STG2_MAX % 1000000 == 0)
	{
		suffix = 'M';
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu%c",POLLARD_STG2_MAX / 1000000, suffix);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu%c",POLLARD_STG2_MAX / 1000000, suffix);
#else
		sprintf(stg2str,"%I64u%c",POLLARD_STG2_MAX / 1000000, suffix);
#endif
	}
	else if (POLLARD_STG2_MAX % 1000 == 0)
	{
		suffix = 'K';
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu%c",POLLARD_STG2_MAX / 1000, suffix);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu%c",POLLARD_STG2_MAX / 1000, suffix);
#else
		sprintf(stg2str,"%I64u%c",POLLARD_STG2_MAX / 1000, suffix);
#endif
	}
	else
	{
#if defined(__unix__) && (BITS_PER_DIGIT == 64)
		sprintf(stg2str,"%lu",POLLARD_STG2_MAX);
#elif defined(__unix__) && (BITS_PER_DIGIT == 32)
		sprintf(stg2str,"%llu",POLLARD_STG2_MAX);
#else
		sprintf(stg2str,"%I64u",POLLARD_STG2_MAX);
#endif
	}

	if (VFLAG >= 0)
	{
		printf("pm1: starting B1 = %s, B2 = %s on C%d",
			stg1str,stg2str,ndigits(n));
		fflush(stdout);	
	}
	logprint(flog,"pm1: starting B1 = %s, B2 = %s on C%d\n",
		stg1str,stg2str,ndigits(n));
	
	return;
}

void pm1exit(int sig)
{
	printf("\nAborting...\n");
	PM1_ABORT = 1;
	return;
}

