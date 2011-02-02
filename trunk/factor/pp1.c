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
implements williams's p+1 algorihm.  inspired by Bressoud's book, and
various papers by montgomery and kruppa
*/

// these are used by the top level function, so both YAFU and GMP-ECM
// paths must use these prototypes
void pp1_init();
void pp1_finalize();
void pp1_print_B1_B2(z *n, FILE *flog);
int mwilliams(z *n, z *f, uint32 rp);
void pp1exit(int sig);

// declarations/definitions when using YAFU's ECM
#if !defined(HAVE_GMP) || !defined(HAVE_GMP_ECM)

	void mNextV(z *h, uint32 j,z *mTwo, z *n);

	void pp1_init()
	{
		if (WILL_STG2_MAX > 0xEE6B2800)
		{
			fprintf(stderr,"primes greater than 4e9 not supported in P+1, reducing.\n");
			WILL_STG2_MAX = 0xEE6B2800;
		}

		PP1_ABORT = 0;
		signal(SIGINT,pp1exit);

		return;
	}

	void pp1_finalize()
	{
		signal(SIGINT,NULL);

		return;
	}

	int mwilliams(z *n, z *f, uint32 rp)
	{
		//run Williams' p+1 method on n returning any factor found or 0 in f
		
		z v, vv, t2,s,v_ad,mTwo, acc,t1;
		z *t;
		int exitcode=0;
		int i = 0, j,charcount;
		int c = 1;
		uint32 q;

		clock_t start, stop;
		double tt1=0;
		clock_t start2, stop2;
		double tt2=0;
		clock_t update_start=0, update_stop;
		double t_update,update_time=5.0;

		int a,b,D;
		int brent_suyama=1;
		int pairing=1,paired=0;
		int resume1 = 0, resume2 = 0, save = 0;
		uint8 marks[2310];
		uint8 nmarks[2310];
		FILE *instate;
		char *str;

		if (WILL_STG1_MAX <= 2310)
		{
			if (WILL_STG1_MAX <= 210)
			{
				printf("WILL_STG1_MAX too low\n");
				return 1;
			}
			D = 210;
		}
		else
			D = 2310;

		zInit(&mTwo);
		zInit(&v);
		zInit(&vv);
		zInit(&s);
		zInit(&v_ad);
		zInit(&acc);
		zInit(&t1);
		zInit(&t2);
		zClear(f);

		//build an array to hold values of f(b)
		t = (z *)malloc(D * sizeof(z));
		for (j=0;j<D;j++)
			zInit(&t[j]);

		start = clock();

		str = (char *)malloc(1024*sizeof(char));

		//check to see if a stage 2 resume file exists for this input
		instate = fopen("pp1_stg2_save.dat","r");
		if (instate != NULL)
		{
			if (!feof(instate))
			{
				//read in n.  this might not read in everything... needs
				//to check if there are more than 1024 characters in n
				fgets(str,1024,instate);
				str2hexz(str,&v);
				if (zCompare(&v,n) == 0)
				{
					printf("restarting stage 2 from saved data set\n");
					resume2 = 1;
				}
			}
			fclose(instate);
		}

		//check to see if a stage 1 resume file exists for this input, 
		//only if there is not a stage 2
		instate = fopen("pp1_stg1_save.dat","r");
		if (instate != NULL && !resume2)
		{
			if (!feof(instate))
			{
				//read in n.  this might not read in everything... needs
				//to check if there are more than 1024 characters in n
				fgets(str,1024,instate);
				str2hexz(str,&v);
				if (zCompare(&v,n) == 0)
				{
					printf("restarting stage 1 from saved data set\n");
					resume1 = 1;
				}
			}
			fclose(instate);
		}

		if (resume1)
		{
			recover_stg1(n,&v,WILLIAMS_METHOD);
			monty_init(n);
			zCopy(&zTwo,&mTwo);
			to_monty(&mTwo,n);
		}
		else if (resume2)
		{
			charcount=0;
			i=0;
			recover_stg2(n,&acc,&v,&vv,WILLIAMS_METHOD);
			monty_init(n);
			zCopy(&zTwo,&mTwo);
			to_monty(&mTwo,n);
			goto resume_stg2;
		}
		else
		{
			monty_init(n);
			v.val[0] = rp;
			to_monty(&v,n);
			zCopy(&zTwo,&mTwo);
			to_monty(&mTwo,n);

			//make sure the first primes up to 10M are computed
			GetPRIMESRange(0,10001000);
		}

		if (resume1)
			printf("resuming at prime %u\n",PRIMES[0]);

		i=0;
		charcount = 0;
		update_start = clock();
		t_update=0;
		while (PRIMES[i] < WILL_STG1_MAX)
		{
			//compute more primes if we need to.  it is more efficient to 
			//get a large chunk, even if we won't use them all.
			if ((uint32)i >= NUM_P)
			{
				GetPRIMESRange(PRIMES[i-1]-1,PRIMES[i-1]+10001000);
				i=0;
			}

			//check for a factor and update the screen (under verbose) every 
			//'update_time' seconds.
			update_stop = clock();
			t_update = (double)(update_stop - update_start)/(double)CLOCKS_PER_SEC;
			if (t_update > update_time)
			{
				if (PP1_ABORT)
				{
					save_state(STG1,n,&v,NULL,NULL,0,WILLIAMS_METHOD);
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

				monty_sub(&v,&mTwo,&t2,n);
				zLEGCD(&t2,n,f);
				if (zCompare(f,&zOne) > 0)
				{
					if (zCompare(f,n) == 0)
						zCopy(&zZero,f);
					exitcode=1;
					stop = clock();
					tt1 = (double)(stop-start)/(double)CLOCKS_PER_SEC;
					goto free;
				}
			}

			q = PRIMES[i];
			c = (int)pow(q,floor(log(WILL_STG1_MAX)/log(q)));

			mNextV(&v,c,&mTwo,n);	
			i++;
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
		stop = clock();
		tt1 = (double)(stop-start)/(double)CLOCKS_PER_SEC;

		monty_sub(&v,&mTwo,&t2,n);
		zLEGCD(&t2,n,f);
		if (zCompare(f,&zOne) > 0)
		{
			if (zCompare(f,n) == 0)
				zCopy(&zZero,f);
			exitcode=1;
			stop = clock();
			tt1 = (double)(stop-start)/(double)CLOCKS_PER_SEC;
			goto free;
		}

		if (save)
			save_state(STG1,n,&v,NULL,NULL,0,WILLIAMS_METHOD);

	resume_stg2:

		if (resume2)
			printf("resuming at prime %u",PRIMES[0]);

		start2 = clock();

		zCopy(&v,&t[1]);
		mNextV(&t[1],1,&mTwo,n);
		a=1;
		
		for (j=2;j<D;j++)
		{
			if (spGCD(j,D) == 1)
			{
				//we should be able to build on previously computed values,
				//but for some reason this breaks stage2 - factors that 
				//should be found are not.
				//zCopy(&t[a],&t[j]);
				//mNextV(&t[j],j-a,&mTwo);
				//a=j;
				
				zCopy(&v,&t[j]);
				mNextV(&t[j],j,&mTwo,n);		//find and store all V_b(m)	
			}
		}

		a = PRIMES[i]/D + (PRIMES[i]%D != 0);
		a *= D;

		//compute V_a
		zCopy(&v,&vv);
		mNextV(&vv,a,&mTwo,n);		

		//initialize info needed for giant step
		zCopy(&v,&s);
		mNextV(&s,D,&mTwo,n);		
		zCopy(&v,&v_ad);
		mNextV(&v_ad,a-D,&mTwo,n);
			
		zCopy(&montyconst.one,&acc);
		
		memset(&marks,0,2310*sizeof(uint8));
		memset(&nmarks,0,2310*sizeof(uint8));
		paired=0;

		//stage 2
		while (PRIMES[i]<WILL_STG2_MAX)
		{
			b = a - PRIMES[i];
			if (brent_suyama && pairing)
			{
				if (!marks[b])
				{
					//not marked, so doesn't have a match on the other side of the previous a
					//accumulate it, and mark the next range of a
					//we accumulate Va - Vb
					monty_sub(&vv,&t[b],&t2,n);
					monty_mul(&t2,&acc,&t1,n);
					zCopy(&t1,&acc);
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
				//accumulate Va - Vb
				monty_sub(&vv,&t[b],&t2,n);
				monty_mul(&t2,&acc,&t1,n);
				zCopy(&t1,&acc);
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
				if (brent_suyama && pairing)
				{
					memcpy(&marks,&nmarks,2310*sizeof(uint8));
					memset(&nmarks,0,2310*sizeof(uint8));
					
					//p+1 appears to not need brent_suyama in order to pair...
					//giant step - use the addition formula for V
					//s is the giant step (Vd)
					zCopy(&vv,&t2);			//t2 is temp
					monty_mul(&vv,&s,&t1,n);
					monty_sub(&t1,&v_ad,&vv,n);		//V_a+d = VaVd - Va-d
					zCopy(&t2,&v_ad);			//v_ad holds the previous vv
				}
				else
				{
					//giant step - use the addition formula for V
					zCopy(&vv,&t2);			//t2 is temp
					monty_mul(&vv,&s,&t1,n);
					monty_sub(&t1,&v_ad,&vv,n);		//V_a+d = VaVd - Va-d
					zCopy(&t2,&v_ad);			//v_ad holds the previous vv
				}
				//next range
				a += D;
			}

			update_stop = clock();
			t_update = (double)(update_stop - update_start)/(double)CLOCKS_PER_SEC;
			if (t_update > update_time)
			{
				if (PP1_ABORT)
				{
					save_state(STG2,n,&v,&vv,&acc,i,WILLIAMS_METHOD);
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

				zLEGCD(&acc,n,f);
				if (zCompare(f,&zOne) > 0)
				{
					if (zCompare(f,n) == 0)
						zCopy(&zZero,f);
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
		stop2 = clock();
		tt2 = (double)(stop2-start2)/(double)CLOCKS_PER_SEC;

		zLEGCD(&acc,n,f);
		if (zCompare(f,&zOne) > 0)
		{
			if (zCompare(f,n) == 0)
				zCopy(&zZero,f);
			exitcode=1;
			stop2 = clock();
			tt2 = (double)(stop2-start2)/(double)CLOCKS_PER_SEC;
			goto free;
		}

		if (save)
			save_state(STG2,n,&v,&vv,&acc,i,WILLIAMS_METHOD);

	free:

		if (VFLAG >= 0)
			printf("\n");

		for (j=0;j<D;j++)
			zFree(&t[j]);
		free(t);

		zFree(&mTwo);
		zFree(&v);
		zFree(&vv);
		zFree(&s);
		zFree(&v_ad);
		zFree(&acc);
		zFree(&t1);
		zFree(&t2);
		free(str);
		monty_free();
		if (exitcode == 0)
			zCopy(&zZero,f);

		return i;
	}

	void mNextV(z *h, uint32 j,z *mTwo, z *n)
	{
		//compute v_j mod p using the duplication and addition formula's
		//for lucas V sequences.  Algorithm described by Mersenne Wiki entry.
		//montgomery ladder powering algorithm...

		uint8 t,bi[32];
		int i;
		z w1,w2,w3;

		zInit(&w1);
		zInit(&w2);
		zInit(&w3);
		
		//x=h=Vn 
		zCopy(h,&w1);

		//y=(B^2-2) mod N
		monty_sqr(h,&w3,n);
		monty_sub(&w3,mTwo,&w2,n);

		//get j in binary form
		i=0;
		while (j>0)
		{
			i++;
			bi[i] = j % 2;
			j = j/2;
		}
		t=i;

		//compute loop
		//for each bit of M to the right of the most significant bit
		for(i=t-1;i>=1;i--)
		{
			//if the bit is 1
			if (bi[i])
			{
				//x=(x*y-B) mod N 
				monty_mul(&w1,&w2,&w3,n);
				monty_sub(&w3,h,&w1,n);
				//y=(y^2-2) mod N 
				monty_sqr(&w2,&w3,n);
				monty_sub(&w3,mTwo,&w2,n);
			}
			else
			{
				//y=(x*y-B) mod N 
				monty_mul(&w1,&w2,&w3,n);
				monty_sub(&w3,h,&w2,n);
				//x=(x^2-2) mod N 
				monty_sqr(&w1,&w3,n);
				monty_sub(&w3,mTwo,&w1,n);
			}
		}
		zCopy(&w1,h);

		zFree(&w1);
		zFree(&w2);
		zFree(&w3);

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

	} ecm_pp1_data_t;

	ecm_pp1_data_t pp1_data;

	void pp1_init()
	{
		mpz_init(pp1_data.gmp_n);
		mpz_init(pp1_data.gmp_factor);
		ecm_init(pp1_data.params);
		//gmp_randseed_ui(tdata->params->rng, 
		//	get_rand(&obj->seed1, &obj->seed2));

		pp1_data.params->method = ECM_PP1;
		//pp1_data.params->verbose = 1;

		TMP_STG2_MAX = WILL_STG2_MAX;

		PP1_ABORT = 0;
		signal(SIGINT,pp1exit);

		return;
	}

	void pp1_finalize()
	{
		ecm_clear(pp1_data.params);
		mpz_clear(pp1_data.gmp_n);
		mpz_clear(pp1_data.gmp_factor);

		WILL_STG2_MAX = TMP_STG2_MAX;
		signal(SIGINT,NULL);

		return;
	}

	int mwilliams(z *n, z *f, uint32 rp)
	{
		int status;
		size_t count;

		pp1_data.params->B1done = 1.0 + floor (1 * 128.) / 134217728.;
		if (VFLAG >= 3)
			pp1_data.params->verbose = VFLAG - 2;		

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		mpz_import(pp1_data.gmp_n, (size_t)(abs(n->size)), -1, sizeof(uint32), 
			0, (size_t)0, n->val);
#else
		//wrapper for YAFU bigints and call to gmp-ecm
		mp2gmp(n, pp1_data.gmp_n);
#endif

		if (PP1_STG2_ISDEFAULT == 0)
		{
			//not default, tell gmp-ecm to use the requested B2
			//printf("using requested B2 value\n");
			sp642z(WILL_STG2_MAX,f);
			mp2gmp(f,pp1_data.params->B2);
			zClear(f);
		}

		status = ecm_factor(pp1_data.gmp_factor, pp1_data.gmp_n,
				WILL_STG1_MAX, pp1_data.params);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		zClear(n);
		mpz_export(n->val, &count, -1, sizeof(uint32),
				0, (size_t)0, pp1_data.gmp_n);
		n->size = count;
#else
		//update n: not sure if gmp-ecm modifies it
		gmp2mp(pp1_data.gmp_n, n);
#endif

		//NOTE: this required a modification to the GMP-ECM source code in pp1.c
		//in order to get the automatically computed B2 value out of the
		//library
		//gmp2mp(pp1_data.params->B2,f);
		//WILL_STG2_MAX = z264(f);

#if defined(_WIN64) && BITS_PER_DIGIT == 32
		zClear(f);
		mpz_export(f->val, &count, -1, sizeof(uint32),
				0, (size_t)0, pp1_data.gmp_factor);
		f->size = count;
#else
		//pull out any factor found
		gmp2mp(pp1_data.gmp_factor,f);
#endif

		//the return value is the stage the factor was found in, if no error
		pp1_data.stagefound = status;

		return status;
	}

#endif

void williams_loop(int trials, fact_obj_t *fobj)
{
	//use william's p+1 algorithm 'trials' times on n.
	//we allow multiple trials since we may not always get a p+1 group
	//with our random choice of base
	z *n = &fobj->pp1_obj.n;
	z d,f,t;
	int i,it;
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

	//initialize the flag to watch for interrupts, and set the
	//pointer to the function to call if we see a user interrupt
	PP1_ABORT = 0;
	signal(SIGINT,pp1exit);

	zInit(&d);
	zInit(&f);
	zInit(&t);

	pp1_init();

	i=0;
	while (i<trials)
	{
		//watch for an abort
		if (PP1_ABORT)
		{
			print_factors(fobj);
			exit(1);
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
			break;
		}
		
		base = spRand(3,MAX_DIGIT);

		pp1_print_B1_B2(n,flog);

#if defined(HAVE_GMP) && defined(HAVE_GMP_ECM)
		//need a new line to make screen output look right, when
		//using GMP-ECM, because the "processed" status is not printed
		if (VFLAG >= 0)
			printf("\n");
#endif

		it = mwilliams(n,&f,base);
		
		if (zCompare(&f,&zOne) > 0 && zCompare(&f,n) < 0)
		{
			//non-trivial factor found
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;
			if (isPrime(&f))
			{
				f.type = PRP;
				add_to_factor_list(fobj, &f);
				logprint(flog,"prp%d = %s\n",
					ndigits(&f),z2decstr(&f,&gstr2));
			}
			else
			{
				f.type = COMPOSITE;
				add_to_factor_list(fobj, &f);
				logprint(flog,"c%d = %s\n",
					ndigits(&f),z2decstr(&f,&gstr2));
			}
			start = clock();

			//reduce input
			zDiv(n,&f,&t,&d);
			zCopy(&t,n);
			i++;
			break;
		}

		i++;
	}

	fclose(flog);

	pp1_finalize();

	signal(SIGINT,NULL);
	zFree(&d);
	zFree(&f);
	zFree(&t);
	return;
}

void pp1_print_B1_B2(z *n, FILE *flog)
{
	char suffix;
	char stg1str[20];
	char stg2str[20];

	if (WILL_STG1_MAX % 1000000000 == 0)
	{
		suffix = 'B';
		sprintf(stg1str,"%u%c",WILL_STG1_MAX / 1000000000, suffix);
	}
	else if (WILL_STG1_MAX % 1000000 == 0)
	{
		suffix = 'M';
		sprintf(stg1str,"%u%c",WILL_STG1_MAX / 1000000, suffix);
	}
	else if (WILL_STG1_MAX % 1000 == 0)
	{
		suffix = 'K';
		sprintf(stg1str,"%u%c",WILL_STG1_MAX / 1000, suffix);
	}
	else
	{
		sprintf(stg1str,"%u",WILL_STG1_MAX);
	}

	if (PP1_STG2_ISDEFAULT == 0)
	{
		if (WILL_STG2_MAX % 1000000000 == 0)
		{
			suffix = 'B';
			sprintf(stg2str,"%" PRIu64 "%c",WILL_STG2_MAX / 1000000000, suffix);
		}
		else if (WILL_STG2_MAX % 1000000 == 0)
		{
			suffix = 'M';
			sprintf(stg2str,"%" PRIu64 "%c",WILL_STG2_MAX / 1000000, suffix);
		}
		else if (WILL_STG2_MAX % 1000 == 0)
		{
			suffix = 'K';
			sprintf(stg2str,"%" PRIu64 "%c",WILL_STG2_MAX / 1000, suffix);
		}
		else
		{
			sprintf(stg2str,"%" PRIu64 "",WILL_STG2_MAX);
		}
	}
	else
		sprintf(stg2str, "gmp-ecm default");

	if (VFLAG >= 0)
	{
		printf("pp1: starting B1 = %s, B2 = %s on C%d"
			,stg1str,stg2str,ndigits(n));
		fflush(stdout);	
	}
	logprint(flog, "pp1: starting B1 = %s, B2 = %s on C%d\n"
		,stg1str,stg2str,ndigits(n));

	return;
}

void pp1exit(int sig)
{
	printf("\nAborting...\n");
	PP1_ABORT = 1;
	return;
}

