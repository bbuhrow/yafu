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
#include "factor.h"
#include "common.h"
#include "util.h"
#include <gmp.h>

#define NUM_SQUFOF_MULT 16

/*
implements shanks's squfof algorihm.  priceless help from the papers
of Stephen McMath, Danial Shanks, and Jason Gower
*/

typedef struct
{
	uint32 mult;
	uint32 valid;
	uint32 P;
	uint32 bn;
	uint32 Qn;
	uint32 Q0;
	uint32 b0;
	uint32 it;
	uint32 imax;	
} mult_t;

#ifdef NOTDEF
typedef struct
{
	// these have to be allocated on the heap to ensure they are aligned
	uint32 *mult;
	uint32 *valid;
	uint32 *P;
	uint32 *bn;
	uint32 *Qn;
	uint32 *Q0;
	uint32 *b0;
	uint32 *it;
	uint32 *imax;	
} mult_2x_t;

double *twoarray, *fudge_array;
uint32 *rarray, *darray, *narray;
uint8 *sqr_tests, *sqr_tests2;

#endif

void shanks_mult_unit(uint64 N, mult_t *mult_save, uint64 *f);

uint64 sp_shanks_loop(mpz_t N, fact_obj_t *fobj)
{
	//call shanks with multiple small multipliers
	const int multipliers[NUM_SQUFOF_MULT] = {1, 3, 5, 7, 
				11, 3*5, 3*7, 3*11, 
				5*7, 5*11, 7*11, 
				3*5*7, 3*5*11, 3*7*11, 
				5*7*11, 3*5*7*11};
				
	int i, rounds,j;
	uint64 n64, nn64, f64, big1, big2;
	mult_t mult_save[NUM_SQUFOF_MULT];
	mpz_t gmptmp;
	
	big1 = 0xFFFFFFFFFFFFFFFFULL;
	big2 = 0x3FFFFFFFFFFFFFFFULL;

	if (mpz_sizeinbase(N,2) > 62)
	{
		printf("N too big (%d bits), exiting...\n", (int)mpz_sizeinbase(N,2));
		return 1;
	}	

	n64 = mpz_get_64(N);

	if (mpz_sizeinbase(N,2) <= 40)
		return LehmanFactor(n64, 3.5, 1, 0.1);

	//default return value
	f64 = 1;

	if (n64 <= 3)
		return n64;

	mpz_init(gmptmp);

	for (i=NUM_SQUFOF_MULT-1;i>=0;i--)
	{
		// can we multiply without overflowing 64 bits?
		if (big2/(uint64)multipliers[i] < n64)
		{
			//this multiplier makes the input bigger than 64 bits
			mult_save[i].mult = multipliers[i];
			mult_save[i].valid = 0;
			continue;
		}

		//form the multiplied input
		nn64 = n64 * (uint64)multipliers[i];

		mult_save[i].mult = multipliers[i];
		mult_save[i].valid = 1;

		//set imax = N^1/4
		mpz_set_64(gmptmp, nn64);
		mpz_sqrt(gmptmp, gmptmp);	
		mult_save[i].b0 = mpz_get_ui(gmptmp);
		mult_save[i].imax = (uint32)sqrt((double)mult_save[i].b0) / 2;

		//set up recurrence
		mult_save[i].Q0 = 1;
		mult_save[i].P = mult_save[i].b0;
		mult_save[i].Qn = (uint32)(nn64 - 
			(uint64)mult_save[i].b0 * (uint64)mult_save[i].b0);
			
		if (mult_save[i].Qn == 0)
		{
			//N is a perfect square
			f64 = (uint64)mult_save[i].b0;
			goto done;
		}
		mult_save[i].bn = (mult_save[i].b0 + mult_save[i].P)
			/ mult_save[i].Qn;
		mult_save[i].it = 0;

	}

	//now process the multipliers a little at a time.  this allows more
	//multipliers to be tried in order to hopefully find one that 
	//factors the input quickly
	rounds = 6;
	for (i = 0; i < rounds; i++)
	{
		for (j=0; j < NUM_SQUFOF_MULT; j++)
		{
			if (mult_save[j].valid == 0)
				continue;

			//form the input
			nn64 = n64 * multipliers[j];
			//try to factor
			shanks_mult_unit(nn64,&mult_save[j],&f64);

			//check the output for a non-trivial factor
			if (f64 == (uint64)-1)
			{
				//this is an error condition, stop processing this multiplier
				mult_save[j].valid = 0;
			}
			else if (f64 > 1)
			{
				if (f64 != multipliers[j])
				{
					//factor found.  check for and remove small multiplier if necessary.
					nn64 = gcd64(f64,multipliers[j]);
					f64 /= nn64;

					if (f64 != 1)
					{
						//found a non-trivial factor, return it;
						goto done;
					}
					else
					{
						//found trivial factor, stop processing this multiplier
						mult_save[j].valid = 0;
					}
				}
				else
				{
					//found trivial factor, stop processing this multiplier
					mult_save[j].valid = 0;
				}
			}
		}
	}
	//default return value
	f64 = 1;

	//if we've got to here, then the number is still unfactored.  returning
	//a value of 1 signifies this
done:

	/*
	gettimeofday(&myTVend, NULL);
	difference = my_difftime (&myTVstart, &myTVend);

	t_time = 
		((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (1)
		printf("Total squfof time = %6.6f seconds.\n",t_time);
		*/

	mpz_clear(gmptmp);

	return f64;
}

void shanks_mult_unit(uint64 N, mult_t *mult_save, uint64 *f)
{
	//use shanks SQUFOF to factor N.  almost all computation can be done with longs
	//input N < 2^63
	//return 1 in f if no factor is found
	uint32 imax,i,Q0,b0,Qn,bn,P,bbn,Ro,S,So,t1,t2;
	int j=0;	

	//initialize output
	*f=0;

	//load previous save point
	P = mult_save->P;
	bn = mult_save->bn;
	Qn = mult_save->Qn;
	Q0 = mult_save->Q0;
	b0 = mult_save->b0;
	i = mult_save->it;
	imax = i + mult_save->imax;
	
	while (1)
	{
		j=0;

		//i must be even on entering the unrolled loop below
		if (i & 0x1)
		{
			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q
			bn = (b0 + P) / Qn;
			i++;
		}

		while (1)
		{
			//at the start of every iteration, we need to know:
			//	P from the previous iteration
			//	bn from the previous iteration
			//	Qn from the previous iteration
			//	Q0 from the previous iteration
			//	iteration count, i
			if (i >= imax)
			{
				//haven't progressed to the next stage yet.  let another
				//multiplier try for awhile.  save state so we
				//know where to start back up in the next round
				mult_save->P = P;
				mult_save->bn = bn;
				mult_save->Qn = Qn;
				mult_save->Q0 = Q0;
				mult_save->it = i;
				//signal to do nothing but continue looking for a factor
				*f = 0;	
				return;
			}

			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q

			// this happens fairly often, but special-casing the division
			// still makes things slower (at least on modern intel chips)
			//if (Qn > ((b0 + P) >> 1))
			//	bn = 1;
			//else
			bn = (b0 + P) / Qn;		
			i++;

			//even iteration
			//check for square Qn = S*S
			t2 = Qn & 31;
			if (t2 == 0 || t2 == 1 || t2 == 4 ||
				t2 == 9 || t2 == 16 || t2 == 17 || t2 == 25)			
			{
				// extra squaritude tests are also slower on modern intel chips
				//t2 = Qn & 63;
				//if (t2 < 32 || t2 == 33 || t2 == 36 || 
				//	t2 == 41 ||t2 == 49 || t2 == 57)
				//{
					t1 = (uint32)sqrt(Qn);
					if (Qn == t1 * t1)
						break;
				//}
			}

			//odd iteration
			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q

			bn = (b0 + P) / Qn;
			i++;
					
		}

		//reduce to G0
		S = (int)sqrt(Qn);
		Ro = P + S*((b0 - P)/S);
		t1 = Ro;
		So = (uint32)(((int64)N - (int64)t1*(int64)t1)/(int64)S);
		bbn = (b0+Ro)/So;

		//search for symmetry point
		while (1)
		{
			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

		}

		*f = gcd64(Ro,N);
		//found a factor - don't know yet if it's trivial or not.
		//we don't need to remember any state info, as one way or the other
		//this multiplier will be invalidated		
		if (*f > 1)
			return;
	}
}

#ifdef NOTDEF
void shanks_mult_unit_2x(uint64 N, mult_2x_t *mult_save, int start_mult, uint64 *f);

uint64 sp_shanks_loop_2x(mpz_t N, fact_obj_t *fobj)
{
	//call shanks with multiple small multipliers
	const int multipliers[NUM_SQUFOF_MULT] = {1, 3, 5, 7, 
				11, 3*5, 3*7, 3*11, 
				5*7, 5*11, 7*11, 
				3*5*7, 3*5*11, 3*7*11, 
				5*7*11, 3*5*7*11};
				
	int i, rounds, j;
	uint64 n64, nn64, f64, big1, big2;
	mult_2x_t mult_save;
	mpz_t gmptmp;
	
	big1 = 0xFFFFFFFFFFFFFFFFULL;
	big2 = 0x3FFFFFFFFFFFFFFFULL;

	if (mpz_sizeinbase(N,2) > 62)
	{
		printf("N too big (%d bits), exiting...\n", (int)mpz_sizeinbase(N,2));
		return 1;
	}

	n64 = mpz_get_64(N);

	//default return value
	f64 = 1;

	if (n64 <= 3)
		return n64;

	mult_save.b0 = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.bn = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.imax = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.it = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.mult = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.P = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.Q0 = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.Qn = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));
	mult_save.valid = (uint32 *)xmalloc_align(NUM_SQUFOF_MULT * sizeof(double));

	mpz_init(gmptmp);

	twoarray = (double *)xmalloc_align(2 * sizeof(double));
	fudge_array = (double *)xmalloc_align(2 * sizeof(double));
	rarray = (uint32 *)xmalloc_align(2 * sizeof(uint32));
	darray = (uint32 *)xmalloc_align(2 * sizeof(uint32));
	narray = (uint32 *)xmalloc_align(2 * sizeof(uint32));
	twoarray[0] = 2.0;
	twoarray[1] = 2.0;
	fudge_array[0] = 1.0/65536.0;
	fudge_array[1] = 1.0/65536.0;

	sqr_tests = (uint8 *)xmalloc_align(16 * sizeof(uint8));
	sqr_tests2 = (uint8 *)xmalloc_align(16 * sizeof(uint8));
	sqr_tests[0] = 0;
	sqr_tests[1] = 1;
	sqr_tests[2] = 4;
	sqr_tests[3] = 9;
	sqr_tests[4] = 0;
	sqr_tests[5] = 1;
	sqr_tests[6] = 4;
	sqr_tests[7] = 9;
	sqr_tests[8] = 0;
	sqr_tests[9] = 1;
	sqr_tests[10] = 4;
	sqr_tests[11] = 9;
	sqr_tests[12] = 0;
	sqr_tests[13] = 1;
	sqr_tests[14] = 4;
	sqr_tests[15] = 9;

	sqr_tests2[0] = 16;
	sqr_tests2[1] = 17;
	sqr_tests2[2] = 25;
	sqr_tests2[3] = 255;
	sqr_tests2[4] = 16;
	sqr_tests2[5] = 17;
	sqr_tests2[6] = 25;
	sqr_tests2[7] = 255;
	sqr_tests2[8] = 16;
	sqr_tests2[9] = 17;
	sqr_tests2[10] = 25;
	sqr_tests2[11] = 255;
	sqr_tests2[12] = 16;
	sqr_tests2[13] = 17;
	sqr_tests2[14] = 25;
	sqr_tests2[15] = 255;	

	j = 0;	//num valid multipliers
	for (i=NUM_SQUFOF_MULT-1;i>=0;i--)
	{
		// can we multiply without overflowing 64 bits?
		if (big2/(uint64)multipliers[i] < n64)
		{
			//this multiplier makes the input bigger than 64 bits
			mult_save.mult[i] = (uint32)multipliers[i];
			mult_save.valid[i] = 0;
			continue;
		}

		//form the multiplied input
		nn64 = n64 * (uint64)multipliers[i];

		mult_save[i].mult = multipliers[i];
		mult_save[i].valid = 1;

		//set imax = N^1/4
		mpz_set_64(gmptmp, nn64);
		mpz_sqrt(gmptmp, gmptmp);	
		mult_save.b0[i] = mpz_get_ui(gmptmp);
		mult_save.imax[i] = (uint32)sqrt((double)mult_save.b0[i]) / 2;

		//set up recurrence
		mult_save.Q0[i] = 1;
		mult_save.P[i] = mult_save.b0[i];
		mult_save.Qn[i] = (uint32)(nn64 - 
			(uint64)mult_save.b0[i] * (uint64)mult_save.b0[i]);
			
		if (mult_save.Qn[i] == 0)
		{
			//N is a perfect square
			f64 = (uint64)mult_save.b0[i];
			goto done;
		}
		mult_save.bn[i] = (mult_save.b0[i] + mult_save.P[i])
			/ mult_save.Qn[i];
		mult_save.it[i] = 0;

		j++;
	}		

	//now process the multipliers a little at a time in batches of 4.  this allows more
	//multipliers to be tried in order to hopefully find one that 
	//factors the input quickly
	rounds = 6;
	for (i = 0; i < rounds; i++)
	{
		for (j=0; j < NUM_SQUFOF_MULT; j++)
		{
			if (mult_save.valid[j] == 0)
				continue;

			//form the input
			nn64 = n64 * multipliers[j];
			//try to factor
			shanks_mult_unit_2x(nn64, &mult_save, j, &f64);

			//check the output for a non-trivial factor
			if (f64 == (uint64)-1)
			{
				//this is an error condition, stop processing this multiplier
				mult_save.valid[j] = 0;
			}
			else if (f64 > 1)
			{
				if (f64 != multipliers[j])
				{
					//factor found.  check for and remove small multiplier if necessary.
					nn64 = gcd64(f64,multipliers[j]);
					f64 /= nn64;

					if (f64 != 1)
					{
						//found a non-trivial factor, return it;
						goto done;
					}
					else
					{
						//found trivial factor, stop processing this multiplier
						mult_save.valid[j] = 0;
					}
				}
				else
				{
					//found trivial factor, stop processing this multiplier
					mult_save.valid[j] = 0;
				}
			}
		}
	}
	//default return value
	f64 = 1;

	//if we've got to here, then the number is still unfactored.  returning
	//a value of 1 signifies this
done:

	mpz_clear(gmptmp);

	align_free(mult_save.b0);
	align_free(mult_save.bn);
	align_free(mult_save.imax);
	align_free(mult_save.it);
	align_free(mult_save.mult);
	align_free(mult_save.P);
	align_free(mult_save.Q0);
	align_free(mult_save.Qn);
	align_free(mult_save.valid);

	align_free(twoarray);
	align_free(fudge_array);
	align_free(rarray);
	align_free(darray);
	align_free(narray);
	align_free(sqr_tests);
	align_free(sqr_tests2);

	return f64;
}

//typedef struct
//{
//	// these have to be allocated on the heap to ensure they are aligned
//	uint32 *mult;		//0
//	uint32 *valid;		//8
//	uint32 *P;			//16
//	uint32 *bn;			//24
//	uint32 *Qn;			//32
//	uint32 *Q0;			//40
//	uint32 *b0;			//48
//	uint32 *it;			//56
//	uint32 *imax;		//64
//} mult_2x_t;

void shanks_mult_unit_2x(uint64 N, mult_2x_t *mult_save, int start_mult, uint64 *f)
{
	//use shanks SQUFOF to factor N.  almost all computation can be done with longs
	//input N < 2^63
	//return 1 in f if no factor is found
	uint32 imax,i,Q0,b0,Qn,bn,P,bbn,Ro,S,So,t1,t2;

	//initialize output
	*f=0;

	// loop through all multipliers, regardless of if they are valid or not.
	// later, add code to get smarter about which multipliers are used.
	asm(
		"movq	16(%0,1), %%rax  \n\t"				/* address of P */
		"movdqa	(%%rax, %%rcx, 4), %%xmm0 \n\t"		/* xmm0 := values of P */
		"movq	24(%0,1), %%rax  \n\t"				/* address of bn */
		"movdqa	(%%rax, %%rcx, 4), %%xmm1 \n\t"		/* xmm1 := values of bn */
		"movq	32(%0,1), %%rax  \n\t"				/* address of Qn */
		"movdqa	(%%rax, %%rcx, 4), %%xmm2 \n\t"		/* xmm2 := values of Qn */
		"movq	40(%0,1), %%rax  \n\t"				/* address of Q0 */
		"movdqa	(%%rax, %%rcx, 4), %%xmm3 \n\t"		/* xmm3 := values of Q0 */
		"movq	48(%0,1), %%rax  \n\t"				/* address of b0 */
		"movdqa	(%%rax, %%rcx, 4), %%xmm4 \n\t"		/* xmm4 := values of b0 */
		"movq	56(%0,1), %%rax  \n\t"				/* address of i */
		"movdqa	(%%rax, %%rcx, 4), %%xmm5 \n\t"		/* xmm5 := values of i */
		"0:	/n/t"									/* top of loop */
		"cmpl	%%ebx, %%ecx \n\t"
		"jge	1f	\n\t"							/* jump out of loop if i >= imax */
		"movdqa	%%xmm0, %%xmm6 \n\t"				/* store tmp P */
		"movdqa	%%xmm2, %%xmm7 \n\t"				/* move Qn for mul */
		"movdqa	%%xmm1, %%xmm9 \n\t"				/* move bn for mul */
		"movdqa	%%xmm2, %%xmm8 \n\t"				/* move Qn for mul */
		"movdqa	%%xmm1, %%xmm10 \n\t"				/* move bn for mul */
		"psrldq	$4, %%xmm8, %%xmm8 \n\t"
		"psrldq	$4, %%xmm10, %%xmm10 \n\t"
		"pmuludq	%%xmm7, %%xmm9 \n\t"			/* words 1 and 3 */
		"pmuludq	%%xmm8, %%xmm10 \n\t"			/* words 2 and 4 */
		"pslldq	$4, %%xmm10, %%xmm10 \n\t"
		"por	%%xmm9, %%xmm10 \n\t"				/* words 1 - 4 */


		"addl	$1, %%ebx \n\t"
		"jmp	0b \n\t"
		"1:	/n/t"
		:
		: "r"(mult_save), "b"(i), "c"(start_mult), "d"(imax), "r"(success)
		: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");

	
	while (1)
	{
		//i must be even on entering the unrolled loop below
		if (i & 0x1)
		{
			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q
			bn = (b0 + P) / Qn;
			i++;
		}

		while (1)
		{
			//at the start of every iteration, we need to know:
			//	P from the previous iteration
			//	bn from the previous iteration
			//	Qn from the previous iteration
			//	Q0 from the previous iteration
			//	iteration count, i
			if (i >= imax)
			{
				//haven't progressed to the next stage yet.  let another
				//multiplier try for awhile.  save state so we
				//know where to start back up in the next round
				mult_save->P = P;
				mult_save->bn = bn;
				mult_save->Qn = Qn;
				mult_save->Q0 = Q0;
				mult_save->it = i;
				//signal to do nothing but continue looking for a factor
				*f = 0;	
				return;
			}

			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q

			// this happens fairly often, but special-casing the division
			// still makes things slower (at least on modern intel chips)
			//if (Qn > ((b0 + P) >> 1))
			//	bn = 1;
			//else
			bn = (b0 + P) / Qn;		
			i++;

			//even iteration
			//check for square Qn = S*S
			//idea: broadcast the results of Qn & x to each of several bytes in
			//an sse2 register.  then pcmpeqb with a pre-set sse2 register containing
			//all of the valid square endings (below) followed by pmovmskb and check 
			//for non-zero.  2 or 3 clocks to accomplish all of the tests below, or more...
#ifdef NOTDEF
			t2 = Qn & 31;
			if (t2 == 0 || t2 == 1 || t2 == 4 ||
				t2 == 9 || t2 == 16 || t2 == 17 || t2 == 25)			
			{
				// extra squaritude tests are also slower on modern intel chips
				//t2 = Qn & 63;
				//if (t2 < 32 || t2 == 33 || t2 == 36 || 
				//	t2 == 41 ||t2 == 49 || t2 == 57)
				//{
					t1 = (uint32)sqrt(Qn);
					if (Qn == t1 * t1)
						break;
				//}
			}
#else
			//this is almost as fast.  the broadcast is the bottleneck.
			t2 = ((Qn & 63) << 8) | (Qn & 63);
			t2 = (t2 << 16) | t2;
			asm(
				"movd	%%eax, %%xmm0 \n\t"
				"pshufd	$0, %%xmm0, %%xmm1 \n\t"
				"pcmpeqb	%%xmm5, %%xmm1 \n\t"
				"pmovmskb	%%xmm1, %0 \n\t"
				: "=r"(t1)
				: "a"(t2)
				: "xmm0", "xmm1", "memory", "cc");

			if (t1 != 0)
			{
				t2 = (uint32)sqrt(Qn);
				if (Qn == t2 * t2)
					break;
			}

#endif

			//odd iteration
			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn*(t1-P);
			Q0 = t2;	//remember last Q

			bn = (b0 + P) / Qn;
			i++;
					
		}

		//reduce to G0
		S = (int)sqrt(Qn);
		Ro = P + S*((b0 - P)/S);
		t1 = Ro;
		So = (uint32)(((int64)N - (int64)t1*(int64)t1)/(int64)S);
		bbn = (b0+Ro)/So;

		//search for symmetry point
		while (1)
		{
			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

			t1 = Ro;		//hold Ro for this iteration
			Ro = bbn*So - Ro;
			t2 = So;		//hold So for this iteration
			So = S + bbn*(t1-Ro);
			S = t2;			//remember last S
			bbn = (b0+Ro)/So;
			
			//check for symmetry point
			if (Ro == t1)
				break;

		}

		*f = gcd64(Ro,N);
		//found a factor - don't know yet if it's trivial or not.
		//we don't need to remember any state info, as one way or the other
		//this multiplier will be invalidated		
		if (*f > 1)
			return;
	}
}
#endif

/****
Code by Warren D. Smith, Dec 2011, to implement Lehman integer-factorization
algorithm.  Rigorous O(N^(1/3)) step factoring
(or prime-proving) algorithm invented by RS Lehman 1974.
Compile:
 gcc LehmanClean.c -Wall -O6 -o lehcl
 gcc LehmanClean.c -Wall -DNDEBUG -O6 -o lehcl
*******/

static double sqr_tab[1024]; 

void make_sqr_tab()
{
	int i;
	for(i=0; i<1024; i++)
		sqr_tab[i] = sqrt((double)i);

	return;
}

static unsigned char issq1024[1024];
static unsigned char issq4199[4199];

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 1
#endif

void MakeIssq()
{
	int i;
	for(i=0; i<1024; i++){ issq1024[i] = 0; }
	for(i=0; i<1024; i++){ issq1024[(i*i)%1024] = 1; }
	for(i=0; i<4199; i++){ issq4199[i] = 0; }
	for(i=0; i<3465; i++){ issq4199[(i*i)%3465] |= 2; }
	for(i=0; i<4199; i++){ issq4199[(i*i)%4199] |= 1; }
	//printf("Square tables built.\n");
}

//the 6542 primes up to 65536=2^16, then sentinel 65535 at end
static uint16 prime[6543]; 

void MakePrimeTable(){
	uint32 i,j,k;
	prime[0]=2;
	prime[1]=3;
	prime[2]=5;
	k=3;
	for(i=7; i<65536; i+=2){
		for(j=0; prime[j]*(uint32)prime[j] <= i; j++){
			if(i%prime[j]==0) goto NONPRIME;
		}
		prime[k]=i;
		k++;
		NONPRIME: ;
	}

	prime[k] = 65535; //sentinel
	//printf("Prime table[0..%d] built: ", k);
	//for(i=0; i<20; i++){ printf("%d,", prime[i]); }
	//printf("%d,...,%d,(%d)\n", prime[20],prime[6541],prime[6542]);
}

/*** LehmanFactor  Testing:
Experimentally: Always returns N if N=prime<65536 and DoTrial=FALSE
and 0.1<=LehTune<=9.6.
  But: If LehTune<=0.8 then the   gcd(a+b, N)   can return N as its
mode of success if N=3.
  So the theorem that the special prime-detect return always will
fire for N prime can thus
  be violated if LehTune<0.8 and N=3 (still correctly returns N though), but
       valid if LehTune>=0.8  and N>2.
  Also valid if LehTune>=0.25 and N>31.
  Also valid if LehTune>=0.10 and N>83.

High end: if N < 2^43.000049 = 8796393022207  then apparently no
overflow issues.
Later note: sorry, the simple search loop I was using was too
simplistic a search approach;
The currently least known overflow-failure also occurs at   N =
5132209842943 = 2^42.2227.

Factor Validity: testing 100 million N with N=A*B with
 2611953 <= A <= B <= 2965920
LehmanFactor always succeeded in factoring.
However, this is deceptive, see above on weird rare failures.  It
appears that Lehman
can be sensitive to the distribution of N and there are "hard locales"
and "easy locales" which
are not obvious.  Also changing the tuning constants can alter the failure set.
***/

uint64 LehmanFactor(uint64 N, double Tune, int DoTrial, double CutFrac)
{
	uint32 b,p,k,r,B,U,Bred,inc,FirstCut,ip = 1;
	uint64 a,c,kN,kN4,B2;
	double Tune2, Tune3, x, sqrtn;
	mpz_t tmpz;

	if((N&1)==0) 
		return(2); //N is even

	if(Tune<0.1)
	{
		printf("Sorry0, Lehman only implemented for Tune>=0.1\n");
		return(1);
	}

	mpz_init(tmpz);
	mpz_set_64(tmpz, N);
	mpz_root(tmpz, tmpz, 3);
	B = Tune * (1+(double)mpz_get_ui(tmpz));

	FirstCut = CutFrac*B;

	//assures prime N will not activate "wrong" Lehman return
	if(FirstCut<84)
		FirstCut=84; 

	if(FirstCut>65535)
		FirstCut = 65535; 

	if(DoTrial)
	{
		for(ip=1; ; ip++)
		{                   
			p = prime[ip];
			if(p>=FirstCut) 
				break;
			if(N%p==0) 
				return(p);
		}
	}

	if(N>=8796393022207ull)
	{
		printf("Sorry1, Lehman only implemented for N<8796393022207\n");
		return(1);
	}

	Tune2 = Tune*Tune;
	Tune3 = Tune2*Tune;
	Bred = B / Tune3;

	B2 = B*B;
	kN = 0;

	//Lehman suggested (to get more average speed) trying highly-divisible k first. However,
	//my experiments on trying to to that have usually slowed things down versus this simple loop:
	sqrtn = sqrt((double)N);
	for(k=1; k<=Bred; k++)
	{
		if(k&1)
		{ 
			inc=4; 
			r=(k+N)%4; 
		}
		else
		{ 
			inc=2; 
			r=1; 
		} 

		kN += N;
		if(kN >= 1152921504606846976ull)
		{
			printf("Sorry2, overflow, N=%" PRIu64 " is too large\n", N);
			return(1);
		}

		//Actually , even if overflow occurs here, one could still use approximate
		//arithmetic to compute kN4's most-signif 64 bits only, then still exactly compute x...
		//With appropriate code alterations is should be possible to extent the range... but
		//I have not tried that idea for trying to "cheat" to gain more precision.
		kN4 = kN*4;
		if (k < 1024)
			x = sqrtn * sqr_tab[k];
		else
			x = sqrt((double)kN);

		a = x;
		if((uint64)a*(uint64)a==kN)
		{ 
			B2 = gcd64((uint64)a, N);
			return(B2);
		}

		x *= 2;
		a = x+0.9999999665; //very carefully chosen.
		//Let me repeat that: a = x+0.9999999665.  Really.
		b=a%inc;  
		b = a + (inc+r-b)%inc;   //b is a but adjusted upward to make b%inc=r.
		c = (uint64)b*(uint64)b - kN4;  //this is the precision bottleneck.
		//At this point, I used to do a test:
		//if( c+kN4 != (uint64)b*(uint64)b ) //overflow-caused failure: exit!
		//	printf("Sorry3, unrepairable overflow, N=%llu is too large\n", N);
		//  return(0);
		//However, I've now reconsidered.  I claim C language computes c mod 2^64 correctly.
		//If overflow happens, this is irrelevant because c is way smaller than 2^64, kN4, and b*b.
		//Hence c should be correct, despite the overflow. Hence I am removing this error-exit.

		U = x + B2/(2*x); 
		//old code was  U = sqrt((real)(B2+kN4+0.99));   and was 4% slower.

		//Below loop is: for(all integers a with 0<=a*a-kN4<=B*B and with a%inc==r)
		for(a=b;  a<=U;  c+=inc*(a+a+inc), a+=inc )
		{
			//again, even though this assert can fail due to overflow, that overflow should not matter:
			/** Programming trick:    naive code:     c = a*a-kN4;
			In the inner loop c is bounded between 0 and T^2*N^(2/3)
			and can be updated additively by adding inc*(anew+aold) to it
			when we update a to anew=a+inc. This saves a multiplication in
			the inner loop and/or allows us to reduce precision. **/
			if(issq1024[c&1023])
			{
				if(issq4199[c%3465]&2)
				{
					if(issq4199[c%4199]&1)
					{
						b = sqrt(c + 0.9);
						if(b*b==c)
						{ 
							//square found
							B2 = gcd64((uint64)(a+b), N);
							if(B2>=N)
								printf("theorem failure: B2=%" PRIu64 " N=%" PRIu64 "\n", B2,N); 
							return(B2);
						}
					}
				}
			}
		}
	}

	//square-finding has failed so resume missing part of trial division: 
	if(DoTrial)
	{ 
		if(B>65535) 
			B = 65535;
		for( ; ; ip++)
		{ 
			p = prime[ip];
			if(p>=B) 
				break;
			if(N%p==0) 
				return(p);
		}
	}

	return(N); //N is prime
}

void init_lehman()
{
	MakeIssq();
	MakePrimeTable();
	make_sqr_tab();

	return;
}

