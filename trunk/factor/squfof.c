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
	//FILE *out;
	
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
			bn = (b0 + P) / Qn;
			i++;

			//even iteration
			//check for square Qn = S*S
			t2 = Qn & 31;
			if (t2 == 0 || t2 == 1 || t2 == 4 ||
				t2 == 9 || t2 == 16 || t2 == 17 || t2 == 25)
			{
				t1 = (uint32)sqrt(Qn);
				if (Qn == t1 * t1)
					break;
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


