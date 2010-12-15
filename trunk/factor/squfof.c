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

#define NUM_SQUFOF_MULT 16
#define SQUFOF_QSIZE 50

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

typedef struct
{
	uint64 mult;
	uint64 valid;
	uint64 P;
	uint64 bn;
	uint64 Qn;
	uint64 Q0;
	uint64 b0;
	uint64 it;
	uint64 imax;
} mult_big_t;

void shanks_mult_unit(uint64 N, mult_t *mult_save, uint64 *f);
void shanks_mult_unit_parallel(uint64 N, mult_t *mult_save, 
							   mult_t *mult_save2, uint64 *f, uint64 *f2);
uint64 sp_shanks_loop_big(z *N, fact_obj_t *fobj);
void shanks_mult_unit_big(z *N, mult_big_t *mult_save, uint64 *f);


uint64 sp_shanks_loop(z *N, fact_obj_t *fobj)
{
	//call shanks with multiple small multipliers
	const int multipliers[NUM_SQUFOF_MULT] = {1, 3, 5, 7, 
				11, 3*5, 3*7, 3*11, 
				5*7, 5*11, 7*11, 
				3*5*7, 3*5*11, 3*7*11, 
				5*7*11, 3*5*7*11};
				
	int i, rounds,j;
	uint64 n64, nn64, f64, big1, big2;
	mult_t mult_save[NUM_SQUFOF_MULT];//, mult_save2;
	z tmp1, tmp2;
	//struct timeval myTVstart, myTVend;
	//TIME_DIFF *	difference;
	//double t_time;

	//otherwise mingw complains
	big1 = ((uint64)0xFFFFFFFF << 32) | (uint64)0xFFFFFFFF;
	big2 = ((uint64)0x3FFFFFFF << 32) | (uint64)0xFFFFFFFF;

	if (zBits(N) > 62)
	{
		if (zBits(N) < 100)
			return sp_shanks_loop_big(N,fobj);

		printf("N too big (%d bits), exiting...\n",zBits(N));
		return 1;
	}

	n64 = z264(N);

	//default return value
	f64 = 1;

	if (n64 <= 3)
		return n64;

	//gettimeofday(&myTVstart, NULL);

	//prepare to factor using squfof
	zInit(&tmp1);
	zInit(&tmp2);

	for (i=NUM_SQUFOF_MULT-1;i>=0;i--)
	{
		// can we multiply without overflowing 64 bits?
		if (big1/(uint64)multipliers[i] < n64)
		{
			//this multiplier makes the input bigger than 64 bits
			mult_save[i].valid = 0;
			continue;
		}

		//form the multiplied input
		nn64 = n64 * (uint64)multipliers[i];

		// check if resulting number is too big
		if (nn64 < big2)
		{
			//this multiplier is ok
			mult_save[i].mult = multipliers[i];
			mult_save[i].valid = 1;
			//set imax = N^1/4
			//b0 = (uint32)sqrt((double)(N));
			sp642z(nn64,&tmp1);
			zNroot(&tmp1,&tmp2,2);
			mult_save[i].b0 = (uint32)tmp2.val[0];
			mult_save[i].imax = (uint32)sqrt(mult_save[i].b0) / 2;

			//set up recurrence
			mult_save[i].Q0 = 1;
			mult_save[i].P = mult_save[i].b0;
			mult_save[i].Qn = (uint32)(nn64 - 
				(uint64)mult_save[i].b0*(uint64)mult_save[i].b0);
			
			if (mult_save[i].Qn == 0)
			{
				//N is a perfect square
				zFree(&tmp1);
				zFree(&tmp2);
				f64 = (uint64)mult_save[i].b0;
				goto done;
			}
			mult_save[i].bn = (mult_save[i].b0 + mult_save[i].P)
				/ mult_save[i].Qn;
			mult_save[i].it = 0;
		}
		else
		{
			//this multiplier makes the input too big
			mult_save[i].mult = multipliers[i];
			mult_save[i].valid = 0;
		}
	}

	zFree(&tmp1);
	zFree(&tmp2);


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

	return f64;
}

void shanks_mult_unit(uint64 N, mult_t *mult_save, uint64 *f)
{
	//use shanks SQUFOF to factor N.  almost all computation can be done with longs
	//input N < 2^63
	//return 1 in f if no factor is found
	uint32 imax,i,Q0,b0,Qn,bn,P,bbn,Ro,S,So,t1,t2;
	int j=0;
	FILE *out;
	
	//initialize output
	*f=0;

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

		//i must be even on entering the loop below
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
				t1 = (uint32)sqrt((double)Qn);
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

			//j+=4;
			//if (j > 100000000)
			//{
			//	out = fopen("squfof_error.log","a");
			//	fprintf(out,"error on N = 0x%x%x\n",
			//		(uint32)(N >> 32),(uint32)(N&0xFFFFFFFF));
			//	fclose(out);
			//	*f=-1;
			//	//this gets stuck very rarely, but it does happen.
			//	//we don't need to remember any state info, as this will
			//	//invalidate this multiplier
			//	return;		
			//}
		}
	
		*f = gcd64(Ro,N);
		//found a factor - don't know yet if it's trivial or not.
		//we don't need to remember any state info, as this will
		//invalidate this multiplier
		if (*f > 1)
			return;
	}
}

/* Implementation of algorithm explained in Gower and Wagstaff paper */
int SQUFOF_alpertron(int64 N, int64 queue[])
{
double sqrtn;
int B, Q, Q1, P, P1, L, S;
int i, r, s, t, q, u;
int queueHead, queueTail, queueIndex;
    /* Step 1: Initialize */
if ((N & 3) == 1)
{
  N <<= 1;
}
sqrtn = sqrt(N);
S = (int)sqrtn;
if ((long)(S+1)*(long)(S+1)<=N)
{
  S++;
}
if ((long)S*(long)S > N)
{
  S--;
}
if ((long)S*(long)S == N)
{
  return S;
}
Q1 = 1;
P = S;
Q = (int)N - P*P;
L = (int)(2*sqrt(2*sqrtn));
B = L << 1;
queueHead = 0;
queueTail = 0;
   /* Step 2: Cycle forward to find a proper square form */
for (i=0; i<=B; i++)
{
  q = (S+P)/Q;
  P1 = q*Q-P;
  if (Q <= L)
  {
    if ((Q & 1) == 0)
    {
      queue[queueHead++] = Q >> 1;
      queue[queueHead++] = P % (Q >> 1);
      if (queueHead == 100)
      {
        queueHead = 0;
      }
    }
    else if (Q+Q<=L)
    {
      queue[queueHead++] = Q;
      queue[queueHead++] = P % Q;
      if (queueHead == 100)
      {
        queueHead = 0;
      }
    }
  }
  t = Q1+q*(P-P1);
  Q1 = Q;
  Q = t;
  P = P1;
  if ((i & 1) == 0 && ((Q & 7) < 2 || (Q & 7) == 4))
  {
    r = (int)sqrt(Q);
    if (r*r == Q)
    {
      queueIndex = queueTail;
      for (;;)
      {
        if (queueIndex == queueHead)
        {
          /* Step 3: Compute inverse square root of the square form */
          Q1 = r;
          u = (S-P)%r;
          u += (u >> 31) & r;
          P = S-u;
          Q = (int)((N-(long)P*(long)P)/Q1);
             /* Step 4: Cycle in the reverse direction to find a factor of N */
          for (;;)
          {
            q = (S+P)/Q;
            P1 = q*Q-P;
            if (P == P1)
            {
              break;
            }
            t = Q1 +q*(P-P1);
            Q1 = Q;
            Q = t;
            P = P1;
          }
          /* Step 5: Get the factor of N */
          if ((Q & 1) == 0)
          {
            return Q >> 1;
          }
          return Q;
        }
        s = queue[queueIndex++];
        t = queue[queueIndex++];
        if (queueIndex == 100)
        {
          queueIndex = 0;
        }
        if ((P-t)%s == 0)
        {
          break;
        }
      }
      if (r > 1)
      {
        queueTail = queueIndex;
      }
      if (r == 1)
      {
        queueIndex = queueTail;
        for (;;)
        {
          if (queueIndex == queueHead)
          {
            break;
          }
          if (queue[queueIndex] == 1)
          {
            return 0;
          }
          queueIndex += 2;
          if (queueIndex == 100)
          {
            queueIndex = 0;
          }
        }
      }
    }
  }
}
return 0;
}

int qqueue[100];
int qpoint;

void enqu(int q,int *iter)
{
qqueue[qpoint] = q;
if (++qpoint > 100) *iter = -1;
}


int squfof_rds(int64 n,int *fac1,int *fac2)
{   /* start of squfof: factor n as fac1*fac2  faster in FP?????*/
int64 temp, temp1;
register int iq,ll,l2,p,pnext,q,qlast,r,s,t,i;
int jter,iter;

qlast = 1;
s = (int)sqrt(n);

p = s;
temp1 = s * s;
temp = n - temp1;                 /* temp = n - floor(sqrt(n))^2   */
if (temp == 0)
   {                                   /* Here n is a square            */
   *fac1 = s;
   *fac2 = s;
   return(1);
   }

q = (int)temp;              /* q = excess of n over next smaller square */
ll = 1 + 2*(int)sqrt((double)(p+p));
l2 = ll/2;
qpoint = 0;

/*   In the loop below, we need to check if q is a square right before   */
/*  the end of the loop.  Is there a faster way? The current way is      */
/*   EXPENSIVE! (many branches and double prec sqrt)                     */

for (jter=0; jter < 800000; jter++)      /* I see no way to speed this   */
   {                                     /*  main loop                   */
   iq = (s + p)/q;   
   pnext = iq*q - p;
   if (q <= ll)
      {
      if ((q & 1) == 0) enqu(q/2,&jter);
      else if (q <= l2) enqu(q,&jter);
      if (jter < 0)
          {                        
          return 0;
          }
      }
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;                          /* check for square; even iter   */
   if (jter & 1) continue;             /* jter is odd:omit square test  */
   r = (int)sqrt((double)q);                 /* r = floor(sqrt(q))      */
   if (q != r*r) continue;
   if (qpoint == 0) goto gotit;
   for (i=0; i<qpoint-1; i+=2)      /* treat queue as list for simplicity*/
      {
      if (r == qqueue[i]) goto contin;
      if (r == qqueue[i+1]) goto contin;
      }
   if (r == qqueue[qpoint-1]) continue;
   goto gotit;
contin:;   
   }   /* end of main loop */

gotit:   ;
qlast = r;
p = p + r*((s - p)/r);
temp = (int64)p * (int64)p;
temp = n - temp;
temp1 = temp / qlast;
q = (int)temp1;					/* q = (n - p*p)/qlast (div is exact)*/                
for (iter=0; iter<40000; iter++)
   {                              /* begin second main loop            */
   iq = (s + p)/q;                /* unroll it, of course              */
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   iq = (s + p)/q;
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   iq = (s + p)/q;
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   iq = (s + p)/q;
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   }


return(0);                               /* this shouldn't happen      */

gotfac:   ; if ((q & 1) == 0) q/=2;      /* q was factor or 2*factor   */
*fac1 = q;
temp = n / q;
*fac2 = (int)temp;
return(1);
}


uint64 sp_shanks_loop_big(z *N, fact_obj_t *fobj)
{
	//call shanks with multiple small multipliers
	const int multipliers[NUM_SQUFOF_MULT] = {1, 3, 5, 7, 
				11, 3*5, 3*7, 3*11, 
				5*7, 5*11, 7*11, 
				3*5*7, 3*5*11, 3*7*11, 
				5*7*11, 3*5*7*11};
				
	int i, rounds,j;
	uint64 nn64, f64;
	mult_big_t mult_save[NUM_SQUFOF_MULT];
	z tmp1, tmp2, tmp3;
	//struct timeval myTVstart, myTVend;
	//TIME_DIFF *	difference;
	//double t_time;

	if (zBits(N) > 100)
	{
		printf("N too big (%d bits), exiting...\n",zBits(N));
		return 1;
	}


	//default return value
	f64 = 1;

	if (zCompare(N,&zTwo) <= 0)
		return (uint64)N->val[0];

	//gettimeofday(&myTVstart, NULL);

	//prepare to factor using squfof
	zInit(&tmp1);
	zInit(&tmp2);
	zInit(&tmp3);

	for (i=NUM_SQUFOF_MULT-1;i>=0;i--)
	{
		//form the multiplied input
		zShortMul(N,(uint64)multipliers[i],&tmp1);

		mult_save[i].mult = multipliers[i];
		mult_save[i].valid = 1;
		//set imax = N^1/4
		//b0 = (uint32)sqrt((double)(N));
		zNroot(&tmp1,&tmp2,2);
		mult_save[i].b0 = z264(&tmp2);
		zNroot(&tmp2,&tmp3,2);
		mult_save[i].imax = (uint64)tmp3.val[0] / 2;

		//set up recurrence
		mult_save[i].Q0 = 1;
		mult_save[i].P = mult_save[i].b0;

		zSqr(&tmp2,&tmp2);
		zSub(&tmp1,&tmp2,&tmp3);

		mult_save[i].Qn = (uint64)tmp3.val[0];
		
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
			zShortMul(N,(uint64)multipliers[j],&tmp1);

			//try to factor
			shanks_mult_unit_big(&tmp1,&mult_save[j],&f64);

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

	zFree(&tmp1);
	zFree(&tmp2);
	zFree(&tmp3);

	/*
	gettimeofday(&myTVend, NULL);
	difference = my_difftime (&myTVstart, &myTVend);

	t_time = 
		((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (1)
		printf("Total squfof time = %6.6f seconds.\n",t_time);
		*/

	return f64;
}

void shanks_mult_unit_big(z *N, mult_big_t *mult_save, uint64 *f)
{
	//use shanks SQUFOF to factor N.  almost all computation can be done with longs
	//input N < 2^63
	//return 1 in f if no factor is found
	uint64 imax,i,Q0,b0,Qn,bn,P,bbn,Ro,S,So,t1,t2;
	int j=0;
	z tmp1, tmp2;
	//FILE *out;

	zInit(&tmp1);
	zInit(&tmp2);
	
	//initialize output
	*f=0;

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
		while (1)
		{
			//at the start of every iteration, we need to know:
			//	P from the previous iteration
			//	bn from the previous iteration
			//	Qn from the previous iteration
			//	Q0 from the previous iteration
			//	iteration count, i

			t1 = P;		//hold Pn for this iteration
			P = bn*Qn - P;
			t2 = Qn;	//hold Qn for this iteration
			Qn = Q0 + bn * (t1 - P);
			Q0 = t2;	//remember last Q
			bn = (b0+P)/Qn;

			if (!(i & 0x1))  //i even
			{
				//check for square Qn = S*S
				t2 = Qn & 31;
				if (t2 == 0 || t2 == 1 || t2 == 4 ||
					t2 == 9 || t2 == 16 || t2 == 17 || t2 == 25)
				{
					//sp642z(Qn,&tmp1);
					//zNroot(&tmp1,&tmp2,2);
					//t1 = (uint64)tmp2.val[0];
					//with a max of 100 bits, Qn should always
					//be less than 53 bits, so we can use doubles
					t1 = (uint64)sqrt((double)Qn);
					if (Qn == t1 * t1)
					{
						i++;
						break;
					}
				}
			}
			i++;

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
				zFree(&tmp1);
				zFree(&tmp2);
				return;
			}
		}

		//reduce to G0
		//sp642z(Qn,&tmp1);
		//zNroot(&tmp1,&tmp2,2);
		//S = (uint64)tmp2.val[0];
		//printf("S = %lu\n",S);
		S = (int64)sqrt((double)Qn);
		Ro = P + S * ((b0 - P)/S);
		t1 = Ro;
		sp642z(t1,&tmp1);
		zSqr(&tmp1,&tmp2);
		zSub(N,&tmp2,&tmp1);
		zShortDiv(&tmp1,(fp_digit)S,&tmp2);
		//So = (uint32)(((int64)N - (int64)t1*(int64)t1)/(int64)S);
		
		So = z264(&tmp2);
		if (tmp2.size < 0)
			So = (uint64)((int64)So * -1);

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


			j++;
			
			if (j > 100000000)
			{
				//out = fopen("squfof_error.log","a");
				printf("error \n");
				//fclose(out);
				*f=-1;
				//this gets stuck very rarely, but it does happen.
				//we don't need to remember any state info, as this will
				//invalidate this multiplier
				zFree(&tmp1);
				zFree(&tmp2);
				return;		
			}

		}
	
		//*f = gcd64(Ro,N);
		sp642z(Ro,&tmp1);
		zBinGCD(&tmp1,N,&tmp2);
		if (tmp2.size == 1)
			*f = tmp2.val[0];
		else
			*f = 1;

		//found a factor - don't know yet if it's trivial or not.
		//we don't need to remember any state info, as this will
		//invalidate this multiplier
		if (*f > 1)
		{
			zFree(&tmp1);
			zFree(&tmp2);
			return;
		}
	}
}

