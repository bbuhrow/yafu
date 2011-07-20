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
#include "soe.h"

void zTrial(fact_obj_t *fobj)
{
	//trial divide n using primes below limit. optionally, print factors found
	uint32 r,k=0;
	uint32 limit = fobj->div_obj.limit;
	int print = fobj->div_obj.print;
	fp_digit q;
	z *n = &fobj->div_obj.n;
	z tmp;

	zInit(&tmp);

	if (P_MAX < limit)
		GetPRIMESRange(0,limit);

	while (!(n->size == 1 && n->val[0] <= 1) && (PRIMES[k] < limit))
	{
		if (k >= NUM_P)
		{
			GetPRIMESRange(PRIMES[k-1]-1,PRIMES[k-1]+10000000);
			k=0;
		}

		q = (fp_digit)PRIMES[k];
		r = zShortMod(n,q);
		
		if (r != 0)
			k++;
		else
		{
			zShortDiv(n,q,n);
			sp2z(q,&tmp);
			tmp.type = PRIME;
			add_to_factor_list(fobj, &tmp);
			if (print && (VFLAG > 0))
				printf("div: found prime factor = %" PRIu64 "\n",q);
		}
	}

	zFree(&tmp);
}

void zFermat(fp_digit limit, fact_obj_t *fobj)
{
 //	  Fermat's factorization method (wikipedia psuedo-code)


	z *n = &fobj->div_obj.n;
	z a, b2, tmp, tmp2, tmp3, maxa;
	int i;
	int numChars;
	fp_digit reportIt, reportInc;
	fp_digit count;

	if ((n->val[0] & 1) == 0)
	{
		zShiftRight(n,n,1);
		isPrime(&zTwo);
		add_to_factor_list(fobj, &zTwo);
		return;
	}

	zInit(&a);

	if (isSquare(n))
	{
		zNroot(n,&a,2);	
		isPrime(&a);
		add_to_factor_list(fobj, &a);
		add_to_factor_list(fobj, &a);
		zCopy(&zOne,n);
		zFree(&a);
		return;
	}

	zInit(&b2);
	zInit(&tmp);
	zInit(&maxa);

	zShortAdd(n,1,&maxa);
	zShiftRight(&maxa,&maxa,1);

	zNroot(n,&a,2);		//floor(sqrt(a))
	zShortAdd(&a,1,&a); //ceil(sqrt(a))
	zSqr(&a,&tmp);
	zSub(&tmp,n,&b2);

	count = 0;
	numChars = 0;
	reportIt = limit / 100;
	reportInc = reportIt;
	while((b2.size > 0) && (!isSquare(&b2)))
	{
		//special case the increment
		if (a.val[0] < MAX_DIGIT)
			a.val[0]++;	
		else
			zShortAdd(&a,1,&a);

		if (zCompare(&maxa,&a) > 0)
			break;	//give up

		//b2 = a*a - N = b2 + 2*a - 1
		zShiftLeft(&tmp,&a,1);
		zAdd(&b2,&tmp,&b2);

		//special case the decrement
		if (b2.val[0] > 0)
			b2.val[0]--;
		else
			zShortSub(&b2,1,&b2);

		count++;
		if (count > limit)
			break;

		//progress report
		if ((count > reportIt) && (VFLAG > 1))
		{
			for (i=0; i< numChars; i++)
				printf("\b");
			numChars = printf("%" PRIu64 "%%",(uint64)((double)count / (double)limit * 100));
			fflush(stdout);
			reportIt += reportInc;
		}
	}

	if ((b2.size > 0) && (isSquare(&b2)))
	{
		//printf("found square at count = %d: a = %s, b2 = %s",count,
		//	z2decstr(&a,&gstr1),z2decstr(&b2,&gstr2));
		zInit(&tmp2);
		zInit(&tmp3);
		zNroot(&b2,&tmp,2);
		zAdd(&a,&tmp,&tmp);
		isPrime(&tmp);	//sets type property of tmp
		add_to_factor_list(fobj, &tmp);
		zDiv(n,&tmp,&tmp2,&tmp3);
		zCopy(&tmp2,n);
		zNroot(&b2,&tmp,2);
		zSub(&a,&tmp,&tmp);
		isPrime(&tmp);
		add_to_factor_list(fobj, &tmp);
		zDiv(n,&tmp,&tmp2,&tmp3);
		zCopy(&tmp2,n);
		zFree(&tmp2);
		zFree(&tmp3);
	}

	zFree(&tmp);
	zFree(&a);
	zFree(&b2);
	zFree(&maxa);

	return;

}

