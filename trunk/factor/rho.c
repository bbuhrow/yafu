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
#include "monty.h"
#include "util.h"
#include "yafu_ecm.h"

int mbrent(fact_obj_t *fobj);

void brent_loop(fact_obj_t *fobj)
{
	//repeatedly use brent's rho on n
	//we always use three constants 'c'.
	//it may be desirable to make the number of different
	//polynomials, and their values, configurable, but for 
	//now it is hardcoded.
	z *n = &fobj->rho_obj.n;
	z *f;

	z d,t;
	int i;
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
	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->flogname);
		return;
	}

	//initialize some local args
	zInit(&d);
	zInit(&t);
	zInit(&fobj->rho_obj.factors[0]);
	f = &fobj->rho_obj.factors[0];

	fobj->rho_obj.curr_poly = 0;
	while(fobj->rho_obj.curr_poly < 3)
	{
		//for each different constant, first check primalty because each
		//time around the number may be different
		start = clock();
		if (isPrime(n))
		{
			n->type = PRP;
			logprint(flog,"prp%d = %s\n",ndigits(n),z2decstr(n,&gstr1));
			add_to_factor_list(fobj, n);
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;
			zCopy(&zOne,n);
			break;
		}

		//verbose: print status to screen
		if (VFLAG >= 0)
			printf("rho: x^2 + %u, starting %d iterations on C%d ",
			fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly], fobj->rho_obj.iterations, ndigits(n));
		logprint(flog, "rho: x^2 + %u, starting %d iterations on C%d\n",
			fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly], fobj->rho_obj.iterations, ndigits(n));
		
		//call brent's rho algorithm, using montgomery arithmetic.
		mbrent(fobj);

		//check to see if 'f' is non-trivial
		if (zCompare(f,&zOne) > 0 && zCompare(f,n) < 0)
		{	
			//non-trivial factor found
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;

			//check if the factor is prime
			if (isPrime(f))
			{
				add_to_factor_list(fobj, f);
				if (VFLAG > 0)
					printf("rho: found prp%d factor = %s\n",ndigits(f),z2decstr(f,&gstr1));
				logprint(flog,"prp%d = %s\n",
					ndigits(f),z2decstr(f,&gstr2));
			}
			else
			{
				add_to_factor_list(fobj, f);
				if (VFLAG > 0)
					printf("rho: found c%d factor = %s\n",ndigits(f),z2decstr(f,&gstr1));
				logprint(flog,"c%d = %s\n",
					ndigits(f),z2decstr(f,&gstr2));
			}
			start = clock();

			//reduce input
			zDiv(n,f,&t,&d);
			zCopy(&t,n);
		}
		else
		{
			//no factor found, log the effort we made.
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;

			fobj->rho_obj.curr_poly++; //try a different function
		}
	}

	fobj->rho_obj.ttime = tt;
	fclose(flog);
	zFree(&d);
	zFree(&t);
	return;
}


int mbrent(fact_obj_t *fobj)
{
	/*
	run pollard's rho algorithm on n with Brent's modification, 
	returning the first factor found in f, or else 0 for failure.
	use f(x) = x^2 + c
	see, for example, bressoud's book.
	use montgomery arithmetic. 
	*/

	z *n = &fobj->rho_obj.n;
	z *f;
	z x,y,q,g,ys,t1,t2,cc;

	uint32 i=0,k,r,m,c;
	int it;
	int imax = fobj->rho_obj.iterations;

	//initialize local arbs
	zInit(&x);
	zInit(&y);
	zInit(&q);
	zInit(&g);
	zInit(&ys);
	zInit(&t1);
	zInit(&t2);
	zInit(&cc);

	// make space for a factor
	fobj->rho_obj.num_factors = 1;
	f = &fobj->rho_obj.factors[0];

	//starting state of algorithm.  
	r = 1;
	m = 10;
	i = 0;
	it = 0;
	c = fobj->rho_obj.curr_poly;
	sp2z(fobj->rho_obj.polynomials[c],&cc);
	q.val[0] = 1;
	y.val[0] = 0;
	g.val[0] = 1;

	//initialize montgomery arithmetic, and a few constants
	monty_init(n);
	zCopy(&montyconst.one,&g);
	zCopy(&montyconst.one,&q);
	to_monty(&cc,n);

	do
	{
		zCopy(&y,&x);	
		for(i=0;i<=r;i++)
		{
			monty_sqr(&y,&t1,n);		//y=(y*y + c)%n
			monty_add(&t1,&cc,&y,n);
		}

		k=0;
		do
		{
			zCopy(&y,&ys);	
			for(i=1;i<=MIN(m,r-k);i++)
			{
				monty_sqr(&y,&t1,n);	//y=(y*y + c)%n
				monty_add(&t1,&cc,&y,n);	
				monty_sub(&x,&y,&t1,n);	//q = q*abs(x-y) mod n
				t1.size = abs(t1.size);
				monty_mul(&q,&t1,&t2,n);
				zCopy(&t2,&q);
			}
			zLEGCD(&q,n,&g);
			k+=m;
			it++;

			if (it>imax)
			{
				zCopy(&zZero,f);
				goto free;
			}
			g.size = abs(g.size);
		} while (k<r && ((g.size == 1) && (g.val[0] == 1)));
		r*=2;
	} while ((g.size == 1) && (g.val[0] == 1));

	if (zCompare(&g,n) == 0)
	{
		//back track
		it=0;
		do
		{
			monty_sqr(&ys,&t1,n);		//ys = (ys*ys + c) mod n
			monty_add(&t1,&cc,&ys,n);
			monty_sub(&ys,&x,&t1,n);
			zLEGCD(&t1,n,&g);
			it++;
			if (it>imax)
			{
				zCopy(&zZero,f);
				goto free;
			}
			g.size = abs(g.size);
		} while ((g.size == 1) && (g.val[0] == 1));
		if (zCompare(&g,n) == 0)
		{
			zCopy(&zZero,f);
			goto free;
		}
		else
		{
			zCopy(&g,f);
			goto free;
		}
	}
	else
	{
		zCopy(&g,f);
		goto free;
	}

free:
	if (VFLAG >= 0)
		printf("\n");

	zFree(&x);
	zFree(&y);
	zFree(&q);
	zFree(&g);
	zFree(&ys);
	zFree(&t1);
	zFree(&t2);
	zFree(&cc);
	monty_free();
	return it;
}
