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

	//z *n = &fobj->rho_obj.n;
	//z *f;
	mpz_t n,f;
	mpz_t x,y,q,g,ys,t1,t2,cc;

	uint32 i=0,k,r,m,c;
	int it;
	int imax = fobj->rho_obj.iterations;

	mpz_init(n);

	mpz_import(n, abs(fobj->rho_obj.n.size), -1, sizeof(fp_digit), 
		0, (size_t)0, fobj->rho_obj.n.val);

	//initialize local arbs
	mpz_init(x);
	mpz_init(y);
	mpz_init(q);
	mpz_init(g);
	mpz_init(ys);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(cc);

	// make space for a factor
	fobj->rho_obj.num_factors = 1;
	zInit(&fobj->rho_obj.factors[0]);
	mpz_init(f);

	//starting state of algorithm.  
	r = 1;
	m = 10;
	i = 0;
	it = 0;
	c = fobj->rho_obj.curr_poly;
	mpz_set_ui(cc, fobj->rho_obj.polynomials[c]);
	mpz_set_ui(q, 1);
	mpz_set_ui(y, 0);
	mpz_set_ui(g, 1); 

	do
	{
		mpz_set(x,y);
		for(i=0;i<=r;i++)
		{
			mpz_mul(t1,y,y);		//y = (y*y + c) mod n
			mpz_add_ui(t1, t1, c);
			mpz_tdiv_r(t1, t1, n);			
		}

		k=0;
		do
		{
			mpz_set(ys, y);
			for(i=1;i<=MIN(m,r-k);i++)
			{
				mpz_mul(t1,y,y); //y=(y*y + c)%n
				mpz_add_ui(t1, t1, c);
				mpz_tdiv_r(y, t1, n);	

				mpz_sub(t1, x, y); //q = q*abs(x-y) mod n
				if (mpz_sgn(t1) < 0)
					mpz_add(t1, t1, n);
				mpz_mul(q, t1, q); 
				mpz_tdiv_r(q, q, n);	
			}
			mpz_gcd(g, q, n);
			k+=m;
			it++;

			if (it>imax)
			{
				mpz_set_ui(f, 0);
				goto free;
			}
			if (mpz_sgn(g) < 0)
				mpz_neg(g, g); 
		} while (k<r && (mpz_get_ui(g) == 1));
		r*=2;
	} while (mpz_get_ui(g) == 1);

	if (mpz_cmp(g,n) == 0)
	{
		//back track
		it=0;
		do
		{
			mpz_mul(t1, ys, ys); //ys = (ys*ys + c) mod n
			mpz_add_ui(t1, t1, c);
			mpz_tdiv_r(ys, t1, n); 

			mpz_sub(t1, ys, x);
			if (mpz_sgn(t1) < 0)
				mpz_add(t1, t1, n);
			mpz_gcd(g, t1, n);
			it++;
			if (it>imax)
			{
				mpz_set_ui(f, 0);
				goto free;
			}
			if (mpz_sgn(g) < 0)
				mpz_neg(g, g); 
		} while ((mpz_size(g) == 1) && (mpz_get_ui(g) == 1));
		if (mpz_cmp(g,n) == 0)
		{
			mpz_set_ui(f, 0);
			goto free;
		}
		else
		{
			mpz_set(f, g);
			goto free;
		}
	}
	else
	{
		mpz_set(f, g);
		goto free;
	}

free:
	if (VFLAG >= 0)
		printf("\n");

	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(q);
	mpz_clear(g);
	mpz_clear(ys);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(cc);
	
	mpz_export(fobj->rho_obj.factors[0].val, &c, -1, sizeof(fp_digit),
		0, (size_t)0, f);
	fobj->rho_obj.factors[0].size = c;

	mpz_clear(f);
	mpz_clear(n);

	return it;
}
