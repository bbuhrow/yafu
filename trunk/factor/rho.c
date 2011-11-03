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
#include "gmp_xface.h"

int mbrent(fact_obj_t *fobj);

void brent_loop(fact_obj_t *fobj)
{
	//repeatedly use brent's rho on n
	//we always use three constants 'c'.
	//it may be desirable to make the number of different
	//polynomials, and their values, configurable, but for 
	//now it is hardcoded.
	mpz_t d,t;
	FILE *flog;
	clock_t start, stop;
	double tt;
		
	//check for trivial cases
	if ((mpz_cmp_ui(fobj->rho_obj.gmp_n, 1) == 0) || (mpz_cmp_ui(fobj->rho_obj.gmp_n, 0) == 0))
		return;

	if (mpz_cmp_ui(fobj->rho_obj.gmp_n, 2) == 0)
		return;

	//open the log file
	flog = fopen(fobj->flogname,"a");
	if (flog == NULL)
	{
		printf("could not open %s for writing\n",fobj->flogname);
		return;
	}

	//initialize some local args
	mpz_init(d);
	mpz_init(t);	

	fobj->rho_obj.curr_poly = 0;
	while(fobj->rho_obj.curr_poly < 3)
	{
		//for each different constant, first check primalty because each
		//time around the number may be different
		start = clock();
		if (mpz_probab_prime_p(fobj->rho_obj.gmp_n, NUM_WITNESSES))
		{
			logprint(flog,"prp%d = %s\n",mpz_sizeinbase(fobj->rho_obj.gmp_n,10),
				mpz_get_str(gstr1.s, 10, fobj->rho_obj.gmp_n));

			add_to_factor_list(fobj, fobj->rho_obj.gmp_n);
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;

			mpz_set_ui(fobj->rho_obj.gmp_n, 1);
			break;
		}

		//verbose: print status to screen
		if (VFLAG >= 0)
			printf("rho: x^2 + %u, starting %d iterations on C%u ",
			fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly], fobj->rho_obj.iterations, 
			(int)mpz_sizeinbase(fobj->rho_obj.gmp_n,10));

		logprint(flog, "rho: x^2 + %u, starting %d iterations on C%u\n",
			fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly], fobj->rho_obj.iterations, 
			(int)mpz_sizeinbase(fobj->rho_obj.gmp_n,10));
		
		//call brent's rho algorithm, using montgomery arithmetic.
		mbrent(fobj);

		//check to see if 'f' is non-trivial
		if ((mpz_cmp_ui(fobj->rho_obj.gmp_f, 1) > 0)
			&& (mpz_cmp(fobj->rho_obj.gmp_f, fobj->rho_obj.gmp_n) < 0))
		{				
			//non-trivial factor found
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;

			//check if the factor is prime
			if (mpz_probab_prime_p(fobj->rho_obj.gmp_f, NUM_WITNESSES))
			{
				add_to_factor_list(fobj, fobj->rho_obj.gmp_f);

				if (VFLAG > 0)
					gmp_printf("rho: found prp%d factor = %Zd\n",
					mpz_sizeinbase(fobj->rho_obj.gmp_f, 10),fobj->rho_obj.gmp_f);

				logprint(flog,"prp%d = %s\n",
					mpz_sizeinbase(fobj->rho_obj.gmp_f, 10),
					mpz_get_str(gstr1.s, 10, fobj->rho_obj.gmp_f));
			}
			else
			{
				add_to_factor_list(fobj, fobj->rho_obj.gmp_f);

				if (VFLAG > 0)
					gmp_printf("rho: found c%d factor = %Zd\n",
					mpz_sizeinbase(fobj->rho_obj.gmp_f, 10),fobj->rho_obj.gmp_f);

				logprint(flog,"c%d = %s\n",
					mpz_sizeinbase(fobj->rho_obj.gmp_f, 10),
					mpz_get_str(gstr1.s, 10, fobj->rho_obj.gmp_f));
			}
			start = clock();

			//reduce input
			mpz_tdiv_q(fobj->rho_obj.gmp_n, fobj->rho_obj.gmp_n, fobj->rho_obj.gmp_f);
			
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
	mpz_clear(d);
	mpz_clear(t);

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

	mpz_t x,y,q,g,ys,t1,t2,cc;

	uint32 i=0,k,r,m,c;
	int it;
	int imax = fobj->rho_obj.iterations;

	//initialize local arbs
	mpz_init(x);
	mpz_init(y);
	mpz_init(q);
	mpz_init(g);
	mpz_init(ys);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(cc);

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
			mpz_tdiv_r(t1, t1, fobj->rho_obj.gmp_n);			
		}

		k=0;
		do
		{
			mpz_set(ys, y);
			for(i=1;i<=MIN(m,r-k);i++)
			{
				mpz_mul(t1,y,y); //y=(y*y + c)%n
				mpz_add_ui(t1, t1, c);
				mpz_tdiv_r(y, t1, fobj->rho_obj.gmp_n);	

				mpz_sub(t1, x, y); //q = q*abs(x-y) mod n
				if (mpz_sgn(t1) < 0)
					mpz_add(t1, t1, fobj->rho_obj.gmp_n);
				mpz_mul(q, t1, q); 
				mpz_tdiv_r(q, q, fobj->rho_obj.gmp_n);	
			}
			mpz_gcd(g, q, fobj->rho_obj.gmp_n);
			k+=m;
			it++;

			if (it>imax)
			{
				mpz_set_ui(fobj->rho_obj.gmp_f, 0);
				goto free;
			}
			if (mpz_sgn(g) < 0)
				mpz_neg(g, g); 
		} while (k<r && (mpz_get_ui(g) == 1));
		r*=2;
	} while (mpz_get_ui(g) == 1);

	if (mpz_cmp(g,fobj->rho_obj.gmp_n) == 0)
	{
		//back track
		it=0;
		do
		{
			mpz_mul(t1, ys, ys); //ys = (ys*ys + c) mod n
			mpz_add_ui(t1, t1, c);
			mpz_tdiv_r(ys, t1, fobj->rho_obj.gmp_n); 

			mpz_sub(t1, ys, x);
			if (mpz_sgn(t1) < 0)
				mpz_add(t1, t1, fobj->rho_obj.gmp_n);
			mpz_gcd(g, t1, fobj->rho_obj.gmp_n);
			it++;
			if (it>imax)
			{
				mpz_set_ui(fobj->rho_obj.gmp_f, 0);
				goto free;
			}
			if (mpz_sgn(g) < 0)
				mpz_neg(g, g); 
		} while ((mpz_size(g) == 1) && (mpz_get_ui(g) == 1));
		if (mpz_cmp(g,fobj->rho_obj.gmp_n) == 0)
		{
			mpz_set_ui(fobj->rho_obj.gmp_f, 0);
			goto free;
		}
		else
		{
			mpz_set(fobj->rho_obj.gmp_f, g);
			goto free;
		}
	}
	else
	{
		mpz_set(fobj->rho_obj.gmp_f, g);
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
	
	return it;
}
