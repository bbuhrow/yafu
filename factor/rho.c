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

#include "arith.h"
#include "factor.h"
#include "ytools.h"
#include "gmp_xface.h"
#include "monty.h"

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
    if (fobj->LOGFLAG)
    {
        flog = fopen(fobj->flogname, "a");
        if (flog == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("could not open %s for writing\n", fobj->flogname);
            return;
        }
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
		if (is_mpz_prp(fobj->rho_obj.gmp_n, fobj->NUM_WITNESSES))
		{
            char* s = mpz_get_str(NULL, 10, fobj->rho_obj.gmp_n);
			logprint(flog,"prp%d = %s\n", gmp_base10(fobj->rho_obj.gmp_n), s);

			add_to_factor_list(fobj->factors, fobj->rho_obj.gmp_n, 
                fobj->VFLAG, fobj->NUM_WITNESSES);
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;

			mpz_set_ui(fobj->rho_obj.gmp_n, 1);
            free(s);
			break;
		}

		//verbose: print status to screen
		if (fobj->VFLAG >= 0)
			printf("rho: x^2 + %u, starting %d iterations on C%u ",
			fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly], fobj->rho_obj.iterations, 
			(int)gmp_base10(fobj->rho_obj.gmp_n));

		logprint(flog, "rho: x^2 + %u, starting %d iterations on C%u\n",
			fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly], fobj->rho_obj.iterations, 
			(int)gmp_base10(fobj->rho_obj.gmp_n));
		
		//call brent's rho algorithm
		mbrent(fobj);

        if (fobj->VFLAG >= 0)
        {
            printf("\n");
        }

		//check to see if 'f' is non-trivial
		if ((mpz_cmp_ui(fobj->rho_obj.gmp_f, 1) > 0)
			&& (mpz_cmp(fobj->rho_obj.gmp_f, fobj->rho_obj.gmp_n) < 0))
		{				
            char* s = mpz_get_str(NULL, 10, fobj->rho_obj.gmp_f);

			//non-trivial factor found
			stop = clock();
			tt = (double)(stop - start)/(double)CLOCKS_PER_SEC;

			//check if the factor is prime
			if (is_mpz_prp(fobj->rho_obj.gmp_f, fobj->NUM_WITNESSES))
			{
				add_to_factor_list(fobj->factors, fobj->rho_obj.gmp_f, 
                    fobj->VFLAG, fobj->NUM_WITNESSES);

				if (fobj->VFLAG > 0)
					gmp_printf("rho: found prp%d factor = %Zd\n",
					gmp_base10(fobj->rho_obj.gmp_f),fobj->rho_obj.gmp_f);

				logprint(flog,"prp%d = %s\n",
					gmp_base10(fobj->rho_obj.gmp_f), s);
			}
			else
			{
				add_to_factor_list(fobj->factors, fobj->rho_obj.gmp_f, 
                    fobj->VFLAG, fobj->NUM_WITNESSES);

				if (fobj->VFLAG > 0)
					gmp_printf("rho: found c%d factor = %Zd\n",
					gmp_base10(fobj->rho_obj.gmp_f),fobj->rho_obj.gmp_f);

				logprint(flog,"c%d = %s\n",
					gmp_base10(fobj->rho_obj.gmp_f), s);
			}
            free(s);
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
    if (fobj->LOGFLAG)
    {
        fclose(flog);
    }
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
	*/

	mpz_t x,y,q,g,ys,t1;

	uint32_t i=0,k,r,m;
	int it;
	int imax = fobj->rho_obj.iterations;

	// initialize local arbs
	mpz_init(x);
	mpz_init(y);
	mpz_init(q);
	mpz_init(g);
	mpz_init(ys);
	mpz_init(t1);

	// starting state of algorithm.  
	r = 1;
	m = 256;
	i = 0;
	it = 0;

	mpz_set_ui(q, 1);
	mpz_set_ui(y, 0);
	mpz_set_ui(g, 1); 

	do
	{
		mpz_set(x,y);
		for(i=0;i<=r;i++)
		{
			mpz_mul(t1,y,y);		//y = (y*y + c) mod n
            mpz_add_ui(t1, t1, fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly]);
			mpz_tdiv_r(y, t1, fobj->rho_obj.gmp_n);			
		}

		k=0;
		do
		{
			mpz_set(ys, y);
			for(i=1;i<=MIN(m,r-k);i++)
			{
				mpz_mul(t1,y,y); //y=(y*y + c)%n
                mpz_add_ui(t1, t1, fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly]);
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

            // abort after the specified number of gcd's
			if (it>imax)
			{
				mpz_set_ui(fobj->rho_obj.gmp_f, 0);
				goto free;
			}

		} while (k<r && (mpz_get_ui(g) == 1));
		r*=2;
	} while (mpz_get_ui(g) == 1);

	if (mpz_cmp(g,fobj->rho_obj.gmp_n) == 0)
	{
		// back track
		do
		{
			mpz_mul(t1, ys, ys); //ys = (ys*ys + c) mod n
            mpz_add_ui(t1, t1, fobj->rho_obj.polynomials[fobj->rho_obj.curr_poly]);
			mpz_tdiv_r(ys, t1, fobj->rho_obj.gmp_n); 

			mpz_sub(t1, ys, x);
			if (mpz_sgn(t1) < 0)
				mpz_add(t1, t1, fobj->rho_obj.gmp_n);
			mpz_gcd(g, t1, fobj->rho_obj.gmp_n);
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
	//if (fobj->VFLAG >= 0)
	//	printf("\n");

	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(q);
	mpz_clear(g);
	mpz_clear(ys);
	mpz_clear(t1);
	
	return it;
}

int montybrent(monty_t *mdata, mpz_t n, mpz_t f, uint32_t a, uint32_t imax)
{
    /*
    run pollard's rho algorithm on n with Brent's modification,
    returning the first factor found in f, or else 0 for failure.
    use f(x) = x^2 + c
    see, for example, bressoud's book.
    */

	mpz_ptr x = mdata->x;
	mpz_ptr y = mdata->y;
	mpz_ptr c = mdata->c;
	mpz_ptr q = mdata->q;
	mpz_ptr g = mdata->g;
	mpz_ptr ys = mdata->ys;
	mpz_ptr t1 = mdata->t1;

    uint32_t i = 0, k, r, m;
    int it;

    // starting state of algorithm.  
    r = 1;
    m = 256;
    i = 0;
    it = 0;

    mpz_set_ui(q, 1);
    mpz_set_ui(y, 0);
    mpz_set_ui(g, 1);
    mpz_set_ui(c, a);

    monty_init(n, mdata);
    to_monty(mdata, c);

    do
    {
        mpz_set(x, y);
        for (i = 0; i <= r; i++)
        {
            monty_mul(mdata, y, y, y);
            monty_add(mdata, y, c, y);
        }

        k = 0;
        do
        {
            mpz_set(ys, y);
            for (i = 1; i <= MIN(m, r - k); i++)
            {
                monty_mul(mdata, y, y, y);
                monty_add(mdata, y, c, y);
                monty_sub(mdata, y, x, t1);
                monty_mul(mdata, t1, q, q);
            }
            mpz_gcd(g, q, mdata->n);
            k += m;
            it++;

            if (it>imax)
            {
                mpz_set_ui(f, 0);
                goto free;
            }

        } while (k<r && (mpz_get_ui(g) == 1));
        r *= 2;
    } while (mpz_get_ui(g) == 1);

    if (mpz_cmp(g, n) == 0)
    {
        // back track
        do
        {
            monty_mul(mdata, ys, ys, ys);
            monty_add(mdata, ys, c, ys);
            monty_sub(mdata, ys, x, t1);
            mpz_gcd(g, t1, mdata->n);
        } while ((mpz_size(g) == 1) && (mpz_get_ui(g) == 1));

        if (mpz_cmp(g, n) == 0)
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

    return it;
}

uint64_t spbrent(uint64_t N, uint64_t c, int imax)
{


    /*
    run pollard's rho algorithm on n with Brent's modification,
    returning the first factor found in f, or else 0 for failure.
    use f(x) = x^2 + c
    see, for example, bressoud's book.
    */
    uint64_t x, y, q, g, ys, t1, f = 0, nhat;
    uint32_t i = 0, k, r, m;
    int it;
    
    // start out checking gcd fairly often
    r = 1;

    // under 48 bits, don't defer gcd quite as long
    i = _trail_zcnt64(N);
    if (i > 20)
        m = 32;
    else if (i > 16)
        m = 160;
    else if (i > 3)
        m = 256;
    else
        m = 384;

    it = 0;
    q = 1;
    g = 1;

    x = (((N + 2) & 4) << 1) + N; // here x*a==1 mod 2**4
    x *= 2 - N * x;               // here x*a==1 mod 2**8
    x *= 2 - N * x;               // here x*a==1 mod 2**16
    x *= 2 - N * x;               // here x*a==1 mod 2**32         
    x *= 2 - N * x;               // here x*a==1 mod 2**64
    nhat = (uint64_t)0 - x;

    // Montgomery representation of c
    c = u64div(c, N);
    y = c;

    do
    {
        x = y;
        for (i = 0; i <= r; i++)
        {
            y = mulredc63(y, y + c, N, nhat);
        }

        k = 0;
        do
        {
            ys = y;
            for (i = 1; i <= MIN(m, r - k); i++)
            {
                y = mulredc63(y, y + c, N, nhat);
                t1 = x > y ? y - x + N : y - x;
                q = mulredc63(q, t1, N, nhat);
            }

            g = bingcd64(N, q);
            k += m;
            it++;

            if (it>imax)
            {
                f = 0;
                goto done;
            }

        } while ((k<r) && (g == 1));
        r *= 2;
    } while (g == 1);

    if (g == N)
    {
        //back track
        do
        {
            ys = mulredc63(ys, ys + c, N, nhat);
            t1 = x > ys ? ys - x + N : ys - x;
            g = bingcd64(N, t1);
        } while (g == 1);

        if (g == N)
        {
            f = 0;
        }
        else
        {
            f = g;
        }
    }
    else
    {
        f = g;
    }

done:

    return f;
}

uint64_t spbrent64(uint64_t N, int imax)
{


	/*
	run pollard's rho algorithm on n with Brent's modification,
	returning the first factor found in f, or else 0 for failure.
	use f(x) = x^2 + c
	see, for example, bressoud's book.
	*/
	uint64_t x, y, q, g, ys, t1, f = 0, nhat;
	uint64_t c = 1;
	uint32_t i = 0, k, r, m;
	int it;

	// start out checking gcd fairly often
	r = 1;

	// under 48 bits, don't defer gcd quite as long
	i = _trail_zcnt64(N);
	if (i > 20)
		m = 32;
	else if (i > 16)
		m = 160;
	else if (i > 3)
		m = 256;
	else
		m = 384;

	it = 0;
	q = 1;
	g = 1;

	x = (((N + 2) & 4) << 1) + N; // here x*a==1 mod 2**4
	x *= 2 - N * x;               // here x*a==1 mod 2**8
	x *= 2 - N * x;               // here x*a==1 mod 2**16
	x *= 2 - N * x;               // here x*a==1 mod 2**32         
	x *= 2 - N * x;               // here x*a==1 mod 2**64
	nhat = (uint64_t)0 - x;

	// Montgomery representation of c
	c = u64div(c, N);
	y = c;

	do
	{
		x = y;
		for (i = 0; i <= r; i++)
		{
			y = mulredc(y, y + c, N, nhat);
		}

		k = 0;
		do
		{
			ys = y;
			for (i = 1; i <= MIN(m, r - k); i++)
			{
				y = mulredc(y, y + c, N, nhat);
				t1 = x > y ? y - x + N : y - x;
				q = mulredc(q, t1, N, nhat);
			}

			g = bingcd64(N, q);
			k += m;
			it++;

			if (it > imax)
			{
				f = 0;
				goto done;
			}

		} while ((k < r) && (g == 1));
		r *= 2;
	} while (g == 1);

	if (g == N)
	{
		//back track
		do
		{
			ys = mulredc(ys, ys + c, N, nhat);
			t1 = x > ys ? ys - x + N : ys - x;
			g = bingcd64(N, t1);
		} while (g == 1);

		if (g == N)
		{
			f = 0;
		}
		else
		{
			f = g;
		}
	}
	else
	{
		f = g;
	}

done:

	return f;
}
