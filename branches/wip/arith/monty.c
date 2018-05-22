/*
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
*/

/* 
I gratefully acknowledge Tom St. Denis's TomsFastMath library, on which
many of the arithmetic routines here are based
*/

#include "yafu.h"
#include "util.h"
#include "monty.h"

/*
implements routines to perform computations with 
montgomery arithmatic
*/

// local routines
void fp_montgomery_calc_normalization(mpz_t r, mpz_t rhat, mpz_t b);
int fp_montgomery_setup(mpz_t n, mpz_t r, mpz_t nhat);


void fp_montgomery_calc_normalization(mpz_t r, mpz_t r2modn, mpz_t n)
{
    int x, bits;
    int nfullwords = mpz_sizeinbase(n, 2) / GMP_LIMB_BITS + 
        ((mpz_sizeinbase(n, 2) % GMP_LIMB_BITS) > 0);

    mpz_set_ui(r, 1);
    mpz_mul_2exp(r, r, nfullwords * GMP_LIMB_BITS);

    /*
    mpz_set(r2modn, r);

    // now compute r^2 % n so that we can do fast 
    // conversions into montgomery representation.
    bits = mpz_sizeinbase(r, 2) - 1;
    while (mpz_cmp(r2modn, n) >= 0)
    {
        mpz_tdiv_q_2exp(r2modn, r2modn, 1);
        bits++;
    }

    for (x = 0; x < bits; x++) {                
        mpz_mul_2exp(r2modn, r2modn, 1);
        if (mpz_cmp(r2modn, n) >= 0)
            mpz_sub(r2modn, r2modn, n);
    }
    */

    return;
}

int fp_montgomery_setup(mpz_t n, mpz_t r, mpz_t nhat)
{
    mpz_invert(nhat, n, r);
    mpz_sub(nhat, r, nhat);

    return 0;
}

void to_monty(monty_t *mdata, mpz_t x)
{
	//given a number x in normal (hexadecimal) representation, 
	//find its montgomery representation

	//this uses some precomputed monty constants
	//xhat = (x * r) mod n
    // = x * R^2 / R mod n
    // = REDC(x * R^2)

    //mpz_mul(x, x, mdata->rhat);
    //monty_redc(mdata, x);

    mpz_mul_2exp(x, x, mpz_sizeinbase(mdata->r,2));
    mpz_tdiv_r(x, x, mdata->n);

	return;
}

monty_t * monty_init(mpz_t n)
{
	//for a input modulus n, initialize constants for 
	//montogomery representation
	//this assumes that n is relatively prime to 2, i.e. is odd.	
	monty_t *mdata;

	mdata = (monty_t *)malloc(sizeof(monty_t));

	// initialize space such that we can use vectorized
    mpz_init(mdata->n);
    mpz_init(mdata->r);
    mpz_init(mdata->rhat);
    mpz_init(mdata->nhat);
    mpz_init(mdata->tmp);

    mpz_set(mdata->n, n);

    fp_montgomery_calc_normalization(mdata->r, mdata->rhat, mdata->n);
    fp_montgomery_setup(mdata->n, mdata->r, mdata->nhat);

    mpz_sub_ui(mdata->r, mdata->r, 1);

	return mdata;
}

void monty_add(monty_t *mdata, mpz_t u, mpz_t v, mpz_t w)
{
    mpz_add(w, u, v);
	if (mpz_cmp(w, mdata->n) >= 0)
        mpz_sub(w, mdata->n, w);

	return;
}

void monty_mul(monty_t *mdata, mpz_t u, mpz_t v, mpz_t w)
{
    mpz_mul(w, u, v);
    monty_redc(mdata, w);
	return;
}

void monty_sub(monty_t *mdata, mpz_t u, mpz_t v, mpz_t w)
{
    if (mpz_cmp(u, v) >= 0)
    {
        mpz_sub(w, u, v);
    }
    else
    {
        mpz_sub(w, u, v);
        mpz_add(w, w, mdata->n);
    }
	
	return;
}

void monty_free(monty_t *mdata)
{
    mpz_clear(mdata->n);
    mpz_clear(mdata->nhat);
    mpz_clear(mdata->rhat);
    mpz_clear(mdata->r);
    mpz_clear(mdata->tmp);

	return;
}

void monty_redc(monty_t *mdata, mpz_t x)
{
    mpz_and(mdata->tmp, x, mdata->r), mdata->tmp->_mp_size = mdata->r->_mp_size;
    mpz_mul(mdata->tmp, mdata->tmp, mdata->nhat);   
    // q'=a*b*(-n^-1) mod r
    mpz_and(mdata->tmp, mdata->tmp, mdata->r), mdata->tmp->_mp_size = mdata->r->_mp_size;
    mpz_mul(mdata->tmp, mdata->tmp, mdata->n);      // q'n
    mpz_add(x, mdata->tmp, x);                      // q'n + ab
    mpz_tdiv_q_2exp(x, x, mpz_sizeinbase(mdata->r, 2));
    if (mpz_cmp(x, mdata->n) >= 0)
        mpz_sub(x, x, mdata->n);
}

