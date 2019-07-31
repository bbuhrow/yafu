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
#include "arith.h"

/********************* arbitrary-precision Montgomery arith **********************/
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

monty_t * monty_alloc()
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
	mpz_init(mdata->x);
	mpz_init(mdata->y);
	mpz_init(mdata->c);
	mpz_init(mdata->q);
	mpz_init(mdata->g);
	mpz_init(mdata->ys);
	mpz_init(mdata->t1);

	return mdata;
}

void monty_init(mpz_t n, monty_t *mdata)
{
	//for a input modulus n, initialize constants for 
	//montogomery representation
	//this assumes that n is relatively prime to 2, i.e. is odd.	

	// initialize space such that we can use vectorized
    mpz_set(mdata->n, n);

    fp_montgomery_calc_normalization(mdata->r, mdata->rhat, mdata->n);
    fp_montgomery_setup(mdata->n, mdata->r, mdata->nhat);

    mpz_sub_ui(mdata->r, mdata->r, 1);

	return;
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
	mpz_clear(mdata->x);
	mpz_clear(mdata->y);
	mpz_clear(mdata->q);
	mpz_clear(mdata->c);
	mpz_clear(mdata->g);
	mpz_clear(mdata->ys);
	mpz_clear(mdata->t1);

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

/********************* 128-bit Montgomery arith **********************/
void to_monty128(monty128_t *mdata, uint64 * x)
{
	//given a number x in normal (hexadecimal) representation, 
	//find its montgomery representation

	//this uses some precomputed monty constants
	//xhat = (x * r) mod n
	// = x * R^2 / R mod n
	// = REDC(x * R^2)
	mpz_t m;
	mpz_t n;

	mpz_init(m);
	mpz_init(n);

	mpz_set_ui(m, x[1]);
	mpz_mul_2exp(m, m, 64);
	mpz_add_ui(m, m, x[0]);
	
	mpz_set_ui(n, mdata->n[1]);
	mpz_mul_2exp(n, n, 64);
	mpz_add_ui(n, n, mdata->n[0]);

	// implied R = 2^128
	mpz_mul_2exp(m, m, 128);
	mpz_mod(m, m, n);

	x[0] = mpz_get_ui(m);
	mpz_tdiv_q_2exp(m, m, 64);
	x[1] = mpz_get_ui(m);

	mpz_clear(m);
	mpz_clear(n);

	return;
}

void monty128_init(monty128_t * mdata, uint64 * n)
{
	//for a input modulus n, initialize constants for 
	//montogomery representation
	//this assumes that n is relatively prime to 2, i.e. is odd.	
	uint64 b = n[0];
	uint64 x;

	mdata->n[0] = n[0];
	mdata->n[1] = n[1];

	// invert (odd) n mod 2^64
	x = (((b + 2) & 4) << 1) + b; // here x*a==1 mod 2**4
	x *= 2 - b * x;               // here x*a==1 mod 2**8
	x *= 2 - b * x;               // here x*a==1 mod 2**16
	x *= 2 - b * x;               // here x*a==1 mod 2**32         
	x *= 2 - b * x;               // here x*a==1 mod 2**64

	mdata->rho = (uint64)((uint64)0 - ((uint64)x));

	mdata->one[0] = 1;
	mdata->one[1] = 0;
	to_monty128(mdata, mdata->one);
	
	return;
}

void ciosFullMul128x(uint64 *u, uint64 *v, uint64 rho, uint64 *n, uint64 *w)
{
	__asm__(
		"movq %0, %%r10	\n\t"			/* u ptr in r10 */
		"movq %2, %%r11	\n\t"			/* w ptr in r11 */
		"movq 0(%1), %%r9	\n\t"		/* v[0] ptr in r9 */

		/* begin s += u * v */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
		"mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(u[0] * v)         */
		"addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
		"mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(u[1] * v)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
		"addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
		"adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */

		"movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
		"mulx %4, %%r9, %%r14	\n\t"   /* m = rho * w[0]         */
		"movq %3, %%r10	\n\t"			/* n ptr in r10 */

		/* begin s = (s + n * m) >> 64 */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
		"mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(n[0] * m)         */
		"addq 0(%%r11), %%r12  \n\t"    /* r12 = w[0] (could be rdx) + LO(n[0] * m) */

		/* r12 should be 0 here */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
		"mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(n[1] * m)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
		"xorq %%r14, %%r14 \n\t"
		"addq 8(%%r11), %%r12  \n\t"
		"movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

		"adcq 16(%%r11), %%r13 \n\t"
		"movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
		"adcq $0, %%r14 \n\t"
		"movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */

		/* round 2 */

		"movq %0, %%r10	\n\t"			/* u ptr in r10 */
		"movq 8(%1), %%r9	\n\t"		/* v[1] ptr in r9 */

		/* begin s += u * v */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
		"mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(u[0] * v)         */
		"addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
		"mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(u[1] * v)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
		"addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
		"adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */

		"movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
		"mulx %4, %4, %%r14	\n\t"       /* m = rho * w[0]         */
		"movq %3, %%r10	\n\t"			/* n ptr in r10 */

		/* begin s = (s + n * m) >> 64 */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
		"mulx %4, %%r12, %%r14	\n\t"   /* r14 = HI(n[0] * m)         */
		"addq 0(%%r11), %%r12  \n\t"    /* r12 = w[0] (could be rdx) + LO(n[0] * m) */

		/* r12 should be 0 here */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
		"mulx %4, %%r12, %%r13	\n\t"   /* r13 = HI(n[1] * m)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
		"xorq %%r14, %%r14 \n\t"
		"addq 8(%%r11), %%r12  \n\t"
		"movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

		"adcq 16(%%r11), %%r13 \n\t"
		"movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
		"adcq $0, %%r14 \n\t"
		"movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */



		:
		: "r"(u), "r"(v), "r"(w), "r"(n), "r"(rho)
		: "r9", "r10", "rdx", "r11", "r12", "r13", "r14", "cc", "memory");

	return;
}


void ciosMul128x(uint64 *u, uint64 v, uint64 rho, uint64 *n, uint64 *w)
{
	__asm__(
		"movq %0, %%r10	\n\t"			/* u ptr in r10 */
		"movq %2, %%r11	\n\t"			/* w ptr in r11 */

		/* begin s += u * v */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
		"mulx %1, %%r12, %%r14	\n\t"   /* r14 = HI(u[0] * v)         */
		"addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
		"mulx %1, %%r12, %%r13	\n\t"   /* r13 = HI(u[1] * v)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
		"addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
		"adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */

		"movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
		"mulx %4, %4, %%r14	\n\t"       /* m = rho * w[0]         */
		"movq %3, %%r10	\n\t"			/* n ptr in r10 */

		/* begin s = (s + n * m) >> 64 */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
		"mulx %4, %%r12, %%r14	\n\t"   /* r14 = HI(n[0] * m)         */
		"addq 0(%%r11), %%r12  \n\t"    /* we just need the potential carry bit */

		/* r12 should be 0 here */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
		"mulx %4, %%r12, %%r13	\n\t"   /* r13 = HI(n[1] * m)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
		"xorq %%r14, %%r14 \n\t"
		"addq 8(%%r11), %%r12  \n\t"    
		"movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

		"adcq 16(%%r11), %%r13 \n\t"    
		"movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
		"adcq $0, %%r14 \n\t"
		"movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */

		:
		: "r"(u), "r"(v), "r"(w), "r"(n), "r"(rho)
		: "r10", "rdx", "r11", "r12", "r13", "r14", "cc", "memory");

	return;
}

void ciosMul128x_1(uint64 *u, uint64 v, uint64 rho, uint64 *n, uint64 *w)
{
	__asm__(
		"movq %0, %%r10	\n\t"			/* u ptr in r10 */
		"movq %2, %%r11	\n\t"			/* w ptr in r11 */

		/* begin s += u * v */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
		"mulx %1, %%r12, %%r14	\n\t"   /* r14 = HI(u[0] * v)         */
		"addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
		"mulx %1, %%r12, %%r13	\n\t"   /* r13 = HI(u[1] * v)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
		"addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
		"adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */
		:
	: "r"(u), "r"(v), "r"(w), "r"(n), "r"(rho)
		: "r10", "rdx", "r11", "r12", "r13", "r14", "cc", "memory");

	return;
}

void ciosMul128x_2(uint64 *u, uint64 v, uint64 rho, uint64 *n, uint64 *w)
{
	__asm__(
		"movq %2, %%r11	\n\t"			/* w ptr in r11 */
		"movq %3, %%r10	\n\t"			/* n ptr in r10 */
		"movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
		"mulx %4, %4, %%r14	\n\t"       /* m = rho * w[0]         */

		/* begin s = (s + n * m) >> 64 */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
		"mulx %4, %%r12, %%r14	\n\t"   /* r14 = HI(n[0] * m)         */
		"addq 0(%%r11), %%r12  \n\t"    /* r12 = w[0] (could be rdx) + LO(n[0] * m) */

		/* r12 should be 0 here */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
		"mulx %4, %%r12, %%r13	\n\t"   /* r13 = HI(n[1] * m)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
		"xorq %%r14, %%r14 \n\t"
		"addq 8(%%r11), %%r12  \n\t"
		"movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

		"adcq 16(%%r11), %%r13 \n\t"
		"movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
		"adcq $0, %%r14 \n\t"
		"movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */

		:
	: "r"(u), "r"(v), "r"(w), "r"(n), "r"(rho)
		: "r10", "rdx", "r11", "r12", "r13", "r14", "cc", "memory");

	return;
}

void mulmod128(uint64 * u, uint64 * v, uint64 * w, monty128_t *mdata)
{
	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64 s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

	//printf("u: %016lx%016lx\n", u[1], u[0]);
	//printf("v: %016lx%016lx\n", v[1], v[0]);
	//printf("n: %016lx%016lx\n", mdata->n[1], mdata->n[0]);

	//ciosMul128x(u, v[0], mdata->rho, mdata->n, s);
	//ciosMul128x(u, v[1], mdata->rho, mdata->n, s);
	ciosFullMul128x(u, v, mdata->rho, mdata->n, s);

	if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
	{
		__asm__(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
			: "r"(s[0]), "r"(s[1]), "r"(mdata->n[0]), "r"(mdata->n[1]), "r"(w)
			: "r11", "cc", "memory");
	}
	else
	{
		w[0] = s[0];
		w[1] = s[1];
	}

	//printf("final result: %016lx%016lx\n", w[1], w[0]);

	return;
}

void sqrmod128(uint64 * u, uint64 * w, monty128_t *mdata)
{
	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64 s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

	//ciosMul128x(u, u[0], mdata->rho, mdata->n, s);
	//ciosMul128x(u, u[1], mdata->rho, mdata->n, s);
	ciosFullMul128x(u, u, mdata->rho, mdata->n, s);

	if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
	{
		__asm__(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
			: "r"(s[0]), "r"(s[1]), "r"(mdata->n[0]), "r"(mdata->n[1]), "r"(w)
			: "r11", "cc", "memory");
	}
	else
	{
		w[0] = s[0];
		w[1] = s[1];
	}

	return;
}

void addmod128(uint64 * a, uint64 * b, uint64 * w, uint64 * n)
{
	w[1] = a[1];
	w[0] = a[0];
	__asm__(
		"movq %0, %%r8 \n\t"
		"movq %1, %%r9 \n\t"
		"subq %4, %%r8 \n\t"		/* t = x - n */
		"sbbq %5, %%r9 \n\t"
		"addq %2, %0 \n\t"			/* x += y */
		"adcq %3, %1 \n\t"
		"addq %2, %%r8 \n\t"		/* t += y */
		"adcq %3, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+r"(w[0]), "+r"(w[1])
		: "r"(b[0]), "r"(b[1]), "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

	return;
}

void submod128(uint64 * a, uint64 * b, uint64 * w, uint64 * n)
{
	__asm__(
		"movq %6, %%r11 \n\t"
		"xorq %%r8, %%r8 \n\t"
		"xorq %%r9, %%r9 \n\t"
		"movq %0, 0(%%r11) \n\t"
		"movq %1, 8(%%r11) \n\t"
		"subq %2, 0(%%r11) \n\t"
		"sbbq %3, 8(%%r11) \n\t"
		"cmovc %4, %%r8 \n\t"
		"cmovc %5, %%r9 \n\t"
		"addq %%r8, 0(%%r11) \n\t"
		"adcq %%r9, 8(%%r11) \n\t"
		"1: \n\t"
		:
	: "r"(a[0]), "r"(a[1]), "r"(b[0]), "r"(b[1]), "r"(n[0]), "r"(n[1]), "r"(w)
		: "r8", "r9", "r11", "cc", "memory");

	return;
}
