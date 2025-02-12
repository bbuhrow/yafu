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

/*

license for the ciosModMul128 code which is slightly faster than my version:


MIT License

Copyright(c) 2024 Pierre

Permission is hereby granted, free of charge, to any person obtaining a copy
of this softwareand associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright noticeand this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "common.h"
//#include "ytools.h"
#include "monty.h"
#include "arith.h"
#include "gmp.h"
#include <stdint.h>
#include <immintrin.h>

#ifdef __GNUC__
#include <stdbool.h>
#endif


#if (defined(GCC_ASM64X) || defined(__MINGW64__)) && !defined(ASM_ARITH_DEBUG)
	
__inline uint64_t _umul128(uint64_t x, uint64_t y, uint64_t* hi)
{
    __asm__(
        "mulq %3	\n\t"
        : "=&a"(x), "=&d"(y)
        : "0"(x), "1"(y)
        : "cc"
    );

    *hi = y;
    return x;
}

#endif

/********************* arbitrary-precision Montgomery arith **********************/
void fp_montgomery_calc_normalization(mpz_t r, mpz_t rhat, mpz_t b);
int fp_montgomery_setup(mpz_t n, mpz_t r, mpz_t nhat);


void fp_montgomery_calc_normalization(mpz_t r, mpz_t r2modn, mpz_t n)
{
    int nfullwords = mpz_sizeinbase(n, 2) / GMP_LIMB_BITS + 
        ((mpz_sizeinbase(n, 2) % GMP_LIMB_BITS) > 0);

    //printf("GMP_LIMB_BITS is %u\n", GMP_LIMB_BITS);
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
void to_monty128(monty128_t *mdata, uint64_t * x)
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

void monty128_init(monty128_t * mdata, uint64_t * n)
{
	//for a input modulus n, initialize constants for 
	//montogomery representation
	//this assumes that n is relatively prime to 2, i.e. is odd.	
	uint64_t b = n[0];
	uint64_t x;

	mdata->n[0] = n[0];
	mdata->n[1] = n[1];

	// invert (odd) n mod 2^64
	x = (((b + 2) & 4) << 1) + b; // here x*a==1 mod 2**4
	x *= 2 - b * x;               // here x*a==1 mod 2**8
	x *= 2 - b * x;               // here x*a==1 mod 2**16
	x *= 2 - b * x;               // here x*a==1 mod 2**32         
	x *= 2 - b * x;               // here x*a==1 mod 2**64

	mdata->rho = (uint64_t)((uint64_t)0 - ((uint64_t)x));

	mdata->one[0] = 1;
	mdata->one[1] = 0;
	to_monty128(mdata, mdata->one);
	
	return;
}

void ciosFullMul128x(uint64_t *u, uint64_t *v, uint64_t rho, uint64_t *n, uint64_t *w)
{
#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax

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

#else


#endif
	return;
}


/********************* start of Perig's 128-bit code **********************/
// Note: slightly modified modular subtract at the end
#ifdef _MSC_VER

#else
#define USE_PERIG_128BIT
#endif

#ifdef USE_PERIG_128BIT

typedef __uint128_t uint128_t;

// borrow::diff = a - b - borrow_in
#define INLINE_ASM 1
static inline uint8_t my_sbb64(uint8_t borrow_in, uint64_t a, uint64_t b, uint64_t* diff)
{
#if (INLINE_ASM && (defined( __INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)))
	return _subborrow_u64(borrow_in, a, b, (unsigned long long*)diff);
#elif defined(__GNUC__)
	bool c;
	c = __builtin_usubll_overflow(a, b, (unsigned long long*)diff);
	c |= __builtin_usubll_overflow(*diff, borrow_in, (unsigned long long*)diff);
	return c;
#elif defined(_MSC_VER)
	return _subborrow_u64(borrow_in, a, b, (unsigned long long*)diff);
#else
	if (__builtin_constant_p(borrow_in) && borrow_in == 0) {
		if (__builtin_constant_p(a) && a == 0) {
			*diff = -b;
			return 1;
		}
		else if (__builtin_constant_p(b) && b == 0) {
			*diff = a;
			return 0;
		}
		else {
			uint64_t tmp = a - b;
			uint8_t borrow = (tmp > a);
			*diff = tmp;
			return borrow;
		}
	}
	else {
		uint64_t tmp1 = a - borrow_in;
		uint8_t borrow = (tmp1 > a);
		if (__builtin_constant_p(b) && b == 0) {
			*diff = tmp1;
			return borrow;
		}
		else {
			uint64_t tmp2 = tmp1 - b;
			borrow |= (tmp2 > tmp1);
			*diff = tmp2;
			return borrow;
		}
	}
#endif
}

// subtract until the result is less than the modulus
static void ciosSubtract128(uint64_t* res_lo, uint64_t* res_hi, uint64_t carries, uint64_t mod_lo, uint64_t mod_hi)
{
	uint64_t n_lo, n_hi;
	uint64_t t_lo, t_hi;
	uint8_t b;
	n_lo = *res_lo;
	n_hi = *res_hi;
	// save, subtract the modulus until a borrows occurs
	do {
		t_lo = n_lo;
		t_hi = n_hi;
		b = my_sbb64(0, n_lo, mod_lo, &n_lo);
		b = my_sbb64(b, n_hi, mod_hi, &n_hi);
#ifdef _MSC_VER
		b = my_sbb64(b, carries, 0, &carries);
#else
		if (__builtin_constant_p(carries) && carries == 0) {
		}
		else {
			b = my_sbb64(b, carries, 0, &carries);
		}
#endif
		
	} while (b == 0);
	// get the saved values when a borrow occurs
	*res_lo = t_lo;
	*res_hi = t_hi;
}

static void ciosModMul128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
	uint64_t mmagic)
{

#ifdef _MSC_VER
	uint64_t a_lo = *res_lo, a_hi = *res_hi;
	uint64_t cshi, cslo, cchi, cclo;
	uint64_t t0, t1, t2, t3, m, ignore;

	//cc = (uint128_t)a_lo * b_lo;	// #1
	//t0 = (uint64_t)cc;
	//cc = cc >> 64;
	cclo = _umul128(a_lo, b_lo, &cchi);
	t0 = cclo;
	cclo = cchi;

	//cc += (uint128_t)a_lo * b_hi;	// #2
	cslo = _umul128(a_lo, b_hi, &cshi);
	cchi = _addcarry_u64(0, cclo, cslo, &cclo);
	cchi += cshi;
	
	//t1 = (uint64_t)cc;
	//cc = cc >> 64;
	//t2 = (uint64_t)cc;
	t1 = cclo;
	t2 = cchi;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	//cs = (uint128_t)m * mod_lo;	// #4
	cslo = _umul128(m, mod_lo, &cshi);

	//cs += t0;
	//cs = cs >> 64;
	cshi += _addcarry_u64(0, t0, cslo, &cslo);
	cslo = cshi;

	//cs += (uint128_t)m * mod_hi;	// #5
	//cs += t1;
	cclo = _umul128(m, mod_hi, &cchi);
	cchi += _addcarry_u64(0, cclo, cslo, &cclo);
	cchi += _addcarry_u64(0, t1, cclo, &cclo);

	//t0 = (uint64_t)cs;
	//cs = cs >> 64;
	t0 = cclo;
	cslo = cchi;
	
	//cs += t2;
	cshi = _addcarry_u64(0, t2, cslo, &cslo);

	//t1 = (uint64_t)cs;
	//cs = cs >> 64;
	//t2 = (uint64_t)cs;
	t1 = cslo;
	t2 = cshi;

#if PARANOID
	assert(cs >> 64 == 0);
#endif

	//cc = (uint128_t)a_hi * b_lo;	// #6
	//cc += t0;
	//t0 = (uint64_t)cc;
	//cc = cc >> 64;
	cclo = _umul128(a_hi, b_lo, &cchi);
	cchi += _addcarry_u64(0, cclo, t0, &cclo);
	t0 = cclo;
	cclo = cchi;
	
	//cc += (uint128_t)a_hi * b_hi;	// #7
	//cc += t1;
	//t1 = (uint64_t)cc;
	//cc = cc >> 64;
	cslo = _umul128(a_hi, b_hi, &cshi);
	cshi += _addcarry_u64(0, cclo, cslo, &cslo);
	cshi += _addcarry_u64(0, cslo, t1, &cslo);
	t1 = cslo;
	cclo = cshi;

	//cc += t2;
	//t2 = (uint64_t)cc;
	//cc = cc >> 64;
	//t3 = (uint64_t)cc;
	cchi = _addcarry_u64(0, cclo, t2, &cclo);
	t2 = cclo;
	t3 = cchi;

#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	//cs = (uint128_t)m * mod_lo;	// #9
	//cs += t0;
	//cs = cs >> 64;
	cslo = _umul128(m, mod_lo, &cshi);
	cshi += _addcarry_u64(0, t0, cslo, &cslo);
	cslo = cshi;

	//cs += (uint128_t)m * mod_hi;	// #10
	//cs += t1;
	cclo = _umul128(m, mod_hi, &cchi);
	cchi += _addcarry_u64(0, cclo, cslo, &cclo);
	cchi += _addcarry_u64(0, t1, cclo, &cclo);

	//t0 = (uint64_t)cs;
	//cs = cs >> 64;
	t0 = cclo;
	cslo = cchi;

	//cs += t2;
	cshi = _addcarry_u64(0, t2, cslo, &cslo);

	//t1 = (uint64_t)cs;
	//cs = cs >> 64;
	t1 = cslo;
	cslo = cshi;

	//cs += t3;
	//t2 = (uint64_t)cs;
	cshi = _addcarry_u64(0, t3, cslo, &cslo);
	t2 = cslo;

	if (t2) {
		unsigned char carry = _subborrow_u64(0, t0, mod_lo, &t0);
		_subborrow_u64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}


#else
	uint64_t a_lo = *res_lo, a_hi = *res_hi;
	uint128_t cs, cc;
	uint64_t t0, t1, t2, t3, m, ignore;

	cc = (uint128_t)a_lo * b_lo;	// #1
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)a_lo * b_hi;	// #2
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	t2 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #5
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	t2 = (uint64_t)cs;
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	cc = (uint128_t)a_hi * b_lo;	// #6
	cc += t0;
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)a_hi * b_hi;	// #7
	cc += t1;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t)cc;
	cc = cc >> 64;
	t3 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	cs = (uint128_t)m * mod_lo;	// #9
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #10
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t)cs;
	if (t2) {
		unsigned char carry = my_sbb64(0, t0, mod_lo, &t0);
		my_sbb64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}

#endif
	*res_lo = t0;
	*res_hi = t1;
}

/********************* end of Perig's 128-bit code **********************/

// modSqr version I wrote based on the modMul
static void ciosModSqr128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
	uint64_t mmagic)
{
#ifdef _MSC_VER
	*res_lo = b_lo;
	*res_hi = b_hi;
	ciosModMul128(res_lo, res_hi, b_lo, b_hi, mod_lo, mod_hi, mmagic);
#else
	uint128_t cs, cc, b_lohi;
	uint64_t t0, t1, t2, t3, m, ignore;

	cc = (uint128_t)b_lo * b_lo;	// #1
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	b_lohi = (uint128_t)b_lo * b_hi;	// #2
	cc += b_lohi;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	t2 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #5
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	t2 = (uint64_t)cs;
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	cc = b_lohi + t0;
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)b_hi * b_hi;	// #6
	cc += t1;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t)cc;
	cc = cc >> 64;
	t3 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	cs = (uint128_t)m * mod_lo;	// #9
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #10
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t)cs;
	if (t2) {
		unsigned char carry = my_sbb64(0, t0, mod_lo, &t0);
		my_sbb64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}

	*res_lo = t0;
	*res_hi = t1;

#endif

}
#endif

// already defined within mingw64/msys2
#if 0 //defined( GCC_ASM64X ) && !defined(__MINGW32__)

__inline uint8_t _addcarry_u64(uint64_t x, uint8_t w, uint64_t y, uint64_t *sum)
{
    uint64_t s, c; 
    s = y;
    c = 0;

    ASM_G("movq %2, %%rax		\n\t"
        "addq %3, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "addq %%rax, %4		\n\t"
        "adcq $0, %5		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(x), "r"((uint64_t)w), "0"(s), "1"(c)
        : "rax", "memory", "cc");

    *sum = s;
    return c;
}

#endif

void mulmod128(uint64_t * u, uint64_t * v, uint64_t * w, monty128_t *mdata)
{
#ifdef USE_PERIG_128BIT
	uint64_t t[2];
	t[0] = v[0];
	t[1] = v[1];
	ciosModMul128(&t[0], &t[1], u[0], u[1], mdata->n[0], mdata->n[1], mdata->rho);
	w[0] = t[0];
	w[1] = t[1];
	return;

#else


	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64_t s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
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
#else
    // TODO: implement portable u128 x u128 modular multiplication
    uint64_t t[3], U, c2, c3;
    uint8_t c1, c4;

    // z = 0
    // for (i = 0; i < t; i++)
    // {
    //     u = (z0 + xi * y0) * -m’ mod b
    //     z = (z + xi * y + u * m) / b
    // }
    // if (x >= m) { x -= m; }

    //printf("nhat = %" PRIx64 "\n", mdata->rho);
    //printf("a = %" PRIx64 ",%" PRIx64 "\n", u[1], u[0]);
    //printf("b = %" PRIx64 ",%" PRIx64 "\n", v[1], v[0]);
    //printf("n = %" PRIx64 ",%" PRIx64 "\n", mdata->n[1], mdata->n[0]);

    //i = 0;
    //u = (z0 + x0 * y0) * nhat mod b;
    t[0] = _umul128(u[0], v[0], &t[1]);
    U = t[0] * mdata->rho;

    //printf("t[0] = %" PRIx64 ", t[1] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], t[1], U);

    //j = 0;
    //z = (z + x0 * y0 + u * m0) / b;
    c1 = _addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    t[2] = _addcarry_u64(c1, t[1], c2, &t[1]);                              // c1 applies to t2
    
    //j = 1;
    //z = (z + x0 * y1 + u * m1) / b;
    c1 = _addcarry_u64(0, _umul128(u[0], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = _addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2
    t[0] = t[1];                                                            // divide by b
    t[1] = t[2];
    t[2] = c1 + c4;

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    //i = 0;
    //u = (z + x1 * y0) * nhat mod b;
    c1 = _addcarry_u64(0, _umul128(u[1], v[0], &c2), t[0], &t[0]);          // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;
    U = t[0] * mdata->rho;

    //printf("t[0] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], U);

    //j = 0;
    //z = (z + x1 * y0 + u * m0) / b;
    c1 = _addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;

    //j = 1;
    //z = (z + x1 * y1 + u * m1) / b;
    c1 = _addcarry_u64(0, _umul128(u[1], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = _addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    w[0] = t[1];
    w[1] = t[2];

    //exit(1);

    return;



#endif

#endif
	return;
}

void sqrmod128(uint64_t * u, uint64_t * w, monty128_t *mdata)
{
#ifdef USE_PERIG_128BIT
	ciosModSqr128(&w[0], &w[1], u[0], u[1], mdata->n[0], mdata->n[1], mdata->rho);
	return;

#else


	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64_t s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
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

#else

    mulmod128(u, u, w, mdata);

#endif

#endif
	return;
}

void addmod128(uint64_t * a, uint64_t * b, uint64_t * w, uint64_t * n)
{
#if defined(GCC_ASM64X)
    // requires GCC_ASM64 syntax
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

#else

    uint8_t c;
    uint64_t t[2];
    c = _addcarry_u64(0, a[0], b[0], &t[0]);
    c = _addcarry_u64(c, a[1], b[1], &t[1]);
    if (c || (t[1] > n[1]) || ((t[1] == n[1]) && (t[0] > n[0])))
    {
        c = _subborrow_u64(0, t[0], n[0], &w[0]);
        c = _subborrow_u64(c, t[1], n[1], &w[1]);
    }
    else
    {
        w[0] = t[0];
        w[1] = t[1];
    }


#endif
	return;
}

void submod128(uint64_t * a, uint64_t * b, uint64_t * w, uint64_t * n)
{
#if defined(GCC_ASM64X)
    // requires GCC_ASM64 syntax
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

#else

    uint8_t c;
    uint64_t t[2];
    c = _subborrow_u64(0, a[0], b[0], &t[0]);
    c = _subborrow_u64(c, a[1], b[1], &t[1]);
    if (c)
    {
        c = _addcarry_u64(0, t[0], n[0], &w[0]);
        c = _addcarry_u64(c, t[1], n[1], &w[1]);
    }
    else
    {
        w[0] = t[0];
        w[1] = t[1];
    }

#endif

	return;
}
