// MIT License
// 
// Copyright (c) 2025 Ben Buhrow
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
#include "monty.h"
#include "arith.h"
#include "gmp.h"
#include <stdint.h>
#include <immintrin.h>
#include <stdio.h>

#ifdef __GNUC__
#include <stdbool.h>
#endif


#if (defined(GCC_ASM64X) || defined(__MINGW64__)) && !defined(ASM_ARITH_DEBUG)
#if !defined(_MSC_VER)
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
#endif

#if defined(_MSC_VER) && defined(__clang__)

uint64_t _udiv128(uint64_t hi, uint64_t lo, uint64_t d, uint64_t* r) // uint64_t c, uint64_t n)
{
	__asm__("divq %4"
		: "=a"(lo), "=d"(hi)
		: "1"(hi), "0"(lo), "r"(d));

	*r = hi;
	return lo;
}

#endif

/********************* arbitrary-precision Montgomery arith **********************/
// generic Montgomery arithmetic using gmp
// specialized things like mpz_powm will be much faster

void fp_montgomery_calc_normalization(mpz_t r, mpz_t rhat, mpz_t b);
int fp_montgomery_setup(mpz_t n, mpz_t r, mpz_t nhat);


void fp_montgomery_calc_normalization(mpz_t r, mpz_t r2modn, mpz_t n)
{
    int nfullwords = mpz_sizeinbase(n, 2) / GMP_LIMB_BITS + 
        ((mpz_sizeinbase(n, 2) % GMP_LIMB_BITS) > 0);

    mpz_set_ui(r, 1);
    mpz_mul_2exp(r, r, nfullwords * GMP_LIMB_BITS);

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
	// given a number x in normal representation, 
	// find its montgomery representation

    mpz_mul_2exp(x, x, mpz_sizeinbase(mdata->r,2));
    mpz_tdiv_r(x, x, mdata->n);

	return;
}

monty_t * monty_alloc()
{
	// for a input modulus n, initialize constants for 
	// montogomery representation
	// this assumes that n is relatively prime to 2, i.e. is odd.	
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
	// for a input modulus n, initialize constants for 
	// montogomery representation
	// this assumes that n is relatively prime to 2, i.e. is odd.	

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
	// given a number x in normal representation, 
	// find its montgomery representation
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
	// for a input modulus n, initialize constants for 
	// montogomery representation
	// this assumes that n is relatively prime to 2, i.e. is odd.	
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

	ASM_G(
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
#ifdef USE_PERIG_128BIT


uint64_t my_ctz64(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#if defined(__BMI1__)
	uint64_t t;
	asm(" tzcntq %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_ctzll(n);
	return 64;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_ctzll(n);
	return 64;
#else
	if (n == 0)
		return 64;
	uint64_t r = 0;
	if ((n & 0xFFFFFFFFull) == 0)
		r += 32, n >>= 32;
	if ((n & 0xFFFFull) == 0)
		r += 16, n >>= 16;
	if ((n & 0xFFull) == 0)
		r += 8, n >>= 8;
	if ((n & 0xFull) == 0)
		r += 4, n >>= 4;
	if ((n & 0x3ull) == 0)
		r += 2, n >>= 2;
	if ((n & 0x1ull) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count trailing zeroes in binary representation 
uint64_t my_ctz52(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#if defined(__BMI1__)
	uint64_t t;
	asm(" tzcntq %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_ctzll(n);
	return 52;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_ctzll(n);
	return 52;
#else
	if (n == 0)
		return 52;
	uint64_t r = 0;
	if ((n & 0xFFFFFFFFull) == 0)
		r += 32, n >>= 32;
	if ((n & 0xFFFFull) == 0)
		r += 16, n >>= 16;
	if ((n & 0xFFull) == 0)
		r += 8, n >>= 8;
	if ((n & 0xFull) == 0)
		r += 4, n >>= 4;
	if ((n & 0x3ull) == 0)
		r += 2, n >>= 2;
	if ((n & 0x1ull) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count trailing zeroes in binary representation 
uint64_t my_ctz104(uint64_t n_lo, uint64_t n_hi)
{
	if (n_lo) {
		return my_ctz52(n_lo);
	}
	return 52 + my_ctz52(n_hi);
}

uint64_t my_ctz128(uint64_t n_lo, uint64_t n_hi)
{
	if (n_lo) {
		return my_ctz64(n_lo);
	}
	return 64 + my_ctz64(n_hi);
}

// count leading zeroes in binary representation
uint64_t my_clz52(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI1__
	uint64_t t;
	asm(" lzcntq %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_clzll(n);
	return 52;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_clzll(n);
	return 52;
#else
	if (n == 0)
		return 52;
	uint64_t r = 0;
	if ((n & (0xFFFFFFFFull << 32)) == 0)
		r += 32, n <<= 32;
	if ((n & (0xFFFFull << 48)) == 0)
		r += 16, n <<= 16;
	if ((n & (0xFFull << 56)) == 0)
		r += 8, n <<= 8;
	if ((n & (0xFull << 60)) == 0)
		r += 4, n <<= 4;
	if ((n & (0x3ull << 62)) == 0)
		r += 2, n <<= 2;
	if ((n & (0x1ull << 63)) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count leading zeroes in binary representation 
uint64_t my_clz104(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return my_clz52(n_hi);
	}
	return 52 + my_clz52(n_lo);
}

// count leading zeroes in binary representation
uint64_t my_clz64(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI1__
	uint64_t t;
	asm(" lzcntq %1, %0\n": "=r"(t) : "r"(n) : "flags");
	return t;
#else
	if (n)
		return __builtin_clzll(n);
	return 64;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_clzll(n);
	return 64;
#else
	if (n == 0)
		return 64;
	uint64_t r = 0;
	if ((n & (0xFFFFFFFFull << 32)) == 0)
		r += 32, n <<= 32;
	if ((n & (0xFFFFull << 48)) == 0)
		r += 16, n <<= 16;
	if ((n & (0xFFull << 56)) == 0)
		r += 8, n <<= 8;
	if ((n & (0xFull << 60)) == 0)
		r += 4, n <<= 4;
	if ((n & (0x3ull << 62)) == 0)
		r += 2, n <<= 2;
	if ((n & (0x1ull << 63)) == 0)
		r += 1;
	return r;
#endif
#endif
}

uint64_t my_clz128(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return my_clz64(n_hi);
	}
	return 64 + my_clz64(n_lo);
}

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

void ciosModMul128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
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
void ciosModSqr128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
	uint64_t mmagic)
{
#if 1 //def _MSC_VER
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
		ASM_G(
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

void mulmod128n(uint64_t* u, uint64_t* v, uint64_t* w, uint64_t* n, uint64_t rho)
{
#ifdef USE_PERIG_128BIT
	uint64_t t[2];
	t[0] = v[0];
	t[1] = v[1];
	ciosModMul128(&t[0], &t[1], u[0], u[1], n[0], n[1], rho);
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
	ciosFullMul128x(u, v, rho, n, s);

	if ((s[2]) || (s[1] > n[1]) || ((s[1] == n[1]) && (s[0] > n[0])))
	{
		ASM_G(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
		: "r"(s[0]), "r"(s[1]), "r"(n[0]), "r"(n[1]), "r"(w)
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
	U = t[0] * rho;

	//printf("t[0] = %" PRIx64 ", t[1] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], t[1], U);

	//j = 0;
	//z = (z + x0 * y0 + u * m0) / b;
	c1 = _addcarry_u64(0, _umul128(U, n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
	t[2] = _addcarry_u64(c1, t[1], c2, &t[1]);                              // c1 applies to t2

	//j = 1;
	//z = (z + x0 * y1 + u * m1) / b;
	c1 = _addcarry_u64(0, _umul128(u[0], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
	c4 = _addcarry_u64(0, _umul128(U, n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
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
	U = t[0] * rho;

	//printf("t[0] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], U);

	//j = 0;
	//z = (z + x1 * y0 + u * m0) / b;
	c1 = _addcarry_u64(0, _umul128(U, n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
	c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
	t[2] += c1;

	//j = 1;
	//z = (z + x1 * y1 + u * m1) / b;
	c1 = _addcarry_u64(0, _umul128(u[1], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
	c4 = _addcarry_u64(0, _umul128(U, n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
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
		ASM_G(
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

void sqrmod128n(uint64_t* u, uint64_t* w, uint64_t* n, uint64_t rho)
{
#ifdef USE_PERIG_128BIT
	ciosModSqr128(&w[0], &w[1], u[0], u[1], n[0], n[1], rho);
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
	ciosFullMul128x(u, u, rho, n, s);

	if ((s[2]) || (s[1] > n[1]) || ((s[1] == n[1]) && (s[0] > n[0])))
	{
		ASM_G(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
		: "r"(s[0]), "r"(s[1]), "r"(n[0]), "r"(n[1]), "r"(w)
			: "r11", "cc", "memory");
	}
	else
	{
		w[0] = s[0];
		w[1] = s[1];
	}

#else

	mulmod128n(u, u, w, n, rho);

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
	ASM_G(
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
	ASM_G(
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

void dblmod128(uint64_t* a, uint64_t* n)
{

#if defined(__x86_64__) && defined(__unix__)
	// requires GCC_ASM64 syntax
	uint64_t t1, t0;

	// do the 2-word variant of this:
	//uint64_t r;
	//uint64_t tmp = x - n;
	//uint8_t c = _addcarry_u64(0, tmp, y, &r);
	//return (c) ? r : x + y;
	t1 = a[1];
	t0 = a[0];

	__asm__ volatile (
		"movq %%rax, %%r8 \n\t"
		"movq %%rdx, %%r9 \n\t"
		"subq %4, %%r8 \n\t"		/* t = x - n */
		"sbbq %5, %%r9 \n\t"
		"addq %%rax, %0 \n\t"		/* x += x */
		"adcq %%rdx, %1 \n\t"
		"addq %%rax, %%r8 \n\t"		/* t = t + x */
		"adcq %%rdx, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+&r"(a[0]), "+&r"(a[1])
		: "a"(t0), "d"(t1), "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

#else

#ifdef USE_PERIG_128BIT

	uint128_t x = ((uint128_t)a[1] << 64) | a[0];
	uint128_t n128 = ((uint128_t)n[1] << 64) | n[0];
	uint128_t r = (x >= n128 - x) ? x - (n128 - x) : x + x;
	a[1] = (uint64_t)(r >> 64);
	a[0] = (uint64_t)r;

#else

	addmod128(a, a, a, n);

#endif

#endif
	return;
}

void chkmod128(uint64_t* a, uint64_t* n)
{

#if defined(__x86_64__) && defined(__unix__)

	// not really faster than the branching version below
	// but probably should be checked on more cpus.
	uint64_t t0, t1;
	t0 = a[0];
	t1 = a[1];
	__asm__ volatile (
		"xorq %%r8, %%r8 \n\t"
		"xorq %%r9, %%r9 \n\t"
		"subq %2, %%r8 \n\t"		/* t = 0 - n */
		"sbbq %3, %%r9 \n\t"
		"addq %0, %%r8 \n\t"		/* t += x */
		"adcq %1, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+&r"(t0), "+&r"(t1)
		: "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

	a[0] = t0;
	a[1] = t1;
#else

	if ((a[1] > n[1]) || ((a[1] == n[1]) && (a[0] >= n[0])))
	{
		uint8_t c1 = _subborrow_u64(0, a[0], n[0], &a[0]);
		_subborrow_u64(c1, a[1], n[1], &a[1]);
	}

#endif
	return;
}

/********************* 52-bit Vector Montgomery arith **********************/
// speed critical functions are declared as static inline in the header
#ifdef USE_AVX512F


void printvec(char* msg, __m512i v)
{
	uint64_t m[8];
	storeu64(m, v);
	int i;
	printf("%s: ", msg);
	for (i = 0; i < 8; i++)
		printf("%016lx ", m[i]);
	printf("\n");
	return;
}

void mulredc52_mask_add_vec(__m512i* c0, __mmask8 addmsk, __m512i a0, __m512i b0, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, C1;

#ifndef IFMA
	__m512d prod1_hd, prod2_hd;
	__m512d prod1_ld, prod2_ld;
	__m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	__m512i one = _mm512_set1_epi64(1);
	__mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = C1 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, t1);

	// propagate t0's carry before we throw it away.
	// it is almost always exactly 1, but not always, for example when m = 0;
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));

	__mmask8 bmsk = _mm512_cmpge_epu64_mask(t0, n0);
	t0 = _mm512_mask_sub_epi64(t0, bmsk, t0, n0);

	// conditional addmod (double result)
	t0 = _mm512_mask_slli_epi64(t0, addmsk, t0, 1);
	bmsk = _mm512_mask_cmpgt_epu64_mask(addmsk, t0, n0);
	*c0 = _mm512_mask_sub_epi64(t0, bmsk & addmsk, t0, n0);
	// _mm512_and_epi64(t0, lo52mask);

	return;
}


/********************* 104-bit Vector Montgomery arith **********************/
// speed critical functions are declared as static inline in the header

uint64_t multiplicative_inverse(uint64_t a)
{
	// compute the 64-bit inverse of a mod 2^64
	//    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
	uint64_t x0 = (3 * a) ^ 2;
	uint64_t y = 1 - a * x0;
	uint64_t x1 = x0 * (1 + y);
	y *= y;
	uint64_t x2 = x1 * (1 + y);
	y *= y;
	uint64_t x3 = x2 * (1 + y);
	y *= y;
	uint64_t x4 = x3 * (1 + y);
	return x4;
}

__m512i multiplicative_inverse104_x8(uint64_t* a)
{
	//    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
	__m512i x0, x1, x2, x3, x4, x5, y, n, i0, i1;
	__m512i three = _mm512_set1_epi64(3), two = _mm512_set1_epi64(2), one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	n = loadu64(a);
	_mm512_mullo_epi52(x0, n, three);
	x0 = _mm512_xor_epi64(x0, two);
	_mm512_mullo_epi52(y, n, x0);
	y = _mm512_sub_epi64(one, y);
	y = _mm512_and_epi64(lo52mask, y);

	x1 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x1, x0, x1);
	_mm512_mullo_epi52(y, y, y);

	x2 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x2, x1, x2);
	_mm512_mullo_epi52(y, y, y);

	x3 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x3, x2, x3);
	_mm512_mullo_epi52(y, y, y);

	x4 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x4, x3, x4);

	return x4;
}


#endif












