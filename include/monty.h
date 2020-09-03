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

#ifndef MONTY_H
#define MONTY_H

#include "yafu.h"
#include "arith.h"
#include "common.h"

/********************* arbitrary-precision Montgomery arith **********************/
typedef struct
{
    mpz_t rhat;
    mpz_t nhat;
    mpz_t n;
    mpz_t r;
    mpz_t tmp;

	mpz_t x;
	mpz_t y;
	mpz_t c; 
	mpz_t q; 
	mpz_t g; 
	mpz_t ys; 
	mpz_t t1;

} monty_t;

monty_t * monty_alloc();
void monty_init(mpz_t n, monty_t *mdata);
void to_monty(monty_t *mdata, mpz_t x);
void monty_add(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_mul(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_sub(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_free(monty_t *mdata);
void monty_redc(monty_t *mdata, mpz_t x);

/********************* 128-bit Montgomery arith **********************/
typedef struct
{
	uint64 r[2];
	uint64 n[2];
	uint64 np[2];
	uint64 nhat[2];
	uint64 rhat[2];
	uint64 rmask[2];
	uint64 one[2];
	uint64 mtmp1[2];
	uint64 mtmp2[2];
	uint64 mtmp3[2];
	uint64 mtmp4[2];
	uint64 rho;
} monty128_t;


void to_monty128(monty128_t *mdata, uint64 * x);
void monty128_init(monty128_t * mdata, uint64 * n);
void mulmod128(uint64 * u, uint64 * v, uint64 * w, monty128_t *mdata);
void sqrmod128(uint64 * u, uint64 * w, monty128_t *mdata);
void addmod128(uint64 * u, uint64 * v, uint64 * w, uint64 * n);
void submod128(uint64 * u, uint64 * v, uint64 * w, uint64 * n);


/********************* 64-bit Montgomery arith **********************/

#if (defined(GCC_ASM64X) || defined(__MINGW64__)) && !defined(ASM_ARITH_DEBUG)

#if defined(USE_AVX512F) || defined(USE_BMI2)
__inline uint64 mulx64(uint64 x, uint64 y, uint64* hi) {
    __asm__(
        "mulx %3, %0, %1	\n\t"
        : "=&d"(x), "=&a"(y)
        : "0"(x), "1"(y)
    );

    *hi = y;
    return x;
}
#endif
__inline uint64 mul64(uint64 x, uint64 y, uint64* hi) {
    __asm__(
        "mulq %3	\n\t"
        : "=&a"(x), "=&d"(y)
        : "0"(x), "1"(y)
        : "cc"
    );

    *hi = y;
    return x;
}



__inline uint64 submod(uint64 a, uint64 b, uint64 n)
{
    __asm__(
        "xorq %%r8, %%r8 \n\t"
        "subq %1, %0 \n\t"
        "cmovc %2, %%r8 \n\t"
        "addq %%r8, %0 \n\t"
        : "+r"(a)
        : "r"(b), "r"(n)
        : "r8", "cc");

    return a;
}

__inline uint64 addmod(uint64 x, uint64 y, uint64 n)
{
    uint64 t = x - n;
    x += y;
    __asm__("add %2, %1\n\t"
        "cmovc %1, %0\n\t"
        :"+r" (x), "+&r" (t)
        : "r" (y)
        : "cc"
        );
    return x;
}

__inline uint64 u64div(uint64 c, uint64 n)
{
    __asm__("divq %4"
        : "=a"(c), "=d"(n)
        : "1"(c), "0ULL"(0), "r"(n));

    return n;
}


#if defined(USE_AVX512F) || defined(USE_BMI2)

__inline uint64 mulredc(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulx %2, %%r10, %%r11	\n\t"
            "movq %%r10, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "xorq %%r12, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=&d"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulx %2, %%r10, %%r11	\n\t"
            "movq %3, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "mulq %%r10 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=d"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64 sqrredc(uint64 x, uint64 n, uint64 nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulx %1, %%r10, %%r11	\n\t"
            "movq %%r10, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "xorq %%r12, %%r12 \n\t"
            "mulq %2 \n\t"
            "mulq %3 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %3, %%r12 \n\t"
            "subq %3, %%rdx \n\t"
            "cmovc %%r12, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=&d"(x)
            : "0"(x), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulx %1, %%r10, %%r11	\n\t"
            "movq %2, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "mulq %%r10 \n\t"
            "mulq %3 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "subq %3, %%rdx \n\t"
            "cmovc %3, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=d"(x)
            : "0"(x), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64 mulredc63(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    __asm__(
        "mulx %2, %%r10, %%r11	\n\t"
        "movq %3, %%rax		\n\t"
        "xorq %%r8, %%r8 \n\t"
        "mulq %%r10 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%r8 \n\t"
        "addq %%r8, %%rdx \n\t"
        : "=d"(x)
        : "0"(x), "r"(y), "r"(nhat), "r"(n)
        : "rax", "r8", "r10", "r11", "cc");

    return x;
}

__inline uint64 sqrredc63(uint64 x, uint64 n, uint64 nhat)
{
    __asm__(
        "mulx %1, %%r10, %%r11	\n\t"
        "movq %2, %%rax		\n\t"
        "xorq %%r8, %%r8 \n\t"
        "mulq %%r10 \n\t"
        "mulq %3 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "subq %3, %%rdx \n\t"
        "cmovc %3, %%r8 \n\t"
        "addq %%r8, %%rdx \n\t"
        : "=d"(x)
        : "0"(x), "r"(nhat), "r"(n)
        : "rax", "r8", "r10", "r11", "cc");

    return x;
}
#else


__inline uint64 mulredc(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "movq $0, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "xorq %%rax, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "movq $0, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64 mulredc63(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    __asm__(
        "mulq %2	\n\t"
        "movq %%rax, %%r10		\n\t"
        "movq %%rdx, %%r11		\n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(x)
        : "0"(x), "r"(y), "r"(nhat), "r"(n)
        : "rdx", "r10", "r11", "cc");

    return x;
}

__inline uint64 sqrredc(uint64 x, uint64 n, uint64 nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "movq $0, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "xorq %%rax, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(x), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "movq $0, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(x), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64 sqrredc63(uint64 x, uint64 n, uint64 nhat)
{
    __asm__(
        "mulq %2	\n\t"
        "movq %%rax, %%r10		\n\t"
        "movq %%rdx, %%r11		\n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(x)
        : "0"(x), "r"(x), "r"(nhat), "r"(n)
        : "rdx", "r10", "r11", "cc");

    return x;
}

#endif


#elif _MSC_VER

// TODO: need something portable to replace 64-bit assembler versions
// of modular multiplication.  This is getting closer, but for now these things
// are spread out locally where they are needed instead of being gathered here,
// apologies for the uglyness.

//https://msdn.microsoft.com/en-us/library/windows/desktop/82cxdw50(v=vs.85).aspx
//// ARBooker's 63-bit C mulredc
//typedef __uint128_t u128;
//
//__inline uint64 mulredc63(x, y, n, nbar)
//{
//    union { u128 d; uint64 l[2]; } t;
//    return (t.d = (u128)(x)*(y), t.d += (t.l[0] * nbar)*(u128)n,
//        t.l[1] -= n, t.l[1] + (n&((int64)t.l[1] >> 63)));
//}

#include <immintrin.h>
#include <intrin.h>

__inline uint64 u64div(uint64 c, uint64 n)
{
    uint64 r;
    //mpz_t a;
    //mpz_init(a);
    //mpz_set_ui(a, c);
    //mpz_mul_2exp(a, a, 64);
    //r = mpz_tdiv_ui(a, n);
    //mpz_clear(a);
    // first available in Visual Studio 2019
    _udiv128(c, 0, n, &r);

    return r;
}

__inline uint64 mulredc(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64(tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

__inline uint64 mulredc63(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64(tl, th, ah, &x);
    return x;
}

__inline uint64 sqrredc(uint64 x, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64(tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

__inline uint64 sqrredc63(uint64 x, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64(tl, th, ah, &x);
    return x;
}

__inline uint64 submod(uint64 a, uint64 b, uint64 n)
{
    uint64 r0;
    if (_subborrow_u64(0, a, b, &r0))
        r0 += n;
    return r0;
}

__inline uint64 addmod(uint64 x, uint64 y, uint64 n)
{
    uint64 r;
    uint8 c = _addcarry_u64(0, x, y, &r);
    
    if (c || (r < x))
        r -= n;
    return r;

    //uint64 r;
    //r = x + y;
    //if ((r >= n) || (r < x)) {
    //    r -= n;
    //}
}


// good to 60 bit inputs
__inline uint64 sqrredc60(uint64 x, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    uint8 c;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// good to 60 bit inputs
__inline uint64 mulredc60(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    uint64 th, tl, u, ah, al;
    uint8 c;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// this works if inputs are 62 bits or less
#define addmod60(x, y, n) ((x) + (y))

#else

/*
__inline uint64 u64div(uint64 c, uint64 n)
{
    uint64 q, r;
    uint64 u[2] = { c, 0 };

    // first available in Visual Studio 2019
    //_udiv128(c, 0, n, &r);

    spDivide(&q, &r, u, n);
    return r;
}

__inline uint64 mulredc(uint64 x, uint64 y, uint64 n, uint64 nhat)
{
    uint64 a_hi = x >> 32;
    uint64 b_hi = y >> 32;
    uint64 a_lo = x & 0xFFFFFFFFULL;
    uint64 b_lo = y & 0xFFFFFFFFULL;

    uint64 lo_prod = a_lo * b_lo;
    uint64 med_prod1 = a_hi * b_lo;
    uint64 med_prod2 = a_lo * b_hi;
    uint64 med_term = med_prod1 + med_prod2;
    uint64 c = 0;
    uint64 hi_prod = a_hi * b_hi;

    uint64 xy_hi = (((lo_prod >> 32) + med_term) >> 32) + hi_prod;
    if ((med_prod1 < 0 && med_prod2 < 0) || ((med_prod1 < 0 || med_prod2 < 0) && med_term >= 0)) xy_hi += 1ULL << 32;
    uint64 xy_lo = (med_term << 32) + lo_prod;

    uint64 u = xy_lo * nhat;

    a_hi = u >> 32;
    b_hi = n >> 32;
    a_lo = u & 0xFFFFFFFFULL;
    b_lo = n & 0xFFFFFFFFULL;

    lo_prod = a_lo * b_lo;
    med_prod1 = a_hi * b_lo;
    med_prod2 = a_lo * b_hi;
    med_term = med_prod1 + med_prod2;
    hi_prod = a_hi * b_hi;

    uint64 un_hi = (((lo_prod >> 32) + med_term) >> 32) + hi_prod;
    if ((med_prod1 < 0 && med_prod2 < 0) || ((med_prod1 < 0 || med_prod2 < 0) && med_term >= 0)) un_hi += 1ULL << 32;

    uint64 un_lo = ((med_term & 0xFFFFFFFFULL) << 32) + lo_prod;
    uint64 r_lo = un_lo + xy_lo;
    c = 0;
    if (r_lo < un_lo)
        c = 1;
    uint64 r_hi = un_hi + xy_hi + c;
    if (r_hi >= n)
        r_hi -= n;
    return r_hi;
}

__inline uint64 sqrredc(uint64 x, uint64 n, uint64 nhat)
{
    uint64 a_hi = x >> 32;
    uint64 a_lo = x & 0xFFFFFFFFULL;

    uint64 lo_prod = a_lo * a_lo;
    uint64 med_prod1 = a_hi * a_lo;
    uint64 med_term = med_prod1 + med_prod1;
    uint64 c = 0;
    uint64 hi_prod = a_hi * a_hi;

    // the medium term could overflow		
    uint64 xy_hi = (((lo_prod >> 32) + med_term) >> 32) + hi_prod;
    if (med_prod1 < 0) xy_hi += 1ULL << 32;
    uint64 xy_lo = (med_term << 32) + lo_prod;

    uint64 u = xy_lo * nhat;

    a_hi = u >> 32;
    uint64 b_hi = n >> 32;
    a_lo = u & 0xFFFFFFFFULL;
    uint64 b_lo = n & 0xFFFFFFFFULL;

    lo_prod = a_lo * b_lo;
    med_prod1 = a_hi * b_lo;
    uint64 med_prod2 = a_lo * b_hi;
    med_term = med_prod1 + med_prod2;
    hi_prod = a_hi * b_hi;

    // the medium term could overflow		
    uint64 un_hi = (((lo_prod >> 32) + med_term) >> 32) + hi_prod;
    if ((med_prod1 < 0 && med_prod2 < 0) || ((med_prod1 < 0 || med_prod2 < 0) && med_term >= 0)) un_hi += 1ULL << 32;

    uint64 un_lo = ((med_term & 0xFFFFFFFFULL) << 32) + lo_prod;
    uint64 r_lo = un_lo + xy_lo;
    c = 0;
    if (r_lo < un_lo)
        c = 1;
    uint64 r_hi = un_hi + xy_hi + c;
    if (r_hi >= n)
        r_hi -= n;
    return r_hi;
}

__inline uint64 submod(uint64 a, uint64 b, uint64 n)
{
    uint64 r0 = a - b;
    return a + (n & ((r0) >> 63));
}

__inline uint64 addmod(uint64 x, uint64 y, uint64 n)
{
    uint64 t = x - n;
    x += y;
    y += t;
    if (y < t) x = y;
    return x;
}

*/
#endif

#endif
