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
	uint64_t r[2];
	uint64_t n[2];
	uint64_t np[2];
	uint64_t nhat[2];
	uint64_t rhat[2];
	uint64_t rmask[2];
	uint64_t one[2];
	uint64_t mtmp1[2];
	uint64_t mtmp2[2];
	uint64_t mtmp3[2];
	uint64_t mtmp4[2];
	uint64_t rho;
} monty128_t;

#ifdef _MSC_VER

#else
#define USE_PERIG_128BIT
#endif

#ifdef USE_PERIG_128BIT

typedef __uint128_t uint128_t;

#endif

void to_monty128(monty128_t *mdata, uint64_t * x);
void monty128_init(monty128_t * mdata, uint64_t * n);
void mulmod128(uint64_t * u, uint64_t * v, uint64_t * w, monty128_t *mdata);
void sqrmod128(uint64_t * u, uint64_t * w, monty128_t *mdata);
void sqrmod128n(uint64_t* u, uint64_t* w, uint64_t *n, uint64_t* nhat);
void addmod128(uint64_t * u, uint64_t * v, uint64_t * w, uint64_t * n);
void submod128(uint64_t * u, uint64_t * v, uint64_t * w, uint64_t * n);
void dblmod128(uint64_t* a, uint64_t* n);
void chkmod128(uint64_t* a, uint64_t* n);
uint64_t my_clz64(uint64_t n);
uint64_t my_clz128(uint64_t n_lo, uint64_t n_hi);
uint64_t my_clz104(uint64_t n_lo, uint64_t n_hi);
uint64_t my_clz52(uint64_t n);
uint64_t my_ctz104(uint64_t n_lo, uint64_t n_hi);
uint64_t my_ctz52(uint64_t n);

/********************* 64-bit Montgomery arith **********************/

#if (defined(GCC_ASM64X) || defined(__MINGW64__)) && !defined(ASM_ARITH_DEBUG)


__inline uint64_t _umul128(uint64_t x, uint64_t y, uint64_t* hi);

#if defined(USE_AVX512F) || defined(USE_BMI2)
__inline uint64_t mulx64(uint64_t x, uint64_t y, uint64_t* hi) {
    __asm__(
        "mulx %3, %0, %1	\n\t"
        : "=&d"(x), "=&a"(y)
        : "0"(x), "1"(y)
    );

    *hi = y;
    return x;
}
#endif
__inline uint64_t mul64(uint64_t x, uint64_t y, uint64_t* hi) {
    __asm__(
        "mulq %3	\n\t"
        : "=&a"(x), "=&d"(y)
        : "0"(x), "1"(y)
        : "cc"
    );

    *hi = y;
    return x;
}

__inline uint64_t submod(uint64_t a, uint64_t b, uint64_t n)
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

__inline uint64_t addmod(uint64_t x, uint64_t y, uint64_t n)
{
    uint64_t t = x - n;
    x += y;
    __asm__("add %2, %1\n\t"
        "cmovc %1, %0\n\t"
        :"+r" (x), "+&r" (t)
        : "r" (y)
        : "cc"
        );
    return x;
}

__inline uint64_t u64div(uint64_t c, uint64_t n)
{
#if 1
    __asm__("divq %4"
        : "=a"(c), "=d"(n)
        : "1"(c), "0"(0ULL), "r"(n));
#else
// this should work if the above won't compile (e.g. on clang)
    uint64_t tmp = 0;
    __asm__("divq %4"
        : "=a"(tmp), "=d"(n)
        : "1"(c), "0"(tmp), "r"(n));
#endif
    return n;
}


#if defined(USE_AVX512F) || defined(USE_BMI2)

__inline uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
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

__inline uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
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

__inline uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
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

__inline uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
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


__inline uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
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

__inline uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
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

__inline uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
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

__inline uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
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

#include <immintrin.h>
#include <intrin.h>

#if defined(_MSC_VER) && defined(__clang__)
extern uint64_t _udiv128(uint64_t hi, uint64_t lo, uint64_t d, uint64_t* r);
#endif

__inline uint64_t u64div(uint64_t c, uint64_t n)
{
    uint64_t r;
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

__inline uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

__inline uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    return x;
}

__inline uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

__inline uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    return x;
}

__inline uint64_t submod(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r0;
    if (_subborrow_u64(0, a, b, &r0))
        r0 += n;
    return r0;
}

__inline uint32_t submod32(uint32_t a, uint32_t b, uint32_t n)
{
    uint32_t r0;
    if (_subborrow_u32(0, a, b, &r0))
        r0 += n;
    return r0;
}

__inline uint64_t addmod(uint64_t x, uint64_t y, uint64_t n)
{
#if 0
    uint64_t r;
    uint64_t tmp = x - n;
    uint8_t c = _addcarry_u64(0, tmp, y, &r);
    return (c) ? r : x + y;
#else
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x>=n-y) ? x-(n-y) : x+y;
#endif
}

__inline uint32_t addmod32(uint32_t x, uint32_t y, uint32_t n)
{
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x >= n - y) ? x - (n - y) : x + y;
}



// good to 60 bit inputs
__inline uint64_t sqrredc60(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    uint8_t c;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// good to 60 bit inputs
__inline uint64_t mulredc60(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    uint8_t c;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// this works if inputs are 62 bits or less
#define addmod60(x, y, n) ((x) + (y))

#endif


#ifdef USE_AVX512F

#include <immintrin.h>


#define and64 _mm512_and_epi64
#define store64 _mm512_store_epi64
#define storeu64 _mm512_storeu_epi64
#define mstoreu64 _mm512_mask_storeu_epi64
#define storeu512 _mm512_storeu_si512
#define add64 _mm512_add_epi64
#define sub64 _mm512_sub_epi64
#define set64 _mm512_set1_epi64
#define srli64 _mm512_srli_epi64
#define load64 _mm512_load_epi64
#define loadu64 _mm512_loadu_epi64
#define loadu512 _mm512_loadu_si512
#define castpd _mm512_castsi512_pd
#define castepu _mm512_castpd_si512

typedef struct
{
    //ALIGNED_MEM
    uint64_t data[2][8];
} vec_u104_t;

/********************* 52-bit vector Montgomery arith **********************/
#define carryprop(lo, hi, mask) \
	{ __m512i carry = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, carry);		\
	lo = _mm512_and_epi64(mask, lo); }

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
#define ROUNDING_MODE (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)
#else
#define ROUNDING_MODE _MM_FROUND_CUR_DIRECTION
#endif


#ifdef IFMA

#define FORCE_INLINE __inline

#define mul52hi(b, c) \
	_mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b)


#define mul52lo(b, c) \
	_mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b)

FORCE_INLINE static void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
	*l = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b);
	*h = _mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b);
	return;
}

#else

static __m512d dbias;
static __m512i vbias1;
static __m512i vbias2;
static __m512i vbias3;

#define mul52lo(b, c) \
	_mm512_and_si512(_mm512_mullo_epi64(b, c), _mm512_set1_epi64(0x000fffffffffffffull))

__inline static __m512i mul52hi(__m512i b, __m512i c)
{
	__m512d prod1_ld = _mm512_cvtepu64_pd(b);
	__m512d prod2_ld = _mm512_cvtepu64_pd(c);
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	return _mm512_sub_epi64(castepu(prod1_ld), vbias1);
}
__inline static void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
	__m512d prod1_ld = _mm512_cvtepu64_pd(b);
	__m512d prod2_ld = _mm512_cvtepu64_pd(c);
	__m512d prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	*h = _mm512_sub_epi64(castepu(prod1_hd), vbias1);
	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	*l = _mm512_castpd_si512(prod1_ld);
	*l = _mm512_and_si512(*l, _mm512_set1_epi64(0x000fffffffffffffull));
	*h = _mm512_and_si512(*h, _mm512_set1_epi64(0x000fffffffffffffull));
	return;
}

#endif



#ifdef IFMA
#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), a, b);

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
    lo = _mm512_madd52lo_epu64(lo, a, b); \
    hi = _mm512_madd52hi_epu64(hi, a, b);
#else

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi = _mm512_add_epi64(hi, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo = _mm512_add_epi64(lo, _mm512_sub_epi64(castepu(prod1_ld), vbias3));

#define VEC_MUL2_ACCUM_LOHI_PD(c, a, b, lo1, hi1, lo2, hi2) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
	prod3_ld = _mm512_cvtepu64_pd(c);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi1 = _mm512_add_epi64(hi1, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
	hi2 = _mm512_add_epi64(hi2, _mm512_sub_epi64(castepu(prod2_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod2_hd = _mm512_sub_pd(castpd(vbias2), prod2_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo1 = _mm512_add_epi64(lo1, _mm512_sub_epi64(castepu(prod1_ld), vbias3)); \
	lo2 = _mm512_add_epi64(lo2, _mm512_sub_epi64(castepu(prod2_ld), vbias3));

#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_and_si512(_mm512_mullo_epi64(a, b), _mm512_set1_epi64(0x000fffffffffffffull));
#endif


// better to #define these?  Or make them static inline here?
//__m512i _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_addsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_subsetc_epi52(__m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_subsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout);
//__m512i _mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout);

// static inline works and is quite a bit faster for tinyecm
__inline static __m512i _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_mask_addsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_subsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_mask_subsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(mask, b, a);
    t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_mask_add_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}
__inline static __m512i _mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_mask_sub_epi64(a, m, a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
    __m512i t2 = _mm512_mask_sub_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}
__inline static __m512i _mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = c & _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

void mulredc52_mask_add_vec(__m512i* c0, __mmask8 addmsk, __m512i a0, __m512i b0, __m512i n0, __m512i vrho);

/********************* 104-bit vector Montgomery arith **********************/
// can we make these static inline in the header?  does it help?
//void mask_mulredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//    __m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0, __m512i vrho);
//void sqrredc104_vec(__m512i* c1, __m512i* c0,
//    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void mask_sqrredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void mask_sqrredc104_vec_pos(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void mask_sqrredc104_exact_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void addmod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
//	__m512i b1, __m512i b0, __m512i n1, __m512i n0);
//void mask_addmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
//	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0);
//void mask_dblmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0);
//void mask_redsub104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0);
//void redsub104_x8(__m512i* c1, __m512i* c0,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0);
//void submod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
//	__m512i b1, __m512i b0, __m512i n1, __m512i n0);

__inline static void mask_mulredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, t3, C1, C2;

#ifndef IFMA
	__m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	__m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	__m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	__m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	__m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	__m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	int biascount = 0;
	__m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	__m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = t3 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b0, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b0, a0, a1, t0, t1, C1, t2);
	//t1 = _mm512_add_epi64(t1, C1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b1, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b1, a0, a1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);
	C1 = C2 = zero;

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);

	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
__inline static void mulredc104_vec(__m512i* c1, __m512i* c0, 
	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, t3, C1, C2;

#ifndef IFMA
	__m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	__m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	__m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	__m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	__m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	__m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	int biascount = 0;
	__m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	__m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = t3 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b0, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b0, a0, a1, t0, t1, C1, t2);
	//t1 = _mm512_add_epi64(t1, C1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b1, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b1, a0, a1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);
	C1 = C2 = zero;

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);

	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);

	return;
}
__inline static void sqrredc104_vec(__m512i* c1, __m512i* c0,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, t3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
	__m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	__m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	__m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	__m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	__m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	__m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	int biascount = 0;
	__m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	__m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = t3 = C1 = C2 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	t1 = sqr_lo;
	t2 = sqr_hi;
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	// t0 = _mm512_add_epi64(t1, one);
	// t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);

	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		__mmask8 bmsk;
		bmsk = _mm512_cmpgt_epi64_mask(n0, t1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, _mm512_mask_set1_epi64(zero, bmsk, 1));
	}
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);

	return;
}
__inline static void mask_sqrredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
	__m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	__m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	__m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	__m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	__m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	__m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	int biascount = 0;
	__m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	__m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	t1 += sqr_lo;
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(sqr_hi, C2);
	// we throw t0 away after this so first propagate its carry.
	// t0 = _mm512_add_epi64(t1, one);
	// t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);

	carryprop(t2, C2, lo52mask);
	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	//scarry = _mm512_cmpge_epu64_mask(t2, n1);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
__inline static void mask_sqrredc104_vec_pos(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with positive variant n'
	__m512i m;
	__m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
	__m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	__m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	__m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	__m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	__m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	__m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	int biascount = 0;
	__m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	__m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	t1 += sqr_lo;

	// note, we leave rho = 0 - rho so that we get -m,
	// and thus the muladd becomes mulsub, since there is
	// no fmsub52 in avx512-ifma.
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_sub_epi64(t1, C1);
	t2 = _mm512_add_epi64(sqr_hi, C2);
	t0 = t1;
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_sub_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);
	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally add
		t1 = _mm512_mask_add_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_add_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_add_epi64(t2, scarry, t2, _mm512_srli_epi64(t1, 52));
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
__inline static void mask_sqrredc104_exact_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#ifndef IFMA
	__m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	__m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	__m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	__m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	__m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	__m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	int biascount = 0;
	__m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	__m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	t1 = sqr_lo;
	t2 = sqr_hi;
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	// t0 = _mm512_add_epi64(t1, one);
	// t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);

	//carryprop(t2, C2, lo52mask);
	//scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	scarry = _mm512_cmpge_epu64_mask(t2, n1);
	//scarry |= (_mm512_cmpeq_epu64_mask(t2, n1) & _mm512_cmpgt_epu64_mask(t1, n0));

	if (scarry > 0) {
		// conditionally subtract when result >= n
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
__inline static void addmod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
	__m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	// add
	__mmask8 bmsk;
	//a0 = _mm512_addsetc_epi52(a0, b0, &bmsk);
	//a1 = _mm512_adc_epi52(a1, bmsk, b1, &bmsk);
	a0 = _mm512_add_epi64(a0, b0);
	a1 = _mm512_add_epi64(a1, b1);
	a1 = _mm512_add_epi64(a1, _mm512_srli_epi64(a0, 52));
	a0 = _mm512_and_epi64(a0, lo52mask);

	// compare
	//__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);
	// *c0 = _mm512_mask_sub_epi64(a0, msk, a0, n0);
	// *c1 = _mm512_mask_sub_epi64(a1, msk, a1, n1);
	// *c1 = _mm512_mask_sub_epi64(*c1, msk, *c1, _mm512_srli_epi64(*c0, 63));
	return;
}
__inline static void mask_addmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// add
	__mmask8 bmsk;
	a0 = _mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
	a1 = _mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

	// compare
	__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void mask_dblmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	// add
	__mmask8 bmsk;
	//a0 = _mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
	//a1 = _mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

	a0 = _mm512_mask_slli_epi64(a0, addmsk, a0, 1);
	a1 = _mm512_mask_slli_epi64(a1, addmsk, a1, 1);
	// when doubling, it is safe to check both carries before adding
	// in the previous carry, because the shift makes room for
	// the previous carry.  So either the upper word shift generates
	// a carry or doesn't, the addition won't cause one.
	a1 = _mm512_add_epi64(a1, _mm512_srli_epi64(a0, 52));
	a0 = _mm512_and_epi64(lo52mask, a0);

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void mask_redsub104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__mmask8 bmsk;

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void redsub104_x8(__m512i* c1, __m512i* c0,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__mmask8 bmsk;

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);

	return;
}
__inline static void submod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
	__m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// compare
	__mmask8 msk = _mm512_cmplt_epu64_mask(a1, b1);
	msk |= _mm512_cmpeq_epu64_mask(a1, b1) & _mm512_cmplt_epu64_mask(a0, b0);

	// subtract
	__mmask8 bmsk;
	a0 = _mm512_subsetc_epi52(a0, b0, &bmsk);
	a1 = _mm512_sbb_epi52(a1, bmsk, b1, &bmsk);

	// conditionally add N
	*c0 = _mm512_mask_addsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_adc_epi52(a1, msk, bmsk, n1, &bmsk);
	return;
}

uint64_t multiplicative_inverse(uint64_t a);
__m512i multiplicative_inverse104_x8(uint64_t* a);

#endif

#endif
