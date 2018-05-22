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
#include "arith.h"
#include "common.h"

typedef struct
{
    mpz_t rhat;
    mpz_t nhat;
    mpz_t n;
    mpz_t r;
    mpz_t tmp;
} monty_t;



/********************* montgomery arith **********************/
monty_t * monty_init(mpz_t n);
void to_monty(monty_t *mdata, mpz_t x);
void monty_add(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_mul(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_sub(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_free(monty_t *mdata);
void monty_redc(monty_t *mdata, mpz_t x);

#if defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

__inline uint64 addmod(uint64 x, uint64 y, uint64 n)
{
    uint64 t = x - n;
    x += y;
    asm("add %2, %1\n\t"
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
        : "1"(c), "0"(0), "r"(n));

    return n;
}

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

#else

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

__inline uint64 mulredc(uint64 x, uint64 y, uint64 n, uint64 nhat) {}
__inline uint64 u64div(uint64 c, uint64 n) {}
__inline uint64 addmod(uint64 x, uint64 y, uint64 n) {}
__inline uint64 mulredc63(uint64 x, uint64 y, uint64 n, uint64 nhat) {}

#endif

