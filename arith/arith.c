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
www.boo.net/~jasonp, Scott Contini's mpqs implementation, Tom St. 
Denis Tom's Fast Math library, and Jeff Hurchalla's Binary GCD.  Many
thanks to their kind donation of code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "arith.h"
#include "ytools.h"
#include "common.h"
#include <math.h>
#include "mpz_aprcl.h"



uint64_t mpz_get_64(mpz_t src)
{
	uint64_t out = mpz_getlimbn(src, 0);

	return out;
}

void mpz_set_64(mpz_t dest, uint64_t src)
{
    mpz_set_ui(dest, src);
}

int ndigits_1(uint64_t n)
{
    int i = 0;
    while (n != 0)
    {
        n /= 10;
        i++;
    }
    if (i == 0)
        i++;
    return i;
}

int gmp_base10(mpz_t x)
{
    mpz_t t;	//temp
    int g;		//guess: either correct or +1

    mpz_init(t);
    g = mpz_sizeinbase(x, 10);
    mpz_set_ui(t, 10);
    mpz_pow_ui(t, t, g - 1);
    g = g - (mpz_cmp(t, x) > 0 ? 1 : 0);
    mpz_clear(t);
    return g;
}

// borrowed from jasonp... 
double zlog(mpz_t x) {

    uint32_t i = mpz_size(x);

    switch (i) {

    case 0:
        return 0;
    case 1:
        return log((double)((uint32_t)mpz_get_ui(x)));
    case 2:

        return log((double)(mpz_getlimbn(x, 0)) +
            18446744073709551616.0 * mpz_getlimbn(x, 1));

    default:

        return 64 * (i - 3) * LN2 +
            log((double)(mpz_getlimbn(x, i - 3)) + 18446744073709551616.0 * (
            ((double)mpz_getlimbn(x, i - 2) + 18446744073709551616.0 *
                mpz_getlimbn(x, i - 1))));

    }
}

uint64_t spBits(uint64_t n)
{
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }

    return i;
}

int bits64(uint64_t n)
{
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }
    return i;
}


// make a _udiv128 for all supported compilers (Windows: msvc, msvc+clang, msvc+intel; Linux: intel, clang, gcc)
// For MSVC, use the intrinsic
#if defined(_MSC_VER) && (!defined(__clang__))
#include <intrin.h>
#pragma intrinsic(_udiv128)
#else
// For other compilers, we'll implement a fallback or use inline assembly
// This is a simplified implementation for demonstration
static __inline uint64_t _udiv128(uint64_t high, uint64_t low, uint64_t divisor, uint64_t* remainder) {

    // Use __int128 if available (GCC/Clang)
#if (defined(__clang__) || defined(__INTEL_COMPILER)) && defined(GCC_ASM64X)

    __asm__("divq %4"
        : "=a"(low), "=d"(high)
        : "1"(high), "0"(low), "r"(divisor));

    *remainder = high;
    return low;

#elif defined(HAS_UINT128)

    __uint128_t dividend = ((__uint128_t)high << 64) | low;
    __uint128_t quotient = dividend / divisor;
    *remainder = dividend % divisor;
    return (uint64_t)quotient;

#else

// Fallback implementation using double-precision division
// This is not as accurate but works for demonstration
    double d_dividend = (double)high * (1.0 * (1ULL << 32) * (1ULL << 32)) + (double)low;
    double d_quotient = d_dividend / (double)divisor;
    uint64_t quotient = (uint64_t)d_quotient;

    // brb added.  I think all compiler cases that I support are covered but in
    // case not, spam this message.
    printf("warning: _udiv128 is not accurate\n");

    // Calculate remainder more accurately
    uint64_t temp_high, temp_low;
    // Multiply quotient * divisor
    temp_low = quotient * divisor; // This may overflow, but that's handled below

    // Subtract from original dividend to get remainder
    if (low >= temp_low) {
        *remainder = low - temp_low;
    }
    else {
        *remainder = (UINT64_MAX - temp_low + 1) + low;
        quotient--;
    }

    return quotient;

#endif
}

#endif


// make a _umul128 for all supported compilers (Windows: msvc, msvc+clang, msvc+intel; Linux: intel, clang, gcc)
#if defined(_MSC_VER) || defined(__MINGW32__)
#include <immintrin.h>
#pragma intrinsic(_umul128)
#else
// For other compilers, we'll implement a fallback or use inline assembly
// This is a simplified implementation for demonstration
static __inline uint64_t _umul128(uint64_t x, uint64_t y, uint64_t* hi) {

#if (defined(__clang__) || defined(__INTEL_COMPILER)) && defined(GCC_ASM64X)

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

// Use __int128 if available (GCC/Clang)
#elif defined(HAS_UINT128)

    __uint128_t prod = (__uint128_t)x * (__uint128_t)y;
    *hi = (uint64_t)(prod >> 64);
    return (uint64_t)prod;

#else

    // Fallback: split into 32-bit parts
    uint64_t x_low = x & 0xFFFFFFFF;
    uint64_t x_high = x >> 32;
    uint64_t y_low = y & 0xFFFFFFFF;
    uint64_t y_high = y >> 32;

    uint64_t p0 = x_low * y_low;
    uint64_t p1 = x_low * y_high;
    uint64_t p2 = x_high * y_low;
    uint64_t p3 = x_high * y_high;

    uint64_t middle = p1 + p2 + (p0 >> 32);
    *hi = p3 + (middle >> 32);
    return p0 + (middle << 32);


#endif
}

#endif

#if defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)


#if defined(__INTEL_COMPILER) || defined(__clang__)
#define ASM_ __asm__
#elif defined(_WIN32)
#define ASM_ ASM_M
#else
#define ASM_ ASM_G
#endif

// uint64_t u64div(uint64_t c, uint64_t n)
// {
// #if 1
//     __asm__("divq %4"
//         : "=a"(c), "=d"(n)
//         : "1"(c), "0"(0ULL), "r"(n));
// #else
//     // this should work if the above won't compile (e.g. on clang)
//     uint64_t tmp = 0;
//     __asm__("divq %4"
//         : "=a"(tmp), "=d"(n)
//         : "1"(c), "0"(tmp), "r"(n));
// #endif
//     return n;
// }


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

void spAdd(uint64_t u, uint64_t v, uint64_t* sum, uint64_t* carry)
{
    //fp_word s,c;
    uint64_t s, c;

    s = v;
    c = 0;

    ASM_("movq %2, %%rax		\n\t"
        "addq %%rax, %3		\n\t"
        "adcq $0, %4		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(u), "0"(s), "1"(c)
        : "rax", "memory", "cc");

    *sum = s;
    *carry = c;

    return;
}

void spAdd3(uint64_t u, uint64_t v, uint64_t w, uint64_t* sum, uint64_t* carry)
{
    //fp_word s,c;
    uint64_t s, c;

    s = v;
    c = 0;

    ASM_("movq %2, %%rax		\n\t"
        "addq %3, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "addq %%rax, %4		\n\t"
        "adcq $0, %5		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(u), "r"(w), "0"(s), "1"(c)
        : "rax", "memory", "cc");

    *sum = s;
    *carry = c;

    return;
}

void spSub3(uint64_t u, uint64_t v, uint64_t w, uint64_t* sub, uint64_t* borrow)
{
    //fp_word s,b;
    uint64_t s, b;

    s = v;
    b = 0;

    ASM_("movq %2, %%rax		\n\t"
        "subq %4, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "subq %3, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "movq %%rax, %4		\n\t"
        : "=r"(s), "=r"(b)
        : "r"(u), "r"(w), "0"(s), "1"(b)
        : "rax", "memory", "cc");

    *sub = s;
    *borrow = b;

    return;
}

void spSub(uint64_t u, uint64_t v, uint64_t* sub, uint64_t* borrow)
{
    //fp_word s,b;
    uint64_t s, b;

    s = v;
    b = 0;

    ASM_("movq %2, %%rax		\n\t"
        "subq %3, %%rax		\n\t"
        "adcq $0, %4		\n\t"
        "movq %%rax, %3		\n\t"
        : "=r"(s), "=r"(b)
        : "r"(u), "0"(s), "1"(b)
        : "rax", "memory", "cc");

    *sub = s;
    *borrow = b;

    return;
}

uint64_t spDivide(uint64_t* q, uint64_t* r, uint64_t u[2], uint64_t v)
{
    *r = u[1];
    *q = u[0];
    ASM_("divq %4"
        : "=a"(*q), "=d"(*r)
        : "1"(*r), "0"(*q), "r"(v));

    return 0;
}

void spMultiply(uint64_t u, uint64_t v, uint64_t* product, uint64_t* carry)
{
    *product = v;
    *carry = u;

    ASM_("movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product)
        : "rax", "rdx", "cc");

    return;
}

uint64_t spPRP2(uint64_t p)
{
    // do a base-2 prp test on the input, where p is greater than 2^32
    // i.e., compute 2^(p-1) % p.
    // since p is more than 32 bits we can do the accumulation division 
    // free for the first 5 iterations.  may not be much, but it's something.
    uint64_t result;

    ASM_(
        "xorq	%%rbx, %%rbx \n\t"
        "xorq	%%rdi, %%rdi \n\t"
        "addq	$1, %%rdi \n\t"		/* n = 1 */
        "0:	\n\t"					/* begin loop */
        "test	$1, %%rcx \n\t"		/* exp & 0x1 */
        "je	2f		\n\t"			/* bit not set, skip accumulation into n */
        "movq	%%rax, %%rsi \n\t"	/* save acc */
        "mulq	%%rdi \n\t"			/* n * acc mod m */
        "movq	%%rax, %%rdi \n\t"	/* save n */
        "movq	%%rsi, %%rax \n\t"	/* restore acc */
        "2:			\n\t"			/* square acc stage */
        "shrq	$1, %%rcx \n\t"		/* base >>= 1 */
        "addq	$1, %%rbx \n\t"
        "mulq	%%rax \n\t"			/* acc = acc * acc*/
        "cmpq	$5, %%rbx \n\t"		/* 5 iterations? */
        "jb 0b \n\t"
        "3:	\n\t"					/* begin loop */
        "test	$1, %%rcx \n\t"		/* exp & 0x1 */
        "je	4f		\n\t"			/* bit not set, skip accumulation into n */
        "movq	%%rax, %%rsi \n\t"	/* save acc */
        "mulq	%%rdi \n\t"			/* n * acc mod m */
        "divq	%3 \n\t"
        "movq	%%rdx, %%rdi \n\t"	/* save n */
        "movq	%%rsi, %%rax \n\t"	/* restore acc */
        "4:			\n\t"			/* square acc stage */
        "shrq	$1, %%rcx \n\t"		/* base >>= 1 */
        "mulq	%%rax \n\t"			/* acc = acc * acc*/
        "divq	%3 \n\t"
        "cmpq	$0, %%rcx \n\t"		/* exp == 0? */
        "movq	%%rdx, %%rax \n\t"	/* mod m */
        "jne 3b \n\t"
        "movq	%%rdi, %0 \n\t"
        : "=r"(result)
        : "a"(2), "c"(p - 1), "r"(p)
        : "rbx", "rdx", "rdi", "rsi", "cc");				/* return result */

    return result;

}

uint64_t spModExp_asm(uint64_t b, uint64_t e, uint64_t m)
{
    uint64_t result;


    ASM_(
        "xorq	%%rdi, %%rdi \n\t"
        "addq	$1, %%rdi \n\t"		/* n = 1 */
        "cmpq	$0, %%rcx \n\t"		/* exp == 0? */
        "je 1f \n\t"
        "0:	\n\t"					/* begin loop */
        "test	$1, %%rcx \n\t"		/* exp & 0x1 */
        "je	2f		\n\t"			/* bit not set, skip accumulation into n */
        "movq	%%rax, %%rsi \n\t"	/* save acc */
        "mulq	%%rdi \n\t"			/* n * acc mod m */
        "divq	%3 \n\t"
        "movq	%%rdx, %%rdi \n\t"	/* save n */
        "movq	%%rsi, %%rax \n\t"	/* restore acc */
        "2:			\n\t"			/* square acc stage */
        "shrq	$1, %%rcx \n\t"		/* base >>= 1 */
        "mulq	%%rax \n\t"			/* acc = acc * acc*/
        "divq	%3 \n\t"
        "cmpq	$0, %%rcx \n\t"		/* exp == 0? */
        "movq	%%rdx, %%rax \n\t"	/* mod m */
        "jne 0b \n\t"
        "1:			\n\t"			/* end loop */
        "movq	%%rdi, %0 \n\t"
        : "=r"(result)
        : "a"(b), "c"(e), "r"(m)
        : "rdx", "rdi", "rsi", "cc");				/* return result */

    return result;
}

void spMulAdd(uint64_t u, uint64_t v, uint64_t w, uint64_t t, uint64_t* lower, uint64_t* carry)
{
    uint64_t k, p;
    spMultiply(u, v, &p, carry);
    spAdd3(p, w, t, lower, &k);
    *carry += k;
    return;
}

void spMulMod(uint64_t u, uint64_t v, uint64_t m, uint64_t* w)
{
    uint64_t p[2];
    uint64_t q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}


#else

uint64_t spDivide(uint64_t * q, uint64_t * r, uint64_t u[2], uint64_t v)
{
    *q = _udiv128(u[1], u[0], v, r);
    return 0;
}

void spMultiply(uint64_t u, uint64_t v, uint64_t* product, uint64_t* carry)
{
    *product = _umul128(u, v, carry);
    return;
}

void spAdd(uint64_t u, uint64_t v, uint64_t* sum, uint64_t* carry)
{
    *carry = _addcarry_u64(0, u, v, sum);
    return;
}

void spAdd3(uint64_t u, uint64_t v, uint64_t w, uint64_t* sum, uint64_t* carry)
{
    unsigned char c;
    *carry = _addcarry_u64(0, u, v, sum);
    c = _addcarry_u64(*carry, *sum, w, sum);
    *carry += c;
    return;
}

void spSub3(uint64_t u, uint64_t v, uint64_t w, uint64_t* sub, uint64_t* borrow)
{
    unsigned char b;
    *borrow = _subborrow_u64(0, u, v, sub);
    b = _subborrow_u64(0, *sub, w, sub);
    *borrow += b;
    return;
}

void spSub(uint64_t u, uint64_t v, uint64_t* sub, uint64_t* borrow)
{
    *borrow = _subborrow_u64(0, u, v, sub);
    return;
}

void spMulAdd(uint64_t u, uint64_t v, uint64_t w, uint64_t t, uint64_t* lower, uint64_t* carry)
{
    uint64_t k, p;
    spMultiply(u, v, &p, carry);
    spAdd3(p, w, t, lower, &k);
    *carry += k;
    return;
}

void spMulMod(uint64_t u, uint64_t v, uint64_t m, uint64_t* w)
{
    uint64_t p[2];
    uint64_t q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}

#endif




// make a 128-bit add/sub for all supported compilers (Windows: msvc, msvc+clang, msvc+intel; Linux: intel, clang, gcc)
// created with help from Claude-code
// 
// Platform detection for compiler intrinsics
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#define COMPILER_MSVC
#elif defined(__GNUC__) || defined(__clang__)
#define COMPILER_GCC_LIKE
#endif


#ifdef HAS_UINT128

static inline int uint128_add(uint128_t* result, const uint128_t a, const uint128_t b)
{
    *result = a + b;
    return (*result < a);
}

static inline int uint128_sub(uint128_t* result, const uint128_t a, const uint128_t b)
{
    *result = a - b;
    return a < b;
}

#else
// 128-bit unsigned integer represented as array of two 64-bit integers
// [0] = low 64 bits, [1] = high 64 bits (little-endian order)
//typedef uint64_t uint128_t[2];
typedef uint64_t uint128_t[2];

/**
 * Add two 128-bit unsigned integers
 * @param result: Output array to store the sum
 * @param a: First operand (128-bit)
 * @param b: Second operand (128-bit)
 * @return: Carry bit (0 or 1)
 */
static inline int uint128_add(uint128_t result, const uint128_t a, const uint128_t b) {
    uint64_t carry = 0;

#if defined(COMPILER_MSVC)
    // Use MSVC intrinsics for optimal performance
    unsigned char c1 = _addcarry_u64(0, a[0], b[0], &result[0]);
    unsigned char c2 = _addcarry_u64(c1, a[1], b[1], &result[1]);
    carry = c2;
#elif defined(COMPILER_GCC_LIKE)
    // Use GCC/Clang builtin overflow detection
    // Use GCC/Clang builtin overflow detection
    carry = __builtin_add_overflow(a[0], b[0], &result[0]);
    result[1] = __builtin_addcll(a[1], b[1], carry, &carry);
#else
    // Portable fallback implementation
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];

    // Check for overflow in low part
    if (result[0] < a[0]) {
        result[1]++;
        // Check for overflow in high part after increment
        if (result[1] < a[1]) {
            carry = 1;
        }
    }
    else if (result[1] < a[1]) {
        // Overflow in high part without low overflow
        carry = 1;
    }
#endif

    return (int)carry;
}

/**
 * Add two 128-bit unsigned integers with input carry
 * @param result: Output array to store the sum
 * @param a: First operand (128-bit)
 * @param b: Second operand (128-bit)
 * @param carry_in: Input carry (0 or 1)
 * @return: Output carry bit (0 or 1)
 */
static inline int uint128_adc(uint128_t result, const uint128_t a, const uint128_t b, int carry_in) {
    uint64_t carry = 0;

#if defined(COMPILER_MSVC)
    unsigned char c1 = _addcarry_u64(carry_in, a[0], b[0], &result[0]);
    unsigned char c2 = _addcarry_u64(c1, a[1], b[1], &result[1]);
    carry = c2;
#else
    // First add the operands
    int c = uint128_add(result, a, b);

    // Then add the carry_in
    if (carry_in) {
        uint128_t one = { 1, 0 };
        c |= uint128_add(result, result, one);
    }
    carry = c;
#endif

    return (int)carry;
}

/**
 * Subtract two 128-bit unsigned integers
 * @param result: Output array to store the difference (a - b)
 * @param a: Minuend (128-bit)
 * @param b: Subtrahend (128-bit)
 * @return: Borrow bit (0 or 1), 1 indicates underflow (b > a)
 */
static inline int uint128_sub(uint128_t result, const uint128_t a, const uint128_t b) {
    uint64_t borrow = 0;

#if defined(COMPILER_MSVC)
    // Use MSVC intrinsics for optimal performance
    unsigned char b1 = _subborrow_u64(0, a[0], b[0], &result[0]);
    unsigned char b2 = _subborrow_u64(b1, a[1], b[1], &result[1]);
    borrow = b2;
#elif defined(COMPILER_GCC_LIKE)
    // Use GCC/Clang builtin overflow detection
    borrow = __builtin_sub_overflow(a[0], b[0], &result[0]);
    result[1] = __builtin_subcll(a[1], b[1], borrow, &borrow);
#else
    // Portable fallback implementation
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];

    // Check for borrow from low part
    if (a[0] < b[0]) {
        result[1]--;
        // Check for underflow in high part after decrement
        if (a[1] < b[1] || result[1] > a[1]) {
            borrow = 1;
        }
    }
    else if (a[1] < b[1]) {
        // Underflow in high part without low borrow
        borrow = 1;
    }
#endif

    return (int)borrow;
}

/**
 * Subtract two 128-bit unsigned integers with input borrow
 * @param result: Output array to store the difference
 * @param a: Minuend (128-bit)
 * @param b: Subtrahend (128-bit)
 * @param borrow_in: Input borrow (0 or 1)
 * @return: Output borrow bit (0 or 1)
 */
static inline int uint128_sbb(uint128_t result, const uint128_t a, const uint128_t b, int borrow_in) {
    uint64_t borrow = 0;

#if defined(COMPILER_MSVC)
    unsigned char b1 = _subborrow_u64(borrow_in, a[0], b[0], &result[0]);
    unsigned char b2 = _subborrow_u64(b1, a[1], b[1], &result[1]);
    borrow = b2;
#else
    // First subtract the operands
    int b = uint128_sub(result, a, b);

    // Then subtract the borrow_in
    if (borrow_in) {
        uint128_t one = { 1, 0 };
        b |= uint128_sub(result, result, one);
    }
    borrow = b;
#endif

    return (int)borrow;
}

/**
 * Compare two 128-bit unsigned integers
 * @param a: First operand
 * @param b: Second operand
 * @return: -1 if a < b, 0 if a == b, 1 if a > b
 */
static inline int uint128_compare(const uint128_t a, const uint128_t b) {
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    return 0;
}

/**
 * Check if 128-bit integer is zero
 * @param a: Operand to check
 * @return: true if zero, false otherwise
 */
static inline int uint128_is_zero(const uint128_t a) {
    return (a[0] == 0) && (a[1] == 0);
}

/**
 * Set 128-bit integer to zero
 * @param result: Output to set to zero
 */
static inline void uint128_zero(uint128_t result) {
    result[0] = 0;
    result[1] = 0;
}

/**
 * Set 128-bit integer from a 64-bit value
 * @param result: Output 128-bit integer
 * @param value: Input 64-bit value
 */
static inline void uint128_set64(uint128_t result, uint64_t value) {
    result[0] = value;
    result[1] = 0;
}

/**
 * Copy 128-bit integer
 * @param dest: Destination
 * @param src: Source
 */
static inline void uint128_copy(uint128_t dest, const uint128_t src) {
    dest[0] = src[0];
    dest[1] = src[1];
}


#endif




/*
Claude-code generated 128-bit / 128-bit division.  Here was the conversation:


===================================================================================================
Me:
I'd like to write a function in C that divides one 128-bit integer, represented as an array of two 64-bit integers,
by another 128-bit integer, represented as an array of two 64-bit integers, returning the quotient and remainder.

Claude:
I'll help you create a function to divide two 128-bit integers represented as arrays of two 64-bit integers.
This requires implementing long division since C doesn't have native 128-bit division on all platforms.

128-bit Integer Division
Code
This implementation provides a complete 128-bit division function with the following key features:

Main Function:

uint128_div() - Takes dividend and divisor as arrays of two 64-bit integers, returns quotient and remainder
Returns 0 on success, -1 on division by zero
Key Implementation Details:

Data Structure: Uses a uint128_t structure internally for easier manipulation
Algorithm: Implements binary long division, similar to how you'd do division by hand but in base 2
Edge Cases: Handles division by zero, zero dividend, and cases where dividend < divisor
Bit Manipulation: Uses efficient bit shifting and comparison operations
Helper Functions:

is_zero() - Check if a 128-bit number is zero
compare() - Compare two 128-bit numbers
shift_left_1()/shift_right_1() - Bit shifting operations
subtract() - 128-bit subtraction
find_msb() - Find most significant bit position
Usage:

uint64_t dividend[2] = {low_bits, high_bits};
uint64_t divisor[2] = {low_bits, high_bits};
uint64_t quotient[2], remainder[2];

int result = uint128_div(dividend, divisor, quotient, remainder);
The algorithm works by:

Aligning the divisor with the most significant bits of the dividend
Performing subtraction when possible and setting quotient bits
Shifting right and repeating until complete
The example code includes test cases demonstrating basic division, large number division, and error handling for
division by zero.

===================================================================================================
Me:
Can you rewrite to make use of Knuth algorithm D?

Claude:
I'll rewrite the division function to use Knuth's Algorithm D from "The Art of Computer Programming" Volume 2. This is
the classic algorithm for multi-precision division.

128-bit Integer Division
Code‚Version 2
This implementation now uses Knuth's Algorithm D for multi-precision division. Here are the key improvements and features:

Knuth Algorithm D Implementation:

Normalization (Step D1): Multiplies both dividend and divisor by a factor to ensure the leading digit of the
divisor is BASE/2, which improves the accuracy of quotient digit estimation.
Quotient Digit Estimation (Step D3): Uses the two leading digits of the current dividend portion and the
leading digit of the divisor to estimate the quotient digit (q).
Adjustment: Refines the quotient digit estimate to prevent overflow.
Multiply and Subtract (Step D4): Performs the trial subtraction.
Correction (Steps D5-D6): If the subtraction goes negative, adds back and decreases the quotient digit.
Unnormalization (Step D8): Removes the normalization factor from the remainder.
Key Features:

32-bit Digit Arithmetic: Works with 32-bit "digits" as Knuth recommends, making the algorithm more manageable
Optimized Single-Digit Case: Special handling for divisors that fit in a single 32-bit digit
Proper Normalization: Implements Knuth's normalization technique for accuracy
Error Handling: Robust handling of edge cases (zero division, zero dividend, etc.)
Algorithm Advantages:

Mathematically Sound: Based on Knuth's rigorous analysis
Efficient: O(mn) complexity where m and n are the number of digits
Accurate: Normalization ensures quotient digit estimates are almost always correct
Well-Tested: This is the standard algorithm used in most multiple-precision libraries
Usage Example:

uint64_t dividend[2] = {0xFFFFFFFFFFFFFFFFULL, 0x1ULL};
uint64_t divisor[2] = {0x100000000ULL, 0};
uint64_t quotient[2], remainder[2];

int result = uint128_div(dividend, divisor, quotient, remainder);
The implementation handles all the tricky aspects of Knuth's algorithm including proper carry propagation,
normalization, and the add-back correction step that's needed when the quotient digit estimate is too large.

===================================================================================================
Me:
ok great.  now please modify this code to use 64-bit dividends and the function _udiv128 to perform
128-bit by 64-bit divisions when necessary


I'll modify the code to use 64-bit digits and leverage the _udiv128 intrinsic for efficient 128-bit by 64-bit
divisions.This updated implementation now uses 64-bit digits and leverages the _udiv128 intrinsic for efficient
divisions. Here are the key changes and improvements:

**Major Changes:**

1. **64-bit Digits**: Now works with 64-bit digits instead of 32-bit, reducing the number of operations needed.

2. **_udiv128 Integration**: Uses the hardware-accelerated _udiv128 intrinsic for 128-bit by 64-bit divisions,
which is much faster than implementing this in software.

3. **Cross-Platform Support**: Includes fallback implementations for non-MSVC compilers using __int128 or manual computation.

**Key Optimizations:**

1. **Single-Digit Divisor Fast Path**: When the divisor fits in a single 64-bit digit, uses _udiv128 directly for maximum efficiency.

2. **Efficient Normalization**: Uses bit-level operations optimized for 64-bit arithmetic.

3. **Hardware Multiplication**: Uses _umul128 intrinsic where available for 64-bit x 64-bit -> 128-bit multiplication.

**Algorithm Flow:**

1. **Input Validation**: Checks for division by zero and handles edge cases.

2. **Single 64-bit Divisor**: If divisor < 2^64, uses _udiv128 for direct hardware division:
   q[1] = _udiv128(0, u[1], v[0], &remainder);
   q[0] = _udiv128(remainder, u[0], v[0], &final_remainder);

3. **Multi-Digit Division**: Uses full Knuth Algorithm D with 64-bit digits for divisors >= 2^64.

**Performance Benefits:**

- **Hardware Acceleration**: _udiv128 provides native CPU support for 128 / 64 operations
- **Reduced Iterations**: 64-bit digits mean fewer loop iterations compared to 32-bit version
- **Better Cache Usage**: Fewer memory operations due to larger digit size
- **Compiler Optimizations**: Modern compilers can better optimize 64-bit arithmetic

**Compiler Support:**

- **MSVC**: Uses _udiv128 and _umul128 intrinsics
- **GCC/Clang**: Uses __int128 builtin type when available
- **Fallback**: Provides software implementation for other cases

This implementation should be significantly faster than the previous 32-bit digit version, especially for cases
where the divisor fits in 64 bits, which can take advantage of the hardware _udiv128 instruction.

BRB:
note that I've modified Claude's code a little to make use of existing low level code.

*/





#include <stdint.h>
#include <stdbool.h>

// Knuth Algorithm D for 128-bit division using 64-bit digits
// Based on "The Art of Computer Programming" Volume 2, Section 4.3.1

#define BASE64 (1ULL << 63) // We'll use 2^63 as our effective base to avoid overflow issues

// Convert 128-bit number (two uint64_t) to array of 64-bit digits (little-endian)
static void to_digits_64(const uint64_t num[2], uint64_t digits[2]) {
    digits[0] = num[0]; // Low 64 bits
    digits[1] = num[1]; // High 64 bits
}

// Convert array of 64-bit digits back to 128-bit number
static void from_digits_64(const uint64_t digits[2], uint64_t num[2]) {
    num[0] = digits[0];
    num[1] = digits[1];
}

// Find the actual length of a number in 64-bit digits (ignoring leading zeros)
static int actual_length_64(const uint64_t digits[], int max_len) {
    for (int i = max_len - 1; i >= 0; i--) {
        if (digits[i] != 0) {
            return i + 1;
        }
    }
    return 0; // All zeros
}

// Shift left by n bits for multi-digit number
static void shift_left_64(uint64_t digits[], int len, int shift) {
    if (shift == 0 || len == 0) return;

    if (shift >= 64) {
        // Shift by whole digits
        int digit_shift = shift / 64;
        shift %= 64;

        for (int i = len - 1; i >= digit_shift; i--) {
            digits[i] = digits[i - digit_shift];
        }
        for (int i = 0; i < digit_shift; i++) {
            digits[i] = 0;
        }
        if (shift == 0) return;
    }

    // Shift by remaining bits
    uint64_t carry = 0;
    for (int i = 0; i < len; i++) {
        uint64_t new_carry = digits[i] >> (64 - shift);
        digits[i] = (digits[i] << shift) | carry;
        carry = new_carry;
    }
}

// Shift right by n bits for multi-digit number
static void shift_right_64(uint64_t digits[], int len, int shift) {
    if (shift == 0 || len == 0) return;

    if (shift >= 64) {
        // Shift by whole digits
        int digit_shift = shift / 64;
        shift %= 64;

        for (int i = 0; i < len - digit_shift; i++) {
            digits[i] = digits[i + digit_shift];
        }
        for (int i = len - digit_shift; i < len; i++) {
            digits[i] = 0;
        }
        if (shift == 0) return;
    }

    // Shift by remaining bits
    uint64_t carry = 0;
    for (int i = len - 1; i >= 0; i--) {
        uint64_t new_carry = digits[i] << (64 - shift);
        digits[i] = (digits[i] >> shift) | carry;
        carry = new_carry;
    }
}

// Subtract: u = u - v * q, return borrow
static bool sub_mul_64(uint64_t u[], const uint64_t v[], int n, uint64_t q) {
    uint64_t borrow = 0;

    for (int i = 0; i < n; i++) {
        // Calculate v[i] * q
        uint64_t high, low;

        low = _umul128(v[i], q, &high);

        // Add previous borrow to the product
        low += borrow;
        if (low < borrow) high++;

        // Subtract from u[i]
        if (u[i] >= low) {
            u[i] -= low;
            borrow = high;
        }
        else {
            u[i] = u[i] - low; // Will underflow correctly
            borrow = high + 1;
        }
    }

    return borrow > 0;
}

// Add: u = u + v, return carry
static uint64_t add_64(uint64_t u[], const uint64_t v[], int n) {
    uint64_t carry = 0;

    for (int i = 0; i < n; i++) {
        uint64_t sum = u[i] + v[i] + carry;
        carry = (sum < u[i]) ? 1 : 0;
        u[i] = sum;
    }

    return carry;
}

// Knuth Algorithm D implementation using 64-bit digits and _udiv128
int knuth_div_64(uint64_t u[], int m_plus_n, const uint64_t v[], int n,
    uint64_t q[], uint64_t r[]) {

    // Step D1: Normalize
    // Find the normalization factor d
    int shift = _lead_zcnt64(v[n - 1]);

    // Shift u and v left by shift bits
    uint64_t u_norm[4] = { 0 }; // Extended for normalization
    uint64_t v_norm[2] = { 0 };

    // Copy and extend u
    for (int i = 0; i < m_plus_n; i++) {
        u_norm[i] = u[i];
    }
    shift_left_64(u_norm, m_plus_n + 1, shift);

    // Copy and normalize v
    for (int i = 0; i < n; i++) {
        v_norm[i] = v[i];
    }
    shift_left_64(v_norm, n, shift);

    int m = m_plus_n - n; // Number of quotient digits

    // Step D2: Initialize j
    for (int j = m; j >= 0; j--) {
        // Step D3: Calculate q‚ (q-hat) using _udiv128
        uint64_t q_hat;
        uint64_t r_hat;

        if (u_norm[j + n] == v_norm[n - 1]) {
            q_hat = 0xffffffffffffffffull;
            r_hat = u_norm[j + n - 1] + v_norm[n - 1];
        }
        else {
            // Use _udiv128 for 128-bit by 64-bit division
            q_hat = _udiv128(u_norm[j + n], u_norm[j + n - 1], v_norm[n - 1], &r_hat);
        }

        // Adjust q_hat if necessary (Step D3 continued)
        while (n >= 2 && q_hat > 0) {
            // Check if q_hat * v_norm[n-2] > (r_hat << 64) + u_norm[j + n - 2]
            uint64_t product_high, product_low;

            product_low = _umul128(q_hat, v_norm[n - 2], &product_high);

            if (product_high > r_hat ||
                (product_high == r_hat && product_low > u_norm[j + n - 2])) {
                q_hat--;
                r_hat += v_norm[n - 1];
                if (r_hat < v_norm[n - 1]) break; // Overflow
            }
            else {
                break;
            }
        }

        // Step D4: Multiply and subtract
        bool negative = sub_mul_64(&u_norm[j], v_norm, n, q_hat);

        // Step D5: Test remainder
        if (negative) {
            // Step D6: Add back
            q_hat--;
            add_64(&u_norm[j], v_norm, n);
        }

        // Store quotient digit
        if (j < 2) q[j] = q_hat;

        // Step D7: Loop on j (handled by for loop)
    }

    // Step D8: Unnormalize remainder
    for (int i = 0; i < n; i++) {
        r[i] = u_norm[i];
    }
    shift_right_64(r, n, shift);

    return 0;
}

// Main 128-bit division function using Knuth Algorithm D with 64-bit digits
int uint128_div(const uint64_t dividend[2], const uint64_t divisor[2],
    uint64_t quotient[2], uint64_t remainder[2]) {

    // Check for division by zero
    if (divisor[0] == 0 && divisor[1] == 0) {
        return -1;
    }

    // Convert to 64-bit digits
    uint64_t u[3] = { 0 }; // dividend + space for normalization
    uint64_t v[2] = { 0 }; // divisor
    uint64_t q[2] = { 0 }; // quotient
    uint64_t r[2] = { 0 }; // remainder

    to_digits_64(dividend, u);
    to_digits_64(divisor, v);

    // Find actual lengths
    int u_len = actual_length_64(u, 2);
    int v_len = actual_length_64(v, 2);

    // Handle special cases
    if (u_len == 0) {
        // Dividend is zero
        quotient[0] = quotient[1] = 0;
        remainder[0] = remainder[1] = 0;
        return 0;
    }

    if (u_len < v_len) {
        // Dividend < divisor
        quotient[0] = quotient[1] = 0;
        remainder[0] = dividend[0];
        remainder[1] = dividend[1];
        return 0;
    }

    // Handle single-digit divisor case using _udiv128
    if (v_len == 1) {
        uint64_t rem;
        if (u_len == 1) {
            q[0] = u[0] / v[0];
            r[0] = u[0] % v[0];
            r[1] = 0;
        }
        else {
            // Use _udiv128 for 128-bit by 64-bit division
            q[1] = _udiv128(0, u[1], v[0], &rem);
            q[0] = _udiv128(rem, u[0], v[0], &r[0]);
            r[1] = 0;
        }
        q[1] = (u_len > 1) ? q[1] : 0;
    }
    else {
        // Multi-digit case: use full Knuth Algorithm D
        knuth_div_64(u, u_len, v, v_len, q, r);
    }

    // Convert results back to 64-bit format
    from_digits_64(q, quotient);
    from_digits_64(r, remainder);

    return 0;
}



uint64_t u64div(uint64_t c, uint64_t n)
{
    uint64_t r;
    _udiv128(c, 0, n, &r);

    return r;
}



void spModExp(uint64_t a, uint64_t b, uint64_t m, uint64_t* u)
{
    //computes a^b mod m = u using the binary method
    //see, for instance, the handbook of applied cryptography
    uint64_t n, bb, aa, t, prod[2];

    n = 1;
    aa = a;
    bb = b;
    while (bb != 0)
    {
        if (bb & 0x1)
        {
            spMultiply(aa, n, &prod[0], &prod[1]);		//n*a
            spDivide(&t, &n, prod, m);					//n*a mod m
        }
        bb >>= 1;
        //compute successive squares of a
        spMultiply(aa, aa, &prod[0], &prod[1]);
        spDivide(&t, &aa, prod, m);
    }
    *u = n;

    return;
}

void ShanksTonelli_1(uint64_t a, uint64_t p, uint64_t* sq)
{
    //a is a quadratic residue mod p
    //p is an odd prime
    //find x where x^2 == a mod p
    //we assume p will always fit into an uint64_t, therefore x will as well.
    //see paper by Ezra Brown
    uint64_t x = 0, b = 0, g = 0, n = 0, s = 0, r = 0, e = 0, b2m = 0, tmp = 0;
    int i;

    //factor p-1 = Q*2^S, where Q is odd and S >= 1.
    s = p - 1;
    e = 0;
    while (!(s & 1))
    {
        s >>= 1;
        e++;
    }

    //find a quadratic non-residue mod p.  keep it small to reduce the work of modexp
    n = 3;
    while (1)
    {
        if (jacobi_1(n, p) < 0)
            break;
        n++;
    }

    //approximate the root x = a^[(s+1)/2] mod p
    spModExp(a, (s + 1) / 2, p, &x);

    //guess at fudge factor b = a^s
    spModExp(a, s, p, &b);

    //initialize g = n^s
    spModExp(n, s, p, &g);

    //initialize r = e
    r = e;

    while (1)
    {
        //find m such that b^(2^m) == 1 mod p with m between 0 and r-1
        b2m = b;
        for (i = 0; i < (int)r; i++)
        {
            if (b2m == 1) break;

            //successivly square b mod p
            spMulMod(b2m, b2m, p, &b2m);
        }

        if (i == 0)
        {
            *sq = x;
            goto free;
        }

        //replace x by x*g^(2^(r-m-1))
        spModExp(g, 1 << (r - i - 1), p, &tmp);
        spMulMod(tmp, x, p, &x);

        //replace g by g^(2^(r-m)) and
        //replace b by b*g
        spModExp(g, 1 << (r - i), p, &g);
        spMulMod(g, b, p, &b);

        r = i;
    }

free:
    //return the smallest solution always
    if (*sq > (p >> 1))
        * sq = p - *sq;

    return;
}

uint32_t modinv_1(uint32_t a, uint32_t p) {

    /* thanks to the folks at www.mersenneforum.org */

    uint32_t ps1, ps2, parity, dividend, divisor, rem, q, t;


    q = 1;
    rem = a;
    dividend = p;
    divisor = a;
    ps1 = 1;
    ps2 = 0;
    parity = 0;

    while (divisor > 1) {
        rem = dividend - divisor;
        t = rem - divisor;
        if (rem >= divisor) {
            q += ps1; rem = t; t -= divisor;
            if (rem >= divisor) {
                q += ps1; rem = t; t -= divisor;
                if (rem >= divisor) {
                    q += ps1; rem = t; t -= divisor;
                    if (rem >= divisor) {
                        q += ps1; rem = t; t -= divisor;
                        if (rem >= divisor) {
                            q += ps1; rem = t; t -= divisor;
                            if (rem >= divisor) {
                                q += ps1; rem = t; t -= divisor;
                                if (rem >= divisor) {
                                    q += ps1; rem = t; t -= divisor;
                                    if (rem >= divisor) {
                                        q += ps1; rem = t;
                                        if (rem >= divisor) {
                                            q = dividend / divisor;
                                            rem = dividend % divisor;
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    if (parity == 0)
        return ps1;
    else
        return p - ps1;
}

uint32_t modinv_1b(uint32_t a, uint32_t p) {

    /* thanks to the folks at www.mersenneforum.org */

    /* modification: p is fixed at 2^32.  a is only valid if odd */

    uint64_t dividend = (uint64_t)0x1 << 32;
    uint32_t ps1, ps2, parity, divisor, rem, q, t;

    q = 1;
    rem = a;
    //dividend = p;
    divisor = a;
    ps1 = 1;
    ps2 = 0;
    parity = 0;

    while (divisor > 1) {
        rem = (uint32_t)(dividend - (uint64_t)divisor);
        t = rem - divisor;
        if (rem >= divisor) {
            q += ps1; rem = t; t -= divisor;
            if (rem >= divisor) {
                q += ps1; rem = t; t -= divisor;
                if (rem >= divisor) {
                    q += ps1; rem = t; t -= divisor;
                    if (rem >= divisor) {
                        q += ps1; rem = t; t -= divisor;
                        if (rem >= divisor) {
                            q += ps1; rem = t; t -= divisor;
                            if (rem >= divisor) {
                                q += ps1; rem = t; t -= divisor;
                                if (rem >= divisor) {
                                    q += ps1; rem = t; t -= divisor;
                                    if (rem >= divisor) {
                                        q += ps1; rem = t;
                                        if (rem >= divisor) {
                                            q = (uint32_t)(dividend / (uint64_t)divisor);
                                            rem = (uint32_t)(dividend % (uint64_t)divisor);
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    if (parity == 0)
        return ps1;
    else
        return 0xFFFFFFFF - ps1 + 1;
}

uint32_t modinv_1c(uint32_t a, uint32_t p) {

    /* thanks to the folks at www.mersenneforum.org */
    // for use when it is known that p >> a, in which case
    // the first set of if/else blocks can be skipped
    uint32_t ps1, ps2, parity, dividend, divisor, rem, q, t;

    q = p / a;
    rem = p % a;
    dividend = a;
    divisor = rem;
    ps1 = q;
    ps2 = 1;
    parity = ~0;

    while (divisor > 1) {
        rem = dividend - divisor;
        t = rem - divisor;
        if (rem >= divisor) {
            q += ps1; rem = t; t -= divisor;
            if (rem >= divisor) {
                q += ps1; rem = t; t -= divisor;
                if (rem >= divisor) {
                    q += ps1; rem = t; t -= divisor;
                    if (rem >= divisor) {
                        q += ps1; rem = t; t -= divisor;
                        if (rem >= divisor) {
                            q += ps1; rem = t; t -= divisor;
                            if (rem >= divisor) {
                                q += ps1; rem = t; t -= divisor;
                                if (rem >= divisor) {
                                    q += ps1; rem = t; t -= divisor;
                                    if (rem >= divisor) {
                                        q += ps1; rem = t;
                                        if (rem >= divisor) {
                                            q = dividend / divisor;
                                            rem = dividend % divisor;
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    if (parity == 0)
        return ps1;
    else
        return p - ps1;
}

int is_mpz_prp(mpz_t n, int num_witnesses)
{
    int i = mpz_probab_prime_p(n, num_witnesses);
    return ((i == 1) || (i == 2)) && mpz_strongbpsw_prp(n);
}

int pull_twos(uint64_t* n, int* j, uint64_t p)
{
    int c = 0;

    while (!(*n & 1))
    {
        *n >>= 1;
        c = 1 - c;
    }
    if ((c * (p * p - 1) % 16) == 8)
        * j *= -1;
    return c;
}

int jacobi_1(uint64_t n, uint64_t p)
{
    //compute the jacobi symbol (n/p) for positive inputs
    //p must be odd
    //based on routine in Bressoud's book

    int j = 1;
    uint64_t t, nn = n;

    //return an error condition if p is even
    if (!(p & 1))
        return -2;

    nn = nn % p;

    //if p divides n then (n/p) = 0
    if (nn == 0)
        return 0;

    pull_twos(&nn, &j, p);
    while (nn > 1)
    {
        if (((nn - 1) * (p - 1)) % 8 == 4)
            j = -1 * j;
        t = nn;
        nn = p % nn;
        p = t;

        pull_twos(&nn, &j, p);
    }
    return j;
}

uint64_t spGCD(uint64_t x, uint64_t y)
{
    uint64_t a, b, c;
    a = x; b = y;
    while (b != 0)
    {
        c = a % b;
        a = b;
        b = c;
    }
    return a;
}

// straight from wikipedia.
uint64_t spBinGCD(uint64_t u, uint64_t v)
{
    // binary GCD for non-zero inputs.
    int shift;

    /* Let shift := lg K, where K is the greatest power of 2
    dividing both u and v. */
    for (shift = 0; ((u | v) & 1) == 0; ++shift) {
        u >>= 1;
        v >>= 1;
    }

    while ((u & 1) == 0)
        u >>= 1;

    /* From here on, u is always odd. */
    do {
        /* remove all factors of 2 in v -- they are not common */
        /*   note: v is not zero, so while will terminate */
        while ((v & 1) == 0)  /* Loop X */
            v >>= 1;

        /* Now u and v are both odd. Swap if necessary so u <= v,
        then set v = v - u (which is even). For bignums, the
        swapping is just pointer movement, and the subtraction
        can be done in-place. */
        if (u > v) {
            uint64_t t = v; v = u; u = t;
        }  // Swap u and v.
        v = v - u;                       // Here v >= u.
    } while (v != 0);

    /* restore common factors of 2 */
    return u << shift;
}

// assume u is odd
uint64_t spBinGCD_odd(uint64_t u, uint64_t v)
{
    /* From here on, u is always odd. */
    do {
        /* remove all factors of 2 in v -- they are not common */
        /*   note: v is not zero, so while will terminate */
        while ((v & 1) == 0)  /* Loop X */
            v >>= 1;

        /* Now u and v are both odd. Swap if necessary so u <= v,
        then set v = v - u (which is even). For bignums, the
        swapping is just pointer movement, and the subtraction
        can be done in-place. */
        if (u > v) {
            uint64_t t = v; v = u; u = t;
        }  // Swap u and v.
        v = v - u;                       // Here v >= u.
    } while (v != 0);

    /* restore common factors of 2 */
    return u;
}

// much faster version: assuming u is odd
uint64_t bingcd64(uint64_t u, uint64_t v)
{
#if 1
    if (u == 0) {
        return v;
    }
    if (v != 0) {
        int j = _trail_zcnt64(v);
        v = (uint64_t)(v >> j);
        while (1) {
            uint64_t tmp = u;
            uint64_t sub1 = (uint64_t)(v - tmp);
            uint64_t sub2 = (uint64_t)(tmp - v);
            if (tmp == v)
                break;
            u = (tmp >= v) ? v : tmp;
            v = (tmp >= v) ? sub2 : sub1;
            // For the line below, the standard way to write this algorithm
            // would have been to use _trail_zcnt64(v)  (instead of
            // _trail_zcnt64(sub1)).  However, as pointed out by
            // https://gmplib.org/manual/Binary-GCD, "in twos complement the
            // number of low zero bits on u-v is the same as v-u, so counting or
            // testing can begin on u-v without waiting for abs(u-v) to be
            // determined."  Hence we are able to use sub1 for the argument.
            // By removing the dependency on abs(u-v), the CPU can execute
            // _trail_zcnt64() at the same time as abs(u-v).
            j = _trail_zcnt64(sub1);
            v = (uint64_t)(v >> j);
        }
    }
    return u;
#else
// For reference, or if in the future we need to allow an even u,
// this version allows u to be even or odd.
    if (u == 0) {
        return v;
    }
    if (v != 0) {
        int i = _trail_zcnt64(u);
        int j = _trail_zcnt64(v);
        u = (uint64_t)(u >> i);
        v = (uint64_t)(v >> j);
        int k = (i < j) ? i : j;
        while (1) {
            uint64_t tmp = u;
            uint64_t sub1 = (uint64_t)(v - tmp);
            uint64_t sub2 = (uint64_t)(tmp - v);
            if (tmp == v)
                break;
            u = (tmp >= v) ? v : tmp;
            v = (tmp >= v) ? sub2 : sub1;
            // For the line below, the standard way to write this algorithm
            // would have been to use _trail_zcnt64(v)  (instead of
            // _trail_zcnt64(sub1)).  However, as pointed out by
            // https://gmplib.org/manual/Binary-GCD, "in twos complement the
            // number of low zero bits on u-v is the same as v-u, so counting or
            // testing can begin on u-v without waiting for abs(u-v) to be
            // determined."  Hence we are able to use sub1 for the argument.
            // By removing the dependency on abs(u-v), the CPU can execute
            // _trail_zcnt64() at the same time as abs(u-v).
            j = _trail_zcnt64(sub1);
            v = (uint64_t)(v >> j);
        }
        u = (uint64_t)(u << k);
    }
    return u;
#endif
}

uint64_t gcd64(uint64_t x, uint64_t y)
{
    uint64_t a, b, c;
    a = x; b = y;
    while (b != 0)
    {
        c = a % b;
        a = b;
        b = c;
    }
    return a;
}

void dblGCD(double x, double y, double* w)
{
    double a, b, c;
    a = x; b = y;
    while (b != 0)
    {
        c = a - b * (floor(a / b));
        a = b;
        b = c;
    }
    *w = a;
    return;
}

int llt(uint32_t exp, int VFLAG)
{
    mpz_t tmp, tmp2, n;
    clock_t start, stop;
    double t;
    uint32_t i, j, nchars;
    uint64_t d;

    mpz_init(tmp);
    mpz_set_ui(tmp, exp);
    if (!mpz_strongbpsw_prp(tmp))
    {
        mpz_clear(tmp);
        printf("exponent is not prime\n");
        return 0;
    }

    start = clock();
    mpz_init(n);
    mpz_setbit(n, exp);
    mpz_sub_ui(n, n, 1);
    //should vary the depth depending on the size of p
    for (i = 1; i < MIN(sqrt(exp) - 1, 1000000); i++)
    {
        d = 2 * i * (uint64_t)exp + 1;
        if (mpz_tdiv_ui(n, d) == 0)
        {
            mpz_clear(n);
            mpz_clear(tmp);
            if (VFLAG > 1)
                printf("2*%d*p+1 is a factor\n", i);
            return 0;
        }
    }
    if (VFLAG > 1)
        printf("trial division to %u bits is complete\n",
        (uint32_t)spBits(2 * 1000000 * exp + 1));

    mpz_init(tmp2);
    mpz_set_ui(tmp, 4);
    nchars = 0;
    //else do the ll test
    for (i = 0; i < exp - 2; i++)
    {
        mpz_mul(tmp, tmp, tmp);
        mpz_sub_ui(tmp, tmp, 2);
        /* Adapted from http://rosettacode.org/wiki/Lucas-Lehmer_test#GMP */
        /* mpz_tdiv_r(tmp, tmp, n); but more efficiently done given mod 2^p-1 */
        if (mpz_sgn(tmp) < 0) mpz_add(tmp, tmp, n);
        /* while (n > mp) { n = (n >> p) + (n & mp) } if (n==mp) n=0 */
        /* but in this case we can have at most one loop plus a carry */
        mpz_tdiv_r_2exp(tmp2, tmp, exp);
        mpz_tdiv_q_2exp(tmp, tmp, exp);
        mpz_add(tmp, tmp, tmp2);
        while (mpz_cmp(tmp, n) >= 0) mpz_sub(tmp, tmp, n);

        if (VFLAG > 1)
        {
            if ((i & 511) == 0)
            {
                for (j = 0; j < nchars; j++)
                    printf("\b");
                nchars = printf("llt iteration %d", i);
                fflush(stdout);
            }
        }
    }
    if (VFLAG > 1)
        printf("\n");

    if (VFLAG > 1)
    {
        stop = clock();
        t = (double)(stop - start) / (double)CLOCKS_PER_SEC;
        printf("elapsed time = %6.4f\n", t);
    }

    mpz_clear(n);
    mpz_clear(tmp2);

    if (mpz_cmp_ui(tmp, 0) == 0)
    {
        mpz_clear(tmp);
        return 1;
    }
    else
    {
        mpz_clear(tmp);
        return 0;
    }
}

void gordon(int bits, mpz_t p, gmp_randstate_t gmp_randstate)
{
    //find a random strong prime of size 'bits'
    //follows the Handbook of applied cryptography
    /*
    SUMMARY: a strong prime p is generated.
    1. Generate two large random primes s and t of roughly equal bitlength (see Note 4.54).
    2. Select an integer i0. Find the first prime in the sequence 2it + 1, for i = i0; i0 +
        1; i0 + 2; : : : (see Note 4.54). Denote this prime by r = 2it+ 1.
    3. Compute p0 = 2(sr-2 mod r)s - 1.
    4. Select an integer j0. Find the first prime in the sequence p0 +2jrs, for j = j0; j0 +
        1; j0 + 2; : : : (see Note 4.54). Denote this prime by p = p0 + 2jrs.
    5. Return(p).

  4.54 Note (implementing Gordon’s algorithm)
    (i) The primes s and t required in step 1 can be probable primes generated by Algorithm
    4.44. TheMiller-Rabin test (Algorithm 4.24) can be used to test each candidate
    for primality in steps 2 and 4, after ruling out candidates that are divisible by a small
    prime less than some boundB. See Note 4.45 for guidance on selecting B. Since the
    Miller-Rabin test is a probabilistic primality test, the output of this implementation
    of Gordon’s algorithm is a probable prime.
    (ii) By carefully choosing the sizes of primes s, t and parameters i0, j0, one can control
    the exact bitlength of the resulting prime p. Note that the bitlengths of r and s will
    be about half that of p, while the bitlength of t will be slightly less than that of r.
    */

    int i, j, s_len, n_words;
    mpz_t s, t, r, tmp, tmp2, p0;

    mpz_init(s);
    mpz_init(t);
    mpz_init(r);
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_init(p0);

    //need to check allocation of tmp vars.  how big do they get?

    //1. s and t should be about half the bitlength of p
    //random s of bitlength bits/2
    s_len = bits / 2 - 4;
    mpz_urandomb(s, gmp_randstate, s_len);
    mpz_setbit(s, s_len);

    //random t of bitlength bits/2
    mpz_urandomb(t, gmp_randstate, s_len);
    mpz_setbit(t, s_len);

    //2. Select an integer i0. Find the first prime in the sequence 2i(t) + 1, for i = i0; i0 +
    //1; i0 + 2; : : : (see Note 4.54). Denote this prime by r = 2i(t)+ 1.
    i = 1;
    mpz_mul_2exp(r, t, 1);
    mpz_add_ui(r, r, 1);
    while (!mpz_probab_prime_p(r, 1))
    {
        i++;
        mpz_mul_2exp(r, t, 1);
        mpz_mul_ui(r, r, i);
        mpz_add_ui(r, r, 1);
    }

    //3. Compute p0 = 2(sr-2 mod r)s - 1.
    //zMul(&s, &r, &tmp);
    //zShortSub(&tmp, 2, &p0);
    //zDiv(&p0, &r, &tmp, &tmp2);
    //zMul(&tmp2, &s, &p0);
    //zShiftLeft(&p0, &p0, 1);
    //zShortSub(&p0, 1, &p0);
    mpz_mul(tmp, s, r);
    mpz_sub_ui(p0, tmp, 2);
    mpz_tdiv_r(tmp2, p0, r);
    mpz_mul(p0, tmp2, s);
    mpz_mul_2exp(p0, p0, 1);
    mpz_sub_ui(p0, p0, 1);

    //4. Select an integer j0. Find the first prime in the sequence p0 +2jrs, for j = j0; j0 +
    //1; j0 + 2; : : : (see Note 4.54). Denote this prime by p = p0 + 2jrs.
    j = 1;
    //zMul(&r, &s, &tmp);
    //zShiftLeft(&tmp, &tmp, 1);
    //zAdd(&p0, &tmp, p);
    mpz_mul(tmp, r, s);
    mpz_mul_2exp(tmp, tmp, 1);
    mpz_add(p, p0, tmp);
    while (!mpz_probab_prime_p(p, 1))
    {
        j++;
        //zMul(&r, &s, &tmp);
        //zShiftLeft(&tmp, &tmp, 1);
        //zShortMul(&tmp, j, &tmp);
        //zAdd(&p0, &tmp, p);
        mpz_mul(tmp, r, s);
        mpz_mul_2exp(tmp, tmp, 1);
        mpz_mul_ui(tmp, tmp, j);
        mpz_add(p, p0, tmp);
    }

    mpz_clear(s);
    mpz_clear(t);
    mpz_clear(r);
    mpz_clear(tmp);
    mpz_clear(tmp2);
    mpz_clear(p0);
    return;
}

void build_RSA(int bits, mpz_t in, gmp_randstate_t gmp_randstate)
{
    int i;
    int words, subwords;
    mpz_t p, q;

    mpz_init(p);
    mpz_init(q);

    if (bits < 65)
    {
        printf("bitlength too small\n");
        return;
    }

    i = 0;
    while (mpz_sizeinbase(in, 2) != bits)
    {
        gordon(bits / 2, p, gmp_randstate);
        gordon(bits / 2, q, gmp_randstate);
        mpz_mul(in, p, q);
        i++;
    }

    mpz_clear(p);
    mpz_clear(q);
    return;
}


