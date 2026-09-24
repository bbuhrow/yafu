/*----------------------------------------------------------------------
limb1.c  --  single-limb (64-bit) arithmetic: implementation

Non-inline single-limb routines: the wide sp* helpers, single-word modular
exponentiation, modular inverse, the Jacobi symbol, Tonelli-Shanks, and the
scalar GCDs.  The hot inline kernels (addmod / submod / mulredc / sqrredc /
bingcd64 / ...) live in limb1.h.

Consolidated from arith.c.  The compiler / OS / ISA detection and the
_umul128 / _udiv128 / _addcarry_u64 / _subborrow_u64 polyfills now come from
mp_platform.h (via limb1.h) instead of the per-file #if ladders that used to
be re-derived (and had drifted) here and in monty.c / microecm.c.
----------------------------------------------------------------------*/

#include "limb1.h"
#include <math.h>       /* floor() for dblGCD */


/* ============================================================
   Wide single-word sp* operations.
   Two mutually-exclusive implementations:
     * GNU-style x86-64 inline asm  -- the fast path on gcc / clang / icc /
       mingw (selected by common.h's GCC_ASM64X)
     * portable, built on the mp_platform.h polyfills -- MSVC, or the
       ASM_ARITH_DEBUG A/B-comparison path
   The portable branch no longer defines its own _addcarry_u64 / rettype /
   COMPILER_MSVC block; it uses the single MSVC-convention polyfills from
   mp_platform.h.
   ============================================================ */
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

#else  /* portable path (MSVC, or ASM_ARITH_DEBUG) */

uint64_t spDivide(uint64_t *q, uint64_t *r, uint64_t u[2], uint64_t v)
{
    *q = _udiv128(u[1], u[0], v, r);
    return 0;
}

void spMultiply(uint64_t u, uint64_t v, uint64_t *product, uint64_t *carry)
{
    *product = _umul128(u, v, carry);
    return;
}

void spAdd(uint64_t u, uint64_t v, uint64_t *sum, uint64_t *carry)
{
    *carry = _addcarry_u64(0, u, v, sum);
    return;
}

void spAdd3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sum, uint64_t *carry)
{
    mp_carry_t c;
    *carry = _addcarry_u64(0, u, v, sum);
    c = _addcarry_u64((mp_carry_t)*carry, *sum, w, sum);
    *carry += c;
    return;
}

void spSub3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sub, uint64_t *borrow)
{
    mp_carry_t b;
    *borrow = _subborrow_u64(0, u, v, sub);
    b = _subborrow_u64(0, *sub, w, sub);
    *borrow += b;
    return;
}

void spSub(uint64_t u, uint64_t v, uint64_t *sub, uint64_t *borrow)
{
    *borrow = _subborrow_u64(0, u, v, sub);
    return;
}

void spMulAdd(uint64_t u, uint64_t v, uint64_t w, uint64_t t,
              uint64_t *lower, uint64_t *carry)
{
    uint64_t k, p;
    spMultiply(u, v, &p, carry);
    spAdd3(p, w, t, lower, &k);
    *carry += k;
    return;
}

void spMulMod(uint64_t u, uint64_t v, uint64_t m, uint64_t *w)
{
    uint64_t p[2];
    uint64_t q;
    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);
    return;
}

/* NOTE: spPRP2() and spModExp_asm() exist only in the asm branch above; there
   is no portable fallback for them, so a plain-MSVC build currently has no
   spPRP2.  If MSVC needs it, say the word and I'll add a portable spPRP2 built
   on spModExp (2^(p-1) mod p). */

#endif /* GCC_ASM64X && !ASM_ARITH_DEBUG */


/* ============================================================
   Small single-word helpers  (bit / digit counts, single divide)
   ============================================================ */
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

uint64_t u64div(uint64_t c, uint64_t n)
{
    uint64_t r;
    _udiv128(c, 0, n, &r);

    return r;
}


/* ============================================================
   Modular exponentiation & square root (Tonelli-Shanks)
   ============================================================ */
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


/* ============================================================
   Modular inverse (mersenneforum lineage)
   ============================================================ */
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


/* ============================================================
   Jacobi symbol
   ============================================================ */
static int pull_twos(uint64_t* n, int* j, uint64_t p)
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


/* ============================================================
   Scalar GCDs
   ============================================================ */
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

// bingcd64 is now the canonical inline definition in limb1.h
// (retired here; arith.c/monty.h/uecm/upm1 copies verified identical).

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

// my_clz64/my_ctz64/my_clz52/my_ctz52 no longer live here: the 64-bit pair
// became _lead_full_zcnt/_trail_full_zcnt in mp_platform.h, and the 52-bit
// pair became inline wrappers in limb1.h.  See the note in limb1.h.
