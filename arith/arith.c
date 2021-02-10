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
#include "ytools.h"
#include "common.h"

uint64 mpz_get_64(mpz_t src)
{

	uint64 out = mpz_getlimbn(src, 0);
#if GMP_LIMB_BITS == 32
	if (mpz_size(src) >= 2)
		out |= ((uint64)mpz_getlimbn(src, 1) << 32ULL);
#endif

	return out;

}

void mpz_set_64(mpz_t dest, uint64 src)
{
    //printf("%d %d %d\n", GMP_LIMB_BITS, sizeof(mp_limb_t), mp_bits_per_limb);
#if GMP_LIMB_BITS == 64
    mpz_set_ui(dest, src);
#else
	/* mpz_import is terribly slow */
    mpz_set_ui(dest, src >> 32);
    mpz_mul_2exp(dest, dest, 32);
    mpz_set_ui(dest, src & 0xffffffff)
#endif

}

int ndigits_1(fp_digit n)
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

double rint(double x)
{
	 double i, r = modf(x, &i);
	 if (r < 0.0) {
		 r += 1.0; 
		 i -= 1.0;
	 }
	 return (r > 0.5 || (r == 0.5 && ((int)i & 1)) ? i + 1.0 : i);
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

#if GMP_LIMB_BITS == 32
    uint32 i = mpz_size(x);

    switch (i) {
    case 0:
        return 0;
    case 1:
        return log((double)((uint32)mpz_get_ui(x)));
    case 2:

        return log((double)(mpz_getlimbn(x, 0)) +
            MP_RADIX * mpz_getlimbn(x, 1));

    default:

        return 32 * (i - 3) * LN2 +
            log((double)(mpz_getlimbn(x, i - 3)) + MP_RADIX * (
            ((double)mpz_getlimbn(x, i - 2) + MP_RADIX *
                mpz_getlimbn(x, i - 1))));

    }

#else

    uint32 i = mpz_size(x);

    switch (i) {

    case 0:
        return 0;
    case 1:
        return log((double)((uint32)mpz_get_ui(x)));
    case 2:

        return log((double)(mpz_getlimbn(x, 0)) +
            18446744073709551616.0 * mpz_getlimbn(x, 1));

    default:

        return 64 * (i - 3) * LN2 +
            log((double)(mpz_getlimbn(x, i - 3)) + 18446744073709551616.0 * (
            ((double)mpz_getlimbn(x, i - 2) + 18446744073709551616.0 *
                mpz_getlimbn(x, i - 1))));

    }


#endif

}

fp_digit spBits(fp_digit n)
{
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }

    return i;
}

int bits64(uint64 n)
{
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }
    return i;
}

#if defined(MSC_ASM32A) && !defined(ASM_ARITH_DEBUG)

void spAdd(fp_digit u, fp_digit v, fp_digit* sum, fp_digit* carry)
{
    fp_digit s, c;

    s = v;
    c = 0;
    ASM_M{
        mov eax, u
        add s, eax
        adc c, 0
    }

    *sum = s;
    *carry = c;

    return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sum, fp_digit* carry)
{
    fp_digit s, c;

    s = v;
    c = 0;

    ASM_M{
        mov eax, u
        add eax, w
        adc c, 0
        add s, eax
        adc c, 0
    }

    *sum = s;
    *carry = c;

    return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sub, fp_digit* borrow)
{
    fp_digit s, b;

    s = v;
    b = 0;

    ASM_M{
        mov eax, u
        sub eax, s
        adc b, 0
        sub eax, w
        adc b, 0
        mov s, eax
    }

    *sub = s;
    *borrow = b;

    return;
}

void spSub(fp_digit u, fp_digit v, fp_digit* sub, fp_digit* borrow)
{
    fp_digit s, b;

    s = v;
    b = 0;

    ASM_M{
        mov eax, u
        sub eax, s
        adc b, 0
        mov s, eax
    }

    *sub = s;
    *borrow = b;

    return;
}

fp_digit spDivide(fp_digit* q, fp_digit* r, fp_digit u[2], fp_digit v)
{
    fp_digit tmpq, tmpr;
    tmpr = *r = u[1];
    tmpq = *q = u[0];

    ASM_M{
        mov eax, tmpq
        mov edx, tmpr
        div v
        mov tmpr, edx
        mov tmpq, eax
    }

    *q = tmpq;
    *r = tmpr;

    return 0;
}

void spMultiply(fp_digit u, fp_digit v, fp_digit* product, fp_digit* carry)
{
    ASM_M{
        mov eax, u
        mul v
        mov v, eax
        mov u, edx
    }

    *product = v;
    *carry = u;

    return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit* lower, fp_digit* carry)
{
    //u*v + (w+t)  used a lot in multiplication
    //fp_word uu;
    //uu = (fp_word)u * (fp_word)v + (fp_word)w + (fp_word)t;
    //*lower = (fp_digit)uu;
    //*carry = (fp_digit)(uu >> BITS_PER_DIGIT);
    fp_digit tmp;

    spMultiply(u, v, lower, carry);
    spAdd3(*lower, w, t, lower, &tmp);
    spAdd(tmp, *carry, carry, &tmp);
    return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit* w)
{
    fp_digit p[2];
    fp_digit q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}

#elif defined(GCC_ASM32A) && !defined(ASM_ARITH_DEBUG)

void spAdd(fp_digit u, fp_digit v, fp_digit* sum, fp_digit* carry)
{
    fp_word s, c;

    s = v;
    c = 0;

    ASM_G("movl %2, %%eax		\n\t"
        "addl %%eax, %3		\n\t"
        "adcl $0, %4		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(u), "0"(s), "1"(c)
        : "eax", "memory", "cc");

    *sum = s;
    *carry = c;

    return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sum, fp_digit* carry)
{
    fp_word s, c;

    s = v;
    c = 0;

    ASM_G("movl %2, %%eax		\n\t"
        "addl %3, %%eax		\n\t"
        "adcl $0, %5		\n\t"
        "addl %%eax, %4		\n\t"
        "adcl $0, %5		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(u), "r"(w), "0"(s), "1"(c)
        : "eax", "memory", "cc");

    *sum = s;
    *carry = c;

    return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sub, fp_digit* borrow)
{
    fp_word s, b;

    s = v;
    b = 0;

    ASM_G("movl %2, %%eax		\n\t"
        "subl %4, %%eax		\n\t"
        "adcl $0, %5		\n\t"
        "subl %3, %%eax		\n\t"
        "adcl $0, %5		\n\t"
        "movl %%eax, %4		\n\t"
        : "=r"(s), "=r"(b)
        : "r"(u), "r"(w), "0"(s), "1"(b)
        : "eax", "memory", "cc");

    *sub = s;
    *borrow = b;

    return;
}

void spSub(fp_digit u, fp_digit v, fp_digit* sub, fp_digit* borrow)
{
    fp_word s, b;

    s = v;
    b = 0;

    ASM_G("movl %2, %%eax		\n\t"
        "subl %3, %%eax		\n\t"
        "adcl $0, %4		\n\t"
        "movl %%eax, %3		\n\t"
        : "=r"(s), "=r"(b)
        : "r"(u), "0"(s), "1"(b)
        : "eax", "memory", "cc");

    *sub = s;
    *borrow = b;

    return;
}


fp_digit spDivide(fp_digit* q, fp_digit* r, fp_digit u[2], fp_digit v)
{
    *r = u[1];
    *q = u[0];
    ASM_G("divl %4"
        : "=a"(*q), "=d"(*r)
        : "1"(*r), "0"(*q), "r"(v));

    return 0;
}

void spMultiply(fp_digit u, fp_digit v, fp_digit* product, fp_digit* carry)
{
    *product = v;
    *carry = u;

    ASM_G("movl %2, %%eax	\n\t"
        "mull %3	\n\t"
        "movl %%eax, %0		\n\t"
        "movl %%edx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product)
        : "eax", "%rdx", "cc");

    return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit* lower, fp_digit* carry)
{
    fp_digit k, p;
    spMultiply(u, v, &p, carry);
    spAdd3(p, w, t, lower, &k);
    *carry += k;
    return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit* w)
{
    fp_digit p[2];
    fp_digit q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}

#elif defined(GCC_ASM64X) && !defined(ASM_ARITH_DEBUG)

void spAdd(fp_digit u, fp_digit v, fp_digit* sum, fp_digit* carry)
{
    //fp_word s,c;
    uint64 s, c;

    s = v;
    c = 0;

    ASM_G("movq %2, %%rax		\n\t"
        "addq %%rax, %3		\n\t"
        "adcq $0, %4		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(u), "0"(s), "1"(c)
        : "rax", "memory", "cc");

    *sum = s;
    *carry = c;

    return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sum, fp_digit* carry)
{
    //fp_word s,c;
    uint64 s, c;

    s = v;
    c = 0;

    ASM_G("movq %2, %%rax		\n\t"
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

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sub, fp_digit* borrow)
{
    //fp_word s,b;
    uint64 s, b;

    s = v;
    b = 0;

    ASM_G("movq %2, %%rax		\n\t"
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

void spSub(fp_digit u, fp_digit v, fp_digit* sub, fp_digit* borrow)
{
    //fp_word s,b;
    uint64 s, b;

    s = v;
    b = 0;

    ASM_G("movq %2, %%rax		\n\t"
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

fp_digit spDivide(fp_digit* q, fp_digit* r, fp_digit u[2], fp_digit v)
{
    *r = u[1];
    *q = u[0];
    ASM_G("divq %4"
        : "=a"(*q), "=d"(*r)
        : "1"(*r), "0"(*q), "r"(v));

    return 0;
}

void spMultiply(fp_digit u, fp_digit v, fp_digit* product, fp_digit* carry)
{
    *product = v;
    *carry = u;

    ASM_G("movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product)
        : "rax", "rdx", "cc");

    return;
}

uint64 spPRP2(uint64 p)
{
    // do a base-2 prp test on the input, where p is greater than 2^32
    // i.e., compute 2^(p-1) % p.
    // since p is more than 32 bits we can do the accumulation division 
    // free for the first 5 iterations.  may not be much, but it's something.
    uint64 result;

    ASM_G(
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

uint64 spModExp_asm(uint64 b, uint64 e, uint64 m)
{
    uint64 result;

    ASM_G(
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

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit* lower, fp_digit* carry)
{
    fp_digit k, p;
    spMultiply(u, v, &p, carry);
    spAdd3(p, w, t, lower, &k);
    *carry += k;
    return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit* w)
{
    fp_digit p[2];
    fp_digit q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}

#elif BITS_PER_DIGIT == 32

fp_digit spDivide(fp_digit * q, fp_digit * r, fp_digit u[2], fp_digit v)
{
    fp_word uu, qq;
    uu = (fp_word)u[1] << BITS_PER_DIGIT | u[0];
    qq = (uu / (fp_word)v);
    *r = (fp_digit)(uu % (fp_word)v);
    *q = (fp_digit)qq;
    return (fp_digit)(qq >> BITS_PER_DIGIT);
}

void spMultiply(fp_digit u, fp_digit v, fp_digit* product, fp_digit* carry)
{
    fp_word uu;
    uu = (fp_word)u * (fp_word)v;
    *product = (fp_digit)uu;
    *carry = (fp_digit)(uu >> BITS_PER_DIGIT);
    return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit* lower, fp_digit* carry)
{
    //u*v + (w+t)  used a lot in multiplication
    //fp_word uu;
    //uu = (fp_word)u * (fp_word)v + (fp_word)w + (fp_word)t;
    //*lower = (fp_digit)uu;
    //*carry = (fp_digit)(uu >> BITS_PER_DIGIT);
    fp_digit tmp;

    spMultiply(u, v, lower, carry);
    spAdd3(*lower, w, t, lower, &tmp);
    spAdd(tmp, *carry, carry, &tmp);
    return;
}

void spAdd(fp_digit u, fp_digit v, fp_digit* sum, fp_digit* carry)
{
    fp_word w;
    w = (fp_word)u + (fp_word)v;
    *sum = (fp_digit)w;
    *carry = (fp_digit)(w >> BITS_PER_DIGIT);
    return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sum, fp_digit* carry)
{
    fp_word uu;
    uu = (fp_word)u + (fp_word)v + (fp_word)w;
    *sum = (fp_digit)uu;
    *carry = (fp_digit)(uu >> BITS_PER_DIGIT);
    return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sub, fp_digit* borrow)
{
    fp_word uu;
    uu = (fp_word)u - (fp_word)v - (fp_word)w;
    *sub = (fp_digit)uu;
    *borrow = (fp_digit)((uu >> BITS_PER_DIGIT) && MAX_DIGIT);
    return;
}

void spSub(fp_digit u, fp_digit v, fp_digit* sub, fp_digit* borrow)
{
    fp_word uu;
    uu = (fp_word)u - (fp_word)v;
    *sub = (fp_digit)uu;
    *borrow = (fp_digit)((uu >> BITS_PER_DIGIT) && MAX_DIGIT);
    return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit* w)
{
    fp_digit p[2];
    fp_digit q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}

#elif BITS_PER_DIGIT == 64

// TBD: needs low level assembly routines similar to mod_64

fp_digit spDivide(fp_digit * q, fp_digit * r, fp_digit u[2], fp_digit v)
{
    *q = _udiv128(u[1], u[0], v, r);
    return 0;
}

void spMultiply(fp_digit u, fp_digit v, fp_digit* product, fp_digit* carry)
{
    *product = _umul128(u, v, carry);
    return;
}

void spAdd(fp_digit u, fp_digit v, fp_digit* sum, fp_digit* carry)
{
    *carry = _addcarry_u64(0, u, v, sum);
    return;
}

void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sum, fp_digit* carry)
{
    unsigned char c;
    *carry = _addcarry_u64(0, u, v, sum);
    c = _addcarry_u64(*carry, *sum, w, sum);
    *carry += c;
    return;
}

void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit* sub, fp_digit* borrow)
{
    unsigned char b;
    *borrow = _subborrow_u64(0, u, v, sub);
    b = _subborrow_u64(0, *sub, w, sub);
    *borrow += b;
    return;
}

void spSub(fp_digit u, fp_digit v, fp_digit* sub, fp_digit* borrow)
{
    *borrow = _subborrow_u64(0, u, v, sub);
    return;
}

void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit* lower, fp_digit* carry)
{
    fp_digit k, p;
    spMultiply(u, v, &p, carry);
    spAdd3(p, w, t, lower, &k);
    *carry += k;
    return;
}

void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit* w)
{
    fp_digit p[2];
    fp_digit q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}

#endif


void spModExp(fp_digit a, fp_digit b, fp_digit m, fp_digit* u)
{
    //computes a^b mod m = u using the binary method
    //see, for instance, the handbook of applied cryptography
    fp_digit n, bb, aa, t, prod[2];

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

void ShanksTonelli_1(fp_digit a, fp_digit p, fp_digit* sq)
{
    //a is a quadratic residue mod p
    //p is an odd prime
    //find x where x^2 == a mod p
    //we assume p will always fit into an fp_digit, therefore x will as well.
    //see paper by Ezra Brown
    fp_digit x = 0, b = 0, g = 0, n = 0, s = 0, r = 0, e = 0, b2m = 0, tmp = 0;
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

uint32 modinv_1(uint32 a, uint32 p) {

    /* thanks to the folks at www.mersenneforum.org */

    uint32 ps1, ps2, parity, dividend, divisor, rem, q, t;


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

uint32 modinv_1b(uint32 a, uint32 p) {

    /* thanks to the folks at www.mersenneforum.org */

    /* modification: p is fixed at 2^32.  a is only valid if odd */

    uint64 dividend = (uint64)0x1 << 32;
    uint32 ps1, ps2, parity, divisor, rem, q, t;

    q = 1;
    rem = a;
    //dividend = p;
    divisor = a;
    ps1 = 1;
    ps2 = 0;
    parity = 0;

    while (divisor > 1) {
        rem = (uint32)(dividend - (uint64)divisor);
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
                                            q = (uint32)(dividend / (uint64)divisor);
                                            rem = (uint32)(dividend % (uint64)divisor);
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

uint32 modinv_1c(uint32 a, uint32 p) {

    /* thanks to the folks at www.mersenneforum.org */
    // for use when it is known that p >> a, in which case
    // the first set of if/else blocks can be skipped
    uint32 ps1, ps2, parity, dividend, divisor, rem, q, t;

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
    return mpz_probab_prime_p(n, num_witnesses) &&
        mpz_strongbpsw_prp(n);
}

int pull_twos(fp_digit* n, int* j, fp_digit p)
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


int jacobi_1(fp_digit n, fp_digit p)
{
    //compute the jacobi symbol (n/p) for positive inputs
    //p must be odd
    //based on routine in Bressoud's book

    int j = 1;
    fp_digit t, nn = n;

    //return an error condition if p is even
    if (!(p & 1))
        return -2;

    nn = nn % p;

    //if p divides n then (n/p) = 0
    if (nn == 0)
        return 0;

    /*
    if (nn<0)
    {
        nn *= -1;
        if (p%4 == 3)
            j = -1;
    }
    */

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

fp_digit spGCD(fp_digit x, fp_digit y)
{
    fp_digit a, b, c;
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
uint64 spBinGCD(uint64 u, uint64 v)
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
            uint64 t = v; v = u; u = t;
        }  // Swap u and v.
        v = v - u;                       // Here v >= u.
    } while (v != 0);

    /* restore common factors of 2 */
    return u << shift;
}

// assume u is odd
uint64 spBinGCD_odd(uint64 u, uint64 v)
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
            uint64 t = v; v = u; u = t;
        }  // Swap u and v.
        v = v - u;                       // Here v >= u.
    } while (v != 0);

    /* restore common factors of 2 */
    return u;
}

// much faster version: assuming x is odd
uint64 bingcd64(uint64 x, uint64 y)
{
    if (y) {
        y >>= _trail_zcnt64(y);
        while (x != y)
            if (x < y)
                y -= x, y >>= _trail_zcnt64(y);
            else
                x -= y, x >>= _trail_zcnt64(x);
    }
    return x;
}

uint64 gcd64(uint64 x, uint64 y)
{
    uint64 a, b, c;
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

int llt(uint32 exp, int VFLAG)
{
    mpz_t tmp, tmp2, n;
    clock_t start, stop;
    double t;
    uint32 i, j, nchars;
    fp_digit d;

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
        d = 2 * i * (fp_digit)exp + 1;
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
        (uint32)spBits(2 * 1000000 * exp + 1));

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


