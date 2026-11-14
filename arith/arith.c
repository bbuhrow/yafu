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
	return mpz_get_u64(src);
}

void mpz_set_64(mpz_t dest, uint64_t src)
{
    // 64-bit-safe on LLP64 (MinGW + GMP), see arith/limb/gmp_u64.h
    mpz_set_u64(dest, src);
}

int gmp_base10(mpz_t x)
{
    mpz_t t;	//temp
    int g;		//guess: either correct or +1

    mpz_init(t);
    g = (int)mpz_sizeinbase(x, 10);
    mpz_set_ui(t, 10);
    mpz_pow_ui(t, t, (unsigned long)(g - 1));
    g = g - (mpz_cmp(t, x) > 0 ? 1 : 0);
    mpz_clear(t);
    return g;
}

// borrowed from jasonp... 
double zlog(mpz_t x) {

    uint32_t i = (uint32_t)mpz_size(x);

    switch (i) {

    case 0:
        return 0;
    case 1:
        return log((double)((uint32_t)mpz_get_ui(x)));
    case 2:

        return log((double)(mpz_getlimbn(x, 0)) +
            18446744073709551616.0 * (double)mpz_getlimbn(x, 1));

    default:

        return 64 * (i - 3) * LN2 +
            log((double)(mpz_getlimbn(x, i - 3)) + 18446744073709551616.0 * (
            ((double)mpz_getlimbn(x, i - 2) + 18446744073709551616.0 *
                (double)mpz_getlimbn(x, i - 1))));

    }
}
/* The single-word (sp*) helpers, modular exponentiation / inverse / sqrt,
   jacobi, the scalar gcds and the 128-bit divide that used to be defined
   here now live in arith/limb/limb1.c and limb2.c. */

int is_mpz_prp(mpz_t n, int num_witnesses)
{
    int i = mpz_probab_prime_p(n, num_witnesses);
    return ((i == 1) || (i == 2)) && mpz_strongbpsw_prp(n);
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
        if (mpz_tdiv_u64(n, d) == 0)
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
                nchars = (uint32_t)printf("llt iteration %d", i);
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

    int i, j, s_len;
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
    mpz_urandomb(s, gmp_randstate, (mp_bitcnt_t)s_len);
    mpz_setbit(s, (mp_bitcnt_t)s_len);

    //random t of bitlength bits/2
    mpz_urandomb(t, gmp_randstate, (mp_bitcnt_t)s_len);
    mpz_setbit(t, (mp_bitcnt_t)s_len);

    //2. Select an integer i0. Find the first prime in the sequence 2i(t) + 1, for i = i0; i0 +
    //1; i0 + 2; : : : (see Note 4.54). Denote this prime by r = 2i(t)+ 1.
    i = 1;
    mpz_mul_2exp(r, t, 1);
    mpz_add_ui(r, r, 1);
    while (!mpz_probab_prime_p(r, 1))
    {
        i++;
        mpz_mul_2exp(r, t, 1);
        mpz_mul_ui(r, r, (unsigned long)i);
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
        mpz_mul_ui(tmp, tmp, (unsigned long)j);
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


