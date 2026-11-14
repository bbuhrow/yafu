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

#ifndef ARITH_H
#define ARITH_H

#include <stdint.h>
#include <sys/types.h>
#include "gmp.h"
#include "gmp_u64.h"	// 64-bit-safe mpz_*_u64 wrappers (LLP64 Windows + GMP)

#define LIMB_BLKSZ 10	
#define MAX_DIGITS 100
#define MP_RADIX 4294967296.0
#define LN2		0.69314718055994530942

// types of numbers
#define PRIME 0
#define PRP 1
#define COMPOSITE 2
#define UNKNOWN 3

/* All one- and two-limb arithmetic (bit scans, carry chains, 64 x 64 -> 128
   multiply, 128 / 64 divide, single-word modular arithmetic, gcd, jacobi,
   128-bit arithmetic, ...) lives in the standalone library in arith/limb.
   This header keeps only the constants and the GMP-based utilities that are
   specific to YAFU. */
#include "limb1.h"
#include "limb2.h"

/********************* a few gmp-based utilities **********************/
double zlog(mpz_t x);
int llt(uint32_t exp, int vflag);
int gmp_base10(mpz_t x);
int is_mpz_prp(mpz_t n, int num_witnesses);
uint64_t mpz_get_64(mpz_t src);
void mpz_set_64(mpz_t dest, uint64_t src);
void build_RSA(int bits, mpz_t n, gmp_randstate_t gmp_randstate);
void gordon(int bits, mpz_t n, gmp_randstate_t gmp_randstate);
void fftmul(mpz_t c, mpz_t a, mpz_t b, int bits_per_word, int fftlen);

#endif
