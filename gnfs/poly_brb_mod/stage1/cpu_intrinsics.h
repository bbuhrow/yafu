/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: cpu_intrinsics.h 680 2011-11-24 16:25:27Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef CPU_INTRINSICS_H
#define CPU_INTRINSICS_H

#include <mp.h>

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef _MSC_VER

#include <intrin.h>
#pragma intrinsic(__emulu)

#define PROD32(hi, lo, a, b)		\
	{	uint64 __t = __emulu(a,b);	\
		hi = (uint32)(__t >> 32);	\
		lo = (uint32)(__t); }

#elif defined(GCC_ASM32X)

#define PROD32(hi, lo, a, b) \
	asm("mull %2  \n\t"      \
	    :"=d"(hi), "=a"(lo)  \
	    :"%rm"(a), "1"(b)    \
	    :"cc")

#else

#define PROD32(hi, lo, a, b) \
	{ uint64 t = (uint64)(a) * (b); \
	  hi = (uint32)(t >> 32);	\
	  lo = (uint32)(t); }

#endif

/*------------------- Montgomery arithmetic --------------------------*/
static INLINE uint64 
montmul64(uint64 a, uint64 b,
		uint64 n, uint32 w) 
{
	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 n0 = (uint32)n;
	uint32 n1 = (uint32)(n >> 32);
	uint32 acc0, acc1, acc2, nmult;
	uint32 prod_lo, prod_hi;
	uint64 prod;

	PROD32(prod_hi, prod_lo, a0, b0);
	acc0 = prod_lo;

	prod = (uint64)prod_hi;
	PROD32(prod_hi, prod_lo, a1, b0);
	prod += ((uint64)prod_hi << 32 | prod_lo);
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32);

	nmult = acc0 * w;

	PROD32(prod_hi, prod_lo, nmult, n0);
	prod = acc0 + ((uint64)prod_hi << 32 | prod_lo);
	prod = prod >> 32;

	PROD32(prod_hi, prod_lo, nmult, n1);
	prod += (uint64)acc1 + ((uint64)prod_hi << 32 | prod_lo);
	acc0 = (uint32)prod;
	prod = (prod >> 32) + (uint64)acc2;
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32);

	PROD32(prod_hi, prod_lo, a0, b1);
	prod = (uint64)acc0 + ((uint64)prod_hi << 32 | prod_lo);
	acc0 = (uint32)prod;
	prod = prod >> 32;

	PROD32(prod_hi, prod_lo, a1, b1);
	prod += (uint64)acc1 + ((uint64)prod_hi << 32 | prod_lo);
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32) + acc2;

	nmult = acc0 * w;

	PROD32(prod_hi, prod_lo, nmult, n0);
	prod = acc0 + ((uint64)prod_hi << 32 | prod_lo);
	prod = prod >> 32;

	PROD32(prod_hi, prod_lo, nmult, n1);
	prod += acc1 + ((uint64)prod_hi << 32 | prod_lo);
	acc0 = (uint32)prod;
	prod = (prod >> 32) + (uint64)acc2;
	acc1 = (uint32)prod;
	acc2 = (uint32)(prod >> 32);

	prod = (uint64)acc1 << 32 | acc0;
	if (acc2 || prod >= n)
		return prod - n;
	else
		return prod;
}

/*------------------ Initializing Montgomery arithmetic -----------------*/
static INLINE uint32 
montmul32_w(uint32 n) {

	uint32 res = 2 + n;
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	return res * (2 + n * res);
}

static INLINE uint64 
montmul64_r(uint64 n) 
{
	mp_t num, den, rem;

	num.val[0] = 0;
	num.val[1] = 0;
	num.val[2] = 0;
	num.val[3] = 0;
	num.val[4] = 1;
	num.nwords = 5;

	den.val[0] = (uint32)n;
	den.val[1] = (uint32)(n >> 32);
	den.nwords = 2;
	if (den.val[1] == 0)
		den.nwords = 1;

	mp_mod(&num, &den, &rem);
	return (uint64)rem.val[1] << 32 | rem.val[0];
}

/*------------------------- Miscellaneous -------------------------------*/

#define HOST_BATCH_SIZE 8192

#define INVERT_BATCH_SIZE 512

#ifdef __cplusplus
}
#endif

#endif /* !CPU_INTRINSICS_H */

