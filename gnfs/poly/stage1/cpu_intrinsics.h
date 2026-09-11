/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: cpu_intrinsics.h 1095 2026-06-15 11:25:15Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef CPU_INTRINSICS_H
#define CPU_INTRINSICS_H

#include <mp.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if defined(_MSC_VER) && defined(_WIN64)

	#include <intrin.h>
	#pragma intrinsic(__umul64)
	#pragma intrinsic(__umul128)
	#pragma intrinsic(__subborrow_u32)
	#pragma intrinsic(__subborrow_u64)

	#define PROD32(hi, lo, a, b) 		\
	{					\
		lo = _umul64(a, b, &(hi));	\
	}

	#define PROD64(hi, lo, a, b) 		\
	{					\
		lo = _umul128(a, b, &(hi));	\
	}

	#define accum32(a0, a1, a2, b0, b1)	\
	{					\
		uint8 cy = _subborrow_u32(0, a0, b0, &(a0));	\
		cy = _subborrow_u32(cy, a1, b1, &(a1));		\
		_subborrow_u32(cy, a2, 0, &(a2));		\
	}

	#define accum64(a0, a1, a2, b0, b1)	\
	{					\
		uint8 cy = _subborrow_u64(0, a0, b0, &(a0));	\
		cy = _subborrow_u64(cy, a1, b1, &(a1));		\
		_subborrow_u64(cy, a2, 0, &(a2));		\
	}

#elif defined(GCC_ASM64X)

	#define PROD32(hi, lo, a, b) \
		asm("mull %2  \n\t"      \
		:"=d"(hi), "=a"(lo)  \
		:"%rm"(a), "1"(b)    \
		:"cc")

	#define PROD64(hi, lo, a, b) \
		asm("mulq %2  \n\t"      \
		:"=d"(hi), "=a"(lo)  \
		:"%rm"(a), "1"(b)    \
		:"cc")

	#define accum32(a0, a1, a2, b0, b1)	\
		asm("subl %3, %0  \n\t"		\
		    "sbbl %4, %1  \n\t"		\
		    "sbbl $0, %2  \n\t"		\
		: "+r"(a0), "+r"(a1), "+r"(a2)	\
		: "g"(b0), "g"(b1) : "cc")

	#define accum64(a0, a1, a2, b0, b1)	\
		asm("subq %3, %0  \n\t"		\
		    "sbbq %4, %1  \n\t"		\
		    "sbbq $0, %2  \n\t"		\
		: "+r"(a0), "+r"(a1), "+r"(a2)	\
		: "g"(b0), "g"(b1) : "cc")

#else
	#error "unsupported compiler combination"
#endif

/*------------------- Montgomery arithmetic --------------------------*/
static INLINE uint32 
montmul32(uint32 a, uint32 b, uint32 n, uint32 w) 
{
	uint32 acc0, acc1, acc2 = 0;
	uint32 q;
	uint32 prod_lo, prod_hi;

	PROD32(acc1, acc0, a, b);
	q = acc0 * w;
	PROD32(prod_hi, prod_lo, q, n);
	accum32(acc0, acc1, acc2, prod_lo, prod_hi);

	if (acc2)
		return acc1 + n;
	else
		return acc1;
}

static INLINE uint64 
montmul64(uint64 a, uint64 b, uint64 n, uint64 w) 
{
	uint64 acc0, acc1, acc2 = 0;
	uint64 q;
	uint64 prod_lo, prod_hi;

	PROD64(acc1, acc0, a, b);
	q = acc0 * w;
	PROD64(prod_hi, prod_lo, q, n);
	accum64(acc0, acc1, acc2, prod_lo, prod_hi);

	if (acc2)
		return acc1 + n;
	else
		return acc1;
}

static INLINE uint64
mod64_32(uint64 num, uint32 t, uint32 n, uint32 w) {

	uint32 x, rem, hi, tmp;

	x = (uint32)num;
	tmp = x * w;
	PROD32(hi, rem, tmp, n);
	x = (uint64)(num >> 32);
	tmp = x - hi;
	tmp = tmp * w;
	if (hi > x)
		tmp++;
	PROD32(hi, rem, tmp, n);

	rem = montmul32(hi, t, n, w);

	if (rem == 0)
		return 0;
	else
		return n - rem;
}

static INLINE uint32
mod128_32(uint128 num, uint32 t, uint32 n, uint32 w) 
{
	uint32 x, tmp, lo, hi;

	x = num.w[0];
	tmp = x * w;
	PROD32(hi, lo, tmp, n);
	x = num.w[1];
	tmp = x - hi;
	tmp = tmp * w;
	if (hi > x)
		tmp++;
	PROD32(hi, lo, tmp, n);
	x = num.w[2];
	tmp = x - hi;
	tmp = tmp * w;
	if (hi > x)
		tmp++;
	PROD32(hi, lo, tmp, n);
	x = num.w[3];
	tmp = x - hi;
	tmp = tmp * w;
	if (hi > x)
		tmp++;
	PROD32(hi, lo, tmp, n);
	hi = montmul32(hi, t, n, w);

	if (hi == 0)
		return 0;
	else
		return n - hi;
}

static INLINE uint64
mod128_64(uint128 num, uint64 t, uint64 n, uint64 w) {

	uint64 x, rem, hi, tmp;

	x = (uint64)num.w[1] << 32 | num.w[0];
	tmp = x * w;
	PROD64(hi, rem, tmp, n);
	x = (uint64)num.w[3] << 32 | num.w[2];
	tmp = x - hi;
	tmp = tmp * w;
	if (hi > x)
		tmp++;
	PROD64(hi, rem, tmp, n);

	rem = montmul64(hi, t, n, w);

	if (rem == 0)
		return 0;
	else
		return n - rem;
}

/*------------------ Initializing Montgomery arithmetic -----------------*/
static INLINE uint64 
montmul64_w(uint64 n) {

	uint64 v = 3 * n ^ 2;
	v = v * ((uint64)2 - n * v);
	v = v * ((uint64)2 - n * v);
	v = v * ((uint64)2 - n * v);
	v = v * ((uint64)2 - n * v);
	return v;
}

static INLINE uint32 
montmul32_w(uint32 n) {

	uint32 v = 3 * n ^ 2;
	v = v * ((uint32)2 - n * v);
	v = v * ((uint32)2 - n * v);
	v = v * ((uint32)2 - n * v);
	return v;
}

#define SPECIALQ_BATCH_SIZE 50


#ifdef __cplusplus
}
#endif

#endif /* !CPU_INTRINSICS_H */

