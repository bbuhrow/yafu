/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#if defined(__CUDACC__) && !defined(CUDA_INTRINSICS_H)
#define CUDA_INTRINSICS_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef int int32;
typedef unsigned int uint32;
typedef unsigned long long uint64;
typedef long long int64;

/*------------------- Low-level functions ------------------------------*/

__device__ void
accum3(uint32 &a0, uint32 &a1, uint32 &a2,
	uint32 b0, uint32 b1) {

	asm("add.cc.u32 %0, %0, %3;   /* inline */   \n\t"
	    "addc.cc.u32 %1, %1, %4;   /* inline */   \n\t"
	    "addc.u32 %2, %2, %5;   /* inline */   \n\t"
		: "+r"(a0), "+r"(a1), "+r"(a2)
		: "r"(b0), "r"(b1), "r"(0) );
}

__device__ void
accum3_shift(uint32 &a0, uint32 &a1, uint32 &a2,
	uint32 b0, uint32 b1) {

	asm("add.cc.u32 %0, %1, %3;   /* inline */   \n\t"
	    "addc.cc.u32 %1, %2, %4;   /* inline */   \n\t"
	    "addc.u32 %2, %5, %5;   /* inline */   \n\t"
		: "=r"(a0), "+r"(a1), "+r"(a2)
		: "r"(b0), "r"(b1), "r"(0) );
}

__device__ void 
add3(uint32 a, uint32 b, uint32 c, uint32* sum, uint32* carry)
{
	uint32 s = c, cs = 0;

	asm("add.cc.u32 %0, %0, %2;   /* inline */   \n\t"
		"addc.u32 %1, 0, 0;   /* inline */   \n\t"
		"add.cc.u32 %0, %0, %3;   /* inline */   \n\t"
		"addc.u32 %1, %1, 0;   /* inline */   \n\t"
		: "+r"(s), "+r"(cs)
		: "r"(a), "r"(b));

	*sum = s;
	*carry = cs;

	return;
}

__device__ void 
add2(uint32 a, uint32 b, uint32* sum, uint32* carry)
{
	asm("add.cc.u32 %0, %2, %3;   /* inline */   \n\t"
		"addc.u32 %1, 0, 0;   /* inline */   \n\t"
		: "=r"(*sum), "=r"(*carry)
		: "r"(a), "r"(b));

	return;
}

__device__ void 
accum3lh(uint32 a, uint32 b, uint32 c, uint32 hi, uint32* sum, uint32* carry)
{
	// sum = a + lo + c
	// carry = sumcarry + hi
	uint32 s = c, cs = 0;

	asm("add.cc.u32 %0, %0, %2;   /* inline */   \n\t"
		"addc.u32 %1, 0, 0;   /* inline */   \n\t"
		"add.cc.u32 %0, %0, %3;   /* inline */   \n\t"
		"addc.u32 %1, %1, %4;   /* inline */   \n\t"
		: "+r"(s), "+r"(cs)
		: "r"(a), "r"(b), "r"(hi));

	*sum = s;
	*carry = cs;

	return;
}

__device__ uint32
innermul(uint32 a, uint32 b, uint32 prevhi, uint32 accum, uint32* carry)
{
	//uint32 lo = a * b;
	//uint32 hi = __umulhi(a, b);
	//accum3lh(accum, lo, prevhi, hi, &accum, carry);

	uint64 s = 0;

	asm("mad.wide.u32 %0, %1, %2, %3;    \n\t"
		"add.u64 %0, %0, %4;    \n\t"
		: "+l"(s)
		: "r"(a), "r"(b), "l"((uint64)accum), "l"((uint64)prevhi));

	*carry = (uint32)(s >> 32);
	return (uint32)s;
}

__device__ void 
accumlh(uint32 lo, uint32 hi, uint32* t, uint32* carry)
{
	asm("add.cc.u32 %0, %0, %2;   /* inline */   \n\t"
		"addc.u32 %1, %3, 0;   /* inline */   \n\t"
		: "+r"(*t), "=r"(*carry)
		: "r"(lo), "r"(hi));

	return;
}

/*----------------- Squaring ----------------------------------------*/

__device__ uint64 
wide_sqr32(uint32 a)
{
	uint32 a0, a1;

	asm("{ .reg .u64 %dprod; \n\t"
	    "mul.wide.u32 %dprod, %2, %2; \n\t"
	    "cvt.u32.u64 %0, %dprod;      \n\t"
	    "shr.u64 %dprod, %dprod, 32;  \n\t"
	    "cvt.u32.u64 %1, %dprod;      \n\t"
	    "}                   \n\t"
	    : "=r"(a0), "=r"(a1)
	    : "r"(a));

	return (uint64)a1 << 32 | a0;
}

/* -------------------- subtraction ------------------------*/
__device__ void
possub96(uint32* a, uint32* b, uint32* c)
{
	// subtract b from a when we know there will be no overflow
	uint32 r0, r1, r2;

	asm("{  \n\t"
		"sub.cc.u32 %0, %3, %6;        \n\t"
		"subc.cc.u32 %1, %4, %7;        \n\t"
		"subc.cc.u32 %2, %5, %8;        \n\t"
		"} \n\t"
		: "=r"(r0), "=r"(r1), "=r"(r2)
		: "r"(a[0]), "r"(a[1]), "r"(a[2]),
		"r"(b[0]), "r"(b[1]), "r"(b[2]));

	c[0] = r0;
	c[1] = r1;
	c[2] = r2;
	return;
}

/* -------------------- comparison ------------------------*/
__device__ void
cmp96(uint32* a, uint32* b, uint32* c)
{
	// subtract b from a when we know there will be no overflow
	uint32 r0, r1, r2;

	asm("{  \n\t"
		"sub.cc.u32 %0, %3, %6;        \n\t"
		"subc.cc.u32 %1, %4, %7;        \n\t"
		"subc.cc.u32 %2, %5, %8;        \n\t"
		"} \n\t"
		: "=r"(r0), "=r"(r1), "=r"(r2)
		: "r"(a[0]), "r"(a[1]), "r"(a[2]),
		"r"(b[0]), "r"(b[1]), "r"(b[2]));

	c[0] = r0;
	c[1] = r1;
	c[2] = r2;
	return;
}

/* -------------------- Modular subtraction ------------------------*/

__device__ uint32 
modsub32(uint32 a, uint32 b, uint32 p) 
{
	uint32 r;

	asm("{  \n\t"
	    ".reg .pred %pborrow;           \n\t"
	    ".reg .u32 %borrow;           \n\t"
	    "mov.b32 %borrow, 0;           \n\t"
	    "sub.cc.u32 %0, %1, %2;        \n\t"
	    "subc.u32 %borrow, %borrow, 0; \n\t"
	    "setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
	    "@%pborrow add.u32 %0, %0, %3; \n\t"
	    "} \n\t"
	    : "=r"(r) : "r"(a), "r"(b), "r"(p) );

	return r;
}

__device__ uint64 
modsub64(uint64 a, uint64 b, uint64 p) 
{
	uint32 r0, r1;
	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 p0 = (uint32)p;
	uint32 p1 = (uint32)(p >> 32);

	asm("{  \n\t"
	    ".reg .pred %pborrow;           \n\t"
	    ".reg .u32 %borrow;           \n\t"
	    "mov.b32 %borrow, 0;           \n\t"
	    "sub.cc.u32 %0, %2, %4;        \n\t"
	    "subc.cc.u32 %1, %3, %5;        \n\t"
	    "subc.u32 %borrow, %borrow, 0; \n\t"
	    "setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
	    "@%pborrow add.cc.u32 %0, %0, %6; \n\t"
	    "@%pborrow addc.u32 %1, %1, %7; \n\t"
	    "} \n\t"
	    : "=r"(r0), "=r"(r1)
	    : "r"(a0), "r"(a1), 
	      "r"(b0), "r"(b1), 
	      "r"(p0), "r"(p1));

	return ((uint64)r1 << 32) | (uint64)r0;
}

__device__ void
modsub96(uint32 *a, uint32 *b, uint32 *c, uint32 *p)
{
	uint32 r0, r1, r2;

	asm("{  \n\t"
		".reg .pred %pborrow;           \n\t"
		".reg .u32 %borrow;           \n\t"
		"mov.b32 %borrow, 0;           \n\t"
		"sub.cc.u32 %0, %3, %6;        \n\t"
		"subc.cc.u32 %1, %4, %7;        \n\t"
		"subc.cc.u32 %2, %5, %8;        \n\t"
		"subc.u32 %borrow, %borrow, 0; \n\t"
		"setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
		"@%pborrow add.cc.u32 %0, %0, %9; \n\t"
		"@%pborrow addc.cc.u32 %1, %1, %10; \n\t"
		"@%pborrow addc.cc.u32 %2, %2, %11; \n\t"
		"} \n\t"
		: "+r"(r0), "+r"(r1), "+r"(r2)
		: "r"(a[0]), "r"(a[1]), "r"(a[2]),
		"r"(b[0]), "r"(b[1]), "r"(b[2]),
		"r"(p[0]), "r"(p[1]), "r"(p[2]));

	c[0] = r0;
	c[1] = r1;
	c[2] = r2;
	return;
}

/* -------------------- Modular addition  ------------------------*/

__device__ uint32
modadd32(uint32 a, uint32 b, uint32 p)
{
	uint32 r;

	asm("{  \n\t"
		".reg .pred %pcarry;           \n\t"
		".reg .u32 %carry;           \n\t"
		"mov.b32 %carry, 0;           \n\t"
		"add.cc.u32 %0, %1, %2;        \n\t"
		"addc.u32 %carry, %carry, 0; \n\t"
		"setp.ne.u32 %pcarry, %carry, 0;  \n\t"
		"@%pcarry sub.u32 %0, %0, %3; \n\t"
		"} \n\t"
		: "=r"(r) : "r"(a), "r"(b), "r"(p));

	return r;
}

__device__ uint64
modadd64(uint64 a, uint64 b, uint64 p)
{
	uint32 r0, r1;
	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 p0 = (uint32)p;
	uint32 p1 = (uint32)(p >> 32);

	asm("{  \n\t"
		".reg .pred %pborrow;           \n\t"
		".reg .u32 %borrow;           \n\t"
		".reg .u32 %carry;           \n\t"
		".reg .u32 %sum0;           \n\t"
		".reg .u32 %sum1;           \n\t"
		"mov.b32 %carry, 0;           \n\t"
		"add.cc.u32 %0, %2, %4;        \n\t"
		"addc.cc.u32 %1, %3, %5;        \n\t"
		"addc.u32 %carry, %carry, 0; \n\t"
		"sub.cc.u32 %sum0, %0, %6;        \n\t"
		"subc.cc.u32 %sum1, %1, %7;        \n\t"
		"subc.u32 %borrow, %carry, 0; \n\t"
		"setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
		"@!%pborrow mov.b32 %0, %sum0; \n\t"
		"@!%pborrow mov.b32 %1, %sum1; \n\t"
		"} \n\t"
		: "=r"(r0), "=r"(r1)
		: "r"(a0), "r"(a1),	"r"(b0), "r"(b1), "r"(p0), "r"(p1));

	return ((uint64)r1 << 32) | r0;
}

__device__ void
modadd96(uint32* a, uint32* b, uint32* c, uint32* p)
{
	uint32 r0, r1, r2;

	asm("{  \n\t"
		".reg .pred %pborrow;           \n\t"
		".reg .u32 %borrow;           \n\t"
		".reg .u32 %carry;           \n\t"
		".reg .u32 %sum0;           \n\t"
		".reg .u32 %sum1;           \n\t"
		".reg .u32 %sum2;           \n\t"
		"mov.b32 %carry, 0;           \n\t"
		"add.cc.u32 %0, %3, %6;        \n\t"
		"addc.cc.u32 %1, %4, %7;        \n\t"
		"addc.cc.u32 %2, %5, %8;        \n\t"
		"addc.u32 %carry, %carry, 0; \n\t"
		"sub.cc.u32 %sum0, %0, %9;        \n\t"
		"subc.cc.u32 %sum1, %1, %10;        \n\t"
		"subc.cc.u32 %sum2, %2, %11;        \n\t"
		"subc.u32 %borrow, %carry, 0; \n\t"
		"setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
		"@!%pborrow mov.b32 %0, %sum0; \n\t"
		"@!%pborrow mov.b32 %1, %sum1; \n\t"
		"@!%pborrow mov.b32 %2, %sum2; \n\t"
		"} \n\t"
		: "=r"(r0), "=r"(r1), "=r"(r2)
		: "r"(a[0]), "r"(a[1]), "r"(a[2]),
		"r"(b[0]), "r"(b[1]), "r"(b[2]),
		"r"(p[0]), "r"(p[1]), "r"(p[2]));

	c[0] = r0;
	c[1] = r1;
	c[2] = r2;

	return;
}

__device__ void
modaddsub96(uint32* a, uint32* b, uint32* s, uint32* d, uint32* p)
{
	// compute (a - b) and (a + b) mod n
	uint32 s0, s1, s2;
	uint32 d0, d1, d2;

	asm("{  \n\t"
		".reg .pred %pborrow;           \n\t"
		".reg .u32 %borrow;           \n\t"
		".reg .u32 %carry;           \n\t"
		".reg .u32 %sum0;           \n\t"
		".reg .u32 %sum1;           \n\t"
		".reg .u32 %sum2;           \n\t"
		"mov.b32 %carry, 0;           \n\t"
		"add.cc.u32 %0, %6, %9;        \n\t"
		"addc.cc.u32 %1, %7, %10;        \n\t"
		"addc.cc.u32 %2, %8, %11;        \n\t"
		"addc.u32 %carry, %carry, 0; \n\t"
		"mov.b32 %borrow, 0;           \n\t"
		"sub.cc.u32 %3, %6, %9;        \n\t"
		"subc.cc.u32 %4, %7, %10;        \n\t"
		"subc.cc.u32 %5, %8, %11;        \n\t"
		"subc.u32 %borrow, %borrow, 0; \n\t"
		"setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
		"@%pborrow add.cc.u32 %3, %3, %12; \n\t"
		"@%pborrow addc.cc.u32 %4, %4, %13; \n\t"
		"@%pborrow addc.cc.u32 %5, %5, %14; \n\t"
		"sub.cc.u32 %sum0, %0, %12;        \n\t"
		"subc.cc.u32 %sum1, %1, %13;        \n\t"
		"subc.cc.u32 %sum2, %2, %14;        \n\t"
		"subc.u32 %borrow, %carry, 0; \n\t"
		"setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
		"@!%pborrow mov.b32 %0, %sum0; \n\t"
		"@!%pborrow mov.b32 %1, %sum1; \n\t"
		"@!%pborrow mov.b32 %2, %sum2; \n\t"
		"} \n\t"
		: "=r"(s0), "=r"(s1), "=r"(s2), "=r"(d0), "=r"(d1), "=r"(d2)
		: "r"(a[0]), "r"(a[1]), "r"(a[2]),
		"r"(b[0]), "r"(b[1]), "r"(b[2]),
		"r"(p[0]), "r"(p[1]), "r"(p[2]));

	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	d[0] = d0;
	d[1] = d1;
	d[2] = d2;

	return;
}

__device__ void 
modadd128(uint64* a, uint64* b, uint64* c, uint64* n)
{
	uint32 aa[4];
	uint32 bb[4];
	uint32 cc[4];
	uint32 nn[4];
	uint32 addcarry;
	uint32 subborrow;

	aa[0] = (uint32)a[0];
	aa[1] = (uint32)(a[0] >> 32);
	aa[2] = (uint32)a[1];
	aa[3] = (uint32)(a[1] >> 32);
	bb[0] = (uint32)b[0];
	bb[1] = (uint32)(b[0] >> 32);
	bb[2] = (uint32)b[1];
	bb[3] = (uint32)(b[1] >> 32);
	nn[0] = (uint32)n[0];
	nn[1] = (uint32)(n[0] >> 32);
	nn[2] = (uint32)n[1];
	nn[3] = (uint32)(n[1] >> 32);

	asm("add.cc.u32 %0, %0, %5;        \n\t"
		"addc.cc.u32 %1, %1, %6; \n\t"
		"addc.cc.u32 %2, %2, %7; \n\t"
		"addc.cc.u32 %3, %3, %8; \n\t"
		"addc.u32 %4, 0, 0; \n\t"
		: "+r"(aa[0]), "+r"(aa[1]), "+r"(aa[2]), "+r"(aa[3]), "=r"(addcarry)
		: "r"(bb[0]), "r"(bb[1]), "r"(bb[2]), "r"(bb[3]));

	cc[0] = aa[0];
	cc[1] = aa[1];
	cc[2] = aa[2];
	cc[3] = aa[3];

	asm("sub.cc.u32 %0, %0, %5;        \n\t"
		"subc.cc.u32 %1, %1, %6; \n\t"
		"subc.cc.u32 %2, %2, %7; \n\t"
		"subc.cc.u32 %3, %3, %8; \n\t"
		"subc.u32 %4, 0, 0; \n\t"
		: "+r"(aa[0]), "+r"(aa[1]), "+r"(aa[2]), "+r"(aa[3]), "=r"(subborrow)
		: "r"(nn[0]), "r"(nn[1]), "r"(nn[2]), "r"(nn[3]));

	if ((addcarry == 1) || (subborrow == 0))
	{
		c[0] = ((uint64)aa[1] << 32) | (uint64)aa[0];
		c[1] = ((uint64)aa[3] << 32) | (uint64)aa[2];
	}
	else
	{
		c[0] = ((uint64)cc[1] << 32) | (uint64)cc[0];
		c[1] = ((uint64)cc[3] << 32) | (uint64)cc[2];
	}

}

/*------------------------------- GCD --------------------------------*/
__device__  uint32
gcd32(uint32 x, uint32 y) {

	/* assumes x and y are odd and nonzero */

	uint32 u = x; 
	uint32 v = y;

	do {
		uint32 shift = 31 - __clz(v & -v);
		v = v >> shift;

		x = min(u, v);
		y = max(u, v);
		u = x;
		v = y - x;
	} while (v != 0);

	return u;
}

/*-------------------------- Modular inverse -------------------------*/

__device__ uint32 
modinv32(uint32 a, uint32 p) {

	uint32 ps1, ps2, dividend, divisor, rem, q, t;
	uint32 parity;

	q = 1; rem = a; dividend = p; divisor = a;
	ps1 = 1; ps2 = 0; parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend - q * divisor;
			q *= ps1;
		} } } } } } } } }

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

__device__ uint64
modinv64(uint64 a, uint64 p, uint64* likely_gcd) {

	uint64 ps1, ps2, dividend, divisor, rem, q, t;
	uint32 parity;

	q = 1; rem = a; dividend = p; divisor = a;
	ps1 = 1; ps2 = 0; parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend - q * divisor;
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}

	if (divisor == 1)
		dividend = divisor;
	*likely_gcd = dividend;
	
	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

__device__ void
modinv96(uint32 *a, uint32 *p, uint32* inv, uint32 *gcd) {

#ifdef __GNUC__

	__uint128_t ps1, ps2, dividend, divisor, rem, q, t;
	uint32 parity;

	q = 1;
	rem = (__uint128_t)a[2]; rem <<= 32;
	rem |= (__uint128_t)a[1];  rem <<= 32;
	rem |= (__uint128_t)a[0];
	dividend = (__uint128_t)p[2]; dividend <<= 32;
	dividend |= (__uint128_t)p[1];  dividend <<= 32;
	dividend |= (__uint128_t)p[0];
	divisor = (__uint128_t)a[2]; divisor <<= 32;
	divisor |= (__uint128_t)a[1];  divisor <<= 32;
	divisor |= (__uint128_t)a[0];
	ps1 = 1;
	ps2 = 0;
	parity = 0;

	while ( divisor > 1) {
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
											rem = dividend - q * divisor;
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

	if (divisor == 1)
		dividend = divisor;

	gcd[0] = (uint32_t)dividend; dividend >>= 32;
	gcd[1] = (uint32_t)dividend; dividend >>= 32;
	gcd[2] = (uint32_t)dividend;

	t = (__uint128_t)p[2]; t <<= 32;
	t |= (__uint128_t)p[1];  t <<= 32;
	t |= (__uint128_t)p[0];

	if (parity == 0)
		ps2 = ps1;
	else
		ps2 = t - ps1;

	inv[0] = (uint32_t)ps2; ps2 >>= 32;
	inv[1] = (uint32_t)ps2; ps2 >>= 32;
	inv[2] = (uint32_t)ps2;

#endif

}

/*------------------- Montgomery arithmetic --------------------------*/
__device__ uint32 
montmul32(uint32 a, uint32 b,
		uint32 n, uint32 w) {

	uint32 acc0, acc1, acc2 = 0;
	uint32 q;
	uint32 prod_lo, prod_hi;

	acc0 = a * b;
	acc1 = __umulhi(a, b);

	q = acc0 * w;

	prod_lo = q * n;
	prod_hi = __umulhi(q, n);

	accum3(acc0, acc1, acc2, prod_lo, prod_hi);

	if (acc2 || acc1 >= n)
		return acc1 - n;
	else
		return acc1;
}

__device__ uint64 
montmul64(uint64 a, uint64 b,
		uint64 n, uint32 w) {

	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 n0 = (uint32)n;
	uint32 n1 = (uint32)(n >> 32);
	uint32 acc0, acc1, acc2 = 0;
	uint32 q0, q1;
	uint32 prod_lo, prod_hi;
	uint64 r;

	acc0 = a0 * b0;
	acc1 = __umulhi(a0, b0);
	q0 = acc0 * w;
	prod_lo = q0 * n0;
	prod_hi = __umulhi(q0, n0);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);
	
	prod_lo = a0 * b1;
	prod_hi = __umulhi(a0, b1);
	accum3_shift(acc0, acc1, acc2, prod_lo, prod_hi);
	prod_lo = a1 * b0;
	prod_hi = __umulhi(a1, b0);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);
	prod_lo = q0 * n1;
	prod_hi = __umulhi(q0, n1);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);
	q1 = acc0 * w;
	prod_lo = q1 * n0;
	prod_hi = __umulhi(q1, n0);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);
	
	prod_lo = a1 * b1;
	prod_hi = __umulhi(a1, b1);
	accum3_shift(acc0, acc1, acc2, prod_lo, prod_hi);
	prod_lo = q1 * n1;
	prod_hi = __umulhi(q1, n1);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);

	r = (uint64)acc1 << 32 | acc0;
	if (acc2 || r >= n)
		return r - n;
	else
		return r;
}

__device__ void montmul96(uint32* a, uint32* b, uint32* c, uint32* n, uint32 w)
{
	// cios approach
	int i, j;
	uint32 t[5];

	for (i = 0; i < 3 + 2; i++)
	{
		t[i] = 0;
	}

	uint32 m;
	uint32 lo, hi;
	uint32 C, S;

	for (i = 0; i < 3; i++)
	{
		// T = T + a * b[i]
		// this works but is slower than a loop,
		// perhaps because the number of I/O causes
		// it to use more resources.  Note: on H200 card,
		// which has pretty good double precision performance.
		//asm("mad.lo.cc.u32 %0, %5, %8, %0;    \n\t"
		//	"madc.lo.cc.u32 %1, %6, %8, %1;    \n\t"
		//	"madc.lo.cc.u32 %2, %7, %8, %2;    \n\t"
		//	"addc.cc.u32 %3, %3, 0;    \n\t"
		//	"mad.hi.cc.u32 %1, %5, %8, %1;    \n\t"
		//	"madc.hi.cc.u32 %2, %6, %8, %2;    \n\t"
		//	"madc.hi.cc.u32 %3, %7, %8, %3;    \n\t"
		//	"addc.u32 %4, 0, 0;    \n\t"
		//	: "+r"(t[0]), "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "=r"(t[4])
		//	: "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(b[i]));

		// on 0th iteration can simplify the
		// accumulation a bit because C = 0;
		//lo = a[0] * b[i];
		//hi = __umulhi(a[0], b[i]);
		//accumlh(lo, hi, &t[0], &C);

		asm("mad.lo.cc.u32 %0, %2, %3, %0;    \n\t"
			"madc.hi.u32 %1, %2, %3, 0;    \n\t"
			: "+r"(t[0]), "=r"(C)
			: "r"(a[0]), "r"(b[i]));

		for (j = 1; j < 3; j++)
		{
			//lo = a[j] * b[i];
			//hi = __umulhi(a[j], b[i]);
			//accum3lh(t[j], lo, C, hi, &t[j], &C);
			uint64 s = 0;
			
			asm("mad.wide.u32 %0, %1, %2, %3;    \n\t"
				"add.u64 %0, %0, %4;    \n\t"
				: "+l"(s)
				: "r"(a[j]), "r"(b[i]), "l"((uint64)t[j]), "l"((uint64)C));
			
			C = (uint32)(s >> 32);
			t[j] = (uint32)s;
		}

		m = t[0] * w;

		add2(t[3], C, &S, &C);
		t[3] = S;
		t[4] = C;

		//lo = n[0] * m;
		//hi = __umulhi(n[0], m);
		//
		//accumlh(lo, hi, &t[0], &C);

		asm("mad.lo.cc.u32 %0, %2, %3, %0;    \n\t"
			"madc.hi.u32 %1, %2, %3, 0;    \n\t"
			: "+r"(t[0]), "=r"(C)
			: "r"(n[0]), "r"(m));

		for (j = 1; j < 3; j++)
		{
			//lo = n[j] * m;
			//hi = __umulhi(n[j], m);
			//accum3lh(t[j], lo, C, hi, &t[j - 1], &C);
			uint64 s = 0;

			asm("mad.wide.u32 %0, %1, %2, %3;    \n\t"
				"add.u64 %0, %0, %4;    \n\t"
				: "+l"(s)
				: "r"(n[j]), "r"(m), "l"((uint64)t[j]), "l"((uint64)C));

			C = (uint32)(s >> 32);
			t[j - 1] = (uint32)s;
		}

		add2(t[3], C, &t[2], &C);
		t[3] = t[4] + C;
		t[4] = 0;
	}

	uint32 cc[3] = { 0,0,0 };

	cc[0] = t[0];
	cc[1] = t[1];
	cc[2] = t[2];

	// non-balanced code paths, usually bad but this is faster... ?
	// AMM: only reduce if result > R
	if (t[3]) // || (carry == 0))
	{
		asm("sub.cc.u32 %0, %0, %3;        \n\t"
			"subc.cc.u32 %1, %1, %4; \n\t"
			"subc.u32 %2, %2, %5; \n\t"
			: "+r"(t[0]), "+r"(t[1]), "+r"(t[2])
			: "r"(n[0]), "r"(n[1]), "r"(n[2]));

		c[0] = t[0];
		c[1] = t[1];
		c[2] = t[2];
	}
	else
	{
		c[0] = cc[0];
		c[1] = cc[1];
		c[2] = cc[2];
	}

	return;
}

__device__ void montsqr96(uint32* a, uint32* c, uint32* n, uint32 w)
{
	// actual squaring enhancement using cios
	int i, j;
	uint32 t[5];

	for (i = 0; i < 3 + 2; i++)
	{
		t[i] = 0;
	}

	uint32 m;
	uint32 lo, hi;
	uint32 C, S;

	for (i = 0; i < 3; i++)
	{
		uint32 hicarry, p;

		//lo = a[i] * a[i];
		//hi = __umulhi(a[i], a[i]);
		//accumlh(lo, hi, &t[i], &C);
		asm("mad.lo.cc.u32 %0, %2, %2, %0;    \n\t"
			"madc.hi.u32 %1, %2, %2, 0;    \n\t"
			: "+r"(t[i]), "=r"(C)
			: "r"(a[i]));

		p = 0;
		for (j = i + 1; j < 3; j++)
		{
			lo = a[j] * a[i];
			hi = __umulhi(a[j], a[i]);

			hicarry = hi >> 31;
			hi = __funnelshift_l(lo, hi, 1);
			lo <<= 1;

			accum3lh(t[j], lo, C, hi, &t[j], &C);
			add2(C, p, &C, &p);
			p += hicarry;
		}

		m = t[0] * w;

		add2(t[3], C, &t[3], &C);
		t[4] = C + p;

		asm("mad.lo.cc.u32 %0, %2, %3, %0;    \n\t"
			"madc.hi.u32 %1, %2, %3, 0;    \n\t"
			: "+r"(t[0]), "=r"(C)
			: "r"(n[0]), "r"(m));

		for (j = 1; j < 3; j++)
		{
			//lo = n[j] * m;
			//hi = __umulhi(n[j], m);
			//accum3lh(t[j], lo, C, hi, &t[j - 1], &C);
			uint64 s = 0;

			asm("mad.wide.u32 %0, %1, %2, %3;    \n\t"
				"add.u64 %0, %0, %4;    \n\t"
				: "+l"(s)
				: "r"(n[j]), "r"(m), "l"((uint64)t[j]), "l"((uint64)C));

			C = (uint32)(s >> 32);
			t[j - 1] = (uint32)s;
		}

		add2(t[3], C, &t[2], &C);
		t[3] = t[4] + C;
		t[4] = 0;

	}

	uint32 cc[3] = { 0,0,0 };

	cc[0] = t[0];
	cc[1] = t[1];
	cc[2] = t[2];

	// non-balanced code paths, usually bad but this is faster... ?
	// AMM: only reduce if result > R
	if (t[3])
	{
		asm("sub.cc.u32 %0, %0, %3;        \n\t"
			"subc.cc.u32 %1, %1, %4; \n\t"
			"subc.u32 %2, %2, %5; \n\t"
			: "+r"(t[0]), "+r"(t[1]), "+r"(t[2])
			: "r"(n[0]), "r"(n[1]), "r"(n[2]));

		c[0] = t[0];
		c[1] = t[1];
		c[2] = t[2];
	}
	else
	{
		c[0] = cc[0];
		c[1] = cc[1];
		c[2] = cc[2];
	}

	return;
}

__device__ void montmul128(uint32* a, uint32* b, uint32* c, uint32* n, uint32 w)
{
	// cios approach
	int i, j;
	uint32 t[6];

	for (i = 0; i < 4 + 2; i++)
	{
		t[i] = 0;
	}

	uint32 m;
	uint32 lo, hi;
	uint32 C, S;

	for (i = 0; i < 4; i++)
	{
		C = 0;
		for (j = 0; j < 4; j++)
		{
			lo = a[j] * b[i];
			hi = __umulhi(a[j], b[i]);
			accum3lh(t[j], lo, C, hi, &t[j], &C);
		}

		m = t[0] * w;

		add2(t[4], C, &S, &C);
		t[4] = S;
		t[5] = C;
		
		lo = n[0] * m;
		hi = __umulhi(n[0], m);

		accumlh(lo, hi, &t[0], &C);

		for (j = 1; j < 4; j++)
		{
			lo = n[j] * m;
			hi = __umulhi(n[j], m);
			accum3lh(t[j], lo, C, hi, &t[j - 1], &C);
		}

		add2(t[4], C, &t[3], &C);
		t[4] = t[5] + C;
		t[5] = 0;
	}

	uint32 cc[4] = { 0,0,0,0 };
	uint32 carry = 0;

	cc[0] = t[0];
	cc[1] = t[1];
	cc[2] = t[2];
	cc[3] = t[3];

	// non-balanced code paths, usually bad but this is faster... ?
	// AMM: only reduce if result > R
	if (t[4])
	{
		asm("sub.cc.u32 %0, %0, %5;        \n\t"
			"subc.cc.u32 %1, %1, %6; \n\t"
			"subc.cc.u32 %2, %2, %7; \n\t"
			"subc.cc.u32 %3, %3, %8; \n\t"
			"subc.u32 %4, 0, 0; \n\t"
			: "+r"(t[0]), "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "=r"(carry)
			: "r"(n[0]), "r"(n[1]), "r"(n[2]), "r"(n[3]));

		c[0] = t[0];
		c[1] = t[1];
		c[2] = t[2];
		c[3] = t[3];
	}
	else
	{
		c[0] = cc[0];
		c[1] = cc[1];
		c[2] = cc[2];
		c[3] = cc[3];
	}

	return;
}

__device__ void montsqr128(uint32* a, uint32* c, uint32* n, uint32 w)
{
	// actual squaring enhancement using cios

	int i, j;
	uint32 t[6];

	for (i = 0; i < 4 + 2; i++)
	{
		t[i] = 0;
	}

	uint32 m;
	uint32 lo, hi;
	uint32 C, S;

	for (i = 0; i < 4; i++)
	{
		lo = a[i] * a[i];
		hi = __umulhi(a[i], a[i]);
		accumlh(lo, hi, &t[i], &C);
		
		uint32 hicarry, p;
		p = 0;

		for (j = i + 1; j < 4; j++)
		{
			lo = a[j] * a[i];
			hi = __umulhi(a[j], a[i]);

			hicarry = hi >> 31;
			hi = __funnelshift_l(lo, hi, 1);
			lo <<= 1;

			accum3lh(t[j], lo, C, hi, &t[j], &C);
			
			//C += p;
			//p = hicarry;
			add2(C, p, &C, &p);
			p += hicarry;
		}

		m = t[0] * w;

		add2(t[4], C, &t[4], &C);
		t[5] = C + p;

		lo = n[0] * m;
		hi = __umulhi(n[0], m);
		accumlh(lo, hi, &t[0], &C);

		for (j = 1; j < 4; j++)
		{
			lo = n[j] * m;
			hi = __umulhi(n[j], m);
			accum3lh(t[j], lo, C, hi, &t[j - 1], &C);
		}

		add2(t[4], C, &t[3], &C);
		t[4] = t[5] + C;
		t[5] = 0;

	}

	uint32 cc[4] = { 0,0,0,0 };
	uint32 carry = 0;

	cc[0] = t[0];
	cc[1] = t[1];
	cc[2] = t[2];
	cc[3] = t[3];

	// non-balanced code paths, usually bad but this is faster... ?
	// AMM: only reduce if result > R
	if (t[4])
	{
		asm("sub.cc.u32 %0, %0, %5;        \n\t"
			"subc.cc.u32 %1, %1, %6; \n\t"
			"subc.cc.u32 %2, %2, %7; \n\t"
			"subc.cc.u32 %3, %3, %8; \n\t"
			"subc.u32 %4, 0, 0; \n\t"
			: "+r"(t[0]), "+r"(t[1]), "+r"(t[2]), "+r"(t[3]), "=r"(carry)
			: "r"(n[0]), "r"(n[1]), "r"(n[2]), "r"(n[3]));

		c[0] = t[0];
		c[1] = t[1];
		c[2] = t[2];
		c[3] = t[3];
	}
	else
	{
		c[0] = cc[0];
		c[1] = cc[1];
		c[2] = cc[2];
		c[3] = cc[3];
	}

	return;
}

#define _BIGWORDS 4

__device__ void bigmontmul(uint32* a, uint32* b, uint32* c, uint32* n, uint32 w)
{
	// cios approach
	int i, j;
	uint32 t[_BIGWORDS + 2];

	for (i = 0; i < _BIGWORDS + 2; i++)
	{
		t[i] = 0;
	}

	uint32 m;
	uint32 lo, hi;
	uint32 C, S;

	for (i = 0; i < _BIGWORDS; i++)
	{
		C = 0;
		for (j = 0; j < _BIGWORDS; j++)
		{
			lo = a[j] * b[i];
			hi = __umulhi(a[j], b[i]);
			accum3lh(t[j], lo, C, hi, &t[j], &C);
		}

		m = t[0] * w;

		add2(t[_BIGWORDS], C, &S, &C);
		t[_BIGWORDS] = S;
		t[_BIGWORDS + 1] = C;

		lo = n[0] * m;
		hi = __umulhi(n[0], m);

		accumlh(lo, hi, &t[0], &C);

		for (j = 1; j < _BIGWORDS; j++)
		{
			lo = n[j] * m;
			hi = __umulhi(n[j], m);
			accum3lh(t[j], lo, C, hi, &t[j - 1], &C);
		}

		add2(t[_BIGWORDS], C, &t[_BIGWORDS - 1], &C);
		t[_BIGWORDS] = t[_BIGWORDS + 1] + C;
		t[_BIGWORDS + 1] = 0;
	}

	uint32 cc[_BIGWORDS];
	uint32 carry = 0;

	for (i = 0; i < _BIGWORDS; i++)
	{
		cc[i] = t[i];
	}

	for (i = 0; i < _BIGWORDS; i += 4)
	{
		asm("sub.cc.u32 %4, 0, %4;        \n\t"
			"subc.cc.u32 %0, %0, %5;        \n\t"
			"subc.cc.u32 %1, %1, %6; \n\t"
			"subc.cc.u32 %2, %2, %7; \n\t"
			"subc.cc.u32 %3, %3, %8; \n\t"
			"subc.u32 %4, 0, 0; \n\t"
			: "+r"(t[i + 0]), "+r"(t[i + 1]), "+r"(t[i + 2]), "+r"(t[i + 3]), "+r"(carry)
			: "r"(n[i + 0]), "r"(n[i + 1]), "r"(n[i + 2]), "r"(n[i + 3]));
	}

	// non-balanced code paths, usually bad but this is faster... ?
	// AMM: only reduce if result > R
	if (t[_BIGWORDS])
	{
		for (i = 0; i < _BIGWORDS; i++)
		{
			c[i] = t[i];
		}
	}
	else
	{
		for (i = 0; i < _BIGWORDS; i++)
		{
			c[i] = cc[i];
		}
	}

	return;
}

__device__ void bigmontsqr(uint32* a, uint32* c, uint32* n, uint32 w)
{
	// actual squaring enhancement using cios

	int i, j;
	uint32 t[_BIGWORDS + 2];

	for (i = 0; i < _BIGWORDS + 2; i++)
	{
		t[i] = 0;
	}

	uint32 m;
	uint32 lo, hi;
	uint32 C, S;

	for (i = 0; i < _BIGWORDS; i++)
	{
		lo = a[i] * a[i];
		hi = __umulhi(a[i], a[i]);
		accumlh(lo, hi, &t[i], &C);

		uint32 hicarry, p;
		p = 0;

		for (j = i + 1; j < _BIGWORDS; j++)
		{
			lo = a[j] * a[i];
			hi = __umulhi(a[j], a[i]);

			hicarry = hi >> 31;
			hi = __funnelshift_l(lo, hi, 1);
			lo <<= 1;

			accum3lh(t[j], lo, C, hi, &t[j], &C);
			//C += p;
			//p = hicarry;
			add2(C, p, &C, &p);
			p += hicarry;
		}

		m = t[0] * w;

		add2(t[_BIGWORDS], C, &t[_BIGWORDS], &C);
		t[_BIGWORDS + 1] = C + p;

		lo = n[0] * m;
		hi = __umulhi(n[0], m);
		accumlh(lo, hi, &t[0], &C);

		for (j = 1; j < _BIGWORDS; j++)
		{
			lo = n[j] * m;
			hi = __umulhi(n[j], m);
			accum3lh(t[j], lo, C, hi, &t[j - 1], &C);
		}

		add2(t[_BIGWORDS], C, &t[_BIGWORDS - 1], &C);
		t[_BIGWORDS] = t[_BIGWORDS + 1] + C;
		t[_BIGWORDS + 1] = 0;

	}

	uint32 cc[_BIGWORDS];
	uint32 carry = 0;

	for (i = 0; i < _BIGWORDS; i++)
	{
		cc[i] = t[i];
	}

	for (i = 0; i < _BIGWORDS; i += 4)
	{
		asm("sub.cc.u32 %4, 0, %4;        \n\t"
			"subc.cc.u32 %0, %0, %5;        \n\t"
			"subc.cc.u32 %1, %1, %6; \n\t"
			"subc.cc.u32 %2, %2, %7; \n\t"
			"subc.cc.u32 %3, %3, %8; \n\t"
			"subc.u32 %4, 0, 0; \n\t"
			: "+r"(t[i + 0]), "+r"(t[i + 1]), "+r"(t[i + 2]), "+r"(t[i + 3]), "+r"(carry)
			: "r"(n[i + 0]), "r"(n[i + 1]), "r"(n[i + 2]), "r"(n[i + 3]));
	}

	// AMM: only reduce if result > R
	if (t[_BIGWORDS])
	{
		for (i = 0; i < _BIGWORDS; i++)
			c[i] = t[i];
	}
	else
	{
		for (i = 0; i < _BIGWORDS; i++)
			c[i] = cc[i];
	}

	return;
}

/*------------------ Initializing Montgomery arithmetic -----------------*/
__device__ uint32 
montmul32_w(uint32 n) {

	uint32 res = 2 + n;
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	return res * (2 + n * res);
}

__device__ uint32 
montmul32_r(uint32 n) {

	uint32 r0 = ((uint64)1 << 63) % n;
	uint32 r1;

	r1 = r0 + r0;

	if (r1 < r0)
		r1 -= n;

	return modsub32(r1, n, n);
}

__device__ uint64 
montmul64_r(uint64 n, uint32 w) {

	uint32 shift;
	uint32 i;
	uint64 shifted_n;
	uint64 res;

	shift = __clzll(n);
	shifted_n = n << shift;
	res = -shifted_n;

	for (i = 64 - shift; i < 72; i++) {
		if (res >> 63)
			res = res + res - shifted_n;
		else
			res = res + res;

		if (res >= shifted_n)
			res -= shifted_n;
	}

	res = res >> shift;
	res = montmul64(res, res, n, w);
	res = montmul64(res, res, n, w);
	return montmul64(res, res, n, w);
}

#ifdef __cplusplus
}
#endif

#endif /* defined(__CUDACC__) && !defined(CUDA_INTRINSICS_H) */

