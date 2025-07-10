// MIT License
// 
// Copyright (c) 2024 Ben Buhrow
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// MIT License
// 
// Copyright (c) 2024 Pierre
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <stdint.h>
#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include "tinyprp.h"
#include "monty.h"			// functions for 52- 64- 104- and 128-bit Montgomery arithmetic
#include "ytools.h"
#include "mpz_aprcl.h"

//#define GMP_CHECK

#ifdef GMP_CHECK
#include "gmp.h"
#endif

int fermat_prp_64x1(uint64_t n)
{
	uint64_t rho = multiplicative_inverse(n);	// pos variant.  neg would be (uint64_t)0ull - multiplicative_inverse(n);
	uint64_t unityval = ((uint64_t)0 - n) % n;   // unityval == R  (mod n)
	uint64_t result = unityval;
	uint64_t e = (n - 1); // / 2;

	// penultimate-hi-bit mask
#if defined(USE_AVX2) || defined(USE_AVX512F)
	// technically need to check the ABM flag, but I don't
	// have that in place anywhere yet.  AVX2 is generally equivalent.

#if defined( __INTEL_COMPILER) || defined(_MSC_VER)

	uint64_t m = 1ULL << (62 - __lzcnt64(n));   // set a mask at the leading bit - 2

#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER)

	uint64_t m = 1ULL << (62 - __builtin_clzll(n));

#endif

#else
	// these builtin functions will have an efficient implementation
	// for the current processor architecture.
#if defined( __INTEL_COMPILER) || defined(_MSC_VER)

	uint32_t pos;
	if (_BitScanReverse64(&pos, n))
		return pos;
	else
		return 64;

	uint64_t m = 1ULL << (62 - pos);   // set a mask at the leading bit - 2

#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER)

	uint64_t m = 1ULL << (62 - __builtin_clzll(n));

#endif

#endif

	result = addmod(result, result, n);

	while (m > 0)
	{
		result = sqrredc_pos(result, n, rho);
		if (e & m) result = addmod(result, result, n);
		m >>= 1;
	}

	// - Fermat primality check:
	//   (2^(n-1) == 1) mod n
	return (int)(result == unityval);
	// 
	// - Euler's criterion 2^(n>>1) == legendre_symbol(2,n) (https://en.wikipedia.org/wiki/Euler%27s_criterion)
	// - Euler primality check:
	//   (2^(n>>1) == 1) mod n
	//   (2^(n>>1) == n-1) mod n
	uint64_t legendre = ((n >> 1) ^ (n >> 2)) & 1;	// shortcut calculation of legendre symbol
	uint64_t m1 = submod(n, unityval, n);
	return ((result == (legendre ? m1 : unityval)));
}

int MR_2sprp_64x1(uint64_t n)
{
	uint64_t rho = multiplicative_inverse(n);	// pos variant.  neg would be (uint64_t)0ull - multiplicative_inverse(n);
	uint64_t one = ((uint64_t)0 - n) % n;		// one == R  (mod n)
	uint64_t result = one;
	uint64_t e = (n - 1);

	// compute d and s
	uint64_t s = my_ctz64(e);
	uint64_t d = e >> s;

	// penultimate-hi-bit mask of d
#if defined(USE_AVX2) || defined(USE_AVX512F)
	// technically need to check the ABM flag, but I don't
	// have that in place anywhere yet.  AVX2 is generally equivalent.

#if defined( __INTEL_COMPILER) || defined(_MSC_VER)

	uint64_t m = 1ULL << (62 - __lzcnt64(d));   // set a mask at the leading bit - 2

#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER)

	uint64_t m = 1ULL << (62 - __builtin_clzll(d));

#endif

#else
	// these builtin functions will have an efficient implementation
	// for the current processor architecture.
#if defined( __INTEL_COMPILER) || defined(_MSC_VER)

	uint32_t pos;
	if (_BitScanReverse64(&pos, d))
		return pos;
	else
		return 64;

	uint64_t m = 1ULL << (62 - pos);   // set a mask at the leading bit - 2

#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER)

	uint64_t m = 1ULL << (62 - __builtin_clzll(d));

#endif

#endif

	// compute b^d using RL binary exponentiation.  RL because then
	// exponent 1-bits are adds instead of multiplies.
	// we know the first bit is set and the first squaring is of unity,
	// so the first iteration is easier (and hence the penultimate mask bit).
	result = addmod(result, result, n);

	while (m > 0)
	{
		result = sqrredc_pos(result, n, rho);
		if (d & m) result = addmod(result, result, n);
		m >>= 1;
	}

	// check current result == 1
	if (result == one) return 1;

	// now compute b^(2^r*d) for 0 <= r < s, 
	// and check for congruence to -1 as we go.
	uint64_t minus1 = n - one;

	while (s > 1)
	{
		if (result == minus1) return 1;
		result = sqrredc_pos(result, n, rho);
		if (result == one) return 0;
		s--;
	}

	if (result == minus1) return 1;
	else return 0;
}

int modexp_128x1b(uint64_t* n, uint64_t* e, uint64_t b)
{
	// compute b^e % n and check if equal to 1 or n-1 (-1)

	uint64_t rho = multiplicative_inverse(n[0]);
	uint64_t one[2]; // = ((uint64_t)0 - n) % n;  // unityval == R  (mod n)
	uint64_t r[2];
	uint64_t d[2];
	uint64_t n1[2];
	uint64_t base[2];

#ifdef POS_VARIANT

#else
	rho = 0ULL - rho;
#endif

	uint128_t n128 = ((uint128_t)n[1] << 64) | n[0];
	uint128_t unity128 = (uint128_t)0 - n128;
	unity128 = unity128 % n128;
	one[1] = (uint64_t)(unity128 >> 64);
	one[0] = (uint64_t)unity128;

	submod128(n, one, n1, n);

	// base to monty-rep
	r[1] = 0;
	r[0] = 0;

	d[0] = one[0];
	d[1] = one[1];

	while (b > 0)
	{
		if (b & 1)
			addmod128(r, d, r, n);
		addmod128(d, d, d, n);
		b >>= 1;
	}

	base[0] = r[0];
	base[1] = r[1];

	//printf("n: %016lx%016lx\n", n[1], n[0]);
	//printf("e: %016lx%016lx\n", e[1], e[0]);
	//printf("1: %016lx%016lx\n", one[1], one[0]);
	//printf("b: %016lx%016lx (%lu)\n", base[1], base[0], b);

	// RL-64x2
	r[1] = one[1];
	r[0] = one[0];
	d[0] = e[0];
	d[1] = e[1];

	int i = 0;
	while (d[0] > 0)
	{
		if (d[0] & 1)
			mulmod128n(r, base, r, n, rho);

		d[0] >>= 1;
		i++;

		sqrmod128n(base, base, n, rho);
	}

	for (; (i < 64) && (d[1] > 0); i++)
	{
		sqrmod128n(base, base, n, rho);
	}

	while (d[1] > 0)
	{
		if (d[1] & 1)
			mulmod128n(r, base, r, n, rho);

		d[1] >>= 1;
		i++;

		sqrmod128n(base, base, n, rho);
	}

	// AMM possibly needs a final correction by n
	chkmod128(r, n);
	chkmod128(r, n);

	//printf("r: %016lx%016lx\n", r[1], r[0]);

	if ((r[0] == one[0]) && (r[1] == one[1])) return 1;
	if ((r[0] == n1[0]) && (r[1] == n1[1])) return -1;
	return 0;
}

int fermat_prp_128x1(uint64_t* n)
{
	// assumes has no small factors.  
	// assumes n has two 64-bit words, n0 and n1.
	// do a base-2 fermat prp test using LR binexp.

	uint64_t rho = multiplicative_inverse(n[0]);
	uint64_t unityval[2]; // = ((uint64_t)0 - n) % n;  // unityval == R  (mod n)
	uint64_t m;
	uint64_t r[2];
	uint64_t e[2];
	uint64_t one[2] = { 1ULL, 0ULL };
	uint64_t zero[2] = { 0ULL, 0ULL };

#ifdef POS_VARIANT

#else
	rho = 0ULL - rho;
#endif

	uint128_t n128 = ((uint128_t)n[1] << 64) | n[0];
	uint128_t unity128 = (uint128_t)0 - n128;
	unity128 = unity128 % n128;
	unityval[1] = (uint64_t)(unity128 >> 64);
	unityval[0] = (uint64_t)unity128;

	e[1] = n[1];
	e[0] = n[0] - 1;	// n odd: won't overflow

	r[1] = unityval[1];
	r[0] = unityval[0];

	int lzcnt = my_clz64(e[1]);

#ifndef POS_VARIANT
	int protect = (lzcnt < 3) ? 1 : 0;
#endif

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	dblmod128(r, n);

#ifndef POS_VARIANT
	if (protect)
	{
		m = (e[1] <= 1) ? 0 : 1ULL << (62 - lzcnt);   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);
			chkmod128(r, n);

			if (e[1] & m) dblmod128(r, n);
			m >>= 1;
		}

		m = (e[1] >= 1) ? 1ULL << 63 : 1ULL << (62 - my_clz64(e[0]));   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);
			chkmod128(r, n);

			if (e[0] & m) dblmod128(r, n);
			m >>= 1;
		}
	}
	else
#endif
	{
		m = (e[1] <= 1) ? 0 : 1ULL << (62 - lzcnt);   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);

			if (e[1] & m) dblmod128(r, n);
			m >>= 1;
		}

		m = (e[1] >= 1) ? 1ULL << 63 : 1ULL << (62 - my_clz64(e[0]));   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);

			if (e[0] & m) dblmod128(r, n);
			m >>= 1;
		}
	}

	chkmod128(r, n);
	return ((r[0] == unityval[0]) && (r[1] == unityval[1]));
}

int MR_2sprp_128x1(uint64_t *n)
{
	// assumes has no small factors.  
	// assumes n has two 64-bit words, n0 and n1.
	// do a base-2 Miller-Rabin sprp test.

	uint64_t rho = multiplicative_inverse(n[0]);
	uint64_t m;
	uint64_t r[2];
	uint64_t e[2];
	uint64_t one[2] = { 1ULL, 0ULL };
	uint64_t zero[2] = { 0ULL, 0ULL };

#ifdef POS_VARIANT

#else
	rho = 0ULL - rho;
#endif

	uint128_t n128 = ((uint128_t)n[1] << 64) | n[0];
	uint128_t unity128 = (uint128_t)0 - n128;
	unity128 = unity128 % n128;
	one[1] = (uint64_t)(unity128 >> 64);
	one[0] = (uint64_t)unity128;

	e[1] = n[1];
	e[0] = n[0] - 1;	// n odd: won't overflow

	r[1] = one[1];
	r[0] = one[0];

	// compute d and s
	uint64_t s = my_ctz128(e[0], e[1]);
	uint64_t d[2];
	d[0] = e[0];
	d[1] = e[1];

	if (s < 64)
	{
		d[0] >>= s;
		d[0] |= (d[1] << (64 - s));
		d[1] >>= s;
	}
	else
	{
		d[0] = d[1] >> (s - 64);
		d[1] = 0;
	}

	int lzcnt = my_clz128(d[0], d[1]);

#ifndef POS_VARIANT
	int protect = ((lzcnt-s) < 3) ? 1 : 0;
#endif

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	dblmod128(r, n);

#ifndef POS_VARIANT
	if (protect)
	{
		m = (d[1] <= 1) ? 0 : 1ULL << (62 - lzcnt);   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);
			chkmod128(r, n);

			if (d[1] & m) dblmod128(r, n);
			m >>= 1;
		}

		m = (d[1] >= 1) ? 1ULL << 63 : 1ULL << (62 - my_clz64(d[0]));   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);
			chkmod128(r, n);

			if (d[0] & m) dblmod128(r, n);
			m >>= 1;
		}
	}
	else
#endif
	{
		m = (d[1] <= 1) ? 0 : 1ULL << (62 - lzcnt);   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);

			if (d[1] & m) dblmod128(r, n);
			m >>= 1;
		}

		m = (d[1] >= 1) ? 1ULL << 63 : 1ULL << (62 - my_clz64(d[0]));   // set a mask at the leading bit - 2
		while (m > 0)
		{
			sqrmod128n(r, r, n, rho);

			if (d[0] & m) dblmod128(r, n);
			m >>= 1;
		}
	}

	chkmod128(r, n);
	if ((r[0] == one[0]) && (r[1] == one[1])) return 1;

	// now compute b^(2^r*d) for 0 <= r < s, 
	// and check for congruence to -1 as we go.
	uint64_t mone[2];
	submod128(n, one, mone, n);

	while (s > 1)
	{
		if ((r[0] == mone[0]) && (r[1] == mone[1])) return 1;
		sqrmod128n(r, r, n, rho);
		if (protect) chkmod128(r, n);
		if ((r[0] == one[0]) && (r[1] == one[1])) return 0;
		s--;
	}

	if ((r[0] == mone[0]) && (r[1] == mone[1])) return 1;
	else return 0;
}

/* *******************************************************************************
 * mpz_lucas_prp:
 * A "Lucas pseudoprime" with parameters (P,Q) is a composite n with D=P^2-4Q,
 * (n,2QD)=1 such that U_(n-(D/n)) == 0 mod n [(D/n) is the Jacobi symbol]
 * *******************************************************************************/
#if 0
int lucas_prp_128x1(uint64_t *n, long int p, long int q)
{
	uint64_t res[2];
	uint64_t index[2];
	uint128_t n128 = ((uint128_t)n[1] << 64) | (uint128_t)n[0];
	int s = 0, j = 0;
	int ret = 0;
	long int d = p * p - 4 * q;
	int sd = d < 0 ? 1 : 0;
	int sq = q < 0 ? 1 : 0;

	if (d == 0) /* Does not produce a proper Lucas sequence */
		return -1;

	if ((n[1] == 0) && (n[0] < 2))
		return 0;

	if ((n[0] & 1 == 0))
		return 0;

	if (sd)	d *= -1;
	if (sq) q *= -1;

	res[0] = (uint64_t)((uint64_t)d * (uint64_t)q * 2);
	res[1] = 0;
	bin_gcd128(res, n, res);

	if (((res[1] == n[1]) && (res[0] == n[0])) ||
		((res[1] == 0) && (res[0] == 1)))
	{

	}
	else
	{
		// if the gcd is anything other than 1 or n, return composite.
		return 0;
	}

	/* index = n-(D/n), where (D/n) is the Jacobi symbol */
	index[0] = n[0]; 
	index[1] = n[1];
	if (jacobi_128(d * sd, n128) < 0)
	{
		uint64_t c = (n[0] == 0xffffffffffffffffull) ? 1 : 0;
		index[0] += 1;
		index[1] += c;
	}
	else
	{
		uint64_t c = (n[0] == 0) ? 1 : 0;
		index[0] -= 1;
		index[1] -= c;
	}
		

	/* mpz_lucasumod(res, p, q, index, n); */
	uint64_t uh[2], vl[2], vh[2], ql[2], qh[2], tmp[2];
	uint64_t rho = multiplicative_inverse(n[0]);

	// signs
	int suh = 0;
	int svl = 0;
	int svh = 0;
	int sql = 0;
	int sqh = 0;
	int st = 0;

	// initialize our lucas variables into Montgomery representation
	uint64_t one[2];
	uint128_t unity128 = (uint128_t)0 - n128;
	unity128 = unity128 % n128;
	one[1] = (uint64_t)(unity128 >> 64);
	one[0] = (uint64_t)unity128;

	uh[0] = one[0];
	vl[0] = one[0];
	vh[0] = one[0]; // p is always 1. p;
	ql[0] = one[0];
	qh[0] = one[0];
	tmp[0] = 0;

	uh[1] = one[1];
	vl[1] = one[1];
	vh[1] = one[1];
	ql[1] = one[1];
	qh[1] = one[1];
	tmp[1] = 0;

	// vl = 2 = one + one
	tmp[0] = vl[0];
	vl[0] += one[0];
	vl[1] += one[1];
	vl[1] += (vl[0] < tmp[0]) ? 1 : 0;

	//s = mpz_scan1(index, 0);	// index of least significant 1-bit, i.e., ctz
	s = my_ctz128(index[0], index[1]);
	int sz = 128 - my_clz128(index[0], index[1]);

	for (j = sz - 1; j >= s + 1; j--)
	{
		/* ql = ql*qh (mod n) */
		mulmod128n(ql, qh, ql, n, rho);
		sql = sql ^ sqh;

		int bit = (j >= 64) ? index[1] & (1ull << (j - 64)) : index[0] & (1ull << j);

		if (bit > 0)
		{
			/* qh = ql*q */
			mpz_mul_si(qh, ql, q);

			/* uh = uh*vh (mod n) */
			mpz_mul(uh, uh, vh);
			mpz_mod(uh, uh, n);

			/* vl = vh*vl - p*ql (mod n) */
			mpz_mul(vl, vh, vl);
			mpz_mul_si(tmp, ql, p);
			mpz_sub(vl, vl, tmp);
			mpz_mod(vl, vl, n);

			/* vh = vh*vh - 2*qh (mod n) */
			mpz_mul(vh, vh, vh);
			mpz_mul_si(tmp, qh, 2);
			mpz_sub(vh, vh, tmp);
			mpz_mod(vh, vh, n);
		}
		else
		{
			/* qh = ql */
			mpz_set(qh, ql);

			/* uh = uh*vl - ql (mod n) */
			mpz_mul(uh, uh, vl);
			mpz_sub(uh, uh, ql);
			mpz_mod(uh, uh, n);

			/* vh = vh*vl - p*ql (mod n) */
			mpz_mul(vh, vh, vl);
			mpz_mul_si(tmp, ql, p);
			mpz_sub(vh, vh, tmp);
			mpz_mod(vh, vh, n);

			/* vl = vl*vl - 2*ql (mod n) */
			mpz_mul(vl, vl, vl);
			mpz_mul_si(tmp, ql, 2);
			mpz_sub(vl, vl, tmp);
			mpz_mod(vl, vl, n);
		}
	}
	/* ql = ql*qh */
	mpz_mul(ql, ql, qh);

	/* qh = ql*q */
	mpz_mul_si(qh, ql, q);

	/* uh = uh*vl - ql */
	mpz_mul(uh, uh, vl);
	mpz_sub(uh, uh, ql);

	/* vl = vh*vl - p*ql */
	mpz_mul(vl, vh, vl);
	mpz_mul_si(tmp, ql, p);
	mpz_sub(vl, vl, tmp);

	/* ql = ql*qh */
	mpz_mul(ql, ql, qh);

	for (j = 1; j <= s; j++)
	{
		/* uh = uh*vl (mod n) */
		mpz_mul(uh, uh, vl);
		mpz_mod(uh, uh, n);

		/* vl = vl*vl - 2*ql (mod n) */
		mpz_mul(vl, vl, vl);
		mpz_mul_si(tmp, ql, 2);
		mpz_sub(vl, vl, tmp);
		mpz_mod(vl, vl, n);

		/* ql = ql*ql (mod n) */
		mpz_mul(ql, ql, ql);
		mpz_mod(ql, ql, n);
	}

	mpz_mod(res, uh, n); /* uh contains our return value */


	if ((res[1] == 0) && (res[0] == 0))
	{
		return 1;
	}
	else
	{
		return 0;
	}

}/* method mpz_lucas_prp */
#endif

int pull_twos_128(int* n, int* j, uint128_t p)
{
	int c = 0;

	while (!(*n & 1))
	{
		*n >>= 1;
		c = 1 - c;
	}
	if ((c * (p * p - 1) % 16) == 8)
		*j *= -1;
	return c;
}

int jacobi_128(int n, uint128_t p)
{
	// compute the jacobi symbol (n/p)
	// p must be odd and positive
	// based on routine in Bressoud's book

	int j = 1;
	uint128_t t;
	int nn;
	int sign;

	// return an error condition if p is even
	if (!(p & 1))
		return -2;

	// pull out the (-1) power if n is negative
	// (-1 / p) = (-1)^((p-1)/2)
	if (n < 0)
	{
		if (((p - 1) / 2) & 1)
		{
			sign = -1;

		}
		else
		{
			sign = 1;
		}
		nn = -1 * n;
	}
	else
	{
		sign = 1;
		nn = n;
	}

	nn = nn % p;

	// if p divides n then (n/p) = 0
	if (nn == 0)
		return 0;

	pull_twos_128(&nn, &j, p);
	while (nn > 1)
	{
		if (((nn - 1) * (p - 1)) % 8 == 4)
			j = -1 * j;
		t = (uint128_t)nn;
		nn = p % (uint128_t)nn;
		p = t;

		pull_twos_128(&nn, &j, p);
	}
	return sign * j;
}

// A lucas PRP with Selfridge parameters on 1 128-bit input (two 64-bit words) (in-progress)
int selfridge_prp_128x1(uint64_t* n)
{
	int result;
	mpz_t gn;
	mpz_init(gn);
	mpz_set_ui(gn, n[1]);
	mpz_mul_2exp(gn, gn, 64);
	mpz_add_ui(gn, gn, n[0]);

#if 1
	long int d = 5, p = 1, q = 0;
	int max_d = 1000000;
	int jacobi = 0;
	uint128_t n128;

	n128 = n[1];
	n128 <<= 64;
	n128 |= (uint128_t)n[0];

	if ((n[1] == 0) && (n[0] < 2))
		return 0;

	if ((n[0] & 1 == 0))
		return 0;

	while (1)
	{
		jacobi = jacobi_128(d, n128);

		/* if jacobi == 0, d is a factor of n, therefore n is composite... */
		/* if d == n, then either n is either prime or 9... */
		if (jacobi == 0)
		{
			if ((d == n128) && (d != 9))
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
		if (jacobi == -1)
			break;

		/* if we get to the 5th d, make sure we aren't dealing with a square... */
		if (d == 13)
		{
			if (mpz_perfect_square_p(gn))
			{
				return 0;
			}
		}

		if (d < 0)
		{
			d *= -1;
			d += 2;
		}
		else
		{
			d += 2;
			d *= -1;
		}

		/* make sure we don't search forever */
		if (d >= max_d)
		{
			return -2;
		}
	}

	q = (1 - d) / 4;

	result = mpz_lucas_prp(gn, p, q); // lucas_prp_128x1(gn, p, q);


#else
	int result = mpz_extrastrongselfridge_prp(gn);
#endif

	mpz_clear(gn);
	return result;
}

// a BPSW test on 1 128-bit input
int bpsw_prp_128x1(uint64_t* n)
{
	if (MR_2sprp_128x1(n) == 0)
		return 0;

	return selfridge_prp_128x1(n);
}

static uint128_t seed = ((uint128_t)0x123456789ull << 92) + ((uint128_t)0xabcdef << 36) + 0x987654321ull;
static uint128_t my_random(void)
{
	// based on linear congruential generator, period = 2^128
	seed = seed * 137 + 13;
	// shuffle
	uint128_t x = seed ^ (seed >> 17) ^ (seed << 13);
	return x;
}
static my_random_reset(void)
{
	seed = ((uint128_t)0x123456789ull << 92) + ((uint128_t)0xabcdef << 36) + 0x987654321ull;
}

static inline uint64_t my_rdtsc(void)
{
#if defined(__x86_64__)
	// supported by GCC and Clang for x86 platform
	return __rdtsc();
#elif INLINE_ASM && defined(__aarch64__)
	// should be a 64 bits wallclock counter
	// document for old/recent architecture and/or BMC chipsets mention it
	// could be a 56 bit counter.
	uint64_t val;

	asm volatile ("mrs %0, cntvct_el0":"=r" (val));

	// I am not sure what the clock unit is, it depends on pre-scaler setup
	// A multiplication by 32 might be needed on my platform 
	return val * 32;	// aarch64 emulation on x86_64 ?
	return ((val / 3) * 25) << 4;	// maybe for ARM M1 ?
	return val;
#else
#error "todo : unsupported _rdtsc implementation\n"
	return 0;
#endif
}

static uint64_t sm_primes[168] = {
	2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
	103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
	199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
	313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
	433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
	563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
	673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
	811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
	941,947,953,967,971,977,983,991,997
};

int test_tinyprp()
{
	uint64_t prp[16];
	int correct = 0;
	int i, k;

#ifdef GMP_CHECK
	mpz_t gmp2, gmpn, gmpn1;
	mpz_init(gmp2);
	mpz_init(gmpn);
	mpz_init(gmpn1);
#endif

	struct timeval start, stop;
	int bits;
	uint64_t elapsed;

	mpz_t gmpb;
	mpz_t gmpe;
	mpz_t gmp1;
	mpz_t gmpm1;
	mpz_t gmpm;
	mpz_t gmpa;
	mpz_t gmpmm1;

	mpz_init(gmpb);
	mpz_init(gmpe);
	mpz_init(gmp1);
	mpz_init(gmpm1);
	mpz_init(gmpmm1);
	mpz_init(gmpm);
	mpz_init(gmpa);
	
	printf("test of mpz_powm +/-1 on random 6k + 1 inputs\n");
	my_random_reset();
	for (bits = 80; bits <= 128; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t totalnum = 0;
		uint32_t num = 100;
		double telapsed = 0;
		int k;

		uint32_t num1 = 0, numm1 = 0;
		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;
			prp[1] = (uint64_t)(x >> 64); // &0xfffffffffffffull;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int y;

			for (y = 0; y < num; y++)
			{
				mpz_set_ui(gmpm, prp[1]);
				mpz_mul_2exp(gmpm, gmpm, 64);
				mpz_add_ui(gmpm, gmpm, prp[0]);
				mpz_sub_ui(gmpmm1, gmpm, 1);
				mpz_tdiv_q_2exp(gmpe, gmpmm1, 1);

				int j;
				for (j = 0; j < 168; j++)
				{
					mpz_set_ui(gmpb, sm_primes[j]);
					mpz_powm(gmpa, gmpb, gmpe, gmpm);
					if (mpz_cmp_ui(gmpa, 1) == 0)
						num1++;
					if (mpz_cmp(gmpa, gmpmm1) == 0)
						numm1++;
				}

				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

			totalnum += num;
		} while (totalnum < 10000); // (0);  (elapsed < (1ull << 30));

		printf("%013lx%013lx\n", prp[1], prp[0]);

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %u +1 and %u -1 tests out of %u %d-bit inputs: %1.2f%%\n",
			num1, numm1, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	printf("test of modexp_104x8b +/-1 on random 6k + 1 inputs\n");
	my_random_reset();
	for (bits = 80; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t totalnum = 0;
		uint32_t num = 100;
		double telapsed = 0;
		int k;

		uint32_t num1 = 0, numm1 = 0;
		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x & 0xfffffffffffffull;
			prp[1] = (uint64_t)(x >> 52) & 0xfffffffffffffull;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int y;

			for (y = 0; y < num; y++)
			{
				int j;

				uint64_t e[2];
				e[0] = prp[0];
				e[1] = prp[1];

				e[0] -= 1;		// is odd: won't carry.
				e[0] >>= 1;
				e[0] |= (e[1] << 51);
				e[1] >>= 1;
				e[0] &= 0xfffffffffffffull;

#ifdef GMP_CHECK
				mpz_set_ui(gmpn, prp[1]);
				mpz_mul_2exp(gmpn, gmpn, 52);
				mpz_add_ui(gmpn, gmpn, prp[0]);
				mpz_sub_ui(gmpn1, gmpn, 1);
				mpz_tdiv_q_2exp(gmpe, gmpn1, 1);
#endif

				uint8_t onemsk = 0;
				uint8_t monemsk = 0;

				for (j = 0; j < 21; j++)
				{
					modexp_104x8b(&onemsk, &monemsk, prp, e, sm_primes + j * 8);
					num1 += _mm_popcnt_u32(onemsk);
					numm1 += _mm_popcnt_u32(monemsk);

#ifdef GMP_CHECK
					uint8_t gonemsk = 0;
					uint8_t gmonemsk = 0;

					int m;
					for (m = 0; m < 8; m++)
					{
						mpz_set_ui(gmpb, sm_primes[j * 8 + m]);
						mpz_powm(gmp2, gmpb, gmpe, gmpn);

						if (mpz_cmp_ui(gmp2, 1) == 0)
							gonemsk |= (1 << m);
						if (mpz_cmp(gmp2, gmpn1) == 0)
							gmonemsk |= (1 << m);
					}


					if (onemsk != gonemsk)
					{
						printf("1 masks not equal: %02x,%02x\n", onemsk, gonemsk);
						gmp_printf("input %Zd\n", gmpn);
						printf("primes: ");
						for (m = 0; m < 8; m++)
						{
							if (((onemsk ^ gonemsk) & (1 << m)) > 0)
							{
								mpz_set_ui(gmpb, sm_primes[j * 8 + m]);
								mpz_powm(gmp2, gmpb, gmpe, gmpn);
								gmp_printf("%u: %Zd ", sm_primes[j * 8 + m], gmp2);
							}
						}
						printf("\n");
					}

					if (monemsk != gmonemsk)
					{
						printf("-1 masks not equal: %02x,%02x\n", monemsk, gmonemsk);
						gmp_printf("input %Zd\n", gmpn);
						printf("primes: ");
						for (m = 0; m < 8; m++)
						{
							if (((monemsk ^ gmonemsk) & (1 << m)) > 0)
							{
								mpz_set_ui(gmpb, sm_primes[j * 8 + m]);
								mpz_powm(gmp2, gmpb, gmpe, gmpn);
								gmp_printf("%u: %Zd ", sm_primes[j * 8 + m], gmp2);
							}
						}
						printf("\n");
					}

					if ((onemsk != gonemsk) || (monemsk != gmonemsk))
						exit(0);
#endif
				}

				prp[0] += inc;
				inc = 6 - inc;
			}

			

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

			totalnum += num;
		} while (totalnum < 10000); // (0);  (elapsed < (1ull << 30));

		printf("%013lx%013lx\n", prp[1], prp[0]);

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %u +1 and %u -1 tests out of %u %d-bit inputs: %1.2f%%\n",
			num1, numm1, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	printf("test of modexp_128x1b on random 6k + 1 inputs\n");
	for (bits = 80; bits <= 128; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t totalnum = 0;
		uint32_t num = 100;
		double telapsed = 0;
		int k;

		uint32_t num1 = 0, numm1 = 0;
		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;
			prp[1] = (uint64_t)(x >> 64);

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int y;

			for (y = 0; y < num; y++)
			{
				int j;

				uint64_t e[2];
				e[0] = prp[0];
				e[1] = prp[1];

				e[0] -= 1;		// is odd: won't carry.
				e[0] >>= 1;
				e[0] |= (e[1] << 63);
				e[1] >>= 1;

				for (j = 0; j < 168; j++)
				{
					int result = modexp_128x1b(prp, e, sm_primes[j]);
					if (result == 1)
						num1++;
					if (result == -1)
						numm1++;
				}

				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

			totalnum += num;
		} while (totalnum < 10000); // (0);  (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %u +1 and %u -1 tests out of %u %d-bit inputs: %1.2f%%\n",
			num1, numm1, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	return 0;

	printf("test of fermat_prp_64x1 on random 6k + 1 inputs\n");
	my_random_reset();
	for (bits = 20; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += fermat_prp_64x1(prp[0]);
				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	printf("test of MR_2sprp_64x1 on random 6k + 1 inputs\n");
	my_random_reset();
	for (bits = 20; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += MR_2sprp_64x1(prp[0]);
				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %d MR_2sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	printf("test of mpz_probab_prime_p on random 6k + 1 inputs\n");
	my_random_reset();
	for (bits = 20; bits <= 0; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;
		mpz_t gmp_prp;

		mpz_init(gmp_prp);

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			mpz_set_ui(gmp_prp, (uint64_t)x);

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += (mpz_probab_prime_p(gmp_prp, 1) > 0);
				mpz_add_ui(gmp_prp, gmp_prp, inc);
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		mpz_clear(gmp_prp);

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %d mpz_probab_prime_p out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	printf("test of fermat_prp_128x1 on random 6k + 1 inputs\n");
	for (bits = 100; bits <= 128; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;
			prp[1] = (uint64_t)(x >> 64);

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += fermat_prp_128x1(prp);
				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	printf("test of MR_2sprp_128x1 on random 6k + 1 inputs\n");
	for (bits = 100; bits <= 128; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {

			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);

			prp[0] = (uint64_t)x;
			prp[1] = (uint64_t)(x >> 64);

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += MR_2sprp_128x1(prp);
				prp[0] += inc;
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num));
		printf("found %d MR-2sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num, bits, 100. * (double)numprp / (double)(k * num));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num));
	}
	printf("\n");

	return 0;
	// test of fermat_prp_52x8 on random 6k+1 inputs
	for (bits = 20; bits <= 52; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;
		if (bits > 52) bits = 52;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint64_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
			}

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += _mm_popcnt_u32(fermat_prp_52x8(prp));
				for (i = 0; i < 8; i++)
				{
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// test of fermat_prp_104x8 on random 6k+1 inputs
	for (bits = 80; bits <= 104; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;
		int k;
		if (bits > 104) bits = 104;

		k = 0;
		int fail = 0;
		elapsed = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint128_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
				prp[i + 8] = (uint64_t)(x >> 52) & 0xfffffffffffffull;
			}

			uint64_t inc = 4;
			int j;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			for (j = 0; j < num; j++)
			{
				uint8_t prpmask = fermat_prp_104x8(prp);
				numprp += _mm_popcnt_u32(prpmask);
				for (i = 0; i < 8; i++)
				{
#ifdef GMP_CHECK
					if ((prpmask & (1 << i)) == 0)
					{
						mpz_set_ui(gmp2, 2);
						mpz_set_ui(gmpn, prp[8 + i]);
						mpz_mul_2exp(gmpn, gmpn, 52);
						mpz_add_ui(gmpn, gmpn, prp[i]);
						mpz_sub_ui(gmpn1, gmpn, 1);
						mpz_powm(gmp2, gmp2, gmpn1, gmpn);
						if (mpz_cmp_ui(gmp2, 1) == 0)
						{
							//printf("prp %016lx%016lx failed in lane %d\n", prp[i+8], prp[i], i);
							//gmp_printf("mpz result = %Zx\n", gmp2);
							//exit(1);
							fail++;
						}
					}
#endif
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);
		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d fermat-prp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
#ifdef GMP_CHECK
		printf("GMP checked number of failures: %d\n", fail);
#endif
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// test of fermat_prp_52x8 on random 6k+1 inputs
	for (bits = 20; bits <= 52; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 1000000;
		double telapsed = 0;
		int k;
		if (bits > 52) bits = 52;

		numprp = 0;
		k = 0;
		elapsed = 0;
		telapsed = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint64_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
			}

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			uint64_t inc = 4;
			int j;

			for (j = 0; j < num; j++)
			{
				numprp += _mm_popcnt_u32(MR_2sprp_52x8(prp));
				for (i = 0; i < 8; i++)
				{
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);

		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d MR-2sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// test of MR_2sprp_104x8 on random 6k+1 inputs
	for (bits = 50; bits <= 104; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;
		int k;
		if (bits > 104) bits = 104;

		k = 0;
		elapsed = 0;
		int fail = 0;
		do {
			for (i = 0; i < 8; i++)
			{
				uint128_t x;
				do {
					x = my_random();
					uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
					uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
					x &= maskAnd;
					x |= maskOr;
					x /= 6;
					x *= 6;	// now a multiple of 6
					x += 1;	// number like 6*k + 1
				} while (x >> (bits - 1) != 1);
				prp[i] = (uint64_t)x & 0xfffffffffffffull;
				prp[i + 8] = (uint64_t)(x >> 52) & 0xfffffffffffffull;
			}

			uint64_t inc = 4;
			int j;

			ticks1 = my_rdtsc();
			gettimeofday(&start, NULL);

			for (j = 0; j < num; j++)
			{
				uint8_t prpmask = MR_2sprp_104x8(prp);
				numprp += _mm_popcnt_u32(prpmask);
				for (i = 0; i < 8; i++)
				{
#ifdef GMP_CHECK
					if ((prpmask & (1 << i)) == 0)
					{
						mpz_set_ui(gmp2, 2);
						mpz_set_ui(gmpn, prp[8 + i]);
						mpz_mul_2exp(gmpn, gmpn, 52);
						mpz_add_ui(gmpn, gmpn, prp[i]);
						mpz_sub_ui(gmpn1, gmpn, 1);
						mpz_powm(gmp2, gmp2, gmpn1, gmpn);
						if (mpz_cmp_ui(gmp2, 1) == 0)
						{
							//printf("prp %016lx%016lx failed in lane %d\n", prp[i + 8], prp[i], i);
							//gmp_printf("mpz result = %Zx\n", gmp2);
							//exit(1);
							fail++;
						}
					}
#endif
					prp[i] += inc;
				}
				inc = 6 - inc;
			}

			k++;
			ticks2 = my_rdtsc();
			elapsed += (ticks2 - ticks1);
			gettimeofday(&stop, NULL);
			telapsed += ytools_difftime(&start, &stop);
		} while (elapsed < (1ull << 30));

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / (k * num * 8));
		printf("found %d MR-2sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, k * num * 8, bits, 100. * (double)numprp / (double)(k * num * 8));
#ifdef GMP_CHECK
		printf("GMP checked number of failures: %d\n", fail);
#endif
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)(k * num * 8));
	}
	printf("\n");

	// bases for MR-sprp check:
	uint64_t bases[24] = { 3, 5, 7, 11,
		13, 17, 19, 23,
		29, 31, 37, 41,
		43, 47, 53, 59,
		61, 67, 71, 73,
		79, 83, 89, 97 };

	// test of MR_sprp_104x8 on PRP 6k+1 inputs
	for (bits = 50; bits <= 104; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1 = my_rdtsc();
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;

		if (bits > 104) bits = 104;

		//printf("commencing test of random 6k+1 %d-bit inputs\n", bits);
		elapsed = 0;

		uint64_t inc[8] = { 4, 4, 4, 4, 4, 4, 4, 4 };

		for (i = 0; i < 8; i++)
		{
			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);
			prp[i] = (uint64_t)x & 0xfffffffffffffull;
			prp[i + 8] = (uint64_t)(x >> 52);
		}

		uint8_t isprp = 0;
		while (isprp != 0xff)
		{
			isprp = fermat_prp_104x8(prp);

			for (i = 0; i < 8; i++)
			{
				if ((isprp & (1 << i)) == 0)
				{
					prp[i] += inc[i];
					inc[i] = 6 - inc[i];
				}
			}
		}

		ticks1 = my_rdtsc();

		uint64_t basecount[8];
		uint64_t maxcount[8];
		uint64_t currentbase[8];
		for (i = 0; i < 8; i++)
		{
			//printf("prp%d = %016lx%016lx\n", i, prp[i+8], prp[i]); fflush(stdout);
			basecount[i] = 0;
			currentbase[i] = bases[0];
			if (bits <= 62)
				maxcount[i] = 8;
			else if (bits <= 82)
				maxcount[i] = 16;	// was 12, but need a multiple of 8 (same cost anyway with avx512)
			else if (bits <= 112)
				maxcount[i] = 16;
			else
				maxcount[i] = 24;	// was 20, but need a multiple of 8 (same cost anyway with avx512)
		}

		uint32_t tested = 0;
		gettimeofday(&start, NULL);

		elapsed = 0;
		int fail = 0;
		while ((elapsed < (1ull << 30)))
		{
			uint8_t prpmask = MR_sprp_104x8(prp, currentbase);
			for (i = 0; i < 8; i++)
			{
#ifdef GMP_CHECK
				if ((prpmask & (1 << i)) == 0)
				{
					mpz_set_ui(gmpn, prp[8 + i]);
					mpz_mul_2exp(gmpn, gmpn, 52);
					mpz_add_ui(gmpn, gmpn, prp[i]);
					if (mpz_probab_prime_p(gmpn, 1) > 0)
					{
						//printf("prp %016lx%016lx failed in lane %d\n", prp[i + 8], prp[i], i);
						//gmp_printf("mpz result = %Zx\n", gmp2);
						//exit(1);
						fail++;
					}
				}
#endif

				if (prpmask & (1 << i))
				{
					// the input in position i could be prp, increment the
					// base if there are more of them
					if (basecount[i] < maxcount[i])
					{
						basecount[i]++;
						currentbase[i] = bases[basecount[i]];
					}
					else
					{
						// we've tested enough bases to know this is prime
						tested++;
						numprp++;
						//prp[i] += inc[i];
						//inc[i] = 6 - inc[i];
						currentbase[i] = bases[0];
						basecount[i] = 0;
					}
				}
				else
				{
					// the input in position i is definitely not prime,
					// replace it and increment num tested
					tested++;
					//prp[i] += inc[i];
					//inc[i] = 6 - inc[i];
					currentbase[i] = bases[0];
					basecount[i] = 0;
				}
			}

			ticks2 = my_rdtsc();
			elapsed = (ticks2 - ticks1);
		}

		gettimeofday(&stop, NULL);
		telapsed += ytools_difftime(&start, &stop);

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / tested);
		printf("found %d MR-sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, tested, bits, 100. * (double)numprp / (double)tested);
#ifdef GMP_CHECK
		printf("GMP checked number of failures: %d\n", fail);
#endif
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)tested);
	}
	printf("\n");

	// test of MR_sprp_104x8base on PRP 6k+1 inputs
	for (bits = 50; bits <= 104; bits += 1)
	{
		uint32_t numprp = 0;
		uint64_t ticks1;
		uint64_t ticks2;
		uint32_t num = 100000;
		double telapsed = 0;
		uint64_t one[16];

		//printf("commencing test of random 6k+1 %d-bit inputs\n", bits);
		elapsed = 0;

		if (bits > 104) bits = 104;

		uint64_t inc[8] = { 4, 4, 4, 4, 4, 4, 4, 4 };

		for (i = 0; i < 8; i++)
		{
			uint128_t x;
			do {
				x = my_random();
				uint128_t maskAnd = ((uint128_t)1 << (bits - 1)) - 1;	// clear msbits
				uint128_t maskOr = ((uint128_t)1 << (bits - 1)) | ((uint128_t)1 << (bits / 2));	// force msb, force another bit
				x &= maskAnd;
				x |= maskOr;
				x /= 6;
				x *= 6;	// now a multiple of 6
				x += 1;	// number like 6*k + 1
			} while (x >> (bits - 1) != 1);
			prp[i] = (uint64_t)x & 0xfffffffffffffull;
			prp[i + 8] = (uint64_t)(x >> 52);
		}

		uint8_t isprp = 0;
		while (isprp != 0xff)
		{
			isprp = fermat_prp_104x8(prp);

			for (i = 0; i < 8; i++)
			{
				//printf("%016lx%016lx : %u (%u)\n", prp[i+8], prp[i], isprp & (1 << i), isprp);
				if ((isprp & (1 << i)) == 0)
				{
					prp[i] += inc[i];
					inc[i] = 6 - inc[i];
				}
			}
		}

		uint128_t o128;
		uint128_t n128;
		for (i = 0; i < 8; i++)
		{
			//printf("prp%d = %016lx%016lx\n", i, prp[i+8], prp[i]); fflush(stdout);
			n128 = ((uint128_t)prp[i + 8] << 52) + prp[i];
			o128 = (uint128_t)1 << 104;
			o128 = o128 % n128;
			one[i] = (uint64_t)o128 & 0xfffffffffffffull;
			one[i + 8] = (uint64_t)(o128 >> 52);
			//printf("one%d = %016lx%016lx\n", i, one[i + 8], one[i]); fflush(stdout);
		}

		uint64_t basecount = 0;
		uint64_t maxcount;
		uint64_t currentbase[8];

		if (bits <= 62)
			maxcount = 8;
		else if (bits <= 82)
			maxcount = 16;	// was 12, but need a multiple of 8 (same cost anyway with avx512)
		else if (bits <= 112)
			maxcount = 16;
		else
			maxcount = 24;	// was 20, but need a multiple of 8 (same cost anyway with avx512)

		for (i = 0; i < 8; i++)
		{
			currentbase[i] = bases[i];
		}

		uint32_t tested = 0;
		gettimeofday(&start, NULL);

		ticks1 = my_rdtsc();

		elapsed = 0;
		int tnum = 0;
		while ((elapsed < (1ull << 30)))
		{
			uint64_t ntest[2], otest[2];
			ntest[1] = prp[tnum + 8];
			ntest[0] = prp[tnum];
			otest[1] = one[tnum + 8];
			otest[0] = one[tnum];

			// so far, simple RL is faster than kary-LR
			uint8_t prpmask = MR_sprp_104x8base(ntest, otest, currentbase);

			if (prpmask == 0xff)
			{
				// the input is prp to all current bases, increment the
				// base if there are more of them
				if ((basecount + 8) <= maxcount)
				{
					for (i = 0; i < 8; i++)
					{
						currentbase[i] = bases[basecount + i];
					}
					basecount += 8;
				}
				else
				{
					// we've tested enough bases to know this is prime
					tested++;
					numprp++;
					basecount = 0;
					for (i = 0; i < 8; i++)
					{
						currentbase[i] = bases[basecount + i];
					}
					tnum = (tnum + 1) & 7;
				}
			}
			else
			{
				// the input in position i is definitely not prime,
				// replace it and increment num tested
				tested++;
				basecount = 0;
				for (i = 0; i < 8; i++)
				{
					currentbase[i] = bases[basecount + i];
				}
				tnum = (tnum + 1) & 7;
			}

			ticks2 = my_rdtsc();
			elapsed = (ticks2 - ticks1);
		}

		gettimeofday(&stop, NULL);
		telapsed += ytools_difftime(&start, &stop);

		printf("total ticks = %lu, ticks per %d-bit input = %lu\n",
			elapsed, bits, (elapsed) / tested);
		printf("found %d MR-sprp out of %u %d-bit inputs: %1.2f%%\n",
			numprp, tested, bits, 100. * (double)numprp / (double)tested);
		printf("elapsed time: %1.4f sec, %1.4f us / input\n", telapsed, 1000000. * telapsed / (double)tested);
	}
	printf("\n");

#ifdef GMP_CHECK
	mpz_clear(gmpn);
	mpz_clear(gmpn1);
	mpz_clear(gmp2);
#endif
	return 0;
}

#ifdef USE_AVX512F


__m512i rem_epu64_x8(__m512i n, __m512i d)
{
	// DANGER: I haven't proven this works for every possible input.
	__m512d d1pd = _mm512_cvtepu64_pd(d);
	__m512d n1pd = _mm512_cvtepu64_pd(n);
	__m512i q, q2, r;

	//n1pd = _mm512_div_pd(n1pd, d1pd);
	//q = _mm512_cvt_roundpd_epu64(n1pd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	n1pd = _mm512_div_round_pd(n1pd, d1pd, ROUNDING_MODE);
	q = _mm512_cvttpd_epu64(n1pd);

	__m512i qd = _mm512_mullox_epi64(q, d);
	r = _mm512_sub_epi64(n, qd);

	// fix q too big by a little, with special check for
	// numerators close to 2^64 and denominators close to 1
	// DANGER: the special check is unused for 64-bits, only for 32-bits.
	// This routine is only used in modmul32 and input numerators
	// shouldn't get that large in normal cases.  The factor base
	// would need to be close to 2^32...
	__mmask8 err = _mm512_cmpgt_epu64_mask(r, n); // |
		//(_mm512_cmpgt_epu64_mask(r, d) & _mm512_cmplt_epu64_mask(
		//	_mm512_sub_epi64(_mm512_set1_epi64(0), r), _mm512_set1_epi64(1024)));
	if (err)
	{
		n1pd = _mm512_cvtepu64_pd(_mm512_sub_epi64(_mm512_set1_epi64(0), r));

		//n1pd = _mm512_div_pd(n1pd, d1pd);
		//q2 = _mm512_add_epi64(_mm512_set1_epi64(1), _mm512_cvt_roundpd_epu64(n1pd,
		//	(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

		n1pd = _mm512_div_round_pd(n1pd, d1pd, ROUNDING_MODE);
		q2 = _mm512_add_epi64(_mm512_set1_epi64(1), _mm512_cvttpd_epu64(n1pd));

		q = _mm512_mask_sub_epi64(q, err, q, q2);
		r = _mm512_mask_add_epi64(r, err, r, _mm512_mullox_epi64(q2, d));
	}

	// fix q too small by a little bit
	err = _mm512_cmpge_epu64_mask(r, d);
	if (err)
	{
		n1pd = _mm512_cvtepu64_pd(r);

		//n1pd = _mm512_div_pd(n1pd, d1pd);
		//q2 = _mm512_cvt_roundpd_epu64(n1pd,
		//	(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

		n1pd = _mm512_div_round_pd(n1pd, d1pd, ROUNDING_MODE);
		q2 = _mm512_cvttpd_epu64(n1pd);

		q = _mm512_mask_add_epi64(q, err, q, q2);
		r = _mm512_mask_sub_epi64(r, err, r, _mm512_mullox_epi64(q2, d));
	}

	return r;
}

// a Fermat PRP test on 8x 52-bit inputs
uint8_t fermat_prp_52x8(uint64_t* n)
{
	// assumes has no small factors.  assumes n <= 52 bits.
	// assumes n is a list of 8 52-bit integers
	// do a base-2 fermat prp test on each using LR binexp.
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__m512i vrho = multiplicative_inverse104_x8(n);
	__m512i unity;
	__m512i r;
	__m512i nvec;
	__m512i evec;
	__m512i m;
	__m512i zero = _mm512_setzero_si512();
	__m512i one = _mm512_set1_epi64(1);

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zero, vrho), lo52mask);
	nvec = loadu64(n);
	evec = _mm512_sub_epi64(nvec, one);

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
	r = _mm512_rem_epu64(_mm512_set1_epi64(1ULL << 52), nvec);
#else
	r = rem_epu64_x8(_mm512_set1_epi64(1ULL << 52), nvec);
#endif

	// penultimate-hi-bit mask
	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	unity = r;

	r = _mm512_add_epi64(r, r);
	__mmask8 ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	while (_mm512_cmpgt_epu64_mask(m, zero))
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec);
		mulredc52_mask_add_vec(&r, bitcmp, r, r, nvec, vrho);
		m = _mm512_srli_epi64(m, 1);
	}

	// AMM possibly needs a final correction by n
	ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	return _mm512_cmpeq_epu64_mask(unity, r);
}

// a Fermat PRP test on 8x 104-bit inputs
uint8_t fermat_prp_104x8(uint64_t* n)
{
	// assumes has no small factors.  assumes n >= 54 bits.
	// assumes n is a list of 8 104-bit integers (16 52-bit words)
	// in the format: 8 lo-words, 8 hi-words.
	// do a base-2 fermat prp test on each using LR binexp.
	__m512i vrho = multiplicative_inverse104_x8(n);
	vec_u104_t unity;
	__m512i nvec[2];
	__m512i evec[2];
	__m512i m;
	__m512i zero = _mm512_setzero_si512();
	__m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	uint64_t tmp = 0;

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zero, vrho), lo52mask);

	nvec[0] = loadu64(&n[0]);
	nvec[1] = loadu64(&n[8]);
	submod104_x8(&evec[1], &evec[0],
		nvec[1], nvec[0], zero, one, nvec[1], nvec[0]);

	// the 128-bit division we do the slow way
	int i;
	for (i = 0; i < 8; i++)
	{
		uint128_t mod = ((uint128_t)n[i + 8] << 52) + n[i];
		uint128_t t = (uint128_t)1 << 104;
		t %= mod;

		unity.data[0][i] = (uint64_t)t & 0xfffffffffffffULL;
		unity.data[1][i] = (uint64_t)(t >> 52);
	}

	// penultimate-hi-bit mask
	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec[1]));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
	m = _mm512_mask_set1_epi64(m, _mm512_cmple_epu64_mask(evec[1], one), 0);

	__mmask8 protect = _mm512_cmpgt_epi64_mask(_mm512_srli_epi64(evec[1], 49), zero);

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	// Note: the first 5 iterations can be done much more cheaply in
	// single precision and then converted into montgomery representation,
	// but that would require a 208-bit division; not worth it.
	__m512i r0 = loadu64(unity.data[0]);
	__m512i r1 = loadu64(unity.data[1]);

	addmod104_x8(&r1, &r0, r1, r0, r1, r0, nvec[1], nvec[0]);

	__mmask8 done;
	if (protect)
	{
		done = _mm512_cmpeq_epu64_mask(m, zero);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[1]);

			mask_sqrredc104_exact_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
			mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zero);
		}
	}
	else
	{

		__mmask8 done = _mm512_cmpeq_epu64_mask(m, zero);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[1]);

			mask_sqrredc104_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
			mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zero);
		}
	}

	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec[0]));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
	m = _mm512_mask_set1_epi64(m, _mm512_cmpge_epu64_mask(evec[1], one), 1ULL << 51);

	if (protect)
	{
		done = _mm512_cmpeq_epu64_mask(m, zero);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[0]);

			mask_sqrredc104_exact_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
			mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zero);
		}
	}
	else
	{

		__mmask8 done = _mm512_cmpeq_epu64_mask(m, zero);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec[0]);

			mask_sqrredc104_vec(&r1, &r0, ~done, r1, r0, nvec[1], nvec[0], vrho);
			mask_dblmod104_x8(&r1, &r0, (~done) & bitcmp, r1, r0, nvec[1], nvec[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zero);
		}
	}

	// AMM possibly needs a final correction by n
	addmod104_x8(&r1, &r0, r1, r0, zero, zero, nvec[1], nvec[0]);

	uint8_t isprp =
		_mm512_cmpeq_epu64_mask(loadu64(unity.data[0]), r0) &
		_mm512_cmpeq_epu64_mask(loadu64(unity.data[1]), r1);

	return isprp;
}

// a Miller-Rabin SPRP test on 8x 52-bit inputs using base 2
uint8_t MR_2sprp_52x8(uint64_t* n)
{
	// assumes has no small factors.  assumes n <= 52 bits.
	// assumes n is a list of 8 52-bit integers
	// TODO: do a base-2 MR sprp test on each using LR binexp.
	// as of now this is Fermat test again...
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	__m512i vrho = multiplicative_inverse104_x8(n);
	__m512i unity;
	__m512i r;
	__m512i nvec;
	__m512i evec;
	__m512i m;
	__m512i zero = _mm512_setzero_si512();
	__m512i one = _mm512_set1_epi64(1);

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zero, vrho), lo52mask);
	nvec = loadu64(n);
	evec = _mm512_sub_epi64(nvec, one);

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
	r = _mm512_rem_epu64(_mm512_set1_epi64(1ULL << 52), nvec);
#else
	r = rem_epu64_x8(_mm512_set1_epi64(1ULL << 52), nvec);
#endif

	// penultimate-hi-bit mask
	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(evec));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	unity = r;

	r = _mm512_add_epi64(r, r);
	__mmask8 ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	while (_mm512_cmpgt_epu64_mask(m, zero))
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, evec);
		mulredc52_mask_add_vec(&r, bitcmp, r, r, nvec, vrho);
		m = _mm512_srli_epi64(m, 1);
	}

	// AMM possibly needs a final correction by n
	ge = _mm512_cmpge_epi64_mask(r, nvec);
	r = _mm512_mask_sub_epi64(r, ge, r, nvec);

	return _mm512_cmpeq_epu64_mask(unity, r);
}

// a Miller-Rabin SPRP test on 8x 104-bit inputs using base 2
uint8_t MR_2sprp_104x8(uint64_t* n)
{
	// assumes has no small factors.  assumes n >= 54 bits.
	// assumes n is a list of 8 104-bit integers (16 52-bit words)
	// in the format: 8 lo-words, 8 hi-words.
	// do a Miller-Rabin sprp test using base 2.
	__m512i vrho = multiplicative_inverse104_x8(n);
	__m512i mone[2];
	vec_u104_t unity;
	__m512i nv[2];
	__m512i dv[2];
	__m512i rv[2];
	__m512i bv[2];
	__m512i n1v[2];
	__m512i tv[2];
	__m512i m;
	__m512i zerov = _mm512_setzero_si512();
	__m512i onev = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	uint64_t tmp = 0;

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zerov, vrho), lo52mask);

	nv[0] = loadu64(&n[0]);
	nv[1] = loadu64(&n[8]);

	// the 128-bit division we do one at a time
	int i;
	for (i = 0; i < 8; i++)
	{
		uint128_t mod = ((uint128_t)n[i + 8] << 52) + n[i];
		uint128_t one = (uint128_t)1 << 104;
		one %= mod;

		unity.data[0][i] = (uint64_t)one & 0xfffffffffffffULL;
		unity.data[1][i] = (uint64_t)(one >> 52) & 0xfffffffffffffULL;
	}

	mone[0] = loadu64(unity.data[0]);
	mone[1] = loadu64(unity.data[1]);

	// compute d and tzcnt
	submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], zerov, onev, nv[1], nv[0]);

	__mmask8 done = 0;
	dv[1] = n1v[1];
	dv[0] = n1v[0];
	__m512i tzcntv = zerov;
	while (done != 0xff)
	{
		__m512i c = _mm512_mask_slli_epi64(dv[1], ~done, dv[1], 51);
		dv[0] = _mm512_mask_srli_epi64(dv[0], ~done, dv[0], 1);
		dv[0] = _mm512_mask_or_epi64(dv[0], ~done, c, dv[0]);
		dv[1] = _mm512_mask_srli_epi64(dv[1], ~done, dv[1], 1);
		tzcntv = _mm512_mask_add_epi64(tzcntv, ~done, tzcntv, onev);
		done = done | _mm512_cmpeq_epi64_mask(_mm512_and_epi64(dv[0], onev), onev);
	}
	dv[0] = _mm512_and_epi64(dv[0], lo52mask);

	// penultimate-hi-bit mask based on d
	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(dv[1]));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
	m = _mm512_mask_set1_epi64(m, _mm512_cmple_epi64_mask(dv[1], onev), 0);

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually (and hence the penultimate mask bit)
	addmod104_x8(&rv[1], &rv[0], mone[1], mone[0], mone[1], mone[0], nv[1], nv[0]);

	__mmask8 protect = _mm512_cmpgt_epi64_mask(_mm512_srli_epi64(n1v[1], 49), zerov);

	// compute b^d
	if (protect)
	{
		done = _mm512_cmpeq_epu64_mask(m, zerov);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[1]);

			mask_sqrredc104_exact_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
			mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zerov);
		}
	}
	else
	{
		done = _mm512_cmpeq_epu64_mask(m, zerov);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[1]);

			mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
			mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zerov);
		}
	}

	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(dv[0]));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
	m = _mm512_mask_set1_epi64(m, _mm512_cmpge_epi64_mask(dv[1], onev), 1ULL << 51);

	if (protect)
	{
		done = _mm512_cmpeq_epu64_mask(m, zerov);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[0]);

			mask_sqrredc104_exact_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
			mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zerov);
		}
	}
	else
	{
		done = _mm512_cmpeq_epu64_mask(m, zerov);
		while (done != 0xff)
		{
			__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[0]);

			mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
			mask_dblmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], nv[1], nv[0]);

			m = _mm512_srli_epi64(m, 1);
			done = _mm512_cmpeq_epu64_mask(m, zerov);
		}
	}

	// AMM possibly needs a final correction by n
	addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

	// check current result == 1
	__mmask8 is1prp = _mm512_cmpeq_epu64_mask(rv[1], mone[1]) &
		_mm512_cmpeq_epu64_mask(rv[0], mone[0]);

	// now compute b^(2^s*d) and check for congruence to -1 as we go.
	// check while tzcnt is > 1 for all inputs or all are already not prp.
	done = is1prp;
	__mmask8 ism1prp = 0;

	submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], mone[1], mone[0], nv[1], nv[0]);

	while (done != 0xff)
	{
		tzcntv = _mm512_mask_sub_epi64(tzcntv, ~done, tzcntv, onev);

		// prp by -1 check
		ism1prp = (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
			_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

		is1prp |= ism1prp;	// stop checking it if we've found a prp criteria.
		done = (is1prp | ism1prp);

		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

		// definitely not prp by 1 check, stop checking
		done |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], zerov) &
			_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], onev));

		done |= _mm512_mask_cmple_epu64_mask(~done, tzcntv, onev);
	}

	addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

	// check current result == m-1
	ism1prp |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
		_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

	return (is1prp | ism1prp);
}

// a Miller-Rabin SPRP test on 8x 104-bit inputs using an
// independent arbitrary 52-bit base on each
uint8_t MR_sprp_104x8(uint64_t* n, uint64_t* bases)
{
	// assumes has no small factors.  assumes n >= 54 bits.
	// assumes n is a list of 8 104-bit integers (16 52-bit words)
	// in the format: 8 lo-words, 8 hi-words.
	// assume bases is a list of 8 small (single-word) bases, one for each input n.
	// do a Miller-Rabin sprp test on each using the supplied bases.
	__m512i vrho = multiplicative_inverse104_x8(n);
	__m512i mone[2];
	vec_u104_t unity;
	__m512i nv[2];
	__m512i dv[2];
	__m512i rv[2];
	__m512i bv[2];
	__m512i n1v[2];
	__m512i tv[2];
	__m512i m;
	__m512i zerov = _mm512_setzero_si512();
	__m512i onev = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	uint64_t tmp = 0;

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zerov, vrho), lo52mask);

	nv[0] = loadu64(&n[0]);
	nv[1] = loadu64(&n[8]);

	// the 128-bit division we do one at a time
	int i;
	for (i = 0; i < 8; i++)
	{
		uint128_t mod = ((uint128_t)n[i + 8] << 52) + n[i];
		uint128_t one = (uint128_t)1 << 104;
		one %= mod;

		unity.data[0][i] = (uint64_t)one & 0xfffffffffffffULL;
		unity.data[1][i] = (uint64_t)(one >> 52);
	}

	mone[0] = loadu64(unity.data[0]);
	mone[1] = loadu64(unity.data[1]);

	// get bases into Monty rep
	bv[0] = loadu64(bases);
	bv[1] = zerov;

	__m512i mpow[2];

	mpow[0] = mone[0];
	mpow[1] = mone[1];

	rv[0] = zerov;
	rv[1] = zerov;

	__mmask8 done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	while (done != 0xff)
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(onev, bv[0]);
		mask_addmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], mpow[1], mpow[0], nv[1], nv[0]);
		addmod104_x8(&mpow[1], &mpow[0], mpow[1], mpow[0], mpow[1], mpow[0], nv[1], nv[0]);

		bv[0] = _mm512_srli_epi64(bv[0], 1);
		done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	}

	bv[0] = rv[0];
	bv[1] = rv[1];

	// compute d and tzcnt
	submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], zerov, onev, nv[1], nv[0]);

	done = 0;
	dv[1] = n1v[1];
	dv[0] = n1v[0];
	__m512i tzcntv = zerov;
	while (done != 0xff)
	{
		__m512i c = _mm512_mask_slli_epi64(dv[1], ~done, dv[1], 51);
		dv[0] = _mm512_mask_srli_epi64(dv[0], ~done, dv[0], 1);
		dv[0] = _mm512_mask_or_epi64(dv[0], ~done, c, dv[0]);
		dv[1] = _mm512_mask_srli_epi64(dv[1], ~done, dv[1], 1);
		tzcntv = _mm512_mask_add_epi64(tzcntv, ~done, tzcntv, onev);
		done = done | _mm512_cmpeq_epi64_mask(_mm512_and_epi64(dv[0], onev), onev);
	}
	dv[0] = _mm512_and_epi64(dv[0], lo52mask);

	// penultimate-hi-bit mask based on d
	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(dv[1]));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
	m = _mm512_mask_set1_epi64(m, _mm512_cmple_epi64_mask(dv[1], onev), 0);

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually (and hence the penultimate mask bit)
	rv[0] = bv[0];
	rv[1] = bv[1];

	// compute b^d
	done = _mm512_cmpeq_epu64_mask(m, zerov);
	while (done != 0xff)
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[1]);

		mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
		mask_mulredc104_vec(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		m = _mm512_srli_epi64(m, 1);
		done = _mm512_cmpeq_epu64_mask(m, zerov);
	}

	m = _mm512_sub_epi64(_mm512_set1_epi64(62), _mm512_lzcnt_epi64(dv[0]));
	m = _mm512_sllv_epi64(_mm512_set1_epi64(1), m);
	m = _mm512_mask_set1_epi64(m, _mm512_cmpge_epi64_mask(dv[1], onev), 1ULL << 51);

	done = _mm512_cmpeq_epu64_mask(m, zerov);
	while (done != 0xff)
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(m, dv[0]);

		mask_sqrredc104_vec(&rv[1], &rv[0], ~done, rv[1], rv[0], nv[1], nv[0], vrho);
		mask_mulredc104_vec(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		m = _mm512_srli_epi64(m, 1);
		done = _mm512_cmpeq_epu64_mask(m, zerov);
	}

	// AMM possibly needs a final correction by n
	addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

	// check current result == 1
	__mmask8 is1prp = _mm512_cmpeq_epu64_mask(rv[1], mone[1]) &
		_mm512_cmpeq_epu64_mask(rv[0], mone[0]);

	// now compute b^(2^s*d) and check for congruence to -1 as we go.
	// check while tzcnt is > 1 for all inputs or all are already not prp.
	done = is1prp;
	__mmask8 ism1prp = 0;

	submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], mone[1], mone[0], nv[1], nv[0]);

	while (done != 0xff)
	{
		tzcntv = _mm512_mask_sub_epi64(tzcntv, ~done, tzcntv, onev);

		// prp by -1 check
		ism1prp = (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
			_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

		is1prp |= ism1prp;	// stop checking it if we've found a prp criteria.
		done = (is1prp | ism1prp);

		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		addmod104_x8(&rv[1], &rv[0], rv[1], rv[0], zerov, zerov, nv[1], nv[0]);

		// definitely not prp by 1 check, stop checking
		done |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], zerov) &
			_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], onev));

		done |= _mm512_mask_cmple_epu64_mask(~done, tzcntv, onev);
	}

	addmod104_x8(&rv[1], &rv[0], rv[1], rv[0], zerov, zerov, nv[1], nv[0]);

	// check current result == m-1
	ism1prp |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
		_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

	return (is1prp | ism1prp);
}

void modexp_104x8b(uint8_t *is1msk, uint8_t *ism1msk, uint64_t* n, uint64_t *e, uint64_t* bases)
{
	// assumes n >= 54 bits.
	// assumes n is 104-bit integer (2 52-bit words)
	// assumes e is 104-bit integer (2 52-bit words)
	// assume bases is a list of 8 small (single-word) bases.
	// compute b^e % n for each base b and compare result to 1 and -1
	__m512i vrho = _mm512_set1_epi64(multiplicative_inverse(n[0]));
	__m512i one[2];
	__m512i nv[2];
	__m512i rv[2];
	__m512i bv[2];
	__m512i n1v[2];
	__m512i m;
	__m512i zerov = _mm512_setzero_si512();
	__m512i onev = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	uint64_t d[2];
	uint64_t tmp = 0;

	// normal (negative) Montgomery mul
	vrho = _mm512_and_epi64(_mm512_sub_epi64(zerov, vrho), lo52mask);

	// Monty 1
	uint128_t mod = ((uint128_t)n[1] << 52) + n[0];
	uint128_t one128 = (uint128_t)1 << 104;
	one128 %= mod;

	one[0] = set64((uint64_t)one128 & 0xfffffffffffffULL);
	one[1] = set64((uint64_t)(one128 >> 52));

#ifdef DEBUG_THIS
	printf("n = %013lx%013lx\n", n[1], n[0]);
	printf("e = %013lx%013lx\n", e[1], e[0]);
	printf("1 = %013lx%013lx\n", (uint64_t)(one128 >> 52), (uint64_t)one128 & 0xfffffffffffffULL);
#endif

	// load N
	nv[0] = set64(n[0]);
	nv[1] = set64(n[1]);

	// and E
	d[0] = e[0];
	d[1] = e[1];

	// get bases into Monty rep
	bv[0] = loadu64(bases);
	bv[1] = zerov;

	__m512i mpow[2];

	mpow[0] = one[0];
	mpow[1] = one[1];

	rv[0] = zerov;
	rv[1] = zerov;

	__mmask8 done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	while (done != 0xff)
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(onev, bv[0]);
		mask_addmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], mpow[1], mpow[0], nv[1], nv[0]);
		addmod104_x8(&mpow[1], &mpow[0], mpow[1], mpow[0], mpow[1], mpow[0], nv[1], nv[0]);

		bv[0] = _mm512_srli_epi64(bv[0], 1);
		done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	}

	bv[0] = rv[0];
	bv[1] = rv[1];

#ifdef DEBUG_THIS
	{
		uint64_t bases0[8];
		uint64_t bases1[8];
		storeu64(bases0, bv[0]);
		storeu64(bases1, bv[1]);
		printf("bases:\n");
		int i;
		for (i = 0; i < 8; i++)
		{
			printf("%013lx%013lx\n", bases1[i], bases0[i]);
		}
	
	}
#endif
	
	// RL-52x2
	rv[0] = one[0];
	rv[1] = one[1];
	int i = 0;
	while (d[0] > 0)
	{
		if (d[0] & 1)
			mulredc104_vec(&rv[1], &rv[0], 
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d[0] >>= 1;
		i++;

		sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

	for (; (i < 52) && (d[1] > 0); i++)
	{
		sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

	while (d[1] > 0)
	{
		if (d[1] & 1)
			mulredc104_vec(&rv[1], &rv[0], 
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d[1] >>= 1;

		sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

	// AMM possibly needs a final correction by n
	chkmod104_x8(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0]);
	chkmod104_x8(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0]);

#ifdef DEBUG_THIS
	{
		uint64_t bases0[8];
		uint64_t bases1[8];
		storeu64(bases0, rv[0]);
		storeu64(bases1, rv[1]);
		printf("results:\n");
		int i;
		for (i = 0; i < 8; i++)
		{
			printf("%013lx%013lx\n", bases1[i], bases0[i]);
		}
	
	}
#endif

	// check current result == 1
	*is1msk = _mm512_cmpeq_epu64_mask(rv[1], one[1]) &
		_mm512_cmpeq_epu64_mask(rv[0], one[0]);

	submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], one[1], one[0], nv[1], nv[0]);

	// check current result == m-1
	*ism1msk = (_mm512_cmpeq_epu64_mask(rv[1], n1v[1]) &
		_mm512_cmpeq_epu64_mask(rv[0], n1v[0]));

	return;
}

// a Miller-Rabin SPRP test on 1 104-bit input using 8x
// different bases
uint8_t MR_sprp_104x8base(uint64_t* n, uint64_t* one, uint64_t* bases)
{
	// assumes has no small factors and is odd.  assumes n >= 54 bits.
	// assumes n is a 104-bit integer with two 52 bit words: [lo,hi].
	// assumes one is a 104-bit integer equal to (1 << 104) mod n.
	// assume bases is a list of 8 small (single-word) bases.
	// do a Miller-Rabin sprp test on the input to each supplied base.
	// uint128_t n128 = ((uint128_t)n[1] << 52) + (uint128_t)n[0];
	__m512i vrho = _mm512_set1_epi64(multiplicative_inverse(n[0]));
	__m512i mone[2];
	__m512i nv[2];
	__m512i dv[2];
	__m512i rv[2];
	__m512i bv[2];
	__m512i n1v[2];
	__m512i tv[2];
	__m512i zerov = _mm512_setzero_si512();
	__m512i onev = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	uint64_t tmp = 0;

	vrho = _mm512_and_epi64(_mm512_sub_epi64(zerov, vrho), lo52mask);

	nv[0] = _mm512_set1_epi64(n[0]);
	nv[1] = _mm512_set1_epi64(n[1]);

	mone[0] = _mm512_set1_epi64(one[0]);
	mone[1] = _mm512_set1_epi64(one[1]);

	// get bases into Monty rep
	bv[0] = loadu64(bases);
	bv[1] = zerov;

	__m512i mpow[2];

	mpow[0] = mone[0];
	mpow[1] = mone[1];

	rv[0] = zerov;
	rv[1] = zerov;

	__mmask8 done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	while (done != 0xff)
	{
		__mmask8 bitcmp = _mm512_test_epi64_mask(onev, bv[0]);
		mask_addmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], mpow[1], mpow[0], nv[1], nv[0]);
		addmod104_x8(&mpow[1], &mpow[0], mpow[1], mpow[0], mpow[1], mpow[0], nv[1], nv[0]);

		bv[0] = _mm512_srli_epi64(bv[0], 1);
		done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	}

	bv[0] = rv[0];
	bv[1] = rv[1];

	// compute d and tzcnt
	uint64_t d[2];
	d[1] = n[1];
	d[0] = n[0] - 1;			// n odd, so this won't carry
	int ntz = my_ctz104(d[0], d[1]);
	n1v[0] = _mm512_set1_epi64(d[0]);
	n1v[1] = _mm512_set1_epi64(d[1]);
	if (ntz < 52)
	{
		uint64_t shift = d[1] & ((1ULL << ntz) - 1);
		d[0] = (d[0] >> ntz) + (shift << (52 - ntz));
		d[1] >>= ntz;
	}
	else
	{
		d[0] = d[1];
		d[1] = 0;
		d[0] >>= (ntz - 52);
	}


#if 0
	// LR
	uint64_t m;
	if (d[1] <= 1) m = 0;
	else m = (1ull << (62 - my_clz52(d[1])));

	while (m > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		if (m & d[1])
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);
		m >>= 1;
	}

	if (d[1] == 0) m = (1ull << (62 - my_clz52(d[0])));
	else m = (1ull << 51);

	while (m > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		if (m & d[0])
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);
		m >>= 1;
	}
#elif 0
	// RL-104
	uint128_t d128 = ((uint128_t)d[1] << 52) + (uint128_t)d[0];
	rv[0] = mone[0];
	rv[1] = mone[1];
	while (d128 > 0)
	{
		if (d128 & 1)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d128 >>= 1;

		if (d128)
			sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}
#elif 0
	// LR-kary
	uint64_t g[256];
	// 
	// precomputation
	rv[0] = bv[0];
	rv[1] = bv[1];

	int i;
	storeu64(&g[0 * 8], mone[0]);		// g0 = 1
	storeu64(&g[1 * 8], mone[1]);		// g0 = 1
	storeu64(&g[2 * 8], rv[0]);			// g1 = g 
	storeu64(&g[3 * 8], rv[1]);			// g1 = g 
	sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
	storeu64(&g[4 * 8], rv[0]);			// g2 = g^2
	storeu64(&g[5 * 8], rv[1]);			// g2 = g^2
	for (i = 3; i < 16; i++)
	{
		mask_mulredc104_vec(&rv[1], &rv[0], 0xff, rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);
		storeu64(&g[(i * 2) * 8], rv[0]);			// gi = g^i
		storeu64(&g[(i * 2 + 1) * 8], rv[1]);		// gi = g^i
	}

	uint128_t d128 = ((uint128_t)d[1] << 52) + (uint128_t)d[0];
	rv[0] = mone[0];
	rv[1] = mone[1];
	int lz = my_clz104(d[0], d[1]);
	int msb = 112 - lz;
	int m;

	m = (d128 >> msb) & 0xf;
	rv[0] = loadu64(&g[(2 * m) * 8]);
	rv[1] = loadu64(&g[(2 * m + 1) * 8]);
	msb -= 4;

	while (msb > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		m = (d128 >> msb) & 0xf;

		if (m > 0)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				loadu64(&g[(2 * m + 1) * 8]), loadu64(&g[(2 * m) * 8]),
				rv[1], rv[0], nv[1], nv[0], vrho);

		msb -= 4;
	}

	msb += 4;
	m = (int)d128 & ((1 << msb) - 1);
	while (msb > 0)
	{
		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		msb--;
	}

	mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
		loadu64(&g[(2 * m + 1) * 8]), loadu64(&g[(2 * m) * 8]),
		rv[1], rv[0], nv[1], nv[0], vrho);
#else
	// RL-52x2
	rv[0] = mone[0];
	rv[1] = mone[1];
	int i = 0;
	while (d[0] > 0)
	{
		if (d[0] & 1)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d[0] >>= 1;
		i++;

		sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

	for (; (i < 52) && (d[1] > 0); i++)
	{
		sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

	while (d[1] > 0)
	{
		if (d[1] & 1)
			mask_mulredc104_vec(&rv[1], &rv[0], 0xff,
				rv[1], rv[0], bv[1], bv[0], nv[1], nv[0], vrho);

		d[1] >>= 1;

		if (d[1])
			sqrredc104_vec(&bv[1], &bv[0], bv[1], bv[0], nv[1], nv[0], vrho);
	}

#endif

	// AMM possibly needs a final correction by n
	addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

	// check current result == 1
	__mmask8 is1prp = _mm512_cmpeq_epu64_mask(rv[1], mone[1]) &
		_mm512_cmpeq_epu64_mask(rv[0], mone[0]);

	// now compute b^(2^s*d) and check for congruence to -1 as we go.
	// check while tzcnt is > 1 for all inputs or all are already not prp.
	done = is1prp;
	__mmask8 ism1prp = 0;

	submod104_x8(&n1v[1], &n1v[0], nv[1], nv[0], mone[1], mone[0], nv[1], nv[0]);

	while ((done != 0xff) && (ntz > 0))
	{
		ntz--;

		// prp by -1 check
		ism1prp = (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
			_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

		is1prp |= ism1prp;	// stop checking it if we've found a prp criteria.
		done = (is1prp | ism1prp);

		sqrredc104_vec(&rv[1], &rv[0], rv[1], rv[0], nv[1], nv[0], vrho);
		addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

		// definitely not prp by 1 check, stop checking
		done |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], zerov) &
			_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], onev));
	}

	addmod104_x8(&rv[1], &rv[0], zerov, zerov, rv[1], rv[0], nv[1], nv[0]);

	// check current result == m-1
	ism1prp |= (_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[1], n1v[1]) &
		_mm512_mask_cmpeq_epu64_mask(~is1prp, rv[0], n1v[0]));

	return (is1prp | ism1prp);
}

#endif

