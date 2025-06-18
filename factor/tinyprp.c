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
#include "monty.h"
#include "ytools.h"

#ifdef GMP_CHECK
#include "gmp.h"
#endif


int fermat_prp_64x1(uint64_t n)
{
	// assumes has no small factors.
	// do a base-2 fermat prp test using LR binexp.

	uint64_t rho = multiplicative_inverse(n);
	uint64_t unityval = ((uint64_t)0 - n) % n;  // unityval == R  (mod n)
	uint64_t m = 1ULL << (62 - __lzcnt64(n));   // set a mask at the leading bit - 2
	uint64_t r = unityval;
	uint64_t e = n - 1;

	// we know the first bit is set and the first squaring is of unity,
	// so we can do the first iteration manually with no squaring.
	r = addmod(r, r, n);

	while (m > 0)
	{
		r = sqrredc(r, n, rho);
		if (e & m) r = addmod(r, r, n);
		m >>= 1;
	}
	return (r == unityval);
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

static uint128_t my_random(void)
{
	// based on linear congruential generator, period = 2^128
	static uint128_t seed = ((uint128_t)0x123456789ull << 92) + ((uint128_t)0xabcdef << 36) + 0x987654321ull;
	seed = seed * 137 + 13;
	// shuffle
	uint128_t x = seed ^ (seed >> 17) ^ (seed << 13);
	return x;
}

static inline uint64_t my_rdtsc(void)
{
#if defined(__x86_64__)
	// supported by GCC and Clang for x86 platform
	return _rdtsc();
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

	// test of fermat_prp_64x1 on random 6k+1 inputs
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

	// test of fermat_prp_128x1 on random 6k+1 inputs
	for (bits = 118; bits <= 0; bits += 1)
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

	// test of fermat_prp_52x8 on random 6k+1 inputs
	for (bits = 20; bits <= 0; bits += 1)
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
	for (bits = 20; bits <= 0; bits += 1)
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
	for (bits = 50; bits <= 0; bits += 1)
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
	for (bits = 50; bits <= 0; bits += 1)
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
	for (bits = 50; bits <= 0; bits += 1)
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
	// do a base-2 MR sprp test on each using LR binexp.
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

	submod104_x8(&n1v[1], &n1v[0], zerov, zerov, mone[1], mone[0], nv[1], nv[0]);

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

	rv[0] = mone[0];
	rv[1] = mone[1];

	bv[0] = _mm512_srli_epi64(bv[0], 1);
	__mmask8 done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	while (done != 0xff)
	{
		addmod104_x8(&mpow[1], &mpow[0], mpow[1], mpow[0], mpow[1], mpow[0], nv[1], nv[0]);
		__mmask8 bitcmp = _mm512_test_epi64_mask(onev, bv[0]);
		mask_addmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], mpow[1], mpow[0], nv[1], nv[0]);

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

	submod104_x8(&n1v[1], &n1v[0], zerov, zerov, mone[1], mone[0], nv[1], nv[0]);

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

	rv[0] = mone[0];
	rv[1] = mone[1];

	bv[0] = _mm512_srli_epi64(bv[0], 1);

	__mmask8 done = _mm512_cmpeq_epi64_mask(bv[0], zerov);
	while (done != 0xff)
	{
		addmod104_x8(&mpow[1], &mpow[0], mpow[1], mpow[0], mpow[1], mpow[0], nv[1], nv[0]);
		__mmask8 bitcmp = _mm512_test_epi64_mask(onev, bv[0]);
		mask_addmod104_x8(&rv[1], &rv[0], (~done) & bitcmp, rv[1], rv[0], mpow[1], mpow[0], nv[1], nv[0]);

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

	submod104_x8(&n1v[1], &n1v[0], zerov, zerov, mone[1], mone[0], nv[1], nv[0]);

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

