/* lasieve-prepn.c
  By Jens Franke.
  6/13/04: Hacked up for use in GGNFS by Chris Monico.
  9/30/22: Vector AVX512 code contributed by Ben Buhrow

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/




/*4:*/
#line 19 "lasieve-prepn.w"

#include <sys/types.h> 
#include <math.h> 
#include <stdlib.h>
#include "asm/siever-config.h"
#include "recurrence6.h"
#include "asm/32bit.h"
#include <stdint.h>
#include <stdio.h>
#include "if.h"
#include <immintrin.h>
#include "avx512_aux.h"

#include <fenv.h>
#pragma STDC FENV_ACCESS ON

#ifdef _MSC_VER
// so that I can read the code in MSVC without it being grayed out.
// It will not build in Visual studio.
#define AVX512_LASIEVE_SETUP
#endif

#ifdef AVX512_LASIEVE_SETUP

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
#define ROUNDING_MODE (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)
#else
#define ROUNDING_MODE _MM_FROUND_CUR_DIRECTION
#endif

__m512i rem_epu32_x16(__m512i n, __m512i d)
{
	// DANGER: I haven't proven this works for every possible input.
	__m512 d1ps = _mm512_cvtepu32_ps(d);
	__m512 n1ps = _mm512_cvtepu32_ps(n);
	__m512i q2, q, r;

	//n1ps = _mm512_div_ps(n1ps, d1ps);
	//q = _mm512_cvt_roundps_epu32(n1ps, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	n1ps = _mm512_div_round_ps(n1ps, d1ps, ROUNDING_MODE);
	q = _mm512_cvttps_epu32(n1ps);

	__m512i qd = _mm512_mullo_epi32(q, d);
	r = _mm512_sub_epi32(n, qd);

	// fix q too big by a little, with special check for
	// numerators close to 2^32 and denominators close to 1
	__mmask16 err = _mm512_cmpgt_epu32_mask(r, n) |
		(_mm512_cmpgt_epu32_mask(r, d) & _mm512_cmplt_epu32_mask(
			_mm512_sub_epi32(_mm512_set1_epi32(0), r), _mm512_set1_epi32(1024)));
	if (err)
	{
		n1ps = _mm512_cvtepu32_ps(_mm512_sub_epi32(_mm512_set1_epi32(0), r));
		
		//n1ps = _mm512_div_ps(n1ps, d1ps);
		//q2 = _mm512_add_epi32(_mm512_set1_epi32(1), tmp_cvtround_ps_epu32(n1ps,
		//    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

		n1ps = _mm512_div_round_ps(n1ps, d1ps, ROUNDING_MODE);
		q2 = _mm512_add_epi32(_mm512_set1_epi32(1), _mm512_cvttps_epu32(n1ps));

		q = _mm512_mask_sub_epi32(q, err, q, q2);
		r = _mm512_mask_add_epi32(r, err, r, _mm512_mullo_epi32(q2, d));
	}

	// fix q too small by a little bit
	err = _mm512_cmpge_epu32_mask(r, d);
	if (err)
	{
		n1ps = _mm512_cvtepu32_ps(r);

		//n1ps = _mm512_div_ps(n1ps, d1ps);
		//q2 = tmp_cvtround_ps_epu32(n1ps,
		//    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
		n1ps = _mm512_div_round_ps(n1ps, d1ps, ROUNDING_MODE);
		q2 = _mm512_cvttps_epu32(n1ps);

		q = _mm512_mask_add_epi32(q, err, q, q2);
		r = _mm512_mask_sub_epi32(r, err, r, _mm512_mullo_epi32(q2, d));
	}

	return r;
}

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
	
	__m512i qd = _mm512_mullo_epi64(q, d);
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
		r = _mm512_mask_add_epi64(r, err, r, _mm512_mullo_epi64(q2, d));
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
		r = _mm512_mask_sub_epi64(r, err, r, _mm512_mullo_epi64(q2, d));
	}

	return r;
}

void divrem_mask_epu32_x16(__m512i* q, __m512i* r, __mmask16 msk, __m512i n, __m512i d)
{
	// DANGER: I haven't proven this works for every possible input.
	__m512 d1ps = _mm512_cvtepu32_ps(d);
	__m512 n1ps = _mm512_cvtepu32_ps(n);
	__m512i q2;

	// n1ps = _mm512_div_ps(n1ps, d1ps);
	// *q = _mm512_mask_cvt_roundps_epu32(*q, msk, n1ps, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

	n1ps = _mm512_div_round_ps(n1ps, d1ps, ROUNDING_MODE);
	*q = _mm512_mask_cvttps_epu32(*q, msk, n1ps);

	__m512i qd = _mm512_mullo_epi32(*q, d);
	*r = _mm512_mask_sub_epi32(*r, msk, n, qd);

	// fix q too big by a little, with special check for
	// numerators close to 2^32 and denominators close to 1
	__mmask16 err = _mm512_mask_cmpgt_epu32_mask(msk, *r, n) |
		(_mm512_mask_cmpgt_epu32_mask(msk, *r, d) & _mm512_mask_cmplt_epu32_mask(msk,
			_mm512_sub_epi32(_mm512_set1_epi32(0), *r), _mm512_set1_epi32(1024)));
	if (err)
	{
		n1ps = _mm512_cvtepu32_ps(_mm512_sub_epi32(_mm512_set1_epi32(0), *r));

		//n1ps = _mm512_div_ps(n1ps, d1ps);
		//q2 = _mm512_add_epi32(_mm512_set1_epi32(1), tmp_cvtround_ps_epu32(n1ps,
		//    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)));

		n1ps = _mm512_div_round_ps(n1ps, d1ps, ROUNDING_MODE);
		q2 = _mm512_add_epi32(_mm512_set1_epi32(1), _mm512_cvttps_epu32(n1ps));

		*q = _mm512_mask_sub_epi32(*q, err, *q, q2);
		*r = _mm512_mask_add_epi32(*r, err, *r, _mm512_mullo_epi32(q2, d));
	}

	// fix q too small by a little bit
	err = _mm512_mask_cmpge_epu32_mask(msk, *r, d);
	if (err)
	{
		n1ps = _mm512_cvtepu32_ps(*r);

		//n1ps = _mm512_div_ps(n1ps, d1ps);
		//q2 = tmp_cvtround_ps_epu32(n1ps,
		//    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
		n1ps = _mm512_div_round_ps(n1ps, d1ps, ROUNDING_MODE);
		q2 = _mm512_cvttps_epu32(n1ps);

		*q = _mm512_mask_add_epi32(*q, err, *q, q2);
		*r = _mm512_mask_sub_epi32(*r, err, *r, _mm512_mullo_epi32(q2, d));
	}

	return;
}

__m512i modmul32_16(__m512i z, __mmask16 ndmsk, __m512i x, __m512i y, __m512i p)
{
	// multiply the 16-element 32-bit vectors a and b to produce two 8-element
	// 64-bit vector products e64 and o64, where e64 is the even elements
	// of a*b and o64 is the odd elements of a*b
	//__m512i t1 = _mm512_shuffle_epi32(a, 0xB1);
	//__m512i t2 = _mm512_shuffle_epi32(b, 0xB1);

	//_mm512_shuffle_epi32(a, 0xB1);
	//_mm512_shuffle_epi32(b, 0xB1);
	__m512i e = _mm512_mul_epu32(x, y);
	__m512i o = _mm512_mul_epu32(_mm512_shuffle_epi32(x, 0xB1), _mm512_shuffle_epi32(y, 0xB1));

	//e = _mm512_rem_epu64(e, _mm512_and_epi64(p, _mm512_set1_epi64(0x00000000ffffffff)));
	//o = _mm512_rem_epu64(o, _mm512_and_epi64(_mm512_shuffle_epi32(p, 0xB1), _mm512_set1_epi64(0x00000000ffffffff)));

	e = rem_epu64_x8(e, _mm512_and_epi64(p, _mm512_set1_epi64(0x00000000ffffffff)));
	o = rem_epu64_x8(o, _mm512_and_epi64(_mm512_shuffle_epi32(p, 0xB1), _mm512_set1_epi64(0x00000000ffffffff)));

	x = _mm512_or_epi64(e, _mm512_shuffle_epi32(o, 0xB1));
	return _mm512_mask_blend_epi32(ndmsk, z, x);
}

__m512i modinv_16(__m512i y, __mmask16 ndmsk, __m512i x, __m512i p)
{
	__m512i ps1, ps2, parity, dividend, divisor, rem, q, t;
	__mmask16 dmsk = 0, lmsk, inmsk = ndmsk;

	__m512i zero = _mm512_setzero_epi32();
	__m512i one = _mm512_set1_epi32(1);
	__m512i mask32 = _mm512_set1_epi64(0xffffffff);

	q = one;
	rem = x;
	dividend = p;
	divisor = rem;
	ps1 = q;

	ps2 = zero;
	parity = zero;

	ndmsk = ndmsk & _mm512_cmpgt_epu32_mask(divisor, one);
	dmsk = ~ndmsk;

	while (ndmsk > 0) {

		rem = _mm512_mask_sub_epi32(rem, ndmsk, dividend, divisor);
		t = _mm512_sub_epi32(rem, divisor);
		lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);

		if (lmsk > 0) {

			q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
			rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

			t = _mm512_sub_epi32(t, divisor);
			lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);


			if (lmsk > 0) {

				q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
				rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

				t = _mm512_sub_epi32(t, divisor);
				lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);


				if (lmsk > 0) {
					q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
					rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

					t = _mm512_sub_epi32(t, divisor);
					lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);


					if (lmsk > 0) {
						q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
						rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

						t = _mm512_sub_epi32(t, divisor);
						lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);


						if (lmsk > 0) {
							q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
							rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

							t = _mm512_sub_epi32(t, divisor);
							lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);


							if (lmsk > 0) {
								q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
								rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

								t = _mm512_sub_epi32(t, divisor);
								lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);


								if (lmsk > 0) {
									q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
									rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

									t = _mm512_sub_epi32(t, divisor);
									lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);

									if (lmsk > 0) {

										q = _mm512_mask_add_epi32(q, lmsk, q, ps1);
										rem = _mm512_mask_add_epi32(rem, lmsk, t, zero);

										lmsk = _mm512_mask_cmpge_epu32_mask(ndmsk, rem, divisor);


										if (lmsk > 0) {

											divrem_mask_epu32_x16(&q, &rem, lmsk, dividend, divisor);
											q = _mm512_mask_mullo_epi32(q, lmsk, q, ps1);

										}
									}
								}
							}
						}
					}
				}
			}
		}

		q = _mm512_mask_add_epi32(q, ndmsk, q, ps2);
		parity = _mm512_mask_andnot_epi32(parity, ndmsk, parity, _mm512_set1_epi32(0xffffffff));
		dividend = _mm512_mask_mov_epi32(dividend, ndmsk, divisor);
		divisor = _mm512_mask_mov_epi32(divisor, ndmsk, rem);
		ps2 = _mm512_mask_mov_epi32(ps2, ndmsk, ps1);
		ps1 = _mm512_mask_mov_epi32(ps1, ndmsk, q);
		ndmsk = _mm512_cmpgt_epu32_mask(divisor, one);
		dmsk = ~ndmsk;

	}

	ps1 = _mm512_mask_blend_epi32(_mm512_cmpeq_epi32_mask(parity, zero), _mm512_sub_epi32(p, ps1), ps1);
	return _mm512_mask_blend_epi32(inmsk, y, ps1);
}

__m512i modsub32_16(__m512i z, __mmask16 ndmsk, __m512i x, __m512i y, __m512i p)
{
	z = _mm512_mask_sub_epi32(z, ndmsk, x, y);
	__mmask16 m = _mm512_cmpgt_epu32_mask(y, x);
	z = _mm512_mask_add_epi32(z, ndmsk & m, z, p);
	return z;
}

#endif


void
lasieve_setup(u32_t* FB, u32_t* proots, u32_t fbsz, i32_t a0, i32_t a1, i32_t b0, i32_t b1,
	u32_t* ri_ptr, u32_t FBsize)
{
	u32_t b0_ul, b1_ul, absa0, absa1;

	if (b0 < 0 || b1 < 0)
		Schlendrian("lasieve_setup called with negative 2-nd coordinate (%d,%d)\n",
			b0, b1);
	if (fbsz <= 0)return;

#if !defined(AVX512_LASIEVE_SETUP) && defined( HAVE_ASM_LASIEVE_SETUP)
	/*7:*/
#line 155 "lasieve-prepn.w"

	if (FB[fbsz - 1] < FLOAT_SETUP_BOUND1 &&
		fabs(a0) + FB[fbsz - 1] * b0 < FLOAT_SETUP_BOUND2 &&
		fabs(a1) + FB[fbsz - 1] * b1 < FLOAT_SETUP_BOUND2) {
		asm_lasieve_setup(FB, proots, fbsz, a0, a1, b0, b1, ri_ptr);
		return;
	}/*:7*/
#line 39 "lasieve-prepn.w"

#endif

	b0_ul = (u32_t)b0;
	b1_ul = (u32_t)b1;
	if (a0 >= 0) {
		absa0 = (u32_t)a0;
#define A0MOD0(p) (absa0%p)

#define A0MOD1(p) absa0

		if (a1 >= 0) {
			absa1 = (u32_t)a1;



#define A1MOD0(p) (absa1%p)
#define A1MOD1(p) absa1
			/*5:*/
#line 95 "lasieve-prepn.w"

			{
				u32_t fbi, fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
				fbp_bound = absa1 < absa0 ? absa0 : absa1;
				if (fbp_bound < b0_ul)fbp_bound = b0_ul;
				if (fbp_bound < b1_ul)fbp_bound = b1_ul;

#ifdef AVX512_LASIEVE_SETUP

				
				for (fbi = 0; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (_mm512_cmpgt_epu32_mask(m32, _mm512_set1_epi32(fbp_bound)) > 0)
						break;

					__m512i vb0 = rem_epu32_x16(_mm512_set1_epi32(b0_ul), m32);
					__m512i vb1 = rem_epu32_x16(_mm512_set1_epi32(b1_ul), m32);
					__m512i vabsa0 = rem_epu32_x16(_mm512_set1_epi32(absa0), m32);
					__m512i vabsa1 = rem_epu32_x16(_mm512_set1_epi32(absa1), m32);
					__mmask16 m1 = _mm512_cmpeq_epu32_mask(pr, m32);

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, vabsa0, t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, vabsa1, m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);

				}

				_mm256_zeroupper();

				for (; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}

#else
				for (fbi = 0; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#endif

				/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)



#ifdef AVX512_LASIEVE_SETUP
				
				for (; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__mmask16 m1 = _mm512_cmpeq_epu32_mask(pr, m32);
					__m512i vb0 = _mm512_set1_epi32(b0_ul);
					__m512i vb1 = _mm512_set1_epi32(b1_ul);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, _mm512_set1_epi32(absa0), t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, _mm512_set1_epi32(absa1), m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);

				}

				_mm256_zeroupper();

				for (; fbi < fbsz; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = b0_ul;
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), b1_ul);
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(absa0, modmul32(proots[fbi], b0_ul));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], b1_ul), absa1));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}

#else
				for (; fbi < fbsz; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						//printf("rare case? proots = %u, modulo32 = %u, b0_ul = %u, ",
						//	proots[fbi], modulo32, b0_ul);
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
						//printf("x = %u\n", x);
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#endif

				/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
			}

			/*:5*/
#line 55 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
		}
		else {
			u32_t aux;

			absa1 = (u32_t)(-a1);
#define A1MOD0(p) ((aux= absa1%p)> 0 ? p-aux : 0 )
#define A1MOD1(p) (p-absa1)
			/*5:*/
#line 95 "lasieve-prepn.w"

			{
				u32_t fbi, fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
				fbp_bound = absa1 < absa0 ? absa0 : absa1;
				if (fbp_bound < b0_ul)fbp_bound = b0_ul;
				if (fbp_bound < b1_ul)fbp_bound = b1_ul;

#ifdef AVX512_LASIEVE_SETUP
				for (fbi = 0; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (_mm512_cmpgt_epu32_mask(m32, _mm512_set1_epi32(fbp_bound)) > 0)
						break;

					__m512i vb0 = rem_epu32_x16(_mm512_set1_epi32(b0_ul), m32);
					__m512i vb1 = rem_epu32_x16(_mm512_set1_epi32(b1_ul), m32);
					__m512i vabsa0 = rem_epu32_x16(_mm512_set1_epi32(absa0), m32);
					__m512i vabsa1 = rem_epu32_x16(_mm512_set1_epi32(absa1), m32); 
					__mmask16 m1 = _mm512_cmpgt_epu32_mask(vabsa1, zero);

					vabsa1 = _mm512_mask_sub_epi32(vabsa1, m1, m32, vabsa1);
					m1 = _mm512_cmpeq_epu32_mask(pr, m32);

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, vabsa0, t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, vabsa1, m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);

				}

				_mm256_zeroupper();

				for (; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}

#else
				for (fbi = 0; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#endif

				/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)

#ifdef AVX512_LASIEVE_SETUP
				for (; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__mmask16 m1 = _mm512_cmpeq_epu32_mask(pr, m32);
					__m512i vb0 = _mm512_set1_epi32(b0_ul);
					__m512i vb1 = _mm512_set1_epi32(b1_ul);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, _mm512_set1_epi32(absa0), t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, _mm512_sub_epi32(m32, _mm512_set1_epi32(absa1)), m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);

				}

				_mm256_zeroupper();

				for (; fbi < fbsz; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = b0_ul;
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), b1_ul);
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(absa0, modmul32(proots[fbi], b0_ul));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], b1_ul), (modulo32 - absa1)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}

#else
				for (; fbi < fbsz; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}

#endif

				/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
			}

			/*:5*/
#line 64 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
#undef A0MOD0
#undef A0MOD1
		}
	}
	else {
		absa0 = (u32_t)(-a0);
#define A0MOD0(p) ((aux= absa0%p)> 0 ? p-aux : 0 )
#define A0MOD1(p) (p-absa0)
		if (a1 >= 0) {
			u32_t aux;
			absa1 = (u32_t)a1;
#define A1MOD0(p) (absa1%p)
#define A1MOD1(p) absa1
			/*5:*/
#line 95 "lasieve-prepn.w"

			{
				u32_t fbi, fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
				fbp_bound = absa1 < absa0 ? absa0 : absa1;
				if (fbp_bound < b0_ul)fbp_bound = b0_ul;
				if (fbp_bound < b1_ul)fbp_bound = b1_ul;

#ifdef AVX512_LASIEVE_SETUP

				for (fbi = 0; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (_mm512_cmpgt_epu32_mask(m32, _mm512_set1_epi32(fbp_bound)) > 0)
						break;

					__m512i vb0 = rem_epu32_x16(_mm512_set1_epi32(b0_ul), m32);
					__m512i vb1 = rem_epu32_x16(_mm512_set1_epi32(b1_ul), m32);
					__m512i vabsa0 = rem_epu32_x16(_mm512_set1_epi32(absa0), m32); 
					__m512i vabsa1 = rem_epu32_x16(_mm512_set1_epi32(absa1), m32);
					__mmask16 m1 = _mm512_cmpgt_epu32_mask(vabsa0, zero);

					vabsa0 = _mm512_mask_sub_epi32(vabsa0, m1, m32, vabsa0);

					m1 = _mm512_cmpeq_epu32_mask(pr, m32);

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, vabsa0, t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, vabsa1, m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);

				}

				_mm256_zeroupper();

				for (; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}

#else
				for (fbi = 0; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#endif

				/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)


#ifdef AVX512_LASIEVE_SETUP
				for (; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__mmask16 m1 = _mm512_cmpeq_epu32_mask(pr, m32);
					__m512i vb0 = _mm512_set1_epi32(b0_ul);
					__m512i vb1 = _mm512_set1_epi32(b1_ul);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, _mm512_sub_epi32(m32, _mm512_set1_epi32(absa0)), t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, _mm512_set1_epi32(absa1), m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);


				}

				_mm256_zeroupper();

				for (; fbi < fbsz; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = b0_ul;
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), b1_ul);
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(modulo32 - absa0, modmul32(proots[fbi], b0_ul));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], b1_ul), absa1));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#else

				for (; fbi < fbsz; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#endif

				/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
			}

			/*:5*/
#line 79 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
		}
		else {
			u32_t aux;
			absa1 = (u32_t)(-a1);
#define A1MOD0(p) ((aux= absa1%p)> 0 ? p-aux : 0 )
#define A1MOD1(p) (p-absa1)
			/*5:*/
#line 95 "lasieve-prepn.w"

			{
				u32_t fbi, fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
				fbp_bound = absa1 < absa0 ? absa0 : absa1;
				if (fbp_bound < b0_ul)fbp_bound = b0_ul;
				if (fbp_bound < b1_ul)fbp_bound = b1_ul;

#ifdef AVX512_LASIEVE_SETUP
				for (fbi = 0; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (_mm512_cmpgt_epu32_mask(m32, _mm512_set1_epi32(fbp_bound)) > 0)
						break;

					__m512i vb0 = rem_epu32_x16(_mm512_set1_epi32(b0_ul), m32);
					__m512i vb1 = rem_epu32_x16(_mm512_set1_epi32(b1_ul), m32);
					__m512i vabsa0 = rem_epu32_x16(_mm512_set1_epi32(absa0), m32);
					__m512i vabsa1 = rem_epu32_x16(_mm512_set1_epi32(absa1), m32);

					__mmask16 m1 = _mm512_cmpgt_epu32_mask(vabsa0, zero);
					__mmask16 m2 = _mm512_cmpgt_epu32_mask(vabsa1, zero);

					vabsa0 = _mm512_mask_sub_epi32(vabsa0, m1, m32, vabsa0);
					vabsa1 = _mm512_mask_sub_epi32(vabsa1, m2, m32, vabsa1);

					m1 = _mm512_cmpeq_epu32_mask(pr, m32);

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, vabsa0, t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, vabsa1, m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);

				}

				_mm256_zeroupper();

				for (; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}

#else
				for (fbi = 0; fbi < fbsz && FB[fbi] <= fbp_bound; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#endif


				/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)


#ifdef AVX512_LASIEVE_SETUP
//u32_t count = fbsz;
				for (; (fbi < (fbsz - 16)) && (fbsz >= 16); fbi += 16)/*6:*/
				{
					__m512i x;
					__m512i m32 = _mm512_loadu_si512(&FB[fbi]);
					__m512i pr = _mm512_loadu_si512(&proots[fbi]);
					__mmask16 m1 = _mm512_cmpeq_epu32_mask(pr, m32);
					__m512i vb0 = _mm512_set1_epi32(b0_ul);
					__m512i vb1 = _mm512_set1_epi32(b1_ul);
					__m512i zero = _mm512_setzero_epi32();
					uint32_t xm[16];

					if (m1 > 0)
					{
						int ii;
						for (ii = 0; ii < 16; ii++)
						{
							if (m1 & (1 << ii))
							{
								modulo32 = FB[fbi + ii];
								u32_t x32 = B0MOD(modulo32);
								if (x32 == 0)x32 = modulo32;
								else {
									x32 = modmul32(modinv32(x32), B1MOD(modulo32));
									if (x32 > 0)x32 = modulo32 - x32;
								}
								xm[ii] = x32;
							}
						}
						x = _mm512_loadu_si512(xm);
					}

					// else
					m1 = ~m1;
					{
						__m512i t = zero;
						t = modmul32_16(t, m1, pr, vb0, m32);
						x = modsub32_16(x, m1, _mm512_sub_epi32(m32, _mm512_set1_epi32(absa0)), t, m32);
						__mmask16 m2 = _mm512_cmpgt_epu32_mask(x, zero);
						m1 = m1 & m2;
						t = modmul32_16(t, m1, pr, vb1, m32);
						t = modsub32_16(t, m1, t, _mm512_sub_epi32(m32, _mm512_set1_epi32(absa1)), m32);
						x = modinv_16(x, m1, x, m32);
						x = modmul32_16(x, m1, x, t, m32);
						x = _mm512_mask_mov_epi32(x, m1 & ~m2, m32);
					}

					ri_ptr += get_recurrence_info_16(ri_ptr, m32, x, FBsize);

				}

				_mm256_zeroupper();

				for (; fbi < fbsz; fbi++)
				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = b0_ul;
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), b1_ul);
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(modulo32 - absa0, modmul32(proots[fbi], b0_ul));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], b1_ul), modulo32 - absa1));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#else


				for (; fbi < fbsz; fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

				{
					u32_t x;
					modulo32 = FB[fbi];

					if (proots[fbi] == modulo32) {
						x = B0MOD(modulo32);
						if (x == 0)x = modulo32;
						else {
							x = modmul32(modinv32(x), B1MOD(modulo32));
							if (x > 0)x = modulo32 - x;
						}
					}
					else {
						x = modsub32(A0MOD(modulo32), modmul32(proots[fbi], B0MOD(modulo32)));
						if (x != 0) {
							x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[fbi], B1MOD(modulo32)),
								A1MOD(modulo32)));
						}
						else {

							x = FB[fbi];
						}
					}
					ri_ptr += get_recurrence_info(ri_ptr, FB[fbi], x, FBsize);
				}
#endif

				/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
			}

			/*:5*/
#line 87 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
		}
	}

}

/*:4*/
