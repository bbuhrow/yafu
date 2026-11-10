#!/usr/bin/env python3
"""
Write tinyecm_ifma104.c (to the directory given): a copy of tinyecm.c whose
tecm_mulredc104_x8 does every 52x52 product with madd52lo / madd52hi (IFMA)
instead of the FP fmadd_round pairs.  The original is kept as
tecm_mulredc104_x8_fp.  Used by tiny_bench104.c -DIFMA104_SRC.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'factor', 'gmp-ecm', 'tinyecm.c')

NEW = r'''
#define LO52(a, b, c) _mm512_madd52lo_epu64((a), (b), (c))
#define HI52(a, b, c) _mm512_madd52hi_epu64((a), (b), (c))
__inline static void tecm_mulredc104_x8(vec_u104_t* p, vec_u104_t* x, vec_u104_t* y, vec_u104_t* N, vec_u104_t* invN)
{
	__m512i z = _mm512_setzero_si512();
	__m512i lomask52 = _mm512_set1_epi64(0x000fffffffffffffULL);
	__m512i x0 = _mm512_load_si512(x->data[0]);
	__m512i x1 = _mm512_load_si512(x->data[1]);
	__m512i y0 = _mm512_load_si512(y->data[0]);
	__m512i y1 = _mm512_load_si512(y->data[1]);
	__m512i invN0 = _mm512_load_si512(invN->data[0]);
	__m512i invN1 = _mm512_load_si512(invN->data[1]);
	__m512i N0 = _mm512_load_si512(N->data[0]);
	__m512i N1 = _mm512_load_si512(N->data[1]);

	/* T = x * y */
	__m512i T0 = LO52(z, x0, y0);
	__m512i T1 = LO52(LO52(HI52(z, x0, y0), x0, y1), x1, y0);
	__m512i T2 = LO52(HI52(HI52(z, x0, y1), x1, y0), x1, y1);
	__m512i T3 = HI52(z, x1, y1);
	carryprop(T1, T2, lomask52);
	carryprop(T2, T3, lomask52);

	/* m = T mod 2^104 * invN mod 2^104 */
	__m512i m0 = LO52(z, T0, invN0);
	__m512i m1 = LO52(LO52(HI52(z, T0, invN0), T0, invN1), T1, invN0);
	m1 = _mm512_and_si512(m1, lomask52);

	/* m * N; only the high half and the carry out of the low half matter */
	__m512i mN1 = LO52(LO52(HI52(z, m0, N0), m0, N1), m1, N0);
	__m512i mN2 = LO52(HI52(HI52(z, m0, N1), m1, N0), m1, N1);
	__m512i mN3 = HI52(z, m1, N1);
	carryprop(mN1, mN2, lomask52);
	carryprop(mN2, mN3, lomask52);

	/* (T - mN) >> 104, + N if negative */
	__mmask8 msk = _mm512_cmplt_epu64_mask(T3, mN3);
	msk |= _mm512_cmpeq_epu64_mask(T3, mN3) & _mm512_cmplt_epu64_mask(T2, mN2);
	__mmask8 bmsk;
	T2 = _mm512_subsetc_epi52(T2, mN2, &bmsk);
	T3 = _mm512_sbb_epi52(T3, bmsk, mN3, &bmsk);
	T2 = _mm512_mask_addsetc_epi52(T2, msk, T2, N0, &bmsk);
	T3 = _mm512_mask_adc_epi52(T3, msk, bmsk, N1, &bmsk);
	_mm512_store_si512(p->data[0], _mm512_and_epi64(T2, lo52mask));
	_mm512_store_si512(p->data[1], _mm512_and_epi64(T3, lo52mask));
}
'''


def main():
    out = sys.argv[1]
    src = open(SRC, newline='').read()
    key = '__inline static void tecm_mulredc104_x8('
    assert src.count(key) == 1
    src = src.replace(key, NEW + '\n__inline static void tecm_mulredc104_x8_fp(')
    open(os.path.join(out, 'tinyecm_ifma104.c'), 'w', newline='').write(src)


if __name__ == '__main__':
    main()
