#pragma once
#include <immintrin.h>


#if defined( __GNUC__ ) && defined(_WIN64)

void _avx512_mask_divrem32(__m512i md, __m512i mr, __mmask16 m, 
	__m512i dividend, __m512i divisor, __m512i* q, __m512i* r);
__m512i _avx512_div32(__m512i dividend, __m512i divisor);
__m512i _mm512_rem_epu64(__m512i dividend, __m512i divisor);
__m512i _mm512_rem_epu32(__m512i dividend, __m512i divisor);

#endif


