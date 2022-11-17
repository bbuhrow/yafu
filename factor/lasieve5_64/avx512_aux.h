#pragma once
#include <immintrin.h>


#if defined( __GNUC__ ) && defined(_WIN64)

__inline void _avx512_mask_divrem32(__m512i md, __m512i mr, __mmask16 m,
	__m512i dividend, __m512i divisor, __m512i* q, __m512i* r);

__inline __m512i _avx512_div32(__m512i dividend, __m512i divisor);

//#define _avx512_div32(result, dividend, divisor) \
//{                                           \
//uint32_t memdiv[16];                         \
//uint32_t memrem[16];                         \
//                                             \
//_mm512_store_epi32(memdiv, dividend);        \
//_mm512_store_epi32(memrem, divisor);         \
//                                             \
//int i;                                       \
//for (i = 0; i < 16; i++)                     \
//{                                            \
//    uint64_t q1 = memdiv[i] / memrem[i];     \
//    memdiv[i] = q1;                          \
//}                                            \
//result = _mm512_load_epi32(memdiv);         \
//} 


__inline __m512i _mm512_rem_epu64(__m512i dividend, __m512i divisor);
__inline __m512i _mm512_rem_epu32(__m512i dividend, __m512i divisor);


#endif


