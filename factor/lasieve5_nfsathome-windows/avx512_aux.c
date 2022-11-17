#include "avx512_aux.h"
#include <immintrin.h>
#include <stdint.h>

#if defined( __GNUC__ ) && defined(_WIN64)

void _avx512_mask_divrem32(__m512i md, __m512i mr, __mmask16 m, 
    __m512i dividend, __m512i divisor, __m512i* q, __m512i* r)
{
    uint32_t memdiv[16];
    uint32_t memrem[16];
    uint32_t memmd[16];
    uint32_t memmr[16];

    _mm512_store_epi32(memdiv, dividend);
    _mm512_store_epi32(memrem, divisor);

    int i;
    for (i = 0; i < 16; i++)
    {
        if (m & (1 << i))
        {
            uint64_t q1 = memdiv[i] / memrem[i];
            uint64_t r1 = memdiv[i] % memrem[i];
            memdiv[i] = q1;
            memrem[i] = r1;
        }
    }

    *q = _mm512_mask_load_epi32(md, m, memdiv);
    *r = _mm512_mask_load_epi32(mr, m, memrem);

    return;
}

__m512i _avx512_div32(__m512i dividend, __m512i divisor)
{
    uint32_t memdiv[16];
    uint32_t memrem[16];

    _mm512_store_epi32(memdiv, dividend);
    _mm512_store_epi32(memrem, divisor);

    int i;
    for (i = 0; i < 16; i++)
    {
        uint64_t q1 = memdiv[i] / memrem[i];
        memdiv[i] = q1;
    }

    return _mm512_load_epi32(memdiv);
}

__m512i _mm512_rem_epu64(__m512i dividend, __m512i divisor)
{
    uint64_t memdiv[8];
    uint64_t memrem[8];

    _mm512_store_epi64(memdiv, dividend);
    _mm512_store_epi64(memrem, divisor);

    int i;
    for (i = 0; i < 8; i++)
    {
        uint64_t q1 = memdiv[i] % memrem[i];
        memrem[i] = q1;
    }

    return _mm512_load_epi64(memrem);
}

__m512i _mm512_rem_epu32(__m512i dividend, __m512i divisor)
{
    uint32_t memdiv[16];
    uint32_t memrem[16];

    _mm512_store_epi32(memdiv, dividend);
    _mm512_store_epi32(memrem, divisor);

    int i;
    for (i = 0; i < 16; i++)
    {
        uint32_t q1 = memdiv[i] % memrem[i];
        memrem[i] = q1;
    }

    return _mm512_load_epi32(memrem);
}

#endif
