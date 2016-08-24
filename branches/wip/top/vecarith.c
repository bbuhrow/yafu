#include <stdint.h>
#include <immintrin.h>


void vec_u32_qr_idiv(uint32_t *q, uint32_t *r, uint32_t *a, uint32_t *b)
{
    // divides each element of the input array 'a' by each element of the 
    // input array 'b' to produce a quotient 'q' and remainder 'r' for
    // each element.

    /*"vmovapd	(%4), %%ymm3 \n\t"       vector of eight 2.0 values */
    /*"vmovapd	(%5), %%ymm6 \n\t"       vector of eight integer 1 values */
    /*"vmovapd	(%6), %%ymm10 \n\t"      vector of eight 1/65536.0 values */
    /*"vmovapd	(%7), %%ymm7 \n\t"       vector of eight magic uint32 -> double conversion values */

    
    
    __m256i vb = _mm256_load_si256((__m256i *)b);
    __m256 vbf;
    __m256i magic;
    __m256i vb64a, vb64b;
    __m128 vbf_128_a, vbf_128_b;

    /* "vcvtdq2ps	    (%1, %%rcx, 4), %%ymm1 \n\t"	    convert denominator uint32 to float */
    vbf = _mm256_cvtepi32_ps(vb);
    /* "vpmovzxdq      16(%1, %%rcx, 4), %%ymm11 \n\t"     zero extend 32-bit to 64-bit */
    vb64a = _mm256_cvtepu32_epi64(_mm_load_si128((__m128i *)b));
    /* "vrcpps	        %%ymm1, %%ymm1 \n\t"	            approx 1/d */
    vbf = _mm256_rcp_ps(vbf);
    /* "vpmovzxdq      (%1, %%rcx, 4), %%ymm0 \n\t"        zero extend 32-bit to 64-bit */
    vb64b = _mm256_cvtepu32_epi64(_mm_load_si128((__m128i *)&b[4]));
    /* "vextracti128   $1,%%ymm1,%%xmm14 \n\t"             grab high half of results */
    vbf_128_a = _mm256_extractf128_ps(vbf, 1);
    /* "vpaddq         %%ymm0, %%ymm7, %%ymm0 \n\t"        add magic constant as integer */
    vb64b = _mm256_add_epi64(vb64b, magic);
    /* "vcvtps2pd      %%xmm1, %%ymm2 \n\t"	            convert back to double precision */
    /* "vpaddq         %%ymm11, %%ymm7, %%ymm11 \n\t"      add magic constant as integer */
    /* "vsubpd         %%ymm7, %%ymm0, %%ymm0 \n\t"        sub magic constant as double */
    /* "vcvtps2pd      %%xmm14, %%ymm14 \n\t"	            convert back to double precision */
    /* "vsubpd         %%ymm7, %%ymm11, %%ymm11 \n\t"      sub magic constant as double */


    /* "vmulpd	        %%ymm2, %%ymm0, %%ymm4 \n\t"	        d * 1/d */
    /* "vmulpd	        %%ymm14, %%ymm11, %%ymm15 \n\t"	        d * 1/d */
    /* "vcvtdq2pd	    (%2, %%rcx, 4), %%ymm1 \n\t"	        convert numerator uint32 to double */
    /* "vmulpd	        %%ymm4, %%ymm2, %%ymm4 \n\t"	        d * 1/d * 1/d */
    /* "vmulpd	        %%ymm15, %%ymm14, %%ymm15 \n\t"	        d * 1/d * 1/d */
    /* "vfmsub132pd	%%ymm3, %%ymm4, %%ymm2 \n\t"	        1/d = 2 * 1/d - d * 1/d^2 */
    /* "vcvtdq2pd	    16(%2, %%rcx, 4), %%ymm12 \n\t"	        convert numerator uint32 to double */
    /* "vfmsub132pd	%%ymm3, %%ymm15, %%ymm14 \n\t"	        1/d = 2 * 1/d - d * 1/d^2 */
    /* "vmulpd	        %%ymm2, %%ymm0, %%ymm4 \n\t"	        d * 1/d */
    /* "vmulpd	        %%ymm14, %%ymm11, %%ymm15 \n\t"	        d * 1/d */
    /* "vmovapd        %%ymm1, %%ymm8 \n\t"
    /* "vmulpd	        %%ymm4, %%ymm2, %%ymm4 \n\t"	        d * 1/d * 1/d */
    /* "vmulpd	        %%ymm15, %%ymm14, %%ymm15 \n\t"	        d * 1/d * 1/d */
    /* "vfmsub132pd	%%ymm3, %%ymm4, %%ymm2 \n\t"	        1/d = 2 * 1/d - d * 1/d^2 */
    /* "vmovapd        %%ymm12, %%ymm13 \n\t"
    /* "vfmsub132pd	%%ymm3, %%ymm15, %%ymm14 \n\t"	        1/d = 2 * 1/d - d * 1/d^2 */
    /* "vfmadd132pd	%%ymm2, %%ymm10, %%ymm1 \n\t"	        n * 1/d + 1/65536 */
    /* "vfmadd132pd	%%ymm14, %%ymm10, %%ymm12 \n\t"	        n * 1/d + 1/65536 */
    /* "vroundpd	    $3, %%ymm1, %%ymm9 \n\t"	            truncate */
    /* "vroundpd	    $3, %%ymm12, %%ymm14 \n\t"	            truncate */
    /* "vmulpd	        %%ymm9, %%ymm0, %%ymm4 \n\t"	        ans * denominators */
    /* "vmulpd	        %%ymm14, %%ymm11, %%ymm15 \n\t"	        ans * denominators */
    /* "vcmppd         $0xe, %%ymm8, %%ymm4, %%ymm5 \n\t"      if (ymm4[j] > numerator[j]) set ymm5[j] = 0xffffffffffffffff, for j=0:3 */
    /* "vcmppd         $0xe, %%ymm13, %%ymm15, %%ymm12 \n\t"   if (ymm4[j] > numerator[j]) set ymm5[j] = 0xffffffffffffffff, for j=0:3 */
    /* "vpand          %%ymm6, %%ymm5, %%ymm5 \n\t"
    /* "vpand          %%ymm6, %%ymm12, %%ymm12 \n\t"
    /* "vsubpd         %%ymm5, %%ymm9, %%ymm1 \n\t"
    /* "vsubpd         %%ymm12, %%ymm14, %%ymm2 \n\t"
    /* "vcvttpd2dq	    %%ymm1, %%xmm4 \n\t"	                truncate to uint32 */
    /* "vcvttpd2dq	    %%ymm2, %%xmm5 \n\t"	                truncate to uint32 */
    /* "vmovntdq       %%xmm4, (%3, %%rcx, 4) \n\t"            move results out (non-temporal) */
    /* "vmovntdq       %%xmm5, 16(%3, %%rcx, 4) \n\t"          move results out (non-temporal) */


    return;
}