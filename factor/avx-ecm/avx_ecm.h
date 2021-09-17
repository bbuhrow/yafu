/*
Copyright (c) 2019, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/

#include "gmp.h"
#include "ytools.h"
#include <stdint.h>
#include <immintrin.h>

#define INV_2_POW_64 5.4210108624275221700372640043497e-20

#define DO_STAGE2_INV
#define DEFINED 1
#define MAX_WINSIZE 10
#define BLOCKWORDS 4
#define strto_uint64_t strtoull

//#ifndef MAXBITS
//#define MAXBITS 1040
//#endif
//#define DIGITBITS 32

#ifndef DIGITBITS
#define DIGITBITS 52
#endif

#if DIGITBITS==52

#define base_t uint64_t
#define base_signed_t int64_t

#define VEC_HALFBITS 26
#define VEC_HALFMASK 0x3ffffff
#define VEC_MAXDIGIT 0xfffffffffffffULL
#define VEC_HIBITMASK 0x8000000000000ULL
#define VECLEN 8

#elif DIGITBITS==32

#define base_t uint32_t
#define base_signed_t int32_t

#define VEC_HALFBITS 16
#define VEC_HALFMASK 0xffff
#define VEC_MAXDIGIT 0xffffffff
#define VEC_HIBITMASK 0x80000000
#define VECLEN 16

#else
#error "DIGITBITS must be either 52 or 32"
#endif

typedef struct
{
    base_t *data;
    int size;
    uint32_t signmask;
    uint32_t WORDS_ALLOC;
} vec_bignum_t;

// a vector math library for montgomery arithmetic using AVX-512
typedef struct
{
    mpz_t nhat;
    mpz_t rhat;
    mpz_t gmp_t1;
    mpz_t gmp_t2;
    vec_bignum_t *r;
    vec_bignum_t *n;
    vec_bignum_t *vnhat;
    vec_bignum_t *vrhat;
    vec_bignum_t *rmask;
    vec_bignum_t *one;
    vec_bignum_t *mtmp1;
    vec_bignum_t *mtmp2;
    vec_bignum_t *mtmp3;
    vec_bignum_t *mtmp4;
    vec_bignum_t **g;             // storage for windowed method precomputation
    base_t *vrho;
    base_t rho;
    int nbits;
    int isMersenne;
    uint32_t MAXBITS;
    uint32_t NWORDS;
    uint32_t NBLOCKS;
    int *vnbits;
    int use_vnbits;
} vec_monty_t;

#ifdef USE_AVX512F

//#define PRINT_DEBUG 2
//#define NPRINT_DEBUG 0
//#define SPRINT_DEBUG 0
//#define DEBUGLANE 0
//#define DEBUG_MERSENNE
//#define DEBUG_VECMUL

#ifndef SKYLAKEX
#define _mm512_mullo_epi64 _mm512_mullox_epi64
// also need an answer for _mm512_cvtepu64_pd, probably among other things,
// because KNL does not support AVX512DQ
#define _mm512_cvtepu64_pd(x) _mm512_sub_pd(_mm512_castsi512_pd(_mm512_add_epi64((x), _mm512_set1_epi64(0x4330000000000000ULL))), _mm512_set1_pd(4503599627370496.))
#define _mm512_cvtpd_epu64(x) _mm512_sub_epi64(_mm512_castpd_si512(_mm512_add_pd((x), _mm512_set1_pd(4503599627370496.))), _mm512_set1_epi64(0x4330000000000000ULL))
#endif

void print_vechex(base_t* a, int v, int n, const char* pre);
void print_regvechex(__m512i a, int v, const char* pre);
void print_regvechex64(__m512i a, int v, const char* pre);
void print_regvechexrange(__m512i a, int v1, int v2, const char* pre);

// ---------------------------------------------------------------------
// emulated instructions
// ---------------------------------------------------------------------
__m512i _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout);
__m512i _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout);
__m512i _mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout);
__m512i _mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout);
__m512i _mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout);
__m512i _mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout);
__m512i _mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout);

// ---------------------------------------------------------------------
// emulated instructions
// ---------------------------------------------------------------------
//#define _mm512_addsetc_epi52(a, b, cout) {          \
//    __m512i t = _mm512_add_epi64(a, b);                                        \
//    (*(cout)) = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL)); \
//    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));            \
//    t;                                                                         \
//}
//
//#define _mm512_adc_epi52(a, c, b, cout) {                                     \
//    __m512i t = _mm512_add_epi64(a, b);                                        \
//    t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(c, 1));                    \
//    (*(cout)) = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL)); \
//    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));            \
//    t;                                                                         \
//}
//
//#define _mm512_mask_adc_epi52(a, m, c, b, cout) {                             \
//    __m512i t = _mm512_add_epi64(a, b);                                        \
//    t = _mm512_mask_add_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));         \
//    (*(cout)) = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL)); \
//    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));            \
//    t;                                                                         \
//}
//
//#define _mm512_addcarry_epi52(a, c, cout) {        \
//    __m512i t = _mm512_add_epi64(a, _mm512_maskz_set1_epi64(c, 1));            \
//    (*(cout)) = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0xfffffffffffffULL)); \
//    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));            \
//    t;                                                                         \
//}
//
//#define _mm512_subborrow_epi52(a, c, cout) {       \
//    __m512i t = _mm512_sub_epi64(a, _mm512_maskz_set1_epi64(c, 1));            \
//    (*(cout)) = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0));                  \
//    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));            \
//    t;                                                                         \
//}
//
//#define _mm512_sbb_epi52(a, c, b, cout) {  \
//__m512i t = _mm512_sub_epi64(a, b);                                            \
//(*(cout)) = _mm512_cmpgt_epu64_mask(b, a);                                         \
//__m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));               \
//(*(cout)) = _mm512_kor((*(cout)), _mm512_cmpgt_epu64_mask(t2, t));                     \
//t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));              \
//t2; }
//
//#define _mm512_mask_sbb_epi52(a, m, c, b, cout) {                             \
//__m512i t = _mm512_mask_sub_epi64(a, m, a, b);                                 \
//(*(cout)) = _mm512_mask_cmpgt_epu64_mask(m, b, a);                                 \
//__m512i t2 = _mm512_mask_sub_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));    \
//(*(cout)) = _mm512_kor((*(cout)), _mm512_mask_cmpgt_epu64_mask(m, t2, t));             \
//t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));              \
//t2; }
//
#define _mm512_adcsbb_epi52(x, c, b, y, lomask, sum, diff) { \
    __m512i ts = _mm512_add_epi64(x, y);                                        \
    __m512i td = _mm512_sub_epi64(x, y);                                  \
    __m512i carry_512i = _mm512_maskz_set1_epi64(c, 1);                       \
    __m512i borrow_512i = _mm512_maskz_set1_epi64(b, 1);                      \
    ts = _mm512_add_epi64(ts, carry_512i);                                  \
    b = _mm512_cmpgt_epu64_mask(y, x);                                \
    __m512i t2 = _mm512_sub_epi64(td, borrow_512i);                                  \
    c = _mm512_cmpgt_epu64_mask(ts, lomask);                         \
    b = _mm512_kor(b, _mm512_cmpgt_epu64_mask(t2, td));           \
    sum = _mm512_and_epi64(ts, lomask);                                 \
    diff = _mm512_and_epi64(t2, lomask);                               \
}

#define _mm512_mask_adcsbb_epi52(x, y, cm, bm, c, b, z, lomask, sum, diff) { \
    __m512i ts = _mm512_mask_add_epi64(x, cm, x, z);                                     \
    __m512i td = _mm512_mask_sub_epi64(y, bm, y, z);                         \
    __m512i carry_512i = _mm512_maskz_set1_epi64(c, 1);                      \
    __m512i borrow_512i = _mm512_maskz_set1_epi64(b, 1);                     \
    ts = _mm512_mask_add_epi64(x, cm, ts, carry_512i);                        \
    b = _mm512_mask_cmpgt_epu64_mask(bm, y, z);                               \
    __m512i t2 = _mm512_mask_sub_epi64(y, bm, td, borrow_512i);               \
    c = _mm512_mask_cmpgt_epu64_mask(cm, ts, lomask);                         \
    b = _mm512_kor(b, _mm512_mask_cmpgt_epu64_mask(bm, t2, td));           \
    sum = _mm512_and_epi64(ts, lomask);                                 \
    diff = _mm512_and_epi64(t2, lomask);                               \
}

#define castpd _mm512_castsi512_pd
#define castepu _mm512_castpd_si512

#ifdef IFMA

#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), a, b);

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
    lo = _mm512_madd52lo_epu64(lo, a, b); \
    hi = _mm512_madd52hi_epu64(hi, a, b);

#define VEC_MUL_LOHI_PD(a, b, lo, hi) \
    lo = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), a, b); \
    hi = _mm512_madd52hi_epu64(_mm512_set1_epi64(0), a, b);

#define VEC_CARRYPROP_LOHI(lo, hi) \
	a0 = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, a0);		\
	lo = _mm512_and_epi64(vlmask, lo);

#define VEC_MUL4_ACCUM(x, b0, b1, b2, b3) \
    te0 = _mm512_madd52lo_epu64(te0, x, b0); \
    te2 = _mm512_madd52lo_epu64(te2, x, b1); \
    te4 = _mm512_madd52lo_epu64(te4, x, b2); \
    te6 = _mm512_madd52lo_epu64(te6, x, b3); \
    te1 = _mm512_madd52hi_epu64(te1, x, b0); \
    te3 = _mm512_madd52hi_epu64(te3, x, b1); \
    te5 = _mm512_madd52hi_epu64(te5, x, b2); \
    te7 = _mm512_madd52hi_epu64(te7, x, b3);

#define SUB_BIAS_HI(bias1, bias2, bias3, bias4) {}
#define SUB_BIAS_LO(bias1, bias2, bias3, bias4) {}

#else

#define _mm512_mullo_epi52(c, a, b) \
    i0 = _mm512_srli_epi64(a, 32); \
    i1 = _mm512_srli_epi64(b, 32); \
    i0 = _mm512_mul_epu32(b, i0); \
    i1 = _mm512_mul_epu32(a, i1); \
    c = _mm512_mul_epu32(a, b); \
    i0 = _mm512_slli_epi64(i0, 32); \
    i1 = _mm512_slli_epi64(i1, 32); \
    c = _mm512_add_epi64(c, i0); \
    c = _mm512_add_epi64(c, i1); \
    c = _mm512_and_epi64(c, vlmask); 

//#define _mm512_mullo_epi52(c, a, b) \
//    prod1_hd = _mm512_cvtepu64_pd(a); \
//    prod2_hd = _mm512_cvtepu64_pd(b); \
//    prod3_hd = _mm512_fmadd_round_pd(prod1_hd, prod2_hd, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
//    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd); \
//    prod3_hd = _mm512_fmadd_round_pd(prod1_hd, prod2_hd, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
//    c = _mm512_sub_epi64(_mm512_castpd_si512(prod3_hd), _mm512_set1_epi64(0x4330000000000000ULL));

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi = _mm512_add_epi64(hi, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo = _mm512_add_epi64(lo, _mm512_sub_epi64(castepu(prod1_ld), vbias3));

#define VEC_MUL_LOHI_PD(a, b, lo, hi) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi = _mm512_sub_epi64(castepu(prod1_hd), vbias1); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo = _mm512_sub_epi64(castepu(prod1_ld), vbias3);

#define VEC_CARRYPROP_LOHI(lo, hi) \
	a0 = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, a0);		\
	lo = _mm512_and_epi64(vlmask, lo);

#define VEC_MUL4_ACCUM(x, b0, b1, b2, b3) \
    prod5_ld = _mm512_cvtepu64_pd(x);                                                                           \
    prod1_ld = _mm512_cvtepu64_pd(b0);                                                                          \
    prod2_ld = _mm512_cvtepu64_pd(b1);                                                                          \
    prod3_ld = _mm512_cvtepu64_pd(b2);                                                                          \
    prod4_ld = _mm512_cvtepu64_pd(b3);                                                                          \
    prod1_hd = _mm512_fmadd_round_pd(prod5_ld, prod1_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));                                                 \
    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));                                                 \
    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));                                                 \
    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));                                                 \
    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                            \
    prod1_ld = _mm512_fmadd_round_pd(prod5_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));                                                 \
    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));                                                 \
    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));                                                 \
    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));

#define VEC_SQR_MUL4_A(t0, t1, t2, t4, t5) \
prod5_ld = _mm512_cvtepu64_pd(t0);                                                                          \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                          \
prod3_ld = _mm512_cvtepu64_pd(t4);                                                                          \
prod4_ld = _mm512_cvtepu64_pd(t5);                                                                          \
prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));                                                 \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));                                                 \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));                                                 \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));                                                 \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                            \
prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));                                                 \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));                                                 \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));                                                 \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));

#define VEC_SQR_MUL4_B(t1, t2, t3, t4, t5) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                          \
prod3_ld = _mm512_cvtepu64_pd(t3);                                                                          \
prod4_ld = _mm512_cvtepu64_pd(t4);                                                                          \
prod5_ld = _mm512_cvtepu64_pd(t5);                                                                          \
prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod5_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));                                                 \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));                                                 \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));                                                 \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));                                                 \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                            \
prod5_ld = _mm512_fmadd_round_pd(prod2_ld, prod5_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));                                                 \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));                                                 \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));                                                 \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod5_ld));

#define VEC_SQR_MUL4_C(t1, t2, t3, t4) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                          \
prod3_ld = _mm512_cvtepu64_pd(t3);                                                                          \
prod4_ld = _mm512_cvtepu64_pd(t4);                                                                          \
prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod1_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));                                                 \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));                                                 \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));                                                 \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));                                                 \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                            \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
prod5_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));                                                 \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));                                                 \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));                                                 \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod5_ld));

#define VEC_SQR_MUL4_D(t1, t2, t3, t4, t5, o1, o2, o3, o4) \
prod5_ld = _mm512_cvtepu64_pd(t1);                                                                         \
prod1_ld = _mm512_cvtepu64_pd(t2);                                                                         \
prod2_ld = _mm512_cvtepu64_pd(t3);                                                                         \
prod3_ld = _mm512_cvtepu64_pd(t4);                                                                         \
prod4_ld = _mm512_cvtepu64_pd(t5);                                                                         \
prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias,                                                \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                 \
prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,                                                \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                 \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,                                                \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                 \
prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias,                                                \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                 \
o2 = _mm512_add_epi64(o2, _mm512_castpd_si512(prod2_hd));                                                \
o4 = _mm512_add_epi64(o4, _mm512_castpd_si512(prod4_hd));                                                \
o2 = _mm512_add_epi64(o2, _mm512_castpd_si512(prod1_hd));                                                \
o4 = _mm512_add_epi64(o4, _mm512_castpd_si512(prod3_hd));                                                \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                           \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                           \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                           \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                           \
prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));  \
prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));  \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));  \
prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));  \
o1 = _mm512_add_epi64(o1, _mm512_castpd_si512(prod2_ld));                                                \
o3 = _mm512_add_epi64(o3, _mm512_castpd_si512(prod4_ld));                                                \
o1 = _mm512_add_epi64(o1, _mm512_castpd_si512(prod1_ld));                                                \
o3 = _mm512_add_epi64(o3, _mm512_castpd_si512(prod3_ld));

#define VEC_MUL_MUL4_A(a0, a1, a2, b2, b3) \
prod5_ld = _mm512_cvtepu64_pd(a0);                                                                        \
prod1_ld = _mm512_cvtepu64_pd(a1);                                                                        \
prod2_ld = _mm512_cvtepu64_pd(a2);                                                                        \
prod3_ld = _mm512_cvtepu64_pd(b2);                                                                        \
prod4_ld = _mm512_cvtepu64_pd(b3);                                                                        \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));                                               \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));                                               \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));                                               \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));                                               \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                          \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                          \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                          \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                          \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));                                               \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));                                               \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));                                               \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));

#define VEC_MUL_MUL4_B(a0, a1, b0, b1, b2) \
prod5_ld = _mm512_cvtepu64_pd(a0);                                                                        \
prod1_ld = _mm512_cvtepu64_pd(a1);                                                                        \
prod2_ld = _mm512_cvtepu64_pd(b0);                                                                        \
prod3_ld = _mm512_cvtepu64_pd(b1);                                                                        \
prod4_ld = _mm512_cvtepu64_pd(b2);                                                                        \
prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));                                               \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));                                               \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));                                               \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));                                               \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                          \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                          \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                          \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                          \
prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));                                               \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));                                               \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));                                               \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));

#define VEC_SQR_MUL2(t1, t2) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                          \
prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));                                                 \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));                                                 \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));                                                 \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));

#define VEC_SQR_MUL2_C(t1, t2) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                          \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));                                                 \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));                                                 \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));                                                 \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));

#define VEC_SQR_MUL2_D(t1, t2, t3) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                          \
prod3_ld = _mm512_cvtepu64_pd(t3);                                                                          \
prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));                                                 \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));                                                 \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));                                                 \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));

#define VEC_SQR_MUL2_B(t1, t2, t3) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                          \
prod3_ld = _mm512_cvtepu64_pd(t3);                                                                          \
prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,                                                 \
(_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));                                                                  \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));                                                 \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));                                                 \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                            \
prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));                                                 \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));

#define VEC_MUL_MUL2_A(a2, a3, b2, b3) \
prod1_ld = _mm512_cvtepu64_pd(a2);                                                                          \
prod2_ld = _mm512_cvtepu64_pd(a3);                                                                          \
prod3_ld = _mm512_cvtepu64_pd(b2);                                                                          \
prod4_ld = _mm512_cvtepu64_pd(b3);                                                                          \
prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));      \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));                                                 \
te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));                                                 \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                            \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                            \
prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));   \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));                                                 \
te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld))

#define VEC_SQR_MUL3_A(t1, t2, t3, t4) \
    prod1_ld = _mm512_cvtepu64_pd(t1);                                                                        \
    prod2_ld = _mm512_cvtepu64_pd(t2);                                                                        \
    prod3_ld = _mm512_cvtepu64_pd(t3);                                                                        \
    prod4_ld = _mm512_cvtepu64_pd(t4);                                                                        \
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
    prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
    prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));    \
    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));                                               \
    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));                                               \
    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));                                               \
    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                          \
    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                          \
    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                          \
    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));                                               \
    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));                                               \
    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));

#define VEC_SQR_MUL3_B(t1, t2, t3, t4) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                            \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                            \
prod3_ld = _mm512_cvtepu64_pd(t3);                                                                            \
prod4_ld = _mm512_cvtepu64_pd(t4);                                                                            \
prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));        \
prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));        \
prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));        \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));                                                   \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));                                                   \
te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));                                                   \
prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);                                              \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                              \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                              \
prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));     \
prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));     \
prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));     \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));                                                   \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));                                                   \
te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));

#define VEC_SQR_MUL3_C(t1, t2, t3, t4) \
prod1_ld = _mm512_cvtepu64_pd(t1);                                                                            \
prod2_ld = _mm512_cvtepu64_pd(t2);                                                                            \
prod3_ld = _mm512_cvtepu64_pd(t3);                                                                            \
prod4_ld = _mm512_cvtepu64_pd(t4);                                                                            \
prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));        \
prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));        \
prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));        \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));                                                   \
te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));                                                   \
te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));                                                   \
prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);                                              \
prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);                                              \
prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);                                              \
prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));     \
prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));     \
prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));     \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));                                                   \
te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));                                                   \
te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));


#define CARRYPROP_ACCUM_0(out) \
acc_e0 = _mm512_add_epi64(acc_e0, te0);                             \
acc_e1 = _mm512_add_epi64(acc_e1, te1);                             \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_srli_epi64(acc_e1, 52);                             \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
_mm512_mullo_epi52(a0, nhatvec_e, acc_e0);                          \
out = a0;                                                           \
b0 = _mm512_load_epi64(n->data + 0 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);                      \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));   \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
acc_e0 = acc_e1;                                                    \
acc_e1 = acc_e2;                                                    \
acc_e2 = zero;

#define CARRYPROP_ACCUM_1(out, in1) \
acc_e0 = _mm512_add_epi64(acc_e0, te2);                             \
acc_e1 = _mm512_add_epi64(acc_e1, te3);                             \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_srli_epi64(acc_e1, 52);                             \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
b1 = _mm512_load_epi64(n->data + 1 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(in1, b1, acc_e0, acc_e1);                     \
_mm512_mullo_epi52(a0, nhatvec_e, acc_e0);                          \
out = a0;                                                           \
b0 = _mm512_load_epi64(n->data + 0 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);                      \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));   \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
acc_e0 = acc_e1;                                                    \
acc_e1 = acc_e2;                                                    \
acc_e2 = zero;

#define CARRYPROP_ACCUM_2(out, in1, in2) \
acc_e0 = _mm512_add_epi64(acc_e0, te4);                             \
acc_e1 = _mm512_add_epi64(acc_e1, te5);                             \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_srli_epi64(acc_e1, 52);                             \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
b1 = _mm512_load_epi64(n->data + 2 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(in1, b1, acc_e0, acc_e1);                     \
b1 = _mm512_load_epi64(n->data + 1 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(in2, b1, acc_e0, acc_e1);                     \
_mm512_mullo_epi52(a0, nhatvec_e, acc_e0);                          \
out = a0;                                                           \
b0 = _mm512_load_epi64(n->data + 0 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);                      \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));   \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
acc_e0 = acc_e1;                                                    \
acc_e1 = acc_e2;                                                    \
acc_e2 = zero;

#define CARRYPROP_ACCUM_3(out, in1, in2, in3) \
acc_e0 = _mm512_add_epi64(acc_e0, te6);                             \
acc_e1 = _mm512_add_epi64(acc_e1, te7);                             \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_srli_epi64(acc_e1, 52);                             \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
b1 = _mm512_load_epi64(n->data + 3 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(in1, b1, acc_e0, acc_e1);                     \
b1 = _mm512_load_epi64(n->data + 2 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(in2, b1, acc_e0, acc_e1);                     \
b1 = _mm512_load_epi64(n->data + 1 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(in3, b1, acc_e0, acc_e1);                     \
_mm512_mullo_epi52(a0, nhatvec_e, acc_e0);                          \
out = a0;                                                           \
b0 = _mm512_load_epi64(n->data + 0 * VECLEN);                       \
VEC_MUL_ACCUM_LOHI_PD(a0, b0, acc_e0, acc_e1);                      \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                                 \
acc_e2 = _mm512_add_epi64(acc_e2, _mm512_srli_epi64(acc_e1, 52));   \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);                          \
acc_e0 = acc_e1;                                                    \
acc_e1 = acc_e2;                                                    \
acc_e2 = zero;

#define ACCUM_SHIFT(out, inlo, inhi) \
acc_e0 = _mm512_add_epi64(acc_e0, inlo);           \
acc_e1 = _mm512_add_epi64(acc_e1, inhi);           \
VEC_CARRYPROP_LOHI(acc_e0, acc_e1);                \
acc_e2 = _mm512_srli_epi64(acc_e1, 52);            \
acc_e1 = _mm512_and_epi64(acc_e1, vlmask);         \
out = acc_e0;                                       \
acc_e0 = acc_e1;                                   \
acc_e1 = acc_e2;                                   \
acc_e2 = zero;

// avoid the slow _mm512_mullo_epi64 by using _mm512_mul_epu32 and a 32-bit shift.
// works because the bias's low 32-bits are all zero.
#define SUB_BIAS_HI(bias1, bias2, bias3, bias4) \
    b0 = _mm512_set1_epi64((bias1) * 0x467);         \
    b1 = _mm512_set1_epi64((bias2) * 0x467);         \
    b2 = _mm512_set1_epi64((bias3) * 0x467);         \
    b3 = _mm512_set1_epi64((bias4) * 0x467);         \
    b0 = _mm512_slli_epi64(b0, 52);            \
    b1 = _mm512_slli_epi64(b1, 52);            \
    b2 = _mm512_slli_epi64(b2, 52);            \
    b3 = _mm512_slli_epi64(b3, 52);            \
    te1 = _mm512_sub_epi64(te1, b0); \
    te3 = _mm512_sub_epi64(te3, b1); \
    te5 = _mm512_sub_epi64(te5, b2); \
    te7 = _mm512_sub_epi64(te7, b3);


#define SUB_BIAS_LO(bias1, bias2, bias3, bias4) \
    b0 = _mm512_set1_epi64((bias1) * 0x433);         \
    b1 = _mm512_set1_epi64((bias2) * 0x433);         \
    b2 = _mm512_set1_epi64((bias3) * 0x433);         \
    b3 = _mm512_set1_epi64((bias4) * 0x433);         \
    b0 = _mm512_slli_epi64(b0, 52);            \
    b1 = _mm512_slli_epi64(b1, 52);            \
    b2 = _mm512_slli_epi64(b2, 52);            \
    b3 = _mm512_slli_epi64(b3, 52);            \
    te0 = _mm512_sub_epi64(te0, b0); \
    te2 = _mm512_sub_epi64(te2, b1); \
    te4 = _mm512_sub_epi64(te4, b2); \
    te6 = _mm512_sub_epi64(te6, b3);

#endif

#endif

void print_hexbignum(vec_bignum_t *a, const char *pre);
void print_vechexbignum(vec_bignum_t* a, const char* pre);
void print_hex(vec_bignum_t *a, const char *pre);
void print_vechex(base_t *a, int v, int n, const char *pre);
vec_monty_t * vec_monty_alloc(uint32_t words);
void vec_monty_free(vec_monty_t *mdata);
void copy_vec_lane(vec_bignum_t *src, vec_bignum_t *dest, int num, int size);
void vecCopy(vec_bignum_t * src, vec_bignum_t * dest);
void vecCopyn(vec_bignum_t * src, vec_bignum_t * dest, int size);
void vecClear(vec_bignum_t *n);
vec_bignum_t * vecInit(uint32_t words);
void vecFree(vec_bignum_t *);

// 52-BIT functions
void vecmulmod52_1(vec_bignum_t *a, base_t *b, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecredc52_base(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecmulmod52_fixed1040_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecmulmod52_fixed832_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecmulmod52_fixed624_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecmulmod52_fixed416_bfips(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsqrmod52_fixed1040_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsqrmod52_fixed832_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsqrmod52_fixed624_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsqrmod52_fixed416_bfips(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecmulmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecsqrmod52(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsubmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata);
void vecsubmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void vecaddmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata);
void vecaddmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void vec_simul_addsub52_fixed1040(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *sum, vec_bignum_t *diff, 
    vec_monty_t* mdata);
void vec_simul_addsub52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
    vec_monty_t* mdata);
void vec_simul_addsub52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
    vec_monty_t* mdata);
void vec_bignum52_mask_sub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, uint32_t wmask);
void vec_bignum52_mask_rshift_1(vec_bignum_t * u, uint32_t wmask);
uint32_t vec_bignum52_mask_lshift_1(vec_bignum_t * u, uint32_t wmask);
uint32_t vec_bignum52_mask_lshift_n(vec_bignum_t * u, int n, uint32_t wmask);
uint32_t vec_eq52(base_t * u, base_t * v, int sz);
uint32_t vec_gte52(vec_bignum_t * u, vec_bignum_t * v);
void vec_bignum52_add_1(vec_bignum_t *a, base_t *b, vec_bignum_t *c);
void vec_bignum52_add(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c);
void vec_bignum52_sub(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c);
void vec_bignum52_mask_rshift_n(vec_bignum_t* u, vec_bignum_t* v, int n, uint32_t wmask);
void vec_bignum52_mask_rshift_vn(vec_bignum_t* u, vec_bignum_t* v, int* n, uint32_t wmask);
void vecmul52_1(vec_bignum_t* a, __m512i b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);


void vecmod_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* s, vec_monty_t* mdata);
void vecmul52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void veckmul(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void vecksqr_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void veckmul_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);


// 32-BIT functions
void vecmulmod(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecsqrmod(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecmulmod_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsqrmod_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsubmod(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata);
void vecaddmod(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata);
void vecaddmod_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void vecsubmod_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void vec_simul_addsub_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
    vec_monty_t* mdata);
void vec_simul_addsub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *sum, vec_bignum_t *diff, vec_monty_t* mdata);
void vec_bignum_mask_sub(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, uint32_t wmask);
void vec_bignum_mask_rshift_1(vec_bignum_t * u, uint32_t wmask);
uint32_t vec_bignum_mask_lshift_1(vec_bignum_t * u, uint32_t wmask);
uint32_t vec_eq(base_t * u, base_t * v, int sz);
uint32_t vec_gte(vec_bignum_t * u, vec_bignum_t * v);

void extract_bignum_from_vec_to_mpz(mpz_t dest, vec_bignum_t *vec_src, int num, int sz);
void broadcast_mpz_to_vec(vec_bignum_t *vec_dest, mpz_t src);
void insert_mpz_to_vec(vec_bignum_t *vec_dest, mpz_t src, int lane);

// declare function pointers to the type of reduction needed
extern void(*vecmulmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
extern void(*vecsqrmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
extern void(*vecaddmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t*);
extern void(*vecsubmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t*);
extern void(*vecaddsubmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t*);

// ecm stuff
typedef struct 
{
	vec_bignum_t *X;
	vec_bignum_t *Z;
} ecm_pt;

typedef struct 
{
	vec_bignum_t *sum1;
	vec_bignum_t *diff1;
	vec_bignum_t *sum2;
	vec_bignum_t *diff2;
	vec_bignum_t *tt1;
	vec_bignum_t *tt2;
	vec_bignum_t *tt3;
	vec_bignum_t *tt4;
	vec_bignum_t *tt5;
	vec_bignum_t *s;
	vec_bignum_t *n;
	ecm_pt pt1;
	ecm_pt pt2;
	ecm_pt pt3;
	ecm_pt pt4;
	ecm_pt pt5;
	uint64_t sigma;

	uint32_t *map;
	ecm_pt *Pa;
	ecm_pt *Pd;
	ecm_pt *Pad;
	vec_bignum_t **Paprod;	
    vec_bignum_t** Pa_inv;
	vec_bignum_t **Pbprod;
	ecm_pt *Pb;
    ecm_pt* Pdnorm;
	vec_bignum_t *stg2acc;
	uint32_t stg1Add;
	uint32_t stg1Doub;
    uint32_t paired;
	uint32_t ptadds;
    uint32_t ptdups;
    uint32_t numinv;
	uint64_t numprimes;
    uint64_t A;
    uint32_t last_pid;
	uint32_t amin;

	uint32_t U;
	uint32_t L;
	uint32_t D;
	uint32_t R;

    uint64_t* primes;
    uint64_t num_p;
    uint64_t min_p;
    uint64_t max_p;

    uint64_t STAGE1_MAX;
    uint64_t STAGE2_MAX;
    uint32_t NWORDS;
    uint32_t NBLOCKS;
    uint32_t MAXBITS;

} ecm_work;

// a wrapper for found factors and information about
// how and where they were found.
typedef struct
{
    mpz_t factor;
    uint64_t sigma;
    int thread_id;
    int vec_id;
    int curve_id;
    int stg_id;
    uint32_t B1;
    uint64_t B2;
} avx_ecm_factor_t;

typedef struct
{
    avx_ecm_factor_t *factors;
    int numfactors;
    uint64_t *sigma;
    vec_monty_t *mdata;
    ecm_work *work;
    uint32_t curves;
    uint32_t b1;
    uint32_t b2;
    ecm_pt *P;
    uint64_t lcg_state;
    uint32_t tid;
    uint32_t total_threads;
    uint32_t phase_done;
    uint32_t ecm_phase;     // 0 == build curve, 1 == stage 1, 2 == stage 2
    uint32_t* pairmap_v;
    uint32_t* pairmap_u;
    uint32_t pairmap_steps;
    uint32_t* Qmap;
    uint32_t* Qrmap;
    Queue_t** Q;
    int verbose;
    int save_b1;

    uint64_t* primes;
    uint64_t nump;
    uint64_t minp;
    uint64_t maxp;

    uint32_t MAXBITS;
    uint32_t NWORDS;
    uint32_t NBLOCKS;
    uint64_t STAGE1_MAX;
    uint64_t STAGE2_MAX;
    uint32_t PRIME_RANGE;
    int DO_STAGE2;
} thread_data_t;

void vececm(thread_data_t *tdata);
void vec_ecm_pt_init(ecm_pt *pt, uint32_t words);
void vec_ecm_pt_free(ecm_pt *pt);
void vec_ecm_work_init(ecm_work *work);
void vec_ecm_work_free(ecm_work *work);




