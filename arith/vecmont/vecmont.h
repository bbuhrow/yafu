/*
Copyright (c) 2026, Ben Buhrow
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
*/

/*
 * Vector Montgomery arithmetic: shared types and the constant-S template
 * mul/sqr API (vecarith52_tmpl.c).  8 numbers per AVX-512 vector (one per
 * lane, radix 2^52, R = 2^(52*NWORDS)), all with the same modulus.
 *
 * Owner of the digit-size macros, vec_bignum_t and vec_monty_t.  AVX-ECM
 * (factor/avx-ecm/avx_ecm.h) includes this header and builds on it.
 */

#ifndef VECMONT_H
#define VECMONT_H

#include "gmp.h"
#include <stdint.h>

#define DIGITBITS 52

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

// Signed-int aliases for VECLEN and BLOCKWORDS.
// Use these in pointer arithmetic expressions that mix int loop variables
// with VECLEN/BLOCKWORDS to avoid -Wsign-conversion warnings.
// e.g.:  a->data + (i * BLOCKWORDS_I + j) * VECLEN_I
// instead of: a->data + (i * BLOCKWORDS + j) * VECLEN
#define VECLEN_I     ((int)VECLEN)

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

typedef void (*vecmont_mul_fn)(vec_bignum_t *a, vec_bignum_t *b,
    vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
typedef void (*vecmont_sqr_fn)(vec_bignum_t *a, vec_bignum_t *c,
    vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);

#if defined(USE_AVX512F) && defined(IFMA)
// Constant-S template mul/sqr cores (vecarith52_tmpl.c), any NWORDS >= 2:
// c = a*b/R mod N in almost-Montgomery form (c < R).  Inputs < R with
// normalized 52-bit words; c may alias a or b; s (>= S words) must not
// alias a, b or c.  N needs >= 3 spare bits below R if the outputs feed
// AVX-ECM's add/sub (see AVXECM_SPARE_BITS).
// Sets *mul / *sqr for NWORDS = S; returns 0 if S is unsupported (S < 2).
int vecmulmod52_tmpl_get(int S, vecmont_mul_fn *mul, vecmont_sqr_fn *sqr);
#endif

#endif /* VECMONT_H */
