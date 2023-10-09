/*
Copyright (c) 2014, Ben Buhrow
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


Copyright (c) 2018 by The Mayo Clinic, though its Special Purpose
 Processor Development Group (SPPDG). All Rights Reserved Worldwide.
 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain
 a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied,
           including conditions of title, non-infringement, merchantability,
          or fitness for a particular purpose
 See the License for the specific language governing permissions and
 limitations under the License.
This file is a snapshot of a work in progress, originated by Mayo
 Clinic SPPDG.
*/

#include "avx_ecm.h"

#ifdef _MSC_VER
//#define USE_AVX512F
#endif

#ifdef USE_AVX512F
#include <immintrin.h>


// ---------------------------------------------------------------------
// emulated instructions
// in speed-critical functions these are not used... noticably faster to macro them
// ---------------------------------------------------------------------
__m512i _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i _mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_mask_add_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i _mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = c & _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i _mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

__m512i _mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}

__m512i _mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_mask_sub_epi64(a, m, a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
    __m512i t2 = _mm512_mask_sub_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}

void print_vechex(base_t* a, int v, int n, const char* pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    int i;
    if (pre != NULL)
        printf("%s: ", pre);
    for (i = n - 1; i >= 0; i--)
    {
        printf("%013"PRIx64"", a[v + i * VECLEN]);
    }
    printf("\n");
    return;
}

void print_regvechex(__m512i a, int v, const char* pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    base_t aa[VECLEN];

    _mm512_store_epi64(aa, a);

    if (pre != NULL)
        printf("%s: ", pre);

    printf("%08"PRIx64"", aa[v]);
    if (pre != NULL)
        printf("\n");
    return;
}

void print_regvechexrange(__m512i a, int v1, int v2, const char* pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    base_t aa[VECLEN];
    int i;

    _mm512_store_epi64(aa, a);

    if (pre != NULL)
        printf("%s: ", pre);

    for (i = v2; i >= v1; i--)
        printf("%08"PRIx64"", aa[i]);

    if (pre != NULL)
        printf("\n");
    return;
}

void print_hexbignum(vec_bignum_t* a, const char* pre)
{
    // print n hexdigits of the bignum a
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    if (pre != NULL)
        printf("%s", pre);

    for (i = NWORDS - 1; i >= 0; i--)
    {
        printf("%013"PRIx64"", a->data[i]);
    }
    printf("\n");
    return;
}

void print_vechexbignum(vec_bignum_t* a, const char* pre)
{
    // print n hexdigits of the bignum a
    int i, j;
    uint32_t NWORDS = a->WORDS_ALLOC;
    if (pre != NULL)
        printf("%s\n", pre);

    for (j = 0; j < VECLEN; j++)
    {
#ifdef DEBUGLANE
        if (j != DEBUGLANE)
            continue;
#endif
        //for (i = 2 * NWORDS - 1; i >= 0; i--)
        for (i = NWORDS - 1; i >= 0; i--)
        {
            printf("%016llx", a->data[i * VECLEN + j]);
        }
        printf("\n");
    }
    return;
}

void print_regvechex64(__m512i a, int v, const char* pre)
{
    // print n hexdigits of the v'th position of vec bignum a
    uint64_t aa[VECLEN];

    _mm512_store_epi64(aa, a);

    if (pre != NULL)
        printf("%s: ", pre);

    printf("%016"PRIx64"", aa[v]);
    if (pre != NULL)
        printf("\n");
    return;
}

uint32_t vec_gte52(vec_bignum_t* u, vec_bignum_t* v)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __mmask8 mdecided = 0;
    __mmask8 mgte = 0;

    for (i = NWORDS - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi64(u->data + i * VECLEN);
        __m512i b = _mm512_load_epi64(v->data + i * VECLEN);

        mgte |= _mm512_mask_cmp_epu64_mask(~mdecided, a, b, _MM_CMPINT_GT);
        mdecided = mdecided | _mm512_mask_cmp_epu64_mask(~mdecided, a, b, _MM_CMPINT_LT);

        if (mdecided == 0xff)
            break;
    }

    //equal if still undecided
    mgte |= ~mdecided;

    return (uint32_t)mgte;
}

uint32_t vec_eq52(base_t* u, base_t* v, int sz)
{
    // decide if each of the bignums in vec 'u' is >=
    // the corresponding bignum in vec 'v'.
    // return a mask of results.
    int i;
    __mmask8 meq = 0xff;

    for (i = sz - 1; i >= 0; --i)
    {
        __m512i a = _mm512_load_epi64(u + i * VECLEN);
        __m512i b = _mm512_load_epi64(v + i * VECLEN);

        meq = _mm512_mask_cmp_epu64_mask(meq, a, b, _MM_CMPINT_EQ);

        if (meq == 0)
            break;
    }

    return (uint32_t)meq;
}

uint32_t vec_bignum52_mask_lshift_1(vec_bignum_t* u, uint32_t wmask)
{
    // return the left shift of bignum u by 1
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i highmask = _mm512_set1_epi64(VEC_MAXDIGIT);
    __m512i word;

    for (i = 0; i < NWORDS; i++)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_srli_epi64(word, (DIGITBITS - 1));
        _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
            _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, 1), carry)));
        carry = nextcarry;
    }

    _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
        _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, 1), carry)));

    // return an overflow mask
    return wmask & _mm512_cmp_epi64_mask(carry, _mm512_set1_epi64(0), _MM_CMPINT_GT);
}

uint32_t vec_bignum52_mask_lshift_n(vec_bignum_t* u, int n, uint32_t wmask)
{
    // return the left shift of bignum u by n bits
    // n is assumed less than DIGITBITS
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i highmask = _mm512_set1_epi64(VEC_MAXDIGIT);
    __m512i word;

    if (n > DIGITBITS)
    {
        printf("error, vec_bignum52_mask_lshift_n expects n < %d\n", DIGITBITS);
        exit(0);
    }

    for (i = 0; i < NWORDS; i++)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_srli_epi64(word, (DIGITBITS - n));
        _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
            _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, n), carry)));
        carry = nextcarry;
    }

    _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
        _mm512_and_epi64(highmask, _mm512_or_epi64(_mm512_slli_epi64(word, n), carry)));

    // return an overflow mask
    return wmask & _mm512_cmp_epi64_mask(carry, _mm512_set1_epi64(0), _MM_CMPINT_GT);
}

void vec_bignum52_mask_rshift_n(vec_bignum_t* u, vec_bignum_t* v, int n, uint32_t wmask)
{
    // return the right shift of bignum u by n bits
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i lowmask;
    __m512i word;
    int wshift = n / 52;
    int bshift = n % 52;

    lowmask = _mm512_set1_epi64((1ULL << (uint64_t)bshift) - 1ULL);

    for (i = 2 * NWORDS - 1; (i - wshift) >= 0; i--)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_slli_epi64(_mm512_and_epi64(word, lowmask), (DIGITBITS - bshift));
        _mm512_mask_store_epi64(v->data + (i - wshift) * VECLEN, (__mmask8)wmask,
            _mm512_or_epi64(_mm512_srli_epi64(word, bshift), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum52_mask_rshift_vn(vec_bignum_t* u, vec_bignum_t* v, int* n, uint32_t wmask)
{
    // return the right shift of bignum u by n bits
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i lowmask, vbshift, vbnshift;
    __m512i word;

    int bshift_a[VECLEN];
    int wshift = n[0] / 52;
    for (i = 0; i < VECLEN; i++)
    {
        bshift_a[i] = n[i] % 52;
    }

    lowmask = _mm512_set_epi64(
        (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);

    vbshift = _mm512_set_epi64(
        bshift_a[7],
        bshift_a[6],
        bshift_a[5],
        bshift_a[4],
        bshift_a[3],
        bshift_a[2],
        bshift_a[1],
        bshift_a[0]);

    vbnshift = _mm512_sub_epi64(_mm512_set1_epi64(DIGITBITS), vbshift);

    for (i = 2 * NWORDS - 1; (i - wshift) >= 0; i--)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_sllv_epi64(_mm512_and_epi64(word, lowmask), vbnshift);
        _mm512_mask_store_epi64(v->data + (i - wshift) * VECLEN, (__mmask8)wmask,
            _mm512_or_epi64(_mm512_srlv_epi64(word, vbshift), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum52_mask_rshift_1(vec_bignum_t* u, uint32_t wmask)
{
    // return the right shift of bignum u by 1
    int i;
    uint32_t NWORDS = u->WORDS_ALLOC;
    __m512i nextcarry;
    __m512i carry = _mm512_set1_epi64(0);
    __m512i lowmask = _mm512_set1_epi64(1);
    __m512i word;

    //carry = 0;
    //for (i = sb - 1; i >= 0; --i)
    //{
    //    nextcarry = (b->data[i] & mask) << y;
    //    a->data[i] = b->data[i] >> x | carry;
    //    carry = nextcarry;
    //}

    for (i = NWORDS - 1; i >= 0; i--)
    {
        word = _mm512_load_epi64(u->data + i * VECLEN);
        nextcarry = _mm512_slli_epi64(_mm512_and_epi64(word, lowmask), (DIGITBITS - 1));
        _mm512_mask_store_epi64(u->data + i * VECLEN, (__mmask8)wmask,
            _mm512_or_epi64(_mm512_srli_epi64(word, 1), carry));
        carry = nextcarry;
    }

    return;
}

void vec_bignum52_mask_sub(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, uint32_t wmask)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    if (wmask == 0)
        return;

    // subtract the selected elements ('1' in the mask)
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_mask_store_epi64(c->data + i * VECLEN, (__mmask8)wmask,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    if (carry)
    {
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi64(c->data + i * VECLEN, (__mmask8)wmask & carry, _mm512_set1_epi64(0));
    }

    return;
}

void vec_bignum52_add(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    uint32_t NWORDS = MAX(a->size, b->size);
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    // subtract the selected elements ('1' in the mask)
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    if (carry)
    {
        c->size = NWORDS + 1;
        // subtract any final borrows that exceed the size of b.
        _mm512_mask_store_epi64(c->data + i * VECLEN, carry, _mm512_set1_epi64(1));
    }

    return;
}

void vec_bignum52_sub(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    uint32_t NWORDS = MAX(a->size, b->size);
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    // subtract the selected elements ('1' in the mask)
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    bvec = _mm512_xor_si512(bvec, bvec);

    while ((carry > 0) && (i < c->WORDS_ALLOC))
    {
        //printf("sub writing word %d\n", i);
        // subtract any final borrows that exceed the size of b.
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN, cvec);
        i++;
    }

    return;
}

void vec_bignum52_add_1(vec_bignum_t* a, base_t* b, vec_bignum_t* c)
{
    // assumptions:
    // a, c are of length VECLEN * NWORDS
    // b is vector of 1-word base_t's
    // s1 is of length VECLEN
    // a, b, c, and s1 are aligned
    // a and b are both positive
    // a >= b
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;

    bvec = _mm512_load_epi64(b);
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64(VEC_MAXDIGIT), cvec));
    }

    _mm512_store_epi64(c->data + i * VECLEN, _mm512_maskz_set1_epi64(carry, 1));

    return;
}







#endif