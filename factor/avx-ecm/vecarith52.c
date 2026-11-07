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
#include "monty.h"
#include <math.h>

// The old fixed416..1040 / BFIPS mul / sqr cores live in
// vecarith52_oldcores.c, outside the builds (kept for A/B tests in
// factor/avx-ecm/bench and arith/rns): AVX-ECM uses the constant-S template
// cores (arith/vecmont/vecarith52_tmpl.c) and the RNS backend.

#ifdef USE_AVX512F
#include <immintrin.h>

#define USE_AMM 1
//#define DEBUG_MUL
//#define DO_INTERMEDIATE_CARRYPROP


void vecmul52_1(vec_bignum_t* a, __m512i b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // unsigned multiply a * b, where b is a single base_t digit.
    // resulting in a NWORDS+1 size result.
    int i;
    int NWORDS = (int)n->WORDS_ALLOC;
    // needed in loops
    __m512i a0, carry;                                     // 4
    UNUSED_VAR __m512d prod1_hd;                 // 23
    UNUSED_VAR __m512d prod1_ld, prod2_ld;        // 28
    UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1;
    __m512i zero = _mm512_set1_epi64(0);

    carry = zero;
    a0 = _mm512_load_epi64(a->data);
    VEC_MUL_LOHI_PD(a0, b, acc_e0, acc_e1);
    carry = acc_e1;
    _mm512_store_epi64(c->data, acc_e0);

    for (i = 1; i < NWORDS; i++)
    {
        a0 = _mm512_load_epi64(a->data + i * VECLEN_I);
        VEC_MUL_LOHI_PD(a0, b, acc_e0, acc_e1);
        acc_e0 = _mm512_add_epi64(acc_e0, carry);
        VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
        carry = acc_e1;
        _mm512_store_epi64(c->data + i * VECLEN_I, acc_e0);
    }
    _mm512_store_epi64(c->data + i * VECLEN_I, carry);

    c->size = NWORDS + 1;
    return;
}

void vecaddmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    int NWORDS = (int)mdata->NWORDS; // a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 mask = 0;
    UNUSED_VAR __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(b->data + i * VECLEN_I);
        
        //cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        __m512i t = _mm512_add_epi64(avec, bvec);
        t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_cmpgt_epu64_mask(t, lomask);
        cvec = _mm512_and_epi64(t, lomask);

        _mm512_store_epi64(c->data + i * VECLEN_I, cvec);
    }

    mask = carry;	// sub mask
    mask2 = mask;	// keep looking mask

#ifndef USE_AMM
    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN_I);
        // compare those that have not already been decided using the mask
        mask |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_GT);
        mask2 |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (mask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    mask |= (~mask2);
#endif

    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN_I);
        
        //bvec = _mm512_mask_sbb_epi52(cvec, mask, carry, nvec, &carry);
        __m512i t = _mm512_mask_sub_epi64(cvec, mask, cvec, nvec);
        carryout = _mm512_mask_cmpgt_epu64_mask(mask, nvec, cvec);
        __m512i t2 = _mm512_mask_sub_epi64(cvec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _kor_mask8(carryout, _mm512_cmpgt_epu64_mask(t2, t));
        bvec = _mm512_and_epi64(t2, lomask);

        _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
    }
    return;
}

void vecsubmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // s1 is of length VECLEN_I
    // a, b, c, n, and s1 are aligned
    // a and b are both positive
    // a >= b
    // n is the montgomery base
    int i;
    int NWORDS = (int)mdata->NWORDS; // a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 mask = 0;
    //__mmask8 mask2 = 0;
    __m512i nvec;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);

    // subtract
    carry = 0;
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(b->data + i * VECLEN_I);
        //cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);

        __m512i t = _mm512_sub_epi64(avec, bvec);
        carryout = _mm512_cmpgt_epu64_mask(bvec, avec);
        __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _kor_mask8(carryout, _mm512_cmpgt_epu64_mask(t2, t));
        cvec = _mm512_and_epi64(t2, lomask);

        _mm512_store_epi64(c->data + i * VECLEN_I, cvec);
    }

    // if we had a final carry, then b was bigger than a so we need to re-add n.
    mask = carry;
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        avec = _mm512_load_epi64(c->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN_I);
        //cvec = _mm512_mask_adc_epi52(avec, mask, carry, nvec, &carry);

        __m512i t = _mm512_mask_add_epi64(avec, mask, avec, nvec);
        t = _mm512_mask_add_epi64(avec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_mask_cmpgt_epu64_mask(mask, t, lomask);
        cvec = _mm512_and_epi64(t, lomask);

        _mm512_store_epi64(c->data + i * VECLEN_I, cvec);
    }

    // if we did not have a final carry, then we still are not in
    // positive territory.  add n again.  Necessary with AMM.
    mask = (~carry) & mask;
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        avec = _mm512_load_epi64(c->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN_I);
        //cvec = _mm512_mask_adc_epi52(avec, mask, carry, nvec, &carry);
    
        __m512i t = _mm512_mask_add_epi64(avec, mask, avec, nvec);
        t = _mm512_mask_add_epi64(avec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_mask_cmpgt_epu64_mask(mask, t, lomask);
        cvec = _mm512_and_epi64(t, lomask);
    
        _mm512_store_epi64(c->data + i * VECLEN_I, cvec);
    }
    return;
}

void vecsignedaddmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_bignum_t *n)
{
    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    int NWORDS = (int)a->WORDS_ALLOC;
    __mmask8 borrow = 0;
    __mmask8 carry = 0;
    __mmask8 mask = 0;
    UNUSED_VAR __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;

    // this is nominally an add, but if a or b are negative then
    // it's a subtract.  Since the vector can contain both cases
    // at the same time we have to perform both operations with 
    // the appropriate masks.
    // if both positive or negative then we add.
    mask = (__mmask8)(a->signmask ^ b->signmask);

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(b->data + i * VECLEN_I);
        // this won't work, because the sbb will replace the lanes
        // where mask=1 with avec.
        cvec = _mm512_mask_adc_epi52(avec, mask, carry, bvec, &carry);
        cvec = _mm512_mask_sbb_epi52(avec, ~mask, borrow, bvec, &borrow);
        _mm512_store_epi64(c->data + i * VECLEN_I, cvec);
    }

    mask = carry;	// sub mask
    mask2 = mask;	// keep looking mask

    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(n->data + i * VECLEN_I);
        // compare those that have not already been decided using the mask
        mask |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_GT);
        mask2 |= _mm512_mask_cmp_epu64_mask(~mask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (mask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    mask |= (~mask2);

    // subtract n from c when c is not less than n, as indicated by a 1 bit in mask
    carry = 0;
    for (i = 0; (i < NWORDS) && (mask > 0); i++)
    {
        cvec = _mm512_load_epi64(c->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(n->data + i * VECLEN_I);
        bvec = _mm512_mask_sbb_epi52(cvec, mask, carry, nvec, &carry);
        _mm512_store_epi64(c->data + i * VECLEN_I, bvec);
    }
    return;
}

void vec_simul_addsub52(vec_bignum_t *a, vec_bignum_t *b, 
    vec_bignum_t *sum, vec_bignum_t *diff, 
    vec_monty_t* mdata)
{
    // awkward to implement the 2*n modsub addition 
    // that is sometimes necessary in AMM.  Just do individual add/sub.
    vecaddmod52(a, b, sum, mdata);
    vecsubmod52(a, b, diff, mdata);
    return;


    // assumptions:
    // a, b, c are of length VECLEN_I * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    // produce sum = a + b and diff = a - b at the same time which
    // saves 3N loads (only have to load a,b, and n once)
    int i;
    int NWORDS = (int)mdata->NWORDS; // a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 borrow = 0;
    __mmask8 cmask = 0;
    UNUSED_VAR __mmask8 cmask2 = 0;
    __mmask8 bmask = 0;
    //__mmask8 bmask2 = 0;
    __m512i avec;
    __m512i bvec;
    //__m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffff);

    for (i = 0; i < NWORDS; i++)
    {
        // add/sub
        avec = _mm512_load_epi64(a->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(b->data + i * VECLEN_I);

        _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, avec, bvec);

        _mm512_store_epi64(sum->data + i * VECLEN_I, avec);
        _mm512_store_epi64(diff->data + i * VECLEN_I, bvec);
    }

    cmask = carry;	    // result too big, need to subtract n
    bmask = borrow;     // result too small, need to add n
    cmask2 = cmask;	    // keep looking mask for add

#ifndef USE_AMM
    // compare.  the initial mask is equal to the addition carry
    // because if the most significant word has a carry, then the
    // result is bigger than n.
    for (i = NWORDS - 1; i >= 0; i--)
    {
        cvec = _mm512_load_epi64(sum->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN_I);
        // compare those that have not already been decided using the mask
        cmask |= _mm512_mask_cmp_epu64_mask(~cmask2, cvec, nvec, _MM_CMPINT_GT);
        cmask2 |= _mm512_mask_cmp_epu64_mask(~cmask2, cvec, nvec, _MM_CMPINT_LT);

        // decided all of them, stop comparing.
        if (cmask2 == 0xff)
        {
            break;
        }
    }

    // check for equal as well by flipping mask bits that have still
    // not been decided (i.e., are equal)
    cmask |= (~cmask2);
#endif

    carry = 0;
    borrow = 0;

    for (i = 0; i < NWORDS; i++)
    {
        // conditional sub/add
        avec = _mm512_load_epi64(sum->data + i * VECLEN_I);
        bvec = _mm512_load_epi64(diff->data + i * VECLEN_I);
        nvec = _mm512_load_epi64(mdata->n->data + i * VECLEN_I);

        _mm512_mask_adcsbb_epi52(bvec, avec, bmask, cmask, carry, borrow, nvec, lomask, bvec, avec);

        _mm512_store_epi64(sum->data + i * VECLEN_I, avec);
        _mm512_store_epi64(diff->data + i * VECLEN_I, bvec);
    }

    return;
}

#else

void vecaddmod52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vecsubmod52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vec_simul_addsub52(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff, vec_monty_t* mdata)
{
    return;
}

void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    return;
}

void vecaddmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vecsubmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    return;
}

void vec_simul_addsub52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff, vec_monty_t* mdata)
{
    return;
}

void vecmodexp52(vec_bignum_t* d, vec_bignum_t* b, vec_bignum_t* e, vec_bignum_t* m,
    vec_bignum_t* s, vec_bignum_t* one, vec_monty_t* mdata)
{

    return;
}


#endif
