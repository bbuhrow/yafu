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

#include "avx_ecm.h"


// define function pointers to the type of reduction needed
void(*vecmulmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecsqrmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecaddmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecsubmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);
void(*vecaddsubmod_ptr)(vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_bignum_t*, vec_monty_t*);


vec_bignum_t * vecInit(uint32_t words)
{
    int i;
    size_t sz = VECLEN * (2 * words + 4);
    vec_bignum_t *n;
    n = (vec_bignum_t *)malloc(sizeof(vec_bignum_t));

    n->data = (base_t *)xmalloc_align(sz * sizeof(base_t));
    if (n->data == NULL)
    {
        printf("could not allocate memory\n");
        exit(2);
    }

    for (i = 0; i < sz; i++)
    {
        n->data[i] = 0;
    }
    n->size = 1;
    n->WORDS_ALLOC = words;
    return n;
}

void vecCopy(vec_bignum_t * src, vec_bignum_t * dest)
{
    //physically copy the digits of u into the digits of v
    int su = VECLEN * (2 * src->WORDS_ALLOC + 1);

    memcpy(dest->data, src->data, su * sizeof(base_t));
    dest->size = src->size; // = NWORDS;
    dest->WORDS_ALLOC = src->WORDS_ALLOC;
    return;
}

void vecCopyn(vec_bignum_t * src, vec_bignum_t * dest, int size)
{
    //physically copy the digits of u into the digits of v
    int su = VECLEN * size;

    memcpy(dest->data, src->data, su * sizeof(base_t));
    dest->size = size;
    dest->WORDS_ALLOC = src->WORDS_ALLOC;
    return;
}

void vecClear(vec_bignum_t *n)
{
    memset(n->data, 0, VECLEN * (2 * n->WORDS_ALLOC + 1) * sizeof(base_t));
    return;
}

void vecFree(vec_bignum_t *n)
{
    align_free(n->data);
    free(n);
}

void copy_vec_lane(vec_bignum_t *src, vec_bignum_t *dest, int num, int size)
{
    int j;

    for (j = 0; j < size; j++)
    {
        dest->data[num + j * VECLEN] = src->data[num + j * VECLEN];
    }

    return;
}

vec_monty_t* vec_monty_alloc(uint32_t words)
{
    int i;
    vec_monty_t *mdata = (vec_monty_t *)malloc(sizeof(vec_monty_t));

    mpz_init(mdata->nhat);
    mpz_init(mdata->rhat);
    mpz_init(mdata->gmp_t1);
    mpz_init(mdata->gmp_t2);
    mdata->r = vecInit(words);
    mdata->n = vecInit(words);
    mdata->vnhat = vecInit(words);
    mdata->vrhat = vecInit(words);
    mdata->rmask = vecInit(words);
    mdata->one = vecInit(words);
    mdata->mtmp1 = vecInit(words);
    mdata->mtmp2 = vecInit(words);
    mdata->mtmp3 = vecInit(words);
    mdata->mtmp4 = vecInit(3*words);

    mdata->g = (vec_bignum_t **)malloc((1 << MAX_WINSIZE) * sizeof(vec_bignum_t *));
    mdata->g[0] = vecInit(words);

    for (i = 1; i < (1 << MAX_WINSIZE); i++)
    {
        mdata->g[i] = vecInit(words);
    }

    mdata->vrho = (base_t *)xmalloc_align(VECLEN * sizeof(base_t));
    mdata->vnbits = (int*)xmalloc_align(VECLEN * sizeof(int));

    return mdata;
}

void vec_monty_free(vec_monty_t *mdata)
{
    int i;

    vecFree(mdata->mtmp1);
    vecFree(mdata->mtmp2);
    vecFree(mdata->mtmp3);
    vecFree(mdata->mtmp4);
    vecFree(mdata->one);
    vecFree(mdata->r);
    vecFree(mdata->n);
    vecFree(mdata->vnhat);
    vecFree(mdata->vrhat);
    vecFree(mdata->rmask);
    mpz_clear(mdata->nhat);
    mpz_clear(mdata->rhat);
    mpz_clear(mdata->gmp_t1);
    mpz_clear(mdata->gmp_t2);

    align_free(mdata->vrho);
    align_free(mdata->vnbits);

    for (i = 0; i < (1 << MAX_WINSIZE); i++)
    {
        vecFree(mdata->g[i]);
    }
    free(mdata->g);

    return;
}


int get_winsize(int bits)
{
    // the window size is based on minimizing the total number of multiplications
    // in the windowed exponentiation.  experiments show that this is best;
    // the growing size of the table doesn't change the calculus, at least
    // on the KNL.
    int size;
    int muls;
    int minmuls = 99999999;
    int minsize = 4;

    for (size = 2; size <= 8; size++)
    {
        muls = (bits / size) + (1 << size);
        if (muls < minmuls)
        {
            minmuls = muls;
            minsize = size;
        }
    }

    return minsize;
}

int get_bitwin(vec_bignum_t* e, int bitloc, int winsize, int lane, int winmask)
{
    int bstr;
    int bitstart = (bitloc - winsize + 1);
    int word = bitloc / DIGITBITS;
    int word2 = bitstart / DIGITBITS;

    bitstart = bitstart % DIGITBITS;

    if (word == word2)
    {
        bstr = (e->data[lane + word * VECLEN] >> bitstart) & winmask;
    }
    else
    {
        int upperbits = (bitloc % DIGITBITS) + 1;

        bstr = (e->data[lane + word2 * VECLEN] >> bitstart);
        bstr |= ((e->data[lane + word * VECLEN]) << (winsize - upperbits));
        bstr &= winmask;
    }

    return bstr;
}


