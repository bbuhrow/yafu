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

vec_bignum_t * vecInit(void)
{
    int i;
    size_t sz = VECLEN * (2 * NWORDS + 4);
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

    return n;
}

void vecCopy(vec_bignum_t * src, vec_bignum_t * dest)
{
    //physically copy the digits of u into the digits of v
    int su = VECLEN * (2 * NWORDS + 1);

    memcpy(dest->data, src->data, su * sizeof(base_t));
    dest->size = src->size; // = NWORDS;
    return;
}

void vecCopyn(vec_bignum_t * src, vec_bignum_t * dest, int size)
{
    //physically copy the digits of u into the digits of v
    int su = VECLEN * size;

    memcpy(dest->data, src->data, su * sizeof(base_t));
    dest->size = size;
    return;
}

void vecClear(vec_bignum_t *n)
{
    memset(n->data, 0, VECLEN * (2 * NWORDS + 1) * sizeof(base_t));
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

vec_monty_t* vec_monty_alloc(void)
{
    int i;
    vec_monty_t *mdata = (vec_monty_t *)malloc(sizeof(vec_monty_t));

    mpz_init(mdata->nhat);
    mpz_init(mdata->rhat);
    mpz_init(mdata->gmp_t1);
    mpz_init(mdata->gmp_t2);
    mdata->r = vecInit();
    mdata->n = vecInit();
    mdata->vnhat = vecInit();
    mdata->vrhat = vecInit();
    mdata->rmask = vecInit();
    mdata->one = vecInit();
    mdata->mtmp1 = vecInit();
    mdata->mtmp2 = vecInit();
    mdata->mtmp3 = vecInit();
    mdata->mtmp4 = vecInit();

    mdata->g = (vec_bignum_t **)malloc((1 << MAX_WINSIZE) * sizeof(vec_bignum_t *));
    mdata->g[0] = vecInit();

    for (i = 1; i < (1 << MAX_WINSIZE); i++)
    {
        mdata->g[i] = vecInit();
    }

    mdata->vrho = (base_t *)xmalloc_align(VECLEN * sizeof(base_t));

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

    for (i = 0; i < (1 << MAX_WINSIZE); i++)
    {
        vecFree(mdata->g[i]);
    }
    free(mdata->g);

    return;
}

