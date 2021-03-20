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

#define INV_2_POW_64 5.4210108624275221700372640043497e-20

#define DO_STAGE2_INV
#define DEFINED 1
#define MAX_WINSIZE 8
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
} vec_monty_t;

void print_hexbignum(vec_bignum_t *a, const char *pre);
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
void vecmulmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecsqrmod52(vec_bignum_t *a, vec_bignum_t *c, vec_bignum_t *n, vec_bignum_t *s, vec_monty_t *mdata);
void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata);
void vecsubmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata);
void vecsubmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void vecaddmod52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *c, vec_monty_t* mdata);
void vecaddmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata);
void vec_simul_addsub52(vec_bignum_t *a, vec_bignum_t *b, vec_bignum_t *sum, vec_bignum_t *diff, 
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

void(*vecmulmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
void(*vecsqrmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t *);
void(*vecaddmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t*);
void(*vecsubmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t*);
void(*vecaddsubmod_ptr)(vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_bignum_t *, vec_monty_t*);

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




