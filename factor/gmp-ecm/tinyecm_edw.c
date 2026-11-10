/*
tinyecm_edw.c: tinyecm with stage 1 on twisted Edwards curves.

Copyright (c) 2014, Ben Buhrow and (c) 2022, Jeff Hurchalla (the parts
taken from tinyecm.c / microecm.c: the 104-bit vector arithmetic, the
x-only stage 2 and its pairing tables).  Dual-licensed, at your option,
under the FreeBSD license or the Mozilla Public License 2.0, exactly as
tinyecm.c; see the full license texts there.

The curve code (GKL family (iv) curves built without inversions, the
small-constant table curves, Edwards stage 1 on precomputed chains, the
hand-off to the Montgomery image, stage 2) is edw_kern.h, instantiated
twice here:
  - scalar, 128-bit residues (portable CIOS Montgomery, R = 2^128);
  - x8 with AVX-512, 104-bit residues as two 52-bit digits (IFMA madd52
    or FP52 products, R = 2^104).
The drivers keep tinyecm's interfaces, B1 levels, curve counts, stage-2
tables (B2 = 25 B1) and stage-2 accumulation rules:
  - scalar: a product that would make acc 0 mod n ends that loop;
  - x8: such a product freezes the lane's acc for the rest of stage 2.
See arith/micro_tiny_edwards_plan.md.
*/

#include "tinyecm_edw.h"
#include "edw_tiny_curves.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

#ifndef ALIGNED_MEM
#if defined(_MSC_VER)
#define ALIGNED_MEM __declspec(align(64))
#else
#define ALIGNED_MEM __attribute__((aligned(64)))
#endif
#endif

#if defined(_MSC_VER)
#define TEDW_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define TEDW_INLINE inline __attribute__((always_inline))
#else
#define TEDW_INLINE __inline
#endif

/* ================================================================ scalar */
typedef struct { uint64_t lo, hi; } tedw_u128;

/* hi:lo = a b + c + d (never overflows 128 bits) */
static TEDW_INLINE uint64_t tedw_mac(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t *hi)
{
#if defined(__SIZEOF_INT128__)
    unsigned __int128 t = (unsigned __int128)a * b + c + d;
    *hi = (uint64_t)(t >> 64);
    return (uint64_t)t;
#else
    uint64_t h;
    uint64_t l = _umul128(a, b, &h);
    l += c;
    h += (l < c);
    l += d;
    h += (l < d);
    *hi = h;
    return l;
#endif
}

typedef struct {
    tedw_u128 n;
    uint64_t rho;                               // -1/n mod 2^64
    tedw_u128 one, Rsqr;                        // R mod n, R^2 mod n (R = 2^128)
    tedw_u128 c6, c10, c16, c32, c48;
} tedw_ctx;

static TEDW_INLINE int tedw_geq(tedw_u128 a, tedw_u128 b)
{
    return (a.hi > b.hi) || ((a.hi == b.hi) && (a.lo >= b.lo));
}

static TEDW_INLINE tedw_u128 tedw_sub_raw(tedw_u128 a, tedw_u128 b)
{
    tedw_u128 r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (uint64_t)(a.lo < b.lo);
    return r;
}

static TEDW_INLINE tedw_u128 tedw_add_raw(tedw_u128 a, tedw_u128 b, int *carry)
{
    tedw_u128 r;
    r.lo = a.lo + b.lo;
    uint64_t c1 = (uint64_t)(r.lo < a.lo);
    uint64_t s = a.hi + b.hi;
    int c2 = (s < a.hi);
    r.hi = s + c1;
    *carry = c2 | (r.hi < s);
    return r;
}

/* a - b with the borrow out (branch-free: the conditional corrections in
   these 2-word routines are unpredictable once n > 2^64) */
static TEDW_INLINE tedw_u128 tedw_sub_b(tedw_u128 a, tedw_u128 b, uint64_t *borrow)
{
    tedw_u128 r;
    uint64_t b1 = (uint64_t)(a.lo < b.lo);
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - b1;
    *borrow = (uint64_t)(a.hi < b.hi) | ((uint64_t)(a.hi == b.hi) & b1);
    return r;
}

/* mask ? a : b, with mask all ones or zero */
static TEDW_INLINE tedw_u128 tedw_sel(uint64_t mask, tedw_u128 a, tedw_u128 b)
{
    tedw_u128 r;
    r.lo = (a.lo & mask) | (b.lo & ~mask);
    r.hi = (a.hi & mask) | (b.hi & ~mask);
    return r;
}

/* CIOS Montgomery product, two words; inputs with u v < n R */
static TEDW_INLINE tedw_u128 tedw_mulredc(tedw_u128 u, tedw_u128 v, const tedw_ctx *c)
{
    uint64_t t0, t1, t2, cy, m, h;
    t0 = tedw_mac(u.lo, v.lo, 0, 0, &cy);
    t1 = tedw_mac(u.lo, v.hi, cy, 0, &t2);
    m = t0 * c->rho;
    (void)tedw_mac(m, c->n.lo, t0, 0, &cy);
    t0 = tedw_mac(m, c->n.hi, t1, cy, &cy);
    t1 = t2 + cy;
    t2 = (uint64_t)(t1 < cy);
    t0 = tedw_mac(u.hi, v.lo, t0, 0, &cy);
    t1 = tedw_mac(u.hi, v.hi, t1, cy, &h);
    t2 += h;
    m = t0 * c->rho;
    (void)tedw_mac(m, c->n.lo, t0, 0, &cy);
    t0 = tedw_mac(m, c->n.hi, t1, cy, &cy);
    t1 = t2 + cy;
    t2 = (uint64_t)(t1 < cy);
    tedw_u128 r;
    r.lo = t0;
    r.hi = t1;
    uint64_t bw;
    tedw_u128 d = tedw_sub_b(r, c->n, &bw);
    return tedw_sel((uint64_t)0 - (t2 | (bw ^ 1)), d, r);   // r - n unless r < n
}

static TEDW_INLINE tedw_u128 tedw_addmod(tedw_u128 x, tedw_u128 y, const tedw_ctx *c)
{
    int carry;
    uint64_t bw;
    tedw_u128 s = tedw_add_raw(x, y, &carry);
    tedw_u128 d = tedw_sub_b(s, c->n, &bw);
    return tedw_sel((uint64_t)0 - ((uint64_t)carry | (bw ^ 1)), d, s);
}

static TEDW_INLINE tedw_u128 tedw_submod(tedw_u128 x, tedw_u128 y, const tedw_ctx *c)
{
    uint64_t bw;
    int carry;
    tedw_u128 r = tedw_sub_b(x, y, &bw);
    tedw_u128 nm;
    nm.lo = c->n.lo & ((uint64_t)0 - bw);
    nm.hi = c->n.hi & ((uint64_t)0 - bw);
    return tedw_add_raw(r, nm, &carry);
}

static TEDW_INLINE tedw_u128 tedw_zero(void)
{
    tedw_u128 z;
    z.lo = z.hi = 0;
    return z;
}

static TEDW_INLINE int tedw_iszero(tedw_u128 a)
{
    return (a.lo | a.hi) == 0;
}

static void tedw_to_mpz(mpz_t out, tedw_u128 a)
{
    mpz_set_ui(out, (unsigned long)(a.hi >> 32));
    mpz_mul_2exp(out, out, 32);
    mpz_add_ui(out, out, (unsigned long)(a.hi & 0xffffffffULL));
    mpz_mul_2exp(out, out, 32);
    mpz_add_ui(out, out, (unsigned long)(a.lo >> 32));
    mpz_mul_2exp(out, out, 32);
    mpz_add_ui(out, out, (unsigned long)(a.lo & 0xffffffffULL));
}

static tedw_u128 tedw_from_mpz(const mpz_t a)
{
    tedw_u128 r;
    mpz_t t;
    mpz_init_set(t, a);
    r.lo = (uint64_t)mpz_get_ui(t) & 0xffffffffULL;
    mpz_tdiv_q_2exp(t, t, 32);
    r.lo |= ((uint64_t)mpz_get_ui(t) & 0xffffffffULL) << 32;
    mpz_tdiv_q_2exp(t, t, 32);
    r.hi = (uint64_t)mpz_get_ui(t) & 0xffffffffULL;
    mpz_tdiv_q_2exp(t, t, 32);
    r.hi |= ((uint64_t)mpz_get_ui(t) & 0xffffffffULL) << 32;
    mpz_clear(t);
    return r;
}

static uint64_t tedw_inverse64(uint64_t a)
{
    uint64_t x0 = (3 * a) ^ 2;
    uint64_t y = 1 - a * x0;
    uint64_t x1 = x0 * (1 + y);
    y *= y;
    uint64_t x2 = x1 * (1 + y);
    y *= y;
    uint64_t x3 = x2 * (1 + y);
    y *= y;
    return x3 * (1 + y);
}

/* n odd, < 2^128 */
static void tedw_ctx_init(tedw_ctx *c, const mpz_t n)
{
    mpz_t t, nn;
    mpz_init(t);
    mpz_init_set(nn, n);
    c->n = tedw_from_mpz(nn);
    c->rho = (uint64_t)0 - tedw_inverse64(c->n.lo);
    mpz_set_ui(t, 1);
    mpz_mul_2exp(t, t, 128);
    mpz_mod(t, t, nn);
    c->one = tedw_from_mpz(t);
    mpz_set_ui(t, 1);
    mpz_mul_2exp(t, t, 256);
    mpz_mod(t, t, nn);
    c->Rsqr = tedw_from_mpz(t);
    tedw_u128 c2 = tedw_addmod(c->one, c->one, c);
    tedw_u128 c4 = tedw_addmod(c2, c2, c);
    tedw_u128 c8 = tedw_addmod(c4, c4, c);
    c->c16 = tedw_addmod(c8, c8, c);
    c->c32 = tedw_addmod(c->c16, c->c16, c);
    c->c48 = tedw_addmod(c->c32, c->c16, c);
    c->c6 = tedw_addmod(c4, c2, c);
    c->c10 = tedw_addmod(c8, c2, c);
    mpz_clear(t);
    mpz_clear(nn);
}

#define EK_T tedw_u128
#define EK_CTX tedw_ctx
#define EK_FN(x) tedw_##x
#define EK_MUL(a, b) tedw_mulredc((a), (b), c)
#define EK_SQR(a) tedw_mulredc((a), (a), c)
#define EK_ADD(a, b) tedw_addmod((a), (b), c)
#define EK_SUB(a, b) tedw_submod((a), (b), c)
#define EK_ZERO tedw_zero()
#define EK_ACC_DECL
#define EK_ACC(acc, t) { tedw_u128 tmp_ = EK_MUL(acc, t); if (tedw_iszero(tmp_)) break; acc = tmp_; }
#include "edw_kern.h"
#undef EK_T
#undef EK_CTX
#undef EK_FN
#undef EK_MUL
#undef EK_SQR
#undef EK_ADD
#undef EK_SUB
#undef EK_ZERO
#undef EK_ACC_DECL
#undef EK_ACC

/* tinyecm.c's pre-paired stage-2 sequences (B2 = 25 B1) */
static const int tedw_numb1_27 = 76;
static const uint8_t tedw_b1_27[76] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7 };
static const int tedw_numb1_47 = 121;
static const uint8_t tedw_b1_47[121] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7,1,19,37,47,
0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,
37,23,19,13,7,1,29,31,41,59,0,59,49,47,41,23,
19,17,13,7,1,43,53,0,59 };
static const int tedw_numb1_70 = 175;
static const uint8_t tedw_b1_70[175] = {
31,41,43,47,49,0,59,47,43,41,37,31,29,19,13,7,
1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,
0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,
41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,
31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,
29,31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,
53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,
53,49,47,43,31,29,23,11,17,0,47,43,41,37,31,23,
19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,
7,17,29,0,53,47,43,41,17,13,11,1,23,31,37 };
static const int tedw_numb1_85 = 225;
static const uint8_t tedw_b1_85[225] = {
	1,53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
	31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
	0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
	31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
	29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
	17,13,11,1,23,31,37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,
	43,41,31,17,7,1,11,13,19,29 };
static const int tedw_numb1_125 = 319;
static const uint8_t tedw_b1_125[319] = {
	23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
	29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,41,31,23,
	17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,31,41,59,0,59,49,
	47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,29,23,11,17,0,47,
	43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,17,13,11,1,23,31,
	37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,43,41,31,17,7,1,
	11,13,19,29,0,59,53,49,47,37,29,11,13,17,23,31,0,59,43,41,37,29,23,17,13,1,31,47,0,59,53,49,47,41,37,
	31,19,13,7,11,17,29,43,0,47,29,19,11,7,1,41,43,59,0,53,49,37,23,13,11,7,1,17,19,29,41,43,59,0,59,49,
	41,37,23,13,1,7,11,29,43,47,53,0,59,53,49,31,23,13,7,1,17,29,43,47,0,59,31,29,19,11,7,37,49,53 };
static const int tedw_numb1_165 = 425;
static const uint8_t tedw_b1_165[425] = {
	13,7,1,11,19,47,59,0,59,49,43,37,31,29,23,19,17,7,11,13,47,53,0,53,47,41,37,31,23,19,11,1,13,29,43,
	59,0,53,49,41,37,31,19,17,1,7,23,29,47,59,0,59,53,47,43,41,29,19,17,13,7,1,23,31,49,0,53,47,41,37,29,
	23,19,11,7,17,31,43,49,59,0,47,43,41,37,23,19,17,13,7,11,29,53,0,53,49,43,37,29,23,11,7,1,13,19,31,41,
	0,53,49,47,43,37,31,23,17,11,13,41,0,59,47,43,37,31,29,23,11,1,17,19,41,0,59,53,19,13,7,1,29,43,47,49,
	0,53,49,47,41,29,19,17,13,11,7,1,23,31,43,59,0,53,49,41,37,23,19,13,11,7,1,17,43,47,0,47,43,41,31,19,
	17,7,1,13,37,49,0,59,49,37,29,13,1,7,11,17,19,41,47,53,0,49,47,31,29,7,1,13,17,19,23,37,59,0,47,37,31,
	19,17,13,11,1,29,41,43,53,0,59,41,17,13,7,1,19,23,31,47,49,53,0,59,53,47,43,31,29,7,1,11,17,37,41,49,
	0,49,43,37,23,19,13,1,7,17,0,59,49,41,37,31,29,23,1,11,13,53,0,53,43,41,37,29,23,17,13,11,7,1,19,31,49,
	0,53,43,31,29,23,19,17,1,13,37,41,59,0,53,43,37,31,23,13,1,17,29,59,0,59,49,41,37,23,19,11,1,7,29,0,59,
	43,17,13,11,1,7,23,29,37,41,49,0,49,47,43,41,29,1,7,13,19,23,31,59,0,59,49,47,31,29,13,7,37,41,43,0,49,
	41,29,23,13,11,7,1,17,19,31,43,53,0,53,47,43,37,29,23,17,1,11,13,31,41,49,59,0,53,47,41,19,13,11,1,17,
	23,43,0,53,49,47,37,23,19,11,7,17,29,31,43,0,53,31,19,17,13,7,1,29,37,59 };
static const int tedw_numb1_205 = 511;
static const uint8_t tedw_b1_205[511] = {
	1,23,41,0,59,53,49,47,37,23,19,17,13,1,7,29,43,0,53,49,41,31,29,19,17,11,7,1,13,37,59,0,49,47,29,23,
	13,7,1,17,31,37,43,0,59,49,47,43,37,31,29,17,13,7,1,11,19,53,0,59,53,49,41,37,23,13,1,11,17,19,29,43,
	47,0,53,49,47,43,23,19,11,1,7,17,37,41,0,59,53,41,37,31,29,19,17,11,1,13,43,47,0,53,47,41,19,17,7,1,
	11,23,31,43,59,0,59,53,41,31,13,11,7,1,17,29,37,0,49,43,37,29,11,1,13,17,19,23,41,0,59,49,47,43,41,37,
	31,19,7,1,13,23,29,53,0,53,49,43,41,37,31,29,23,13,7,17,19,47,59,0,49,47,37,29,23,17,11,7,13,19,31,41,
	53,0,59,43,29,23,19,17,13,11,1,41,0,59,37,31,23,17,13,11,7,1,19,29,43,53,0,49,47,43,41,31,19,17,1,7,11,
	13,23,0,47,43,37,29,13,11,7,1,17,19,23,31,59,0,59,37,31,29,23,19,13,1,7,11,41,47,53,0,53,49,43,31,23,
	17,13,41,59,0,59,53,31,19,17,1,7,11,23,37,47,49,0,59,53,47,43,41,37,31,23,19,17,11,1,0,59,53,49,47,31,
	17,13,7,1,11,29,37,0,53,43,31,17,13,7,1,29,41,49,0,53,49,41,29,23,11,7,1,19,31,47,0,47,43,41,29,23,19,
	7,1,11,49,0,59,31,29,23,17,11,7,1,13,41,43,0,59,43,37,17,1,7,11,13,19,41,49,0,59,53,43,41,37,31,29,23,
	13,11,1,47,0,59,53,47,31,19,17,13,1,7,11,29,37,43,49,0,49,43,41,31,17,13,7,11,23,37,53,0,53,49,41,23,
	19,13,11,7,1,17,37,59,0,49,47,43,37,31,29,23,1,7,41,0,59,43,41,37,31,17,13,11,7,47,49,0,59,49,47,37,31,
	29,19,17,7,1,0,53,47,37,19,13,1,11,31,41,0,49,47,37,23,17,13,11,7,19,31,53,0,59,53,47,29,13,11,7,1,23,
	41,0,49,47,41,37,19,11,13,17,23,29,31,43,0,59,29,19,13,1,41,43,47,53,0,59,53,43,41,37,23,17,11,7,1,13,
	29,49 };

static int tedw_s2table(uint32_t B1, const uint8_t **barray, int *numb)
{
    switch (B1) {
    case 27: *barray = tedw_b1_27; *numb = tedw_numb1_27; return 1;
    case 47: *barray = tedw_b1_47; *numb = tedw_numb1_47; return 1;
    case 70: *barray = tedw_b1_70; *numb = tedw_numb1_70; return 1;
    case 85: *barray = tedw_b1_85; *numb = tedw_numb1_85; return 1;
    case 125: *barray = tedw_b1_125; *numb = tedw_numb1_125; return 1;
    case 165: *barray = tedw_b1_165; *numb = tedw_numb1_165; return 1;
    case 205: *barray = tedw_b1_205; *numb = tedw_numb1_205; return 1;
    default: return 0;
    }
}

/* table constant (|v| < 2^63) to Montgomery form */
static tedw_u128 tedw_const(int64_t v, const tedw_ctx *c)
{
    tedw_u128 a;
    a.lo = (v < 0) ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    a.hi = 0;
    a = tedw_mulredc(a, c->Rsqr, c);            // |v| R^2 / R, |v| R^2 < n R
    return (v < 0) ? tedw_submod(tedw_zero(), a, c) : a;
}

static void tedw_table_curve(int i, tedw_curve *cv, const tedw_ctx *c)
{
    tedw_u128 v[6];
    int k;
    for (k = 0; k < 6; k++)
        v[k] = tedw_const(edw_tiny_curves[i][k], c);
    tedw_curve_from(v, cv, c);
}

/* [m]Q on E', m >= 1 */
static void tedw_w_start(tedw_wpt *R, uint32_t m, const tedw_ctx *c)
{
    R->X = tedw_const(5, c);
    R->Y = tedw_const(8, c);
    R->Z = c->one;
    int top = 31;
    while (top > 0 && !((m >> top) & 1))
        top--;
    for (int i = top - 1; i >= 0; i--) {
        tedw_w_dbl(R, c);
        if ((m >> i) & 1)
            tedw_w_addq(R, c);
    }
}

/* the curve source of one input: table curves, then GKL m0, m0 + 1, ... */
typedef struct {
    tedw_wpt R;
    uint32_t m0;
    int started, tab, ntab;
} tedw_src;

#ifdef TEDW_STATS
uint64_t tedw_stat_curves, tedw_stat_gkl, tedw_stat_starts;
#endif

static void tedw_next_curve(tedw_src *s, tedw_curve *cv, const tedw_ctx *c)
{
#ifdef TEDW_STATS
    tedw_stat_curves++;
    if (s->tab >= s->ntab) tedw_stat_gkl++;
    if (s->tab >= s->ntab && !s->started) tedw_stat_starts++;
#endif
    if (s->tab < s->ntab) {
        tedw_table_curve(s->tab++, cv, c);
        return;
    }
    if (!s->started) {
        tedw_w_start(&s->R, s->m0, c);
        s->started = 1;
    }
    tedw_gkl_curve(&s->R, cv, c);
    tedw_w_addq(&s->R, c);
}

static int tedw_check(tedw_u128 x, const mpz_t n, mpz_t f)
{
    mpz_t g;
    int r = 0;
    mpz_init(g);
    tedw_to_mpz(g, x);
    mpz_gcd(g, g, n);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
        mpz_set(f, g);
        r = 1;
    }
    mpz_clear(g);
    return r;
}

/* curves from the source; returns 1 if a factor was found (in f) */
static int tedw_curves(const mpz_t n, mpz_t f, uint32_t B1, uint32_t B2, uint32_t curves,
    const tedw_ctx *c, tedw_src *src, int *stage)
{
    edw_chain_t ch;
    const uint8_t *barray;
    int numb;
    uint32_t i;
    *stage = 0;
    if (!edw_chain_for(B1, &ch) || !tedw_s2table(B1, &barray, &numb))
        return 0;
    for (i = 0; i < curves; i++) {
        tedw_curve cv;
        tedw_pt P;
        tedw_next_curve(src, &cv, c);
        tedw_stage1(&cv, &ch, &P, c);
        if (tedw_check(P.Z, n, f)) {
            *stage = 1;
            return 1;
        }
        if (B2 > B1) {
            tedw_u128 acc = tedw_stage2(&P, &cv, B1, barray, numb, c);
            if (tedw_check(acc, n, f)) {
                *stage = 2;
                return 1;
            }
        }
    }
    return 0;
}

static uint32_t tedw_lcg_m0(uint64_t *pran)
{
    *pran = 6364136223846793005ULL * (*pran) + 1442695040888963407ULL;
    return 3 + (uint32_t)((*pran >> 40) & 0xff);
}

int getfactor_tecm_edw(mpz_t n, mpz_t f, int target_bits, uint64_t *pran)
{
    tedw_ctx ctx;
    tedw_src src;
    int bits, arbitrary, stage, found = 0;

    mpz_set_ui(f, 1);
    if (mpz_even_p(n)) {
        mpz_set_ui(f, 2);
        return 1;
    }
    if (mpz_sizeinbase(n, 2) > 128 || mpz_cmp_ui(n, 3) < 0)
        return 0;
    if (target_bits == 0) {
        bits = (int)mpz_sizeinbase(n, 2);
        arbitrary = 1;
    } else {
        bits = target_bits * 2;
        arbitrary = 0;
    }
    tedw_ctx_init(&ctx, n);
    src.m0 = tedw_lcg_m0(pran);
    src.started = 0;
    src.tab = 0;
    src.ntab = EDW_TINY_NCURVES;

    if (arbitrary) {
        // fast attempts at possible small factors (as tecm_dispatch)
        if (tedw_curves(n, f, 70, 25 * 70, 1, &ctx, &src, &stage)) return 1;
        if (bits > 58)
            if (tedw_curves(n, f, 125, 25 * 125, 1, &ctx, &src, &stage)) return 1;
    }
    if (bits <= 48) found = tedw_curves(n, f, 70, 25 * 70, 32, &ctx, &src, &stage);
    else if (bits <= 52) found = tedw_curves(n, f, 85, 25 * 85, 32, &ctx, &src, &stage);
    else if (bits <= 58) found = tedw_curves(n, f, 125, 25 * 125, 32, &ctx, &src, &stage);
    else if (bits <= 62) found = tedw_curves(n, f, 165, 25 * 165, 42, &ctx, &src, &stage);
    else found = tedw_curves(n, f, 205, 25 * 205, 42, &ctx, &src, &stage);
    return found;
}

static void tedw_curve_by_id(uint32_t id, tedw_curve *cv, const tedw_ctx *c)
{
    if (id & TECM_EDW_TABLE) {
        tedw_table_curve((int)(id & ~TECM_EDW_TABLE) % EDW_TINY_NCURVES, cv, c);
    } else {
        tedw_wpt R;
        tedw_w_start(&R, id, c);
        tedw_gkl_curve(&R, cv, c);
    }
}

int tecm_edw_one_curve(mpz_t n, uint32_t id, uint32_t B1, uint32_t B2, mpz_t f, int *stage)
{
    tedw_ctx ctx;
    tedw_src src;
    tedw_ctx_init(&ctx, n);
    mpz_set_ui(f, 1);
    if (id & TECM_EDW_TABLE) {
        src.tab = (int)(id & ~TECM_EDW_TABLE) % EDW_TINY_NCURVES;
        src.ntab = src.tab + 1;
        src.m0 = 3;
    } else {
        src.tab = src.ntab = 0;
        src.m0 = id;
    }
    src.started = 0;
    return tedw_curves(n, f, B1, B2, 1, &ctx, &src, stage);
}

void tecm_edw_stage1_point(mpz_t n, uint32_t id, uint32_t B1, mpz_t X, mpz_t Z)
{
    tedw_ctx ctx;
    tedw_curve cv;
    edw_chain_t ch;
    tedw_pt P;
    tedw_u128 one;
    tedw_ctx_init(&ctx, n);
    tedw_curve_by_id(id, &cv, &ctx);
    mpz_set_ui(X, 0);
    mpz_set_ui(Z, 0);
    if (!edw_chain_for(B1, &ch))
        return;
    tedw_stage1(&cv, &ch, &P, &ctx);
    one.lo = 1;
    one.hi = 0;
    tedw_to_mpz(X, tedw_mulredc(P.X, one, &ctx));
    tedw_to_mpz(Z, tedw_mulredc(P.Z, one, &ctx));
}

/* ================================================================== x8 */
#ifdef USE_AVX512F
#include <immintrin.h>
#include "monty.h"          /* the 52-bit carry helpers */

#define TEDW8_M52 0x000fffffffffffffULL
#define carryprop_tedw(lo, hi, mask) { __m512i carry_ = _mm512_srli_epi64(lo, 52); \
    hi = _mm512_add_epi64(hi, carry_); lo = _mm512_and_si512(mask, lo); }

typedef struct { __m512i d0, d1; } tedw_v;

#ifndef IFMA
static __m512d tedw8_dbias;
static __m512i tedw8_vbias1, tedw8_vbias2;

static TEDW_INLINE void tedw8_lohi(__m512i a, __m512i b, __m512i *l, __m512i *h)
{
    __m512d a_ld = _mm512_cvtepu64_pd(a);
    __m512d b_ld = _mm512_cvtepu64_pd(b);
    __m512d hd = _mm512_fmadd_round_pd(a_ld, b_ld, tedw8_dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    *h = _mm512_and_si512(_mm512_sub_epi64(_mm512_castpd_si512(hd), tedw8_vbias1),
        _mm512_set1_epi64(TEDW8_M52));
    hd = _mm512_sub_pd(_mm512_castsi512_pd(tedw8_vbias2), hd);
    *l = _mm512_and_si512(_mm512_castpd_si512(_mm512_fmadd_round_pd(a_ld, b_ld, hd,
        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC))), _mm512_set1_epi64(TEDW8_M52));
}
#endif

/* 104-bit Montgomery product (R = 2^104, positive inverse invN mod 2^104) */
static TEDW_INLINE tedw_v tedw8_mulredc(tedw_v x, tedw_v y, tedw_v N, tedw_v invN)
{
    const __m512i mask = _mm512_set1_epi64(TEDW8_M52);
    __m512i T0, T1, T2, T3, m0, m1, mN1, mN2, mN3;
#ifdef IFMA
    const __m512i z = _mm512_setzero_si512();
#define TL(a, b, cc) _mm512_madd52lo_epu64((a), (b), (cc))
#define TH(a, b, cc) _mm512_madd52hi_epu64((a), (b), (cc))
    T0 = TL(z, x.d0, y.d0);
    T1 = TL(TL(TH(z, x.d0, y.d0), x.d0, y.d1), x.d1, y.d0);
    T2 = TL(TH(TH(z, x.d0, y.d1), x.d1, y.d0), x.d1, y.d1);
    T3 = TH(z, x.d1, y.d1);
    carryprop_tedw(T1, T2, mask);
    carryprop_tedw(T2, T3, mask);
    m0 = TL(z, T0, invN.d0);
    m1 = _mm512_and_si512(TL(TL(TH(z, T0, invN.d0), T0, invN.d1), T1, invN.d0), mask);
    mN1 = TL(TL(TH(z, m0, N.d0), m0, N.d1), m1, N.d0);
    mN2 = TL(TH(TH(z, m0, N.d1), m1, N.d0), m1, N.d1);
    mN3 = TH(z, m1, N.d1);
#undef TL
#undef TH
#else
    __m512i l, h, l2, h2, l3, h3, l4, h4;
    tedw8_lohi(x.d0, y.d0, &T0, &T1);
    tedw8_lohi(x.d0, y.d1, &l2, &h2);
    tedw8_lohi(x.d1, y.d0, &l3, &h3);
    tedw8_lohi(x.d1, y.d1, &l4, &h4);
    T1 = _mm512_add_epi64(_mm512_add_epi64(T1, l2), l3);
    T2 = _mm512_add_epi64(_mm512_add_epi64(h2, h3), l4);
    T3 = h4;
    carryprop_tedw(T1, T2, mask);
    carryprop_tedw(T2, T3, mask);
    tedw8_lohi(T0, invN.d0, &m0, &m1);
    tedw8_lohi(T0, invN.d1, &l, &h);
    m1 = _mm512_add_epi64(m1, l);
    tedw8_lohi(T1, invN.d0, &l, &h);
    m1 = _mm512_and_si512(_mm512_add_epi64(m1, l), mask);
    __m512i mN0;
    tedw8_lohi(m0, N.d0, &mN0, &mN1);
    tedw8_lohi(m0, N.d1, &l2, &h2);
    tedw8_lohi(m1, N.d0, &l3, &h3);
    tedw8_lohi(m1, N.d1, &l4, &h4);
    mN1 = _mm512_add_epi64(_mm512_add_epi64(mN1, l2), l3);
    mN2 = _mm512_add_epi64(_mm512_add_epi64(h2, h3), l4);
    mN3 = h4;
    (void)mN0;
#endif
    carryprop_tedw(mN1, mN2, mask);
    carryprop_tedw(mN2, mN3, mask);
    // (T - mN) / 2^104, + N if T_hi < mN_hi (the low halves are equal)
    __mmask8 msk = _mm512_cmplt_epu64_mask(T3, mN3);
    msk |= _mm512_cmpeq_epu64_mask(T3, mN3) & _mm512_cmplt_epu64_mask(T2, mN2);
    __mmask8 bmsk;
    T2 = _mm512_subsetc_epi52(T2, mN2, &bmsk);
    T3 = _mm512_sbb_epi52(T3, bmsk, mN3, &bmsk);
    T2 = _mm512_mask_addsetc_epi52(T2, msk, T2, N.d0, &bmsk);
    T3 = _mm512_mask_adc_epi52(T3, msk, bmsk, N.d1, &bmsk);
    tedw_v r;
    r.d0 = _mm512_and_si512(T2, mask);
    r.d1 = _mm512_and_si512(T3, mask);
    return r;
}

static TEDW_INLINE tedw_v tedw8_addmod(tedw_v x, tedw_v y, tedw_v N)
{
    __mmask8 bmsk;
    __m512i x0 = _mm512_addsetc_epi52(x.d0, y.d0, &bmsk);
    __m512i x1 = _mm512_adc_epi52(x.d1, bmsk, y.d1, &bmsk);
    __mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(x1, N.d1);
    msk |= (_mm512_cmpeq_epu64_mask(x1, N.d1) & _mm512_cmpge_epu64_mask(x0, N.d0));
    x0 = _mm512_mask_subsetc_epi52(x0, msk, x0, N.d0, &bmsk);
    x1 = _mm512_mask_sbb_epi52(x1, msk, bmsk, N.d1, &bmsk);
    tedw_v r;
    r.d0 = _mm512_and_si512(x0, _mm512_set1_epi64(TEDW8_M52));
    r.d1 = _mm512_and_si512(x1, _mm512_set1_epi64(TEDW8_M52));
    return r;
}

static TEDW_INLINE tedw_v tedw8_submod(tedw_v x, tedw_v y, tedw_v N)
{
    __mmask8 msk = _mm512_cmplt_epu64_mask(x.d1, y.d1);
    msk |= _mm512_cmpeq_epu64_mask(x.d1, y.d1) & _mm512_cmplt_epu64_mask(x.d0, y.d0);
    __mmask8 bmsk;
    __m512i x0 = _mm512_subsetc_epi52(x.d0, y.d0, &bmsk);
    __m512i x1 = _mm512_sbb_epi52(x.d1, bmsk, y.d1, &bmsk);
    x0 = _mm512_mask_addsetc_epi52(x0, msk, x0, N.d0, &bmsk);
    x1 = _mm512_mask_adc_epi52(x1, msk, bmsk, N.d1, &bmsk);
    tedw_v r;
    r.d0 = _mm512_and_si512(x0, _mm512_set1_epi64(TEDW8_M52));
    r.d1 = _mm512_and_si512(x1, _mm512_set1_epi64(TEDW8_M52));
    return r;
}

static TEDW_INLINE tedw_v tedw8_zero(void)
{
    tedw_v z;
    z.d0 = z.d1 = _mm512_setzero_si512();
    return z;
}

static TEDW_INLINE __mmask8 tedw8_zmask(tedw_v a)
{
    return _mm512_cmpeq_epu64_mask(_mm512_or_si512(a.d0, a.d1), _mm512_setzero_si512());
}

static TEDW_INLINE tedw_v tedw8_blend(__mmask8 m, tedw_v a, tedw_v b)
{
    tedw_v r;
    r.d0 = _mm512_mask_mov_epi64(a.d0, m, b.d0);
    r.d1 = _mm512_mask_mov_epi64(a.d1, m, b.d1);
    return r;
}

typedef struct {
    tedw_v n, rho, one, Rsqr;
    tedw_v c6, c10, c16, c32, c48;
} tedw_ctx8;

#define EK_T tedw_v
#define EK_CTX tedw_ctx8
#define EK_FN(x) tedw8_##x
#define EK_MUL(a, b) tedw8_mulredc((a), (b), c->n, c->rho)
#define EK_SQR(a) tedw8_mulredc((a), (a), c->n, c->rho)
#define EK_ADD(a, b) tedw8_addmod((a), (b), c->n)
#define EK_SUB(a, b) tedw8_submod((a), (b), c->n)
#define EK_ZERO tedw8_zero()
#define EK_ACC_DECL __mmask8 ek_msk = 0xff;
#define EK_ACC(acc, t) { tedw_v tmp_ = EK_MUL(acc, t); \
    ek_msk &= (__mmask8)~tedw8_zmask(tmp_); acc = tedw8_blend(ek_msk, acc, tmp_); }
#include "edw_kern.h"
#undef EK_T
#undef EK_CTX
#undef EK_FN
#undef EK_MUL
#undef EK_SQR
#undef EK_ADD
#undef EK_SUB
#undef EK_ZERO
#undef EK_ACC_DECL
#undef EK_ACC

static void tedw8_init_consts(void)
{
#ifndef IFMA
    tedw8_dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    tedw8_vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);
    tedw8_vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);
#endif
}

/* per-lane modulus data */
typedef struct {
    ALIGNED_MEM uint64_t n[2][8], rho[2][8], one[2][8], Rsqr[2][8];
} tedw_lanes;

static tedw_v tedw8_load(uint64_t a[2][8])
{
    tedw_v r;
    r.d0 = _mm512_load_si512(a[0]);
    r.d1 = _mm512_load_si512(a[1]);
    return r;
}

static void tedw8_mpz_to_lane(uint64_t a[2][8], int lane, const mpz_t x)
{
    mpz_t t;
    mpz_init_set(t, x);
    a[0][lane] = (uint64_t)(tedw_from_mpz(t).lo & TEDW8_M52);
    mpz_tdiv_q_2exp(t, t, 52);
    a[1][lane] = (uint64_t)(tedw_from_mpz(t).lo & TEDW8_M52);
    mpz_clear(t);
}

static void tedw8_lane_to_mpz(mpz_t x, const uint64_t *d0, const uint64_t *d1, int lane)
{
    tedw_u128 v;
    v.lo = d0[lane] | (d1[lane] << 52);
    v.hi = d1[lane] >> 12;
    tedw_to_mpz(x, v);
}

/* lane setup for modulus n (odd, < 2^104) */
static void tedw8_set_lane(tedw_lanes *L, int lane, const mpz_t n)
{
    mpz_t t, R;
    mpz_init(t);
    mpz_init(R);
    tedw8_mpz_to_lane(L->n, lane, n);
    mpz_set_ui(R, 1);
    mpz_mul_2exp(R, R, 104);
    mpz_mod(t, R, n);
    tedw8_mpz_to_lane(L->one, lane, t);
    mpz_mul(t, t, t);
    mpz_mod(t, t, n);
    tedw8_mpz_to_lane(L->Rsqr, lane, t);
    mpz_invert(t, n, R);
    tedw8_mpz_to_lane(L->rho, lane, t);
    mpz_clear(t);
    mpz_clear(R);
}

static void tedw8_ctx_load(tedw_ctx8 *c, tedw_lanes *L)
{
    c->n = tedw8_load(L->n);
    c->rho = tedw8_load(L->rho);
    c->one = tedw8_load(L->one);
    c->Rsqr = tedw8_load(L->Rsqr);
    tedw_v c2 = tedw8_addmod(c->one, c->one, c->n);
    tedw_v c4 = tedw8_addmod(c2, c2, c->n);
    tedw_v c8 = tedw8_addmod(c4, c4, c->n);
    c->c16 = tedw8_addmod(c8, c8, c->n);
    c->c32 = tedw8_addmod(c->c16, c->c16, c->n);
    c->c48 = tedw8_addmod(c->c32, c->c16, c->n);
    c->c6 = tedw8_addmod(c4, c2, c->n);
    c->c10 = tedw8_addmod(c8, c2, c->n);
}

/* a signed constant per lane (|v| < 2^63) to Montgomery form */
static tedw_v tedw8_const(const int64_t *v, const tedw_ctx8 *c)
{
    ALIGNED_MEM uint64_t a[2][8];
    __mmask8 neg = 0;
    int i;
    for (i = 0; i < 8; i++) {
        uint64_t u = (v[i] < 0) ? (uint64_t)0 - (uint64_t)v[i] : (uint64_t)v[i];
        a[0][i] = u & TEDW8_M52;
        a[1][i] = u >> 52;
        if (v[i] < 0) neg |= (__mmask8)(1 << i);
    }
    tedw_v r = tedw8_mulredc(tedw8_load(a), c->Rsqr, c->n, c->rho);
    return tedw8_blend(neg, r, tedw8_submod(tedw8_zero(), r, c->n));
}

static void tedw8_table_curve(const int *k, tedw8_curve *cv, const tedw_ctx8 *c)
{
    tedw_v w[6];
    int64_t v[8];
    int i, j;
    for (j = 0; j < 6; j++) {
        for (i = 0; i < 8; i++)
            v[i] = edw_tiny_curves[k[i]][j];
        w[j] = tedw8_const(v, c);
    }
    tedw8_curve_from(w, cv, c);
}

static void tedw8_blend_curve(__mmask8 m, tedw8_curve *cv, const tedw8_curve *src)
{
    cv->P.X = tedw8_blend(m, cv->P.X, src->P.X);
    cv->P.Y = tedw8_blend(m, cv->P.Y, src->P.Y);
    cv->P.Z = tedw8_blend(m, cv->P.Z, src->P.Z);
    cv->P.T = tedw8_blend(m, cv->P.T, src->P.T);
    cv->dn2 = tedw8_blend(m, cv->dn2, src->dn2);
    cv->dd = tedw8_blend(m, cv->dd, src->dd);
    cv->dd2 = tedw8_blend(m, cv->dd2, src->dd2);
    cv->a24n = tedw8_blend(m, cv->a24n, src->a24n);
    cv->a24d = tedw8_blend(m, cv->a24d, src->a24d);
}

static void tedw8_blend_w(__mmask8 m, tedw8_wpt *R, const tedw8_wpt *S)
{
    R->X = tedw8_blend(m, R->X, S->X);
    R->Y = tedw8_blend(m, R->Y, S->Y);
    R->Z = tedw8_blend(m, R->Z, S->Z);
}

/* the constant [3]Q = (1205/49, -40664/343) = (413315 : -1992536 : 16807) */
static void tedw8_q3(tedw8_wpt *R, const tedw_ctx8 *c)
{
    int64_t v[3][8];
    int i;
    for (i = 0; i < 8; i++) {
        v[0][i] = 413315;
        v[1][i] = -1992536;
        v[2][i] = 16807;
    }
    R->X = tedw8_const(v[0], c);
    R->Y = tedw8_const(v[1], c);
    R->Z = tedw8_const(v[2], c);
}

/* stage 1 / stage 2 checks of lane l (mod the lane's n) */
static int tedw8_lane_check(const tedw_v *x, tedw_lanes *L, int lane, mpz_t f)
{
    ALIGNED_MEM uint64_t d0[8], d1[8];
    mpz_t g, n;
    int r = 0;
    _mm512_store_si512(d0, x->d0);
    _mm512_store_si512(d1, x->d1);
    mpz_init(g);
    mpz_init(n);
    tedw8_lane_to_mpz(g, d0, d1, lane);
    tedw8_lane_to_mpz(n, L->n[0], L->n[1], lane);
    mpz_gcd(g, g, n);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
        mpz_set(f, g);
        r = 1;
    }
    mpz_clear(g);
    mpz_clear(n);
    return r;
}

static void tedw8_store_factor(uint64_t *f64, uint32_t i, const mpz_t f)
{
    tedw_u128 v = tedw_from_mpz(f);
    f64[2 * i] = v.lo;
    f64[2 * i + 1] = v.hi;
}

/* the x8 list: inputs ix[0..num-1] (indices into q64 / f64), curves kstart
   .. kstart + curves - 1 of each input's sequence (table curves, then GKL
   m = 3, 4, ...); returns the curves run */
static uint64_t tedw8_list(const uint64_t *q64, uint64_t *f64, const uint32_t *ix, uint32_t num,
    uint32_t B1, uint32_t B2, int kstart, uint32_t curves)
{
    tedw_lanes L;
    tedw_ctx8 ctx;
    const tedw_ctx8 *c = &ctx;
    edw_chain_t ch;
    const uint8_t *barray;
    int numb;
    int kcur[8], k8[8], lane;
    uint32_t idx[8], j = 0;
    uint8_t msk = 0;
    uint64_t total = 0;
    const int ntab = EDW_TINY_NCURVES;
    const int kend = kstart + (int)curves;
    tedw8_wpt R;
    mpz_t n, f;

    if (num == 0 || !edw_chain_for(B1, &ch) || !tedw_s2table(B1, &barray, &numb))
        return 0;
    mpz_init(n);
    mpz_init(f);
    memset(&R, 0, sizeof(R));
    // every lane gets a valid modulus, so idle lanes compute harmlessly
    for (lane = 0; lane < 8; lane++) {
        mpz_set_ui(n, (unsigned long)(q64[2 * ix[0] + 1] >> 26));
        mpz_mul_2exp(n, n, 26);
        mpz_add_ui(n, n, (unsigned long)(q64[2 * ix[0] + 1] & 0x3ffffff));
        mpz_mul_2exp(n, n, 26);
        mpz_add_ui(n, n, (unsigned long)(q64[2 * ix[0]] >> 26));
        mpz_mul_2exp(n, n, 26);
        mpz_add_ui(n, n, (unsigned long)(q64[2 * ix[0]] & 0x3ffffff));
        tedw8_set_lane(&L, lane, n);
        kcur[lane] = kstart;
        idx[lane] = ix[0];
    }

    while (1) {
        uint8_t lmsk = 0;
        for (lane = 0; lane < 8; lane++) {
            if (!(msk & (1 << lane)) && j < num) {
                uint32_t i = ix[j++];
                mpz_set_ui(n, (unsigned long)(q64[2 * i + 1] >> 26));
                mpz_mul_2exp(n, n, 26);
                mpz_add_ui(n, n, (unsigned long)(q64[2 * i + 1] & 0x3ffffff));
                mpz_mul_2exp(n, n, 26);
                mpz_add_ui(n, n, (unsigned long)(q64[2 * i] >> 26));
                mpz_mul_2exp(n, n, 26);
                mpz_add_ui(n, n, (unsigned long)(q64[2 * i] & 0x3ffffff));
                tedw8_set_lane(&L, lane, n);
                kcur[lane] = kstart;
                idx[lane] = i;
                lmsk |= (uint8_t)(1 << lane);
            }
        }
        msk |= lmsk;
        if (msk == 0)
            break;
        if (lmsk)
            tedw8_ctx_load(&ctx, &L);

        uint8_t tabm = 0, newm = 0;
        for (lane = 0; lane < 8; lane++) {
            k8[lane] = (kcur[lane] < ntab) ? kcur[lane] : 0;
            if (kcur[lane] < ntab) tabm |= (uint8_t)(1 << lane);
            if (kcur[lane] == ntab) newm |= (uint8_t)(1 << lane);
        }
        tedw8_curve cv, cvg;
        if (tabm)
            tedw8_table_curve(k8, &cv, c);
        if ((uint8_t)~tabm) {
            __mmask8 gm = (__mmask8)(uint8_t)~tabm;
            if (newm) {
                tedw8_wpt Q3;
                tedw8_q3(&Q3, c);
                tedw8_blend_w(newm, &R, &Q3);
            }
            tedw8_gkl_curve(&R, &cvg, c);
            tedw8_wpt R2 = R;
            tedw8_w_addq(&R2, c);
            tedw8_blend_w(gm, &R, &R2);
            if (!tabm)
                cv = cvg;
            else
                tedw8_blend_curve(gm, &cv, &cvg);
        }
        for (lane = 0; lane < 8; lane++)
            kcur[lane]++;

        tedw8_pt P;
        tedw8_stage1(&cv, &ch, &P, c);
        for (lane = 0; lane < 8; lane++) {
            if (msk & (1 << lane)) {
                total++;
                if (tedw8_lane_check(&P.Z, &L, lane, f)) {
                    tedw8_store_factor(f64, idx[lane], f);
                    msk &= (uint8_t)~(1 << lane);
                }
            }
        }
        if (B2 > B1 && msk) {
            tedw_v acc = tedw8_stage2(&P, &cv, B1, barray, numb, c);
            for (lane = 0; lane < 8; lane++) {
                if ((msk & (1 << lane)) && tedw8_lane_check(&acc, &L, lane, f)) {
                    tedw8_store_factor(f64, idx[lane], f);
                    msk &= (uint8_t)~(1 << lane);
                }
            }
        }
        for (lane = 0; lane < 8; lane++)
            if (kcur[lane] >= kend)
                msk &= (uint8_t)~(1 << lane);
    }
    mpz_clear(n);
    mpz_clear(f);
    return total;
}

/* one n in all lanes, lane i running curve base + i of the sequence: table
   curves 0..4, then GKL m = idx - 2 (3, 4, 5, ...).  Lanes start their
   [m]Q lazily, from a running [k]Q (k = 3, 4, ...), and then step it by
   [8]Q each round (a doubling instead for the lane at m = 8). */
typedef struct {
    tedw_lanes L;
    tedw_ctx8 ctx;
    tedw8_wpt R, S, cur;
    int mcur[8];            // the lane's current GKL m, 0 if not started
    int kcur;               // cur = [kcur]Q, 0 if not started
    int have_s;
    int base;
} tedw8_one_t;

static void tedw8_one_init(tedw8_one_t *st, const mpz_t n)
{
    int lane;
    tedw8_set_lane(&st->L, 0, n);           // one mpz setup, copied to all lanes
    for (lane = 1; lane < 8; lane++) {
        st->L.n[0][lane] = st->L.n[0][0];
        st->L.n[1][lane] = st->L.n[1][0];
        st->L.rho[0][lane] = st->L.rho[0][0];
        st->L.rho[1][lane] = st->L.rho[1][0];
        st->L.one[0][lane] = st->L.one[0][0];
        st->L.one[1][lane] = st->L.one[1][0];
        st->L.Rsqr[0][lane] = st->L.Rsqr[0][0];
        st->L.Rsqr[1][lane] = st->L.Rsqr[1][0];
    }
    tedw8_ctx_load(&st->ctx, &st->L);
    for (lane = 0; lane < 8; lane++)
        st->mcur[lane] = 0;
    st->kcur = 0;
    st->have_s = 0;
    st->base = 0;
}

/* lanes in gm get their GKL curve of this round in cvg, and step */
static void tedw8_one_gkl(tedw8_one_t *st, __mmask8 gm, tedw8_curve *cvg)
{
    const tedw_ctx8 *c = &st->ctx;
    int lane;
    // start lanes that need their first [m]Q, m = idx - 2, in increasing m
    for (lane = 0; lane < 8; lane++) {
        int m;
        if (!(gm & (1 << lane)) || st->mcur[lane])
            continue;
        m = st->base + lane - 2;
        if (st->kcur == 0) {
            tedw8_q3(&st->cur, c);
            st->kcur = 3;
        }
        while (st->kcur < m) {
            tedw8_w_addq(&st->cur, c);
            st->kcur++;
        }
        tedw8_blend_w((__mmask8)(1 << lane), &st->R, &st->cur);
        st->mcur[lane] = m;
    }
    tedw8_gkl_curve(&st->R, cvg, c);
    // step by [8]Q
    if (!st->have_s) {
        int64_t v[3][8];
        for (lane = 0; lane < 8; lane++) {
            v[0][lane] = 5;
            v[1][lane] = 8;
            v[2][lane] = 1;
        }
        st->S.X = tedw8_const(v[0], c);
        st->S.Y = tedw8_const(v[1], c);
        st->S.Z = tedw8_const(v[2], c);
        tedw8_w_dbl(&st->S, c);
        tedw8_w_dbl(&st->S, c);
        tedw8_w_dbl(&st->S, c);
        st->have_s = 1;
    }
    __mmask8 dm = 0;
    for (lane = 0; lane < 8; lane++)
        if ((gm & (1 << lane)) && st->mcur[lane] == 8)
            dm |= (__mmask8)(1 << lane);
    tedw8_wpt R2 = st->R;
    tedw8_w_addp(&R2, &st->S, c);
    if (dm) {
        tedw8_wpt R3 = st->R;
        tedw8_w_dbl(&R3, c);
        tedw8_blend_w(dm, &R2, &R3);
    }
    tedw8_blend_w(gm, &st->R, &R2);
    for (lane = 0; lane < 8; lane++)
        if (gm & (1 << lane))
            st->mcur[lane] += 8;
}

/* rounds of 8 curves; returns 1 with the factor in f */
static int tedw8_one_rounds(tedw8_one_t *st, const mpz_t n, mpz_t f, uint32_t B1, int rounds,
    int *done_rounds)
{
    const tedw_ctx8 *c = &st->ctx;
    edw_chain_t ch;
    const uint8_t *barray;
    int numb, r, lane;
    const int ntab = EDW_TINY_NCURVES;
    (void)n;
    *done_rounds = 0;
    if (!edw_chain_for(B1, &ch) || !tedw_s2table(B1, &barray, &numb))
        return 0;
    for (r = 0; r < rounds; r++) {
        int k8[8];
        uint8_t tabm = 0;
        for (lane = 0; lane < 8; lane++) {
            int id = st->base + lane;
            k8[lane] = (id < ntab) ? id : 0;
            if (id < ntab) tabm |= (uint8_t)(1 << lane);
        }
        tedw8_curve cv, cvg;
        if (tabm)
            tedw8_table_curve(k8, &cv, c);
        if ((uint8_t)~tabm) {
            __mmask8 gm = (__mmask8)(uint8_t)~tabm;
            tedw8_one_gkl(st, gm, &cvg);
            if (!tabm)
                cv = cvg;
            else
                tedw8_blend_curve(gm, &cv, &cvg);
        }
        st->base += 8;
        (*done_rounds)++;
        tedw8_pt P;
        tedw8_stage1(&cv, &ch, &P, c);
        for (lane = 0; lane < 8; lane++)
            if (tedw8_lane_check(&P.Z, &st->L, lane, f))
                return 1;
        tedw_v acc = tedw8_stage2(&P, &cv, B1, barray, numb, c);
        for (lane = 0; lane < 8; lane++)
            if (tedw8_lane_check(&acc, &st->L, lane, f))
                return 1;
    }
    return 0;
}
#endif /* USE_AVX512F */

int getfactor_tecm_edw_x8(mpz_t n, mpz_t f, int target_bits, uint64_t *pran)
{
#ifdef USE_AVX512F
    tedw8_one_t *st;
    int dr, found = 0;
    uint32_t B1;
    int rounds;
    mpz_set_ui(f, 1);
    if (mpz_even_p(n)) {
        mpz_set_ui(f, 2);
        return 1;
    }
    if (mpz_sizeinbase(n, 2) > 104 || mpz_cmp_ui(n, 3) < 0)
        return getfactor_tecm_edw(n, f, target_bits, pran);
    tedw8_init_consts();
    st = (tedw8_one_t *)malloc(sizeof(tedw8_one_t) + 64);
    tedw8_one_t *s = (tedw8_one_t *)(((uintptr_t)st + 63) & ~(uintptr_t)63);
    tedw8_one_init(s, n);
    if (target_bits == 0) {
        // fast attempts (as tecm_dispatch_x8): one round each
        if (tedw8_one_rounds(s, n, f, 47, 1, &dr) || tedw8_one_rounds(s, n, f, 70, 1, &dr) ||
            tedw8_one_rounds(s, n, f, 125, 1, &dr)) {
            free(st);
            return 1;
        }
    }
    if (target_bits <= 20) { B1 = 27; rounds = 4; }
    else if (target_bits <= 22) { B1 = 47; rounds = 4; }
    else if (target_bits <= 24) { B1 = 70; rounds = 4; }
    else if (target_bits <= 26) { B1 = 85; rounds = 4; }
    else if (target_bits <= 29) { B1 = 125; rounds = 4; }
    else if (target_bits <= 31) { B1 = 165; rounds = 6; }
    else { B1 = 205; rounds = 6; }
    found = tedw8_one_rounds(s, n, f, B1, rounds, &dr);
    free(st);
    return found;
#else
    return getfactor_tecm_edw(n, f, target_bits, pran);
#endif
}

void getfactor_tecm_edw_x8_list(uint64_t *q64, uint64_t *f64, int target_bits,
    uint32_t num_in, uint64_t *pran)
{
    uint32_t i;
    for (i = 0; i < num_in; i++) {
        f64[2 * i] = 1;
        f64[2 * i + 1] = 0;
    }
#ifdef USE_AVX512F
    uint32_t *ix = (uint32_t *)malloc((num_in + 1) * sizeof(uint32_t));
    uint32_t num = num_in, B1, curves;
    int kstart = 0;
    tedw8_init_consts();
    for (i = 0; i < num_in; i++)
        ix[i] = i;
    (void)pran;
    if (target_bits == 0) {
        // fast attempts (as tecm_dispatch_x8_list), one curve each, then
        // only the inputs not yet factored go on
        static const uint32_t pre[3] = { 47, 70, 125 };
        int p;
        for (p = 0; p < 3; p++) {
            tedw8_list(q64, f64, ix, num, pre[p], 25 * pre[p], kstart, 1);
            kstart++;
            uint32_t k = 0;
            for (i = 0; i < num; i++)
                if (f64[2 * ix[i]] == 1 && f64[2 * ix[i] + 1] == 0)
                    ix[k++] = ix[i];
            num = k;
        }
    }
    if (target_bits <= 20) { B1 = 27; curves = 32; }
    else if (target_bits <= 22) { B1 = 47; curves = 32; }
    else if (target_bits <= 24) { B1 = 70; curves = 32; }
    else if (target_bits <= 26) { B1 = 85; curves = 32; }
    else if (target_bits <= 29) { B1 = 125; curves = 32; }
    else if (target_bits <= 31) { B1 = 165; curves = 42; }
    else { B1 = 205; curves = 42; }
    tedw8_list(q64, f64, ix, num, B1, 25 * B1, kstart, curves);
    free(ix);
#else
    mpz_t n, f;
    mpz_init(n);
    mpz_init(f);
    for (i = 0; i < num_in; i++) {
        mpz_set_ui(n, (unsigned long)(q64[2 * i + 1] >> 26));
        mpz_mul_2exp(n, n, 26);
        mpz_add_ui(n, n, (unsigned long)(q64[2 * i + 1] & 0x3ffffff));
        mpz_mul_2exp(n, n, 26);
        mpz_add_ui(n, n, (unsigned long)(q64[2 * i] >> 26));
        mpz_mul_2exp(n, n, 26);
        mpz_add_ui(n, n, (unsigned long)(q64[2 * i] & 0x3ffffff));
        if (getfactor_tecm_edw(n, f, target_bits, pran)) {
            tedw_u128 v = tedw_from_mpz(f);
            f64[2 * i] = v.lo;
            f64[2 * i + 1] = v.hi;
        }
    }
    mpz_clear(n);
    mpz_clear(f);
#endif
}

uint64_t tecm_edw_x8_run(uint64_t *q64, uint64_t *f64, uint32_t num, uint32_t B1, uint32_t curves)
{
#ifdef USE_AVX512F
    uint32_t i;
    uint64_t total;
    uint32_t *ix = (uint32_t *)malloc((num + 1) * sizeof(uint32_t));
    tedw8_init_consts();
    for (i = 0; i < num; i++) {
        ix[i] = i;
        f64[2 * i] = 1;
        f64[2 * i + 1] = 0;
    }
    total = tedw8_list(q64, f64, ix, num, B1, 25 * B1, 0, curves);
    free(ix);
    return total;
#else
    (void)q64; (void)f64; (void)num; (void)B1; (void)curves;
    return 0;
#endif
}

int tecm_edw_x8_one(mpz_t n, mpz_t f, uint32_t B1, uint32_t curves, int *rounds)
{
#ifdef USE_AVX512F
    int found;
    tedw8_one_t *st = (tedw8_one_t *)malloc(sizeof(tedw8_one_t) + 64);
    tedw8_one_t *s = (tedw8_one_t *)(((uintptr_t)st + 63) & ~(uintptr_t)63);
    tedw8_init_consts();
    mpz_set_ui(f, 1);
    tedw8_one_init(s, n);
    found = tedw8_one_rounds(s, n, f, B1, (int)((curves + 7) / 8), rounds);
    free(st);
    return found;
#else
    (void)n; (void)f; (void)B1; (void)curves;
    *rounds = 0;
    return 0;
#endif
}
