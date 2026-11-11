/* ===========================================================================
   limb2.h -- double-limb arithmetic: 128-bit scalar, and the 104-bit
              (2 x 52-bit) AVX-512 vector layer.

   Companion to limb1.h (single-limb 64-bit scalar + 52-bit vector), which it
   includes, and which in turn pulls the compiler / ISA layer and intrinsic
   polyfills from mp_platform.h.  Nothing outside arith/limb is needed.

   Conventions
     - a two-limb value is a little-endian pair: x[0] = low 64 bits, x[1] = high
     - Montgomery routines use R = 2^128 and the NEGATIVE inverse
       rho = -n^-1 mod 2^64 (monty128_t.rho); the modulus must be odd

   Sections
     1. plain two-limb add / sub / compare / shift / multiply / square
        (inline) and divide (limb2.c)
     2. 128-bit Montgomery arithmetic (limb2.c: the CIOS cores are large, so
        they stay out of line)
     3. 104 / 128-bit zero-safe bit scans
     4. 128-bit gcd
     5. 104-bit vector Montgomery arithmetic, 8 lanes (AVX-512, optional IFMA)

   Opt-outs: define LIMB2_NO_VECTOR (or LIMB1_SCALAR_ONLY) before including
   to take only the scalar layers.
   =========================================================================== */

#ifndef LIMB2_H
#define LIMB2_H

#include "limb1.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================
   1. Plain two-limb arithmetic
   =========================================================================== */

/* r = a + b, returns the carry out.  r may alias a or b. */
MP_FORCE_INLINE mp_carry_t add128(const uint64_t *a, const uint64_t *b, uint64_t *r)
{
    mp_carry_t c = mp_addcarry_u64(0, a[0], b[0], &r[0]);
    return mp_addcarry_u64(c, a[1], b[1], &r[1]);
}

/* r = a - b, returns the borrow out.  r may alias a or b. */
MP_FORCE_INLINE mp_carry_t sub128(const uint64_t *a, const uint64_t *b, uint64_t *r)
{
    mp_carry_t c = mp_subborrow_u64(0, a[0], b[0], &r[0]);
    return mp_subborrow_u64(c, a[1], b[1], &r[1]);
}

/* -1, 0, 1 for a < b, a == b, a > b */
MP_FORCE_INLINE int cmp128(const uint64_t *a, const uint64_t *b)
{
    if (a[1] != b[1]) return (a[1] < b[1]) ? -1 : 1;
    if (a[0] != b[0]) return (a[0] < b[0]) ? -1 : 1;
    return 0;
}

/* r = a << s and r = a >> s, for 0 <= s < 128.  r may alias a. */
MP_FORCE_INLINE void shl128(const uint64_t *a, int s, uint64_t *r)
{
    uint64_t lo = a[0], hi = a[1];
    if (s >= 64)     { r[1] = lo << (s - 64); r[0] = 0; }
    else if (s > 0)  { r[1] = (hi << s) | (lo >> (64 - s)); r[0] = lo << s; }
    else             { r[1] = hi; r[0] = lo; }
}

MP_FORCE_INLINE void shr128(const uint64_t *a, int s, uint64_t *r)
{
    uint64_t lo = a[0], hi = a[1];
    if (s >= 64)     { r[0] = hi >> (s - 64); r[1] = 0; }
    else if (s > 0)  { r[0] = (lo >> s) | (hi << (64 - s)); r[1] = hi >> s; }
    else             { r[0] = lo; r[1] = hi; }
}

/* r[0..2] = a[0..1] * b : 128 x 64 -> 192 bits.  r must not alias a. */
MP_FORCE_INLINE void mul128x64(const uint64_t *a, uint64_t b, uint64_t *r)
{
    uint64_t h0, h1, l1;
    mp_carry_t c;
    r[0] = _umul128(a[0], b, &h0);
    l1 = _umul128(a[1], b, &h1);
    c = mp_addcarry_u64(0, l1, h0, &r[1]);
    r[2] = h1 + c;
}

/* r[0..3] = a[0..1] * b[0..1] : 128 x 128 -> 256 bits.  r must not alias. */
MP_FORCE_INLINE void mul128(const uint64_t *a, const uint64_t *b, uint64_t *r)
{
    uint64_t h00, l01, h01, l10, h10, l11, h11, t;
    mp_carry_t c;
    r[0] = _umul128(a[0], b[0], &h00);
    l01 = _umul128(a[0], b[1], &h01);
    l10 = _umul128(a[1], b[0], &h10);
    l11 = _umul128(a[1], b[1], &h11);
    c = mp_addcarry_u64(0, h00, l01, &t);
    c = mp_addcarry_u64(c, h01, l11, &r[2]);
    r[3] = h11 + c;
    c = mp_addcarry_u64(0, t, l10, &r[1]);
    c = mp_addcarry_u64(c, r[2], h10, &r[2]);
    r[3] += c;
}

/* r[0..3] = a[0..1]^2 : 128 -> 256 bits.  r must not alias a. */
MP_FORCE_INLINE void sqr128(const uint64_t *a, uint64_t *r)
{
    uint64_t h00, l01, h01, l11, h11;
    mp_carry_t c;
    r[0] = _umul128(a[0], a[0], &h00);
    l01 = _umul128(a[0], a[1], &h01);
    l11 = _umul128(a[1], a[1], &h11);
    /* add the cross product twice: (h01:l01) << 65 */
    r[3] = h11 + (h01 >> 63);
    h01 = (h01 << 1) | (l01 >> 63);
    l01 <<= 1;
    c = mp_addcarry_u64(0, h00, l01, &r[1]);
    c = mp_addcarry_u64(c, l11, h01, &r[2]);
    r[3] += c;
}

/* q = n / d and r = n % d (128-bit operands); q or r may be NULL.
   Returns 0, or -1 on division by zero.  Knuth D with 64-bit digits. */
int uint128_div(const uint64_t n[2], const uint64_t d[2],
    uint64_t q[2], uint64_t r[2]);

/* c = a mod b (128-bit), b != 0 */
void mod_128(const uint64_t *a, const uint64_t *b, uint64_t *c);


/* ===========================================================================
   2. 128-bit Montgomery arithmetic
   =========================================================================== */
typedef struct
{
	uint64_t r[2];
	uint64_t n[2];
	uint64_t np[2];
	uint64_t nhat[2];
	uint64_t rhat[2];
	uint64_t rmask[2];
	uint64_t one[2];
	uint64_t mtmp1[2];
	uint64_t mtmp2[2];
	uint64_t mtmp3[2];
	uint64_t mtmp4[2];
	uint64_t rho;
} monty128_t;

#ifdef HAS_UINT128
#define USE_PERIG_128BIT
#endif

/* setup: n odd.  monty128_init fills in n, rho and one (= R mod n);
   to_monty128 converts x in place (x need not be reduced). */
void monty128_init(monty128_t * mdata, uint64_t * n);
void to_monty128(monty128_t *mdata, uint64_t * x);

/* w = u * v * R^-1 mod n and w = u^2 * R^-1 mod n.  The ...n variants take
   the modulus and rho directly instead of a monty128_t. */
void mulmod128(uint64_t * u, uint64_t * v, uint64_t * w, monty128_t *mdata);
void mulmod128n(uint64_t* u, uint64_t* v, uint64_t* w, uint64_t* n, uint64_t nhat);
void sqrmod128(uint64_t * u, uint64_t * w, monty128_t *mdata);
void sqrmod128n(uint64_t* u, uint64_t* w, uint64_t *n, uint64_t nhat);
/* the unreduced core: w[0..2] must be ZERO on entry; on return it holds
   u*v*R^-1 mod n, possibly plus n (the caller does the final subtraction) */
void ciosFullMul128x(uint64_t *u, uint64_t *v, uint64_t rho, uint64_t *n, uint64_t *w);

/* w = u +/- v mod n (inputs reduced); a = 2a mod n; a = a - n if a >= n */
void addmod128(uint64_t * u, uint64_t * v, uint64_t * w, uint64_t * n);
void submod128(uint64_t * u, uint64_t * v, uint64_t * w, uint64_t * n);
void dblmod128(uint64_t* a, uint64_t* n);
void chkmod128(uint64_t* a, uint64_t* n);


/* ===========================================================================
   3. 104/128-bit bit scans.  COUNT semantics, and SAFE at zero (they return
   the width) -- unlike the _lead_zcnt64 / _trail_zcnt64 pair in
   mp_platform.h, which are undefined at zero on the gcc/clang builtin paths.
   The two families are NOT interchangeable; see the note in limb2.c.
   =========================================================================== */
uint64_t my_clz128(uint64_t n_lo, uint64_t n_hi);
uint64_t my_clz104(uint64_t n_lo, uint64_t n_hi);
uint64_t my_ctz104(uint64_t n_lo, uint64_t n_hi);
uint64_t my_ctz128(uint64_t nlo, uint64_t nhi);

/* legacy single-word names (these were functions in monty.c) */
#define my_clz64 _lead_full_zcnt
#define my_ctz64 _trail_full_zcnt
#define my_clz52 _lead_full_zcnt52
#define my_ctz52 _trail_full_zcnt52


/* ===========================================================================
   4. 128-bit gcd.
   gcd128 (Euclid) works on every compiler: it divides with uint128_div
   rather than the native 128-bit `%`, which clang-cl and icx-cl accept but
   cannot link (no __umodti3 in the MSVC runtime) and MSVC does not have.
   bin_gcd128 (binary, u odd, clobbers u and v) still needs a native 128-bit
   integer type for its adds and compares, so it is gcc / clang / icc only.
   =========================================================================== */
static UNUSED_FUNC void gcd128(uint64_t* u, uint64_t* v, uint64_t* w);
static UNUSED_FUNC void gcd128(uint64_t* u, uint64_t* v, uint64_t* w)
{
	uint64_t a[2], b[2], c[2];
	a[0] = u[0]; a[1] = u[1];
	b[0] = v[0]; b[1] = v[1];
	while ((b[0] | b[1]) != 0)
	{
		(void)uint128_div(a, b, (uint64_t *)0, c);
		a[0] = b[0]; a[1] = b[1];
		b[0] = c[0]; b[1] = c[1];
	}
	w[1] = a[1];
	w[0] = a[0];
	return;
}

#if defined(__SIZEOF_INT128__) && (__SIZEOF_INT128__ == 16)

static UNUSED_FUNC void bin_gcd128(uint64_t *u, uint64_t *v, uint64_t *w);
static UNUSED_FUNC void bin_gcd128(uint64_t *u, uint64_t *v, uint64_t *w)
{
	//w = gcd(u, v);
	if ((u[1] == 0) && (u[0] == 0))
	{
		w[1] = v[1];
		w[0] = v[0];
		return;
	}
	if ((v[1] != 0) || (v[0] != 0)) {
		int j = (int)my_ctz128(v[0], v[1]);
		//v = (uint64_t)(v >> j);
		//if (j > 64) printf("j overflow\n");
		// v >>= j for any j in [0, 127].  (The open-coded two-limb shift that
		// used to be here was wrong for j == 0, a shift by 64, and for j >= 64.)
		shr128(v, j, v);
		while (1) {

			//uint64_t tmp = u;
			//uint64_t sub1 = (uint64_t)(v - tmp);
			//uint64_t sub2 = (uint64_t)(tmp - v);
			//if (tmp == v)
			//	break;
			//u = (tmp >= v) ? v : tmp;
			//v = (tmp >= v) ? sub2 : sub1;


			uint128_t t = (uint128_t)u[1] << 64 | (uint128_t)u[0];
			uint128_t v128 = ((uint128_t)v[1] << 64 | (uint128_t)v[0]);
			uint128_t s1 = v128 - t;
			uint128_t s2 = t - v128;

			if (t == v128)
				break;


			u[0] = (t >= v128) ? (uint64_t)v128 : (uint64_t)t;
			u[1] = (t >= v128) ? (uint64_t)(v128 >> 64) : (uint64_t)(t >> 64);

			v[0] = (t >= v128) ? (uint64_t)s2 : (uint64_t)s1;
			v[1] = (t >= v128) ? (uint64_t)(s2 >> 64) : (uint64_t)(s1 >> 64);

			// For the line below, the standard way to write this algorithm
			// would have been to use _trail_zcnt64(v)  (instead of
			// _trail_zcnt64(sub1)).  However, as pointed out by
			// https://gmplib.org/manual/Binary-GCD, "in twos complement the
			// number of low zero bits on u-v is the same as v-u, so counting or
			// testing can begin on u-v without waiting for abs(u-v) to be
			// determined."  Hence we are able to use sub1 for the argument.
			// By removing the dependency on abs(u-v), the CPU can execute
			// _trail_zcnt64() at the same time as abs(u-v).
			j = (int)my_ctz128((uint64_t)s1, (uint64_t)(s1 >> 64));
			//v = (uint64_t)(v >> j);
			//if (j > 64) printf("j overflow\n");
			shr128(v, j, v);
		}
	}
	w[1] = u[1];
	w[0] = u[0];
}
#endif


/* ===========================================================================
   5. 104-bit vector Montgomery arithmetic (8 lanes of 2 x 52 bits, AVX-512)
   Uses the 52-bit primitives from limb1.h (vec_u104_t, carryprop, the
   _mm512_*_epi52 carry ops, VEC_MUL_ACCUM_LOHI_PD).
   =========================================================================== */
#if MP_USE_AVX512F && !defined(LIMB2_NO_VECTOR) && !defined(LIMB1_SCALAR_ONLY)

// can we make these static inline in the header?  does it help?
//void mask_mulredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//    __m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0, __m512i vrho);
//void sqrredc104_vec(__m512i* c1, __m512i* c0,
//    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void mask_sqrredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void mask_sqrredc104_vec_pos(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//    __m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void mask_sqrredc104_exact_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho);
//void addmod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
//	__m512i b1, __m512i b0, __m512i n1, __m512i n0);
//void mask_addmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
//	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0);
//void mask_dblmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0);
//void mask_redsub104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0);
//void redsub104_x8(__m512i* c1, __m512i* c0,
//	__m512i a1, __m512i a0, __m512i n1, __m512i n0);
//void submod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
//	__m512i b1, __m512i b0, __m512i n1, __m512i n0);

void init_monty104(void);

#if defined(__GNUC__) && !MP_USE_IFMA
// the non-IFMA paths below intentionally declare local copies of the
// dbias/vbias constants used by the mul52 macros, shadowing the statics
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

__inline static void mask_mulredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, t3, C1, C2;

#if !MP_USE_IFMA
	UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	UNUSED_VAR int biascount = 0;
	UNUSED_VAR __m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	UNUSED_VAR __mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = t3 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b0, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b0, a0, a1, t0, t1, C1, t2);
	//t1 = _mm512_add_epi64(t1, C1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b1, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b1, a0, a1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);
	C1 = C2 = zero;

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);

	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
__inline static void mulredc104_vec(__m512i* c1, __m512i* c0, 
	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, t3, C1, C2;

#if !MP_USE_IFMA
	UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	UNUSED_VAR int biascount = 0;
	UNUSED_VAR __m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	UNUSED_VAR __mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = t3 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b0, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b0, a0, a1, t0, t1, C1, t2);
	//t1 = _mm512_add_epi64(t1, C1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b1, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(a1, b1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(b1, a0, a1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);
	C1 = C2 = zero;

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);
	//VEC_MUL2_ACCUM_LOHI_PD(m, n0, n1, t0, C1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);

	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);

	return;
}
__inline static void sqrredc104_vec(__m512i* c1, __m512i* c0,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, t3, C1, C2, sqr_lo, sqr_hi;

#if !MP_USE_IFMA
	UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	UNUSED_VAR int biascount = 0;
	UNUSED_VAR __m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	UNUSED_VAR __mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = t3 = C1 = C2 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	t1 = sqr_lo;
	t2 = sqr_hi;
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	// t0 = _mm512_add_epi64(t1, one);
	// t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);

	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		__mmask8 bmsk;
		bmsk = _mm512_cmpgt_epi64_mask(n0, t1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, _mm512_mask_set1_epi64(zero, bmsk, 1));
	}
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);

	return;
}
__inline static void mask_sqrredc104_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#if !MP_USE_IFMA
	UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	UNUSED_VAR int biascount = 0;
	UNUSED_VAR __m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	UNUSED_VAR __mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	t1 = _mm512_add_epi64(t1, sqr_lo);
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(sqr_hi, C2);
	// we throw t0 away after this so first propagate its carry.
	// t0 = _mm512_add_epi64(t1, one);
	// t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);

	carryprop(t2, C2, lo52mask);
	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	//scarry = _mm512_cmpge_epu64_mask(t2, n1);

	if (scarry > 0) {
		// conditionally subtract when needed (AMM - only on overflow)
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
__inline static void mask_sqrredc104_vec_pos(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with positive variant n'
	__m512i m;
	__m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#if !MP_USE_IFMA
	UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	UNUSED_VAR int biascount = 0;
	UNUSED_VAR __m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	UNUSED_VAR __mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	t1 = _mm512_add_epi64(t1, sqr_lo);

	// note, we leave rho = 0 - rho so that we get -m,
	// and thus the muladd becomes mulsub, since there is
	// no fmsub52 in avx512-ifma.
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_sub_epi64(t1, C1);
	t2 = _mm512_add_epi64(sqr_hi, C2);
	t0 = t1;
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_sub_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);
	carryprop(t2, C2, lo52mask);
	scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	if (scarry > 0) {
		// conditionally add
		t1 = _mm512_mask_add_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_add_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_add_epi64(t2, scarry, t2, _mm512_srli_epi64(t1, 52));
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
__inline static void mask_sqrredc104_exact_vec(__m512i* c1, __m512i* c0, __mmask8 mulmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1, t2, C3, C1, C2, sqr_lo, sqr_hi;

#if !MP_USE_IFMA
	UNUSED_VAR __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
	UNUSED_VAR __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
	UNUSED_VAR __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	UNUSED_VAR __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
	UNUSED_VAR __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
	UNUSED_VAR __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
	UNUSED_VAR int biascount = 0;
	UNUSED_VAR __m512i i0, i1;
#endif

	__m512i zero = _mm512_set1_epi64(0);
	UNUSED_VAR __m512i one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);
	UNUSED_VAR __mmask8 scarry2;
	__mmask8 scarry;

	t0 = t1 = t2 = C1 = C2 = C3 = sqr_lo = sqr_hi = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a0, sqr_lo, sqr_hi);
	t1 = sqr_lo;
	t2 = sqr_hi;
	VEC_MUL_ACCUM_LOHI_PD(a0, a0, t0, t1);

	// m0
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, C2);

	t1 = _mm512_add_epi64(t1, C1);
	t2 = _mm512_add_epi64(t2, C2);
	// we throw t0 away after this so first propagate its carry.
	// t0 = _mm512_add_epi64(t1, one);
	// t0 = _mm512_mask_add_epi64(t1, _mm512_cmpgt_epu64_mask(m, zero), t1, one);
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
	t1 = t2;
	t2 = C1 = C2 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a1, a1, t1, t2);

	t0 = _mm512_add_epi64(t0, sqr_lo);
	t1 = _mm512_add_epi64(t1, sqr_hi);

	// m1
	m = mul52lo(t0, vrho);

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, C1);
	VEC_MUL_ACCUM_LOHI_PD(m, n1, t1, t2);

	t1 = _mm512_add_epi64(t1, C1);

	// final carryprop
	carryprop(t0, t1, lo52mask);
	carryprop(t1, t2, lo52mask);

	//carryprop(t2, C2, lo52mask);
	//scarry = _mm512_cmp_epu64_mask(C2, zero, _MM_CMPINT_GT);

	// result >= n, as a two-word comparison.  (This used to test the high
	// words only -- t2 >= n1 -- which also subtracts when the high words are
	// equal and the low word is SMALLER, giving a wrong, wrapped result.  That
	// is a 2^-52 event for moduli near 104 bits but frequent when n1 is only a
	// few bits wide.)
	scarry = _mm512_cmpgt_epu64_mask(t2, n1) |
		(_mm512_cmpeq_epu64_mask(t2, n1) & _mm512_cmpge_epu64_mask(t1, n0));
	//scarry |= (_mm512_cmpeq_epu64_mask(t2, n1) & _mm512_cmpgt_epu64_mask(t1, n0));

	if (scarry > 0) {
		// conditionally subtract when result >= n
		C1 = _mm512_mask_set1_epi64(zero, _mm512_cmpgt_epi64_mask(n0, t1), 1);
		t1 = _mm512_mask_sub_epi64(t1, scarry, t1, n0);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, n1);
		t2 = _mm512_mask_sub_epi64(t2, scarry, t2, C1);
	}

	// on Zen4-epyc it is slower to do this:
	// *c0 = _mm512_mask_and_epi64(a0, mulmsk, lo52mask, t1);
	// *c1 = _mm512_mask_and_epi64(a1, mulmsk, lo52mask, t2);

	// than this:
	*c0 = _mm512_and_epi64(lo52mask, t1);
	*c1 = _mm512_and_epi64(lo52mask, t2);
	*c0 = _mm512_mask_mov_epi64(*c0, ~mulmsk, a0);
	*c1 = _mm512_mask_mov_epi64(*c1, ~mulmsk, a1);

	return;
}
#if defined(__GNUC__) && !MP_USE_IFMA
#pragma GCC diagnostic pop
#endif

__inline static void addmod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
	__m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	// add
	__mmask8 bmsk;
	//a0 = _mm512_addsetc_epi52(a0, b0, &bmsk);
	//a1 = _mm512_adc_epi52(a1, bmsk, b1, &bmsk);
	a0 = _mm512_add_epi64(a0, b0);
	a1 = _mm512_add_epi64(a1, b1);
	a1 = _mm512_add_epi64(a1, _mm512_srli_epi64(a0, 52));
	a0 = _mm512_and_epi64(a0, lo52mask);

	// compare
	//__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);
	// *c0 = _mm512_mask_sub_epi64(a0, msk, a0, n0);
	// *c1 = _mm512_mask_sub_epi64(a1, msk, a1, n1);
	// *c1 = _mm512_mask_sub_epi64(*c1, msk, *c1, _mm512_srli_epi64(*c0, 63));
	return;
}
__inline static void chkmod104_x8(__m512i* c1, __m512i* c0, 
	__m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// check if larger than n and reduce if so.
	__mmask8 bmsk = _mm512_cmpgt_epu64_mask(n0, b0);

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(b1, n1);
	__mmask8 msk2 = _mm512_cmpge_epu64_mask(b0, n0);
	msk |= (_mm512_cmpeq_epu64_mask(b1, n1) & msk2);

	// conditionally subtract N
	// *c0 = _mm512_mask_subsetc_epi52(b0, msk, b0, n0, &bmsk);
	// *c1 = _mm512_mask_sbb_epi52(b1, msk, bmsk, n1, &bmsk);
	*c0 = _mm512_mask_sub_epi64(b0, msk, b0, n0);
	*c1 = _mm512_mask_sub_epi64(b1, msk, b1, n1);
	*c0 = _mm512_and_epi64(*c0, _mm512_set1_epi64(0xfffffffffffffull));
	*c1 = _mm512_mask_sub_epi64(*c1, msk & bmsk, *c1, _mm512_set1_epi64(1));
	return;
}
__inline static void mask_addmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// add
	__mmask8 bmsk;
	a0 = _mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
	a1 = _mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

	// compare
	__mmask8 msk = bmsk | _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void mask_dblmod104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	// add
	__mmask8 bmsk;
	//a0 = _mm512_mask_addsetc_epi52(a0, addmsk, a0, b0, &bmsk);
	//a1 = _mm512_mask_adc_epi52(a1, addmsk, bmsk, b1, &bmsk);

	a0 = _mm512_mask_slli_epi64(a0, addmsk, a0, 1);
	a1 = _mm512_mask_slli_epi64(a1, addmsk, a1, 1);
	// when doubling, it is safe to check both carries before adding
	// in the previous carry, because the shift makes room for
	// the previous carry.  So either the upper word shift generates
	// a carry or doesn't, the addition won't cause one.
	a1 = _mm512_add_epi64(a1, _mm512_srli_epi64(a0, 52));
	a0 = _mm512_and_epi64(lo52mask, a0);

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void mask_redsub104_x8(__m512i* c1, __m512i* c0, __mmask8 addmsk,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__mmask8 bmsk;

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, addmsk & msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, addmsk & msk, bmsk, n1, &bmsk);
	return;
}
__inline static void redsub104_x8(__m512i* c1, __m512i* c0,
	__m512i a1, __m512i a0, __m512i n1, __m512i n0)
{
	__mmask8 bmsk;

	// compare
	__mmask8 msk = _mm512_cmpgt_epu64_mask(a1, n1);
	msk |= (_mm512_cmpeq_epu64_mask(a1, n1) & _mm512_cmpge_epu64_mask(a0, n0));

	// conditionally subtract N
	*c0 = _mm512_mask_subsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_sbb_epi52(a1, msk, bmsk, n1, &bmsk);

	return;
}
__inline static void submod104_x8(__m512i* c1, __m512i* c0, __m512i a1, __m512i a0,
	__m512i b1, __m512i b0, __m512i n1, __m512i n0)
{
	// compare
	__mmask8 msk = _mm512_cmplt_epu64_mask(a1, b1);
	msk |= _mm512_cmpeq_epu64_mask(a1, b1) & _mm512_cmplt_epu64_mask(a0, b0);

	// subtract
	__mmask8 bmsk;
	a0 = _mm512_subsetc_epi52(a0, b0, &bmsk);
	a1 = _mm512_sbb_epi52(a1, bmsk, b1, &bmsk);

	// conditionally add N
	*c0 = _mm512_mask_addsetc_epi52(a0, msk, a0, n0, &bmsk);
	*c1 = _mm512_mask_adc_epi52(a1, msk, bmsk, n1, &bmsk);
	return;
}

/* 8 lanes of the 52-bit inverse of a[i] mod 2^52 (a[i] odd) */
__m512i multiplicative_inverse104_x8(uint64_t* a);

#endif /* MP_USE_AVX512F && !LIMB2_NO_VECTOR && !LIMB1_SCALAR_ONLY */

#ifdef __cplusplus
}
#endif

#endif /* LIMB2_H */
