/* ===========================================================================
   limb2.c -- double-limb arithmetic: 128-bit scalar Montgomery.

   Moved verbatim from monty.c.  What changed in the move:

     - The two SETUP routines did NOT move.  to_monty128() needs GMP (mpz_t,
       17 references) to build the Montgomery constants, and monty128_init()
       calls it -- so although monty128_init has no mpz_ references of its own,
       the whole setup path is GMP-bound.  Both stay in monty.c with the rest
       of the GMP glue, which leaves limb2.o with NO undefined references and
       no GMP dependency: it is pure double-limb arithmetic.  Their
       declarations remain in limb2.h, since they are part of the 128-bit
       interface.
     - monty.c's local multiplicative_inverse() is GONE: limb1.h now holds the
       single canonical inline definition (verified byte-identical to the
       monty.c / uecm_ / upm1_ copies before consolidating).
     - monty.c's local _umul128 and its (already #if 0'd) local _addcarry_u64
       are GONE: mp_platform.h supplies both, with one uniform convention.
     - my_clz64 / my_ctz64 became _lead_full_zcnt / _trail_full_zcnt in
       mp_platform.h, and my_clz52 / my_ctz52 became _lead_full_zcnt52 /
       _trail_full_zcnt52 in limb1.h -- all four are now thin zero-safe
       wrappers over the single _lead_zcnt64 / _trail_zcnt64 ladder instead of
       four separate ~37-line copies of it.  Only the 104- and 128-bit
       wrappers remain here, retargeted to the new names.
     - The BMI asm branches in the bit-scan functions are DEAD, for TWO
       independent reasons, and moving the code changes neither:
         (1) `#if (INLINE_ASM && ...)` -- INLINE_ASM was #define'd at
             monty.c:609, *after* its uses at 421/460/516/564, and an
             undefined identifier is 0 in #if.  mp_platform.h now defines
             INLINE_ASM up front, so this reason alone would have gone away.
         (2) the inner `#ifdef __BMI1__` -- that macro does not exist.  gcc
             and clang define __BMI__ and __BMI2__; __BMI1__ is never set by
             any of them, so the branch is unreachable even with -mbmi or
             -march=native.  It looks like a typo for __BMI__.
       Verified empirically: building limb1.c with -mbmi, with and without the
       opt-in, produces byte-identical assembly.  The functions have always
       compiled their __builtin_ctzll/__builtin_clzll branch, which gcc already
       lowers to tzcnt/lzcnt anyway, so there is likely nothing to win here --
       but if you want to test it, fix the macro name first.
   =========================================================================== */

#include "limb2.h"
#include <stdbool.h>   /* my_sbb64 uses bool; monty.c included this at its top */

/* ---- Montgomery setup ------------------------------------------------ */
void ciosFullMul128x(uint64_t *u, uint64_t *v, uint64_t rho, uint64_t *n, uint64_t *w)
{
#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax

	ASM_G(
		"movq %0, %%r10	\n\t"			/* u ptr in r10 */
		"movq %2, %%r11	\n\t"			/* w ptr in r11 */
		"movq 0(%1), %%r9	\n\t"		/* v[0] ptr in r9 */

		/* begin s += u * v */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
		"mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(u[0] * v)         */
		"addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
		"mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(u[1] * v)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
		"addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
		"adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */

		"movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
		"mulx %4, %%r9, %%r14	\n\t"   /* m = rho * w[0]         */
		"movq %3, %%r10	\n\t"			/* n ptr in r10 */

		/* begin s = (s + n * m) >> 64 */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
		"mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(n[0] * m)         */
		"addq 0(%%r11), %%r12  \n\t"    /* r12 = w[0] (could be rdx) + LO(n[0] * m) */

		/* r12 should be 0 here */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
		"mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(n[1] * m)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
		"xorq %%r14, %%r14 \n\t"
		"addq 8(%%r11), %%r12  \n\t"
		"movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

		"adcq 16(%%r11), %%r13 \n\t"
		"movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
		"adcq $0, %%r14 \n\t"
		"movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */

		/* round 2 */

		"movq %0, %%r10	\n\t"			/* u ptr in r10 */
		"movq 8(%1), %%r9	\n\t"		/* v[1] ptr in r9 */

		/* begin s += u * v */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by u[0]  */
		"mulx %%r9, %%r12, %%r14 \n\t"  /* r14 = HI(u[0] * v)         */
		"addq %%r12, 0(%%r11) \n\t"     /* w[0] = w[0] + LO(u[0] * v) */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by u[1]  */
		"mulx %%r9, %%r12, %%r13 \n\t"  /* r13 = HI(u[1] * v)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(u[0] * v) + LO(u[1] * v) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(u[1] * v) + prevcarry                */
		"addq %%r12, 8(%%r11) \n\t"		/* w[1] = w[1] + HI(u[0] * v) + LO(u[1] * v)*/
		"adcq %%r13, 16(%%r11) \n\t"	/* w[2] = w[2] + HI(u[1] * v) + prevcarry */

		"movq 0(%%r11), %%rdx	\n\t"   /* ready to multiply by w[0]  */
		"mulx %4, %4, %%r14	\n\t"       /* m = rho * w[0]         */
		"movq %3, %%r10	\n\t"			/* n ptr in r10 */

		/* begin s = (s + n * m) >> 64 */
		"movq 0(%%r10), %%rdx	\n\t"   /* ready to multiply by n[0]  */
		"mulx %4, %%r12, %%r14	\n\t"   /* r14 = HI(n[0] * m)         */
		"addq 0(%%r11), %%r12  \n\t"    /* r12 = w[0] (could be rdx) + LO(n[0] * m) */

		/* r12 should be 0 here */

		"movq 8(%%r10), %%rdx	\n\t"   /* ready to multiply by n[1]  */
		"mulx %4, %%r12, %%r13	\n\t"   /* r13 = HI(n[1] * m)         */
		"adcq %%r14, %%r12 \n\t"        /* r12 = HI(n[0] * m) + LO(n[1] * m) + prevcarry */
		"adcq $0, %%r13 \n\t"           /* r13 = HI(n[1] * m) + prevcarry                */
		"xorq %%r14, %%r14 \n\t"
		"addq 8(%%r11), %%r12  \n\t"
		"movq %%r12, 0(%%r11)	\n\t"	/* w[0] = w[1] + HI(n[0] * m) + LO(n[1] * m) + prevcarry */

		"adcq 16(%%r11), %%r13 \n\t"
		"movq %%r13, 8(%%r11)	\n\t"	/* w[1] = w[2] + HI(n[1] * m) + prevcarry */
		"adcq $0, %%r14 \n\t"
		"movq %%r14, 16(%%r11)	\n\t"   /* w[2] = carry out */



		:
		: "r"(u), "r"(v), "r"(w), "r"(n), "r"(rho)
		: "r9", "r10", "rdx", "r11", "r12", "r13", "r14", "cc", "memory");

#else


#endif
	return;
}

/* ---- 104/128-bit bit scan --------------------------------------------
   These delegate to the zero-safe COUNT wrappers: _lead_full_zcnt /
   _trail_full_zcnt (mp_platform.h) and _lead_full_zcnt52 /
   _trail_full_zcnt52 (limb1.h).  Do NOT be tempted to route them through
   _lead_zcnt64: that one is undefined at zero on the builtin paths, and on
   the _BitScanReverse64 branches (MSVC, clang-cl, icc without BMI2) it
   returns the INDEX of the high bit rather than a leading-zero count. */

uint64_t my_ctz104(uint64_t n_lo, uint64_t n_hi)
{
	if (n_lo) {
		return _trail_full_zcnt52(n_lo);
	}
	return 52 + _trail_full_zcnt52(n_hi);
}

uint64_t my_ctz128(uint64_t n_lo, uint64_t n_hi)
{
	if (n_lo) {
		return _trail_full_zcnt(n_lo);
	}
	return 64 + _trail_full_zcnt(n_hi);
}

uint64_t my_clz104(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return _lead_full_zcnt52(n_hi);
	}
	return 52 + _lead_full_zcnt52(n_lo);
}

uint64_t my_clz128(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return _lead_full_zcnt(n_hi);
	}
	return 64 + _lead_full_zcnt(n_lo);
}

/* ---- CIOS cores ------------------------------------------------------ */

static inline uint8_t my_sbb64(uint8_t borrow_in, uint64_t a, uint64_t b, uint64_t* diff)
{
#if (INLINE_ASM && (defined( __INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)))
	return _subborrow_u64(borrow_in, a, b, (unsigned long long*)diff);
#elif defined(__GNUC__)
	bool c;
	c = __builtin_usubll_overflow(a, b, (unsigned long long*)diff);
	c |= __builtin_usubll_overflow(*diff, borrow_in, (unsigned long long*)diff);
	return c;
#elif defined(_MSC_VER)
	return _subborrow_u64(borrow_in, a, b, (unsigned long long*)diff);
#else
	if (__builtin_constant_p(borrow_in) && borrow_in == 0) {
		if (__builtin_constant_p(a) && a == 0) {
			*diff = -b;
			return 1;
		}
		else if (__builtin_constant_p(b) && b == 0) {
			*diff = a;
			return 0;
		}
		else {
			uint64_t tmp = a - b;
			uint8_t borrow = (tmp > a);
			*diff = tmp;
			return borrow;
		}
	}
	else {
		uint64_t tmp1 = a - borrow_in;
		uint8_t borrow = (tmp1 > a);
		if (__builtin_constant_p(b) && b == 0) {
			*diff = tmp1;
			return borrow;
		}
		else {
			uint64_t tmp2 = tmp1 - b;
			borrow |= (tmp2 > tmp1);
			*diff = tmp2;
			return borrow;
		}
	}
#endif
}

static void ciosSubtract128(uint64_t* res_lo, uint64_t* res_hi, uint64_t carries, uint64_t mod_lo, uint64_t mod_hi)
{
	uint64_t n_lo, n_hi;
	uint64_t t_lo, t_hi;
	uint8_t b;
	n_lo = *res_lo;
	n_hi = *res_hi;
	// save, subtract the modulus until a borrows occurs
	do {
		t_lo = n_lo;
		t_hi = n_hi;
		b = my_sbb64(0, n_lo, mod_lo, &n_lo);
		b = my_sbb64(b, n_hi, mod_hi, &n_hi);
#ifdef _MSC_VER
		b = my_sbb64(b, carries, 0, &carries);
#else
		if (__builtin_constant_p(carries) && carries == 0) {
		}
		else {
			b = my_sbb64(b, carries, 0, &carries);
		}
#endif

	} while (b == 0);
	// get the saved values when a borrow occurs
	*res_lo = t_lo;
	*res_hi = t_hi;
}

#if 1 //def USE_PERIG_128BIT



void ciosModMul128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
	uint64_t mmagic)
{

#ifdef _MSC_VER
	uint64_t a_lo = *res_lo, a_hi = *res_hi;
	uint64_t cshi, cslo, cchi, cclo;
	uint64_t t0, t1, t2, t3, m;

	//cc = (uint128_t)a_lo * b_lo;	// #1
	//t0 = (uint64_t)cc;
	//cc = cc >> 64;
	cclo = _umul128(a_lo, b_lo, &cchi);
	t0 = cclo;
	cclo = cchi;

	//cc += (uint128_t)a_lo * b_hi;	// #2
	cslo = _umul128(a_lo, b_hi, &cshi);
	cchi = _addcarry_u64(0, cclo, cslo, &cclo);
	cchi += cshi;
	
	//t1 = (uint64_t)cc;
	//cc = cc >> 64;
	//t2 = (uint64_t)cc;
	t1 = cclo;
	t2 = cchi;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	//cs = (uint128_t)m * mod_lo;	// #4
	cslo = _umul128(m, mod_lo, &cshi);

	//cs += t0;
	//cs = cs >> 64;
	cshi += _addcarry_u64(0, t0, cslo, &cslo);
	cslo = cshi;

	//cs += (uint128_t)m * mod_hi;	// #5
	//cs += t1;
	cclo = _umul128(m, mod_hi, &cchi);
	cchi += _addcarry_u64(0, cclo, cslo, &cclo);
	cchi += _addcarry_u64(0, t1, cclo, &cclo);

	//t0 = (uint64_t)cs;
	//cs = cs >> 64;
	t0 = cclo;
	cslo = cchi;
	
	//cs += t2;
	cshi = _addcarry_u64(0, t2, cslo, &cslo);

	//t1 = (uint64_t)cs;
	//cs = cs >> 64;
	//t2 = (uint64_t)cs;
	t1 = cslo;
	t2 = cshi;

#if PARANOID
	assert(cs >> 64 == 0);
#endif

	//cc = (uint128_t)a_hi * b_lo;	// #6
	//cc += t0;
	//t0 = (uint64_t)cc;
	//cc = cc >> 64;
	cclo = _umul128(a_hi, b_lo, &cchi);
	cchi += _addcarry_u64(0, cclo, t0, &cclo);
	t0 = cclo;
	cclo = cchi;
	
	//cc += (uint128_t)a_hi * b_hi;	// #7
	//cc += t1;
	//t1 = (uint64_t)cc;
	//cc = cc >> 64;
	cslo = _umul128(a_hi, b_hi, &cshi);
	cshi += _addcarry_u64(0, cclo, cslo, &cslo);
	cshi += _addcarry_u64(0, cslo, t1, &cslo);
	t1 = cslo;
	cclo = cshi;

	//cc += t2;
	//t2 = (uint64_t)cc;
	//cc = cc >> 64;
	//t3 = (uint64_t)cc;
	cchi = _addcarry_u64(0, cclo, t2, &cclo);
	t2 = cclo;
	t3 = cchi;

#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	//cs = (uint128_t)m * mod_lo;	// #9
	//cs += t0;
	//cs = cs >> 64;
	cslo = _umul128(m, mod_lo, &cshi);
	cshi += _addcarry_u64(0, t0, cslo, &cslo);
	cslo = cshi;

	//cs += (uint128_t)m * mod_hi;	// #10
	//cs += t1;
	cclo = _umul128(m, mod_hi, &cchi);
	cchi += _addcarry_u64(0, cclo, cslo, &cclo);
	cchi += _addcarry_u64(0, t1, cclo, &cclo);

	//t0 = (uint64_t)cs;
	//cs = cs >> 64;
	t0 = cclo;
	cslo = cchi;

	//cs += t2;
	cshi = _addcarry_u64(0, t2, cslo, &cslo);

	//t1 = (uint64_t)cs;
	//cs = cs >> 64;
	t1 = cslo;
	cslo = cshi;

	//cs += t3;
	//t2 = (uint64_t)cs;
	cshi = _addcarry_u64(0, t3, cslo, &cslo);
	t2 = cslo;

	if (t2) {
		unsigned char carry = _subborrow_u64(0, t0, mod_lo, &t0);
		_subborrow_u64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}


#else
	uint64_t a_lo = *res_lo, a_hi = *res_hi;
	uint128_t cs, cc;
	uint64_t t0, t1, t2, t3, m;

	cc = (uint128_t)a_lo * b_lo;	// #1
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)a_lo * b_hi;	// #2
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	t2 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #5
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	t2 = (uint64_t)cs;
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	cc = (uint128_t)a_hi * b_lo;	// #6
	cc += t0;
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)a_hi * b_hi;	// #7
	cc += t1;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t)cc;
	cc = cc >> 64;
	t3 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	cs = (uint128_t)m * mod_lo;	// #9
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #10
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t)cs;
	ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);

	//if (t2) {
	//	unsigned char carry = my_sbb64(0, t0, mod_lo, &t0);
	//	my_sbb64(carry, t1, mod_hi, &t1);
	//	ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	//}

#endif
	*res_lo = t0;
	*res_hi = t1;
}

/********************* end of Perig's 128-bit code **********************/

// modSqr version I wrote based on the modMul
void ciosModSqr128(uint64_t* res_lo, uint64_t* res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
	uint64_t mmagic)
{
#if 1 //def _MSC_VER
	*res_lo = b_lo;
	*res_hi = b_hi;
	ciosModMul128(res_lo, res_hi, b_lo, b_hi, mod_lo, mod_hi, mmagic);
#else
	uint128_t cs, cc, b_lohi;
	uint64_t t0, t1, t2, t3, m;

	cc = (uint128_t)b_lo * b_lo;	// #1
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	b_lohi = (uint128_t)b_lo * b_hi;	// #2
	cc += b_lohi;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	t2 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t)m * mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #5
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	t2 = (uint64_t)cs;
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	cc = b_lohi + t0;
	t0 = (uint64_t)cc;
	cc = cc >> 64;
	cc += (uint128_t)b_hi * b_hi;	// #6
	cc += t1;
	t1 = (uint64_t)cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t)cc;
	cc = cc >> 64;
	t3 = (uint64_t)cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	cs = (uint128_t)m * mod_lo;	// #9
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t)m * mod_hi;	// #10
	cs += t1;
	t0 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t)cs;
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t)cs;
	if (t2) {
		unsigned char carry = my_sbb64(0, t0, mod_lo, &t0);
		my_sbb64(carry, t1, mod_hi, &t1);
		//ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}

	*res_lo = t0;
	*res_hi = t1;

#endif

}
#endif

/* ---- 128-bit modular interface --------------------------------------- */

void mulmod128(uint64_t * u, uint64_t * v, uint64_t * w, monty128_t *mdata)
{
#ifdef USE_PERIG_128BIT
	uint64_t t[2];
	t[0] = v[0];
	t[1] = v[1];
	ciosModMul128(&t[0], &t[1], u[0], u[1], mdata->n[0], mdata->n[1], mdata->rho);
	w[0] = t[0];
	w[1] = t[1];
	return;

#else


	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64_t s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, v, mdata->rho, mdata->n, s);

	if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
	{
		ASM_G(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
			: "r"(s[0]), "r"(s[1]), "r"(mdata->n[0]), "r"(mdata->n[1]), "r"(w)
			: "r11", "cc", "memory");
	}
	else
	{
		w[0] = s[0];
		w[1] = s[1];
	}
#else
    // TODO: implement portable u128 x u128 modular multiplication
    uint64_t t[3], U, c2, c3;
    uint8_t c1, c4;

    // z = 0
    // for (i = 0; i < t; i++)
    // {
    //     u = (z0 + xi * y0) * -m’ mod b
    //     z = (z + xi * y + u * m) / b
    // }
    // if (x >= m) { x -= m; }

    //printf("nhat = %" PRIx64 "\n", mdata->rho);
    //printf("a = %" PRIx64 ",%" PRIx64 "\n", u[1], u[0]);
    //printf("b = %" PRIx64 ",%" PRIx64 "\n", v[1], v[0]);
    //printf("n = %" PRIx64 ",%" PRIx64 "\n", mdata->n[1], mdata->n[0]);

    //i = 0;
    //u = (z0 + x0 * y0) * nhat mod b;
    t[0] = _umul128(u[0], v[0], &t[1]);
    U = t[0] * mdata->rho;

    //printf("t[0] = %" PRIx64 ", t[1] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], t[1], U);

    //j = 0;
    //z = (z + x0 * y0 + u * m0) / b;
    c1 = _addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    t[2] = _addcarry_u64(c1, t[1], c2, &t[1]);                              // c1 applies to t2
    
    //j = 1;
    //z = (z + x0 * y1 + u * m1) / b;
    c1 = _addcarry_u64(0, _umul128(u[0], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = _addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2
    t[0] = t[1];                                                            // divide by b
    t[1] = t[2];
    t[2] = c1 + c4;

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    //i = 0;
    //u = (z + x1 * y0) * nhat mod b;
    c1 = _addcarry_u64(0, _umul128(u[1], v[0], &c2), t[0], &t[0]);          // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;
    U = t[0] * mdata->rho;

    //printf("t[0] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], U);

    //j = 0;
    //z = (z + x1 * y0 + u * m0) / b;
    c1 = _addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;

    //j = 1;
    //z = (z + x1 * y1 + u * m1) / b;
    c1 = _addcarry_u64(0, _umul128(u[1], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = _addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    w[0] = t[1];
    w[1] = t[2];

    //exit(1);

    return;



#endif

#endif
	return;
}

void mulmod128n(uint64_t* u, uint64_t* v, uint64_t* w, uint64_t* n, uint64_t rho)
{
#ifdef USE_PERIG_128BIT
	uint64_t t[2];
	t[0] = v[0];
	t[1] = v[1];
	ciosModMul128(&t[0], &t[1], u[0], u[1], n[0], n[1], rho);
	w[0] = t[0];
	w[1] = t[1];
	return;

#else


	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64_t s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
	// requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, v, rho, n, s);

	if ((s[2]) || (s[1] > n[1]) || ((s[1] == n[1]) && (s[0] > n[0])))
	{
		ASM_G(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
		: "r"(s[0]), "r"(s[1]), "r"(n[0]), "r"(n[1]), "r"(w)
			: "r11", "cc", "memory");
	}
	else
	{
		w[0] = s[0];
		w[1] = s[1];
	}
#else
	// TODO: implement portable u128 x u128 modular multiplication
	uint64_t t[3], U, c2, c3;
	uint8_t c1, c4;

	// z = 0
	// for (i = 0; i < t; i++)
	// {
	//     u = (z0 + xi * y0) * -m’ mod b
	//     z = (z + xi * y + u * m) / b
	// }
	// if (x >= m) { x -= m; }

	//printf("nhat = %" PRIx64 "\n", mdata->rho);
	//printf("a = %" PRIx64 ",%" PRIx64 "\n", u[1], u[0]);
	//printf("b = %" PRIx64 ",%" PRIx64 "\n", v[1], v[0]);
	//printf("n = %" PRIx64 ",%" PRIx64 "\n", mdata->n[1], mdata->n[0]);

	//i = 0;
	//u = (z0 + x0 * y0) * nhat mod b;
	t[0] = _umul128(u[0], v[0], &t[1]);
	U = t[0] * rho;

	//printf("t[0] = %" PRIx64 ", t[1] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], t[1], U);

	//j = 0;
	//z = (z + x0 * y0 + u * m0) / b;
	c1 = _addcarry_u64(0, _umul128(U, n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
	t[2] = _addcarry_u64(c1, t[1], c2, &t[1]);                              // c1 applies to t2

	//j = 1;
	//z = (z + x0 * y1 + u * m1) / b;
	c1 = _addcarry_u64(0, _umul128(u[0], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
	c4 = _addcarry_u64(0, _umul128(U, n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
	c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
	c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2
	t[0] = t[1];                                                            // divide by b
	t[1] = t[2];
	t[2] = c1 + c4;

	//printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

	//i = 0;
	//u = (z + x1 * y0) * nhat mod b;
	c1 = _addcarry_u64(0, _umul128(u[1], v[0], &c2), t[0], &t[0]);          // c1,c2 apply to t1
	c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
	t[2] += c1;
	U = t[0] * rho;

	//printf("t[0] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], U);

	//j = 0;
	//z = (z + x1 * y0 + u * m0) / b;
	c1 = _addcarry_u64(0, _umul128(U, n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
	c1 = _addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
	t[2] += c1;

	//j = 1;
	//z = (z + x1 * y1 + u * m1) / b;
	c1 = _addcarry_u64(0, _umul128(u[1], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
	c4 = _addcarry_u64(0, _umul128(U, n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
	c1 = _addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
	c4 = _addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2

	//printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

	w[0] = t[1];
	w[1] = t[2];

	//exit(1);

	return;



#endif

#endif
	return;
}

void sqrmod128(uint64_t * u, uint64_t * w, monty128_t *mdata)
{
#ifdef USE_PERIG_128BIT
	ciosModSqr128(&w[0], &w[1], u[0], u[1], mdata->n[0], mdata->n[1], mdata->rho);
	return;

#else


	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64_t s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, u, mdata->rho, mdata->n, s);

	if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
	{
		ASM_G(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
			: "r"(s[0]), "r"(s[1]), "r"(mdata->n[0]), "r"(mdata->n[1]), "r"(w)
			: "r11", "cc", "memory");
	}
	else
	{
		w[0] = s[0];
		w[1] = s[1];
	}

#else

    mulmod128(u, u, w, mdata);

#endif

#endif
	return;
}

void sqrmod128n(uint64_t* u, uint64_t* w, uint64_t* n, uint64_t rho)
{
#ifdef USE_PERIG_128BIT
	ciosModSqr128(&w[0], &w[1], u[0], u[1], n[0], n[1], rho);
	return;

#else


	// integrate multiply and reduction steps, alternating
	// between iterations of the outer loops.
	uint64_t s[3];

	s[0] = 0;
	s[1] = 0;
	s[2] = 0;

#if defined( USE_AVX2 ) && defined(GCC_ASM64X)
	// requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, u, rho, n, s);

	if ((s[2]) || (s[1] > n[1]) || ((s[1] == n[1]) && (s[0] > n[0])))
	{
		ASM_G(
			"movq %4, %%r11 \n\t"
			"movq %0, 0(%%r11) \n\t"
			"movq %1, 8(%%r11) \n\t"
			"subq %2, 0(%%r11) \n\t"
			"sbbq %3, 8(%%r11) \n\t"
			:
		: "r"(s[0]), "r"(s[1]), "r"(n[0]), "r"(n[1]), "r"(w)
			: "r11", "cc", "memory");
	}
	else
	{
		w[0] = s[0];
		w[1] = s[1];
	}

#else

	mulmod128n(u, u, w, n, rho);

#endif

#endif
	return;
}

void addmod128(uint64_t * a, uint64_t * b, uint64_t * w, uint64_t * n)
{
#if defined(GCC_ASM64X)
    // requires GCC_ASM64 syntax
	w[1] = a[1];
	w[0] = a[0];
	ASM_G(
		"movq %0, %%r8 \n\t"
		"movq %1, %%r9 \n\t"
		"subq %4, %%r8 \n\t"		/* t = x - n */
		"sbbq %5, %%r9 \n\t"
		"addq %2, %0 \n\t"			/* x += y */
		"adcq %3, %1 \n\t"
		"addq %2, %%r8 \n\t"		/* t += y */
		"adcq %3, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+r"(w[0]), "+r"(w[1])
		: "r"(b[0]), "r"(b[1]), "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

#else

    uint8_t c;
    uint64_t t[2];
    c = _addcarry_u64(0, a[0], b[0], &t[0]);
    c = _addcarry_u64(c, a[1], b[1], &t[1]);
    if (c || (t[1] > n[1]) || ((t[1] == n[1]) && (t[0] > n[0])))
    {
        c = _subborrow_u64(0, t[0], n[0], &w[0]);
        c = _subborrow_u64(c, t[1], n[1], &w[1]);
    }
    else
    {
        w[0] = t[0];
        w[1] = t[1];
    }


#endif
	return;
}

void submod128(uint64_t * a, uint64_t * b, uint64_t * w, uint64_t * n)
{
#if defined(GCC_ASM64X)
    // requires GCC_ASM64 syntax
	ASM_G(
		"movq %6, %%r11 \n\t"
		"xorq %%r8, %%r8 \n\t"
		"xorq %%r9, %%r9 \n\t"
		"movq %0, 0(%%r11) \n\t"
		"movq %1, 8(%%r11) \n\t"
		"subq %2, 0(%%r11) \n\t"
		"sbbq %3, 8(%%r11) \n\t"
		"cmovc %4, %%r8 \n\t"
		"cmovc %5, %%r9 \n\t"
		"addq %%r8, 0(%%r11) \n\t"
		"adcq %%r9, 8(%%r11) \n\t"
		"1: \n\t"
		:
	: "r"(a[0]), "r"(a[1]), "r"(b[0]), "r"(b[1]), "r"(n[0]), "r"(n[1]), "r"(w)
		: "r8", "r9", "r11", "cc", "memory");

#else

    uint8_t c;
    uint64_t t[2];
    c = _subborrow_u64(0, a[0], b[0], &t[0]);
    c = _subborrow_u64(c, a[1], b[1], &t[1]);
    if (c)
    {
        c = _addcarry_u64(0, t[0], n[0], &w[0]);
        c = _addcarry_u64(c, t[1], n[1], &w[1]);
    }
    else
    {
        w[0] = t[0];
        w[1] = t[1];
    }

#endif

	return;
}

void dblmod128(uint64_t* a, uint64_t* n)
{

#if defined(__x86_64__) && defined(__GNUC__)
	// requires GCC_ASM64 syntax
	uint64_t t1, t0;

	// do the 2-word variant of this:
	//uint64_t r;
	//uint64_t tmp = x - n;
	//uint8_t c = _addcarry_u64(0, tmp, y, &r);
	//return (c) ? r : x + y;
	t1 = a[1];
	t0 = a[0];

	__asm__ volatile (
		"movq %%rax, %%r8 \n\t"
		"movq %%rdx, %%r9 \n\t"
		"subq %4, %%r8 \n\t"		/* t = x - n */
		"sbbq %5, %%r9 \n\t"
		"addq %%rax, %0 \n\t"		/* x += x */
		"adcq %%rdx, %1 \n\t"
		"addq %%rax, %%r8 \n\t"		/* t = t + x */
		"adcq %%rdx, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+&r"(a[0]), "+&r"(a[1])
		: "a"(t0), "d"(t1), "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

#else

#ifdef USE_PERIG_128BIT

	uint128_t x = ((uint128_t)a[1] << 64) | a[0];
	uint128_t n128 = ((uint128_t)n[1] << 64) | n[0];
	uint128_t r = (x >= n128 - x) ? x - (n128 - x) : x + x;
	a[1] = (uint64_t)(r >> 64);
	a[0] = (uint64_t)r;

#else

	addmod128(a, a, a, n);

#endif

#endif
	return;
}

void chkmod128(uint64_t* a, uint64_t* n)
{

#if defined(__x86_64__) && defined(__unix__)

	// not really faster than the branching version below
	// but probably should be checked on more cpus.
	uint64_t t0, t1;
	t0 = a[0];
	t1 = a[1];
	__asm__ volatile (
		"xorq %%r8, %%r8 \n\t"
		"xorq %%r9, %%r9 \n\t"
		"subq %2, %%r8 \n\t"		/* t = 0 - n */
		"sbbq %3, %%r9 \n\t"
		"addq %0, %%r8 \n\t"		/* t += x */
		"adcq %1, %%r9 \n\t"
		"cmovc %%r8, %0 \n\t"
		"cmovc %%r9, %1 \n\t"
		: "+&r"(t0), "+&r"(t1)
		: "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc", "memory");

	a[0] = t0;
	a[1] = t1;
#else

	if ((a[1] > n[1]) || ((a[1] == n[1]) && (a[0] >= n[0])))
	{
		uint8_t c1 = _subborrow_u64(0, a[0], n[0], &a[0]);
		_subborrow_u64(c1, a[1], n[1], &a[1]);
	}

#endif
	return;
}

