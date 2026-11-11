/* ===========================================================================
   limb2.c -- double-limb arithmetic: out-of-line implementation.

     - 128-bit Montgomery multiply / square (three interchangeable cores:
       mulx asm, Perig Aubry's __int128 CIOS, and a portable CIOS built on
       the mp_platform.h polyfills), modular add / sub / double
     - 104 / 128-bit zero-safe bit scans
     - plain 128 / 128 division (Knuth D) and reduction
     - 128-bit Montgomery setup (no GMP)
     - the two AVX-512 helpers that are too big or too cold to inline

   The Montgomery cores were moved verbatim from monty.c.  Notes from that
   move that still matter:

     - my_clz64 / my_ctz64 became _lead_full_zcnt / _trail_full_zcnt in
       mp_platform.h, and my_clz52 / my_ctz52 became _lead_full_zcnt52 /
       _trail_full_zcnt52 in limb1.h -- thin zero-safe wrappers over the one
       _lead_zcnt64 / _trail_zcnt64 ladder.  limb2.h keeps the old names as
       macros.
     - monty.c's BMI asm branches in the bit-scan functions were dead code
       (gated on an INLINE_ASM that was defined after its uses, and on
       __BMI1__, a macro no compiler defines).  They were dropped; the
       builtins already lower to tzcnt / lzcnt.
   =========================================================================== */

#include "limb2.h"
#include <stddef.h>
#include <stdbool.h>   /* my_sbb64 uses bool */

/* ---- Montgomery setup ------------------------------------------------ */
void ciosFullMul128x(uint64_t *u, uint64_t *v, uint64_t rho, uint64_t *n, uint64_t *w)
{
#if (MP_USE_AVX2 || MP_USE_BMI2) && MP_HAS_GNU_ASM_X64
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax

	__asm__(
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

	// No mulx / GNU asm on this build.  This branch used to be empty, so the
	// function silently did nothing; keep the contract instead (w[0..2] is
	// u*v*R^-1 mod n, possibly plus n) by going through the general
	// multiplier, which returns a fully reduced two-limb value.
	mulmod128n(u, v, w, n, rho);
	w[2] = 0;

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
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
	return mp_subborrow_u64(borrow_in, a, b, diff);
#elif defined(__GNUC__)
	bool c;
	c = __builtin_usubll_overflow(a, b, (unsigned long long*)diff);
	c |= __builtin_usubll_overflow(*diff, borrow_in, (unsigned long long*)diff);
	return c;
#elif defined(_MSC_VER)
	return mp_subborrow_u64(borrow_in, a, b, diff);
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

static UNUSED_FUNC void ciosSubtract128(uint64_t* res_lo, uint64_t* res_hi, uint64_t carries, uint64_t mod_lo, uint64_t mod_hi)
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
	cchi = mp_addcarry_u64(0, cclo, cslo, &cclo);
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
	cshi += mp_addcarry_u64(0, t0, cslo, &cslo);
	cslo = cshi;

	//cs += (uint128_t)m * mod_hi;	// #5
	//cs += t1;
	cclo = _umul128(m, mod_hi, &cchi);
	cchi += mp_addcarry_u64(0, cclo, cslo, &cclo);
	cchi += mp_addcarry_u64(0, t1, cclo, &cclo);

	//t0 = (uint64_t)cs;
	//cs = cs >> 64;
	t0 = cclo;
	cslo = cchi;
	
	//cs += t2;
	cshi = mp_addcarry_u64(0, t2, cslo, &cslo);

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
	cchi += mp_addcarry_u64(0, cclo, t0, &cclo);
	t0 = cclo;
	cclo = cchi;
	
	//cc += (uint128_t)a_hi * b_hi;	// #7
	//cc += t1;
	//t1 = (uint64_t)cc;
	//cc = cc >> 64;
	cslo = _umul128(a_hi, b_hi, &cshi);
	cshi += mp_addcarry_u64(0, cclo, cslo, &cslo);
	cshi += mp_addcarry_u64(0, cslo, t1, &cslo);
	t1 = cslo;
	cclo = cshi;

	//cc += t2;
	//t2 = (uint64_t)cc;
	//cc = cc >> 64;
	//t3 = (uint64_t)cc;
	cchi = mp_addcarry_u64(0, cclo, t2, &cclo);
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
	cshi += mp_addcarry_u64(0, t0, cslo, &cslo);
	cslo = cshi;

	//cs += (uint128_t)m * mod_hi;	// #10
	//cs += t1;
	cclo = _umul128(m, mod_hi, &cchi);
	cchi += mp_addcarry_u64(0, cclo, cslo, &cclo);
	cchi += mp_addcarry_u64(0, t1, cclo, &cclo);

	//t0 = (uint64_t)cs;
	//cs = cs >> 64;
	t0 = cclo;
	cslo = cchi;

	//cs += t2;
	cshi = mp_addcarry_u64(0, t2, cslo, &cslo);

	//t1 = (uint64_t)cs;
	//cs = cs >> 64;
	t1 = cslo;
	cslo = cshi;

	//cs += t3;
	//t2 = (uint64_t)cs;
	cshi = mp_addcarry_u64(0, t3, cslo, &cslo);
	t2 = cslo;

	if (t2) {
		unsigned char carry = mp_subborrow_u64(0, t0, mod_lo, &t0);
		mp_subborrow_u64(carry, t1, mod_hi, &t1);
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

#if (MP_USE_AVX2 || MP_USE_BMI2) && MP_HAS_GNU_ASM_X64
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, v, mdata->rho, mdata->n, s);

	if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
	{
		__asm__(
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
    //     u = (z0 + xi * y0) * -m' mod b
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
    c1 = mp_addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    t[2] = mp_addcarry_u64(c1, t[1], c2, &t[1]);                              // c1 applies to t2
    
    //j = 1;
    //z = (z + x0 * y1 + u * m1) / b;
    c1 = mp_addcarry_u64(0, _umul128(u[0], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = mp_addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = mp_addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = mp_addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2
    t[0] = t[1];                                                            // divide by b
    t[1] = t[2];
    t[2] = c1 + c4;

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    //i = 0;
    //u = (z + x1 * y0) * nhat mod b;
    c1 = mp_addcarry_u64(0, _umul128(u[1], v[0], &c2), t[0], &t[0]);          // c1,c2 apply to t1
    c1 = mp_addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;
    U = t[0] * mdata->rho;

    //printf("t[0] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], U);

    //j = 0;
    //z = (z + x1 * y0 + u * m0) / b;
    c1 = mp_addcarry_u64(0, _umul128(U, mdata->n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
    c1 = mp_addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
    t[2] += c1;

    //j = 1;
    //z = (z + x1 * y1 + u * m1) / b;
    c1 = mp_addcarry_u64(0, _umul128(u[1], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
    c4 = mp_addcarry_u64(0, _umul128(U, mdata->n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
    c1 = mp_addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
    c4 = mp_addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2

    //printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

    // The result is c1 + c4 : t[2] : t[1], in [0, 2n).  This path used to drop
    // the top carry and skip the final subtraction, which is wrong for moduli
    // near 2^128 and left results unreduced otherwise; the asm and __int128
    // paths have always returned a value in [0, n).  Do the same here.
    if ((c1 | c4) || (t[2] > mdata->n[1]) || ((t[2] == mdata->n[1]) && (t[1] >= mdata->n[0])))
    {
        c1 = mp_subborrow_u64(0, t[1], mdata->n[0], &t[1]);
        (void)mp_subborrow_u64(c1, t[2], mdata->n[1], &t[2]);
    }

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

#if (MP_USE_AVX2 || MP_USE_BMI2) && MP_HAS_GNU_ASM_X64
	// requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, v, rho, n, s);

	if ((s[2]) || (s[1] > n[1]) || ((s[1] == n[1]) && (s[0] > n[0])))
	{
		__asm__(
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
	//     u = (z0 + xi * y0) * -m' mod b
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
	c1 = mp_addcarry_u64(0, _umul128(U, n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
	t[2] = mp_addcarry_u64(c1, t[1], c2, &t[1]);                              // c1 applies to t2

	//j = 1;
	//z = (z + x0 * y1 + u * m1) / b;
	c1 = mp_addcarry_u64(0, _umul128(u[0], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
	c4 = mp_addcarry_u64(0, _umul128(U, n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
	c1 = mp_addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
	c4 = mp_addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2
	t[0] = t[1];                                                            // divide by b
	t[1] = t[2];
	t[2] = c1 + c4;

	//printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

	//i = 0;
	//u = (z + x1 * y0) * nhat mod b;
	c1 = mp_addcarry_u64(0, _umul128(u[1], v[0], &c2), t[0], &t[0]);          // c1,c2 apply to t1
	c1 = mp_addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
	t[2] += c1;
	U = t[0] * rho;

	//printf("t[0] = %" PRIx64 ", u = %" PRIx64 "\n", t[0], U);

	//j = 0;
	//z = (z + x1 * y0 + u * m0) / b;
	c1 = mp_addcarry_u64(0, _umul128(U, n[0], &c2), t[0], &t[0]);      // c1,c2 apply to t1
	c1 = mp_addcarry_u64(c1, t[1], c2, &t[1]);                                // c1 applies to t2
	t[2] += c1;

	//j = 1;
	//z = (z + x1 * y1 + u * m1) / b;
	c1 = mp_addcarry_u64(0, _umul128(u[1], v[1], &c2), t[1], &t[1]);          // c1,c2 apply to t1
	c4 = mp_addcarry_u64(0, _umul128(U, n[1], &c3), t[1], &t[1]);      // c1,c2 apply to t1
	c1 = mp_addcarry_u64(c1, t[2], c2, &t[2]);                                // c1 applies to t2
	c4 = mp_addcarry_u64(c4, t[2], c3, &t[2]);                                // c4 applies to t2

	//printf("t = %" PRIx64 ",%" PRIx64 ",%" PRIx64 "\n", t[2], t[1], t[0]);

	// The result is c1 + c4 : t[2] : t[1], in [0, 2n).  This path used to drop
	// the top carry and skip the final subtraction, which is wrong for moduli
	// near 2^128 and left results unreduced otherwise; the asm and __int128
	// paths have always returned a value in [0, n).  Do the same here.
	if ((c1 | c4) || (t[2] > n[1]) || ((t[2] == n[1]) && (t[1] >= n[0])))
	{
	    c1 = mp_subborrow_u64(0, t[1], n[0], &t[1]);
	    (void)mp_subborrow_u64(c1, t[2], n[1], &t[2]);
	}

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

#if (MP_USE_AVX2 || MP_USE_BMI2) && MP_HAS_GNU_ASM_X64
    // requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, u, mdata->rho, mdata->n, s);

	if ((s[2]) || (s[1] > mdata->n[1]) || ((s[1] == mdata->n[1]) && (s[0] > mdata->n[0])))
	{
		__asm__(
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

#if (MP_USE_AVX2 || MP_USE_BMI2) && MP_HAS_GNU_ASM_X64
	// requires mulx in BMI2 (via the AVX2 macro) and GCC_ASM64 syntax
	ciosFullMul128x(u, u, rho, n, s);

	if ((s[2]) || (s[1] > n[1]) || ((s[1] == n[1]) && (s[0] > n[0])))
	{
		__asm__(
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
#if MP_HAS_GNU_ASM_X64
    // requires GCC_ASM64 syntax
	// work in locals and store at the end: w may alias a, b, or both.  (This
	// used to copy a into w first, which destroyed b when w == b.)  The
	// outputs are early-clobber so an input can never share their register.
	uint64_t t0 = a[0], t1 = a[1];
	__asm__(
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
		: "+&r"(t0), "+&r"(t1)
		: "r"(b[0]), "r"(b[1]), "r"(n[0]), "r"(n[1])
		: "r8", "r9", "cc");
	w[0] = t0;
	w[1] = t1;

#else

    uint8_t c;
    uint64_t t[2];
    c = mp_addcarry_u64(0, a[0], b[0], &t[0]);
    c = mp_addcarry_u64(c, a[1], b[1], &t[1]);
    if (c || (t[1] > n[1]) || ((t[1] == n[1]) && (t[0] >= n[0])))     // >=: a sum equal to n must reduce to 0
    {
        c = mp_subborrow_u64(0, t[0], n[0], &w[0]);
        c = mp_subborrow_u64(c, t[1], n[1], &w[1]);
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
#if MP_HAS_GNU_ASM_X64
    // requires GCC_ASM64 syntax
	__asm__(
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
    c = mp_subborrow_u64(0, a[0], b[0], &t[0]);
    c = mp_subborrow_u64(c, a[1], b[1], &t[1]);
    if (c)
    {
        c = mp_addcarry_u64(0, t[0], n[0], &w[0]);
        c = mp_addcarry_u64(c, t[1], n[1], &w[1]);
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
	//uint8_t c = mp_addcarry_u64(0, tmp, y, &r);
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
		uint8_t c1 = mp_subborrow_u64(0, a[0], n[0], &a[0]);
		mp_subborrow_u64(c1, a[1], n[1], &a[1]);
	}

#endif
	return;
}



/* ===========================================================================
   Plain (non-Montgomery) two-limb division and reduction.
   All operands are little-endian limb pairs: x[0] = low 64 bits, x[1] = high.
   =========================================================================== */

/* q = n / d, r = n % d for 128-bit n, d.  Either of q, r may be NULL.
   Returns 0 on success, -1 on division by zero.

   Knuth algorithm D with 64-bit digits.  When d fits in one limb this is two
   chained 128/64 divides; otherwise the quotient is a single 64-bit digit:
   one 128/64 divide for the estimate, then at most two corrections.  This
   replaces the ~650-line generic version that used to live in arith.c, and
   the mod-only variant in mod128.c. */
int uint128_div(const uint64_t n[2], const uint64_t d[2],
    uint64_t q[2], uint64_t r[2])
{
    uint64_t q0, q1, r0, r1;

    if ((d[0] | d[1]) == 0)
        return -1;

    if (d[1] == 0)
    {
        /* single-limb divisor: n[1] % d[0] < d[0], so the second divide's
           quotient fits in one limb */
        q1 = n[1] / d[0];
        r1 = n[1] % d[0];
        q0 = _udiv128(r1, n[0], d[0], &r0);
        r1 = 0;
    }
    else if ((n[1] < d[1]) || ((n[1] == d[1]) && (n[0] < d[0])))
    {
        q1 = q0 = 0;
        r0 = n[0];
        r1 = n[1];
    }
    else
    {
        /* two-limb divisor: normalize so the divisor's top bit is set */
        int s = (int)_lead_zcnt64(d[1]);            /* d[1] != 0 here */
        uint64_t v1, v0, u2, u1, u0, rhat, phi, plo, t0, t1;
        mp_carry_t c;

        if (s)
        {
            v1 = (d[1] << s) | (d[0] >> (64 - s));
            v0 = d[0] << s;
            u2 = n[1] >> (64 - s);
            u1 = (n[1] << s) | (n[0] >> (64 - s));
            u0 = n[0] << s;
        }
        else
        {
            v1 = d[1]; v0 = d[0];
            u2 = 0; u1 = n[1]; u0 = n[0];
        }

        /* estimate: u2 < 2^s <= 2^63 <= v1, so this cannot overflow */
        q0 = _udiv128(u2, u1, v1, &rhat);

        /* correct while q0 * v0 > rhat:u0 (at most twice) */
        plo = _umul128(q0, v0, &phi);
        while ((phi > rhat) || ((phi == rhat) && (plo > u0)))
        {
            q0--;
            c = mp_subborrow_u64(0, plo, v0, &plo);
            phi -= c;                       /* product shrinks by v0 */
            rhat += v1;
            if (rhat < v1)                  /* rhat overflowed 2^64: done */
                break;
        }

        /* remainder = (u2:u1:u0) - q0 * (v1:v0), which now fits in 128 bits */
        plo = _umul128(q0, v0, &phi);
        t0 = q0 * v1 + phi;                 /* low limb of the high partial */
        c = mp_subborrow_u64(0, u0, plo, &r0);
        (void)mp_subborrow_u64(c, u1, t0, &t1);
        r1 = t1;

        if (s)
        {
            r0 = (r0 >> s) | (r1 << (64 - s));
            r1 >>= s;
        }
        q1 = 0;
    }

    if (q) { q[0] = q0; q[1] = q1; }
    if (r) { r[0] = r0; r[1] = r1; }
    return 0;
}

/* c = a mod b (128-bit).  b must be non-zero. */
void mod_128(const uint64_t *a, const uint64_t *b, uint64_t *c)
{
    (void)uint128_div(a, b, NULL, c);
}


/* ===========================================================================
   128-bit Montgomery setup.  No GMP: the library stands alone.
   =========================================================================== */

/* x = x * 2^128 mod n  (x into Montgomery representation, R = 2^128).
   x need not be reduced on entry.  n must be odd and non-zero. */
void to_monty128(monty128_t *mdata, uint64_t *x)
{
    uint64_t r[2];
    int i;

    /* x mod n, then 128 modular doublings.  ~1 ns per doubling; this replaces
       an mpz_init / mul_2exp / mod / clear sequence that cost more than that
       in allocator traffic alone. */
    if ((x[1] > mdata->n[1]) || ((x[1] == mdata->n[1]) && (x[0] >= mdata->n[0])))
        (void)uint128_div(x, mdata->n, NULL, r);
    else
    {
        r[0] = x[0];
        r[1] = x[1];
    }

    for (i = 0; i < 128; i++)
        dblmod128(r, mdata->n);

    x[0] = r[0];
    x[1] = r[1];
}

void monty128_init(monty128_t * mdata, uint64_t * n)
{
	// for a input modulus n, initialize constants for
	// montogomery representation
	// this assumes that n is relatively prime to 2, i.e. is odd.
	uint64_t zero_minus_n[2];
	mp_carry_t c;

	mdata->n[0] = n[0];
	mdata->n[1] = n[1];

	// negative inverse of (odd) n mod 2^64
	mdata->rho = (uint64_t)0 - multiplicative_inverse(n[0]);

	// one = R mod n = (2^128 - n) mod n: a single 128-bit reduction
	c = mp_subborrow_u64(0, 0, n[0], &zero_minus_n[0]);
	(void)mp_subborrow_u64(c, 0, n[1], &zero_minus_n[1]);
	(void)uint128_div(zero_minus_n, n, NULL, mdata->one);

	return;
}


/* ===========================================================================
   AVX-512 helpers (52-bit and 104-bit lanes) that stay out of line.
   =========================================================================== */
#if MP_USE_AVX512F && !defined(LIMB2_NO_VECTOR) && !defined(LIMB1_SCALAR_ONLY)

#if defined(__GNUC__) && !MP_USE_IFMA
// the non-IFMA path declares local copies of the dbias/vbias constants used
// by the mul52 macros, shadowing the (per-TU, otherwise unset) statics
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

void mulredc52_mask_add_vec(__m512i* c0, __mmask8 addmsk, __m512i a0, __m512i b0, __m512i n0, __m512i vrho)
{
	// CIOS modular multiplication with normal (negative) single-word nhat
	__m512i m;
	__m512i t0, t1;

#if !MP_USE_IFMA
	__m512d prod1_hd;
	__m512d prod1_ld, prod2_ld;
	// local constants: the file-scope statics of the same names in limb1.h
	// are per-translation-unit and nothing initializes this unit's copies.
	__m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
	__m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);
	__m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);
	__m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
#endif

	__m512i zero = _mm512_set1_epi64(0);

	t0 = t1 = zero;

	VEC_MUL_ACCUM_LOHI_PD(a0, b0, t0, t1);

	// m0
#if MP_USE_IFMA
	m = mul52lo(t0, vrho);
#else
	m = _mm512_and_si512(_mm512_mullo_epi64(t0, vrho), _mm512_set1_epi64(0x000fffffffffffffull));
#endif

	VEC_MUL_ACCUM_LOHI_PD(m, n0, t0, t1);

	// propagate t0's carry before we throw it away.
	// it is almost always exactly 1, but not always, for example when m = 0;
	t0 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));

	__mmask8 bmsk = _mm512_cmpge_epu64_mask(t0, n0);
	t0 = _mm512_mask_sub_epi64(t0, bmsk, t0, n0);

	// conditional addmod (double result)
	t0 = _mm512_mask_slli_epi64(t0, addmsk, t0, 1);
	bmsk = _mm512_mask_cmpgt_epu64_mask(addmsk, t0, n0);
	*c0 = _mm512_mask_sub_epi64(t0, bmsk & addmsk, t0, n0);

	return;
}

#if defined(__GNUC__) && !MP_USE_IFMA
#pragma GCC diagnostic pop
#endif

__m512i multiplicative_inverse104_x8(uint64_t* a)
{
	//    assert(a%2 == 1);  // the inverse (mod 2<<64) only exists for odd values
	__m512i x0, x1, x2, x3, x4, y, n;
	__m512i three = _mm512_set1_epi64(3), two = _mm512_set1_epi64(2), one = _mm512_set1_epi64(1);
	__m512i lo52mask = _mm512_set1_epi64(0x000fffffffffffffull);

	n = loadu64(a);
	_mm512_mullo_epi52(x0, n, three);
	x0 = _mm512_xor_epi64(x0, two);
	_mm512_mullo_epi52(y, n, x0);
	y = _mm512_sub_epi64(one, y);
	y = _mm512_and_epi64(lo52mask, y);

	x1 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x1, x0, x1);
	_mm512_mullo_epi52(y, y, y);

	x2 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x2, x1, x2);
	_mm512_mullo_epi52(y, y, y);

	x3 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x3, x2, x3);
	_mm512_mullo_epi52(y, y, y);

	x4 = _mm512_add_epi64(y, one);
	_mm512_mullo_epi52(x4, x3, x4);

	return x4;
}

void init_monty104(void)
{
#if !MP_USE_IFMA
	// sets THIS translation unit's copies of the limb1.h bias statics.  Kept
	// for source compatibility; the kernels above and in limb2.h use local
	// constants and do not depend on it.
	dbias = _mm512_castsi512_pd(set64(0x4670000000000000ULL));
	vbias1 = set64(0x4670000000000000ULL);
	vbias2 = set64(0x4670000000000001ULL);
	vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);
#endif
}

#endif /* MP_USE_AVX512F vector helpers */
