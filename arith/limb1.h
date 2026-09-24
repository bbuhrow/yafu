/*----------------------------------------------------------------------
limb1.h  --  single-limb arithmetic: 64-bit scalar + 52-bit vector

Hot inline kernels live here (static inline via MP_FORCE_INLINE); the
heavier / non-inline single-limb routines are declared here and defined in
limb1.c.  The 64-bit modular + redc family and the 52-bit vector kernels
were moved down from monty.h; the compiler/OS/ISA detection and the
_umul128 / __umulh / _addcarry_u64 / _subborrow_u64 polyfills come from
mp_platform.h, so monty.h's per-branch re-definitions of those are gone.

Linkage note: the modular block (submod/addmod/mulredc/sqrredc/...) was
plain `__inline` (external linkage) in monty.h -- a one-definition-rule
hazard across the compiler matrix.  It is now `static inline` via
MP_FORCE_INLINE, which also force-inlines (always_inline / __forceinline),
matching how microecm already treats these.  The Hurchalla functions and
the 52-bit vector kernels were already `__inline static`, so they moved
verbatim.
----------------------------------------------------------------------*/

#ifndef LIMB1_H
#define LIMB1_H

#include "mp_platform.h"

#ifdef _MSC_VER
#pragma warning(disable: 4505)   /* unreferenced static function removed */
#endif


/* ====================================================================
   64-bit scalar: Hurchalla positive-inverse redc, binary gcd, and the
   gated modular-arith family (moved from monty.h)
   ==================================================================== */

/* --- The following functions are written by Jeff Hurchalla, Copyright 2022 --- */

// Using the inline asm in this file can increase performance by ~20-25%
// (surprisingly).  Hence these macros are defined by default.
#if defined(__x86_64__) || defined(_M_X64)
#  define ALT_MULREDC_USE_INLINE_ASM_X86
#endif


// for this algorithm, see https://jeffhurchalla.com/2022/04/28/montgomery-redc-using-the-positive-inverse-mod-r/
__inline static uint64_t mulredc_pos_alt(uint64_t x, uint64_t y, uint64_t N, uint64_t invN)
{
#if defined(_MSC_VER)
	uint64_t T_hi;
	uint64_t T_lo = _umul128(x, y, &T_hi);
	uint64_t m = T_lo * invN;
	uint64_t mN_hi = __umulh(m, N);
#else
	__uint128_t prod = (__uint128_t)x * y;
	uint64_t T_hi = (uint64_t)(prod >> 64);
	uint64_t T_lo = (uint64_t)(prod);
	uint64_t m = T_lo * invN;
	__uint128_t mN = (__uint128_t)m * N;
	uint64_t mN_hi = (uint64_t)(mN >> 64);
#endif
	uint64_t tmp = T_hi + N;
#if defined(ALT_MULREDC_USE_INLINE_ASM_X86) && !defined(_MSC_VER)
	__asm__(
		"subq %[mN_hi], %[tmp] \n\t"    /* tmp = T_hi + N - mN_hi */
		"subq %[mN_hi], %[T_hi] \n\t"   /* T_hi = T_hi - mN_hi */
		"cmovaeq %[T_hi], %[tmp] \n\t"  /* tmp = (T_hi >= mN_hi) ? T_hi : tmp */
		: [tmp] "+&r"(tmp), [T_hi]"+&r"(T_hi)
		: [mN_hi] "r"(mN_hi)
		: "cc");
	uint64_t result = tmp;
#else
	tmp = tmp - mN_hi;
	uint64_t result = T_hi - mN_hi;
	result = (T_hi < mN_hi) ? tmp : result;
#endif
	return result;
}

/* Canonical single-limb binary GCD (gmplib twos-complement trick).  This one
   definition replaces arith.c's bingcd64, monty.h's bin_gcd64, and (when the
   micro files migrate) uecm_bingcd64 / upm1_bingcd64 -- all four were verified
   byte-identical.  Left as plain `static` rather than force-inlined: it has a
   loop, so inlining stays the compiler's call, matching the original monty.h. */
static UNUSED_FUNC uint64_t bingcd64(uint64_t u, uint64_t v)
{
#if 1
	if (u == 0) {
		return v;
	}
	if (v != 0) {
		int j = _trail_zcnt64(v);
		v = (uint64_t)(v >> j);
		while (1) {
			uint64_t tmp = u;
			uint64_t sub1 = (uint64_t)(v - tmp);
			uint64_t sub2 = (uint64_t)(tmp - v);
			if (tmp == v)
				break;
			u = (tmp >= v) ? v : tmp;
			v = (tmp >= v) ? sub2 : sub1;
			// For the line below, the standard way to write this algorithm
			// would have been to use _trail_zcnt64(v)  (instead of
			// _trail_zcnt64(sub1)).  However, as pointed out by
			// https://gmplib.org/manual/Binary-GCD, "in twos complement the
			// number of low zero bits on u-v is the same as v-u, so counting or
			// testing can begin on u-v without waiting for abs(u-v) to be
			// determined."  Hence we are able to use sub1 for the argument.
			// By removing the dependency on abs(u-v), the CPU can execute
			// _trail_zcnt64(void) at the same time as abs(u-v).
			j = _trail_zcnt64(sub1);
			v = (uint64_t)(v >> j);
		}
	}
	return u;
#else
	// For reference, or if in the future we need to allow an even u,
	// this version allows u to be even or odd.
	if (u == 0) {
		return v;
	}
	if (v != 0) {
		int i = _trail_zcnt64(u);
		int j = _trail_zcnt64(v);
		u = (uint64_t)(u >> i);
		v = (uint64_t)(v >> j);
		int k = (i < j) ? i : j;
		while (1) {
			uint64_t tmp = u;
			uint64_t sub1 = (uint64_t)(v - tmp);
			uint64_t sub2 = (uint64_t)(tmp - v);
			if (tmp == v)
				break;
			u = (tmp >= v) ? v : tmp;
			v = (tmp >= v) ? sub2 : sub1;
			// For the line below, the standard way to write this algorithm
			// would have been to use _trail_zcnt64(v)  (instead of
			// _trail_zcnt64(sub1)).  However, as pointed out by
			// https://gmplib.org/manual/Binary-GCD, "in twos complement the
			// number of low zero bits on u-v is the same as v-u, so counting or
			// testing can begin on u-v without waiting for abs(u-v) to be
			// determined."  Hence we are able to use sub1 for the argument.
			// By removing the dependency on abs(u-v), the CPU can execute
			// _trail_zcnt64() at the same time as abs(u-v).
			j = _trail_zcnt64(sub1);
			v = (uint64_t)(v >> j);
		}
		u = (uint64_t)(u << k);
	}
	return u;
#endif
}

/* back-compat alias: monty / ECM call sites use bin_gcd64 */
#define bin_gcd64 bingcd64

// full strength mul/sqr redc
__inline static uint64_t mulredc_pos(uint64_t x, uint64_t y, uint64_t n, uint64_t inv)
{
#if defined(__unix__) && defined(__x86_64__)
	// On Intel Skylake: 9 cycles latency, 7 fused uops.
	__uint128_t prod = (__uint128_t)x * y;
	uint64_t Thi = (uint64_t)(prod >> 64);
	uint64_t rrax = (uint64_t)(prod);
	__asm__(
		"imulq %[inv], %%rax \n\t"        /* m = T_lo * invN */
		"mulq %[n] \n\t"                  /* mN = m * N */
		"leaq (%[Thi], %[n]), %%rax \n\t" /* rax = T_hi + N */
		"subq %%rdx, %%rax \n\t"          /* rax = rax - mN_hi */
		"subq %%rdx, %[Thi] \n\t"         /* t_hi = T_hi - mN_hi */
		"cmovbq %%rax, %[Thi] \n\t"       /* t_hi = (T_hi<mN_hi) ? rax : t_hi */
		: [Thi] "+&bcSD"(Thi), "+&a"(rrax)
		: [n] "r"(n), [inv]"r"(inv)
		: "rdx", "cc");
	uint64_t result = Thi;
	return result;

#else
	return mulredc_pos_alt(x, y, n, inv);
#endif
}
__inline static uint64_t sqrredc_pos(uint64_t x, uint64_t n, uint64_t inv)
{
#if defined(__unix__) && defined(__x86_64__)
	// On Intel Skylake: 9 cycles latency, 7 fused uops.
	__uint128_t prod = (__uint128_t)x * x;
	uint64_t Thi = (uint64_t)(prod >> 64);
	uint64_t rrax = (uint64_t)(prod);
	__asm__(
		"imulq %[inv], %%rax \n\t"        /* m = T_lo * invN */
		"mulq %[n] \n\t"                  /* mN = m * N */
		"leaq (%[Thi], %[n]), %%rax \n\t" /* rax = T_hi + N */
		"subq %%rdx, %%rax \n\t"          /* rax = rax - mN_hi */
		"subq %%rdx, %[Thi] \n\t"         /* t_hi = T_hi - mN_hi */
		"cmovbq %%rax, %[Thi] \n\t"       /* t_hi = (T_hi<mN_hi) ? rax : t_hi */
		: [Thi] "+&bcSD"(Thi), "+&a"(rrax)
		: [n] "r"(n), [inv]"r"(inv)
		: "rdx", "cc");
	uint64_t result = Thi;
	return result;

#else
	return mulredc_pos_alt(x, x, n, inv);
#endif
}
__inline static uint64_t mfma64(uint64_t x, uint64_t y, uint64_t c, uint64_t N, uint64_t invN)
{
#if defined(_MSC_VER)
	uint64_t T_hi;
	uint64_t T_lo = _umul128(x, y, &T_hi);
	uint64_t m = T_lo * invN;
	uint64_t mN_hi = __umulh(m, N);
#else
	__uint128_t z = (__uint128_t)x * y;
	uint64_t u = (uint64_t)(z >> 64);
	uint64_t v = (uint64_t)z;
	uint64_t w = (u < N - c) ? u + c : u + c - N;  // modular add
	uint64_t T_hi = w;
	uint64_t T_lo = v;

	uint64_t m = T_lo * invN;
	__uint128_t mN = (__uint128_t)m * N;
	uint64_t mN_hi = (uint64_t)(mN >> 64);
	uint64_t tmp = T_hi + N;
	tmp = tmp - mN_hi;
	uint64_t result = T_hi - mN_hi;
	result = (T_hi < mN_hi) ? tmp : result;
	return result;

#endif

}

/* --- end Hurchalla functions --- */


/********************* 64-bit modular arith **********************/

#if (defined(GCC_ASM64X) || defined(__MINGW64__)) && !defined(ASM_ARITH_DEBUG)

MP_FORCE_INLINE uint64_t submod(uint64_t a, uint64_t b, uint64_t n)
{
    __asm__(
        "xorq %%r8, %%r8 \n\t"
        "subq %1, %0 \n\t"
        "cmovc %2, %%r8 \n\t"
        "addq %%r8, %0 \n\t"
        : "+r"(a)
        : "r"(b), "r"(n)
        : "r8", "cc");

    return a;
}

MP_FORCE_INLINE uint64_t addmod(uint64_t x, uint64_t y, uint64_t n)
{
    uint64_t t = x - n;
    x += y;
    __asm__("add %2, %1\n\t"
        "cmovc %1, %0\n\t"
        :"+r" (x), "+&r" (t)
        : "r" (y)
        : "cc"
        );
    return x;
}

#if defined(USE_AVX512F) || defined(USE_BMI2)

MP_FORCE_INLINE uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulx %2, %%r10, %%r11	\n\t"
            "movq %%r10, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "xorq %%r12, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=&d"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulx %2, %%r10, %%r11	\n\t"
            "movq %3, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "mulq %%r10 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=d"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "cc");

    }
    return x;
}

MP_FORCE_INLINE uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulx %1, %%r10, %%r11	\n\t"
            "movq %%r10, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "xorq %%r12, %%r12 \n\t"
            "mulq %2 \n\t"
            "mulq %3 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %3, %%r12 \n\t"
            "subq %3, %%rdx \n\t"
            "cmovc %%r12, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=&d"(x)
            : "0"(x), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulx %1, %%r10, %%r11	\n\t"
            "movq %2, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "mulq %%r10 \n\t"
            "mulq %3 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "subq %3, %%rdx \n\t"
            "cmovc %3, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=d"(x)
            : "0"(x), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "cc");

    }
    return x;
}

MP_FORCE_INLINE uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulx %2, %%r10, %%r11	\n\t"
        "movq %3, %%rax		\n\t"
        "xorq %%r8, %%r8 \n\t"
        "mulq %%r10 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%r8 \n\t"
        "addq %%r8, %%rdx \n\t"
        : "=d"(x)
        : "0"(x), "r"(y), "r"(nhat), "r"(n)
        : "rax", "r8", "r10", "r11", "cc");

    return x;
}

MP_FORCE_INLINE uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulx %1, %%r10, %%r11	\n\t"
        "movq %2, %%rax		\n\t"
        "xorq %%r8, %%r8 \n\t"
        "mulq %%r10 \n\t"
        "mulq %3 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "subq %3, %%rdx \n\t"
        "cmovc %3, %%r8 \n\t"
        "addq %%r8, %%rdx \n\t"
        : "=d"(x)
        : "0"(x), "r"(nhat), "r"(n)
        : "rax", "r8", "r10", "r11", "cc");

    return x;
}


#else // !USE_AVX512F && !USE_BMI2


MP_FORCE_INLINE uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "movq $0, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "xorq %%rax, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "movq $0, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "cc");

    }
    return x;
}

MP_FORCE_INLINE uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulq %2	\n\t"
        "movq %%rax, %%r10		\n\t"
        "movq %%rdx, %%r11		\n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(x)
        : "0"(x), "r"(y), "r"(nhat), "r"(n)
        : "rdx", "r10", "r11", "cc");

    return x;
}

MP_FORCE_INLINE uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "movq $0, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "xorq %%rax, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(x), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "movq $0, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(x), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "cc");

    }
    return x;
}

MP_FORCE_INLINE uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulq %2	\n\t"
        "movq %%rax, %%r10		\n\t"
        "movq %%rdx, %%r11		\n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(x)
        : "0"(x), "r"(x), "r"(nhat), "r"(n)
        : "rdx", "r10", "r11", "cc");

    return x;
}

#endif // !USE_AVX512F && !USE_BMI2


#else // (!GCC_ASM64X && !__MINGW64__) || ASM_ARITH_DEBUG

// TODO: need something portable to replace 64-bit assembler versions
// of modular multiplication.  This is getting closer, but for now these things
// are spread out locally where they are needed instead of being gathered here,
// apologies for the uglyness.


MP_FORCE_INLINE uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

MP_FORCE_INLINE uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    return x;
}

MP_FORCE_INLINE uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

MP_FORCE_INLINE uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    return x;
}

MP_FORCE_INLINE uint64_t submod(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r0;
    if (_subborrow_u64(0, a, b, &r0))
        r0 += n;
    return r0;
}

MP_FORCE_INLINE uint32_t submod32(uint32_t a, uint32_t b, uint32_t n)
{
    uint32_t r0;
    if (_subborrow_u32(0, a, b, &r0))
        r0 += n;
    return r0;
}

MP_FORCE_INLINE uint64_t addmod(uint64_t x, uint64_t y, uint64_t n)
{
#if 0
    uint64_t r;
    uint64_t tmp = x - n;
    uint8_t c = _addcarry_u64(0, tmp, y, &r);
    return (c) ? r : x + y;
#else
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x>=n-y) ? x-(n-y) : x+y;
#endif
}

MP_FORCE_INLINE uint32_t addmod32(uint32_t x, uint32_t y, uint32_t n)
{
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x >= n - y) ? x - (n - y) : x + y;
}



// good to 60 bit inputs
MP_FORCE_INLINE uint64_t sqrredc60(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    uint8_t c;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// good to 60 bit inputs
MP_FORCE_INLINE uint64_t mulredc60(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    uint8_t c;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// this works if inputs are 62 bits or less
#define addmod60(x, y, n) ((x) + (y))

#endif // (!GCC_ASM64X && !__MINGW64__) || ASM_ARITH_DEBUG


/* ====================================================================
   52-bit vector single-limb (AVX-512).  These primitives (carryprop,
   mul52*, VEC_MUL_ACCUM_LOHI_PD, the _mm512_*_epi52 carry ops, the
   and64/store64/... aliases, and vec_u104_t) are shared with the 104-bit
   double-limb code, so limb2.h includes limb1.h to reuse them.  In monty.h
   this lived inside one big `#ifdef USE_AVX512F` that also wrapped the
   104-bit code; here limb1.h closes it and limb2.h reopens it.
   ==================================================================== */
/* A consumer that has its own 52-bit vector code (microecm / micropm1) can
   `#define LIMB1_SCALAR_ONLY` before including this header to take just the
   scalar 64-bit layer and skip the vector primitives below -- avoiding a
   collision on the shared dbias/vbias globals. */
#if defined(USE_AVX512F) && !defined(LIMB1_SCALAR_ONLY)

#define and64 _mm512_and_si512
#define store64 _mm512_store_si512
#define storeu64 _mm512_storeu_si512
#define mstoreu64 _mm512_mask_storeu_si512
#define storeu512 _mm512_storeu_si512
#define add64 _mm512_add_epi64
#define sub64 _mm512_sub_epi64
#define set64 _mm512_set1_epi64
#define srli64 _mm512_srli_epi64
#define load64 _mm512_load_si512
#define loadu64 _mm512_loadu_si512
#define loadu512 _mm512_loadu_si512
#define castpd _mm512_castsi512_pd
#define castepu _mm512_castpd_si512

typedef struct
{
    //ALIGNED_MEM
    uint64_t data[2][8];
} vec_u104_t;

/********************* 52-bit vector Montgomery arith **********************/
#define carryprop(lo, hi, mask) \
	{ __m512i carry = _mm512_srli_epi64(lo, 52);	\
	hi = _mm512_add_epi64(hi, carry);		\
	lo = _mm512_and_epi64(mask, lo); }

#if defined(INTEL_COMPILER) || defined(INTEL_LLVM_COMPILER)
#define ROUNDING_MODE (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)
#else
#define ROUNDING_MODE _MM_FROUND_CUR_DIRECTION
#endif


#ifdef IFMA

#define FORCE_INLINE __inline

#define mul52hi(b, c) \
	_mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b)


#define mul52lo(b, c) \
	_mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b)

FORCE_INLINE static void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
	*l = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), c, b);
	*h = _mm512_madd52hi_epu64(_mm512_set1_epi64(0), c, b);
	return;
}

#else

static __m512d dbias;
static __m512i vbias1;
static __m512i vbias2;
static __m512i vbias3;

#define mul52lo(b, c) \
	_mm512_and_si512(_mm512_mullo_epi64(b, c), _mm512_set1_epi64(0x000fffffffffffffull))

__inline static __m512i mul52hi(__m512i b, __m512i c)
{
	__m512d prod1_ld = _mm512_cvtepu64_pd(b);
	__m512d prod2_ld = _mm512_cvtepu64_pd(c);
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	return _mm512_sub_epi64(castepu(prod1_ld), vbias1);
}
__inline static void mul52lohi(__m512i b, __m512i c, __m512i* l, __m512i* h)
{
	__m512d prod1_ld = _mm512_cvtepu64_pd(b);
	__m512d prod2_ld = _mm512_cvtepu64_pd(c);
	__m512d prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	*h = _mm512_sub_epi64(castepu(prod1_hd), vbias1);
	prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
	*l = _mm512_castpd_si512(prod1_ld);
	*l = _mm512_and_si512(*l, _mm512_set1_epi64(0x000fffffffffffffull));
	*h = _mm512_and_si512(*h, _mm512_set1_epi64(0x000fffffffffffffull));
	return;
}

#endif



#ifdef IFMA
#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_madd52lo_epu64(_mm512_set1_epi64(0), a, b);

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
    lo = _mm512_madd52lo_epu64(lo, a, b); \
    hi = _mm512_madd52hi_epu64(hi, a, b);
#else

#define VEC_MUL_ACCUM_LOHI_PD(a, b, lo, hi) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi = _mm512_add_epi64(hi, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo = _mm512_add_epi64(lo, _mm512_sub_epi64(castepu(prod1_ld), vbias3));

#define VEC_MUL2_ACCUM_LOHI_PD(c, a, b, lo1, hi1, lo2, hi2) \
	prod1_ld = _mm512_cvtepu64_pd(a);		\
	prod2_ld = _mm512_cvtepu64_pd(b);		\
	prod3_ld = _mm512_cvtepu64_pd(c);		\
    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
    hi1 = _mm512_add_epi64(hi1, _mm512_sub_epi64(castepu(prod1_hd), vbias1)); \
	hi2 = _mm512_add_epi64(hi2, _mm512_sub_epi64(castepu(prod2_hd), vbias1)); \
    prod1_hd = _mm512_sub_pd(castpd(vbias2), prod1_hd); \
	prod2_hd = _mm512_sub_pd(castpd(vbias2), prod2_hd); \
	prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); \
	lo1 = _mm512_add_epi64(lo1, _mm512_sub_epi64(castepu(prod1_ld), vbias3)); \
	lo2 = _mm512_add_epi64(lo2, _mm512_sub_epi64(castepu(prod2_ld), vbias3));

#define _mm512_mullo_epi52(c, a, b) \
    c = _mm512_and_si512(_mm512_mullo_epi64(a, b), _mm512_set1_epi64(0x000fffffffffffffull));
#endif


// better to #define these?  Or make them static inline here?
//__m512i _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_addsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_subsetc_epi52(__m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_subsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout);
//__m512i _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout);
//__m512i _mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout);
//__m512i _mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout);

// static inline works and is quite a bit faster for tinyecm
__inline static __m512i _mm512_addsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_mask_addsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_subsetc_epi52(__m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_mask_subsetc_epi52(__m512i c, __mmask8 mask, __m512i a, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(mask, b, a);
    t = _mm512_mask_and_epi64(c, mask, t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_adc_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_mask_adc_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, b);
    t = _mm512_mask_add_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpgt_epu64_mask(t, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_sbb_epi52(__m512i a, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, b);
    *cout = _mm512_cmpgt_epu64_mask(b, a);
    __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(c, 1));
    *cout = (__mmask8)_mm512_kor(*cout, _mm512_cmpgt_epu64_mask(t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}
__inline static __m512i _mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
{
    __m512i t = _mm512_mask_sub_epi64(a, m, a, b);
    *cout = _mm512_mask_cmpgt_epu64_mask(m, b, a);
    __m512i t2 = _mm512_mask_sub_epi64(a, m, t, _mm512_maskz_set1_epi64(c, 1));
    *cout = (__mmask8)_mm512_kor(*cout, _mm512_mask_cmpgt_epu64_mask(m, t2, t));
    t2 = _mm512_and_epi64(t2, _mm512_set1_epi64(0xfffffffffffffULL));
    return t2;
}
__inline static __m512i _mm512_addcarry_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_add_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = c & _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0xfffffffffffffULL));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}
__inline static __m512i _mm512_subborrow_epi52(__m512i a, __mmask8 c, __mmask8* cout)
{
    __m512i t = _mm512_sub_epi64(a, _mm512_maskz_set1_epi64(c, 1));
    *cout = _mm512_cmpeq_epu64_mask(a, _mm512_set1_epi64(0));
    t = _mm512_and_epi64(t, _mm512_set1_epi64(0xfffffffffffffULL));
    return t;
}

void mulredc52_mask_add_vec(__m512i* c0, __mmask8 addmsk, __m512i a0, __m512i b0, __m512i n0, __m512i vrho);

#endif /* USE_AVX512F (52-bit vector) */

/* ====================================================================
   Small inline helpers
   ==================================================================== */

/* Montgomery inverse via Newton iteration (Dumas): 64-bit inverse mod 2^64.
   Canonical inline definition -- verified byte-identical to monty.c's
   multiplicative_inverse and the uecm_/upm1_ copies, which are retired as those
   files migrate (monty.c -> limb2.c; the micro files in their own pass). */
MP_FORCE_INLINE uint64_t multiplicative_inverse(uint64_t a)
{
    uint64_t x0 = (3 * a) ^ 2, y = 1 - a * x0;
    uint64_t x1 = x0 * (1 + y); y *= y;
    uint64_t x2 = x1 * (1 + y); y *= y;
    uint64_t x3 = x2 * (1 + y); y *= y;
    return x3 * (1 + y);
}

/* count leading / trailing zeros: the zero-safe 64-bit pair is now
   _lead_full_zcnt / _trail_full_zcnt in mp_platform.h, and the 52-bit pair is
   _lead_full_zcnt52 / _trail_full_zcnt52 above.  Nothing to declare here. */


/* ====================================================================
   Heavier single-limb routines      [DECLARATIONS -> limb1.c]
   (consolidated from arith.c / arith.h)
   ==================================================================== */

/* wide sp* single-word ops (one asm and one portable definition in limb1.c) */
void     spAdd(uint64_t u, uint64_t v, uint64_t *sum, uint64_t *carry);
void     spAdd3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sum, uint64_t *carry);
void     spSub(uint64_t u, uint64_t v, uint64_t *sub, uint64_t *borrow);
void     spSub3(uint64_t u, uint64_t v, uint64_t w, uint64_t *sub, uint64_t *borrow);
void     spMultiply(uint64_t u, uint64_t v, uint64_t *product, uint64_t *carry);
void     spMulAdd(uint64_t u, uint64_t v, uint64_t w, uint64_t t,
                  uint64_t *lower, uint64_t *carry);
void     spMulMod(uint64_t u, uint64_t v, uint64_t m, uint64_t *w);
uint64_t spDivide(uint64_t *q, uint64_t *r, uint64_t u[2], uint64_t v);
uint64_t u64div(uint64_t c, uint64_t n);
uint64_t spPRP2(uint64_t p);                    /* asm branch only -- see note in limb1.c */

/* bit / digit counts */
uint64_t spBits(uint64_t n);
int      bits64(uint64_t n);
int      ndigits_1(uint64_t n);

/* modular exponentiation / inverse / sqrt / jacobi */
void     spModExp(uint64_t a, uint64_t b, uint64_t m, uint64_t *u);
uint32_t modinv_1(uint32_t a, uint32_t p);
uint32_t modinv_1b(uint32_t a, uint32_t p);
uint32_t modinv_1c(uint32_t a, uint32_t p);
int      jacobi_1(uint64_t n, uint64_t p);
void     ShanksTonelli_1(uint64_t a, uint64_t p, uint64_t *sq);

/* scalar GCDs.  bingcd64 is now a single inline definition above (canonical);
   arith.c's and monty.h's copies are retired.  The uecm_/upm1_ copies come out
   during the micro-file migration.  The spGCD / gcd64 / spBinGCD family below
   stay external in limb1.c. */
/* ---- zero-safe bit scan, 52-bit field ---------------------------------
   The 64-bit pair lives in mp_platform.h as _lead_full_zcnt/_trail_full_zcnt.
   These are the 52-bit-width equivalents of monty.c's my_clz52/my_ctz52, and
   they reproduce that behaviour exactly: for a non-zero argument both return
   the full 64-bit scan, and only the zero case returns 52.

   That is worth reading twice for the leading case.  _lead_full_zcnt52 does
   NOT return a count relative to a 52-bit field -- it does not subtract the
   12 high padding bits -- so for n != 0 it is identical to _lead_full_zcnt.
   That is monty.c's long-standing behaviour, preserved verbatim here rather
   than silently "corrected", since my_clz104 is built on it and changing it
   would shift that result by 12.  If a true 52-bit-relative count was the
   intent, this is the place to fix it. */
MP_FORCE_INLINE uint64_t _lead_full_zcnt52(uint64_t n) {
    if (n)
        return (uint64_t)_lead_zcnt64(n);
    return 52;
}
MP_FORCE_INLINE uint64_t _trail_full_zcnt52(uint64_t n) {
    if (n)
        return (uint64_t)_trail_zcnt64(n);
    return 52;
}

uint64_t spGCD(uint64_t x, uint64_t y);
uint64_t gcd64(uint64_t x, uint64_t y);
uint64_t spBinGCD(uint64_t u, uint64_t v);
uint64_t spBinGCD_odd(uint64_t u, uint64_t v);
void     dblGCD(double x, double y, double *w);


#endif /* LIMB1_H */
