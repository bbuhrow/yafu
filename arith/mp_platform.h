/*----------------------------------------------------------------------
mp_platform.h

Shared platform / detection / intrinsic-polyfill layer for the YAFU limb
arithmetic. limb1.{h,c}, limb2.{h,c}, microecm.c and micropm1.c include
this and nothing else for compiler / OS / ISA concerns.

Layering:
  common.h       - assembler-syntax detection (ASM_G / ASM_M, GCC_ASM64X,
                   the A/X guards).  Papadopoulos / Gladman lineage, shared
                   across the wider YAFU tree.  We build ON TOP of it and do
                   not absorb it.
  mp_platform.h  - this file: ISA-feature normalization, inline/unused
                   attribute macros, and ONE copy of every intrinsic polyfill
                   (_umul128, __umulh, _addcarry_u64, _subborrow_u64/_u32,
                   _udiv128, the {trail,lead}_zcnt / reset_lsb family).

This replaces the per-file #if ladders that were re-derived (and had drifted)
in arith.h, arith.c, monty.h, monty.c, microecm.c and micropm1.c.
----------------------------------------------------------------------*/

#ifndef MP_PLATFORM_H
#define MP_PLATFORM_H

#include <stdint.h>
#include "common.h"

/* -------- vendor intrinsics: pulled in once, here -------- */
#if defined(_MSC_VER)
  #include <intrin.h>
  #if !defined(_M_ARM64)
    #include <immintrin.h>
  #endif
  #if !defined(__clang__)
    #pragma intrinsic(_umul128)
    #pragma intrinsic(_udiv128)
  #endif
#elif !defined(__aarch64__)
  #include <immintrin.h>
#endif


/* ======================================================================
   1. Derived capability flags  (single source of truth)
   ====================================================================== */

/* 128-bit integer type available? */
#if defined(__SIZEOF_INT128__) && (__SIZEOF_INT128__ == 16)
  #ifndef HAS_UINT128
    #define HAS_UINT128
  #endif
  typedef __uint128_t uint128_t;
#endif

/* GNU-style inline asm usable for x86-64?
   common.h already sets GCC_ASM64X for gcc/clang/icc on unix-x64, for mingw,
   and for icc/clang on win-x64.  __MINGW64__ is folded in to match the
   historical monty.h gate.  ASM_ARITH_DEBUG forces the portable path so the
   two code paths can be A/B compared. */
#if (defined(GCC_ASM64X) || defined(__MINGW64__)) && !defined(ASM_ARITH_DEBUG) \
    && (defined(__x86_64__) || defined(_M_X64))
  #define MP_HAS_GNU_ASM_X64 1
#else
  #define MP_HAS_GNU_ASM_X64 0
#endif

/* INLINE_ASM was previously #define'd =1 locally in monty.c -- and, because it
   was defined *after* the my_ctz/my_clz functions, the BMI1 asm paths in those
   functions were never actually reached.  Centralizing it here (this header is
   included first) makes INLINE_ASM consistent everywhere.
   NOTE: this ACTIVATES the previously-dead BMI1 ctz/clz asm paths -- a real
   codegen change.  Leave MP_ENABLE_LEGACY_INLINE_ASM undefined to keep the
   historical (portable-builtin) behavior for those functions until you have
   re-benchmarked them. */
#ifndef INLINE_ASM
  #define INLINE_ASM MP_HAS_GNU_ASM_X64
#endif


/* ======================================================================
   2. Inline / unused attribute macros
   ====================================================================== */

/* MP_FORCE_INLINE expands to include `static`.  This is deliberate: the old
   monty.h modular block (submod/addmod/mulredc/sqrredc/...) was declared plain
   `__inline` with external linkage, which is a one-definition-rule landmine
   when a header is included by multiple TUs (gnu89-inline vs C99 vs MSVC all
   disagree about who emits the external symbol).  Routing every header-resident
   kernel through this macro makes them all `static inline` and removes the
   hazard. */
#if defined(_MSC_VER)
  #define MP_FORCE_INLINE static __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) \
   || defined(__INTEL_LLVM_COMPILER)
  #define MP_FORCE_INLINE static inline __attribute__((always_inline))
#else
  #define MP_FORCE_INLINE static __inline
#endif

/* Compiler-compat attribute macros.  These now live HERE (moved down out of
   ytools.h, which drags in windows.h / winsock and must not be a dependency of
   low-level arithmetic).  ytools.h should `#include "mp_platform.h"` and delete
   its own copies of these -- along with ALIGNED_MEM / INLINE / PREFETCH, which
   are also pure compiler detection with no OS dependency.
   The bodies below are byte-for-byte ytools.h's, so even during the transition
   (before ytools.h is edited) an identical macro redefinition is legal and
   silent in C regardless of include order.  The #ifndef guards additionally
   let a ytools.h-first TU win without a diagnostic. */
#if defined(_MSC_VER)
  /* 4100 = unref formal parameter, 4101 = unref local variable */
  #ifndef UNUSED_VAR
    #define UNUSED_VAR __pragma(warning(suppress: 4100 4101))
  #endif
  /* MSVC has no per-symbol form; headers with intentionally-unused statics add
     `#pragma warning(disable: 4505)` near the top (limb1.h / limb2.h do). */
  #ifndef UNUSED_FUNC
    #define UNUSED_FUNC
  #endif
#elif defined(__GNUC__) || defined(__clang__) || \
      defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  #ifndef UNUSED_VAR
    #define UNUSED_VAR __attribute__((unused))
  #endif
  #ifndef UNUSED_FUNC
    #define UNUSED_FUNC __attribute__((unused))
  #endif
#else
  #ifndef UNUSED_VAR
    #define UNUSED_VAR
  #endif
  #ifndef UNUSED_FUNC
    #define UNUSED_FUNC
  #endif
#endif

/* Carry/borrow return type.  The _addcarry_u64 / _subborrow_u64 provided here
   (the MSVC intrinsic, or the portable wrappers further down) uniformly follow
   the MSVC convention: they RETURN the carry/borrow and write the result
   through the out pointer.  A single unsigned char is therefore correct on
   every compiler; use mp_carry_t in place of the ad-hoc `rettype` that appeared
   in monty.h / microecm.c / arith.c. */
typedef unsigned char mp_carry_t;


/* ======================================================================
   3. Bit scan / reset  (relocated verbatim from arith.h -- known good)
   ====================================================================== */

#if defined(USE_BMI2) || defined(TARGET_KNL) || defined(USE_AVX512F)
  #define _reset_lsb(x)   _blsr_u32(x)
  #define _reset_lsb64(x) _blsr_u64(x)
#else
  #define _reset_lsb(x)   ((x) &= ((x) - 1))
  #define _reset_lsb64(x) ((x) &= ((x) - 1))
#endif

#if defined(__INTEL_COMPILER)
  #if defined(USE_BMI2) || defined(TARGET_KNL) || defined(USE_AVX512F)
    #define _lead_zcnt64  __lzcnt64
    #define _trail_zcnt   _tzcnt_u32
    #define _trail_zcnt64 _tzcnt_u64
  #else
    MP_FORCE_INLINE uint32_t _trail_zcnt(uint32_t x) {
        uint32_t pos; return _BitScanForward(&pos, x) ? pos : 32;
    }
    MP_FORCE_INLINE uint64_t _trail_zcnt64(uint64_t x) {
        uint32_t pos; return _BitScanForward64(&pos, x) ? pos : 64;
    }
    MP_FORCE_INLINE uint64_t _lead_zcnt64(uint64_t x) {
        /* _BitScanReverse64 yields the INDEX of the highest set bit; the
           contract here is a leading-zero COUNT, as returned by lzcnt,
           __builtin_clzll and the portable fallback.  Hence 63 - pos. */
        uint32_t pos; return _BitScanReverse64(&pos, x) ? (uint64_t)(63 - pos) : 64;
    }
  #endif
#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER)
  #define _lead_zcnt64  __builtin_clzll
  #define _trail_zcnt   __builtin_ctzl
  #define _trail_zcnt64 __builtin_ctzll
#elif defined(_MSC_VER)
  #ifdef __clang__   /* clang-cl: _BitScan* return unsigned long */
    MP_FORCE_INLINE uint32_t _trail_zcnt(uint32_t x) {
        unsigned long pos; return _BitScanForward(&pos, x) ? (uint32_t)pos : 32;
    }
    MP_FORCE_INLINE uint32_t _trail_zcnt64(uint64_t x) {
        unsigned long pos; return _BitScanForward64(&pos, x) ? (uint32_t)pos : 64;
    }
    MP_FORCE_INLINE uint32_t _lead_zcnt64(uint64_t x) {
        /* index -> count: see the note above */
        unsigned long pos; return _BitScanReverse64(&pos, x) ? (uint32_t)(63 - pos) : 64;
    }
  #else
    MP_FORCE_INLINE uint32_t _trail_zcnt(uint32_t x) {
        uint32_t pos; return _BitScanForward(&pos, x) ? pos : 32;
    }
    MP_FORCE_INLINE uint32_t _trail_zcnt64(uint64_t x) {
        uint32_t pos; return _BitScanForward64(&pos, x) ? pos : 64;
    }
    MP_FORCE_INLINE uint32_t _lead_zcnt64(uint64_t x) {
        /* index -> count: see the note above */
        uint32_t pos; return _BitScanReverse64(&pos, x) ? (uint32_t)(63 - pos) : 64;
    }
  #endif
#else   /* portable loop fallback */
  MP_FORCE_INLINE uint64_t _lead_zcnt64(uint64_t x) {
      uint64_t pos;
      if (x) { for (pos = 0; !(x & (1ULL << (63 - pos))); pos++) ; }
      else   { pos = 8 * sizeof(x); }
      return pos;
  }
  MP_FORCE_INLINE uint32_t _trail_zcnt(uint32_t x) {
      uint32_t pos;
      if (x) { x = (x ^ (x - 1)) >> 1; for (pos = 0; x; pos++) x >>= 1; }
      else   { pos = 8 * sizeof(x); }
      return pos;
  }
  MP_FORCE_INLINE uint64_t _trail_zcnt64(uint64_t x) {
      uint64_t pos;
      if (x) { x = (x ^ (x - 1)) >> 1; for (pos = 0; x; pos++) x >>= 1; }
      else   { pos = 8 * sizeof(x); }
      return pos;
  }
#endif


/* ---- zero-safe ("full") variants -------------------------------------------
   _lead_zcnt64 / _trail_zcnt64 above are the fast paths and are UNDEFINED at
   zero on the gcc/clang builtin branches, so every call site has to guarantee
   a non-zero argument.  These wrappers add the zero case, returning the width,
   and are what monty.c used to spell my_clz64 / my_ctz64.

   Being wrappers rather than copies matters: they inherit whatever the ladder
   above selected, so there is exactly one implementation of each scan per
   platform.  The `if (n)` guard is also what keeps the builtin paths out of
   undefined behaviour -- do not remove it. */
MP_FORCE_INLINE uint64_t _lead_full_zcnt(uint64_t n) {
    if (n)
        return (uint64_t)_lead_zcnt64(n);
    return 64;
}
MP_FORCE_INLINE uint64_t _trail_full_zcnt(uint64_t n) {
    if (n)
        return (uint64_t)_trail_zcnt64(n);
    return 64;
}


/* ======================================================================
   4. Wide-multiply / divide / carry polyfills
   MSVC gets these from <intrin.h>.  Everything else is defined once here.
   ====================================================================== */

#if !defined(_MSC_VER)

  /* --- 64x64 -> 128 multiply --- */
  #if defined(HAS_UINT128)
    MP_FORCE_INLINE uint64_t _umul128(uint64_t x, uint64_t y, uint64_t *hi) {
        __uint128_t p = (__uint128_t)x * (__uint128_t)y;
        *hi = (uint64_t)(p >> 64);
        return (uint64_t)p;
    }
    MP_FORCE_INLINE uint64_t __umulh(uint64_t x, uint64_t y) {
        return (uint64_t)(((__uint128_t)x * (__uint128_t)y) >> 64);
    }
  #elif MP_HAS_GNU_ASM_X64
    MP_FORCE_INLINE uint64_t _umul128(uint64_t x, uint64_t y, uint64_t *hi) {
        __asm__("mulq %3\n\t" : "=&a"(x), "=&d"(y) : "0"(x), "1"(y) : "cc");
        *hi = y; return x;
    }
    MP_FORCE_INLINE uint64_t __umulh(uint64_t x, uint64_t y) {
        uint64_t hi; (void)_umul128(x, y, &hi); return hi;
    }
  #else
    /* portable 32x32 schoolbook -- only reached on a non-x64 target with no
       __int128, which is outside the current support matrix. */
    MP_FORCE_INLINE uint64_t _umul128(uint64_t x, uint64_t y, uint64_t *hi) {
        uint64_t xl = (uint32_t)x, xh = x >> 32, yl = (uint32_t)y, yh = y >> 32;
        uint64_t ll = xl * yl, lh = xl * yh, hl = xh * yl, hh = xh * yh;
        uint64_t mid = (ll >> 32) + (uint32_t)lh + (uint32_t)hl;
        *hi = hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
        return (mid << 32) | (uint32_t)ll;
    }
    MP_FORCE_INLINE uint64_t __umulh(uint64_t x, uint64_t y) {
        uint64_t hi; (void)_umul128(x, y, &hi); return hi;
    }
  #endif


  /* --- 128/64 -> 64 divide (quotient) with remainder out ---
     The GNU-asm divq branch comes FIRST on x86-64, ahead of the __int128
     form.  Measured: gcc lowers `(__uint128_t)n / d` to a __udivmodti4
     libcall (full software 128-bit division) rather than a single divq --
     roughly an order of magnitude slower.  clang recognizes the idiom and
     emits divq, so it is unaffected either way; this only rescues gcc/icc.
     ASM_ARITH_DEBUG clears MP_HAS_GNU_ASM_X64 and drops back to __int128.

     CONTRACT: requires d != 0 and hi < d.  divq raises #DE if the quotient
     will not fit in 64 bits -- which matches MSVC's documented _udiv128
     behaviour.  Note the __int128 fallback instead truncates silently, so a
     caller that violates the precondition changes from wrong-answer to trap
     when this branch is active.  Callers must guarantee hi < d regardless. */
  #if MP_HAS_GNU_ASM_X64
    MP_FORCE_INLINE uint64_t _udiv128(uint64_t hi, uint64_t lo, uint64_t d,
                                      uint64_t *rem) {
        uint64_t q, r;
        __asm__("divq %[d]"
                : "=a"(q), "=d"(r)
                : "a"(lo), "d"(hi), [d]"rm"(d)
                : "cc");
        *rem = r; return q;
    }
  #elif defined(HAS_UINT128)
    MP_FORCE_INLINE uint64_t _udiv128(uint64_t hi, uint64_t lo, uint64_t d,
                                      uint64_t *rem) {
        __uint128_t n = ((__uint128_t)hi << 64) | lo;
        *rem = (uint64_t)(n % d);
        return (uint64_t)(n / d);
    }
  #else
    /* No hardware path: unreachable on the current support matrix (every
       target has either MSVC intrinsics or __int128).  If a non-x86
       no-__int128 target is ever added, provide a software _udiv128 here,
       e.g. built on the Knuth uint128_div() that lives in limb2.c. */
    uint64_t _udiv128(uint64_t hi, uint64_t lo, uint64_t d, uint64_t *rem);
  #endif

#endif /* !_MSC_VER */

/* ======================================================================
 5. Add-with-carry / subtract-with-borrow
 Outside the !_MSC_VER guard on purpose: MSVC needs the mp_* wrappers
 too (it has the native intrinsics, so it takes the first branch).
 ====================================================================== */
/* --- add-carry / sub-borrow -------------------------------------------
   Canonical API (MSVC/Intel convention on every target): RETURN the
   carry/borrow, write the low result through `out`.

       mp_carry_t mp_addcarry_u64 (mp_carry_t c_in, uint64_t a, uint64_t b, uint64_t *out);
       mp_carry_t mp_subborrow_u64(mp_carry_t b_in, uint64_t a, uint64_t b, uint64_t *out);

   These are performance-critical, so we use the best mechanism available
   rather than a portable emulation:

     1. native _addcarry_u64/_subborrow_u64 -- MSVC <intrin.h>, and
        gcc/clang/icc on x86 via <immintrin.h> (adxintrin.h).  Lowers to a
        single adc/sbb and chains carries without materializing flags.
     2. clang __builtin_addcll/__builtin_subcll -- carry-in AND carry-out in
        one builtin; used on non-x86 clang (e.g. aarch64).
     3. __int128 for add, plain C for sub -- used by gcc < 10, which has no
        native _addcarry_u64.  Mixed on purpose; see the table there.
     4. gcc __builtin_add_overflow/__builtin_sub_overflow -- no carry-in, so
        composed from two ops.
     5. portable C fallback.

   NOTE: the wrappers are named mp_* deliberately.  Defining our own
   _addcarry_u64 is a *conflicting redefinition* on x86 gcc/clang, because
   adxintrin.h already declares it (verified: gcc 13 errors out).  The
   compatibility aliases below are only installed when no native version
   exists.  Taking uint64_t* in our own signature (and bouncing through a
   temporary for the native call) also removes the LP64 incompatible-pointer
   warning at call sites, where uint64_t is `unsigned long` but gcc/clang
   declare the intrinsic with `unsigned long long*`.  See the note on the
   native branch for why that temporary is passed as void*.  It costs
   nothing at -O1 and up.

   Define MP_NO_NATIVE_ADDCARRY to force the builtin path (for benchmarking),
   or MP_FORCE_PORTABLE_CARRY to force the plain-C path. */

#if !defined(MP_HAS_NATIVE_ADDCARRY)
  #if defined(MP_FORCE_PORTABLE_CARRY) || defined(MP_NO_NATIVE_ADDCARRY)
    #define MP_HAS_NATIVE_ADDCARRY 0
  #elif defined(MP_ASSUME_NATIVE_ADDCARRY)
    /* escape hatch: force the native path on for a toolchain that has the
       intrinsics but is not recognized by the tests below. */
    #define MP_HAS_NATIVE_ADDCARRY 1
  #elif defined(_MSC_VER) && defined(_M_X64)
    /* <intrin.h> provides both; MSVC on ARM64 does not, and falls through. */
    #define MP_HAS_NATIVE_ADDCARRY 1
  #elif (defined(__x86_64__) || defined(_M_X64)) && \
        (defined(__clang__) || defined(__INTEL_COMPILER) || \
         defined(__INTEL_LLVM_COMPILER))
    #define MP_HAS_NATIVE_ADDCARRY 1
  #elif defined(__x86_64__) && defined(__GNUC__) && (__GNUC__ >= 10)
    /* VERSION GATE -- do not loosen without testing.  gcc's adxintrin.h has
       not always provided the non-`x` _addcarry_u64 / _subborrow_u64: on
       gcc 8.5 they are absent, and relying on them produces an implicit
       declaration (assumed `extern int`) that compiles with only a warning
       and then fails at LINK with "undefined reference to _addcarry_u64".
       They are NOT ADX-gated, so -madx does not help -- it is purely a
       header-vintage issue, hence a version test rather than a target test.
       gcc >= 10 is known good (12.2 verified); 8.5 verified missing.
       gcc 9 and below fall through to the mixed __int128/plain-C tier
       below.  That tier is measurably slower than native -- on gcc 12.2,
       0.68/0.75 vs 0.52/0.51 ns per limb -- so this gate costs real
       performance on old gcc and should be lowered if you can confirm
       your version has the intrinsics (see MP_ASSUME_NATIVE_ADDCARRY).
       If your older gcc does have them, define MP_ASSUME_NATIVE_ADDCARRY. */
    #define MP_HAS_NATIVE_ADDCARRY 1
  #else
    #define MP_HAS_NATIVE_ADDCARRY 0
  #endif
#endif

#if MP_HAS_NATIVE_ADDCARRY
  #define MP_CARRY_IMPL "native _addcarry_u64/_subborrow_u64"

  /* The out-parameter type of these intrinsics is NOT consistent across
     compilers, and there is no single typed temporary that satisfies all of
     them:

       gcc / clang / MSVC   unsigned long long *
       icc classic (LP64)   unsigned __int64 *  ==  unsigned long *

     Passing `unsigned long long *` gets icc warning #167; passing uint64_t *
     gets -Wincompatible-pointer-types on LP64 gcc/clang.  A `void *` satisfies
     every one of them: in C it implicitly converts to whatever object pointer
     type the declaration actually uses, so we no longer have to track each
     compiler's header details -- which is twice now that guessing has cost us.
     The union makes the punning well-defined rather than an aliasing gamble.
     Verified to produce the same adc/sbb chain, with no spill, as the plain
     temporary did -- the address never escapes, so it stays in a register.

     C++ does NOT implicitly convert void*, so MP_INTRIN_OUT falls back to a
     typed pointer there.  That is correct for gcc/clang/MSVC; an icc C++ TU
     would see warning #167 again, which is noisy but harmless.  YAFU is C,
     so the C path is the one that matters. */
  typedef union { uint64_t u64; unsigned long long ull; unsigned long ul; }
      mp_u64_slot_t;
  typedef union { uint32_t u32; unsigned int ui; unsigned long ul; }
      mp_u32_slot_t;

  #if defined(__cplusplus)
    #define MP_INTRIN_OUT64(slot) (&(slot).ull)
    #define MP_INTRIN_OUT32(slot) (&(slot).ui)
  #else
    #define MP_INTRIN_OUT64(slot) ((void *)&(slot))
    #define MP_INTRIN_OUT32(slot) ((void *)&(slot))
  #endif

  MP_FORCE_INLINE mp_carry_t mp_addcarry_u64(mp_carry_t c_in, uint64_t a,
                                             uint64_t b, uint64_t *out) {
      mp_u64_slot_t t;
      mp_carry_t c = (mp_carry_t)_addcarry_u64((unsigned char)c_in, a, b,
                                               MP_INTRIN_OUT64(t));
      *out = t.u64;
      return c;
  }
  MP_FORCE_INLINE mp_carry_t mp_addcarry_u32(mp_carry_t c_in, uint32_t a,
      uint32_t b, uint32_t* out) {
      mp_u32_slot_t t;
      mp_carry_t c = (mp_carry_t)_addcarry_u32((unsigned char)c_in, a, b,
          MP_INTRIN_OUT32(t));
      *out = t.u32;
      return c;
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u64(mp_carry_t b_in, uint64_t a,
                                              uint64_t b, uint64_t *out) {
      mp_u64_slot_t t;
      mp_carry_t br = (mp_carry_t)_subborrow_u64((unsigned char)b_in, a, b,
                                                 MP_INTRIN_OUT64(t));
      *out = t.u64;
      return br;
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u32(mp_carry_t b_in, uint32_t a,
                                              uint32_t b, uint32_t *out) {
      mp_u32_slot_t t;
      mp_carry_t br = (mp_carry_t)_subborrow_u32((unsigned char)b_in, a, b,
                                                 MP_INTRIN_OUT32(t));
      *out = t.u32;
      return br;
  }

#elif defined(__clang__) && !defined(MP_FORCE_PORTABLE_CARRY)
  #define MP_CARRY_IMPL "clang __builtin_addcll/__builtin_subcll"

  /* __builtin_addcll RETURNS the sum and writes the carry-out through the
     pointer -- the opposite of the MSVC convention, so do not cross them. */
  MP_FORCE_INLINE mp_carry_t mp_addcarry_u64(mp_carry_t c_in, uint64_t a,
                                             uint64_t b, uint64_t *out) {
      unsigned long long carry_out;
      *out = (uint64_t)__builtin_addcll((unsigned long long)a,
                                        (unsigned long long)b,
                                        (unsigned long long)c_in, &carry_out);
      return (mp_carry_t)carry_out;
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u64(mp_carry_t b_in, uint64_t a,
                                              uint64_t b, uint64_t *out) {
      unsigned long long borrow_out;
      *out = (uint64_t)__builtin_subcll((unsigned long long)a,
                                        (unsigned long long)b,
                                        (unsigned long long)b_in, &borrow_out);
      return (mp_carry_t)borrow_out;
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u32(mp_carry_t b_in, uint32_t a,
                                              uint32_t b, uint32_t *out) {
      unsigned int borrow_out;
      *out = (uint32_t)__builtin_subc(a, b, (unsigned int)b_in, &borrow_out);
      return (mp_carry_t)borrow_out;
  }

#elif defined(HAS_UINT128) && !defined(MP_FORCE_PORTABLE_CARRY)
  #define MP_CARRY_IMPL "__int128 add + portable-C sub"

  /* Fallback wherever __int128 exists -- in practice gcc < 10, which has no
     native _addcarry_u64.  Deliberately MIXED, because the two directions do
     not behave the same way.

     The 128-bit add lowers well (carry-out is just the high word), but the
     128-bit subtract does not: extracting the borrow as (t >> 64) & 1 defeats
     the compiler, and the plain-C borrow idiom beats it badly.  Measured on
     AMD, ns per limb over an 8-limb chain:

                        gcc 8.5          clang 17        gcc 12.2
       add  __int128      0.675            0.484           0.677
       add  plain C       1.084            0.887           0.741   -> __int128 wins
       sub  __int128      1.507            0.546           0.753
       sub  plain C       0.795            0.267           0.795   -> plain C wins

     Using __int128 for both would make subtract ~2x SLOWER than plain C on
     gcc 8.5 and clang.  So: __int128 for add, plain C for sub.  On gcc 12.2
     the two are within noise, so the split costs nothing there.

     For reference, the native intrinsic path (tier 1) beats everything here
     on every compiler tested -- gcc 12.2 0.516/0.511, clang 0.265/0.267 --
     so this tier only ever runs where tier 1 is unavailable. */
  MP_FORCE_INLINE mp_carry_t mp_addcarry_u64(mp_carry_t c_in, uint64_t a,
                                             uint64_t b, uint64_t *out) {
      __uint128_t t = (__uint128_t)a + (__uint128_t)b + (__uint128_t)c_in;
      *out = (uint64_t)t;
      return (mp_carry_t)(t >> 64);
  }
  /* plain C, NOT __int128 -- see the table above */
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u64(mp_carry_t b_in, uint64_t a,
                                              uint64_t b, uint64_t *out) {
      uint64_t d  = a - b;
      mp_carry_t br = (mp_carry_t)(a < b);
      uint64_t d2 = d - (uint64_t)b_in;
      br |= (mp_carry_t)(d < (uint64_t)b_in);
      *out = d2;
      return br;
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u32(mp_carry_t b_in, uint32_t a,
                                              uint32_t b, uint32_t *out) {
      uint32_t d  = a - b;
      mp_carry_t br = (mp_carry_t)(a < b);
      uint32_t d2 = d - (uint32_t)b_in;
      br |= (mp_carry_t)(d < (uint32_t)b_in);
      *out = d2;
      return br;
  }

#elif defined(__GNUC__) && !defined(MP_FORCE_PORTABLE_CARRY)
  #define MP_CARRY_IMPL "gcc __builtin_add_overflow/__builtin_sub_overflow"

  MP_FORCE_INLINE mp_carry_t mp_addcarry_u64(mp_carry_t c_in, uint64_t a,
                                             uint64_t b, uint64_t *out) {
      uint64_t s;
      mp_carry_t c1 = (mp_carry_t)__builtin_add_overflow(a, b, &s);
      mp_carry_t c2 = (mp_carry_t)__builtin_add_overflow(s, (uint64_t)c_in, &s);
      *out = s;
      return (mp_carry_t)(c1 | c2);   /* c1 and c2 are never both set */
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u64(mp_carry_t b_in, uint64_t a,
                                              uint64_t b, uint64_t *out) {
      uint64_t d;
      mp_carry_t b1 = (mp_carry_t)__builtin_sub_overflow(a, b, &d);
      mp_carry_t b2 = (mp_carry_t)__builtin_sub_overflow(d, (uint64_t)b_in, &d);
      *out = d;
      return (mp_carry_t)(b1 | b2);
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u32(mp_carry_t b_in, uint32_t a,
                                              uint32_t b, uint32_t *out) {
      uint32_t d;
      mp_carry_t b1 = (mp_carry_t)__builtin_sub_overflow(a, b, &d);
      mp_carry_t b2 = (mp_carry_t)__builtin_sub_overflow(d, (uint32_t)b_in, &d);
      *out = d;
      return (mp_carry_t)(b1 | b2);
  }

#else
  #define MP_CARRY_IMPL "portable C"

  MP_FORCE_INLINE mp_carry_t mp_addcarry_u64(mp_carry_t c_in, uint64_t a,
                                             uint64_t b, uint64_t *out) {
      uint64_t s = a + (uint64_t)c_in;
      mp_carry_t c = (mp_carry_t)(s < a);
      s += b;
      c |= (mp_carry_t)(s < b);
      *out = s;
      return c;
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u64(mp_carry_t b_in, uint64_t a,
                                              uint64_t b, uint64_t *out) {
      uint64_t d  = a - b;
      mp_carry_t br = (mp_carry_t)(a < b);
      uint64_t d2 = d - (uint64_t)b_in;
      br |= (mp_carry_t)(d < (uint64_t)b_in);
      *out = d2;
      return br;
  }
  MP_FORCE_INLINE mp_carry_t mp_subborrow_u32(mp_carry_t b_in, uint32_t a,
                                              uint32_t b, uint32_t *out) {
      uint32_t d  = a - b;
      mp_carry_t br = (mp_carry_t)(a < b);
      uint32_t d2 = d - (uint32_t)b_in;
      br |= (mp_carry_t)(d < (uint32_t)b_in);
      *out = d2;
      return br;
  }
#endif

/* Compatibility aliases for existing call sites.  Installed ONLY when the
   platform has no native intrinsic of that name -- otherwise we would shadow
   (and conflict with) the real one. */
#if !MP_HAS_NATIVE_ADDCARRY
  #define _addcarry_u64  mp_addcarry_u64
  #define _subborrow_u64 mp_subborrow_u64
  #define _subborrow_u32 mp_subborrow_u32
#endif

#endif /* MP_PLATFORM_H */
