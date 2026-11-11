/*----------------------------------------------------------------------
mp_platform.h

Shared platform / detection / intrinsic-polyfill layer for the YAFU limb
arithmetic. limb1.{h,c}, limb2.{h,c}, microecm.c and micropm1.c include
this and nothing else for compiler / OS / ISA concerns.

Self-contained: includes only <stdint.h> and the compiler's own intrinsic
headers.  It depends on nothing else in YAFU, so the arith/limb directory
can be copied into another project as-is.

Contents:
  - target / compiler / ISA-feature normalization (MP_ARCH_*, MP_USE_*,
    MP_HAS_GNU_ASM_X64)
  - inline / unused attribute macros
  - ONE copy of every intrinsic polyfill: _umul128, __umulh, _udiv128,
    _addcarry_u64/_u32, _subborrow_u64/_u32, the {trail,lead}_zcnt /
    reset_lsb family, popcount.

This replaces the per-file #if ladders that were re-derived (and had drifted)
in arith.h, arith.c, monty.h, monty.c, microecm.c and micropm1.c.
----------------------------------------------------------------------*/

#ifndef MP_PLATFORM_H
#define MP_PLATFORM_H

#include <stdint.h>

/* -------- target architecture -------- */
#if defined(__x86_64__) || defined(_M_X64)
  #define MP_ARCH_X64 1
#else
  #define MP_ARCH_X64 0
#endif
#if defined(__aarch64__) || defined(_M_ARM64)
  #define MP_ARCH_ARM64 1
#else
  #define MP_ARCH_ARM64 0
#endif

/* -------- vendor intrinsics: pulled in once, here -------- */
#if defined(_MSC_VER)
  #include <intrin.h>
  #if MP_ARCH_X64
    #include <immintrin.h>
    #if !defined(__clang__)
      #pragma intrinsic(_umul128)
      #pragma intrinsic(_udiv128)
    #endif
  #endif
#elif MP_ARCH_X64 || defined(__i386__)
  #include <immintrin.h>
#endif

/* -------- ISA feature selection --------
   One MP_USE_* flag per feature the limb code cares about.  Each is on when
   the host project asks for it explicitly (YAFU's -DUSE_BMI2 / -DUSE_AVX2 /
   -DUSE_AVX512F / -DIFMA build flags) or when the compiler has the feature
   enabled (-mbmi2, -march=..., /arch:AVX512, ...), so the library works both
   inside YAFU and dropped into another project with no -D flags at all.
   Define MP_NO_AUTO_ISA to disable the compiler-macro detection and rely on
   the explicit flags only.

   If the host passes ANY of the explicit flags, auto-detection turns itself
   off: a project that selects its ISA with -D flags has code of its own
   gated on those same flags, and the library must not quietly pick a
   different path than its caller (e.g. -march=native -DUSE_AVX512F without
   -DIFMA on an IFMA machine means "build the non-IFMA path"). */
#if !defined(MP_NO_AUTO_ISA) && (defined(USE_SSE41) || defined(USE_AVX2) \
    || defined(USE_BMI2) || defined(USE_AVX512F) || defined(IFMA) || defined(TARGET_KNL))
  #define MP_NO_AUTO_ISA
#endif
#if !MP_ARCH_X64
  #define MP_USE_BMI2    0
  #define MP_USE_AVX2    0
  #define MP_USE_AVX512F 0
  #define MP_USE_IFMA    0
#else
  #if defined(USE_BMI2) || defined(TARGET_KNL) || defined(USE_AVX512F) \
      || (!defined(MP_NO_AUTO_ISA) && defined(__BMI__) && defined(__BMI2__))
    #define MP_USE_BMI2 1
  #else
    #define MP_USE_BMI2 0
  #endif
  #if defined(USE_AVX2) || defined(USE_AVX512F) \
      || (!defined(MP_NO_AUTO_ISA) && defined(__AVX2__))
    #define MP_USE_AVX2 1
  #else
    #define MP_USE_AVX2 0
  #endif
  #if defined(USE_AVX512F) \
      || (!defined(MP_NO_AUTO_ISA) && defined(__AVX512F__) && defined(__AVX512DQ__))
    #define MP_USE_AVX512F 1
  #else
    #define MP_USE_AVX512F 0
  #endif
  #if MP_USE_AVX512F && (defined(IFMA) \
      || (!defined(MP_NO_AUTO_ISA) && defined(__AVX512IFMA__)))
    #define MP_USE_IFMA 1
  #else
    #define MP_USE_IFMA 0
  #endif
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
   True for gcc, mingw, clang (including clang-cl) and both Intel compilers on
   x86-64; false for native MSVC and for every non-x86 target.  This is the
   same set YAFU's common.h calls GCC_ASM64X (plus __MINGW64__, which the old
   monty.h gate added), derived here directly from compiler macros so this
   header needs nothing from the host project.  ASM_ARITH_DEBUG forces the
   portable path so the two code paths can be A/B compared. */
#if MP_ARCH_X64 && !defined(ASM_ARITH_DEBUG) \
    && (defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) \
        || defined(__INTEL_LLVM_COMPILER))
  #define MP_HAS_GNU_ASM_X64 1
#else
  #define MP_HAS_GNU_ASM_X64 0
#endif

/* Note: this header deliberately does NOT define a generic INLINE_ASM macro.
   Other YAFU files test that name with their own meaning, and a library
   header must not change what they compile. */


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
#if defined(_MSC_VER) && !defined(__clang__)   /* clang-cl takes the attribute form */
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
   3. Bit scan / reset : Make every path a value-returning inline
   ====================================================================== */

#if MP_USE_BMI2
    /* BMI branch (where _blsr_u64/_blsr_u32 are available) */
    MP_FORCE_INLINE uint64_t _reset_lsb64(uint64_t x) { return _blsr_u64(x); }
    MP_FORCE_INLINE uint32_t _reset_lsb(uint32_t x) { return _blsr_u32(x); }
#else
    /* portable fallback branch */
    MP_FORCE_INLINE uint64_t _reset_lsb64(uint64_t x) { return x & (x - 1); }
    MP_FORCE_INLINE uint32_t _reset_lsb(uint32_t x) { return x & (x - 1); }
#endif

#if defined(__INTEL_COMPILER)
  #if MP_USE_BMI2
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

/* Software 128/64 -> 64 divide: Knuth algorithm D specialised to a two-digit
   quotient in base 2^32 (Hacker's Delight "divlu").  Contract as for _udiv128:
   d != 0 and hi < d.  This is the _udiv128 of targets with neither a hardware
   128/64 divide nor __int128 (MSVC on ARM64); it is defined everywhere so it
   can be tested everywhere.  Not force-inlined -- it is big. */
static UNUSED_FUNC uint64_t mp_udiv128_soft(uint64_t hi, uint64_t lo, uint64_t d,
                                         uint64_t *rem) {
    const uint64_t b = (uint64_t)1 << 32;
    uint64_t un32, un10, un1, un0, vn1, vn0, q1, q0, rhat;
    int s = 0;
    while (!(d >> 63)) { d <<= 1; s++; }        /* normalize divisor */
    vn1 = d >> 32;
    vn0 = d & 0xffffffffULL;
    un32 = s ? ((hi << s) | (lo >> (64 - s))) : hi;
    un10 = lo << s;
    un1 = un10 >> 32;
    un0 = un10 & 0xffffffffULL;
    q1 = un32 / vn1;
    rhat = un32 - q1 * vn1;
    while (q1 >= b || q1 * vn0 > b * rhat + un1) {
        q1--; rhat += vn1;
        if (rhat >= b) break;
    }
    un32 = un32 * b + un1 - q1 * d;              /* mod 2^64 on purpose */
    q0 = un32 / vn1;
    rhat = un32 - q0 * vn1;
    while (q0 >= b || q0 * vn0 > b * rhat + un0) {
        q0--; rhat += vn1;
        if (rhat >= b) break;
    }
    *rem = (un32 * b + un0 - q0 * d) >> s;
    return q1 * b + q0;
}

#if !(defined(_MSC_VER) && MP_ARCH_X64)

  #if defined(__MINGW32__)
    /* mingw-w64's <windows.h> (winnt.h -> psdk_inc/intrin-impl.h) declares
       __umulh and declares+defines _umul128 with external linkage, which
       collides with our static versions.  Route ours through private names so
       either include order works, and tell intrin-impl.h its _umul128 already
       exists so a later <windows.h> doesn't define it again.  (mingw has no
       __umulh body at all, so we can't just defer to theirs.) */
    #define _umul128 mp_umul128
    #define __umulh  mp_umulh
    #define __INTRINSIC_DEFINED__umul128
  #endif

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
  #elif defined(_MSC_VER) && MP_ARCH_ARM64
    /* MSVC on ARM64: <intrin.h> has __umulh (one UMULH instruction) but no
       _umul128, so build the latter from it. */
    MP_FORCE_INLINE uint64_t _umul128(uint64_t x, uint64_t y, uint64_t *hi) {
        *hi = __umulh(x, y);
        return x * y;
    }
  #else
    /* portable 32x32 schoolbook -- only reached on a target with neither
       __int128 nor a 64x64->128 intrinsic. */
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
    /* no hardware or __int128 path (MSVC on ARM64, or an unknown compiler) */
    #define _udiv128 mp_udiv128_soft
  #endif

#endif /* !(_MSC_VER && x64) */

/* clang-cl on x64: its <intrin.h> supplies _umul128 and __umulh but not
   _udiv128 (verified missing in the LLVM shipped with VS 2022), so calling it
   is an implicit-declaration error.  Supply the divq form under our own name
   and route _udiv128 to it; if a later clang-cl grows a real _udiv128 the
   macro simply shadows it, which is harmless. */
#if defined(_MSC_VER) && (defined(__clang__) || defined(__INTEL_COMPILER)) \
    && MP_ARCH_X64      /* clang-cl, icx-cl, and classic Intel icl */
  MP_FORCE_INLINE uint64_t mp_udiv128_clangcl(uint64_t hi, uint64_t lo, uint64_t d,
                                              uint64_t *rem) {
  #if MP_HAS_GNU_ASM_X64
      uint64_t q, r;
      __asm__("divq %[d]"
              : "=a"(q), "=d"(r)
              : "a"(lo), "d"(hi), [d]"rm"(d)
              : "cc");
      *rem = r; return q;
  #else     /* ASM_ARITH_DEBUG */
      return mp_udiv128_soft(hi, lo, d, rem);
  #endif
  }
  #undef _udiv128
  #define _udiv128 mp_udiv128_clangcl
#endif

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
  /* 32-bit add-with-carry: one definition for every non-native tier (a
     64-bit add is the best any of them can do). */
  MP_FORCE_INLINE mp_carry_t mp_addcarry_u32(mp_carry_t c_in, uint32_t a,
                                             uint32_t b, uint32_t *out) {
      uint64_t s = (uint64_t)a + (uint64_t)b + (uint64_t)c_in;
      *out = (uint32_t)s;
      return (mp_carry_t)(s >> 32);
  }

  #undef _addcarry_u64
  #undef _addcarry_u32
  #undef _subborrow_u64
  #undef _subborrow_u32
  #define _addcarry_u64  mp_addcarry_u64
  #define _addcarry_u32  mp_addcarry_u32
  #define _subborrow_u64 mp_subborrow_u64
  #define _subborrow_u32 mp_subborrow_u32
#endif


/* ======================================================================
   6. Population count
   ====================================================================== */
#if (defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) \
     || defined(__INTEL_LLVM_COMPILER)) \
    && !(defined(_MSC_VER) && defined(__INTEL_COMPILER))  /* icl on Windows has no __builtin_popcount */
  MP_FORCE_INLINE uint32_t mp_popcount64(uint64_t x) {
      return (uint32_t)__builtin_popcountll(x);
  }
  MP_FORCE_INLINE uint32_t mp_popcount32(uint32_t x) {
      return (uint32_t)__builtin_popcount(x);
  }
#elif defined(_MSC_VER) && MP_ARCH_X64
  /* popcnt is not baseline x86-64, but every CPU this code targets has it */
  MP_FORCE_INLINE uint32_t mp_popcount64(uint64_t x) {
      return (uint32_t)__popcnt64(x);
  }
  MP_FORCE_INLINE uint32_t mp_popcount32(uint32_t x) {
      return (uint32_t)__popcnt(x);
  }
#else
  MP_FORCE_INLINE uint32_t mp_popcount64(uint64_t x) {
      x = x - ((x >> 1) & 0x5555555555555555ULL);
      x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
      x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
      return (uint32_t)((x * 0x0101010101010101ULL) >> 56);
  }
  MP_FORCE_INLINE uint32_t mp_popcount32(uint32_t x) {
      return mp_popcount64((uint64_t)x);
  }
#endif

#endif /* MP_PLATFORM_H */
