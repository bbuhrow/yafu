/*
Portable replacements for the x86 _addcarry_u64 / _subborrow_u64 /
_subborrow_u32 intrinsics, for non-x86 gcc/clang builds (e.g. aarch64).

These keep the x86 convention: the carry/borrow is RETURNED and the
result is written through the pointer.  Note that the gcc/clang builtins
__builtin_addcll / __builtin_subcll use the opposite convention (they
return the sum and write the carry through the pointer), so they cannot
be substituted with a simple argument-reordering macro.

Stopgap until mp_platform.h is adopted everywhere: it provides the same
aliases with the same convention.  If mp_platform.h is already included
this header does nothing; if it is included afterward, it #undefs and
replaces these aliases.  Once mp_platform.h is used throughout, replace
#include "portable_carry.h" with #include "mp_platform.h" and delete this
file.
*/

#ifndef YAFU_PORTABLE_CARRY_H
#define YAFU_PORTABLE_CARRY_H

#ifndef MP_PLATFORM_H

#include <stdint.h>

static inline uint8_t yafu_addcarry_u64(uint8_t c_in, uint64_t a, uint64_t b, uint64_t* out)
{
    __uint128_t s = (__uint128_t)a + b + c_in;
    *out = (uint64_t)s;
    return (uint8_t)(s >> 64);
}

static inline uint8_t yafu_subborrow_u64(uint8_t b_in, uint64_t a, uint64_t b, uint64_t* out)
{
    uint64_t d = a - b;
    uint8_t br = (a < b);
    *out = d - b_in;
    br |= (d < b_in);
    return br;
}

static inline uint8_t yafu_subborrow_u32(uint8_t b_in, uint32_t a, uint32_t b, uint32_t* out)
{
    uint32_t d = a - b;
    uint8_t br = (a < b);
    *out = d - b_in;
    br |= (d < b_in);
    return br;
}

#undef _addcarry_u64
#undef _subborrow_u64
#undef _subborrow_u32
#define _addcarry_u64  yafu_addcarry_u64
#define _subborrow_u64 yafu_subborrow_u64
#define _subborrow_u32 yafu_subborrow_u32

#endif // !MP_PLATFORM_H

#endif
