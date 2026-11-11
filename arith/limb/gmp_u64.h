/*
 * gmp_u64.h -- 64-bit-safe wrappers for GMP's mpz_*_ui functions.
 *
 * GMP's "_ui" functions take/return a C `unsigned long`.  That is 64 bits on
 * LP64 platforms (Linux, macOS) but only 32 bits on LLP64 Windows, so with
 * real GMP under MinGW/MSYS2 any uint64_t passed to mpz_set_ui, mpz_add_ui,
 * mpz_mul_ui, mpz_cmp_ui, mpz_tdiv_q_ui, ... is silently truncated to its low
 * 32 bits, and mpz_get_ui only returns the low 32 bits.  MPIR (used by the
 * MSVC builds) declares these with its 64-bit `mpir_ui` type on Win64, which
 * is why the problem only shows up in MinGW + GMP builds.
 *
 * Use the mpz_*_u64 forms below whenever the value may exceed 32 bits.  When
 * the GMP "ui" type is already 64 bits they compile to the plain _ui call.
 * Otherwise values that fit in 32 bits still take the _ui fast path and
 * larger ones go through a read-only one-limb mpz (no allocation).
 */

#ifndef YAFU_GMP_U64_H
#define YAFU_GMP_U64_H

#include <stdint.h>
#include <limits.h>
#include <gmp.h>

#if defined(__MPIR_VERSION) || (ULONG_MAX > 0xffffffffUL)
#define GMP_UI_IS_64 1
#else
#define GMP_UI_IS_64 0
#if GMP_LIMB_BITS != 64
#error "gmp_u64.h: 32-bit 'ui' with non-64-bit limbs is not supported"
#endif
#endif

#ifdef _MSC_VER
#define GMP_U64_INLINE static __inline
#else
#define GMP_U64_INLINE static inline
#endif

#if GMP_UI_IS_64

GMP_U64_INLINE void mpz_set_u64(mpz_ptr z, uint64_t v) { mpz_set_ui(z, v); }
GMP_U64_INLINE uint64_t mpz_get_u64(mpz_srcptr z) { return (uint64_t)mpz_getlimbn(z, 0); }
GMP_U64_INLINE void mpz_add_u64(mpz_ptr r, mpz_srcptr a, uint64_t v) { mpz_add_ui(r, a, v); }
GMP_U64_INLINE void mpz_sub_u64(mpz_ptr r, mpz_srcptr a, uint64_t v) { mpz_sub_ui(r, a, v); }
GMP_U64_INLINE void mpz_mul_u64(mpz_ptr r, mpz_srcptr a, uint64_t v) { mpz_mul_ui(r, a, v); }
GMP_U64_INLINE void mpz_addmul_u64(mpz_ptr r, mpz_srcptr a, uint64_t v) { mpz_addmul_ui(r, a, v); }
GMP_U64_INLINE void mpz_submul_u64(mpz_ptr r, mpz_srcptr a, uint64_t v) { mpz_submul_ui(r, a, v); }
GMP_U64_INLINE int mpz_cmp_u64(mpz_srcptr a, uint64_t v) { return mpz_cmp_ui(a, v); }
GMP_U64_INLINE uint64_t mpz_tdiv_q_u64(mpz_ptr q, mpz_srcptr a, uint64_t d) { return mpz_tdiv_q_ui(q, a, d); }
GMP_U64_INLINE uint64_t mpz_tdiv_r_u64(mpz_ptr r, mpz_srcptr a, uint64_t d) { return mpz_tdiv_r_ui(r, a, d); }
GMP_U64_INLINE uint64_t mpz_tdiv_u64(mpz_srcptr a, uint64_t d) { return mpz_tdiv_ui(a, d); }
GMP_U64_INLINE uint64_t mpz_fdiv_u64(mpz_srcptr a, uint64_t d) { return mpz_fdiv_ui(a, d); }
GMP_U64_INLINE void mpz_divexact_u64(mpz_ptr q, mpz_srcptr a, uint64_t d) { mpz_divexact_ui(q, a, d); }
GMP_U64_INLINE int mpz_divisible_u64_p(mpz_srcptr a, uint64_t d) { return mpz_divisible_ui_p(a, d); }

#else /* 32-bit 'unsigned long' (LLP64 Windows with GMP) */

/* read-only mpz view of a 64-bit value; v must outlive the view */
#define GMP_U64_VIEW(name, v) \
    mp_limb_t name##_limb = (mp_limb_t)(v); \
    mpz_t name; \
    mpz_roinit_n(name, &name##_limb, (name##_limb != 0) ? 1 : 0)

#define GMP_U64_FITS32(v) ((v) <= 0xffffffffULL)

GMP_U64_INLINE void mpz_set_u64(mpz_ptr z, uint64_t v)
{
    if (GMP_U64_FITS32(v)) { mpz_set_ui(z, (unsigned long)v); return; }
    { GMP_U64_VIEW(t, v); mpz_set(z, t); }
}

GMP_U64_INLINE uint64_t mpz_get_u64(mpz_srcptr z)
{
    return (uint64_t)mpz_getlimbn(z, 0);
}

GMP_U64_INLINE void mpz_add_u64(mpz_ptr r, mpz_srcptr a, uint64_t v)
{
    if (GMP_U64_FITS32(v)) { mpz_add_ui(r, a, (unsigned long)v); return; }
    { GMP_U64_VIEW(t, v); mpz_add(r, a, t); }
}

GMP_U64_INLINE void mpz_sub_u64(mpz_ptr r, mpz_srcptr a, uint64_t v)
{
    if (GMP_U64_FITS32(v)) { mpz_sub_ui(r, a, (unsigned long)v); return; }
    { GMP_U64_VIEW(t, v); mpz_sub(r, a, t); }
}

GMP_U64_INLINE void mpz_mul_u64(mpz_ptr r, mpz_srcptr a, uint64_t v)
{
    if (GMP_U64_FITS32(v)) { mpz_mul_ui(r, a, (unsigned long)v); return; }
    { GMP_U64_VIEW(t, v); mpz_mul(r, a, t); }
}

GMP_U64_INLINE void mpz_addmul_u64(mpz_ptr r, mpz_srcptr a, uint64_t v)
{
    if (GMP_U64_FITS32(v)) { mpz_addmul_ui(r, a, (unsigned long)v); return; }
    { GMP_U64_VIEW(t, v); mpz_addmul(r, a, t); }
}

GMP_U64_INLINE void mpz_submul_u64(mpz_ptr r, mpz_srcptr a, uint64_t v)
{
    if (GMP_U64_FITS32(v)) { mpz_submul_ui(r, a, (unsigned long)v); return; }
    { GMP_U64_VIEW(t, v); mpz_submul(r, a, t); }
}

GMP_U64_INLINE int mpz_cmp_u64(mpz_srcptr a, uint64_t v)
{
    if (GMP_U64_FITS32(v)) return mpz_cmp_ui(a, (unsigned long)v);
    { GMP_U64_VIEW(t, v); return mpz_cmp(a, t); }
}

/* the remainder-returning forms mirror the _ui semantics: |remainder| */
GMP_U64_INLINE uint64_t mpz_tdiv_q_u64(mpz_ptr q, mpz_srcptr a, uint64_t d)
{
    if (GMP_U64_FITS32(d)) return mpz_tdiv_q_ui(q, a, (unsigned long)d);
    {
        GMP_U64_VIEW(t, d);
        mpz_t r;
        uint64_t rem;
        mpz_init(r);
        mpz_tdiv_qr(q, r, a, t);
        rem = (uint64_t)mpz_getlimbn(r, 0);
        mpz_clear(r);
        return rem;
    }
}

GMP_U64_INLINE uint64_t mpz_tdiv_r_u64(mpz_ptr r, mpz_srcptr a, uint64_t d)
{
    if (GMP_U64_FITS32(d)) return mpz_tdiv_r_ui(r, a, (unsigned long)d);
    { GMP_U64_VIEW(t, d); mpz_tdiv_r(r, a, t); return (uint64_t)mpz_getlimbn(r, 0); }
}

GMP_U64_INLINE uint64_t mpz_tdiv_u64(mpz_srcptr a, uint64_t d)
{
    if (GMP_U64_FITS32(d)) return mpz_tdiv_ui(a, (unsigned long)d);
    {
        GMP_U64_VIEW(t, d);
        mpz_t r;
        uint64_t rem;
        mpz_init(r);
        mpz_tdiv_r(r, a, t);
        rem = (uint64_t)mpz_getlimbn(r, 0);
        mpz_clear(r);
        return rem;
    }
}

GMP_U64_INLINE uint64_t mpz_fdiv_u64(mpz_srcptr a, uint64_t d)
{
    if (GMP_U64_FITS32(d)) return mpz_fdiv_ui(a, (unsigned long)d);
    {
        GMP_U64_VIEW(t, d);
        mpz_t r;
        uint64_t rem;
        mpz_init(r);
        mpz_fdiv_r(r, a, t);
        rem = (uint64_t)mpz_getlimbn(r, 0);
        mpz_clear(r);
        return rem;
    }
}

GMP_U64_INLINE void mpz_divexact_u64(mpz_ptr q, mpz_srcptr a, uint64_t d)
{
    if (GMP_U64_FITS32(d)) { mpz_divexact_ui(q, a, (unsigned long)d); return; }
    { GMP_U64_VIEW(t, d); mpz_divexact(q, a, t); }
}

GMP_U64_INLINE int mpz_divisible_u64_p(mpz_srcptr a, uint64_t d)
{
    if (GMP_U64_FITS32(d)) return mpz_divisible_ui_p(a, (unsigned long)d);
    { GMP_U64_VIEW(t, d); return mpz_divisible_p(a, t); }
}

#endif /* GMP_UI_IS_64 */

/* signed 64-bit set, built on mpz_set_u64 */
GMP_U64_INLINE void mpz_set_s64(mpz_ptr z, int64_t v)
{
    if (v < 0) { mpz_set_u64(z, (uint64_t)0 - (uint64_t)v); mpz_neg(z, z); }
    else mpz_set_u64(z, (uint64_t)v);
}

#endif /* YAFU_GMP_U64_H */
