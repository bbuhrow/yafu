#ifndef _COFACTORIZE_H_
#define _COFACTORIZE_H_

#include <gmp.h>

typedef signed char s8;
typedef unsigned char u8;
typedef signed short s16;
typedef unsigned short u16;
typedef signed int s32;
typedef unsigned int u32;

#ifdef _MSC_VER
typedef signed __int64 s64;
typedef unsigned __int64 u64;
#define INLINE _inline
#else
typedef long long s64;
typedef unsigned long long u64;
#define INLINE inline
#endif

static INLINE u64 gmp2u64(mpz_t src) 
{
  /* mpz_export is terribly slow */
  u64 ans = mpz_getlimbn(src, 0);
#if GMP_LIMB_BITS == 32
  if (mpz_size(src) >= 2)
    ans |= (u64)mpz_getlimbn(src, 1) << 32;
#endif
  return ans;
}

static INLINE void u64_2gmp(u64 src, mpz_t dest) 
{
#if GMP_LIMB_BITS == 64
  dest->_mp_d[0] = src;
  dest->_mp_size = (src ? 1 : 0);
#else
  /* mpz_import is terribly slow */
  mpz_set_ui(dest, (u32)(src >> 32));
  mpz_mul_2exp(dest, dest, 32);
  mpz_add_ui(dest, dest, (u32)src);
#endif
}


u32 tinyqs(mpz_t n, mpz_t factor1, mpz_t factor2);

#endif /* !_COFACTORIZE_H_ */
