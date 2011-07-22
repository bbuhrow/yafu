/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: gmp_xface.h 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#ifdef HAVE_GMP

#ifndef _GMP_XFACE_H_
#define _GMP_XFACE_H_

#include "util.h"
#include <gmp.h>
#include "yafu.h"

#ifdef __cplusplus
extern "C" {
#endif

	/* Note that when GMP_LIMB_BITS == 64 it is possible
	   to use mpz_set_{ui|si}, except that 64-bit
	   MSVC forces the input argument in these calls to
	   be 32 bits in size and not 64 */

/*--------------------------------------------------------------------*/
static INLINE void mp2gmp(z *src, mpz_t dest) {

	mpz_import(dest, (size_t)(abs(src->size)), -1, sizeof(fp_digit), 
			0, (size_t)0, src->val);

	if (src->size < 0)
		mpz_neg(dest, dest);
}

/*--------------------------------------------------------------------*/
static INLINE void gmp2mp(mpz_t src, z *dest) {

	size_t count;

	zClear(dest);
	mpz_export(dest->val, &count, -1, sizeof(fp_digit),
			0, (size_t)0, src);
	dest->size = count;
}

/*--------------------------------------------------------------------*/
static INLINE void uint64_2gmp(uint64 src, mpz_t dest) {

#if GMP_LIMB_BITS == 64
	dest->_mp_d[0] = src;
	dest->_mp_size = (src ? 1 : 0);
#else
	/* mpz_import is terribly slow */
	mpz_set_ui(dest, (uint32)(src >> 32));
	mpz_mul_2exp(dest, dest, 32);
	mpz_add_ui(dest, dest, (uint32)src);
#endif
}

/*--------------------------------------------------------------------*/
static INLINE void int64_2gmp(int64 src, mpz_t dest) {

	if (src < 0) {
		uint64_2gmp((uint64)(-src), dest);
		mpz_neg(dest, dest);
	}
	else {
		uint64_2gmp((uint64)src, dest);
	}
}

/*--------------------------------------------------------------------*/
static INLINE uint64 gmp2uint64(mpz_t src) {

	/* mpz_export is terribly slow */
	uint64 ans = mpz_getlimbn(src, 0);
#if GMP_LIMB_BITS == 32
	if (mpz_size(src) >= 2)
		ans |= (uint64)mpz_getlimbn(src, 1) << 32;
#endif
        return ans;
}

/*--------------------------------------------------------------------*/
static INLINE int64 gmp2int64(mpz_t src) {

	if (mpz_cmp_ui(src, 0) < 0) {
		return -gmp2uint64(src);
	}
	else {
       		return gmp2uint64(src);
	}
}


#ifdef __cplusplus
}
#endif

#endif /* _GMP_XFACE_H_ */

#endif /* HAVE_GMP */
