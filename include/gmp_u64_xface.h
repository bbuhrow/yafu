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

#ifndef _GMP_U64_XFACE_H_
#define _GMP_U64_XFACE_H_

#include "ytools.h"
#include <gmp.h>
#include "gmp_u64.h"

#ifdef __cplusplus
extern "C" {
#endif

	/* Note that when GMP_LIMB_BITS == 64 it is possible
	   to use mpz_set_{ui|si}, except that 64-bit
	   MSVC forces the input argument in these calls to
	   be 32 bits in size and not 64 */

/*--------------------------------------------------------------------*/
static INLINE void uint64_2gmp(uint64_t src, mpz_t dest) {

	// 64-bit-safe on LLP64 (MinGW + GMP), see gmp_u64.h
	mpz_set_u64(dest, src);
}

/*--------------------------------------------------------------------*/
static INLINE void int64_2gmp(int64_t src, mpz_t dest) {

	if (src < 0) {
		uint64_2gmp((uint64_t)(-src), dest);
		mpz_neg(dest, dest);
	}
	else {
		uint64_2gmp((uint64_t)src, dest);
	}
}

/*--------------------------------------------------------------------*/
static INLINE uint64_t gmp2uint64(mpz_t src) {

	// 64-bit-safe on LLP64 (MinGW + GMP), see gmp_u64.h
	return mpz_get_u64(src);
}

/*--------------------------------------------------------------------*/
static INLINE int64_t gmp2int64(mpz_t src) {

	if (mpz_cmp_ui(src, 0) < 0) {
		return -(int64_t)gmp2uint64(src);
	}
	else {
       	return (int64_t)gmp2uint64(src);
	}
}


#ifdef __cplusplus
}
#endif

#endif //  _GMP_U64_XFACE_H_ 

