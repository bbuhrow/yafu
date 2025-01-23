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

#ifndef _GMP_XFACE_H_
#define _GMP_XFACE_H_

#include "msieve_common.h"
#include "gmp_u64_xface.h"
#include "ytools.h"
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

	/* Note that when GMP_LIMB_BITS == 64 it is possible
	   to use mpz_set_{ui|si}, except that 64-bit
	   MSVC forces the input argument in these calls to
	   be 32 bits in size and not 64 */

static INLINE void mp_t2gmp(mp_t *src, mpz_t dest) {

	mpz_import(dest, (size_t)(src->nwords), -1, sizeof(uint32_t),
			0, (size_t)0, src->val);
}

/*--------------------------------------------------------------------*/
static INLINE void gmp2mp_t(mpz_t src, mp_t *dest) {

	size_t count;

	memset(dest->val, 0, MAX_MP_WORDS * sizeof(uint32_t)); //mp_clear(dest);
	mpz_export(dest->val, &count, -1, sizeof(uint32_t),
			0, (size_t)0, src);
	dest->nwords = count;
}


#ifdef __cplusplus
}
#endif

#endif //  _GMP_XFACE_H_ 

