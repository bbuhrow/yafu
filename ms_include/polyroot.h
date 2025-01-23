/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: polyroot.h 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _POLYROOT_H_
#define _POLYROOT_H_

#include <ms_common.h>
#include <dd.h>
#include <ddcomplex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* extended-precision polynomial rootfinder */

#define MAX_ROOTFINDER_DEGREE 10

uint32 find_poly_roots(dd_t *poly, uint32 degree, dd_complex_t *roots);

#ifdef __cplusplus
}
#endif

#endif /* _POLYROOT_H_ */
