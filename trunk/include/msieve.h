/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/24/08

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.    Much of the functionality in here
			is ignored, as I'm only using the bits pertaining to
			block lanczos and certian other cpu utilities.
--------------------------------------------------------------------*/

#ifndef _MSIEVE_H_
#define _MSIEVE_H_

#ifdef __cplusplus
extern "C" {
#endif

	/* Lightweight factoring API */

#include "util.h"
#include "arith.h"
#include "yafu.h"

/* version info */




//void msieve_run(msieve_obj *obj);

#define MSIEVE_DEFAULT_LOGFILE "msieve.log"
#define MSIEVE_DEFAULT_SAVEFILE "msieve.dat"
#define MSIEVE_DEFAULT_NFS_FBFILE "msieve.fb"

#ifdef __cplusplus
}
#endif

#endif /* _MSIEVE_H_ */

