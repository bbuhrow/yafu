/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@gmail.com
 

 * Modified:	Ben Buhrow
 * Date:		11/24/09
 * Purpose:		Port into Yafu-1.14.


*/

#define TFM_DEFINES
#include "fp_sqr_comba.c"
#include "yafu.h"
#include "util.h"

/* generic comba squarer */

void fp_sqr_comba(z *A, z *B)
{
  int       pa, ix, iz, sA;
  fp_digit  c0, c1, c2;
  z    *dst;
  z loc;
#ifdef TFM_ISO
  uint64   tt;
#endif    

  /* get size of output and trim */
  sA = abs(A->size);
   pa = sA + sA;

  /* number of output digits to produce */
  COMBA_START;
  CLEAR_CARRY;

  if (A == B) {
     //zClear(&atmp1);
	  zInit(&loc);
     //dst = &atmp1;
	  dst = &loc;
  } else {
     zClear(B);
     dst = B;
  }

	if (dst->alloc < pa)
	{
		zGrow(dst,pa + LIMB_BLKSZ);
	}
	zClear(dst);
	
  for (ix = 0; ix < pa; ix++) { 
      int      tx, ty, iy;
      fp_digit *tmpy, *tmpx;

      /* get offsets into the two bignums */
      ty = MIN(sA-1, ix);
      tx = ix - ty;

      /* setup temp aliases */
      tmpx = A->val + tx;
      tmpy = A->val + ty;

      /* this is the number of times the loop will iterrate,
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MIN(sA-tx, ty+1);

      /* now for squaring tx can never equal ty 
       * we halve the distance since they approach 
       * at a rate of 2x and we have to round because 
       * odd cases need to be executed
       */
      iy = MIN(iy, (ty-tx+1)>>1);

      /* forward carries */
      CARRY_FORWARD;

      /* execute loop */
      for (iz = 0; iz < iy; iz++) {
		  SQRADD2(*tmpx++, *tmpy--);
      }

      /* even columns have the square term in them */
      if ((ix&1) == 0) {
		  SQRADD(A->val[ix>>1],A->val[ix>>1]);
      }

      /* store it */
      COMBA_STORE(dst->val[ix]);
  }

  COMBA_FINI;

  /* setup dest */
  dst->size = pa;
  fp_clamp (dst);
  if (dst != B) {
     zCopy(dst, B);
	 zFree(&loc);
  }
}



/* $Source: /cvs/libtom/tomsfastmath/src/sqr/Attic/fp_sqr_comba_generic.c,v $ */
/* $Revision: 1.3 $ */
/* $Date: 2007/02/15 00:31:32 $ */
