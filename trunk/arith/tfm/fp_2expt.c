/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@gmail.com
 *
 * Modified:	Ben Buhrow
 * Date:		11/24/09
 * Purpose:		Port into Yafu-1.14.  
 */
#include "tfm.h"
#include "yafu.h"

/* computes a = 2**b */
void fp_2expt(z *a, int b)
{
   int     k;

   /* zero a as per default */
   zClear(a);

   if (b < 0) { 
      return;
   }

   k = b / DIGIT_BIT;

  /* set the used count of where the bit will go */
  a->size = k + 1;
  if (a->size >= a->alloc)
	  zGrow(a,a->size + 2);

  /* put the single bit in its place */
  a->val[k] = ((fp_digit)1) << (b % DIGIT_BIT);
}


/* $Source: /cvs/libtom/tomsfastmath/src/exptmod/fp_2expt.c,v $ */
/* $Revision: 1.1 $ */
/* $Date: 2006/12/31 21:25:53 $ */
