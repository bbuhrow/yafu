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
#include "tfm.h"
#include "yafu.h"

int fp_cmp_mag(z *a, z *b)
{
   int x;
	int sa = abs(a->size);
	int sb = abs(b->size);

   if (sa > sb) {
      return FP_GT;
   } else if (sa < sb) {
      return FP_LT;
   } else {
      for (x = sa - 1; x >= 0; x--) {
          if (a->val[x] > b->val[x]) {
             return FP_GT;
          } else if (a->val[x] < b->val[x]) {
             return FP_LT;
          }
      }
   }
   return FP_EQ;
}

/* $Source: /cvs/libtom/tomsfastmath/src/addsub/fp_cmp_mag.c,v $ */
/* $Revision: 1.1 $ */
/* $Date: 2006/12/31 21:25:53 $ */
