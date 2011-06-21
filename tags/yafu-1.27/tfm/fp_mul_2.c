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

void fp_mul_2(z * a, z * b)
{
  int     x, oldused;

  int sb = abs(b->size);
  int sa = abs(a->size);
  int signa = a->size < 0;
//  int signb = b->size < 0;
   
  oldused = sb;
  b->size = sa;

  {
    fp_digit r, rr, *tmpa, *tmpb;

    /* alias for source */
    tmpa = a->val;
    
    /* alias for dest */
    tmpb = b->val;

    /* carry */
    r = 0;
    for (x = 0; x < sa; x++) {
    
      /* get what will be the *next* carry bit from the 
       * MSB of the current digit 
       */
      rr = *tmpa >> ((fp_digit)(DIGIT_BIT - 1));
      
      /* now shift up this digit, add in the carry [from the previous] */
      *tmpb++ = ((*tmpa++ << ((fp_digit)1)) | r);
      
      /* copy the carry that would be from the source 
       * digit into the next iteration 
       */
      r = rr;
    }

    /* new leading digit? */
    if (r != 0 && b->size != (FP_SIZE-1)) {
      /* add a MSB which is always 1 at this point */
      *tmpb = 1;
      ++(b->size);
    }

    /* now zero any excess digits on the destination 
     * that we didn't write to 
     */
    tmpb = b->val + b->size;
    for (x = b->size; x < oldused; x++) {
      *tmpb++ = 0;
    }
  }
  if (signa)
	b->size *= -1;
}


/* $Source: /cvs/libtom/tomsfastmath/src/mul/fp_mul_2.c,v $ */
/* $Revision: 1.1 $ */
/* $Date: 2006/12/31 21:25:53 $ */
