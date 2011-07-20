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

/* unsigned subtraction ||a|| >= ||b|| ALWAYS! */

void s_fp_sub(z *a, z *b, z *c)
{
//  int      x, oldbused, oldused;
//  fp_word  t;

  
	//if (BITS_PER_DIGIT == 64)
	//{
		if (b->size < 0)
			printf("error, negative number in s_fp_sub\n");
		zSub(a,b,c);
		return;
	//}
	/*

  oldused  = c->size;
  oldbused = b->size;
  c->size  = a->size;
  t       = 0;
  for (x = 0; x < oldbused; x++) {
     t         = ((fp_word)a->val[x]) - (((fp_word)b->val[x]) + t);
     c->val[x]  = (fp_digit)t;
     t         = (t >> DIGIT_BIT)&1;
  }
  for (; x < a->size; x++) {
     t         = ((fp_word)a->val[x]) - t;
     c->val[x]  = (fp_digit)t;
     t         = (t >> DIGIT_BIT);
   }
  for (; x < oldused; x++) {
     c->val[x] = 0;
  }
  fp_clamp(c);
  */
}


/* $Source: /cvs/libtom/tomsfastmath/src/addsub/s_fp_sub.c,v $ */
/* $Revision: 1.1 $ */
/* $Date: 2006/12/31 21:25:53 $ */
