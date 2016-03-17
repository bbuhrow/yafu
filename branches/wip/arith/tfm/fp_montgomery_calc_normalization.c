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

/* computes a = B**n mod b without division or multiplication useful for
 * normalizing numbers in a Montgomery system.
 */
void fp_montgomery_calc_normalization(z *a, z *b)
{
  int     x, bits;

  /* how many bits of last digit does b use */
  bits = zBits(b) % BITS_PER_DIGIT;
  if (!bits) bits = BITS_PER_DIGIT;

  /* compute A = B^(n-1) * 2^(bits-1) */
  if (b->size > 1) {
     fp_2expt (a, (b->size - 1) * BITS_PER_DIGIT + bits - 1);
  } else {
	  //printf("b.size == 1\n");
     sp2z((fp_digit)1,a);
     bits = 1;
  }

  /* now compute C = A * B mod b */
  for (x = bits - 1; x < (int)BITS_PER_DIGIT; x++) {
    fp_mul_2 (a, a);
	//  zShiftLeft(a,a,1);
    if (fp_cmp_mag (a, b) != FP_LT) {
      s_fp_sub (a, b, a);
    }
	  //if (zCompare(a,b) > 0)
		//  zSub(a,b,a);
  }
}

/* $Source: /cvs/libtom/tomsfastmath/src/mont/fp_montgomery_calc_normalization.c,v $ */
/* $Revision: 1.1 $ */
/* $Date: 2006/12/31 21:25:53 $ */
