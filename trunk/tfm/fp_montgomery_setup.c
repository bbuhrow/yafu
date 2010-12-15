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

/* setups the montgomery reduction */
int fp_montgomery_setup(z *a, fp_digit *rho)
{
  fp_digit x, b;

/* fast inversion mod 2**k
 *
 * Based on the fact that
 *
 * XA = 1 (mod 2**n)  =>  (X(2-XA)) A = 1 (mod 2**2n)
 *                    =>  2*X*A - X*X*A*A = 1
 *                    =>  2*(1) - (1)     = 1
 */
  b = a->val[0];

  if ((b & 1) == 0) {
    return FP_VAL;
  }

  x = (((b + 2) & 4) << 1) + b; /* here x*a==1 mod 2**4 */
  x *= 2 - b * x;               /* here x*a==1 mod 2**8 */
  x *= 2 - b * x;               /* here x*a==1 mod 2**16 */
  x *= 2 - b * x;               /* here x*a==1 mod 2**32 */
  //printf("rho32 = %llu\n",(((fp_word) 1 << ((fp_word) 32)) - ((fp_word)x)));
#ifdef FP_64BIT
  x *= 2 - b * x;               /* here x*a==1 mod 2**64 */
#endif

  
  /* rho = -1/m mod b */
	*rho = (((fp_word) 1 << ((fp_word) DIGIT_BIT)) - ((fp_word)x));
	//printf("rho = %llu\n",*rho);

  return FP_OKAY;
}


/* $Source: /cvs/libtom/tomsfastmath/src/mont/fp_montgomery_setup.c,v $ */
/* $Revision: 1.1 $ */
/* $Date: 2006/12/31 21:25:53 $ */
