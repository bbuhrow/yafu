/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "arith.h"

/********************* montgomery arith **********************/
void monty_init(z *n);
void to_monty(z *x, z *n);
void monty_add(z *u, z *v, z *w, z *n);
void monty_mul(z *u, z *v, z *w, z *n);
void monty_sqr(z *x, z *w, z *n);
void monty_sub(z *u, z *v, z *w, z *n);
void monty_free(void);
void zREDC(z *T, z *n);
//monty exponentiation
void zmModExp(z *a, z *b, z *u, z *n);
void zmModExpw(z *a, z *e, z *u, z *n, int k);

monty montyconst;
