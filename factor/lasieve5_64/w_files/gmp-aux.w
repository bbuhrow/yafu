@* Auxilliary stuff for gmp.

Copyright (C) 2002 Jens Franke, T. Kleinjung
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@(gmp-aux.h@>=
void adjust_mpz_bufsize(mpz_t **x,size_t *alloc_ptr,size_t size,size_t increment);
int string2mpz(mpz_t rop,char *x,int base);
#if __GNU_MP_VERSION < 3 || (__GNU_MP_VERSION == 3 && __GNU_MP_VERSION_MINOR == 0)
// SMJS #error NEED_MPZ_MUL_SI
// SMJS
//#define NEED_MPZ_MUL_SI
//void mpz_mul_si(mpz_t x,mpz_t y,long int z);
#endif

#define mpz_add_si(r,o,s) \
  ({ int _s; _s=(s); \
     _s>=0 ? mpz_add_ui(r,o,(mp_limb_t)(_s)) : mpz_sub_ui(r,o,(mp_limb_t)(-_s)); })

void mpz_set_ull(mpz_t targ,ullong src);
ullong mpz_get_ull(mpz_t src);
int mpz_cmp_ull(mpz_t op1,ullong op2);
#ifdef ULL_NO_UL
void mpz_ull_init();
void mpz_mul_ull(mpz_t rop,mpz_t op1,ullong op2);
int mpz_fits_sllong_p(mpz_t);
int mpz_fits_ullong_p(mpz_t);
long long int mpz_get_sll(mpz_t);
void mpz_set_sll(mpz_t,long long int);
unsigned long long mpz_get_ull(mpz_t);
void mpz_set_ull(mpz_t,unsigned long long);
#else
#define mpz_ull_init()
#define mpz_mul_ull mpz_mul_ui
#define mpz_mul_sll mpz_mul_si
#define mpz_fits_sllong_p mpz_fits_slong_p
#define mpz_fits_ullong_p mpz_fits_ulong_p
#define mpz_get_sll mpz_get_si
#define mpz_set_sll mpz_set_si
#define mpz_get_ull mpz_get_ui
#define mpz_set_ull mpz_set_ui
#endif
@
@c
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <gmp.h>

#include "asm/siever-config.h"
#include "if.h"
#include "gmp-aux.h"

/* SMJS Commented out because doesn't seem to be used and complains about 
   pointer types - I don't want to think about how to fix it if its not used
void
adjust_mpz_bufsize(mpz_t **x,size_t *alloc_ptr,size_t size,size_t increment)
{
  size_t old_alloc;

  old_alloc=*alloc_ptr;
  adjust_bufsize(x,alloc_ptr,size,increment,sizeof(**x));
  while(old_alloc<*alloc_ptr) mpz_init((*x)[old_alloc++]);
}
*/

@ In some versions of gmp, |mpz_set_str(x,"0\n",10)| sets |x| to
a value different from zero. Also, a leading '+' may not be recognised.
The following function is inteded as an interface to |mpz_set_str|
which produces the expected behaviour on all valid input strings,
not as a general purpose function. For this
reason, we do not bother what |string2mpz(x," +-1 23")| does.

I do not know what happens if a |'-'| is separated from the first digit
by a whitespace.
@c
int
string2mpz(mpz_t rop,char *x,int base)
{
  size_t l;
  char *y;
  int rv;

  x+=strspn(x," \t+");
  if(strlen(x)==0) mpz_set_ui(rop,0);
  y=strdup(x);
  for(l=strlen(y)-1;y[l]=='\n';l--) {
    y[l]='\0';
    if(l==0) break;
  }
  rv=mpz_set_str(rop,y,base);
  free(y);
  return rv;
}

@
@c
#ifdef NEED_MPZ_MUL_SI
void mpz_mul_si(mpz_t x,mpz_t y,long int z)
{
  if(z<0) {
    mpz_mul_ui(x,y,(unsigned long)(-z));
    mpz_neg(x,x);
  } else {
    mpz_mul_ui(x,y,(unsigned long)z);
  }
}
#endif
