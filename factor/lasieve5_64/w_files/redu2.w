@* Two dimensional lattice reduction.

@*3 Copying.
Copyright (C) 2001 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@ Note that we dont use lattice reduction to organize sieving with the factor
base primes which are larger than the line length. Instead, we use a method
based on the relation with continued fractions and diophantine approximations.
Therefore, the code contained in this file is not time crtitical.

@(redu2.h@>=
int reduce2(i32_t*,i32_t*,i32_t*,i32_t*,i64_t,i64_t,i64_t,i64_t,double);
int reduce2gmp(i64_t*,i64_t*,i64_t*,i64_t*,mpz_t,mpz_t,double);
extern i32_t n_iter;

@ This is done in the straigthforward way. A basis is given by
the fifth to eighth argument (|a0|,|b0|), (|a1|,|b1|). The
reduced basis is written to the addresses given by the first four arguments.
The ninth argument will often be larger than 1 an punishes large |b|.

As a safeguard against numerical instability, the scalar products
are recalculated from the modified vectors after each reduction step and
not updated in a way which would save multiplications. Note again that this
function is not likely time critical for lattice sieving.

The scalar products of the two vector with themselves are called
|a0sq| and |a1sq|. Their scalar product is called |s|.
@c
#include <math.h>
#include <sys/types.h>
#include "gmp.h"
#include "asm/siever-config.h"
#include "if.h"
#include "gmp-aux.h"
#include "redu2.h"

i32_t n_iter=0;
int
reduce2(i32_t *a0_ptr,i32_t *b0_ptr,i32_t *a1_ptr,i32_t *b1_ptr,
	i64_t a0,i64_t b0,i64_t a1,i64_t b1,double sigma)
{
  double a0sq,a1sq,s;
  i32_t j=0;

  a0sq=((double)a0)*a0;
  a0sq+=sigma*((double)b0)*b0;
  a1sq=((double)a1)*a1;
  a1sq+=sigma*((double)b1)*b1;

  for(;;) {
    s=((double)a0)*a1;
    s+=sigma*((double)b0)*b1;

    if(++j>1024) break; /* n_iter++; */
    if(a0sq<a1sq) {
      i64_t k;

      k=rint(s/a0sq);
      if(k==0) break;
      a1-=k*a0;
      b1-=k*b0;
      a1sq=((double)a1)*a1;
      a1sq+=sigma*((double)b1)*b1;
    } else {
      i64_t k;

      k=rint(s/a1sq);
      if(k==0) break;
      a0-=k*a1;
      b0-=k*b1;
      a0sq=((double)a0)*a0;
      a0sq+=sigma*((double)b0)*b0;
    }
  }

  n_iter += j;

  if(b0<0) {
    b0=-b0;
    a0=-a0;
  }
  if(b1<0) {
    b1=-b1;
    a1=-a1;
  }
  if(a0sq>a1sq) {
    *a0_ptr=(i32_t)a0; if ((i64_t)(*a0_ptr)!=a0) return 1;
    *b0_ptr=(i32_t)b0; if ((i64_t)(*b0_ptr)!=b0) return 1;
    *a1_ptr=(i32_t)a1; if ((i64_t)(*a1_ptr)!=a1) return 1;
    *b1_ptr=(i32_t)b1; if ((i64_t)(*b1_ptr)!=b1) return 1;
  } else {
    *a0_ptr=(i32_t)a1; if ((i64_t)(*a0_ptr)!=a1) return 1;
    *b0_ptr=(i32_t)b1; if ((i64_t)(*b0_ptr)!=b1) return 1;
    *a1_ptr=(i32_t)a0; if ((i64_t)(*a1_ptr)!=a0) return 1;
    *b1_ptr=(i32_t)b0; if ((i64_t)(*b1_ptr)!=b0) return 1;
  }    
  return 0;
}

#ifdef VERY_LARGE_Q
static mpz_t r2_a0, r2_a1, r2_b0, r2_b1;
static mpz_t r2_a0sq, r2_a1sq, r2_s, r2_aux, r2_sigma;
static int reduce2gmp_isinit=0;

/*
CAVE for large SNFS projects: sigma will be rounded to an integer
*/
int
reduce2gmp(i64_t *a0_ptr,i64_t *b0_ptr,i64_t *a1_ptr,i64_t *b1_ptr,
        mpz_t q, mpz_t r,double sigma)
{
  double a0sq,a1sq,s;

  if (!reduce2gmp_isinit) {
    mpz_init(r2_a0);
    mpz_init(r2_b0);
    mpz_init(r2_a1);
    mpz_init(r2_b1);
    mpz_init(r2_a0sq);
    mpz_init(r2_a1sq);
    mpz_init(r2_s);
    mpz_init(r2_aux);
    mpz_init(r2_sigma);
  }
  mpz_set(r2_a0,q); mpz_set_ui(r2_b0,0);
  mpz_set(r2_a1,r); mpz_set_ui(r2_b1,1);
  mpz_set_d(r2_sigma,sigma);
  if (mpz_sgn(r2_sigma)<=0) mpz_set_ui(r2_sigma,1);

  mpz_mul(r2_a0sq,r2_a0,r2_a0);
  mpz_mul(r2_aux,r2_b0,r2_b0);
  mpz_addmul(r2_a0sq,r2_aux,r2_sigma);
  mpz_mul(r2_a1sq,r2_a1,r2_a1);
  mpz_mul(r2_aux,r2_b1,r2_b1);
  mpz_addmul(r2_a1sq,r2_aux,r2_sigma);

  for(;;) {
    mpz_mul(r2_s,r2_a0,r2_a1);
    mpz_mul(r2_aux,r2_b0,r2_b1);
    mpz_addmul(r2_s,r2_aux,r2_sigma);

    n_iter++;
    if(mpz_cmp(r2_a0sq,r2_a1sq)<0) {
      mpz_fdiv_qr(r2_s,r2_aux,r2_s,r2_a0sq);
      mpz_add(r2_aux,r2_aux,r2_aux);
      if (mpz_cmp(r2_aux,r2_a0sq)>0) mpz_add_ui(r2_s,r2_s,1);
      if (mpz_sgn(r2_s)==0) break;
      mpz_submul(r2_a1,r2_s,r2_a0);
      mpz_submul(r2_b1,r2_s,r2_b0);
      mpz_mul(r2_a1sq,r2_a1,r2_a1);
      mpz_mul(r2_aux,r2_b1,r2_b1);
      mpz_addmul(r2_a1sq,r2_aux,r2_sigma);
    } else {
      mpz_fdiv_qr(r2_s,r2_aux,r2_s,r2_a1sq);
      mpz_add(r2_aux,r2_aux,r2_aux);
      if (mpz_cmp(r2_aux,r2_a1sq)>0) mpz_add_ui(r2_s,r2_s,1);
      if (mpz_sgn(r2_s)==0) break;
      mpz_submul(r2_a0,r2_s,r2_a1);
      mpz_submul(r2_b0,r2_s,r2_b1);
      mpz_mul(r2_a0sq,r2_a0,r2_a0);
      mpz_mul(r2_aux,r2_b0,r2_b0);
      mpz_addmul(r2_a0sq,r2_aux,r2_sigma);
    }
  }
  if(mpz_sgn(r2_b0)<0) {
    mpz_neg(r2_b0,r2_b0);
    mpz_neg(r2_a0,r2_a0);
  }
  if(mpz_sgn(r2_b1)<0) {
    mpz_neg(r2_b1,r2_b1);
    mpz_neg(r2_a1,r2_a1);
  }
  if (mpz_cmp(r2_a0sq,r2_a1sq)<=0) {
    mpz_swap(r2_a0,r2_a1);
    mpz_swap(r2_b0,r2_b1);
  }
#if 0
  if (!mpz_sizeinbase(r2_a0,2)>31) return 1;
  if (!mpz_sizeinbase(r2_b0,2)>31) return 1;
  if (!mpz_sizeinbase(r2_a1,2)>31) return 1;
  if (!mpz_sizeinbase(r2_b1,2)>31) return 1;
  *a0_ptr=(i32_t)mpz_get_si(r2_a0);
  *b0_ptr=(i32_t)mpz_get_si(r2_b0);
  *a1_ptr=(i32_t)mpz_get_si(r2_a1);
  *b1_ptr=(i32_t)mpz_get_si(r2_b1);
#else
  if (!mpz_sizeinbase(r2_a0,2)>63) return 1;
  if (!mpz_sizeinbase(r2_b0,2)>63) return 1;
  if (!mpz_sizeinbase(r2_a1,2)>63) return 1;
  if (!mpz_sizeinbase(r2_b1,2)>63) return 1;
  *a0_ptr=(i64_t)mpz_get_sll(r2_a0);
  *b0_ptr=(i64_t)mpz_get_sll(r2_b0);
  *a1_ptr=(i64_t)mpz_get_sll(r2_a1);
  *b1_ptr=(i64_t)mpz_get_sll(r2_b1);
#endif
  return 0;
}
#endif

