//@f mpz int
//@f ullong int

//Combining |mpz| and |ullong|.
//
//Copyright (C) 2002 Jens Franke, T. Kleinjung
//This file is part of gnfs4linux, distributed under the terms of the 
//GNU General Public Licence and WITHOUT ANY WARRANTY.
//
//You should have received a copy of the GNU General Public License along
//with this program; see the file COPYING.  If not, write to the Free
//Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
//02111-1307, USA.

#include <sys/types.h>
#include <limits.h>
typedef unsigned long long ullong;
#include <gmp.h>
#include "asm/siever-config.h"
#include "gmp-aux.h"


#ifdef ULL_NO_UL
static unsigned int have_init=0;
static mpz_t auxz,auxz2;

#define ULLONG_MAX 0xffffffffffffffffULL

void
mpz_ull_init()
{
  if(have_init!=0) return;
  mpz_init(auxz);
  mpz_init(auxz2);
  have_init=1;
}
#endif

#ifndef BITS_PER_GMP_ULONG
#define BITS_PER_GMP_ULONG (sizeof(ulong)*CHAR_BIT)
#endif

#ifdef ULL_NO_UL
//#include <stdio.h>
void mpz_set_ull(mpz_t targ,unsigned long long src)
{
    //printf("sizeof(ulong) = %u\n", sizeof(ulong));
    //printf("sizeof(ullong) = %u\n", sizeof(ullong));
    //printf("sizeof(mp_limb_t) = %u\n", sizeof(mp_limb_t));
    //printf("ulong_max = %lx\n", ULONG_MAX);
    //printf("BITS_PER_GMP_ULONG = %lu\n", BITS_PER_GMP_ULONG);
    // printf("src = %llx\n", src);
    // printf("ulongmax<<BITS_PER_GMP_ULONG = %llx\n", ((unsigned long long)ULONG_MAX << BITS_PER_GMP_ULONG));
    // printf("src&prev = %llx\n", (src & ((unsigned long long)ULONG_MAX << BITS_PER_GMP_ULONG)));
    // printf("prev>>BITS_PER_GMP_ULONG = %lx\n", 
    //     (ulong)((src & ((unsigned long long)ULONG_MAX << BITS_PER_GMP_ULONG)) >> BITS_PER_GMP_ULONG));
    // printf("src&ulong_max = %lx\n", (ulong)(src & ULONG_MAX));
    mpz_set_ui(targ,
	     (ulong)((src&((unsigned long long)ULONG_MAX<<BITS_PER_GMP_ULONG))>>BITS_PER_GMP_ULONG));
    mpz_mul_2exp(targ,targ,BITS_PER_GMP_ULONG);
    mpz_add_ui(targ,targ,(ulong)(src&ULONG_MAX));
    // gmp_printf("result = %Zx\n", targ);
}
#endif

#ifdef ULL_NO_UL
ullong
mpz_get_ull(mpz_t src)
{
  ullong res;
  if(sizeof(ullong)==2*(BITS_PER_GMP_ULONG/CHAR_BIT)) { //sizeof(ulong)) {
    mpz_fdiv_q_2exp(auxz,src, BITS_PER_GMP_ULONG);
    res=mpz_get_ui(auxz);
    res<<= BITS_PER_GMP_ULONG;
    res|=mpz_get_ui(src);
  } else {
    /* CAVE: This has not been checked! */
    if(sizeof(ulong)==sizeof(ullong)) return mpz_get_ui(src);
    else {
      ulong i;
      res=mpz_get_ui(src);
      mpz_fdiv_q_2exp(auxz,src,CHAR_BIT*sizeof(ulong));
      res|=((ullong)mpz_get_ui(auxz))<<(sizeof(ulong)*CHAR_BIT);
      for(i=2;i*sizeof(ulong)<sizeof(ullong);i++) {
	mpz_fdiv_q_2exp(auxz,src,CHAR_BIT*sizeof(ulong));
	res|=((ullong)mpz_get_ui(auxz))<<(i*sizeof(ulong)*CHAR_BIT);
      }
    }
  }
  return res;
}
#endif

#ifdef ULL_NO_UL
int
mpz_cmp_ull(mpz_t op1,ullong op2)
{
  mpz_set_ull(auxz,op2);
  return mpz_cmp(op1,auxz);
}
#endif

#ifdef ULL_NO_UL
void
mpz_add_ull(mpz_t rop, mpz_t op1, ullong op2)
{
    mpz_set_ull(auxz, op2);
    mpz_add(rop, op1, auxz);
}
#endif

#ifdef ULL_NO_UL
void
mpz_tdiv_q_ull(mpz_t rop, mpz_t op1, ullong op2)
{
    mpz_set_ull(auxz, op2);
    mpz_tdiv_q(rop, op1, auxz);
}
#endif

#ifdef ULL_NO_UL
void
mpz_mul_ull(mpz_t rop,mpz_t op1,ullong op2)
{
  mpz_set_ull(auxz,op2);
  mpz_mul(rop,op1,auxz);
}

void
mpz_mul_sll(mpz_t rop,mpz_t op1,long long int op2)
{
  mpz_set_sll(auxz,op2);
  mpz_mul(rop,op1,auxz);
}
#endif

#ifdef ULL_NO_UL
long long int
mpz_get_sll(mpz_t x)
{
  if(mpz_sgn(x)<0) {
    mpz_neg(auxz2,x);
    return -((long long int)mpz_get_ull(auxz2));
  }
  else return mpz_get_ull(x);
}
#endif

#ifdef ULL_NO_UL
void
mpz_set_sll(mpz_t x,long long int src)
{
  if(src<0) {
    mpz_set_ull(x,(ullong)(-src));
    mpz_neg(x,x);
  } else mpz_set_ull(x,(ullong)src);
}
#endif

#ifdef ULL_NO_UL
int
mpz_fits_ullong_p(mpz_t x)
{
  if(mpz_sgn(x)<0) return 0;
  if(mpz_sizeinbase(x,2)>CHAR_BIT*sizeof(ullong)) return 0;
  return 1;
}
#endif

#ifdef ULL_NO_UL
int
mpz_fits_sllong_p(mpz_t x)
{
  if(mpz_sgn(x)>0) return mpz_get_ull(x)<ULLONG_MAX/2;
  else {
    mpz_neg(auxz2,x);
    return mpz_get_ull(auxz2)<=ULLONG_MAX/2;
  }
}
#endif
