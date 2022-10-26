#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>
#include "siever-config.h"
#include "../if.h"

#include "montgomery_mul.h"

#define uchar   unsigned char

ulong montgomery_inv_n;
ulong *montgomery_modulo_n;
ulong montgomery_modulo_R2[NMAX_ULONGS], montgomery_modulo_R4[NMAX_ULONGS];
size_t montgomery_ulongs;
mpz_t montgomery_gmp_help;
static int montgomery_multiplication_is_init=0;


/* function pointers */
void ASM_ATTR (*asm_mulmod)(ulong *,ulong *,ulong *)=NULL;
void ASM_ATTR (*asm_zero)(ulong *)=NULL;
void ASM_ATTR (*asm_copy)(ulong *,ulong *)=NULL;
void ASM_ATTR (*asm_half)(ulong *)=NULL;
void ASM_ATTR (*asm_sub)(ulong *,ulong *,ulong *)=NULL;
void ASM_ATTR (*asm_add2)(ulong *,ulong *)=NULL;

void (*asm_sub_n)(ulong *,ulong *)=NULL;
void (*asm_squmod)(ulong *,ulong *)=NULL;
void (*asm_diff)(ulong *,ulong *,ulong *)=NULL;
void (*asm_add2_ui)(ulong *,ulong)=NULL;
int (*asm_inv)(ulong *,ulong *)=NULL;



/* Moved to header file
extern void ASM_ATTR asm_mulm64(ulong *,ulong *,ulong *);
extern void ASM_ATTR asm_zero64(ulong *);
extern void ASM_ATTR asm_sub64_3(ulong *,ulong *,ulong *);
extern void ASM_ATTR asm_copy64(ulong *,ulong *);
extern void ASM_ATTR asm_half64(ulong *);
extern void ASM_ATTR asm_add64(ulong *,ulong *);

extern void asm_sub_n64(ulong *,ulong *);
extern void asm_add64_ui(ulong *,ulong);
extern void asm_sqm64(ulong *,ulong *);
extern void asm_diff64(ulong *,ulong *,ulong *);
extern int asm_inv64(ulong *,ulong *);
*/

#ifdef ULONG_HAS_32BIT
extern void asm_mulm96(ulong *,ulong *,ulong *);
extern void asm_sqm96(ulong *,ulong *);
extern void asm_add96(ulong *,ulong *);
extern void asm_diff96(ulong *,ulong *,ulong *);
extern void asm_sub96_3(ulong *,ulong *,ulong *);
extern void asm_add96_ui(ulong *,ulong);
extern void asm_copy96(ulong *,ulong *);
extern void asm_zero96(ulong *);
extern void asm_sub_n96(ulong *,ulong *);
extern void asm_half96(ulong *);
extern int asm_inv96(ulong *,ulong *);
#endif

/* Moved to header file
extern void ASM_ATTR asm_sub128_3(ulong *,ulong *,ulong *);
extern void ASM_ATTR asm_zero128(ulong *);
extern void ASM_ATTR asm_copy128(ulong *,ulong *);
extern void ASM_ATTR asm_half128(ulong *);
extern void ASM_ATTR asm_mulm128(ulong *,ulong *,ulong *);
extern void ASM_ATTR asm_add128(ulong *,ulong *);

extern void asm_sub_n128(ulong *,ulong *);
void asm_sqm128(ulong *,ulong *);
extern void asm_diff128(ulong *,ulong *,ulong *);
extern void asm_add128_ui(ulong *,ulong);
extern int asm_inv128(ulong *,ulong *);
*/

#ifdef ULONG_HAS_32BIT
extern void asm_mulm160(ulong *,ulong *,ulong *);
void asm_sqm160(ulong *,ulong *);
extern void asm_diff160(ulong *,ulong *,ulong *);
extern void asm_add160(ulong *,ulong *);
extern void asm_sub160_3(ulong *,ulong *,ulong *);
extern void asm_add160_ui(ulong *,ulong);
extern void asm_zero160(ulong *);
extern void asm_copy160(ulong *,ulong *);
extern void asm_sub_n160(ulong *,ulong *);
extern void asm_half160(ulong *);
extern int asm_inv160(ulong *,ulong *);
#endif

/* Moved to header file
extern void ASM_ATTR asm_zero192(ulong *);
extern void ASM_ATTR asm_mulm192(ulong *,ulong *,ulong *);
extern void ASM_ATTR asm_sub192_3(ulong *,ulong *,ulong *);
extern void ASM_ATTR asm_add192(ulong *,ulong *);
extern void ASM_ATTR asm_copy192(ulong *,ulong *);
extern void ASM_ATTR asm_half192(ulong *);

extern void asm_sub_n192(ulong *,ulong *);
void asm_sqm192(ulong *,ulong *);
extern void asm_diff192(ulong *,ulong *,ulong *);
extern void asm_add192_ui(ulong *,ulong);
extern int asm_inv192(ulong *,ulong *);
*/


extern int asm_invert(ulong *,ulong *);



ulong montgomery_inv_add_table[32][8];  /* assuming NMAX_ULONGS<=8 */
ulong montgomery_inv_shift_table[16]={
4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
};

/* ----------------------------------------- */

void asm_sqm128(ulong *x,ulong *y)
{
  asm_mulm128(x,y,y);
}

#ifdef ULONG_HAS_32BIT
void asm_sqm160(ulong *x,ulong *y)
{
  asm_mulm160(x,y,y);
}
#endif

void asm_sqm192(ulong *x,ulong *y)
{
  asm_mulm192(x,y,y);
}

#ifndef HAVE_ASM_INV
int asm_inv64(ulong *res, ulong *b)
{
  return asm_invert(res,b);
}

#ifdef ULONG_HAS_32BIT
int asm_inv96(ulong *res, ulong *b)
{
  return asm_invert(res,b);
}
#endif

int asm_inv128(ulong *res, ulong *b)
{
  return asm_invert(res,b);
}

#ifdef ULONG_HAS_32BIT
int asm_inv160(ulong *res, ulong *b)
{
  return asm_invert(res,b);
}
#endif

int asm_inv192(ulong *res, ulong *b)
{
  return asm_invert(res,b);
}

#endif

/* ----------------------------------------- */

void init_montgomery_R2(mpz_t n)
{
  size_t i;

  mpz_set_ui(montgomery_gmp_help,1);
  mpz_mul_2exp(montgomery_gmp_help,montgomery_gmp_help,2*mp_bits_per_limb*montgomery_ulongs);
  mpz_fdiv_r(montgomery_gmp_help,montgomery_gmp_help,n);
  montgomery_modulo_R2[0]=mpz_get_ui(montgomery_gmp_help);
  for (i=1; i<montgomery_ulongs; i++) {
    mpz_fdiv_q_2exp(montgomery_gmp_help,montgomery_gmp_help,mp_bits_per_limb);
    montgomery_modulo_R2[i]=mpz_get_ui(montgomery_gmp_help);
  }
}


uchar montgomery_inv_table[128]={
1, 171, 205, 183, 57, 163, 197, 239,
241, 27, 61, 167, 41, 19, 53, 223,
225, 139, 173, 151, 25, 131, 165, 207,
209, 251, 29, 135, 9, 243, 21, 191,
193, 107, 141, 119, 249, 99, 133, 175,
177, 219, 253, 103, 233, 211, 245, 159,
161, 75, 109, 87, 217, 67, 101, 143,
145, 187, 221, 71, 201, 179, 213, 127,
129, 43, 77, 55, 185, 35, 69, 111,
113, 155, 189, 39, 169, 147, 181, 95,
97, 11, 45, 23, 153, 3, 37, 79,
81, 123, 157, 7, 137, 115, 149, 63,
65, 235, 13, 247, 121, 227, 5, 47,
49, 91, 125, 231, 105, 83, 117, 31,
33, 203, 237, 215, 89, 195, 229, 15,
17, 59, 93, 199, 73, 51, 85, 255 };

ulong montgomery_inverse()
{
  ulong inv, h, a;

  a=montgomery_modulo_n[0];
  inv=(ulong)montgomery_inv_table[(a&0xff)>>1];
  h=(ulong)a*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=(ulong)a*inv; h&=0xffff0000;
  h*=inv; inv-=h;
  if (mp_bits_per_limb==64) {
    h=(ulong)a*inv; h&=0xffffffff00000000ULL;
    h*=inv; inv-=h;
  }
  if (inv*(ulong)a!=1ULL) Schlendrian("montgomery_inverse\n");
  return (-inv);
}


void asm_invert_init()
{
  u32_t i, j, r, s;

  asm_zero(montgomery_inv_add_table[1]);
  for (i=1,r=1,s=2; i<=4; i++,r<<=1,s<<=1) {
    for (j=0; j<r; j++)
      asm_copy(montgomery_inv_add_table[s+2*j],montgomery_inv_add_table[r+j]);
    for (j=1; j<r; j+=2) {
      asm_copy(montgomery_inv_add_table[s+j],montgomery_inv_add_table[r+j]);
      asm_half(montgomery_inv_add_table[s+j]);
    }
    for (; j<s; j+=2) {
      asm_copy(montgomery_inv_add_table[s+j],montgomery_inv_add_table[j]);
      asm_add2_ui(montgomery_inv_add_table[s+j],1);
      asm_half(montgomery_inv_add_table[s+j]);
    }
  }
}


int set_montgomery_multiplication(mpz_t n)
{
  size_t new, old;
  int j;

  if (!montgomery_multiplication_is_init) init_montgomery_multiplication();
  old=montgomery_ulongs;
  new=mpz_size(n);
  if (new>NMAX_ULONGS) return 0;
  for (j=1; j<NMAX_ULONGS; j++) montgomery_modulo_n[j]=0;
  montgomery_modulo_n[0]=mpz_get_ui(n);
  montgomery_ulongs=1;
  if (new>1) mpz_fdiv_q_2exp(montgomery_gmp_help,n,mp_bits_per_limb);
  while (montgomery_ulongs<new) {
    montgomery_modulo_n[montgomery_ulongs]=mpz_get_ui(montgomery_gmp_help);
    mpz_fdiv_q_2exp(montgomery_gmp_help,montgomery_gmp_help,mp_bits_per_limb);
    montgomery_ulongs++;
  }
/* modular multiplication may give wrong results if most significant limb is
   2^mp_bits_per_limb-1, therefore: */
  if (montgomery_modulo_n[montgomery_ulongs-1]+1==0) {
    if (montgomery_ulongs>=NMAX_ULONGS) return 0;
    montgomery_modulo_n[montgomery_ulongs]=0;
    montgomery_ulongs++;
  }
  if ((mp_bits_per_limb==32) && (montgomery_ulongs<2)) montgomery_ulongs=2;
       /* have no 32 bit functions */
  if (montgomery_ulongs>NMAX_ULONGS)
    complain("set_montgomery_multiplication\n");
  if (montgomery_ulongs*mp_bits_per_limb>192) {
    montgomery_ulongs=old;
    return 0;
  }
  if (old!=montgomery_ulongs) {
    if (montgomery_ulongs*mp_bits_per_limb==64) {
      asm_mulmod=asm_mulm64;
      asm_squmod=asm_sqm64;
      asm_add2=asm_add64;
      asm_diff=asm_diff64;
      asm_sub=asm_sub64_3;
      asm_add2_ui=asm_add64_ui;
      asm_zero=asm_zero64;
      asm_copy=asm_copy64;
      asm_sub_n=asm_sub_n64;
      asm_half=asm_half64;
      asm_inv=asm_inv64;
#ifdef ULONG_HAS_32BIT
    } else if (montgomery_ulongs*mp_bits_per_limb==96) {
      asm_mulmod=asm_mulm96;
      asm_squmod=asm_sqm96;
      asm_add2=asm_add96;
      asm_diff=asm_diff96;
      asm_sub=asm_sub96_3;
      asm_add2_ui=asm_add96_ui;
      asm_zero=asm_zero96;
      asm_copy=asm_copy96;
      asm_sub_n=asm_sub_n96;
      asm_half=asm_half96;
      asm_inv=asm_inv96;
#endif
    } else if (montgomery_ulongs*mp_bits_per_limb==128) {
      asm_mulmod=asm_mulm128;
      asm_squmod=asm_sqm128;
      asm_add2=asm_add128;
      asm_diff=asm_diff128;
      asm_sub=asm_sub128_3;
      asm_add2_ui=asm_add128_ui;
      asm_zero=asm_zero128;
      asm_copy=asm_copy128;
      asm_sub_n=asm_sub_n128;
      asm_half=asm_half128;
      asm_inv=asm_inv128;
#ifdef ULONG_HAS_32BIT
    } else if (montgomery_ulongs*mp_bits_per_limb==160) {
      asm_mulmod=asm_mulm160;
      asm_squmod=asm_sqm160;
      asm_add2=asm_add160;
      asm_diff=asm_diff160;
      asm_sub=asm_sub160_3;
      asm_add2_ui=asm_add160_ui;
      asm_zero=asm_zero160;
      asm_copy=asm_copy160;
      asm_sub_n=asm_sub_n160;
      asm_half=asm_half160;
      asm_inv=asm_inv160;
#endif
    } else if (montgomery_ulongs*mp_bits_per_limb==192) {
      asm_mulmod=asm_mulm192;
      asm_squmod=asm_sqm192;
      asm_add2=asm_add192;
      asm_diff=asm_diff192;
      asm_sub=asm_sub192_3;
      asm_add2_ui=asm_add192_ui;
      asm_zero=asm_zero192;
      asm_copy=asm_copy192;
      asm_sub_n=asm_sub_n192;
      asm_half=asm_half192;
      asm_inv=asm_inv192;
    } else return 0;
  }
  montgomery_inv_n=montgomery_inverse();
  if (montgomery_inv_n*montgomery_modulo_n[0]+1) {
    fprintf(stderr,"init_montgomery_multiplication failed %lu %lu\n",montgomery_inv_n,montgomery_modulo_n[0]);
    return 0;
  }
  init_montgomery_R2(n);
  asm_mulmod(montgomery_modulo_R4,montgomery_modulo_R2,montgomery_modulo_R2);
  asm_invert_init();
  return 1;
}


void init_montgomery_multiplication()
{
  mpz_init(montgomery_gmp_help);
  montgomery_modulo_n=(ulong *)xmalloc(NMAX_ULONGS*sizeof(ulong));
  montgomery_ulongs=0;
  montgomery_multiplication_is_init=1;
}

