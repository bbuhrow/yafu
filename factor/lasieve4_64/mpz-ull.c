/*1:*/
#line 15 "mpz-ull.w"

#include <sys/types.h> 
#include <limits.h> 
typedef unsigned long long ullong;
#include <gmp.h> 
#include "asm/siever-config.h"
#include "gmp-aux.h"
#ifdef ULL_NO_UL
static ulong have_init= 0;
static mpz_t auxz,auxz2;

#define ULLONG_MAX 0xffffffffffffffffULL

void
mpz_ull_init()
{
if(have_init!=0)return;
mpz_init(auxz);
mpz_init(auxz2);
have_init= 1;
}
#endif

/*:1*//*2:*/
#line 40 "mpz-ull.w"

#define BITS_PER_ULONG (sizeof(ulong)*CHAR_BIT)
#ifdef ULL_NO_UL
void mpz_set_ull(mpz_t targ,unsigned long long src)
{
mpz_set_ui(targ,
(ulong)((src&((unsigned long long)ULONG_MAX<<BITS_PER_ULONG))>>BITS_PER_ULONG));
mpz_mul_2exp(targ,targ,BITS_PER_ULONG);
mpz_add_ui(targ,targ,(ulong)(src&ULONG_MAX));
}
#endif

/*:2*//*3:*/
#line 53 "mpz-ull.w"

#ifdef ULL_NO_UL
ullong
mpz_get_ull(mpz_t src)
{
ullong res;

if(sizeof(ullong)==2*sizeof(ulong)){
mpz_fdiv_q_2exp(auxz,src,sizeof(ulong)*CHAR_BIT);
res= mpz_get_ui(auxz);
res<<= sizeof(ulong)*CHAR_BIT;
res|= mpz_get_ui(src);
}else{

if(sizeof(ulong)==sizeof(ullong))return mpz_get_ui(src);
else{
ulong i;
res= mpz_get_ui(src);
mpz_fdiv_q_2exp(auxz,src,CHAR_BIT*sizeof(ulong));
res|= ((ullong)mpz_get_ui(auxz))<<(sizeof(ulong)*CHAR_BIT);
for(i= 2;i*sizeof(ulong)<sizeof(ullong);i++){
mpz_fdiv_q_2exp(auxz,src,CHAR_BIT*sizeof(ulong));
res|= ((ullong)mpz_get_ui(auxz))<<(i*sizeof(ulong)*CHAR_BIT);
}
}
}
return res;
}
#endif

/*:3*//*4:*/
#line 84 "mpz-ull.w"

#ifdef ULL_NO_UL
int
mpz_cmp_ull(mpz_t op1,ullong op2)
{
mpz_set_ull(auxz,op2);
return mpz_cmp(op1,auxz);
}
#endif

/*:4*//*5:*/
#line 95 "mpz-ull.w"

#ifdef ULL_NO_UL
void
mpz_mul_ull(mpz_t rop,mpz_t op1,ullong op2)
{
mpz_set_ull(auxz,op2);
mpz_mul(rop,op1,auxz);
}
#endif

/*:5*//*6:*/
#line 106 "mpz-ull.w"

#ifdef ULL_NO_UL
long long int
mpz_get_sll(mpz_t x)
{
if(mpz_sgn(x)<0){
mpz_neg(auxz2,x);
return-((long long int)mpz_get_ull(auxz2));
}
else return mpz_get_ull(x);
}
#endif

/*:6*//*7:*/
#line 120 "mpz-ull.w"

#ifdef ULL_NO_UL
void
mpz_set_sll(mpz_t x,long long int src)
{
if(src<0){
mpz_set_ull(x,(ullong)(-src));
mpz_neg(x,x);
}else mpz_set_ull(x,(ullong)src);
}
#endif

/*:7*//*8:*/
#line 133 "mpz-ull.w"

#ifdef ULL_NO_UL
int
mpz_fits_ullong_p(mpz_t x)
{
if(mpz_sgn(x)<0)return 0;
if(mpz_sizeinbase(x,2)> CHAR_BIT*sizeof(ullong))return 0;
return 1;
}
#endif

/*:8*//*9:*/
#line 145 "mpz-ull.w"

#ifdef ULL_NO_UL
int
mpz_fits_sllong_p(mpz_t x)
{
if(mpz_sgn(x)> 0)return mpz_get_ull(x)<ULLONG_MAX/2;
else{
mpz_neg(auxz2,x);
return mpz_get_ull(auxz2)<=ULLONG_MAX/2;
}
}
#endif/*:9*/
