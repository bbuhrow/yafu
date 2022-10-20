/*2:*/
#line 50 "gmp-aux.w"

#include <stdlib.h> 
#include <sys/types.h> 
#include <string.h> 
#include <gmp.h> 

#include "asm/siever-config.h"
#include "if.h"
#include "gmp-aux.h"














/*:2*//*3:*/
#line 82 "gmp-aux.w"

int
string2mpz(mpz_t rop,char*x,int base)
{
size_t l;
char*y;
int rv;

x+= strspn(x," \t+");
if(strlen(x)==0)mpz_set_ui(rop,0);
y= strdup(x);
for(l= strlen(y)-1;y[l]=='\n';l--){
y[l]= '\0';
if(l==0)break;
}
rv= mpz_set_str(rop,y,base);
free(y);
return rv;
}

/*:3*//*4:*/
#line 103 "gmp-aux.w"

#ifdef NEED_MPZ_MUL_SI
void mpz_mul_si(mpz_t x,mpz_t y,long int z)
{
if(z<0){
mpz_mul_ui(x,y,(unsigned long)(-z));
mpz_neg(x,x);
}else{
mpz_mul_ui(x,y,(unsigned long)z);
}
}
#endif/*:4*/
