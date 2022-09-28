/*1:*/
#line 5 "mpz-trialdiv.w"

#include <stdlib.h> 
#include <gmp.h> 

#include "siever-config.h"
#include "../if.h"

mp_limb_t mpz_asm_td(mp_limb_t,mp_limb_t,mp_limb_t*,mp_size_t);

static mp_limb_t*lbuf,*mibuf;
static size_t lbuf_alloc= 0,mibuf_alloc= 0;
extern unsigned char mpqs_256_inv_table[128];

u32_t*
mpz_trialdiv(mpz_t N,u32_t*pbuf,u32_t ncp,char*errmsg)
{
u32_t np,np1,i,e2,nlimb;

if(mpz_sgn(N)<=0){
if(mpz_sgn(N)==0)return pbuf;
mpz_neg(N,N);
}

if(mibuf_alloc<ncp){
if(mibuf_alloc> 0)free(mibuf);
mibuf_alloc= ncp;
mibuf= xmalloc(mibuf_alloc*sizeof(*mibuf));
}

e2= 0;
while((mpz_get_ui(N)%2)==0){
mpz_fdiv_q_2exp(N,N,1);
e2++;
}

nlimb= N[0]._mp_size;
if(lbuf_alloc<nlimb){
if(lbuf_alloc> 0)free(lbuf);
lbuf_alloc= nlimb;
lbuf= xmalloc(lbuf_alloc*sizeof(*lbuf));
}

memcpy(lbuf,N[0]._mp_d,nlimb*sizeof(*lbuf));
np= 0;
for(i= 0;i<ncp;i++){
mp_limb_t x,p,r;

p= pbuf[i];
x= mpqs_256_inv_table[(p&255)/2];
x= 2*x-p*x*x;
x= 2*x-p*x*x;
x= 2*x-p*x*x;
r= mpz_asm_td(p,x,lbuf,nlimb);
if(r!=0){
if(errmsg!=NULL)
Schlendrian("%s : %u does not divide\n",errmsg,pbuf[i]);
memcpy(lbuf,N[0]._mp_d,nlimb*sizeof(*lbuf));
continue;
}
while(lbuf[nlimb-1]==0)nlimb--;
if(errmsg==NULL)memcpy(N[0]._mp_d,lbuf,nlimb*sizeof(*lbuf));
mibuf[np]= x;
pbuf[np++]= p;
}
np1= np;
if(errmsg!=NULL)memcpy(N[0]._mp_d,lbuf,nlimb*sizeof(*lbuf));
for(i= 0;i<np1;i++){
mp_limb_t x,p,r;

p= pbuf[i];
x= mibuf[i];
for(;;){
r= mpz_asm_td(p,x,lbuf,nlimb);
if(r!=0){
memcpy(lbuf,N[0]._mp_d,nlimb*sizeof(*lbuf));
break;
}
while(lbuf[nlimb-1]==0)nlimb--;
memcpy(N[0]._mp_d,lbuf,nlimb*sizeof(*lbuf));
pbuf[np++]= p;
}
}
N[0]._mp_size= nlimb;
for(i= 0;i<e2;i++)
pbuf[np++]= 2;
return pbuf+np;
}/*:1*/
