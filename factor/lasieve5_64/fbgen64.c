/* Originally written by F. Bahr. For use in the quality estimator for GNFS
   polynomials, adjustments (sometimes crude adjustments) have been made by
   J. Franke. */

/*
  Copyright (C) 2001 Friedrich Bahr

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

/* very crude adjustment for 63 bit by Thorsten Kleinjung, very slow */

#include <sys/types.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>

#include "asm/siever-config.h"
#include "if.h"
#include "gmp-aux.h"
#include "fbgen64.h"

static u64_t mod64;
static u64_t modulo64hbit;
static i64_t modulo64bit[64],modulo64ebits;

static u64_t * polr,t,polr_alloc=0;



static u64_t ma64(u64_t x, u64_t y)
{
  u64_t z;

  z=x+y; if (z>=mod64) z-=mod64;
  return z;
}

static u64_t ms64(u64_t x, u64_t y)
{
  u64_t z;

  if (y>=mod64) Schlendrian("ms64\n");
  z=x+(mod64-y); if (z>=mod64) z-=mod64;
  return z;
}

static mpz_t mm64_aux1, mm64_aux2;
static int mm64_init=0;

static u64_t mm64(u64_t x, u64_t y)
{
  if (!mm64_init) {
    mpz_init(mm64_aux1);
    mpz_init(mm64_aux2);
    mm64_init=1;
  }
  mpz_set_ull(mm64_aux1,x);
  mpz_set_ull(mm64_aux2,y);
  mpz_mul(mm64_aux1,mm64_aux1,mm64_aux2);
  mpz_set_ull(mm64_aux2,mod64);
  mpz_mod(mm64_aux1,mm64_aux1,mm64_aux2);
  return mpz_get_ull(mm64_aux1);
}


static u64_t modsq64(u64_t x)
{
  return mm64(x,x);
}

/* uses same auxiliary variables as mm64 */
static u64_t mi64(u64_t x)
{
  if (!mm64_init) {
    mpz_init(mm64_aux1);
    mpz_init(mm64_aux2);
    mm64_init=1;
  }
  mpz_set_ull(mm64_aux1,x);
  mpz_set_ull(mm64_aux2,mod64);
  if (!mpz_invert(mm64_aux1,mm64_aux1,mm64_aux2))
    Schlendrian("mi64: inverse does not exist\n");
  return mpz_get_ull(mm64_aux1);
}

static mpz_t fdiv_ull_aux;
static int fdiv_ull_init=0;

static u64_t mpz_fdiv_ull(mpz_t op, u64_t div)
{
  if (!fdiv_ull_init) {
    mpz_init(fdiv_ull_aux);
    fdiv_ull_init=1;
  }
  mpz_set_ull(fdiv_ull_aux,div);
  mpz_fdiv_r(fdiv_ull_aux,op,fdiv_ull_aux);
  return mpz_get_ull(fdiv_ull_aux);
}


static inline u64_t polmodsq64(u64_t * P,u64_t dP){
   u64_t i,j,a,m;
   m=P[dP];
   P[dP<<1]=modsq64(m);
   for(i=0;i<dP;i++){
     a=mm64(m,P[i]);
     P[dP+i]=ma64(a,a);
   }
   for(j=dP;j>1;){
     m=P[--j];
     P[j<<1]=ma64(P[j<<1],modsq64(m));
     a=mm64(m,P[0]);
     P[j]=ma64(a,a);
     for(i=1;i<j;i++){
       a=mm64(m,P[i]);
       P[j+i]=ma64(P[j+i],ma64(a,a));
     }
   }
   if(dP!=0)
     P[0]=modsq64(P[0]);
   dP=dP<<1;
   return dP;
}

static u64_t poldivmod64(u64_t *P,u64_t dP,u64_t *D,u64_t dD){
   u64_t i,a,b,j;
   /* D!=0 ! vorausgesetzt */
   while(dP>=dD){
     a=P[dP],b=D[dD];
     for(i=j=dP-dD;i<dP;i++)
       P[i]=ms64(mm64(P[i],b),mm64(D[i-j],a));
     for(i=0;i<j;i++)
       P[i]=mm64(P[i],b);
     for(dP--;P[dP]==0 && dP>0;dP--);
   }
   return dP;
}

static u64_t polcpymod64(u64_t * Q,u64_t * P,u64_t dP){
   memcpy(Q, P, (dP+1)*sizeof(u64_t));
   return dP;
}

static u64_t polgcdmod64(u64_t * P,u64_t dP, u64_t * Q, u64_t dQ){
   /* hier wird dP > dQ angenommen */
   /* stets ist Q das Ergebnis */
   u64_t * S,d,*R=Q;
   while(dQ){
     d=poldivmod64(P,dP,Q,dQ);
     S=P,P=Q,Q=S,dP=dQ,dQ=d;
   }
   if(Q[0]==0){
     if(P!=R)
       polcpymod64(R,P,dP);
     return dP;
   }
   R[0]=1;
   return 0;
}

static inline u64_t poldivnormmod64(u64_t *P,u64_t dP,u64_t *D, u64_t dD){ /* D 
normiert */
   u64_t i,a,j;
   while(dP>=dD){
     a=P[dP];
     for(i=j=dP-dD;i<dP;i++)
       P[i]=ms64(P[i],mm64(D[i-j],a));
     for(dP--;P[dP]==0 && dP>0;dP--);
   }
   return dP;
}

static u64_t polnormmod64(u64_t *P,u64_t dP){
   u64_t a,i;
   if(P[dP]!=1 && P[dP]!=0){
     a=mi64(P[dP]);
     for(i=0;i<dP;i++)
       P[i]=mm64(P[i],a);
     P[dP]=1;
   }
   return dP;
}

/* CAVE Was passiert, falls die Eingabe durch |mod64| teilbar ist? */

static u64_t polredmod64(u64_t * T, mpz_t *A, u64_t dA){
   u64_t i;
   for(i=0;i<=dA;i++)
     T[i] = mpz_fdiv_ull(A[i],mod64);
   while(T[dA]==0) dA--;
   return polnormmod64(T,dA);;
  }

static u64_t polmod64(u64_t * T, mpz_t *A, u64_t dA){
  u64_t i;
  for(i=0;i<=dA;i++)
    T[i] = mpz_fdiv_ull(A[i],mod64);
  while(dA>0 && T[dA]==0) dA--;
  return dA;
}

static u64_t polXpotmodPmod64(u64_t * Q,u64_t * P, u64_t dP){
  u64_t dQ=modulo64bit[modulo64hbit],k=0;
   /* Q = X hoch (mod64-1)/2 - 1 modulo P */
  Q+=modulo64ebits;
  for(;k<dQ;k++)
    Q[k]=0;
  Q[dQ]=1;
  if (dQ>=dP)
    dQ=poldivnormmod64(Q,dQ,P,dP);
  for(k=modulo64hbit;k>0;){
    k--;
    dQ=polmodsq64(Q,dQ);
    if(modulo64bit[k]){ /* Multiplikation mit X */
      Q--;
      dQ++;
      Q[0]=0;
    }
    if (dQ>=dP)
      dQ=poldivnormmod64(Q,dQ,P,dP);
  }
  return dQ;
}

static void polchcomod64(u64_t *P,u64_t dP){
  u64_t i=dP,j;
  for(;i>0;){
    for(i--,j=i;j<dP;j++)
      P[j]=ma64(P[j],P[j+1]);
  }
}

static u64_t polrootrecmod64(u64_t * T, u64_t dT, u64_t * r, u64_t * Q){
   u64_t * P, dP, dQ,i=0,b;
   /* Ersetze T(X) durch T(X+1) */
   if(dT>1){
     polchcomod64(T,dT);
     dQ=polXpotmodPmod64(Q,T,dT);
     Q[0]=ma64(Q[0],1);
     P=Q+dT+1;
     dP=polcpymod64(P,T,dT);
     dQ=polgcdmod64(P,dP,Q,dQ);
     b=r[0];
     r[0]=ma64(b,1);
     if(dQ>0 && dQ<dT){
       dQ=polnormmod64(Q,dQ);
       dP=polcpymod64(P,T,dT);
       dP=poldivnormmod64(P,dP,Q,dQ);
       dP=dT-dQ,P+=dQ;
       i=polrootrecmod64(Q,dQ,r,P+dP+1);  /*i=dQ*/
       r[i]=ma64(b,1);
       i+=polrootrecmod64(P,dP,r+i,P+dP+1); /*i=dT*/
       return i;
     }
     else
       return polrootrecmod64(T,dT,r,Q);
   }
   else {
     r[0]=ms64(r[0],T[0]);
     return 1;
   }
}

static u64_t polvalmod64(u64_t * P, u64_t dP, u64_t a){
  u64_t i=dP,r=P[dP];
  for(;i>0;)
    r=ma64(P[--i],mm64(a,r));
  return r;
}

static void polprintmod64(u64_t *P, u64_t dP, char * c){
  u64_t i;
/* SMJS
  printf("\n%s %llu\n",c,dP);
  for(i=0;i<=dP;i++)
    printf("%llu ",P[i]);
*/
  printf("\n%s "UL_FMTSTR"\n",c,dP);
  for(i=0;i<=dP;i++)
    printf(UL_FMTSTR" ",P[i]);
  printf("\n");
}

static void mod64multab(dT){
  u64_t i,c;
  modulo64ebits=0;
  for(i=0,c=mod64;(c>>=1)>dT;i++){
    if(c&1){
      modulo64bit[i]=1;
      modulo64ebits++;
    } else modulo64bit[i]=0;
  }
  modulo64bit[i]=c;
  modulo64hbit=i;/* Anzahl der notwendigen Quadratbildungen */
}


static u64_t polrootmod64(mpz_t * A, u64_t dT,u64_t *r){
   u64_t i=0,j,*P,dP,*Q,dQ,*S,dS,*T=polr;
   dT=polredmod64(T,A,dT);
   if(T[0]==0){  /* Nullstelle 0 */
     for(i=1;T[i]==0 && i <dT;i++);
     T+=i;dT-=i;
     r[0]=0,i=1;
   }
   if(dT==0)
     return i;
   mod64multab(dT);
   /* Zerlegung Q = (P,X^(mod64-1)/2+1); S = (P,X^(mod64-1)/2+1)*/
   P=T+dT+1;
   dP=polcpymod64(P,T,dT);
   Q=P+dP+1;
   dQ=polXpotmodPmod64(Q,P,dP);
   S=Q+dP+1;/* die ggT-Berechnung kann P ergeben */
   dS=polcpymod64(S,Q,dQ);
   Q[0]=ma64(Q[0],1);
   dQ=polgcdmod64(P,dP,Q,dQ);
   dP=polcpymod64(P,T,dT);
   S[0]=ms64(S[0],1);
   dS=polgcdmod64(P,dP,S,dS);
   if(dQ>0){
     dQ=polnormmod64(Q,dQ);
     r[i]=0;
     j=polrootrecmod64(Q,dQ,r+i,S+dS+1);  /*j=dQ*/
     // SMJSif(j!=dQ) fprintf(stderr,"Falsche Nullstellenanzahl Q, P %llu\n",mod64);
     if(j!=dQ) fprintf(stderr,"Falsche Nullstellenanzahl Q, P "UL_FMTSTR"\n",mod64);
     i+=j;
   }
   if(dS>0){
     dS=polnormmod64(S,dS);
     r[i]=0;
     j=polrootrecmod64(S,dS,r+i,S+dS+1);  /*j=dS*/
     // SMJS if(j!=dS) fprintf(stderr,"Falsche Nullstellenanzahl S, P %llu\n",mod64);
     if(j!=dS) fprintf(stderr,"Falsche Nullstellenanzahl S, P "UL_FMTSTR"\n",mod64);
     i+=j;
   }
   qsort(r,i,sizeof(u64_t),u64_cmp012);
   return i;
}

u32_t
root_finder64(u64_t *root_buf,mpz_t *A,u32_t adeg,u64_t p)
{
  u32_t res;
  if(polr_alloc<300*adeg*adeg) {
    /* CAVE falls adeg==0? */
    if(polr_alloc>0) free(polr);
    polr_alloc=300*adeg*adeg;
    polr=xmalloc(polr_alloc*sizeof(*polr));
  }
  if(p==2) {
    u32_t i,pv;

    res=0;
    if(mpz_fdiv_ull(A[0],p)==0) root_buf[res++]=0;
    for(i=0,pv=0;i<=adeg;i++)
      if(mpz_fdiv_ull(A[i],p)!=0) pv^=1;
    if(pv==0) root_buf[res++]=1;
    if(mpz_fdiv_ull(A[adeg],p)==0)
      root_buf[res++]=2;
    return res;
  }
  mod64=p;
  res=(u32_t)polrootmod64(A,(u64_t)adeg,root_buf);
  /* CAVE: Diese Division erfolgt wohl ein zweites Mal. */
  if(mpz_fdiv_ull(A[adeg],p)==0) root_buf[res++]=p;
  return res;
}
