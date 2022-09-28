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
#include "asm/32bit.h"
#include "primgen32.h"
#include "fbgen.h"

volatile u32_t modulo32;
u32_t modulo32hbit;
i32_t modulo32bit[32],modulo32ebits;

u32_t * polr,t,polr_alloc=0;

inline u32_t polmodsq32(u32_t * P,u32_t dP){
   u32_t i,j,a,m;
   m=P[dP];
   P[dP<<1]=modsq32(m);
   for(i=0;i<dP;i++){
     a=modmul32(m,P[i]);
     P[dP+i]=modadd32(a,a);
   }
   for(j=dP;j>1;){
     m=P[--j];
     P[j<<1]=modadd32(P[j<<1],modsq32(m));
     a=modmul32(m,P[0]);
     P[j]=modadd32(a,a);
     for(i=1;i<j;i++){
       a=modmul32(m,P[i]);
       P[j+i]=modadd32(P[j+i],modadd32(a,a));
     }
   }
   if(dP!=0)
     P[0]=modsq32(P[0]);
   dP=dP<<1;
   return dP;
}

u32_t poldivmod32(u32_t *P,u32_t dP,u32_t *D,u32_t dD){
   u32_t i,a,b,j;
   /* D!=0 ! vorausgesetzt */
   while(dP>=dD){
     a=P[dP],b=D[dD];
     for(i=j=dP-dD;i<dP;i++)
       P[i]=modsub32(modmul32(P[i],b),modmul32(D[i-j],a));
     for(i=0;i<j;i++)
       P[i]=modmul32(P[i],b);
     for(dP--;P[dP]==0 && dP>0;dP--);
   }
   return dP;
}

u32_t polcpymod32(u32_t * Q,u32_t * P,u32_t dP){
   memcpy(Q, P, (dP+1)*sizeof(u32_t));
   return dP;
}

u32_t polgcdmod32(u32_t * P,u32_t dP, u32_t * Q, u32_t dQ){
   /* hier wird dP > dQ angenommen */
   /* stets ist Q das Ergebnis */
   u32_t * S,d,*R=Q;
   while(dQ){
     d=poldivmod32(P,dP,Q,dQ);
     S=P,P=Q,Q=S,dP=dQ,dQ=d;
   }
   if(Q[0]==0){
     if(P!=R)
       polcpymod32(R,P,dP);
     return dP;
   }
   R[0]=1;
   return 0;
}

inline u32_t poldivnormmod32(u32_t *P,u32_t dP,u32_t *D, u32_t dD){ /* D 
normiert */
   u32_t i,a,j;
   while(dP>=dD){
     a=P[dP];
     for(i=j=dP-dD;i<dP;i++)
       P[i]=modsub32(P[i],modmul32(D[i-j],a));
     for(dP--;P[dP]==0 && dP>0;dP--);
   }
   return dP;
}

u32_t polnormmod32(u32_t *P,u32_t dP){
   u32_t a,i;
   if(P[dP]!=1 && P[dP]!=0){
     a=modinv32(P[dP]);
     for(i=0;i<dP;i++)
       P[i]=modmul32(P[i],a);
     P[dP]=1;
   }
   return dP;
}

/* CAVE Was passiert, falls die Eingabe durch |modulo32| teilbar ist? */

u32_t polredmod32(u32_t * T, mpz_t *A, u32_t dA){
   u32_t i;
   for(i=0;i<=dA;i++)
     T[i] = mpz_fdiv_ui(A[i],modulo32);
   while(T[dA]==0) dA--;
   return polnormmod32(T,dA);;
  }

u32_t polmod32(u32_t * T, mpz_t *A, u32_t dA){
   u32_t i;
   for(i=0;i<=dA;i++)
     T[i] = mpz_fdiv_ui(A[i],modulo32);
   while(dA>0 && T[dA]==0) dA--;
   return dA;
  }

u32_t polXpotmodPmod32(u32_t * Q,u32_t * P, u32_t dP){
   u32_t dQ=modulo32bit[modulo32hbit],k=0;
    /* Q = X hoch (modulo32-1)/2 - 1 modulo P */
   Q+=modulo32ebits;
   for(;k<dQ;k++)
     Q[k]=0;
   Q[dQ]=1;
   if (dQ>=dP)
     dQ=poldivnormmod32(Q,dQ,P,dP);
   for(k=modulo32hbit;k>0;){
     k--;
     dQ=polmodsq32(Q,dQ);
     if(modulo32bit[k]){ /* Multiplikation mit X */
       Q--;
       dQ++;
       Q[0]=0;
     }
     if (dQ>=dP)
       dQ=poldivnormmod32(Q,dQ,P,dP);
   }
   return dQ;
}

void polchcomod32(u32_t *P,u32_t dP){
   u32_t i=dP,j;
   for(;i>0;){
     for(i--,j=i;j<dP;j++)
       P[j]=modadd32(P[j],P[j+1]);
   }
}

u32_t polrootrecmod32(u32_t * T, u32_t dT, u32_t * r, u32_t * Q){
   u32_t * P, dP, dQ,i=0,b;
   /* Ersetze T(X) durch T(X+1) */
   if(dT>1){
     polchcomod32(T,dT);
     dQ=polXpotmodPmod32(Q,T,dT);
     Q[0]=modadd32(Q[0],1);
     P=Q+dT+1;
     dP=polcpymod32(P,T,dT);
     dQ=polgcdmod32(P,dP,Q,dQ);
     b=r[0];
     r[0]=modadd32(b,1);
     if(dQ>0 && dQ<dT){
       dQ=polnormmod32(Q,dQ);
       dP=polcpymod32(P,T,dT);
       dP=poldivnormmod32(P,dP,Q,dQ);
       dP=dT-dQ,P+=dQ;
       i=polrootrecmod32(Q,dQ,r,P+dP+1);  /*i=dQ*/
       r[i]=modadd32(b,1);
       i+=polrootrecmod32(P,dP,r+i,P+dP+1); /*i=dT*/
       return i;
     }
     else
       return polrootrecmod32(T,dT,r,Q);
   }
   else {
     r[0]=modsub32(r[0],T[0]);
     return 1;
   }
}

u32_t polvalmod32(u32_t * P, u32_t dP, u32_t a){
   u32_t i=dP,r=P[dP];
   for(;i>0;)
     r=modadd32(P[--i],modmul32(a,r));
   return r;
}

void polprintmod32(u32_t *P, u32_t dP, char * c){
   u32_t i;
   printf("\n%s %u\n",c,dP);
   for(i=0;i<=dP;i++)
     printf("%u ",P[i]);
   printf("\n");
}

void mod32multab(dT){
   u32_t i,c;
   modulo32ebits=0;
   for(i=0,c=modulo32;(c>>=1)>dT;i++){
     if(c&1){
       modulo32bit[i]=1;
       modulo32ebits++;
     }
         else
           modulo32bit[i]=0;
   }
   modulo32bit[i]=c;
   modulo32hbit=i;/* Anzahl der notwendigen Quadratbildungen */
  }

int uintcmp(const void *x,const void *y){
   const u32_t ux =*(const u32_t *) x;
   const u32_t uy =*(const u32_t *) y;
   return ux>=uy?1:-1;
}

u32_t polrootmod32(mpz_t * A, u32_t dT,u32_t * r){
   u32_t i=0,j,*P,dP,*Q,dQ,*S,dS,*T=polr;
   dT=polredmod32(T,A,dT);
   if(T[0]==0){  /* Nullstelle 0 */
     for(i=1;T[i]==0 && i <dT;i++);
     T+=i;dT-=i;
     r[0]=0,i=1;
   }
   if(dT==0)
     return i;
   mod32multab(dT);
   /* Zerlegung Q = (P,X^(modulo32-1)/2+1); S = (P,X^(modulo32-1)/2+1)*/
   P=T+dT+1;
   dP=polcpymod32(P,T,dT);
   Q=P+dP+1;
   dQ=polXpotmodPmod32(Q,P,dP);
   S=Q+dP+1;/* die ggT-Berechnung kann P ergeben */
   dS=polcpymod32(S,Q,dQ);
   Q[0]=modadd32(Q[0],1);
   dQ=polgcdmod32(P,dP,Q,dQ);
   dP=polcpymod32(P,T,dT);
   S[0]=modsub32(S[0],1);
   dS=polgcdmod32(P,dP,S,dS);
   if(dQ>0){
     dQ=polnormmod32(Q,dQ);
     r[i]=0;
     j=polrootrecmod32(Q,dQ,r+i,S+dS+1);  /*j=dQ*/
     if(j!=dQ) fprintf(stderr,"Falsche Nullstellenanzahl Q, P %u\n",modulo32);
     i+=j;
   }
   if(dS>0){
     dS=polnormmod32(S,dS);
     r[i]=0;
     j=polrootrecmod32(S,dS,r+i,S+dS+1);  /*j=dS*/
     if(j!=dS) fprintf(stderr,"Falsche Nullstellenanzahl S, P %u\n",modulo32);
     i+=j;
   }
   qsort(r,i,sizeof(u32_t),uintcmp);
   return i;
}

u32_t
root_finder(u32_t *root_buf,mpz_t *A,u32_t adeg,u32_t p)
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
    if(mpz_fdiv_ui(A[0],p)==0) root_buf[res++]=0;
    for(i=0,pv=0;i<=adeg;i++)
      if(mpz_fdiv_ui(A[i],p)!=0) pv^=1;
    if(pv==0) root_buf[res++]=1;
    if(mpz_fdiv_ui(A[adeg],p)==0)
      root_buf[res++]=2;
    return res;
  }
  modulo32=p;
  res=polrootmod32(A,adeg,root_buf);
  /* CAVE: Diese Division erfolgt wohl ein zweites Mal. */
  if(mpz_fdiv_ui(A[adeg],p)==0) root_buf[res++]=p;
  return res;
}
