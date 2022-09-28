/*
  Copyright (C) 2000 Jens Franke
  This file is part of mpqs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <sys/types.h>
#include <stdio.h>
#include <gmp.h>

#include "siever-config.h"
#include "32bit.h"

volatile u32_t modulo32;

void Schlendrian(char *fmt,...);

u32_t modpow32(u32_t x,u32_t a)
{ /*  returns x^a mod modulo32  */
  u32_t potenz;
  x%=modulo32;
  for(potenz=1;;) {
    if(a&1) potenz=modmul32(x,potenz);
      a>>=1;
    if (a==0) return potenz;
    x=modsq32(x);
  }
}

u32_t modsqrt32(u32_t x)
{
  u32_t q,e,g,i,j,k,l,r,keinePZ;
  x%=modulo32;
  if(x==0) return x;
  if(!(modulo32&1)) Schlendrian("modsqrt: %u not a prime!\n",modulo32);
  /*===================================================\
  | Fuer Zahlen der Form 4n+3 ist das Problem einfach. |
  | Resultat = x^{modulo32+1\over4}                      |
  \===================================================*/
  if(modulo32&2) return modpow32(x,(modulo32+1)>>2);
  /*=======================\
  | Zerlege modulo32-1=q2^e. |
  \=======================*/
  for(e=0,q=modulo32-1;!(q&1);e++,q>>=1);
  /*==========================================================================\
  | Finde Erzeuger g der 2-Torsion der primen Restklassengruppe modulo modulo32 |
  \==========================================================================*/
  for(i=2;i<modulo32;i++) {
    {
      g=modpow32(i,q);
      if(g==1) continue;
      j=g;
      keinePZ=1;
      for(k=1;k<e;k++) {
	if(j==(modulo32-1)) keinePZ=0;
	j=modsq32(j);
	if(j==1 && keinePZ) Schlendrian("modsqrt: %u not a prime!\n",modulo32);
      }
      if(j==modulo32-1) break;
      if(keinePZ) Schlendrian("modsqrt: %u not a prime!\n",modulo32);
    }
  }
  /*====================================================\
  | Setze r=x^{q+1\over2},k=x^q.                        |
  | Bestimme eine Potenz m von g mit km^2=1. r<--rm ist |
  | die Quadratwurzel.                                  |
  \====================================================*/
  r=modpow32(x,(q+1)>>1);
  k=modpow32(x,q);
  while(k!=1)
    {
      j=k;
      for(i=0;i<e && j!=1;i++) j=modsq32(j);
      if(i==e || j!=1)
	Schlendrian("modsqrt: %u not a prime!\n",modulo32);
      for(l=i,l++;l<e;l++) g=modsq32(g);
      r=modmul32(r,g);
      g=modsq32(g);
      k=modmul32(k,g);
      e=i;
    }
  return r;
}

static u32_t x_o;

#if 0
u32_t modinv32(u32_t x)
#else
static u32_t modinv32_eucl(u32_t x)
#endif
{
  /*===========================================================\
  | Anfang: (x,y,z)<--(x,1,0)                                  |
  |         (m,a,b)<--(modulo32,0,1)                             |
  | Iteration: q=m/x,r=m%x;                                    |
  | ((x,y,z),(m,a,b))<--((r;a-qy,b-qz),(x,y,z))                |
  | Ein Tripel (i,j,k) erfuellt stets i=xj+mk.                 |
  | Wir brauchen nur die ersten beiden Komponenten der Tripel, |
  | und auch das nur modulo modulo32.                            |
  \===========================================================*/
  u32_t q,a,r,y,m;
  x_o=x;
  m=modulo32;
  a=0;
  y=1;
  while(x!=0) {
#if 0
    {register u32_t *j=&q;
    asm("divl %2; movl %4,(%3)" : "=edx" (r) :
	"0" (0), "r" (x), "r" (j), "eax" (m));}
#else
    q=m/x;
    r=m%x;
#endif
    m=x;
    x=r;
    r=y;
    y=modsub32(a,modmul32(q,y));
    a=r;
  }
  if(m!=1) Schlendrian("D=%u and M=%u have ggT %u\n",x_o,modulo32,m);
  return a;
}

u32_t modinv32(u32_t x)
{
  /*===========================================================\
  | Anfang: (x,y,z)<--(x,1,0)
  |         (m,a,b)<--(modulo32,0,1) 
  | Iteration: 
            x gerade:
                y gerade: dann auch z gerade
		(x,y,z)<---(x/2,y/2,z/2)
                y ungerade:
                (x,y,z)<---(x/2,(y+modulo32)/2,(z-x_o)/2) oder
                (x,y,z)<---(x/2,(y-modulo32)/2,(z+x_o)/2)
            x ungerade:
	        x>m:
		((x,y,z),(m,a,b))<---((x-m,y-a,z-b),(m,a,b))
		x<m:
		((x,y,z),(m,a,b))<---((m-x,a-y,z-b),(x,y,z))
		x=m: Beide 1: return a; else Krach.
 Ein Tripel (i,j,k) erfuellt stets i=xj+mk. 
 Wir brauchen nur die ersten beiden Komponenten der Tripel, 
 und auch das nur modulo modulo32. */
  u32_t m,y,a,xx,yy;

  if(modulo32%2==0) return modinv32_eucl(x);
  x_o=x;
  m=modulo32;
  a=0;
  y=1;
  if(x==0) Schlendrian("Modular inversion of 0\n");
  for(;;) { /* m ist stets ungerade */
    while(!(x&1)) {
      x>>=1;
      if(y&1) {
	if(y>=modulo32) y=(y-modulo32)>>2;
	else y=y/2+modulo32/2+1;
      } else y>>=1;
    }
    if(x>m) {
      x-=m;
      if(y<a) y=y+(modulo32-a);
      else y-=a;
    } else if(x<m) {
      xx=x;
      yy=y;
      x=m-x;
      if(a<y) y=a+(modulo32-y);
      else y=a-y;
      m=xx;
      a=yy;
    } else { /* x=m */
      if(x==1) return a;
      else Schlendrian("D=%u and M=%u have ggT %u\n",x_o,modulo32,m);
    }
  }
}

