/* recurrence6.c
  By Jens Franke.
  6/13/04: Hacked up for use in GGNFS by Chris Monico.
                                                                                                                                                                                                             
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
#include <limits.h>

#include "asm/siever-config.h"
#include "if.h"
#include "recurrence6.h"

static u32_t A,A_bits,ub;

void rec_info_init(u32_t A1,u32_t ub1)
{
  u32_t aa;

  if(A1> USHRT_MAX)complain("Recurrence init: A=%u exceeds %u\n",A1,USHRT_MAX);
  A= A1;
  if(A%2!=0)complain("rec_info_init with odd range %u\n",A1);
  if(ub1> (USHRT_MAX / 4 + 1))
    complain("Recurrence init: ub=%u exceeds %u\n",ub1,(USHRT_MAX / 4 + 1));
  if(ub1<2)
    complain("Recurrence time %u too small\n",ub1);
  ub= ub1;
  for(aa= 1,A_bits= 0;A_bits<=CHAR_BIT*sizeof(u32_t);A_bits++,aa*= 2)
    if(aa==A)break;
  if(A_bits> CHAR_BIT*sizeof(u32_t))
    complain("rec_info_init: A=%u not a power of two\n",A);
}

u32_t get_recurrence_info(u32_t*res_ptr,u32_t p,u32_t r)
{
  u32_t b,c,s,t;

  if(r == 0) {
    b= 1;
    s= 2*ub+1;
    c= A - 3;
    t= 2*ub;

    goto done;
  }
  if(r >= p) {

    b= 1;
    s= 0;
    c= A-1;
    t= 2*ub+1;

    goto done;
}

#ifdef HAVE_ASM_GETBC
asm_getbc(r,p,A,&b,&s,&c,&t);
#else
{
b= r;
s= 1;
c= p;
t= 0;
for(;;){
u32_t k;

if(b<A){

k= (c-A)/b+1;
t+= k*s;
c-= k*b;
break;
}
k= c/b;
c= c%b;
t+= k*s;
if(c<A){

k= (b-A)/c+1;
s+= k*t;
b-= k*c;
break;
}
k= b/c;
b= b%c;
s+= k*t;
}
}
#endif

#if 0
{
  u32_t b1,c1,s1,t1;

  asm_getbc(r,p,A,&b1,&s1,&c1,&t1);
  if(b1!=b||c1!=c||s1!=s||t1!=t)Schlendrian("Asm ri\n");
}
#endif

  {
    i32_t d;

    d= t-s;
    s= (s<=2*ub?s:2*ub+2-(s&1));
    if(t> 2*ub){
      if(d<0||d> 2*ub)d= 2*ub+2+(s&1)-(t&1);
      t= s+d;
    }
  }

done:
  {
    res_ptr[0]= (s<<A_bits)+b;
    res_ptr[1]= (t<<A_bits)-c;
  }

#if 0
{
u32_t z[2];
u32_t have_it[3]= {0,0,0};

z[0]= A+b;
z[1]= s;
{
u32_t ot;
#ifdef ULONG_RI
u32_t*xx;
#else
u16_t*xx;
#endif

ot= (z[0]&1)+2*(z[1]&1);
if(ot==1)xx= x1;
else if(ot==2)xx= x2;
else xx= x3;
have_it[ot-1]= 1;
#ifdef ULONG_RI
*xx= ((z[1]/2)<<A_bits)|(z[0]/2);
#else
xx[0]= z[0]/2;
xx[1]= z[1]/2;
#endif
}

if(b+c<=A&&s<=t){
if(b+c<A){
z[0]= A;
z[1]= 1;
}else{
z[0]= 0;
z[1]= t-s;
}
}else{
z[0]= A+b-c;
z[1]= t+s;
}
#if 0
z[1]= (z[1]<=2*ub?z[1]:2*ub+2-(z[1]&1));
#else
z[1]= (z[1]<=USHRT_MAX?z[1]:2*ub+2-(z[1]&1));
#endif
{
u32_t ot;
#ifdef ULONG_RI
u32_t*xx;
#else
u16_t*xx;
#endif

ot= (z[0]&1)+2*(z[1]&1);
if(ot==1)xx= x1;
else if(ot==2)xx= x2;
else xx= x3;
have_it[ot-1]= 1;
#ifdef ULONG_RI
*xx= ((z[1]/2)<<A_bits)|(z[0]/2);
#else
xx[0]= z[0]/2;
xx[1]= z[1]/2;
#endif

z[0]= A-c;
z[1]= t;

{
u32_t ot;
#ifdef ULONG_RI
u32_t*xx;
#else
u16_t*xx;
#endif

ot= (z[0]&1)+2*(z[1]&1);
if(ot==1)xx= x1;
else if(ot==2)xx= x2;
else xx= x3;
have_it[ot-1]= 1;
#ifdef ULONG_RI
*xx= ((z[1]/2)<<A_bits)|(z[0]/2);
#else
xx[0]= z[0]/2;
xx[1]= z[1]/2;
#endif

if(have_it[0]==0)Schlendrian("???");
}

#endif
return 2;
}
