/*2:*/
#line 14 "medsched.w"

#include <sys/types.h> 
#include <limits.h> 

#include "siever-config.h"
#include "medsched.h"
#include "../if.h"

#define L1_SIZE (1<<L1_BITS)
#define i_bits (I_bits-1)
#define n_i (1<<i_bits)

#define U16_SHIFT (CHAR_BIT*sizeof(u16_t))

u32_t*medsched0(u32_t*,u32_t*,u32_t*,u32_t**,u32_t);

u32_t*
medsched(ri,ij_ptr,ij_ptr_ub,sched_ptr,fbi_offs,ot)
u32_t*ri,*ij_ptr,*ij_ptr_ub,**sched_ptr,fbi_offs,ot;
{
u32_t ij;
u32_t ot_mask,ot_tester;
u32_t*sched;

sched= *sched_ptr;

if(ot!=0){
ot_tester= (ot&1)|((ot&2)<<(i_bits-1));
ot_mask= n_i|1;
}else return medsched0(ri,ij_ptr,ij_ptr_ub,sched_ptr,fbi_offs);

while(ij_ptr<ij_ptr_ub){
u16_t a,b;

a= n_i-(ri[0]&(n_i-1));
b= n_i-(ri[1]&(n_i-1));
if(ot==0)ij= *ij_ptr;
else/*3:*/
#line 69 "medsched.w"

{
ij= 0;
if((ri[0]&ot_mask)==ot_tester)ij= ri[0];
else{
if((ri[1]&ot_mask)==(ot_tester^n_i))ij= ri[1];
else{
if((ri[0]&(n_i-1))<=(ri[1]&(n_i-1))&&ri[0]<=ri[1]){


if((ri[0]&(n_i-1))==(ri[1]&(n_i-1)))ij= ri[1]-ri[0];
else ij= n_i;
if(ot!=2)
Schlendrian("Exceptional situation for oddness type %u ?\n",
ot);
}
else ij= ri[0]+ri[1];
}
}
ij= (ij+((~ot_tester)&n_i))/2;
}/*:3*/
#line 51 "medsched.w"

while(ij<L1_SIZE){
u16_t i;

*(sched++)= (fbi_offs<<U16_SHIFT)|ij;
i= ij&(n_i-1);
if(i<b)ij+= ri[0];
if(i>=a)ij+= ri[1];
}
ri+= 2;
*(ij_ptr++)= ij-L1_SIZE;
fbi_offs++;
}
*sched_ptr= sched;
return ri;
}

/*:2*/
