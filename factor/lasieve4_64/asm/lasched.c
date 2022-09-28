/*2:*/
#line 14 "lasched.w"

#include <sys/types.h> 
#include <limits.h> 

#include "siever-config.h"
#include "lasched.h"
#include "../if.h"

#define L1_SIZE (1<<L1_BITS)
#define i_bits (I_bits-1)
#define n_i (1<<i_bits)




#define U16_SHIFT (CHAR_BIT*sizeof(u16_t))

u32_t*lasched0(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t*lasched1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t*lasched2(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t*lasched3(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t*lasched0nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t*lasched1nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t*lasched2nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t*lasched3nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);

u32_t*
lasched(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs,ot)
u32_t*ri,*ij_ptr,*ij_ptr_ub,n1_j,**sched_ptr,fbi_offs,ot;
{
u32_t ij,ij_ub;
u32_t ot_mask,ot_tester;

ij_ub= n1_j<<i_bits;

#ifdef SCHED_NT_BOUND
if(ij_ub> SCHED_NT_BOUND)
switch(ot){
case 0:
return lasched0nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
case 1:
return lasched1nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
case 2:
return lasched2nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
default:
return lasched3nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
}
#endif

switch(ot){
case 0:
return lasched0(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
case 1:
return lasched1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
case 2:
return lasched2(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
default:
return lasched3(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
}

if(ot!=0){
ot_tester= (ot&1)|((ot&2)<<(i_bits-1));
ot_mask= n_i|1;
}

while(ij_ptr<ij_ptr_ub){
u16_t a,b;

a= n_i-(ri[0]&(n_i-1));
b= n_i-(ri[1]&(n_i-1));
if(ot==0)ij= *ij_ptr;
else/*3:*/
#line 102 "lasched.w"

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
#line 85 "lasched.w"

while(ij<ij_ub){
u16_t i;

*(sched_ptr[ij>>L1_BITS]++)= (fbi_offs<<U16_SHIFT)|(ij&(L1_SIZE-1));
i= ij&(n_i-1);
if(i<b)ij+= ri[0];
if(i>=a)ij+= ri[1];
}
ri+= 2;
*(ij_ptr++)= ij-ij_ub;
fbi_offs++;
}
return ri;
}

/*:2*/
