@* Scheduling function for the lattice siever.
Copyright (C) 2002 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
@(lasched.h@>=
u32_t *lasched(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t,u32_t);
u32_t *lasched_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t,u32_t);

@
@c
#include <sys/types.h>
#include <limits.h>

#include "siever-config.h"
#include "lasched.h"
#include "../if.h"

#define L1_SIZE (1<<L1_BITS)
#define i_bits (I_bits-1)
#define n_i (1<<i_bits)

/* Amount by which to shift a larger number to provide space for a
   16-bit number. */

#define U16_SHIFT (CHAR_BIT*sizeof(u16_t))

u32_t * ASM_ATTR lasched0(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched2(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched3(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched0nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched1nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched2nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched3nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);

u32_t * ASM_ATTR lasched0_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched1_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched2_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched3_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);

u32_t *
lasched_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs,ot)
     u32_t *ri,*ij_ptr,*ij_ptr_ub,n1_j,**sched_ptr,fbi_offs,ot;
{
  u32_t ij,ij_ub;
  u32_t ot_mask,ot_tester;

  ij_ub=n1_j<<i_bits;

  switch(ot) {
  case 0:
    return lasched0_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  case 1:
    return lasched1_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  case 2:
    return lasched2_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  default:
    return lasched3_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  }
  return NULL;
}






u32_t *
lasched(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs,ot)
     u32_t *ri,*ij_ptr,*ij_ptr_ub,n1_j,**sched_ptr,fbi_offs,ot;
{
  u32_t ij,ij_ub;
  u32_t ot_mask,ot_tester;

  ij_ub=n1_j<<i_bits;

#ifdef SCHED_NT_BOUND
  if(ij_ub>SCHED_NT_BOUND)
    switch(ot) {
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

  switch(ot) {
  case 0:
    return lasched0(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  case 1:
    return lasched1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  case 2:
    return lasched2(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  default:
    return lasched3(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  }

  if(ot!=0) {
    ot_tester=(ot&1)|((ot&2)<<(i_bits-1));
    ot_mask=n_i|1;
  } 

  while ( ij_ptr < ij_ptr_ub) {
    u16_t a,b;

    a=n_i-(ri[0]&(n_i-1));
    b=n_i-(ri[1]&(n_i-1));
    if(ot == 0) ij=*ij_ptr;
    else @<Calculate first sieving event from |ri|@>@;
    while(ij<ij_ub) {
      u16_t i;

      *(sched_ptr[ij>>L1_BITS]++)=(fbi_offs<<U16_SHIFT) | (ij&(L1_SIZE-1));
      i=ij&(n_i-1);
      if(i<b) ij+=ri[0];
      if(i>=a) ij+=ri[1];
    }
    ri+=2;
    *(ij_ptr++)=ij-ij_ub;
    fbi_offs++;
  }
  return ri;
}

@
@<Calculate first sieving event from |ri|@>=
{
  ij=0;
  if( (ri[0]&ot_mask) == ot_tester ) ij=ri[0];
  else {
    if( (ri[1]&ot_mask) == (ot_tester^n_i) ) ij=ri[1];
    else {
      if((ri[0]&(n_i-1))<=(ri[1]&(n_i-1)) && ri[0]<=ri[1]) {
	/* This corresponds to the line
	   |if(b+c<=A && s<=t)| in recurrence6.w */
	if((ri[0]&(n_i-1))==(ri[1]&(n_i-1))) ij=ri[1]-ri[0];
	else ij=n_i;
	if(ot != 2)
	  Schlendrian("Exceptional situation for oddness type %u ?\n",
		      ot);
      }
      else ij=ri[0]+ri[1];
    }
  }
  ij=(ij+((~ot_tester)&n_i))/2;
}
