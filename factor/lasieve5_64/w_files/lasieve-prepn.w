@* Preparation of lattice sieving.

@*3 Copying.
Copyright (C) 2001 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@ These functions should be somewhat quicker than the straightforward thing
since unnecessary divisions and jumps in the inner loop are eliminated.
@(lasieve-prepn.h@>=
void lasieve_setup(u32_t*,u32_t*,u32_t,i32_t,i32_t,i32_t,i32_t,u32_t*,i32_t);
void lasieve_setup64(u32_t*,u32_t*,u32_t,i64_t,i64_t,i64_t,i64_t,u32_t*,i32_t);

@
@c
#include <sys/types.h>
#include <math.h>

#include "if.h"
#include "asm/siever-config.h"
#include "recurrence6.h"
#include "asm/32bit.h"

#include <stdio.h>
#include <stdlib.h>

void
lasieve_setup(FB,proots,fbsz,a0, a1, b0, b1,ri_ptr,d)
  u32_t *FB,*proots,fbsz,*ri_ptr;
i32_t a0,a1,b0,b1,d;
{
  u32_t b0_ul,b1_ul,absa0,absa1;

  if(b0<0 || b1<0)
    Schlendrian("lasieve_setup called with negative 2-nd coordinate (%d,%d)\n",
               b0,b1);
  if(fbsz<=0) return;
#ifdef HAVE_ASM_LASIEVE_SETUP
  @<Use floating point arithmetic if possible@>@;
#endif
  b0_ul=(u32_t)b0;
  b1_ul=(u32_t)b1;
#ifdef HAVE_MM_LASIEVE_SETUP
  @<Use Montgomery multiplication arithmetic if possible@>@;
#endif
  if(a0>=0) {
    absa0=(u32_t)a0;
#define A0MOD0(p) (absa0%p)
    /* The next macro applies when a0 is smaller than p. */
#define A0MOD1(p) absa0
    if(a1>=0) {
      absa1=(u32_t)a1;
      /* Both a-coordinates non-negative. */
      /* |#define| how to obtain the remainder of a0 modulo a
         factor base element. */
#define A1MOD0(p) (absa1%p)
#define A1MOD1(p) absa1
      @<Do the setup32@>@;
#undef A1MOD0
#undef A1MOD1
    } else {
      u32_t aux;
      /* |a0>=0 && a1<0 */
      absa1=(u32_t)(-a1);
#define A1MOD0(p) ((aux=absa1%p)>0 ? p-aux : 0 )
#define A1MOD1(p) (p-absa1)
      @<Do the setup32@>@;
#undef A1MOD0
#undef A1MOD1
#undef A0MOD0
#undef A0MOD1
    }
  } else {
    absa0=(u32_t)(-a0);
#define A0MOD0(p) ((aux=absa0%p)>0 ? p-aux : 0 )
#define A0MOD1(p) (p-absa0)
    if(a1>=0) {
      u32_t aux;
      absa1=(u32_t)a1;
#define A1MOD0(p) (absa1%p)
#define A1MOD1(p) absa1
      @<Do the setup32@>@;
#undef A1MOD0
#undef A1MOD1
    } else {
      u32_t aux;
      absa1=(u32_t)(-a1);
#define A1MOD0(p) ((aux=absa1%p)>0 ? p-aux : 0 )
#define A1MOD1(p) (p-absa1)
      @<Do the setup32@>@;
#undef A1MOD0
#undef A1MOD1
    }
  }
}

@
@<Do the setup32@>=
{
  u32_t fbi,fbp_bound;

  /* Do the stuff for the small primes for which taking residues mod p
     of the $a_i$ or $b_i$ still is a non-trivial operation. */
#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
  fbp_bound= absa1<absa0 ? absa0 : absa1;
  if(fbp_bound<b0_ul) fbp_bound=b0_ul;
  if(fbp_bound<b1_ul) fbp_bound=b1_ul;
  for(fbi=0;fbi<fbsz && FB[fbi]<=fbp_bound;fbi++) @<Setup loop32@>@;
#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)
  /* Now, we are in a situation where the absolute values of all four
     numbers |a0|, |a1|, |b0| and |b1| are smaller than the factor base
     element we are considering. We |#define| all macros to do the
     simplest possible things. */
  for(;fbi<fbsz;fbi++) @<Setup loop32@>@;
#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
}

@
@<Setup loop32@>=
{
  u32_t x;
  modulo32=FB[fbi];

  if(proots[fbi]==modulo32) {
    x=B0MOD(modulo32);
    if(x==0) x=modulo32;
    else {
      x=modmul32(modinv32(x),B1MOD(modulo32));
      if(x>0) x=modulo32-x;
    }
  } else {
    x=modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
    if(x!=0) {
      x=modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
                                          A1MOD(modulo32)));
    } else {
      /* Case $\infty$. */
      x=FB[fbi];
    }
  }
  ri_ptr+=get_recurrence_info(ri_ptr,FB[fbi],x);
}

@
@<Use floating point arithmetic if possible@>=
if(FB[fbsz-1]<FLOAT_SETUP_BOUND1 &&
   fabs(a0)+FB[fbsz-1]*b0<FLOAT_SETUP_BOUND2 &&
   fabs(a1)+FB[fbsz-1]*b1<FLOAT_SETUP_BOUND2) {
     asm_lasieve_setup(FB,proots,fbsz,a0,a1,b0,b1,ri_ptr);
     return;
   }

@
@<Use Montgomery multiplication arithmetic if possible@>=
{
  if(d==1) {
    if(a0>=0) {
      absa0=(u32_t)a0;
      if(a1>=0) {
	absa1=(u32_t)a1;
	/* Both a-coordinates non-negative. */
	ri_ptr=asm_lasieve_mm_setup0(FB,proots,fbsz,
				     absa0,absa1,b0_ul,b1_ul,ri_ptr);
      } else {
	u32_t aux;
	/* |a0>=0 && a1<0 */
	absa1=(u32_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup1(FB,proots,fbsz,
				     absa0,absa1,b0_ul,b1_ul,ri_ptr);
      }
    } else {
      absa0=(u32_t)(-a0);
      if(a1>=0) {
	u32_t aux;
	absa1=(u32_t)a1;
	ri_ptr=asm_lasieve_mm_setup2(FB,proots,fbsz,
				     absa0,absa1,b0_ul,b1_ul,ri_ptr);
      } else {
	u32_t aux;
	absa1=(u32_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup3(FB,proots,fbsz,
				     absa0,absa1,b0_ul,b1_ul,ri_ptr);
      }
    }
    return;
  }


#ifdef HAVE_MM_LASIEVE_SETUP2
  if(d==2) {
    if(a0>=0) {
      absa0=(u32_t)a0;
      if(a1>=0) {
	absa1=(u32_t)a1;
	/* Both a-coordinates non-negative. */
	ri_ptr=asm_lasieve_mm_setup20(FB,proots,fbsz,
				      absa0,absa1,b0_ul,b1_ul,ri_ptr);
      } else {
	u32_t aux;
	/* |a0>=0 && a1<0 */
	absa1=(u32_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup21(FB,proots,fbsz,
				      absa0,absa1,b0_ul,b1_ul,ri_ptr);
      }
    } else {
      absa0=(u32_t)(-a0);
      if(a1>=0) {
	u32_t aux;
	absa1=(u32_t)a1;
	ri_ptr=asm_lasieve_mm_setup22(FB,proots,fbsz,
				      absa0,absa1,b0_ul,b1_ul,ri_ptr);
      } else {
	u32_t aux;
	absa1=(u32_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup23(FB,proots,fbsz,
				      absa0,absa1,b0_ul,b1_ul,ri_ptr);
      }
    }
    return;
  }
#endif
  Schlendrian("No special setup function for d=%d\n",(int)d);
}


@
@c
#ifdef VERY_LARGE_Q
void
lasieve_setup64(u32_t *FB, u32_t *proots, u32_t fbsz, i64_t a0, i64_t a1,
                i64_t b0, i64_t b1, u32_t *ri_ptr, i32_t d)
{
  u64_t b0_ull,b1_ull,absa0,absa1;

  if(b0<0 || b1<0)
    Schlendrian("lasieve_setup called with negative 2-nd coordinate (%lld,%lld)\n",
               b0,b1);
  if(fbsz<=0) return;
  b0_ull=(u64_t)b0;
  b1_ull=(u64_t)b1;
#ifdef HAVE_MM_LASIEVE_SETUP_64
  @<Use Montgomery multiplication arithmetic for 64 bit if possible@>@;
#endif
  if(a0>=0) {
    absa0=(u64_t)a0;
    if(a1>=0) {
      absa1=(u64_t)a1;
      /* Both a-coordinates non-negative. */
      @<Do the setup64pp@>@;
    } else {
      u64_t aux;
      /* |a0>=0 && a1<0 */
      absa1=(u64_t)(-a1);
      @<Do the setup64pn@>@;
    }
  } else {
    absa0=(u64_t)(-a0);
    if(a1>=0) {
      absa1=(u64_t)a1;
      @<Do the setup64np@>@;
    } else {
      absa1=(u64_t)(-a1);
      @<Do the setup64nn@>@;
    }
  }
}
#endif

@
@<Do the setup64pp@>=
{
  u64_t fbi;

  for(fbi=0;fbi<fbsz;fbi++) {
    u32_t x;
    u64_t p64, h;

    modulo32=FB[fbi];
    p64=(u64_t)modulo32;

    if(proots[fbi]==modulo32) {
      x=((u32_t)(b0_ull%p64));
      if(x==0) x=modulo32;
      else {
        x=modmul32(modinv32(x),((u32_t)(b1_ull%p64)));
        if(x>0) x=modulo32-x;
      }
    } else {
      h=absa0+((u64_t)(modulo32-proots[fbi]))*(b0_ull%p64);
      x=(u32_t)(h%p64);
      if(x!=0) {
        h=absa1+((u64_t)(modulo32-proots[fbi]))*(b1_ull%p64);
        x=modulo32-x;
        x=modmul32(asm_modinv32(x),(u32_t)(h%p64));
      } else {
        /* Case $\infty$. */
        x=FB[fbi];
      }
    }
    ri_ptr+=get_recurrence_info(ri_ptr,FB[fbi],x);
  }
}

@
@<Do the setup64pn@>=
{
  u64_t fbi;

  for(fbi=0;fbi<fbsz;fbi++) {
    u32_t x;
    u64_t p64, h;

    modulo32=FB[fbi];
    p64=(u64_t)modulo32;

    if(proots[fbi]==modulo32) {
      x=((u32_t)(b0_ull%p64));
      if(x==0) x=modulo32;
      else {
        x=modmul32(modinv32(x),((u32_t)(b1_ull%p64)));
        if(x>0) x=modulo32-x;
      }
    } else {
      h=absa0+((u64_t)(modulo32-proots[fbi]))*(b0_ull%p64);
      x=(u32_t)(h%p64);
      if(x!=0) {
        h=absa1+((u64_t)proots[fbi])*(b1_ull%p64);
        x=modmul32(asm_modinv32(x),(u32_t)(h%p64));
      } else {
        /* Case $\infty$. */
        x=FB[fbi];
      }
    }
    ri_ptr+=get_recurrence_info(ri_ptr,FB[fbi],x);
  }
}

@
@<Do the setup64np@>=
{
  u64_t fbi;

  for(fbi=0;fbi<fbsz;fbi++) {
    u32_t x;
    u64_t p64, h;

    modulo32=FB[fbi];
    p64=(u64_t)modulo32;

    if(proots[fbi]==modulo32) {
      x=((u32_t)(b0_ull%p64));
      if(x==0) x=modulo32;
      else {
        x=modmul32(modinv32(x),((u32_t)(b1_ull%p64)));
        if(x>0) x=modulo32-x;
      }
    } else {
      h=absa0+((u64_t)proots[fbi])*(b0_ull%p64);
      x=(u32_t)(h%p64);
      if(x!=0) {
        h=absa1+((u64_t)(modulo32-proots[fbi]))*(b1_ull%p64);
        x=modmul32(asm_modinv32(x),(u32_t)(h%p64));
      } else {
        /* Case $\infty$. */
        x=FB[fbi];
      }
    }
    ri_ptr+=get_recurrence_info(ri_ptr,FB[fbi],x);
  }
}

@
@<Do the setup64nn@>=
{
  u64_t fbi;

  for(fbi=0;fbi<fbsz;fbi++) {
    u32_t x;
    u64_t p64, h;

    modulo32=FB[fbi];
    p64=(u64_t)modulo32;

    if(proots[fbi]==modulo32) {
      x=((u32_t)(b0_ull%p64));
      if(x==0) x=modulo32;
      else {
        x=modmul32(modinv32(x),((u32_t)(b1_ull%p64)));
        if(x>0) x=modulo32-x;
      }
    } else {
      h=absa0+((u64_t)proots[fbi])*(b0_ull%p64);
      x=(u32_t)(h%p64);
      if(x!=0) {
        x=modulo32-x;
        h=absa1+((u64_t)proots[fbi])*(b1_ull%p64);
        x=modmul32(asm_modinv32(x),(u32_t)(h%p64));
      } else {
        /* Case $\infty$. */
        x=FB[fbi];
      }
    }
    ri_ptr+=get_recurrence_info(ri_ptr,FB[fbi],x);
  }
}

@
@<Use Montgomery multiplication arithmetic for 64 bit if possible@>=
{
  if(d==1) {
    if(a0>=0) {
      absa0=(u64_t)a0;
      if(a1>=0) {
	absa1=(u64_t)a1;
	/* Both a-coordinates non-negative. */
	ri_ptr=asm_lasieve_mm_setup0_64(FB,proots,fbsz,
				     absa0,absa1,b0_ull,b1_ull,ri_ptr);
      } else {
	/* |a0>=0 && a1<0 */
	absa1=(u64_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup1_64(FB,proots,fbsz,
				     absa0,absa1,b0_ull,b1_ull,ri_ptr);
      }
    } else {
      absa0=(u64_t)(-a0);
      if(a1>=0) {
	absa1=(u64_t)a1;
	ri_ptr=asm_lasieve_mm_setup2_64(FB,proots,fbsz,
				     absa0,absa1,b0_ull,b1_ull,ri_ptr);
      } else {
	absa1=(u64_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup3_64(FB,proots,fbsz,
				     absa0,absa1,b0_ull,b1_ull,ri_ptr);
      }
    }
    return;
  }

#ifdef HAVE_MM_LASIEVE_SETUP2
  if(d==2) {
    if(a0>=0) {
      absa0=(u64_t)a0;
      if(a1>=0) {
	absa1=(u64_t)a1;
	/* Both a-coordinates non-negative. */
	ri_ptr=asm_lasieve_mm_setup20_64(FB,proots,fbsz,
				      absa0,absa1,b0_ull,b1_ull,ri_ptr);
      } else {
	/* |a0>=0 && a1<0 */
	absa1=(u64_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup21_64(FB,proots,fbsz,
				      absa0,absa1,b0_ull,b1_ull,ri_ptr);
      }
    } else {
      absa0=(u64_t)(-a0);
      if(a1>=0) {
	absa1=(u64_t)a1;
	ri_ptr=asm_lasieve_mm_setup22_64(FB,proots,fbsz,
				      absa0,absa1,b0_ull,b1_ull,ri_ptr);
      } else {
	absa1=(u64_t)(-a1);
	ri_ptr=asm_lasieve_mm_setup23_64(FB,proots,fbsz,
				      absa0,absa1,b0_ull,b1_ull,ri_ptr);
      }
    }
    return;
  }
#endif
  Schlendrian("No special setup function for d=%d\n",(int)d);
}
