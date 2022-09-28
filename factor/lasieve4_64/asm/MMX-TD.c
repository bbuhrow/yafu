#include <sys/types.h>
#include <malloc.h>
#include <limits.h>

#include "siever-config.h"
#include "32bit.h"
#include "../if.h"

/*
 * Auxilliary information for trial division using MMX instructions,
 * stored in the form (root0,root1,root2,root3,prime0,
 *                     prime0,prime1,prime2,prime3,
 *                     mi0,mi1,mi2,mi3),
 * where the mi are modular Inverses to the primes.
 */

static u16_t *(MMX_TdAux[2]),*(MMX_TdBound[2]);

/*
 * The projective roots used to update MMX_TdAux when the line is changed.
 */

static u16_t **(MMX_TdPr[2]);

static size_t MMX_TdAlloc[2]={0,0};

static int jps=0;

#ifndef MMX_REGW
#define MMX_REGW 4
#endif

/* Read-Ahead safety. */
#define RAS MMX_REGW

static u16_t *
mmx_xmalloc(size_t n)
{
  u16_t *r;
  r=(u16_t*)memalign(MMX_REGW*sizeof(u16_t),n*sizeof(u16_t));
  if(r==NULL)
    complain("mmx_malloc(%Lu) failed: %m\n",(u64_t)n);
  return r;
}

void
MMX_TdAllocate(int jps_arg,size_t s0,size_t s1)
{
  int side;

  MMX_TdAlloc[0]=MMX_REGW*((s0+MMX_REGW-1)/MMX_REGW);
  MMX_TdAlloc[1]=MMX_REGW*((s1+MMX_REGW-1)/MMX_REGW);
  jps=jps_arg;

  for(side=0;side<2;side++) {
    int i;

    MMX_TdAux[side]=mmx_xmalloc(3*(MMX_TdAlloc[side]+RAS));
    MMX_TdPr[side]=xmalloc(jps*sizeof(*(MMX_TdPr[side])));
    for(i=0;i<jps;i++)
      MMX_TdPr[side][i]=mmx_xmalloc(MMX_TdAlloc[side]+RAS);
  }
}

u16_t*
MMX_TdInit(int side,u16_t *x,u16_t *x_ub,u32_t *pbound_ptr,
	   int initialize)
{
  u16_t *y,*z,*u;
  u32_t p_bound;

  if(initialize==1) {
    int i;
    if(x_ub>x+4*MMX_TdAlloc[side])
      Schlendrian("Buffer overflow in MMX_TdInit\n");
    z=MMX_TdAux[side]+MMX_REGW;
    y=x;
    i=0;
    while(y+4*MMX_REGW<x_ub) {
      int j;
      for(j=0;j<MMX_REGW;j++,y+=4,z++,i++) {
	u16_t mi;
	u32_t r,rr;
	int k;

	modulo32=*y;
	mi=*y;
	*z=mi;
	mi=2*mi-mi*mi*modulo32;
	mi=2*mi-mi*mi*modulo32;
	mi=2*mi-mi*mi*modulo32;
	*(z+MMX_REGW)=mi;
	r=y[1];
	rr=r;
	for(k=0;k<jps;k++) {
	  MMX_TdPr[side][k][i]=rr;
	  rr=modadd32(rr,r);
	}
      }
      z+=2*MMX_REGW;
    }
  }
  z=MMX_TdAux[side];
  u=MMX_TdPr[side][jps-1];
  p_bound=*pbound_ptr;
  for(y=x;y<x_ub-4*MMX_REGW;y=y+4*MMX_REGW,z+=3*MMX_REGW,u+=MMX_REGW) {
    int i;
    if(z[MMX_REGW]>p_bound) break;
    for(i=0;i<MMX_REGW;i++) {
      modulo32=y[4*i];
      z[i]=modsub32(0,y[4*i+3]);
      y[4*i+3]=modadd32(y[4*i+3],u[i]);
    }
  }
  *pbound_ptr=z[-MMX_REGW-1];
  MMX_TdBound[side]=z;
  return y;
}

void
MMX_TdUpdate(int side,int j_step)
{
#if 0
  u16_t *x,*y;

  y=MMX_TdPr[side][j_step-1];
  for(x=MMX_TdAux[side];x<MMX_TdBound[side];x+=3*MMX_REGW,y+=MMX_REGW) {
    int i;
    for(i=0;i<MMX_REGW;i++) {
      modulo32=x[MMX_REGW+i];
      x[i]=modsub32(x[i],y[i]);
    }
  }
#else
#if MMX_REGW == 4
  asm_TdUpdate4(MMX_TdAux[side],MMX_TdBound[side],MMX_TdPr[side][j_step-1]);
#elif MMX_REGW == 8
  asm_TdUpdate8(MMX_TdAux[side],MMX_TdBound[side],MMX_TdPr[side][j_step-1]);
#else
#error "No asm for this MMX_REGW"
#endif
#endif
}

u32_t *asm_MMX_Td4(u32_t*,u32_t,u16_t*,u16_t*);
u32_t *asm_MMX_Td8(u32_t*,u32_t,u16_t*,u16_t*);

#ifdef MMX_TDBENCH
u64_t MMX_TdNloop=0;
#endif

u32_t *
MMX_Td(u32_t *pbuf,int side,u16_t strip_i)
{
#if 1
#ifdef MMX_TDBENCH
  MMX_TdNloop+=(MMX_TdBound[side]-MMX_TdAux[side])/MMX_REGW;
#endif
#if MMX_REGW == 4
  return asm_MMX_Td4(pbuf,strip_i,MMX_TdAux[side],MMX_TdBound[side]);
#elif MMX_REGW == 8
  return asm_MMX_Td8(pbuf,strip_i,MMX_TdAux[side],MMX_TdBound[side]);
#else
#error "No asm for this MMX_REGW"
#endif
#else
  u16_t* x;

  for(x=MMX_TdAux[side];x<MMX_TdBound[side];x+=2*MMX_REGW) {
    int i;

    for(i=0;i<MMX_REGW;i++,x++) {
      u16_t t;

      modulo32=x[MMX_REGW];
    
      t=strip_i+x[0];
      t*=x[2*MMX_REGW];
      if(((modulo32*(u32_t)t)&0xffff0000)==0)
	*(pbuf++)=modulo32;
    }
  }
  return pbuf;
#endif
}
