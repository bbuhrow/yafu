#include <stdio.h>
#include <limits.h>
#include <sys/types.h>
#include <math.h>

#include "asm/siever-config.h"
#include "if.h"
#include "real-poly-aux.h"

static double *aux;
static size_t aux_alloc=0;

/* Action of GL(2,Z) on homogeneous polys. */

void
tpol(double *rop,double *op,i32_t deg,i32_t x0,i32_t x1,i32_t y0,i32_t y1)
{
  i32_t d;

  if(deg==U32_MAX)
    complain("Degree too large\n");
  if(deg==0) {
    rop[0]=op[0];
    return;
  }
  if(aux_alloc<deg+1) {
    if(aux_alloc>0) free(aux);
    aux_alloc=deg+1;
    aux=xmalloc(aux_alloc*sizeof(*aux));
  }
  rop[0]=op[0]*y1;
  rop[1]=op[0]*y0;
  aux[0]=x1;
  aux[1]=x0;
  for(d=1;;d++) {
    i32_t i;
    for(i=0;i<=d;i++) rop[i]+=op[d]*aux[i];
    if(d==deg) return;
    /* Multiply rop by (y0*T+y1), where T is the free variable. */
    rop[d+1]=y0*rop[d];
    for(i=d;i>0;i--)
      rop[i]=y1*rop[i]+y0*rop[i-1];
    rop[0]*=y1;
    /* Multiply aux by (x0*T+x1). */
    aux[d+1]=x0*aux[d];
    for(i=d;i>0;i--)
      aux[i]=x1*aux[i]+x0*aux[i-1];
    aux[0]*=x1;
  }
}

double
rpol_eval(double *p,i32_t d,double x,double y)
{
  double v,z;
  i32_t i;
  for(i=0,v=p[0],z=x;i<d;) {
    v=y*v+p[++i]*z;
    z*=x;
  }
  return v;
}

/* Find roots and optima in [-1,1]. */

#define EPS 1e-12
static double parg_lb,parg_ub;

static void
rpol_prepare_lb(double *p,i32_t d,i32_t *nr,i32_t *no,double *r,double *o)
{
  i32_t i,ndr,ndo,n;
  double *deriv,*o1;

  if(d==1) {
    o[0]=parg_lb;
    o[1]=parg_ub;
    *no=2;
    if(r==NULL || nr == NULL) return;
    r[0]=-p[0]/p[1];
    if(r[0]<=parg_ub && r[0]>=parg_lb) *nr=1;
    else *nr=0;
    return;
  }
  o[0]=parg_lb;
  deriv=alloca(d*sizeof(*deriv));
  o1=alloca(d*sizeof(*o1));
  for(i=0;i<d;i++) deriv[i]=(i+1)*p[i+1];
  rpol_prepare_lb(deriv,d-1,&ndr,&ndo,o+1,o1);
  *no=ndr+2;
  o[ndr+1]=parg_ub;
  if(r==NULL || nr==NULL) return;
  for(i=0,n=0;i+1<*no;i++) {
    double pv1,pv2;
    pv1=rpol_eval(p,d,o[i],1.0);
    pv2=rpol_eval(p,d,o[i+1],1.0);
    if(pv1==0.0) {
      r[n++]=o[i];
      continue;
    }
    if(pv2==0.0) {
      if(i+1==*no) r[n++]=o[i+1];
      continue;
    }
    if(pv1<0.0) {
      double x,y;

      if(pv2<0.0) continue;
      x=o[i];
      y=o[i+1];
      while(y>x+EPS) {
	double z,pv;

	z=(x+y)/2;
	pv=rpol_eval(p,d,z,1.0);
	if(pv<0.0) {
	  if(z==x) break;
	  x=z;
	} else {
	  if(z==y) break;
	  y=z;
	}
      }
      r[n++]=x;
    } else {
      double x,y;

      if(pv2>0.0) continue;
      x=o[i];
      y=o[i+1];
      while(y>x+EPS) {
	double z,pv;

	z=(x+y)/2;
	pv=rpol_eval(p,d,z,1.0);
	if(pv>0.0) {
	  if(z==x) break;
	  x=z;
	} else {
	  if(z==y) break;
	  y=z;
	}
      }
      r[n++]=x;
    }
  }
  *nr=n;
}

void
get_sieve_report_bounds(sr_bounds,poly,d,a,b,steps,smult,lps)
     double *poly,smult,lps;
     i32_t d,a,b,steps;
     unsigned char **sr_bounds;
{
  double *r1,*o1; /* Roots and local optima with |x|<=|y|. */
  double *o2; /* Local optima with |y|<=|x|. */
  i32_t i,j,nr1,no1,no2;

  if(d<=0) complain("Trying to get sieve report bounds for a poldeg %d\n",d);
  if(a%2!=0)
    complain("Odd number %d of steps in i-direction for poly min\n",a);
  r1=alloca(d*sizeof(*r1));
  o1=alloca((d+1)*sizeof(*o1));
  o2=alloca((d+1)*sizeof(*o2));
  /* Calculate o2, using o1 as a polynomial buffer. */
  for(i=0;i<=d;i++) o1[d-i]=poly[i];
  parg_lb=0;
  parg_ub=steps*b;
  rpol_prepare_lb(o1,d,NULL,&no2,NULL,o2);
  /* Calculate r1 and o1. */
  parg_ub=steps*a/2;
  parg_lb=-(steps*a)/2;
  rpol_prepare_lb(poly,d,&nr1,&no1,r1,o1);
  for(i=0;i<a;i++) {
    for(j=0;j<b;j++) {
      /* Treat a rectangle. */
      i32_t o,u,l,r,k;
      double min_pv,pv,ll,rr;
      /* Untere Ecke. */
      u=j*steps;
      /* Left corner. */
      l=(i-a/2)*steps;
      /* Oben. */
      o=u+steps;
      /* Right corner. */
      r=l+steps;
      /* First, the minimum of the polynomial values at the corners,
         excluding the zero corner. */
      min_pv=fabs(rpol_eval(poly,d,l,o));
      pv=fabs(rpol_eval(poly,d,r,o));
      if(pv<min_pv) min_pv=pv;
      if(u==0) {
	if(l<=0 && r>=0) {
	  pv=fabs(rpol_eval(poly,d,1.0,0.0));
	  if(pv<min_pv) min_pv=pv;
	}
	u=1;
      }
      if(l==0) {
	if(u<=1) {
	  pv=fabs(rpol_eval(poly,d,0.0,1.0));
	  if(pv<min_pv) min_pv=pv;
	}
	l=1;
      }
      if(r==0) {
	if(u<=1) {
	  pv=fabs(rpol_eval(poly,d,0.0,1.0));
	  if(pv<min_pv) min_pv=pv;
	}
	r=-1;
      }
      /* In ll and rr, we hold the edges of the projection of the square
	 to the line $R\times\{1\}$ with projection center 0. This is used when
	 we look at the roots. */
      if(l<0) ll=(double)l/(double)u;
      else ll=(double)l/(double)o;
      if(r>0) rr=(double)r/(double)u;
      else rr=(double)r/(double)o;
      pv=fabs(rpol_eval(poly,d,l,u));
      if(pv<min_pv) min_pv=pv;
      pv=fabs(rpol_eval(poly,d,r,u));
      if(pv<min_pv) min_pv=pv;
      /* Next, the contribution of local minima to min of poly. */
      for(k=0;k<no1;k++) {
	double x;

	x=u*o1[k];
	if(x<(double)l || x>(double)r) continue;
	pv=fabs(rpol_eval(poly,d,x,u));
	if(pv<min_pv) min_pv=pv;
      }
      for(k=0;k<no2;k++) {
	if(l>0) {
	  /* Does the local optimum enter the square trough its left edge? */
	  double y;

	  y=l*o2[k];
	  if(y<(double)u || y>(double)o) continue;
	  pv=fabs(rpol_eval(poly,d,l,y));
	  if(pv<min_pv) min_pv=pv;
	}
	if(r<0) {
	  /* Perhaps the local optimum enters trough right edge of square. */
	  double y;

	  y=r*o2[k];
	  if(y<(double)u || y>(double)o) continue;
	  pv=fabs(rpol_eval(poly,d,r,y));
	  if(pv<min_pv) min_pv=pv;
	}
      }
      /* Finally, the roots. */
      for(k=0;k<nr1;k++) {
	i32_t m;
	if(r1[k]<ll || r1[k]>rr) continue;
	for(m=u;m<o;m++) {
	  int n;
	  n=floor(r1[k]*m);
	  if(n>=l && n<=r) {
	    pv=fabs(rpol_eval(poly,d,n,m));
	    if(pv<min_pv) min_pv=pv;
	  }
	  n++;
	  if(n>=l && n<=r) {
	    pv=fabs(rpol_eval(poly,d,n,m));
	    if(pv<min_pv) min_pv=pv;
	  }
	}
      }
      if(min_pv<=0.0) {
	sr_bounds[j][i]=0;
	continue;
      }
      pv=smult*log(min_pv)-lps;
      if(pv<=0.0) {
	sr_bounds[j][i]=0;
	continue;
      }
      if(pv>=UCHAR_MAX) {
	errprintf("Warning: Found sieve report bound %g\n",pv);
	sr_bounds[j][i]=UCHAR_MAX;
	continue;
      }
      sr_bounds[j][i]=rint(pv);
    }
  }
}
