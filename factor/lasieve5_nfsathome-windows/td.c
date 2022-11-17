
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#include <sys/time.h>
#include "gmp.h"
#include "asm/siever-config.h"
#include "if.h"
#include "primgen32.h"
#include "primgen64.h"


static u32_t td_n, td_nl, td_nbl, td_nb;
static mpz_t *td_tree, td_rem, td_aux0, td_aux1, td_aux2;
static mpz_t td_aux_rec0, td_aux_rec1, td_aux_rec2;
static mpz_t td_aux_r, td_aux_y, td_aux_z;
static size_t td_tree_alloc=0;

static size_t *td_prec, td_prec_alloc=0;

/* ---------------------------------------- */

static void td_adjust(u32_t nf)
{
  size_t na;
  u32_t h;

  td_n=nf;
  for (td_nl=0,h=td_n-1; h; td_nl++) h>>=1;
  na=(size_t)(1<<td_nl);
  if (na<2) na=2; /* case nf<2 */
  if (td_tree_alloc<na) {
    if (td_tree_alloc) td_tree=(mpz_t *)xrealloc(td_tree,na*sizeof(mpz_t));
    else { /* init */
      td_tree=(mpz_t *)xmalloc(na*sizeof(mpz_t));
      mpz_init(td_aux0);
      mpz_init(td_aux1);
      mpz_init(td_aux2);
      mpz_init(td_rem);
      mpz_init(td_aux_r);
      mpz_init(td_aux_y);
      mpz_init(td_aux_z);
      mpz_init(td_aux_rec0);
      mpz_init(td_aux_rec1);
      mpz_init(td_aux_rec2);
    }
    for (; td_tree_alloc<na; td_tree_alloc++) mpz_init(td_tree[td_tree_alloc]);
  }

  if (td_prec_alloc<=(size_t)td_nl) {
    na=(size_t)(td_nl+1);
    if (td_prec_alloc)
      td_prec=(size_t *)xrealloc(td_prec,na*sizeof(size_t));
    else
      td_prec=(size_t *)xmalloc(na*sizeof(size_t));
    td_prec_alloc=na;
  }
}

/* ----------------------- public ----------------------------- */

void td_pprod_get(mpz_t rop, u64_t p0, u64_t p1)
{
  u32_t i, l, p, b;
  mpz_t *td_pprod_aux;
  double db;
  pr32_struct pgs;
  pr64_struct pgs64;
  u64_t p64;

  db=(double)(p1-p0); db/=M_LN2; db/=(double)(mp_bits_per_limb);
  if (db>4294967295.) complain("td_pprod_get: more than 2^32 limbs\n");
  b=(u32_t)db; printf("expsize: %u limbs\n",b);
  for (l=1; (l<32) && (b>16); l++) b=(b+1)/2;
  td_pprod_aux=(mpz_t *)xmalloc(l*sizeof(mpz_t));
  for (i=0; i<l; i++) mpz_init_set_ui(td_pprod_aux[i],1);
  if (p1>>32) {
    initprime64(&pgs64);
    if (mp_bits_per_limb<64) {
      mpz_t td_ull32_aux;

      mpz_init(td_ull32_aux);
      if (mp_bits_per_limb!=32)
        complain("td_pprod_get: mp_bits_per_limb neither 32 nor 64\n");
      for (p64=pr64_seek(&pgs64,p0); (p64<p1) && (p64!=0);
           p64=nextprime64(&pgs64)) {
        mpz_set_ui(td_ull32_aux,(u32_t)(p64>>32));
        mpz_mul_2exp(td_ull32_aux,td_ull32_aux,32);
        mpz_add_ui(td_ull32_aux,td_ull32_aux,(u32_t)p64);
        mpz_mul(td_pprod_aux[0],td_pprod_aux[0],td_ull32_aux);
        if (mpz_size(td_pprod_aux[0])>b) {
          for (i=1; i<l; i++) {
            mpz_mul(td_pprod_aux[i],td_pprod_aux[i],td_pprod_aux[i-1]);
            mpz_set_ui(td_pprod_aux[i-1],1);
            if ((mpz_size(td_pprod_aux[i])>>i)<=b) break;
          }
        }
      }
      mpz_clear(td_ull32_aux);
    } else {
      for (p64=pr64_seek(&pgs64,p0); (p64<p1) && (p64!=0);
           p64=nextprime64(&pgs64)) {
        mpz_mul_ui(td_pprod_aux[0],td_pprod_aux[0],p64);
        if (mpz_size(td_pprod_aux[0])>b) {
          for (i=1; i<l; i++) {
            mpz_mul(td_pprod_aux[i],td_pprod_aux[i],td_pprod_aux[i-1]);
            mpz_set_ui(td_pprod_aux[i-1],1);
            if ((mpz_size(td_pprod_aux[i])>>i)<=b) break;
          }
        }
      }
    }
    clearprime64(&pgs64);
  } else {
    initprime32(&pgs);
    for (p=pr32_seek(&pgs,(u32_t)p0); (p<(u32_t)p1) && (p!=0);
         p=nextprime32(&pgs)) {
      mpz_mul_ui(td_pprod_aux[0],td_pprod_aux[0],p);
      if (mpz_size(td_pprod_aux[0])>b) {
        for (i=1; i<l; i++) {
          mpz_mul(td_pprod_aux[i],td_pprod_aux[i],td_pprod_aux[i-1]);
          mpz_set_ui(td_pprod_aux[i-1],1);
          if ((mpz_size(td_pprod_aux[i])>>i)<=b) break;
        }
      }
    }
    clearprime32(&pgs);
  }

  for (i=1; i<l; i++)
    mpz_mul(td_pprod_aux[i],td_pprod_aux[i],td_pprod_aux[i-1]);
  mpz_set(rop,td_pprod_aux[l-1]);

  for (i=0; i<l; i++) mpz_clear(td_pprod_aux[i]);
  free(td_pprod_aux);
}


void td_pprod_store(mpz_t op, u64_t p0, u64_t p1)
{
  struct stat statbuf;
  char *fname;
  FILE *fi;

  // SMJS asprintf(&fname,"pprod/pprod.%llu-%llu",p0,p1-1);
  asprintf(&fname,"pprod/pprod."UL_FMTSTR"-"UL_FMTSTR,p0,p1-1);
  if (stat(fname,&statbuf)==0)
    complain("td_pprod_store: file already exists: %s\n",fname);
  fi=fopen(fname,"w");
  if (fi==NULL) complain("cannot open %s\n",fname);
  if (!mpz_out_raw(fi,op)) complain("mpz-error writing file\n");
  fclose(fi);
  free(fname);
}


int td_pprod_load(mpz_t rop, u64_t p0, u64_t p1)
{
  struct stat statbuf;
  char *fname;
  FILE *fi;

  // SMJSasprintf(&fname,"pprod/pprod.%llu-%llu",p0,p1-1);
  asprintf(&fname,"pprod/pprod."UL_FMTSTR"-"UL_FMTSTR,p0,p1-1);
  if (stat(fname,&statbuf)) { mpz_set_ui(rop,0); return -1; }
//    complain("td_pprod_load: file does not exist: %s\n",fname);
  fi=fopen(fname,"r");
  if (fi==NULL) complain("cannot open %s\n",fname);
  if (!mpz_inp_raw(rop,fi)) complain("mpz-error reading file\n");
  fclose(fi);
  free(fname);
  return 0;
}

#define THRESHOLD   256

void td_init(mpz_t *op, u32_t nf)
{
  u32_t l, n, nn, i, d;

  td_adjust(nf);
  d=1<<td_nl;
  if (td_n<THRESHOLD) {
    d/=2;
    nn=td_n>>1;
    for (i=0; i<nn; i++)
      mpz_mul(td_tree[d+i],op[2*i],op[2*i+1]);
    if (td_n&1) { mpz_set(td_tree[d+nn],op[td_n-1]); nn++; }
    n=nn;
    for (; nn<d; nn++) mpz_set_ui(td_tree[d+nn],0);
    for (l=1; l<td_nl; l++) {
      d>>=1;
      nn=n>>1;
      for (i=0; i<nn; i++)
        mpz_mul(td_tree[d+i],td_tree[2*d+2*i],td_tree[2*d+2*i+1]);
      if (n&1) { mpz_set(td_tree[d+nn],td_tree[2*d+n-1]); nn++; }
      n=nn;
      for (; nn<d; nn++) mpz_set_ui(td_tree[d+nn],0);
    }
    if (nf==1) { mpz_set(td_tree[1],op[0]); d=1; n=1; }
  } else {
    u32_t d6, ii, in, il, id;

    td_nbl=td_n>>(td_nl-6); if ((td_nbl<<(td_nl-6))<td_n) td_nbl++;
    d>>=1;
    d6=d>>5; ii=0; nn=td_n;
    for (il=0; il<d6; il++) {
      if (!nn) {
        mpz_set_ui(td_tree[d6+il],1);
        continue;
      }
      id=(d6+il)<<5;
      if (nn>=td_nbl) n=td_nbl; else n=nn;
      for (i=0; 2*i+1<n; i++) {
        mpz_mul(td_tree[id+i],op[ii],op[ii+1]); ii+=2;
      }
      if (2*i<n) mpz_set(td_tree[id+i++],op[ii++]);
      for (; i<(1<<5); i++) mpz_set_ui(td_tree[id+i],1);
      for (l=0,in=(1<<5); l<5; l++) {
        id>>=1; in>>=1;
        for (i=0; i<in; i++)
          mpz_mul(td_tree[id+i],td_tree[2*id+2*i],td_tree[2*id+2*i+1]);
      }
      if (id!=d6+il) Schlendrian("td_inita\n");
      nn-=n;
    }
    if (ii!=td_n) Schlendrian("td_inite\n");

    d=d6;
    for (l=6; l<td_nl; l++) {
      d>>=1;
      for (i=0; i<d; i++)
        mpz_mul(td_tree[d+i],td_tree[2*d+2*i],td_tree[2*d+2*i+1]);
    }
    n=1;
  }
  if ((d!=1) || (n!=1)) complain("td_initc\n");
  if (!mpz_sgn(td_tree[1])) complain("td_initd\n");
  mpz_set_ui(td_rem,1);
}


void td_execute(mpz_t pprod)
{
  size_t sp, st, s, sb;

  sp=mpz_size(pprod); st=mpz_size(td_tree[1]);
printf("sizes %zu %zu\n",sp,st);
  if (sp>2*st) {
    sb=st*mp_bits_per_limb;
    mpz_set_ui(td_aux_r,1);
    mpz_mul_2exp(td_aux_r,td_aux_r,2*sb);
    mpz_mod(td_aux_r,td_aux_r,td_tree[1]);
    s=sp/st; s--;
    mpz_fdiv_q_2exp(td_aux0,pprod,s*sb);
    s--;
    while (1) {
      mpz_fdiv_r_2exp(td_aux_z,td_aux0,sb);
      mpz_fdiv_q_2exp(td_aux_y,td_aux0,sb);
      mpz_mul(td_aux0,td_aux_y,td_aux_r);
      if (mpz_sgn(td_aux_z)) {
        mpz_mul_2exp(td_aux_z,td_aux_z,sb);
        memcpy(td_aux_z[0]._mp_d,pprod[0]._mp_d+s*st,st*sizeof(mp_limb_t));
      } else { /* extremly rare */
        mpz_set_ui(td_aux_z,1);
        mpz_mul_2exp(td_aux_z,td_aux_z,sb);
        mpz_set(td_aux_y,td_aux_z);
        memcpy(td_aux_z[0]._mp_d,pprod[0]._mp_d+s*st,st*sizeof(mp_limb_t));
        mpz_sub(td_aux_z,td_aux_z,td_aux_y);
      }
      mpz_add(td_aux0,td_aux0,td_aux_z);
      if (!s) break;
      s--;
    }
    mpz_mod(td_aux1,td_aux0,td_tree[1]);
  } else {
    mpz_fdiv_r(td_aux1,pprod,td_tree[1]);
  }
  mpz_mul(td_rem,td_rem,td_aux1);
  mpz_mod(td_rem,td_rem,td_tree[1]);
}


void td_end(mpz_t *rop, mpz_t *op)
{
  u32_t l, i, d;

  mpz_set(td_tree[1],td_rem);
  if (td_n<THRESHOLD) {
    for (l=1,d=1; l<td_nl; l++) {
      for (i=0; i<d; i++) {
        if (!mpz_sgn(td_tree[2*d+2*i+1])) break;
        mpz_mod(td_tree[2*d+2*i],td_tree[d+i],td_tree[2*d+2*i]);
        mpz_mod(td_tree[2*d+2*i+1],td_tree[d+i],td_tree[2*d+2*i+1]);
      }
      if ((i<d) && (mpz_sgn(td_tree[2*d+2*i])))
        mpz_mod(td_tree[2*d+2*i],td_tree[d+i],td_tree[2*d+2*i]);
      d<<=1;
    }
    if (td_nl) {
      if (d!=(1<<(td_nl-1))) complain("td_end\n");
      for (i=0; 2*i+1<td_n; i++) {
        mpz_mod(rop[2*i],td_tree[d+i],op[2*i]);
        mpz_mod(rop[2*i+1],td_tree[d+i],op[2*i+1]);
      }
      if (2*i+1==td_n) mpz_mod(rop[2*i],td_tree[d+i],op[2*i]);
    } else mpz_mod(rop[0],td_tree[1],op[0]);
  } else {
    u32_t ii, in, il, id, n, nn;
//  clock_t time_begin;

    td_nbl=td_n>>(td_nl-6); if ((td_nbl<<(td_nl-6))<td_n) td_nbl++;

//time_begin=clock();

    for (l=1,d=1; l+5<td_nl; l++) {
      for (i=0; i<d; i++) {
        mpz_mod(td_tree[2*d+2*i],td_tree[d+i],td_tree[2*d+2*i]);
        mpz_mod(td_tree[2*d+2*i+1],td_tree[d+i],td_tree[2*d+2*i+1]);
      }
      d<<=1;
    }
//printf("endtime %.2f\n",((double)(clock()-time_begin))/CLOCKS_PER_SEC);

    ii=0; nn=td_n;
    for (il=0; il<d; il++) {
      if (!nn) break;
      id=d+il;
      if (nn>=td_nbl) n=td_nbl; else n=nn;
      for (l=0,in=1; l<5; l++) {
        for (i=0; i<in; i++) {
          mpz_mod(td_tree[2*id+2*i],td_tree[id+i],td_tree[2*id+2*i]);
          mpz_mod(td_tree[2*id+2*i+1],td_tree[id+i],td_tree[2*id+2*i+1]);
        }
        id<<=1; in<<=1;
      }

      for (i=0; 2*i+1<n; i++) {
        mpz_mod(rop[ii],td_tree[id+i],op[ii]);
        mpz_mod(rop[ii+1],td_tree[id+i],op[ii+1]); ii+=2;
      }
      if (2*i+1==n) { mpz_mod(rop[ii],td_tree[id+i],op[ii]); ii++; }

      nn-=n;
    }
  }
  for (i=0; i<td_n; i++) mpz_gcd(rop[i],rop[i],op[i]);
}

/* ------ scaled remainder trees ---------- */

void td2_init(mpz_t *op, u32_t nf)
{
  u32_t l, n, nn, i, d;
  size_t m, nb;

  td_adjust(nf);
  d=1<<td_nl;
  if (td_n<THRESHOLD) {
    d/=2;
    nn=td_n>>1;
    for (i=0; i<nn; i++)
      mpz_mul(td_tree[d+i],op[2*i],op[2*i+1]);
    if (td_n&1) { mpz_set(td_tree[d+nn],op[td_n-1]); nn++; }
    n=nn;
    for (; nn<d; nn++) mpz_set_ui(td_tree[d+nn],0);
    for (l=1; l<td_nl; l++) {
      d>>=1;
      nn=n>>1;
      for (i=0; i<nn; i++)
        mpz_mul(td_tree[d+i],td_tree[2*d+2*i],td_tree[2*d+2*i+1]);
      if (n&1) { mpz_set(td_tree[d+nn],td_tree[2*d+n-1]); nn++; }
      n=nn;
      for (; nn<d; nn++) mpz_set_ui(td_tree[d+nn],0);
    }
    if (nf==1) { mpz_set(td_tree[1],op[0]); d=1; n=1; }
  } else {
    u32_t d6, ii, in, il, id;

    td_nbl=td_n>>(td_nl-6); if ((td_nbl<<(td_nl-6))<td_n) td_nbl++;
    d>>=1;
    d6=d>>5; ii=0; nn=td_n;
    for (il=0; il<d6; il++) {
      if (!nn) break;
      id=(d6+il)<<5;
      if (nn>=td_nbl) n=td_nbl; else n=nn;
      for (i=0; 2*i+1<n; i++) {
        mpz_mul(td_tree[id+i],op[ii],op[ii+1]); ii+=2;
      }
      if (2*i<n) mpz_set(td_tree[id+i++],op[ii++]);
      for (; i<(1<<5); i++) mpz_set_ui(td_tree[id+i],1);
      for (l=0,in=(1<<5); l<5; l++) {
        id>>=1; in>>=1;
        for (i=0; i<in; i++)
          mpz_mul(td_tree[id+i],td_tree[2*id+2*i],td_tree[2*id+2*i+1]);
      }
      if (id!=d6+il) Schlendrian("td_inita\n");
      nn-=n;
    }
    if (ii!=td_n) Schlendrian("td_inite\n");
    if (il<d6) {
      u32_t i0, i1;

      for (l=0,i0=il+d6,i1=d6+d6; l<6; l++,i0<<=1,i1<<=1)
        for (i=i0; i<i1; i++)
          mpz_set_ui(td_tree[i],1);
    }

    d=d6;
    for (l=6; l<td_nl; l++) {
      d>>=1;
      for (i=0; i<d; i++)
        mpz_mul(td_tree[d+i],td_tree[2*d+2*i],td_tree[2*d+2*i+1]);
    }
    n=1;
  }
  if ((d!=1) || (n!=1)) complain("td_initc\n");
  if (!mpz_sgn(td_tree[1])) complain("td_initd\n");
  mpz_set_ui(td_rem,1);
  for (l=0,d=1; l<td_nl; l++) {
    for (i=0,m=0; i<d; i++) {
      nb=mpz_sizeinbase(td_tree[d+i],2);
      if (nb>m) m=nb;
    }
    td_prec[l]=m;
//printf("level %u: %zu\n",l,m);
    d<<=1;
  }
  for (l=td_nl,m=0; l; m++) l>>=1;
  for (l=0; l<td_nl; l++) td_prec[l]+=(3+m+100);  /* ???? */
}


void td2_reciprocal(mpz_t quot, mpz_t rem, mpz_t op, size_t nb)
/* 2^(2*nb)=quot*op+rem, 0<=rem<op, condition op<2^nb */
{
  clock_t time_begin;
#if 1
  if (nb&31) Schlendrian("td2_reciprocal\n");
time_begin=clock();
  mpz_set_ui(td_aux_rec0,1);
  mpz_mul_2exp(td_aux_rec0,td_aux_rec0,nb);
  mpz_fdiv_q_2exp(td_aux0,op,nb/2); mpz_add_ui(td_aux0,td_aux0,1);
  mpz_fdiv_qr(td_aux_rec1,td_aux_rec2,td_aux_rec0,td_aux0);
//printf("%.2f\n",((double)(clock()-time_begin))/CLOCKS_PER_SEC);
  mpz_mul_2exp(td_aux1,td_aux0,nb/2); mpz_sub(td_aux1,td_aux1,op);
  mpz_mul(td_aux1,td_aux1,td_aux_rec1);
  mpz_fdiv_q_2exp(td_aux1,td_aux1,nb/2);
  mpz_add(td_aux1,td_aux1,td_aux_rec2);
  mpz_mul(td_aux1,td_aux1,td_aux_rec1);
  mpz_fdiv_q_2exp(td_aux1,td_aux1,nb/2);
  mpz_mul_2exp(td_aux_rec1,td_aux_rec1,nb/2);
  mpz_add(td_aux_rec1,td_aux_rec1,td_aux1); /* almost quot */

  mpz_set_ui(td_aux_rec0,1);
  mpz_mul_2exp(td_aux_rec0,td_aux_rec0,2*nb);
  mpz_mul(td_aux0,td_aux_rec1,op);
  mpz_sub(td_aux_rec0,td_aux_rec0,td_aux0);
  mpz_fdiv_qr(quot,rem,td_aux_rec0,op);
  mpz_add(quot,quot,td_aux_rec1);

#else
  mpz_set_ui(td_aux_rec0,1);
  mpz_mul_2exp(td_aux_rec0,td_aux_rec0,2*nb);
  mpz_fdiv_qr(quot,rem,td_aux_rec0,op);
#endif
}


void td2_execute(mpz_t pprod)
{
  size_t sp, st, s, sb;
  clock_t time_begin;

  sp=mpz_size(pprod); st=mpz_size(td_tree[1]);
//printf("sizes %zu %zu\n",sp,st);
  if (td_prec[0]%mp_bits_per_limb)
    td_prec[0]+=(mp_bits_per_limb-(td_prec[0]%mp_bits_per_limb));
  if (sp>2*st) {
    if (st*mp_bits_per_limb>td_prec[0]) Schlendrian("td2_execute: prec 0\n");
    sb=st*mp_bits_per_limb;
sb=td_prec[0];
time_begin=clock();
    td2_reciprocal(td_aux2,td_aux_r,td_tree[1],sb);
//printf("%.2f\n",((double)(clock()-time_begin))/CLOCKS_PER_SEC);
/*
    mpz_set_ui(td_aux_r,1);
    mpz_mul_2exp(td_aux_r,td_aux_r,2*sb);
    mpz_fdiv_qr(td_aux2,td_aux_r,td_aux_r,td_tree[1]);
*/
    s=sp/st; s--;
    mpz_fdiv_q_2exp(td_aux0,pprod,s*sb);
    s--;
    while (1) {
      mpz_fdiv_r_2exp(td_aux_z,td_aux0,sb);
      mpz_fdiv_q_2exp(td_aux_y,td_aux0,sb);
      mpz_mul(td_aux0,td_aux_y,td_aux_r);
      if (mpz_sgn(td_aux_z)) {
        mpz_mul_2exp(td_aux_z,td_aux_z,sb);
        memcpy(td_aux_z[0]._mp_d,pprod[0]._mp_d+s*st,st*sizeof(mp_limb_t));
      } else { /* extremly rare */
        mpz_set_ui(td_aux_z,1);
        mpz_mul_2exp(td_aux_z,td_aux_z,sb);
        mpz_set(td_aux_y,td_aux_z);
        memcpy(td_aux_z[0]._mp_d,pprod[0]._mp_d+s*st,st*sizeof(mp_limb_t));
        mpz_sub(td_aux_z,td_aux_z,td_aux_y);
      }
      mpz_add(td_aux0,td_aux0,td_aux_z);
      if (!s) break;
      s--;
    }
    mpz_fdiv_q_2exp(td_aux1,td_aux0,td_prec[0]);
    mpz_mul(td_aux1,td_aux1,td_aux2);
    mpz_fdiv_q_2exp(td_aux1,td_aux1,td_prec[0]);
    mpz_mul(td_aux1,td_aux1,td_tree[1]);
    mpz_sub(td_aux0,td_aux0,td_aux1);

    mpz_mod(td_aux1,td_aux0,td_tree[1]); /* replace by subtractions ??? */



    mpz_mul(td_aux1,td_aux1,td_aux2);
    mpz_fdiv_q_2exp(td_tree[1],td_aux1,td_prec[0]);


  } else {
complain("todo\n");
    mpz_fdiv_r(td_aux1,pprod,td_tree[1]);
  mpz_mul(td_rem,td_rem,td_aux1);
  mpz_mod(td_rem,td_rem,td_tree[1]);
if (mpz_sizeinbase(td_rem,2)>td_prec[0]-3) Schlendrian("td2_execute\n");
  mpz_mul_2exp(td_aux1,td_rem,td_prec[0]);
  mpz_fdiv_q(td_tree[1],td_aux1,td_tree[1]);

  }
}


void td2_end(mpz_t *rop, mpz_t *op)
{
  u32_t l, i, d;
  int c;
  size_t sh;

  if (td_n<THRESHOLD) {
complain("todo td2_end\n");
    for (l=1,d=1; l<td_nl; l++) {
      for (i=0; i<d; i++) {
        if (!mpz_sgn(td_tree[2*d+2*i+1])) break;
        mpz_mod(td_tree[2*d+2*i],td_tree[d+i],td_tree[2*d+2*i]);
        mpz_mod(td_tree[2*d+2*i+1],td_tree[d+i],td_tree[2*d+2*i+1]);
      }
      if ((i<d) && (mpz_sgn(td_tree[2*d+2*i])))
        mpz_mod(td_tree[2*d+2*i],td_tree[d+i],td_tree[2*d+2*i]);
      d<<=1;
    }
    if (td_nl) {
      if (d!=(1<<(td_nl-1))) complain("td_end\n");
      for (i=0; 2*i+1<td_n; i++) {
        mpz_mod(rop[2*i],td_tree[d+i],op[2*i]);
        mpz_mod(rop[2*i+1],td_tree[d+i],op[2*i+1]);
      }
      if (2*i+1==td_n) mpz_mod(rop[2*i],td_tree[d+i],op[2*i]);
    } else mpz_mod(rop[0],td_tree[1],op[0]);
  } else {
    u32_t ii, in, il, id, n, nn;
    size_t sh;

    td_nbl=td_n>>(td_nl-6); if ((td_nbl<<(td_nl-6))<td_n) td_nbl++;

{
  clock_t time_begin;

  time_begin=clock();
    for (l=1,d=1; l+5<td_nl; l++) {
      sh=td_prec[l-1]-td_prec[l];
      for (i=0; i<d; i++) {
        mpz_mul(td_aux0,td_tree[d+i],td_tree[2*d+2*i+1]);
        mpz_mul(td_aux1,td_tree[d+i],td_tree[2*d+2*i]);
        mpz_fdiv_r_2exp(td_tree[2*d+2*i],td_aux0,td_prec[l-1]);
        mpz_fdiv_q_2exp(td_tree[2*d+2*i],td_tree[2*d+2*i],sh);
        mpz_fdiv_r_2exp(td_tree[2*d+2*i+1],td_aux1,td_prec[l-1]);
        mpz_fdiv_q_2exp(td_tree[2*d+2*i+1],td_tree[2*d+2*i+1],sh);
      }
      d<<=1;
    }
//printf("%.2f\n",((double)(clock()-time_begin))/CLOCKS_PER_SEC);
}

    ii=0; nn=td_n;
    for (il=0; il<d; il++) {
      if (!nn) break;
      id=d+il;
      if (nn>=td_nbl) n=td_nbl; else n=nn;
      for (l=0,in=1; l<5; l++) {
        sh=td_prec[td_nl-5+l-1]-td_prec[td_nl-5+l];
        for (i=0; i<in; i++) {
          mpz_set(td_aux0,td_tree[2*id+2*i]);
          mpz_mul(td_tree[2*id+2*i],td_tree[id+i],td_tree[2*id+2*i+1]);
          mpz_fdiv_r_2exp(td_tree[2*id+2*i],td_tree[2*id+2*i],
                          td_prec[td_nl-5+l-1]);
          mpz_fdiv_q_2exp(td_tree[2*id+2*i],td_tree[2*id+2*i],sh);

          mpz_mul(td_tree[2*id+2*i+1],td_tree[id+i],td_aux0);
          mpz_fdiv_r_2exp(td_tree[2*id+2*i+1],td_tree[2*id+2*i+1],
                          td_prec[td_nl-5+l-1]);
          mpz_fdiv_q_2exp(td_tree[2*id+2*i+1],td_tree[2*id+2*i+1],sh);
        }
        id<<=1; in<<=1;
      }

      sh=td_prec[td_nl-1];
      for (i=0; 2*i+1<n; i++) {
        mpz_mul(td_aux0,td_tree[id+i],op[ii]);
        mpz_mul(td_aux0,td_aux0,op[ii+1]);
        c=mpz_tstbit(td_aux0,td_prec[td_nl-1]-1);
        mpz_fdiv_q_2exp(td_aux0,td_aux0,td_prec[td_nl-1]);
        if (c) mpz_add_ui(td_aux0,td_aux0,1);
        mpz_mod(rop[ii],td_aux0,op[ii]);
        mpz_mod(rop[ii+1],td_aux0,op[ii+1]);

        ii+=2;
      }
      if (2*i+1==n) {
        mpz_mul(td_aux0,td_tree[id+i],op[ii]);
        c=mpz_tstbit(td_aux0,td_prec[td_nl-1]-1);
        mpz_fdiv_q_2exp(td_aux0,td_aux0,td_prec[td_nl-1]);
        if (c) mpz_add_ui(td_aux0,td_aux0,1);
        mpz_mod(rop[ii],td_aux0,op[ii]);

        ii++;
      }

      nn-=n;
    }
  }
  for (i=0; i<td_n; i++) mpz_gcd(rop[i],rop[i],op[i]);

}


void mpout(mpz_t x)
{
  gmp_printf("%Zu\n",x);
}
