/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/


#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifndef _WIN64 // SMJS
#include <sys/times.h>
#endif
#include <sys/time.h>

#define uchar         unsigned char

#include "gmp.h"
#include "asm/siever-config.h"
#include "if.h"
#include "asm/montgomery_mul.h"
#include "pm1.h"


#define BUF_INC    256

extern ulong montgomery_inv_n;
extern ulong *montgomery_modulo_n;
extern ulong montgomery_modulo_R2[NMAX_ULONGS];
extern size_t montgomery_ulongs;


extern void gcd(ulong *,ulong *,ulong *);
extern int asm_invert(ulong *,ulong *);



static int pm1_is_init;
mpz_t gmp_f;

static uchar *pm1_prime_bit;
static u32_t pm1_prime_max=0;

static u32_t pm1_index1, pm1_index2;

typedef struct {
  u32_t B1;
  u32_t len;
  uchar *ex;
} pm1_scheme1_t;

static pm1_scheme1_t *B1_scheme;
static size_t B1_scheme_len=0, B1_scheme_alloc=0;

typedef struct {
  u32_t B1;
  u32_t B2;
  u32_t D;
  uchar *rpc;
  uchar *tests;
  u32_t tablen;
  u32_t *tab;
  u32_t beginD;
  u32_t endD;
} pm1_scheme2_t;

static pm1_scheme2_t *B2_scheme;
static size_t B2_scheme_len=0, B2_scheme_alloc=0;

static ulong *mm_stack;
static size_t mm_stack_alloc=0;

static ulong **mm_B2_tab1, **mm_B2_tab2;
static u32_t mm_B2_tab1_len=0, mm_B2_tab2_len=0;


/* ------------------- conversion ------------------- */

static void mpz_set_v(mpz_t targ, ulong *v)
{
  size_t i;

  mpz_set_ui(targ,0);
  for (i=0; i<montgomery_ulongs; i++) {
    mpz_mul_2exp(targ,targ,mp_bits_per_limb);
    mpz_add_ui(targ,targ,v[montgomery_ulongs-1-i]);
  }
}


static void mpz_get_v(ulong *v, mpz_t src)
{
  size_t i;

  asm_zero(v);
  for (i=0; i<mpz_size(src); i++) v[i]=mpz_getlimbn(src,i);
}

/* ------------------- bit manipulation ------------------- */

static int get_bit(uchar *rop, u32_t pos)
{
  return (int)((rop[pos>>3])&(1<<(pos&7)));
}


static void set_bit(uchar *rop, u32_t pos)
{
  uchar x, y;

  x=1<<(pos&7);
  y=rop[pos>>3];
  if (!(y&x)) {
    y^=x;
    rop[pos>>3]=y;
  }
}


static void clear_bit(uchar *rop, u32_t pos)
{
  uchar x, y;

  x=1<<(pos&7);
  y=rop[pos>>3];
  if ((y&x)) {
    y^=x;
    rop[pos>>3]=y;
  }
}

static int get_bit64(u64_t *rop, u32_t pos)
{
  return (int)((rop[pos>>3])&(1<<(pos&7)));
}


static void set_bit64(u64_t *rop, u32_t pos)
{
  uchar x, y;

  x=1<<(pos&7);
  y=rop[pos>>3];
  if (!(y&x)) {
    y^=x;
    rop[pos>>3]=y;
  }
}


static void clear_bit64(u64_t *rop, u32_t pos)
{
  uchar x, y;

  x=1<<(pos&7);
  y=rop[pos>>3];
  if ((y&x)) {
    y^=x;
    rop[pos>>3]=y;
  }
}
static uchar pop8_tab[256]={
0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

static u32_t popcount(u64_t *v, u32_t len)
{
  u32_t res, i;
  uchar *u;

  u=(uchar *)v;
  for (i=0,res=0; i<8*len; i++) res+=(u32_t)(pop8_tab[u[i]]);
  return res;
}

/* ------------------------ init primes ------------------------- */

/* ok, if limit is <2^32-2^17 or something like that */
static void pm1_init_primes(u32_t limit)
{
  u32_t i, j, l;

  if (limit>(1UL<<31)) complain("pm1_init_primes: too big\n");
  if (pm1_prime_max) free(pm1_prime_bit);
  l=(limit+15)>>4;
  pm1_prime_max=l<<4;
/* pm1_prime_bit[i] stores odd primes in [16i+1,16i+15] */
  pm1_prime_bit=(uchar *)xmalloc(l*sizeof(uchar));
  memset(pm1_prime_bit,255,l*sizeof(uchar));
  pm1_prime_bit[0]=254; /* 1 is not prime */
  j=3;
  do {
    for (i=j*j; i<pm1_prime_max; i+=2*j) clear_bit(pm1_prime_bit,i>>1);
    while (j<pm1_prime_max) {
      j+=2;
      if (get_bit(pm1_prime_bit,j>>1)) break;
    }
    if (j>65535) j=65535;
  } while (j*j<pm1_prime_max);
}

/* ------------------------ init phase 1 ---------------------- */

static void pm1_init_B1(u32_t B1)
{
  u32_t i, p, q, len, l8;
  mpz_t B1_aux;

  if (pm1_prime_max<B1) pm1_init_primes(B1);
  adjust_bufsize((void **)(&B1_scheme),&B1_scheme_alloc,
                 B1_scheme_len+1,BUF_INC,sizeof(pm1_scheme1_t));
  B1_scheme[B1_scheme_len].B1=B1;
  mpz_init_set_ui(B1_aux,1);
  for (q=2; q<B1; q*=2) mpz_add(B1_aux,B1_aux,B1_aux);
  for (i=1; 2*i+1<B1; i++)
    if (get_bit(pm1_prime_bit,i)) {
      p=2*i+1;
      for (q=p; q<B1; q*=p) mpz_mul_ui(B1_aux,B1_aux,(ulong)p);
    }
  len=(u32_t)mpz_sizeinbase(B1_aux,2)-1;
  l8=(len+7)/8;
  B1_scheme[B1_scheme_len].len=len;
  B1_scheme[B1_scheme_len].ex=(uchar *)xcalloc((size_t)l8,sizeof(uchar));
  for (i=0; i<len; i++)
    if (mpz_tstbit(B1_aux,(ulong)(len-1-i)))
      set_bit(B1_scheme[B1_scheme_len].ex,i);
  mpz_clear(B1_aux);

#ifdef PM1_STAT	
/* STAT */
  printf("# STAT: step1 %u muls\n",B1_scheme[B1_scheme_len].len);
#endif
  B1_scheme_len++;
}

/* ------------------------ init phase 2 ---------------------- */

static u32_t find_triples(uchar **rop, uchar **pr, u32_t dim1, u32_t dim2)
{
  u32_t i, j, d64, *cnt, ind, res, l;
  u64_t **bit_p, *r1, *r2;
  uchar *tmp;

  cnt=(u32_t *)xmalloc(dim1*sizeof(u32_t));
  d64=(dim2+63)/64;
  bit_p=(u64_t **)xmalloc(dim1*sizeof(u64_t *));
  for (i=0; i<dim1; i++) {
    bit_p[i]=(u64_t *)xcalloc((size_t)d64,sizeof(u64_t));
    for (j=0; j<dim2; j++) if (pr[i][j]) set_bit64(bit_p[i],j);
    cnt[i]=popcount(bit_p[i],d64);
  }
  r1=(u64_t *)xcalloc((size_t)d64,sizeof(u64_t));
  r2=(u64_t *)xcalloc((size_t)d64,sizeof(u64_t));

  res=dim2; /* squarings */
  tmp=(uchar *)xmalloc(2*dim1*dim2*sizeof(uchar *)); tmp[0]=0;
  ind=1; l=0;
  while (1) {
    u32_t bi0, bi1, bi2, bn, n;

/* todo: improve triple search algorithm */
    bi0=0; bn=cnt[0];
    for (i=1; i<dim1; i++) if (cnt[i]>bn) { bi0=i; bn=cnt[i]; }
    bn=0; bi1=bi0;
    for (i=0; i<dim1; i++)
      if (i!=bi0) {
        for (j=0; j<d64; j++) r1[j]=bit_p[bi0][j]&bit_p[i][j];
        n=popcount(r1,d64);
        if (n>bn) { bi1=i; bn=n; }
      }
    for (j=0; j<d64; j++) r1[j]=bit_p[bi0][j]&bit_p[bi1][j];
    bn=0; bi2=bi0;
    for (i=0; i<dim1; i++)
      if ((i!=bi0) && (i!=bi1)) {
        for (j=0; j<d64; j++) r2[j]=r1[j]&bit_p[i][j];
        n=popcount(r2,d64);
        if (n>bn) { bi2=i; bn=n; }
      }
    if (bn<5) break;

    if (dim1>254) {
      tmp[ind++]=(uchar)(bi0>>8); tmp[ind++]=(uchar)(bi0&0xff);
      tmp[ind++]=(uchar)(bi1>>8); tmp[ind++]=(uchar)(bi1&0xff);
      tmp[ind++]=(uchar)(bi2>>8); tmp[ind++]=(uchar)(bi2&0xff);
    } else {
      tmp[ind++]=(uchar)bi0; tmp[ind++]=(uchar)bi1; tmp[ind++]=(uchar)bi2;
    }
    res+=4;
    for (j=0; j<dim2; j++)
      if (pr[bi0][j]&pr[bi1][j]&pr[bi2][j]) {
        if (!(get_bit64(bit_p[bi0],j)&get_bit64(bit_p[bi1],j)&get_bit64(bit_p[bi2],j)))
          complain("bits\n");
        pr[bi0][j]=0; pr[bi1][j]=0; pr[bi2][j]=0;
        cnt[bi0]--; cnt[bi1]--; cnt[bi2]--;
        clear_bit64(bit_p[bi0],j);
        clear_bit64(bit_p[bi1],j);
        clear_bit64(bit_p[bi2],j);
        if (dim2>254) {
          tmp[ind++]=(uchar)(j>>8); tmp[ind++]=(uchar)(j&0xff);
        } else tmp[ind++]=(uchar)j;
        res+=2;
      }
    if (dim2>254) tmp[ind++]=255;
    tmp[ind++]=255;
  }
  if (dim1>254) tmp[ind++]=255;
  tmp[ind++]=255;
  tmp=xrealloc(tmp,ind*sizeof(uchar));
  *rop=tmp;
  for (i=0; i<dim1; i++) res+=cnt[i]; /* remaining tests */
  free(cnt);
  for (i=0; i<dim1; i++) free(bit_p[i]);
  free(bit_p);
  free(r1);
  free(r2);
  return res;
}


static u32_t create_B2_scheme(pm1_scheme2_t *s, u32_t B1, u32_t B2, u32_t d, int verb)
{
  u32_t nmul, n, i, b0, b1, p;
  u32_t dim1, dim2, i1, i2, l, *tabnr;
  u32_t n0, n1;
  uchar **pr, **prt, *tmp, *tmpt;

  if ((!d) || (d%30)) Schlendrian("create_B2_scheme\n");
  if (s->D) { free(s->rpc); free(s->tests); free(s->tab); }
  s->D=d;
  nmul=19; /* table up to 30 */
  if (d>30)
    nmul+=(19+16*(d/60-1)); /* table up to D */
  n=(d*d)>>1; while (n) { n>>=1; nmul++; }  /* exponentiation D^2 */

  if (2*B1<d) b0=0; else b0=(B1-d/2)/d;
  n=b0>>1; while (n) { n>>=1; nmul+=2; }  /* 2 exponentiations b0 */
  b1=1+(B2+d/2)/d;
  nmul+=2*(b1-b0); /* table of multiples of D */
  s->beginD=b0;
  s->endD=b1;

  dim1=b1-b0; dim2=d/4;
  pr=(uchar **)xmalloc(dim1*sizeof(uchar *));
  for (i=0; i<dim1; i++) pr[i]=(uchar *)xcalloc((size_t)dim2,sizeof(uchar));
  for (i=B1/2; 2*i+1<B2; i++)
    if (get_bit(pm1_prime_bit,i)) {
      n=2*i+1-b0*d;
      i1=n/d; i2=n%d;
      if (2*i2>d) { i1++; i2=d-i2; }
      i2/=2;
      if ((i1>=dim1) || (i2>=dim2)) Schlendrian("create_B2_scheme %u\n",2*i+1);
      pr[i1][i2]=1;
    }
  tabnr=(u32_t *)xmalloc(dim2*sizeof(u32_t));
  for (i2=0,n=0; i2<dim2; i2++) {
    for (i1=0; i1<dim1; i1++) if (pr[i1][i2]) break;
    if (i1<dim1) tabnr[i2]=n++; else tabnr[i2]=U32_MAX;
  }
  s->tablen=n;
  s->tab=(u32_t *)xmalloc(n*sizeof(u32_t));
  for (i=0; i<dim2; i++) if (tabnr[i]<U32_MAX) s->tab[tabnr[i]]=2*i+1;

/* abort if table lengths are unbalanced or >=2^16-1: */
  if ((n>4*dim1) || (dim1>4*n) || (n>=0xffff) || (dim1>=0xffff)) {
    free(tabnr);
    for (i=0; i<dim1; i++) free(pr[i]);
    free(pr);
    s->rpc=xmalloc(1);
    s->tests=xmalloc(1);
    return 0;
  }

/* compress pr */
  for (i2=0,n=0; i2<dim2; i2++) {
    if (tabnr[i2]==U32_MAX) continue;
    if (i2==n) { n++; continue; }
    for (i1=0; i1<dim1; i1++) pr[i1][n]=pr[i1][i2];
    n++;
  }
  if (n!=s->tablen) Schlendrian("create_B2_scheme\n");
  dim2=n;

/* rpc */
  prt=(uchar **)xmalloc(dim2*sizeof(uchar *));
  for (i=0; i<dim2; i++) prt[i]=(uchar *)xcalloc((size_t)dim1,sizeof(uchar));
  for (i1=0; i1<dim1; i1++)
    for (i2=0; i2<dim2; i2++)
      prt[i2][i1]=pr[i1][i2];
  n0=find_triples(&tmp,pr,dim1,dim2);
  n1=find_triples(&tmpt,prt,dim2,dim1);
  if (n0<n1) {
    nmul+=n0;
    s->rpc=tmp;
    free(tmpt);
    s->rpc[0]=0;
  } else {
    nmul+=n1;
    s->rpc=tmpt;
    free(tmp);
    s->rpc[0]=1;
    for (i1=0; i1<dim1; i1++)
      for (i2=0; i2<dim2; i2++)
        pr[i1][i2]=prt[i2][i1];
  }

  for (i=0; i<dim2; i++) free(prt[i]);
  free(prt);

  tmp=(uchar *)xmalloc((dim1*dim2+dim1)*sizeof(uchar));
  n=0;
  for (i1=0; i1<dim1; i1++) {
    l=0;
    for (i2=0; i2<dim2; i2++)
      if (pr[i1][i2]) {
        i=i2+1-l;
/* for B2<2^16 dim2 should always be <255 */
        while (i>255) { tmp[n++]=255; i-=255; nmul++; }
        tmp[n++]=(uchar)i; /* nmul++; */
        l=i2+1;
      }
    tmp[n++]=0;
  }
  s->tests=(uchar *)xmalloc(n*sizeof(uchar));
  memcpy(s->tests,tmp,n*sizeof(uchar));
  free(tmp);

  free(tabnr);
  for (i=0; i<dim1; i++) free(pr[i]);
  free(pr);
  return nmul;
}


static void pm1_init_B2(u32_t B1, u32_t B2)
{
  u32_t d, nmul, best, bestd, len;

  if (pm1_prime_max<B2) pm1_init_primes(B2);
  adjust_bufsize((void **)(&B2_scheme),&B2_scheme_alloc,
                 B2_scheme_len+1,BUF_INC,sizeof(pm1_scheme2_t));
  B2_scheme[B2_scheme_len].B1=B1;
  B2_scheme[B2_scheme_len].B2=B2;
  B2_scheme[B2_scheme_len].D=0;
  bestd=0; best=U32_MAX;
  for (d=30; d*d<60*B2; d+=30) {
    nmul=create_B2_scheme(B2_scheme+B2_scheme_len,B1,B2,d,0);
    if (!nmul) continue;
    if (nmul<best) { bestd=d; best=nmul; }
//printf("S %u muls, D: %u (%u * %u)\n",nmul,d,B2_scheme[B2_scheme_len].tablen,B2_scheme[B2_scheme_len].endD-B2_scheme[B2_scheme_len].beginD);
  }
  if (!bestd)
    complain("pm1_init_B2: parameters B1=%u, B2=%u too weird\n",B1,B2);
  d=bestd;
  create_B2_scheme(B2_scheme+B2_scheme_len,B1,B2,d,0);

#ifdef PM1_STAT
/* STAT */
  printf("# STAT: step2 %u muls, D: %u (%u * %u)\n",best,bestd,
         B2_scheme[B2_scheme_len].tablen,
         B2_scheme[B2_scheme_len].endD-B2_scheme[B2_scheme_len].beginD);
#endif

  len=16*(d/60)+16;
  if (len<32) len=32;
  if (mm_B2_tab1_len<len) {
    if (mm_B2_tab1_len) free(mm_B2_tab1);
    mm_B2_tab1_len=len;
    mm_B2_tab1=(ulong **)xmalloc(mm_B2_tab1_len*sizeof(ulong *));
  }
  len=B2_scheme[B2_scheme_len].endD-B2_scheme[B2_scheme_len].beginD;
  len*=2; len+=2;
  if (len<12) len=12;
  if (mm_B2_tab2_len<len) {
    if (mm_B2_tab2_len) free(mm_B2_tab2);
    mm_B2_tab2_len=len;
    mm_B2_tab2=(ulong **)xmalloc(mm_B2_tab2_len*sizeof(ulong *));
  }
  len=B2_scheme[B2_scheme_len].endD-B2_scheme[B2_scheme_len].beginD;
  len+=B2_scheme[B2_scheme_len].tablen;
  B2_scheme_len++;
}
 
/* ------------------------ init ---------------------- */

static void pm1_init_pointers()
{
  u32_t i;
  ulong *ptr;
  size_t mem;

/* check stack space */
  mem=(size_t)(mm_B2_tab1_len+mm_B2_tab2_len);
  adjust_bufsize((void **)(&mm_stack),&mm_stack_alloc,
                 mem*montgomery_ulongs,BUF_INC,sizeof(ulong *));

  ptr=mm_stack;
  for (i=0; i<mm_B2_tab1_len; i++) {
    mm_B2_tab1[i]=ptr; ptr+=montgomery_ulongs;
  }
  for (i=0; i<mm_B2_tab2_len; i++) {
    mm_B2_tab2[i]=ptr; ptr+=montgomery_ulongs;
  }
}

/* ------------------ step 1 --------------------- */

ulong mm_b[NMAX_ULONGS];
ulong mm_one[NMAX_ULONGS];
ulong mm_a[NMAX_ULONGS];

ulong mm_u[NMAX_ULONGS];
ulong mm_w[NMAX_ULONGS];
ulong mm_prod[NMAX_ULONGS];

ulong mm_A[NMAX_ULONGS];
ulong mm_B[NMAX_ULONGS];
ulong mm_C[NMAX_ULONGS];

static int pm1_step1()
{
  u32_t i;

  asm_zero(mm_one); mm_one[0]=1;
  asm_mulmod(mm_one,montgomery_modulo_R2,mm_one); /* one=1 */

  asm_copy(mm_a,mm_one);
  asm_add2(mm_a,mm_a);  /* a=2 */
  for (i=0; i<B1_scheme[pm1_index1].len; i++) {
    asm_squmod(mm_a,mm_a);
    if (get_bit(B1_scheme[pm1_index1].ex,i)) asm_add2(mm_a,mm_a);
  }

  asm_sub(mm_a,mm_a,mm_one);
  if (asm_invert(mm_b,mm_a)==0) {
    gcd(mm_b,mm_a,montgomery_modulo_n);
    mpz_set_v(gmp_f,mm_b);
    return 1;
  }
  asm_add2(mm_a,mm_one);
  return 0;
}

/* ------------------ step 2 --------------------- */

/* uses mm_u */
static void pm1_power(ulong *rop, ulong *op, u32_t ex)
{
  u32_t i;

  if (!ex) { asm_copy(rop,mm_one); return; }
  if (ex==1) { asm_copy(rop,op); return; }
  for (i=31; i; i--) if (ex&(1<<i)) break;
  asm_copy(mm_u,op);
  for (; i; i--) {
    asm_squmod(mm_u,mm_u);
    if (ex&(1<<(i-1))) asm_mulmod(mm_u,mm_u,op);
  }
  asm_copy(rop,mm_u);
}


static int pm1_step2()
{
  u32_t i, j, d, b0, b1, n, p, len1, len2;
  ulong **t1, **t2, **tab1;
  u32_t ind;
  uchar *te;

  d=B2_scheme[pm1_index2].D;
  b0=B2_scheme[pm1_index2].beginD;
  b1=B2_scheme[pm1_index2].endD;
  len1=B2_scheme[pm1_index2].tablen;
  len2=b1-b0;
  tab1=(ulong **)xmalloc(len1*sizeof(ulong *));
  t1=mm_B2_tab1;
  t2=mm_B2_tab2;
/* compute first table: */
  asm_copy(t1[0],mm_a);
  asm_squmod(t2[0],t1[0]); /* 2 */
  asm_mulmod(t2[0],t2[0],t1[0]); /* 3 */
  asm_squmod(t2[0],t2[0]); /* 6 */
  asm_squmod(t2[0],t2[0]); /* 12 */
  asm_squmod(t2[0],t2[0]); /* 24 */
  asm_squmod(t2[1],t2[0]); /* 48 */
  asm_mulmod(t2[2],t2[1],t2[0]); /* 72 */
  asm_mulmod(t2[3],t2[2],t2[1]); /* 120 */

  asm_mulmod(t1[2],t2[1],t1[0]); /* 49 */
  asm_mulmod(t1[6],t2[3],t1[2]); /* 169 */
  asm_mulmod(t2[4],t2[3],t2[2]); /* 192 */
  asm_mulmod(t1[10],t2[4],t1[6]); /* 361 */
  asm_mulmod(t1[4],t2[3],t1[0]); /* 121 */
  asm_mulmod(t2[5],t2[3],t2[1]); /* 168 */
  asm_mulmod(t1[8],t2[5],t1[4]); /* 289 */
  asm_mulmod(t2[6],t2[5],t2[2]); /* 240 */
  asm_mulmod(t1[12],t2[6],t1[8]); /* 529 */
  asm_mulmod(t2[7],t2[6],t2[2]); /* 312 */
  asm_mulmod(t1[14],t2[7],t1[12]); /* 841 */

/* extend to 30+{1,7,11,13,17,19,23,29} */
  if (d>30) {
    asm_mulmod(t2[9],t2[6],t2[3]); /* 360 */
    asm_squmod(t2[0],t2[6]);
    asm_squmod(t2[0],t2[0]);  /* 960 */
    asm_mulmod(t2[1],t2[0],t2[9]); /* 1320 */
    asm_mulmod(t2[2],t2[1],t2[6]); /* 1560 */
    asm_mulmod(t2[8],t2[2],t2[6]); /* 1800 */
    asm_mulmod(t2[3],t2[1],t2[9]); /* 1680 */
    asm_mulmod(t2[5],t2[3],t2[9]); /* 2040 */
    asm_mulmod(t2[4],t2[2],t2[9]); /* 1920 */
    asm_mulmod(t2[6],t2[4],t2[9]); /* 2280 */
    asm_mulmod(t2[7],t2[6],t2[9]); /* 2640 */
    for (i=0; i<8; i++)
      asm_mulmod(t1[16+2*i],t1[2*i],t2[i]);
  }

/* extend to [60,d/2] */
  for (j=2; 60*j<d; j++)
    for (i=0; i<8; i++) {
      asm_mulmod(t2[i],t2[i],t2[8]);
      asm_mulmod(t1[16*j+2*i],t1[16*(j-1)+2*i],t2[i]);
    }

  for (i=0; i<len1; i++) {
    u32_t t[8]={1,7,11,13,17,19,23,29}, m, k;

    n=B2_scheme[pm1_index2].tab[i]; m=n/30; k=n%30;
    for (j=0; j<8; j++) if (t[j]==k) break;
    if (j>=8) Schlendrian("pm1_step2\n");
    tab1[i]=t1[16*m+2*j];
  }

/* compute second table */
  pm1_power(t2[0],mm_a,d*d);
  pm1_power(t2[3],t2[0],b0);
  pm1_power(t2[2],t2[3],b0);  /* b0^2 */
  asm_squmod(t2[3],t2[3]); /* 2*b0 */
  asm_squmod(t2[1],t2[0]);  /* 2 */
  asm_mulmod(t2[0],t2[0],t2[3]); /* 2*b0+1 */
  asm_mulmod(t2[4],t2[2],t2[0]); /* (b0+1)^2 */
  for (i=2; i+b0<b1; i++) {
    asm_mulmod(t2[0],t2[0],t2[1]);
    asm_mulmod(t2[2*i+2],t2[2*i],t2[0]);
  }

/* rpc */
  asm_zero(mm_prod); mm_prod[0]=1;
  ind=1;
  te=B2_scheme[pm1_index2].rpc;

/* todo: implement len1 or len2>254, probably not needed in the near future */
if (len2>254) complain("to implement\n");
if (len1>254) complain("to implement\n");
  if (te[0]) {
    for (i=1; i<=len2; i++) asm_squmod(t2[2*i+1],t2[2*i]);
    while (te[ind]!=255) {
      u32_t ii[3];

      ii[0]=(u32_t)(te[ind++]);
      ii[1]=(u32_t)(te[ind++]);
      ii[2]=(u32_t)(te[ind++]);
      asm_copy(mm_A,tab1[ii[0]]);
      asm_add2(mm_A,tab1[ii[1]]);
      asm_mulmod(mm_u,mm_A,tab1[ii[2]]);
      asm_add2(mm_A,tab1[ii[2]]);
      asm_mulmod(mm_B,tab1[ii[0]],tab1[ii[1]]);
      asm_mulmod(mm_C,mm_B,tab1[ii[2]]);
      asm_add2(mm_B,mm_u);
      asm_mulmod(mm_u,mm_A,mm_B);
      asm_sub(mm_C,mm_C,mm_u);
      while (te[ind]!=255) {
        i=(u32_t)(te[ind++]);
        if (i>len2) Schlendrian("pm1_step2 %u %u %u\n",i,len2,ind);
/* test tab1[ii[0,1,2]] against tab2[i] */
        asm_copy(mm_w,t2[2*i+3]);
        asm_add2(mm_w,mm_B);
        asm_sub(mm_u,t2[2*i+2],mm_A);
        asm_mulmod(mm_u,mm_u,mm_w);
        asm_sub(mm_u,mm_u,mm_C);
        asm_mulmod(mm_prod,mm_prod,mm_u);
      }
      ind++;
    }
  } else {
    for (i=0; i<len1; i++) asm_squmod(t1[2*i+1],tab1[i]);
    while (te[ind]!=255) {
      u32_t ii[3];

      ii[0]=2+2*((u32_t)(te[ind++]));
      ii[1]=2+2*((u32_t)(te[ind++]));
      ii[2]=2+2*((u32_t)(te[ind++]));
      asm_copy(mm_A,t2[ii[0]]);
      asm_add2(mm_A,t2[ii[1]]);
      asm_mulmod(mm_u,mm_A,t2[ii[2]]);
      asm_add2(mm_A,t2[ii[2]]);
      asm_mulmod(mm_B,t2[ii[0]],t2[ii[1]]);
      asm_mulmod(mm_C,mm_B,t2[ii[2]]);
      asm_add2(mm_B,mm_u);
      asm_mulmod(mm_u,mm_A,mm_B);
      asm_sub(mm_C,mm_C,mm_u);
      while (te[ind]!=255) {
        i=(u32_t)(te[ind++]);
        if (i>len1) Schlendrian("pm1_step2\n");
/* test tab2[ii[0,1,2]] against tab1[i] */
        asm_copy(mm_w,t1[2*i+1]);
        asm_add2(mm_w,mm_B);
        asm_sub(mm_u,tab1[i],mm_A);
        asm_mulmod(mm_u,mm_u,mm_w);
        asm_sub(mm_u,mm_u,mm_C);
        asm_mulmod(mm_prod,mm_prod,mm_u);
      }
      ind++;
    }
  }

/* tests */
  te=B2_scheme[pm1_index2].tests;
  ind=0;
  for (i=0; i<len2; i++) {
    j=0;
    while (te[ind]) {
      j+=(u32_t)(te[ind++]);
      if (j>len1) Schlendrian("pm1_step2\n");
/* test tab1[j-1] against tab2[i] */
      asm_sub(mm_w,tab1[j-1],t2[2*i+2]);
      asm_mulmod(mm_prod,mm_prod,mm_w);
    }
    ind++;
  }
  free(tab1);
  if (asm_invert(mm_w,mm_prod)==0) {
    gcd(mm_b,mm_prod,montgomery_modulo_n);
    mpz_set_v(gmp_f,mm_b);
    return 1;
  }

  return 0;
}

/* ------------------ main functions --------------------- */

static int pm1_init(mpz_t N, u32_t B1, u32_t B2)
{
  u32_t i;

  if (!set_montgomery_multiplication(N)) return -1;
  if (!pm1_is_init) {
    mpz_init(gmp_f);
    pm1_is_init=1;
  }
  if (B2<=B1) return -1;
//  if (B2<=B1) B2=2*B1; /* at least one prime between B1 and B2 */
  for (i=0; i<(u32_t)B1_scheme_len; i++)
    if (B1_scheme[i].B1==B1)
      break;
  if (i==(u32_t)B1_scheme_len) pm1_init_B1(B1);
  pm1_index1=i;

  for (i=0; i<(u32_t)B2_scheme_len; i++)
    if ((B2_scheme[i].B1==B1) && (B2_scheme[i].B2==B2))
      break;
  if (i==(u32_t)B2_scheme_len) pm1_init_B2(B1,B2);
  pm1_index2=i;

  pm1_init_pointers();
  return 0;
}


int pm1(mpz_t **fptr)
{
  if (pm1_step1()) {  *fptr=&gmp_f; return 1; }
  if (pm1_step2()) {  *fptr=&gmp_f; return 1; }
  return 0;
}


int pm1_factor(mpz_t N, u32_t B1, u32_t B2, mpz_t **fptr)
{
  int suc;

  if (pm1_init(N,B1,B2)) return -1;
  suc=pm1(fptr);
  return suc;
}

