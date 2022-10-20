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
#include "ecm.h"


#define BUF_INC    256

extern ulong montgomery_inv_n;
extern ulong *montgomery_modulo_n;
extern ulong montgomery_modulo_R2[NMAX_ULONGS];
extern size_t montgomery_ulongs;


extern void gcd(ulong *,ulong *,ulong *);
extern int asm_invert(ulong *,ulong *);



static int ecm_is_init=0;
mpz_t gmp_f, gmp_aux0, gmp_aux1, gmp_aux2, gmp_aux3;

static uchar *ecm_prime_bit;
static u32_t ecm_prime_max=0;

static uchar *addition_chain;
static u32_t *B1_ac, *B1_prime, B1_len, B1_ac_maxlen=0, B1_max=0;


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
} scheme2_t;

static scheme2_t *B2_scheme;
static size_t B2_scheme_len=0, B2_scheme_alloc=0;

static ulong *mm_stack;
static size_t mm_stack_alloc=0;

static ulong **mm_B1_x, **mm_B1_z;
static ulong **mm_B2_inv, **mm_B2_tab1, **mm_B2_tab2;
static u32_t mm_B2_inv_len=0, mm_B2_tab1_len=0, mm_B2_tab2_len=0;


/* ------------------- conversion ------------------- */

void mpz_set_v(mpz_t targ, ulong *v)
{
  size_t i;

  mpz_set_ui(targ,0);
  for (i=0; i<montgomery_ulongs; i++) {
    mpz_mul_2exp(targ,targ,mp_bits_per_limb);
    mpz_add_ui(targ,targ,v[montgomery_ulongs-1-i]);
  }
}


void mpz_get_v(ulong *v, mpz_t src)
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

uchar pop8_tab[256]={
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

u32_t popcount(u64_t *v, u32_t len)
{
  u32_t res, i;
  uchar *u;

  u=(uchar *)v;
  for (i=0,res=0; i<8*len; i++) res+=(u32_t)(pop8_tab[u[i]]);
  return res;
}

/* ------------------------ init primes ------------------------- */

/* ok, if limit is <2^32-2^17 or something like that */
static void ecm_init_primes(u32_t limit)
{
  u32_t i, j, l;

  if (limit>(1UL<<31)) complain("ecm_init_primes: too big\n");
  if (ecm_prime_max) free(ecm_prime_bit);
  l=(limit+15)>>4;
  ecm_prime_max=l<<4;
/* ecm_prime_bit[i] stores odd primes in [16i+1,16i+15] */
  ecm_prime_bit=(uchar *)xmalloc(l*sizeof(uchar));
  memset(ecm_prime_bit,255,l*sizeof(uchar));
  ecm_prime_bit[0]=254; /* 1 is not prime */
  j=3;
  do {
    for (i=j*j; i<ecm_prime_max; i+=2*j) clear_bit(ecm_prime_bit,i>>1);
    while (j<ecm_prime_max) {
      j+=2;
      if (get_bit(ecm_prime_bit,j>>1)) break;
    }
    if (j>65535) j=65535;
  } while (j*j<ecm_prime_max);
}

/* ------------------------ init phase 1 ---------------------- */

static u32_t find_addition_chain(uchar **rop, u32_t p)
{
  u32_t i, j, *seq, r, best, nmul;
  // SMJS Initialise ind here to stop compiler complaint
  u32_t diff, l, ll, ind = 0;
  uchar *tmp;

  seq=(u32_t *)xmalloc(p*sizeof(u32_t)); /* ok for p<10000 */
  tmp=(uchar *)xmalloc((1+4*p)*sizeof(uchar));
  best=U32_MAX;
  for (r=(p+1)/2; r<p; r++) {
    seq[0]=p; seq[1]=r; i=2;
    while (1) {
      l=seq[i-1];
      if (l<2) break;
      if (l==2) { seq[i]=1; i++; break; }
      ll=seq[i-2];
      diff=ll-l;
      if (diff&1) {
        if ((l&1) && (3*diff<l)) {
          seq[i]=ll/2; seq[i+1]=ll/2-diff;
          if (diff>seq[i+1]) seq[i+1]=diff;
          i+=2;
          continue;
        }
        if ((ll&1) && (4*diff<l)) {
          seq[i]=l/2+diff; seq[i+1]=l/2;
          if (diff>seq[i+1]) seq[i+1]=diff;
          i+=2;
          continue;
        }
      } else {
        if ((((l+ll)&3)==0) && (17*diff<2*l)) {
          seq[i]=diff; seq[i+1]=(3*ll-l)/4;
          seq[i+2]=(ll+l)/4; seq[i+3]=(3*l-ll)/4;
          if (diff>seq[i+3]) seq[i+3]=diff;
          i+=4;
          continue;
        }
/* todo: cases 4,8,16.. divides diff */
      }
      if (2*diff<l) seq[i]=l-diff; else seq[i]=diff;
      i++;
    }
    nmul=0;
    for (j=0; j<i-1; j++)
      if (seq[j]&1) nmul+=6;
      else {
        for (l=j+1; l<i; l++) if (2*seq[l]==seq[j]) break;
        if (l<i) nmul+=5; else nmul+=6;
      }
//for (j=0; j<i; j++) printf("%u ",seq[j]); printf("   %u\n",nmul);
    if (nmul<best) {
      best=nmul;
      ind=0; j=i-1;
      while (j) {
        j--;
        ll=seq[j];
        if ((ll&1)==0) { /* search for duplication */
          for (l=j+1; l<i; l++) if (2*seq[l]==ll) break;
          if (l<i) {
            tmp[ind++]=2; tmp[ind++]=(uchar)(i-1-l);
            continue;
          }
        }
/* search for addition */
        diff=0;
        for (l=j+1; l<i; l++) {
          u32_t i0, i1;
          for (i0=j+1; i0<i; i0++) {
            if (seq[l]+seq[i0]==ll) {
              for (i1=j+1; i1<i; i1++)
                if ((seq[l]==seq[i0]+seq[i1]) || (seq[i0]==seq[l]+seq[i1])) {
                  tmp[ind++]=4; tmp[ind++]=(uchar)(i-1-l);
                  tmp[ind++]=(uchar)(i-1-i0); tmp[ind++]=(uchar)(i-1-i1);
                  diff=1;
                  break;
                }
            }
            if (diff) break;
            if (seq[l]==seq[i0]+ll) {
              for (i1=j+1; i1<i; i1++)
                if (seq[l]+seq[i0]==seq[i1]) {
                  tmp[ind++]=4; tmp[ind++]=(uchar)(i-1-l);
                  tmp[ind++]=(uchar)(i-1-i0); tmp[ind++]=(uchar)(i-1-i1);
                  diff=1;
                  break;
                }
            }
            if (diff) break;
          }
          if (diff) break;
        }
        if (diff) continue;
        Schlendrian("find_addition_chain for p=%u\n",p);
      }
      tmp[ind++]=0;
    }
  }
//  printf("%u: %u  ",p,best);
  free(seq);
  *rop=tmp;
  return ind;
}


static u32_t ac_mul_count(uchar *rop)
{
  u32_t res=0, l;
  uchar *ptr=rop;

  while (*ptr) {
    l=(u32_t)(*ptr);
    res+=(4+l/2);
    ptr+=l;
  }
  return res;
}


static u32_t ac_add_count(uchar *rop)
{
  u32_t res=0, l;
  uchar *ptr=rop;

  while (*ptr) {
    l=(u32_t)(*ptr);
    res+=(2+l);
    ptr+=l;
  }
  return res;
}


static u32_t ac_copy_count(uchar *rop)
{
  u32_t res=0, l;
  uchar *ptr=rop;

  while (*ptr) {
    l=(u32_t)(*ptr);
    res+=(3*(l/2)-2);
    ptr+=l;
  }
  return res;
}


static void ecm_init_B1(u32_t B1, u32_t B2)
{
  u32_t i, ind, n, len, m;
  uchar **ac_tmp;

  if (B1_max) {
    free(addition_chain);
    free(B1_ac);
    free(B1_prime);
    free(mm_B1_x);
    free(mm_B1_z);
  }
  if (ecm_prime_max<B1) ecm_init_primes(B1);
  for (i=1,n=0; 2*i+1<B1; i++) if (get_bit(ecm_prime_bit,i)) n++;
  B1_prime=(u32_t *)xmalloc(n*sizeof(*B1_prime));
  for (i=1,ind=0; 2*i+1<B1; i++)
    if (get_bit(ecm_prime_bit,i))
      B1_prime[ind++]=2*i+1;
  if (ind!=n) Schlendrian("ecm_init_B1\n");
  B1_ac=(u32_t *)xmalloc(n*sizeof(*B1_prime));
  B1_len=n;
/* search for good addition chains */
  ac_tmp=(uchar **)xmalloc(B1_len*sizeof(uchar *));
  for (i=0,n=0; i<B1_len; i++)
    n+=find_addition_chain(ac_tmp+i,B1_prime[i]);

  addition_chain=(uchar *)xmalloc(n*sizeof(uchar));
  for (i=0,ind=0; i<B1_len; i++) {
    B1_ac[i]=ind;
    len=0; m=1;
    while (ac_tmp[i][len]) { len+=((u32_t)ac_tmp[i][len]); m++; }
    len++;
    memcpy(addition_chain+ind,ac_tmp[i],len*sizeof(uchar));
    ind+=len;
    if (B1_ac_maxlen<m) B1_ac_maxlen=m;
  }
  if (n!=ind) Schlendrian("ecm_init_B1 %u %u\n",n,ind);

  for (i=0; i<B1_len; i++) free(ac_tmp[i]);
  free(ac_tmp);
  mm_B1_x=(ulong **)xmalloc(B1_ac_maxlen*sizeof(ulong *));
  mm_B1_z=(ulong **)xmalloc(B1_ac_maxlen*sizeof(ulong *));
  B1_max=B1;

#ifdef ECM_STAT	
{
  u32_t na, nc, ma, mc;

/* STAT */
  n=0; na=0; nc=0;
  for (i=2; i<B2; i*=2) { n+=5; na+=4; nc+=1; }
  for (i=3; i<B2; i*=3) { n+=11; na+=10; nc+=5; }
  for (ind=1; ind<B1_len; ind++) {
    len=B1_ac[ind];
    m=ac_mul_count(addition_chain+len);
    ma=ac_add_count(addition_chain+len);
    mc=ac_copy_count(addition_chain+len);
    for (i=B1_prime[ind]; i<B1; i*=B1_prime[ind]) { n+=m; na+=ma; nc+=mc; }
  }
  printf("# STAT: step1 %u muls  %u adds  %u copys\n",n,na,nc);
}
#endif
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


static u32_t create_B2_scheme(scheme2_t *s, u32_t B1, u32_t B2, u32_t d, int verb)
{
  u32_t nmul, n, i, b0, b1, p;
  u32_t dim1, dim2, i1, i2, l, *tabnr;
  u32_t n0, n1;
  uchar **pr, **prt, *tmp, *tmpt;

  if ((!d) || (d%30)) Schlendrian("create_B2_scheme\n");
  if (s->D) { free(s->rpc); free(s->tests); free(s->tab); }
  s->D=d;
  nmul=80; /* table up to 30 and 30Q */
  nmul+=24*(d/30-1); /* table up to D */
/* 30*Q -> D*Q: */
  n=d/30; while (!(n&1)) { n>>=1; nmul+=5; }
  for (i=0; i<B1_len; i++) {
    p=B1_prime[i];
    if (n%p==0) {
      n/=p;
      nmul+=ac_mul_count(addition_chain+B1_ac[i]);
      while (n%p==0) {
        n/=p;
        nmul+=ac_mul_count(addition_chain+B1_ac[i]);
      }
      if (n==1) break;
    }
  }
/* avoid cases where D/30 is not B1-smooth */
  if (n>1) { s->D=0; return 0; }

  if (2*B1<d) b0=0; else b0=(B1-d/2)/d;
  n=b0; while (n>2) { n>>=1; nmul+=11; }
  if (n) nmul+=5; if (n>1) nmul+=6; /* computing b0*D*Q, (b0+1)*D*Q */
  b1=1+(B2+d/2)/d;
  nmul+=6*(b1-b0-2); if (!b0) nmul--; /* table of multiples of D */
  s->beginD=b0;
  s->endD=b1;

  dim1=b1-b0; dim2=d/4;
  pr=(uchar **)xmalloc(dim1*sizeof(uchar *));
  for (i=0; i<dim1; i++) pr[i]=(uchar *)xcalloc((size_t)dim2,sizeof(uchar));
  for (i=B1/2; 2*i+1<B2; i++)
    if (get_bit(ecm_prime_bit,i)) {
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
  nmul+=3*(n+b1-b0); /* Inversion */

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

#ifdef ECM_STAT
static void ecm_step2_count(u32_t *naptr, u32_t *ncptr)
{
  u32_t i, j, d, b0, b1, n, p, len1, len2;
  u32_t ind;
  uchar *te;
  scheme2_t sch;
  u32_t na, nc, nadd, ndup;

  na=0; nc=0; nadd=0; ndup=0;
  sch=B2_scheme[B2_scheme_len];
  d=sch.D;
  b0=sch.beginD;
  b1=sch.endD;
  len1=sch.tablen;
  len2=b1-b0;
/* Q=(mm_x,mm_z), compute table: */
  ndup++; /* 2Q */
  ndup++; /* 4Q */
  nadd++; /* 3Q */
  ndup++; /* 6Q */
  ndup++; /* 12Q */
  nadd++; /* 18Q */

  nc+=2; /* Q */
  nadd++; /* 7Q */
  nadd++; /* 13Q */
  nadd++; /* 19Q */
  nadd++; /* 11Q */
  nadd++; /* 17Q */
  nadd++; /* 23Q */
  nadd++; /* 29Q */
  nadd++; /* 30Q */

/* extend to 30+{1,7,11,13,17,19,23,29} */
  if (d>30)
    for (i=0; i<8; i++)
      nadd++;
/* extend to [60,d/2] */
  for (j=2; 60*j<d; j++)
    for (i=0; i<8; i++)
      nadd++;

/* compute D*Q */
  n=d/30;
  while (!(n&1)) {
    n>>=1;
    ndup++;
  }
  for (i=0; i<B1_len; i++) {
    p=B1_prime[i];
    if (n%p==0) {
      n/=p;
      na+=ac_add_count(addition_chain+B1_ac[i]);
      while (n%p==0) {
        n/=p;
        na+=ac_add_count(addition_chain+B1_ac[i]);
      }
      if (n==1) break;
    }
  }
  if (n>1) complain("ecm_step2: unusual B1,B2-parameters\n");

/* compute beginD*D*Q, (beginD+1)*D*Q */
  if (!b0) {
    nc+=2; /* 0Q */
    nc+=2;   /* 1Q */
  } else {
    n=b0; i=0; while (n>2) { n>>=1; i++; }
    if (n==1) {
      nc+=2;   /* 1Q */
      ndup++;  /* 2Q */
    } else {
      ndup++; /* 2Q */
      nadd++; /* 3Q */
    }
    for (; i; i--) {
      if ((b0>>(i-1))&1) {
        nadd++;
        ndup++;
      } else {
        nadd++;
        ndup++;
      }
    }
  }
/* compute second table */
  if (!b0) { /* we cannot compute 2=1+1 by addition formula */
    ndup++;
    i=3;
  } else i=2;
  for (; i+b0<b1; i++)
    nadd++;

  nc++;

/* rpc */
  nc++;
  ind=1;
  te=sch.rpc;

/* todo: implement len1 or len2>254, probably not needed in the near future */
if (len2>254) complain("to implement\n");
if (len1>254) complain("to implement\n");
  if (te[0]) {
    while (te[ind]!=255) {
      u32_t ii[3];

      ii[0]=(u32_t)(te[ind++]);
      ii[1]=(u32_t)(te[ind++]);
      ii[2]=(u32_t)(te[ind++]);
      nc++;
      na++;
      na++;
      na++;
      na++;
      while (te[ind]!=255) {
        i=(u32_t)(te[ind++]);
        if (i>len2) Schlendrian("ecm_step2 %u %u %u\n",i,len2,ind);
/* test tab1[ii[0,1,2]] against tab2[i] */
        nc++;
        na++;
        na++;
        na++;
      }
      ind++;
    }
  } else {
    while (te[ind]!=255) {
      u32_t ii[3];

      ii[0]=2+2*((u32_t)(te[ind++]));
      ii[1]=2+2*((u32_t)(te[ind++]));
      ii[2]=2+2*((u32_t)(te[ind++]));
      nc++;
      na++;
      na++;
      na++;
      na++;
      while (te[ind]!=255) {
        i=(u32_t)(te[ind++]);
        if (i>len1) Schlendrian("ecm_step2\n");
/* test tab2[ii[0,1,2]] against tab1[i] */
        nc++;
        na++;
        na++;
        na++;
      }
      ind++;
    }
  }

/* tests */
  te=sch.tests;
  ind=0;
  for (i=0; i<len2; i++) {
    j=0;
    while (te[ind]) {
      j+=(u32_t)(te[ind++]);
      if (j>len1) Schlendrian("ecm_step2\n");
/* test tab1[j-1] against tab2[i] */
      na++;
    }
    ind++;
  }

  na+=ndup*4; nc+=ndup;
  na+=nadd*6; nc+=nadd*4;
  *naptr=na;
  *ncptr=nc;
}
#endif


static void ecm_init_B2(u32_t B1, u32_t B2)
{
  u32_t d, nmul, best, bestd, len;

  if (ecm_prime_max<B2) ecm_init_primes(B2);
  adjust_bufsize((void **)(&B2_scheme),&B2_scheme_alloc,
                 B2_scheme_len+1,BUF_INC,sizeof(scheme2_t));
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
    complain("ecm_init_B2: parameters B1=%u, B2=%u too weird\n",B1,B2);
  d=bestd;
  create_B2_scheme(B2_scheme+B2_scheme_len,B1,B2,d,0);

#ifdef ECM_STAT
{
  u32_t na, nc;

/* STAT */
  ecm_step2_count(&na,&nc);
  printf("# STAT: step2 %u muls  %u adds  %u copys, D: %u (%u * %u)\n",
         best,na,nc,bestd,
         B2_scheme[B2_scheme_len].tablen,
         B2_scheme[B2_scheme_len].endD-B2_scheme[B2_scheme_len].beginD);
}
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
  if (mm_B2_inv_len<len) {
    if (mm_B2_inv_len) free(mm_B2_inv);
    mm_B2_inv_len=len;
    mm_B2_inv=(ulong **)xmalloc(mm_B2_inv_len*sizeof(ulong *));
  }
  B2_scheme_len++;
}

/* ------------------------ init ---------------------- */

static void ecm_init()
{
  mpz_init(gmp_f);
  mpz_init(gmp_aux0);
  mpz_init(gmp_aux1);
  mpz_init(gmp_aux2);
  mpz_init(gmp_aux3);

  ecm_is_init=1;
}


static void ecm_init_pointers(ecm_t e)
{
  u32_t i;
  ulong *ptr;
  size_t mem;

/* check stack space */
  mem=(size_t)(2*B1_ac_maxlen+mm_B2_tab1_len+mm_B2_tab2_len+mm_B2_inv_len);
  adjust_bufsize((void **)(&mm_stack),&mm_stack_alloc,
                 mem*montgomery_ulongs,BUF_INC,sizeof(ulong *));

  ptr=mm_stack;
  for (i=0; i<B1_ac_maxlen; i++) {
    mm_B1_x[i]=ptr; ptr+=montgomery_ulongs;
    mm_B1_z[i]=ptr; ptr+=montgomery_ulongs;
  }
  for (i=0; i<mm_B2_tab1_len; i++) {
    mm_B2_tab1[i]=ptr; ptr+=montgomery_ulongs;
  }
  for (i=0; i<mm_B2_tab2_len; i++) {
    mm_B2_tab2[i]=ptr; ptr+=montgomery_ulongs;
  }
  for (i=0; i<mm_B2_inv_len; i++) {
    mm_B2_inv[i]=ptr; ptr+=montgomery_ulongs;
  }
}

/* ------------------- choose elliptic curve ---------------------- */

/*
Choose curves as in Montgomery: Speeding the Pollard and EC Methods of fact.,
Math. Comp. p. 263:
E: y^2=x^3-12x, P=[-2,-4], aP=[x_a,y_a] for a=2,3,...
store x_a mod N in curve_param
*/
static int ecm_next_curve(ecm_t e)
{
#ifdef ECM_ZEIT
zeitA(10);
#endif

  if (!mpz_sgn(e->curve_x)) { /* first curve, 2P=[4,4] */
    mpz_set_ui(e->curve_x,4);
    mpz_set_ui(e->curve_y,4);
    asm_zero(e->curve_param);
    e->curve_param[0]=4;
#ifdef ECM_ZEIT
zeitB(10);
#endif
    return 0;
  }
/* add P to [e->curve_x,e->curve_y] */
  mpz_add_ui(gmp_aux0,e->curve_x,2);
  mpz_add_ui(gmp_aux1,e->curve_y,4);
  if (!mpz_invert(gmp_aux3,gmp_aux0,e->N)) { /* extremely rare */
    mpz_gcd(gmp_f,gmp_aux0,e->N);
#ifdef ECM_ZEIT
zeitB(10);
#endif
    return 1;
  }
  mpz_sub(gmp_aux2,gmp_aux1,gmp_aux0);
  mpz_sub(gmp_aux2,gmp_aux2,gmp_aux0);
  mpz_mul(gmp_aux2,gmp_aux2,gmp_aux3);
  mpz_mul_2exp(gmp_aux2,gmp_aux2,1);   /* nu=2*(-2*x+y)/(x+2) */
  mpz_mul(gmp_aux3,gmp_aux3,gmp_aux1);  /* lambda=(y+4)/(x+2) */
  mpz_sub_ui(gmp_aux0,gmp_aux0,4);
  mpz_neg(gmp_aux0,gmp_aux0);
  mpz_addmul(gmp_aux0,gmp_aux3,gmp_aux3); /* rx=lambda^2-x+2 */
  mpz_mul(gmp_aux1,gmp_aux0,gmp_aux3);
  mpz_neg(gmp_aux1,gmp_aux1);
  mpz_sub(gmp_aux1,gmp_aux1,gmp_aux2);  /* ry=-lambda*rx-nu */
  mpz_mod(e->curve_x,gmp_aux0,e->N);
  mpz_mod(e->curve_y,gmp_aux1,e->N);
  mpz_get_v(e->curve_param,e->curve_x);
#ifdef ECM_ZEIT
zeitB(10);
#endif
  return 0;
}

/* ------------------ elliptic curve operations --------------------- */

ulong mm_b[NMAX_ULONGS];
ulong mm_one[NMAX_ULONGS];
ulong mm_u[NMAX_ULONGS];
ulong mm_v[NMAX_ULONGS];
ulong mm_w[NMAX_ULONGS];
ulong mm_prod[NMAX_ULONGS];

ulong mm_a[NMAX_ULONGS];
ulong mm_x[NMAX_ULONGS];
ulong mm_z[NMAX_ULONGS];
  ulong mm_x1[NMAX_ULONGS];
  ulong mm_z1[NMAX_ULONGS];

ulong mm_A[NMAX_ULONGS];
ulong mm_B[NMAX_ULONGS];
ulong mm_C[NMAX_ULONGS];

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 1 copy
     Uses the following global variables:
     - n : number to factor
     - b : (a+2)/4 mod n
     - u, v, w : auxiliary variables
Modifies: x2, z2, u, v, w
*/
static void ecm_duplicate(ulong *x2, ulong *z2, ulong *x1, ulong *z1)
{
  asm_copy(mm_w,x1); asm_add2(mm_w,z1);
  asm_squmod(mm_u,mm_w); /* u=(x1+z1)^2 */
  asm_sub(mm_w,x1,z1);
  asm_squmod(mm_v,mm_w); /* v=(x1-z1)^2 */
  asm_mulmod(x2,mm_u,mm_v); /* x2=(u*v) */
  asm_sub(mm_w,mm_u,mm_v); /* w=u-v=4*x1*z1 */
  asm_mulmod(mm_u,mm_b,mm_w);
  asm_add2(mm_u,mm_v); /* u=(v+b*w) */
  asm_mulmod(z2,mm_w,mm_u); /* z2=(w*u) */
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
     using 6 mul, 6 add/sub, 4 copy.
     One assumes that Q-R=P or R-Q=P where P=(x:z).
     Uses the following global variables:
     - n : number to factor
     - x, z : coordinates of P
     - u, v, w : auxiliary variables
Modifies: x3, z3, u, v, w.
(x3,z3) may be identical to (x2,z2), (x1,z1) and to (x,z)
*/
/* returns garbage if two points are equal */
static void ecm_add(ulong *x3, ulong *z3, ulong *x2, ulong *z2, ulong *x1, ulong *z1, ulong *x, ulong *z)
{
  asm_sub(mm_u,x2,z2);
  asm_copy(mm_v,x1); asm_add2(mm_v,z1);
  asm_mulmod(mm_u,mm_u,mm_v); /* u=(x2-z2)(x1+z1) */
  asm_copy(mm_w,x2); asm_add2(mm_w,z2);
  asm_sub(mm_v,x1,z1);
  asm_mulmod(mm_v,mm_w,mm_v); /* v=(x2+z2)(x1-z1) */
  asm_copy(mm_w,mm_u); asm_add2(mm_w,mm_v);
  asm_sub(mm_u,mm_u,mm_v);
  asm_squmod(mm_w,mm_w); /* w=(u+v)^2 */
  asm_mulmod(mm_w,mm_w,z);
  asm_squmod(mm_u,mm_u); /* v=(u-v)^2 */
  asm_mulmod(z3,x,mm_u); /* z3=(xw) */
  asm_copy(x3,mm_w); /* x3=(zx3) */
}


static void ecm_mul(ulong *ropx, ulong *ropz, ulong *opx, ulong *opz, uchar *cmd)
{
  u32_t i, j, k, len;
  ulong *rx, *rz, *o1x, *o1z, *o2x, *o2z, *o3x, *o3z;

#if 1
  if (!cmd[0]) Schlendrian("ecm_mul: nothing to do\n");
  asm_copy(mm_B1_x[0],opx);
  asm_copy(mm_B1_z[0],opz);
  i=0; j=1;
  while (cmd[i]) {
    len=(u32_t)(cmd[i]);
    rx=mm_B1_x[j]; rz=mm_B1_z[j];
    if (len==2) {
      k=(u32_t)(cmd[i+1]);
      o1x=mm_B1_x[k]; o1z=mm_B1_z[k];
      ecm_duplicate(rx,rz,o1x,o1z);
    } else if (len==4) {
      k=(u32_t)(cmd[i+1]);
      o1x=mm_B1_x[k]; o1z=mm_B1_z[k];
      k=(u32_t)(cmd[i+2]);
      o2x=mm_B1_x[k]; o2z=mm_B1_z[k];
      k=(u32_t)(cmd[i+3]);
      o3x=mm_B1_x[k]; o3z=mm_B1_z[k];
      ecm_add(rx,rz,o1x,o1z,o2x,o2z,o3x,o3z);
    } else Schlendrian("ecm_mul\n");
    i+=len; j++;
  }
  asm_copy(ropx,mm_B1_x[j-1]);
  asm_copy(ropz,mm_B1_z[j-1]);
#else
  if (!cmd[0]) Schlendrian("ecm_mul: nothing to do\n");
  i=0; j=0;
  while (cmd[i]) {
    len=(u32_t)(cmd[i]);
    if (cmd[i+len]) { rx=mm_B1_x[j]; rz=mm_B1_z[j]; }
    else { rx=ropx; rz=ropz; }
    if (len==2) {
      k=(u32_t)(cmd[i+1]);
      if (k) { o1x=mm_B1_x[k-1]; o1z=mm_B1_z[k-1]; }
      else { o1x=opx; o1z=opz; }
      ecm_duplicate(rx,rz,o1x,o1z);
    } else if (len==4) {
      k=(u32_t)(cmd[i+1]);
      if (k) { o1x=mm_B1_x[k-1]; o1z=mm_B1_z[k-1]; }
      else { o1x=opx; o1z=opz; }
      k=(u32_t)(cmd[i+2]);
      if (k) { o2x=mm_B1_x[k-1]; o2z=mm_B1_z[k-1]; }
      else { o2x=opx; o2z=opz; }
      k=(u32_t)(cmd[i+3]);
      if (k) { o3x=mm_B1_x[k-1]; o3z=mm_B1_z[k-1]; }
      else { o3x=opx; o3z=opz; }
      ecm_add(rx,rz,o1x,o1z,o2x,o2z,o3x,o3z);
    } else Schlendrian("ecm_mul\n");
    i+=len; j++;
  }
#endif
}

/* ------------------ step 1 --------------------- */

static int ecm_step1(ecm_t e)
{
  u32_t i, q, p;

#ifdef ECM_ZEIT
zeitA(3);
#endif

 /* use Montgomerys initialisation */
  asm_zero(mm_one); mm_one[0]=1;
  asm_mulmod(mm_one,montgomery_modulo_R2,mm_one); /* one=1 */

  asm_mulmod(mm_a,montgomery_modulo_R2,e->curve_param);

  asm_copy(mm_v,mm_one); asm_add2(mm_v,mm_one); asm_add2(mm_v,mm_one); /* 3 */
  asm_add2(mm_v,mm_v); asm_add2(mm_v,mm_v); /* 12 */
  asm_squmod(mm_u,mm_a); asm_sub(mm_u,mm_u,mm_v); /* a^2-12 */
  asm_copy(mm_v,mm_u);
  asm_copy(mm_w,mm_a); asm_add2(mm_w,mm_w); asm_add2(mm_w,mm_w); /* 4a */
  asm_sub(mm_u,mm_u,mm_w); /* u=a^2-4a-12 */
  asm_add2(mm_v,mm_w); asm_add2(mm_v,mm_w); asm_add2(mm_v,mm_w); /* v=a^2+12a-12 */

  asm_squmod(mm_x,mm_v); asm_squmod(mm_w,mm_u);
  asm_add2(mm_x,mm_w); asm_add2(mm_x,mm_w); asm_add2(mm_x,mm_w);/* x=3u^2+v^2 */
  asm_mulmod(mm_z,mm_u,mm_v);
  asm_add2(mm_z,mm_z); asm_add2(mm_z,mm_z); /* z=4uv */

  asm_mulmod(mm_b,mm_z,mm_w);  /* b=4u^3v */

  asm_sub(mm_a,mm_v,mm_u); asm_squmod(mm_w,mm_a);
  asm_mulmod(mm_w,mm_w,mm_a);  /* (u-v)^3 */

  asm_copy(mm_a,mm_u);
  asm_add2(mm_a,mm_a); asm_add2(mm_a,mm_u); asm_add2(mm_a,mm_v); /* 3u+v */

  asm_mulmod(mm_a,mm_w,mm_a); /* a=(v-u)^3*(3u+v) */
  if (asm_invert(mm_u,mm_b)==0) {
    gcd(mm_a,mm_b,montgomery_modulo_n);
    mpz_set_v(gmp_f,mm_a);
    return 1;
  }
  asm_mulmod(mm_a,mm_a,mm_u);
  asm_sub(mm_a,mm_a,mm_one); asm_sub(mm_a,mm_a,mm_one); /* a=(v-u)^3*(3u+v)/(4u^3v)-2 */
  asm_copy(mm_b,mm_a);
  asm_half(mm_b);
  asm_add2(mm_b,mm_one);
  asm_half(mm_b); /* now b=(a+2)/4 */

#ifdef ECM_ZEIT
zeitB(3);
#endif

/* use B2 as upper limit for powers of 2 and 3 to take into account that
   there is a torsion group of order divisible by 12 or 24 */
  for (q=2; q<B2_scheme[e->B2_index].B2; q*=2)
    ecm_duplicate(mm_x,mm_z,mm_x,mm_z);
  for (q=3; q<B2_scheme[e->B2_index].B2; q*=3) {
    ecm_duplicate(mm_x1,mm_z1,mm_x,mm_z);
    ecm_add(mm_x,mm_z,mm_x,mm_z,mm_x1,mm_z1,mm_x,mm_z);
  }
  for (i=1; i<B1_len; i++) {
    p=B1_prime[i]; if (p>=e->B1) break;
    for (q=p; q<e->B1; q*=p)  /* upper limit B1 for other primes */
      ecm_mul(mm_x,mm_z,mm_x,mm_z,addition_chain+B1_ac[i]);
  }

#if 0
/* if inversion is very very fast: */
  if (asm_invert(mm_w,mm_z)==0) {
    gcd(mm_a,mm_z,montgomery_modulo_n);
    mpz_set_v(gmp_f,mm_a);
    return 1;
  }
#endif

  return 0;
}

/* ------------------ step 2 --------------------- */


static int ecm_step2(ecm_t e)
{
  u32_t i, j, d, b0, b1, n, p, len1, len2;
  ulong **t1, **t2, **invptr, **invarray, **tab1;
  u32_t ind;
  uchar *te;

  d=B2_scheme[e->B2_index].D;
  b0=B2_scheme[e->B2_index].beginD;
  b1=B2_scheme[e->B2_index].endD;
  len1=B2_scheme[e->B2_index].tablen;
  len2=b1-b0;
  invptr=(ulong **)xmalloc((len1+len2)*sizeof(ulong *));
  tab1=(ulong **)xmalloc(len1*sizeof(ulong *));
  invarray=mm_B2_inv;
  t1=mm_B2_tab1;
  t2=mm_B2_tab2;
/* Q=(mm_x,mm_z), compute table: */
  ecm_duplicate(t2[0],t2[1],mm_x,mm_z); /* 2Q */
  ecm_duplicate(t2[2],t2[3],t2[0],t2[1]); /* 4Q */
  ecm_add(t2[4],t2[5],t2[0],t2[1],mm_x,mm_z,mm_x,mm_z); /* 3Q */
  ecm_duplicate(t2[6],t2[7],t2[4],t2[5]); /* 6Q */
  ecm_duplicate(t2[8],t2[9],t2[6],t2[7]); /* 12Q */
  ecm_add(t2[10],t2[11],t2[8],t2[9],t2[6],t2[7],t2[6],t2[7]); /* 18Q */

  asm_copy(t1[0],mm_x); asm_copy(t1[1],mm_z); /* Q */
  ecm_add(t1[2],t1[3],t2[2],t2[3],t2[4],t2[5],mm_x,mm_z); /* 7Q */
  ecm_add(t1[6],t1[7],t1[2],t1[3],t2[6],t2[7],mm_x,mm_z); /* 13Q */
  ecm_add(t1[10],t1[11],t1[6],t1[7],t2[6],t2[7],t1[2],t1[3]); /* 19Q */
  ecm_add(t1[4],t1[5],t2[8],t2[9],mm_x,mm_z,t1[6],t1[7]); /* 11Q */
  ecm_add(t1[8],t1[9],t2[10],t2[11],mm_x,mm_z,t1[10],t1[11]); /* 17Q */
  ecm_add(t1[12],t1[13],t1[8],t1[9],t2[6],t2[7],t1[4],t1[5]); /* 23Q */
  ecm_add(t1[14],t1[15],t1[12],t1[13],t2[6],t2[7],t1[8],t1[9]); /* 29Q */

  ecm_add(t2[0],t2[1],t2[10],t2[11],t2[8],t2[9],t2[6],t2[7]); /* 30Q */

/* extend to 30+{1,7,11,13,17,19,23,29} */
  if (d>30)
    for (i=0; i<8; i++)
      ecm_add(t1[16+2*i],t1[17+2*i],t2[0],t2[1],t1[2*i],t1[1+2*i],
              t1[14-2*i],t1[15-2*i]);
/* extend to [60,d/2] */
  for (j=2; 60*j<d; j++)
    for (i=0; i<8; i++)
      ecm_add(t1[16*j+2*i],t1[16*j+2*i+1],t2[0],t2[1],
              t1[16*(j-1)+2*i],t1[16*(j-1)+2*i+1],
              t1[16*(j-2)+2*i],t1[16*(j-2)+2*i+1]);

  for (i=0; i<len1; i++) {
    u32_t t[8]={1,7,11,13,17,19,23,29}, m, k;

    n=B2_scheme[e->B2_index].tab[i]; m=n/30; k=n%30;
    for (j=0; j<8; j++) if (t[j]==k) break;
    if (j>=8) Schlendrian("ecm_step2\n");
    tab1[i]=t1[16*m+2*j];
    invptr[i]=t1[16*m+2*j+1];
  }

/* compute D*Q */
  n=d/30;
  while (!(n&1)) {
    n>>=1;
    ecm_duplicate(t2[0],t2[1],t2[0],t2[1]);
  }
  for (i=0; i<B1_len; i++) {
    p=B1_prime[i];
    if (n%p==0) {
      n/=p;
      ecm_mul(t2[0],t2[1],t2[0],t2[1],addition_chain+B1_ac[i]);
      while (n%p==0) {
        n/=p;
        ecm_mul(t2[0],t2[1],t2[0],t2[1],addition_chain+B1_ac[i]);
      }
      if (n==1) break;
    }
  }
  if (n>1) complain("ecm_step2: unusual B1,B2-parameters\n");

/* compute beginD*D*Q, (beginD+1)*D*Q */
  if (!b0) {
    asm_copy(t2[2],mm_one); asm_copy(t2[3],mm_one); /* 0Q */
    asm_copy(t2[4],t2[0]); asm_copy(t2[5],t2[1]);   /* 1Q */
  } else {
    n=b0; i=0; while (n>2) { n>>=1; i++; }
    if (n==1) {
      asm_copy(t2[2],t2[0]); asm_copy(t2[3],t2[1]);   /* 1Q */
      ecm_duplicate(t2[4],t2[5],t2[0],t2[1]);         /* 2Q */
    } else {
      ecm_duplicate(t2[2],t2[3],t2[0],t2[1]);         /* 2Q */
      ecm_add(t2[4],t2[5],t2[2],t2[3],t2[0],t2[1],t2[0],t2[1]); /* 3Q */
    }
    for (; i; i--) {
      if ((b0>>(i-1))&1) {
        ecm_add(t2[2],t2[3],t2[4],t2[5],t2[2],t2[3],t2[0],t2[1]);
        ecm_duplicate(t2[4],t2[5],t2[4],t2[5]);
      } else {
        ecm_add(t2[4],t2[5],t2[4],t2[5],t2[2],t2[3],t2[0],t2[1]);
        ecm_duplicate(t2[2],t2[3],t2[2],t2[3]);
      }
    }
  }
/* compute second table */
  if (!b0) { /* we cannot compute 2=1+1 by addition formula */
    ecm_duplicate(t2[6],t2[7],t2[4],t2[5]);
    i=3;
  } else i=2;
  for (; i+b0<b1; i++)
    ecm_add(t2[2*i+2],t2[2*i+3],t2[2*i],t2[2*i+1],t2[0],t2[1],t2[2*i-2],t2[2*i-1]);

  for (i=0; i+b0<b1; i++) {
    invptr[len1+i]=t2[2*i+3];
  }

/* normalise */
  asm_mulmod(invarray[0],invptr[0],invptr[1]);
  for (i=1; i<len1+len2-1; i++) asm_mulmod(invarray[i],invarray[i-1],invptr[i+1]);
  if (asm_invert(invarray[len1+len2-1],invarray[len1+len2-2])==0) {
    gcd(mm_a,invarray[len1+len2-2],montgomery_modulo_n);
    mpz_set_v(gmp_f,mm_a);
    free(tab1);
    free(invptr);
    return 1;
  }
  for (i=len1+len2-1; i>1; i--) {
    asm_mulmod(invarray[i-1],invarray[i],invptr[i]);
    asm_mulmod(invptr[i],invarray[i],invarray[i-2]);
  }
  asm_mulmod(invarray[0],invarray[1],invptr[1]);
  asm_mulmod(invptr[1],invptr[0],invarray[1]);
  asm_copy(invptr[0],invarray[0]);

  for (i=0; i<len1; i++) asm_mulmod(tab1[i],tab1[i],invptr[i]);
  for (i=0; i<len2; i++) asm_mulmod(t2[2*i+2],t2[2*i+2],invptr[len1+i]);

/* rpc */
//  asm_zero(mm_prod); mm_prod[0]=1;
  asm_copy(mm_prod,mm_z);
  ind=1;
  te=B2_scheme[e->B2_index].rpc;

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
        if (i>len2) Schlendrian("ecm_step2 %u %u %u\n",i,len2,ind);
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
        if (i>len1) Schlendrian("ecm_step2\n");
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
  te=B2_scheme[e->B2_index].tests;
  ind=0;
  for (i=0; i<len2; i++) {
    j=0;
    while (te[ind]) {
      j+=(u32_t)(te[ind++]);
      if (j>len1) Schlendrian("ecm_step2\n");
/* test tab1[j-1] against tab2[i] */
      asm_sub(mm_w,tab1[j-1],t2[2*i+2]);
      asm_mulmod(mm_prod,mm_prod,mm_w);
    }
    ind++;
  }
  free(tab1);
  free(invptr);
  if (asm_invert(mm_w,mm_prod)==0) {
    gcd(mm_a,mm_prod,montgomery_modulo_n);
    mpz_set_v(gmp_f,mm_a);
    return 1;
  }
  return 0;
}

/* ------------------ main functions --------------------- */

void ecm_set_params(ecm_t e, u32_t B1, u32_t B2)
{
  u32_t i;
  size_t mem;

/* todo: avoid multiple set_montgomery_multiplication calls */
  if (!set_montgomery_multiplication(e->N))
    Schlendrian("ecm_set_params without prior ecm_curve_init?\n");
  if (B2<=B1) B2=2*B1; /* at least one prime between B1 and B2 */
  if ((e->B1==B1) && (e->B2_index<(u32_t)B2_scheme_len) &&
      (B2_scheme[e->B2_index].B1==B1) &&
      (B2_scheme[e->B2_index].B2==B2)) return;

  if (B1_max<B1) ecm_init_B1(B1,B2);
  e->B1=B1;

  for (i=0; i<(u32_t)B2_scheme_len; i++)
    if ((B2_scheme[i].B1==B1) && (B2_scheme[i].B2==B2))
      break;
  if (i==(u32_t)B2_scheme_len) ecm_init_B2(B1,B2);
  e->B2_index=i;
}


void ecm_curve_init(ecm_t e)
{
  if (!ecm_is_init) ecm_init();
  mpz_init(e->N);
  mpz_init(e->curve_x);
  mpz_init(e->curve_y);
}


int ecm_curve_set(ecm_t e, mpz_t N, u32_t B1, u32_t B2)
{
  if (!set_montgomery_multiplication(N)) return -1;
  if (B2<=B1) return -1;
#ifdef ECM_ZEIT
zeita(12);
#endif
  mpz_set(e->N,N);
  mpz_set_ui(e->curve_x,0);
  ecm_set_params(e,B1,B2);
#ifdef ECM_ZEIT
zeitb(12);
#endif
  return 0;
}


int ecm_reset_n(ecm_t e, mpz_t N)
{
#if 1
  if (!mpz_divisible_p(e->N,N)) return -1;
#endif
  if (!set_montgomery_multiplication(N)) return -1;
  mpz_set(e->N,N);
  return 0;
}


void ecm_curve_clear(ecm_t e)
{
  mpz_clear(e->N);
  mpz_clear(e->curve_x);
  mpz_clear(e->curve_y);
}


int ecm(ecm_t e, mpz_t **fptr)
{
  if (ecm_next_curve(e)) { *fptr=&gmp_f; return 1; }
  ecm_init_pointers(e);
#ifdef ECM_ZEIT
zeitA(1);
#endif
  if (ecm_step1(e)) {
    *fptr=&gmp_f;
#ifdef ECM_ZEIT
zeitB(1);
#endif
    return 1;
  }
#ifdef ECM_ZEIT
zeitB(1);
zeitA(2);
#endif
  if (ecm_step2(e)) {
    *fptr=&gmp_f;
#ifdef ECM_ZEIT
zeitB(2);
#endif
    return 1;
  }
#ifdef ECM_ZEIT
zeitB(2);
#endif
  return 0;
}


int ecm_factor(mpz_t N, u32_t B1, u32_t B2, mpz_t **fptr, u32_t ncurves)
{
  u32_t i;
  int suc;
  ecm_t e;

  ecm_curve_init(e);
  if (ecm_curve_set(e,N,B1,B2)) return -1;
  for (i=0; i<ncurves; i++) {
    suc=ecm(e,fptr);
    if (suc) {
      ecm_curve_clear(e);
      return (int)(i+1); /* number of processed curves */
    }
  }
  ecm_curve_clear(e);
  return 0;
}

