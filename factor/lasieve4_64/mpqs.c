/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/
/* Written by T. Kleinjung, modified by J. Franke */
/* CJM, 11/30/04: This is a hack of generic/mpqs.c for easy compilation
   and inclusion in GGNFS.
*/
/*
#define xmalloc malloc
#define MPQS_STAT
#define MPQS_ZEIT
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#if defined( _MSC_VER ) && defined( _WIN64 )
#include "mpir.h"
#include "asm/if.h"
#else
#include <gmp.h>
#include "if.h"
#endif

#include "asm/siever-config.h"

#ifdef MPQS_STAT
extern u64_t stat_asm_eval, stat_asm_td;
extern u32_t stat_asm_div, stat_final_mulmod;
#endif

#define mpz_get_ull mpz_get_ui
#define mpz_set_ull mpz_set_ui

#define MPQS_MAX_FBSIZE      512
#define MPQS_MIN_EXCESS      10
#define MPQS_MAX_ADIV_ALL    6
#define MPQS_SIEVELEN        (1<<(L1_BITS<15?L1_BITS:15))
#define MPQS_REL_ENTRIES     32
#define MPQS_TD_MAX_NDIV     (MPQS_REL_ENTRIES-5)
#define MPQS_MAX_NRELS       256  /* ??? */
#define MPQS_MAX_NPRELS      2048  /* ??? */
#define MPQS_MAX_NCRELS      256  /* ??? */
#define MPQS_MAX_NRELS_BUF   1024  /* ??? */
#define MPQS_MIN_RELBUFFER   16
#define MPQS_GAUSS_MAX       512
#define MPQS_MAX_FACTORS     5
#define MPQS_MAX_NPRIMES     1024
#define MPQS_HASH_OVERFLOW   65535
#define MPQS_MULT_NTESTPR    25       /* prime(MPQS_MULT_NTESTPR)<2^8 */
#define MPQS_MULT_NCAND      8        /* largest candidate<2^8 */
#define MPQS_FB_MAXPRIME     4096     /* used in sieve and final */

/* common with asm functions */
u16_t mpqs_nFBk_1;
u16_t mpqs_td_begin, mpqs_sievebegin;
u16_t mpqs_FB_inv_info[4*MPQS_MAX_NPRIMES];
u16_t mpqs_nFBk, mpqs_FBk[3];  /* 3*5*7*11 too big */
unsigned char mpqs_FB_log[MPQS_MAX_FBSIZE];
u16_t mpqs_nFB, mpqs_nAdiv_total;
u16_t mpqs_Adiv_all[MPQS_MAX_ADIV_ALL];
u32_t mpqs_FBk_inv[3], mpqs_FB_A_inv[MPQS_MAX_ADIV_ALL];
/* end */

static int mpqs_isinit=0;
static mpz_t mpqs_gmp_help, mpqs_dummy, mpqs_help;
static mpz_t mpqs_N, mpqs_kN, mpqs_factors[MPQS_MAX_FACTORS];
static unsigned char mpqs_f_status[MPQS_MAX_FACTORS];
static u32_t mpqs_multiplier, mpqs_nfactors;
static double mpqs_complexity;
static u16_t mpqs_prime_table[MPQS_MAX_NPRIMES];
static u16_t mpqs_prime_sqrthelp[MPQS_MAX_NPRIMES];
static u16_t mpqs_nAdiv, mpqs_accept;
static unsigned char mpqs_Adiv_log[MPQS_MAX_ADIV_ALL];;
static u16_t mpqs_pmax;
static u16_t mpqs_Adiv[MPQS_MAX_ADIV_ALL];
static u16_t mpqs_Adiv_sqrt[MPQS_MAX_ADIV_ALL], mpqs_Adiv_all_sqrt[MPQS_MAX_ADIV_ALL];
static u16_t mpqs_Adiv_start1[MPQS_MAX_ADIV_ALL], mpqs_Adiv_start2[MPQS_MAX_ADIV_ALL];
static unsigned char mpqs_Adiv_active[MPQS_MAX_ADIV_ALL];
static unsigned char mpqs_A_mask[(1<<MPQS_MAX_ADIV_ALL)];
static u16_t mpqs_nA;
static i16_t mpqs_A_index, mpqs_B_index;
static u16_t mpqs_SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_FBSIZE];
static u16_t mpqs_2A_all_inv[MPQS_MAX_FBSIZE];
static u16_t mpqs_Adiv_SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_ADIV_ALL];
static i64_t mpqs_A, mpqs_B, mpqs_C, mpqs_2B, mpqs_Bi[MPQS_MAX_ADIV_ALL];
static u32_t mpqs_disp;
static u16_t mpqs_nsurvivors;
static u32_t mpqs_nrels, mpqs_nlp, mpqs_nsp, mpqs_excess;
static unsigned char mpqs_logbound;
static u16_t mpqs_nprels, mpqs_ncrels;
static u16_t **mpqs_relations, **mpqs_part_rels, **mpqs_comb_rels;
static u16_t **mpqs_rel_buffer;
static u16_t mpqs_rel_hash_table[256][16];
static u16_t mpqs_hash_table[128][16];
static u32_t mpqs_hash_index[128][16];
static u32_t mpqs_lp[128*15];
static u32_t **mpqs_gauss_row;
static unsigned char mpqs_sol[MPQS_GAUSS_MAX];
static i16_t mpqs_gauss_c[MPQS_GAUSS_MAX], mpqs_gauss_d[MPQS_GAUSS_MAX];
static long mpqs_gauss_m, mpqs_gauss_n, mpqs_gauss_n32, mpqs_gauss_k;
static i16_t mpqs_exp[128*16+MPQS_MAX_FBSIZE+3+MPQS_MAX_ADIV_ALL];
static mpz_t mpqs_sq1, mpqs_sq2;
#if 0
static u16_t mpqs_zero1, mpqs_zero2;
#endif
static double mpqs_kN_dbl;
static u64_t mpqs_kN_64, mpqs_A_inv_64;

u16_t mpqs_FB[2*MPQS_MAX_FBSIZE];
u16_t mpqs_FB_start[2*MPQS_MAX_FBSIZE];
u32_t mpqs_sievelen;
u32_t mpqs_FB_inv[MPQS_MAX_FBSIZE];
unsigned char *mpqs_sievearray;

unsigned char mpqs_256_inv_table[128]={
1, 171, 205, 183, 57, 163, 197, 239, 
241, 27, 61, 167, 41, 19, 53, 223, 
225, 139, 173, 151, 25, 131, 165, 207, 
209, 251, 29, 135, 9, 243, 21, 191, 
193, 107, 141, 119, 249, 99, 133, 175, 
177, 219, 253, 103, 233, 211, 245, 159, 
161, 75, 109, 87, 217, 67, 101, 143, 
145, 187, 221, 71, 201, 179, 213, 127, 
129, 43, 77, 55, 185, 35, 69, 111, 
113, 155, 189, 39, 169, 147, 181, 95, 
97, 11, 45, 23, 153, 3, 37, 79, 
81, 123, 157, 7, 137, 115, 149, 63, 
65, 235, 13, 247, 121, 227, 5, 47, 
49, 91, 125, 231, 105, 83, 117, 31, 
33, 203, 237, 215, 89, 195, 229, 15, 
17, 59, 93, 199, 73, 51, 85, 255 };

#ifdef MPQS_STAT
u32_t stat_mpqs_nsieves, stat_mpqs_nsurvivors, stat_mpqs_ntrials, stat_mpqs_ndiv;
#endif



static u64_t mpqs_inv_64(u64_t a)
{
  u64_t inv, h;

  inv=(u64_t)mpqs_256_inv_table[(a&0xff)>>1];
  h=a*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffff0000U;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffffffff00000000ULL;
  h*=inv; inv-=h;
  /*if (inv*a!=1ULL) Schlendrian("mpqs_inv_64");*/
  return inv;
}


static unsigned char mpqs_jacobi_tab2[8]={0,0,0,1,0,1,0,0};

static long mpqs_jacobi(u16_t a, u16_t b) /* always gcd(a,b)=1 and 0<a<b */
{
  long e, l, m, n, r;

  m=(long)a; n=(long)b;
  if (!(n&1)) Schlendrian("mpqs_jacobi: b even\n");
  e=0;
  if (!(m&1)) {
    l=mpqs_jacobi_tab2[n&7];
    do {
      m/=2; e+=l;
    } while (!(m&1));
  }
  while (m>1) {
    r=n-m;
    if (r>0) {
      l=m&n;
      l>>=1;
      e+=l;
      n=m; m=r;
    } else {
      m=-r;
    }
    if (m&1) continue;
    l=mpqs_jacobi_tab2[n&7];
    do {
      m/=2; e+=l;
    } while (!(m&1));
  }
  if (!m) Schlendrian("jacobi\n");
  e=e&1; e=1-2*e;
  return e;
}


static long mpqs_jacobi0(u16_t a, u16_t b) /* always gcd(a,b)=1 and 0<a<b */
{
  long e, l, m, n;

  m=(long)a; n=(long)b;
  if (!(n&1)) Schlendrian("mpqs_jacobi: b even\n");
  e=0;
  if (!(m&1)) {
    l=mpqs_jacobi_tab2[n&7];
    do {
      m/=2; e+=l;
    } while (!(m&1));
  }
  while (m>1) {
    l=m; m=n; n=l;
    l=m&n;
    l>>=1;
    e+=l;
    l=m/n; m-=l*n;
    if (m&1) continue;
    l=mpqs_jacobi_tab2[n&7];
    do {
      m/=2; e+=l;
    } while (!(m&1));
  }
  if (!m) Schlendrian("jacobi\n");
  e=e&1; e=1-2*e;
  return e;
}


static u16_t mpqs_invert(u16_t a, u16_t p)
{
  u32_t v1=0, v2=1;
  u32_t b=(u32_t)a, pp=(u32_t)p, q;

  while (b>1) {
    pp-=b; v1+=v2;
    if (pp>=b) {
      pp-=b; v1+=v2;
      if (pp>=b) {
        q=pp/b; pp%=b;
        v1+=q*v2;
      }
    }
    if (pp<=1) { v2=p-v1; break; }
    b-=pp; v2+=v1;
    if (b>=pp) {
      b-=pp; v2+=v1;
      if (b>=pp) {
        q=b/pp; b%=pp;
        v2+=q*v1;
      }
    }
  }
/*  if ((v2*a-1)%(u32_t)p) complain("invert");*/
  return v2;
}


static u16_t mpqs_powmod(u16_t a, u16_t e, u16_t p)
{
  u16_t ex=e;
  u32_t aa, res;

  if (!ex) return 1;
  aa=(u32_t)a;
  res=1;
  while (1) {
    if (ex&1) { res*=aa; res%=(u32_t)p; }
    ex>>=1;
    if (!ex) return (u16_t)res;
    aa*=aa; aa%=(u32_t)p;
  }
  return 0; /* never reached */
}


static u16_t mpqs_sqrt_init(u16_t p)
{
  u16_t e, u, i, j;
  u32_t g=0, b;

  if (p&2) return 0;
  if (p&4) return mpqs_powmod(2,(p-1)>>2,p);
  e=0; u=p-1;
  while (!(u&1)) { e++; u>>=1; }
  for (i=2; i<p; i++) {
    g=mpqs_powmod(i,u,p);
    if (g==1) continue;
    b=g;
    for (j=0; j<e-1; j++) { b*=b; b%=(u32_t)p; }
    if (b==p-1) break;
  }
  if (i>=p) Schlendrian("mpqs_sqrt_init\n");
  return g;
}


static u16_t mpqs_sqrt(u16_t a, u16_t p, u16_t help)
                 /* a!=0, a<p, p prime, p!=2, (a/p)=1 */
{
  u16_t e, u, i, l;
  u32_t b, g, r, k;

  if (p&2) return mpqs_powmod(a,(p+1)>>2,p);
  if (p&4) {
    b=mpqs_powmod(a,(p+3)>>3,p);
    g=b*b; g%=(u32_t)p;
    if (g==a) return b;
    b*=(u32_t)help; b%=(u32_t)p;
    return b;
  }
  e=0; u=p-1;
  while (!(u&1)) { e++; u>>=1; }
  g=help;
#if 1
  r=mpqs_powmod(a,u>>1,p);
  k=r*r; k%=(u32_t)p;
  k*=a; k%=(u32_t)p;
  r*=a; r%=(u32_t)p;
#else
  r=mpqs_powmod(a,(u+1)/2,p);
  k=mpqs_powmod(a,u,p);
#endif
  while (k!=1) {
    b=k;
    for (i=0; i<e && b!=1; i++) { b*=b; b%=(u32_t)p; }
    for (l=i+1; l<e; l++) { g*=g; g%=(u32_t)p; }
    r*=g; r%=(u32_t)p;
    g*=g; g%=(u32_t)p;
    k*=g; k%=(u32_t)p;
    e=i;
  }
  return (u16_t)r;
}

#if 0
static u32_t mpqs_inv(u16_t i)
{
  u32_t v1=0, v2=1, q, a, p;

  p=0; a=i;
  while (a>1) {
    p-=a; v1+=v2;
    if (p>=a) {
      p-=a; v1+=v2;
      if (p>=a) {
        q=p/a; p%=a;
        v1+=q*v2;
      }
    }
    if (p<=1) { v2=-v1; break; }
    a-=p; v2+=v1;
    if (a>=p) {
      a-=p; v2+=v1;
      if (a>=p) {
        q=a/p; a%=p;
        v2+=q*v1;
      }
    }
  }
/*  if (v2*i!=1) complain("inv: %lu %lu %lu %lu %d ",v1,v2,a,p,i);
  if (v2*i!=1) Schlendrian("mpqs_inv");*/
  return v2;
}
#else
static u32_t mpqs_inv(u16_t a)
{
  u32_t inv, h;

  inv=(u32_t)mpqs_256_inv_table[(a&0xff)>>1];
  h=(u32_t)a*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=(u32_t)a*inv; h&=0xffff0000;
  h*=inv; inv-=h;
  if (inv*(u32_t)a!=1ULL) Schlendrian("mpqs_inv");
  return inv;
}
#endif

/* ---------------------------------------------------------- */

void mpqs_init()
{
  long i, j;
  u32_t p, add, d;

  mpz_init(mpqs_N);
  mpz_init(mpqs_kN);
  mpz_init(mpqs_dummy);
  mpz_init(mpqs_gmp_help);
  mpz_init(mpqs_help);
  mpz_init(mpqs_sq1);
  mpz_init(mpqs_sq2);
  for (i=0; i<MPQS_MAX_FACTORS; i++) mpz_init(mpqs_factors[i]);
  mpqs_sievearray=(unsigned char *)xmalloc(MPQS_SIEVELEN*sizeof(unsigned char));
  mpqs_relations=(u16_t **)xmalloc(MPQS_MAX_NRELS*sizeof(u16_t *));
  mpqs_relations[0]=(u16_t *)xmalloc(MPQS_MAX_NRELS*MPQS_REL_ENTRIES*sizeof(u16_t));
  for (i=1; i<MPQS_MAX_NRELS; i++)
    mpqs_relations[i]=mpqs_relations[i-1]+MPQS_REL_ENTRIES;
  mpqs_part_rels=(u16_t **)xmalloc(MPQS_MAX_NPRELS*sizeof(u16_t *));
  mpqs_part_rels[0]=(u16_t *)xmalloc(MPQS_MAX_NPRELS*MPQS_REL_ENTRIES*sizeof(u16_t));
  for (i=1; i<MPQS_MAX_NPRELS; i++)
    mpqs_part_rels[i]=mpqs_part_rels[i-1]+MPQS_REL_ENTRIES;
  mpqs_comb_rels=(u16_t **)xmalloc(MPQS_MAX_NCRELS*sizeof(u16_t *));
  mpqs_comb_rels[0]=(u16_t *)xmalloc(MPQS_MAX_NCRELS*2*MPQS_REL_ENTRIES*sizeof(u16_t));
  for (i=1; i<MPQS_MAX_NCRELS; i++)
    mpqs_comb_rels[i]=mpqs_comb_rels[i-1]+2*MPQS_REL_ENTRIES;
  mpqs_rel_buffer=(u16_t **)xmalloc(MPQS_MAX_NRELS_BUF*sizeof(u16_t *));
  mpqs_rel_buffer[0]=(u16_t *)xmalloc(MPQS_MAX_NRELS_BUF*MPQS_REL_ENTRIES*sizeof(u16_t));
  for (i=1; i<MPQS_MAX_NRELS_BUF; i++)
    mpqs_rel_buffer[i]=mpqs_rel_buffer[i-1]+MPQS_REL_ENTRIES;

  mpqs_gauss_row=(u32_t **)xmalloc(MPQS_GAUSS_MAX*sizeof(u32_t *));
  mpqs_gauss_row[0]=(u32_t *)xmalloc(MPQS_GAUSS_MAX*MPQS_GAUSS_MAX/32*sizeof(u32_t));
#if 0
  prime_table_init();
  for (i=0; i<MPQS_MAX_NPRIMES; i++) {
    mpqs_prime_table[i]=(u16_t)get_next_prime();
    mpqs_prime_sqrthelp[i]=mpqs_sqrt_init(mpqs_prime_table[i]);
  }
#else
  mpqs_prime_table[0]=2;
  mpqs_prime_table[1]=3;
  mpqs_prime_table[2]=5; mpqs_prime_sqrthelp[2]=mpqs_sqrt_init(5);
  i=3; p=7; add=4;
  while (i<MPQS_MAX_NPRIMES) {
    for (j=2; j<i; j++) {
      d=(u32_t)mpqs_prime_table[j];
      if (d*d>=p) break;
      if (p%d==0) break;
    }
    if (d*d>p) {
      mpqs_prime_table[i]=(u16_t)p;
      mpqs_prime_sqrthelp[i]=mpqs_sqrt_init(p);
      i++;
    }
    p+=add; add=6-add;
  }
#endif

  mpqs_isinit=1;
}


static u16_t mpqs_multiplier_cand[4][MPQS_MULT_NCAND]={
 { 1, 17, 33, 41, 57, 65, 73, 89 },
 { 3, 11, 19, 35, 43, 51, 59, 67 },
 { 5, 13, 21, 29, 37, 53, 61, 69 },
 { 7, 15, 23, 31, 39, 47, 55, 71 },
};

static void mpqs_choose_multiplier()
{
  u16_t Nmod8, mult, p, mm;
  u16_t residue[MPQS_MULT_NTESTPR];
  double value[MPQS_MULT_NCAND], v, vp, vmin, dn;
  long i, j;  

  for (i=0; i<MPQS_MULT_NTESTPR; i++)
    residue[i]=(u16_t)mpz_mod_ui(mpqs_dummy,mpqs_N,(u32_t)mpqs_prime_table[i]);
  Nmod8=(u16_t)mpz_mod_ui(mpqs_dummy,mpqs_N,8);
  dn=log(mpz_get_d(mpqs_N))/log(2.);
  dn=100.9-dn;
  if (dn>7) mm=128; else mm=(u16_t)(exp(dn*log(2)));
  for (j=0; j<MPQS_MULT_NCAND; j++) {
    mult=mpqs_multiplier_cand[Nmod8/2][j];
    if (mult<=mm) {
      v=log((double)mult)-1.5*log(2);
      for (i=1; i<MPQS_MULT_NTESTPR; i++) {
        p=mpqs_prime_table[i]; vp=log((double)p)/((double)p);
        if (mult%p) {
          if (mpqs_jacobi((mult*residue[i])%p,p)==1) v-=vp;
          else v+=vp;
        }
      }
    } else v=1000;
    value[j]=v;
  }
  vmin=value[0];
  for (j=1; j<MPQS_MULT_NCAND; j++)
    if (value[j]<vmin) vmin=value[j];
  for (j=0; j<MPQS_MULT_NCAND; j++)
    if (value[j]==vmin) break;
  if (j>=MPQS_MULT_NCAND) complain("mpqs_choose_multiplier.vmin\n");
  mpqs_complexity=(vmin+log(mpz_get_d(mpqs_N)))/log(2.);
  mpqs_multiplier=mpqs_multiplier_cand[Nmod8/2][j];
#ifdef MPQS_STAT
printf("%ld ",mpqs_multiplier);
#endif
  mpz_mul_ui(mpqs_kN,mpqs_N,mpqs_multiplier);
#if 0
  mpz_fdiv_q_2exp(mpqs_dummy,mpqs_kN,32);
  mpqs_kN_64=(u64_t)mpz_get_ui(mpqs_dummy); mpqs_kN_64<<=32;
  mpqs_kN_64+=(u64_t)mpz_get_ui(mpqs_kN);
#else
  mpqs_kN_64=mpz_get_ull(mpqs_kN);
#endif
}


static void mpqs_choose_multiplier0()
{
  u32_t Nmod8;

  Nmod8=mpz_mod_ui(mpqs_dummy,mpqs_N,8);
  mpqs_multiplier=Nmod8;  /* vorläufig !!!!!!!!!!!!!!!!!!!!!! */
  mpz_mul_ui(mpqs_kN,mpqs_N,mpqs_multiplier);
  mpz_fdiv_q_2exp(mpqs_dummy,mpqs_kN,32);
  mpqs_kN_64=(u64_t)mpz_get_ui(mpqs_dummy); mpqs_kN_64<<=32;
  mpqs_kN_64+=(u64_t)mpz_get_ui(mpqs_kN);
}

#if 1
static u16_t mpqs_param[12][6]={
 { 40, 3, 2, 4, 11, 16},
 { 50, 3, 2, 4, 12, 16},
 { 60, 3, 3, 4, 15, 16},
 { 70, 3, 3, 5, 14, 16},
 { 80, 3, 3, 5, 14, 16},
 { 90, 3, 3, 5, 15, 20},
 { 110, 3, 3, 5, 17, 20},
 { 120, 3, 3, 5, 19, 20},
 { 140, 3, 4, 6, 18, 30},
 { 140, 3, 4, 6, 20, 40},
 { 160, 3, 4, 6, 21, 50},
 { 180, 4, 4, 6, 23, 70}
};
/* vor retry-implementierung
static u16_t mpqs_param[12][6]={
 { 40, 3, 2, 4, 11, 16},
 { 50, 3, 2, 4, 12, 16},
 { 60, 3, 3, 4, 15, 16},
 { 70, 3, 3, 5, 14, 16},
 { 80, 3, 3, 5, 14, 16},
 { 90, 3, 3, 5, 15, 20},
 { 110, 3, 3, 5, 17, 20},
 { 120, 3, 3, 5, 19, 20},
 { 140, 3, 4, 6, 18, 30},
 { 140, 3, 4, 6, 20, 40},
 { 160, 3, 4, 6, 21, 50},
 { 180, 4, 4, 6, 23, 70}
};
*/
#else
static u16_t mpqs_param[10][6]={
 { 60, 3, 3, 4, 13, 16},
 { 70, 3, 3, 5, 14, 16},
 { 80, 3, 3, 5, 14, 16},
 { 90, 3, 3, 5, 15, 20},
 { 110, 3, 3, 5, 16, 20},
 { 120, 3, 3, 5, 18, 20},
 { 130, 3, 4, 6, 18, 30},
 { 130, 3, 4, 6, 19, 30},
 { 150, 3, 4, 6, 20, 50},
 { 180, 4, 4, 6, 23, 70}
};
#endif

static void mpqs_choose_parameter(long retry)
{
  long n, r;

#if 1
  n=(long)(mpqs_complexity+0.5);
  if (n>96) n=96;
  if (n<51) n=51;
  n=(n+3)/4;
  n-=13;  /* 0<=n<=11 */
#else
  n=(long)(mpqs_complexity+0.5);
  if (n>96) n=96;
  if (n<59) n=59;
  n=(n+3)/4;
  n-=15;  /* 0<=n<=9 */
#endif

  r=retry;
  if ((r) && (n<9)) { n++; r--; }
  mpqs_nFB=mpqs_param[n][0];
  if (r) mpqs_nFB+=r*mpqs_nFB/4;
  if (mpqs_nFB>=MPQS_MAX_FBSIZE) mpqs_nFB=MPQS_MAX_FBSIZE-1;
  mpqs_sievebegin=mpqs_param[n][1];
  if (r) mpqs_sievebegin+=r;
  mpqs_nAdiv=mpqs_param[n][2];
  mpqs_nAdiv_total=mpqs_param[n][3];
  mpqs_accept=mpqs_param[n][4];
  mpqs_td_begin=mpqs_param[n][5];

  mpqs_sievelen=MPQS_SIEVELEN;  /* !!! */
  mpqs_disp=mpqs_sievelen/2;
  mpqs_nfactors=0; mpqs_nrels=0; mpqs_nlp=0; mpqs_nprels=0; mpqs_ncrels=0;
  for (n=0; n<128; n++) mpqs_hash_table[n][0]=0;
  for (n=0; n<256; n++) mpqs_rel_hash_table[n][0]=0;
}


static void mpqs_generate_FB()
{
  u16_t *fb, p, rest, i, nfb;

  fb=mpqs_FB; nfb=0; i=0; mpqs_nFBk=0;
  *fb++=2; *fb++=1; nfb++;
  for (i=1; i<MPQS_MAX_NPRIMES; i++) {
    p=mpqs_prime_table[i];
    if (p>MPQS_FB_MAXPRIME) break;
    rest=(u16_t)mpz_mod_ui(mpqs_dummy,mpqs_kN,p);  /* ca 5M */
    if (rest) {
      if (mpqs_jacobi(rest,p)==1) { /* ca 10M */
	    assert(p);
        *fb++=p;
        *fb++=mpqs_sqrt(rest,p,mpqs_prime_sqrthelp[i]); /* ca 13M */
        nfb++;
        if (nfb>=mpqs_nFB) { break; }
      }
    } else {
      if (mpqs_multiplier%(u32_t)p)
        complain("mpqs: N has small divisor: %u\n",p);
      mpqs_FBk[mpqs_nFBk++]=p;
    }
  }
  mpqs_nFB=nfb;
  mpqs_nFBk_1=mpqs_nFBk+1;
  mpqs_nsp=1+mpqs_nFB+mpqs_nFBk;
  mpqs_pmax=mpqs_FB[2*(nfb-1)];
  mpqs_FB[2*nfb]=mpqs_sievelen;
}


static int mpqs_SI_init()
{
  double d;
  long a, i, j, k, l, n;
  i64_t prod;
  u16_t inv, p, *fb=mpqs_FB;
  double A_div_log[MPQS_MAX_ADIV_ALL]; /* maximal number of A's */
  double v;

  d=mpz_get_d(mpqs_kN); mpqs_kN_dbl=d;
  d*=8.; d=sqrt(d); d/=(double)mpqs_sievelen;
  d=log(d); d/=(double)mpqs_nAdiv; d=exp(d);
  a=(long)d;
/*  if (a>=mpqs_pmax) { printf("%ld %ld ",a,mpqs_pmax); return -1; }*/
  if (a>=mpqs_pmax) a=mpqs_pmax-1;
  for (i=0; i<mpqs_nFB; i++) if (mpqs_FB[2*i]>a) break;
  i-=mpqs_nAdiv_total/2;
  if (i<1) i=1; /* first prime >2 */
  if (mpqs_FB[2*i]<3) return -2;
  if (i+mpqs_nAdiv_total>mpqs_nFB) /*return -3;*/ i=mpqs_nFB-mpqs_nAdiv_total-1;
  for (j=0; j<mpqs_nAdiv_total; j++) {
    mpqs_Adiv_all[j]=mpqs_FB[2*(i+j)];
    mpqs_Adiv_all_sqrt[j]=mpqs_FB[2*(i+j)+1];
  }
  for (j=i+mpqs_nAdiv_total; j<mpqs_nFB; j++) {
    mpqs_FB[2*(j-mpqs_nAdiv_total)]=mpqs_FB[2*j];
    mpqs_FB[2*(j-mpqs_nAdiv_total)+1]=mpqs_FB[2*j+1];
  }
  mpqs_nFB-=mpqs_nAdiv_total;
  if (i<mpqs_sievebegin) mpqs_sievebegin=i;
#if 1
  p=mpqs_FB[2];
  mpqs_FB_inv[1]=mpqs_inv(p);
  mpqs_FB_inv_info[2]=p; mpqs_FB_inv_info[3]=p;
  mpqs_FB_inv_info[6]=(u16_t)mpqs_FB_inv[1];
  mpqs_FB_inv_info[7]=(u16_t)mpqs_FB_inv[1];
  if (mpqs_td_begin&1) complain("mpqs_td_begin is odd\n");
  for (j=2; j<mpqs_td_begin; j+=2) {
    p=mpqs_FB[2*j];
    mpqs_FB_inv[j]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j]=p; mpqs_FB_inv_info[4*j+1]=p;
    mpqs_FB_inv_info[4*j+4]=(u16_t)mpqs_FB_inv[j];
    mpqs_FB_inv_info[4*j+5]=(u16_t)mpqs_FB_inv[j];
    p=mpqs_FB[2*(j+1)];
    mpqs_FB_inv[j+1]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j+2]=p; mpqs_FB_inv_info[4*j+3]=p;
    mpqs_FB_inv_info[4*j+6]=(u16_t)mpqs_FB_inv[j+1];
    mpqs_FB_inv_info[4*j+7]=(u16_t)mpqs_FB_inv[j+1];
  }
  for (; j<mpqs_nFB; j++) mpqs_FB_inv[j]=mpqs_inv(mpqs_FB[2*j]);
  for (j=0; j<mpqs_nFBk; j++) mpqs_FBk_inv[j]=mpqs_inv(mpqs_FBk[j]);
  for (j=0; j<mpqs_nAdiv_total; j++)
    mpqs_FB_A_inv[j]=mpqs_inv(mpqs_Adiv_all[j]);
#endif
/* compute log-approximations */
  d=mpz_get_d(mpqs_kN);
  d/=8.; d=sqrt(d); d*=(double)mpqs_sievelen;
  d=log(d);
  d-=log((double)((u64_t)1<<mpqs_accept));
  mpqs_logbound=128;  /* ??? */
  d=(double)(mpqs_logbound)/d;
  for (i=0; i<mpqs_nFB; i++) {
    p=mpqs_FB[2*i];
    mpqs_FB_log[i]=(unsigned char)(0.5+d*log((double)p));
  }
  for (i=0; i<mpqs_nAdiv_total; i++) {
    p=mpqs_Adiv_all[i];
    mpqs_Adiv_log[i]=(unsigned char)(0.5+d*log((double)p));
  }

/* compute mask for choice of A-divisors */
  n=1<<mpqs_nAdiv_total;
  for (i=0; i<n; i++) mpqs_A_mask[i]=0;
  for (i=0; i<mpqs_nAdiv_total; i++) {
    k=1<<i;
    for (j=k; j<n; j+=(2*k))
      for (l=0; l<k; l++) mpqs_A_mask[j+l]++;
  } /* now mpqs_A_mask[i] contains the number of ones of i (binary) */
  j=0;
  for (i=0; i<n; i++)
    if (mpqs_A_mask[i]==mpqs_nAdiv) {
      mpqs_A_mask[j]=i; j++;
    }
/* ensure that |C|<2^64 */
  for (i=0; i<mpqs_nAdiv_total; i++)
    A_div_log[i]=log((double)(mpqs_Adiv_all[i]))/log(2.);
  d=log(mpqs_kN_dbl)/log(2.);
  for (i=0; i<j; i++) {
    v=0;
    for (k=0; k<mpqs_nAdiv_total; k++)
      if (mpqs_A_mask[i]&(1<<k)) v+=A_div_log[k];
    if (d-v>63.95) {
      j--; /* printf(":");*/
      mpqs_A_mask[i]=mpqs_A_mask[j];
      i--;
    }
  }

  mpqs_nA=j; mpqs_A_index=-1; mpqs_B_index=(1<<mpqs_nAdiv)-1;

  prod=1; d=1.;
  for (i=0; i<mpqs_nAdiv_total; i++) {
    prod*=(i64_t)(mpqs_Adiv_all[i]);
    d*=(double)(mpqs_Adiv_all[i]);
  }
  if (d>9223372036854775807.) return -4;
  for (i=1; i<mpqs_nFB; i++) {
    p=fb[2*i];
    inv=(u16_t)(prod%(i64_t)p);
    inv=mpqs_invert(inv,p);
    inv<<=1; if (inv>=p) inv-=p;
    mpqs_2A_all_inv[i]=inv;
  }
  mpqs_FB[2*mpqs_nFB]=mpqs_sievelen;
  return 1;
}


static int mpqs_next_pol()
{
  long i, ind, j;
  i64_t add, bi, prod_other;
  u16_t p, s1, s2, inv, *fb=mpqs_FB, *fbs=mpqs_FB_start;
  i16_t sh;
  u32_t bb, cc, c1, c2;
  unsigned char mask;
#if 0
/* Used for calculation of zeros. */
  double d;
#endif
  u16_t aa;
#define GCC_BUG
#ifndef GCC_BUG
  u32_t al, *ul;
#endif

  mpqs_B_index++;
  if (mpqs_B_index>=1<<(mpqs_nAdiv-1)) {
/*zeitA(7);*/
    mpqs_A_index++;
    if (mpqs_A_index>=mpqs_nA) return 0;
    mask=mpqs_A_mask[mpqs_A_index];
    for (i=0; i<mpqs_nAdiv_total; i++)
      if (mask & (1<<i)) mpqs_Adiv_active[i]=1;
      else mpqs_Adiv_active[i]=0;
    j=0;
    for (i=0; i<mpqs_nAdiv_total; i++)
      if (mpqs_Adiv_active[i]) {
        mpqs_Adiv[j]=mpqs_Adiv_all[i];
        mpqs_Adiv_sqrt[j]=mpqs_Adiv_all_sqrt[i];
        j++;
      }
    mpqs_A=1;
    for (i=0; i<mpqs_nAdiv; i++) mpqs_A*=(i64_t)(mpqs_Adiv[i]);
#if 0
/* compute zeros */
    d=sqrt(mpqs_kN_dbl)/(double)mpqs_A;
    if (d<(double)mpqs_disp) {
      mpqs_zero1=mpqs_disp-(u16_t)d;
      mpqs_zero2=mpqs_disp+(u16_t)d;
    } else {
      mpqs_zero1=mpqs_sievelen;
      mpqs_zero2=mpqs_sievelen;
    }
#endif
/* compute B_i's */
    for (i=0; i<mpqs_nAdiv; i++) {
      p=mpqs_Adiv[i];
      bi=mpqs_A/(i64_t)(p);
      inv=(u16_t)(bi%(i64_t)p); /* bi>0 */
      inv=mpqs_invert(inv,p);
      bi*=(i64_t)inv; bi%=mpqs_A;
      bi*=(i64_t)mpqs_Adiv_sqrt[i]; bi%=mpqs_A;
      mpqs_Bi[i]=bi;
    }
    mpqs_B=0;
    for (i=0; i<mpqs_nAdiv; i++) mpqs_B+=mpqs_Bi[i];
    prod_other=1;
    for (i=0; i<mpqs_nAdiv_total; i++)
      if (!mpqs_Adiv_active[i]) prod_other*=(i64_t)(mpqs_Adiv_all[i]);
    for (i=1; i<mpqs_nFB; i++) {
      p=fb[2*i];
	  assert(p);
#if 0
      bb=(u32_t)(prod_other%(i64_t)p);
      bb*=(u32_t)(mpqs_2A_all_inv[i]); bb%=(u32_t)p;
#else
bb=(u32_t)(((i64_t)(mpqs_2A_all_inv[i])*prod_other)%(u32_t)p);
#endif
      for (j=0; j<mpqs_nAdiv; j++) {
#if 1
        cc=bb*(u32_t)(mpqs_Bi[j]%(i64_t)p); cc%=(u32_t)p;
#else
cc=(u32_t)(((i64_t)bb*mpqs_Bi[j])%(u32_t)p);
#endif
        mpqs_SI_add[j][i]=(u16_t)cc;
      }
      if (bb&1) bb+=p;
      bb>>=1;  /* now bb=1/A mod p */
      cc=p-(u32_t)(mpqs_B%(i64_t)p);  /* mpqs_B>0 */
      c1=cc+fb[2*i+1]; c1*=bb; c1+=mpqs_disp; c1%=(u32_t)p;
      c2=cc+(p-fb[2*i+1]); c2*=bb; c2+=mpqs_disp; c2%=(u32_t)p;
#if 0
      if (c1<c2) { fbs[2*i]=(u16_t)c1; fbs[2*i+1]=(u16_t)c2; }
      else { fbs[2*i]=(u16_t)c2; fbs[2*i+1]=(u16_t)c1; }
#else
      fbs[2*i]=(u16_t)c1; fbs[2*i+1]=(u16_t)c2;
#endif
    }

    mpqs_2B=2*mpqs_B;
    mpqs_A_inv_64=mpqs_inv_64(mpqs_A);
#if 0
    mpz_set_sll(mpqs_help,mpqs_B);
    mpz_mul(mpqs_help,mpqs_help,mpqs_help);
    mpz_sub(mpqs_help,mpqs_help,mpqs_kN);
    mpz_set_sll(mpqs_dummy,mpqs_A);
    mpz_fdiv_qr(mpqs_help,mpqs_dummy,mpqs_help,mpqs_dummy);
    if (mpz_sgn(mpqs_dummy)) {
      mpz_out_str(stdout,10,mpqs_help); printf("\n");
      mpz_out_str(stdout,10,mpqs_dummy); printf("\n");
      Schlendrian("mpqs_next_pol\n");
    }
    mpz_get_sll(&mpqs_C,mpqs_help);
#else
    mpqs_C=mpqs_B*mpqs_B;
    mpqs_C-=mpqs_kN_64;
    mpqs_C*=mpqs_A_inv_64;
#endif

    for (i=0; i<mpqs_nAdiv_total; i++)
      if (!mpqs_Adiv_active[i]) {
        p=mpqs_Adiv_all[i];
        bb=(u32_t)(mpqs_A%(i64_t)p);
        bb=mpqs_invert(bb,p);
        bb<<=1; if (bb>=p) bb-=p;
        for (j=0; j<mpqs_nAdiv; j++) {
          cc=bb*(u32_t)(mpqs_Bi[j]%(i64_t)p); cc%=(u32_t)p;
          mpqs_Adiv_SI_add[j][i]=(u16_t)cc;
        }
        if (bb&1) bb+=p;
        bb>>=1;  /* now bb=1/A mod p */
        cc=p-(u32_t)(mpqs_B%(i64_t)p);  /* mpqs_B>0 */
        c1=cc+mpqs_Adiv_all_sqrt[i]; c1*=bb; c1+=mpqs_disp; c1%=(u32_t)p;
        c2=cc+(p-mpqs_Adiv_all_sqrt[i]); c2*=bb; c2+=mpqs_disp; c2%=(u32_t)p;
        if (c1<c2) {
          mpqs_Adiv_start1[i]=(u16_t)c1; mpqs_Adiv_start2[i]=(u16_t)c2;
        } else {
          mpqs_Adiv_start1[i]=(u16_t)c2; mpqs_Adiv_start2[i]=(u16_t)c1;
        }
      } else {
        p=mpqs_Adiv_all[i];
        sh=(i16_t)(mpqs_2B%(i64_t)p);
        if (sh>0) bb=(u32_t)sh; else bb=(u32_t)((i16_t)p+sh);
        bb=(u32_t)mpqs_invert(bb,p);
        sh=(i16_t)(mpqs_C%(i64_t)p); sh=-sh;
        if (sh>0) cc=(u32_t)sh; else cc=(u32_t)((i16_t)p+sh);
        cc*=bb; cc%=(u32_t)p;  /* now p|2*B*cc+C */
        cc+=mpqs_disp; cc%=(u32_t)p;
        mpqs_Adiv_start1[i]=(u16_t)cc;
      }
    mpqs_B_index=0;
/*zeitB(7);*/
    return 1;
  }
/*zeitA(7);*/
  ind=mpqs_B_index; i=0;
  while (1) {
    if (ind&1) break;
    i++; ind>>=1;
  }
  ind>>=1; ind++;
  add=2*mpqs_Bi[i];
  if (ind&1) {
    mpqs_B-=add;
    for (j=1; j<mpqs_nFB; j++) {
#if 0
      p=fb[2*j]; s1=fbs[2*j]; s2=fbs[2*j+1];
      s1+=mpqs_SI_add[i][j]; if (s1>=p) s1-=p;
      s2+=mpqs_SI_add[i][j]; if (s2>=p) s2-=p;
/*      if (s1<s2) {*/
        fbs[2*j]=s1; fbs[2*j+1]=s2;
/*      } else {
        fbs[2*j]=s2; fbs[2*j+1]=s1;
      }*/
#else
      p=fb[2*j]; /*s1=fbs[2*j]; s2=fbs[2*j+1];*/
      aa=mpqs_SI_add[i][j];
#ifndef GCC_BUG
      al=(((u32_t)aa)<<16)+(u32_t)aa;
      ul=(u32_t *)(fbs+2*j);
      *ul+=al;
#else
      fbs[2*j]+=aa; fbs[2*j+1]+=aa;
#endif
      if (fbs[2*j]>=p) fbs[2*j]-=p;
      if (fbs[2*j+1]>=p) fbs[2*j+1]-=p;
/*      s1+=aa; if (s1>=p) s1-=p;
      s2+=aa; if (s2>=p) s2-=p;
      fbs[2*j]=s1; fbs[2*j+1]=s2;*/
#endif
    }
  } else {
    mpqs_B+=add;
    for (j=1; j<mpqs_nFB; j++) {
#if 0
      p=fb[2*j]; s1=fbs[2*j]; s2=fbs[2*j+1];
      s1+=(p-mpqs_SI_add[i][j]); if (s1>=p) s1-=p;
      s2+=(p-mpqs_SI_add[i][j]); if (s2>=p) s2-=p;
/*      if (s1<s2) {*/
        fbs[2*j]=s1; fbs[2*j+1]=s2;
/*      } else {
        fbs[2*j]=s2; fbs[2*j+1]=s1;
      }*/
#else
      p=fb[2*j]; /*s1=fbs[2*j]; s2=fbs[2*j+1];*/
      aa=p-mpqs_SI_add[i][j];
#ifndef GCC_BUG
      al=(((u32_t)aa)<<16)+(u32_t)aa;
      ul=(u32_t *)(fbs+2*j);
      *ul+=al;
#else
      fbs[2*j]+=aa; fbs[2*j+1]+=aa;
#endif
      if (fbs[2*j]>=p) fbs[2*j]-=p;
      if (fbs[2*j+1]>=p) fbs[2*j+1]-=p;
/*      s1+=aa; if (s1>=p) s1-=p;
      s2+=aa; if (s2>=p) s2-=p;
      fbs[2*j]=s1; fbs[2*j+1]=s2;*/
#endif
    }
  }
/*zeitB(7);*/
  mpqs_2B=2*mpqs_B;
#if 0
  mpz_set_sll(mpqs_help,mpqs_B);
  mpz_mul(mpqs_help,mpqs_help,mpqs_help);
  mpz_sub(mpqs_help,mpqs_help,mpqs_kN);
  mpz_set_sll(mpqs_dummy,mpqs_A);
#if 0
  mpz_fdiv_qr(mpqs_help,mpqs_dummy,mpqs_help,mpqs_dummy);
  if (mpz_sgn(mpqs_dummy)) Schlendrian("mpqs_next_pol\n");
#else
  mpz_divexact(mpqs_help,mpqs_help,mpqs_dummy);
#endif
  mpz_get_sll(&mpqs_C,mpqs_help);
#else
  mpqs_C=mpqs_B*mpqs_B;
  mpqs_C-=mpqs_kN_64;
  mpqs_C*=mpqs_A_inv_64;
#endif
/*printf("pol: %Ld %Ld %Ld \n",mpqs_A,mpqs_B,mpqs_C);*/
  for (j=0; j<mpqs_nAdiv_total; j++)
    if (!mpqs_Adiv_active[j]) {
      p=mpqs_Adiv_all[j];
      s1=mpqs_Adiv_start1[j]; s2=mpqs_Adiv_start2[j];
      if (ind&1) {
        s1+=mpqs_Adiv_SI_add[i][j]; if (s1>=p) s1-=p;
        s2+=mpqs_Adiv_SI_add[i][j]; if (s2>=p) s2-=p;
      } else {
        s1+=(p-mpqs_Adiv_SI_add[i][j]); if (s1>=p) s1-=p;
        s2+=(p-mpqs_Adiv_SI_add[i][j]); if (s2>=p) s2-=p;
      }
      mpqs_Adiv_start1[j]=s1; mpqs_Adiv_start2[j]=s2;
    } else {
      p=mpqs_Adiv_all[i];
      sh=(i16_t)(mpqs_2B%(i64_t)p);
      if (sh>0) bb=(u32_t)sh; else bb=(u32_t)((i16_t)p+sh);
      bb=(u32_t)mpqs_invert(bb,p);
      sh=(i16_t)(mpqs_C%(i64_t)p); sh=-sh;
      if (sh>0) cc=(u32_t)sh; else cc=(u32_t)((i16_t)p+sh);
      cc*=bb; cc%=(u32_t)p;  /* now p|2*B*cc+C */
      cc+=mpqs_disp; cc%=(u32_t)p;
      mpqs_Adiv_start1[i]=(u16_t)cc;
    }
/*zeitB(7);*/
  return 1;
}


static void mpqs_sieve()
{
  unsigned char *sv=mpqs_sievearray, *fblog=mpqs_FB_log, *svend;
  u16_t *fb=mpqs_FB, *fbs=mpqs_FB_start;
  u16_t p, s1, s2, pb;
  long i;
  unsigned char lo;
  u32_t *ulsv, *ulsvend, mask;

  ulsv=(u32_t *)mpqs_sievearray;
  ulsvend=ulsv+mpqs_sievelen/4;
  if (mpqs_B&1) {
    mask=(u32_t)(mpqs_FB_log[0]); mask*=0x00040004;
  } else {
    mask=(u32_t)(mpqs_FB_log[0]); mask*=0x04000400;
  }
  while (ulsv<ulsvend) {
    *ulsv++=mask;
    *ulsv++=mask;
    *ulsv++=mask;
    *ulsv++=mask;
  }
/*printf("%lu %lu \n",(u32_t)(mpqs_sievearray[0]),(u32_t)(mpqs_sievearray[1]));*/
/*  mask*=0x00000101;
  if (mpqs_zero1<mpqs_sievelen) {
    s1=mpqs_zero1/4; if (s1<128) s1=0; else s1-=128;
    ulsv=(u32_t *)mpqs_sievearray; ulsv+=s1; ulsvend=ulsv+256;
    while (ulsv<ulsvend) { *ulsv++=mask; *ulsv++=mask; *ulsv++=mask; *ulsv++=mask; }
  }
  if (mpqs_zero2<mpqs_sievelen) {
    s1=mpqs_zero2/4; if (s1>(mpqs_sievelen/4)-128) s1=(mpqs_sievelen/4)-256; else s1-=128;
    ulsv=(u32_t *)mpqs_sievearray; ulsv+=s1; ulsvend=ulsv+256;
    while (ulsv<ulsvend) { *ulsv++=mask; *ulsv++=mask; *ulsv++=mask; *ulsv++=mask; }
  }*/
  pb=mpqs_sievelen>>2;
  i=mpqs_sievebegin; fb+=2*i; fbs+=2*i;

#ifdef ASM_MPQS_SIEVER
  asm_sieve();
#else
  for (;;) {
    p=*fb; if (p>pb) break;
    lo=fblog[i]; sv=mpqs_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++; 
    svend=sv+mpqs_sievelen-4*p;
    while (sv<svend) {
      sv[s1]+=lo; sv[s2]+=lo; sv+=p;
      sv[s1]+=lo; sv[s2]+=lo; sv+=p;
      sv[s1]+=lo; sv[s2]+=lo; sv+=p;
      sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    }
    svend+=(p+p);
    if (sv<svend) {
      sv[s1]+=lo; sv[s2]+=lo; sv+=p;
      sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    }
    svend+=p;
    if (sv<svend) {
      sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    }
    svend+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }
#endif
  sv=mpqs_sievearray;
#ifdef MPQS_ZEIT
zeita(6);
#endif
  for (i=0; i<mpqs_nAdiv_total; i++)
    if (!mpqs_Adiv_active[i]) {
      p=mpqs_Adiv_all[i]; lo=mpqs_Adiv_log[i];
      s1=mpqs_Adiv_start1[i];
      s2=mpqs_Adiv_start2[i];
      while (s1<mpqs_sievelen) { sv[s1]+=lo; s1+=p; }
      while (s2<mpqs_sievelen) { sv[s2]+=lo; s2+=p; }
    } else {
      p=mpqs_Adiv_all[i]; lo=mpqs_Adiv_log[i];
      s1=mpqs_Adiv_start1[i];
      while (s1<mpqs_sievelen) { sv[s1]+=lo; s1+=p; }
    }
#ifdef MPQS_ZEIT
zeitb(6);
#endif
}


static u16_t mpqs_evaluate()
{
  u32_t *ulsv, *ulsvend, h;
  unsigned char *sv;
  u16_t **rels, buffer[256];
  long i, nmaxsurv;

  mpqs_nsurvivors=0;
  nmaxsurv=MPQS_MAX_NRELS-mpqs_nrels;
  if (MPQS_MAX_NPRELS-mpqs_nprels<nmaxsurv)
    nmaxsurv=MPQS_MAX_NPRELS-mpqs_nprels;
  if (MPQS_MAX_NCRELS-mpqs_ncrels<nmaxsurv)
    nmaxsurv=MPQS_MAX_NCRELS-mpqs_ncrels;
  if (nmaxsurv<MPQS_MIN_RELBUFFER) return 0;
  if (nmaxsurv>255) nmaxsurv=255;  /* 1 weniger als buffer-Länge !!! */
  ulsv=(u32_t *)mpqs_sievearray;
  ulsvend=ulsv+(mpqs_sievelen>>2);
  rels=mpqs_rel_buffer;
#ifndef ASM_MPQS_EVAL
  while (ulsv<ulsvend) {
    h=*ulsv;
    if (h&0x80808080) {
      sv=(unsigned char *)ulsv;
      for (i=0; i<4; i++)
        if (*sv++&0x80)
          buffer[mpqs_nsurvivors++]=(u16_t)(sv-mpqs_sievearray-1);
      /* CJM: test added here. */
      if (mpqs_nsurvivors>nmaxsurv-4) {
        while (ulsv<ulsvend) *ulsv++=0;
        break;
      }
    }
    *ulsv++=0;
    h=*ulsv;
    if (h&0x80808080) {
      sv=(unsigned char *)ulsv;
      for (i=0; i<4; i++)
        if (*sv++&0x80)
          buffer[mpqs_nsurvivors++]=(u16_t)(sv-mpqs_sievearray-1);
      if (mpqs_nsurvivors>nmaxsurv-4) {
        while (ulsv<ulsvend) *ulsv++=0;
        break;
      }
    }
    *ulsv++=0;
  }
#else
  mpqs_nsurvivors=asm_evaluate(ulsv,ulsvend,buffer,nmaxsurv); /*printf("%ld ",mpqs_nsurvivors);*/
#ifdef MPQS_STAT
  stat_asm_eval++;
#endif
#endif

  for (i=0; i<mpqs_nsurvivors; i++) {
    mpqs_sievearray[buffer[i]]=(unsigned char)(i+1);
    rels[i][0]=buffer[i];
  }
  return mpqs_nsurvivors;
}


/* change hash_table size to 256*15 ??? */
static inline u16_t mpqs_hash_lp(u32_t lp)
{
  long j;
  u32_t lp1;
  u16_t lphash, nh;

  lphash=(u16_t)(lp&0x00ff); lphash>>=1;
  lp1=lp>>8;
  nh=mpqs_hash_table[lphash][0];
  if (nh>=15) { /*printf(".");*/ return MPQS_HASH_OVERFLOW; }
  for (j=1; j<=nh; j++)
    if (mpqs_hash_table[lphash][j]==lp1) break;
  if (j>nh) { /* new lp */
    mpqs_hash_table[lphash][j]=lp1;
    mpqs_hash_table[lphash][0]++;
    mpqs_nlp++;
    mpqs_hash_index[lphash][j]=mpqs_nprels;
    mpqs_lp[mpqs_nprels]=lp;
  }
  return mpqs_hash_index[lphash][j];
}


static inline int mpqs_hash_rel(i64_t axb)
{
  long j;
  u32_t hash, hash2, nh;

  hash=(u32_t)(axb&0xffffff);
  hash2=(u16_t)(hash>>8); hash&=0xff;
  nh=mpqs_rel_hash_table[hash][0];
  if (nh>=15) return 1;
  for (j=1; j<=nh; j++)
    if (mpqs_rel_hash_table[hash][j]==hash2) return 1;
  mpqs_rel_hash_table[hash][nh+1]=hash2;
  mpqs_rel_hash_table[hash][0]++;
  return 0;
}


static int mpqs_decompose()
{
  long i, j;
  u16_t *fb=mpqs_FB, *fbs=mpqs_FB_start;
  unsigned char *sv, *svend, i1, i2;
  u16_t ind, p, s1, s2, nr, lpi, ii;
  u16_t **rels;
  i64_t axb, *llp;
  u64_t qx, pr;
  u16_t minus;
  double dbl_qx, dx;
  i16_t x;
  u32_t lp, ulqx;
  u32_t inv, ls1, ls2;
  u32_t ax, ay, az, at;

#ifdef MPQS_ZEIT
zeitA(11);
#endif
  rels=mpqs_rel_buffer-1;
  for (i=1; i<=mpqs_nsurvivors; i++) rels[i][4]=0;
  i=mpqs_td_begin; fb+=2*i; fbs+=2*i;
  for (; i<mpqs_nFB; i++) {
    p=*fb; fb+=2; sv=mpqs_sievearray; svend=sv+mpqs_sievelen-2*p;
    s1=*fbs++; s2=*fbs++;
    while (sv<svend) {
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=mpqs_nFBk_1+i; rels[i1][4]++; }
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=mpqs_nFBk_1+i; rels[i2][4]++; }
      }
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=mpqs_nFBk_1+i; rels[i1][4]++; }
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=mpqs_nFBk_1+i; rels[i2][4]++; }
      }
    }
    svend+=p;
    while (sv<svend) {
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=mpqs_nFBk_1+i; rels[i1][4]++; }
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=mpqs_nFBk_1+i; rels[i2][4]++; }
      }
    }
    svend+=p;
    if (sv+s1<svend) {
      i1=sv[s1];
      if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=mpqs_nFBk_1+i; rels[i1][4]++; }
    }
    if (sv+s2<svend) {
      i2=sv[s2];
      if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=mpqs_nFBk_1+i; rels[i2][4]++; }
    }
  }
#ifdef MPQS_ZEIT
zeitB(11);
#endif
  for (i=1; i<=mpqs_nsurvivors; i++) {
    nr=rels[i][4];
    ind=rels[i][0];
    x=(i16_t)(ind)-mpqs_disp;

    axb=mpqs_A*(i64_t)x+mpqs_B;
    if (axb<0) axb=-axb;
    if (mpqs_hash_rel(axb)) goto next;

    qx=mpqs_A*(u64_t)x+mpqs_2B;
    qx*=(u64_t)x;
    qx+=(u64_t)mpqs_C;

/* das folgende besser durch ungefähre Berechnung der Nullstellen
   und Vergleiche ersetzen ?? */
    minus=0; dx=(double)x;
    dbl_qx=(double)mpqs_A*dx+(double)mpqs_2B;
    dbl_qx*=dx;
    dx=(double)mpqs_C; if (dx>0.) dx-=18446744073709551616.;
    dbl_qx+=dx;
    if (fabs(dbl_qx)>=18446744073709551616.) { /*printf(";");*/ goto next; }
    if (dbl_qx<0.) { minus=1; qx=-qx; }

#ifdef ASM_MPQS_TD
#ifdef MPQS_ZEIT
zeitA(12);
#endif
#ifdef MPQS_STAT
stat_asm_td++;
#endif
if (asm_td(rels[i],minus,&qx,&ulqx)) {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
  goto next;
}
/*
if (nr>=40) {
complain("%Lu %ld, %ld, ",qx,ulqx,nr);
}
*/
/*printf(",");*/
/*if (nr>27) {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
  goto next;
}*/
/*printf("%lu ",ulqx);*/
/*      for (j=0; j<mpqs_nFBk; j++) {
        p=mpqs_FBk[j];
        if (ulqx%(u32_t)p==0) {
          complain("%Lu %lu %lu ",qx,ulqx,p);
        }
      }*/
/*      for (j=0; j<mpqs_nFBk; j++) {
        p=mpqs_FBk[j];
        if (ulqx%(u32_t)p==0) {
          complain("%Lu %lu %lu ",qx,ulqx,p);
        }
      }*/

#else
    for (j=1; j<mpqs_td_begin; j++) {
      p=mpqs_FB[2*j];
      inv=mpqs_FB_inv[j];
      ls1=ind+(p-mpqs_FB_start[2*j]); ls2=ind+(p-mpqs_FB_start[2*j+1]);
      ls1*=inv; ls2*=inv;
      ls1&=0xffff0000; ls2&=0xffff0000;
      if (!ls1 || !ls2) {
        rels[i][5+nr]=mpqs_nFBk_1+j; nr++; /*printf("%d.",p);*/
      }
    }

    pr=1;
    for (j=0; j<nr; j++) pr*=(u64_t)(mpqs_FB[2*(rels[i][5+j]-mpqs_nFBk_1)]);
    rels[i][4]=nr;
    if (minus) {
      rels[i][5+rels[i][4]]=0; rels[i][4]++;
    }
    while (!(qx&1)) {
      if (rels[i][4]<MPQS_TD_MAX_NDIV) {
        rels[i][5+rels[i][4]]=mpqs_nFBk_1; rels[i][4]++;
        qx>>=1;
      } else {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
        goto next;
      }
    }
    if (pr<4294967296ULL) {
      if (qx>(pr<<32)) {
/*        printf(",");*/
#ifdef MPQS_ZEIT
zeitB(12);
#endif
        goto next;
      }
    }
    ax=(u32_t)qx; ay=(u32_t)pr; /*az=0;*/
    inv=(u32_t)mpqs_256_inv_table[(ay&0xff)>>1];
    at=ay*inv; at&=0x0000ff00;
    at*=inv; inv-=at;
    at=ay*inv; at&=0xffff0000;
    at*=inv; inv-=at;
/*    if (ay*inv!=1) complain("dec.neu8");*/
    az=ax*inv;
/*if (qx%pr) complain("a");
if (qx/pr!=az) complain("b");*/
ulqx=az;

    for (j=0; j<nr; j++) {
      ii=rels[i][5+j];
      p=mpqs_FB[2*(ii-mpqs_nFBk_1)];
      while (ulqx%(u32_t)p==0) {
        if (rels[i][4]<MPQS_TD_MAX_NDIV) {
          rels[i][5+rels[i][4]]=ii; rels[i][4]++;
          ulqx/=(u32_t)p;
        } else {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
          goto next;
        }
      }
    }


    if (ulqx>1)
      for (j=0; j<mpqs_nFBk; j++) {
        p=mpqs_FBk[j];
        while (ulqx%(u32_t)p==0) {
          if (rels[i][4]<MPQS_TD_MAX_NDIV) {
            rels[i][5+rels[i][4]]=1+j; rels[i][4]++;
            ulqx/=(u32_t)p;
          } else {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
            goto next;
          }
        }
        if (ulqx==1) break;
      }
    if (ulqx>1)
      for (j=0; j<mpqs_nAdiv_total; j++) {
        p=mpqs_Adiv_all[j];
        while (ulqx%(u32_t)p==0) {
          if (rels[i][4]<MPQS_TD_MAX_NDIV) {
            rels[i][5+rels[i][4]]=1+mpqs_nFB+mpqs_nFBk+j; rels[i][4]++;
            ulqx/=(u32_t)p;
          } else {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
            goto next;
          }
        }
        if (ulqx==1) break;
      }

#endif

#if 0
    for (j=0; j<nr; j++) {
      u32_t jj;

      ii=rels[i][5+j]; jj=ii-mpqs_nFBk_1;
      p=mpqs_FB[2*jj];
      inv=mpqs_FB_inv[jj];
      while (1) {
	u32_t rr;
        rr=ulqx*inv;
        if (((u64_t)rr*(u64_t)p) & 0xffffffff00000000ULL) {
/*          if (ulqx%(u32_t)p==0) complain("%lu %lu %lu %lu y",rr,p,ulqx,inv);*/
          break;
        }
/*        if (ulqx%(u32_t)p) complain("z");*/
        if (rels[i][4]<MPQS_TD_MAX_NDIV) {
          rels[i][5+rels[i][4]]=ii; rels[i][4]++;
 /*         ulqx/=(u32_t)p;*/
          ulqx=rr;
        } else {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
          goto next;
        }
      }
    }
#endif

    if (rels[i][4]<MPQS_TD_MAX_NDIV-mpqs_nAdiv) {
      nr=rels[i][4]; ind=0;
      for (j=0; j<mpqs_nAdiv_total; j++)
        if (mpqs_Adiv_active[j]) {
          rels[i][5+nr+ind]=1+mpqs_nFB+mpqs_nFBk+j;
          ind++;
        }
      rels[i][4]+=mpqs_nAdiv;
    } else {
#ifdef MPQS_ZEIT
 zeitB(12);
#endif
      goto next;
    }
#ifdef MPQS_ZEIT
zeitB(12);
#endif
    if (ulqx==1) {
      llp=(i64_t *)(mpqs_relations[mpqs_nrels]);
      *llp=axb;
      nr=rels[i][4];
      for (j=0; j<=nr; j++)
        mpqs_relations[mpqs_nrels][j+4]=rels[i][j+4];
      mpqs_nrels++;
      goto next;
    }
    if (ulqx<1048576) { /*printf("%d ",x);*/  /* !!! */
      if (ulqx<=mpqs_pmax) {
        /*printf("%d %lu %u ",x,qx,ulqx);*/ /* removes garbage output from 64bit siever */
return -1; /* CJM, 11/30/04. */
        complain("dec.2\n");
      }
      lp=ulqx;
      lpi=mpqs_hash_lp(lp);
      if (lpi==MPQS_HASH_OVERFLOW) goto next;
      if (lpi<mpqs_nprels) {
        if (mpqs_lp[lpi]!=lp) Schlendrian("lp");
        llp=(i64_t *)(mpqs_comb_rels[mpqs_ncrels]);
        *llp=axb;
        llp=(i64_t *)(mpqs_part_rels[lpi]);
        axb=*llp;
        llp=(i64_t *)(mpqs_comb_rels[mpqs_ncrels]); llp++;
        *llp=axb;
        nr=rels[i][4];
        for (j=0; j<nr; j++)
          mpqs_comb_rels[mpqs_ncrels][j+9]=rels[i][j+5];
        for (j=0; j<mpqs_part_rels[lpi][4]; j++)
          mpqs_comb_rels[mpqs_ncrels][j+9+nr]=mpqs_part_rels[lpi][j+5];
        mpqs_comb_rels[mpqs_ncrels][8]=nr+256*mpqs_part_rels[lpi][4];
        mpqs_ncrels++;
      } else {
        llp=(i64_t *)(mpqs_part_rels[mpqs_nprels]);
        *llp=axb;
        nr=rels[i][4];
        for (j=0; j<=nr; j++)
          mpqs_part_rels[mpqs_nprels][j+4]=rels[i][j+4];
        mpqs_nprels++;
      }
    } /*else printf("<");*/
next:
    ; /* Null statement to hush compiler warning. */
  }
  if (mpqs_nsp<mpqs_nrels+mpqs_ncrels)
    mpqs_excess=mpqs_nrels+mpqs_ncrels-mpqs_nsp;
  else mpqs_excess=0;
  return 0;
}


/*
order of the primes:
-1, primes in mpqs_FBk, primes in mpqs_FB, primes in mpqs_Adiv_total
*/

/*
static u32_t mpqs_gauss_mask0[]  =
{
  0x80000000UL, 0x40000000UL, 0x20000000UL, 0x10000000UL,
  0x08000000UL, 0x04000000UL, 0x02000000UL, 0x01000000UL,
  0x00800000UL, 0x00400000UL, 0x00200000UL, 0x00100000UL,
  0x00080000UL, 0x00040000UL, 0x00020000UL, 0x00010000UL,
  0x00008000UL, 0x00004000UL, 0x00002000UL, 0x00001000UL,
  0x00000800UL, 0x00000400UL, 0x00000200UL, 0x00000100UL,
  0x00000080UL, 0x00000040UL, 0x00000020UL, 0x00000010UL,
  0x00000008UL, 0x00000004UL, 0x00000002UL, 0x00000001UL
};
*/
static u32_t mpqs_gauss_mask[]  =
{
  0x00000001UL, 0x00000002UL, 0x00000004UL, 0x00000008UL,
  0x00000010UL, 0x00000020UL, 0x00000040UL, 0x00000080UL,
  0x00000100UL, 0x00000200UL, 0x00000400UL, 0x00000800UL,
  0x00001000UL, 0x00002000UL, 0x00004000UL, 0x00008000UL,
  0x00010000UL, 0x00020000UL, 0x00040000UL, 0x00080000UL,
  0x00100000UL, 0x00200000UL, 0x00400000UL, 0x00800000UL,
  0x01000000UL, 0x02000000UL, 0x04000000UL, 0x08000000UL,
  0x10000000UL, 0x20000000UL, 0x40000000UL, 0x80000000UL
};

#define PRUNING

static void mpqs_matrix_init()
{
  long i, j;
  u16_t nr, ii, pi;
  u32_t mask;
#ifdef PRUNING
  u16_t mpqs_gauss_wt[300];  /* verbessern !!!! */

  for (i=0; i<mpqs_nsp; i++) mpqs_gauss_wt[i]=0;
#endif

  mpqs_gauss_n=mpqs_nrels+mpqs_ncrels;
  if (mpqs_gauss_n>=MPQS_GAUSS_MAX) return;
  mpqs_gauss_m=mpqs_nsp;
  if (mpqs_gauss_m>=mpqs_gauss_n) Schlendrian("gauss: no excess\n");
/* build matrix */
  mpqs_gauss_n32=(mpqs_gauss_n+31)/32;
  for (i=1; i<mpqs_gauss_m; i++)
    mpqs_gauss_row[i]=mpqs_gauss_row[i-1]+mpqs_gauss_n32;
  for (i=0; i<mpqs_gauss_m; i++)
    for (j=0; j<mpqs_gauss_n32; j++)
      mpqs_gauss_row[i][j]=0;
  for (i=0; i<mpqs_nrels; i++) {
    nr=mpqs_relations[i][4]; mask=mpqs_gauss_mask[i%32];
    for (j=0; j<nr; j++) {
      pi=mpqs_relations[i][5+j];
      mpqs_gauss_row[pi][i/32]^=mask;
#ifdef PRUNING
      mpqs_gauss_wt[pi]++;
#endif
    }
  }
  for (i=0; i<mpqs_ncrels; i++) {
    nr=mpqs_comb_rels[i][8]; nr=(nr&255)+(nr>>8);
    ii=i+mpqs_nrels; mask=mpqs_gauss_mask[ii%32];
    for (j=0; j<nr; j++) {
      pi=mpqs_comb_rels[i][9+j];
      mpqs_gauss_row[pi][ii/32]^=mask;
#ifdef PRUNING
      mpqs_gauss_wt[pi]++;
#endif
    }
  }
  for (i=0; i<mpqs_gauss_m; i++) mpqs_gauss_c[i]=-1;
  mpqs_gauss_k=0;
  for (i=0; i<mpqs_gauss_n; i++) mpqs_gauss_d[i]=-1;

#ifdef PRUNING
/* simple pruning */
{
  u32_t k, k32;

#if 1
/*  u32_t anz[5];
  for (j=0; j<5; j++) anz[j]=0;
  for (j=0; j<mpqs_gauss_m; j++)
    if (mpqs_gauss_wt[j]<5) anz[mpqs_gauss_wt[j]]++;
  printf("| %lu: ",mpqs_gauss_m);
  for (j=0; j<5; j++) printf("%lu ",anz[j]); printf("| ");
  for (j=0; j<5; j++) stat_gauss_anz[j]+=anz[j];*/

  for (j=0; j<mpqs_gauss_m; j++)
    if (mpqs_gauss_wt[j]==1) { /*printf("!");*/
      k32=0;
      while (mpqs_gauss_row[j][k32]==0) k32++;
      for (k=0; k<32; k++)
        if (mpqs_gauss_row[j][k32]&mpqs_gauss_mask[k]) break;
      k+=(k32<<5);
      if (k>=mpqs_gauss_n) Schlendrian("pruning\n");

      mpqs_gauss_d[k]=j; mpqs_gauss_c[j]=k;
/* No need to clear the k-th column since it will not appear in a solution. */
    }
#else    /* unsinn!!! */
  for (j=0; j<mpqs_gauss_m; j++)
    if (mpqs_gauss_wt[j]==1) { /*printf("!");*/
      mpqs_gauss_m--;
      for (k=0; k<mpqs_gauss_n32; k++)
        mpqs_gauss_row[j][k]=mpqs_gauss_row[mpqs_gauss_m][k];
      mpqs_gauss_wt[j]=mpqs_gauss_wt[mpqs_gauss_m]; break;
    }
#endif
}
#endif
}


static int mpqs_matrix()
{
  long i, j, l, t, k32;
  u32_t mask;
 
#ifdef MPQS_ZEIT
zeita(8);
#endif
/* solve matrix */
  while (mpqs_gauss_k<mpqs_gauss_n) {
#ifdef PRUNING
    if (mpqs_gauss_d[mpqs_gauss_k]>=0) { mpqs_gauss_k++; continue; }
#endif
    mask=mpqs_gauss_mask[mpqs_gauss_k%32]; k32=mpqs_gauss_k/32;
    j=0;
    while ((j<mpqs_gauss_m) && (mpqs_gauss_c[j]>=0 ||
       ((mpqs_gauss_row[j][k32] & mask)==0))) j++;
    if (j==mpqs_gauss_m) { mpqs_gauss_d[mpqs_gauss_k]=-1; break; }
    mpqs_gauss_d[mpqs_gauss_k]=j; mpqs_gauss_c[j]=mpqs_gauss_k;
    for (t=0; t<j; t++)
      if (mpqs_gauss_row[t][k32] & mask)
        for (l=k32; l<mpqs_gauss_n32; l++)
          mpqs_gauss_row[t][l]^=mpqs_gauss_row[j][l];
    for (t=j+1; t<mpqs_gauss_m; t++)
      if (mpqs_gauss_row[t][k32] & mask)
        for (l=k32; l<mpqs_gauss_n32; l++)
          mpqs_gauss_row[t][l]^=mpqs_gauss_row[j][l];
    mpqs_gauss_k++;
  }
  if (mpqs_gauss_k>=mpqs_gauss_n) {
#ifdef MPQS_ZEIT
zeitb(8);
#endif
    return 0;
  }
  for (i=0; i<mpqs_gauss_n; i++) mpqs_sol[i]=0;
  if (mpqs_gauss_d[mpqs_gauss_k]!=-1) Schlendrian("gauss1\n");
  for (i=0; i<mpqs_gauss_k; i++)
    if (mpqs_gauss_d[i]!=-1)
      if (mpqs_gauss_row[mpqs_gauss_d[i]][k32] & mask)
        mpqs_sol[i]=1;
  mpqs_sol[mpqs_gauss_k]=1;
  mpqs_gauss_k++;
#ifdef MPQS_ZEIT
zeitb(8);
#endif
  return 1;
}


static int mpqs_final()
{
  long j, k;
  u16_t nr, nr1, nr2;
  u32_t p, prod1, prod2;
  u64_t *up;
  int split;

/*test_rels(); test_combs();*/
  if (!mpqs_nfactors) {
    mpqs_nfactors=1;
    mpz_set(mpqs_factors[0],mpqs_N);
    mpqs_f_status[0]=0;
  }
  while (mpqs_matrix()) {
#ifdef MPQS_STAT
stat_mpqs_ntrials++;
#endif
    for (j=0; j<mpqs_nsp; j++) mpqs_exp[j]=0;
    for (j=0; j<mpqs_nrels; j++)
      if (mpqs_sol[j]) {
        nr=mpqs_relations[j][4];
        if (j&1)     /* 'random' distribution of full relations */
          for (k=0; k<nr; k++) mpqs_exp[mpqs_relations[j][5+k]]--;
        else
          for (k=0; k<nr; k++) mpqs_exp[mpqs_relations[j][5+k]]++;
      }
    for (j=0; j<mpqs_ncrels; j++)
      if (mpqs_sol[mpqs_nrels+j]) {
        nr=mpqs_comb_rels[j][8];
        nr1=nr&255; nr2=nr>>8;
        for (k=0; k<nr1; k++) mpqs_exp[mpqs_comb_rels[j][9+k]]++;
        for (k=0; k<nr2; k++) mpqs_exp[mpqs_comb_rels[j][9+nr1+k]]--;
      }
    for (j=0; j<mpqs_nsp; j++)
      if (mpqs_exp[j]&1) Schlendrian("final.odd\n");
      else mpqs_exp[j]>>=1;



    mpz_set_ui(mpqs_sq1,1); prod1=1;
    mpz_set_ui(mpqs_sq2,1); prod2=1;
    for (j=1; j<mpqs_nsp; j++)
      if (mpqs_exp[j]) {
        if (j<1+mpqs_nFBk) {
          p=(u32_t)(mpqs_FBk[j-1]);
        } else if (j<1+mpqs_nFB+mpqs_nFBk) {
          k=j-mpqs_nFBk-1;
          p=(u32_t)(mpqs_FB[2*k]);
        } else {
          k=j-mpqs_nFB-mpqs_nFBk-1;
          p=(u32_t)(mpqs_Adiv_all[k]);
        }
        if (mpqs_exp[j]>0) {
          for (k=0; k<mpqs_exp[j]; k++) {
            prod1*=p;
            if (prod1&0xfff00000) {  /* p<4096 */
              mpz_mul_ui(mpqs_sq1,mpqs_sq1,prod1);
              mpz_mod(mpqs_sq1,mpqs_sq1,mpqs_N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              prod1=1;
            }
          }
        } else {
          for (k=0; k<-mpqs_exp[j]; k++) {
            prod2*=p;
            if (prod2&0xfff00000) {  /* p<4096 */
              mpz_mul_ui(mpqs_sq2,mpqs_sq2,prod2);
              mpz_mod(mpqs_sq2,mpqs_sq2,mpqs_N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              prod2=1;
            }
          }
        }
      }
    if (prod1>1) {
      mpz_mul_ui(mpqs_sq1,mpqs_sq1,prod1);
      mpz_mod(mpqs_sq1,mpqs_sq1,mpqs_N);
    }
    if (prod2>1) {
      mpz_mul_ui(mpqs_sq2,mpqs_sq2,prod2);
      mpz_mod(mpqs_sq2,mpqs_sq2,mpqs_N);
    }

    for (j=0; j<mpqs_nrels; j++)
      if (mpqs_sol[j]) {
        up=(u64_t *)(mpqs_relations[j]); /*printf("Q: %Lu  ",*up);*/
        mpz_set_ull(mpqs_gmp_help,*up);
        if (j&1) {
          mpz_mul(mpqs_sq1,mpqs_sq1,mpqs_gmp_help);
          mpz_mod(mpqs_sq1,mpqs_sq1,mpqs_N);
        } else {
          mpz_mul(mpqs_sq2,mpqs_sq2,mpqs_gmp_help);
          mpz_mod(mpqs_sq2,mpqs_sq2,mpqs_N);
        }
      }

    for (j=0; j<mpqs_ncrels; j++)
      if (mpqs_sol[mpqs_nrels+j]) {
        up=(u64_t *)(mpqs_comb_rels[j]); /*printf("Q: %Lu  ",*up);*/
        mpz_set_ull(mpqs_gmp_help,*up);
        mpz_mul(mpqs_sq2,mpqs_sq2,mpqs_gmp_help);
        mpz_mod(mpqs_sq2,mpqs_sq2,mpqs_N);
        up++;
        mpz_set_ull(mpqs_gmp_help,*up);
        mpz_mul(mpqs_sq1,mpqs_sq1,mpqs_gmp_help);
        mpz_mod(mpqs_sq1,mpqs_sq1,mpqs_N);
      }

    mpz_add(mpqs_gmp_help,mpqs_sq1,mpqs_sq2);
    mpz_gcd(mpqs_gmp_help,mpqs_gmp_help,mpqs_N);
    if (mpz_cmp_ui(mpqs_gmp_help,1)==0) continue;
    if (mpz_cmp(mpqs_gmp_help,mpqs_N)==0) continue;

    for (j=0; j<mpqs_nfactors; j++)
      if (mpqs_f_status[j]==0) {
        mpz_gcd(mpqs_sq1,mpqs_gmp_help,mpqs_factors[j]);
        if (mpz_cmp_ui(mpqs_sq1,1)==0) continue;
        if (mpz_cmp(mpqs_sq1,mpqs_factors[j])==0) continue;
        mpz_set(mpqs_factors[mpqs_nfactors],mpqs_sq1);
        mpz_divexact(mpqs_factors[j],mpqs_factors[j],mpqs_sq1);
        mpqs_f_status[j]=psp(mpqs_factors[j]);
        mpqs_f_status[mpqs_nfactors]=psp(mpqs_factors[mpqs_nfactors]);
        mpqs_nfactors++;
        break;
      }

    split=1;
    for (j=0; j<mpqs_nfactors; j++) if (mpqs_f_status[j]==0) split=0;
    if (split) return 1;
    if (mpqs_nfactors>=MPQS_MAX_FACTORS)
      Schlendrian("final: too many factors\n");
  }
  return 0;
}



static long mpqs_factor0(mpz_t N, long max_bits, mpz_t **factors, long retry)
{
  long nbits, i, err;
  int ev;

#ifdef MPQS_ZEIT
zeita(9);
#endif
#ifdef MPQS_STAT
  stat_mpqs_nsieves=0; stat_mpqs_nsurvivors=0; stat_mpqs_ntrials=0; stat_mpqs_ndiv=0;
#endif
  if (!mpqs_isinit) mpqs_init();
#ifdef MPQS_ZEIT
zeitb(9);
#endif
  nbits=mpz_sizeinbase(N,2);
  if (nbits>96) {
    fprintf(stderr,"warning: mpqs called with >96 Bit\n");
    return -2;
  }
#ifdef MPQS_ZEIT
zeita(1);
#endif
  mpz_set(mpqs_N,N);
  mpqs_choose_multiplier();
  mpqs_choose_parameter(retry);
  mpqs_generate_FB();
  if ((err=mpqs_SI_init())<1) {
    printf("%ld ",err);
    fprintf(stderr,"warning: mpqs: error in self-initialization ");
    mpz_out_str(stderr,10,N); printf("\n");
#ifdef MPQS_ZEIT
zeitb(1);
#endif
    return -3;
  }
#ifdef MPQS_ZEIT
zeitb(1);
#endif
  while (1) { /*printf("%ld,%ld; ",stat_mpqs_nsieves,mpqs_nrels);*/
#ifdef MPQS_ZEIT
zeita(2);
#endif
    if (!mpqs_next_pol()) {     
    /* Silence the frequent, but harmless warning below */
    /* if (retry) {
        fprintf(stderr,"warning: not enough polynomials in mpqs with ");
        mpz_out_str(stderr,10,N); fprintf(stderr,"\n");
      }
    */
#ifdef MPQS_ZEIT
zeitb(2);
#endif
      return -4;
    }
#ifdef MPQS_ZEIT
zeitb(2); zeita(3);
#endif
    mpqs_sieve();
#ifdef MPQS_STAT
stat_mpqs_nsieves++;
#endif
#ifdef MPQS_ZEIT
zeitb(3);
zeita(4);
#endif
    ev=mpqs_evaluate();
#ifdef MPQS_STAT
stat_mpqs_nsurvivors+=ev; /*printf("%d %Ld  ",ev,mpqs_C); */
#endif
#ifdef MPQS_ZEIT
zeitb(4);
#endif
    if (ev) {
#ifdef MPQS_ZEIT
zeita(5);
#endif
    { int mpqs_dec_res; /* CJM, 11/30/04. */

      if ((mpqs_dec_res = mpqs_decompose()))
        return -1;
    }
#ifdef MPQS_ZEIT
zeitb(5);
#endif
      if (mpqs_nsurvivors && (mpqs_excess>MPQS_MIN_EXCESS)) {
#ifdef MPQS_ZEIT
zeita(8);
#endif
        mpqs_matrix_init();
#ifdef MPQS_ZEIT
zeitb(8);
zeita(10);
#endif
        if (mpqs_final()) {
#ifdef MPQS_ZEIT
zeitb(10);
#endif
#ifdef MPQS_STAT
printf("Stat: %lu %lu %lu %lu %lu %lu %lu \n",stat_mpqs_nsieves,stat_mpqs_nsurvivors,stat_mpqs_ntrials,stat_mpqs_ndiv,mpqs_nrels,mpqs_nprels,mpqs_ncrels);
#endif
          for (i=0; i<mpqs_nfactors; i++)
            if (mpz_sizeinbase(mpqs_factors[i],2)>max_bits) return 0;
          *factors=mpqs_factors;
          return mpqs_nfactors;
        }
#ifdef MPQS_ZEIT
zeitb(10);
#endif
      }
    }
    if (mpqs_nrels>MPQS_MAX_NRELS-MPQS_MIN_RELBUFFER) {
      /* Silence the frequent, but harmless warning */
      /* fprintf(stderr,"warning: too many relations in mpqs\n"); */
      return -5;
    }
  }
  return -1;
}


long mpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
  long err;

  err=mpqs_factor0(N,max_bits,factors,0);
  if (err==-4) err=mpqs_factor0(N,max_bits,factors,1);
  if (err==-4) err=mpqs_factor0(N,max_bits,factors,3);
  return err;
}

