/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/
/* Written by T. Kleinjung, with some modifications by J. Franke. */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include "gmp.h"
#include "asm/siever-config.h"
#include "if.h"
#include "gmp-aux.h"
#include "asm/mpqs-config.h"


#define uchar unsigned char
#ifdef TOTAL_STAT
#define MPQS_STAT
ulong td_total=0,tds_total=0,total_nfbp=0;
ulong total_ulqx=0;
double sieve_events=0,total_sieve_length=0;
#endif

#ifdef MPQS_STAT
extern u64_t stat_asm_eval, stat_asm_td;
extern u64_t stat_td_cand,stat_td_surv;
extern u64_t stat_ff,stat_pf,stat_comb;
extern u32_t stat_asm_div, stat_final_mulmod;
extern u32_t stat_counter0, stat_retry;
extern u32_t stat_size[12];
#endif


#define L1BITS               14

#define MPQS_MAX_FBSIZE      512
#define MPQS_MIN_EXCESS      10
#define MPQS_MAX_ADIV_ALL    10
#define MPQS_SIEVELENMAX     (1<<L1BITS)
#define MPQS_REL_ENTRIES     32
#define MPQS_TD_MAX_NDIV     (MPQS_REL_ENTRIES-5)
#define MPQS_MAX_NRELS       384  /* ??? */
#define MPQS_MAX_NPRELS      2048  /* ??? */
#define MPQS_MAX_NCRELS      256  /* ??? */
#define MPQS_MAX_NRELS_BUF   256  /* ??? */
#define MPQS_MIN_RELBUFFER   16
#define MPQS_GAUSS_MAX       512
#define MPQS_MAX_FACTORS     5
#define MPQS_MAX_NPRIMES     1024
#define MPQS_HASH_OVERFLOW   65535
#define MPQS_MULT_NTESTPR    25       /* prime(MPQS_MULT_NTESTPR)<2^8 */
#define MPQS_MULT_NCAND      8        /* largest candidate<2^8 */
#define MPQS_FB_MAXPRIME     4096     /* used in nextpol, sieve and final */
#define MPQS_SI_NVALUES      4
#define MPQS_MAX_TINYPROD    128

/* common with asm functions */
ushort mpqs_nFBk_1;
ushort mpqs_td_begin, mpqs_sievebegin, mpqs_sievebegink;
ushort mpqs_FB_inv_info[4*MPQS_MAX_NPRIMES];
ushort mpqs_nFBk, mpqs_FBk[3];  /* 3*5*7*11 too big */
uchar mpqs_FB_log[MPQS_MAX_FBSIZE], mpqs_FBk_log[3];
ushort mpqs_nFB, mpqs_nAdiv_total;
ushort mpqs_Adiv_all[MPQS_MAX_ADIV_ALL];
u32_t mpqs_FBk_inv[3], mpqs_FB_A_inv[MPQS_MAX_ADIV_ALL];

ushort mpqs_Adiv_disp[MPQS_MAX_ADIV_ALL];
ushort mpqs_FB_mm_inv[MPQS_MAX_FBSIZE], mpqs_FB_disp[MPQS_MAX_FBSIZE];
ushort mpqs_FB_np_p[2*MPQS_MAX_FBSIZE];
ushort mpqs_FB_np_px[2*MPQS_MAX_FBSIZE];
/* end */

static int mpqs_isinit=0;
static mpz_t mpqs_gmp_help, mpqs_dummy, mpqs_help;
static mpz_t mpqs_N, mpqs_kN, mpqs_factors[MPQS_MAX_FACTORS];
static uchar mpqs_f_status[MPQS_MAX_FACTORS];
static u32_t mpqs_multiplier, mpqs_nfactors;
static double mpqs_complexity;
static ushort mpqs_prime_table[MPQS_MAX_NPRIMES];
static ushort mpqs_prime_sqrthelp[MPQS_MAX_NPRIMES];
static ushort mpqs_nAdiv, mpqs_accept;
static uchar mpqs_Adiv_log[MPQS_MAX_ADIV_ALL];;
static ushort mpqs_pmax;
static ushort mpqs_Adiv[MPQS_MAX_ADIV_ALL];
static ushort mpqs_Adiv_sqrt[MPQS_MAX_ADIV_ALL], mpqs_Adiv_all_sqrt[MPQS_MAX_ADIV_ALL];
static ushort mpqs_Adiv_start1[MPQS_MAX_ADIV_ALL], mpqs_Adiv_start2[MPQS_MAX_ADIV_ALL];
static uchar mpqs_Adiv_active[MPQS_MAX_ADIV_ALL];
static ushort mpqs_A_mask[(1<<MPQS_MAX_ADIV_ALL)];
static ushort mpqs_nA;
static short mpqs_A_index, mpqs_B_index;
static ushort mpqs_SI_inv_table[MPQS_MAX_ADIV_ALL][MPQS_MAX_FBSIZE];
static ushort mpqs_Adiv_SI_inv_table[MPQS_MAX_ADIV_ALL][MPQS_MAX_ADIV_ALL];
static ushort mpqs_SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_FBSIZE], mpqs_SIk_add[MPQS_MAX_ADIV_ALL][3];
static ushort mpqs_Adiv_SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_ADIV_ALL];
static long long mpqs_A, mpqs_B, mpqs_C, mpqs_2B, mpqs_Bi[MPQS_MAX_ADIV_ALL];
static u32_t mpqs_disp;
static ushort mpqs_nsurvivors;
static u32_t mpqs_nrels, mpqs_nlp, mpqs_nsp, mpqs_excess;
static uchar mpqs_logbound;
static ushort mpqs_nprels, mpqs_ncrels;
static ushort **mpqs_relations, **mpqs_part_rels, **mpqs_comb_rels;
static ushort **mpqs_rel_buffer;
static ushort mpqs_rel_hash_table[256][16];
static ushort mpqs_hash_table[128][16];
static u32_t mpqs_hash_index[128][16];
static u32_t mpqs_lp[128*15];
static uchar *mpqs_rel_active; /* auf bitarray umstellen !!! */
static ushort mpqs_nrels_active, mpqs_nlp_active;
static ushort mpqs_lp_count[128*16+MPQS_MAX_FBSIZE+3+MPQS_MAX_ADIV_ALL];
u32_t **mpqs_gauss_row, *mpqs_gauss_mat;
static u64_t **mpqs_gauss_row64;
static uchar mpqs_sol[MPQS_GAUSS_MAX];
short mpqs_gauss_c[MPQS_GAUSS_MAX], mpqs_gauss_d[MPQS_GAUSS_MAX];
short mpqs_gauss_col[MPQS_GAUSS_MAX];
u32_t mpqs_gauss_m, mpqs_gauss_n, mpqs_gauss_n32, mpqs_gauss_k, mpqs_gauss_j;
u32_t mpqs_gauss_n64;
static short mpqs_exp[128*16+MPQS_MAX_FBSIZE+3+MPQS_MAX_ADIV_ALL];
static mpz_t mpqs_sq1, mpqs_sq2;
static ushort mpqs_zero1, mpqs_zero2;
static double mpqs_kN_dbl, mpqs_maxval_dbl;
static u64_t mpqs_kN_64, mpqs_A_inv_64;

ushort mpqs_FB[2*MPQS_MAX_FBSIZE], mpqs_FB0[MPQS_MAX_FBSIZE+2*256];
ushort mpqs_FB_start[2*MPQS_MAX_FBSIZE], mpqs_FBk_start[3];
u32_t mpqs_sievelen;
u32_t mpqs_FB_inv[MPQS_MAX_FBSIZE];
uchar *mpqs_sievearray, *mpqs_tinyarray;
u32_t mpqs_tiny_prod, mpqs_ntiny;
ushort mpqs_sinit_tab[8*MPQS_SI_NVALUES];

uchar mpqs_256_inv_table[128]={
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


uchar mpqs_A_table[256]={
0, 1, 2, 4, 8, 16, 32, 64,
128, 3, 5, 6, 9, 10, 12, 17,
18, 20, 24, 33, 34, 36, 40, 48,
65, 66, 68, 72, 80, 96, 129, 130,
132, 136, 144, 160, 192, 7, 11, 13,
14, 19, 21, 22, 25, 26, 28, 35,
37, 38, 41, 42, 44, 49, 50, 52,
56, 67, 69, 70, 73, 74, 76, 81,
82, 84, 88, 97, 98, 100, 104, 112,
131, 133, 134, 137, 138, 140, 145, 146,
148, 152, 161, 162, 164, 168, 176, 193,
194, 196, 200, 208, 224, 15, 23, 27,
29, 30, 39, 43, 45, 46, 51, 53,
54, 57, 58, 60, 71, 75, 77, 78,
83, 85, 86, 89, 90, 92, 99, 101,
102, 105, 106, 108, 113, 114, 116, 120,
135, 139, 141, 142, 147, 149, 150, 153,
154, 156, 163, 165, 166, 169, 170, 172,
177, 178, 180, 184, 195, 197, 198, 201,
202, 204, 209, 210, 212, 216, 225, 226,
228, 232, 240, 31, 47, 55, 59, 61,
62, 79, 87, 91, 93, 94, 103, 107,
109, 110, 115, 117, 118, 121, 122, 124,
143, 151, 155, 157, 158, 167, 171, 173,
174, 179, 181, 182, 185, 186, 188, 199,
203, 205, 206, 211, 213, 214, 217, 218,
220, 227, 229, 230, 233, 234, 236, 241,
242, 244, 248, 63, 95, 111, 119, 123,
125, 126, 159, 175, 183, 187, 189, 190,
207, 215, 219, 221, 222, 231, 235, 237,
238, 243, 245, 246, 249, 250, 252, 127,
191, 223, 239, 247, 251, 253, 254, 255
};

u32_t mpqs_A_table_n[10]={
0, 1, 9, 37, 93, 163, 219, 247, 255, 256
};

#ifdef MPQS_STAT
u32_t stat_mpqs_nsieves, stat_mpqs_nsurvivors, stat_mpqs_ntrials, stat_mpqs_ndiv;
#endif

static u64_t mpqs_inv_64(u64_t a)
{
  u64_t inv, h;

  inv=(u64_t)mpqs_256_inv_table[(a&0xff)>>1];
  h=a*inv; h&=0xff00ULL;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffff0000ULL;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffffffff00000000ULL;
  h*=inv; inv-=h;
  /*if (inv*a!=1ULL) Schlendrian("mpqs_inv_64");*/
  return inv;
}


static uchar mpqs_jacobi_tab2[8]={0,0,0,1,0,1,0,0};

static int mpqs_jacobi(ushort a, ushort b) /* always gcd(a,b)=1 and 0<a<b */
{
  int e, l, m, n, r;

//zeitA(26);
  m=(int)a; n=(int)b;
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
//zeitB(26);
  return e;
}


static ushort mpqs_invert(ushort a, ushort p)
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


#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p32; h32=h32>>16


static ushort mpqs_powmod0(ushort a, ushort e, ushort p)
{
  ushort ex=e;
  u32_t aa, res;

  if (!ex) return 1;
  aa=(u32_t)a;
  res=1;
  while (1) {
    if (ex&1) { res*=aa; res%=(u32_t)p; }
    ex>>=1;
    if (!ex) return (ushort)res;
    aa*=aa; aa%=(u32_t)p;
  }
  return 0; /* never reached */
}

static ushort mpqs_powmod(ushort a, ushort e, ushort p)
{
  ushort ex=e, rrr;
  u32_t aa, res, p32, inv, h, h1;

//rrr=mpqs_powmod0(a,e,p);
  if (!ex) return 1;
  aa=(u32_t)a;
  p32=(u32_t)p;

  inv=(u32_t)mpqs_256_inv_table[(p&0xff)>>1];
  h=p32*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=p32*inv; h&=0xffff0000;
  h*=inv; inv-=h;
  if (inv*p32!=1) Schlendrian("mpqs_powmod_inv");
  inv=-inv;

  res=0x10000%p32;
  aa<<=16; aa%=p32;
  while (1) {
    if (ex&1) {
      res*=aa;
      h1=(res&0x0000ffff)*inv; h1&=0x0000ffff; res+=h1*p32; res>>=16;
    }
    ex>>=1;
    if (!ex) {
      h1=(res&0x0000ffff)*inv; h1&=0x0000ffff; res+=h1*p32; res>>=16;
/*
if (((ushort)res)!=rrr)
  complain("mpqs_powmod %u %u %u  %u %u\n",(u32_t)a,(u32_t)e,p32,(u32_t)rrr,res);
*/
      return (ushort)res;
    }
    aa*=aa;
    h1=(aa&0x0000ffff)*inv; h1&=0x0000ffff; aa+=h1*p32; aa>>=16;
  }
  return 0; /* never reached */
}


static u32_t mpqs_powmod_mm(u32_t a, u32_t e, u32_t p32, u32_t mmi, u32_t one)
{
  u32_t aa, res, ex=e, h32, hh32;

  res=one;
  if (!ex) return res;

  aa=a;
  while (1) {
    if (ex&1) {
      h32=res*aa;
      MMREDUCE;
      res=h32;
    }
    ex>>=1;
    if (!ex) return res;
    h32=aa*aa;
    MMREDUCE;
    aa=h32;
  }
  return 0; /* never reached */
}


static ushort mpqs_sqrt_init(ushort p)
{
  ushort e, u, i, j;
  // SMJS initialised g to stop compiler warning
  u32_t g = 0, b;

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


static ushort mpqs_sqrt0(ushort a, ushort p, ushort help)
                 /* 0<a<p, p prime, p!=2, (a/p)=1 */
{
  ushort e, u, i, l;
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
  r=mpqs_powmod(a,u>>1,p);
  k=r*r; k%=(u32_t)p;
  k*=a; k%=(u32_t)p;
  r*=a; r%=(u32_t)p;
  while (k!=1) {
    b=k;
    for (i=0; i<e && b!=1; i++) { b*=b; b%=(u32_t)p; }
    for (l=i+1; l<e; l++) { g*=g; g%=(u32_t)p; }
    r*=g; r%=(u32_t)p;
    g*=g; g%=(u32_t)p;
    k*=g; k%=(u32_t)p;
    e=i;
  }
  return (ushort)r;
}


static ushort mpqs_sqrt(ushort a, ushort p, ushort help)
                 /* 0<a<p, p prime, p!=2, (a/p)=1 */
{
  u32_t e, u, i, l;
  u32_t b, g, r, k,aa;
  u32_t p32, mmi, h, h32, hh32, one, mminit;

  p32=(u32_t)p;
  mmi=(u32_t)mpqs_256_inv_table[(p&0xff)>>1];
  h=p32*mmi; h&=0xff00;
  h*=mmi; mmi-=h;
  h=p32*mmi; h&=0xffff0000;
  h*=mmi; mmi-=h;
//  if (mmi*p32!=1) Schlendrian("mpqs_sqrt_inv");
  mmi=-mmi;
  mminit=(0xffffffff)%p32; mminit++; /* mminit!=p32 */

  aa=(u32_t)a;
  h32=aa*mminit; MMREDUCE; aa=h32; if (aa>=p32) aa-=p32;
  h32=mminit; MMREDUCE; one=h32; if (one>=p32) one-=p32;

  if (p&2) {
    aa=mpqs_powmod_mm(aa,(p32+1)>>2,p32,mmi,one);
    h32=aa;
    MMREDUCE;
    if (h32>=p32) h32-=p32;
    return (ushort)h32;
  }
  if (p&4) {
    b=mpqs_powmod_mm(aa,(p32+3)>>3,p32,mmi,one);
    h32=b*b; MMREDUCE; g=h32; if (g>=p32) g-=p32;
    h32=b;
    if (g!=aa) h32*=help; 
    MMREDUCE; b=h32; if (b>=p32) b-=p32;
    return (ushort)b;
  }
  e=0; u=p32-1;
  while (!(u&1)) { e++; u>>=1; }
  h32=((u32_t)help)*mminit; MMREDUCE; g=h32;
  r=mpqs_powmod_mm(aa,u>>1,p32,mmi,one);
  h32=r*r; MMREDUCE; k=h32;
  h32=r*aa; MMREDUCE; r=h32;
  h32=k*aa; MMREDUCE; k=h32; if (k>=p32) k-=p32;
  while (k!=one) {
    b=k;
    for (i=0; i<e && b!=one; i++) {
      h32=b*b; MMREDUCE; b=h32; if (b>=p32) b-=p32;
    }

    for (l=i+1; l<e; l++) {
      h32=g*g; MMREDUCE; g=h32;
    }
    h32=r*g; MMREDUCE; r=h32;
    h32=g*g; MMREDUCE; g=h32;
    h32=k*g; MMREDUCE; k=h32; if (k>=p32) k-=p32;
    e=i;
  }
  h32=r; MMREDUCE; r=h32; if (r>=p32) r-=p32;
  return (ushort)r;
}


static u32_t mpqs_inv(ushort a)
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

/* ---------------------------------------------------------- */


static ushort mpqs_multiplier_cand[4][MPQS_MULT_NCAND]={
 { 1, 17, 33, 41, 57, 65, 73, 89 },
 { 3, 11, 19, 35, 43, 51, 59, 67 },
 { 5, 13, 21, 29, 37, 53, 61, 69 },
 { 7, 15, 23, 31, 39, 47, 55, 71 },
};

static double mpqs_multiplier_cand_log[4][MPQS_MULT_NCAND];

static char mpqs_multiplier_cand_jacobi[4][MPQS_MULT_NCAND][MPQS_MULT_NTESTPR];

void mpqs_init()
{
  u32_t i, j, k;
  u32_t p, add, d;

  mpz_init(mpqs_N);
  mpz_init(mpqs_kN);
  mpz_init(mpqs_dummy);
  mpz_init(mpqs_gmp_help);
  mpz_init(mpqs_help);
  mpz_init(mpqs_sq1);
  mpz_init(mpqs_sq2);
  for (i=0; i<MPQS_MAX_FACTORS; i++) mpz_init(mpqs_factors[i]);
  mpqs_sievearray=(uchar *)xmalloc((1+MPQS_SIEVELENMAX)*sizeof(uchar));
  mpqs_tinyarray=(uchar *)xmalloc(32*MPQS_MAX_TINYPROD*sizeof(uchar));
  mpqs_relations=(ushort **)xmalloc(MPQS_MAX_NRELS*sizeof(ushort *));
  mpqs_relations[0]=(ushort *)xmalloc(MPQS_MAX_NRELS*MPQS_REL_ENTRIES*sizeof(ushort));
  for (i=1; i<MPQS_MAX_NRELS; i++)
    mpqs_relations[i]=mpqs_relations[i-1]+MPQS_REL_ENTRIES;
  mpqs_part_rels=(ushort **)xmalloc(MPQS_MAX_NPRELS*sizeof(ushort *));
  mpqs_part_rels[0]=(ushort *)xmalloc(MPQS_MAX_NPRELS*MPQS_REL_ENTRIES*sizeof(ushort));
  for (i=1; i<MPQS_MAX_NPRELS; i++)
    mpqs_part_rels[i]=mpqs_part_rels[i-1]+MPQS_REL_ENTRIES;
  mpqs_comb_rels=(ushort **)xmalloc(MPQS_MAX_NCRELS*sizeof(ushort *));
  mpqs_comb_rels[0]=(ushort *)xmalloc(MPQS_MAX_NCRELS*2*MPQS_REL_ENTRIES*sizeof(ushort));
  for (i=1; i<MPQS_MAX_NCRELS; i++)
    mpqs_comb_rels[i]=mpqs_comb_rels[i-1]+2*MPQS_REL_ENTRIES;
  mpqs_rel_buffer=(ushort **)xmalloc(MPQS_MAX_NRELS_BUF*sizeof(ushort *));
  mpqs_rel_buffer[0]=(ushort *)xmalloc(MPQS_MAX_NRELS_BUF*MPQS_REL_ENTRIES*sizeof(ushort));
  for (i=1; i<MPQS_MAX_NRELS_BUF; i++)
    mpqs_rel_buffer[i]=mpqs_rel_buffer[i-1]+MPQS_REL_ENTRIES;

  mpqs_rel_active=(uchar *)xmalloc(MPQS_MAX_NRELS*sizeof(uchar));
  mpqs_gauss_row64=(u64_t **)xmalloc(MPQS_GAUSS_MAX*sizeof(u64_t *));
  mpqs_gauss_row64[0]=(u64_t *)xmalloc(MPQS_GAUSS_MAX*MPQS_GAUSS_MAX/64*sizeof(u64_t));
  mpqs_gauss_row=(u32_t **)xmalloc(MPQS_GAUSS_MAX*sizeof(u32_t *));
  mpqs_gauss_row[0]=(u32_t *)xmalloc(MPQS_GAUSS_MAX*MPQS_GAUSS_MAX/32*sizeof(u32_t));
  mpqs_gauss_mat=mpqs_gauss_row[0];

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
      mpqs_prime_table[i]=(ushort)p;
      mpqs_prime_sqrthelp[i]=mpqs_sqrt_init(p);
      i++;
    }
    p+=add; add=6-add;
  }

  for (i=0; i<4; i++)
    for (j=0; j<MPQS_MULT_NCAND; j++) {
      d=mpqs_multiplier_cand[i][j];
      for (k=1; k<MPQS_MULT_NTESTPR; k++) {
        p=(u32_t)mpqs_prime_table[k];
        if (d%p) {
          mpqs_multiplier_cand_jacobi[i][j][k]=
            (char)(mpqs_jacobi((ushort)(d%p),(ushort)p));
        } else mpqs_multiplier_cand_jacobi[i][j][k]=0;
      }
      mpqs_multiplier_cand_log[i][j]=log((double)d)-1.5*log(2);
    }

  mpqs_isinit=1;
}


static void mpqs_choose_multiplier()
{
  ushort Nmod8, mult, p, mm, r;
  ushort residue[MPQS_MULT_NTESTPR];
  double value[MPQS_MULT_NCAND], v, v0, vp, vmin, dn;
  double plog[MPQS_MULT_NTESTPR];
  char jac[MPQS_MULT_NTESTPR], c;
  u32_t i, j;

  for (i=1; i<MPQS_MULT_NTESTPR; i++) {
    p=mpqs_prime_table[i];
    r=(ushort)mpz_mod_ui(mpqs_dummy,mpqs_N,(ulong)p);
    jac[i]=(char)mpqs_jacobi(r,p);
  }
  Nmod8=(ushort)mpz_mod_ui(mpqs_dummy,mpqs_N,8);
  dn=log(mpz_get_d(mpqs_N))/log(2.);
  dn=100.9-dn;
  if (dn>7) mm=128; else mm=(ushort)(exp(dn*log(2)));
  for (i=1; i<MPQS_MULT_NTESTPR; i++) {
    p=mpqs_prime_table[i];
    plog[i]=log((double)p)/((double)p);
  }
  for (j=0; j<MPQS_MULT_NCAND; j++) {
    mult=mpqs_multiplier_cand[Nmod8/2][j];
    if (mult<=mm) {
      v=mpqs_multiplier_cand_log[Nmod8/2][j];
      for (i=1; i<MPQS_MULT_NTESTPR; i++) {
        c=mpqs_multiplier_cand_jacobi[Nmod8/2][j][i];
        if (c) {
          vp=plog[i];
          if (c==jac[i]) v-=vp; else v+=vp;
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
// printf("%u ",mpqs_multiplier);
#endif
  mpz_mul_ui(mpqs_kN,mpqs_N,(ulong)mpqs_multiplier);
  mpz_fdiv_q_2exp(mpqs_dummy,mpqs_kN,mp_bits_per_limb);
  mpqs_kN_64=(u64_t)mpz_get_ui(mpqs_dummy); mpqs_kN_64<<=mp_bits_per_limb;
  if (mp_bits_per_limb>=64) mpqs_kN_64=0;
  mpqs_kN_64+=(u64_t)mpz_get_ui(mpqs_kN);
}


static void mpqs_choose_parameter(int retry)
{
  int n, r;

  n=(int)(mpqs_complexity+0.5);
  if (n>96) n=96;
  if (n<43) n=43;
  n=(n+3)/4;
  n-=11;  /* 0<=n<=13 */

  r=retry;
  if ((r) && (n<13)) n++;
//  if ((r) && (n<11)) { n++; r--; }
#ifdef MPQS_STAT
stat_size[n]++;
#endif
  mpqs_nFB=mpqs_param[n][0];
  if (r) mpqs_nFB+=r*mpqs_nFB/4;
  if (mpqs_nFB>=MPQS_MAX_FBSIZE) mpqs_nFB=MPQS_MAX_FBSIZE-1;
  mpqs_sievebegin=mpqs_param[n][1];
  mpqs_nAdiv=mpqs_param[n][2];
  mpqs_nAdiv_total=mpqs_param[n][3];
  mpqs_accept=mpqs_param[n][4];
  mpqs_td_begin=mpqs_param[n][5];

  mpqs_sievelen=mpqs_param[n][6];
  if (r) {
    if (mpqs_sievelen<MPQS_SIEVELENMAX>>r) mpqs_sievelen<<=r;
    else mpqs_sievelen=MPQS_SIEVELENMAX;
  }
  mpqs_disp=mpqs_sievelen/2;
  mpqs_nfactors=0; mpqs_nrels=0; mpqs_nlp=0; mpqs_nprels=0; mpqs_ncrels=0;
  for (n=0; n<128; n++) mpqs_hash_table[n][0]=0;
  for (n=0; n<256; n++) mpqs_rel_hash_table[n][0]=0;
}


static void mpqs_generate_FB()
{
  ushort *fb, p, rest, i, nfb;

  fb=mpqs_FB; nfb=0; i=0; mpqs_nFBk=0;
  *fb++=2; *fb++=1; nfb++;
  for (i=1; i<MPQS_MAX_NPRIMES; i++) {
    p=mpqs_prime_table[i];
    if (p>MPQS_FB_MAXPRIME) break;
    rest=(ushort)mpz_mod_ui(mpqs_dummy,mpqs_kN,(ulong)p);
    if (rest) {
//zeitA(25);
      if (mpqs_jacobi(rest,p)==1) {
        *fb++=p;
        *fb++=mpqs_sqrt(rest,p,mpqs_prime_sqrthelp[i]);
        nfb++;
        if (nfb>=mpqs_nFB) { /*zeitB(25);*/ break; }
      }
//zeitB(25);
    } else {
      if (mpqs_multiplier%(u32_t)p)
        complain("mpqs: N has small divisor: %u\n",p);
      mpqs_FBk[mpqs_nFBk++]=p;
    }
  }
  nfb-=((nfb-mpqs_nAdiv_total)&3);
  mpqs_nFB=nfb;
  mpqs_nFBk_1=mpqs_nFBk+1;
  mpqs_nsp=1+mpqs_nFB+mpqs_nFBk;
  mpqs_pmax=mpqs_FB[2*(nfb-1)];
  mpqs_FB[2*nfb]=mpqs_sievelen;
  mpqs_FB[2*nfb]=0xffff;
}


static int mpqs_SI_init()
{
  double d;
  u32_t a, i, j, k, l, n;
  long long prod;
  ushort inv, p, *fb=mpqs_FB, pk, tp_max;
  double A_value[20], A_div_log[MPQS_MAX_ADIV_ALL]; /* maximal number of A's */
  double v, vmax;
  uchar lo;
  u16_t pb;

  d=mpz_get_d(mpqs_kN); mpqs_kN_dbl=d;
  d*=8.; d=sqrt(d); d/=(double)mpqs_sievelen;
  d=log(d); d/=(double)mpqs_nAdiv; d=exp(d);
  a=(u32_t)d;
/*  if (a>=mpqs_pmax) { printf("%u %u ",a,mpqs_pmax); return -1; }*/
  if (a>=mpqs_pmax) a=mpqs_pmax-1;
  for (i=0; i<mpqs_nFB; i++) if (mpqs_FB[2*i]>a) break;
  if (i>mpqs_nAdiv_total/2) i-=mpqs_nAdiv_total/2; else i=1;
  if (i<1) i=1; /* first prime >2 */
  if (mpqs_FB[2*i]<3) return -2;
  if (i+mpqs_nAdiv_total>mpqs_nFB) {
    if (mpqs_nFB<=mpqs_nAdiv_total+1) return -3;
    i=mpqs_nFB-mpqs_nAdiv_total-1;
  }
  for (j=0; j<mpqs_nAdiv_total; j++) {
    mpqs_Adiv_all[j]=mpqs_FB[2*(i+j)];
    mpqs_Adiv_all_sqrt[j]=mpqs_FB[2*(i+j)+1];
  }
  for (j=i+mpqs_nAdiv_total; j<mpqs_nFB; j++) {
    mpqs_FB[2*(j-mpqs_nAdiv_total)]=mpqs_FB[2*j];
    mpqs_FB[2*(j-mpqs_nAdiv_total)+1]=mpqs_FB[2*j+1];
  }
  mpqs_nFB-=mpqs_nAdiv_total;
/* tiny */
  mpqs_sievebegin=1; mpqs_sievebegink=0; mpqs_tiny_prod=1; mpqs_ntiny=0;
  tp_max=mpqs_sievelen>>4;
  if (tp_max>MPQS_MAX_TINYPROD) tp_max=MPQS_MAX_TINYPROD;
  while (1) {
    if (mpqs_sievebegin<mpqs_nFB) p=mpqs_FB[2*mpqs_sievebegin];
    else p=tp_max+1;
    if (mpqs_sievebegink<mpqs_nFBk) pk=mpqs_FBk[mpqs_sievebegink];
    else pk=tp_max+1;
    if (p<pk) {
      if (mpqs_tiny_prod*(u32_t)p>=tp_max) break;
      mpqs_sievebegin++; mpqs_ntiny++;
      mpqs_tiny_prod*=(u32_t)p;
    } else {
      if (mpqs_tiny_prod*(u32_t)pk>=tp_max) break;
      mpqs_sievebegink++; mpqs_ntiny++;
      mpqs_tiny_prod*=(u32_t)pk;
    }
  }
#ifdef MPQS_STAT
//for (j=0; j<mpqs_sievebegin; j++) printf("%u ",mpqs_FB[2*j]);
//printf("\n");
#endif

/*for (j=0; j<mpqs_sievebegin; j++) printf("%u ",mpqs_FB[2*j]);
printf("\n");*/
  if (mpqs_nFB<mpqs_td_begin) {
    if (mpqs_nFB&3) complain("mpqs_nFB\n");
//printf("%u\n",mpqs_nFB);
    mpqs_td_begin=mpqs_nFB;
  }

#ifdef HAVE_XMM_MUL
  while (mpqs_td_begin&3) mpqs_td_begin++;
  p=mpqs_FB[2];
  mpqs_FB_inv[1]=mpqs_inv(p);
  mpqs_FB_inv_info[2]=p; mpqs_FB_inv_info[3]=p;
  mpqs_FB_inv_info[10]=(ushort)mpqs_FB_inv[1];
  mpqs_FB_inv_info[11]=(ushort)mpqs_FB_inv[1];
  p=mpqs_FB[4];
  mpqs_FB_inv[2]=mpqs_inv(p);
  mpqs_FB_inv_info[4]=p; mpqs_FB_inv_info[5]=p;
  mpqs_FB_inv_info[12]=(ushort)mpqs_FB_inv[2];
  mpqs_FB_inv_info[13]=(ushort)mpqs_FB_inv[2];
  p=mpqs_FB[6];
  mpqs_FB_inv[3]=mpqs_inv(p);
  mpqs_FB_inv_info[6]=p; mpqs_FB_inv_info[7]=p;
  mpqs_FB_inv_info[14]=(ushort)mpqs_FB_inv[3];
  mpqs_FB_inv_info[15]=(ushort)mpqs_FB_inv[3];
  if (mpqs_td_begin&3) complain("mpqs_td_begin is not divisible by 4\n");
  for (j=4; j<mpqs_td_begin; j+=4) {
    p=mpqs_FB[2*j];
    mpqs_FB_inv[j]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j]=p; mpqs_FB_inv_info[4*j+1]=p;
    mpqs_FB_inv_info[4*j+8]=(ushort)mpqs_FB_inv[j];
    mpqs_FB_inv_info[4*j+9]=(ushort)mpqs_FB_inv[j];
    p=mpqs_FB[2*(j+1)];
    mpqs_FB_inv[j+1]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j+2]=p; mpqs_FB_inv_info[4*j+3]=p;
    mpqs_FB_inv_info[4*j+10]=(ushort)mpqs_FB_inv[j+1];
    mpqs_FB_inv_info[4*j+11]=(ushort)mpqs_FB_inv[j+1];
    p=mpqs_FB[2*(j+2)];
    mpqs_FB_inv[j+2]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j+4]=p; mpqs_FB_inv_info[4*j+5]=p;
    mpqs_FB_inv_info[4*j+12]=(ushort)mpqs_FB_inv[j+2];
    mpqs_FB_inv_info[4*j+13]=(ushort)mpqs_FB_inv[j+2];
    p=mpqs_FB[2*(j+3)];
    mpqs_FB_inv[j+3]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j+6]=p; mpqs_FB_inv_info[4*j+7]=p;
    mpqs_FB_inv_info[4*j+14]=(ushort)mpqs_FB_inv[j+3];
    mpqs_FB_inv_info[4*j+15]=(ushort)mpqs_FB_inv[j+3];
  }
#else
  if (mpqs_td_begin&1) mpqs_td_begin++;
  p=mpqs_FB[2];
  mpqs_FB_inv[1]=mpqs_inv(p);
  mpqs_FB_inv_info[2]=p; mpqs_FB_inv_info[3]=p;
  mpqs_FB_inv_info[6]=(ushort)mpqs_FB_inv[1];
  mpqs_FB_inv_info[7]=(ushort)mpqs_FB_inv[1];
  if (mpqs_td_begin&1) complain("mpqs_td_begin is odd\n");
  for (j=2; j<mpqs_td_begin; j+=2) {
    p=mpqs_FB[2*j];
    mpqs_FB_inv[j]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j]=p; mpqs_FB_inv_info[4*j+1]=p;
    mpqs_FB_inv_info[4*j+4]=(ushort)mpqs_FB_inv[j];
    mpqs_FB_inv_info[4*j+5]=(ushort)mpqs_FB_inv[j];
    p=mpqs_FB[2*(j+1)];
    mpqs_FB_inv[j+1]=mpqs_inv(p);
    mpqs_FB_inv_info[4*j+2]=p; mpqs_FB_inv_info[4*j+3]=p;
    mpqs_FB_inv_info[4*j+6]=(ushort)mpqs_FB_inv[j+1];
    mpqs_FB_inv_info[4*j+7]=(ushort)mpqs_FB_inv[j+1];
  }
#endif
  for (; j<mpqs_nFB; j++) mpqs_FB_inv[j]=mpqs_inv(mpqs_FB[2*j]);
  for (j=0; j<mpqs_nFBk; j++) mpqs_FBk_inv[j]=mpqs_inv(mpqs_FBk[j]);
  for (j=0; j<mpqs_nAdiv_total; j++)
    mpqs_FB_A_inv[j]=mpqs_inv(mpqs_Adiv_all[j]);

/* for next_pol: */
{
/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/

  u32_t mmi, h32, hh32, invh, p32;
  ushort invhelp[MPQS_MAX_ADIV_ALL], hhh;

  for (j=1; j<mpqs_nFB; j++) {
    p=mpqs_FB[2*j]; p32=(u32_t)p;
    mmi=-mpqs_FB_inv[j];
    invhelp[0]=mpqs_Adiv_all[0]%p;
    h32=(u32_t)(invhelp[0]);
    for (i=1; i<mpqs_nAdiv_total; i++) {
      h32*=(u32_t)(mpqs_Adiv_all[i]);
      MMREDUCE;
      invhelp[i]=(ushort)h32;
    }
    MMREDUCE; if (h32>=p) h32-=p; /* divide by 2^16 */
    invh=(u32_t)mpqs_invert((ushort)h32,p);
    invh<<=1;

    for (i=mpqs_nAdiv_total-1; i; i--) {
      h32=invh; h32*=(u32_t)(invhelp[i-1]);
      MMREDUCE; if (h32>=p) h32-=p;
      mpqs_SI_inv_table[i][j]=(ushort)h32;
      h32=invh; h32*=(u32_t)mpqs_Adiv_all[i];
      MMREDUCE; invh=h32;
    }
    if (invh>=p) invh-=p;
    mpqs_SI_inv_table[0][j]=(ushort)invh;
  }
  for (j=0; j<mpqs_nAdiv_total; j++) {
    u32_t h;
    ushort pa;

    p=mpqs_Adiv_all[j];
    for (i=0; i<mpqs_nAdiv_total; i++) {
      pa=mpqs_Adiv_all[i]%p;
      if (pa) {
        h=(u32_t)mpqs_invert(pa,p); //h*=2; if (h>=p) h-=p;
        h<<=17; h%=((u32_t)p);
      } else h=0;
      mpqs_Adiv_SI_inv_table[i][j]=(ushort)h;
    }
    mpqs_Adiv_disp[j]=(ushort)(mpqs_disp%((u32_t)p));
  }
}

  for (j=1; j<mpqs_nFB; j++) {
    mpqs_FB_mm_inv[j]=-mpqs_FB_inv[j];
    p=mpqs_FB[2*j];
    mpqs_FB_disp[j]=(ushort)(mpqs_disp%((u32_t)p));
  }
#ifdef HAVE_XMM_MUL
  for (j=0; j<mpqs_nFB; j+=8) {
    for (i=0; i<8; i++) mpqs_FB_np_px[2*j+i]=mpqs_FB[2*j+2*i];
    for (i=0; i<8; i++) mpqs_FB_np_px[2*j+8+i]=mpqs_FB[2*j+1+2*i];
  }
#else
  for (j=0; j<mpqs_nFB; j+=4) {
    mpqs_FB_np_p[2*j]=mpqs_FB[2*j];
    mpqs_FB_np_p[2*j+1]=mpqs_FB[2*j+2];
    mpqs_FB_np_p[2*j+2]=mpqs_FB[2*j+4];
    mpqs_FB_np_p[2*j+3]=mpqs_FB[2*j+6];
    mpqs_FB_np_p[2*j+4]=mpqs_FB[2*j+1];
    mpqs_FB_np_p[2*j+5]=mpqs_FB[2*j+3];
    mpqs_FB_np_p[2*j+6]=mpqs_FB[2*j+5];
    mpqs_FB_np_p[2*j+7]=mpqs_FB[2*j+7];
  }
#endif

/* compute log-approximations */
  d=mpz_get_d(mpqs_kN);
  d/=8.; d=sqrt(d); d*=(double)mpqs_sievelen;
  mpqs_maxval_dbl=d;
  d=log(d);
  d-=log((double)(1<<mpqs_accept));
  mpqs_logbound=128;  /* ??? */
  d=(double)(mpqs_logbound)/d;
  for (i=0; i<mpqs_nFB; i++) {
    p=mpqs_FB[2*i];
    mpqs_FB_log[i]=(uchar)(0.5+d*log((double)p));
  }
  for (i=0; i<mpqs_nFBk; i++) {
    p=mpqs_FBk[i];
    mpqs_FBk_log[i]=(uchar)(0.5+d*log((double)p));
  }
  for (i=0; i<mpqs_nAdiv_total; i++) {
    p=mpqs_Adiv_all[i];
    mpqs_Adiv_log[i]=(uchar)(0.5+d*log((double)p));
  }

/* fill mpqs_FB0 */
  lo=0; n=0; i=mpqs_sievebegin;
#if 1
  pb=mpqs_sievelen>>3;
  for (; i<mpqs_nFB; i++) {
    p=mpqs_FB[2*i]; if (p>pb) break;
    if (lo<mpqs_FB_log[i]) {
      mpqs_FB0[n++]=0;
      lo=mpqs_FB_log[i];
      mpqs_FB0[n++]=(u16_t)lo;
    }
    mpqs_FB0[n++]=p;
  }
  mpqs_FB0[n++]=0; mpqs_FB0[n++]=0;
#endif
  for (j=4; j; j--) {
    pb=mpqs_sievelen/j;
    for (; i<mpqs_nFB; i++) {
      p=mpqs_FB[2*i]; if (p>pb) break;
      if (lo<mpqs_FB_log[i]) {
        mpqs_FB0[n++]=0;
        lo=mpqs_FB_log[i];
        mpqs_FB0[n++]=(u16_t)lo;
      }
      mpqs_FB0[n++]=p;
    }
    mpqs_FB0[n++]=0; mpqs_FB0[n++]=0;
  }
  for (; i<mpqs_nFB; i++) {
    p=mpqs_FB[2*i];
    if (lo<mpqs_FB_log[i]) {
      mpqs_FB0[n++]=0;
      lo=mpqs_FB_log[i];
      mpqs_FB0[n++]=(u16_t)lo;
    }
    mpqs_FB0[n++]=p;
  }
  mpqs_FB0[n++]=0; mpqs_FB0[n++]=0;

#if 0
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
  i=0;
  while (i<j) {
    v=0;
    for (k=0; k<mpqs_nAdiv_total; k++)
      if (mpqs_A_mask[i]&(1<<k)) v+=A_div_log[k];
    if (d-v>63.95) {
      j--; /* printf(":");*/
      mpqs_A_mask[i]=mpqs_A_mask[j];
    } else i++;
  }
#else
{
/* experimental */
  double *alog;
  ushort tmp;
  u32_t *asort;

/* compute mask for choice of A-divisors */
  j=0;
  if (mpqs_nAdiv_total>8) {
    k=0; if (mpqs_nAdiv>8) k=mpqs_nAdiv-8;
    for (; k<=mpqs_nAdiv; k++) {
      u32_t i0, i1;
      ushort m16;

      for (i0=mpqs_A_table_n[k]; i0<mpqs_A_table_n[k+1]; i0++) {
        m16=(ushort)(mpqs_A_table[i0]); m16<<=8;
        if (m16>>mpqs_nAdiv_total) break;
        for (i1=mpqs_A_table_n[mpqs_nAdiv-k];
             i1<mpqs_A_table_n[mpqs_nAdiv-k+1]; i1++) {
          mpqs_A_mask[j]=m16+(ushort)(mpqs_A_table[i1]);
          j++;
        }
      }
    }
  } else {
    for (i=mpqs_A_table_n[mpqs_nAdiv];
         i<mpqs_A_table_n[mpqs_nAdiv+1]; i++) {
      ushort m16;

      m16=(ushort)(mpqs_A_table[i]);
      if (m16>>mpqs_nAdiv_total) break;
      mpqs_A_mask[j]=m16;
      j++;
    }
  }
/* sort mpqs_A_mask */
  alog=(double *)xmalloc(j*sizeof(double));

/* ensure that |C|<2^96 */
  for (i=0; i<mpqs_nAdiv_total; i++)
    A_div_log[i]=log((double)(mpqs_Adiv_all[i]))/log(2.);
  d=log(mpqs_kN_dbl)/log(2.);
  for (i=0; i<j; i++) {
    v=0;
    for (k=0; k<mpqs_nAdiv_total; k++)
      if (mpqs_A_mask[i]&(1<<k)) v+=A_div_log[k];
    alog[i]=v;
    if (d-v>62.95) {
      j--; /* printf(":");*/
      mpqs_A_mask[i]=mpqs_A_mask[j];
      i--;
    }
  }

  d=mpz_get_d(mpqs_kN); mpqs_kN_dbl=d;
  d*=8.; d=sqrt(d); d/=(double)mpqs_sievelen;
  d=log(d)/log(2.);
  for (i=0; i<j; i++)
    if (d>alog[i]) alog[i]=d-alog[i];
    else alog[i]-=d;

{
  u32_t *a0, *a1, *atmp, zp, i0, i1, n0, n1;

  asort=(u32_t *)xmalloc(2*j*sizeof(*asort));
  for (i=0; i<j; i++) asort[i]=i;
  a0=asort; a1=asort+j;
  for (i=0; i+1<j; i+=2) {
    if (alog[a0[i]]<alog[a0[i+1]]) { a1[i]=a0[i]; a1[i+1]=a0[i+1]; }
    else { a1[i]=a0[i+1]; a1[i+1]=a0[i]; }
  }
  if (i<j) a1[i]=a0[i];
  zp=2;
  while (zp<j) {
    atmp=a0; a0=a1; a1=atmp;
    for (i=0; i<j;) {
      i0=i; i1=i+zp; n0=i1; n1=n0+zp;
      if (n0>j) n0=j;
      if (n1>j) n1=j;
      while (1) {
        if (i0>=n0) {
          while (i1<n1) a1[i++]=a0[i1++];
          break;
        }
        if (i1>=n1) {
          while (i0<n0) a1[i++]=a0[i0++];
          break;
        }
        if (alog[a0[i0]]<alog[a0[i1]]) a1[i++]=a0[i0++];
        else a1[i++]=a0[i1++];
      }
    }
    zp+=zp;
  }
  for (i=0; i<j; i++) a1[i]=(u32_t)(mpqs_A_mask[a1[i]]);
  for (i=0; i<j; i++) mpqs_A_mask[i]=(ushort)(a1[i]);
  free(asort);
}
  free(alog);
}
#endif

  mpqs_nA=j; mpqs_A_index=-1; mpqs_B_index=(1<<mpqs_nAdiv)-1;

  mpqs_FB[2*mpqs_nFB]=mpqs_sievelen;
  mpqs_FB[2*mpqs_nFB]=0xffff;
  return 1;
}


static int mpqs_next_pol()
{
  u32_t i, ind, j;
  long long add, bi;
  ushort p, s1, s2, *fb=mpqs_FB, *fbs=mpqs_FB_start;
  u32_t mod, inv;
  short sh;
  u32_t bb, cc, c1, c2;
  ushort mask;
  double d;
  ushort aa;
  u32_t al, *ul;
  u32_t B_modp;

  mpqs_B_index++;
  if (mpqs_B_index>=1<<(mpqs_nAdiv-1)) {
    u32_t Bi_mul[MPQS_MAX_ADIV_ALL], Adiv_nr[MPQS_MAX_ADIV_ALL];

    mpqs_A_index++;
    if (mpqs_A_index>=mpqs_nA) return 0;
#ifdef MPQS_ZEIT
  zeitA(7);
#endif
    mask=mpqs_A_mask[mpqs_A_index];
    for (i=0; i<mpqs_nAdiv_total; i++)
      if (mask & (1<<i)) mpqs_Adiv_active[i]=1;
      else mpqs_Adiv_active[i]=0;
    j=0;
    for (i=0; i<mpqs_nAdiv_total; i++)
      if (mpqs_Adiv_active[i]) {
        mpqs_Adiv[j]=mpqs_Adiv_all[i];
        mpqs_Adiv_sqrt[j]=mpqs_Adiv_all_sqrt[i];
        Adiv_nr[j]=i;
        j++;
      }
    mpqs_A=1;
    for (i=0; i<mpqs_nAdiv; i++) mpqs_A*=(long long)(mpqs_Adiv[i]);
/* precomputations for sieve init */
    mpqs_sinit_tab[0]=0; mpqs_sinit_tab[1]=0;
    i=1;
    for (j=1; j<MPQS_SI_NVALUES; j++) {
      d=sqrt(mpqs_kN_dbl+mpqs_maxval_dbl*((double)(mpqs_A>>j)))/((double)mpqs_A);
      if (d<(double)mpqs_disp) s1=mpqs_disp-(ushort)d; else s1=0;
      s1-=(s1&63);
      if (s1>mpqs_sinit_tab[2*i-2]) {
        mpqs_sinit_tab[2*i]=s1; mpqs_sinit_tab[2*i+1]=j;
        i++;
      }
    }
    for (j=1; j<MPQS_SI_NVALUES; j++) {
      d=mpqs_kN_dbl-mpqs_maxval_dbl*((double)(mpqs_A>>(MPQS_SI_NVALUES-j)));
      if (d<0.) break;
      d=sqrt(d)/((double)mpqs_A);
      if (d<(double)mpqs_disp) s1=mpqs_disp-(ushort)d; else s1=0;
      s1-=(s1&63);
      if (s1>mpqs_sinit_tab[2*i-2]) {
        mpqs_sinit_tab[2*i]=s1; mpqs_sinit_tab[2*i+1]=MPQS_SI_NVALUES-1-j;
        i++;
      }
    }

    if (mpqs_sinit_tab[2*i-2]==mpqs_disp) i--;
    for (j=i; j<2*i; j++) {
      mpqs_sinit_tab[2*j]=mpqs_sievelen-mpqs_sinit_tab[2*(2*i-1-j)];
      if (j<2*i-1) mpqs_sinit_tab[2*j+1]=mpqs_sinit_tab[2*(2*i-2-j)+1];
    }
/* compute B_i's */
    for (i=0; i<mpqs_nAdiv; i++) {
      p=mpqs_Adiv[i];
      bi=mpqs_A/(long long)(p);
      inv=(u32_t)(bi%(long long)p); /* bi>0 */
      inv=mpqs_invert((ushort)inv,p);
      mod=(inv*((u32_t)mpqs_Adiv_sqrt[i]))%((u32_t)p);
      Bi_mul[i]=mod;
      bi*=(long long)mod;
      mpqs_Bi[i]=bi;  /* >=0 */
    }
    mpqs_B=0;
    for (i=0; i<mpqs_nAdiv; i++) mpqs_B+=mpqs_Bi[i];
#ifdef ASM_MPQS_NEXT_POL
{
    u32_t jj;
    u64_t pi0_64, *fbs64;

    pi0_64=1ULL<<(16-mpqs_nAdiv); pi0_64+=(pi0_64<<16); pi0_64+=(pi0_64<<32);
    fbs64=(u64_t *)fbs;
#ifdef HAVE_XMM_MUL
    for (i=0; i<mpqs_nFB; i+=8) {
      *fbs64++=pi0_64;
      *fbs64++=pi0_64;
      *fbs64++=0;
      *fbs64++=0;
    }
    for (j=0; j<mpqs_nAdiv; j++) {
      jj=Adiv_nr[j];
      asm_next_pol10_xmm((mpqs_nFB+7)/8,&(mpqs_SI_inv_table[jj][0]),
                         &(mpqs_SI_add[j][0]),Bi_mul[j]);
    }
    asm_next_pol11_xmm((mpqs_nFB+7)/8);
#else
    for (i=0; i<mpqs_nFB; i+=4) {
      *fbs64++=pi0_64;
      *fbs64++=0;
    }
    for (j=0; j<mpqs_nAdiv; j++) {
      jj=Adiv_nr[j];
      asm_next_pol10((mpqs_nFB+3)/4,&(mpqs_SI_inv_table[jj][0]),
                      &(mpqs_SI_add[j][0]),Bi_mul[j]);
    }
    asm_next_pol11((mpqs_nFB+3)/4);
#endif
}
#else

/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/
{
    u32_t mmi, h32, hh32, p32;
    u32_t cc1, cc2, bbb;
    u32_t jj, pi, pi0;

    pi0=1<<(16-mpqs_nAdiv);
    for (i=1; i<mpqs_nFB; i++) {
      p=fb[2*i]; p32=(u32_t)p;
      mmi=mpqs_FB_mm_inv[i];

      bbb=0; pi=pi0;
      for (j=0; j<mpqs_nAdiv; j++) {
        jj=Adiv_nr[j];
        cc=mpqs_SI_inv_table[jj][i];
        cc*=Bi_mul[j]; h32=cc;
        MMREDUCE; cc=h32; if (cc>=p) cc-=p;
        mpqs_SI_add[j][i]=(ushort)cc;
        bbb+=mpqs_SI_add[j][i]; if (bbb>=p) bbb-=p;
        pi*=((u32_t)mpqs_SI_inv_table[jj][i]); h32=pi;
        MMREDUCE; pi=h32;
      }
      cc=bbb;
      if (cc&1) cc+=p; cc>>=1;
      cc=mpqs_FB_disp[i]+(p-cc); if (cc>=p) cc-=p;
      cc1=fb[2*i+1];
      h32=cc1*pi; MMREDUCE; cc1=h32;
      if (cc1>=p) cc1-=p;
      cc2=p-cc1;
      cc1+=cc; if (cc1>=p) cc1-=p;
      cc2+=cc; if (cc2>=p) cc2-=p;
      fbs[2*i]=(ushort)cc1; fbs[2*i+1]=(ushort)cc2;
    }
}
#endif

    mpqs_2B=2*mpqs_B;
    mpqs_A_inv_64=mpqs_inv_64(mpqs_A);
    mpqs_C=mpqs_B*mpqs_B;
    mpqs_C-=mpqs_kN_64;
    mpqs_C*=mpqs_A_inv_64;

/* divisors of A */
    for (i=0; i<mpqs_nAdiv_total; i++)
      if (!mpqs_Adiv_active[i]) {
/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/

        u32_t mmi, h32, hh32, pi, jj, p32;

        p=mpqs_Adiv_all[i]; p32=(u32_t)p;
        mmi=-mpqs_inv(p);

        bb=0; pi=1<<(16-mpqs_nAdiv);
        for (j=0; j<mpqs_nAdiv; j++) {
          jj=Adiv_nr[j];
          cc=mpqs_Adiv_SI_inv_table[jj][i];
          cc*=Bi_mul[j]; h32=cc;
          MMREDUCE; cc=h32; if (cc>=p) cc-=p;
          mpqs_Adiv_SI_add[j][i]=(ushort)cc;
          bb+=mpqs_Adiv_SI_add[j][i]; if (bb>=p) bb-=p;
          pi*=((u32_t)mpqs_Adiv_SI_inv_table[jj][i]); h32=pi;
          MMREDUCE; pi=h32;
        }
        cc=bb;
        if (cc&1) cc+=p; cc>>=1;
        cc=mpqs_Adiv_disp[i]+(p-cc); if (cc>=p) cc-=p;
        c1=mpqs_Adiv_all_sqrt[i];
        h32=c1*pi; MMREDUCE; c1=h32;
        if (c1>=p) c1-=p;
        c2=p-c1;
        c1+=cc; if (c1>=p) c1-=p;
        c2+=cc; if (c2>=p) c2-=p;
        mpqs_Adiv_start1[i]=(ushort)c1; mpqs_Adiv_start2[i]=(ushort)c2;
      } else {
/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/
// SMJS Initialise j0 to stop compiler warning, not sure about this one

        u32_t mmi, h32, hh32, jj, cc0, j0 = 0, p32;

        p=mpqs_Adiv_all[i]; p32=(u32_t)p;
        sh=(short)(mpqs_2B%(long long)p);
        if (sh>0) bb=(u32_t)sh; else bb=(u32_t)((short)p+sh);
        bb=(u32_t)mpqs_invert(bb,p);
        sh=(short)(mpqs_C%(long long)p); sh=-sh;
        if (sh>0) cc=(u32_t)sh; else cc=(u32_t)((short)p+sh);

        cc*=bb; cc%=(u32_t)p;  /* now p|2*B*cc+C */
        cc0=2*cc; if (cc0>=p) cc0-=p;
        cc+=mpqs_Adiv_disp[i]; if (cc>=p) cc-=p;
        mpqs_Adiv_start1[i]=(ushort)cc;

        mmi=-mpqs_inv(p);
        for (j=0; j<mpqs_nAdiv; j++) {
          jj=Adiv_nr[j];
          if (jj!=i) {
            cc=mpqs_Adiv_SI_inv_table[jj][i];
            cc*=Bi_mul[j]; h32=cc;
            MMREDUCE; cc=h32; if (cc>=p) cc-=p;
            mpqs_Adiv_SI_add[j][i]=(ushort)cc;
            cc0+=cc; if (cc0>=p) cc0-=p;
          } else j0=j;
        }
        if (cc0) cc0=p-cc0;
/* using relation: 2*start_value=-sum of change_values */
        mpqs_Adiv_SI_add[j0][i]=cc0;
      }
/* divisors of k */
    for (i=0; i<mpqs_nFBk; i++) {
      p=mpqs_FBk[i]; B_modp=0;
      bb=(u32_t)(((u64_t)mpqs_A)%((u64_t)p));
      bb=mpqs_invert(bb,p);
      for (j=0; j<mpqs_nAdiv; j++) {
        cc=(u32_t)((u64_t)(mpqs_Bi[j])%((u64_t)p));
        B_modp+=((u32_t)p-cc);
        cc=2*bb*cc; cc%=(u32_t)p;
        mpqs_SIk_add[j][i]=(ushort)cc;
      }
      c1=B_modp*bb+mpqs_disp; c1%=(u32_t)p;
      mpqs_FBk_start[i]=(ushort)c1;
    }
    mpqs_B_index=0;
#ifdef MPQS_ZEIT
    zeitB(7);
#endif
//printf("pol1: %lld %lld %lld \n",mpqs_A,mpqs_B,mpqs_C);
    return 1;
  }
#ifdef MPQS_ZEIT
  zeitA(7);
#endif
  ind=mpqs_B_index; i=0;
  while (1) {
    if (ind&1) break;
    i++; ind>>=1;
  }
  ind>>=1; ind++;
  add=2*mpqs_Bi[i];
  if (ind&1) {
    mpqs_B-=add;
#ifdef ASM_MPQS_NEXT_POL
#ifdef HAVE_XMM_MUL
    asm_next_pol3plus_xmm((mpqs_nFB+7)/8,&(mpqs_SI_add[i][0]));
#else
    asm_next_pol3plus((mpqs_nFB+3)/4,&(mpqs_SI_add[i][0]));
#endif
#else
    for (j=1; j<mpqs_nFB; j++) {
      p=fb[2*j]; aa=mpqs_SI_add[i][j];
      fbs[2*j]+=aa; fbs[2*j+1]+=aa;
      if (fbs[2*j]>=p) fbs[2*j]-=p;
      if (fbs[2*j+1]>=p) fbs[2*j+1]-=p;
    }
#endif
    for (j=0; j<mpqs_nFBk; j++) {
      p=mpqs_FBk[j];
      aa=mpqs_FBk_start[j]+mpqs_SIk_add[i][j];
      if (aa>=p) aa-=p;
      mpqs_FBk_start[j]=aa;
    }
  } else {
    mpqs_B+=add;
#ifdef ASM_MPQS_NEXT_POL
#ifdef HAVE_XMM_MUL
    asm_next_pol3minus_xmm((mpqs_nFB+7)/8,&(mpqs_SI_add[i][0]));
#else
    asm_next_pol3minus((mpqs_nFB+3)/4,&(mpqs_SI_add[i][0]));
#endif
#else
    for (j=1; j<mpqs_nFB; j++) {
      p=fb[2*j]; aa=p-mpqs_SI_add[i][j];
      fbs[2*j]+=aa; fbs[2*j+1]+=aa;
      if (fbs[2*j]>=p) fbs[2*j]-=p;
      if (fbs[2*j+1]>=p) fbs[2*j+1]-=p;
    }
#endif
    for (j=0; j<mpqs_nFBk; j++) {
      p=mpqs_FBk[j];
      aa=mpqs_FBk_start[j]+(p-mpqs_SIk_add[i][j]);
      if (aa>=p) aa-=p;
      mpqs_FBk_start[j]=aa;
    }
  }
/*zeitB(7);*/
  mpqs_2B=2*mpqs_B;
  mpqs_C=mpqs_B*mpqs_B;
  mpqs_C-=mpqs_kN_64;
  mpqs_C*=mpqs_A_inv_64;
//printf("pol2: %lld %lld %lld \n",mpqs_A,mpqs_B,mpqs_C);
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
      p=mpqs_Adiv_all[j];
      s1=mpqs_Adiv_start1[j];
      if (ind&1) {
        s1+=mpqs_Adiv_SI_add[i][j]; if (s1>=p) s1-=p;
      } else {
        s1+=(p-mpqs_Adiv_SI_add[i][j]); if (s1>=p) s1-=p;
      }
      mpqs_Adiv_start1[j]=s1;
    }
#ifdef MPQS_ZEIT
  zeitB(7);
#endif
  return 1;
}


static void mpqs_sieve()
{
  uchar *sv=mpqs_sievearray, *fblog=mpqs_FB_log, *svend;
  ushort *fb=mpqs_FB, *fbs=mpqs_FB_start;
  ushort p, s1, s2, end, pb;
  u32_t i;
  uchar lo;
  u32_t *ulsv, *ulsvend, m, m1, m2, *ulsv0;
  u64_t *ullsv, *ullsvend, mask0, mask1, m64[3];
#ifdef TOTAL_STAT
  double sl2se=0;
#endif

#ifdef MPQS_ZEIT
  zeitA(12);
#endif
#ifdef MPQS_ZEIT
zeitA(13);
#endif
/* initialize odd tiny primes */
  ulsv=(u32_t *)mpqs_tinyarray;
  if (mpqs_ntiny==0) {
    *ulsv=0;
  } else {
    u32_t j;
    u32_t h, *ulsv0;

    for (i=0; i<mpqs_tiny_prod; i++) *ulsv++=0;
    for (i=1; i<mpqs_sievebegin; i++) {
      p=fb[2*i]; h=mpqs_tiny_prod/((u32_t)p);
      lo=fblog[i]; sv=mpqs_tinyarray;
      s1=fbs[2*i]; s2=fbs[2*i+1];
      for (j=0; j<h; j++) {
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
      }
    }
    for (i=0; i<mpqs_sievebegink; i++) {
      p=mpqs_FBk[i]; h=mpqs_tiny_prod/((u32_t)p);
      lo=mpqs_FBk_log[i]; sv=mpqs_tinyarray;
      s1=mpqs_FBk_start[i];
      for (j=0; j<h; j++) {
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
      }
    }
  }
  memcpy(mpqs_tinyarray+4*mpqs_tiny_prod,mpqs_tinyarray,4*mpqs_tiny_prod);
  memcpy(mpqs_tinyarray+8*mpqs_tiny_prod,mpqs_tinyarray,8*mpqs_tiny_prod);
#ifdef ASM_MPQS_SIEVE_INIT16
  memcpy(mpqs_tinyarray+16*mpqs_tiny_prod,mpqs_tinyarray,16*mpqs_tiny_prod);
#endif

/* add contributions of powers of 2 */
  lo=mpqs_FB_log[0];
  m=((u32_t)mpqs_A)*(-1-mpqs_disp)+(u32_t)mpqs_2B;
  m*=(-1-mpqs_disp);
  m+=(u32_t)mpqs_C;
  m1=((u32_t)mpqs_A)*(-3-2*mpqs_disp)+(u32_t)mpqs_2B;
  m2=2*((u32_t)mpqs_A);

  m64[0]=0; m64[1]=0;
  sv=(uchar *)(&(m64[0]));
  for (i=0; i<16; i++) {
    m1+=m2; m+=m1;
    if (m&1) continue;
    if (m&6) complain("si mod 8\n");
    if (m&8) { sv[i]+=3*lo; continue; }
    if (m&16) { sv[i]+=4*lo; continue; }
    sv[i]+=6*lo;
  }

  mask0=m64[0]; mask1=m64[1];
  ullsv=(u64_t *)mpqs_tinyarray;
#ifdef ASM_MPQS_SIEVE_INIT16
  ullsvend=ullsv+4*mpqs_tiny_prod;
#else
  ullsvend=ullsv+2*mpqs_tiny_prod;
#endif
  while (ullsv<ullsvend) {
    *ullsv+++=mask0;
    *ullsv+++=mask1;
  }
#ifdef MPQS_ZEIT
zeitB(13);
#endif

/* contribution of log |q(x)| */
#ifdef MPQS_ZEIT
zeitA(14);
#endif
#ifdef ASM_MPQS_SIEVE_INIT
  mask0=0x0101010101010101ULL*(u64_t)(lo);
#ifdef ASM_MPQS_SIEVE_INIT16
  asm_sieve_init16(mpqs_sievearray,mpqs_sievelen,mpqs_sinit_tab,&mask0,
                   mpqs_tinyarray,mpqs_tiny_prod);
#else
  asm_sieve_init(mpqs_sievearray,mpqs_sievelen,mpqs_sinit_tab,&mask0,
                 mpqs_tinyarray,mpqs_tiny_prod);
#endif
#else
{
  u64_t *tptr0, *tptr1, *tptr, mask;

  tptr0=(u64_t *)mpqs_tinyarray;
  tptr1=tptr0+2*mpqs_tiny_prod;
  tptr=tptr0;
  i=0;
  ullsv=(u64_t *)mpqs_sievearray;
  while (1) {
    if (mpqs_sinit_tab[2*i]==mpqs_sievelen) break;
    ullsvend=((u64_t *)mpqs_sievearray)+mpqs_sinit_tab[2*i+2]/8;
    mask=0x0101010101010101ULL*(u64_t)(lo*mpqs_sinit_tab[2*i+1]);

    while (ullsv<ullsvend) {
      *ullsv++=*tptr+++mask;
      *ullsv++=*tptr+++mask;
      if (tptr==tptr1) tptr=tptr0;
    }
    i++;
  }
}
#endif
#ifdef MPQS_ZEIT
zeitB(14);
#endif

#ifdef MPQS_ZEIT
  zeitB(12);
#endif

#ifdef ASM_MPQS_SIEVE
#ifdef ASM_MPQS_SIEVE0
  asm_sievea();
#else
  asm_sieve();
#endif
#else
  pb=mpqs_sievelen>>2;
  i=mpqs_sievebegin; fb+=2*i; fbs+=2*i;
  for (;;) {
    p=*fb; if (p>pb) break;
#ifdef TOTAL_STAT
    sl2se+=2/(double)p;
#endif
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
  svend=mpqs_sievearray+mpqs_sievelen;
  pb=mpqs_sievelen/3;
  for (;;) {
    p=*fb; if (p>pb) break;
#ifdef TOTAL_STAT
    sl2se+=2/(double)p;
#endif
    lo=fblog[i]; sv=mpqs_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++; 
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }
  pb=mpqs_sievelen>>1;
  for (;;) {
    p=*fb; if (p>pb) break;
#ifdef TOTAL_STAT
    sl2se+=2/(double)p;
#endif
    lo=fblog[i]; sv=mpqs_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }
  pb=mpqs_sievelen;
  for (;;) {
    p=*fb; if (p>pb) break;
#ifdef TOTAL_STAT
    sl2se+=2/(double)p;
#endif
    lo=fblog[i]; sv=mpqs_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }

  for (; i<mpqs_nFB;) {
    p=*fb;
#ifdef TOTAL_STAT
    sl2se+=2/(double)p;
#endif
    lo=fblog[i]; sv=mpqs_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }
#endif
  sv=mpqs_sievearray;
#ifdef TOTAL_STAT
  sieve_events+=sl2se*mpqs_sievelen;
  total_sieve_length+=mpqs_sievelen;
#endif
#ifdef MPQS_ZEIT
zeitA(6);
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
zeitB(6);
#endif
  for (i=mpqs_sievebegink; i<mpqs_nFBk; i++) {
    p=mpqs_FBk[i];
    lo=mpqs_FBk_log[i];
    s1=mpqs_FBk_start[i];
#if 1
    sv=mpqs_sievearray;
    svend=sv+mpqs_sievelen-4*p;
    while (sv<svend) {
      sv[s1]+=lo; sv+=p;
      sv[s1]+=lo; sv+=p;
      sv[s1]+=lo; sv+=p;
      sv[s1]+=lo; sv+=p;
    }
    svend+=(p+p);
    if (sv<svend) {
      sv[s1]+=lo; sv+=p;
      sv[s1]+=lo; sv+=p;
    }
    svend+=p;
    if (sv<svend) {
      sv[s1]+=lo; sv+=p;
    }
    svend+=p;
    if (sv+s1<svend) sv[s1]+=lo;
#else
    while (s1<mpqs_sievelen) { sv[s1]+=lo; s1+=p; }
#endif
  }
}


static ushort mpqs_evaluate()
{
  u32_t *ulsv, *ulsvend, h;
  uchar *sv;
  ushort **rels, buffer[256];
  u32_t i, nmaxsurv;

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
#ifdef MPQS_STAT
  stat_asm_eval++;
#endif
#ifdef ASM_MPQS_EVAL
#ifdef HAVE_XMM_MUL
  if (mpqs_td_begin<mpqs_nFB)
    mpqs_nsurvivors=asm_evaluate0_xmm(ulsv,ulsvend,buffer,nmaxsurv);
  else
    mpqs_nsurvivors=asm_evaluate_xmm(ulsv,ulsvend,buffer,nmaxsurv); 
#else
  if (mpqs_td_begin<mpqs_nFB)
    mpqs_nsurvivors=asm_evaluate0(ulsv,ulsvend,buffer,nmaxsurv);
  else
    mpqs_nsurvivors=asm_evaluate(ulsv,ulsvend,buffer,nmaxsurv);
#endif
#else
  while (ulsv<ulsvend) {
    h=*ulsv;
    if (h&0x80808080) {
      sv=(uchar *)ulsv;
      for (i=0; i<4; i++)
        if (*sv++&0x80)
          buffer[mpqs_nsurvivors++]=(ushort)(sv-mpqs_sievearray-1);
      if (mpqs_nsurvivors>nmaxsurv-8) {
        while (ulsv<ulsvend) *ulsv++=0;
        break;
      }
    }
    *ulsv++=0;
    h=*ulsv;
    if (h&0x80808080) {
      sv=(uchar *)ulsv;
      for (i=0; i<4; i++)
        if (*sv++&0x80)
          buffer[mpqs_nsurvivors++]=(ushort)(sv-mpqs_sievearray-1);
      if (mpqs_nsurvivors>nmaxsurv-8) {
        while (ulsv<ulsvend) *ulsv++=0;
        break;
      }
    }
    *ulsv++=0;
  }
#endif

  for (i=0; i<mpqs_nsurvivors; i++) {
    mpqs_sievearray[buffer[i]]=(uchar)(i+1);
    rels[i][0]=buffer[i];
  }
  return mpqs_nsurvivors;
}


/* change hash_table size to 256*15 ??? */
static inline ushort mpqs_hash_lp(u32_t lp)
{
  u32_t j;
  u32_t lp1;
  ushort lphash, nh;

  lphash=(ushort)(lp&0x00ff); lphash>>=1;
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


static inline int mpqs_hash_rel(long long axb)
{
  u32_t j;
  u32_t hash, hash2, nh;

  hash=(u32_t)(axb&0xffffff);
  hash2=(ushort)(hash>>8); hash&=0xff;
  nh=mpqs_rel_hash_table[hash][0];
  if (nh>=15) return 1;
  for (j=1; j<=nh; j++)
    if (mpqs_rel_hash_table[hash][j]==hash2) return 1;
  mpqs_rel_hash_table[hash][nh+1]=hash2;
  mpqs_rel_hash_table[hash][0]++;
  return 0;
}


static void mpqs_decompose()
{
  u32_t i, j;
  ushort *fb=mpqs_FB, *fbs=mpqs_FB_start;
  uchar *sv, *svend, i1, i2;
  ushort ind, p, s1, s2, nr, d, lpi, r, ii;
  ushort **rels;
  long long axb, *llp;
  u64_t qx, pr, prasm;
  ushort d2, ii2, minus;
  double dbl_qx, dx;
  short x;
  u32_t *ultarg, *ulsrc, lp, ulqx;
  u32_t inv, h, hi, ls1, ls2;
  u32_t ax, ay, az, at, asmh;
  u32_t rr, jj;

#ifdef MPQS_ZEIT
zeitA(11);
#endif
  rels=mpqs_rel_buffer-1;
  for (i=1; i<=mpqs_nsurvivors; i++) rels[i][4]=0;
  i=mpqs_td_begin; fb+=2*i; fbs+=2*i;
#ifdef ASM_MPQS_TDSIEVE
{
  ushort *buf[256];
  u32_t n;

  for (ii=1; ii<=mpqs_nsurvivors; ii++) buf[ii]=rels[ii]+5+rels[ii][4];
  mpqs_sievearray[mpqs_sievelen]=0;
  n=asm_tdsieve(fb,fbs,&(buf[0]),mpqs_nFBk_1+i-1);
  fb+=n; fbs+=n;
  for (ii=1; ii<=mpqs_nsurvivors; ii++)
    rels[ii][4]=(ushort)(buf[ii]-rels[ii]-5);
}
#else
  for (; i<mpqs_nFB; i++) {
    p=*fb; fb+=2; sv=mpqs_sievearray; svend=sv+mpqs_sievelen-2*p;
    s1=*fbs++; s2=*fbs++;
    while (sv<svend) {
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=mpqs_nFBk_1+i; rels[i1][4]++;}
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=mpqs_nFBk_1+i; rels[i2][4]++;}      }
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=mpqs_nFBk_1+i; rels[i1][4]++;}
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=mpqs_nFBk_1+i; rels[i2][4]++;}
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
#endif
#ifdef MPQS_ZEIT
zeitB(11);
#endif
#ifdef TOTAL_STAT
 td_total+=mpqs_nsurvivors*mpqs_td_begin;
#endif
#ifdef MPQS_STAT
 stat_td_cand+=mpqs_nsurvivors;
#endif
 for (i=1; i<=mpqs_nsurvivors; i++) {
   nr=rels[i][4];
   ind=rels[i][0];
   x=(short)(ind)-mpqs_disp;

#ifdef TOTAL_STAT
    tds_total+=nr;
#endif

   axb=mpqs_A*(long long)x+mpqs_B;
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
#ifdef MPQS_STAT
   stat_asm_td++;
#endif

#ifdef MPQS_ZEIT
   zeitA(11);
#endif
#ifdef ASM_MPQS_TD
#ifndef A64_STYLE_TD
   if (asm_td(rels[i],minus,&qx,&ulqx)) {
#ifdef MPQS_ZEIT
     zeitB(11);
#endif
     goto next;
   }
#else
   if((ulqx=asm_td(rels[i],minus,qx))==0) {
#ifdef MPQS_ZEIT
     zeitB(11);
#endif
     goto next;
   }
#endif
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
 zeitB(11);
#endif
        goto next;
      }
    }
    if (pr<4294967296ULL) {
      if (qx>(pr<<32)) {
/*        printf(",");*/
#ifdef MPQS_ZEIT
zeitB(11);
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
    az=ax*inv;
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
 zeitB(11);
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
 zeitB(11);
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
 zeitB(11);
#endif
            goto next;
          }
        }
        if (ulqx==1) break;
      }

#endif

#if 0
    for (j=0; j<nr; j++) {
      ii=rels[i][5+j]; jj=ii-mpqs_nFBk_1;
      p=mpqs_FB[2*jj];
      inv=mpqs_FB_inv[jj];
      while (1) {
        rr=ulqx*inv;
        if (((u64_t)rr*(u64_t)p) & 0xffffffff00000000ULL) {
/*          if (ulqx%(u32_t)p==0) complain("%u %u %u %u y",rr,p,ulqx,inv);*/
          break;
        }
/*        if (ulqx%(u32_t)p) complain("z");*/
        if (rels[i][4]<MPQS_TD_MAX_NDIV) {
          rels[i][5+rels[i][4]]=ii; rels[i][4]++;
 /*         ulqx/=(u32_t)p;*/
          ulqx=rr;
        } else {
#ifdef MPQS_ZEIT
 zeitB(11);
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
#ifdef MPQS_ZEIT1
 zeitB(11);
#endif
      goto next;
    }
#ifdef MPQS_ZEIT1
zeitB(11);
#endif
    if (ulqx==1) {
#ifdef MPQS_STAT
      stat_td_surv++;
#endif
      llp=(long long *)(mpqs_relations[mpqs_nrels]);
      *llp=axb;
      nr=rels[i][4];
#ifdef USE_MEMCPY
      memcpy(mpqs_relations[mpqs_nrels]+4,rels[i]+4,(1+nr)*sizeof(ushort));
#else
      for (j=0; j<=nr; j++)
        mpqs_relations[mpqs_nrels][j+4]=rels[i][j+4];
#endif
      mpqs_nrels++;
      goto next;
    }
    if (ulqx<1048576) { /*printf("%d ",x);*/  /* !!! */
      if (ulqx<=mpqs_pmax) {
        gmp_fprintf(stderr,"error in mpqs_decompose for %Zd: %lld %u\n",
                    mpqs_N,qx,ulqx);
        goto next;
//        complain("dec.2\n");
      }
      lp=ulqx;
      lpi=mpqs_hash_lp(lp);
      if (lpi==MPQS_HASH_OVERFLOW) goto next;
#ifdef MPQS_STAT
      stat_td_surv++;
#endif
      if (lpi<mpqs_nprels) {
        if (mpqs_lp[lpi]!=lp) Schlendrian("lp");
        llp=(long long *)(mpqs_comb_rels[mpqs_ncrels]);
        *llp=axb;
        llp=(long long *)(mpqs_part_rels[lpi]);
        axb=*llp;
        llp=(long long *)(mpqs_comb_rels[mpqs_ncrels]); llp++;
        *llp=axb;
        nr=rels[i][4];
#ifdef USE_MEMCPY
        memcpy(mpqs_comb_rels[mpqs_ncrels]+9,rels[i]+5,nr*sizeof(ushort));
        memcpy(mpqs_comb_rels[mpqs_ncrels]+9+nr,mpqs_part_rels[lpi]+5,
               mpqs_part_rels[lpi][4]*sizeof(ushort));
#else
        for (j=0; j<nr; j++)
          mpqs_comb_rels[mpqs_ncrels][j+9]=rels[i][j+5];
        for (j=0; j<mpqs_part_rels[lpi][4]; j++)
          mpqs_comb_rels[mpqs_ncrels][j+9+nr]=mpqs_part_rels[lpi][j+5];
#endif
        mpqs_comb_rels[mpqs_ncrels][8]=nr+256*mpqs_part_rels[lpi][4];
        mpqs_ncrels++;
      } else {
        llp=(long long *)(mpqs_part_rels[mpqs_nprels]);
        *llp=axb;
        nr=rels[i][4];
#ifdef USE_MEMCPY
        memcpy(mpqs_part_rels[mpqs_nprels]+4,rels[i]+4,
               (1+nr)*sizeof(ushort));
#else
        for (j=0; j<=nr; j++)
          mpqs_part_rels[mpqs_nprels][j+4]=rels[i][j+4];
#endif
        mpqs_nprels++;
      }
    } /*else stat_counter0++;*/ /*else printf("<");*/
next: ;
  }
  if (mpqs_nsp<mpqs_nrels+mpqs_ncrels)
    mpqs_excess=mpqs_nrels+mpqs_ncrels-mpqs_nsp;
  else mpqs_excess=0;
}


/*
order of the primes:
-1, primes in mpqs_FBk, primes in mpqs_FB, primes in mpqs_Adiv_total
*/

#ifdef ASM_MPQS_GAUSS
u32_t mpqs_gauss_mask[]  =
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

static int mpqs_matrix_init()
{
  u32_t i, j, n32;
  ushort nr, ii, pi;
  u32_t mask;

  mpqs_gauss_n=mpqs_nrels+mpqs_ncrels;
  if (mpqs_gauss_n>=MPQS_GAUSS_MAX) {
    return -1;
  }
  mpqs_gauss_m=mpqs_nsp;
  if (mpqs_gauss_m>=mpqs_gauss_n) Schlendrian("gauss: no excess\n");
/* build matrix */
  mpqs_gauss_n32=(mpqs_gauss_n+31)/32;
  if (mpqs_gauss_n32&1) mpqs_gauss_n32++;
  for (i=1; i<mpqs_gauss_m; i++)
    mpqs_gauss_row[i]=mpqs_gauss_row[i-1]+mpqs_gauss_n32;
  for (i=0; i<mpqs_gauss_m; i++)
    for (j=0; j<mpqs_gauss_n32; j++)
      mpqs_gauss_row[i][j]=0;
#if 0
  for (i=0; i<mpqs_nrels; i++) {
    nr=mpqs_relations[i][4];
    ii=i; mask=mpqs_gauss_mask[ii%32];
    for (j=0; j<nr; j++) {
      pi=mpqs_relations[ii][5+j];
      mpqs_gauss_row[pi][ii/32]^=mask;
    }
  }
  for (i=0; i<mpqs_ncrels; i++) {
    nr=mpqs_comb_rels[i][8]; nr=(nr&255)+(nr>>8);
    ii=i+mpqs_nrels; mask=mpqs_gauss_mask[ii%32];
    for (j=0; j<nr; j++) {
      pi=mpqs_comb_rels[i][9+j];
      mpqs_gauss_row[pi][ii/32]^=mask;
    }
  }
#else
  for (i=0; i<mpqs_nrels; i++) {
    nr=mpqs_relations[i][4];
    ii=i+mpqs_ncrels; mask=mpqs_gauss_mask[ii%32];
    for (j=0; j<nr; j++) {
      pi=mpqs_relations[i][5+j];
      mpqs_gauss_row[pi][ii/32]^=mask;
    }
  }
  for (i=0; i<mpqs_ncrels; i++) {
    nr=mpqs_comb_rels[i][8]; nr=(nr&255)+(nr>>8);
    ii=i; mask=mpqs_gauss_mask[ii%32];
    for (j=0; j<nr; j++) {
      pi=mpqs_comb_rels[i][9+j];
      mpqs_gauss_row[pi][ii/32]^=mask;
    }
  }
#endif
  for (i=0; i<mpqs_gauss_m+1; i++) mpqs_gauss_c[i]=-1;
  mpqs_gauss_k=mpqs_gauss_n;
  for (i=0; i<mpqs_gauss_n; i++) mpqs_gauss_d[i]=-1;
  mpqs_gauss_j=0;
  return 0;
}
#else

static int mpqs_matrix_init()
{
  u32_t i, j, n64, i64;
  ushort nr, ii, pi;
  u64_t mask;

  mpqs_gauss_n=mpqs_nrels+mpqs_ncrels;
  if (mpqs_gauss_n>=MPQS_GAUSS_MAX) {
    return -1;
  }
  mpqs_gauss_m=mpqs_nsp;
  if (mpqs_gauss_m>=mpqs_gauss_n) Schlendrian("gauss: no excess\n");
/* build matrix */
  mpqs_gauss_n64=(mpqs_gauss_n+63)/64;
mpqs_gauss_n32=2*mpqs_gauss_n64;
  for (i=1; i<mpqs_gauss_m; i++)
    mpqs_gauss_row64[i]=mpqs_gauss_row64[i-1]+mpqs_gauss_n64;
for (i=0; i<mpqs_gauss_m; i++)
  mpqs_gauss_row[i]=(u32_t *)(mpqs_gauss_row64[i]);
  for (i=0; i<mpqs_gauss_m; i++)
    for (j=0; j<mpqs_gauss_n64; j++)
      mpqs_gauss_row64[i][j]=0;
#if 0
  for (i=0; i<mpqs_nrels; i++) {
    nr=mpqs_relations[i][4];
    ii=i; mask=1ULL<<(i&63); i64=ii>>6;
    for (j=0; j<nr; j++) {
      pi=mpqs_relations[ii][5+j];
      mpqs_gauss_row64[pi][i64]^=mask;
    }
  }
  for (i=0; i<mpqs_ncrels; i++) {
    nr=mpqs_comb_rels[i][8]; nr=(nr&255)+(nr>>8);
    ii=i+mpqs_nrels; mask=1ULL<<(ii&63); i64=ii>>6;
    for (j=0; j<nr; j++) {
      pi=mpqs_comb_rels[i][9+j];
      mpqs_gauss_row64[pi][i64]^=mask;
    }
  }
#else
  for (i=0; i<mpqs_nrels; i++) {
    nr=mpqs_relations[i][4];
    ii=i+mpqs_ncrels; mask=1ULL<<(ii&63); i64=ii>>6;
    for (j=0; j<nr; j++) {
      pi=mpqs_relations[i][5+j];
      mpqs_gauss_row64[pi][i64]^=mask;
    }
  }
  for (i=0; i<mpqs_ncrels; i++) {
    nr=mpqs_comb_rels[i][8]; nr=(nr&255)+(nr>>8);
    ii=i; mask=1ULL<<(ii&63); i64=ii>>6;
    for (j=0; j<nr; j++) {
      pi=mpqs_comb_rels[i][9+j];
      mpqs_gauss_row64[pi][i64]^=mask;
    }
  }
#endif
  for (i=0; i<mpqs_gauss_m+1; i++) mpqs_gauss_c[i]=-1;
  mpqs_gauss_k=mpqs_gauss_n;
  for (i=0; i<mpqs_gauss_n; i++) mpqs_gauss_d[i]=-1;
  mpqs_gauss_j=0;
  return 0;
}
#endif


#ifdef ASM_MPQS_GAUSS
static int mpqs_rowechelon()
{
// SMJS Was: return;
  return 0;
}

static int mpqs_matrix()
{
  u32_t i, j, l, t, k32, len;
  u32_t mask;
  u64_t h;
 
#ifdef MPQS_ZEIT
zeitA(8);
#endif
/* solve matrix */
#ifdef ASM_MPQS_GAUSS
  asm_gauss();
  if (!(mpqs_gauss_k+1)) {
#ifdef MPQS_ZEIT
zeitB(8);
#endif
    return 0;
  }
#else
  while (1) {
    if (!mpqs_gauss_k) {
#ifdef MPQS_ZEIT
zeitB(8);
#endif
      return 0;
    }
    mpqs_gauss_k--;
    mask=(1UL<<(mpqs_gauss_k&31)); k32=mpqs_gauss_k/32;
    j=mpqs_gauss_j;
    while ((j<mpqs_gauss_m) && (mpqs_gauss_c[j]>=0 ||
       ((mpqs_gauss_row[j][k32] & mask)==0))) j++;
    if (j==mpqs_gauss_m) { mpqs_gauss_d[mpqs_gauss_k]=-1; break; }
/* exchange j and mpqs_gauss_j */
    if (j!=mpqs_gauss_j) {
      for (l=0; l<mpqs_gauss_n32; l++) {
        h=mpqs_gauss_row[mpqs_gauss_j][l];
        mpqs_gauss_row[mpqs_gauss_j][l]=mpqs_gauss_row[j][l];
        mpqs_gauss_row[j][l]=h;
      }
      j=mpqs_gauss_j;
    }
    mpqs_gauss_d[mpqs_gauss_k]=j; mpqs_gauss_c[j]=mpqs_gauss_k;
    while (mpqs_gauss_c[mpqs_gauss_j]>=0) mpqs_gauss_j++;
    for (t=0; t<j; t++)
      if (mpqs_gauss_row[t][k32] & mask)
        for (l=0; l<=k32; l++)
          mpqs_gauss_row[t][l]^=mpqs_gauss_row[j][l];
    for (t=j+1; t<mpqs_gauss_m; t++)
      if (mpqs_gauss_row[t][k32] & mask)
        for (l=0; l<=k32; l++)
          mpqs_gauss_row[t][l]^=mpqs_gauss_row[j][l];
  }
#endif
  mask=1ULL<<(mpqs_gauss_k&31);
  k32=mpqs_gauss_k/32;

  for (i=0; i<mpqs_gauss_n; i++) mpqs_sol[i]=0;
  if (mpqs_gauss_d[mpqs_gauss_k]!=-1) Schlendrian("gauss1 %u\n",mpqs_gauss_k);
  for (i=mpqs_gauss_k+1; i<mpqs_gauss_n; i++)
    if (mpqs_gauss_d[i]!=-1)
      if (mpqs_gauss_row[mpqs_gauss_d[i]][k32] & mask)
        mpqs_sol[i]=1;
  mpqs_sol[mpqs_gauss_k]=1;
#ifdef MPQS_ZEIT
zeitB(8);
#endif
  return 1;
}
#else

static int mpqs_rowechelon()
{
  u32_t i, j, k, l, t, k64, col, row, r0, c0, zz;
  u64_t tmp, tab[16];
  uchar ucm[1024], mm[4], m, z;

//zeitA(8);
  col=0; row=0;
  while (1) {
    if (row>=mpqs_gauss_m) break;
    if (col>=mpqs_gauss_n) break;

    k64=col>>6;
    r0=row; c0=col;
    for (i=0; i<4; i++,col++) {
//zeitA(26);
      for (j=row; j<mpqs_gauss_m; j++) {
        m=(uchar)(0xf&mpqs_gauss_row64[j][k64]>>(c0&63));
        if (!m) continue;
        for (k=0,z=1; k<i; k++,z+=z) if (m&z) m^=mm[k];
        if (m&z) break;
      }
//zeitB(26);
      if (j>=mpqs_gauss_m) {
        mpqs_gauss_d[col]=-1;
        continue;
       }
      if (j>row) {
        for (l=k64; l<mpqs_gauss_n64; l++) {
          tmp=mpqs_gauss_row64[j][l];
          mpqs_gauss_row64[j][l]=mpqs_gauss_row64[row][l];
          mpqs_gauss_row64[row][l]=tmp;
        }
      }
      m=(uchar)(mpqs_gauss_row64[row][k64]>>(c0&63));
      for (k=0,z=1; k<i; k++,z+=z) {
        if (m&z) {
          m^=mm[k];
          for (l=k64; l<mpqs_gauss_n64; l++)
            mpqs_gauss_row64[row][l]^=
                mpqs_gauss_row64[mpqs_gauss_d[c0+k]][l];
        }
      }
      mm[i]=m;
      for (k=0; k<i; k++) {
        if (mpqs_gauss_d[c0+k]==-1) continue;
        if (mm[k]&z) {
          mm[k]^=mm[i];
          for (l=k64; l<mpqs_gauss_n64; l++)
            mpqs_gauss_row64[mpqs_gauss_d[c0+k]][l]^=
                mpqs_gauss_row64[row][l];
        }
      }
      mpqs_gauss_d[col]=row;
      mpqs_gauss_c[row]=col;
      row++;
      if (row>=mpqs_gauss_m) break;
    }
//zeitA(25);
    for (j=0; j<mpqs_gauss_m; j++)
      ucm[j]=(uchar)((mpqs_gauss_row64[j][k64]>>(c0&63))&0xf);
    for (j=r0; j<row; j++) ucm[j]=0;
    for (l=k64; l<mpqs_gauss_n64; l++) {
      tab[0]=0ULL;
      for (j=0,zz=1; j<4; j++,zz+=zz) {
        if (mpqs_gauss_d[c0+j]==-1) tab[zz]=0ULL;
        else tab[zz]=mpqs_gauss_row64[mpqs_gauss_d[c0+j]][l];
        for (k=1; k<zz; k++) tab[zz+k]=tab[k]^tab[zz];
      }
      for (t=0; t<mpqs_gauss_m; t++)
        mpqs_gauss_row64[t][l]^=tab[ucm[t]];
    }
//zeitB(25);
  }
  mpqs_gauss_k=0;
//zeitB(8);
}


static int mpqs_matrix()
{
  u32_t i, k64;
  int notzero;
  u64_t mask;

  for (i=0; i<mpqs_gauss_n; i++) mpqs_sol[i]=0;
  for (; mpqs_gauss_k<mpqs_gauss_n; mpqs_gauss_k++) {
    if (mpqs_gauss_d[mpqs_gauss_k]>=0) continue;
    notzero=0;
    mask=1ULL<<(mpqs_gauss_k&63); k64=mpqs_gauss_k>>6;
    for (i=0; i<mpqs_gauss_k; i++)
      if (mpqs_gauss_d[i]!=-1)
        if (mpqs_gauss_row64[mpqs_gauss_d[i]][k64] & mask) {
          mpqs_sol[i]=1;
          notzero=1;
        }

    if (notzero) {
      mpqs_sol[mpqs_gauss_k]=1;
      mpqs_gauss_k++;
      return 1;
    }
  }
  return 0;
}
#endif

static u32_t mpqs_is_power(mpz_t op) /* very sloppy */
{
  u32_t res;

  for (res=2; res<12; res++)
    if (mpz_root(mpqs_gmp_help,op,(ulong)res))
      return res;
  return 0;
}


static int mpqs_final()
{
  u32_t i, j, k;
  ushort nr, nr1, nr2;
  u32_t p, lp, prod1, prod2;
  u64_t *up;
  int split;
  size_t modthres;

  modthres=mpz_sizeinbase(mpqs_N,2)*4;
  if (!mpqs_nfactors) {
    mpqs_nfactors=1;
    mpz_set(mpqs_factors[0],mpqs_N);
    mpqs_f_status[0]=0;
  }
  mpqs_rowechelon();
  while (mpqs_matrix()) {
#ifdef MPQS_STAT
stat_mpqs_ntrials++;
#endif
    memset(mpqs_exp,0,mpqs_nsp*sizeof(short));
    mpz_set_ui(mpqs_sq1,1); prod1=1;
    mpz_set_ui(mpqs_sq2,1); prod2=1;

    for (j=0; j<mpqs_nrels; j++)
      if (mpqs_sol[mpqs_ncrels+j]) {
        nr=mpqs_relations[j][4];
        up=(u64_t *)(mpqs_relations[j]); /*printf("Q: %llu  ",*up);*/
        mpz_set_ull(mpqs_gmp_help,*up);
        if (j&1) {   /* 'random' distribution of full relations */
          for (k=0; k<nr; k++) mpqs_exp[mpqs_relations[j][5+k]]--;
          mpz_mul(mpqs_sq1,mpqs_sq1,mpqs_gmp_help);
          if (mpz_sizeinbase(mpqs_sq1,2)>modthres)
            mpz_mod(mpqs_sq1,mpqs_sq1,mpqs_N);
        } else {
          for (k=0; k<nr; k++) mpqs_exp[mpqs_relations[j][5+k]]++;
          mpz_mul(mpqs_sq2,mpqs_sq2,mpqs_gmp_help);
          if (mpz_sizeinbase(mpqs_sq2,2)>modthres)
            mpz_mod(mpqs_sq2,mpqs_sq2,mpqs_N);
        }
      }

    for (j=0; j<mpqs_ncrels; j++)
      if (mpqs_sol[j]) {
        nr=mpqs_comb_rels[j][8];
        nr1=nr&255; nr2=nr>>8;
        for (k=0; k<nr1; k++) mpqs_exp[mpqs_comb_rels[j][9+k]]++;
        for (k=0; k<nr2; k++) mpqs_exp[mpqs_comb_rels[j][9+nr1+k]]--;

        up=(u64_t *)(mpqs_comb_rels[j]); /*printf("Q: %llu  ",*up);*/
        mpz_set_ull(mpqs_gmp_help,*up);
        mpz_mul(mpqs_sq2,mpqs_sq2,mpqs_gmp_help);
        if (mpz_sizeinbase(mpqs_sq2,2)>modthres)
          mpz_mod(mpqs_sq2,mpqs_sq2,mpqs_N);
        up++;
        mpz_set_ull(mpqs_gmp_help,*up);
        mpz_mul(mpqs_sq1,mpqs_sq1,mpqs_gmp_help);
        if (mpz_sizeinbase(mpqs_sq1,2)>modthres)
          mpz_mod(mpqs_sq1,mpqs_sq1,mpqs_N);
      }

    for (j=0; j<mpqs_nsp; j++)
      if (mpqs_exp[j]&1) Schlendrian("final.odd\n");
      else mpqs_exp[j]>>=1;

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
              mpz_mul_ui(mpqs_sq1,mpqs_sq1,(ulong)prod1);
              if (mpz_sizeinbase(mpqs_sq1,2)>modthres) {
                mpz_mod(mpqs_sq1,mpqs_sq1,mpqs_N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              }
              prod1=1;
            }
          }
        } else {
          for (k=0; k<-mpqs_exp[j]; k++) {
            prod2*=p;
            if (prod2&0xfff00000) {  /* p<4096 */
              mpz_mul_ui(mpqs_sq2,mpqs_sq2,(ulong)prod2);
              if (mpz_sizeinbase(mpqs_sq2,2)>modthres) {
                mpz_mod(mpqs_sq2,mpqs_sq2,mpqs_N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              }
              prod2=1;
            }
          }
        }
      }
    if (prod1>1) mpz_mul_ui(mpqs_sq1,mpqs_sq1,(ulong)prod1);
    if (prod2>1) mpz_mul_ui(mpqs_sq2,mpqs_sq2,(ulong)prod2);

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
#if 1
        mpqs_f_status[j]=psp(mpqs_factors[j]);
        mpqs_f_status[mpqs_nfactors]=psp(mpqs_factors[mpqs_nfactors]);
#else
        mpqs_f_status[j]=mpz_probab_prime_p1(mpqs_factors[j],1);
        mpqs_f_status[mpqs_nfactors]=
          mpz_probab_prime_p1(mpqs_factors[mpqs_nfactors],1);
#endif
        mpqs_nfactors++;
        break;
      }

    split=1;
    for (j=0; j<mpqs_nfactors; j++) if (mpqs_f_status[j]==0) split=0;
    if (split) return 1;
    if (mpqs_nfactors>=MPQS_MAX_FACTORS)
      Schlendrian("final: too many factors\n");
  }
  for (j=0; j<mpqs_nfactors; j++)
    while (mpqs_f_status[j]==0) {
      if ((k=mpqs_is_power(mpqs_factors[j]))) {
        if (mpqs_nfactors+k-1>MPQS_MAX_FACTORS)
          Schlendrian("final: too many factors\n");
        if (!mpz_root(mpqs_factors[j],mpqs_factors[j],(ulong)k))
          Schlendrian("error in power detection\n");
        mpqs_f_status[j]=psp(mpqs_factors[j]);
        for (k--; k; k--) {
          mpz_set(mpqs_factors[mpqs_nfactors],mpqs_factors[j]);
          mpqs_f_status[mpqs_nfactors]=mpqs_f_status[j];
          mpqs_nfactors++;
        }
      } else return 0;
    }
  if (j>=mpqs_nfactors) return 1;
  return 0;
}



static int mpqs_factor0(mpz_t N, size_t max_bits, mpz_t **factors, int retry)
{
  size_t nbits;
  int err;
  u32_t ev, i;

#ifdef MPQS_ZEIT
zeitA(9);
#endif
#ifdef MPQS_STAT
  stat_mpqs_nsieves=0; stat_mpqs_nsurvivors=0; stat_mpqs_ntrials=0; stat_mpqs_ndiv=0;
  if (retry) stat_retry++;
#endif
  if (!mpqs_isinit) mpqs_init();
#ifdef MPQS_ZEIT
zeitB(9);
#endif
  nbits=mpz_sizeinbase(N,2);
  if (nbits>96) {
    fprintf(stderr,"warning: mpqs called with >96 Bit\n");
    return -2;
  }
#ifdef MPQS_ZEIT
zeitA(1);
#endif
//zeitA(20);
  mpz_set(mpqs_N,N);
//zeitB(20);
//zeitA(21);
  mpqs_choose_multiplier();
//zeitB(21);
//zeitA(22);
  mpqs_choose_parameter(retry);
//zeitB(22);
//zeitA(23);
  mpqs_generate_FB();
//zeitB(23);
//zeitA(24);
  if ((err=mpqs_SI_init())<1) {
    printf("%d ",err);
    fprintf(stderr,"warning: mpqs: error in self-initialization ");
    mpz_out_str(stderr,10,N); printf("\n");
#ifdef MPQS_ZEIT
zeitB(1);
#endif
    return -3;
  }
//zeitB(24);
#ifdef MPQS_ZEIT
zeitB(1);
#endif
  while (1) { // printf("%u,%u; ",stat_mpqs_nsieves,mpqs_nrels);
#ifdef MPQS_ZEIT
zeitA(2);
#endif
    if (!mpqs_next_pol()) {
      if (retry) {
        fprintf(stderr,"warning: not enough polynomials in mpqs for ");
        mpz_out_str(stderr,10,N); fprintf(stderr,"\n");
      }
#ifdef MPQS_ZEIT
zeitB(2);
#endif
      return -4;
    }
#ifdef MPQS_ZEIT
zeitB(2); zeitA(3);
#endif
    mpqs_sieve();
#ifdef MPQS_STAT
stat_mpqs_nsieves++;
#endif
#ifdef MPQS_ZEIT
zeitB(3);
zeitA(4);
#endif
    ev=mpqs_evaluate();
#ifdef MPQS_STAT
stat_mpqs_nsurvivors+=ev; // printf("%u %lld  \n",ev,mpqs_C);
#endif
#ifdef MPQS_ZEIT
zeitB(4);
#endif
    if (ev) {
#ifdef MPQS_ZEIT
zeitA(5);
#endif
      mpqs_decompose(); /*test_rels();*/
#ifdef MPQS_ZEIT
zeitB(5);
#endif
      if (mpqs_nsurvivors && (mpqs_excess>MPQS_MIN_EXCESS)) {
#ifdef MPQS_ZEIT
zeitA(11);
#endif
        if (mpqs_matrix_init()) return -6;
#ifdef MPQS_ZEIT
zeitB(11);
zeitA(10);
#endif
        if (mpqs_final()) {
#ifdef MPQS_ZEIT
zeitB(10);
#endif
#ifdef MPQS_STAT
// printf("Stat: %u %u %u %u %u %u %u \n",stat_mpqs_nsieves,stat_mpqs_nsurvivors,stat_mpqs_ntrials,stat_mpqs_ndiv,mpqs_nrels,mpqs_nprels,mpqs_ncrels);
stat_ff+=mpqs_nrels;
stat_pf+=(mpqs_nprels+mpqs_ncrels);
stat_comb+=mpqs_ncrels;
#endif
          for (i=0; i<mpqs_nfactors; i++)
            if (mpz_sizeinbase(mpqs_factors[i],2)>max_bits) return 0;
          *factors=mpqs_factors;
          return mpqs_nfactors;
        }
#ifdef MPQS_ZEIT
zeitB(10);
#endif
      }
    }
    if (mpqs_nrels>(MPQS_MAX_NRELS-MPQS_MIN_RELBUFFER)) {
      fprintf(stderr,"warning: too many relations in mpqs\n");
      return -5;
    }
  }
  return -1;
}


int mpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
  int err;

//gmp_printf("%Zd\n",N);
  err=mpqs_factor0(N,max_bits,factors,0);
  if (err==-4) err=mpqs_factor0(N,max_bits,factors,1);
  if (err==-4) err=mpqs_factor0(N,max_bits,factors,2);
  return err;
}


void
mpqs_total_stat()
{
#ifdef TOTAL_STAT
  printf("%u -> %u trial division candidates\n",
	 stat_td_cand,stat_td_surv);
  printf("%u divisions by sieving, true trial divisions: %u\n",
	 tds_total,td_total);
  printf("sieve events: %.3g, total sieve length: %.3g\n",
	 sieve_events,total_sieve_length);
  printf("total fb primes: %u sum of ulqx after td: %u\n",
	 total_nfbp,total_ulqx);
#endif
}
