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
#include <gmp.h>
#include "asm/siever-config.h"
#include "if.h"
#include "gmp-aux.h"
#ifdef VARIANT2
#include "asm/mpqs-config2.h"
#else
#include "asm/mpqs-config.h"
#endif


#ifdef MPQS3_STAT
extern u64_t stat_asm_eval, stat_asm_td;
extern u64_t stat_td_cand,stat_td_surv;
extern u64_t stat_ff,stat_pf,stat_comb,stat_pp;
extern u32_t stat_asm_div, stat_final_mulmod;
extern u32_t stat_retry;
extern u32_t stat_size[15];
#endif

/* 2^-64 */
#define DBL_CONVERT   5.421010862427522170037264004e-20

#define uchar unsigned char

#define MPQS3_NULONG           3

#define L1BITS               16

#define MPQS3_MAX_FBSIZE      1024
#define MPQS3_MIN_EXCESS      10
#define MPQS3_MAX_ADIV_ALL    12   /* at most 16 */
#define MPQS3_MAX_FBXSIZE     (MPQS3_MAX_FBSIZE+32)
#define MPQS3_SIEVELEN        (1<<L1BITS)
#define MPQS3_REL_ENTRIES     36
#define MPQS3_TD_MAX_NDIV     (MPQS3_REL_ENTRIES-7)
#define MPQS3_MAX_NRELS       512  /* ??? */
#define MPQS3_MAX_NPRELS      8192  /* ??? */
#define MPQS3_MAX_NCRELS      1024  /* ??? */
#define MPQS3_MAX_NRELS_BUF   256  /* ??? */
#define MPQS3_MIN_RELBUFFER   16
#define MPQS3_GAUSS_MAX       1024
#define MPQS3_MAX_FACTORS     6
#define MPQS3_MAX_NPRIMES     2048   /* largest prime MUST be <16384 (final) */
#define MPQS3_HASH_OVERFLOW   65535
#define MPQS3_MULT_NTESTPR    25       /* prime(MPQS3_MULT_NTESTPR)<2^8 */
#define MPQS3_MULT_NCAND      8        /* largest candidate<2^8 */
#define MPQS3_FB_MAXPRIME     16384     /* used in nextpol, sieve and final */
#define MPQS3_SI_NVALUES      4
#define MPQS3_MAX_TINYPROD    128

#define MPQS3_NHASHBITS       8
#define MPQS3_NHASH           (1<<(MPQS3_NHASHBITS-1))
#define MPQS3_NHASHMASK       (2*MPQS3_NHASH-1)
#define MPQS3_HASH_LEN        64
#define MPQS3_RELHASH_LEN     512

#define MPQS3_NRELHASH        256


/* common with asm functions */
ushort mpqs3_nFBk_1;
ushort mpqs3_td_begin, mpqs3_sievebegin, mpqs3_sievebegink;
ushort mpqs3_FB_inv_info[4*MPQS3_MAX_NPRIMES];
ushort mpqs3_nFBk, mpqs3_FBk[3];  /* 3*5*7*11 too big */
uchar mpqs3_FB_log[MPQS3_MAX_FBSIZE], mpqs3_FBk_log[3];
ushort mpqs3_nFB, mpqs3_nFBx, mpqs3_nAdiv_total;
ushort mpqs3_Adiv_all[MPQS3_MAX_ADIV_ALL];
u32_t mpqs3_FBk_inv[3], mpqs3_FB_A_inv[MPQS3_MAX_ADIV_ALL];

ushort mpqs3_Adiv_disp[MPQS3_MAX_ADIV_ALL];
ushort mpqs3_FB_mm_inv[MPQS3_MAX_FBXSIZE], mpqs3_FB_disp[MPQS3_MAX_FBXSIZE];
ushort mpqs3_FB_np_p[2*MPQS3_MAX_FBXSIZE];
ushort mpqs3_FB_np_px[2*MPQS3_MAX_FBXSIZE];
/* end */

static int mpqs3_isinit=0;
static mpz_t mpqs3_gmp_help, mpqs3_dummy, mpqs3_help;
static mpz_t mpqs3_N, mpqs3_kN, mpqs3_factors[MPQS3_MAX_FACTORS];
static uchar mpqs3_f_status[MPQS3_MAX_FACTORS];
static u32_t mpqs3_multiplier, mpqs3_nfactors;
static double mpqs3_complexity;
static ushort mpqs3_prime_table[MPQS3_MAX_NPRIMES];
static ushort mpqs3_prime_sqrthelp[MPQS3_MAX_NPRIMES];
static ushort mpqs3_nAdiv, mpqs3_accept;
static uchar mpqs3_Adiv_log[MPQS3_MAX_ADIV_ALL];;
static ushort mpqs3_pmax;
static ushort mpqs3_Adiv[MPQS3_MAX_ADIV_ALL];
static ushort mpqs3_Adiv_sqrt[MPQS3_MAX_ADIV_ALL], mpqs3_Adiv_all_sqrt[MPQS3_MAX_ADIV_ALL];
static ushort mpqs3_Adiv_start1[MPQS3_MAX_ADIV_ALL], mpqs3_Adiv_start2[MPQS3_MAX_ADIV_ALL];
static uchar mpqs3_Adiv_active[MPQS3_MAX_ADIV_ALL];
static ushort mpqs3_A_mask[(1<<MPQS3_MAX_ADIV_ALL)];
static ushort mpqs3_nA;
static short mpqs3_A_index, mpqs3_B_index;
static ushort mpqs3_SI_inv_table[MPQS3_MAX_ADIV_ALL][MPQS3_MAX_FBSIZE];
static ushort mpqs3_Adiv_SI_inv_table[MPQS3_MAX_ADIV_ALL][MPQS3_MAX_ADIV_ALL];
static ushort mpqs3_SI_add[MPQS3_MAX_ADIV_ALL][MPQS3_MAX_FBXSIZE];
static ushort mpqs3_SIk_add[MPQS3_MAX_ADIV_ALL][3];
static ushort mpqs3_Adiv_SI_add[MPQS3_MAX_ADIV_ALL][MPQS3_MAX_ADIV_ALL];
static u32_t mpqs3_Adiv_all_sq_mod_p2[MPQS3_MAX_ADIV_ALL];
static ushort mpqs3_Adiv_all_sq_inv[MPQS3_MAX_ADIV_ALL];

static double mpqs3_A_dbl, mpqs3_B_dbl, mpqs3_2B_dbl, mpqs3_C_dbl;
static double mpqs3_Bi_dbl[MPQS3_MAX_ADIV_ALL];
static u64_t mpqs3_A_64, mpqs3_B_64, mpqs3_2B_64, mpqs3_C_64;
static u64_t mpqs3_Bi_64[MPQS3_MAX_ADIV_ALL];


static u32_t mpqs3_disp;
static ushort mpqs3_nsurvivors;
static u32_t mpqs3_nrels, mpqs3_nlp, mpqs3_nsp, mpqs3_excess;
static uchar mpqs3_logbound;
static ushort mpqs3_nprels, mpqs3_ncrels;
static ushort **mpqs3_relations, **mpqs3_part_rels;
static ushort *mpqs3_comb_index;
static ushort **mpqs3_rel_buffer;
static ushort mpqs3_rel_hash_table[MPQS3_NRELHASH][MPQS3_RELHASH_LEN];
static ushort mpqs3_hash_table[MPQS3_NHASH][MPQS3_HASH_LEN];
static u32_t mpqs3_hash_index[MPQS3_NHASH][MPQS3_HASH_LEN];
static u32_t mpqs3_lp[MPQS3_NHASH*(MPQS3_HASH_LEN-1)];
static u64_t *mpqs3_gauss_row64;
static uchar mpqs3_sol[MPQS3_GAUSS_MAX];
static short mpqs3_gauss_d[MPQS3_GAUSS_MAX];
static u32_t mpqs3_gauss_m, mpqs3_gauss_n, mpqs3_gauss_n64, mpqs3_gauss_k;
static short mpqs3_exp[128*16+MPQS3_MAX_FBSIZE+3+MPQS3_MAX_ADIV_ALL]; /* ??? */
static mpz_t mpqs3_sq1, mpqs3_sq2;
#if 0
static ushort mpqs_zero1, mpqs_zero2;
#endif
static double mpqs3_kN_dbl, mpqs3_maxval_dbl, mpqs3_A_inv_dbl;
static u64_t mpqs3_kN_64, mpqs3_A_inv_64;

ushort mpqs3_FB[2*MPQS3_MAX_FBXSIZE], mpqs3_FB0[MPQS3_MAX_FBSIZE+2*256];
ushort mpqs3_FB_start[2*MPQS3_MAX_FBXSIZE], mpqs3_FBk_start[3];
u32_t mpqs3_sievelen;
u32_t mpqs3_FB_inv[MPQS3_MAX_FBXSIZE];
uchar *mpqs3_sievearray, *mpqs3_tinyarray;
u32_t mpqs3_tiny_prod, mpqs3_ntiny;
ushort mpqs3_sinit_tab[8*MPQS3_SI_NVALUES];

uchar mpqs3_256_inv_table[128]={
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


/* 0x00, ..., 0xff sorted w.r.t. number of ones: */
uchar mpqs3_A_table[256]={
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

u32_t mpqs3_A_table_n[10]={
0, 1, 9, 37, 93, 163, 219, 247, 255, 256
};

#ifdef MPQS3_STAT
u32_t stat_mpqs_nsieves, stat_mpqs_nsurvivors, stat_mpqs_ntrials, stat_mpqs_ndiv;
#endif

#include "asm/mpqs3arith.c"

void mpqs3_convert(u32_t *rop, double op_dbl, u64_t op_64)
{
  double d;

  rop[0]=(u32_t)op_64;
  rop[1]=(u32_t)(op_64>>32);
  d=op_dbl-(double)op_64;
  d*=DBL_CONVERT;
  if (d<0.) d+=4294967296.;
  rop[2]=(u32_t)(d+0.5);
}

ushort mpqs3_dbl64mod16(double op_dbl, u64_t op_64, ushort p)
{
  u32_t h[3];
  ushort res;

  if (op_dbl>=0.) {
    mpqs3_convert(h,op_dbl,op_64);
    return (ushort)(mpqs3_mod3(h,(u32_t)p));
  }
  mpqs3_convert(h,-op_dbl,-op_64);
  res=(ushort)(mpqs3_mod3(h,(u32_t)p));
  if (res) res=p-res;
  return res;
}

#if 1
u32_t mpqs3_dbl64mod32(double op_dbl, u64_t op_64, u32_t p)
{
  u32_t h[3], res;

  if (op_dbl>=0.) {
    mpqs3_convert(h,op_dbl,op_64);
    return mpqs3_mod3(h,p);
  }
  mpqs3_convert(h,-op_dbl,-op_64);
  res=mpqs3_mod3(h,p);
  if (res) res=p-res;
  return res;
}
#endif

static void mpz_set_v(mpz_t gmp, u32_t *ulp)
{
  mpz_set_ui(gmp,0);
  if (ulp[2]&0x80000000) {
    mpz_set_ui(gmp,1);
    mpz_mul_2exp(gmp,gmp,32);
    mpz_neg(gmp,gmp);
  }
  mpz_add_ui(gmp,gmp,(ulong)(ulp[2]));
  mpz_mul_2exp(gmp,gmp,32);
  mpz_add_ui(gmp,gmp,(ulong)(ulp[1]));
  mpz_mul_2exp(gmp,gmp,32);
  mpz_add_ui(gmp,gmp,(ulong)(ulp[0]));
}



static u64_t mpqs3_inv_64(u64_t a)
{
  u64_t inv, h;

  inv=(u64_t)mpqs3_256_inv_table[(a&0xff)>>1];
  h=a*inv; h&=0xff00ULL;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffff0000ULL;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffffffff00000000ULL;
  h*=inv; inv-=h;
  if (inv*a!=1ULL) Schlendrian("mpqs3_inv_64");
  return inv;
}


static uchar mpqs3_jacobi_tab2[8]={0,0,0,1,0,1,0,0};

static int mpqs3_jacobi(ushort a, ushort b) /* always gcd(a,b)=1 and 0<a<b */
{
  int e, l, m, n, r;

  m=(int)a; n=(int)b;
  if (!(n&1)) Schlendrian("mpqs3_jacobi: b even\n");
  e=0;
  if (!(m&1)) {
    l=mpqs3_jacobi_tab2[n&7];
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
    l=mpqs3_jacobi_tab2[n&7];
    do {
      m/=2; e+=l;
    } while (!(m&1));
  }
  if (!m) Schlendrian("jacobi\n");
  e=e&1; e=1-2*e;
  return e;
}


static ushort mpqs3_invert(ushort a, ushort p)
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


static ushort mpqs3_powmod0(ushort a, ushort e, ushort p)
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

static ushort mpqs3_powmod(ushort a, ushort e, ushort p)
{
  ushort ex=e, rrr;
  u32_t aa, res, p32, inv, h, h1;

//rrr=mpqs3_powmod0(a,e,p);
  if (!ex) return 1;
  aa=(u32_t)a;
  p32=(u32_t)p;

  inv=(u32_t)mpqs3_256_inv_table[(p&0xff)>>1];
  h=p32*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=p32*inv; h&=0xffff0000;
  h*=inv; inv-=h;
  if (inv*p32!=1) Schlendrian("mpqs3_powmod_inv");
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
  complain("mpqs3_powmod %u %u %u  %u %u\n",(u32_t)a,(u32_t)e,p32,(u32_t)rrr,res);
*/
      return (ushort)res;
    }
    aa*=aa;
    h1=(aa&0x0000ffff)*inv; h1&=0x0000ffff; aa+=h1*p32; aa>>=16;
  }
  return 0; /* never reached */
}


static u32_t mpqs3_powmod_mm(u32_t a, u32_t e, u32_t p32, u32_t mmi, u32_t one)
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


static ushort mpqs3_sqrt_init(ushort p)
{
  ushort e, u, i, j;
  // SMJS Initialised g to stop compiler warning (probably if p was 0 which shouldn't happen)
  u32_t g = 0, b;

  if (p&2) return 0;
  if (p&4) return mpqs3_powmod(2,(p-1)>>2,p);
  e=0; u=p-1;
  while (!(u&1)) { e++; u>>=1; }
  for (i=2; i<p; i++) {
    g=mpqs3_powmod(i,u,p);
    if (g==1) continue;
    b=g;
    for (j=0; j<e-1; j++) { b*=b; b%=(u32_t)p; }
    if (b==p-1) break;
  }
  if (i>=p) Schlendrian("mpqs3_sqrt_init\n");
  return g;
}


static ushort mpqs3_sqrt0(ushort a, ushort p, ushort help)
                 /* a!=0, a<p, p prime, p!=2, (a/p)=1 */
{
  ushort e, u, i, l;
  u32_t b, g, r, k;

  if (p&2) return mpqs3_powmod(a,(p+1)>>2,p);
  if (p&4) {
    b=mpqs3_powmod(a,(p+3)>>3,p);
    g=b*b; g%=(u32_t)p;
    if (g==a) return b;
    b*=(u32_t)help; b%=(u32_t)p;
    return b;
  }
  e=0; u=p-1;
  while (!(u&1)) { e++; u>>=1; }
  g=help;
  r=mpqs3_powmod(a,u>>1,p);
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


static ushort mpqs3_sqrt(ushort a, ushort p, ushort help)
                 /* 0<a<p, p prime, p!=2, (a/p)=1 */
{
  u32_t e, u, i, l;
  u32_t b, g, r, k,aa;
  u32_t p32, mmi, h, h32, hh32, one, mminit;

  p32=(u32_t)p;
  mmi=(u32_t)mpqs3_256_inv_table[(p&0xff)>>1];
  h=p32*mmi; h&=0xff00;
  h*=mmi; mmi-=h;
  h=p32*mmi; h&=0xffff0000;
  h*=mmi; mmi-=h;
//  if (mmi*p32!=1) Schlendrian("mpqs3_sqrt_inv");
  mmi=-mmi;
  mminit=(0xffffffff)%p32; mminit++; /* mminit!=p32 */

  aa=(u32_t)a;
  h32=aa*mminit; MMREDUCE; aa=h32; if (aa>=p32) aa-=p32;
  h32=mminit; MMREDUCE; one=h32; if (one>=p32) one-=p32;

  if (p&2) {
    aa=mpqs3_powmod_mm(aa,(p32+1)>>2,p32,mmi,one);
    h32=aa;
    MMREDUCE;
    if (h32>=p32) h32-=p32;
    return (ushort)h32;
  }
  if (p&4) {
    b=mpqs3_powmod_mm(aa,(p32+3)>>3,p32,mmi,one);
    h32=b*b; MMREDUCE; g=h32; if (g>=p32) g-=p32;
    h32=b;
    if (g!=aa) h32*=help;
    MMREDUCE; b=h32; if (b>=p32) b-=p32;
    return (ushort)b;
  }
  e=0; u=p32-1;
  while (!(u&1)) { e++; u>>=1; }
  h32=((u32_t)help)*mminit; MMREDUCE; g=h32;
  r=mpqs3_powmod_mm(aa,u>>1,p32,mmi,one);
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


static u32_t mpqs3_inv(ushort a)
{
  u32_t inv, h;

  inv=(u32_t)mpqs3_256_inv_table[(a&0xff)>>1];
  h=(u32_t)a*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=(u32_t)a*inv; h&=0xffff0000;
  h*=inv; inv-=h;
  if (inv*(u32_t)a!=1) Schlendrian("mpqs3_inv");
  return inv;
}

/* ---------------------------------------------------------- */

static ushort mpqs3_multiplier_cand[4][MPQS3_MULT_NCAND]={
 { 1, 17, 33, 41, 57, 65, 73, 89 },
 { 3, 11, 19, 35, 43, 51, 59, 67 },
 { 5, 13, 21, 29, 37, 53, 61, 69 },
 { 7, 15, 23, 31, 39, 47, 55, 71 },
};

static double mpqs3_multiplier_cand_log[4][MPQS3_MULT_NCAND];

static char mpqs3_multiplier_cand_jacobi[4][MPQS3_MULT_NCAND][MPQS3_MULT_NTESTPR];


void mpqs3_init()
{
  u32_t i, j, k;
  u32_t p, add, d;

  mpz_init(mpqs3_N);
  mpz_init(mpqs3_kN);
  mpz_init(mpqs3_dummy);
  mpz_init(mpqs3_gmp_help);
  mpz_init(mpqs3_help);
  mpz_init(mpqs3_sq1);
  mpz_init(mpqs3_sq2);
  for (i=0; i<MPQS3_MAX_FACTORS; i++) mpz_init(mpqs3_factors[i]);
  mpqs3_sievearray=(uchar *)xmalloc((1+3*MPQS3_SIEVELEN)*sizeof(uchar));
  mpqs3_tinyarray=(uchar *)xmalloc(16*MPQS3_MAX_TINYPROD*sizeof(uchar));
  mpqs3_relations=(ushort **)xmalloc(MPQS3_MAX_NRELS*sizeof(ushort *));
  mpqs3_relations[0]=(ushort *)xmalloc(MPQS3_MAX_NRELS*MPQS3_REL_ENTRIES*sizeof(ushort));
  for (i=1; i<MPQS3_MAX_NRELS; i++)
    mpqs3_relations[i]=mpqs3_relations[i-1]+MPQS3_REL_ENTRIES;
  mpqs3_part_rels=(ushort **)xmalloc(MPQS3_MAX_NPRELS*sizeof(ushort *));
  mpqs3_part_rels[0]=(ushort *)xmalloc(MPQS3_MAX_NPRELS*MPQS3_REL_ENTRIES*sizeof(ushort));
  for (i=1; i<MPQS3_MAX_NPRELS; i++)
    mpqs3_part_rels[i]=mpqs3_part_rels[i-1]+MPQS3_REL_ENTRIES;
  mpqs3_comb_index=(ushort *)xmalloc(2*MPQS3_MAX_NCRELS*sizeof(ushort));
  mpqs3_rel_buffer=(ushort **)xmalloc(MPQS3_MAX_NRELS_BUF*sizeof(ushort *));
  mpqs3_rel_buffer[0]=(ushort *)xmalloc(MPQS3_MAX_NRELS_BUF*MPQS3_REL_ENTRIES*sizeof(ushort));
  for (i=1; i<MPQS3_MAX_NRELS_BUF; i++)
    mpqs3_rel_buffer[i]=mpqs3_rel_buffer[i-1]+MPQS3_REL_ENTRIES;

  mpqs3_gauss_row64=(u64_t *)xmalloc(MPQS3_GAUSS_MAX*MPQS3_GAUSS_MAX/64*sizeof(u64_t));

  mpqs3_prime_table[0]=2;
  mpqs3_prime_table[1]=3;
  mpqs3_prime_table[2]=5; mpqs3_prime_sqrthelp[2]=mpqs3_sqrt_init(5);
  i=3; p=7; add=4;
  while (i<MPQS3_MAX_NPRIMES) {
    for (j=2; j<i; j++) {
      d=(u32_t)mpqs3_prime_table[j];
      if (d*d>=p) break;
      if (p%d==0) break;
    }
    if (d*d>p) {
      mpqs3_prime_table[i]=(ushort)p;
      mpqs3_prime_sqrthelp[i]=mpqs3_sqrt_init(p);
      i++;
    }
    p+=add; add=6-add;
  }

  for (i=0; i<4; i++)
    for (j=0; j<MPQS3_MULT_NCAND; j++) {
      d=mpqs3_multiplier_cand[i][j];
      for (k=1; k<MPQS3_MULT_NTESTPR; k++) {
        p=(u32_t)mpqs3_prime_table[k];
        if (d%p) {
          mpqs3_multiplier_cand_jacobi[i][j][k]=
            (char)(mpqs3_jacobi((ushort)(d%p),(ushort)p));
        } else mpqs3_multiplier_cand_jacobi[i][j][k]=0;
      }
      mpqs3_multiplier_cand_log[i][j]=log((double)d)-1.5*log(2);
    }

  mpqs3_isinit=1;
}


static void mpqs3_choose_multiplier()
{
  ushort Nmod8, mult, p, mm, r;
  double value[MPQS3_MULT_NCAND], v, vp, vmin, dn;
  double plog[MPQS3_MULT_NTESTPR];
  char jac[MPQS3_MULT_NTESTPR], c;
  u32_t i, j;  

  for (i=1; i<MPQS3_MULT_NTESTPR; i++) {
    p=mpqs3_prime_table[i];
    r=(ushort)mpz_mod_ui(mpqs3_dummy,mpqs3_N,(ulong)p);
    jac[i]=(char)mpqs3_jacobi(r,p);
  }
  Nmod8=(ushort)mpz_mod_ui(mpqs3_dummy,mpqs3_N,8);
  dn=log(mpz_get_d(mpqs3_N))/log(2.);
  dn=164.9-dn;
  if (dn>7) mm=128; else mm=(ushort)(exp(dn*log(2)));
  for (i=1; i<MPQS3_MULT_NTESTPR; i++) {
    p=mpqs3_prime_table[i];
    plog[i]=log((double)p)/((double)p);
  }
  for (j=0; j<MPQS3_MULT_NCAND; j++) {
    mult=mpqs3_multiplier_cand[Nmod8/2][j];
    if (mult<=mm) {
      v=mpqs3_multiplier_cand_log[Nmod8/2][j];
      for (i=1; i<MPQS3_MULT_NTESTPR; i++) {
        c=mpqs3_multiplier_cand_jacobi[Nmod8/2][j][i];
        if (c) {
          vp=plog[i];
          if (c==jac[i]) v-=vp; else v+=vp;
        }
      }
    } else v=1000;
    value[j]=v;
  }
  vmin=value[0];
  for (j=1; j<MPQS3_MULT_NCAND; j++)
    if (value[j]<vmin) vmin=value[j];
  for (j=0; j<MPQS3_MULT_NCAND; j++)
    if (value[j]==vmin) break;
  if (j>=MPQS3_MULT_NCAND) complain("mpqs3_choose_multiplier.vmin\n");
  mpqs3_complexity=(vmin+log(mpz_get_d(mpqs3_N)))/log(2.);
  mpqs3_multiplier=mpqs3_multiplier_cand[Nmod8/2][j];
#ifdef MPQS3_STAT
// printf("%u ",mpqs3_multiplier);
#endif
  mpz_mul_ui(mpqs3_kN,mpqs3_N,(ulong)mpqs3_multiplier);
  mpqs3_kN_dbl=mpz_get_d(mpqs3_kN);
  mpz_fdiv_r_2exp(mpqs3_dummy,mpqs3_kN,64);
  mpqs3_kN_64=mpz_get_ull(mpqs3_dummy);
}


static void mpqs3_choose_parameter(int retry)
{
  int n, r;

  n=(int)(mpqs3_complexity+0.5);
  if (n>148) n=148;
  if (n<91) n=91;
  n=(n+3)/4;
  n-=23;  /* 0<=n<=14 */

  r=retry;
  if ((r) && (n<14)) { n++; r--; }
#ifdef MPQS3_STAT
stat_size[n]++;
#endif
  mpqs3_nFB=mpqs3_param[n][0];
  if (r) mpqs3_nFB+=r*mpqs3_nFB/4;
  if (mpqs3_nFB>=MPQS3_MAX_FBSIZE) mpqs3_nFB=MPQS3_MAX_FBSIZE-1;
  mpqs3_sievebegin=mpqs3_param[n][1];
  if (r) mpqs3_sievebegin+=r;
  mpqs3_nAdiv=mpqs3_param[n][2];
  mpqs3_nAdiv_total=mpqs3_param[n][3];
  if (mpqs3_nAdiv_total>MPQS3_MAX_ADIV_ALL)
    Schlendrian("mpqs3_choose_parameter: mpqs3_nAdiv_total too big\n");
  mpqs3_accept=mpqs3_param[n][4];
  mpqs3_td_begin=mpqs3_param[n][5];
  mpqs3_sievelen=mpqs3_param[n][6];
  if (mpqs3_sievelen>MPQS3_SIEVELEN) mpqs3_sievelen=MPQS3_SIEVELEN;

  mpqs3_disp=mpqs3_sievelen/2;
  mpqs3_nfactors=0; mpqs3_nrels=0; mpqs3_nlp=0; mpqs3_nprels=0; mpqs3_ncrels=0;
  for (n=0; n<MPQS3_NHASH; n++) mpqs3_hash_table[n][0]=0;
  for (n=0; n<MPQS3_NRELHASH; n++) mpqs3_rel_hash_table[n][0]=0;
}


static void mpqs3_generate_FB()
{
  ushort *fb, p, rest, i, nfb;

  fb=mpqs3_FB; nfb=0; i=0; mpqs3_nFBk=0;
  *fb++=2; *fb++=1; nfb++;
  for (i=1; i<MPQS3_MAX_NPRIMES; i++) {
    p=mpqs3_prime_table[i];
    if (p>MPQS3_FB_MAXPRIME) break;
    rest=(ushort)mpz_mod_ui(mpqs3_dummy,mpqs3_kN,(ulong)p);
    if (rest) {
      if (mpqs3_jacobi(rest,p)==1) {
        *fb++=p;
        *fb++=mpqs3_sqrt(rest,p,mpqs3_prime_sqrthelp[i]);
        nfb++;
        if (nfb>=mpqs3_nFB) { break; }
      }
    } else {
      if (mpqs3_multiplier%(u32_t)p)
        complain("mpqs3: N has small divisor: %u\n",p);
      mpqs3_FBk[mpqs3_nFBk++]=p;
    }
  }
  nfb-=((nfb-mpqs3_nAdiv_total)&3);
  mpqs3_nFB=nfb;
  mpqs3_nFBk_1=mpqs3_nFBk+1;
  mpqs3_nsp=1+mpqs3_nFB+mpqs3_nFBk;
  mpqs3_pmax=mpqs3_FB[2*(nfb-1)];
  mpqs3_FB[2*nfb]=mpqs3_sievelen;
  mpqs3_FB[2*nfb]=0xffff;
//printf("pmax: %u\n",mpqs3_pmax);
}

// SMJS Moved outside func because qsort comp func needs it
static double *alog;

static int mpqs3_SI_init()
{
  double d, prod_dbl;
  u64_t prod_64;
  u32_t a, i, j, k, l, n;
  ushort inv, p, *fb=mpqs3_FB, pk, tp_max;
  double A_div_log[MPQS3_MAX_ADIV_ALL]; /* maximal number of A's */
  double v, vmax;
  uchar lo;
  u16_t pb;

  d=mpz_get_d(mpqs3_kN); mpqs3_kN_dbl=d;
  d*=8.; d=sqrt(d); d/=(double)mpqs3_sievelen;
  d=log(d); d/=(double)mpqs3_nAdiv; d=exp(d);
  a=(u32_t)d;
//a=972;
/*  if (a>=mpqs3_pmax) { printf("%u %u ",a,mpqs3_pmax); return -1; }*/
  if (a>=mpqs3_pmax) a=mpqs3_pmax-1;
  for (i=0; i<mpqs3_nFB; i++) if (mpqs3_FB[2*i]>a) break;
  if (i>mpqs3_nAdiv_total/2) i-=mpqs3_nAdiv_total/2; else i=1;
  if (i<1) i=1; /* first prime >2 */
  if (mpqs3_FB[2*i]<3) return -2;
  if (i+mpqs3_nAdiv_total>mpqs3_nFB) {
    if (mpqs3_nFB<=mpqs3_nAdiv_total+1) return -3;
    i=mpqs3_nFB-mpqs3_nAdiv_total-1;
  }
  for (j=0; j<mpqs3_nAdiv_total; j++) {
    mpqs3_Adiv_all[j]=mpqs3_FB[2*(i+j)];
//printf("%u ",mpqs3_Adiv_all[j]);
    mpqs3_Adiv_all_sqrt[j]=mpqs3_FB[2*(i+j)+1];

    mpqs3_Adiv_all_sq_inv[j]=mpqs3_invert(mpqs3_Adiv_all_sqrt[j],mpqs3_Adiv_all[j]);
/* TODO optimise: */
{
  u32_t pp, p2, r0, rr, s, inv;

  pp=(u32_t)mpqs3_Adiv_all[j];
  p2=pp*pp;
  r0=mpz_mod_ui(mpqs3_dummy,mpqs3_kN,p2);
  s=(u32_t)mpqs3_Adiv_all_sqrt[j];
  inv=(u32_t)mpqs3_invert((ushort)(r0%pp),(ushort)pp);
  if (inv&1) inv+=pp; inv>>=1;
  rr=r0+(p2-(s*s)%p2);
  if (rr%pp) complain("sqrt2a\n");
  rr/=pp;
  rr*=s; rr%=pp;
  rr*=inv; rr%=pp;
  s+=rr*pp;
  mpqs3_Adiv_all_sq_mod_p2[j]=p2-s;
/* check 
  s=(u32_t)((((u64_t)s)*((u64_t)s))%((u64_t)p2));
  if (s!=r0) complain("sqrt2b\n");
*/
s=mpqs3_Adiv_all_sq_inv[j];
for (rr=0; rr<33; rr++) { s<<=1; if (s>=pp) s-=pp; }
mpqs3_Adiv_all_sq_inv[j]=s;
}
  }
//printf("\n");
  for (j=i+mpqs3_nAdiv_total; j<mpqs3_nFB; j++) {
    mpqs3_FB[2*(j-mpqs3_nAdiv_total)]=mpqs3_FB[2*j];
    mpqs3_FB[2*(j-mpqs3_nAdiv_total)+1]=mpqs3_FB[2*j+1];
  }
  mpqs3_nFB-=mpqs3_nAdiv_total;

/* put A divisors and k divisors at end of FB, used in next_pol */
  for (j=0; j<mpqs3_nAdiv_total; j++)
    mpqs3_FB[2*(mpqs3_nFB+j)]=mpqs3_Adiv_all[j];
  for (j=0; j<mpqs3_nFBk; j++)
    mpqs3_FB[2*(mpqs3_nFB+mpqs3_nAdiv_total+j)]=mpqs3_FBk[j];
  mpqs3_nFBx=mpqs3_nFB+mpqs3_nAdiv_total+mpqs3_nFBk;

/* tiny */
  mpqs3_sievebegin=1; mpqs3_sievebegink=0; mpqs3_tiny_prod=1; mpqs3_ntiny=0;
  tp_max=mpqs3_sievelen>>4;
  if (tp_max>MPQS3_MAX_TINYPROD) tp_max=MPQS3_MAX_TINYPROD;
  while (1) {
    if (mpqs3_sievebegin<mpqs3_nFB) p=mpqs3_FB[2*mpqs3_sievebegin];
    else p=tp_max+1;
    if (mpqs3_sievebegink<mpqs3_nFBk) pk=mpqs3_FBk[mpqs3_sievebegink];
    else pk=tp_max+1;
    if (p<pk) {
      if (mpqs3_tiny_prod*(u32_t)p>=tp_max) break;
      mpqs3_sievebegin++; mpqs3_ntiny++;
      mpqs3_tiny_prod*=(u32_t)p;
    } else {
      if (mpqs3_tiny_prod*(u32_t)pk>=tp_max) break;
      mpqs3_sievebegink++; mpqs3_ntiny++;
      mpqs3_tiny_prod*=(u32_t)pk;
    }
  }
#ifdef MPQS3_STAT
// for (j=0; j<mpqs3_sievebegin; j++) printf("%u ",mpqs3_FB[2*j]);
// printf("\n");
#endif
  if (mpqs3_nFB<mpqs3_td_begin) {
    if (mpqs3_nFB&3) complain("mpqs3_nFB\n");
//printf("%u\n",mpqs3_nFB);
    mpqs3_td_begin=mpqs3_nFB;
  }

#if 1
#ifdef HAVE_XMM_MUL
  p=mpqs3_FB[2];
  mpqs3_FB_inv[1]=mpqs3_inv(p);
  mpqs3_FB_inv_info[2]=p; mpqs3_FB_inv_info[3]=p;
  mpqs3_FB_inv_info[10]=(ushort)mpqs3_FB_inv[1];
  mpqs3_FB_inv_info[11]=(ushort)mpqs3_FB_inv[1];
  p=mpqs3_FB[4];
  mpqs3_FB_inv[2]=mpqs3_inv(p);
  mpqs3_FB_inv_info[4]=p; mpqs3_FB_inv_info[5]=p;
  mpqs3_FB_inv_info[12]=(ushort)mpqs3_FB_inv[2];
  mpqs3_FB_inv_info[13]=(ushort)mpqs3_FB_inv[2];
  p=mpqs3_FB[6];
  mpqs3_FB_inv[3]=mpqs3_inv(p);
  mpqs3_FB_inv_info[6]=p; mpqs3_FB_inv_info[7]=p;
  mpqs3_FB_inv_info[14]=(ushort)mpqs3_FB_inv[3];
  mpqs3_FB_inv_info[15]=(ushort)mpqs3_FB_inv[3];
  if (mpqs3_td_begin&3) complain("mpqs3_td_begin is not divisible by 4\n");
  for (j=4; j<mpqs3_td_begin; j+=4) {
    p=mpqs3_FB[2*j];
    mpqs3_FB_inv[j]=mpqs3_inv(p);
    mpqs3_FB_inv_info[4*j]=p; mpqs3_FB_inv_info[4*j+1]=p;
    mpqs3_FB_inv_info[4*j+8]=(ushort)mpqs3_FB_inv[j];
    mpqs3_FB_inv_info[4*j+9]=(ushort)mpqs3_FB_inv[j];
    p=mpqs3_FB[2*(j+1)];
    mpqs3_FB_inv[j+1]=mpqs3_inv(p);
    mpqs3_FB_inv_info[4*j+2]=p; mpqs3_FB_inv_info[4*j+3]=p;
    mpqs3_FB_inv_info[4*j+10]=(ushort)mpqs3_FB_inv[j+1];
    mpqs3_FB_inv_info[4*j+11]=(ushort)mpqs3_FB_inv[j+1];
    p=mpqs3_FB[2*(j+2)];
    mpqs3_FB_inv[j+2]=mpqs3_inv(p);
    mpqs3_FB_inv_info[4*j+4]=p; mpqs3_FB_inv_info[4*j+5]=p;
    mpqs3_FB_inv_info[4*j+12]=(ushort)mpqs3_FB_inv[j+2];
    mpqs3_FB_inv_info[4*j+13]=(ushort)mpqs3_FB_inv[j+2];
    p=mpqs3_FB[2*(j+3)];
    mpqs3_FB_inv[j+3]=mpqs3_inv(p);
    mpqs3_FB_inv_info[4*j+6]=p; mpqs3_FB_inv_info[4*j+7]=p;
    mpqs3_FB_inv_info[4*j+14]=(ushort)mpqs3_FB_inv[j+3];
    mpqs3_FB_inv_info[4*j+15]=(ushort)mpqs3_FB_inv[j+3];
  }
#else
  p=mpqs3_FB[2];
  mpqs3_FB_inv[1]=mpqs3_inv(p);
  mpqs3_FB_inv_info[2]=p; mpqs3_FB_inv_info[3]=p;
  mpqs3_FB_inv_info[6]=(ushort)mpqs3_FB_inv[1];
  mpqs3_FB_inv_info[7]=(ushort)mpqs3_FB_inv[1];
  if (mpqs3_td_begin&1) complain("mpqs3_td_begin is odd\n");
  for (j=2; j<mpqs3_td_begin; j+=2) {
    p=mpqs3_FB[2*j];
    mpqs3_FB_inv[j]=mpqs3_inv(p);
    mpqs3_FB_inv_info[4*j]=p; mpqs3_FB_inv_info[4*j+1]=p;
    mpqs3_FB_inv_info[4*j+4]=(ushort)mpqs3_FB_inv[j];
    mpqs3_FB_inv_info[4*j+5]=(ushort)mpqs3_FB_inv[j];
    p=mpqs3_FB[2*(j+1)];
    mpqs3_FB_inv[j+1]=mpqs3_inv(p);
    mpqs3_FB_inv_info[4*j+2]=p; mpqs3_FB_inv_info[4*j+3]=p;
    mpqs3_FB_inv_info[4*j+6]=(ushort)mpqs3_FB_inv[j+1];
    mpqs3_FB_inv_info[4*j+7]=(ushort)mpqs3_FB_inv[j+1];
  }
#endif
  for (; j<mpqs3_nFBx; j++) mpqs3_FB_inv[j]=mpqs3_inv(mpqs3_FB[2*j]);
  for (j=0; j<mpqs3_nFBk; j++) mpqs3_FBk_inv[j]=mpqs3_inv(mpqs3_FBk[j]);
  for (j=0; j<mpqs3_nAdiv_total; j++)
    mpqs3_FB_A_inv[j]=mpqs3_inv(mpqs3_Adiv_all[j]);
#endif

//zeitA(35);
/* for next_pol: */
{
  u32_t mmi, h32, hh32, invh, p32;
  ushort invhelp[MPQS3_MAX_ADIV_ALL], hhh;

  for (j=1; j<mpqs3_nFB; j++) {
    p=mpqs3_FB[2*j]; p32=(u32_t)p;
    mmi=-mpqs3_FB_inv[j];
    invhelp[0]=mpqs3_Adiv_all[0]%p;
    h32=(u32_t)(invhelp[0]);
    for (i=1; i<mpqs3_nAdiv_total; i++) {
      h32*=(u32_t)(mpqs3_Adiv_all[i]);
      MMREDUCE;
      invhelp[i]=(ushort)h32;
    }
    MMREDUCE; if (h32>=p) h32-=p; /* divide by 2^16 */
    invh=(u32_t)mpqs3_invert((ushort)h32,p);
    invh<<=1;
    
    for (i=mpqs3_nAdiv_total-1; i; i--) {
      h32=invh; h32*=(u32_t)(invhelp[i-1]);
      MMREDUCE; if (h32>=p) h32-=p;
      mpqs3_SI_inv_table[i][j]=(ushort)h32;
      h32=invh; h32*=(u32_t)mpqs3_Adiv_all[i];
      MMREDUCE; invh=h32;
    }
    if (invh>=p) invh-=p;
    mpqs3_SI_inv_table[0][j]=(ushort)invh;
  }
  for (j=0; j<mpqs3_nAdiv_total; j++) {
    u32_t h;
    ushort pa;

    p=mpqs3_Adiv_all[j];
    for (i=0; i<mpqs3_nAdiv_total; i++) {
      pa=mpqs3_Adiv_all[i]%p;
      if (pa) {
        h=(u32_t)mpqs3_invert(pa,p); //h*=2; if (h>=p) h-=p;
        h<<=17; h%=((u32_t)p);
      } else h=0;
      mpqs3_Adiv_SI_inv_table[i][j]=(ushort)h;
    }
    mpqs3_Adiv_disp[j]=(ushort)(mpqs3_disp%((u32_t)p));
  }
}
//zeitB(35);

  for (j=1; j<mpqs3_nFBx; j++) {
    mpqs3_FB_mm_inv[j]=-mpqs3_FB_inv[j];
    p=mpqs3_FB[2*j];
    mpqs3_FB_disp[j]=(ushort)(mpqs3_disp%((u32_t)p));
  }

#ifdef HAVE_XMM_MUL
  for (j=0; j<mpqs3_nFBx; j+=8) {
    for (i=0; i<8; i++) mpqs3_FB_np_px[2*j+i]=mpqs3_FB[2*j+2*i];
    for (i=0; i<8; i++) mpqs3_FB_np_px[2*j+8+i]=mpqs3_FB[2*j+1+2*i];
  }
#else
  for (j=0; j<mpqs3_nFBx; j+=4) {
    mpqs3_FB_np_p[2*j]=mpqs3_FB[2*j];
    mpqs3_FB_np_p[2*j+1]=mpqs3_FB[2*j+2];
    mpqs3_FB_np_p[2*j+2]=mpqs3_FB[2*j+4];
    mpqs3_FB_np_p[2*j+3]=mpqs3_FB[2*j+6];
    mpqs3_FB_np_p[2*j+4]=mpqs3_FB[2*j+1];
    mpqs3_FB_np_p[2*j+5]=mpqs3_FB[2*j+3];
    mpqs3_FB_np_p[2*j+6]=mpqs3_FB[2*j+5];
    mpqs3_FB_np_p[2*j+7]=mpqs3_FB[2*j+7];
  }
#endif

/* compute log-approximations */
  d=mpz_get_d(mpqs3_kN);
  d/=8.; d=sqrt(d); d*=(double)mpqs3_sievelen;
  mpqs3_maxval_dbl=d;
  d=log(d);
  d-=log((double)(1ULL<<mpqs3_accept));
  mpqs3_logbound=128;  /* ??? */
  d=(double)(mpqs3_logbound)/d;
#if 0
  for (i=0; i<mpqs3_nFB; i++) {
    p=mpqs3_FB[2*i];
    mpqs3_FB_log[i]=(uchar)(0.5+d*log((double)p));
  }
#else
  if (mpqs3_nFB<8) return -7;
{
  uchar curlog;
  u32_t i0, i1, i2;
  ushort pt;
  double dp;

  for (i=0; i<8; i++) {
    p=mpqs3_FB[2*i];
    mpqs3_FB_log[i]=(uchar)(0.5+d*log((double)p));
  }
  curlog=mpqs3_FB_log[i-1]-1;  
  while (i<mpqs3_nFB) {
    if (curlog>128) return -8;
    curlog++;
    dp=exp((0.5+(double)curlog)/d);
    if (dp>65535.) pt=65535; else pt=(ushort)dp;
    i0=i-1; i1=mpqs3_nFB-1;
    if (mpqs3_FB[2*i1]>pt) {
      while (i1>i0+1) {
        i2=((i0+i1)>>1);
        if (mpqs3_FB[2*i2]>pt) i1=i2; else i0=i2;
      }
      i1--;
    }
    for (; i<=i1; i++) mpqs3_FB_log[i]=curlog;
  }
/*
  for (i=0; i<mpqs3_nFB; i++) {
    p=mpqs3_FB[2*i];
    if (mpqs3_FB_log[i]!=(uchar)(0.5+d*log((double)p))) complain("diff %u %u   %u %u\n",i,p,mpqs3_FB_log[i],(u32_t)(0.5+d*log((double)p)));
  }
*/
}
#endif
  for (i=0; i<mpqs3_nFBk; i++) {
    p=mpqs3_FBk[i];
    mpqs3_FBk_log[i]=(uchar)(0.5+d*log((double)p));
  }
  for (i=0; i<mpqs3_nAdiv_total; i++) {
    p=mpqs3_Adiv_all[i];
    mpqs3_Adiv_log[i]=(uchar)(0.5+d*log((double)p));
  }
#ifdef ASM_MPQS3_SIEVE0
/* fill mpqs3_FB0 */
  lo=0; n=0; i=mpqs3_sievebegin;
#if 1
  pb=mpqs3_sievelen>>3;
  for (; i<mpqs3_nFB; i++) {
    p=mpqs3_FB[2*i]; if (p>pb) break;
    if (lo<mpqs3_FB_log[i]) {
      mpqs3_FB0[n++]=0;
      lo=mpqs3_FB_log[i];
      mpqs3_FB0[n++]=(u16_t)lo;
    }
    mpqs3_FB0[n++]=p;
  }
  mpqs3_FB0[n++]=0; mpqs3_FB0[n++]=0;
#endif
  for (j=4; j; j--) {
    pb=mpqs3_sievelen/j;
    for (; i<mpqs3_nFB; i++) {
      p=mpqs3_FB[2*i]; if (p>pb) break;
      if (lo<mpqs3_FB_log[i]) {
        mpqs3_FB0[n++]=0;
        lo=mpqs3_FB_log[i];
        mpqs3_FB0[n++]=(u16_t)lo;      
      }
      mpqs3_FB0[n++]=p;
    }
    mpqs3_FB0[n++]=0; mpqs3_FB0[n++]=0;
  }
  for (; i<mpqs3_nFB; i++) {
    p=mpqs3_FB[2*i];
    if (lo<mpqs3_FB_log[i]) {
      mpqs3_FB0[n++]=0;
      lo=mpqs3_FB_log[i];
      mpqs3_FB0[n++]=(u16_t)lo;
    }
    mpqs3_FB0[n++]=p;
  }
  mpqs3_FB0[n++]=0; mpqs3_FB0[n++]=0;
#endif

{
/* experimental */
// SMJS Make static outside func because qsort comp func needs it  double *alog;
  ushort tmp;
  u32_t *asort;

{

/* compute mask for choice of A-divisors */
  j=0;
  if (mpqs3_nAdiv_total>8) {
    k=0; if (mpqs3_nAdiv>8) k=mpqs3_nAdiv-8;
    for (; k<=mpqs3_nAdiv; k++) {
      u32_t i0, i1;
      ushort m16;

      for (i0=mpqs3_A_table_n[k]; i0<mpqs3_A_table_n[k+1]; i0++) {
        m16=(ushort)(mpqs3_A_table[i0]); m16<<=8;
        if (m16>>mpqs3_nAdiv_total) break;
        for (i1=mpqs3_A_table_n[mpqs3_nAdiv-k];
             i1<mpqs3_A_table_n[mpqs3_nAdiv-k+1]; i1++) {
          mpqs3_A_mask[j]=m16+(ushort)(mpqs3_A_table[i1]);
          j++;
        }
      }
    }
  } else {
    for (i=mpqs3_A_table_n[mpqs3_nAdiv];
         i<mpqs3_A_table_n[mpqs3_nAdiv+1]; i++) {
      ushort m16;

      m16=(ushort)(mpqs3_A_table[i]);
      if (m16>>mpqs3_nAdiv_total) break;
      mpqs3_A_mask[j]=m16;
      j++;
    }
  }
/* sort mpqs3_A_mask */
  alog=(double *)xmalloc(j*sizeof(double));

/* ensure that |C|<2^96 */
  for (i=0; i<mpqs3_nAdiv_total; i++)
    A_div_log[i]=log((double)(mpqs3_Adiv_all[i]))/log(2.);
  d=log(mpqs3_kN_dbl)/log(2.);
  for (i=0; i<j; i++) {
    v=0;
    for (k=0; k<mpqs3_nAdiv_total; k++)
      if (mpqs3_A_mask[i]&(1<<k)) v+=A_div_log[k];
    alog[i]=v;
    if (d-v>95.95) {
      j--; /* printf(":");*/
      mpqs3_A_mask[i]=mpqs3_A_mask[j];
      i--;
    }
  }

  d=mpz_get_d(mpqs3_kN); mpqs3_kN_dbl=d;
  d*=8.; d=sqrt(d); d/=(double)mpqs3_sievelen;
  d=log(d)/log(2.);
  for (i=0; i<j; i++)
    if (d>alog[i]) alog[i]=d-alog[i];
    else alog[i]-=d;

//zeitA(43);
#if 1
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
/*
  for (i=0; i+1<j; i++)
    if (alog[a1[i]]>alog[a1[i+1]])
      complain("sort %u  %g %g\n",i,alog[a1[i]],alog[a1[i+1]]);
*/
  for (i=0; i<j; i++) a1[i]=(u32_t)(mpqs3_A_mask[a1[i]]);
  for (i=0; i<j; i++) mpqs3_A_mask[i]=(ushort)(a1[i]);
  free(asort);
}
#else
  asort=(u32_t *)xmalloc(j*sizeof(u32_t));
  for (i=0; i<j; i++) asort[i]=i;
  qsort(asort,(size_t)j,sizeof(*asort),A_index_cmp);
  for (i=0; i<j; i++) asort[i]=(u32_t)(mpqs3_A_mask[asort[i]]);
  for (i=0; i<j; i++) mpqs3_A_mask[i]=(ushort)(asort[i]);
  free(asort);
#endif
//zeitB(43);
  free(alog);
//printf("a-dev: %g %g %g\n",alog[0],alog[j/2],alog[j-1]);
}
}

  mpqs3_nA=j; mpqs3_A_index=-1; mpqs3_B_index=(1<<mpqs3_nAdiv)-1;

//  mpqs3_FB[2*mpqs3_nFB]=mpqs3_sievelen;
//  mpqs3_FB[2*mpqs3_nFB]=0xffff;
  return 1;
}

// SMJS Moved from inside above function
int A_index_cmp(const void *a, const void *b)
{
  if (alog[*((u32_t *)(a))]>=alog[*((u32_t *)(b))]) return 1;
  return -1;
}

static int mpqs3_next_pol()
{
  u32_t i, ind, j, inv;
//  u32_t add[3], bi[3], prod_other[3], inv;
  ushort p, s1, s2, *fb=mpqs3_FB, *fbs=mpqs3_FB_start;
  short sh;
  u32_t bb, cc, c1, c2;
  ushort mask;
  double d, A2n;
  ushort aa;
  u32_t al, *ul;
  u32_t mod, B_modp;
  double add_dbl, bi_dbl;
  u64_t add_64, bi_64;

  mpqs3_B_index++;
  if (mpqs3_B_index>=1<<(mpqs3_nAdiv-1)) {
    u32_t Bi_mul[MPQS3_MAX_ADIV_ALL], Adiv_nr[MPQS3_MAX_ADIV_ALL];

#ifdef MPQS3_ZEIT
zeitA(7);
#endif
//zeitA(30);
//zeitA(26);
//zeitA(60);
    mpqs3_A_index++;
    if (mpqs3_A_index>=mpqs3_nA) return 0;
    mask=mpqs3_A_mask[mpqs3_A_index];
    for (i=0; i<mpqs3_nAdiv_total; i++)
      if (mask & (1<<i)) mpqs3_Adiv_active[i]=1;
      else mpqs3_Adiv_active[i]=0;
    j=0;
    for (i=0; i<mpqs3_nAdiv_total; i++)
      if (mpqs3_Adiv_active[i]) {
        mpqs3_Adiv[j]=mpqs3_Adiv_all[i];
        mpqs3_Adiv_sqrt[j]=mpqs3_Adiv_all_sqrt[i];
        Adiv_nr[j]=i;
        j++;
      }
    mpqs3_A_64=1;
    mpqs3_A_dbl=1.;
    for (i=0; i<mpqs3_nAdiv; i++) {
      mpqs3_A_dbl*=(double)(mpqs3_Adiv[i]);
      mpqs3_A_64*=(u64_t)(mpqs3_Adiv[i]);
    }
/* precomputations for sieve init */
    mpqs3_sinit_tab[0]=0; mpqs3_sinit_tab[1]=0;
    i=1; A2n=mpqs3_A_dbl*0.5;
    for (j=1; j<MPQS3_SI_NVALUES; j++) {
      d=sqrt(mpqs3_kN_dbl+mpqs3_maxval_dbl*A2n)/mpqs3_A_dbl;
      A2n*=0.5;
      if (d<(double)mpqs3_disp) s1=mpqs3_disp-(ushort)d; else s1=0;
      s1-=(s1&63);
      if (s1>mpqs3_sinit_tab[2*i-2]) {
        mpqs3_sinit_tab[2*i]=s1; mpqs3_sinit_tab[2*i+1]=j;
        i++;
      }
    }
    A2n*=2.;
    for (j=1; j<MPQS3_SI_NVALUES; j++) {
      d=mpqs3_kN_dbl-mpqs3_maxval_dbl*A2n; A2n*=2.;
      if (d<0.) d=0.;
      d=sqrt(d)/mpqs3_A_dbl; if (d<1.) d=1.;
      if (d<(double)mpqs3_disp) s1=mpqs3_disp-(ushort)d; else s1=0;
      s1-=(s1&63);
      if (s1>mpqs3_sinit_tab[2*i-2]) {
        mpqs3_sinit_tab[2*i]=s1; mpqs3_sinit_tab[2*i+1]=MPQS3_SI_NVALUES-1-j;
        i++;
      }
    }
    if (A2n!=mpqs3_A_dbl) complain("precomp sieve init %lf %lf\n",A2n,mpqs3_A_dbl);

    for (j=i; j<2*i; j++) {
      mpqs3_sinit_tab[2*j]=mpqs3_sievelen-mpqs3_sinit_tab[2*(2*i-1-j)];
      if (j<2*i-1) mpqs3_sinit_tab[2*j+1]=mpqs3_sinit_tab[2*(2*i-2-j)+1];
    }
/* compute B_i's */ /* TODO optimise */
    for (i=0; i<mpqs3_nAdiv; i++) {
      p=mpqs3_Adiv[i];
      bi_64=mpqs3_A_64*mpqs3_inv_64((u64_t)p);
      bi_dbl=mpqs3_A_dbl/(double)p;
      inv=mpqs3_dbl64mod16(bi_dbl,bi_64,p);
      inv=mpqs3_invert((ushort)inv,p);
      mod=(inv*((u32_t)mpqs3_Adiv_sqrt[i]))%((u32_t)p);
      Bi_mul[i]=mod;
      mpqs3_Bi_dbl[i]=bi_dbl*(double)mod;
      mpqs3_Bi_64[i]=bi_64*(u64_t)mod;
    }
    mpqs3_B_dbl=0.; mpqs3_B_64=0;
    for (i=0; i<mpqs3_nAdiv; i++) {
      mpqs3_B_dbl+=mpqs3_Bi_dbl[i];
      mpqs3_B_64+=mpqs3_Bi_64[i];
    }
    mpqs3_2B_dbl=2.*mpqs3_B_dbl;
    mpqs3_2B_64=2*mpqs3_B_64;
//zeitA(63);
//zeitB(26);
#ifdef MPQS3_ZEIT
zeitA(13);
#endif
#ifdef ASM_MPQS3_NEXT_POL
//zeitA(36);
{
    u32_t jj;
    u64_t pi0_64, *fbs64;

    pi0_64=1ULL<<(16-mpqs3_nAdiv); pi0_64+=(pi0_64<<16); pi0_64+=(pi0_64<<32);
    fbs64=(u64_t *)fbs;
#ifdef HAVE_XMM_MUL
    for (i=0; i<mpqs3_nFB; i+=8) {
      *fbs64++=pi0_64;
      *fbs64++=pi0_64;
      *fbs64++=0;
      *fbs64++=0;
    }
    for (j=0; j<mpqs3_nAdiv; j++) {
      jj=Adiv_nr[j];
      asm3_next_pol10_xmm((mpqs3_nFB+7)/8,&(mpqs3_SI_inv_table[jj][0]),
                          &(mpqs3_SI_add[j][0]),Bi_mul[j]);
    }
    asm3_next_pol11_xmm((mpqs3_nFB+7)/8);
#else
    for (i=0; i<mpqs3_nFB; i+=4) {
      *fbs64++=pi0_64;
      *fbs64++=0;
    }
    for (j=0; j<mpqs3_nAdiv; j++) {
      jj=Adiv_nr[j];
      asm3_next_pol10((mpqs3_nFB+3)/4,&(mpqs3_SI_inv_table[jj][0]),
                      &(mpqs3_SI_add[j][0]),Bi_mul[j]);
    }
    asm3_next_pol11((mpqs3_nFB+3)/4);
#endif
}
//zeitB(36);
#else

/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/
{
    u32_t mmi, h32, hh32, p32;
    u32_t cc1, cc2, bbb;
    u32_t jj, pi, pi0;

    pi0=1<<(16-mpqs3_nAdiv);
    for (i=1; i<mpqs3_nFB; i++) {
      p=fb[2*i]; p32=(u32_t)p;
      mmi=mpqs3_FB_mm_inv[i];

      bbb=0; pi=pi0;
      for (j=0; j<mpqs3_nAdiv; j++) {
        jj=Adiv_nr[j];
        cc=mpqs3_SI_inv_table[jj][i];
        cc*=Bi_mul[j]; h32=cc;
        MMREDUCE; cc=h32; if (cc>=p) cc-=p;
        mpqs3_SI_add[j][i]=(ushort)cc;
        bbb+=mpqs3_SI_add[j][i]; if (bbb>=p) bbb-=p;
        pi*=((u32_t)mpqs3_SI_inv_table[jj][i]); h32=pi;
        MMREDUCE; pi=h32;
      }
      cc=bbb;
      if (cc&1) cc+=p; cc>>=1;
      cc=mpqs3_FB_disp[i]+(p-cc); if (cc>=p) cc-=p;
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
#ifdef MPQS3_ZEIT
zeitB(13);
#endif
//zeitB(63);
//zeitA(64);
    mpqs3_C_dbl=(mpqs3_B_dbl*mpqs3_B_dbl-mpqs3_kN_dbl)/mpqs3_A_dbl;
    mpqs3_A_inv_64=mpqs3_inv_64(mpqs3_A_64);
    mpqs3_C_64=mpqs3_B_64*mpqs3_B_64-mpqs3_kN_64;
    mpqs3_C_64*=mpqs3_A_inv_64;

/* divisors of A */
    for (i=0; i<mpqs3_nAdiv_total; i++)
      if (!mpqs3_Adiv_active[i]) {
        u32_t mmi, h32, hh32, pi, jj, p32;

        p=mpqs3_Adiv_all[i]; p32=(u32_t)p;
        mmi=-mpqs3_inv(p);

        bb=0; pi=1<<(16-mpqs3_nAdiv);
        for (j=0; j<mpqs3_nAdiv; j++) {
          jj=Adiv_nr[j];
          cc=mpqs3_Adiv_SI_inv_table[jj][i];
          cc*=Bi_mul[j]; h32=cc;
          MMREDUCE; cc=h32; if (cc>=p) cc-=p;
          mpqs3_Adiv_SI_add[j][i]=(ushort)cc;
mpqs3_SI_add[j][mpqs3_nFB+i]=(ushort)cc;
          bb+=mpqs3_Adiv_SI_add[j][i]; if (bb>=p) bb-=p;
          pi*=((u32_t)mpqs3_Adiv_SI_inv_table[jj][i]); h32=pi;
          MMREDUCE; pi=h32;
        }
        cc=bb;
        if (cc&1) cc+=p; cc>>=1;
        cc=mpqs3_Adiv_disp[i]+(p-cc); if (cc>=p) cc-=p;
        c1=mpqs3_Adiv_all_sqrt[i];
        h32=c1*pi; MMREDUCE; c1=h32;
        if (c1>=p) c1-=p;
        c2=p-c1;
        c1+=cc; if (c1>=p) c1-=p;
        c2+=cc; if (c2>=p) c2-=p;
        mpqs3_Adiv_start1[i]=(ushort)c1; mpqs3_Adiv_start2[i]=(ushort)c2;
fbs[2*(mpqs3_nFB+i)]=(ushort)c1; fbs[2*(mpqs3_nFB+i)+1]=(ushort)c2;
      } else {
        u32_t mmi, h32, hh32, pi, jj, cc0, j0, p32, cc1;

        p=mpqs3_Adiv_all[i]; p32=(u32_t)p;
        cc0=0;
        mmi=-mpqs3_inv(p);
        for (j=0; j<mpqs3_nAdiv; j++) {
          jj=Adiv_nr[j];
          if (jj!=i) {
            cc=mpqs3_Adiv_SI_inv_table[jj][i];
            cc*=Bi_mul[j]; h32=cc;
            MMREDUCE; cc=h32; if (cc>=p) cc-=p;
            mpqs3_Adiv_SI_add[j][i]=(ushort)cc;
mpqs3_SI_add[j][mpqs3_nFB+i]=(ushort)cc;
          } else {
            u32_t k;

#if 1
            cc=mpqs3_dbl64mod32(mpqs3_Bi_dbl[j],mpqs3_Bi_64[j],p32*p32)+mpqs3_Adiv_all_sq_mod_p2[jj];
#else
/* include modmul32 etc */
            cc=Bi_mul[j]; modulo32=p32*p32;
            for (k=0; k<mpqs3_nAdiv; k++) {
              if (k==j) continue;
              cc=modmul32(cc,mpqs3_Adiv_all[Adiv_nr[k]]);
            }
            cc+=mpqs3_Adiv_all_sq_mod_p2[jj];
#endif
            if (cc%p32) complain("not div by p\n");
            cc/=p32;
            cc*=Bi_mul[j]; h32=cc;
            MMREDUCE; cc=h32;
            cc*=mpqs3_Adiv_all_sq_inv[jj]; h32=cc;
            MMREDUCE; cc=h32; if (cc>=p32) cc-=p32;
            mpqs3_Adiv_SI_add[j][i]=(ushort)cc;
mpqs3_SI_add[j][mpqs3_nFB+i]=(ushort)cc;
          }
          cc0+=cc; if (cc0>=p32) cc0-=p32;
        }
        if (cc0) cc0=p32-cc0;
/* using relation: 2*start_value=-sum of change_values */
        if (cc0&1) cc0+=p32; cc0>>=1;
        cc0+=mpqs3_Adiv_disp[i]; if (cc0>=p) cc0-=p;
        mpqs3_Adiv_start1[i]=(ushort)cc0;
//        if (mpqs3_Adiv_start1[i]!=(ushort)cc0) complain("neu %u %u  %u\n",cc0,mpqs3_Adiv_start1[i],p32);
fbs[2*(mpqs3_nFB+i)]=(ushort)cc0;
      }
/* divisors of k */
    for (i=0; i<mpqs3_nFBk; i++) { /* slow, but seems to be fast enough */
      p=mpqs3_FBk[i]; B_modp=0;
      bb=(u32_t)mpqs3_dbl64mod16(mpqs3_A_dbl,mpqs3_A_64,p);
      bb=mpqs3_invert(bb,p);
      for (j=0; j<mpqs3_nAdiv; j++) {
        cc=(u32_t)mpqs3_dbl64mod16(mpqs3_Bi_dbl[j],mpqs3_Bi_64[j],p);
        B_modp+=((u32_t)p-cc);
        cc=2*bb*cc; cc%=(u32_t)p;
        mpqs3_SIk_add[j][i]=(ushort)cc;
mpqs3_SI_add[j][mpqs3_nFB+mpqs3_nAdiv_total+i]=(ushort)cc;
      }
      c1=B_modp*bb+mpqs3_disp; c1%=(u32_t)p;
      mpqs3_FBk_start[i]=(ushort)c1;
fbs[2*(mpqs3_nFB+mpqs3_nAdiv_total+i)]=(ushort)c1;
    }
    mpqs3_B_index=0;
//zeitB(64);
//zeitB(60);
//zeitB(30);
#ifdef MPQS3_ZEIT
zeitB(7);
#endif
    return 1;
  }
//zeitA(61);
  ind=mpqs3_B_index; i=0;
  while (1) {
    if (ind&1) break;
    i++; ind>>=1;
  }
  ind>>=1; ind++;
  add_dbl=2.*mpqs3_Bi_dbl[i];
  add_64=2*mpqs3_Bi_64[i];

//zeitA(62);
  if (ind&1) {
    mpqs3_B_dbl=mpqs3_B_dbl-add_dbl;
    mpqs3_B_64=mpqs3_B_64-add_64;
#ifdef ASM_MPQS3_NEXT_POL
#ifdef HAVE_XMM_MUL
    asm3_next_pol3plus_xmm((mpqs3_nFBx+7)/8,&(mpqs3_SI_add[i][0]));
#else
    asm3_next_pol3plus((mpqs3_nFBx+3)/4,&(mpqs3_SI_add[i][0]));
#endif
#else
    for (j=1; j<mpqs3_nFBx; j++) {
      p=fb[2*j]; /*s1=fbs[2*j]; s2=fbs[2*j+1];*/
      aa=mpqs3_SI_add[i][j];
      fbs[2*j]+=aa; fbs[2*j+1]+=aa;
      if (fbs[2*j]>=p) fbs[2*j]-=p;
      if (fbs[2*j+1]>=p) fbs[2*j+1]-=p;
    }
#endif
  } else {
    mpqs3_B_dbl=mpqs3_B_dbl+add_dbl;
    mpqs3_B_64=mpqs3_B_64+add_64;
#ifdef ASM_MPQS3_NEXT_POL
#ifdef HAVE_XMM_MUL
    asm3_next_pol3minus_xmm((mpqs3_nFBx+7)/8,&(mpqs3_SI_add[i][0]));
#else
    asm3_next_pol3minus((mpqs3_nFBx+3)/4,&(mpqs3_SI_add[i][0]));
#endif
#else
    for (j=1; j<mpqs3_nFBx; j++) {
      p=fb[2*j]; /*s1=fbs[2*j]; s2=fbs[2*j+1];*/
      aa=p-mpqs3_SI_add[i][j];
      fbs[2*j]+=aa; fbs[2*j+1]+=aa;
      if (fbs[2*j]>=p) fbs[2*j]-=p;
      if (fbs[2*j+1]>=p) fbs[2*j+1]-=p;
    }
#endif
  }
//zeitB(62);
  for (j=0; j<mpqs3_nFBk; j++)
    mpqs3_FBk_start[j]=fbs[2*(mpqs3_nFB+mpqs3_nAdiv_total+j)];

  for (j=0; j<mpqs3_nAdiv_total; j++) {
    mpqs3_Adiv_start1[j]=fbs[2*(mpqs3_nFB+j)];
    mpqs3_Adiv_start2[j]=fbs[2*(mpqs3_nFB+j)+1];
  }

  mpqs3_2B_dbl=2.*mpqs3_B_dbl;
  mpqs3_2B_64=2*mpqs3_B_64;

  mpqs3_C_dbl=(mpqs3_B_dbl*mpqs3_B_dbl-mpqs3_kN_dbl)/mpqs3_A_dbl;
  mpqs3_C_64=mpqs3_B_64*mpqs3_B_64-mpqs3_kN_64;
  mpqs3_C_64*=mpqs3_A_inv_64;
//zeitB(61);
  return 1;
}


static void mpqs3_sieve()
{
  uchar *sv=mpqs3_sievearray, *fblog=mpqs3_FB_log, *svend;
  ushort *fb=mpqs3_FB, *fbs=mpqs3_FB_start, *fb0=mpqs3_FB0;
  ushort p, s1, s2, end, pb, save;
  u32_t i;
  uchar lo;
  u32_t *ulsv, *ulsvend, m, m1, m2, *ulsv0;
  u64_t *ullsv, *ullsvend, mask0, mask1, m64[2];

/* initialize odd tiny primes */
  ulsv=(u32_t *)mpqs3_tinyarray;
  if (mpqs3_ntiny==0) {
    *ulsv=0;
  } else {
    u32_t j, p4;
    u32_t h, *ulsv0;

    for (i=0; i<mpqs3_tiny_prod; i++) *ulsv++=0;
    for (i=1; i<mpqs3_sievebegin; i++) {
      p=fb[2*i]; lo=fblog[i]; p4=((u32_t)p<<2);
      sv=mpqs3_tinyarray;
      s1=fbs[2*i]; s2=fbs[2*i+1];
      for (j=0; j<mpqs3_tiny_prod; j+=p4) {
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
      }
    }
    for (i=0; i<mpqs3_sievebegink; i++) {
      p=mpqs3_FBk[i]; lo=mpqs3_FBk_log[i]; p4=((u32_t)p<<2);
      sv=mpqs3_tinyarray;
      s1=mpqs3_FBk_start[i];
      for (j=0; j<mpqs3_tiny_prod; j+=p4) {
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
      }
    }
  }

/* add contributions of powers of 2 */
  lo=mpqs3_FB_log[0];
  m=((u32_t)mpqs3_A_64)*(-1-mpqs3_disp)+(u32_t)mpqs3_2B_64;
  m*=(-1-mpqs3_disp);
  m+=(u32_t)mpqs3_C_64;
  m1=((u32_t)mpqs3_A_64)*(-3-2*mpqs3_disp)+(u32_t)mpqs3_2B_64;
  m2=2*((u32_t)mpqs3_A_64);

  m64[0]=0; m64[1]=0;
  sv=(uchar *)(&(m64[0]));
#if 1
  if (m&1) {
    i=0;
    m-=m1; m1-=m2;
  } else i=1;
/* update (m1,m2) for steps of 2 */
  m1+=m1; m1-=m2; m2<<=2;
  for (; i<16; i+=2) {
    m1+=m2; m+=m1;
//    if (m&6) complain("si mod 8\n");
    if (m&8) { sv[i]+=3*lo; continue; }
    if (m&16) { sv[i]+=4*lo; continue; }
    sv[i]+=6*lo;
  }
#else
  for (i=0; i<16; i++) {
    m1+=m2; m+=m1;
    if (m&1) continue;
    if (m&6) complain("si mod 8\n");
    if (m&8) { sv[i]+=3*lo; continue; }
    if (m&16) { sv[i]+=4*lo; continue; }
    sv[i]+=6*lo;
  }
#endif

//zeitA(70);
#ifndef ASM_MPQS3_SIEVE_INIT0
  memcpy(mpqs3_tinyarray+mpqs3_tiny_prod,mpqs3_tinyarray,mpqs3_tiny_prod);
  memcpy(mpqs3_tinyarray+2*mpqs3_tiny_prod,mpqs3_tinyarray,2*mpqs3_tiny_prod);
  memcpy(mpqs3_tinyarray+4*mpqs3_tiny_prod,mpqs3_tinyarray,4*mpqs3_tiny_prod);
  memcpy(mpqs3_tinyarray+8*mpqs3_tiny_prod,mpqs3_tinyarray,8*mpqs3_tiny_prod);
  memcpy(mpqs3_tinyarray+16*mpqs3_tiny_prod,mpqs3_tinyarray,16);
  mask0=m64[0]; mask1=m64[1];
  ullsv=(u64_t *)mpqs3_tinyarray;
  ullsvend=ullsv+2*mpqs3_tiny_prod+2;
  while (ullsv<ullsvend) {
    *ullsv+++=mask0;
    *ullsv+++=mask1;
  }
#else
  asm_sieve_init0(mpqs3_tinyarray,&(m64[0]),mpqs3_tiny_prod);
#endif
//zeitB(70);

//zeitA(71);
/* contribution of log |q(x)| */
#ifdef ASM_MPQS3_SIEVE_INIT
  mask0=0x0101010101010101ULL*(u64_t)(lo);
  asm_sieve_init(mpqs3_sievearray,mpqs3_sievelen,mpqs3_sinit_tab,&mask0,mpqs3_tinyarray,mpqs3_tiny_prod);
#else
{
  u64_t *tptr0, *tptr1, *tptr, mask;

  tptr0=(u64_t *)mpqs3_tinyarray;
  tptr1=tptr0+2*mpqs3_tiny_prod;
  tptr=tptr0;
  i=0;
  ullsv=(u64_t *)mpqs3_sievearray;
  while (1) {
    if (mpqs3_sinit_tab[2*i]==mpqs3_sievelen) break;
    ullsvend=((u64_t *)mpqs3_sievearray)+mpqs3_sinit_tab[2*i+2]/8;
//if ((mpqs3_sinit_tab[2*i+2]-mpqs3_sinit_tab[2*i])&31) complain("32\n");
    mask=0x0101010101010101ULL*(u64_t)(lo*mpqs3_sinit_tab[2*i+1]);

    while (ullsv<ullsvend) {
      *ullsv++=*tptr+++mask;
      *ullsv++=*tptr+++mask;
      *ullsv++=*tptr+++mask;
      *ullsv++=*tptr+++mask;
      if (tptr>=tptr1) tptr-=(2*mpqs3_tiny_prod);
    }
    i++;
  }
}
#endif
//zeitB(71);

  save=mpqs3_FB[2*mpqs3_nFB];
  mpqs3_FB[2*mpqs3_nFB]=0xffff;

#ifdef ASM_MPQS3_SIEVE
//zeitA(20);
//zeitA(72);
#ifdef ASM_MPQS3_SIEVE0
  asm3_sievea();
#else
  asm3_sieve();
#endif
//zeitB(72);
//zeitB(20);
#else
#if 1
  pb=mpqs3_sievelen>>2;
  i=mpqs3_sievebegin; fb+=2*i; fbs+=2*i;
  for (;;) {
    p=*fb; if (p>pb) break;
    lo=fblog[i]; sv=mpqs3_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++; 
    svend=sv+mpqs3_sievelen-4*p;
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
  pb=mpqs3_sievelen/3;
  svend=mpqs3_sievearray+mpqs3_sievelen;
  for (;;) {
    p=*fb; if (p>pb) break;
    lo=fblog[i]; sv=mpqs3_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }
  pb=mpqs3_sievelen>>1;
  for (;;) {
    p=*fb; if (p>pb) break;
    lo=fblog[i]; sv=mpqs3_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++; 
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }
  pb=mpqs3_sievelen;
  for (;;) {
    p=*fb; if (p>=pb) break;
    lo=fblog[i]; sv=mpqs3_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }

  for (;;) {
    p=*fb; if (p==0xffff) break;
    lo=fblog[i]; sv=mpqs3_sievearray;
    fb+=2;
    s1=*fbs++; s2=*fbs++;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
    i++;
  }
#else
  lo=0;
  fbs+=2*mpqs3_sievebegin;
  for (;;) {
    p=*fb0++;
    if (!p) {
      if (!*fb0) break;
      lo=(uchar)(*fb0++);
      p=*fb0++;
    }
    sv=mpqs3_sievearray;
    s1=*fbs++; s2=*fbs++; 
    svend=sv+mpqs3_sievelen-4*p;
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
  }
  fb0++;
  for (;;) {
    p=*fb0++;
    if (!p) {
      if (!*fb0) break;
      lo=(uchar)(*fb0++);
      p=*fb0++;
    }
    sv=mpqs3_sievearray;
    s1=*fbs++; s2=*fbs++; 
    svend=sv+mpqs3_sievelen-4*p;
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
  }
  fb0++;
  svend=mpqs3_sievearray+mpqs3_sievelen;
  for (;;) {
    p=*fb0++;
    if (!p) {
      if (!*fb0) break;
      lo=(uchar)(*fb0++);
      p=*fb0++;
    }
    sv=mpqs3_sievearray;
    s1=*fbs++; s2=*fbs++;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
  }
  fb0++;
  for (;;) {
    p=*fb0++;
    if (!p) {
      if (!*fb0) break;
      lo=(uchar)(*fb0++);
      p=*fb0++;
    }
    sv=mpqs3_sievearray;
    s1=*fbs++; s2=*fbs++; 
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
  }
  fb0++;
  for (;;) {
    p=*fb0++;
    if (!p) {
      if (!*fb0) break;
      lo=(uchar)(*fb0++);
      p=*fb0++;
    }
    sv=mpqs3_sievearray;
    s1=*fbs++; s2=*fbs++;
    sv[s1]+=lo; sv[s2]+=lo; sv+=p;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
  }
  fb0++;

  for (;;) {
    p=*fb0++;
    if (!p) {
      if (!*fb0) break;
      lo=(uchar)(*fb0++);
      p=*fb0++;
    }
    sv=mpqs3_sievearray;
    s1=*fbs++; s2=*fbs++;
    if (sv+s1<svend) sv[s1]+=lo;
    if (sv+s2<svend) sv[s2]+=lo;
  }
#endif

#endif
  mpqs3_FB[2*mpqs3_nFB]=save;

  sv=mpqs3_sievearray;
#ifdef MPQS3_ZEIT
//zeitA(6);
#endif
//zeitA(73);
  for (i=0; i<mpqs3_nAdiv_total; i++)
    if (!mpqs3_Adiv_active[i]) {
      p=mpqs3_Adiv_all[i]; lo=mpqs3_Adiv_log[i];
      s1=mpqs3_Adiv_start1[i];
      s2=mpqs3_Adiv_start2[i];
#if 0
      while (s1<mpqs3_sievelen) { sv[s1]+=lo; s1+=p; }
      while (s2<mpqs3_sievelen) { sv[s2]+=lo; s2+=p; }
#else
      if ((p<<1)<mpqs3_sievelen) {
        svend=sv+mpqs3_sievelen-2*p;
        while (sv<svend) {
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
      } else {
        while (s1<mpqs3_sievelen) { sv[s1]+=lo; s1+=p; }
        while (s2<mpqs3_sievelen) { sv[s2]+=lo; s2+=p; }
      }
      sv=mpqs3_sievearray;
#endif
    } else {
      p=mpqs3_Adiv_all[i]; lo=mpqs3_Adiv_log[i];
      s1=mpqs3_Adiv_start1[i];
      while (s1<mpqs3_sievelen) { sv[s1]+=lo; s1+=p; }
    }
//zeitB(73);
#ifdef MPQS3_ZEIT
//zeitB(6);
#endif
//zeitA(74);
  for (i=mpqs3_sievebegink; i<mpqs3_nFBk; i++) {
    p=mpqs3_FBk[i];
    lo=mpqs3_FBk_log[i];
    s1=mpqs3_FBk_start[i];
#if 1
    sv=mpqs3_sievearray;
    svend=sv+mpqs3_sievelen-4*p;
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
    while (s1<mpqs3_sievelen) { sv[s1]+=lo; s1+=p; }
#endif
  }
//zeitB(74);
}


static ushort mpqs3_evaluate()
{
  u32_t *ulsv, *ulsvend, h;
  uchar *sv;
  ushort **rels, buffer[256];
  u32_t i, nmaxsurv;

  mpqs3_nsurvivors=0;
  nmaxsurv=MPQS3_MAX_NRELS-mpqs3_nrels;
  if (MPQS3_MAX_NPRELS-mpqs3_nprels<nmaxsurv)
    nmaxsurv=MPQS3_MAX_NPRELS-mpqs3_nprels;
  if (MPQS3_MAX_NCRELS-mpqs3_ncrels<nmaxsurv)
    nmaxsurv=MPQS3_MAX_NCRELS-mpqs3_ncrels;
  if (nmaxsurv<MPQS3_MIN_RELBUFFER) return 0;
  if (nmaxsurv>255) nmaxsurv=255;  /* 1 weniger als buffer-Laenge !!! */
  ulsv=(u32_t *)mpqs3_sievearray;
  ulsvend=ulsv+(mpqs3_sievelen>>2);
  rels=mpqs3_rel_buffer;
#ifdef MPQS3_STAT
  stat_asm_eval++;
#endif
#ifdef ASM_MPQS3_EVAL
  mpqs3_nsurvivors=asm_evaluate0(ulsv,ulsvend,buffer,nmaxsurv);
#else
  while (ulsv<ulsvend) {
    h=*ulsv;
    if (h&0x80808080) {
      sv=(uchar *)ulsv;
      for (i=0; i<4; i++)
        if (*sv++&0x80) {
          buffer[mpqs3_nsurvivors++]=(ushort)(sv-mpqs3_sievearray-1);
//          printf("(%u %u) ",(u32_t)buffer[mpqs3_nsurvivors-1],(u32_t)sv[-1]);
        }
      if (mpqs3_nsurvivors>nmaxsurv-8) {
        while (ulsv<ulsvend) *ulsv++=0;
        break;
      }
    }
    *ulsv++=0;
    h=*ulsv;
    if (h&0x80808080) {
      sv=(uchar *)ulsv;
      for (i=0; i<4; i++)
        if (*sv++&0x80) {
          buffer[mpqs3_nsurvivors++]=(ushort)(sv-mpqs3_sievearray-1);
//          printf("(%u %u) ",(u32_t)buffer[mpqs3_nsurvivors-1],(u32_t)sv[-1]);
        }
      if (mpqs3_nsurvivors>nmaxsurv-8) {
        while (ulsv<ulsvend) *ulsv++=0;
        break;
      }
    }
    *ulsv++=0;
  }
#endif

  for (i=0; i<mpqs3_nsurvivors; i++) {
    mpqs3_sievearray[buffer[i]]=(uchar)(i+1);
    rels[i][0]=buffer[i];
  }
  return mpqs3_nsurvivors;
}



static inline ushort mpqs3_hash_lp(u32_t lp)
{
  u32_t j;
  u32_t lp1;
  ushort lphash, nh;

  lphash=(ushort)(lp&MPQS3_NHASHMASK); lphash>>=1;
  lp1=lp>>MPQS3_NHASHBITS;
  nh=mpqs3_hash_table[lphash][0];
  if (nh>=MPQS3_HASH_LEN-1) { printf("."); return MPQS3_HASH_OVERFLOW; }
  for (j=1; j<=nh; j++)
    if (mpqs3_hash_table[lphash][j]==lp1) break;
  if (j>nh) { /* new lp */
    mpqs3_hash_table[lphash][j]=lp1;
    mpqs3_hash_table[lphash][0]++;
    mpqs3_nlp++;
    mpqs3_hash_index[lphash][j]=mpqs3_nprels;
    mpqs3_lp[mpqs3_nprels]=lp;
  }
  return mpqs3_hash_index[lphash][j];
}


static inline int mpqs3_hash_rel(u32_t axb0)
{
  u32_t j;
  u32_t hash, hash2, nh;

  hash=axb0&0x00ffffff;
  hash2=(ushort)(hash>>8); hash&=0xff;
  nh=mpqs3_rel_hash_table[hash][0];
  if (nh>=MPQS3_RELHASH_LEN-1) { printf(":"); return 1; }
  for (j=1; j<=nh; j++)
    if (mpqs3_rel_hash_table[hash][j]==hash2) return 1;
  mpqs3_rel_hash_table[hash][nh+1]=hash2;
  mpqs3_rel_hash_table[hash][0]++;
  return 0;
}


static void mpqs3_decompose()
{
  u32_t i, j;
  ushort *fb=mpqs3_FB, *fbs=mpqs3_FB_start;
  uchar *sv, *svend, i1, i2, i3, i4;
  ushort ind, p, s1, s2, nr, d, lpi, r, ii;
  ushort **rels;
  u32_t axb[3], qx[3];
  double qx_dbl, axb_dbl, pr_dbl;
  u64_t qx_64, axb_64, pr_64;
  ushort d2, ii2, minus, save;
  double dx;
  short x;
  u32_t *ultarg, *ulsrc, lp, ulqx;
  u32_t inv, h, hi, ls1, ls2;
  u32_t ax, ay, az, at, asmh;
  u32_t rr, jj;
  u32_t *ulp, *ulphelp, pr[3];
  u32_t pend;

#ifdef MPQS3_ZEIT
zeitA(11);
#endif
  rels=mpqs3_rel_buffer-1;
  for (i=1; i<=mpqs3_nsurvivors; i++) rels[i][6]=0;
  i=mpqs3_td_begin; fb+=2*i; fbs+=2*i;

  save=mpqs3_FB[2*mpqs3_nFB];
  mpqs3_FB[2*mpqs3_nFB]=0xffff;

#ifdef ASM_MPQS3_TDSIEVE
//zeitA(42);
{
#error ASM_MPQS3_TDSIEVE
  ushort *buf[256];
  u32_t n;

  for (ii=1; ii<=mpqs3_nsurvivors; ii++) buf[ii]=rels[ii]+7+rels[ii][6];
  mpqs3_sievearray[mpqs3_sievelen]=0;
  n=asm3_tdsieve(fb,fbs,&(buf[0]),mpqs3_nFBk_1+i-1);
  fb+=n; fbs+=n;
  for (ii=1; ii<=mpqs3_nsurvivors; ii++)
    rels[ii][6]=(ushort)(buf[ii]-rels[ii]-7);
//zeitB(42);
}
#else
//zeitA(40);
  pend=mpqs3_sievelen/4;
  for (; i<mpqs3_nFB; i++) {
    p=*fb; if (p>pend) break;
    fb+=2; sv=mpqs3_sievearray; svend=sv+mpqs3_sievelen-2*p;
    s1=*fbs++; s2=*fbs++; i1=0;
    while (sv<svend) {
      i1=0;
      i1|=sv[s1]; i1|=sv[s2]; sv+=p;
      i1|=sv[s1]; i1|=sv[s2]; sv+=p;
      if (i1) {
        sv-=(2*p);
        i1=sv[s1]; i2=sv[s2]; sv+=p;
        if (i1||i2) {
          if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
          if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        }
        i1=sv[s1]; i2=sv[s2]; sv+=p;
        if (i1||i2) {
          if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
          if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        }
      }
    }
    svend+=p;
    if (sv<svend) {
      i1=sv[s1]; i1|=sv[s2]; sv+=p;
      if (i1) {
        sv-=p;
        i1=sv[s1]; i2=sv[s2]; sv+=p;
        if (i1||i2) {
          if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
          if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        }
      }
    }
    svend+=p; i1=0; i2=0;
    if (sv+s1<svend) i1=sv[s1];
    if (sv+s2<svend) i2=sv[s2];

    if (i1||i2) {
      if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
      if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
    }
  }
//zeitB(40);

//zeitA(41);
  pend=mpqs3_sievelen/3;
  for (; i<mpqs3_nFB; i++) {
    p=*fb; if (p>pend) break;
    fb+=2; sv=mpqs3_sievearray; svend=sv+mpqs3_sievelen;
    s1=*fbs++; s2=*fbs++;
    i1=sv[s1]; i2=sv[s2]; sv+=p;
    i1|=sv[s1]; i2|=sv[s2]; sv+=p;
    i1|=sv[s1]; i2|=sv[s2]; sv+=p;
    if (sv+s1<svend) i1|=sv[s1];
    if (sv+s2<svend) i2|=sv[s2];

    if (i1||i2) {
      if (i1) {
        sv-=(3*p); i1=sv[s1]; sv+=p;
        if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
        i1=sv[s1]; sv+=p;
        if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
        i1=sv[s1]; sv+=p;
        if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
        if (sv+s1<svend) {
          i1=sv[s1];
          if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
        }
      }
      if (i2) {
        sv-=(3*p); i2=sv[s2]; sv+=p;
        if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        i2=sv[s2]; sv+=p;
        if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        i2=sv[s2]; sv+=p;
        if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        if (sv+s2<svend) {
          i2=sv[s2];
          if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        }
      }
    }
  }
//zeitB(41);

//zeitA(42);
  pend=mpqs3_sievelen/2;
  for (; i<mpqs3_nFB; i++) {
    p=*fb; if (p>pend) break;
    fb+=2; sv=mpqs3_sievearray; svend=sv+mpqs3_sievelen;
    s1=*fbs++; s2=*fbs++;
    i1=sv[s1]; i2=sv[s2]; sv+=p;
    i1|=sv[s1]; i2|=sv[s2]; sv+=p;
    if (sv+s1<svend) i1|=sv[s1];
    if (sv+s2<svend) i2|=sv[s2];

    if (i1||i2) {
      if (i1) {
        sv-=(p+p); i1=sv[s1]; sv+=p;
        if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
        i1=sv[s1]; sv+=p;
        if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
        if (sv+s1<svend) {
          i1=sv[s1];
          if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
        }
      }
      if (i2) {
        sv-=(p+p); i2=sv[s2]; sv+=p;
        if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        i2=sv[s2]; sv+=p;
        if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        if (sv+s2<svend) {
          i2=sv[s2];
          if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
        }
      }
    }
  }
//zeitB(42);

//zeitA(43);
  pend=mpqs3_sievelen;
  for (; i<mpqs3_nFB; i++) {
    p=*fb; if (p>pend) break;
    fb+=2; sv=mpqs3_sievearray; svend=sv+mpqs3_sievelen;
    s1=*fbs++; s2=*fbs++; i3=0; i4=0;
    i1=sv[s1]; i2=sv[s2]; sv+=p;
    if (sv+s1<svend) i3=sv[s1];
    if (sv+s2<svend) i4=sv[s2];

    if (i1||i2||i3||i4) {
      if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
      if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
      if (i3) {nr=rels[i3][6];rels[i3][nr+7]=mpqs3_nFBk_1+i;rels[i3][6]++;}
      if (i4) {nr=rels[i4][6];rels[i4][nr+7]=mpqs3_nFBk_1+i;rels[i4][6]++;}
    }
  }
//zeitB(43);

//zeitA(44);
  for (; i<mpqs3_nFB; i++) {
    p=*fb; if (p==0xffff) break;
    fb+=2; sv=mpqs3_sievearray; svend=sv+mpqs3_sievelen;
    s1=*fbs++; s2=*fbs++; i1=0; i2=0;
    if (sv+s1<svend) i1=sv[s1];
    if (sv+s2<svend) i2=sv[s2];

    if (i1||i2) {
      if (i1) {nr=rels[i1][6];rels[i1][nr+7]=mpqs3_nFBk_1+i;rels[i1][6]++;}
      if (i2) {nr=rels[i2][6];rels[i2][nr+7]=mpqs3_nFBk_1+i;rels[i2][6]++;}
    }
  }
 if (*fb!=0xffff) complain("s %u %u\n",i,*fb);
//zeitB(44);

#if 0
  for (; i<mpqs3_nFB; i++) {
    p=*fb;
    fb+=2; sv=mpqs3_sievearray; svend=sv+mpqs3_sievelen-2*p;
    s1=*fbs++; s2=*fbs++;
    while (sv<svend) {
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][6]; rels[i1][nr+7]=mpqs3_nFBk_1+i; rels[i1][6]++; }
        if (i2) { nr=rels[i2][6]; rels[i2][nr+7]=mpqs3_nFBk_1+i; rels[i2][6]++; }
      }
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][6]; rels[i1][nr+7]=mpqs3_nFBk_1+i; rels[i1][6]++; }
        if (i2) { nr=rels[i2][6]; rels[i2][nr+7]=mpqs3_nFBk_1+i; rels[i2][6]++; }
      }
    }
    svend+=p;
    if (sv<svend) {
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][6]; rels[i1][nr+7]=mpqs3_nFBk_1+i; rels[i1][6]++; }
        if (i2) { nr=rels[i2][6]; rels[i2][nr+7]=mpqs3_nFBk_1+i; rels[i2][6]++; }
      }
    }
    svend+=p;
    if (sv+s1<svend) {
      i1=sv[s1];
      if (i1) { nr=rels[i1][6]; rels[i1][nr+7]=mpqs3_nFBk_1+i; rels[i1][6]++; }
    }
    if (sv+s2<svend) {
      i2=sv[s2];
      if (i2) { nr=rels[i2][6]; rels[i2][nr+7]=mpqs3_nFBk_1+i; rels[i2][6]++; }
    }
  }
#endif
#endif
  mpqs3_FB[2*mpqs3_nFB]=save;
#ifdef MPQS3_ZEIT
zeitB(11);
#endif
#ifdef MPQS3_STAT
 stat_td_cand+=mpqs3_nsurvivors;
#endif
  for (i=1; i<=mpqs3_nsurvivors; i++) {
    nr=rels[i][6];
    ind=rels[i][0];
    x=(short)(ind)-mpqs3_disp;

    axb_64=mpqs3_A_64*((u64_t)x)+mpqs3_B_64;
    axb_dbl=mpqs3_A_dbl*((double)x)+mpqs3_B_dbl;
    if (axb_dbl<0.) { axb_dbl=-axb_dbl; axb_64=-axb_64; }
    if (mpqs3_hash_rel((u32_t)axb_64)) goto next;
    mpqs3_convert(axb,axb_dbl,axb_64);

    qx_64=mpqs3_A_64*((u64_t)x)+mpqs3_2B_64;
    qx_64*=((u64_t)x);
    qx_64+=mpqs3_C_64;

    minus=0; dx=(double)x;
    qx_dbl=mpqs3_A_dbl*dx+mpqs3_2B_dbl;
    qx_dbl*=dx;
    qx_dbl+=mpqs3_C_dbl;
    if (qx_dbl<0.) { minus=1; qx_64=-qx_64; qx_dbl=-qx_dbl; }
    qx[0]=(u32_t)qx_64;
    qx[1]=(u32_t)(qx_64>>32);
//    qx_dbl-=(double)qx_64;
    qx[2]=(u32_t)(0.5+(qx_dbl-(double)qx_64)*DBL_CONVERT);

#ifdef MPQS3_STAT
stat_asm_td++;
#endif

#ifdef ASM_MPQS3_TD
#ifdef MPQS3_ZEIT
zeitA(12); zeitA(29);
#endif
#ifndef A64_STYLE_TD
    if (asm3_td(rels[i],minus,qx,&ulqx)) {
#ifdef MPQS3_ZEIT
 zeitB(12); zeitB(29);
#endif
      goto next;
    }
#else
    if((ulqx=asm3_td(rels[i],minus,qx))==0) {
#ifdef MPQS3_ZEIT
      zeitB(12); zeitB(29);
#endif
      goto next;
    }
#endif
#ifdef MPQS3_ZEIT
zeitB(29);
#endif
#else
    for (j=1; j<mpqs3_td_begin; j++) {
      p=mpqs3_FB[2*j];
      inv=mpqs3_FB_inv[j];
      ls1=ind+(p-mpqs3_FB_start[2*j]); ls2=ind+(p-mpqs3_FB_start[2*j+1]);
      ls1*=inv; ls2*=inv;
      ls1&=0xffff0000; ls2&=0xffff0000;
      if (!ls1 || !ls2) {
        rels[i][7+nr]=mpqs3_nFBk_1+j; nr++; /*printf("%d.",p);*/
      }
    }

    pr_dbl=1.; pr_64=1;
    for (j=0; j<nr; j++) {
      pr_dbl*=(double)(mpqs3_FB[2*(rels[i][7+j]-mpqs3_nFBk_1)]);
      pr_64*=(u64_t)(mpqs3_FB[2*(rels[i][7+j]-mpqs3_nFBk_1)]);
    }
    rels[i][6]=nr;

    if (minus) {
      rels[i][7+rels[i][6]]=0; rels[i][6]++;
    }
    while (!(qx[0]&1)) {
      if (rels[i][6]<MPQS3_TD_MAX_NDIV) {
        rels[i][7+rels[i][6]]=mpqs3_nFBk_1; rels[i][6]++;
        qx[0]=(qx[0]>>1)|(qx[1]<<31);
        qx[1]=(qx[1]>>1)|(qx[2]<<31);
        qx[2]>>=1;
        qx_dbl/=2.;
      } else {
#ifdef MPQS3_ZEIT
 zeitB(12);
#endif
        goto next;
      }
    }

#if 0
    mpqs3_convert(pr,pr_dbl,pr_64);
    if (pr[2]==0) {
      if (qx[2]>=pr[1]) {
        if (qx[2]>pr[1]) {
#ifdef MPQS3_ZEIT
zeitB(12);
#endif
          goto next;
        }
        if (qx[1]>=pr[0]) {

/*          printf(",");*/
#ifdef MPQS3_ZEIT
zeitB(12);
#endif
          goto next;
        }
      }
    }
#else
    if (qx_dbl>4294967295.5L*pr_dbl) {
#ifdef MPQS3_ZEIT
zeitB(12);
#endif
      goto next;
    }
    pr[0]=(u32_t)pr_64;
#endif

    ax=qx[0]; ay=pr[0]; /*az=0;*/
    inv=(u32_t)mpqs3_256_inv_table[(ay&0xff)>>1];
    at=ay*inv; at&=0x0000ff00;
    at*=inv; inv-=at;
    at=ay*inv; at&=0xffff0000;
    at*=inv; inv-=at;
    az=ax*inv;
    ulqx=az;

    for (j=0; j<nr; j++) {
      ii=rels[i][7+j];
      p=mpqs3_FB[2*(ii-mpqs3_nFBk_1)];
      while (ulqx%(u32_t)p==0) {
        if (rels[i][6]<MPQS3_TD_MAX_NDIV) {
          rels[i][7+rels[i][6]]=ii; rels[i][6]++;
          ulqx/=(u32_t)p;
        } else {
#ifdef MPQS3_ZEIT
 zeitB(12);
#endif
          goto next;
        }
      }
    }

    if (ulqx>1)
      for (j=0; j<mpqs3_nFBk; j++) {
        p=mpqs3_FBk[j];
        while (ulqx%(u32_t)p==0) {
          if (rels[i][6]<MPQS3_TD_MAX_NDIV) {
            rels[i][7+rels[i][6]]=1+j; rels[i][6]++;
            ulqx/=(u32_t)p;
          } else {
#ifdef MPQS3_ZEIT
 zeitB(12);
#endif
            goto next;
          }
        }
        if (ulqx==1) break;
      }
    if (ulqx>1)
      for (j=0; j<mpqs3_nAdiv_total; j++) {
        p=mpqs3_Adiv_all[j];
        while (ulqx%(u32_t)p==0) {
          if (rels[i][6]<MPQS3_TD_MAX_NDIV) {
            rels[i][7+rels[i][6]]=1+mpqs3_nFB+mpqs3_nFBk+j; rels[i][6]++;
            ulqx/=(u32_t)p;
          } else {
#ifdef MPQS3_ZEIT
 zeitB(12);
#endif
            goto next;
          }
        }
        if (ulqx==1) break;
      }
#endif

    if (rels[i][6]<MPQS3_TD_MAX_NDIV-mpqs3_nAdiv) {
      nr=rels[i][6]; ind=0;
      for (j=0; j<mpqs3_nAdiv_total; j++)
        if (mpqs3_Adiv_active[j]) {
          rels[i][7+nr+ind]=1+mpqs3_nFB+mpqs3_nFBk+j;
          ind++;
        }
      rels[i][6]+=mpqs3_nAdiv;
    } else {
#ifdef MPQS3_ZEIT
 zeitB(12);
#endif
      goto next;
    }
#ifdef MPQS3_ZEIT
zeitB(12);
#endif
    if (ulqx==1) {
#ifdef MPQS3_STAT
      stat_td_surv++;
#endif
      ulp=(u32_t *)(mpqs3_relations[mpqs3_nrels]);
      *ulp++=axb[0]; *ulp++=axb[1]; *ulp++=axb[2];
      nr=rels[i][6];
#ifdef USE_MEMCPY
      memcpy(mpqs3_relations[mpqs3_nrels]+6,rels[i]+6,(1+nr)*sizeof(ushort));
#else
      for (j=0; j<=nr; j++)
        mpqs3_relations[mpqs3_nrels][j+6]=rels[i][j+6];
#endif
      mpqs3_nrels++;
      goto next;
    }
//    if ((ulqx>mpqs3_pmax*mpqs3_pmax) && (ulqx<(1UL<<30))) {
/*
    if (ulqx>mpqs3_pmax*mpqs3_pmax) {
      mpz_set_ui(mpqs3_sq1,ulqx);
      if (psp(mpqs3_sq1)==0) stat_pp++;
    }
*/
#if 1
    if (ulqx<1048576) { /*printf("%d ",x);*/  /* !!! */
#else
    if (ulqx<4194304) { /*printf("%d ",x);*/  /* !!! */
#endif
      if (ulqx<=mpqs3_pmax) {
        gmp_fprintf(stderr,"error in mpqs_decompose for %Zd: %u %u %u   %u\n",
                    mpqs3_N,qx[0],qx[1],qx[2],ulqx);
        goto next;
        printf("%d  %u %u %u  %u ",x,qx[0],qx[1],qx[2],ulqx);
        complain("dec.2\n");
      }
      lp=ulqx;
      lpi=mpqs3_hash_lp(lp);
      if (lpi==MPQS3_HASH_OVERFLOW) { /*printf(".");*/ goto next; }
#ifdef MPQS3_STAT
      stat_td_surv++;
#endif
      if (lpi<mpqs3_nprels) {
        if (mpqs3_lp[lpi]!=lp) Schlendrian("lp");
//        mpqs3_lp[mpqs3_nprels]=lp;

        mpqs3_comb_index[2*mpqs3_ncrels]=lpi;
        mpqs3_comb_index[2*mpqs3_ncrels+1]=mpqs3_nprels;
        mpqs3_ncrels++;
      }
      ulp=(u32_t *)(mpqs3_part_rels[mpqs3_nprels]);
      *ulp++=axb[0]; *ulp++=axb[1]; *ulp++=axb[2];
      nr=rels[i][6];
#ifdef USE_MEMCPY
      memcpy(mpqs3_part_rels[mpqs3_nprels]+6,rels[i]+6,
             (1+nr)*sizeof(ushort));
#else
      for (j=0; j<=nr; j++)
        mpqs3_part_rels[mpqs3_nprels][j+6]=rels[i][j+6];
#endif
      mpqs3_nprels++;
    } /*else printf("<");*/
next: ;
  }
  if (mpqs3_nsp<mpqs3_nrels+mpqs3_ncrels)
    mpqs3_excess=mpqs3_nrels+mpqs3_ncrels-mpqs3_nsp;
  else mpqs3_excess=0;
}


/*
order of the primes:
-1, primes in mpqs3_FBk, primes in mpqs3_FB, primes in mpqs3_Adiv_total
*/


static int mpqs3_matrix_init()
{
  u32_t i, ii, j, n64, i64;
  ushort nr, pi, r;
  u64_t mask;

#ifdef MPQS3_ZEIT
zeitA(8);
#endif
  mpqs3_gauss_n=mpqs3_nrels+mpqs3_ncrels;
  if (mpqs3_gauss_n>=MPQS3_GAUSS_MAX) return -1;
  mpqs3_gauss_m=mpqs3_nsp;
  if (mpqs3_gauss_m>=mpqs3_gauss_n) Schlendrian("gauss: no excess\n");
/* build matrix */
  mpqs3_gauss_n64=(mpqs3_gauss_n+63)/64;
  memset(mpqs3_gauss_row64,0,mpqs3_gauss_m*mpqs3_gauss_n64*sizeof(u64_t));
  for (i=0; i<mpqs3_nrels; i++) {
    nr=mpqs3_relations[i][6];
    mask=1ULL<<(i&63); i64=i>>6;
    for (j=0; j<nr; j++) {
      pi=mpqs3_relations[i][7+j];
      mpqs3_gauss_row64[pi+mpqs3_gauss_m*i64]^=mask;
    }
  }
  for (i=0; i<mpqs3_ncrels; i++) {
    ii=i+mpqs3_nrels; mask=1ULL<<(ii&63); i64=ii>>6;

    r=mpqs3_comb_index[2*i];
    nr=mpqs3_part_rels[r][6];
    for (j=0; j<nr; j++) {
      pi=mpqs3_part_rels[r][7+j];
      mpqs3_gauss_row64[pi+mpqs3_gauss_m*i64]^=mask;
    }

    r=mpqs3_comb_index[2*i+1];
    nr=mpqs3_part_rels[r][6];
    for (j=0; j<nr; j++) {
      pi=mpqs3_part_rels[r][7+j];
      mpqs3_gauss_row64[pi+mpqs3_gauss_m*i64]^=mask;
    }
  }
  mpqs3_gauss_k=0;
  for (i=0; i<mpqs3_gauss_n; i++) mpqs3_gauss_d[i]=-1;
#ifdef MPQS3_ZEIT
zeitB(8);
#endif
  return 0;
}


//SMJS Changed to void to stop compiler complaint, no return statement but called as void static int mpqs3_rowechelon()
static void mpqs3_rowechelon()
{
  u32_t i, j, k, l, t, k64, col, row, r0, c0, zz;
  u64_t tmp, tab[16];
  uchar ucm[1024], mm[4], m, z;  /* TODO replace 1024 */

#ifdef MPQS3_ZEIT
zeitA(8);
#endif
  col=0; row=0;
  while (1) {
    if (row>=mpqs3_gauss_m) break;
    if (col>=mpqs3_gauss_n) break;

    k64=col>>6;
    r0=row; c0=col;
    for (i=0; i<4; i++,col++) {
      for (j=row; j<mpqs3_gauss_m; j++) {
        m=(uchar)(0xf&mpqs3_gauss_row64[j+mpqs3_gauss_m*k64]>>(c0&63));
        if (!m) continue;
        for (k=0,z=1; k<i; k++,z+=z) if (m&z) m^=mm[k];
        if (m&z) break;
      }
      if (j>=mpqs3_gauss_m) {
        mpqs3_gauss_d[col]=-1;
        continue;
      }
      if (j>row) {
        for (l=k64; l<mpqs3_gauss_n64; l++) {
          tmp=mpqs3_gauss_row64[j+mpqs3_gauss_m*l];
          mpqs3_gauss_row64[j+mpqs3_gauss_m*l]=
            mpqs3_gauss_row64[row+mpqs3_gauss_m*l];
          mpqs3_gauss_row64[row+mpqs3_gauss_m*l]=tmp;
        }
      }
      m=(uchar)(mpqs3_gauss_row64[row+mpqs3_gauss_m*k64]>>(c0&63));
      for (k=0,z=1; k<i; k++,z+=z) {
        if (m&z) {
          m^=mm[k];
          for (l=k64; l<mpqs3_gauss_n64; l++)
            mpqs3_gauss_row64[row+mpqs3_gauss_m*l]^=
                mpqs3_gauss_row64[mpqs3_gauss_d[c0+k]+mpqs3_gauss_m*l];
        }
      }
      mm[i]=m;
      for (k=0; k<i; k++) {
        if (mpqs3_gauss_d[c0+k]==-1) continue;
        if (mm[k]&z) {
          mm[k]^=mm[i];
          for (l=k64; l<mpqs3_gauss_n64; l++)
            mpqs3_gauss_row64[mpqs3_gauss_d[c0+k]+mpqs3_gauss_m*l]^=
                mpqs3_gauss_row64[row+mpqs3_gauss_m*l];
        }
      }
      mpqs3_gauss_d[col]=row;
      row++;
      if (row>=mpqs3_gauss_m) break;
    }
    for (j=0; j<mpqs3_gauss_m; j++)
      ucm[j]=(uchar)((mpqs3_gauss_row64[j+mpqs3_gauss_m*k64]>>(c0&63))&0xf);
    for (j=r0; j<row; j++) ucm[j]=0;
#ifdef ASM_MPQS3_GAUSS
    for (l=k64; l<mpqs3_gauss_n64; l++)
      asm_re_strip(mpqs3_gauss_row64+mpqs3_gauss_m*l,mpqs3_gauss_m,
                   mpqs3_gauss_d+c0,ucm);
#else
    for (l=k64; l<mpqs3_gauss_n64; l++) {
      tab[0]=0ULL;
      for (j=0,zz=1; j<4; j++,zz+=zz) {
        if (mpqs3_gauss_d[c0+j]==-1) tab[zz]=0ULL;
        else tab[zz]=mpqs3_gauss_row64[mpqs3_gauss_d[c0+j]+mpqs3_gauss_m*l];
        for (k=1; k<zz; k++) tab[zz+k]=tab[k]^tab[zz];
      }
      for (t=0; t<mpqs3_gauss_m; t++)
        mpqs3_gauss_row64[t+mpqs3_gauss_m*l]^=tab[ucm[t]];
    }
#endif
  }
  mpqs3_gauss_k=0;
#ifdef MPQS3_ZEIT
zeitB(8);
#endif
}


static int mpqs3_matrix()
{
  u32_t i, k64;
  int notzero;
  u64_t mask;

  for (i=0; i<mpqs3_gauss_n; i++) mpqs3_sol[i]=0;
  for (; mpqs3_gauss_k<mpqs3_gauss_n; mpqs3_gauss_k++) {
    if (mpqs3_gauss_d[mpqs3_gauss_k]>=0) continue;
    notzero=0;
    mask=1ULL<<(mpqs3_gauss_k&63); k64=mpqs3_gauss_k>>6;
    for (i=0; i<mpqs3_gauss_k; i++)
      if (mpqs3_gauss_d[i]!=-1)
        if (mpqs3_gauss_row64[mpqs3_gauss_d[i]+mpqs3_gauss_m*k64] & mask) {
          mpqs3_sol[i]=1;
          notzero=1;
        }

    if (notzero) {
      mpqs3_sol[mpqs3_gauss_k]=1;
      mpqs3_gauss_k++;
      return 1;
    }
  }
  return 0;
}


static u32_t mpqs3_is_power(mpz_t op) /* very sloppy */
{
  u32_t res;

  for (res=2; res<12; res++)
    if (mpz_root(mpqs3_gmp_help,op,(ulong)res))
      return res;
  return 0;
}


static int mpqs3_final(size_t max_bits)
{
  u32_t i, j, k;
  ushort nr, nr1, nr2, r1, r2;
  u32_t p, lp, prod1, prod2;
  u32_t *ulp;
  int split;
  size_t modthres;

  modthres=mpz_sizeinbase(mpqs3_N,2)*4;
  if (!mpqs3_nfactors) {
    mpqs3_nfactors=1;
    mpz_set(mpqs3_factors[0],mpqs3_N);
    mpqs3_f_status[0]=0;
  }
  mpqs3_rowechelon();
  while (mpqs3_matrix()) {
#ifdef MPQS3_STAT
stat_mpqs_ntrials++;
#endif
    memset(mpqs3_exp,0,mpqs3_nsp*sizeof(short));

    mpz_set_ui(mpqs3_sq1,1); prod1=1;
    mpz_set_ui(mpqs3_sq2,1); prod2=1;

    for (j=0; j<mpqs3_nrels; j++)
      if (mpqs3_sol[j]) {
        ulp=(u32_t *)(mpqs3_relations[j]);
        mpz_set_v(mpqs3_gmp_help,ulp);
        nr=mpqs3_relations[j][6];
        if (j&1) {
          for (k=0; k<nr; k++) mpqs3_exp[mpqs3_relations[j][7+k]]--;
          mpz_mul(mpqs3_sq1,mpqs3_sq1,mpqs3_gmp_help);
          if (mpz_sizeinbase(mpqs3_sq1,2)>modthres)
            mpz_mod(mpqs3_sq1,mpqs3_sq1,mpqs3_N);
        } else {
          for (k=0; k<nr; k++) mpqs3_exp[mpqs3_relations[j][7+k]]++;
          mpz_mul(mpqs3_sq2,mpqs3_sq2,mpqs3_gmp_help);
          if (mpz_sizeinbase(mpqs3_sq2,2)>modthres)
            mpz_mod(mpqs3_sq2,mpqs3_sq2,mpqs3_N);
        }
      }

    for (j=0; j<mpqs3_ncrels; j++)
      if (mpqs3_sol[mpqs3_nrels+j]) {
        r1=mpqs3_comb_index[2*j];
        r2=mpqs3_comb_index[2*j+1];

        nr1=mpqs3_part_rels[r1][6];
        nr2=mpqs3_part_rels[r2][6];
        for (k=0; k<nr1; k++) mpqs3_exp[mpqs3_part_rels[r1][7+k]]++;
        for (k=0; k<nr2; k++) mpqs3_exp[mpqs3_part_rels[r2][7+k]]--;

        ulp=(u32_t *)(mpqs3_part_rels[r1]);
        mpz_set_v(mpqs3_gmp_help,ulp);
        mpz_mul(mpqs3_sq2,mpqs3_sq2,mpqs3_gmp_help);
        if (mpz_sizeinbase(mpqs3_sq2,2)>modthres)
          mpz_mod(mpqs3_sq2,mpqs3_sq2,mpqs3_N);
        ulp=(u32_t *)(mpqs3_part_rels[r2]);
        mpz_set_v(mpqs3_gmp_help,ulp);
        mpz_mul(mpqs3_sq1,mpqs3_sq1,mpqs3_gmp_help);
        if (mpz_sizeinbase(mpqs3_sq1,2)>modthres)
          mpz_mod(mpqs3_sq1,mpqs3_sq1,mpqs3_N);
      }

    for (j=0; j<mpqs3_nsp; j++)
      if (mpqs3_exp[j]&1) Schlendrian("final.odd %u\n",j);
      else mpqs3_exp[j]>>=1;

    for (j=1; j<mpqs3_nsp; j++)
      if (mpqs3_exp[j]) {
        if (j<1+mpqs3_nFBk) {
          p=(u32_t)(mpqs3_FBk[j-1]);
        } else if (j<1+mpqs3_nFB+mpqs3_nFBk) {
          k=j-mpqs3_nFBk-1;
          p=(u32_t)(mpqs3_FB[2*k]);
        } else {
          k=j-mpqs3_nFB-mpqs3_nFBk-1;
          p=(u32_t)(mpqs3_Adiv_all[k]);
        }
        if (mpqs3_exp[j]>0) {
          for (k=0; k<mpqs3_exp[j]; k++) {
            prod1*=p;
            if (prod1&0xfffc0000) {  /* p<16384 */
              mpz_mul_ui(mpqs3_sq1,mpqs3_sq1,(ulong)prod1);
              if (mpz_sizeinbase(mpqs3_sq1,2)>modthres) {
                mpz_mod(mpqs3_sq1,mpqs3_sq1,mpqs3_N);
#ifdef MPQS3_STAT
stat_final_mulmod++;
#endif
              }
              prod1=1;
            }
          }
        } else {
          for (k=0; k<-mpqs3_exp[j]; k++) {
            prod2*=p;
            if (prod2&0xfffc0000) {  /* p<16384 */
              mpz_mul_ui(mpqs3_sq2,mpqs3_sq2,(ulong)prod2);
              if (mpz_sizeinbase(mpqs3_sq2,2)>modthres) {
                mpz_mod(mpqs3_sq2,mpqs3_sq2,mpqs3_N);
#ifdef MPQS3_STAT
stat_final_mulmod++;
#endif
              }
              prod2=1;
            }
          }
        }
      }
    if (prod1>1) mpz_mul_ui(mpqs3_sq1,mpqs3_sq1,(ulong)prod1);
    if (prod2>1) mpz_mul_ui(mpqs3_sq2,mpqs3_sq2,(ulong)prod2);

    mpz_add(mpqs3_gmp_help,mpqs3_sq1,mpqs3_sq2);
    mpz_gcd(mpqs3_gmp_help,mpqs3_gmp_help,mpqs3_N);
    if (mpz_cmp_ui(mpqs3_gmp_help,1)==0) continue;
    if (mpz_cmp(mpqs3_gmp_help,mpqs3_N)==0) continue;

    for (j=0; j<mpqs3_nfactors; j++)
      if (mpqs3_f_status[j]==0) {
        mpz_gcd(mpqs3_sq1,mpqs3_gmp_help,mpqs3_factors[j]);
        if (mpz_cmp_ui(mpqs3_sq1,1)==0) continue;
        if (mpz_cmp(mpqs3_sq1,mpqs3_factors[j])==0) continue;
        mpz_set(mpqs3_factors[mpqs3_nfactors],mpqs3_sq1);
        mpz_divexact(mpqs3_factors[j],mpqs3_factors[j],mpqs3_sq1);
#if 1
        mpqs3_f_status[j]=psp(mpqs3_factors[j]);
        if (mpqs3_f_status[j])
          if (mpz_sizeinbase(mpqs3_factors[j],2)>max_bits)
            return 1;
        mpqs3_f_status[mpqs3_nfactors]=psp(mpqs3_factors[mpqs3_nfactors]);
#else
        mpqs3_f_status[j]=mpz_probab_prime_p1(mpqs3_factors[j],1);
        if (mpqs3_f_status[j])
          if (mpz_sizeinbase(mpqs3_factors[j],2)>max_bits)
            return 1;
        mpqs3_f_status[mpqs3_nfactors]=mpz_probab_prime_p1(mpqs3_factors[mpqs3_nfactors],1);
#endif
        mpqs3_nfactors++;
        if (mpqs3_f_status[mpqs3_nfactors-1])
          if (mpz_sizeinbase(mpqs3_factors[mpqs3_nfactors-1],2)>max_bits)
            return 1;
        break;
      }

    split=1;
    for (j=0; j<mpqs3_nfactors; j++) if (mpqs3_f_status[j]==0) split=0;
    if (split) return 1;
    if (mpqs3_nfactors>=MPQS3_MAX_FACTORS)
      Schlendrian("final: too many factors\n");
  }
  for (j=0; j<mpqs3_nfactors; j++)
    while (mpqs3_f_status[j]==0) {
      if ((k=mpqs3_is_power(mpqs3_factors[j]))) {
        if (mpqs3_nfactors+k-1>MPQS3_MAX_FACTORS)
          Schlendrian("final: too many factors\n");
        if (!mpz_root(mpqs3_factors[j],mpqs3_factors[j],(ulong)k))
          Schlendrian("error in power detection\n");
        mpqs3_f_status[j]=psp(mpqs3_factors[j]);
        for (k--; k; k--) {
          mpz_set(mpqs3_factors[mpqs3_nfactors],mpqs3_factors[j]);
          mpqs3_f_status[mpqs3_nfactors]=mpqs3_f_status[j];
          mpqs3_nfactors++;
        }
      } else return 0;
    }
  if (j>=mpqs3_nfactors) return 1;
  return 0;
}


static int mpqs3_factor0(mpz_t N, size_t max_bits, mpz_t **factors, int retry)
{
  size_t nbits;
  int err;
  u32_t i, ev;

#ifdef MPQS3_ZEIT
zeitA(9);
#endif
#ifdef MPQS3_STAT
  stat_mpqs_nsieves=0; stat_mpqs_nsurvivors=0; stat_mpqs_ntrials=0; stat_mpqs_ndiv=0;
  if (retry) stat_retry++;
#endif
  if (!mpqs3_isinit) mpqs3_init();
#ifdef MPQS3_ZEIT
zeitB(9);
#endif
  nbits=mpz_sizeinbase(N,2);
  if (nbits>148) {
    fprintf(stderr,"warning: mpqs3 called with >148 Bit\n");
    return -2;
  }
  if (nbits<80) {
    fprintf(stderr,"warning: mpqs3 called with <80 Bit\n");
    return -2;
  }
#ifdef MPQS3_ZEIT
zeitA(1);
#endif
  mpz_set(mpqs3_N,N);
  mpqs3_choose_multiplier();
  mpqs3_choose_parameter(retry);
  mpqs3_generate_FB();
  if ((err=mpqs3_SI_init())<1) {
    printf("%d ",err);
    fprintf(stderr,"warning: mpqs3: error in self-initialization ");
    mpz_out_str(stderr,10,N); printf("\n");
#ifdef MPQS3_ZEIT
zeitB(1);
#endif
    return -3;
  }
#ifdef MPQS3_ZEIT
zeitB(1);
#endif
  while (1) { // printf("%u,%u; ",stat_mpqs_nsieves,mpqs3_nrels);
//printf("S: %u %u %u %u %u %u %u \n",stat_mpqs_nsieves,stat_mpqs_nsurvivors,stat_mpqs_ntrials,stat_mpqs_ndiv,mpqs3_nrels,mpqs3_nprels,mpqs3_ncrels);
#ifdef MPQS3_ZEIT
zeitA(2);
#endif
    if (!mpqs3_next_pol()) {
      if (retry) {
        fprintf(stderr,"warning: not enough polynomials in mpqs3 with ");
        mpz_out_str(stderr,10,N); fprintf(stderr,"\n");
      }
#ifdef MPQS3_ZEIT
zeitB(2);
#endif
      return -4;
    }
#ifdef MPQS3_ZEIT
zeitB(2); zeitA(3);
#endif
    mpqs3_sieve();
#ifdef MPQS3_STAT
stat_mpqs_nsieves++;
#endif
#ifdef MPQS3_ZEIT
zeitB(3);
zeitA(4);
#endif
    ev=mpqs3_evaluate();
#ifdef MPQS3_STAT
stat_mpqs_nsurvivors+=ev;
#endif
#ifdef MPQS3_ZEIT
zeitB(4);
#endif
    if (ev) {
#ifdef MPQS3_ZEIT
zeitA(5);
#endif
      mpqs3_decompose();
#ifdef MPQS3_ZEIT
zeitB(5);
#endif
      if (mpqs3_nsurvivors && (mpqs3_excess>MPQS3_MIN_EXCESS)) {
#ifdef MPQS3_ZEIT
zeitA(14);
#endif
        if (mpqs3_matrix_init()) return -6;
#ifdef MPQS3_ZEIT
zeitB(14);
zeitA(10);
#endif
        if (mpqs3_final(max_bits)) {
#ifdef MPQS3_ZEIT
zeitB(10);
#endif
#ifdef MPQS3_STAT
// printf("Stat: %u %u %u %u %u %u %u \n",stat_mpqs_nsieves,stat_mpqs_nsurvivors,stat_mpqs_ntrials,stat_mpqs_ndiv,mpqs3_nrels,mpqs3_nprels,mpqs3_ncrels);
stat_ff+=mpqs3_nrels;
stat_pf+=mpqs3_nprels;
stat_comb+=mpqs3_ncrels;
#endif
          for (i=0; i<mpqs3_nfactors; i++)
            if (mpz_sizeinbase(mpqs3_factors[i],2)>max_bits) return 0;
          *factors=mpqs3_factors;
          return mpqs3_nfactors;
        }
#ifdef MPQS3_ZEIT
zeitB(10);
#endif
      }
    }
    if (mpqs3_nrels>MPQS3_MAX_NRELS-MPQS3_MIN_RELBUFFER) {
      fprintf(stderr,"warning: too many relations in mpqs3\n");
      return -5;
    }
  }
  return -1;
}


int mpqs3_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
  int err;

  err=mpqs3_factor0(N,max_bits,factors,0);
  if (err==-4) err=mpqs3_factor0(N,max_bits,factors,1);
  if (err==-4) err=mpqs3_factor0(N,max_bits,factors,3);

  return err;
}

