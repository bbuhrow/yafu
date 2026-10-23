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

/*
kmpqs1: lasieve5_64/mpqs.c (MPQS for N <= 96 bits) imported into yafu.
The asm kernels (asm/mpqs_*.asm) are replaced by lasieve's own C fallbacks;
timing hooks are dropped.  Configuration is lasieve's TINY variant.
*/

#include "kmpqs_common.h"
#include "kmpqs_kernels.h"

/* every function that can fail has the context in scope as 'qs' */
#define KMPQS_JMPBUF (qs->jmpbuf)

/* from asm/mpqs-config.h (TINY): nFB, sievebegin (unused in TINY), nAdiv,
   nAdiv_total, accept, td_begin, sievelen; one row per 4 bits, 44..96 */
static uint16_t mpqs_param[14][7]={
 { 30, 3, 2, 5, 7, 30, 768},      /* 44 */
 { 30, 3, 3, 5, 9, 30, 1024},     /* 48 */
 { 35, 3, 3, 5, 10, 40, 1280},    /* 52 */
 { 40, 3, 3, 6, 10, 50, 1536},    /* 56 */
 { 50, 3, 3, 6, 10, 50, 1536},    /* 60 */
 { 60, 3, 3, 6, 11, 60, 1536},    /* 64 */
 { 60, 3, 3, 7, 11, 60, 1792},    /* 68 */
 { 80, 3, 3, 7, 12, 80, 2048},    /* 72 */
 { 90, 3, 3, 7, 13, 90, 2048},    /* 76 */
 { 110, 3, 3, 7, 14, 110, 2560},   /* 80 */
 { 120, 3, 4, 7, 14, 120, 2560},   /* 84 */
 { 140, 3, 4, 7, 15, 140, 3072},   /* 88 */
 { 150, 3, 4, 7, 16, 150, 3072},   /* 92 */
 { 160, 4, 4, 8, 17, 160, 3072}    /* 96 */
};


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

/* all per-factorization state; one context per thread */
typedef struct {
  jmp_buf jmpbuf;
  uint16_t nFBk_1;
  uint16_t td_begin, sievebegin, sievebegink;
  uint16_t FB_inv_info[4*MPQS_MAX_NPRIMES];
  uint16_t nFBk, FBk[3]; /* 3*5*7*11 too big */
  uint8_t FB_log[MPQS_MAX_FBSIZE], FBk_log[3];
  uint16_t nFB, nAdiv_total;
  uint16_t Adiv_all[MPQS_MAX_ADIV_ALL];
  uint32_t FBk_inv[3], FB_A_inv[MPQS_MAX_ADIV_ALL];
  uint16_t Adiv_disp[MPQS_MAX_ADIV_ALL];
  uint16_t FB_mm_inv[MPQS_MAX_FBSIZE], FB_disp[MPQS_MAX_FBSIZE];
  uint16_t FB_np_p[2*MPQS_MAX_FBSIZE];
  uint16_t FB_np_px[2*MPQS_MAX_FBSIZE];
  mpz_t gmp_help, dummy, help;
  mpz_t N, kN, factors[MPQS_MAX_FACTORS];
  uint8_t f_status[MPQS_MAX_FACTORS];
  uint32_t multiplier, nfactors;
  double complexity;
  uint16_t prime_table[MPQS_MAX_NPRIMES];
  uint16_t prime_sqrthelp[MPQS_MAX_NPRIMES];
  uint16_t nAdiv, accept;
  uint8_t Adiv_log[MPQS_MAX_ADIV_ALL];;
  uint16_t pmax;
  uint16_t Adiv[MPQS_MAX_ADIV_ALL];
  uint16_t Adiv_sqrt[MPQS_MAX_ADIV_ALL], Adiv_all_sqrt[MPQS_MAX_ADIV_ALL];
  uint16_t Adiv_start1[MPQS_MAX_ADIV_ALL], Adiv_start2[MPQS_MAX_ADIV_ALL];
  uint8_t Adiv_active[MPQS_MAX_ADIV_ALL];
  uint16_t A_mask[(1<<MPQS_MAX_ADIV_ALL)];
  uint16_t nA;
  short A_index, B_index;
  uint16_t SI_inv_table[MPQS_MAX_ADIV_ALL][MPQS_MAX_FBSIZE];
  uint16_t Adiv_SI_inv_table[MPQS_MAX_ADIV_ALL][MPQS_MAX_ADIV_ALL];
  uint16_t SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_FBSIZE], SIk_add[MPQS_MAX_ADIV_ALL][3];
  uint16_t Adiv_SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_ADIV_ALL];
  long long A, B, C, twoB, Bi[MPQS_MAX_ADIV_ALL];
  uint32_t disp;
  uint16_t nsurvivors;
  uint32_t nrels, nlp, nsp, excess;
  uint8_t logbound;
  uint16_t nprels, ncrels;
  uint16_t **relations, **part_rels, **comb_rels;
  uint16_t **rel_buffer;
  uint16_t rel_hash_table[256][16];
  uint16_t hash_table[128][16];
  uint32_t hash_index[128][16];
  uint32_t lp[128*15];
  uint8_t *rel_active; /* auf bitarray umstellen !!! */
  uint16_t nrels_active, nlp_active;
  uint16_t lp_count[128*16+MPQS_MAX_FBSIZE+3+MPQS_MAX_ADIV_ALL];
  uint64_t *gauss_row64;  /* word-major, see kmpqs_rowechelon */
  uint8_t sol[MPQS_GAUSS_MAX];
  short gauss_d[MPQS_GAUSS_MAX];
  uint32_t gauss_m, gauss_n, gauss_k;
  uint32_t gauss_n64;
  short exp[128*16+MPQS_MAX_FBSIZE+3+MPQS_MAX_ADIV_ALL];
  mpz_t sq1, sq2;
  uint16_t zero1, zero2;
  double kN_dbl, maxval_dbl;
  uint64_t kN_64, A_inv_64;
  uint16_t FB[2*MPQS_MAX_FBSIZE], FB0[MPQS_MAX_FBSIZE+2*256];
  uint16_t FB_start[2*MPQS_MAX_FBSIZE], FBk_start[3];
  uint32_t sievelen;
  uint32_t FB_inv[MPQS_MAX_FBSIZE];
  uint8_t *sievearray, *tinyarray;
  uint32_t tiny_prod, ntiny;
  uint16_t sinit_tab[8*MPQS_SI_NVALUES];
  double multiplier_cand_log[4][MPQS_MULT_NCAND];
  char multiplier_cand_jacobi[4][MPQS_MULT_NCAND][MPQS_MULT_NTESTPR];
  uint32_t kN_w[8]; int kN_nw;  /* kN as 32-bit words, for kN mod p */
  /* for the kernels in kmpqs_kernels.h */
  uint16_t FB_p2[2*MPQS_MAX_FBSIZE+KMPQS_VPAD];
  uint16_t FB_inv2[2*MPQS_MAX_FBSIZE+KMPQS_VPAD];
  uint16_t SI_add2[MPQS_MAX_ADIV_ALL][2*MPQS_MAX_FBSIZE+KMPQS_VPAD];
  uint16_t td_list[MPQS_MAX_FBSIZE+KMPQS_VPAD];  /* +slack for vector stores */
  /* floor((2^32-1)/p): x*inv(p) mod 2^32 <= lim  <=>  p | x */
  uint32_t FB_lim[MPQS_MAX_FBSIZE], FBk_lim[3], FB_A_lim[MPQS_MAX_ADIV_ALL];
} mpqs_ctx_t;





static uint8_t mpqs_256_inv_table[128]={
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


static uint8_t mpqs_A_table[256]={
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

static uint32_t mpqs_A_table_n[10]={
0, 1, 9, 37, 93, 163, 219, 247, 255, 256
};

#ifdef MPQS_STAT
static uint32_t stat_mpqs_nsieves, stat_mpqs_nsurvivors, stat_mpqs_ntrials, stat_mpqs_ndiv;
#endif

static uint64_t mpqs_inv_64(uint64_t a)
{
  uint64_t inv, h;

  inv=(uint64_t)mpqs_256_inv_table[(a&0xff)>>1];
  h=a*inv; h&=0xff00ULL;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffff0000ULL;
  h*=inv; inv-=h;
  h=a*inv; h&=0xffffffff00000000ULL;
  h*=inv; inv-=h;
  /*if (inv*a!=1ULL) KMPQS_FAIL("mpqs_inv_64");*/
  return inv;
}


static uint8_t mpqs_jacobi_tab2[8]={0,0,0,1,0,1,0,0};

static int mpqs_jacobi(mpqs_ctx_t *qs, uint16_t a, uint16_t b) /* always gcd(a,b)=1 and 0<a<b */
{
  if (!(b&1)) KMPQS_FAIL("mpqs_jacobi: b even\n");
  return kmpqs_jacobi(a,b);
}


static uint16_t mpqs_invert(uint16_t a, uint16_t p)
{
  uint32_t v1=0, v2=1;
  uint32_t b=(uint32_t)a, pp=(uint32_t)p, q;

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
/*  if ((v2*a-1)%(uint32_t)p) KMPQS_FAIL("invert");*/
  return v2;
}


#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p32; h32=h32>>16


static uint16_t mpqs_powmod(mpqs_ctx_t *qs, uint16_t a, uint16_t e, uint16_t p)
{
  uint16_t ex=e, rrr;
  uint32_t aa, res, p32, inv, h, h1;

//rrr=mpqs_powmod0(a,e,p);
  if (!ex) return 1;
  aa=(uint32_t)a;
  p32=(uint32_t)p;

  inv=(uint32_t)mpqs_256_inv_table[(p&0xff)>>1];
  h=p32*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=p32*inv; h&=0xffff0000;
  h*=inv; inv-=h;
  if (inv*p32!=1) KMPQS_FAIL("mpqs_powmod_inv");
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
if (((uint16_t)res)!=rrr)
  KMPQS_FAIL("mpqs_powmod %u %u %u  %u %u\n",(uint32_t)a,(uint32_t)e,p32,(uint32_t)rrr,res);
*/
      return (uint16_t)res;
    }
    aa*=aa;
    h1=(aa&0x0000ffff)*inv; h1&=0x0000ffff; aa+=h1*p32; aa>>=16;
  }
  return 0; /* never reached */
}


static uint32_t mpqs_powmod_mm(uint32_t a, uint32_t e, uint32_t p32, uint32_t mmi, uint32_t one)
{
  uint32_t aa, res, ex=e, h32, hh32;

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


static uint16_t mpqs_sqrt_init(mpqs_ctx_t *qs, uint16_t p)
{
  uint16_t e, u, i, j;
  // SMJS initialised g to stop compiler warning
  uint32_t g = 0, b;

  if (p&2) return 0;
  if (p&4) return mpqs_powmod(qs, 2,(p-1)>>2,p);
  e=0; u=p-1;
  while (!(u&1)) { e++; u>>=1; }
  for (i=2; i<p; i++) {
    g=mpqs_powmod(qs, i,u,p);
    if (g==1) continue;
    b=g;
    for (j=0; j<e-1; j++) { b*=b; b%=(uint32_t)p; }
    if (b==p-1) break;
  }
  if (i>=p) KMPQS_FAIL("mpqs_sqrt_init\n");
  return g;
}


static uint16_t mpqs_sqrt(uint16_t a, uint16_t p, uint16_t help)
                 /* 0<a<p, p prime, p!=2, (a/p)=1 */
{
  uint32_t e, u, i, l;
  uint32_t b, g, r, k,aa;
  uint32_t p32, mmi, h, h32, hh32, one, mminit;

  p32=(uint32_t)p;
  mmi=(uint32_t)mpqs_256_inv_table[(p&0xff)>>1];
  h=p32*mmi; h&=0xff00;
  h*=mmi; mmi-=h;
  h=p32*mmi; h&=0xffff0000;
  h*=mmi; mmi-=h;
//  if (mmi*p32!=1) KMPQS_FAIL("mpqs_sqrt_inv");
  mmi=-mmi;
  mminit=(0xffffffff)%p32; mminit++; /* mminit!=p32 */

  aa=(uint32_t)a;
  h32=aa*mminit; MMREDUCE; aa=h32; if (aa>=p32) aa-=p32;
  h32=mminit; MMREDUCE; one=h32; if (one>=p32) one-=p32;

  if (p&2) {
    aa=mpqs_powmod_mm(aa,(p32+1)>>2,p32,mmi,one);
    h32=aa;
    MMREDUCE;
    if (h32>=p32) h32-=p32;
    return (uint16_t)h32;
  }
  if (p&4) {
    b=mpqs_powmod_mm(aa,(p32+3)>>3,p32,mmi,one);
    h32=b*b; MMREDUCE; g=h32; if (g>=p32) g-=p32;
    h32=b;
    if (g!=aa) h32*=help;
    MMREDUCE; b=h32; if (b>=p32) b-=p32;
    return (uint16_t)b;
  }
  e=0; u=p32-1;
  while (!(u&1)) { e++; u>>=1; }
  h32=((uint32_t)help)*mminit; MMREDUCE; g=h32;
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
  return (uint16_t)r;
}


static uint32_t mpqs_inv(mpqs_ctx_t *qs, uint16_t a)
{
  uint32_t inv, h;

  inv=(uint32_t)mpqs_256_inv_table[(a&0xff)>>1];
  h=(uint32_t)a*inv; h&=0xff00;
  h*=inv; inv-=h;
  h=(uint32_t)a*inv; h&=0xffff0000;
  h*=inv; inv-=h;
  if (inv*(uint32_t)a!=1ULL) KMPQS_FAIL("mpqs_inv");
  return inv;
}

/* ---------------------------------------------------------- */


static uint16_t mpqs_multiplier_cand[4][MPQS_MULT_NCAND]={
 { 1, 17, 33, 41, 57, 65, 73, 89 },
 { 3, 11, 19, 35, 43, 51, 59, 67 },
 { 5, 13, 21, 29, 37, 53, 61, 69 },
 { 7, 15, 23, 31, 39, 47, 55, 71 },
};



static void mpqs_init(mpqs_ctx_t *qs)
{
  uint32_t i, j, k;
  uint32_t p, add, d;

  mpz_init(qs->N);
  mpz_init(qs->kN);
  mpz_init(qs->dummy);
  mpz_init(qs->gmp_help);
  mpz_init(qs->help);
  mpz_init(qs->sq1);
  mpz_init(qs->sq2);
  for (i=0; i<MPQS_MAX_FACTORS; i++) mpz_init(qs->factors[i]);
  qs->sievearray=(uint8_t *)kmpqs_xmalloc((1+MPQS_SIEVELENMAX)*sizeof(uint8_t));
  qs->tinyarray=(uint8_t *)kmpqs_xmalloc(32*MPQS_MAX_TINYPROD*sizeof(uint8_t)+64);
  qs->relations=(uint16_t **)kmpqs_xmalloc(MPQS_MAX_NRELS*sizeof(uint16_t *));
  qs->relations[0]=(uint16_t *)kmpqs_xmalloc(MPQS_MAX_NRELS*MPQS_REL_ENTRIES*sizeof(uint16_t));
  for (i=1; i<MPQS_MAX_NRELS; i++)
    qs->relations[i]=qs->relations[i-1]+MPQS_REL_ENTRIES;
  qs->part_rels=(uint16_t **)kmpqs_xmalloc(MPQS_MAX_NPRELS*sizeof(uint16_t *));
  qs->part_rels[0]=(uint16_t *)kmpqs_xmalloc(MPQS_MAX_NPRELS*MPQS_REL_ENTRIES*sizeof(uint16_t));
  for (i=1; i<MPQS_MAX_NPRELS; i++)
    qs->part_rels[i]=qs->part_rels[i-1]+MPQS_REL_ENTRIES;
  qs->comb_rels=(uint16_t **)kmpqs_xmalloc(MPQS_MAX_NCRELS*sizeof(uint16_t *));
  qs->comb_rels[0]=(uint16_t *)kmpqs_xmalloc(MPQS_MAX_NCRELS*2*MPQS_REL_ENTRIES*sizeof(uint16_t));
  for (i=1; i<MPQS_MAX_NCRELS; i++)
    qs->comb_rels[i]=qs->comb_rels[i-1]+2*MPQS_REL_ENTRIES;
  qs->rel_buffer=(uint16_t **)kmpqs_xmalloc(MPQS_MAX_NRELS_BUF*sizeof(uint16_t *));
  qs->rel_buffer[0]=(uint16_t *)kmpqs_xmalloc(MPQS_MAX_NRELS_BUF*MPQS_REL_ENTRIES*sizeof(uint16_t));
  for (i=1; i<MPQS_MAX_NRELS_BUF; i++)
    qs->rel_buffer[i]=qs->rel_buffer[i-1]+MPQS_REL_ENTRIES;

  qs->rel_active=(uint8_t *)kmpqs_xmalloc(MPQS_MAX_NRELS*sizeof(uint8_t));
  qs->gauss_row64=(uint64_t *)kmpqs_xmalloc(MPQS_GAUSS_MAX*MPQS_GAUSS_MAX/64*sizeof(uint64_t));

  qs->prime_table[0]=2;
  qs->prime_table[1]=3;
  qs->prime_table[2]=5; qs->prime_sqrthelp[2]=mpqs_sqrt_init(qs, 5);
  i=3; p=7; add=4;
  while (i<MPQS_MAX_NPRIMES) {
    for (j=2; j<i; j++) {
      d=(uint32_t)qs->prime_table[j];
      if (d*d>=p) break;
      if (p%d==0) break;
    }
    if (d*d>p) {
      qs->prime_table[i]=(uint16_t)p;
      qs->prime_sqrthelp[i]=mpqs_sqrt_init(qs, p);
      i++;
    }
    p+=add; add=6-add;
  }

  for (i=0; i<4; i++)
    for (j=0; j<MPQS_MULT_NCAND; j++) {
      d=mpqs_multiplier_cand[i][j];
      for (k=1; k<MPQS_MULT_NTESTPR; k++) {
        p=(uint32_t)qs->prime_table[k];
        if (d%p) {
          qs->multiplier_cand_jacobi[i][j][k]=
            (char)(mpqs_jacobi(qs, (uint16_t)(d%p),(uint16_t)p));
        } else qs->multiplier_cand_jacobi[i][j][k]=0;
      }
      qs->multiplier_cand_log[i][j]=log((double)d)-1.5*log(2);
    }

}


static void mpqs_choose_multiplier(mpqs_ctx_t *qs)
{
  uint16_t Nmod8, mult, p, mm, r;
  uint16_t residue[MPQS_MULT_NTESTPR];
  double value[MPQS_MULT_NCAND], v, v0, vp, vmin, dn;
  double plog[MPQS_MULT_NTESTPR];
  char jac[MPQS_MULT_NTESTPR], c;
  uint32_t i, j;

  uint32_t Nw[8];
  int nNw=kmpqs_mpz_words(qs->N,Nw,8);

  for (i=1; i<MPQS_MULT_NTESTPR; i++) {
    p=qs->prime_table[i];
    r=(uint16_t)kmpqs_mod_words(Nw,nNw,p);
    jac[i]=(char)mpqs_jacobi(qs, r,p);
  }
  Nmod8=(uint16_t)(Nw[0]&7);
  dn=log(mpz_get_d(qs->N))/log(2.);
  dn=100.9-dn;
  if (dn>7) mm=128; else mm=(uint16_t)(exp(dn*log(2)));
  for (i=1; i<MPQS_MULT_NTESTPR; i++) {
    p=qs->prime_table[i];
    plog[i]=log((double)p)/((double)p);
  }
  for (j=0; j<MPQS_MULT_NCAND; j++) {
    mult=mpqs_multiplier_cand[Nmod8/2][j];
    if (mult<=mm) {
      v=qs->multiplier_cand_log[Nmod8/2][j];
      for (i=1; i<MPQS_MULT_NTESTPR; i++) {
        c=qs->multiplier_cand_jacobi[Nmod8/2][j][i];
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
  if (j>=MPQS_MULT_NCAND) KMPQS_FAIL("mpqs_choose_multiplier.vmin\n");
  qs->complexity=(vmin+log(mpz_get_d(qs->N)))/log(2.);
  qs->multiplier=mpqs_multiplier_cand[Nmod8/2][j];
#ifdef MPQS_STAT
// printf("%u ",mpqs_multiplier);
#endif
  mpz_mul_ui(qs->kN,qs->N,(unsigned long)qs->multiplier);
  qs->kN_nw=kmpqs_mpz_words(qs->kN,qs->kN_w,8);
  qs->kN_64=kmpqs_mpz_get_u64(qs->kN);  /* kN mod 2^64 */
}


static void mpqs_choose_parameter(mpqs_ctx_t *qs, int retry)
{
  int n, r;

  n=(int)(qs->complexity+0.5);
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
  qs->nFB=mpqs_param[n][0];
  if (r) qs->nFB+=r*qs->nFB/4;
  if (qs->nFB>=MPQS_MAX_FBSIZE) qs->nFB=MPQS_MAX_FBSIZE-1;
  qs->sievebegin=mpqs_param[n][1];
  qs->nAdiv=mpqs_param[n][2];
  qs->nAdiv_total=mpqs_param[n][3];
  qs->accept=mpqs_param[n][4];
  qs->td_begin=mpqs_param[n][5];

  qs->sievelen=mpqs_param[n][6];
  if (r) {
    if (qs->sievelen<MPQS_SIEVELENMAX>>r) qs->sievelen<<=r;
    else qs->sievelen=MPQS_SIEVELENMAX;
  }
  qs->disp=qs->sievelen/2;
  qs->nfactors=0; qs->nrels=0; qs->nlp=0; qs->nprels=0; qs->ncrels=0;
  for (n=0; n<128; n++) qs->hash_table[n][0]=0;
  for (n=0; n<256; n++) qs->rel_hash_table[n][0]=0;
}


static void mpqs_generate_FB(mpqs_ctx_t *qs)
{
  uint16_t *fb, p, rest, i, nfb;

  fb=qs->FB; nfb=0; i=0; qs->nFBk=0;
  *fb++=2; *fb++=1; nfb++;
  for (i=1; i<MPQS_MAX_NPRIMES; i++) {
    p=qs->prime_table[i];
    if (p>MPQS_FB_MAXPRIME) break;
    rest=(uint16_t)kmpqs_mod_words(qs->kN_w,qs->kN_nw,p);
    if (rest) {
//zeitA(25);
      if (mpqs_jacobi(qs, rest,p)==1) {
        *fb++=p;
        *fb++=mpqs_sqrt(rest,p,qs->prime_sqrthelp[i]);
        nfb++;
        if (nfb>=qs->nFB) { /*zeitB(25);*/ break; }
      }
//zeitB(25);
    } else {
      if (qs->multiplier%(uint32_t)p)
        KMPQS_FAIL("mpqs: N has small divisor: %u\n",p);
      qs->FBk[qs->nFBk++]=p;
    }
  }
  nfb-=((nfb-qs->nAdiv_total)&3);
  qs->nFB=nfb;
  qs->nFBk_1=qs->nFBk+1;
  qs->nsp=1+qs->nFB+qs->nFBk;
  qs->pmax=qs->FB[2*(nfb-1)];
  qs->FB[2*nfb]=qs->sievelen;
  qs->FB[2*nfb]=0xffff;
}


static int mpqs_SI_init(mpqs_ctx_t *qs)
{
  double d;
  uint32_t a, i, j, k, l, n;
  long long prod;
  uint16_t inv, p, *fb=qs->FB, pk, tp_max;
  double A_value[20], A_div_log[MPQS_MAX_ADIV_ALL]; /* maximal number of A's */
  double v, vmax;
  uint8_t lo;
  uint16_t pb;

  d=mpz_get_d(qs->kN); qs->kN_dbl=d;
  d*=8.; d=sqrt(d); d/=(double)qs->sievelen;
  d=log(d); d/=(double)qs->nAdiv; d=exp(d);
  a=(uint32_t)d;
/*  if (a>=mpqs_pmax) { printf("%u %u ",a,mpqs_pmax); return -1; }*/
  if (a>=qs->pmax) a=qs->pmax-1;
  for (i=0; i<qs->nFB; i++) if (qs->FB[2*i]>a) break;
  if (i>qs->nAdiv_total/2) i-=qs->nAdiv_total/2; else i=1;
  if (i<1) i=1; /* first prime >2 */
  if (qs->FB[2*i]<3) return -2;
  if (i+qs->nAdiv_total>qs->nFB) {
    if (qs->nFB<=qs->nAdiv_total+1) return -3;
    i=qs->nFB-qs->nAdiv_total-1;
  }
  for (j=0; j<qs->nAdiv_total; j++) {
    qs->Adiv_all[j]=qs->FB[2*(i+j)];
    qs->Adiv_all_sqrt[j]=qs->FB[2*(i+j)+1];
  }
  for (j=i+qs->nAdiv_total; j<qs->nFB; j++) {
    qs->FB[2*(j-qs->nAdiv_total)]=qs->FB[2*j];
    qs->FB[2*(j-qs->nAdiv_total)+1]=qs->FB[2*j+1];
  }
  qs->nFB-=qs->nAdiv_total;
/* tiny */
  qs->sievebegin=1; qs->sievebegink=0; qs->tiny_prod=1; qs->ntiny=0;
  tp_max=qs->sievelen>>4;
  if (tp_max>MPQS_MAX_TINYPROD) tp_max=MPQS_MAX_TINYPROD;
  while (1) {
    if (qs->sievebegin<qs->nFB) p=qs->FB[2*qs->sievebegin];
    else p=tp_max+1;
    if (qs->sievebegink<qs->nFBk) pk=qs->FBk[qs->sievebegink];
    else pk=tp_max+1;
    if (p<pk) {
      if (qs->tiny_prod*(uint32_t)p>=tp_max) break;
      qs->sievebegin++; qs->ntiny++;
      qs->tiny_prod*=(uint32_t)p;
    } else {
      if (qs->tiny_prod*(uint32_t)pk>=tp_max) break;
      qs->sievebegink++; qs->ntiny++;
      qs->tiny_prod*=(uint32_t)pk;
    }
  }
#ifdef MPQS_STAT
//for (j=0; j<mpqs_sievebegin; j++) printf("%u ",mpqs_FB[2*j]);
//printf("\n");
#endif

/*for (j=0; j<mpqs_sievebegin; j++) printf("%u ",mpqs_FB[2*j]);
printf("\n");*/
  if (qs->nFB<qs->td_begin) {
    if (qs->nFB&3) KMPQS_FAIL("mpqs_nFB\n");
//printf("%u\n",mpqs_nFB);
    qs->td_begin=qs->nFB;
  }

  if (qs->td_begin&1) qs->td_begin++;
  p=qs->FB[2];
  qs->FB_inv[1]=mpqs_inv(qs, p);
  qs->FB_inv_info[2]=p; qs->FB_inv_info[3]=p;
  qs->FB_inv_info[6]=(uint16_t)qs->FB_inv[1];
  qs->FB_inv_info[7]=(uint16_t)qs->FB_inv[1];
  if (qs->td_begin&1) KMPQS_FAIL("mpqs_td_begin is odd\n");
  for (j=2; j<qs->td_begin; j+=2) {
    p=qs->FB[2*j];
    qs->FB_inv[j]=mpqs_inv(qs, p);
    qs->FB_inv_info[4*j]=p; qs->FB_inv_info[4*j+1]=p;
    qs->FB_inv_info[4*j+4]=(uint16_t)qs->FB_inv[j];
    qs->FB_inv_info[4*j+5]=(uint16_t)qs->FB_inv[j];
    p=qs->FB[2*(j+1)];
    qs->FB_inv[j+1]=mpqs_inv(qs, p);
    qs->FB_inv_info[4*j+2]=p; qs->FB_inv_info[4*j+3]=p;
    qs->FB_inv_info[4*j+6]=(uint16_t)qs->FB_inv[j+1];
    qs->FB_inv_info[4*j+7]=(uint16_t)qs->FB_inv[j+1];
  }
  for (; j<qs->nFB; j++) qs->FB_inv[j]=mpqs_inv(qs, qs->FB[2*j]);
  for (j=0; j<qs->nFBk; j++) qs->FBk_inv[j]=mpqs_inv(qs, qs->FBk[j]);
  for (j=0; j<qs->nAdiv_total; j++)
    qs->FB_A_inv[j]=mpqs_inv(qs, qs->Adiv_all[j]);

/* for next_pol: */
{
/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/

  uint32_t mmi, h32, hh32, invh, p32;
  uint16_t invhelp[MPQS_MAX_ADIV_ALL], hhh;

  for (j=1; j<qs->nFB; j++) {
    p=qs->FB[2*j]; p32=(uint32_t)p;
    mmi=-qs->FB_inv[j];
    invhelp[0]=qs->Adiv_all[0]%p;
    h32=(uint32_t)(invhelp[0]);
    for (i=1; i<qs->nAdiv_total; i++) {
      h32*=(uint32_t)(qs->Adiv_all[i]);
      MMREDUCE;
      invhelp[i]=(uint16_t)h32;
    }
    MMREDUCE; if (h32>=p) h32-=p; /* divide by 2^16 */
    invh=(uint32_t)mpqs_invert((uint16_t)h32,p);
    invh<<=1;

    for (i=qs->nAdiv_total-1; i; i--) {
      h32=invh; h32*=(uint32_t)(invhelp[i-1]);
      MMREDUCE; if (h32>=p) h32-=p;
      qs->SI_inv_table[i][j]=(uint16_t)h32;
      h32=invh; h32*=(uint32_t)qs->Adiv_all[i];
      MMREDUCE; invh=h32;
    }
    if (invh>=p) invh-=p;
    qs->SI_inv_table[0][j]=(uint16_t)invh;
  }
  for (j=0; j<qs->nAdiv_total; j++) {
    uint32_t h;
    uint16_t pa;

    p=qs->Adiv_all[j];
    for (i=0; i<qs->nAdiv_total; i++) {
      pa=qs->Adiv_all[i]%p;
      if (pa) {
        h=(uint32_t)mpqs_invert(pa,p); //h*=2; if (h>=p) h-=p;
        h<<=17; h%=((uint32_t)p);
      } else h=0;
      qs->Adiv_SI_inv_table[i][j]=(uint16_t)h;
    }
    qs->Adiv_disp[j]=(uint16_t)(qs->disp%((uint32_t)p));
  }
}

  for (j=1; j<qs->nFB; j++) {
    qs->FB_mm_inv[j]=-qs->FB_inv[j];
    p=qs->FB[2*j];
    qs->FB_disp[j]=(uint16_t)(qs->disp%((uint32_t)p));
  }
  for (j=0; j<qs->nFB; j+=4) {
    qs->FB_np_p[2*j]=qs->FB[2*j];
    qs->FB_np_p[2*j+1]=qs->FB[2*j+2];
    qs->FB_np_p[2*j+2]=qs->FB[2*j+4];
    qs->FB_np_p[2*j+3]=qs->FB[2*j+6];
    qs->FB_np_p[2*j+4]=qs->FB[2*j+1];
    qs->FB_np_p[2*j+5]=qs->FB[2*j+3];
    qs->FB_np_p[2*j+6]=qs->FB[2*j+5];
    qs->FB_np_p[2*j+7]=qs->FB[2*j+7];
  }

/* compute log-approximations */
  d=mpz_get_d(qs->kN);
  d/=8.; d=sqrt(d); d*=(double)qs->sievelen;
  qs->maxval_dbl=d;
  d=log(d);
  d-=log((double)(1<<qs->accept));
  qs->logbound=128;  /* ??? */
  d=(double)(qs->logbound)/d;
  for (i=0; i<qs->nFB; i++) {
    p=qs->FB[2*i];
    qs->FB_log[i]=(uint8_t)(0.5+d*log((double)p));
  }
  for (i=0; i<qs->nFBk; i++) {
    p=qs->FBk[i];
    qs->FBk_log[i]=(uint8_t)(0.5+d*log((double)p));
  }
  for (i=0; i<qs->nAdiv_total; i++) {
    p=qs->Adiv_all[i];
    qs->Adiv_log[i]=(uint8_t)(0.5+d*log((double)p));
  }

/* fill mpqs_FB0 */
  lo=0; n=0; i=qs->sievebegin;
  pb=qs->sievelen>>3;
  for (; i<qs->nFB; i++) {
    p=qs->FB[2*i]; if (p>pb) break;
    if (lo<qs->FB_log[i]) {
      qs->FB0[n++]=0;
      lo=qs->FB_log[i];
      qs->FB0[n++]=(uint16_t)lo;
    }
    qs->FB0[n++]=p;
  }
  qs->FB0[n++]=0; qs->FB0[n++]=0;
  for (j=4; j; j--) {
    pb=qs->sievelen/j;
    for (; i<qs->nFB; i++) {
      p=qs->FB[2*i]; if (p>pb) break;
      if (lo<qs->FB_log[i]) {
        qs->FB0[n++]=0;
        lo=qs->FB_log[i];
        qs->FB0[n++]=(uint16_t)lo;
      }
      qs->FB0[n++]=p;
    }
    qs->FB0[n++]=0; qs->FB0[n++]=0;
  }
  for (; i<qs->nFB; i++) {
    p=qs->FB[2*i];
    if (lo<qs->FB_log[i]) {
      qs->FB0[n++]=0;
      lo=qs->FB_log[i];
      qs->FB0[n++]=(uint16_t)lo;
    }
    qs->FB0[n++]=p;
  }
  qs->FB0[n++]=0; qs->FB0[n++]=0;

{
/* experimental */
  double *alog;
  uint16_t tmp;
  uint32_t *asort;

/* compute mask for choice of A-divisors */
  j=0;
  if (qs->nAdiv_total>8) {
    k=0; if (qs->nAdiv>8) k=qs->nAdiv-8;
    for (; k<=qs->nAdiv; k++) {
      uint32_t i0, i1;
      uint16_t m16;

      for (i0=mpqs_A_table_n[k]; i0<mpqs_A_table_n[k+1]; i0++) {
        m16=(uint16_t)(mpqs_A_table[i0]); m16<<=8;
        if (m16>>qs->nAdiv_total) break;
        for (i1=mpqs_A_table_n[qs->nAdiv-k];
             i1<mpqs_A_table_n[qs->nAdiv-k+1]; i1++) {
          qs->A_mask[j]=m16+(uint16_t)(mpqs_A_table[i1]);
          j++;
        }
      }
    }
  } else {
    for (i=mpqs_A_table_n[qs->nAdiv];
         i<mpqs_A_table_n[qs->nAdiv+1]; i++) {
      uint16_t m16;

      m16=(uint16_t)(mpqs_A_table[i]);
      if (m16>>qs->nAdiv_total) break;
      qs->A_mask[j]=m16;
      j++;
    }
  }
/* sort mpqs_A_mask */
  alog=(double *)kmpqs_xmalloc(j*sizeof(double));

/* ensure that |C|<2^96 */
  for (i=0; i<qs->nAdiv_total; i++)
    A_div_log[i]=log((double)(qs->Adiv_all[i]))/log(2.);
  d=log(qs->kN_dbl)/log(2.);
  for (i=0; i<j; i++) {
    v=0;
    for (k=0; k<qs->nAdiv_total; k++)
      if (qs->A_mask[i]&(1<<k)) v+=A_div_log[k];
    alog[i]=v;
    if (d-v>62.95) {
      j--; /* printf(":");*/
      qs->A_mask[i]=qs->A_mask[j];
      i--;
    }
  }

  d=mpz_get_d(qs->kN); qs->kN_dbl=d;
  d*=8.; d=sqrt(d); d/=(double)qs->sievelen;
  d=log(d)/log(2.);
  for (i=0; i<j; i++)
    if (d>alog[i]) alog[i]=d-alog[i];
    else alog[i]-=d;

{
  uint32_t *a0, *a1, *atmp, zp, i0, i1, n0, n1;

  asort=(uint32_t *)kmpqs_xmalloc(2*j*sizeof(*asort));
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
  for (i=0; i<j; i++) a1[i]=(uint32_t)(qs->A_mask[a1[i]]);
  for (i=0; i<j; i++) qs->A_mask[i]=(uint16_t)(a1[i]);
  free(asort);
}
  free(alog);
}

  qs->nA=j; qs->A_index=-1; qs->B_index=(1<<qs->nAdiv)-1;

  qs->FB[2*qs->nFB]=qs->sievelen;
  qs->FB[2*qs->nFB]=0xffff;

/* duplicated primes and 16-bit inverses for the kernels */
  memset(qs->FB_p2,0,sizeof(qs->FB_p2));
  memset(qs->FB_inv2,0,sizeof(qs->FB_inv2));
  for (j=1; j<qs->nFB; j++) {
    qs->FB_p2[2*j]=qs->FB_p2[2*j+1]=qs->FB[2*j];
    qs->FB_inv2[2*j]=qs->FB_inv2[2*j+1]=(uint16_t)qs->FB_inv[j];
  }
  for (j=1; j<qs->nFB; j++) qs->FB_lim[j]=0xffffffffU/qs->FB[2*j];
  for (j=0; j<qs->nFBk; j++) qs->FBk_lim[j]=0xffffffffU/qs->FBk[j];
  for (j=0; j<qs->nAdiv_total; j++) qs->FB_A_lim[j]=0xffffffffU/qs->Adiv_all[j];
/* kmpqs_td_roots needs ind+p < 2^16; lasieve's parameters keep it far below */
  if (qs->sievelen+qs->pmax>=65536) KMPQS_FAIL("sievelen+pmax too large\n");
  return 1;
}


static int mpqs_next_pol(mpqs_ctx_t *qs)
{
  uint32_t i, ind, j;
  long long add, bi;
  uint16_t p, s1, s2, *fb=qs->FB, *fbs=qs->FB_start;
  uint32_t mod, inv;
  short sh;
  uint32_t bb, cc, c1, c2;
  uint16_t mask;
  double d;
  uint16_t aa;
  uint32_t al, *ul;
  uint32_t B_modp;

  qs->B_index++;
  if (qs->B_index>=1<<(qs->nAdiv-1)) {
    uint32_t Bi_mul[MPQS_MAX_ADIV_ALL], Adiv_nr[MPQS_MAX_ADIV_ALL];

    qs->A_index++;
    if (qs->A_index>=qs->nA) return 0;
    mask=qs->A_mask[qs->A_index];
    for (i=0; i<qs->nAdiv_total; i++)
      if (mask & (1<<i)) qs->Adiv_active[i]=1;
      else qs->Adiv_active[i]=0;
    j=0;
    for (i=0; i<qs->nAdiv_total; i++)
      if (qs->Adiv_active[i]) {
        qs->Adiv[j]=qs->Adiv_all[i];
        qs->Adiv_sqrt[j]=qs->Adiv_all_sqrt[i];
        Adiv_nr[j]=i;
        j++;
      }
    qs->A=1;
    for (i=0; i<qs->nAdiv; i++) qs->A*=(long long)(qs->Adiv[i]);
/* precomputations for sieve init */
    qs->sinit_tab[0]=0; qs->sinit_tab[1]=0;
    i=1;
    for (j=1; j<MPQS_SI_NVALUES; j++) {
      d=sqrt(qs->kN_dbl+qs->maxval_dbl*((double)(qs->A>>j)))/((double)qs->A);
      if (d<(double)qs->disp) s1=qs->disp-(uint16_t)d; else s1=0;
      s1-=(s1&63);
      if (s1>qs->sinit_tab[2*i-2]) {
        qs->sinit_tab[2*i]=s1; qs->sinit_tab[2*i+1]=j;
        i++;
      }
    }
    for (j=1; j<MPQS_SI_NVALUES; j++) {
      d=qs->kN_dbl-qs->maxval_dbl*((double)(qs->A>>(MPQS_SI_NVALUES-j)));
      if (d<0.) break;
      d=sqrt(d)/((double)qs->A);
      if (d<(double)qs->disp) s1=qs->disp-(uint16_t)d; else s1=0;
      s1-=(s1&63);
      if (s1>qs->sinit_tab[2*i-2]) {
        qs->sinit_tab[2*i]=s1; qs->sinit_tab[2*i+1]=MPQS_SI_NVALUES-1-j;
        i++;
      }
    }

    if (qs->sinit_tab[2*i-2]==qs->disp) i--;
    for (j=i; j<2*i; j++) {
      qs->sinit_tab[2*j]=qs->sievelen-qs->sinit_tab[2*(2*i-1-j)];
      if (j<2*i-1) qs->sinit_tab[2*j+1]=qs->sinit_tab[2*(2*i-2-j)+1];
    }
/* compute B_i's */
    for (i=0; i<qs->nAdiv; i++) {
      p=qs->Adiv[i];
      bi=qs->A/(long long)(p);
      inv=(uint32_t)(bi%(long long)p); /* bi>0 */
      inv=mpqs_invert((uint16_t)inv,p);
      mod=(inv*((uint32_t)qs->Adiv_sqrt[i]))%((uint32_t)p);
      Bi_mul[i]=mod;
      bi*=(long long)mod;
      qs->Bi[i]=bi;  /* >=0 */
    }
    qs->B=0;
    for (i=0; i<qs->nAdiv; i++) qs->B+=qs->Bi[i];

/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/
{
/* roots and root increments of all FB primes (lasieve's asm_next_pol10/11) */
    const uint16_t *inv_rows[MPQS_MAX_ADIV_ALL];
    uint16_t *add_rows[MPQS_MAX_ADIV_ALL], *add2_rows[MPQS_MAX_ADIV_ALL];

    for (j=0; j<qs->nAdiv; j++) {
      inv_rows[j]=qs->SI_inv_table[Adiv_nr[j]];
      add_rows[j]=qs->SI_add[j];
      add2_rows[j]=qs->SI_add2[j];
    }
    kmpqs_np_newA(qs->nFB,fb,qs->FB_mm_inv,qs->FB_disp,qs->nAdiv,inv_rows,
                  add_rows,add2_rows,Bi_mul,1u<<(16-qs->nAdiv),fbs);
}

    qs->twoB=2*qs->B;
    qs->A_inv_64=mpqs_inv_64(qs->A);
    /* C=(B^2-kN)/A, exact mod 2^64 (unsigned: the products wrap) */
    qs->C=(long long)(((uint64_t)qs->B*(uint64_t)qs->B-qs->kN_64)*qs->A_inv_64);

/* divisors of A */
    for (i=0; i<qs->nAdiv_total; i++)
      if (!qs->Adiv_active[i]) {
/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/

        uint32_t mmi, h32, hh32, pi, jj, p32;

        p=qs->Adiv_all[i]; p32=(uint32_t)p;
        mmi=-mpqs_inv(qs, p);

        bb=0; pi=1<<(16-qs->nAdiv);
        for (j=0; j<qs->nAdiv; j++) {
          jj=Adiv_nr[j];
          cc=qs->Adiv_SI_inv_table[jj][i];
          cc*=Bi_mul[j]; h32=cc;
          MMREDUCE; cc=h32; if (cc>=p) cc-=p;
          qs->Adiv_SI_add[j][i]=(uint16_t)cc;
          bb+=qs->Adiv_SI_add[j][i]; if (bb>=p) bb-=p;
          pi*=((uint32_t)qs->Adiv_SI_inv_table[jj][i]); h32=pi;
          MMREDUCE; pi=h32;
        }
        cc=bb;
        if (cc&1) cc+=p; cc>>=1;
        cc=qs->Adiv_disp[i]+(p-cc); if (cc>=p) cc-=p;
        c1=qs->Adiv_all_sqrt[i];
        h32=c1*pi; MMREDUCE; c1=h32;
        if (c1>=p) c1-=p;
        c2=p-c1;
        c1+=cc; if (c1>=p) c1-=p;
        c2+=cc; if (c2>=p) c2-=p;
        qs->Adiv_start1[i]=(uint16_t)c1; qs->Adiv_start2[i]=(uint16_t)c2;
      } else {
/*
#define MMREDUCE hh32=(h32&0x0000ffff)*mmi; hh32&=0x0000ffff; \
   h32+=hh32*p; h32=h32>>16
*/
// SMJS Initialise j0 to stop compiler warning, not sure about this one

        uint32_t mmi, h32, hh32, jj, cc0, j0 = 0, p32;

        p=qs->Adiv_all[i]; p32=(uint32_t)p;
        sh=(short)(qs->twoB%(long long)p);
        if (sh>0) bb=(uint32_t)sh; else bb=(uint32_t)((short)p+sh);
        bb=(uint32_t)mpqs_invert(bb,p);
        sh=(short)(qs->C%(long long)p); sh=-sh;
        if (sh>0) cc=(uint32_t)sh; else cc=(uint32_t)((short)p+sh);

        cc*=bb; cc%=(uint32_t)p;  /* now p|2*B*cc+C */
        cc0=2*cc; if (cc0>=p) cc0-=p;
        cc+=qs->Adiv_disp[i]; if (cc>=p) cc-=p;
        qs->Adiv_start1[i]=(uint16_t)cc;

        mmi=-mpqs_inv(qs, p);
        for (j=0; j<qs->nAdiv; j++) {
          jj=Adiv_nr[j];
          if (jj!=i) {
            cc=qs->Adiv_SI_inv_table[jj][i];
            cc*=Bi_mul[j]; h32=cc;
            MMREDUCE; cc=h32; if (cc>=p) cc-=p;
            qs->Adiv_SI_add[j][i]=(uint16_t)cc;
            cc0+=cc; if (cc0>=p) cc0-=p;
          } else j0=j;
        }
        if (cc0) cc0=p-cc0;
/* using relation: 2*start_value=-sum of change_values */
        qs->Adiv_SI_add[j0][i]=cc0;
      }
/* divisors of k */
    for (i=0; i<qs->nFBk; i++) {
      p=qs->FBk[i]; B_modp=0;
      bb=(uint32_t)(((uint64_t)qs->A)%((uint64_t)p));
      bb=mpqs_invert(bb,p);
      for (j=0; j<qs->nAdiv; j++) {
        cc=(uint32_t)((uint64_t)(qs->Bi[j])%((uint64_t)p));
        B_modp+=((uint32_t)p-cc);
        cc=2*bb*cc; cc%=(uint32_t)p;
        qs->SIk_add[j][i]=(uint16_t)cc;
      }
      c1=B_modp*bb+qs->disp; c1%=(uint32_t)p;
      qs->FBk_start[i]=(uint16_t)c1;
    }
    qs->B_index=0;
    return 1;
  }
  ind=qs->B_index; i=0;
  while (1) {
    if (ind&1) break;
    i++; ind>>=1;
  }
  ind>>=1; ind++;
  add=2*qs->Bi[i];
  if (ind&1) {
    qs->B-=add;
    kmpqs_np_update(fbs,qs->FB_p2,qs->SI_add2[i],2*qs->nFB,0);
    for (j=0; j<qs->nFBk; j++) {
      p=qs->FBk[j];
      aa=qs->FBk_start[j]+qs->SIk_add[i][j];
      if (aa>=p) aa-=p;
      qs->FBk_start[j]=aa;
    }
  } else {
    qs->B+=add;
    kmpqs_np_update(fbs,qs->FB_p2,qs->SI_add2[i],2*qs->nFB,1);
    for (j=0; j<qs->nFBk; j++) {
      p=qs->FBk[j];
      aa=qs->FBk_start[j]+(p-qs->SIk_add[i][j]);
      if (aa>=p) aa-=p;
      qs->FBk_start[j]=aa;
    }
  }
/*zeitB(7);*/
  qs->twoB=2*qs->B;
  qs->C=(long long)(((uint64_t)qs->B*(uint64_t)qs->B-qs->kN_64)*qs->A_inv_64);
//printf("pol2: %lld %lld %lld \n",mpqs_A,mpqs_B,mpqs_C);
  for (j=0; j<qs->nAdiv_total; j++)
    if (!qs->Adiv_active[j]) {
      p=qs->Adiv_all[j];
      s1=qs->Adiv_start1[j]; s2=qs->Adiv_start2[j];
      if (ind&1) {
        s1+=qs->Adiv_SI_add[i][j]; if (s1>=p) s1-=p;
        s2+=qs->Adiv_SI_add[i][j]; if (s2>=p) s2-=p;
      } else {
        s1+=(p-qs->Adiv_SI_add[i][j]); if (s1>=p) s1-=p;
        s2+=(p-qs->Adiv_SI_add[i][j]); if (s2>=p) s2-=p;
      }
      qs->Adiv_start1[j]=s1; qs->Adiv_start2[j]=s2;
    } else {
      p=qs->Adiv_all[j];
      s1=qs->Adiv_start1[j];
      if (ind&1) {
        s1+=qs->Adiv_SI_add[i][j]; if (s1>=p) s1-=p;
      } else {
        s1+=(p-qs->Adiv_SI_add[i][j]); if (s1>=p) s1-=p;
      }
      qs->Adiv_start1[j]=s1;
    }
  return 1;
}


static void mpqs_sieve(mpqs_ctx_t *qs)
{
  uint8_t *sv=qs->sievearray, *fblog=qs->FB_log, *svend;
  uint16_t *fb=qs->FB, *fbs=qs->FB_start;
  uint16_t p, s1, s2, end, pb;
  uint32_t i;
  uint8_t lo;
  uint32_t *ulsv, *ulsvend, m, m1, m2, *ulsv0;
  uint64_t *ullsv, *ullsvend, mask0, mask1, m64[3];

/* initialize odd tiny primes */
  ulsv=(uint32_t *)qs->tinyarray;
  if (qs->ntiny==0) {
    *ulsv=0;
  } else {
    uint32_t j;
    uint32_t h, *ulsv0;

    for (i=0; i<qs->tiny_prod; i++) *ulsv++=0;
    for (i=1; i<qs->sievebegin; i++) {
      p=fb[2*i]; h=qs->tiny_prod/((uint32_t)p);
      lo=fblog[i]; sv=qs->tinyarray;
      s1=fbs[2*i]; s2=fbs[2*i+1];
      for (j=0; j<h; j++) {
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
        sv[s1]+=lo; sv[s2]+=lo; sv+=p;
      }
    }
    for (i=0; i<qs->sievebegink; i++) {
      p=qs->FBk[i]; h=qs->tiny_prod/((uint32_t)p);
      lo=qs->FBk_log[i]; sv=qs->tinyarray;
      s1=qs->FBk_start[i];
      for (j=0; j<h; j++) {
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
        sv[s1]+=lo; sv+=p;
      }
    }
  }
  memcpy(qs->tinyarray+4*qs->tiny_prod,qs->tinyarray,4*qs->tiny_prod);
  memcpy(qs->tinyarray+8*qs->tiny_prod,qs->tinyarray,8*qs->tiny_prod);

/* add contributions of powers of 2 */
  lo=qs->FB_log[0];
  m=((uint32_t)qs->A)*(-1-qs->disp)+(uint32_t)qs->twoB;
  m*=(-1-qs->disp);
  m+=(uint32_t)qs->C;
  m1=((uint32_t)qs->A)*(-3-2*qs->disp)+(uint32_t)qs->twoB;
  m2=2*((uint32_t)qs->A);

  m64[0]=0; m64[1]=0;
  sv=(uint8_t *)(&(m64[0]));
  for (i=0; i<16; i++) {
    m1+=m2; m+=m1;
    if (m&1) continue;
    if (m&6) KMPQS_FAIL("si mod 8\n");
    if (m&8) { sv[i]+=3*lo; continue; }
    if (m&16) { sv[i]+=4*lo; continue; }
    sv[i]+=6*lo;
  }

  mask0=m64[0]; mask1=m64[1];
  ullsv=(uint64_t *)qs->tinyarray;
  ullsvend=ullsv+2*qs->tiny_prod;
  while (ullsv<ullsvend) {
    *ullsv+++=mask0;
    *ullsv+++=mask1;
  }

/* contribution of log |q(x)| */
  kmpqs_sieve_init(qs->sievearray,qs->sievelen,qs->sinit_tab,qs->tinyarray,
                   16*qs->tiny_prod,lo);


  kmpqs_sieve_fb(qs->sievearray,qs->sievelen,qs->FB,qs->FB_start,qs->FB_log,qs->sievebegin);
  sv=qs->sievearray;
  for (i=0; i<qs->nAdiv_total; i++)
    if (!qs->Adiv_active[i]) {
      p=qs->Adiv_all[i]; lo=qs->Adiv_log[i];
      s1=qs->Adiv_start1[i];
      s2=qs->Adiv_start2[i];
      while (s1<qs->sievelen) { sv[s1]+=lo; s1+=p; }
      while (s2<qs->sievelen) { sv[s2]+=lo; s2+=p; }
    } else {
      p=qs->Adiv_all[i]; lo=qs->Adiv_log[i];
      s1=qs->Adiv_start1[i];
      while (s1<qs->sievelen) { sv[s1]+=lo; s1+=p; }
    }
  for (i=qs->sievebegink; i<qs->nFBk; i++)
    kmpqs_sieve_one(qs->sievearray,qs->sievelen,qs->FBk[i],qs->FBk_start[i],qs->FBk_log[i]);
}


static uint16_t mpqs_evaluate(mpqs_ctx_t *qs)
{
  uint16_t **rels, buffer[256];
  uint32_t i, nmaxsurv;

  qs->nsurvivors=0;
  nmaxsurv=MPQS_MAX_NRELS-qs->nrels;
  if (MPQS_MAX_NPRELS-qs->nprels<nmaxsurv)
    nmaxsurv=MPQS_MAX_NPRELS-qs->nprels;
  if (MPQS_MAX_NCRELS-qs->ncrels<nmaxsurv)
    nmaxsurv=MPQS_MAX_NCRELS-qs->ncrels;
  if (nmaxsurv<MPQS_MIN_RELBUFFER) return 0;
  if (nmaxsurv>255) nmaxsurv=255;  /* 1 weniger als buffer-Länge !!! */
  rels=qs->rel_buffer;
#ifdef MPQS_STAT
  stat_asm_eval++;
#endif
  qs->nsurvivors=(uint16_t)kmpqs_eval_scan(qs->sievearray,qs->sievelen,buffer,nmaxsurv);

  for (i=0; i<qs->nsurvivors; i++) {
    qs->sievearray[buffer[i]]=(uint8_t)(i+1);
    rels[i][0]=buffer[i];
  }
  return qs->nsurvivors;
}


/* change hash_table size to 256*15 ??? */
static inline uint16_t mpqs_hash_lp(mpqs_ctx_t *qs, uint32_t lp)
{
  uint32_t j;
  uint32_t lp1;
  uint16_t lphash, nh;

  lphash=(uint16_t)(lp&0x00ff); lphash>>=1;
  lp1=lp>>8;
  nh=qs->hash_table[lphash][0];
  if (nh>=15) { /*printf(".");*/ return MPQS_HASH_OVERFLOW; }
  for (j=1; j<=nh; j++)
    if (qs->hash_table[lphash][j]==lp1) break;
  if (j>nh) { /* new lp */
    qs->hash_table[lphash][j]=lp1;
    qs->hash_table[lphash][0]++;
    qs->nlp++;
    qs->hash_index[lphash][j]=qs->nprels;
    qs->lp[qs->nprels]=lp;
  }
  return qs->hash_index[lphash][j];
}


static inline int mpqs_hash_rel(mpqs_ctx_t *qs, long long axb)
{
  uint32_t j;
  uint32_t hash, hash2, nh;

  hash=(uint32_t)(axb&0xffffff);
  hash2=(uint16_t)(hash>>8); hash&=0xff;
  nh=qs->rel_hash_table[hash][0];
  if (nh>=15) return 1;
  for (j=1; j<=nh; j++)
    if (qs->rel_hash_table[hash][j]==hash2) return 1;
  qs->rel_hash_table[hash][nh+1]=hash2;
  qs->rel_hash_table[hash][0]++;
  return 0;
}


static void mpqs_decompose(mpqs_ctx_t *qs)
{
  uint32_t i, j;
  uint16_t *fb=qs->FB, *fbs=qs->FB_start;
  uint8_t *sv, *svend, i1, i2;
  uint16_t ind, p, s1, s2, nr, d, lpi, r, ii;
  uint16_t **rels;
  long long axb, *llp;
  uint64_t qx, pr, prasm;
  uint16_t d2, ii2, minus;
  double dbl_qx, dx;
  short x;
  uint32_t *ultarg, *ulsrc, lp, ulqx;
  uint32_t inv, h, hi, ls1, ls2;
  uint32_t ax, ay, az, at, asmh;
  uint32_t rr, jj;

  rels=qs->rel_buffer-1;
  for (i=1; i<=qs->nsurvivors; i++) rels[i][4]=0;
  i=qs->td_begin; fb+=2*i; fbs+=2*i;
  for (; i<qs->nFB; i++) {
    p=*fb; fb+=2; sv=qs->sievearray; svend=sv+qs->sievelen-2*p;
    s1=*fbs++; s2=*fbs++;
    while (sv<svend) {
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=qs->nFBk_1+i; rels[i1][4]++;}
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=qs->nFBk_1+i; rels[i2][4]++;}      }
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=qs->nFBk_1+i; rels[i1][4]++;}
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=qs->nFBk_1+i; rels[i2][4]++;}
      }
    }
    svend+=p;
    while (sv<svend) {
      i1=sv[s1]; i2=sv[s2]; sv+=p;
      if (i1||i2) {
        if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=qs->nFBk_1+i; rels[i1][4]++; }
        if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=qs->nFBk_1+i; rels[i2][4]++; }
      }
    }
    svend+=p;
    if (sv+s1<svend) {
      i1=sv[s1];
      if (i1) { nr=rels[i1][4]; rels[i1][nr+5]=qs->nFBk_1+i; rels[i1][4]++; }
    }
    if (sv+s2<svend) {
      i2=sv[s2];
      if (i2) { nr=rels[i2][4]; rels[i2][nr+5]=qs->nFBk_1+i; rels[i2][4]++; }
    }
  }
#ifdef MPQS_STAT
 stat_td_cand+=qs->nsurvivors;
#endif
 for (i=1; i<=qs->nsurvivors; i++) {
   nr=rels[i][4];
   ind=rels[i][0];
   x=(short)(ind)-qs->disp;


   axb=qs->A*(long long)x+qs->B;
   if (axb<0) axb=-axb;
   if (mpqs_hash_rel(qs, axb)) goto next;

   qx=qs->A*(uint64_t)x+qs->twoB;
   qx*=(uint64_t)x;
   qx+=(uint64_t)qs->C;

/* das folgende besser durch ungefähre Berechnung der Nullstellen
   und Vergleiche ersetzen ?? */
   minus=0; dx=(double)x;
   dbl_qx=(double)qs->A*dx+(double)qs->twoB;
   dbl_qx*=dx;
   dx=(double)qs->C; if (dx>0.) dx-=18446744073709551616.;
   dbl_qx+=dx;
   if (fabs(dbl_qx)>=18446744073709551616.) { /*printf(";");*/ goto next; }
   if (dbl_qx<0.) { minus=1; qx=-qx; }
#ifdef MPQS_STAT
   stat_asm_td++;
#endif

    {
      uint32_t ntd=kmpqs_td_roots(ind,qs->FB_start,qs->FB_p2,qs->FB_inv2,
                                  qs->td_begin,qs->td_list);
      for (j=0; j<ntd; j++) rels[i][5+nr+j]=qs->nFBk_1+qs->td_list[j];
      nr+=ntd;
    }

    pr=1;
    for (j=0; j<nr; j++) pr*=(uint64_t)(qs->FB[2*(rels[i][5+j]-qs->nFBk_1)]);
    rels[i][4]=nr;
    if (minus) {
      rels[i][5+rels[i][4]]=0; rels[i][4]++;
    }
    while (!(qx&1)) {
      if (rels[i][4]<MPQS_TD_MAX_NDIV) {
        rels[i][5+rels[i][4]]=qs->nFBk_1; rels[i][4]++;
        qx>>=1;
      } else {
        goto next;
      }
    }
    if (pr<4294967296ULL) {
      if (qx>(pr<<32)) {
/*        printf(",");*/
        goto next;
      }
    }
    ax=(uint32_t)qx; ay=(uint32_t)pr; /*az=0;*/
    inv=(uint32_t)mpqs_256_inv_table[(ay&0xff)>>1];
    at=ay*inv; at&=0x0000ff00;
    at*=inv; inv-=at;
    at=ay*inv; at&=0xffff0000;
    at*=inv; inv-=at;
    az=ax*inv;
    ulqx=az;

    for (j=0; j<nr; j++) {
      uint32_t q;

      ii=rels[i][5+j]; jj=ii-qs->nFBk_1;
      while ((q=ulqx*qs->FB_inv[jj])<=qs->FB_lim[jj]) {  /* p | ulqx */
        if (rels[i][4]<MPQS_TD_MAX_NDIV) {
          rels[i][5+rels[i][4]]=ii; rels[i][4]++;
          ulqx=q;
        } else {
          goto next;
        }
      }
    }


    if (ulqx>1)
      for (j=0; j<qs->nFBk; j++) {
        uint32_t q;

        while ((q=ulqx*qs->FBk_inv[j])<=qs->FBk_lim[j]) {
          if (rels[i][4]<MPQS_TD_MAX_NDIV) {
            rels[i][5+rels[i][4]]=1+j; rels[i][4]++;
            ulqx=q;
          } else {
            goto next;
          }
        }
        if (ulqx==1) break;
      }
    if (ulqx>1)
      for (j=0; j<qs->nAdiv_total; j++) {
        uint32_t q;

        while ((q=ulqx*qs->FB_A_inv[j])<=qs->FB_A_lim[j]) {
          if (rels[i][4]<MPQS_TD_MAX_NDIV) {
            rels[i][5+rels[i][4]]=1+qs->nFB+qs->nFBk+j; rels[i][4]++;
            ulqx=q;
          } else {
            goto next;
          }
        }
        if (ulqx==1) break;
      }



    if (rels[i][4]<MPQS_TD_MAX_NDIV-qs->nAdiv) {
      nr=rels[i][4]; ind=0;
      for (j=0; j<qs->nAdiv_total; j++)
        if (qs->Adiv_active[j]) {
          rels[i][5+nr+ind]=1+qs->nFB+qs->nFBk+j;
          ind++;
        }
      rels[i][4]+=qs->nAdiv;
    } else {
      goto next;
    }
    if (ulqx==1) {
#ifdef MPQS_STAT
      stat_td_surv++;
#endif
      llp=(long long *)(qs->relations[qs->nrels]);
      *llp=axb;
      nr=rels[i][4];
      memcpy(qs->relations[qs->nrels]+4,rels[i]+4,(1+nr)*sizeof(uint16_t));
      qs->nrels++;
      goto next;
    }
    if (ulqx<1048576) { /*printf("%d ",x);*/  /* !!! */
      if (ulqx<=qs->pmax) {
        if (kmpqs_verbose)
          gmp_fprintf(stderr,"kmpqs: error in mpqs_decompose for %Zd: %llu %u\n",
                      qs->N,(unsigned long long)qx,ulqx);
        goto next;
//        KMPQS_FAIL("dec.2\n");
      }
      lp=ulqx;
      lpi=mpqs_hash_lp(qs, lp);
      if (lpi==MPQS_HASH_OVERFLOW) goto next;
#ifdef MPQS_STAT
      stat_td_surv++;
#endif
      if (lpi<qs->nprels) {
        if (qs->lp[lpi]!=lp) KMPQS_FAIL("lp");
        llp=(long long *)(qs->comb_rels[qs->ncrels]);
        *llp=axb;
        llp=(long long *)(qs->part_rels[lpi]);
        axb=*llp;
        llp=(long long *)(qs->comb_rels[qs->ncrels]); llp++;
        *llp=axb;
        nr=rels[i][4];
        memcpy(qs->comb_rels[qs->ncrels]+9,rels[i]+5,nr*sizeof(uint16_t));
        memcpy(qs->comb_rels[qs->ncrels]+9+nr,qs->part_rels[lpi]+5,
               qs->part_rels[lpi][4]*sizeof(uint16_t));
        qs->comb_rels[qs->ncrels][8]=nr+256*qs->part_rels[lpi][4];
        qs->ncrels++;
      } else {
        llp=(long long *)(qs->part_rels[qs->nprels]);
        *llp=axb;
        nr=rels[i][4];
        memcpy(qs->part_rels[qs->nprels]+4,rels[i]+4,
               (1+nr)*sizeof(uint16_t));
        qs->nprels++;
      }
    } /*else stat_counter0++;*/ /*else printf("<");*/
next: ;
  }
  if (qs->nsp<qs->nrels+qs->ncrels)
    qs->excess=qs->nrels+qs->ncrels-qs->nsp;
  else qs->excess=0;
}


/*
order of the primes:
-1, primes in mpqs_FBk, primes in mpqs_FB, primes in mpqs_Adiv_total
*/


static int mpqs_matrix_init(mpqs_ctx_t *qs)
{
  uint32_t i, j, i64, m;
  uint16_t nr, ii, pi;
  uint64_t mask, *M=qs->gauss_row64;

  qs->gauss_n=qs->nrels+qs->ncrels;
  if (qs->gauss_n>=MPQS_GAUSS_MAX) {
    return -1;
  }
  qs->gauss_m=m=qs->nsp;
  if (qs->gauss_m>=qs->gauss_n) KMPQS_FAIL("gauss: no excess\n");
/* build matrix, word-major: bit c of row r in M[r+m*(c>>6)] */
  qs->gauss_n64=(qs->gauss_n+63)/64;
  memset(M,0,m*qs->gauss_n64*sizeof(uint64_t));
  for (i=0; i<qs->nrels; i++) {
    nr=qs->relations[i][4];
    ii=i+qs->ncrels; mask=1ULL<<(ii&63); i64=ii>>6;
    for (j=0; j<nr; j++) {
      pi=qs->relations[i][5+j];
      M[pi+m*i64]^=mask;
    }
  }
  for (i=0; i<qs->ncrels; i++) {
    nr=qs->comb_rels[i][8]; nr=(nr&255)+(nr>>8);
    ii=i; mask=1ULL<<(ii&63); i64=ii>>6;
    for (j=0; j<nr; j++) {
      pi=qs->comb_rels[i][9+j];
      M[pi+m*i64]^=mask;
    }
  }
  qs->gauss_k=qs->gauss_n;
  for (i=0; i<qs->gauss_n; i++) qs->gauss_d[i]=-1;
  return 0;
}



static int mpqs_rowechelon(mpqs_ctx_t *qs)
{
  kmpqs_rowechelon(qs->gauss_row64,qs->gauss_m,qs->gauss_n,qs->gauss_n64,qs->gauss_d);
  qs->gauss_k=0;
  return 0;
}


static int mpqs_matrix(mpqs_ctx_t *qs)
{
  uint32_t i, k64;
  int notzero;
  uint64_t mask;

  for (i=0; i<qs->gauss_n; i++) qs->sol[i]=0;
  for (; qs->gauss_k<qs->gauss_n; qs->gauss_k++) {
    if (qs->gauss_d[qs->gauss_k]>=0) continue;
    notzero=0;
    mask=1ULL<<(qs->gauss_k&63); k64=qs->gauss_k>>6;
    for (i=0; i<qs->gauss_k; i++)
      if (qs->gauss_d[i]!=-1)
        if (qs->gauss_row64[qs->gauss_d[i]+qs->gauss_m*k64] & mask) {
          qs->sol[i]=1;
          notzero=1;
        }

    if (notzero) {
      qs->sol[qs->gauss_k]=1;
      qs->gauss_k++;
      return 1;
    }
  }
  return 0;
}

static uint32_t mpqs_is_power(mpqs_ctx_t *qs, mpz_t op) /* very sloppy */
{
  uint32_t res;

  for (res=2; res<12; res++)
    if (mpz_root(qs->gmp_help,op,(unsigned long)res))
      return res;
  return 0;
}


static int mpqs_final(mpqs_ctx_t *qs)
{
  uint32_t i, j, k;
  uint16_t nr, nr1, nr2;
  uint32_t p, lp, prod1, prod2;
  uint64_t *up;
  int split;
  size_t modthres;

  modthres=mpz_sizeinbase(qs->N,2)*4;
  if (!qs->nfactors) {
    qs->nfactors=1;
    mpz_set(qs->factors[0],qs->N);
    qs->f_status[0]=0;
  }
  mpqs_rowechelon(qs);
  while (mpqs_matrix(qs)) {
#ifdef MPQS_STAT
stat_mpqs_ntrials++;
#endif
    memset(qs->exp,0,qs->nsp*sizeof(short));
    mpz_set_ui(qs->sq1,1); prod1=1;
    mpz_set_ui(qs->sq2,1); prod2=1;

    for (j=0; j<qs->nrels; j++)
      if (qs->sol[qs->ncrels+j]) {
        nr=qs->relations[j][4];
        up=(uint64_t *)(qs->relations[j]); /*printf("Q: %llu  ",*up);*/
        kmpqs_mpz_set_u64(qs->gmp_help,*up);
        if (j&1) {   /* 'random' distribution of full relations */
          for (k=0; k<nr; k++) qs->exp[qs->relations[j][5+k]]--;
          mpz_mul(qs->sq1,qs->sq1,qs->gmp_help);
          if (mpz_sizeinbase(qs->sq1,2)>modthres)
            mpz_mod(qs->sq1,qs->sq1,qs->N);
        } else {
          for (k=0; k<nr; k++) qs->exp[qs->relations[j][5+k]]++;
          mpz_mul(qs->sq2,qs->sq2,qs->gmp_help);
          if (mpz_sizeinbase(qs->sq2,2)>modthres)
            mpz_mod(qs->sq2,qs->sq2,qs->N);
        }
      }

    for (j=0; j<qs->ncrels; j++)
      if (qs->sol[j]) {
        nr=qs->comb_rels[j][8];
        nr1=nr&255; nr2=nr>>8;
        for (k=0; k<nr1; k++) qs->exp[qs->comb_rels[j][9+k]]++;
        for (k=0; k<nr2; k++) qs->exp[qs->comb_rels[j][9+nr1+k]]--;

        up=(uint64_t *)(qs->comb_rels[j]); /*printf("Q: %llu  ",*up);*/
        kmpqs_mpz_set_u64(qs->gmp_help,*up);
        mpz_mul(qs->sq2,qs->sq2,qs->gmp_help);
        if (mpz_sizeinbase(qs->sq2,2)>modthres)
          mpz_mod(qs->sq2,qs->sq2,qs->N);
        up++;
        kmpqs_mpz_set_u64(qs->gmp_help,*up);
        mpz_mul(qs->sq1,qs->sq1,qs->gmp_help);
        if (mpz_sizeinbase(qs->sq1,2)>modthres)
          mpz_mod(qs->sq1,qs->sq1,qs->N);
      }

    for (j=0; j<qs->nsp; j++)
      if (qs->exp[j]&1) KMPQS_FAIL("final.odd\n");
      else qs->exp[j]>>=1;

    for (j=1; j<qs->nsp; j++)
      if (qs->exp[j]) {
        if (j<1+qs->nFBk) {
          p=(uint32_t)(qs->FBk[j-1]);
        } else if (j<1+qs->nFB+qs->nFBk) {
          k=j-qs->nFBk-1;
          p=(uint32_t)(qs->FB[2*k]);
        } else {
          k=j-qs->nFB-qs->nFBk-1;
          p=(uint32_t)(qs->Adiv_all[k]);
        }
        if (qs->exp[j]>0) {
          for (k=0; k<qs->exp[j]; k++) {
            prod1*=p;
            if (prod1&0xfff00000) {  /* p<4096 */
              mpz_mul_ui(qs->sq1,qs->sq1,(unsigned long)prod1);
              if (mpz_sizeinbase(qs->sq1,2)>modthres) {
                mpz_mod(qs->sq1,qs->sq1,qs->N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              }
              prod1=1;
            }
          }
        } else {
          for (k=0; k<-qs->exp[j]; k++) {
            prod2*=p;
            if (prod2&0xfff00000) {  /* p<4096 */
              mpz_mul_ui(qs->sq2,qs->sq2,(unsigned long)prod2);
              if (mpz_sizeinbase(qs->sq2,2)>modthres) {
                mpz_mod(qs->sq2,qs->sq2,qs->N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              }
              prod2=1;
            }
          }
        }
      }
    if (prod1>1) mpz_mul_ui(qs->sq1,qs->sq1,(unsigned long)prod1);
    if (prod2>1) mpz_mul_ui(qs->sq2,qs->sq2,(unsigned long)prod2);

    mpz_add(qs->gmp_help,qs->sq1,qs->sq2);
    mpz_gcd(qs->gmp_help,qs->gmp_help,qs->N);
    if (mpz_cmp_ui(qs->gmp_help,1)==0) continue;
    if (mpz_cmp(qs->gmp_help,qs->N)==0) continue;

    for (j=0; j<qs->nfactors; j++)
      if (qs->f_status[j]==0) {
        mpz_gcd(qs->sq1,qs->gmp_help,qs->factors[j]);
        if (mpz_cmp_ui(qs->sq1,1)==0) continue;
        if (mpz_cmp(qs->sq1,qs->factors[j])==0) continue;
        mpz_set(qs->factors[qs->nfactors],qs->sq1);
        mpz_divexact(qs->factors[j],qs->factors[j],qs->sq1);
        qs->f_status[j]=kmpqs_psp(qs->factors[j]);
        qs->f_status[qs->nfactors]=kmpqs_psp(qs->factors[qs->nfactors]);
        qs->nfactors++;
        break;
      }

    split=1;
    for (j=0; j<qs->nfactors; j++) if (qs->f_status[j]==0) split=0;
    if (split) return 1;
    if (qs->nfactors>=MPQS_MAX_FACTORS)
      KMPQS_FAIL("final: too many factors\n");
  }
  for (j=0; j<qs->nfactors; j++)
    while (qs->f_status[j]==0) {
      if ((k=mpqs_is_power(qs, qs->factors[j]))) {
        if (qs->nfactors+k-1>MPQS_MAX_FACTORS)
          KMPQS_FAIL("final: too many factors\n");
        if (!mpz_root(qs->factors[j],qs->factors[j],(unsigned long)k))
          KMPQS_FAIL("error in power detection\n");
        qs->f_status[j]=kmpqs_psp(qs->factors[j]);
        for (k--; k; k--) {
          mpz_set(qs->factors[qs->nfactors],qs->factors[j]);
          qs->f_status[qs->nfactors]=qs->f_status[j];
          qs->nfactors++;
        }
      } else return 0;
    }
  if (j>=qs->nfactors) return 1;
  return 0;
}



static int mpqs_factor0(mpqs_ctx_t *qs, mpz_t N, size_t max_bits, mpz_t **factors, int retry)
{
  size_t nbits;
  int err;
  uint32_t ev, i;

#ifdef MPQS_STAT
  stat_mpqs_nsieves=0; stat_mpqs_nsurvivors=0; stat_mpqs_ntrials=0; stat_mpqs_ndiv=0;
  if (retry) stat_retry++;
#endif
  nbits=mpz_sizeinbase(N,2);
  if (nbits>96) {
    kmpqs_report("warning: mpqs called with >96 Bit\n");
    return -2;
  }
  mpz_set(qs->N,N);
  mpqs_choose_multiplier(qs);
  mpqs_choose_parameter(qs, retry);
  mpqs_generate_FB(qs);
  if ((err=mpqs_SI_init(qs))<1) {
    if (kmpqs_verbose)
      gmp_fprintf(stderr,"kmpqs: warning: mpqs: error %d in self-initialization %Zd\n",err,N);
    return -3;
  }
  while (1) {
    if (!mpqs_next_pol(qs)) {
      if (retry && kmpqs_verbose)
        gmp_fprintf(stderr,"kmpqs: warning: not enough polynomials in mpqs for %Zd\n",N);
      return -4;
    }
    mpqs_sieve(qs);
#ifdef MPQS_STAT
stat_mpqs_nsieves++;
#endif
    ev=mpqs_evaluate(qs);
#ifdef MPQS_STAT
stat_mpqs_nsurvivors+=ev; // printf("%u %lld  \n",ev,mpqs_C);
#endif
    if (ev) {
      mpqs_decompose(qs); /*test_rels();*/
      if (qs->nsurvivors && (qs->excess>MPQS_MIN_EXCESS)) {
        if (mpqs_matrix_init(qs)) return -6;
        if (mpqs_final(qs)) {
#ifdef MPQS_STAT
// printf("Stat: %u %u %u %u %u %u %u \n",stat_mpqs_nsieves,stat_mpqs_nsurvivors,stat_mpqs_ntrials,stat_mpqs_ndiv,mpqs_nrels,mpqs_nprels,mpqs_ncrels);
stat_ff+=qs->nrels;
stat_pf+=(qs->nprels+qs->ncrels);
stat_comb+=qs->ncrels;
#endif
          for (i=0; i<qs->nfactors; i++)
            if (mpz_sizeinbase(qs->factors[i],2)>max_bits) return 0;
          *factors=qs->factors;
          return qs->nfactors;
        }
      }
    }
    if (qs->nrels>(MPQS_MAX_NRELS-MPQS_MIN_RELBUFFER)) {
      kmpqs_report("warning: too many relations in mpqs\n");
      return -5;
    }
  }
  return -1;
}


void *kmpqs1_ctx_new(void)
{
  mpqs_ctx_t *qs;

  qs=(mpqs_ctx_t *)calloc(1,sizeof(mpqs_ctx_t));
  if (qs==NULL) return NULL;
  if (setjmp(qs->jmpbuf)) return NULL;  /* internal error in init: cannot happen */
  mpqs_init(qs);
  return qs;
}


void kmpqs1_ctx_free(void *ctx)
{
  mpqs_ctx_t *qs=(mpqs_ctx_t *)ctx;
  uint32_t i;

  if (qs==NULL) return;
  mpz_clear(qs->N); mpz_clear(qs->kN); mpz_clear(qs->dummy);
  mpz_clear(qs->gmp_help); mpz_clear(qs->help);
  mpz_clear(qs->sq1); mpz_clear(qs->sq2);
  for (i=0; i<MPQS_MAX_FACTORS; i++) mpz_clear(qs->factors[i]);
  free(qs->sievearray); free(qs->tinyarray);
  free(qs->relations[0]); free(qs->relations);
  free(qs->part_rels[0]); free(qs->part_rels);
  free(qs->comb_rels[0]); free(qs->comb_rels);
  free(qs->rel_buffer[0]); free(qs->rel_buffer);
  free(qs->rel_active);
  free(qs->gauss_row64);
  free(qs);
}


int kmpqs1_run(void *ctx, mpz_t N, size_t max_bits, mpz_t **factors)
{
  mpqs_ctx_t *qs=(mpqs_ctx_t *)ctx;
  int err;

  if (setjmp(qs->jmpbuf)) return KMPQS_ERR_INTERNAL;
  err=mpqs_factor0(qs, N,max_bits,factors,0);
  if (err==-4) err=mpqs_factor0(qs, N,max_bits,factors,1);
  if (err==-4) err=mpqs_factor0(qs, N,max_bits,factors,2);
  return err;
}
