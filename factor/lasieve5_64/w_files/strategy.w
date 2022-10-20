@* Doing the cofactorisation.

@*3 Copying.
Copyright (C) 2001,2002 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@(strategy.h@>=
#include <stdarg.h>
#include <stdio.h>

typedef struct {
  u32_t **stindex; /* nr of strategy to use for (n0,n1) */
  u32_t **stlist;  /* list of strategies */
  u16_t mc[2];     /* maximal number of bits for composites */
  unsigned char **bit;
} strat_t;

void print_strategy(strat_t s);
void read_strategy(strat_t *s, u16_t *maxcomp, char *basename, u16_t *maxpr);
int cofactorisation(strat_t *st, mpz_t **large_primes, mpz_t *large_factors,
                      u16_t *max_primebits, u32_t *nlp, mpz_t *FBb_sq,
                      mpz_t *FBb_cu);
void print_strategy_stat();


@
@c
#include <stdio.h>
#include <sys/types.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "asm/siever-config.h"
#include "if.h"
#include "primgen32.h"
#include "asm/32bit.h"
#include "strategy.h"
#include "asm/montgomery_mul.h"
#include "ecm.h"
#include "pm1.h"

#include "mpqs.h"
#include "mpqs3.h"

extern char *input_line;
extern size_t input_line_alloc;

static mpz_t ecm_f[2], ecm_aux;
static u32_t cf_maxcomp[2];

static u32_t nfecm=0, nsecm=0, nfpm1=0, nspm1=0;
static i64_t mpqsaux_clock;

#define CF_STAT
//#define CF_STAT_EXACT

#ifdef CF_STAT
static double **stat_cost;
static double **stat_yield;
static double cost, yield;
static u64_t cf_n=0, cf_necm=0, cf_naux=0, cf_nauxmpqs=0, cf_nauxmpqs3=0;
static u64_t cf_nauxmpqstoobig=0, cf_nauxecm=0;
#endif

#ifdef CF_STAT_EXACT
static u32_t **stat_cand;
static u32_t **stat_success;
static u32_t *stat_mpqsaux, *stat_aux;
#endif


@
@c
void print_strategy(strat_t s)
{
  u32_t i, j, k, l, t;

  printf("Strategy:\n");
  for (i=0; i<=s.mc[0]; i++)
    for (j=0; j<=s.mc[1]; j++)
      if (t=s.stindex[i][j]) {
        printf("%u %u:",i,j);
        for (k=0; l=s.stlist[t][k]; k++) {
          if (l&1) printf(" Y"); else printf(" X");
          if ((l>>2)==1) printf("M");
          else {
            if ((l>>2)&1) printf("P%u:%u",(l>>3)&0x1fff,4*(l>>16));
            else printf("E%u:%u",(l>>3)&0x7ff,l>>14);
          }
          if (l&2) printf("|");
        }
        printf("\n");
      }
  printf("\n");
  printf("bit:\n");
  for (i=0; i<=s.mc[0]; i++) {
    printf("%3u: ",i);
    for (j=0; j<=s.mc[1]; j++) {
      if (s.bit[i][j]) printf("X"); else printf(".");
    }
    printf("\n");
  }
  printf("\n");
}


@ P-1: We use 13 bit for B1 and 16 for B2/4.
@c
static u32_t get_pm1_type(u32_t B1, u32_t B2)
{
  if (B1<2) return 0;
  if ((B1>>13)) return 0;
  if (B2<2*B1) return 0;
  if ((B2>>18)) return 0;
  return 1+2*B1+((B2>>2)<<14);
}


@
@c
static void get_pm1_param(u32_t *b1ptr, u32_t *b2ptr, u32_t type)
{
  *b1ptr=(type>>3)&0x1fff;
  *b2ptr=4*(type>>16);
}


@ ECM: We use 11 bit for B1/2 and 18 for B2/4.
@c
static u32_t get_ecm_type(u32_t B1, u32_t B2)
{
  if (B1<2) return 0;
  if ((B1>>12)) return 0;
  if (B2<2*B1) return 0;
  if ((B2>>20)) return 0;
  return 2*((B1>>1)<<1)+((B2>>2)<<12);
}


@
@c
static void get_ecm_param(u32_t *b1ptr, u32_t *b2ptr, u32_t type)
{
  *b1ptr=(type>>2)&0xffe;
  *b2ptr=(type>>14)<<2;
}


@
@c
u32_t get_fm_type(char **lptr) /* only mpqs, ecm and pm1 at the moment */
{
  char *line, *tail;
  u32_t t, B1, B2;

  line=*lptr;
  if (*line=='M') { /* mpqs */
    t=1; line++;
  } else if (*line=='E') { /* ecm, 11 bit for B1, 18 for B2 */
    line++;
    B1=(u32_t)strtoul(line,&tail,10);
    if ((tail==NULL) || (line==tail))
      complain("st file contains corrupt line:\n%s",line-1);
    line=tail+1;
    B2=(u32_t)strtoul(line,&tail,10);
    if ((tail==NULL) || (line==tail))
      complain("st file contains corrupt line:\n%s",line-1);
    line=tail;
    t=get_ecm_type(B1,B2);
    if (!t) complain("bad ecm parameter in line:\n%s\n",*lptr);
  } else if (*line=='P') {
    line++;
    B1=(u32_t)strtoul(line,&tail,10);
    if ((tail==NULL) || (line==tail))
      complain("st file contains corrupt line:\n%s",line-1);
    line=tail+1;
    B2=(u32_t)strtoul(line,&tail,10);
    if ((tail==NULL) || (line==tail))
      complain("st file contains corrupt line:\n%s",line-1);
    line=tail;
    t=get_pm1_type(B1,B2);
    if (!t) complain("bad pm1 parameter in line:\n%s\n",*lptr);
  } else complain("unknown factorisation type\n");
  *lptr=line;
  return t;
}


@
@c
void read_strategy(strat_t *s, u16_t *maxcomp, char *basename, u16_t *maxpr)
{
  char *ifn;
  FILE *ifile;
  char *line, *tail;
  u32_t n0, n1, len;
  u32_t i, j ,k ,t, tt, b;
  size_t alloc=0;
  u32_t *s_tmp;

  cf_maxcomp[0]=(u32_t)maxcomp[0];
  cf_maxcomp[1]=(u32_t)maxcomp[1];
  mpz_init(ecm_f[0]);
  mpz_init(ecm_f[1]);
  mpz_init(ecm_aux);
  init_montgomery_multiplication();
  memcpy(s->mc,maxcomp,2*sizeof(u16_t));
  s->stindex=(u32_t **)xmalloc((1+maxcomp[0])*sizeof(u32_t *));
  for (i=0; i<=maxcomp[0]; i++)
    s->stindex[i]=(u32_t *)xcalloc((size_t)(1+maxcomp[1]),sizeof(u32_t));
  s->stindex[0][0]=1;
  adjust_bufsize((void**)&(s->stlist),&alloc,1,16,sizeof(u32_t *));
  s->stlist[0]=(u32_t *)xcalloc(1,sizeof(u32_t));
  s->stlist[0][0]=0; /* dummy strategy: do nothing */
  len=1;
#ifdef CF_STAT
  stat_cost=(double **)xmalloc((1+maxcomp[0])*sizeof(double *));
  for (i=0; i<=maxcomp[0]; i++)
    stat_cost[i]=(double *)xcalloc((size_t)(1+maxcomp[1]),sizeof(double));
  stat_yield=(double **)xmalloc((1+maxcomp[0])*sizeof(double *));
  for (i=0; i<=maxcomp[0]; i++)
    stat_yield[i]=(double *)xcalloc((size_t)(1+maxcomp[1]),sizeof(double));
  cost=0.; yield=0.;
#endif
#ifdef CF_STAT_EXACT
  stat_cand=(u32_t **)xmalloc((1+maxcomp[0])*sizeof(u32_t *));
  for (i=0; i<=maxcomp[0]; i++)
    stat_cand[i]=(u32_t *)xcalloc((size_t)(1+maxcomp[1]),sizeof(u32_t));
  stat_success=(u32_t **)xmalloc((1+maxcomp[0])*sizeof(u32_t *));
  for (i=0; i<=maxcomp[0]; i++)
    stat_success[i]=(u32_t *)xcalloc((size_t)(1+maxcomp[1]),sizeof(u32_t));
  i=(maxcomp[0]<maxcomp[1] ? maxcomp[1] : maxcomp[0]);
  stat_aux=(u32_t *)xcalloc((size_t)(1+i),sizeof(u32_t));
  stat_mpqsaux=(u32_t *)xcalloc((size_t)(1+i),sizeof(u32_t));
#endif

  asprintf(&ifn,"%s.st",basename);
  if((ifile=fopen(ifn,"r"))!=0) {
    while (1) {
      if (skip_blanks_comments(&input_line,&input_line_alloc,ifile)==0) break;
      if(input_line == NULL) break;
      line=input_line;
      n0=(u32_t)strtol(line,&tail,10);
      if ((tail==NULL) || (line==tail))
        complain("file %s contains corrupt line:\n%s",ifn,input_line);
      line=tail;
      n1=(u32_t)strtol(line,&tail,10);
      if ((tail==NULL) || (line==tail))
        complain("file %s contains corrupt line:\n%s",ifn,input_line);
      line=tail;
      if (n0>maxcomp[0]) continue;
      if (n1>maxcomp[1]) continue;
#ifdef CF_STAT
      line++;
      stat_cost[n0][n1]=strtod(line,&tail);
      if ((tail==NULL) || (line==tail))
        complain("file %s contains corrupt line:\n%s",ifn,input_line);
      line=tail+1;
      stat_yield[n0][n1]=strtod(line,&tail);
      if ((tail==NULL) || (line==tail))
        complain("file %s contains corrupt line:\n%s",ifn,input_line);
      line=tail;
#endif
      while (*line) {
        if (*line=='X') break;
        if (*line=='Y') break;
        line++;
      }
      s->stindex[n0][n1]=len;
      for (tail=line,i=1; *tail; tail++) if (*tail==',') i++;
      adjust_bufsize((void**)&(s->stlist),&alloc,(size_t)(len+1),
                     16,sizeof(u32_t *));
      s_tmp=(u32_t *)xcalloc((size_t)(i+1),sizeof(u32_t));
      for (j=0; j<i; j++) {
        t=0;
        if (*line=='Y') t++; else if (*line!='X') complain("no XY\n");
        line++;
        tt=get_fm_type(&line);
        t+=4*tt;
        s_tmp[j]=t;
        if (j+1<i) if (*line!=',') complain("no ,  %u %u\n",i,j);
        line++;
      }
      if (*line) complain("st-file corrupt\n");
      s_tmp[j++]=0;
/* mark last factorisation attempt for each side */
      for (i=j-1; i; i--) if ((s_tmp[i-1]&1)==0) { s_tmp[i-1]^=2; break; }
      for (i=j-1; i; i--) if ((s_tmp[i-1]&1)==1) { s_tmp[i-1]^=2; break; }
/* check whether strategy is already stored */
      for (i=0; i<len; i++) {
        for (k=0; k<j; k++) if (s_tmp[k]!=s->stlist[i][k]) break;
        if (k==j) break;
      }
      if (i<len) s->stindex[n0][n1]=i;
      else {
        s->stlist[len]=(u32_t *)xmalloc(j*sizeof(u32_t));
        memcpy(s->stlist[len],s_tmp,j*sizeof(u32_t));
        len++;
      }
      free(s_tmp);
    }
    fclose(ifile);
  } else {
/* create default strategy: bad for >2LP */
    if ((3*maxpr[0]+3<maxcomp[0]) || (3*maxpr[1]+3<maxcomp[1])) {
      logbook(0,"Warning: >2LP but no .st file\n");
    }
/* add two strategies doing only mpqs: */
    adjust_bufsize((void**)&(s->stlist),&alloc,(size_t)(len+2),
                   16,sizeof(u32_t *));
    s->stlist[len]=(u32_t *)xmalloc(3*sizeof(u32_t));
    s->stlist[len][0]=4+2; s->stlist[len][1]=4+2+1; s->stlist[len][2]=0;
    len++;
    s->stlist[len]=(u32_t *)xmalloc(3*sizeof(u32_t));
    s->stlist[len][0]=4+2+1; s->stlist[len][1]=4+2; s->stlist[len][2]=0;
    len++;
    for (n0=0; n0<=maxcomp[0]; n0++) {
      for (n1=0; n1<=maxcomp[1]; n1++) {
        if (n0<n1) s->stindex[n0][n1]=len-1;
        else s->stindex[n0][n1]=len-2;
/* first do mpqs for bigger number then for smaller number, ok for
   <=2LP, probably bad for >2LP */
      }
    }
  }


  s->bit=(unsigned char **)xmalloc((maxcomp[0]+1)*sizeof(unsigned char *));
  for (i=0; i<=maxcomp[0]; i++)
    s->bit[i]=(unsigned char *)xcalloc((size_t)(maxcomp[1]+1),sizeof(unsigned char));
{
  u32_t i0, i1, j0, j1;

  for (i=0; i<=maxcomp[0]; i++) {
    if (i>3) i0=i-3; else i0=0;
    if (i==0) i1=maxpr[0]; else i1=i;
    i1+=1; if (i1>=maxcomp[0]) i1=maxcomp[0];
    for (j=0; j<=maxcomp[1]; j++) {
      if (j>3) j0=j-2; else j0=0;
      if (j==0) j1=maxpr[1]; else j1=j;
      j1+=0; if (j1>=maxcomp[1]) j1=maxcomp[1];
      if (s->stindex[i][j])
        for (k=i0; k<=i1; k++)
          for (t=j0; t<=j1; t++)
            s->bit[k][t]=1;
    }
  }
}
}


@
@c
int cofactorisation(strat_t *st, mpz_t **large_primes, mpz_t *large_factors,
                      u16_t *max_primebits, u32_t *nlp, mpz_t *FBb_sq,
                      mpz_t *FBb_cu)
{
  u32_t s, nb[2];
  u32_t m, *fm, t, j, done[2], B1, B2, pm1done[2];
  static ecm_t e[2];
  static int cf_ecm_init=0;
  clock_t cl;

  for(s=0;s<2;s++) {
    if(mpz_sgn(large_factors[s])>0) {
      if(mpz_cmp_ui(large_factors[s],1)==0)
        nlp[s]=0;
      else {
        nlp[s]=1;
        mpz_set(large_primes[s][0],large_factors[s]);
      }
      nb[s]=0;
    } else {
      mpz_neg(large_factors[s],large_factors[s]);
      nb[s]=(u32_t)(mpz_sizeinbase(large_factors[s],2));
      nlp[s]=2;
    }
  }
  if ((nlp[0]<2) && (nlp[1]<2)) {
#ifdef CF_STAT
    yield+=1.;
#endif
#ifdef CF_STAT_EXACT
    stat_cand[0][0]++;
    stat_success[0][0]++;
#endif
    return 0;
  }

#ifdef CF_STAT
  yield+=stat_yield[nb[0]][nb[1]];
  cost+=stat_cost[nb[0]][nb[1]];
#endif
#ifdef CF_STAT_EXACT
  stat_cand[nb[0]][nb[1]]++;
#endif

  if (cf_ecm_init==0) {
    ecm_curve_init(e[0]);
    ecm_curve_init(e[1]);
    cf_ecm_init=1;
  }
  m=st->stindex[nb[0]][nb[1]];
  if (!m) return 1;
  fm=st->stlist[m];
#ifdef CF_STAT
  cf_n++;
#endif
  for (s=0; s<2; s++)
    if (nlp[s]==2) { nlp[s]=0; done[s]=0; pm1done[s]=0; } else done[s]=2;
  for (j=0; t=fm[j]; j++) {
    i32_t nf, i;
    size_t sf[2];
    mpz_t *fac;

    s=t&1;
    if (done[s]>1) continue;
    if ((t>>2)==1) { /* mpqs */
      if(mpz_sizeinbase(large_factors[s],2)>96)
        nf=mpqs3_factor(large_factors[s],max_primebits[s],&fac);
      else
        nf=mpqs_factor(large_factors[s],max_primebits[s],&fac);
      if(nf<0) return -2;
      if(!nf) return 1;
      for(i=0;i<nf;i++)
        mpz_set(large_primes[s][nlp[s]+i],fac[i]);
      nlp[s]+=nf;
      done[s]=2;
    } else { /* ecm or p-1 */
      if ((t>>2)&1) { /* p-1 */
        get_pm1_param(&B1,&B2,t);
        nf=pm1_factor(large_factors[s],B1,B2,&fac);
        pm1done[s]=1;
nfpm1++; if (nf>0) nspm1++;
      } else { /* ecm */
        get_ecm_param(&B1,&B2,t);
        if (!done[s]) {
          if (ecm_curve_set(e[s],large_factors[s],B1,B2)) return -3;
          done[s]=1;
        } else {
          ecm_set_params(e[s],B1,B2);
        }
        nf=ecm(e[s],&fac);
nfecm++; if (nf>0) nsecm++;
      }
      if (nf<0) return -3;
      if (nf) {
        @<ECM factor found@>@;
      }
    }
    if (t&2) { /* last test for this side */
      if (done[s]<2) return 1;
    }
  }
  if ((done[0]!=2) || (done[1]!=2)) return -1;

#ifdef CF_STAT_EXACT
  stat_success[nb[0]][nb[1]]++;
#endif
  return 0;
}

@ If a factor was found by ECM or p-1, we try to decompose the factor and
the cofactor completely. First some size checks and psp tests are done to
avoid unnecessary computations. If both (factor and cofactor) are composite
we first consider the smaller number.
TODO: If there is a strategy for the bitsize of the number we use it possibly
avoiding a second p-1 test, otherwise an ad hoc strategy is used.

@<ECM factor found@>=
{
  int es, need_test[2], order[2], o;

#ifdef CF_STAT
  cf_necm++;
#endif
  mpz_set(ecm_f[0],fac[0]);
  sf[0]=mpz_sizeinbase(ecm_f[0],2);
  if (sf[0]<=1) return -4;
  if (sf[0]<nb[s]) {
    mpz_fdiv_qr(ecm_f[1],ecm_aux,large_factors[s],ecm_f[0]);
    if (mpz_sgn(ecm_aux)) Schlendrian("ecm found non-divisor\n");
    sf[1]=mpz_sizeinbase(ecm_f[1],2);
/* check of need for psp-test */
    for (es=0; es<2; es++) {
      need_test[es]=0;
      if (sf[es]<=max_primebits[s]) continue;
      if (mpz_cmp(ecm_f[es],FBb_sq[s])<0) return 1; /* prime, too big */
      if (sf[es]<=2*max_primebits[s])  { need_test[es]=1; continue; }
      if (mpz_cmp(ecm_f[es],FBb_cu[s])<0) return 1; /* not smooth */
      need_test[es]=1;
    }
/* psp tests */
    for (es=0; es<2; es++)
      if (need_test[es])
// SMJS psp only has one arg        if (psp(ecm_f[es],1)==1) return 1;
        if (psp(ecm_f[es])==1) return 1;
/* p-1, ecm and mpqs */
    order[0]=0; order[1]=1;
/* if first number is bigger change order of factorisation */
    if ((need_test[0]&need_test[1]) && (sf[0]>sf[1])) {
      order[0]=1; order[1]=0;
    }
    for (o=0; o<2; o++) {
      es=order[o];
      if (!need_test[es]) {
        mpz_set(large_primes[s][nlp[s]],ecm_f[es]); nlp[s]++;
      } else {
#ifdef CF_STAT_EXACT
        stat_aux[mpz_sizeinbase(ecm_f[es],2)]++;
#endif
#ifdef CF_STAT
  cf_naux++;
#endif
        @<Auxiliary factorisation@>;
      }
    }
    done[s]=2;
  }
/* CAVE: shall we try mpqs if all factors were found simultaneously by ecm or p-1? */
}

@
@<Auxiliary factorisation@>=
{
  u32_t ne;
  size_t sz;

cl=clock();
  ne=0;
  while (1) {
    sz=mpz_sizeinbase(ecm_f[es],2);
    if(sz<=96) {
#ifdef CF_STAT_EXACT
        stat_mpqsaux[sz]++;
#endif
#ifdef CF_STAT
  cf_nauxmpqs++;
#endif
      nf=mpqs_factor(ecm_f[es],max_primebits[s],&fac);
    } else {
      if (5*ne>(sz-90)) {
#ifdef CF_STAT
  cf_nauxmpqs3++;
#endif
#ifdef CF_STAT_EXACT
        stat_mpqsaux[sz]++;
#endif
        if (sz>128) {
#ifdef CF_STAT
  cf_nauxmpqstoobig++;
#endif
          mpqsaux_clock+=clock()-cl; return -1;
        }
        nf=mpqs3_factor(ecm_f[es],max_primebits[s],&fac);
      } else { /* ecm */
        B1=500; B2=36000;  /* TODO: more intelligent choice of parameters */
        if (!done[s]) {
          if (ecm_curve_set(e[s],ecm_f[es],B1,B2)) {
            mpqsaux_clock+=clock()-cl;
            return -3;
          }
          done[s]=1;
        } else {
          if (!ne) {
            if (ecm_reset_n(e[s],ecm_f[es])) {
              mpqsaux_clock+=clock()-cl;
              return -3;
            }
          }
          ecm_set_params(e[s],B1,B2);
        }
#ifdef CF_STAT
  cf_nauxecm++;
#endif
        nf=ecm(e[s],&fac);
nfecm++; if (nf>0) nsecm++;
        ne++;
        if (nf<0) { mpqsaux_clock+=clock()-cl; return -3; }
        if (nf==0) continue;
/* TODO: at the moment we ignore the case that ecm finds a composite factor */
        if (mpz_sizeinbase(fac[0],2)<=max_primebits[s]) {
          mpz_set(large_primes[s][nlp[s]],fac[0]); nlp[s]++;
          mpz_fdiv_qr(ecm_f[es],ecm_aux,ecm_f[es],fac[0]);
          if (mpz_sgn(ecm_aux)) {
            gmp_printf("%Zd %Zd\n%Zd\n",ecm_aux,fac[0],large_factors[s]);
            Schlendrian("ecm found non-divisor\n");
          }
          if (mpz_sizeinbase(ecm_f[es],2)<=max_primebits[s]) {
            mpz_set(large_primes[s][nlp[s]],ecm_f[es]); nlp[s]++;
            break;
// SMJS psp only has one arg          } else if (psp(ecm_f[es],1)==1) {
          } else if (psp(ecm_f[es])==1) {
            mpqsaux_clock+=clock()-cl; return 1;
          }
          if (ecm_reset_n(e[s],ecm_f[es])) {
            mpqsaux_clock+=clock()-cl;
            return -3;
          }
          continue;
        }
        // SMJS psp only has one arg if (psp(fac[0],1)==1) { mpqsaux_clock+=clock()-cl; return 1; }
        if (psp(fac[0])==1) { mpqsaux_clock+=clock()-cl; return 1; }
        continue;
      }
    }
    if (nf>0) {
      for(i=0;i<nf;i++)
        mpz_set(large_primes[s][nlp[s]+i],fac[i]);
      nlp[s]+=nf;
      break;
    }
mpqsaux_clock+=clock()-cl;
    if(nf<0) return -4;
    if(!nf) return 1;
  }
mpqsaux_clock+=clock()-cl;
}

@
@c
void print_strategy_stat()
{
#ifdef CF_STAT
  logbook(0,"Expected yield/cost: %.4g  %.4g\n",yield,cost);
  logbook(0,"p-1: %u tests, %u successes  ecm: %u tests, %u successes\n",
  nfpm1,nspm1,nfecm,nsecm);
  mpqsaux_clock=rint((1000.0*mpqsaux_clock)/CLOCKS_PER_SEC);
  logbook(0,"MPQS-AUX %d\n",(int)mpqsaux_clock);
  logbook(0,"COF: %llu tests, %llu ecm, %llu aux:\n",cf_n,cf_necm,cf_naux);
  logbook(0,"       %llu mpqs, %llu mpqs3, %llu ecm, %llu too big\n",
          cf_nauxmpqs,cf_nauxmpqs3,cf_nauxecm,cf_nauxmpqstoobig);
#endif

#ifdef CF_STAT_EXACT
{
  u32_t i, j;

  logbook(0,"\nExact stat:\n");
  for (i=0; i<=cf_maxcomp[0]; i++)
    for (j=0; j<=cf_maxcomp[1]; j++)
      if (stat_cand[i][j]) {
        logbook(0,"%u %u: %u -> %u (%g)\n",i,j,stat_cand[i][j],
               stat_success[i][j],((double)stat_cand[i][j])*stat_yield[i][j]);
      }
  logbook(0,"\nmpqs-aux:\n");
  j=(cf_maxcomp[0]<cf_maxcomp[1] ? cf_maxcomp[1] : cf_maxcomp[0]);
  for (i=0; i<=j; i++)
    if (stat_aux[i])
      logbook(0,"%u: %u (%u)\n",i,stat_aux[i],stat_mpqsaux[i]);
}
#endif
}
