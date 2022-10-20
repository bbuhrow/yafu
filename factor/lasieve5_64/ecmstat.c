/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include "asm/siever-config.h"
#include "if.h"

#include "asm/montgomery_mul.h"

char *input_line=NULL;
size_t input_line_alloc=0;

u32_t n_candidates, n_candidates_alloc=0;
mpz_t *candidates;

#define LEVELS    5
mpz_t N[LEVELS];

/* test numbers for timing */
char *Nstr[13]={
"14394188920690575371","65097776422525603088779139339",
"323859163652677805317143139631058956209",
"1404237468360679436977411373853436211577920596387",
"5678098302081870360364729480037798069525196901951417477721",
"25886921661627941599294334476556349009848626550231384466267535602247",
"110473738115768982925210705854910904125868138683739705729281479130276469317583",
"469736367450300187692528692182364707539915644830763665728621573603706346038874705666503",
"2002487885206564031098488248680740725155169377541287493645480003088415373772205622248958419047821",
"9125591701571767882124324233767455967094296334863356572252095877545380507490417396510477260255681740229773",
"37660282491378331507924244662752019105836222235084270321956451125435577876606070892328222568982128127514311765370523",
"161547511231014199097427087354022911848391497194991225279484107820479699022447775908906825862530739466405749362039858987767323",
"719059028897467194366537545109481277788043269900307547668915242986831748660467425380743813044660741722327545692170682557050086721349973"
};


void random_prime(mpz_t rop, gmp_randstate_t rs, u32_t nb)
{
  while (1) {
    mpz_urandomb(rop,rs,nb);
    if (mpz_sizeinbase(rop,2)<nb) continue;
    if (mpz_probab_prime_p(rop,10)) break;
  }
}


void stat_arith(int argc, char *argv[])
{
  double freq=2.2e9;
  ulong aa[LEVELS+1], bb[LEVELS+1], cc[LEVELS+1], dd[LEVELS+1];
  u32_t i, n, smin, smax, mode, nnum;
  double cl;

  if (argc<6) complain("not enough parameters\n");
  smin=strtoul(argv[2],NULL,10);
  smax=strtoul(argv[3],NULL,10);
  mode=strtoul(argv[4],NULL,10);
  nnum=strtoul(argv[5],NULL,10);
  if (argc>6) freq=strtod(argv[6],NULL);
  nnum>>=4; if (!nnum) nnum=1;

  srand(17);
  for (i=smin; i<smax; i++) {
    if (i>=LEVELS) break;
    if (!set_montgomery_multiplication(N[i])) complain("mm\n");
    printf("# %u bit:",32*(i+2));
    for (n=0; n<LEVELS+1; n++) { aa[i]=(ulong)rand(); bb[i]=(ulong)rand(); }
    asm_mulmod(cc,aa,bb); asm_mulmod(aa,aa,cc); asm_mulmod(bb,bb,cc);
    if (mode&1) {
      cl=-(double)clock();
      for (n=0; n<nnum; n++) {
        asm_copy(dd,bb); asm_copy(dd,cc);
        asm_copy(dd,bb); asm_copy(dd,cc);
        asm_copy(dd,bb); asm_copy(dd,cc);
        asm_copy(dd,bb); asm_copy(dd,cc);
        asm_copy(dd,bb); asm_copy(dd,cc);
        asm_copy(dd,bb); asm_copy(dd,cc);
        asm_copy(dd,bb); asm_copy(dd,cc);
        asm_copy(dd,bb); asm_copy(dd,cc);
      }
      cl+=(double)clock();
      cl/=(double)(16*n); cl*=freq;
      cl/=CLOCKS_PER_SEC;
      printf("  copy %.8f",cl);
    }

    for (n=0; n<LEVELS+1; n++) { aa[i]=(ulong)rand(); bb[i]=(ulong)rand(); }
    asm_mulmod(cc,aa,bb); asm_mulmod(aa,aa,cc); asm_mulmod(bb,bb,cc);
    if (mode&2) {
      cl=-(double)clock();
      for (n=0; n<nnum; n++) {
        asm_add2(aa,bb); asm_add2(bb,cc);
        asm_add2(aa,bb); asm_add2(bb,cc);
        asm_add2(aa,bb); asm_add2(bb,cc);
        asm_add2(aa,bb); asm_add2(bb,cc);
        asm_add2(aa,bb); asm_add2(bb,cc);
        asm_add2(aa,bb); asm_add2(bb,cc);
        asm_add2(aa,bb); asm_add2(bb,cc);
        asm_add2(aa,bb); asm_add2(bb,cc);
      }
      cl+=(double)clock();
      cl/=(double)(16*n); cl*=freq;
      cl/=CLOCKS_PER_SEC;
      printf("  add %.8f",cl);
    }

    for (n=0; n<LEVELS+1; n++) { aa[i]=(ulong)rand(); bb[i]=(ulong)rand(); }
    asm_mulmod(cc,aa,bb); asm_mulmod(aa,aa,cc); asm_mulmod(bb,bb,cc);
    if (mode&4) {
      cl=-(double)clock();
      for (n=0; n<nnum; n++) {
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
        asm_sub(aa,bb,cc); asm_sub(bb,bb,cc);
      }
      cl+=(double)clock();
      cl/=(double)(16*n); cl*=freq;
      cl/=CLOCKS_PER_SEC;
      printf("  sub %.8f",cl);
    }

    for (n=0; n<LEVELS+1; n++) { aa[i]=(ulong)rand(); bb[i]=(ulong)rand(); }
    asm_mulmod(cc,aa,bb); asm_mulmod(aa,aa,cc); asm_mulmod(bb,bb,cc);
    if (mode&8) {
      cl=-(double)clock();
      for (n=0; n<nnum; n++) {
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
        asm_mulmod(aa,bb,cc); asm_mulmod(bb,bb,cc);
      }
      cl+=(double)clock();
      cl/=(double)(16*n); cl*=freq;
      cl/=CLOCKS_PER_SEC;
      printf("  mul %.8f",cl);
    }

    printf("\n");
    if (mp_bits_per_limb==64) i++;
  }
  printf("\n");
  exit(0);
}


int main(int argc, char *argv[])
{
  mpz_t *num, *f;
  u32_t B1, B2;
  u32_t nbit, nbitmin, nbitmax, nnum, nn, i, n;
  i32_t b1;
  gmp_randstate_t rs;
  double prob, cl;

  setbuf(stdout,NULL);
  initzeit(13);
  for (i=0; i<LEVELS; i++) mpz_init_set_str(N[i],Nstr[i],10);
  if (argc<5) {
    printf("Usage: ecmstat nbitmin nbitmax nnumbers B1 [ B2 (def=60B1) ]\n");
    complain("  or:  ecmstat 0 sizemin sizemax mode nnumbers [ freq ]\n");
  }
  nbitmin=strtoul(argv[1],NULL,10);
  if (!nbitmin) stat_arith(argc,argv);

  nbitmax=strtoul(argv[2],NULL,10);
  nnum=strtoul(argv[3],NULL,10);
  b1=strtol(argv[4],NULL,10);
  if (b1<0) B1=(u32_t)(-b1); else B1=(u32_t)b1;
  if (argc>5) B2=strtoul(argv[5],NULL,10); else B2=60*B1;
  nn=nnum; if (nn<100) nn=100;

  num=(mpz_t *)xmalloc(nnum*sizeof(mpz_t));
  for (i=0; i<nnum; i++) mpz_init(num[i]);
  printf("#");
  for (i=0; i<argc; i++) printf(" %s",argv[i]);
  printf("\n");
  gmp_randinit(rs,GMP_RAND_ALG_DEFAULT,128);

  for (nbit=nbitmin; nbit<=nbitmax; nbit++) {
/* generate random numbers of nbit bit */
    for (i=0; i<nnum; i++) random_prime(num[i],rs,nbit);

    n=0;
    for (i=0; i<nnum; i++)
      if (ecm_factor(num[i],B1,B2,&f,1)>0) n++;
    prob=((double)n)/((double)nnum);
    printf("E %u  %u %u %.6f",nbit,B1,B2,prob);
    if ((nbit==nbitmin) && (b1>0)) {
      printf("  ");
      for (i=0; i<LEVELS; i++) {
        cl=-(double)clock();
        n=ecm_factor(N[i],B1,B2,&f,nn);
        cl+=(double)clock();
        if (n>0) cl/=((double)n);
        else if (n==0) cl/=((double)nn);
        else complain("ecm failed for %u bit test number\n",32*(2+i));
        cl/=CLOCKS_PER_SEC;
        printf(" %u:%.8f",64+32*i,cl);
        if (mp_bits_per_limb==64) i++;
      }
    }
    printf("\n");
  }
  exit(0);
}

