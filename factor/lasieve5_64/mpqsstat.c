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
#include "gmp-aux.h"

u64_t stat_td_cand=0,stat_td_surv=0;
u64_t stat_ff=0,stat_pf=0,stat_comb=0;

u64_t stat_asm_eval=0, stat_asm_td=0;
u32_t stat_asm_div=0, stat_final_mulmod=0;
u32_t stat_counter0=0, stat_counter1=0, stat_retry=0;
u32_t stat_size[12]={0,0,0,0,0,0,0,0,0,0,0,0};

char *input_line=NULL;
size_t input_line_alloc=0;

u32_t n_candidates, n_candidates_alloc=0;
mpz_t *candidates, aux0, aux1, aux2;


void random_prime(mpz_t rop, gmp_randstate_t rs, u32_t nb)
{
  while (1) {
    mpz_urandomb(rop,rs,nb);
    if (mpz_sizeinbase(rop,2)<nb) continue;
    if (mpz_probab_prime_p(rop,10)) break;
  }
}


void random_composite(mpz_t rop, gmp_randstate_t rs, u32_t nb)
{
  u32_t r, n;

  if (nb<=40) {
    while (1) {
      random_prime(aux0,rs,nb/2);
      random_prime(aux1,rs,nb-nb/2+(rand()&1));
      mpz_mul(rop,aux0,aux1);
      if (mpz_sizeinbase(rop,2)==(size_t)nb) break;
    }
    return;
  }
  while (1) {
    n=nb; mpz_set_ui(rop,1);
    while (n>40) {
      r=rand()%(n-40);
      random_prime(aux0,rs,r+20);
      mpz_mul(rop,rop,aux0);
      n=nb-(u32_t)mpz_sizeinbase(rop,2);
    }
    random_prime(aux0,rs,n+(rand()&1));
    mpz_mul(rop,rop,aux0);
    if (mpz_sizeinbase(rop,2)==(size_t)nb) break;
  }
}


int main(int argc, char *argv[])
{
  mpz_t *num, *f;
  u32_t nbit0, nbit1;
  u32_t nb, nnum, nn, i;
  int n;
  gmp_randstate_t rs;
  double cl;

  setbuf(stdout,NULL);
  mpz_ull_init();
//  init_montgomery_multiplication();
//  initzeit(123);
  if (argc<4)
    complain("Usage: mpqsstat nbitmin nbitmax nnumbers\n");
  nbit0=strtoul(argv[1],NULL,10);
  if (nbit0<20) complain("nbitmin<20\n");
  nbit1=strtoul(argv[2],NULL,10);
  nnum=strtoul(argv[3],NULL,10);
  num=(mpz_t *)xmalloc(nnum*sizeof(mpz_t));
  for (i=0; i<nnum; i++) mpz_init(num[i]);
  mpz_init(aux0);
  mpz_init(aux1);
  mpz_init(aux2);

  printf("#");
  for (i=0; i<argc; i++) printf(" %s",argv[i]);
  printf("\n");

  gmp_randinit(rs,GMP_RAND_ALG_DEFAULT,128);
  for (nb=nbit0; nb<=nbit1; nb++) {
    for (i=0; i<nnum; i++) random_composite(num[i],rs,nb);
    nn=0;
    cl=-(double)clock();
    if (nb<=96) {
      for (i=0; i<nnum; i++) {
        n=mpqs_factor(num[i],nb,&f);
        if (n>=1) nn++;
      }
    } else {
      for (i=0; i<nnum; i++) {
        n=mpqs3_factor(num[i],nb,&f);
        if (n>=1) nn++;
      }
    }
    cl+=(double)clock();
    cl/=CLOCKS_PER_SEC;
    cl/=((double)nnum);
    printf("M %u  %.8f\n",nb,cl);
    if (nn<nnum) {
      printf("#mpqs failures for %u bit: %u out of %u\n",nb,nnum-nn,nnum);
      fprintf(stderr,"%u bit: %u failures, %u numbers\n",nb,nnum-nn,nnum);
    }
  }
  exit(0);
}

