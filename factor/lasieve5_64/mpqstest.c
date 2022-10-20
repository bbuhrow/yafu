/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/


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

int iter=0;
u64_t stat_asm_eval=0, stat_asm_td=0;
u64_t stat_td_cand=0,stat_td_surv=0;
u64_t stat_ff=0,stat_pf=0,stat_comb=0;
u32_t stat_asm_div=0, stat_final_mulmod=0;
u32_t stat_counter0=0, stat_counter1=0, stat_retry=0;
u32_t stat_size[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};


int main(int argc, char **argv)
{
  mpz_t N, *f;
  int n, i;

  setbuf(stdout,NULL);
  initzeit(55); zeita(0);
  mpz_ull_init();
//  init_montgomery_multiplication();
  mpz_init(N);

#if 1
{
  FILE *fi=NULL;
  int suc=0, fail=0;
  char *input_line=NULL;
  size_t input_line_alloc=0;

  if (argc>1) {
    if(strcmp(argv[1],"-")) fi=fopen(argv[1],"r");
    else fi=stdin;
  } else fi=fopen("comp.96","r");
  if (fi==NULL) complain("cannot open file\n");
  while (getline(&input_line,&input_line_alloc,fi)>0) { iter++;
    if (mpz_set_str(N,input_line,10)) break;
    n=mpqs_factor(N,80,&f);
    if (n>=1) suc++; else fail++;
  }
  if(fi!=stdin) fclose(fi);
  printf("success: %d  failure: %d\n",suc,fail);
}
#else
  mpz_set_str(N,"1208925819614629174706951",10);
  mpz_set_str(N,"1237940039285380274899124231",10);
  mpz_set_str(N,"39614081257132168796771916653",10);
  mpz_set_str(N,"39614081257132168796771916997",10);
  mpz_set_str(N,"1152921504606847177",10);
  mpz_set_str(N,"3963679138604680183237",10);
  mpz_set_str(N,"9094508562016968604909695137",10);
  mpz_set_str(N,"40240008760233165180623235631",10);
  mpz_set_str(N,"183495473821062038788219961",10);
  mpz_set_str(N,"301136548239571127",10);
  mpz_set_str(N,"43131998832981427596515945471",10);
  mpz_set_str(N,"1532124483140137",10);
  mpz_set_str(N,"2439308476736629369",10);
  mpz_set_str(N,"709745647298040929",10);
  for (i=0; i<100; i++)  n=mpqs_factor(N,70,&f);
  for (i=0; i<n; i++) {
    mpz_out_str(stdout,10,f[i]); printf(" ");
  }
  printf("\n");
#endif
  zeitb(0);
  if (stat_retry) printf("Warning: %u retries\n",stat_retry);
  printf("Stat: sieves %Lu, td %Lu->%Lu->%Lu\n",
	 stat_asm_eval,stat_td_cand,stat_asm_td,stat_td_surv);
  printf("  ff %llu, pf %llu->%llu\n",stat_ff,stat_pf,stat_comb);
  for (i=0; i<14; i++) if (stat_size[i])
    printf("size %d: %u  ",i,stat_size[i]); printf("\n");

  printf("timing:  total: "); printzeit(0); printf("\n");
  printf("  init       : "); printzeit(9); printf("\n");
  printf("  init fb    : "); printzeit(1); printf("\n");
  printf("  nextpol    : "); printzeit(2); printf("\n");
  printf("  sieve      : "); printzeit(3); printf("\n");
  printf("  eval       : "); printzeit(4); printf("\n");
  printf("  decompose  : "); printzeit(5); printf("\n");
  printf("  final      :"); printzeit(10); printf("\n");
  printf("    matrix   : "); printzeit(8); printf("\n");
printzeit(20); printzeit(21); printzeit(22); printzeit(23); printzeit(24); printzeit(25); printzeit(26); printf("\n");
}

