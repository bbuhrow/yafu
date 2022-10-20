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
#include <gmp.h>
#include "asm/siever-config.h"
#include "if.h"
#include "gmp-aux.h"

int iter=0;
u64_t stat_asm_eval=0, stat_asm_td=0;
u64_t stat_td_cand=0,stat_td_surv=0;
u64_t stat_ff=0,stat_pf=0,stat_comb=0,stat_pp=0;
u32_t stat_asm_div=0, stat_final_mulmod=0;
u32_t stat_counter0=0, stat_counter1=0, stat_retry=0;
u32_t stat_size[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


int main(int argc, char **argv)
{
  mpz_t N, *f;
  int n, i;

  setbuf(stdout,NULL);
  initzeit(100); zeita(0);
  mpz_init(N);
#ifdef ULL_NO_UL
  mpz_ull_init();
#endif
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
    n=mpqs3_factor(N,200,&f);
    if (n>=1) suc++; else fail++;
  }
  if(fi!=stdin) fclose(fi);
  printf("success: %d  failure: %d\n",suc,fail);
}
#else
  mpz_set_str(N,"709745647298040929",10);
  for (i=0; i<100; i++)  n=mpqs3_factor(N,200,&f);
  for (i=0; i<n; i++) {
    mpz_out_str(stdout,10,f[i]); printf(" ");
  }
  printf("\n");
#endif
  zeitb(0);
  if (stat_retry) printf("Warning: %u retries\n",stat_retry);
  printf("Stat: sieves %Lu, td %Lu->%Lu->%Lu\n",
	 stat_asm_eval,stat_td_cand,stat_asm_td,stat_td_surv);
  printf("  ff %llu, pf %llu->%llu (%llu)\n",stat_ff,stat_pf,stat_comb,stat_pp);
  for (i=0; i<15; i++) if (stat_size[i])
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
//printzeit(20); printzeit(21); printzeit(22); printzeit(23); printzeit(24); printzeit(25); printzeit(26); printf("\n");
//printzeit(30); printzeit(31); printzeit(32); printzeit(33); printzeit(34); printzeit(35); printzeit(36); printf("\n");
printzeit(11); printzeit(29); printzeit(40); printzeit(41); printzeit(42); printzeit(43); printzeit(44); printzeit(50);
// printzeit(32); printzeit(33); printzeit(34); printzeit(35); printzeit(36);
printf("\n");
}

