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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "asm/siever-config.h"
#include "if.h"

#include "asm/montgomery_mul.h"



int main(int argc, char *argv[])
{
  mpz_t N;
  u32_t n, i, necms=0, iter=0, B1, B2;
  FILE *fi;
  mpz_t *f;
  char *input_line=NULL;
  size_t input_line_alloc=0;
  u32_t imax;
  int nc;

  setbuf(stdout,NULL);
  if (argc<4) complain("Usage: ecmtest B1 filename ncurves [ B2 ]\n");
  B1=(u32_t)atoi(argv[1]);
  initzeit(13); zeita(0);
  mpz_init(N); n=0;
  fi=fopen(argv[2],"r");
  imax=strtoul(argv[3],NULL,10);
  if (argc>4) B2=strtoul(argv[4],NULL,10); else B2=80*B1;
  while (getline(&input_line,&input_line_alloc,fi)>0) { iter++;
    if (mpz_set_str(N,input_line,10)) break;
    nc=ecm_factor(N,B1,B2,&f,imax);
    if (nc>0) { necms+=(u32_t)nc; n++; /* gmp_printf("%Zd  ",*f); */ } else necms+=imax;
    printf(".");
  }
  fclose(fi);
  if (n) {
    printf("ECM-tests: %u  Factors: %u\n",necms,n);
  } else printf("No factor found!\n");
  printf("\n");
  zeitb(0);
  for (i=0; i<13; i++) printzeit(i); printf("\n");
  exit(0);
}

