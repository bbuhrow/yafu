@* Input of a pair of NFS polynomials.
@*3 Copying.
Copyright (C) 2001 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@
@(input-poly.h@>=
  void input_poly(mpz_t,mpz_t**,i32_t*,mpz_t**,i32_t*,mpz_t,FILE*);

@
@c
#include <stdio.h>
#include <sys/types.h>
#include <string.h>
#include <gmp.h>
#include <stdlib.h>

#include "asm/siever-config.h"
#include "input-poly.h"
#include "if.h"
#define MAX(_a, _b) ((_a) > (_b) ? (_a) : (_b))

void
input_poly(mpz_t N,mpz_t **A,i32_t *adeg,mpz_t **B,i32_t *bdeg,mpz_t m,
	FILE *fp)
{ char  token[256], value[512], thisLine[1024];
  int   i, j, cont=1;
  mpz_t tmp, tmp2, mpow;

  *adeg = *bdeg = 0;                                                    
  *A = xmalloc(9*sizeof(**A)); /* plenty o' room. */
  *B = xmalloc(9*sizeof(**B));
  for (i=0; i<9; i++) {
    mpz_init_set_ui((*A)[i], 0);
    mpz_init_set_ui((*B)[i], 0);
  }
  while (cont) {
    thisLine[0] = 0;
    fgets(thisLine, 1023, fp);
    if ((sscanf(thisLine, "%255s %511s", token, value)==2) &&
                (thisLine[0] != '#')) {
	  token[sizeof(token)-1] = 0;
      if (strncmp(token, "n:", 2)==0) {
        mpz_set_str(N, value, 10);
      } else if (strncmp(token, "m:", 2)==0) {
        mpz_set_str(m, value, 10);
      } else if ((token[0]=='c') && (token[1] >= '0') && (token[1] <= '8')) {
        mpz_set_str((*A)[token[1]-'0'], value, 10);
        *adeg = MAX(*adeg, token[1]-'0');
      } else if ((token[0]=='Y') && (token[1] >= '0') && (token[1] <= '8')) {
        mpz_set_str((*B)[token[1]-'0'], value, 10);
        *bdeg = MAX(*bdeg, token[1]-'0');
      } else if (strncmp(token, "END_POLY", 8)==0) {
        cont=0;
      }
    }
    if (feof(fp)) cont=0;
  }
  if (*bdeg == 0) {
    mpz_set_ui((*B)[1], 1);
    mpz_neg((*B)[0], m);
    *bdeg=1;
  }

  /* Verify the polynomials: */
  mpz_init(tmp); mpz_init(tmp2); mpz_init(mpow);
  mpz_set_ui(tmp, 0);
  mpz_set_ui(mpow, 1);
  
  if (*bdeg) {
    mpz_neg(m, (*B)[0]);  /* m for temporary use */
    for(i=0; i<=*adeg; i++) {
      mpz_mul(tmp2, mpow, (*A)[i]);
        for(j=i; j <= *adeg; j++)
          mpz_mul(tmp2, tmp2, (*B)[1]);
      mpz_add(tmp, tmp, tmp2);
      mpz_mul(mpow, mpow, m);
    }
  } else {
    for (i=*adeg; i>=0; i--) {
      mpz_mul(tmp, tmp, m);
      mpz_add(tmp, tmp, (*A)[i]);
    }
  }
  mpz_mod(tmp, tmp, N);
  if (mpz_sgn(tmp)) {
    printf("Error: the polynomials don't have a common root:\n");
    for (i=0; i<=*adeg; i++) 
      printf("c%d: %s\n", i, mpz_get_str(token, 10, (*A)[i]));
    for (i=0; i<=*bdeg; i++) 
      printf("Y%d: %s\n", i, mpz_get_str(token, 10, (*B)[i]));
    printf("n: ");
    mpz_out_str(NULL,10,N);
    exit(-1);
  }
  mpz_clear(tmp2); mpz_clear(mpow); mpz_clear(tmp);
}
