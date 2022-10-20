/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/

#include <stdlib.h>
#include <gmp.h>
#include "siever-config.h"
#include "../if.h"

u64_t ASM_ATTR pt64(u64_t);

#include "montgomery_mul.h"

extern ulong *montgomery_modulo_n;
extern ulong montgomery_modulo_R2[NMAX_ULONGS];
extern ulong montgomery_ulongs;


int psp2(mpz_t n)
{
  ulong x[NMAX_ULONGS], ex[NMAX_ULONGS], one[NMAX_ULONGS], k;
  u32_t i, j, e, s;

  if (!set_montgomery_multiplication(n)) return -1;

  if (!(montgomery_modulo_n[0]&1)) return 0;  /* number is even */
  asm_copy(ex,montgomery_modulo_n); ex[0]--;
  for (i=0; i<montgomery_ulongs; i++) if (ex[i]) break;
  if (i>=montgomery_ulongs) return 0;  /* number is 1 */

  e=0;
  while (!(ex[0]&1)) {
    for (i=0; i<montgomery_ulongs-1; i++) ex[i]=(ex[i]>>1)|(ex[i+1]<<63);
    ex[montgomery_ulongs-1]>>=1;
    e++;
  }
  asm_zero(one); one[0]=1;
  asm_mulmod(one,montgomery_modulo_R2,one);

  for (i=1,j=0; i<montgomery_ulongs; i++) if (ex[i]) j=i;
  if ((j==0) && (ex[0]==0)) complain("psp\n");
  for (i=1,k=0; i<64; i++) if (ex[j]&(1ULL<<i)) k=1ULL<<i;

  s=1;
  while (j|(k-1)) {
    k>>=1;
    if (!k) { j--; k=1ULL<<63; }
    s<<=1;
    if (ex[j]&k) s++;
    if (s>32*(montgomery_ulongs-1)) break;
  }
  asm_zero(x);
  x[s>>6]=1ULL<<(s&63);
  asm_mulmod(x,montgomery_modulo_R2,x);

  while (j|(k-1)) {
    k>>=1;
    if (!k) { j--; k=1ULL<<63; }
#if 1
    asm_squmod(x,x);
#else
    asm_mulmod(x,x,x);
#endif
    if (ex[j]&k) asm_add2(x,x);
  }

  if (asm_cmp(x,one)==0) return 1;
  asm_diff(one,montgomery_modulo_n,one);
  if (asm_cmp(x,one)==0) return 1;
  for (i=0; i<e-1; i++) {
#if 1
    asm_squmod(x,x);
#else
    asm_mulmod(x,x,x);
#endif
    if (asm_cmp(x,one)==0) return 1;
  }
  return 0;
}


int psp(mpz_t n)
{
  int res0, res1;

  if(mpz_size(n)<2) {
    if(mpz_sgn(n)!=0) return pt64(n[0]._mp_d[0]);
    return 0;
  }
/*
  if (mpz_size(n)==2) {
    res0=psp2(n);
    res1=mpz_probab_prime_p(n,1);
    if (res0!=res1) {
      printf("psp: %d %d\n",res0,res1);
      gmp_printf("%Zd\n",n);
      complain("");
    }
    return psp2(n);
  }
*/
  if (mpz_size(n)==2) return psp2(n);
  return mpz_probab_prime_p(n,1);
}

