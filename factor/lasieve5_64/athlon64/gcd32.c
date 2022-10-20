/*
  Copyright (C) 2002 Jens Franke
  This file is part of gnfs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/
#include <sys/types.h>

#include "siever-config.h"

u32_t gcd32(u32_t x, u32_t y)
{
  u32_t r;
  if(!x) return y;
  if(!y) return x;
  __asm__ volatile ( "movl %%esi,%%eax\n"
                     "bsfl %%esi,%%ecx\n"
                     "bsfl %%edi,%%edx\n"
                     "shrl %%cl,%%eax\n"
                     "xchgl %%edx,%%ecx\n"
                     "shrl %%cl,%%edi\n"
                     "subl %%edx,%%ecx\n"
                     "sbbl %%esi,%%esi\n"
                     "andl %%ecx,%%esi\n"
                     "addl %%esi,%%edx\n"
                     "subl %%eax,%%edi\n"
                     "jz 2f\n"
                     "1:\n"
                     "sbbl %%ecx,%%ecx\n"
                     "movl %%ecx,%%esi\n"
                     "andl %%edi,%%esi\n"
                     "xorl %%ecx,%%edi\n"
                     "addl %%esi,%%eax\n"
                     "subl %%ecx,%%edi\n"
                     "bsfl %%edi,%%ecx\n"
                     "shrl %%cl,%%edi\n"
                     "subl %%eax,%%edi\n"
                     "jnz 1b\n"
                     "2:\n"
                     "movl %%edx,%%ecx\n"
                     "shll %%cl,%%eax\n" : "=a"(r) : "S"(y), "D"(x) : "rcx", "rdx");
  return r;
}
