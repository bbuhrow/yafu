/*
  Copyright (C) 2004 Jens Franke, Torsten Kleinjung
  This file is part of mpqs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

volatile extern u32_t modulo32;
u32_t gcd32(u32_t x, u32_t y);
int jac32(u32_t x,u32_t y);
u32_t modpow32(u32_t x,u32_t a);
u32_t modsqrt32(u32_t x);
u32_t ASM_ATTR asm_modinv32(u32_t x);

#define modinv32(x) asm_modinv32(x)

#if defined(__APPLE__)
#define modulo32 "_modulo32"
#else
#define modulo32 "modulo32"
#endif

static inline u32_t modsq32(u32_t x)
{
  u32_t res,clobber;
  __asm__ volatile ("mull %%eax\n"
	   "divl "modulo32"(%%rip)" : "=d" (res), "=a" (clobber) : "a" (x) : "cc" );
  return res;
}

static inline u32_t modmul32(u32_t x,u32_t y)
{
  u32_t res,clobber;
  __asm__ volatile ("mull %%ecx\n"
	   "divl "modulo32"(%%rip)" : "=d" (res), "=a" (clobber) : "a" (x), "c" (y) :
	   "cc");
  return res;
}

static inline u32_t modadd32(u32_t x,u32_t y)
{
  u32_t res;
#ifdef HAVE_CMOV
  __asm__ volatile ("xorl %%edx,%%edx\n"
	   "addl %%eax,%%ecx\n"
	   "cmovc "modulo32"(%%rip),%%edx\n"
	   "cmpl "modulo32"(%%rip),%%ecx\n"
	   "cmovae "modulo32"(%%rip),%%edx\n"
	   "subl %%edx,%%ecx\n"
	   "2:\n" : "=c" (res) : "a" (x), "c" (y) : "%edx", "cc");
#else
  __asm__ volatile ("addl %%eax,%%ecx\n"
	   "jc 1f\n"
	   "cmpl "modulo32"(%%rip),%%ecx\n"
	   "jb 2f\n"
	   "1:\n"
	   "subl "modulo32"(%%rip),%%ecx\n"
	   "2:\n" : "=c" (res) : "a" (x), "c" (y) : "cc");
#endif
  return res;
}

static inline u32_t modsub32(u32_t subtrahend,u32_t minuend)
{
  u32_t res;
#ifdef HAVE_CMOV
  __asm__ volatile ("xorl %%edx,%%edx\n"
	   "subl %%eax,%%ecx\n"
	   "cmovbl "modulo32"(%%rip),%%edx\n"
	   "addl %%edx,%%ecx\n"
	   "1:" : "=c" (res) : "a" (minuend), "c" (subtrahend) : "%edx", "cc");
#else
  __asm__ volatile ("subl %%eax,%%ecx\n"
	   "jae 1f\n"
	   "addl "modulo32"(%%rip),%%ecx\n"
	   "1:" : "=c" (res) : "a" (minuend), "c" (subtrahend) : "cc" );
#endif
  return res;
}

#undef modulo32
