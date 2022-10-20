/*
  Copyright (C) 2000 Jens Franke
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
u32_t modinv32(u32_t x);
#define asm_modinv32(x) modinv32(x)

static inline u32_t modsq32(u32_t x)
{
  return ((unsigned long long)x*(unsigned long long)x)%modulo32;
}

static u32_t modmul32(u32_t x,u32_t y)
{
  return ((unsigned long long)x*(unsigned long long)y)%modulo32;
}

static u32_t modsub32(u32_t minuend,u32_t subtrahend)
{
  if(subtrahend<=minuend) return minuend-subtrahend;
  else return modulo32+minuend-subtrahend;
}

static inline u32_t modadd32(u32_t x,u32_t y)
{
  return modsub32(x,modulo32-y);
}
