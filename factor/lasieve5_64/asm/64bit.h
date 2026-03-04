/*
  Copyright (C) 2000,2006 Jens Franke, Thorsten Kleinjung
  This file is part of gnfs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

volatile extern u64_t modulo64;

u64_t modmul64(u64_t x,u64_t y);

#if 1
#define modmul64 asm_modmul64
#define modadd64 asm_modadd64
#else
u64_t modmul64(u64_t x,u64_t y);

static inline u64_t modadd64(u64_t x,u64_t y)
{
  u64_t res;

  res=x+y; if (res<x) res-=modulo64;
  if (res>=modulo64) res-=modulo64;
  return res;
}
#endif
