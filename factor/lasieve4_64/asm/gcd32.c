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
#include "32bit.h"

u32_t gcd32(u32_t x, u32_t y)
{
  u32_t f;
  if(!x) return y;
  if(!y) return x;
  for(f=0;!(x&1) && !(y&1);x>>=1,y>>=1) f++;
  for(;!(x&1);x>>=1);
  for(;!(y&1);y>>=1);
  while(x!=y) {
    if(x<y) {
      y-=x;
      for(;!(y&1);y>>=1);
    } else {
      x-=y;
      for(;!(x&1);x>>=1);
    }
  }
  return x<<f;
}
