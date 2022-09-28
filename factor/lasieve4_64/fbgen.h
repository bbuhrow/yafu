/*
  Copyright (C) 2001 Jens Franke.
  This file is part of gnfs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

u32_t root_finder(u32_t *root_buf,mpz_t *A,u32_t adeg,u32_t p);
u32_t polvalmod32(u32_t * P, u32_t dP, u32_t a);
u32_t polmod32(u32_t * T, mpz_t *A, u32_t dA);
extern u32_t * polr,polr_alloc;

