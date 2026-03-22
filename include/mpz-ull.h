//@f mpz int
//@f ullong int

//Combining |mpz| and |ullong|.
//
//Copyright (C) 2002 Jens Franke, T. Kleinjung
//This file is part of gnfs4linux, distributed under the terms of the 
//GNU General Public Licence and WITHOUT ANY WARRANTY.
//
//You should have received a copy of the GNU General Public License along
//with this program; see the file COPYING.  If not, write to the Free
//Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
//02111-1307, USA.

#include <sys/types.h>
#include <limits.h>
typedef unsigned long ulong;
typedef unsigned long long ullong;
#include <gmp.h>

#define mpz_add_si(r,o,s) \
  ({ int _s; _s= (s); \
     _s>=0 ? mpz_add_ui(r,o,(mp_limb_t)(_s)) : mpz_sub_ui(r,o,(mp_limb_t)(-_s)); })


void mpz_set_ull(mpz_t targ, ullong src);
ullong mpz_get_ull(mpz_t src);
int mpz_cmp_ull(mpz_t op1, ullong op2);

#ifdef ULL_NO_UL
void mpz_ull_init();
void mpz_mul_ull(mpz_t rop, mpz_t op1, ullong op2);
int mpz_fits_sllong_p(mpz_t);
int mpz_fits_ullong_p(mpz_t);
long long int mpz_get_sll(mpz_t);
void mpz_set_sll(mpz_t, long long int);
unsigned long long mpz_get_ull(mpz_t);
void mpz_set_ull(mpz_t, unsigned long long);
void mpz_add_ull(mpz_t rop, mpz_t op1, ullong op2);
void mpz_tdiv_q_ull(mpz_t rop, mpz_t op1, ullong op2);
#else

#define mpz_ull_init()
#define mpz_mul_ull mpz_mul_ui
#define mpz_mul_sll mpz_mul_si
#define mpz_fits_sllong_p mpz_fits_slong_p
#define mpz_fits_ullong_p mpz_fits_ulong_p
#define mpz_get_sll mpz_get_si
#define mpz_set_sll mpz_set_si
#define mpz_get_ull mpz_get_ui
#define mpz_set_ull mpz_set_ui
#define mpz_add_ull mpz_add_ui
#define mpz_tdiv_q_ull mpz_tdiv_q_ui
#endif

