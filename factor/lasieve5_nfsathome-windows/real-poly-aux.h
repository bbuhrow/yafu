/*
  Copyright (C) 2001 Jens Franke.
  This file is part of gnfs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

void tpol(double *rop,double *op,i32_t deg,i32_t x0,i32_t x1,i32_t y0,i32_t y1);
void tpol64(double *rop,double *op,i32_t deg,i64_t x0,i64_t x1,
            i64_t y0,i64_t y1);
double rpol_eval(double *p,i32_t d,double x,double y);
void rpol_eval_clear(void);
void get_sieve_report_bounds(unsigned char**,double*,i32_t,
			     i32_t,i32_t,i32_t,double,double);

double rpol_eval0(double *p,i32_t d,i32_t x0,u16_t y0);

