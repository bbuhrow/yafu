/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/

/* ECM on number N using parameters B1, B2 and curves number
   0,1,...,ncurves-1. Returns 0 iff no factor was found, -1 on
   error and on success the number of curves processed. The factor
   is put at position 0 of the array of size 1 in *fptr. The factor
   is a divisor of N and not 1, but it might be N, and it might be
   composite. */
int ecm_factor(mpz_t N, u32_t B1, u32_t B2, mpz_t **fptr, u32_t ncurves);

/* If you do not want to process all curves at once use the following
   functions: CAVE they will change */

typedef struct {
  u32_t B1;
  u32_t B2_index;
  mpz_t N;
  mpz_t curve_x;
  mpz_t curve_y;
  ulong curve_param[NMAX_ULONGS];

/* todo mm-data? ... */
} ecm_struct;

typedef ecm_struct ecm_t[1];


/* construct structure for doing ECM */
void ecm_curve_init(ecm_t e);

/* set input number to N and parameters to B1, B2, reset curve number.
   Returns -1 on error, otherwise 0. */
int ecm_curve_set(ecm_t e, mpz_t N, u32_t B1, u32_t B2);

/* change parameters B1,B2 but do not reset curve number */
void ecm_set_params(ecm_t e, u32_t B1, u32_t B2);

/* replace input number by N. N must be a divisor of input number.
   B1, B2 and curve number are unchanged. Returns -1 on error, otherwise 0. */
int ecm_reset_n(ecm_t e, mpz_t N);

/* process next curve, return -1 on error, 1 on success, 0 otherwise */
int ecm(ecm_t e, mpz_t **fptr);

/* free allocated memory: */
void ecm_curve_clear(ecm_t e);



