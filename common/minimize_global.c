/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: minimize_global.c 296 2010-06-07 01:29:37Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>

/*-------------------------------------------------------------------------*/
#define IMPULSE() (-impulse * log(((double)get_rand(seed1, \
				seed2) + 1) / MP_RADIX))

static double simplex_core(double p[MAX_VARS+1][MAX_VARS], 
				double y[MAX_VARS+1], double psum[MAX_VARS], 
				double pbest[MAX_VARS], double *ybest, 
				double *yhi, int ndim, int ihi, double impulse, 
				uint32 *seed1, uint32 *seed2, double fac,
				objective_func callback, void *extra)
{
	int i;
	double fac1, fac2, ytry, yimpulse; 
	double ptry[MAX_VARS];

	fac1 = (1.0 - fac) / ndim;
	fac2 = fac1 - fac;

	for (i = 0; i < ndim; i++)
		ptry[i] = psum[i] * fac1 - p[ihi][i] * fac2;

	ytry = callback(ptry, extra);
	yimpulse = ytry - IMPULSE();

	if (ytry <= *ybest) {
		*ybest = ytry;
		for (i = 0; i < ndim; i++)
			pbest[i] = ptry[i];
	}
	if (yimpulse < *yhi) {
		*yhi = yimpulse;
		y[ihi] = ytry;
		for (i = 0; i < ndim; i++) {
			psum[i] += ptry[i] - p[ihi][i];
			p[ihi][i] = ptry[i];
		}
	}
	return ytry;
}

/*-------------------------------------------------------------------------*/
static int32 simplex(double p[MAX_VARS+1][MAX_VARS], double y[MAX_VARS+1], 
			double pbest[MAX_VARS], double *ybest,
			double impulse, uint32 *seed1, 
			uint32 *seed2, int ndim, double tol,
			int iter, objective_func callback, void *extra)
{
	double psum[MAX_VARS];
	double rtol, ytry, yhi, ynhi, ylo;
	int32 i, j;
	int32 ihi, ilo;

	for (i = 0; i < ndim; i++) {
		double sum = 0.0;
		for (j = 0; j < ndim + 1; j++)
			sum += p[j][i];
		psum[i] = sum;
	}

	while (iter > 0) {

		ilo = 0;
		ihi = 1;
		ynhi = ylo = y[0] + IMPULSE();
	        yhi = y[1] + IMPULSE();	
		if (ylo > yhi) {
			ihi = 0; 
			ilo = 1;
			ynhi = yhi;
			yhi = ylo;
			ylo = ynhi;
		}
		for (i = 2; i < ndim + 1; i++) {
			double yt = y[i] + IMPULSE();
			if (yt <= ylo) {
				ilo = i;
				ylo = yt;
			}
			if (yt > yhi) {
				ynhi = yhi;
				ihi = i;
				yhi = yt;
			}
			else if (yt > ynhi) {
				ynhi = yt;
			}
		}

		rtol = 2.0 * fabs(yhi - ylo) /
				(fabs(yhi) + fabs(ylo));
		if (rtol < tol)
			return 0;

		ytry = simplex_core(p, y, psum, pbest, ybest, &yhi,
					ndim, ihi, impulse, seed1,
					seed2, -1.0, callback, extra);
		iter--;

		if (ytry <= ylo) { 
			ytry = simplex_core(p, y, psum, pbest, ybest, &yhi,
						ndim, ihi, impulse, seed1,
						seed2, 2.0, callback, extra);
			iter--;
		}
		else if (ytry >= ynhi) {

			double ysave = yhi;
			ytry = simplex_core(p, y, psum, pbest, ybest, &yhi,
						ndim, ihi, impulse, seed1,
						seed2, 0.5, callback, extra);
			iter--;

			if (ytry >= ysave) {
				for (i = 0; i < ndim + 1; i++) {
					if (i == ilo)
						continue;
					for (j = 0; j < ndim; j++) {
						p[i][j] = 0.5 * (p[i][j] +
								p[ilo][j]);
					}
					y[i] = callback(p[i], extra);
				}
				iter -= ndim;

				for (i = 0; i < ndim; i++) {
					double sum = 0.0;
					for (j = 0; j < ndim + 1; j++)
						sum += p[j][i];
					psum[i] = sum;
				}
			}
		}
	}

	return -1;
}

/*-------------------------------------------------------------------------*/
double minimize_global(double p[MAX_VARS], uint32 ndim,
			double limits[MAX_VARS][2],
			double tol, uint32 iter_limit,
			objective_func callback, 
			void *extra)
{
	uint32 i, j;
	double simp[MAX_VARS+1][MAX_VARS];
	double scores[MAX_VARS+1];
	double impulse, best_score;
	uint32 seed1 = 1111111;
	uint32 seed2 = 2222222;
	uint32 iter = 0;

	best_score = 1e100;

	for (i = 0; i < ndim + 1; i++) {
		for (j = 0; j < ndim; j++) {
			uint32 val = get_rand(&seed1, &seed2);
			simp[i][j] = limits[j][0] + 
				     (double)val / 4294967296.0 *
				     (limits[j][1] - limits[j][0]);
		}
		scores[i] = callback(simp[i], extra);
		if (scores[i] < best_score)
			best_score = scores[i];
	}

	impulse = best_score * 0.1;

	while (iter < iter_limit) {

		iter += 80;

		if (simplex(simp, scores, p, &best_score, 
			impulse, &seed1, &seed2, ndim, tol, 80,
			callback, extra) == 0)
			break;

		impulse *= 0.99;
	}

	return best_score;
}

