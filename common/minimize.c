/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: minimize.c 734 2012-08-04 15:13:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include "ms_common.h"

#define PHI 1.618033988749
#define COMP_PHI (2.0 - PHI)

//#define MSIEVE_STATIC static
#define MSIEVE_STATIC 

/*-------------------------------------------------------------------------*/
MSIEVE_STATIC double evaluate(double *base, double dist, 
			double *search_dir, uint32 ndim,
			objective_func callback, void *extra)
{
	uint32 i;
	double curr_pt[MAX_VARS];

	for (i = 0; i < ndim; i++)
		curr_pt[i] = base[i] + dist * search_dir[i];

	return callback(curr_pt, extra);
}

/*-------------------------------------------------------------------------*/
MSIEVE_STATIC void bracket_min(double *base, double *search_dir,
			double *a_out, double *b_out, double *c_out,
			double tol, double *fb_out, uint32 ndim, 
			objective_func callback, void *extra)
{
	double a = 0;
	double fa = evaluate(base, a, search_dir, ndim, callback, extra);
	double b = 1;
	double fb = evaluate(base, b, search_dir, ndim, callback, extra);
	double c;
	double fc; 
	double t;

	while (fabs(fb - fa) < 10 * tol * fabs(fa)) {
		b *= 2;
		fb = evaluate(base, b, search_dir, ndim, callback, extra);
	}

	if (fb > fa) {
		t = a; a = b; b = t;
		t = fa; fa = fb; fb = t;
	}
	c = b + PHI * (b - a);
	fc = evaluate(base, c, search_dir, ndim, callback, extra);

	while (fb > fc) {
		double r = (b - a) * (fb - fc);
		double q = (b - c) * (fb - fa);
		double u, max_u, fu;

		t = fabs(q - r);
		t = MAX(t, 1e-20);
		if (q - r < 0.0)
			t = -t;

		u = b - ((b - c) * q - (b - a) * r) / (2.0 * t);
		max_u = b + 100.0 * (c - b);

		if ((b - u) * (u - c) > 0.0) {
			fu = evaluate(base, u, search_dir, ndim, 
					callback, extra);

			if (fu < fc) {
				*a_out = b; *b_out = u; *c_out = c;
				*fb_out = fu;
				return;
			}
			else if (fu > fb) {
				*a_out = a; *b_out = b; *c_out = u;
				*fb_out = fb;
				return;
			}

			u = c + PHI * (c - b);
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);
		}
		else if ((c - u) * (u - max_u) > 0.0) {
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);

			if (fu < fc) {
				b = c; c = u; u = c + PHI * (c - b);
				fb = fc; 
				fc = fu; 
				fu = evaluate(base, u, search_dir, ndim,
						callback, extra);
			}
		}
		else if ((u - max_u) * (max_u - c) >= 0.0) {
			u = max_u;
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);
		}
		else {
			u = c + PHI * (c - b);
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);
		}

		a = b; b = c; c = u;
		fa = fb; fb = fc; fc = fu;
	}

	*a_out = a; *b_out = b; *c_out = c;
	*fb_out = fb;
}

/*-------------------------------------------------------------------------*/
MSIEVE_STATIC double minimize_line_core(double *base, double *search_dir,
			double a_in, double b_in, double c_in, 
			double fb_in, double tol, double *min_out, 
			uint32 ndim, int32 *status,
			objective_func callback, void *extra)
{
	int32 i;
	double a, b, x, w, v, u;
	double fx, fw, fv, fu;
	double e = 0.0;

	a = MIN(a_in, c_in);
	b = MAX(a_in, c_in);
	x = w = v = b_in;
	fx = fw = fv = fb_in;

	for (i = 0; i < 300; i++) {
		double xm = 0.5 * (a + b);
		double tol1 = tol * fabs(x) + 1e-20;
		double tol2 = 2.0 * tol1;
		double d = 0, etemp;

		if (fabs(x - xm) <= tol2 - 0.5 * (b - a)) {
			*min_out = x;
			return fx;
		}

		if (fabs(e) > tol1) {
			double r = (x - w) * (fx - fv);
			double q = (x - v) * (fx - fw);
			double p = (x - v) * q - (x - w) * r;
			
			q = 2.0 * (q - r);
			if (q > 0.0)
				p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) ||
			    p <= q * (a - x) || p >= q * (b - x)) {
	
				e = b - x;
				if (x >= xm)
					e = a - x;				
				d = COMP_PHI * e;
			}
			else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = (xm >= x) ? tol1 : -tol1;
			}
		}
		else {
			e = (x >= xm) ? a - x : b - x;
			d = COMP_PHI * e;
		}

		if (fabs(d) < tol1)
			u = x + ((d >= 0) ? tol1 : -tol1);
		else
			u = x + d;
		fu = evaluate(base, u, search_dir, ndim, callback, extra);

		if (fu <= fx) {
			if (u >= x)
				a = x;
			else
				b = x;

			v = w; w = x; x = u;
			fv = fw; fw = fx; fx = fu;
		}
		else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x) {
				v = w; w = u;
				fv = fw; fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}

	printf(":"); // "too many line iterations\n");
	*min_out = x;
	*status = 1;
	return fx;
}

/*-------------------------------------------------------------------------*/
MSIEVE_STATIC double minimize_line(double *base, double *search_dir, 
				uint32 ndim, double ftol,
				int32 *status,
				objective_func callback,
				void *extra)
{
	uint32 i;
	double a;
	double b, fb;
	double c;
	double min_dist, fmin_dist;
	double max_dir;
	double normalized[MAX_VARS];
	double backup[MAX_VARS];

	for (i = 0; i < ndim; i++)
		backup[i] = base[i];

	max_dir = fabs(search_dir[0]);
	for (i = 1; i < ndim; i++)
		max_dir = MAX(max_dir, fabs(search_dir[i]));
	for (i = 0; i < ndim; i++)
		normalized[i] = search_dir[i] / max_dir;

	bracket_min(base, normalized, &a, &b, &c, 
			ftol, &fb, ndim, callback, extra);
	fmin_dist = minimize_line_core(base, normalized, a, b, c, 
					fb, ftol, &min_dist, ndim,
					status, callback, extra);
	if (*status) {
		for (i = 0; i < ndim; i++)
			base[i] = backup[i];
	}
	else {
		for (i = 0; i < ndim; i++) {
			search_dir[i] = normalized[i] * min_dist;
			base[i] += search_dir[i];
		}
	}

	return fmin_dist;
}

/*-------------------------------------------------------------------------*/
double minimize(double p[MAX_VARS], uint32 ndim,
	double ftol, uint32 max_iter,
	objective_func callback, void* extra)
{
	uint32 i, j;
	double best_f, curr_f;
	double p_old[MAX_VARS], p_extrap[MAX_VARS], dir_extrap[MAX_VARS];
	double directions[MAX_VARS][MAX_VARS];
	int32 status = 0;

	memset(directions, 0, sizeof(directions));
	for (i = 0; i < ndim; i++) {
		p_old[i] = p[i];
		directions[i][i] = 1.0;
	}
	best_f = callback(p, extra);

	if (ndim == 1) {
		curr_f = minimize_line(p, directions[0], ndim,
			MAX(ftol, 1e-7), &status,
			callback, extra);
		if (status)
			return best_f;
		return curr_f;
	}

	for (i = 0; i < max_iter; i++) {

		int32 ibig = 0;
		double del = 0.0;
		double start_f = best_f;

		for (j = 0; j < ndim; j++) {
			// when called from poly_sizeopt, ndim is the following:
			// uint32 rotate_dim = deg - 4;
			// uint32 num_vars = rotate_dim + 3;
			curr_f = minimize_line(p, directions[j], ndim,
				MAX(ftol, 1e-7), &status,
				callback, extra);
			if (status)
				return best_f;

			if (best_f - curr_f > del) {
				ibig = j;
				del = best_f - curr_f;
			}
			best_f = curr_f;
		}

		if (2.0 * fabs(start_f - best_f) <=
			ftol * (fabs(start_f) + fabs(best_f)) + 1e-20)
			return best_f;

		for (j = 0; j < ndim; j++) {
			p_extrap[j] = 2.0 * p[j] - p_old[j];
			dir_extrap[j] = p[j] - p_old[j];
			p_old[j] = p[j];
		}

		curr_f = callback(p_extrap, extra);
		if (curr_f < start_f) {

			double t1 = start_f - best_f - del;
			double t2 = start_f - curr_f;
			double t = 2.0 * (start_f - 2.0 * best_f + curr_f) *
				t1 * t1 - del * t2 * t2;

			if (t < 0.0) {
				best_f = minimize_line(p, dir_extrap, ndim,
					MAX(ftol, 1e-7), &status,
					callback, extra);
				if (status)
					return best_f;

				for (j = 0; j < ndim; j++) {
					directions[ibig][j] =
						directions[ndim - 1][j];
					directions[ndim - 1][j] =
						dir_extrap[j];
				}
			}
		}
	}

	return best_f;
}

/*-------------------------------------------------------------------------*/
double minimize_grad(double p[MAX_VARS], uint32 ndim, 
			double ftol, uint32 max_iter,
			objective_func callback, 
			objective_func_grad callback_grad,
			void *extra)
{
	uint32 i, j;
	double g[MAX_VARS];
	double h[MAX_VARS];
	double grad[MAX_VARS];
	double best_f;

	best_f = callback_grad(p, grad, extra);
	for (i = 0; i < ndim; i++)
		g[i] = h[i] = grad[i] = -grad[i];

	for (i = 0; i < max_iter; i++) {

		double gam;
		int32 status = 0;
		double gg = 0;
		double dgg = 0;
		double gmax, gradmax;
		double curr_f = minimize_line(p, grad, ndim,
					MAX(ftol, 1e-7), &status, 
					callback, extra);

		if (status)
			return curr_f;
		if (2.0 * fabs(best_f - curr_f) <= ftol *
				(fabs(best_f) + fabs(curr_f) + 1e-10))
			return curr_f;

		best_f = callback_grad(p, grad, extra);

		gmax = fabs(g[0]);
		gradmax = fabs(grad[0]);
		for (j = 1; j < ndim; j++) {
			gmax = MAX(gmax, fabs(g[j]));
			gradmax = MAX(gradmax, fabs(grad[j]));
		}
		if (gmax == 0)
			return best_f;

		for (j = 0; j < ndim; j++) {
			gg += g[j] * (g[j] / gmax);
			dgg += (grad[j] + g[j]) * (grad[j] / gradmax);
		}

		gam = (dgg / gg) * (gmax / gradmax);
		for (j = 0; j < ndim; j++) {
			g[j] = -grad[j];
			grad[j] = h[j] = -grad[j] + gam * h[j];
		}
	}

	return best_f;
}

/*-------------------------------------------------------------------------*/
#define SWAP(type, a, b) {type tmp = (a); (a) = (b); (b) = tmp; }

void 
solve_dmatrix(double matrix[MAX_VARS][MAX_VARS], 
		double x[MAX_VARS], 
		double b[MAX_VARS], 
		uint32 n) {

	int32 i, j, k;
	uint32 permute[MAX_VARS];
	double pivot;

	for (i = 0; i < n; i++)
		permute[i] = i;

	for (i = 0; i < n - 1; i++) {

		uint32 pivot_idx = i;
		double *pivot_row;

		pivot = fabs(matrix[permute[i]][i]);

		for (j = i + 1; j < n; j++) {
			uint32 new_idx = j;
			double new_pivot = fabs(matrix[permute[j]][i]);

			if (new_pivot > pivot) {
				pivot_idx = new_idx;
				pivot = new_pivot;
			}
		}
		SWAP(int, permute[i], permute[pivot_idx]);
		pivot_row = matrix[permute[i]];
		pivot = 1.0 / pivot_row[i];

		for (j = i + 1; j < n; j++) {
			double *curr_row = matrix[permute[j]];
			double mult = curr_row[i] * pivot;

			for (k = i + 1; k < n; k++) {
				curr_row[k] -= mult * pivot_row[k];
			}
			b[permute[j]] -= mult * b[permute[i]];
		}
	}

	for (i = n - 1; i >= 0; i--) {

		double *curr_row = matrix[permute[i]];
		double sum = b[permute[i]];

		for (j = i + 1; j < n; j++) {
			sum -= x[j] * curr_row[j];
		}

		x[i] = sum / curr_row[i];
	}
}

/*-------------------------------------------------------------------------*/
double minimize_hess(double p[MAX_VARS], uint32 ndim, 
			double ftol, uint32 max_iter,
			objective_func callback, 
			objective_func_hess callback_hess,
			void *extra)
{
	uint32 i, j;
	double grad[MAX_VARS];
	double pnew[MAX_VARS];
	double search_dir[MAX_VARS];
	double hess[MAX_VARS][MAX_VARS];
	double best_f, curr_f;

	for (i = 0, curr_f = 0.0; i < max_iter; i++) {

		int32 status = 0;

		best_f = callback_hess(p, grad, hess, extra);

		for (j = 0; j < ndim; j++) {
			grad[j] = -grad[j];
			pnew[j] = p[j];
		}

		solve_dmatrix(hess, search_dir, grad, ndim);

		curr_f = minimize_line(pnew, search_dir, ndim,
					MAX(ftol, 1e-7), &status, 
					callback, extra);

		if (curr_f < best_f) {
			for (j = 0; j < ndim; j++)
				p[j] = pnew[j];
		}
		else {
			printf("line minimize failed\n");
			curr_f = best_f;
		}

		if (status)
			return curr_f;

		if (2.0 * fabs(best_f - curr_f) <= ftol *
				(fabs(best_f) + fabs(curr_f) + 1e-10))
			break;
	}

	return curr_f;
}

/*-------------------------------------------------------------------------*/
MSIEVE_STATIC double minimize_line_quasinewt(
		double base[MAX_VARS], 
		uint32 ndim,
		double fbase, 
		double gradbase[MAX_VARS],
		double search_dir[MAX_VARS], 
		double new_base[MAX_VARS],
		double ftol,
		double max_step, 
		uint32 *status,
		objective_func callback, 
		void *extra)
{
	uint32 i;
	double sum;
	double slope;
	double test;
	double alamin, alam, alam2 = 0;
	double curr_f;
	double f2 = 0;

	*status = 0;

	for (i = 0, sum = 0.0; i < ndim; i++)
		sum += search_dir[i] * search_dir[i];
	sum = sqrt(sum);

	if (sum > max_step) {
		for (i = 0; i < ndim; i++)
			search_dir[i] *= max_step / sum;
	}

	for (i = 0, slope = 0; i < ndim; i++)
		slope += gradbase[i] * search_dir[i];
	if (slope >= 0) {
		for (i = 0; i < ndim; i++)
			new_base[i] = base[i];
		*status = 1;
		return fbase;
	}

	for (i = 0, test = 0; i < ndim; i++) {
		double t = fabs(search_dir[i]) / 
				MAX(1.0, fabs(base[i]));
		test = MAX(test, t);
	}

	alamin = ftol / test;
	alam = 1.0;
	while (1) {

		double tmplam; 
		double min_decrease = 1e-4 * alam * slope;

		for (i = 0; i < ndim; i++)
			new_base[i] = base[i] + alam * search_dir[i];

		curr_f = callback(new_base, extra);
		if (alam < alamin) {
			*status = 2;
			return curr_f;
		}
		else if (fabs(fbase) < fabs(min_decrease) &&
			 fabs(curr_f) < 1e-4 * fabs(fbase)) {
			return curr_f;
		}
		else if (curr_f < fbase + min_decrease) {
			return curr_f;
		}
		else {
			if (alam == 1.0) {
				tmplam = -slope / (2 * 
					(curr_f - fbase - slope));
			}
			else {
				double disc, rhs1, rhs2, a, b;

				rhs1 = curr_f - fbase - alam * slope;
				rhs2 = f2 - fbase - alam2 * slope;
				a = (rhs1 / (alam * alam) -
				     rhs2 / (alam2 * alam2)) / 
				    (alam - alam2);
				b = (-alam2 * rhs1 / (alam * alam) +
				     alam * rhs2 / (alam2 * alam2)) /
				    (alam - alam2);

				if (a == 0.0) {
					tmplam = -slope / (2.0 * b);
				}
				else {
					disc = b * b - 3.0 * a * slope;
					if (disc < 0.0) {
						tmplam = 0.5 * alam;
					}
					else if (b <= 0.0) {
						tmplam = (-b + sqrt(disc)) / 
							(3.0 * a);
					}
					else {
						tmplam = -slope / 
							(b + sqrt(disc));
					}

					tmplam = MIN(tmplam, 0.5 * alam);
				}
			}
		}

		alam2 = alam;
		f2 = curr_f;
		alam = MAX(tmplam, 0.1 * alam);
	}

	return curr_f;
}

/*-------------------------------------------------------------------------*/
double minimize_grad_bfgs(double p[MAX_VARS], uint32 ndim, 
			double ftol, uint32 max_iter,
			objective_func callback, 
			objective_func_grad callback_grad,
			void *extra)
{
	uint32 i, j, k;
	double pnew[MAX_VARS];
	double g[MAX_VARS];
	double dg[MAX_VARS];
	double hdg[MAX_VARS];
	double xi[MAX_VARS];
	double hess[MAX_VARS][MAX_VARS];
	double curr_f, best_f;
	double sum, max_step;
	uint32 status;

	curr_f = best_f = callback_grad(p, g, extra);

	memset(hess, 0, sizeof(hess));
	for (i = 0, sum = 0.0; i < ndim; i++) {
		hess[i][i] = 1.0;
		xi[i] = -g[i];
		sum += p[i] * p[i];
	}

	max_step = 100.0 * MAX(sqrt(sum), (double)ndim);

	for (i = 0; i < max_iter; i++) {

		double test, den;
		double fac, fad, fae, sumdg, sumxi;

		curr_f = minimize_line_quasinewt(p, ndim, best_f, g,
						xi, pnew, ftol, max_step,
						&status, callback, extra);
		for (j = 0; j < ndim; j++) {
			xi[j] = pnew[j] - p[j];
			p[j] = pnew[j];
		}

		if (status != 0)
			return curr_f;

		for (j = 0, test = 0; j < ndim; j++) {
			double t = fabs(xi[j]) / MAX(fabs(p[j]), 1.0);
			test = MAX(test, t);
		}
		if (test < 4 * ftol)
			return curr_f;

		for (j = 0; j < ndim; j++)
			dg[j] = g[j];
		best_f = callback_grad(p, g, extra);

		den = MAX(best_f, 1.0);
		for (j = 0, test = 0.0; j < ndim; j++) {
			double t = fabs(g[j]) * MAX(p[j], 1.0) / den;
			test = MAX(test, t);
		}
		if (test < ftol)
			return best_f;

		for (j = 0; j < ndim; j++)
			dg[j] = g[j] - dg[j];

		for (j = 0; j < ndim; j++) {
			for (k = 0, hdg[j] = 0; k < ndim; k++) {
				hdg[j] += hess[j][k] * dg[k];
			}
		}

		fac = fae = sumdg = sumxi = 0.0;
		for (j = 0; j < ndim; j++) {
			fac += dg[j] * xi[j];
			fae += dg[j] * hdg[j];
			sumdg += dg[j] * dg[j];
			sumxi += xi[j] * xi[j];
		}

		if (fac > sqrt(1e-15 * sumdg * sumxi)) {
			fac = 1.0 / fac;
			fad = 1.0 / fae;
			for (j = 0; j < ndim; j++)
				dg[j] = fac * xi[j] - fad * hdg[j];

			for (j = 0; j < ndim; j++) {
				for (k = j; k < ndim; k++) {
					hess[j][k] += fac * xi[j] * xi[k] -
						      fad * hdg[j] * hdg[k] +
						      fae * dg[j] * dg[k];
					hess[k][j] = hess[j][k];
				}
			}
		}

		for (j = 0; j < ndim; j++) {
			for (k = 0, xi[j] = 0.0; k < ndim; k++) {
				xi[j] -= hess[j][k] * g[k];
			}
		}
	}

	return curr_f;
}

