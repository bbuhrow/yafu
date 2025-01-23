/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: polyroot.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include <polyroot.h>

/* This rootfinder uses a simplified version of the all-complex
   Jenkins-Traub algorithm, then switches to Newton's method
   in extended precision to polish the roots that are found. In
   practice we cannot use only Newton's method, because this
   behaves essentially randomly unless given a very good initial
   root approximation. Jenkins-Traub is much more complex, and 
   converges just about unconditionally.

   The J-T code was kindly converted by Brian Gladman from the
   awful Fortran code of Algorithm 419 from comm. acm, vol. 15, 
   no. 0. I've removed a lot of the gold-plating of the full 
   implementation that is unlikely to be necessary for NFS 
   polynomials, and converted the paired-array-of-doubles code
   to much simpler array-of-complex code */

/* Jenkins-Traub rootfinder */

typedef struct {
	double r, i;
} complex_t;

static const double are = DBL_EPSILON;
static const double mre = 2.0 * M_SQRT2 * DBL_EPSILON;
static const complex_t czero = { 0, 0 };

typedef struct {
	complex_t poly[MAX_ROOTFINDER_DEGREE + 1]; 
	complex_t poly_aux[MAX_ROOTFINDER_DEGREE + 1];
	complex_t hpoly[MAX_ROOTFINDER_DEGREE + 1]; 
	complex_t hpoly_aux[MAX_ROOTFINDER_DEGREE + 1];

	complex_t angle, angle_inc;
	uint32 hpoly_root_found;
	uint32 degree;
} jt_t;

/*-----------------------------------------------------------------------*/
static complex_t complex(double re, double im) {

	complex_t res;

	res.r = re;
	res.i = im;
	return res;
}

static complex_t cadd(complex_t x, complex_t y) {
	return complex(x.r + y.r, x.i + y.i);
}

static complex_t cneg(complex_t x) {
	return complex(-x.r, -x.i);
}

static complex_t cscale(double a, complex_t b) {
	return complex(a * b.r, a * b.i);
}

static complex_t cmul(complex_t a, complex_t b) {
	return complex(a.r * b.r - a.i * b.i, 
			a.r * b.i + a.i * b.r);
}

static complex_t cmac(complex_t a, complex_t x, complex_t b) {
	return complex(a.r * x.r - a.i * x.i + b.r,
		       a.r * x.i + a.i * x.r + b.i);
}

static double cmod(complex_t x) {

	/* modulus of a complex number avoiding overflow */

	double re = fabs(x.r); 
	double im = fabs(x.i); 

	if (re == im)
		return re * M_SQRT2;

	if (re < im) {
		double t = re / im;
		return im * sqrt(1.0 + t * t);
	}
	else {
		double t = im / re;
		return re * sqrt(1.0 + t * t);
	}
}

static complex_t cdiv(complex_t x, complex_t y) {

	/* complex division avoiding overflow */

	double zr, zi;
	double xr = x.r;
	double xi = x.i;
	double yr = y.r;
	double yi = y.i;

	if (yr == 0 && yi == 0) {
		zr = zi = DBL_MAX;
	}
	else if (fabs(yr) < fabs(yi)) {
		double u = yr / yi;
		double v = yi + u * yr;

		zr = (xr * u + xi) / v;
		zi = (xi * u - xr) / v;
	}
	else {
		double u = yi / yr;
		double v = yr + u * yi;

		zr = (xr + xi * u) / v;
		zi = (xi - xr * u) / v;
	}

	return complex(zr, zi);
}

/*-----------------------------------------------------------------------*/
static double cauchy_bound(uint32 n, const complex_t p[]) {

	/* computes a lower bound on the moduli of the zeros of a
	   polynomial p(x). norms(x) is a polynomial whose i_th coefficient
	   is the modulus of the i_th coeffcient of p(x) but
	   whose constant term is negative. The lower bound is the 
	   (unique) positive root x of norms(x) */

	double x, xmax, f, dx, df;
	double norms[MAX_ROOTFINDER_DEGREE + 1];
	uint32 i;

	for (i = 0; i < n; i++)
		norms[i] = cmod(p[i]);
	norms[i] = -cmod(p[i]);

	/* compute upper estimate of bound: assume all the
	   middle terms of norms(x) are zero */

	xmax = exp((log(-norms[n]) - log(norms[0])) / (double) n);

	/* if ignoring the nonlinear terms of norms(x) produces
	   a smaller root, use that instead */

	if (norms[n - 1] != 0.0) {
		x = -norms[n] / norms[n - 1];
		xmax = MIN(x, xmax);
	}

	/* chop the interval (0, x) until until x is about
	   to make norms(x) change sign */

	do {
		x = xmax;
		xmax = 0.1 * x;

		f = norms[0];
		for (i = 1; i <= n; i++)
			f = f * xmax + norms[i];
	} while (f > 0.0);

	/* do newton iteration until x converges to two decimal places */

	dx = x;
	while (fabs(dx / x) > 0.005) {
		df = 0;
		f = norms[0];
		for (i = 1; i <= n; i++) {
			df = df * x + f;
			f = f * x + norms[i];
		}
		dx = f / df;
		x -= dx;
	}

	return x;
}

/*-----------------------------------------------------------------------*/
static complex_t poly_val(uint32 n, complex_t s, 
			const complex_t p[], complex_t q[]) {

	/* evaluates a polynomial p at s by the horner 
	   recurrence, placing the partial sums in q and 
	   returning the computed value */

	complex_t pv;
	uint32 i;

	pv = q[0] = p[0];
	for (i = 1; i <= n; i++)
		pv = q[i] = cmac(pv, s, p[i]);

	return pv;
}

/*-----------------------------------------------------------------------*/
static void next_hpoly(complex_t correction, jt_t * w) {

	/* calculates the next shifted h polynomial */

	uint32 i;
	complex_t *poly_aux = w->poly_aux;
	complex_t *hpoly_aux = w->hpoly_aux;
	complex_t *hpoly = w->hpoly;

	if (w->hpoly_root_found == 0) {
		hpoly[0] = poly_aux[0];
		for (i = 1; i < w->degree; i++) {
			hpoly[i] = cmac(hpoly_aux[i - 1], correction, 
						poly_aux[i]);
		}
	}
	else {
		/* we are essentially at a root of h(x); remove 
		   it by deflating the polynomial. Calling code 
		   always expects h(x) to have the same degree, 
		   so the high-order coefficient becomes zero */

		hpoly[0] = czero;
		for (i = 1; i < w->degree; i++)
			hpoly[i] = hpoly_aux[i - 1];
	}
}

/*-----------------------------------------------------------------------*/
static complex_t next_correction(complex_t pval, complex_t curr_root,
				jt_t * w) {

	/* computes -pval / hpoly(curr_root)
	   sets flag to true if hval is essentially zero. */

	complex_t *hpoly = w->hpoly;
	complex_t *hpoly_aux = w->hpoly_aux;
	complex_t hval = poly_val(w->degree - 1, curr_root, hpoly, hpoly_aux);

	if (cmod(hval) <= 10.0 * DBL_EPSILON * cmod(hpoly[w->degree - 1])) {
		w->hpoly_root_found = 1;
		return czero;
	}
	else {
		w->hpoly_root_found = 0;
		return cdiv(cneg(pval), hval);
	}
}

/*-----------------------------------------------------------------------*/
#define STAGE3_ITER 10

static uint32 stage3(complex_t *root, jt_t *w) {

	/* carries out the third stage iteration,
	   returns 1 if iteration converges */

	double mp, ms, tp;
	uint32 i, j;
	complex_t pval;
	complex_t correction;
	complex_t curr_root = *root;
	complex_t *poly = w->poly;
	complex_t *poly_aux = w->poly_aux;

	for (i = 0; i < STAGE3_ITER; i++) {

		/* evaluate poly at current root value */

		pval = poly_val(w->degree, curr_root, poly, poly_aux);

		/* calculate bound on the error in evaluating the polynomial 
		   by the horner recurrence */

		mp = cmod(pval);
		ms = cmod(curr_root);
		tp = cmod(poly_aux[0]) * mre / (DBL_EPSILON + mre);
		for (j = 0; j <= w->degree; j++)
			tp = tp * ms + cmod(poly_aux[j]);
		tp = tp * (DBL_EPSILON + mre) - mp * mre;

		if (mp <= 20.0 * tp) {
			/* polynomial value is smaller in value 
			   than a bound on the error in evaluating p, 
			   terminate the iteration */
			*root = curr_root;
			return 1;
		}

		/* calculate next h polynomial */

		correction = next_correction(pval, curr_root, w);
		next_hpoly(correction, w);

		/* use the next h polynomial to calculate the next
		   root estimate, using the current root estimate */

		correction = next_correction(pval, curr_root, w);
		curr_root = cadd(curr_root, correction);
	}

	return 0;
}

/*-----------------------------------------------------------------------*/
static uint32 stage2(uint32 stage2_iter, complex_t *root, jt_t *w) {

	uint32 i;
	complex_t curr_root;
	complex_t correction;
	complex_t pval;

	/* calculate first correction */

	curr_root = *root;
	pval = poly_val(w->degree, curr_root, w->poly, w->poly_aux);
	correction = next_correction(pval, curr_root, w);

	for (i = 0; i < stage2_iter; i++) {

		/* compute next h polynomial and new correction;
		   note that the fixed-shift iteration never changes
		   the value of curr_root, only the h polynomial */

		next_hpoly(correction, w);
		correction = next_correction(pval, curr_root, w);

		if (w->hpoly_root_found == 1)
			break;
	}

	/* attempt stage 3 with the final h polynomial and
	   final correction */

	*root = cadd(curr_root, correction);
	return stage3(root, w);
}

/*-----------------------------------------------------------------------*/
#define STAGE1_ITER 5

static void stage1(uint32 n, complex_t p[], complex_t h[]) {

	uint32 i, j;

	/* the initial h polynomial is a scaled version of the
	   derivative of the input polynomial p(x) */

	for (i = 0; i < n; i++)
		h[i] = cscale((double) (n - i) / n, p[i]);

	/* compute a series of no-shift h polynomials */

	for (i = 0; i < STAGE1_ITER; i++) {

		/* if the constant term is essentially zero, 
		   shift the h coefficients */

		if (cmod(h[n-1]) <= 10.0 * DBL_EPSILON * cmod(p[n-1])) {

			for (j = n - 1; j; j--)
				h[j] = h[j-1];
			h[j] = czero;
		}
		else {

			complex_t tmp = cdiv(cneg(p[n]), h[n-1]);
			for (j = n - 1; j; j--)
				h[j] = cmac(h[j-1], tmp, p[j]);
			h[j] = p[0];
		}
	}
}

/*-----------------------------------------------------------------------*/
static int find_one_root(complex_t *root, jt_t *w) {

	uint32 i, j, k;
	double bound;
	complex_t hpoly_start[MAX_ROOTFINDER_DEGREE + 1];

	/* find linear roots immediately */

	if (w->degree <= 1) {
		*root = cdiv(cneg(w->poly[1]), w->poly[0]);
		return 1;
	}

	/* calculate a lower bound on the modulus of the zeros */

	bound = cauchy_bound(w->degree, w->poly);

	/* stage 1 sets up the initial h polynomial only */

	stage1(w->degree, w->poly, hpoly_start);

	/* try the fixed-shift sequence twice */

	for (i = 0; i < 2; i++) {

		/* inner loop to select a shift */

		for (j = 1; j < 10; j++) {

			/* start point is chosen with modulus 'bound'
			   and a pseudo-random angle. In practice
			   we don't want to repeat previous work,
			   so the starting angle is rotated a fixed 
			   amount (94 degrees) from the previous 
			   start point */

			w->angle = cmul(w->angle, w->angle_inc);
			*root = cscale(bound, w->angle);

			/* do the second stage, with a varying
			   number of iterations.
			   
			   Note that every starting point uses the same
			   h polynomial. This is a change from all other
			   cpoly() versions, including the original 1972 
			   fortran, which uses a global h array that is
			   not reinitialized when a new start point is
			   chosen (I think erroneously) */

			for (k = 0; k < w->degree; k++)
				w->hpoly[k] = hpoly_start[k];

			if (stage2(10 * j, root, w) == 1)
				return 1;
		}
	}

	return 0;
}

/*-----------------------------------------------------------------------*/
static uint32 jenkins_traub(complex_t poly[], 
			uint32 degree, complex_t roots[]) {

	/* main Jenkins-Traub driver; returns number 
	   of roots found */

	uint32 i; 
	uint32 roots_found;
	jt_t w;

	/* remove any zeros at the origin */

	for (i = degree, roots_found = 0; i; i--, roots_found++) {
		if (poly[i].r != 0.0 || poly[i].i != 0.0)
			break;

		roots[roots_found] = czero;
	}
	w.degree = i;

	/* initialize */

	for (i = 0; i <= w.degree; i++)
		w.poly[i] = poly[i];

	w.angle = complex(M_SQRT1_2, -M_SQRT1_2);
	w.angle_inc = complex(cos(94 * M_PI / 180), 
			   sin(94 * M_PI / 180));

	/* loop to find roots */

	for (; roots_found < degree; roots_found++) {

		if (find_one_root(roots + roots_found, &w) == 0)
			break;

		/* deflate the polynomial */

		(w.degree)--;
		for (i = 0; i <= w.degree; i++)
			w.poly[i] = w.poly_aux[i];
	}

	return roots_found;
}

/*------------------------------------------------------------------*/
#define NEWTON_ITER 10

static uint32 polish_root(dd_complex_t *poly, uint32 degree,
			dd_complex_t x, dd_complex_t *root,
			double eps) {

	uint32 i = 0;
	double eps2 = eps * eps;

	for (i = 0; i < NEWTON_ITER; i++) {

		uint32 j = degree;
		dd_complex_t f = poly[j];
		dd_complex_t df = cplx_set_d(0.0, 0.0);
		dd_complex_t dx;
		dd_t abs_x, abs_dx;

		for (j--; (int32)j >= 0; j--) {
			df = cplx_add(cplx_mul(df, x), f);
			f = cplx_add(cplx_mul(f, x), poly[j]);
		}
		dx = cplx_div(f, df);
		x = cplx_sub(x, dx);

		abs_x = dd_add_dd(dd_mul_dd(x.r, x.r),
				  dd_mul_dd(x.i, x.i));
		abs_dx = dd_add_dd(dd_mul_dd(dx.r, dx.r),
				  dd_mul_dd(dx.i, dx.i));

		if (dd_cmp_dd(abs_dx, dd_mul_d(abs_x, eps2)) <= 0)
			break;
	}

	*root = x;
	return 0;
}

/*------------------------------------------------------------------*/
uint32 find_poly_roots(dd_t *poly, uint32 degree, dd_complex_t *roots) {

	uint32 i;
	dd_complex_t ddcoeffs[MAX_ROOTFINDER_DEGREE + 1];
	complex_t dcoeffs[MAX_ROOTFINDER_DEGREE + 1];
	complex_t droots[MAX_ROOTFINDER_DEGREE + 1];

	if (degree == 1) {
		roots[0].r = dd_div_dd(dd_neg(poly[0]), poly[1]);
		roots[0].i = dd_set_d(0.0);
		return 0;
	}

	for (i = 0; i <= degree; i++) {
		ddcoeffs[i].r = poly[i];
		ddcoeffs[i].i = dd_set_d(0.0);
		dcoeffs[degree - i] = complex(poly[i].hi, 0.0);
	}

	/* find the roots to a relative error close to the
	   double-precision limit */

	if (jenkins_traub(dcoeffs, degree, droots) != degree)
		return 1;

	/* polish each root */

	for (i = 0; i < degree; i++) {

		if (polish_root(ddcoeffs, degree,
				cplx_set_d(droots[i].r, droots[i].i),
				roots + i, 1e-30) != 0)
			return 2;

		/* change roots with very small imaginary part to
		   be explicitly real roots */

		if (dd_cmp_dd(dd_fabs(roots[i].i),
			      dd_mul_d(dd_fabs(roots[i].r), 1e-30)) <= 0) {
			roots[i].i = dd_set_d(0.0);
		}
	}

	return 0;
}
