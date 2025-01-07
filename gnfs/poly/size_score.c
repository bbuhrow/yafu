/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: size_score.c 987 2015-07-11 16:13:18Z jasonp_sf $
--------------------------------------------------------------------*/

#include "poly.h"

#if MAX_POLY_DEGREE > 8
#error "polynomial degree > 8 not supported"
#endif

typedef struct {
	double power;
	ddpoly_t *rpoly;
	ddpoly_t *apoly;
} dparam_t;

typedef struct {
	dickman_t *dickman_aux;
	double root_score_r;
	double root_score_a;
	double skew_x;
	double skew_y;
	double rfb_limit;
	double afb_limit;
	double log_rfb_limit;
	double log_afb_limit;
	double sieve_area;
	dpoly_t rpoly;
	dpoly_t apoly;
} murphy_param_t;

static const double recip_factorial[] = {
  5.00000000000000000e-01,
  1.66666666666666657e-01,
  4.16666666666666644e-02,
  8.33333333333333322e-03,
  1.38888888888888894e-03,
  1.98412698412698413e-04,
  2.48015873015873016e-05,
  2.75573192239858925e-06,
  2.75573192239858883e-07,
  2.50521083854417202e-08,
  2.08767569878681002e-09,
  1.60590438368216133e-10,
  1.14707455977297245e-11,
  7.64716373181981641e-13,
  4.77947733238738525e-14,
  2.81145725434552060e-15,
};

#define NUM_RECIP_FACTORIALS (sizeof(recip_factorial) / sizeof(double))

/*------------------------------------------------------------------*/
static int compare_double(const void *x, const void *y) {
	double *xx = (double *)x;
	double *yy = (double *)y;

	if (*xx > *yy)
		return 1;
	if (*xx < *yy)
		return -1;
	return 0;
}

/*------------------------------------------------------------------*/
static double get_polyval(ddpoly_t *poly, double x, double h) {

	/* evaluate poly(x+h) using a taylor series centered
	   at poly(x). We compute poly(x) and the correction
	   from poly(x) separately, because in practice x and h
	   can have such different magnitudes that x+h could
	   equal x to machine precision. */

	double base; 
	double off;
	double hpow;
	dd_t *p = poly->coeff;
	double p0, p1, p2, p3, p4, p5, p6, p7, p8;

	switch (poly->degree) {
	case 0:
		p0 = p[0].hi;

		base = p0;
		off = 0;
		break;
	case 1:
		p0 = p[0].hi;
		p1 = p[1].hi;

		base = p1*x+p0;
		off = p1*h;
		break;
	case 2:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;

		base = (p2*x+p1)*x+p0;
		hpow = h;
		off = hpow * (2*p2*x+p1);
		hpow *= h;
		off += hpow * p2;
		break;
	case 3:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;

		base = ((p3*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * ((3*p3*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * (3*p3*x+p2);
		hpow *= h;
		off += hpow * p3;
		break;
	case 4:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;

		base = (((p4*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * (((4*p4*x+3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * ((6*p4*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * (4*p4*x+p3);
		hpow *= h;
		off += hpow * p4;
		break;
	case 5:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;

		base = ((((p5*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * ((((5*p5*x+4*p4)*x+3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * (((10*p5*x+6*p4)*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * ((10*p5*x+4*p4)*x+p3);
		hpow *= h;
		off += hpow * (5*p5*x+p4);
		hpow *= h;
		off += hpow * p5;
		break;
	case 6:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;
		p6 = p[6].hi;

		base = (((((p6*x+p5)*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * (((((6*p6*x+5*p5)*x+4*p4)*x+3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * ((((15*p6*x+10*p5)*x+6*p4)*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * (((20*p6*x+10*p5)*x+4*p4)*x+p3);
		hpow *= h;
		off += hpow * ((15*p6*x+5*p5)*x+p4);
		hpow *= h;
		off += hpow * (6*p6*x+p5);
		hpow *= h;
		off += hpow * p6;
		break;
 	case 7:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;
		p6 = p[6].hi;
		p7 = p[7].hi;

 		base = ((((((p7*x+p6)*x+p5)*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
 		hpow = h;
 		off = hpow * ((((((7*p7*x+6*p6)*x+5*p5)*x+4*p4)*x+
 						3*p3)*x+2*p2)*x+p1);
 		hpow *= h;
 		off += hpow * (((((21*p7*x+15*p6)*x+10*p5)*x+
 						6*p4)*x+3*p3)*x+p2);
 		hpow *= h;
 		off += hpow * ((((35*p7*x+20*p6)*x+10*p5)*x+4*p4)*x+p3);
 		hpow *= h;
 		off += hpow * (((35*p7*x+15*p6)*x+5*p5)*x+p4);
 		hpow *= h;
 		off += hpow * ((21*p7*x+6*p6)*x+p5);
 		hpow *= h;
 		off += hpow * (7*p7*x+p6);
 		hpow *= h;
 		off += hpow * p7;
 		break;
	case 8:
		p0 = p[0].hi;
		p1 = p[1].hi;
		p2 = p[2].hi;
		p3 = p[3].hi;
		p4 = p[4].hi;
		p5 = p[5].hi;
		p6 = p[6].hi;
		p7 = p[7].hi;
		p8 = p[8].hi;

		base = (((((((p8*x+p7)*x+p6)*x+p5)*x+p4)*x+p3)*x+p2)*x+p1)*x+p0;
		hpow = h;
		off = hpow * (((((((8*p8*x+7*p7)*x+6*p6)*x+5*p5)*x+4*p4)*x+
						3*p3)*x+2*p2)*x+p1);
		hpow *= h;
		off += hpow * ((((((28*p8*x+21*p7)*x+15*p6)*x+10*p5)*x+
						6*p4)*x+3*p3)*x+p2);
		hpow *= h;
		off += hpow * (((((56*p8*x+35*p7)*x+20*p6)*x+10*p5)*x+4*p4)*x+p3);
		hpow *= h;
		off += hpow * ((((70*p8*x+35*p7)*x+15*p6)*x+5*p5)*x+p4);
		hpow *= h;
		off += hpow * (((56*p8*x+21*p7)*x+6*p6)*x+p5);
		hpow *= h;
		off += hpow * ((28*p8*x+7*p7)*x+p6);
		hpow *= h;
		off += hpow * (8*p8*x+p7);
		hpow *= h;
		off += hpow * p8;
		break;
	default:
		base = off = 0;
		break;
	}
	return base + off;
}

/*------------------------------------------------------------------*/
static double integrand(double x, double h, void *params) {

	/* callback for numerical integration */

	dparam_t *aux = (dparam_t *)params;
	double polyval;

	polyval = get_polyval(aux->rpoly, x, h) *
		  get_polyval(aux->apoly, x, h);
	return pow(fabs(polyval), aux->power);
}

/*------------------------------------------------------------------*/
#define INTEGRATE_LIMIT 1e12

uint32 analyze_poly_size(integrate_t *integ_aux,
			ddpoly_t *rpoly, ddpoly_t *apoly, 
			double *result) {

	uint32 i, j;
	uint32 rdeg, adeg;
	dparam_t params;
	dd_complex_t roots[2*MAX_POLY_DEGREE];
	uint32 num_roots;
	double endpoints[2*MAX_POLY_DEGREE+2];
	uint32 num_endpoints;

	rdeg = rpoly->degree;
	adeg = apoly->degree;
	params.power = -2.0 / (rdeg + adeg);
	params.rpoly = rpoly;
	params.apoly = apoly;

	/* Rather than use Murphy's method to rate the average
	   size of polynomial values, we use a result of D.J.
	   Bernstein: For rational poly R(x) and algebraic poly
	   A(x), the number of polynomial values where x is rational
	   and the size of A(x)*R(x) achieves a given value is 
	   proportional to the following superelliptic integral:

	   	dx / ((R(x) * A(x))^2) ^ (1/(deg(R(x))+deg(A(x))))

	   for x from -infinity to +infinity. Larger values of
	   this integral imply more values that get found by
	   a siever. Bernstein writes that this estimate is 
	   "extremely accurate", but has not elaborated to date.
	   The integration variable x refers to the a/b values used
	   by a siever, and it is both more realistic and simpler to
	   make the integration interval finite but large 
	   
	   The integrand takes on a wide range of values depending
	   on how close x is to a root of R(x)*A(x), so we have to
	   partition the integral */


	/* find the roots of R(x)*A(x) */

	if (find_poly_roots(params.rpoly->coeff, rdeg, roots)) {
		printf("rational poly rootfinder failed\n");
		return 1;
	}
	num_roots = rdeg;
	if (find_poly_roots(params.apoly->coeff, adeg, roots + num_roots)) {
		printf("algebraic poly rootfinder failed\n");
		return 1;
	}
	num_roots += adeg;

	/* squeeze out roots whose real part would not lie
	   in the integration interval, along with all complex 
	   conjugates of complex roots.

	   We cannot remove all complex roots because if the
	   imaginary part of the root is small then the value
	   of A(x)*R(x) is 'almost' a singularity at this point,
	   and ignoring it would throw off the numerical integrator */
	   
	for (i = j = 0; i < num_roots; i++) {
		if (roots[i].i.hi < 0.0 || 
		    roots[i].r.hi <= -INTEGRATE_LIMIT ||
		    roots[i].r.hi >= INTEGRATE_LIMIT)
			continue;
		endpoints[j++] = roots[i].r.hi;
	}
	endpoints[j] = -INTEGRATE_LIMIT;
	endpoints[j+1] = INTEGRATE_LIMIT;
	num_endpoints = j + 2;

	/* sort the endpoints into increasing order */

	qsort(endpoints, (size_t)num_endpoints,
			sizeof(double), compare_double);

	/* integrate */

	*result = 0.0;

	integrate_run(integ_aux, integrand,
			&params, endpoints, num_endpoints);

	/* test for convergence */

	if (integ_aux->result > 1 ||
	    integ_aux->error > SIZE_EPS * fabs(integ_aux->result)) {
		printf("integrator failed %e %e\n", 
				integ_aux->error,
				SIZE_EPS * fabs(integ_aux->result));
		return 2;
	}
	*result = integ_aux->result;
	return 0;
}

/*------------------------------------------------------------------*/
static double eval_dpoly(dpoly_t *poly,
			double x0, double xh, 
			double y0, double yh) {

	/* evaluate poly(x0+xh, y0+yh) using a Taylor series
	   centered at poly(x0,y0). We have to do this because
	   the numerical integrator will sample the integrand 
	   extremely close to poly(x,y)=0, and in that case
	   xh and yh will be of such different magnitude compared
	   to x0 and y0 that we will run into numerical trouble */

	uint32 i, j, k;
	uint32 deg = poly->degree;
	double *coeff = poly->coeff;
	double xbinom[MAX_POLY_DEGREE+1][MAX_POLY_DEGREE+1];
	double ybinom[MAX_POLY_DEGREE+1][MAX_POLY_DEGREE+1];
	double sum[MAX_POLY_DEGREE+1];
	double *xrow, *yrow;
	double res;

	xbinom[1][0] = xh;
	xbinom[1][1] = x0;
	for (i = 2; i <= deg; i++) {

		double *polyj = xbinom[i-1];
		double *polyk = xbinom[1];
		double *row = xbinom[i];

		for (j = 0; j <= i; j++)
			row[j] = 0;

		for (j = 0; j < i; j++) {
			row[j+0] += polyj[j] * polyk[0];
			row[j+1] += polyj[j] * polyk[1];
		}
	}

	ybinom[1][0] = yh;
	ybinom[1][1] = y0;
	for (i = 2; i <= deg; i++) {

		double *polyj = ybinom[i-1];
		double *polyk = ybinom[1];
		double *row = ybinom[i];

		for (j = 0; j <= i; j++)
			row[j] = 0;

		for (j = 0; j < i; j++) {
			row[j+0] += polyj[j] * polyk[0];
			row[j+1] += polyj[j] * polyk[1];
		}
	}

	xrow = xbinom[deg];
	yrow = ybinom[deg];
	for (i = 0; i <= deg; i++)
		sum[i] = xrow[i] * coeff[deg] + yrow[i] * coeff[0];

	for (i = 1; i < deg; i++) {

		uint32 xdeg = i;
		uint32 ydeg = deg - i;

		xrow = xbinom[xdeg];
		yrow = ybinom[ydeg];
		for (j = 0; j <= xdeg; j++) {
			for (k = 0; k <= ydeg; k++) {
				sum[j+k] += coeff[i] * xrow[j] * yrow[k];
			}
		}
	}

	res = sum[0];
	for (i = 1; i <= deg; i++)
		res += sum[i];

	return res;
}

/*------------------------------------------------------------------*/
static double murphy_integrand(double r, double h, void *params) {

	murphy_param_t *p = (murphy_param_t *)params;
	double x0, xh, y0, yh;
	double polyval_r, polyval_a;

	/* we need to convert the angular measure (r+h), with
	   r and h of possibly very different magnitudes, into
	   cartesian measures (x0+xh, y0+yh) so that we can use
	   the homogeneous form of the NFS polynomials */

	if (fabs(h) > 0.1) {
		x0 = cos(r + h);
		y0 = sin(r + h);
		xh = yh = 0;
	}
	else {
		uint32 i;
		double term, hpow;

		/* sum the Taylor series for cos(r+h) and sin(r+h)
		   simultaneously */

		x0 = cos(r);
		y0 = sin(r);
		hpow = h;
		xh = -y0 * h;
		yh = x0 * h;

		for (i = 0; i < NUM_RECIP_FACTORIALS; i += 2) {

			hpow *= -h;
			term = hpow * recip_factorial[i];

			xh += x0 * term;
			yh += y0 * term;

			hpow *= h;
			term = hpow * recip_factorial[i+1];

			xh -= y0 * term;
			yh += x0 * term;

			if (fabs(term) <= 1e-16 * (fabs(xh) + fabs(yh)))
				break;
		}
	}

	/* skew the coordinates */

	x0 *= p->skew_x;
	xh *= p->skew_x;
	y0 *= p->skew_y;
	yh *= p->skew_y;

	/* evaluate the NFS polynomials at the skewed coordinates */

	polyval_r = eval_dpoly(&p->rpoly, x0, xh, y0, yh);
	polyval_a = eval_dpoly(&p->apoly, x0, xh, y0, yh);

	/* polyval_[ra] give a measure of the size of sieve values
	   contained in a ray emanating from the origin at an angle 
	   of (r+h) radians. The integrand is the probability that
	   sieve values on this ray are [ra]fb-limit-smooth, after 
	   small primes are removed */

	return dickman(p->dickman_aux,
			(log(fabs(polyval_r)) + p->root_score_r) / 
				p->log_rfb_limit) *
	       dickman(p->dickman_aux,
			(log(fabs(polyval_a)) + p->root_score_a) / 
				p->log_afb_limit);
}

/*------------------------------------------------------------------*/
uint32 analyze_poly_murphy(integrate_t *integ_aux, dickman_t *dickman_aux,
			ddpoly_t *rpoly, double root_score_r,
			ddpoly_t *apoly, double root_score_a,
			double skewness, double *result,
			uint32 *num_real_roots) {

	/* Given the skewness and root score for an NFS polynomial
	   pair, calculate the probability that an average sieve 
	   value in the sieving region has all rational (resp. algebraic)
	   factors less than rfb_limit (resp. afb_limit) 
	 
	   Ideally the sieving area and factor base limits should vary
	   with the size of the NFS input, but we fix them here to
	   be compatible with the code in pol51. That code uses a
	   trapezoidal approximation to the integral and computes
	   Dickman's function via linear interpolation from pre-tabulated
	   values. The following uses a full numerical integrator and 
	   the classical Dickman series instead, but the integrand is
	   so smooth most of the time that the end effect is the same.
	   This code uses about 90% fewer integrand evaluations though */

	uint32 i, j;
	murphy_param_t params;
	dd_complex_t roots[2*MAX_POLY_DEGREE];
	double angles[2*MAX_POLY_DEGREE+2];
	uint32 num_roots, num_angles;
	uint32 rdeg = rpoly->degree;
	uint32 adeg = apoly->degree;

	const double rfb_limit = 5000000;
	const double afb_limit = 10000000;
	const double sieve_area = 1e16;

	params.dickman_aux = dickman_aux;
	params.root_score_r = root_score_r;
	params.root_score_a = root_score_a;
	params.rfb_limit = rfb_limit;
	params.afb_limit = afb_limit;
	params.log_rfb_limit = log(rfb_limit);
	params.log_afb_limit = log(afb_limit);
	params.skew_x = sqrt(sieve_area * skewness);
	params.skew_y = params.skew_x / skewness;

	params.rpoly.degree = rdeg;
	for (i = 0; i <= rdeg; i++)
		params.rpoly.coeff[i] = rpoly->coeff[i].hi;

	params.apoly.degree = adeg;
	for (i = 0; i <= adeg; i++)
		params.apoly.coeff[i] = apoly->coeff[i].hi;

	/* find the roots of rpoly * apoly */

	if (find_poly_roots(rpoly->coeff, rdeg, roots)) {
		printf("rational poly rootfinder failed\n");
		return 1;
	}
	num_roots = rdeg;
	if (find_poly_roots(apoly->coeff, adeg, roots + num_roots)) {
		printf("algebraic poly rootfinder failed\n");
		return 1;
	}
	num_roots += adeg;

	*num_real_roots = 0;
	for (i = 0; i < adeg; i++) {
		if (dd_cmp_d(roots[rdeg+i].i, 0.0) == 0)
			(*num_real_roots)++;
	}

	/* convert the roots to angles between 0 and pi. Since the
	   integrator will skew values of x and y derived from these
	   angles, skew the roots the other way to compensate */

	for (i = j = 0; i < num_roots; i++) {

		if (roots[i].i.hi < 0.0)
			continue;

		angles[j++] = atan2(1.0, roots[i].r.hi / skewness);
	}
	angles[j] = 0.0;
	angles[j+1] = M_PI;
	num_angles = j + 2;

	/* sort in order of increasing angle */

	qsort(angles, (size_t)num_angles, sizeof(double), compare_double);

	integrate_run(integ_aux, murphy_integrand, 
			&params, angles, num_angles);

	*result = integ_aux->result / M_PI;
	return 0;
}
