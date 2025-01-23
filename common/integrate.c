/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: integrate.c 23 2009-07-20 02:59:07Z jasonp_sf $
--------------------------------------------------------------------*/

#include <integrate.h>

/* This code implements numerical integration via double-exponential (DE)
   transformation. Basically the integrand is evaluated at a set of
   points that approach the limits of the integration interval faster
   than exponentially, with weights that tend to zero faster than
   exponentially. The result is that even if the integrand has singularities
   at both endpoints, as long as it is analytic and doesn't vary by 
   much in the interior the interval, this method converges to the numerical
   value of the integral with comparatively very few function evaluations.
   Basically we use the trapezoid rule in transformed coordinates and
   reuse previous computations when cutting the stepsize in half. The
   code computes the integral value, the residual, and an error estimate
   at every step. The DE transformation is comparatively unknown (I used
   to be very interested in numerical integration and had never heard of
   it before Brian Gladman brought it to my attention), and is applicable
   to a wide variety of different integration-type problems. For further 
   explanations, history and a bunch of references see M. Mori, 'The 
   Discovery of the Double Exponential Transformation and its Developments'.

   The code here is based on code from 1996 by Takuya Ooura, which was a
   C port from the original painfully obscure spaghetti fortran. I've cleaned 
   it up a lot and documented how I think things work, but it's always 
   possible I misinterpreted what the original code was doing */


/* bound on the amount of memory allocated */

#define MAX_TRAP_LEVELS 7

/* macro for evaluating the integrand. base is always
   the left or right integration endpoint, and the integrand
   is evaluated at (base+off). The two must be provided separately,
   because base and off can be of such different magnitudes that
   adding them together could cause numerical instability. This
   happens when base+off == base to the limit of double precision.
   When the endpoints are singularities it is important to calculate
   the contribution of off separately */

#define FUNC_EVAL(base, off) func(base, off, params)

/*---------- double exponential integrator, definite interval ------------*/

/* describes the integration of one subinterval */

typedef struct {
	double left;          /* left endpoint */
	double right;         /* right endpoint */
	double result;        /* numerical value of the integral */
	double error;         /* estimated error */
	uint32 level;
} subinterval_t;

/* describes one point in the interval to evaluate */

typedef struct {
	double abscissa;        /* distance from the left or right endpoint
				   where the integrand is evaluated */
	double weight;          /* weight applied to the integrand at
				   'abscissa' before adding to integral value */
	double residual_weight; /* weight applied to the integrand at
				   'abscissa' before adding to residual */
} de_point_t;

/* structure controlling numerical integration. Note that introductory
   numerical methods describe the trapezoid rule as having equal-size
   panels, with the rule at half the step size containing twice as many
   half-size panels. The DE transform is different: we compute the
   integral from -1 to 1 via a transformation to an integral from -infinity
   to +infinity, then evaluate the latter integral using a series of
   trapezoid panels that increase super-exponentially in width. The 'step-size'
   refers to the distance of the first sample point from the origin, and
   halving the step-size here means doubling the number of trapezoid panels. 
   While there are really an infinite number of trapezoid panels to sum, 
   we truncate the sum when the size of the abscissas become too large to 
   be represented accurately. Because the abscissas increase in size so 
   quickly, it only takes a few samples to reach this point.  The de_point_t 
   structures store the value of the abscissa after conversion back to the 
   (-1,1) interval, as offsets from +-1 */

typedef struct {
	uint32 num_points;         /* maximum number of trapezoid panels
				      allowed, beyond which the abscissas
				      are of size close to the limit of double 
				      precision */
	double relative_error;     /* the allowed relative error */
	de_point_t center_point;   /* initial sample point of the interval */
	de_point_t *initial_rule;  /* further trapezoid sample points at the 
				      coarsest step size (1.0). There are 
				      num_points of these, and since the rule 
				      is symmetric about the origin we have
				      up to 2*num_points samples at this
				      resolution */
	de_point_t *refined_rules; /* list of groups of num_points samples
				      that implement the trapezoid rule with
				      successively halved step sizes. For
				      level i, 0 <= i < num_trap_levels, the
				      step size (i.e. the distance of the
				      first abscissa to the origin) is
				      2^-(i+1) and there are 2^i groups of
				      up to num_points samples to add. Sample
				      groups are stored in order of increasing 
				      i, then increasing value of initial
				      abscissa. Each group of samples is also
				      symmetric about the origin.  */
	subinterval_t *heap;       /* used to store integrals of subintervals */
	uint32 num_heap;           /* number of subintervals */
	uint32 num_heap_alloc;     /* memory allocated for subintervals */
} de_t;

/*-----------------------------------------------------------------------*/
static void de_fill_sample_group(de_point_t *rule, 
			uint32 num_points, double first_sample, 
			double abscissa_range,
			double expo_large, double expo_small) {

	/* fill the information for the num_points trapezoid
	   panels, with the first abscissa at (first_sample * 
	   abscissa_range), 0 < first_sample <= 1. The sample 
	   points march off to infinity in transformed coordinates;
	   they are stored in untransformed coordinates, as offsets 
	   from +-abscissa_range, converging super-exponentially 
	   quickly to zero */

	uint32 i;
	double curr_expo = exp(first_sample * abscissa_range);
	double curr_expo_large = M_PI_2 * curr_expo;
	double curr_expo_small = M_PI_2 / curr_expo;

	for (i = 0; i < num_points; i++) {
		de_point_t *curr_point = rule + i;
		double abscissa = 1.0 / (1.0 + exp(curr_expo_large -
						   curr_expo_small));
		double weight = abscissa * (1.0 - abscissa) * abscissa_range;

		curr_point->abscissa = abscissa;
		curr_point->weight = weight * (curr_expo_large +
						curr_expo_small);
		curr_point->residual_weight = 4.0 * weight;

		curr_expo_large *= expo_large;
		curr_expo_small *= expo_small;
	}
}

/*-----------------------------------------------------------------------*/
static void de_init(integrate_t *aux, double relative_error) {

	/* initialize for all future integrations using aux */

	/* parameters of the algorithm: 
	   
	   min_val is the smallest function value we expect to
	   encounter

	   abscissa_scale expands the DE region from (-1,1) to 
	   (-abscissa_scale,abscissa_scale), presumably to avoid 
	   truncation error */

	const double min_val = 1e-30;
	const double abscissa_scale = 8.5;

	uint32 i, j;
	de_t *integ;
	de_point_t *curr_samples;
	uint32 num_points;
	double log_min_val = -log(min_val);
	double log_relative_error = 1.0 - log(relative_error);
	double abscissa_range = abscissa_scale / log_relative_error;
	double expo_large = exp(abscissa_range);
	double expo_small = 1.0 / expo_large;

	/* fill error bounds */

	integ = (de_t *)xmalloc(sizeof(de_t));
	integ->relative_error = relative_error;

	/* describe the initial sample point, in the middle 
	   of the interval */

	integ->center_point.weight = M_PI_2 * abscissa_range * 0.5;
	integ->center_point.residual_weight = abscissa_range;

	/* figure out the maximum number of samples needed */

	num_points = integ->num_points = (uint32)(log(log_min_val / M_PI_2) /
						 abscissa_range) + 1;

	curr_samples = integ->initial_rule = (de_point_t *)xmalloc(
					(num_points << MAX_TRAP_LEVELS) *
					sizeof(de_point_t));

	/* create the initial set of trapezoid sample points, with
	   step size of 1.0 * abscissa_range */

	de_fill_sample_group(curr_samples, num_points, 1.0,
				abscissa_range, expo_large, expo_small);

	/* now create trapezoid sample points in order of increasing
	   level. Level i has 2^i groups of samples; the j_th group
	   contains samples that start a distance of 
	   abscissa_range*j*2^-(i+1) from the origin and march away 
	   to inifinity in transformed coordinates at a given rate
	   (identical for all groups). Combining all the samples 
	   at level i, along with initial_rule and all samples at 
	   levels < i, produces the complete set of samples needed 
	   to implement the trapezoid rule with step size 2^-(i+1) */

	curr_samples += num_points;
	integ->refined_rules = curr_samples;

	for (i = 0; i < MAX_TRAP_LEVELS; i++) {
		uint32 num_groups = 1 << i;
		double step_size = 1.0 / num_groups;
		double sample_val = 0.5 * step_size;

		for (j = 0; j < num_groups; j++) {
			de_fill_sample_group(curr_samples, num_points, 
					sample_val, abscissa_range, 
					expo_large, expo_small);
			sample_val += step_size;
			curr_samples += num_points;
		}
	}

	/* set up the heap of subintervals */

	integ->num_heap = 0;
	integ->num_heap_alloc = 100;
	integ->heap = (subinterval_t *)xmalloc(integ->num_heap_alloc *
						sizeof(subinterval_t));

	aux->internal = (void *)integ;
}

/*-----------------------------------------------------------------------*/
static void de_free(integrate_t *aux) {

	de_t *i = (de_t *)(aux->internal);
	free(i->initial_rule);
	free(i->heap);
	free(i);
	aux->internal = NULL;
}

/*-----------------------------------------------------------------------*/
static void de_run_core(integrate_t *aux, integrand_t func, 
			void *params, subinterval_t *range) {

	/* compute the integral of func(x, params) in the
	   range defined by 'range' */

	uint32 i, j, k;
	uint32 num_groups;
	de_t *integ = (de_t *)aux->internal;
	de_point_t *points;
	uint32 num_points = integ->num_points;
	double lower_limit = range->left;
	double upper_limit = range->right;
	double interval = upper_limit - lower_limit;
	double step_size = 1.0;
	double result;
	double residual;
	double curr_error;
	double target_error;
	double left_val, right_val;

	/* evaluate func at the middle of the interval */

	result = FUNC_EVAL(0.5 * (lower_limit + upper_limit), 0.0);
	residual = result * integ->center_point.residual_weight;
	result *= integ->center_point.weight;

	/* compute the trapezoid rule approximation at the coarsest
	   step size. Sample at abscissa values that are symmetric 
	   about the origin */

	points = integ->initial_rule;

	for (i = 0; i < num_points; i++) {
		de_point_t *p = points + i;
		double abscissa = interval * p->abscissa;
		left_val = FUNC_EVAL(lower_limit, abscissa);
		right_val = FUNC_EVAL(upper_limit, -abscissa);
		result += (left_val + right_val) * p->weight;
		residual += (left_val + right_val) * p->residual_weight;
	}

	/* now compute trapezoid rules corresponding to successively
	   halved step sizes, until the estimate of the error indicates
	   convergence to the desired error tolerance */

	points = integ->refined_rules;
	i = 0;
	num_groups = 1;

	do {
		double old_result = result;
		double old_residual = residual;

		/* compute the trapezoid rule with step size 2^-(i+1) by
		   adding in 2^i groups of additional sample points.
		   Each sample falls somewhere in between two previous
	           samples, and all groups march from near the middle
		   of the interval out to the integration endpoints */

		for (j = 0; j < num_groups; j++, points += num_points) {

			/* update the integral value  and residual */

			for (k = 0; k < num_points; k++) {
				de_point_t *p = points + k;
				double abscissa = interval * p->abscissa;
				left_val = FUNC_EVAL(lower_limit, abscissa);
				right_val = FUNC_EVAL(upper_limit, -abscissa);
				result += (left_val + right_val) * p->weight;
				residual += (left_val + right_val) * 
							p->residual_weight;
			}
		}

		/* use the computed integral value and residual, at
		   the two different step sizes, to update the error 
		   estimate */

		curr_error = fabs(result - 2.0 * old_result) +
			     fabs(residual - 2.0 * old_residual);

		target_error = 0.1 * integ->relative_error * fabs(result);
		step_size *= 0.5;
		num_groups *= 2;

	} while (++i < MAX_TRAP_LEVELS && curr_error > target_error);

	/* compute the final integral value */

	range->result = result * interval * step_size;
	range->error = 2 * curr_error * interval;

#if 0
	printf("done: %le %le res %le err %le lvl %u\n", 
			lower_limit, upper_limit, 
			range->result, range->error, i-1);
#endif
}

/*-----------------------------------------------------------------------*/
/* boilerplate code for managing heaps */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void heapify(subinterval_t *h, uint32 index, uint32 size) {

	uint32 c;
	subinterval_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].error < h[c+1].error)
			c++;

		if (h[index].error < h[c].error) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].error < h[c].error) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void heap_insert(subinterval_t *h, uint32 size) {
	uint32 i = size - 1;
	subinterval_t tmp = h[i];

	while ((i > 0) && (h[HEAP_PARENT(i)].error < tmp.error) ) {
		h[i] = h[HEAP_PARENT(i)];
		i = HEAP_PARENT(i);
	}

	h[i] = tmp;
}

static void make_heap(subinterval_t *h, uint32 size) {

	int32 i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify(h, (uint32)i, size);
}

/*-----------------------------------------------------------------------*/
#define MAX_LEVELS 30

static uint32 de_run(integrate_t *aux, integrand_t func, 
			void *params, double *endpoints,
			uint32 num_endpoints) {

	uint32 i;
	de_t *integ = (de_t *)aux->internal;
	double result = 0;
	double error = 0;
	double rel_err = integ->relative_error;

	if (num_endpoints < 2)
		return 1;

	if (num_endpoints >= integ->num_heap_alloc) {
		integ->num_heap_alloc = MAX(num_endpoints + 100,
						2 * integ->num_heap_alloc);
		integ->heap = (subinterval_t *)xrealloc(integ->heap,
						integ->num_heap_alloc *
						sizeof(subinterval_t));
	}

	for (i = 0; i < num_endpoints - 1; i++) {
		subinterval_t *s = integ->heap + i;

		s->left = endpoints[i];
		s->right = endpoints[i+1];
		s->level = 0;
		de_run_core(aux, func, params, s);
		result += s->result;
		error += s->error;
	}
	integ->num_heap = i;

	if (error < rel_err * fabs(result)) {
		aux->result = result;
		aux->error = error;
		return 0;
	}

	make_heap(integ->heap, integ->num_heap);

	while (error > rel_err * fabs(result)) {

		subinterval_t *s = integ->heap;
		double left = s->left;
		double right = s->right;
		uint32 level = s->level;

		if (level > MAX_LEVELS)
			break;

		result -= s->result;
		error -= s->error;
		s->left = left;
		s->right = 0.5 * (left + right);
		s->level = level + 1;
		de_run_core(aux, func, params, s);
		result += s->result;
		error += s->error;
		heapify(integ->heap, 0, integ->num_heap);

		if (integ->num_heap == integ->num_heap_alloc) {
			integ->num_heap_alloc *= 2;
			integ->heap = (subinterval_t *)xrealloc(integ->heap,
						integ->num_heap_alloc *
						sizeof(subinterval_t));
		}

		s = integ->heap + integ->num_heap++;
		s->left = 0.5 * (left + right);
		s->right = right;
		s->level = level + 1;
		de_run_core(aux, func, params, s);
		result += s->result;
		error += s->error;
		heap_insert(integ->heap, integ->num_heap);
	}

	aux->result = result;
	aux->error = error;
	return 0;
}

/*------------------------- external interface --------------------------*/

void integrate_init(integrate_t *aux, 
			double relative_error,
			enum integrator_type type) {

	aux->type = type;
	if (type == double_exponential)
		de_init(aux, relative_error);
}

void integrate_free(integrate_t *aux) {

	if (aux->type == double_exponential)
		de_free(aux);
}

uint32 integrate_run(integrate_t *aux, 
			integrand_t func, void *params,
			double *endpoints, uint32 num_endpoints) {

	aux->result = 0.0;
	aux->error = 0.0;
	if (aux->type == double_exponential)
		return de_run(aux, func, params, 
				endpoints, num_endpoints);
	return 1;
}
