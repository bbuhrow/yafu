/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: dd.h 849 2013-03-09 08:02:27Z brgladman $
--------------------------------------------------------------------*/

#ifndef _DD_H_
#define _DD_H_

/* Basic arithmetic on double-double operands.
   See the QD package of D.H. Bailey and the
   doubledouble package of Keith Briggs. */

#include <util.h>
#include <mp.h>
#include <ms_gmp_xface.h>

#ifdef __cplusplus
extern "C" {
#endif

   /* These routines *require* IEEE 53-bit double precision,
      even on x86 processors that support higher precision */

#if defined(WIN32) && !defined(_WIN64) 
	#include <float.h>
	typedef uint32 dd_precision_t;
#elif (defined(__GNUC__) || defined(__ICL)) && \
		(defined(__i386__) || defined(__x86_64__))
	#include <float.h>
	typedef uint16 dd_precision_t;
#else
	typedef uint32 dd_precision_t;
#endif

static INLINE dd_precision_t dd_set_precision_ieee(void) {
#if defined(WIN32) && !defined(_WIN64)
	dd_precision_t old_prec = _control87(0, 0);
	_control87(_PC_53, _MCW_PC);
	return old_prec;

#elif defined(GCC_ASM32X) || defined(GCC_ASM64X)
	dd_precision_t old_prec, new_prec;
	ASM_G volatile ("fnstcw %0":"=m"(old_prec));
	new_prec = (old_prec & ~0x0300) | 0x0200;
	ASM_G volatile ("fldcw %0": :"m"(new_prec));
	return old_prec;
	
#else
	return 0;
#endif
}

static INLINE void dd_clear_precision(dd_precision_t old_prec) {
#if defined(WIN32) && !defined(_WIN64)
	_control87(old_prec, 0xffffffff);
#elif defined(GCC_ASM32X) || defined(GCC_ASM64X)
	ASM_G volatile ("fldcw %0": :"m"(old_prec));
#endif
}

static INLINE uint32 dd_precision_is_ieee(void) {
#if defined(WIN32) && !defined(_WIN64)
	dd_precision_t prec = _control87(0, 0);
	return  ((prec & _MCW_PC) == _PC_53) ? 1 : 0;

#elif defined(GCC_ASM32X) || defined(GCC_ASM64X)
	dd_precision_t prec;
	ASM_G volatile ("fnstcw %0":"=m"(prec));
	return ((prec & ~0x0300) == 0x0200) ? 1 : 0;
#else
	return 1;
#endif
}

/* basic structure */

typedef struct {
	double hi;
	double lo;
} dd_t;

/* assignment */

static INLINE dd_t dd_set_d(double x) {
	dd_t dd;
	dd.hi = x;
	dd.lo = 0.0;
	return dd;
}

static INLINE dd_t dd_set_dd(double hi, double lo) {
	dd_t dd;
	dd.hi = hi;
	dd.lo = lo;
	return dd;
}

/* addition */

#define QUICK_TWO_SUM(a, b, sum, err) {		\
	double u = (a) + (b);			\
	(err) = (b) - (u - (a));		\
	(sum) = u;				\
}

#define TWO_SUM(a, b, sum, err) { 		\
	double u = (a) + (b);			\
	double bb = u - (a);			\
	(sum) = u;				\
	(err) = ((a) - (u - bb)) + ((b) - bb);	\
}

static INLINE dd_t dd_add_d(dd_t a, double b) {
	
	double s1, s2;

	TWO_SUM(a.hi, b, s1, s2);
	s2 += a.lo;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	return dd_set_dd(s1, s2);
}

static INLINE dd_t dd_add_dd(dd_t a, dd_t b) {

	double s1, s2, t1, t2;

	TWO_SUM(a.hi, b.hi, s1, s2);
	TWO_SUM(a.lo, b.lo, t1, t2);
	s2 += t1;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	s2 += t2;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	return dd_set_dd(s1, s2);
}

/* subtraction */

#define TWO_DIFF(a, b, sum, err) { 		\
	double u = (a) - (b);			\
	double bb = u - (a);			\
	(sum) = u;				\
	(err) = ((a) - (u - bb)) - ((b) + bb);	\
}

static INLINE dd_t dd_sub_d(dd_t a, double b) {
	
	double s1, s2;

	TWO_DIFF(a.hi, b, s1, s2);
	s2 += a.lo;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	return dd_set_dd(s1, s2);
}

static INLINE dd_t dd_sub_dd(dd_t a, dd_t b) {

	double s1, s2, t1, t2;

	TWO_DIFF(a.hi, b.hi, s1, s2);
	TWO_DIFF(a.lo, b.lo, t1, t2);
	s2 += t1;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	s2 += t2;
	QUICK_TWO_SUM(s1, s2, s1, s2);
	return dd_set_dd(s1, s2);
}

/* multiplication */

#if defined(__GNUC__) && (defined(__powerpc__) || defined(__ppc64__))

static INLINE double __fmsub (double a, double c, double b) {
	double result;
	ASM_G ("fmsub %0, %1, %2, %3"
	    : "=f" (result) 
	    : "f" (a), "f" (c), "f" (b));
	return result;
}

#define TWO_PROD(a, b, prod, err) {	\
	double u = (a) * (b);		\
	(prod) = u;			\
	(err) = __fmsub((a), (b), u);	\
}

#else

#define SPLIT_CONST 134217729.0               /* = 2^27 + 1 */

#define SPLIT(a, hi, lo) {		\
	double u = SPLIT_CONST * (a);\
	(hi) = u - (u - (a));	\
	(lo) = (a) - (hi);		\
}

#define TWO_PROD(a, b, prod, err) {		\
	double a_hi, a_lo, b_hi, b_lo;		\
	double u = (a) * (b);			\
	SPLIT((a), a_hi, a_lo);			\
	SPLIT((b), b_hi, b_lo);			\
	(prod) = u;				\
	(err) = ((a_hi * b_hi - u) + 		\
		a_hi * b_lo + a_lo * b_hi) + 	\
		a_lo * b_lo;			\
}

#endif


static INLINE dd_t dd_mul_d(dd_t a, double b) {

	double p1, p2;

	TWO_PROD(a.hi, b, p1, p2);
	p2 += (a.lo * b);
	QUICK_TWO_SUM(p1, p2, p1, p2);
	return dd_set_dd(p1, p2);
}

static INLINE dd_t dd_mul_dxd(double a, double b) {

	double p, e;
	TWO_PROD(a, b, p, e);
	return dd_set_dd(p, e);
}

static INLINE dd_t dd_mul_dpow2(dd_t a, double b) {

	return dd_set_dd(a.hi * b, a.lo * b);
}

static INLINE dd_t dd_mul_dd(dd_t a, dd_t b) {

	double p1, p2;

	TWO_PROD(a.hi, b.hi, p1, p2);
	p2 += (a.hi * b.lo);
	p2 += (a.lo * b.hi);
	QUICK_TWO_SUM(p1, p2, p1, p2);
	return dd_set_dd(p1, p2);
}

/* division */

static INLINE dd_t dd_div_dxd(double a, double b) {

	double q1, q2, p1, p2;
	double s, e;

	q1 = a / b;

	/* Compute  a - q1 * b */
	TWO_PROD(q1, b, p1, p2);
	TWO_DIFF(a, p1, s, e);
	e -= p2;

	/* get next approximation */
	q2 = (s + e) / b;

	/* renormalize */
	QUICK_TWO_SUM(q1, q2, s, e);
	return dd_set_dd(s, e);
}

static INLINE dd_t dd_div_d(dd_t a, double b) {

	double q1, q2, p1, p2;
	double s, e;
	dd_t r;

	q1 = a.hi / b;  /* approximate quotient */

	/* Compute  this - q1 * b */
	TWO_PROD(q1, b, p1, p2);
	TWO_DIFF(a.hi, p1, s, e);
	e += a.lo;
	e -= p2;

	/* get next approximation */
	q2 = (s + e) / b;

	/* renormalize */
	QUICK_TWO_SUM(q1, q2, r.hi, r.lo);
	return r;
}

static INLINE dd_t dd_div_dd(dd_t a, dd_t b) {

	double q1, q2, q3;
	dd_t r;

	q1 = a.hi / b.hi;
	r = dd_sub_dd(a, dd_mul_d(b, q1));

	q2 = r.hi / b.hi;
	r = dd_sub_dd(r, dd_mul_d(b, q2));

	q3 = r.hi / b.hi;
	QUICK_TWO_SUM(q1, q2, q1, q2);
	return dd_add_d(dd_set_dd(q1, q2), q3);
}

/* miscellaneous functions */

static INLINE dd_t dd_neg(dd_t a) {

	return dd_set_dd(-a.hi, -a.lo);
}

static INLINE int32 dd_cmp_d(dd_t a, double b) {
	if (a.hi < b)
		return -1;
	if (a.hi > b)
		return 1;

	if (a.lo < 0)
		return -1;
	if (a.lo > 0)
		return 1;
	return 0;
}

static INLINE int32 dd_cmp_dd(dd_t a, dd_t b) {
	if (a.hi < b.hi)
		return -1;
	if (a.hi > b.hi)
		return 1;

	if (a.lo < b.lo)
		return -1;
	if (a.lo > b.lo)
		return 1;
	return 0;
}

static INLINE dd_t dd_fabs(dd_t a) {
	if (a.hi < 0 || (a.hi == 0 && a.lo < 0))
		return dd_set_dd(-a.hi, -a.lo);
	else
		return a;
}

static INLINE dd_t dd_mp2dd(mp_t *x) {

	/* convert a multiple-precision number to 
	   a double-double */

	int32 i;
	dd_t d;

	if (mp_is_zero(x))
		return dd_set_d(0.0);
	
	i = x->nwords - 1;
	d = dd_set_d((double)x->val[i]);
	for (i--; i >= 0; i--) {
		d = dd_mul_dpow2(d, MP_RADIX);
		d = dd_add_d(d, (double)x->val[i]);
	}

	return d;
}

static INLINE dd_t dd_signed_mp2dd(signed_mp_t *x) {

	if (x->sign == NEGATIVE)
		return dd_neg(dd_mp2dd(&x->num));
	else
		return dd_mp2dd(&x->num);
}

static INLINE dd_t dd_gmp2dd(mpz_t x) {

	signed_mp_t mpx;

	gmp2mp(x, &mpx.num);

	mpx.sign = POSITIVE;
	if (mpz_sgn(x) < 0)
		mpx.sign = NEGATIVE;

	return dd_signed_mp2dd(&mpx);
}

static INLINE void dd_dd2mp(dd_t d, mp_t *x) {

	/* convert a double-double to multiple-precison
	   format. If both halves of d are integers then
	   the conversion is exact */

	mp_t x1, x2;

	mp_d2mp(&d.hi, &x1);
	if (d.lo < 0) {
		double dneg = -d.lo;
		mp_d2mp(&dneg, &x2);
		mp_sub(&x1, &x2, x);
	}
	else {
		mp_d2mp(&d.lo, &x2);
		mp_add(&x1, &x2, x);
	}
}

#ifdef __cplusplus
}
#endif

#endif /* _DD_H_ */
