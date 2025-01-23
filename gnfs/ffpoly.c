/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.

$Id: ffpoly.c 967 2014-06-26 02:48:22Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>
#include "gnfs.h"

#if MAX_POLY_DEGREE > 8
#error "factor base generation requires poly degree 8 or less"
#endif

/* representation of polynomials with finite-field coefficients */

typedef struct {
	uint32 coef[2 * MAX_POLY_DEGREE + 1];
	uint32 degree;
} _poly_t;
typedef _poly_t poly_t[1];

/*------------------------------------------------------------------*/
static void poly_cp(poly_t dest, poly_t src) {

	dest[0] = src[0];
}

/*------------------------------------------------------------------*/
static void poly_fix_degree(poly_t op) { 

	int32 i = op->degree;

	while ((i > 0) && (op->coef[i] == 0))
		i--;
	op->degree = i;
}

/*------------------------------------------------------------------*/
static void poly_make_monic(poly_t res, poly_t a, uint32 p) {

	uint32 i;
	uint32 d = a->degree;
	uint32 msw = a->coef[d];

	if (msw != 1) {
		msw = mp_modinv_1(msw, p);
		res->degree = d;
		res->coef[d] = 1;
		for (i = 0; i < d; i++)
			res->coef[i] = mp_modmul_1(msw, a->coef[i], p);
	}
	else {
		poly_cp(res, a);
	}
}

/*------------------------------------------------------------------*/
static void poly_add(poly_t res, poly_t a, poly_t b, uint32 p) {

	uint32 i;

	if (a->degree < b->degree) {
		for (i = 0; i <= a->degree; i++)
			res->coef[i] = mp_modadd_1(a->coef[i], b->coef[i], p);
		for (; i <= b->degree; i++)
			res->coef[i] = b->coef[i];
		res->degree = b->degree;
	}
	else {
		for (i = 0; i <= b->degree; i++)
			res->coef[i] = mp_modadd_1(a->coef[i], b->coef[i], p);
		for (; i <= a->degree; i++)
			res->coef[i] = a->coef[i];
		res->degree = a->degree;
	}
	poly_fix_degree(res);
}

/*------------------------------------------------------------------*/
static void poly_mod(poly_t res, poly_t op, poly_t _mod, uint32 p) { 

	/* divide the polynomial 'op' by the polynomial '_mod'
	   and write the remainder to 'res'. All polynomial
	   coefficients are reduced modulo 'p' */

	int32 i;
	uint32 msw;
	poly_t tmp, mod;

	if (_mod->degree == 0) {
		memset(res, 0, sizeof(res[0]));
		return;
	}
	poly_cp(tmp, op);
	poly_make_monic(mod, _mod, p);

	while (tmp->degree >= mod->degree) {

		/* tmp <-- tmp - msw * mod * x^{deg(tmp)- deg(mod)} */

		msw = tmp->coef[tmp->degree];

		tmp->coef[tmp->degree] = 0;
		for (i = mod->degree-1; i >= 0; i--) {
			uint32 c = mp_modmul_1(msw, mod->coef[i], p);
			uint32 j = tmp->degree - (mod->degree - i);
			tmp->coef[j] = mp_modsub_1(tmp->coef[j], c, p);
		}
		poly_fix_degree(tmp);
	}
	poly_cp(res, tmp);
}

/*------------------------------------------------------------------*/
static void poly_modmul(poly_t res, poly_t a, poly_t b, 
			poly_t mod, uint32 p) { 

	uint32 i, j;
	poly_t prod;

	for (i = 0; i <= a->degree; i++)
		prod->coef[i] = mp_modmul_1(a->coef[i], b->coef[0], p);

	for (i = 1; i <= b->degree; i++) {
		for (j = 0; j < a->degree; j++) {
			uint32 c = mp_modmul_1(a->coef[j], b->coef[i], p);
			prod->coef[i+j] = mp_modadd_1(prod->coef[i+j], c, p);
		}
		prod->coef[i+j] = mp_modmul_1(a->coef[j], b->coef[i], p);
	}

	prod->degree = a->degree + b->degree;
	poly_fix_degree(prod);
	poly_mod(res, prod, mod, p);
}

/*------------------------------------------------------------------*/
	/* The following routines are highly performance-critical 
	   for factor base generation. 
	
	   These deal with modular multiplication and squaring 
	   of polynomials. The operands and the modulus are 
	   all packed contiguously into an array of uint32's, 
	   to minimize pointer arithmetic */

#define NUM_POLY_COEFFS (MAX_POLY_DEGREE+1)
#define OP1(i) buf[i]
#define MOD(i) buf[NUM_POLY_COEFFS+i]

	/* The following implements word-based modular multiplication
	   and squaring of polynomials, with each coefficient reduced 
	   modulo a prime p. Because there are no carries from one word 
	   to the next, many simplifications are possible. To multiply 
	   A[] and B[] mod N[], with each quantity a polynomial containing
	   d coefficients, the ordinary word-based method computes:
	  
	   for (words B[i] starting at d and working backwards) {
	  	 accum[] = (accum[] << (one word)) mod N[]
	  	 accum[] = (accum[] + A[] * B[i]) mod N[]  (*)
	   }
	  
	   In the finite field case, A[] * B[i] still has d words, and
	   the shifted up quantity has d+1 words. N[] is assumed monic 
	   so the remainder operation amounts to subtracting N[] * (top
	   word of accumulator). This means that going into the multiply-
	   accumulate operation (*) we know the factor by which N[] is 
	   multiplied, so the multiply-accumulate and the subtraction 
	   of N[] can happen at the same time.
	  
	   Another optimization is to allow the coefficents of the 
	   accumulator to reach p^2 in size instead of p; this reduces the
	   number of 32-bit mod operations to only 2*d instead of d*(d+2).
	   There is a mod operation every time the top word of accum[]
	   is converted from mod p^2 to mod p, and there are d mods when
	   the final answer is copied from the accumulator */

/*------------------------------------------------------------------*/
static INLINE uint32 mul_mac(uint32 a, uint32 b, uint32 c, 
				uint32 q, uint32 _mod, 
				uint32 p, uint64 psq) {

	/* For 32-bit a, b, q, c, _mod and 64-bit psq compute

	   (a * b + c - q * _mod) % p

	   psq is the square of p, and the 32-bit quantities 
	   are all assumed less than p
	
	   Note that this is intended for a high-performance
	   processor with conditional move instructions. If the compiler 
	   does not inline, or use these instructions, or realize that
	   many inputs are cached from previous calls then the overall
	   performance will be dismal. */

#if defined(GCC_ASM32A) && \
	!(defined(__GNUC__) && __GNUC__ < 3 ) && \
	defined(NDEBUG)

	uint32 ans;
	ASM_G(
	    "movl %2, %%eax             \n\t"
	    "mull %1                    \n\t"
	    "addl %3, %%eax             \n\t"
	    "adcl $0, %%edx             \n\t"
	    "movl %%eax, %%esi          \n\t"
	    "movl %%edx, %%edi          \n\t"
	    "movl %4, %%eax             \n\t"
	    "mull %5                    \n\t"
	    "subl %%eax, %%esi          \n\t"
	    "sbbl %%edx, %%edi          \n\t"
	    "movl $0, %%eax             \n\t"
	    "movl $0, %%edx             \n\t"
	    "cmovbl %7, %%eax           \n\t"
	    "cmovbl 4+%7, %%edx         \n\t"
	    "addl %%esi, %%eax          \n\t"
	    "adcl %%edi, %%edx          \n\t"
	    "divl %6                    \n\t"
	    "movl %%edx, %0             \n"
	    :"=g"(ans)
	    :"g"(a), "g"(b), "g"(c), "g"(_mod), 
	     "g"(q), "g"(p), "g"(psq)
	    :"%eax", "%esi", "%edi", "%edx", "cc");

	return ans;

#elif defined(MSC_ASM32A)

	uint32 ans;
	ASM_M
	{
		mov	eax,b
		mul	a
		add	eax,c
		adc	edx,0
		mov	esi,eax
		mov	edi,edx
		mov	eax,_mod
		mul	q
		sub	esi,eax
		sbb	edi,edx
		mov	eax,0
		mov	edx,0
		lea	ecx,[psq]
		cmovb	eax,[ecx]
		cmovb	edx,[ecx+4]
		add	eax,esi
		adc	edx,edi
		div	p
		mov	ans,edx
	}

	return ans;

#else
	uint64 ans, tmp;

	/* for a,b,c < p the following will never exceed psq,
	   so no correction is needed */
	ans = (uint64)a * (uint64)b + (uint64)c;

	tmp = (uint64)q * (uint64)_mod;
	ans = ans - tmp + (tmp > ans ? psq : 0);
	return mp_mod64(ans, p);
#endif
}

/*------------------------------------------------------------------*/
static INLINE uint64 sqr_mac(uint32 a, uint32 b, uint64 c, 
				uint32 q, uint32 _mod, uint64 psq) {

	/* For 32-bit a, b, q, _mod and 64-bit c, psq compute

	   (a * b + c - q * _mod) % psq

	   psq is the square of a prime, c is less than psq 
	   and the 32-bit quantities are all assumed less 
	   than the prime */
	
	uint64 ans;

#if defined(GCC_ASM32A) && \
	!(defined(__GNUC__) && __GNUC__ < 3 ) && \
	defined(NDEBUG)

	ASM_G(
	    "movl %1, %%eax               \n\t"
	    "mull %2                      \n\t"
	    "movl %6, %%esi               \n\t"
	    "movl 4+%6, %%edi             \n\t"
	    "subl %3, %%esi               \n\t"
	    "sbbl 4+%3, %%edi             \n\t"
	    "subl %%esi, %%eax            \n\t"
	    "sbbl %%edi, %%edx            \n\t"
	    "movl $0, %%esi               \n\t"
	    "movl $0, %%edi               \n\t"
	    "cmovbl %6, %%esi             \n\t"
	    "cmovbl 4+%6, %%edi           \n\t"
	    "addl %%eax, %%esi            \n\t"
	    "adcl %%edx, %%edi            \n\t"
	    "movl %5, %%eax               \n\t"
	    "mull %4                      \n\t"
	    "subl %%eax, %%esi            \n\t"
	    "sbbl %%edx, %%edi            \n\t"
	    "movl $0, %%eax               \n\t"
	    "movl $0, %%edx               \n\t"
	    "cmovbl %6, %%eax             \n\t"
	    "cmovbl 4+%6, %%edx           \n\t"
	    "addl %%esi, %%eax            \n\t"
	    "adcl %%edi, %%edx            \n\t"
	    "movl %%eax, %0               \n\t"
#if !defined(_ICL_WIN_)
	    "movl %%edx, 4+%0             \n"
#else   /* possible bug in Intel compiler ? */
	    "addl $4, %0                  \n\t"
	    "movl %%edx, %0               \n"
#endif
	    :"=g"(ans)
	    :"g"(a), "g"(b), "g"(c), "g"(q),
	     "g"(_mod), "g"(psq)
	    :"%eax", "%esi", "%edi", "%edx", "cc");

#elif defined(MSC_ASM32A)
	ASM_M
	{
		mov	eax, a
		mul	b
		lea	ecx,[psq]
		mov	esi,[ecx]
		mov	edi,[ecx+4]
		lea	ecx,[c]
		sub	esi,[ecx]
		sbb	edi,[ecx+4]
		sub	eax,esi
		sbb	edx,edi
		mov	esi,0
		mov	edi,0
		lea	ecx,[psq]
		cmovb	esi,[ecx]
		cmovb	edi,[ecx+4]
		add	esi,eax
		adc	edi,edx
		mov	eax,_mod
		mul	q
		sub	esi,eax
		sbb	edi,edx
		mov	eax,0
		mov	edx,0
		cmovb	eax,[ecx]
		cmovb	edx,[ecx+4]
		add	eax,esi
		adc	edx,edi
		lea	ecx,[ans]
		mov	[ecx],eax
		mov	[ecx+4],edx
	}
#else
	uint64 tmp;
	ans = (uint64)a * (uint64)b;
	tmp = psq - c;
	ans = ans - tmp + (tmp > ans ? psq : 0);
	tmp = (uint64)q * (uint64)_mod;
	ans = ans - tmp + (tmp > ans ? psq : 0);
#endif
	return ans;
}

/*------------------------------------------------------------------*/
static INLINE uint64 sqr_mac0(uint32 a, uint32 b,
			uint32 q, uint32 _mod, uint64 psq) {

	/* For 32-bit a, b, q, _mod, compute

	   (a * b - q * _mod) % psq

	   psq is the square of a prime and the 32-bit quantities 
	   are all assumed less than the prime */

	uint64 ans;

#if defined(GCC_ASM32A) && \
	!(defined(__GNUC__) && __GNUC__ < 3 ) && \
	defined(NDEBUG)

	ASM_G(
	    "movl %1, %%eax               \n\t"
	    "mull %2                      \n\t"
	    "movl %%eax, %%esi            \n\t"
	    "movl %%edx, %%edi            \n\t"
	    "movl %4, %%eax               \n\t"
	    "mull %3                      \n\t"
	    "subl %%eax, %%esi            \n\t"
	    "sbbl %%edx, %%edi            \n\t"
	    "movl $0, %%eax               \n\t"
	    "movl $0, %%edx               \n\t"
	    "cmovbl %5, %%eax             \n\t"
	    "cmovbl 4+%5, %%edx           \n\t"
	    "addl %%esi, %%eax            \n\t"
	    "adcl %%edi, %%edx            \n\t"
	    "movl %%eax, %0               \n\t"
	    "movl %%edx, 4+%0             \n"
	    :"=g"(ans)
	    :"g"(a), "g"(b), "g"(q), "g"(_mod), "g"(psq)
	    :"%eax", "%esi", "%edi", "%edx", "cc");

#elif defined(MSC_ASM32A)
	ASM_M
	{
		mov	eax,a
		mul	b
		mov	esi,eax
		mov	edi,edx
		mov	eax,_mod
		mul	q
		sub	esi,eax
		sbb	edi,edx
		mov	eax,0
		mov	edx,0
		lea	ecx,[psq]
		cmovb	eax,[ecx]
		cmovb	edx,[ecx+4]
		add	eax,esi
		adc	edx,edi
		lea	ecx,[ans]
		mov	[ecx],eax
		mov	[ecx+4],edx
	}
#else
	uint64 tmp;
	ans = (uint64)a * (uint64)b;
	tmp = (uint64)q * (uint64)_mod;
	ans = ans - tmp + (tmp > ans ? psq : 0);
#endif
	return ans;
}

/*------------------------------------------------------------------*/
static void poly_expo_modmul(uint32 *buf, uint32 dm, uint32 shift,
			   uint32 p, uint64 psq) { 

	/* OP1 = OP1 * (x - shift) mod MOD
	   OP1 and MOD are of degree dm */

	uint32 q;
	uint32 zero = 0;

	q = OP1(dm-1);
	switch(dm-1) {
	case 7: OP1(7) = mul_mac(OP1(7), shift, OP1(6), q, MOD(7), p, psq);
	case 6: OP1(6) = mul_mac(OP1(6), shift, OP1(5), q, MOD(6), p, psq);
	case 5: OP1(5) = mul_mac(OP1(5), shift, OP1(4), q, MOD(5), p, psq);
	case 4: OP1(4) = mul_mac(OP1(4), shift, OP1(3), q, MOD(4), p, psq);
	case 3: OP1(3) = mul_mac(OP1(3), shift, OP1(2), q, MOD(3), p, psq);
	case 2: OP1(2) = mul_mac(OP1(2), shift, OP1(1), q, MOD(2), p, psq);
	case 1: OP1(1) = mul_mac(OP1(1), shift, OP1(0), q, MOD(1), p, psq);
	case 0: OP1(0) = mul_mac(OP1(0), shift,   zero, q, MOD(0), p, psq);
		break;
	}
}

/*------------------------------------------------------------------*/
static void poly_expo_square(uint32 *buf, uint32 dm, uint32 p, uint64 psq) { 

	/* OP1 = OP1 * OP1 mod MOD
	   OP1 and MOD are both of degree dm */

	uint32 i;
	uint32 q;
	uint64 acc[NUM_POLY_COEFFS];

	for (i = 0; i < dm; i++)
		acc[i] = (uint64)(OP1(i)) * (uint64)(OP1(dm-1));

	for (i = dm - 2; (int32)i >= 0; i--) {
		q = mp_mod64(acc[dm-1], p);
		switch(dm-1) {
  		case 7: acc[7] = sqr_mac(OP1(7), OP1(i), acc[6], 
  							q, MOD(7), psq);
		case 6: acc[6] = sqr_mac(OP1(6), OP1(i), acc[5], 
							q, MOD(6), psq);
		case 5: acc[5] = sqr_mac(OP1(5), OP1(i), acc[4], 
							q, MOD(5), psq);
		case 4: acc[4] = sqr_mac(OP1(4), OP1(i), acc[3], 
							q, MOD(4), psq);
		case 3: acc[3] = sqr_mac(OP1(3), OP1(i), acc[2], 
							q, MOD(3), psq);
		case 2: acc[2] = sqr_mac(OP1(2), OP1(i), acc[1], 
							q, MOD(2), psq);
		case 1: acc[1] = sqr_mac(OP1(1), OP1(i), acc[0], 
							q, MOD(1), psq);
		case 0: acc[0] = sqr_mac0(OP1(0), OP1(i), q, MOD(0), psq);
			break;
		}
	}

	for (i = 0; i < dm; i++)
		OP1(i) = mp_mod64(acc[i], p);
}

/*------------------------------------------------------------------*/
static void poly_xpow(poly_t res, uint32 shift, uint32 n, 
			poly_t mod, uint32 p) { 

	/* Modular exponentiation of polynomials with 
	   finite-field coefficients, i.e. res = (x-shift) ^ n % mod
	   with all polynomial coefficients reduced modulo p.
	   n is assumed nonzero */

	poly_t modnorm;
	uint32 msw;
	uint32 i, d;
	uint32 buf[2*NUM_POLY_COEFFS] = {0};
	uint64 psq;

	poly_make_monic(modnorm, mod, p);
	d = modnorm->degree;

	OP1(0) = shift;
	OP1(1) = 1;
	for (i = 0; i <= d; i++)
		MOD(i) = modnorm->coef[i];

	msw = 0x80000000;
	while (!(n & msw)) {
		msw >>= 1;
	}
	msw >>= 1;
	 
	psq = (uint64)p * (uint64)p;

	/* use left-to-right binary exponentiation, not
	   the right-to-left variety. For factor base generation
	   the base always has degree less than modnorm, and the 
	   left-to-right method preserves that, saving time during
	   modular multiplication */

	while (msw) {
		poly_expo_square(buf, d, p, psq);
		if (n & msw) {
			poly_expo_modmul(buf, d, shift, p, psq);
		}
		msw >>= 1;
	}

	res->degree = d;
	for (i = 0; i <= d; i++)
		res->coef[i] = OP1(i);
	poly_fix_degree(res);
}

/*------------------------------------------------------------------*/
static void poly_gcd(poly_t g_in, poly_t h_in, uint32 p) { 

	poly_t g, h;

	/* make sure the first GCD iteration actually
	   performs useful work */

	if (g_in->degree > h_in->degree) {
		poly_cp(g, g_in);
		poly_cp(h, h_in);
	}
	else {
		poly_cp(h, g_in);
		poly_cp(g, h_in);
	}

	while ((h->degree > 0) || (h->coef[h->degree])) {
		poly_t r;
		poly_mod(r, g, h, p);
		poly_cp(g, h);
		poly_cp(h, r);
	}
	if (g->degree == 0)
		g->coef[0] = 1;
	poly_cp(g_in, g);
}

/*------------------------------------------------------------------*/
static void get_zeros_rec(uint32 *zeros, uint32 shift, 
			uint32 *num_zeros, poly_t f, uint32 p) {

	/* get the zeros of a poly, f, that is known to split
	   completely over Z/pZ. Many thanks to Bob Silverman 
	   for a neat implementation of Cantor-Zassenhaus splitting */

	poly_t g, xpow;
	uint32 degree1, degree2;

	/* base cases of the recursion: we can find the roots
	   of linear and quadratic polynomials immediately */

	if (f->degree == 1) {
		uint32 w = f->coef[1];
		if (w != 1) {
			w = mp_modinv_1(w, p);
			zeros[(*num_zeros)++] = mp_modmul_1(p - f->coef[0],w,p);
		}
		else {
			zeros[(*num_zeros)++] = (f->coef[0] == 0 ? 0 : 
							p - f->coef[0]);
		}
		return;
	}
	else if (f->degree == 2) {

		/* if f is a quadratic polynomial, then it will 
		   always have two distinct nonzero roots or else
		   we wouldn't have gotten to this point. The two 
		   roots are the solution of a general quadratic 
		   equation, mod p */

		uint32 d = mp_modmul_1(f->coef[0], f->coef[2], p);
		uint32 root1 = p - f->coef[1];
		uint32 root2 = root1;
		uint32 ainv = mp_modinv_1(
				mp_modadd_1(f->coef[2], f->coef[2], p),
				p);

		d = mp_modsub_1(mp_modmul_1(f->coef[1], f->coef[1], p),
				mp_modmul_1(4, d, p),
				p);
		d = mp_modsqrt_1(d, p);

		root1 = mp_modadd_1(root1, d, p);
		root2 = mp_modsub_1(root2, d, p);
		zeros[(*num_zeros)++] = mp_modmul_1(root1, ainv, p);
		zeros[(*num_zeros)++] = mp_modmul_1(root2, ainv, p);
		return;
	}

	/* For an increasing sequence of integers 's', compute 
	   the polynomial gcd((x-s)^(p-1)/2 - 1, f). If the result is
	   not g = 1 or g = f, this is a nontrivial splitting 
	   of f. References require choosing s randomly, but however
	   s is chosen there is a 50% chance that it will split f.
	   Since only 0 <= s < p is valid, we choose each s in turn;
	   choosing random s allows the possibility that the same
	   s gets chosen twice (mod p), which would waste time */

	while (shift < p) {
		poly_xpow(xpow, shift, (p-1)/2, f, p);

		poly_cp(g, xpow);
		g->coef[0] = mp_modsub_1(g->coef[0], 1, p);
		poly_fix_degree(g);

		poly_gcd(g, f, p);

		if (g->degree > 0 && g->degree < f->degree)
			break;
		shift++;
	}

	/* f was split; repeat the splitting process on
	   the two halves of f. The linear factors of f are
	   either somewhere in x^((p-1)/2) - 1, in 
	   x^((p-1)/2) + 1, or 'shift' itself is a linear
	   factor. Test each of these possibilities in turn.
	   In the first two cases, begin trying values of s
	   strictly greater than have been tried thus far */

	degree1 = g->degree;
	get_zeros_rec(zeros, shift + 1, num_zeros, g, p);

	poly_cp(g, xpow);
	g->coef[0] = mp_modadd_1(g->coef[0], 1, p);
	poly_fix_degree(g);
	poly_gcd(g, f, p);
	degree2 = g->degree;

	if (degree2 > 0)
		get_zeros_rec(zeros, shift + 1, num_zeros, g, p);

	if (degree1 + degree2 < f->degree)
		zeros[(*num_zeros)++] = (shift == 0 ? 0 : p - shift);
}

/*------------------------------------------------------------------*/
static void poly_reduce_mod_p(poly_t res, mpz_poly_t *_f, uint32 p) {

	uint32 i;

	res->degree = _f->degree;
	for (i = 0; i <= _f->degree; i++)
		res->coef[i] = mpz_fdiv_ui(_f->coeff[i], p);
	poly_fix_degree(res);
}

/*------------------------------------------------------------------*/
uint32 poly_get_zeros(uint32 *zeros, mpz_poly_t *_f, 
			uint32 p, uint32 *high_coeff,
			uint32 count_only) { 

        /* Find all roots of multiplicity 1 for polynomial _f,
	   when the coefficients of _f are reduced mod p. 
	   The leading coefficient of _f mod p is returned
	   
	   Make count_only nonzero if only the number of roots
	   and not their identity matters; this is much faster */

	poly_t g, f;
	uint32 i, j, num_zeros;

	/* reduce the coefficients mod p */

	poly_reduce_mod_p(f, _f, p);
	*high_coeff = f->coef[_f->degree];

	/* bail out if the polynomial is zero */

	if (f->degree == 0)
		return 0;

	/* pull out roots of zero. We do this early to
	   avoid having to handle degree-1 polynomials
	   in later code */

	num_zeros = 0;
	if (f->coef[0] == 0) {
		for (i = 1; i <= f->degree; i++) {
			if (f->coef[i])
				break;
		}
		for (j = i; i <= f->degree; i++) {
			f->coef[i - j] = f->coef[i];
		}
		f->degree = i - j - 1;
		zeros[num_zeros++] = 0;
	}

	/* handle trivial cases */

	if (f->degree == 0) {
		return num_zeros;
	}
	else if (f->degree == 1) {
		uint32 w = f->coef[1];

		if (count_only)
			return num_zeros + 1;

		if (w != 1) {
			w = mp_modinv_1(w, p);
			zeros[num_zeros++] = mp_modmul_1(p - f->coef[0], 
						w, p);
		}
		else {
			zeros[num_zeros++] = (f->coef[0] == 0 ? 
						0 : p - f->coef[0]);
		}
		return num_zeros;
	}

	/* the rest of the algorithm assumes p is odd, which
	   will not work for p=2. Fortunately, in that case
	   there are only two possible roots, 0 and 1. The above
	   already tried 0, so try 1 here */

	if (p == 2) {
		uint32 parity = 0;
		for (i = 0; i <= f->degree; i++)
			parity ^= f->coef[i];
		if (parity == 0)
			zeros[num_zeros++] = 1;
		return num_zeros;
	}
	 
	/* Compute g = gcd(f, x^(p-1) - 1). The result is
	   a polynomial that is the product of all the linear
	   factors of f. A given factor only occurs once in
	   this polynomial */

	poly_xpow(g, 0, p-1, f, p);
	g->coef[0] = mp_modsub_1(g->coef[0], 1, p);
	poly_fix_degree(g);
	poly_gcd(g, f, p);

	/* no linear factors, no service */

	if (g->degree < 1 || count_only)
		return num_zeros + g->degree;

	/* isolate the linear factors */

	get_zeros_rec(zeros, 0, &num_zeros, g, p);
	return num_zeros;
}

/*------------------------------------------------------------------*/
uint32 poly_get_zeros_and_mult(uint32 *zeros, uint32 *mult,
				mpz_poly_t *_f, uint32 p,
				uint32 *high_coeff) {

	uint32 i;
	uint32 num_roots;
	poly_t f;

	num_roots = poly_get_zeros(zeros, _f, p, high_coeff, 0);
	if (num_roots == 0)
		return num_roots;

	poly_reduce_mod_p(f, _f, p);
	for (i = 0; i < num_roots; i++)
		mult[i] = 0;
	if (f->degree == num_roots)
		return num_roots;

	for (i = 0; i < num_roots; i++) {

		poly_t g, r;
		uint32 root = zeros[i];

		g->degree = 2;
		g->coef[0] = mp_modmul_1(root, root, p);
		g->coef[1] = p - mp_modadd_1(root, root, p);
		g->coef[2] = 1;

		poly_mod(r, f, g, p);
		if (r->degree == 0)
			mult[i] = 1;
	}
	return num_roots;
}

/*------------------------------------------------------------------*/
static void poly_xpow_pd(poly_t res, uint32 p, uint32 d, poly_t f) {

	/* compute x^(p^d) mod f */

	uint32 i;
	mpz_t exponent;
	poly_t x;

	mpz_init_set_ui(exponent, p);
	mpz_pow_ui(exponent, exponent, d);

	x->degree = 1;
	x->coef[0] = 0;
	x->coef[1] = 1;

	poly_cp(res, x);
	for (i = mpz_sizeinbase(exponent, 2) - 2; (int32)i >= 0; i--) {
		poly_modmul(res, res, res, f, p);
		if (mpz_tstbit(exponent, i))
			poly_modmul(res, res, x, f, p);
	}

	mpz_clear(exponent);
}

/*------------------------------------------------------------------*/
uint32 is_irreducible(mpz_poly_t *poly, uint32 p) {

	/* this uses Proposition 3.4.4 of H. Cohen, "A Course
	   in Computational Algebraic Number Theory". The tests
	   below are much simpler than trying to factor 'poly' */

	uint32 i;
	poly_t f, tmp;

	poly_reduce_mod_p(f, poly, p);
	poly_make_monic(f, f, p);

	/* in practice, the degree of f will be 8 or less,
	   and we want to compute GCDs for all prime numbers
	   that divide the degree. For this limited range
	   the loop below avoids duplicated code */

	for (i = 2; i < f->degree; i++) {
		if (f->degree % i)
			continue;

		/* for degree d, compute x^(p^(d/i)) - x */

		poly_xpow_pd(tmp, p, f->degree / i, f);
		if (tmp->degree == 0) {
			tmp->degree = 1;
			tmp->coef[1] = p - 1;
		}
		else {
			tmp->coef[1] = mp_modsub_1(tmp->coef[1], 
						(uint32)1, p);
			poly_fix_degree(tmp);
		}

		/* this must be relatively prime to f */

		poly_gcd(tmp, f, p);
		if (tmp->degree > 0 || tmp->coef[0] != 1) {
			return 0;
		}
	}

	/* final test: x^(p^d) mod f must equal x */

	poly_xpow_pd(tmp, p, f->degree, f);
	if (tmp->degree == 1 && tmp->coef[0] == 0 && tmp->coef[1] == 1)
		return 1;
	return 0;
}

/*------------------------------------------------------------------*/
#define NUM_ISQRT_RETRIES 1000

uint32 inv_sqrt_mod_q(mpz_poly_t *res, mpz_poly_t *s_in, mpz_poly_t *f_in,
			uint32 q, uint32 *rand_seed1, uint32 *rand_seed2) {

	/* find a polynomial res(x) such that (res * res * s_in(x)) == 1 
	   mod f_in(x). The algorithm used is from Per Leslie Jensen's
	   thesis 'Integer Factorization', though I haven't been able
	   to find a reference to it anywhere else. Hendrik Lenstra writes
	   that it is a variation on Cantor-Zassenhaus */

	uint32 i, j;
	mpz_t exponent;
	poly_t f, s, y0, y1;

	/* initialize */

	poly_reduce_mod_p(f, f_in, q);
	poly_reduce_mod_p(s, s_in, q);
	poly_make_monic(f, f, q);

	/* none of this will work if s(x) has zero degree */

	if (s->degree == 0)
		return 0;

	/* compute q ^ (degree(f)) */

	mpz_init_set_ui(exponent, q);
	mpz_pow_ui(exponent, exponent, f->degree);

	/* the algorithm is probabilistic; try a few different
	   initial values, then give up if no answer is found */

	for (i = 0; i < NUM_ISQRT_RETRIES; i++) {
		poly_t r0, r1;

		/* form the polynomial r0(x) + r1(x)*Y, where Y is a
		   symbolic variable, r0(x) is randomly chosen and
		   r1(x) = -1 */

		for (j = 0; j < f->degree; j++)
			r0->coef[j] = get_rand(rand_seed1, rand_seed2) % q;
		r0->degree = f->degree - 1;
		poly_fix_degree(r0);

		r1->degree = 0;
		r1->coef[0] = q - 1;

		/* compute 
		   [ (r0(x) - r1(x)*Y) ^ ((exponent-1)/2) ] mod (Y^2 - s(x))

		   The result is a polynomial y0(x) - y1(x) * Y,
		   and y1(x) will be the desired inverse square 
		   root (or its negative) about 50% of the time.
		   Polynomials in x are reduced mod f(x) and their
		   coefficients are reduced mod q */

		poly_cp(y0, r0);
		poly_cp(y1, r1);
		for (j = mpz_sizeinbase(exponent, 2) - 2; j; j--) {

			poly_t p0, p1, p2;

			/* square (y0, y1) */

			poly_modmul(p2, y1, y1, f, q);
			poly_modmul(p2, p2, s, f, q);
			poly_modmul(p1, y0, y1, f, q);
			poly_add(p1, p1, p1, q);
			poly_modmul(p0, y0, y0, f, q);
			poly_add(p0, p0, p2, q);
			poly_cp(y0, p0);
			poly_cp(y1, p1);

			if (!mpz_tstbit(exponent, j))
				continue;

			/* multiply (y0, y1) by (r0, r1) */

			poly_modmul(p2, y1, r0, f, q);
			poly_modmul(p1, y0, r1, f, q);
			poly_add(p1, p1, p2, q);
			poly_modmul(p2, y1, r1, f, q);
			poly_modmul(p2, p2, s, f, q);
			poly_modmul(p0, y0, r0, f, q);
			poly_add(p0, p0, p2, q);
			poly_cp(y0, p0);
			poly_cp(y1, p1);
		}
		
		/* check if y1(x) is a valid inverse square root */

		poly_modmul(y0, y1, y1, f, q);
		poly_modmul(y0, y0, s, f, q);
		if (y0->degree == 0 && y0->coef[0] == 1)
			break;
	}

	mpz_clear(exponent);

	/* if no inverse square root was found and q is small enough,
	   attempt to find an inverse square root by brute force,
	   trying all q^d elements of the finite field.

	   We can save half the time by avoiding polynomials that 
	   are the negative of polynomials already searched */
	  
	if (i == NUM_ISQRT_RETRIES && q < 150) {
		uint32 c0, c1, c2, c3, c4, c5, c6, c7;
		uint32 start[MAX_POLY_DEGREE];
 
		for (i = 0; i < f->degree; i++)
			start[i] = q - 1;
		for (; i < MAX_POLY_DEGREE; i++)
			start[i] = 0;
		y1->degree = 7;
		start[y1->degree] /= 2;
	
		for (c7 = start[7]; (int32)c7 >= 0; c7--) {
			y1->coef[7] = c7;
			if (c7 == 0) 
				start[--y1->degree] /= 2;

		for (c6 = start[6]; (int32)c6 >= 0; c6--) {
			y1->coef[6] = c6;
			if (c6 == 0 && y1->degree == 6) 
				start[--y1->degree] /= 2; 

		for (c5 = start[5]; (int32)c5 >= 0; c5--) {
			y1->coef[5] = c5;
			if (c5 == 0 && y1->degree == 5) 
				start[--y1->degree] /= 2; 

		for (c4 = start[4]; (int32)c4 >= 0; c4--) {
			y1->coef[4] = c4;
			if (c4 == 0 && y1->degree == 4) 
				start[--y1->degree] /= 2; 

		for (c3 = start[3]; (int32)c3 >= 0; c3--) {
			y1->coef[3] = c3;
			if (c3 == 0 && y1->degree == 3) 
				start[--y1->degree] /= 2; 

		for (c2 = start[2]; (int32)c2 >= 0; c2--) {
			y1->coef[2] = c2;
			if (c2 == 0 && y1->degree == 2) 
				start[--y1->degree] /= 2; 

		for (c1 = start[1]; (int32)c1 >= 0; c1--) {
			y1->coef[1] = c1;
			if (c1 == 0 && y1->degree == 1) 
				start[--y1->degree] /= 2;

		for (c0 = start[0]; (int32)c0 >= 0; c0--) {
			y1->coef[0] = c0;
			poly_modmul(y0, y1, y1, f, q);
			poly_modmul(y0, y0, s, f, q);
			if (y0->degree == 0 && y0->coef[0] == 1)
				goto finished;
		}}}}}}}}
	}

finished:
	if (y0->degree == 0 && y0->coef[0] == 1) {
		res->degree = y1->degree;
		for (i = 0; i <= y1->degree; i++)
			mpz_set_ui(res->coeff[i], y1->coef[i]);
		return 1;
	}

	/* no luck; give up */
	return 0;
}
