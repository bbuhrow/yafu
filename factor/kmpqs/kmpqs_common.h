/*
kmpqs: Kleinjung's small-number MPQS from lasieve5 (gnfs4linux), ported
into yafu.

The algorithms are Copyright (C) 2001 Jens Franke, T. Kleinjung, and this
module is distributed under the terms of the GNU General Public Licence and
WITHOUT ANY WARRANTY.  See factor/lasieve5_64/COPYING.

Internal glue: replaces lasieve's if.h / gmp-aux.h / psp.c / siever-config.h.
*/

#ifndef KMPQS_COMMON_H
#define KMPQS_COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <setjmp.h>
#include <math.h>
#include "gmp.h"
#include "kmpqs.h"
#include "../../top/aprcl/tinyprp.h"

/*
 * lasieve's Schlendrian()/complain() call exit().  Here an internal error
 * longjmps back to the entry point of the current factorization, which
 * then returns KMPQS_ERR_INTERNAL.  Each source file defines KMPQS_JMPBUF.
 */
#define KMPQS_ERR_INTERNAL (-100)

static void kmpqs_report(const char *fmt, ...)
{
	va_list ap;

	if (!kmpqs_verbose)
		return;
	va_start(ap, fmt);
	fprintf(stderr, "kmpqs: ");
	vfprintf(stderr, fmt, ap);
	va_end(ap);
}

#define KMPQS_FAIL(...) \
	do { kmpqs_report(__VA_ARGS__); longjmp(KMPQS_JMPBUF, 1); } while (0)

static void *kmpqs_xmalloc(size_t size)
{
	void *p = malloc(size);

	if (p == NULL) {
		fprintf(stderr, "kmpqs: out of memory (%zu bytes)\n", size);
		exit(1);
	}
	return p;
}

/* low 64 bits of |z| */
static inline uint64_t kmpqs_mpz_get_u64(mpz_t z)
{
#if GMP_LIMB_BITS == 64
	return (uint64_t)mpz_getlimbn(z, 0);
#else
	return (uint64_t)mpz_getlimbn(z, 0) | ((uint64_t)mpz_getlimbn(z, 1) << 32);
#endif
}

static inline void kmpqs_mpz_set_u64(mpz_t z, uint64_t v)
{
	mpz_set_ui(z, (unsigned long)(v >> 32));
	mpz_mul_2exp(z, z, 32);
	mpz_add_ui(z, z, (unsigned long)(v & 0xffffffffULL));
}

/* Base-2 strong probable prime test (lasieve's psp()); returns 0 or 1. */
static inline int kmpqs_psp(mpz_t n)
{
	size_t nbits;

	if (mpz_sgn(n) <= 0)
		return 0;
	nbits = mpz_sizeinbase(n, 2);
	if (nbits <= 2)
		return mpz_cmp_ui(n, 1) > 0;
	if (mpz_even_p(n))
		return 0;
	if (nbits <= 64)
		return MR_2sprp_64x1(kmpqs_mpz_get_u64(n)) ? 1 : 0;
	if (nbits <= 128) {
		uint64_t v[2];
		mpz_t t;

		v[0] = kmpqs_mpz_get_u64(n);
		mpz_init(t);
		mpz_tdiv_q_2exp(t, n, 64);
		v[1] = kmpqs_mpz_get_u64(t);
		mpz_clear(t);
		return MR_2sprp_128x1(v) ? 1 : 0;
	}
	return mpz_probab_prime_p(n, 1) ? 1 : 0;
}

/* Counters lasieve's mpqstest/mpqs3test print; defined by the test driver. */
#if defined(MPQS_STAT) || defined(MPQS3_STAT)
extern uint64_t stat_asm_eval, stat_asm_td;
extern uint64_t stat_td_cand, stat_td_surv;
extern uint64_t stat_ff, stat_pf, stat_comb, stat_pp;
extern uint32_t stat_asm_div, stat_final_mulmod;
extern uint32_t stat_counter0, stat_counter1, stat_retry;
extern uint32_t stat_size[15];
#endif

#endif
