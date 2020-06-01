/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.

$Id: batch_factor.c 638 2011-09-11 15:31:19Z jasonp_sf $
--------------------------------------------------------------------*/

#include <batch_factor.h>
#include <stdint.h>
#include "monty.h"
#include "prime_sieve.h"

/*------------------------------------------------------------------*/

#define BREAKOVER_WORDS 50
void multiply_primes(uint32 first, uint32 last,
    prime_sieve_t *sieve, mpz_t prod) {

    /* recursive routine to multiply all the elements of a
       list of consecutive primes. The current invocation
       multiplies elements first to last of the list, inclusive.
       The primes are created on demand */

    mpz_t half_prod;
    uint32 mid = (last + first) / 2;

    /* base case; accumulate a few primes */

    if (last < first + BREAKOVER_WORDS) {
        mpz_set_ui(prod, (unsigned long)get_next_prime(sieve));
        while (++first <= last) {
            mpz_mul_ui(prod, prod,
                (unsigned long)get_next_prime(sieve));
        }
        return;
    }

    /* recursively handle the left and right halves of the list */

    mpz_init(half_prod);
    multiply_primes(first, mid, sieve, prod);
    multiply_primes(mid + 1, last, sieve, half_prod);

    /* multiply them together. We can take advantage of
       fast multiplication since in general the two half-
       products are about the same size*/

    mpz_mul(prod, prod, half_prod);
    mpz_clear(half_prod);
}



/*------------------------------------------------------------------*/
void relation_to_gmp(relation_batch_t *rb,
				uint32 index, mpz_t out) {

	/* multiply together the unfactored parts of an
	   NFS relation, then convert to an mpz_t */

    int j;
	uint32 i, nwords;
	cofactor_t *c = rb->relations + index;
	uint32 *f = rb->factors + c->factor_list_word + 
			c->num_factors_r + c->num_factors_a;
	
    mpz_ptr num1 = rb->f1a, num2 = rb->f2a, num3 = rb->n;

	if (c->lp_r_num_words > 1 && c->lp_a_num_words > 1) {

		/* rational and algebraic parts need to be
		   multiplied together first */

		//mp_clear(&num1);
		//mp_clear(&num2);
		//num1.nwords = c->lp_r_num_words;
		//for (i = 0; i < c->lp_r_num_words; i++)
		//	num1.val[i] = f[i];

        mpz_set_ui(num1, f[0]);
        for (i = 1; i < c->lp_r_num_words; i++)
        {
            mpz_mul_2exp(num1, num1, 32);
            mpz_add_ui(num1, num1, f[i]);
        }

		//num2.nwords = c->lp_a_num_words;
		f += c->lp_r_num_words;
		//for (i = 0; i < c->lp_a_num_words; i++)
		//	num2.val[i] = f[i];
        mpz_set_ui(num2, f[0]);
        for (i = 1; i < c->lp_a_num_words; i++)
        {
            mpz_mul_2exp(num2, num2, 32);
            mpz_add_ui(num2, num2, f[i]);
        }

		mpz_mul(num3, num2, num1);
	}
	else {
		/* no multiply necesary */

		//mp_clear(&num3);
		nwords = c->lp_r_num_words;

		if (c->lp_a_num_words > 1) {
			nwords = c->lp_a_num_words;
			f += c->lp_r_num_words;
		}

		//num3.nwords = nwords;
		//for (i = 0; i < nwords; i++)
		//	num3.val[i] = f[i];

        mpz_set_ui(num3, f[nwords - 1]);
        for (j = nwords - 2; j >= 0; j--)
        {
            mpz_mul_2exp(num3, num3, 32);
            mpz_add_ui(num3, num3, f[j]);
        }
	}

    mpz_set(out, num3);

    return;
}

/*------------------------------------------------------------------*/
void multiply_relations(uint32 first, uint32 last, 
				relation_batch_t *rb,
				mpz_t prod) {
	mpz_t half_prod;

	/* recursive routine to multiply together a list of
	   relations. The current invocation multiplies relations
	   first to last of the list, inclusive */

	/* base case of recursion */

	if (first == last) {
		relation_to_gmp(rb, first, prod);
		return;
	}

	/* recurse on the left and right halves */

	mpz_init(half_prod);
	if (last == first + 1) {
		relation_to_gmp(rb, first, half_prod);
		relation_to_gmp(rb, last, prod);
	}
	else {
		uint32 mid = (last + first) / 2;
		multiply_relations(first, mid, rb, prod);
		multiply_relations(mid + 1, last, rb, half_prod);
	}

	/* multiply the halves */

	mpz_mul(prod, prod, half_prod);
	mpz_clear(half_prod);
}

/*------------------------------------------------------------------*/

uint64_t pow2m(uint64_t b, uint64_t n)
{
    // compute 2^b mod n
    uint64_t acc, x, rho, am, g[8], mask;
    int i;
    int j;
    int bstr;
    int bits = 64 - __builtin_clzll(b);

    x = (((n + 2) & 4) << 1) + n;   // here x*a==1 mod 2**4
    x *= 2 - n * x;                 // here x*a==1 mod 2**8
    x *= 2 - n * x;                 // here x*a==1 mod 2**16
    x *= 2 - n * x;                 // here x*a==1 mod 2**32         
    x *= 2 - n * x;                 // here x*a==1 mod 2**64

    rho = (uint64_t)0 - x;
    am = u64div(2, n);              // put 2 into Monty rep.

    // precomputations, b^i for 0 <= i < 2^k
    g[1] = am;
    g[2] = mulredc(g[1], am, n, rho);
    g[3] = mulredc(g[2], am, n, rho);
    g[4] = sqrredc(g[2], n, rho);
    g[5] = mulredc(g[4], am, n, rho);
    g[6] = sqrredc(g[3], n, rho);
    g[7] = mulredc(g[6], am, n, rho);

    // L-R windowed exponentiation.
    mask = 0x7ULL << (bits - 3);

    // the first 5 iterations can be done with no mod if BITS==64,
    // before we put the accumulator into monty rep.
    // iter1: sqr 1 and mul
    acc = 1;
    bstr = (b & mask) >> (bits - 3);
    if (bstr > 0)
        acc = acc * (1ULL << bstr);
    mask >>= 3;
    bits -= 3;
    acc *= acc;
    acc *= acc;
    bstr = (b & mask) >> (bits - 2);
    if (bstr > 0)
        acc = acc * (1ULL << bstr);
    mask >>= 2;
    bits -= 2;

    acc = u64div(acc, n);             // put acc into Monty rep.
    //acc = u64div(1, n);

    while (bits > 3)
    {
        bstr = (b & mask) >> (bits - 3);
        acc = sqrredc(acc, n, rho);
        acc = sqrredc(acc, n, rho);
        acc = sqrredc(acc, n, rho);
        mask >>= 3;
        if (bstr > 0)
            acc = mulredc(acc, g[bstr], n, rho);
        bits -= 3;
    }

    bstr = (b & mask);
    for (j = 0; mask > 0; j++)
    {
        acc = sqrredc(acc, n, rho);
        mask >>= 1;
    }

    if (bstr > 0)
        acc = mulredc(acc, g[bstr], n, rho);

    // take result out of monty rep.
    acc = mulredc(acc, 1, n, rho);

    // final check to ensure c < N
    if (acc >= n)
    {
        acc -= n;
    }

    return acc;
}

/* the following recursion base-case is specialized for relations
   containing <= 3 rational and/or algebraic large primes. The rest
   of the batch factoring handles an arbitrary number of large primes,
   so only this routine needs to change when more advanced factoring
   code becomes available. */

void check_batch_relation(relation_batch_t *rb,
			uint32 index,
			mpz_t prime_product) {

	uint32 i;
    int j;
	cofactor_t *c = rb->relations + index;
	uint32 *f = rb->factors + c->factor_list_word;
	uint32 *lp1 = f + c->num_factors_r + c->num_factors_a;
	uint32 *lp2 = lp1 + c->lp_r_num_words;
    mpz_ptr f1r = rb->f1r;
    mpz_ptr f1a = rb->f1a;
    mpz_ptr f2r = rb->f2r;
    mpz_ptr f2a = rb->f2a;
    mpz_ptr small = rb->small;
    mpz_ptr large = rb->large;
    mpz_ptr n = rb->n;
	uint32 lp_r[MAX_LARGE_PRIMES];
	uint32 lp_a[MAX_LARGE_PRIMES];
	uint32 num_r, num_a;

	/* first compute gcd(prime_product, rational large cofactor).
	   The rational part will split into a cofactor with all 
	   factors <= the largest prime in rb->prime_product (stored
	   in f1r) and a cofactor with all factors larger (stored in f2r) */

    // the larger we make lp_cutoff_r, the more TLP's we will find, at
    // the cost of having to split more TLP candidates in f1r.

    c->success = 0;
    c->lp_r[0] = c->lp_r[1] = c->lp_r[2] = 1;

	if (c->lp_r_num_words) {

        mpz_set_ui(n, lp1[c->lp_r_num_words - 1]);
        for (j = c->lp_r_num_words - 2; j >= 0; j--)
        {
            mpz_mul_2exp(n, n, 32);
            mpz_add_ui(n, n, lp1[j]);
        }

        //gmp_printf("processing n = %Zd: ", n); // , prime_product has %u bits\n", n, mpz_sizeinbase(prime_product, 2));

		if (mpz_sizeinbase(n, 2) < 32) {
            mpz_set(f1r, n);
            mpz_set_ui(f2r, 1);
		}
		else {
			//mp_gcd(prime_product, &n, &f1r);
            mpz_gcd(f1r, prime_product, n);

            if (mpz_cmp_ui(f1r, 1) == 0)
            {
                mpz_set(f2r, n);
            }
            else
            {
                mpz_tdiv_q(f2r, n, f1r);
            }

            //if (mpz_cmp(f1r, n) == 0)
            //{
            //    gmp_printf("all factors found in gcd, n = %Zd, f1r = %Zd, (%u bits)\n",
            //        n, f1r, mpz_sizeinbase(f1r, 2));
            //}

			//if (mp_is_one(&f1r))
			//	mp_copy(&n, &f2r);
			//else
			//	mp_div(&n, &f1r, &f2r);
		}

        //gmp_printf("f1r = %Zx\n", f1r);
        //gmp_printf("f2r = %Zx\n", f2r);

		/* give up on this relation if
		     - f1r has a single factor, and that factor
		       exceeds the rational large prime cutoff
		     - f1r is more than 3 words long */

		//if (f1r.nwords > 3 ||
		//    (f1r.nwords == 1 && f1r.val[0] > rb->lp_cutoff_r))
		//	return;
        if ((mpz_sizeinbase(f1r, 2) > 96) || 
            ((mpz_sizeinbase(f1r, 2) < 32) && (mpz_get_ui(f1r) > rb->lp_cutoff_r)))
        {
            //printf("abort 1\n");
            return;
        }

		/* give up on this relation if
		     - f2r has a single factor, and that factor
		       exceeds the rational large prime cutoff
		     - f2r has two words and exceeds the square of the
		       cutoff (meaning at least one of the two factors in
		       f2r would exceed the large prime cutoff)
		     - f2r has two words and is smaller than the square
		       of the largest prime in rb->prime_product, meaning
		       that f2r is prime and way too large
		     - f2r has 3+ words. We don't know anything about
		       f2r in this case, except that all factors exceed the
		       largest prime in rb->prime_product, and it would be
		       much too expensive to find out anything more */

		//if (f2r.nwords >= 3 ||
		//    (f2r.nwords == 2 && 
		//     		(mp_cmp(&f2r, &rb->max_prime2) <= 0 ||
		//     		 mp_cmp(&f2r, &rb->lp_cutoff_r2) > 0)) ||
		//    (f2r.nwords == 1 && f2r.val[0] > rb->lp_cutoff_r))
		//	return;

        //if (mpz_cmp(f1r, f2r) > 0)
        //{
        //    mpz_set(n, f1r);
        //    mpz_set(f1r, f2r);
        //    mpz_set(f2r, n);
        //}

        uint32 r_cutoff_bits = spBits(rb->lp_cutoff_r) * 2;

        if ((mpz_sizeinbase(f2r, 2) > r_cutoff_bits) || ((mpz_get_ui(f2r) > 1) &&
            (((mpz_sizeinbase(f2r, 2) < 32) && (mpz_get_ui(f2r) > rb->lp_cutoff_r)) ||
            ((mpz_sizeinbase(f2r, 2) < 64) && 
                ((mpz_cmp(f2r, rb->max_prime2) <= 0) ||
                 (mpz_cmp(f2r, rb->lp_cutoff_r2) > 0))))))
        {
            //gmp_printf("abort 2: f1r = %Zd, f2r = %Zd\n", f1r, f2r);
            return;
        }

        //if ((mpz_sizeinbase(f2r, 2) >= 64) && ((mpz_sizeinbase(f2r, 2) < 96)))
        if (0)
        {
            // experiment with the last case where f2r possibly splits
            // into a TLP but we know that all of the factors are
            // larger than primes in our GCD.  It is expensive to learn more
            // but for huge jobs the expense may be worth it.
            // We expend a small amount of ECM and hope to get lucky.
            if (mpz_probab_prime_p(f2r, 1))
            {
                //printf("abort 3\n");
                return;
            }

            for (i = num_r = num_a = 0; i < MAX_LARGE_PRIMES; i++)
                lp_r[i] = lp_a[i] = 1;

            tinyecm(f2r, small, 70, 70 * 25, 16, 0);

            if ((mpz_sizeinbase(small, 2) > 32) || (mpz_get_ui(small) > rb->lp_cutoff_r))
            {
                //printf("abort 8\n");
                return;
            }

            // small is acceptable, get large part.
            mpz_tdiv_q(large, f2r, small);
            lp_r[num_r++] = mpz_get_ui(small);

            if (mpz_sizeinbase(large, 2) > 64)
            {
                //printf("abort 8\n");
                return;
            }

            uint64 e = mpz_get_ui(large);
            if (pow2m(e - 1, e) == 1)
            {
                //printf("abort 3\n");
                return;
            }

            // large still on track, try to split it
            uint64 f64 = do_uecm(mpz_get_ui(large));

            // see if any factors found are acceptable
            if (f64 <= 1 || f64 > rb->lp_cutoff_r)
            {
                //printf("abort 9\n");
                return;
            }

            lp_r[num_r++] = f64;
            mpz_tdiv_q_ui(large, large, f64);

            if ((mpz_sizeinbase(large, 2) > 32) || (mpz_get_ui(large) > rb->lp_cutoff_r))
            {
                //printf("abort 10\n");
                return;
            }

            lp_r[num_r++] = mpz_get_ui(large);

            /* yay! Another relation found */

            rb->num_success++;

            for (i = 0; i < MIN(3, num_r); i++)
            {
                c->lp_r[i] = lp_r[i];
            }
            c->success = num_r;
            return;
        }


        //if (mpz_sizeinbase(f2r, 2) > 64)
        //{
        //    //gmp_printf("abort 2: f2r > 64 bits, f2r = %Zd\n", f2r);
        //    return;
        //}
        //
        //if (mpz_get_ui(f2r) > 1)
        //{
        //    if ((mpz_sizeinbase(f2r, 2) < 32) && (mpz_get_ui(f2r) > rb->lp_cutoff_r))
        //    {
        //        //gmp_printf("abort 2: f2r > lp_cutoff_r (%u), f2r = %Zd\n", rb->lp_cutoff_r, f2r);
        //        return;
        //    }
        //
        //    if (mpz_cmp(f2r, rb->max_prime2) <= 0)
        //    {
        //        //gmp_printf("abort 2: f2r < max_prime2 (%Zd), f2r = %Zd\n", rb->max_prime2, f2r);
        //        return;
        //    }
        //
        //    if (mpz_cmp(f2r, rb->lp_cutoff_r2) > 0)
        //    {
        //        //gmp_printf("abort 2: f2r > lp_cutoff_r2 (%Zd), f2r = %Zd\n", rb->lp_cutoff_r2, f2r);
        //        return;
        //    }
        //}

	}
	else {
        mpz_set_ui(f1r, 0);
        mpz_set_ui(f2r, 0);
	}

	/* repeat with the algebraic unfactored part, if any */

	if (c->lp_a_num_words) {
        mpz_set_ui(n, lp2[c->lp_a_num_words - 1]);
        for (j = c->lp_a_num_words - 2; j >= 0; j--)
        {
            mpz_mul_2exp(n, n, 32);
            mpz_add_ui(n, n, lp2[j]);
        }

        if (mpz_sizeinbase(n, 2) < 32) {
            mpz_set(f1a, n);
            mpz_set_ui(f2a, 1);
        }
        else {
            mpz_gcd(f1a, prime_product, n);

            if (mpz_cmp_ui(f1a, 1) == 1)
            {
                mpz_set(f2a, n);
            }
            else
            {
                mpz_tdiv_q(f2a, n, f1a);
            }

        }

        if ((mpz_sizeinbase(f1a, 2) >= 96) ||
            ((mpz_sizeinbase(f1a, 2) < 32) && (mpz_get_ui(f1a) > rb->lp_cutoff_a)))
        {
            return;
        }

        if ((mpz_sizeinbase(f2a, 2) >= 64) ||
            ((mpz_sizeinbase(f2a, 2) < 32) && (mpz_get_ui(f2a) > rb->lp_cutoff_a)) ||
            ((mpz_sizeinbase(f2a, 2) < 64) &&
            ((mpz_cmp(f2a, rb->max_prime2) <= 0) ||
                (mpz_cmp(f2a, rb->lp_cutoff_a2) > 0))))
        {
            return;
        }
	}
	else {
        mpz_set_ui(f1a, 0);
        mpz_set_ui(f2a, 0);
	}

	/* the relation isn't obviously bad; do more work
	   trying to factor everything. Note that when relations 
	   are expected to have three large primes then ~98% of 
	   relations do not make it to this point
	
	   Begin by performing compositeness tests on f2r and f2a,
	   which are necessary if they are two words in size and
	   f1r or f1a is not one (the latter being true means
	   that the sieving already performed the compositeness test) */

	for (i = num_r = num_a = 0; i < MAX_LARGE_PRIMES; i++)
		lp_r[i] = lp_a[i] = 1;

    if ((mpz_sizeinbase(f2r, 2) > 32) && (mpz_sizeinbase(f2r, 2) <= 64))
    {
        uint64 e = mpz_get_ui(f2r);
        if (pow2m(e - 1, e) == 1)
        {
            //printf("abort 3\n");
            return;
        }
    }
    if ((mpz_sizeinbase(f2a, 2) > 32) && (mpz_sizeinbase(f2a, 2) <= 64))
    {
        uint64 e = mpz_get_ui(f2a);
        if (pow2m(e - 1, e) == 1)
        {
            return;
        }
    }

	/* now perform all the factorizations that
	   require SQUFOF, since it is much faster than the
	   QS code. We have to check all of f[1|2][r|a] but
	   for relations with three large primes then at most
	   two of the four choices need factoring */

    if (mpz_sizeinbase(f1r, 2) <= 32)
    {
        if (mpz_get_ui(f1r) > 1)
            lp_r[num_r++] = mpz_get_ui(f1r);
    }
    else if (mpz_sizeinbase(f1r, 2) <= 64)
    {
        uint64 f64 = do_uecm(mpz_get_ui(f1r));

        if (f64 <= 1 || f64 > rb->lp_cutoff_r)
        {
            //printf("abort 4\n");
            return;
        }
        lp_r[num_r++] = f64;
        mpz_tdiv_q_ui(f1r, f1r, f64);

        if ((mpz_sizeinbase(f1r, 2) > 32) || (mpz_get_ui(f1r) > rb->lp_cutoff_r))
        {
            //printf("abort 5\n");
            return;
        }

        lp_r[num_r++] = mpz_get_ui(f1r);
    }

    if (mpz_sizeinbase(f2r, 2) <= 32)
    {
        if (mpz_get_ui(f2r) > 1)
            lp_r[num_r++] = mpz_get_ui(f2r);
    }
    else if (mpz_sizeinbase(f2r, 2) <= 64)
    {
        uint64 f64 = do_uecm(mpz_get_ui(f2r));

        if (f64 <= 1 || f64 > rb->lp_cutoff_r)
        {
            //printf("abort 6\n");
            return;
        }

        lp_r[num_r++] = f64;
        mpz_tdiv_q_ui(f2r, f2r, f64);

        if ((mpz_sizeinbase(f2r, 2) > 32) || (mpz_get_ui(f2r) > rb->lp_cutoff_r))
        {
            //printf("abort 7\n");
            return;
        }

        lp_r[num_r++] = mpz_get_ui(f2r);
    }

    if (mpz_sizeinbase(f1a, 2) <= 32)
    {
        if (mpz_get_ui(f1a) > 1)
            lp_a[num_a++] = mpz_get_ui(f1a);
    }
    else if (mpz_sizeinbase(f1a, 2) <= 64)
    {
        uint64 f64 = do_uecm(mpz_get_ui(f1a));

        if (f64 <= 1 || f64 > rb->lp_cutoff_a)
            return;
        lp_a[num_a++] = f64;
        //mp_divrem_1(&f1r, i, &f1r);
        mpz_tdiv_q_ui(f1a, f1a, f64);

        if ((mpz_sizeinbase(f1a, 2) > 32) || (mpz_get_ui(f1a) > rb->lp_cutoff_a))
            return;

        //if (f1r.nwords > 1 || f1r.val[0] > rb->lp_cutoff_r)
        //    return;

        lp_a[num_a++] = mpz_get_ui(f1a);
    }


    if (mpz_sizeinbase(f2a, 2) <= 32)
    {
        if (mpz_get_ui(f2a) > 1)
            lp_a[num_a++] = mpz_get_ui(f2a);
    }
    else if (mpz_sizeinbase(f2a, 2) <= 64)
    {
        uint64 f64 = do_uecm(mpz_get_ui(f2a));

        if (f64 <= 1 || f64 > rb->lp_cutoff_a)
            return;
        lp_a[num_a++] = f64;
        //mp_divrem_1(&f1r, i, &f1r);
        mpz_tdiv_q_ui(f2a, f2a, f64);

        if ((mpz_sizeinbase(f2a, 2) > 32) || (mpz_get_ui(f2a) > rb->lp_cutoff_a))
            return;

        //if (f1r.nwords > 1 || f1r.val[0] > rb->lp_cutoff_r)
        //    return;

        lp_a[num_a++] = mpz_get_ui(f2a);
    }

	/* only use expensive factoring methods when we know 
	   f1r and/or f1a splits into three large primes and 
	   we know all three primes are smaller than the 
	   largest prime in rb->prime_product. When the latter 
	   is a good deal smaller than the large prime cutoff
	   this happens extremely rarely */

	//if (f1r.nwords == 3) {
	//	if (tinyqs(&f1r, &t0, &t1) == 0)
	//		return;
    //
	//	small = &t0;
	//	large = &t1;
	//	if (mp_cmp(small, large) > 0) {
	//		small = &t1;
	//		large = &t0;
	//	}
    //
	//	if (small->nwords > 1 || small->val[0] > rb->lp_cutoff_r)
	//		return;
	//	lp_r[num_r++] = small->val[0];
	//	i = squfof(large);
	//	if (i <= 1 || i > rb->lp_cutoff_r)
	//		return;
	//	lp_r[num_r++] = i;
	//	mp_divrem_1(large, i, large);
	//	if (large->nwords > 1 || large->val[0] > rb->lp_cutoff_r)
	//		return;
	//	lp_r[num_r++] = large->val[0];
	//}

    if (mpz_sizeinbase(f1r, 2) > 64) {
        //gmp_printf("attempting to factor %u-bit n = %Zx\n", mpz_sizeinbase(f1r, 2), f1r);

        //if (tinyqs(qs_params, f1r, small, large) == 0)
        //    goto done;
        //
        //if (mpz_cmp(small, large) > 0) {
        //    mpz_set(n, small);
        //    mpz_set(small, large);
        //    mpz_set(large, n);
        //}

        int B1, B2, curves, bits = mpz_sizeinbase(f1r, 2);
        int targetBits = bits / 3 + 1;
        if (targetBits <= 25)
        {
            B1 = 70;
            curves = 16;
        }
        else if (targetBits <= 26)
        {
            B1 = 85;
            curves = 16;
        }
        else if (targetBits <= 29)
        {
            B1 = 125;
            curves = 16;
        }
        else if (targetBits <= 31)
        {
            B1 = 165;
            curves = 24;
        }
        else if (targetBits <= 32)
        {
            B1 = 205;
            curves = 24;
        }
        else
        {
            printf("something's wrong, bits = %u, targetBits = %u\n", bits, targetBits);
        }

        tinyecm(f1r, small, B1, B1 * 25, curves, 0);

        if ((mpz_sizeinbase(small, 2) > 32) || (mpz_get_ui(small) > rb->lp_cutoff_r))
        {
            //printf("abort 8\n");
            return;
        }

        mpz_tdiv_q(large, f1r, small);
        lp_r[num_r++] = mpz_get_ui(small);

        uint64 f64 = do_uecm(mpz_get_ui(large));

        if (f64 <= 1 || f64 > rb->lp_cutoff_r)
        {
            //printf("abort 9\n");
            return;
        }

        lp_r[num_r++] = f64;
        mpz_tdiv_q_ui(large, large, f64);

        if ((mpz_sizeinbase(large, 2) > 32) || (mpz_get_ui(large) > rb->lp_cutoff_r))
        {
            //printf("abort 10\n");
            return;
        }

        lp_r[num_r++] = mpz_get_ui(large);
    }

	//if (f1a.nwords == 3) {
	//	if (tinyqs(&f1a, &t0, &t1) == 0)
	//		return;
    //
	//	small = &t0;
	//	large = &t1;
	//	if (mp_cmp(small, large) > 0) {
	//		small = &t1;
	//		large = &t0;
	//	}
    //
	//	if (small->nwords > 1 || small->val[0] > rb->lp_cutoff_a)
	//		return;
	//	lp_a[num_a++] = small->val[0];
	//	i = squfof(large);
	//	if (i <= 1 || i > rb->lp_cutoff_a)
	//		return;
	//	lp_a[num_a++] = i;
	//	mp_divrem_1(large, i, large);
	//	if (large->nwords > 1 || large->val[0] > rb->lp_cutoff_a)
	//		return;
	//	lp_a[num_a++] = large->val[0];
	//}

    if (mpz_sizeinbase(f1a, 2) > 64) {
        //if (tinyqs(qs_params, f1a, small, large) == 0)
        //    return;
        //
        //if (mpz_cmp(small, large) > 0) {
        //    mpz_set(n, small);
        //    mpz_set(small, large);
        //    mpz_set(large, n);
        //}

        int B1, B2, curves, bits = mpz_sizeinbase(f1a, 2);
        int targetBits = bits / 3 + 1;
        if (targetBits <= 25)
        {
            B1 = 70;
            curves = 16;
        }
        else if (targetBits <= 26)
        {
            B1 = 85;
            curves = 16;
        }
        else if (targetBits <= 29)
        {
            B1 = 125;
            curves = 16;
        }
        else if (targetBits <= 31)
        {
            B1 = 165;
            curves = 24;
        }
        else if (targetBits <= 32)
        {
            B1 = 205;
            curves = 24;
        }
        else
        {
            printf("something's wrong, bits = %u, targetBits = %u\n", bits, targetBits);
        }

        tinyecm(f1a, small, B1, B1 * 25, curves, 0);

        if ((mpz_sizeinbase(small, 2) > 32) || (mpz_get_ui(small) > rb->lp_cutoff_a))
            return;

        mpz_tdiv_q(large, f1a, small);
        lp_a[num_a++] = mpz_get_ui(small);

        uint64 f64 = do_uecm(mpz_get_ui(large));

        if (f64 <= 1 || f64 > rb->lp_cutoff_a)
            return;

        lp_a[num_a++] = f64;
        mpz_tdiv_q_ui(large, large, f64);

        if ((mpz_sizeinbase(large, 2) > 32) || (mpz_get_ui(large) > rb->lp_cutoff_a))
            return;

        lp_a[num_a++] = mpz_get_ui(large);
    }

	/* yay! Another relation found */

	rb->num_success++;
    // it is complicated to save into the siqs relation buffer
    // from here, so we just mark this cofactor as a success and
    // calling code will put it into the buffer.
    //printf("batch index %d success %u, %u factors: %u,%u,%u\n", 
    //    index, rb->num_success, num_r, lp_r[0], lp_r[1], lp_r[2]);

    for (i = 0; i < MIN(3, num_r); i++)
    {
        c->lp_r[i] = lp_r[i];
    }
    c->success = num_r;

	//rb->print_relation(rb->savefile, c->a, c->b,
	//		f, c->num_factors_r, lp_r,
	//		f + c->num_factors_r, c->num_factors_a, lp_a);
    //buffer_relation(c->offset, lp_r, c->num_factors_r, f, 
    //    c->a, c->b, c->parity, NULL, NULL, 0, 1);

    //printf("success\n");
    return;
}

/*------------------------------------------------------------------*/
void compute_remainder_tree(uint32 first, uint32 last,
				relation_batch_t *rb,
				mpz_t numerator) {

	/* recursively compute numerator % (each relation in 
	   rb->relations) */

	uint32 mid = (first + last) / 2;
	mpz_t relation_prod, remainder;

	/* recursion base case: numerator already fits in
	   an mp_t, so manually compute the remainder and
	   postprocess each relation */

	if (mpz_sizeinbase(numerator, 2) <= (MAX_MP_WORDS * 32)) {
		if (mpz_sgn(numerator) > 0) {
			while (first <= last)
                check_batch_relation(rb, first++, numerator);
		}
		return;
	}

	/* multiply together the unfactored parts of all the
	   relations from first to last */

	mpz_init(relation_prod);
	mpz_init(remainder);
	multiply_relations(first, last, rb, relation_prod);

	/* use the remainder to deal with the left and right
	   halves of the relation list */

	if (mpz_cmp(numerator, relation_prod) < 0) {
		mpz_clear(relation_prod);
		compute_remainder_tree(first, mid, rb, numerator);
		compute_remainder_tree(mid + 1, last, rb, numerator);
	}
	else {
		mpz_tdiv_r(remainder, numerator, relation_prod);
		mpz_clear(relation_prod);
		compute_remainder_tree(first, mid, rb, remainder);
		compute_remainder_tree(mid + 1, last, rb, remainder);
	}
	mpz_clear(remainder);
}

/*------------------------------------------------------------------*/
void relation_batch_init(FILE *logfile, relation_batch_t *rb,
			uint32 min_prime, uint32 max_prime,
			uint32 lp_cutoff_r, uint32 lp_cutoff_a, 
			qs_savefile_t *savefile,
			print_relation_t print_relation,
            int do_prime_product) {

    mpz_init(rb->prime_product);

    if (do_prime_product)
    {
        prime_sieve_t sieve;
        uint32 num_primes, p;

        /* count the number of primes to multiply. Knowing this
           in advance makes the recursion a lot easier, at the cost
           of a small penalty in runtime */

        init_prime_sieve(&sieve, min_prime + 1, max_prime);
        p = min_prime;
        num_primes = 0;
        while (p < max_prime) {
        	p = get_next_prime(&sieve);
        	num_primes++;
        }
        free_prime_sieve(&sieve);

        ///* compute the product of primes */

        logprint(logfile, "multiplying %u primes from %u to %u\n",
        		num_primes, min_prime, max_prime);

        init_prime_sieve(&sieve, min_prime, max_prime);
        multiply_primes(0, num_primes - 2, &sieve, rb->prime_product);
        free_prime_sieve(&sieve);

        logprint(logfile, "multiply complete, product has %u bits\n",
        			(uint32)mpz_sizeinbase(rb->prime_product, 2));
    }
					
	rb->savefile = savefile;
	rb->print_relation = print_relation;
    rb->conversion_ratio = 0.0;

	/* compute the cutoffs used by the recursion base-case. Large
	   primes have a maximum size specified as input arguments, 
	   but numbers that can be passed to the SQUFOF routine are
	   limited to size 2^62 */

	rb->lp_cutoff_r = lp_cutoff_r;
	lp_cutoff_r = MIN(lp_cutoff_r, 0x7fffffff);
	mpz_init(rb->lp_cutoff_r2);
    mpz_set_ui(rb->lp_cutoff_r2, lp_cutoff_r);
    mpz_mul_ui(rb->lp_cutoff_r2, rb->lp_cutoff_r2, lp_cutoff_r);

	rb->lp_cutoff_a = lp_cutoff_a;
	lp_cutoff_a = MIN(lp_cutoff_a, 0x7fffffff);
    mpz_init(rb->lp_cutoff_a2);
    mpz_set_ui(rb->lp_cutoff_a2, lp_cutoff_a);
    mpz_mul_ui(rb->lp_cutoff_a2, rb->lp_cutoff_a2, lp_cutoff_a);

    mpz_init(rb->max_prime2);
    mpz_set_ui(rb->max_prime2, max_prime);
    mpz_mul_ui(rb->max_prime2, rb->max_prime2, max_prime);

	/* allocate lists for relations and their factors */

	rb->target_relations = 500000;
	rb->num_relations = 0;
	rb->num_relations_alloc = 1000;
	rb->relations = (cofactor_t *)xmalloc(rb->num_relations_alloc *
						sizeof(cofactor_t));

	rb->num_factors = 0;
	rb->num_factors_alloc = 10000;
	rb->factors = (uint32 *)xmalloc(rb->num_factors_alloc *
						sizeof(uint32));

    mpz_init(rb->n);
    mpz_init(rb->f1r);
    mpz_init(rb->f2r);
    mpz_init(rb->f1a);
    mpz_init(rb->f2a);
    mpz_init(rb->t0);
    mpz_init(rb->t1);
    mpz_init(rb->small);
    mpz_init(rb->large);
}

/*------------------------------------------------------------------*/
void relation_batch_free(relation_batch_t *rb) {

	mpz_clear(rb->prime_product);
    mpz_clear(rb->max_prime2);
    mpz_clear(rb->lp_cutoff_a2);
    mpz_clear(rb->lp_cutoff_r2);
	free(rb->relations);
	free(rb->factors);

    mpz_clear(rb->n);
    mpz_clear(rb->f1r);
    mpz_clear(rb->f2r);
    mpz_clear(rb->f1a);
    mpz_clear(rb->f2a);
    mpz_clear(rb->t0);
    mpz_clear(rb->t1);
    mpz_clear(rb->small);
    mpz_clear(rb->large);
}

/*------------------------------------------------------------------*/
void relation_batch_add(uint32 a, uint32 b, int32 offset,
			uint32 *factors_r, uint32 num_factors_r, 
			mpz_t unfactored_r_in,
			uint32 *factors_a, uint32 num_factors_a, 
			mpz_t unfactored_a_in,
            mpz_t tmp_in,
			relation_batch_t *rb) {

	uint32 i;
	uint32 *f;
	cofactor_t *c;

	/* add one relation to the batch */
	if (rb->num_relations == rb->num_relations_alloc) {
		rb->num_relations_alloc *= 2;
		rb->relations = (cofactor_t *)xrealloc(rb->relations,
					rb->num_relations_alloc *
					sizeof(cofactor_t));
	}
	c = rb->relations + rb->num_relations++;
	c->a = a;
	c->b = b;
    c->signed_offset = offset;
	c->num_factors_r = num_factors_r;
	c->num_factors_a = num_factors_a;
	c->lp_r_num_words = 0;
	c->lp_a_num_words = 0;
	c->factor_list_word = rb->num_factors;
    c->lp_r[0] = c->lp_r[1] = c->lp_r[2] = 0;
    c->success = 0;

	/* add its small factors */

	if (rb->num_factors + num_factors_r + num_factors_a +
			2 * MAX_LARGE_PRIMES >= rb->num_factors_alloc) {
		rb->num_factors_alloc *= 2;
		rb->factors = (uint32 *)xrealloc(rb->factors,
					rb->num_factors_alloc *
					sizeof(uint32));
	}
	f = rb->factors + rb->num_factors;

	for (i = 0; i < num_factors_r; i++)
		f[i] = factors_r[i];
	f += i;
	rb->num_factors += i;

	for (i = 0; i < num_factors_a; i++)
		f[i] = factors_a[i];
	f += i;
	rb->num_factors += i;

	/* add its large factors (if any) */

	if (mpz_cmp_ui(unfactored_r_in, 1) > 0) {

        mpz_set(rb->t0, unfactored_r_in);
        i = 0;
        while (mpz_cmp_ui(rb->t0, 0) > 0)
        {
            f[i++] = (uint32)mpz_get_ui(rb->t0);
            mpz_tdiv_q_2exp(rb->t0, rb->t0, 32);
        }

		f += i;
		rb->num_factors += i;
		c->lp_r_num_words = i;
	}

    if (mpz_cmp_ui(unfactored_a_in, 1) > 0) {
        mpz_set(rb->t0, unfactored_a_in);
        i = 0;
        while (mpz_cmp_ui(rb->t0, 0) > 0)
        {
            f[i++] = (uint32)mpz_get_ui(rb->t0);
            mpz_tdiv_q_2exp(rb->t0, rb->t0, 32);
        }
    
        f += i;
        rb->num_factors += i;
        c->lp_a_num_words = i;
    }
}
	
/*------------------------------------------------------------------*/
uint32 relation_batch_run(relation_batch_t *rb) {

	rb->num_success = 0;
	if (rb->num_relations > 0) {
		compute_remainder_tree(0, rb->num_relations - 1,
					rb, rb->prime_product);
	}

	/* wipe out batched relations */

	//rb->num_relations = 0;
	//rb->num_factors = 0;
	return rb->num_success;
}


