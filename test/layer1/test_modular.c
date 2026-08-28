/*----------------------------------------------------------------------
 Layer 1 -- Montgomery modular arithmetic (monty.h) vs GMP.

 Two interfaces are exercised:
   * The mpz Montgomery context (monty_alloc/monty_init/to_monty/monty_mul/
     monty_add/monty_sub/monty_redc) over odd moduli, checked against GMP.
   * The fixed 128-bit interface (monty128_init/to_monty128/mulmod128/
     sqrmod128/addmod128/submod128) over odd 128-bit moduli, checked
     against GMP.

 Montgomery form recap (R = 2^k > n, gcd(R,n)=1, so n must be ODD):
     to_monty(x)      -> x*R mod n
     monty_mul(X,Y)   -> x*y*R mod n      (stays in the Montgomery domain)
     monty_redc(Z)    -> Z*R^{-1} mod n   (leaves the domain)
 so redc(mul(to(x),to(y))) == x*y mod n, which is what we assert.
 add/sub are domain-linear, so redc(add(to(x),to(y))) == (x+y) mod n.

 For the 128-bit interface there is no explicit "from domain" call, so we
 leave the domain by Montgomery-multiplying by literal one:
     mulmod128(Z, {1,0}) == Z*R^{-1} mod n.

 NOTE (verify on first build): monty128's decls live in monty.h in the
 public tree; in the refactored tree they may have moved to limb2.h -- if
 so, change the second include below. Also confirmed against monty.h:
 to_monty/to_monty128 and monty_redc act IN PLACE on their mpz/uint64[2]
 argument; monty_mul/add/sub write to a separate res.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include "test_data.h"
#include <gmp.h>
#include "monty.h"          /* if monty128* moved to limb2.h in your tree, include that too */

/* set an mpz from a uint64[2] little-endian pair (lo, hi) */
static void set_mpz_u128(mpz_t out, const uint64_t v[2])
{
    mpz_set_ui(out, v[1]);
    mpz_mul_2exp(out, out, 64);
    mpz_add_ui(out, out, v[0]);
}
/* read an mpz (assumed < 2^128) back into a uint64[2] little-endian pair */
static void get_mpz_u128(uint64_t v[2], const mpz_t in)
{
    mpz_t t; mpz_init(t);
    v[0] = mpz_get_ui(in) & 0xFFFFFFFFULL;               /* build 64b safely */
    mpz_tdiv_q_2exp(t, in, 32); v[0] |= (mpz_get_ui(t) & 0xFFFFFFFFULL) << 32;
    mpz_tdiv_q_2exp(t, in, 64);
    v[1] = mpz_get_ui(t) & 0xFFFFFFFFULL;
    mpz_tdiv_q_2exp(t, in, 96); v[1] |= (mpz_get_ui(t) & 0xFFFFFFFFULL) << 32;
    mpz_clear(t);
}

/* ================================================================== *
 * mpz Montgomery vs GMP
 * ================================================================== */
static void t_monty_mpz(tk_ctx *tk)
{
    tk_rng *rng = tk_rng_of(tk);
    monty_t *md = monty_alloc();
    mpz_t n, x, y, X, Y, Z, ref, want;
    gmp_randstate_t gmp_rng_state;
    int k;

    mpz_inits(n, x, y, X, Y, Z, ref, want, NULL);
    /* GMP's own randomness, seeded from this test's tk stream so runs remain
       reproducible (for a given GMP version) from the printed seed. */
    gmp_randinit_default(gmp_rng_state);
    gmp_randseed_ui(gmp_rng_state, (unsigned long)tk_rng_u64(rng));

    for (k = 0; k < 2000; k++) {
        int bits = 40 + (int)tk_rng_range(rng, 200);   /* up to 239-bit odd n */
        mpz_urandomb(n, gmp_rng_state, (mp_bitcnt_t)bits);
        mpz_setbit(n, 0);                              /* force odd */
        if (mpz_cmp_ui(n, 3) < 0) mpz_set_ui(n, 3);    /* enforce minimum */

        /* random residues x, y < n */
        mpz_urandomm(x, gmp_rng_state, n);
        mpz_urandomm(y, gmp_rng_state, n);

        monty_init(n, md);

        mpz_set(X, x); to_monty(md, X);
        mpz_set(Y, y); to_monty(md, Y);

        /* multiply */
        monty_mul(md, X, Y, Z);
        mpz_set(ref, Z); monty_redc(md, ref);
        mpz_mul(want, x, y); mpz_mod(want, want, n);
        TK_CHECKF(tk, mpz_cmp(ref, want) == 0, "monty_mul mismatch (k=%d)", k);

        /* add */
        monty_add(md, X, Y, Z);
        mpz_set(ref, Z); monty_redc(md, ref);
        mpz_add(want, x, y); mpz_mod(want, want, n);
        TK_CHECKF(tk, mpz_cmp(ref, want) == 0, "monty_add mismatch (k=%d)", k);

        /* sub */
        monty_sub(md, X, Y, Z);
        mpz_set(ref, Z); monty_redc(md, ref);
        mpz_sub(want, x, y); mpz_mod(want, want, n);   /* mpz_mod is non-negative */
        TK_CHECKF(tk, mpz_cmp(ref, want) == 0, "monty_sub mismatch (k=%d)", k);
    }

    gmp_randclear(gmp_rng_state);
    mpz_clears(n, x, y, X, Y, Z, ref, want, NULL);
    monty_free(md);
}

/* ================================================================== *
 * 128-bit Montgomery vs GMP
 * ================================================================== */
static void t_monty128(tk_ctx *tk)
{
    tk_rng *rng = tk_rng_of(tk);
    mpz_t N, mx, my, mprod, msum, mdiff, mref;
    gmp_randstate_t gmp_rng_state;
    int k;

    mpz_inits(N, mx, my, mprod, msum, mdiff, mref, NULL);
    gmp_randinit_default(gmp_rng_state);
    gmp_randseed_ui(gmp_rng_state, (unsigned long)tk_rng_u64(rng));

    for (k = 0; k < 20000; k++) {
        monty128_t md;
        uint64_t n[2], x[2], y[2], X[2], Y[2], Z[2], out[2];
        uint64_t one[2] = { 1ULL, 0ULL };

        mpz_urandomb(N, gmp_rng_state, 128);
        mpz_setbit(N, 0);                           /* force odd */
        if (mpz_cmp_ui(N, 3) < 0) mpz_set_ui(N, 3); /* enforce minimum */
        get_mpz_u128(n, N);

        /* residues < n */
        mpz_urandomm(mx, gmp_rng_state, N); get_mpz_u128(x, mx);
        mpz_urandomm(my, gmp_rng_state, N); get_mpz_u128(y, my);

        monty128_init(&md, n);

        /* --- mulmod128 --- */
        X[0]=x[0]; X[1]=x[1]; to_monty128(&md, X);
        Y[0]=y[0]; Y[1]=y[1]; to_monty128(&md, Y);
        mulmod128(X, Y, Z, &md);
        mulmod128(Z, one, out, &md);      /* leave the domain */
        set_mpz_u128(mref, out);
        mpz_mul(mprod, mx, my); mpz_mod(mprod, mprod, N);
        TK_CHECKF(tk, mpz_cmp(mref, mprod) == 0, "mulmod128 mismatch (k=%d)", k);

        /* --- sqrmod128 --- */
        sqrmod128(X, Z, &md);
        mulmod128(Z, one, out, &md);
        set_mpz_u128(mref, out);
        mpz_mul(mprod, mx, mx); mpz_mod(mprod, mprod, N);
        TK_CHECKF(tk, mpz_cmp(mref, mprod) == 0, "sqrmod128 mismatch (k=%d)", k);

        /* --- addmod128 / submod128 (plain modular, not domain) --- */
        addmod128(x, y, out, n);
        set_mpz_u128(mref, out);
        mpz_add(msum, mx, my); mpz_mod(msum, msum, N);
        TK_CHECKF(tk, mpz_cmp(mref, msum) == 0, "addmod128 mismatch (k=%d)", k);

        submod128(x, y, out, n);
        set_mpz_u128(mref, out);
        mpz_sub(mdiff, mx, my); mpz_mod(mdiff, mdiff, N);
        TK_CHECKF(tk, mpz_cmp(mref, mdiff) == 0, "submod128 mismatch (k=%d)", k);
    }

    gmp_randclear(gmp_rng_state);
    mpz_clears(N, mx, my, mprod, msum, mdiff, mref, NULL);
}

static const tk_test tk__modular_tests[] = {
    { "monty_mpz",  t_monty_mpz,  "slow gmp modular" },
    { "monty128",   t_monty128,   "slow gmp modular" }
};

const tk_module tk_module_modular = {
    "modular",
    "Montgomery arithmetic (mpz and 128-bit) vs GMP",
    tk__modular_tests,
    (int)(sizeof tk__modular_tests / sizeof tk__modular_tests[0])
};
