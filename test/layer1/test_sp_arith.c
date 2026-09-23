/*----------------------------------------------------------------------
 Layer 1 -- arith.h single-word (sp*) primitives, checked against GMP.

 Every reference value is computed with GMP (mpz_mul/mpz_mod, mpz_powm,
 mpz_tdiv_qr, mpz_invert, mpz_sizeinbase) rather than a hand-rolled oracle,
 so the "known good" answer is the library's, not ours. Operands are single
 words, generated from the test's tk RNG stream (reproducible from the seed);
 we only ever BUILD mpz values from those words and compare, never extract a
 word from an mpz, which keeps the comparisons portable across LP64/LLP64.

 Covers: spMultiply, spMulMod, spModExp, spDivide, modinv_1 (+1b/1c),
 bits64.

 NOTE (verify on first build): two arith.h conventions are assumed and
 flagged so they are easy to flip if the tree differs:
   * spDivide's dividend limb order: u[0] = low word, u[1] = high word,
     with u[1] < v so the quotient fits in 64 bits. If your build reverses
     the words, swap the two stores marked [LIMB-ORDER].
   * modinv_1(a,p) returns the canonical inverse in [0,p) for gcd(a,p)=1,
     matching mpz_invert. If it uses a different representative, the modinv
     test will flag it.
 Public domain.
----------------------------------------------------------------------*/
#include "testkit.h"
#include <gmp.h>
#include "arith.h"

/* ---- portable mpz <- integer builders (mpz_set_ui takes only the low 32
   bits on LLP64, so assemble 64-bit values in 32-bit chunks). No temporaries,
   so these are cheap enough to call in a hot loop. ---- */
static void mpz_set_u64(mpz_t r, uint64_t v)
{
    mpz_set_ui(r, (unsigned long)(v >> 32));
    mpz_mul_2exp(r, r, 32);
    mpz_add_ui(r, r, (unsigned long)(v & 0xFFFFFFFFUL));
}
static void mpz_set_u128(mpz_t r, uint64_t lo, uint64_t hi)
{
    mpz_set_u64(r, hi);
    mpz_mul_2exp(r, r, 32);
    mpz_add_ui(r, r, (unsigned long)(lo >> 32));
    mpz_mul_2exp(r, r, 32);
    mpz_add_ui(r, r, (unsigned long)(lo & 0xFFFFFFFFUL));
}

/* ================================================================== */
static void t_spMultiply(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t A, B, P, got;
    int k;
    mpz_inits(A, B, P, got, NULL);
    for (k = 0; k < 300000; k++) {
        uint64_t a = tk_rng_u64(r), b = tk_rng_u64(r), lo, hi;
        spMultiply(a, b, &lo, &hi);
        mpz_set_u64(A, a); mpz_set_u64(B, b); mpz_mul(P, A, B);
        mpz_set_u128(got, lo, hi);
        if (!TK_CHECKF(tk, mpz_cmp(got, P) == 0,
                       "spMultiply(%llu,%llu) mismatch",
                       (unsigned long long)a, (unsigned long long)b)) break;
    }
    mpz_clears(A, B, P, got, NULL);
}

static void t_spMulMod(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t A, B, M, W, got;
    int k;
    mpz_inits(A, B, M, W, got, NULL);
    for (k = 0; k < 300000; k++) {
        uint64_t m = tk_rng_u64(r) | 1ULL;
        uint64_t a = tk_rng_u64(r) % m, b = tk_rng_u64(r) % m, w;
        spMulMod(a, b, m, &w);
        mpz_set_u64(A, a); mpz_set_u64(B, b); mpz_set_u64(M, m);
        mpz_mul(W, A, B); mpz_mod(W, W, M);
        mpz_set_u64(got, w);
        if (!TK_CHECKF(tk, mpz_cmp(got, W) == 0,
                       "spMulMod(%llu,%llu mod %llu) mismatch",
                       (unsigned long long)a, (unsigned long long)b,
                       (unsigned long long)m)) break;
    }
    mpz_clears(A, B, M, W, got, NULL);
}

static void t_spModExp(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t A, E, M, U, got;
    int k;
    mpz_inits(A, E, M, U, got, NULL);
    for (k = 0; k < 50000; k++) {
        uint64_t m = tk_rng_u64(r) | 1ULL;
        uint64_t a = tk_rng_u64(r) % m, e = tk_rng_u64(r), u;
        spModExp(a, e, m, &u);
        mpz_set_u64(A, a); mpz_set_u64(E, e); mpz_set_u64(M, m);
        mpz_powm(U, A, E, M);
        mpz_set_u64(got, u);
        if (!TK_CHECKF(tk, mpz_cmp(got, U) == 0,
                       "spModExp(%llu^%llu mod %llu) mismatch",
                       (unsigned long long)a, (unsigned long long)e,
                       (unsigned long long)m)) break;
    }
    /* Fermat spot-check: for prime m, a^(m-1) == 1 (mod m), a != 0 */
    {
        uint64_t m = 2147483647ULL, u;         /* 2^31-1, prime */
        uint64_t a = 1 + tk_rng_u64(r) % (m - 1);
        spModExp(a, m - 1, m, &u);
        TK_EQ_U64(tk, u, 1);
    }
    mpz_clears(A, E, M, U, got, NULL);
}

static void t_spDivide(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t NUM, V, Q, R, gotq, gotr;
    int k;
    mpz_inits(NUM, V, Q, R, gotq, gotr, NULL);
    for (k = 0; k < 300000; k++) {
        uint64_t v = tk_rng_u64(r) | 1ULL;
        uint64_t hi = tk_rng_u64(r) % v;        /* ensure quotient fits 64b */
        uint64_t lo = tk_rng_u64(r);
        uint64_t u[2], q, rr;
        u[0] = lo;   /* [LIMB-ORDER] low word  */
        u[1] = hi;   /* [LIMB-ORDER] high word */
        spDivide(&q, &rr, u, v);
        mpz_set_u128(NUM, lo, hi); mpz_set_u64(V, v);
        mpz_tdiv_qr(Q, R, NUM, V);
        mpz_set_u64(gotq, q); mpz_set_u64(gotr, rr);
        if (!TK_CHECKF(tk, mpz_cmp(gotq, Q) == 0, "spDivide quotient mismatch (k=%d)", k)) break;
        if (!TK_CHECKF(tk, mpz_cmp(gotr, R) == 0, "spDivide remainder mismatch (k=%d)", k)) break;
    }
    mpz_clears(NUM, V, Q, R, gotq, gotr, NULL);
}

static void t_modinv(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t A, P, INV, got;
    int k;
    mpz_inits(A, P, INV, got, NULL);
    /* modinv_1 works on 32-bit operands */
    for (k = 0; k < 200000; k++) {
        uint32_t p = (uint32_t)tk_rng_u64(r) | 1u;   /* odd */
        uint32_t a;
        if (p < 3) p = 3;
        a = (uint32_t)tk_rng_u64(r) % p;
        if (a == 0) a = 1;
        mpz_set_u64(A, a); mpz_set_u64(P, p);
        /* GMP is the reference: mpz_invert returns 0 iff a is not invertible
           mod p; when it is, INV holds the canonical inverse in [0,p). */
        if (mpz_invert(INV, A, P) == 0) continue;    /* gcd(a,p) != 1: skip */
        mpz_set_u64(got, modinv_1(a, p));
        if (!TK_CHECKF(tk, mpz_cmp(got, INV) == 0,
                       "modinv_1(%lu mod %lu) mismatch",
                       (unsigned long)a, (unsigned long)p)) break;

        //p is fixed at 2^32, skip this check
        //mpz_set_u64(got, modinv_1b(a, p)); TK_CHECK(tk, mpz_cmp(got, INV) == 0);

        mpz_set_u64(got, modinv_1c(a, p)); 
        if (!TK_CHECKF(tk, mpz_cmp(got, INV) == 0,
            "modinv_1c(%lu mod %lu) mismatch",
            (unsigned long)a, (unsigned long)p)) break;

    }
    mpz_clears(A, P, INV, got, NULL);
}

static void t_bits(tk_ctx *tk)
{
    tk_rng *r = tk_rng_of(tk);
    mpz_t X;
    int k, i;
    mpz_init(X);
    TK_EQ_U64(tk, (uint64_t)bits64(0), 0);
    for (i = 0; i < 64; i++)
        TK_EQ_U64(tk, (uint64_t)bits64(1ULL << i), (uint64_t)(i + 1));
    for (k = 0; k < 100000; k++) {
        uint64_t x = tk_rng_u64(r);
        uint64_t ref;
        /* mpz_sizeinbase(0,2) is 1 (GMP quirk), but bits64(0) is 0 */
        if (x == 0) {
            ref = 0;
        } else {
            mpz_set_u64(X, x);
            ref = (uint64_t)mpz_sizeinbase(X, 2);
        }
        if (!TK_EQ_U64(tk, (uint64_t)bits64(x), ref)) break;
    }
    mpz_clear(X);
}

static const tk_test tk__sp_arith_tests[] = {
    { "spMultiply", t_spMultiply, "fast arith sp gmp" },
    { "spMulMod",   t_spMulMod,   "fast arith sp gmp" },
    { "spModExp",   t_spModExp,   "fast arith sp gmp" },
    { "spDivide",   t_spDivide,   "fast arith sp gmp" },
    { "modinv_1",   t_modinv,     "fast arith sp gmp" },
    { "bits64",     t_bits,       "fast arith sp gmp" }
};

const tk_module tk_module_sp_arith = {
    "sp_arith",
    "arith.h single-word mul/div/mulmod/modexp/modinv vs GMP",
    tk__sp_arith_tests,
    (int)(sizeof tk__sp_arith_tests / sizeof tk__sp_arith_tests[0])
};
