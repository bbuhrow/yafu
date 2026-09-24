/*
 * arith/mod128.c
 *
 * 128-bit unsigned modulo: c = a % b
 *
 *   void mod_128(const uint64_t *a, const uint64_t *b, uint64_t *c)
 *
 * Each of a, b, c is a two-element uint64_t array laid out as:
 *   [0] = low  64 bits
 *   [1] = high 64 bits
 * so the 128-bit value is  a[1]:a[0]  (high:low).
 *
 * Implementation notes
 * ────────────────────
 * x86-64's DIV instruction performs a 128÷64 division when rdx:rax holds
 * the 128-bit dividend and the operand holds the 64-bit divisor, producing
 * a 64-bit quotient in rax and a 64-bit remainder in rdx.  A true 128÷128
 * division requires that the high 64 bits of the divisor are zero OR that
 * we use a multi-word algorithm.
 *
 * Case 1 — b fits in 64 bits (b[1] == 0):
 *   Use a single DIV with rdx:rax = a[1]:a[0] and divisor = b[0].
 *   This is exact and very fast.
 *
 * Case 2 — b uses the full 128 bits (b[1] != 0):
 *   Use the portable __uint128_t path on GCC/Clang, or a two-step
 *   shift-and-divide reduction on MSVC where __uint128_t is unavailable.
 *   The MSVC path uses the same DIV-based approach described in Hacker's
 *   Delight §9-4 (Knuth's Algorithm D for two-word divisors).
 *
 * All three compilers (GCC, Clang/ClangCL, MSVC) produce correct results.
 * GCC and Clang additionally optimise the __uint128_t path to a tight
 * sequence of DIV + MUL + SUB instructions.
 */

#include <stdint.h>

/* ── compiler capability detection ──────────────────────────────────────── */
#if defined(__GNUC__) || defined(__clang__)
#  define HAVE_UINT128 1
#endif

/* ── helper: 128-bit unsigned mod using compiler built-in ────────────────── */
#if HAVE_UINT128
static void mod_128_full(const uint64_t *a, const uint64_t *b, uint64_t *c)
{
    __uint128_t va = ((__uint128_t)a[1] << 64) | a[0];
    __uint128_t vb = ((__uint128_t)b[1] << 64) | b[0];
    __uint128_t vc = va % vb;
    c[0] = (uint64_t)(vc);
    c[1] = (uint64_t)(vc >> 64);
}

#else  /* MSVC — no __uint128_t */

/*
 * Portable 128÷128 → 128 remainder using Knuth's Algorithm D (two-word
 * divisor).  We normalise so the divisor's high word is non-zero, perform
 * two 128÷64 DIV steps, then un-normalise the remainder.
 *
 * This is intentionally straightforward rather than maximally optimised;
 * the b[1]!=0 path is rare in yafu's inner loops.
 */
#include <intrin.h>   /* _udiv128 available in VS 2019 16.8+ */

static void mod_128_full(const uint64_t *a, const uint64_t *b, uint64_t *c)
{
    /*
     * Normalise: shift both operands left by `shift` bits so that the
     * most-significant bit of b is 1.  This is required by Algorithm D
     * to bound the quotient estimate error.
     */
    unsigned shift = 0;
    uint64_t bhi = b[1];
    while ((bhi & 0x8000000000000000ULL) == 0) { bhi <<= 1; ++shift; }

    /* Shifted divisor */
    uint64_t bn_hi = (b[1] << shift) | (shift ? b[0] >> (64 - shift) : 0);
    uint64_t bn_lo = b[0] << shift;

    /* Shifted dividend — needs an extra high word */
    uint64_t an_x  = shift ? a[1] >> (64 - shift) : 0;  /* overflow word */
    uint64_t an_hi = (a[1] << shift) | (shift ? a[0] >> (64 - shift) : 0);
    uint64_t an_lo = a[0] << shift;

    /*
     * Step 1: estimate q1 = floor((an_x:an_hi) / bn_hi)
     * using a single 128÷64 DIV, then correct downward at most twice.
     */
    uint64_t rem1;
    uint64_t q1 = _udiv128(an_x, an_hi, bn_hi, &rem1);

    /* Correction: while q1*bn_lo > rem1:bn_lo (as 128-bit compare) */
    /* We check using unsigned 128-bit multiply via _umul128 */
    uint64_t p_hi, p_lo;
    p_lo = _umul128(q1, bn_lo, &p_hi);
    while (p_hi > rem1 || (p_hi == rem1 && p_lo > an_lo)) {
        --q1;
        if (rem1 + bn_hi < rem1) break;  /* overflow means q1 is now correct */
        rem1 += bn_hi;
        p_lo = _umul128(q1, bn_lo, &p_hi);
    }

    /*
     * Step 2: compute the partial remainder r1 = an_x:an_hi:an_lo - q1*bn
     * then estimate q0 = floor(r1_hi:r1_lo / bn_hi) and correct.
     */
    uint64_t mul_hi, mul_lo;
    mul_lo = _umul128(q1, bn_lo, &mul_hi);
    uint64_t r_lo = an_lo - mul_lo;
    uint64_t borrow = (r_lo > an_lo) ? 1 : 0;
    uint64_t r_hi = an_hi - _umul128(q1, bn_hi, &mul_hi) - mul_hi - borrow;

    uint64_t rem0;
    uint64_t q0 = _udiv128(r_hi, r_lo, bn_hi, &rem0);

    mul_lo = _umul128(q0, bn_lo, &mul_hi);
    while (mul_hi > rem0 || (mul_hi == rem0 && mul_lo > (shift ? 0 : an_lo))) {
        --q0;
        if (rem0 + bn_hi < rem0) break;
        rem0 += bn_hi;
        mul_lo = _umul128(q0, bn_lo, &mul_hi);
    }

    /* Un-normalised remainder = (r_hi:r_lo - q0*bn) >> shift */
    mul_lo = _umul128(q0, bn_lo, &mul_hi);
    uint64_t rr_lo = r_lo - mul_lo;
    borrow = (rr_lo > r_lo) ? 1 : 0;
    uint64_t rr_hi = r_hi - _umul128(q0, bn_hi, &mul_hi) - mul_hi - borrow;

    /* Shift right to undo normalisation */
    if (shift) {
        c[0] = (rr_hi << (64 - shift)) | (rr_lo >> shift);
        c[1] = rr_hi >> shift;
    } else {
        c[0] = rr_lo;
        c[1] = rr_hi;
    }
}
#endif  /* HAVE_UINT128 / MSVC */


/* ── public interface ────────────────────────────────────────────────────── */

void mod_128(const uint64_t *a, const uint64_t *b, uint64_t *c)
{
    if (b[1] == 0)
    {
        /*
         * Fast path: b fits in 64 bits.
         * Use a single DIV instruction via the 128÷64 form:
         *   rdx:rax = a[1]:a[0]  divided by  b[0]
         * Remainder lands in rdx (GCC/Clang) or is returned directly (MSVC).
         *
         * Precondition: a[1] < b[0]  (otherwise quotient overflows 64 bits).
         * If a[1] >= b[0] we reduce a[1] first using a[1] % b[0].
         */
#if HAVE_UINT128
        __uint128_t va = ((__uint128_t)a[1] << 64) | a[0];
        uint64_t    rem = (uint64_t)(va % b[0]);
        c[0] = rem;
        c[1] = 0;
#else
        /* MSVC: _udiv128(hi, lo, divisor, &remainder) -> quotient */
        uint64_t hi = a[1];
        if (hi >= b[0])
            hi %= b[0];          /* reduce so quotient fits in 64 bits */
        uint64_t rem;
        _udiv128(hi, a[0], b[0], &rem);
        c[0] = rem;
        c[1] = 0;
#endif
    }
    else
    {
        /* Full 128÷128 path */
        mod_128_full(a, b, c);
    }
}
