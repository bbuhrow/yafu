/*----------------------------------------------------------------------
limb_selftest.c -- stand-alone self-test and micro-benchmark for the
one/two-limb library (mp_platform.h, limb1.{h,c}, limb2.{h,c}).

Needs nothing but the library and libc -- no GMP, no YAFU -- so it is the
quickest way to check a new compiler / OS / CPU:

    gcc   -O2 [-march=native] limb_selftest.c limb1.c limb2.c -lm -o limb_selftest
    clang -O2 [-march=native] limb_selftest.c limb1.c limb2.c -lm -o limb_selftest
    cl    /O2 [/arch:AVX512]  limb_selftest.c limb1.c limb2.c
    clang-cl /O2 [-march=native] limb_selftest.c limb1.c limb2.c

    ./limb_selftest            run the checks (exit status 0 = all passed)
    ./limb_selftest bench      also time the main kernels

Every routine is checked against an independent, deliberately slow
reference built from bit-at-a-time shifts, adds and compares.
Public domain.
----------------------------------------------------------------------*/
#include "limb1.h"
#include "limb2.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

static int nfail = 0;
static long nchecks = 0;

static void fail(const char *name, const char *file, int line)
{
    if (nfail < 20)
        printf("  FAIL %s (%s:%d)\n", name, file, line);
    nfail++;
}

#define CHECK(cond, name) \
    do { nchecks++; if (!(cond)) { fail(name, __FILE__, __LINE__); } } while (0)

/* ---- xorshift RNG ---- */
static uint64_t rs = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void)
{
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
    return rs;
}
/* random value of a random bit length in [1, maxbits] */
static uint64_t rndbits(int maxbits)
{
    int b = 1 + (int)(rnd() % (uint64_t)maxbits);
    uint64_t x = rnd() >> (64 - b);
    return x | ((uint64_t)1 << (b - 1));
}

/* ---- slow references on little-endian limb arrays ---- */
typedef struct { uint64_t w[4]; } u256;

static int ref_cmp(const uint64_t *a, const uint64_t *b, int n)
{
    int i;
    for (i = n - 1; i >= 0; i--)
        if (a[i] != b[i]) return (a[i] < b[i]) ? -1 : 1;
    return 0;
}
static uint64_t ref_add(uint64_t *r, const uint64_t *a, const uint64_t *b, int n)
{
    uint64_t c = 0; int i;
    for (i = 0; i < n; i++) {
        uint64_t s = a[i] + c, c1 = (s < c);
        r[i] = s + b[i]; c = c1 | (r[i] < s);
    }
    return c;
}
static uint64_t ref_sub(uint64_t *r, const uint64_t *a, const uint64_t *b, int n)
{
    uint64_t c = 0; int i;
    for (i = 0; i < n; i++) {
        uint64_t d = a[i] - c, c1 = (a[i] < c);
        r[i] = d - b[i]; c = c1 | (d < b[i]);
    }
    return c;
}
static uint64_t ref_shl1(uint64_t *a, int n)
{
    uint64_t c = 0; int i;
    for (i = 0; i < n; i++) { uint64_t t = a[i] >> 63; a[i] = (a[i] << 1) | c; c = t; }
    return c;
}
/* 64 x 64 -> 128 from 32-bit pieces */
static void ref_mul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi)
{
    uint64_t al = (uint32_t)a, ah = a >> 32, bl = (uint32_t)b, bh = b >> 32;
    uint64_t ll = al * bl, lh = al * bh, hl = ah * bl, hh = ah * bh;
    uint64_t mid = (ll >> 32) + (uint32_t)lh + (uint32_t)hl;
    *lo = (mid << 32) | (uint32_t)ll;
    *hi = hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
}
/* n-limb x m-limb schoolbook product into r[n+m] */
static void ref_mul(uint64_t *r, const uint64_t *a, int n, const uint64_t *b, int m)
{
    int i, j;
    memset(r, 0, (size_t)(n + m) * sizeof(uint64_t));
    for (i = 0; i < n; i++) {
        uint64_t carry = 0;
        for (j = 0; j < m; j++) {
            uint64_t lo, hi, s;
            ref_mul64(a[i], b[j], &lo, &hi);
            s = r[i + j] + lo; hi += (s < lo);
            r[i + j] = s + carry; hi += (r[i + j] < s);
            carry = hi;
        }
        r[i + m] = carry;
    }
}
/* x[0..n) mod d[0..m), m <= 2, restoring division one bit at a time; the
   remainder comes back in rem[0..2) */
static void ref_mod(uint64_t *rem, const uint64_t *x, int n, const uint64_t *d, int m)
{
    uint64_t r[3] = { 0, 0, 0 }, dd[3] = { 0, 0, 0 };
    int i;
    for (i = 0; i < m; i++) dd[i] = d[i];
    for (i = 64 * n - 1; i >= 0; i--) {
        ref_shl1(r, 3);
        r[0] |= (x[i / 64] >> (i % 64)) & 1;
        if (ref_cmp(r, dd, 3) >= 0) ref_sub(r, r, dd, 3);
    }
    rem[0] = r[0]; rem[1] = r[1];
}

/* ================================================================== */
static void test_primitives(void)
{
    int k, i;
    for (k = 0; k < 200000; k++) {
        uint64_t a = rnd(), b = rnd(), lo, hi, rlo, rhi, d, q, r, q2, r2, t[2], rem[2];
        if ((k & 7) == 0) a = rndbits(64);
        if ((k & 15) == 0) b = ~(uint64_t)0;

        lo = _umul128(a, b, &hi);
        ref_mul64(a, b, &rlo, &rhi);
        CHECK(lo == rlo && hi == rhi, "_umul128");
        CHECK(__umulh(a, b) == rhi, "__umulh");

        d = rndbits(64);
        hi = rnd() % d; lo = rnd();
        q = _udiv128(hi, lo, d, &r);
        q2 = mp_udiv128_soft(hi, lo, d, &r2);
        CHECK(q == q2 && r == r2, "_udiv128 vs soft");
        t[0] = lo; t[1] = hi;
        ref_mod(rem, t, 2, &d, 1);
        CHECK(r == rem[0], "_udiv128 remainder");
        ref_mul64(q, d, &rlo, &rhi);                 /* q*d + r == hi:lo */
        rlo += r; rhi += (rlo < r);
        CHECK(rlo == lo && rhi == hi, "_udiv128 reconstruction");

        {   /* carry chains */
            uint64_t x[4], y[4], z[4], e[4];
            mp_carry_t c = 0;
            for (i = 0; i < 4; i++) { x[i] = rnd(); y[i] = (k & 3) ? rnd() : ~x[i]; }
            for (i = 0; i < 4; i++) c = mp_addcarry_u64(c, x[i], y[i], &z[i]);
            CHECK(c == ref_add(e, x, y, 4) && !ref_cmp(z, e, 4), "mp_addcarry_u64");
            c = 0;
            for (i = 0; i < 4; i++) c = mp_subborrow_u64(c, x[i], y[i], &z[i]);
            CHECK(c == ref_sub(e, x, y, 4) && !ref_cmp(z, e, 4), "mp_subborrow_u64");
        }
        {   /* 32-bit carry ops */
            uint32_t x = (uint32_t)a, y = (uint32_t)b, z;
            mp_carry_t cin = (mp_carry_t)(k & 1), c;
            c = mp_addcarry_u32(cin, x, y, &z);
            CHECK(z == (uint32_t)(x + y + cin) && c == (((uint64_t)x + y + cin) >> 32), "mp_addcarry_u32");
            c = mp_subborrow_u32(cin, x, y, &z);
            CHECK(z == (uint32_t)(x - y - cin) && c == ((uint64_t)x < (uint64_t)y + cin), "mp_subborrow_u32");
        }
        {   /* bit scans, popcount */
            uint64_t x = rndbits(64) << (rnd() % 64);
            int tz = 0, lz = 0, pc = 0;
            if (x == 0) x = 1;
            while (!((x >> tz) & 1)) tz++;
            while (!((x << lz) >> 63)) lz++;
            for (i = 0; i < 64; i++) pc += (int)((x >> i) & 1);
            CHECK((int)_trail_zcnt64(x) == tz, "_trail_zcnt64");
            CHECK((int)_lead_zcnt64(x) == lz, "_lead_zcnt64");
            CHECK((int)_trail_full_zcnt(x) == tz && (int)_lead_full_zcnt(x) == lz, "full zcnt");
            CHECK((int)mp_popcount64(x) == pc, "mp_popcount64");
            CHECK(_reset_lsb64(x) == (x & (x - 1)), "_reset_lsb64");
            if ((uint32_t)x) {
                int t32 = 0;
                while (!(((uint32_t)x >> t32) & 1)) t32++;
                CHECK((int)_trail_zcnt((uint32_t)x) == t32, "_trail_zcnt");
            }
        }
    }
    CHECK(_trail_full_zcnt(0) == 64 && _lead_full_zcnt(0) == 64, "full zcnt at zero");
    CHECK(my_clz128(0, 0) == 128 && my_ctz128(0, 0) == 128, "128-bit scans at zero");
    CHECK(my_clz128(1, 0) == 127 && my_ctz128(0, 1) == 64, "128-bit scans");
}

static void test_limb1(void)
{
    int k;
    for (k = 0; k < 200000; k++) {
        uint64_t n = rndbits(64) | 1, x, y, c, inv, e[2], p[2], t[3], r, g, a, b;
        if (n < 3) n = 3;
        x = rnd() % n; y = rnd() % n; c = rnd() % n;
        inv = multiplicative_inverse(n);
        CHECK(inv * n == 1, "multiplicative_inverse");
        CHECK(((multiplicative_inverse52(n) * n) & 0x000fffffffffffffULL) == 1, "multiplicative_inverse52");

        /* REDC: r = x*y / 2^64 mod n  <=>  r * 2^64 == x*y  (mod n) */
        ref_mul64(x, y, &p[0], &p[1]);
        ref_mod(e, p, 2, &n, 1);                        /* e = x*y mod n */
#define REDC_OK(rr) (t[0] = 0, t[1] = (rr), ref_mod(p, t, 2, &n, 1), (rr) < n && p[0] == e[0])
        r = mulredc_pos(x, y, n, inv);       CHECK(REDC_OK(r), "mulredc_pos");
        r = mulredc_pos_alt(x, y, n, inv);   CHECK(REDC_OK(r), "mulredc_pos_alt");
        r = mulredc(x, y, n, 0 - inv);       CHECK(REDC_OK(r), "mulredc");
        if (n < ((uint64_t)1 << 63)) {
            r = mulredc63(x, y, n, 0 - inv); CHECK(REDC_OK(r), "mulredc63");
        }
        if (n < ((uint64_t)1 << 60)) {                  /* lazy: result < 2n */
            r = mulredc60(x, y, n, 0 - inv);
            if (r >= n) r -= n;
            CHECK(REDC_OK(r), "mulredc60");
        }
        ref_mul64(x, x, &p[0], &p[1]);
        ref_mod(e, p, 2, &n, 1);
        r = sqrredc_pos(x, n, inv);          CHECK(REDC_OK(r), "sqrredc_pos");
        r = sqrredc(x, n, 0 - inv);          CHECK(REDC_OK(r), "sqrredc");
        if (n < ((uint64_t)1 << 63)) {
            r = sqrredc63(x, n, 0 - inv);    CHECK(REDC_OK(r), "sqrredc63");
        }

        /* mfma64 = x*y/2^64 + c (mod n) */
        r = mfma64(x, y, c, n, inv);
        CHECK(r == addmod(mulredc_pos(x, y, n, inv), c, n), "mfma64");

        /* addmod / submod */
        t[0] = x + y; t[1] = (t[0] < x);
        ref_mod(e, t, 2, &n, 1);
        CHECK(addmod(x, y, n) == e[0], "addmod");
        CHECK(submod(addmod(x, y, n), y, n) == x, "submod");
        CHECK(submod(x, y, n) == ((x >= y) ? x - y : x + (n - y)), "submod value");
        if (n <= 0xffffffffULL) {
            CHECK(addmod32((uint32_t)x, (uint32_t)y, (uint32_t)n) == (x + y) % n, "addmod32");
            CHECK(submod32((uint32_t)x, (uint32_t)y, (uint32_t)n) == (x + n - y) % n, "submod32");
        }

        /* gcds */
        a = rndbits(64) | 1; b = rndbits(64);
        g = gcd64(a, b);
        CHECK(bingcd64(a, b) == g, "bingcd64");
        CHECK(spGCD(a, b) == g, "spGCD");
        CHECK(spBinGCD(a, b) == g, "spBinGCD");
        CHECK(spBinGCD_odd(a, b) == g, "spBinGCD_odd");
        CHECK(bingcd64(a, 0) == a, "bingcd64 v=0");

        /* jacobi for a composite modulus is multiplicative over its prime
           factors, and is 0 (not a hang) when n shares a factor with it */
        {
            uint64_t nn = b % 100000;
            CHECK(jacobi_1(nn, 3465) == jacobi_1(nn, 3) * jacobi_1(nn, 3) * jacobi_1(nn, 5)
                                        * jacobi_1(nn, 7) * jacobi_1(nn, 11), "jacobi_1 composite");
        }
        {   /* 128-bit gcds: gcd(f*p, f*q) = f * gcd(p, q); u odd, v with any power of 2 */
            uint64_t gf[2] = { 0, 0 }, gp[2] = { 0, 0 }, gq[2] = { 0, 0 }, gt[4], uu[2], vv[2], w1[2], ex[2];
            gf[0] = rndbits(50) | 1; gp[0] = rndbits(60) | 1; gq[0] = rndbits(64);
            mul128(gf, gp, gt); uu[0] = gt[0]; uu[1] = gt[1];
            mul128(gf, gq, gt); vv[0] = gt[0]; vv[1] = gt[1];
            ref_mul64(gf[0], gcd64(gp[0], gq[0]), &ex[0], &ex[1]);
            gcd128(uu, vv, w1);
            CHECK(w1[0] == ex[0] && w1[1] == ex[1], "gcd128");
#if defined(__SIZEOF_INT128__) && (__SIZEOF_INT128__ == 16)
            bin_gcd128(uu, vv, w1);                     /* clobbers its inputs */
            CHECK(w1[0] == ex[0] && w1[1] == ex[1], "bin_gcd128");
#endif
        }

        /* sp* helpers */
        {
            uint64_t lo, hi, s, cy, w;
            spMultiply(x, y, &lo, &hi);
            ref_mul64(x, y, &p[0], &p[1]);
            CHECK(lo == p[0] && hi == p[1], "spMultiply");
            spAdd(x, y, &s, &cy);
            CHECK(s == x + y && cy == (s < x), "spAdd");
            spSub(x, y, &s, &cy);
            CHECK(s == x - y && cy == (x < y), "spSub");
            spMulMod(x, y, n, &w);
            ref_mod(e, p, 2, &n, 1);
            CHECK(w == e[0], "spMulMod");
            CHECK((int)spBits(n) == 64 - (int)_lead_zcnt64(n), "spBits");
            CHECK(bits64(n) == 64 - (int)_lead_zcnt64(n), "bits64");
            CHECK(u64div(x, n) == (t[0] = 0, t[1] = x, ref_mod(e, t, 2, &n, 1), e[0]), "u64div");
        }
    }

    /* modular exponentiation, PRP, inverse, jacobi, sqrt on primes */
    {
        static const uint64_t primes[] = { 3, 5, 7, 11, 13, 101, 257, 65537, 1000003,
            4294967291ULL, 4294967311ULL, 1000000007ULL, 0xFFFFFFFFFFFFFFC5ULL,
            9223372036854775783ULL, 2305843009213693951ULL, 1152921504606846883ULL };
        int np = (int)(sizeof primes / sizeof primes[0]), i;
        for (i = 0; i < np; i++) {
            uint64_t p = primes[i], u, acc = 1, base, ex;
            for (k = 0; k < 200; k++) {
                uint64_t a = rnd() % p, sq, w;
                int j;
                if (a == 0) a = 1;
                /* spModExp against square-and-multiply on spMulMod */
                ex = rndbits(64); base = a; acc = 1;
                {
                    uint64_t e2 = ex, bb = base;
                    while (e2) { if (e2 & 1) spMulMod(acc, bb, p, &acc); spMulMod(bb, bb, p, &bb); e2 >>= 1; }
                }
                spModExp(a, ex, p, &u);
                CHECK(u == acc, "spModExp");
                spModExp(a, p - 1, p, &u);
                CHECK(u == 1, "spModExp fermat");
                /* jacobi vs Euler's criterion */
                spModExp(a, (p - 1) / 2, p, &u);
                j = jacobi_1(a, p);
                CHECK((u == 1 && j == 1) || (u == p - 1 && j == -1), "jacobi_1");
                if (j == 1) {
                    ShanksTonelli_1(a, p, &sq);
                    spMulMod(sq, sq, p, &w);
                    CHECK(w == a, "ShanksTonelli_1");
                }
                if (p <= 0xffffffffULL) {
                    uint32_t iv = modinv_1((uint32_t)a, (uint32_t)p);
                    CHECK(((uint64_t)iv * a) % p == 1, "modinv_1");
                    iv = modinv_1b((uint32_t)a | 1, 0);     /* inverse mod 2^32 */
                    CHECK((uint32_t)(iv * ((uint32_t)a | 1)) == 1, "modinv_1b");
                }
            }
            if (p > 0xffffffffULL) {
                CHECK(spPRP2(p) == 1, "spPRP2 prime");
            }
        }
        for (k = 0; k < 20000; k++) {
            uint64_t p = (rndbits(31) << 33) | (rnd() >> 31) | 1, u;   /* > 2^32, odd */
            spModExp(2, p - 1, p, &u);
            CHECK(spPRP2(p) == u, "spPRP2");
        }
    }
}

/* reference 128-bit modular multiply: a*b mod n by schoolbook + bitwise mod */
static void ref_mulmod128(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *n)
{
    uint64_t p[4];
    ref_mul(p, a, 2, b, 2);
    ref_mod(r, p, 4, n, 2);
}

static void test_limb2(void)
{
    int k, i;
    for (k = 0; k < 60000; k++) {
        uint64_t a[2], b[2], n[2], z[4], e[4], q[2], r[2], t[4];
        int s = (int)(rnd() % 128);
        a[0] = rnd(); a[1] = (k & 3) ? rndbits(64) : 0;
        b[0] = rnd(); b[1] = (k & 5) ? rndbits(64) : 0;
        if ((k % 9) == 0) { a[0] = a[1] = ~(uint64_t)0; }
        if ((k % 11) == 0) { b[1] |= (uint64_t)1 << 63; b[0] = ~(uint64_t)0; }

        /* plain two-limb arithmetic */
        CHECK(add128(a, b, z) == ref_add(e, a, b, 2) && !ref_cmp(z, e, 2), "add128");
        CHECK(sub128(a, b, z) == ref_sub(e, a, b, 2) && !ref_cmp(z, e, 2), "sub128");
        CHECK(cmp128(a, b) == ref_cmp(a, b, 2), "cmp128");
        shl128(a, s, z);
        e[0] = a[0]; e[1] = a[1];
        for (i = 0; i < s; i++) ref_shl1(e, 2);
        CHECK(!ref_cmp(z, e, 2), "shl128");
        shr128(z, s, t);                                /* (a << s) >> s */
        shl128(t, s, q);
        CHECK(!ref_cmp(q, z, 2), "shr128");
        mul128(a, b, z);   ref_mul(e, a, 2, b, 2);  CHECK(!ref_cmp(z, e, 4), "mul128");
        sqr128(a, z);      ref_mul(e, a, 2, a, 2);  CHECK(!ref_cmp(z, e, 4), "sqr128");
        mul128x64(a, b[0], z); ref_mul(e, a, 2, b, 1); CHECK(!ref_cmp(z, e, 3), "mul128x64");

        /* division: n = q*d + r, r < d */
        if ((b[0] | b[1]) == 0) b[0] = 1;
        CHECK(uint128_div(a, b, q, r) == 0, "uint128_div status");
        ref_mod(e, a, 2, b, 2);
        CHECK(!ref_cmp(r, e, 2), "uint128_div remainder");
        ref_mul(t, q, 2, b, 2);
        ref_add(t, t, r, 2);
        CHECK(!ref_cmp(t, a, 2) && t[2] == 0 && t[3] == 0, "uint128_div quotient");
        mod_128(a, b, z);
        CHECK(!ref_cmp(z, e, 2), "mod_128");

        /* 128-bit Montgomery: n odd, x, y < n */
        n[0] = rnd() | 1; n[1] = (k & 7) ? rndbits(64) : 0;
        if (n[1] == 0 && n[0] < 3) n[0] = 3;
        {
            monty128_t m;
            uint64_t x[2], y[2], w[2], R[3] = { 0, 0, 1 }, one[2], back[2], uno[2] = { 1, 0 };
            ref_mod(x, a, 2, n, 2);
            ref_mod(y, b, 2, n, 2);
            monty128_init(&m, n);
            CHECK(m.rho * n[0] == ~(uint64_t)0, "monty128_init rho");
            ref_mod(one, R, 3, n, 2);
            CHECK(!ref_cmp(m.one, one, 2), "monty128_init one");

            w[0] = a[0]; w[1] = a[1];                   /* unreduced input */
            to_monty128(&m, w);
            t[0] = 0; t[1] = 0; t[2] = x[0]; t[3] = x[1];   /* x * 2^128 */
            ref_mod(e, t, 4, n, 2);
            CHECK(!ref_cmp(w, e, 2), "to_monty128");

            /* r = REDC(x*y)  <=>  r * 2^128 == x*y (mod n) */
            ref_mulmod128(e, x, y, n);
#define REDC128_OK(rr) (chkmod128((rr), n), t[0] = 0, t[1] = 0, t[2] = (rr)[0], t[3] = (rr)[1], \
                        ref_mod(back, t, 4, n, 2), !ref_cmp(back, e, 2))
            mulmod128(x, y, r, &m);          CHECK(REDC128_OK(r), "mulmod128");
            mulmod128n(x, y, r, n, m.rho);   CHECK(REDC128_OK(r), "mulmod128n");
            ref_mulmod128(e, x, x, n);
            sqrmod128(x, r, &m);             CHECK(REDC128_OK(r), "sqrmod128");
            sqrmod128n(x, r, n, m.rho);      CHECK(REDC128_OK(r), "sqrmod128n");

            /* out of Montgomery form again */
            mulmod128(w, uno, back, &m);
            chkmod128(back, n);
            CHECK(!ref_cmp(back, x, 2), "monty128 round trip");

            addmod128(x, y, r, n);
            ref_add(t, x, y, 2); t[2] = (ref_cmp(t, x, 2) < 0);
            ref_mod(e, t, 3, n, 2);
            CHECK(!ref_cmp(r, e, 2), "addmod128");
            {   /* the result may alias either operand, or both */
                uint64_t al[2];
                al[0] = x[0]; al[1] = x[1]; addmod128(al, y, al, n);
                CHECK(!ref_cmp(al, r, 2), "addmod128 w == a");
                al[0] = y[0]; al[1] = y[1]; addmod128(x, al, al, n);
                CHECK(!ref_cmp(al, r, 2), "addmod128 w == b");
                al[0] = y[0]; al[1] = y[1]; submod128(r, al, al, n);
                CHECK(!ref_cmp(al, x, 2), "submod128 w == b");
                al[0] = x[0]; al[1] = x[1]; mulmod128(al, y, al, &m);
                mulmod128(x, y, z, &m);
                CHECK(!ref_cmp(al, z, 2), "mulmod128 w == u");
            }
            submod128(r, y, z, n);
            CHECK(!ref_cmp(z, x, 2), "submod128");
            r[0] = x[0]; r[1] = x[1];
            dblmod128(r, n);
            addmod128(x, x, z, n);
            CHECK(!ref_cmp(r, z, 2), "dblmod128");
        }
    }
    {
        uint64_t z[2] = { 0, 0 }, one[2] = { 1, 0 };
        CHECK(uint128_div(one, z, NULL, NULL) == -1, "uint128_div by zero");
    }
}

#if MP_USE_AVX512F && !defined(LIMB1_NO_VECTOR)
static void test_vec52(void)
{
    const uint64_t M52 = 0x000fffffffffffffULL;
    int k, i;
    for (k = 0; k < 40000; k++) {
        uint64_t n[8], x[8], y[8], inv[8], out[8];
        __m512i vn, vx, vy, vinv, v;
        for (i = 0; i < 8; i++) {
            n[i] = rndbits(52) | 1;
            if (n[i] < 3) n[i] = 3;
            x[i] = rnd() % n[i]; y[i] = rnd() % n[i];
            inv[i] = multiplicative_inverse52(n[i]);
        }
        vn = _mm512_loadu_si512((const void *)n);
        vx = _mm512_loadu_si512((const void *)x);
        vy = _mm512_loadu_si512((const void *)y);
        vinv = _mm512_loadu_si512((const void *)inv);

        v = mulredc52_x8(vx, vy, vn, vinv);
        _mm512_storeu_si512((void *)out, v);
        for (i = 0; i < 8; i++) {
            /* out * 2^52 == x*y (mod n) */
            uint64_t p[2], e[2], t[2], b[2];
            ref_mul64(x[i], y[i], &p[0], &p[1]);
            ref_mod(e, p, 2, &n[i], 1);
            t[0] = out[i] << 52; t[1] = out[i] >> 12;
            ref_mod(b, t, 2, &n[i], 1);
            CHECK(out[i] < n[i] && b[0] == e[0], "mulredc52_x8");
        }
        v = sqrredc52_x8(vx, vn, vinv);
        _mm512_storeu_si512((void *)out, v);
        for (i = 0; i < 8; i++) {
            uint64_t p[2], e[2], t[2], b[2];
            ref_mul64(x[i], x[i], &p[0], &p[1]);
            ref_mod(e, p, 2, &n[i], 1);
            t[0] = out[i] << 52; t[1] = out[i] >> 12;
            ref_mod(b, t, 2, &n[i], 1);
            CHECK(out[i] < n[i] && b[0] == e[0], "sqrredc52_x8");
        }
        v = addmod52_x8(vx, vy, vn);
        _mm512_storeu_si512((void *)out, v);
        for (i = 0; i < 8; i++) CHECK(out[i] == (x[i] + y[i]) % n[i], "addmod52_x8");
        v = submod52_x8(vx, vy, vn);
        _mm512_storeu_si512((void *)out, v);
        for (i = 0; i < 8; i++) CHECK(out[i] == (x[i] + n[i] - y[i]) % n[i], "submod52_x8");

        {   /* product halves */
            __m512i l, h;
            uint64_t lo[8], hi[8];
            mp_mul52lohi_x8(vx, vy, &l, &h);
            _mm512_storeu_si512((void *)lo, l);
            _mm512_storeu_si512((void *)hi, h);
            for (i = 0; i < 8; i++) {
                uint64_t p0, p1;
                ref_mul64(x[i], y[i], &p0, &p1);
                CHECK((lo[i] & M52) == (p0 & M52) && hi[i] == ((p1 << 12) | (p0 >> 52)), "mp_mul52lohi_x8");
            }
            _mm512_storeu_si512((void *)lo, mp_mul52lo_x8(vx, vy));
            _mm512_storeu_si512((void *)hi, mp_mul52hi_x8(vx, vy));
            for (i = 0; i < 8; i++) {
                uint64_t p0, p1;
                ref_mul64(x[i], y[i], &p0, &p1);
                CHECK(lo[i] == (p0 & M52) && hi[i] == ((p1 << 12) | (p0 >> 52)), "mp_mul52lo/hi_x8");
            }
        }
    }
}
#endif

/* ---- micro-benchmarks ---- */
static volatile uint64_t sink;

#define BENCH(label, iters, per, BODY)                                         \
    do { clock_t t0 = clock(); long it_;                                       \
        for (it_ = 0; it_ < (iters); it_++) { BODY; }                          \
        printf("  %-34s %7.2f ns\n", label,                                    \
            1e9 * (double)(clock() - t0) / CLOCKS_PER_SEC / (double)(iters) / (per)); \
    } while (0)

static void bench(void)
{
    uint64_t n = rnd() | 0x8000000000000001ULL, x = rnd() % n, inv = multiplicative_inverse(n);
    uint64_t n62 = n >> 2 | 1, x62 = x % n62, inv62 = multiplicative_inverse(n62);
    uint64_t nn[2], xx[2], yy[2], q[2], r[2];
    monty128_t m;
    const long N = 50000000;

    printf("timings (latency chains unless noted)\n");
    BENCH("mulredc_pos", N, 1, x = mulredc_pos(x, x, n, inv));
    BENCH("mulredc_pos_alt", N, 1, x = mulredc_pos_alt(x, x, n, inv));
    BENCH("mulredc (neg inverse, 64-bit n)", N, 1, x = mulredc(x, x, n, 0 - inv));
    BENCH("mulredc63", N, 1, x62 = mulredc63(x62, x62, n62, 0 - inv62));
    BENCH("addmod + submod", N, 2, x = submod(addmod(x, inv, n), n >> 1, n));
    BENCH("_umul128", N, 1, { uint64_t h; x = _umul128(x | 1, inv, &h) ^ h; });
    BENCH("_udiv128", N / 5, 1, { uint64_t rr; x = _udiv128(x % n, inv, n, &rr) + rr; });
    BENCH("bingcd64", N / 20, 1, x += bingcd64(n, x | 1));
    sink += x + x62;

    nn[0] = rnd() | 1; nn[1] = rnd() | ((uint64_t)1 << 63);
    xx[0] = rnd(); xx[1] = rnd() >> 1; yy[0] = rnd(); yy[1] = rnd() >> 1;
    monty128_init(&m, nn);
    BENCH("mulmod128", N / 2, 1, mulmod128(xx, yy, xx, &m));
    BENCH("sqrmod128", N / 2, 1, sqrmod128(xx, xx, &m));
    BENCH("addmod128", N / 2, 1, addmod128(xx, yy, xx, nn));
    BENCH("uint128_div", N / 10, 1, { xx[0] += 77; uint128_div(nn, xx, q, r); sink += q[0]; });
    BENCH("monty128_init", N / 50, 1, { nn[0] += 2; monty128_init(&m, nn); });
    sink += xx[0];

#if MP_USE_AVX512F && !defined(LIMB1_NO_VECTOR)
    {
        uint64_t a8[8], n8[8], i8[8];
        __m512i vx, vn, vi;
        int i;
        for (i = 0; i < 8; i++) {
            n8[i] = rndbits(52) | ((uint64_t)1 << 51) | 1;
            a8[i] = rnd() % n8[i];
            i8[i] = multiplicative_inverse52(n8[i]);
        }
        vx = _mm512_loadu_si512((const void *)a8);
        vn = _mm512_loadu_si512((const void *)n8);
        vi = _mm512_loadu_si512((const void *)i8);
        BENCH("mulredc52_x8 (per lane)", N, 8, vx = mulredc52_x8(vx, vx, vn, vi));
        _mm512_storeu_si512((void *)a8, vx);
        sink += a8[0];
    }
#endif
}

int main(int argc, char **argv)
{
    printf("limb library self-test\n");
    printf("  compiler path : %s\n",
#if defined(__INTEL_LLVM_COMPILER)
        "Intel icx"
#elif defined(__INTEL_COMPILER)
        "Intel icc"
#elif defined(__clang__) && defined(_MSC_VER)
        "clang-cl"
#elif defined(__clang__)
        "clang"
#elif defined(__GNUC__)
        "gcc"
#elif defined(_MSC_VER)
        "MSVC"
#else
        "unknown"
#endif
    );
    printf("  arch          : %s\n", MP_ARCH_X64 ? "x86-64" : (MP_ARCH_ARM64 ? "arm64" : "other"));
    printf("  gnu asm       : %d   uint128: %d   carry: %s\n", MP_HAS_GNU_ASM_X64,
#ifdef HAS_UINT128
        1,
#else
        0,
#endif
        MP_CARRY_IMPL);
    printf("  BMI2 %d  AVX2 %d  AVX512F %d  IFMA %d\n",
        MP_USE_BMI2, MP_USE_AVX2, MP_USE_AVX512F, MP_USE_IFMA);

    test_primitives();   printf("  primitives   done (%ld checks)\n", nchecks);
    test_limb1();        printf("  limb1        done (%ld checks)\n", nchecks);
    test_limb2();        printf("  limb2        done (%ld checks)\n", nchecks);
#if MP_USE_AVX512F && !defined(LIMB1_NO_VECTOR)
    test_vec52();        printf("  52-bit x8    done (%ld checks)\n", nchecks);
#endif

    if (argc > 1 && strcmp(argv[1], "bench") == 0)
        bench();

    printf("%s: %ld checks, %d failures\n", nfail ? "FAILED" : "PASSED", nchecks, nfail);
    return nfail ? 1 : 0;
}
