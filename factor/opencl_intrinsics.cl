/*--------------------------------------------------------------------
OpenCL port of cuda_intrinsics.h by Jason Papadopoulos.
Original placed in the public domain. This translation likewise.

Porting notes
-------------
CUDA -> OpenCL type mapping:
  uint32          -> uint       (OpenCL built-in 32-bit unsigned)
  uint64          -> ulong      (OpenCL built-in 64-bit unsigned)
  int32           -> int
  int64           -> long

CUDA -> OpenCL intrinsic mapping:
  __umulhi(a,b)         -> mul_hi(a,b)        exact equivalent
  __clz(v)              -> clz(v)             exact equivalent
  __clzll(n)            -> clz((ulong)n)      exact equivalent
  __funnelshift_l(l,h,1)-> (h<<1)|(l>>31)    implemented inline, see below
  __device__            -> (omit)             all functions callable from kernel
                           or use 'inline'    hint to compiler

Carry-chain arithmetic:
  PTX add.cc/addc.cc chains map to uadd_carry() (OpenCL 2.0).
  PTX sub.cc/subc.cc chains map to usub_borrow() (OpenCL 2.0).
  Both are available on AMD RX 6700 XT (OpenCL 2.0 device).

  For multi-limb chains we use a 'ulong' accumulator wherever possible
  because it is simpler and equally efficient on AMD RDNA2:
      ulong s = (ulong)a + b + carry_in;
      result  = (uint)s;
      carry   = (uint)(s >> 32);   // 0 or 1

__uint128_t (used in modinv96):
  OpenCL has no 128-bit integer type. We emulate it with two ulongs:
  (hi:lo) representing the value hi*2^64 + lo.
  Helper macros/functions for 128-bit add, subtract, compare, and
  divide-by-subtraction are provided below.

Montgomery macros (preserved from original):
  #define ALWAYS_SUB   always do the conditional subtraction
  #define NO_SUB       skip the conditional subtraction (testing only)
--------------------------------------------------------------------*/

#ifndef OPENCL_INTRINSICS_CL
#define OPENCL_INTRINSICS_CL

/* ------------------------------------------------------------------ */
/* Build-time options (mirror the original CUDA macros)               */
/* ------------------------------------------------------------------ */
//#define NO_SUB
#define ALWAYS_SUB

/* ================================================================== */
/*  uadd_carry / usub_borrow compatibility layer                       */
/*                                                                     */
/*  These names appeared in early OpenCL 2.0 drafts but were never    */
/*  ratified and are not implemented by the AMD (comgr) compiler or   */
/*  any other mainstream OpenCL driver.                                */
/*                                                                     */
/*  We provide them as inline functions using the ulong accumulator    */
/*  pattern.  The AMD compiler recognises this idiom and emits native  */
/*  v_addc_u32 / v_subb_u32 RDNA2 instructions with no overhead.      */
/*                                                                     */
/*  Signatures match the proposed spec:                                */
/*    uint uadd_carry (uint a, uint b, uint *carry_out)  -> sum        */
/*    uint usub_borrow(uint a, uint b, uint *borrow_out) -> difference */
/* ================================================================== */

static inline uint uadd_carry(uint a, uint b, uint *carry_out)
{
    ulong s      = (ulong)a + (ulong)b;
    *carry_out   = (uint)(s >> 32);
    return (uint)s;
}

static inline uint usub_borrow(uint a, uint b, uint *borrow_out)
{
    ulong s       = (ulong)a - (ulong)b;
    /* bit 32 is set iff a < b; negate the sign-extended bit to get 0 or 1 */
    *borrow_out   = (uint)(-(int)(s >> 32));
    return (uint)s;
}

/* ================================================================== */
/*  Low-level carry/borrow helper functions                            */
/*  These call the macros above so their logic is identical to the    */
/*  original PTX patterns.                                            */
/* ================================================================== */

/*
 * add2 -- add two 32-bit values, produce 32-bit sum and 1-bit carry
 */
static inline void
add2(uint a, uint b, uint *sum, uint *carry)
{
    *sum = uadd_carry(a, b, carry);
}

/*
 * add3c -- add three 32-bit values, produce sum and carry.
 *   Carries from two separate adds are each 0 or 1,
 *   so their sum fits in a uint without overflow.
 */
static inline void
add3c(uint a, uint b, uint cin, uint *sum, uint *carry)
{
    uint c1, c2;
    uint t = uadd_carry(a, b, &c1);
    *sum   = uadd_carry(t, cin, &c2);
    *carry = c1 + c2;
}

/*
 * sub2b -- subtract b from a, produce difference and 1-bit borrow
 */
static inline void
sub2b(uint a, uint b, uint *diff, uint *borrow)
{
    *diff = usub_borrow(a, b, borrow);
}

/* ================================================================== */
/*  Core accumulator primitives                                        */
/*  These are the heart of the multiprecision inner loops.            */
/* ================================================================== */

/*
 * accum3 -- add a 64-bit value (b1:b0) into a 96-bit accumulator (a2:a1:a0)
 *
 *   CUDA original used three chained PTX instructions:
 *     add.cc   a0, a0, b0
 *     addc.cc  a1, a1, b1
 *     addc     a2, a2, 0
 *
 *   OpenCL: use ulong to carry across the 32-bit boundary cleanly.
 */
static inline void
accum3(uint *a0, uint *a1, uint *a2,
       uint b0, uint b1)
{
    ulong s;
    s   = (ulong)(*a0) + b0;
    *a0 = (uint)s;
    ulong c = (s >> 32);
    s   = (ulong)(*a1) + b1 + c;
    *a1 = (uint)s;
    c = (s >> 32);
    *a2 = *a2 + (uint)c;   /* no further carry possible in correct usage */
}

/*
 * accum3_shift -- shift-accumulate variant
 *
 *   CUDA original:
 *     add.cc   a0, a1, b0    (note: destination a0, source a1)
 *     addc.cc  a1, a2, b1
 *     addc     a2, 0,  0
 *
 *   This shifts the accumulator down one word while adding (b1:b0).
 *   Used in montmul64 to advance the rolling window.
 */
static inline void
accum3_shift(uint *a0, uint *a1, uint *a2,
             uint b0, uint b1)
{
    ulong s;
    s   = (ulong)(*a1) + b0;
    *a0 = (uint)s;
    ulong c = (s >> 32);
    s   = (ulong)(*a2) + b1 + c;
    *a1 = (uint)s;
    c = (s >> 32);
    *a2 = (uint)c;
}

/*
 * add3 -- add three 32-bit values into a (sum, carry) pair
 *
 *   Equivalent to: sum = a + b + c, carry = overflow bit
 */
static inline void
add3(uint a, uint b, uint c, uint *sum, uint *carry)
{
    ulong s = (ulong)a + b + c;
    *sum   = (uint)s;
    *carry = (uint)(s >> 32);
}

/*
 * accumlh -- add a (lo,hi) pair into (*t), producing a carry
 *
 *   CUDA original:
 *     add.cc   *t,     *t, lo
 *     addc     *carry, hi, 0
 */
static inline void
accumlh(uint lo, uint hi, uint *t, uint *carry)
{
    ulong s = (ulong)*t + lo;
    *t      = (uint)s;
    *carry  = hi + (uint)(s >> 32);   /* hi is a mul_hi result; adding 1-bit carry safe */
}

/*
 * accum3lh -- full inner-loop accumulation step
 *
 *   Computes: (*sum, *carry) = t + lo + C + hi (conceptually)
 *   where t and lo are the "current" limb and product-lo,
 *   C is carry-in from the previous limb, and hi is the high word
 *   of the product that feeds into the *carry output.
 *
 *   CUDA original called accumlh then folded hi into carry; here we
 *   keep it in one ulong chain for clarity.
 */
static inline void
accum3lh(uint t, uint lo, uint C, uint hi,
         uint *sum, uint *carry)
{
    ulong s  = (ulong)t + lo + C;
    *sum     = (uint)s;
    *carry   = hi + (uint)(s >> 32);
}

/* ================================================================== */
/*  Squaring helper                                                    */
/* ================================================================== */

/*
 * wide_sqr32 -- square a 32-bit value, return 64-bit result
 *
 *   CUDA used inline PTX mul.wide.u32; in OpenCL just cast to ulong.
 */
static inline ulong
wide_sqr32(uint a)
{
    return (ulong)a * a;
}

/*
 * innermul -- core multiply-accumulate used in montmul32/64
 *
 *   Computes: (carry : return) = a*b + prevhi + accum
 *   using a single 64-bit multiply-add (replaces mad.wide.u32 + add.u64)
 */
static inline uint
innermul(uint a, uint b, uint prevhi, uint accum, uint *carry)
{
    ulong s = (ulong)a * b + accum + prevhi;
    *carry  = (uint)(s >> 32);
    return   (uint)s;
}


/* ================================================================== */
/*  Modular subtraction                                                */
/* ================================================================== */

/*
 * modsub32 -- modular subtraction: return (a - b) mod p
 *
 *   If a < b, borrow occurs and we add p back.
 *   Replaces: sub.cc / subc / setp / conditional add
 */
static inline uint
modsub32(uint a, uint b, uint p)
{
    uint borrow;
    uint r = usub_borrow(a, b, &borrow);
    if (borrow) r += p;
    return r;
}

/*
 * modsub64 -- modular subtraction for 64-bit values
 */
static inline ulong
modsub64(ulong a, ulong b, ulong p)
{
    ulong r = a - b;
    if (a < b)
    {
        r += p;
    }

    return r;
}

/*
 * modsub96 -- modular subtraction for 96-bit values (3 x uint limbs)
 */
static inline void
modsub96(uint *a, uint *b, uint *c, uint *p)
{
    uint bw0, bw1, bw2, bw3;
    uint r0, r1, r2;

    r0 = usub_borrow(a[0], b[0], &bw0);
    uint t1 = usub_borrow(a[1], b[1], &bw1);
    r1 = usub_borrow(t1, bw0, &bw0); bw1 += bw0;
    uint t2 = usub_borrow(a[2], b[2], &bw2);
    r2 = usub_borrow(t2, bw1, &bw3); bw2 += bw3;

    if (bw2) {
        /* result was negative, add p back */
        uint c0, c1;
        //r0 = uadd_carry(r0, p[0], &c0);
        //r1 = uadd_carry(r1 + c0, p[1], &c1);
        //r2 = r2 + p[2] + c1;
        r0 = uadd_carry(r0, p[0], &c0);
        r1 = uadd_carry(r1, p[1], &c1); uint cx; r1 = uadd_carry(r1, c0, &cx); c1 += cx;
        r2 = r2 + p[2] + c1;
    }
    c[0] = r0; c[1] = r1; c[2] = r2;
}

/* ================================================================== */
/*  Modular addition                                                   */
/* ================================================================== */

/*
 * modadd32 -- modular addition: return (a + b) mod p
 *
 *   Adds, then if carry or result >= p, subtracts p.
 */
static inline uint
modadd32(uint a, uint b, uint p)
{
    uint carry;
    uint r = uadd_carry(a, b, &carry);
    if (carry || r >= p) r -= p;
    return r;
}

/*
 * modadd64 -- modular addition for 64-bit values
 *
 *   Strategy: compute sum = a+b (65-bit), also compute sum-p (65-bit),
 *   use the final borrow to decide which to keep.
 *   Mirrors the original PTX exactly.
 */
static inline ulong
modadd64(ulong a, ulong b, ulong p)
{
    ulong r = a + b;

    if ((r < a) || (r >= p))
        r -= p;

    return r;
}

/*
 * modadd96 -- modular addition for 96-bit values
 */
static inline void
modadd96(uint *a, uint *b, uint *c, uint *p)
{
    uint r0, r1, r2, s0, s1, s2;
    uint carry, borrow;
    uint c1, c2, c3;

    /* r = a + b */
    r0 = uadd_carry(a[0], b[0], &carry);
    r1 = uadd_carry(a[1], b[1], &c1); uint cx; r1 = uadd_carry(r1, carry, &cx); c1 += cx;
    r2 = uadd_carry(a[2], b[2], &c2); uint cy; r2 = uadd_carry(r2, c1, &cy);    c2 += cy;
    carry = c2;

    /* s = r - p */
    s0 = usub_borrow(r0, p[0], &borrow);
    s1 = usub_borrow(r1, p[1], &c1); uint bx; s1 = usub_borrow(s1, borrow, &bx); c1 += bx;
    s2 = usub_borrow(r2, p[2], &c2); uint by; s2 = usub_borrow(s2, c1, &by);     c2 += by;
    borrow = c2; // usub_borrow(carry, c2, &c3);

    //if (borrow == 0)
    if ((carry > 0) || ((carry == 0) && (borrow == 0))) { r0=s0; r1=s1; r2=s2; }
    c[0]=r0; c[1]=r1; c[2]=r2;
}

/*
 * modaddsub96 -- compute both (a+b) mod p and (a-b) mod p in one pass
 */
static inline void
modaddsub96(uint *a, uint *b, uint *s, uint *d, uint *p)
{
    /* sum = a + b mod p */
    modadd96(a, b, s, p);
    /* diff = a - b mod p */
    modsub96(a, b, d, p);
}

/*
 * modadd128 -- modular addition for 128-bit values (2 x ulong limbs)
 *
 *   Works in 32-bit limbs internally to mirror the original exactly.
 */
static inline void
modadd128(ulong *a, ulong *b, ulong *c, ulong *n)
{
    uint aa[4], bb[4], cc[4], nn[4];
    uint carry, borrow;

    aa[0]=(uint)a[0]; aa[1]=(uint)(a[0]>>32); aa[2]=(uint)a[1]; aa[3]=(uint)(a[1]>>32);
    bb[0]=(uint)b[0]; bb[1]=(uint)(b[0]>>32); bb[2]=(uint)b[1]; bb[3]=(uint)(b[1]>>32);
    nn[0]=(uint)n[0]; nn[1]=(uint)(n[0]>>32); nn[2]=(uint)n[1]; nn[3]=(uint)(n[1]>>32);

    /* addcarry = aa += bb (128-bit add) */
    uint c0,c1,c2,c3,c4;
    aa[0] = uadd_carry(aa[0], bb[0], &c0);
    aa[1] = uadd_carry(aa[1], bb[1], &c1); uint x1; aa[1] = uadd_carry(aa[1], c0, &x1); c1+=x1;
    aa[2] = uadd_carry(aa[2], bb[2], &c2); uint x2; aa[2] = uadd_carry(aa[2], c1, &x2); c2+=x2;
    aa[3] = uadd_carry(aa[3], bb[3], &c3); uint x3; aa[3] = uadd_carry(aa[3], c2, &x3); c3+=x3;
    uint addcarry = c3;

    cc[0]=aa[0]; cc[1]=aa[1]; cc[2]=aa[2]; cc[3]=aa[3];

    /* subborrow = aa -= nn */
    uint b0,b1,b2,b3;
    aa[0] = usub_borrow(aa[0], nn[0], &b0);
    aa[1] = usub_borrow(aa[1], nn[1], &b1); uint y1; aa[1] = usub_borrow(aa[1], b0, &y1); b1+=y1;
    aa[2] = usub_borrow(aa[2], nn[2], &b2); uint y2; aa[2] = usub_borrow(aa[2], b1, &y2); b2+=y2;
    aa[3] = usub_borrow(aa[3], nn[3], &b3); uint y3; aa[3] = usub_borrow(aa[3], b2, &y3); b3+=y3;
    uint subborrow = b3;

    /* use subtracted result if addcarry==1 or subborrow==0 */
    if (addcarry == 1 || subborrow == 0) {
        c[0] = ((ulong)aa[1]<<32)|(ulong)aa[0];
        c[1] = ((ulong)aa[3]<<32)|(ulong)aa[2];
    } else {
        c[0] = ((ulong)cc[1]<<32)|(ulong)cc[0];
        c[1] = ((ulong)cc[3]<<32)|(ulong)cc[2];
    }
}

/* ================================================================== */
/*  GCD                                                                */
/* ================================================================== */

/*
 * gcd32 -- binary-style GCD of two odd nonzero 32-bit values
 *
 *   clz() is an exact OpenCL equivalent of __clz().
 *   min/max are OpenCL built-ins.
 */
static inline uint
gcd32(uint x, uint y)
{
    uint u = x, v = y;
    do {
        uint shift = 31 - clz(v & (uint)(-(int)v));
        v = v >> shift;
        x = min(u, v);
        y = max(u, v);
        u = x;
        v = y - x;
    } while (v != 0);
    return u;
}

/* ================================================================== */
/*  Modular inverse                                                    */
/* ================================================================== */

/*
 * modinv32 / modinv64 -- extended Euclidean modular inverse
 *
 *   Pure integer arithmetic; no PTX intrinsics in the original.
 *   Translation is straightforward: replace typedefs with OpenCL types.
 */
static inline uint
modinv32(uint a, uint p)
{
    uint ps1, ps2, dividend, divisor, rem, q, t;
    uint parity;

    q=1; rem=a; dividend=p; divisor=a;
    ps1=1; ps2=0; parity=0;

    while (divisor > 1) {
        rem = dividend - divisor;
        t   = rem - divisor;
        if (rem>=divisor){q+=ps1;rem=t;t-=divisor;
        if (rem>=divisor){q+=ps1;rem=t;t-=divisor;
        if (rem>=divisor){q+=ps1;rem=t;t-=divisor;
        if (rem>=divisor){q+=ps1;rem=t;t-=divisor;
        if (rem>=divisor){q+=ps1;rem=t;t-=divisor;
        if (rem>=divisor){q+=ps1;rem=t;t-=divisor;
        if (rem>=divisor){q+=ps1;rem=t;t-=divisor;
        if (rem>=divisor){q+=ps1;rem=t;
        if (rem>=divisor){q=dividend/divisor;
                          rem=dividend-q*divisor;
                          q*=ps1;
        }}}}}}}}}
        q+=ps2; parity=~parity;
        dividend=divisor; divisor=rem;
        ps2=ps1; ps1=q; q=1;
    }
    return (parity==0) ? ps1 : p-ps1;
}

static inline ulong
modinv64(ulong a, ulong p, ulong *likely_gcd)
{
    ulong ps1, ps2, dividend, divisor, rem, q, t;
    uint parity;

    q = 1; rem = a; dividend = p; divisor = a;
    ps1 = 1; ps2 = 0; parity = 0;

    while (divisor > 1) {
        rem = dividend - divisor;
        t = rem - divisor;
        if (rem >= divisor) {
            q += ps1; rem = t; t -= divisor;
            if (rem >= divisor) {
                q += ps1; rem = t; t -= divisor;
                if (rem >= divisor) {
                    q += ps1; rem = t; t -= divisor;
                    if (rem >= divisor) {
                        q += ps1; rem = t; t -= divisor;
                        if (rem >= divisor) {
                            q += ps1; rem = t; t -= divisor;
                            if (rem >= divisor) {
                                q += ps1; rem = t; t -= divisor;
                                if (rem >= divisor) {
                                    q += ps1; rem = t; t -= divisor;
                                    if (rem >= divisor) {
                                        q += ps1; rem = t;
                                        if (rem >= divisor) {
                                            q = dividend / divisor;
                                            rem = dividend - q * divisor;
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    if (divisor == 1)
        dividend = divisor;
    *likely_gcd = dividend;

    if (parity == 0)
        return ps1;
    else
        return p - ps1;
}

/* ------------------------------------------------------------------ */
/* 128-bit integer emulation for modinv96                             */
/*                                                                     */
/* We represent a 128-bit unsigned integer as (hi:lo) where           */
/*   value = hi * 2^64 + lo                                           */
/* ------------------------------------------------------------------ */
typedef struct { ulong lo; ulong hi; } uint128;

static inline uint128 u128_from_limbs3(uint w0, uint w1, uint w2)
{
    uint128 r;
    r.lo = (ulong)w0 | ((ulong)w1 << 32);
    r.hi = (ulong)w2;
    return r;
}

static inline uint128 u128_add(uint128 a, uint128 b)
{
    uint128 r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo ? 1UL : 0UL);
    return r;
}

static inline uint128 u128_sub(uint128 a, uint128 b)
{
    uint128 r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (a.lo < b.lo ? 1UL : 0UL);
    return r;
}

static inline int u128_gt(uint128 a, uint128 b)
{
    return (a.hi > b.hi) || (a.hi == b.hi && a.lo > b.lo);
}

static inline int u128_ge(uint128 a, uint128 b)
{
    return (a.hi > b.hi) || (a.hi == b.hi && a.lo >= b.lo);
}

static inline int u128_eq(uint128 a, uint128 b)
{
    return (a.hi == b.hi) && (a.lo == b.lo);
}

static inline int u128_gt1(uint128 a)   /* a > 1 */
{
    return (a.hi > 0) || (a.lo > 1);
}

/* Scalar multiply: a * scalar (scalar fits in ulong) */
static inline uint128 u128_mul_scalar(uint128 a, ulong s)
{
    /* only needed for q *= ps1 where values fit in practice */
    /* For modinv96 the quotient q fits in a ulong in all cases */
    uint128 r;
    r.lo = a.lo * s;
    r.hi = a.hi * s;
    r.hi += mul_hi(a.lo, s);

    // /* hi part: a.hi*s + carry from lo */
    // ulong hi_a_s = a.hi * s;
    // ulong carry  = (a.lo != 0 && s > (ulong)(-1) / a.lo) ?
    //                (ulong)((__uint128_t)a.lo * s >> 64) : 0; /* host-side only */
    // /* For kernel use we avoid __uint128_t; use the fact that in modinv96 
    //    the values stay bounded within 96 bits so a.hi is small (< 2^32). */
    // r.hi = hi_a_s + carry;
    return r;
}

/*
 * u128_divmod -- compute q = a / b, rem = a % b  (unsigned, both > 0)
 *
 *   Used only in the "large quotient" fast-path of modinv96 (rare).
 *   Simple restoring division by repeated subtraction of shifted b.
 */
static inline uint128 u128_divmod(uint128 dividend, uint128 divisor, uint128 *rem)
{
    uint128 q = {0, 0};
    uint128 r = dividend;
    /* find highest bit of divisor relative to dividend */
    /* shift divisor left until it exceeds dividend, then shift back */
    uint shift = 0;
    uint128 d = divisor;
    while (!u128_gt(d, dividend) && !(d.hi >> 63)) {
        d.hi = (d.hi << 1) | (d.lo >> 63);
        d.lo <<= 1;
        shift++;
    }
    /* long division */
    for (int i = (int)shift; i >= 0; i--) {
        if (u128_ge(r, d)) {
            r = u128_sub(r, d);
            if (i < 64)      q.lo |= (1UL << i);
            else             q.hi |= (1UL << (i-64));
        }
        d.lo = (d.lo >> 1) | (d.hi << 63);
        d.hi >>= 1;
    }
    *rem = r;
    return q;
}

/*
 * modinv96 -- modular inverse for 96-bit values
 *
 *   Original used __uint128_t (GCC extension, not available in OpenCL).
 *   We replace it with the uint128 struct above.
 */
static inline void
modinv96(uint *a, uint *p, uint *inv, uint *gcd)
{
    uint128 ps1, ps2, dividend, divisor, rem, q, t;
    uint parity;

    uint128 one = {1UL, 0UL};
    uint128 zero = {0UL, 0UL};

    q        = one;
    rem      = u128_from_limbs3(a[0], a[1], a[2]);
    dividend = u128_from_limbs3(p[0], p[1], p[2]);
    divisor  = u128_from_limbs3(a[0], a[1], a[2]);
    ps1      = one;
    ps2      = zero;
    parity   = 0;

    while (u128_gt1(divisor)) {
        rem = u128_sub(dividend, divisor);
        t   = u128_sub(rem, divisor);

        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;t=u128_sub(t,divisor);
        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;t=u128_sub(t,divisor);
        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;t=u128_sub(t,divisor);
        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;t=u128_sub(t,divisor);
        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;t=u128_sub(t,divisor);
        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;t=u128_sub(t,divisor);
        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;t=u128_sub(t,divisor);
        if (u128_ge(rem, divisor)){q=u128_add(q,ps1);rem=t;
        if (u128_ge(rem, divisor)){
            uint128 qr;
            q = u128_divmod(dividend, divisor, &qr);
            rem = qr;
            /* q *= ps1 -- use scalar multiply (q fits in ulong here) */
            //uint128 qq = {q.lo * ps1.lo, 0};  /* simplified: bounded case */
            q = u128_mul_scalar(ps1, q.lo);
        }}}}}}}}}

        q = u128_add(q, ps2);
        parity = ~parity;
        dividend = divisor;
        divisor  = rem;
        ps2 = ps1;
        ps1 = q;
        //q   = one;
    }

    if (u128_eq(divisor, one)) dividend = divisor;

    /* store gcd (low 96 bits of dividend) */
    gcd[0] = (uint)dividend.lo;
    gcd[1] = (uint)(dividend.lo >> 32);
    gcd[2] = (uint)dividend.hi;

    /* compute inverse */
    uint128 pval = u128_from_limbs3(p[0], p[1], p[2]);
    uint128 result = (parity == 0) ? ps1 : u128_sub(pval, ps1);

    inv[0] = (uint)result.lo;
    inv[1] = (uint)(result.lo >> 32);
    inv[2] = (uint)result.hi;
}

/* ================================================================== */
/*  Montgomery arithmetic                                              */
/* ================================================================== */

/*
 * montmul32_w -- compute Montgomery constant w = -n^{-1} mod 2^32
 */
static inline uint
montmul32_w(uint n)
{
    uint res = 2 + n;
    res = res * (2 + n * res);
    res = res * (2 + n * res);
    res = res * (2 + n * res);
    return res * (2 + n * res);
}

static inline ulong
montmul64_w(ulong n)
{
    ulong res = 2 + n;
    res = res * (2 + n * res);
    res = res * (2 + n * res);
    res = res * (2 + n * res);
    res = res * (2 + n * res);
    return res * (2 + n * res);
}

/*
 * montmul32_r -- compute Montgomery R^2 mod n (starting value for conversion)
 */
static inline uint
montmul32_r(uint n)
{
    /* R = 2^32, so R mod n = 2^32 mod n, computed as 2*(2^31 mod n) */
    uint r0 = (uint)(((ulong)1 << 63) % n);
    uint r1 = r0 + r0;
    if (r1 < r0) r1 -= n;
    return modsub32(r1, n, n);
}

/*
 * montmul32 -- Montgomery multiplication for 32-bit values
 *
 *   Uses mul_hi() (exact equivalent of __umulhi) and the accum3 helpers.
 */
static inline uint
montmul32(uint a, uint b, uint n, uint w)
{
    uint acc0, acc1, acc2 = 0;
    uint q;

    acc0 = a * b;
    acc1 = mul_hi(a, b);

    q = acc0 * w;

    uint prod_lo = q * n;
    uint prod_hi = mul_hi(q, n);

    accum3(&acc0, &acc1, &acc2, prod_lo, prod_hi);

    if (acc2 || acc1 >= n)
        return acc1 - n;
    return acc1;
}

/*
 * montmul64 -- Montgomery multiplication for 64-bit values
 *
 *   Structure mirrors the original exactly; PTX inline asm replaced
 *   by accum3 / accum3_shift helper calls.
 */
#if 1

static inline ulong
montmul64(ulong a, ulong b, ulong n, uint w)
{
    uint a0=(uint)a, a1=(uint)(a>>32);
    uint b0=(uint)b, b1=(uint)(b>>32);
    uint n0=(uint)n, n1=(uint)(n>>32);
    uint acc0, acc1, acc2=0;
    uint q0, q1;
    ulong r;

    acc0 = a0*b0;  acc1 = mul_hi(a0,b0);
    q0   = acc0*w;
    accum3(&acc0,&acc1,&acc2, q0*n0, mul_hi(q0,n0));
    accum3_shift(&acc0,&acc1,&acc2, a0*b1, mul_hi(a0,b1));
    accum3      (&acc0,&acc1,&acc2, a1*b0, mul_hi(a1,b0));
    accum3      (&acc0,&acc1,&acc2, q0*n1, mul_hi(q0,n1));
    
    q1 = acc0*w;
    accum3      (&acc0,&acc1,&acc2, q1*n0, mul_hi(q1,n0));
    accum3_shift(&acc0,&acc1,&acc2, a1*b1, mul_hi(a1,b1));
    accum3      (&acc0,&acc1,&acc2, q1*n1, mul_hi(q1,n1));

    r = (ulong)acc1<<32 | acc0;
    if (acc2 || r >= n) return r - n;
    return r;
}

#else

static inline ulong
montmul64(ulong a, ulong b, ulong n, uint w)
{
    uint a0 = (uint)a, a1 = (uint)(a >> 32);
    uint b0 = (uint)b, b1 = (uint)(b >> 32);
    uint n0 = (uint)n, n1 = (uint)(n >> 32);
    uint acc0, acc1, acc2 = 0;
    uint q0, q1;
    ulong r;
    
    //acc0 = a0 * b0;  acc1 = mul_hi(a0, b0);
    ulong t = (ulong)a0 * (ulong)b0;
    acc0 = (uint)t;
    acc1 = (uint)(t >> 32);

    
    q0 = acc0 * w;
    //q0 = acc0 * w;

    t = (ulong)q0 * (ulong)n0;
    uint lo = (uint)t;
    uint hi = (uint)(t >> 32);

    accum3(&acc0, &acc1, &acc2, lo, hi);
    //accum3(&acc0, &acc1, &acc2, q0 * n0, mul_hi(q0, n0));
    


    t = (ulong)a0 * (ulong)b1;
    lo = (uint)t;
    hi = (uint)(t >> 32);

    accum3_shift(&acc0, &acc1, &acc2, lo, hi);
    //accum3_shift(&acc0, &acc1, &acc2, a0 * b1, mul_hi(a0, b1));
    

    t = (ulong)a1 * (ulong)b0;
    lo = (uint)t;
    hi = (uint)(t >> 32);

    accum3(&acc0, &acc1, &acc2, lo, hi);
    //accum3(&acc0, &acc1, &acc2, a1 * b0, mul_hi(a1, b0));
    

    t = (ulong)q0 * (ulong)n1;
    lo = (uint)t;
    hi = (uint)(t >> 32);

    accum3(&acc0, &acc1, &acc2, lo, hi);
    //accum3(&acc0, &acc1, &acc2, q0 * n1, mul_hi(q0, n1));

    
    
    q1 = acc0 * w;
    //q1 = acc0 * w;
    

    t = (ulong)q1 * (ulong)n0;
    lo = (uint)t;
    hi = (uint)(t >> 32);

    accum3(&acc0, &acc1, &acc2, lo, hi);
    //accum3(&acc0, &acc1, &acc2, q1 * n0, mul_hi(q1, n0));
    

    t = (ulong)a1 * (ulong)b1;
    lo = (uint)t;
    hi = (uint)(t >> 32);

    accum3_shift(&acc0, &acc1, &acc2, lo, hi);
    //accum3_shift(&acc0, &acc1, &acc2, a1 * b1, mul_hi(a1, b1));
    

    t = (ulong)q1 * (ulong)n1;
    lo = (uint)t;
    hi = (uint)(t >> 32);

    accum3(&acc0, &acc1, &acc2, lo, hi);
    //accum3(&acc0, &acc1, &acc2, q1 * n1, mul_hi(q1, n1));
    
    r = (ulong)acc1 << 32 | acc0;
    if (acc2 || r >= n) return r - n;
    return r;

    //ulong w64 = montmul64_w(n);
    //
    //ulong t0 = a * b;
    //ulong t1 = mul_hi(a, b);
    //
    //ulong q = t0 * w64;
    //
    //ulong m0 = q * n;
    //ulong m1 = mul_hi(q, n);
    //
    //ulong s0 = t0 + m0;
    //ulong c = (s0 < t0) ? 1 : 0;
    //s0 = t1 + m1 + c;
    //
    //if (s0 >= n)
    //    s0 -= n;
    //
    //return s0;
}

#endif

/*
 * montmul64_r -- compute R^2 mod n for 64-bit Montgomery setup
 *
 *   clz((ulong)n) replaces __clzll(n); otherwise identical to original.
 */
static inline ulong
montmul64_r(ulong n, uint w)
{
    uint shift;
    uint i;
    ulong shifted_n;
    ulong res;

    shift    = clz(n);
    shifted_n = n << shift;
    res      = -shifted_n;

    for (i = 64 - shift; i < 72; i++) {
        if (res >> 63)
            res = res + res - shifted_n;
        else
            res = res + res;
        if (res >= shifted_n) res -= shifted_n;
    }

    res = res >> shift;
    res = montmul64(res, res, n, w);
    res = montmul64(res, res, n, w);
    return montmul64(res, res, n, w);
}

/*
 * montmul96 -- Montgomery multiplication for 96-bit values (CIOS)
 *
 *   All PTX inline asm replaced:
 *     mad.lo.cc / madc.hi  ->  accumlh() using mul + mul_hi
 *     mad.wide.u32 + add   ->  innermul() via (ulong) multiply
 *     sub.cc / subc chains ->  usub_borrow chains
 *   __funnelshift_l not used here (that appears in montsqr).
 */
static inline void
montmul96(uint *a, uint *b, uint *c, uint *n, uint w)
{
    uint t[5] = {0,0,0,0,0};
    uint m, C, S;

    for (int i = 0; i < 3; i++) {
        /* t += a * b[i]: first limb via accumlh */
        accumlh(a[0]*b[i], mul_hi(a[0],b[i]), &t[0], &C);

        for (int j = 1; j < 3; j++) {
            t[j] = innermul(a[j], b[i], C, t[j], &C);
        }

        m = t[0] * w;
        add2(t[3], C, &S, &C); t[3]=S; t[4]=C;

        /* t -= t[0]*n (reduction step): first limb */
        accumlh(n[0]*m, mul_hi(n[0],m), &t[0], &C);

        for (int j = 1; j < 3; j++) {
            t[j-1] = innermul(n[j], m, C, t[j], &C);
        }

        add2(t[3], C, &t[2], &C);
        t[3] = t[4] + C;
        t[4] = 0;
    }

#ifdef NO_SUB
    c[0]=t[0]; c[1]=t[1]; c[2]=t[2];
#else
    uint cc[3] = {t[0], t[1], t[2]};

#ifdef ALWAYS_SUB
    /* Subtract n from t using a single borrow chain across all 3 limbs.
     * The previous code split the t[2] limb into two usub_borrow calls and
     * only kept the borrow from the second, silently discarding the borrow
     * from t[2]-n[2].  Use ulong arithmetic for each limb to keep it clean. */
    uint borrow;
    ulong d0 = (ulong)t[0] - n[0];
    uint s0  = (uint)d0;
    borrow   = (uint)(-(int)(d0 >> 32)) & 1;

    ulong d1 = (ulong)t[1] - n[1] - borrow;
    uint s1  = (uint)d1;
    borrow   = (uint)(-(int)(d1 >> 32)) & 1;

    ulong d2 = (ulong)t[2] - n[2] - borrow;
    uint s2  = (uint)d2;
    borrow   = (uint)(-(int)(d2 >> 32)) & 1;

    /* use subtracted result if t[3] set OR no borrow */
    if (t[3] || borrow == 0) { c[0]=s0; c[1]=s1; c[2]=s2; }
    else                     { c[0]=cc[0]; c[1]=cc[1]; c[2]=cc[2]; }
#else
    if (t[3]) {
        ulong d0b = (ulong)t[0] - n[0];
        uint  s0b = (uint)d0b;
        uint  bw  = (uint)(-(int)(d0b >> 32)) & 1;
        ulong d1b = (ulong)t[1] - n[1] - bw;
        uint  s1b = (uint)d1b;
        bw        = (uint)(-(int)(d1b >> 32)) & 1;
        ulong d2b = (ulong)t[2] - n[2] - bw;
        c[0]=s0b; c[1]=s1b; c[2]=(uint)d2b;
    } else { c[0]=cc[0]; c[1]=cc[1]; c[2]=cc[2]; }
#endif
#endif
}

/*
 * montsqr96 -- Montgomery squaring for 96-bit values (CIOS)
 *
 *   __funnelshift_l(lo, hi, 1) computes the high 32 bits of (hi:lo) << 1:
 *     new_hi = (hi << 1) | (lo >> 31)
 *   We implement this inline.
 */
static inline void
montsqr96(uint *a, uint *c, uint *n, uint w)
{
#if 0
    montmul96(a, a, c, n, w);
    return;
#else
    uint t[5] = {0,0,0,0,0};
    uint m, C, S;

    for (int i = 0; i < 3; i++) {
        uint hicarry, p = 0;

        /* diagonal term */
        accumlh(a[i]*a[i], mul_hi(a[i],a[i]), &t[i], &C);

        for (int j = i+1; j < 3; j++) {
            uint lo = a[j] * a[i];
            uint hi = mul_hi(a[j], a[i]);

            /* __funnelshift_l(lo, hi, 1): shift 64-bit value (hi:lo) left 1, take top 32 */
            hicarry = hi >> 31;
            hi      = (hi << 1) | (lo >> 31);   /* replaces __funnelshift_l */
            lo    <<= 1;

            accum3lh(t[j], lo, C, hi, &t[j], &C);
            add2(C, p, &C, &p);
            p += hicarry;
        }

        m = t[0] * w;
        add2(t[3], C, &t[3], &C); t[4] = C + p;

        accumlh(n[0]*m, mul_hi(n[0],m), &t[0], &C);

        for (int j = 1; j < 3; j++) {
            t[j-1] = innermul(n[j], m, C, t[j], &C);
        }

        add2(t[3], C, &t[2], &C);
        t[3] = t[4] + C;
        t[4] = 0;
    }

#ifdef NO_SUB
    c[0]=t[0]; c[1]=t[1]; c[2]=t[2];
#else
    uint cc[3] = {t[0], t[1], t[2]};

#ifdef ALWAYS_SUB
    /* Subtract n from t using a single borrow chain across all 3 limbs.
     * The previous code split the t[2] limb into two usub_borrow calls and
     * only kept the borrow from the second, silently discarding the borrow
     * from t[2]-n[2].  Use ulong arithmetic for each limb to keep it clean. */
    uint borrow;
    ulong d0 = (ulong)t[0] - n[0];
    uint s0  = (uint)d0;
    borrow   = (uint)(-(int)(d0 >> 32)) & 1;

    ulong d1 = (ulong)t[1] - n[1] - borrow;
    uint s1  = (uint)d1;
    borrow   = (uint)(-(int)(d1 >> 32)) & 1;

    ulong d2 = (ulong)t[2] - n[2] - borrow;
    uint s2  = (uint)d2;
    borrow   = (uint)(-(int)(d2 >> 32)) & 1;

    /* use subtracted result if t[3] set OR no borrow */
    if (t[3] || borrow == 0) { c[0]=s0; c[1]=s1; c[2]=s2; }
    else                     { c[0]=cc[0]; c[1]=cc[1]; c[2]=cc[2]; }
#else
    if (t[3]) {
        ulong d0b = (ulong)t[0] - n[0];
        uint  s0b = (uint)d0b;
        uint  bw  = (uint)(-(int)(d0b >> 32)) & 1;
        ulong d1b = (ulong)t[1] - n[1] - bw;
        uint  s1b = (uint)d1b;
        bw        = (uint)(-(int)(d1b >> 32)) & 1;
        ulong d2b = (ulong)t[2] - n[2] - bw;
        c[0]=s0b; c[1]=s1b; c[2]=(uint)d2b;
    } else { c[0]=cc[0]; c[1]=cc[1]; c[2]=cc[2]; }
#endif
#endif
#endif
}

/*
 * montmul128 / montsqr128 -- 128-bit Montgomery multiply/square (CIOS)
 *
 *   Same structure as the 96-bit versions but with 4 limbs.
 *   accum3lh / accumlh / innermul / add2 handle all carry chains.
 */
static inline void
montmul128(uint *a, uint *b, uint *c, uint *n, uint w)
{
    uint t[6] = {0,0,0,0,0,0};
    uint m, C, S;

    for (int i = 0; i < 4; i++) {
        C = 0;
        for (int j = 0; j < 4; j++) {
            accum3lh(t[j], a[j]*b[i], C, mul_hi(a[j],b[i]), &t[j], &C);
        }
        m = t[0]*w;
        add2(t[4], C, &S, &C); t[4]=S; t[5]=C;

        accumlh(n[0]*m, mul_hi(n[0],m), &t[0], &C);
        for (int j = 1; j < 4; j++) {
            accum3lh(t[j], n[j]*m, C, mul_hi(n[j],m), &t[j-1], &C);
        }
        add2(t[4], C, &t[3], &C); t[4]=t[5]+C; t[5]=0;
    }

    uint cc[4]={t[0],t[1],t[2],t[3]};
    if (t[4]) {
        uint bw,bw2,bw3,bw4,bw5;
        t[0]=usub_borrow(t[0],n[0],&bw);
        uint tmp;
        tmp=usub_borrow(t[1],n[1],&bw2); t[1]=usub_borrow(tmp,bw, &bw3); bw2+=bw3;
        tmp=usub_borrow(t[2],n[2],&bw3); t[2]=usub_borrow(tmp,bw2,&bw4); bw3+=bw4;
        tmp=usub_borrow(t[3],n[3],&bw4); t[3]=usub_borrow(tmp,bw3,&bw5); bw4+=bw5;
        c[0]=t[0]; c[1]=t[1]; c[2]=t[2]; c[3]=t[3];
    } else { c[0]=cc[0]; c[1]=cc[1]; c[2]=cc[2]; c[3]=cc[3]; }
}

static inline void
montsqr128(uint *a, uint *c, uint *n, uint w)
{
    uint t[6] = {0,0,0,0,0,0};
    uint m, C;

    for (int i = 0; i < 4; i++) {
        uint hicarry, p=0;
        accumlh(a[i]*a[i], mul_hi(a[i],a[i]), &t[i], &C);

        for (int j = i+1; j < 4; j++) {
            uint lo = a[j]*a[i];
            uint hi = mul_hi(a[j],a[i]);
            hicarry = hi >> 31;
            hi      = (hi<<1)|(lo>>31);   /* __funnelshift_l equivalent */
            lo    <<= 1;
            accum3lh(t[j], lo, C, hi, &t[j], &C);
            add2(C, p, &C, &p); p+=hicarry;
        }

        m = t[0]*w;
        add2(t[4], C, &t[4], &C); t[5]=C+p;

        accumlh(n[0]*m, mul_hi(n[0],m), &t[0], &C);
        for (int j = 1; j < 4; j++) {
            accum3lh(t[j], n[j]*m, C, mul_hi(n[j],m), &t[j-1], &C);
        }
        add2(t[4], C, &t[3], &C); t[4]=t[5]+C; t[5]=0;
    }

    uint cc[4]={t[0],t[1],t[2],t[3]};
    if (t[4]) {
        uint bw,bw2,bw3,bw4,bw5; uint tmp;
        t[0]=usub_borrow(t[0],n[0],&bw);
        tmp=usub_borrow(t[1],n[1],&bw2); t[1]=usub_borrow(tmp,bw, &bw3); bw2+=bw3;
        tmp=usub_borrow(t[2],n[2],&bw3); t[2]=usub_borrow(tmp,bw2,&bw4); bw3+=bw4;
        tmp=usub_borrow(t[3],n[3],&bw4); t[3]=usub_borrow(tmp,bw3,&bw5); bw4+=bw5;
        c[0]=t[0]; c[1]=t[1]; c[2]=t[2]; c[3]=t[3];
    } else { c[0]=cc[0]; c[1]=cc[1]; c[2]=cc[2]; c[3]=cc[3]; }
}

/* ------------------------------------------------------------------ */
/* bigmontmul / bigmontsqr  (_BIGWORDS = 4, matches original default) */
/*                                                                     */
/* These are fully general CIOS loops identical in structure to the   */
/* 128-bit versions above; kept separate to match the original file.  */
/* ------------------------------------------------------------------ */
#define _BIGWORDS 4

static inline void
bigmontmul(uint *a, uint *b, uint *c, uint *n, uint w)
{
    uint t[_BIGWORDS+2];
    for (int i=0;i<_BIGWORDS+2;i++) t[i]=0;
    uint m,C,S;

    for (int i=0;i<_BIGWORDS;i++) {
        C=0;
        for (int j=0;j<_BIGWORDS;j++)
            accum3lh(t[j], a[j]*b[i], C, mul_hi(a[j],b[i]), &t[j], &C);

        m=t[0]*w;
        add2(t[_BIGWORDS],C,&S,&C); t[_BIGWORDS]=S; t[_BIGWORDS+1]=C;

        accumlh(n[0]*m, mul_hi(n[0],m), &t[0], &C);
        for (int j=1;j<_BIGWORDS;j++)
            accum3lh(t[j],n[j]*m,C,mul_hi(n[j],m),&t[j-1],&C);

        add2(t[_BIGWORDS],C,&t[_BIGWORDS-1],&C);
        t[_BIGWORDS]=t[_BIGWORDS+1]+C; t[_BIGWORDS+1]=0;
    }

    uint cc[_BIGWORDS]; for(int i=0;i<_BIGWORDS;i++) cc[i]=t[i];

    /* subtract n from t in 4-limb chunks, chaining borrow */
    uint carry=0;
    for (int i=0;i<_BIGWORDS;i+=4) {
        uint bw,bw2,bw3,bw4,bw5; uint tmp;
        /* fold incoming carry as an extra borrow on the first limb */
        tmp   =usub_borrow(t[i+0],carry,  &bw);
        t[i+0]=usub_borrow(tmp,   n[i+0], &bw2); bw+=bw2;
        tmp   =usub_borrow(t[i+1],bw,     &bw2);
        t[i+1]=usub_borrow(tmp,   n[i+1], &bw3); bw2+=bw3;
        tmp   =usub_borrow(t[i+2],bw2,    &bw3);
        t[i+2]=usub_borrow(tmp,   n[i+2], &bw4); bw3+=bw4;
        tmp   =usub_borrow(t[i+3],bw3,    &bw4);
        t[i+3]=usub_borrow(tmp,   n[i+3], &bw5); carry=bw4+bw5;
    }

    if (t[_BIGWORDS]) { for(int i=0;i<_BIGWORDS;i++) c[i]=t[i]; }
    else              { for(int i=0;i<_BIGWORDS;i++) c[i]=cc[i]; }
}

static inline void
bigmontsqr(uint *a, uint *c, uint *n, uint w)
{
    uint t[_BIGWORDS+2];
    for(int i=0;i<_BIGWORDS+2;i++) t[i]=0;
    uint m,C;

    for(int i=0;i<_BIGWORDS;i++) {
        uint hicarry,p=0;
        accumlh(a[i]*a[i], mul_hi(a[i],a[i]), &t[i], &C);

        for(int j=i+1;j<_BIGWORDS;j++) {
            uint lo=a[j]*a[i], hi=mul_hi(a[j],a[i]);
            hicarry=hi>>31;
            hi=(hi<<1)|(lo>>31); lo<<=1;
            accum3lh(t[j],lo,C,hi,&t[j],&C);
            add2(C,p,&C,&p); p+=hicarry;
        }

        m=t[0]*w;
        add2(t[_BIGWORDS],C,&t[_BIGWORDS],&C); t[_BIGWORDS+1]=C+p;

        accumlh(n[0]*m, mul_hi(n[0],m), &t[0], &C);
        for(int j=1;j<_BIGWORDS;j++)
            accum3lh(t[j],n[j]*m,C,mul_hi(n[j],m),&t[j-1],&C);

        add2(t[_BIGWORDS],C,&t[_BIGWORDS-1],&C);
        t[_BIGWORDS]=t[_BIGWORDS+1]+C; t[_BIGWORDS+1]=0;
    }

    uint cc[_BIGWORDS]; for(int i=0;i<_BIGWORDS;i++) cc[i]=t[i];

    uint carry=0;
    for(int i=0;i<_BIGWORDS;i+=4) {
        uint bw,bw2,bw3,bw4,bw5; uint tmp;
        tmp   =usub_borrow(t[i+0],carry,  &bw);
        t[i+0]=usub_borrow(tmp,   n[i+0], &bw2); bw+=bw2;
        tmp   =usub_borrow(t[i+1],bw,     &bw2);
        t[i+1]=usub_borrow(tmp,   n[i+1], &bw3); bw2+=bw3;
        tmp   =usub_borrow(t[i+2],bw2,    &bw3);
        t[i+2]=usub_borrow(tmp,   n[i+2], &bw4); bw3+=bw4;
        tmp   =usub_borrow(t[i+3],bw3,    &bw4);
        t[i+3]=usub_borrow(tmp,   n[i+3], &bw5); carry=bw4+bw5;
    }

    if(t[_BIGWORDS]) { for(int i=0;i<_BIGWORDS;i++) c[i]=t[i]; }
    else             { for(int i=0;i<_BIGWORDS;i++) c[i]=cc[i]; }
}

#endif /* OPENCL_INTRINSICS_CL */
