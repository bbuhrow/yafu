# AVX-ECM stage 1 on twisted Edwards curves: background and plan

Status (2026-11-09, branch `wip_edwards`): Phase 0 done, the gate passes
(Part C). Phase 1 done (Part D): `-avxecm_curve 1` builds family (iv)
curves and runs their Montgomery image through the existing PRAC stage 1
and both stage 2s, lane-exact against the oracle on IFMA / AVX-512F / AVX2.
Part A is the
background (what the three papers say, how Edwards curves work, which
curves to use, and how they connect to the existing stage 2). Part B is the
implementation plan for IFMA, AVX-512F and AVX2 (and the RNS backend above
the crossover, which comes for free).

Sources (PDFs in the repo root):

* [EECM] Bernstein, Birkner, Lange, Peters, *ECM using Edwards curves*
  (preprint 2008-01-09; `eecm-20080109.pdf`). The journal version (Math.
  Comp. 82, 2013) adds more material, including the a = -1 curves; this
  preprint does not have it.
* [BK] Bos, Kleinjung, *ECM at Work* (Asiacrypt 2012; `bos_ecm_at_work.pdf`).
* [GKL] Gélin, Kleinjung, Lenstra, *Parametrizations for families of
  ECM-friendly curves* (`kleinjung_parameterizations_for_families_of_ecm_friendly_curves.pdf`).
* Cited by them and used here: Hisil, Wong, Carter, Dawson, *Twisted Edwards
  curves revisited* (2008: extended coordinates); Bernstein, Birkner, Lange,
  *Starfish on strike* (2010: a = -1 curves for ECM); Barbulescu, Bos,
  Bouvier, Kleinjung, Montgomery, *Finding ECM-friendly curves through a
  study of Galois properties* (ANTS 2012).

---

# Part A: background

## A1. What decides ECM's speed

Stage 1 computes Q = [s]P with s = prod q^k (q^k < B1, as yafu builds it).
log2 s ~ 1.4427 B1. Two things decide the time per factor found:

1. **Cost per bit of s.** This is the number of modular multiplications per
   bit. The curve shape and the coordinates decide it.
2. **Probability per curve.** This is the chance that #E(F_p) (more
   precisely ord P) is B1-smooth apart from one prime below B2. The
   rational torsion of the curve family decides it, along with subtler
   Galois effects that raise the average power of 2 (and 3) dividing
   #E(F_p).

Edwards curves win on (1). The best Edwards families roughly tie or slightly
beat the Suyama curves on (2).

## A2. Where yafu is today (Montgomery, x-only, PRAC)

* Curve: By^2 = x^3 + Ax^2 + x from Suyama's sigma (`vec_build_one_curve`).
  Group order is divisible by 12. A24 = (A+2)/4 is a full-size residue.
* Point: (X : Z) only, so P and -P look the same. There is no general
  addition, only *differential* addition (P+Q needs P-Q).
* `vec_add` (differential add) = 4M + 2S; `vec_duplicate` = 3M + 2S
  (2M + 2S + 1D when A24 is small, which it isn't for Suyama).
* Scalar multiplication per prime with PRAC (Montgomery's Lucas-chain
  heuristic). GMP-ECM counted 1.39 differential adds + 0.135 doublings per
  bit of s at B1 = 1e6 [EECM §3.4]. That is ~9.0 multiplications per bit, or
  ~8.3 "mul-equivalents" if a template squaring costs 0.78 of a
  multiplication (IFMA, `avxecm_ifma_analysis.md` §4).
* Measured: c1 = 14.2-15.5 template mulmods per unit of B1
  (`avxecm_stage2_fft_design.md` §13.3), i.e. ~10-10.7 per bit. The
  difference from 8.3 is add/sub, copies and PRAC bookkeeping.

The x-only representation is the reason PRAC is needed: without y you
cannot add arbitrary multiples, so windowing (precompute 3P, 5P, ... and
add them) is impossible. Differential chains need ~1.4 additions per bit no
matter how clever they are.

## A3. Edwards curves

**Edwards form** [EECM §2]: x^2 + y^2 = 1 + d x^2 y^2, d not in {0, 1}.

* Addition law: (x1,y1) + (x2,y2) = ((x1y2 + y1x2)/(1 + d x1x2y1y2),
  (y1y2 - x1x2)/(1 - d x1x2y1y2)). It is the same formula for doubling
  ("unified"). If d is a non-square it is *complete* (no exceptions at all).
* Neutral element (0, 1). Negation -(x, y) = (-x, y): **y is invariant
  under negation**, which plays the role Montgomery's x plays.
* (0, -1) has order 2; (±1, 0) have order 4. So every Edwards curve has a
  rational point of order 4, and hence 4 | #E(F_p).

**Twisted Edwards** [EECM §2.3]: a x^2 + y^2 = 1 + d x^2 y^2. The addition
law has a x1x2 in place of x1x2. Twisted Edwards curves are exactly the
curves that are birationally equivalent to Montgomery curves:

    A = 2(a + d)/(a - d),  B = 4/(a - d),   u = (1 + y)/(1 - y),  v = u/x.

This map is a group isomorphism. That is the key fact for the stage-2
hand-off (§A7).

**Coordinates and costs** (M = multiplication, S = squaring, D =
multiplication by a curve constant):

| system | add | doubling | notes |
|---|---|---|---|
| Montgomery (X:Z), yafu | 4M+2S (differential) | 3M+2S (full A24) | no windows possible |
| projective Edwards (X:Y:Z) [EECM §2.1] | 10M+1S+1D | 3M+4S | |
| inverted Edwards (Z/X, Z/Y) [EECM §2.2] | 9M+1S+1D | 3M+4S+1D | GMP-EECM 2008 |
| **extended, a = -1** (X:Y:Z:T), T = XY/Z [Hisil et al.] | **8M** (cached operand), 7M without output T | **3M+4S**, +1M for T | Starfish, EECM-MPFQ, [BK], [GKL] |

The a = -1 extended formulas, as EFD writes them (dbl-2008-hwcd,
add-2008-hwcd-3):

    dbl:  A = X^2, B = Y^2, C = 2Z^2, D = -A
          E = (X+Y)^2 - A - B, G = D + B, F = G - C, H = D - B
          X3 = E F, Y3 = G H, Z3 = F G   [T3 = E H only if an add follows]
    add (operand 2 cached as (Y2-X2, Y2+X2, 2d T2, 2 Z2)):
          A = (Y1-X1)(Y2-X2), B = (Y1+X1)(Y2+X2), C = T1 (2d T2), D = Z1 (2 Z2)
          E = B - A, F = D - C, G = D + C, H = B + A
          X3 = E F, Y3 = G H, Z3 = F G   [T3 = E H only if an add follows]

Hisil's mixing trick [BK §5]: a doubling followed by a doubling skips T
(3M+4S). A doubling followed by an addition computes T (4M+4S). An addition
is always followed by a doubling, so it skips T3 (7M). So **each addition
costs 8M in total** (7M + the 1M of T in the preceding doubling). a = -1
makes D = -A free, and multiplication by 2d is absorbed into the cached
table entry, so d may be a full-size residue at no cost.

## A4. Scalar multiplication: why Edwards wins

With full additions, stage 1 can use **signed sliding windows / wNAF**:
precompute the odd multiples P, 3P, ..., (2^(w-1) - 1)P (negatives are free:
negate x, i.e. swap Y-X and Y+X and negate 2dT), then per bit do one
doubling and, for a fraction eps ~ 1/(w+1) of bits, one addition.

    Edwards cost per bit ~ (3M + 4S) + eps (8M)

GMP-EECM's other change [EECM §3.5] is to **batch primes**: multiply many
q^k into one big integer t and recode t, instead of running one short chain
per prime. Windows only pay on long scalars. The table costs 2^(w-2)
additions per batch, so the batches must be long (the optimum w grows with
the batch length L: eps(w) ~ 1/(w+1) + 2^(w-2)/L).

Per-bit estimate, mul-equivalents, sqr = 0.78 M (template, IFMA) or 1.0 M
(RNS, where a squaring costs about the same as a multiplication):

| method | sqr = 0.78 | sqr = 1.0 |
|---|---:|---:|
| PRAC on Montgomery, full A24 (yafu now) | 8.33 | 9.0 |
| Montgomery ladder, small A24 and small x0 (GMP-ECM param 1 batch) | ~8.1 | ~9.0 |
| Edwards a = -1 extended, w = 6 (eps ~ 0.14) | 7.25 | 8.13 |
| Edwards a = -1 extended, w = 8 (eps ~ 0.11) | 7.01 | 7.89 |
| Edwards a = -1 extended, w = 10 (eps ~ 0.095) | 6.88 | 7.76 |

**Expected stage-1 gain: ~1.15-1.2x on the radix (template) path, ~1.1-1.15x
on RNS**, before add/sub overhead (Phase 0 measured composed operations
instead: 1.20-1.30x on the radix path, Part C2) (similar for both: ~7 add/sub per Edwards
doubling, ~5-6 per PRAC step). [BK] Table 3 agrees: at B1 = 8192, EECM-MPFQ
uses 93040 M+S against GMP-ECM's 104428 (1.12x counting S = M).

A property that matters for SIMD: **the scalar is the same for all 8 lanes
(and all threads)**. Recoding is done once, and the control flow is
identical across lanes. Only the data differs. Window lookups are
uniform-index loads, not gathers.

**Bos-Kleinjung chains** [BK]: for a *fixed* small B1 they searched over
10^12 integers that have cheap addition-subtraction chains (D doublings
plus very few additions of the input point only, the "no-storage" form
2^D ± 2^n1 ± ...) for B1-powersmoothness. They then greedily covered
lcm(1..B1) with them. Results vs windowed EECM: 2-9% fewer M+S (9% at
B1 = 256, 2% at 8192), and 10 residues of memory instead of 38-550. That
matters on GPUs (55x less memory at B1 = 8192, 2.2x faster than windows on
a GTX 580). On a CPU with a big L2 and B1 >= 1e4 the gain shrinks to a
couple of percent, and it needs a precomputed chain table per B1.
Dixon-Lenstra-style batching (low-weight products of primes) is the
light-weight version of the same idea.

## A5. Torsion: which curves to use

Mazur: E_tors(Q) is Z/m (m = 1..10, 12) or Z/2 x Z/2m (m = 1..4). Rational
torsion of order T maps into E(F_p) for good p, so T | #E(F_p). Bigger T
means smoother group orders.

* Edwards curves always have Z/4, so the possible groups are Z/4, Z/8, Z/12,
  Z/2 x Z/4, Z/2 x Z/8 [EECM §4].
* **a = 1 families with large torsion** [EECM §4-6]:
  - Z/2 x Z/8: x8 = (u^2 + 2u + 2)/(u^2 - 2), d = (2x8^2 - 1)/x8^4 (Thm 4.5);
    non-torsion points from the Atkin-Morain construction (Thm 5.2) or from
    a search for small ones (§6.1).
  - Z/12: x3 = (u^2 - 1)/(u^2 + 1), y3 = -(u - 1)^2/(u^2 + 1),
    d = (u^2 + 1)^3 (u^2 - 4u + 1)/((u - 1)^6 (u + 1)^2) (Thm 4.9).
  - These need a = 1, so the 8M a = -1 addition is not available: there is
    an extra D, or 9M.
* **a = -1 restricts the torsion** to Z/6, Z/8 or Z/2 x Z/4 among the
  interesting ones (order > 4). *Starfish* compared 1000 curves of each type
  and found that the faster a = -1 arithmetic outweighs the smaller torsion.
  Its best curve was a Z/6 curve "C6": -x^2 + y^2 = 1 - (13312/18225) x^2 y^2.
* **Galois families** (Barbulescu et al., via [GKL] §2): the a = -1,
  Z/2 x Z/4 curves are -x^2 + y^2 = 1 - e^4 x^2 y^2 (e not in {0, ±1}). For
  special e the *average* 2-adic valuation of #E(F_p) is higher than the
  generic 14/3:
  - (i) e = g^2, (ii) e = g^2/2, (iii) e^2 - 1 a square: 29/6;
  - **(iv) e^2 + 1 a square, i.e. e = (g - 1/g)/2: 16/3.** Here
    16 | #E(F_p) always, and 32 | #E(F_p) when p = 1 mod 4 and g(g-1)(g+1)
    is a QR mod p.
* **[GKL] gives infinite families** for (i)-(iv). For the best family (iv)
  (Corollary 3.6, j = 4): take points (x, y) = [m]Q on the rank-1 curve
  E': y^2 = x^3 - x^2 - 9x + 9, Q = (5, 8); t = (4x + 4)/(y - 4); then

      e  = (t^2 + 4t)/(t^2 - 4)
      x1 = (2t^3 + 2t^2 - 8t - 8)/(t^4 + 6t^3 + 12t^2 + 16t)
      y1 = (t^6 + 6t^5 + 10t^4 - 16t^3 - 48t^2 - 32t - 32)
           /(t^6 + 6t^5 + 10t^4 + 16t^3 + 48t^2 + 64t)

  gives -x^2 + y^2 = 1 - e^4 x^2 y^2 with the non-torsion point (x1, y1).
  m = 1 gives e = 15/8, point (28/195, 3152/3495). (The PDF text was
  extracted by machine. Recheck these rational functions symbolically before
  coding: Phase 0.)
* **Yields** ([GKL] Tables 1-2, EECM's test: count the b-bit primes each
  curve finds, b = 15..63): family (iv) curves C[m] beat C6 by ~1.0% on
  average, the Z/12 a = 1 curve C12 by ~1.4%, and the older family (iv) curve
  C_{2x4} by ~0.6%. The 100 curves C[1..100] behave alike (best/worst within
  1-9% at tiny b, 0.5% at b = 26). **So m works like sigma: any m >= 1 gives
  an equally good curve.**
* **The simplest generic option** (Thm 3.1, j = 1, any rational t not in
  {0, ±1, ±3, ±1/3}): e = 3(t^2 - 1)/(8t), x1 = 128t^3/(27t^6 + 63t^4 -
  63t^2 - 27), y1 = (9t^4 - 2t^2 + 9)/(9t^4 - 9). It needs no auxiliary
  elliptic curve, but it gives the generic 14/3 2-adic valuation.

Compared with what yafu runs now: Suyama curves (12 | #E) are in the same
class as the Z/12 curves. Expect family (iv) to be **about equal, maybe 1-2%
better, per curve**. The gain is the cheaper arithmetic, not the
probability. (Phase 0 measured more: family (iv) finds 6.5% ± 1.0% more
factors per curve than yafu's Suyama curves, j = 1 3.6% fewer; Part C1.)

Note: GMP-EECM's other 2008 trick, *small* a, d and base-point
coordinates, saves multiplications by tiny constants. It does not carry
over here: each lane is a different curve, so constants are full-size,
lane-varying residues. (Montgomery param 1 is different: A24 = d/2^64 is
small, and a lane-varying single-word constant still allows `vecmulmod_1`.)

## A6. Exceptions

The 8M extended addition is not unified with doubling. It fails when the
two inputs are equal (or differ by certain 2/4-torsion points). In a wNAF
chain that happens modulo p only if a small multiple of the running point
is the identity mod p. The failed formula then produces zero coordinates mod
p, which the final gcd catches. The effect is an early stage-1 find (or N if
it happens mod every prime, which is vanishingly rare). The lane-exact
oracle (Phase 0) will show whether it ever matters.

The x-only PRAC stage 1 has an exceptional case of its own (found in
Phase 1, Part D2). A differential addition is only correct while the
difference point is not the identity. Modulo p, once a point's remaining
order r is a small prime, a PRAC chain for a later prime q > r can meet a
difference [k]P with r | k. The chain then collapses to Z = 0, and stage 1
finds p although [s]P != O mod p. Example: r = 61 with 61^2 | ord P, B1 =
3000, so s contains only 61^1, collapsing in the chain for q = 647. This
happens on any curve family, Suyama included. It costs nothing and is
occasionally a free find. An exact oracle has to replay the chains
(`edw_ref.py prac_stage1_z`); the order alone does not predict it.

## A7. Hand-off to the existing stage 2

Stage 2 (pairing and FFT, with or without Brent-Suyama) needs only:

* the Montgomery x-coordinate of Q as (X : Z) in `P->X`, `P->Z`, and
* A24 = (A+2)/4 in `work->s`.

It never uses B or y (Brent-Suyama rebuilds a y-carrying model as
g y^2 = x^3 + A x^2 + x, which is the curve or its twist, so B is irrelevant).

Through the birational map of §A3, with Q = (X : Y : Z : T) in extended
coordinates:

    X_M = Z + Y,   Z_M = Z - Y            (u = (1 + y)/(1 - y))
    A24 = (A + 2)/4 = a/(a - d)
        = 1/(1 + d)       for a = -1
        = 1/(1 - e^4)     for the d = -e^4 families

* That is 2 add/sub per lane at the end of stage 1 plus one constant per
  curve, computed with mpz at build time. **Stage 2 does not change.** The
  FFT continuation, its memory model, auto-B2 and Brent-Suyama all apply
  as they are.
* It is exact: the map is a group isomorphism between E and the Montgomery
  curve, so x_M([k]Q) = x_M([k]Q') exactly when y([k]Q) = y([k]Q'). The
  stage-2 condition (jD = ±k mod ord Q) is unchanged. Edwards' own y-only
  arithmetic would just be Montgomery's x-only arithmetic in other
  coordinates (u = (1+y)/(1-y)), so nothing is gained by a native Edwards
  stage 2. [EECM §3] also expects stage 2's curve operations to be minor;
  ours are ~10% of the FFT stage 2 (roots).
* Stage-1 gcd: the identity is (0 : 1 : 1) and the 2-torsion point is
  (0 : -1 : 1), so gcd(X, N) catches both. gcd(Z_M, N) = gcd(Z - Y, N)
  catches only the identity. Take the gcd of X: a point of order 2 at the
  end of stage 1 is a find.
* 1 + d must be invertible mod N. If it isn't, gcd(1 + d, N) is a factor.
  The same holds for every inversion in curve construction.

---

# Part B: implementation plan

All curve arithmetic goes through the existing pointer API
(`vecmulmod_ptr`, `vecsqrmod_ptr`, `vecaddmod_ptr`, `vecsubmod_ptr`,
`vecaddsubmod_ptr`). So the template cores (IFMA, AVX-512F FP, AVX2 FP),
exact word sizing and the RNS backend all come for free. No new kernels are
needed in the first version.

## Phase 0: reference model and decision gate (Python, no C)

1. **Symbolic check** (sympy) of every formula used: the a = -1 extended
   dbl/add, the cached-operand form, the Edwards -> Montgomery map and
   A24 = 1/(1+d), and the [GKL] j = 1 and j = 4 parametrizations (point on
   the curve, e^2 + 1 a square for j = 4, m = 1 gives e = 15/8).
2. **Exact Edwards reference** (extend `factor/avx-ecm/bench/ref/ecm_ref.py`
   and `arith/vecpoly/ecm_tests/bs_oracle.py`): build the curve from m (or
   t) mod p, compute ord P mod p by BSGS on the Montgomery image, apply
   yafu's s (q^k < B1), and predict stage-1 and stage-2 finds lane by lane
   for given D, dF, blocks and Brent-Suyama polynomial. That gives the same
   lane-for-lane test that validated the FFT stage 2 (§14.3 of the stage-2
   doc).
3. **Op-count model:** for B1 = 1e3..1e8, the batch length L, window w,
   chunking and table cost, count M, S and add/sub per bit. Price them with
   measured per-ISA template mul/sqr/add ns (`ae_bench -b`, plus the AVX2
   and AVX-512F sqr/mul ratios, which are not in the docs yet) and the RNS
   ratios. Predict stage-1 time vs PRAC per ISA and size.
4. **Yield sanity check:** with the oracle, count finds over ~10^4 random
   (p, curve) pairs at p ~ 2^30..2^40: Suyama vs family (iv) vs j = 1.
   Expect family (iv) >= Suyama within noise. This confirms the curve choice
   on yafu's own s and B2.

**Gate:** proceed if the predicted stage-1 gain is >= ~8% on IFMA radix and
the yields are not worse.

## Phase 1: curve construction (mpz, per lane)

* `vec_build_one_curve_edwards(tdata, X, Y, Z, T, dd, a24, m)` next to
  `vec_build_one_curve`. Family (iv): compute [m]Q on E' mod N (affine,
  mpz, ~2 log2 m inversions; microseconds per curve), then t, e, d = -e^4,
  (x1, y1) mod N, then (X : Y : Z : T) = (x1, y1, 1, x1 y1), 2d, and
  A24 = 1/(1 + d). Convert everything to the backend's Montgomery form with
  `vec_monty_mulR` / `vec_lane_set`.
* Any failed inversion is gcd > 1 with N: report it as a find (stage 0) and
  give the lane a dummy curve, as the stage-2 code already does with
  `found[]`.
* "sigma" becomes m (family (iv)) or t (j = 1); random m from the existing
  LCG; fixed m for tests (consecutive per lane, like sigma today). Also keep
  j = 1 (generic, cheapest to build) and possibly C6-style Z/6 for A/B.
* Option, e.g. `-avxecm_curve` (cmdOptions; `ecm_obj` field): 0 = Suyama
  Montgomery (today), 1 = Edwards family (iv), 2 = Edwards j = 1. Default
  stays 0 until Phase 6. (Changing structs in factor.h / cmdOptions.h needs
  `make clean`; edit CRLF files with the Edit tool.)

## Phase 2: point arithmetic

* `ecm_ept { X, Y, Z, T }` and a cached form `{ YmX, YpX, T2d, Z2 }`.
  Allocate in `vec_ecm_work_init` with `max(NWORDS, RNS words)` like the
  other vectors.
* `edw_dbl` (3M+4S), `edw_dbl_t` (4M+4S), `edw_add_cached` (7M, no T3;
  8M variant for table building), `edw_to_cached` (1M for 2d T plus
  add/sub), `edw_neg_cached` (swap YmX and YpX, negate T2d: free when
  selecting).
* **Bounds audit** (both backends):
  - Radix: `vecaddmod52` / `vecsubmod52` assume AMM outputs < ~2N with
    >= 3 spare bits. E = (X+Y)^2 - A - B, F = G - C and H = -(A + B) are
    two-deep add/sub chains. Check them with an `ae_bench -a`-style test on
    real Edwards sequences.
  - RNS: inputs to a mul must be < 8kN, and every submod adds 2kN. The
    worst chain (E: sqr output minus two mul outputs) is < 5kN. Write it
    down the way `handoff_vector_arith.md` §4.1 does for `CROSS_PRODUCT`.
  - Prefer forms that need no explicit negation (fold -1 signs into which
    of G/H is subtracted from which; projective sign changes are free).
* Optional later: dual-product kernels. The doubling's 4 squarings / 3
  final products are independent pairs, which helps the latency-bound
  small sizes (S <= 7 words, `avxecm_ifma_analysis.md` §5).

## Phase 3: stage-1 driver (batched scalar, signed windows)

* Per prime range (the existing `PRIME_RANGE` loop in `vececm`): build the
  scalar in **chunks** of ~L bits (multiply q^k into an mpz until it
  reaches L; L ~ 2^14..2^20, from the model). Recode each chunk once, in
  the main thread, to a signed-window digit string (w-NAF or signed sliding
  window), shared read-only by all threads. Powers of 2 become plain
  doublings.
* Per chunk and thread: build the table of odd multiples P, 3P, ...,
  (2^(w-1) - 1)P in cached form (1 dbl + 2^(w-2) adds), then walk the
  digits: dbl or dbl_t depending on whether the next step adds, then
  add_cached with ±table[digit]. Software-prefetch the next digit's entry,
  which is known in advance.
* Choose w per run from eps(w) = 1/(w+1) + 2^(w-2)/L under a table-memory
  cap: bytes = 2^(w-2) x 4 x NWORDS x 64. At 1024 bits, w = 9 is 128 x 5 KB
  = 640 KB per thread, and a 1-2 MB per-thread cap is reasonable (more
  with RNS's larger vectors).
* Keep the large-B1 `save_b1_intermediate` checkpoints (convert to
  Montgomery, §A7, at each checkpoint), the progress output and
  `stg1Add` / `stg1Doub`-style counters.
* `vec_ecm_stage1` dispatches on the curve type; the PRAC path stays
  intact.

## Phase 4: end of stage 1 and stage-2 hand-off

* gcd on X per lane (§A7) for stage-1 finds, using the existing
  `vec_check_factor` loop on the Edwards X.
* Convert: `P->X = Z + Y`, `P->Z = Z - Y`, `work->s = A24` (set at build
  time). The pairing stage 2, the FFT stage 2 and Brent-Suyama then run
  unchanged.
* Reporting and resume:
  - Finds print the curve as e.g. `param=edw4, m=...`, enough to rebuild
    it.
  - `save_b1` lines: GMP-ECM's resume needs the Montgomery A (or a param
    it knows) and X. Write `A=` (A = 4 A24 - 2) plus the converted X.
  - Check the exact resume syntax GMP-ECM accepts against
    `ecm-master.tar.gz` (in the repo root) before relying on it.
* The auto-B2 model uses T1 = c1 B1 mu (`s2m.c1`). Edwards needs its own
  c1 per ISA (measured in Phase 6). Otherwise `-avxecm_stg2 2` picks B2
  for the old stage-1 speed.

## Phase 5: verification

* Phase 0 oracle vs yafu, fixed m, lane for lane: stage-1 and stage-2 finds
  at 400, 520, 1024, 1250, 2048 (RNS) and 4000 bits, on IFMA, AVX-512F and
  AVX2 clean builds (check with objdump: no zmm in AVX2, no vpmadd52 in
  AVX-512F). Plain, FFT and Brent-Suyama stage 2.
* p15 / p20 end-to-end finds (64 curves, B1 3000 / 11000) at the §4.6 size
  list of `handoff_vector_arith.md`, including word boundaries and the RNS
  crossover. 4 threads.
* `ae_bench -a`-style add/sub range test on Edwards chains at 1-26 spare
  bits (it must fail only below 3).
* MSVC / clang-cl compile of the new code.

## Phase 6: tuning, measurement, default

* Stage-1 time, Edwards vs PRAC, per ISA, at 200-4000 bits and
  B1 = 1e4..1e7. Tune w, L and the memory cap. Measure with 1 thread and
  with all threads (table memory and bandwidth).
* Refit c1 per ISA for the auto-B2 model.
* Optional (expensive): curves-to-find a known p35/p40, Suyama vs family
  (iv). Only the Phase 0 yield check is required.
* Decide the default (`-avxecm_curve 1` if it wins everywhere; keep
  Suyama for any sizes where it doesn't).

## Phase 7 (optional, later)

* **Bos-Kleinjung no-storage chains** for small fixed B1 (the yafu
  `factor()` early B1 levels; tinyecm/microecm are separate). [BK] shows
  2-9% over windows at B1 <= 8192 with tiny memory. It needs an offline
  chain search or their published chains [BK ref 7]. Low priority on CPUs.
* **Montgomery re-parameterization** (the user's original second item):
  GMP-ECM's param 1/3 (small d, so A24 is a one-word constant, x0 = 2, and
  a batch Montgomery ladder over all of s) saves about 1M per doubling over
  Suyama. yafu has an unfinished `PARAM1` path (`build_one_curve_param1`,
  `vec_ecm_stage1_batch`, `vecmulmod_1`; radix only). Per the table in
  §A4 it does not beat Edwards. It is only worth finishing if Edwards loses
  somewhere (e.g. RNS sizes, where squarings are not cheaper). Check the
  torsion / probability of each GMP-ECM param in `ecm-master.tar.gz`
  first.
* Double-base chains with Edwards tripling: the journal EECM found little
  gain over signed sliding windows. Not planned.

## Expected outcome

Stage 1 ~1.15-1.2x faster on the radix path (IFMA, AVX-512F, AVX2), and
~1.1x above the RNS crossover. Per-curve probability equal or up to ~1-2%
better (family (iv)). Stage 2 is unchanged. Since stage 1 is 60-75% of the
time with the FFT stage 2 at its best B2, time per factor drops ~10-15%,
and the auto-B2 choice shifts slightly toward a larger B2 (stage 2 becomes
relatively more expensive).

(Updated by Phase 0, Part C3: stage 1 1.20-1.30x, +6.5% finds per curve,
~20-25% less time per factor.)

---

# Part C: Phase 0 results (2026-11-09)

Tools, all in `arith/edwards/` (Python 3, no yafu changes):

* `edw_symbolic.py`: sympy identities (needs sympy).
* `edw_ref.py`: exact arithmetic mod p. `selftest`; `lanes fam p B1 D dF
  blocks S m0 nbatch`, the lane oracle (same arguments and output as
  `bs_oracle.py lanes`, plus a stage-0 column for curve-build finds; `fam`
  = suyama | iv | j1); `yield bits B1 B2 ncurves [seed]`.
* `edw_model.py [B1 ...]`: stage-1 cost model from `ops_<isa>.txt`.
* `ops_ifma.txt`, `ops_noifma.txt`, `ops_avx2.txt`: `ae_bench -e` output.
  The new `-e` mode in `factor/avx-ecm/bench/ae_bench.c` times mul, sqr,
  add, sub and the composed point operations, ns per 8-lane call.

## C0. Formula checks

* `edw_symbolic.py`, all identities hold:
  - the Edwards -> Montgomery map, and (A+2)/4 = a/(a-d) = 1/(1-e^4) for
    a = -1, d = -e^4;
  - GKL Thm 3.5 j = 1 and j = 4 points on their curves for all t;
  - family (iv): e^2 + 1 is the paper's square on E'
    (y^2 = x^3 - x^2 - 9x + 9, t = (4x+4)/(y-4));
  - family (i): e = (3x/y)^2;
  - m = 1 gives e = 15/8 and the point (28/195, 3152/3495).
  So the machine-extracted formulas in §A5 are right.
* `edw_ref.py selftest` (2379 checks over 40- and 61-bit p):
  - the extended a = -1 doubling (signs folded, see below) and the cached
    addition, run inside a width-w NAF scalar multiply (w = 2..6), agree
    with the affine Edwards law;
  - T = XY/Z holds whenever T is produced;
  - the Montgomery x-only ladder with A24 = 1/(1+d) started from
    u = (1+y)/(1-y) gives exactly (Z+Y : Z-Y) of the Edwards multiple, so
    the stage-2 hand-off of §A7 is exact;
  - point counts on small p: 16 | #E for family (iv), 8 | #E for j1.
* The doubling as checked (and as timed below) folds the a = -1 signs:
  H' = A + B, E' = H' - (X+Y)^2, G' = A - B, F' = 2Z^2 + G'. These are
  E, F, G, H negated, so X3 = E'F', Y3 = G'H', Z3 = F'G', T3 = E'H' are
  unchanged. That is 6 add/sub and no negations.
* The lane oracle reproduces `bs_oracle.py` exactly for Suyama curves.

## C1. Yield: finds per curve, yafu's s and plain stage 2

One random curve of each family per random prime p (paired), exact ord P
mod p. Stage 1 finds when [s]P = O (Edwards: or (0,-1), X = 0); stage 2
when the remaining order is a prime <= B2 (no D-structure).

| bits | B1 | B2 | primes | Suyama | family (iv) | j1 | (iv) - Suyama |
|---:|---:|---:|---:|---:|---:|---:|---|
| 32 | 400 | 4e4 | 40000 | 9876 | 10409 (1.054) | 9109 (0.922) | +533 ± 123 |
| 40 | 2000 | 2e5 | 20000 | 3551 | 3838 (1.081) | 3572 (1.006) | +287 ± 77 |
| 48 | 1e4 | 1e6 | 12000 | 1664 | 1825 (1.097) | 1568 (0.942) | +161 ± 54 |

Pooled: **family (iv) +6.5% ± 1.0%** over yafu's Suyama curves; j1 -3.6%
± 1.5%. That is more than the literature suggested: [GKL] measured +1.4%
against an Edwards Z/12 curve, not against Suyama.

Average valuations of ord P: v2 2.33 (iv) vs 1.44 (Suyama), v3 0.49 vs
1.17. Family (iv) trades 3s for 2s and comes out ahead. j1 (generic
Z/2 x Z/4) is worse than Suyama, so the family matters. **Use family (iv)
only.**

## C2. Stage-1 cost per bit of s

Composed-operation timings (`ae_bench -e`, Rocket Lake, `taskset -c 3`,
ns per 8-lane call; selected rows):

| ISA | S | mul | sqr | add | sub | PRAC add (4M+2S+6a) | PRAC dbl (3M+2S+4a) | Edw dbl (3M+4S+6a) | Edw dbl+T | Edw add (7M+6a) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IFMA | 8 | 57.7 | 47.4 | 7.3 | 23.3 | 430 | 343 | 443 | 502 | 503 |
| IFMA | 20 | 344 | 265 | 23.1 | 70.9 | 2222 | 1777 | 2336 | 2686 | 2716 |
| AVX-512F | 20 | 702 | 495 | 23.2 | 71.6 | 4115 | 3278 | 4337 | 5052 | 5229 |
| AVX2 | 20 | 997 | 880 | 24.3 | 63.5 | 6136 | 5050 | 6840 | 7882 | 7274 |

`edw_model.py` counts the operations on yafu's real scalar:

* **PRAC:** it replays `vec_prac` (`ORIG_PRAC` is off, so conditions 3, 4,
  5, 9; v chosen by `lucas_cost`). At B1 = 1e6 that is 1.374 additions and
  0.151 doublings per bit (B1 = 11000: 1.322 / 0.185).
* **Edwards:** a width-w NAF of the whole scalar, which gives
  1/(w+1) additions per bit, plus the table.

**Calibration (IFMA):** the PRAC model gives 643 / 3302 ns per bit at
S = 8 / 20, B1 = 1e6. yafu measured 689 / 3398 (stage 1 0.994 s / 4.90 s,
8 curves). So the model is 3-7% low; the difference is loop overhead,
copies and pointer swaps, which an Edwards loop has too.

Edwards vs PRAC, ns per bit at B1 = 1e6 (ratio = PRAC / Edwards):

| ISA | S | PRAC | Edw w=8 | Edw w=10 | ratio w=8 | ratio w=10 | table at w=10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| IFMA | 8 | 643 | 554 | 542 | 1.16 | 1.19 | 0.5 MB |
| IFMA | 12 | 1300 | 1042 | 1018 | 1.25 | 1.28 | 0.75 MB |
| IFMA | 20 | 3302 | 2698 | 2636 | 1.22 | 1.25 | 1.25 MB |
| IFMA | 32 | 8149 | 6700 | 6545 | 1.22 | 1.25 | 2 MB |
| AVX-512F | 8 | 1050 | 865 | 846 | 1.21 | 1.24 | 0.5 MB |
| AVX-512F | 20 | 6113 | 4984 | 4863 | 1.23 | 1.26 | 1.25 MB |
| AVX2 | 8 | 2094 | 1773 | 1734 | 1.18 | 1.21 | 0.5 MB |
| AVX2 | 20 | 9178 | 7810 | 7643 | 1.18 | 1.20 | 1.25 MB |

Every size from S = 4 to 40 is in the script output: ratios are 1.15-1.35
at w = 8-10. At B1 = 11000, w = 10 is still the best (ratios 1.15-1.28).

Notes:

* Larger w keeps winning on paper (w = 13 gives 1.21-1.35), but the table
  becomes 2-20 MB per thread. The timings assume a hot table operand, so
  w will be capped by a per-thread memory limit (~1-2 MB, Phase 3) and
  tuned for real.
* Edwards gains more than the plan's M/S count suggested (1.20-1.30x, not
  1.15-1.2x), for two reasons:
  - yafu's simplified PRAC needs 1.37 additions per bit (GMP-ECM's full
    PRAC: 1.39 additions + 0.135 doublings);
  - PRAC also pays 6 add/sub per addition.
* AVX2 gains the least: its squaring is relatively expensive (0.88 of a
  mul at S = 20, vs 0.77 on IFMA).
* RNS (above ~1170 / 1490 / 1540 bits) was not measured. With sqr ~ mul
  the count model gives ~1.15x there.

## C3. Gate

Both conditions pass by a wide margin:

* predicted stage-1 gain >= 8%: measured model **1.20-1.30x** (IFMA
  1.19-1.28x, AVX-512F 1.24-1.26x, AVX2 1.20-1.21x at w = 10);
* yields not worse: **+6.5%**.

Time per factor, with T2 ~ 0.45 T1 (the FFT stage 2 near its best B2):
(1 + 0.45) / (1/1.25 + 0.45) = 1.16x from speed, times 1.065 from yield,
**~1.24x less time per factor** on IFMA. A bit less on AVX2, and once
auto-B2 is re-optimized for the faster stage 1, a bit more.

Changes to the plan from Phase 0:

* Use family (iv) only. j1 stays as a test curve, not an option.
* Window width w ~ 8-10 under a per-thread table memory cap.
* Phase 1 must skip m whose t is in the excluded set {-2, -1, ±4} or hits
  y = 4 on E' (t infinite). Those are curve-build "finds" mod p, which the
  oracle already reports as stage 0.
* Phase 5 runs the lane-exact comparison with `edw_ref.py lanes iv ...`
  (an `oracle_cmp.sh` variant for `-avxecm_curve 1`).

---

# Part D: Phase 1 results (2026-11-09)

## D1. What was built

Phase 1 builds the family (iv) curves and hands their Montgomery image
(§A7) to the *existing* x-only code. That is a complete, testable step on
its own: stage 1 is still PRAC, but the curves, and so the +6.5% yield of
C1, are already in use. Phases 2-3 then only replace the stage-1
arithmetic.

* Option `-avxecm_curve 0|1` (`OptionArray[144]`, `options->avxecm_curve`,
  `ecm_obj.avx_curve`, `thread_data_t.curve_type`); also readable from
  `yafu.ini`.
  - Default `AVXECM_CURVE_DEFAULT` = 0 (Suyama) in `include/factor.h`.
  - Constants `AVXECM_CURVE_SUYAMA` / `AVXECM_CURVE_EDW4` in `avx_ecm.h`.
* `vec_build_one_curve_edwards` (`avxecm.c`, beside `vec_build_one_curve`),
  per lane, mpz mod N:
  - [m]Q on E' left to right (`edw_wstep`, affine; one inversion per step);
  - t, e, y1 (Horner);
  - A24 = 1/(1 - e^4) and u = (1 + y1)/(1 - y1);
  - then X = u, Z = 1, `work->s` = A24 in Montgomery form, exactly where
    the Suyama builder puts its values.
* Failed inversions:
  - If a denominator shares a factor g with N (1 < g < N), the lane becomes
    (1 : g). That is the identity modulo g's primes, it stays there through
    PRAC, and the stage-1 gcd reports it.
  - If a denominator is 0 mod N, m moves to m + 1. The reported sigma is
    the m actually used.
* Random m: the LCG, in [2, 2^31). `-sigma s` gives lanes s, s+1, ..., as
  for Suyama.
* `-v` prints the curve family.
* `save_b1` lines for Edwards curves write `A=` (the Montgomery
  A = 4 A24 - 2) instead of `SIGMA=`, plus X as before
  (`vec_ecm_save_curve`). GMP-ECM's resume reads `A=`
  (`ecm-master/resume.c`).
* x1 is not computed yet (the x-only image does not need it); Phase 2 adds
  it.

## D2. Verification

`arith/edwards/edw_oracle_cmp.sh` (like `oracle_cmp.sh`: yafu
`-avxecm_curve 1 -avxecm_stg2 1` vs `edw_ref.py lanes iv-mont`, lane for
lane, stage 1 and stage 2). Every row has 0 batches differing, on IFMA,
AVX-512F (no `vpmadd52` in the binary) and AVX2 (no zmm):

| input | B1 / B2 | batches | stage-1 / stage-2 finds |
|---|---|---:|---|
| 400-bit, p40, plain | 3000 / 231284 | 30 | 6 / 50 |
| same, Brent-Suyama D_6 | | 15 | 4 / 31-33 (dF differs per ISA) |
| 400-bit, p32 | | 15 | 24 / 68 |
| 1024-bit, p40 | | 12 | 4 / 21 |
| 2048-bit (RNS backend), p40 | | 8 | 2 / 5 |

Also:

* The Suyama default still matches `bs_oracle.py`.
* Random m with 4 threads finds a p16 at 520 bits, with the pairing and the
  FFT stage 2.
* GMP-ECM 7.0.6 `-resume` of yafu's 8 `save_b1` lines (B1 3000, B2 231284)
  finds p on exactly the lanes whose remaining order is a prime <= B2
  (lanes 1, 4, 5). So `A=` and X are right.

One oracle fix came out of this: **PRAC can find p in stage 1 when the
order alone says it can't.**

* The case: the point's remaining order r is a small prime (here 61, with
  61^2 | ord P but only 61 in s).
* What happens: a later PRAC chain (q = 647) meets an intermediate
  difference that is a multiple of r, i.e. the identity mod p. The x-only
  differential addition then collapses to Z = 0.
* The fix: `edw_ref.py` now replays yafu's PRAC chains mod p
  (`prac_stage1_z`) whenever 1 < r < B1, for Suyama and `iv-mont` alike.
* This is a property of x-only PRAC, not of the curves: a ladder by the same
  s does not find it. The Edwards stage 1 of Phase 3 has no differential
  additions, so it will not have it (its own exceptional cases, §A6, are
  different).

## D3. Next

Phase 2 (extended-coordinate point arithmetic, add/sub bounds audit) and
Phase 3 (signed-window driver). Phase 4's `save_b1` and stage-2 hand-off
parts are already done: stage 1 will convert back with X_M = Z + Y,
Z_M = Z - Y and keep `work->s` = A24.
