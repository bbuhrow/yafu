# AVX-ECM stage 2: FFT continuation (polynomial multipoint evaluation)

Status (2026-11-08, `wip_arith`, uncommitted): **first implementation
working** on IFMA, AVX-512F and AVX2; opt-in with `-avxecm_stg2 1`.
Results and open items in §10.
References: Montgomery, *An FFT Extension of the Elliptic Curve Method of
Factorization* (dissertation, 1992), ch. 3-5; Montgomery, *Speeding the
Pollard and Elliptic Curve Methods of Factorization* (1987), §4-5;
GMP-ECM 7 (`stage2.c`, `polyeval.c`, `ecm_ntt.c`, `mpzspv.c`, `bestd.c`).

## 1. Goal and constraints

AVX-ECM's stage 2 today is the baby-step/giant-step *pairing* continuation
(`vec_ecm_stage2_init` / `vec_ecm_stage2_pair` in `avxecm.c`): about one
mulmod per prime in (B1, B2], so B2 is held to 100*B1. The FFT continuation
tests d1*d2 differences with O(d log^2 d) work, which makes GMP-ECM-sized
B2 (B2 ~ (11/6 B1)^1.43, about 1000*B1 at B1 = 1e6) affordable and raises
the per-curve success probability.

Constraints (from the user):
* Keep AVX-ECM's paradigm: **one curve per vector lane**, 8 curves per
  `vec_bignum_t` (`VECLEN` = 8 in every build; AVX2 runs two 4-lane halves).
* Targets: AVX-512 with IFMA, AVX-512F without IFMA, AVX2 (+FMA).
* Reuse the new base arithmetic: the constant-S template Montgomery cores
  (`arith/vecmont`) and the 52-bit product layer (`arith/rns/rns_vec52.h`:
  `LO52`/`HI52`, native IFMA or exact FP emulation, 512 or 256 bits).

## 2. The algorithm

Q = stage-1 output on a Montgomery curve (x-only, X:Z). For p | N let q be
the order of Q mod p. Pick D (even, product of small primes) and

* S1 = { k : 1 <= k < D/2, gcd(k, D) = 1 },  |S1| = phi(D)/2 <= dF = 2^t,
* baby roots  a_k = x([k]Q)  (affine),   F(X) = prod_k (X - a_k), padded to
  degree dF with duplicate roots (harmless: a duplicate only repeats a test),
* giant roots b_j = x([jD]Q), j = j0, j0+1, ...; block i uses dF of them:
  G_i(X) = prod_{j in block i} (X - b_j),
* H = prod_i G_i  mod F,
* result: V = prod_k H(a_k) (mod N), one value per lane; gcd(V, N).

x([jD]Q) = x([k]Q) mod p iff q | jD - k or q | jD + k, so every prime
q in (B1, B2] coprime to D is covered by j = round(q/D) and k = |q - jD|.
Coverage per block = dF * D; blocks = ceil((B2 - B1) / (dF*D)), j0 =
round(B1/D). D's prime factors are < 30 < B1 so every q > B1 is coprime.

No Brent-Suyama extension in v1 (exponent 1). Dickson / x^E roots are a
later option (§9).

### 2.1 Roots

Computed in the radix-2^52 Montgomery domain with the template cores
(independent of the curve backend, which may be RNS above the crossover):

* Q, A24 = (A+2)/4 come from the backend through `vec_lane_get` (backend
  Montgomery form, R_b), converted per lane to the radix domain
  (x * R_b^-1 * R mod N, R_b read off the backend's Montgomery one).
* Baby: [2]Q by xDBL, then [k+2]Q = xADD([k]Q, [2]Q; [k-2]Q) for odd k up to
  D/2: D/4 xADDs (6 mulmod each), keep k in S1.
* Giant: [D]Q (ladder), [j0 D]Q and [(j0+1) D]Q (ladders), then
  [(j+1)D]Q = xADD([jD]Q, [D]Q; [(j-1)D]Q): 6 mulmod per root.
* Normalization: affine x = X/Z by Montgomery's batch-inversion trick per
  lane over a block (3 mulmod per root + 1 mpz_invert per lane per block),
  then X * Zinv (1 mulmod). A failed inversion means gcd(Z, N) > 1: that lane
  has found a factor; its non-invertible product is folded into the lane's
  result (so the final gcd reports it) and the lane continues with garbage.

Per root: ~10 mulmods. Per block: dF of them.

### 2.2 Polynomial work per block (dominant cost)

1. G_i from its dF roots: product tree, log2(dF) levels (§4.2).
2. H := H * G_i mod F: one product of length 2dF, then a Barrett-style
   reduction with precomputed transforms of invF = floor(X^(2dF-1)/F) and F.
First block: H := G_1 mod F (G_1 - F, both monic of degree dF).

### 2.3 Once per curve batch

* F from the baby roots, **keeping the product tree** (log2(dF) levels of dF
  coefficients) for the evaluation at the end.
* invF by Newton iteration; forward transforms of F and invF kept.
* Evaluation of H at the dF roots of F by the transposed (Tellegen /
  "scaled remainder tree") algorithm, §4.4: one middle product with invF
  at the root, two middle products per tree node going down, no other
  inversions. Product of the dF values into the lane result.

## 3. Arithmetic: polynomials over Z/NZ, one curve per lane

### 3.1 Coefficients

A coefficient is 8 lanes x n words (n = NWORDS of the radix domain, radix
2^52, Montgomery form wrt R = 2^(52n)), exactly the `vec_bignum_t` layout:
word w of lane l at `c[w*8 + l]`. A polynomial is an array of coefficients
(stride 8n u64), so a coefficient pointer can be wrapped in a `vec_bignum_t`
and passed to the template mul. Invariant: **every stored coefficient is in
[0, 2N)** (template mul maps [0,2N)^2 -> [0,2N) because N < R/8; the
coefficient add/sub here reduce mod 2N; the CRT output is < 2N).

### 3.2 Products: NTT over K word-size primes, then CRT (Montgomery ch. 3.4)

Coefficient products of polynomials of length L are computed as cyclic
convolutions modulo K primes p_i in (2^49, 2^50), p_i = 1 mod 2^20, and
recombined by the CRT. P = prod p_i > 4 L (2N)^2, so K = ceil((2*52n + 2 +
log2 L + 2) / 49).

Every lane holds a different curve but uses the **same primes, twiddles and
constants**: the NTT of an 8-lane polynomial is 8 independent NTTs done with
broadcast constants and no shuffles, at every stage. NTT-form layout: prime
i, index t at `v[(i*L + t)*8 + lane]`.

Residue arithmetic uses only `LO52`/`HI52` (52x52-bit products, low/high 52
bits; IFMA `vpmadd52luq/huq` natively, exact FP emulation otherwise), so one
kernel source (`vp_kern.h`) builds for IFMA / AVX-512F / AVX2 like
`rns_yaob_kern.h`. Values are lazy in [0, 2p) or [0, 4p) (< 2^52).

* Shoup mulmod by a constant w (w' = floor(w 2^52 / p)), x < 2^52:
  q = hi52(x w'), r = lo52(x w) - lo52(q p) mod 2^52, r in [0, 2p).
* Montgomery mulmod (pointwise products), a, b < 2p: t = a b < 2^100,
  m = lo52(lo52(t) * (-p^-1 mod 2^52)), r = hi52(t) + hi52(m p) + [lo52(t) != 0]
  = a b 2^-52 mod p, r in [0, 2p).
* Forward NTT: decimation in frequency, Gentleman-Sande butterflies (Harvey
  lazy form), natural order in, bit-reversed out, values in [0, 2p).
* Inverse NTT: decimation in time, Cooley-Tukey butterflies, bit-reversed
  in, natural out, values in [0, 4p). The 1/L, the 2^52 of the Montgomery
  pointwise product and the CRT's (P/p_i)^-1 are one Shoup constant per
  (prime, L), applied in the CRT.
* Twiddles: per prime, per level 2m (m = 1, 2, ..., Lmax/2), m entries
  w_{2m}^j and their Shoup companions; any L <= Lmax uses the same levels.

**Conversion to residues** (x < 2^(52n), digits d_k): with
c_{k,i} = 2^(52k+104) mod p_i, accumulate L = sum lo52(d_k c), H = sum
hi52(d_k c) (exact in 64 bits for n < 2^11), then a 3-limb Montgomery REDC
by 2^104 gives x mod p_i in [0, 2p): 2n + 6 products per prime.

**CRT back** (Montgomery 3.4.4): y_i = r_i * (L^-1 2^52 (P/p_i)^-1) mod p_i
(Shoup), kappa = floor(sum y_i / p_i + 1/8) in FP (h/P < 1/4 gives the
margin), then S = sum_i y_i C_i + kappa C_K with C_i = (P/p_i) R^-1 2^104
mod N and C_K = -P R^-1 2^104 mod N, accumulated in n+2 words with
LO52/HI52, and two Montgomery word-reductions (divide by 2^104): output
h R^-1 mod N, < 2N. So an NTT product of Montgomery-form inputs is a
Montgomery-form output, exactly like the template mul. 2(K+1)n + 4n + 3K
products per coefficient.

### 3.3 Small products

Below a threshold degree (tunable, start 16), polynomial products use
schoolbook multiplication with the template mul (d^2 mulmods). Tree leaves
(degree 1 x degree 1) are one mul and one add.

## 4. Polynomial algorithms (all lengths powers of 2)

Monic polynomials are stored without their leading 1.

### 4.1 Product of two monic polynomials of degree d
(X^d + a)(X^d + b) = X^2d + X^d (a + b) + a b, deg(ab) <= 2d - 2: one cyclic
convolution of length 2d.

### 4.2 Product tree
Level 0: the roots, negated (X - r). Level s: pairs of monic degree-2^(s-1)
polynomials multiplied (§4.1). F keeps level 0 and every k-th level after
it (`tree[s]`, dF coefficients each, k = 1..4 chosen by the memory
model, §11); a level that is not kept is rebuilt from the nearest kept
level below during the evaluation. G keeps none and is built in place.
§4.1 works in place (c = a, b = a + d): a + b is saved first, the
product's inputs are consumed by the forward transforms before the CRT
writes the output.

### 4.3 Reduction modulo F
For P = H G (deg <= 2dF - 1): Q = floor(P_hi * invF / X^(dF-1)) with
P_hi = floor(P / X^dF), invF = floor(X^(2dF-1)/F); R = P_lo - (Q f mod
X^dF) where F = X^dF + f. Verified against a scalar reference (§7).
As built, R is computed modulo X^dF - 1 so the product by f is a wrapped
(cyclic) product of length dF instead of 2dF:
P = H G + X^dF H = Q F + R gives, mod X^dF - 1,
R = (P_lo + P_hi) - Q - (Q f mod X^dF - 1). Steps: transform H (all K
primes) and G (4 primes at a time, multiplied in), inverse, CRT P_hi into
a scratch polynomial and P_lo over H; then Q from P_hi with the length-2dF
transform of invF; then the length-dF cyclic product Q f. The transforms of
invF (length 2dF) and f (length dF) are stored, or recomputed per block in
the lean mode (§11).
invF by Newton: I_{2k} from I_k, k = 1, 2, ..., dF (Montgomery fig. 3.5.1 /
GMP-ECM `PolyInvert`); the product (A mod X^2k) S_k is only needed in
coefficients k..2k-1, which a cyclic product of length 2k gives exactly.

### 4.4 Evaluation (transposed / scaled remainder tree)
For a node M = A B (degrees dA, dB), let U_M[s] (s = 1..deg M) be the
coefficient of X^-s in the Laurent series H/M. Then
* root: U_F[s] = (H * invF)[2dF - 1 - s], s = 1..dF (a middle product:
  1/F = sum_{r >= dF} f'_r X^-r and invF holds f'_dF .. f'_{2dF-1}
  reversed; implemented as a product of length 2dF with the stored
  transform of invF and a reversed slice);
* children: U_A[s] = U_M[s + dB] + sum_{r<dB} b_r U_M[s + r] (s = 1..dA),
  a middle product with the sibling B (cyclic length deg M suffices), and
  U_B likewise with A;
* leaf X - a: U[1] = H(a).
Cost ~ 2 products per node (vs 1 for the tree build), no inversions.
As built the evaluation runs in place in H: a node transforms U_M once
(first half of the NTT buffer), then forms U_A and U_B one after the other
in the second half; U_M's upper half is the only extra input and stays in
place until U_B is done.

## 5. Parameters

dF = 2^t; D from GMP-ECM's power-of-2 table (phi(D)/2 <= dF, densest):

| t | dF | D | phi(D)/2 | D/dF |
|---:|---:|---:|---:|---:|
| 8 | 256 | 2310 | 240 | 9.0 |
| 9 | 512 | 4620 | 480 | 9.0 |
| 10 | 1024 | 9240 | 960 | 9.0 |
| 11 | 2048 | 19110 | 2016 | 9.3 |
| 12 | 4096 | 39270 | 3840 | 9.6 |
| 13 | 8192 | 79170 | 8064 | 9.7 |
| 14 | 16384 | 158340 | 16128 | 9.7 |

Block coverage dF*D (t = 12: 1.6e8). Default B2 as GMP-ECM:
(11/6 * B1)^1.43; t chosen from B2 (about sqrt((B2-B1)/9.5 / blocks_target))
and capped by a per-thread memory budget. Memory: §11 (as built, about
64 dF [(t/k + 5.5) n + 5K] bytes per thread).

The pairing stage 2 stays for B2 <= 100*B1-ish and as the fallback
(option to force either).

## 6. Code layout

| path | contents |
|---|---|
| `arith/vecpoly/vecpoly.h` | API: context (primes, twiddles, CRT tables; shared, read-only), per-thread work, polynomial ops |
| `arith/vecpoly/vp_setup.c` | prime search, tables, scalar references for tests |
| `arith/vecpoly/vp_kern.h` | ISA-generic kernels on `rns_vec52.h`: to/from residues, NTT, pointwise, coefficient add/sub |
| `arith/vecpoly/vp_avx512.c`, `vp_avx2.c` | kernel instantiations (512 / 256 bits) |
| `arith/vecpoly/vp_poly.c` | products, trees, reciprocal, mulmod F, evaluation |
| `arith/vecpoly/vp_test.c`, `Makefile` | standalone tests vs GMP and benchmarks |
| `factor/avx-ecm/avxecm_stage2fft.c` | ECM side: parameters, roots, driver, hook in `vececm` |

## 7. Verification plan

* Kernels vs scalar GMP: conversion, NTT round trip, cyclic products of
  random polynomials (every lane, several n and L), CRT exactness at the
  extreme input bound (all coefficients 2N-1).
* Polynomial ops vs GMP per lane: product tree root vs direct expansion,
  reciprocal, mulmod F, evaluation vs Horner at every root.
* ECM: the found-factor oracle. Known group orders: pick sigma, p with
  known largest/second-largest factors of #E(F_p), check stage 2 finds p iff
  the largest factor q of the order (after B1) lies in (B1, B2]; compare
  lane-for-lane with the pairing stage 2 at B2 = 100*B1 (FFT stage 2 must
  find a superset), and with GMP-ECM's `-sigma 0:x` / `-param` runs where
  parameterizations match.
* All three ISAs (IFMA, `NOIFMA`, `AVX2`) on Rocket Lake.

## 8. Cost model (512-bit N, n = 10, K = 22, IFMA)

Mulmod (template) ~ 300 vector ops; Shoup butterfly ~ 6; conversion in
~ 2n+15 per prime; CRT ~ 2(K+1)n. Per coefficient per tree level:
~3 log2(2s)/2 K butterflies + conversions ~ 10-15 mulmod-equivalents.
Per block (dF = 4096): roots ~ 41k mulmods, G tree ~ 0.5M, mulmod F ~ 0.2M
mulmod-equivalents for dF*D = 1.6e8 of B2 range, vs ~1 mulmod per prime
(~8.6e6 per 1.6e8) for pairing: ~10x cheaper per unit of B2.

## 9. Later

* Brent-Suyama (Dickson polynomials) for more chances per root.
* Native-FP Shoup butterflies for non-IFMA targets (the LO52/HI52 emulation
  costs ~3x the native ops).
* Cache blocking of the NTT for 2dF*64 B > L2; truncated / wrapped products
  (length dF instead of 2dF in the reduction); keep transforms between tree
  levels where possible.
* RNS-domain roots above the crossover (roots are ~10% of the work).
* Multithreading one batch across threads (today each thread owns 8 curves).

## 10. Status and results (2026-11-08)

### Code (as built)

* `arith/vecpoly/`: `vecpoly.h`, `vp_setup.c`, `vp_kern.c` (one source, the
  width picked from `RNS_HAVE_AVX512` / `RNS_HAVE_AVX2`), `vp_poly.c`,
  `vp_test.c`, `Makefile` (`make`, `make NOIFMA=1`, `make AVX2=1`).
* `factor/avx-ecm/avxecm_stage2fft.[ch]`: parameters, x-only point
  arithmetic in the vecpoly domain, batch normalization, the driver.
* Hook: thread phase 4 in `vececm` (`vec_ecm_stage2_fft_work_fcn`); shared
  setup / free in `vec_ecm_main`; `thread_data_t.s2fft`.
* Options: `avxecm_stg2` (0 pairing, 1 FFT, 2 FFT with B2 from the memory
  limit (§11.5); `ecm_obj.avx_stg2_fft`, default
  `AVXECM_STG2_FFT_DEFAULT` = 0) and `avxecm_stg2mem` (MB per thread,
  `ecm_obj.avx_stg2_mem`, default 256). With the FFT stage 2 and no
  explicit B2, `ecm.c` uses GMP-ECM's default B2 = (11/6 B1)^1.43. The
  progress line prints the real B2 (the covered B2 is rounded up to whole
  blocks).
* Makefile `ECM_SRCS` (+ `vp_kern` in `FP_KERNEL_OBJS`), both VS ecm
  projects (+ filters).

### Verification

* `vp_test`: every polynomial operation, every lane, against GMP at 200,
  400, 520, 1024, 2048 bits (products up to 512 x 512 including the
  all-(2N-1) bound, tree, reciprocal, mulmod, evaluation for m = 1..512,
  schoolbook and NTT paths): pass with IFMA, `NOIFMA=1` (no `vpmadd52` in
  the kernels) and `AVX2=1` (no zmm).
* ECM, fixed sigmas, lane for lane (`arith/vecpoly/ecm_tests/cmp_s2.sh`: pairing vs FFT at the same
  requested B2; every pairing stage-2 find must be an FFT find):

| input | backend | IFMA | AVX-512F | AVX2 |
|---|---|---|---|---|
| 400-bit, p40, B1 3000, B2 231284, 240 curves | radix | 47/47 (+1) | 47/47 (+1) | 47/47 (+1) |
| 1024-bit, p40, 96 curves | radix | | | 23/23 (+1) |
| 2048-bit, p40, 64 curves | RNS | 7/7 | | 7/7 |
| 2^439-1 (p27), B1 1000, B2 60000, 160 curves | radix | | | 93/93 (+4) |

  (+k: FFT-only finds; its covered B2 is a little larger.) The AVX2
  column was run with a binary that had AVX-512 vecpoly kernels linked in
  (see the correction under Speed); §11.4 repeats every row with clean
  builds of all three ISAs. 4 threads:
  identical finds. Mersenne arithmetic is never selected in this tree
  (`forceNoMersenne` = 1); the code reads the modulus and Montgomery one
  generically, so it would work there too, but that is untested.
* yafu builds clean (new files warning-free under `-Wall -Wconversion`) for
  avx2 / avx512 / ifma; clang-cl and MSVC v143 `cl` compile the new files
  for avx2 / avx512 / avx512ifma (the full VS ecm project stops earlier on
  a pre-existing missing `ecm.h` include path in `avxppm1.c`).

### Speed (Rocket Lake, 1 thread, 8 curves, B1 = 1e6, semiprime inputs; `ecm_tests/bench_s2.sh`, inputs from `ecm_tests/mkn.py`)

After the memory trim (§11; re-measured, all binaries from clean builds):

| bits | ISA | stage 1 | pairing, B2 = 1e8 | FFT, B2 = 1.2e8 | FFT, B2 = 9.7e8 (default) |
|---:|---|---:|---:|---:|---:|
| 400 | IFMA | 1.06 s | 0.90 s | 0.13 s | 0.49 s |
| 400 | AVX-512F | 1.62 s | 1.02 s | 0.19 s | 0.71 s |
| 400 | AVX2 | 3.24 s | 1.59 s | 0.30 s | 1.13 s |
| 1024 | IFMA | 4.96 s | 2.52 s | 0.41 s | 1.55 s (dF 4096, 138 MB) |

So at the same B2 the FFT stage 2 is 5-7x faster, and GMP-ECM's ~10x
larger default B2 costs less than the old 100*B1 pairing stage 2 (AVX2:
about 0.7x of it).

Correction: the first version of this table gave AVX2 as 0.16 s / 0.55 s.
That "AVX2" yafu binary was not a clean build and had AVX-512 (IFMA)
vecpoly kernels linked in (`objdump` shows zmm / `vpmadd52` in `vpk_*`);
the AVX2 rows of the ECM verification table below were also run with it.
`vp_test` built with `make AVX2=1` was always genuine AVX2. The numbers
above and the §11 ECM verification use a clean AVX2 build (no zmm in the
binary).
`vp_test -b` (520 bits, IFMA, dF = 4096): G tree 42 ms, mulmod 37 ms,
F + 1/F ~100 ms, evaluation 86 ms; ~93 mulmod-equivalents per root for the
tree.

### Open items / next steps

1. **Default on?** `AVXECM_STG2_FFT_DEFAULT` is 0 for now. yafu's
   `factor()` curve counts follow GMP-ECM's default B2, so turning it on
   should make AVX-ECM's work per level match them; needs a decision and a
   t-level check.
2. **Memory.** Done (§11): 1024 bits at dF = 4096 went from ~290 MB to
   143 MB per thread (77 MB lean); the setup picks dF, the tree levels
   kept and the lean mode within the limit; `-avxecm_stg2mem` takes
   fractional MB; §11.3 gives the best B2 per memory limit.
3. **Speed.** Native FP butterflies for non-IFMA; cache-blocked / fused
   NTT passes for 2dF*64 B > L2; keep a transform between tree levels;
   tune `school_max` (8) and the cost model in `s2_tcost`; RNS-domain roots
   above the crossover.
4. **Brent-Suyama** (Dickson polynomials, GMP-ECM `-power`/`-dickson`) for
   extra chances per root.
5. Stage 2 for P-1/P+1 (`avxppm1.c`) could reuse vecpoly (P-1 roots are
   just powers); unreachable from yafu today.
6. `ecmtune()` could measure the stage-2 cost ratio and choose B2; see
   §11.5 for choosing B2 from the memory limit (done: `-avxecm_stg2 2`).

## 11. Memory, time and the choice of B2 (2026-11-08)

### 11.1 What uses memory

Per thread, with dF = 2^t, n words per coefficient, c = 64n bytes per
coefficient (8 lanes x n words), 64 bytes per NTT vector and
K = `vp_nprimes(n, t+1)` = floor((104n + t + 3) / 50) + 1 ~ 2.1n primes:

| buffer | bytes | notes |
|---|---|---|
| NTT buffer `w->ntt` | K * 2dF * 64 | the only full-size transform buffer |
| second operand `w->grp` | 4 * 2dF * 64 | 4 primes at a time (`VP_GRP`) |
| scratch `w->tmp[0..1]` | 1.5 * dF * c | P_hi / Q in the reduction; a+b, reversed children |
| G, H | 2 * dF * c | F lives in H until block 0 (H = G_0 - F) |
| F's tree | ceil(t/k) * dF * c | levels 0, k, 2k, ...; k = 1..4 |
| 1/F | dF * c | |
| transforms of 1/F (length 2dF) and F (length dF) | 3K * dF * 64 | full mode only |
| F | dF * c | lean mode only (its transform is recomputed) |
| root batch X, Z, prefix | 3 * 256 * c | batch inversion per 256 roots |
| shared: twiddles | 2 * K * 2dF * 8 | one table of w_2dF^j per prime, levels stride through it |

So, per thread,

    full:  M ~ 64 dF [(t/k + 5.5) n + 5K + 8]  ~ 64 dF n (t/k + 16)
    lean:  M ~ 64 dF [(t/k + 6.5) n + 2K + 8]  ~ 64 dF n (t/k + 11)

It is linear in the modulus size (n) and grows like dF log dF in the
polynomial length. It does not depend on B1 or B2 directly: B2 only
decides how many blocks of dF giant steps are run (time, not memory). The
dependence on B2 (and so on B1, through the default
B2 = (11/6 B1)^1.43) is through the choice of dF: one block covers
dF * D ~ 9.5 dF^2 of the B2 range, so a larger B2 makes a larger dF
worthwhile (the unconstrained optimum is a few blocks, dF ~ sqrt(B2/40)),
and the memory limit caps dF. Also D(t) <= B1 is required (t = 15 needs
B1 >= 324870).

Memory per thread in MB (`ecm_tests/s2model.py`; full with k = 1 / lean
with k = 3, shared tables included):

| bits | dF=512 | 1024 | 2048 | 4096 | 8192 | 16384 | 32768 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 200 | 4/2 | 7/4 | 15/9 | 31/17 | 64/36 | 132/73 | 271/145 |
| 400 | 7/4 | 14/8 | 29/16 | 59/32 | 121/68 | 250/135 | 526/274 |
| 700 | 12/7 | 24/14 | 49/28 | 102/55 | 210/116 | 433/231 | 893/461 |
| 1024 | 17/10 | 34/20 | 70/39 | 143/77 | 296/163 | 610/324 | 1259/647 |
| 2048 | 34/20 | 68/40 | 138/77 | 285/152 | 587/321 | 1212/640 | 2502/1278 |

Before the trim (§10), 1024 bits at dF = 4096 took ~290 MB; now 143 MB
(full) or 77 MB (lean, every 3rd tree level). What changed:

1. One NTT buffer instead of four. The second operand of a full-length
   product is transformed 4 primes at a time into `w->grp` and multiplied
   into the first; products of length <= Lmax/2 use the two halves of the
   buffer. (The per-prime-group conversion re-reads the coefficients; no
   measurable cost.)
2. The reduction's product by f is wrapped to length dF (§4.3): the stored
   transform of F is half as long, and the mulmod got faster (1024 bits,
   dF = 4096: 83 -> 68 ms).
3. Everything in place: G's tree, the evaluation (H becomes the values),
   F built into H; the evaluation's root step uses the stored transform of
   1/F instead of a separate product; the Newton reciprocal uses length-2k
   cyclic products for its middle part (§4.3). Evaluation 203 -> 187 ms,
   reciprocal + transforms 119 -> 84 ms at the same size.
4. Scratch sized to need: 1.5 dF coefficients instead of 6 x 2dF.
5. Roots are normalized 256 at a time (8 extra `mpz_invert` per 256 roots,
   well under 1%) instead of keeping X, Z, prefix products and roots for a
   whole block.
6. Twiddles: one table of Lmax/2 entries per prime instead of one per
   level (shared memory halved).
7. Optional, chosen by the cost model when the limit requires it: keep
   only every k-th tree level (the evaluation rebuilds the others: about
   (k-1)/2 extra tree builds, once per curve batch), and the lean mode (no
   stored transforms of 1/F and F: +2-10% per block, 3K dF 64 bytes
   less).

Not done: splitting the length-2dF products into cyclic + negacyclic
halves would halve `w->ntt` (another ~K dF 64 bytes); keeping the tree on
disk (GMP-ECM `-treefile`).

### 11.2 Time

In template mulmods (mu(n) = one 8-lane mulmod; IFMA: 17 / 59 / 110 / 352 /
1430 ns at n = 4 / 8 / 11 / 20 / 40):

    T2 = mu(n) [ blocks * Cb + Cs ]
    Cb = dF (10 + 1.78 tc)                       giant roots, G's tree, H G mod F
    Cs = dF (17 + 3.7 tc) + rebuilt tree levels  baby roots, F's tree, 1/F, evaluation
    tc(t, n) ~ 8t - 42.8 + 4.05 t^2 / n          G's tree, mulmods per root

(tc fit to `vp_test -b`, n = 4..40, t = 8..13, within ~10%; the lean mode
adds ~0.12 tc per block. `s2_cost` in `avxecm_stage2fft.c` is this model;
it overestimates yafu's stage 2 by ~7%.) One block covers dF * D ~ 9.5 dF^2
of B2, so the cost per unit of B2 is

    mu(n) (10 + 1.78 tc) / D   ~   t / dF

Doubling dF (one more t, about twice the memory) halves the cost per unit
of B2, less ~10%. That is why stage 2 runs faster with more memory even
when the covered B2 is larger: the time is blocks x (cost per block), and
the blocks fall like 1/dF^2 while the cost per block grows like dF log dF.
Measured, 1024 bits, B1 = 1e6, default B2 (~9.4e8), IFMA, 1 thread:

| limit (MB) | dF | blocks | tree levels kept | stage 2 |
|---:|---:|---:|---:|---:|
| 16 | 512 | 382 | 5 of 9 | 5.90 s |
| 32 | 1024 | 96 | 5 of 10 | 3.48 s |
| 64 | 2048 | 24 | 6 of 11 | 2.20 s |
| 128 | 4096 | 6 | 6 of 12 | 1.58 s |
| 256 | 4096 | 6 | 12 of 12 | 1.76 s |

(stage 1: 4.9 s; the 128 / 256 MB difference is noise.) The covered B2 is
also rounded up to a whole block, so a request just above a block boundary
costs a whole extra block.

### 11.3 The most efficient B2 for a memory limit

The quantity to minimize is the expected time per factor,
(T1 + T2) / P(B1, B2), with P the probability that the group order of a
curve modulo p (taken as a random number of size p/23, Suyama torsion) is
B1-smooth apart from one prime in (B1, B2] (Dickman rho). For a fixed
limit, T2 is (nearly) linear in B2 at the largest dF that fits;
`ecm_tests/s2model.py` scans every configuration (t, k, lean, blocks) that
fits, keeps the cheapest stage 2 for each covered B2, and finds the
optimum (T1 = c1 B1 mu(n), c1 = 16.9 at 400 bits and 13.9 at 1024 bits,
measured):

| bits | B1 | p | limit MB | best B2 | T2/T1 | dF | blocks | GMP default B2 | its T2/T1 | its loss | B2 within 2% |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 400 | 1e6 | 35 | 16 | 4.7e8 | 0.35 | 2048 | 12 | 9.4e8 | 0.64 | 4.8% | 3.1e8 .. 7.1e8 |
| 400 | 1e6 | 35 | 64 | 8.1e8 | 0.37 | 4096 | 5 | 9.7e8 | 0.43 | 0.0% | 6.4e8 .. 1.3e9 |
| 400 | 1e6 | 35 | 256 | 2.0e9 | 0.62 | 8192 | 3 | 9.7e8 | 0.43 | 1.3% | 8.1e8 .. 2.6e9 |
| 400 | 3e6 | 40 | 16 | 1.3e9 | 0.28 | 2048 | 33 | 4.4e9 | 0.91 | 15.9% | 8.6e8 .. 2.0e9 |
| 400 | 3e6 | 40 | 64 | 3.9e9 | 0.37 | 8192 | 6 | 4.5e9 | 0.42 | 0.1% | 2.6e9 .. 6.5e9 |
| 400 | 3e6 | 40 | 256 | 7.8e9 | 0.47 | 16384 | 3 | 4.5e9 | 0.37 | 3.3% | 5.2e9 .. 1.0e10 |
| 1024 | 2.5e5 | 30 | 16 | 6.4e7 | 0.37 | 512 | 27 | 1.3e8 | 0.69 | 6.3% | 4.3e7 .. 7.8e7 |
| 1024 | 2.5e5 | 30 | 64 | 1.6e8 | 0.46 | 2048 | 4 | 1.6e8 | 0.46 | 0.0% | 1.2e8 .. 2.4e8 |
| 1024 | 1e6 | 35 | 16 | 2.3e8 | 0.31 | 512 | 97 | 9.1e8 | 1.20 | 24.7% | 1.5e8 .. 2.9e8 |
| 1024 | 1e6 | 35 | 64 | 5.9e8 | 0.32 | 2048 | 15 | 9.4e8 | 0.48 | 2.0% | 3.9e8 .. 8.6e8 |
| 1024 | 1e6 | 35 | 256 | 2.0e9 | 0.51 | 8192 | 3 | 9.7e8 | 0.35 | 2.5% | 1.3e9 .. 2.6e9 |
| 1024 | 3e6 | 40 | 64 | 1.6e9 | 0.27 | 2048 | 41 | 4.4e9 | 0.70 | 9.4% | 1.3e9 .. 2.5e9 |
| 1024 | 3e6 | 40 | 256 | 5.2e9 | 0.34 | 8192 | 8 | 4.5e9 | 0.31 | 0.1% | 3.2e9 .. 7.8e9 |
| 1024 | 3e6 | 40 | 1024 | 1.0e10 | 0.45 | 16384 | 4 | 4.5e9 | 0.30 | 4.8% | 7.8e9 .. 1.3e10 |
| 1024 | 1.1e7 | 45 | 64 | 5.9e9 | 0.26 | 2048 | 151 | 2.8e10 | 1.19 | 28.9% | 3.8e9 .. 9.2e9 |
| 1024 | 1.1e7 | 45 | 256 | 1.8e10 | 0.26 | 8192 | 27 | 2.9e10 | 0.41 | 2.2% | 1.2e10 .. 2.7e10 |
| 1024 | 1.1e7 | 45 | 1024 | 5.3e10 | 0.33 | 32768 | 5 | 3.2e10 | 0.24 | 2.2% | 4.3e10 .. 8.5e10 |

What this says:

1. **It is predictable.** Given n and the limit, the largest dF that fits
   fixes the cost per block (Cb) and per batch (Cs); the best B2 is where
   T2 ~ 0.3-0.6 T1, i.e. blocks ~ (0.4 T1 / mu - Cs) / Cb, rounded to a
   whole block. The optimum is flat (B2 can vary by ~2x within 2% of the
   best) and hardly depends on the factor size targeted: for the same B1,
   p30..p45 give the same B2 to within a block or two.
2. **The best B2 grows with memory**, about 2-3x per 4x memory (when dF
   doubles, a block covers 4x more for 2x the cost).
3. **GMP-ECM's default B2** is within ~5% of the best whenever the limit
   lets dF reach roughly sqrt(B2/40) (at 1024 bits: 64 MB for B1 <= 1e6,
   256 MB up to B1 = 1.1e7). Under a tight limit it is too large: 10-100%
   more time per factor (1024 bits, 16-64 MB, B1 >= 1e6), because T2 then
   exceeds T1. There the right B2 is ~4-5x smaller.
4. With a generous limit the best B2 is larger than GMP-ECM's default
   (2x at 1024 bits, B1 = 1e6, 256 MB), but that only gains 2-5%.

`python3 ecm_tests/s2model.py bits B1 MB` lists, for one case, the
cheapest stage 2 for every covered B2 and its expected time per p35 / p40
factor relative to the best.

### 11.4 Verification of the trimmed code

* `vp_test` (200, 400, 520, 1024, 2048 bits; IFMA, `NOIFMA=1`, `AVX2=1`):
  pass, now also for the lean mode (every 3rd tree level, no stored
  transforms), the in-place evaluation, and the in-place tree without a
  tree.
* ECM lane for lane (`cmp_s2.sh`), clean builds of all three ISAs (the
  AVX2 binary has no zmm instructions); `AVXECM_S2FFT_CFG="t,k,lean"`
  forces a configuration (t = log2 dF, every k-th tree level, lean 0/1):

| input | config | IFMA | AVX-512F | AVX2 |
|---|---|---|---|---|
| 400-bit, p40, B1 3000, B2 231284, 240 curves | model's choice | 47/47 (+1) | 47/47 (+1) | 47/47 (+1) |
| same | t=8, k=3, lean | 47/47 (+7) | 47/47 (+7) | 47/47 (+7) |
| same | t=7, k=2, full | 47/47 (+2) | 47/47 (+2) | 47/47 (+2) |
| same | t=6, k=4, lean | 47/47 (+1) | 47/47 (+1) | 47/47 (+1) |
| 1024-bit, p40, same B1/B2, 96 curves | model's choice | 23/23 (+1) | 23/23 (+1) | 23/23 (+1) |
| same | t=8, k=2, lean | 23/23 (+5) | 23/23 (+5) | 23/23 (+5) |
| 2048-bit (RNS backend), p40, 64 curves | model's choice | 7/7 | | 7/7 |
| same | t=8, k=3, lean | 7/7 (+5) | | 7/7 (+5) |
| 2^439-1, B1 1000, B2 60000, 160 curves | model's choice | | | 93/93 (+4) |

  (Forced t = 8 covers B2 up to ~6e5, hence the extra FFT-only finds.)
  At 1024 bits, B1 = 1e6, the six memory limits of §11.2 found the same
  lanes.

### 11.5 B2 from the memory limit: `-avxecm_stg2 2`

`-avxecm_stg2 2` is the FFT stage 2 with B2 chosen by `avxecm_s2fft_setup`
from the memory limit when no B2 is given (an explicit `-B2ecm` is still
used as is; modes 0 and 1 are unchanged). `s2_auto_B2` does in C what
`s2model.py` does: for every configuration (t, k, lean) that fits the limit
and every block count up to T2 = 4 T1 (1..8, then +1/8 per step), it
evaluates P(B1, covered B2) / (T1 + T2) and keeps the best.

* T1 is predicted, not measured: stage 2 is set up before stage 1 runs.
  T1 = 15 B1 template mulmods (measured 13.9-16.9 with the radix backend;
  the RNS backend's stage 1 is cheaper per mulmod, so T1 is overestimated
  there and B2 comes out somewhat large, which the flat optimum tolerates).
  T2 = 0.93 `s2_cost`.
* The target factor size comes from B1 through GMP-ECM's B1 table
  (11e3 -> p20, 5e4 -> p25, 25e4 -> p30, 1e6 -> p35, 3e6 -> p40, ...,
  interpolated in log B1); the optimum hardly depends on it (§11.3).
* P: Dickman rho table (step 1e-3 up to 40) and the one-large-prime
  integral, group order size p/23. The search takes a few ms.
* `-v` prints the chosen B2 with the target and the predicted T2/T1.

Measured (IFMA, 1 thread, B1 = 1e6):

| bits | limit | chosen B2 | dF | blocks | predicted T2/T1 | measured T2/T1 |
|---:|---:|---:|---:|---:|---:|---:|
| 1024 | 16 MB | 2.4e8 | 512 | 102 | 0.30 | 0.32 (1.61 s) |
| 1024 | 64 MB | 6.3e8 | 2048 | 16 | 0.31 | 0.30 (1.52 s) |
| 1024 | 256 MB | 1.9e9 | 8192 (lean) | 3 | 0.47 | 0.49 (2.42 s) |
| 400 | 16 MB | 4.3e8 | 2048 (lean, 3 of 11 levels) | 11 | 0.37 | 0.32 |
| 400 | 256 MB | 1.3e9 | 8192 | 2 | 0.56 | 0.52 |

For comparison, mode 1 at 1024 bits / 16 MB runs B2 = 9e8 in 5.96 s; by
§11.3 that is ~25% more time per p35 factor than the automatic choice.

Caveat: `factor()`'s curve counts and t-level accounting assume
GMP-ECM's B2 for a given B1. In mode 2 the progress output and log show
the actual B2, but the work bookkeeping is not adjusted, so mode 2 is
meant for direct `ecm()` runs; making it the default would need a
t-level calculation that uses the actual B2.
