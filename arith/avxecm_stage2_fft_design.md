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
3. **Speed.** Native FP butterflies for non-IFMA and radix-4 passes: done
   (§12). Remaining: cache-blocked / fused
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

## 12. Native floating-point NTT for non-IFMA targets; radix-4 passes (2026-11-09)

### 12.1 Design

Without IFMA, the integer NTT emulated every 52 x 52 -> 104-bit product
with double-precision FMAs (`rns_vec52.h`, ~10 vector ops per LO52/HI52;
~25 per Shoup twiddle product). The FP NTT keeps the NTT domain in doubles
instead (selected at compile time: `VP_FPNTT` = non-IFMA build; `-DVP_INT_NTT`
or `make ... INTNTT=1` forces the integer NTT for A/B tests):

* Values are doubles holding integers x with |x| <= p (p < 2^50) after
  every kernel call.
* Product modulo p (as FLINT's fft_small): h = a b (rounded),
  l = fma(a, b, -h) (exact), q = round(h / p), r = fnmadd(q, p, h) + l
  = a b - q p. Every step is exact for p < 2^50 (h - q p and r are
  integers below 2^53 in magnitude); |a|, |b| <= p gives |r| <= 0.875 p.
  Twiddles are (w, w/p) pairs with |w| <= p/2, and then
  q = round(x (w/p)) works for |x| up to 4p with |r| <= 0.75 p.
* Reduction of sums: x - p round(x/p) (3 ops).
* The pointwise product is plain a b mod p (no Montgomery 2^-52), so the
  CRT scales lose their 2^52 factor; `ctx->fpntt` (from `vpk_fp_ntt()`)
  tells `vp_ctx_init` which table format to build (same sizes).
* Conversions stay as they were (FP-emulated exact 104-bit
  multiply-accumulate) with an FP boundary: to_ntt centers its [0, 2p)
  integer residue and stores a double; from_ntt multiplies by the CRT
  scale in FP, maps to [0, p) and feeds the double straight into the
  multiply-accumulate.
* No MXCSR switch in the FP transforms (every step is exact or explicitly
  rounded); the AVX2 conversions still run with round-toward-zero, which
  the FP steps there tolerate.

Two changes for every ISA, IFMA included:

* **Radix-4 passes**: two radix-2 levels per sweep over a prime's L
  values (plus one radix-2 level when log2 L is odd), half the memory
  traffic. At L = 4096 a prime's transform is 256 KB, larger than L1, and
  the transforms were L2-bandwidth bound (AVX2 and AVX-512F FP ran at the
  same speed). In the FP path reductions sit only where the exactness
  bounds need them: DIF reduces the second level's sum (<= 4p), DIT the
  second level's outputs (<= 2p); per two levels 4 twiddle products +
  2 (DIF) or 4 (DIT) reductions instead of 4 + 4 / 4 + 8.
* **Conversion constants as doubles** (`ctx->cvtd`, `ctx->crtd`): the
  FP-emulated multiply-accumulate converted each 52-bit table entry from
  integer inside the inner loop; now it broadcasts a stored double.

### 12.2 Kernel timings

`vp_test -k -l 12` (L = 4096, best of 7, Rocket Lake, 1 thread), ns per
element-prime (transforms, pointwise) or per coefficient-prime
(conversions):

| bits | build | to_ntt | fwd | pmul | inv | from_ntt |
|---:|---|---:|---:|---:|---:|---:|
| 400 | IFMA, before radix-4 | 7.4 | 10.3 | 2.5 | 10.2 | 9.3 |
| 400 | IFMA | 7.2 | 7.7 | 2.5 | 7.8 | 9.5 |
| 400 | AVX-512F, integer NTT (before) | 16.0 | 19.3 | 4.5 | 18.3 | 17.6 |
| 400 | AVX-512F, FP NTT | 13.9 | 14.4 | 2.4 | 17.2 | 14.1 |
| 400 | AVX2, integer NTT (before) | 26.6 | 31.4 | 6.5 | 30.9 | 30.9 |
| 400 | AVX2, FP NTT | 25.7 | 14.9 | 2.8 | 21.1 | 23.5 |
| 1024 | IFMA | 14.6 | 9.1 | 4.7 | 10.4 | 19.0 |
| 1024 | AVX-512F, FP NTT | 24.2 | 15.4 | 5.1 | 19.2 | 25.2 |
| 1024 | AVX2, FP NTT | 48.3 | 16.0 | 6.1 | 24.3 | 43.4 |
| 2048 | IFMA | 29.3 | 11.2 | 6.3 | 12.7 | 35.3 |
| 2048 | AVX-512F, FP NTT | 40.8 | 16.8 | 6.1 | 20.3 | 45.6 |
| 2048 | AVX2, FP NTT | 81.5 | 17.3 | 6.7 | 24.3 | 74.7 |

* The FP transforms are 1.5-2x faster than the emulated integer ones on
  AVX2 and ~1.3x on AVX-512F; pointwise products 2x. AVX2 FP transforms
  run as fast as AVX-512F ones: Rocket Lake has one 512-bit FMA port
  (ports 0+1 fused), so native FP gets no width advantage from AVX-512.
* FP transforms are within 1.5-2x of IFMA's integer ones, and FP pointwise
  products match IFMA.
* What is left between non-IFMA and IFMA is mostly the conversions (O(K n)
  emulated 104-bit multiply-accumulates per coefficient): at 1024 bits
  they are ~60% of an AVX2 product (to_ntt 3.3x, from_ntt 2.3x IFMA).

### 12.3 ECM stage 2 (yafu, 1 thread, 8 curves, B1 = 1e6, `-avxecm_stg2 1`)

| bits | build | stage 1 | stage 2, B2 = 1.2e8 | stage 2, default B2 (9.7e8) | before §12 (default B2) |
|---:|---|---:|---:|---:|---:|
| 400 | IFMA | 1.04 s | 0.10 s | 0.38 s | 0.49 s |
| 400 | AVX-512F | 1.63 s | 0.16 s | 0.60 s | 0.71 s |
| 400 | AVX2 | 3.4 s | 0.25 s | 0.88 s | 1.13 s |
| 1024 | IFMA | 5.0 s | 0.38 s | 1.41 s | 1.55 s |
| 1024 | AVX-512F | 9.4 s | 0.65 s | 2.34 s | |
| 1024 | AVX2 | 15.0 s | 0.88 s | 3.17 s | |

Relative to IFMA, stage 2 now costs 1.6-1.7x on AVX-512F and 2.2-2.3x on
AVX2, while stage 1 costs 1.6-1.9x and 3.0-3.3x: on non-IFMA machines
stage 2 is relatively cheaper than on IFMA ones, which favors a larger B2
there (§11.3; the cost model and the automatic B2 still use constants
fit on IFMA).

### 12.4 Verification

* `vp_test` (200, 400, 520, 1024, 2048 bits, every polynomial operation
  and lane against GMP, including the all-(2N-1) bound): pass for IFMA,
  AVX-512F FP, AVX-512F integer (`INTNTT=1`), AVX2 FP and AVX2 integer.
* ECM lane for lane (`cmp_s2.sh`), clean yafu builds (AVX2 binary: no zmm;
  AVX-512F: no `vpmadd52`): 400-bit p40 (240 curves), model's choice and
  forced t=8 k=3 lean: 47/47 (+1) and 47/47 (+7); 1024-bit p40 (96
  curves), forced t=7 k=2 full: 23/23 (+1); on IFMA, AVX-512F and AVX2.
* clang-cl compiles `vp_kern.c` / `vp_setup.c` / `vp_poly.c` warning-free
  for avx2 / avx512 / avx512ifma; gcc warning-free for all five variants.

### 12.5 Next

(All three items below were taken up in §13.)

* Conversions for non-IFMA (the remaining gap): the emulated multiply-
  accumulate costs ~5 ops per 52 x 52 product; options are a native-FP
  formulation with fewer, wider partial products, or reducing K n (larger
  or fewer primes per conversion).
* Refit `s2_tcost` / the `s2_cost` ratios and S2_C1 per ISA (or measure
  the stage-1 / stage-2 cost ratio at run time) so `-avxecm_stg2 2` picks
  B2 for the machine it runs on.
* Radix-4 halved the transforms' memory traffic; a cache-blocked
  (four-step / two-pass) transform would take L = 4096..32768 out of L2
  entirely.

## 13. Conversions, cache-blocked transforms, per-ISA cost model (2026-11-09)

### 13.1 Conversions

* **from_ntt, loops interchanged** (all ISAs). The CRT multiply-accumulate
  ran prime-outer / limb-inner into n + 2 column accumulators, which do
  not fit in registers (n = 20 at 1024 bits): every multiply-accumulate
  loaded and stored two spilled accumulators. Now the K scaled residues
  y_i (and kappa as the (K+1)-th) are computed first into a stack array,
  then each column k = sum_i y_i C_i[k] is accumulated limb-outer /
  prime-inner in two registers; the high half carries into column k + 1.
  Stack: `VP_MAXK` (= 2.1 `VP_MAXW` + 4) vectors.
* **to_ntt, native-FP final reduction** (non-IFMA). The accumulated
  X = h1 2^104 + h0 2^52 + lo (X = x 2^104 mod p by the conversion table)
  is reduced as h1 + h0 (2^-52 mod p) + lo (2^-104 mod p) with two FP
  constant products and one reduction (`ctx->rdcd`), instead of the
  emulated two-word Montgomery REDC (~22 ops instead of ~45). Bound: each
  product <= 1.25 p for inputs < 2^52 even with round-toward-zero (the
  AVX2 conversions run in that mode), sum < 2.6 p, reduced to <= p/2.

ns per coefficient-prime, `vp_test -k -l 12`, before -> after:

| | from_ntt 400 / 1024 / 2048 bits | to_ntt 400 / 1024 / 2048 bits |
|---|---|---|
| IFMA | 9.5 -> 6.2 / 19.0 -> 14.2 / 35.3 -> 30.0 | (integer path, unchanged) |
| AVX-512F | 14.1 -> 9.4 / 25.2 -> 18.5 / 45.6 -> 32.7 | 13.9 -> 12.4 / 24.2 -> 21.5 / 40.8 -> 37.6 |
| AVX2 | 23.5 -> 18.6 / 43.4 -> 35.9 / 74.7 -> 64.7 | 25.7 -> 19.5 / 48.3 -> 41.3 / 81.5 -> 73.0 |

### 13.2 Cache-blocked transforms

The transforms are now recursive and depth first (`vp_fwd_rec`,
`vp_inv_rec`): above the block size (`ctx->blk`, default `VP_BLK` = 512
elements = 32 KB per prime) one radix-4 pass (or a radix-2 pass to even
out the level count), then the sub-blocks; a block of `blk` elements
then finishes all its remaining levels while it is in L1. Levels up to
`VP_TWS_M` = 256 read their twiddles from a contiguous per-level table
(`ctx->tws` / `itws`, K (2 * 256 - 1) * 16 bytes each, shared) instead of
striding through the w_Lmax table, so the in-L1 levels' twiddles stay in
L1 too. The pass functions are per representation (FP / integer); the
recursion is shared. `vp_test -B blk` sets the block size (tested down
to 4 and 8, which exercise every recursion path).

* One thread: no measurable difference, even at L = 16384 (1 MB per
  prime, twice L2): Rocket Lake's prefetchers stream from L3 well enough.
* 8 concurrent copies (as yafu runs ECM, one batch per thread), L = 16384,
  ns per element-prime, blocked / unblocked: IFMA fwd 82 / 98, inv
  65 / 106; AVX2 fwd 67 / 93, inv 66 / 117 (20-45% faster). L = 4096:
  neutral. Kept on, B = 512.
* Note that 8 concurrent copies run the transforms ~4x slower per copy
  than one copy (shared L3 / memory bandwidth); stage 2 throughput with
  many threads is bandwidth-limited, which favors smaller dF per thread
  more than the single-thread model says.

### 13.3 Per-ISA cost model

`s2_cost` (§11.2) is now parameterized per build (`s2_model_t s2m` in
`avxecm_stage2fft.c`, selected by `USE_AVX512IFMA` / `USE_AVX512F` /
AVX2): tree cost tc = a t + b + c t^2/n and the per-tree ratios fit by
`ecm_tests/fit_s2.py` to `vp_test -b` per ISA (n = 4..40, t = 8..13; max
error 14-28%), the stage-1 cost c1 (mulmods per unit of B1, from yafu
stage-1 times) and a scale fit to yafu stage-2 times (measured / model
0.84-1.09 over 5 configurations per ISA):

| ISA | a, b, c | mm, mi, ev (full) | lean | c1 | scale | template mulmod ns, n = 8 / 20 |
|---|---|---|---|---|---|---|
| IFMA | 6.82, -37.0, 3.34 | 0.82, 0.79, 1.84 | 0.94, 0.57, 1.91 | 15.5 | 0.90 | 61 / 359 |
| AVX-512F | 3.45, -7.6, 2.69 | 0.72, 0.76, 1.92 | 0.84, 0.55, 1.87 | 14.2 | 0.92 | 105 / 732 |
| AVX2 | 4.22, -14.7, 1.40 | 0.67, 0.68, 1.83 | 0.73, 0.49, 1.81 | 14.5 | 1.00 | 223 / 1032 |

In units of each build's own template mulmod, the polynomial work is 2-3x
cheaper on AVX-512F / AVX2 than on IFMA (their template mulmod is
relatively slower than their native-FP transforms), so stage 2 is
relatively cheaper there and the best B2 is larger. `s2model.py --isa
ifma|f|avx2` uses the same constants; e.g. 1024 bits, B1 = 1e6, 256 MB:
best B2 2.6e9 (IFMA) / 2.6e9 (AVX-512F) / 3.2e9 (AVX2); 400 bits, 64 MB:
9.7e8 / 1.1e9 / 2.6e9. GMP-ECM's default B2 with ample memory loses up to
~5% on IFMA and up to ~12% on AVX2 against the best.

`-avxecm_stg2 2` (B1 = 1e6, 1 thread), predicted vs measured T2/T1:

| bits | limit | IFMA | AVX-512F | AVX2 |
|---:|---:|---|---|---|
| 400 | 256 MB | B2 1.9e9: 0.55 / 0.55 | B2 1.9e9: 0.47 / 0.52 | B2 2.6e9: 0.43 / 0.43 |
| 400 | 64 MB | 9.7e8: 0.38 / 0.36 | 1.1e9: 0.37 / 0.40 | 2.6e9: 0.49 / 0.50 |
| 1024 | 256 MB | 2.6e9: 0.45 / 0.49 | 2.6e9: 0.37 / 0.39 | 3.2e9: 0.39 / 0.40 |
| 1024 | 64 MB | 7.8e8: 0.31 / 0.35 | 8.6e8: 0.29 / 0.30 | 9.4e8: 0.29 / 0.29 |

### 13.4 ECM stage 2 now (yafu, 1 thread, 8 curves, B1 = 1e6, default B2 9.7e8, 256 MB)

| bits | ISA | stage 1 | stage 2 | stage 2 / IFMA | before §13 | stage 1 / IFMA |
|---:|---|---:|---:|---:|---:|---:|
| 400 | IFMA | 1.02 s | 0.35 s | 1 | 0.38 s | 1 |
| 400 | AVX-512F | 1.65 s | 0.57 s | 1.63 | 0.60 s | 1.6 |
| 400 | AVX2 | 3.27 s | 0.81 s | 2.31 | 0.88 s | 3.2 |
| 1024 | IFMA | 5.05 s | 1.39 s | 1 | 1.41 s | 1 |
| 1024 | AVX-512F | 9.36 s | 2.06 s | 1.48 | 2.34 s | 1.85 |
| 1024 | AVX2 | 14.73 s | 2.95 s | 2.12 | 3.17 s | 2.9 |

(For reference, the pairing stage 2 at B2 = 1e8: IFMA 0.99 / 2.96 s, AVX2
1.93 / 6.07 s at 400 / 1024 bits; the FFT stage 2 at B2 = 1.2e8: 0.11 /
0.39 s and 0.23 / 0.90 s.)

### 13.5 Verification

* `vp_test` 200..2048 bits, every operation and lane vs GMP: pass for
  IFMA, AVX-512F FP / integer, AVX2 FP / integer, with the default block
  and with blocks of 4 and 8 elements.
* ECM lane for lane (`cmp_s2.sh`), clean builds of the three ISAs (AVX2:
  no zmm; AVX-512F: no `vpmadd52`), every row 0 pairing finds missed:
  400-bit 47/47 (+1), 1024-bit 23/23 (+1), 2048-bit (RNS) 7/7, forced
  t=8 k=3 lean 47/47 (+7) and 7/7 (+5), forced t=7 k=2 full 23/23 (+1);
  AVX2 2^439-1 93/93 (+4).
* gcc warning-free for the five standalone variants; clang-cl compiles
  `vp_setup.c`, `vp_kern.c`, `vp_poly.c`, `avxecm_stage2fft.c` for avx2 /
  avx512 / avx512ifma without warnings.

### 13.6 Next

* The conversions are now the largest share on non-IFMA (to_ntt +
  from_ntt ~60% of an AVX2 product at 1024 bits): their emulated 104-bit
  multiply-accumulate is ~5 ops per 52 x 52 product, close to the floor
  for exact FP products; the remaining levers are fewer products (larger
  digits / fewer primes) or IFMA.
* Multithreaded stage 2 is bandwidth-bound (§13.2); a per-thread memory
  model that accounts for shared L3 / DRAM bandwidth, or sharing one
  batch's polynomial work across threads, would help large thread counts.

## 14. Brent-Suyama extension (2026-11-09)

### 14.1 Method

With a root polynomial f, the FFT stage 2 uses the roots x([f(k)]Q)
(baby, k in S1) and x([f(jD)]Q) (giant). The product of the differences
vanishes mod p when f(jD) = +-f(k) (mod ord Q). For f(x) = x^S or the
Dickson polynomial D_S(x, -1) (D_0 = 2, D_1 = x, D_n = x D_{n-1} + D_{n-2}),
f(jD) - f(k) contains (jD - k)(jD + k), so every prime the plain stage 2
finds is still found (q = jD +- k); a larger q is found when it divides
one of the other algebraic factors, i.e. when f(jD) = +-f(k) mod q for a
non-trivial pair. S = 1 is the existing stage 2 and stays the default.

* **Points with y.** Finite differences need general point additions,
  which x-only Montgomery arithmetic does not have. Q = (x0 : 1) is put on
  the curve g y^2 = x^3 + A x^2 + x with g = x0^3 + A x0^2 + x0 and
  y = 1: the stage-1 curve or its quadratic twist, either way with the
  same x-coordinates for every multiple of Q (which is all stage 2 uses).
  Affine addition: lambda = (y2 - y1)/(x2 - x1), x3 = g lambda^2 - A - x1
  - x2, y3 = lambda (x1 - x3) - y1; doubling lambda = (3x^2 + 2Ax + 1)/(2gy).
* **Tables.** For an arithmetic progression i_m = i0 + m delta the
  finite-difference table [Delta^e f(i_0)]Q, e = 0..deg, advances by
  E_e <- E_e + E_{e+1} (deg additions per root; E_0 is the root). Initial
  entries: x-only ladder for [c]Q and [c+1]Q (c = Delta^e f(i0), an
  integer of ~deg log2(i) bits), then Okeya-Sakurai y-recovery
  y_c = [Z2((x0 X1 + Z1)(x0 Z1 + X1 + 2A Z1) - 2A Z1^2) - (x0 Z1 - X1)^2 X2]
  / (2g Z1^2 Z2) in projective form, all normalized with one batch
  inversion.
* **M progressions side by side** so each step's M deg denominators share
  one batch inversion (one `mpz_invert` per lane): giant j = j0 + r + M m
  (r < M, i0 = (j0 + r) D, delta = M D); baby k = c + 6 (r + (M/2) m),
  c in {1, 5} (S1 consists of such k because 6 | D), kept when k < D/2 and
  gcd(k, D) = 1. M (a power of 2, `s2_bs_m`) balances table setup
  (M (deg + 1) ladders) against inversions (one per M roots): ~32 for
  degree 3, 16 for 6, 8 for 12, 2 for 30 at dF = 4096.
* **Integer coincidences.** The table multipliers are positive integers
  (positive coefficients and arguments) and are tracked exactly (mpz).
  When two adjacent ones are equal the two points are equal and the step
  is a doubling (found by the oracle: x^2 with j0 + r = M/2 makes
  Delta f = Delta^2 f exactly, which collapsed every lane before the fix).
  Coincidences modulo ord Q (T_e = +-T_{e+1} or T_e = 0 mod q) make a
  denominator vanish mod p: the batch inversion folds it into the lane's
  result, a genuine extra find (seen with D_12 / D_30 at tiny q; verified
  case by case).

### 14.2 Code and options

* `-avxecm_bs S`: 0 or 1 off, `auto` (§14.9; the default with
  `-avxecm_stg2 2`, else off), `S` or `pS` x^S, `dS` Dickson
  D_S(x, -1), S <= 32 (a string, since the option parser takes `-` as the
  next switch). `ecm_obj.avx_stg2_bs` (> 1 power, < -1 Dickson);
  `avxecm_s2fft_setup(..., bs, ...)`.
* `avxecm_stage2fft.c`: `s2_fpoly`, `s2_findiff`, `s2_inv_batch`,
  `s2_mul2` (ladder with an mpz scalar returning [c]Q and [c+1]Q),
  `s2_bs_init` (tables, y-recovery), `s2_bs_step`, `s2_bs_baby`; the
  giant loop emits M roots per step. (First version; §14.7 describes the
  current tables, setup and cost terms.)
* Test aid: `AVXECM_S2FFT_BSM=M` forces the number of progressions.
* `-v` prints the polynomial, degree and M.
* Cost model (`s2_cost`): per root deg (7 + 12/n) + 1 + 8 (4 + 160/n)/M
  mulmods, ~1.6 dF baby candidates, table setup M (deg + 1) 10 deg
  (log2(giant arg) + log2(D/2)); measured / model 0.88-1.15 on all three
  ISAs at 400 and 1024 bits (1.23 for D_30 at 400 bits on IFMA).

### 14.3 Verification: an exact oracle

`ecm_tests/bs_oracle.py` decides exactly which curves stage 2 must find
for a small p: Suyama curve from sigma (as yafu builds it), order of the
starting point by baby-step giant-step on the Hasse interval, yafu's
stage-1 exponent (2^e and p^e < B1 for p < B1), then the condition
f(jD) = +-f(k) (mod ord Q) over the actual D, dF, blocks; stage-1 lanes
also count as stage-2 finds (Q is the identity mod p). Self-test: the
x-only ladder against the affine group law on the twist and the
y-recovery formula. `ecm_tests/oracle_cmp.sh` compares it with yafu lane
by lane.

* Plain stage 2: 0 differing batches of 42 (336 curves, 400 and 1024 bits)
  - an exact check of the existing FFT stage 2 too.
* Brent-Suyama, 400-bit p40 (240 curves, B1 3000, B2 231284) and 1024-bit
  (96 curves): x^2, x^3, x^6, x^12, D_2, D_3, D_6, D_12, D_30, M = 4 / 16 /
  adaptive, on IFMA, AVX-512F and AVX2: every lane the oracle predicts is
  found; the only differences are 1-2 extra yafu finds for D_12 / D_30,
  each verified to be a table collision mod ord Q (e.g. q = 20635493:
  T_22 = +-T_23 in a giant table; q = 552833: T_9 = +-T_10 in a baby
  table) - genuine factors.

### 14.4 Cost (first version; §14.7 has the current one)

Stage 2 time relative to the plain stage 2 (yafu, 1 thread, B1 1e6, B2 9.7e8,
dF 4096, 6 blocks):

| ISA, bits | D_3 | D_6 | x^6 | D_12 | D_30 |
|---|---:|---:|---:|---:|---:|
| IFMA 400 / 1024 | 1.17 / 1.10 | 1.41 / 1.34 | 1.55 / 1.35 | 1.87 / 1.88 | 3.55 / 3.49 |
| AVX-512F 400 / 1024 | 1.14 / 1.20 | 1.35 / 1.56 | 1.37 / 1.52 | 1.77 / 2.14 | 3.24 / 4.05 |
| AVX2 400 / 1024 | 1.20 / 1.21 | 1.49 / 1.54 | 1.49 / 1.58 | 2.16 / 2.29 | 3.92 / 4.10 (1) |

(1) AVX2 1024-bit D_30 picked dF = 8192, 2 blocks. The extension's work
is mostly template mulmods (7 per affine addition), so relative to the
NTT-heavy plain stage 2 it costs more on AVX-512F / AVX2.

### 14.5 Probability

Expected number of non-trivial pairs for a prime q > B2:
E(q) = nbaby ngiant r / q, r the number of non-trivial x with
f(x) = +-f(y) mod q. r depends on q's residue class (x^6: gcd(6, q - 1)
etc.), so the chance 1 - exp(-E) is averaged over a sample of r values
(`bs_model.py`, 1200 samples over random q): means x^2 1.00, D_3 1.94,
x^6 3.77, D_6 3.78, D_12 5.66, D_30 9.58 (powers and Dickson polynomials of
the same degree have nearly the same mean; Dickson spreads it over all q).
Oracle-exact check on 3000 random curves (p 40 bits, B1 3000, D 510,
dF 64, 7 blocks; 1893 curves with one prime q > B2 left): caught by
D_6 49 (model 47.4), D_30 95 (87.6), x^12 57 (62.8) - within Poisson
noise. Stage-2 finds grow +11% (D_6) / +21% (D_30) at this tiny size.

### 14.6 Is it worth it? (first version; §14.8 has the current answer)

Least expected time per factor, (T1 + T2)/P, each option with its own best
B2 and configuration within the memory limit, relative to the plain
stage 2 (`bs_model.py table`, `--isa`; per-ISA cost constants of §13.3 and
the Brent-Suyama terms of §14.2):

| bits | B1 | p | MB | IFMA best | AVX-512F best | AVX2 best | D_12 (IFMA / AVX2) | D_30 (IFMA / AVX2) |
|---:|---:|---:|---:|---|---|---|---|---|
| 400 | 2.5e5 | 30 | 256 | plain | plain | plain | 1.14 / 1.26 | 1.39 / 1.63 |
| 400 | 1e6 | 35 | 256 | D_3 0.997 | plain | plain | 1.06 / 1.13 | 1.23 / 1.41 |
| 400 | 3e6 | 40 | 64 | D_6 0.981 | D_3 0.990 | D_3 1.000 | 1.01 / 1.07 | 1.10 / 1.20 |
| 1024 | 1e6 | 35 | 256 | D_3 0.994 | plain | plain | 1.07 / 1.13 | 1.24 / 1.37 |
| 1024 | 3e6 | 40 | 256 | D_3 0.987 | D_3 0.994 | D_3 0.998 | 1.02 / 1.07 | 1.12 / 1.21 |
| 1024 | 1.1e7 | 45 | 256 | D_6 0.978 | D_3 0.989 | D_3 0.991 | 1.00 / 1.04 | 1.07 / 1.14 |
| 1024 | 4.3e7 | 50 | 1024 | D_6 0.973 | | D_3 0.989 | 0.99 / 1.02 | 1.03 / 1.10 |
| 1024 | 1.1e8 | 55 | 1024 | D_6 0.971 | | D_3/D_6 0.988 | 0.98 / 1.01 | 1.02 / 1.08 |
| 2048 | 4.3e7 | 50 | 1024 | D_6 0.977 | | D_3 0.989 | 1.00 / 1.03 | 1.05 / 1.11 |

Conclusions:

1. **Small gain at best.** Brent-Suyama with D_3 / D_6 lowers the expected
   time per factor by 1-3% on IFMA from p40 up (2.5-3% at p50-p55), ~1% on
   AVX-512F / AVX2; at p30-p35 it is break-even or a loss. Degree 12 and
   30 lose almost everywhere (they only approach break-even at p55).
2. **Why so little.** At these dF the extra chance comes from q a few
   times above B2 (E(q) ~ pairs r / q with pairs ~ dF^2 blocks), while
   each root costs 7 deg template mulmods on top of ~150-300 for the
   polynomial work; GMP-ECM's larger dF (2^16 and up) puts far more pairs
   behind each root. The best B2 with the extension is somewhat smaller
   than without it.
3. **Non-IFMA gains least**: the root work is template mulmods, which are
   relatively expensive there (§13.3).
4. **Recommendation:** keep it off by default; `-avxecm_bs d6` on IFMA
   (or `d3` elsewhere) is a ~2-3% improvement for p45+ work. The
   probability model is validated only at 40-bit p (exact oracle); at
   p45+ it is the standard rho-based estimate.

### 14.7 Cheaper roots (2026-11-09, second pass)

Â§14.4 - 14.6 describe the first version. The second pass makes each root
cheaper without changing which roots are produced:

* **Scaled Weierstrass coordinates.** With X = g x, Y = g^2 y the curve
  g y^2 = x^3 + A x^2 + x becomes Y^2 = X^3 + a2 X^2 + a4 X (a2 = g A,
  a4 = g^2): the addition x3 = g lambda^2 - ... loses its multiplication by
  g, 6 mulmods per affine addition (3 batch inversion, lambda, lambda^2,
  lambda (X1 - X3)) instead of 7. All roots (baby and giant) are scaled by
  the same g per lane, so the product of differences only gains a factor
  g^(dF ...), invertible unless Q has order 2 mod p (a find anyway).
* **Entry-major tables, vector calls.** Entry (e, r) lives at e M + r. A
  step reads entries e < deg and e + 1 > 0, two contiguous ranges offset by
  M, so each operation of the step is one call over deg M coefficients
  (`s2_aff_add`) instead of deg M single-coefficient calls; the M roots of
  a step are one contiguous copy. The batch inversion (`s2_inv_batch`)
  runs S2_NCH = 4 interleaved chains so its products are vector calls too,
  at 3 S2_NCH extra mulmods per inversion.
* **2-D table setup.** For one family of C progressions i0 + r' sigma
  (step delta along each), row e of the tables over r' is
  g_e(i0 + r' sigma) with g_e = Delta_delta^e f a polynomial of degree
  deg - e, i.e. itself a finite-difference table in steps of sigma. The
  setup ladders only the (deg + 1)(deg + 2)/2 entries of that triangle
  (Delta_sigma^j Delta_delta^e f(i0), j + e <= deg, all positive) and then
  steps the triangle C - 1 times (deg (deg + 1)/2 additions per step, one
  batch inversion), copying its column j = 0 into progression r'. The
  1-D setup (C (deg + 1) ladders) is kept for C <= ~(deg + 2)/2;
  `s2_bs_famcost` picks the cheaper. Giant tables are one family
  (i0 = j0 D, sigma = D, delta = M D); baby tables two (k = 1 + 6 r' and
  5 + 6 r', delta = 3 M).
* **More progressions.** With the setup nearly independent of M, M now
  balances one batch inversion (8 lanes x mpz_invert) per M roots against
  the C - 1 extra triangle steps: `s2_bs_m` evaluates the cost for M = 2,
  4, ..., 128 (<= dF) and takes the least. Typical: 128 for D_3 / D_6 at
  400 bits, 64 at 1024 bits (D_6) and for D_12, 2-4 for D_30 (still 1-D).
* **Multipliers mod 2^61 - 1.** The integer multipliers that detect equal
  points (doubling, Â§14.1) are tracked modulo the prime 2^61 - 1 instead of
  as mpz. Not mod 2^64: the high differences (deg! delta^deg and their
  neighbours) are divisible by large powers of 2, so mod 2^64 many unequal
  multipliers look equal (found by the oracle: D_20 / D_30 lost lanes).
* Scratch per thread: 2 (deg + 1) M + 7 max((deg + 1) M,
  (deg + 1)(deg + 2)/2) + 2 (deg + 1)(deg + 2)/2 + 36 coefficients
  (D_6, M = 64, 1024 bits: ~5.5 MB). `-v -v` prints the setup time.

Cost model (`s2_cost`, mirrored in `bs_model.py`): per root
deg (6 + 12/n) + 1 + (8 (4 + 160/n) + 12) / M template mulmods; setup per
family as above with a ladder at 10 deg log2(argument) + 20 mulmods.
Measured / model stage-2 ratio 0.93 - 1.04 for D_3 / D_6 / x^6 on all three
ISAs at 400 and 1024 bits, 0.92 - 1.16 for D_12 / D_30. Measured per
affine addition ~7 template mulmods (D_6, IFMA, 1024 bits); table setups
0.013 s (D_3), 0.038 s (D_6), 0.21 s (D_12) of a 1.4 - 2.4 s stage 2.

Verification: `oracle_cmp.sh` on IFMA, AVX-512F and AVX2 for plain, x^2,
D_3, x^3, x^6, D_6, D_12, D_20, D_30, adaptive M and M = 4, 400 and
1024 bits: every predicted lane found; the only differences are 1-2 extra
D_12 / D_30 finds per 30 batches, each a table collision mod ord Q
(`ecm_tests/coll2.py` simulates the new 1-D / 2-D tables and lists them;
e.g. ord Q = 2871359, giant table step 50; ord Q = 676057, baby table).

Stage 2 time relative to the plain stage 2 (B1 1e6, B2 9.7e8, dF 4096,
6 blocks, 1 thread; first version in parentheses):

| ISA, bits | D_3 | D_6 | x^6 | D_12 | D_30 |
|---|---:|---:|---:|---:|---:|
| IFMA 400 | 1.07 (1.17) | 1.25 (1.41) | 1.34 (1.55) | 1.72 (1.87) | 2.92 (3.55) |
| IFMA 1024 | 1.02 (1.10) | 1.23 (1.34) | 1.23 (1.35) | 1.73 (1.88) | 3.16 (3.49) |
| AVX-512F 400 | 1.13 (1.14) | 1.23 (1.35) | 1.32 (1.37) | 1.70 (1.77) | 2.75 (3.24) |
| AVX-512F 1024 | 1.12 (1.20) | 1.37 (1.56) | 1.38 (1.52) | 2.29 (2.14) (1) | 3.70 (4.05) |
| AVX2 400 | 1.11 (1.20) | 1.32 (1.49) | 1.32 (1.49) | 1.87 (2.16) | 3.36 (3.92) |
| AVX2 1024 | 1.13 (1.21) | 1.41 (1.54) | 1.40 (1.58) | 2.15 (2.29) | (2) |

(1) a noisy run (a rerun gave giant roots 2.6 s instead of 3.0 s, i.e.
~2.0). (2) picked dF = 8192 with 2 blocks, not comparable.

### 14.8 Is it worth it now?

`bs_model.py table --isa ...` and the larger cases, time per factor
relative to the plain stage 2, each with its own best B2 and configuration
(first version, Â§14.6, in parentheses):

| bits | B1 | p | MB | IFMA best | AVX-512F best | AVX2 best | D_12 IFMA / AVX2 |
|---:|---:|---:|---:|---|---|---|---|
| 400 | 2.5e5 | 30 | 256 | D_3 0.990 (plain) | D_3 0.999 (plain) | plain (plain) | 1.12 / 1.23 |
| 400 | 1e6 | 35 | 256 | D_3/D_6 0.984 (0.997) | D_3 0.988 (plain) | D_3 0.996 (plain) | 1.03 / 1.10 |
| 400 | 3e6 | 40 | 64 | D_6 0.965 (0.981) | D_6 0.975 (0.990) | D_3 0.987 (1.000) | 0.99 / 1.03 |
| 1024 | 1e6 | 35 | 256 | D_3 0.981 (0.994) | D_3 0.988 (plain) | D_3 0.991 (plain) | 1.04 / 1.09 |
| 1024 | 3e6 | 40 | 256 | D_6 0.973 (0.987) | D_3 0.983 (0.994) | D_3 0.985 (0.998) | 1.00 / 1.04 |
| 1024 | 1.1e7 | 45 | 256 | D_6 0.968 (0.978) | D_6 0.980 (0.989) | D_3 0.983 (0.991) | 0.99 / 1.02 |
| 1024 | 4.3e7 | 50 | 1024 | D_6 0.966 (0.973) | D_6 0.976 | D_6 0.982 (0.989) | 0.98 / 1.01 |
| 1024 | 1.1e8 | 55 | 1024 | D_6 0.965 (0.971) | D_6 0.975 | D_6 0.980 (0.988) | 0.97 / 1.00 |
| 2048 | 4.3e7 | 50 | 1024 | D_6 0.969 (0.977) | D_6 0.981 | D_3 0.983 (0.989) | 0.98 / 1.01 |

* The gain roughly doubled: IFMA 2-3.5% (D_6 from p40 up, D_3 below),
  AVX-512F 1.5-2.5%, AVX2 1-2%. Brent-Suyama now wins almost everywhere;
  plain is best only for small B1 on AVX2.
* D_3 and D_6 are within ~0.5% of each other almost everywhere (D_6 ahead
  for larger p and on IFMA). D_12 approaches break-even only for p50+ on
  IFMA and never wins; D_30 never wins; the powers x^S at best tie with
  the Dickson polynomial of the same degree (same cost, nearly the same
  r mean).
* The remaining limit is the ~6 deg mulmods per root (the polynomial work
  per root is ~150-300), not setup or inversions.

### 14.9 Automatic choice

`-avxecm_bs auto` chooses among plain, D_3, D_6, D_12 the least predicted
time per factor; it is the default with `-avxecm_stg2 2` (automatic B2),
where polynomial and B2 are chosen together. An explicit `-avxecm_bs`
(including `0` for plain) overrides it; with `-avxecm_stg2 1` the default
stays plain (`auto` must be asked for; B2 is then fixed and only the
polynomial is chosen, `s2_auto_bs`).

* `s2_prob(..., pairs, bs)` adds the q > B2 term of Â§14.5 with the
  r distributions measured by `bs_model.py rsamples` (1200 samples),
  stored as `s2_rhist` for x^2, x^3, x^6, D_3, D_6, D_12 (r: weight):
  D_3 {4: 0.485}, D_6 {2: 0.251, 4: 0.262, 10: 0.223}, D_12 {2: 0.147,
  4: 0.262, 6: 0.104, 10: 0.102, 22: 0.122} (r = 0 otherwise).
  pairs = nbaby x blocks x dF. With an explicit polynomial the automatic
  B2 now also counts its extra chance.
* `s2_auto_B2` loops over the candidates (`s2_bs_cand`) as well as t, k,
  lean and the block count; `-v` prints the choice
  ("... (chosen automatically)" / "no Brent-Suyama (chosen
  automatically)").
* Option values: `AVXECM_BS_AUTO` = 64 (`auto`), `AVXECM_BS_DEFAULT` = 65
  (option not given) in `avxecm_stage2fft.h`; `avx_ecm_main.c` maps the
  default to auto or off.
* Choices on Rocket Lake (1 thread, 256 MB): IFMA D_3 at B1 2.5e5 (400
  bits) and 1e6 (1024), D_6 from B1 3e6 / 1.1e7 up; AVX-512F the same;
  AVX2 plain at B1 2.5e5 (400 bits), D_3 elsewhere up to 4.3e7 - matching
  the tables above.

### 14.10 Next

* Validate the time-per-factor gain end to end (curves to find a known
  p40-p45 with plain vs auto) - expensive: ~1000s of curves per point for
  a 2-3% effect.
* The per-root cost is now the affine additions themselves; fewer
  mulmods per root would need a different root generation (e.g. GMP-ECM's
  polynomial-evaluation roots for larger dF).
