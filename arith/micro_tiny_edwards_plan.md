# Edwards curves and AVX2 for microecm and tinyecm: review and plan

Status (2026-11-09, branch `wip_edwards`): Phase 0 done, the gate passes
for microecm (scalar and x8) and, more narrowly, for tinyecm scalar (Part
C). No production code yet.

Goal: new, standalone implementations of microecm (n <= 64 bits) and
tinyecm (n <= 128 bits) that run stage 1 on twisted Edwards curves:
* single lane;
* AVX-512 vectorized (x8, IFMA and AVX-512F);
* AVX2 vectorized (x4, FP52 arithmetic from `arith/avx2mont`).

The existing `microecm.c` and `tinyecm.c` stay as they are. The new code
gets its own files and entry points, so the two can be compared head to
head and the callers switched later.

Part A is the review: what the existing code does, and why the AVX-ECM
results (`avxecm_edwards_plan.md`) do not carry over one to one. Part B is
the plan.

---

# Part A: review

## A1. The existing code

| | microecm scalar | microecm x8 | tinyecm scalar | tinyecm x8 |
|---|---|---|---|---|
| n | <= 64 bits | < 2^52 | <= 128 bits | <= 104 bits |
| residue | 1 x 64 | 1 x 52, 8 lanes | 2 x 64 (`monty128`) | 2 x 52, 8 lanes |
| mulredc | `mulx` + asm (gcc) | IFMA, or FP52 via `fmadd_round` | `mulmod128` | always FP52, even with IFMA |
| sqr | = mul | = mul | `sqrmod128` | = mul |
| add/sub | exact | exact | exact | exact |
| curve | Suyama sigma (first 8 prebuilt) | Suyama, random sigma | Suyama | Suyama |
| curve build | ~15M + 1 `uecm_modinv_64` | same, per-lane modinv | ~20M + `mpz_invert` + mpz conversions | per-lane `mpz_invert` |
| stage 1 | PRAC per prime, deep-searched step tables (`uprac70`, `uprac85`) | same | same | same |
| stage 2 | x-only pairing, D = 60, B2 = 25 B1 (50 B1 at B1 >= 125) | same | same, D = 120 tables | same |
| entry | `getfactor_uecm` | `getfactor_uecm_x8_list` | `getfactor_tecm` | `getfactor_tecm_x8`, `_x8_list` |

* B1 levels: 27, 47, 70, 85, 125, 165, 205, with 32-48 curves, chosen by
  the input's bit size (`uecm_dispatch`, `tecm_dispatch`, `*_x8_list`).
  Optional "arbitrary" pre-passes at B1 = 47 / 70 / 125.
* Stage 1 also doubles until 2^k >= 4 B1 (Jeff's ~1% tweak) and merges
  some primes into composites with cheaper chains (37 x 53, 61 x 127, ...).
* Callers: `getfactor_uecm` and `getfactor_tecm(_x8)` from
  `factor/batch_factor.c`, `factor/autofactor.c`, `factor/qs/tdiv.c`,
  `top/cmdParser/calc.c`. The `_x8_list` entry points are used only by
  `top/test.c` benchmarks today. `test/layer2/test_ecm.c` has miss-rate
  tests. The lasieve5 directories carry their own copies (out of scope).
* Side observations (they apply to the new code; the old files are not
  touched):
  - microecm x8 computes m = T_lo * invN mod 2^52 with
    `_mm512_mullo_epi64` (3 uops, ~15 cycles latency) plus a mask. With
    IFMA, `madd52lo(T_lo, invN)` gives the same 52 bits in one 4-cycle uop.
  - tinyecm x8 never uses IFMA for its 104-bit products.
  - tinyecm scalar builds every curve through mpz.

## A2. Why the AVX-ECM result does not transfer directly

AVX-ECM got 1.16-1.29x in stage 1 from three things: cheap squarings
(S ~ 0.78 M), long scalars (64k-bit chunks, w = 8-10), and +6.5% yield
from family (iv). Here:

**1. The scalars are short.** log2 s per level, as the code builds it:

| B1 | 27 | 47 | 70 | 85 | 125 | 165 | 205 |
|---|---:|---:|---:|---:|---:|---:|---:|
| bits of s | 38 | 66 | 99 | 126 | 176 | 236 | 300 |

Windows barely amortize: w = 3-5, and the table is a visible share of the
cost. This is the regime of Bos-Kleinjung chains ([BK]: 9% fewer M+S
than windows at B1 = 256). Since B1 is fixed per level, the chain for
each level can be searched offline and stored as a digit table, as
`uprac70` / `uprac85` already do for PRAC.

**2. With one word, a squaring costs a multiplication.** Edwards' main
saving (4 of the doubling's 7 products are squarings) is gone in
microecm, both scalar and x8. Counting M only:

| | M per bit |
|---|---:|
| PRAC, x-only (1.25-1.37 adds/bit x 6M + ~0.15 dup x 5M) | ~8.2-9.0 |
| Edwards a = -1, w = 4, 126 bits (dbl 7M, +1M T before an add, add 7M, table 3 adds) | ~8.8 |
| Edwards, offline-optimized chain (BK-style, -5..9%) | ~8.1-8.4 |

So the operation count is roughly a tie, and the 2^k doublings get worse
(Montgomery dup 5M, Edwards dbl 7M). For tinyecm's 2-word arithmetic a
squaring saves one of four partial products (S ~ 0.8 M), which helps a
little.

**3. These kernels are latency-bound, and Edwards has more ILP.** This is
the argument for doing it anyway, and the first thing to measure:
* One curve (scalar) or one vector of curves (x8) is a single dependency
  chain. The scalar `mulredc` is ~9 cycles latency (code comment,
  Skylake) against a throughput of ~2-3 cycles. The x8 IFMA `mulredc` is
  madd52 (4) -> vpmullq (15) -> madd52hi (4) -> compare/blend: ~25+
  cycles latency against ~4 cycles throughput.
* PRAC's differential add is three multiplications deep (2 muls ->
  2 sqrs -> 2 muls) and the chain is serial, ~4.5 mul latencies per bit.
* The Edwards doubling is two deep (4 independent squarings, then 3-4
  independent products), ~2 latencies per bit plus a little for adds.
* If the kernels really are latency-bound, Edwards could be 1.3-2x
  faster in stage 1 even at equal M counts.
* The fair control is PRAC with two interleaved curves (or two vectors),
  which also fills the pipeline. That doubles the batch size, though,
  which the scalar `getfactor_uecm` (one n, curves until a find) can't
  use.

**4. Curve construction is a real cost here.**
* Family (iv) built as in AVX-ECM ([m]Q on E', one inversion per step)
  would cost more than a whole B1 = 27 curve.
* Even today's single inversion per curve is ~15-25% of a small-B1
  microecm curve, and much more in tinyecm scalar (mpz).
* Fix: a precomputed table of small-height curves, loaded in projective
  form with no inversion (x1 = xn/xd, y1 = yn/yd gives X = xn yd,
  Y = yn xd, Z = xd yd, T = xn yn; d = dn/dd is folded into the cached
  table entries by scaling them with dd), and a projective A24 = num/den
  in stage 2's few doublings (+2M each).
* That makes **zero inversions per curve**.

**5. Yield.** Family (iv) found +5.4% (32-bit p, B1 = 400) to +9.7%
(48-bit) more than Suyama in Part C1. At 16-32-bit factors and B1 =
27-205 this must be re-measured. Because additions are only ~20% of the
operations at these scalars, a = +1 curves become candidates too: EECM's
Z/2 x Z/8 and Z/12 families have small-height members, and their 9M
addition (vs 8M) costs only ~2%. Higher torsion matters more at small B1.

**6. Stage 2 stays.** The hand-off of AVX-ECM §A7 (X_M = Z + Y,
Z_M = Z - Y, A24 = 1/(1+d), or 1/(1-d) for a = 1) works unchanged, so the
existing x-only pairing stage 2 is reused. Stage 2 is roughly 35-50% of a
curve here (B2 = 25-50 B1), so a stage-1 speedup is diluted about by
half in the total.

**Realistic expectations:**
* Time per factor 0-10% better in the scalar code if the latency argument
  is weak, 10-30% if it holds.
* The vector versions probably gain the most (longest latencies).
* AVX2 x4 is new capability. Its value is against scalar `mulx` code,
  which is strong at 64 bits.

Phase 0 decides with measurements before any production code is written.

---

# Part B: plan

## B0. Design

* **New files, old ones untouched:**

  | file | content |
  |---|---|
  | `factor/gmp-ecm/edw_small_kern.h` | curve code, written once against an arithmetic macro layer (the `rns_yaob_kern.h` / `rns_vec52.h` pattern): table curve load, Edwards dbl/add/cache, chain walker, stage-1 gcd, Montgomery hand-off, x-only stage 2 (copied from the old files, with projective A24), and a PRAC stage 1 as a control |
  | `factor/gmp-ecm/edw_small_tables.h` | generated: curve table(s), per-B1 digit chains, stage-2 pair tables (copied) |
  | `factor/gmp-ecm/microecm_edw.c` / `.h` | 64-bit scalar arithmetic layer + kernel + dispatch, and the list entry that picks the vector backend |
  | `factor/gmp-ecm/uecm_edw_avx512.c` | 52-bit x8 layer (IFMA, AVX-512F FP52) + kernel |
  | `factor/gmp-ecm/uecm_edw_avx2.c` | 52-bit x4 layer (AVX2 FP52) + kernel |
  | `factor/gmp-ecm/tinyecm_edw.c` / `.h` | 128-bit scalar layer + kernel + dispatch |
  | `factor/gmp-ecm/tecm_edw_avx512.c`, `tecm_edw_avx2.c` | 104-bit x8 / x4 layers + kernel |
  | `arith/edwards/tiny_*.py` | offline tools: chain search, curve search, yield, oracle |

* **One kernel, six backends** (scalar64, x8-52, x4-52, scalar128,
  x8-104, x4-104). The residue type is a value type (`uint64_t`,
  `__m512i`, `__m256i`, or a small struct for the 2-word types). If 2-word
  struct codegen is poor, split into a 1-word and a 2-word kernel.
* **API** (new names, same shapes as today):
  - `getfactor_uecm_edw(n, is_arbitrary, pran)`
  - `getfactor_uecm_edw_list(q, f, num, pran)`: x8 / x4 / scalar by ISA,
    vectors only for n < 2^52
  - `getfactor_tecm_edw(n, f, target_bits, pran)`
  - `getfactor_tecm_edw_x(n, f, target_bits, pran)`: one n, 8 or 4 curves
  - `getfactor_tecm_edw_list(...)`

  Plus a compile-time or runtime switch for the PRAC control path.
* **Curves:** a table of K >= 64 small-height curves per family.
  - Lane / curve c of an input uses table entry (c + offset) mod K. The
    offset comes from `*pran`, so repeated calls on the same n still see
    different curves.
  - Every input needs at most 48 curves, so K = 64 never repeats a curve
    within one input.
* **AVX2 rounding:** FP52 needs exact hi / lo products. Evaluate both:
  - (a) MXCSR round-toward-zero around the batch (`a2m` approach;
    `-frounding-math`, `fenv_access` on MSVC);
  - (b) round-to-nearest hi with a signed lo correction (~3 extra integer
    ops, no global state).

  Prefer (b) if it is within ~5%: these kernels sit inside SIQS / NFS
  cofactorization threads, where a global rounding mode is a trap.
  Nothing double-based (PRAC's `v`, the LCG) may run under RZ.

## Phase 0: model, tools, microbenchmarks (gate)

1. **Exact scalars:** `tiny_chains.py` builds each level's s exactly as
   `uecm_stage1` does (doublings to 4 B1, the composite merges).
2. **Chain search** (`tiny_chains.py`): per level, the cheapest
   signed-digit chain under the cost model.
   - Costs: dbl 3M + 4S, dbl + T 4M + 4S, add 7M, cached table entry
     ~8M + 1M. Price with S = M (micro) and S = 0.8 M (tiny), and later
     with the measured ns.
   - Candidates: w-NAF and sliding windows (w = 2..6); fractional
     windows; reordered / regrouped prime powers; a small BK-style search
     for low-weight multiples (2^a +- 2^b +- 1) covering parts of s.
   - Emit C digit tables. Report M per bit against PRAC's counted chains
     (instrument the old code's `stg1Add` / `stg1Doub`).
3. **Curve search** (`tiny_curves.py`):
   - family (iv) curves with small g = a/b (e = (g - 1/g)/2) and a small
     non-torsion point, by direct search;
   - the a = 1 Z/2 x Z/8 and Z/12 small-parameter curves of EECM §6.
   - Verify every curve symbolically (on the curve, torsion structure,
     point not torsion), the `edw_symbolic.py` way. Keep numerators and
     denominators <= 64 bits where possible (multiword is fine: a few
     extra M to reduce).
4. **Yield** (extend `edw_ref.py yield`): 16-32-bit p, each level's B1 and
   B2, and the real pairing coverage (parse the `b1_*` tables). Compare
   Suyama random sigma, table family (iv), table Z/2 x Z/8 and table
   Z/12. Paired primes, ~10^4 per point.
5. **C microbenchmarks** (standalone, `factor/gmp-ecm/bench/` or
   `arith/edwards/`):
   - latency and throughput of `mulredc` per backend;
   - old udup / uadd vs Edwards dbl / dbl+T / add;
   - **full stage 1 per level**: the old PRAC code (copied), PRAC with 2
     interleaved curves, and an Edwards prototype on the Phase-0 chain;
   - scalar64 first, then x8 IFMA (with the madd52lo m) and scalar128.

**Gate** (per backend): go if Edwards stage 1 is >= ~10% faster, or if
the predicted time per factor (stage 1 + stage 2 + build, times yield) is
>= 8% better. A backend that fails the gate keeps its PRAC stage 1 but can
still use the table curves and inversion-free build if those win alone.

## Phase 1: scalar microecm_edw (64-bit)

* Arithmetic layer from `microecm.c` (`uecm_mulredc` / asm, `addmod`,
  `submod`, `bingcd64`).
* Table curve load, inversion-free (§A2.4). Degenerate curves mod p
  (denominator, d or 1 + d = 0 mod p) show up as zero coordinates and are
  caught by the gcd as a stage-0 / stage-1 find, like AVX-ECM.
* Chain walker over the digit table. T only before an addition; the
  unified add-2008-hwcd-3, so the identity and equal inputs are safe.
* Stage-1 gcd on X (catches the 2-torsion point too). Hand-off (Z + Y,
  Z - Y). Stage 2 with projective A24.
* Debug hook: per curve, (table index, stage of find), for the oracle.

## Phase 2: microecm_edw x8 (AVX-512, n < 2^52)

* 52-bit layer: IFMA (`madd52lo/hi`, madd52lo for m) and AVX-512F FP52.
* List driver with lane refill (as `microecm_x8_list`), per-lane curve
  index.
* Interleave independent dbl products explicitly if the compiler does not.

## Phase 3: tinyecm_edw scalar (128-bit) and x8 (104-bit)

* Scalar: `mulmod128` / `sqrmod128` (check that sqr is really cheaper),
  no mpz anywhere in the curve loop.
* x8: port `tecm_mulredc104_x8` and add an IFMA version (4 madd52 pairs +
  reduction); a dedicated sqr (3 partial products).
* Both `_x` (one n, 8 curves; the `getfactor_tecm_x8` use in
  batch_factor) and `_list`.

## Phase 4: AVX2 x4 (52-bit and 104-bit)

* FP52 layer from `a2m_fp52.c`: bias-trick conversions, compare-based
  masks (no `vpcmpuq` / mask registers: use `vpcmpgtq` with a sign flip,
  and blends), rounding variant (a) or (b) per B0.
* microecm x4 and tinyecm x4 on the same kernel, plus the PRAC + Suyama
  control on the same layer. That is also the plain "AVX2 port of today's
  microecm / tinyecm", so both questions get answered.
* Check with objdump: no zmm in the AVX2 objects.

## Phase 5: verification

* **Oracle:** `edw_ref.py` tiny mode predicts per (n, table index) the
  stage of the find (0 / 1 / 2 / none) from exact ord P mod p, the
  level's s and the pairing coverage. Compare against the debug hook for
  every backend: 0 mismatches required.
  - The PRAC control is checked against the same oracle with
    `prac_stage1_z` (the x-only collapse of Part D2 exists here too).
* **Arithmetic:** random and chained mul / sqr / add / sub against GMP per
  backend, at n near 2^52, 2^64, 2^104, 2^128 and small n.
* **Miss rates and correctness** with the `test/layer2/test_ecm.c` inputs
  (40 / 50 / 60 / 64-bit semiprimes, 3LP splits up to 105 bits): equal or
  fewer misses than the old code.

## Phase 6: tuning and measurement

* Per level: chain (from measured ns, not the model), and K.
* **Re-tune the dispatch:** with cheaper curves and better yield, the best
  B1 and curve count per bit size shift. Fit the expected time per factor
  over the levels (and B2 = 25 vs 50 B1).
* Benchmarks: the `top/test.c` markers (`microecm`, `uecm_52_list`,
  `tinyecm_128`, `tinyecm_104_list`) duplicated for the new entry points,
  1 thread and all threads. Report old vs new, per ISA build (IFMA,
  AVX-512F, AVX2, and no-AVX scalar).
* Then real workloads: SIQS DLP/TLP batch cofactorization and `factor()`
  on a set of composites, total time.

## Phase 7: integration (after review)

* A switch in the callers (`batch_factor.c`, `tdiv.c`, `autofactor.c`,
  `calc.c`), default decided from Phase 6.
* MSVC project files (`build.vc22`), `Makefile` sources, the `-frounding-math`
  / `fenv_access` handling if variant (a) wins.
* Tests in `test/layer2/test_ecm.c`; results written up in this file.
* Later, separately: the lasieve5 copies, and GPU `cuda_tinyecm.cu`.

## Decisions taken as defaults (say if you want otherwise)

1. The AVX2 work includes both the Edwards version and the PRAC / Suyama
   control on the same x4 layer.
2. New versions use a fixed curve table with a `*pran`-derived offset,
   not random sigma.
3. The new versions may re-tune B1 levels and curve counts. Their tables
   start as copies of today's, so the first comparisons are like for like.
4. The tinyecm x8 IFMA 104-bit product and the madd52lo m-trick go only
   into the new code. Back-porting them to the old files would be a
   separate small change.

## Risks

* **The latency argument may not hold for scalar64.** Out-of-order cores
  may already overlap PRAC chains well enough. Phase 0 tells; scalar
  microecm is the strongest incumbent.
* **Small-height table curves might yield less** than random family (iv)
  members. The yield step measures it, and the build trick also works
  with mod-n-computed curves (one batched inversion).
* **2-word struct codegen** in one shared kernel. Fallback: split
  kernels.
* **AVX2 x4 vs scalar mulx at <= 52 bits.** It may only tie per core. It
  still frees the scalar ports in the hyperthreaded sieve workloads, but
  that is hard to measure. Judge it on total cofactorization time.

---

# Part C: Phase 0 results (2026-11-09)

Tools, all in `arith/edwards/` (nothing in yafu's build changed):

| file | what |
|---|---|
| `tiny_chains.py` | replays the existing x-only stage 1 on integer multiples (doublings to 4 B1, `uprac`, the `uprac70` / `uprac85` step tables parsed from `microecm.c`): exact s and PRAC counts per level; DP search for the cheapest Edwards chain; `emit FILE.h` writes the chains |
| `edw_tiny_chains.h` | the emitted chains (digit tables, table size m, top digit, tz) |
| `tiny_curves.py` | small-height curve search (`search`) and the yield comparison at the dispatch sizes (`yield`, `yieldall`), with stage 2 = the exact (A, b) pairs of microecm's tables |
| `tiny_bench.c` | microecm 64-bit scalar: PRAC (microecm.c's own) vs Edwards, x1 and 2-way interleaved, curve build, stage 2 |
| `tiny_bench8.c` | microecm x8 (AVX-512, n < 2^52), same; IFMA, AVX-512F, and IFMA with `madd52lo` for m |
| `tiny_bench128.c` | tinyecm scalar 128-bit: PRAC (tinyecm.c's own, monty.c arithmetic) vs Edwards; `-DINLINE128` swaps both to an inline `__int128` CIOS |
| `timing_phase0.txt`, `yield_tiny.txt` | raw output |

Every bench has `-t`. The check: run the Edwards chain on the odd part of
s, hand off to (Z + Y : Z - Y), and do tz x-only doublings with projective
A24 = dd/(dd + dn). The result must give the same Montgomery x as the
existing PRAC stage 1 on the image curve.
* 0 mismatches everywhere: 2100 scalar64 curves, 3360 x8 lanes per ISA,
  700 x 128-bit.
* That also confirms the replayed scalars are exactly the ones the code
  uses. tinyecm's step tables are identical to microecm's.

Build: WSL gcc 11.4, `-O2 -march=native` (Rocket Lake i5-11600K), 1
thread, `taskset -c 3`, min of 3.

## C1. Scalars and multiplication counts

| B1 | bits of s | PRAC adds / dups | PRAC M | Edwards m | dbls / adds | Edwards M | PRAC / Edw (M) |
|---:|---:|---|---:|---:|---|---:|---:|
| 27 | 39 | 30 / 21 | 285 | 7 | 35 / 7 | 341 | 0.84 |
| 47 | 67 | 60 / 32 | 520 | 7 | 64 / 12 | 582 | 0.89 |
| 70 | 98 | 95 / 43 | 785 | 7 | 95 / 18 | 845 | 0.93 |
| 85 | 125 | 128 / 48 | 1008 | 7 | 123 / 23 | 1081 | 0.93 |
| 125 | 178 | 188 / 65 | 1453 | 9 | 174 / 29 | 1498 | 0.97 |
| 165 | 236 | 262 / 76 | 1952 | 7 | 232 / 44 | 2010 | 0.97 |
| 205 | 300 | 341 / 90 | 2496 | 13 | 296 / 49 | 2534 | 0.98 |

Cost model (S = M):
* PRAC: add 6M, dup 5M.
* Edwards: dbl 7M, +1M for T before an add, add 7M; a cached table entry
  costs 4M including the dd scaling.
* The power of 2 in s is done after the hand-off as x-only doublings: 5M
  each instead of 7M.

The deep-searched PRAC chains are very good here (~1.0 add + 0.38 dup per
bit, ~8 M/bit). **Edwards needs 2-16% more multiplications.** That is
worse than the Part A estimate: short scalars make the table expensive,
and S = M.

## C2. Stage 1 time anyway: Edwards wins on ILP

ns per curve (x8: per vector of 8 curves), stage 1 only:

| B1 | scalar64 PRAC | Edw | ratio | PRAC x2 | Edw x2 | ratio | x8 IFMA PRAC | Edw | ratio | x8 IFMA+m52 PRAC | Edw | ratio | x8 AVX-512F PRAC | Edw | ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 27 | 495 | 438 | 1.13 | 402 | 388 | 1.04 | 1165 | 926 | 1.26 | 765 | 625 | 1.23 | 1746 | 1346 | 1.30 |
| 47 | 883 | 734 | 1.20 | 744 | 647 | 1.15 | 2073 | 1573 | 1.32 | 1357 | 1075 | 1.26 | 3113 | 2273 | 1.37 |
| 70 | 1286 | 1058 | 1.22 | 1009 | 935 | 1.08 | 3063 | 2270 | 1.35 | 1997 | 1544 | 1.29 | 4606 | 3272 | 1.41 |
| 85 | 1667 | 1346 | 1.24 | 1307 | 1209 | 1.08 | 3930 | 2893 | 1.36 | 2566 | 1962 | 1.31 | 5956 | 4173 | 1.43 |
| 125 | 2408 | 1856 | 1.30 | 1940 | 1649 | 1.18 | 5636 | 3973 | 1.42 | 3670 | 2691 | 1.36 | 8494 | 5726 | 1.48 |
| 165 | 3271 | 2495 | 1.31 | 2667 | 2232 | 1.19 | 7598 | 5350 | 1.42 | 4928 | 3597 | 1.37 | 11354 | 7669 | 1.48 |
| 205 | 4213 | 3133 | 1.34 | 3458 | 2803 | 1.23 | 9730 | 6705 | 1.45 | 6300 | 4517 | 1.40 | 14620 | 9608 | 1.52 |

(scalar64 at 64-bit n; 50-bit n is within 1%.)

* PRAC does fewer multiplications, yet it is 1.13-1.34x slower in the
  scalar code and 1.23-1.52x slower in x8: **the latency argument holds.**
  The Edwards doubling's independent squarings and products overlap;
  PRAC's chain doesn't.
* 2-way interleaving helps PRAC more than Edwards (it has more idle
  latency), but Edwards x2 still wins 1.04-1.23x. Edwards x1 is about as
  fast as PRAC x2, without needing two curves in flight.
* **x8 IFMA with `madd52lo` for m** (instead of `vpmullq` + mask):
  - the existing PRAC stage 1 gets 1.52x faster, and so does stage 2
    (1372 -> 902 ns at B1 = 27);
  - with it, the Edwards stage 1 is 1.86-2.15x faster than today's x8
    stage 1.
* **tinyecm scalar 128-bit:**
  - with monty.c's arithmetic, Edwards is 0.87x at B1 = 27, ~0.95-0.99x
    at 47-85, and 1.04-1.06x at 125-205;
  - with an inline `__int128` CIOS it is about the same (0.83-1.04x);
  - the 2-word multiply is throughput-bound enough that the extra
    multiplications cost what the ILP saves.
* Side findings in monty.c:
  - `addmod128` (asm path) is wrong when the output aliases its second
    operand. tinyecm doesn't appear to hit that case today.
  - The bench built with `-flto` gave mismatches everywhere. The cause
    was not investigated (likely the inline asm once it is inlined).
  - So the new code should not depend on monty.c.

## C3. Curve construction

| | Suyama build (existing) | load of precomputed constants |
|---|---:|---:|
| scalar64 | 189-236 ns (1 `uecm_modinv_64`) | 7 ns |
| x8 IFMA, per 8 curves | 1349-1382 ns (8 scalar modinv) | 6-9 ns |
| scalar128 | 1400-1680 ns (`mpz_invert` + mpz conversions) | ~690 ns with monty.c's mpz `to_monty128`; ~6 mulmods with R^2 |

The Suyama build is ~50% of a B1 = 27 stage 1 in scalar64, and **more
than a whole B1 = 27 stage 1 in x8.** At small B1, inversion-free curve
construction is worth as much as the Edwards arithmetic.

**But a small rational table is not available.** Family (iv) curves with a
non-torsion point of small height are rare:
* Direct search over g = a/b with a, b <= 100 and points x = xn/xd <=
  1000 finds only 8 curves. All have constants <= 29 bits; none more turn
  up even when constants up to 127 bits are allowed.
* GKL's [m]Q (with -[m]Q and the torsion translates on E') gives 2 curves
  at 20 bits, 2 at 59, then 117+ bits.
* The a = +1 families have more: 30 Z/12 curves and 13 Z/2 x Z/8 curves
  with constants <= 63 bits (`curves_z12.json`, `curves_z2x8.json`).

So microecm and tinyecm need ~48 curves per input from a construction
mod n, not a table (§C6).

## C4. Yield at the dispatch sizes

Setup: one random curve per random prime p, paired, 40000 primes per row.
Exact ord P mod p, microecm's s and its exact stage-2 (A, b) pairs. Finds
per curve relative to Suyama (std error of the paired difference ~1%):

| p bits | B1 | Suyama finds | GKL (iv), random m | table (iv), 8 curves | table Z/12 | table Z/2 x Z/8 |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 27 | 10713 | **1.141** | 1.122 | 1.085 | 1.097 |
| 22 | 47 | 10380 | **1.111** | 1.100 | 1.068 | 1.082 |
| 24 | 70 | 8783 | **1.142** | 1.092 | 1.077 | 1.101 |
| 26 | 85 | 6657 | **1.119** | 1.106 | 1.074 | 1.119 |
| 29 | 125 | 5879 | **1.098** | 1.124 | 1.080 | 1.085 |
| 31 | 165 | 4873 | **1.091** | 1.086 | 1.069 | 1.077 |
| 32 | 205 | 4753 | **1.125** | 1.122 | 1.090 | 1.103 |

* Every Edwards family beats Suyama by 7-14% at these sizes. That is more
  than at AVX-ECM sizes (+5-10% in `avxecm_edwards_plan.md` Part C1).
* GKL family (iv) with random m is the best or tied everywhere:
  **+9-14%**.
* Z/12 is the weakest (+7-9%). Z/2 x Z/8 is in between and would need the
  a = +1 addition. **Family (iv) it is.**

## C5. Time per factor, and the gate

Per curve: build + stage 1 + stage 2 (stage 2 unchanged; ns):

| | B1 | today | Edwards (+ m52 for x8) | speed | x yield (iv) | **per factor** |
|---|---:|---:|---:|---:|---:|---:|
| scalar64 | 27 | 236 + 495 + 644 = 1375 | 75 + 438 + 644 = 1157 | 1.19 | 1.14 | **1.36x** |
| scalar64 | 85 | 236 + 1667 + 1176 = 3079 | 75 + 1346 + 1176 = 2597 | 1.19 | 1.12 | **1.33x** |
| scalar64 | 205 | 236 + 4213 + 4416 = 8865 | 75 + 3133 + 4416 = 7624 | 1.16 | 1.13 | **1.31x** |
| x8 IFMA (per 8) | 27 | 1382 + 1165 + 1372 = 3919 | 150 + 625 + 902 = 1677 | 2.34 | 1.14 | **2.7x** |
| x8 IFMA (per 8) | 85 | 1382 + 3930 + 2409 = 7721 | 150 + 1962 + 1629 = 3741 | 2.06 | 1.12 | **2.3x** |
| x8 IFMA (per 8) | 205 | 1382 + 9730 + 4644 = 15756 | 150 + 4517 + 3138 = 7805 | 2.02 | 1.13 | **2.3x** |
| x8 AVX-512F (per 8) | 85 | 1492 + 5956 + 3721 = 11169 | 300 + 4173 + 3721 = 8194 | 1.36 | 1.12 | **1.5x** |
| scalar128 | 85 | 1448 + 10237 + 8254 = 19939 | 600 + 10353 + 8254 = 19207 | 1.04 | 1.12 | **1.16x** |
| scalar128 | 27 | 1448 + 2836 + 4322 = 8606 | 600 + 3259 + 4322 = 8181 | 1.05 | 1.14 | **1.20x** |

Notes:
* The Edwards build is an estimate: ~50 multiplications per curve for the
  inversion-free GKL construction of §C6. That is 75 ns in scalar64, ~20
  ns per lane in x8, and ~600 ns in scalar128.
* The 8 stage-2 doublings with projective A24 add ~1%.
* Of the x8 IFMA gain, `madd52lo` alone (keeping PRAC and Suyama) is
  ~1.3x.

**Gate: passed.**
* microecm scalar64: ~1.3x less time per factor, 1.16-1.19x of it from
  speed.
* microecm x8: ~2.3-2.7x on IFMA (1.3x of it from `madd52lo`), ~1.5x on
  AVX-512F.
* tinyecm scalar128: ~1.15-1.2x, almost all from yield and the cheaper
  build; stage 1 itself is a wash. Use Edwards anyway, to keep one code
  path, but the curve construction is what matters there.
* tinyecm x8 104-bit: ~2.2-2.6x on IFMA with an IFMA 104-bit multiply,
  ~1.35x on AVX-512F (§C7).
* Not measured yet (its own phase): AVX2 x4.

## C6. Changes to the plan

1. **Curves: inversion-free GKL family (iv) mod n**, not a table:
   - keep [m]Q on E' in projective coordinates and step m -> m + 1 with
     one projective addition per curve (~14M);
   - get t = tn/td = (4X + 4Z)/(Y - 4Z), then e, x1, y1 as homogeneous
     polynomials in (tn, td), and from them the Edwards point
     (X, Y, Z, T) and d = dn/dd, all without dividing (~35M);
   - a random start m0 per input (a short projective scalar multiple once
     per n, e.g. m0 < 2^8), so retries see other curves;
   - optional: the 8 small-height curves (`curves_iv.json`) as the first
     curves, loaded for free.
2. **x8 IFMA `madd52lo` m** is mandatory in the new code. It is also a
   one-line, ~1.3x win for the *existing* `microecm.c` x8 (stage 1 and
   stage 2) that doesn't change any results, so it is worth doing
   separately.
3. **Power of 2 last, in Montgomery form:**
   - the chain covers only the odd part of s;
   - after the hand-off come tz x-only doublings with projective A24,
     which are cheaper than Edwards doublings;
   - the stage-1 gcd is then on Z_M. A 2-torsion point becomes the
     identity after the doublings, so it is still caught.
4. **tinyecm scalar128:**
   - use Edwards for a single code path, but don't expect a stage-1 gain;
   - spend the effort on the build (no mpz in the curve loop) and an
     inline 128-bit layer;
   - re-check per level once that layer is tuned (sqr with 3 partial
     products, branch-free reductions).
5. **Chains:** DP-optimal signed digits are within 1-2% of each other for
   m = 5..13, so the choice per level is not critical. The tables in
   `edw_tiny_chains.h` are ready to use.
6. **tinyecm x8: an all-IFMA 104-bit multiply** (§C7) is mandatory in the
   new code. Like `madd52lo` in microecm x8, it would also be a
   result-neutral ~1.5x for the existing `tinyecm.c` x8.

## C7. tinyecm x8 (104-bit, AVX-512)

Added after the first write-up. `tiny_bench104.c` measures it the same way
(`-t` gives 0 of 1680 lanes mismatching for each build).

`tecm_mulredc104_x8` today does all its 52x52 products as FP
`fmadd_round` pairs, even in IFMA builds. Only two of the m cross terms
use `madd52lo`. `tiny_ifma104.py` writes a copy of tinyecm.c whose
multiply uses `madd52lo/hi` for every product (19 madd52 instructions
per mulredc, same reduction and final correction).

100-bit n, ns per vector of 8 curves:

| B1 | today: PRAC | Edw | ratio | stage 2 | IFMA-104: PRAC | Edw | ratio | stage 2 | AVX-512F: PRAC | Edw | ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 27 | 5660 | 5836 | 0.97 | 8316 | 3884 | 3340 | 1.16 | 5568 | 6031 | 6087 | 0.99 |
| 47 | 10391 | 9969 | 1.04 | 10859 | 7126 | 5630 | 1.27 | 7288 | 11032 | 10314 | 1.07 |
| 70 | 14519 | 14319 | 1.01 | 13934 | 9560 | 8090 | 1.18 | 9370 | 15448 | 14907 | 1.04 |
| 85 | 18609 | 18112 | 1.03 | 15940 | 11972 | 10326 | 1.16 | 10608 | 19788 | 19058 | 1.04 |
| 125 | 27661 | 25002 | 1.11 | 21413 | 18289 | 14319 | 1.28 | 14287 | 29394 | 26552 | 1.11 |
| 165 | 37976 | 33787 | 1.12 | 27587 | 25450 | 19168 | 1.33 | 18323 | 40245 | 35282 | 1.14 |
| 205 | 49292 | 42137 | 1.17 | 32611 | 33280 | 23998 | 1.39 | 21710 | 52345 | 44518 | 1.18 |

The Suyama build (`build_curves_104_x8`: 8 `mpz_invert` + mpz Montgomery
conversions) takes **7.5-7.9 us per vector**. That is more than a whole
stage 1 at B1 <= 70 on IFMA.

* **The FP 104-bit multiply is throughput-bound**: with it, Edwards is a
  wash at small B1 and 1.1-1.17x at B1 >= 125. That is the same picture as
  scalar128, and the same on AVX-512F.
* **The IFMA 104-bit multiply** makes the existing PRAC stage 1 1.46-1.53x
  faster and stage 2 ~1.5x faster. With it the multiply is fast enough
  that latency matters again: Edwards then wins 1.16-1.39x.
* Edwards + IFMA-104 against today's stage 1: **1.69x (B1 = 27) to 2.05x
  (B1 = 205).**

Time per factor (build + stage 1 + stage 2, x yield). The Edwards build is
estimated as ~50 multiplications per curve, ~0.5 us per vector on IFMA and
~0.9 us on AVX-512F.

| | B1 | today | new | speed | **per factor** |
|---|---:|---:|---:|---:|---:|
| IFMA | 27 | 7890 + 5660 + 8316 = 21866 | 475 + 3340 + 5568 = 9383 | 2.33 | **2.65x** |
| IFMA | 85 | 7890 + 18609 + 15940 = 42439 | 475 + 10326 + 10608 = 21409 | 1.98 | **2.2x** |
| IFMA | 205 | 7890 + 49292 + 32611 = 89793 | 475 + 23998 + 21710 = 46183 | 1.94 | **2.2x** |
| AVX-512F | 85 | 7808 + 19788 + 16870 = 44466 | 880 + 19058 + 16870 = 36808 | 1.21 | **1.35x** |
| AVX-512F | 205 | 7808 + 52345 + 34294 = 94447 | 880 + 44518 + 34294 = 79692 | 1.19 | **1.34x** |

Of the IFMA gain, the IFMA multiply alone (keeping PRAC and the Suyama
build) is ~1.4x at B1 = 85.

---

# Part D: optimization pass of the existing microecm / tinyecm (2026-11-10)

Before Phase 1, the arithmetic improvements Phase 0 found went into the
existing files directly. Same curves, same chains, same results.

| file | change |
|---|---|
| `microecm.c` x8 | m = T_lo * invN mod 2^52 with `madd52lo` (IFMA) instead of `vpmullq` + mask, in `uecm_mulredc_x8` and `uecm_sqrredc_x8` |
| `microecm.c` scalar | the gcc inline-asm `uecm_mulredc` / `uecm_sqrredc` are used on any x86-64 gcc / clang except clang-cl. They required `__unix__`, so MinGW / MSYS2 builds used `uecm_mulredc_alt`, ~5% slower in stage 1 and 5-10% in stage 2 (measured on Linux by forcing that path) |
| `tinyecm.c` x8 | with IFMA, `tecm_mulredc104_x8` does all its 52x52 products with `madd52lo/hi` (19 instructions) instead of FP `fmadd_round` pairs; the FP version stays for AVX-512F |
| `tinyecm.c` x8 build | R^2 mod n per lane in `monty104_x8_t.rsqr`; Montgomery conversions are mulredcs by R^2 instead of mpz per lane; one `mpz_invert` per lane, no mpz products or conversions after it |
| `tinyecm.c` scalar | inline `__int128` CIOS `mulmod128` / `sqrmod128` / `addmod128` / `submod128` (where `__SIZEOF_INT128__`; MSVC keeps monty.c; `-DTECM_NO_INLINE128` forces it), replacing out-of-line calls into monty.c |
| `tinyecm.c` scalar build | R^2 in the work struct, conversions by mulmod; one `mpz_invert` |
| `monty.c` | bug fix: the asm `addmod128` clobbered its 2nd operand when the output aliased it (`arith/limb2.c` has the same code but isn't built) |

Tried and dropped: `vpmuludq` pieces for the AVX-512F m (no change).

**Verification:**
* The Phase 0 exact checks pass for all four benches on IFMA and
  AVX-512F builds.
* `tiny_e2e.c` hashes every factor returned through the four public
  entry points (18 rows, same seeds). The hashes are identical before and
  after on both builds, so results are unchanged to the bit.

**Benchmarks** (Rocket Lake, 1 thread, `taskset`, WSL gcc 11.4; the
"before" binaries are built from the previous commit by
`arith/edwards/tiny_build.sh`; raw output in the session scratchpad):

Stage level, ns, before -> after:

| B1 | x8 52-bit IFMA stage 1 | stage 2 | x8 104-bit IFMA stage 1 | stage 2 | scalar128 stage 1 | stage 2 |
|---:|---|---|---|---|---|---|
| 27 | 1151 -> 775 (1.49x) | 1358 -> 896 (1.52x) | 5640 -> 3889 (1.45x) | 8326 -> 5759 (1.45x) | 2896 -> 2290 (1.26x) | 4316 -> 4032 (1.07x) |
| 47 | 2086 -> 1345 (1.55x) | 1699 -> 1128 (1.51x) | 10423 -> 7208 (1.45x) | 10883 -> 7697 (1.41x) | 5382 -> 4275 (1.26x) | 6045 -> 5299 (1.14x) |
| 70 | 3061 -> 2002 (1.53x) | 2123 -> 1416 (1.50x) | 14462 -> 9596 (1.51x) | 13907 -> 9837 (1.41x) | 7931 -> 6278 (1.26x) | 7244 -> 6857 (1.06x) |
| 85 | 3955 -> 2546 (1.55x) | 2417 -> 1601 (1.51x) | 18681 -> 12100 (1.54x) | 15978 -> 11363 (1.41x) | 10208 -> 8266 (1.23x) | 8150 -> 7763 (1.05x) |
| 125 | 5647 -> 3666 (1.54x) | 3147 -> 2129 (1.48x) | 27696 -> 18498 (1.50x) | 21548 -> 15227 (1.42x) | 14773 -> 12029 (1.23x) | 10989 -> 10533 (1.04x) |
| 165 | 7579 -> 4940 (1.53x) | 4055 -> 2738 (1.48x) | 38355 -> 25453 (1.51x) | 27710 -> 19762 (1.40x) | 20081 -> 16084 (1.25x) | 14114 -> 13706 (1.03x) |
| 205 | 9773 -> 6295 (1.55x) | 4788 -> 3127 (1.53x) | 49168 -> 33286 (1.48x) | 32678 -> 23262 (1.40x) | 25970 -> 20624 (1.26x) | 16621 -> 15970 (1.04x) |

Curve builds:
* x8 104-bit: 7652 -> 479 ns per 8 curves IFMA (16x); 7834 -> 661 ns
  AVX-512F (12x).
* scalar128: 1418 -> 163 ns (8.7x).

Unchanged within noise (±2%): scalar64 on Linux (no code change there),
and the AVX-512F stage-level numbers except the x8 104-bit build.

End to end through the public entry points, us per input (10000 random
inputs per row; uecm: balanced semiprimes; tecm: p of t bits times q of
2t bits, target t). Found rates are identical before and after.

| level | uecm x8 list, IFMA | tecm scalar, IFMA | tecm x8 list, IFMA | tecm x8 list, AVX-512F |
|---|---|---|---|---|
| B1 27 (40-bit / t = 20) | 1.23 -> 0.97 (1.27x) | 30.4 -> 26.0 (1.17x) | 9.99 -> 6.51 (1.53x) | 10.28 -> 8.98 (1.14x) |
| B1 47 (44 / 22) | 1.64 -> 1.26 (1.30x) | 47.4 -> 41.9 (1.13x) | 15.46 -> 10.18 (1.52x) | 15.99 -> 14.15 (1.13x) |
| B1 70 (48 / 24) | 2.49 -> 1.81 (1.38x) | 87.9 -> 66.1 (1.33x)* | 22.84 -> 15.10 (1.51x) | 23.60 -> 21.45 (1.10x) |
| B1 85 (52 / 26) | 3.73 -> 2.70 (1.38x) | 126.2 -> 104.7 (1.20x) | 34.18 -> 22.62 (1.51x) | 35.97 -> 33.21 (1.08x) |
| B1 125 (t = 29) | | 252.0 -> 208.9 (1.21x) | 67.68 -> 44.01 (1.54x) | 68.23 -> 65.18 (1.05x) |
| B1 165 (t = 31) | | 404.5 -> 331.9 (1.22x) | 116.74 -> 70.67 (1.65x) | 107.86 -> 102.65 (1.05x) |
| B1 205 (t = 32) | | 488.8 -> 407.8 (1.20x) | 123.91 -> 85.72 (1.45x) | 130.15 -> 125.70 (1.04x) |

(* the 87.9 before is a noisy run; the AVX-512F build's scalar tecm, same
code, went 76.6 -> 66.8, 1.15x. In general scalar tecm is 1.13-1.22x on
both builds. scalar uecm: 0.98-1.03x, noise.)

**Effect on the Edwards plan:**
* The incumbents are now faster, so Part C5's per-factor gains shrink to
  the Edwards-only part. Against the optimized code (same method as C5:
  build + stage 1 + stage 2, times yield):
  - microecm x8 IFMA: ~1.55x (B1 = 205) to ~2.0x (B1 = 27). Its Suyama
    build (8 scalar modinv, 1351 ns) is untouched and is now larger than
    the whole stage 1 at B1 = 27.
  - tinyecm x8 IFMA: ~1.2-1.4x.
  - tinyecm scalar: ~1.05-1.1x. Its build is now cheap, so this is mostly
    yield.
  - microecm scalar: unchanged, ~1.3x.
* Phase 0's Edwards prototypes already used these arithmetic forms, so
  the Edwards stage-1 numbers stand.

---

# Part E: Phase 1, scalar microecm_edw (2026-11-10)

## E1. What was built

`factor/gmp-ecm/microecm_edw.c` / `.h` is a standalone 64-bit scalar
microecm on Edwards curves. `microecm.c` is untouched by it. It is in the
Makefile and both VS `factor_common` projects; no caller uses it yet.

* **API:**
  - `getfactor_uecm_edw(q64, is_arbitrary, pran)`: same contract, B1
    levels, curve counts, arbitrary pre-passes and B2 as `getfactor_uecm`;
  - test hooks `uecm_edw_one_curve`, `uecm_edw_curve_params` and
    `uecm_edw_stage1_point`.
* **Arithmetic:** microecm's (asm `mulredc`, `addmod`, `submod`, binary
  gcd), copied with `uedw_` names; portable / MSVC fallbacks.
* **Curves:** GKL family (iv) without inversions (§C6):
  - [m]Q on E' in projective coordinates: a mixed addition of Q = (5, 8)
    per curve, ~10M; a projective doubling for the start;
  - t = tn/td, then e, x1, y1 as homogeneous polynomials and the Edwards
    point (X : Y : Z : T) with d = dn/dd, ~36M;
  - the formulas were checked in Python against `edw_ref.py curve_iv`
    before coding;
  - the start m0 is in [2, 257], from the caller's `*pran`, computed once
    per input (double-and-add) and shared by the arbitrary pre-passes and
    the main run.
* **Stage 1:** the `edw_tiny_chains.h` chain on the odd part of s
  (extended coordinates, cached table scaled by dd), the hand-off (Z + Y :
  Z - Y), and tz x-only doublings with projective A24 = dd / (dd + dn).
  The gcd is taken on Z.
* **Stage 2:** microecm's `uecm_stage2`, D = 60, with the same tables
  (copied), with the 4-5 doublings using the projective A24.

## E2. Verification

`arith/edwards/uedw_test.c` (stdin driver) and `uedw_oracle.py`:

* `params`: the curve for m = 1..40 mod random 40-64-bit primes, out of
  Montgomery form, against `edw_ref.py curve_iv`. **0 of 2000 mismatch.**
* `curves`: one curve per random (n = p q, m) at every level. For each
  prime the oracle takes:
  - the exact stage-1 point (x-only ladder by s);
  - every stage-2 cross term mod p, from an exact Python replay of
    microecm's stage 2;
  - then it replays the C accumulator loops. Each loop stops multiplying
    when a term would make acc 0 mod n.

  That control flow matters, and the first oracle got it wrong. When both
  primes hit, the earlier term wins. And the `break` in the B1 <= 70
  pre-step loop leaves only that loop, so the main loop can still find one
  prime after both hit together. The replay also covers the x-only
  collapse when a leftover order divides a difference used (r | 30 etc.).

  **0 mismatches in 37122 curves** (two seeds, all seven levels: 1253
  stage-1 finds, 11436 stage-2 finds, 24433 none). **3 extra stage-1
  finds** come from exceptional cases of the Edwards formulas: the
  stage-1 point is (0 : 0) mod that prime, found by the gcd, a correct
  factor (§A6 of the AVX-ECM plan). The oracle recognizes and reports
  these.

## E3. Speed

`tiny_e2e.c -DWITH_EDW`: the same 10000 balanced semiprimes per row
through `getfactor_uecm` (with the Part D code) and `getfactor_uecm_edw`.
Found rates were 1.0 for both (`getfactor_uecm` missed 0.02% at 58 and 64
bits). us per input:

| bits | B1 | uecm, IFMA build | uecm_edw | speedup | uecm, AVX-512F build | uecm_edw | speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 27 | 2.99 | 2.81 | 1.07 | 2.97 | 2.62 | 1.13 |
| 44 | 47 | 4.03 | 3.63 | 1.11 | 4.16 | 3.96 | 1.05 |
| 48 | 70 | 6.42 | 5.33 | 1.20 | 7.07 | 5.71 | 1.24 |
| 52 | 85 | 10.17 | 8.24 | 1.24 | 10.91 | 8.66 | 1.26 |
| 58 | 125 | 21.28 | 16.96 | 1.25 | 22.40 | 18.22 | 1.23 |
| 62 | 165 | 33.29 | 26.20 | 1.27 | 34.71 | 26.26 | 1.32 |
| 64 | 205 | 40.66 | 31.68 | 1.28 | 42.53 | 32.18 | 1.32 |

(The two builds run the same scalar code; the spread between them is run
to run noise, ~3-5% per row.)

**1.2-1.3x at B1 >= 70, as predicted in C5. Only 1.05-1.13x at B1 =
27 / 47.** At the small levels the per-input start ([m0]Q by
double-and-add, up to 8 doublings and additions on E') and the per-curve
GKL build (~46M) weigh more against a ~300M stage 1. Follow-ups for those
two levels:
* a cheaper start (a smaller m0 range, or small precomputed start points);
* fewer per-curve multiplications in the build (shared powers, smaller
  constants).

## E4. Small B1: the small-constant curves first

Revisiting §C3: the search found few small-constant curves, but most tiny
inputs need only a few curves.

**Curves needed per input** (`arith/edwards/uedw_stats.c`, microecm_edw
built with `-DUEDW_STATS`, 10000 balanced semiprimes per level; curves
run, the finding one included):

| bits | B1 | mean | P(1) | P(<= 2) | P(<= 3) | P(<= 5) | P(<= 8) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 27 | 1.96 | 0.521 | 0.763 | 0.880 | 0.969 | 0.997 |
| 44 | 47 | 2.10 | 0.482 | 0.732 | 0.856 | 0.960 | 0.994 |
| 48 | 70 | 2.35 | 0.438 | 0.675 | 0.808 | 0.935 | 0.987 |
| 52 | 85 | 3.07 | 0.329 | 0.549 | 0.699 | 0.860 | 0.955 |
| 58 | 125 | 3.40 | 0.294 | 0.503 | 0.646 | 0.822 | 0.940 |
| 62 | 165 | 4.12 | 0.246 | 0.429 | 0.564 | 0.749 | 0.891 |
| 64 | 205 | 4.15 | 0.247 | 0.428 | 0.567 | 0.747 | 0.888 |

So a handful of cheap first curves covers most inputs at small B1.

**What the table holds:**
* The 8 curves of `curves_iv.json` are only 4 curves: for d = -e^4, the
  curves with d and 1/d are isomorphic (same group order) and come in
  pairs (g = 4 / 5/3, 9/2 / 11/7, 13/3 / 8/5, 17 / 9/8).
* g = 4 is GKL's m = 1.
* GKL's m = 2 fits in 59 bits.
* `arith/edwards/tiny_table.py` writes these 5 distinct curves to
  `factor/gmp-ecm/edw_tiny_curves.h` as signed 64-bit (X, Y, Z, T, dn,
  dd).

**What changed in `microecm_edw.c`:**
* A curve source per input: the 5 table curves first (6 mulredcs each to
  load), then GKL m = m0, m0 + 1, ... with m0 in [3, 258].
* [m0]Q is computed only when the table runs out, i.e. for ~3-25% of the
  inputs, depending on the level.
* Test hooks take a curve id (GKL m, or `UECM_EDW_TABLE | i`). There is
  also `uecm_edw_set_table` (experiments) and `uecm_edw_hist` (with
  `UEDW_STATS`).

**Verification:** `uedw_oracle.py` now also covers table curves:
* `params`: 0 of 1350 mismatch, table and GKL curves;
* `curves` (30% table curves): 0 mismatches in 37133 curves, two seeds,
  plus 1 extra find from an Edwards exception.

**Speed**, us per input, same inputs and run (IFMA build, 1 thread):

| bits | B1 | getfactor_uecm | edw, GKL only | edw, table first | speedup (table) | speedup (GKL only) |
|---:|---:|---:|---:|---:|---:|---:|
| 40 | 27 | 2.77 | 2.51 | 2.26 | **1.23** | 1.11 |
| 44 | 47 | 4.17 | 3.63 | 3.38 | **1.24** | 1.15 |
| 48 | 70 | 6.41 | 5.34 | 5.10 | **1.26** | 1.20 |
| 52 | 85 | 10.16 | 8.23 | 7.99 | **1.27** | 1.23 |
| 58 | 125 | 20.88 | 16.68 | 16.43 | **1.27** | 1.25 |
| 62 | 165 | 32.75 | 27.61 | 25.61 | **1.28** | 1.19 |
| 64 | 205 | 39.62 | 31.50 | 31.45 | **1.26** | 1.26 |

microecm_edw is now 1.23-1.28x faster than the optimized `getfactor_uecm`
at every level. More small curves (the a = +1 Z/12 and Z/2 x Z/8 tables)
would only help the few inputs needing more than 5 curves, and would need
the a = +1 addition. Not pursued.

---

# Part F: Phase 2, microecm_edw x8 (2026-11-10)

## F1. What was built

In `microecm_edw.c`, under `USE_AVX512F`:
`getfactor_uecm_edw_x8_list(q64, f64, num, pran)`, plus the test hook
`uecm_edw_x8_run`. Eight inputs (n < 2^52) at a time, one per lane, with
lane refills, as `getfactor_uecm_x8_list` does. It uses that function's
B1 / curve-count table.

* **Arithmetic:** R = 2^52, microecm's x8 Montgomery multiply. With IFMA
  m comes from `madd52lo`; AVX-512F uses the FP52 path.
* **Curves per lane:** the 5 table curves, then GKL m = 3, 4, ...
  - The start [3]Q = (1205/49, -40664/343) = (413315 : -1992536 : 16807)
    is a constant, so no lane needs a scalar multiplication on E'. The
    price is that the x8 sequence is the same for every input (no `*pran`
    randomization).
  - Table constants (up to 59 bits) load per lane as two 52-bit digits:
    a0 R^2/R + a1 R^3/R.
  - A vector with any lane past the table builds GKL curves for all lanes
    and blends.
* **Stage 1:** the scalar code's chains and formulas, vectorized.
* **Stage 2:** microecm's `uecm_stage2_x8`, with projective A24. Two
  differences from the scalar stage 2 are kept as they are:
  - it uses the B2 = 25 B1 tables at B1 >= 125;
  - a lane's accumulator *skips* a product that would be 0 mod n instead
    of breaking out of the loop.

## F2. Verification

`uedw_oracle.py x8` predicts, for each input, the first curve of the x8
sequence that finds a prime and which prime, under the x8 stage-2
semantics, plus the total number of curves the list runs. It feeds the
`x` command of `uedw_test` (built with `USE_AVX512F`).

| build | inputs | wrong factors | curve totals |
|---|---:|---:|---|
| IFMA | 21000 + 1400 | 0 | equal at every level but 2, off by 2 |
| AVX-512F | 21000 | 0 | equal at every level but 2, off by 1-2 |

* The input that found the *other* prime (IFMA, B1 = 27) is an Edwards
  exception: the first table curve's stage-1 point is (0 : 0) mod that
  prime. The oracle now checks for that and counts such finds separately.
* The curve totals differ by the same effect: a correct find a curve or
  two early.

## F3. Speed

`tiny_e2e.c -DWITH_EDW`: 20000 balanced semiprimes per row; us per input:

| bits | B1 | x8 list, IFMA | x8 edw, IFMA | speedup | x8 list, AVX-512F | x8 edw, AVX-512F | speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 27 | 0.97 | 0.58 | **1.66** | 1.56 | 1.05 | **1.48** |
| 44 | 47 | 1.27 | 0.79 | **1.60** | 2.26 | 1.55 | **1.46** |
| 48 | 70 | 1.87 | 1.19 | **1.57** | 3.44 | 2.30 | **1.49** |
| 52 | 85 | 2.75 | 1.77 | **1.56** | 5.01 | 3.50 | **1.43** |

* These are against the Part D (optimized) x8 code.
* Against the x8 code before Part D: IFMA 1.23 -> 0.58 us at 40 bits
  (2.1x) and 3.73 -> 1.77 us at 52 bits (2.1x).
* The same run measured scalar microecm_edw at 1.17-1.33x over
  `getfactor_uecm`.

## F4. Next

* Phase 3: tinyecm_edw, scalar 128-bit and x8 104-bit.
* Phase 4: AVX2 x4.
* Integration into callers: the user will do it later.

## F5. To do at the end of this effort: B1 / B2 choices (user request, 2026-11-10)

The levels B1 = 27, 47, 70, 85, 125, 165, 205 and the stage-2 multipliers
B2 = 25 B1 / 50 B1 are legacy values that were found to work well. They
were never optimized systematically for a given factor size. Investigate
and optimize:

* **The question:** for an expected N-bit factor, which B1 minimizes the
  expected time per factor, statistically? N = 16..32 for scalar 64-bit
  microecm_edw, and the corresponding ranges for tinyecm_edw and the x8 /
  x4 versions.
* **The same for the B2 multiplier.** 25 and 50 may not be best; the
  stage-2 cost per pair is cheaper relative to stage 1 now.
* **Method:**
  - exact success probability per curve from the yield model
    (`tiny_curves.py`: exact ord P, s, pairing coverage), for family (iv)
    including the table curves;
  - measured cost per curve per B1 / B2 from the C code;
  - minimize cost / P(success), or the full curve-count distribution with
    the per-input start costs;
  - check against end-to-end runs.
* **New B1 values need new chains and stage-2 tables:**
  - chains: `tiny_chains.py` builds them for any s;
  - stage-2 pair tables: generate them, or replace the precomputed tables
    with a small generator.
* Also revisit the dispatch thresholds (target bits -> level) and the curve
  counts.

---

# Part G: Phase 3, tinyecm_edw (2026-11-10)

## G1. What was built

* **`factor/gmp-ecm/edw_kern.h`:** the curve code, written once against a
  small arithmetic interface (value type, mul / sqr / add / sub macros, a
  context, and the stage-2 accumulation rule) and included once per
  arithmetic type. It holds:
  - E' steps: + Q, doubling, general projective addition;
  - the GKL curve build;
  - the Edwards doubling, addition and cached form;
  - the chain walker with the hand-off and the x-only power-of-2 tail;
  - the x-only stage 2.
* **`factor/gmp-ecm/tinyecm_edw.c` / `.h`** instantiate it twice.

**Scalar, 128-bit:**
* `getfactor_tecm_edw`: portable 2-word CIOS (`__int128` or `_umul128`),
  R = 2^128.
* The modular corrections are branch-free. With branches, inputs above
  2^64 ran 30% slower: the conditional subtractions mispredict once n's
  high word is nonzero.
* tinyecm's dispatch, curve counts and its own B2 = 25 B1 tables (its
  b1_165 has 425 entries, microecm's 435).
* Stage 2: a product that would make acc 0 mod n ends that loop.
* The curve source per input is the same as microecm_edw: 5 table curves,
  then GKL from a random m0.

**x8, 104-bit (AVX-512):**
* Two 52-bit digits, R = 2^104. With IFMA, all products are `madd52lo/hi`;
  AVX-512F uses the FP52 path. The carry helpers come from `monty.h`.
* `getfactor_tecm_edw_x8_list` (8 inputs, one per lane): table curves,
  then GKL m = 3, 4, ... from the constant [3]Q. With target_bits = 0 the
  pre-passes run first, and only unfactored inputs go on to the main pass
  (the old code reran every input).
* `getfactor_tecm_edw_x8` (one n, 8 curves per round): lane i of round r
  runs curve 8r + i (table, then GKL m = id - 2).
  - Each lane starts its [m]Q lazily from a running [k]Q and steps it by
    [8]Q (a doubling for the lane at m = 8).
  - The modulus is set up once with mpz and copied to all lanes.
  - The first version did 8 mpz setups and an eager 14-step start; it was
    0.74-0.95x at small t.
* Stage 2 keeps tinyecm x8's rule: a product that would make acc 0 mod n
  freezes that lane for the rest of stage 2.

## G2. Verification

`arith/edwards/tedw_test.c` (driver) and `uedw_oracle.py`. The oracle's
stage-2 replay now has three rules (scalar break, microecm x8 skip,
tinyecm x8 freeze) and reads either file's tables.

| check | what | result |
|---|---|---|
| `tcurves` | scalar, one curve per (n = p q, id), n up to 128 bits | 0 mismatches in 2100 + 1050 + 13999 (`tedw_oracle_run.txt`) |
| `tx8` | x8 list: first finding curve and prime per input, total curves run | 0 of 2100 + 14000 on IFMA and AVX-512F; curve totals equal |
| `tone` | x8 single n: first (round, lane) that finds, and the prime | 0 of 420 + 2800 on IFMA and AVX-512F (also after the lazy-start rewrite) |

## G3. Speed

`tiny_e2e.c -DWITH_EDW`, 8000 inputs per row (the single-n rows: 2000).
Inputs: p of t bits times q of 2t bits, target t. us per input, old ->
new (speedup):

| t | B1 | scalar (IFMA build) | x8 list, IFMA | x8 single n, IFMA | x8 list, AVX-512F | x8 single n, AVX-512F |
|---:|---:|---|---|---|---|---|
| 20 | 70 / 27 | 26.3 -> 24.6 (1.07) | 6.56 -> 5.76 (1.14) | 13.5 -> 13.5 (1.00) | 9.05 -> 8.97 (1.01) | 18.9 -> 19.1 (0.99) |
| 22 | 70 / 47 | 40.3 -> 35.8 (1.13) | 10.3 -> 7.98 (1.29) | 20.4 -> 17.2 (1.19) | 14.2 -> 12.0 (1.18) | 28.3 -> 26.1 (1.08) |
| 24 | 70 | 67.1 -> 58.6 (1.15) | 14.9 -> 12.6 (1.18) | 27.3 -> 23.4 (1.17) | 21.3 -> 18.8 (1.13) | 39.1 -> 35.9 (1.09) |
| 26 | 85 | 104 -> 93.8 (1.11) | 23.3 -> 17.9 (1.30) | 36.3 -> 30.9 (1.18) | 33.3 -> 27.7 (1.20) | 52.2 -> 46.7 (1.12) |
| 29 | 125 | 209 -> 183 (1.14) | 44.1 -> 32.8 (1.34) | 61.9 -> 48.4 (1.28) | 64.0 -> 52.0 (1.23) | 88.8 -> 75.4 (1.18) |
| 31 | 165 | 330 -> 305 (1.08) | 69.0 -> 53.2 (1.30) | 93.0 -> 72.6 (1.28) | 101 -> 82.8 (1.22) | 135 -> 112 (1.20) |
| 32 | 205 | 404 -> 362 (1.12) | 85.6 -> 62.4 (1.37) | 112 -> 86.6 (1.29) | 126 -> 99.7 (1.27) | 162 -> 135 (1.20) |

* The B1 column gives the scalar level / x8 level where they differ; the
  dispatches use different thresholds.
* All against the Part D (already optimized) tinyecm.
* Find rates are equal or higher everywhere (e.g. t = 32 x8 list 0.9825 ->
  0.9906).
* The scalar build ran the same code on both ISAs (1.07-1.15x on the
  AVX-512F build too).
* As Phase 0 predicted, the scalar gain is mostly yield: the Edwards stage
  1 is about even with PRAC at 128 bits.
* The x8 versions gain more, from Edwards' ILP with the fast IFMA
  multiply.
* The smallest level (t = 20) is dominated by the per-input setup; it is
  about even.

## G4. Next

* Phase 4: AVX2 x4. Add an FP52 arithmetic type to edw_kern.h: 52-bit for
  microecm, 104-bit for tinyecm.
* Possibly move microecm_edw onto edw_kern.h too (it has the same code
  written out).
* Then the B1 / B2 study (F5).

---

# Part H: Phase 4a, microecm_edw on AVX2 x4 (2026-11-09)

`factor/gmp-ecm/microecm_edw_x4.c`: the x8 list's algorithm (same curves,
chains, B2 = 25 B1 stage 2 with the "skip" accumulator, same results) on
AVX2, 4 inputs n < 2^52 per ymm, using the shared kernel `edw_kern.h`.
`getfactor_uecm_edw_x8_list` calls it on AVX2-only builds;
`getfactor_uecm_edw_x4_list` / `uecm_edw_x4_run` are the direct entry
points.  Built with `-mfma` (Makefile `EDW_X4_OBJS`); no MXCSR change.

## H1. Arithmetic

Residues are doubles holding integers in [0, n), R = 2^52.  The product
uses FMA with the default round-to-nearest (variant (b) of B0):
* `t = fma(x, y, 2^104)`, `hd = t - 2^104`, `l1 = fma(x, y, -hd)`: x y =
  hd + l1 exactly, hd = H1 2^52, |l1| <= 2^51 (signed low half).  The same
  split with bias 3 2^103 handles signed products.
* `m = l1 rho mod 2^52` (rho = -1/n), signed, |m| <= 2^51; `v = fma(m, n,
  3 2^103)` gives H2.
* When |l1| < 2^51 the low half of m n is exactly -l1, so r = H1 + H2
  straight from the high halves: `r = fma(v, 2^-52, H1 - 3 2^51)`, with
  r - n and r + n formed by two more FMAs in parallel and two sign-bit
  blends to land in [0, n) (r is in [-n/2, 3n/2)).
* A tie |l1| = 2^51 (probability ~2^-52 per lane) needs the carry
  (l1 + l2) / 2^52: a predictable branch takes the full exact form
  (`uedw4_mulredc_full`).  The check costs within noise (0-3%).
* add / sub: one FP add, one FP sub, one sign blend.

The first version (full form always) had a ~58-cycle dependent latency
per product, against ~12 for scalar mulx code.  The Edwards formulas offer
only 2-4 independent products per step, so it was latency-bound and lost
to scalar microecm_edw (0.77-0.89x).  The direct-from-high-halves form cut
the latency to ~40 cycles: microbenchmark 0.52 ns per lane-product with 8
independent chains (scalar 0.77 ns), 1.15 ns with 2 chains.

Interleaving two ymm per residue (`UEDW4_W = 2`, 8 lanes) was slower (register
spills), so W = 1 is the default; the option is kept.

## H2. Verification

| check | result |
|---|---|
| `uedw_oracle.py x8` on the x4 build (`uedw_test.c -DUEDW_TEST_X4`), AVX2-only | 0 of 21000 + 21000 mismatches (1-2 Edwards-exception early finds, allowed) |
| same, AVX-512 build (gcc may use VL forms) | 0 of 21000 |
| `u4_arith.c`: mul / add / sub vs __int128, n near 2^52 and small, x, y near n, forced ties | 0 of 16M lanes bad (1.08M ties) |
| `objdump` of the AVX2-only build | no zmm |

Run log: `arith/edwards/uedw4_oracle_run.txt`.

## H3. Speed (i5-11600K, one core, `x4_build.sh`, `tiny_e2e` 4000 inputs, two seeds)

AVX2-only build (no x8 there: the old `getfactor_uecm_x8_list` falls back
to scalar microecm):

| bits | B1 | old microecm us | microecm_edw us | x4 list us | x4 vs old | x4 vs edw scalar |
|---|---|---|---|---|---|---|
| 40 | 27 | 2.70-2.83 | 2.27-2.32 | 2.18-2.19 | 1.23-1.30x | 1.04-1.07x |
| 44 | 47 | 4.14-4.21 | 3.37 | 3.11 | 1.33-1.35x | 1.08x |
| 48 | 70 | 6.37-6.39 | 5.03-5.24 | 4.54-4.55 | 1.40x | 1.11-1.15x |
| 52 | 85 | 10.0-10.4 | 8.19-8.31 | 7.01-7.33 | 1.37-1.48x | 1.12-1.19x |

Found rates are 1.0000 throughout.  On AVX-512 the x8 list is 3-5x
faster than x4, as expected; x4 there is only a fallback.

## H4. Next

* Phase 4b: tinyecm_edw x4 (104-bit as two 52-bit FP digits).  Its
  product has 11 partial products, so latency will matter more; the
  high-half shortcut needs care with the carries between digits.
* Possible x4 refinements: precompute y rho for constant multipliers
  (cuts m's latency); a cheaper add/sub.
* Then the B1 / B2 study (F5).
