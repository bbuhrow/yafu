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
