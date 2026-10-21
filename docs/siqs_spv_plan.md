# SIQS Small-Prime Variation (SPV) Tuning Plan

**Repo:** yafu.git · **Target:** AVX512BW · **Job sizes:** c60–c90 · **Status:** Approved plan — implementation

## Objective

Determine empirically whether raising `small_limit` (the number of FB primes not sieved,
currently 256, set at `factor/qs/SIQS.c:4196`) improves total time-to-factorization for
c60–c90 inputs on AVX512BW, vectorize the small-prime divisor check in `filter_SPV` to
make higher limits viable, and improve the online tuner. Adopt code changes only where
measurements justify them; "keep 256" is an acceptable outcome.

## Governing convention: `TRY_SMALLP_OPT` guard

- Add `#define TRY_SMALLP_OPT 1` near the top of `include/qs_impl.h`. Commenting out
  that one line reverts the qs library to pre-experiment behavior with **zero code
  deletion**.
- Every experimental edit inside the qs library is wrapped
  `#ifdef TRY_SMALLP_OPT / <new code> / #else / <original code preserved verbatim> / #endif`.
- **Guarded (qs purview):**
  - `include/qs_impl.h`: the define; new `tiny_fb_element_siqs` members (`invmod`,
    `thresh`); any new instrumentation counters in `static_conf_t`/`dynamic_conf_t`.
    Struct layout byte-identical when define is off.
  - `factor/qs/siqs_aux.c`: `invmod`/`thresh` computation in `make_fb_siqs()`
    (adjacent to lines 209–223, 286–300).
  - `factor/qs/tdiv_small.c`: vectorized SPV loop replacing lines 216–308 (original
    batch-of-4 + tail loops preserved in `#else`).
  - `factor/qs/SIQS.c`: small_limit override application in `siqs_static_init`; sized
    `tinylist` allocation + start-index assert (vs fixed 512 at lines 4499–4506); all
    `siqs_sync` tuner changes (lines 182–289); chosen-value print; any Phase E
    default/DLP-branch `closnuf` change.
- **Unguarded (outside qs purview):** `factor.h` / `factor_common.c`
  (`gbl_override_spv_limit` field + init + copy), `cmdOptions.c/.h` parse + help
  (both copies: `factor/avx-ecm/demo/` and `factor/qs/qs_demo/`). With the define off,
  `-siqsSPV` is parsed but ignored.
- **Untouched:** `td_small_p` SS-search path (SIQS.c:1467–1553), GPU (`siqs_gpu*.cu`),
  `tinySIQS.c`, `smallmpqs.c`.

## Background (verified code facts)

- `sieve_small_fb_start` = first FB index with prime > `small_limit` (SIQS.c:4641–4646);
  all consumers key off it and shift together when the limit changes.
- `make_fb_siqs()` fills `tinylist` only for `prime < small_limit` (siqs_aux.c:209–223,
  286–300); primes ≥ limit get `inv32/corr32/shift32` precompute (siqs_aux.c:258).
- `tf_small_cutoff = avg + 2.5·sd`, `sd = √28` magic (SIQS.c:5328–5345); `avg` = mean
  per-location log contribution of unsieved primes. Non-DLP path subtracts it from
  `closnuf` (SIQS.c:5377); DLP path (`use_dlp==1`, SIQS.c:5351–5368) uses a digits-based
  `closnuf` and does **not** subtract — key asymmetry for c85/c90.
- `filter_SPV` (tdiv_small.c) divides out unsieved primes per candidate, accumulating
  `logp` into `bits`; abandons candidate if `bits < tf_closnuf + tf_small_cutoff`
  (line 316).
- Online tuner (`siqs_sync`, SIQS.c:182–289): 3 measurements (`orig`, `+5`, `−5`),
  `num_avg` polys each (10/5/2 by bits), 2% (or 5% DLP) superiority bias. Disabled by
  `-noopt`, `-siqsTFSm` (SIQS.c:158–164), and restart (SIQS.c:3157).
- **Constraint:** `tinylist` is a fixed 512-entry allocation (SIQS.c:4499–4506).
  π(2048)=309 fits; π(4096)=564 silently overflows — latent bug, must be fixed before
  limits ≥ 4096.
- Compressed roots are `uint16_t` (qs_impl.h:305–315) — safe to 65535. Divisionless
  32-bit mod exact for offsets < 2³²; actual offsets < ~1.3M.
- Author's note (SIQS.c:4194–4195): raising the limit requires co-raising TFSm /
  lowering TF bounds; prior experiments found ~neutral gains at 256.

## Phase B — `-siqsSPV` knob + allocation fix (guarded) — Commit 1

1. New CLI/.ini option `-siqsSPV <limit>`, plumbed exactly like `siqsTFSm`
   (see unguarded files above); consumed in `siqs_static_init` **before** the FB loop at
   SIQS.c:4475. Default 256 = zero behavior change.
2. Sized `tinylist` allocation (derived from limit) + assert
   `sieve_small_fb_start < alloc`; update stale comment at SIQS.c:4194–4196 (guarded).
3. **Acceptance:** splash prints expected SPV bound (SIQS.c:4911–4932);
   `-siqsSPV 256` matches baseline relations/timing; define-off build compiles clean
   and is behaviorally identical to pre-experiment.

## Phase C — Vectorized small-prime divisor check (guarded) — Commit 2

Lemire-style divisibility test, 16 primes per zmm (AVX512F 32-bit lanes; all AVX512BW
CPUs have them):

- New per-prime constants in `make_fb_siqs()`: `invmod[i] = p⁻¹ mod 2³²` (Newton
  iteration; p odd — p=2 already handled by the powers-of-two loop at tdiv_small.c:193),
  `thresh[i] = (2³²−1)/p`.
- Hit test per root: `p | (offset − root)` ⟺ `(offset + p − root) · invmod ≡ ≤ thresh
  (mod 2³²)` — one `mullo` per root, no 64-bit splitting. The `p−root` form naturally
  covers the `root==prime`/`root==0` corner cases the current code special-cases
  (tdiv_small.c:299–303).
- Per batch of 16 primes: load 16×u16 `fbc->prime/root1/root2`, widen (`vpmovzxwd`),
  `t = offset_bcast + pvec − rootvec`, `mullo` by `invmod`, `cmpule` vs `thresh`, OR
  the two masks, iterate hits with `tzcnt` into the existing scalar `DIVIDE_ONE_PRIME`
  (mpz divide + repeat-while-divisible + `bits += logp`). Hits are rare (0–2/candidate);
  scalar division stays.
- Replace both loops (batch-of-4 and tail, tdiv_small.c:216–308) with one masked loop
  over `[2, sieve_small_fb_start)`; keep existing `small_inv`/`correction` (still used
  by poly_roots_32k.c:2044–2054); `#else` keeps the scalar path; add temporary
  assert-both-paths mode.
- ~10–12 vector uops per 16 primes vs ~30+ scalar uops per 4 primes (~4–6× fewer).
  Interleave two 16-prime groups to hide `vpmulld` latency (~10 cyc).
- **Acceptance:** identical `valid_Qs`/`Qvals` vs scalar path on fixed c60/c80 inputs;
  stage timers show `t_tdiv` reduction.

## Phase D — Sweep (vectorized filter on) — Commit 4

1. Per size point: `small_limit ∈ {256, 512, 1024}` × `TFSm ∈ {−10, −5, 0, +5, +10}`
   relative to each run's formula value (`avg` shifts with the limit), using
   `-noopt -siqsTFSm`.
2. DLP regime (c85, c90): test the co-tuning hypothesis — raised limit may need
   explicit TFSm co-tuning, or a guarded one-line change to subtract
   `tf_small_cutoff` in the DLP `closnuf` branch too. Measure both.
3. Decide on **total time**, not rels/sec; watch `t_tdiv` and filter time for
   SPARSE_STORE work migration (filter.c:289).

## Phase E — Tuner + defaults (guarded) — Commit 3

1. Two-stage ±5→±2 refinement; median-of-polys; optional extra early round trying
   limit 256 vs Phase D winner — **gated to bits ≥ 230, off by default behind a flag**.
2. Print chosen `tf_small_cutoff` each run (guarded).
3. Adopt size-dependent defaults (or DLP-branch change) only on repeated, clear
   Phase D wins; "keep 256" is an acceptable outcome.

## Verification checklist

- Define off: bit-identical relations and timing vs pre-experiment build on all 5
  benchmarks.
- Define on, defaults (`-siqsSPV 256`): same relation sets as baseline within tuner
  noise.
- Vectorized vs scalar assert mode: identical per-report decisions.
- All 5 benchmark factorizations produce correct factors.
- Both `cmdOptions` copies updated consistently; verify which links into `yafu` via the
  Makefile; smoke-test `-siqsSPV` in both `yafu` and the standalone `siqs` binary.
- OPT_DEBUG (SIQS.c:56) stays commented in committed code.

## Risks

- `tinylist` overflow at limit > ~3600 (fixed in Commit 1).
- DLP threshold asymmetry silently cutting relation yield at c85+ without co-tuning
  (Commit 4.2).
- `vpmulld` latency (mitigated by interleaving two 16-prime groups).
- Guard drift: any later qs-library edit must also be wrapped — enforced at every commit.

## Commit sequencing

1. Guard scaffold + knob + alloc fix
2. Vectorized filter
3. Tuner changes
4. Sweep results → default decision

Each commit keeps define-off builds clean, enabling bisect/revert at any boundary.

## Build notes

- Use the repo Makefile; targets `make yafu` and `make siqs`.
- Timed runs must use fresh savefiles — a savefile restart silently disables the tuner.
- OPT_DEBUG traces are for local experiments only, not committed enabled.
