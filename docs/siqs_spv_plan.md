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

## Results (Commit 4 sweep — measured Oct 2026, Xeon Gold 6254, AVX512BW, 1 thread)

The sweep harness lives in `docs/spv_sweep.sh`; re-run it on other hardware to
check whether these conclusions hold (it rediscovers the formula
`tf_small_cutoff` per (input, limit) itself):

    make siqs USER_CFLAGS=-DSIQS_PROFILE     # fills the stage-timer columns
    cd <scratch dir>
    docs/spv_sweep.sh <path>/siqs_demo input.txt "256 512 1024" "-10 -5 0 5 10"

Setup: c70 (232-bit) and c80 (265-bit, DLP) semiprimes; fresh savefiles every run;
`-noopt -siqsSPV <L> -siqsTFSm <formula+offset>`; formula `tf_small_cutoff` measured
per config (c70: 19/20/21 at L=256/512/1024; c80: 21/22/23). Single runs — treat
deltas below ~5% as noise. Stage timers via `-DSIQS_PROFILE` builds.

c70 (non-DLP), total seconds / t_tdiv:

| L   | -10       | -5        | 0 (formula) | +5        | +10       |
|-----|-----------|-----------|-------------|-----------|-----------|
| 256 | 13.1/5.56 | 12.5/4.98 | 13.0/4.47   | 12.5/4.28 | 14.3/3.98 |
| 512 | 13.0/5.27 | 12.5/4.72 | 12.9/4.36   | 13.2/4.09 | 13.6/4.03 |
| 1024| 13.6/5.34 | 13.3/4.85 | 13.9/4.68   | 14.1/4.56 | 13.6/4.13 |

c80 (DLP), total seconds / t_tdiv:

| L   | -10       | -5        | 0 (formula) | +5        | +10       |
|-----|-----------|-----------|-------------|-----------|-----------|
| 256 | 77.7/40.2 | 74.0/37.0 | 70.4/32.3   | 80.4/34.7 | 78.0/31.1 |
| 512 | 75.5/38.3 | 72.4/35.0 | 70.5/31.8   | 73.6/30.9 | 75.8/30.8 |
| 1024| 76.5/36.1 | 76.7/35.3 | 72.1/32.3   | 73.3/31.1 | 71.8/29.4 |

Findings:

1. **Keep 256.** Raising the limit is neutral-to-negative on total time at both
   sizes. The vectorized filter cuts t_tdiv (~4-5% at L256), which offsets the
   extra SPV work as the limit rises (tdiv flat, sieve load unchanged), but the
   freed time does not translate into faster factorizations. Relation yield
   rises with the limit (c80: 1490/1599/1839 rels at formula) without helping
   total time.
2. **TFSm formula value is optimal at both sizes.** -10 sharply increases tdiv
   (more surviving candidates); +5/+10 shift work into the sieve without
   finishing sooner. This held for every limit tried.
3. **DLP co-tuning rejected.** Subtracting tf_small_cutoff from closnuf in the
   DLP branch (mirroring the non-DLP path) measured ~2x total time on c80
   (70.4 -> 139.0 s at L256/TFSm=21; t_tdiv 32 -> 106 s). The one-line change
   is left in the code, commented out, with the measurement noted.
4. No size-dependent default is adopted; `-siqsSPV` remains an override knob
   (default 256 = stock behavior) and the Commit-3 tuner (two-stage, median)
   picks tf_small_cutoff per run.

Bugs found and fixed during the sweep (all guarded):

- At raised limits the fb index of the first prime > limit approaches the
  10-bit boundary and can land unaligned; the stock realignment loop in
  tdiv_medprimes_32k_avx2 then cannot complete before bound10 and the SIMD
  groups run misaligned (segfault). siqs_static_init now rounds the SPV
  boundary down to a 16-aligned index in exactly that case and adjusts
  small_limit so the unsieved range and tinylist stay consistent (no change
  at the default limit).
- siqs_start did not initialize the Commit-3 tuner fields of the
  stack-allocated siqs_userdata_t (num_samples), crashing yafu (caught by the
  yafu smoke test, not by siqs_demo where the stack happened to be zeroish).
