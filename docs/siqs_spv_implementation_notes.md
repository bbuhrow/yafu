# SIQS SPV Experiment — Implementation Notes

Post-mortem/context for future sessions working on the code introduced by the
four SPV commits (see `docs/siqs_spv_plan.md` for the plan and results, and
`docs/spv_sweep.sh` for the measurement harness).  Focus: what bit us, why,
and how it was fixed — so the same class of bug isn't rediscovered.

## 0. The guard, and how to keep it honest

- Everything inside the qs library sits behind `#define TRY_SMALLP_OPT 1`
  at the top of `include/qs_impl.h`, wrapped `#ifdef / #else <original
  verbatim> / #endif`.  Driver-side plumbing (`factor.h`,
  `factor_common.c`, the three `cmdOptions` copies, `driver.c`,
  `siqs_demo.c`, `avx-ecm/demo/main.c`) is intentionally unguarded; with the
  define off, `-siqsSPV` parses and is ignored.
- Every commit was verified in BOTH define states before staging.  The
  define-off build caught one real breakage (see 3.2) and the define-on
  state was exercised by actual factorizations, not just compilation.
- When adding any later qs-library edit during follow-up work, wrap it the
  same way; several of the gotchas below came from edits that straddled a
  preprocessor boundary.

## 1. Repo hygiene gotchas (cost real time)

### 1.1 Mixed line-ending styles
- The repo has **no `.gitattributes` and `core.autocrlf` unset**, and files
  are committed inconsistently: qs-library files (`qs_impl.h`, `SIQS.c`,
  `siqs_aux.c`, `tdiv_small.c`, `factor.h`, `factor_common.c`,
  `qs_demo/*`, `driver.c`) are **CRLF in HEAD**, while `top/cmdParser/*`
  and `avx-ecm/demo/*` are **LF in HEAD**.  The working tree was CRLF
  almost everywhere.
- Editing a file whose endings don't match your editor's insertions fails
  exact-match edits (the `edit` tool matches `\n`, CRLF files don't contain
  it) — symptom: "Could not find oldString" on text you can plainly see.
- The fix used: convert the target file to LF (`sed -i 's/\r$//'`), edit,
  then convert back to the file's **HEAD-native** ending before
  `git add` (`sed -i 's/$/\r/'` for the CRLF-in-HEAD files).  Committing
  the wrong direction makes the diff a whole-file rewrite (4000+ lines of
  noise on `driver.c` alone).

### 1.2 Three cmdOptions copies, not two
- The plan said "both copies" (avx-ecm/demo, qs_demo), but the **Makefile
  links `top/cmdParser/cmdOptions.c` into yafu** — a third, fuller copy
  (NUMOPTIONS 133 vs 114/31).  Always verify which copy actually links
  before assuming the driver plumbing is done; `driver.c`'s
  `options_to_factobj` is yafu's, `siqs_demo.c` and `avx-ecm/demo/main.c`
  have their own.
- Each copy's parallel arrays (`OptionArray`, `OptionHelp`, `needsArg`,
  `LongOptionAliases`) must all grow together; the avx-ecm aliases array
  was found miscounted after editing and had to be regenerated to exactly
  one entry per option.  Count-check scripts (python over the array
  blocks) caught it; do that before building.

### 1.3 Build environment
- This machine has no `../gmp-install/6.2.0-gcc` despite `config.mk`
  pointing there; builds need `GMP_PREFIX=/opt/opencode/tmp/795601290/opencode/gmp`.
- **Changing compile flags does not trigger rebuilds**: adding
  `USER_CFLAGS=-DSIQS_PROFILE` to an up-to-date tree silently rebuilds
  nothing.  `touch` the sources that consume the flag (or `make clean`)
  or you'll measure a binary without the instrumentation you think you
  added.

## 2. Commit 1 (knob + allocation) gotchas

- The `tinylist` was sized from a count over `fobj->primes` (all primes to
  10^7, present in every `fact_obj_t`), plus slack for the dummy index-0
  entry.  The consistency check (start < alloc) is a plain check-and-exit,
  **not** `assert()` — the release build defines `NDEBUG`, which would
  compile real asserts out.
- `siqs_static_init` computes `tiny_alloc` inside the FB-construction
  `while(1)` loop but checks it after the loop; the variable therefore had
  to be declared at function top (guarded) to stay in scope.

## 3. Commit 2 (vectorized filter) gotchas

### 3.1 Lemire test details worth remembering
- Constants: `invmod = p^-1 mod 2^32` via Newton iteration starting from
  `x = p` (exact mod 8 for odd p; 5 doublings reach 2^32), `thresh =
  0xffffffffu / p`.  Hit ⟺ `(offset + p − root) * invmod <= thresh`
  (mod 2^32).  The `+ p` form keeps the intermediate non-negative and
  covers the `root == 0` / `root == prime` corner cases the scalar loop
  special-cased — the assert mode confirmed exact hit-set equality, so the
  corner-case reasoning held.
- `prime == 2` is FB index 1 and never reaches this code (powers-of-2 loop
  handles it); the vector path runs on odd primes only.  Guard on
  `USE_AVX512BW`, not `USE_AVX512F`: the masked u16 tail loads
  (`_mm256_maskz_loadu_epi16`) need AVX512BW+VL.  AVX512F-only builds fall
  back to the preserved scalar path.
- Tail groups must use **masked loads for every array** (fbc prime/roots
  AND tinylist invmod/thresh): for small inputs `B` can be barely larger
  than `sieve_small_fb_start`, so an unmasked 256-bit load of 16 u16 runs
  off the end of the array even though the *used* lanes would be fine.
- Correctness check design: don't compare `bits` or Q values (destructive,
  order-dependent); compare the **hit index sets** (vectorized hits are
  processed in ascending index order, so a memcmp works).  That found zero
  mismatches on c70+c80 full factorizations and is cheap enough to enable
  for a whole run (`SPV_ASSERT_BOTH` in tdiv_small.c).

### 3.2 The #ifdef/#else wrap that ate a function
- Rewriting `siqs_sync`'s tuner block by slicing the original at the
  `"// free sieving structure"` marker cut the block **midway** — that
  comment sits *inside* `if (t[tid].dconf->buffered_rels)`, before the
  counter resets and the block's closing brace.  The result compiled with
  the define ON (the stray remainder happened to be syntactically
  absorbable) but failed define-off with confusing errors far away
  (`invalid storage class for function 'now_sec'` = a function definition
  swallowed by an unclosed brace).
- Lessons: (a) never slice a guarded region at a comment marker — verify
  the extracted region is brace-balanced on its own; (b) an unbalanced
  #else branch surfaces as errors at the *next function definition*, not
  at the edit site; (c) checking `cpp -fpreprocessed` brace balance of the
  affected function in both define states finds this in seconds.
- The final structure guards only the tuner `if`-block; the
  `if (buffered_rels)` opener, merge lines, counter resets and closing
  brace stay shared outside the guard (identical in both paths).

## 4. Commit 3 (tuner) gotchas

### 4.1 Stack-allocated state, initialized in a different function
- `siqs_userdata_t udata` is a **stack local of `SIQS()`**; its fields are
  initialized in `siqs_start()`, which only touched the pre-existing ones.
  The new `num_samples` was left as stack garbage → in yafu it was huge →
  `rels_samples[num_samples++]` overflow → SIGSEGV.  siqs_demo didn't
  crash because its stack happened to be zeroish — **the standalone binary
  is not a sufficient smoke test for shared-state lifetime bugs; run the
  same input through yafu too.**
- Fix: guarded init of all new udata fields in `siqs_start` next to the
  existing ones.

## 5. Commit 4 (sweep) gotchas

### 5.1 FB index vs prime value (the alignment segfault)
- The factor base holds only primes with a root of `n mod p` (~half of all
  primes, plus multiplier adjustments), so **FB indices are ~2× denser in
  prime value than natural prime indices**: index ~30 ≈ prime 251, index
  ~87 ≈ prime ~1024.  Any reasoning of the form "limit 1024 ⇒ start index
  π(1024) = 172" is wrong by ~2×.
- With limit 1024, `sieve_small_fb_start` (~87) lands **within 8 of
  `fb_10bit_B`** and unaligned.  `tdiv_medprimes_32k_avx2`'s realignment
  loop (scalar-test until `i & 15 == 0`) can then exit at `bound10`
  itself, and the SIMD groups (which require 16-byte = 8-element alignment
  of the compressed fb arrays, because the `small_inv` shift changes on
  8-aligned index boundaries) run misaligned → `vmovdqa` segfault inside
  `MOD_CMP_8X_vec`.  The stock code only ever worked because at limit 256
  the boundary sits ~58 indices below the 10-bit boundary.
- Fix (guarded, in `siqs_static_init`): only when the realignment cannot
  complete (next 16-aligned index would exceed the would-be `bound10`),
  round `sieve_small_fb_start` **down** to 16-aligned and adjust
  `small_limit = prime[start-1] + 1`.  Rounding down (not up) matters:
  rounding up would put primes > limit into the unsieved range where
  `make_fb_siqs` (fill condition `prime < small_limit`, unguarded original)
  has no tinylist entries.  Adjusting `small_limit` keeps
  `prime < small_limit ⟺ index < start` exactly true, which also closes
  the pre-existing edge case `limit == prime` (unsieved prime with no
  tinylist entry).  At the default limit the trigger condition is never
  met → zero behavior change (verified: start 30, limit 256 unchanged).
- The same crash reproduced under c70 but not c80 at the same limit —
  different `n` ⇒ different jacobi subset ⇒ different boundary alignment.
  A single passing config does not clear an alignment hazard.

### 5.2 Measurement harness gotchas (encoded into docs/spv_sweep.sh)
- `siqs_demo` requires a positional argument (`numRequired = 1`) and the
  `@`-substitution form: `echo N | siqs_demo "siqs(@)" -flag value ...`.
  Options must be **separate argv entries** — embedding them in the
  expression string makes the calc parser eat them
  (`unrecognized variable or function 'siqsSPV256'`).
- The "trial factoring cutoff at N bits" splash line prints `tf_closnuf`,
  **not** `tf_small_cutoff`; the formula value needs the guarded
  `SPV formula tf_small_cutoff = N` diagnostic (added, prints at
  `-v -v`).  TFSm offsets are relative to the formula value, which shifts
  ~1 bit per limit doubling.
- Stage timers (`Tdiv time` etc.) print only under `SIQS_PROFILE` (never
  defined by the Makefile; pass `USER_CFLAGS=-DSIQS_PROFILE` and touch
  sources) **and** at `VFLAG > 0`.
- Noise discipline: single runs vary ±3-5%; conclusions were drawn only
  from deltas well outside that (the DLP variant's 2× blowup) or from
  consistent direction across the whole matrix (TFSm formula optimal).

### 5.3 Measurement outcomes (see plan doc for tables)
- Raising the limit: t_tdiv flat, relation yield up, total time flat —
  the vectorized filter pays for the extra SPV work but buys nothing.
- DLP co-tuning (subtracting `tf_small_cutoff` from `closnuf` in the DLP
  branch): ~2× total time (t_tdiv 32→106 s) — the lowered report
  threshold floods trial division.  Line left in code, commented out.
- Mid-run `small_limit` switching is structurally impossible (it defines
  the fb/sieve layout); the "tuner tries two limits" idea from the plan
  was dropped for this reason, with the sweep answering the same question.

## 6. Verification summary (what was actually run)

- Define-off: clean build of `siqs_demo` + `yafu`; c70 factorization with
  identical timing; tuner reverts to the original 3-measurement/mean path.
- Define-on: `-siqsSPV 256` matches baseline; 512/1024 factor correctly in
  both binaries (including through the alignment fix); assert-both-paths
  clean on full c70+c80 runs; `SPV_ASSERT_BOTH` left commented off.
- `OPT_DEBUG` remains commented in committed code; `-DSIQS_PROFILE` was
  used only in throwaway measurement builds and the released binaries are
  plain `-O2` release builds.
