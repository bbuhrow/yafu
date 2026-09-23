# collharness -- Phase 0 test harness for the `gpu_gerbicz` OpenCL port

Portable (gcc / clang / MSVC), plain C99, zero GMP/CUDA/OpenCL
dependencies. Lets later phases (the OpenCL collision engine, sieve
kernels, host integration) be validated bit-for-bit without an NVIDIA
GPU.

## Building

### gcc / clang (Linux, macOS, MinGW)

```
make                # uses cc by default
make CC=clang
```

Produces `./collharness`. Built and tested clean under both gcc 13 and
clang 18 with `-std=c99 -Wall -Wextra -pedantic`.

### MSVC (Windows, `cl.exe`)

No Makefile is provided for MSVC; a plain command line works because
nothing here uses VLAs, `__int128`, or any GNU extension:

```
cl /std:c11 /W4 /I include src\main.c src\collcase_io.c src\gen.c ^
   src\cpuref_exact.c src\cpuref_stats.c src\cmp.c /Fe:collharness.exe
```

(`/std:c11` rather than a C99 flag because MSVC has no dedicated C99
switch and its C11 mode is the closest strict superset it offers; the
source itself is C99 -- it was not exercised under MSVC in this
sandbox, since no Windows toolchain was available here. If a future
phase touches this harness, verifying an actual MSVC build is a
worthwhile trust-but-verify step.)

## CLI

```
collharness gen --out FILE [options]     generate a synthetic collcase
collharness ref --in FILE --out FILE     (re)compute expected output
collharness cmp OBSERVED REFERENCE       compare two result files
collharness suite --out-dir DIR [--large] generate the named suite
collharness info FILE                    print a file's header
```

`gen --help`-equivalent options (see also `collharness --help`):

| option | default | meaning |
|---|---|---|
| `--n N` | 10000 | element count |
| `--root-bytes 4\|8` | 4 | key width |
| `--key-bits K` | 24 | signed key range is +/-2^(K-1) |
| `--shift S` | 20 | bits allocated to `p` in each value |
| `--bucket-hash 0\|1` | 0 | 0 = mask hash, 1 = multiplicative |
| `--hash-word-cap W` | 0 (= 4096) | filter hash-table size override |
| `--seed S` | 1 | PRNG seed (splitmix64 -- fully deterministic) |
| `--density D` | 0.01 | filler collision density, in [0,1] |
| `--skew N` | 0 | force N extra distinct keys into one bucket |
| `--no-edge-cases` | | omit the fixed edge-case suite |
| `--no-expected` | | skip computing the expected-output section |

## Format

See `include/collcase.h` for the full, byte-precise collcase v1
layout (header fields, array layout, optional expected-output
section with `found_t`/`specialq_t`-shaped entries and the 424-byte
stats block). Two points worth calling out:

- All integers are little-endian and are read/written field-by-field
  (`include/byteio.h`), never as a raw `fwrite` of a struct -- so the
  format is identical regardless of host endianness or compiler
  padding/alignment choices.
- Keys are stored as the **raw two's-complement bit pattern** in
  their `root_bytes`-wide slot, not sign-extended to 64 bits on disk.
  In memory (`collcase_t.keys`, always `uint64_t`), a 4-byte key is
  zero-extended (matching what the real engine's `scatter_roots_kernel`
  does), and sign-extension back to a signed offset only happens when
  an entry is emitted -- see `sign_extend_key()` in
  `src/cpuref_exact.c`.

## What's modeled, and one deliberate shortcut

`src/cpuref_exact.c` is the ground-truth found-set (0.4): sort by raw
key, group, and for every pair in a group with matching `q_index` and
`gcd(p1,p2)==1`, emit. Independent of any hash filter, so this is what
every engine (CUDA, OpenCL, this harness) must agree with.

`src/cpuref_stats.c` (0.5) faithfully replays
`filter_per_bucket_kernel`'s multi-round hash-collision filter
round-for-round (same shift schedule, same ping-pong tables, same
convergence/cap/zero stop conditions, same 102-slot histogram
layout) to reproduce `candidate_count` and `filter_iters_hist`
bit-for-bit. This is a genuine simulation, not a shortcut, because the
filter is a probabilistic/approximate Bloom-style construction with
real false positives.

`dedup_count` and `value_match_count`, however, are **not** re-derived
by also simulating the engine's D/S/X secondary hash tables. Those
tables are an exact, collision-free lookup (no false positives) --
an implementation detail for GPU parallelism, not a source of
approximation -- so as long as no capacity/overflow path is hit
(assumed throughout; see the plan's "Open issues"), the secondary
hash's output is mathematically identical to grouping the filter's
candidate output by exact value and counting true original
multiplicities directly, which is what `cpuref_stats()` does. This
is documented in `include/cpuref.h` and was cross-checked, not just
asserted: every generated case in this delivery has
`dedup_count == (# distinct duplicated keys)` and
`value_match_count == sum of their true multiplicities` verified by
hand against the deterministic edge-case suite (see STATUS.md).

## Verification performed

- `tiny_basic4`/`tiny_basic8`/`tiny_hashmode1` (edge cases only,
  n=64): found_count/dedup_count/value_match_count hand-verified
  against the 0.3 edge-case design (50 / 8 / 24 in every case,
  independent of root_bytes or bucket_hash, as expected since none of
  those change the *key values themselves*, only their encoding or
  bucketing).
- Round-trip: for every suite case, stripping the embedded expected
  output and recomputing it via `collharness ref` reproduces it
  exactly (`collharness cmp` reports MATCH), including one case
  (`med_hashmode1`) that exercises the saturation/subset-check branch
  of the comparator (found_count > 999).
- Performance: n=2,000,000 generates (with full exact reference +
  stats) in ~1.1s; the optional n=30,000,000 case in ~19s. Both are
  included in the delivered suite.

## Suite contents (`suite_out/`, from `collharness suite --out-dir suite_out --large`)

| file | n | root_bytes | key_bits | bucket_hash | notes |
|---|---|---|---|---|---|
| tiny_basic4 | 64 | 4 | 24 | 0 | edge cases only |
| tiny_basic8 | 64 | 8 | 40 | 0 | edge cases, 8-byte keys |
| tiny_hashmode1 | 64 | 4 | 24 | 1 | edge cases, multiplicative hash |
| tiny_skew_mask | 329 | 4 | 24 | 0 | +300-key bucket skew, mask hash |
| tiny_skew_mix | 329 | 4 | 24 | 1 | +300-key bucket skew, mult. hash |
| med_default4 | 2,000,000 | 4 | 24 | 0 | dominated by *accidental* birthday collisions (small 24-bit key space at 2M scale) -- deliberately left this way; see note below |
| med_default8 | 2,000,000 | 8 | 40 | 0 | 40-bit key space -- collisions are almost entirely the intentional density-driven filler pairs, not accidental |
| med_hashmode1 | 2,000,000 | 4 | 24 | 1 | same as med_default4 but multiplicative hash; exercises the comparator's saturation path (found_count > 999) |
| large_default4 | 30,000,000 | 4 | 27 | 0 | optional large case (0.3/0.7) |

**Note on med_default4/med_hashmode1:** with `key_bits=24` (a 2^24 ≈
16.8M-value key space) and n=2,000,000, random filler singleton keys
collide with each other by chance far more often than the 1% filler
`--density` alone would produce (birthday-paradox effect) -- that's
why `dedup_count` for these two is ~125k rather than the ~20k seen in
`med_default8` (40-bit key space, collisions negligible except by
design). Both are valid, useful, and different test data -- one
stresses "lots of small accidental duplicate groups", the other
"collisions strictly by density knob" -- so both are kept in the
suite rather than treated as a bug. Widen `--key-bits` when generating
custom cases if a lower accidental-collision rate is wanted.

## End-to-end real-engine recipe (plan task 0.8)

This harness produces and checks purely synthetic data. Producing and
comparing an actual `test_ad` cell's hits from YAFU itself is a
separate, smaller recipe, using the engine's own test-mode dump
(`stage1.c`'s `test_%s.hits` files, one line per hit: `a_d p m`,
post-filter):

1. Run YAFU's existing (CUDA) `gpu_gerbicz` engine with
   `nfs_args="test_ad=<AD> test_pmin=<PMIN> test_pmax=<PMAX> test_qmin=<QMIN> test_qmax=<QMAX>"`
   for a chosen cell. This produces `test_gerbicz.hits` (the file name
   is `test_<engine-name>.hits`, taken from the selected engine's
   `v->name` -- see `stage1_engine_select()`).
2. Once the OpenCL engine exists (Phase 3+) and is registered, rerun
   the identical `nfs_args` selecting the OpenCL engine; it will
   produce its own `test_<opencl-engine-name>.hits`.
3. Compare: sort both files' lines and diff them. Because a hit's
   fields are exact integers (`a_d`, `p`, `m` as decimal `mpz`
   values) textual sort+diff is sufficient; no floating-point
   tolerance is involved at this level (the floating-point coefficient
   filtering already happened inside `check_found_array` before a hit
   is ever written).
4. A proposed cell list, to cover the four regimes named in the plan
   (small/<480 bits, larger, pp32 p_max<65536, pp64 p_max>=65536), is
   an open item for the person running this recipe to supply real
   `test_ad`/`N` values for -- see STATUS.md "Inputs and blockers".
   This harness's synthetic suite is not a substitute for this step;
   it validates the *collision engine's internal logic* in isolation,
   while this recipe validates the *whole pipeline* (trans kernels +
   collision engine + `handle_collision`) end to end.

## Known limitations / assumptions

- No bucket-capacity overflow or `CANDIDATE_CAP`/`VALUE_MATCH_CAP`
  overflow paths are modeled -- `cpuref_stats()` assumes an
  effectively unbounded bucket/candidate/value-match capacity
  (equivalent to the real engine's `ensure_capacity` retry loop
  having already converged). Stress-testing those overflow paths is
  called out as future work (matches the plan's existing "Open
  issues": hash-cap and candidate-overflow behavior belongs to
  Phase 3/7, not this harness).
- `cpuref_exact`'s O(n log n) sort via `qsort` and `cpuref_stats`'s
  per-key multiplicity lookup via binary search are adequate for this
  delivery's sizes (30M elements in under 20 seconds) but are not
  algorithmically optimal; fine for an offline test harness, flagged
  here in case a much larger stress case is wanted later.
