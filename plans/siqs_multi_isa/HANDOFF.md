# Handoff: SIQS large-prime strategies, kernels and multi-ISA port

*Written 2026-10-22 at commit `ce6f673` (branch `L2sieve`, pushed to
`origin/L2sieve`). This is the context document for continuing this work
in another session. Line numbers are as of `ce6f673`. Search by symbol
name if they drift.*

Related documents:
- `direct_sieve.md`: the direct interval sieve design and its history.
- `Sink-Stripe Sieve Stores.md`: the sink-stripe idea from MPU-GMP.
- `plans/siqs_multi_isa/plan.md`: the plan. Section 0 holds the final
  status; section 5.1 records the correctness-gate change.
- `L2_SIEVE_DEBUG_NOTES.md`: older notes on the L2 / sparse sieve.

---

## 1. What this effort was

YAFU's SIQS sieves the factor base in ranges. This work concerns the two
ranges above `med_B`:

| range (fb index) | primes | strategies |
|---|---|---|
| `[med_B, large_B)` ("mid") | ~32k .. sieve interval | **bucket** (stock medbucket) or **sparse** (per-block sparse sieve) |
| `[large_B, B)` ("large") | > sieve interval | **bucket** (stock largebucket) or **direct** (direct interval sieve) |

The work happened in two stages.

### 1.1 Kernel experiments (AVX-512 only; commits `5b9841d`, `6ca8388`, `46de3d0`)

Earlier in this session, the direct sieve's scatter of roots larger than
the interval was rewritten and benchmarked. All timings are QS seconds,
i5-11600K, `-threads 1`, fixed seed.

**Scatter kernels:**
- **bitloop** (original): compare, then a `tzcnt`/`blsr` loop over the
  mask bits.
- **sink stripe**: every lane does a read-modify-write; misses are blended
  to a rotating 1 KB stripe in the second tail guard block.
- **creg**: compress in a register, unaligned store, batched.
- **cmem**: `vpcompressd` straight to memory, batched.

| input | NB | bitloop | sink | creg | cmem | compress, unbatched |
|---|---|---|---|---|---|---|
| c77 | 2 | 26.40 | 23.75 | 20.95 | **20.67** | 24.08 |
| c77 | 4 | 23.02 | 21.22 | **20.06** | 20.11 | 22.67 |
| c77 | 8 | 23.50 | 22.75 | **22.61** | 22.98 | 24.28 |
| c79 | 2 | 47.18 | 42.42 | 35.42 | **35.03** | 41.45 |
| c79 | 4 | 39.27 | 38.43 | 34.71 | **33.97** | 38.97 |
| c79 | 8 | 38.55 | 38.16 | 37.13 | **36.66** | 40.22 |

The gain comes from batching: one loop exit per about 500 hits instead of
four per 16 primes. Compressing without batching is no better than the
sink stripe. The batch size (`DIRECT_CBUF_FLUSH`) is flat from 256 to
4096; it is set to 512.

**Compress stores in the sparse sieve and its tdiv.** With these, the
sparse mid range roughly ties the bucket mid range (within ±2% at NB 2
and 4; mixed at NB 8). Bit loops were 3–16% slower than compress.

### 1.2 Multi-ISA effort (this handoff's main subject; commits `edbaf44`..`ce6f673`)

**Goal.** Make the strategies and kernels selectable at run time and port
them to AVX2, SSE4.1, SSE2, generic C and NEON. Then pick defaults per
CPU and ISA from measurements.

**Done.** All of that, plus build and benchmark tooling, a tuner
extension, an NB cap, and two bug fixes found along the way.

---

## 2. Commits (oldest first)

| commit | what |
|---|---|
| `5b9841d` | sink store vs compress-store kernels for the direct sieve (compile-time `DIRECT_SCATTER_MODE`) |
| `6ca8388` | compress stores in the sparse sieve and tdiv_large; `DIRECT_MEDBUCKET 0` default |
| `46de3d0` | plan doc (plus your tuning script) |
| `edbaf44` | **phase 0**: ISA fences, per-core L2 detection, `tools/` build and bench scripts |
| `3b4157d` | `.gitattributes`: LF line endings for `*.sh` and `tools/*.py` |
| `72b8c9f` | **phase 1**: strategies and kernels chosen at run time (AVX-512) |
| `5a053fd` | **phase 2**: scan hook, scalar kernels, portable root updates; strategies at every ISA |
| `77fed87` | **phase 3**: AVX2 kernels, plus the fix for phase 2's sparse tdiv order bug |
| `26830c0` | **phases 4+5**: 4-lane template with SSE4.1 and NEON backends, plus NEON emulation |
| `ce6f673` | **phase 6**: defaults table, NB cap, `ytools_cpu_is_amd`, `siqs_tune.py` params, docs |

---

## 3. How it works now

### 3.1 Compile-time structure

- **`include/qs_impl.h:94` `SIQS_EXPERIMENTAL_SIEVES`** is defined for every
  ISA unless `-DSIQS_NO_EXPERIMENTAL_SIEVES` is given (which builds the
  stock-only reference). It switches on `TRY_DIRECT_SIEVE`, the one macro
  for both the sparse and the direct sieve. `TRY_L2_SIEVE`,
  `TRY_SPARSE_SIEVE` and `TRY_LIST_SIEVE` were folded into it or removed
  after `858a047`, together with the whole-interval L2 sieve and the
  linked-list sieve.
- **Each kernel ISA's code lives in its own file**, compiled when the
  build has that ISA:

| file | content | compiled when |
|---|---|---|
| `factor/qs/direct_sieve_avx512.c` | AVX-512 kernels: helpers `direct_cstore16` :36 / `direct_sink_scatter16` :77 / `direct_scatter16` :128; pass 1 :207; fused :284 (prefilter :389, bitmap gather :438); sparse :596; tdiv :760 | `USE_AVX512BW` |
| `factor/qs/direct_sieve_avx2.c` | AVX2 kernels: scatter :72, resieve test `direct_rs8` :243, fused :266, sparse :467, tdiv :572, sparse root update :679 | `USE_AVX2` (AVX2 and AVX-512 builds) |
| `factor/qs/direct_sieve_v4.h` | 4-lane template: scatter :72, fused :201, sparse :363, tdiv :470, root update :576 | included by the next two |
| `factor/qs/direct_sieve_sse41.c` | SSE4.1 backend (pshufb compress) | `USE_SSE41` |
| `factor/qs/direct_sieve_neon.c` | NEON backend (vqtbl1q compress, vaddvq masks) | AArch64 `__ARM_NEON`, or `-DSIQS_NEON_EMU` on x86 |
| `factor/qs/arm_neon_emu.h` | scalar emulation of the dozen NEON intrinsics used | `SIQS_NEON_EMU` only |
| `factor/qs/direct_sieve.c` | portable code: survivor pool :42-200, `direct_block_sub4` :205, scalar kernels :244-660, portable root updates :669-815 | always (with `TRY_DIRECT_SIEVE`) |

- **Kernel pattern.** Each kernel is a `SIQS_FORCEINLINE` body taking
  `const int kmode`, plus thin per-kernel wrappers (`*_bitloop`, `*_sink`,
  `*_creg`, `*_cmem`). The compiler specializes each wrapper, so the hot
  loops have no branches on the kernel choice. Shared helpers in
  `qs_impl.h`:
  - `SIQS_FORCEINLINE` :1004
  - `direct_compress_flush` :1014
  - `siqs_ctz32` / `siqs_popcnt32` (portable; MSVC has no builtins)
- **Kernel names** are `<kernel>_<isa>_<mode>`, e.g.
  `direct_resieve_update_avx512_cmem`, `tdiv_large_avx2_creg` (4-lane:
  `<kernel>_<mode>_<sfx>`). The AVX-512 kernels were moved out of
  `med_sieve_32k_avx2.c` / `tdiv_med_32k_avx2.c` and renamed from
  `*_32k_*` / `tdiv_large_32k_avx2_*` after `ce6f673`.
- **4-lane prototypes** are declared through `SIQS_DECLARE_V4_KERNELS(sfx)`
  (`qs_impl.h:1067`). `SIQS_HAVE_NEON_KERNELS` is set at :1087.

### 3.2 Run time

**Options** (`top/cmdParser/cmdOptions.c:1349`..:1385; values in
`include/siqs_strategy.h`; copied to `fobj->qs_obj.gbl_siqs_*` in
`top/driver.c`). All work in yafu.ini too.

| option | values | meaning |
|---|---|---|
| `-siqsMid` | `bucket` \| `sparse` | strategy for `[med_B, large_B)` |
| `-siqsLarge` | `bucket` \| `direct` | strategy for `[large_B, B)` |
| `-siqsKernel` | `bitloop` \| `sink` \| `creg` \| `cmem` | direct scatter kernel |
| `-siqsSparseKernel` | `bitloop` \| `creg` \| `cmem` | sparse sieve / tdiv_large kernel |
| `-siqsLPisa` | `scalar` \| `sse41` \| `avx2` \| `avx512` \| `neon` | kernel ISA (default: best available) |

On 4/8-lane ISAs `cmem` means `creg`. In scalar code, `creg`/`cmem` mean
the branch-free append.

**`siqs_resolve_strategy`** (`factor/qs/SIQS.c:4409`; called at :5737,
after `scan_ptr` is chosen):
1. Picks the kernel ISA from the build and CPU flags.
2. Fills in defaults from the table at :4447 (section 5).
3. Sets the fb bounds with `siqs_lp_boundary` (:4393). The direct sieve
   rounds `large_B` up, so every prime in its range is larger than the
   interval; buckets round down. `fb_blk8i = med_B` when the mid range uses
   buckets, which empties every sparse loop.
4. Sets `sconf->lp_direct`, `mid_bucket`, `lp_kernel`, `sparse_kernel`,
   `lp_isa` and `lp_portable_roots` (`qs_impl.h:639`ff).
5. Sets the kernel function pointers (`SIQS.c:74`ff):
   `large_sieveblock_ptr`, `tdiv_large_ptr`, `direct_sieve_largeprimes_ptr`,
   `direct_resieve_update_ptr` and `nextRoots_lp_ptr`. The 4-lane selection
   is the `SIQS_SELECT_V4` macro at :4578.

**NB cap** (`SIQS.c:5214`): with AVX2 or AVX-512 kernels, no `-siqsNB`, and
not both ranges forced to buckets, NB is capped at L2/(4·32 KB). That is 4
on a 512 KB L2.

**`-v`** prints a `strategy:` line (from `print_siqs_experiments`,
`SIQS.c:3399`), for example:
`strategy: mid = sparse, large = direct, kernel = cmem, sparse kernel = cmem (avx512 kernels)`.
`print_fb_bounds` (:4647) now runs after the strategy is resolved
(`-v -v`).

### 3.3 Per-poly flow (`process_poly`, `SIQS.c:1779`)

- **Direct pre-checks** (:1861): a-primes must lie below `large_B`, or
  below `med_B` when the mid range uses buckets.
- **Direct pre-block** (:2448):
  - For the first poly of an 'a', memset both intervals and run pass 1
    (`direct_sieve_largeprimes_ptr`, :2462).
  - Apply the minblock correction with `direct_block_sub4`.
- **Block loops** (positive side, then negative):
  - `blk_sieve` is `sieve + i*32768` for the direct sieve, else the single
    block window.
  - The med sieve runs, then the sparse sieve (`large_sieveblock_ptr`,
    :2577 / :2600 / :2753 / :2776).
  - `lp_sieveblock_ptr` runs only if the large range uses buckets or the
    mid range does.
  - `blk_sieve[65535]` is saved and restored around the block
    (`direct_dummy_byte`, :2555 / :2731), because the med sieve dumps
    a-prime logs there.
  - `dconf->sieve` is pointed at the block for the scan.
- **Fused pass** (:2827): when another poly follows, memset both
  intervals, then run `direct_resieve_update_ptr` (:2843), then
  `direct_finalize_survivors` (:2845).
- **nextRoots** (:3047ff):
  - The AVX-512 path (`nextRoots_32k_avx2_small` + knl medbucket /
    largebucket) runs unless `lp_portable_roots` is set.
  - Otherwise small primes go through `nextRoots_32k_avx2_small` (AVX2
    builds; it stops at med_B without AVX-512) or
    `nextRoots_32k_generic_small`, then `nextRoots_lp_ptr` (:3129 / :3149).
  - Stock bucket+bucket on lower ISAs uses the stock `nextRoots_ptr`.
- **firstRoots** (`factor/qs/poly_roots_32k.c`): with the direct sieve, no
  slice bookkeeping or bucket fill for `[large_B, B)` (:2447, :2484).

**Scans** (`factor/qs/tdiv_scan.c`):
- Every non-AVX-512 scan ends each report with `siqs_report_lp` (:30). It
  runs tdiv_LP for bucketed primes and, with the direct sieve, parks the
  survivor through `direct_save_survivor`.
- Those scans call `SIQS_SPARSE_TDIV` (:50) **after** `resieve_med_ptr`;
  see bug 1 in section 7.
- The AVX-512 scan (`check_relations_siqs_16_avx512`, :995) has its own
  equivalent at :1176.

**Survivor pool** (`direct_sieve.c`):
- `direct_pool_alloc` :42, `direct_save_survivor` :99,
  `direct_resieve_hit` :126 (bitmap, then binary search, then a guarded
  division), `direct_finalize_survivors` :170.
- Allocated at `SIQS.c:3870`. The sieve buffer is `(2·NB+4)` blocks
  (:3951). The sink stripe uses the second tail guard block.

**Portable root updates** (`direct_sieve.c`):
- `nextRoots_32k_lp_sparse_scalar` :669 updates the sparse roots in
  min/max order into `sieve_p` / `sieve_n`.
- `nextRoots_32k_lp_buckets` :713 fills the mid and large buckets in
  firstRoots' element format and slice bookkeeping.
- `nextRoots_32k_lp_portable` :810 calls both.
- The AVX2 and 4-lane versions replace only the sparse part.

### 3.4 Other code touched

- **`ytools/ytools.c`:**
  - `ytools_get_l2_core_size` :866 returns the true per-core L2: CPUID
    leaf 4 on Intel, 0x80000006 on AMD. The existing "L2" value is
    max(L2, L3) on purpose, because msieve's linear algebra uses it, and
    0x80000006 reports a stale 256 KB on Rocket Lake.
  - `ytools_cpu_is_amd` :900.
- **`factor/qs/update_poly_roots_32k_avx2.c:3263, :3423`:** the AVX-512 L2
  root-update sections are now guarded by `USE_AVX512F`.
- **Build lists:** `Makefile`, `Makefile.gcc` and
  `build.vc22/siqs{,-clang}/siqs.vcxproj` list the four new `.c` files.
  The MSVC projects are edited but **not build-tested**.
- **`siqs_tune.py`:** `Mid`, `Large` and `Kernel` are string-categorical
  parameters (`STR_CHOICES`), parsed from the `strategy:` line.

---

## 4. Tooling and environment

**Environment:**
- Builds and runs are in WSL Ubuntu 22.04, gcc 11.4.
- GMP is at `../gmp-install/6.2.1` and ECM at `../ecm-install/7.0.6-wsl`
  (relative to the repo; `build_matrix.sh` makes these absolute).
- `config.mk` sets `USE_AVX512IFMA := 1`, so the in-tree `make yafu` is an
  IFMA build.

**`tools/build_matrix.sh [tags]`:**
- One rsync'd copy and build per ISA level, in `$BUILD_ROOT/<tag>/src`;
  the binary goes to `$BUILD_ROOT/<tag>/yafu` (default `~/yafu_matrix`).
- Tags: `avx512 avx2 sse41 sse2 generic avx512stock avx512neonemu`.
- Pass extra make arguments with `MAKE_ARGS="USER_CFLAGS=-D..."`.
- **Caveat:** make doesn't track `USER_CFLAGS` changes. Delete
  `$BUILD_ROOT/<tag>` when switching `-D` flags (this bit me once).

**`tools/bench_matrix.sh [-i inputs] [-n "NB list|def"] [-r reps] [-o csv] [-p] label=bin[,args] ...`:**
- Runs one job at a time with `-seed 777 -noopt` (`-s none` turns that
  off), deletes `siqs.dat` before each run, and refuses to start if a yafu
  is already running.
- Writes a CSV and keeps the last log of each variant as
  `~/siqsbench/last_<label>_<input>_<nb>.txt`.
- `-p` wraps runs in `perf stat`.

**`tools/bench_summary.py csv [--base L] [--groups "a,b;c"]`:**
- Prints a table of means with % against the base variant.
- Checks that work counters (polys, touched, blocks, full, partial) are
  identical within each group, and flags runs that didn't factor.

**Also in `tools/`:**
- `tools/bench_inputs.txt`: c70, c77 (256-bit), c79 (261-bit), c85 and
  c95, each with its factors.
- `tools/bench_ref.txt`: the reference counters.

**Profiling** (WSL): perf doesn't resolve symbols for this binary. Map
addresses with `nm -n` instead; a small bisect script is in the session
history. The pattern:

```bash
perf record -o p.data ./yafu "siqs(@)" ...; nm -n yafu | grep " t "
```

**Quick check that everything is correct:**

```bash
tools/build_matrix.sh avx512 avx512neonemu
```

```bash
tools/bench_matrix.sh -i c77 -o /tmp/chk.csv a=~/yafu_matrix/avx512/yafu "v2=~/yafu_matrix/avx512/yafu,-siqsLPisa avx2" "v4=~/yafu_matrix/avx512/yafu,-siqsLPisa sse41 -siqsMid sparse" "s=~/yafu_matrix/avx512/yafu,-siqsLPisa scalar -siqsMid sparse" "n=~/yafu_matrix/avx512neonemu/yafu,-siqsLPisa neon -siqsMid sparse"
```

```bash
python3 tools/bench_summary.py /tmp/chk.csv
```

These pass the `-siqsNB 4` default, so the NB cap doesn't apply. All
variants must show 89994 polys and 4698192 touched. The 4-lane and scalar
variants need `-siqsMid sparse`, because their default mid range is bucket
(87317 polys).

---

## 5. Defaults table (what auto chooses)

From `SIQS.c:4447`:

| kernel ISA | `[med_B, large_B)` | `[large_B, B)` | kernel | NB |
|---|---|---|---|---|
| AVX-512 | sparse if 2·NB·32 KB ≤ L2, else bucket | direct | cmem (Intel), creg (AMD) | capped at L2/128 KB |
| AVX2 | sparse if it fits in L2, else bucket | direct | creg | capped at L2/128 KB |
| SSE4.1 / NEON / scalar | bucket | direct | creg / append | yafu's table |

`-DDIRECT_DEFAULT_MID` / `_LARGE` / `_KERNEL` / `_SPARSE_KERNEL` override
the table at compile time.

---

## 6. Final analysis (i5-11600K, 512 KB L2; WSL; `-threads 1`; single runs)

### 6.1 Correctness gate (all pass)

- **Stock builds at different ISA levels don't do identical work.** Their
  med-sieve kernels approximate logp differently: on c77 at NB 4, the
  polynomial counts are 93182 (AVX2), 91142 (SSE4.1), 90887 (SSE2) and
  90248 (generic). So cross-ISA bit-identity is only possible for the
  large-prime kernels.
- **Inside one AVX-512 binary,** `-siqsLPisa avx512|avx2|sse41|scalar|neon`
  (neon from the `avx512neonemu` build), with every kernel, gives
  identical counters for sparse+direct: c77 89994/4698192, c79
  150705/5895674.
- **bucket+bucket equals the stock build** (`avx512stock`) exactly, and on
  each lower-ISA build equals that build's stock counters.
- **Mid buckets via the portable root update** differ by about 0.01% from
  the knl path (different slice logp bookkeeping) and factor fine.
- **Every build factors c70, c77, c79 and c85** under every strategy
  tested.

### 6.2 New defaults against each build's own stock path, at yafu's NB choice (QS seconds)

| build | c77 | c79 | c85 |
|---|---|---|---|
| AVX-512 | 20.2 vs 21.5 (−6%) | 33.6 vs 36.8 (−9%) | 110.2 vs 141.3 (−22%) |
| AVX2 | 25.1 vs 28.1 (−11%) | 44.1 vs 48.2 (−9%) | 145.1 vs 175.2 (−17%) |
| SSE4.1 | 34.7 vs 35.6 (−3%) | 60.5 vs 62.6 (−3%) | — |
| SSE2 | 36.8 vs 39.3 (−7%) | 65.1 vs 69.7 (−7%) | — |
| generic | 45.2 vs 48.3 (−6%) | 76.8 vs 82.8 (−7%) | — |

The c85 rows use the NB cap (NB 4), against stock at its table NB (7 and
11). Run-to-run noise was several percent on this machine; one unchanged
SSE4.1 configuration varied by 9% between two sessions.

### 6.3 What the measurements say

1. **The direct sieve wins for primes above the interval at every ISA and
   NB tested.** Large buckets (`-siqsLarge bucket`) were the slowest option
   everywhere: +5–7% against stock at NB 12/16 on AVX-512, +35–75% on AVX2.
2. **The sparse mid range needs wide vectors and an interval that fits in
   L2.**
   - AVX-512 at NB 12/16 (768 KB / 1 MB): sparse+direct was +3–4% against
     stock, while bucket+direct tied stock.
   - With 4-lane kernels (SSE4.1), bucket+direct beat sparse+direct by
     2–4%.
   - Scalar sparse is far behind: its tdiv walks every sparse prime for
     every report.
3. **NB: smaller is better once the direct sieve is used.** yafu's NB
   table (6 on AVX-512, 9–11 elsewhere) was tuned for buckets.
   - AVX-512 / AVX2: NB 4 beat the table's NB by 4–20%.
   - 4-lane and scalar builds still prefer NB 9.
4. **Batching is what makes compress kernels win.** On wide ISAs: AVX-512
   cmem ≈ creg > sink > bitloop; AVX2 creg > sink (by about 10%) >
   bitloop (+12 to +41% against stock). The same order held at 4 lanes,
   with smaller gaps.
5. **Remaining cost.** The fused pass (resieve, root update, next-poly
   scatter) is the largest remaining cost: about 29% of the AVX2 profile
   and 58% of the scalar one. The med sieve is next at about 28% (AVX2).

### 6.4 Performance notes on the scalar and 4-lane code

- **gcc turns branch-free ternaries back into branches.** Use arithmetic
  masks, as in `DIRECT_SCALAR_HIT` (`direct_sieve.c:384`) and the root updates. That change alone took the scalar fused pass
  from 64.8 s to 44.3 s on c77.
- **The 4-lane resieve test is scalar and branch-free**, with no
  prefilter. AVX2 uses a 512-bit register prefilter plus masked gathers;
  AVX-512 uses a 1024-bit prefilter plus gathers.

---

## 7. Bugs and lessons

1. **Sparse tdiv order (phase 2 bug, fixed in `77fed87`).**
   - **Bug:** the new lower-ISA scans ran `tdiv_large_ptr` before
     `resieve_med_ptr`, so resieved primes (below med_B) were appended
     after larger sparse primes.
   - **Why it mattered:** `filter.c:302` merges each relation's
     `fb_offsets` with the a-prime list assuming both are sorted.
   - **Symptom:** "odd exponent found" in the sqrt; c79 failed on the AVX2
     build at NB ≥ 4. The relations still passed `check_relation`. The
     SSE4.1 build got away with it.
   - **Rule:** factor indices must be appended in ascending fb order.
2. **The a-prime hypothesis was wrong.** a-primes were all below med_B
   here. The direct pre-check still guards against them.
3. **Small fb primes left in Q.** About 0.2% of reports leave one in Q,
   identically for stock and all strategies (`DEBUG_MISSED_FB_PRIMES`,
   `tdiv.c:254`). It predates this work.
4. **Line endings.**
   - Git Bash `sed -i` on Windows rewrites CRLF files as LF; it happened
     twice and was reverted. Edit with the Edit tool or Python preserving
     `\r\n`.
   - Running git from WSL sees every file as modified (autocrlf and
     filemode differ). Run git from the Windows side.
5. **Shell escaping.** Heredocs through this tool layer mangle `\\n` and
   backslash-newline. Write scripts to files rather than inline heredocs.
6. **Splitting phases from a mixed working tree.** The rsync'd build
   copies under `~/yafu_matrix/<tag>/src` are exact snapshots of what was
   tested. They were used to commit phases 1, 2 and 3 separately (excluding
   `build.vc22`, which isn't rsync'd).

---

## 8. Open items and suggested next steps

1. **Other hardware.**
   - Run `tools/build_matrix.sh` + `bench_matrix.sh` on AMD Zen 4/5, to
     check cmem against creg and the table's AMD rule, which is
     unmeasured.
   - Also on a newer Intel core with a 1–2 MB L2 (tests the L2 rule and NB
     cap), and on an older AVX2-only Intel and Zen 2/3.
2. **ARM.** Build natively on AArch64 (Apple M or Graviton):
   - Use a plain build with no `USE_*` x86 flags and the native compiler,
     so `__aarch64__ && __ARM_NEON` enables `direct_sieve_neon.c`. It is
     untested whether the Makefile needs ARM-specific adjustments.
   - Check the counters against `bench_ref.txt` with
     `-siqsLPisa scalar` vs `neon`, then time them.
   - The NEON backend is only verified through `arm_neon_emu.h`.
3. **Larger inputs** (c95+) and **multi-threaded** runs (`-t N`): neither
   is measured. The NB cap and the L2 rule assume one thread per core.
4. **NB policy.**
   - The cap is a simple rule fitted to c77–c85.
   - `siqs_tune.py` can now tune `Mid` / `Large` / `Kernel` together with
     NB; build a proper per-size NB table for the direct sieve.
   - The 4-lane / scalar builds may want a different NB too.
5. **Speed.**
   - The fused pass is the top remaining cost on AVX2 and scalar.
   - Ideas: an AVX2 gather-free survivor test with a larger register
     prefilter; vectorize the scalar resieve test with SSE2 where
     available.
   - `tdiv_large_scalar` walks every sparse prime per report, which is
     why scalar builds use mid buckets.
6. **MSVC / Windows builds.** The `.vcxproj` files list the new sources,
   but they are untested. Watch for `SIQS_FORCEINLINE`, the intrinsics in
   `direct_sieve_avx2.c`, and `siqs_ctz32`.
7. **Possible cleanup.**
   - Done: the `TRY_LIST_SIEVE` and whole-interval `TRY_L2_SIEVE` code is
     removed, and `TRY_L2_SIEVE` / `TRY_SPARSE_SIEVE` are folded into
     `TRY_DIRECT_SIEVE`.
8. **Calibration.** Bucket and direct use different logp conventions (per
   slice vs per 16-group; `DIRECT_LOGP_BIAS`), which mixes yield into the
   speed comparisons. A common convention would make strategy comparisons
   cleaner.
