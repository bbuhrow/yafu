# OpenCL `gpu_gerbicz` port — STATUS

This file travels from phase to phase. Attach it (with `OCL_GERBICZ_PLAN.md`) to every new conversation.
**Update protocol:** at the end of each conversation, rewrite the sections marked *(update)*. Keep it short: fold old detail into one-line summaries. Do not delete "Locked decisions" or "Facts" without recording why.

- **Last updated:** 2026-09-23, end of Phase 0
- **Current phase:** 1 (not started; prompt ready — see `OCL_GERBICZ_PHASE1_PROMPT.md`)
- **Target:** AMD GPUs without CUDA (dev card RX 6700 XT, OpenCL 2.0). Gerbicz engine only.

## Phase table *(update)*

| Phase | Title | Status | Conversation notes |
|---|---|---|---|
| 0 | Baseline and harness | **Done** | Portable C99 `collharness` tool built, tested (gcc 13 + clang 18, `-Wall -Wextra -pedantic`, zero warnings). All tasks 0.1–0.10 complete except 0.9 (no NVIDIA machine available) and part of 0.8 (real `test_ad`/N values still needed from the user — see Inputs). |
| 1 | OpenCL foundation | Not started | Needs `ocl_xface.h` |
| 2 | Primitives (scan, reduce-max, fill, radix sort) | Not started | |
| 3 | Collision engine (3a–3d) | Not started | |
| 4 | Sieve (trans) kernels | Not started | |
| 5 | Host integration and registry | Not started | |
| 6 | Validation | Not started | |
| 7 | Performance and hardening | Not started | |

Status values: Not started / In progress / Done / Done with caveats.

## Locked decisions
(Full text in plan §2.)
1. OpenCL 2.0 minimum; sub-groups optional and gated.
2. Static link, runtime-loaded `.cl` files, no DSO.
3. `ocl_` symbol prefix; CUDA and OpenCL builds mutually exclusive for now.
4. Base-offset kernel arguments instead of sub-buffers or SVM.
5. Shared context and program per device; per-thread queue and kernel objects.
6. No-sub-group baseline collision kernels.
7. Bit-exact first, then optimize; behavior changes logged here.
8. Registry: `STAGE1_ENGINE_OCL_GERBICZ` / `ocl_gerbicz` / `HAVE_OCL_POLY` (proposed; confirm in Phase 5).

## Facts (verified from source before Phase 0)
(Detail and line references in plan §3.)
- Keys are signed (two's complement); sort full 32/64-bit width, not `key_bits`.
- Key 0 = empty slot. Packed value = `(q_index << shift) | p`; `q_index` indexes `q_batch`.
- `found_array[0].p1` counts all attempted stores; only 999 entries are stored.
- `STAGE1_OVERFLOW_SKIP` is declared but unused; CUDA engine `exit(-1)`s on candidate and value-match overflow.
- CUDA hash-word cap comes from the shared-memory opt-in limit; OpenCL local memory (64 KB on the dev card) caps at 4096 words.
- Filter output is order-independent, so a sequential CPU model can match it exactly.
- `stage1_core.cu:502` truncates a 64-bit value to `uint32` (pp64_r64). Port verbatim.
- Intrinsics port exists; `modinv64` (extra out-param) and `montmul64` (`uint w`) signatures differ.
- Cofactorization cache key is device name only; needs a source hash and driver version.
- **(Phase 0 addition)** `test_%s.hits` dump filenames are `test_<engine v->name>.hits` (from `stage1_engine_select()`), not the literal `test_cpu_gerbicz.hits` / `test_ocl_gerbicz.hits` named in the original Phase 0 prompt — that prompt's task 0.8 wording was a simplification. No correction needed to any locked decision; just note the real naming convention when writing the actual end-to-end recipe (done, see harness README's "End-to-end real-engine recipe" section).
- **(Phase 0 addition)** The engine's D/S/X secondary hash tables (used to find, for every original element, which deduped candidate key it matches) are an exact, collision-free lookup — not a source of approximation, unlike the per-bucket filter. This let the CPU stats model compute `dedup_count`/`value_match_count` directly from the candidate set's exact-duplicate structure instead of also simulating open addressing, with no loss of fidelity (assuming no capacity overflow — see Open issues). Cross-checked, not just asserted: verified by hand against the deterministic edge-case suite and by independent round-trip (`gen` → strip → `ref` → `cmp`) on every suite case.

## Inputs and blockers *(update)*

| Item | State |
|---|---|
| `ocl_xface.h` (OpenCL `gpu_info_t`, `gpu_launch_t`, `OCL_TRY`, `clGetErrorString`) | **Missing.** Needed for Phase 1. |
| Sample YAFU invocation / N / `nfs_args` for end-to-end `test_ad` cells | **Still not provided.** The harness's recipe (README "End-to-end real-engine recipe") is ready to use once these are supplied: need `test_ad`, `test_pmin`, `test_pmax`, `test_qmin`, `test_qmax` for 4 regimes (small <480 bits, larger, pp32 p_max<65536, pp64 p_max>=65536). |
| NVIDIA machine for CUDA golden dumps | Unknown (optional, task 0.9 not attempted). |
| Project files present | `collision_engine.cu/.h`, `collision_bucket.h`, `cuda_intrinsics.h`, `stage1_core.cu/.h`, `stage1_engine.c/.h`, `stage1.c`, `stage1_sieve_gpu.c`, `gpu_cofactorization_cl.c/.h`, `opencl_intrinsics.cl` |
| Harness repo path convention in the user's actual repo | Not asked/answered this phase; harness delivered standalone, needs a home directory decided before Phase 1 assumes one. |

## Interfaces and formats *(update — Phase 0 fills the collcase spec)*
- **collcase v1 file format:** Done — see `include/collcase.h` in the delivered harness for the byte-precise spec (64-byte header, `keys[n]`/`values[n]`/`q_batch[num_q]` arrays, optional expected-output section with `found_t`-shaped entries + 424-byte stats block). All integers little-endian, written/read field-by-field (never a raw struct `fwrite`), so it's endianness- and padding-independent.
- **Harness tool CLI:** Done — single binary `collharness` with subcommands `gen`, `ref`, `cmp`, `suite`, `info`. See harness README.
- **Canonical result form and comparison rules:** Done — `collharness cmp` canonicalizes (p1<p2 within an entry; sort by (q,qroot,offset,p1,p2)) and applies a saturation-aware rule: exact set match below the 999-entry cap, subset check plus found_count equality at/above it. Verified against a real saturated case (med_hashmode1, found_count=19954).
- **Stats-model parameters (hash-word cap override etc.):** Done — `--hash-word-cap` CLI flag / `collcase_t.hash_word_cap` field, 0 = default 4096 words (derived from the 64 KiB / 3-tables-of-uint32 AMD local-memory limit, matching the Facts entry above).

## Artifacts produced *(update)*

| Phase | File | Purpose | Location |
|---|---|---|---|
| 0 | `collharness/` (full source tree: `include/`, `src/`, `Makefile`, `README.md`, `run_tests.sh`) | Portable C99 test harness: collcase v1 I/O, generator, CPU exact reference, CPU stats model, comparator, CLI | delivered as `collharness.zip` |
| 0 | `collharness/suite_out_tiny/*.collcase` | The 5 tiny/skew named-suite cases (cheap, <6 KB each), bundled directly | inside `collharness.zip` |
| 0 | `med_default4.collcase`, `med_default8.collcase`, `med_hashmode1.collcase` (~19–27 MB each), `large_default4.collcase` (~250 MB, n=30,000,000) | The med/large named-suite cases | **not bundled** (too large to ship inline) — fully reproducible in ~1s (med) / ~19s (large) via `collharness suite --out-dir DIR --large`; see harness README's suite table for exact parameters/seeds (seeded by case name via FNV-1a, so regeneration is byte-identical) |
| 0 | `OCL_GERBICZ_PHASE1_PROMPT.md` | Phase 1 kickoff prompt | delivered alongside this file |

## Open issues and risks *(update)*
- Candidate overflow on the OpenCL path may be more likely if the local-memory cap is lower than CUDA's; needs a recoverable skip (Phase 3).
- `MAX_MEM_ALLOC_SIZE` may be smaller than the two bucket arrays or the root arrays at the current sizing (Phase 5).
- AMD program compile time is long; the binary cache must be robust (Phase 1).
- Sorting fewer than the full key bits is unproven; do not do it before checking (Phase 2 / 7).
- **(Phase 0 addition)** The CPU stats model (`cpuref_stats()`) assumes bucket/candidate/value-match capacity is never exceeded (no `ensure_capacity`-retry-loop equivalent, no `CANDIDATE_CAP`/`VALUE_MATCH_CAP` overflow simulated). This is fine for correctness validation of the filter/dedup/match *logic*, but a dedicated overflow-path stress test (forcing single buckets or single dedup runs past those caps) is still open — natural fit for Phase 3b or Phase 7 per the plan's existing exit-criteria note.
- **(Phase 0 addition)** `cpuref_exact`/`cpuref_stats` use `qsort` and binary search rather than a radix sort or hash map; adequate up to the delivered 30M-element case (~19s) but not algorithmically optimal. Only worth revisiting if a much larger stress case is needed later.
- **(Phase 0 addition)** The harness was built and tested only under gcc 13 and clang 18 (both installed in the Phase 0 sandbox). No MSVC toolchain was available to actually compile it; the source deliberately avoids VLAs, `__int128`, and GNU extensions to stay MSVC-compatible, and the README gives an exact `cl.exe` invocation, but this has not been executed. Worth a real MSVC build the first time this harness is touched from a Windows environment.

## Handoff notes for the next phase *(update)*
- Next phase: 1. Use `OCL_GERBICZ_PHASE1_PROMPT.md`.
- The Phase 0 harness is fully self-contained and does not need CUDA, OpenCL, or GMP — it can be built and run on any machine, including one with neither GPU vendor's SDK installed. Good for CI.
- Before Phase 1 starts in earnest, the user should: (a) decide/tell the next conversation where in their repo `collharness/` should live (it wasn't asked this phase — see Inputs), (b) upload `ocl_xface.h` to the Project if it exists yet, (c) supply real `test_ad`/N/`nfs_args` values for the four cell regimes so the end-to-end recipe in the harness README can actually be exercised (this is independent of Phase 1's OpenCL-foundation work and can happen in parallel/later).
- The collcase v1 format, the exact-reference semantics, and the filter-pipeline stats semantics are all locked in as of this phase (see "Interfaces and formats" and the two new Facts entries above) — Phase 3 (collision engine port) should treat `include/collcase.h`, `cpuref_exact()`, and `cpuref_stats()` as the ground truth to validate the real OpenCL kernels against, via `collharness cmp` once the OpenCL engine can dump its own collcase-shaped output (a small dump shim, analogous to the CUDA-side task 0.9 that was skipped, is likely worth doing in Phase 3's own scope instead).
- One naming correction worth remembering: real engine hit-dump files are `test_<engine-name>.hits` (e.g. `test_gerbicz.hits`), not `test_cpu_gerbicz.hits`/`test_ocl_gerbicz.hits` as originally written in the Phase 0 kickoff prompt — harmless, already reflected in the harness README, no other document needed correcting.
