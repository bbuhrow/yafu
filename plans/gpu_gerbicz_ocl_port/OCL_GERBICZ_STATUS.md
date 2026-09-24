# OpenCL `gpu_gerbicz` port — STATUS

This file travels from phase to phase. Attach it (with `OCL_GERBICZ_PLAN.md`) to every new conversation.
**Update protocol:** at the end of each conversation, rewrite the sections marked *(update)*. Keep it short: fold old detail into one-line summaries. Do not delete "Locked decisions" or "Facts" without recording why.

- **Last updated:** 2026-09-24, end of Phase 5
- **Current phase:** 6 (not started; prompt ready — see `OCL_GERBICZ_PHASE6_PROMPT.md`)
- **Target:** AMD GPUs without CUDA (dev card RX 6700 XT, OpenCL 2.0). Gerbicz engine only.

## Phase table *(update)*

| Phase | Title | Status | Conversation notes |
|---|---|---|---|
| 0 | Baseline and harness | **Done** | Portable C99 `collharness`; full source in `/mnt/project`. |
| 1 | OpenCL foundation | **Done** | `ocl_shared.*`, `ocl_gerbicz_ctx.*`. Phase 4 found and fixed a `modinv32` bug in `opencl_intrinsics.cl` (see Facts). Phase 5: `ocl_thread_init` gained a queue parameter; dead `program_*` fields deleted. |
| 2 | Primitives | **Done** | fill / reduce-max / scan / radix sort. 41/41 tests (Phase 2). |
| 3 | Collision engine | **Done** | 10 kernels, no-sub-group baseline; 8/8 on real 6700 XT. Cross-queue race fixed (Facts). |
| 4 | Sieve (trans) kernels | **Done** | 3 kernels, 17 args each; 10/10 on real 6700 XT. |
| 5 | Host integration and registry | **Done with caveats** | `ocl_gerbicz_driver.*` (msieve-independent core), `stage1_sieve_ocl_gerbicz.c` (msieve glue), registry row `ocl_gerbicz`, `ocl_thread_init` queue parameter, overflow → skip. **19/19 new tests + Phase 1 smoke + Phase 4's 10/10 pass — on PoCL 5.0 CPU only.** Caveats: real-hardware run outstanding; glue/registry compiled only against a hand-written `shim/stage1.h` (real `stage1.h` is not in the Project); Phase 3's own `test_ocl_collision.c` could not be re-run (not readable in the Project). gcc + clang clean (`-Wall -Wextra -pedantic`); **MSVC unreachable (sixth phase in a row)**, mingw-w64 compile-only passed for all Phase 5 sources. |
| 6 | Validation | Not started | |
| 7 | Performance and hardening | Not started | |

## Locked decisions
1. OpenCL 2.0 minimum; sub-groups optional and gated.
2. Static link, runtime-loaded `.cl` files, no DSO.
3. `ocl_` symbol prefix; CUDA and OpenCL builds mutually exclusive for now.
4. Base-offset kernel arguments instead of sub-buffers or SVM.
5. Shared context and program per device; per-thread queue and kernel objects. **(Phase 5: context shared, programs NOT — see Open issues.)**
6. No-sub-group baseline collision kernels.
7. Bit-exact first, then optimize; behavior changes logged here. **Phase 5 logged two deliberate deviations from CUDA: pp64 `root_bytes` forced to 8, batch-level overflow skip instead of `exit(-1)` (Facts).**
8. Registry: `STAGE1_ENGINE_OCL_GERBICZ` / token `ocl_gerbicz` / `HAVE_OCL_POLY` — **confirmed**, row mirrors `GPU_GERBICZ` (`{ENV_TRUNK_Q, ENV_TRUNK_P, is_gpu=1}`, `STAGE1_OVERFLOW_SKIP`, `max_threads=4`). Nothing OpenCL-specific changed it.

## Facts
(Phase 1-4 facts folded to one line each; full text in earlier revisions.)
- P1-P3: `modinv32` had a stray `q=1;` (fixed, patch already in Project copy of `opencl_intrinsics.cl`); `modinv64` takes an extra `*gcd` out-param; `montmul32_w` is deliberately called with a truncated `pp` in pp64; trans-kernel scratch is per-variant (`roots_out` for pp32_r32/pp64_r64, `p_out` for pp32_r64); `stage1_core.cu:502` truncating read ported verbatim (not fixed); collision keys are always zero-extended u64 internally; `d_D`/`d_value_counts` scanned in place; cross-queue race ⇒ **engines that hand buffers to each other must share one in-order queue**.
- **(P5) Layout at the seam, confirmed against source and by test.** Trans output index = `row*num_entries + soa_offset + p_offset + m*num_p (+ plane*num_entries*num_specialq)`. `p_out` is always u32 (`(row<<shift)|p`). Roots width per variant: PP32_R32 start u32 / out u32; **PP32_R64 start u32 / out u64 (zero-extended; `p_out` is scratch)**; PP64_R64 start u64 / out i64. `keys_in = d_root_array`, `data_in = d_p_array`, both offset 0 (the per-soa offset `j` goes to the *sieve* `p_out_off`/`roots_out_off`, in elements of that buffer's type). Same cl_mem, no copy.
- **(P5) `ocl_thread_init(..., cl_command_queue existing_queue)`**: NULL ⇒ create + own (`owns_queue=1`); non-NULL ⇒ borrow, never released by `ocl_thread_free`. `ocl_primitives_init` passes NULL; collision and sieve pass `prim->th.queue`. The three create-then-discard hacks are gone.
- **(P5) Dead fields resolved by DELETION**: `ocl_gerbicz_device_t` is now `{ ocl_device_t dev; }`; `program_collision`/`program_sieve` and the unused `ocl_gerbicz_thread_t` removed. Why not wire them up: each engine builds a different program (different sources and `-D` options, already disk-cached), so a device-level copy is a duplicate handle with a second lifetime.
- **(P5) Overflow**: `ocl_gerbicz_run_batch` maps `OCL_COLLISION_CANDIDATE_OVERFLOW` / `_VALUE_MATCH_OVERFLOW` → `OCL_GERBICZ_OVERFLOW_SKIP`; the worker drops that **batch** (not the whole cell), counts it, logs `N batch(es) skipped on collision overflow`, and continues. Tested with a real 6M-key overflow (status 1), then a normal batch on the same driver matches the reference.
- **(P5) `stage1.c` never reads `->overflow`** — the "skip" mechanism the registry comment describes did not exist anywhere; it is implemented in the worker, not the caller.
- **(P5) CUDA bug at p_max == 65536** (not fixed in CUDA, fixed in port): `pp_is_64 = p_max >= 65536`, but `root_bytes = key_bits > 32 ? 8 : 4`; at exactly `p_max == 65536` `key_bits == 32` ⇒ `root_bytes == 4` while the pp64 kernel writes 8-byte roots. Port forces `root_bytes = 8` whenever `pp_is_64`. Test `PP64_R64 at p_max=65536` fails without it (mutation-verified).
- **(P5) CUDA quirk, ported verbatim**: in the small-batch aprog logic, once `key_bits > 32` the line `num_aprog_vals = MIN(max_batch, max_batch64)/batch_size` *replaces* the `max_aprog_vals` cap — e.g. 235 aprog planes in the worker test. Looks unintended but changes nothing about correctness; flagged for Phase 7.
- **(P5) Other CUDA-driver differences (deliberate)**: `sieve_data_init` sizing uses the already-capped `num_threads` parameter (CUDA re-read uncapped `obj->num_threads`); the sieve_fb factories are reused from `d->threads[i]` (CUDA allocated its own duplicates); the dead "randomize special-q range" branch (`num_pieces > 51`, unreachable because `num_pieces ≤ 50`) is dropped; elapsed time is wall-clock (CUDA: cpu + event time); `emergency_gpu_cleanup`/atexit not ported; collision stats log is condensed (no histogram).
- **(P5) Kernel arg limit**: the trans kernels need 17 args; the **Project's** `ocl_xface.h` still had `GPU_MAX_KERNEL_ARGS 15` and no `GPU_ARG_LOCAL` — a 17-arg descriptor would overrun `arg_type[15]`. The Project's real `ocl_xface.c` also had no `GPU_ARG_LOCAL` case (it would `exit(-1)` on the filter kernel), was guarded by `HAVE_OCL_BATCH_FACTOR` only (empty TU under `HAVE_OCL_POLY`), and defined `clGetErrorString` a second time next to `ocl_shared.c`. All patched (see Artifacts). The user's mention of `gpu_shared_init` — no such function exists in that file.
- **(P5) Real `gpu_init` vs sandbox stub**: real one enumerates `CL_DEVICE_TYPE_GPU` only (PoCL CPU is invisible to it; the stub uses `ALL`), and fills `warp_size` from `CL_DEVICE_PREFERRED_WORK_GROUP_SIZE_MULTIPLE` (an OpenCL 3.0 *device* query; returns 0 / error on 2.x drivers). The driver treats `warp_size <= 0` as 32.
- **(P5) Seam-test design**: expected = `cpuref_exact(collcase built from cpuref_trans output)` (independent of both engines); buffers pre-poisoned with *structured* stale data (every slot key 0x5555, valid-looking value) so an un-cleared slot shows as a false hit; planted coprime (p1,p2) pairs constructed to collide at a chosen row/offset (both signs); a "bad q" (= pa·pb) makes zero-qq_prod slots. **Mutation-verified**: dropping the aprog pre-clear, ignoring the soa offset on `p_out` or `roots_out`, dropping the pp64 `root_bytes` fix, or skipping the found-count reset each makes tests fail.
- **(P5) Cache race**: `save_binary_to_cache` writes non-atomically; every worker thread builds its own three programs, so `ocl_gerbicz_sieve_data_init` does one serial warm-up build first (driver init + free with degree 4) so threads only ever *load* a finished cache.
- **(P5) `test_ocl_foundation.c` was a 4th `ocl_thread_init` call site** (and used `gd->program_collision`) — the plan said three. Updated.

## Inputs and blockers *(update)*

| Item | State |
|---|---|
| Patches to the real tree | **User action:** replace/patch the 12 Project files listed under Artifacts (all patches verified to reproduce the delivered files byte-for-byte from the Project copies; CRLF preserved on the CRLF files). `opencl_intrinsics.cl` needs nothing (Project copy already has the `modinv32` fix). |
| Real `gpu_init` etc. | Provided as `ocl_xface.c`; patched this phase. Sandbox tests still use `ocl_xface_stub.c`; on hardware run `make XFACE=ocl_xface.c check`. |
| `stage1.h` | **Not in the Project.** `stage1_sieve_ocl_gerbicz.c`/`stage1_engine.c` were compiled only against `shim/stage1.h` (field layouts and prototypes inferred from `stage1.c`/`stage1_sieve_gpu.c`). Expect small prototype mismatches on the first real build (`sieve_fb_*`, `poly_stats_add_qdone`, `task->d->stats`). |
| `test_ocl_collision.c` | Listed in the Project but `project_read` returns "No doc". Phase 3's 8 cases were not re-run; the collision path is exercised only through this phase's pipeline tests. Re-run on your tree. |
| Real AMD RX 6700 XT | Phase 3 and 4 confirmed. **Phase 5 has not been run on hardware.** |
| `handle_collision` | Same signature as the CUDA path (`task, threadid, p1*p2, q, uint128 root, int64 res`); worker calls it identically. |
| File location | "Alongside" (default assumed). |
| Benchmarks / sample YAFU invocation | No update. |
| MSVC | Unreachable, sixth phase running. mingw-w64 compile-only substituted (`make mingw-check CL_HEADERS_DIR=...`): all Phase 5 sources pass. |

## Interfaces and formats *(update)*
- Foundation / primitives / collision / sieve APIs: unchanged except `ocl_thread_init` (+ `existing_queue`) and `ocl_gerbicz_device_t` (`program_*` removed).
- **(P5) Driver API** (`ocl_gerbicz_driver.h`, no msieve types): `ocl_gerbicz_driver_init(drv, gd, degree, ocl_dir, cache_dir, max_entries32, max_entries64)` / `_free`; `p_reset`, `p_add(drv, p, num_roots, uint64 roots[])`, `p_start(drv, pp_is_64)` (transposes + uploads, drops arrays with `num_p*num_roots < 50`, sets `num_arrays`/`num_entries`); `q_reset`, `q_add`, `q_nextbatch`; `unused_bits`, `max_batch`, `plan_batch` (pure batch arithmetic); `run_batch(drv, batch, found_out[1000], &elapsed)` → `OK | OVERFLOW_SKIP | INVALID | CL_ERROR`. One driver = one host thread (own queue/programs/buffers).
- **(P5) Glue / registry**: `ocl_gerbicz_sieve_data_init/free`, `ocl_gerbicz_thread_data_init/free`, `stage1_specialq_ocl_gerbicz` (vtable signature). Declared in `stage1_engine.c` under `#ifdef HAVE_OCL_POLY` (no `stage1.h` edit needed). Select with `stage1_engine=ocl_gerbicz`. Runtime args: `ocldir=<dir of .cl files>` (default cwd), `ocl_cache=<dir>` (default `.`), `gpu_mem_mb=`, `collhash=`, `collstats=`, `colldebug=`.
- The `.cl` files loaded at run time: `ocl_primitives_kernels.cl`, `ocl_collision_kernels.cl`, `opencl_intrinsics.cl`, `stage1_trans_kernels.cl`.

## Artifacts produced *(update)*
Delivered as `ocl_gerbicz_phase5.zip` (full buildable tree; unchanged Project files included byte-identical so `make check` works standalone).

| File | Purpose |
|---|---|
| `ocl_gerbicz_driver.h/.c` | **New.** Core driver (p/q staging, trans launches, seam into collision, batch planning, overflow mapping). |
| `stage1_sieve_ocl_gerbicz.c` | **New.** msieve glue + registry entry points. |
| `test_ocl_gerbicz_driver.c` | **New.** 19 tests: 6 full-pipeline seam cases, overflow + recovery (2), registry/planner (9), worker-through-vtable (2). |
| `shim/stage1.h` | **New.** Compile-check stand-in for msieve's `stage1.h`. Never link a real build against it. |
| `patches/*.patch` | Diffs vs the Project files (ocl_shared.h/.c, ocl_primitives.c, ocl_collision.c, ocl_sieve.c, ocl_gerbicz_ctx.h/.c, ocl_xface.h/.c, stage1_engine.h/.c, test_ocl_foundation.c) + the patched files themselves at the top level. |
| `Makefile` | Builds/runs foundation, sieve and driver tests; `XFACE=` selects stub vs real xface; `mingw-check`. |
| `test_output.txt` | Real pass/fail output of `make check` (PoCL 5.0, CPU device). |

## Open issues and risks *(update)*
- **Real-hardware run of Phase 5** (and re-run of Phase 3's `test_ocl_collision`) — outstanding. PoCL runs the kernels but says nothing about occupancy, local-memory caps (Phase 3 risk, still open), or RDNA2 timing.
- **First real build against msieve's `stage1.h`** may need small fixes to the glue (see Inputs).
- **Decision #5 not fully met:** one program per engine *per thread* (cache-backed), not one shared program per device. Move to a shared program in Phase 7 if build time or memory matters.
- **Overflow policy is batch-level.** A pathological cell can skip many batches and lose coverage silently (only a log line). Phase 6 should measure how often it happens on real `test_ad` regimes.
- **CUDA aprog cap bypass** (Facts) — decide in Phase 7 whether to keep bit-exact behavior.
- **Blocking writes / per-soa `clFinish`** in the driver and `ocl_sieve_run` (kept for safety); pure Phase 7 material along with all other performance work.
- Still open from earlier: `clGetErrorString`'s non-thread-safe fallback buffer; MSVC never run; scan/radix-sort and no-sub-group compaction unoptimized (decision #7).
- Cache write is non-atomic (mitigated by the warm-up build; a second concurrent *process* could still race).

## Handoff notes for the next phase *(update)*
- Next phase: 6 (Validation). Use `OCL_GERBICZ_PHASE6_PROMPT.md`.
- **Before Phase 6 (user):** (1) apply the Phase 5 patches / replace the 12 files in the real tree; (2) build msieve/YAFU with `-DHAVE_OCL_POLY` against the real `stage1.h` and fix whatever the shim guessed wrong; (3) on the 6700 XT run `make XFACE=ocl_xface.c check` (Phase 1, 4 and 5 tests) and Phase 3's own `test_ocl_collision`; (4) one real run with `stage1_engine=ocl_gerbicz` and paste any errors/log.
- Validation approach worth keeping: independent references (`cpuref_trans` + `cpuref_exact`), planted known collisions, structured poison in the seam buffers, and **mutation checks** (break the code on purpose, confirm a test fails) — the aprog pre-clear was only caught once the poison was structured, not constant.
