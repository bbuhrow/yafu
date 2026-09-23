# OpenCL `gpu_gerbicz` port — STATUS

This file travels from phase to phase. Attach it (with `OCL_GERBICZ_PLAN.md`) to every new conversation.
**Update protocol:** at the end of each conversation, rewrite the sections marked *(update)*. Keep it short: fold old detail into one-line summaries. Do not delete "Locked decisions" or "Facts" without recording why.

- **Last updated:** 2026-09-23, end of Phase 3
- **Current phase:** 4 (not started; prompt ready — see `OCL_GERBICZ_PHASE4_PROMPT.md`)
- **Target:** AMD GPUs without CUDA (dev card RX 6700 XT, OpenCL 2.0). Gerbicz engine only.

## Phase table *(update)*

| Phase | Title | Status | Conversation notes |
|---|---|---|---|
| 0 | Baseline and harness | **Done** | Portable C99 `collharness` tool built and tested. Its full source is now present in `/mnt/project` (`collcase.h/_io.c`, `cpuref*.c`, `cmp.c`, `gen.c`, `main.c`, etc.) — used for real in Phase 3. |
| 1 | OpenCL foundation | **Done** | `ocl_shared.h/.c` + `ocl_gerbicz_ctx.h/.c`. Three discrepancies found (header guard, missing `MAX_MEM_ALLOC_SIZE`, device-version-vs-`-cl-std` gap). MSVC unreachable; mingw-w64 substituted. |
| 2 | Primitives (scan, reduce-max, fill, radix sort) | **Done** | `ocl_primitives.h/.c` + kernels. 41/41 tests pass. Two real bugs found and fixed (missing `block_sums` write; private-vs-local array crash). |
| 3 | Collision engine (3a–3d) | **Done** | `ocl_collision.h/.c` + `ocl_collision_kernels.cl`: all 10 kernels (scatter, filter, dedup, secondary hash x2, match x3, emit x2) ported, no-sub-group baseline. **8/8 validation tests pass against the REAL Phase 0 harness** (`gen_generate`/`cpuref_exact`/`cpuref_stats`/`cmp_run`, linked directly — not reimplemented), including two cases that hit `FOUND_ARRAY_SIZE`'s 999-cap and exercise the harness's saturation-aware comparator. One real cross-queue race bug found and fixed (see Facts) — significant, affected the entire pipeline. Two new `ocl_xface.h` patch additions (`GPU_ARG_LOCAL`, `GPU_MAX_KERNEL_ARGS` 15→20). Overflow handling returns status codes, not wired to `STAGE1_OVERFLOW_SKIP` yet (Phase 5's job). MSVC still unreachable; mingw-w64 check passed. |
| 4 | Sieve (trans) kernels | Not started | |
| 5 | Host integration and registry | Not started | |
| 6 | Validation | Not started | |
| 7 | Performance and hardening | Not started | |

Status values: Not started / In progress / Done / Done with caveats.

## Locked decisions
(Full text in plan §2. No changes this phase.)
1. OpenCL 2.0 minimum; sub-groups optional and gated.
2. Static link, runtime-loaded `.cl` files, no DSO.
3. `ocl_` symbol prefix; CUDA and OpenCL builds mutually exclusive for now.
4. Base-offset kernel arguments instead of sub-buffers or SVM.
5. Shared context and program per device; per-thread queue and kernel objects.
6. No-sub-group baseline collision kernels.
7. Bit-exact first, then optimize; behavior changes logged here.
8. Registry: `STAGE1_ENGINE_OCL_GERBICZ` / `ocl_gerbicz` / `HAVE_OCL_POLY` (proposed; confirm in Phase 5).

## Facts
(Full detail and line references in plan §3 and prior STATUS revisions. New this phase:)
- **(Phase 3 addition)** `GPU_MAX_KERNEL_ARGS=15` was insufficient for `ocl_count_and_store_matched_values` (needs 16). Bumped to `20` (matching `cuda_xface.h`'s own already-bumped value) via `ocl_xface_phase3.h.patch`. This is the same discrepancy Phase 2 flagged as a *Phase 4* concern (for the fused trans+scatter kernel) — it simply arrived one phase earlier than expected. Applying this patch also fixes it for Phase 4.
- **(Phase 3 addition)** OpenCL has no `gpu_arg_type_t` case for CUDA's dynamic shared memory (`cudaFuncSetAttribute`/`MaxDynamicSharedMemorySize`) — every existing case sets an argument's *value*; a dynamically-sized `__local` kernel parameter needs its *size* set instead (`clSetKernelArg(kernel, idx, size, NULL)`). Added `GPU_ARG_LOCAL` (same patch file), with `gpu_arg_t.uint32_arg` repurposed to carry the byte size for that one slot. Needed for `filter_per_bucket_kernel`'s runtime-sized hash table (`3 * max_tsize_words * 4` bytes, exactly mirroring the CUDA original's own runtime sizing).
- **(Phase 3 addition)** `d_D`/`d_value_counts` don't need CUDA's separate `_scan`/`_offsets` buffers — CUDA only has them because `cub::DeviceScan::ExclusiveSum` always writes to a distinct output; neither pre-scan buffer is read again once the scanned version exists. This port scans both **in place** via Phase 2's `ocl_scan_exclusive_u32`, one buffer fewer per pair.
- **(Phase 3 addition — real bug, now fixed)** `ocl_primitives_t` (Phase 2) and `ocl_collision_engine_t` (Phase 3) each create their own OpenCL command queue via Phase 1's `ocl_thread_init`. Since the two constantly hand buffers back and forth (fill before scatter, reduce-max after scatter, sort then dedup, scan then scatter, ...), two independent queues gave OpenCL no ordering guarantee between them — a real race. Symptom: `candidate_count`/`dedup_count` looked plausible but `value_match_count` was 0 in every real test, traced to the secondary hash's bitmask table (`S`) reading back stale/garbage instead of freshly zeroed, because the zeroing fill and the kernel that read it ran on different queues. **Fixed** by discarding the collision engine's own queue after `ocl_thread_init` and reusing the primitives engine's queue instead (decision #5's "one queue per thread" taken literally — not one per *kind* of kernel). **Not fixed at the root**: `ocl_thread_init()` itself should accept an optional pre-existing queue rather than always creating one; the create-then-discard workaround is fine for now (two callers, one clearly "owns" the queue) but won't generalize cleanly once more concurrent engines exist. Flagged for Phase 5 or whenever this recurs.

## Inputs and blockers *(update)*

| Item | State |
|---|---|
| `ocl_xface.h` patches | Phase 1's guard patch confirmed applied to the real tree. Phase 3's `ocl_xface_phase3.h.patch` (GPU_ARG_LOCAL, arg-count bump) is new this phase, not yet applied to the real tree — **user action needed** before Phase 4 builds against the real `ocl_xface.h`. |
| Real `gpu_init`/`gpu_launch_init`/`gpu_launch_set` implementation | **Still the Phase 1 sandbox stub**, now three-going-on-four phases running. `gpu_launch_set`'s `GPU_ARG_LOCAL` case (this phase's addition) is only implemented in the stub — if/when a real `ocl_xface.c` shows up, it needs this case too. |
| Real AMD RX 6700 XT hardware | Still reachable to the user, not this sandbox. Phase 2's benchmark request is still outstanding — no update this phase. |
| File location | Still "alongside for now" per the user; Phase 3 followed suit (`harness/` also added alongside, copied from what's now in `/mnt/project`). |
| Sample YAFU invocation / real `test_ad` regimes | No update this phase — still just the one example command line + A100 timing data point from Phase 2. |
| Real MSVC (cl.exe) toolchain | Still not available, fourth phase running into this. Mingw-w64 compile-only check substituted again. |
| Project files present | Now also includes the full Phase 0 harness source (`collcase.h/_io.c`, `cpuref.h`, `cpuref_exact.c`, `cpuref_stats.c`, `cmp.h/.c`, `gen.h/.c`, `main.c`, `byteio.h`, `bucket_hash.h`, `prng.h`, `run_tests.sh`) — used directly in Phase 3's validation. |

## Interfaces and formats *(update)*
- OpenCL program/cache/thread lifecycle, base-offset convention, sub-group build-time gate: unchanged from Phase 1.
- Primitives API (fill/reduce-max/scan/sort): unchanged from Phase 2.
- **(Phase 3) Collision engine API:** `ocl_collision_init(engine, gd, prim, cl_source_path, cache_dir)` / `_free()` / `ocl_collision_run(engine, data)`. `data` (`ocl_collision_data_t`) mirrors `collision_data_t` with `cl_mem`+`cl_ulong offset` pairs for the four external buffers (`keys_in`, `data_in`, `q_batch`, `found_array`); the engine's own scratch buffers carry no offset (decision #4 rule 5 — see `ocl_collision.h`'s header comment for the exact reasoning). Returns `ocl_collision_status_t`, not `exit(-1)`, on the two overflow paths (task 3.4).
- **(Phase 3) `ocl_collision_init` requires a shared queue.** Callers must pass an `ocl_primitives_t*` that has already been initialized (`ocl_primitives_init`) on the *same* `ocl_gerbicz_device_t` — the collision engine borrows that queue (see the cross-queue-race Fact above). Passing primitives initialized against a different device, or calling `ocl_collision_run` before `ocl_primitives_init`, is undefined.
- **(Phase 3) New `-D` flags:** `COLL_BLOCK_THREADS=128`, `COLL_MAX_FILTER_ITERS=20`, `COLL_MATCH_ARENA_WIDTH=8`, `FOUND_ARRAY_SIZE=1000`, plus Phase 1's `LOG2_NUM_BUCKETS`/`BUCKET_HASH_MIX`. All baked in from the same header constants used host-side — never hardcoded twice.

## Artifacts produced *(update)*

| Phase | File | Purpose | Location |
|---|---|---|---|
| 0 | `collharness/` (now also directly in `/mnt/project`) | Portable C99 test harness | `collharness.zip` (prior), `/mnt/project` |
| 1 | `ocl_shared.*`, `ocl_gerbicz_ctx.*`, `ocl_xface.h.patch`, `ocl_xface_stub.c` | OpenCL foundation | `ocl_gerbicz_phase1.zip` |
| 2 | `ocl_primitives.*` | Fill / reduce-max / scan / radix sort | `ocl_gerbicz_phase2.zip` |
| 3 | `ocl_collision.h`, `ocl_collision.c`, `ocl_collision_kernels.cl` | The collision engine itself | `ocl_gerbicz_phase3.zip` |
| 3 | `test_ocl_collision.c` | Validation against the real Phase 0 harness (8 cases) | `ocl_gerbicz_phase3.zip` |
| 3 | `ocl_xface_phase3.h.patch` | `GPU_ARG_LOCAL` + `GPU_MAX_KERNEL_ARGS` bump | `ocl_gerbicz_phase3.zip` |
| 3 | `harness/` | Copy of the Phase 0 harness sources, for standalone building | `ocl_gerbicz_phase3.zip` |
| 3 | `Makefile`, `README.md` | Build/run instructions, discrepancies, the bug found, test table | `ocl_gerbicz_phase3.zip` |
| 3 | `OCL_GERBICZ_PHASE4_PROMPT.md` | Phase 4 kickoff prompt | delivered alongside this file |

## Open issues and risks *(update)*
- Candidate overflow risk from a lower OpenCL local-memory cap: still open, Phase 3 didn't hit it in testing (all cases stayed within `max_tsize_words`'s cap derived from `local_mem_size`) but hasn't been deliberately stress-tested either.
- Real MSVC still not run on anything in this project — fourth phase in a row.
- `clGetErrorString()`'s non-thread-safe fallback buffer (Phase 1 finding) — still relevant, still unaddressed, still only matters once Phase 5 runs multiple threads concurrently.
- The scan/radix-sort performance non-optimizations (Phase 2) are unchanged; Phase 3's own kernels have their own decision-#7 simplifications too (see README: no-sub-group compaction throughout, matching CUDA's own `SCATTER_DIRECT_ATOMIC` branch as "the model" per the plan).
- **(Phase 3 addition)** The `ocl_thread_init`-creates-its-own-queue design (Phase 1) doesn't generalize to two components that need to share state — worked around this phase (see Facts), not fixed at the root. Whoever touches `ocl_shared.c` next (likely Phase 5, wiring multiple engines/threads together for real) should consider giving `ocl_thread_init` an optional pre-existing-queue parameter instead of patching around it again.
- **(Phase 3 addition)** `GPU_ARG_LOCAL`'s host-side handling exists only in the Phase 1 sandbox stub (`ocl_xface_stub.c`) — if/when the real `gpu_launch_set` implementation is found or written, it needs this case too, or every kernel using `GPU_ARG_LOCAL` (currently just `filter_per_bucket_kernel`) will silently fail to size its local memory correctly.
- **(Phase 3 addition)** Overflow status codes (`OCL_COLLISION_CANDIDATE_OVERFLOW`/`_VALUE_MATCH_OVERFLOW`) exist but aren't wired to anything yet — no test in this phase actually exercises the overflow paths (would need a case with a pathologically dense duplicate structure at CANDIDATE_CAP scale, which wasn't attempted given the sandbox's CPU-emulated OpenCL device's speed). Worth a dedicated overflow-path test before or during Phase 5, matching Phase 0's own flagged-but-not-yet-done "dedicated overflow-path stress test" open issue.

## Handoff notes for the next phase *(update)*
- Next phase: 4 (Sieve/trans kernels). Use `OCL_GERBICZ_PHASE4_PROMPT.md`.
- Phase 4 should apply `ocl_xface_phase3.h.patch` to the real tree (on top of Phase 1's patch, already applied) before starting, and reuse `GPU_ARG_LOCAL` if the trans kernels need any dynamically-sized local memory of their own.
- Phase 4 is the first phase to actually use the arg-count bump for its originally-intended purpose (the fused trans+scatter kernel, 17 args) — the bump to 20 already covers it, no further patching needed there.
- The cross-queue-race lesson generalizes: **any new component that calls into both `ocl_primitives_t` and `ocl_collision_engine_t` (or introduces a third component) must share one queue among all of them**, not create its own. Phase 4's sieve kernels, if they end up calling any Phase 2 primitives or feeding directly into Phase 3's collision engine, need the same treatment.
- Validation approach worth keeping: Phase 3's real win was linking directly against the actual Phase 0 harness source instead of hand-rolling a comparator — caught bugs (the saturation cases) that a simpler "does it run" smoke test wouldn't have. Phase 4 should do the same wherever the harness (or a small addition to it) can express what "correct" means for the trans kernels.
