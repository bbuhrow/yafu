# Context
I'm porting the `gpu_gerbicz` CUDA stage-1 NFS polynomial-selection collision engine (YAFU/msieve) to OpenCL, targeting AMD GPUs that can't run CUDA. The work is split into phases, one per conversation. THIS conversation is **Phase 3 only** (the collision engine itself: scatter, per-bucket filter, sort/dedup, secondary hash, match, emit — plan sub-phases 3a–3d). Don't port the sieve/trans kernels (Phase 4) or do host integration (Phase 5), and don't modify the Project's source files.

Attached: `OCL_GERBICZ_PLAN.md` and `OCL_GERBICZ_STATUS.md` (progress tracker — **read the Phase 1 and Phase 2 additions to Facts, Interfaces, and Open issues before starting**; Phase 2 in particular found two real bugs from insufficiently-scaled testing — read what triggered each before writing this phase's own tests). Also attach Phase 1's and Phase 2's deliverables (`ocl_gerbicz_phase1.zip`, `ocl_gerbicz_phase2.zip`): this phase calls `ocl_scan_exclusive_u32`/`ocl_reduce_max_u32`/`ocl_radix_sort_u64`/`ocl_fill_u32`/`ocl_fill_u64` directly rather than reimplementing any of them — read `ocl_primitives.h`'s header comment first (in particular, the case for a single `ulong`-only radix sort — you'll want the same reasoning for anything else in this phase that seems to call for a 32-/64-bit split). The source files are in `/mnt/project/`, and now include `cuda_xface.h` (the CUDA original) alongside `ocl_xface.h`.

# Working rules
- Keep explanations brief and to the point.
- Deliver work as files, plus one zip of everything. No CUDA dependencies in anything new. Build under gcc, clang, and MSVC if reachable this session; if not, say so explicitly in STATUS (this has now recurred in Phases 0, 1, and 2) rather than silently skipping it again.
- Follow the base-offset kernel-argument convention and the `HAVE_SUBGROUPS` build-time gate exactly as documented in `ocl_gerbicz_ctx.h`.
- Use Phase 2's primitives (`ocl_primitives.h`) for scan/reduce-max/sort/fill rather than writing new versions of any of them.
- **Test at the scale that actually exercises multi-workgroup code paths, not just a single large smoke run.** Phase 2's two real bugs (a kernel not writing what its own design comment said it would; a private array that only crashed at `num_wg >= ~25`) both hid behind small/single-workgroup test cases. Use the same "small edge case → one-workgroup boundary → several workgroups → real CANDIDATE_CAP scale" progression for every new kernel here, not just the ones that seem likely to need it.
- If a fact in the plan turns out wrong or incomplete, say so and record the correction in STATUS (Phase 1 found three, Phase 2 found two more plus two real bugs — check whether anything similar shows up here).
- Ask me, in a single message at the start, for anything missing (listed under "Inputs" below), and proceed with sensible defaults meanwhile.

# Phase 3 goal
Port `collision_engine.cu`'s actual collision-search logic to OpenCL: scatter roots into buckets, per-bucket hash filter, sort+dedup the candidate keys (via Phase 2's radix sort), secondary hash for exact matching, and emit found pairs — bit-exact against the CPU reference, per decision #7 ("bit-exact first, then optimize"; do not attempt sub-group or other optimizations this phase beyond what decision #6 already requires as baseline).

# Tasks (draft — confirm/adjust against the plan's §2 locked decisions, and Phase 1/2's actual code, once reread)
**3.1 Confirm understanding.** Skim `collision_engine.cu` in full (not just the previously-read line ranges) and list any discrepancy against the plan, Phase 1, or Phase 2 (short).

**3.2 Data layout.** Design the OpenCL equivalent of `struct collision_engine` (`collision_engine.cu` 600-805): device buffers (`cl_mem` instead of raw pointers), `ensure_capacity`-style growth, `init`/`free_device`/`free` lifecycle mirroring the CUDA original and Phase 1's `ocl_gerbicz_thread_t` shape. Reuse `ocl_primitives_t`'s scratch-buffer pattern as a model, but this engine's own buffers (`d_arr_a`/`d_arr_b`/`d_candidate_keys`/etc.) are logically distinct from Phase 2's internal scratch space — don't share cl_mem objects between the two.

**3.3a Scatter + reduce-max.** Port `scatter_roots_kernel` (111-153) without the `__match_any_sync` sub-group path (decision #6: no-sub-group baseline; `SCATTER_DIRECT_ATOMIC`-equivalent global-atomic-only version). Overflow handling: `atomicOr` on an overflow flag → OpenCL atomic equivalent. Reduce-max over `bucket_count` via Phase 2's `ocl_reduce_max_u32`. Host grow-and-retry loop (886-926) if overflow flag set.

**3.3b Per-bucket filter.** Port `filter_per_bucket_kernel` (155-336): 128-thread groups (`BLOCK_THREADS`), the local hash table sized from `LOCAL_MEM_SIZE` (Phase 1's `ocl_device_t.local_mem_size`, not a hardcoded CUDA shared-memory constant), survivor compaction via a local atomic counter, `MAX_FILTER_ITERS` iteration cap. Collect the `filter_iters_hist`/bucket-size stats (`collision_data_t`, `collision_engine.h`) only when `collect_stats != 0`, matching the CUDA original's "zero atomic ops when off" behavior.

**3.3c Sort, dedup, secondary hash.** Port `dedup`/`count_secondary`/`scatter_secondary`/`find_candidate_slot` (338-407) and the D/S/X secondary hash tables (Phase 0's Fact: these are an exact, collision-free lookup, not an approximation — this phase should preserve that property, not approximate it). Sort candidates via Phase 2's `ocl_radix_sort_u64` (not a new sort). `MAX_DSIZE`/`MAX_SSIZE` buffer sizing per the corrected Phase 2 fact (`MAX_DSIZE` ≈ 1.05M, not ~4M — don't reintroduce that confusion here).

**3.3d Match + emit.** Port the three match kernels (409-524) and `emit_found_kernel`/`emit_found_arena_kernel` (526-598). `found_array[0].p1` atomic-counts all attempts; only `FOUND_ARRAY_SIZE-1` (998) entries actually stored (`stage1_core.h`) — preserve this saturation behavior exactly, and make it observable in the test harness (i.e., test a case that actually saturates, the way Phase 0's harness README documented for its own `med_hashmode1` case).

**3.4 Overflow handling.** Locked decision #7 + the plan's open issue: replace CUDA's `exit(-1)` on candidate overflow (1005-1013) and value-match overflow (1108-1113) with a return code the (future, Phase 5) driver can treat as a skip — design the return-code contract now even though Phase 5 wires it up, and register it under the existing but currently-unused `STAGE1_OVERFLOW_SKIP` (`stage1_engine.h`/`.c`) rather than inventing a new mechanism.

**3.5 Tests.** Build a small collcase-shaped dump path (Phase 0's STATUS already flagged this as "likely worth doing in Phase 3's own scope") so `collharness cmp` can validate this engine's output against the Phase 0 harness's exact/stats CPU references — this is the real validation, not just a hand-rolled comparison. Use the harness's own suite cases (tiny/skew, and regenerate med/large via `collharness suite` per its README) rather than inventing new ones. Follow the "small → boundary → several workgroups → real scale" progression from the working rules above for anything Phase 0's suite doesn't already cover at the right granularity (e.g. the local hash table's `MAX_FILTER_ITERS` cap, the `FOUND_ARRAY_SIZE` saturation case).

**3.6 Wrap-up.** Update STATUS and draft the Phase 4 prompt using this file's structure.

# Source sections to read (from `/mnt/project/` and the Phase 1/2 deliverables)
- `collision_engine.cu`, whole file (1165 lines) — this phase ports nearly all of it; the plan's line-referenced map (§4) is a starting index, not a substitute for reading it
- `collision_bucket.h`, `collision_engine.h`, `stage1_core.h` (`FOUND_ARRAY_SIZE`, `found_t`, `specialq_t`) — whole files
- `stage1_engine.h`/`.c` — `STAGE1_OVERFLOW_SKIP` and the registry shape task 3.4/decision #8 will eventually plug into
- `ocl_primitives.h` (whole file) and `ocl_gerbicz_ctx.h` (whole file) — the two things this phase builds directly on
- The Phase 0 harness's `include/collcase.h`, `cpuref_exact()`, `cpuref_stats()` — ground truth for 3.5

# Inputs (ask me in one message)
1. Any update on real AMD benchmark numbers for Phase 2's primitives (still outstanding — see STATUS)?
2. Any update on the real `gpu_init`/`gpu_launch_init`/`gpu_launch_set` implementation, or should Phase 3 keep building against the Phase 1 stub?
3. Is the Phase 0 `collharness` tool available to attach this session (needed for 3.5's real validation), or does it need regenerating from `OCL_GERBICZ_PHASE0_PROMPT.md`'s original spec?
4. Same file-location answer as before (alongside Phase 1/2's files), or has a repo path been decided?

# End-of-conversation protocol
Deliver: the collision engine's OpenCL source (kernels + host struct/lifecycle), the collcase dump shim and its `collharness cmp` validation output, a build script/Makefile addition, a zip of everything, the updated `OCL_GERBICZ_STATUS.md`, and the Phase 4 prompt. Then give me a 3–5 line summary of what was done, what was deferred, and anything I need to do before Phase 4.
