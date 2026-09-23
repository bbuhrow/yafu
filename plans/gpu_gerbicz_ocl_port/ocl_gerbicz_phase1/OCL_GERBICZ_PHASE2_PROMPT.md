# Context
I'm porting the `gpu_gerbicz` CUDA stage-1 NFS polynomial-selection collision engine (YAFU/msieve) to OpenCL, targeting AMD GPUs that can't run CUDA. The work is split into phases, one per conversation. THIS conversation is **Phase 2 only** (primitives: exclusive scan, reduce-max, buffer fill, LSD radix sort — replacing CUB). Don't write collision-engine logic kernels yet (scatter/filter/dedup/match/emit are Phase 3), and don't modify the Project's source files.

Attached: `OCL_GERBICZ_PLAN.md` (full plan, locked decisions, verified facts, line-referenced source map) and `OCL_GERBICZ_STATUS.md` (progress tracker — **read its Phase 1 additions to the Facts, Interfaces, and Open issues sections before starting**, especially the base-offset argument convention and the sub-group build-time gate, both fixed in Phase 1 and binding on every kernel from here on). Also attach Phase 1's deliverable (`ocl_gerbicz_phase1.zip`): `ocl_shared.h/.c` and `ocl_gerbicz_ctx.h/.c` are the foundation this phase's kernels build on — read `ocl_gerbicz_ctx.h`'s header comment block first, it is the base-offset/sub-group contract. The source files are in `/mnt/project/`.

# Working rules
- Keep explanations brief and to the point.
- Deliver work as files, plus one zip of everything. No CUDA dependencies in anything new. Build under gcc, clang, and MSVC if reachable this session; if not, say so explicitly in STATUS (this has recurred every phase so far — Phase 0 and Phase 1 both lacked it) rather than silently skipping it again.
- Follow the base-offset kernel-argument convention and the `HAVE_SUBGROUPS` build-time gate exactly as documented in `ocl_gerbicz_ctx.h` (Phase 1) — don't redesign either.
- Use `ocl_build_program_cached()` / `ocl_thread_init()` from `ocl_shared.h` for this phase's program and kernel objects rather than writing new build/cache logic.
- If a fact in the plan turns out wrong or incomplete, say so and record the correction in STATUS (Phase 1 found three: `ocl_xface.h`'s guard, the missing `MAX_MEM_ALLOC_SIZE` field, and the device-version-vs-`-cl-std` gap — check whether any similar surprises show up here).
- Ask me, in a single message at the start, for anything missing (listed under "Inputs" below), and proceed with sensible defaults meanwhile.

# Phase 2 goal
Implement the four primitives Phase 3's collision engine needs in place of CUB (`cub::DeviceReduce::Max`, `cub::DeviceRadixSort::SortKeys`, `cub::DeviceScan::ExclusiveSum` — see `collision_engine.cu` 741–757 for exactly which CUB calls these replace): exclusive scan (up to ~4M elements), reduce-max (16384 elements), a buffer-fill kernel, and an LSD radix sort of 32- or 64-bit keys at **full key width** (not `key_bits` — Section 3 of the plan is explicit that masking to `key_bits` is unproven and is a later, separate optimization).

# Tasks (draft — confirm/adjust against the plan's §2 locked decisions and Phase 1's actual foundation code once reread)
**2.1 Confirm understanding.** Skim Phase 1's `ocl_gerbicz_ctx.h` comment block and `ocl_shared.h`'s public API; list any discrepancy or anything Phase 2 needs that Phase 1 didn't anticipate (short).

**2.2 Buffer fill.** A simple `__global` fill kernel (uint32/uint64 variants, or a single kernel parameterized appropriately) — used throughout Phase 3 to zero/reset device buffers between batches without a round-trip through `clEnqueueFillBuffer` if that's not portable enough across the target driver, or using it directly if it is (check and note which).

**2.3 Reduce-max.** Up to 16384 elements (`NUM_BUCKETS` from `collision_bucket.h` — this is exactly `d_bucket_count`'s size in `collision_engine.cu` 741–743). Block-local reduction + a second small pass, or single-pass if the work-group size covers it at this size. No sub-groups required for correctness (decision #6) — sub-group variant only as the optional fast path if time allows, gated per Phase 1's `HAVE_SUBGROUPS` convention.

**2.4 Exclusive scan.** Up to ~4M elements (`d_D`/`d_D_scan` are `MAX_DSIZE`-sized per `collision_engine.cu` 723–725; check `MAX_DSIZE`'s actual value in the source and confirm ~4M is right rather than assuming). Per the plan: block-local scan plus a block-sums pass — explicitly **no decoupled look-back** (that's the CUB single-pass trick; skip it, two-pass is fine and simpler to get right first per decision #7's "bit-exact first, then optimize").

**2.5 LSD radix sort.** Sort 32- or 64-bit keys (`d_candidate_keys` → `d_sorted_keys`, `CANDIDATE_CAP` = 2^22 elements per `collision_engine.cu` 31–37) at full key width per the plan's Section 3 signed-key fact — do not mask to `key_bits` bits. Built from 2.4's scan as the per-digit counting/scatter primitive, the standard way. Radix width (4 vs 8 bits per pass) is an open choice — pick one, document why, leave retuning for Phase 7.

**2.6 Tests.** CPU reference implementations of all four primitives (small, deliberately simple — not performance code) and a test harness comparing GPU output against them: randomized inputs plus the edge sizes (0, 1, exactly one work-group, one more than a work-group boundary, the full CANDIDATE_CAP-scale size for the sort). Run and show actual pass/fail output, not just "should work."

**2.7 Benchmarks.** Rough throughput numbers on whatever OpenCL device is available this session (see Inputs #1) at a couple of representative sizes (small, and CANDIDATE_CAP-scale for the sort) — order-of-magnitude only, not a tuning exercise (that's Phase 7).

**2.8 Wrap-up.** Update STATUS (phase table, interfaces and formats, artifacts with filenames, open issues, handoff notes) and draft the Phase 3 prompt using this file's structure.

# Source sections to read (from `/mnt/project/` and Phase 1's deliverable)
- `collision_engine.cu` 30–37 (constants: `CANDIDATE_CAP`, `MAX_C_ILOG2`, etc.) and 741–757 (the exact CUB calls being replaced, with their buffer sizes)
- `collision_bucket.h` (whole file — `NUM_BUCKETS`, already used by Phase 1's `-D` flags)
- `ocl_gerbicz_ctx.h` (whole file — base-offset convention, sub-group gate; every new kernel signature follows this)
- `ocl_shared.h` (whole file — `ocl_build_program_cached`, `ocl_thread_init`, `ocl_device_t` fields Phase 2 should query rather than re-derive, e.g. `local_mem_size`, `max_work_group_size`, `pref_wg_multiple`)
- `opencl_intrinsics.cl` (skim — existing conventions Phase 2's new `.cl` file(s) should match)

# Inputs (ask me in one message)
1. Same OpenCL device situation as Phase 1, or has anything changed (e.g. is a real AMD card reachable now)? Affects whether benchmarks (2.7) mean anything beyond "runs correctly."
2. Any preference on radix width (4-bit vs 8-bit passes) for 2.5, or leave it to my judgement (documented either way)?
3. Should Phase 2's new files live alongside Phase 1's (`ocl_shared.*`, `ocl_gerbicz_ctx.*`) in the same directory, or is there a repo location decided yet (Phase 1's Handoff notes flagged this as still unresolved)?
4. Any update on: real `test_ad`/N/`nfs_args` values (four regimes), the real `gpu_init`/`gpu_launch_init`/`gpu_launch_set` implementation question (Phase 1 used a stand-in stub), or MSVC reachability? Not required for Phase 2 itself, but flagged every phase until answered, per STATUS.

# End-of-conversation protocol
Deliver: the four primitives' source files, their CPU-reference test harness with real pass/fail output shown, a build script/Makefile addition (building on Phase 1's), a zip of everything, the updated `OCL_GERBICZ_STATUS.md`, and the Phase 3 prompt. Then give me a 3–5 line summary of what was done, what was deferred, and anything I need to do before Phase 3.
