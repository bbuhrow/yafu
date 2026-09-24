# OpenCL `gpu_gerbicz` port — STATUS

This file travels from phase to phase. Attach it (with `OCL_GERBICZ_PLAN.md`) to every new conversation.
**Update protocol:** at the end of each conversation, rewrite the sections marked *(update)*. Keep it short: fold old detail into one-line summaries. Do not delete "Locked decisions" or "Facts" without recording why.

- **Last updated:** 2026-09-23, end of Phase 2
- **Current phase:** 3 (not started; prompt ready — see `OCL_GERBICZ_PHASE3_PROMPT.md`)
- **Target:** AMD GPUs without CUDA (dev card RX 6700 XT, OpenCL 2.0). Gerbicz engine only.

## Phase table *(update)*

| Phase | Title | Status | Conversation notes |
|---|---|---|---|
| 0 | Baseline and harness | **Done** | Portable C99 `collharness` tool built and tested. |
| 1 | OpenCL foundation | **Done** | `ocl_shared.h/.c` + `ocl_gerbicz_ctx.h/.c`; base-offset convention and sub-group gate designed; smoke test passes gcc+clang, cache-invalidation verified. Three discrepancies found (header guard, missing `MAX_MEM_ALLOC_SIZE`, device-version-vs-`-cl-std` gap). MSVC unreachable; mingw-w64 compile-only check substituted. |
| 2 | Primitives (scan, reduce-max, fill, radix sort) | **Done** | `ocl_primitives.h/.c` + `ocl_primitives_kernels.cl`: fill, reduce-max, exclusive scan (two-level/three-kernel), LSD radix sort over `ulong` keys (8×8-bit passes). 41/41 correctness tests pass (gcc+clang) up to and including `CANDIDATE_CAP`-scale (4,194,304) sort and the real `VALUE_MATCH_CAP+1` (4,194,305) scan. Two real bugs found and fixed during testing (see Facts). Benchmarks run only on the sandbox's CPU-emulated OpenCL device — real AMD numbers still needed from the user (see Inputs). MSVC still unreachable; mingw-w64 compile-only check passed. |
| 3 | Collision engine (3a–3d) | Not started | |
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

## Facts (verified from source before Phase 0)
(Detail and line references in plan §3. Unchanged facts omitted here — see prior STATUS revisions for the full Phase 0/1 list. New this phase:)
- Keys are signed (two's complement); sort full 32/64-bit width, not `key_bits`. Key 0 = empty slot. `found_array[0].p1` counts all attempted stores; only 999 entries stored. Filter output is order-independent. `stage1_core.cu:502` truncates a 64-bit value to `uint32` (port verbatim). Intrinsics `modinv64`/`montmul64` signatures differ from CUDA. (Full detail: plan §3.)
- **(Phase 1)** `ocl_xface.h`'s guard, missing `MAX_MEM_ALLOC_SIZE`, device-version-vs-`-cl-std` gap — see prior STATUS revision. `ocl_xface.h.patch` has since been applied to the real project tree (confirmed: `/mnt/project/ocl_xface.h` now shows the widened guard).
- **(Phase 2 addition)** `MAX_DSIZE` (`collision_engine.cu:36`, `(1u<<20)+64`) is **1,048,640** (~1.05M), not ~4M. The plan/Phase-2-prompt's "~4M" scan-size figure refers to a *different* `cub::DeviceScan::ExclusiveSum` call three lines later (`collision_engine.cu:754-756`, over `d_value_counts`/`d_value_offsets`, size `VALUE_MATCH_CAP+1` = `CANDIDATE_CAP+1` = **4,194,305**). Both sizes are now covered by Phase 2's test suite. No locked decision affected — just a disambiguation between two source lines.
- **(Phase 2 addition)** Only one radix sort width is needed. `scatter_roots_kernel` (`collision_engine.cu:122-126`) zero-extends a 4-byte root into a `uint64 k` (sign-extension happens later, only on emit); `d_candidate_keys`/`d_sorted_keys` are `uint64*` unconditionally regardless of `root_bytes`. So "sort the full 32- or 64-bit key" (plan §3) is about never truncating to `key_bits`, not about needing two sort implementations — there's only one storage width at the point sorting happens. Phase 2 implements a single `ocl_radix_sort_u64`, plain unsigned (no sign-bit-flip), matching `cub::DeviceRadixSort::SortKeys`'s own default behavior on a `uint64_t*` buffer.
- **(Phase 2 addition)** `GPU_MAX_KERNEL_ARGS` is `15` in `ocl_xface.h` but was bumped to `20` in `cuda_xface.h` on 2026-05-25 for the fused trans+scatter kernel (17 args) — the OpenCL header wasn't updated to match. Not a problem for Phase 2 (largest kernel here has 9 args) but **Phase 4 will need this bumped** before porting the fused trans+scatter kernel. No patch written this phase (out of Phase 2's scope) — flagging so Phase 4 doesn't discover it mid-port.
- **(Phase 2 addition — real bug, now fixed)** `ocl_scan_local_u32`'s `block_sums[wg] = acc;` write was present in the design comment but missing from the actual kernel body. Silent with exactly one scan workgroup (`n <= 4096`, where the "block sum" is scanned down to 0 regardless), so it passed at `n=4096` and failed immediately at `n=4097`. Found by seeding the block-sums buffer with a sentinel before the kernel ran and confirming the kernel never touched it. Fixed; the multi-workgroup path is now explicitly covered by the test suite (`n=4097, 8193, 100000`, plus both ~1M/~4M real sizes).
- **(Phase 2 addition — real bug, now fixed)** The radix scatter kernel's per-workgroup bin-cursor array (`uint cursor[256]`) was originally declared `private`. Even though only `lid==0` uses it, PoCL's CPU backend segfaulted with a corrupt stack once `num_wg` got large enough (first reproduced at `n=100000`, `num_wg=25`; small sizes with `num_wg<25` didn't trigger it) — apparently sizing private memory per-work-item across the whole workgroup regardless of which lanes' control flow actually reaches the declaration. Moved to `__local` (one workgroup-shared allocation) and the crash disappeared, including at full `CANDIDATE_CAP` scale (`num_wg=1024`). Worth remembering generally: a moderately large array declared inside an early-`return`-gated branch should probably be `__local`, not `private`, unless every lane truly needs its own copy.

## Inputs and blockers *(update)*

| Item | State |
|---|---|
| `ocl_xface.h` guard patch | **Applied to the real tree** (confirmed via `/mnt/project/ocl_xface.h` this phase). |
| Real implementation of `gpu_init`/`gpu_launch_init`/`gpu_launch_set` | **Still the Phase 1 sandbox stub.** User attached `cuda_xface.h` (the CUDA original) this phase as a reference for the arg-count bump (see Facts), but no real `ocl_xface.c`/equivalent was supplied — `ocl_xface_stub.c` is still what's actually running. |
| Real AMD RX 6700 XT hardware | **Reachable to the user, not to this sandbox.** This container still has no `/dev/dri`/ROCm and no route to the user's machine. Phase 2's benchmarks (see Phase 2's own table) are the same PoCL CPU-emulated device as Phase 1 — order-of-magnitude sanity only. **Action needed:** user runs `make smoke` from Phase 2's deliverable on the real card and pastes back the benchmark lines. |
| Radix width preference | User: no preference — Phase 2 picked 8-bit passes (8 passes for a 64-bit key), documented as an open Phase 7 tuning target. |
| File location | User: alongside Phase 1's files for now (same directory), will distribute into the real tree themselves. Phase 2 followed suit. |
| Sample YAFU invocation / N / `nfs_args` | **Partially answered.** User gave one example command line and a real timing data point: `./yafu "nfs(rsa(420))" -v -np -nfs_stage1_args "stage1_engine=gpu_gerbicz stage2_threads=15"` searched `a_d` 120–5100 in ~1200s on one A100 GPU. This is useful context (order of magnitude for a real run) but is not yet the four `test_ad`/`test_pmin`/`test_pmax`/`test_qmin`/`test_qmax` regime values the harness's end-to-end recipe needs (small/<480 bits, larger, pp32, pp64) — still open. |
| NVIDIA machine for CUDA golden dumps | Unknown (optional, task 0.9 not attempted). |
| Real MSVC (cl.exe) toolchain | **Still not available**, third phase running into this. Mingw-w64 compile-only cross-check substituted again (passed, zero warnings) — still not equivalent to real MSVC. |
| Project files present | Now also includes `cuda_xface.h` (CUDA original, reference-only) and `ocl_xface_stub.c` (the Phase 1 stub, apparently committed into the tree by the user). |
| Where Phase 1/2's files should ultimately live in the repo | Still not decided (user: "I will distribute to the real tree" — no path given yet). |

## Interfaces and formats *(update)*
- OpenCL program/cache/thread lifecycle, base-offset convention, sub-group build-time gate, `collision_bucket.h`-derived `-D` flags: unchanged from Phase 1, see prior STATUS revision.
- **(Phase 2) Primitives API:** `ocl_primitives_init()`/`_free()` (one program, all 8 kernels); `ocl_fill_u32/u64`, `ocl_reduce_max_u32`, `ocl_scan_exclusive_u32` (in-place, `n <= OCL_PRIM_SCAN_MAX_N` ≈ 16.7M), `ocl_radix_sort_u64` (`keys_in`/`keys_out` must be distinct allocations; internally manages its own scratch/ping-pong buffers). All follow the base-offset convention (`cl_ulong` element offset immediately after each buffer arg). See `ocl_primitives.h`'s header comment for the full design rationale, including why there's only one sort width.
- **(Phase 2) Chunking constants:** `OCL_PRIM_LOCAL_SIZE=256`, `OCL_PRIM_ELEMS_PER_THREAD=16`, `OCL_PRIM_ELEMS_PER_WG=4096`, `OCL_PRIM_RADIX_BITS=8`, `OCL_PRIM_RADIX_BINS=256`, `OCL_PRIM_RADIX_PASSES_U64=8`. Baked into the `.cl` source as `-D` flags from the same header macros — never hardcoded a second time.

## Artifacts produced *(update)*

| Phase | File | Purpose | Location |
|---|---|---|---|
| 0 | `collharness/`, suite cases | Portable C99 test harness | `collharness.zip` (prior phase) |
| 1 | `ocl_shared.h/.c`, `ocl_gerbicz_ctx.h/.c`, `ocl_xface.h.patch`, `ocl_xface_stub.c`, smoke test, `patched/` | OpenCL foundation | `ocl_gerbicz_phase1.zip` |
| 2 | `ocl_primitives.h`, `ocl_primitives.c`, `ocl_primitives_kernels.cl` | Fill / reduce-max / exclusive scan / LSD radix sort (`ulong` keys) | `ocl_gerbicz_phase2.zip` |
| 2 | `test_ocl_primitives.c` | CPU references + correctness tests (41 cases) + rough benchmarks | `ocl_gerbicz_phase2.zip` |
| 2 | `Makefile`, `README.md` | Build/run instructions, discrepancies, the two bugs found, benchmark caveat | `ocl_gerbicz_phase2.zip` |
| 2 | `OCL_GERBICZ_PHASE3_PROMPT.md` | Phase 3 kickoff prompt | delivered alongside this file |

## Open issues and risks *(update)*
- Candidate overflow on the OpenCL path may be more likely if the local-memory cap is lower than CUDA's; needs a recoverable skip (Phase 3).
- `MAX_MEM_ALLOC_SIZE` now actually queried (`ocl_device_t.max_mem_alloc_size`, Phase 1) — Phase 5's sizing check has real data to check against.
- AMD program compile time is long; cache now keyed on source hash (Phase 1), verified to invalidate correctly.
- Sorting fewer than the full key bits is unproven; do not do it before checking (Phase 2 / 7). **(Phase 2 note: not attempted — full 64-bit width used throughout, per the zero-extension fact above.)**
- Real MSVC still not run on anything in this project (harness, Phase 1 foundation, or Phase 2 primitives) — recurring every phase. Would be good to close out for good if a Windows machine becomes reachable.
- `clGetErrorString()`'s fallback path uses a non-thread-safe `static char buf[32]` (Phase 1 finding) — still relevant, becomes a real risk once Phase 5's per-thread queues run concurrently.
- **(Phase 2 addition)** The scan primitive's two-level design caps at `OCL_PRIM_ELEMS_PER_WG²` ≈ 16.7M elements (block-sums array must itself fit one workgroup's scan). Comfortably covers every current use (`CANDIDATE_CAP+1` ≈ 4.19M is the largest), but if a future phase needs a bigger scan, this cap needs lifting (a third scan level, or decoupled look-back) — `ocl_scan_exclusive_u32` returns `-1` rather than silently misbehaving if `n` exceeds it, so this will surface as an explicit error, not silent corruption.
- **(Phase 2 addition)** Performance is intentionally not addressed: the scan's block-sums combine and the radix scatter are both single-thread-serial within their scope (see Facts and README for why, and decision #7 for the rationale). Real throughput numbers await the user running the benchmark on real AMD hardware (see Inputs) — Phase 7 is where actual tuning happens.
- **(Phase 2 addition)** `GPU_MAX_KERNEL_ARGS=15` vs CUDA's `20` (see Facts) — Phase 4 needs to bump this (in `ocl_xface.h`, via a small patch like Phase 1's guard fix) before porting the fused trans+scatter kernel's larger arg list.

## Handoff notes for the next phase *(update)*
- Next phase: 3 (Collision engine — scatter, per-bucket filter, sort/dedup/secondary-hash, match, emit). Use `OCL_GERBICZ_PHASE3_PROMPT.md`.
- Phase 3 is the first phase where correctness actually matters end-to-end against real collision-engine semantics — it should validate against the Phase 0 `collharness` suite (via a small OpenCL-side dump shim, as Phase 0's STATUS already flagged as "likely worth doing in Phase 3's own scope"), not just against Phase 2's primitive-level CPU references.
- Phase 3 will call `ocl_scan_exclusive_u32`/`ocl_reduce_max_u32`/`ocl_radix_sort_u64`/`ocl_fill_*` directly rather than reimplementing any of them — the collision engine's own `d_bucket_count` (16384, reduce-max), `d_candidate_keys`→`d_sorted_keys` (up to `CANDIDATE_CAP`, radix sort), and the two `ExclusiveSum` call sites (`MAX_DSIZE` and `VALUE_MATCH_CAP+1`, both now confirmed and tested) map directly onto Phase 2's primitives.
- Before Phase 3 (or at latest Phase 4) starts in earnest: (a) resolve the real `gpu_init`/`gpu_launch_init`/`gpu_launch_set` question (still open, three phases running — user has only supplied `cuda_xface.h` as CUDA-side reference so far, not an OpenCL implementation), (b) get real AMD benchmark numbers for at least Phase 2's primitives (user has the hardware; this sandbox does not), (c) decide the eventual repo location for `ocl_shared.*`/`ocl_gerbicz_ctx.*`/`ocl_primitives.*` (still "alongside for now" per the user).
- The two bugs found this phase (missing `block_sums` write, private-vs-local array) are both the kind that only show up at a scale smaller test runs don't reach — Phase 3's own collision-engine kernels should get the same "small edge case, then boundary case, then real scale" test progression Phase 2 used, not just a single large-scale smoke test.
