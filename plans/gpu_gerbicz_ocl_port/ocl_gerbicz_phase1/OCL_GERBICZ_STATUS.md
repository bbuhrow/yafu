# OpenCL `gpu_gerbicz` port — STATUS

This file travels from phase to phase. Attach it (with `OCL_GERBICZ_PLAN.md`) to every new conversation.
**Update protocol:** at the end of each conversation, rewrite the sections marked *(update)*. Keep it short: fold old detail into one-line summaries. Do not delete "Locked decisions" or "Facts" without recording why.

- **Last updated:** 2026-09-23, end of Phase 1
- **Current phase:** 2 (not started; prompt ready — see `OCL_GERBICZ_PHASE2_PROMPT.md`)
- **Target:** AMD GPUs without CUDA (dev card RX 6700 XT, OpenCL 2.0). Gerbicz engine only.

## Phase table *(update)*

| Phase | Title | Status | Conversation notes |
|---|---|---|---|
| 0 | Baseline and harness | **Done** | Portable C99 `collharness` tool built and tested. See prior STATUS revision for detail. |
| 1 | OpenCL foundation | **Done** | `ocl_shared.h/.c` + `ocl_gerbicz_ctx.h/.c` built, base-offset convention and sub-group gate designed and documented, smoke test passes on gcc+clang with zero warnings, disk cache verified to actually invalidate on source change. Two real discrepancies found and fixed (see Facts). MSVC not run (no toolchain reachable); mingw-w64 compile-only cross-check substituted and passed. |
| 2 | Primitives (scan, reduce-max, fill, radix sort) | Not started | |
| 3 | Collision engine (3a–3d) | Not started | |
| 4 | Sieve (trans) kernels | Not started | |
| 5 | Host integration and registry | Not started | |
| 6 | Validation | Not started | |
| 7 | Performance and hardening | Not started | |

Status values: Not started / In progress / Done / Done with caveats.

## Locked decisions
(Full text in plan §2. No changes this phase — see Facts for two clarifications that don't require changing any locked decision, just how it's implemented.)
1. OpenCL 2.0 minimum; sub-groups optional and gated.
2. Static link, runtime-loaded `.cl` files, no DSO.
3. `ocl_` symbol prefix; CUDA and OpenCL builds mutually exclusive for now.
4. Base-offset kernel arguments instead of sub-buffers or SVM.
5. Shared context and program per device; per-thread queue and kernel objects.
6. No-sub-group baseline collision kernels.
7. Bit-exact first, then optimize; behavior changes logged here.
8. Registry: `STAGE1_ENGINE_OCL_GERBICZ` / `ocl_gerbicz` / `HAVE_OCL_POLY` (proposed; confirm in Phase 5).

## Facts (verified from source before Phase 0)
(Detail and line references in plan §3. Unchanged this phase.)
- Keys are signed (two's complement); sort full 32/64-bit width, not `key_bits`.
- Key 0 = empty slot. Packed value = `(q_index << shift) | p`; `q_index` indexes `q_batch`.
- `found_array[0].p1` counts all attempted stores; only 999 entries are stored.
- `STAGE1_OVERFLOW_SKIP` is declared but unused; CUDA engine `exit(-1)`s on candidate and value-match overflow.
- CUDA hash-word cap comes from the shared-memory opt-in limit; OpenCL local memory (64 KB on the dev card) caps at 4096 words.
- Filter output is order-independent, so a sequential CPU model can match it exactly.
- `stage1_core.cu:502` truncates a 64-bit value to `uint32` (pp64_r64). Port verbatim.
- Intrinsics port exists; `modinv64` (extra out-param) and `montmul64` (`uint w`) signatures differ.
- Cofactorization cache key is device name only; needs a source hash and driver version. **(Phase 1: fixed for the new gerbicz build path — see below. The existing `gpu_cofactorization_cl.c` cache itself is untouched; only the new shared module fixes this.)**
- **(Phase 0 addition)** `test_%s.hits` dump filenames are `test_<engine v->name>.hits`, not the literal `test_cpu_gerbicz.hits`/`test_ocl_gerbicz.hits` named in the original Phase 0 prompt. No correction needed to any locked decision.
- **(Phase 0 addition)** The engine's D/S/X secondary hash tables are an exact, collision-free lookup, not a source of approximation. Cross-checked against the deterministic edge-case suite and independent round-trip.
- **(Phase 1 addition)** `ocl_xface.h`'s feature guard (`#ifdef HAVE_OCL_BATCH_FACTOR`) would silently empty the header under a gerbicz-only build (`HAVE_OCL_POLY` alone, per locked decision #8's proposed guard name). Fixed by `ocl_xface.h.patch` (widens the guard to `defined(HAVE_OCL_BATCH_FACTOR) || defined(HAVE_OCL_POLY)`). Not applied to the real repo by this phase (project sources aren't modified in place) — apply the patch before Phase 2/3 code is compiled against the real tree.
- **(Phase 1 addition)** `gpu_info_t` (in `ocl_xface.h`) has no field for `CL_DEVICE_MAX_MEM_ALLOC_SIZE`, though the plan's Phase 1 task list calls for exposing it. `local_mem_size` and `max_work_group_size` are already covered by `shared_mem_size`/`max_threads_per_block`. Not worth a struct change — the new `ocl_device_t` (Phase 1) queries it directly from `cl_device_id` alongside the `gpu_info_t*` it wraps.
- **(Phase 1 addition)** A device's `CL_DEVICE_VERSION` is not proof a given `-cl-std=CLx.y` build option will be accepted. Observed directly in the Phase 1 sandbox: PoCL 5.0 reports device version "OpenCL 3.0" but rejects `-cl-std=CL2.0` — its `CL_DEVICE_OPENCL_C_ALL_VERSIONS` list is `{1.0, 1.1, 1.2, 3.0}` with no 2.0 entry at all, so 1.2 and 3.0 both build but 2.0 does not. `ocl_device_query()` now queries that list (falling back to the single-string `CL_DEVICE_OPENCL_C_VERSION` on older ICDs) and stores the actual best-supported value in `ocl_device_t::cl_c_std`, used by `ocl_build_program_cached()` instead of a hardcoded `"CL2.0"`. Locked decision #1's *floor* (refuse device version < 2.0) is unchanged and still enforced via `ocl_check_min_version()` against the device/platform version — this fact only changes how the build-option string is *picked* once that floor is cleared. Worth remembering for the real AMD target too: never hardcode `"-cl-std=CL2.0"` in Phase 2+ kernels' build calls.
- **(Phase 1 addition)** `gpu_init()`, `gpu_launch_init()`, and `gpu_launch_set()` are declared in `ocl_xface.h` and already called by `gpu_cofactorization_cl.c`, but no implementation (`ocl_xface.c` or equivalent) was among the files handed to Phase 1. Phase 1 shipped a sandbox-only reference implementation (`ocl_xface_stub.c`, clearly marked as such) purely so its own smoke test could build and run. **Action needed before Phase 3 depends on this:** confirm whether a real implementation already exists elsewhere in the repo; if so, reconcile it against the stub (particularly `gpu_launch_set`'s per-arg-type `clSetKernelArg` dispatch and whatever `warp_size`/preferred-multiple `gpu_init` reports) and drop the stub. If no real implementation exists yet, the stub is a reasonable starting point but wasn't written with production hardening in mind (no input validation beyond what OpenCL itself enforces).

## Inputs and blockers *(update)*

| Item | State |
|---|---|
| `ocl_xface.h` (OpenCL `gpu_info_t`, `gpu_launch_t`, `OCL_TRY`, `clGetErrorString`) | **Now attached and read** (Phase 1). Needs `ocl_xface.h.patch` applied before compiling anything against it under `HAVE_OCL_POLY` alone — see Facts. |
| Real implementation of `gpu_init`/`gpu_launch_init`/`gpu_launch_set` (declared in `ocl_xface.h`, not defined in any file handed to this phase) | **Still missing / unconfirmed.** Phase 1 used a sandbox-only stub (`ocl_xface_stub.c`) to build+run its smoke test. See Facts. |
| Sample YAFU invocation / N / `nfs_args` for end-to-end `test_ad` cells | **Still not provided.** Independent of Phase 1/2's own scope; the harness's recipe is ready once these are supplied. |
| NVIDIA machine for CUDA golden dumps | Unknown (optional, task 0.9 not attempted). |
| Real MSVC (cl.exe) toolchain | **Still not available.** Phase 0's sandbox lacked it; Phase 1's did too (Linux container, no Microsoft domains in the allowed network egress list). Substituted a compile-only mingw-w64 cross-check for Phase 1's two new files (`ocl_shared.c`, `ocl_gerbicz_ctx.c`) — passed with zero warnings, but this only catches gross portability breaks, not MSVC-specific ones. A real MSVC build is still owed. |
| Project files present | `collision_engine.cu/.h`, `collision_bucket.h`, `cuda_intrinsics.h`, `stage1_core.cu/.h`, `stage1_engine.c/.h`, `stage1.c`, `stage1_sieve_gpu.c`, `gpu_cofactorization_cl.c/.h`, `opencl_intrinsics.cl`, `ocl_xface.h` (now present) |
| Harness repo path convention in the user's actual repo | Not asked/answered this phase either; still needs a home directory decided. |
| Where Phase 1's new files (`ocl_shared.*`, `ocl_gerbicz_ctx.*`) should live in the repo relative to `collision_engine.*`/`stage1_sieve_gpu.c` | Not asked this phase — see Handoff notes. |

## Interfaces and formats *(update — Phase 1 adds the OpenCL foundation's own interfaces)*
- **collcase v1 file format:** Done in Phase 0 — see `include/collcase.h` in the delivered harness.
- **Harness tool CLI:** Done in Phase 0 — `collharness` with `gen`/`ref`/`cmp`/`suite`/`info`.
- **Canonical result form and comparison rules:** Done in Phase 0.
- **Stats-model parameters:** Done in Phase 0.
- **(Phase 1) OpenCL program/cache/thread lifecycle:** Done — see `ocl_shared.h`'s function-level doc comments. Summary: `ocl_device_init()` (query + context) → `ocl_build_program_cached()` (per program, keyed on device+driver+build-opts+source-hash) → `ocl_thread_init()` (per thread: queue + one `cl_kernel` per entry point via the existing `gpu_launch_init`). Two programs per device for gerbicz specifically (`ocl_gerbicz_device_t.program_collision` / `.program_sieve`), per the plan's "separate programs for sieve vs collision kernels" bullet.
- **(Phase 1) Base-offset kernel-argument convention:** Done — full rules written up in `ocl_gerbicz_ctx.h`'s header comment (element-based `cl_ulong`/`ulong` offsets, argument-pairing order, when an offset is/isn't needed). This is the convention every Phase 2/3/4 kernel signature must follow.
- **(Phase 1) Sub-group build-time gate:** Done — `-DHAVE_SUBGROUPS=0/1` decided from a real `cl_khr_subgroups`/`cl_intel_subgroups` extension-string query at program-build time (not a runtime kernel branch). `.cl` sources guard with `#if HAVE_SUBGROUPS`. Verified for real in the smoke test: the sandbox device reports `cl_khr_subgroups`, so `HAVE_SUBGROUPS=1` was actually compiled and exercised, not merely declared.
- **(Phase 1) `-D` flags from `collision_bucket.h`:** Done — `ocl_gerbicz_build_opts()` bakes `LOG2_NUM_BUCKETS`/`BUCKET_HASH_MIX` in as `-D` flags at build time rather than letting `.cl` sources redefine them (avoids the drift `collision_bucket.h`'s own comment warns about).

## Artifacts produced *(update)*

| Phase | File | Purpose | Location |
|---|---|---|---|
| 0 | `collharness/` (full tree), suite cases, `OCL_GERBICZ_PHASE1_PROMPT.md` | Portable C99 test harness | delivered as `collharness.zip` (prior phase) |
| 1 | `ocl_shared.h`, `ocl_shared.c` | Generic OpenCL device/context/program-cache/thread infrastructure (usable by cofactorization too) | `ocl_gerbicz_phase1.zip` |
| 1 | `ocl_gerbicz_ctx.h`, `ocl_gerbicz_ctx.c` | Gerbicz-specific layer: two-program split, base-offset convention doc, sub-group build flag, `collision_bucket.h`-derived `-D` flags | `ocl_gerbicz_phase1.zip` |
| 1 | `ocl_xface.h.patch` | Fixes the `ocl_xface.h` guard discrepancy (see Facts) | `ocl_gerbicz_phase1.zip` |
| 1 | `ocl_xface_stub.c` | **Sandbox-only** reference impl of `gpu_init`/`gpu_launch_init`/`gpu_launch_set`, not for production use without reconciling against the real implementation (see Facts) | `ocl_gerbicz_phase1.zip` |
| 1 | `ocl_gerbicz_smoke.cl`, `test_ocl_foundation.c` | Task 1.7 smoke test (build→cache→launch→verify, twice) | `ocl_gerbicz_phase1.zip` |
| 1 | `patched/ocl_xface.h`, `patched/collision_bucket.h` | Copies (one patched, one not) so `make smoke` works standalone without touching the real repo | `ocl_gerbicz_phase1.zip` |
| 1 | `Makefile`, `README.md` | Build/run instructions, discrepancies list, smoke-test results | `ocl_gerbicz_phase1.zip` |
| 1 | `OCL_GERBICZ_PHASE2_PROMPT.md` | Phase 2 kickoff prompt | delivered alongside this file |

## Open issues and risks *(update)*
- Candidate overflow on the OpenCL path may be more likely if the local-memory cap is lower than CUDA's; needs a recoverable skip (Phase 3).
- `MAX_MEM_ALLOC_SIZE` may be smaller than the two bucket arrays or the root arrays at the current sizing (Phase 5). **(Phase 1 note: now actually queried and stored in `ocl_device_t.max_mem_alloc_size`, so Phase 5's sizing check has real data to check against, not just a plan-level flag.)**
- AMD program compile time is long; the binary cache must be robust (Phase 1). **Addressed this phase** — cache key now includes a source hash, so a changed `.cl` file can't silently serve a stale binary (verified: deliberately breaking the smoke kernel's math forced a real recompile and the resulting mismatch was correctly detected, not masked by a stale cache hit).
- Sorting fewer than the full key bits is unproven; do not do it before checking (Phase 2 / 7).
- **(Phase 0 addition)** The CPU stats model assumes capacity is never exceeded; a dedicated overflow-path stress test is still open (Phase 3b/7).
- **(Phase 0 addition)** `cpuref_exact`/`cpuref_stats` use `qsort`/binary search, not a radix sort or hash map; fine up to 30M elements (~19s), revisit only if a much larger stress case is needed.
- **(Phase 0 addition)** The harness was only actually compiled under gcc 13 and clang 18; no MSVC toolchain was available to test it, though the source avoids VLAs/`__int128`/GNU extensions for MSVC-compatibility.
- **(Phase 1 addition)** Real MSVC still not run on anything in this project (Phase 0's harness or Phase 1's foundation code). Both were designed to be MSVC-compatible (no VLAs, no GNU extensions, C89-style block-scoped declarations throughout Phase 1's new files) and Phase 1's two new `.c` files passed a compile-only mingw-w64 cross-check with zero warnings, but neither substitutes for the real thing. If a Windows/MSVC machine becomes available before Phase 2, running both the harness and Phase 1's smoke test there would close this out for good instead of it recurring every phase.
- **(Phase 1 addition)** `clGetErrorString()`'s fallback path (unrecognized error codes) uses a function-local `static char buf[32]`, which is not thread-safe. No effect through Phase 1 (single-threaded), but Phase 5+ (per-thread queues actually running concurrently) should either make this thread-local or switch unrecognized codes to a `snprintf` into a caller-supplied buffer.
- **(Phase 1 addition)** `ocl_xface_stub.c` hardcodes `gi->warp_size = 32` as a placeholder inside `gpu_init()`, since the real preferred-work-group-multiple needs a built kernel to query. `ocl_device_t.pref_wg_multiple` (queried for real once a kernel is built, in `ocl_thread_init()`) is the value later phases should actually use for launch geometry — not `gpu_info_t.warp_size`, which stays a rough placeholder in the stub. If/when the stub is replaced by a real `ocl_xface.c`, check what that one does for this field instead of assuming 32.

## Handoff notes for the next phase *(update)*
- Next phase: 2 (Primitives — exclusive scan, reduce-max, buffer fill, LSD radix sort). Use `OCL_GERBICZ_PHASE2_PROMPT.md`.
- Phase 1's foundation is exercised only by a trivial vector-add smoke kernel; Phase 2 is the first phase to write kernels that do real work; it should build its primitives directly on `ocl_gerbicz_ctx.h`'s conventions (base-offset args, `ocl_build_program_cached`, `ocl_thread_init`) rather than reinventing any of that.
- Before Phase 2 (or definitely before Phase 3) starts in earnest: (a) apply `ocl_xface.h.patch` to the real project tree, (b) resolve the `gpu_init`/`gpu_launch_init`/`gpu_launch_set` implementation question (real one vs. reconciling the Phase 1 stub) — Phase 2's own primitives don't strictly need `gpu_launch_t` (they can use `clSetKernelArg` directly if simpler), but Phase 3's collision engine will, (c) decide where `ocl_shared.*`/`ocl_gerbicz_ctx.*` live in the repo relative to `collision_engine.*` (not asked this phase).
- The `collharness` Phase 0 deliverable is unrelated to Phase 1's own testing (Phase 1's smoke test is pure OpenCL pipeline, not a gerbicz-logic test) — Phase 3 is where the two connect, once the OpenCL collision engine can dump collcase-shaped output for `collharness cmp` to check.
- One naming/config reminder carried over from Phase 0: real engine hit-dump files are `test_<engine-name>.hits`, not `test_cpu_gerbicz.hits`/`test_ocl_gerbicz.hits`.
