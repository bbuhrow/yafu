# Phase 1 kickoff prompt

**How to use:** start a new conversation in the same Claude Project (so the source files are available under `/mnt/project/`), attach `OCL_GERBICZ_PLAN.md`, `OCL_GERBICZ_STATUS.md`, and `ocl_xface.h` (if it now exists), and paste everything inside the fence below.

````
# Context
I'm porting the `gpu_gerbicz` CUDA stage-1 NFS polynomial-selection collision engine (YAFU/msieve) to OpenCL, targeting AMD GPUs that can't run CUDA. The work is split into phases, one per conversation. THIS conversation is **Phase 1 only** (OpenCL foundation). Don't write collision-engine kernels yet, and don't modify the Project's source files.

Attached: `OCL_GERBICZ_PLAN.md` (full plan, locked decisions, verified facts, line-referenced source map) and `OCL_GERBICZ_STATUS.md` (progress tracker — read its Phase 0 additions to the Facts, Interfaces, and Open issues sections before starting). The source files are in `/mnt/project/`. The Phase 0 test harness (`collharness/`) is attached or already in the repo at <PATH FROM STATUS.md's Inputs section — fill in once known>; its `include/collcase.h`, `cpuref_exact()`, and `cpuref_stats()` are the ground truth this and later phases validate against.

# Working rules
- Keep explanations brief and to the point.
- Deliver work as files, plus one zip of everything. No CUDA dependencies in anything new. Build under gcc, clang, and (if available this session) MSVC; if MSVC still isn't available, say so explicitly in STATUS rather than silently skipping it again.
- If a fact in the plan turns out wrong or incomplete, say so and record the correction in STATUS.
- Ask me, in a single message at the start, for anything missing (listed under "Inputs" below), and proceed with sensible defaults meanwhile.

# Phase 1 goal
Stand up the OpenCL context/program/build infrastructure (`ocl_` prefix, per the locked decisions) that every later phase's kernels will run inside, without yet porting any actual collision-engine logic.

# Tasks (draft — confirm/adjust against the plan's §2 locked decisions and this phase's actual `ocl_xface.h` once seen)
**1.1 Confirm understanding.** Skim `ocl_xface.h` (now attached) and the source sections below; list any discrepancy with the plan or with Phase 0's STATUS additions (short).

**1.2 Context/program lifecycle.** One `cl_context` and one built `cl_program` per device (per locked decision 5), runtime-loaded `.cl` source files (locked decision 2 — no DSO). Design the load/build/cache path: where compiled binaries are cached, how a cache is invalidated (device name alone is known to be insufficient per Phase 0's cofactorization-cache Fact — use a source hash + driver version instead, as already flagged).

**1.3 Per-thread queue/kernel objects.** Per locked decision 5: shared context/program, but per-thread `cl_command_queue` and `cl_kernel` instances. Design the struct(s) and lifecycle (init/free) mirroring `collision_engine`'s own init/free/ensure_capacity shape in `collision_engine.cu`, so Phase 3 can reuse the pattern directly.

**1.4 Error handling.** An `OCL_TRY`-style macro (per `ocl_xface.h` if it defines one; otherwise propose one matching `CUDA_TRY`'s shape in `collision_engine.cu`) plus `clGetErrorString`.

**1.5 Base-offset kernel-argument convention.** Per locked decision 4 (no sub-buffers, no SVM): settle the exact calling convention for kernels that need to operate on a sub-range of a larger buffer (base-offset arguments), since this affects every kernel signature from Phase 2 onward. Write it down precisely enough that Phase 2/3 don't have to re-derive it.

**1.6 Sub-group gating.** Per locked decision 1 (OpenCL 2.0 minimum, sub-groups optional and gated): design the capability check and the compile-time or runtime switch between sub-group and no-sub-group kernel variants (Phase 3's baseline collision kernels are no-sub-group per locked decision 6, but the gating mechanism itself belongs here).

**1.7 Build and smoke-test.** A minimal `.cl` kernel (e.g. a trivial vector-add or the eventual scan/reduce-max primitives' skeleton, whichever is more useful to prove the pipeline) compiled and run through the new infrastructure, on whatever OpenCL device is available in this sandbox (note in STATUS if none is — the actual AMD RX 6700 XT target hardware is presumably not available here either, matching Phase 0's "no NVIDIA machine" situation on the CUDA side).

**1.8 Wrap-up.** Update STATUS (phase table, interfaces and formats, artifacts with filenames, open issues, handoff notes) and draft the Phase 2 prompt using this file's structure.

# Source sections to read (from `/mnt/project/`)
- `collision_engine.h` (whole file — the `collision_data_t` DSO interface this phase's infrastructure must eventually plug into)
- `collision_engine.cu` 600–805 (`collision_engine` struct: init/free/ensure_capacity/free_device — the lifecycle shape to mirror)
- `gpu_cofactorization_cl.c/.h` (whole files — the existing CUDA-to-OpenCL precedent for context/program setup conventions already used elsewhere in this codebase; follow its style where it doesn't conflict with a locked decision)
- `opencl_intrinsics.cl` (skim — existing OpenCL kernel-source conventions, includes, macros already in use)
- `ocl_xface.h` (whole file, once attached)

# Inputs (ask me in one message)
1. Where should `collharness/` (the Phase 0 test harness) live in the repo? (Needed so later phases can reference it by a stable path instead of "wherever it was attached this session".)
2. Is `ocl_xface.h` attached to this conversation? If not, what should Phase 1 assume/stub for its contents (`gpu_info_t`, `gpu_launch_t`, `OCL_TRY`, `clGetErrorString`)?
3. Is there an OpenCL-capable device (any vendor) available in this session's sandbox to actually build/run kernels against, or should Phase 1's "smoke test" be build-only (compile the `.cl` source, skip execution) until real AMD hardware is available?
4. Any update on real `test_ad`/N/`nfs_args` values for the four cell regimes (small/<480 bits, larger, pp32, pp64)? Not required for Phase 1 itself, but flagged every phase per STATUS until answered.

# End-of-conversation protocol
Deliver: the OpenCL foundation source files, a build script/Makefile addition, a zip of everything, the updated `OCL_GERBICZ_STATUS.md`, and the Phase 2 prompt. Then give me a 3–5 line summary of what was done, what was deferred, and anything I need to do before Phase 2.
````
