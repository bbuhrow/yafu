# Context
I'm porting the `gpu_gerbicz` CUDA stage-1 NFS polynomial-selection collision engine (YAFU/msieve) to OpenCL, targeting AMD GPUs that can't run CUDA. The work is split into phases, one per conversation. THIS conversation is **Phase 4 only** (the sieve/trans kernels: `sieve_kernel_trans_pp32_r32`, `pp32_r64`, `pp64_r64`). Don't do host integration or registry wiring (Phase 5), and don't modify the Project's source files.

Attached: `OCL_GERBICZ_PLAN.md` and `OCL_GERBICZ_STATUS.md` (progress tracker — **read the Phase 1/2/3 additions before starting**, especially Phase 3's cross-queue-race lesson: any new component touching Phase 2's primitives or Phase 3's collision engine must share their queue, not create its own). Also attach the Phase 1/2/3 deliverables (`ocl_gerbicz_phase1.zip`, `ocl_gerbicz_phase2.zip`, `ocl_gerbicz_phase3.zip`) and apply `ocl_xface_phase3.h.patch` to the real project tree before building against it (on top of Phase 1's already-applied guard patch). The source files are in `/mnt/project/`.

# Working rules
- Keep explanations brief and to the point.
- Deliver work as files, plus one zip of everything. No CUDA dependencies in anything new. Build under gcc, clang, and MSVC if reachable this session; if not, say so explicitly in STATUS (this has now recurred in every phase so far) rather than silently skipping it again.
- Follow the base-offset kernel-argument convention and the `HAVE_SUBGROUPS` build-time gate exactly as documented in `ocl_gerbicz_ctx.h`.
- **Share one OpenCL command queue across every component this phase touches** (Phase 2's `ocl_primitives_t`, Phase 3's `ocl_collision_engine_t`, and whatever this phase adds) — Phase 3 lost real time to a cross-queue race from not doing this; don't repeat it. If this phase's sieve kernels don't call into Phase 2/3 at all, this may not apply — check first.
- **Validate against the real Phase 0 harness**, the way Phase 3 did (linking `gen.c`/`cpuref_exact.c`/`cpuref_stats.c`/`cmp.c` directly), not a hand-rolled comparator — if the harness doesn't yet have a way to express "correct sieve/trans output," say so and propose the smallest addition to it rather than building a separate, parallel validation path.
- If a fact in the plan turns out wrong or incomplete, say so and record the correction in STATUS (every phase so far has found at least one real discrepancy or bug — check whether this one does too).
- Ask me, in a single message at the start, for anything missing (listed under "Inputs" below), and proceed with sensible defaults meanwhile.

# Phase 4 goal
Port the three CUDA trans kernels (`stage1_core.cu` — `pp32_r32` lines 22-215, `pp32_r64` lines 218-415, `pp64_r64` lines 418-613) to OpenCL, using `opencl_intrinsics.cl`'s already-ported `modinv32`/`modinv64`/`montmul32`/`montmul64` (note the signature differences from CUDA already flagged in the plan's Facts — `modinv64` takes an extra out-param, `montmul64` takes a `uint w`). Base-offset arguments throughout, per the convention. Verbatim first, including the known `stage1_core.cu:502` truncation quirk in `pp64_r64` — port the bug faithfully, don't fix it (decision #7).

# Tasks (draft — confirm/adjust against the plan's §2 locked decisions and the actual source once reread)
**4.1 Confirm understanding.** Skim `stage1_core.cu` in full and `opencl_intrinsics.cl`'s actual `modinv32`/`modinv64`/`montmul32`/`montmul64` signatures (don't trust the plan's summary alone — Phases 1-3 all found real discrepancies against source that looked settled). List anything that doesn't match (short).

**4.2 Base-offset arguments for the trans kernels.** Per the plan's own note: "`roots_out` (r32, pp64_r64) or `p_out` (pp32_r64) is used as scratch for qq_prod, then overwritten. For `num_aprog_vals > 1` the host pre-clears the roots array." Work out exactly how this scratch-reuse interacts with the base-offset convention (a scratch buffer that's also an output buffer, sliced per-batch) before writing any kernel signatures — this is exactly the kind of thing that bit Phase 3 (the D/D_scan buffer-reuse question) and is worth getting right on paper first.

**4.3 Port `pp32_r32`, `pp32_r64`, `pp64_r64`.** Using the already-ported `opencl_intrinsics.cl`. Verbatim, including the `pp64_r64:502` truncation.

**4.4 `GPU_ARG_LOCAL` / arg-count check.** Phase 3 needed both a new arg type and a higher `GPU_MAX_KERNEL_ARGS` for kernels far smaller than the trans kernels look like they'll be (17 args per the plan's own note about the retired fused-kernel spike). Check whether either limit needs touching again before assuming they don't.

**4.5 Tests.** CPU reference values for a handful of real (p, q) pairs (the plan's own exit criterion: "roots match the reference for all three variants, including `num_aprog_vals > 1`"), plus edge sizes the way every prior phase's tests did (small, one-workgroup boundary, several workgroups). Validate through the real Phase 0 harness if it can express this; if not, say what the smallest addition to it would look like rather than building a parallel comparator.

**4.6 Wrap-up.** Update STATUS and draft the Phase 5 prompt using this file's structure.

# Source sections to read (from `/mnt/project/` and the Phase 1/2/3 deliverables)
- `stage1_core.cu`, whole file (pp32_r32 22-215, pp32_r64 218-415, pp64_r64 418-613)
- `opencl_intrinsics.cl`, whole file — this phase's kernels are built directly on it
- `cuda_intrinsics.h` — the CUDA original, for comparison against the OpenCL port's already-noted signature differences
- `ocl_gerbicz_ctx.h`, `ocl_primitives.h`, `ocl_collision.h` — the three things this phase may build on or need to coexist with

# Inputs (ask me in one message)
1. Any update on real AMD benchmark numbers (outstanding since Phase 2)?
2. Any update on the real `gpu_init`/`gpu_launch_init`/`gpu_launch_set` implementation (outstanding since Phase 1; this phase's `GPU_ARG_LOCAL` case is currently only implemented in the sandbox stub)?
3. Same file-location answer as before, or has a repo path been decided?
4. Does the Phase 0 harness need a small addition to validate sieve/trans output, or does one already exist that wasn't mentioned yet?

# End-of-conversation protocol
Deliver: the trans kernels' OpenCL source (+ host launch code), their test harness with real pass/fail output, a build script/Makefile addition, a zip of everything, the updated `OCL_GERBICZ_STATUS.md`, and the Phase 5 prompt. Then give me a 3–5 line summary of what was done, what was deferred, and anything I need to do before Phase 5.
