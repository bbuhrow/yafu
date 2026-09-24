
# Context
I'm porting the `gpu_gerbicz` CUDA stage-1 NFS polynomial-selection collision engine (YAFU/msieve) to OpenCL, targeting AMD GPUs that can't run CUDA. The work is split into phases, one per conversation. THIS conversation is **Phase 6 only** (validation). Don't start Phase 7 (performance/hardening). Don't modify the Project's source files directly — only the working tree.

Attached: `OCL_GERBICZ_PLAN.md` (read its Phase 6 section and reconcile this prompt with it; this draft was written from STATUS alone) and `OCL_GERBICZ_STATUS.md` (**read the Phase 5 additions first**, especially: the seam layout fact, the two deliberate deviations from CUDA — pp64 `root_bytes` forced to 8 and batch-level overflow skip — the CUDA aprog-cap quirk, and the "Inputs and blockers" table). Also attach `ocl_gerbicz_phase5.zip` (full buildable tree). Source files are in `/mnt/project/`.

# Working rules
- Keep explanations brief and to the point.
- Deliver work as files, plus one zip of everything. No CUDA dependencies in anything new. Build under gcc, clang, and MSVC if reachable; otherwise say so explicitly in STATUS (recurring — six phases).
- Test against **independent** references (as Phases 4-5 did) and **mutation-check** any new test (break the code on purpose; confirm the test fails).
- If a plan fact turns out wrong or incomplete, record it in STATUS.
- Ask me in one message at the start for missing inputs (below) and proceed with defaults.

# Phase 6 goal
Establish that the OpenCL engine (`stage1_engine=ocl_gerbicz`) finds the same stage-1 hits as the reference engines on real workloads, and characterize its failure modes, on real AMD hardware.

# Tasks (draft — confirm against the plan's Phase 6 section)
**6.1 Confirm state.** Report which of my pre-Phase-6 actions are done (patches applied, real-`stage1.h` build, `make XFACE=ocl_xface.c check` on the 6700 XT, Phase 3's `test_ocl_collision` re-run, one real `ocl_gerbicz` run). Fix whatever the first real build against msieve's `stage1.h` exposes (the glue was only shim-compiled).
**6.2 Hardware regression.** Interpret the real-hardware output of Phases 1/3/4/5 tests; diagnose any failure that PoCL didn't show (local-memory caps, work-group limits, `warp_size` = 0 on 2.x drivers, `CL_KERNEL_WORK_GROUP_SIZE` below 128 for the trans kernels).
**6.3 End-to-end equivalence.** Using the registry's `test_ad=`/`test_pmin`/`test_pmax`/`test_qmin`/`test_qmax` test mode (writes `test_<engine>.hits`), run `cpu_gerbicz` (and `gpu_gerbicz` if a CUDA box is available) vs `ocl_gerbicz` on the same cells; compare hit sets (canonicalize; account for the 999-entry found cap and the deliberate deviations). Cover: pp32 and pp64 (`p_max` around 65536, incl. exactly 65536), aprog>1 (small special-q counts), degrees 4/5/6, the trivial q=1 row, multi-thread (`num_threads` 2-4).
**6.4 Overflow and saturation.** Find real cells that hit `OCL_GERBICZ_OVERFLOW_SKIP` and found-array saturation; measure frequency and lost coverage; decide whether batch-level skip needs a retry with half the batch.
**6.5 Robustness.** Soft/hard stop mid-batch, deadline handling (wall-clock vs CUDA's cpu+event time), repeated init/free, cache warm-up race with 4 threads on a cold cache, `gpu_mem_mb=` small values.
**6.6 Wrap-up.** Update STATUS, draft the Phase 7 prompt. List anything Phase 7 should optimize, from measurements, not guesses.

# Inputs (ask me in one message)
1. Results of the pre-Phase-6 actions (STATUS "Handoff notes").
2. A sample YAFU invocation and real `test_ad` regimes to validate on.
3. Is a CUDA machine available for a `gpu_gerbicz` cross-check?
4. Any real benchmark numbers yet?

# End-of-conversation protocol
Deliver: new/changed test scripts and sources, real pass/fail output, a zip, the updated `OCL_GERBICZ_STATUS.md`, and the Phase 7 prompt. Then a 3-5 line summary: done / deferred / what I must do before Phase 7.
