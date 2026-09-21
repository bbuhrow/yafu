# OpenCL `gpu_gerbicz` port — STATUS

This file travels from phase to phase. Attach it (with `OCL_GERBICZ_PLAN.md`) to every new conversation.
**Update protocol:** at the end of each conversation, rewrite the sections marked *(update)*. Keep it short: fold old detail into one-line summaries. Do not delete "Locked decisions" or "Facts" without recording why.

- **Last updated:** (none yet — pre-Phase 0)
- **Current phase:** 0 (not started; prompt ready)
- **Target:** AMD GPUs without CUDA (dev card RX 6700 XT, OpenCL 2.0). Gerbicz engine only.

## Phase table *(update)*

| Phase | Title | Status | Conversation notes |
|---|---|---|---|
| 0 | Baseline and harness | Not started | |
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

## Inputs and blockers *(update)*

| Item | State |
|---|---|
| `ocl_xface.h` (OpenCL `gpu_info_t`, `gpu_launch_t`, `OCL_TRY`, `clGetErrorString`) | **Missing.** Needed for Phase 1. |
| Sample YAFU invocation / N / `nfs_args` for end-to-end `test_ad` cells | Not provided yet (Phase 0 asks). |
| NVIDIA machine for CUDA golden dumps | Unknown (optional). |
| Project files present | `collision_engine.cu/.h`, `collision_bucket.h`, `cuda_intrinsics.h`, `stage1_core.cu/.h`, `stage1_engine.c/.h`, `stage1.c`, `stage1_sieve_gpu.c`, `gpu_cofactorization_cl.c/.h`, `opencl_intrinsics.cl` |

## Interfaces and formats *(update — Phase 0 fills the collcase spec)*
- **collcase v1 file format:** (TBD in Phase 0)
- **Harness tool CLI:** (TBD)
- **Canonical result form and comparison rules:** (TBD)
- **Stats-model parameters (hash-word cap override etc.):** (TBD)

## Artifacts produced *(update)*

| Phase | File | Purpose | Location |
|---|---|---|---|
| — | — | — | — |

## Open issues and risks *(update)*
- Candidate overflow on the OpenCL path may be more likely if the local-memory cap is lower than CUDA's; needs a recoverable skip (Phase 3).
- `MAX_MEM_ALLOC_SIZE` may be smaller than the two bucket arrays or the root arrays at the current sizing (Phase 5).
- AMD program compile time is long; the binary cache must be robust (Phase 1).
- Sorting fewer than the full key bits is unproven; do not do it before checking (Phase 2 / 7).

## Handoff notes for the next phase *(update)*
- Next phase: 0. Use `OCL_GERBICZ_PHASE0_PROMPT.md`.
- (The ending assistant writes here: what the next conversation needs to know that is not in the tables above.)
