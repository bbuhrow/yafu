# Dynamic Stage-1 / Stage-2 Load Balancing — Plan Summary

Goal: replace the fixed `stage2_threads=S` split with elastic balancing —
CPU-bound stage-1 producers can temporarily act as stage-2 consumers when the
stage-2 queue backs up, while GPU-bound stage-1 engines never do (one GPU
producer is already enough to saturate many consumers). Chosen approach:
**Option A, elastic borrowing** — keep the existing two-pool structure from
`parallel_stage2.md`; add a non-blocking way for an idle-ish producer to steal
and run one queued stage-2 hit.

This also fixes a real, separately-discovered bug: the per-coefficient
deadline (`coeff_deadline`) is stored on `task_data_t` but never read anywhere
in `search_coeff_core` — stage 1 can currently sit indefinitely on one `a_d`.
Chunking stage-1 work (needed for borrowing checkpoints anyway) fixes this for
free.

---

## 1. Core mechanism

- `search_coeff_core` currently makes **one** blocking call to
  `v->specialq(...)` covering the whole `[special_q_min2, special_q_max2)`
  window for the chosen `a_d`. Nothing checks time or backlog until that call
  returns.
- Change: slice that window into resumable sub-chunks. After each chunk:
  1. Check `task->coeff_deadline` against elapsed time — **actually enforce
     it** (currently dead code).
  2. If the stage-2 backlog is high **and** the current engine is CPU-bound,
     try to steal and run one queued stage-2 task via a new non-blocking
     threadpool primitive.
- GPU engines (`gpu_cubsort`, `gpu_gerbicz`) skip step 2 entirely via a new
  `cpu_bound` flag on the engine vtable — they're latency-bound on the device,
  not consuming a CPU scheduling slot, so they should never be pulled into
  stage-2 work.
- Backlog signal reuses the **existing** `poly_stage_stats_t` counters
  (`hits - stage2_done`, already computed for the roll-up line) via a new thin
  accessor `poly_stats_backlog()`. No new source of truth.

## 2. Single-thread scaling (`-t 1`)

The scheme above does **not** degrade safely to single-threaded operation as
originally specified. Two distinct issues:

1. **Thread-count overshoot.** `num_stage2_dedicated` defaults to 1
   independent of the producer count, so `-t 1` would actually spawn 2
   threads (1 producer + 1 dedicated consumer) — violating "single-threaded"
   before borrowing logic even runs.
   - Fix: `num_stage2_dedicated = 0` whenever `obj->num_threads <= 1`,
     overriding the `stage2_threads=` parse in that case. `stage2_workers`
     sizing (`num_stage2_dedicated + obj->num_threads`) already collapses
     correctly to a single bundle with no further change.

2. **Potential deadlock at N=1.** Borrowing only drains the stage-2 queue at
   chunk *boundaries* (between `v->specialq(...)` calls). `handle_collision`
   currently submits each hit with `threadpool_add_task(pool, &tc, 1)` —
   **blocking**. If one chunk produces enough hits to fill the queue before
   the call returns, and nothing else drains it (dedicated=0, and the lone
   producer is busy inside `specialq`, not at a checkpoint), that blocking
   add waits forever.
   - Fix: make the submission itself deadlock-proof, independent of thread
     count. In `handle_collision`, switch to non-blocking submission with a
     synchronous inline fallback:
     ```c
     if (threadpool_add_task(d->stage2_threadpool, &task_control, 0) == -2) {
         /* queue (or free-task pool) full and nobody's draining fast enough —
            run it in place, on this thread, right now */
         stage1_hit_run(hit_data,
                        d->num_stage2_dedicated + threadid);
     }
     ```
     This is a general safety valve, not just a single-thread patch — at any
     thread count, if backlog spikes faster than the periodic
     borrow-checkpoints keep up, this catches it instead of relying on chunk
     size being tuned exactly right.

Net effect at N=1: no extra thread spawned; every hit either enqueues
(harmlessly — the same thread pops it moments later at its own checkpoint) or,
once the small queue fills, runs synchronously inline immediately. This
degenerates to the original sequential stage-1/stage-2 behavior through the
*same* code paths, not a special case — deadline enforcement at chunk
boundaries still applies.

## 3. Threadpool primitive (`thread.c` / `thread.h`)

Investigated `thread.c`: no non-blocking dequeue exists. `threadpool_add_task`
can be non-blocking (enqueue side only). The only dequeue path,
`threadpool_task_get_task`, is `static` and always blocks on `new_tasks_cond`.

Required addition:
- Extract the existing "run task, return it to free-queue, broadcast" block
  out of `worker_thr_routine` into a shared static helper
  `threadpool_execute_and_release()` — pure refactor, no behavior change.
- Add a new public function:
  ```c
  int threadpool_try_run_task(struct threadpool *pool, int thread_num);
  /* returns 0 = ran a task, 1 = queue empty (no-op), -1 = error */
  ```
  Non-blocking: locks `pool->mutex`, checks `threadpool_queue_is_empty`,
  dequeues if non-empty, unlocks, then executes via the shared helper on the
  *calling* thread (no new pthread spawned).
- Same lock (`pool->mutex`) and same execute/release path the normal worker
  loop already uses — no new lock ordering introduced.

## 4. Bundle-index partitioning (avoiding a collision hazard)

`stage1_hit_run` resolves per-worker private state via
`stage2_workers[thread_num]`. If a borrowing producer used its own `threadid`
directly, it could collide with a dedicated stage-2 consumer thread that
happens to share the same numeric id in its own pool.

Fix: partition the index space.
```
stage2_workers[]  size = num_stage2_dedicated + num_threads (producers)
  [0 .. num_stage2_dedicated)                    -> dedicated pool workers
  [num_stage2_dedicated .. +num_threads)          -> producer i borrows index
                                                      num_stage2_dedicated + i
```
Each producer always borrows into the same private slot, so no contention
even though the *task* came from the shared queue.

## 5. Code changes

### `thread.c` / `thread.h`
- New `threadpool_execute_and_release()` (extracted, static).
- New public `threadpool_try_run_task(pool, thread_num)`.

### `stage1.c` — `search_coeff_core`
- Wrap the single `v->specialq(...)` call in a chunked loop over
  `[special_q_min2, special_q_max2)`, chunk size `window / QCHUNK_DIVISOR`
  capped at `QCHUNK_MAX` (both new tunables, need benchmarking).
- `stage1_engine_cell_fits` check moves inside the loop (per-chunk clamp
  instead of once for the whole window).
- After each chunk: enforce `task->coeff_deadline` (fixes the dead-code bug);
  then, if `v->cpu_bound` and `poly_stats_backlog(d->stats) >
  STAGE2_HIGH_WATERMARK`, call
  `threadpool_try_run_task(d->stage2_threadpool, d->num_stage2_dedicated + threadid)`.

### `search_coeffs` (same file)
- Parse `coeff_deadline=` from `obj->nfs_args` instead of the hardcoded
  `deadline_per_coeff = 8640000`.

### `stage1_engine.h` (vtable)
- Add `cpu_bound` flag so GPU engines opt out of borrowing.

### `poly_stats.{c,h}`
- Add `poly_stats_backlog()` accessor wrapping the existing mutex-guarded
  `hits - stage2_done` read. No new state.

### `poly_skew.c`
- `stage1_sieve_data_t` gains `num_stage2_dedicated` (copied through next to
  `sieve_data.stats = data->stats;`).
- `num_stage2_dedicated` computation: 0 whenever `obj->num_threads <= 1`
  (overrides `stage2_threads=` in that case), otherwise the existing
  `stage2_threads=` parse.
- `stage2_workers` array sized `num_stage2_dedicated + obj->num_threads`
  instead of just the dedicated count; `build_stage2_workers()` itself is
  unchanged (it already just builds N identical independent bundles). At
  `num_stage2_dedicated = 0` this collapses to exactly one bundle.
- Dedicated `d->stage2_threadpool` still sized `num_stage2_dedicated` only —
  only the bundle *array* grows to cover borrowers. At `num_stage2_dedicated
  = 0`, no dedicated pool thread is spawned at all.
- Teardown loop bound changes from `num_stage2_dedicated` to
  `num_stage2_dedicated + obj->num_threads`.

### `stage1.c` — `handle_collision`
- Stage-2 submission changes from blocking (`threadpool_add_task(pool, &tc,
  1)`) to non-blocking with a synchronous inline fallback: on a `-2` return
  (queue/free-task pool full), run the hit immediately on the calling thread
  via `stage1_hit_run(hit_data, d->num_stage2_dedicated + threadid)` instead
  of enqueuing. Prevents deadlock at `num_stage2_dedicated = 0` and acts as a
  general backlog safety valve at any thread count.

## 6. What doesn't change

- The `file_lock`/`save_poly` shared-write serialization from
  `parallel_stage2.md` §2c is untouched — borrowed tasks run the identical
  `stage1_hit_run` → callback chain, just dispatched from a different call
  site (a producer's chunk-boundary checkpoint instead of a dedicated
  consumer's blocking `get_task`).
- `build_stage2_workers()` logic, `stage2_worker_t` layout, and the
  `sizeopt_callback`/`rootopt_callback` bodies are all unchanged.

---

## TODO / open items

1. **Confirm the actual msieve timing call** to replace the placeholder
   `read_clock()` in the chunk loop (likely in `common/util.c`).
2. **Pick `QCHUNK_DIVISOR` / `QCHUNK_MAX`** — needs benchmarking; too fine
   over-slices small windows, too coarse defeats the checkpoint's purpose.
3. **Pick `STAGE2_HIGH_WATERMARK`** — threshold on `poly_stats_backlog()` that
   triggers borrowing.
4. **Verify `stage1_engine_cell_fits` is safe to call per-chunk** — confirm
   its cost is negligible and its clamp behavior is correct when applied
   repeatedly to sub-windows rather than once for the whole range.
5. **Add `cpu_bound` to all four engine vtables** (`cpu_gerbicz`,
   `cpu_hashtable` → true; `gpu_cubsort`, `gpu_gerbicz` → false).
6. **`best_saved_combined_e` semantics** (carried over from
   `parallel_stage2.md`'s own verify list, item 3) — if it's a save-gate
   rather than a stat, per-worker copies (now including borrower slots) mean
   more independent local-bests; decide whether to leave as-is or move under
   the shared lock.
7. **Write/port the `s2_pattern_test.c`-style standalone concurrency test**
   for the borrowing path specifically: producer threads intermixing chunk
   work and borrowed stage-2 tasks, dedicated consumers running concurrently,
   verifying no bundle-index collisions and no interleaved `.p` file writes.
8. **Decide dedicated-pool floor for multi-threaded runs** — keep
   `num_stage2_dedicated` small (e.g. 1–2) so the queue never fully starves
   during a moment when every producer happens to be mid-chunk and hasn't
   reached a checkpoint yet. (At `-t 1` this floor is forced to 0 regardless
   — see §2.)
9. **Verify the synchronous fallback in `handle_collision` doesn't re-enter
   any lock the caller already holds** — `stage1_hit_run` ultimately reaches
   `rootopt_callback`, which takes `file_lock`; confirm `handle_collision`'s
   call stack holds no conflicting lock when this inline path fires.
10. **Test `-t 1` end-to-end** — confirm no thread-count regression and no
    stall, ideally as part of the same standalone concurrency test in item 7
    (add a dedicated single-thread case).