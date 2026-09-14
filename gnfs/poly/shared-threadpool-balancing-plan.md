# Dynamic Stage-1 / Stage-2 Load Balancing — Shared Threadpool (Option B)

Alternative to `dynamic-stage2-balancing-plan.md`'s elastic-borrowing approach
(Option A). Where Option A keeps separate stage-1/stage-2 pools and adds a
non-blocking steal primitive, Option B removes the pool split entirely: one
`N`-worker threadpool, one priority-ordered task queue, stage-2 tasks always
win. GPU stage-1 engines are unaffected either way — they run outside this
pool and feed hits in via the existing `handle_collision` path.

---

## 1. Task model

One queue, two task kinds:

```c
enum task_kind { TASK_STAGE1_CHUNK, TASK_STAGE2_HIT };
```

- `TASK_STAGE1_CHUNK` — search one q-sub-window for one `a_d`, then re-enqueue
  its own continuation (or move to the next `a_d`). Replaces the persistent
  per-producer loop in `search_coeff_core`.
- `TASK_STAGE2_HIT` — same payload as today's `stage1_hit_data_t`.

**Priority rule:** stage-2 tasks always dequeue before stage-1 tasks. This is
the entire balancing mechanism — no watermark, no borrow trigger, no fairness/
aging logic needed, because stage-2 priority is the explicit design goal
(finish polynomials), not a load-balance compromise.

## 2. Stage-1 progression is queue-resident, not thread-resident

`a_d` ownership becomes state that travels through the queue instead of
living on a producer thread's stack:

```c
typedef struct {
    mp_t     ad;
    uint32_t q_cursor;      /* resume point in [special_q_min2, special_q_max2) */
    uint64_t deadline_end;  /* absolute wall-clock deadline, set once per a_d */
} stage1_chunk_state_t;
```

Worker executing a `TASK_STAGE1_CHUNK`:
1. Run `v->specialq()` over one chunk starting at `q_cursor`.
2. Any collisions → `handle_collision` → enqueue `TASK_STAGE2_HIT` (high
   priority).
3. Check `deadline_end` (wall-clock — see §5); if exceeded or window
   exhausted, call `find_next_ad` and reset cursor; else advance `q_cursor`.
4. Re-enqueue itself as a new `TASK_STAGE1_CHUNK` (low priority) rather than
   looping in place, so a pending stage-2 task can jump ahead on the next
   dequeue.

`num_stage1_slots` in-flight continuations are seeded once at startup and
perpetuate 1:1 — no separate slot-throttling logic (see §4: not needed, since
workers should never idle).

## 3. Stage-2 bundles: bundle pool, not per-thread bundle

`stage2_worker_t` (private `sizeopt_data`/`rootopt_data`) can no longer be
indexed by thread identity — any worker can run any `TASK_STAGE2_HIT` at any
time. Replace with a fixed-size free-list:

```c
stage2_worker_t *bundle_pool[N];  /* N = total worker count; mutex-guarded free-list */
```

Sized exactly `N` — a worker can only run one task at a time, so concurrent
stage-2 executions are bounded by `N` automatically. No tunable cap, no
exhaustion/requeue handling: worst case is all `N` workers on stage-2 at once,
which is exactly `N` bundles in use, never more.

`file_lock`/stats-mutex usage inside `rootopt_callback` is unchanged — bundles
stay private per in-flight task; only the pool they're drawn from is now
shared/dynamic rather than a static array indexed by `threadid`.

## 4. Scheduling policy

- Dequeue order: check stage-2 lane first; if empty, check stage-1 lane.
- Workers never idle: stage-1 continuations are always available to fall back
  to (until the search space is genuinely exhausted), so no
  backlog-triggered throttle on stage-1 slot count is needed — holding a
  worker back from stage-1 would only produce idling, never a benefit.
- No fairness/aging rule for stage-1: starvation under sustained stage-2
  volume is working as intended, not a bug.

## 5. Timing: wall-clock, not `get_cpu_time()`

msieve's existing `get_cpu_time()` (via `getrusage(RUSAGE_SELF, ...)` on
non-Windows) measures **process-wide** user CPU time, not per-thread or
per-task. In the shared pool, an `a_d`'s successive chunks can land on
different threads interleaved with unrelated stage-1 and stage-2 work, so a
`get_cpu_time()` delta doesn't isolate this `a_d`'s time. Even
`RUSAGE_THREAD` wouldn't fix it, since it's current-thread-only and chunks
move between threads.

Use wall-clock elapsed time instead:

```c
double get_wall_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}
```

Add as a new function alongside `get_cpu_time()` in `common/util.c` (don't
change `get_cpu_time()`'s semantics — other callers may depend on it). Set
`deadline_end = get_wall_time() + coeff_deadline` once per `a_d`; compare
against `get_wall_time()` at each chunk boundary.

## 6. `-t 1` degeneration

Free by construction: one worker, one queue, same priority rule. It runs any
pending `TASK_STAGE2_HIT` before the next `TASK_STAGE1_CHUNK`. No dedicated-
pool-size-0 special case, no inline-fallback deadlock avoidance, no bundle-
index partitioning — all of which Option A needs specifically to handle
`N=1` safely.

## 7. Early-abort interaction

Once a stage-2 task reports Murphy-E above the abort threshold: stop pulling
new tasks from the stage-1 lane (just don't dequeue/re-enqueue further
`TASK_STAGE1_CHUNK`s) and let already-in-flight stage-2 tasks finish. Because
stage-1 state lives entirely in queued task payloads rather than producer
thread stacks, "drop remaining stage-1 work" is just clearing that lane —
no thread cancellation or join logic, and no need to interrupt a chunk mid-
`specialq()`; in-flight chunks simply finish and their continuations aren't
re-queued.

## 8. GPU path

Unchanged. GPU engines run on their own thread(s) outside the shared pool and
call `handle_collision` → enqueues `TASK_STAGE2_HIT` into the same shared
queue. The `cpu_bound` vtable flag from Option A becomes unnecessary — GPU
engines were never pool workers under this design in the first place.

## 9. Comparison to Option A (elastic borrowing)

| | Option A | Option B (this plan) |
|---|---|---|
| Producer identity | fixed per OS thread | none — `a_d` state travels as task payload |
| Stage-2 concurrency bound | `stage2_threads` (thread count) | fixed at `N`, no tunable |
| Balancing mechanism | watermark-triggered steal at chunk checkpoints | implicit — priority dequeue |
| New primitive needed | `threadpool_try_run_task` | priority-aware dequeue |
| `N=1` handling | special-cased (dedicated=0, inline fallback, index partitioning) | free by construction |
| Early abort | not addressed | trivial — clear stage-1 lane, let stage-2 drain |
| Code churn | moderate, additive | larger — bundle-per-thread → bundle pool; producer loop → continuation tasks |

---

## Proposed code changes

### `thread.c` / `thread.h`
- Replace single FIFO with a two-lane (or single-list-with-priority-field)
  queue; dequeue always checks the stage-2 lane first.
- `worker_thr_routine` becomes a flat "dequeue highest-priority task, run it,
  repeat" loop — no fixed producer/consumer roles.
- No `threadpool_try_run_task` needed (that primitive is Option-A-specific).

### `stage1.c` — `search_coeff_core`
- Replace the single blocking `v->specialq()` call over the whole window with
  chunked execution driven by `TASK_STAGE1_CHUNK` payloads (`stage1_chunk_state_t`
  above): one chunk per task invocation, self-re-enqueuing continuation.
- `stage1_engine_cell_fits` check moves inside the per-chunk path.
- Deadline check per chunk boundary using `get_wall_time()` (§5) against
  `deadline_end`, computed once when an `a_d` is picked.

### `search_coeffs` (same file)
- Parse `coeff_deadline=` from `obj->nfs_args` (replacing the hardcoded
  `deadline_per_coeff = 8640000`); seed `num_stage1_slots` continuation tasks
  at startup.

### `common/util.c`
- Add `get_wall_time()` (`CLOCK_MONOTONIC`-based), alongside the existing
  `get_cpu_time()` — don't modify `get_cpu_time()` itself.

### `poly_skew.c`
- Replace `stage2_workers[]` static per-thread array with `bundle_pool[N]`
  free-list (§3): acquire on `TASK_STAGE2_HIT` dequeue, release after
  `stage1_hit_run`/callback chain completes.
- `stage2_threads=` parameter is superseded — no dedicated-pool sizing logic
  needed; pool is just `N` workers total.

### `stage1_engine.h` (vtable)
- `cpu_bound` flag from Option A is not needed here — omit.

### `poly_stats.{c,h}`
- `poly_stats_backlog()` accessor is optional under this design (no
  watermark-driven throttle), but may still be useful for the roll-up
  display line. Keep if cheap; not required for the balancing logic itself.

### Early-abort path (wherever it currently lives / will be added per phase 3)
- On abort trigger: stop dequeuing/re-enqueuing `TASK_STAGE1_CHUNK`; let
  in-flight `TASK_STAGE2_HIT` executions drain naturally.

---

## TODO / open items

1. Priority queue implementation — two lanes (separate lists + condvars) vs.
   one list with a priority field and a single condvar.
2. Confirm `stage1_engine_cell_fits` is safe/cheap to call per-chunk rather
   than once per full window (carried over from Option A's equivalent item).
3. Pick chunk-size tunables analogous to Option A's `QCHUNK_DIVISOR`/
   `QCHUNK_MAX` — needs benchmarking; same tradeoff (too fine over-slices,
   too coarse defeats checkpointing).
4. Decide whether "clear the stage-1 lane" on early abort should also try to
   interrupt an in-flight chunk mid-`specialq()`, or just prevent its
   continuation from being re-queued (current lean: the latter — simpler,
   and the in-flight chunk finishes quickly anyway since it's bounded by
   chunk size).
5. Standalone concurrency test: sustained stage-2 volume with stage-2-lane
   priority, confirm bundle pool never exceeds `N` concurrent uses, confirm
   `-t 1` degenerates correctly, confirm early-abort drain behavior.
6. Wire `get_wall_time()` into msieve's build (confirm `clock_gettime`
   availability/portability matches existing YAFU platform support matrix).
