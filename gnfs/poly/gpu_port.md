# GPU engine port — kyleaskine's `stage1_sieve_gpu.c` → the stage-1 contract

kyleaskine's file already unifies **#3 (CUB sort)** and **#4 (Gerbicz DSO)** via
`use_collision_engine`. The port reshapes the *host boundary* to the contract and
keeps everything below it — the `.cu` kernels, the CUB/Gerbicz DSO interfaces,
the SoA machinery, and (for the first cut) the engine's own factory. All of it
lives under `#ifdef HAVE_CUDA`.

**Can't be built here** (no `nvcc`, no CUDA/DSO headers), so this is an
apply-and-build guide keyed to kyleaskine's line numbers, like the earlier
wiring docs. Verify against your CUDA build; §8 is the checklist.

A design note up front: I'm keeping the engine's **own** `sieve_p_fb`/`sieve_q_fb`
for this first working port (self-contained, known-good population), *not*
switching to the driver's factory. Population parity with #1/#2 is a clean
follow-up (§7); getting a correct GPU path first is worth more than parity.

---

## 1. Registry change (small — needed so the GPU knows #3 vs #4)

`gpu_data_init` must pick the DSO (`load_sort_engine` for #3 vs
`load_collision_engine` for #4) at init, but the contract's
`sieve_data_init(obj, num_threads)` doesn't carry the engine id. Add it.

**`stage1_engine.h`** — the vtable slot gains the id:
```c
	void * (*sieve_data_init)(msieve_obj *obj, uint32 num_threads,
				stage1_engine_id id);   /* +id */
```
**`stage1.c`** — the driver's `stage1_sieve_data_init` passes it:
```c
	if (v->sieve_data_init)
		d->hw_data = v->sieve_data_init(obj, num_threads, v->id);   /* +v->id */
```
The CPU rows keep `sieve_data_init = NULL`, so nothing else changes. Both GPU
rows (`gpu_cubsort`, `gpu_gerbicz`) point at the *same* `gpu_sieve_data_init`
and the *same* `stage1_specialq_gpu`; the id is the only thing that differs.

---

## 2. The contract functions (in `stage1_sieve_gpu.c`)

### `gpu_sieve_data_init` — from `gpu_data_init` (line 1778)

Keep the body **except**: it now takes `num_threads` + `id`, drops `poly`, drops
the `MIN(4, …)` cap (the driver already capped by the envelope `max_threads=4`),
sets `use_collision_engine` from the **id**, and creates **no pools**.
```c
void *
gpu_sieve_data_init(msieve_obj *obj, uint32 num_threads, stage1_engine_id id)
{
	device_data_t *d;
	gpu_config_t gpu_config;
	/* … gpu_init, GPU-exists checks, alloc d, gpu_info, logprintf … (unchanged) */

	read_collision_engine_args(obj, d);                  /* keep: collhash/stats/debug */
	d->use_collision_engine = (id == STAGE1_ENGINE_GPU_GERBICZ);   /* registry wins */
	if (d->use_collision_engine)
		load_collision_engine(obj, d);
	else
		load_sort_engine(obj, d);

	/* … gpu_mem + max_sort_entries32/64 sizing … (unchanged; uses sizeof, not degree) */

	d->num_threads = num_threads;                        /* param, already capped */
	d->max_sort_entries32 /= num_threads;
	d->max_sort_entries64 /= num_threads;
	d->threads = (device_thread_data_t *)xcalloc(num_threads,
					sizeof(device_thread_data_t));

	/* DELETE the two threadpool_init blocks and the thread_control setup —
	   the polysize driver owns both pools. DELETE `d->poly = poly;`. */

	return d;                                            /* → stored in the driver's d->hw_data */
}
```
In `read_collision_engine_args` (line 1461), **stop setting the enable bit** —
keep only the `collhash=`/`collstats=`/`colldebug=` reads. The id drives the
engine choice now; leave `collengine=` parsing in as a harmless no-op or remove
it.

### `gpu_sieve_data_free(void *hw_data)`

Take whatever `gpu_data_init`'s teardown did, minus the pools: free the
per-thread GPU state, the loaded engine (`*_free` + `dlclose`), `gpu_info`, and
`d` itself. (If teardown was inline in the old `sieve_lattice_gpu`, lift it here.)

### `gpu_thread_data_init(void *data, int threadid)` — retarget (line 1646)

Currently `data` is the GPU `device_data_t`. Under the contract `data` is the
**polysize** `stage1_sieve_data_t *d`, so reach the GPU context through it and
publish the per-thread slot back:
```c
void
gpu_thread_data_init(void *data, int threadid)
{
	stage1_sieve_data_t *d = (stage1_sieve_data_t *)data;
	device_data_t *gd = (device_data_t *)d->hw_data;
	device_thread_data_t *t = gd->threads + threadid;

	/* … the existing per-thread setup: CUcontext, module, cuMemAlloc of the
	   p/root/found arrays, stream, per-thread sort/collision engine, events,
	   and the engine's own fb (t->sieve_p_fb = sieve_fb_alloc(); …) … */

	d->threads[threadid].hw_thread_data = t;   /* so the worker finds it by threadid */
}
```
`gpu_thread_data_free` mirrors it (reach `gd` via `d->hw_data`).

### `stage1_specialq_gpu` — from `sieve_lattice_gpu_core` (line 1338)

This is the per-`a_d` worker. Reshape the entry to the contract signature, force
the ranges, pull the GPU thread state by `threadid`:
```c
void
stage1_specialq_gpu(task_data_t *task, uint32 threadid,
		uint64 special_q_min, uint64 special_q_max,
		uint32 p_min, uint32 p_max)
{
	poly_coeff_t *c = task->c;
	msieve_obj *obj = task->obj;
	device_data_t *gd = (device_data_t *)task->d->hw_data;
	device_thread_data_t *t = gd->threads + threadid;

	/* keep sieve_lattice_gpu_core's body, but:
	   - use special_q_min/max, p_min/p_max instead of deriving them
	   - sieve_fb_init(t->sieve_p_fb, c, …) as before (engine's own fb)
	   - pass `task` down to sieve_specialq → check_found_array (see §3) */

	sieve_specialq(obj, c, gd, t, task, special_q_min, special_q_max,
			p_min, p_max, (double)task->coeff_deadline);
}
```
`sieve_specialq` (line 968) and `handle_special_q_batch` (line 732) are internal
— thread a `task_data_t *task` parameter through them so it reaches
`check_found_array`. The `use_collision_engine` branch inside
`handle_special_q_batch` (line 843) is untouched: it's exactly the #3/#4 split,
now driven by the id set in §2.

**Delete** the old top-level `sieve_lattice_gpu` (line 1973) and any
`poly_stage1_run`/coefficient-dispatch loop in this file — the polysize driver
owns that now.

---

## 3. `check_found_array` — retarget to the shared handler (lines 637–730)

Add a `task_data_t *task` parameter (thread it in from `stage1_specialq_gpu`).
The clamp/read stays; the per-hit block collapses to one call:
```c
		if (coeff <= c->coeff_max)
			handle_collision(task, (uint64)p1 * p2, (uint64)q,
					gpu_promote128(qroot), offset);
```
**Delete** the `status`/`if (status==1){ build hit_data; submit to
d->stage2_threadpool }`/`else if (status==2) crap++` machinery — the shared
`handle_collision` does the Kleinjung check, submits to the driver's stage-2
pool, and bumps `hits`. Also delete `stage1_hit_data_t`, `stage1_hit_run`,
`stage1_hit_free`.

`c->found_count++` (the per-coefficient tally) has no post-check status to gate
on anymore; either drop it or bump it right before the call as an "attempts"
counter — the authoritative counts are the global stats now.

Add the promote helper near the top of the file (same as #1's `ht_promote128`):
```c
static uint128
gpu_promote128(uint64 r)
{
	uint128 u;
	u.w[0] = (uint32)r; u.w[1] = (uint32)(r >> 32); u.w[2] = 0; u.w[3] = 0;
	return u;
}
```

---

## 4. The two fb callbacks — ABI adapt (same as #1)

The polysize factory calls back with `(uint64 p, mpz_t *roots)`, not
`(uint32 p, uint64 *roots)`. Convert with `gmp2uint64`.

**`store_specialq` (line 428)** — signature + the one root read:
```c
store_specialq(uint64 q, uint32 num_roots, mpz_t *roots, void *extra)
{
	uint64 q2 = (uint64)q * q;
	…
	s->p    = (uint32)q;
	s->pp   = q2;
	s->root = gmp2uint64(roots[i]);        /* was roots[i] as uint64 */
}
```
**`store_p_soa` (line 314)** — signature `(uint64 p, uint32 num_roots,
mpz_t *roots, void *extra)`; cast `p` to `uint32` where it's stored into the SoA;
and where the body walks `roots` as a `uint64*` (`rs[m]`), read
`gmp2uint64(roots[<running index>])` instead. It's the same shape as #1's
`store_p_packed` conversion — convert at the point of use, keep the SoA packing
untouched.

---

## 5. `device_data_t` / removals

- `device_data_t` (line 533): delete `gpu_threadpool` and `stage2_threadpool`.
- `stage1_hit_data_t` / `stage1_hit_run` / `stage1_hit_free` (605–635): delete.
- `read_collision_engine_args`: keep, minus the enable bit (§2).
- `g_active_gpu_device` + `emergency_gpu_cleanup`: **keep** — the signal-handler
  GPU reset is harmless and useful with one instance (the override guarantees one).

---

## 6. The `#3/#4` selection, end to end

`stage1_engine=gpu_cubsort` → registry id `GPU_CUBSORT` → `use_collision_engine=0`
→ `load_sort_engine` → CUB path in `handle_special_q_batch`.
`stage1_engine=gpu_gerbicz` → id `GPU_GERBICZ` → `use_collision_engine=1` →
`load_collision_engine` → Gerbicz DSO path. One worker, one file, the id decides.
`collengine=` is retired as the selector; `collhash=`/`collstats=`/`colldebug=`
remain as Gerbicz tuning knobs (ride `nfs_stage1_args`).

---

## 7. Follow-up: population parity (use the driver's factory)

To make the GPU search the identical `(p,q)` population #1/#2/#3 do — drop the
engine's own `sieve_fb_alloc`/`sieve_fb_init`, and in `stage1_specialq_gpu` use
`task->d->threads[threadid].sieve_p_fb` / `sieve_q_fb` (already reset+init'd by
the driver's `search_coeff_core`). The catch is that the driver currently inits
the factory with CPU-side bounds; the GPU's `sieve_fb_init` (line 1410) may use
different ones, so parity means unifying those bounds in the driver. Worth doing
once the engine runs, not before.

---

## 8. Verify on build

1. `HAVE_CUDA` defined; file compiled by `nvcc` (or the host `.c` by the C
   compiler with the CUDA driver API headers), the `.cu` unchanged.
2. Both DSOs resolvable at run time: `cub/sort_engine.{so,dll}` (#3) and the
   collision-engine lib (#4). The registry only offers the engine that its DSO
   loads; a missing DSO should log + fall back, not crash.
3. Envelope `{max_special_q=2^32−1, max_p=2^27−1}`, `max_threads=4` — already the
   registry values; the driver clamps `q` / skips over-cap `a_d` before dispatch.
4. Found-array overflow is a clamp (line 655), not an `exit` — good. Confirm the
   collision-engine DSO path also skips gracefully on its own overflow rather
   than `exit(-1)`.
5. `thread_num` is pool-local `0..num_threads-1` (confirmed for the msieve pool),
   so `gd->threads[threadid]` is in range.
6. Smoke test both: `stage1_engine=gpu_cubsort` and `stage1_engine=gpu_gerbicz`
   over one forced cell; hits should flow to the shared stage 2, `q` should move,
   the `.p` should fill — same status line as the CPU engines, because it's the
   same shared path from `handle_collision` on down.
7. Then the real comparison: same N / a_d range / `stage2_threads`, swap
   `stage1_engine=` across `cpu_gerbicz` / `gpu_cubsort` / `gpu_gerbicz`, and read
   `poly_score.py` on the one `.p`. This is where GPU stage 1 should finally let
   you push `stage2_threads` toward full core count (producer moved off-CPU).
