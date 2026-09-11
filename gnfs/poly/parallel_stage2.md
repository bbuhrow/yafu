# Single-instance override + parallel stage 2

Two changes that land together (the override alone would run ~THREADS× slower —
it collapses your across-instance stage-2 workers into one; the parallel stage 2
restores and then exceeds that).

The concurrency shape is validated standalone in `s2_pattern_test.c` (4 producers
/ 3 consumers / 1.2M tasks: zero cross-worker corruption, zero interleaved file
lines, exact funnel). What can't be compiled here is the msieve-struct wiring —
see the verify list at the end.

---

## Part 1 — `nfs_poly.c`: one msieve instance owns the threading

In `do_msieve_polyselect` (≈ line 897) the machinery spawns `fobj->THREADS`
YAFU workers, each building a `num_threads`-wide msieve pool. Force **one** YAFU
worker and hand the real count to msieve. Range dispatch is dynamic, so the lone
worker still covers the whole a_d space.

**A — force one worker** (top of the function, before the `thread_data` malloc /
spawn loop):
```c
	int saved_threads = fobj->THREADS;   /* the real thread count */
	fobj->THREADS = 1;                    /* one YAFU poly worker */
```

**B — restore on every exit path** (Path-A abort, Path-B, normal end — after the
workers are torn down, before each `return`):
```c
	fobj->THREADS = saved_threads;
```

**C — give msieve the real count.** `init_poly_threaddata` passes
`(uint32_t)fobj->THREADS` (now 1) to `msieve_obj_new`. Instead:
- add a param `int num_msieve_threads` to `init_poly_threaddata`,
- at the call site pass `saved_threads`,
- use `(uint32_t)num_msieve_threads` in the `num_threads` slot of **both**
  `msieve_obj_new(...)` calls (CUDA + non-CUDA).

Then `obj->num_threads = saved_threads` natively and the driver's stage-1 pool is
sized correctly with no further change.

Caveats: `deadline /= fobj->THREADS` comes out right at 1 (the single instance
runs the full deadline). The Path-B test-sieve estimate (`… / fobj->THREADS`)
would be off, but it's the embryonic path we're relocating. Restore promptly so
nothing downstream sees THREADS=1.

---

## Part 2 — parallel stage 2

Per-worker private optimization state; the file, best-E, stats, and the poly heap
are shared under one lock. `thread_num` is pool-local `0..S-1`
(`common/thread.c:481`), so indexing `stage2_workers[thread_num]` is safe.

### 2a `poly_skew.h`

```c
/* one per stage-2 worker: private sizeopt/rootopt state, shared file/best/stats */
typedef struct {
	poly_sizeopt_t          sizeopt_data;
	poly_rootopt_t          rootopt_data;
	sizeopt_callback_data_t sizeopt_callback_data;
	rootopt_callback_data_t rootopt_callback_data;
} stage2_worker_t;
```
`rootopt_callback_data_t` (already gained `stats`, `all_poly_file`) gains:
```c
	mutex_t *file_lock;   /* shared; serializes all_poly_file + save_poly */
```
`poly_stage1_t` gains:
```c
	stage2_worker_t *stage2_workers;   /* S bundles */
	uint32 num_stage2_workers;
```

### 2b `poly_skew.c` — build S private bundles

Extract the (currently single) sizeopt/rootopt setup into a loop. Copy the exact
field-fill from the original, including the per-degree `min_e_bernstein` block:
```c
static stage2_worker_t *
build_stage2_workers(uint32 s, msieve_obj *obj, mpz_t n, uint32 degree,
		<params_type> *params, poly_config_t *config,
		FILE *shared_poly_file, mutex_t *file_lock,
		poly_stage_stats_t *stats)
{
	uint32 i;
	stage2_worker_t *w = (stage2_worker_t *)xcalloc(s, sizeof(stage2_worker_t));

	for (i = 0; i < s; i++) {
		stage2_worker_t *k = w + i;

		/* size opt — PRIVATE */
		poly_sizeopt_init(&k->sizeopt_data, sizeopt_callback,
				&k->sizeopt_callback_data);
		mpz_set(k->sizeopt_data.gmp_N, n);
		k->sizeopt_data.degree           = degree;
		k->sizeopt_data.max_stage1_norm  = params->stage1_norm;
		k->sizeopt_data.max_sizeopt_norm = params->stage2_norm;
		k->sizeopt_data.best_saved_combined_e = 0.0;
		k->sizeopt_data.num_rootopt = 0;
		k->sizeopt_data.num_saved   = 0;

		/* root opt — PRIVATE */
		poly_rootopt_init(&k->rootopt_data, obj, rootopt_callback,
				&k->rootopt_callback_data);
		mpz_set(k->rootopt_data.gmp_N, n);
		k->rootopt_data.degree           = degree;
		k->rootopt_data.max_sizeopt_norm = params->stage2_norm;
		k->rootopt_data.min_e            = params->final_norm;
		k->rootopt_data.min_e_bernstein  = 0;
		/* ---- copy the original degree==4 / degree==5 min_e_bernstein block here ---- */

		/* link + SHARED */
		k->sizeopt_callback_data.rootopt          = &k->rootopt_data;
		k->sizeopt_callback_data.rootopt_callback = &k->rootopt_callback_data;
		k->sizeopt_callback_data.stats            = stats;
		k->rootopt_callback_data.config           = config;
		k->rootopt_callback_data.all_poly_file    = shared_poly_file; /* SHARED */
		k->rootopt_callback_data.file_lock        = file_lock;        /* SHARED */
		k->rootopt_callback_data.stats            = stats;
	}
	return w;
}
```
Orchestrator (full-pipeline branch) — open the file once, one lock, parse S,
build, point stage 1 at worker[0] as the default `callback_data`:
```c
	uint32 num_s2 = 1;
	mutex_t s2_file_lock;
	stage2_worker_t *workers;
	FILE *shared_poly_file;

	if (obj->nfs_args) {
		const char *tmp = strstr(obj->nfs_args, "stage2_threads=");
		if (tmp) num_s2 = MAX(1, atoi(tmp + 15));
	}
	sprintf(buf, "%s.p", obj->savefile.name);
	shared_poly_file = fopen(buf, "a");        /* the ONE all-poly file */
	mutex_init(&s2_file_lock);

	workers = build_stage2_workers(num_s2, obj, n, degree, params, config,
			shared_poly_file, &s2_file_lock, &poly_stats);

	poly_stage1_init(&stage1_data, stage1_callback, &workers[0].sizeopt_data);
	stage1_data.stage2_workers     = workers;
	stage1_data.num_stage2_workers = num_s2;
```
This replaces the single-instance POLYSIZE/POLYROOT setup + the per-thread `.p`
open. After `poly_stage1_run` and the stats report/free: free each worker
(`poly_sizeopt_free` / `poly_rootopt_free`), `fclose(shared_poly_file)`,
`mutex_free(&s2_file_lock)`, `free(workers)`.

### 2c `poly_skew.c` — serialize the two shared writes

`rootopt_callback` writes the shared file **and** calls `save_poly(config, …)`
(shared heap). Put both under the lock:
```c
	if (data->file_lock) mutex_lock(data->file_lock);

	fprintf(data->all_poly_file, "# norm %le alpha %lf e %.3le rroots %u\n"
			"skew: %.2lf\n", size_score, root_score, combined_score,
			num_real_roots, skewness);
	for (i = 0; i <= degree; i++)
		gmp_fprintf(data->all_poly_file, "c%u: %Zd\n", i, coeff1[i]);
	for (i = 0; i <= 1; i++)
		gmp_fprintf(data->all_poly_file, "Y%u: %Zd\n", i, coeff2[i]);
	fflush(data->all_poly_file);

	save_poly(config, &poly);              /* shared heap -> inside the lock */

	if (data->file_lock) mutex_unlock(data->file_lock);

	/* stats hook (from the wiring patch) stays here, OUTSIDE the file lock */
	if (data->stats) {
		char adbuf[64];
		gmp_snprintf(adbuf, sizeof(adbuf), "%Zd", coeff1[degree]);
		poly_stats_add_rootopt(data->stats, combined_score, adbuf);
	}
```
`sizeopt_callback` needs no lock — it only touches its worker's private
`sizeopt_data`/`rootopt_data` (plus the stats bump, which is internally locked).

### 2d `stage1.c` — route each hit to its worker's bundle

`stage1_sieve_data_t` (beside `engine`, `stats`):
```c
	stage2_worker_t *stage2_workers;
	uint32 num_stage2_workers;
```
`stage1_hit_data_t` gains a back-pointer so the worker is reachable by `thread_num`:
```c
	stage1_sieve_data_t *d;
```
`handle_collision` — set it (next to `hit_data->stats = d->stats;`):
```c
	hit_data->d = d;
```
`poly_stage1_run` — copy through (next to `sieve_data.stats = data->stats;`):
```c
	sieve_data.stage2_workers     = data->stage2_workers;
	sieve_data.num_stage2_workers = data->num_stage2_workers;
```
Size the pool by S:
```c
	d->stage2_threadpool = threadpool_init(
			MAX(1, d->num_stage2_workers), 1000, &thread_control);
```
Route in `stage1_hit_run` (`thread_num` is pool-local 0..S-1):
```c
static void
stage1_hit_run(void *data, int thread_num)
{
	stage1_hit_data_t *hit_data = (stage1_hit_data_t *)data;
	void *extra = hit_data->callback_data;          /* single-instance fallback */

	if (hit_data->d && hit_data->d->num_stage2_workers > 0)
		extra = &hit_data->d->stage2_workers[thread_num].sizeopt_data;

	hit_data->callback(hit_data->ad, hit_data->p, hit_data->m, extra);

	if (hit_data->stats)
		poly_stats_stage2_done(hit_data->stats);
}
```

### Config

`stage2_threads=S` rides the `nfs_stage1_args` passthrough (default 1 = today).
Tune S up watching `q`: when it stops pegging at 1000, stage 2 keeps up; total
worker threads ≈ (stage-1 `num_threads`) + S, so keep the sum near your core
count.

---

## Verify on build (highest-risk change — the pattern is proven, the wiring isn't)

1. `thread_num` pool-local 0..S-1 — confirmed (`common/thread.c:481`); the
   `stage2_workers[thread_num]` index is safe.
2. `poly_sizeopt_init`/`poly_rootopt_init` have no static globals (confirmed by
   scan) — S copies are independent.
3. **`best_saved_combined_e` semantics.** If root-opt uses it as a *save gate*
   (only save when better than best-so-far), per-worker copies mean each worker
   saves its own local-bests → more polys emitted. If it's only a *stat*, no
   behavior change. Check `stage2.c` / `optimize.c`; if it's a gate and you want
   the old global behavior, move it into the shared struct under the lock. The
   job-wide best-E on the status line is already correct (from `poly_stats`).
4. **`save_poly` / `poly_config_t` heap** is shared — kept inside the file lock
   above. Confirm `save_poly` touches nothing else cross-worker.
5. Runtime: funnel `hits ≥ sizeopt ≥ rootopt`; `q` drains to 0 at end; no crash
   at S>1; CPU ≈ (num_threads + S) with the Part-1 override; the single `.p`
   file's lines stay intact (no interleave — the whole point of the lock).
6. Teardown frees every worker's sizeopt/rootopt, closes the one file, frees the
   lock and the array.
