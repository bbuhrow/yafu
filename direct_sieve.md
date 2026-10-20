# Direct Interval Sieve Design (TRY_DIRECT_SIEVE) — design notes, not yet implemented

*Written 2026-10-20 to preserve a full design pass made in session. Companion to
`L2_SIEVE_DEBUG_NOTES.md` (read that first: sections 1, 3, 4, 5 contain the
component map, the bug lessons, and tooling). Companion to `direct_sieve.md`'s
sibling experiment TRY_LIST_SIEVE (section 6 of the debug notes).*

## 1. Goal

New experiment guarded by `#define TRY_DIRECT_SIEVE` (add to `include/qs_impl.h`
near the other TRY_ defines, line ~79-94). Uses the whole-interval sieve idea of
TRY_L2_SIEVE but generalizes it:

- **any number of blocks** (8 gets replaced with `sconf->num_blocks`; the
  interval is `num_blocks * 32768` = `sconf->sieve_interval`).
- **both sides' intervals stay in memory simultaneously** (pos and neg), so the
  large primes `[large_B, B)` can be sieved for both signs in one pass and the
  trial-division scan becomes a separate, deferred pass. This requires 2x the
  `dconf->sieve` storage in `SIQS.c->siqs_dynamic_init()` (~line 3425).
- **bucket machinery for [large_B, B) is retired** (no scatter, no
  lp_sieveblock_ptr, no tdiv_LP). Instead: one direct sieve pass + one fused
  resieve/poly-update pass over `[large_B, B)`.
- **[med_B, large_B) keeps the TRY_SPARSE_SIEVE machinery**
  (`large_sieveblock_32k_avx512bw`, `tdiv_large_32k_avx2`), just with the
  boundary moved from `fb_blk8i` to `large_B`.
- **trial division is deferred**: the block loops sieve only; the scan
  (`check_relations_siqs_16_avx512`) runs block-by-block afterwards, and
  relation finalization runs once after the large-prime resieve.

## 2. Key facts established this session (don't re-derive)

1. `sconf->sieve_interval = qs_blocksize * num_blocks` (SIQS.c:4345).
   The fb boundary computation at SIQS.c:4499-4515 already finds the first
   16-aligned index with `prime > sieve_interval` and stores it in
   `fb->large_B` — **but only in the non-L2 branch**. Under
   `#if defined(TRY_L2_SIEVE)` line 4513 forces `large_B = fb_blk8i`.
   Under TRY_DIRECT_SIEVE we want the computed value (for NB=8 it equals
   fb_blk8i exactly, since 8*32768 = 262144 = the L2 boundary; for other NB it
   is the correct "prime > interval" boundary). Guard precedence needed:
   `#if defined(TRY_DIRECT_SIEVE)` keep computed `i`, `#elif defined(TRY_L2_SIEVE)`
   force fb_blk8i. Also enforce the existing clamp at 4598-4599
   (`large_B < med_B -> med_B`).

2. `firstRoots_32k` (poly_roots_32k.c:1939) already computes, for the first
   poly of each 'a':
   - L2 section (2279-2345): firstroots1/2 (min/max sorted), sieve_p/n->root1/2,
     and rootupdates for [med_B, bound) where bound is fb_blk8i (or B under
     TRY_LIST_SIEVE).
   - Middle bucket loop (2350): scatter + firstroots (UNSORTED) + rootupdates
     for [fb_blk8i, large_B) — EMPTY under current TRY_L2 config, and must be
     SKIPPED under TRY_DIRECT_SIEVE (no buckets; range covered by L2 section).
   - Third bucket loop (2442): scatter + firstroots (UNSORTED) + rootupdates
     for [large_B, B) — under TRY_DIRECT_SIEVE keep firstroots + rootupdates,
     skip the FILL_ONE_PRIME_P/N scatters. The unsortedness does not matter:
     the direct sieve treats {root1, root2} as an unordered pair (each root is
     checked against the interval independently, both sides).

3. `nextRoots_32k_avx2_small` (update_poly_roots_32k_avx2.c:3147) is called
   AFTER `nextB` at SIQS.c:2700. `nextB` (new_poly.c:1058) only updates
   mpz_poly_b/c; the transition parameters it applies are
   `curr_poly->nu[numB]` / `curr_poly->gray[numB]` and the precomputed
   `rootupdates[(v-1)*bound + j]` (filled once per 'a' by firstRoots for ALL
   j < s). Therefore a combined resieve+update pass can run BEFORE `nextB`
   inside process_poly and compute the identical next-poly roots that
   nextRoots would compute after it: read `v = curr_poly->nu[dconf->numB]`,
   `sign = curr_poly->gray[dconf->numB]`, `ptr = &rootupdates[(v-1)*bound + ...]`.
   Under TRY_DIRECT_SIEVE, nextRoots' L2 sections get
   `l2upd_bound = fb->large_B` (they keep updating firstroots + sieve_p/n roots
   for the sparse-sieve range [med_B, large_B) only).

4. **Root state per procedure (critical for deferral):**
   - `med_sieve_ptr` advances `comp_sieve_p/n->root1/2` (uint16, all med_B
     entries; small primes below sieve_small_fb_start are static per poly) by
     one block per call.
   - `large_sieveblock_32k_avx512bw` (med_sieve_32k_avx2.c:1128) advances and
     wraps `dconf->sieve_p/n->root1/2` over [med_B, bound) per block.
   - `tdiv_med_ptr` / `resieve_med_ptr` / `tdiv_large_32k_avx2` only READ these
     root arrays (no mutation), but they expect the state as of *that block*:
     resieve_med's "corrections" trick checks roots already advanced past the
     block; tdiv_large checks roots wrapped to the next interval.
   - => if the scan is deferred out of the block loops, the per-block root
     state must be **snapshot after sieving block i, restored before scanning
     block i**. Nothing else mutates these arrays between, so snapshot/restore
     of the 4 arrays suffices:
     `comp_sieve_p->root1/root2` (med_B x uint16, full arrays — small-prime
     entries are identical across blocks so including them is harmless),
     `comp_sieve_n->root1/root2`, and `dconf->sieve_p->root1/root2` +
     `sieve_n->root1/root2` over [med_B, large_B) (uint32).
   - `set_aprime_roots(..., 1)` invalidates a-prime roots before med_sieve;
     `set_aprime_roots(..., 0)` restores them before the scan. Order in the
     scan loop: restore snapshot FIRST, then set_aprime_roots(0) (snapshot
     holds invalid markers for a-prime slots), then scan_ptr.

5. `check_relations_siqs_16_avx512` (tdiv_scan.c:947) and `filter_SPV`
   (tdiv_small.c:75) are the only scan-side consumers of sieve BYTES.
   filter_SPV reads `sieve[reports[rn]]` from the pointer passed in — so under
   TRY_DIRECT_SIEVE it must be passed `sieve + blocknum*32768`. The AVX512BW
   branch of the scan (use_dlp >= 2, tdiv_scan.c:964-993) also needs the
   blocknum offset (currently latent bug, debug notes section 4.3).
   tdiv_med/resieve_med/tdiv_large/tdiv_LP never touch dconf->sieve bytes.

6. Relation finalization lives in `trial_divide_Q_siqs` (tdiv.c:212): it uses
   `dconf->curr_poly` (a-poly factors, `curr_poly->index`) and calls
   `buffer_relation`. Therefore **all of pass 2 (resieve+update), pass 3
   (finalization) must run before `nextB`** inside process_poly. `nextB`
   must not run until every report of the interval has been finalized.

7. `dconf->Qvals`, `valid_Qs`, `smooth_num`, `reports` are already POINTERS
   (qs_impl.h:719-722); only `fb_offsets` is a fixed member
   `uint32_t fb_offsets[MAX_SIEVE_REPORTS][MAX_SMOOTH_PRIMES]`. Under the
   guard, change it to `uint32_t (*fb_offsets)[MAX_SMOOTH_PRIMES];` (all
   existing `dconf->fb_offsets[r][k]` indexing stays valid) and allocate
   `2 * num_blocks * MAX_SIEVE_REPORTS` rows. All five arrays get
   `2 * num_blocks * MAX_SIEVE_REPORTS` entries; global report id
   `gid = side * (num_blocks * MAX_SIEVE_REPORTS) + block * MAX_SIEVE_REPORTS
   + local_report`. Per-block code keeps using LOCAL indices by repointing the
   five dconf pointers per block (see 4.4 below).

8. Allocation sites (SIQS.c): reports :3785, Qvals :3787-3792 (mpz_init each),
   valid_Qs :3794, smooth_num :3795, fb_offsets is a struct member (no
   allocation today). Teardown: frees at :6283-6291 (loop bound
   MAX_SIEVE_REPORTS must become the new count). Sieve allocation :3426-3433;
   sieve_p/n fb arrays :3318-3340 (sievefb_alloc); the prime/logp copy loops
   at :3450-3490 (extend the med_B..fb_blk8i copy to large_B under the guard).

9. Guard precedence design: everywhere the new behavior conflicts with
   TRY_LIST_SIEVE, write `#if defined(TRY_DIRECT_SIEVE) ... #elif
   defined(TRY_LIST_SIEVE) ... #else/#endif`. Then both defines may be active
   at once (DIRECT wins; the list sieve's init is never called so its B-sized
   ctx never touches the large_B-sized sieve_p/n arrays). No need to comment
   out TRY_LIST_SIEVE when testing, though it may be commented for cleanliness.

## 3. Per-poly pipeline (process_poly, SIQS.c ~2180-2740)

```
for (numB = 1; numB < maxB; numB++, tot_poly++)     // existing loop
{
    // ---- stage A: sieve positive side (deferred scan) ----
    memset(dconf->sieve, blockinit, num_blocks * 32768);            // 2248 area
    memset(dconf->sieve + minblock*32768, blockinit - 4, 32768);    // keep minblock semantics
    for (i = 0; i < num_blocks; i++)                                // 2253
    {
        blk_sieve = dconf->sieve + i * 32768;                       // 2255 condition: add TRY_DIRECT_SIEVE
        set_aprime_roots(..., fb_sieve_p, 1);
        // per-block memsets at 2283/2303 DISABLED (guard them off)
        med_sieve_ptr(blk_sieve, ..., i == minblock ? blockinit-4 : blockinit);
        large_sieveblock_32k_avx512bw(blk_sieve, ...);              // 2295/2315, bound -> large_B
        // NO lp_sieveblock_ptr, NO scan_ptr
        // NEW: snapshot comp_sieve_p root1/2 + sieve_p root1/2 [med_B, large_B) -> slot [side=0][i]
    }

    // ---- stage B: sieve negative side (same, mirrored) ----
    //   interval init at 2395 area: memset(dconf->sieve + num_blocks*32768, ...)
    //   blk_sieve = dconf->sieve + (num_blocks << 15) + i*32768
    //   snapshot -> slot [side=1][i]

    // ---- stage C: pass 1 — direct interval sieve for [large_B, B) ----
    direct_sieve_largeprimes_32k(sconf, dconf);      // both sides, one prime loop

    // ---- stage D: deferred block-by-block scan, pos then neg ----
    for (i = 0; i < num_blocks; i++)
    {
        restore snapshot [0][i]; set_aprime_roots(..., fb_sieve_p, 0);
        repoint Qvals/fb_offsets/smooth_num/valid_Qs/reports to gid base (side 0, block i)
        scan_ptr(i, 0, sconf, dconf);                // check_relations: scan+filter+med+resieve_med+tdiv_large only
        direct_block_reports[i] = dconf->num_reports;
    }
    dconf->sieve += (num_blocks << 15);              // point scan at neg interval
    for (i = 0; i < num_blocks; i++) { ... side 1 ... }
    dconf->sieve -= (num_blocks << 15);              // restore

    // ---- stage E: pass 2 — fused resieve + poly-update for [large_B, B) ----
    direct_resieve_update_largeprimes_32k(sconf, dconf);   // BEFORE nextB

    // ---- stage F: pass 3 — deferred relation finalization ----
    for (side = 0; side < 2; side++)
        for (b = 0; b < num_blocks; b++)
            for (j = 0; j < direct_block_reports[side*num_blocks + b]; j++)
            {
                gid = (side*num_blocks + b) * MAX_SIEVE_REPORTS + j;
                if (dconf->valid_Qs_base[gid])
                    trial_divide_Q_siqs(gid, side, dconf->numB - 1, b, sconf, dconf);
            }
    // (keep total_surviving_reports counting in check_relations as today)

    nextB(dconf, sconf, 1);                          // 2621, unchanged
    nextRoots_32k_avx2_small(sconf, dconf);          // 2700, L2 bound -> large_B
    // knl medbucket/largebucket calls (2705-2720) SKIPPED under the guard
    //   (same as the existing TRY_LIST_SIEVE skip at 2702-2704)
}
```

minblock is computed at SIQS.c:2183-2195 before the p-side loop — available
for the interval-init memsets.

## 4. New/changed code, file by file

### 4.1 include/qs_impl.h
- Add `#define TRY_DIRECT_SIEVE` near line 94 with a comment block describing
  the scheme and its requirements.
- Under `#if defined(TRY_DIRECT_SIEVE)`:
  - `uint32_t (*fb_offsets)[MAX_SMOOTH_PRIMES];` instead of the 2D array
    (qs_impl.h:721).
  - New dynamic_conf_t fields: base pointers `Qvals_base` (mpz_t*),
    `valid_Qs_base`, `smooth_num_base`, `reports_base`,
    `fb_offsets_base` (row pointer), `direct_flags_p` / `direct_flags_n`
    (uint32_t*, interval-sized), snapshot pointers `direct_comp_snap`
    (uint16_t*) and `direct_fbroots_snap` (uint32_t*), and
    `uint32_t direct_block_reports[128];` (num_blocks is <= 22 in the param
    table; note the limit in a comment).
  - Prototypes for the two new functions.

### 4.2 SIQS.c — siqs_dynamic_init (~3290-3830) and teardown (~6280-6300)
- `sievefb_alloc`: `#if defined(TRY_DIRECT_SIEVE)` = `fb->large_B`
  (#elif TRY_LIST_SIEVE = B, #else fb_blk8i).
- Extend the prime/logp copy loop (3450-3490, currently [med_B, fb_blk8i)) to
  large_B under the guard.
- Sieve allocation (3426): under the guard
  `xmalloc_align(qs_blocksize * (2 * num_blocks + 4))`, then
  `dconf->sieve = &base[2 * qs_blocksize]` — 2 guard blocks below, 2 above;
  pos interval = [dconf->sieve, +NB*32768), neg = next NB*32768.
- Allocate the five report arrays with
  `nreports = 2 * num_blocks * MAX_SIEVE_REPORTS` entries (reports, Qvals +
  mpz_init loop, valid_Qs, smooth_num, fb_offsets rows of MAX_SMOOTH_PRIMES);
  store the base pointers. NOTE: fb_offsets alone is
  2*8*2048*200*4 = ~26 MB per thread at NB=8 — acceptable for the experiment,
  but flag it; a shared arena could shrink it later.
- Allocate flags: `2 * (num_blocks << 15) * sizeof(uint32_t)` and MEMSET TO
  ZERO (the TRY_LIST_SIEVE segfault lesson: never leave flags uninitialized —
  debug notes section 6 "Segfault during integration").
- Allocate snapshots: comp 2*NB*2*med_B uint16; fb-roots 2*NB*2*(large_B -
  med_B) uint32.
- Teardown: free the new allocations; mpz_clear loop bound = nreports; free
  fb_offsets (now heap).

### 4.3 process_poly (SIQS.c 2243-2492) — per section 3 above
Guard details:
- Interval-init points (2246 / 2394): `#if defined(TRY_DIRECT_SIEVE)` branch
  FIRST (before the TRY_LIST_SIEVE / L2 branches): whole-interval memset +
  minblock blockinit-4 memset (minblock memset is NOT applied per block
  anymore).
- blk_sieve selection (2255/2402): `#if defined(TRY_DIRECT_SIEVE) ||
  (defined(TRY_L2_SIEVE) && !defined(TRY_SPARSE_SIEVE))` -> offset form; neg
  side adds `(num_blocks << 15)`.
- Per-block memsets (2282-2284, 2302-2304, 2427-2429, 2447-2449): add
  `&& !defined(TRY_DIRECT_SIEVE)` to the keep-condition.
- lp_sieveblock_ptr (2352/2474): skip under DIRECT (like the LIST branch).
- Remove scan_ptr calls from the block loops (2364, 2486) under the guard;
  set_aprime_roots(...,1) stays in the sieve loop.
- Insert stages C-F after the neg-side loop (after ~2492, before the timing
  code at 2495 and nextB at 2621).
- Snapshot insert points: end of each block's sieve work in both loops
  (after med+sparse sieve, before loop end).
- USE_SS_SEARCH blocks repoint dconf->sieve — the experiment config has SS
  off; do not attempt to make DIRECT+SS coexist (same stance as TRY_LIST_SIEVE).

### 4.4 tdiv_scan.c — check_relations_siqs_16_avx512 (947-1114)
- Scan base: under DIRECT use `dconf->sieve + blocknum * 32768` in BOTH the
  use_dlp>=2 BW branch (971-972) and the normal branch (1013-1019 condition).
- filter_SPV call (1062-1070): pass `dconf->sieve + blocknum * 32768` under
  DIRECT (this is the "block offsets" requirement).
- Keep tdiv_med_ptr / resieve_med_ptr; call tdiv_large_32k_avx2 under DIRECT
  (extend the SPARSE branch condition).
- Per-report loop (1088-1110): under DIRECT keep the valid_Qs /
  total_surviving_reports accounting but SKIP tdiv_LP_ptr and
  trial_divide_Q_siqs (deferred). (tdiv_LP is bucket-only; no buckets exist.)

### 4.5 med_sieve_32k_avx2.c
- `large_sieveblock_32k_avx512bw` line 1173: `bound = full_fb->large_B` under
  the guard (fb_blk8i otherwise). Scalar tail comment mentions 8*32768 — the
  L2-semantics tail is unchanged and correct for any interval (roots wrapped
  by `interval` would be wrong here — NOTE the tail subtracts 32768 per block
  semantics; the sparse sieve wraps per BLOCK (32768), which is correct for
  primes > 32768 covering [med_B, large_B); do not confuse with the
  whole-interval L2 tail).
- NEW `direct_sieve_largeprimes_32k(sconf, dconf)`:
  - `interval = sconf->num_blocks << 15; large_B = fb->large_B;`
    `sievep = dconf->sieve; sieven = dconf->sieve + interval;`
  - 16-wide loop j in [large_B, B) (B and large_B are both 16-aligned):
    load prime, fr1 = firstroots1[j], fr2 = firstroots2[j]; logp approximated
    as `full_fb->list->logprime[j+15]` (same approximation style as
    large_sieveblock); masks `fr < interval` for both roots;
    neg roots `prime - fr1`, `prime - fr2` with their own masks
    (nr in (0, p]; nr == p when fr == 0 is auto-masked since p > interval).
  - Scalar hit loops over set bits (large_sieveblock style):
    `sievep[hit] -= logp;` / `sieven[hit] -= logp;`. Do NOT wrap or store
    roots here (pass 2 re-derives everything from firstroots).
  - First poly of an 'a': firstroots come from firstRoots' third loop;
    later polys: from pass 2 of the previous poly. No other init needed.
- Optionally put pass 2 here too (or in tdiv_med_32k_avx2.c next to
  resieve_large_32k_avx2).

### 4.6 NEW `direct_resieve_update_largeprimes_32k(sconf, dconf)`
- Reads: `v = dconf->curr_poly->nu[dconf->numB]`,
  `sign = dconf->curr_poly->gray[dconf->numB]`,
  `ptr = &rootupdates[(v - 1) * bound + large_B]` (bound = fb->B),
  `tag = dconf->tot_poly & 0x7FFFF` (per-poly tag for flags; tot_poly is
  per-thread and unique per process_poly body — incremented at loop increment,
  i.e. after the body runs).
- 16-wide loop j in [large_B, B):
  - Resieve: fr1/fr2 (current poly, NOT yet updated — "firstroots values are
    not updated yet, so we load them"), pos hits: for each root < interval,
    `f = direct_flags_p[hit]`; if `(f >> 12) == tag` and `(f & 0xFFF)`:
    `rpt = (f & 0xFFF) - 1; gid = (hit >> 15) * MAX_SIEVE_REPORTS + rpt;`
    GUARDED division into `Qvals_base[gid]`: only if
    `mpz_tdiv_ui(Q, prime) == 0` then
    `do { fb_offsets_base[gid][++smooth_num_base[gid]] = j;
    mpz_tdiv_q_ui(Q, Q, prime); } while (mpz_tdiv_ui(Q, prime) == 0);`
    (guarded first division — the Bug 2 / Bug 4 failure class must be
    impossible: a stale/aliased flag can then only waste a tdiv_ui, never
    corrupt Q; the list sieve used the same guard).
    Neg side identical with `direct_flags_n`, roots `prime - fr`,
    `gid = (num_blocks + (hit >> 15)) * MAX_SIEVE_REPORTS + rpt`.
  - Update (mirror nextRoots' L2 vector loop): sign > 0:
    `r = fr - u; if ((int)r < 0) r += prime;` (COMPUTE_NEXT_ROOTS_P
    semantics / _mm512_mask_add on the > mask); sign < 0:
    `r = fr + u; if (r >= prime) r -= prime;`. Store back into
    firstroots1/2[j]. No sorting needed (unordered pair semantics).
  - ptr advances 16 per iteration. Scalar tail for the last < 16 primes.
- flags encoding written by the scan: during check_relations' report loop
  (or in check_relations after num_reports is final), for each report of
  block b: `direct_flags_p[parity-independent pos]`:
  `flags[side][b*32768 + reports[j]] = (tag << 12) | (j + 1);`
  (needs the tag available in tdiv_scan.c — pass it via a dynamic_conf_t
  field set in process_poly stage C, e.g. `dconf->direct_poly_tag`, or
  recompute `tot_poly & 0x7FFFF` there). rpt+1 <= 2048 fits 12 bits because
  num_reports is clamped to MAX_SIEVE_REPORTS-1 = 2047.

### 4.7 poly_roots_32k.c — firstRoots_32k (1939+)
- L2 section bound (2285-2288): `#if defined(TRY_DIRECT_SIEVE)` ->
  `[med_B, large_B)`; #elif TRY_LIST_SIEVE -> B; #else fb_blk8i.
- Middle bucket loop (2350): skip entirely under DIRECT
  (`#if !defined(TRY_DIRECT_SIEVE) ... #endif` around it, keeping the
  `bound_val/check_bound/logp` setup consistent).
- Third bucket loop (2442): under DIRECT keep firstroots1/2 + rootupdates +
  CHECK_NEW_SLICE, but compile out FILL_ONE_PRIME_P / FILL_ONE_PRIME_N.

### 4.8 tdiv_med_32k_avx2.c — tdiv_large_32k_avx2 (916)
- `bound` at :950 and the scalar tail loop at :975: use `fb->large_B` instead
  of `fb->fb_blk8i` under the guard.

## 5. Invariants / correctness reasoning (as with section 6 of the debug notes)

- Every sieve byte of both intervals is final before any scan: med (per block
  during stage A/B), sparse [med_B, large_B) (per block), direct large primes
  (stage C, whole interval both sides). The scan's hi-bit report detection
  therefore sees complete bytes — the same property the L2 scheme relies on.
- Every (prime, root) hit of [large_B, B) is visited exactly once by pass 1
  and at most once by pass 2 (p > interval => at most one hit per root);
  both passes read the SAME root values (firstroots are only advanced in
  pass 2's update half), so the resieve covers exactly the sieve's hit set.
- All divisions into report Qvals are guarded (divisibility checked first);
  a wrong root/flag can lose a factor but cannot corrupt Q (the Q->0 /
  smooth_num-explosion cascade from Bugs 2/4 cannot occur).
- Relation finalization (trial_divide_Q_siqs -> buffer_relation) runs before
  nextB, so curr_poly context is intact.
- Snapshot/restore keeps comp_sieve and sparse root state exactly as the
  interleaved (stock) pipeline would have had it at each scan_ptr call.
- num_blocks generality: allocation is 2*NB + 4 blocks (no more hardcoded 8 /
  20-block assumptions on this path; debug notes 4.4 does not apply here).

## 6. Integration steps remaining

1. Apply 4.1-4.8 (single cohesive edit set; guards written with DIRECT
   precedence per 2.9).
2. Build (Makefile + config.mk, SKYLAKEX / USE_AVX512F + USE_AVX512BW config
   from the debug notes; no new files, so no Makefile changes needed — the
   two new functions go into existing files).
3. Test: 256-bit input in N.txt (factors in factors.txt) and the 261-bit
   79-digit input from debug notes section 3 Bug 4, `-siqsNB 8 -threads 1`,
   expecting full factorization, zero "does not divide" errors, no crash.
4. Test num_blocks != 8 (e.g. `-siqsNB 12` and `-siqsNB 16`) — the whole
   point of this experiment; large_B then differs from fb_blk8i.
5. Compare relation counts / timing against the TRY_LIST_SIEVE and
   stock-TRY_SPARSE builds on the same inputs.
6. If relations double-record a-primes or filter/LA misbehaves, check the
   a-prime caveat from TRY_LIST_SIEVE (debug notes section 6 end) first.

## 7. Assumptions made in this design (confirm before/while implementing)

1. **Deferred finalization**: the per-block report loop in check_relations
   must NOT call trial_divide_Q_siqs under the guard; finalization moves to
   stage F. This is forced by "large primes are divided after all
   block-by-block steps" — otherwise relations needing large primes would be
   misclassified at block time.
2. **Snapshot/restore** of comp_sieve and sparse roots is required (see 2.4);
   the strategy as stated implies it but does not spell it out.
3. **Memory**: report storage grows to 2*NB*MAX_SIEVE_REPORTS rows
   (~26 MB/thread fb_offsets at NB=8) plus ~2 MB flags and ~1.7 MB snapshots.
   Fine for an experiment; can be arena-compressed later.
4. TRY_DIRECT_SIEVE is incompatible (untested) with USE_SS_SEARCH,
   USE_BATCHPOLY / USE_XLBATCHPOLY, and TLP (use_dlp >= 2) paths — same
   stance as the previous experiments; the BW scan branch gets the blocknum
   offset anyway per 4.4.
