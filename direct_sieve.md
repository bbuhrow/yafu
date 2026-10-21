# Direct Interval Sieve Design (TRY_DIRECT_SIEVE)

*First design pass 2026-10-20; revised the same day after review and a
first (unfused) implementation. This file describes the design as built,
plus the fused pass (section 6) that is being implemented next. Companion
to `L2_SIEVE_DEBUG_NOTES.md` (component map, bug lessons, tooling).*

## 1. Goal

`#define TRY_DIRECT_SIEVE` (include/qs_impl.h, after TRY_LIST_SIEVE) sieves
the primes larger than the sieve interval, fb indices **[large_B, B)**,
directly into whole intervals, with no bucket sieve:

- **any number of blocks**: interval = `num_blocks * 32768` =
  `sconf->sieve_interval`.
- **both sides' intervals in memory at once**: pos at
  `dconf->sieve[0, interval)`, neg at `[interval, 2*interval)`. The
  allocation is `2*num_blocks + 4` blocks with 2 guard blocks on each end.
- **[med_B, large_B) keeps the TRY_SPARSE_SIEVE machinery**
  (`large_sieveblock_32k_avx512bw` per block, `tdiv_large_32k_avx2` in
  the scan), unchanged.
- **only relation finalization is deferred**. Sieving and the scan stay
  block by block.

Header guards: TRY_DIRECT_SIEVE `#undef`s TRY_LIST_SIEVE (DIRECT wins, so
the list code compiles out), and `#error`s unless TRY_L2_SIEVE +
TRY_SPARSE_SIEVE + USE_AVX512F are defined, or if USE_SS_SEARCH /
USE_BATCHPOLY / USE_XLBATCHPOLY are defined. At runtime it exits if
scan_ptr is not `check_relations_siqs_16_avx512` (the only scan that defers
finalization), or if an a-prime index is >= large_B (the direct sieve never
consults set_aprime_roots).

## 2. Key facts (verified against the code; don't re-derive)

1. **Boundary.** `fb_blk8i` is computed against `sieve_interval`
   (SIQS.c fb bounds), not 8*32768, so under TRY_L2_SIEVE
   `large_B == fb_blk8i` for every NB. The sparse sieve, tdiv_large and
   nextRoots' L2 sections all stop at fb_blk8i, so they already end at
   large_B. **Under DIRECT the boundary rounds UP**: it is the first 16-aligned
   index with prime > interval (stock rounds down via `i -= 16`). So every
   prime in [large_B, B) is > interval, which means at most one hit per root
   per side. The <= 15 primes > interval that land below the boundary are
   handled correctly by the sparse sieve.
2. **Roots.** `update_data.firstroots1/2[j]` (int, sized B) hold the
   current poly's roots in [0, p) as an **unordered pair**. The pos hit is
   r, the neg hit is p - r (r == 0 gives p > interval, so it is masked out
   automatically). firstRoots_32k's third loop computes them, plus
   rootupdates, for the first poly of each 'a'. Under DIRECT its
   FILL_ONE_PRIME_P/N bucket scatters are compiled out. The middle bucket
   loop [fb_blk8i, large_B) is empty.
3. **Root updates.** nextB (new_poly.c) and nextRoots_32k_avx2_small both
   read `nu[dconf->numB]` / `gray[dconf->numB]`, and neither changes numB
   (the loop increment does). So code running before nextB can advance the
   roots to the next poly with
   `ptr = &rootupdates[(nu[numB]-1)*B + j]`: for sign > 0,
   `r = r - u (+p if u > r)`; for sign < 0, `r = r + u (-p if >= p)`. The
   knl medbucket/largebucket calls after nextRoots_32k_avx2_small are
   skipped under DIRECT.
4. **The med sieve only subtracts** (`s_init` is unused in
   med_sieveblock_32k_avx512bw), and so does the sparse sieve. Sieve bytes
   can therefore be pre-loaded with large-prime hits before the block loops,
   and the minblock "blockinit - 4" init can be applied as a `-= 4` at any
   time (uint8 arithmetic mod 256 is order-independent). **But it is not
   confined to its block:** it dumps a-prime logs at `blk_sieve[65535]`
   (see 6.3), which the contiguous layout has to protect.
5. **The scan reads sieve bytes only through dconf->sieve**: the sparse
   scan branch, the BW (use_dlp >= 2) branch and filter_SPV. process_poly
   sets `dconf->sieve = blk_sieve` around each scan_ptr call and then
   restores it, so tdiv_scan.c needed no offset changes.
6. **Finalization.** `trial_divide_Q_siqs(report_num, parity, poly_id,
   bnum, ...)` uses `Qvals/fb_offsets/smooth_num/reports[report_num]` plus
   bnum and curr_poly, so it must run before nextB.

## 3. Per-poly pipeline as built (unfused), SIQS.c process_poly

```
for (numB = 1; numB < maxB; numB++, tot_poly++)
{
    memset(sieve, blockinit, 2 * NB * 32768);         // both sides
    if (minblock < NB) memset both sides' minblock to blockinit - 4;
    dpool->num = 0;
    direct_sieve_largeprimes_32k(sieve, ...);          // pass 1, both sides

    for (i = 0; i < NB; i++)                           // pos side
    {
        blk_sieve = sieve + i*32768;                   // neg: sieve + (NB+i)*32768
        set_aprime_roots(..., 1);
        med_sieve_ptr(blk_sieve, ...);                 // no per-block memset
        large_sieveblock_32k_avx512bw(blk_sieve, ...); // sparse [med_B, large_B)
        // no lp_sieveblock_ptr
        set_aprime_roots(..., 0);
        dconf->sieve = blk_sieve; scan_ptr(i, side); dconf->sieve = sieve;
    }
    // neg side identical

    direct_resieve_update_largeprimes_32k(...);        // pass 2
    direct_finalize_survivors(...);
    nextB; nextRoots_32k_avx2_small;                   // no knl bucket calls
}
```

Inside check_relations_siqs_16_avx512, the per-report loop calls
`direct_save_survivor(j, parity, blocknum, ...)` for each valid report
instead of tdiv_LP_ptr + trial_divide_Q_siqs.

## 4. Data structures

### 4.1 Survivor pool (`direct_pool_t`, `dconf->dpool`)
Holds only reports that survive filter_SPV / tdiv_med / resieve_med /
tdiv_large, not all reports.
- `pos[k] = side*interval + block*32768 + block_loc`. Survivors are appended
  in scan order (pos side blocks 0..NB-1, then neg side), so pos[] is
  sorted ascending.
- `Q[k]` gets the report's Qval through `mpz_swap` (O(1); allocations
  circulate between the scan arrays and the pool), along with
  `smooth_num[k]` and `fb_offsets[k*MAX_SMOOTH_PRIMES ...]` (smooth_num+1
  entries copied).
- It starts at 1024 entries and doubles with xrealloc (new mpz's are
  initialized).
- `bitmap`: 2*interval bits (64 KB at NB=8), one bit per survivor position.
  It starts zeroed and is kept zeroed: finalize clears the bits it set. So
  there is no per-poly tag and no stale-flag aliasing.
- The struct members `fb_offsets[MAX_SIEVE_REPORTS][...]` and the scan
  arrays are unchanged.

### 4.2 Resieve hit test (`direct_resieve_hit`)
For a root hit at interval position pos of fb prime index j:
1. Test the bitmap bit (the common case is a miss).
2. On a set bit, binary-search pool->pos.
3. Do a **guarded** division: only if `mpz_tdiv_ui(Q, p) == 0`, then divide
   repeatedly and append j to fb_offsets (capped at MAX_SMOOTH_PRIMES - 1).
   A Q of 0 is skipped.

A wrong hit can only waste a tdiv_ui, never corrupt Q, so the Bug 2/4
failure class cannot occur.

### 4.3 Finalization (`direct_finalize_survivors`)
For each survivor k: clear its bitmap bit, `mpz_swap` Q into `Qvals[0]`, copy
smooth_num / fb_offsets / block_loc into slot 0, set valid_Qs[0] = 1, and
call `trial_divide_Q_siqs(0, side, numB-1, block, ...)`. Then set num = 0.

## 5. Invariants

- Every sieve byte is final before its block is scanned: [large_B, B) is
  sieved before the block loops (pass 1), and med + sparse are sieved per
  block before the scan of that block.
- Pass 1 and the resieve read the same root values (only the update step
  changes them), so the resieve visits exactly the sieve's hit set.
- All divisions in the direct resieve are guarded.
- Finalization runs before nextB.

Results of the unfused build (WSL gcc 11, i5-11600K, -threads 1, runs 3-6
at a time in parallel, so the timings are noisy):

| input | NB | sparse+buckets | list | direct (unfused) |
|---|---|---|---|---|
| 261-bit | 8 | 46.4 s | 134.3 s | 58.3 s |
| 256-bit | 8 | 37.1 s | 104.0 s | 45.1 s |
| 261-bit | 16 | 54.2 s | — | 61.6 s |
| 261-bit | default (6) | 38.8 s | — | 48.5 s |

All runs factored correctly with zero errors (NB 12 also checked). Relation
yield per poly matches stock (261-bit NB=8: 62985 vs 62475 polys, 14093 vs
14015 fulls), so the 15-25% gap is per-poly cost: two full passes over
[large_B, B) per poly against the bucket path's one.

## 6. Fused pass (implemented 2026-10-20)

Observation: once both sides of poly k are scanned, the sieve bytes are
dead, since finalization only reads the pool. So the interval can be
re-initialized for poly k+1 right away, and one loop over [large_B, B) can
do, per prime:

1. **resieve** poly k's roots against the pool (as pass 2);
2. **update** the roots to poly k+1 (as pass 2);
3. **sieve** poly k+1's roots into the freshly memset intervals (as pass 1).

New per-poly pipeline:

```
direct_primed = 0;                                   // per process_poly call (per 'a')
for (numB = 1; numB < maxB; numB++, tot_poly++)
{
    if (!direct_primed)                              // first poly of this 'a'
    {
        memset(sieve, blockinit, 2 * NB * 32768);
        direct_sieve_largeprimes_32k(sieve, ...);    // standalone pass 1
    }
    if (minblock < NB) both sides' minblock -= 4;    // fact 2.4: order-independent
    dpool->num = 0;

    ... pos / neg block loops unchanged ...

    next = (numB + 1 < maxB);
    if (next) memset(sieve, blockinit, 2 * NB * 32768);   // bytes are dead now
    direct_resieve_update_largeprimes_32k(sconf, dconf, next ? sieve : NULL);
    direct_primed = next;
    direct_finalize_survivors(...);                  // pool only; sieve untouched
    nextB; nextRoots_32k_avx2_small;
}
```

`direct_resieve_update_largeprimes_32k(sconf, dconf, uint8_t* sieve)`:
- `sieve != NULL`: resieve, then update, then sieve the new roots
  (pos r, neg p - r, same masks and logp approximation as pass 1).
- `sieve == NULL` (last poly of the 'a'): resieve only. Skip the update and
  the sieve, since nothing consumes them: process_poly returns and the next
  'a' starts with firstRoots.

Notes:
- minblock depends on poly k+1's b, which is only known after nextB, so it
  is applied at the top of the next poly as `-= 4` over that block on both
  sides (fact 2.4) instead of in the memset.
- Nothing between the fused pass and poly k+1's block loops touches the
  sieve buffer (finalize, nextB and nextRoots only touch the pool, poly and
  roots; SS search is excluded).
- The resieve for poly k reads the roots before the update overwrites them
  in the same iteration, so the hit-set invariant (section 5) holds.
- Per poly, [large_B, B) is walked once (plus one extra pass 1 per 'a'),
  matching the bucket path's single root-update pass.

### 6.1 Implementation notes
- The minblock `-= 4` is an explicit AVX512BW loop (`_mm512_sub_epi8`).
  gcc 11 at -O2 does not auto-vectorize; the scalar byte loop cost ~4%
  of total time. So TRY_DIRECT_SIEVE now also requires USE_AVX512BW.
- **Bitmap pre-filter by gather.** A perf profile of the fused build showed
  `direct_resieve_hit` at ~10% of runtime: nearly every lane hits the
  interval (~28 root hits per 16 primes), and each hit was a scalar call
  that almost always missed the bitmap. The fused loop now gathers the
  32-bit bitmap word for each lane's position
  (`_mm512_mask_i32gather_epi32`, index pos>>5) and tests bit pos&31.
  Neg-side positions are computed as `p - r + interval`, masked
  `< 2*interval`. Only lanes with a set bit reach the scalar
  binary-search + guarded division. This was the big win.

### 6.2 Results
Sequential single runs (one job at a time), WSL gcc 11, i5-11600K
(Rocket Lake: 48 KB L1D, 512 KB L2/core), -threads 1. QS elapsed seconds:

| input | NB | sparse+buckets | unfused | fused | fused + gather |
|---|---|---|---|---|---|
| 256-bit | 8 | 26.7 / 26.6 / 27.2 / 26.8 | 31.2 / 31.3 / 31.6 / 31.8 | 32.7 / 32.4 (scalar minblock), 31.9 / 31.7 | **27.4 / 27.4** |
| 261-bit | 8 | 35.1 | — | — | 36.2 |
| 261-bit | 12 | 40.0 | — | — | 43.3 |
| 261-bit | 16 | 46.2 | — | — | 49.8 |

- Fusing by itself bought ~nothing (unfused ~= fused): the second walk over
  [large_B, B) was cheap. The cost was the per-hit scalar bitmap probe.
- With the gather filter, direct is within ~1-3% of sparse+buckets at
  NB=8 and ~8% behind at NB=12/16. Larger NB means a larger scatter target:
  both sides are 2*NB*32 KB = 512 KB / 768 KB / 1 MB, against a 512 KB L2.
  That fits the L2-residency hypothesis: the bucket path's writes are
  sequential, and its per-block application stays in L1.
- Correctness: all of 256-bit NB8, 261-bit NB 6/8/12/16 factor with the
  correct P39s and zero errors. The poly count is identical to the unfused
  build (261-bit NB8: 62985), so the fused pass's hit set matches.

### 6.3 Block-edge bug (found at small NB, fixed)

At NB=2/4, direct needed noticeably more polys than stock: 261-bit NB2
partials/poly were −7% and DLP attempts −13%, with far more "outside
range" residues. Diagnosis:
- `DEBUG_MISSED_FB_PRIMES` showed no missing divisions.
- `DEBUG_DIRECT_ROOTS` showed every `[large_B, B)` root was correct.
- A report dump with a fixed seed (`DEBUG_DUMP_REPORTS`, `-seed 12345
  -noopt`), diffed between builds on identical polys, found the junk
  reports: almost all at **location 32767 of every block that has a
  predecessor**, with bytes ~162-166 (a large bogus subtraction).

Cause: for a-primes <= 8192, `set_aprime_roots` keeps the prime and marks
the roots 0xffff. The unrolled small-prime med-sieve paths then subtract
logp at `blk_sieve[65535]`. In stock's single-block window that byte is
guard space. In the contiguous direct layout it is the **next block's last
byte** (for the last neg block it is guard space). Fix (SIQS.c, both
block loops): save `blk_sieve[65535]` before the block's med/sparse sieve
and restore it afterwards.

(The med sieve also dumps logs for a-primes > 8192 at `blk_sieve[0]`;
that happens in both builds, inside the block, and is unchanged.)

### 6.4 Results after the fix, and logp calibration

Yield per poly after the fix (261-bit, `-noopt` so both builds use the same
tf_small_cutoff; per poly):

| build | NB | polys | partials | survivors | outside range | prp + DLP att |
|---|---|---|---|---|---|---|
| sparse | 2 | 203.5k | 1.133 | 10.14 | 6.24 | 3.02 |
| direct | 2 | 210.6k | 1.036 | 9.47 | 6.24 | 2.34 |
| direct, bias 1 | 2 | 186.9k | 1.466 | 12.48 | 6.52 | 5.07 |
| sparse | 4 | 120.4k | 1.853 | 17.83 | 11.66 | 4.65 |
| direct | 4 | 123.9k | 1.739 | 16.93 | 11.49 | 3.93 |
| direct, bias 1 | 4 | 111.9k | 2.311 | 21.03 | 11.99 | 7.51 |

- The junk is gone: "outside range" now matches stock.
- The remaining gap is **sieve-credit calibration, not correctness**. On
  identical polys, stock's bytes are on average ~1 lower (more credit),
  because each bucket slice applies the logp of the slice's *last* prime to
  every prime in the slice, while direct uses the 16-group's top logp.
  Stock's extra near-threshold reports (bytes 254/255) are mostly DLP
  candidates.
- `DIRECT_LOGP_BIAS` (qs_impl.h, default 0) adds a constant to the direct
  range's logp. With bias 1, direct overtakes stock in yield per poly by
  25-30%. The yield is very sensitive to large-prime credit (equivalently,
  to the report threshold), which is a lever for the parameter exploration.

Wall time (sequential single runs, normal cutoff tuning on, QS seconds):

| input | NB | sparse+buckets | direct | direct bias 1 |
|---|---|---|---|---|
| 261-bit | 2 | 32.8 | 33.8 | 33.4 |
| 261-bit | 4 | 30.1 | 29.8 | 29.7 |
| 261-bit | 8 | 34.6 | 34.3 | 34.4 |
| 256-bit | 2 | 25.6 | 25.9 | 25.9 |
| 256-bit | 4 | 23.5 | 23.3 | 23.6 |
| 256-bit | 8 | 26.7 | 26.4 | 27.0 |

- With the fix, direct is even with or slightly faster than stock at NB 4
  and 8, and ~1-3% slower at NB2. Before the fix it was +3% at NB8 on the
  261-bit input (36.2 s).
- **NB=4 is the fastest setting for both builds** on these 77-79 digit
  inputs: ~12% faster than NB=8 and ~8% faster than NB=2.
- Bias 1's higher yield per poly is cancelled out in wall time by the extra
  trial division (and by the adaptive tf_small_cutoff). Calibration should
  be explored together with the thresholds, not on its own.

Diagnostics left in the code (all off by default):
- `DEBUG_MISSED_FB_PRIMES` (tdiv.c): counts fb primes left in finalized Qs,
  by range. Both builds leave a few (~0.2-0.4% of reports). This predates
  TRY_DIRECT_SIEVE.
- `DEBUG_DIRECT_ROOTS` (SIQS.c): checks Q(r) == 0 mod p for every direct
  root, each poly.
- `DEBUG_DUMP_REPORTS` (tdiv_scan.c): dumps every report with its sieve
  byte and survival flag for the first 20000 scan calls, keyed by call
  index. Use with `-seed N -noopt` to diff builds on identical polys.

### 6.5 Possible next steps
- Profile fused+gather again to split the fused loop's cost between the
  root update, the next-poly scatter and the gathers.
- Shrink the scatter target at larger NB: e.g. sieve the next poly's pos
  and neg halves in two sub-passes (two 256 KB targets at NB=8), or cap the
  direct range's NB dependence.
- Compare on a CPU with a larger L2 (e.g. 1-2 MB/core) where 2*interval
  fits.
- Parameter exploration at small NB (2/4): report threshold / tf_closnuf,
  tf_small_cutoff, DIRECT_LOGP_BIAS, and the large-prime multipliers (the
  single large prime cutoff is currently 90 * pmax).

## 7. Not done / caveats

- USE_SS_SEARCH, batch polys and non-AVX512 builds are excluded at compile
  time; TLP (use_dlp >= 2) is untested (the BW scan branch works through
  the dconf->sieve repointing).
- Pool and bitmap memory: ~MAX_SMOOTH_PRIMES*4 bytes per survivor (grows on
  demand) + interval/4 bytes of bitmap. That is small next to the original
  design's ~26 MB per thread.
- Sieve memory: (2*NB + 4) * 32 KB per thread.
