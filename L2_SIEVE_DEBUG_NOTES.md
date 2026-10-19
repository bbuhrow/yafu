# L2 Sieve Debugging Notes (SIQS, TRY_L2_SIEVE)

*Status as of 2026-10-19. Written to allow resuming the debugging session later.*

## 1. What the L2 sieve is

The modified SIQS in `factor/qs/` adds a "L2 sieve" class of primes under
`TRY_L2_SIEVE` (defined in `include/qs_impl.h:79`). Instead of sieving every
prime block-by-block, primes in the range **(med_B, fb_blk8i)** — approximately
(49152, 262144) — are sieved across the **whole 8-block interval**
(8 x 32768 = 262144 bytes) at once, once per polynomial per sign side. This is
an L2-cache-sized direct interval sieve. All other primes keep the original
block-by-block treatment.

Build config used throughout: `USE_AVX512F + USE_AVX512BW` (SKYLAKEX),
`USE_BATCHPOLY`/`USE_BATCHPOLY_X2` **not** defined (plain per-poly dispatch),
`TRY_L2_SIEVE` defined. Run with `-siqsNB 8` to force exactly 8 blocks.
Inputs: 256-bit and 261-bit (79 digit) semiprimes, `use_dlp = 1`,
`scan_ptr = check_relations_siqs_16_avx512`.

### Component map (with the boundaries that matter)

From a live run (79-digit input, `-v -v`), fb boundaries:

```
        small: 1024   10bit: 88    11bit: 160   12bit: 296   13bit: 528
        14bit: 944    15bit: 1760  med: 2576    L2: 11520    large: 11520
        large_x2: 36176  all: 36176
start primes:
        SPV: 241  10bit: 1033  15bit: 32749  med: 49367
        L2: 262109  large: 262109  large_x2: 907657  all: 907657
```

i.e. `med_B = 2576` (primes < 49367), `fb_blk8i = large_B = 11520`
(primes < 262109; note with NB=8, sieve_interval = 262144 so `large_B ==
fb_blk8i` and `nextRoots_32k_knl_medbucket`'s bucket loop is EMPTY),
`B = 36176`. The fb contains only QR primes, so fb indices are ~half the
integer-prime rank.

Data flow per polynomial per side:

1. **Roots for the L2 range** are computed for the first poly of each 'a' by
   the L2 section in `factor/qs/poly_roots_32k.c:2273-2331` (firstroots +
   rootupdates), and updated per-poly by the L2 sections in
   `nextRoots_32k_avx2_small()` (`factor/qs/update_poly_roots_32k_avx2.c:3262+
   sign>0, 3385+ sign<0`) using `COMPUTE_NEXT_ROOTS_P/N`, which dereference the
   running `ptr` into `dconf->rootupdates` (32-bit int array).
2. **L2 sieve** `L2_sieveblock_32k_avx512bw()` (`factor/qs/med_sieve_32k_avx2.c:968-1118`):
   called once per side from `process_poly()` (`SIQS.c:2243-2248` p-side,
   `2376-2381` n-side) after a single `memset(sieve, blockinit, 32768*8)`.
   Vector loop over 16 primes at a time; fast loop sieves all lanes while
   `root2 < 8*32768` (`vblock`), masked loop handles lanes dropping out, tail
   handles the final root1 hit, then wraps both roots (subtract 8*32768, store
   min/max back into `dconf->sieve_p/n->root1/root2`).
   NOTE: the per-block memsets in `process_poly` are compiled OUT under
   `TRY_L2_SIEVE` (`SIQS.c:2279-2281, 2296-2298, 2410-2412, 2427-2429`).
3. Per-block processing is unchanged: `med_sieve_ptr` (primes < med_B),
   `lp_sieveblock_ptr` (bucket primes), `scan_ptr` =
   `check_relations_siqs_16_avx512` (`factor/qs/tdiv_scan.c:947-1094`, which
   correctly uses `dconf->sieve + blocknum*32768` in its non-BW branch under
   TRY_L2_SIEVE).
4. **tdiv for L2 primes**: `tdiv_L2primes_32k_avx2()`
   (`factor/qs/tdiv_med_32k_avx2.c:854-913`, called from `tdiv_scan.c:1065-1067`).
   For each report at interval position `block_loc = reports[rn] + bnum*32768`,
   it advances the position past the interval end
   (`tmp = 1 + divide(8*32768 - block_loc, inv32/corr32/shift32)`,
   `tmp = block_loc + tmp*prime - 8*32768`) and compares against the **wrapped
   roots** stored by the L2 sieve. Match => `DIVIDE_ONE_PRIME` (unchecked first
   division). This "trust the progression match" design mirrors the bucket tdiv.
5. **Bucket primes** `[fb_blk8i, B)` are bucket-sieved as before; with NB=8 the
   `nextRoots_32k_knl_medbucket` loop is empty (its TRY_L2 branch at
   `update_poly_roots_32k_knl.c:81-85` sets `med_B = fb_blk8i`) and
   `nextRoots_32k_knl_largebucket` continues from the state medbucket hands it.

Key struct facts (`include/qs_impl.h`):
- `dynamic_conf_t`: `Qvals` ptr (657), `valid_Qs` ptr (658),
  `fb_offsets[MAX_SIEVE_REPORTS][MAX_SMOOTH_PRIMES]` = `[2048][200]` (659),
  **`smooth_num` ptr immediately after (660)**.
- `MAX_SIEVE_REPORTS = 2048`, `MAX_SMOOTH_PRIMES = 200`, `BUCKET_ALLOC = 2048`.
- `dconf->sieve` = `xmalloc_align(qs_blocksize * 20)` offset by +2 blocks under
  TRY_L2_SIEVE (`SIQS.c:3388-3395`) — 18 blocks usable above `dconf->sieve`.
- `dconf->sieve_p/n` arrays are sized `fb_blk8i` (`SIQS.c:3295-3311`).
- `DIVIDE_ONE_PRIME` (tdiv_macros_common.h:3-8 and tdiv_large.c:342-347):
  `do { fb_offsets[++smooth_num] = i; mpz_tdiv_q_ui(Q,Q,prime); } while
  (mpz_tdiv_ui(Q,prime)==0);` — **no divisibility check before the first
  division, no bound check on smooth_num.**

## 2. Tooling built for this investigation (all under /opt/opencode/tmp/795601290/opencode/)

- `guardalloc.c` / `guardalloc.so` — LD_PRELOAD allocator giving every heap
  chunk its own mmap with a PROT_NONE page immediately after. Converts any
  page-crossing OOB heap write into an immediate SIGSEGV at the writing
  instruction. (Smoke test: `gt3.c`; note -O2 dead-store-eliminates deliberate
  OOB writes in tests — use -O0.)
- `gdbcmds2.txt` — breakpoint
  `break tdiv_large.c:639 if __gmpz_tdiv_ui(dconf->Qvals[report_num], prime) != 0`
  catches the **first division of Q by a non-divisor** (the corruption entry
  point). Note: use `__gmpz_tdiv_ui` (gdb can't resolve `mpz_tdiv_ui`).
- `gdbcmds6.txt` — richer capture: Q sign/position variant tests via
  `__gmpz_*` calls in gdb (needs explicit return-type casts), bucket dump for
  an element across all blocks.
- `N.txt` — 256-bit test semiprime
  82681505245784387573134107388964238632969984682532062611426448367569436877569
  (factors in `factors.txt`). Crash reproduced deterministically in ~15 s with
  `-threads 1`.
- Repro command:
  `echo $N | env LD_PRELOAD=.../guardalloc.so gdb -batch -x gdbcmds2.txt --args ./yafu "siqs(@)" -siqsNB 8 -threads 1`
- Build notes: repo builds with `Makefile` + `config.mk`. GMP 6.2.0 can be
  built out-of-tree from the read-only clone `.slim/clonedeps/repos/gnu__gmp`.
  No ASan/valgrind on this box (no libasan.so.5 runtime); the guard allocator
  substitutes. `timeout` on PATH is a broken python shim — use /usr/bin/timeout.

## 3. Bugs found and their status

### Bug 1 — poly_roots_32k.c first-poly bucket setup not L2-aware  [FIXED, verified]

`factor/qs/poly_roots_32k.c:1964` (`bound_val = fb->med_B`) and `:1987`
(`lp_bucket_p->fb_bounds[0] = fb->med_B`) initialized the first-poly/new-'a'
bucket bookkeeping from the non-L2-aware `med_B`, while the L2-adjusted bucket
loop (line 2333+) scatters with `bound_val = fb->fb_blk8i`.

Effect: slice-0 bucket elements encoded relative to fb_blk8i but decoded via
`i = fb_bound + (e>>16)` with `fb_bound = med_B` — every slice-0 entry decoded
to a wrong factor-base index (off by fb_blk8i - med_B). `tdiv_LP`'s unchecked
DIVIDE_ONE_PRIME then divided Q by non-factors.

Fix (applied by user): use `fb->fb_blk8i` under TRY_L2_SIEVE for both lines.
Verified by capture: `fb_bounds[0] = 11456` after fix.

### Bug 2 — stale `ptr` in nextRoots_32k_avx2_small's L2 sections  [FIXED, verified]

`factor/qs/update_poly_roots_32k_avx2.c` — `nextRoots_32k_avx2_small()`: `ptr`
(= `&rootupdates[(v-1)*bound + startprime]`, line 3174) is advanced in
lockstep with `j` only by the two scalar loops (3179, 3205); the second breaks
at the first 16-aligned j WITHOUT advancing ptr; the AVX2 16x block then
processes [j, med_B) using the separate `sm_ptr` and only `j` jumps to med_B.
In stock yafu nothing afterwards used `ptr`. The new L2 sections
(line ~3262 sign>0, ~3385 sign<0) run
`for (j = med_B; j < fb_blk8i; j++, ptr++)` with `COMPUTE_NEXT_ROOTS_P/N`,
which dereference `*ptr` (poly_macros_common.h:708+; all active variants use
`*ptr`, the 32-bit rootupdates array).

Effect: for every poly >= 1, every L2-range prime got its root updated with
the update value of a different, much smaller prime (~2400 entries off) ->
wrong L2 sieve positions AND wrong wrapped roots -> `tdiv_L2primes` false
progression matches -> unchecked division of Q by non-factors -> Q reaches 0
-> every later prime "divides" (0 % p == 0) -> smooth_num explodes ->
`fb_offsets[++smooth_num]` overflows the report's row, then the whole
fb_offsets array inside dynamic_conf_t, smashing the adjacent `smooth_num`
pointer (qs_impl.h:660) and cascading into `Qvals` -> SIGSEGV (in
`__gmpn_divrem_1` via corrupted Qvals in the original backtrace; at the
fb_offsets store in the guard-page build, `mov %ebp,(%rax,%r14,4)` with
smooth_num = 409969).

Fix (applied by user): re-point before both L2 loops:
`ptr = &rootupdates[(v - 1) * bound + sconf->factor_base->med_B];`

After fixes 1-3 the run gets much farther (20234 relations found on the
79-digit input before the next crash — see Bug 4).

### Bug 3 — fb_blk8i rounds UP, admitting primes > 262144 into the L2 range  [FIXED, verified]

The fb_blk8i loop (`SIQS.c:4437-4456`) breaks at the first index that is both
16-aligned AND has prime > 8*blocksize — rounding the boundary UP. Live data:
primes 262147, 262303, 262459 at indices fb_blk8i-15 .. fb_blk8i-1 were inside
the L2 range. For those, roots can exceed 262144, and
`L2_sieveblock_32k_avx512bw`'s fast loop initializes `valid_mask = 0xffff`
unconditionally (med_sieve_32k_avx2.c:997) — its first iteration sieves all 32
root positions with no range check, writing up to ~315 bytes past the 8-block
region (currently lands in the unused tail of the 20-block sieve allocation,
i.e. silent, but latent OOB), and the wrapped roots for those primes violate
the routine's own < p contract going into tdiv_L2primes.

Fix (applied by user): round the boundary DOWN
(`fb_blk8i = <first index with prime > 8*blocksize> & ~15;`), plus harden the
fast loop: `valid_mask = _mm512_cmplt_epu32_mask(vroot2, vblock);` instead of
0xffff. (The invariants hold: with root1 <= root2 < prime, the masked loop +
tail handle partial validity and the final wrap cannot underflow. Verified by
the post-fix fb bounds printout: L2 ends at 262109 < 262144.) The 15 displaced
primes correctly fall into the bucket range (p > interval => <= 1 hit/root).

### Bug 4 — NEW, UNRESOLVED: crash at `dconf->smooth_num[report_num] = smooth_num;`  [OPEN]

After fixes 1-3, on the 79-digit input
86601586107718322788075198193322239264919417729992471551587341952662571831281
(261 bits, FB 36176 primes, 8 blocks, 1 thread), the run found **20234
relations** (10195 full + 10039 from 224165 partial) before crashing:

```
Program received signal SIGSEGV, Segmentation fault.
0x0000000000497c7b in tdiv_LP_avx2 (report_num=2, parity=<optimized out>,
    bnum=<optimized out>, sconf=0x30f0d50, dconf=0x3121710) at factor/qs/tdiv_large.c:670
670             dconf->smooth_num[report_num] = smooth_num;
#1 check_relations_siqs_16_avx512 (blocknum=6, ...) at factor/qs/tdiv_scan.c:1080
#2 process_poly (...) at factor/qs/SIQS.c:2350
#3 SIQS (...) at factor/qs/SIQS.c:946
```

Analysis (not yet confirmed): line 670 is a plain store through the
`dconf->smooth_num` POINTER. Faulting there means the pointer field itself is
corrupt — and `smooth_num` (qs_impl.h:660) is the field **immediately after**
`fb_offsets[2048][200]` (line 659) in dynamic_conf_t. So the leading
hypothesis is a **residual fb_offsets overflow** (a report whose smooth_num
exceeds MAX_SMOOTH_PRIMES=200, or row-2047 overflow past the array) — i.e.
some remaining unchecked-division / runaway path, of the same class as bugs
1-2 but from a different trigger. Candidate sources to investigate when
resuming:

1. Rerun the first-bogus-division breakpoint (gdbcmds2.txt) on the current
   binary — if it fires, capture which stage is dividing Q by a non-divisor
   now (tdiv_L2primes with the new boundary? resieve? tdiv_med?).
2. Audit every `fb_offsets[++smooth_num]` site for missing
   `smooth_num < MAX_SMOOTH_PRIMES` guards (tdiv_macros_common.h, tdiv_large.c,
   tdiv_med*, tdiv_resieve*, tdiv_small.c).
3. Defensive fix worth making regardless: bound smooth_num in all
   DIVIDE_ONE_PRIME / DIVIDE_RESIEVED_PRIME macros.
4. Check whether the crash reproduces with the guard allocator (which turns
   the fb_offsets array-end overflow into a precise fault at the writing
   instruction) and/or with MALLOC_PERTURB_ set.
5. Note the crash blocknum=6 and report_num=2; check dconf->numB at capture.

## 4. Secondary issues (real, but not crashers — fix eventually)

1. **tdiv_L2primes exact-boundary edge**: when a root lands exactly on
   8*32768, the L2 sieve correctly stores wrapped root 0 (position 0 of the
   next interval), but tdiv_L2primes' `1 + divide(8*32768 - block_loc, ...)`
   computes `prime` instead of 0 for reports where (8*32768 - block_loc) is a
   multiple of the prime -> that prime silently isn't divided out -> occasional
   lost relations. (No false positives: match still implies loc ≡ root mod p.)
2. **filter_SPV** (`factor/qs/tdiv_small.c:189`): reads
   `sieve[dconf->reports[report_num]]` with no `blocknum * 32768` offset —
   under L2 it validates every report against block 0's sieve bytes. The
   blocknum param is available; fix is `sieve + blocknum*32768` at the call
   site in tdiv_scan.c or inside filter_SPV.
3. **check_relations_siqs_16_avx512 AVX512BW branch** (`tdiv_scan.c:964-993`,
   taken when use_dlp >= 2, i.e. TLP inputs): scans `dconf->sieve` with NO
   blocknum offset under TRY_L2_SIEVE (the non-BW else branch at 1013-1019 has
   it). Latent until a triple-large-prime run.
4. **Hardcoded 8 blocks**: `memset(sieve, blockinit, 32768*8)` (SIQS.c:2245,
   2378), `vblock = 8*32768` in L2_sieveblock, and the L2 comments assume 8.
   With `-siqsNB` != 8 or the default num_blocks for larger inputs (param
   table goes up to 22), blocks >= 8 are never initialized/L2-sieved (per-block
   memsets are compiled out under TRY_L2_SIEVE) and num_blocks > 18 would write
   past the 20-block sieve allocation. Should be `sconf->num_blocks`-driven or
   guarded.
5. `L2_sieveblock_32k_avx512bw`'s scalar tail guard
   (`if ((fb_blk8i - i) < 16) bound = fb_blk8i;`) is fine, but note the loop
   assumes `fb_blk8i > med_B` and `med_B` 16-aligned (SIQS.c:4429, 4448 —
   still true after the round-down change).

## 5. Misc facts established (so we don't re-derive)

- `use_dlp = 1` at >= 70 digits (SIQS.c:4665); scan_ptr = check_relations_siqs_16_avx512.
- num_blocks defaults: param table (siqs_aux.c:458+); NB=8 forced via -siqsNB 8.
- blockinit = tf_closnuf (SIQS.c:5037); for 77-79 digits = digits_n+7-2 (AVX512F)
  = 82-84. Sieve bytes wrap past 0 into 128-255 when hit hard; the scan's
  hi-bit test detects reports this way. filter_SPV: `bits = (255 - bits) +
  tf_closnuf + 1`.
- medbucket's bucket loop under L2+NB=8 is EMPTY (large_B == fb_blk8i); the
  largebucket "at most one hit per interval" method covers all bucket primes —
  consistent, since all bucket primes > 262144 = interval.
- inv32/corr32/shift32 are initialized for ALL fb primes (siqs_aux.c:258,337).
- sm_rootupdates is uint16 and only valid below med_B; the L2 sections must use
  the 32-bit rootupdates (they do, via *ptr — once ptr is correct).
- The tdiv_LP bucket element decode is `i = slicebound[k] + (e >> 16)`, encode
  `(j - bound_val) << 16 | (root & 32767)`; CHECK_NEW_SLICE records bounds at
  index >= 1 (poly_macros_common.h:55-96); fb_bounds[0] belongs to the bucket
  function's init.
- Heap layout relevant to the corruption cascade: dconf->reports (8KB, malloc
  at SIQS.c:3747) immediately precedes dconf->Qvals (32KB, :3749). The scan's
  report loop (tdiv_scan.c:1021-1052) is unbounded, clamped only after
  (:1055) — another latent overflow if a block ever floods.
