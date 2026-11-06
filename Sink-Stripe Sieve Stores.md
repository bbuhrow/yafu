# Sink-Stripe Sieve Stores

The "miss stripe" from Math-Prime-Util-GMP's SIQS, written up so it can be tried in YAFU's large-prime sieve loop.

Nothing here has been benchmarked in YAFU. Section 6 suggests an experiment.

## 1The current YAFU pattern

YAFU compares a vector of roots (16 with AVX-512, 8 with AVX2) against the interval size, then loops over the set bits:

```
while (m1 > 0) {
    pos = _trail_zcnt(m1);
    sieve[r_1[pos]] -= logp;
    m1 = _reset_lsb(m1);
}
```

That costs one read-modify-write per hit, plus loop control. The trip count is `popcount(m1)`, which depends on the data. When roots hit at a moderate rate (roughly 10–90%), the count changes from vector to vector and the loop-exit branch mispredicts often, at about 15–20 cycles each on current x86.

## 2The MPU-GMP technique

From `siqs_sieve_kernel_inc.c`: `siqs_sieve_one_hit_local` and the one-hit tier of `siqs_run_sieve_kernel`, near line 216.

```
uint32_t sink = length;                  /* first byte past the interval */
for (; i < end; i++) {
    pos = root1[i] < length ? root1[i] : sink;      /* cmov, not a branch */
    sieve[pos] += logp[i];
    if (root1[i] != root2[i]) {
        pos = root2[i] < length ? root2[i] : sink + 1;
        sieve[pos] += logp[i];
    }
    sink += 2;
    if (sink == length + 4096) sink = length;       /* rotate through 4 KiB */
}
```

- **Every root is written.** The compare chooses the address, not whether a store happens. A root that misses by any amount gets the sink address. Its value is never used as an address, so no padding as large as pmax is needed.
- **No data-dependent branches.** gcc 11.4 at `-O3 -march=native` compiles each root to `cmp` → `cmovbe` → `addb (sieve,idx)`. The only branch left is the `root1 != root2` check, which is almost always taken the same way.
- **Why the sink rotates.** If every miss used one byte, each read-modify-write would wait on the previous one through store-to-load forwarding, so the misses would run as one serial chain. Rotating keeps consecutive miss stores on different bytes, and the stripe stays in L1. MPU uses 4096 bytes; a few hundred would probably also break the chain.
- **Sink bytes are garbage and never read.** They sit past the scanned interval and are never cleared. Overflow or wrap there doesn't matter.
- **When MPU uses it.** Only for p > interval length, and only when there are ≥ 5200 such primes. Below that it uses a blind store into a buffer padded to pmax. The 5200 threshold is a safety argument about allocation size, not a tuned speed crossover. The kernels were tuned on an Apple M1 Pro, not x86.

The technique does nothing for hits, which still land at random positions in the interval. It only makes misses cheap and branch-free.

## 3Adapting it to YAFU's vector compare

YAFU already has the compare mask, so the scalar `cmov` becomes a vector blend. The loop over set bits becomes a fixed loop over all lanes, fully unrolled.

```
/* Assumes: interval length L in bytes; sieve buffer has >= L + SINK bytes,
 * with [L, L+SINK) never scanned or relied on.  SINK is a multiple of 32. */
#define SINK 1024
const __m512i lane = _mm512_setr_epi32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
const __m512i vL   = _mm512_set1_epi32(L);
uint32_t sink = L;
ALIGNED64 uint32_t idx[16];

for (i = start; i < end; i += 16) {
    __m512i r1 = _mm512_load_epi32(root1 + i);
    __mmask16 m1 = _mm512_cmplt_epu32_mask(r1, vL);    /* existing compare */
    __m512i vs = _mm512_add_epi32(_mm512_set1_epi32(sink), lane);
    _mm512_store_epi32(idx, _mm512_mask_blend_epi32(m1, vs, r1));
    /* fixed 16 RMWs, no data-dependent branch */
    sieve[idx[0]]  -= logp; sieve[idx[1]]  -= logp; sieve[idx[2]]  -= logp;
    sieve[idx[3]]  -= logp; sieve[idx[4]]  -= logp; sieve[idx[5]]  -= logp;
    sieve[idx[6]]  -= logp; sieve[idx[7]]  -= logp; sieve[idx[8]]  -= logp;
    sieve[idx[9]]  -= logp; sieve[idx[10]] -= logp; sieve[idx[11]] -= logp;
    sieve[idx[12]] -= logp; sieve[idx[13]] -= logp; sieve[idx[14]] -= logp;
    sieve[idx[15]] -= logp;
    sink += 16;
    if (sink >= L + SINK) sink = L;     /* well predicted */
    /* same for root2 with its own sink offsets, e.g. sink+16 .. sink+31 */
}
```

- **Don't use a scatter.** AVX-512 has no byte scatter, and a dword scatter loses updates when two lanes hit the same position. Sequential scalar read-modify-writes handle duplicate positions correctly.
- **Reading the indices.** Reloading the indices after a 64-byte vector store normally forwards fine on Intel and AMD. Extracting lanes with `_mm512_extracti32x4_epi32` and `_mm_extract_epi32` is the alternative if it measures better.
- **Early-out.** Keep an `m1 == 0` early-out only for prime ranges where all-miss vectors are common and predictable. The goal is to remove unpredictable branches.
- **AVX2.** The same pattern works with 8 lanes and `_mm256_blendv_epi8`. AVX2 has only signed compares, which is fine as long as roots and L stay below 2³¹.
- **Padding.** The extra SINK bytes must be allocated and must lie outside any SIMD scan, memset or bucket range that matters.

## 4When it should win and when it shouldn't

Let h be the fraction of roots that hit, about (interval length)/p for primes larger than the interval.

| hit rate h | set-bit loop | sink stripe |
| --- | --- | --- |
| high (≳ 0.5) | \~16h stores, frequent exit mispredicts | 16 stores, no mispredicts; likely better |
| moderate (0.1–0.5) | popcount varies most; mispredicts dominate | likely better |
| low (≲ 1/16) | mostly 0–1 bits; `m == 0` early-out predicts well | 16 wasted stores per vector; likely worse |

The break-even depends on store throughput (1–2 stores per cycle on recent x86) against the mispredict penalty. Using one method for some prime ranges and the other for the rest is reasonable. The sink is most useful for primes just above the interval size, where hits are common but not certain. It is probably not useful for the largest primes, where almost every root misses.

## 5An alternative to compare against

The sink fixes the branch but adds stores. You can instead compact hits from many vectors first: `_mm512_mask_compressstoreu_epi32(buf + n, m1, r1); n += popcount(m1);` over a batch of primes, then one loop over `buf[0..n)`. That leaves one loop-exit mispredict per batch rather than per vector, with no wasted stores. It is the same idea a bucket sieve uses, so it may fit YAFU's existing structure better. It's worth benchmarking it alongside the sink.

## 6Suggested experiment

1. Instrument the current loop: average popcount per vector, broken down by prime range.
2. Implement the sink variant behind a compile-time switch.
3. Time per-polynomial sieve cost at a few sizes (C80, C95, C100) with `perf stat -e branch-misses,cycles`, sink on and off.
4. If the result is mixed, choose per prime range using the measured hit rate.