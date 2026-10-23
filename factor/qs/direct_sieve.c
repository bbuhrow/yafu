/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// Portable (any ISA) parts of the large-prime strategies (see
// direct_sieve.md and plans/siqs_multi_isa/plan.md):
//  - the direct sieve's survivor pool, resieve hit test and finalization
//  - scalar kernels for direct pass 1, the fused resieve/update/scatter
//    pass, the sparse sieve and its tdiv, with the same kernel choices as
//    the AVX-512 versions (bitloop = branchy, sink = cmov to a sink stripe,
//    creg/cmem = branch-free scalar append + batched sieve).  They update
//    the sieve and divide in the same order as the AVX-512 kernels, so all
//    kernels do bit-identical work.
//  - nextRoots_32k_lp_portable: per-poly root updates for [med_B, B) under
//    any strategy (sparse roots, mid/large bucket fills, nothing for the
//    direct range, whose roots the fused pass advances).

#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "tdiv_macros_common.h"
#include "poly_macros_32k.h"
#include "poly_macros_common.h"

#if defined(USE_SSE41) || defined(USE_AVX2) || (defined(__SSE2__) && !defined(FORCE_GENERIC))
#include <immintrin.h>
#endif

#if defined(TRY_DIRECT_SIEVE)

// ----------------------------------------------------------------------
// survivor pool
// ----------------------------------------------------------------------
// Scan survivors of every block of both sides are parked here until the
// [large_B, B) resieve (pass 2) has divided them, then finalized.

direct_pool_t* direct_pool_alloc(uint32_t num_blocks)
{
    direct_pool_t* pool = (direct_pool_t*)xmalloc(sizeof(direct_pool_t));
    uint32_t i;

    pool->num = 0;
    pool->active = 0;
    pool->ab_saved = 0;
    pool->ab_miss = 0;
    pool->ex_polys = pool->ex_resolved = 0;
    memset(pool->ex_hist, 0, sizeof(pool->ex_hist));
    memset(pool->ex_qbits, 0, sizeof(pool->ex_qbits));
    pool->ex_loc0 = 0;
    pool->ex_skip = pool->ex_groups = 0;
    pool->ab_samp = pool->ab_sampdiv = 0;
    pool->ab_num[0] = pool->ab_num[1] = 0;
    pool->ab_wrong[0] = pool->ab_wrong[1] = 0;
    pool->alloc = 1024;
    pool->interval = num_blocks << 15;
    pool->pos = (uint32_t*)xmalloc(pool->alloc * sizeof(uint32_t));
    pool->Q = (mpz_t*)xmalloc(pool->alloc * sizeof(mpz_t));
    for (i = 0; i < pool->alloc; i++)
        mpz_init(pool->Q[i]);
    pool->smooth_num = (int*)xmalloc(pool->alloc * sizeof(int));
    pool->act = (uint8_t*)xmalloc(pool->alloc * sizeof(uint8_t));
    pool->ab_est = (int16_t*)xmalloc(pool->alloc * sizeof(int16_t));
    pool->ab_rem = (int16_t*)xmalloc(pool->alloc * sizeof(int16_t));
    pool->ab_byte = (uint8_t*)xmalloc(pool->alloc * sizeof(uint8_t));
    memset(pool->ab_bytes, 0, sizeof(pool->ab_bytes));
    memset(pool->ab_allq, 0, sizeof(pool->ab_allq));
    pool->ab_junk = pool->ab_junk_last = pool->ab_junk_first = 0;
    memset(pool->ab_junk_blk, 0, sizeof(pool->ab_junk_blk));
    memset(pool->ab_junk_loc, 0, sizeof(pool->ab_junk_loc));
    pool->ab_junk_side[0] = pool->ab_junk_side[1] = 0;
    pool->ab_sn0 = (int16_t*)xmalloc(pool->alloc * sizeof(int16_t));
    memset(pool->ab_hist, 0, sizeof(pool->ab_hist));
    pool->fb_offsets = (uint32_t*)xmalloc(
        (size_t)pool->alloc * MAX_SMOOTH_PRIMES * sizeof(uint32_t));

    // one bit per position of both intervals; must start (and is kept) zero
    pool->bitmap = (uint64_t*)xmalloc(2 * pool->interval / 8);
    memset(pool->bitmap, 0, 2 * pool->interval / 8);

    // the hit list is sized from the factor base (direct_hits_reserve)
    pool->hits = NULL;
    pool->hit_gend = NULL;
    pool->hits_alloc = 0;

    return pool;
}

// make room for the hit list of [large_B, B): at most 4 entries per prime,
// plus room for the compress kernels' full-vector stores past the end
void direct_hits_reserve(direct_pool_t* pool, static_conf_t* sconf)
{
    fb_list* fb = sconf->factor_base;
    uint32_t need = 4 * (fb->B - fb->large_B) + 96;

    if (pool->hits_alloc >= need)
        return;

    if (pool->hits_alloc > 0)
    {
        align_free(pool->hits);
        free(pool->hit_gend);
    }
    pool->hits = (uint32_t*)xmalloc_align(need * sizeof(uint32_t));
    pool->hit_gend = (uint32_t*)xmalloc(((fb->B - fb->large_B) / 16 + 1) * sizeof(uint32_t));
    pool->hits_alloc = need;
}

void direct_pool_free(direct_pool_t* pool)
{
    uint32_t i;

#if defined(DIRECT_EXIT_STATS)
    if (pool->ex_polys > 0)
        printf("direct exit: %" PRIu64 " polys, %" PRIu64 " with all survivors resolved\n",
            pool->ex_polys, pool->ex_resolved);
    if (pool->ex_polys > 0)
    {
        for (i = 0; i < 16; i++)
            printf("  exhist %2u: %" PRIu64 "\n", i, pool->ex_hist[i]);
        for (i = 0; i < 256; i++)
            if (pool->ex_qbits[i])
                printf("  exqbits %3u: %" PRIu64 "\n", i, pool->ex_qbits[i]);
        printf("  exloc0 %" PRIu64 "\n", pool->ex_loc0);
        printf("  exskip %" PRIu64 " of %" PRIu64 " resieve groups skipped\n", pool->ex_skip, pool->ex_groups);
    }
#endif
#ifdef DIRECT_ABANDON_STATS
    if (pool->ab_saved > 0)
        printf("direct abandon: %" PRIu64 " survivors; nolarge %" PRIu64
            " (%" PRIu64 " wrong), hopeless %" PRIu64 " (%" PRIu64 " wrong)\n",
            pool->ab_saved, pool->ab_num[0], pool->ab_wrong[0],
            pool->ab_num[1], pool->ab_wrong[1]);
#if defined(DIRECT_ABANDON_NOLARGE)
    if (pool->ab_saved > 0)
    {
        printf("direct abandon: %" PRIu64 " bitmap hits not found in pos[]\n", pool->ab_miss);
        printf("direct abandon: sampled %" PRIu64 " big-error no-hit survivors, %" PRIu64 " divisible by a range prime\n",
            pool->ab_samp, pool->ab_sampdiv);
        for (i = 0; i < 256; i++)
            if (pool->ab_bytes[i])
                printf("  sminusq %4d: %" PRIu64 "\n", (int)i - 128, pool->ab_bytes[i]);
        for (i = 0; i < 256; i++)
            if (pool->ab_allq[i])
                printf("  allsminusq %4d: %" PRIu64 "\n", (int)i - 128, pool->ab_allq[i]);
        printf("  blockinit %u\n", (uint32_t)pool->ab_init);
        printf("  junk %" PRIu64 " loc32767 %" PRIu64 " loc0 %" PRIu64 " side0 %" PRIu64 " side1 %" PRIu64 "\n",
            pool->ab_junk, pool->ab_junk_last, pool->ab_junk_first, pool->ab_junk_side[0], pool->ab_junk_side[1]);
        for (i = 0; i < 16; i++)
            printf("  junkblk %2u: %" PRIu64 "  junkloc %5u: %" PRIu64 "\n", i, pool->ab_junk_blk[i],
                i << 11, pool->ab_junk_loc[i]);
        printf("direct abandon: est large logp histogram (est: no-resieve-hit / hit)\n");
        for (i = 0; i < 64; i++)
            if (pool->ab_hist[0][i] | pool->ab_hist[1][i])
                printf("  est %2u: %" PRIu64 " / %" PRIu64 "\n", i,
                    pool->ab_hist[0][i], pool->ab_hist[1][i]);
    }
#endif
#endif

    for (i = 0; i < pool->alloc; i++)
        mpz_clear(pool->Q[i]);
    free(pool->Q);
    free(pool->pos);
    free(pool->smooth_num);
    free(pool->act);
    free(pool->ab_est);
    free(pool->ab_rem);
    free(pool->ab_byte);
    free(pool->ab_sn0);
    free(pool->fb_offsets);
    free(pool->bitmap);
    if (pool->hits_alloc > 0)
    {
        align_free(pool->hits);
        free(pool->hit_gend);
    }
    free(pool);

    return;
}

static void direct_pool_grow(direct_pool_t* pool)
{
    uint32_t i, newalloc = 2 * pool->alloc;

    pool->pos = (uint32_t*)xrealloc(pool->pos, newalloc * sizeof(uint32_t));
    pool->Q = (mpz_t*)xrealloc(pool->Q, newalloc * sizeof(mpz_t));
    for (i = pool->alloc; i < newalloc; i++)
        mpz_init(pool->Q[i]);
    pool->smooth_num = (int*)xrealloc(pool->smooth_num, newalloc * sizeof(int));
    pool->act = (uint8_t*)xrealloc(pool->act, newalloc * sizeof(uint8_t));
    pool->ab_est = (int16_t*)xrealloc(pool->ab_est, newalloc * sizeof(int16_t));
    pool->ab_rem = (int16_t*)xrealloc(pool->ab_rem, newalloc * sizeof(int16_t));
    pool->ab_byte = (uint8_t*)xrealloc(pool->ab_byte, newalloc * sizeof(uint8_t));
    pool->ab_sn0 = (int16_t*)xrealloc(pool->ab_sn0, newalloc * sizeof(int16_t));
    pool->fb_offsets = (uint32_t*)xrealloc(pool->fb_offsets,
        (size_t)newalloc * MAX_SMOOTH_PRIMES * sizeof(uint32_t));
    pool->alloc = newalloc;

    return;
}

// 1 if Q >= p^2.  Q has no fb factors below p when this is asked, so
// otherwise Q is 1 or prime: dividing it by p or a larger fb prime can
// only turn a large prime into a full relation.  p^2 needs 64 bits and
// unsigned long may be 32, hence the double (exact: p < 2^26).
static int direct_q_open(mpz_t Q, uint32_t p)
{
    return mpz_cmp_d(Q, (double)p * (double)p) >= 0;
}

#if defined(DIRECT_ABANDON_NOLARGE)
// 1 if the sieve byte of this survivor shows no [large_B, B) prime hit:
// the total logp sieved here, less that of the fb primes below large_B
// that tdiv found dividing Q (once per prime, as sieved; tiny primes
// aren't sieved), is less than one range prime's logp.  The minblock's
// lower init, the per-16 logp approximation and powers only inflate the
// remainder, i.e. err toward keeping the survivor active.
static int direct_large_logp(uint32_t report_num, static_conf_t* sconf,
    dynamic_conf_t* dconf)
{
    fb_list* fb = sconf->factor_base;
    uint32_t* offs = dconf->fb_offsets[report_num];
    int smooth_num = dconf->smooth_num[report_num];
    int init = sconf->blockinit;
#if defined(DIRECT_ABANDON_STATS)
    dconf->dpool->ab_init = sconf->blockinit;
#endif
    int v = dconf->sieve[dconf->reports[report_num]];
    int sieved, j;
    uint32_t prev = 0xffffffff;

    // reports have the high bit set, i.e. init - sieved wrapped below 0
    if (init >= 128)
        return 255;
    if (v > init)
        v -= 256;
    sieved = init - v;

    for (j = 0; j <= smooth_num; j++)
    {
        uint32_t i = offs[j];
        if ((i != prev) && (i >= sconf->sieve_small_fb_start) && (i < fb->large_B))
            sieved -= (int)fb->list->logprime[i];
        prev = i;
    }

    return sieved;
}

static int direct_no_large_hit(int large_logp, static_conf_t* sconf)
{
    return large_logp < ((int)sconf->factor_base->list->logprime[
        sconf->factor_base->large_B] - DIRECT_ABANDON_MARGIN);
}
#endif

#if defined(DIRECT_ABANDON_HOPELESS)
// largest cofactor finalization keeps under this LP variation
static double direct_cofactor_max(static_conf_t* sconf)
{
    if (sconf->use_dlp >= 3) return sconf->large_prime_max4;
    if (sconf->use_dlp >= 2) return sconf->large_prime_max3;
    if (sconf->use_dlp >= 1) return (double)sconf->large_prime_max2;
    return (double)sconf->large_prime_max;
}
#endif

// move scan survivor report_num of block bnum, side parity, into the pool.
// Called in scan order, so pool->pos stays sorted.
void direct_save_survivor(uint32_t report_num, uint8_t parity, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_pool_t* pool = dconf->dpool;
    fb_list* fb = sconf->factor_base;
    uint32_t k, pos;
    int smooth_num = dconf->smooth_num[report_num];

    if (pool->num == pool->alloc)
        direct_pool_grow(pool);

    k = pool->num++;
    pos = parity * pool->interval + (bnum << 15) + dconf->reports[report_num];
    pool->pos[k] = pos;
    mpz_swap(pool->Q[k], dconf->Qvals[report_num]);
    pool->smooth_num[k] = smooth_num;
    if (smooth_num >= 0)
        memcpy(&pool->fb_offsets[(size_t)k * MAX_SMOOTH_PRIMES],
            dconf->fb_offsets[report_num], (smooth_num + 1) * sizeof(uint32_t));
    pool->bitmap[pos >> 6] |= (1ULL << (pos & 63));

    // all fb primes below large_B are divided out by now: count the
    // survivors the [large_B, B) resieve can still divide
    pool->act[k] = (fb->large_B < fb->B) &&
        direct_q_open(pool->Q[k], fb->list->prime[fb->large_B]);

#if defined(DIRECT_ABANDON_NOLARGE) || defined(DIRECT_ABANDON_HOPELESS)
    pool->ab_saved++;
#endif
#if defined(DIRECT_ABANDON_NOLARGE)
    int large_logp = direct_large_logp(report_num, sconf, dconf);
#if defined(DIRECT_ABANDON_EXHAUSTED)
    // range logp still unaccounted for; the resieve hits count it down
    pool->ab_rem[k] = (int16_t)large_logp;
    pool->logprime = fb->list->logprime;
#endif
#if defined(DIRECT_ABANDON_STATS)
    // calibration: histogrammed at finalization by whether the resieve
    // divided the survivor
    pool->ab_est[k] = (int16_t)large_logp;
    {
        // diagnostic: sieved logp S against the size of Q(x) itself
        mpz_t qx, t;
        uint32_t off = (bnum << sconf->qs_blockbits) + dconf->reports[report_num];
        int v = dconf->sieve[dconf->reports[report_num]], d;
        mpz_init(qx);
        mpz_init(t);
        if (sconf->knmod8 == 1)
        {
            mpz_mul_ui(t, dconf->curr_poly->mpz_poly_a, off);
            if (parity) mpz_sub(qx, t, dconf->curr_poly->mpz_poly_b);
            else mpz_add(qx, t, dconf->curr_poly->mpz_poly_b);
        }
        else
        {
            mpz_mul_ui(t, dconf->curr_poly->mpz_poly_a, off);
            if (parity) mpz_submul_ui(t, dconf->curr_poly->mpz_poly_b, 2);
            else mpz_addmul_ui(t, dconf->curr_poly->mpz_poly_b, 2);
            mpz_set(qx, t);
        }
        mpz_mul_ui(qx, qx, off);
        mpz_add(qx, qx, dconf->curr_poly->mpz_poly_c);
        if (v > (int)sconf->blockinit) v -= 256;
        d = ((int)sconf->blockinit - v) - (int)mpz_sizeinbase(qx, 2);
        d += 128;
        pool->ab_byte[k] = (uint8_t)((d < 0) ? 0 : ((d > 255) ? 255 : d));
        pool->ab_allq[pool->ab_byte[k]]++;
        if (d > 128 + 4)
        {
            uint32_t loc = dconf->reports[report_num];
            pool->ab_junk++;
            if (loc == 32767) pool->ab_junk_last++;
            if (loc == 0) pool->ab_junk_first++;
            pool->ab_junk_blk[(bnum < 15) ? bnum : 15]++;
            pool->ab_junk_side[parity]++;
            pool->ab_junk_loc[loc >> 11]++;
        }
        mpz_clear(qx);
        mpz_clear(t);
    }
    pool->ab_sn0[k] = (int16_t)smooth_num;
#endif
    if (pool->act[k] && direct_no_large_hit(large_logp, sconf))
    {
        pool->act[k] = 2;
        pool->ab_num[0]++;
    }
#endif
#if defined(DIRECT_ABANDON_HOPELESS)
    if ((pool->act[k] == 1) && (mpz_cmp_d(pool->Q[k], direct_cofactor_max(sconf) *
        pow((double)sconf->pmax, DIRECT_ABANDON_K)) > 0))
    {
        pool->act[k] = 3;
        pool->ab_num[1]++;
    }
#endif

#if defined(DIRECT_ABANDON_STATS)
    // abandoned survivors stay counted, so resieving sees every one
    pool->active += (pool->act[k] != 0);
#else
    pool->active += (pool->act[k] == 1);
#endif

    //uint32_t p = sconf->factor_base->list->prime[sconf->factor_base->large_B];
    //if (mpz_cmp_ui(dconf->Qvals[report_num], p * p) >= 0)
    //    pool->active++;

    return;
}

// a root of fb prime fbi hit interval position pos: if a survivor sits
// there, divide it out (guarded: a wrong hit can only waste a tdiv_ui).
void direct_resieve_hit(direct_pool_t* pool, uint32_t pos,
    uint32_t fbi, uint32_t prime)
{
    int lo, hi, mid, k;
    int smooth_num;
    uint32_t* fb_offsets;

    if ((pool->bitmap[pos >> 6] & (1ULL << (pos & 63))) == 0)
        return;

    // binary search the sorted pool for this position
    lo = 0;
    hi = (int)pool->num - 1;
    k = -1;
    while (lo <= hi)
    {
        mid = (lo + hi) >> 1;
        if (pool->pos[mid] == pos) { k = mid; break; }
        if (pool->pos[mid] < pos) lo = mid + 1;
        else hi = mid - 1;
    }
    if (k < 0)
    {
#if defined(DIRECT_ABANDON_STATS)
        pool->ab_miss++;
#endif
        return;
    }

    if ((mpz_sgn(pool->Q[k]) == 0) || (mpz_tdiv_ui(pool->Q[k], prime) != 0))
        return;

    smooth_num = pool->smooth_num[k];
    fb_offsets = &pool->fb_offsets[(size_t)k * MAX_SMOOTH_PRIMES];
    do
    {
        if (smooth_num < (MAX_SMOOTH_PRIMES - 1))
            fb_offsets[++smooth_num] = fbi;
        mpz_tdiv_q_ui(pool->Q[k], pool->Q[k], prime);
    } while (mpz_tdiv_ui(pool->Q[k], prime) == 0);
    pool->smooth_num[k] = smooth_num;

#if defined(DIRECT_ABANDON_STATS)
    // an abandoned survivor that the resieve divides was abandoned wrongly;
    // it stays counted in active (see direct_save_survivor)
    if ((pool->act[k] == 2) || (pool->act[k] == 3))
    {
        pool->ab_wrong[pool->act[k] - 2]++;
        pool->act[k] = 4;     // counted once; stays in active
    }
#endif

#if defined(DIRECT_ABANDON_EXHAUSTED)
    // the range primes divided so far account for all of the range logp
    // in this survivor's sieve byte: no more resieve hits can come (e.g.
    // Q = q * L1 * L2 is left as L1 * L2 > prime^2, never hit again)
    pool->ab_rem[k] = (int16_t)(pool->ab_rem[k] - (int)pool->logprime[fbi]);
    if ((pool->act[k] == 1) &&
        (pool->ab_rem[k] < ((int)pool->logprime[fbi] - DIRECT_ABANDON_MARGIN)))
    {
        pool->act[k] = 2;
        pool->ab_num[0]++;
#if !defined(DIRECT_ABANDON_STATS)
        pool->active--;
#endif
        return;
    }
#endif

    // resieving is done in ascending prime order (to within a 16-group),
    // so once Q < prime^2 it is 1 or prime.  The flag keeps a survivor
    // that is hit again (Q itself a larger fb prime) from being
    // uncounted twice.
    if ((pool->act[k] == 1) && !direct_q_open(pool->Q[k], prime))
    
    {
        pool->act[k] = 0;
        pool->active--;
    }

    //if (mpz_cmp_ui(pool->Q[k], prime * prime) < 0)
    //    pool->active--;

    return;
}

// resieve one 16-group (fb index j) of [large_B, B) against the survivors,
// in the AVX-512 order (pos root1 lanes, pos root2 lanes, neg root1, neg
// root2), so survivors' fb_offsets come out identical on every ISA.
// Branch-free in-range + survivor-bitmap test (in-range is a coin flip, so
// a branch on it mispredicts constantly); only the rare survivor hit
// branches to the out-of-line handler.  Out-of-range lanes test bit 0 of a
// zero word instead.
SIQS_FORCEINLINE void direct_resieve_group(direct_pool_t* pool, uint32_t j,
    const uint32_t* p, const uint32_t* r1, const uint32_t* r2)
{
    const uint64_t* bitmap = pool->bitmap;
    uint32_t interval = pool->interval;
    int k;

#define DIRECT_SCALAR_HIT(r, off, k) \
    do { uint32_t _r = (r); \
        uint32_t _m = 0 - (uint32_t)(_r < interval); \
        uint32_t _q = (_r + (off)) & _m; \
        uint64_t _w = bitmap[_q >> 6] & (uint64_t)(int64_t)(int32_t)_m; \
        if ((_w >> (_q & 63)) & 1) \
            direct_resieve_hit(pool, _q, j + (k), p[k]); } while (0)

    for (k = 0; k < 16; k++) DIRECT_SCALAR_HIT(r1[k], 0, k);
    for (k = 0; k < 16; k++) DIRECT_SCALAR_HIT(r2[k], 0, k);
    for (k = 0; k < 16; k++) DIRECT_SCALAR_HIT(p[k] - r1[k], interval, k);
    for (k = 0; k < 16; k++) DIRECT_SCALAR_HIT(p[k] - r2[k], interval, k);
#undef DIRECT_SCALAR_HIT
}

// resieve the 16-groups of [large_B, B) from the hit list the compress
// kernels left for this poly: the list holds exactly the in-range roots,
// so a group can have a survivor hit only if one of its entries' positions
// is marked in the bitmap.  Such a group is resieved in full with
// direct_resieve_group, so the divisions are the same, in the same order,
// as testing every root.  Must run before the roots are advanced.
void direct_resieve_hits(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    fb_list* fb = sconf->factor_base;
    direct_pool_t* pool = dconf->dpool;
    const uint64_t* bitmap = pool->bitmap;
    const uint32_t* hits = pool->hits;
    const uint32_t* gend = pool->hit_gend;
    uint32_t* primes = dconf->update_data.prime;
    uint32_t* fr1 = (uint32_t*)dconf->update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)dconf->update_data.firstroots2;
    uint32_t ngroups = (fb->B - fb->large_B) / 16;
    uint32_t g, k = 0;

    // straight bitmap tests: a register-sized prefilter in front of them
    // passes about half the groups here (20-30 entries each), which made
    // it slower
    for (g = 0; g < ngroups; g++)
    {
        uint32_t end = gend[g];
        uint64_t any = 0;

        for (; k < end; k++)
        {
            uint32_t pos = hits[k] & 0xffffff;
            any |= bitmap[pos >> 6] >> (pos & 63);
        }

        if (any & 1)
        {
            uint32_t j = fb->large_B + 16 * g;
            direct_resieve_group(pool, j, primes + j, fr1 + j, fr2 + j);
            if (pool->active == 0)
                break;
        }
    }
}

// for vectorized walks of the hit list: if hit list entry k sits on a
// survivor, resieve its 16-group in full (as direct_resieve_hits) and
// return the index one past the group's entries, else return k + 1.
// Entries must be visited in ascending order, skipping to the returned
// index, so that each group is resieved at most once.
uint32_t direct_resieve_hit_entry(static_conf_t* sconf, dynamic_conf_t* dconf, uint32_t k)
{
    fb_list* fb = sconf->factor_base;
    direct_pool_t* pool = dconf->dpool;
    const uint32_t* gend = pool->hit_gend;
    uint32_t pos = pool->hits[k] & 0xffffff;
    uint32_t lo, hi, j;

    if (((pool->bitmap[pos >> 6] >> (pos & 63)) & 1) == 0)
        return k + 1;

    // the group holding entry k: the first with gend > k
    lo = 0;
    hi = (fb->B - fb->large_B) / 16 - 1;
    while (lo < hi)
    {
        uint32_t mid = (lo + hi) >> 1;
        if (gend[mid] > k) hi = mid;
        else lo = mid + 1;
    }

    j = fb->large_B + 16 * lo;
    direct_resieve_group(pool, j, dconf->update_data.prime + j,
        (uint32_t*)dconf->update_data.firstroots1 + j,
        (uint32_t*)dconf->update_data.firstroots2 + j);
    return gend[lo];
}

// finalize every survivor of this poly (both sides): copy it into report
// slot 0 of the scan arrays and run the normal finalization.  Must run
// before nextB (trial_divide_Q_siqs uses curr_poly).  Leaves the pool
// empty and its bitmap zeroed.
void direct_finalize_survivors(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_pool_t* pool = dconf->dpool;
    uint32_t interval = pool->interval;
    uint32_t k;

#if defined(DIRECT_EXIT_STATS)
    for (k = 0; k < pool->num; k++)
    {
        // the survivors still active at the end: what are they?
        if (pool->act[k] == 1)
        {
            size_t b = mpz_sizeinbase(pool->Q[k], 2);
            pool->ex_qbits[(b < 255) ? b : 255]++;
            pool->ex_loc0 += ((pool->pos[k] & 32767) == 0);
        }
    }
#endif

    for (k = 0; k < pool->num; k++)
    {
        uint32_t pos = pool->pos[k];
        uint8_t side = (pos >= interval) ? 1 : 0;
        uint32_t loc = pos - side * interval;
        int smooth_num = pool->smooth_num[k];

        pool->bitmap[pos >> 6] &= ~(1ULL << (pos & 63));

#if defined(DIRECT_ABANDON_STATS) && defined(DIRECT_ABANDON_NOLARGE)
        {
            // estimate error: est less the exact logp of the range primes
            // the resieve divided (once per prime, as sieved)
            uint32_t* offs = &pool->fb_offsets[(size_t)k * MAX_SMOOTH_PRIMES];
            uint32_t prev = 0xffffffff;
            int j, e = pool->ab_est[k];
            for (j = pool->ab_sn0[k] + 1; j <= smooth_num; j++)
            {
                if (offs[j] != prev)
                    e -= (int)sconf->factor_base->list->logprime[offs[j]];
                prev = offs[j];
            }
            // sample the big-error no-hit survivors: does any range prime
            // still divide Q (a resieve miss)?
            if ((e >= 31) && (smooth_num == pool->ab_sn0[k]) &&
                ((pool->ab_hist[0][63] % 10) == 0))
            {
                uint32_t pi;
                pool->ab_samp++;
                pool->ab_bytes[pool->ab_byte[k]]++;
                for (pi = sconf->factor_base->large_B; pi < sconf->factor_base->B; pi++)
                {
                    if (mpz_tdiv_ui(pool->Q[k], sconf->factor_base->list->prime[pi]) == 0)
                    {
                        pool->ab_sampdiv++;
                        break;
                    }
                }
            }
            e += 32;
            e = (e < 0) ? 0 : ((e > 63) ? 63 : e);
            pool->ab_hist[smooth_num > pool->ab_sn0[k]][e]++;
        }
#endif

        mpz_swap(dconf->Qvals[0], pool->Q[k]);
        dconf->smooth_num[0] = smooth_num;
        if (smooth_num >= 0)
            memcpy(dconf->fb_offsets[0], &pool->fb_offsets[(size_t)k * MAX_SMOOTH_PRIMES],
                (smooth_num + 1) * sizeof(uint32_t));
        dconf->reports[0] = loc & 32767;
        dconf->valid_Qs[0] = 1;

        trial_divide_Q_siqs(0, side, dconf->numB - 1, loc >> 15, sconf, dconf);
    }

#if defined(DIRECT_EXIT_STATS)
    // polys whose survivors all resolved (so the resieve could stop early)
    pool->ex_polys++;
    pool->ex_resolved += ((pool->num > 0) && (pool->active == 0));
    pool->ex_hist[(pool->active < 15) ? pool->active : 15]++;
#endif
    pool->num = 0;
    pool->active = 0;

    return;
}

// subtract 4 from each byte of a 32 KB block (the minblock init, applied
// after the fact: the sieve is subtractive, so order doesn't matter).
// Explicit SIMD: gcc -O2 does not vectorize the byte loop (the scalar loop
// cost ~4% of total time).
void direct_block_sub4(uint8_t* blk)
{
    int k;
#if defined(USE_AVX512BW)
    __m512i v4 = _mm512_set1_epi8(4);
    for (k = 0; k < 32768; k += 64)
        _mm512_store_si512((__m512i*)(blk + k),
            _mm512_sub_epi8(_mm512_load_si512((__m512i*)(blk + k)), v4));
#elif defined(USE_AVX2)
    __m256i v4 = _mm256_set1_epi8(4);
    for (k = 0; k < 32768; k += 32)
        _mm256_store_si256((__m256i*)(blk + k),
            _mm256_sub_epi8(_mm256_load_si256((__m256i*)(blk + k)), v4));
#elif defined(USE_SSE41) || (defined(__SSE2__) && !defined(FORCE_GENERIC))
    __m128i v4 = _mm_set1_epi8(4);
    for (k = 0; k < 32768; k += 16)
        _mm_store_si128((__m128i*)(blk + k),
            _mm_sub_epi8(_mm_load_si128((__m128i*)(blk + k)), v4));
#else
    // SWAR per-byte subtract: with each byte's top bit forced on, no
    // borrow crosses a byte boundary; then fix each top bit
    uint64_t* w = (uint64_t*)blk;
    const uint64_t hi = 0x8080808080808080ULL;
    const uint64_t four = 0x0404040404040404ULL;
    for (k = 0; k < 32768 / 8; k++)
    {
        uint64_t x = w[k];
        w[k] = ((x | hi) - four) ^ ((x ^ ~four) & hi);
    }
#endif
}

// ----------------------------------------------------------------------
// scalar kernels
// ----------------------------------------------------------------------

// scalar counterpart of direct_scatter16 (direct_sieve_avx512.c): one
// 16-prime group of [large_B, B) into both intervals.  p, r1, r2 are the
// group's primes and roots (pos-side roots; neg side is p - r).  CREG /
// CMEM append to the hit list hits[*nc...] and sieve the entries past
// *flushed once DIRECT_CBUF_FLUSH of them are pending; the caller flushes
// the rest after its loop.
SIQS_FORCEINLINE void direct_scatter16_scalar(uint8_t* sieve, uint32_t interval,
    const uint32_t* p, const uint32_t* r1, const uint32_t* r2, uint8_t logp,
    const int kmode, uint32_t* hits, uint32_t* nc, uint32_t* flushed, uint32_t* sink)
{
    uint8_t* sieven = sieve + interval;
    int k;

    if ((kmode == SIQS_K_CREG) || (kmode == SIQS_K_CMEM))
    {
        // branch-free append, in the AVX-512 compress order
        uint32_t n = *nc;
        uint32_t tag = (uint32_t)logp << 24;
        for (k = 0; k < 16; k++) { uint32_t r = r1[k]; hits[n] = r | tag; n += (r < interval); }
        for (k = 0; k < 16; k++) { uint32_t r = r2[k]; hits[n] = r | tag; n += (r < interval); }
        for (k = 0; k < 16; k++) { uint32_t r = p[k] - r1[k]; hits[n] = (r + interval) | tag; n += (r < interval); }
        for (k = 0; k < 16; k++) { uint32_t r = p[k] - r2[k]; hits[n] = (r + interval) | tag; n += (r < interval); }
        if ((n - *flushed) >= DIRECT_CBUF_FLUSH)
        {
            direct_compress_flush(sieve, hits + *flushed, n - *flushed);
            *flushed = n;
        }
        *nc = n;
    }
    else if (kmode == SIQS_K_SINK)
    {
        // misses go to a rotating stripe in the second tail guard block
        uint32_t s = 2 * interval + 32768 + *sink;
        for (k = 0; k < 16; k++) { uint32_t r = r1[k]; sieve[(r < interval) ? r : s + k] -= logp; }
        for (k = 0; k < 16; k++) { uint32_t r = r2[k]; sieve[(r < interval) ? r : s + 16 + k] -= logp; }
        for (k = 0; k < 16; k++) { uint32_t r = p[k] - r1[k]; sieve[(r < interval) ? interval + r : s + 32 + k] -= logp; }
        for (k = 0; k < 16; k++) { uint32_t r = p[k] - r2[k]; sieve[(r < interval) ? interval + r : s + 48 + k] -= logp; }
        *sink += 64;
        if (*sink >= DIRECT_SINK_BYTES)
            *sink = 0;
    }
    else
    {
        for (k = 0; k < 16; k++) if (r1[k] < interval) sieve[r1[k]] -= logp;
        for (k = 0; k < 16; k++) if (r2[k] < interval) sieve[r2[k]] -= logp;
        for (k = 0; k < 16; k++) if ((p[k] - r1[k]) < interval) sieven[p[k] - r1[k]] -= logp;
        for (k = 0; k < 16; k++) if ((p[k] - r2[k]) < interval) sieven[p[k] - r2[k]] -= logp;
    }
}

// direct pass 1 (see direct_sieve_largeprimes_body in direct_sieve_avx512.c)
SIQS_FORCEINLINE void direct_sieve_largeprimes_scalar_body(uint8_t* sieve,
    static_conf_t* sconf, dynamic_conf_t* dconf, const int kmode)
{
    fb_list* fb = sconf->factor_base;
    update_t update_data = dconf->update_data;
    uint32_t interval = sconf->num_blocks << 15;
    uint8_t* sieven = sieve + interval;
    uint32_t* primes = update_data.prime;
    uint32_t* fr1 = (uint32_t*)update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)update_data.firstroots2;
    uint32_t j;
    uint8_t logp;
    direct_pool_t* pool = dconf->dpool;
    uint32_t* hits;
    uint32_t nc = 0, flushed = 0, sink = 0;

    direct_hits_reserve(pool, sconf);
    hits = pool->hits;

    for (j = fb->large_B; (j + 16) <= fb->B; j += 16)
    {
        // approximate the next 16 logp's as equal to the largest
        logp = update_data.logp[j + 15] + DIRECT_LOGP_BIAS;
        direct_scatter16_scalar(sieve, interval, primes + j, fr1 + j, fr2 + j,
            logp, kmode, hits, &nc, &flushed, &sink);
        pool->hit_gend[(j - fb->large_B) / 16] = nc;
    }
    direct_compress_flush(sieve, hits + flushed, nc - flushed);

    // scalar tail (empty when B - large_B is a multiple of 16)
    for (; j < fb->B; j++)
    {
        uint32_t prime = primes[j];
        logp = update_data.logp[j] + DIRECT_LOGP_BIAS;

        if (fr1[j] < interval) sieve[fr1[j]] -= logp;
        if (fr2[j] < interval) sieve[fr2[j]] -= logp;
        if ((prime - fr1[j]) < interval) sieven[prime - fr1[j]] -= logp;
        if ((prime - fr2[j]) < interval) sieven[prime - fr2[j]] -= logp;
    }
}

void direct_sieve_largeprimes_scalar_bitloop(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_scalar_body(sieve, sconf, dconf, SIQS_K_BITLOOP);
}

void direct_sieve_largeprimes_scalar_sink(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_scalar_body(sieve, sconf, dconf, SIQS_K_SINK);
}

void direct_sieve_largeprimes_scalar_append(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    direct_sieve_largeprimes_scalar_body(sieve, sconf, dconf, SIQS_K_CMEM);
}

// fused resieve / root update / next-poly scatter (see
// direct_resieve_update_body in direct_sieve_avx512.c).  Hits are tested in
// the AVX-512 order (per 16-group: pos root1 lanes, pos root2 lanes, neg
// root1, neg root2), so survivors' fb_offsets come out identical.  The
// compress kernels resieve from the hit list their previous scatter left
// (direct_resieve_hits) and leave the next poly's.
SIQS_FORCEINLINE void direct_resieve_update_scalar_body(static_conf_t* sconf,
    dynamic_conf_t* dconf, uint8_t* sieve, const int kmode)
{
    fb_list* fb = sconf->factor_base;
    direct_pool_t* pool = dconf->dpool;
    update_t update_data = dconf->update_data;
    uint32_t interval = pool->interval;
    uint32_t large_B = fb->large_B;
    uint32_t bound = fb->B;
    uint32_t* primes = update_data.prime;
    uint32_t* fr1 = (uint32_t*)update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)update_data.firstroots2;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    uint32_t* ptr = (uint32_t*)&dconf->rootupdates[(v - 1) * bound + large_B];
    int resieve = (pool->active > 0);
    uint8_t* sieven = (sieve == NULL) ? NULL : sieve + interval;
    const int hitlist = ((kmode == SIQS_K_CREG) || (kmode == SIQS_K_CMEM));
    uint32_t j;
    int k;
    uint8_t logp;
    uint32_t* hits;
    uint32_t nc = 0, flushed = 0, sink = 0;

    if ((sieve == NULL) && (!resieve))
        return;

    direct_hits_reserve(pool, sconf);
    hits = pool->hits;

    if (resieve && hitlist)
    {
        direct_resieve_hits(sconf, dconf);
        resieve = (pool->active > 0);
    }

    for (j = large_B; (j + 16) <= bound; j += 16, ptr += 16)
    {
        uint32_t* p = primes + j;
        uint32_t* r1 = fr1 + j;
        uint32_t* r2 = fr2 + j;

        if (resieve && !hitlist)
        {
            direct_resieve_group(pool, j, p, r1, r2);
            if (pool->active == 0)
                resieve = 0;
        }

        if (sieve == NULL)
            continue;

        // advance to the next poly; the pair stays unordered.  Branch-free:
        // roots and primes are < 2^31, so the sign bit flags the wrap.
        if (sign > 0)
        {
            for (k = 0; k < 16; k++)
            {
                uint32_t u = ptr[k];
                uint32_t a = r1[k] - u, b = r2[k] - u;
                r1[k] = a + (p[k] & (0 - (a >> 31)));
                r2[k] = b + (p[k] & (0 - (b >> 31)));
            }
        }
        else
        {
            for (k = 0; k < 16; k++)
            {
                uint32_t u = ptr[k];
                uint32_t a = r1[k] + u - p[k], b = r2[k] + u - p[k];
                r1[k] = a + (p[k] & (0 - (a >> 31)));
                r2[k] = b + (p[k] & (0 - (b >> 31)));
            }
        }

        // sieve the next poly's roots, both sides (as pass 1)
        logp = update_data.logp[j + 15] + DIRECT_LOGP_BIAS;
        direct_scatter16_scalar(sieve, interval, p, r1, r2, logp, kmode,
            hits, &nc, &flushed, &sink);
        pool->hit_gend[(j - large_B) / 16] = nc;
    }

    // nc stays 0 when sieve == NULL
    direct_compress_flush(sieve, hits + flushed, nc - flushed);

    // scalar tail (empty when B - large_B is a multiple of 16)
    for (; j < bound; j++, ptr++)
    {
        uint32_t prime = primes[j];
        uint32_t root1 = fr1[j];
        uint32_t root2 = fr2[j];
        uint32_t u = *ptr;

        if (resieve)
        {
            if (root1 < interval) direct_resieve_hit(pool, root1, j, prime);
            if (root2 < interval) direct_resieve_hit(pool, root2, j, prime);
            if ((prime - root1) < interval)
                direct_resieve_hit(pool, interval + prime - root1, j, prime);
            if ((prime - root2) < interval)
                direct_resieve_hit(pool, interval + prime - root2, j, prime);
        }

        if (sieve == NULL)
            continue;

        if (sign > 0)
        {
            root1 = (u > root1) ? root1 + prime - u : root1 - u;
            root2 = (u > root2) ? root2 + prime - u : root2 - u;
        }
        else
        {
            root1 += u;
            root2 += u;
            if (root1 >= prime) root1 -= prime;
            if (root2 >= prime) root2 -= prime;
        }
        fr1[j] = root1;
        fr2[j] = root2;

        logp = update_data.logp[j] + DIRECT_LOGP_BIAS;
        if (root1 < interval) sieve[root1] -= logp;
        if (root2 < interval) sieve[root2] -= logp;
        if ((prime - root1) < interval) sieven[prime - root1] -= logp;
        if ((prime - root2) < interval) sieven[prime - root2] -= logp;
    }
}

void direct_resieve_update_scalar_bitloop(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    direct_resieve_update_scalar_body(sconf, dconf, sieve, SIQS_K_BITLOOP);
}

void direct_resieve_update_scalar_sink(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    direct_resieve_update_scalar_body(sconf, dconf, sieve, SIQS_K_SINK);
}

void direct_resieve_update_scalar_append(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    direct_resieve_update_scalar_body(sconf, dconf, sieve, SIQS_K_CMEM);
}

// sparse sieve of one 32k block for fb [med_B, fb_blk8i) (see
// large_sieveblock_body in direct_sieve_avx512.c): every prime is > 32768
// and root1 <= root2, so each root hits the block at most once.
SIQS_FORCEINLINE void large_sieveblock_scalar_body(uint8_t* sieve, sieve_fb* fb,
    fb_list* full_fb, const int kmode)
{
    uint32_t i;
    int k;
    uint32_t prime, root1, root2;
    uint8_t logp;
    uint32_t bound = full_fb->fb_blk8i;
    ALIGNED_MEM uint32_t cbuf[DIRECT_CBUF_SIZE];
    uint32_t nc = 0;

    for (i = full_fb->med_B; (i + 16) <= bound; i += 16)
    {
        // approximate the next 16 logp's as equal to the largest
        logp = fb->logp[i + 15];

        if (kmode == SIQS_K_BITLOOP)
        {
            for (k = 0; k < 16; k++)
            {
                root1 = fb->root1[i + k];
                root2 = fb->root2[i + k];
                if (root1 < 32768) { sieve[root1] -= logp; root1 += fb->prime[i + k]; }
                if (root2 < 32768) { sieve[root2] -= logp; root2 += fb->prime[i + k]; }
                root1 -= 32768;
                root2 -= 32768;
                fb->root1[i + k] = (root1 < root2) ? root1 : root2;
                fb->root2[i + k] = (root1 < root2) ? root2 : root1;
            }
        }
        else
        {
            uint32_t tag = (uint32_t)logp << 24;
            for (k = 0; k < 16; k++)
            {
                prime = fb->prime[i + k];
                root1 = fb->root1[i + k];
                root2 = fb->root2[i + k];
                cbuf[nc] = root1 | tag;
                nc += (root1 < 32768);
                cbuf[nc] = root2 | tag;
                nc += (root2 < 32768);
                root1 += (root1 < 32768) ? prime : 0;
                root2 += (root2 < 32768) ? prime : 0;
                root1 -= 32768;
                root2 -= 32768;
                fb->root1[i + k] = (root1 < root2) ? root1 : root2;
                fb->root2[i + k] = (root1 < root2) ? root2 : root1;
            }
            if (nc >= DIRECT_CBUF_FLUSH)
            {
                direct_compress_flush(sieve, cbuf, nc);
                nc = 0;
            }
        }
    }
    direct_compress_flush(sieve, cbuf, nc);

    // primes left over after the 16-groups (none when med_B and fb_blk8i
    // are both 16-aligned); same as the AVX-512 tail
    for (; i < bound; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        while (root2 < (uint32_t)(32768))
        {
            sieve[root1] -= logp;
            sieve[root2] -= logp;
            root1 += prime;
            root2 += prime;
        }

        if (root1 < (uint32_t)(32768))
        {
            sieve[root1] -= logp;
            root1 += prime;
        }

        root1 -= 32768;
        root2 -= 32768;
        fb->root1[i] = (root1 < root2) ? root1 : root2;
        fb->root2[i] = (root1 < root2) ? root2 : root1;
    }
}

void large_sieveblock_scalar_bitloop(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    large_sieveblock_scalar_body(sieve, fb, full_fb, SIQS_K_BITLOOP);
}

void large_sieveblock_scalar_append(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    large_sieveblock_scalar_body(sieve, fb, full_fb, SIQS_K_CMEM);
}

// the sparse-sieve tdiv for [med_B, fb_blk8i) (see tdiv_large_body in
// direct_sieve_avx512.c): guarded divisions for the 16-groups, in ascending
// fb index order, then the AVX-512 version's scalar tail.
void tdiv_large_scalar(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    uint32_t i;
    int k;
    uint32_t tmp, prime, root1, root2, report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb* fbc = parity ? dconf->sieve_n : dconf->sieve_p;
    uint32_t block_loc;
    uint32_t bound = sconf->factor_base->fb_blk8i;

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {
        if (!dconf->valid_Qs[report_num])
            continue;

        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        for (i = sconf->factor_base->med_B; (i + 16) <= bound; i += 16)
        {
            for (k = 0; k < 16; k++)
            {
                tmp = block_loc + fbc->prime[i + k] - 32768;
                if ((tmp == fbc->root1[i + k]) || (tmp == fbc->root2[i + k]))
                {
                    uint32_t j = i + k;
                    DIVIDE_RESIEVED_PRIME_2(j);
                }
            }
        }

        for (; i < bound; i++)
        {
            prime = fbc->prime[i];
            root1 = fbc->root1[i];
            root2 = fbc->root2[i];

            // these primes are all > blocksize, so we only need one step
            tmp = block_loc + prime - 32768;
            if (tmp == root1 || tmp == root2)
            {
                DIVIDE_ONE_PRIME;
            }
        }

        dconf->smooth_num[report_num] = smooth_num;
    }

    return;
}

// ----------------------------------------------------------------------
// portable root updates for [med_B, B)
// ----------------------------------------------------------------------
// Runs after a small-prime update (nextRoots_32k_generic_small or an ISA
// version of it), in place of the ISA-specific medbucket/largebucket
// passes:
//  [med_B, fb_blk8i)    sparse sieve: roots to sieve_p/n in min/max order
//  [fb_blk8i, large_B)  mid buckets (when fb_blk8i == med_B)
//  [large_B, B)         large buckets, unless the direct sieve owns it
// The bucket fills use the same element format and slice bookkeeping as
// firstRoots_32k, which every ISA's lp_sieveblock / tdiv_LP reads.
// nextRoots_32k_lp_portable = nextRoots_32k_lp_sparse_scalar (the sparse
// range) + nextRoots_32k_lp_buckets (the bucket ranges); the AVX2 version
// swaps in a vector sparse update.

// [med_B, fb_blk8i): sparse-sieve roots to sieve_p/n in min/max order
void nextRoots_32k_lp_sparse_scalar(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    fb_list* fb = sconf->factor_base;
    update_t update_data = dconf->update_data;
    uint32_t bound = fb->B;
    uint32_t med_B = fb->med_B;
    uint32_t sparse_B = fb->fb_blk8i;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    uint32_t* ptr = (uint32_t*)&dconf->rootupdates[(v - 1) * bound + med_B];
    uint32_t j, root1, root2, prime;

    for (j = med_B; j < sparse_B; j++, ptr++)
    {
        uint32_t u = *ptr;
        prime = update_data.prime[j];
        root1 = update_data.firstroots1[j];
        root2 = update_data.firstroots2[j];

        if (sign > 0)
        {
            root1 = (u > root1) ? root1 + prime - u : root1 - u;
            root2 = (u > root2) ? root2 + prime - u : root2 - u;
        }
        else
        {
            root1 += u;
            root2 += u;
            if (root1 >= prime) root1 -= prime;
            if (root2 >= prime) root2 -= prime;
        }

        if (root2 < root1) { uint32_t t = root1; root1 = root2; root2 = t; }
        update_data.firstroots1[j] = root1;
        update_data.firstroots2[j] = root2;
        dconf->sieve_p->root1[j] = root1;
        dconf->sieve_p->root2[j] = root2;
        dconf->sieve_n->root1[j] = prime - root2;
        dconf->sieve_n->root2[j] = prime - root1;
    }
}

// bucket ranges [fb_blk8i, large_B) (mid buckets) and [large_B, B) (large
// buckets, unless the direct sieve owns that range)
void nextRoots_32k_lp_buckets(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    fb_list* fb = sconf->factor_base;
    update_t update_data = dconf->update_data;
    lp_bucket* lp_bucket_p = dconf->buckets;
    uint32_t bound = fb->B;
    uint32_t sparse_B = fb->fb_blk8i;
    uint32_t large_B = fb->large_B;
    uint32_t stop = sconf->lp_direct ? large_B : bound;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    uint32_t* ptr = (uint32_t*)&dconf->rootupdates[(v - 1) * bound + sparse_B];
    uint32_t j, root1, root2, prime, nroot1, nroot2, interval;
    int numblocks = sconf->num_blocks;
    uint32_t* bptr;
    int bnum, room, k;
    int bound_index = 0;
    int check_bound;
    uint32_t bound_val = sparse_B;
    uint32_t* numptr_p, * numptr_n, * sliceptr_p, * sliceptr_n;
    uint8_t* slicelogp_ptr;
    uint32_t* slicebound_ptr;
    uint8_t logp = 0;

    interval = numblocks << 15;

    if ((lp_bucket_p->alloc_slices == 0) || (sparse_B >= stop))
    {
        if (lp_bucket_p->alloc_slices != 0)
            lp_bucket_p->num_slices = 0;
        return;
    }

    // bucket ranges [sparse_B, stop)
    lp_bucket_p->fb_bounds[0] = sparse_B;
    sliceptr_p = lp_bucket_p->list;
    sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);
    numptr_p = lp_bucket_p->num;
    numptr_n = lp_bucket_p->num + numblocks;
    for (j = 0; j < lp_bucket_p->list_size; j++)
        numptr_p[j] = 0;
    lp_bucket_p->num_slices = 0;
    slicelogp_ptr = lp_bucket_p->logp;
    slicebound_ptr = lp_bucket_p->fb_bounds;
    check_bound = sparse_B + BUCKET_ALLOC / 2;

    logp = update_data.logp[sparse_B - 1];
    for (j = sparse_B; j < stop; j++, ptr++)
    {
        uint32_t u = *ptr;

        CHECK_NEW_SLICE(j);

        prime = update_data.prime[j];
        root1 = update_data.firstroots1[j];
        root2 = update_data.firstroots2[j];

        if (sign > 0)
        {
            root1 = (u > root1) ? root1 + prime - u : root1 - u;
            root2 = (u > root2) ? root2 + prime - u : root2 - u;
        }
        else
        {
            root1 += u;
            root2 += u;
            if (root1 >= prime) root1 -= prime;
            if (root2 >= prime) root2 -= prime;
        }

        update_data.firstroots1[j] = root1;
        update_data.firstroots2[j] = root2;
        nroot1 = prime - root1;
        nroot2 = prime - root2;

        if (j < large_B)
        {
            // p < interval: possibly several hits per root
            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);
        }
        else
        {
            FILL_ONE_PRIME_P(j);
            root1 = nroot1;
            root2 = nroot2;
            FILL_ONE_PRIME_N(j);
        }
    }

    lp_bucket_p->num_slices = bound_index + 1;
    lp_bucket_p->logp[bound_index] = sconf->lp_direct ?
        fb->list->logprime[large_B - 1] : logp;

    return;
}

void nextRoots_32k_lp_portable(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    nextRoots_32k_lp_sparse_scalar(sconf, dconf);
    nextRoots_32k_lp_buckets(sconf, dconf);
}

#endif
