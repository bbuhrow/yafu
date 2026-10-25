/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.
----------------------------------------------------------------------*/

// 4-lane kernels for the large-prime strategies, written once against a
// small vector interface and included by direct_sieve_sse41.c (SSE4.1) and
// direct_sieve_neon.c (NEON).  The including file defines:
//
//   v4t                    4 x uint32 vector type
//   V4_LOAD(p) V4_STORE(p, v)          unaligned load / store of 4 uint32
//   V4_SET1(x)  V4_SETR(a, b, c, d)
//   V4_ADD V4_SUB V4_AND V4_OR         lane-wise
//   V4_ANDNOT(m, x)                    ~m & x
//   V4_LT(a, b)  V4_EQ(a, b)           all-ones lanes where a < b / a == b
//                                      (all values here are < 2^31)
//   V4_MASK(m)                         4-bit mask of all-ones lanes
//   V4_BLEND(a, b, m)                  m ? b : a
//   V4_MIN V4_MAX                      unsigned
//   V4_COMPRESS(v, m)                  set lanes of v (mask m) packed low,
//                                      using v4_cmp_lut[m] (16-byte table)
//   V4_FN(name)                        name##_<isa>
//
// Each 16-prime group is handled as four quarters, so the logp
// approximation, fb boundaries and the order in which hits and divisions
// are visited match the AVX-512 / AVX2 / scalar kernels: all do
// bit-identical work.  CREG and CMEM are the same here (compress through
// the table, then a full store).

// byte shuffle per 4-bit lane mask: the set lanes' bytes packed to the low
// end, the rest 0x80 (zero for both pshufb and vqtbl1q)
static uint8_t V4_FN(v4_cmp_lut)[16][16];
static int V4_FN(v4_lut_ready) = 0;

void V4_FN(direct_v4_init)(void)
{
    uint32_t m, b, n, k;

    if (V4_FN(v4_lut_ready))
        return;

    for (m = 0; m < 16; m++)
    {
        for (k = 0; k < 16; k++)
            V4_FN(v4_cmp_lut)[m][k] = 0x80;
        for (b = 0, n = 0; b < 4; b++)
        {
            if ((m >> b) & 1)
            {
                for (k = 0; k < 4; k++)
                    V4_FN(v4_cmp_lut)[m][4 * n + k] = (uint8_t)(4 * b + k);
                n++;
            }
        }
    }
    V4_FN(v4_lut_ready) = 1;
}

// append the lanes of v selected by the 4-bit mask m to buf[n...] (buf
// needs 4 entries of room past n); returns the new count
SIQS_FORCEINLINE uint32_t V4_FN(v4_cstore)(uint32_t* buf, uint32_t n, uint32_t m, v4t v)
{
    V4_STORE(buf + n, V4_COMPRESS(v, m));
    // popcount of the 4-bit mask from a packed nibble table
    return n + (uint32_t)((0x4332322132212110ull >> (4 * m)) & 15);
}

// ----------------------------------------------------------------------
// direct scatter: one 16-prime group of [large_B, B) into both intervals.
// CREG / CMEM append to the hit list hits[*nc...] and sieve the entries
// past *flushed once DIRECT_CBUF_FLUSH of them are pending.
// ----------------------------------------------------------------------
SIQS_FORCEINLINE void V4_FN(direct_scatter16)(uint8_t* sieve, uint32_t interval,
    const uint32_t* p, const uint32_t* r1, const uint32_t* r2, uint8_t logp,
    const int kmode, uint32_t* hits, uint32_t* nc, uint32_t* flushed, uint32_t* sink)
{
    v4t vint = V4_SET1(interval);
    v4t vtag = V4_SET1((uint32_t)logp << 24);
    v4t pos[16], msk[16];
    int q, k;

    // order: pos root1 quarters, pos root2, neg root1, neg root2
    for (q = 0; q < 4; q++)
    {
        v4t vp = V4_LOAD(p + 4 * q);
        v4t a = V4_LOAD(r1 + 4 * q);
        v4t b = V4_LOAD(r2 + 4 * q);
        v4t na = V4_SUB(vp, a);
        v4t nb = V4_SUB(vp, b);
        pos[q] = a;                      msk[q] = V4_LT(a, vint);
        pos[4 + q] = b;                  msk[4 + q] = V4_LT(b, vint);
        pos[8 + q] = V4_ADD(na, vint);   msk[8 + q] = V4_LT(na, vint);
        pos[12 + q] = V4_ADD(nb, vint);  msk[12 + q] = V4_LT(nb, vint);
    }

    if ((kmode == SIQS_K_CREG) || (kmode == SIQS_K_CMEM))
    {
        uint32_t n = *nc;
        for (k = 0; k < 16; k++)
            n = V4_FN(v4_cstore)(hits, n, V4_MASK(msk[k]), V4_OR(pos[k], vtag));
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
        ALIGNED_MEM uint32_t idx[64];
        v4t vs = V4_ADD(V4_SET1(2 * interval + 32768 + *sink), V4_SETR(0, 1, 2, 3));
        v4t v4 = V4_SET1(4);

        for (k = 0; k < 16; k++)
        {
            V4_STORE(idx + 4 * k, V4_BLEND(vs, pos[k], msk[k]));
            vs = V4_ADD(vs, v4);
        }
        for (k = 0; k < 64; k += 8)
        {
            sieve[idx[k + 0]] -= logp; sieve[idx[k + 1]] -= logp;
            sieve[idx[k + 2]] -= logp; sieve[idx[k + 3]] -= logp;
            sieve[idx[k + 4]] -= logp; sieve[idx[k + 5]] -= logp;
            sieve[idx[k + 6]] -= logp; sieve[idx[k + 7]] -= logp;
        }
        *sink += 64;
        if (*sink >= DIRECT_SINK_BYTES)
            *sink = 0;
    }
    else
    {
        ALIGNED_MEM uint32_t t[4];
        for (k = 0; k < 16; k++)
        {
            uint32_t m = V4_MASK(msk[k]);
            if (m)
            {
                V4_STORE(t, pos[k]);
                while (m) { sieve[t[siqs_ctz32(m)]] -= logp; m &= m - 1; }
            }
        }
    }
}

// direct pass 1
SIQS_FORCEINLINE void V4_FN(direct_sieve_largeprimes_body)(uint8_t* sieve,
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
        V4_FN(direct_scatter16)(sieve, interval, primes + j, fr1 + j, fr2 + j,
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

void V4_FN(direct_sieve_largeprimes_bitloop)(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    V4_FN(direct_sieve_largeprimes_body)(sieve, sconf, dconf, SIQS_K_BITLOOP);
}

void V4_FN(direct_sieve_largeprimes_sink)(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    V4_FN(direct_sieve_largeprimes_body)(sieve, sconf, dconf, SIQS_K_SINK);
}

void V4_FN(direct_sieve_largeprimes_creg)(uint8_t* sieve, static_conf_t* sconf, dynamic_conf_t* dconf)
{
    V4_FN(direct_sieve_largeprimes_body)(sieve, sconf, dconf, SIQS_K_CREG);
}

// ----------------------------------------------------------------------
// fused resieve / root update / next-poly scatter
// ----------------------------------------------------------------------
SIQS_FORCEINLINE void V4_FN(direct_resieve_update_body)(static_conf_t* sconf,
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
    const uint64_t* bitmap = pool->bitmap;
    // CREG resieves the 16-groups from the hit list its previous scatter
    // left (direct_resieve_hits); the other kernels test every root
    const int hitlist = (kmode == SIQS_K_CREG) || (kmode == SIQS_K_CMEM);
    uint32_t j;
    int q, k;
    uint8_t logp;
    uint32_t* hits;
    ALIGNED_MEM uint32_t qpos[64];
    ALIGNED_MEM uint32_t qin[64];
    uint32_t nc = 0, flushed = 0, sink = 0;
    v4t vint = V4_SET1(interval);

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
        if (resieve && !hitlist)
        {
            // positions (masked to 0 when out of range) and in-range masks
            // for the 4 root kinds, lanes 0..15 each: the AVX-512 order
            for (q = 0; q < 4; q++)
            {
                v4t vp = V4_LOAD(primes + j + 4 * q);
                v4t a = V4_LOAD(fr1 + j + 4 * q);
                v4t b = V4_LOAD(fr2 + j + 4 * q);
                v4t na = V4_SUB(vp, a);
                v4t nb = V4_SUB(vp, b);
                v4t m;
                m = V4_LT(a, vint);
                V4_STORE(qin + 4 * q, m);       V4_STORE(qpos + 4 * q, V4_AND(a, m));
                m = V4_LT(b, vint);
                V4_STORE(qin + 16 + 4 * q, m);  V4_STORE(qpos + 16 + 4 * q, V4_AND(b, m));
                m = V4_LT(na, vint);
                V4_STORE(qin + 32 + 4 * q, m);  V4_STORE(qpos + 32 + 4 * q, V4_AND(V4_ADD(na, vint), m));
                m = V4_LT(nb, vint);
                V4_STORE(qin + 48 + 4 * q, m);  V4_STORE(qpos + 48 + 4 * q, V4_AND(V4_ADD(nb, vint), m));
            }

            // branch-free survivor bitmap test; only the rare hit branches
            for (k = 0; k < 64; k++)
            {
                uint32_t qq = qpos[k];
                uint64_t w = bitmap[qq >> 6] & (uint64_t)(int64_t)(int32_t)qin[k];
                if ((w >> (qq & 63)) & 1)
                    direct_resieve_hit(pool, qq, j + (k & 15), primes[j + (k & 15)]);
            }

            if (pool->active == 0)
                resieve = 0;
        }

        if (sieve == NULL)
            continue;

        // advance to the next poly; the pair stays unordered
        for (q = 0; q < 4; q++)
        {
            v4t vp = V4_LOAD(primes + j + 4 * q);
            v4t u = V4_LOAD(ptr + 4 * q);
            v4t a = V4_LOAD(fr1 + j + 4 * q);
            v4t b = V4_LOAD(fr2 + j + 4 * q);

            if (sign > 0)
            {
                a = V4_ADD(V4_SUB(a, u), V4_AND(vp, V4_LT(a, u)));
                b = V4_ADD(V4_SUB(b, u), V4_AND(vp, V4_LT(b, u)));
            }
            else
            {
                a = V4_ADD(a, u);
                b = V4_ADD(b, u);
                a = V4_SUB(a, V4_ANDNOT(V4_LT(a, vp), vp));
                b = V4_SUB(b, V4_ANDNOT(V4_LT(b, vp), vp));
            }
            V4_STORE(fr1 + j + 4 * q, a);
            V4_STORE(fr2 + j + 4 * q, b);
        }

        // sieve the next poly's roots, both sides (as pass 1)
        logp = update_data.logp[j + 15] + DIRECT_LOGP_BIAS;
        V4_FN(direct_scatter16)(sieve, interval, primes + j, fr1 + j, fr2 + j,
            logp, kmode, hits, &nc, &flushed, &sink);
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

void V4_FN(direct_resieve_update_bitloop)(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    V4_FN(direct_resieve_update_body)(sconf, dconf, sieve, SIQS_K_BITLOOP);
}

void V4_FN(direct_resieve_update_sink)(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    V4_FN(direct_resieve_update_body)(sconf, dconf, sieve, SIQS_K_SINK);
}

void V4_FN(direct_resieve_update_creg)(static_conf_t* sconf, dynamic_conf_t* dconf, uint8_t* sieve)
{
    V4_FN(direct_resieve_update_body)(sconf, dconf, sieve, SIQS_K_CREG);
}

// ----------------------------------------------------------------------
// sparse sieve and its tdiv, [med_B, fb_blk8i)
// ----------------------------------------------------------------------
SIQS_FORCEINLINE void V4_FN(large_sieveblock_body)(uint8_t* sieve, sieve_fb* fb,
    fb_list* full_fb, const int kmode)
{
    uint32_t i;
    uint32_t prime, root1, root2;
    uint8_t logp;
    uint32_t bound = full_fb->fb_blk8i;
    v4t vblock = V4_SET1(32768);
    ALIGNED_MEM uint32_t cbuf[DIRECT_CBUF_SIZE];
    ALIGNED_MEM uint32_t t[4];
    uint32_t nc = 0;

    for (i = full_fb->med_B; (i + 16) <= bound; i += 16)
    {
        int h;
        v4t vtag;

        // approximate the next 16 logp's as equal to the largest
        logp = fb->logp[i + 15];
        vtag = V4_SET1((uint32_t)logp << 24);

        for (h = 0; h < 16; h += 4)
        {
            v4t vp = V4_LOAD(fb->prime + i + h);
            v4t v1 = V4_LOAD(fb->root1 + i + h);
            v4t v2 = V4_LOAD(fb->root2 + i + h);
            v4t m1 = V4_LT(v1, vblock);
            v4t m2 = V4_LT(v2, vblock);

            if (kmode == SIQS_K_BITLOOP)
            {
                uint32_t m;
                m = V4_MASK(m1);
                if (m)
                {
                    V4_STORE(t, v1);
                    while (m) { sieve[t[siqs_ctz32(m)]] -= logp; m &= m - 1; }
                }
                m = V4_MASK(m2);
                if (m)
                {
                    V4_STORE(t, v2);
                    while (m) { sieve[t[siqs_ctz32(m)]] -= logp; m &= m - 1; }
                }
            }
            else
            {
                nc = V4_FN(v4_cstore)(cbuf, nc, V4_MASK(m1), V4_OR(v1, vtag));
                nc = V4_FN(v4_cstore)(cbuf, nc, V4_MASK(m2), V4_OR(v2, vtag));
            }

            // advance the roots that hit, reduce both for the next block
            v1 = V4_SUB(V4_ADD(v1, V4_AND(vp, m1)), vblock);
            v2 = V4_SUB(V4_ADD(v2, V4_AND(vp, m2)), vblock);
            V4_STORE(fb->root1 + i + h, V4_MIN(v1, v2));
            V4_STORE(fb->root2 + i + h, V4_MAX(v1, v2));
        }

        if (nc >= DIRECT_CBUF_FLUSH)
        {
            direct_compress_flush(sieve, cbuf, nc);
            nc = 0;
        }
    }
    direct_compress_flush(sieve, cbuf, nc);

    // primes left over after the 16-groups; same as the AVX-512 tail
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

void V4_FN(large_sieveblock_bitloop)(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    V4_FN(large_sieveblock_body)(sieve, fb, full_fb, SIQS_K_BITLOOP);
}

void V4_FN(large_sieveblock_creg)(uint8_t* sieve, sieve_fb* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    V4_FN(large_sieveblock_body)(sieve, fb, full_fb, SIQS_K_CREG);
}

SIQS_FORCEINLINE void V4_FN(tdiv_large_body)(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf, const int kmode)
{
    uint32_t i;
    uint32_t tmp, prime, root1, root2, report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb* fbc = parity ? dconf->sieve_n : dconf->sieve_p;
    uint32_t block_loc;
    uint32_t bound = sconf->factor_base->fb_blk8i;

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {
        v4t vbl;
        ALIGNED_MEM uint32_t mbuf[64 + 16];
        uint32_t nm = 0, k;

        if (!dconf->valid_Qs[report_num])
            continue;

        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];
        vbl = V4_SET1(block_loc - 32768);

        for (i = sconf->factor_base->med_B; (i + 16) <= bound; i += 16)
        {
            uint32_t msk = 0;
            int h;

            for (h = 0; h < 4; h++)
            {
                v4t t = V4_ADD(V4_LOAD(fbc->prime + i + 4 * h), vbl);
                v4t e = V4_OR(V4_EQ(t, V4_LOAD(fbc->root1 + i + 4 * h)),
                    V4_EQ(t, V4_LOAD(fbc->root2 + i + 4 * h)));
                msk |= V4_MASK(e) << (4 * h);
            }

            if (kmode == SIQS_K_BITLOOP)
            {
                while (msk)
                {
                    uint32_t j = i + siqs_ctz32(msk);
                    DIVIDE_RESIEVED_PRIME_2(j);
                    msk &= msk - 1;
                }
            }
            else if (msk)
            {
                v4t vidx = V4_ADD(V4_SET1(i), V4_SETR(0, 1, 2, 3));
                for (h = 0; h < 4; h++)
                {
                    nm = V4_FN(v4_cstore)(mbuf, nm, (msk >> (4 * h)) & 15, vidx);
                    vidx = V4_ADD(vidx, V4_SET1(4));
                }
                if (nm >= 64)
                {
                    for (k = 0; k < nm; k++)
                    {
                        uint32_t j = mbuf[k];
                        DIVIDE_RESIEVED_PRIME_2(j);
                    }
                    nm = 0;
                }
            }
        }

        for (k = 0; k < nm; k++)
        {
            uint32_t j = mbuf[k];
            DIVIDE_RESIEVED_PRIME_2(j);
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
}

void V4_FN(tdiv_large_bitloop)(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    V4_FN(tdiv_large_body)(parity, poly_id, bnum, sconf, dconf, SIQS_K_BITLOOP);
}

void V4_FN(tdiv_large_creg)(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    V4_FN(tdiv_large_body)(parity, poly_id, bnum, sconf, dconf, SIQS_K_CREG);
}

// ----------------------------------------------------------------------
// per-poly roots for [med_B, B): vector sparse update + portable buckets
// ----------------------------------------------------------------------
void V4_FN(nextRoots_32k_lp)(static_conf_t* sconf, dynamic_conf_t* dconf)
{
    fb_list* fb = sconf->factor_base;
    update_t update_data = dconf->update_data;
    uint32_t bound = fb->B;
    uint32_t med_B = fb->med_B;
    uint32_t sparse_B = fb->fb_blk8i;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    uint32_t* ptr = (uint32_t*)&dconf->rootupdates[(v - 1) * bound + med_B];
    uint32_t* fr1 = (uint32_t*)update_data.firstroots1;
    uint32_t* fr2 = (uint32_t*)update_data.firstroots2;
    uint32_t j, root1, root2, prime;

    for (j = med_B; (j + 4) <= sparse_B; j += 4, ptr += 4)
    {
        v4t vp = V4_LOAD(update_data.prime + j);
        v4t r1 = V4_LOAD(fr1 + j);
        v4t r2 = V4_LOAD(fr2 + j);
        v4t u = V4_LOAD(ptr);
        v4t lo, hi;

        if (sign > 0)
        {
            r1 = V4_ADD(V4_SUB(r1, u), V4_AND(vp, V4_LT(r1, u)));
            r2 = V4_ADD(V4_SUB(r2, u), V4_AND(vp, V4_LT(r2, u)));
        }
        else
        {
            r1 = V4_ADD(r1, u);
            r2 = V4_ADD(r2, u);
            r1 = V4_SUB(r1, V4_ANDNOT(V4_LT(r1, vp), vp));
            r2 = V4_SUB(r2, V4_ANDNOT(V4_LT(r2, vp), vp));
        }

        lo = V4_MIN(r1, r2);
        hi = V4_MAX(r1, r2);
        V4_STORE(fr1 + j, lo);
        V4_STORE(fr2 + j, hi);
        V4_STORE(dconf->sieve_p->root1 + j, lo);
        V4_STORE(dconf->sieve_p->root2 + j, hi);
        V4_STORE(dconf->sieve_n->root1 + j, V4_SUB(vp, hi));
        V4_STORE(dconf->sieve_n->root2 + j, V4_SUB(vp, lo));
    }

    for (; j < sparse_B; j++, ptr++)
    {
        uint32_t uu = *ptr;
        prime = update_data.prime[j];
        root1 = fr1[j];
        root2 = fr2[j];

        if (sign > 0)
        {
            root1 = (uu > root1) ? root1 + prime - uu : root1 - uu;
            root2 = (uu > root2) ? root2 + prime - uu : root2 - uu;
        }
        else
        {
            root1 += uu;
            root2 += uu;
            if (root1 >= prime) root1 -= prime;
            if (root2 >= prime) root2 -= prime;
        }

        if (root2 < root1) { uint32_t t = root1; root1 = root2; root2 = t; }
        fr1[j] = root1;
        fr2[j] = root2;
        dconf->sieve_p->root1[j] = root1;
        dconf->sieve_p->root2[j] = root2;
        dconf->sieve_n->root1[j] = prime - root2;
        dconf->sieve_n->root2[j] = prime - root1;
    }

    nextRoots_32k_lp_buckets(sconf, dconf);
}
