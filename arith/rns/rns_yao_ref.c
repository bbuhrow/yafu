/*
 * RNS Montgomery modular multiplication, Yao regime: setup and scalar
 * reference.  See rns_yao.h and arith/rns_design.md §0b/§0c.
 *
 * The reference does each channel with GMP but follows the exact reduction
 * sequence of the vector kernels (fold until < 2^w), so it defines the
 * residues those kernels reproduce bit for bit when built with
 * RNSY_COLUMN_FOLD = 0.
 */

#include "rns_yao.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NCH RNSY_NCH
#define MAXN RNSY_MAXN

static void die(const char *msg)
{
    fprintf(stderr, "rnsy: %s\n", msg);
    exit(1);
}

// ---------------------------------------------------------------------------
// channel access and the fold reduction
// ---------------------------------------------------------------------------

static void get_ch(rnsy_fixed_t *fx, mpz_t X, uint64_t *x, int base, int ch)
{
    int t;
    mpz_set_ui(X, 0);
    for (t = fx->s - 1; t >= 0; t--)
    {
        mpz_mul_2exp(X, X, RNS_WBITS);
        mpz_add_ui(X, X, x[(base * fx->s + t) * fx->n + ch]);
    }
}

static void put_ch(rnsy_fixed_t *fx, uint64_t *x, int base, int ch, mpz_t X)
{
    mpz_t v;
    int t;
    if (mpz_sgn(X) < 0 || mpz_cmp(X, fx->two_w) >= 0)
        die("channel value out of range");
    mpz_init_set(v, X);
    for (t = 0; t < fx->s; t++)
    {
        x[(base * fx->s + t) * fx->n + ch] = mpz_get_ui(v) & RNS_MASK52;
        mpz_tdiv_q_2exp(v, v, RNS_WBITS);
    }
    mpz_clear(v);
}

// X -> (X mod 2^w) + (X >> w) * d until X < 2^w
static void fold(rnsy_fixed_t *fx, mpz_t X, uint64_t d, mpz_t tmp)
{
    while (mpz_cmp(X, fx->two_w) >= 0)
    {
        mpz_tdiv_q_2exp(tmp, X, (mp_bitcnt_t)fx->w);
        mpz_tdiv_r_2exp(X, X, (mp_bitcnt_t)fx->w);
        mpz_addmul_ui(X, tmp, d);
    }
}

static uint64_t mpz_low52(mpz_t x, mpz_t t)
{
    mpz_fdiv_r_2exp(t, x, RNS_WBITS);
    return mpz_get_ui(t);
}

// out += sum_k digit_k(a) * tab[k], a < 2^w: a product by a per-N constant
// through its rotations tab[k] = f * 2^(52k) mod m (the vector form)
static void cmul_acc(rnsy_fixed_t *fx, mpz_t out, mpz_t a, mpz_t *tab, mpz_t tmp)
{
    int k;
    for (k = 0; k < fx->s; k++)
    {
        mpz_tdiv_q_2exp(tmp, a, (mp_bitcnt_t)(RNS_WBITS * k));
        mpz_addmul_ui(out, tab[k], mpz_get_ui(tmp) & RNS_MASK52);
    }
}

// tab[k] = f * 2^(52k) mod m
static void rotations(rnsy_fixed_t *fx, mpz_t *tab, mpz_t f, mpz_t m)
{
    int k;
    for (k = 0; k < fx->s; k++)
    {
        mpz_mul_2exp(tab[k], f, (mp_bitcnt_t)(RNS_WBITS * k));
        mpz_mod(tab[k], tab[k], m);
    }
}

// ---------------------------------------------------------------------------
// fixed setup
// ---------------------------------------------------------------------------

static int capacity_s(int maxbits, int n)
{
    int s;
    for (s = 1; ; s++)
        if ((double)n * 52 * s - 0.01 >= maxbits + 2.0 * log2((double)RNS_HEADROOM * (n + 2)) + 1)
            return s;
}

// fit the constants t[a][b] (a = summed index, b = output channel) into 52
// bits: shift out the common trailing zeros of each output when needed
static int fit_constants(uint64_t out[MAXN][MAXN], int *z, mpz_t full[MAXN][MAXN],
    int n, int *vfull, int *vstored)
{
    int a, b, ok = 1;
    for (b = 0; b < n; b++)
    {
        int zmin = 1 << 30, need = 0;
        for (a = 0; a < n; a++)
        {
            int bits = (int)mpz_sizeinbase(full[a][b], 2);
            int tz = (int)mpz_scan1(full[a][b], 0);
            if (bits > *vfull)
                *vfull = bits;
            if (bits > 52)
                need = 1;
            if (tz < zmin)
                zmin = tz;
        }
        z[b] = need ? zmin : 0;
        for (a = 0; a < n; a++)
        {
            mpz_t v;
            mpz_init(v);
            mpz_tdiv_q_2exp(v, full[a][b], (mp_bitcnt_t)z[b]);
            if (mpz_sizeinbase(v, 2) > 52)
                ok = 0;
            else
                out[a][b] = mpz_get_ui(v);
            if ((int)mpz_sizeinbase(v, 2) > *vstored)
                *vstored = (int)mpz_sizeinbase(v, 2);
            mpz_clear(v);
        }
    }
    return ok;
}

rnsy_fixed_t *rnsy_fixed_init_n(int maxbits, int n)
{
    rnsy_fixed_t *fx;
    uint64_t dpool[2 * MAXN];
    mpz_t pool[2 * MAXN], full[MAXN][MAXN], t, u, bound;
    int s, i, j, k, cnt, ok = 1;
    uint64_t d;

    if (n < 2 || n > MAXN)
        return NULL;
    s = capacity_s(maxbits, n);
    if (s > RNSY_MAXS)
        return NULL;

    fx = (rnsy_fixed_t *)calloc(1, sizeof(rnsy_fixed_t));
    for (i = 0; i < MAXN; i++)
    {
        mpz_init(fx->b[i]); mpz_init(fx->c[i]);
        mpz_init(fx->negMC_b[i]); mpz_init(fx->K_b[i]);
        mpz_init(fx->g1_c[i]); mpz_init(fx->Cjinv_c[i]);
        mpz_init(fx->Bi[i]); mpz_init(fx->Biinv_b[i]);
        for (j = 0; j < MAXN; j++)
            mpz_init(full[i][j]);
    }
    mpz_init(fx->MB); mpz_init(fx->MC); mpz_init(fx->two_w);
    mpz_init(t); mpz_init(u); mpz_init(bound);

    fx->maxbits = maxbits;
    fx->n = n;
    fx->k = n + 2;
    fx->s = s;
    fx->w = 52 * s;
    fx->nwords = 2 * n * s + 1;
    fx->rword = 2 * n * s;
    mpz_setbit(fx->two_w, (mp_bitcnt_t)fx->w);

    // FCFS: 2^w - d, d odd ascending, pairwise coprime
    for (cnt = 0, d = 1; cnt < 2 * n; d += 2)
    {
        int good = 1;
        mpz_sub_ui(t, fx->two_w, d);
        for (i = 0; i < cnt && good; i++)
        {
            mpz_gcd(u, t, pool[i]);
            good = (mpz_cmp_ui(u, 1) == 0);
        }
        if (good)
        {
            mpz_init_set(pool[cnt], t);
            dpool[cnt++] = d;
        }
    }

    // split: B gets the smaller offsets, so every extension constant has
    // a uniform sign
    mpz_set_ui(fx->MB, 1);
    mpz_set_ui(fx->MC, 1);
    for (i = 0; i < n; i++)
    {
        fx->db[i] = dpool[i];
        fx->dc[i] = dpool[n + i];
        mpz_set(fx->b[i], pool[i]);
        mpz_set(fx->c[i], pool[n + i]);
        mpz_mul(fx->MB, fx->MB, fx->b[i]);
        mpz_mul(fx->MC, fx->MC, fx->c[i]);
    }
    for (i = 0; i < 2 * n; i++)
        mpz_clear(pool[i]);

    mpz_set_ui(bound, (unsigned long)(RNS_HEADROOM * fx->k));
    mpz_mul(bound, bound, bound);
    mpz_mul_2exp(bound, bound, (mp_bitcnt_t)maxbits);
    if (mpz_cmp(fx->MB, bound) < 0)
        die("M_B bound check failed");
    mpz_set_ui(bound, (unsigned long)(2 * fx->k));
    mpz_mul_2exp(bound, bound, (mp_bitcnt_t)maxbits);
    if (mpz_cmp(fx->MC, bound) < 0)
        die("M_C bound check failed");

    // first extension constants, full[i][j] = prod_{k!=i} (dc_j - db_k),
    // checked against the true CRT constants
    for (i = 0; i < n; i++)
    {
        mpz_divexact(fx->Bi[i], fx->MB, fx->b[i]);
        for (j = 0; j < n; j++)
        {
            mpz_set_ui(full[i][j], 1);
            for (k = 0; k < n; k++)
                if (k != i)
                    mpz_mul_ui(full[i][j], full[i][j], fx->dc[j] - fx->db[k]);
            mpz_sub(u, full[i][j], fx->Bi[i]);
            if (!mpz_divisible_p(u, fx->c[j]))
                die("first extension constant mismatch");
        }
    }
    fx->vB = fx->vBt = 0;
    ok &= fit_constants(fx->Bt, fx->zB, full, n, &fx->vB, &fx->vBt);

    // second extension, full[j][i] = prod_{k!=j} (dc_k - db_i): n-1 factors
    // (c_k - b_i) = -(dc_k - db_i), so C_j == (-1)^(n-1) full (mod b_i)
    fx->cneg = (n - 1) & 1;
    for (j = 0; j < n; j++)
    {
        mpz_t Cj;
        mpz_init(Cj);
        mpz_divexact(Cj, fx->MC, fx->c[j]);
        for (i = 0; i < n; i++)
        {
            mpz_set_ui(full[j][i], 1);
            for (k = 0; k < n; k++)
                if (k != j)
                    mpz_mul_ui(full[j][i], full[j][i], fx->dc[k] - fx->db[i]);
            if (fx->cneg)
                mpz_add(u, full[j][i], Cj);
            else
                mpz_sub(u, full[j][i], Cj);
            if (!mpz_divisible_p(u, fx->b[i]))
                die("second extension constant mismatch");
        }
        fx->Cr[j] = mpz_low52(Cj, t);
        mpz_invert(fx->Cjinv_c[j], Cj, fx->c[j]);
        mpz_invert(fx->g1_c[j], fx->MB, fx->c[j]);
        mpz_mul(fx->g1_c[j], fx->g1_c[j], Cj);
        mpz_mod(fx->g1_c[j], fx->g1_c[j], fx->c[j]);
        mpz_clear(Cj);
    }
    fx->vC = fx->vCt = 0;
    // full[j][i] is indexed (summed j, output i): store as Ct[j][i]
    ok &= fit_constants(fx->Ct, fx->zC, full, n, &fx->vC, &fx->vCt);

    fx->shifted = 0;
    for (i = 0; i < n; i++)
        if (fx->zB[i] || fx->zC[i])
            fx->shifted = 1;

    for (i = 0; i < n; i++)
    {
        fx->Br[i] = mpz_low52(fx->Bi[i], t);
        mpz_invert(fx->Biinv_b[i], fx->Bi[i], fx->b[i]);

        mpz_neg(t, fx->MC);
        mpz_mod(fx->negMC_b[i], t, fx->b[i]);

        // K_i = ceil(2^(w + vC + ceil(log2 n)) / b_i) * b_i > sum_j xi'_j |C_j|~
        mpz_set_ui(t, 0);
        mpz_setbit(t, (mp_bitcnt_t)(fx->w + fx->vC + (int)ceil(log2((double)n))));
        mpz_cdiv_q(fx->K_b[i], t, fx->b[i]);
        mpz_mul(fx->K_b[i], fx->K_b[i], fx->b[i]);
    }

    mpz_set_ui(u, 0);
    mpz_setbit(u, RNS_WBITS);
    mpz_invert(t, fx->MB, u);
    fx->MBinv_r = mpz_get_ui(t);
    mpz_invert(t, fx->MC, u);
    fx->MCinv_r = mpz_get_ui(t);

    for (i = 0; i < MAXN; i++)
        for (j = 0; j < MAXN; j++)
            mpz_clear(full[i][j]);
    mpz_clear(t); mpz_clear(u); mpz_clear(bound);

    if (!ok)
    {
        rnsy_fixed_free(fx);
        return NULL;
    }
    return fx;
}

rnsy_fixed_t *rnsy_fixed_init(int maxbits)
{
    rnsy_fixed_t *fx = rnsy_fixed_init_n(maxbits, RNSY_NCH);
    if (fx == NULL)
        die("no 8-channel setup for this size");
    return fx;
}

// Batch-kernel cost per 8 numbers, in vector instructions: 5 channel
// products plus folds per channel pair, two n x n extensions of s lo/hi
// pairs.  The linear coefficients were fit to a measured n sweep (Rocket
// Lake, 512-4096 bits): the pick is within ~5% (noise level) of the
// fastest n at every size.  The trailing-zero shift costs too little to
// show up, so shifted setups get no penalty.
static double batch_cost(int n, int s, int shifted)
{
    (void)shifted;
    return n * (5.0 * s * s + 4.0 * s + 14.0) + 4.0 * n * n * s;
}

int rnsyb_choose_n(int maxbits)
{
    int n, best = 0;
    double bestc = 1e30;
    for (n = 2; n <= MAXN; n++)
    {
        rnsy_fixed_t *fx;
        int s = capacity_s(maxbits, n);
        double c;
        if (s > RNSY_MAXS)
            continue;
        c = batch_cost(n, s, 0);
        if (c >= bestc)
            continue;           // can't win even without shifts
        if ((fx = rnsy_fixed_init_n(maxbits, n)) == NULL)
            continue;
        c = batch_cost(n, s, fx->shifted);
        if (c < bestc)
        {
            bestc = c;
            best = n;
        }
        rnsy_fixed_free(fx);
    }
    return best;
}

void rnsy_fixed_free(rnsy_fixed_t *fx)
{
    int i;
    for (i = 0; i < MAXN; i++)
    {
        mpz_clear(fx->b[i]); mpz_clear(fx->c[i]);
        mpz_clear(fx->negMC_b[i]); mpz_clear(fx->K_b[i]);
        mpz_clear(fx->g1_c[i]); mpz_clear(fx->Cjinv_c[i]);
        mpz_clear(fx->Bi[i]); mpz_clear(fx->Biinv_b[i]);
    }
    mpz_clear(fx->MB); mpz_clear(fx->MC); mpz_clear(fx->two_w);
    free(fx);
}

// ---------------------------------------------------------------------------
// variable setup
// ---------------------------------------------------------------------------

uint64_t *rnsy_alloc(rnsy_ctx_t *ctx)
{
    return rns_aligned_alloc((size_t)ctx->fx->nwords);
}

// digits of X (radix 2^52, ndig of them) into v[k * stride], k < ndig
static void digits_mpz(uint64_t *v, int ndig, int stride, mpz_t X)
{
    mpz_t t;
    int k;
    mpz_init_set(t, X);
    for (k = 0; k < ndig; k++)
    {
        v[k * stride] = mpz_get_ui(t) & RNS_MASK52;
        mpz_tdiv_q_2exp(t, t, RNS_WBITS);
    }
    if (mpz_sgn(t) != 0)
        die("vector constant does not fit");
    mpz_clear(t);
}

// one-number kernel constants (n = 8, unshifted only)
static void build_vconst(rnsy_ctx_t *ctx)
{
    rnsy_fixed_t *fx = ctx->fx;
    int s = fx->s, nvec = 7 * s + 23 + 2 * s + 3 * s * s, i, j, k;
    uint64_t *p;
    mpz_t t;

    if (fx->n != NCH || fx->shifted || !fx->cneg)
    {
        ctx->vblock = NULL;
        return;
    }

    mpz_init(t);
    p = ctx->vblock = rns_aligned_alloc((size_t)nvec * NCH);
    ctx->vdb = p; p += NCH;
    ctx->vdc = p; p += NCH;
    ctx->vf1 = p; p += s * NCH;
    ctx->vg1 = p; p += s * NCH;
    ctx->vg2 = p; p += s * NCH;
    ctx->vnegMC = p; p += s * NCH;
    ctx->vK = p; p += (s + 1) * NCH;
    ctx->vEb = p; p += (s + 1) * NCH;
    ctx->vEc = p; p += (s + 1) * NCH;
    ctx->vBrN = p; p += NCH;
    ctx->vCr = p; p += NCH;
    ctx->vbt = p; p += NCH * NCH;
    ctx->vct = p; p += NCH * NCH;
    ctx->vscratch = p; p += 2 * s * NCH;
    ctx->vF1 = p; p += s * s * NCH;
    ctx->vG1 = p; p += s * s * NCH;
    ctx->vG2 = p; p += s * s * NCH;

    for (i = 0; i < NCH; i++)
    {
        ctx->vdb[i] = fx->db[i];
        ctx->vdc[i] = fx->dc[i];
        digits_mpz(ctx->vf1 + i, s, NCH, ctx->f1_b[i]);
        digits_mpz(ctx->vg1 + i, s, NCH, fx->g1_c[i]);
        digits_mpz(ctx->vg2 + i, s, NCH, ctx->g2_c[i]);
        digits_mpz(ctx->vnegMC + i, s, NCH, fx->negMC_b[i]);
        digits_mpz(ctx->vK + i, s + 1, NCH, fx->K_b[i]);

        // E = subconst + 2m, per channel, for x + E - y in submod
        get_ch(fx, t, ctx->subconst, 0, i);
        mpz_addmul_ui(t, fx->b[i], 2);
        digits_mpz(ctx->vEb + i, s + 1, NCH, t);
        get_ch(fx, t, ctx->subconst, 1, i);
        mpz_addmul_ui(t, fx->c[i], 2);
        digits_mpz(ctx->vEc + i, s + 1, NCH, t);

        for (k = 0; k < s; k++)
        {
            digits_mpz(ctx->vF1 + k * s * NCH + i, s, NCH, ctx->F1[i][k]);
            digits_mpz(ctx->vG1 + k * s * NCH + i, s, NCH, ctx->G1[i][k]);
            digits_mpz(ctx->vG2 + k * s * NCH + i, s, NCH, ctx->G2[i][k]);
        }

        ctx->vBrN[i] = ctx->BrN[i];
        ctx->vCr[i] = fx->Cr[i];
        for (j = 0; j < NCH; j++)
        {
            ctx->vbt[i * NCH + j] = fx->Bt[i][j];
            ctx->vct[i * NCH + j] = fx->Ct[i][j];   // vct[j] lane i = Ct[j][i]
        }
    }
    mpz_clear(t);
}

// batch kernel constants
static void build_bconst(rnsy_ctx_t *ctx)
{
    rnsy_fixed_t *fx = ctx->fx;
    int n = fx->n, s = fx->s, i, k, a;
    const size_t zn = (size_t)n, zs = (size_t)s;
    size_t words = zn * (zs + 1) * 16 + zn * zs * NCH +
        3 * zn * zs * zs + zn * (zs + 2) + 2 * zn * (zs + 1);
    uint64_t *p;
    mpz_t t;

    mpz_init(t);
    p = ctx->bblock = rns_aligned_alloc(words);
    ctx->balpha = p; p += zn * (zs + 1) * 16;     // keeps 64-byte alignment
    ctx->bscratch = p; p += zn * zs * NCH;
    ctx->bF1 = p; p += zn * zs * zs;
    ctx->bG1 = p; p += zn * zs * zs;
    ctx->bG2 = p; p += zn * zs * zs;
    ctx->bK = p; p += zn * (zs + 2);
    ctx->bEb = p; p += zn * (zs + 1);
    ctx->bEc = p; p += zn * (zs + 1);

    for (i = 0; i < n; i++)
    {
        for (k = 0; k < s; k++)
        {
            digits_mpz(ctx->bF1 + (i * s + k) * s, s, 1, ctx->F1[i][k]);
            digits_mpz(ctx->bG1 + (i * s + k) * s, s, 1, ctx->G1[i][k]);
            digits_mpz(ctx->bG2 + (i * s + k) * s, s, 1, ctx->G2[i][k]);
        }
        digits_mpz(ctx->bK + i * (s + 2), s + 2, 1, fx->K_b[i]);
        get_ch(fx, t, ctx->subconst, 0, i);
        mpz_addmul_ui(t, fx->b[i], 2);
        digits_mpz(ctx->bEb + i * (s + 1), s + 1, 1, t);
        get_ch(fx, t, ctx->subconst, 1, i);
        mpz_addmul_ui(t, fx->c[i], 2);
        digits_mpz(ctx->bEc + i * (s + 1), s + 1, 1, t);

        // alpha' <= n + 1 < 16: table of alpha * |-M_C|_{b_i}, s+1 digits
        for (a = 0; a < 16; a++)
        {
            mpz_mul_ui(t, fx->negMC_b[i], (unsigned long)((a <= n + 1) ? a : 0));
            digits_mpz(ctx->balpha + (i * (s + 1)) * 16 + a, s + 1, 16, t);
        }
    }
    mpz_clear(t);
}

rnsy_ctx_t *rnsy_ctx_init(rnsy_fixed_t *fx, mpz_t N)
{
    rnsy_ctx_t *ctx;
    mpz_t t, u;
    int i, j;

    if (mpz_even_p(N) || (int)mpz_sizeinbase(N, 2) > fx->maxbits)
    {
        fprintf(stderr, "rnsy: N must be odd and < 2^%d\n", fx->maxbits);
        return NULL;
    }
    mpz_init(t); mpz_init(u);
    mpz_gcd(t, N, fx->MB);
    mpz_gcd(u, N, fx->MC);
    if (mpz_cmp_ui(t, 1) != 0 || mpz_cmp_ui(u, 1) != 0)
    {
        mpz_clear(t); mpz_clear(u);
        return NULL;
    }

    ctx = (rnsy_ctx_t *)calloc(1, sizeof(rnsy_ctx_t));
    ctx->fx = fx;
    mpz_init_set(ctx->N, N);
    for (i = 0; i < RNSY_MAXN; i++)
    {
        mpz_init(ctx->f1_b[i]);
        mpz_init(ctx->g2_c[i]);
        for (j = 0; j < RNSY_MAXS; j++)
        {
            mpz_init(ctx->F1[i][j]);
            mpz_init(ctx->G1[i][j]);
            mpz_init(ctx->G2[i][j]);
        }
    }

    mpz_set_ui(u, 0);
    mpz_setbit(u, RNS_WBITS);
    mpz_invert(t, fx->MB, u);
    mpz_mul(t, t, N);
    for (i = 0; i < fx->n; i++)
    {
        mpz_mul_ui(u, t, fx->Br[i]);
        ctx->BrN[i] = mpz_low52(u, u);
    }

    for (i = 0; i < fx->n; i++)
    {
        mpz_invert(ctx->f1_b[i], N, fx->b[i]);
        mpz_mul(ctx->f1_b[i], ctx->f1_b[i], fx->Biinv_b[i]);
        mpz_neg(ctx->f1_b[i], ctx->f1_b[i]);
        mpz_mod(ctx->f1_b[i], ctx->f1_b[i], fx->b[i]);

        mpz_invert(ctx->g2_c[i], fx->MB, fx->c[i]);
        mpz_mul(ctx->g2_c[i], ctx->g2_c[i], N);
        mpz_mul(ctx->g2_c[i], ctx->g2_c[i], fx->Cjinv_c[i]);
        mpz_mod(ctx->g2_c[i], ctx->g2_c[i], fx->c[i]);

        rotations(fx, ctx->F1[i], ctx->f1_b[i], fx->b[i]);
        rotations(fx, ctx->G1[i], fx->g1_c[i], fx->c[i]);
        rotations(fx, ctx->G2[i], ctx->g2_c[i], fx->c[i]);
    }

    mpz_init(ctx->MBmodN);
    mpz_init(ctx->MBinvmodN);
    mpz_mod(ctx->MBmodN, fx->MB, N);
    mpz_invert(ctx->MBinvmodN, fx->MB, N);

    ctx->subconst = rnsy_alloc(ctx);
    mpz_mul_ui(t, N, (unsigned long)(RNS_SUBMULT * fx->k));
    rnsy_from_int(ctx, ctx->subconst, t);
    ctx->one = rnsy_alloc(ctx);
    mpz_set_ui(t, 1);
    rnsy_to_monty(ctx, ctx->one, t);
    build_vconst(ctx);
    build_bconst(ctx);

    mpz_clear(t); mpz_clear(u);
    return ctx;
}

void rnsy_ctx_free(rnsy_ctx_t *ctx)
{
    int i, j;
    for (i = 0; i < RNSY_MAXN; i++)
    {
        mpz_clear(ctx->f1_b[i]);
        mpz_clear(ctx->g2_c[i]);
        for (j = 0; j < RNSY_MAXS; j++)
        {
            mpz_clear(ctx->F1[i][j]);
            mpz_clear(ctx->G1[i][j]);
            mpz_clear(ctx->G2[i][j]);
        }
    }
    mpz_clear(ctx->N); mpz_clear(ctx->MBmodN); mpz_clear(ctx->MBinvmodN);
    rns_vfree(ctx->subconst);
    rns_vfree(ctx->one);
    if (ctx->vblock)
        rns_vfree(ctx->vblock);
    rns_vfree(ctx->bblock);
    free(ctx);
}

// ---------------------------------------------------------------------------
// conversions
// ---------------------------------------------------------------------------

void rnsy_from_int(rnsy_ctx_t *ctx, uint64_t *x, mpz_t X)
{
    rnsy_fixed_t *fx = ctx->fx;
    mpz_t t;
    int i;

    mpz_init(t);
    for (i = 0; i < fx->n; i++)
    {
        mpz_mod(t, X, fx->b[i]);
        put_ch(fx, x, 0, i, t);
        mpz_mod(t, X, fx->c[i]);
        mpz_mul(t, t, fx->Cjinv_c[i]);
        mpz_mod(t, t, fx->c[i]);
        put_ch(fx, x, 1, i, t);
    }
    x[fx->rword] = mpz_low52(X, t);
    mpz_clear(t);
}

void rnsy_to_int(rnsy_ctx_t *ctx, mpz_t X, uint64_t *x)
{
    rnsy_fixed_t *fx = ctx->fx;
    mpz_t t;
    int i;

    mpz_init(t);
    mpz_set_ui(X, 0);
    for (i = 0; i < fx->n; i++)
    {
        get_ch(fx, t, x, 0, i);
        mpz_mul(t, t, fx->Biinv_b[i]);
        mpz_mod(t, t, fx->b[i]);
        mpz_addmul(X, t, fx->Bi[i]);
    }
    mpz_mod(X, X, fx->MB);
    mpz_clear(t);
}

void rnsy_to_monty(rnsy_ctx_t *ctx, uint64_t *x, mpz_t a)
{
    mpz_t t;
    mpz_init(t);
    mpz_mul(t, a, ctx->MBmodN);
    mpz_mod(t, t, ctx->N);
    rnsy_from_int(ctx, x, t);
    mpz_clear(t);
}

void rnsy_from_monty(rnsy_ctx_t *ctx, mpz_t a, uint64_t *x)
{
    rnsy_to_int(ctx, a, x);
    mpz_mul(a, a, ctx->MBinvmodN);
    mpz_mod(a, a, ctx->N);
}

// batch layout: word (base, ch, t) of lane l at ((base*n + ch)*s + t)*8 + l
uint64_t *rnsyb_alloc(rnsy_ctx_t *ctx)
{
    return rns_aligned_alloc((size_t)ctx->fx->nwords * NCH);
}

void rnsyb_vset_lane(rnsy_ctx_t *ctx, uint64_t *v, int lane, uint64_t *x)
{
    rnsy_fixed_t *fx = ctx->fx;
    int n = fx->n, s = fx->s, base, ch, t;
    for (base = 0; base < 2; base++)
        for (ch = 0; ch < n; ch++)
            for (t = 0; t < s; t++)
                v[((base * n + ch) * s + t) * NCH + lane] = x[(base * s + t) * n + ch];
    v[fx->rword * NCH + lane] = x[fx->rword];
}

void rnsyb_vget_lane(rnsy_ctx_t *ctx, uint64_t *x, uint64_t *v, int lane)
{
    rnsy_fixed_t *fx = ctx->fx;
    int n = fx->n, s = fx->s, base, ch, t;
    for (base = 0; base < 2; base++)
        for (ch = 0; ch < n; ch++)
            for (t = 0; t < s; t++)
                x[(base * s + t) * n + ch] = v[((base * n + ch) * s + t) * NCH + lane];
    x[fx->rword] = v[fx->rword * NCH + lane];
}

// ---------------------------------------------------------------------------
// scalar reference
// ---------------------------------------------------------------------------

void rnsy_mulmod_ref(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    rnsy_fixed_t *fx = ctx->fx;
    int n = fx->n;
    mpz_t xi[MAXN], xr[MAXN], q, a, b, tmp;
    uint64_t ur, sr, alpha, r_x = x[fx->rword], r_y = y[fx->rword];
    int i, j;

    mpz_init(q); mpz_init(a); mpz_init(b); mpz_init(tmp);
    for (i = 0; i < n; i++)
    {
        mpz_init(xi[i]);
        mpz_init(xr[i]);
    }

    // 1: xi_i = fold(fold(x_i y_i) f1_i), in B
    for (i = 0; i < n; i++)
    {
        get_ch(fx, a, x, 0, i);
        get_ch(fx, b, y, 0, i);
        mpz_mul(a, a, b);
        fold(fx, a, fx->db[i], tmp);
        mpz_set_ui(xi[i], 0);
        cmul_acc(fx, xi[i], a, ctx->F1[i], tmp);
        fold(fx, xi[i], fx->db[i], tmp);
    }

    // 2-3: Q_j = fold(sum_i xi_i |B_i|~);  xi'_j = fold(T'_j g1_j + Q_j g2_j)
    for (j = 0; j < n; j++)
    {
        mpz_set_ui(q, 0);
        for (i = 0; i < n; i++)
            mpz_addmul_ui(q, xi[i], fx->Bt[i][j]);
        mpz_mul_2exp(q, q, (mp_bitcnt_t)fx->zB[j]);
        fold(fx, q, fx->dc[j], tmp);

        get_ch(fx, a, x, 1, j);
        get_ch(fx, b, y, 1, j);
        mpz_mul(a, a, b);
        fold(fx, a, fx->dc[j], tmp);
        mpz_set_ui(xr[j], 0);
        cmul_acc(fx, xr[j], a, ctx->G1[j], tmp);
        cmul_acc(fx, xr[j], q, ctx->G2[j], tmp);
        fold(fx, xr[j], fx->dc[j], tmp);
    }

    // 4: redundant channel and alpha'
    ur = r_x * r_y * fx->MBinv_r;
    for (i = 0; i < n; i++)
        ur += mpz_getlimbn(xi[i], 0) * ctx->BrN[i];
    ur &= RNS_MASK52;
    sr = 0;
    for (j = 0; j < n; j++)
        sr += mpz_getlimbn(xr[j], 0) * fx->Cr[j];
    alpha = ((sr - ur) * fx->MCinv_r) & RNS_MASK52;
    if (alpha > (uint64_t)n + 1)
        die("alpha' out of range");

    // 5: U_i = fold(+-sum_j xi'_j |C_j|~ + alpha' |-M_C|_{b_i}), the minus
    // sign (n even) taken as K_i - sum
    for (i = 0; i < n; i++)
    {
        mpz_set_ui(b, 0);
        for (j = 0; j < n; j++)
            mpz_addmul_ui(b, xr[j], fx->Ct[j][i]);
        mpz_mul_2exp(b, b, (mp_bitcnt_t)fx->zC[i]);
        if (fx->cneg)
            mpz_sub(q, fx->K_b[i], b);
        else
            mpz_set(q, b);
        if (mpz_sgn(q) < 0)
            die("K_b too small");
        mpz_addmul_ui(q, fx->negMC_b[i], alpha);
        fold(fx, q, fx->db[i], tmp);
        put_ch(fx, z, 0, i, q);
    }
    for (j = 0; j < n; j++)
        put_ch(fx, z, 1, j, xr[j]);
    z[fx->rword] = ur;

    for (i = 0; i < n; i++)
    {
        mpz_clear(xi[i]);
        mpz_clear(xr[i]);
    }
    mpz_clear(q); mpz_clear(a); mpz_clear(b); mpz_clear(tmp);
}

static void addsub_ref(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y, int sub)
{
    rnsy_fixed_t *fx = ctx->fx;
    mpz_t a, b, tmp;
    int base, ch;

    mpz_init(a); mpz_init(b); mpz_init(tmp);
    for (base = 0; base < 2; base++)
    {
        for (ch = 0; ch < fx->n; ch++)
        {
            uint64_t d = base ? fx->dc[ch] : fx->db[ch];
            get_ch(fx, a, x, base, ch);
            get_ch(fx, b, y, base, ch);
            if (sub)
            {
                // x + SUBMULT*k*N + 2m - y, 2m > 2^w > y
                get_ch(fx, tmp, ctx->subconst, base, ch);
                mpz_add(a, a, tmp);
                mpz_addmul_ui(a, base ? fx->c[ch] : fx->b[ch], 2);
                mpz_sub(a, a, b);
            }
            else
                mpz_add(a, a, b);
            fold(fx, a, d, tmp);
            put_ch(fx, z, base, ch, a);
        }
    }
    if (sub)
        z[fx->rword] = (x[fx->rword] + ctx->subconst[fx->rword] - y[fx->rword]) & RNS_MASK52;
    else
        z[fx->rword] = (x[fx->rword] + y[fx->rword]) & RNS_MASK52;
    mpz_clear(a); mpz_clear(b); mpz_clear(tmp);
}

void rnsy_addmod_ref(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    addsub_ref(ctx, z, x, y, 0);
}

void rnsy_submod_ref(rnsy_ctx_t *ctx, uint64_t *z, uint64_t *x, uint64_t *y)
{
    addsub_ref(ctx, z, x, y, 1);
}
