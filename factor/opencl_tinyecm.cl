/*============================================================================
 * opencl_tinyecm.cl  -- OpenCL kernel translation of cuda_tinyecm.cu
 *
 * Translation notes
 * -----------------
 * CUDA -> OpenCL concept mapping used here:
 *   __device__          -> (inline function, no annotation needed)
 *   __global__ void     -> __kernel void
 *   __constant T arr[]  -> __constant T arr[]   (same keyword, same meaning)
 *   blockIdx.x          -> get_group_id(0)
 *   threadIdx.x         -> get_local_id(0)
 *   blockDim.x          -> get_local_size(0)
 *   blockIdx.x*blockDim.x+threadIdx.x -> get_global_id(0)
 *   uint / ulong -> uint / ulong  (OpenCL built-in types)
 *   uint32 (bare)       -> uint
 *
 * This file includes opencl_intrinsics.cl for all the multiprecision
 * helpers (montmul96, modsub96, etc.) and adds the ECM-specific functions.
 *
 * The 64-bit ECM kernel (gbl_ecm) and the 96-bit ECM kernel (gbl_ecm96)
 * are both translated here.  The P-1 kernel (gbl_pm196) is also included.
 *
 * The #if 0 blocks in the original (dead code) are omitted for clarity.
 * USE_D30 is the active stage-2 variant (matching the original default).
 *============================================================================*/

#include "opencl_intrinsics.cl"

#define USE_D30

/* -------------------------------------------------------------------------
 * Shared types
 * ------------------------------------------------------------------------- */

typedef struct { ulong X; ulong Z; } uecm_pt;

typedef struct {
    uint X[3];
    uint Z[3];
} uecm96_pt;

#define INV_2_POW_32 2.3283064365386962890625e-10

/* =========================================================================
 * Utility / RNG
 * ========================================================================= */

static inline uint
lcg_rand_32b(uint lower, uint upper, ulong *ploc_lcg)
{
    *ploc_lcg = 6364136223846793005UL * (*ploc_lcg) + 1442695040888963407UL;
    return lower + (uint)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}

/* =========================================================================
 * 64-bit ECM point operations
 * ========================================================================= */

static inline void
uaddxz(uint rho, ulong n,
        ulong p1x, ulong p1z, ulong p2x, ulong p2z,
        ulong pix, ulong piz,
        ulong *pox, ulong *poz)
{
    ulong diff1 = modsub64(p1x, p1z, n);
    ulong sum1  = modadd64(p1x, p1z, n);
    ulong diff2 = modsub64(p2x, p2z, n);
    ulong sum2  = modadd64(p2x, p2z, n);

    ulong tt1 = montmul64(diff1, sum2, n, rho);   // U
    ulong tt2 = montmul64(sum1,  diff2, n, rho);  // V

    ulong tt3 = modadd64(tt1, tt2, n);
    ulong tt4 = modsub64(tt1, tt2, n);
    tt1 = montmul64(tt3, tt3, n, rho);   // (U+V)^2
    tt2 = montmul64(tt4, tt4, n, rho);   // (U-V)^2

    *pox = montmul64(tt1, piz, n, rho);  // Z*(U+V)^2
    *poz = montmul64(tt2, pix, n, rho);  // x*(U-V)^2
}

static inline void
uadd(uint rho, ulong n, uecm_pt p1, uecm_pt p2, uecm_pt pi, uecm_pt *po)
{
    ulong diff1 = modsub64(p1.X, p1.Z, n);
    ulong sum1 = modadd64(p1.X, p1.Z, n);
    ulong diff2 = modsub64(p2.X, p2.Z, n);
    ulong sum2 = modadd64(p2.X, p2.Z, n);

    ulong tt1 = montmul64(diff1, sum2, n, rho); //U
    ulong tt2 = montmul64(sum1, diff2, n, rho); //V

    ulong tt3 = modadd64(tt1, tt2, n);
    ulong tt4 = modsub64(tt1, tt2, n);
    tt1 = montmul64(tt3, tt3, n, rho);   //(U + V)^2
    tt2 = montmul64(tt4, tt4, n, rho);   //(U - V)^2

    ulong tmpx = montmul64(tt1, pi.Z, n, rho);     //Z * (U + V)^2
    ulong tmpz = montmul64(tt2, pi.X, n, rho);     //x * (U - V)^2

    po->X = tmpx;
    po->Z = tmpz;
}

static inline void
udup(ulong s, uint rho, ulong n,
     ulong insum, ulong indiff, uecm_pt *P)
{
    ulong tt1 = montmul64(indiff, indiff, n, rho);  // U=(x1-z1)^2
    ulong tt2 = montmul64(insum,  insum,  n, rho);  // V=(x1+z1)^2
    P->X = montmul64(tt1, tt2, n, rho);             // x=U*V

    ulong tt3 = modsub64(tt2, tt1, n);              // w=V-U
    tt2 = montmul64(tt3, s,  n, rho);               // w=(A+2)/4*w
    tt2 = modadd64(tt2, tt1, n);                    // w=w+U

    P->Z = montmul64(tt2, tt3, n, rho);             // Z=w*(V-U)
}

static inline void
uprac(uint rho, ulong n, uecm_pt *P, ulong c, double v, ulong s)
{
    // require postive odd c
    ulong d, e, r;

    d = c;
    r = (ulong)((double)d * v + 0.5);
    d = c - r;
    e = 2 * r - c;

    ulong s1, d1;
    ulong swp;
    uecm_pt pt1, pt2, pt3;

    // the first one is always a doubling
    // point1 is [1]P
    pt1.X = pt2.X = pt3.X = P->X;
    pt1.Z = pt2.Z = pt3.Z = P->Z;

    d1 = modsub64(pt1.X, pt1.Z, n);
    s1 = modadd64(pt1.X, pt1.Z, n);

    // point2 is [2]P
    udup(s, rho, n, s1, d1, &pt1);

    while (d != e)
    {
        if (d < e)
        {
            r = d;
            d = e;
            e = r;
            swp = pt1.X;
            pt1.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt2.Z;
            pt2.Z = swp;
        }

        // in our small-B1 cases these are the only
        // PRAC conditions used.  Need to verify that
        // continues to be the case if B1 grows.
        if ((d + 3) / 4 <= e)
        {
            d -= e;

            uecm_pt pt4;
            uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = swp;
        }
        else //if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;

            d1 = modsub64(pt1.X, pt1.Z, n);
            s1 = modadd64(pt1.X, pt1.Z, n);

            uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            udup(s, rho, n, s1, d1, &pt1);        // A = 2A
        }
    }

    uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
    return;
}

static inline ulong
uecm_build(uecm_pt *P, uint rho, ulong n,
           uint sigma, ulong *ps, ulong one, ulong Rsqr)
{
    ulong t1, t2, t3, t4, t5;
    ulong u, v;

    t2 = modadd64(one, one, n);
    t4 = modadd64(t2, t2, n);
    t5 = modadd64(one, t4, n);

    u = montmul64((ulong)sigma, Rsqr, n, rho);  // to_monty(sigma)
    v = modadd64(u, u, n);
    v = modadd64(v, v, n);            // 4*sigma
    u = montmul64(u, u, n, rho);
    u = modsub64(u, t5, n);           // sigma^2 - 5
    t1 = montmul64(u, u, n, rho);
    ulong tmpx = montmul64(t1, u, n, rho);  // u^3

    t2 = modadd64(v, v, n);             // 2*v
    t2 = modadd64(t2, t2, n);           // 4*v
    t2 = modadd64(t2, t2, n);           // 8*v
    t2 = modadd64(t2, t2, n);           // 16*v
    t5 = montmul64(t2, tmpx, n, rho);    // 16*u^3*v
    t1 = montmul64(v, v, n, rho);
    ulong tmpz = montmul64(t1, v, n, rho);  // v^3

    //compute parameter A
    t1 = modsub64(v, u, n);           // (v - u)
    t2 = montmul64(t1, t1, n, rho);
    t4 = montmul64(t2, t1, n, rho);   // (v - u)^3
    t1 = modadd64(u, u, n);           // 2u
    t2 = modadd64(u, v, n);           // u + v
    t3 = modadd64(t1, t2, n);         // 3u + v
    t1 = montmul64(t3, t4, n, rho);   // a = (v-u)^3 * (3u + v)

    // accomplish the division by multiplying by the modular inverse
    t2 = 1;
    t5 = montmul64(t5, t2, n, rho);   // take t5 out of monty rep
    t3 = modinv64(t5, n, &t2);
    t3 = montmul64(t3, Rsqr, n, rho); // to_monty(t3)
    *ps = montmul64(t3, t1, n, rho);

    P->X = tmpx;
    P->Z = tmpz;

    return t2;
}

static inline void
uecm_stage1(uint rho, ulong n, uecm_pt *P, uint stg1, ulong s)
{
    uint q = 2;
    while (q < (stg1 * 4)) {
        ulong diff1 = modsub64(P->X, P->Z, n);
        ulong sum1  = modadd64(P->X, P->Z, n);
        udup(s, rho, n, sum1, diff1, P);
        q *= 2;
    }

    if (stg1 > 205)
    {
        // anything greater than 205 will use B1=250
        uprac(rho, n, P, 211, 0.612429949509495031, s);
        uprac(rho, n, P, 223, 0.625306711365725132, s);
        uprac(rho, n, P, 227, 0.580178728295464130, s);
        uprac(rho, n, P, 229, 0.548409048446403258, s);
        uprac(rho, n, P, 233, 0.618033988749894903, s);
        uprac(rho, n, P, 239, 0.618033988749894903, s);
        uprac(rho, n, P, 241, 0.625306711365725132, s);
        uprac(rho, n, P, 3, 0.618033988749894903, s);
    }

    if (stg1 > 175)
    {
        // anything greater than 175 will use B1=200
        uprac(rho, n, P, 3, 0.618033988749894903, s);
        uprac(rho, n, P, 3, 0.618033988749894903, s);
        uprac(rho, n, P, 3, 0.618033988749894903, s);
        uprac(rho, n, P, 3, 0.618033988749894903, s);
        uprac(rho, n, P, 5, 0.618033988749894903, s);
        uprac(rho, n, P, 5, 0.618033988749894903, s);
        uprac(rho, n, P, 5, 0.618033988749894903, s);
        uprac(rho, n, P, 7, 0.618033988749894903, s);
        uprac(rho, n, P, 7, 0.618033988749894903, s);
        uprac(rho, n, P, 11, 0.580178728295464130, s);
        uprac(rho, n, P, 11, 0.580178728295464130, s);
        uprac(rho, n, P, 13, 0.618033988749894903, s);
        uprac(rho, n, P, 13, 0.618033988749894903, s);
        uprac(rho, n, P, 17, 0.618033988749894903, s);
        uprac(rho, n, P, 19, 0.618033988749894903, s);
        uprac(rho, n, P, 3427, 0.618033988749894903, s);
        uprac(rho, n, P, 29, 0.548409048446403258, s);
        uprac(rho, n, P, 31, 0.618033988749894903, s);
        uprac(rho, n, P, 4181, 0.618033988749894903, s);
        uprac(rho, n, P, 2173, 0.618033988749894903, s);
        uprac(rho, n, P, 43, 0.618033988749894903, s);
        uprac(rho, n, P, 47, 0.548409048446403258, s);
        uprac(rho, n, P, 8909, 0.580178728295464130, s);
        uprac(rho, n, P, 12017, 0.643969713705029423, s);
        uprac(rho, n, P, 67, 0.580178728295464130, s);
        uprac(rho, n, P, 71, 0.591965645556728037, s);
        uprac(rho, n, P, 73, 0.618033988749894903, s);
        uprac(rho, n, P, 79, 0.618033988749894903, s);
        uprac(rho, n, P, 8051, 0.632839806088706269, s);
        uprac(rho, n, P, 89, 0.618033988749894903, s);
        uprac(rho, n, P, 101, 0.556250337855490828, s);
        uprac(rho, n, P, 103, 0.632839806088706269, s);
        uprac(rho, n, P, 107, 0.580178728295464130, s);
        uprac(rho, n, P, 109, 0.548409048446403258, s);
        uprac(rho, n, P, 17653, 0.586779411332316370, s);
        uprac(rho, n, P, 131, 0.618033988749894903, s);
        uprac(rho, n, P, 22879, 0.591384619013580526, s);
        uprac(rho, n, P, 157, 0.640157392785047019, s);
        uprac(rho, n, P, 163, 0.551390822543526449, s);
        uprac(rho, n, P, 173, 0.612429949509495031, s);
        uprac(rho, n, P, 179, 0.618033988749894903, s);
        uprac(rho, n, P, 181, 0.551390822543526449, s);
        uprac(rho, n, P, 191, 0.618033988749894903, s);
        uprac(rho, n, P, 193, 0.618033988749894903, s);
        uprac(rho, n, P, 199, 0.551390822543526449, s);
    }
}

static inline ulong
gcd64(ulong u, ulong v)
{
    ulong a = u, b = v, c;
    while (b != 0) { c = a % b; a = b; b = c; }
    return a;
}

static inline ulong
crossprodxz(ulong ptx, ulong ptz,
            ulong Pbx, ulong Pbz, ulong Pbprod,
            ulong ptprod, ulong acc, uint rho, ulong n)
{
    ulong tt1 = modsub64(ptx,  Pbx,  n);
    ulong tt2 = modadd64(ptz,  Pbz,  n);
    ulong tt3 = montmul64(tt1, tt2,  n, rho);
    tt1 = modadd64(tt3,  Pbprod, n);
    tt2 = modsub64(tt1,  ptprod, n);
    return montmul64(acc, tt2, n, rho);
}

static inline ulong
uecm_stage2_D30(uecm_pt *P, uint rho, ulong n,
                uint B1, ulong s, ulong unityval)
{
    uecm_pt Pa, Pstep, Pdiff, pt5, pt6;
    ulong result;

    ulong Pbprod[8], PbX[8], PbZ[8];
    ulong diff1, sum1;

    PbX[0] = P->X;  PbZ[0] = P->Z;

    diff1 = modsub64(P->X, P->Z, n);
    sum1  = modadd64(P->X, P->Z, n);
    udup(s, rho, n, sum1, diff1, &Pa);
    PbX[1] = Pa.X;  PbZ[1] = Pa.Z;

    uaddxz(rho, n, PbX[0], PbZ[0], PbX[1], PbZ[1], PbX[0], PbZ[0], &PbX[3], &PbZ[3]);

    diff1 = modsub64(PbX[3], PbZ[3], n);
    sum1  = modadd64(PbX[3], PbZ[3], n);
    udup(s, rho, n, sum1, diff1, &pt6);

    uaddxz(rho, n, PbX[3], PbZ[3], PbX[1], PbZ[1], PbX[0], PbZ[0], &pt5.X, &pt5.Z);
    uaddxz(rho, n, pt6.X, pt6.Z, pt5.X, pt5.Z, PbX[0], PbZ[0], &PbX[2], &PbZ[2]);
    uaddxz(rho, n, PbX[2], PbZ[2], pt6.X, pt6.Z, pt5.X, pt5.Z, &PbX[4], &PbZ[4]);
    uaddxz(rho, n, PbX[4], PbZ[4], pt6.X, pt6.Z, PbX[2], PbZ[2], &PbX[6], &PbZ[6]);
    uaddxz(rho, n, PbX[6], PbZ[6], pt6.X, pt6.Z, PbX[4], PbZ[4], &PbX[7], &PbZ[7]);
    uaddxz(rho, n, pt6.X, pt6.Z, PbX[0], PbZ[0], pt5.X, pt5.Z, &PbX[1], &PbZ[1]);
    uaddxz(rho, n, PbX[1], PbZ[1], pt6.X, pt6.Z, PbX[0], PbZ[0], &PbX[3], &PbZ[3]);
    uaddxz(rho, n, PbX[3], PbZ[3], pt6.X, pt6.Z, PbX[1], PbZ[1], &PbX[5], &PbZ[5]);

    diff1 = modsub64(PbX[0], PbZ[0], n);
    sum1  = modadd64(PbX[0], PbZ[0], n);
    udup(s, rho, n, sum1, diff1, &pt5);
    diff1 = modsub64(pt5.X, pt5.Z, n);
    sum1  = modadd64(pt5.X, pt5.Z, n);
    udup(s, rho, n, sum1, diff1, &pt5);   // pt5=[4]Q

    uaddxz(rho, n, PbX[4], PbZ[4], PbX[3], PbZ[3], pt5.X, pt5.Z, &Pdiff.X, &Pdiff.Z);

    Pstep.X = Pdiff.X;  Pstep.Z = Pdiff.Z;
    diff1 = modsub64(Pstep.X, Pstep.Z, n);
    sum1  = modadd64(Pstep.X, Pstep.Z, n);
    udup(s, rho, n, sum1, diff1, &Pstep);   // [60]Q

    Pdiff.X = Pstep.X;  Pdiff.Z = Pstep.Z;
    Pa.X    = Pstep.X;  Pa.Z    = Pstep.Z;
    diff1 = modsub64(Pa.X, Pa.Z, n);
    sum1  = modadd64(Pa.X, Pa.Z, n);
    udup(s, rho, n, sum1, diff1, &Pa);   // [120]Q

    for (int j = 0; j < 8; j++)
        Pbprod[j] = montmul64(PbX[j], PbZ[j], n, rho);

    uint aval = 120;
    while (aval < B1) {
        ulong tx = Pa.X, tz = Pa.Z;
        uaddxz(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, &Pa.X, &Pa.Z);
        Pdiff.X = tx;  Pdiff.Z = tz;
        aval += 60;
    }

    ulong Paprod = montmul64(Pa.X, Pa.Z, n, rho);
    ulong acc = unityval;
    result = 0;

    while (aval < (25 * B1)) {
        for (int j = 0; j < 8; j++) {
            ulong tmp = crossprodxz(Pa.X, Pa.Z, PbX[j], PbZ[j],
                                    Pbprod[j], Paprod, acc, rho, n);
            if ((tmp == 0) && (result == 0)) result = acc;
            acc = tmp;
        }
        ulong tx = Pa.X, tz = Pa.Z;
        uaddxz(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, &Pa.X, &Pa.Z);
        Pdiff.X = tx;  Pdiff.Z = tz;
        aval += 60;
        Paprod = montmul64(Pa.X, Pa.Z, n, rho);
    }

    return (result == 0) ? acc : result;
}

/* =========================================================================
 * 64-bit ECM kernel
 * =========================================================================
 *
 * CUDA original: gbl_ecm(num, n_in, rho_in, one_in, rsq_in, sigma_in,
 *                         f_out, stg1, curve)
 * OpenCL:        __kernel with get_global_id(0) replacing blockIdx*blockDim+threadIdx
 */
__kernel void gbl_ecm(
    int num,
    __global ulong  *n_in,
    __global uint   *rho_in,
    __global ulong  *one_in,
    __global ulong  *rsq_in,
    __global uint   *sigma_in,
    __global ulong  *f_out,
    uint stg1,
    int  curve)
{
    int idx = (int)get_global_id(0);
    if (idx >= num) return;

    uint  rho      = rho_in[idx];
    ulong n        = n_in[idx];
    ulong unityval = one_in[idx];
    ulong Rsqr     = rsq_in[idx];
    uecm_pt P;
    ulong s;

    f_out[idx] = 1;

    ulong gcd_val;
    gcd_val = uecm_build(&P, rho, n, sigma_in[idx], &s, unityval, Rsqr);


    if (gcd_val > 1 && gcd_val < n) {
        if ((f_out[idx] == 1) || (f_out[idx] == n))
            f_out[idx] = gcd_val;
    }

    uecm_stage1(rho, n, &P, stg1, s);
    gcd_val = gcd64(n, P.Z);
    if ((f_out[idx] == 1) || (f_out[idx] == n))
        f_out[idx] = gcd_val;

    ulong stg2acc = uecm_stage2_D30(&P, rho, n, stg1, s, unityval);
    gcd_val = gcd64(stg2acc, n);
    if ((f_out[idx] == 1) || (f_out[idx] == n))
        f_out[idx] = gcd_val;


    ulong sigma64 = sigma_in[idx];
    sigma_in[idx] = lcg_rand_32b(7, (uint)-1, &sigma64);


}

/* =========================================================================
 * 96-bit ECM point operations
 * ========================================================================= */

static inline void
uaddxz96(uint rho, uint *n,
         uint *p1x, uint *p1z, uint *p2x, uint *p2z,
         uint *pix,  uint *piz,
         uint *pox,  uint *poz)
{
    uint diff1[3], sum1[3], diff2[3], sum2[3];
    uint tt1[3], tt2[3];

    modaddsub96(p1x, p1z, sum1, diff1, n);
    modaddsub96(p2x, p2z, sum2, diff2, n);

    montmul96(diff1, sum2, tt1, n, rho);
    montmul96(sum1,  diff2, tt2, n, rho);

    modaddsub96(tt1, tt2, sum1, diff1, n);

    montsqr96(sum1,  tt1, n, rho);
    montsqr96(diff1, tt2, n, rho);

    montmul96(tt1, piz, tt1, n, rho);
    montmul96(tt2, pix, tt2, n, rho);

    pox[0] = tt1[0]; pox[1] = tt1[1]; pox[2] = tt1[2];
    poz[0] = tt2[0]; poz[1] = tt2[1]; poz[2] = tt2[2];
}

static inline void
uadd96(uint rho, uint *n, uecm96_pt p1, uecm96_pt p2, uecm96_pt pi, uecm96_pt *po)
{
    uaddxz96(rho, n, p1.X, p1.Z, p2.X, p2.Z, pi.X, pi.Z, po->X, po->Z);
}

static inline void
udup96(uint *s, uint rho, uint *n,
       uint *insum, uint *indiff, uecm96_pt *P)
{
    uint tt1[3], tt2[3], tt3[3];
    montsqr96(indiff, tt1, n, rho);
    montsqr96(insum,  tt2, n, rho);
    montmul96(tt1, tt2, P->X, n, rho);
    modsub96(tt2, tt1, tt3, n);
    montmul96(tt3, s, tt2, n, rho);
    modadd96(tt2, tt1, tt2, n);
    montmul96(tt2, tt3, P->Z, n, rho);
}

static inline void
udup96as(uint *s, uint rho, uint *n, uecm96_pt *P)
{
    uint tt1[3], tt2[3], tt3[3], insum[3], indiff[3];
    modsub96(P->X, P->Z, indiff, n);
    modadd96(P->X, P->Z, insum,  n);
    montsqr96(indiff, tt1, n, rho);
    montsqr96(insum,  tt2, n, rho);
    montmul96(tt1, tt2, P->X, n, rho);
    modsub96(tt2, tt1, tt3, n);
    montmul96(tt3, s, tt2, n, rho);
    modadd96(tt2, tt1, tt2, n);
    montmul96(tt2, tt3, P->Z, n, rho);
}

static inline void swap96(uint *a, uint *b)
{
    uint t;
    t=a[0]; a[0]=b[0]; b[0]=t;
    t=a[1]; a[1]=b[1]; b[1]=t;
    t=a[2]; a[2]=b[2]; b[2]=t;
}

static inline void threeswap96(uint *a, uint *b, uint *c)
{
    uint t;
    t=a[0]; a[0]=b[0]; b[0]=c[0]; c[0]=t;
    t=a[1]; a[1]=b[1]; b[1]=c[1]; c[1]=t;
    t=a[2]; a[2]=b[2]; b[2]=c[2]; c[2]=t;
}

static inline void
uprac96(uint rho, uint *n, uecm96_pt *P, ulong c, double v, uint *s)
{
    ulong d, e, r;
    d = c;
    r = (ulong)((double)d * v + 0.5);
    d = c - r;
    e = 2 * r - c;

    uecm96_pt pt1, pt2, pt3, pt4;

    pt1.X[0]=pt2.X[0]=pt3.X[0]=P->X[0];
    pt1.Z[0]=pt2.Z[0]=pt3.Z[0]=P->Z[0];
    pt1.X[1]=pt2.X[1]=pt3.X[1]=P->X[1];
    pt1.Z[1]=pt2.Z[1]=pt3.Z[1]=P->Z[1];
    pt1.X[2]=pt2.X[2]=pt3.X[2]=P->X[2];
    pt1.Z[2]=pt2.Z[2]=pt3.Z[2]=P->Z[2];

    udup96as(s, rho, n, &pt1);

    while (d != e) {
        if (d < e) {
            r = d; d = e; e = r;
            swap96(pt1.X, pt2.X);
            swap96(pt1.Z, pt2.Z);
        }
        if ((d + 3) / 4 <= e) {
            d -= e;
            uadd96(rho, n, pt2, pt1, pt3, &pt4);
            threeswap96(pt2.X, pt4.X, pt3.X);
            threeswap96(pt2.Z, pt4.Z, pt3.Z);
        } else {
            d = (d - e) / 2;
            uadd96(rho, n, pt2, pt1, pt3, &pt2);
            udup96as(s, rho, n, &pt1);
        }
    }
    uadd96(rho, n, pt1, pt2, pt3, P);
}

static inline void
uecm96_build(uecm96_pt *P, uint rho, uint *n,
             uint sigma, uint *ps, uint *one, uint *Rsqr, uint *gcd_out)
{
    uint t1[3], t2[3], t3[3], t4[3], t5[3];
    uint u[3], v[3], x[3], z[3];

    modadd96(one, one, t2, n);
    modadd96(t2,  t2,  t4, n);
    modadd96(one, t4,  t5, n);

    t1[0]=sigma; t1[1]=0; t1[2]=0;
    montmul96(t1, Rsqr, u, n, rho);   // to_monty(sigma)

    modadd96(u, u, v, n);
    modadd96(v, v, v, n);              // 4*sigma

    montmul96(u, u, u, n, rho);
    modsub96(u, t5, u, n);             // sigma^2 - 5

    montsqr96(u, t1, n, rho);
    montmul96(t1, u, x, n, rho);      // u^3

    modadd96(v, v, t2, n);
    modadd96(t2, t2, t2, n);
    modadd96(t2, t2, t2, n);
    modadd96(t2, t2, t2, n);          // 16*v
    montmul96(t2, x, t5, n, rho);    // 16*u^3*v

    montsqr96(v, t1, n, rho);
    montmul96(t1, v, z, n, rho);      // v^3

    modsub96(v, u, t1, n);
    montsqr96(t1, t2, n, rho);
    montmul96(t2, t1, t4, n, rho);    // (v-u)^3

    modadd96(u, u, t1, n);
    modadd96(t1, u, t2, n);           // 3u (note: t2 = 2u+u = 3u)
    modadd96(v, t2, t3, n);           // 3u+v
    montmul96(t3, t4, t1, n, rho);    // a=(v-u)^3*(3u+v)

    t2[0]=1; t2[1]=0; t2[2]=0;
    montmul96(t5, t2, t5, n, rho);    // take t5 out of monty
    modinv96(t5, n, t3, gcd_out);
    montmul96(t3, Rsqr, t3, n, rho);  // to_monty(t3)
    montmul96(t3, t1, ps, n, rho);

    P->X[0]=x[0]; P->X[1]=x[1]; P->X[2]=x[2];
    P->Z[0]=z[0]; P->Z[1]=z[1]; P->Z[2]=z[2];
}

static inline void
uecm96_stage1(uint rho, uint *n, uecm96_pt *P, uint stg1, uint *s)
{
    uint q = 2;
    while (q < stg1) {
        udup96as(s, rho, n, P);
        q *= 2;
    }

    if (stg1 > 600) {
        uprac96(rho, n, P, 601, 0.538760479780959312, s);
        uprac96(rho, n, P, 607, 0.620181980807415711, s);
        uprac96(rho, n, P, 613, 0.554014025731451420, s);
        uprac96(rho, n, P, 617, 0.524531838023777786, s);
        uprac96(rho, n, P, 619, 0.591965645556728037, s);
        uprac96(rho, n, P, 631, 0.634253156188120615, s);
        uprac96(rho, n, P, 641, 0.595297340080313875, s);
        uprac96(rho, n, P, 643, 0.612429949509495031, s);
        uprac96(rho, n, P, 647, 0.612429949509495031, s);
        uprac96(rho, n, P, 653, 0.632839806088706269, s);
        uprac96(rho, n, P, 659, 0.537850274464956368, s);
        uprac96(rho, n, P, 661, 0.553431118763021646, s);
        uprac96(rho, n, P, 673, 0.617214616534403904, s);
        uprac96(rho, n, P, 677, 0.552936068843375761, s);
        uprac96(rho, n, P, 683, 0.551033585330114040, s);
        uprac96(rho, n, P, 691, 0.553431118763021646, s);
    }
    if (stg1 > 500) {
        uprac96(rho, n, P, 503, 0.618033988749894903, s);
        uprac96(rho, n, P, 509, 0.704785919459383070, s);
        uprac96(rho, n, P, 521, 0.553431118763021646, s);
        uprac96(rho, n, P, 523, 0.655268038071952219, s);
        uprac96(rho, n, P, 541, 0.553757286936507276, s);
        uprac96(rho, n, P, 547, 0.618033988749894903, s);
        uprac96(rho, n, P, 557, 0.612429949509495031, s);
        uprac96(rho, n, P, 563, 0.525973719485773428, s);
        uprac96(rho, n, P, 569, 0.580178728295464130, s);
        uprac96(rho, n, P, 571, 0.620181980807415711, s);
        uprac96(rho, n, P, 577, 0.704785919459383070, s);
        uprac96(rho, n, P, 587, 0.518139744387644541, s);
        uprac96(rho, n, P, 593, 0.567795980221515451, s);
        uprac96(rho, n, P, 599, 0.553431118763021646, s);
        uprac96(rho, n, P, 5,   0.618033988749894903, s);
        uprac96(rho, n, P, 23,  0.522786351415446049, s);
    }
    if (stg1 > 400) {
        uprac96(rho, n, P, 184861, 0.628445063812485882, s);
        uprac96(rho, n, P, 409,    0.551390822543526449, s);
        uprac96(rho, n, P, 419,    0.524531838023777786, s);
        uprac96(rho, n, P, 421,    0.643616254685781319, s);
        uprac96(rho, n, P, 431,    0.551390822543526449, s);
        uprac96(rho, n, P, 433,    0.553431118763021646, s);
        uprac96(rho, n, P, 439,    0.618033988749894903, s);
        uprac96(rho, n, P, 217513, 0.586747080461318182, s);
        uprac96(rho, n, P, 449,    0.612429949509495031, s);
        uprac96(rho, n, P, 457,    0.612429949509495031, s);
        uprac96(rho, n, P, 463,    0.553431118763021646, s);
        uprac96(rho, n, P, 467,    0.553431118763021646, s);
        uprac96(rho, n, P, 479,    0.618033988749894903, s);
        uprac96(rho, n, P, 487,    0.591965645556728037, s);
        uprac96(rho, n, P, 499,    0.618033988749894903, s);
    }
    if (stg1 > 300) {
        uprac96(rho, n, P, 307,    0.580178728295464130, s);
        uprac96(rho, n, P, 311,    0.618033988749894903, s);
        uprac96(rho, n, P, 313,    0.618033988749894903, s);
        uprac96(rho, n, P, 317,    0.618033988749894903, s);
        uprac96(rho, n, P, 337,    0.618033988749894903, s);
        uprac96(rho, n, P, 124573, 0.552705982126542983, s);
        uprac96(rho, n, P, 349,    0.632839806088706269, s);
        uprac96(rho, n, P, 95663,  0.618033988749894903, s);
        uprac96(rho, n, P, 142763, 0.524817056539543136, s);
        uprac96(rho, n, P, 373,    0.524531838023777786, s);
        uprac96(rho, n, P, 145157, 0.580178728295464130, s);
        uprac96(rho, n, P, 397,    0.580178728295464130, s);
        uprac96(rho, n, P, 2317,   0.618033988749894903, s);
        uprac96(rho, n, P, 19,     0.618033988749894903, s);
    }
    if (stg1 > 250) {
        uprac96(rho, n, P, 251, 0.541554796058780874, s);
        uprac96(rho, n, P, 263, 0.612429949509495031, s);
        uprac96(rho, n, P, 269, 0.618033988749894903, s);
        if (stg1 <= 300)
            uprac96(rho, n, P, 271, 0.618033988749894903, s);
        uprac96(rho, n, P, 277, 0.618033988749894903, s);
        uprac96(rho, n, P, 281, 0.580178728295464130, s);
        uprac96(rho, n, P, 283, 0.580178728295464130, s);
        uprac96(rho, n, P, 293, 0.551390822543526449, s);
        uprac96(rho, n, P, 17,  0.618033988749894903, s);
    }
    if (stg1 > 200) {
        uprac96(rho, n, P, 48319, 0.552793637425581075, s);
        if (stg1 > 250) {
            uprac96(rho, n, P, 57311, 0.552740754023503311, s);
            uprac96(rho, n, P, 241,   0.625306711365725132, s);
        } else {
            uprac96(rho, n, P, 53743, 0.718522647487825128, s);
        }
        uprac96(rho, n, P, 54253, 0.580178728295464130, s);
        uprac96(rho, n, P, 233,   0.618033988749894903, s);
        uprac96(rho, n, P, 3,     0.618033988749894903, s);
    }

//#define FULLY_UNROLL
    if (stg1 > 175)
    {
        // anything greater than 175 gets B1=200.
        // here we have pair-optimized some primes
#ifdef FULLY_UNROLL
        // this removes hundreds of swaps, and completely removes
        // the prac overhead (all of the if/else stuff).
        // But the crazy thing
        // about GPUs is that it actually runs slightly slower
        // on a H200.  The swaps and switches are apparently invisible given
        // the latency hiding intrinsic to the threadblock architecture,
        // or compiler optimizations, or both. And this is obviously
        // lots more lines of code, which also take effort and cache 
        // to store, load, and decode.
        uecm96_pt p1, p2, p3, p4;

        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p1, p2, p3, P); // 3 = 2 + 1 (1)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p1, p2, p3, P); // 3 = 2 + 1 (1)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p1, p2, p3, P); // 3 = 2 + 1 (1)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p1, p2, p3, P); // 3 = 2 + 1 (1)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, P); // 5 = 2 + 3 (1)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, P); // 5 = 2 + 3 (1)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, P); // 5 = 2 + 3 (1)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, P); // 7 = 2 + 5 (3)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, P); // 7 = 2 + 5 (3)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, P); // 11 = 4 + 7 (3)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, P); // 11 = 4 + 7 (3)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, P); // 13 = 5 + 8 (3)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, P); // 13 = 5 + 8 (3)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p1); // 5 = 2 + 3 (1)
        udup96as(s, rho, n, &p4); // 6 = 2 * 3
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 5 + 6 (1)
        uadd96(rho, n, p4, p3, p1, P); // 17 = 6 + 11 (5)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p3, p4, p1, &p2); // 8 = 5 + 3 (2)
        uadd96(rho, n, p4, p2, p3, &p1); // 11 = 3 + 8 (5)
        uadd96(rho, n, p2, p1, p4, P); // 19 = 8 + 11 (3)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p4); // 21 = 8 + 13 (5)
        uadd96(rho, n, p1, p4, p2, &p3); // 34 = 13 + 21 (8)
        uadd96(rho, n, p4, p3, p1, &p2); // 55 = 21 + 34 (13)
        uadd96(rho, n, p3, p2, p4, &p1); // 89 = 34 + 55 (21)
        uadd96(rho, n, p2, p1, p3, &p4); // 144 = 55 + 89 (34)
        uadd96(rho, n, p1, p4, p2, &p3); // 233 = 89 + 144 (55)
        uadd96(rho, n, p4, p3, p1, &p2); // 377 = 144 + 233 (89)
        uadd96(rho, n, p3, p2, p4, &p1); // 610 = 233 + 377 (144)
        uadd96(rho, n, p2, p1, p3, &p2); // 987 = 377 + 610 (233)
        udup96as(s, rho, n, &p1); // 1220 = 2 * 610
        uadd96(rho, n, p2, p1, p3, &p4); // 2207 = 987 + 1220 (233)
        uadd96(rho, n, p1, p4, p2, P); // 3427 = 1220 + 2207 (987)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 4 + 7 (3)
        uadd96(rho, n, p4, p3, p1, &p2); // 18 = 7 + 11 (4)
        uadd96(rho, n, p3, p2, p4, P); // 29 = 11 + 18 (7)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p2, p3, p4, &p1); // 13 = 8 + 5 (3)
        uadd96(rho, n, p3, p1, p2, &p4); // 18 = 5 + 13 (8)
        uadd96(rho, n, p1, p4, p3, P); // 31 = 13 + 18 (5)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p4); // 21 = 8 + 13 (5)
        uadd96(rho, n, p1, p4, p2, &p3); // 34 = 13 + 21 (8)
        uadd96(rho, n, p4, p3, p1, &p2); // 55 = 21 + 34 (13)
        uadd96(rho, n, p3, p2, p4, &p1); // 89 = 34 + 55 (21)
        uadd96(rho, n, p2, p1, p3, &p4); // 144 = 55 + 89 (34)
        uadd96(rho, n, p1, p4, p2, &p3); // 233 = 89 + 144 (55)
        uadd96(rho, n, p4, p3, p1, &p2); // 377 = 144 + 233 (89)
        uadd96(rho, n, p3, p2, p4, &p1); // 610 = 233 + 377 (144)
        uadd96(rho, n, p2, p1, p3, &p4); // 987 = 377 + 610 (233)
        uadd96(rho, n, p1, p4, p2, &p3); // 1597 = 610 + 987 (377)
        uadd96(rho, n, p4, p3, p1, &p2); // 2584 = 987 + 1597 (610)
        uadd96(rho, n, p3, p2, p4, P); // 4181 = 1597 + 2584 (987)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p4); // 21 = 8 + 13 (5)
        uadd96(rho, n, p1, p4, p2, &p3); // 34 = 13 + 21 (8)
        uadd96(rho, n, p4, p3, p1, &p2); // 55 = 21 + 34 (13)
        uadd96(rho, n, p3, p2, p4, &p1); // 89 = 34 + 55 (21)
        uadd96(rho, n, p2, p1, p3, &p4); // 144 = 55 + 89 (34)
        uadd96(rho, n, p1, p4, p2, &p3); // 233 = 89 + 144 (55)
        uadd96(rho, n, p3, p4, p1, &p2); // 377 = 233 + 144 (89)
        uadd96(rho, n, p4, p2, p3, &p1); // 521 = 144 + 377 (233)
        uadd96(rho, n, p1, p2, p4, &p3); // 898 = 521 + 377 (144)
        uadd96(rho, n, p2, p3, p1, &p4); // 1275 = 377 + 898 (521)
        uadd96(rho, n, p3, p4, p2, P); // 2173 = 898 + 1275 (377)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p3, p4, p1, &p2); // 8 = 5 + 3 (2)
        uadd96(rho, n, p4, p2, p3, &p4); // 11 = 3 + 8 (5)
        udup96as(s, rho, n, &p2); // 16 = 2 * 8
        uadd96(rho, n, p4, p2, p3, &p1); // 27 = 11 + 16 (5)
        uadd96(rho, n, p2, p1, p4, P); // 43 = 16 + 27 (11)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 4 + 7 (3)
        uadd96(rho, n, p4, p3, p1, &p2); // 18 = 7 + 11 (4)
        uadd96(rho, n, p3, p2, p4, &p1); // 29 = 11 + 18 (7)
        uadd96(rho, n, p2, p1, p3, P); // 47 = 18 + 29 (11)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, &p2); // 7 = 2 + 5 (3)
        uadd96(rho, n, p3, p2, p1, &p4); // 12 = 5 + 7 (2)
        uadd96(rho, n, p2, p4, p3, &p1); // 19 = 7 + 12 (5)
        uadd96(rho, n, p4, p1, p2, &p3); // 31 = 12 + 19 (7)
        uadd96(rho, n, p1, p3, p4, &p2); // 50 = 19 + 31 (12)
        uadd96(rho, n, p3, p2, p1, &p4); // 81 = 31 + 50 (19)
        uadd96(rho, n, p2, p4, p3, &p1); // 131 = 50 + 81 (31)
        uadd96(rho, n, p1, p4, p2, &p3); // 212 = 131 + 81 (50)
        uadd96(rho, n, p4, p3, p1, &p2); // 293 = 81 + 212 (131)
        uadd96(rho, n, p3, p2, p4, &p1); // 505 = 212 + 293 (81)
        uadd96(rho, n, p2, p1, p3, &p4); // 798 = 293 + 505 (212)
        uadd96(rho, n, p1, p4, p2, &p3); // 1303 = 505 + 798 (293)
        uadd96(rho, n, p4, p3, p1, &p2); // 2101 = 798 + 1303 (505)
        uadd96(rho, n, p3, p2, p4, &p1); // 3404 = 1303 + 2101 (798)
        uadd96(rho, n, p2, p1, p3, &p4); // 5505 = 2101 + 3404 (1303)
        uadd96(rho, n, p1, p4, p2, P); // 8909 = 3404 + 5505 (2101)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p1); // 5 = 2 + 3 (1)
        udup96as(s, rho, n, &p4); // 6 = 2 * 3
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 5 + 6 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 17 = 6 + 11 (5)
        uadd96(rho, n, p3, p2, p4, &p1); // 28 = 11 + 17 (6)
        uadd96(rho, n, p2, p1, p3, &p4); // 45 = 17 + 28 (11)
        uadd96(rho, n, p1, p4, p2, &p3); // 73 = 28 + 45 (17)
        uadd96(rho, n, p4, p3, p1, &p2); // 118 = 45 + 73 (28)
        uadd96(rho, n, p3, p2, p4, &p1); // 191 = 73 + 118 (45)
        uadd96(rho, n, p2, p1, p3, &p4); // 309 = 118 + 191 (73)
        uadd96(rho, n, p1, p4, p2, &p3); // 500 = 191 + 309 (118)
        uadd96(rho, n, p4, p3, p1, &p2); // 809 = 309 + 500 (191)
        uadd96(rho, n, p3, p2, p4, &p3); // 1309 = 500 + 809 (309)
        udup96as(s, rho, n, &p2); // 1618 = 2 * 809
        uadd96(rho, n, p3, p2, p4, &p1); // 2927 = 1309 + 1618 (309)
        uadd96(rho, n, p2, p1, p3, &p4); // 4545 = 1618 + 2927 (1309)
        uadd96(rho, n, p1, p4, p2, &p3); // 7472 = 2927 + 4545 (1618)
        uadd96(rho, n, p4, p3, p1, P); // 12017 = 4545 + 7472 (2927)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, &p2); // 7 = 2 + 5 (3)
        uadd96(rho, n, p3, p2, p1, &p4); // 12 = 5 + 7 (2)
        uadd96(rho, n, p2, p4, p3, &p2); // 19 = 7 + 12 (5)
        udup96as(s, rho, n, &p4); // 24 = 2 * 12
        uadd96(rho, n, p2, p4, p3, &p1); // 43 = 19 + 24 (5)
        uadd96(rho, n, p4, p1, p2, P); // 67 = 24 + 43 (19)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, &p1); // 7 = 2 + 5 (3)
        udup96as(s, rho, n, &p3); // 10 = 2 * 5
        uadd96(rho, n, p1, p3, p4, &p2); // 17 = 7 + 10 (3)
        uadd96(rho, n, p3, p2, p1, &p4); // 27 = 10 + 17 (7)
        uadd96(rho, n, p2, p4, p3, &p1); // 44 = 17 + 27 (10)
        uadd96(rho, n, p4, p1, p2, P); // 71 = 27 + 44 (17)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p2); // 21 = 8 + 13 (5)
        udup96as(s, rho, n, &p1); // 26 = 2 * 13
        uadd96(rho, n, p2, p1, p3, &p4); // 47 = 21 + 26 (5)
        uadd96(rho, n, p1, p4, p2, P); // 73 = 26 + 47 (21)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p1, p2, p3, &p4); // 21 = 13 + 8 (5)
        uadd96(rho, n, p2, p4, p1, &p3); // 29 = 8 + 21 (13)
        uadd96(rho, n, p4, p3, p2, &p1); // 50 = 21 + 29 (8)
        uadd96(rho, n, p3, p1, p4, P); // 79 = 29 + 50 (21)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p3, p4, p1, &p2); // 8 = 5 + 3 (2)
        uadd96(rho, n, p4, p2, p3, &p1); // 11 = 3 + 8 (5)
        uadd96(rho, n, p2, p1, p4, &p3); // 19 = 8 + 11 (3)
        uadd96(rho, n, p1, p3, p2, &p4); // 30 = 11 + 19 (8)
        uadd96(rho, n, p3, p4, p1, &p2); // 49 = 19 + 30 (11)
        uadd96(rho, n, p4, p2, p3, &p1); // 79 = 30 + 49 (19)
        uadd96(rho, n, p2, p1, p4, &p3); // 128 = 49 + 79 (30)
        uadd96(rho, n, p1, p3, p2, &p4); // 207 = 79 + 128 (49)
        uadd96(rho, n, p3, p4, p1, &p2); // 335 = 128 + 207 (79)
        uadd96(rho, n, p4, p2, p3, &p1); // 542 = 207 + 335 (128)
        uadd96(rho, n, p2, p1, p4, &p2); // 877 = 335 + 542 (207)
        udup96as(s, rho, n, &p1); // 1084 = 2 * 542
        uadd96(rho, n, p2, p1, p4, &p3); // 1961 = 877 + 1084 (207)
        uadd96(rho, n, p1, p3, p2, &p4); // 3045 = 1084 + 1961 (877)
        uadd96(rho, n, p3, p4, p1, &p2); // 5006 = 1961 + 3045 (1084)
        uadd96(rho, n, p4, p2, p3, P); // 8051 = 3045 + 5006 (1961)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p4); // 21 = 8 + 13 (5)
        uadd96(rho, n, p1, p4, p2, &p3); // 34 = 13 + 21 (8)
        uadd96(rho, n, p4, p3, p1, &p2); // 55 = 21 + 34 (13)
        uadd96(rho, n, p3, p2, p4, P); // 89 = 34 + 55 (21)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 4 + 7 (3)
        uadd96(rho, n, p4, p3, p1, &p2); // 18 = 7 + 11 (4)
        uadd96(rho, n, p3, p2, p4, &p3); // 29 = 11 + 18 (7)
        udup96as(s, rho, n, &p2); // 36 = 2 * 18
        uadd96(rho, n, p3, p2, p4, &p1); // 65 = 29 + 36 (7)
        uadd96(rho, n, p2, p1, p3, P); // 101 = 36 + 65 (29)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p3, p4, p1, &p2); // 8 = 5 + 3 (2)
        uadd96(rho, n, p4, p2, p3, &p1); // 11 = 3 + 8 (5)
        uadd96(rho, n, p1, p2, p4, &p3); // 19 = 11 + 8 (3)
        uadd96(rho, n, p2, p3, p1, &p2); // 27 = 8 + 19 (11)
        udup96as(s, rho, n, &p3); // 38 = 2 * 19
        uadd96(rho, n, p2, p3, p1, &p4); // 65 = 27 + 38 (11)
        uadd96(rho, n, p3, p4, p2, P); // 103 = 38 + 65 (27)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, &p2); // 7 = 2 + 5 (3)
        uadd96(rho, n, p3, p2, p1, &p4); // 12 = 5 + 7 (2)
        uadd96(rho, n, p2, p4, p3, &p1); // 19 = 7 + 12 (5)
        uadd96(rho, n, p4, p1, p2, &p4); // 31 = 12 + 19 (7)
        udup96as(s, rho, n, &p1); // 38 = 2 * 19
        uadd96(rho, n, p4, p1, p2, &p3); // 69 = 31 + 38 (7)
        uadd96(rho, n, p1, p3, p4, P); // 107 = 38 + 69 (31)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 4 + 7 (3)
        uadd96(rho, n, p4, p3, p1, &p2); // 18 = 7 + 11 (4)
        uadd96(rho, n, p2, p3, p4, &p1); // 29 = 18 + 11 (7)
        uadd96(rho, n, p3, p1, p2, &p4); // 40 = 11 + 29 (18)
        uadd96(rho, n, p1, p4, p3, &p2); // 69 = 29 + 40 (11)
        uadd96(rho, n, p4, p2, p1, P); // 109 = 40 + 69 (29)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, &p2); // 7 = 2 + 5 (3)
        uadd96(rho, n, p2, p3, p1, &p4); // 12 = 7 + 5 (2)
        uadd96(rho, n, p3, p4, p2, &p1); // 17 = 5 + 12 (7)
        uadd96(rho, n, p4, p1, p3, &p2); // 29 = 12 + 17 (5)
        uadd96(rho, n, p1, p2, p4, &p3); // 46 = 17 + 29 (12)
        uadd96(rho, n, p2, p3, p1, &p4); // 75 = 29 + 46 (17)
        uadd96(rho, n, p3, p4, p2, &p1); // 121 = 46 + 75 (29)
        uadd96(rho, n, p4, p1, p3, &p2); // 196 = 75 + 121 (46)
        uadd96(rho, n, p2, p1, p4, &p3); // 317 = 196 + 121 (75)
        uadd96(rho, n, p1, p3, p2, &p4); // 438 = 121 + 317 (196)
        uadd96(rho, n, p3, p4, p1, &p2); // 755 = 317 + 438 (121)
        uadd96(rho, n, p4, p2, p3, &p1); // 1193 = 438 + 755 (317)
        uadd96(rho, n, p2, p1, p4, &p3); // 1948 = 755 + 1193 (438)
        uadd96(rho, n, p1, p3, p2, &p4); // 3141 = 1193 + 1948 (755)
        uadd96(rho, n, p3, p4, p1, &p3); // 5089 = 1948 + 3141 (1193)
        udup96as(s, rho, n, &p4); // 6282 = 2 * 3141
        uadd96(rho, n, p3, p4, p1, &p2); // 11371 = 5089 + 6282 (1193)
        uadd96(rho, n, p4, p2, p3, P); // 17653 = 6282 + 11371 (5089)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p4); // 21 = 8 + 13 (5)
        uadd96(rho, n, p1, p4, p2, &p3); // 34 = 13 + 21 (8)
        uadd96(rho, n, p3, p4, p1, &p2); // 55 = 34 + 21 (13)
        uadd96(rho, n, p4, p2, p3, &p1); // 76 = 21 + 55 (34)
        uadd96(rho, n, p2, p1, p4, P); // 131 = 55 + 76 (21)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p4, p1, p2, &p3); // 5 = 3 + 2 (1)
        uadd96(rho, n, p1, p3, p4, &p1); // 7 = 2 + 5 (3)
        udup96as(s, rho, n, &p3); // 10 = 2 * 5
        uadd96(rho, n, p1, p3, p4, &p2); // 17 = 7 + 10 (3)
        uadd96(rho, n, p3, p2, p1, &p4); // 27 = 10 + 17 (7)
        uadd96(rho, n, p2, p4, p3, &p1); // 44 = 17 + 27 (10)
        uadd96(rho, n, p4, p1, p2, &p3); // 71 = 27 + 44 (17)
        uadd96(rho, n, p1, p3, p4, &p2); // 115 = 44 + 71 (27)
        uadd96(rho, n, p3, p2, p1, &p4); // 186 = 71 + 115 (44)
        uadd96(rho, n, p2, p4, p3, &p1); // 301 = 115 + 186 (71)
        uadd96(rho, n, p4, p1, p2, &p3); // 487 = 186 + 301 (115)
        uadd96(rho, n, p1, p3, p4, &p2); // 788 = 301 + 487 (186)
        uadd96(rho, n, p3, p2, p1, &p4); // 1275 = 487 + 788 (301)
        uadd96(rho, n, p2, p4, p3, &p1); // 2063 = 788 + 1275 (487)
        uadd96(rho, n, p4, p1, p2, &p3); // 3338 = 1275 + 2063 (788)
        uadd96(rho, n, p1, p3, p4, &p2); // 5401 = 2063 + 3338 (1275)
        uadd96(rho, n, p3, p2, p1, &p4); // 8739 = 3338 + 5401 (2063)
        uadd96(rho, n, p2, p4, p3, &p1); // 14140 = 5401 + 8739 (3338)
        uadd96(rho, n, p4, p1, p2, P); // 22879 = 8739 + 14140 (5401)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p1); // 5 = 2 + 3 (1)
        udup96as(s, rho, n, &p4); // 6 = 2 * 3
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 5 + 6 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 17 = 6 + 11 (5)
        uadd96(rho, n, p3, p2, p4, &p1); // 28 = 11 + 17 (6)
        uadd96(rho, n, p2, p1, p3, &p2); // 45 = 17 + 28 (11)
        udup96as(s, rho, n, &p1); // 56 = 2 * 28
        uadd96(rho, n, p2, p1, p3, &p4); // 101 = 45 + 56 (11)
        uadd96(rho, n, p1, p4, p2, P); // 157 = 56 + 101 (45)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 4 + 7 (3)
        uadd96(rho, n, p4, p3, p1, &p2); // 18 = 7 + 11 (4)
        uadd96(rho, n, p3, p2, p4, &p1); // 29 = 11 + 18 (7)
        uadd96(rho, n, p2, p1, p3, &p2); // 47 = 18 + 29 (11)
        udup96as(s, rho, n, &p1); // 58 = 2 * 29
        uadd96(rho, n, p2, p1, p3, &p4); // 105 = 47 + 58 (11)
        uadd96(rho, n, p1, p4, p2, P); // 163 = 58 + 105 (47)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p2, p3, p4, &p1); // 13 = 8 + 5 (3)
        uadd96(rho, n, p3, p1, p2, &p4); // 18 = 5 + 13 (8)
        uadd96(rho, n, p1, p4, p3, &p2); // 31 = 13 + 18 (5)
        uadd96(rho, n, p4, p2, p1, &p4); // 49 = 18 + 31 (13)
        udup96as(s, rho, n, &p2); // 62 = 2 * 31
        uadd96(rho, n, p4, p2, p1, &p3); // 111 = 49 + 62 (13)
        uadd96(rho, n, p2, p3, p4, P); // 173 = 62 + 111 (49)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p1, p2, p3, &p4); // 21 = 13 + 8 (5)
        uadd96(rho, n, p2, p4, p1, &p3); // 29 = 8 + 21 (13)
        uadd96(rho, n, p4, p3, p2, &p1); // 50 = 21 + 29 (8)
        uadd96(rho, n, p3, p1, p4, &p2); // 79 = 29 + 50 (21)
        uadd96(rho, n, p2, p1, p3, &p4); // 129 = 79 + 50 (29)
        uadd96(rho, n, p1, p4, p2, P); // 179 = 50 + 129 (79)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 4 + 7 (3)
        uadd96(rho, n, p4, p3, p1, &p2); // 18 = 7 + 11 (4)
        uadd96(rho, n, p3, p2, p4, &p1); // 29 = 11 + 18 (7)
        uadd96(rho, n, p2, p1, p3, &p4); // 47 = 18 + 29 (11)
        uadd96(rho, n, p4, p1, p2, &p3); // 76 = 47 + 29 (18)
        uadd96(rho, n, p1, p3, p4, &p2); // 105 = 29 + 76 (47)
        uadd96(rho, n, p3, p2, p1, P); // 181 = 76 + 105 (29)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p4); // 21 = 8 + 13 (5)
        uadd96(rho, n, p1, p4, p2, &p3); // 34 = 13 + 21 (8)
        uadd96(rho, n, p4, p3, p1, &p4); // 55 = 21 + 34 (13)
        udup96as(s, rho, n, &p3); // 68 = 2 * 34
        uadd96(rho, n, p4, p3, p1, &p2); // 123 = 55 + 68 (13)
        uadd96(rho, n, p3, p2, p4, P); // 191 = 68 + 123 (55)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p4); // 3 = 1 + 2 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 5 = 2 + 3 (1)
        uadd96(rho, n, p4, p3, p1, &p2); // 8 = 3 + 5 (2)
        uadd96(rho, n, p3, p2, p4, &p1); // 13 = 5 + 8 (3)
        uadd96(rho, n, p2, p1, p3, &p2); // 21 = 8 + 13 (5)
        udup96as(s, rho, n, &p1); // 26 = 2 * 13
        uadd96(rho, n, p2, p1, p3, &p4); // 47 = 21 + 26 (5)
        uadd96(rho, n, p1, p4, p2, &p3); // 73 = 26 + 47 (21)
        uadd96(rho, n, p4, p3, p1, &p2); // 120 = 47 + 73 (26)
        uadd96(rho, n, p3, p2, p4, P); // 193 = 73 + 120 (47)
        p1.X[0] = p2.X[0] = p3.X[0] = P->X[0];
        p1.Z[0] = p2.Z[0] = p3.Z[0] = P->Z[0];
        p1.X[1] = p2.X[1] = p3.X[1] = P->X[1];
        p1.Z[1] = p2.Z[1] = p3.Z[1] = P->Z[1];
        p1.X[2] = p2.X[2] = p3.X[2] = P->X[2];
        p1.Z[2] = p2.Z[2] = p3.Z[2] = P->Z[2];
        udup96as(s, rho, n, &p1); // 2 = 2 * 1
        uadd96(rho, n, p2, p1, p3, &p2); // 3 = 1 + 2 (1)
        udup96as(s, rho, n, &p1); // 4 = 2 * 2
        uadd96(rho, n, p2, p1, p3, &p4); // 7 = 3 + 4 (1)
        uadd96(rho, n, p1, p4, p2, &p3); // 11 = 4 + 7 (3)
        uadd96(rho, n, p4, p3, p1, &p2); // 18 = 7 + 11 (4)
        uadd96(rho, n, p3, p2, p4, &p1); // 29 = 11 + 18 (7)
        uadd96(rho, n, p2, p1, p3, &p4); // 47 = 18 + 29 (11)
        uadd96(rho, n, p1, p4, p2, &p3); // 76 = 29 + 47 (18)
        uadd96(rho, n, p4, p3, p1, &p2); // 123 = 47 + 76 (29)
        uadd96(rho, n, p3, p2, p4, P); // 199 = 76 + 123 (47)

#else

        uprac96(rho, n, P, 3, 0.618033988749894903, s);
        uprac96(rho, n, P, 3, 0.618033988749894903, s);
        uprac96(rho, n, P, 3, 0.618033988749894903, s);
        uprac96(rho, n, P, 3, 0.618033988749894903, s);
        uprac96(rho, n, P, 5, 0.618033988749894903, s);
        uprac96(rho, n, P, 5, 0.618033988749894903, s);
        uprac96(rho, n, P, 5, 0.618033988749894903, s);
        uprac96(rho, n, P, 7, 0.618033988749894903, s);
        uprac96(rho, n, P, 7, 0.618033988749894903, s);
        uprac96(rho, n, P, 11, 0.580178728295464130, s);
        uprac96(rho, n, P, 11, 0.580178728295464130, s);
        uprac96(rho, n, P, 13, 0.618033988749894903, s);
        uprac96(rho, n, P, 13, 0.618033988749894903, s);
        uprac96(rho, n, P, 17, 0.618033988749894903, s);
        uprac96(rho, n, P, 19, 0.618033988749894903, s);
        uprac96(rho, n, P, 3427, 0.618033988749894903, s); //[23,149], saving 3
        uprac96(rho, n, P, 29, 0.548409048446403258, s);
        uprac96(rho, n, P, 31, 0.618033988749894903, s);
        uprac96(rho, n, P, 4181, 0.618033988749894903, s); //[37,113], saving 3
        uprac96(rho, n, P, 2173, 0.618033988749894903, s); //[41,53], saving 3
        uprac96(rho, n, P, 43, 0.618033988749894903, s);
        uprac96(rho, n, P, 47, 0.548409048446403258, s);
        uprac96(rho, n, P, 8909, 0.580178728295464130, s); //[59,151], saving 2
        uprac96(rho, n, P, 12017, 0.643969713705029423, s); //[61,197], saving 3
        uprac96(rho, n, P, 67, 0.580178728295464130, s);
        uprac96(rho, n, P, 71, 0.591965645556728037, s);
        uprac96(rho, n, P, 73, 0.618033988749894903, s);
        uprac96(rho, n, P, 79, 0.618033988749894903, s);
        uprac96(rho, n, P, 8051, 0.632839806088706269, s); //[83,97], saving 3
        uprac96(rho, n, P, 89, 0.618033988749894903, s);
        uprac96(rho, n, P, 101, 0.556250337855490828, s);
        uprac96(rho, n, P, 103, 0.632839806088706269, s);
        uprac96(rho, n, P, 107, 0.580178728295464130, s);
        uprac96(rho, n, P, 109, 0.548409048446403258, s);
        uprac96(rho, n, P, 17653, 0.586779411332316370, s); //[127,139], saving 3
        uprac96(rho, n, P, 131, 0.618033988749894903, s);
        uprac96(rho, n, P, 22879, 0.591384619013580526, s); //[137,167], saving 3
        uprac96(rho, n, P, 157, 0.640157392785047019, s);
        uprac96(rho, n, P, 163, 0.551390822543526449, s);
        uprac96(rho, n, P, 173, 0.612429949509495031, s); //[173,347], savings = 2
        uprac96(rho, n, P, 179, 0.618033988749894903, s); //[179,379], savings = 1, [179,461], savings = 3
        uprac96(rho, n, P, 181, 0.551390822543526449, s);
        uprac96(rho, n, P, 191, 0.618033988749894903, s);
        uprac96(rho, n, P, 193, 0.618033988749894903, s);
        uprac96(rho, n, P, 199, 0.551390822543526449, s);
#endif

    }

    return;
}

/* -------------------------------------------------------------------------
 * gcd96: 96-bit GCD using Euclidean algorithm via uint128 emulation
 * (replaces __uint128_t which is unavailable in OpenCL)
 * ------------------------------------------------------------------------- */
static inline void
gcd96(uint *u, uint *v, uint *gcd_out)
{
    uint128 a = u128_from_limbs3(u[0], u[1], u[2]);
    uint128 b = u128_from_limbs3(v[0], v[1], v[2]);
    uint128 c;
    uint128 zero = {0UL, 0UL};

    while (!u128_eq(b, zero)) {
        u128_divmod(a, b, &c);   /* c = remainder of a / b */
        a = b;
        b = c;
    }

    gcd_out[0] = (uint)a.lo;
    gcd_out[1] = (uint)(a.lo >> 32);
    gcd_out[2] = (uint)a.hi;
}

/* -------------------------------------------------------------------------
 * Stage 2 D=30
 * ------------------------------------------------------------------------- */
#ifdef USE_D30
static inline void
uecm96_stage2_D30(uecm96_pt *P, uint rho, uint *n,
                  uint B1, uint B2, uint *s, uint *unityval, uint *result)
{
    uecm96_pt Pa, Pstep, Pdiff, pt5, pt6;
    uint Pbprod[8][3], PbX[8][3], PbZ[8][3];
    uint diff1[3], sum1[3];

    PbX[0][0]=P->X[0]; PbX[0][1]=P->X[1]; PbX[0][2]=P->X[2];
    PbZ[0][0]=P->Z[0]; PbZ[0][1]=P->Z[1]; PbZ[0][2]=P->Z[2];

    modsub96(P->X, P->Z, diff1, n);
    modadd96(P->X, P->Z, sum1,  n);
    udup96(s, rho, n, sum1, diff1, &Pa);
    PbX[1][0]=Pa.X[0]; PbX[1][1]=Pa.X[1]; PbX[1][2]=Pa.X[2];
    PbZ[1][0]=Pa.Z[0]; PbZ[1][1]=Pa.Z[1]; PbZ[1][2]=Pa.Z[2];

    uaddxz96(rho, n, PbX[0], PbZ[0], PbX[1], PbZ[1], PbX[0], PbZ[0], PbX[3], PbZ[3]);

    modsub96(PbX[3], PbZ[3], diff1, n);
    modadd96(PbX[3], PbZ[3], sum1,  n);
    udup96(s, rho, n, sum1, diff1, &pt6);

    uaddxz96(rho, n, PbX[3], PbZ[3], PbX[1], PbZ[1], PbX[0], PbZ[0], pt5.X, pt5.Z);
    uaddxz96(rho, n, pt6.X, pt6.Z, pt5.X, pt5.Z, PbX[0], PbZ[0], PbX[2], PbZ[2]);
    uaddxz96(rho, n, PbX[2], PbZ[2], pt6.X, pt6.Z, pt5.X, pt5.Z, PbX[4], PbZ[4]);
    uaddxz96(rho, n, PbX[4], PbZ[4], pt6.X, pt6.Z, PbX[2], PbZ[2], PbX[6], PbZ[6]);
    uaddxz96(rho, n, PbX[6], PbZ[6], pt6.X, pt6.Z, PbX[4], PbZ[4], PbX[7], PbZ[7]);
    uaddxz96(rho, n, pt6.X, pt6.Z, PbX[0], PbZ[0], pt5.X, pt5.Z, PbX[1], PbZ[1]);
    uaddxz96(rho, n, PbX[1], PbZ[1], pt6.X, pt6.Z, PbX[0], PbZ[0], PbX[3], PbZ[3]);
    uaddxz96(rho, n, PbX[3], PbZ[3], pt6.X, pt6.Z, PbX[1], PbZ[1], PbX[5], PbZ[5]);

    modsub96(PbX[0], PbZ[0], diff1, n);
    modadd96(PbX[0], PbZ[0], sum1,  n);
    udup96(s, rho, n, sum1, diff1, &pt5);
    modsub96(pt5.X, pt5.Z, diff1, n);
    modadd96(pt5.X, pt5.Z, sum1,  n);
    udup96(s, rho, n, sum1, diff1, &pt5);  // [4]Q

    uaddxz96(rho, n, PbX[4], PbZ[4], PbX[3], PbZ[3], pt5.X, pt5.Z, Pdiff.X, Pdiff.Z);

    Pstep.X[0]=Pdiff.X[0]; Pstep.X[1]=Pdiff.X[1]; Pstep.X[2]=Pdiff.X[2];
    Pstep.Z[0]=Pdiff.Z[0]; Pstep.Z[1]=Pdiff.Z[1]; Pstep.Z[2]=Pdiff.Z[2];
    modsub96(Pstep.X, Pstep.Z, diff1, n);
    modadd96(Pstep.X, Pstep.Z, sum1,  n);
    udup96(s, rho, n, sum1, diff1, &Pstep);  // [60]Q

    Pdiff.X[0]=Pstep.X[0]; Pdiff.X[1]=Pstep.X[1]; Pdiff.X[2]=Pstep.X[2];
    Pdiff.Z[0]=Pstep.Z[0]; Pdiff.Z[1]=Pstep.Z[1]; Pdiff.Z[2]=Pstep.Z[2];
    Pa.X[0]=Pdiff.X[0]; Pa.X[1]=Pdiff.X[1]; Pa.X[2]=Pdiff.X[2];
    Pa.Z[0]=Pdiff.Z[0]; Pa.Z[1]=Pdiff.Z[1]; Pa.Z[2]=Pdiff.Z[2];
    modsub96(Pa.X, Pa.Z, diff1, n);
    modadd96(Pa.X, Pa.Z, sum1,  n);
    udup96(s, rho, n, sum1, diff1, &Pa);  // [120]Q

    for (int j = 0; j < 8; j++)
        montmul96(PbX[j], PbZ[j], Pbprod[j], n, rho);

    uint aval = 120;
    while (aval < B1) {
        pt5.X[0]=Pa.X[0]; pt5.X[1]=Pa.X[1]; pt5.X[2]=Pa.X[2];
        pt5.Z[0]=Pa.Z[0]; pt5.Z[1]=Pa.Z[1]; pt5.Z[2]=Pa.Z[2];
        uaddxz96(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, Pa.X, Pa.Z);
        Pdiff.X[0]=pt5.X[0]; Pdiff.X[1]=pt5.X[1]; Pdiff.X[2]=pt5.X[2];
        Pdiff.Z[0]=pt5.Z[0]; Pdiff.Z[1]=pt5.Z[1]; Pdiff.Z[2]=pt5.Z[2];
        aval += 60;
    }

    uint Paprod[3];
    montmul96(Pa.X, Pa.Z, Paprod, n, rho);

    uint acc[3];
    acc[0]=unityval[0]; acc[1]=unityval[1]; acc[2]=unityval[2];

    while (aval < B2) {
        uint tmp[3];
        for (int j = 0; j < 8; j++) {
            modsub96(Pa.X, PbX[j], diff1, n);
            modadd96(Pa.Z, PbZ[j], sum1,  n);
            montmul96(diff1, sum1, tmp, n, rho);
            modadd96(tmp, Pbprod[j], sum1, n);
            modsub96(sum1, Paprod, diff1, n);
            montmul96(acc, diff1, acc, n, rho);
        }
        pt5.X[0]=Pa.X[0]; pt5.X[1]=Pa.X[1]; pt5.X[2]=Pa.X[2];
        pt5.Z[0]=Pa.Z[0]; pt5.Z[1]=Pa.Z[1]; pt5.Z[2]=Pa.Z[2];
        uaddxz96(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, Pa.X, Pa.Z);
        Pdiff.X[0]=pt5.X[0]; Pdiff.X[1]=pt5.X[1]; Pdiff.X[2]=pt5.X[2];
        Pdiff.Z[0]=pt5.Z[0]; Pdiff.Z[1]=pt5.Z[1]; Pdiff.Z[2]=pt5.Z[2];
        aval += 60;
        montmul96(Pa.X, Pa.Z, Paprod, n, rho);
    }

    result[0]=acc[0]; result[1]=acc[1]; result[2]=acc[2];
}
#endif /* USE_D30 */

/* =========================================================================
 * 96-bit ECM kernel
 *
 * CUDA original: gbl_ecm96(num, n_in, rho_in, one_in, rsq_in, sigma_in,
 *                            f_out, stg1, stg2, curve)
 * ========================================================================= */
__kernel void gbl_ecm96(
    int num,
    __global uint *n_in,
    __global uint *rho_in,
    __global uint *one_in,
    __global uint *rsq_in,
    __global uint *sigma_in,
    __global uint *f_out,
    uint stg1,
    uint stg2,
    int  curve)
{
    int idx = (int)get_global_id(0);
    if (idx >= num) return;

    uint rho = rho_in[idx];
    uint n[3], unityval[3], Rsqr[3];
    uecm96_pt P;
    uint s[3];

    n[0]        = n_in[idx*3+0];       n[1]        = n_in[idx*3+1];       n[2]        = n_in[idx*3+2];
    unityval[0] = one_in[idx*3+0];     unityval[1] = one_in[idx*3+1];     unityval[2] = one_in[idx*3+2];
    Rsqr[0]     = rsq_in[idx*3+0];     Rsqr[1]     = rsq_in[idx*3+1];     Rsqr[2]     = rsq_in[idx*3+2];

    f_out[3*idx+0] = 1;
    f_out[3*idx+1] = 0;
    f_out[3*idx+2] = 0;

    uint gcd_val[3];
    uecm96_build(&P, rho, n, sigma_in[idx], s, unityval, Rsqr, gcd_val);

    if (gcd_val[0] > 1) {
        if ((gcd_val[0] != n[0]) || (gcd_val[1] != n[1]) || (gcd_val[2] != n[2])) {
            if (((f_out[3*idx+0]==1) && (f_out[3*idx+1]==0) && (f_out[3*idx+2]==0)) ||
                ((f_out[3*idx+0]==n[0]) && (f_out[3*idx+1]==n[1]) && (f_out[3*idx+2]==n[2])))
            {
                f_out[3*idx+0] = gcd_val[0];
                f_out[3*idx+1] = gcd_val[1];
                f_out[3*idx+2] = gcd_val[2];
            }
        }
    }

    uecm96_stage1(rho, n, &P, stg1, s);
    gcd96(n, P.Z, gcd_val);

    if (((f_out[3*idx+0]==1) && (f_out[3*idx+1]==0) && (f_out[3*idx+2]==0)) ||
        ((f_out[3*idx+0]==n[0]) && (f_out[3*idx+1]==n[1]) && (f_out[3*idx+2]==n[2])))
    {
        f_out[3*idx+0] = gcd_val[0];
        f_out[3*idx+1] = gcd_val[1];
        f_out[3*idx+2] = gcd_val[2];
    }

    if (1) {
        uint stg2acc[3];
#ifdef USE_D30
        uecm96_stage2_D30(&P, rho, n, stg1, stg2, s, unityval, stg2acc);
#else
        stg2acc[0] = 1; stg2acc[1] = 0; stg2acc[2] = 0;
#endif

        gcd96(n, stg2acc, gcd_val);
        if (((f_out[3 * idx + 0] == 1) && (f_out[3 * idx + 1] == 0) && (f_out[3 * idx + 2] == 0)) ||
            ((f_out[3 * idx + 0] == n[0]) && (f_out[3 * idx + 1] == n[1]) && (f_out[3 * idx + 2] == n[2])))
        {
            f_out[3 * idx + 0] = gcd_val[0];
            f_out[3 * idx + 1] = gcd_val[1];
            f_out[3 * idx + 2] = gcd_val[2];
        }
    }

    ulong sigma64 = sigma_in[idx];
    sigma_in[idx] = lcg_rand_32b(7, (uint)-1, &sigma64);
}

/* =========================================================================
 * P-1 kernel support (pm196)
 * ========================================================================= */

/* __constant__ tables -- same data as original, same keyword in OpenCL */
__constant int tpm1_ewin100[34] = {
    12, 12, 14, 3, 12, 7, 13, 6, 12, 0, 15, 4, 1, 8, 7, 3, 0, 14, 13, 6,
    5, 6, 15, 13, 0, 11, 0, 14, 13, 3, 8, 8, 12, 0 };

__constant int tpm1_ewin333[119] = {
    2, 5, 1, 6, 12, 5, 1, 5, 0, 15, 3, 13, 2, 4, 2, 0, 4, 13, 11, 9, 4, 5,
    4, 13, 7, 15, 0, 11, 10, 7, 5, 4, 7, 0, 14, 11, 12, 10, 12, 4, 11, 2,
    5, 2, 10, 10, 7, 3, 14, 11, 8, 0, 15, 2, 2, 3, 10, 11, 6, 9, 2, 8, 15,
    4, 12, 13, 14, 13, 0, 7, 3, 3, 12, 3, 9, 8, 4, 6, 15, 0, 3, 9, 11, 14,
    5, 7, 4, 4, 11, 14, 8, 6, 11, 14, 0, 12, 13, 4, 12, 12, 11, 6, 13, 0,
    10, 8, 13, 6, 13, 15, 2, 5, 14, 13, 11, 11, 15, 0, 0 };

__constant uint tpm1_ewin500[23] = {
    0xd474c2d4, 0xb7330cfe, 0xb00f3a15, 0x74d81ab8, 0x23102ec9,
    0x1693c48d, 0x9845eb35, 0x9c0da860, 0x9477df49, 0x598d1d83,
    0x8c8bb315, 0xa5add55b, 0xb193f4f7, 0x90e6ec89, 0x2e477998,
    0x6fefd0b8, 0xce5273a3, 0x8952ee05, 0x0f1dc5dd, 0xa487cf80,
    0x484dd0af, 0x21644a26, 0xfd100000 };

__constant uint tpm1_map[60] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

static inline void
tpm1_expRL(uint *P, uint *in, uint *n, uint *one,
           uint m, uint rho)
{
    uint s[3];
    P[0]=one[0]; P[1]=one[1]; P[2]=one[2];
    s[0]=in[0];  s[1]=in[1];  s[2]=in[2];
    while (m > 0) {
        if (m & 1) montmul96(P, s, P, n, rho);
        montsqr96(s, s, n, rho);
        m >>= 1;
    }
}

static inline void
pm196_stage1(uint *P, uint *N, uint rho, uint *one, uint stg1)
{
    int i;
    uint g[16][3];

    g[0][0]=P[0]=one[0]; g[0][1]=P[1]=one[1]; g[0][2]=P[2]=one[2];
    for (i = 1; i < 16; i++)
        modadd96(g[i-1], g[i-1], g[i], N);

    switch (stg1) {
    case 100:
        for (i = 0; i < 34; i++) {
            montsqr96(P, P, N, rho);
            montsqr96(P, P, N, rho);
            montsqr96(P, P, N, rho);
            montsqr96(P, P, N, rho);
            if (tpm1_ewin100[i] > 0) montmul96(P, g[tpm1_ewin100[i]], P, N, rho);
        }
        break;
    case 333:
        for (i = 0; i < 119; i++) {
            montsqr96(P, P, N, rho);
            montsqr96(P, P, N, rho);
            montsqr96(P, P, N, rho);
            montsqr96(P, P, N, rho);
            if (tpm1_ewin333[i] > 0) montmul96(P, g[tpm1_ewin333[i]], P, N, rho);
        }
        break;
    case 500:
        for (i = 0; i < 23; i++) {
            int j;
            for (j = 7; j >= 0; j--) {
                montsqr96(P, P, N, rho);
                montsqr96(P, P, N, rho);
                montsqr96(P, P, N, rho);
                montsqr96(P, P, N, rho);
                montmul96(P, g[(tpm1_ewin500[i] >> (j*4)) & 0xf], P, N, rho);
            }
        }
        break;
    }
}

static inline void
pm196_stage2_pair(uint *P, uint *acc, uint *n, uint rho, uint b1)
{
    int w = 60;
    uint d[18][3], six[3], x12[3], xmid[3], x24[3], x25[3];
    uint x36[3], x60[3], x72[3], five[3], pw[3], pgiant[3];
    int i, j;
    uint b2 = 50 * b1;

    d[1][0]=P[0]; d[1][1]=P[1]; d[1][2]=P[2];
    montsqr96(P,    d[2],  n, rho);
    montmul96(P,    d[2],  six, n, rho);
    montsqr96(six,  six,   n, rho);
    montsqr96(six,  x12,   n, rho);
    montsqr96(x12,  x24,   n, rho);
    montmul96(x24,  x12,   x36, n, rho);
    montmul96(x24,  x36,   x60, n, rho);
    montsqr96(x36,  x72,   n, rho);
    montsqr96(d[2], five,  n, rho);
    montmul96(five, d[1],  five, n, rho);
    montmul96(x24,  d[1],  x25, n, rho);

    j = 1;
    xmid[0]=x12[0]; xmid[1]=x12[1]; xmid[2]=x12[2];
    while ((j + 6) < 50) {
        montmul96(d[tpm1_map[j]],   x36,  d[tpm1_map[j+6]], n, rho);
        montmul96(d[tpm1_map[j+6]], xmid, d[tpm1_map[j+6]], n, rho);
        montmul96(xmid, x72, xmid, n, rho);
        j += 6;
    }

    montmul96(x25,              x36,  d[tpm1_map[11]], n, rho);
    montmul96(d[tpm1_map[11]], x60,  d[tpm1_map[11]], n, rho);
    montmul96(x60, x72, xmid, n, rho);
    j = 11;
    while ((j + 6) < 60) {
        montmul96(d[tpm1_map[j]],   x36,  d[tpm1_map[j+6]], n, rho);
        montmul96(d[tpm1_map[j+6]], xmid, d[tpm1_map[j+6]], n, rho);
        montmul96(xmid, x72, xmid, n, rho);
        j += 6;
    }

    tpm1_expRL(pw, P, n, acc, 120*120, rho);

    uint x14400[3];
    x14400[0]=pw[0]; x14400[1]=pw[1]; x14400[2]=pw[2];

    uint x240x120[3];
    montsqr96(pw, x240x120, n, rho);
    xmid[0]=x240x120[0]; xmid[1]=x240x120[1]; xmid[2]=x240x120[2];
    pgiant[0]=pw[0]; pgiant[1]=pw[1]; pgiant[2]=pw[2];

    i = 2*w;
    while (i < (int)b1) {
        montmul96(pgiant, x14400,   pgiant, n, rho);
        montmul96(pgiant, xmid,     pgiant, n, rho);
        montmul96(xmid,   x240x120, xmid,   n, rho);
        i += 2*w;
    }

    while (i < (int)b2) {
        uint tmp[3], sub[3];

        modsub96(pgiant, d[1], sub, n);
        montmul96(acc, sub, acc, n, rho);
        for (j = 3; j < 18; j++) {
            modsub96(pgiant, d[j], sub, n);
            montmul96(acc, sub, acc, n, rho);
        }
        montmul96(pgiant, x14400,   pgiant, n, rho);
        montmul96(pgiant, xmid,     pgiant, n, rho);
        montmul96(xmid,   x240x120, xmid,   n, rho);
        i += 120;
    }
    /* suppress unused warning for tmp */
    // (void)tmp;
}

/* =========================================================================
 * P-1 kernel
 *
 * CUDA original: gbl_pm196(num, n_in, rho_in, one_in, f_out, B1, B2)
 * ========================================================================= */
__kernel void gbl_pm196(
    int num,
    __global uint *n_in,
    __global uint *rho_in,
    __global uint *one_in,
    __global uint *f_out,
    uint B1,
    uint B2)
{
    int idx = (int)get_global_id(0);
    if (idx >= num) return;

    uint rho = rho_in[idx];
    uint n[3], unityval[3], P[3], acc[3], gcd_val[3];

    n[0]        = n_in[idx*3+0];   n[1]        = n_in[idx*3+1];   n[2]        = n_in[idx*3+2];
    unityval[0] = one_in[idx*3+0]; unityval[1] = one_in[idx*3+1]; unityval[2] = one_in[idx*3+2];

    f_out[3*idx+0] = 1;
    f_out[3*idx+1] = 0;
    f_out[3*idx+2] = 0;

    pm196_stage1(P, n, rho, unityval, B1);
    modsub96(P, unityval, acc, n);
    gcd96(n, acc, gcd_val);

    if (((f_out[3*idx+0]==1) && (f_out[3*idx+1]==0) && (f_out[3*idx+2]==0)) ||
        ((f_out[3*idx+0]==n[0]) && (f_out[3*idx+1]==n[1]) && (f_out[3*idx+2]==n[2])))
    {
        f_out[3*idx+0] = gcd_val[0];
        f_out[3*idx+1] = gcd_val[1];
        f_out[3*idx+2] = gcd_val[2];
    }

    acc[0]=unityval[0]; acc[1]=unityval[1]; acc[2]=unityval[2];
    pm196_stage2_pair(P, acc, n, rho, B1);
    gcd96(n, acc, gcd_val);

    if (((f_out[3*idx+0]==1) && (f_out[3*idx+1]==0) && (f_out[3*idx+2]==0)) ||
        ((f_out[3*idx+0]==n[0]) && (f_out[3*idx+1]==n[1]) && (f_out[3*idx+2]==n[2])))
    {
        f_out[3*idx+0] = gcd_val[0];
        f_out[3*idx+1] = gcd_val[1];
        f_out[3*idx+2] = gcd_val[2];
    }
}
