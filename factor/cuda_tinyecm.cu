// ============================================================================
// cuda_ecm64.cu - CUDA kernel implementation (compile with nvcc)
// ============================================================================
#include <cuda_runtime.h>
#include <stdint.h>

#ifdef __CUDACC__
#include "cuda_intrinsics.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    uint64_t X;
    uint64_t Z;
} uecm_pt;

#define INV_2_POW_32 2.3283064365386962890625e-10

__device__ uint32_t lcg_rand_32b(uint32_t lower, uint32_t upper, uint64_t* ploc_lcg)
{
    *ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}

__device__ void uaddxz(uint32_t rho, uint64_t n, uint64_t p1x, uint64_t p1z,
    uint64_t p2x, uint64_t p2z, uint64_t pix, uint64_t piz, 
    uint64_t* pox, uint64_t* poz) {

    uint64_t diff1 = modsub64(p1x, p1z, n);
    uint64_t sum1 = modadd64(p1x, p1z, n);
    uint64_t diff2 = modsub64(p2x, p2z, n);
    uint64_t sum2 = modadd64(p2x, p2z, n);

    uint64_t tt1 = montmul64(diff1, sum2, n, rho); //U
    uint64_t tt2 = montmul64(sum1, diff2, n, rho); //V

    uint64_t tt3 = modadd64(tt1, tt2, n);
    uint64_t tt4 = modsub64(tt1, tt2, n);
    tt1 = montmul64(tt3, tt3, n, rho);   //(U + V)^2
    tt2 = montmul64(tt4, tt4, n, rho);   //(U - V)^2

    uint64_t tmpx = montmul64(tt1, piz, n, rho);     //Z * (U + V)^2
    uint64_t tmpz = montmul64(tt2, pix, n, rho);     //x * (U - V)^2

    *pox = tmpx;
    *poz = tmpz;
}

__device__ void uadd(uint32_t rho, uint64_t n, uecm_pt p1,
    uecm_pt p2, uecm_pt pi, uecm_pt *po) {

    uint64_t diff1 = modsub64(p1.X, p1.Z, n);
    uint64_t sum1 = modadd64(p1.X, p1.Z, n);
    uint64_t diff2 = modsub64(p2.X, p2.Z, n);
    uint64_t sum2 = modadd64(p2.X, p2.Z, n);

    uint64_t tt1 = montmul64(diff1, sum2, n, rho); //U
    uint64_t tt2 = montmul64(sum1, diff2, n, rho); //V

    uint64_t tt3 = modadd64(tt1, tt2, n);
    uint64_t tt4 = modsub64(tt1, tt2, n);
    tt1 = montmul64(tt3, tt3, n, rho);   //(U + V)^2
    tt2 = montmul64(tt4, tt4, n, rho);   //(U - V)^2

    uint64_t tmpx = montmul64(tt1, pi.Z, n, rho);     //Z * (U + V)^2
    uint64_t tmpz = montmul64(tt2, pi.X, n, rho);     //x * (U - V)^2

    po->X = tmpx;
    po->Z = tmpz;

}

__device__ void udup(uint64_t s, uint32_t rho, uint64_t n,
    uint64_t insum, uint64_t indiff, uecm_pt *P) {

    uint64_t tt1 = montmul64(indiff, indiff, n, rho);          // U=(x1 - z1)^2
    uint64_t tt2 = montmul64(insum, insum, n, rho);           // V=(x1 + z1)^2
    P->X = montmul64(tt1, tt2, n, rho);         // x=U*V

    uint64_t tt3 = modsub64(tt2, tt1, n);          // w = V-U
    tt2 = montmul64(tt3, s, n, rho);      // w = (A+2)/4 * w
    tt2 = modadd64(tt2, tt1, n);          // w = w + U

    P->Z = montmul64(tt2, tt3, n, rho);         // Z = w*(V-U)
}

__device__ void uprac(uint32_t rho, uint64_t n, uecm_pt *P,
    uint64_t c, double v, uint64_t s)
{
    uint64_t d, e, r;

    // we require c != 0
    int shift = __ffsll(c) - 1;

    d = c;
    r = (uint64_t)((double)d * v + 0.5);

    d = c - r;
    e = 2 * r - c;

    uint64_t s1, s2, d1, d2;
    uint64_t swp;
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
        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        {
            d = (2 * d - e) / 3;
            e = (e - d) / 2;

            uecm_pt pt4;
            uadd(rho, n, pt1, pt2, pt3, &pt4); // T = A + B (C)
            uecm_pt pt5;
            uadd(rho, n, pt4, pt1, pt2, &pt5); // T2 = T + A (B)
            uadd(rho, n, pt2, pt4, pt1, &pt2); // B = B + T (A)

            swp = pt1.X;
            pt1.X = pt5.X;
            pt5.X = swp;
            swp = pt1.Z;
            pt1.Z = pt5.Z;
            pt5.Z = swp;
        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        {
            d = (d - e) / 2;

            d1 = modsub64(pt1.X, pt1.Z, n);
            s1 = modadd64(pt1.X, pt1.Z, n);

            uadd(rho, n, pt1, pt2, pt3, &pt2);        // B = A + B (C)
            udup(s, rho, n, s1, d1, &pt1);        // A = 2A
        }
        else if ((d + 3) / 4 <= e)
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
        else if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;

            d2 = modsub64(pt1.X, pt1.Z, n);
            s2 = modadd64(pt1.X, pt1.Z, n);

            uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 2 == 0)
        {
            d /= 2;

            d2 = modsub64(pt1.X, pt1.Z, n);
            s2 = modadd64(pt1.X, pt1.Z, n);

            uadd(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)
            udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 3 == 0)
        {
            d = d / 3 - e;

            d1 = modsub64(pt1.X, pt1.Z, n);
            s1 = modadd64(pt1.X, pt1.Z, n);

            uecm_pt pt4;
            udup(s, rho, n, s1, d1, &pt4);        // T = 2A
            uecm_pt pt5;
            uadd(rho, n, pt1, pt2, pt3, &pt5);        // T2 = A + B (C)
            uadd(rho, n, pt4, pt1, pt1, &pt1);        // A = T + A (A)
            uadd(rho, n, pt4, pt5, pt3, &pt4);        // T = T + T2 (C)

            swp = pt3.X;
            pt3.X = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

        }
        else if ((d + e) % 3 == 0)
        {
            d = (d - 2 * e) / 3;

            uecm_pt pt4;
            uadd(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)


            d2 = modsub64(pt1.X, pt1.Z, n);
            s2 = modadd64(pt1.X, pt1.Z, n);
            uadd(rho, n, pt4, pt1, pt2, &pt2);        // B = T + A (B)
            udup(s, rho, n, s2, d2, &pt4);        // T = 2A
            uadd(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else if ((d - e) % 3 == 0)
        {
            d = (d - e) / 3;

            uecm_pt pt4;
            uadd(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)

            d2 = modsub64(pt1.X, pt1.Z, n);
            s2 = modadd64(pt1.X, pt1.Z, n);
            uadd(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)

            swp = pt2.X;
            pt2.X = pt4.X;
            pt4.X = swp;
            swp = pt2.Z;
            pt2.Z = pt4.Z;
            pt4.Z = swp;

            udup(s, rho, n, s2, d2, &pt4);        // T = 2A
            uadd(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else
        {
            e /= 2;

            d2 = modsub64(pt2.X, pt2.Z, n);
            s2 = modadd64(pt2.X, pt2.Z, n);

            uadd(rho, n, pt3, pt2, pt1, &pt3);        // C = C + B (A)
            udup(s, rho, n, s2, d2, &pt2);        // B = 2B
        }
    }
        
    uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)

    int i;
    for (i = 0; i < shift; i++)
    {
        d1 = modsub64(P->X, P->Z, n);
        s1 = modadd64(P->X, P->Z, n);
        udup(s, rho, n, s1, d1, P);     // P = 2P
    }


    return;
}

__device__ void uprac85(uint32_t rho, uint64_t n, uecm_pt* P, uint64_t s)
{
    uint64_t s1, s2, d1, d2;
    uint64_t swp;
    int i;
    uint8_t steps[146] = {
        0,6,0,6,0,6,0,6,0,4,
        6,0,4,6,0,4,4,6,0,4,
        4,6,0,5,4,6,0,3,3,4,
        6,0,3,5,4,6,0,3,4,3,
        4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,
        0,4,3,4,3,5,3,3,3,3,
        3,3,3,3,4,6,0,3,3,3,
        3,3,3,3,3,3,4,3,4,3,
        4,6,0,3,4,3,5,4,6,0,
        5,3,3,3,4,6,0,5,4,3,
        5,4,6,0,4,3,3,3,5,4,
        6,0,4,3,5,3,3,4,6,0,
        3,3,3,3,5,4,6,0,3,3,
        3,4,3,3,4,6 };

    uecm_pt pt1, pt2, pt3;
    for (i = 0; i < 146; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X = pt2.X = pt3.X = P->X;
            pt1.Z = pt2.Z = pt3.Z = P->Z;

            d1 = modsub64(pt1.X, pt1.Z, n);
            s1 = modadd64(pt1.X, pt1.Z, n);
            udup(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uecm_pt pt4;
            uadd(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            swp = pt1.X;
            pt1.X = pt4.X;
            pt4.X = pt3.X;
            pt3.X = pt2.X;
            pt2.X = swp;
            swp = pt1.Z;
            pt1.Z = pt4.Z;
            pt4.Z = pt3.Z;
            pt3.Z = pt2.Z;
            pt2.Z = swp;
        }
        else if (steps[i] == 4)
        {
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
        else if (steps[i] == 5)
        {
            d2 = modsub64(pt1.X, pt1.Z, n);
            s2 = modadd64(pt1.X, pt1.Z, n);

            uadd(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            udup(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uadd(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
}

__device__ uint64_t uecm_build(uecm_pt* P, uint32_t rho, uint64_t n,
    uint32_t sigma, uint64_t* ps, uint64_t five, uint64_t Rsqr)
{
    uint64_t t1, t2, t3, t4;
    uint64_t u, v;

    u = montmul64((uint64_t)sigma, Rsqr, n, rho);  // to_monty(sigma)

    //printf("sigma = %" PRIu64 ", u = %" PRIu64 ", n = %" PRIu64 "\n", sigma, u, n);

    v = modadd64(u, u, n);
    v = modadd64(v, v, n);            // 4*sigma

    //printf("v = %" PRIu64 "\n", v);

    u = montmul64(u, u, n, rho);
    t1 = five;

    //printf("monty(5) = %" PRIu64 "\n", t1);

    u = modsub64(u, t1, n);           // sigma^2 - 5

    //printf("u = %" PRIu64 "\n", u);

    t1 = montmul64(u, u, n, rho);
    uint64_t tmpx = montmul64(t1, u, n, rho);  // u^3

    uint64_t v2 = modadd64(v, v, n);             // 2*v
    uint64_t v4 = modadd64(v2, v2, n);           // 4*v
    uint64_t v8 = modadd64(v4, v4, n);           // 8*v
    uint64_t v16 = modadd64(v8, v8, n);          // 16*v
    uint64_t t5 = montmul64(v16, tmpx, n, rho);    // 16*u^3*v

    t1 = montmul64(v, v, n, rho);
    uint64_t tmpz = montmul64(t1, v, n, rho);  // v^3

    //compute parameter A
    t1 = modsub64(v, u, n);           // (v - u)
    t2 = montmul64(t1, t1, n, rho);
    t4 = montmul64(t2, t1, n, rho);   // (v - u)^3

    t1 = modadd64(u, u, n);           // 2u
    t2 = modadd64(u, v, n);           // u + v
    t3 = modadd64(t1, t2, n);         // 3u + v

    t1 = montmul64(t3, t4, n, rho);   // a = (v-u)^3 * (3u + v)

    // u holds the denom (jeff note: isn't it t5 that has the denom?)
    // t1 holds the numer
    // accomplish the division by multiplying by the modular inverse
    t2 = 1;
    t5 = montmul64(t5, t2, n, rho);   // take t5 out of monty rep

    uint64_t gcd;
    t3 = modinv64(t5, n, (unsigned long long *)&gcd);

    t3 = montmul64(t3, Rsqr, n, rho); // to_monty(t3)
    *ps = montmul64(t3, t1, n, rho);

    P->X = tmpx;
    P->Z = tmpz;

    return gcd;
}

__device__ void uecm_stage1(uint32_t rho, uint64_t n, uecm_pt *P,
    uint32_t stg1, uint64_t s)
{
    uint32_t q;

    // handle the only even case
    q = 2;
    while (q < (stg1 * 4))  // jeff: multiplying by 4 improves perf ~1%
    {
        uint64_t diff1 = modsub64(P->X, P->Z, n);
        uint64_t sum1 = modadd64(P->X, P->Z, n);
        udup(s, rho, n, sum1, diff1, P);
        q *= 2;
    }

    if (stg1 == 27)
    {
        uprac(rho, n, P, 3, 0.61803398874989485, s);
        uprac(rho, n, P, 3, 0.61803398874989485, s);
        uprac(rho, n, P, 3, 0.61803398874989485, s);
        uprac(rho, n, P, 5, 0.618033988749894903, s);
        uprac(rho, n, P, 5, 0.618033988749894903, s);
        uprac(rho, n, P, 7, 0.618033988749894903, s);
        uprac(rho, n, P, 11, 0.580178728295464130, s);
        uprac(rho, n, P, 13, 0.618033988749894903, s);
        uprac(rho, n, P, 17, 0.618033988749894903, s);
        uprac(rho, n, P, 19, 0.618033988749894903, s);
        uprac(rho, n, P, 23, 0.522786351415446049, s);
    }
    else 
    {
        uprac85(rho, n, P, s);

        if (stg1 == 85)
        {
            uprac(rho, n, P, 61, 0.522786351415446049, s);
        }
        else
        {
            uprac(rho, n, P, 5, 0.618033988749894903, s);
            uprac(rho, n, P, 11, 0.580178728295464130, s);
            uprac(rho, n, P, 89, 0.618033988749894903, s);
            uprac(rho, n, P, 97, 0.723606797749978936, s);
            uprac(rho, n, P, 101, 0.556250337855490828, s);
            uprac(rho, n, P, 107, 0.580178728295464130, s);
            uprac(rho, n, P, 109, 0.548409048446403258, s);
            uprac(rho, n, P, 113, 0.618033988749894903, s);

            if (stg1 == 125)
            {
                // jeff: moved 61 to here
                uprac(rho, n, P, 61, 0.522786351415446049, s);
                uprac(rho, n, P, 103, 0.632839806088706269, s);
            }
            else
            {
                uprac(rho, n, P, 7747, 0.552188778811121, s); // 61 x 127
                uprac(rho, n, P, 131, 0.618033988749894903, s);
                uprac(rho, n, P, 14111, 0.632839806088706, s);  // 103 x 137
                uprac(rho, n, P, 20989, 0.620181980807415, s);  // 139 x 151
                uprac(rho, n, P, 157, 0.640157392785047019, s);
                uprac(rho, n, P, 163, 0.551390822543526449, s);

                if (stg1 == 165)
                {
                    uprac(rho, n, P, 149, 0.580178728295464130, s);
                }
                else
                {
                    uprac(rho, n, P, 13, 0.618033988749894903, s);
                    uprac(rho, n, P, 167, 0.580178728295464130, s);
                    uprac(rho, n, P, 173, 0.612429949509495031, s);
                    uprac(rho, n, P, 179, 0.618033988749894903, s);
                    uprac(rho, n, P, 181, 0.551390822543526449, s);
                    uprac(rho, n, P, 191, 0.618033988749894903, s);
                    uprac(rho, n, P, 193, 0.618033988749894903, s);
                    uprac(rho, n, P, 29353, 0.580178728295464, s);  // 149 x 197
                    uprac(rho, n, P, 199, 0.551390822543526449, s);
                }
            }
        }
    }

#if 0
        else if (stg1 == 47)
        {
            // jeff: improved perf slightly by using one more uprac for 3,
            // and removing uprac for 47.
            uprac(rho, n, P, 3, 0.618033988749894903, s, size);
            uprac(rho, n, P, 3, 0.618033988749894903, s, size);
            uprac(rho, n, P, 3, 0.618033988749894903, s, size);
            uprac(rho, n, P, 3, 0.618033988749894903, s, size);
            uprac(rho, n, P, 5, 0.618033988749894903, s, size);
            uprac(rho, n, P, 5, 0.618033988749894903, s, size);
            uprac(rho, n, P, 7, 0.618033988749894903, s, size);
            uprac(rho, n, P, 11, 0.580178728295464130, s, size);
            uprac(rho, n, P, 13, 0.618033988749894903, s, size);
            uprac(rho, n, P, 17, 0.618033988749894903, s, size);
            uprac(rho, n, P, 19, 0.618033988749894903, s, size);
            uprac(rho, n, P, 23, 0.522786351415446049, s, size);
            uprac(rho, n, P, 29, 0.548409048446403258, s, size);
            uprac(rho, n, P, 31, 0.618033988749894903, s, size);
            uprac(rho, n, P, 37, 0.580178728295464130, s, size);
            uprac(rho, n, P, 41, 0.548409048446403258, s, size);
            uprac(rho, n, P, 43, 0.618033988749894903, s, size);
            //        uprac(rho, n, P, 47, 0.548409048446403258, s, size);
        }
        else if (stg1 == 59)
        {   // jeff: probably stg1 of 59 would benefit from similar changes
            // as stg1 of 47 above, but I didn't bother. Stg1 of 59 seems to
            // always perform worse than stg1 of 47, so there doesn't seem
            // to be any reason to ever use stg1 of 59.
            uprac(rho, n, P, 3, 0.61803398874989485, s, size);
            uprac(rho, n, P, 3, 0.61803398874989485, s, size);
            uprac(rho, n, P, 3, 0.61803398874989485, s, size);
            uprac(rho, n, P, 5, 0.618033988749894903, s, size);
            uprac(rho, n, P, 5, 0.618033988749894903, s, size);
            uprac(rho, n, P, 7, 0.618033988749894903, s, size);
            uprac(rho, n, P, 7, 0.618033988749894903, s, size);
            uprac(rho, n, P, 11, 0.580178728295464130, s, size);
            uprac(rho, n, P, 13, 0.618033988749894903, s, size);
            uprac(rho, n, P, 17, 0.618033988749894903, s, size);
            uprac(rho, n, P, 19, 0.618033988749894903, s, size);
            uprac(rho, n, P, 23, 0.522786351415446049, s, size);
            uprac(rho, n, P, 29, 0.548409048446403258, s, size);
            uprac(rho, n, P, 31, 0.618033988749894903, s, size);
            uprac(rho, n, P, 1961, 0.552936068843375, s, size);   // 37 * 53
            uprac(rho, n, P, 41, 0.548409048446403258, s, size);
            uprac(rho, n, P, 43, 0.618033988749894903, s, size);
            uprac(rho, n, P, 47, 0.548409048446403258, s, size);
            uprac(rho, n, P, 59, 0.548409048446403258, s, size);
        }
        else if (stg1 == 70)
        {
            // call prac with best ratio found in deep search.
            // some composites are cheaper than their
            // constituent primes.
            //uprac70(rho, n, P, s, size);
        }
        else // if (stg1 >= 85)
        {
            //uprac85(rho, n, P, s, size);

            if (stg1 == 85)
            {
                uprac(rho, n, P, 61, 0.522786351415446049, s, size);
            }
            else
            {
                uprac(rho, n, P, 5, 0.618033988749894903, s, size);
                uprac(rho, n, P, 11, 0.580178728295464130, s, size);
                //            uprac(rho, n, P, 61, 0.522786351415446049, s, size);
                uprac(rho, n, P, 89, 0.618033988749894903, s, size);
                uprac(rho, n, P, 97, 0.723606797749978936, s, size);
                uprac(rho, n, P, 101, 0.556250337855490828, s, size);
                uprac(rho, n, P, 107, 0.580178728295464130, s, size);
                uprac(rho, n, P, 109, 0.548409048446403258, s, size);
                uprac(rho, n, P, 113, 0.618033988749894903, s, size);

                if (stg1 == 125)
                {
                    // jeff: moved 61 to here
                    uprac(rho, n, P, 61, 0.522786351415446049, s, size);
                    uprac(rho, n, P, 103, 0.632839806088706269, s, size);
                }
                else
                {
                    uprac(rho, n, P, 7747, 0.552188778811121, s, size); // 61 x 127
                    uprac(rho, n, P, 131, 0.618033988749894903, s, size);
                    uprac(rho, n, P, 14111, 0.632839806088706, s, size);  // 103 x 137
                    uprac(rho, n, P, 20989, 0.620181980807415, s, size);  // 139 x 151
                    uprac(rho, n, P, 157, 0.640157392785047019, s, size);
                    uprac(rho, n, P, 163, 0.551390822543526449, s, size);

                    if (stg1 == 165)
                    {
                        uprac(rho, n, P, 149, 0.580178728295464130, s, size);
                    }
                    else
                    {
                        uprac(rho, n, P, 13, 0.618033988749894903, s, size);
                        uprac(rho, n, P, 167, 0.580178728295464130, s, size);
                        uprac(rho, n, P, 173, 0.612429949509495031, s, size);
                        uprac(rho, n, P, 179, 0.618033988749894903, s, size);
                        uprac(rho, n, P, 181, 0.551390822543526449, s, size);
                        uprac(rho, n, P, 191, 0.618033988749894903, s, size);
                        uprac(rho, n, P, 193, 0.618033988749894903, s, size);
                        uprac(rho, n, P, 29353, 0.580178728295464, s, size);  // 149 x 197
                        uprac(rho, n, P, 199, 0.551390822543526449, s, size);
                    }
                }
            }
        }
#endif

    return;
}

__device__ uint64_t crossprodxz(uint64_t ptx, uint64_t ptz, 
    uint64_t Pbx, uint64_t Pbz, uint64_t Pbprod,
    uint64_t ptprod, uint64_t acc, uint32_t rho, uint64_t n)
{
    // accumulate the cross product  (zimmerman syntax).
    // page 342 in C&P
    uint64_t tt1 = modsub64(ptx, Pbx, n);
    uint64_t tt2 = modadd64(ptz, Pbz, n);
    uint64_t tt3 = montmul64(tt1, tt2, n, rho);
    tt1 = modadd64(tt3, Pbprod, n);
    tt2 = modsub64(tt1, ptprod, n);

    acc = montmul64(acc, tt2, n, rho);

    return acc;
}

__device__ uint64_t uecm_stage2_D30(uecm_pt* P, uint32_t rho, uint64_t n,
    uint32_t B1, uint64_t s, uint64_t unityval)
{
    int b;
    int i, j, k;
    uecm_pt Pa;
    uecm_pt Pstep;
    uecm_pt Pdiff;
    uint64_t result;
    uecm_pt pt5, pt6;

    int idx = 1; // threadIdx.x;

    // Q = P = result of stage 1
    // Compute small differences: 1, 7, 11, 13, 17, 19, 23, 29
    //__shared__ uint64_t Pbprod[8 * 128];
    //__shared__ uint64_t PbX[8 * 128];   
    //__shared__ uint64_t PbZ[8 * 128];
    uint64_t Pbprod[8];
    uint64_t PbX[8];
    uint64_t PbZ[8];
    uint64_t diff1;
    uint64_t sum1;

    // common init
    {
        // [1]Q
        PbX[0 * idx] = P->X;
        PbZ[0 * idx] = P->Z;

        // [2]Q
        diff1 = modsub64(P->X, P->Z, n);
        sum1 = modadd64(P->X, P->Z, n);
        udup(s, rho, n, sum1, diff1, &Pa);
        PbX[1 * idx] = Pa.X;
        PbZ[1 * idx] = Pa.Z;

        // [2]Q + [1]Q([1]Q) = [3]Q
        uaddxz(rho, n,
            PbX[0 * idx], PbZ[0 * idx],
            PbX[1 * idx], PbZ[1 * idx],
            PbX[0 * idx], PbZ[0 * idx],
            &PbX[3 * idx], &PbZ[3 * idx]);  // <-- temporary

        // 2*[3]Q = [6]Q
        diff1 = modsub64(PbX[3 * idx], PbZ[3 * idx], n);
        sum1 = modadd64(PbX[3 * idx], PbZ[3 * idx], n);
        udup(s, rho, n, sum1, diff1, &pt6);   // pt6 = [6]Q

        // [3]Q + [2]Q([1]Q) = [5]Q
        uaddxz(rho, n, PbX[3 * idx], PbZ[3 * idx],
            PbX[1 * idx], PbZ[1 * idx],
            PbX[0 * idx], PbZ[0 * idx],
            &pt5.X, &pt5.Z);    // <-- pt5 = [5]Q

        // [6]Q + [5]Q([1]Q) = [11]Q
        uaddxz(rho, n, pt6.X, pt6.Z,
            pt5.X, pt5.Z,
            PbX[0 * idx], PbZ[0 * idx],
            &PbX[2 * idx], &PbZ[2 * idx]);    // <-- [11]Q

        // [11]Q + [6]Q([5]Q) = [17]Q
        uaddxz(rho, n,
            PbX[2 * idx], PbZ[2 * idx],
            pt6.X, pt6.Z, pt5.X, pt5.Z,
            &PbX[4 * idx], &PbZ[4 * idx]);    // <-- [17]Q

        // [17]Q + [6]Q([11]Q) = [23]Q
        uaddxz(rho, n, PbX[4 * idx], PbZ[4 * idx],
            pt6.X, pt6.Z, PbX[2 * idx], PbZ[2 * idx],
            &PbX[6 * idx], &PbZ[6 * idx]);    // <-- [23]Q

        // [23]Q + [6]Q([17]Q) = [29]Q
        uaddxz(rho, n,
            PbX[6 * idx], PbZ[6 * idx],
            pt6.X, pt6.Z, PbX[4 * idx], PbZ[4 * idx],
            &PbX[7 * idx], &PbZ[7 * idx]);    // <-- [29]Q

        // [6]Q + [1]Q([5]Q) = [7]Q
        uaddxz(rho, n, pt6.X, pt6.Z,
            PbX[0 * idx], PbZ[0 * idx], pt5.X, pt5.Z,
            &PbX[1 * idx], &PbZ[1 * idx]);    // <-- [7]Q

        // [7]Q + [6]Q([1]Q) = [13]Q
        uaddxz(rho, n,
            PbX[1 * idx], PbZ[1 * idx], pt6.X, pt6.Z,
            PbX[0 * idx], PbZ[0 * idx],
            &PbX[3 * idx], &PbZ[3 * idx]);    // <-- [13]Q

        // [13]Q + [6]Q([7]Q) = [19]Q
        uaddxz(rho, n,
            PbX[3 * idx], PbZ[3 * idx], pt6.X, pt6.Z,
            PbX[1 * idx], PbZ[1 * idx],
            &PbX[5 * idx], &PbZ[5 * idx]);    // <-- [19]Q

        // 4*[1]Q = [4]Q
        diff1 = modsub64(PbX[0 * idx], PbZ[0 * idx], n);
        sum1 = modadd64(PbX[0 * idx], PbZ[0 * idx], n);
        udup(s, rho, n, sum1, diff1, &pt5);   // pt5 = [2]Q

        diff1 = modsub64(pt5.X, pt5.Z, n);
        sum1 = modadd64(pt5.X, pt5.Z, n);
        udup(s, rho, n, sum1, diff1, &pt5);   // pt5 = [4]Q

        // Pd = [2w]Q
        // [17]Q + [13]Q([4]Q) = [30]Q
        uaddxz(rho, n, PbX[4 * idx], PbZ[4 * idx],
            PbX[3 * idx], PbZ[3 * idx], pt5.X, pt5.Z,
            &Pdiff.X, &Pdiff.Z);   // <-- [30]Q

        Pstep.X = Pdiff.X;
        Pstep.Z = Pdiff.Z;
        diff1 = modsub64(Pstep.X, Pstep.Z, n);
        sum1 = modadd64(Pstep.X, Pstep.Z, n);
        udup(s, rho, n, sum1, diff1, &Pstep);   // Pstep = [60]Q

        Pdiff.X = Pstep.X;
        Pdiff.Z = Pstep.Z;

        Pa.X = Pdiff.X;
        Pa.Z = Pdiff.Z;
        diff1 = modsub64(Pa.X, Pa.Z, n);
        sum1 = modadd64(Pa.X, Pa.Z, n);
        udup(s, rho, n, sum1, diff1, &Pa);   // Pa = [120]Q

        // Now we have Pa = 120, Pstep = 60 and Pdiff = 60.  ready for giant step.

        // make all of the Pbprod's
        Pbprod[0 * idx] = montmul64(PbX[0 * idx], PbZ[0 * idx], n, rho);
        Pbprod[1 * idx] = montmul64(PbX[1 * idx], PbZ[1 * idx], n, rho);
        Pbprod[2 * idx] = montmul64(PbX[2 * idx], PbZ[2 * idx], n, rho);
        Pbprod[3 * idx] = montmul64(PbX[3 * idx], PbZ[3 * idx], n, rho);
        Pbprod[4 * idx] = montmul64(PbX[4 * idx], PbZ[4 * idx], n, rho);
        Pbprod[5 * idx] = montmul64(PbX[5 * idx], PbZ[5 * idx], n, rho);
        Pbprod[6 * idx] = montmul64(PbX[6 * idx], PbZ[6 * idx], n, rho);
        Pbprod[7 * idx] = montmul64(PbX[7 * idx], PbZ[7 * idx], n, rho);

    }

    // __syncthreads();

    // ---------------------------------------------------------------------
    // to here is the same for any B1.
    // ---------------------------------------------------------------------

    // advance giant step to one step beyond B1.
    uint32_t aval = 120;
    while (aval < B1)
    {
        pt5.X = Pa.X;
        pt5.Z = Pa.Z;
        uaddxz(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, &Pa.X, &Pa.Z);
        Pdiff.X = pt5.X;
        Pdiff.Z = pt5.Z;
        aval += 60;
    }

    //initialize Paprod
    uint64_t Paprod = montmul64(Pa.X, Pa.Z, n, rho);

    //initialize accumulator
    uint64_t acc = unityval;
    result = 0;

    while (aval < (25 * B1))
    {
        int j;
        uint64_t tmp;

        for (j = 0; j < 8; j++)
        {
            tmp = crossprodxz(Pa.X, Pa.Z, PbX[j * idx], PbZ[j * idx],
                Pbprod[j * idx], Paprod, acc, rho, n);
            if ((tmp == 0) && (result == 0)) result = acc;
            acc = tmp;
        }

        // giant step - use the addition formula for ECM
        pt5.X = Pa.X;
        pt5.Z = Pa.Z;
        uaddxz(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, &Pa.X, &Pa.Z);
        Pdiff.X = pt5.X;
        Pdiff.Z = pt5.Z;
        aval += 60;
        Paprod = montmul64(Pa.X, Pa.Z, n, rho);

    }

    if (result == 0)
        result = acc;

    // __syncthreads();

    return result;
}

__device__ uint64_t gcd64(uint64_t u, uint64_t v)
{
    uint64_t a, b, c;
    a = u; b = v;
    while (b != 0)
    {
        c = a % b;
        a = b;
        b = c;
    }

    return a;

    // this has problems on the gpu as-is.  for now the
    // division based one is working just fine.
    // 
    //uint64_t retval = 1;
    //if (u == 0) {
    //    retval = v;
    //}
    //else if (v != 0) {
    //    int j = __ffsll(v) - 1;
    //    v = (uint64_t)(v >> j);
    //    while (1) {
    //        uint64_t tmp = u;
    //        uint64_t sub1 = (uint64_t)(v - tmp);
    //        uint64_t sub2 = (uint64_t)(tmp - v);
    //        if (tmp == v)
    //            break;
    //        u = (tmp >= v) ? v : tmp;
    //        v = (tmp >= v) ? sub2 : sub1;
    //        // For the line below, the standard way to write this algorithm
    //        // would have been to use _trail_zcnt64(v)  (instead of
    //        // _trail_zcnt64(sub1)).  However, as pointed out by
    //        // https://gmplib.org/manual/Binary-GCD, "in twos complement the
    //        // number of low zero bits on u-v is the same as v-u, so counting or
    //        // testing can begin on u-v without waiting for abs(u-v) to be
    //        // determined."  Hence we are able to use sub1 for the argument.
    //        // By removing the dependency on abs(u-v), the CPU can execute
    //        // _trail_zcnt64() at the same time as abs(u-v).
    //        j = __ffsll(sub1) - 1;
    //        v = (uint64_t)(v >> j);
    //    }
    //    retval = u;
    //}
    //
    //return retval;
}

__global__ void gbl_init(int num, uint64_t* n,
    uint32_t* rho, uint64_t* unity, uint64_t *rsq)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < num) {
        uint64_t tmp1;

        rho[idx] = montmul32_w(n[idx] & 0xffffffff);
        uint32_t locrho = rho[idx];

        // uint64_t unityval = u64div(1, n);
        // Let R = 2^64.  We can see R%n == (R-n)%n  (mod n)
        unity[idx] = ((uint64_t)0 - n[idx]) % n[idx];   // unityval == R  (mod n)
        uint64_t unityval = unity[idx];

        uint64_t two = modadd64(unityval, unityval, n[idx]);
        uint64_t four = modadd64(two, two, n[idx]);
        uint64_t five = modadd64(unityval, four, n[idx]);
        uint64_t eight = modadd64(four, four, n[idx]);
        uint64_t sixteen = modadd64(eight, eight, n[idx]);
        uint64_t two_8 = montmul64(sixteen, sixteen, n[idx], locrho);   // R*2^8          (mod n)
        uint64_t two_16 = montmul64(two_8, two_8, n[idx], locrho);      // R*2^16         (mod n)
        uint64_t two_32 = montmul64(two_16, two_16, n[idx], locrho);    // R*2^32         (mod n)
        rsq[idx] = montmul64(two_32, two_32, n[idx], locrho);      // R*2^64 == R*R  (mod n)
    }

    return;
}


__global__ void gbl_ecm(int num, uint64_t *n_in, uint32_t* rho_in, uint64_t* one_in, 
    uint64_t* rsq_in, uint32_t *sigma_in, uint64_t* f_out, uint32_t stg1, int curve)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
    int numf = 0;

    if (idx < num)
    {
        uint64_t rho = rho_in[idx];
        uint64_t n = n_in[idx];
        uint64_t unityval = one_in[idx];
        uint64_t Rsqr = rsq_in[idx];

        uecm_pt P;
        uint64_t s;
        uint64_t likely_gcd = 1;

        //uint32_t smallsigma[8] = { 11, 61, 56, 81, 83, 7, 30, 51 };
        uint64_t Z[8] = { 85184, 14526784, 11239424, 34012224, 36594368, 21952, 1728000, 8489664 };
        uint64_t X[8] = { 1560896ULL, 51312965696ULL, 30693697091ULL, 281784327616ULL,
            326229015104ULL, 85184ULL, 716917375ULL, 17495004736ULL };
        uint64_t negt1[8] = { 146313216ULL, 476803160866816ULL, 236251574395731ULL, 4838810084806656ULL,
                              5902145938870272ULL, 655360ULL, 1305683671875ULL, 109380272541696ULL };
        uint64_t d[8] = { 1098870784ULL, 200325818077184ULL, 110006210374144ULL, 1460769954361344ULL,
            1732928528232448ULL, 38162432ULL, 1376481360000ULL, 57103695458304ULL };

        f_out[idx] = 1;
        //int lastcurve = curve + 2;
        //for (; curve < lastcurve; curve++)
        {
            uint64_t gcd;

            if (curve < 8)
            {
                // lookup point
                P.X = X[curve];
                P.Z = Z[curve];

                // some computation left to do for S parameter for this 'n'
                uint64_t num;

                uint64_t dem = d[curve];
                num = montmul64(negt1[curve], Rsqr, n, rho);     // to Monty rep.

                // The mulredc postcondition guarantees  num < n.
                num = n - num;

                dem = modinv64(dem, n, (unsigned long long*) &gcd);
                dem = montmul64(dem, Rsqr, n, rho);              // to Monty rep.
                s = montmul64(num, dem, n, rho);
                P.X = montmul64(P.X, Rsqr, n, rho);              // to Monty rep.
                P.Z = montmul64(P.Z, Rsqr, n, rho);              // to Monty rep.
            }
            else
            {
                uint64_t two = modadd64(unityval, unityval, n);
                uint64_t four = modadd64(two, two, n);
                uint64_t five = modadd64(unityval, four, n);

                gcd = uecm_build(&P, rho, n, sigma_in[idx], &s, five, Rsqr);
            }

            if (likely_gcd > 1)
            {
                // If the gcd gave us a factor, we're done.  If not, since gcd != 1
                // the inverse calculated in uecm_build would have bogus, and so this
                // curve is probably set up for failure (hence we continue).
                if (likely_gcd < n) // || n % likely_gcd != 0)
                {
                    if ((f_out[idx] == 1) || (f_out[idx] == n))
                    {
                        // if we haven't already found a factor, assign this result.
                        f_out[idx] = likely_gcd;
                    }
                }
            }

            uecm_stage1(rho, n, &P, stg1, s);

            uint64_t result = gcd64(n, P.Z);

            if ((f_out[idx] == 1) || (f_out[idx] == n))
            {
                // if we haven't already found a factor, assign this result.
                f_out[idx] = result;
            }

            if (1)
            {
                uint64_t stg2acc = uecm_stage2_D30(&P, rho, n, stg1, s, unityval);

                result = gcd64(stg2acc, n);

                if ((f_out[idx] == 1) || (f_out[idx] == n))
                {
                    // if we haven't already found a factor, assign this result.
                    f_out[idx] = result;
                }
            }

            // new curve
            //sigma_in[idx]++;
            uint64_t sigma64 = sigma_in[idx];
            sigma_in[idx] = lcg_rand_32b(7, (uint32_t)-1, &sigma64);
        }
    }

    return;
}

__global__ void gbl_ecm_xy(int num, uint64_t* n_in, uint32_t* rho_in, uint64_t* one_in,
    uint64_t* rsq_in, uint32_t* sigma_in, uint64_t* f_out, uint32_t stg1, int curve)
{
    int gidx = blockIdx.x * blockDim.x * blockDim.y + threadIdx.x;
    int tid = threadIdx.x;
    int numf = 0;

    //if (gidx < num) {
    {
        // get a number to work on from global memory
        uint64_t tmp1;

        uint64_t n = n_in[gidx];
        uint32_t rho = montmul32_w(n & 0xffffffff);
        uint64_t unityval = ((uint64_t)0 - n) % n;   // unityval == R  (mod n)
        uint64_t two = modadd64(unityval, unityval, n);
        uint64_t four = modadd64(two, two, n);
        uint64_t five = modadd64(unityval, four, n);
        uint64_t eight = modadd64(four, four, n);
        uint64_t sixteen = modadd64(eight, eight, n);
        uint64_t two_8 = montmul64(sixteen, sixteen, n, rho);   // R*2^8          (mod n)
        uint64_t two_16 = montmul64(two_8, two_8, n, rho);      // R*2^16         (mod n)
        uint64_t two_32 = montmul64(two_16, two_16, n, rho);    // R*2^32         (mod n)
        uint64_t Rsqr = montmul64(two_32, two_32, n, rho);      // R*2^64 == R*R  (mod n)

        uecm_pt P;
        uint64_t s;
        uint64_t likely_gcd = 1;

        //uint32_t smallsigma[8] = { 11, 61, 56, 81, 83, 7, 30, 51 };
        uint64_t Z[8] = { 85184, 14526784, 11239424, 34012224, 36594368, 21952, 1728000, 8489664 };
        uint64_t X[8] = { 1560896ULL, 51312965696ULL, 30693697091ULL, 281784327616ULL,
            326229015104ULL, 85184ULL, 716917375ULL, 17495004736ULL };
        uint64_t negt1[8] = { 146313216ULL, 476803160866816ULL, 236251574395731ULL, 4838810084806656ULL,
                              5902145938870272ULL, 655360ULL, 1305683671875ULL, 109380272541696ULL };
        uint64_t d[8] = { 1098870784ULL, 200325818077184ULL, 110006210374144ULL, 1460769954361344ULL,
            1732928528232448ULL, 38162432ULL, 1376481360000ULL, 57103695458304ULL };

        f_out[gidx] = 1;
        //int lastcurve = curve + 2;
        for (curve = 0; curve < 40; curve++)
        {
            uint64_t gcd;

            if (curve < 8)
            {
                //ubuild(&P, rho, &work, goodsigma[curve]); // sigma);
                //sigma = smallsigma[curve];
                // lookup point
                P.X = X[curve];
                P.Z = Z[curve];
                // some computation left to do for S parameter for this 'n'

                uint64_t num, uvc, u3vc;

                uint64_t dem = d[curve];
                // jeff note:  uecm_modinv_64(dem, n) appears to require dem < n,
                // so I now take the remainder to achieve dem < n.

                // This is a faster way to compute dem = dem % n, even if the CPU
                // has extremely fast division (as present in many new CPUs).
                dem = montmul64(dem, unityval, n, rho);

                num = montmul64(negt1[curve], Rsqr, n, rho);     // to Monty rep.
                // The mulredc postcondition guarantees  num < n.
                num = n - num;

                dem = modinv64(dem, n, (unsigned long long*) &gcd);

                dem = montmul64(dem, Rsqr, n, rho);              // to Monty rep.
                s = montmul64(num, dem, n, rho);

                P.X = montmul64(P.X, Rsqr, n, rho);              // to Monty rep.
                P.Z = montmul64(P.Z, Rsqr, n, rho);              // to Monty rep.
            }
            else
            {
                uint64_t two = modadd64(unityval, unityval, n);
                uint64_t four = modadd64(two, two, n);
                uint64_t five = modadd64(unityval, four, n);

                gcd = uecm_build(&P, rho, n, sigma_in[gidx], &s, five, Rsqr);
            }

            if (likely_gcd > 1)
            {
                // If the gcd gave us a factor, we're done.  If not, since gcd != 1
                // the inverse calculated in uecm_build would have bogus, and so this
                // curve is probably set up for failure (hence we continue).
                if (likely_gcd < n) // || n % likely_gcd != 0)
                {
                    if ((f_out[gidx] == 1) || (f_out[gidx] == n))
                    {
                        // if we haven't already found a factor, assign this result.
                        f_out[gidx] = likely_gcd;
                    }
                }
            }

            uecm_stage1(rho, n, &P, stg1, s);

            uint64_t result = gcd64(n, P.Z);

            if ((f_out[gidx] == 1) || (f_out[gidx] == n))
            {
                // if we haven't already found a factor, assign this result.
                f_out[gidx] = result;
            }

            if (1)
            {
                uint64_t stg2acc = uecm_stage2_D30(&P, rho, n, stg1, s, unityval);

                result = gcd64(stg2acc, n);

                if ((f_out[gidx] == 1) || (f_out[gidx] == n))
                {
                    // if we haven't already found a factor, assign this result.
                    f_out[gidx] = result;
                }
            }

            // new curve
            sigma_in[gidx]++;

            __syncthreads();

            if (tid == 0)
            {
                // nominate the first thread in the block to 
                // swap all factored results
                // to the end of the block's list.  if we have
                // fewer than threads/2 left, exit so more 
                // efficient kernels can finish the job.
                int i;
                int gbl_idx_start = blockIdx.x * blockDim.x * blockDim.y;
                int gbl_idx_end = gbl_idx_start + blockDim.x * blockDim.y - numf;
                for (i = 0; i < blockDim.x; i++)
                {
                    if ((f_out[gbl_idx_start + i] > 1) &&
                        (f_out[gbl_idx_start + i] < n_in[gbl_idx_start + i]))
                    {
                        // swap this factored modulus to the end of the list
                        uint64_t m = n_in[gbl_idx_start + i];
                        n_in[gbl_idx_start + i] = n_in[gbl_idx_end - 1];
                        n_in[gbl_idx_end - 1] = m;

                        m = rsq_in[gbl_idx_start + i];
                        rsq_in[gbl_idx_start + i] = rsq_in[gbl_idx_end - 1];
                        rsq_in[gbl_idx_end - 1] = m;

                        m = one_in[gbl_idx_start + i];
                        one_in[gbl_idx_start + i] = one_in[gbl_idx_end - 1];
                        one_in[gbl_idx_end - 1] = m;

                        uint32_t r = rho_in[gbl_idx_start + i];
                        rho_in[gbl_idx_start + i] = rho_in[gbl_idx_end - 1];
                        rho_in[gbl_idx_end - 1] = r;

                        gbl_idx_end--;
                        numf++;
                    }
                }
            }

            __syncthreads();

            n = n_in[gidx];
            rho = montmul32_w(n & 0xffffffff);
            unityval = ((uint64_t)0 - n) % n;   // unityval == R  (mod n)
            two = modadd64(unityval, unityval, n);
            four = modadd64(two, two, n);
            five = modadd64(unityval, four, n);
            eight = modadd64(four, four, n);
            sixteen = modadd64(eight, eight, n);
            two_8 = montmul64(sixteen, sixteen, n, rho);   // R*2^8          (mod n)
            two_16 = montmul64(two_8, two_8, n, rho);      // R*2^16         (mod n)
            two_32 = montmul64(two_16, two_16, n, rho);    // R*2^32         (mod n)
            Rsqr = montmul64(two_32, two_32, n, rho);      // R*2^64 == R*R  (mod n)
        }
    }

    return;
}


typedef struct
{
    uint32_t X[3];
    uint32_t Z[3];
} uecm96_pt;

__device__ void uaddxz96(uint32_t rho, uint32_t* n, uint32_t* p1x, uint32_t* p1z,
    uint32_t* p2x, uint32_t* p2z, uint32_t* pix, uint32_t* piz,
    uint32_t* pox, uint32_t* poz) {

    uint32 diff1[3], sum1[3], diff2[3], sum2[3];
    uint32 tt1[3], tt2[3], tt3[3], tt4[3];

    modsub96(p1x, p1z, diff1, n);
    modadd96(p1x, p1z, sum1, n);
    modsub96(p2x, p2z, diff2, n);
    modadd96(p2x, p2z, sum2, n);

    montmul96(diff1, sum2, tt1, n, rho);    // U
    montmul96(sum1, diff2, tt2, n, rho);    // V

    modadd96(tt1, tt2, tt3, n);
    modsub96(tt1, tt2, tt4, n);
    montmul96(tt3, tt3, tt1, n, rho);       // (U + V)^2
    montmul96(tt4, tt4, tt2, n, rho);       // (U - V)^2

    montmul96(tt1, piz, tt3, n, rho);       // Z * (U + V)^2
    montmul96(tt2, pix, tt4, n, rho);       // x * (U - V)^2

    pox[0] = tt3[0];
    pox[1] = tt3[1];
    pox[2] = tt3[2];
    poz[0] = tt4[0];
    poz[1] = tt4[1];
    poz[2] = tt4[2];

    return;
}

__device__ void uadd96(uint32 rho, uint32* n, uecm96_pt p1,
    uecm96_pt p2, uecm96_pt pi, uecm96_pt* po) {

    uaddxz96(rho, n, p1.X, p1.Z, p2.X, p2.Z, pi.X, pi.Z, po->X, po->Z);
    return;
}

__device__ void udup96(uint32* s, uint32 rho, uint32* n,
    uint32* insum, uint32* indiff, uecm96_pt* P) {

    uint32 tt1[3], tt2[3], tt3[3];

    montmul96(indiff, indiff, tt1, n, rho);     // U=(x1 - z1)^2
    montmul96(insum, insum, tt2, n, rho);       // V=(x1 + z1)^2
    montmul96(tt1, tt2, P->X, n, rho);          // x=U*V

    modsub96(tt2, tt1, tt3, n);                 // w = V-U
    montmul96(tt3, s, tt2, n, rho);             // w = (A+2)/4 * w
    modadd96(tt2, tt1, tt2, n);                 // w = w + U

    montmul96(tt2, tt3, P->Z, n, rho);          // Z = w*(V-U)
}

__device__ void swap96(uint32* a, uint32* b)
{
    uint32 t;
    t = a[0];
    a[0] = b[0];
    b[0] = t;

    t = a[1];
    a[1] = b[1];
    b[1] = t;

    t = a[2];
    a[2] = b[2];
    b[2] = t;
}

__device__ void threeswap96(uint32* a, uint32* b, uint32* c)
{
    uint32 t;
    t = a[0];
    a[0] = b[0];
    b[0] = c[0];
    c[0] = t;

    t = a[1];
    a[1] = b[1];
    b[1] = c[1];
    c[1] = t;

    t = a[2];
    a[2] = b[2];
    b[2] = c[2];
    c[2] = t;
}

__device__ void fourswap96(uint32* a, uint32* b, uint32* c, uint32* d)
{
    uint32 t;
    t = a[0];
    a[0] = b[0];
    b[0] = c[0];
    c[0] = d[0];
    d[0] = t;

    t = a[1];
    a[1] = b[1];
    b[1] = c[1];
    c[1] = d[1];
    d[1] = t;

    t = a[2];
    a[2] = b[2];
    b[2] = c[2];
    c[2] = d[2];
    d[2] = t;
}

__device__ void uprac96(uint32 rho, uint32* n, uecm96_pt* P,
    uint64_t c, double v, uint32* s)
{
    uint64_t d, e, r;

    // we require c != 0
    int shift = __ffsll(c) - 1;

    d = c;
    r = (uint64_t)((double)d * v + 0.5);

    d = c - r;
    e = 2 * r - c;

    uint32_t s1[3], s2[3], d1[3], d2[3];
    uecm96_pt pt1, pt2, pt3, pt4, pt5;

    // the first one is always a doubling
    // point1 is [1]P
    pt1.X[0] = pt2.X[0] = pt3.X[0] = P->X[0];
    pt1.Z[0] = pt2.Z[0] = pt3.Z[0] = P->Z[0];
    pt1.X[1] = pt2.X[1] = pt3.X[1] = P->X[1];
    pt1.Z[1] = pt2.Z[1] = pt3.Z[1] = P->Z[1];
    pt1.X[2] = pt2.X[2] = pt3.X[2] = P->X[2];
    pt1.Z[2] = pt2.Z[2] = pt3.Z[2] = P->Z[2];

    modsub96(pt1.X, pt1.Z, d1, n);
    modadd96(pt1.X, pt1.Z, s1, n);

    // point2 is [2]P
    udup96(s, rho, n, s1, d1, &pt1);

    while (d != e)
    {
        if (d < e)
        {
            r = d;
            d = e;
            e = r;
            swap96(pt1.X, pt2.X);
            swap96(pt1.Z, pt2.Z);
        }
        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        {
            d = (2 * d - e) / 3;
            e = (e - d) / 2;

            uadd96(rho, n, pt1, pt2, pt3, &pt4); // T = A + B (C)
            uadd96(rho, n, pt4, pt1, pt2, &pt5); // T2 = T + A (B)
            uadd96(rho, n, pt2, pt4, pt1, &pt2); // B = B + T (A)

            swap96(pt1.X, pt5.X);
            swap96(pt1.Z, pt5.Z);
        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        {
            d = (d - e) / 2;

            modsub96(pt1.X, pt1.Z, d1, n);
            modadd96(pt1.X, pt1.Z, s1, n);

            uadd96(rho, n, pt1, pt2, pt3, &pt2);        // B = A + B (C)
            udup96(s, rho, n, s1, d1, &pt1);        // A = 2A
        }
        else if ((d + 3) / 4 <= e)
        {
            d -= e;

            uadd96(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            threeswap96(pt2.X, pt4.X, pt3.X);
            threeswap96(pt2.Z, pt4.Z, pt3.Z);
        }
        else if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;

            modsub96(pt1.X, pt1.Z, d2, n);
            modadd96(pt1.X, pt1.Z, s2, n);

            uadd96(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            udup96(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 2 == 0)
        {
            d /= 2;

            modsub96(pt1.X, pt1.Z, d2, n);
            modadd96(pt1.X, pt1.Z, s2, n);

            uadd96(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)
            udup96(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (d % 3 == 0)
        {
            d = d / 3 - e;

            modsub96(pt1.X, pt1.Z, d1, n);
            modadd96(pt1.X, pt1.Z, s1, n);

            udup96(s, rho, n, s1, d1, &pt4);        // T = 2A
            uadd96(rho, n, pt1, pt2, pt3, &pt5);        // T2 = A + B (C)
            uadd96(rho, n, pt4, pt1, pt1, &pt1);        // A = T + A (A)
            uadd96(rho, n, pt4, pt5, pt3, &pt4);        // T = T + T2 (C)

            threeswap96(pt3.X, pt2.X, pt4.X);
            threeswap96(pt3.Z, pt2.Z, pt4.Z);
        }
        else if ((d + e) % 3 == 0)
        {
            d = (d - 2 * e) / 3;

            uadd96(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)

            modsub96(pt1.X, pt1.Z, d2, n);
            modadd96(pt1.X, pt1.Z, s2, n);
            uadd96(rho, n, pt4, pt1, pt2, &pt2);        // B = T + A (B)
            udup96(s, rho, n, s2, d2, &pt4);        // T = 2A
            uadd96(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else if ((d - e) % 3 == 0)
        {
            d = (d - e) / 3;

            uadd96(rho, n, pt1, pt2, pt3, &pt4);        // T = A + B (C)

            modsub96(pt1.X, pt1.Z, d2, n);
            modadd96(pt1.X, pt1.Z, s2, n);
            uadd96(rho, n, pt3, pt1, pt2, &pt3);        // C = C + A (B)

            swap96(pt2.X, pt4.X);
            swap96(pt2.Z, pt4.Z);

            udup96(s, rho, n, s2, d2, &pt4);        // T = 2A
            uadd96(rho, n, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else
        {
            e /= 2;

            modsub96(pt2.X, pt2.Z, d2, n);
            modadd96(pt2.X, pt2.Z, s2, n);

            uadd96(rho, n, pt3, pt2, pt1, &pt3);        // C = C + B (A)
            udup96(s, rho, n, s2, d2, &pt2);        // B = 2B
        }
    }

    uadd96(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)

    int i;
    for (i = 0; i < shift; i++)
    {
        modsub96(P->X, P->Z, d1, n);
        modadd96(P->X, P->Z, s1, n);
        udup96(s, rho, n, s1, d1, P);     // P = 2P
    }

    return;
}

__device__ void uprac96_85(uint32_t rho, uint32* n, uecm96_pt* P, uint32* s)
{
    uint32 s1[3], s2[3], d1[3], d2[3];
    int i;
    uint8_t steps[146] = {
        0,6,0,6,0,6,0,6,0,4,
        6,0,4,6,0,4,4,6,0,4,
        4,6,0,5,4,6,0,3,3,4,
        6,0,3,5,4,6,0,3,4,3,
        4,6,0,5,5,4,6,0,5,3,
        3,4,6,0,3,3,4,3,4,6,
        0,4,3,4,3,5,3,3,3,3,
        3,3,3,3,4,6,0,3,3,3,
        3,3,3,3,3,3,4,3,4,3,
        4,6,0,3,4,3,5,4,6,0,
        5,3,3,3,4,6,0,5,4,3,
        5,4,6,0,4,3,3,3,5,4,
        6,0,4,3,5,3,3,4,6,0,
        3,3,3,3,5,4,6,0,3,3,
        3,4,3,3,4,6 };

    uecm96_pt pt1, pt2, pt3, pt4;
    for (i = 0; i < 146; i++)
    {
        if (steps[i] == 0)
        {
            pt1.X[0] = pt2.X[0] = pt3.X[0] = P->X[0];
            pt1.Z[0] = pt2.Z[0] = pt3.Z[0] = P->Z[0];
            pt1.X[1] = pt2.X[1] = pt3.X[1] = P->X[1];
            pt1.Z[1] = pt2.Z[1] = pt3.Z[1] = P->Z[1];
            pt1.X[2] = pt2.X[2] = pt3.X[2] = P->X[2];
            pt1.Z[2] = pt2.Z[2] = pt3.Z[2] = P->Z[2];

            modsub96(pt1.X, pt1.Z, d1, n);
            modadd96(pt1.X, pt1.Z, s1, n);
            udup96(s, rho, n, s1, d1, &pt1);
        }
        else if (steps[i] == 3)
        {
            // integrate step 4 followed by swap(1,2)
            uadd96(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            fourswap96(pt1.X, pt4.X, pt3.X, pt2.X);
            fourswap96(pt1.Z, pt4.Z, pt3.Z, pt2.Z);
        }
        else if (steps[i] == 4)
        {
            uadd96(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            threeswap96(pt2.X, pt4.X, pt3.X);
            threeswap96(pt2.Z, pt4.Z, pt3.Z);
        }
        else if (steps[i] == 5)
        {
            modsub96(pt1.X, pt1.Z, d2, n);
            modadd96(pt1.X, pt1.Z, s2, n);

            uadd96(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            udup96(s, rho, n, s2, d2, &pt1);        // A = 2A
        }
        else if (steps[i] == 6)
        {
            uadd96(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)
        }
    }
}

__device__ void uecm96_build(uecm96_pt* P, uint32_t rho, uint32* n,
    uint32_t sigma, uint32* ps, uint32* five, uint32* Rsqr)
{
    uint32 t1[3], t2[3], t3[3], t4[3], t5[3];
    uint32 u[3], v[3], x[3], z[3];

    t1[0] = sigma;
    t1[1] = 0;
    t1[2] = 0;
    montmul96(t1, Rsqr, u, n, rho);  // to_monty(sigma)

    modadd96(u, u, v, n);
    modadd96(v, v, v, n);            // 4*sigma

    montmul96(u, u, u, n, rho);
    modsub96(u, five, u, n);           // sigma^2 - 5

    montmul96(u, u, t1, n, rho);
    montmul96(t1, u, x, n, rho);  // u^3, preserving u

    modadd96(v, v, t2, n);             // 2*v
    modadd96(t2, t2, t2, n);           // 4*v
    modadd96(t2, t2, t2, n);           // 8*v
    modadd96(t2, t2, t2, n);          // 16*v
    montmul96(t2, x, t5, n, rho);    // 16*u^3*v

    montmul96(v, v, t1, n, rho);
    montmul96(t1, v, z, n, rho);  // v^3, preserving v

    // compute parameter A
    modsub96(v, u, t1, n);           // (v - u)
    montmul96(t1, t1, t2, n, rho);
    montmul96(t2, t1, t4, n, rho);   // (v - u)^3

    modadd96(u, u, t1, n);           // 2u
    modadd96(t1, u, t2, n);          // 3u
    modadd96(v, t2, t3, n);          // 3u + v

    montmul96(t3, t4, t1, n, rho);   // a = (v-u)^3 * (3u + v)

    // u holds the denom (jeff note: isn't it t5 that has the denom?)
    // t1 holds the numer
    // accomplish the division by multiplying by the modular inverse
    t2[0] = 1;
    t2[1] = 0;
    t2[2] = 0;
    montmul96(t5, t2, t5, n, rho);   // take t5 out of monty rep
    
    uint32 gcd[3];
    modinv96(t5, n, t3, gcd);
    
    montmul96(t3, Rsqr, t3, n, rho); // to_monty(t3)
    montmul96(t3, t1, ps, n, rho);

    P->X[0] = x[0];
    P->X[1] = x[1];
    P->X[2] = x[2];
    P->Z[0] = z[0];
    P->Z[1] = z[1];
    P->Z[2] = z[2];

    return;
}

__device__ void uecm96_stage1(uint32_t rho, uint32* n, uecm96_pt* P,
    uint32_t stg1, uint32* s)
{
    uint32_t q;
    uint32 diff1[3], sum1[3];

    // handle the only even case
    q = 2;
    while (q < (stg1 * 4))  // jeff: multiplying by 4 improves perf ~1%
    {
        modsub96(P->X, P->Z, diff1, n);
        modadd96(P->X, P->Z, sum1, n);
        udup96(s, rho, n, sum1, diff1, P);
        q *= 2;
    }

    if (stg1 == 27)
    {
        uprac96(rho, n, P, 3, 0.61803398874989485, s);
        uprac96(rho, n, P, 3, 0.61803398874989485, s);
        uprac96(rho, n, P, 3, 0.61803398874989485, s);
        uprac96(rho, n, P, 5, 0.618033988749894903, s);
        uprac96(rho, n, P, 5, 0.618033988749894903, s);
        uprac96(rho, n, P, 7, 0.618033988749894903, s);
        uprac96(rho, n, P, 11, 0.580178728295464130, s);
        uprac96(rho, n, P, 13, 0.618033988749894903, s);
        uprac96(rho, n, P, 17, 0.618033988749894903, s);
        uprac96(rho, n, P, 19, 0.618033988749894903, s);
        uprac96(rho, n, P, 23, 0.522786351415446049, s);
    }
    else
    {
        uprac96_85(rho, n, P, s);

        if (stg1 == 85)
        {
            uprac96(rho, n, P, 61, 0.522786351415446049, s);
        }
        else
        {
            uprac96(rho, n, P, 5, 0.618033988749894903, s);
            uprac96(rho, n, P, 11, 0.580178728295464130, s);
            uprac96(rho, n, P, 89, 0.618033988749894903, s);
            uprac96(rho, n, P, 97, 0.723606797749978936, s);
            uprac96(rho, n, P, 101, 0.556250337855490828, s);
            uprac96(rho, n, P, 107, 0.580178728295464130, s);
            uprac96(rho, n, P, 109, 0.548409048446403258, s);
            uprac96(rho, n, P, 113, 0.618033988749894903, s);

            if (stg1 == 125)
            {
                // jeff: moved 61 to here
                uprac96(rho, n, P, 61, 0.522786351415446049, s);
                uprac96(rho, n, P, 103, 0.632839806088706269, s);
            }
            else
            {
                uprac96(rho, n, P, 7747, 0.552188778811121, s); // 61 x 127
                uprac96(rho, n, P, 131, 0.618033988749894903, s);
                uprac96(rho, n, P, 14111, 0.632839806088706, s);  // 103 x 137
                uprac96(rho, n, P, 20989, 0.620181980807415, s);  // 139 x 151
                uprac96(rho, n, P, 157, 0.640157392785047019, s);
                uprac96(rho, n, P, 163, 0.551390822543526449, s);

                if (stg1 == 165)
                {
                    uprac96(rho, n, P, 149, 0.580178728295464130, s);
                }
                else
                {
                    uprac96(rho, n, P, 13, 0.618033988749894903, s);
                    uprac96(rho, n, P, 167, 0.580178728295464130, s);
                    uprac96(rho, n, P, 173, 0.612429949509495031, s);
                    uprac96(rho, n, P, 179, 0.618033988749894903, s);
                    uprac96(rho, n, P, 181, 0.551390822543526449, s);
                    uprac96(rho, n, P, 191, 0.618033988749894903, s);
                    uprac96(rho, n, P, 193, 0.618033988749894903, s);
                    uprac96(rho, n, P, 29353, 0.580178728295464, s);  // 149 x 197
                    uprac96(rho, n, P, 199, 0.551390822543526449, s);
                }
            }
        }
    }

    return;
}

__device__ void gcd96(uint32_t* u, uint32_t* v, uint32_t* gcd)
{
#ifdef __GNUC__
    //uint64_t a, b, c;
    //a = u; b = v;
    //while (b != 0)
    //{
    //    c = a % b;
    //    a = b;
    //    b = c;
    //}
    //
    //return a;

    __uint128_t a, b, c;
    a = (__uint128_t)u[2];
    a <<= 32;
    a |= (__uint128_t)u[1];
    a <<= 32;
    a |= (__uint128_t)u[0];

    b = (__uint128_t)v[2];
    b <<= 32;
    b |= (__uint128_t)v[1];
    b <<= 32;
    b |= (__uint128_t)v[0];

    while (b != 0)
    {
        c = a % b;
        a = b;
        b = c;
    }

    gcd[0] = (uint32_t)a;
    a >>= 32;
    gcd[1] = (uint32_t)a;
    a >>= 32;
    gcd[2] = (uint32_t)a;

#else
    // need alternative 96-bit gcd on MSVC.


#endif

    return;
}

#if 1

__device__ void crossprodxz96(uint32_t* res, uint32_t *ptx, uint32_t *ptz,
    uint32_t *Pbx, uint32_t *Pbz, uint32_t *Pbprod,
    uint32_t *ptprod, uint32_t *acc, uint32_t rho, uint32_t *n)
{
    // accumulate the cross product  (zimmerman syntax).
        // page 342 in C&P
    uint32_t tt1[3];
    uint32_t tt2[3];
    uint32_t tt3[3];

    modsub96(ptx, Pbx, tt1, n);
    modadd96(ptz, Pbz, tt2, n);
    montmul96(tt1, tt2, tt3, n, rho);

    modadd96(tt3, Pbprod, tt1, n);
    modsub96(tt1, ptprod, tt2, n);

    montmul96(acc, tt2, res, n, rho);

    return;
}

__device__ void uecm96_stage2_D30(uecm96_pt* P, uint32_t rho, uint32_t *n,
    uint32_t B1, uint32_t B2, uint32_t *s, uint32_t *unityval, uint32_t *result)
{
    int b;
    int i, j, k;
    uecm96_pt Pa;
    uecm96_pt Pstep;
    uecm96_pt Pdiff;
    uecm96_pt pt5, pt6;

    //int idx = 1; // threadIdx.x;

    // Q = P = result of stage 1
    // Compute small differences: 1, 7, 11, 13, 17, 19, 23, 29
    //__shared__ uint64_t Pbprod[8 * 128];
    //__shared__ uint64_t PbX[8 * 128];   
    //__shared__ uint64_t PbZ[8 * 128];
    uint32_t Pbprod[8][3];
    uint32_t PbX[8][3];
    uint32_t PbZ[8][3];
    uint32_t diff1[3];
    uint32_t sum1[3];

    // common init
    {
        // [1]Q
        PbX[0][0] = P->X[0];
        PbZ[0][0] = P->Z[0];
        PbX[0][1] = P->X[1];
        PbZ[0][1] = P->Z[1];
        PbX[0][2] = P->X[2];
        PbZ[0][2] = P->Z[2];

        // [2]Q
        modsub96(P->X, P->Z, diff1, n);
        modadd96(P->X, P->Z, sum1, n);
        udup96(s, rho, n, sum1, diff1, &Pa);
        PbX[1][0] = Pa.X[0];
        PbZ[1][0] = Pa.Z[0];
        PbX[1][1] = Pa.X[1];
        PbZ[1][1] = Pa.Z[1];
        PbX[1][2] = Pa.X[2];
        PbZ[1][2] = Pa.Z[2];

        // [2]Q + [1]Q([1]Q) = [3]Q
        uaddxz96(rho, n,
            PbX[0], PbZ[0],
            PbX[1], PbZ[1],
            PbX[0], PbZ[0],
            PbX[3], PbZ[3]);  // <-- temporary

        // 2*[3]Q = [6]Q
        modsub96(PbX[3], PbZ[3], diff1, n);
        modadd96(PbX[3], PbZ[3], sum1, n);
        udup96(s, rho, n, sum1, diff1, &pt6);   // pt6 = [6]Q

        // [3]Q + [2]Q([1]Q) = [5]Q
        uaddxz96(rho, n, PbX[3], PbZ[3],
            PbX[1], PbZ[1],
            PbX[0], PbZ[0],
            pt5.X, pt5.Z);    // <-- pt5 = [5]Q

        // [6]Q + [5]Q([1]Q) = [11]Q
        uaddxz96(rho, n, pt6.X, pt6.Z,
            pt5.X, pt5.Z,
            PbX[0], PbZ[0],
            PbX[2], PbZ[2]);    // <-- [11]Q

        // [11]Q + [6]Q([5]Q) = [17]Q
        uaddxz96(rho, n,
            PbX[2], PbZ[2],
            pt6.X, pt6.Z, pt5.X, pt5.Z,
            PbX[4], PbZ[4]);    // <-- [17]Q

        // [17]Q + [6]Q([11]Q) = [23]Q
        uaddxz96(rho, n, PbX[4], PbZ[4],
            pt6.X, pt6.Z, PbX[2], PbZ[2],
            PbX[6], PbZ[6]);    // <-- [23]Q

        // [23]Q + [6]Q([17]Q) = [29]Q
        uaddxz96(rho, n,
            PbX[6], PbZ[6],
            pt6.X, pt6.Z, PbX[4], PbZ[4],
            PbX[7], PbZ[7]);    // <-- [29]Q

        // [6]Q + [1]Q([5]Q) = [7]Q
        uaddxz96(rho, n, pt6.X, pt6.Z,
            PbX[0], PbZ[0], pt5.X, pt5.Z,
            PbX[1], PbZ[1]);    // <-- [7]Q

        // [7]Q + [6]Q([1]Q) = [13]Q
        uaddxz96(rho, n,
            PbX[1], PbZ[1], pt6.X, pt6.Z,
            PbX[0], PbZ[0],
            PbX[3], PbZ[3]);    // <-- [13]Q

        // [13]Q + [6]Q([7]Q) = [19]Q
        uaddxz96(rho, n,
            PbX[3], PbZ[3], pt6.X, pt6.Z,
            PbX[1], PbZ[1],
            PbX[5], PbZ[5]);    // <-- [19]Q

        // 4*[1]Q = [4]Q
        modsub96(PbX[0], PbZ[0], diff1, n);
        modadd96(PbX[0], PbZ[0], sum1, n);
        udup96(s, rho, n, sum1, diff1, &pt5);   // pt5 = [2]Q

        modsub96(pt5.X, pt5.Z, diff1, n);
        modadd96(pt5.X, pt5.Z, sum1, n);
        udup96(s, rho, n, sum1, diff1, &pt5);   // pt5 = [4]Q

        // Pd = [2w]Q
        // [17]Q + [13]Q([4]Q) = [30]Q
        uaddxz96(rho, n, PbX[4], PbZ[4],
            PbX[3], PbZ[3], pt5.X, pt5.Z,
            Pdiff.X, Pdiff.Z);   // <-- [30]Q

        Pstep.X[0] = Pdiff.X[0];
        Pstep.Z[0] = Pdiff.Z[0];
        Pstep.X[1] = Pdiff.X[1];
        Pstep.Z[1] = Pdiff.Z[1];
        Pstep.X[2] = Pdiff.X[2];
        Pstep.Z[2] = Pdiff.Z[2];
        modsub96(Pstep.X, Pstep.Z, diff1, n);
        modadd96(Pstep.X, Pstep.Z, sum1, n);
        udup96(s, rho, n, sum1, diff1, &Pstep);   // Pstep = [60]Q

        Pdiff.X[0] = Pstep.X[0];
        Pdiff.Z[0] = Pstep.Z[0];
        Pdiff.X[1] = Pstep.X[1];
        Pdiff.Z[1] = Pstep.Z[1];
        Pdiff.X[2] = Pstep.X[2];
        Pdiff.Z[2] = Pstep.Z[2];

        Pa.X[0] = Pdiff.X[0];
        Pa.Z[0] = Pdiff.Z[0];
        Pa.X[1] = Pdiff.X[1];
        Pa.Z[1] = Pdiff.Z[1];
        Pa.X[2] = Pdiff.X[2];
        Pa.Z[2] = Pdiff.Z[2];
        modsub96(Pa.X, Pa.Z, diff1, n);
        modadd96(Pa.X, Pa.Z, sum1, n);
        udup96(s, rho, n, sum1, diff1, &Pa);   // Pa = [120]Q

        // Now we have Pa = 120, Pstep = 60 and Pdiff = 60.  ready for giant step.

        // make all of the Pbprod's
        montmul96(PbX[0], PbZ[0], Pbprod[0], n, rho);
        montmul96(PbX[1], PbZ[1], Pbprod[1], n, rho);
        montmul96(PbX[2], PbZ[2], Pbprod[2], n, rho);
        montmul96(PbX[3], PbZ[3], Pbprod[3], n, rho);
        montmul96(PbX[4], PbZ[4], Pbprod[4], n, rho);
        montmul96(PbX[5], PbZ[5], Pbprod[5], n, rho);
        montmul96(PbX[6], PbZ[6], Pbprod[6], n, rho);
        montmul96(PbX[7], PbZ[7], Pbprod[7], n, rho);

    }

    // __syncthreads();

    // ---------------------------------------------------------------------
    // to here is the same for any B1.
    // ---------------------------------------------------------------------

    // advance giant step to one step beyond B1.
    uint32_t aval = 120;
    while (aval < B1)
    {
        pt5.X[0] = Pa.X[0];
        pt5.Z[0] = Pa.Z[0];
        pt5.X[1] = Pa.X[1];
        pt5.Z[1] = Pa.Z[1];
        pt5.X[2] = Pa.X[2];
        pt5.Z[2] = Pa.Z[2];
        uaddxz96(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, Pa.X, Pa.Z);
        Pdiff.X[0] = pt5.X[0];
        Pdiff.Z[0] = pt5.Z[0];
        Pdiff.X[1] = pt5.X[1];
        Pdiff.Z[1] = pt5.Z[1];
        Pdiff.X[2] = pt5.X[2];
        Pdiff.Z[2] = pt5.Z[2];
        aval += 60;
    }

    //initialize Paprod
    uint32 Paprod[3];
    montmul96(Pa.X, Pa.Z, Paprod, n, rho);

    //initialize accumulator
    uint32 acc[3];
    acc[0] = unityval[0];
    acc[1] = unityval[1];
    acc[2] = unityval[2];
    result[0] = 0;
    result[1] = 0;
    result[2] = 0;

    while (aval < B2)
    {
        int j;
        uint32 tmp[3];

        for (j = 0; j < 8; j++)
        {
            crossprodxz96(tmp, Pa.X, Pa.Z, PbX[j], PbZ[j],
                Pbprod[j], Paprod, acc, rho, n);

            if ((tmp[0] == 0) && (tmp[1] == 0) && (tmp[2] == 0) &&
                (result[0] == 0) && (result[1] == 0) && (result[2] == 0))
            {
                result[0] = acc[0];
                result[1] = acc[1];
                result[2] = acc[2];
            }
            acc[0] = tmp[0];
            acc[1] = tmp[1];
            acc[2] = tmp[2];
        }

        // giant step - use the addition formula for ECM
        pt5.X[0] = Pa.X[0];
        pt5.Z[0] = Pa.Z[0];
        pt5.X[1] = Pa.X[1];
        pt5.Z[1] = Pa.Z[1];
        pt5.X[2] = Pa.X[2];
        pt5.Z[2] = Pa.Z[2];
        uaddxz96(rho, n, Pa.X, Pa.Z, Pstep.X, Pstep.Z, Pdiff.X, Pdiff.Z, Pa.X, Pa.Z);
        Pdiff.X[0] = pt5.X[0];
        Pdiff.Z[0] = pt5.Z[0];
        Pdiff.X[1] = pt5.X[1];
        Pdiff.Z[1] = pt5.Z[1];
        Pdiff.X[2] = pt5.X[2];
        Pdiff.Z[2] = pt5.Z[2];
        aval += 60;
        montmul96(Pa.X, Pa.Z, Paprod, n, rho);

    }

    if ((result[0] == 0) && (result[1] == 0) && (result[2] == 0))
    {
        result[0] = acc[0];
        result[1] = acc[1];
        result[2] = acc[2];
    }

    // __syncthreads();

    return;
}

#endif


__global__ void gbl_ecm96(int num, uint32_t* n_in, uint32_t* rho_in, uint32_t* one_in,
    uint32_t* rsq_in, uint32_t* sigma_in, uint32_t* f_out, 
    uint32_t stg1, uint32_t stg2, int curve)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
    int numf = 0;

    if (idx < num)
    {
        uint32_t rho = rho_in[idx];
        uint32_t* n = &n_in[idx * 3];
        uint32_t* unityval = &one_in[idx * 3];
        uint32_t* Rsqr = &rsq_in[idx * 3];

        uecm96_pt P;
        uint32_t s[3];

        //uint32_t smallsigma[8] = { 11, 61, 56, 81, 83, 7, 30, 51 };
        uint64_t Z[8] = { 85184, 14526784, 11239424, 34012224, 36594368, 21952, 1728000, 8489664 };
        uint64_t X[8] = { 1560896ULL, 51312965696ULL, 30693697091ULL, 281784327616ULL,
            326229015104ULL, 85184ULL, 716917375ULL, 17495004736ULL };
        uint64_t negt1[8] = { 146313216ULL, 476803160866816ULL, 236251574395731ULL, 4838810084806656ULL,
                              5902145938870272ULL, 655360ULL, 1305683671875ULL, 109380272541696ULL };
        uint64_t d[8] = { 1098870784ULL, 200325818077184ULL, 110006210374144ULL, 1460769954361344ULL,
            1732928528232448ULL, 38162432ULL, 1376481360000ULL, 57103695458304ULL };

        f_out[3 * idx + 0] = 1;
        f_out[3 * idx + 1] = 0;
        f_out[3 * idx + 2] = 0;

        //int lastcurve = curve + 2;
        //for (; curve < lastcurve; curve++)
        {
            uint32_t gcd[3];

            if (curve < 8)
            {
                //ubuild(&P, rho, &work, goodsigma[curve]); // sigma);
                //sigma = smallsigma[curve];
                // lookup point
                P.X[0] = X[curve] & 0xffffffff;
                P.Z[0] = Z[curve] & 0xffffffff;
                P.X[1] = X[curve] >> 32;
                P.Z[1] = Z[curve] >> 32;
                P.X[2] = 0;
                P.Z[2] = 0;

                // some computation left to do for S parameter for this 'n'
                uint32_t num[3];
                uint32_t neg[3];
                uint32_t dem[3];

                dem[0] = d[curve] & 0xffffffff;
                dem[1] = d[curve] >> 32;
                dem[2] = 0;
                neg[0] = negt1[curve] & 0xffffffff;
                neg[1] = negt1[curve] >> 32;
                neg[2] = 0;

                montmul96(neg, Rsqr, num, n, rho);     // to Monty rep.

                // The mulredc postcondition guarantees  num < n.
                // not true with modmul96: uses AMM currently.
                // do a modsub first to guarantee num < n.
                modsub96(num, n, num, n);
                possub96(n, num, num); // num = n - num;

                modinv96(dem, n, dem, gcd);

                montmul96(dem, Rsqr, dem, n, rho);              // to Monty rep.
                montmul96(num, dem, s, n, rho);
                montmul96(P.X, Rsqr, P.X, n, rho);              // to Monty rep.
                montmul96(P.Z, Rsqr, P.Z, n, rho);              // to Monty rep.
            }
            else
            {
                uint32_t two[3], four[3], five[3];

                modadd96(unityval, unityval, two, n);
                modadd96(two, two, four, n);
                modadd96(unityval, four, five, n);

                uecm96_build(&P, rho, n, sigma_in[idx], s, five, Rsqr);
            }

            //if (likely_gcd > 1)
            //{
            //    // If the gcd gave us a factor, we're done.  If not, since gcd != 1
            //    // the inverse calculated in uecm_build would have bogus, and so this
            //    // curve is probably set up for failure (hence we continue).
            //    if (likely_gcd < n) // || n % likely_gcd != 0)
            //    {
            //        if ((f_out[idx] == 1) || (f_out[idx] == n))
            //        {
            //            // if we haven't already found a factor, assign this result.
            //            f_out[idx] = likely_gcd;
            //        }
            //    }
            //}

            uecm96_stage1(rho, n, &P, stg1, s);

            uint32_t result[3];
            gcd96(n, P.Z, result);
            
            if (((f_out[3 * idx + 0] == 1) &&
                (f_out[3 * idx + 1] == 0) &&
                (f_out[3 * idx + 2] == 0)) ||
                ((f_out[3 * idx + 0] == n[0]) &&
                (f_out[3 * idx + 1] == n[1]) &&
                (f_out[3 * idx + 2] == n[2])))
            {
                // if we haven't already found a factor, assign this result.
                f_out[3 * idx + 0] = result[0];
                f_out[3 * idx + 1] = result[1];
                f_out[3 * idx + 2] = result[2];
            }

#if 1
            uint32_t stg2acc[3];
                
            uecm96_stage2_D30(&P, rho, n, stg1, stg2, s, unityval, stg2acc);

            gcd96(n, stg2acc, result);

            if (((f_out[3 * idx + 0] == 1) &&
                (f_out[3 * idx + 1] == 0) &&
                (f_out[3 * idx + 2] == 0)) ||
                ((f_out[3 * idx + 0] == n[0]) &&
                (f_out[3 * idx + 1] == n[1]) &&
                (f_out[3 * idx + 2] == n[2])))
            {
                // if we haven't already found a factor, assign this result.
                f_out[3 * idx + 0] = result[0];
                f_out[3 * idx + 1] = result[1];
                f_out[3 * idx + 2] = result[2];
            }
#endif

            // new curve
            uint64_t sigma64 = sigma_in[idx];
            sigma_in[idx] = lcg_rand_32b(7, (uint32_t)-1, &sigma64);

            //sigma_in[idx]++;
        }
    }

    return;
}


#ifdef __cplusplus
}
#endif
