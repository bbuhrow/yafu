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
    // require postive odd c
    uint64_t d, e, r;

    d = c;
    r = (uint64_t)((double)d * v + 0.5);
    d = c - r;
    e = 2 * r - c;

    uint64_t s1, d1;
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

__device__ uint64_t uecm_build(uecm_pt* P, uint32_t rho, uint64_t n,
    uint32_t sigma, uint64_t* ps, uint64_t one, uint64_t Rsqr)
{
    uint64_t t1, t2, t3, t4, t5;
    uint64_t u, v;

    t2 = modadd64(one, one, n);
    t4 = modadd64(t2, t2, n);
    t5 = modadd64(one, t4, n);

    u = montmul64((uint64_t)sigma, Rsqr, n, rho);  // to_monty(sigma)
    v = modadd64(u, u, n);
    v = modadd64(v, v, n);            // 4*sigma
    u = montmul64(u, u, n, rho);
    u = modsub64(u, t5, n);           // sigma^2 - 5
    t1 = montmul64(u, u, n, rho);
    uint64_t tmpx = montmul64(t1, u, n, rho);  // u^3

    t2 = modadd64(v, v, n);             // 2*v
    t2 = modadd64(t2, t2, n);           // 4*v
    t2 = modadd64(t2, t2, n);           // 8*v
    t2 = modadd64(t2, t2, n);           // 16*v
    t5 = montmul64(t2, tmpx, n, rho);    // 16*u^3*v
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

    // accomplish the division by multiplying by the modular inverse
    t2 = 1;
    t5 = montmul64(t5, t2, n, rho);   // take t5 out of monty rep
    t3 = modinv64(t5, n, (unsigned long long *)&t2);
    t3 = montmul64(t3, Rsqr, n, rho); // to_monty(t3)
    *ps = montmul64(t3, t1, n, rho);

    P->X = tmpx;
    P->Z = tmpz;

    return t2;
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




#if 0
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
#endif

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

#if 0
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
#endif

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

        f_out[idx] = 1;
        //int lastcurve = curve + 2;
        //for (; curve < lastcurve; curve++)
        {
            uint64_t gcd;
            gcd = uecm_build(&P, rho, n, sigma_in[idx], &s, unityval, Rsqr);

            if (gcd > 1)
            {
                // If the gcd gave us a factor, we're done.  If not, since gcd != 1
                // the inverse calculated in uecm_build would have bogus, and so this
                // curve is probably set up for failure (hence we continue).
                if (gcd < n) // || n % likely_gcd != 0)
                {
                    if ((f_out[idx] == 1) || (f_out[idx] == n))
                    {
                        // if we haven't already found a factor, assign this result.
                        f_out[idx] = gcd;
                    }
                }
            }

            uecm_stage1(rho, n, &P, stg1, s);

            gcd = gcd64(n, P.Z);

            if ((f_out[idx] == 1) || (f_out[idx] == n))
            {
                // if we haven't already found a factor, assign this result.
                f_out[idx] = gcd;
            }

            if (1)
            {
                uint64_t stg2acc = uecm_stage2_D30(&P, rho, n, stg1, s, unityval);

                gcd = gcd64(stg2acc, n);

                if ((f_out[idx] == 1) || (f_out[idx] == n))
                {
                    // if we haven't already found a factor, assign this result.
                    f_out[idx] = gcd;
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

typedef struct
{
    uint32_t X[3];
    uint32_t Z[3];
} uecm96_pt;

__device__ void uaddxz96(uint32_t rho, uint32_t* n, uint32_t* p1x, uint32_t* p1z,
    uint32_t* p2x, uint32_t* p2z, uint32_t* pix, uint32_t* piz,
    uint32_t* pox, uint32_t* poz) {

    uint32 diff1[3], sum1[3], diff2[3], sum2[3];
    uint32 tt1[3], tt2[3];

    //modsub96(p1x, p1z, diff1, n);
    //modadd96(p1x, p1z, sum1, n);
    //modsub96(p2x, p2z, diff2, n);
    //modadd96(p2x, p2z, sum2, n);

    modaddsub96(p1x, p1z, sum1, diff1, n);
    modaddsub96(p2x, p2z, sum2, diff2, n);

    montmul96(diff1, sum2, tt1, n, rho);    // U
    montmul96(sum1, diff2, tt2, n, rho);    // V

    //modadd96(tt1, tt2, tt3, n);
    //modsub96(tt1, tt2, tt4, n);
    modaddsub96(tt1, tt2, sum1, diff1, n);

    montsqr96(sum1, tt1, n, rho);       // (U + V)^2
    montsqr96(diff1, tt2, n, rho);       // (U - V)^2

    montmul96(tt1, piz, tt1, n, rho);       // Z * (U + V)^2
    montmul96(tt2, pix, tt2, n, rho);       // x * (U - V)^2

    pox[0] = tt1[0];
    pox[1] = tt1[1];
    pox[2] = tt1[2];
    poz[0] = tt2[0];
    poz[1] = tt2[1];
    poz[2] = tt2[2];

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

    montsqr96(indiff, tt1, n, rho);     // U=(x1 - z1)^2
    montsqr96(insum, tt2, n, rho);       // V=(x1 + z1)^2
    montmul96(tt1, tt2, P->X, n, rho);          // x=U*V

    modsub96(tt2, tt1, tt3, n);                 // w = V-U
    montmul96(tt3, s, tt2, n, rho);             // w = (A+2)/4 * w
    modadd96(tt2, tt1, tt2, n);                 // w = w + U

    montmul96(tt2, tt3, P->Z, n, rho);          // Z = w*(V-U)
}

__device__ void udup96as(uint32* s, uint32 rho, uint32* n, uecm96_pt* P) {

    uint32 tt1[3], tt2[3], tt3[3];
    uint32 insum[3], indiff[3];

    modsub96(P->X, P->Z, indiff, n);
    modadd96(P->X, P->Z, insum, n);

    montsqr96(indiff, tt1, n, rho);     // U=(x1 - z1)^2
    montsqr96(insum, tt2, n, rho);       // V=(x1 + z1)^2
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

__device__ void uprac96(uint32 rho, uint32* n, uecm96_pt* P,
    uint64_t c, double v, uint32* s)
{
    // require positive odd c
    uint64_t d, e, r;

    d = c;
    r = (uint64_t)((double)d * v + 0.5);
    d = c - r;
    e = 2 * r - c;

    uecm96_pt pt1, pt2, pt3, pt4, pt5;

    // the first one is always a doubling
    // point1 is [1]P
    pt1.X[0] = pt2.X[0] = pt3.X[0] = P->X[0];
    pt1.Z[0] = pt2.Z[0] = pt3.Z[0] = P->Z[0];
    pt1.X[1] = pt2.X[1] = pt3.X[1] = P->X[1];
    pt1.Z[1] = pt2.Z[1] = pt3.Z[1] = P->Z[1];
    pt1.X[2] = pt2.X[2] = pt3.X[2] = P->X[2];
    pt1.Z[2] = pt2.Z[2] = pt3.Z[2] = P->Z[2];

    // point2 is [2]P
    udup96as(s, rho, n, &pt1);

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
        // in our small-B1 cases these are the only
        // PRAC conditions used.  Need to verify that
        // continues to be the case if B1 grows.
        if ((d + 3) / 4 <= e)
        {
            d -= e;

            uadd96(rho, n, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            threeswap96(pt2.X, pt4.X, pt3.X);
            threeswap96(pt2.Z, pt4.Z, pt3.Z);
        }
        else
        {
            d = (d - e) / 2;

            uadd96(rho, n, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            udup96as(s, rho, n, &pt1);        // A = 2A
        }
    }

    uadd96(rho, n, pt1, pt2, pt3, P);     // A = A + B (C)

    return;
}

__device__ void uecm96_build(uecm96_pt* P, uint32_t rho, uint32* n,
    uint32_t sigma, uint32* ps, uint32* one, uint32* Rsqr, uint32_t *gcd)
{
    uint32 t1[3], t2[3], t3[3], t4[3], t5[3];
    uint32 u[3], v[3], x[3], z[3];

    modadd96(one, one, t2, n);
    modadd96(t2, t2, t4, n);
    modadd96(one, t4, t5, n);

    t1[0] = sigma;
    t1[1] = 0;
    t1[2] = 0;
    montmul96(t1, Rsqr, u, n, rho);  // to_monty(sigma)

    modadd96(u, u, v, n);
    modadd96(v, v, v, n);            // 4*sigma

    montmul96(u, u, u, n, rho);
    modsub96(u, t5, u, n);           // sigma^2 - 5

    montsqr96(u, t1, n, rho);
    montmul96(t1, u, x, n, rho);  // u^3, preserving u

    modadd96(v, v, t2, n);             // 2*v
    modadd96(t2, t2, t2, n);           // 4*v
    modadd96(t2, t2, t2, n);           // 8*v
    modadd96(t2, t2, t2, n);          // 16*v
    montmul96(t2, x, t5, n, rho);    // 16*u^3*v

    montsqr96(v, t1, n, rho);
    montmul96(t1, v, z, n, rho);  // v^3, preserving v

    // compute parameter A
    modsub96(v, u, t1, n);           // (v - u)
    montsqr96(t1, t2, n, rho);
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
    
    modinv96(t5, n, t3, gcd);
    
    montmul96(t3, Rsqr, t3, n, rho);    // to_monty(t3)
    montmul96(t3, t1, ps, n, rho);      // S = (v-u)^3 * (3u + v) / 16*u^3*v

    P->X[0] = x[0];
    P->X[1] = x[1];
    P->X[2] = x[2];
    P->Z[0] = z[0];
    P->Z[1] = z[1];
    P->Z[2] = z[2];

    return;
}

#if 0
__device__ void uecm96_build_param3(uecm96_pt* P, uint32_t rho, uint32* n,
    uint32_t sigma, uint32* ps, uint32* one, uint32* Rsqr, uint32_t* gcd)
{
    uint32 t1[4], nn[4], a[3], x[3], z[3];
    int i;

    t1[0] = sigma;
    t1[1] = 0;
    t1[2] = 0;
    t1[3] = 0;

    nn[0] = n[0];
    nn[1] = n[1];
    nn[2] = n[2];
    nn[3] = 0;
    
    /* A=4*d-2 with d = sigma/2^GMP_NUMB_BITS*/
    /* Compute d = sigma/2^GMP_NUMB_BITS */
    int carry = 0;
    for (i = 0; i < 32; i++)
    {
        carry = 0;
        if ((t1[0] & 1) == 1)
        {
            //carry = add96(t1, t1, n);
        }
        //rshift96_1(t1);
        if (carry)
        {
            t1[2] |= 0x80000000;
        }
    }

    //while (cmp96(t1, n) > 0)
    //{
    //    sub96(t1, t1, n);
    //}

    //lshift128(t1, t1, 2);
    //sub128_1(t1, t1, 2);

    //while (cmp128(t1, nn) > 0)
    //{
    //    sub128(t1, t1, nn);
    //}

    montmul96(t1, Rsqr, a, n, rho);    // a = to_monty(4*d-2)

    t1[0] = 2;
    t1[1] = 0;
    t1[2] = 0;

    montmul96(t1, Rsqr, x, n, rho);      // x = to_monty(2)

    P->X[0] = x[0];
    P->X[1] = x[1];
    P->X[2] = x[2];
    P->Z[0] = z[0];
    P->Z[1] = z[1];
    P->Z[2] = z[2];

    return;
}
#endif

__device__ void uecm96_stage1(uint32_t rho, uint32* n, uecm96_pt* P,
    uint32_t stg1, uint32* s)
{
    uint32_t q;

    // handle the only even case
    q = 2;
    while (q < (stg1))  // jeff: multiplying by 4 improves perf ~1%
    {
        udup96as(s, rho, n, P);
        q *= 2;
    }

    if (stg1 > 600)
    {
        // ./cuda_3lp -b1 700 -b2 100000 -m 4
        // not pair optimized
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

    if (stg1 > 500)
    {
        // ./cuda_3lp -b1 600 -b2 100000 -m 4
        // not pair optimized
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
        uprac96(rho, n, P, 5, 0.618033988749894903, s);
        uprac96(rho, n, P, 23, 0.522786351415446049, s);    // <-- pairs well
    }

    if (stg1 > 400)
    {
        // ./cuda_3lp -b1 500 -b2 100000 -m 4
        // anything greater than 400 gets B1=500
        // pair-optimized only within the set 400 < p < 500
        uprac96(rho, n, P, 184861, 0.628445063812485882, s);//[401,461] is 155.0, savings = 2
        uprac96(rho, n, P, 409, 0.551390822543526449, s);
        uprac96(rho, n, P, 419, 0.524531838023777786, s);
        uprac96(rho, n, P, 421, 0.643616254685781319, s);
        uprac96(rho, n, P, 431, 0.551390822543526449, s);
        uprac96(rho, n, P, 433, 0.553431118763021646, s);
        uprac96(rho, n, P, 439, 0.618033988749894903, s);
        uprac96(rho, n, P, 217513, 0.586747080461318182, s);//[443,491] is 156.0, savings = 1
        uprac96(rho, n, P, 449, 0.612429949509495031, s);
        uprac96(rho, n, P, 457, 0.612429949509495031, s);
        //uprac96(rho, n, P, 461, 0.520959819992697026, s);
        //[401,461] is 155.0, savings = 2
        //[43, 461] is 126.0, savings = 2
        //[367,461] is 156.0, savings = 4
        //[379,461] is 156.0, savings = 3
        uprac96(rho, n, P, 463, 0.553431118763021646, s);
        uprac96(rho, n, P, 467, 0.553431118763021646, s);
        uprac96(rho, n, P, 479, 0.618033988749894903, s);
        uprac96(rho, n, P, 487, 0.591965645556728037, s);
        //[431,491] is 155.0, savings = 1
        //[11,491] is 108.0, savings = 2
        //[443,491] is 156.0, savings = 1
        //[457,491] is 156.0, savings = 1
        //uprac96(rho, n, P, 491, 0.618347119656228017, s);
        uprac96(rho, n, P, 499, 0.618033988749894903, s);
    }

    if (stg1 > 300)
    {
        // anything greater than 300 gets B1=400
        // pair-optimized within the set 250 < p < 400
        uprac96(rho, n, P, 307, 0.580178728295464130, s);
        uprac96(rho, n, P, 311, 0.618033988749894903, s);
        uprac96(rho, n, P, 313, 0.618033988749894903, s);
        uprac96(rho, n, P, 317, 0.618033988749894903, s);
        //uprac96(rho, n, P, 331, 0.543797196679826511, s);
        uprac96(rho, n, P, 337, 0.618033988749894903, s);
        uprac96(rho, n, P, 124573, 0.552705982126542983, s); //[347,359], savings = 3
        uprac96(rho, n, P, 349, 0.632839806088706269, s);
        uprac96(rho, n, P, 95663, 0.618033988749894903, s);// [271,353], savings = 2
        //uprac96(rho, n, P, 359, 0.612429949509495031, s);
        uprac96(rho, n, P, 142763, 0.524817056539543136, s);// [367,389], savings = 1
        uprac96(rho, n, P, 373, 0.524531838023777786, s);
        uprac96(rho, n, P, 145157, 0.580178728295464130, s);//[379,383], savings = 2
        //uprac96(rho, n, P, 383, 0.537965503694917802, s);
        //uprac96(rho, n, P, 389, 0.632839806088706269, s);
        uprac96(rho, n, P, 397, 0.580178728295464130, s);
        // one more factor of 7 will fit - and it pairs with 331, savings 2
        uprac96(rho, n, P, 2317, 0.618033988749894903, s);  // 7 x 331, savings 2
        uprac96(rho, n, P, 19, 0.618033988749894903, s);
    }

    if (stg1 > 250)
    {
        // anything greater than 250 gets B1=300
        // pair-optimized within the set 200 < p < 300
        uprac96(rho, n, P, 251, 0.541554796058780874, s);
        //uprac96(rho, n, P, 257, 0.551390822543526449, s); // paired with 223 below
        uprac96(rho, n, P, 263, 0.612429949509495031, s);
        uprac96(rho, n, P, 269, 0.618033988749894903, s);
        if (stg1 > 300)
        {
            // paired with 353, above
        }
        else
        {
            uprac96(rho, n, P, 271, 0.618033988749894903, s);
        }
        uprac96(rho, n, P, 277, 0.618033988749894903, s);
        uprac96(rho, n, P, 281, 0.580178728295464130, s);
        uprac96(rho, n, P, 283, 0.580178728295464130, s);
        uprac96(rho, n, P, 293, 0.551390822543526449, s);
        uprac96(rho, n, P, 17, 0.618033988749894903, s);
    }

    if (stg1 > 200)
    {
        // anything greater than 200 gets B1=250
        // pair-optimized only within the set 200 < p < 250
        uprac96(rho, n, P, 48319, 0.552793637425581075, s);     // 211 x 229 savings 3
        if (stg1 > 250)
        {
            uprac96(rho, n, P, 57311, 0.552740754023503311, s); // 223 x 257 savings 3
            uprac96(rho, n, P, 241, 0.625306711365725132, s);   // [241,347], savings = 2
        }
        else
        {
            uprac96(rho, n, P, 53743, 0.718522647487825128, s); // 223 x 241 savings 2
        }
        uprac96(rho, n, P, 54253, 0.580178728295464130, s);     // 227 x 239 savings 2
        uprac96(rho, n, P, 233, 0.618033988749894903, s);
        // one more factor of 3 will fit
        uprac96(rho, n, P, 3, 0.618033988749894903, s);         
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

#if 0
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
#endif

    return;
}

#if 0
__device__ void uecm96_stage1_ladder(uint32_t rho, uint32* n, uecm96_pt* P,
    uint32_t stg1, uint32* s)
{
    // prac outperforms this, slightly, for B1=250, but for large B1 this
    // could be the way to go with the full scalar multiplier passed in.
    uecm96_pt R0;
    uecm96_pt R1;
    int startbit;
    int numwords;
    int pow2;
    uint32_t e[11];

    // R0 = P
    // R1 = [2]P
    R0.Z[0] = R1.Z[0] = P->Z[0];
    R0.Z[1] = R1.Z[1] = P->Z[1];
    R0.Z[2] = R1.Z[2] = P->Z[2];
    R0.X[0] = R1.X[0] = P->X[0];
    R0.X[1] = R1.X[1] = P->X[1];
    R0.X[2] = R1.X[2] = P->X[2];

    udup96as(s, rho, n, &R1);

    if (stg1 > 400)
    {
        // exponent for B1=500:
        // 32b08645 16c9a10a 28bd7eee 08c84eb3 
        // 35fe7448 8e2381be cf794bb0 81128015 
        // db8056fb 67bf4817 a181d61f
        // we start at bit t-2, so the top word
        // has its first set bit removed.
        e[10] = 0x12b08645;
        e[9] = 0x16c9a10a;
        e[8] = 0x28bd7eee;
        e[7] = 0x08c84eb3;
        e[6] = 0x35fe7448;
        e[5] = 0x8e2381be;
        e[4] = 0xcf794bb0;
        e[3] = 0x81128015;
        e[2] = 0xdb8056fb;
        e[1] = 0x67bf4817;
        e[0] = 0xa181d61f;

        startbit = 28;
        numwords = 11;
        pow2 = 8;
    }
    else if (stg1 > 200)
    {
        // exponent for B1=250:
        // 32b08645 16c9a10a 28bd7eee 08c84eb3 
        // 35fe7448 8e2381be cf794bb0 81128015 
        // db8056fb 67bf4817 a181d61f
        // we start at bit t-2, so the top word
        // has its first set bit removed.
        e[10] = 0x12b08645;
        e[9] = 0x16c9a10a;
        e[8] = 0x28bd7eee;
        e[7] = 0x08c84eb3;
        e[6] = 0x35fe7448;
        e[5] = 0x8e2381be;
        e[4] = 0xcf794bb0;
        e[3] = 0x81128015;
        e[2] = 0xdb8056fb;
        e[1] = 0x67bf4817;
        e[0] = 0xa181d61f;

        startbit = 28;
        numwords = 11;
        pow2 = 7;
    }

    int i;
    int j;

    for (j = startbit; j >= 0; j--)
    {
        if (((1 << j) & e[10]) == 0)
        {
            uaddxz96(rho, n, R1.X, R1.Z, R0.X, R0.Z, P->X, P->Z, R1.X, R1.Z);
            udup96as(s, rho, n, &R0);
        }
        else
        {
            uaddxz96(rho, n, R1.X, R1.Z, R0.X, R0.Z, P->X, P->Z, R0.X, R0.Z);
            udup96as(s, rho, n, &R1);
        }
    }

    for (i = numwords - 2; i >= 0; i--)
    {
        int j;
        for (j = 31; j >= 0; j--)
        {
            if (((1 << j) & e[i]) == 0)
            {
                uaddxz96(rho, n, R1.X, R1.Z, R0.X, R0.Z, P->X, P->Z, R1.X, R1.Z);
                udup96as(s, rho, n, &R0);
            }
            else
            {
                uaddxz96(rho, n, R1.X, R1.Z, R0.X, R0.Z, P->X, P->Z, R0.X, R0.Z);
                udup96as(s, rho, n, &R1);
            }
        }
    }

    // and the last 7 '0' bits
    for (i = 0; i < pow2; i++)
    {
        udup96as(s, rho, n, &R0);
    }

    P->Z[0] = R0.Z[0];
    P->Z[1] = R0.Z[1];
    P->Z[2] = R0.Z[2];
    P->X[0] = R0.X[0];
    P->X[1] = R0.X[1];
    P->X[2] = R0.X[2];
    
    return;
}
#endif

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

#define USE_D30

#ifdef USE_D30

__device__ void uecm96_stage2_D30(uecm96_pt* P, uint32_t rho, uint32_t *n,
    uint32_t B1, uint32_t B2, uint32_t *s, uint32_t *unityval, uint32_t *result)
{
    int b;
    int i, j, k;
    uecm96_pt Pa;
    uecm96_pt Pstep;
    uecm96_pt Pdiff;
    uecm96_pt pt5, pt6;

    // Q = P = result of stage 1
    // With D=30, compute small differences: 1, 7, 11, 13, 17, 19, 23, 29.
    // at 96-bit modulus size these can be put into registers.
    // Does it help or hurt to use a larger D?  D=60 requires more precomputation
    // and more storage, perhaps pushing the function beyond what can be
    // accomodated with registers or reducing the number of threads in the
    // block, either of which are likely a net reduction in performance.
    // D=60 would need to compute and store info for 7 more primes:
    // 31, 37, 41, 43, 47, 53, 59.
    // for B1=500, B2=50*500, D=60 would reduce the total number of multiplications
    // in the main loop by about 1800:
    // (25000 - 500) / 120 = 208, (25000 - 500) / 60 = 408 giant steps.
    // 408*7+408*8*2 = 9384
    // 204*7+204*15*2 = 7548
    // The total stage 1 cost of B1 = 500 is 6047 muls, so already with B2=50 we
    // are greater than the cost of stage 1.  D=60 brings it closer.
    uint32_t Pbprod[8][3];
    uint32_t PbX[8][3];
    uint32_t PbZ[8][3];
    uint32_t diff1[3];
    uint32_t sum1[3];

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
        uaddxz96(rho, n, 
            PbX[3], PbZ[3],
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

    while (aval < B2)
    {
        int j;
        uint32 tmp[3];

        for (j = 0; j < 8; j++)
        {
            // cross product
            modsub96(Pa.X, PbX[j], diff1, n);
            modadd96(Pa.Z, PbZ[j], sum1, n);
            montmul96(diff1, sum1, tmp, n, rho);

            modadd96(tmp, Pbprod[j], sum1, n);
            modsub96(sum1, Paprod, diff1, n);
            montmul96(acc, diff1, acc, n, rho);
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

    result[0] = acc[0];
    result[1] = acc[1];
    result[2] = acc[2];

    return;
}

#elif defined(USE_D60)

__device__ void uecm96_stage2_D60(uecm96_pt* P, uint32_t rho, uint32_t* n,
    uint32_t B1, uint32_t B2, uint32_t* s, uint32_t* unityval, uint32_t* result)
{
    int b;
    int i, j, k;
    uecm96_pt Pa;
    uecm96_pt Pstep;
    uecm96_pt Pdiff;
    uecm96_pt pt5, pt6;

    // Q = P = result of stage 1
    // With D=30, compute small differences: 1, 7, 11, 13, 17, 19, 23, 29.
    // at 96-bit modulus size these can be put into registers.
    // Does it help or hurt to use a larger D?  D=60 requires more precomputation
    // and more storage, perhaps pushing the function beyond what can be
    // accomodated with registers or reducing the number of threads in the
    // block, either of which are likely a net reduction in performance.
    // UPDATE: on the A100 we have to go from 384 threads/block to <= 256 threads/block,
    // with 128 threads/block slightly faster than 256.
    // D=60 would need to compute and store info for 7 more primes:
    // 31, 37, 41, 43, 47, 53, 59.
    // for B1=500, B2=50*500, D=60 would reduce the total number of multiplications
    // in the main loop by about 1800:
    // (25000 - 500) / 120 = 208, (25000 - 500) / 60 = 408 giant steps.
    // 408*7+408*8*2 = 9384
    // 204*7+204*15*2 = 7548
    // The total stage 1 cost of B1 = 500 is 6047 muls, so already with B2=50 we
    // are greater than the cost of stage 1.  D=60 brings it closer.
    uint32_t Pbprod[15][3];
    uint32_t PbX[15][3];
    uint32_t PbZ[15][3];
    uint32_t diff1[3];
    uint32_t sum1[3];

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

        // ------------------------------------------------------------
        // 5 progression
        // ------------------------------------------------------------
        // [3]Q + [2]Q([1]Q) = [5]Q
        uaddxz96(rho, n,
            PbX[3], PbZ[3],
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

        // [29]Q + [6]Q([23]Q) = [35]Q
        uaddxz96(rho, n,
            PbX[7], PbZ[7],
            pt6.X, pt6.Z, PbX[6], PbZ[6],
            Pstep.X, Pstep.Z);    // <-- [35]Q (temporary)

        // [35]Q + [6]Q([29]Q) = [41]Q
        uaddxz96(rho, n,
            Pstep.X, Pstep.Z,
            pt6.X, pt6.Z, PbX[7], PbZ[7],
            PbX[8], PbZ[8]);    // <-- [41]Q

        // [41]Q + [6]Q([35]Q) = [47]Q
        uaddxz96(rho, n,
            PbX[8], PbZ[8],
            pt6.X, pt6.Z, Pstep.X, Pstep.Z,
            PbX[9], PbZ[9]);    // <-- [47]Q

        // [47]Q + [6]Q([41]Q) = [53]Q
        uaddxz96(rho, n,
            PbX[9], PbZ[9],
            pt6.X, pt6.Z, PbX[8], PbZ[8],
            PbX[10], PbZ[10]);    // <-- [53]Q

        // [53]Q + [6]Q([47]Q) = [59]Q
        uaddxz96(rho, n,
            PbX[10], PbZ[10],
            pt6.X, pt6.Z, PbX[9], PbZ[9],
            PbX[11], PbZ[11]);    // <-- [59]Q

        // ------------------------------------------------------------
        // 1 progression
        // ------------------------------------------------------------

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

        // [19]Q + [6]Q([13]Q) = [25]Q
        uaddxz96(rho, n,
            PbX[5], PbZ[5], pt6.X, pt6.Z,
            PbX[3], PbZ[3],
            Pstep.X, Pstep.Z);    // <-- [25]Q (temporary)

        // [25]Q + [6]Q([19]Q) = [31]Q
        uaddxz96(rho, n,
            Pstep.X, Pstep.Z, pt6.X, pt6.Z,
            PbX[5], PbZ[5],
            PbX[12], PbZ[12]);    // <-- [31]Q

        // [31]Q + [6]Q([25]Q) = [37]Q
        uaddxz96(rho, n,
            PbX[12], PbZ[12], pt6.X, pt6.Z,
            Pstep.X, Pstep.Z,
            PbX[13], PbZ[13]);    // <-- [37]Q

        // [37]Q + [6]Q([31]Q) = [43]Q
        uaddxz96(rho, n,
            PbX[13], PbZ[13], pt6.X, pt6.Z,
            PbX[12], PbZ[12],
            PbX[14], PbZ[14]);    // <-- [43]Q

        // Pd = [2w]Q
        // [31]Q + [29]Q([2]Q) = [60]Q
        uaddxz96(rho, n, PbX[12], PbZ[12],
            PbX[7], PbZ[7], Pa.X, Pa.Z,
            Pstep.X, Pstep.Z);   // <-- [60]Q

        modsub96(Pstep.X, Pstep.Z, diff1, n);
        modadd96(Pstep.X, Pstep.Z, sum1, n);
        udup96(s, rho, n, sum1, diff1, &Pstep);   // Pstep = [120]Q

        Pdiff.X[0] = Pstep.X[0];
        Pdiff.Z[0] = Pstep.Z[0];
        Pdiff.X[1] = Pstep.X[1];
        Pdiff.Z[1] = Pstep.Z[1];
        Pdiff.X[2] = Pstep.X[2];
        Pdiff.Z[2] = Pstep.Z[2];

        modsub96(Pdiff.X, Pdiff.Z, diff1, n);
        modadd96(Pdiff.X, Pdiff.Z, sum1, n);
        udup96(s, rho, n, sum1, diff1, &Pa);   // Pa = [240]Q

        // Now we have Pa = 240, Pstep = 120 and Pdiff = 120.  ready for giant step.

        // make all of the Pbprod's
        montmul96(PbX[0], PbZ[0], Pbprod[0], n, rho);
        montmul96(PbX[1], PbZ[1], Pbprod[1], n, rho);
        montmul96(PbX[2], PbZ[2], Pbprod[2], n, rho);
        montmul96(PbX[3], PbZ[3], Pbprod[3], n, rho);
        montmul96(PbX[4], PbZ[4], Pbprod[4], n, rho);
        montmul96(PbX[5], PbZ[5], Pbprod[5], n, rho);
        montmul96(PbX[6], PbZ[6], Pbprod[6], n, rho);
        montmul96(PbX[7], PbZ[7], Pbprod[7], n, rho);
        montmul96(PbX[8], PbZ[8], Pbprod[8], n, rho);
        montmul96(PbX[9], PbZ[9], Pbprod[9], n, rho);
        montmul96(PbX[10], PbZ[10], Pbprod[10], n, rho);
        montmul96(PbX[11], PbZ[11], Pbprod[11], n, rho);
        montmul96(PbX[12], PbZ[12], Pbprod[12], n, rho);
        montmul96(PbX[13], PbZ[13], Pbprod[13], n, rho);
        montmul96(PbX[14], PbZ[14], Pbprod[14], n, rho);
    }

    // advance giant step to one step beyond B1.
    uint32_t aval = 240;
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
        aval += 120;
    }

    //initialize Paprod
    uint32 Paprod[3];
    montmul96(Pa.X, Pa.Z, Paprod, n, rho);

    //initialize accumulator
    uint32 acc[3];
    acc[0] = unityval[0];
    acc[1] = unityval[1];
    acc[2] = unityval[2];

    // commence main loop
    while (aval < B2)
    {
        int j;
        uint32 tmp[3];

        for (j = 0; j < 15; j++)
        {
            // cross product
            modsub96(Pa.X, PbX[j], diff1, n);
            modadd96(Pa.Z, PbZ[j], sum1, n);
            montmul96(diff1, sum1, tmp, n, rho);

            modadd96(tmp, Pbprod[j], sum1, n);
            modsub96(sum1, Paprod, diff1, n);
            montmul96(acc, diff1, acc, n, rho);
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
        aval += 120;
        montmul96(Pa.X, Pa.Z, Paprod, n, rho);
    }

    result[0] = acc[0];
    result[1] = acc[1];
    result[2] = acc[2];

    return;
}

#else


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
        uint32_t n[3];
        uint32_t unityval[3];
        uint32_t Rsqr[3];
        uecm96_pt P;
        uint32_t s[3];

        n[0] = n_in[idx * 3 + 0];
        n[1] = n_in[idx * 3 + 1];
        n[2] = n_in[idx * 3 + 2];
        unityval[0] = one_in[idx * 3 + 0];
        unityval[1] = one_in[idx * 3 + 1];
        unityval[2] = one_in[idx * 3 + 2];
        Rsqr[0] = rsq_in[idx * 3 + 0];
        Rsqr[1] = rsq_in[idx * 3 + 1];
        Rsqr[2] = rsq_in[idx * 3 + 2];
        f_out[3 * idx + 0] = 1;
        f_out[3 * idx + 1] = 0;
        f_out[3 * idx + 2] = 0;

        {
            uint32_t gcd[3];

            uecm96_build(&P, rho, n, sigma_in[idx], s, unityval, Rsqr, gcd);

            if (gcd[0] > 1)
            {
                // If the gcd gave us a factor, we're done.  We don't have a 
                // good way of ignoring a bad inverse, so if the gcd is n
                // just continue.
                if ((gcd[0] != n[0]) || (gcd[1] != n[1]) || (gcd[2] != n[2]))
                {
                    if (((f_out[3 * idx + 0] == 1) &&
                        (f_out[3 * idx + 1] == 0) &&
                        (f_out[3 * idx + 2] == 0)) ||
                        ((f_out[3 * idx + 0] == n[0]) &&
                        (f_out[3 * idx + 1] == n[1]) &&
                        (f_out[3 * idx + 2] == n[2])))
                    {
                        // if we haven't already found a factor, assign this result.
                        f_out[3 * idx + 0] = gcd[0];
                        f_out[3 * idx + 1] = gcd[1];
                        f_out[3 * idx + 2] = gcd[2];
                    }
                }
            }

            //uecm96_stage1_ladder(rho, n, &P, stg1, s);
            uecm96_stage1(rho, n, &P, stg1, s);

            gcd96(n, P.Z, gcd);
            
            if (((f_out[3 * idx + 0] == 1) &&
                (f_out[3 * idx + 1] == 0) &&
                (f_out[3 * idx + 2] == 0)) ||
                ((f_out[3 * idx + 0] == n[0]) &&
                (f_out[3 * idx + 1] == n[1]) &&
                (f_out[3 * idx + 2] == n[2])))
            {
                // if we haven't already found a factor, assign this result.
                f_out[3 * idx + 0] = gcd[0];
                f_out[3 * idx + 1] = gcd[1];
                f_out[3 * idx + 2] = gcd[2];
            }

            uint32_t stg2acc[3];
#ifdef USE_D30
            uecm96_stage2_D30(&P, rho, n, stg1, stg2, s, unityval, stg2acc);
#elif defined(USE_D60)
            uecm96_stage2_D60(&P, rho, n, stg1, stg2, s, unityval, stg2acc);
#else
            stg2acc[0] = 1;
            stg2acc[1] = 0;
            stg2acc[2] = 0;
#endif
            
            gcd96(n, stg2acc, gcd);

            if (((f_out[3 * idx + 0] == 1) &&
                (f_out[3 * idx + 1] == 0) &&
                (f_out[3 * idx + 2] == 0)) ||
                ((f_out[3 * idx + 0] == n[0]) &&
                (f_out[3 * idx + 1] == n[1]) &&
                (f_out[3 * idx + 2] == n[2])))
            {
                // if we haven't already found a factor, assign this result.
                f_out[3 * idx + 0] = gcd[0];
                f_out[3 * idx + 1] = gcd[1];
                f_out[3 * idx + 2] = gcd[2];
            }

            // new curve
            uint64_t sigma64 = sigma_in[idx];
            sigma_in[idx] = lcg_rand_32b(7, (uint32_t)-1, &sigma64);
        }
    }

    return;
}


// cce3c7d6c0f418730ed656fd0b0ed388c0
__constant__ int tpm1_ewin100[34] = { // 170 muls
        12, 12, 14, 3, 12, 7, 13, 6, 12, 0, 15, 4, 1, 8, 7, 3, 0, 14, 13, 6,
        5, 6, 15, 13, 0, 11, 0, 14, 13, 3, 8, 8, 12, 0 };


__constant__ int tpm1_ewin333[119] = { // 595 muls
    2, 5, 1, 6, 12, 5, 1, 5, 0, 15, 3, 13, 2, 4, 2, 0, 4, 13, 11, 9, 4, 5,
    4, 13, 7, 15, 0, 11, 10, 7, 5, 4, 7, 0, 14, 11, 12, 10, 12, 4, 11, 2,
    5, 2, 10, 10, 7, 3, 14, 11, 8, 0, 15, 2, 2, 3, 10, 11, 6, 9, 2, 8, 15,
    4, 12, 13, 14, 13, 0, 7, 3, 3, 12, 3, 9, 8, 4, 6, 15, 0, 3, 9, 11, 14,
    5, 7, 4, 4, 11, 14, 8, 6, 11, 14, 0, 12, 13, 4, 12, 12, 11, 6, 13, 0,
    10, 8, 13, 6, 13, 15, 2, 5, 14, 13, 11, 11, 15, 0, 0 };

__constant__ uint32_t tpm1_ewin500[23] = { 
    0xd474c2d4,
    0xb7330cfe,
    0xb00f3a15,
    0x74d81ab8,
    0x23102ec9,
    0x1693c48d,
    0x9845eb35,
    0x9c0da860,
    0x9477df49,
    0x598d1d83,
    0x8c8bb315,
    0xa5add55b,
    0xb193f4f7,
    0x90e6ec89,
    0x2e477998,
    0x6fefd0b8,
    0xce5273a3,
    0x8952ee05,
    0x0f1dc5dd,
    0xa487cf80,
    0x484dd0af,
    0x21644a26,
    0xfd100000 };
    
__device__ void pm196_stage1(uint32_t* P, uint32_t* N, 
    uint32_t rho, uint32_t* one, uint32_t stg1)
{
    int i;
    uint32_t g[16][3];

    g[0][0] = P[0] = one[0];
    g[0][1] = P[1] = one[1];
    g[0][2] = P[2] = one[2];

    for (i = 1; i < 16; i++)
    {
        modadd96(g[i - 1], g[i - 1], g[i], N);
    }

    switch (stg1)
    {
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
            for (j = 7; j >= 0; j--)
            {
                montsqr96(P, P, N, rho);
                montsqr96(P, P, N, rho);
                montsqr96(P, P, N, rho);
                montsqr96(P, P, N, rho);
                montmul96(P, g[(tpm1_ewin500[i] >> (j * 4)) & 0xf], P, N, rho);
            }
        }
        break;
    }
    return;
}


__constant__ uint32_t tpm1_map[60] = {
    0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
    0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
    0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
    0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
    0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
    0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

__device__ void tpm1_expRL(uint32_t* P, uint32_t* in, uint32_t* n, uint32_t* one, 
    uint32_t m, uint32_t rho)
{
    uint32_t s[3];
    P[0] = one[0];
    P[1] = one[1];
    P[2] = one[2];
    s[0] = in[0];
    s[1] = in[1];
    s[2] = in[2];

    while (m > 0)
    {
        if (m & 1)
            montmul96(P, s, P, n, rho);
        montsqr96(s, s, n, rho);
        m >>= 1;
    }
    return;
}

__device__ void pm196_stage2_pair(uint32_t* P, uint32_t* acc, uint32_t *n, 
    uint32_t rho, uint32_t b1)
{
    int w = 60;
    uint32_t d[18][3], six[3], x12[3], xmid[3], x24[3], x25[3];
    uint32_t x36[3], x60[3], x72[3], five[3], pw[3], pgiant[3];
    int i, j;
    uint32_t b2 = 50 * b1;

    // we accumulate f(vw) - f(u) where f(n) = n^2, so that
    // we can pair together primes vw+/-u
    // see: P. L. Montgomery, "Speeding the Pollard and Elliptic Curve Methods
    // of Factorization," Mathematics of Computation, Vol 48, No. 177, 1987
    // 
    // b^(n+h)^2 = b^(n^2 + 2h(n) + h^2) = (b^(n^2))*(b^(2hn))*(b^(h^2))

    // u=1, u^2=1, b^u^2 = P
    d[1][0] = P[0];
    d[1][1] = P[1];
    d[1][2] = P[2];
    montsqr96(P, d[2], n, rho);
    montmul96(P, d[2], six, n, rho);
    montsqr96(six, six, n, rho);                // P^6
    montsqr96(six, x12, n, rho);                // P^12
    montsqr96(x12, x24, n, rho);                // P^24
    montmul96(x24, x12, x36, n, rho);           // P^36
    montmul96(x24, x36, x60, n, rho);           // P^60
    montsqr96(x36, x72, n, rho);                // P^72
    montsqr96(d[2], five, n, rho);
    montmul96(five, d[1], five, n, rho);        // P^5
    montmul96(x24, d[1], x25, n, rho);          // P^25


    // P^7^2 = P^(1+6)^2 = P^(1^2) * P^(12*1) * P^36
    // P^13^2 = P^(7+6)^2 = P^(7^2) * P^(12*7) * P^36
    // P^19^2 = P^(13+6)^2 = P^(13^2) * P^(12*13) * P^36
    // ...

    // 1, 7, 13, 19, 25, 31, 37, 43, 49
    // unnecessary powers will be mapped to scratch d[0].
    j = 1;
    //copy128(x12, xmid);
    xmid[0] = x12[0];
    xmid[1] = x12[1];
    xmid[2] = x12[2];
    while ((j + 6) < 50)
    {
        montmul96(d[tpm1_map[j]], x36, d[tpm1_map[j + 6]], n, rho);
        montmul96(d[tpm1_map[j + 6]], xmid, d[tpm1_map[j + 6]], n, rho);
        montmul96(xmid, x72, xmid, n, rho);
        j += 6;
    }

    // P^11^2 = P^(5+6)^2 = P^(5^2) * P^(12*5) * P^36
    // P^17^2 = P^(11+6)^2 = P^(11^2) * P^(12*11) * P^36
    // P^23^2 = P^(17+6)^2 = P^(17^2) * P^(12*17) * P^36
    // ...

    // 11, 17, 23, 29, 35, 41, 47, 53, 59
    // unnecessary powers will be mapped to scratch d[0].
    montmul96(x25, x36, d[tpm1_map[11]], n, rho);
    montmul96(d[tpm1_map[11]], x60, d[tpm1_map[11]], n, rho);
    montmul96(x60, x72, xmid, n, rho);
    j = 11;
    while ((j + 6) < 60)
    {
        montmul96(d[tpm1_map[j]], x36, d[tpm1_map[j + 6]], n, rho);
        montmul96(d[tpm1_map[j + 6]], xmid, d[tpm1_map[j + 6]], n, rho);
        montmul96(xmid, x72, xmid, n, rho);
        j += 6;
    }

    // P^(2w)^2, assumes w=60
    // P^(120^2) = P^(14440)
    tpm1_expRL(pw, P, n, acc, 120 * 120, rho);

    uint32_t x14400[3];

    //copy128(pw, x14400);
    x14400[0] = pw[0];
    x14400[1] = pw[1];
    x14400[2] = pw[2];

    uint32_t x240x120[3];
    montsqr96(pw, x240x120, n, rho);

    //copy128(x240x120, xmid);
    xmid[0] = x240x120[0];
    xmid[1] = x240x120[1];
    xmid[2] = x240x120[2];

    // P^(240^2) = P^(120+120)^2 = P^(120^2) * P^(240*120) * P^(14400)
    // P^(360^2) = P^(240+120)^2 = P^(240^2) * P^(240*240) * P^(14400)
    // P^(480^2) = P^(360+120)^2 = P^(360^2) * P^(240*360) * P^(14400)
    // ...

    //copy128(pw, pgiant);
    pgiant[0] = pw[0];
    pgiant[1] = pw[1];
    pgiant[2] = pw[2];

    i = 2 * w;
    while (i < b1)
    {
        montmul96(pgiant, x14400, pgiant, n, rho);
        montmul96(pgiant, xmid, pgiant, n, rho);
        montmul96(xmid, x240x120, xmid, n, rho);
        i += 2 * w;
    }

    // acc is monty(one) on input
    while (i < b2)
    {
        uint32_t tmp[3], sub[3];
        int j;

        modsub96(pgiant, d[1], sub, n);
        montmul96(acc, sub, acc, n, rho);

        for (j = 3; j < 18; j++)
        {
            modsub96(pgiant, d[j], sub, n);
            montmul96(acc, sub, acc, n, rho);
        }

        montmul96(pgiant, x14400, pgiant, n, rho);
        montmul96(pgiant, xmid, pgiant, n, rho);
        montmul96(xmid, x240x120, xmid, n, rho);

        i += 120;
    }

    return;
}

__global__ void gbl_pm196(int num, uint32_t* n_in, uint32_t* rho_in, uint32_t* one_in, 
    uint32_t* f_out, uint32_t B1, uint32_t B2)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
    int numf = 0;

    if (idx < num)
    {
        uint32_t rho = rho_in[idx];
        uint32_t n[3];
        uint32_t unityval[3];
        uint32_t P[3];

        n[0] = n_in[idx * 3 + 0];
        n[1] = n_in[idx * 3 + 1];
        n[2] = n_in[idx * 3 + 2];
        unityval[0] = one_in[idx * 3 + 0];
        unityval[1] = one_in[idx * 3 + 1];
        unityval[2] = one_in[idx * 3 + 2];
        f_out[3 * idx + 0] = 1;
        f_out[3 * idx + 1] = 0;
        f_out[3 * idx + 2] = 0;

        uint32_t acc[3];
        pm196_stage1(P, n, rho, unityval, B1);
        modsub96(P, unityval, acc, n);  // test stage1 P-1

        uint32_t gcd[3];

        gcd96(n, acc, gcd);

        if (((f_out[3 * idx + 0] == 1) &&
            (f_out[3 * idx + 1] == 0) &&
            (f_out[3 * idx + 2] == 0)) ||
            ((f_out[3 * idx + 0] == n[0]) &&
            (f_out[3 * idx + 1] == n[1]) &&
            (f_out[3 * idx + 2] == n[2])))
        {
            // if we haven't already found a factor, assign this result.
            f_out[3 * idx + 0] = gcd[0];
            f_out[3 * idx + 1] = gcd[1];
            f_out[3 * idx + 2] = gcd[2];
        }

        acc[0] = unityval[0];
        acc[1] = unityval[1];
        acc[2] = unityval[2];

        pm196_stage2_pair(P, acc, n, rho, B1);

        gcd96(n, acc, gcd);

        if (((f_out[3 * idx + 0] == 1) &&
            (f_out[3 * idx + 1] == 0) &&
            (f_out[3 * idx + 2] == 0)) ||
            ((f_out[3 * idx + 0] == n[0]) &&
            (f_out[3 * idx + 1] == n[1]) &&
            (f_out[3 * idx + 2] == n[2])))
        {
            // if we haven't already found a factor, assign this result.
            f_out[3 * idx + 0] = gcd[0];
            f_out[3 * idx + 1] = gcd[1];
            f_out[3 * idx + 2] = gcd[2];
        }
    }

    return;
}

#ifdef __cplusplus
}
#endif
