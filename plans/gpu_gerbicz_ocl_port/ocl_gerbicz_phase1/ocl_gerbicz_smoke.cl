/* ocl_gerbicz_smoke.cl -- Phase 1 pipeline smoke test.
 *
 * Deliberately trivial (vector add), but exercises the two conventions
 * that every later kernel must follow:
 *   - base-offset arguments (see ocl_gerbicz_ctx.h)
 *   - the HAVE_SUBGROUPS build-time gate (see ocl_gerbicz_ctx.h)
 * so a failure here means the *pipeline*, not kernel logic, is broken.
 */

#ifndef HAVE_SUBGROUPS
#define HAVE_SUBGROUPS 0
#endif

__kernel void smoke_vecadd(__global const uint *a, ulong a_off,
                            __global const uint *b, ulong b_off,
                            __global uint *out,      ulong out_off,
                            uint n)
{
    a   += a_off;
    b   += b_off;
    out += out_off;

    uint gid = (uint)get_global_id(0);
    if (gid >= n)
        return;

#if HAVE_SUBGROUPS
    /* Not a meaningful use of sub-groups -- just proves the gated path
     * compiles and runs when HAVE_SUBGROUPS=1 is passed as a build
     * option. Real sub-group kernels arrive in Phase 3/7. */
    uint sg_sum = sub_group_reduce_add(1u);
    out[gid] = a[gid] + b[gid] + (sg_sum - sg_sum); /* sg_sum cancels out */
#else
    out[gid] = a[gid] + b[gid] + 1000u;
#endif
}
