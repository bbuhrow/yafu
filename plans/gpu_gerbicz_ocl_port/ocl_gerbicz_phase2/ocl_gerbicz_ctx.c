#include <stdio.h>
#include <string.h>
#include "ocl_gerbicz_ctx.h"
#include "collision_bucket.h"   /* LOG2_NUM_BUCKETS, BUCKET_HASH_MIX -- single
                                   source of truth, see that header's own
                                   drift warning */

int
ocl_gerbicz_device_init(ocl_gerbicz_device_t *gd, gpu_info_t *info)
{
    memset(gd, 0, sizeof(*gd));
    if (ocl_device_init(&gd->dev, info) != 0)
        return -1;

    if (!ocl_check_min_version(&gd->dev, 2, 0)) {
        fprintf(stderr,
            "ocl_gerbicz_device_init: device '%s' reports OpenCL %d.%d, "
            "need >= 2.0 (locked decision #1)\n",
            info->name, gd->dev.version_major, gd->dev.version_minor);
        ocl_device_free(&gd->dev);
        return -1;
    }
    return 0;
}

void
ocl_gerbicz_device_free(ocl_gerbicz_device_t *gd)
{
    if (gd->program_collision) clReleaseProgram(gd->program_collision);
    if (gd->program_sieve)     clReleaseProgram(gd->program_sieve);
    gd->program_collision = gd->program_sieve = NULL;
    ocl_device_free(&gd->dev);
}

void
ocl_gerbicz_build_opts(const ocl_gerbicz_device_t *gd, char *buf, size_t buf_sz)
{
    snprintf(buf, buf_sz,
             "-DHAVE_SUBGROUPS=%d -DLOG2_NUM_BUCKETS=%uu -DBUCKET_HASH_MIX=0x%llxULL",
             gd->dev.has_subgroups ? 1 : 0,
             (unsigned)LOG2_NUM_BUCKETS,
             (unsigned long long)BUCKET_HASH_MIX);
}
