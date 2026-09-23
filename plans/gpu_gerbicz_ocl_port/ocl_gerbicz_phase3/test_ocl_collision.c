/* test_ocl_collision.c -- Phase 3 task 3.5 validation.
 *
 * Uses the REAL Phase 0 harness (gen.c/cpuref_exact.c/cpuref_stats.c/
 * cmp.c/collcase_io.c, linked directly, not shelled out to) to:
 *   1. generate a collcase_t (gen_generate, with_expected=1 -- this
 *      computes the CPU reference via cpuref_exact/cpuref_stats),
 *   2. upload it, run this phase's OpenCL collision engine,
 *   3. read back found_array + stats, write an "observed" collcase
 *      file (same input arrays, our own expected-output section),
 *   4. call cmp_run() -- the harness's own comparator, not a hand
 *      rolled one -- against a reference file written from the
 *      already-computed expected section.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ocl_collision.h"
#include "collcase.h"
#include "gen.h"
#include "cpuref.h"
#include "cmp.h"

/* Host-side mirrors of the device's found_t/specialq_t (ocl_collision_kernels.cl)
 * -- identical layout to collcase.h's own cc_found_t/cc_specialq_t (both are
 * plain structs of fixed-width scalars with natural alignment; collcase.h's own
 * comment already requires this to track stage1_core.h exactly), so reuse them
 * by alias rather than redeclaring. */
typedef cc_found_t found_t;
typedef cc_specialq_t specialq_t;

static int g_pass = 0, g_fail = 0;

static void
report(const char *name, int ok)
{
    if (ok) { g_pass++; printf("  [PASS] %s\n", name); }
    else    { g_fail++; printf("  [FAIL] %s\n", name); }
}

/* Runs the OpenCL engine on collcase `c`, fills `out` with the same
 * shape as `c` but with an expected-output section built from the
 * engine's own found_array + stats. Returns the ocl_collision_status_t. */
static ocl_collision_status_t
run_case(ocl_collision_engine_t *e, ocl_gerbicz_device_t *gd,
         const collcase_t *c, collcase_t *out)
{
    cl_int err;
    cl_mem d_keys, d_values, d_qbatch, d_found;
    ocl_collision_data_t data;
    ocl_collision_status_t st;
    found_t *h_found;
    uint32_t found_count, entry_count, i;

    memset(out, 0, sizeof(*out));
    out->n = c->n; out->root_bytes = c->root_bytes; out->key_bits = c->key_bits;
    out->shift = c->shift; out->bucket_hash = c->bucket_hash; out->num_q = c->num_q;
    out->hash_word_cap = c->hash_word_cap;
    out->keys = (uint64_t *)malloc((size_t)c->n * sizeof(uint64_t));
    out->values = (uint32_t *)malloc((size_t)c->n * sizeof(uint32_t));
    out->q_batch = (cc_specialq_t *)malloc((size_t)c->num_q * sizeof(cc_specialq_t));
    memcpy(out->keys, c->keys, (size_t)c->n * sizeof(uint64_t));
    memcpy(out->values, c->values, (size_t)c->n * sizeof(uint32_t));
    memcpy(out->q_batch, c->q_batch, (size_t)c->num_q * sizeof(cc_specialq_t));

    /* Upload. Keys on disk/in memory are already the zero-extended
     * 64-bit pattern (collcase.h's own doc comment); the engine reads
     * them at root_bytes width via a byte pointer + cast (see
     * ocl_scatter_roots), so upload at native width to match exactly
     * what a real per-batch buffer would contain. */
    if (c->root_bytes == 4) {
        uint32_t *tmp = (uint32_t *)malloc((size_t)c->n * sizeof(uint32_t));
        for (i = 0; i < c->n; i++) tmp[i] = (uint32_t)c->keys[i];
        d_keys = clCreateBuffer(gd->dev.context, CL_MEM_READ_ONLY,
                (size_t)c->n * sizeof(uint32_t), NULL, &err);
        clEnqueueWriteBuffer(e->th.queue, d_keys, CL_TRUE, 0,
                (size_t)c->n * sizeof(uint32_t), tmp, 0, NULL, NULL);
        free(tmp);
    } else {
        d_keys = clCreateBuffer(gd->dev.context, CL_MEM_READ_ONLY,
                (size_t)c->n * sizeof(uint64_t), NULL, &err);
        clEnqueueWriteBuffer(e->th.queue, d_keys, CL_TRUE, 0,
                (size_t)c->n * sizeof(uint64_t), c->keys, 0, NULL, NULL);
    }
    d_values = clCreateBuffer(gd->dev.context, CL_MEM_READ_ONLY,
            (size_t)c->n * sizeof(uint32_t), NULL, &err);
    clEnqueueWriteBuffer(e->th.queue, d_values, CL_TRUE, 0,
            (size_t)c->n * sizeof(uint32_t), c->values, 0, NULL, NULL);

    d_qbatch = clCreateBuffer(gd->dev.context, CL_MEM_READ_ONLY,
            (size_t)(c->num_q ? c->num_q : 1) * sizeof(specialq_t), NULL, &err);
    if (c->num_q) {
        /* cc_specialq_t and the kernel's specialq_t share layout
         * (collcase.h's own comment: field order/types must track
         * stage1_core.h exactly) -- direct byte copy is valid. */
        clEnqueueWriteBuffer(e->th.queue, d_qbatch, CL_TRUE, 0,
                (size_t)c->num_q * sizeof(specialq_t), c->q_batch, 0, NULL, NULL);
    }

    h_found = (found_t *)calloc(OCL_COLL_FOUND_ARRAY_SIZE, sizeof(found_t));
    d_found = clCreateBuffer(gd->dev.context, CL_MEM_READ_WRITE,
            OCL_COLL_FOUND_ARRAY_SIZE * sizeof(found_t), NULL, &err);
    clEnqueueWriteBuffer(e->th.queue, d_found, CL_TRUE, 0,
            OCL_COLL_FOUND_ARRAY_SIZE * sizeof(found_t), h_found, 0, NULL, NULL);

    memset(&data, 0, sizeof(data));
    data.keys_in = d_keys; data.keys_in_off = 0;
    data.data_in = d_values; data.data_in_off = 0;
    data.q_batch = d_qbatch; data.q_batch_off = 0;
    data.found_array = d_found; data.found_array_off = 0;
    data.num_elements = c->n;
    data.key_bits = c->key_bits;
    data.root_bytes = c->root_bytes;
    data.shift = c->shift;
    data.bucket_hash = (int)c->bucket_hash;
    data.debug = 0;
    data.collect_stats = 1;

    st = ocl_collision_run(e, &data);

    clEnqueueReadBuffer(e->th.queue, d_found, CL_TRUE, 0,
            OCL_COLL_FOUND_ARRAY_SIZE * sizeof(found_t), h_found, 0, NULL, NULL);

    found_count = h_found[0].p1;
    entry_count = found_count < (OCL_COLL_FOUND_ARRAY_SIZE - 1u) ?
            found_count : (OCL_COLL_FOUND_ARRAY_SIZE - 1u);

    out->has_expected = 1;
    out->found_count = found_count;
    out->entry_count = entry_count;
    out->entries = (cc_found_t *)malloc((size_t)(entry_count ? entry_count : 1) * sizeof(cc_found_t));
    for (i = 0; i < entry_count; i++) {
        found_t *f = &h_found[i + 1];
        out->entries[i].p1 = f->p1; out->entries[i].p2 = f->p2;
        out->entries[i].q = f->q; out->entries[i].pad = 0;
        out->entries[i].qroot = f->qroot; out->entries[i].offset = f->offset;
    }
    out->stats.candidate_count = data.candidate_count;
    out->stats.dedup_count = data.dedup_count;
    out->stats.value_match_count = data.value_match_count;
    out->stats.bucket_max = data.bucket_max;
    memcpy(out->stats.filter_iters_hist, data.filter_iters_hist,
           sizeof(out->stats.filter_iters_hist));

    free(h_found);
    clReleaseMemObject(d_keys); clReleaseMemObject(d_values);
    clReleaseMemObject(d_qbatch); clReleaseMemObject(d_found);
    return st;
}

static void
test_one(ocl_collision_engine_t *e, ocl_gerbicz_device_t *gd,
         const char *name, gen_params_t *gp)
{
    collcase_t c, observed;
    ocl_collision_status_t st;
    char obs_path[128], ref_path[128];
    int rc;

    collcase_init(&c);
    gp->with_expected = 1;
    if (gen_generate(gp, &c) != 0) {
        report(name, 0);
        return;
    }

    st = run_case(e, gd, &c, &observed);
    if (st != OCL_COLLISION_OK) {
        fprintf(stderr, "  (%s: ocl_collision_run returned status %d)\n", name, (int)st);
        report(name, 0);
        collcase_free(&c); collcase_free(&observed);
        return;
    }

    snprintf(obs_path, sizeof(obs_path), "/tmp/occ_%s_obs.collcase", name);
    snprintf(ref_path, sizeof(ref_path), "/tmp/occ_%s_ref.collcase", name);
    collcase_write(obs_path, &observed);
    collcase_write(ref_path, &c); /* c already has_expected=1 from gen_generate */

    rc = cmp_run(obs_path, ref_path);
    report(name, rc == 0);

    collcase_free(&c);
    collcase_free(&observed);
}

int
main(void)
{
    gpu_config_t config;
    ocl_gerbicz_device_t gd;
    ocl_primitives_t prim;
    ocl_collision_engine_t eng;
    gen_params_t gp;

    gpu_init(&config);
    if (config.num_gpu == 0) { fprintf(stderr, "no OpenCL devices found\n"); return 1; }
    if (ocl_gerbicz_device_init(&gd, &config.info[0]) != 0) return 1;
    printf("device: %s  cl_c_std=%s local_mem=%llu\n", config.info[0].name,
           gd.dev.cl_c_std, (unsigned long long)gd.dev.local_mem_size);

    if (ocl_primitives_init(&prim, &gd, "ocl_primitives_kernels.cl", ".") != 0) return 1;
    if (ocl_collision_init(&eng, &gd, &prim, "ocl_collision_kernels.cl", ".") != 0) return 1;

    printf("\n== small edge cases (harness's own 0.3 edge suite included) ==\n");
    gen_params_defaults(&gp);
    gp.n = 200; gp.seed = 1;
    test_one(&eng, &gd, "tiny_n200", &gp);

    gen_params_defaults(&gp);
    gp.n = 1; gp.seed = 2; gp.include_edge_cases = 0;
    test_one(&eng, &gd, "n1_noedge", &gp);

    printf("\n== one-workgroup boundary (COLL_BLOCK_THREADS=128 buckets) ==\n");
    gen_params_defaults(&gp);
    gp.n = 5000; gp.seed = 3; gp.collision_density = 0.05;
    test_one(&eng, &gd, "n5000_density", &gp);

    printf("\n== bucket-grow loop (forces ensure_capacity's retry path) ==\n");
    gen_params_defaults(&gp);
    gp.n = 2000; gp.seed = 4; gp.bucket_skew_count = 500;
    test_one(&eng, &gd, "skew500", &gp);

    printf("\n== multiplicative bucket hash ==\n");
    gen_params_defaults(&gp);
    gp.n = 3000; gp.seed = 5; gp.bucket_hash = 1;
    test_one(&eng, &gd, "bhash1", &gp);

    printf("\n== root_bytes=8 (r64) ==\n");
    gen_params_defaults(&gp);
    gp.n = 3000; gp.seed = 6; gp.root_bytes = 8; gp.key_bits = 40;
    test_one(&eng, &gd, "root8", &gp);

    printf("\n== several workgroups, real-ish scale ==\n");
    gen_params_defaults(&gp);
    gp.n = 200000; gp.seed = 7; gp.collision_density = 0.02;
    test_one(&eng, &gd, "n200000", &gp);

    printf("\n== hash_word_cap override (small, to force filter iteration) ==\n");
    gen_params_defaults(&gp);
    gp.n = 20000; gp.seed = 8; gp.hash_word_cap = 32; gp.collision_density = 0.1;
    test_one(&eng, &gd, "smallcap", &gp);

    printf("\n%d passed, %d failed\n", g_pass, g_fail);

    ocl_collision_free(&eng);
    ocl_primitives_free(&prim);
    ocl_gerbicz_device_free(&gd);
    return g_fail ? 1 : 0;
}
