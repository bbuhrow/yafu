/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#include "lanczos_gpu.h"
#include "lanczos_gpu_core.h"

static const char * gpu_kernel_names[] = 
{
	"lanczos_kernel_mask",
	"lanczos_kernel_xor",
	"lanczos_kernel_inner_prod",
	"lanczos_kernel_outer_prod",
};
 
typedef struct {
	uint32 row_off;
	uint32 col_off;
} entry_idx_t;

#if 0
// Tried using compressible memory on an A100. Did not help 
static CUresult setProp(CUmemAllocationProp *prop, int UseCompressibleMemory)
{
    CUdevice currentDevice;
    CUDA_TRY(cuCtxGetDevice(&currentDevice))

    memset(prop, 0, sizeof(CUmemAllocationProp));
    prop->type = CU_MEM_ALLOCATION_TYPE_PINNED;
    prop->location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    prop->location.id = currentDevice;

    if (UseCompressibleMemory)
        prop->allocFlags.compressionType = CU_MEM_ALLOCATION_COMP_GENERIC;

    return CUDA_SUCCESS;
}

CUresult allocateCompressible(void **adr, size_t size, int UseCompressibleMemory)
{
    CUmemAllocationProp prop = {};
    setProp(&prop, UseCompressibleMemory);

    size_t granularity = 0;
    CUDA_TRY(cuMemGetAllocationGranularity(&granularity, &prop, CU_MEM_ALLOC_GRANULARITY_MINIMUM))
    size = ((size - 1) / granularity + 1) * granularity;
    CUdeviceptr dptr;
    CUDA_TRY(cuMemAddressReserve(&dptr, size, 0, 0, 0))
    
	CUmemGenericAllocationHandle allocationHandle;
    CUDA_TRY(cuMemCreate(&allocationHandle, size, &prop, 0))

    // Check if cuMemCreate was able to allocate compressible memory.
    if (UseCompressibleMemory) {
        CUmemAllocationProp allocationProp = {};
        cuMemGetAllocationPropertiesFromHandle(&allocationProp, allocationHandle);
        if (allocationProp.allocFlags.compressionType != CU_MEM_ALLOCATION_COMP_GENERIC) {
            printf("Could not allocate compressible memory...\n");
            exit(-1);
        }
    }

    CUDA_TRY(cuMemMap(dptr, size, 0, allocationHandle, 0))
    CUDA_TRY(cuMemRelease(allocationHandle))

    CUmemAccessDesc accessDescriptor;
    accessDescriptor.location.id = prop.location.id;
    accessDescriptor.location.type = prop.location.type;
    accessDescriptor.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;

    CUDA_TRY(cuMemSetAccess(dptr, size, &accessDescriptor, 1))

    *adr = (void *)dptr;
    return CUDA_SUCCESS;
}

CUresult freeCompressible(void *ptr, size_t size, int UseCompressibleMemory)
{
    CUmemAllocationProp prop = {};
    setProp(&prop, UseCompressibleMemory);

    size_t granularity = 0;
    CUDA_TRY(cuMemGetAllocationGranularity(&granularity, &prop, CU_MEM_ALLOC_GRANULARITY_MINIMUM))
    size = ((size - 1) / granularity + 1) * granularity;

    if (ptr == NULL) return CUDA_SUCCESS;
    if (cuMemUnmap((CUdeviceptr)ptr, size) != CUDA_SUCCESS ||
        cuMemAddressFree((CUdeviceptr)ptr, size) != CUDA_SUCCESS)
        return CUDA_ERROR_INVALID_VALUE;
    return CUDA_SUCCESS;
}
#endif

/*-------------------------------------------------------------------*/
static void copy_dense(packed_matrix_t *p) 
{
	/* copy the dense arrays to device memory */

	uint32 i, j, k;
	uint32 ncols = p->ncols;
	gpudata_t *d = (gpudata_t *)p->extra;
	uint32 num_dense_blocks = (p->num_dense_rows + VBITS - 1) / VBITS;
	v_t *tmp = (v_t *)xmalloc(ncols * sizeof(v_t));

	d->dense_blocks = (CUdeviceptr *)xmalloc(num_dense_blocks *
						sizeof(CUdeviceptr));

	for (i = 0; i < num_dense_blocks; i++) {

		for (j = 0; j < ncols; j++) {
			la_col_t *col = p->unpacked_cols + j;
			uint32 *src = col->data + col->weight;
			for (k = 0; k < VWORDS; k++) {
				uint32 t = i * VWORDS + k;
				tmp[j].w[k] = (uint64)src[2 * t + 1] << 32 |
					(uint64)src[2 * t];
			}
		}

		if (d->use_cudamanaged) {
			CUDA_TRY(cuMemAllocManaged(&d->dense_blocks[i],
				ncols * sizeof(v_t),
				CU_MEM_ATTACH_GLOBAL))
			CUDA_TRY(cuMemcpy(d->dense_blocks[i],
				(CUdeviceptr) tmp,
				ncols * sizeof(v_t)))
			CUDA_TRY(cuMemAdvise(d->dense_blocks[i],
				ncols * sizeof(v_t),
				CU_MEM_ADVISE_SET_READ_MOSTLY,
				d->gpu_info->device_handle))
		} else {
			/* CUDA_TRY(allocateCompressible((void **)&d->dense_blocks[i], ncols * sizeof(v_t), 1)) */
			CUDA_TRY(cuMemAlloc(&d->dense_blocks[i],
					ncols * sizeof(v_t)))
			CUDA_TRY(cuMemcpyHtoD(d->dense_blocks[i], tmp,
					ncols * sizeof(v_t)))
		}
	}

	free(tmp);
}

/*-------------------------------------------------------------------*/
static uint32 extract_block(la_col_t *cols,
			uint32 row_min, uint32 row_max,
			uint32 col_min, uint32 col_max,
			uint32 nnz, uint32 *blocksize,
			entry_idx_t **entries_in, 
			uint32 *max_entries_in)
{
	uint32 i, j;
	uint32 num_entries = 0;
	entry_idx_t *entries = *entries_in;
	uint32 max_entries = *max_entries_in;

	for (i = col_min; (i < col_max) && (num_entries < nnz); i++) {

		la_col_t *col = cols + i;

		for (j = 0; j < col->weight; j++) {
			uint32 idx = col->data[j];

			if (idx >= row_max)
				break;

			if (idx >= row_min) {

				entry_idx_t *e;

				if (num_entries == max_entries) {
					max_entries *= 2;
					entries = (entry_idx_t *)xrealloc(
							entries, 
							max_entries *
							sizeof(entry_idx_t));
				}

				e = entries + num_entries++;
				e->row_off = idx;
				e->col_off = i;
			}
		}
	}

	*blocksize = i - col_min;
	*entries_in = entries;
	*max_entries_in = max_entries;
	return num_entries;
}

/*-------------------------------------------------------------------*/
static uint32 extract_block_trans(la_col_t *cols,
			uint32 row_min, uint32 row_max,
			uint32 col_min, uint32 col_max,
			uint32 nnz, uint32 num_dense_rows,
			uint32 *blocksize,
			entry_idx_t **entries_in,
			uint32 *max_entries_in)
{
	uint32 i, j;
	uint32 num_entries = 0;
	uint32 my_row_max, my_blocksize;
	entry_idx_t *entries = *entries_in;
	uint32 max_entries = *max_entries_in;
	uint32 min_nnz = 9 * (nnz / 10);
	uint32 max_nnz = 11 * (nnz / 10);

	/* Need to figure out what my_row_max to use to get about nnz nonzeros */

	my_blocksize = *blocksize;
	my_row_max = row_min + my_blocksize;
	if (my_row_max > row_max) {
		my_row_max = row_max;
		my_blocksize = row_max - row_min;
	}

	while (1) {
		num_entries = 0;
		for (i = col_min; i < col_max; i++) {
			la_col_t *col = cols + i;
			for (j = 0; j < col->weight; j++) {
				uint32 idx = col->data[j];
				if (idx >= my_row_max) break;
				if (idx >= row_min) num_entries++;
			}
		}
		if (num_entries > max_nnz) {
			my_blocksize = 4 * (my_blocksize / 5);
			if ((row_min == 0) && (my_blocksize <= num_dense_rows)) {
				/* just grab a few and continue */
				my_row_max = num_dense_rows + 10;
				break;
			}
			my_row_max = row_min + my_blocksize;
			if (my_blocksize == 2) break;
			min_nnz = 0;
			continue;
		}
		if (num_entries < min_nnz) {
			my_blocksize = 5 * (my_blocksize / 4);
			my_row_max = row_min + my_blocksize;
			if (my_row_max >= row_max) {
				my_row_max = row_max;
				break;
			}
			max_nnz = (uint32)(-1);
			continue;
		}
		break;
	}

	if (num_entries == 0) /* shouldn't happen */
		my_row_max = MIN(my_row_max + 10, row_max);
	num_entries = 0;
	for (i = col_min; i < col_max; i++) {

		la_col_t *col = cols + i;

		for (j = 0; j < col->weight; j++) {
			uint32 idx = col->data[j];

			if (idx >= my_row_max)
				break;

			if (idx >= row_min) {

				entry_idx_t *e;

				if (num_entries == max_entries) {
					max_entries *= 2;
					entries = (entry_idx_t *)xrealloc(
							entries, 
							max_entries *
							sizeof(entry_idx_t));
				}

				e = entries + num_entries++;
				e->row_off = idx;
				e->col_off = i;
			}
		}
	}

	*blocksize = my_row_max - row_min;
	*entries_in = entries;
	*max_entries_in = max_entries;
	return num_entries;
}

/*-------------------------------------------------------------------*/
static int compare_row_off(const void *x, const void *y) {
	entry_idx_t *xx = (entry_idx_t *)x;
	entry_idx_t *yy = (entry_idx_t *)y;

	if (xx->row_off > yy->row_off)
		return 1;
	if (xx->row_off < yy->row_off)
		return -1;

	return (int)xx->col_off - (int)yy->col_off;
}

/*-------------------------------------------------------------------*/
static void radix_sort(entry_idx_t *arr, uint32 n) {

	/* simple radix sort, much faster than qsort() */

	uint32 i, pass, skip;
	uint64 *a, *b, *from, *to, *temp;

	a = (uint64 *) malloc(n * sizeof(uint64));
	b = (uint64 *) malloc(n * sizeof(uint64));

	for (i = 0; i < n; i++) {
		entry_idx_t *e = arr + i;
		a[i] = ((uint64)(e->row_off) << 32) | (uint64)(e->col_off);
	}

	from = a;
	to = b;
	skip = 0;
	for (pass = 0; pass < 8; pass++)  {
		uint32 box[256] = { 0 };

		for (i = 0; i < n; i++) box[ (from[i] >> (8*pass)) & 255]++;
		if (box[0] == n) { /* this word is all 0's, don't need to sort */
			skip++;
			continue;
		}
		for (i = 1; i < 256; i++) box[i] += box[i-1];
		for (i = n - 1; i != (uint32)(-1); i--) to[--box[(from[i] >> (8*pass)) & 255]] = from[i];

		temp = from;
		from = to;
		to = temp;
	}

	if (skip & 1) to = b;
	else to = a;
	for (i = 0; i < n; i++) {
		entry_idx_t *e = arr + i;
		e->row_off = (uint32)(to[i] >> 32);
		e->col_off = (uint32)(to[i]);
	}

	free(a);
	free(b);
}

/*-------------------------------------------------------------------*/
static void pack_matrix_block(gpudata_t *d, block_row_t *b,
			entry_idx_t *entries, uint32 num_entries,
			uint32 row_min, uint32 row_max, 
			uint32 col_min, uint32 col_max,
			uint32 is_trans)
{

	uint32 i, j;
	uint32 num_rows = row_max - row_min;

	/* convert a block of matrix rows from COO to CSR format */

	uint32 *col_entries = (uint32 *)xmalloc(num_entries * 
					sizeof(uint32));
	uint32 *row_entries = (uint32 *)xcalloc(num_rows + 1,
					sizeof(uint32));

	if (is_trans) {
		for (i = 0; i < num_entries; i++) {
			entry_idx_t *e = entries + i;
			j = e->row_off;
			e->row_off = e->col_off;
			e->col_off = j;
		}
	}
	else {
		/* qsort(entries, num_entries, sizeof(entry_idx_t),
				compare_row_off); */
		radix_sort(entries, num_entries);
	}

	for (i = j = 0; i < num_entries; i++, j++) {

		entry_idx_t *e = entries + i;

		col_entries[i] = e[0].col_off - col_min;

		if (i > 0 && e[0].row_off != e[-1].row_off) {
			row_entries[e[-1].row_off - row_min] = j;
			j = 0;
		}
	}
	row_entries[entries[i-1].row_off - row_min] = j;

	for (i = j = 0; i < num_rows; i++) {
		uint32 t = row_entries[i];
		row_entries[i] = j;
		j += t;
	}
	row_entries[num_rows] = num_entries;

	b->num_rows = num_rows;
	b->num_cols = col_max - col_min;
	b->num_col_entries = num_entries;
	printf("%u %u %u\n", num_entries, num_rows, b->blocksize);

	if (d->use_cudamanaged) {
		CUDA_TRY(cuMemAllocManaged(&b->col_entries,
				num_entries * sizeof(uint32),
				CU_MEM_ATTACH_GLOBAL))
		CUDA_TRY(cuMemcpy(b->col_entries,
				(CUdeviceptr) col_entries,
				num_entries * sizeof(uint32)))
		CUDA_TRY(cuMemAdvise(b->col_entries,
				num_entries * sizeof(uint32),
				CU_MEM_ADVISE_SET_READ_MOSTLY,
				d->gpu_info->device_handle))

		CUDA_TRY(cuMemAllocManaged(&b->row_entries,
				(num_rows + 1) * sizeof(uint32),
				CU_MEM_ATTACH_GLOBAL))
		CUDA_TRY(cuMemcpy(b->row_entries,
				(CUdeviceptr) row_entries,
				(num_rows + 1) * sizeof(uint32)))
		CUDA_TRY(cuMemAdvise(b->row_entries,
				(num_rows + 1) * sizeof(uint32),
				CU_MEM_ADVISE_SET_READ_MOSTLY,
				d->gpu_info->device_handle))
	} else {
		/* CUDA_TRY(allocateCompressible((void **)&b->col_entries, num_entries * sizeof(uint32), 1)) */
		CUDA_TRY(cuMemAlloc(&b->col_entries,
				num_entries * sizeof(uint32)))
		CUDA_TRY(cuMemcpyHtoD(b->col_entries,
				col_entries,
				num_entries * sizeof(uint32)))

		/* CUDA_TRY(allocateCompressible((void **)&b->row_entries, (num_rows + 1) * sizeof(uint32), 1)) */
		CUDA_TRY(cuMemAlloc(&b->row_entries,
				(num_rows + 1) * sizeof(uint32)))
		CUDA_TRY(cuMemcpyHtoD(b->row_entries,
				row_entries,
				(num_rows + 1) * sizeof(uint32)))
	}

	free(col_entries);
	free(row_entries);
}

/*-------------------------------------------------------------------*/
static void gpu_matrix_init(packed_matrix_t *p) {

	uint32 start_row = 0;
	uint32 start_col = 0;
	uint32 blocksize;
	gpudata_t *d = (gpudata_t *)p->extra;

	uint32 num_block_rows = 0;
	uint32 num_trans_block_rows = 0;
	uint32 num_block_rows_alloc = 100;
	uint32 num_trans_block_rows_alloc = 100;
	block_row_t *block_rows = (block_row_t *)xmalloc(
					num_block_rows_alloc *
					sizeof(block_row_t));
	block_row_t *trans_block_rows = (block_row_t *)xmalloc(
					num_trans_block_rows_alloc *
					sizeof(block_row_t));

	uint32 num_entries_alloc = 10000;
	entry_idx_t *entries = (entry_idx_t *)xmalloc(
					num_entries_alloc *
					sizeof(entry_idx_t));

	/* deal with the dense rows */

	copy_dense(p);

	/* deal with the sparse rows */

	printf("converting matrix to CSR and copying it onto the GPU\n");

	while (start_col < p->ncols) {

		block_row_t *b;
		uint32 num_entries;

		num_entries = extract_block(p->unpacked_cols,
					0, p->nrows,
					start_col, p->ncols,
				 	p->block_nnz,
					&blocksize,
					&entries,
					&num_entries_alloc);

		if (num_entries > 2147483647) {
			printf("max column entries is 2147483647\n");
			printf("adjust preferred block to compensate\n");
			exit(42);
		}

		if (num_block_rows == num_block_rows_alloc) {
			num_block_rows_alloc *= 2;
			block_rows = (block_row_t *)xrealloc(
					block_rows,
					num_block_rows_alloc *
					sizeof(block_row_t));
		}

		b = block_rows + num_block_rows++;
		b->blocksize = blocksize;
		pack_matrix_block(d, b, entries, num_entries,
				0, p->nrows, 
				start_col,
				start_col + b->blocksize,
				0);

		start_col += b->blocksize;
	}

	d->num_block_rows = num_block_rows;
	d->block_rows = block_rows;

	/* handle the transpose of the matrix */

	/* First rows are heavy so suggest a small initial blocksize */
	blocksize = p->block_nnz / 10000;
	while (start_row < p->nrows) {

		block_row_t *b;
		uint32 num_entries;

		num_entries = extract_block_trans(p->unpacked_cols,
					start_row,
					p->nrows,
					0, p->ncols,
					p->block_nnz,
					p->num_dense_rows,
					&blocksize,
					&entries,
					&num_entries_alloc);

		if (num_entries > 2147483647) {
			printf("max column entries is 2147483647\n");
			printf("adjust preferred transpose block to compensate\n");
			exit(42);
		}

		if (num_trans_block_rows == num_trans_block_rows_alloc) {
			num_trans_block_rows_alloc *= 2;
			trans_block_rows = (block_row_t *)xrealloc(
					trans_block_rows,
					num_trans_block_rows_alloc *
					sizeof(block_row_t));
		}

		b = trans_block_rows + num_trans_block_rows++;
		b->blocksize = blocksize;
		pack_matrix_block(d, b, entries, num_entries,
				0, p->ncols, 
				start_row,
				start_row + b->blocksize,
				1);

		start_row += b->blocksize;
	}

	d->num_trans_block_rows = num_trans_block_rows;
	d->trans_block_rows = trans_block_rows;

	free(entries);
}

/*-------------------------------------------------------------------*/
static void gpu_matrix_free(packed_matrix_t *p) {

	uint32 i;
	gpudata_t *d = (gpudata_t *)p->extra;

	for (i = 0; i < d->num_block_rows; i++) {
		block_row_t *b = d->block_rows + i;

		CUDA_TRY(cuMemFree(b->row_entries))
		CUDA_TRY(cuMemFree(b->col_entries))
	}
	free(d->block_rows);

	for (i = 0; i < d->num_trans_block_rows; i++) {
		block_row_t *b = d->trans_block_rows + i;

		CUDA_TRY(cuMemFree(b->row_entries))
		CUDA_TRY(cuMemFree(b->col_entries))
	}
	free(d->trans_block_rows);

	for (i = 0; i < (p->num_dense_rows + VBITS - 1) / VBITS; i++)
		CUDA_TRY(cuMemFree(d->dense_blocks[i]))
	free(d->dense_blocks);
}

/*------------------------------------------------------------------------*/
static void
load_spmv_engine(msieve_obj *obj, gpudata_t *d)
{
	char libname[256];
	#if defined(WIN32) || defined(_WIN64)
	const char *suffix = ".dll";
	#else
	const char *suffix = ".so";
	#endif

	if (d->gpu_info->compute_version_major < 2) {
		printf("error: GPU compute capability >= 2.0 required\n");
		exit(-1);
	}

	sprintf(libname, "cub/spmv_engine%s", suffix);

	/* override from input args */

	if (obj->nfs_args != NULL) {
		char *tmp = strstr(obj->nfs_args, "spmvlib=");

		if (tmp != NULL) {
			uint32 i;
			for (i = 0, tmp += 8; i < sizeof(libname) - 1; i++) {
				if (*tmp == 0 || isspace(*tmp))
					break;

				libname[i] = *tmp++;
			}
			libname[i] = 0;
		}
	}

	d->spmv_engine_handle = load_dynamic_lib(libname);
	if (d->spmv_engine_handle == NULL) {
		printf("error: failed to load GPU matrix multiply engine\n");
		exit(-1);
	}

	/* the spmv engine uses the same CUDA context */

	d->spmv_engine_init = get_lib_symbol(
					d->spmv_engine_handle,
					"spmv_engine_init");
	d->spmv_engine_free = get_lib_symbol(
					d->spmv_engine_handle,
					"spmv_engine_free");
	d->spmv_engine_run = get_lib_symbol(
					d->spmv_engine_handle,
					"spmv_engine_run");
	if (d->spmv_engine_init == NULL ||
	    d->spmv_engine_free == NULL ||
	    d->spmv_engine_run == NULL) {
		printf("error: cannot find GPU matrix multiply function\n");
		exit(-1);
	}
}

/*-------------------------------------------------------------------*/
void matrix_extra_init(msieve_obj *obj, packed_matrix_t *p,
			uint32 first_block_size) {

	uint32 i;
	int check_vbits;
	gpudata_t *d;
	gpu_config_la_t gpu_config;
	gpu_info_la_t *gpu_info;

	/* select card, save info struct */

	gpu_init_la(&gpu_config);
	if (gpu_config.num_gpu == 0) {
		printf("error: no CUDA-enabled GPUs found\n");
		exit(-1);
	}
	if (obj->which_gpu >= (uint32)gpu_config.num_gpu) {
		printf("error: GPU %u does not exist "
			"or is not CUDA-enabled\n", obj->which_gpu);
		exit(-1);
	}

	p->extra = d = (gpudata_t *)xcalloc(1, sizeof(gpudata_t));

	d->gpu_info = gpu_info = (gpu_info_la_t *)xmalloc(sizeof(gpu_info_la_t));
	memcpy(gpu_info, gpu_config.info + obj->which_gpu,
			sizeof(gpu_info_la_t)); 

	logprintf(obj, "using GPU %u (%s)\n", obj->which_gpu, gpu_info->name);
	logprintf(obj, "selected card has CUDA arch %d.%d\n",
			gpu_info->compute_version_major,
			gpu_info->compute_version_minor);

 	/* CUDA_TRY(cuDevicePrimaryCtxSetFlags(d->gpu_info->device_handle,
			CU_CTX_SCHED_BLOCKING_SYNC)) */

	/* initialize context */

	CUDA_TRY(cuCtxCreate(&d->gpu_context,
			CU_CTX_SCHED_BLOCKING_SYNC,
			d->gpu_info->device_handle))

	/* CUDA_TRY(cuDevicePrimaryCtxRetain(&d->gpu_context,
			d->gpu_info->device_handle)) */

	load_spmv_engine(obj, d);
	d->spmv_engine = d->spmv_engine_init(&check_vbits);
	if (check_vbits != VBITS) {
                printf("error: SpMV library compiled for VBITS=%d\n", check_vbits);
                exit(-1);
	}

	/* load kernels */

	CUDA_TRY(cuModuleLoad(&d->gpu_module, "lanczos_kernel.ptx"))

	d->launch = (gpu_launch_la_t *)xmalloc(NUM_GPU_FUNCTIONS *
				sizeof(gpu_launch_la_t));

	for (i = 0; i < NUM_GPU_FUNCTIONS; i++) {
		gpu_launch_la_t *launch = d->launch + i;

		gpu_launch_init_la(d->gpu_module, gpu_kernel_names[i],
				launch);

		launch->threads_per_block = MIN(256, 
				launch->threads_per_block);
	}

	/* allocate scratch arrays */

	CUDA_TRY(cuMemAlloc(&d->gpu_scratch, VBITS * sizeof(v_t)))

	/* Set preferred nonzeros per matrix block */

	p->block_nnz = 1750000000;
	if (obj->nfs_args != NULL) {
		const char *tmp;
		tmp = strstr(obj->nfs_args, "block_nnz=");
		if (tmp != NULL) p->block_nnz = (uint32)atoi(tmp + 10);
		if (p->block_nnz < 100000) p->block_nnz = 100000;
		if (p->block_nnz > 1750000000) p->block_nnz = 1750000000; 
	}
	printf("Nonzeros per block: %u\n", p->block_nnz);

	/* should we used CUDA managed memory to store the matrix */

	d->use_cudamanaged = 0;
	if (obj->nfs_args != NULL) {
		const char *tmp;
		tmp = strstr(obj->nfs_args, "use_managed=1");
		if (tmp != NULL) {
			if (d->gpu_info->concurrent_managed_access)
				d->use_cudamanaged = 2; /* can prefetch */
			else d->use_cudamanaged = 1;
			printf("Storing matrix in managed memory\n");
		}
	}
	
	/* Adjust L2 fetch granularity. Default is 128. Tried 32 for VBITS=256, but makes no difference */
	/* if (gpu_info->compute_version_major >= 8) CUDA_TRY(cuCtxSetLimit(CU_LIMIT_MAX_L2_FETCH_GRANULARITY, 32)) */

	/* set up the matrix on the card */

	gpu_matrix_init(p);
}

/*-------------------------------------------------------------------*/
void matrix_extra_free(packed_matrix_t *p) {

	gpudata_t *d = (gpudata_t *)p->extra;

	gpu_matrix_free(p);

	CUDA_TRY(cuMemFree(d->gpu_scratch))

	free(d->launch);

	d->spmv_engine_free(d->spmv_engine);
	unload_dynamic_lib(d->spmv_engine_handle);

	CUDA_TRY(cuCtxDestroy(d->gpu_context))
	/* CUDA_TRY(cuDevicePrimaryCtxRelease(d->gpu_info->device_handle)) */

	free(d->gpu_info);
	free(d);
}

/*-------------------------------------------------------------------*/
static void mul_packed_gpu(packed_matrix_t *p, 
				gpuvec_t *x, gpuvec_t *b) {

	uint32 i;
	uint32 start_col = 0;
	gpudata_t *d = (gpudata_t *)p->extra;

	CUDA_TRY(cuMemsetD8(b->gpu_vec, 0, 
			p->nrows * sizeof(v_t)));

	/* sweep through the matrix a block col at a time */

	for (i = 0; i < d->num_block_rows; i++) {

		block_row_t *blk = d->block_rows + i;
		spmv_data_t spmv_data;

		if (d->use_cudamanaged == 2) {
			CUDA_TRY(cuMemPrefetchAsync(blk->col_entries,
				blk->num_col_entries * sizeof(uint32),
				d->gpu_info->device_handle, 0))
			CUDA_TRY(cuMemPrefetchAsync(blk->row_entries,
				(blk->num_rows + 1) * sizeof(uint32),
				d->gpu_info->device_handle, 0))
		}
		spmv_data.num_rows = blk->num_rows;
		spmv_data.num_cols = blk->num_cols;
		spmv_data.num_col_entries = blk->num_col_entries;
		spmv_data.col_entries = blk->col_entries;
		spmv_data.row_entries = blk->row_entries;
		spmv_data.vector_in = (CUdeviceptr)((v_t *)x->gpu_vec + start_col);
		spmv_data.vector_out = b->gpu_vec;

		d->spmv_engine_run(d->spmv_engine, &spmv_data);

		start_col += blk->blocksize;
	}

	/* handle dense rows */

	for (i = 0; i < (p->num_dense_rows + VBITS - 1) / VBITS; i++) {
		if (d->use_cudamanaged == 2) {
			CUDA_TRY(cuMemPrefetchAsync(d->dense_blocks[i],
				p->ncols * sizeof(v_t),
				d->gpu_info->device_handle, 0))
		}
		mul_BxN_NxB_gpu(p, 
			d->dense_blocks[i], 
			x->gpu_vec, 
			(CUdeviceptr)((v_t *)b->gpu_vec + VBITS * i), 
			p->ncols);
	}
}

/*-------------------------------------------------------------------*/
static void mul_packed_trans_gpu(packed_matrix_t *p, 
				gpuvec_t *x, gpuvec_t *b) {

	uint32 i;
	uint32 start_row = 0;
	gpudata_t *d = (gpudata_t *)p->extra;

	CUDA_TRY(cuMemsetD8(b->gpu_vec, 0, 
			p->ncols * sizeof(v_t)));

	/* sweep through the matrix a block row at a time */

	for (i = 0; i < d->num_trans_block_rows; i++) {

		block_row_t *blk = d->trans_block_rows + i;
		spmv_data_t spmv_data;

		if (d->use_cudamanaged == 2) {
			CUDA_TRY(cuMemPrefetchAsync(blk->col_entries,
				blk->num_col_entries * sizeof(uint32),
				d->gpu_info->device_handle, 0))
			CUDA_TRY(cuMemPrefetchAsync(blk->row_entries,
				(blk->num_rows + 1) * sizeof(uint32),
				d->gpu_info->device_handle, 0))
		}
		spmv_data.num_rows = blk->num_rows;
		spmv_data.num_cols = blk->num_cols;
		spmv_data.num_col_entries = blk->num_col_entries;
		spmv_data.col_entries = blk->col_entries;
		spmv_data.row_entries = blk->row_entries;
		spmv_data.vector_in = (CUdeviceptr)((v_t *)x->gpu_vec + start_row);
		spmv_data.vector_out = b->gpu_vec;

		d->spmv_engine_run(d->spmv_engine, &spmv_data);

		start_row += blk->blocksize;
	}

	/* handle dense rows */

	for (i = 0; i < (p->num_dense_rows + VBITS - 1) / VBITS; i++) {
		if (d->use_cudamanaged == 2) {
			CUDA_TRY(cuMemPrefetchAsync(d->dense_blocks[i],
				p->ncols * sizeof(v_t),
				d->gpu_info->device_handle, 0))
		}
		mul_NxB_BxB_acc_gpu(p,
			d->dense_blocks[i],
			(CUdeviceptr)((v_t *)x->gpu_vec + VBITS * i),
			(CUdeviceptr)((v_t *)b->gpu_vec + VBITS * i),
			p->ncols);
	}
}

/*-------------------------------------------------------------------*/
void mul_core(packed_matrix_t *A, void *x_in, void *b_in) {
    
	gpuvec_t *x = (gpuvec_t *)x_in;
	gpuvec_t *b = (gpuvec_t *)b_in;

	mul_packed_gpu(A, x, b);

#ifdef LANCZOS_GPU_DEBUG
	{
		uint32 i, j;
		v_t *tmp = (v_t *) xmalloc(A->ncols * 
						sizeof(v_t));

		CUDA_TRY(cuMemcpyDtoH(tmp, b->gpu_vec, 
					A->nrows * sizeof(v_t)))
		CUDA_TRY(cuMemcpyDtoH(x->host_vec, x->gpu_vec, 
					A->ncols * sizeof(v_t)))

		mul_unpacked(A, x->host_vec, b->host_vec);

		for (i = 0; i < MIN(A->ncols, A->nrows); i++) {
			for (j = 0; j < VWORDS; j++) {				
				if (tmp[i].w[j] != b->host_vec[i].w[j]) { 
					printf("m error %u %" PRIx64 " %" PRIx64 "\n", 
							i, b->host_vec[i].w[j], tmp[i].w[j]);
					exit(-1);
				}
			}
		}

		free(tmp);
	}
#endif
}

/*-------------------------------------------------------------------*/
void mul_trans_core(packed_matrix_t *A, void *x_in, void *b_in) {
    
	gpuvec_t *x = (gpuvec_t *)x_in;
	gpuvec_t *b = (gpuvec_t *)b_in;

	mul_packed_trans_gpu(A, x, b);

#ifdef LANCZOS_GPU_DEBUG
	{
		uint32 i, j;
		v_t *tmp = (v_t *)xmalloc(A->ncols * 
						sizeof(v_t));

		CUDA_TRY(cuMemcpyDtoH(tmp, b->gpu_vec, 
					A->ncols * sizeof(v_t)))
		CUDA_TRY(cuMemcpyDtoH(x->host_vec, x->gpu_vec, 
					A->nrows * sizeof(v_t)))

		mul_trans_unpacked(A, x->host_vec, b->host_vec);

		for (i = 0; i < A->ncols; i++) {
			for (j = 0; j < VWORDS; j++) {				
				if (tmp[i].w[j] != b->host_vec[i].w[j]) { 
					printf("tr error %u %" PRIx64 " %" PRIx64 "\n", 
							i, b->host_vec[i].w[j], tmp[i].w[j]);
					exit(-1);
				}
			}
		}

		free(tmp);
	}
#endif
}

/*-------------------------------------------------------------------*/
size_t packed_matrix_sizeof(packed_matrix_t *p) {

	uint32 i;
	size_t mem_use, tot_mem_use;
	gpudata_t *d = (gpudata_t*) p->extra;

	/* account for the vectors used in the lanczos iteration */

#ifdef HAVE_MPI
	mem_use = (6 * p->nsubcols + 2 * 
			MAX(p->nrows, p->ncols)) * sizeof(v_t);
#else
	mem_use = 7 * p->max_ncols * sizeof(v_t);
#endif

	/* and for the vv kernel scratch array */

	mem_use += VBITS * sizeof(v_t);

	tot_mem_use = mem_use;
	printf("vector memory use: %.1f MB\n", (double)mem_use/1048576);

	/* and for the matrix */

	/* dense rows */
	mem_use = ((p->num_dense_rows + VBITS - 1) / VBITS) * p->ncols * sizeof(v_t);

	tot_mem_use += mem_use;
	printf("dense rows memory use: %.1f MB\n", (double)mem_use/1048576);

	mem_use = 0;

	/* matrix in CSR format */
	for (i = 0; i < d->num_block_rows; i++) {
		block_row_t *b = d->block_rows + i;
		mem_use += (b->num_rows + 1 + b->num_col_entries) * sizeof(uint32);
	}

	/* transpose matrix in CSR format */
	for (i = 0; i < d->num_trans_block_rows; i++) {
		block_row_t *b = d->trans_block_rows + i;
		mem_use += (b->num_rows + 1 + b->num_col_entries) * sizeof(uint32);
	}

	tot_mem_use += mem_use;
	printf("sparse matrix memory use: %.1f MB\n", (double)mem_use/1048576);

	return tot_mem_use;
}
