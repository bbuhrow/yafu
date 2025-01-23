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

#include "lanczos_gpu_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------*/
__global__ void
lanczos_kernel_mask(v_t *x, v_t mask, uint32 n)
{
	uint32 i;
	uint32 num_threads = gridDim.x * blockDim.x;
	uint32 grid_id = blockIdx.x * blockDim.x + threadIdx.x;

	for (i = grid_id; i < n; i += num_threads)
		x[i] = v_and(x[i], mask);
}

/*------------------------------------------------------------------------*/
__global__ void
lanczos_kernel_xor(v_t *dest, v_t *src, uint32 n)
{
	uint32 i;
	uint32 num_threads = gridDim.x * blockDim.x;
	uint32 grid_id = blockIdx.x * blockDim.x + threadIdx.x;

	for (i = grid_id; i < n; i += num_threads)
		dest[i] = v_xor(dest[i], src[i]);
}

/*------------------------------------------------------------------------*/
__global__ void
lanczos_kernel_inner_prod(v_t *y, v_t *v,
			v_t *x, uint32 n)
{
	uint32 i, j;
	uint32 num_threads = gridDim.x * blockDim.x;
	uint32 grid_id = blockIdx.x * blockDim.x + threadIdx.x;
	v_t acc;
	__shared__ v_t c[32*VWORDS][3];

	for (i = threadIdx.x; i < 32 * VWORDS; i += blockDim.x) {
		acc = x[2 * i];
		c[i][0] = acc;

		acc = v_xor(acc, x[2 * i + 1]);
		c[i][2] = acc;

		acc = v_xor(acc, x[2 * i]);
		c[i][1] = acc;
	}

	__syncthreads();

	for (i = grid_id; i < n; i += num_threads) {
		v_t vi = v[i];
		for (j = 0; j < VWORDS; j++) acc.w[j] = 0;

		for (j = 0; j < 32 * VWORDS; j++) {
			uint32 k = (vi.w[j >> 5] >> (2*(j & 31))) & 3;
			if (k != 0) acc = v_xor(acc, c[j][k-1]);
		}
		y[i] = v_xor(y[i], acc);
	}
}

/*------------------------------------------------------------------------*/

/* thanks to Patrick Stach for ideas on this */

#define MAX_OUTER_THREADS 256

__global__ void
lanczos_kernel_outer_prod(v_t *x, v_t *y,
			v_t *xy, uint32 n) 
{
	uint32 i, w_x, w_y;
	uint32 num_threads = gridDim.x * blockDim.x;
	uint32 grid_id = blockIdx.x * blockDim.x + threadIdx.x;
	uint32 block_id = threadIdx.x;
	__shared__ uint64 scratch[3 * MAX_OUTER_THREADS];

	for (w_x = 0; w_x < VWORDS; w_x++) {
		for (w_y = 0; w_y < VWORDS; w_y++) {
			uint64 *s = scratch + (block_id & ~0x1f);
			scratch[block_id + 0*MAX_OUTER_THREADS] = 0;
			scratch[block_id + 1*MAX_OUTER_THREADS] = 0;
			scratch[block_id + 2*MAX_OUTER_THREADS] = 0;

			for (i = grid_id; i < n; i += num_threads) {
				uint32 j;
				uint32 k = block_id & 0x1f;
				uint64 xi = x[i].w[w_x];
				uint64 yi = y[i].w[w_y];

				if (k != 0)
					xi = (xi >> (2 * k)) | (xi << (64 - (2 * k)));

#pragma unroll
				for (j = 0; j < 32; j++) {
					uint32 off = bfe(xi, 2 * j, 2);
					uint64 tmp = yi;

					if (off == 0) {
						tmp = 0;
						off = 1;
					}

					s[((k + j) & 0x1f) + 
						MAX_OUTER_THREADS * (off - 1)] ^= tmp;
				}
			}

			s = scratch + block_id;
			__syncthreads();
			s[0*MAX_OUTER_THREADS] ^= s[2*MAX_OUTER_THREADS];
			s[1*MAX_OUTER_THREADS] ^= s[2*MAX_OUTER_THREADS];
			__syncthreads();

			for (i = MAX_OUTER_THREADS / 2; i >= 32; i >>= 1) {
				if (block_id < i) {
					s[0*MAX_OUTER_THREADS] ^= s[0*MAX_OUTER_THREADS + i];
					s[1*MAX_OUTER_THREADS] ^= s[1*MAX_OUTER_THREADS + i];
				}
				__syncthreads();
			}

			if (block_id < 32) {
				uint64 res = scratch[block_id];
				i = 2 * block_id;
				atomicXor(&xy[64 * w_x + i].w[w_y], res);
			}
			else if (block_id < 64) {
				uint64 res = scratch[MAX_OUTER_THREADS + block_id - 32];
				i = 2 * block_id - 64 + 1;
				atomicXor(&xy[64 * w_x + i].w[w_y], res);
			}
			__syncthreads();
		}
	}
}

#ifdef __cplusplus
}
#endif
