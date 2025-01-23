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

#ifndef _COMMON_LANCZOS_GPU_LANCZOS_GPU_CORE_H_
#define _COMMON_LANCZOS_GPU_LANCZOS_GPU_CORE_H_

#if defined(__CUDACC__) /*------------- device code -------------*/

typedef short int16;
typedef unsigned short uint16;
typedef int int32;
typedef unsigned int uint32;
typedef unsigned long long uint64;
typedef long long int64;

#if defined(_WIN64) || defined(__LP64__)
	#define PTR_CONSTRAINT(x) "l"(x)
#else
	#define PTR_CONSTRAINT(x) "r"(x)
#endif

#define VWORDS ((VBITS + 63) / 64)

typedef struct {
	uint64 w[VWORDS];
} v_t;

__device__ v_t v_and(v_t a, v_t b) {
	v_t res;

	res.w[0] = a.w[0] & b.w[0];
	#if VWORDS > 1
	res.w[1] = a.w[1] & b.w[1];
	#if VWORDS > 2
	res.w[2] = a.w[2] & b.w[2];
	#if VWORDS > 3
	res.w[3] = a.w[3] & b.w[3];
	#if VWORDS > 4
	res.w[4] = a.w[4] & b.w[4];
	#if VWORDS > 5
	res.w[5] = a.w[5] & b.w[5];
	#if VWORDS > 6
	res.w[6] = a.w[6] & b.w[6];
	#if VWORDS > 7
	res.w[7] = a.w[7] & b.w[7];
	#endif
	#endif
	#endif
	#endif
	#endif
	#endif
	#endif
	return res;
}

__device__ v_t v_xor(v_t a, v_t b) {
	v_t res;

	res.w[0] = a.w[0] ^ b.w[0];
	#if VWORDS > 1
	res.w[1] = a.w[1] ^ b.w[1];
	#if VWORDS > 2
	res.w[2] = a.w[2] ^ b.w[2];
	#if VWORDS > 3
	res.w[3] = a.w[3] ^ b.w[3];
	#if VWORDS > 4
	res.w[4] = a.w[4] ^ b.w[4];
	#if VWORDS > 5
	res.w[5] = a.w[5] ^ b.w[5];
	#if VWORDS > 6
	res.w[6] = a.w[6] ^ b.w[6];
	#if VWORDS > 7
	res.w[7] = a.w[7] ^ b.w[7];
	#endif
	#endif
	#endif
	#endif
	#endif
	#endif
	#endif
	return res;
}

__device__ void v_atomicxor(v_t * a, v_t b) {
	atomicXor(&(a->w[0]), b.w[0]);
	#if VWORDS > 1
	atomicXor(&(a->w[1]), b.w[1]);
	#if VWORDS > 2
	atomicXor(&(a->w[2]), b.w[2]);
	#if VWORDS > 3
	atomicXor(&(a->w[3]), b.w[3]);
	#if VWORDS > 4
	atomicXor(&(a->w[4]), b.w[4]);
	#if VWORDS > 5
	atomicXor(&(a->w[5]), b.w[5]);
	#if VWORDS > 6
	atomicXor(&(a->w[6]), b.w[6]);
	#if VWORDS > 7
	atomicXor(&(a->w[7]), b.w[7]);
	#endif
	#endif
	#endif
	#endif
	#endif
	#endif
	#endif
	return;
}

static const v_t v_zero = {{0}};

__device__ uint32
bfe(uint64 x, uint32 pos, uint32 bits)
{
#if __CUDA_ARCH__ >= 200

	uint32 res;
	uint32 hi = (uint32)(x >> 32);
	uint32 lo = (uint32)x;

	if (pos < 32) {
	       if (pos + bits > 32) {
			res = ((lo >> pos) | (hi << (32 - pos))) &
				((1 << bits) - 1);
	       }
	       else {
			asm("bfe.u32 %0, %1, %2, %3; \n\t"
				: "=r"(res) : "r"(lo), "r"(pos), "r"(bits));
	       }
	}
	else {
		asm("bfe.u32 %0, %1, %2, %3; \n\t"
			: "=r"(res) : "r"(hi), "r"(pos - 32), "r"(bits));
	}

	return res;

#else

	return (uint32)(x >> pos) & ((1 << bits) - 1);
#endif
}

__device__ uint64
load_bypassL1(uint64 *addr)
{
#if __CUDA_ARCH__ >= 200 && VBITS == 64

	uint64 res;

	asm("ld.global.cg.u64 %0, [%1]; \n\t"
		: "=l"(res) : PTR_CONSTRAINT(addr));

	return res;
#else
	return addr[0];
#endif
}

__device__ void
store_bypassL1(uint64 x, uint64 *addr)
{
#if __CUDA_ARCH__ >= 200 && VBITS == 64

	asm("st.global.cg.u64 [%0], %1; \n\t"
		: : PTR_CONSTRAINT(addr), "l"(x));
#else
	addr[0] = x;
#endif
}

#endif  /*--------------------------- device code -----------------*/

#ifdef __cplusplus
extern "C" {
#endif

#define MATMUL_THREADS 256

typedef union {

	uint32 w;

	/* offset is a column offset unless head is nonzero,
	   in which case it is a row offset (i.e. a trigger that
	   a new row is starting) */

	struct {
		uint32 offset : 31;
		uint32 head : 1;
	} d; 
} gpu_entry_idx_t;

#ifdef __cplusplus
}
#endif

#endif /* !_COMMON_LANCZOS_GPU_LANCZOS_GPU_CORE_H_ */
