#include "common.h"
#include <immintrin.h>

#ifndef _POLY_COMMON_H_
#define _POLY_COMMON_H_

#define USE_SSE2

typedef struct 
{
	//read/write data inputs
	uint32_t *numptr_n;
	uint32_t *numptr_p;
	uint32_t *sliceptr_n;
	uint32_t *sliceptr_p;
	uint32_t *update_data_prime;
	int *update_data_root1;
	int *update_data_root2;
	uint8_t *update_data_logp;
	lp_bucket *lp_bucket_p;
	int *ptr;

	//read only inputs:
	uint32_t large_B;
	uint32_t B;
	uint32_t interval;
	int numblocks;

	//read/write words
	uint32_t bound_val;
	int bound_index;
	int check_bound;
	uint8_t logp;

	uint32_t intervalm1;
    uint32_t filler;
    uint32_t *scratch;

} polysieve_t;

typedef struct 
{
	uint16_t *first_r1;
	uint16_t *first_r2;
	uint16_t *fbp1;
	uint16_t *fbp2;
	uint16_t *fbn1;
	uint16_t *fbn2;
	uint16_t *primes;
	uint16_t *updates;
	uint32_t start;
	uint32_t stop;
} small_update_t;

#define CHECK_NEW_SLICE(j)									\
	if (j >= check_bound)							\
	{														\
		room = 0;											\
        /* find the most filled bucket */ \
		for (k=0;k<numblocks;k++)							\
		{													\
			if (*(numptr_p + k) > room)						\
				room = *(numptr_p + k);						\
			if (*(numptr_n + k) > room)						\
				room = *(numptr_n + k);						\
		}													\
		room = BUCKET_ALLOC - room;							\
        /* if it is filled close to the allocation, start recording in a new set of buckets */ \
		if (room < 32)										\
		{													\
			logp = update_data.logp[j];						\
			slicelogp_ptr[bound_index] = logp;			\
			bound_index++;									\
			slicebound_ptr[bound_index] = j;		\
			bound_val = j;									\
			sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
			sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
			numptr_p += (numblocks << 1);							\
			numptr_n += (numblocks << 1);							\
			check_bound += BUCKET_ALLOC >> 1;					\
		}													\
		else												\
			check_bound += room >> 1;						\
	}										\
	else if ((j - bound_val) >= 65536)		\
	{										\
		slicelogp_ptr[bound_index] = logp;			\
		bound_index++;									\
		slicebound_ptr[bound_index] = j;		\
		bound_val = j;									\
		sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
		sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
		numptr_p += (numblocks << 1);							\
		numptr_n += (numblocks << 1);							\
		check_bound += BUCKET_ALLOC >> 1;					\
	}

#define CHECK_NEW_SLICE_BATCH(j)									\
	if (j >= check_bound)							\
        	{														\
		room = 0;											\
		for (k = 0; k < numblocks; k++) { \
            for (p = 0; p < dconf->poly_batchsize; p++) { \
                if (numptr_p[k + p * poly_offset] > room) { \
                    room = numptr_p[k + p * poly_offset]; \
                                                                } } } \
		room = BUCKET_ALLOC - room;							\
		if (room < 32)										\
                		{													\
			logp = update_data.logp[j];						\
			lp_bucket_p->logp[bound_index] = logp;			\
			bound_index++;									\
			lp_bucket_p->fb_bounds[bound_index] = j;		\
			bound_val = j;									\
			sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
			sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
			numptr_p += (numblocks << 1);							\
			numptr_n += (numblocks << 1);							\
			check_bound += BUCKET_ALLOC >> 1;					\
                		}													\
                        		else												\
			check_bound += room >> 1;						\
        	}										\
            	else if ((j - bound_val) >= 65536)		\
    	{										\
		lp_bucket_p->logp[bound_index] = logp;			\
		bound_index++;									\
		lp_bucket_p->fb_bounds[bound_index] = j;		\
		bound_val = j;									\
		sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
		sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
		numptr_p += (numblocks << 1);							\
		numptr_n += (numblocks << 1);							\
		check_bound += BUCKET_ALLOC >> 1;					\
    	}

#define CHECK_NEW_SLICE_BATCH_2(j)									\
	if (j >= cb[p])							\
	{														\
		room = 0;											\
        /* find the most filled bucket */ \
		for (k = 0; k < numblocks; k++)							\
		{													\
			if (numptr_p[k] > room)						\
				room = numptr_p[k];						\
			if (numptr_n[k] > room)						\
				room = numptr_n[k];						\
		}													\
		room = BUCKET_ALLOC - room;							\
        /* if it is filled close to the allocation, start recording in a new set of buckets */ \
		if (room < 32)										\
		{													\
			logp = update_data.logp[j];						\
			slicelogp_ptr[p * lp_bucket_p->alloc_slices + bi[p]] = logp;			\
			bi[p]++;									\
            /*printf("new slice bound = %u for poly %d\n", j, p); */\
			slicebound_ptr[p * lp_bucket_p->alloc_slices + bi[p]] = j;		\
			bv[p] = j;									\
			sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
			sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
			numptr_p += (numblocks << 1);							\
			numptr_n += (numblocks << 1);							\
			cb[p] += BUCKET_ALLOC >> 1;					\
		}													\
		else												\
			cb[p] += room >> 1;						\
	}										\
	else if ((j - bv[p]) >= 65536)		\
	{										\
		slicelogp_ptr[p * lp_bucket_p->alloc_slices + bi[p]] = logp;			\
		bi[p]++;									\
		slicebound_ptr[p * lp_bucket_p->alloc_slices + bi[p]] = j;		\
		bv[p] = j;									\
		sliceptr_p += (numblocks << (BUCKET_BITS + 1));		\
		sliceptr_n += (numblocks << (BUCKET_BITS + 1));		\
		numptr_p += (numblocks << 1);							\
		numptr_n += (numblocks << 1);							\
		cb[p] += BUCKET_ALLOC >> 1;					\
	}


#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)

	#define COMPUTE_4_PROOTS(j)								\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				ptrs = _mm_load_si128((__m128i *)(&rootupdates[(v-1) * bound + j])); \
				root1s = _mm_load_si128((__m128i *)(update_data.firstroots1 + j)); \
				root1s = _mm_sub_epi32(root1s, ptrs); 	 					/* root1 -= ptr */ \
				root2s = _mm_load_si128((__m128i *)(update_data.firstroots2 + j)); \
				root2s = _mm_sub_epi32(root2s, ptrs); 	 					/* root2 -= ptr */ \
				tmp1 = _mm_xor_si128(tmp1, tmp1); 							/* zero xmm4 */ \
				tmp2 = _mm_xor_si128(tmp2, tmp2);							/* zero xmm5 */ \
				primes = _mm_load_si128((__m128i *)(update_data.prime + j)); \
				tmp1 = _mm_cmpgt_epi32(tmp1,root1s); 						/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
				tmp2 = _mm_cmpgt_epi32(tmp2,root2s); 						/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
				tmp1 = _mm_and_si128(tmp1, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				tmp2 = _mm_and_si128(tmp2, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				root1s = _mm_add_epi32(root1s, tmp1); 						/* selectively add back prime (modular subtract) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots1 + j),root1s);		/* save new root1 values */ \
				root2s = _mm_add_epi32(root2s, tmp2); 						/* selectively add back prime (modular subtract) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots2 + j),root2s); 		/* save new root2 values */ \
			} while (0);

	#define COMPUTE_4_NROOTS(j)								\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				ptrs = _mm_load_si128((__m128i *)(&rootupdates[(v-1) * bound + j])); \
				root1s = _mm_load_si128((__m128i *)(update_data.firstroots1 + j)); \
				root1s = _mm_add_epi32(root1s, ptrs); 	 					/* root1 += ptr */ \
				root2s = _mm_load_si128((__m128i *)(update_data.firstroots2 + j)); \
				root2s = _mm_add_epi32(root2s, ptrs); 	 					/* root2 += ptr */ \
				tmp1 = _mm_shuffle_epi32(root1s, 0xe4); 					/* copy root1 to xmm4 */ \
				tmp2 = _mm_shuffle_epi32(root2s, 0xe4);						/* copy root2 to xmm5 */ \
				primes = _mm_load_si128((__m128i *)(update_data.prime + j)); \
				tmp1 = _mm_cmpgt_epi32(tmp1,primes); 						/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
				tmp2 = _mm_cmpgt_epi32(tmp2,primes); 						/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
				tmp1 = _mm_and_si128(tmp1, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				tmp2 = _mm_and_si128(tmp2, primes); 						/* copy prime to overflow locations (are set to 1) */ \
				root1s = _mm_sub_epi32(root1s, tmp1); 						/* selectively sub back prime (modular addition) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots1 + j),root1s);		/* save new root1 values */ \
				root2s = _mm_sub_epi32(root2s, tmp2); 						/* selectively sub back prime (modular addition) */ \
				_mm_store_si128((__m128i *)(update_data.firstroots2 + j),root2s); 		/* save new root2 values */ \
			} while (0);

	#define COMPUTE_8X_SMALL_PROOTS	\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				for (j = h.start; j < h.stop; j += 8) \
				{	\
					ptrs = _mm_load_si128((__m128i *)(h.updates + j)); \
					root1s = _mm_load_si128((__m128i *)(h.first_r1 + j)); \
					root2s = _mm_load_si128((__m128i *)(h.first_r2 + j)); \
					primes = _mm_load_si128((__m128i *)(h.primes + j)); \
					root1s = _mm_sub_epi16(root1s, ptrs); \
					root2s = _mm_sub_epi16(root2s, ptrs); \
					tmp1 = _mm_xor_si128(tmp1, tmp1); \
					tmp2 = _mm_xor_si128(tmp2, tmp2); \
					tmp1 = _mm_cmpgt_epi16(tmp1,root1s);	\
					tmp2 = _mm_cmpgt_epi16(tmp2,root2s);	\
					tmp1 = _mm_and_si128(tmp1, primes);	\
					tmp2 = _mm_and_si128(tmp2, primes);	\
					root1s = _mm_add_epi16(root1s, tmp1);	\
					root2s = _mm_add_epi16(root2s, tmp2);	\
					tmp1 = root2s;	\
					tmp1 = _mm_max_epi16(tmp1, root1s); \
					root2s = _mm_min_epi16(root2s, root1s); \
					tmp2 = primes; \
					_mm_store_si128((__m128i *)(h.fbp1 + j), root2s);	\
					primes = _mm_sub_epi16(primes, root2s); \
					_mm_store_si128((__m128i *)(h.first_r1 + j), root2s);	\
					_mm_store_si128((__m128i *)(h.fbp2 + j), tmp1);	\
					tmp2 = _mm_sub_epi16(tmp2, tmp1); \
					_mm_store_si128((__m128i *)(h.first_r2 + j), tmp1);	\
					_mm_store_si128((__m128i *)(h.fbn1 + j), tmp2);	\
					_mm_store_si128((__m128i *)(h.fbn2 + j), primes);	\
				}	\
			} while(0);

	#define COMPUTE_8X_SMALL_NROOTS	\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				for (j = h.start; j < h.stop; j += 8) \
				{	\
					ptrs = _mm_load_si128((__m128i *)(h.updates + j)); \
					root1s = _mm_load_si128((__m128i *)(h.first_r1 + j)); \
					root2s = _mm_load_si128((__m128i *)(h.first_r2 + j)); \
					primes = _mm_load_si128((__m128i *)(h.primes + j)); \
					root1s = _mm_add_epi16(root1s, ptrs); \
					root2s = _mm_add_epi16(root2s, ptrs); \
					tmp1 = _mm_xor_si128(tmp1, tmp1); \
					tmp2 = _mm_xor_si128(tmp2, tmp2); \
					root1s = _mm_sub_epi16(root1s, primes); \
					root2s = _mm_sub_epi16(root2s, primes); \
					tmp1 = _mm_cmpgt_epi16(tmp1,root1s);	\
					tmp2 = _mm_cmpgt_epi16(tmp2,root2s);	\
					tmp1 = _mm_and_si128(tmp1, primes);	\
					tmp2 = _mm_and_si128(tmp2, primes);	\
					root1s = _mm_add_epi16(root1s, tmp1);	\
					root2s = _mm_add_epi16(root2s, tmp2);	\
					tmp1 = root2s;	\
					tmp1 = _mm_max_epi16(tmp1, root1s); \
					root2s = _mm_min_epi16(root2s, root1s); \
					tmp2 = primes; \
					_mm_store_si128((__m128i *)(h.fbp1 + j), root2s);	\
					primes = _mm_sub_epi16(primes, root2s); \
					_mm_store_si128((__m128i *)(h.first_r1 + j), root2s);	\
					_mm_store_si128((__m128i *)(h.fbp2 + j), tmp1);	\
					tmp2 = _mm_sub_epi16(tmp2, tmp1); \
					_mm_store_si128((__m128i *)(h.first_r2 + j), tmp1);	\
					_mm_store_si128((__m128i *)(h.fbn1 + j), tmp2);	\
					_mm_store_si128((__m128i *)(h.fbn2 + j), primes);	\
				}	\
			} while(0);

#elif defined(GCC_ASM64X) && !defined(FORCE_GENERIC)

#ifdef USE_SSE2

	#define COMPUTE_8X_SMALL_PROOTS	\
		__asm (	\
			"movq   %0, %%r13 \n\t"	\
			"xorq	%%r15, %%r15 \n\t"	\
			"xorq	%%rax, %%rax \n\t"	\
			"movl   68(%%r13), %%r15d \n\t"	/* r15d = stop */	\
			"movl   64(%%r13), %%eax \n\t"	/* eax = start */	\
			"movq   0(%%r13), %%rbx \n\t"		/* rbx = first_r1 */	\
			"movq   8(%%r13), %%rcx \n\t"		/* rcx = first_r2 */	\
			"movq   48(%%r13), %%rdx \n\t"	/* rdx = primes */	\
			"movq   56(%%r13), %%r8 \n\t"		/* r8 = updates */	\
			"movq   16(%%r13), %%r9 \n\t"		/* r9 = fbp1 */	\
			"movq   24(%%r13), %%r10 \n\t"	/* r10 = fbp2 */	\
			"movq   32(%%r13), %%r11 \n\t"	/* r11 = fbn1 */	\
			"movq   40(%%r13), %%r12 \n\t"	/* r12 = fbn2 */	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the P side */	\
			"movdqa	(%%r8, %%rax, 2), %%xmm3 \n\t"			/* xmm3 = ptr */	\
			"movdqa (%%rbx, %%rax, 2), %%xmm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"movdqa (%%rcx, %%rax, 2), %%xmm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"psubw	%%xmm3, %%xmm1 \n\t"					/* root1 -= ptr */	\
			"psubw	%%xmm3, %%xmm2 \n\t"					/* root2 -= ptr */	\
			"movdqa (%%rdx, %%rax, 2), %%xmm0 \n\t"			/* xmm0 = next 8 primes */	\
			"pxor	%%xmm4, %%xmm4 \n\t"					/* zero xmm4 */	\
			"pxor	%%xmm5, %%xmm5 \n\t"					/* zero xmm5 */	\
			"pcmpgtw	%%xmm1, %%xmm4 \n\t"				/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */	\
			"pcmpgtw	%%xmm2, %%xmm5 \n\t"				/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */	\
			"pand	%%xmm0, %%xmm4 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"pand	%%xmm0, %%xmm5 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"paddw	%%xmm4, %%xmm1 \n\t"					/* selectively add back prime (modular subtract) */	\
			"paddw	%%xmm5, %%xmm2 \n\t"					/* selectively add back prime (modular subtract) */	\
			"movdqa %%xmm2, %%xmm5 \n\t"					/* xmm5 = root2 copy */	\
			"pmaxsw	%%xmm1, %%xmm5 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"pminsw	%%xmm1, %%xmm2 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
															/* now copy results to appropriate data structures */	\
			"movdqa	%%xmm0, %%xmm4 \n\t"					/* copy primes */	\
															/* root1p always gets the smaller roots (LT) */	\
			"movdqa	%%xmm2, (%%r9, %%rax, 2) \n\t"			/* update root1p */	\
			"psubw	%%xmm2, %%xmm0 \n\t"					/* prime - LT roots */	\
			"movdqa	%%xmm2, (%%rbx, %%rax, 2) \n\t"			/* update firstroots1 */	\
															/* root2p always gets the bigger roots (GT) */	\
			"movdqa	%%xmm5, (%%r10, %%rax, 2) \n\t"			/* update root2p */	\
			"psubw	%%xmm5, %%xmm4 \n\t"					/* prime - GT roots */	\
			"movdqa	%%xmm5, (%%rcx, %%rax, 2) \n\t"			/* update firstroots2 */	\
															/* root1n always gets prime - bigger roots (LT) */	\
			"movdqa	%%xmm4, (%%r11, %%rax, 2) \n\t"			/* update root1n */	\
															/* root2n always gets prime - smaller roots (GT) */	\
			"movdqa	%%xmm0, (%%r12, %%rax, 2) \n\t"			/* update root2n */	\
			"addl	$8, %%eax \n\t"	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "g"(&h)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "rax", "rbx", "rcx", "rdx",	\
			"r8", "r9", "r10", "r11", "r12", "r13", "r15", "cc", "memory");

	#define COMPUTE_8X_SMALL_NROOTS	\
		ASM_G (	\
			"movq   %0, %%rsi \n\t"	\
			"xorq	%%r15, %%r15 \n\t"	\
			"xorq	%%rax, %%rax \n\t"	\
			"movl   68(%%rsi), %%r15d \n\t"	/* r15d = stop */	\
			"movl   64(%%rsi), %%eax \n\t"	/* eax = start */	\
			"movq   0(%%rsi), %%rbx \n\t"		/* rbx = first_r1 */	\
			"movq   8(%%rsi), %%rcx \n\t"		/* rcx = first_r2 */	\
			"movq   48(%%rsi), %%rdx \n\t"	/* rdx = primes */	\
			"movq   56(%%rsi), %%r8 \n\t"		/* r8 = updates */	\
			"movq   16(%%rsi), %%r9 \n\t"		/* r9 = fbp1 */	\
			"movq   24(%%rsi), %%r10 \n\t"	/* r10 = fbp2 */	\
			"movq   32(%%rsi), %%r11 \n\t"	/* r11 = fbn1 */	\
			"movq   40(%%rsi), %%r12 \n\t"	/* r12 = fbn2 */	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the N side */	\
			"movdqa (%%r8, %%rax, 2), %%xmm3 \n\t"			/* xmm3 = next 8 updates */	\
			"movdqa (%%rbx, %%rax, 2), %%xmm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"movdqa (%%rcx, %%rax, 2), %%xmm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"paddw	%%xmm3, %%xmm1 \n\t"					/* root1 += ptr */	\
			"paddw	%%xmm3, %%xmm2 \n\t"					/* root2 += ptr */	\
			"movdqa (%%rdx, %%rax, 2), %%xmm0 \n\t"			/* xmm0 = next 8 primes */	\
			"pxor	%%xmm4, %%xmm4 \n\t"					/* zero xmm4 */	\
			"pxor	%%xmm5, %%xmm5 \n\t"					/* zero xmm5 */	\
			"psubw	%%xmm0, %%xmm1 \n\t"					/* (modular subtract) */		\
			"psubw	%%xmm0, %%xmm2 \n\t"					/* (modular subtract) */	\
			"pcmpgtw	%%xmm1, %%xmm4 \n\t"				/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */	\
			"pcmpgtw	%%xmm2, %%xmm5 \n\t"				/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */	\
			"pand	%%xmm0, %%xmm4 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"pand	%%xmm0, %%xmm5 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"paddw	%%xmm4, %%xmm1 \n\t"					/* selectively add back prime (modular subtract) */	\
			"paddw	%%xmm5, %%xmm2 \n\t"					/* selectively add back prime (modular subtract) */	\
			"movdqa %%xmm2, %%xmm5 \n\t"					/* xmm5 = root2 copy */	\
			"pmaxsw	%%xmm1, %%xmm5 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"pminsw	%%xmm1, %%xmm2 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
															/* now copy results to appropriate data structures */	\
			"movdqa	%%xmm0, %%xmm4 \n\t"					/* copy primes */	\
															/* root1p always gets the smaller roots (LT) */	\
			"movdqa	%%xmm2, (%%r9, %%rax, 2) \n\t"			/* update root1p */	\
			"psubw	%%xmm2, %%xmm0 \n\t"					/* prime - LT roots */	\
			"movdqa	%%xmm2, (%%rbx, %%rax, 2) \n\t"			/* update firstroots1 */	\
															/* root2p always gets the bigger roots (GT) */	\
			"movdqa	%%xmm5, (%%r10, %%rax, 2) \n\t"			/* update root2p */	\
			"psubw	%%xmm5, %%xmm4 \n\t"					/* prime - GT roots */	\
			"movdqa	%%xmm5, (%%rcx, %%rax, 2) \n\t"			/* update firstroots2 */	\
															/* root1n always gets prime - bigger roots (LT) */	\
			"movdqa	%%xmm4, (%%r11, %%rax, 2) \n\t"			/* update root1n */	\
															/* root2n always gets prime - smaller roots (GT) */	\
			"movdqa	%%xmm0, (%%r12, %%rax, 2) \n\t"			/* update root2n */	\
			"addl	$8, %%eax \n\t"	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "g"(&h)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "rax", "rsi", "rbx", "rcx", "rdx",	\
			"r8", "r9", "r10", "r11", "r12", "r15", "cc", "memory");

	#define COMPUTE_4_PROOTS(j)								\
		ASM_G (											\
			"movdqa (%%rax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%rcx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%rdx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%rbx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%rcx) \n\t"			/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%rdx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	

	#define COMPUTE_4_NROOTS(j)								\
		ASM_G (											\
			"movdqa (%%rax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%rcx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%rdx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%rbx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%rcx) \n\t"			/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%rdx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	


#endif

#elif defined(GCC_ASM32X)

	#define COMPUTE_4_PROOTS(j)								\
		ASM_G (											\
			"movdqa (%%eax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%ecx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%edx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%ebx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%ecx) \n\t"			/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%edx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	

	#define COMPUTE_4_NROOTS(j)								\
		ASM_G (											\
			"movdqa (%%eax), %%xmm3 \n\t"			/* xmm3 = next 4 values of rootupdates */ \
			"movdqa (%%ecx), %%xmm1 \n\t"			/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%edx), %%xmm2 \n\t"			/* xmm2 = next 4 values of root2 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%ebx), %%xmm0 \n\t"			/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%ecx) \n\t"			/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%edx) \n\t"			/* save new root2 values */ \
			: \
			: "a"(&rootupdates[(v-1) * bound + j]), "b"(update_data.prime + j), "c"(update_data.firstroots1 + j), "d"(update_data.firstroots2 + j) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "cc");	


#endif

#if defined (GCC_ASM64X) && !defined(FORCE_GENERIC)
	#define CHECK_NEW_SLICE_ASM \
		"cmpl   104(%%rsi),%%r15d	\n\t"		/* compare j with check_bound */ \
			/* note this is the counter j, not the byte offset j */ \
		"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
			/* else, this is the "else-if" check */ \
		"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
		"subl   96(%%rsi),%%ebx \n\t"			/* ebx = j - bound_val */ \
		"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
		"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
			/* now we are in the else-if block of CHECK_NEW_SLICE */ \
		"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
		"movl   100(%%rsi),%%edx \n\t"		/* move bound_index into rdx */ \
		"movq   64(%%rsi),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
		"movq	16(%%r9),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
		"movq	56(%%rsi),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
		"movzbl (%%r14,%%r15),%%ebx \n\t"		/* bring in logp */ \
		"movb	%%bl, 108(%%rsi) \n\t"		/* shove logp into output */ \
		"movb   %%bl,(%%r8,%%rdx) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
		"incq   %%rdx \n\t"						/* increment bound_index locally */ \
		"movl   %%edx,100(%%rsi) \n\t"		/* copy bound_index back to structure */ \
		"movq	8(%%r9),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
		"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
			/* note this is the counter j, not the byte offset j */ \
		"movl   %%r15d,96(%%rsi) \n\t"		/* bound_val = j */ \
		"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
		"movl   92(%%rsi),%%ebx \n\t"			/* put numblocks into ebx */ \
		"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
		"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
		"addq   %%rbx,8(%%rsi) \n\t"			/* numptr_p += (numblocks << 1) */ \
		"addq   %%rbx,0(%%rsi) \n\t"			/* numptr_n += (numblocks << 1) */ \
		"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
			/* note also, this works because we've already left shifted by 1 */ \
		"addq   %%rbx,24(%%rsi) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
		"addq   %%rbx,16(%%rsi) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
		"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
		"cmp	%%rax,%%rax \n\t"				/* force jump */ \
		"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
		"1:		\n\t"									\
			/* now we are in the if block of CHECK_NEW_SLICE */ \
		"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
		"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
		"cmpl   92(%%rsi),%%ebx \n\t"			/* compare with numblocks */ \
		"jae    3f \n\t"						/* jump past loop if condition met */ \
			/* condition not met, put a couple things in registers */ \
		"movq	8(%%rsi),%%r10 \n\t"			/* numptr_p into r10 */ \
		"movq	0(%%rsi),%%r11 \n\t"			/* numptr_n into r11 */ \
		"5:		\n\t"							\
			/* now we are in the room loop */ \
			/* room is in register ecx */ \
		"movl   (%%r10,%%rbx,4),%%eax \n\t"		/* value at numptr_p + k */ \
		"movl   (%%r11,%%rbx,4),%%edx \n\t"		/* value at numptr_n + k */ \
		"cmpl   %%ecx,%%eax \n\t"				/* *(numptr_p + k) > room ? */ \
		"cmova  %%eax,%%ecx \n\t"				/* new value of room if so */ \
		"cmpl   %%ecx,%%edx \n\t"				/* *(numptr_p + k) > room ? */ \
		"cmova  %%edx,%%ecx \n\t"				/* new value of room if so */ \
		"incq   %%rbx \n\t"						/* increment counter */ \
		"cmpl   92(%%rsi),%%ebx \n\t"			/* compare to numblocks */ \
		"jl     5b \n\t"						/* iterate loop if condition met */ \
		"3:		\n\t"							\
		"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
		"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
		"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
		"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
		"jle    4f \n\t"						/* jump if less than */ \
		"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
		"addl   %%ebx,104(%%rsi) \n\t"		/* add (room >> 1) to check_bound */ \
		"cmpq	%%rax,%%rax \n\t"				/* force jump */ \
		"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
		"4:		\n\t"							\
			/* now we are inside the (room < 2) block */ \
		"xorq	%%rax, %%rax \n\t" \
		"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
		"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
		"xorq	%%rdx, %%rdx \n\t" \
		"movl   100(%%rsi),%%edx \n\t"		/* move bound_index into rdx */ \
		"movq   64(%%rsi),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
		"movq	16(%%r9),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
		"movq	56(%%rsi),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
		"movzbl (%%r14,%%r15),%%ebx \n\t"		/* bring in logp */ \
		"movb	%%bl, 108(%%rsi) \n\t"		/* shove logp into output */ \
		"movb   %%bl,(%%r8,%%rdx) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
		"incq   %%rdx \n\t"						/* increment bound_index locally */ \
		"movl   %%edx,100(%%rsi) \n\t"		/* copy bound_index back to structure */ \
		"movq	8(%%r9),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
		"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
			/* note this is the counter j, not the byte offset j */ \
		"movl   %%r15d,96(%%rsi) \n\t"		/* bound_val = j */ \
		"xorq	%%rbx, %%rbx \n\t" \
		"movl   92(%%rsi),%%ebx \n\t"			/* put numblocks into ebx */ \
		"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
		"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
		"addq   %%rbx,8(%%rsi) \n\t"			/* numptr_p += (numblocks << 1) */ \
		"addq   %%rbx,0(%%rsi) \n\t"			/* numptr_n += (numblocks << 1) */ \
		"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
			/* note also, this works because we've already left shifted by 1 */ \
		"addq   %%rbx,24(%%rsi) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
		"addq   %%rbx,16(%%rsi) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
		"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
		"2:		\n\t"

#endif


#if defined(MSC_ASM32A) && !defined(FORCE_GENERIC) && !defined(__INTEL_COMPILER)
	#define COMPUTE_NEXT_ROOTS_P	\
	do {	\
		uint32_t update = *ptr;	\
		ASM_M {					\
			ASM_M xor ebx, ebx	\
			ASM_M xor ecx, ecx	\
			ASM_M mov eax, root1	\
			ASM_M mov edx, root2	\
			ASM_M sub eax, update	\
			ASM_M cmovc ebx, prime	\
			ASM_M sub edx, update	\
			ASM_M cmovc ecx, prime	\
			ASM_M add eax, ebx	\
			ASM_M add edx, ecx	\
			ASM_M mov root1, eax	\
			ASM_M mov root2, edx}	\
		} while (0);		

	#define COMPUTE_NEXT_ROOTS_N	\
	do {	\
		uint32_t update = *ptr;	\
		ASM_M {					\
			ASM_M mov eax, root1	\
			ASM_M mov edx, root2	\
			ASM_M mov ebx, eax		\
			ASM_M add ebx, update	\
			ASM_M mov ecx, edx		\
			ASM_M add ecx, update	\
			ASM_M sub eax, prime	\
			ASM_M sub edx, prime	\
			ASM_M add eax, update	\
			ASM_M cmovae eax, ebx	\
			ASM_M add edx, update	\
			ASM_M cmovae edx, ecx	\
			ASM_M mov root1, eax	\
			ASM_M mov root2, edx}	\
		} while (0);

	#define COMPUTE_FIRST_ROOTS	\
	do {	\
		ASM_M {					\
			ASM_M xor ebx, ebx	\
			ASM_M xor ecx, ecx	\
			ASM_M mov eax, root1	\
			ASM_M mov edx, root2	\
			ASM_M sub eax, bmodp	\
			ASM_M cmovc ebx, prime	\
			ASM_M sub edx, bmodp	\
			ASM_M cmovc ecx, prime	\
			ASM_M add eax, ebx	\
			ASM_M add edx, ecx	\
			ASM_M mov root1, eax	\
			ASM_M mov root2, edx}	\
		} while (0);
	

#elif defined(GCC_ASM64X) && !defined(FORCE_GENERIC)

#ifdef _WIN32
#define ASM_ ASM_M
#else
#define ASM_ ASM_G
#endif

	#define COMPUTE_NEXT_ROOTS_P						\
		ASM_ (											\
			"xorl %%r8d, %%r8d		\n\t"	/*r8d = 0*/	\
			"xorl %%r9d, %%r9d		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - ptr*/	\
			"cmovc %3, %%r8d		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - ptr*/	\
			"cmovc %3, %%r9d		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%r8d, %%eax		\n\t"		\
			"addl %%r9d, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "g"(*ptr), "g"(prime)		\
			: "r8", "r9", "cc");

	#define COMPUTE_NEXT_ROOTS_N		\
		ASM_ (							\
			"movl %%eax, %%r8d		\n\t"	\
			"addl %2, %%r8d			\n\t"	/*r8d = root1 + ptr*/		\
			"movl %%edx, %%r9d		\n\t"								\
			"addl %2, %%r9d			\n\t"	/*r9d = root2 + ptr*/		\
			"subl %3, %%eax			\n\t"	/*root1 = root1 - prime*/	\
			"subl %3, %%edx			\n\t"	/*root2 = root2 - prime*/	\
			"addl %2, %%eax			\n\t"	/*root1 + ptr*/				\
			"cmovae %%r8d, %%eax	\n\t"	/*other caluclation if no overflow*/	\
			"addl %2, %%edx			\n\t"	/*root2 + ptr*/							\
			"cmovae %%r9d, %%edx	\n\t"	/*other caluclation if no overflow*/	\
			: "+a"(root1), "+d"(root2)		\
			: "g"(*ptr), "g"(prime)			\
			: "r8", "r9", "cc");

	#define COMPUTE_FIRST_ROOTS			\
		ASM_ (											\
			"xorl %%r8d, %%r8d		\n\t"	/*r8d = 0*/	\
			"xorl %%r9d, %%r9d		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - bmodp*/	\
			"cmovc %3, %%r8d		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - bmodp*/	\
			"cmovc %3, %%r9d		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%r8d, %%eax		\n\t"		\
			"addl %%r9d, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "r"(bmodp), "g"(prime)		\
			: "r8", "r9", "cc");

#elif defined(GCC_ASM32X) && !defined(FORCE_GENERIC)

	#define COMPUTE_NEXT_ROOTS_P						\
		ASM_G (											\
			"xorl %%ecx, %%ecx		\n\t"	/*r8d = 0*/	\
			"xorl %%edi, %%edi		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - ptr*/	\
			"cmovc %3, %%ecx		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - ptr*/	\
			"cmovc %3, %%edi		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%ecx, %%eax		\n\t"		\
			"addl %%edi, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "g"(*ptr), "g"(prime)		\
			: "ecx", "edi", "cc");	

	#define COMPUTE_NEXT_ROOTS_N		\
		ASM_G (							\
			"movl %%eax, %%ecx		\n\t"	\
			"addl %2, %%ecx			\n\t"	/*r8d = root1 + ptr*/		\
			"movl %%edx, %%edi		\n\t"								\
			"addl %2, %%edi			\n\t"	/*r9d = root2 + ptr*/		\
			"subl %3, %%eax			\n\t"	/*root1 = root1 - prime*/	\
			"subl %3, %%edx			\n\t"	/*root2 = root2 - prime*/	\
			"addl %2, %%eax			\n\t"	/*root1 + ptr*/				\
			"cmovae %%ecx, %%eax	\n\t"	/*other caluclation if no overflow*/	\
			"addl %2, %%edx			\n\t"	/*root2 + ptr*/							\
			"cmovae %%edi, %%edx	\n\t"	/*other caluclation if no overflow*/	\
			: "+a"(root1), "+d"(root2)		\
			: "g"(*ptr), "g"(prime)			\
			: "ecx", "edi", "cc");

	#define COMPUTE_FIRST_ROOTS			\
		ASM_G (											\
			"xorl %%ecx, %%ecx		\n\t"	/*r8d = 0*/	\
			"xorl %%edi, %%edi		\n\t"	/*r9d = 0*/	\
			"subl %2, %%eax			\n\t"	/*root1 - bmodp*/	\
			"cmovc %3, %%ecx		\n\t"	/*prime into r8 if overflow*/	\
			"subl %2, %%edx			\n\t"	/*root2 - bmodp*/	\
			"cmovc %3, %%edi		\n\t"	/*prime into r9 if overflow*/	\
			"addl %%ecx, %%eax		\n\t"		\
			"addl %%edi, %%edx		\n\t"		\
			: "+a"(root1), "+d"(root2)			\
			: "r"(bmodp), "g"(prime)		\
			: "ecx", "edi", "cc");	


#else

	#define COMPUTE_NEXT_ROOTS_P		\
		root1 = (int)root1 - *ptr;		\
		root2 = (int)root2 - *ptr;		\
		root1 += ((root1 >> 31) * prime);			\
		root2 += ((root2 >> 31) * prime);	

	#define COMPUTE_NEXT_ROOTS_N		\
		root1 = (int)root1 + *ptr;		\
		root2 = (int)root2 + *ptr;		\
		root1 -= ((root1 >= prime) * prime);	\
		root2 -= ((root2 >= prime) * prime);	

	#define COMPUTE_FIRST_ROOTS		\
		root1 = (int)root1 - bmodp;		\
		root2 = (int)root2 - bmodp;		\
		if (root1 < 0) root1 += prime;			\
		if (root2 < 0) root2 += prime;

#define COMPUTE_NEXT_ROOTS_BATCH(i) \
        if (gray[numB + i] > 0) { \
            root1 = (int)root1 - rootupdates[(nu[numB + i] - 1) * bound + j]; \
            root2 = (int)root2 - rootupdates[(nu[numB + i] - 1) * bound + j]; \
            root1 += ((root1 >> 31) * prime); \
            root2 += ((root2 >> 31) * prime); \
        } else { \
            root1 = (int)root1 + rootupdates[(nu[numB + i] - 1) * bound + j]; \
            root2 = (int)root2 + rootupdates[(nu[numB + i] - 1) * bound + j]; \
            root1 -= ((root1 >= prime) * prime); \
            root2 -= ((root2 >= prime) * prime); \
        }

#define COMPUTE_NEXT_ROOTS_BATCH_P(i) \
        root1 = (int)root1 - rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root2 = (int)root2 - rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root1 += ((root1 >> 31) * prime); \
        root2 += ((root2 >> 31) * prime);

#define COMPUTE_NEXT_ROOTS_BATCH_N(i) \
        root1 = (int)root1 + rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root2 = (int)root2 + rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
        root1 -= ((root1 >= prime) * prime); \
        root2 -= ((root2 >= prime) * prime); \

#endif

#endif

