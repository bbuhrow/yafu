/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"
#include "util.h"
#include "common.h"

typedef struct 
{
	//read/write data inputs
	uint32 *numptr_n;
	uint32 *numptr_p;
	uint32 *sliceptr_n;
	uint32 *sliceptr_p;
	uint32 *update_data_prime;
	int *update_data_root1;
	int *update_data_root2;
	uint8 *update_data_logp;
	lp_bucket *lp_bucket_p;
	int *ptr;

	//read only inputs:
	uint32 large_B;
	uint32 B;
	uint32 interval;
	int numblocks;

	//read/write words
	uint32 bound_val;
	int bound_index;
	int check_bound;
	uint8 logp;

} polysieve_t;

typedef struct 
{
	uint16 *first_r1;
	uint16 *first_r2;
	uint16 *fbp1;
	uint16 *fbp2;
	uint16 *fbn1;
	uint16 *fbn2;
	uint16 *primes;
	uint16 *updates;
	uint32 start;
	uint32 stop;
} small_update_t;

#define FILL_ONE_PRIME_P(i)					\
	if (root1 < interval)					\
	{										\
		bnum = root1 >> BLOCKBITS;			\
		bptr = sliceptr_p +					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
		numptr_p[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = sliceptr_p +					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_p[bnum]++;					\
	}

#define FILL_ONE_PRIME_N(i)						\
	if (root1 < interval)						\
	{											\
		bnum = root1 >> BLOCKBITS;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
		numptr_n[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> BLOCKBITS;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_n[bnum]++;					\
	}			

#define FILL_ONE_PRIME_LOOP_P(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
		numptr_p[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> BLOCKBITS;				\
	}											\
	bnum = root2 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_p[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> BLOCKBITS;				\
	} 

#define FILL_ONE_PRIME_LOOP_N(i)				\
	bnum = root1 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_n +						\
			(bnum << BUCKET_BITS) +				\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & BLOCKSIZEm1)); \
		numptr_n[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> BLOCKBITS;				\
	}											\
	bnum = root2 >> BLOCKBITS;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_n +						\
			(bnum << BUCKET_BITS) +				\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & BLOCKSIZEm1)); \
		numptr_n[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> BLOCKBITS;				\
	} 

#define CHECK_NEW_SLICE(j)									\
	if (j >= check_bound)							\
	{														\
		room = 0;											\
		for (k=0;k<numblocks;k++)							\
		{													\
			if (*(numptr_p + k) > room)						\
				room = *(numptr_p + k);						\
			if (*(numptr_n + k) > room)						\
				room = *(numptr_n + k);						\
		}													\
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


#if defined(_MSC_VER)

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

#endif

#if defined(MSC_ASM32A)
	#define COMPUTE_NEXT_ROOTS_P	\
	do {	\
		uint32 update = *ptr;	\
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
		uint32 update = *ptr;	\
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
	

#elif defined(GCC_ASM64X)

	#define COMPUTE_NEXT_ROOTS_P						\
		ASM_G (											\
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

	#define COMPUTE_8X_SMALL_PROOTS	\
		ASM_G (	\
			"movq   %0, %%rsi \n\t"	\
			"xorq	%%r15, %%r15 \n\t"	\
			"xorq	%%rax, %%rax \n\t"	\
			"movl   68(%%rsi,1), %%r15d \n\t"	/* r15d = stop */	\
			"movl   64(%%rsi,1), %%eax \n\t"	/* eax = start */	\
			"movq   0(%%rsi,1), %%rbx \n\t"		/* rbx = first_r1 */	\
			"movq   8(%%rsi,1), %%rcx \n\t"		/* rcx = first_r2 */	\
			"movq   48(%%rsi,1), %%rdx \n\t"	/* rdx = primes */	\
			"movq   56(%%rsi,1), %%r8 \n\t"		/* r8 = updates */	\
			"movq   16(%%rsi,1), %%r9 \n\t"		/* r9 = fbp1 */	\
			"movq   24(%%rsi,1), %%r10 \n\t"	/* r10 = fbp2 */	\
			"movq   32(%%rsi,1), %%r11 \n\t"	/* r11 = fbn1 */	\
			"movq   40(%%rsi,1), %%r12 \n\t"	/* r12 = fbn2 */	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the P side */	\
			"movdqa	(%%r8, %%rax, 2), %%xmm3 \n\t"			/* xmm3 = ptr */	\
			"movdqa (%%rbx, %%rax, 2), %%xmm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"movdqa (%%rcx, %%rax, 2), %%xmm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"movdqa (%%rdx, %%rax, 2), %%xmm0 \n\t"			/* xmm0 = next 8 primes */	\
			"psubw	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */	\
			"psubw	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */	\
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */	\
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */	\
			"pcmpgtw	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */	\
			"pcmpgtw	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */	\
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */	\
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */	\
			"paddw	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */	\
			"paddw	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */	\
			"movdqa %%xmm2, %%xmm5 \n\t"			/* xmm5 = root2 copy */	\
			"pmaxsw	%%xmm1, %%xmm5 \n\t"		/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"pminsw	%%xmm1, %%xmm2 \n\t"		/* xmm2 = root2 < root1 ? root2 : root1 */	\
			/* now copy results to appropriate data structures */	\
			"movdqa	%%xmm0, %%xmm4 \n\t"			/* copy primes */	\
			/* root1p always gets the smaller roots (LT) */	\
			"movdqa	%%xmm2, (%%r9, %%rax, 2) \n\t"				/* update root1p */	\
			"psubw	%%xmm2, %%xmm0 \n\t"			/* prime - LT roots */	\
			"movdqa	%%xmm2, (%%rbx, %%rax, 2) \n\t"				/* update firstroots1 */	\
			/* root2p always gets the bigger roots (GT) */	\
			"movdqa	%%xmm5, (%%r10, %%rax, 2) \n\t"				/* update root2p */	\
			"psubw	%%xmm5, %%xmm4 \n\t"			/* prime - GT roots */	\
			"movdqa	%%xmm5, (%%rcx, %%rax, 2) \n\t"				/* update firstroots2 */	\
			/* root1n always gets prime - bigger roots (LT) */	\
			"movdqa	%%xmm4, (%%r11, %%rax, 2) \n\t"				/* update root1n */	\
			/* root2n always gets prime - smaller roots (GT) */	\
			"movdqa	%%xmm0, (%%r12, %%rax, 2) \n\t"				/* update root2n */	\
			"addl	$8, %%eax \n\t"	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "g"(&h)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "rax", "rsi", "rbx", "rcx", "rdx",	\
			"r8", "r9", "r10", "r11", "r12", "r15", "cc", "memory");

#define COMPUTE_8X_SMALL_NROOTS	\
		ASM_G (	\
			"movq   %0, %%rsi \n\t"	\
			"xorq	%%r15, %%r15 \n\t"	\
			"xorq	%%rax, %%rax \n\t"	\
			"movl   68(%%rsi,1), %%r15d \n\t"	/* r15d = stop */	\
			"movl   64(%%rsi,1), %%eax \n\t"	/* eax = start */	\
			"movq   0(%%rsi,1), %%rbx \n\t"		/* rbx = first_r1 */	\
			"movq   8(%%rsi,1), %%rcx \n\t"		/* rcx = first_r2 */	\
			"movq   48(%%rsi,1), %%rdx \n\t"	/* rdx = primes */	\
			"movq   56(%%rsi,1), %%r8 \n\t"		/* r8 = updates */	\
			"movq   16(%%rsi,1), %%r9 \n\t"		/* r9 = fbp1 */	\
			"movq   24(%%rsi,1), %%r10 \n\t"	/* r10 = fbp2 */	\
			"movq   32(%%rsi,1), %%r11 \n\t"	/* r11 = fbn1 */	\
			"movq   40(%%rsi,1), %%r12 \n\t"	/* r12 = fbn2 */	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the N side */	\
			"movdqa (%%r8, %%rax, 2), %%xmm3 \n\t"			/* xmm3 = next 8 updates */	\
			"movdqa (%%rdx, %%rax, 2), %%xmm0 \n\t"			/* xmm0 = next 8 primes */	\
			"movdqa (%%rbx, %%rax, 2), %%xmm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"movdqa (%%rcx, %%rax, 2), %%xmm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"paddw	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */	\
			"paddw	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */	\
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */	\
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */	\
			"psubw	%%xmm0, %%xmm1 \n\t"			/* (modular subtract) */		\
			"psubw	%%xmm0, %%xmm2 \n\t"			/* (modular subtract) */	\
			"pcmpgtw	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */	\
			"pcmpgtw	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */	\
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */	\
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */	\
			"paddw	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */	\
			"paddw	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */	\
			"movdqa %%xmm2, %%xmm5 \n\t"			/* xmm5 = root2 copy */	\
			"pmaxsw	%%xmm1, %%xmm5 \n\t"		/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"pminsw	%%xmm1, %%xmm2 \n\t"		/* xmm2 = root2 < root1 ? root2 : root1 */	\
			/* now copy results to appropriate data structures */	\
			"movdqa	%%xmm0, %%xmm4 \n\t"			/* copy primes */	\
			/* root1p always gets the smaller roots (LT) */	\
			"movdqa	%%xmm2, (%%r9, %%rax, 2) \n\t"				/* update root1p */	\
			"psubw	%%xmm2, %%xmm0 \n\t"			/* prime - LT roots */	\
			"movdqa	%%xmm2, (%%rbx, %%rax, 2) \n\t"				/* update firstroots1 */	\
			/* root2p always gets the bigger roots (GT) */	\
			"movdqa	%%xmm5, (%%r10, %%rax, 2) \n\t"				/* update root2p */	\
			"psubw	%%xmm5, %%xmm4 \n\t"			/* prime - GT roots */	\
			"movdqa	%%xmm5, (%%rcx, %%rax, 2) \n\t"				/* update firstroots2 */	\
			/* root1n always gets prime - bigger roots (LT) */	\
			"movdqa	%%xmm4, (%%r11, %%rax, 2) \n\t"				/* update root1n */	\
			/* root2n always gets prime - smaller roots (GT) */	\
			"movdqa	%%xmm0, (%%r12, %%rax, 2) \n\t"				/* update root2n */	\
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

	#define COMPUTE_NEXT_ROOTS_N		\
		ASM_G (							\
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

	// The assembly for putting a prime into a bucket is fairly regular, so we break 
	// macros out to make the inline loop shorter within nextRoots().  The only things
	// that change between root1 and root2 are the input root register.  we use two
	// in order to take advantage of the fractional clock latency of movd, and to
	// break up a dependency bottleneck between movd and the cmpl.

	// macro for adding an element (specified by r8d) to the end of a bucket list
	#define UPDATE_ROOT1(it) \
		"movl   %%r8d,%%ebx \n\t"				/* ebx becomes bnum */ \
		"movl   %%r15d,%%edi \n\t"				/* edi becomes fb offset */ \
		"movl   %%r8d,%%eax \n\t"				/* eax becomes block location */ \
		"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root by blksize  = bnum */ \
		"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root & BLOCKSIZEm1 */ \
		"movl   %%ebx,%%ecx \n\t"				/* ecx becomes bucket address offset */ \
		"addl	$" it ", %%edi \n\t"			/* add iteration number to j */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"subl   %%r12d,%%edi \n\t"				/* j - bound_val */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
		"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */

	// macro for adding an element (specified by r9d) to the end of a bucket list
	#define UPDATE_ROOT2(it) \
		"movl   %%r9d,%%ebx \n\t"				/* ebx becomes bnum */ \
		"movl   %%r15d,%%edi \n\t"				/* edi becomes fb offset */ \
		"movl   %%r9d,%%eax \n\t"				/* eax becomes block location */ \
		"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root by blksize  = bnum */ \
		"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root & BLOCKSIZEm1 */ \
		"movl   %%ebx,%%ecx \n\t"				/* ecx becomes bucket address offset */ \
		"addl	$" it ", %%edi \n\t"			/* add iteration number to j */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"subl   %%r12d,%%edi \n\t"				/* j - bound_val */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
		"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */

	// macro for iteratively adding an element (specified by r8d) to the end of a bucket list
	#define UPDATE_ROOT1_LOOP(it) \
		"1:		\n\t"	 						/* beginning of loop */ \
		"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
		"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
		"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
		"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
		"addl	$" it ", %%edi \n\t"			/* add iteration number to j */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"subl   %%r12d,%%edi \n\t"				/* j - bound_val */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
		"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"movl   %%r8d,%%eax \n\t"				/* mov root1 to eax */ \
		"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
		"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
		"addl	%%edx,%%r8d \n\t"				/* increment root by prime */ \
		"cmpl   %%r13d,%%r8d \n\t"				/* root > interval? */ \
		"jb		1b \n\t"						/* repeat if necessary */

	// macro for iteratively adding an element (specified by r9d) to the end of a bucket list
	#define UPDATE_ROOT2_LOOP(it) \
		"1:		\n\t"	 						/* beginning of loop */ \
		"movl   %%r9d,%%ebx \n\t"				/* copy root to ebx */ \
		"shrl   $" BLOCKBITStxt ",%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
		"movl   %%ebx,%%ecx \n\t"				/* make into 64 bit register for addressing */ \
		"movl   %%r15d,%%edi \n\t"				/* copy j to edi */ \
		"addl	$" it ", %%edi \n\t"			/* add iteration number to j */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"subl   %%r12d,%%edi \n\t"				/* j - bound_val */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
		"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"movl   %%r9d,%%eax \n\t"				/* mov root1 to eax */ \
		"andl   $" BLOCKSIZEm1txt ",%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
		"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
		"addl	%%edx,%%r9d \n\t"				/* increment root by prime */ \
		"cmpl   %%r13d,%%r9d \n\t"				/* root > interval? */ \
		"jb		1b \n\t"						/* repeat if necessary */



#elif defined(GCC_ASM32X)

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

#endif

//this is in the poly library, even though the bulk of the time is spent
//bucketizing large primes, because it's where the roots of a poly are updated
void nextRoots(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//update the roots 
	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	int *rootupdates = dconf->rootupdates;

	update_t update_data = dconf->update_data;

	uint32 startprime = 2;
	uint32 bound = sconf->factor_base->B;

	char v = dconf->curr_poly->nu[dconf->numB];
	char sign = dconf->curr_poly->gray[dconf->numB];
	int *ptr;
	uint16 *sm_ptr;

	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32 med_B = sconf->factor_base->med_B;
	uint32 large_B = sconf->factor_base->large_B;

	uint32 j, interval; //, fb_offset;
	int k,numblocks;
	uint32 root1, root2, prime;

	int bound_index=0;
	int check_bound = BUCKET_ALLOC/2 - 1;
	uint32 bound_val = med_B;
	uint32 *numptr_p, *numptr_n, *sliceptr_p,*sliceptr_n;
	
#if !defined(USE_POLY_SSE2_ASM) || defined(PROFILING)
	uint32 *bptr;
	int bnum, room;
#endif

	uint8 logp=0;
	polysieve_t helperstruct;

	numblocks = sconf->num_blocks;
	interval = numblocks << BLOCKBITS;
	
	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->fb_bounds[0] = med_B;

		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		
		//reuse this for a sec...
		prime = 2*numblocks*lp_bucket_p->alloc_slices;

		//reset lp_buckets
		for (j=0;j<prime;j++)
			numptr_p[j] = 0;
	
		lp_bucket_p->num_slices = 0;

	}
	else
	{
		sliceptr_p = NULL;
		sliceptr_n = NULL;
		numptr_p = NULL;
		numptr_n = NULL;
	}

	k=0;
	ptr = &rootupdates[(v-1) * bound + startprime];	

	if (sign > 0)
	{
#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=startprime;j<sconf->sieve_small_fb_start;j++,ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			//we don't sieve these, so ordering doesn't matter
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			fb_p->root1[j] = (uint16)root1;
			fb_p->root2[j] = (uint16)root2;
			fb_n->root1[j] = (uint16)(prime - root2);
			fb_n->root2[j] = (uint16)(prime - root1);
			if (fb_n->root1[j] == prime)
				fb_n->root1[j] = 0;
			if (fb_n->root2[j] == prime)
				fb_n->root2[j] = 0;

		}

		// do one at a time up to the 10bit boundary, where
		// we can start doing things 8 at a time and be
		// sure we can use aligned moves (static_data_init).		
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint32)update_data.sm_firstroots1[j];
			root2 = (uint32)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
		
#if defined(GCC_ASM64X) || defined(_MSC_VER) //NOTDEF //GCC_ASM64X
		
		// update 8 at a time using SSE2 and no branching
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound];
		{
			// 8x root updates in parallel with SSE2, for use only with 32k
			// versions because of the signed comparisons.
			// stuff of potential use for full 16 bit updates
			// 8x unsigned gteq emulation
			// "psubusw	%%xmm1, %%xmm4 \n\t"	/* xmm2 := orig root - new root */
			// "pcmpeqw	%%xmm0, %%xmm4 \n\t"	/* xmm2 := a >= b ? 1 : 0 */
			// 8x MIN/MAX swap in SSE2 */
			// "movdqa	%%xmm2, %%xmm4 \n\t"			/* copy root2 */
			// "pcmpltd	%%xmm1, %%xmm4 \n\t"		/* root2 < root1? root2 --> 1's:root2 --> 0 */
			// "movdqa	%%xmm4, %%xmm5 \n\t"			/* copy ans */
			// "movdqa	%%xmm4, %%xmm6 \n\t"			/* copy ans */
			// "movdqa	%%xmm4, %%xmm7 \n\t"			/* copy ans */
			// "pand	%%xmm1, %%xmm4 \n\t"			/* copy root1 to where root1 is GT */
			// "pandn	%%xmm2, %%xmm5 \n\t"			/* copy root2 to where root2 is GT */
			// "pandn	%%xmm1, %%xmm6 \n\t"			/* copy root1 to where root1 is LT */
			// "pand	%%xmm2, %%xmm7 \n\t"			/* copy root2 to where root2 is LT */
			// "por	%%xmm4, %%xmm5 \n\t"			/* combine GT */
			// "por	%%xmm6, %%xmm7 \n\t"			/* combine LT */
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = sconf->factor_base->fb_10bit_B;		// 64
			h.stop = sconf->factor_base->fb_15bit_B;		// 68
			if ((h.stop - 8) > h.start)
				h.stop -= 8;

			COMPUTE_8X_SMALL_PROOTS;
			
			j = h.stop;
		}	

#else
		ptr = &dconf->rootupdates[(v-1) * bound + sconf->factor_base->fb_10bit_B];
		for (j=sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
#endif		

		// assembly code may not get all the way to 15 bits since we 
		// do things in blocks of 8 there.  Make sure we are at the 15 bit
		// boundary before we switch to using update_data.firstroots1/2.
		// this should only run a few iterations, if any.
		ptr = &dconf->rootupdates[(v-1) * bound + j];
		for ( ; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint16)update_data.sm_firstroots1[j];
			root2 = (uint16)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}	

		// continue one at a time once we exceed 15 bits, because the 8x SSE2
		// code has a hard time with unsigned 16 bit comparisons
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.firstroots1[j] = root2;
				update_data.firstroots2[j] = root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.firstroots1[j] = root1;
				update_data.firstroots2[j] = root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}	

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		

		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[med_B-1];

		if (med_B % 16 != 0)
		{
			printf("med_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((large_B - med_B) % 16 != 0)
		{
			printf("med range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];						//72
		helperstruct.large_B = med_B;				//80
		helperstruct.B = large_B;					//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
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
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"	/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm8 \n\t"			/* copy of prime for neg root loops */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "memory", "cc");

		// refresh local pointers and constants before entering the next loop
		numptr_n = helperstruct.numptr_n;
		numptr_p = helperstruct.numptr_p;
		sliceptr_n = helperstruct.sliceptr_n;
		sliceptr_p = helperstruct.sliceptr_p;
		ptr = helperstruct.ptr;	
		bound_val = helperstruct.bound_val;	
		check_bound = helperstruct.check_bound;
		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;


#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_PROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;
		}


#else
		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;
			//fb_offset = (j - bound_val) << 16;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif
			
		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[large_B-1];

		if (large_B % 16 != 0)
		{
			printf("large_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((bound - large_B) % 16 != 0)
		{
			printf("large range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = large_B;				//80
		helperstruct.B = bound;						//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
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
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"psubd	%%xmm3, %%xmm1 \n\t"			/* root1 -= ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"psubd	%%xmm3, %%xmm2 \n\t"			/* root2 -= ptr */ \
			"pxor	%%xmm4, %%xmm4 \n\t"			/* zero xmm4 */ \
			"pxor	%%xmm5, %%xmm5 \n\t"			/* zero xmm5 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"pcmpgtd	%%xmm1, %%xmm4 \n\t"		/* signed comparison: 0 > root1? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm2, %%xmm5 \n\t"		/* signed comparison: 0 > root2? if so, set xmm5 dword to 1's */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy prime to compute neg roots */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy prime to compute neg roots */ \
			"paddd	%%xmm4, %%xmm1 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"paddd	%%xmm5, %%xmm2 \n\t"			/* selectively add back prime (modular subtract) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* edx is free at this point... can it be used? */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");

		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;


#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_PROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);
			
			j++;
		}


#else
		logp = update_data.logp[j-1];
		for (j=large_B;j<bound;j++,ptr++)				
		{				
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_P;		

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
			FILL_ONE_PRIME_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

	}
	else
	{

#ifdef QS_TIMING
		gettimeofday(&qs_timing_start, NULL);
#endif

		for (j=startprime;j<sconf->sieve_small_fb_start;j++,ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			//we don't sieve these, so ordering doesn't matter
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			fb_p->root1[j] = (uint16)root1;
			fb_p->root2[j] = (uint16)root2;
			fb_n->root1[j] = (uint16)(prime - root2);
			fb_n->root2[j] = (uint16)(prime - root1);
			if (fb_n->root1[j] == prime)
				fb_n->root1[j] = 0;
			if (fb_n->root2[j] == prime)
				fb_n->root2[j] = 0;

		}

		// do one at a time up to the 10bit boundary, where
		// we can start doing things 8 at a time and be
		// sure we can use aligned moves (static_data_init).	
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint32)update_data.sm_firstroots1[j];
			root2 = (uint32)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
		
#if defined(GCC_ASM64X) || defined(_MSC_VER) //NOTDEF //GCC_ASM64X
		// update 8 at a time using SSE2 and no branching		
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound];
		{
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = sconf->factor_base->fb_10bit_B;		// 64
			h.stop = sconf->factor_base->fb_15bit_B;		// 68
			if ((h.stop - 8) > h.start)
				h.stop -= 8;

			COMPUTE_8X_SMALL_NROOTS;
			
			j = h.stop;
		}
		sm_ptr = &dconf->sm_rootupdates[(v-1) * bound + j];
		

#else
		ptr = &dconf->rootupdates[(v-1) * bound + sconf->factor_base->fb_10bit_B];
		for (j=sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}
#endif		

		// assembly code may not get all the way to 15 bits since we 
		// do things in blocks of 8 there.  Make sure we are at the 15 bit
		// boundary before we switch to using update_data.firstroots1/2.
		// this should only run a few iterations, if any.
		ptr = &dconf->rootupdates[(v-1) * bound + j];
		for ( ; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint16)update_data.sm_firstroots1[j];
			root2 = (uint16)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16)root2;
				update_data.sm_firstroots2[j] = (uint16)root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16)root1;
				update_data.sm_firstroots2[j] = (uint16)root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}	

		// continue one at a time once we exceed 15 bits, because the 8x SSE2
		// code has a hard time with unsigned 16 bit comparisons
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.firstroots1[j] = root2;
				update_data.firstroots2[j] = root1;

				fb_p->root1[j] = (uint16)root2;
				fb_p->root2[j] = (uint16)root1;
				fb_n->root1[j] = (uint16)(prime - root1);
				fb_n->root2[j] = (uint16)(prime - root2);
			}
			else
			{
				update_data.firstroots1[j] = root1;
				update_data.firstroots2[j] = root2;

				fb_p->root1[j] = (uint16)root1;
				fb_p->root2[j] = (uint16)root2;
				fb_n->root1[j] = (uint16)(prime - root2);
				fb_n->root2[j] = (uint16)(prime - root1);
			}
		}	

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		
		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[med_B-1];

		if (med_B % 16 != 0)
		{
			printf("med_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((large_B - med_B) % 16 != 0)
		{
			printf("med range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = med_B;				//80
		helperstruct.B = large_B;					//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
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
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t"	/* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy of prime for neg root calculation */ \
			"movdqa %%xmm0, %%xmm8 \n\t"			/* copy of prime for neg root loops */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm0 \n\t" 				/* next prime */ \
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
			"movd	%%xmm0,%%edx \n\t"				/* else, extract prime from xmm0 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"psrldq $4,%%xmm8 \n\t" 				/* next prime */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
			"movd	%%xmm8,%%edx \n\t"				/* else, extract prime from xmm8 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2_LOOP("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "memory", "cc");

		// refresh local pointers and constants before entering the next loop
		numptr_n = helperstruct.numptr_n;
		numptr_p = helperstruct.numptr_p;
		sliceptr_n = helperstruct.sliceptr_n;
		sliceptr_p = helperstruct.sliceptr_p;
		ptr = helperstruct.ptr;	
		bound_val = helperstruct.bound_val;	
		check_bound = helperstruct.check_bound;
		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_NROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);

			j++;
		}


#else

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			FILL_ONE_PRIME_LOOP_P(j);

			root1 = (prime - update_data.firstroots1[j]);
			root2 = (prime - update_data.firstroots2[j]);

			FILL_ONE_PRIME_LOOP_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG3 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif

		
#if defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		logp = update_data.logp[large_B-1];
		
		if (large_B % 16 != 0)
		{
			printf("large_B must be divisible by 16!\n");
			exit(-1);
		}
		if ((bound - large_B) % 16 != 0)
		{
			printf("large range must be divisible by 16!\n");
			exit(-1);
		}

		//load up our helper struct, so we don't have a list
		//of a 100 things to stick into our asm block
		helperstruct.numptr_n = numptr_n;			//0
		helperstruct.numptr_p = numptr_p;			//8
		helperstruct.sliceptr_n = sliceptr_n;		//16
		helperstruct.sliceptr_p = sliceptr_p;		//24
		helperstruct.update_data_prime = update_data.prime;	//32
		helperstruct.update_data_root1 = update_data.firstroots1;	//40
		helperstruct.update_data_root2 = update_data.firstroots2;	//48
		helperstruct.update_data_logp = update_data.logp;	//56
		helperstruct.lp_bucket_p = lp_bucket_p;		//64
		helperstruct.ptr = &rootupdates[(v-1) * bound];  //72;						//72
		helperstruct.large_B = large_B;				//80
		helperstruct.B = bound;						//84
		helperstruct.interval = interval;			//88
		helperstruct.numblocks = numblocks;			//92
		helperstruct.bound_val = bound_val;			//96
		helperstruct.bound_index = bound_index;		//100
		helperstruct.check_bound = check_bound;		//104
		helperstruct.logp = logp;					//108

		ASM_G (		\
			"movq	%0,%%rsi \n\t"					/* move helperstruct into rsi */ \
			"movl   80(%%rsi,1),%%r15d \n\t"		/* large_B = j = r15d */ \
													/* do the loop comparison */ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j >= bound ? */ \
			"jae    9f	\n\t"						/* jump to end of loop, if test fails */ \
			"8: \n\t"	\
				/* ================================================ */	\
				/* ========== BEGIN CHECK_NEW_SLICE BLOCK ========= */	\
				/* ================================================ */	\
			"cmpl   104(%%rsi,1),%%r15d	\n\t"		/* compare j with check_bound */ \
				/* note this is the counter j, not the byte offset j */ \
			"jge     1f \n\t"						/* jump into "if" code if comparison works */ \
				/* else, this is the "else-if" check */ \
			"movl   %%r15d,%%ebx \n\t"				/* copy j into ebx */ \
			"subl   96(%%rsi,1),%%ebx \n\t"			/* ebx = j - bound_val */ \
			"cmpl   $0xffff,%%ebx \n\t"				/* compare to 2^16 */ \
			"jbe    2f \n\t"						/* exit CHECK_NEW_SLICE if this comparison fails too */ \
				/* now we are in the else-if block of CHECK_NEW_SLICE */ \
			"xorq	%%rdx, %%rdx \n\t"				/* clear rdx */ \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t"				/* clear rbx */ \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (translate to bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 (negative blocks are contiguous) */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shlq   $" BUCKET_BITStxt ",%%rbx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"cmp	%%rax,%%rax \n\t"				/* force jump */
			"je		2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"1:		\n\t"									\
				/* now we are in the if block of CHECK_NEW_SLICE */ \
			"xorl   %%ecx,%%ecx \n\t"				/* ecx = room  = 0 */ \
			"xorq	%%rbx, %%rbx \n\t"				/* loop counter = 0 */ \
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare with numblocks */ \
			"jae    3f \n\t"						/* jump past loop if condition met */ \
				/* condition not met, put a couple things in registers */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	0(%%rsi,1),%%r11 \n\t"			/* numptr_n into r11 */ \
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
			"cmpl   92(%%rsi,1),%%ebx \n\t"			/* compare to numblocks */ \
			"jl     5b \n\t"						/* iterate loop if condition met */ \
			"3:		\n\t"							\
			"movl   $" BUCKET_ALLOCtxt ",%%ebx \n\t" /* move bucket allocation into register for subtraction */ \
			"subl   %%ecx,%%ebx \n\t"				/* room = bucket_alloc - room */ \
			"cmpl   $31,%%ebx \n\t"					/* answer less than 32? */ \
			"movl   %%ebx,%%ecx \n\t"				/* copy answer back to room register */ \
			"jle    4f \n\t"						/* jump if less than */ \
			"sarl   %%ebx	\n\t"					/* room >> 1 (copy of room) */ \
			"addl   %%ebx,104(%%rsi,1) \n\t"		/* add (room >> 1) to check_bound */ \
			"cmpq	%%rax,%%rax \n\t"				/* force jump */
			"je     2f \n\t"						/* jump out of CHECK_NEW_SLICE */ \
			"4:		\n\t"							\
				/* now we are inside the (room < 2) block */ \
			"xorq	%%rax, %%rax \n\t" \
			"movl   %%r15d,%%eax \n\t"				/* copy j to scratch reg */ \
			"shll   $0x4,%%eax \n\t"				/* multiply by 16 bytes per j */ \
			"xorq	%%rdx, %%rdx \n\t" \
			"movl   100(%%rsi,1),%%edx \n\t"		/* move bound_index into rdx */ \
			"movq   64(%%rsi,1),%%r9 \n\t"			/* move lp_bucket_p ptr into r9 */ \
			"movq	16(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->logp ptr into r8 */ \
			"movq	56(%%rsi,1),%%r14 \n\t"			/* move updata_data.logp pointer into r14 */ \
			"movzbl (%%r14,%%r15,1),%%ebx \n\t"		/* bring in logp */ \
			"movb	%%bl, 108(%%rsi,1) \n\t"		/* shove logp into output */ \
			"movb   %%bl,(%%r8,%%rdx,1) \n\t"		/* mov logp into lp_bucket_p->logp[bound_index] */ \
			"incq   %%rdx \n\t"						/* increment bound_index locally */ \
			"movl   %%edx,100(%%rsi,1) \n\t"		/* copy bound_index back to structure */ \
			"movq	8(%%r9,1),%%r8 \n\t"			/* move lp_bucket_p->fb_bounds ptr into r8 */ \
			"movl   %%r15d,(%%r8,%%rdx,4) \n\t"		/* mov j into lp_bucket_p->fb_bounds[bound_index] */ \
				/* note this is the counter j, not the byte offset j */ \
			"movl   %%r15d,96(%%rsi,1) \n\t"		/* bound_val = j */ \
			"xorq	%%rbx, %%rbx \n\t" \
			"movl   92(%%rsi,1),%%ebx \n\t"			/* put numblocks into ebx */ \
			"shll	$2,%%ebx \n\t"					/* numblocks * 4 (bytes) */ \
			"shll	$1,%%ebx \n\t"					/* numblocks << 1 */ \
			"addq   %%rbx,8(%%rsi,1) \n\t"			/* numptr_p += (numblocks << 1) */ \
			"addq   %%rbx,0(%%rsi,1) \n\t"			/* numptr_n += (numblocks << 1) */ \
			"shll   $" BUCKET_BITStxt ",%%ebx \n\t"	/* numblocks << (BUCKET_BITS + 1) */ \
				/* note also, this works because we've already left shifted by 1 */ \
			"addq   %%rbx,24(%%rsi,1) \n\t"			/* sliceptr_p += (numblocks << 11) */ \
			"addq   %%rbx,16(%%rsi,1) \n\t"			/* sliceptr_n += (numblocks << 11) */ \
			"addl   $" HALFBUCKET_ALLOCtxt ",104(%%rsi,1) \n\t"		/* add 2^(BUCKET_BITS-1) to check_bound */ \
			"2:		\n\t"						\
				/* ================================================ */	\
				/* ============ BEGIN - GET NEW ROOTS 1 =========== */	\
				/* ================================================ */	\
			"movq	72(%%rsi,1),%%rdi \n\t"			/* edi = ptr */ \
			"movq	40(%%rsi,1),%%r14 \n\t"			/* move updata_data.root1 pointer into r14 */ \
			"movdqa (%%rdi,%%r15,4), %%xmm3 \n\t"	/* xmm3 = next 4 values of rootupdates */ \
			"movq	48(%%rsi,1),%%r13 \n\t"			/* move updata_data.root2 pointer into r13 */ \
			"movdqa (%%r14,%%r15,4), %%xmm1 \n\t"	/* xmm1 = next 4 values of root1 */ \
			"paddd	%%xmm3, %%xmm1 \n\t"			/* root1 += ptr */ \
			"movdqa (%%r13,%%r15,4), %%xmm2 \n\t"	/* xmm2 = next 4 values of root2 */ \
			"movq	32(%%rsi,1),%%r12 \n\t"			/* move updata_data.prime pointer into r12 */ \
			"paddd	%%xmm3, %%xmm2 \n\t"			/* root2 += ptr */ \
			"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1 to xmm4 */ \
			"movdqa (%%r12,%%r15,4), %%xmm0 \n\t"	/* xmm0 = next 4 primes */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2 to xmm5 */ \
			"pcmpgtd	%%xmm0, %%xmm4 \n\t"		/* signed comparison: root1 > p? if so, set xmm4 dword to 1's */ \
			"pcmpgtd	%%xmm0, %%xmm5 \n\t"		/* signed comparison: root2 > p? if so, set xmm5 dword to 1's */ \
			"movdqa %%xmm0, %%xmm6 \n\t"			/* copy prime to xmm6 */ \
			"pand	%%xmm0, %%xmm4 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"movdqa %%xmm0, %%xmm7 \n\t"			/* copy prime to xmm7 */ \
			"pand	%%xmm0, %%xmm5 \n\t"			/* copy prime to overflow locations (are set to 1) */ \
			"psubd	%%xmm4, %%xmm1 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm1, (%%r14,%%r15,4) \n\t"	/* save new root1 values */ \
			"psubd	%%xmm5, %%xmm2 \n\t"			/* selectively sub back prime (modular addition) */ \
			"movdqa %%xmm2, (%%r13,%%r15,4) \n\t"	/* save new root2 values */ \
			"psubd	%%xmm1, %%xmm6 \n\t"			/* form negative root1's; prime - root1 */ \
			"psubd	%%xmm2, %%xmm7 \n\t"			/* form negative root2's; prime - root2 */ \
			"movq	8(%%rsi,1),%%r10 \n\t"			/* numptr_p into r10 */ \
			"movq	24(%%rsi,1),%%r11 \n\t"			/* sliceptr_p into r11 */ \
			"movl	88(%%rsi,1),%%r13d \n\t"		/* interval */ \
			"movl	96(%%rsi,1),%%r12d \n\t"		/* store bound_val in a register */ \
				/* edx is free at this point... can it be used? */ \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,1 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,1 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,2 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,2 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,3 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,3 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm1 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm2 \n\t" 				/* next root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm1,%%r8d \n\t"				/* else, extract root1,4 from xmm1 */ \
			"movd	%%xmm2,%%r9d \n\t"				/* else, extract root2,4 from xmm2 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE POS BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,1 ========== */	\
				/* ================================================ */	\
			"movq	0(%%rsi,1),%%r10 \n\t"			/* numptr_n into r10 */ \
			"movq	16(%%rsi,1),%%r11 \n\t"			/* sliceptr_n into r10 */ \
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,1 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,1 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,1 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    1f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("0") \
			"1:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,2 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,2 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,2 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,2 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    3f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("1") \
			"3:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,3 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,3 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,3 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"psrldq $4,%%xmm6 \n\t" 				/* nextn root1 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,3 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"psrldq $4,%%xmm7 \n\t" 				/* nextn nroot2 */ \
			"jae    5f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("2") \
			"5:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT1,4 ========== */	\
				/* ================================================ */	\
			"movd	%%xmm6,%%r8d \n\t"				/* else, extract nroot1,4 from xmm6 */ \
			"movd	%%xmm7,%%r9d \n\t"				/* else, extract nroot2,4 from xmm7 */ \
			"cmpl	%%r13d,%%r8d \n\t"				/* root1 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT1("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* =========== UPDATE NEG BUCKET - ROOT2,4 ========== */	\
				/* ================================================ */	\
			"cmpl	%%r13d,%%r9d \n\t"				/* root2 > interval? */ \
			"jae    7f \n\t" 						/* jump if CF = 1 */ \
			UPDATE_ROOT2("3") \
			"7:		\n\t" \
				/* ================================================ */	\
				/* ======== END OF LOOP - UPDATE AND CHECK ======== */	\
				/* ================================================ */	\
			"addl   $4,%%r15d \n\t"					/* increment j by 1*/ \
			"cmpl   84(%%rsi,1),%%r15d \n\t"		/* j < bound ? */ \
			"jb     8b \n\t"	\
			"9:		\n\t"				\
			"movl	%%r15d, %%eax \n\t" \
			:  \
			: "g"(&helperstruct) \
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "memory", "cc");

		bound_index = helperstruct.bound_index;
		logp = helperstruct.logp;

#elif defined(HAS_SSE2)

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound; )
		{
			CHECK_NEW_SLICE(j);

			COMPUTE_4_NROOTS(j);

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);

			j++;

			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];
			prime = update_data.prime[j];

			FILL_ONE_PRIME_P(j);

			root1 = (prime - root1);
			root2 = (prime - root2);

			FILL_ONE_PRIME_N(j);
			
			j++;
		}


#else

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound;j++,ptr++)				
		{				
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_N;		

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
			FILL_ONE_PRIME_N(j);
		}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		POLY_STG4 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

	}

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->num_slices = bound_index + 1;
		lp_bucket_p->logp[bound_index] = logp;
	}

	return;
}
