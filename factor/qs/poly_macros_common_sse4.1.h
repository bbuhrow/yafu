#include "common.h"

#if defined(GCC_ASM64X)

#if defined( _WIN32)
#define ASM_ ASM_M
#else
#define ASM_ ASM_G
#endif

	#define COMPUTE_8X_SMALL_PROOTS_SSE41	\
		ASM_ (	\
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
			"pxor	%%xmm6, %%xmm6 \n\t" \
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the P side */	\
			"movdqa	(%%r8, %%rax, 2), %%xmm3 \n\t"			/* xmm3 = ptr */	\
			"movdqa (%%rbx, %%rax, 2), %%xmm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"movdqa (%%rcx, %%rax, 2), %%xmm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"movdqa	%%xmm1, %%xmm4 \n\t"					/* copy r1 */ \
			"movdqa	%%xmm2, %%xmm5 \n\t"					/* copy r2 */ \
			"movdqa	%%xmm1, %%xmm7 \n\t"					/* copy r1 */ \
			"movdqa	%%xmm2, %%xmm8 \n\t"					/* copy r2 */ \
			"psubw	%%xmm3, %%xmm1 \n\t"					/* root1 -= ptr */	\
			"psubw	%%xmm3, %%xmm2 \n\t"					/* root2 -= ptr */	\
			"psubusw	%%xmm3, %%xmm4 \n\t"					/* root1 -= ptr */	\
			"psubusw	%%xmm3, %%xmm5 \n\t"					/* root2 -= ptr */	\
			"movdqa (%%rdx, %%rax, 2), %%xmm0 \n\t"			/* xmm0 = next 8 primes */	\
			"pcmpeqw	%%xmm3, %%xmm7 \n\t"				/* xmm4 := root1 == ptr ? 1 : 0 */ \
			"pcmpeqw	%%xmm3, %%xmm8 \n\t"				/* xmm5 := root2 == ptr ? 1 : 0 */ \
			"pcmpeqw	%%xmm6, %%xmm4 \n\t"				/* xmm4 := ptr >= root1 ? 1 : 0 */ \
			"pcmpeqw	%%xmm6, %%xmm5 \n\t"				/* xmm5 := ptr >= root2 ? 1 : 0 */ \
			"pandn	%%xmm4, %%xmm7 \n\t"					/* ptr > root1 (greater than and not equal) */	\
			"pandn	%%xmm5, %%xmm8 \n\t"					/* ptr > root1 (greater than and not equal) */	\
			"pand	%%xmm0, %%xmm7 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"pand	%%xmm0, %%xmm8 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"paddw	%%xmm7, %%xmm1 \n\t"					/* selectively add back prime (modular sub) */	\
			"paddw	%%xmm8, %%xmm2 \n\t"					/* selectively add back prime (modular sub) */	\
			"movdqa %%xmm2, %%xmm5 \n\t"					/* xmm5 = root2 copy */	\
			"pmaxuw	%%xmm1, %%xmm5 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"pminuw	%%xmm1, %%xmm2 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
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
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "rax", "rsi", "rbx", "rcx", "rdx",	\
			"r8", "r9", "r10", "r11", "r12", "r15", "cc", "memory");

	#define COMPUTE_8X_SMALL_NROOTS_SSE41	\
		ASM_ (	\
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
			"pxor	%%xmm6, %%xmm6 \n\t" \
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the N side */	\
			"movdqa (%%r8, %%rax, 2), %%xmm3 \n\t"			/* xmm3 = next 8 updates */	\
			"movdqa (%%rbx, %%rax, 2), %%xmm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"movdqa (%%rcx, %%rax, 2), %%xmm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"movdqa	%%xmm1, %%xmm7 \n\t"					/* copy root1 */ \
			"movdqa	%%xmm2, %%xmm8 \n\t"					/* copy root2 */ \
			"paddw	%%xmm3, %%xmm1 \n\t"					/* root1 += ptr */	\
			"paddw	%%xmm3, %%xmm2 \n\t"					/* root2 += ptr */	\
			"paddusw	%%xmm3, %%xmm7 \n\t"					/* root1 += ptr */	\
			"paddusw	%%xmm3, %%xmm8 \n\t"					/* root2 += ptr */	\
			"movdqa (%%rdx, %%rax, 2), %%xmm0 \n\t"			/* xmm0 = next 8 primes */	\
			"movdqa	%%xmm0, %%xmm4 \n\t"					/* copy primes before comparison */ \
			"movdqa	%%xmm0, %%xmm5 \n\t"					/* copy primes before comparison */ \
			"psubusw	%%xmm7, %%xmm4 \n\t"				/* xmm4 := prime - root1 */ \
			"psubusw	%%xmm8, %%xmm5 \n\t"				/* xmm5 := prime - root2 */ \
			"pcmpeqw	%%xmm6, %%xmm4 \n\t"				/* xmm4 := root1 >= prime ? 1 : 0 */ \
			"pcmpeqw	%%xmm6, %%xmm5 \n\t"				/* xmm5 := root2 >= prime ? 1 : 0 */ \
			"pand	%%xmm0, %%xmm4 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"pand	%%xmm0, %%xmm5 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"psubw	%%xmm4, %%xmm1 \n\t"					/* selectively sub back prime (modular add) */	\
			"psubw	%%xmm5, %%xmm2 \n\t"					/* selectively sub back prime (modular add) */	\
			"movdqa %%xmm2, %%xmm5 \n\t"					/* xmm5 = root2 copy */	\
			"pmaxuw	%%xmm1, %%xmm5 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"pminuw	%%xmm1, %%xmm2 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
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
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "rax", "rsi", "rbx", "rcx", "rdx",	\
			"r8", "r9", "r10", "r11", "r12", "r15", "cc", "memory");

#elif defined (_MSC_VER) && defined(_WIN64)

#define COMPUTE_8X_SMALL_PROOTS_SSE41	\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i tmp3;	\
				__m128i tmp4;	\
				__m128i zeros; \
				zeros = _mm_xor_si128(zeros, zeros); \
				for (j = h.start; j < h.stop; j += 8) \
				{	\
					ptrs = _mm_load_si128((__m128i *)(h.updates + j)); \
					root1s = _mm_load_si128((__m128i *)(h.first_r1 + j)); \
					root2s = _mm_load_si128((__m128i *)(h.first_r2 + j)); \
					primes = _mm_load_si128((__m128i *)(h.primes + j)); \
					tmp1 = root1s; \
					tmp2 = root2s; \
					tmp3 = root1s; \
					tmp4 = root2s; \
					root1s = _mm_sub_epi16(root1s, ptrs); \
					root2s = _mm_sub_epi16(root2s, ptrs); \
					tmp1 = _mm_subs_epu16(tmp1, ptrs); \
					tmp2 = _mm_subs_epu16(tmp2, ptrs); \
					tmp3 = _mm_cmpeq_epi16(tmp3,ptrs);	\
					tmp4 = _mm_cmpeq_epi16(tmp4,ptrs);	\
					tmp1 = _mm_cmpeq_epi16(tmp1,zeros);	\
					tmp2 = _mm_cmpeq_epi16(tmp2,zeros);	\
					tmp3 = _mm_andnot_si128(tmp3, tmp1);	\
					tmp4 = _mm_andnot_si128(tmp4, tmp2);	\
					tmp3 = _mm_and_si128(tmp3, primes);	\
					tmp4 = _mm_and_si128(tmp4, primes);	\
					root1s = _mm_add_epi16(root1s, tmp3);	\
					root2s = _mm_add_epi16(root2s, tmp4);	\
					tmp1 = root2s;	\
					tmp1 = _mm_max_epu16(tmp1, root1s); \
					root2s = _mm_min_epu16(root2s, root1s); \
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

	#define COMPUTE_8X_SMALL_NROOTS_SSE41	\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i tmp3;	\
				__m128i tmp4;	\
				__m128i zeros; \
				zeros = _mm_xor_si128(zeros, zeros); \
				for (j = h.start; j < h.stop; j += 8) \
				{	\
					ptrs = _mm_load_si128((__m128i *)(h.updates + j)); \
					root1s = _mm_load_si128((__m128i *)(h.first_r1 + j)); \
					root2s = _mm_load_si128((__m128i *)(h.first_r2 + j)); \
					primes = _mm_load_si128((__m128i *)(h.primes + j)); \
					tmp3 = root1s; \
					tmp4 = root2s; \
					root1s = _mm_add_epi16(root1s, ptrs); \
					root2s = _mm_add_epi16(root2s, ptrs); \
					tmp3 = _mm_adds_epu16(tmp3, ptrs); \
					tmp4 = _mm_adds_epu16(tmp4, ptrs); \
					tmp1 = primes; \
					tmp2 = primes; \
					tmp1 = _mm_subs_epu16(tmp1, tmp3); \
					tmp2 = _mm_subs_epu16(tmp2, tmp4); \
					tmp1 = _mm_cmpeq_epi16(tmp1,zeros);	\
					tmp2 = _mm_cmpeq_epi16(tmp2,zeros);	\
					tmp1 = _mm_and_si128(tmp1, primes);	\
					tmp2 = _mm_and_si128(tmp2, primes);	\
					root1s = _mm_sub_epi16(root1s, tmp1);	\
					root2s = _mm_sub_epi16(root2s, tmp2);	\
					tmp1 = root2s;	\
					tmp1 = _mm_max_epu16(tmp1, root1s); \
					root2s = _mm_min_epu16(root2s, root1s); \
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

