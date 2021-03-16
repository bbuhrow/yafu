
#if defined(GCC_ASM64X)


#define COMPUTE_16X_SMALL_PROOTS_AVX2	\
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
			"vpxor	%%ymm6, %%ymm6, %%ymm6 \n\t" \
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the P side */	\
			"vmovdqa	(%%r8, %%rax, 2), %%ymm3 \n\t"			/* xmm3 = ptr */	\
			"vmovdqa    (%%rbx, %%rax, 2), %%ymm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"vmovdqa    (%%rcx, %%rax, 2), %%ymm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"vmovdqa	%%ymm1, %%ymm4 \n\t"					/* copy r1 */ \
			"vmovdqa	%%ymm2, %%ymm5 \n\t"					/* copy r2 */ \
			"vmovdqa	%%ymm1, %%ymm7 \n\t"					/* copy r1 */ \
			"vmovdqa	%%ymm2, %%ymm8 \n\t"					/* copy r2 */ \
			"vpsubw	    %%ymm3, %%ymm1, %%ymm1 \n\t"					/* root1 -= ptr */	\
			"vpsubw	    %%ymm3, %%ymm2, %%ymm2 \n\t"					/* root2 -= ptr */	\
			"vpsubusw	%%ymm3, %%ymm4, %%ymm4 \n\t"					/* root1 -= ptr */	\
			"vpsubusw	%%ymm3, %%ymm5, %%ymm5 \n\t"					/* root2 -= ptr */	\
			"vmovdqa    (%%rdx, %%rax, 2), %%ymm0 \n\t"			/* xmm0 = next 8 primes */	\
			"vpcmpeqw	%%ymm3, %%ymm7, %%ymm7 \n\t"				/* xmm4 := root1 == ptr ? 1 : 0 */ \
			"vpcmpeqw	%%ymm3, %%ymm8, %%ymm8 \n\t"				/* xmm5 := root2 == ptr ? 1 : 0 */ \
			"vpcmpeqw	%%ymm6, %%ymm4, %%ymm4 \n\t"				/* xmm4 := ptr >= root1 ? 1 : 0 */ \
			"vpcmpeqw	%%ymm6, %%ymm5, %%ymm5 \n\t"				/* xmm5 := ptr >= root2 ? 1 : 0 */ \
			"vpandn	    %%ymm4, %%ymm7, %%ymm7 \n\t"					/* ptr > root1 (greater than and not equal) */	\
			"vpandn	    %%ymm5, %%ymm8, %%ymm8 \n\t"					/* ptr > root1 (greater than and not equal) */	\
			"vpand	    %%ymm0, %%ymm7, %%ymm7 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpand	    %%ymm0, %%ymm8, %%ymm8 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpaddw	    %%ymm7, %%ymm1, %%ymm1 \n\t"					/* selectively add back prime (modular sub) */	\
			"vpaddw	    %%ymm8, %%ymm2, %%ymm2 \n\t"					/* selectively add back prime (modular sub) */	\
			"vmovdqa    %%ymm2, %%ymm5 \n\t"					/* xmm5 = root2 copy */	\
			"vpmaxuw	%%ymm1, %%ymm5, %%ymm5 \n\t"					/* xmm5 = root2 > root1 ? root2 : root1 */	\
			"vpminuw	%%ymm1, %%ymm2, %%ymm2 \n\t"					/* xmm2 = root2 < root1 ? root2 : root1 */	\
															/* now copy results to appropriate data structures */	\
			"vmovdqa	%%ymm0, %%ymm4 \n\t"					/* copy primes */	\
															/* root1p always gets the smaller roots (LT) */	\
			"vmovdqa	%%ymm2, (%%r9, %%rax, 2) \n\t"			/* update root1p */	\
			"vpsubw	    %%ymm2, %%ymm0, %%ymm0 \n\t"					/* prime - LT roots */	\
			"vmovdqa	%%ymm2, (%%rbx, %%rax, 2) \n\t"			/* update firstroots1 */	\
															/* root2p always gets the bigger roots (GT) */	\
			"vmovdqa	%%ymm5, (%%r10, %%rax, 2) \n\t"			/* update root2p */	\
			"vpsubw	    %%ymm5, %%ymm4, %%ymm4 \n\t"					/* prime - GT roots */	\
			"vmovdqa	%%ymm5, (%%rcx, %%rax, 2) \n\t"			/* update firstroots2 */	\
															/* root1n always gets prime - bigger roots (LT) */	\
			"vmovdqa	%%ymm4, (%%r11, %%rax, 2) \n\t"			/* update root1n */	\
															/* root2n always gets prime - smaller roots (GT) */	\
			"vmovdqa	%%ymm0, (%%r12, %%rax, 2) \n\t"			/* update root2n */	\
			"addl	$16, %%eax \n\t"	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "g"(&h)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "rax", "rsi", "rbx", "rcx", "rdx",	\
			"r8", "r9", "r10", "r11", "r12", "r15", "cc", "memory");

#define COMPUTE_16X_SMALL_NROOTS_AVX2	\
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
			"vpxor	%%ymm6, %%ymm6, %%ymm6 \n\t" \
			"cmpl	%%r15d, %%eax \n\t"	\
			"jge	1f \n\t"	\
			"0: \n\t"	\
			/* compute 8 new roots on the N side */	\
			"vmovdqa    (%%r8, %%rax, 2), %%ymm3 \n\t"			/* xmm3 = next 8 updates */	\
			"vmovdqa    (%%rbx, %%rax, 2), %%ymm1 \n\t"			/* xmm1 = next 8 values of root1 */	\
			"vmovdqa    (%%rcx, %%rax, 2), %%ymm2 \n\t"			/* xmm2 = next 8 values of root2 */	\
			"vmovdqa	%%ymm1, %%ymm7 \n\t"					/* copy root1 */ \
			"vmovdqa	%%ymm2, %%ymm8 \n\t"					/* copy root2 */ \
			"vpaddw	    %%ymm3, %%ymm1, %%ymm1 \n\t"					/* root1 += ptr */	\
			"vpaddw	    %%ymm3, %%ymm2, %%ymm2 \n\t"					/* root2 += ptr */	\
			"vpaddusw	%%ymm3, %%ymm7, %%ymm7 \n\t"					/* root1 += ptr */	\
			"vpaddusw	%%ymm3, %%ymm8, %%ymm8 \n\t"					/* root2 += ptr */	\
			"vmovdqa    (%%rdx, %%rax, 2), %%ymm0 \n\t"			/* ymm0 = next 8 primes */	\
			"vmovdqa	%%ymm0, %%ymm4 \n\t"					/* copy primes before comparison */ \
			"vmovdqa	%%ymm0, %%ymm5 \n\t"					/* copy primes before comparison */ \
			"vpsubusw	%%ymm7, %%ymm4, %%ymm4 \n\t"				/* ymm4 := prime - root1 */ \
			"vpsubusw	%%ymm8, %%ymm5, %%ymm5 \n\t"				/* ymm5 := prime - root2 */ \
			"vpcmpeqw	%%ymm6, %%ymm4, %%ymm4 \n\t"				/* ymm4 := root1 >= prime ? 1 : 0 */ \
			"vpcmpeqw	%%ymm6, %%ymm5, %%ymm5 \n\t"				/* ymm5 := root2 >= prime ? 1 : 0 */ \
			"vpand	    %%ymm0, %%ymm4, %%ymm4 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpand	    %%ymm0, %%ymm5, %%ymm5 \n\t"					/* copy prime to overflow locations (are set to 1) */	\
			"vpsubw	    %%ymm4, %%ymm1, %%ymm1 \n\t"					/* selectively sub back prime (modular add) */	\
			"vpsubw	    %%ymm5, %%ymm2, %%ymm2 \n\t"					/* selectively sub back prime (modular add) */	\
			"vmovdqa    %%ymm2, %%ymm5 \n\t"					/* ymm5 = root2 copy */	\
			"vpmaxuw	%%ymm1, %%ymm5, %%ymm5 \n\t"					/* ymm5 = root2 > root1 ? root2 : root1 */	\
			"vpminuw	%%ymm1, %%ymm2, %%ymm2 \n\t"					/* ymm2 = root2 < root1 ? root2 : root1 */	\
															/* now copy results to appropriate data structures */	\
			"vmovdqa	%%ymm0, %%ymm4 \n\t"					/* copy primes */	\
															/* root1p always gets the smaller roots (LT) */	\
			"vmovdqa	%%ymm2, (%%r9, %%rax, 2) \n\t"			/* update root1p */	\
			"vpsubw	    %%ymm2, %%ymm0, %%ymm0 \n\t"					/* prime - LT roots */	\
			"vmovdqa	%%ymm2, (%%rbx, %%rax, 2) \n\t"			/* update firstroots1 */	\
															/* root2p always gets the bigger roots (GT) */	\
			"vmovdqa	%%ymm5, (%%r10, %%rax, 2) \n\t"			/* update root2p */	\
			"vpsubw	    %%ymm5, %%ymm4, %%ymm4 \n\t"					/* prime - GT roots */	\
			"vmovdqa	%%ymm5, (%%rcx, %%rax, 2) \n\t"			/* update firstroots2 */	\
															/* root1n always gets prime - bigger roots (LT) */	\
			"vmovdqa	%%ymm4, (%%r11, %%rax, 2) \n\t"			/* update root1n */	\
															/* root2n always gets prime - smaller roots (GT) */	\
			"vmovdqa	%%ymm0, (%%r12, %%rax, 2) \n\t"			/* update root2n */	\
			"addl	$16, %%eax \n\t"	\
			"cmpl	%%r15d, %%eax \n\t"	\
			"jb		0b \n\t"	\
			"1: \n\t"	\
			:	\
			: "g"(&h)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "rax", "rsi", "rbx", "rcx", "rdx",	\
			"r8", "r9", "r10", "r11", "r12", "r15", "cc", "memory");

#define CLEAN_AVX2 __asm__ volatile ("vzeroupper   \n\t");

#elif defined (_MSC_VER) && defined(_WIN64)

#define COMPUTE_16X_SMALL_PROOTS_AVX2	\
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

	#define COMPUTE_16X_SMALL_NROOTS_AVX2	\
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

#ifdef USE_AVX512BW
#include <immintrin.h>

#define COMPUTE_32X_SMALL_PROOTS_AVX512	\
		do {	\
				__m512i primes;	\
				__m512i root1s;	\
				__m512i root2s;	\
				__m512i ptrs;	\
				__m512i tmp1;	\
				__m512i tmp2;	\
				__m512i tmp3;	\
				__m512i tmp4;	\
                __mmask32 m1, m2; \
				for (j = h.start; j < h.stop - 32; j += 32) \
				{	\
					ptrs = _mm512_loadu_si512((__m512i *)(h.updates + j)); \
					root1s = _mm512_loadu_si512((__m512i *)(h.first_r1 + j)); \
					root2s = _mm512_loadu_si512((__m512i *)(h.first_r2 + j)); \
					primes = _mm512_loadu_si512((__m512i *)(h.primes + j)); \
					m1 = _mm512_cmpgt_epu16_mask(ptrs, root1s);	\
					m2 = _mm512_cmpgt_epu16_mask(ptrs, root2s);	\
					root1s = _mm512_sub_epi16(root1s, ptrs); \
					root2s = _mm512_sub_epi16(root2s, ptrs); \
					root1s = _mm512_mask_add_epi16(root1s, m1, root1s, primes);	\
					root2s = _mm512_mask_add_epi16(root2s, m2, root2s, primes);	\
					tmp1 = _mm512_max_epu16(root1s, root2s); \
					tmp2 = _mm512_min_epu16(root1s, root2s); \
                    tmp3 = _mm512_sub_epi16(primes, tmp1); \
                    tmp4 = _mm512_sub_epi16(primes, tmp2); \
					_mm512_storeu_si512((__m512i *)(h.fbp1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.first_r1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.fbp2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.first_r2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.fbn1 + j), tmp3);	\
					_mm512_storeu_si512((__m512i *)(h.fbn2 + j), tmp4);	\
				}	\
                h.stop = j - 32; \
			} while(0);


#define COMPUTE_32X_SMALL_NROOTS_AVX512	\
		do {	\
				__m512i primes;	\
				__m512i root1s;	\
				__m512i root2s;	\
				__m512i ptrs;	\
				__m512i tmp1;	\
				__m512i tmp2;	\
				__m512i tmp3;	\
				__m512i tmp4;	\
                __mmask32 m1, m2; \
				for (j = h.start; j < h.stop - 32; j += 32) \
				{	\
					ptrs = _mm512_loadu_si512((__m512i *)(h.updates + j)); \
					root1s = _mm512_loadu_si512((__m512i *)(h.first_r1 + j)); \
					root2s = _mm512_loadu_si512((__m512i *)(h.first_r2 + j)); \
					primes = _mm512_loadu_si512((__m512i *)(h.primes + j)); \
					tmp1 = _mm512_add_epi16(root1s, ptrs); \
					tmp2 = _mm512_add_epi16(root2s, ptrs); \
                    m1 = _mm512_cmplt_epu16_mask(tmp1, root1s);	\
					m2 = _mm512_cmplt_epu16_mask(tmp2, root2s);	\
                    m1 |= _mm512_cmpge_epu16_mask(tmp1, primes);	\
					m2 |= _mm512_cmpge_epu16_mask(tmp2, primes);	\
					root1s = _mm512_mask_sub_epi16(tmp1, m1, tmp1, primes);	\
					root2s = _mm512_mask_sub_epi16(tmp2, m2, tmp2, primes);	\
					tmp1 = _mm512_max_epu16(root1s, root2s); \
					tmp2 = _mm512_min_epu16(root1s, root2s); \
                    tmp3 = _mm512_sub_epi16(primes, tmp1); \
                    tmp4 = _mm512_sub_epi16(primes, tmp2); \
					_mm512_storeu_si512((__m512i *)(h.fbp1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.first_r1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.fbp2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.first_r2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.fbn1 + j), tmp3);	\
					_mm512_storeu_si512((__m512i *)(h.fbn2 + j), tmp4);	\
				}	\
                h.stop = j - 32; \
			} while(0);


#define COMPUTE_16X_SMALL_PROOTS_AVX512	\
		do {	\
				__m256i primes;	\
				__m256i root1s;	\
				__m256i root2s;	\
				__m256i ptrs;	\
				__m256i tmp1;	\
				__m256i tmp2;	\
				__m256i tmp3;	\
				__m256i tmp4;	\
                __mmask16 m1, m2; \
				for (j = h.start; j < h.stop; j += 16) \
				{	\
					ptrs = _mm256_load_si256((__m256i *)(h.updates + j)); \
					root1s = _mm256_load_si256((__m256i *)(h.first_r1 + j)); \
					root2s = _mm256_load_si256((__m256i *)(h.first_r2 + j)); \
					primes = _mm256_load_si256((__m256i *)(h.primes + j)); \
					m1 = _mm256_cmpgt_epu16_mask(ptrs, root1s);	\
					m2 = _mm256_cmpgt_epu16_mask(ptrs, root2s);	\
					root1s = _mm256_sub_epi16(root1s, ptrs); \
					root2s = _mm256_sub_epi16(root2s, ptrs); \
					root1s = _mm256_mask_add_epi16(root1s, m1, root1s, primes);	\
					root2s = _mm256_mask_add_epi16(root2s, m2, root2s, primes);	\
					tmp1 = _mm256_max_epu16(root1s, root2s); \
					tmp2 = _mm256_min_epu16(root1s, root2s); \
                    tmp3 = _mm256_sub_epi16(primes, tmp1); \
                    tmp4 = _mm256_sub_epi16(primes, tmp2); \
					_mm256_store_si256((__m256i *)(h.fbp1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.first_r1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.fbp2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.first_r2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.fbn1 + j), tmp3);	\
					_mm256_store_si256((__m256i *)(h.fbn2 + j), tmp4);	\
				}	\
			} while(0);

#define COMPUTE_16X_SMALL_NROOTS_AVX512	\
		do {	\
				__m256i primes;	\
				__m256i root1s;	\
				__m256i root2s;	\
				__m256i ptrs;	\
				__m256i tmp1;	\
				__m256i tmp2;	\
				__m256i tmp3;	\
				__m256i tmp4;	\
                __mmask16 m1, m2; \
				for (j = h.start; j < h.stop; j += 16) \
				{	\
					ptrs = _mm256_load_si256((__m256i *)(h.updates + j)); \
					root1s = _mm256_load_si256((__m256i *)(h.first_r1 + j)); \
					root2s = _mm256_load_si256((__m256i *)(h.first_r2 + j)); \
					primes = _mm256_load_si256((__m256i *)(h.primes + j)); \
					tmp1 = _mm256_add_epi16(root1s, ptrs); \
					tmp2 = _mm256_add_epi16(root2s, ptrs); \
                    m1 = _mm256_cmplt_epu16_mask(tmp1, root1s);	\
					m2 = _mm256_cmplt_epu16_mask(tmp2, root2s);	\
                    m1 |= _mm256_cmpge_epu16_mask(tmp1, primes);	\
					m2 |= _mm256_cmpge_epu16_mask(tmp2, primes);	\
					root1s = _mm256_mask_sub_epi16(tmp1, m1, tmp1, primes);	\
					root2s = _mm256_mask_sub_epi16(tmp2, m2, tmp2, primes);	\
					tmp1 = _mm256_max_epu16(root1s, root2s); \
					tmp2 = _mm256_min_epu16(root1s, root2s); \
                    tmp3 = _mm256_sub_epi16(primes, tmp1); \
                    tmp4 = _mm256_sub_epi16(primes, tmp2); \
					_mm256_store_si256((__m256i *)(h.fbp1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.first_r1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.fbp2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.first_r2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.fbn1 + j), tmp3);	\
					_mm256_store_si256((__m256i *)(h.fbn2 + j), tmp4);	\
				}	\
                h.stop = j - 32; \
			} while(0);

#endif
