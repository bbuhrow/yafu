
#if defined(GCC_ASM64X) || defined(__MINGW64__)

#define _FINALIZE_SORT_UPDATE_SSE41 \
	"paddw	%%xmm4, %%xmm1 \n\t"			/* r1 = r1 + (p - b) */ \
	"paddw	%%xmm5, %%xmm2 \n\t"			/* r2 = r2 + (p - b) */ \
	"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1s */ \
	"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2s */ \
	"pmaxuw	%%xmm4, %%xmm2 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	"pminuw	%%xmm5, %%xmm1 \n\t"			/* replace xmm1 with min of root1 and root2 */ \
	"movdqa %%xmm1, (%2) \n\t"				/* write new root1's */ \
	"movdqa %%xmm2, (%3) \n\t"				/* write new root2's */

#define _SSE41_SMALL_PRIME_SIEVE \
	for (; i < med_B; i += 8) \
		asm (			\
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$15, %%rsi \n\t"				/* logp's range from 15 to 16... call 'em = 15 */ \
			"movdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"movdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"movdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
			\
			_8P_FINAL_STEP_SIEVE					\
			_FINALIZE_SORT_UPDATE_SSE41 \
			\
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx", "rdx", \
				"rsi", "rdi", "cc", "memory");

#elif defined(_MSC_VER)


#define _FINALIZE_SORT_UPDATE_SSE41 \
	root1s = _mm_add_epi16(root1s, tmp1); \
	root2s = _mm_add_epi16(root2s, tmp2); \
	tmp1 = root1s; \
	tmp2 = root2s; \
	root2s = _mm_max_epu16(root2s, tmp1); \
	root1s = _mm_min_epu16(root1s, tmp2); \
	_mm_store_si128((__m128i *)(fb->root1 + i), root1s); \
	_mm_store_si128((__m128i *)(fb->root2 + i), root2s);

#define _SSE41_SMALL_PRIME_SIEVE \
		for (; i < med_B; i += 8) { \
			logp = 15; \
			primes = _mm_load_si128((__m128i *)(fb->prime + i)); \
			root1s = _mm_load_si128((__m128i *)(fb->root1 + i)); \
			root2s = _mm_load_si128((__m128i *)(fb->root2 + i)); \
			_8P_FINAL_STEP_SIEVE	 \
			_FINALIZE_SORT_UPDATE_SSE41 \
		} \
	} while (0);


#endif


