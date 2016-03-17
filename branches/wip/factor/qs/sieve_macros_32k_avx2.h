
#if defined(GCC_ASM64X) || defined(__MINGW64__)



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


