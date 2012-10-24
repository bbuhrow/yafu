
#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)

		// the original roots are the current roots + BLOCKSIZE
		// to test if this block_loc is on the root progression, we first
		// advance the current block_loc one step past blocksize and
		// then test if this value is equal to either root

		// to advance the block_loc, we compute 
		// steps = 1 + (BLOCKSIZE - block_loc) / prime

		// to do the div quickly using precomputed values, do:
		// tmp = ((BLOCKSIZE - block_loc) + correction) * inv >> shift
		// shift is either 24 or 26 bits, depending on the size of prime
		// in order to keep enough precision.  See Agner Fog's optimization
		// manuals.
		// also note that with 64k versions, the final addition will overflow,
		// but since the addition does not saturate and since the final 
		// subtraction always yeilds a number less than 2^16, the overflow
			// does not hurt.

		// will miss cases like this:
		//	index = 56, prime = 587, inv = 28581, corr = 1, root1 = 0, root2 = 322, 
		// loc = 1070, result = 0
		// which we can catch by comparing to prime as well as root1/2; but it probably isn't
		// worth it.

				/*"movdqa (%6), %%xmm0 \n\t"		 move in BLOCKSIZE */ \
				/*"movdqa (%7), %%xmm1 \n\t"		 move in block_loc */ \

	#define MOD_INIT_8X												\
		ASM_G (														\
			"movdqa (%0), %%xmm0 \n\t"		/* move in BLOCKSIZE */ \
			"movdqa (%1), %%xmm1 \n\t"		/* move in block_loc */ \
			:														\
			: "r" (bl_sizes), "r" (bl_locs)							\
			: "xmm0", "xmm1");

	#define MOD_CMP_8X(xtra_bits)																		\
		ASM_G (																				\
			"movdqa %%xmm1, %%xmm7 \n\t"	/* copy block_loc */ \
			"movdqa	%%xmm0, %%xmm4 \n\t"	/* copy BLOCKSIZE */ \
			"movdqa (%1), %%xmm3 \n\t"		/* move in primes */							\
			"psubw	%%xmm1, %%xmm4 \n\t"	/* BLOCKSIZE - block_loc */						\
			"paddw	(%2), %%xmm4 \n\t"		/* apply corrections */							\
			"movdqa (%4), %%xmm6 \n\t"		/* move in root1s */							\
			"pmulhuw	(%3), %%xmm4 \n\t"	/* (unsigned) multiply by inverses */		\
			"movdqa (%5), %%xmm2 \n\t"		/* move in root2s */							\
			"psrlw	$" xtra_bits ", %%xmm4 \n\t"		/* to get to total shift of 24/26/28 bits */			\
			"paddw	%%xmm3, %%xmm7 \n\t"	/* add primes and block_loc */					\
			"pmullw	%%xmm3, %%xmm4 \n\t"	/* (signed) multiply by primes */				\
			"psubw	%%xmm0, %%xmm7 \n\t"	/* substract blocksize */						\
			"paddw	%%xmm7, %%xmm4 \n\t"	/* add in block_loc + primes - blocksize */		\
			"pcmpeqw	%%xmm4, %%xmm6 \n\t"	/* compare to root1s */						\
			"pcmpeqw	%%xmm4, %%xmm2 \n\t"	/* compare to root2s */						\
			"por	%%xmm6, %%xmm2 \n\t"	/* combine compares */							\
			"pmovmskb %%xmm2, %0 \n\t"		/* export to result */							\
			: "=r" (tmp3)																	\
			: "r" (fbc->prime + i), "r" (fullfb_ptr->correction + i), \
				"r" (fullfb_ptr->small_inv + i), "r" (fbc->root1 + i), \
					"r" (fbc->root2 + i) \
			: "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7");

	#define STEP_COMPARE_COMBINE \
		"psubw %%xmm1, %%xmm2 \n\t"		/* subtract primes from root1s */ \
		"psubw %%xmm1, %%xmm3 \n\t"		/* subtract primes from root2s */ \
		"pcmpeqw %%xmm2, %%xmm5 \n\t"	/* root1s ?= 0 */ \
		"pcmpeqw %%xmm3, %%xmm6 \n\t"	/* root2s ?= 0 */ \
		"por %%xmm5, %%xmm7 \n\t"		/* combine results */ \
		"por %%xmm6, %%xmm0 \n\t"		/* combine results */

	#define INIT_RESIEVE \
		"movdqa (%4), %%xmm4 \n\t"		/* bring in corrections to roots */				\
		"pxor %%xmm0, %%xmm0 \n\t"		/* zero xmm8 */ \
		"movdqa (%2), %%xmm2 \n\t"		/* bring in 8 root1s */ \
		"paddw %%xmm4, %%xmm2 \n\t"		/* correct root1s */ \
		"movdqa (%3), %%xmm3 \n\t"		/* bring in 8 root2s */ \
		"paddw %%xmm4, %%xmm3 \n\t"		/* correct root2s */ \
		"movdqa (%1), %%xmm1 \n\t"		/* bring in 8 primes */ \
		"pxor %%xmm7, %%xmm7 \n\t"		/* zero xmm7 */ \
		"pxor %%xmm5, %%xmm5 \n\t"		/* zero xmm5 */ \
		"pxor %%xmm6, %%xmm6 \n\t"		/* zero xmm6 */

	#define RESIEVE_8X_14BIT_MAX \
		asm ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"por	%%xmm0, %%xmm7 \n\t" \
			"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
			: "=r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
			: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

	#define RESIEVE_8X_15BIT_MAX \
		asm ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"por	%%xmm0, %%xmm7 \n\t" \
			"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
			: "=r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
			: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

	#define RESIEVE_8X_16BIT_MAX \
		asm ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			"por	%%xmm0, %%xmm7 \n\t" \
			"pmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
			: "=r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
			: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);



#elif defined(_MSC_VER)

	#define MOD_INIT_8X	

	#define STEP_COMPARE_COMBINE \
		root1s = _mm_sub_epi16(root1s, primes); \
		root2s = _mm_sub_epi16(root2s, primes); \
		tmp1 = _mm_cmpeq_epi16(tmp1, root1s); \
		tmp2 = _mm_cmpeq_epi16(tmp2, root2s); \
		combine = _mm_xor_si128(combine, tmp1); \
		combine = _mm_xor_si128(combine, tmp2);

	#define INIT_RESIEVE \
		c = _mm_load_si128((__m128i *)corrections); \
		root1s = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
		root1s = _mm_add_epi16(root1s, c); \
		root2s = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
		root2s = _mm_add_epi16(root2s, c); \
		primes = _mm_load_si128((__m128i *)(fbc->prime + i)); \
		combine = _mm_xor_si128(combine, combine); \
		tmp1 = _mm_xor_si128(tmp1, tmp1); \
		tmp2 = _mm_xor_si128(tmp2, tmp2);

	#define RESIEVE_8X_14BIT_MAX \
		do { \
			__m128i tmp1;	\
			__m128i tmp2;	\
			__m128i root1s;	\
			__m128i root2s;	\
			__m128i primes;	\
			__m128i c;	\
			__m128i combine;	\
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			result = _mm_movemask_epi8(combine); \
		} while (0);

	#define RESIEVE_8X_15BIT_MAX \
		do { \
			__m128i tmp1;	\
			__m128i tmp2;	\
			__m128i root1s;	\
			__m128i root2s;	\
			__m128i primes;	\
			__m128i c;	\
			__m128i combine;	\
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			result = _mm_movemask_epi8(combine); \
		} while (0);

	#define RESIEVE_8X_16BIT_MAX \
		do { \
			__m128i tmp1;	\
			__m128i tmp2;	\
			__m128i root1s;	\
			__m128i root2s;	\
			__m128i primes;	\
			__m128i c;	\
			__m128i combine;	\
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			result = _mm_movemask_epi8(combine); \
		} while (0);


	#define MOD_CMP_8X(xtra_bits)						\
		do { \
				__m128i t1; \
				__m128i t2;	\
				__m128i lr1;	\
				__m128i lr2;	\
				__m128i lp;	\
				__m128i c;	\
				__m128i sinv; \
				__m128i blksz; \
				__m128i blkloc; \
			blksz = _mm_load_si128((__m128i *)(bl_sizes)); \
			blkloc = _mm_load_si128((__m128i *)(bl_locs)); \
			t1 = blksz; \
			lp = _mm_load_si128((__m128i *)(fbc->prime + i)); \
			t1 = _mm_sub_epi16(t1, blkloc); \
			c = _mm_load_si128((__m128i *)(fullfb_ptr->correction + i));	\
			t2 = blkloc; \
			sinv = _mm_load_si128((__m128i *)(fullfb_ptr->small_inv + i));	\
			c = _mm_add_epi16(c, t1); \
			lr1 = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
			c = _mm_mulhi_epu16(c, sinv); \
			lr2 = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
			c = _mm_srli_epi16(c, xtra_bits); \
			t2 = _mm_add_epi16(t2, lp); \
			c = _mm_mullo_epi16(c, lp); \
			c = _mm_add_epi16(c, t2); \
			c = _mm_sub_epi16(c, blksz); \
			lr1 = _mm_cmpeq_epi16(lr1, c); \
			lr2 = _mm_cmpeq_epi16(lr2, c); \
			lr2 = _mm_or_si128(lr2, lr1); \
			tmp3 = _mm_movemask_epi8(lr2); \
		} while (0);

#else	/* compiler not recognized*/

	#define MOD_INIT_8X	

	#define COMPARE_RESIEVE_VALS(x)	\
		if (r1 == 0) result |= x; \
		if (r2 == 0) result |= x;

	#define STEP_RESIEVE \
		r1 -= p;	\
		r2 -= p;

	#define RESIEVE_1X_14BIT_MAX(x)	\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);

	#define RESIEVE_1X_15BIT_MAX(x)	\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);			\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);

	#define RESIEVE_1X_16BIT_MAX(x)	\
		STEP_RESIEVE;						\
		COMPARE_RESIEVE_VALS(x);

	#define RESIEVE_8X_14BIT_MAX \
		do { \
			int p = (int)fbc->prime[i];										\
			int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
			int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
			RESIEVE_1X_14BIT_MAX(0x2);										\
			\
			p = (int)fbc->prime[i+1];										\
			r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_14BIT_MAX(0x8);										\
			\
			p = (int)fbc->prime[i+2];										\
			r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_14BIT_MAX(0x20);									\
			\
			p = (int)fbc->prime[i+3];										\
			r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_14BIT_MAX(0x80);									\
			\
			p = (int)fbc->prime[i+4];										\
			r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_14BIT_MAX(0x200);									\
			\
			p = (int)fbc->prime[i+5];										\
			r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_14BIT_MAX(0x800);									\
			\
			p = (int)fbc->prime[i+6];										\
			r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_14BIT_MAX(0x2000);									\
			\
			p = (int)fbc->prime[i+7];										\
			r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_14BIT_MAX(0x8000);									\
		} while (0); 

	#define RESIEVE_8X_15BIT_MAX \
		do { \
			int p = (int)fbc->prime[i];										\
			int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
			int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
			RESIEVE_1X_15BIT_MAX(0x2);										\
			\
			p = (int)fbc->prime[i+1];										\
			r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_15BIT_MAX(0x8);										\
			\
			p = (int)fbc->prime[i+2];										\
			r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_15BIT_MAX(0x20);									\
			\
			p = (int)fbc->prime[i+3];										\
			r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_15BIT_MAX(0x80);									\
			\
			p = (int)fbc->prime[i+4];										\
			r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_15BIT_MAX(0x200);									\
			\
			p = (int)fbc->prime[i+5];										\
			r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_15BIT_MAX(0x800);									\
			\
			p = (int)fbc->prime[i+6];										\
			r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_15BIT_MAX(0x2000);									\
			\
			p = (int)fbc->prime[i+7];										\
			r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_15BIT_MAX(0x8000);									\
		} while (0); 

	#define RESIEVE_8X_16BIT_MAX \
		do { \
			int p = (int)fbc->prime[i];										\
			int r1 = (int)fbc->root1[i] + BLOCKSIZE - block_loc;			\
			int r2 = (int)fbc->root2[i] + BLOCKSIZE - block_loc;			\
			RESIEVE_1X_16BIT_MAX(0x2);										\
			\
			p = (int)fbc->prime[i+1];										\
			r1 = (int)fbc->root1[i+1] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+1] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_16BIT_MAX(0x8);										\
			\
			p = (int)fbc->prime[i+2];										\
			r1 = (int)fbc->root1[i+2] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+2] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_16BIT_MAX(0x20);									\
			\
			p = (int)fbc->prime[i+3];										\
			r1 = (int)fbc->root1[i+3] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+3] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_16BIT_MAX(0x80);									\
			\
			p = (int)fbc->prime[i+4];										\
			r1 = (int)fbc->root1[i+4] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+4] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_16BIT_MAX(0x200);									\
			\
			p = (int)fbc->prime[i+5];										\
			r1 = (int)fbc->root1[i+5] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+5] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_16BIT_MAX(0x800);									\
			\
			p = (int)fbc->prime[i+6];										\
			r1 = (int)fbc->root1[i+6] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+6] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_16BIT_MAX(0x2000);									\
			\
			p = (int)fbc->prime[i+7];										\
			r1 = (int)fbc->root1[i+7] + BLOCKSIZE - block_loc;				\
			r2 = (int)fbc->root2[i+7] + BLOCKSIZE - block_loc;				\
			RESIEVE_1X_16BIT_MAX(0x8000);									\
		} while (0); 

#endif

