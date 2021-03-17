

#define DIVIDE_ONE_PRIME \
	do \
	{						\
		fb_offsets[++smooth_num] = i;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime); \
	} while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0); 
#define DIVIDE_ONE_PRIME_2(j) \
	do \
    	{						\
		fb_offsets[++smooth_num] = (j);	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[j]); \
        } while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[j]) == 0); 


#define DIVIDE_RESIEVED_PRIME(j) \
	while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[i+j]) == 0) \
	{						\
		fb_offsets[++smooth_num] = i+j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[i+j]);		\
	}
#define DIVIDE_RESIEVED_PRIME_2(j) \
    while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[j]) == 0) \
    {						\
	    fb_offsets[++smooth_num] = j;	\
	    mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[j]);		\
    }

#define DIVIDE_LARGE_PRIME(j) \
    while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0) \
    {						\
		fb_offsets[++smooth_num] = j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime);		\
    }

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	#define TDIV_MED_CLEAN asm volatile("emms");
#elif defined(MSC_ASM32A)
	#define TDIV_MED_CLEAN ASM_M {emms};
#else
	#define TDIV_MED_CLEAN /*nothing*/
#endif




#ifdef USE_AVX512BW



#define MOD_CMP_32X_vec(v_bits)				{ \
        __m512i v_primes, v_y4, v_y7, v_r1, v_r2, v_corr, v_inv;       \
        __mmask32 m_r1, m_r2; \
        uint32_t msk32; \
        v_corr = _mm512_load_si512((__m512i *)(fullfb_ptr->correction + i)); \
        v_inv = _mm512_load_si512((__m512i *)(fullfb_ptr->small_inv + i)); \
        v_primes = _mm512_load_si512((__m512i *)(fbc->prime + i)); \
        v_r1 = _mm512_load_si512((__m512i *)(fbc->root1 + i)); \
        v_r2 = _mm512_load_si512((__m512i *)(fbc->root2 + i)); \
        v_y4 = _mm512_sub_epi16(v_blksz, v_blkloc); \
        v_y4 = _mm512_add_epi16(v_y4, v_corr); \
        v_y4 = _mm512_mulhi_epu16(v_y4, v_inv); \
        v_y4 = _mm512_srlv_epi16(v_y4, v_bits); \
        v_y7 = _mm512_add_epi16(v_blkloc, v_primes); \
        v_y4 = _mm512_mullo_epi16(v_y4, v_primes); \
        v_y7 = _mm512_sub_epi16(v_y7, v_blksz); \
        v_y4 = _mm512_add_epi16(v_y7, v_y4); \
        m_r1 = _mm512_cmpeq_epi16_mask(v_y4, v_r1); \
        m_r2 = _mm512_cmpeq_epi16_mask(v_y4, v_r2); \
        msk32 = m_r1 | m_r2; \
        while ((pos = _trail_zcnt(msk32)) < 32) { \
            buffer[tmp3++] = pos + i; \
            msk32 = _reset_lsb(msk32); \
        }}

#define MOD_CMP_8X_vec_bw(xtra_bits)				{ \
        __m128i v_primes, v_y4, v_y7, v_r1, v_r2;       \
        __mmask8 m_r1, m_r2; \
        uint32_t msk32; \
        v_primes = _mm_load_si128((__m128i *)(fbc->prime + i)); \
        v_y4 = _mm_mask_sub_epi16(v_blksz128, 0xff, v_blksz128, v_blkloc128); \
        v_y4 = _mm_mask_add_epi16(v_y4, 0xff, v_y4, _mm_load_si128((__m128i *)(fullfb_ptr->correction + i))); \
        v_r1 = _mm_load_si128((__m128i *)(fbc->root1 + i)); \
        v_y4 = _mm_mask_mulhi_epu16(v_y4, 0xff, v_y4, _mm_load_si128((__m128i *)(fullfb_ptr->small_inv + i))); \
        v_r2 = _mm_load_si128((__m128i *)(fbc->root2 + i)); \
        v_y4 = _mm_mask_srli_epi16 (v_y4, 0xff, v_y4, xtra_bits); \
        v_y7 = _mm_mask_add_epi16(v_blkloc128, 0xff, v_blkloc128, v_primes); \
        v_y4 = _mm_mask_mullo_epi16(v_y4, 0xff, v_y4, v_primes); \
        v_y7 = _mm_mask_sub_epi16(v_y7, 0xff, v_y7, v_blksz128); \
        v_y4 = _mm_mask_add_epi16(v_y7, 0xff, v_y7, v_y4); \
        m_r1 = _mm_cmpeq_epu16_mask(v_y4, v_r1); \
        m_r2 = _mm_cmpeq_epu16_mask(v_y4, v_r2); \
        msk32 = m_r1 | m_r2; \
        while ((pos = _trail_zcnt(msk32)) < 32) { \
            buffer[tmp3++] = pos + i; \
            msk32 = _reset_lsb(msk32); \
        }}

#define MOD_INIT_32X									\
        __m512i v_blksz = _mm512_load_si512((__m512i *)bl_sizes);  \
        __m512i v_blkloc = _mm512_load_si512((__m512i *)bl_locs);  \
        __m128i v_blksz128 = _mm_load_si128((__m128i *)bl_sizes);  \
        __m128i v_blkloc128 = _mm_load_si128((__m128i *)bl_locs);  \
        uint32_t msk32, pos;



#define INIT_CORRECTIONS_32 \
	corrections[0] = 32768 - block_loc; \
	corrections[1] = 32768 - block_loc; \
	corrections[2] = 32768 - block_loc; \
	corrections[3] = 32768 - block_loc; \
	corrections[4] = 32768 - block_loc; \
	corrections[5] = 32768 - block_loc; \
	corrections[6] = 32768 - block_loc; \
	corrections[7] = 32768 - block_loc; \
	corrections[8] = 32768 - block_loc; \
	corrections[9] = 32768 - block_loc; \
	corrections[10] = 32768 - block_loc; \
	corrections[11] = 32768 - block_loc; \
	corrections[12] = 32768 - block_loc; \
	corrections[13] = 32768 - block_loc; \
	corrections[14] = 32768 - block_loc; \
	corrections[15] = 32768 - block_loc; \
    corrections[16] = 32768 - block_loc; \
    corrections[17] = 32768 - block_loc; \
    corrections[18] = 32768 - block_loc; \
    corrections[19] = 32768 - block_loc; \
    corrections[20] = 32768 - block_loc; \
    corrections[21] = 32768 - block_loc; \
    corrections[22] = 32768 - block_loc; \
    corrections[23] = 32768 - block_loc; \
    corrections[24] = 32768 - block_loc; \
    corrections[25] = 32768 - block_loc; \
    corrections[26] = 32768 - block_loc; \
    corrections[27] = 32768 - block_loc; \
    corrections[28] = 32768 - block_loc; \
    corrections[29] = 32768 - block_loc; \
    corrections[30] = 32768 - block_loc; \
    corrections[31] = 32768 - block_loc;

#define STEP_COMPARE_COMBINE_AVX512 \
	v512_y2 = _mm512_sub_epi16(v512_y2, v512_p); \
    v512_y3 = _mm512_sub_epi16(v512_y3, v512_p); \
    m1 |= _mm512_cmpeq_epi16_mask(v512_y2, v512_y5); \
    m2 |= _mm512_cmpeq_epi16_mask(v512_y3, v512_y6);

#define INIT_RESIEVE_AVX512 \
    __m512i v512_p, v512_y0, v512_y2, v512_y3, v512_corr, v512_y5, v512_y6, v512_y7; \
    uint32_t msk32, pos; \
    __mmask32 m1 = 0, m2 = 0; \
    v512_p = v512_y0 = v512_y2 = v512_y3 = v512_corr = v512_y5 = v512_y6 = v512_y7 = _mm512_set1_epi16(0); \
    v512_corr = _mm512_load_si512(corrections);                   \
    v512_y0 = _mm512_xor_si512(v512_y0, v512_y0);               \
    v512_y2 = _mm512_load_si512(fbc->root1 + i);                \
    v512_y3 = _mm512_load_si512(fbc->root2 + i);                \
    v512_y2 = _mm512_add_epi16(v512_corr, v512_y2);               \
    v512_y3 = _mm512_add_epi16(v512_corr, v512_y3);               \
    v512_p = _mm512_load_si512(fbc->prime + i);                 \
    v512_y5 = _mm512_xor_si512(v512_y5, v512_y5);               \
    v512_y6 = _mm512_xor_si512(v512_y6, v512_y6);               \
    v512_y7 = _mm512_xor_si512(v512_y7, v512_y7);


#define GATHER_BIT_INDICES_AVX512 \
        while ((pos = _trail_zcnt(msk32)) < 32) { \
            buffer[result++] = pos + i; \
            msk32 = _reset_lsb(msk32); \
        }


#define RESIEVE_32X_14BIT_MAX_VEC_AVX512 \
			INIT_RESIEVE_AVX512 \
			STEP_COMPARE_COMBINE_AVX512	\
			STEP_COMPARE_COMBINE_AVX512	\
			STEP_COMPARE_COMBINE_AVX512	\
			STEP_COMPARE_COMBINE_AVX512	\
			msk32 = m1 ^ m2; \
            GATHER_BIT_INDICES_AVX512

#define RESIEVE_32X_15BIT_MAX_VEC_AVX512 \
			INIT_RESIEVE_AVX512 \
			STEP_COMPARE_COMBINE_AVX512	\
			STEP_COMPARE_COMBINE_AVX512	\
			msk32 = m1 ^ m2; \
            GATHER_BIT_INDICES_AVX512

#define RESIEVE_32X_16BIT_MAX_VEC_AVX512 \
			INIT_RESIEVE_AVX512 \
			STEP_COMPARE_COMBINE_AVX512	\
			msk32 = m1 ^ m2; \
            GATHER_BIT_INDICES_AVX512


#define STEP_COMPARE_COMBINE_8X_AVX512 \
        v128_x2 = _mm_sub_epi16(v128_x2, v128_p); \
        v128_x3 = _mm_sub_epi16(v128_x3, v128_p); \
        m1 |= _mm_cmpeq_epi16_mask(v128_x2, v128_x5); \
        m2 |= _mm_cmpeq_epi16_mask(v128_x3, v128_x6);


#define INIT_RESIEVE_8X_AVX512 \
__m128i v128_p, v128_x0, v128_x2, v128_x3, v128_x4, v128_x5, v128_x6, v128_x7; \
uint32_t msk32, pos; \
__mmask8 m1 = 0, m2 = 0; \
v128_x4 = _mm_load_si128((__m128i *)corrections);                  \
v128_x0 = _mm_xor_si128(v128_x0, v128_x0);              \
v128_x2 = _mm_load_si128((__m128i *)(fbc->root1 + i));               \
v128_x3 = _mm_load_si128((__m128i *)(fbc->root2 + i));               \
v128_x2 = _mm_add_epi16(v128_x4, v128_x2);              \
v128_x3 = _mm_add_epi16(v128_x4, v128_x3);              \
v128_p = _mm_load_si128((__m128i *)(fbc->prime + i));                \
v128_x5 = _mm_xor_si128(v128_x5, v128_x5);              \
v128_x6 = _mm_xor_si128(v128_x6, v128_x6);              \
v128_x7 = _mm_xor_si128(v128_x7, v128_x7);              


#define RESIEVE_8X_14BIT_MAX_VEC_AVX512 \
			INIT_RESIEVE_8X_AVX512 \
			STEP_COMPARE_COMBINE_8X_AVX512	\
			STEP_COMPARE_COMBINE_8X_AVX512	\
			STEP_COMPARE_COMBINE_8X_AVX512	\
			STEP_COMPARE_COMBINE_8X_AVX512	\
			msk32 = m1 ^ m2; \
            GATHER_BIT_INDICES_AVX512


#define RESIEVE_8X_16BIT_MAX_VEC_AVX512 \
			INIT_RESIEVE_8X_AVX512 \
			STEP_COMPARE_COMBINE_8X_AVX512	\
			msk32 = m1 ^ m2; \
            GATHER_BIT_INDICES_AVX512

#endif
