
#ifdef USE_YAFU_TDIV
#define DIVIDE_ONE_PRIME \
	do	\
	{	\
		fb_offsets[++smooth_num] = i;	\
		zShortDiv32(tmp32, prime, tmp32);	\
	} while (zShortMod32(tmp32, prime) == 0);
#else
#define DIVIDE_ONE_PRIME \
	do \
	{						\
		fb_offsets[++smooth_num] = i;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], prime); \
	} while (mpz_tdiv_ui(dconf->Qvals[report_num], prime) == 0); 
#endif

#ifdef USE_YAFU_TDIV
#define DIVIDE_RESIEVED_PRIME(j) \
	while (zShortMod32(tmp32, fbc->prime[i+j]) == 0)	\
	{	\
		fb_offsets[++smooth_num] = i+j;	\
		zShortDiv32(tmp32, fbc->prime[i+j], tmp32);	\
	}
#else
#define DIVIDE_RESIEVED_PRIME(j) \
	while (mpz_tdiv_ui(dconf->Qvals[report_num], fbc->prime[i+j]) == 0) \
	{						\
		fb_offsets[++smooth_num] = i+j;	\
		mpz_tdiv_q_ui(dconf->Qvals[report_num], dconf->Qvals[report_num], fbc->prime[i+j]);		\
	}
#endif

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	#define TDIV_MED_CLEAN asm volatile("emms");
#elif defined(MSC_ASM32A)
	#define TDIV_MED_CLEAN ASM_M {emms};
#else
	#define TDIV_MED_CLEAN /*nothing*/
#endif
