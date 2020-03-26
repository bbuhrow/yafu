

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

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	#define TDIV_MED_CLEAN asm volatile("emms");
#elif defined(MSC_ASM32A)
	#define TDIV_MED_CLEAN ASM_M {emms};
#else
	#define TDIV_MED_CLEAN /*nothing*/
#endif
