

#define FILL_ONE_PRIME_P(i)					\
	if (root1 < interval)					\
	{										\
		bnum = root1 >> 15;			\
		bptr = sliceptr_p +					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & 32767)); \
		numptr_p[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> 15;			\
		bptr = sliceptr_p +					\
			(bnum << BUCKET_BITS) +			\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & 32767)); \
		numptr_p[bnum]++;					\
	}

#define FILL_ONE_PRIME_N(i)						\
	if (root1 < interval)						\
	{											\
		bnum = root1 >> 15;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & 32767)); \
		numptr_n[bnum]++;					\
	}										\
	if (root2 < interval)					\
	{										\
		bnum = root2 >> 15;			\
		bptr = sliceptr_n +					\
			(bnum << BUCKET_BITS) +			\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & 32767)); \
		numptr_n[bnum]++;					\
	}			

#define FILL_ONE_PRIME_LOOP_P(i)				\
	bnum = root1 >> 15;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & 32767)); \
		numptr_p[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> 15;				\
	}											\
	bnum = root2 >> 15;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_p +						\
			(bnum << BUCKET_BITS) +				\
			numptr_p[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & 32767)); \
		numptr_p[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> 15;				\
	} 

#define FILL_ONE_PRIME_LOOP_N(i)				\
	bnum = root1 >> 15;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_n +						\
			(bnum << BUCKET_BITS) +				\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root1 & 32767)); \
		numptr_n[bnum]++;					\
		root1 += prime;							\
		bnum = root1 >> 15;				\
	}											\
	bnum = root2 >> 15;					\
	while (bnum < numblocks)					\
	{											\
		bptr = sliceptr_n +						\
			(bnum << BUCKET_BITS) +				\
			numptr_n[bnum];					\
		*bptr = ((i - bound_val) << 16 | (root2 & 32767)); \
		numptr_n[bnum]++;					\
		root2 += prime;							\
		bnum = root2 >> 15;				\
	} 

#if defined(GCC_ASM64X) || defined(__MINGW64__)
	// The assembly for putting a prime into a bucket is fairly regular, so we break 
	// macros out to make the inline loop shorter within nextRoots().  The only things
	// that change between root1 and root2 are the input root register.  we use two
	// in order to take advantage of the fractional clock latency of movd, and to
	// break up a dependency bottleneck between movd and the cmpl.

    // macro for adding an element (specified by r8d) to the end of a bucket list
    #define UPDATE_ROOT1(it) \
        "movl   %%r8d,%%ebx \n\t"				/* ebx becomes bnum */ \
		"movl   %%r15d,%%edi \n\t"				/* edi becomes fb offset */ \
		"shrl   $15,%%ebx \n\t"	/* right shift root by blksize  = bnum */ \
		"andl   $32767,%%r8d \n\t"	/* root & BLOCKSIZEm1 */ \
		"movl   %%ebx,%%ecx \n\t"				/* ecx becomes bucket address offset */ \
		"addl	$" it ", %%edi \n\t"			/* add iteration number to j */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"subl   %%r12d,%%edi \n\t"				/* j - bound_val */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
		"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"orl	%%r8d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */

	// macro for adding an element (specified by r9d) to the end of a bucket list
	#define UPDATE_ROOT2(it) \
		"movl   %%r9d,%%ebx \n\t"				/* ebx becomes bnum */ \
		"movl   %%r15d,%%edi \n\t"				/* edi becomes fb offset */ \
		"shrl   $15,%%ebx \n\t"	/* right shift root by blksize  = bnum */ \
		"andl   $32767,%%r9d \n\t"	/* root & BLOCKSIZEm1 */ \
		"movl   %%ebx,%%ecx \n\t"				/* ecx becomes bucket address offset */ \
		"addl	$" it ", %%edi \n\t"			/* add iteration number to j */ \
		"movl   (%%r10,%%rbx,4),%%r14d \n\t"	/* numptr_p[bnum] */ \
		"subl   %%r12d,%%edi \n\t"				/* j - bound_val */ \
		"shll   $" BUCKET_BITStxt ",%%ecx \n\t"	/* bnum << BUCKET_BITS */ \
		"shll	$16,%%edi \n\t"					/* move (j - bound_val) to upper word */ \
		"addl   %%r14d,%%ecx \n\t"				/* (bnum << 11) + numptr_p[bnum] */ \
		"addl   $1,(%%r10,%%rbx,4) \n\t"		/* store new numptr to memory */ \
		"orl	%%r9d,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */

	// macro for iteratively adding an element (specified by r8d) to the end of a bucket list
	#define UPDATE_ROOT1_LOOP(it) \
		"1:		\n\t"	 						/* beginning of loop */ \
		"movl   %%r8d,%%ebx \n\t"				/* copy root1 to ebx */ \
		"shrl   $15,%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
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
		"andl   $32767,%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
		"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
		"addl	%%edx,%%r8d \n\t"				/* increment root by prime */ \
		"cmpl   %%r13d,%%r8d \n\t"				/* root > interval? */ \
		"jb		1b \n\t"						/* repeat if necessary */

	// macro for iteratively adding an element (specified by r9d) to the end of a bucket list
	#define UPDATE_ROOT2_LOOP(it) \
		"1:		\n\t"	 						/* beginning of loop */ \
		"movl   %%r9d,%%ebx \n\t"				/* copy root to ebx */ \
		"shrl   $15,%%ebx \n\t"	/* right shift root1 by 15  = bnum */ \
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
		"andl   $32767,%%eax \n\t"	/* root1 & BLOCKSIZEm1 */ \
		"orl	%%eax,%%edi \n\t"				/* combine two words to reduce write port pressure */ \
		"movl   %%edi,(%%r11,%%rcx,4) \n\t"		/* store new fb_index/loc to memory */ \
		"addl	%%edx,%%r9d \n\t"				/* increment root by prime */ \
		"cmpl   %%r13d,%%r9d \n\t"				/* root > interval? */ \
		"jb		1b \n\t"						/* repeat if necessary */

#endif

