
#if defined(GCC_ASM64X) || defined(__MINGW64__)

#define _AVX2_SMALL_PRIME_SIEVE_32k_DIV3 \
	for (; i < full_fb->fb_32k_div3-16; i += 16) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$13, %%rsi \n\t"				/* logps = 13 in this range */ \
			"vmovdqa	(%1), %%ymm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%ymm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%ymm2 \n\t"				/* bring in 8 root2's */ \
			 \
			_16P_STEP_SIEVE_AVX2 \
			_16P_STEP_SIEVE_AVX2 \
			_16P_STEP_SIEVE_AVX2 \
			_16P_FINAL_STEP_SIEVE_AVX2	 \
			_FINALIZE_SORT_UPDATE_AVX2 \
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", \
				"rcx", "rdx", "rsi", "rdi", "cc", "memory");
	
#define _AVX2_SMALL_PRIME_SIEVE_14b \
	for (; i < full_fb->fb_14bit_B-16; i += 16) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$14, %%rsi \n\t"				/* logp's range from 13 to 14... call 'em = 14 */ \
			"vmovdqa	(%1), %%ymm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%ymm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%ymm2 \n\t"				/* bring in 8 root2's */ \
			 \
			_16P_STEP_SIEVE_AVX2 \
			_16P_STEP_SIEVE_AVX2 \
			_16P_FINAL_STEP_SIEVE_AVX2	\
			_FINALIZE_SORT_UPDATE_AVX2 \
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", \
				"rcx", "rdx", "rsi", "rdi", "cc", "memory");

#define _AVX2_SMALL_PRIME_SIEVE_15b \
	for (; i < full_fb->fb_15bit_B-16; i += 16) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$15, %%rsi \n\t"				/* logp's range from 14 to 15... call 'em = 15 */ \
			"vmovdqa	(%1), %%ymm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%ymm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%ymm2 \n\t"				/* bring in 8 root2's */ \
				\
			_16P_STEP_SIEVE_AVX2 \
			_16P_FINAL_STEP_SIEVE_AVX2			\
			_FINALIZE_SORT_UPDATE_AVX2 \
			\
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx", "rdx", "rsi", \
				"rdi", "cc", "memory");

#define _16P_STEP_SIEVE_AVX2 \
	"vpextrw	$0, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$0, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$1, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$1, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$2, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$2, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$3, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$3, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$4, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$4, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$5, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$5, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$6, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$6, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$7, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$7, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"vpsrldq	$16, %%ymm1, %%ymm4 \n\t"		/* next 8 words */ \
	"vpsrldq	$16, %%ymm2, %%ymm5 \n\t"		/* next 8 words */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$0, %%xmm4, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$0, %%xmm5, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$1, %%xmm4, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$1, %%xmm5, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$2, %%xmm4, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$2, %%xmm5, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$3, %%xmm4, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$3, %%xmm5, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$4, %%xmm4, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$4, %%xmm5, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$5, %%xmm4, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$5, %%xmm5, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$6, %%xmm4, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$6, %%xmm5, %%ebx \n\t"			/* extract root2 */ \
	"vpextrw	$7, %%xmm4, %%ecx \n\t"			/* extract root1 */ \
	"vpextrw	$7, %%xmm5, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpaddw	%%ymm0, %%ymm1, %%ymm1 \n\t"			/* increment root1's by primes */ \
	"vpaddw	%%ymm0, %%ymm2, %%ymm2 \n\t"			/* increment root2's by primes */

#define _16P_FINAL_STEP_SIEVE_AVX2 \
	"vmovdqa	%%ymm1, %%ymm4 \n\t"			/* copy root1's before comparison */ \
	"vmovdqa	%%ymm2, %%ymm5 \n\t"			/* copy root1's before comparison */ \
	/* workaround for lack of unsigned greater than... */ \
	"vpsubusw	%%ymm3, %%ymm4, %%ymm4 \n\t"		/* ymm4 := root1 - 32767 */ \
	"vpcmpeqw	%%ymm6, %%ymm4, %%ymm4 \n\t"		/* ymm4 := root1 <= 32767 ? 1 : 0 */ \
	"vpsubusw	%%ymm3, %%ymm5, %%ymm5 \n\t"		/* ymm5 := root2 - 32767 */ \
	"vpcmpeqw	%%ymm6, %%ymm5, %%ymm5 \n\t"		/* ymm5 := root2 <= 32767 ? 1 : 0 */ \
	"vpmovmskb	%%ymm4, %%ecx \n\t" \
	"vpmovmskb	%%ymm5, %%edi \n\t" \
	"testl	$0x2, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* if bits are equal to zero then root1 ! <= blocksize-1; jump */ \
	"vpextrw	$0, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$0, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$1, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$1, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x20, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$2, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x20, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$2, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x80, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$3, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x80, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$3, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x200, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$4, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x200, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$4, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x800, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$5, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x800, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$5, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x2000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$6, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$6, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$7, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$7, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	"vpsrldq	$16, %%ymm1, %%ymm6 \n\t"		/* next 8 words */ \
	"vpsrldq	$16, %%ymm2, %%ymm3 \n\t"		/* next 8 words */ \
	"testl	$0x20000, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* if bits are equal to zero then root1 ! <= blocksize-1; jump */ \
	"vpextrw	$0, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x20000, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$0, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x80000, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$1, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x80000, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$1, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x200000, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"vpand	%%ymm0, %%ymm4, %%ymm4 \n\t"			/* clear primes whose root1's are >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$2, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x200000, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$2, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x800000, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"vpand	%%ymm0, %%ymm5, %%ymm5 \n\t"			/* clear primes whose root2's are >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$3, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x800000, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$3, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x2000000, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"vpsubw	%%ymm7, %%ymm4, %%ymm4 \n\t"			/* advance root1's still below blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$4, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2000000, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$4, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8000000, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"vpsubw	%%ymm7, %%ymm5, %%ymm5 \n\t"			/* advance root1's still below blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$5, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8000000, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$5, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x20000000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$6, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x20000000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$6, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x80000000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$7, %%xmm6, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x80000000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$7, %%xmm3, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */

#define _INIT_AVX2_SMALL_PRIME_SIEVE \
	asm (	\
		"movl	$0x7fff7fff, %%ecx \n\t"		/* load 2 copies of blocksize-1 in a 32bit reg */ \
		"movl	$0x80008000, %%edi \n\t"		/* load 2 copies of blocksize in a 32bit reg */ \
		"movd	%%ecx, %%xmm3 \n\t"				/* we need 32k-1 because we do signed comparisons */ \
		"movd	%%edi, %%xmm7 \n\t"				\
		"vpxor	%%ymm6, %%ymm6, %%ymm6 \n\t"			/* ymm6 := 0 */ \
		"vpshufd	$0, %%ymm3, %%ymm3 \n\t"		/* broadcast blocksize-1 to all words of ymm3 */ \
		"vpshufd	$0, %%ymm7, %%ymm7 \n\t"		/* broadcast blocksize to all words of ymm7 */ \
		: \
		:  \
		: "ecx", "edi", "xmm3", "xmm6", "xmm7");

#define _FINALIZE_SORT_UPDATE_AVX2 \
	"vpaddw	%%ymm4, %%ymm1, %%ymm1 \n\t"			/* r1 = r1 + (p - b) */ \
	"vpaddw	%%ymm5, %%ymm2, %%ymm2 \n\t"			/* r2 = r2 + (p - b) */ \
	"vmovdqa	%%ymm1, %%ymm4 \n\t"			/* copy root1s */ \
	"vmovdqa	%%ymm2, %%ymm5 \n\t"			/* copy root2s */ \
	"vpmaxuw	%%ymm4, %%ymm2, %%ymm2 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	"vpminuw	%%ymm5, %%ymm1, %%ymm1\n\t"			/* replace xmm1 with min of root1 and root2 */ \
	"vmovdqa %%ymm1, (%2) \n\t"				/* write new root1's */ \
	"vmovdqa %%ymm2, (%3) \n\t"				/* write new root2's */

#define _AVX2_SMALL_PRIME_SIEVE \
	for (; i < med_B; i += 16) \
		asm (			\
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$15, %%rsi \n\t"				/* logp's range from 15 to 16... call 'em = 15 */ \
			"vmovdqa	(%1), %%ymm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%ymm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%ymm2 \n\t"				/* bring in 8 root2's */ \
			\
			_16P_FINAL_STEP_SIEVE_AVX2					\
			_FINALIZE_SORT_UPDATE_AVX2 \
			\
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx", "rdx", \
				"rsi", "rdi", "cc", "memory");

#elif defined(_MSC_VER)

#define _SSE2_AVX2_SMALL_PRIME_SIEVE_32k_DIV3 \
		for (; i < full_fb->fb_32k_div3-16; i += 8) { \
			logp = 13; \
			primes = _mm_load_si128((__m128i *)(fb->prime + i)); \
			root1s = _mm_load_si128((__m128i *)(fb->root1 + i)); \
			root2s = _mm_load_si128((__m128i *)(fb->root2 + i)); \
			_8P_STEP_SIEVE \
			_8P_STEP_SIEVE \
			_8P_STEP_SIEVE \
			_8P_FINAL_STEP_SIEVE	 \
			_FINALIZE_SORT_UPDATE_SSE2 \
		} \
	} while (0);


#define _SSE2_AVX2_SMALL_PRIME_SIEVE_14b \
		for (; i < full_fb->fb_14bit_B-16; i += 8) { \
			logp = 14; \
			primes = _mm_load_si128((__m128i *)(fb->prime + i)); \
			root1s = _mm_load_si128((__m128i *)(fb->root1 + i)); \
			root2s = _mm_load_si128((__m128i *)(fb->root2 + i)); \
			_8P_STEP_SIEVE \
			_8P_STEP_SIEVE \
			_8P_FINAL_STEP_SIEVE	 \
			_FINALIZE_SORT_UPDATE_SSE2 \
		} \
	} while (0);

#define _SSE2_AVX2_SMALL_PRIME_SIEVE_15b \
		for (; i < full_fb->fb_15bit_B-16; i += 8) { \
			logp = 15; \
			primes = _mm_load_si128((__m128i *)(fb->prime + i)); \
			root1s = _mm_load_si128((__m128i *)(fb->root1 + i)); \
			root2s = _mm_load_si128((__m128i *)(fb->root2 + i)); \
			_8P_STEP_SIEVE \
			_8P_FINAL_STEP_SIEVE	 \
			_FINALIZE_SORT_UPDATE_SSE2 \
		} \
	} while (0);


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


