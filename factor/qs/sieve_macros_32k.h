
#define BLOCK_INIT \
	memset(sieve,s_init,32768);

#ifdef ASM_SIEVING

#define SIEVE_MED_P_ASM \
	__asm__ ( \
		"movq	%0,%%r12 \n\t"					/* move helperstruct into rsi */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* move startprime (i) into r8d */ \
		"movq	(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"movl   44(%%r12,1),%%r15d \n\t"		/* move med_B into r15d */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
  		"movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
		"cmpl   %%r15d,%%r8d \n\t"				/* i >= med_B? */ \
		"jae    9f \n\t"						/* jump to exit if so */ \
		"7: \n\t"								/* start of 2x sieving loop */  \
		"movq   8(%%r12,1),%%rdx \n\t"			/* rdx holds prime pointer */ \
		"movl   %%r8d,%%r14d \n\t"				/* copy i to ecx */ \
		"movq   32(%%r12,1),%%rax \n\t"			/* rax holds logp pointer */ \
  		"movzwl	(%%r13,%%r8,2),%%r9d \n\t"		/* bring in root1 */ \
  		"movzwl (%%r11,%%r8,2),%%edi \n\t"		/* bring in root2 */	 \
  		"movzwl (%%rdx,%%r8,2),%%r10d \n\t"		/* bring in prime */	 \
  		"movzbl (%%rax,%%r8,2),%%esi \n\t"		/* bring in logp */ \
  		"cmpl   $0x4000,%%r10d \n\t"			/* prime > 16384 */ \
  		"ja     8f \n\t"						/* jump to 1x sieving if so */ \
			/* ==================================================================== */ \
			/* = 2x sieving											              = */ \
			/* ==================================================================== */ \
  		"movl   $32768,%%r13d \n\t"			/* copy blocksize ; root1 ptr overwritten */	 \
  		"movl   %%r10d,%%edx \n\t"				/* copy prime to edx; prime ptr overwritten */ \
  		"subl   %%r10d,%%r13d	 \n\t"			/* stop = blocksize - prime */ \
  		"leaq   (%%rdx,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
  		"cmpl   %%r13d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
  		"jae    1f \n\t"						/* jump past loop if so */ \
  		"leal   (%%r10,%%r10,1),%%r11d \n\t"	/* 2x prime in r11; root2 prime overwritten */ \
		"0:  \n\t"								/* sieve to "stop"(r13d) */ \
  		"movl   %%r9d,%%edx \n\t" \
  		"movl   %%edi,%%eax \n\t"				/* logp pointer overwritten */ \
  		"addl   %%r11d,%%edi \n\t" \
  		"subb   %%sil,(%%rbx,%%rdx,1)	 \n\t"	/* rbx holds sieve */ \
  		"addl   %%r11d,%%r9d \n\t" \
  		"subb   %%sil,(%%rbx,%%rax,1) \n\t" \
  		"subb   %%sil,(%%rcx,%%rdx,1)	 \n\t"	/* rcx holds sieve2 */ \
  		"subb   %%sil,(%%rcx,%%rax,1) \n\t" \
  		"cmpl   %%r13d,%%edi \n\t" \
  		"jb     0b \n\t"						/* repeat */ \
		"1:  \n\t" \
  		"cmpl   $32767,%%edi \n\t"	/* root2 >= blocksize? */ \
		"ja     1f  \n\t"						/* jump to extra root1 check */ \
				/* data16 and nop */  \
		"0:  \n\t"								/* sieve to "blocksize" */ \
  		"movl   %%r9d,%%ecx \n\t" \
  		"movl   %%edi,%%edx \n\t" \
  		"addl   %%r10d,%%edi \n\t" \
  		"subb   %%sil,(%%rbx,%%rcx,1) \n\t" \
  		"addl   %%r10d,%%r9d \n\t" \
  		"subb   %%sil,(%%rbx,%%rdx,1) \n\t" \
  		"cmpl   $32767,%%edi \n\t" \
  		"jbe    0b \n\t"						/* repeat */ \
		"1:  \n\t" \
  		"cmpl   $32767,%%r9d \n\t"	/* root1 >= blocksize? */ \
  		"ja     1f \n\t"						/* jump past extra root1 block if so */ \
  		"movl   %%r9d,%%r13d \n\t" \
  		"subb   %%sil,(%%rbx,%%r13,1) \n\t" \
  		"movl   %%edi,%%esi \n\t" \
  		"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
  		"movl   %%esi,%%r9d \n\t" \
		"1:  \n\t" \
  		"movq   16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
  		"movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
  		"incl   %%r8d \n\t" \
 		"leaq   0xffffffffffff8000(%%r9),%%r10 \n\t" \
 		"cmpl   %%r15d,%%r8d \n\t" \
 		"leaq   0xffffffffffff8000(%%rdi),%%r9 \n\t" \
 		"movw   %%r10w,0x0(%%r13,%%r14,2) \n\t" /* write out new root1 */ \
 		"movw   %%r9w,(%%r11,%%r14,2) \n\t"		/* write out new root2 */ \
 		"jb     7b \n\t"						/* repeat 2x sieving loop */ \
		"cmpl   %%r15d,%%r8d \n\t" \
 		"jae    9f	 \n\t"						/* exit to lp sieving */ \
		"movq   8(%%r12,1),%%rdx \n\t"			/* rdx holds prime pointer */ \
		"movq   32(%%r12,1),%%rax \n\t"			/* rax holds logp pointer */ \
		"8: \n\t"													\
			/* ==================================================================== */ \
			/* = 1x sieving											              = */ \
			/* ==================================================================== */ \
		"cmpl	%%r15d,%%r8d \n\t"										\
		"jae    9f	 \n\t"						/* exit to lp sieving */ \
		"3: \n\t"								/* top of 1x sieving loop */ \
  		"movzwl (%%rdx,%%r8,2),%%r10d \n\t"		/* bring in prime */	 \
		"movzwl	(%%r13,%%r8,2),%%r9d \n\t"		/* bring in root1 */ \
  		"movzwl (%%r11,%%r8,2),%%edi \n\t"		/* bring in root2 */	 \
  		"movzbl (%%rax,%%r8,2),%%esi \n\t"		/* bring in logp */ \
 		"cmpl    $32768,%%r10d \n\t"										\
 		"ja     2f \n\t"						/* if prime > blocksize, exit loop */ \
 		"cmpl    $32767,%%edi \n\t"	/* if root2 > blocksize, skip to extra root1 check */ \
 		"ja     1f \n\t"																\
		"0: \n\t"								/* top of 1x unrolled loop */ \
 		"movl   %%r9d,%%r13d \n\t"										\
 		"movl   %%edi,%%r11d \n\t"												\
 		"addl   %%r10d,%%edi \n\t"												\
 		"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t"							\
 		"addl   %%r10d,%%r9d \n\t"										\
 		"subb   %%sil,(%%r11,%%rbx,1) \n\t"								\
 		"cmpl   $32767,%%edi \n\t"							\
 		"jbe    0b \n\t"						/* back to top of 1x unrolled loop */ \
		"1:  \n\t"																\
 		"cmpl   $32767,%%r9d \n\t"				/* extra root1 check */ \
 		"ja     1f \n\t"						/* if root1 > blocksize, skip past extra root1 check */ \
 		"movl   %%r9d,%%r13d \n\t"										\
 		"movl   %%edi,%%edx \n\t"										\
 		"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
 		"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t"							\
 		"movl   %%edx,%%r9d \n\t"				/* root1 = old root2 (swap) */ \
		"1:  \n\t"								/* end of extra root1 check */ \
		"movq   16(%%r12,1),%%r13 \n\t"									\
 		"movq   24(%%r12,1),%%r11 \n\t"									\
		"incl   %%r8d \n\t"												\
 		"leaq   0xffffffffffff8000(%%r9),%%r10 \n\t"							\
 		"leaq   0xffffffffffff8000(%%rdi),%%r9 \n\t"						\
 		"cmpl   %%r15d,%%r8d \n\t"										\
 		"movw   %%r10w,0x0(%%r13,%%r14,2) \n\t"							\
 		"movw   %%r9w,(%%r11,%%r14,2) \n\t"									\
 		"jae    9f \n\t"												\
 		"movq   8(%%r12,1),%%rdx \n\t"									\
 		"movq   32(%%r12,1),%%rax \n\t"									\
		"movl   %%r8d,%%r14d \n\t"				/* copy i to r14 */ \
 		"jmp	3b \n\t"												\
		"2:  \n\t"														\
			/* ==================================================================== */ \
			/* = prime > blocksize sieving							              = */ \
			/* ==================================================================== */ \
		"cmpl   %%r15d,%%r8d \n\t" \
		"jae    9f \n\t" \
		"4: \n\t"								/* top of > blocksize sieving loop */ \
		"movzwl (%%rdx,%%r8,2),%%r10d \n\t"		/* bring in prime */	 \
		"movzwl	(%%r13,%%r8,2),%%r9d \n\t"		/* bring in root1 */ \
  		"movzwl (%%r11,%%r8,2),%%edi \n\t"		/* bring in root2 */	 \
  		"movzbl (%%rax,%%r8,2),%%esi \n\t"		/* bring in logp */ \
 		"cmpl   $32767,%%r9d \n\t"	/* root1 > blocksize, skip to root update */ \
 		"ja     1f \n\t" \
 		"movl   %%r9d,%%r11d \n\t" \
 		"addl   %%r10d,%%r9d \n\t" \
 		"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 		"cmpl   $32767,%%edi \n\t"	/* root2 > blocksize, skip to swap roots */ \
 		"ja     2f \n\t" \
 		"movl   %%edi,%%r13d \n\t" \
 		"addl   %%r10d,%%edi \n\t" \
 		"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
		"3: \n\t" \
 		"movq   16(%%r12,1),%%r13 \n\t" \
 		"movq   24(%%r12,1),%%r11 \n\t" \
		"1:  \n\t"								/* update roots */ \
 		"incl   %%r8d \n\t" \
 		"leaq   0xffffffffffff8000(%%r9),%%rsi \n\t" \
 		"leaq   0xffffffffffff8000(%%rdi),%%r9 \n\t" \
 		"cmpl   %%r15d,%%r8d \n\t" \
 		"movw   %%si,0x0(%%r13,%%r14,2) \n\t" \
 		"movw   %%r9w,(%%r11,%%r14,2) \n\t" \
 		"jae    9f \n\t" \
 		"movl   %%r8d,%%r14d \n\t" \
 		"movq   8(%%r12,1),%%rdx \n\t" \
 		"movq   32(%%r12,1),%%rax \n\t" \
 		"jmp    4b \n\t"						/* back up to top of bigprime loop */ \
		"2: \n\t"								/* swap roots */ \
 		"movl   %%edi,%%edx \n\t" \
 		"movl   %%r9d,%%edi \n\t" \
 		"movl   %%edx,%%r9d \n\t" \
 		"jmp    3b \n\t"						/* jump to update roots */ \
 		"9: \n\t"								/* exit flag */ \
		"movl	%%r8d, 40(%%r12,1) \n\t"		/* copy out final value of i */ \
		:																\
		: "g"(&asm_input)												\
		: "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "memory", "cc");


#else

#define CHECK_2X_DONE \
	if (prime > 16384) \
		break;

#define CHECK_1X_DONE \
	if (prime > 32768) \
		break;

#define UPDATE_ROOTS \
	fb->root1[i] = (uint16)(root1 - 32768); \
	fb->root2[i] = (uint16)(root2 - 32768);

#define SIEVE_2X \
	stop = 32768 - prime;	\
	s2 = sieve + prime;	\
	while (root2 < stop) { \
		sieve[root1] -= logp; \
		sieve[root2] -= logp; \
		s2[root1] -= logp; \
		s2[root2] -= logp; \
		root1 += (prime << 1); \
		root2 += (prime << 1); \
	}
	
#define SIEVE_1X \
	while (root2 < 32768) { \
		sieve[root1] -= logp; \
		sieve[root2] -= logp; \
		root1 += prime; \
		root2 += prime; \
	}

#define SIEVE_LAST \
	if (root1 < 32768) { \
		sieve[root1] -= logp; \
		root1 += prime; \
		tmp = root2; \
		root2 = root1; \
		root1 = tmp; \
	}

#define SIEVE_BIG \
	if (root1 < 32768) { \
		sieve[root1] -= logp; \
		root1 += prime; \
		if (root2 < 32768) { \
			sieve[root2] -= logp; \
			root2 += prime; \
		} \
		else { \
			tmp=root2; \
			root2=root1; \
			root1=tmp; \
		} \
	}
	

#endif

