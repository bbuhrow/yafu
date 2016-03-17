
#define BLOCK_INIT \
	memset(sieve,s_init,32768);

#if defined(GCC_ASM64X) || defined(__MINGW64__)

#define SIEVE_GT_BLOCKSIZE_ASM \
	__asm__ ( \
		"movq	%0,%%r12 \n\t"					/* move helperstruct into r12 */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* move startprime (i) into r8d */ \
		"movq	(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"movl   44(%%r12,1),%%r15d \n\t"		/* move med_B into r15d */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
  		"movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
		"cmpl   44(%%r12,1),%%r8d \n\t"				/* i >= med_B? */ \
		"jae    9f \n\t"						/* jump to exit if so */ \
		"movq   8(%%r12,1),%%rdx \n\t"			/* rdx holds prime pointer */ \
		"movl   %%r8d,%%r14d \n\t"				/* copy i to ecx */ \
		"movq   32(%%r12,1),%%rax \n\t"			/* rax holds logp pointer */ \
			/* ==================================================================== */ \
			/* = prime > blocksize sieving							              = */ \
			/* ==================================================================== */ \
		"cmpl   44(%%r12,1),%%r8d \n\t" \
		"jae    9f \n\t" \
		"4: \n\t"								/* top of > blocksize sieving loop */ \
		"movzwl (%%rdx,%%r8,2),%%r10d \n\t"		/* bring in prime */	 \
		"movzwl	(%%r13,%%r8,2),%%r9d \n\t"		/* bring in root1 */ \
  		"movzwl (%%r11,%%r8,2),%%edi \n\t"		/* bring in root2 */	 \
  		"movzbl (%%rax,%%r8,2),%%esi \n\t"		/* bring in logp */ \
 		"cmpl   $32767,%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
 		"ja     1f \n\t" \
 		"movl   %%r9d,%%r11d \n\t" \
 		"addl   %%r10d,%%r9d \n\t" \
 		"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 		"cmpl   $32767,%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
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
 		"cmpl   44(%%r12,1),%%r8d \n\t" \
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
		: "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "memory", "cc");
		

#define SIEVE_13b_ASM \
	__asm__ ( \
		"movq	%0,%%r12 \n\t"					/* move helperstruct into r12 */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* move startprime (i) into r8d */ \
		"movq	(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
  		"movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
		"cmpl   44(%%r12,1),%%r8d \n\t"				/* i >= med_B? */ \
		"jae    8f \n\t"						/* jump to exit if so */ \
		"7: \n\t"								/* start of 2x sieving loop */  \
		"movq   8(%%r12,1),%%rdx \n\t"			/* rdx holds prime pointer */ \
		"movl   %%r8d,%%r14d \n\t"				/* copy i to ecx */ \
		"movq   32(%%r12,1),%%rax \n\t"			/* rax holds logp pointer */ \
  		"movzwl	(%%r13,%%r8,2),%%r9d \n\t"		/* bring in root1 */ \
  		"movzwl (%%r11,%%r8,2),%%edi \n\t"		/* bring in root2 */	 \
  		"movzwl (%%rdx,%%r8,2),%%r10d \n\t"		/* bring in prime */	 \
  		"movzbl (%%rax,%%r8,2),%%esi \n\t"		/* bring in logp */ \
			/* ==================================================================== */ \
			/* = 2x sieving											              = */ \
			/* ==================================================================== */ \
  		"movl   $32768,%%r13d \n\t"				/* copy blocksize ; root1 ptr overwritten */	 \
  		"movl   %%r10d,%%edx \n\t"				/* copy prime to edx; prime ptr overwritten */ \
  		"subl   %%r10d,%%r13d	 \n\t"			/* stop = blocksize - prime */ \
  		"leaq   (%%rdx,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
  		"cmpl   %%r13d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
  		"jae    1f \n\t"						/* jump past loop if so */ \
  		"leal   (%%r10,%%r10,1),%%r11d \n\t"	/* 2x prime in r11; root2 prime overwritten */ \
		"0:  \n\t"								/* sieve to "stop"(r13d) */ \
  		"movl   %%r9d,%%edx \n\t" \
  		"movl   %%edi,%%eax \n\t"				/* logp pointer overwritten */ \  		
  		"subb   %%sil,(%%rbx,%%rdx,1)	 \n\t"	/* rbx holds sieve */ \
  		"addl   %%r11d,%%edi \n\t" \
  		"subb   %%sil,(%%rbx,%%rax,1) \n\t" \
  		"subb   %%sil,(%%rcx,%%rdx,1)	 \n\t"	/* rcx holds sieve2 */ \
        "addl   %%r11d,%%r9d \n\t" \
  		"subb   %%sil,(%%rcx,%%rax,1) \n\t" \
  		"cmpl   %%r13d,%%edi \n\t" \
  		"jb     0b \n\t"						/* repeat */ \
		"1:  \n\t" \
  		"cmpl   $32767,%%edi \n\t"				/* root2 >= blocksize? */ \
		"ja     1f  \n\t"						/* jump to extra root1 check */ \
				 \
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
  		"cmpl   $32767,%%r9d \n\t"				/* root1 >= blocksize? */ \
  		"ja     1f \n\t"						/* jump past extra root1 block if so */ \
  		"movl   %%r9d,%%r13d \n\t" \
  		"subb   %%sil,(%%rbx,%%r13,1) \n\t" \
  		"movl   %%edi,%%esi \n\t" \
  		"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
  		"movl   %%esi,%%r9d \n\t" \
		"1:  \n\t" \
  		"movq   16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
 		"leaq   0xffffffffffff8000(%%r9),%%r10 \n\t" \
		"incl   %%r8d \n\t" \
 		"leaq   0xffffffffffff8000(%%rdi),%%r9 \n\t" \
 		"movw   %%r10w,0x0(%%r13,%%r14,2) \n\t" /* write out new root1 */ \
		"movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
		"cmpl   44(%%r12,1),%%r8d \n\t" \
		"movw   %%r9w,(%%r11,%%r14,2) \n\t"		/* write out new root2 */ \
 		"jb     7b \n\t"						/* repeat 2x sieving loop */ \
		"8: \n\t"													\
		"movl	%%r8d, 40(%%r12,1) \n\t"		/* copy out final value of i */ \
		:																\
		: "g"(&asm_input)												\
		: "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "memory", "cc");


#define _8P_STEP_SIEVE \
	"pextrw	$0, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"pextrw	$0, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"pextrw	$1, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"pextrw	$1, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"pextrw	$2, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"pextrw	$2, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"pextrw	$3, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"pextrw	$3, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"pextrw	$4, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"pextrw	$4, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"pextrw	$5, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"pextrw	$5, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"pextrw	$6, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"pextrw	$6, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
	"pextrw	$7, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
	"pextrw	$7, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"paddw	%%xmm0, %%xmm1 \n\t"			/* increment root1's by primes */ \
	"paddw	%%xmm0, %%xmm2 \n\t"			/* increment root2's by primes */

	// alternate formulation
	//"pextrw	$0, %%xmm1, %%eax \n\t"			/* extract root1 */
	//"pextrw	$0, %%xmm2, %%ebx \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */
	//"pextrw	$1, %%xmm1, %%ecx \n\t"			/* extract root1 */
	//"pextrw	$1, %%xmm2, %%edi \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */
	//"pextrw	$2, %%xmm1, %%eax \n\t"			/* extract root1 */
	//"pextrw	$2, %%xmm2, %%ebx \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */
	//"pextrw	$3, %%xmm1, %%ecx \n\t"			/* extract root1 */
	//"pextrw	$3, %%xmm2, %%edi \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */
	//"pextrw	$4, %%xmm1, %%eax \n\t"			/* extract root1 */
	//"pextrw	$4, %%xmm2, %%ebx \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */
	//"pextrw	$5, %%xmm1, %%ecx \n\t"			/* extract root1 */
	//"pextrw	$5, %%xmm2, %%edi \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */
	//"pextrw	$6, %%xmm1, %%eax \n\t"			/* extract root1 */
	//"pextrw	$6, %%xmm2, %%ebx \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */
	//"pextrw	$7, %%xmm1, %%ecx \n\t"			/* extract root1 */
	//"pextrw	$7, %%xmm2, %%edi \n\t"			/* extract root2 */
	//"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */
	//"subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */

#define _8P_FINAL_STEP_SIEVE \
	"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1's before comparison */ \
	"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root1's before comparison */ \
	/* workaround for lack of unsigned greater than... */ \
	"psubusw	%%xmm3, %%xmm4 \n\t"		/* xmm4 := root1 - 32767 */ \
	"pcmpeqw	%%xmm6, %%xmm4 \n\t"		/* xmm4 := root1 <= 32767 ? 1 : 0 */ \
	"psubusw	%%xmm3, %%xmm5 \n\t"		/* xmm5 := root2 - 32767 */ \
	"pcmpeqw	%%xmm6, %%xmm5 \n\t"		/* xmm5 := root2 <= 32767 ? 1 : 0 */ \
	"pmovmskb	%%xmm4, %%ecx \n\t" \
	"pmovmskb	%%xmm5, %%edi \n\t" \
	"testl	$0x2, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* if bits are equal to zero then root1 ! <= blocksize-1; jump */ \
	"pextrw	$0, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$0, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$1, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$1, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x20, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"pand	%%xmm0, %%xmm4 \n\t"			/* clear primes whose root1's are >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$2, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x20, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$2, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x80, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"pand	%%xmm0, %%xmm5 \n\t"			/* clear primes whose root2's are >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$3, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x80, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$3, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x200, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"psubw	%%xmm7, %%xmm4 \n\t"			/* advance root1's still below blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$4, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x200, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$4, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x800, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"psubw	%%xmm7, %%xmm5 \n\t"			/* advance root1's still below blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$5, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x800, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$5, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x2000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$6, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$6, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$7, %%xmm1, %%eax \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"pextrw	$7, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */

#define _FINALIZE_SORT_UPDATE_SSE2 \
	"paddw	%%xmm4, %%xmm1 \n\t"			/* r1 = r1 + (p - b) */ \
	"paddw	%%xmm5, %%xmm2 \n\t"			/* r2 = r2 + (p - b) */ \
	"movdqa	%%xmm1, %%xmm4 \n\t"			/* copy root1s */ \
	"movdqa	%%xmm2, %%xmm5 \n\t"			/* copy root2s */ \
	"pmaxsw	%%xmm4, %%xmm2 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	"pminsw	%%xmm5, %%xmm1 \n\t"			/* replace xmm1 with min of root1 and root2 */ \
	"movdqa %%xmm1, (%2) \n\t"				/* write new root1's */ \
	"movdqa %%xmm2, (%3) \n\t"				/* write new root2's */

#define _INIT_SSE2_SMALL_PRIME_SIEVE \
	asm (	\
		"movl	$0x7fff7fff, %%ecx \n\t"		/* load 2 copies of blocksize-1 in a 32bit reg */ \
		"movl	$0x80008000, %%edi \n\t"		/* load 2 copies of blocksize in a 32bit reg */ \
		"movd	%%ecx, %%xmm3 \n\t"				/* we need 32k-1 because we do signed comparisons */ \
		"movd	%%edi, %%xmm7 \n\t"				\
		"pxor	%%xmm6, %%xmm6 \n\t"			/* xmm6 := 0 */ \
		"pshufd	$0, %%xmm3, %%xmm3 \n\t"		/* broadcast blocksize-1 to all words of xmm3 */ \
		"pshufd	$0, %%xmm7, %%xmm7 \n\t"		/* broadcast blocksize to all words of xmm7 */ \
		: \
		:  \
		: "ecx", "edi", "xmm3", "xmm6", "xmm7");

#define _SSE2_SMALL_PRIME_SIEVE_32k_DIV3 \
	for (; i < full_fb->fb_32k_div3-8; i += 8) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$13, %%rsi \n\t"				/* logps = 13 in this range */ \
			"movdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"movdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"movdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
			 \
			_8P_STEP_SIEVE \
			_8P_STEP_SIEVE \
			_8P_STEP_SIEVE \
			_8P_FINAL_STEP_SIEVE	 \
			_FINALIZE_SORT_UPDATE_SSE2 \
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", \
				"rcx", "rdx", "rsi", "rdi", "cc", "memory");
	
#define _SSE2_SMALL_PRIME_SIEVE_14b \
	for (; i < full_fb->fb_14bit_B-8; i += 8) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$14, %%rsi \n\t"				/* logp's range from 13 to 14... call 'em = 14 */ \
			"movdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"movdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"movdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
			 \
			_8P_STEP_SIEVE \
			_8P_STEP_SIEVE \
			_8P_FINAL_STEP_SIEVE	\
			_FINALIZE_SORT_UPDATE_SSE2 \
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", \
				"rcx", "rdx", "rsi", "rdi", "cc", "memory");

#define _SSE2_SMALL_PRIME_SIEVE_15b \
	for (; i < full_fb->fb_15bit_B-8; i += 8) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$15, %%rsi \n\t"				/* logp's range from 14 to 15... call 'em = 15 */ \
			"movdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"movdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"movdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
				\
			_8P_STEP_SIEVE \
			_8P_FINAL_STEP_SIEVE			\
			_FINALIZE_SORT_UPDATE_SSE2 \
			\
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx", "rdx", "rsi", \
				"rdi", "cc", "memory");



#elif defined(_MSC_VER)


#define _8P_STEP_SIEVE \
	sieve[_mm_extract_epi16 (root1s, 0)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 0)] -= logp; \
	sieve[_mm_extract_epi16 (root1s, 1)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 1)] -= logp; \
	sieve[_mm_extract_epi16 (root1s, 2)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 2)] -= logp; \
	sieve[_mm_extract_epi16 (root1s, 3)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 3)] -= logp; \
	sieve[_mm_extract_epi16 (root1s, 4)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 4)] -= logp; \
	sieve[_mm_extract_epi16 (root1s, 5)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 5)] -= logp; \
	sieve[_mm_extract_epi16 (root1s, 6)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 6)] -= logp; \
	sieve[_mm_extract_epi16 (root1s, 7)] -= logp; \
	sieve[_mm_extract_epi16 (root2s, 7)] -= logp; \
	root1s = _mm_add_epi16(root1s, primes); \
	root2s = _mm_add_epi16(root2s, primes);


#define _8P_FINAL_STEP_SIEVE \
	tmp1 = root1s; \
	tmp2 = root2s; \
	tmp1 = _mm_subs_epu16 (tmp1, blocksize_m1); \
	tmp2 = _mm_subs_epu16 (tmp2, blocksize_m1); \
	tmp1 = _mm_cmpeq_epi16(tmp1, zeros); \
	tmp2 = _mm_cmpeq_epi16(tmp2, zeros); \
	root1mask = _mm_movemask_epi8(tmp1); \
	root2mask = _mm_movemask_epi8(tmp2); \
	if (root1mask & 0x2) { \
		sieve[_mm_extract_epi16 (root1s, 0)] -= logp; \
		if (root2mask & 0x2) \
			sieve[_mm_extract_epi16 (root2s, 0)] -= logp; \
	} \
	if (root1mask & 0x8) { \
		sieve[_mm_extract_epi16 (root1s, 1)] -= logp; \
		if (root2mask & 0x8) \
			sieve[_mm_extract_epi16 (root2s, 1)] -= logp; \
	} \
	tmp1 = _mm_and_si128(tmp1, primes); \
	if (root1mask & 0x20) { \
		sieve[_mm_extract_epi16 (root1s, 2)] -= logp; \
		if (root2mask & 0x20) \
			sieve[_mm_extract_epi16 (root2s, 2)] -= logp; \
	} \
	tmp2 = _mm_and_si128(tmp2, primes); \
	if (root1mask & 0x80) { \
		sieve[_mm_extract_epi16 (root1s, 3)] -= logp; \
		if (root2mask & 0x80) \
			sieve[_mm_extract_epi16 (root2s, 3)] -= logp; \
	} \
	tmp1 = _mm_sub_epi16(tmp1, blocksize); \
	if (root1mask & 0x200) { \
		sieve[_mm_extract_epi16 (root1s, 4)] -= logp; \
		if (root2mask & 0x200) \
			sieve[_mm_extract_epi16 (root2s, 4)] -= logp; \
	} \
	tmp2 = _mm_sub_epi16(tmp2, blocksize); \
	if (root1mask & 0x800) { \
		sieve[_mm_extract_epi16 (root1s, 5)] -= logp; \
		if (root2mask & 0x800) \
			sieve[_mm_extract_epi16 (root2s, 5)] -= logp; \
	} \
	if (root1mask & 0x2000) { \
		sieve[_mm_extract_epi16 (root1s, 6)] -= logp; \
		if (root2mask & 0x2000) \
			sieve[_mm_extract_epi16 (root2s, 6)] -= logp; \
	} \
	if (root1mask & 0x8000) { \
		sieve[_mm_extract_epi16 (root1s, 7)] -= logp; \
		if (root2mask & 0x8000) \
			sieve[_mm_extract_epi16 (root2s, 7)] -= logp; \
	}

#define _FINALIZE_SORT_UPDATE_SSE2 \
	root1s = _mm_add_epi16(root1s, tmp1); \
	root2s = _mm_add_epi16(root2s, tmp2); \
	tmp1 = root1s; \
	tmp2 = root2s; \
	root2s = _mm_max_epi16(root2s, tmp1); \
	root1s = _mm_min_epi16(root1s, tmp2); \
	_mm_store_si128((__m128i *)(fb->root1 + i), root1s); \
	_mm_store_si128((__m128i *)(fb->root2 + i), root2s);

#define _INIT_SSE2_SMALL_PRIME_SIEVE \
	do { \
		__m128i blocksize_m1 = _mm_set_epi32(0x7fff7fff, 0x7fff7fff, 0x7fff7fff, 0x7fff7fff); \
		__m128i blocksize = _mm_set_epi32(0x80008000, 0x80008000, 0x80008000, 0x80008000); \
		__m128i zeros = _mm_set_epi32(0, 0, 0, 0); \
		__m128i tmp1; \
		__m128i tmp2; \
		__m128i primes; \
		__m128i root1s; \
		__m128i root2s; \
		uint32 root1mask; \
		uint32 root2mask; \
		uint8 logp; 

#define _SSE2_SMALL_PRIME_SIEVE_32k_DIV3 \
		for (; i < full_fb->fb_32k_div3-8; i += 8) { \
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


#define _SSE2_SMALL_PRIME_SIEVE_14b \
		for (; i < full_fb->fb_14bit_B-8; i += 8) { \
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

#define _SSE2_SMALL_PRIME_SIEVE_15b \
		for (; i < full_fb->fb_15bit_B-8; i += 8) { \
			logp = 15; \
			primes = _mm_load_si128((__m128i *)(fb->prime + i)); \
			root1s = _mm_load_si128((__m128i *)(fb->root1 + i)); \
			root2s = _mm_load_si128((__m128i *)(fb->root2 + i)); \
			_8P_STEP_SIEVE \
			_8P_FINAL_STEP_SIEVE	 \
			_FINALIZE_SORT_UPDATE_SSE2 \
		} \
	} while (0);

#endif


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


