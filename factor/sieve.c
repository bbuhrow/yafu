/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"

#if defined(_MSC_VER)
	#include <mmintrin.h>
#endif

typedef struct
{
	uint8 *sieve;					//0
	uint16 *primeptr;				//8
	uint16 *root1ptr;				//16
	uint16 *root2ptr;				//24
	uint16 *logptr;					//32
	uint32 startprime;				//40
	uint32 med_B;					//44
} helperstruct_t;

//#undef SSE2_ASM_SIEVING
//#define ASM_SIEVING 1
//#define SSE2_ASM_SIEVING 1

#if defined(SSE2_ASM_SIEVING)

	#define SIEVE_ONE_PRIME_X2_LOOP(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"				/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"				/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"				/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"				/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r10d \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
  		"movl   $" BLOCKSIZEtxt ",%%r13d \n\t"	/* copy blocksize ; root1 ptr overwritten */	 \
  		"movl   %%r10d,%%edx \n\t"				/* copy prime to edx; prime ptr overwritten */ \
		"psrldq	$2,%%xmm0 \n\t" \
  		"subl   %%r10d,%%r13d	 \n\t"			/* stop = blocksize - prime */ \
  		"leaq   (%%rdx,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
  		"cmpl   %%r13d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
		"psrldq	$2,%%xmm1 \n\t" \
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
  		"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"	/* root2 >= blocksize? */ \
		"ja     1f  \n\t"						/* jump to extra root1 check */ \
		"0:  \n\t"								/* sieve to "blocksize" */ \
  		"movl   %%r9d,%%ecx \n\t" \
  		"movl   %%edi,%%edx \n\t" \
  		"addl   %%r10d,%%edi \n\t" \
  		"subb   %%sil,(%%rbx,%%rcx,1) \n\t" \
  		"addl   %%r10d,%%r9d \n\t" \
  		"subb   %%sil,(%%rbx,%%rdx,1) \n\t" \
  		"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
  		"jbe    0b \n\t"						/* repeat */ \
		"1:  \n\t" \
  		"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"	/* root1 >= blocksize? */ \
  		"ja     1f \n\t"						/* jump past extra root1 block if so */ \
  		"movl   %%r9d,%%r13d \n\t" \
  		"subb   %%sil,(%%rbx,%%r13,1) \n\t" \
  		"movl   %%edi,%%esi \n\t" \
  		"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
  		"movl   %%esi,%%r9d \n\t" \
		"1:  \n\t" \
		"psrldq	$2,%%xmm2 \n\t" \
 		"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
		"psrldq	$2,%%xmm3 \n\t" \
 		"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			/* ==================================================================== */ \
			/* = put new roots, in r10 and r9, back into xmm registers            = */ \
			/* ==================================================================== */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X2_LOOP_1(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r10d \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
  			"movl   $" BLOCKSIZEtxt ",%%r13d \n\t"			/* copy blocksize ; root1 ptr overwritten */	 \
  			"movl   %%r10d,%%edx \n\t"				/* copy prime to edx; prime ptr overwritten */ \
			"psrldq	$2,%%xmm0 \n\t" \
  			"subl   %%r10d,%%r13d	 \n\t"			/* stop = blocksize - prime */ \
  			"leaq   (%%rdx,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
  			"cmpl   %%r13d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
			"psrldq	$2,%%xmm1 \n\t" \
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
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 >= blocksize? */ \
			"ja     1f  \n\t"						/* jump to extra root1 check */ \
			"0:  \n\t"								/* sieve to "blocksize" */ \
  			"movl   %%r9d,%%ecx \n\t" \
  			"movl   %%edi,%%edx \n\t" \
  			"addl   %%r10d,%%edi \n\t" \
  			"subb   %%sil,(%%rbx,%%rcx,1) \n\t" \
  			"addl   %%r10d,%%r9d \n\t" \
  			"subb   %%sil,(%%rbx,%%rdx,1) \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
  			"jbe    0b \n\t"						/* repeat */ \
			"1:  \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 >= blocksize? */ \
  			"ja     1f \n\t"						/* jump past extra root1 block if so */ \
  			"movl   %%r9d,%%r13d \n\t" \
  			"subb   %%sil,(%%rbx,%%r13,1) \n\t" \
  			"movl   %%edi,%%esi \n\t" \
  			"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
  			"movl   %%esi,%%r9d \n\t" \
			"1:  \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
		"movq	8(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm6 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X2_LOOP_2(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r10d \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
  			"movl   $" BLOCKSIZEtxt ",%%r13d \n\t"			/* copy blocksize ; root1 ptr overwritten */	 \
  			"movl   %%r10d,%%edx \n\t"				/* copy prime to edx; prime ptr overwritten */ \
			"psrldq	$2,%%xmm0 \n\t" \
  			"subl   %%r10d,%%r13d	 \n\t"			/* stop = blocksize - prime */ \
  			"leaq   (%%rdx,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
  			"cmpl   %%r13d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
			"psrldq	$2,%%xmm1 \n\t" \
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
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 >= blocksize? */ \
			"ja     1f  \n\t"						/* jump to extra root1 check */ \
			"0:  \n\t"								/* sieve to "blocksize" */ \
  			"movl   %%r9d,%%ecx \n\t" \
  			"movl   %%edi,%%edx \n\t" \
  			"addl   %%r10d,%%edi \n\t" \
  			"subb   %%sil,(%%rbx,%%rcx,1) \n\t" \
  			"addl   %%r10d,%%r9d \n\t" \
  			"subb   %%sil,(%%rbx,%%rdx,1) \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
  			"jbe    0b \n\t"						/* repeat */ \
			"1:  \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 >= blocksize? */ \
  			"ja     1f \n\t"						/* jump past extra root1 block if so */ \
  			"movl   %%r9d,%%r13d \n\t" \
  			"subb   %%sil,(%%rbx,%%r13,1) \n\t" \
  			"movl   %%edi,%%esi \n\t" \
  			"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
  			"movl   %%esi,%%r9d \n\t" \
			"1:  \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
		"movq	16(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm7 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X2_LOOP_3(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r10d \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
  			"movl   $" BLOCKSIZEtxt ",%%r13d \n\t"			/* copy blocksize ; root1 ptr overwritten */	 \
  			"movl   %%r10d,%%edx \n\t"				/* copy prime to edx; prime ptr overwritten */ \
			"psrldq	$2,%%xmm0 \n\t" \
  			"subl   %%r10d,%%r13d	 \n\t"			/* stop = blocksize - prime */ \
  			"leaq   (%%rdx,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
  			"cmpl   %%r13d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
			"psrldq	$2,%%xmm1 \n\t" \
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
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 >= blocksize? */ \
			"ja     1f  \n\t"						/* jump to extra root1 check */ \
			"0:  \n\t"								/* sieve to "blocksize" */ \
  			"movl   %%r9d,%%ecx \n\t" \
  			"movl   %%edi,%%edx \n\t" \
  			"addl   %%r10d,%%edi \n\t" \
  			"subb   %%sil,(%%rbx,%%rcx,1) \n\t" \
  			"addl   %%r10d,%%r9d \n\t" \
  			"subb   %%sil,(%%rbx,%%rdx,1) \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
  			"jbe    0b \n\t"						/* repeat */ \
			"1:  \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 >= blocksize? */ \
  			"ja     1f \n\t"						/* jump past extra root1 block if so */ \
  			"movl   %%r9d,%%r13d \n\t" \
  			"subb   %%sil,(%%rbx,%%r13,1) \n\t" \
  			"movl   %%edi,%%esi \n\t" \
  			"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
  			"movl   %%esi,%%r9d \n\t" \
			"1:  \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
		"movq	24(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm8 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X2_LOOP_4(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r10d \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
  			"movl   $" BLOCKSIZEtxt ",%%r13d \n\t"			/* copy blocksize ; root1 ptr overwritten */	 \
  			"movl   %%r10d,%%edx \n\t"				/* copy prime to edx; prime ptr overwritten */ \
			"psrldq	$2,%%xmm0 \n\t" \
  			"subl   %%r10d,%%r13d	 \n\t"			/* stop = blocksize - prime */ \
  			"leaq   (%%rdx,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
  			"cmpl   %%r13d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
			"psrldq	$2,%%xmm1 \n\t" \
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
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 >= blocksize? */ \
			"ja     1f  \n\t"						/* jump to extra root1 check */ \
			"0:  \n\t"								/* sieve to "blocksize" */ \
  			"movl   %%r9d,%%ecx \n\t" \
  			"movl   %%edi,%%edx \n\t" \
  			"addl   %%r10d,%%edi \n\t" \
  			"subb   %%sil,(%%rbx,%%rcx,1) \n\t" \
  			"addl   %%r10d,%%r9d \n\t" \
  			"subb   %%sil,(%%rbx,%%rdx,1) \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
  			"jbe    0b \n\t"						/* repeat */ \
			"1:  \n\t" \
  			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 >= blocksize? */ \
  			"ja     1f \n\t"						/* jump past extra root1 block if so */ \
  			"movl   %%r9d,%%r13d \n\t" \
  			"subb   %%sil,(%%rbx,%%r13,1) \n\t" \
  			"movl   %%edi,%%esi \n\t" \
  			"leal   (%%r10,%%r9,1),%%edi \n\t"		/* root2 = root1 + prime */ \
  			"movl   %%esi,%%r9d \n\t" \
			"1:  \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
				"movq	32(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm9 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X1_LOOP(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
			"cmpl    $" BLOCKSIZEm1txt ",%%edi \n\t"			/* if root2 > blocksize, skip to extra root1 check */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
			"0: \n\t"								/* top of 1x unrolled loop */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%r11d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
 			"jbe    0b \n\t"						/* back to top of 1x unrolled loop */ \
			"1:  \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* extra root1 check */ \
 			"ja     1f \n\t"						/* if root1 > blocksize, skip past extra root1 check */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%edx \n\t" \
 			"leal   (%%r10,%%r9,1),%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
			"1:  \n\t"								/* end of extra root1 check */ \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X1_LOOP_1(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
			"cmpl    $" BLOCKSIZEm1txt ",%%edi \n\t"			/* if root2 > blocksize, skip to extra root1 check */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
			"0: \n\t"								/* top of 1x unrolled loop */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%r11d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
 			"jbe    0b \n\t"						/* back to top of 1x unrolled loop */ \
			"1:  \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* extra root1 check */ \
 			"ja     1f \n\t"						/* if root1 > blocksize, skip past extra root1 check */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%edx \n\t" \
 			"leal   (%%r10,%%r9,1),%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
			"1:  \n\t"								/* end of extra root1 check */ \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
				"movq	8(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm6 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X1_LOOP_2(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
			"cmpl    $" BLOCKSIZEm1txt ",%%edi \n\t"			/* if root2 > blocksize, skip to extra root1 check */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
			"0: \n\t"								/* top of 1x unrolled loop */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%r11d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
 			"jbe    0b \n\t"						/* back to top of 1x unrolled loop */ \
			"1:  \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* extra root1 check */ \
 			"ja     1f \n\t"						/* if root1 > blocksize, skip past extra root1 check */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%edx \n\t" \
 			"leal   (%%r10,%%r9,1),%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
			"1:  \n\t"								/* end of extra root1 check */ \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
				"movq	16(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm7 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X1_LOOP_3(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
			"cmpl    $" BLOCKSIZEm1txt ",%%edi \n\t"			/* if root2 > blocksize, skip to extra root1 check */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
			"0: \n\t"								/* top of 1x unrolled loop */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%r11d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
 			"jbe    0b \n\t"						/* back to top of 1x unrolled loop */ \
			"1:  \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* extra root1 check */ \
 			"ja     1f \n\t"						/* if root1 > blocksize, skip past extra root1 check */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%edx \n\t" \
 			"leal   (%%r10,%%r9,1),%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
			"1:  \n\t"								/* end of extra root1 check */ \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
				"movq	24(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm8 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_X1_LOOP_4(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%edi \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%r9d \n\t" \
			"cmpl    $" BLOCKSIZEm1txt ",%%edi \n\t"			/* if root2 > blocksize, skip to extra root1 check */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
			"0: \n\t"								/* top of 1x unrolled loop */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%r11d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
 			"jbe    0b \n\t"						/* back to top of 1x unrolled loop */ \
			"1:  \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* extra root1 check */ \
 			"ja     1f \n\t"						/* if root1 > blocksize, skip past extra root1 check */ \
 			"movl   %%r9d,%%r13d \n\t" \
 			"movl   %%edi,%%edx \n\t" \
 			"leal   (%%r10,%%r9,1),%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
			"1:  \n\t"								/* end of extra root1 check */ \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
				"movq	32(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm9 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_TWO_PRIME_1  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r14 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r15 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"shll	$16,%%r10d \n\t" 			/* move 2nd word up */ \
			"shll	$16,%%r9d \n\t" 			/* move 2nd word up */ \
			"orl	%%r14d, %%r10d \n\t"		/* combine roots */ \
			"orl	%%r15d, %%r9d \n\t"			/* combine roots */ \
			"movq	8(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"movd	%%r10d,%%xmm11 \n\t"				/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm6 \n\t"	/* xmm0 = next 8 prime */ \
		"movd	%%r9d,%%xmm10 \n\t"				/* insert root1 */

	#define SIEVE_TWO_PRIME_2  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r14 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r15 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"shll	$16,%%r10d \n\t" 			/* move 2nd word up */ \
			"shll	$16,%%r9d \n\t" 			/* move 2nd word up */ \
			"orl	%%r14d, %%r10d \n\t"		/* combine roots */ \
			"orl	%%r15d, %%r9d \n\t"			/* combine roots */ \
			"movq	16(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"movd	%%r10d,%%xmm5 \n\t"				/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm7 \n\t"	/* xmm0 = next 8 prime */ \
		"movd	%%r9d,%%xmm4 \n\t"				/* insert root1 */

	#define SIEVE_TWO_PRIME_3  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"pslldq	$4,%%xmm4 \n\t" \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"pslldq	$4,%%xmm5 \n\t" \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"por	%%xmm4, %%xmm10 \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
		"por	%%xmm5, %%xmm11 \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r14 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r15 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"shll	$16,%%r10d \n\t" 			/* move 2nd word up */ \
			"shll	$16,%%r9d \n\t" 			/* move 2nd word up */ \
			"orl	%%r14d, %%r10d \n\t"		/* combine roots */ \
			"orl	%%r15d, %%r9d \n\t"			/* combine roots */ \
			"movq	24(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"movd	%%r10d,%%xmm5 \n\t"				/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm8 \n\t"	/* xmm0 = next 8 prime */ \
		"movd	%%r9d,%%xmm4 \n\t"				/* insert root1 */

	#define SIEVE_TWO_PRIME_4  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"pslldq	$8,%%xmm4 \n\t" \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"pslldq	$8,%%xmm5 \n\t" \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"por	%%xmm4, %%xmm10 \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
		"por	%%xmm5, %%xmm11 \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r14 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r15 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"shll	$16, %%r10d \n\t" 			/* move 2nd word up */ \
			"shll	$16, %%r9d \n\t" 			/* move 2nd word up */ \
			"orl	%%r14d, %%r10d \n\t"		/* combine roots */ \
			"orl	%%r15d, %%r9d \n\t"			/* combine roots */ \
			"movq	32(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"movd	%%r10d,%%xmm5 \n\t"				/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm9 \n\t"	/* xmm0 = next 8 prime */ \
		"movd	%%r9d,%%xmm4 \n\t"				/* insert root1 */


	#define SIEVE_ONE_PRIME(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_1(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"movq	8(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm6 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

	#define SIEVE_ONE_PRIME_2(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"movq	16(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm7 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */


	#define SIEVE_ONE_PRIME_3(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"movq	24(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm8 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */


	#define SIEVE_ONE_PRIME_4(x)  \
			/* ============================================================ */ \
			/* =            bring in info for this iteration              = */ \
			/* ============================================================ */ \
		"movd	%%xmm0,%%r10d \n\t"		/* extract prime,x from xmm0 */ \
		"movd	%%xmm3,%%esi \n\t"		/* extract logp,x from xmm0 */ \
		"movd	%%xmm2,%%edi \n\t"		/* extract root2,x from xmm2 */ \
		"movd	%%xmm1,%%r9d \n\t"		/* extract root1,x from xmm1 */ \
		"andl	$0xffff, %%r9d \n\t" \
		"andl	$0xffff, %%r10d \n\t" \
		"psrldq	$2,%%xmm0 \n\t" \
		"andl	$0xff, %%esi \n\t" \
		"andl	$0xffff, %%edi \n\t" \
			"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* root1 > blocksize, skip to root update */ \
			"psrldq	$2,%%xmm1 \n\t" \
 			"ja     1f \n\t" \
 			"movl   %%r9d,%%r11d \n\t" \
 			"addl   %%r10d,%%r9d \n\t" \
 			"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 			"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"				/* root2 > blocksize, skip to swap roots */ \
 			"ja     2f \n\t" \
 			"movl   %%edi,%%r13d \n\t" \
 			"addl   %%r10d,%%edi \n\t" \
 			"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
			"1:  \n\t"								/* update roots */ \
 			"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
			"psrldq	$2,%%xmm2 \n\t" \
 			"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
			"psrldq	$2,%%xmm3 \n\t" \
			"jmp	3f \n\t" \
			"2: \n\t"								/* swap roots */ \
 			"movl   %%edi,%%edx \n\t" \
 			"movl   %%r9d,%%edi \n\t" \
 			"movl   %%edx,%%r9d \n\t" \
 			"jmp    1b \n\t"						/* jump to update roots */ \
				/* ==================================================================== */ \
				/* = put new roots, in r10 and r9, back into xmm registers            = */ \
				/* ==================================================================== */ \
			"3: \n\t" \
			"movq	32(%%r12,1),%%r14 \n\t"			/* move prime pointer into rcx */ \
		"pinsrw	$" x ",%%r10d,%%xmm5 \n\t"		/* insert root2 */ \
		"movdqa 16(%%r14,%%r8,2), %%xmm9 \n\t"	/* xmm0 = next 8 prime */ \
		"pinsrw	$" x ",%%r9d,%%xmm4 \n\t"		/* insert root1 */

#endif

#if defined(_MSC_VER)
	#define ADDRESS_SIEVE(j) (sieve[j])
	#define ADDRESS_LOC(j) (bptr[j].loc)
#else
	#define ADDRESS_SIEVE(j) (*(sieve + (j)))
	#define ADDRESS_LOC(j) ((bptr + j)->loc)
#endif

void lp_sieveblock(uint8 *sieve, sieve_fb_compressed *fb, fb_list *full_fb, uint32 bnum, uint32 numblocks,
							 lp_bucket *lp, uint32 start_prime, uint8 s_init, int side)
{
	uint32 i,j,lpnum,basebucket;
	uint32 med_B, fb_14bit_B, fb_15bit_B;
	uint8 logp;
#if !defined(SSE2_ASM_SIEVING) && !defined(ASM_SIEVING)
	sieve_fb_compressed *fbptr;
	uint32 prime, root1, root2, tmp, stop, p16;
#endif
	uint32 *bptr;
	helperstruct_t asm_input;

	med_B = full_fb->med_B;
	fb_14bit_B = full_fb->fb_14bit_B;
	fb_15bit_B = full_fb->fb_15bit_B;
	
#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	//initialize the block
	memset(sieve,s_init,BLOCKSIZE);

#ifdef SSE2_ASM_SIEVING
	//this currently does work - but med_B is loaded with magic number for the test c75
	//would need to load in sconf to make it general purpose.  Also it could be optimized
	//by doing two primes at a time, combining the new roots from both primes in one dword,
	//and hence being able to remove the pinsrw's in favor of movd's and pslldq's

	stop = start_prime + (8 - start_prime % 8);

	//printf("starting x2 sieving at offset %u, prime %u\n",start_prime,fb->prime[start_prime]);

	//get to the first multiple of 8
	for (i=start_prime;i<stop;i++)
	{	
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];

		//if we are past the blocksize, bail out because there are faster methods
		if (prime > BLOCKSIZE)
			break;

		while (root2 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			ADDRESS_SIEVE(root2) -= logp;
			root1 += prime;
			root2 += prime;
		}

		//don't forget the last proot1[i], and compute the roots for the next block
		if (root1 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			root1 += prime;
			//root1 will be bigger on the next iteration, switch them now
			tmp = root2;
			root2 = root1;
			root1 = tmp;
		}
			
		fb->root1[i] = (uint16)(root1 - BLOCKSIZE);
		fb->root2[i] = (uint16)(root2 - BLOCKSIZE);
	}	

	//load up the helperstruct
	asm_input.logptr = fb->logp;
	asm_input.primeptr = fb->prime;
	asm_input.root1ptr = fb->root1;
	asm_input.root2ptr = fb->root2;
	asm_input.sieve = sieve;
	asm_input.startprime = i;
	asm_input.med_B = fb_14bit_B - (fb_14bit_B % 8);

	//printf("starting x2 sieving at offset %u, prime %u\n",i,fb->prime[i]);

	__asm__ (		\
		"movq	%0,%%r12 \n\t"					/* move helperstruct into rsi */ \
		"movq	0(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* start_prime = i = r15d */ \
		"movq	8(%%r12,1),%%rcx \n\t"			/* move prime pointer into rcx */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* move root1 pointer into r13 */ \
		"movq	24(%%r12,1),%%r14 \n\t"			/* move root2 pointer into r14 */ \
		"movq	32(%%r12,1),%%r11 \n\t"			/* move logp pointer into r11 */ \
		"movdqa (%%rcx,%%r8,2), %%xmm0 \n\t"	/* xmm0 = next 8 prime */ \
		"movdqa (%%r13,%%r8,2), %%xmm1 \n\t"	/* xmm1 = next 8 root1 start locations */ \
		"movdqa (%%r14,%%r8,2), %%xmm2 \n\t"	/* xmm2 = next 8 root2 start locations */ \
		"movdqa (%%r11,%%r8,2), %%xmm3 \n\t"	/* xmm3 = next 8 logp */ \
		"cmpl   44(%%r12,1),%%r8d \n\t"			/* i >= med_B ? */ \
		"jae    6f	\n\t"						/* jump to end of loop, if test fails */ \
		"5: \n\t"	\
		SIEVE_ONE_PRIME_X2_LOOP_1("0")	\
		SIEVE_ONE_PRIME_X2_LOOP("1")	\
		SIEVE_ONE_PRIME_X2_LOOP_2("2")	\
		SIEVE_ONE_PRIME_X2_LOOP("3")	\
		SIEVE_ONE_PRIME_X2_LOOP_3("4")	\
		SIEVE_ONE_PRIME_X2_LOOP("5")	\
		SIEVE_ONE_PRIME_X2_LOOP_4("6")	\
		SIEVE_ONE_PRIME_X2_LOOP("7")	 \
		"movq	16(%%r12,1),%%r13 \n\t"			/* move root1 pointer into r13 */ \
		"movq	24(%%r12,1),%%r14 \n\t"			/* move root2 pointer into r14 */ \
		"movdqa %%xmm5, (%%r13,%%r8,2) \n\t"	/* write out new root1s */ \
		"movdqa %%xmm4, (%%r14,%%r8,2) \n\t"	/* write out new root2s */ \
		"addl	$8, %%r8d \n\t"					/* increment by 8 */ \
		"movdqa %%xmm6, %%xmm0 \n\t"			/* copy over new primes */ \
		"movdqa %%xmm7, %%xmm1 \n\t"			/* copy over new root1s */ \
		"movdqa %%xmm8, %%xmm2 \n\t"			/* copy over new root2s */ \
		"movdqa %%xmm9, %%xmm3 \n\t"			/* copy over new logps */ \
		"cmpl	44(%%r12,1),%%r8d \n\t"			/* i < med_B ? */ \
		"jb		5b \n\t"						/* do next 8 primes */ \
		"6:		\n\t"							/* end of loop */ \
		"movl	%%r8d, 40(%%r12,1) \n\t"		/* copy out final value of i */ \
		:	\
		: "g"(&asm_input) \
		: "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "xmm1", "xmm2", "xmm0", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "memory", "cc");

	i = asm_input.startprime;

	//load up the helperstruct
	asm_input.logptr = fb->logp;
	asm_input.primeptr = fb->prime;
	asm_input.root1ptr = fb->root1;
	asm_input.root2ptr = fb->root2;
	asm_input.sieve = sieve;
	asm_input.startprime = i;
#ifdef YAFU_64K
	asm_input.med_B = med_B;
#else
	asm_input.med_B = fb_15bit_B + (8 - fb_15bit_B % 8);
#endif

	//printf("starting x1 sieving at offset %u, prime %u\n",i,fb->prime[i]);

	__asm__ (		\
		"movq	%0,%%r12 \n\t"					/* move helperstruct into rsi */ \
		"movq	0(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* start_prime = i = r15d */ \
		"cmpl   44(%%r12,1),%%r8d \n\t"			/* i >= med_B ? */ \
		"jae    6f	\n\t"						/* jump to end of loop, if test fails */ \
		"movq	8(%%r12,1),%%rcx \n\t"			/* move prime pointer into rcx */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* move root1 pointer into r13 */ \
		"movq	24(%%r12,1),%%r14 \n\t"			/* move root2 pointer into r14 */ \
		"movq	32(%%r12,1),%%r11 \n\t"			/* move logp pointer into r11 */ \
		"movdqa (%%rcx,%%r8,2), %%xmm0 \n\t"	/* xmm0 = next 8 prime */ \
		"movdqa (%%r13,%%r8,2), %%xmm1 \n\t"	/* xmm1 = next 8 root1 start locations */ \
		"movdqa (%%r14,%%r8,2), %%xmm2 \n\t"	/* xmm2 = next 8 root2 start locations */ \
		"movdqa (%%r11,%%r8,2), %%xmm3 \n\t"	/* xmm3 = next 8 logp */ \		
		"5: \n\t"	\
		SIEVE_ONE_PRIME_X1_LOOP_1("0")	\
		SIEVE_ONE_PRIME_X1_LOOP("1")	\
		SIEVE_ONE_PRIME_X1_LOOP_2("2")	\
		SIEVE_ONE_PRIME_X1_LOOP("3")	\
		SIEVE_ONE_PRIME_X1_LOOP_3("4")	\
		SIEVE_ONE_PRIME_X1_LOOP("5")	\
		SIEVE_ONE_PRIME_X1_LOOP_4("6")	\
		SIEVE_ONE_PRIME_X1_LOOP("7")	\
		"movq	16(%%r12,1),%%r13 \n\t"			/* move root1 pointer into r13 */ \
		"movq	24(%%r12,1),%%r14 \n\t"			/* move root2 pointer into r14 */ \
		"movdqa %%xmm5, (%%r13,%%r8,2) \n\t"	/* write out new root1s */ \
		"movdqa %%xmm4, (%%r14,%%r8,2) \n\t"	/* write out new root2s */ \
		"addl	$8, %%r8d \n\t"					/* increment by 8 */ \
		"movdqa %%xmm6, %%xmm0 \n\t"			/* copy over new primes */ \
		"movdqa %%xmm7, %%xmm1 \n\t"			/* copy over new root1s */ \
		"movdqa %%xmm8, %%xmm2 \n\t"			/* copy over new root2s */ \
		"movdqa %%xmm9, %%xmm3 \n\t"			/* copy over new logps */ \
		"cmpl	44(%%r12,1),%%r8d \n\t"			/* i < med_B ? */ \
		"jb		5b \n\t"						/* do next 8 primes */ \
		"6:		\n\t"							/* end of loop */ \
		"movl	%%r8d, 40(%%r12,1) \n\t"		/* copy out final value of i */ \
		:	\
		: "g"(&asm_input) \
		: "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "xmm1", "xmm2", "xmm0", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "memory", "cc");

	i = asm_input.startprime;

#if !defined(YAFU_64K)
	//load up the helperstruct
	asm_input.logptr = fb->logp;
	asm_input.primeptr = fb->prime;
	asm_input.root1ptr = fb->root1;
	asm_input.root2ptr = fb->root2;
	asm_input.sieve = sieve;
	asm_input.startprime = i;
	asm_input.med_B = med_B;

	__asm__ (		\
		"movq	%0,%%r12 \n\t"					/* move helperstruct into rsi */ \
		"movq	0(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* start_prime = i = r15d */ \
		"cmpl   44(%%r12,1),%%r8d \n\t"			/* i >= med_B ? */ \
		"jae    6f	\n\t"						/* jump to end of loop, if test fails */ \
		"movq	8(%%r12,1),%%rcx \n\t"			/* move prime pointer into rcx */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* move root1 pointer into r13 */ \
		"movq	24(%%r12,1),%%r14 \n\t"			/* move root2 pointer into r14 */ \
		"movq	32(%%r12,1),%%r11 \n\t"			/* move logp pointer into r11 */ \
		"movdqa (%%rcx,%%r8,2), %%xmm0 \n\t"	/* xmm0 = next 8 prime */ \
		"movdqa (%%r13,%%r8,2), %%xmm1 \n\t"	/* xmm1 = next 8 root1 start locations */ \
		"movdqa (%%r14,%%r8,2), %%xmm2 \n\t"	/* xmm2 = next 8 root2 start locations */ \
		"movdqa (%%r11,%%r8,2), %%xmm3 \n\t"	/* xmm3 = next 8 logp */ \
		"5: \n\t"	\
		SIEVE_ONE_PRIME_1("0")	\
		SIEVE_ONE_PRIME("1")	\
		SIEVE_ONE_PRIME_2("2")	\
		SIEVE_ONE_PRIME("3")	\
		SIEVE_ONE_PRIME_3("4")	\
		SIEVE_ONE_PRIME("5")	\
		SIEVE_ONE_PRIME_4("6")	\
		SIEVE_ONE_PRIME("7")	 \
		/*SIEVE_TWO_PRIME_1	\
		SIEVE_TWO_PRIME_2	\
		SIEVE_TWO_PRIME_3	\
		SIEVE_TWO_PRIME_4	\
		"pslldq	$12,%%xmm4 \n\t" \		
		"pslldq	$12,%%xmm5 \n\t" \
		"por	%%xmm4, %%xmm10 \n\t" \
		"por	%%xmm5, %%xmm11 \n\t" */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* move root1 pointer into r13 */ \
		"movq	24(%%r12,1),%%r14 \n\t"			/* move root2 pointer into r14 */ \
		"movdqa %%xmm5, (%%r13,%%r8,2) \n\t"	/* write out new root1s */ \
		"movdqa %%xmm4, (%%r14,%%r8,2) \n\t"	/* write out new root2s */ \
		"addl	$8, %%r8d \n\t"					/* increment by 8 */ \
		"movdqa %%xmm6, %%xmm0 \n\t"			/* copy over new primes */ \
		"movdqa %%xmm7, %%xmm1 \n\t"			/* copy over new root1s */ \
		"movdqa %%xmm8, %%xmm2 \n\t"			/* copy over new root2s */ \
		"movdqa %%xmm9, %%xmm3 \n\t"			/* copy over new logps */ \
		"cmpl	44(%%r12,1),%%r8d \n\t"			/* i < med_B ? */ \
		"jb		5b \n\t"						/* do next 8 primes */ \
		"6:		\n\t"							/* end of loop */ \
		"movl	%%r8d, 40(%%r12,1) \n\t"		/* copy out final value of i */ \
		:	\
		: "g"(&asm_input) \
		: "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "xmm1", "xmm2", "xmm0", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm8", "xmm9", "xmm10", "xmm11", "memory", "cc");

	i = asm_input.startprime;

#endif

#elif ASM_SIEVING

	asm_input.logptr = fb->logp;
	asm_input.primeptr = fb->prime;
	asm_input.root1ptr = fb->root1;
	asm_input.root2ptr = fb->root2;
	asm_input.sieve = sieve;
	asm_input.startprime = start_prime;
	asm_input.med_B = med_B;

	__asm__ ( \
		"movq	%0,%%r12 \n\t"					/* move helperstruct into rsi */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* move startprime (i) into r8d */ \
		"movq	(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"movl   44(%%r12,1),%%r15d \n\t"		/* move med_B into r15d */ \
		"movq	16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
  		"movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
		"cmpl   %%r15d,%%r8d \n\t"				/* i >= med_B? */ \
		"jae    9f \n\t"						/* jump to exit if so */ \
		"7: \n\t"								/* start of 2x sieving loop */  		
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
  		"movl   $" BLOCKSIZEtxt ",%%r13d \n\t"			/* copy blocksize ; root1 ptr overwritten */	 \
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
  		"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"	/* root2 >= blocksize? */ \
		"ja     1f  \n\t"						/* jump to extra root1 check */ \
				/* data16 and nop */  \
		"0:  \n\t"								/* sieve to "blocksize" */ \
  		"movl   %%r9d,%%ecx \n\t" \
  		"movl   %%edi,%%edx \n\t" \
  		"addl   %%r10d,%%edi \n\t" \
  		"subb   %%sil,(%%rbx,%%rcx,1) \n\t" \
  		"addl   %%r10d,%%r9d \n\t" \
  		"subb   %%sil,(%%rbx,%%rdx,1) \n\t" \
  		"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t" \
  		"jbe    0b \n\t"						/* repeat */ \
		"1:  \n\t" \
  		"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"	/* root1 >= blocksize? */ \
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
 		"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t" \
 		"cmpl   %%r15d,%%r8d \n\t" \
 		"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
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
 		"cmpl    $" BLOCKSIZEtxt ",%%r10d \n\t"										\
 		"ja     2f \n\t"						/* if prime > blocksize, exit loop */ \
 		"cmpl    $" BLOCKSIZEm1txt ",%%edi \n\t"	/* if root2 > blocksize, skip to extra root1 check */ \
 		"ja     1f \n\t"																\
		"0: \n\t"								/* top of 1x unrolled loop */ \
 		"movl   %%r9d,%%r13d \n\t"										\
 		"movl   %%edi,%%r11d \n\t"												\
 		"addl   %%r10d,%%edi \n\t"												\
 		"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t"							\
 		"addl   %%r10d,%%r9d \n\t"										\
 		"subb   %%sil,(%%r11,%%rbx,1) \n\t"								\
 		"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"							\
 		"jbe    0b \n\t"						/* back to top of 1x unrolled loop */ \
		"1:  \n\t"																\
 		"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"				/* extra root1 check */ \
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
 		"leaq   " negBLOCKSIZE "(%%r9),%%r10 \n\t"							\
 		"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t"						\
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
 		"cmpl   $" BLOCKSIZEm1txt ",%%r9d \n\t"	/* root1 > blocksize, skip to root update */ \
 		"ja     1f \n\t" \
 		"movl   %%r9d,%%r11d \n\t" \
 		"addl   %%r10d,%%r9d \n\t" \
 		"subb   %%sil,(%%r11,%%rbx,1) \n\t" \
 		"cmpl   $" BLOCKSIZEm1txt ",%%edi \n\t"	/* root2 > blocksize, skip to swap roots */ \
 		"ja     2f \n\t" \
 		"movl   %%edi,%%r13d \n\t" \
 		"addl   %%r10d,%%edi \n\t" \
 		"subb   %%sil,0x0(%%r13,%%rbx,1) \n\t" \
		"3: \n\t" \
 		"movq   16(%%r12,1),%%r13 \n\t" \
 		"movq   24(%%r12,1),%%r11 \n\t" \
		"1:  \n\t"								/* update roots */ \
 		"incl   %%r8d \n\t" \
 		"leaq   " negBLOCKSIZE "(%%r9),%%rsi \n\t" \
 		"leaq   " negBLOCKSIZE "(%%rdi),%%r9 \n\t" \
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
 		"jmp    3b \n\t"						/* jump to update roots */
 		"9: \n\t"								/* exit flag */ \
		"movl	%%r8d, 40(%%r12,1) \n\t"		/* copy out final value of i */ \
		:																\
		: "g"(&asm_input)												\
		: "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "memory", "cc");
	
	i = asm_input.startprime;

#else

	p16 = BLOCKSIZE >> 1;
	for (i=start_prime;i<med_B;i++)
	{	
		uint8 *s2;		

#ifdef USE_COMPRESSED_FB
		fbptr = fb + i;
		prime = fbptr->prime_and_logp & 0xFFFF;
		root1 = fbptr->roots & 0xFFFF;
		root2 = fbptr->roots >> 16;
		logp = fbptr->prime_and_logp >> 16;
#else
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];
#endif

		//if we are past the blocksize, bail out because there are faster methods
		if (prime > p16)
			break;

		stop = BLOCKSIZE - prime;
		s2 = sieve + prime;

		while (root2 < stop)
		{
			sieve[root1] -= logp;
			sieve[root2] -= logp;
			s2[root1] -= logp;
			s2[root2] -= logp;
			root1 += (prime << 1);
			root2 += (prime << 1);
		}

		while (root2 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			ADDRESS_SIEVE(root2) -= logp;
			root1 += prime;
			root2 += prime;
		}

		//don't forget the last proot1[i], and compute the roots for the next block
		if (root1 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			root1 += prime;
			//root1 will be bigger on the next iteration, switch them now
			tmp = root2;
			root2 = root1;
			root1 = tmp;
		}
			
#ifdef USE_COMPRESSED_FB
		fbptr->roots = (uint32)(((root2 - BLOCKSIZE) << 16) | (root1 - BLOCKSIZE));
#else
		fb->root1[i] = (uint16)(root1 - BLOCKSIZE);
		fb->root2[i] = (uint16)(root2 - BLOCKSIZE);
#endif

	}

	for (;i<med_B;i++)
	{	
#ifdef USE_COMPRESSED_FB
		fbptr = fb + i;
		prime = fbptr->prime_and_logp & 0xFFFF;
		root1 = fbptr->roots & 0xFFFF;
		root2 = fbptr->roots >> 16;
		logp = fbptr->prime_and_logp >> 16;
#else
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];
#endif

		//if we are past the blocksize, bail out because there are faster methods
		if (prime > BLOCKSIZE)
			break;

		while (root2 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			ADDRESS_SIEVE(root2) -= logp;
			root1 += prime;
			root2 += prime;
		}

		//don't forget the last proot1[i], and compute the roots for the next block
		if (root1 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			root1 += prime;
			//root1 will be bigger on the next iteration, switch them now
			tmp = root2;
			root2 = root1;
			root1 = tmp;
		}
			
#ifdef USE_COMPRESSED_FB
		fbptr->roots = (uint32)(((root2 - BLOCKSIZE) << 16) | (root1 - BLOCKSIZE));
#else
		fb->root1[i] = (uint16)(root1 - BLOCKSIZE);
		fb->root2[i] = (uint16)(root2 - BLOCKSIZE);
#endif
	}

	//if there are primes left bigger than the blocksize, this will take
	//care of them.  if not, it doesn't run at all.
	for (;i<med_B;i++)
	{	
#ifdef USE_COMPRESSED_FB
		fbptr = fb + i;
		prime = fbptr->prime_and_logp & 0xFFFF;
		root1 = fbptr->roots & 0xFFFF;
		root2 = fbptr->roots >> 16;
		logp = fbptr->prime_and_logp >> 16;
#else
		prime = fb->prime[i];
		root1 = fb->root1[i];
		root2 = fb->root2[i];
		logp = fb->logp[i];
#endif

		if (root1 < BLOCKSIZE)
		{
			ADDRESS_SIEVE(root1) -= logp;
			root1 += prime;

			if (root2 < BLOCKSIZE)
			{
				ADDRESS_SIEVE(root2) -= logp;
				root2 += prime;
			}
			else
			{
				tmp=root2;
				root2=root1;
				root1=tmp;
			}
		}

#ifdef USE_COMPRESSED_FB
		fbptr = fb + i;
		prime = fbptr->prime_and_logp & 0xFFFF;
		root1 = fbptr->roots & 0xFFFF;
		root2 = fbptr->roots >> 16;
		logp = fbptr->prime_and_logp >> 16;
#else
		fb->root1[i] = (uint16)(root1 - BLOCKSIZE);
		fb->root2[i] = (uint16)(root2 - BLOCKSIZE);
#endif

	}

#endif	


#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		SIEVE_STG1 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);

		gettimeofday(&qs_timing_start, NULL);
#endif


	//finally, dump the buckets into the now cached 
	//sieve block in prefetched batches
	//the starting address for this slice and bucket can be computed as an offset from
	//lp->list, the list of bucket hits.  each block is a bucket, and each bucket gets
	//2^BUCKET_BITS bytes of storage.  negative buckets are arranged immediately after
	//positive buckets, thus we have numblocks positive buckets followed by numblocks
	//negative buckets in every slice.  each slice is a complete set of such buckets.  slices
	//are contiguous in memory, so the whole thing is physically one giant linear list of
	//bucket hits which we have subdivided.  
	bptr = lp->list + (bnum << BUCKET_BITS);
	if (side)
	{
		bptr += (numblocks << BUCKET_BITS);
		basebucket = numblocks;
	}
	else
		basebucket = 0;

	//use x8 when cache line has 32 bytes
	//use x16 when chache line has 64 bytes
#if defined(CACHE_LINE_64)
	//printf("cache_line_64\n");
	for (j=0;j<lp->num_slices;j++)
	{
		lpnum = *(lp->num + bnum + basebucket);
		//printf("dumping %d primes from slice %d, bucket %d\n",lpnum, j, bnum);
		logp = *(lp->logp + j);
		for (i = 0; i < (lpnum & (uint32)(~15)); i += 16)
		{
#if defined(MANUAL_PREFETCH)
#if defined(GCC_ASM64X) || defined (__MINGW32__)
			//_mm_prefetch(bptr + i + 16,_MM_HINT_NTA);
			__asm ("prefetchnta %0	\n\t"
				: 
				:"m"(bptr + i + 16));
#else
			_mm_prefetch(bptr + i + 16,_MM_HINT_NTA);
#endif
#endif

			sieve[bptr[i  ] & 0x0000ffff] -= logp;
			sieve[bptr[i+1] & 0x0000ffff] -= logp;
			sieve[bptr[i+2] & 0x0000ffff] -= logp;
			sieve[bptr[i+3] & 0x0000ffff] -= logp;
			sieve[bptr[i+4] & 0x0000ffff] -= logp;
			sieve[bptr[i+5] & 0x0000ffff] -= logp;
			sieve[bptr[i+6] & 0x0000ffff] -= logp;
			sieve[bptr[i+7] & 0x0000ffff] -= logp;
			sieve[bptr[i+8] & 0x0000ffff] -= logp;
			sieve[bptr[i+9] & 0x0000ffff] -= logp;
			sieve[bptr[i+10] & 0x0000ffff] -= logp;
			sieve[bptr[i+11] & 0x0000ffff] -= logp;
			sieve[bptr[i+12] & 0x0000ffff] -= logp;
			sieve[bptr[i+13] & 0x0000ffff] -= logp;
			sieve[bptr[i+14] & 0x0000ffff] -= logp;
			sieve[bptr[i+15] & 0x0000ffff] -= logp;
		}

		for (;i<lpnum;i++)
			sieve[bptr[i] & 0x0000ffff] -= logp;
		
		//point to the next slice of primes
		bptr += (numblocks << (BUCKET_BITS + 1));
		basebucket += (numblocks << 1);
	}
#else

	for (j=0;j<lp->num_slices;j++)
	{
		lpnum = *(lp->num + bnum + basebucket);
		//printf("dumping %d primes from slice %d, bucket %d\n",lpnum, j, bnum);
		logp = *(lp->logp + j);
		for (i = 0; i < (lpnum & (uint32)(~7)); i += 8)
		//the slices can be considered stacks; the highest indices are put in last.
		//so it makes sense to take these out first, as they are more likely to 
		//still be in cache.
		//for (i = lpnum; i > (lpnum & 15); i -=16)
		{
			
#if defined(MANUAL_PREFETCH)
#if defined(GCC_ASM64X) || defined (__MINGW32__)
			//_mm_prefetch(bptr + i + 16,_MM_HINT_NTA);
			
#else
			//_mm_prefetch(bptr + i + 8,_MM_HINT_NTA);
#endif
#endif

			sieve[bptr[i  ] & 0x0000ffff] -= logp;
			sieve[bptr[i+1] & 0x0000ffff] -= logp;
			sieve[bptr[i+2] & 0x0000ffff] -= logp;
			sieve[bptr[i+3] & 0x0000ffff] -= logp;
			sieve[bptr[i+4] & 0x0000ffff] -= logp;
			sieve[bptr[i+5] & 0x0000ffff] -= logp;
			sieve[bptr[i+6] & 0x0000ffff] -= logp;
			sieve[bptr[i+7] & 0x0000ffff] -= logp;
		}

		for (;i<lpnum;i++)
			sieve[bptr[i] & 0x0000ffff] -= logp;
		
		//point to the next slice of primes
		bptr += (numblocks << (BUCKET_BITS + 1));
		basebucket += (numblocks << 1);
	}

#endif

#ifdef QS_TIMING
		gettimeofday (&qs_timing_stop, NULL);
		qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

		SIEVE_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
		free(qs_timing_diff);
#endif

}

void test_block_siqs(uint8 *sieve, sieve_fb *fb, uint32 start_prime)
{
	//sieve the small primes over a block, in order to estimate their average
	//contribution to a sieve location.  there are probably analytical methods
	//to do this, but this is fast, easy, and gives decent results.
	uint32 prime, root1, root2;
	uint32 i,j;
	uint8 *inner_sieve;
	uint8 logp;
	sieve_fb *fbptr;
	
	//initialize block
	memset(sieve,0,BLOCKSIZE);

	inner_sieve = sieve;

	//break block into sections to keep in L1 cache
	//blocksize should be a multiple of inner_blocksize
	for (i=0;i<BLOCKSIZE;i+=INNER_BLOCKSIZE)
	{
		//for some small number of primes (those that fit in L1 cache/2
		for (j=2; j<start_prime; j++)
		{
			fbptr = fb + j;
			prime = fbptr->prime;
			root1 = fbptr->root1;
			root2 = fbptr->root2;
			logp = fbptr->logprime;

			while (root2 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				inner_sieve[root2] += logp;
				root1 += prime;
				root2 += prime;
			}

			//don't forget the last proot1[i], and compute the roots for the next block
			if (root1 < INNER_BLOCKSIZE)
			{
				inner_sieve[root1] += logp;
				root1 += prime;
				//don't update the root pointers... this is just a test.
			}
		}
		//move inner_sieve to the right to paint a new stripe of the small fb primes in sieve
		inner_sieve += INNER_BLOCKSIZE;
	}
	return;
}
