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

#include "common.h"

#ifdef USE_AVX512F
#include <immintrin.h>
#endif

#ifdef _MSC_VER
#define USE_AVX2
#define GCC_ASM64X
#endif

// protect avx2 code under MSVC builds.  USE_AVX2 should be manually
// enabled at the top of qs.h for MSVC builds on supported hardware
#if defined( USE_AVX2 ) && defined (GCC_ASM64X)


#include "yafu.h"
#include "qs.h"
#include "sieve_macros_32k.h"
#include "sieve_macros_32k_avx2.h"

// vpext uses p0 and p5
// subb uses 2p0156 2p237 p4
// vpext avoids read port conflicts; subb is free to use p237
// subb can be using p1 or p6 while vpext uses p0 and p5
// thus hopefully once things get rolling we have an 
// effective reciprocal throughput of 1.
// short of gather/scatter, not sure if I can do any better.
#define _8P_STEP_SIEVE_AVX2 \
	"vpextrw	$0, %%xmm1, %%eax \n\t"			/* extract root1 */ \
	"vpextrw	$0, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
    "subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$1, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
    "subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$1, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$2, %%xmm1, %%eax \n\t"			/* extract root1 */ \
    "subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$2, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
    "subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$3, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
    "subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$3, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$4, %%xmm1, %%eax \n\t"			/* extract root1 */ \
    "subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$4, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
    "subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$5, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
    "subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$5, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$6, %%xmm1, %%eax \n\t"			/* extract root1 */ \
    "subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$6, %%xmm2, %%ebx \n\t"			/* extract root2 */ \
    "subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$7, %%xmm1, %%ecx \n\t"			/* extract root1 */ \
    "subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpextrw	$7, %%xmm2, %%edi \n\t"			/* extract root2 */ \
	"subb	%%sil, (%%rdx, %%rcx, 1) \n\t"	/* read/modify/write sieve */ \
	"vpaddw	%%xmm0, %%xmm1, %%xmm1 \n\t"			/* increment root1's by primes */ \
	"vpaddw	%%xmm0, %%xmm2, %%xmm2 \n\t"			/* increment root2's by primes */ \
    "subb	%%sil, (%%rdx, %%rdi, 1) \n\t"	/* read/modify/write sieve */ 


#define _8P_FINAL_STEP_SIEVE_AVX2 \
	/* workaround for lack of unsigned greater than... */ \
	"vpsubusw	%%xmm3, %%xmm1, %%xmm4 \n\t"		/* xmm4 := root1 - 32767 */ \
    "vpsubusw	%%xmm3, %%xmm2, %%xmm5 \n\t"		/* xmm5 := root2 - 32767 */ \
	"vpcmpeqw	%%xmm6, %%xmm4, %%xmm4 \n\t"		/* xmm4 := root1 <= 32767 ? 1 : 0 */ \
	"vpcmpeqw	%%xmm6, %%xmm5, %%xmm5 \n\t"		/* xmm5 := root2 <= 32767 ? 1 : 0 */ \
	"vpmovmskb	%%xmm4, %%ecx \n\t" \
	"vpmovmskb	%%xmm5, %%edi \n\t" \
	"testl	$0x2, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* if bits are equal to zero then root1 ! <= blocksize-1; jump */ \
	"vpextrw	$0, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$0, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$1, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$1, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpand	%%xmm0, %%xmm4, %%xmm4 \n\t"			/* clear primes whose root1's are >= blocksize */ \
	"testl	$0x20, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$2, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$2, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x20, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpand	%%xmm0, %%xmm5, %%xmm5 \n\t"			/* clear primes whose root2's are >= blocksize */ \
	"testl	$0x80, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$3, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$3, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x80, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpsubw	%%xmm7, %%xmm4, %%xmm4 \n\t"			/* advance root1's still below blocksize */ \
    "testl	$0x200, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$4, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$4, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x200, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpsubw	%%xmm7, %%xmm5, %%xmm5 \n\t"			/* advance root1's still below blocksize */ \
	"testl	$0x800, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$5, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$5, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x800, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x2000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$6, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$6, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$7, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$7, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */



#ifdef NOTDEF
#define _16P_FINAL_STEP_SIEVE_AVX2 \
	/* workaround for lack of unsigned greater than... */ \
	"vpsubusw	%%ymm3, %%ymm1, %%ymm4 \n\t"		/* xmm4 := root1 - 32767 */ \
    "vpsubusw	%%ymm3, %%ymm2, %%ymm5 \n\t"		/* xmm5 := root2 - 32767 */ \
	"vpcmpeqw	%%ymm6, %%ymm4, %%ymm4 \n\t"		/* xmm4 := root1 <= 32767 ? 1 : 0 */ \
	"vpcmpeqw	%%ymm6, %%ymm5, %%ymm5 \n\t"		/* xmm5 := root2 <= 32767 ? 1 : 0 */ \
	"vpmovmskb	%%ymm4, %%ecx \n\t" \
	"vpmovmskb	%%ymm5, %%edi \n\t" \
	"testl	$0x2, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* if bits are equal to zero then root1 ! <= blocksize-1; jump */ \
	"vpextrw	$0, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$0, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$1, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$1, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpand	%%ymm0, %%ymm4, %%ymm4 \n\t"			/* clear primes whose root1's are >= blocksize */ \
	"testl	$0x20, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$2, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$2, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x20, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpand	%%ymm0, %%ymm5, %%ymm5 \n\t"			/* clear primes whose root2's are >= blocksize */ \
	"testl	$0x80, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$3, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$3, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x80, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpsubw	%%ymm7, %%ymm4, %%ymm4 \n\t"			/* advance root1's still below blocksize */ \
    "testl	$0x200, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$4, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$4, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x200, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "vpsubw	%%ymm7, %%ymm5, %%ymm5 \n\t"			/* advance root1's still below blocksize */ \
	"testl	$0x800, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$5, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$5, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x800, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x2000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$6, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$6, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$7, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$7, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
    "vextracti128 $1,%%ymm1, %%xmm1 \n\t" /* now the high half */ \
    "vextracti128 $1,%%ymm2, %%xmm2 \n\t" /* now the high half */ \
    "sarl   $16, %%ecx \n\t" \
    "sarl   $16, %%edi \n\t" \
    "testl	$0x2, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* if bits are equal to zero then root1 ! <= blocksize-1; jump */ \
	"vpextrw	$0, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$0, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$1, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$1, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x20, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$2, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$2, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x20, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x80, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$3, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$3, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x80, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
    "testl	$0x200, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$4, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$4, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x200, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x800, %%ecx \n\t"				/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$5, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$5, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x800, %%edi \n\t"				/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x2000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$6, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$6, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x2000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */ \
	 \
	"testl	$0x8000, %%ecx \n\t"			/* test if root1 >= blocksize (thus so is root2) */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"vpextrw	$7, %%xmm1, %%eax \n\t"			/* else extract root */ \
    "vpextrw	$7, %%xmm2, %%ebx \n\t"			/* else extract root */ \
	"subb	%%sil, (%%rdx, %%rax, 1) \n\t"	/* read/modify/write sieve */ \
	"testl	$0x8000, %%edi \n\t"			/* test if root2 >= blocksize */ \
	"je 1f \n\t"							/* jump past R/M/W if so */ \
	"subb	%%sil, (%%rdx, %%rbx, 1) \n\t"	/* read/modify/write sieve */ \
	"1: \n\t"								/* both root1 and root2 were >= blocksize */

#endif

#define _FINALIZE_SORT_UPDATE_AVX2 \
	"vpaddw	    %%xmm4, %%xmm1, %%xmm1 \n\t"			/* r1 = r1 + (p - b) */ \
	"vpaddw	    %%xmm5, %%xmm2, %%xmm2 \n\t"			/* r2 = r2 + (p - b) */ \
	"vpmaxsw	%%xmm1, %%xmm2, %%xmm5 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	"vpminsw	%%xmm2, %%xmm1, %%xmm4 \n\t"			/* replace xmm1 with min of root1 and root2 */ \
	"vmovdqa    %%xmm4, (%2) \n\t"				/* write new root1's */ \
	"vmovdqa    %%xmm5, (%3) \n\t"				/* write new root2's */

#define _INIT_AVX2_SMALL_PRIME_SIEVE \
	asm (	\
		"movl	$0x7fff7fff, %%ecx \n\t"		/* load 2 copies of blocksize-1 in a 32bit reg */ \
		"movl	$0x80008000, %%edi \n\t"		/* load 2 copies of blocksize in a 32bit reg */ \
		"vmovd	%%ecx, %%xmm3 \n\t"				/* we need 32k-1 because we do signed comparisons */ \
		"vmovd	%%edi, %%xmm7 \n\t"				\
		"vpxor	    %%xmm6, %%xmm6, %%xmm6 \n\t"			/* xmm6 := 0 */ \
		"vpshufd	$0, %%xmm3, %%xmm3 \n\t"		/* broadcast blocksize-1 to all words of xmm3 */ \
		"vpshufd	$0, %%xmm7, %%xmm7 \n\t"		/* broadcast blocksize to all words of xmm7 */ \
		: \
		:  \
		: "ecx", "edi", "xmm3", "xmm6", "xmm7");

#define _INIT_AVX2_SMALL_PRIME_SIEVE_x16 \
	asm (	\
		"movl	$0x7fff7fff, %%ecx \n\t"		/* load 2 copies of blocksize-1 in a 32bit reg */ \
		"movl	$0x80008000, %%edi \n\t"		/* load 2 copies of blocksize in a 32bit reg */ \
		"vmovd	%%ecx, %%xmm3 \n\t"				/* we need 32k-1 because we do signed comparisons */ \
		"vmovd	%%edi, %%xmm7 \n\t"				\
		"vpxor	    %%ymm6, %%ymm6, %%ymm6 \n\t"			/* xmm6 := 0 */ \
		"vpshufd	$0, %%xmm3, %%xmm3 \n\t"		/* broadcast blocksize-1 to all words of xmm3 */ \
		"vpshufd	$0, %%xmm7, %%xmm7 \n\t"		/* broadcast blocksize to all words of xmm7 */ \
        "vinserti128	$1, %%xmm3, %%ymm3, %%ymm3 \n\t"		/* broadcast blocksize-1 to all words of xmm3 */ \
		"vinserti128	$1, %%xmm7, %%ymm7, %%ymm7 \n\t"		/* broadcast blocksize to all words of xmm7 */ \
		: \
		:  \
		: "ecx", "edi", "xmm3", "xmm6", "xmm7");


#define _FINALIZE_SORT_UPDATE_AVX2_MM \
	"vpaddw	    %%xmm4, %%xmm1, %%xmm1 \n\t"			/* r1 = r1 + (p - b) */ \
	"vpaddw	    %%xmm5, %%xmm2, %%xmm2 \n\t"			/* r2 = r2 + (p - b) */ \
	"vpmaxuw	%%xmm1, %%xmm2, %%xmm5 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	"vpminuw	%%xmm2, %%xmm1, %%xmm4 \n\t"			/* replace xmm1 with min of root1 and root2 */ \
	"vmovdqa    %%xmm4, (%2) \n\t"				/* write new root1's */ \
	"vmovdqa    %%xmm5, (%3) \n\t"				/* write new root2's */

#define _FINALIZE_SORT_UPDATE_AVX2_x16 \
	"vpaddw	    %%ymm4, %%ymm1, %%ymm1 \n\t"			/* r1 = r1 + (p - b) */ \
	"vpaddw	    %%ymm5, %%ymm2, %%ymm2 \n\t"			/* r2 = r2 + (p - b) */ \
	"vpmaxuw	%%ymm1, %%ymm2, %%ymm5 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	"vpminuw	%%ymm2, %%ymm1, %%ymm4 \n\t"			/* replace xmm1 with min of root1 and root2 */ \
	"vmovdqa    %%ymm4, (%2) \n\t"				/* write new root1's */ \
	"vmovdqa    %%ymm5, (%3) \n\t"				/* write new root2's */

#define _AVX2_SMALL_PRIME_SIEVE \
	for (; i < med_B; i += 8) \
		asm (			\
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$15, %%rsi \n\t"				/* logp's range from 15 to 16... call 'em = 15 */ \
			"vmovdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
			\
			_8P_FINAL_STEP_SIEVE_AVX2					\
			_FINALIZE_SORT_UPDATE_AVX2 \
			\
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx", "rdx", \
				"rsi", "rdi", "cc", "memory");

#define _AVX2_SMALL_PRIME_SIEVE_x16 \
	for (; i < med_B; i += 16) \
		asm (			\
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$15, %%rsi \n\t"				/* logp's range from 15 to 16... call 'em = 15 */ \
			"vmovdqa	(%1), %%xmm0 \n\t"				/* bring in 16 primes */ \
			"vmovdqa	(%2), %%xmm1 \n\t"				/* bring in 16 root1's */ \
			"vmovdqa	(%3), %%xmm2 \n\t"				/* bring in 16 root2's */ \
			\
			_8P_FINAL_STEP_SIEVE_AVX2					\
			_FINALIZE_SORT_UPDATE_AVX2 \
			"vmovdqa	16(%1), %%xmm0 \n\t"				/* bring in 16 primes */ \
			"vmovdqa	16(%2), %%xmm1 \n\t"				/* bring in 16 root1's */ \
			"vmovdqa	16(%3), %%xmm2 \n\t"				/* bring in 16 root2's */ \
			\
			_8P_FINAL_STEP_SIEVE_AVX2					\
			"vpaddw	    %%xmm4, %%xmm1, %%xmm1 \n\t"			/* r1 = r1 + (p - b) */ \
	        "vpaddw	    %%xmm5, %%xmm2, %%xmm2 \n\t"			/* r2 = r2 + (p - b) */ \
	        "vpmaxuw	%%xmm1, %%xmm2, %%xmm5 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	        "vpminuw	%%xmm2, %%xmm1, %%xmm4 \n\t"			/* replace xmm1 with min of root1 and root2 */ \
	        "vmovdqa    %%xmm4, 16(%2) \n\t"				/* write new root1's */ \
	        "vmovdqa    %%xmm5, 16(%3) \n\t"				/* write new root2's */ \
			\
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx", "rdx", \
				"rsi", "rdi", "cc", "memory");

#define _AVX2_SMALL_PRIME_SIEVE_32k_DIV3 \
	for (; i < full_fb->fb_32k_div3-8; i += 8) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$13, %%rsi \n\t"				/* logps = 13 in this range */ \
			"vmovdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
			 \
			_8P_STEP_SIEVE_AVX2 \
			_8P_STEP_SIEVE_AVX2 \
			_8P_STEP_SIEVE_AVX2 \
			_8P_FINAL_STEP_SIEVE_AVX2	 \
			_FINALIZE_SORT_UPDATE_AVX2 \
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", \
				"rcx", "rdx", "rsi", "rdi", "cc", "memory");

#define _AVX2_SMALL_PRIME_SIEVE_14b \
	for (; i < full_fb->fb_14bit_B-8; i += 8) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$14, %%rsi \n\t"				/* logp's range from 13 to 14... call 'em = 14 */ \
			"vmovdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
			 \
			_8P_STEP_SIEVE_AVX2 \
			_8P_STEP_SIEVE_AVX2 \
			_8P_FINAL_STEP_SIEVE_AVX2	\
			_FINALIZE_SORT_UPDATE_AVX2 \
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", \
				"rcx", "rdx", "rsi", "rdi", "cc", "memory");

#define _AVX2_SMALL_PRIME_SIEVE_15b \
	for (; i < full_fb->fb_15bit_B-8; i += 8) \
		asm ( \
			"movq	%0,	%%rdx \n\t"					/* sieve array address */ \
			"movq	$15, %%rsi \n\t"				/* logp's range from 14 to 15... call 'em = 15 */ \
			"vmovdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */ \
			"vmovdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */ \
			"vmovdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */ \
				\
			_8P_STEP_SIEVE_AVX2 \
			_8P_FINAL_STEP_SIEVE_AVX2			\
			_FINALIZE_SORT_UPDATE_AVX2 \
			\
			: \
			: "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i) \
			: "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx", "rdx", "rsi", \
				"rdi", "cc", "memory");

#define _AVX512_SMALL_PRIME_SIEVE_15b \
    for (k = i; k < full_fb->fb_15bit_B - 32; k += 32) \
    asm(                                                                                         \
        "movq	%0,	%%rdx \n\t"					/* sieve array address */                             \
        "movq	$15, %%rsi \n\t"				/* logp's range from 14 to 15... call 'em = 15 */     \
        "vmovdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */                           \
        "vmovdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                          \
        "vmovdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                          \
        _8P_STEP_SIEVE_AVX2                                                                           \
        "vmovdqa	16(%1), %%xmm0 \n\t"				/* bring in 8 primes */                       \
        "vmovdqa	16(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                      \
        "vmovdqa	16(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                      \
        _8P_STEP_SIEVE_AVX2                                                                           \
        "vmovdqa	32(%1), %%xmm0 \n\t"				/* bring in 8 primes */                       \
        "vmovdqa	32(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                      \
        "vmovdqa	32(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                      \
        _8P_STEP_SIEVE_AVX2                                                                           \
        "vmovdqa	48(%1), %%xmm0 \n\t"				/* bring in 8 primes */                       \
        "vmovdqa	48(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                      \
        "vmovdqa	48(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                      \
        _8P_STEP_SIEVE_AVX2                                                                           \
        : \
        : "r"(sieve), "r"(fb->prime + k), "r"(fb->root1 + k), "r"(fb->root2 + k)                  \
        : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx",    \
            "rdx", "rsi", "rdi", "cc", "memory"); \

#define _AVX512_15b \
    asm(                                                                                         \
        "movq	%0,	%%rdx \n\t"					/* sieve array address */                             \
        "movq	$15, %%rsi \n\t"				/* logp's range from 14 to 15... call 'em = 15 */     \
        "vmovdqa	(%1), %%xmm0 \n\t"				/* bring in 8 primes */                           \
        "vmovdqa	(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                          \
        "vmovdqa	(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                          \
        _8P_STEP_SIEVE_AVX2                                                                           \
        "vmovdqa	16(%1), %%xmm0 \n\t"				/* bring in 8 primes */                       \
        "vmovdqa	16(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                      \
        "vmovdqa	16(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                      \
        _8P_STEP_SIEVE_AVX2                                                                           \
        "vmovdqa	32(%1), %%xmm0 \n\t"				/* bring in 8 primes */                       \
        "vmovdqa	32(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                      \
        "vmovdqa	32(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                      \
        _8P_STEP_SIEVE_AVX2                                                                           \
        "vmovdqa	48(%1), %%xmm0 \n\t"				/* bring in 8 primes */                       \
        "vmovdqa	48(%2), %%xmm1 \n\t"				/* bring in 8 root1's */                      \
        "vmovdqa	48(%3), %%xmm2 \n\t"				/* bring in 8 root2's */                      \
        _8P_STEP_SIEVE_AVX2                                                                           \
        : \
        : "r"(sieve), "r"(fb->prime + i), "r"(fb->root1 + i), "r"(fb->root2 + i)                  \
        : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "rax", "rbx", "rcx",    \
            "rdx", "rsi", "rdi", "cc", "memory");

// registers clobbered in this loop:
// r9,r10,r14,rax,rcx,rdx,rsi,rdi
// protected registers:
// rbx,r8,r12,r15,r11,r13
// free registers:
// none


#define SIEVE_2X_BLOCK_ASM(id) \
    "vpextrw $" id ", %%xmm2, %%r10d \n\t"	    /* bring in prime */	 \
    "vpextrw $" id ", %%xmm5, %%r14d \n\t"		/* bring in stop value */ \
    "leaq   (%%r10,%%rbx,1),%%rcx	 \n\t"	/* sieve2 = sieve + prime */ \
    "vpextrw $" id ", %%xmm0, %%r9d \n\t"		/* bring in root1 */ \
  	"vpextrw $" id ", %%xmm1, %%edi \n\t"		/* bring in root2 */	 \
  	"vpextrw $" id ", %%xmm3, %%esi \n\t"		/* bring in logp */ \
    "leal   (%%r10,%%r10,1),%%r10d \n\t"	/* 2x prime in r11; root2 prime overwritten */ \
  	"cmpl   %%r14d,%%edi \n\t"				/* root2 >= blocksize-prime? */ \
  	"jae    1f \n\t"						/* jump past loop if so */ \
	"0:  \n\t"								/* sieve to "stop"(r14d) */ \
    "leaq   (%%rbx,%%r9,1),%%rdx \n\t" \
    "leaq   (%%rbx,%%rdi,1),%%rax \n\t"			 \
    "leaq   (%%rcx,%%r9,1),%%r11 \n\t" \
    "leaq   (%%rcx,%%rdi,1),%%r13 \n\t"			 \
    "subb   %%sil,(%%rdx)	 \n\t" \
    "subb   %%sil,(%%rax) \n\t" \
    "prefetcht0 (%%r11) \n\t" \
    "prefetcht0 (%%r13) \n\t" \
    "addl   %%r10d,%%edi \n\t" \
    "addl   %%r10d,%%r9d \n\t" \
    "leaq   (%%rbx,%%r9,1),%%rdx \n\t" \
    "leaq   (%%rbx,%%rdi,1),%%rax \n\t"			 \
    "shrl   $1, %%r10d \n\t"                 /* done with 2x prime */ \
    "subb   %%sil,(%%r11)	 \n\t"	 \
    "subb   %%sil,(%%r13) \n\t" \
    "addl   %%r10d,%%edi \n\t" \
    "addl   %%r10d,%%r9d \n\t" \
    "subb   %%sil,(%%rdx)	 \n\t" \
    "subb   %%sil,(%%rax) \n\t" \
    "shll   $1, %%r10d \n\t"                 /* 2x prime */ \
    "cmpl   %%r14d,%%edi \n\t" \
    "jb     0b \n\t"						/* repeat */ \
    "1:  \n\t" \
    "shrl   $1, %%r10d \n\t"                 /* done with 2x prime */ \
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
    "vpinsrw $" id ", %%edi, %%xmm1, %%xmm1 \n\t"		/* put back new root2 */ \
    "cmpl   $32767,%%r9d \n\t"				/* root1 >= blocksize? */ \
    "ja     1f \n\t"						/* jump past extra root1 block if so */ \
    "subb   %%sil,(%%rbx,%%r9,1) \n\t"      /* need one last write in this block from root 1 */ \
    "addl   %%r10d,%%r9d \n\t"              /* and one last increment by prime */ \
    "1:  \n\t" \
    "vpinsrw $" id ", %%r9d, %%xmm0, %%xmm0 \n\t"		/* put back new root1 */


#define SIEVE_13b_ASM_AVX2 \
	__asm__ ( \
		"movq	%0,%%r12 \n\t"					/* move helperstruct into r12 */ \
		"movl   40(%%r12,1),%%r8d \n\t"			/* move startprime (i) into r8d */ \
		"movq	(%%r12,1),%%rbx \n\t"			/* move sieve into rbx */ \
		"cmpl   44(%%r12,1),%%r8d \n\t"				/* i >= bound? */ \
		"jae    8f \n\t"						/* jump to exit if so */ \
        "movq	16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
        "movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
        "movq   8(%%r12,1),%%rdx \n\t"			/* rdx holds prime pointer */ \
        "movq   32(%%r12,1),%%rax \n\t"			/* rax holds logp pointer */ \
        "movl   $0x80008000,%%ecx \n\t"		    /* blocksize temporarily in rcx */ \
        "vmovd      %%ecx, %%xmm4 \n\t"         /* broadcast blocksize to xmm4 */ \
        "vpshufd	$0, %%xmm4, %%xmm4 \n\t"    /* broadcast blocksize to xmm4 */ \
            /* ==================================================================== */ \
			/* = 2x sieving	loop									              = */ \
			/* ==================================================================== */ \
		"7: \n\t"								/* start of 2x sieving loop */  \
        "vmovdqa (%%r13, %%r8, 2), %%xmm0 \n\t"          /* load 8 root1's */ \
        "vmovdqa (%%r11, %%r8, 2), %%xmm1 \n\t"          /* load 8 root2's */ \
        "vmovdqa (%%rdx, %%r8, 2), %%xmm2 \n\t"          /* load 8 primes's */ \
        "vmovdqa (%%rax, %%r8, 2), %%xmm3 \n\t"          /* load 8 logp's */ \
        "vpsubw     %%xmm2, %%xmm4, %%xmm5 \n\t"    /* xmm5 = blocksize - prime */ \
        "vpsubw     %%xmm2, %%xmm5, %%xmm5 \n\t"    /* xmm5 = blocksize - 2*prime */ \
        /* registers clobbered in this loop: */ \
        /* r9,r10,r14,rax,rcx,rdx,rsi,rdi */ \
  		SIEVE_2X_BLOCK_ASM("0") \
        SIEVE_2X_BLOCK_ASM("1") \
        SIEVE_2X_BLOCK_ASM("2") \
        SIEVE_2X_BLOCK_ASM("3") \
        SIEVE_2X_BLOCK_ASM("4") \
        SIEVE_2X_BLOCK_ASM("5") \
        SIEVE_2X_BLOCK_ASM("6") \
        SIEVE_2X_BLOCK_ASM("7") \
            /* ==================================================================== */ \
            /* = cleanup        									              = */ \
            /* ==================================================================== */ \
        "vpsubw     %%xmm4, %%xmm0, %%xmm0 \n\t"    /* root1 -= blocksize */ \
        "vpsubw     %%xmm4, %%xmm1, %%xmm1 \n\t"    /* root2 -= blocksize */ \
        "movq	16(%%r12,1),%%r13 \n\t"			/* r13 holds root1 pointer */ \
        "movq   24(%%r12,1),%%r11 \n\t"			/* r11 holds root2 pointer */ \
        "vpmaxuw	%%xmm0, %%xmm1, %%xmm5 \n\t"	/* replace xmm2 with max of root1 and root2 */ \
        "vpminuw	%%xmm0, %%xmm1, %%xmm6 \n\t"	/* replace xmm1 with min of root1 and root2 */ \
        "movq   %%r8, %%r14 \n\t" \
        "addq   $8, %%r8 \n\t"                      /* get ready for next iteration */ \    
        "vmovdqa %%xmm6, (%%r13, %%r14, 2) \n\t"    /* store 8 new root1's */ \
        "movq   8(%%r12,1),%%rdx \n\t"			    /* rdx holds prime pointer */ \
        "vmovdqa %%xmm5, (%%r11, %%r14, 2) \n\t"    /* store 8 new root2's */ \
        "movq   32(%%r12,1),%%rax \n\t"			    /* rax holds logp pointer */ \
        "cmpl   44(%%r12,1),%%r8d \n\t" \
        "jb     7b \n\t"						/* repeat 2x sieving loop */ \
        "8: \n\t"													\
        "movl	%%r8d, 40(%%r12,1) \n\t"		/* copy out final value of i */ \
        :																\
        : "g"(&asm_input)												\
        : "rax", "rbx", "rcx", "rdx", "rdi", "rsi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "memory", "cc");


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

void med_sieveblock_32k_avx2(uint8* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32 start_prime, uint8 s_init)
{
    uint32 i;
    uint32 med_B;

    uint32 prime, root1, root2, tmp, stop;
    uint8 logp;

#ifdef USE_AVX512BW
    __m512i vblock = _mm512_set1_epi16(BLOCKSIZE);
#endif

#if defined( TARGET_KNC ) || defined(noUSE_AVX512BW)
    __m512i vlomask = _mm512_set1_epi32(0x000000ff);
    __m512i vhimask = _mm512_set1_epi32(0xffffff00);
    __m512i vlosieve1, vhisieve1, vlosieve2, vhisieve2;
    __m512i vpmul = _mm512_setr_epi32(
        0, 1, 2, 3,
        4, 5, 6, 7,
        8, 9, 10, 11,
        12, 13, 14, 15);

    __m512i vblock = _mm512_set1_epi32(32768);
    __m512i vzero = _mm512_setzero_epi32();
#endif


    helperstruct_t asm_input;

    med_B = full_fb->med_B;

    CLEAN_AVX2;

#ifdef QS_TIMING
    gettimeofday(&qs_timing_start, NULL);
#endif

    //initialize the block
    BLOCK_INIT;

    CLEAN_AVX2;

#if defined(noUSE_AVX512BW)
    // beyond 11 bit we don't need to loop... 16 steps takes
    // us past blocksize all at once when prime > 2048
    for (i = start_prime; i < full_fb->fb_11bit_B; i++)
    {
        __m512i vprime, vroot1, vroot2, vlogp, vidx1, vidx2, v16p;
        __m512i vnextprime, vnextroot1, vnextroot2, vnextidx1, vnextidx2;
        __mmask16 mask1, mask2;
        int steps1 = 0, steps2 = 0;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        vprime = _mm512_set1_epi32(prime);
        vroot1 = _mm512_set1_epi32(root1);
        vroot2 = _mm512_set1_epi32(root2);
        vlogp = _mm512_set1_epi32(logp);
        v16p = _mm512_slli_epi32(vprime, 4);
        vidx2 = _mm512_mullo_epi32(vprime, vpmul);
        vidx1 = _mm512_add_epi32(vroot1, vidx2);
        vidx2 = _mm512_add_epi32(vroot2, vidx2);

        mask2 = _mm512_cmp_epu32_mask(vidx2, vblock, _MM_CMPINT_LT);

        //printf("prime = %u, root1 = %u, root2 = %u, initial mask2 = %08x\n", prime, root1, root2, mask2);

        // while all 16 steps hit the block, loop using only mask2 (larger offset):
        while (mask2 == 0xffff)
        {
            vhisieve2 = _mm512_i32gather_epi32(vidx2, sieve, _MM_SCALE_1);
            vhisieve2 = _mm512_sub_epi8(vhisieve2, vlogp);
            _mm512_i32scatter_epi32(sieve, vidx2, vhisieve2, _MM_SCALE_1);

            vnextidx1 = _mm512_add_epi32(vidx1, v16p);
            vnextidx2 = _mm512_add_epi32(vidx2, v16p);
            mask2 = _mm512_cmp_epu32_mask(vnextidx2, vblock, _MM_CMPINT_LT);

            vhisieve1 = _mm512_i32gather_epi32(vidx1, sieve, _MM_SCALE_1);
            vhisieve1 = _mm512_sub_epi8(vhisieve1, vlogp);
            _mm512_i32scatter_epi32(sieve, vidx1, vhisieve1, _MM_SCALE_1);

            vidx1 = vnextidx1;
            vidx2 = vnextidx2;
            steps1 += 16;
            steps2 += 16;

            //printf("steps = %d, mask2 = %08x\n", steps1, mask2);
        }

        // last iteration using separate mask1 and mask2
        mask1 = _mm512_cmp_epu32_mask(vidx1, vblock, _MM_CMPINT_LT);
        //printf("steps1 = %d, mask1 = %08x, mask2 = %08x\n", steps1, mask1, mask2);
        vhisieve2 = _mm512_mask_i32gather_epi32(vzero, mask2, vidx2, sieve, _MM_SCALE_1);
        vhisieve2 = _mm512_sub_epi8(vhisieve2, vlogp);
        _mm512_mask_i32scatter_epi32(sieve, mask2, vidx2, vhisieve2, _MM_SCALE_1);

        vhisieve1 = _mm512_mask_i32gather_epi32(vzero, mask1, vidx1, sieve, _MM_SCALE_1);
        vhisieve1 = _mm512_sub_epi8(vhisieve1, vlogp);
        _mm512_mask_i32scatter_epi32(sieve, mask1, vidx1, vhisieve1, _MM_SCALE_1);

        // test to see if lower offset (mask1) took more steps
        // and determine the new roots.
        steps1 += _mm_popcnt_u32(mask1);
        steps2 += _mm_popcnt_u32(mask2);

        //printf("steps1 = %d, steps2 = %d, mask1 = %08x, mask2 = %08x\n", steps1, steps2, mask1, mask2);

        if (steps1 > steps2)
        {
            //printf("swapped\n");
            fb->root2[i] = (uint16)(root1 + steps1 * prime - 32768);
            fb->root1[i] = (uint16)(root2 + steps2 * prime - 32768);
        }
        else
        {
            fb->root1[i] = (uint16)(root1 + steps1 * prime - 32768);
            fb->root2[i] = (uint16)(root2 + steps2 * prime - 32768);
        }
    }

    asm_input.logptr = fb->logp;
    asm_input.primeptr = fb->prime;
    asm_input.root1ptr = fb->root1;
    asm_input.root2ptr = fb->root2;
    asm_input.sieve = sieve;
    asm_input.startprime = i;
    asm_input.med_B = full_fb->fb_13bit_B - 8;

    SIEVE_13b_ASM_AVX2;

    i = asm_input.startprime;

    //printf("small prime loop done\n");
    //exit(1);

#elif defined(USE_ASM_SMALL_PRIME_SIEVING)
    // sieve primes less than 2^13 using optimized loops: it becomes
    // inefficient to do fully unrolled sse2 loops as the number of
    // steps through the block increases.  While I didn't specifically
    // test whether 2^13 is the best point to start using loops, it seemed
    // good enough.  any gains if it is not optmial will probably be minimal.

    asm_input.logptr = fb->logp;
    asm_input.primeptr = fb->prime;
    asm_input.root1ptr = fb->root1;
    asm_input.root2ptr = fb->root2;
    asm_input.sieve = sieve;
    asm_input.startprime = start_prime;
    asm_input.med_B = full_fb->fb_10bit_B;

    SIEVE_13b_ASM;

    i = asm_input.startprime;

    asm_input.logptr = fb->logp;
    asm_input.primeptr = fb->prime;
    asm_input.root1ptr = fb->root1;
    asm_input.root2ptr = fb->root2;
    asm_input.sieve = sieve;
    asm_input.startprime = i;
    asm_input.med_B = full_fb->fb_13bit_B - 8;

    SIEVE_13b_ASM_AVX2;

    i = asm_input.startprime;

#else
    for (i = start_prime; i < full_fb->fb_13bit_B - 8; i++)
    {
        uint8* s2;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }
#endif

    // the small prime sieve stops just before prime exceeds 2^13
    // the next sse2 sieve assumes primes exceed 2^13.  since
    // some of the primes in the next set of 8 primes could be less
    // than the cutoff and some are greater than, we have to do this
    // small set of crossover primes manually, one at a time.
    for (; i < full_fb->fb_13bit_B; i++)
    {
        uint8* s2;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    // sieve primes 8 at a time, where 8192 < p < blocksize/3 (10923)
    _INIT_AVX2_SMALL_PRIME_SIEVE;
    _AVX2_SMALL_PRIME_SIEVE_32k_DIV3;

    // the sse2 sieve stops just before prime exceeds blocksize/3
    // the next sse2 sieve assumes primes exceed blocksize/3.  since
    // some of the primes in the next set of 8 primes could be less
    // than the cutoff and some are greater than, we have to do this
    // small set of crossover primes manually, one at a time.
    for (; i < full_fb->fb_32k_div3; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;

        UPDATE_ROOTS;
    }

    // sieve primes 8 at a time, where blocksize/3 < p < 2^14
    _INIT_AVX2_SMALL_PRIME_SIEVE;
    _AVX2_SMALL_PRIME_SIEVE_14b;

    // do this small set of crossover primes manually, one at a time,
    // this time for the 14 bit crossover.
    for (; i < full_fb->fb_14bit_B; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;

        UPDATE_ROOTS;
    }


#ifdef USE_AVX512BW
    
    __m512i vp;
    __m512i vr1;
    __m512i vr2;
    __m512i vi;
    __mmask32 result2;
    __mmask32 result1;
    __mmask32 mfinal;
    uint32 res2;
    uint32 res1;
    uint16* r1;
    uint16* r2;
    int k;

    // sieve primes 8 at a time, 2^14 < p < 2^15
#if 1
    _INIT_AVX2_SMALL_PRIME_SIEVE;
    _AVX2_SMALL_PRIME_SIEVE_15b;
#else
    __m128i vp_128;
    __m128i vr1_128;
    __m128i vr2_128;
    __m128i vi_128;
    __m128i vblock_128 = _mm_set1_epi16(BLOCKSIZE);
    __mmask8 result2_128;
    __mmask8 result1_128;
    __mmask8 mfinal_128;

    logp = 15;
    for (; i < full_fb->fb_15bit_B - 8; i += 8) {
        //printf("loading from index %d\n", i); fflush(stdout);
        r1 = fb->root1 + i;
        r2 = fb->root2 + i;

        // all hit at least once.
        sieve[r1[0]] -= logp;
        sieve[r2[0]] -= logp;
        sieve[r1[1]] -= logp;
        sieve[r2[1]] -= logp;
        sieve[r1[2]] -= logp;
        sieve[r2[2]] -= logp;
        sieve[r1[3]] -= logp;
        sieve[r2[3]] -= logp;
        sieve[r1[4]] -= logp;
        sieve[r2[4]] -= logp;
        sieve[r1[5]] -= logp;
        sieve[r2[5]] -= logp;
        sieve[r1[6]] -= logp;
        sieve[r2[6]] -= logp;
        sieve[r1[7]] -= logp;
        sieve[r2[7]] -= logp;

#ifdef __CYGWIN__
        vp_128 = _mm_load_si128((__m128i *)(fb->prime + i));
        vr1_128 = _mm_load_si128((__m128i *)(fb->root1 + i));
        vr2_128 = _mm_load_si128((__m128i *)(fb->root2 + i));
#else
        vp_128 = _mm_load_epi32((fb->prime + i));
        vr1_128 = _mm_load_epi32((fb->root1 + i));
        vr2_128 = _mm_load_epi32((fb->root2 + i));
#endif

        vr1_128 = _mm_add_epi16(vr1_128, vp_128);
        vr2_128 = _mm_add_epi16(vr2_128, vp_128);

        result2_128 = _mm_cmp_epu16_mask(vr2_128, vblock_128, _MM_CMPINT_LT);
        res2 = result2_128;

        while (res2 > 0) {
            int idx = _trail_zcnt(res2);
            sieve[fb->root2[i + idx]] -= logp;
            sieve[fb->root1[i + idx]] -= logp;
            res2 = _reset_lsb(res2);
        }

        // res1 will have fewer set bits this way, so we
        // have fewer overall loop iterations
        result1_128 = _mm_cmp_epu16_mask(vr1_128, vblock_128, _MM_CMPINT_LT);
        res1 = result1_128 & (~result2_128);

        while (res1 > 0) {
            int idx = _trail_zcnt(res1);
            sieve[fb->root1[i + idx]] -= logp;
            res1 = _reset_lsb(res1);
        }

        vr1_128 = _mm_mask_add_epi16(vr1_128, result1_128, vr1_128, vp_128);
        vr2_128 = _mm_mask_add_epi16(vr2_128, result2_128, vr2_128, vp_128);
        vr1_128 = _mm_sub_epi16(vr1_128, vblock_128);
        vr2_128 = _mm_sub_epi16(vr2_128, vblock_128);

#ifdef __CYGWIN__
        _mm_store_si128((__m128i *)(fb->root1 + i), _mm_min_epu16(vr1_128, vr2_128));
        _mm_store_si128((__m128i *)(fb->root2 + i), _mm_max_epu16(vr1_128, vr2_128));
#else
        _mm_store_epi32(fb->root1 + i, _mm_min_epu16(vr1_128, vr2_128));
        _mm_store_epi32(fb->root2 + i, _mm_max_epu16(vr1_128, vr2_128));
#endif
    }
#endif

    // do this small set of crossover primes manually, one at a time,
    // this time for the 15 bit crossover.
    for (i = full_fb->fb_15bit_B - 8; i < med_B; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        if ((prime > 32768) && ((i & 15) == 0))
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;

        UPDATE_ROOTS;
    }

    // sieve primes 32 at a time, 2^15 < p < med_B
    logp = 15;
__assume_aligned(fb->prime, 64); 
__assume_aligned(fb->root1, 64); 
__assume_aligned(fb->root2, 64); 
    for (; i < med_B - 32; i += 32) {
        //printf("loading from index %d\n", i); fflush(stdout);
        vp = _mm512_load_epi32((fb->prime + i));
        vr1 = _mm512_load_epi32((fb->root1 + i));
        vr2 = _mm512_load_epi32((fb->root2 + i));

        result2 = _mm512_cmp_epu16_mask(vr2, vblock, _MM_CMPINT_LT);
        res2 = result2;

        while (res2 > 0) {
            int idx = _trail_zcnt(res2);
            sieve[fb->root2[i + idx]] -= logp;
            sieve[fb->root1[i + idx]] -= logp;
            res2 = _reset_lsb(res2);
        }

        // res1 will have fewer set bits this way, so we
        // have fewer overall loop iterations
        result1 = _mm512_cmp_epu16_mask(vr1, vblock, _MM_CMPINT_LT);
        res1 = result1 & (~result2);

        while (res1 > 0) {
            int idx = _trail_zcnt(res1);
            sieve[fb->root1[i + idx]] -= logp;
            res1 = _reset_lsb(res1);
        }

        vr1 = _mm512_mask_add_epi16(vr1, result1, vr1, vp);
        vr2 = _mm512_mask_add_epi16(vr2, result2, vr2, vp);
        vr1 = _mm512_sub_epi16(vr1, vblock);
        vr2 = _mm512_sub_epi16(vr2, vblock);
        _mm512_store_epi32(fb->root1 + i, _mm512_min_epu16(vr1, vr2));
        _mm512_store_epi32(fb->root2 + i, _mm512_max_epu16(vr1, vr2));
    }

    // do the final iteration with any excessive indices masked.
    vi = _mm512_set1_epi16(i);
    vi = _mm512_add_epi16(vi, _mm512_set_epi16(0, 1, 2, 3, 4, 5, 6, 7,
        8, 9, 10, 11, 12, 13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23,
        24, 25, 26, 27, 28, 29, 30, 31));
    mfinal = _mm512_cmp_epi16_mask(vi, _mm512_set1_epi16(med_B), _MM_CMPINT_LT);

    vp = _mm512_maskz_load_epi32(mfinal, (fb->prime + i));
    vr1 = _mm512_maskz_load_epi32(mfinal, (fb->root1 + i));
    vr2 = _mm512_maskz_load_epi32(mfinal, (fb->root2 + i));

    result2 = _mm512_mask_cmp_epu16_mask(mfinal, vr2, vblock, _MM_CMPINT_LT);
    res2 = result2;

    while (res2 > 0) {
        int idx = _trail_zcnt(res2);
        sieve[fb->root2[i + idx]] -= logp;
        sieve[fb->root1[i + idx]] -= logp;
        res2 = _reset_lsb(res2);
    }

    // res1 will have fewer set bits this way, so we
    // have fewer overall loop iterations
    result1 = _mm512_mask_cmp_epu16_mask(mfinal, vr1, vblock, _MM_CMPINT_LT);
    res1 = result1 & (~result2);

    while (res1 > 0) {
        int idx = _trail_zcnt(res1);
        sieve[fb->root1[i + idx]] -= logp;
        res1 = _reset_lsb(res1);
    }

    vr1 = _mm512_mask_add_epi16(vr1, result1, vr1, vp);
    vr2 = _mm512_mask_add_epi16(vr2, result2, vr2, vp);
    vr1 = _mm512_sub_epi16(vr1, vblock);
    vr2 = _mm512_sub_epi16(vr2, vblock);
    _mm512_mask_store_epi32(fb->root1 + i, mfinal, _mm512_min_epu16(vr1, vr2));
    _mm512_mask_store_epi32(fb->root2 + i, mfinal, _mm512_max_epu16(vr1, vr2));

#else

    // sieve primes 8 at a time, 2^14 < p < 2^15
    _INIT_AVX2_SMALL_PRIME_SIEVE;
    _AVX2_SMALL_PRIME_SIEVE_15b;


    // do this small set of crossover primes manually, one at a time,
    // this time for the 15 bit crossover.
    for (i = full_fb->fb_15bit_B - 8; i < med_B; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        if ((prime > 32768) && ((i & 15) == 0))
            break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;

        UPDATE_ROOTS;
    }

	// sieve primes 8 at a time, 2^15 < p < med_B
	_INIT_AVX2_SMALL_PRIME_SIEVE;
	_AVX2_SMALL_PRIME_SIEVE;

#endif
    CLEAN_AVX2;


#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
    SIEVE_STG1 += yafu_difftime (&qs_timing_start, &qs_timing_stop);
	gettimeofday(&qs_timing_start, NULL);
#endif

	return;

}

#endif // USE_AVX2

