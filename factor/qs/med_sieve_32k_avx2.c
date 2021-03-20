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

#include "qs_impl.h"
#include "sieve_macros_32k.h"
#include "sieve_macros_32k_avx2.h"

// asm and sieve routines for linux/mingw
#if defined( USE_AVX2 ) && defined (GCC_ASM64X)

#include <immintrin.h>

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
    for (k = i; k < full_fb->fb_15bit_B - 32; k += 32) { \
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
                                                                                      \
        vp = _mm512_load_epi32((fb->prime + k));                                      \
        vr1 = _mm512_load_epi32((fb->root1 + k));                                     \
        vr2 = _mm512_load_epi32((fb->root2 + k));                                     \
                                                                                        \
        vr1 = _mm512_add_epi16(vr1, vp);                                                \
        vr2 = _mm512_add_epi16(vr2, vp);                                                \
                                                                                      \
        result2 = _mm512_cmp_epu16_mask(vr2, vblock, _MM_CMPINT_LT);                  \
        res2 = result2;                                                               \
                                                                                      \
        while (res2 > 0) {                                                            \
            int idx = _trail_zcnt(res2);                                              \
            sieve[fb->root2[k + idx]] -= logp;                                        \
            sieve[fb->root1[k + idx]] -= logp;                                        \
            res2 = _reset_lsb(res2);                                                  \
        }                                                                             \
                                                                                      \
        /* res1 will have fewer set bits this way, so we       */                     \
        /* have fewer overall loop iterations                  */                     \
        result1 = _mm512_cmp_epu16_mask(vr1, vblock, _MM_CMPINT_LT);                  \
        res1 = result1 & (~result2);                                                  \
                                                                                      \
        while (res1 > 0) {                                                            \
            int idx = _trail_zcnt(res1);                                              \
            sieve[fb->root1[k + idx]] -= logp;                                        \
            res1 = _reset_lsb(res1);                                                  \
        }                                                                             \
                                                                                      \
        vr1 = _mm512_mask_add_epi16(vr1, result1, vr1, vp);                           \
        vr2 = _mm512_mask_add_epi16(vr2, result2, vr2, vp);                           \
        vr1 = _mm512_sub_epi16(vr1, vblock);                                          \
        vr2 = _mm512_sub_epi16(vr2, vblock);                                          \
        _mm512_store_epi32(fb->root1 + k, _mm512_min_epu16(vr1, vr2));                \
        _mm512_store_epi32(fb->root2 + k, _mm512_max_epu16(vr1, vr2));                \
    } 
    
//#define _AVX512_SMALL_PRIME_SIEVE_15b_2 
    


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


        typedef struct
        {
            uint8_t* sieve;					//0
            uint16_t* primeptr;				//8
            uint16_t* root1ptr;				//16
            uint16_t* root2ptr;				//24
            uint16_t* logptr;					//32
            uint32_t startprime;				//40
            uint32_t med_B;					//44
        } helperstruct_t;


void med_sieveblock_32k_avx2(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    uint32_t i;
    uint32_t med_B;

    uint32_t prime, root1, root2, tmp, stop;
    uint8_t logp;

#ifdef USE_AVX512BW
    __m512i vblock = _mm512_set1_epi16(BLOCKSIZE);
#endif

#if defined( TARGET_KNC ) || defined(USE_AVX512BW)
    __m512i vzero = _mm512_setzero_epi32();

    ALIGNED_MEM uint16_t r_id1[32];
    ALIGNED_MEM uint16_t r_id2[32];

    uint32_t bound = full_fb->fb_15bit_B - 32;

    if (full_fb->med_B > (full_fb->fb_15bit_B + 32))
        bound += 32;
#else
    uint32_t bound = full_fb->fb_13bit_B - 8;
#endif


    helperstruct_t asm_input;

    med_B = full_fb->med_B;

    //initialize the block
    BLOCK_INIT;

#if defined(USE_AVX512BW)

    for (i = start_prime; i < bound; i += 32)
    {
        __m512i vprime, vroot1, vroot2;
        __mmask32 valid_mask_1, valid_mask_2, initial_mask;
        uint32_t msk_2;
        int pos;

        vprime = _mm512_loadu_si512(fb->prime + i);
        vroot1 = _mm512_loadu_si512(fb->root1 + i);
        vroot2 = _mm512_loadu_si512(fb->root2 + i);
        logp = fb->logp[i]; // approximate the next 32 logp's as equal to this one.

        // we don't sieve primes that are part of the poly
        valid_mask_1 = initial_mask = valid_mask_2 = _mm512_cmpgt_epu16_mask(vprime, vzero);

        // make it so we write to a dummy sieve location for non-sieved primes
        vroot1 = _mm512_mask_add_epi16(vroot1, ~initial_mask, vroot1, vblock);
        vroot2 = _mm512_mask_add_epi16(vroot2, ~initial_mask, vroot2, vblock);

        // until things start to drop off the end of the interval, 
        // simply dump in all logs.
        while (valid_mask_2 == initial_mask)
        {
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            sieve[r_id1[0]] -= logp;
            sieve[r_id1[1]] -= logp;
            sieve[r_id1[2]] -= logp;
            sieve[r_id1[3]] -= logp;
            sieve[r_id1[4]] -= logp;
            sieve[r_id1[5]] -= logp;
            sieve[r_id1[6]] -= logp;
            sieve[r_id1[7]] -= logp;
            sieve[r_id1[8]] -= logp;
            sieve[r_id1[9]] -= logp;
            sieve[r_id1[10]] -= logp;
            sieve[r_id1[11]] -= logp;
            sieve[r_id1[12]] -= logp;
            sieve[r_id1[13]] -= logp;
            sieve[r_id1[14]] -= logp;
            sieve[r_id1[15]] -= logp;
            sieve[r_id1[16]] -= logp;
            sieve[r_id1[17]] -= logp;
            sieve[r_id1[18]] -= logp;
            sieve[r_id1[19]] -= logp;
            sieve[r_id1[20]] -= logp;
            sieve[r_id1[21]] -= logp;
            sieve[r_id1[22]] -= logp;
            sieve[r_id1[23]] -= logp;
            sieve[r_id1[24]] -= logp;
            sieve[r_id1[25]] -= logp;
            sieve[r_id1[26]] -= logp;
            sieve[r_id1[27]] -= logp;
            sieve[r_id1[28]] -= logp;
            sieve[r_id1[29]] -= logp;
            sieve[r_id1[30]] -= logp;
            sieve[r_id1[31]] -= logp;
            sieve[r_id2[0]] -= logp;
            sieve[r_id2[1]] -= logp;
            sieve[r_id2[2]] -= logp;
            sieve[r_id2[3]] -= logp;
            sieve[r_id2[4]] -= logp;
            sieve[r_id2[5]] -= logp;
            sieve[r_id2[6]] -= logp;
            sieve[r_id2[7]] -= logp;
            sieve[r_id2[8]] -= logp;
            sieve[r_id2[9]] -= logp;
            sieve[r_id2[10]] -= logp;
            sieve[r_id2[11]] -= logp;
            sieve[r_id2[12]] -= logp;
            sieve[r_id2[13]] -= logp;
            sieve[r_id2[14]] -= logp;
            sieve[r_id2[15]] -= logp;
            sieve[r_id2[16]] -= logp;
            sieve[r_id2[17]] -= logp;
            sieve[r_id2[18]] -= logp;
            sieve[r_id2[19]] -= logp;
            sieve[r_id2[20]] -= logp;
            sieve[r_id2[21]] -= logp;
            sieve[r_id2[22]] -= logp;
            sieve[r_id2[23]] -= logp;
            sieve[r_id2[24]] -= logp;
            sieve[r_id2[25]] -= logp;
            sieve[r_id2[26]] -= logp;
            sieve[r_id2[27]] -= logp;
            sieve[r_id2[28]] -= logp;
            sieve[r_id2[29]] -= logp;
            sieve[r_id2[30]] -= logp;
            sieve[r_id2[31]] -= logp;

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }
        
        // as roots start to exceed the block size, selectively 
        // dump in logs
        while (valid_mask_2 > 0)
        {
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            msk_2 = valid_mask_2;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos]] -= logp;
                sieve[r_id1[pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // now all larger roots are invalid.  Last iteration for 
        // possibly still valid root1s.  If they are still valid, 
        // record the sieve hit, advance them, and swap with the
        // other root
        _mm512_store_si512(r_id1, vroot1);
        valid_mask_2 = valid_mask_1 & _mm512_cmplt_epu16_mask(vroot1, vblock);
        msk_2 = valid_mask_2;

        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }

        // reduce both roots and store back for the next block
        vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
        vroot1 = _mm512_sub_epi16(vroot1, vblock);
        vroot2 = _mm512_sub_epi16(vroot2, vblock);
        _mm512_storeu_si512(fb->root1 + i, _mm512_min_epu16(vroot1, vroot2));
        _mm512_storeu_si512(fb->root2 + i, _mm512_max_epu16(vroot1, vroot2));
    }


    if ((med_B - i) < 32)
    {
        bound = med_B;
    }

    for (; i < bound; i++)
    {
        uint8_t* s2;

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

#elif defined(USE_AVX2) && defined(USE_BMI2)

//#define TEST_SIEVE

    bound = full_fb->fb_15bit_B - 16;

    //if (full_fb->med_B > (full_fb->fb_15bit_B + 16))
    //bound += 16;

    //printf("bound = %u, startprime = %u\n", bound, start_prime);

    for (i = start_prime; i < bound; i++)
    {
        uint8_t* s2;

        if ((i & 15) == 0)
            break;

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

    __m256i vzero = _mm256_setzero_si256();
    __m256i vblock = _mm256_set1_epi16(BLOCKSIZE);
    __m256i vprime, vroot1, vroot2, vtmp1, vtmp2;
    __m256i valid_mask_1, valid_mask_2, initial_mask;
    ALIGNED_MEM uint16_t r_id1[16];
    ALIGNED_MEM uint16_t r_id2[16];

    for (; i < bound; i += 16)
    {
        uint32_t msk_2;
        int pos;
        int initial_mask_32;

        vprime = _mm256_load_si256((__m256i *)(fb->prime + i));
        vroot1 = _mm256_load_si256((__m256i *)(fb->root1 + i));
        vroot2 = _mm256_load_si256((__m256i *)(fb->root2 + i));
        logp = fb->logp[i]; // approximate the next 32 logp's as equal to this one.

#ifdef TEST_SIEVE
        int j;
        printf("root1s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", fb->root1[i + j]);
        }
        printf("\n");
        printf("root2s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", fb->root2[i + j]);
        }
        printf("\n");
        printf("primes @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", fb->prime[i + j]);
        }
        printf("\n");
        _mm256_store_si256((__m256i*)r_id1, vblock);
        printf("vblock @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id1[j]);
        }
        printf("\n");
#endif

        // we don't sieve primes that are part of the poly
        valid_mask_1 = initial_mask = valid_mask_2 = _mm256_cmpgt_epi16(vprime, vzero);

        // make it so we write to a dummy sieve location for non-sieved primes
        vtmp1 = _mm256_andnot_si256(valid_mask_2, vblock);
        vtmp2 = _mm256_andnot_si256(valid_mask_2, vblock);
        vroot1 = _mm256_add_epi16(vtmp1, vroot1);
        vroot2 = _mm256_add_epi16(vtmp2, vroot2);
        initial_mask_32 = _mm256_movemask_epi8(initial_mask);

#ifdef TEST_SIEVE
        _mm256_store_si256((__m256i*)r_id1, initial_mask);
        printf("initial mask\n");
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id1[j]);
        }
        printf("\n");
#endif

        // until things start to drop off the end of the interval, 
        // simply dump in all logs.
        do
        {
            _mm256_store_si256((__m256i *)r_id1, vroot1);
            _mm256_store_si256((__m256i *)r_id2, vroot2);

            sieve[r_id1[0]] -= logp;
            sieve[r_id1[1]] -= logp;
            sieve[r_id1[2]] -= logp;
            sieve[r_id1[3]] -= logp;
            sieve[r_id1[4]] -= logp;
            sieve[r_id1[5]] -= logp;
            sieve[r_id1[6]] -= logp;
            sieve[r_id1[7]] -= logp;
            sieve[r_id1[8]] -= logp;
            sieve[r_id1[9]] -= logp;
            sieve[r_id1[10]] -= logp;
            sieve[r_id1[11]] -= logp;
            sieve[r_id1[12]] -= logp;
            sieve[r_id1[13]] -= logp;
            sieve[r_id1[14]] -= logp;
            sieve[r_id1[15]] -= logp;
            sieve[r_id2[0]] -= logp;
            sieve[r_id2[1]] -= logp;
            sieve[r_id2[2]] -= logp;
            sieve[r_id2[3]] -= logp;
            sieve[r_id2[4]] -= logp;
            sieve[r_id2[5]] -= logp;
            sieve[r_id2[6]] -= logp;
            sieve[r_id2[7]] -= logp;
            sieve[r_id2[8]] -= logp;
            sieve[r_id2[9]] -= logp;
            sieve[r_id2[10]] -= logp;
            sieve[r_id2[11]] -= logp;
            sieve[r_id2[12]] -= logp;
            sieve[r_id2[13]] -= logp;
            sieve[r_id2[14]] -= logp;
            sieve[r_id2[15]] -= logp;

            vroot1 = _mm256_add_epi16(vroot1, vprime);
            vroot2 = _mm256_add_epi16(vroot2, vprime);

            vtmp2 = _mm256_srli_epi16(vroot2, 15);
            vtmp2 = _mm256_or_si256(_mm256_cmpgt_epi16(vtmp2, vzero),
                _mm256_cmpeq_epi16(vroot2, vblock));

            valid_mask_2 = _mm256_andnot_si256(vtmp2, valid_mask_2);
        } while (_mm256_movemask_epi8(valid_mask_2) == initial_mask_32);

        // zero out the primes where roots have exceeded the block
        vprime = _mm256_and_si256(valid_mask_2, vprime);

#ifdef TEST_SIEVE
        _mm256_store_si256((__m256i*)r_id1, vroot1);
        _mm256_store_si256((__m256i*)r_id2, vroot2);
        printf("after first loop\n");
        printf("root1s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id1[j]);
        }
        printf("\n");
        printf("root2s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id2[j]);
        }
        printf("\n");
        _mm256_store_si256((__m256i*)r_id2, vprime);
        printf("primes @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id2[j]);
        }
        printf("\n");
#endif

        // as roots start to exceed the block size, selectively 
        // dump in logs
        while ((msk_2 = _mm256_movemask_epi8(valid_mask_2)) > 0)
        {
            _mm256_store_si256((__m256i *)r_id1, vroot1);
            _mm256_store_si256((__m256i *)r_id2, vroot2);

            msk_2 &= 0xaaaaaaaa;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos >> 1]] -= logp;
                sieve[r_id1[pos >> 1]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }

            vroot1 = _mm256_add_epi16(vroot1, vprime);
            vroot2 = _mm256_add_epi16(vroot2, vprime);

            vtmp2 = _mm256_srli_epi16(vroot2, 15);
            vtmp2 = _mm256_or_si256(_mm256_cmpgt_epi16(vtmp2, vzero),
                _mm256_cmpeq_epi16(vroot2, vblock));

            valid_mask_2 = _mm256_andnot_si256(vtmp2, valid_mask_2);

            // zero out the primes where roots have exceeded the block
            vprime = _mm256_and_si256(valid_mask_2, vprime);
        }

#ifdef TEST_SIEVE
        _mm256_store_si256((__m256i*)r_id1, vroot1);
        _mm256_store_si256((__m256i*)r_id2, vroot2);
        printf("after second loop\n");
        printf("root1s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id1[j]);
        }
        printf("\n");
        printf("root2s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id2[j]);
        }
        printf("\n");
        _mm256_store_si256((__m256i*)r_id2, vprime);
        printf("primes @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id2[j]);
        }
        printf("\n");
#endif

        // now all larger roots are invalid.  Last iteration for 
        // possibly still valid root1s.  If they are still valid, 
        // record the sieve hit, advance them, and swap with the
        // other root
        _mm256_store_si256((__m256i*)r_id1, vroot1);

        vtmp2 = _mm256_srli_epi16(vroot1, 15);
        vtmp2 = _mm256_or_si256(_mm256_cmpgt_epi16(vtmp2, vzero),
            _mm256_cmpeq_epi16(vroot1, vblock));
        valid_mask_2 = _mm256_andnot_si256(vtmp2, valid_mask_1);
        msk_2 = _mm256_movemask_epi8(valid_mask_2);

        msk_2 &= 0xaaaaaaaa;
        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos >> 1]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }

        // reload the primes, then zero out the ones where root1 has
        // already exceeded the block.
        vprime = _mm256_load_si256((__m256i*)(fb->prime + i));
        vprime = _mm256_and_si256(valid_mask_2, vprime);
        vprime = _mm256_and_si256(initial_mask, vprime);

        // reduce both roots and store back for the next block
        vroot1 = _mm256_add_epi16(vroot1, vprime);

#ifdef TEST_SIEVE
        _mm256_store_si256((__m256i*)r_id1, vroot1);
        _mm256_store_si256((__m256i*)r_id2, vroot2);
        printf("after last loop\n");
        printf("root1s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id1[j]);
        }
        printf("\n");
        printf("root2s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id2[j]);
        }
        printf("\n");
        _mm256_store_si256((__m256i*)r_id2, vprime);
        printf("primes @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id2[j]);
        }
        printf("\n");
#endif


        vroot1 = _mm256_sub_epi16(vroot1, vblock);
        vroot2 = _mm256_sub_epi16(vroot2, vblock);
        _mm256_store_si256((__m256i *)(fb->root1 + i), _mm256_min_epu16(vroot1, vroot2));
        _mm256_store_si256((__m256i *)(fb->root2 + i), _mm256_max_epu16(vroot1, vroot2));

#ifdef TEST_SIEVE
        printf("after storage\n");
        printf("root1s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", fb->root1[i+j]);
        }
        printf("\n");
        printf("root2s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", fb->root1[i+j]);
        }
        printf("\n");
        //exit(1);
#endif
    }

    if ((med_B - i) < 32)
    {
        bound = med_B;
    }

    for (; i < bound; i++)
    {
        uint8_t* s2;

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

    for (; i < full_fb->fb_13bit_B; i++)
    {
        uint8_t* s2;

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

#else
    for (i = start_prime; i < full_fb->fb_13bit_B - 8; i++)
    {
        uint8_t* s2;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    // the small prime sieve stops just before prime exceeds 2^13
    // the next sse2 sieve assumes primes exceed 2^13.  since
    // some of the primes in the next set of 8 primes could be less
    // than the cutoff and some are greater than, we have to do this
    // small set of crossover primes manually, one at a time.
    for (; i < full_fb->fb_13bit_B; i++)
    {
        uint8_t* s2;

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
#endif

    


#ifdef USE_AVX512BW
    
    __m512i vp;
    __m512i vr1;
    __m512i vr2;
    __mmask32 result2;
    __mmask32 result1;
    uint32_t res2;
    uint32_t res1;

    //for (; i < full_fb->fb_15bit_B - 32; i += 32)
    //{
    //    __m512i vprime, vroot1, vroot2;
    //    __mmask32 valid_mask_1, valid_mask_2, initial_mask;
    //    uint32_t msk_2;
    //    int pos;
    //
    //    vprime = _mm512_loadu_si512(fb->prime + i);
    //    vroot1 = _mm512_loadu_si512(fb->root1 + i);
    //    vroot2 = _mm512_loadu_si512(fb->root2 + i);
    //    logp = fb->logp[i]; // approximate the next 32 logp's as equal to this one.
    //
    //    // we don't sieve primes that are part of the poly
    //    valid_mask_1 = initial_mask = valid_mask_2 = _mm512_cmpgt_epu16_mask(vprime, vzero);
    //
    //    // as roots start to exceed the block size, selectively 
    //    // dump in logs
    //    while (valid_mask_2 > 0)
    //    {
    //        _mm512_store_si512(r_id1, vroot1);
    //        _mm512_store_si512(r_id2, vroot2);
    //
    //        msk_2 = valid_mask_2;
    //        while (msk_2 > 0) {
    //            pos = _trail_zcnt(msk_2);
    //            sieve[r_id2[pos]] -= logp;
    //            sieve[r_id1[pos]] -= logp;
    //            msk_2 = _reset_lsb(msk_2);
    //        }
    //
    //        vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
    //        vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);
    //
    //        valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
    //    }
    //
    //    // now all larger roots are invalid.  Last iteration for 
    //    // possibly still valid root1s.  If they are still valid, 
    //    // record the sieve hit, advance them, and swap with the
    //    // other root
    //    _mm512_store_si512(r_id1, vroot1);
    //    valid_mask_2 = valid_mask_1 & _mm512_cmplt_epu16_mask(vroot1, vblock);
    //    msk_2 = valid_mask_2;
    //
    //    while (msk_2 > 0) {
    //        pos = _trail_zcnt(msk_2);
    //        sieve[r_id1[pos]] -= logp;
    //        msk_2 = _reset_lsb(msk_2);
    //    }
    //
    //    // reduce both roots and store back for the next block
    //    vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
    //    vroot1 = _mm512_sub_epi16(vroot1, vblock);
    //    vroot2 = _mm512_sub_epi16(vroot2, vblock);
    //    _mm512_storeu_si512(fb->root1 + i, _mm512_min_epu16(vroot1, vroot2));
    //    _mm512_storeu_si512(fb->root2 + i, _mm512_max_epu16(vroot1, vroot2));
    //}

    // sieve primes 32 at a time, 2^15 < p < med_B
    logp = 15;
    for (; i < med_B - 32; i += 32) {
        //printf("med_B 32x loop @ %d\n\n", i); fflush(stdout);
        vp = _mm512_loadu_si512((fb->prime + i));
        vr1 = _mm512_loadu_si512((fb->root1 + i));
        vr2 = _mm512_loadu_si512((fb->root2 + i));

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
        _mm512_storeu_si512(fb->root1 + i, _mm512_min_epu16(vr1, vr2));
        _mm512_storeu_si512(fb->root2 + i, _mm512_max_epu16(vr1, vr2));
    }


    for (; i < med_B; i++)
    {
        //printf("med_B loop @ %d\n\n", i); fflush(stdout);
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;

        UPDATE_ROOTS;
    }

#elif defined(USE_AVX2) && defined(USE_BMI2)

    // do this small set of crossover primes manually, one at a time,
    // this time for the 15 bit crossover.
    for (; i < med_B; i++)
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

#else
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

    return;

}



#endif // USE_AVX2

// asm and sieve routines for msvc
#if defined( USE_AVX2 ) && defined(_MSC_VER)

	#include <immintrin.h>

void med_sieveblock_32k_avx2(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    uint32_t i;
    uint32_t med_B;

    uint32_t prime, root1, root2, tmp, stop;
    uint8_t logp;

#ifdef USE_AVX512BW
    __m512i vblock = _mm512_set1_epi16(BLOCKSIZE);
#endif

#if defined( TARGET_KNC ) || defined(USE_AVX512BW)
    __m512i vzero = _mm512_setzero_epi32();

    ALIGNED_MEM uint16_t r_id1[32];
    ALIGNED_MEM uint16_t r_id2[32];

    uint32_t bound = full_fb->fb_15bit_B - 32;

    if (full_fb->med_B > (full_fb->fb_15bit_B + 32))
        bound += 32;
#else
    uint32_t bound = full_fb->fb_13bit_B - 8;
#endif

    med_B = full_fb->med_B;

    //initialize the block
    BLOCK_INIT;

#if defined(USE_AVX512BW)

    for (i = start_prime; i < bound; i += 32)
    {
        __m512i vprime, vroot1, vroot2;
        __mmask32 valid_mask_1, valid_mask_2, initial_mask;
        uint32_t msk_2;
        int pos;

        vprime = _mm512_loadu_si512(fb->prime + i);
        vroot1 = _mm512_loadu_si512(fb->root1 + i);
        vroot2 = _mm512_loadu_si512(fb->root2 + i);
        logp = fb->logp[i]; // approximate the next 32 logp's as equal to this one.

        // we don't sieve primes that are part of the poly
        valid_mask_1 = initial_mask = valid_mask_2 = _mm512_cmpgt_epu16_mask(vprime, vzero);

        // make it so we write to a dummy sieve location for non-sieved primes
        vroot1 = _mm512_mask_add_epi16(vroot1, ~initial_mask, vroot1, vblock);
        vroot2 = _mm512_mask_add_epi16(vroot2, ~initial_mask, vroot2, vblock);

        // until things start to drop off the end of the interval, 
        // simply dump in all logs.
        while (valid_mask_2 == initial_mask)
        {
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            sieve[r_id1[0]] -= logp;
            sieve[r_id1[1]] -= logp;
            sieve[r_id1[2]] -= logp;
            sieve[r_id1[3]] -= logp;
            sieve[r_id1[4]] -= logp;
            sieve[r_id1[5]] -= logp;
            sieve[r_id1[6]] -= logp;
            sieve[r_id1[7]] -= logp;
            sieve[r_id1[8]] -= logp;
            sieve[r_id1[9]] -= logp;
            sieve[r_id1[10]] -= logp;
            sieve[r_id1[11]] -= logp;
            sieve[r_id1[12]] -= logp;
            sieve[r_id1[13]] -= logp;
            sieve[r_id1[14]] -= logp;
            sieve[r_id1[15]] -= logp;
            sieve[r_id1[16]] -= logp;
            sieve[r_id1[17]] -= logp;
            sieve[r_id1[18]] -= logp;
            sieve[r_id1[19]] -= logp;
            sieve[r_id1[20]] -= logp;
            sieve[r_id1[21]] -= logp;
            sieve[r_id1[22]] -= logp;
            sieve[r_id1[23]] -= logp;
            sieve[r_id1[24]] -= logp;
            sieve[r_id1[25]] -= logp;
            sieve[r_id1[26]] -= logp;
            sieve[r_id1[27]] -= logp;
            sieve[r_id1[28]] -= logp;
            sieve[r_id1[29]] -= logp;
            sieve[r_id1[30]] -= logp;
            sieve[r_id1[31]] -= logp;
            sieve[r_id2[0]] -= logp;
            sieve[r_id2[1]] -= logp;
            sieve[r_id2[2]] -= logp;
            sieve[r_id2[3]] -= logp;
            sieve[r_id2[4]] -= logp;
            sieve[r_id2[5]] -= logp;
            sieve[r_id2[6]] -= logp;
            sieve[r_id2[7]] -= logp;
            sieve[r_id2[8]] -= logp;
            sieve[r_id2[9]] -= logp;
            sieve[r_id2[10]] -= logp;
            sieve[r_id2[11]] -= logp;
            sieve[r_id2[12]] -= logp;
            sieve[r_id2[13]] -= logp;
            sieve[r_id2[14]] -= logp;
            sieve[r_id2[15]] -= logp;
            sieve[r_id2[16]] -= logp;
            sieve[r_id2[17]] -= logp;
            sieve[r_id2[18]] -= logp;
            sieve[r_id2[19]] -= logp;
            sieve[r_id2[20]] -= logp;
            sieve[r_id2[21]] -= logp;
            sieve[r_id2[22]] -= logp;
            sieve[r_id2[23]] -= logp;
            sieve[r_id2[24]] -= logp;
            sieve[r_id2[25]] -= logp;
            sieve[r_id2[26]] -= logp;
            sieve[r_id2[27]] -= logp;
            sieve[r_id2[28]] -= logp;
            sieve[r_id2[29]] -= logp;
            sieve[r_id2[30]] -= logp;
            sieve[r_id2[31]] -= logp;

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // as roots start to exceed the block size, selectively 
        // dump in logs
        while (valid_mask_2 > 0)
        {
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            msk_2 = valid_mask_2;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos]] -= logp;
                sieve[r_id1[pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // now all larger roots are invalid.  Last iteration for 
        // possibly still valid root1s.  If they are still valid, 
        // record the sieve hit, advance them, and swap with the
        // other root
        _mm512_store_si512(r_id1, vroot1);
        valid_mask_2 = valid_mask_1 & _mm512_cmplt_epu16_mask(vroot1, vblock);
        msk_2 = valid_mask_2;

        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }

        // reduce both roots and store back for the next block
        vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
        vroot1 = _mm512_sub_epi16(vroot1, vblock);
        vroot2 = _mm512_sub_epi16(vroot2, vblock);
        _mm512_storeu_si512(fb->root1 + i, _mm512_min_epu16(vroot1, vroot2));
        _mm512_storeu_si512(fb->root2 + i, _mm512_max_epu16(vroot1, vroot2));
    }
    
    if ((med_B - i) < 32)
    {
        bound = med_B;
    }

    for (; i < bound; i++)
    {
        uint8_t* s2;

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

#elif defined(USE_AVX2) && defined(USE_BMI2)

    bound = full_fb->fb_15bit_B - 16;

    for (i = start_prime; i < bound; i++)
    {
        uint8_t* s2;

        if ((i & 15) == 0)
            break;

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

    __m256i vzero = _mm256_setzero_si256();
    __m256i vblock = _mm256_set1_epi16(BLOCKSIZE);
    __m256i vprime, vroot1, vroot2, vtmp1, vtmp2;
    __m256i valid_mask_1, valid_mask_2, initial_mask;
    ALIGNED_MEM uint16_t r_id1[16];
    ALIGNED_MEM uint16_t r_id2[16];

    for (; i < bound; i += 16)
    {
        uint32_t msk_2;
        int pos;
        int initial_mask_32;

        vprime = _mm256_load_si256((__m256i*)(fb->prime + i));
        vroot1 = _mm256_load_si256((__m256i*)(fb->root1 + i));
        vroot2 = _mm256_load_si256((__m256i*)(fb->root2 + i));
        logp = fb->logp[i]; // approximate the next 32 logp's as equal to this one.

        // we don't sieve primes that are part of the poly
        valid_mask_1 = initial_mask = valid_mask_2 = _mm256_cmpgt_epi16(vprime, vzero);

        // make it so we write to a dummy sieve location for non-sieved primes
        vtmp1 = _mm256_andnot_si256(valid_mask_2, vblock);
        vtmp2 = _mm256_andnot_si256(valid_mask_2, vblock);
        vroot1 = _mm256_add_epi16(vtmp1, vroot1);
        vroot2 = _mm256_add_epi16(vtmp2, vroot2);
        initial_mask_32 = _mm256_movemask_epi8(initial_mask);

        // until things start to drop off the end of the interval, 
        // simply dump in all logs.
        do
        {
            _mm256_store_si256((__m256i*)r_id1, vroot1);
            _mm256_store_si256((__m256i*)r_id2, vroot2);

            sieve[r_id1[0]] -= logp;
            sieve[r_id1[1]] -= logp;
            sieve[r_id1[2]] -= logp;
            sieve[r_id1[3]] -= logp;
            sieve[r_id1[4]] -= logp;
            sieve[r_id1[5]] -= logp;
            sieve[r_id1[6]] -= logp;
            sieve[r_id1[7]] -= logp;
            sieve[r_id1[8]] -= logp;
            sieve[r_id1[9]] -= logp;
            sieve[r_id1[10]] -= logp;
            sieve[r_id1[11]] -= logp;
            sieve[r_id1[12]] -= logp;
            sieve[r_id1[13]] -= logp;
            sieve[r_id1[14]] -= logp;
            sieve[r_id1[15]] -= logp;
            sieve[r_id2[0]] -= logp;
            sieve[r_id2[1]] -= logp;
            sieve[r_id2[2]] -= logp;
            sieve[r_id2[3]] -= logp;
            sieve[r_id2[4]] -= logp;
            sieve[r_id2[5]] -= logp;
            sieve[r_id2[6]] -= logp;
            sieve[r_id2[7]] -= logp;
            sieve[r_id2[8]] -= logp;
            sieve[r_id2[9]] -= logp;
            sieve[r_id2[10]] -= logp;
            sieve[r_id2[11]] -= logp;
            sieve[r_id2[12]] -= logp;
            sieve[r_id2[13]] -= logp;
            sieve[r_id2[14]] -= logp;
            sieve[r_id2[15]] -= logp;

            vroot1 = _mm256_add_epi16(vroot1, vprime);
            vroot2 = _mm256_add_epi16(vroot2, vprime);

            vtmp2 = _mm256_srli_epi16(vroot2, 15);
            vtmp2 = _mm256_or_si256(_mm256_cmpgt_epi16(vtmp2, vzero),
                _mm256_cmpeq_epi16(vroot2, vblock));

            valid_mask_2 = _mm256_andnot_si256(vtmp2, valid_mask_2);
        } while (_mm256_movemask_epi8(valid_mask_2) == initial_mask_32);

        // zero out the primes where roots have exceeded the block
        vprime = _mm256_and_si256(valid_mask_2, vprime);

        // as roots start to exceed the block size, selectively 
        // dump in logs
        while ((msk_2 = _mm256_movemask_epi8(valid_mask_2)) > 0)
        {
            _mm256_store_si256((__m256i*)r_id1, vroot1);
            _mm256_store_si256((__m256i*)r_id2, vroot2);

            msk_2 &= 0xaaaaaaaa;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos >> 1]] -= logp;
                sieve[r_id1[pos >> 1]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }

            vroot1 = _mm256_add_epi16(vroot1, vprime);
            vroot2 = _mm256_add_epi16(vroot2, vprime);

            vtmp2 = _mm256_srli_epi16(vroot2, 15);
            vtmp2 = _mm256_or_si256(_mm256_cmpgt_epi16(vtmp2, vzero),
                _mm256_cmpeq_epi16(vroot2, vblock));

            valid_mask_2 = _mm256_andnot_si256(vtmp2, valid_mask_2);

            // zero out the primes where roots have exceeded the block
            vprime = _mm256_and_si256(valid_mask_2, vprime);
        }

        // now all larger roots are invalid.  Last iteration for 
        // possibly still valid root1s.  If they are still valid, 
        // record the sieve hit, advance them, and swap with the
        // other root
        _mm256_store_si256((__m256i*)r_id1, vroot1);

        vtmp2 = _mm256_srli_epi16(vroot1, 15);
        vtmp2 = _mm256_or_si256(_mm256_cmpgt_epi16(vtmp2, vzero),
            _mm256_cmpeq_epi16(vroot1, vblock));
        valid_mask_2 = _mm256_andnot_si256(vtmp2, valid_mask_1);
        msk_2 = _mm256_movemask_epi8(valid_mask_2);

        msk_2 &= 0xaaaaaaaa;
        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos >> 1]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }

        // reload the primes, then zero out the ones where root1 has
        // already exceeded the block.
        vprime = _mm256_load_si256((__m256i*)(fb->prime + i));
        vprime = _mm256_and_si256(valid_mask_2, vprime);
        vprime = _mm256_and_si256(initial_mask, vprime);

        // reduce both roots and store back for the next block
        vroot1 = _mm256_add_epi16(vroot1, vprime);


        vroot1 = _mm256_sub_epi16(vroot1, vblock);
        vroot2 = _mm256_sub_epi16(vroot2, vblock);
        _mm256_store_si256((__m256i*)(fb->root1 + i), _mm256_min_epu16(vroot1, vroot2));
        _mm256_store_si256((__m256i*)(fb->root2 + i), _mm256_max_epu16(vroot1, vroot2));

    }

    for (; i < bound; i++)
    {
        uint8_t* s2;

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

#else
    for (i = start_prime; i < full_fb->fb_13bit_B - 8; i++)
    {
        uint8_t* s2;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    // the small prime sieve stops just before prime exceeds 2^13
    // the next sse2 sieve assumes primes exceed 2^13.  since
    // some of the primes in the next set of 8 primes could be less
    // than the cutoff and some are greater than, we have to do this
    // small set of crossover primes manually, one at a time.
    for (; i < full_fb->fb_13bit_B; i++)
    {
        uint8_t* s2;

        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];

        //// special exit condition: when prime > 8192 and i % 8 is 0;
        //if ((prime > 8192) && ((i&7) == 0))
        //	break;

        // invalid root (part of poly->a)
        if (prime == 0)
            continue;

        SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

#endif

    

    

#ifdef USE_AVX512BW

    __m512i vp;
    __m512i vr1;
    __m512i vr2;
    __mmask32 result2;
    __mmask32 result1;
    uint32_t res2;
    uint32_t res1;

    for (; i < full_fb->fb_15bit_B - 32; i += 32)
    {
        __m512i vprime, vroot1, vroot2;
        __mmask32 valid_mask_1, valid_mask_2, initial_mask;
        uint32_t msk_2;
        int pos;

        vprime = _mm512_loadu_si512(fb->prime + i);
        vroot1 = _mm512_loadu_si512(fb->root1 + i);
        vroot2 = _mm512_loadu_si512(fb->root2 + i);
        logp = fb->logp[i]; // approximate the next 32 logp's as equal to this one.

        // we don't sieve primes that are part of the poly
        valid_mask_1 = initial_mask = valid_mask_2 = _mm512_cmpgt_epu16_mask(vprime, vzero);

        // as roots start to exceed the block size, selectively 
        // dump in logs
        while (valid_mask_2 > 0)
        {
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            msk_2 = valid_mask_2;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos]] -= logp;
                sieve[r_id1[pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // now all larger roots are invalid.  Last iteration for 
        // possibly still valid root1s.  If they are still valid, 
        // record the sieve hit, advance them, and swap with the
        // other root
        _mm512_store_si512(r_id1, vroot1);
        valid_mask_2 = valid_mask_1 & _mm512_cmplt_epu16_mask(vroot1, vblock);
        msk_2 = valid_mask_2;

        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }

        // reduce both roots and store back for the next block
        vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
        vroot1 = _mm512_sub_epi16(vroot1, vblock);
        vroot2 = _mm512_sub_epi16(vroot2, vblock);
        _mm512_storeu_si512(fb->root1 + i, _mm512_min_epu16(vroot1, vroot2));
        _mm512_storeu_si512(fb->root2 + i, _mm512_max_epu16(vroot1, vroot2));
    }

    // sieve primes 32 at a time, 2^15 < p < med_B
    logp = 15;
    for (; i < med_B - 32; i += 32) {
        //printf("loading from index %d\n", i); fflush(stdout);
        vp = _mm512_loadu_si512((fb->prime + i));
        vr1 = _mm512_loadu_si512((fb->root1 + i));
        vr2 = _mm512_loadu_si512((fb->root2 + i));

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
        _mm512_storeu_si512(fb->root1 + i, _mm512_min_epu16(vr1, vr2));
        _mm512_storeu_si512(fb->root2 + i, _mm512_max_epu16(vr1, vr2));
    }

    for (; i < med_B; i++)
    {
        prime = fb->prime[i];
        root1 = fb->root1[i];
        root2 = fb->root2[i];
        logp = fb->logp[i];
    
        if (prime == 0)
            continue;

        SIEVE_1X;
        SIEVE_LAST;
    
        UPDATE_ROOTS;
    }

#elif defined(USE_AVX2) && defined(USE_BMI2)

    // do this small set of crossover primes manually, one at a time,
    // this time for the 15 bit crossover.
    for (; i < med_B; i++)
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
    _INIT_SSE2_SMALL_PRIME_SIEVE;
    _SSE41_SMALL_PRIME_SIEVE;



#else

    // sieve primes 8 at a time, where 8192 < p < blocksize/3
    _INIT_SSE2_SMALL_PRIME_SIEVE;
    _SSE2_SMALL_PRIME_SIEVE_32k_DIV3;

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
    _INIT_SSE2_SMALL_PRIME_SIEVE;
    _SSE2_SMALL_PRIME_SIEVE_14b;

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


    // sieve primes 8 at a time, 2^14 < p < 2^15
    _INIT_SSE2_SMALL_PRIME_SIEVE;
    _SSE2_SMALL_PRIME_SIEVE_15b;


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
    _INIT_SSE2_SMALL_PRIME_SIEVE;
    _SSE41_SMALL_PRIME_SIEVE;


#endif

    return;

}



#endif
