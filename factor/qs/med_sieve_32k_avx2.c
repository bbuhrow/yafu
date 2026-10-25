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
#include "sieve_macros_32k_sse4.1.h"

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


#define _FINALIZE_SORT_UPDATE_AVX2 \
	"vpaddw	    %%xmm4, %%xmm1, %%xmm1 \n\t"			/* r1 = r1 + (p - b) */ \
	"vpaddw	    %%xmm5, %%xmm2, %%xmm2 \n\t"			/* r2 = r2 + (p - b) */ \
	"vpmaxsw	%%xmm1, %%xmm2, %%xmm5 \n\t"			/* replace xmm2 with max of root1 and root2 */ \
	"vpminsw	%%xmm2, %%xmm1, %%xmm4 \n\t"			/* replace xmm1 with min of root1 and root2 */ \
	"vmovdqa    %%xmm4, (%2) \n\t"				/* write new root1's */ \
	"vmovdqa    %%xmm5, (%3) \n\t"				/* write new root2's */

#define _INIT_AVX2_SMALL_PRIME_SIEVE \
	ASM_G (	\
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


#define _AVX2_SMALL_PRIME_SIEVE \
	for (; i < med_B; i += 8) \
		ASM_G (			\
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


#endif // USE_AVX2


#ifdef USE_AVX2

void med_sieveblock_32k_avx2(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    uint32_t i;
    uint32_t med_B;

    uint32_t prime, root1, root2, tmp, stop;
    uint8_t logp;

    uint32_t bound = full_fb->fb_13bit_B - 8;

    med_B = full_fb->med_B;

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
#if (defined(_MSC_VER) && defined(__clang__)) || (defined(__MINGW32__))
    // clang-cl defines _mm256_set1_epi16 as taking a short, and therefore sets 32768 to -32768.
    // so we have to do this instead:
    __m256i vblock = _mm256_set1_epi16(32767);
    vblock = _mm256_add_epi16(vblock, _mm256_set1_epi16(1));
#else
    __m256i vblock = _mm256_set1_epi16(BLOCKSIZE);
#endif
    //__m256i vinterval = _mm256_set1_epi16(4 * BLOCKSIZE);
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
        printf("primes that hit block @ i = %u: ", i);
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
        printf("primes that hit block @ i = %u: ", i);
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
        printf("primes where root1 hits block @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", r_id2[j]);
        }
        printf("\n");
#endif


        vroot1 = _mm256_sub_epi16(vroot1, vblock);
        vroot2 = _mm256_sub_epi16(vroot2, vblock);
        _mm256_store_si256((__m256i*)(fb->root1 + i), _mm256_min_epu16(vroot1, vroot2));
        _mm256_store_si256((__m256i*)(fb->root2 + i), _mm256_max_epu16(vroot1, vroot2));

#ifdef TEST_SIEVE
        printf("after storage\n");
        printf("root1s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", fb->root1[i + j]);
        }
        printf("\n");
        printf("root2s @ i = %u: ", i);
        for (j = 0; j < 16; j++)
        {
            printf("%04x, ", fb->root1[i + j]);
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

        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

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
#if defined(GCC_ASM64X)
    _INIT_AVX2_SMALL_PRIME_SIEVE;
    _AVX2_SMALL_PRIME_SIEVE;
#else
    _INIT_SSE2_SMALL_PRIME_SIEVE;
    _SSE41_SMALL_PRIME_SIEVE;
#endif

    CLEAN_AVX2;

    return;

}

#endif

// intrinsics and sieve routines for processors with AVX512BW
#ifdef USE_AVX512BW

// sink-strip stores (see Sink-Stripe Sieve Stores.md): the selective loops
// write every lane, blending missed lanes to a rotating strip instead of
// looping over the set bits of the hit mask.  The strip is the last
// MED_SINK_BYTES of the previous block, sieve[-MED_SINK_BYTES, 0), so the
// RMWs index from sbase = sieve - MED_SINK_BYTES and hits are biased by
// MED_SINK_BYTES.  That block is finished (sieved, scanned, tdiv'd) before
// this one is sieved, in both the single-block window (where it is the
// leading guard) and the whole-interval direct layout (where pos block 0
// uses the leading guard and neg block 0 uses the last pos block), so
// nothing needs saving or restoring.
// bit 0: 13-15 bit selective loops, bit 1: 15 bit < p < med_B loop
#ifndef MED_SINK_STRIP
#define MED_SINK_STRIP 3
#endif

// power of 2, >= 64.  each use takes 64 bytes (32 lanes x 2 roots)
#ifndef MED_SINK_BYTES
#define MED_SINK_BYTES 256
#endif

#define MED_RMW4(r, k) \
    sbase[r[k]] -= logp; sbase[r[k + 1]] -= logp; \
    sbase[r[k + 2]] -= logp; sbase[r[k + 3]] -= logp;
#define MED_RMW32(r) \
    MED_RMW4(r, 0) MED_RMW4(r, 4) MED_RMW4(r, 8) MED_RMW4(r, 12) \
    MED_RMW4(r, 16) MED_RMW4(r, 20) MED_RMW4(r, 24) MED_RMW4(r, 28)

// sbase-relative index: hit lanes (root < 32768) get root + MED_SINK_BYTES,
// missed lanes get their strip offset
#define MED_SINK_IDX(m, vs, vr) \
    _mm512_mask_blend_epi16(m, vs, _mm512_add_epi16(vr, vsink_bias))

#define MED_SINK_NEXT \
    vsink1 = _mm512_and_si512(_mm512_add_epi16(vsink1, vsink_step), vsink_wrap); \
    vsink2 = _mm512_add_epi16(vsink1, vsink_half);

static const uint16_t med_sink_lanes[32] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 };

void med_sieveblock_32k_avx512bw(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    uint32_t i;
    uint32_t med_B;

    uint32_t prime, root1, root2, tmp, stop;
    uint8_t logp;
#if defined(_MSC_VER) && defined(__clang__)
    __m512i vblock = _mm512_set1_epi16(32767);
    vblock = _mm512_add_epi16(vblock, _mm512_set1_epi16(1));
#else
    __m512i vblock = _mm512_set1_epi16(BLOCKSIZE);
#endif
    __m512i vzero = _mm512_setzero_epi32();

    ALIGNED_MEM uint16_t r_id1[32];
    ALIGNED_MEM uint16_t r_id2[32];

#if MED_SINK_STRIP
    uint8_t* sbase = sieve - MED_SINK_BYTES;
    __m512i vsink_bias = _mm512_set1_epi16(MED_SINK_BYTES);
    __m512i vsink1 = _mm512_loadu_si512(med_sink_lanes);
    __m512i vsink_half = _mm512_set1_epi16(32);
    __m512i vsink2 = _mm512_add_epi16(vsink1, vsink_half);
    __m512i vsink_step = _mm512_set1_epi16(64);
    __m512i vsink_wrap = _mm512_set1_epi16(MED_SINK_BYTES - 1);
#endif

    uint32_t bound = full_fb->fb_15bit_B - 32;

    if (full_fb->med_B > (full_fb->fb_15bit_B + 32))
        bound += 32;

    med_B = full_fb->med_B;

    //printf("start: %d, stop: %d\n", start_prime, bound);

    for (i = start_prime; i < bound; i++)
    {
        uint8_t* s2;

        if ((i & 31) == 0)
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

    //printf("start: %d, stop: %d\n", i, bound);

    for ( ; i < bound; i += 32)
    {
        __m512i vprime, vroot1, vroot2;
        __mmask32 valid_mask_1, valid_mask_2, initial_mask;
        uint32_t msk_2;
        int pos;

        vprime = _mm512_load_si512(fb->prime + i);
        vroot1 = _mm512_load_si512(fb->root1 + i);
        vroot2 = _mm512_load_si512(fb->root2 + i);
        logp = fb->logp[i+16]; // approximate the next 32 logp's as equal to the midpoint

        // we don't sieve primes that are part of the poly
        valid_mask_1 = initial_mask = valid_mask_2 = _mm512_cmpgt_epu16_mask(vprime, vzero);

        // make it so we write to a dummy sieve location for non-sieved primes
        // NOTE: this doesn't work in subset-sum variations where we sieve 
        // a whole interval - end up overwriting things in subsequent blocks.
        // plus the interval could be larger than 16 bits.  we could maybe save and
        // restore sieve indices that are impacted by this.
        //vroot1 = _mm512_mask_add_epi16(vroot1, ~initial_mask, vroot1, vinterval); 
        //vroot2 = _mm512_mask_add_epi16(vroot2, ~initial_mask, vroot2, vinterval); 
        vroot1 = _mm512_mask_set1_epi16(vroot1, ~initial_mask, 0); 
        vroot2 = _mm512_mask_set1_epi16(vroot2, ~initial_mask, 0); 

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
#if MED_SINK_STRIP & 1
            _mm512_store_si512(r_id1, MED_SINK_IDX(valid_mask_2, vsink1, vroot1));
            _mm512_store_si512(r_id2, MED_SINK_IDX(valid_mask_2, vsink2, vroot2));
            MED_RMW32(r_id1);
            MED_RMW32(r_id2);
            MED_SINK_NEXT;
#else
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            msk_2 = valid_mask_2;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos]] -= logp;
                sieve[r_id1[pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }
#endif

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // now all larger roots are invalid.  Last iteration for
        // possibly still valid root1s.  If they are still valid,
        // record the sieve hit, advance them, and swap with the
        // other root
        valid_mask_2 = valid_mask_1 & _mm512_cmplt_epu16_mask(vroot1, vblock);
#if MED_SINK_STRIP & 1
        _mm512_store_si512(r_id1, MED_SINK_IDX(valid_mask_2, vsink1, vroot1));
        MED_RMW32(r_id1);
        MED_SINK_NEXT;
#else
        _mm512_store_si512(r_id1, vroot1);
        msk_2 = valid_mask_2;

        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }
#endif

        // reduce both roots and store back for the next block
        vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
        vroot1 = _mm512_sub_epi16(vroot1, vblock);
        vroot2 = _mm512_sub_epi16(vroot2, vblock);
        _mm512_store_si512(fb->root1 + i, _mm512_min_epu16(vroot1, vroot2));
        _mm512_store_si512(fb->root2 + i, _mm512_max_epu16(vroot1, vroot2));
    }

    if ((med_B - i) < 32)
    {
        bound = med_B;
    }

    //printf("start: %d, stop: %d\n", i, bound);

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

        // medB can be greater than 32k, so can't do sieve2x
        //SIEVE_2X;
        SIEVE_1X;
        SIEVE_LAST;
        UPDATE_ROOTS;
    }

    __m512i vp;
    __m512i vr1;
    __m512i vr2;
    __mmask32 result2;
    __mmask32 result1;
    uint32_t res2;
    uint32_t res1;

    //printf("start: %d, stop: %d\n", i, full_fb->fb_15bit_B - 32);

    for (; i < full_fb->fb_15bit_B - 32; i += 32)
    {
        __m512i vprime, vroot1, vroot2;
        __mmask32 valid_mask_1, valid_mask_2, initial_mask;
        uint32_t msk_2;
        int pos;

        vprime = _mm512_load_si512(fb->prime + i);
        vroot1 = _mm512_load_si512(fb->root1 + i);
        vroot2 = _mm512_load_si512(fb->root2 + i);
        logp = fb->logp[i]; // approximate the next 32 logp's as equal to this one.

        // we don't sieve primes that are part of the poly
        valid_mask_1 = initial_mask = valid_mask_2 = _mm512_cmpgt_epu16_mask(vprime, vzero);

        // as roots start to exceed the block size, selectively 
        // dump in logs
        while (valid_mask_2 > 0)
        {
#if MED_SINK_STRIP & 1
            _mm512_store_si512(r_id1, MED_SINK_IDX(valid_mask_2, vsink1, vroot1));
            _mm512_store_si512(r_id2, MED_SINK_IDX(valid_mask_2, vsink2, vroot2));
            MED_RMW32(r_id1);
            MED_RMW32(r_id2);
            MED_SINK_NEXT;
#else
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            msk_2 = valid_mask_2;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos]] -= logp;
                sieve[r_id1[pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }
#endif

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // now all larger roots are invalid.  Last iteration for
        // possibly still valid root1s.  If they are still valid,
        // record the sieve hit, advance them, and swap with the
        // other root
        valid_mask_2 = valid_mask_1 & _mm512_cmplt_epu16_mask(vroot1, vblock);
#if MED_SINK_STRIP & 1
        _mm512_store_si512(r_id1, MED_SINK_IDX(valid_mask_2, vsink1, vroot1));
        MED_RMW32(r_id1);
        MED_SINK_NEXT;
#else
        _mm512_store_si512(r_id1, vroot1);
        msk_2 = valid_mask_2;

        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }
#endif

        // reduce both roots and store back for the next block
        vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
        vroot1 = _mm512_sub_epi16(vroot1, vblock);
        vroot2 = _mm512_sub_epi16(vroot2, vblock);
        _mm512_store_si512(fb->root1 + i, _mm512_min_epu16(vroot1, vroot2));
        _mm512_store_si512(fb->root2 + i, _mm512_max_epu16(vroot1, vroot2));
    }

    //printf("start: %d, stop: %d\n", i, med_B - 32);
    // sieve primes 32 at a time, 2^15 < p < med_B
    logp = 15;
    for (; i < med_B - 32; i += 32) {
        //printf("loading from index %d\n", i); fflush(stdout);
        vp = _mm512_load_si512((fb->prime + i));
        vr1 = _mm512_load_si512((fb->root1 + i));
        vr2 = _mm512_load_si512((fb->root2 + i));

        result2 = _mm512_cmp_epu16_mask(vr2, vblock, _MM_CMPINT_LT);
        result1 = _mm512_cmp_epu16_mask(vr1, vblock, _MM_CMPINT_LT);

#if MED_SINK_STRIP & 2
        // root1 <= root2, so result2 is a subset of result1
        _mm512_store_si512(r_id1, MED_SINK_IDX(result1, vsink1, vr1));
        _mm512_store_si512(r_id2, MED_SINK_IDX(result2, vsink2, vr2));
        MED_RMW32(r_id1);
        MED_RMW32(r_id2);
        MED_SINK_NEXT;
#else
        res2 = result2;

        while (res2 > 0) {
            int idx = _trail_zcnt(res2);
            sieve[fb->root2[i + idx]] -= logp;
            sieve[fb->root1[i + idx]] -= logp;
            res2 = _reset_lsb(res2);
        }

        // res1 will have fewer set bits this way, so we
        // have fewer overall loop iterations
        res1 = result1 & (~result2);

        while (res1 > 0) {
            int idx = _trail_zcnt(res1);
            sieve[fb->root1[i + idx]] -= logp;
            res1 = _reset_lsb(res1);
        }
#endif

        vr1 = _mm512_mask_add_epi16(vr1, result1, vr1, vp);
        vr2 = _mm512_mask_add_epi16(vr2, result2, vr2, vp);
        vr1 = _mm512_sub_epi16(vr1, vblock);
        vr2 = _mm512_sub_epi16(vr2, vblock);
        _mm512_store_si512(fb->root1 + i, _mm512_min_epu16(vr1, vr2));
        _mm512_store_si512(fb->root2 + i, _mm512_max_epu16(vr1, vr2));
    }

    //printf("start: %d, stop: %d\n", i, med_B);
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

    //exit(1);

#ifdef USE_SS_SEARCH
    // restore sieve locations associated with roots we didn't sieve
    sieve[0] = 0x7f;
#endif

    return;

}

#define MED_DUMP4(r, k) \
    sieve[r[k]] -= logp; sieve[r[k + 1]] -= logp; \
    sieve[r[k + 2]] -= logp; sieve[r[k + 3]] -= logp;
#define MED_DUMP32(r) \
    MED_DUMP4(r, 0) MED_DUMP4(r, 4) MED_DUMP4(r, 8) MED_DUMP4(r, 12) \
    MED_DUMP4(r, 16) MED_DUMP4(r, 20) MED_DUMP4(r, 24) MED_DUMP4(r, 28)

// one fixed-step pass over a 32-prime group: lanes in mask mv with a root
// in the block (r < 32768) take a hit at sieve[r] and advance by p, the rest
// write the sink strip (sbase is sieve - MED_SINK_BYTES, the previous
// block's tail, see MED_SINK_STRIP).
#define MED_FIXED_PASS(mv) \
    m1 = (mv) & _mm512_cmplt_epu16_mask(vroot1, vblock); \
    m2 = (mv) & _mm512_cmplt_epu16_mask(vroot2, vblock); \
    _mm512_store_si512(r_id1, MED_SINK_IDX(m1, vsink1, vroot1)); \
    _mm512_store_si512(r_id2, MED_SINK_IDX(m2, vsink2, vroot2)); \
    MED_RMW32(r_id1); \
    MED_RMW32(r_id2); \
    MED_SINK_NEXT; \
    vroot1 = _mm512_mask_add_epi16(vroot1, m1, vroot1, vprime); \
    vroot2 = _mm512_mask_add_epi16(vroot2, m2, vroot2, vprime);

// sink-strip stores in the selective loops below fb_32k_div3 too.  Off by
// default: they remove the branch misses but measured slightly slower here
#ifndef MED_SPLIT_SINK_SELECTIVE
#define MED_SPLIT_SINK_SELECTIVE 0
#endif

// start the vector loop at the aligned group holding start_prime instead of
// sieving up to the first aligned index with scalar code.  Off: the smallest
// primes span too wide a range per 32-lane group, so most of their
// iterations end up in the selective bit loop (measured 3-25% slower)
#ifndef MED_SPLIT_VECTOR_HEAD
#define MED_SPLIT_VECTOR_HEAD 0
#endif

// medium-prime sieve with fixed trip counts.  Below fb_32k_div3 this is the
// avx512bw loop of med_sieveblock_32k_avx512bw.  Above it a root hits at
// most K = ceil(32768 / p) times per block, so each 32-prime group runs K
// MED_FIXED_PASSes, with K from the group's smallest prime:
//   10922 < p < 16384: K = 3 (2 certain hits + at most 1)
//   16384 < p < 32768: K = 2 (1 certain hit + at most 1)
//   32768 < p < med_B: K = 1 (at most 1)
// no loop count depends on the data.  r + p can pass 2^16 only for
// p > 32768, so those lanes (which have at most one hit) sit out passes
// after the first, where a wrapped root would look like another hit.
// a-primes have prime = root = 0 and so only hit sieve[0], as in the
// vector loop.
void med_sieveblock_32k_avx512bw_split(uint8_t* sieve, sieve_fb_compressed* fb, fb_list* full_fb,
    uint32_t start_prime, uint8_t s_init)
{
    uint32_t i;
    uint32_t med_B = full_fb->med_B;
    uint32_t* plist = full_fb->list->prime;
    uint32_t prime, root1, root2, tmp, stop;
    uint32_t b3, vbound;
    uint8_t logp;
    uint8_t* sbase = sieve - MED_SINK_BYTES;
    uint16_t* fprime = fb->prime;
    uint16_t* froot1 = fb->root1;
    uint16_t* froot2 = fb->root2;
#if defined(_MSC_VER) && defined(__clang__)
    __m512i vblock = _mm512_set1_epi16(32767);
    vblock = _mm512_add_epi16(vblock, _mm512_set1_epi16(1));
#else
    __m512i vblock = _mm512_set1_epi16(BLOCKSIZE);
#endif
    __m512i vzero = _mm512_setzero_epi32();
    __m512i vffff = _mm512_set1_epi16((short)0xffff);
    __m512i vsink_bias = _mm512_set1_epi16(MED_SINK_BYTES);
    __m512i vsink1 = _mm512_loadu_si512(med_sink_lanes);
    __m512i vsink_half = _mm512_set1_epi16(32);
    __m512i vsink2 = _mm512_add_epi16(vsink1, vsink_half);
    __m512i vsink_step = _mm512_set1_epi16(64);
    __m512i vsink_wrap = _mm512_set1_epi16(MED_SINK_BYTES - 1);

    ALIGNED_MEM uint16_t r_id1[32];
    ALIGNED_MEM uint16_t r_id2[32];

    b3 = (full_fb->fb_32k_div3 < med_B) ? full_fb->fb_32k_div3 : med_B;
    vbound = b3 & ~31u;

#if MED_SPLIT_VECTOR_HEAD
    // the vector loop starts at the aligned group holding start_prime, with
    // the lanes below it masked off and their roots left alone
    i = start_prime & ~31u;
    if (i >= vbound)
        i = start_prime;
#else
    i = start_prime;
#endif
    for (; i < b3; i++)
    {
        uint8_t* s2;

        if ((i & 31) == 0)
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

    for (; i < vbound; i += 32)
    {
        __m512i vprime, vroot1, vroot2;
        __mmask32 valid_mask_1, valid_mask_2, initial_mask;
        __mmask32 keep_mask;
#if !MED_SPLIT_SINK_SELECTIVE
        uint32_t msk_2;
        int pos;
#endif

        vprime = _mm512_load_si512(fb->prime + i);
        vroot1 = _mm512_load_si512(fb->root1 + i);
        vroot2 = _mm512_load_si512(fb->root2 + i);
        logp = fb->logp[i + 16]; // approximate the next 32 logp's as equal to the midpoint

        // we don't sieve primes that are part of the poly (prime 0, or roots
        // marked 0xffff for a-primes <= 8192) or lanes below start_prime
        // (MED_SPLIT_VECTOR_HEAD).  Their roots are parked at 0xffff, so the
        // unconditional dumps write sieve[65535] (guard space, or saved and
        // restored by process_poly in the direct layout) and, unlike a
        // 0xffff root plus p, never wrap into the block.
        keep_mask = (i < start_prime) ? ~(__mmask32)_bzhi_u32(0xffffffff, start_prime - i) :
            0xffffffff;
        valid_mask_1 = initial_mask = valid_mask_2 = keep_mask &
            _mm512_cmpgt_epu16_mask(vprime, vzero) & _mm512_cmpneq_epu16_mask(vroot1, vffff);
        vroot1 = _mm512_mask_set1_epi16(vroot1, ~initial_mask, 0xffff);
        vroot2 = _mm512_mask_set1_epi16(vroot2, ~initial_mask, 0xffff);

        // until things start to drop off the end of the interval,
        // simply dump in all logs.
        while (valid_mask_2 == initial_mask)
        {
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);
            MED_DUMP32(r_id1);
            MED_DUMP32(r_id2);

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // as roots start to exceed the block size, selectively
        // dump in logs
        while (valid_mask_2 > 0)
        {
#if MED_SPLIT_SINK_SELECTIVE
            _mm512_store_si512(r_id1, MED_SINK_IDX(valid_mask_2, vsink1, vroot1));
            _mm512_store_si512(r_id2, MED_SINK_IDX(valid_mask_2, vsink2, vroot2));
            MED_RMW32(r_id1);
            MED_RMW32(r_id2);
            MED_SINK_NEXT;
#else
            _mm512_store_si512(r_id1, vroot1);
            _mm512_store_si512(r_id2, vroot2);

            msk_2 = valid_mask_2;
            while (msk_2 > 0) {
                pos = _trail_zcnt(msk_2);
                sieve[r_id2[pos]] -= logp;
                sieve[r_id1[pos]] -= logp;
                msk_2 = _reset_lsb(msk_2);
            }
#endif

            vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, valid_mask_2, vroot2, vprime);

            valid_mask_2 &= _mm512_cmplt_epu16_mask(vroot2, vblock);
        }

        // last iteration for possibly still valid root1s
        valid_mask_2 = valid_mask_1 & _mm512_cmplt_epu16_mask(vroot1, vblock);
#if MED_SPLIT_SINK_SELECTIVE
        _mm512_store_si512(r_id1, MED_SINK_IDX(valid_mask_2, vsink1, vroot1));
        MED_RMW32(r_id1);
        MED_SINK_NEXT;
#else
        _mm512_store_si512(r_id1, vroot1);
        msk_2 = valid_mask_2;

        while (msk_2 > 0) {
            pos = _trail_zcnt(msk_2);
            sieve[r_id1[pos]] -= logp;
            msk_2 = _reset_lsb(msk_2);
        }
#endif

        // reduce both roots and store back for the next block
        vroot1 = _mm512_mask_add_epi16(vroot1, valid_mask_2, vroot1, vprime);
        vroot1 = _mm512_sub_epi16(vroot1, vblock);
        vroot2 = _mm512_sub_epi16(vroot2, vblock);
        _mm512_mask_storeu_epi16(fb->root1 + i, keep_mask, _mm512_min_epu16(vroot1, vroot2));
        _mm512_mask_storeu_epi16(fb->root2 + i, keep_mask, _mm512_max_epu16(vroot1, vroot2));
    }

    for (; i < b3; i++)
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

    for (; i < med_B; i += 32)
    {
        __m512i vprime, vroot1, vroot2;
        __mmask32 lm, nbig, m1, m2;

        lm = ((med_B - i) >= 32) ? 0xffffffff : (__mmask32)_bzhi_u32(0xffffffff, med_B - i);
        vprime = _mm512_maskz_loadu_epi16(lm, fprime + i);
        vroot1 = _mm512_maskz_loadu_epi16(lm, froot1 + i);
        vroot2 = _mm512_maskz_loadu_epi16(lm, froot2 + i);
        logp = (uint8_t)fb->logp[i];
        nbig = lm & _mm512_cmple_epu16_mask(vprime, vblock);

        MED_FIXED_PASS(lm);
        if (plist[i] < 32768)
        {
            MED_FIXED_PASS(nbig);
            if (plist[i] < 16384)
            {
                MED_FIXED_PASS(nbig);
            }
        }

        vroot1 = _mm512_sub_epi16(vroot1, vblock);
        vroot2 = _mm512_sub_epi16(vroot2, vblock);
        _mm512_mask_storeu_epi16(froot1 + i, lm, _mm512_min_epu16(vroot1, vroot2));
        _mm512_mask_storeu_epi16(froot2 + i, lm, _mm512_max_epu16(vroot1, vroot2));
    }

#ifdef USE_SS_SEARCH
    // restore sieve locations associated with roots we didn't sieve
    sieve[0] = 0x7f;
#endif

    return;
}

#endif
