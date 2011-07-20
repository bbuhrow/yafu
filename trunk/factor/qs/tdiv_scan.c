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
#include "factor.h"
#include "util.h"
#include "common.h"

//#define SIQSDEBUG 1

/*
We are given an array of bytes that has been sieved.  The basic trial 
division strategy is as follows:

1) Scan through the array and 'mark' locations that meet criteria 
indicating they may factor completely over the factor base.  

2) 'Filter' the marked locations by trial dividing by small primes
that we did not sieve.  These primes are all less than 256.  If after
removing small primes the location does not meet another set of criteria,
remove it from the 'marked' list (do not subject it to further trial
division).

3) Divide out primes from the factor base between 256 and 2^13 or 2^14, 
depending on the version (2^13 for 32k version, 2^14 for 64k).  

4) Resieve primes between 2^{13|14} and 2^16, max.  

5) Primes larger than 2^16 will have been bucket sieved.  Remove these
by scanning the buckets for sieve hits equal to the current block location.

6) If applicable/appropriate, factor a remaining composite with squfof

this file contains code implementing 1)


*/


#if defined(GCC_ASM32X) || defined(GCC_ASM64X) || defined(__MINGW32__)
	//these compilers support SIMD 
	#define SIMD_SIEVE_SCAN 1
	#define SCAN_CLEAN asm volatile("emms");	

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			asm volatile (		\
				"movdqa (%1), %%xmm0   \n\t"		\
				"orpd 16(%1), %%xmm0    \n\t"		\
				"pmovmskb %%xmm0, %0   \n\t"		\
				: "=r"(result)						\
				: "r"(sieveblock + j), "0"(result)	\
				: "%xmm0");

		#define SIEVE_SCAN_64		\
			asm volatile (							\
				"movdqa (%1), %%xmm0   \n\t"		\
				"orpd 16(%1), %%xmm0    \n\t"		\
				"orpd 32(%1), %%xmm0    \n\t"		\
				"orpd 48(%1), %%xmm0    \n\t"		\
				"pmovmskb %%xmm0, %0   \n\t"		\
				: "=r"(result)						\
				: "r"(sieveblock + j), "0"(result)	\
				: "%xmm0");

		#define SIEVE_SCAN_128		\
			asm volatile (			\
				"movdqa (%1), %%xmm0   \n\t"		\
				"orpd 16(%1), %%xmm0    \n\t"		\
				"orpd 32(%1), %%xmm0    \n\t"		\
				"orpd 48(%1), %%xmm0    \n\t"		\
				"orpd 64(%1), %%xmm0    \n\t"		\
				"orpd 80(%1), %%xmm0    \n\t"		\
				"orpd 96(%1), %%xmm0    \n\t"		\
				"orpd 112(%1), %%xmm0    \n\t"		\
				"pmovmskb %%xmm0, %0   \n\t"		\
				: "=r"(result)						\
				: "r"(sieveblock + j), "0"(result)	\
				: "%xmm0");

	#elif defined(HAS_MMX)
		#define SIEVE_SCAN_32		\
			asm volatile (			\
				"movq (%1), %%mm0     \n\t"		\
				"por 8(%1), %%mm0     \n\t"		\
				"por 16(%1), %%mm0    \n\t"		\
				"por 24(%1), %%mm0    \n\t"		\
				"pmovmskb %%mm0, %0   \n\t"		\
				: "=r"(result)					\
				: "r"(sieveblock + j), "0"(result)	\
				: "%mm0");	

		#define SIEVE_SCAN_64		\
			asm volatile (			\
				"movq (%1), %%mm0     \n\t"		\
				"por 8(%1), %%mm0     \n\t"		\
				"por 16(%1), %%mm0    \n\t"		\
				"por 24(%1), %%mm0    \n\t"		\
				"por 32(%1), %%mm0    \n\t"		\
				"por 40(%1), %%mm0    \n\t"		\
				"por 48(%1), %%mm0    \n\t"		\
				"por 56(%1), %%mm0    \n\t"		\
				"pmovmskb %%mm0, %0   \n\t"		\
				: "=r"(result)					\
				: "r"(sieveblock + j), "0"(result)	\
				: "%mm0");

		#define SIEVE_SCAN_128		\
			asm volatile (				\
				"movq (%1), %%mm0     \n\t"		\
				"por 8(%1), %%mm0     \n\t"		\
				"por 16(%1), %%mm0    \n\t"		\
				"por 24(%1), %%mm0    \n\t"		\
				"por 32(%1), %%mm0    \n\t"		\
				"por 40(%1), %%mm0    \n\t"		\
				"por 48(%1), %%mm0    \n\t"		\
				"por 56(%1), %%mm0    \n\t"		\
				"por 64(%1), %%mm0     \n\t"	\
				"por 72(%1), %%mm0    \n\t"		\
				"por 80(%1), %%mm0    \n\t"		\
				"por 88(%1), %%mm0    \n\t"		\
				"por 96(%1), %%mm0    \n\t"		\
				"por 104(%1), %%mm0    \n\t"	\
				"por 112(%1), %%mm0    \n\t"	\
				"por 120(%1), %%mm0		\n\t"	\
				"pmovmskb %%mm0, %0   \n\t"		\
				: "=r"(result)					\
				: "r"(sieveblock + j), "0"(result)	\
				: "%mm0");

	#else
		#undef SIMD_SIEVE_SCAN
	#endif

#elif defined(MSC_ASM32A)
	#define SIMD_SIEVE_SCAN 1
	#define SCAN_CLEAN ASM_M {emms};

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
					ASM_M pmovmskb ecx, xmm0			\
					ASM_M mov result, ecx};			\
			} while (0);


		#define SIEVE_SCAN_64	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 32]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 48]	\
					ASM_M pmovmskb ecx, xmm0			\
					ASM_M mov result, ecx};			\
			} while (0);

		#define SIEVE_SCAN_128	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 32]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 48]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 64]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 80]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 96]	\
					ASM_M por xmm0, XMMWORD PTR [edi + 112]	\
					ASM_M pmovmskb ecx, xmm0			\
					ASM_M mov result, ecx};			\
			} while (0);

	#elif defined(HAS_MMX)

		#define SIEVE_SCAN_32	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movq mm0, QWORD PTR [edi]	\
					ASM_M por mm0, QWORD PTR [edi + 8]	\
					ASM_M por mm0, QWORD PTR [edi + 16]	\
					ASM_M por mm0, QWORD PTR [edi + 24]	\
					ASM_M pmovmskb ecx, mm0			\
					ASM_M mov result, ecx};			\
			} while (0);

		#define SIEVE_SCAN_64	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movq mm0, QWORD PTR [edi]	\
					ASM_M por mm0, QWORD PTR [edi + 8]	\
					ASM_M por mm0, QWORD PTR [edi + 16]	\
					ASM_M por mm0, QWORD PTR [edi + 24]	\
					ASM_M por mm0, QWORD PTR [edi + 32]	\
					ASM_M por mm0, QWORD PTR [edi + 40]	\
					ASM_M por mm0, QWORD PTR [edi + 48]	\
					ASM_M por mm0, QWORD PTR [edi + 56]	\
					ASM_M pmovmskb ecx, mm0			\
					ASM_M mov result, ecx};			\
			} while (0);

		#define SIEVE_SCAN_128	\
			do	{						\
				uint64 *localblock = sieveblock + j;	\
				ASM_M  {			\
					ASM_M mov edi, localblock			\
					ASM_M movq mm0, QWORD PTR [edi]	\
					ASM_M por mm0, QWORD PTR [edi + 8]	\
					ASM_M por mm0, QWORD PTR [edi + 16]	\
					ASM_M por mm0, QWORD PTR [edi + 24]	\
					ASM_M por mm0, QWORD PTR [edi + 32]	\
					ASM_M por mm0, QWORD PTR [edi + 40]	\
					ASM_M por mm0, QWORD PTR [edi + 48]	\
					ASM_M por mm0, QWORD PTR [edi + 56]	\
					ASM_M por mm0, QWORD PTR [edi + 64]	\
					ASM_M por mm0, QWORD PTR [edi + 72]	\
					ASM_M por mm0, QWORD PTR [edi + 80]	\
					ASM_M por mm0, QWORD PTR [edi + 88]	\
					ASM_M por mm0, QWORD PTR [edi + 96]	\
					ASM_M por mm0, QWORD PTR [edi + 104]	\
					ASM_M por mm0, QWORD PTR [edi + 112]	\
					ASM_M por mm0, QWORD PTR [edi + 120]	\
					ASM_M pmovmskb ecx, mm0			\
					ASM_M mov result, ecx};			\
			} while (0);

	#else
		#undef SIMD_SIEVE_SCAN

	#endif

#elif defined(_WIN64)

	#define SIMD_SIEVE_SCAN 1
	#define SCAN_CLEAN /*nothing*/

	#if defined(HAS_SSE2)
		//top level sieve scanning with SSE2
		#define SIEVE_SCAN_32	\
			do	{						\
				__m128i local_block;	\
				__m128i local_block2;	\
				local_block = _mm_load_si128(sieveblock + j); \
				local_block2 = _mm_load_si128(sieveblock + j + 2); \
				local_block = _mm_or_si128(local_block, local_block2); \
				result = _mm_movemask_epi8(local_block); \
			} while (0);


		#define SIEVE_SCAN_64	\
			do	{				  		\
				__m128i local_block;	\
				__m128i local_block2;	\
				__m128i local_block3;	\
				__m128i local_block4;	\
				local_block = _mm_load_si128(sieveblock + j); \
				local_block2 = _mm_load_si128(sieveblock + j + 2); \
				local_block3 = _mm_load_si128(sieveblock + j + 4); \
				local_block = _mm_or_si128(local_block, local_block2); \
				local_block = _mm_or_si128(local_block, local_block3); \
				local_block4 = _mm_load_si128(sieveblock + j + 6); \
				local_block = _mm_or_si128(local_block, local_block4); \
				result = _mm_movemask_epi8(local_block); \
			} while (0);

		#define SIEVE_SCAN_128	\
			do	{						\
				__m128i local_block;	\
				__m128i local_block2;	\
				__m128i local_block3;	\
				__m128i local_block4;	\
				__m128i local_block5;	\
				__m128i local_block6;	\
				__m128i local_block7;	\
				__m128i local_block8;	\
				local_block = _mm_load_si128(sieveblock + j); \
				local_block2 = _mm_load_si128(sieveblock + j + 2); \
				local_block3 = _mm_load_si128(sieveblock + j + 4); \
				local_block = _mm_or_si128(local_block, local_block2); \
				local_block4 = _mm_load_si128(sieveblock + j + 6); \
				local_block = _mm_or_si128(local_block, local_block3); \
				local_block5 = _mm_load_si128(sieveblock + j + 8); \
				local_block = _mm_or_si128(local_block, local_block4); \
				local_block6 = _mm_load_si128(sieveblock + j + 10); \
				local_block = _mm_or_si128(local_block, local_block5); \
				local_block7 = _mm_load_si128(sieveblock + j + 12); \
				local_block = _mm_or_si128(local_block, local_block6); \
				local_block8 = _mm_load_si128(sieveblock + j + 14); \
				local_block = _mm_or_si128(local_block, local_block7); \
				local_block = _mm_or_si128(local_block, local_block8); \
				result = _mm_movemask_epi8(local_block); \
			} while (0);

	#else
		#undef SIMD_SIEVE_SCAN
	#endif

#else	/* compiler not recognized*/

	#define SCAN_CLEAN /*nothing*/

#endif

#define SCAN_MASK 0x8080808080808080ULL

	//when we compress small primes into 16 bits of a 32 bit field, the
	//trick of fooling the sieve routine to not sieve those roots which
	//divide poly_a fails when the blocksize is 2^16, because we're doing this:
	//root1 = fbptr->roots & 0xFFFF;
	//root2 = fbptr->roots >> 16;
	//set_aprime_roots sets roots to all 1's, which then results in root1 and 
	//root2 being set to 65535 in the sieve routine.  this, of course, isn't right
	//so the sieve location 65535 is corrupted by many small prime hits when it
	//shouldn't be, and thus we might end up here more often then we should for
	//offset 65535.
	//
	//even if we do end up here when we shouldn't, often we'll fail to find many
	//small primes which actually divide this location, and we'll bail anyway.  this
	//is safe because we explicitly trial divide by these small primes.  
	//if we make it past the small prime test and go to check the progression
	//of a prime which divides poly_a then the roots we arrive at are false (65535 again)
	//but our computation of the progression will always be 65535 + prime - blocksize,
	//since we set the root to 65535 during the sieve step as well.  
	//65535 != 65535 + prime - blocksize, so we are safe here as well.
	//we may incur more trial division than necessary - is that better than always
	//throwing away block location 65535 - NO, empirically it is much better to just
	//always bail for location 65535 when the blocksize is 65536.

int check_relations_siqs_1(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//not unrolled; for small inputs

	uint32 j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;
	uint64 mask = SCAN_MASK;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j++)
	{
		//check 8 locations simultaneously
		if ((sieveblock[j] & mask) == (uint64)(0))
			continue;

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (k=0;k<8;k++)
		{
			thisloc = (j<<3) + k;
			if ((dconf->sieve[thisloc] & 0x80) == 0)			
				continue;

#ifdef YAFU_64K
			//see discussion near line 377
			if (thisloc ==	65535)
				continue;
#endif
			// log this report
			if (dconf->num_reports < MAX_SIEVE_REPORTS)
				dconf->reports[dconf->num_reports++] = thisloc;			
		}
	}

	if (dconf->num_reports > MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1,blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

int check_relations_siqs_4(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x32; for medium inputs

	uint32 i,j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j+=4)	
	{

#ifdef SIMD_SIEVE_SCAN

		uint32 result = 0;

		SIEVE_SCAN_32;

		if (result == 0)
			continue;

#else
		uint64 mask = SCAN_MASK;

		if (((sieveblock[j] | sieveblock[j+1] | 
			sieveblock[j+2] | sieveblock[j+3]) & 
		      mask) == (uint64)(0))
			continue;
#endif
	
#if defined(SIMD_SIEVE_SCAN)
		// make it safe to perform floating point
		SCAN_CLEAN;
#endif

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<4; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64)(0))
				continue;

			//at least one passed the check, find which one(s) and pass to 
			//trial division stage
			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

#ifdef YAFU_64K
				//see discussion near line 377
				if (thisloc == 65535)
					continue;
#endif
				// log this report
				if (dconf->num_reports < MAX_SIEVE_REPORTS)
					dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

#if defined(SIMD_SIEVE_SCAN)
	// make it safe to perform floating point
	SCAN_CLEAN;
#endif

	if (dconf->num_reports > MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve,dconf->numB-1,blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}

int check_relations_siqs_8(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x64; for large inputs
	uint32 i,j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j+=8)
	{

#ifdef SIMD_SIEVE_SCAN

		uint32 result = 0;

		SIEVE_SCAN_64;

		if (result == 0)
			continue;

#else
		uint64 mask = SCAN_MASK;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3] |
		      sieveblock[j+4] | sieveblock[j+5] | sieveblock[j+6] | sieveblock[j+7]) 
			  & mask) == (uint64)(0))
			continue;
#endif
		
#if defined(SIMD_SIEVE_SCAN)
		// make it safe to perform floating point
		SCAN_CLEAN;
#endif

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<8; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

#ifdef YAFU_64K
				//see discussion near line 377
				if (thisloc == 65535)
					continue;
#endif

				// log this report
				dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

#if defined(SIMD_SIEVE_SCAN)
	// make it safe to perform floating point
	SCAN_CLEAN;
#endif

	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//printf("block %d found %d reports\n", blocknum, dconf->num_reports);

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1, blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}


int check_relations_siqs_16(uint32 blocknum, uint8 parity, 
						   static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//unrolled x128; for large inputs
	uint32 i,j,k,it=BLOCKSIZE>>3;
	uint32 thisloc;
	uint64 *sieveblock;

	sieveblock = (uint64 *)dconf->sieve;
	dconf->num_reports = 0;

	//check for relations
	for (j=0;j<it;j+=16)
	{

#ifdef SIMD_SIEVE_SCAN

		uint32 result = 0;

		SIEVE_SCAN_128;

		if (result == 0)
			continue;

#else
		uint64 mask = SCAN_MASK;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3] |
		      sieveblock[j+4] | sieveblock[j+5] | sieveblock[j+6] | sieveblock[j+7] |
			  sieveblock[j+8] | sieveblock[j+9] | sieveblock[j+10] | sieveblock[j+11] |
		      sieveblock[j+12] | sieveblock[j+13] | sieveblock[j+14] | sieveblock[j+15]
				) & mask) == (uint64)(0))
			continue;
#endif
		

#if defined(SIMD_SIEVE_SCAN)
		// make it safe to perform floating point
		SCAN_CLEAN;
#endif

		//at least one passed the check, find which one(s) and pass to 
		//trial division stage
		for (i=0; i<16; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SCAN_MASK) == (uint64)(0))
				continue;

			for (k=0;k<8;k++)
			{
				thisloc = ((j+i)<<3) + k;
				if ((dconf->sieve[thisloc] & 0x80) == 0)
					continue;

#ifdef YAFU_64K
				//see discussion near line 377
				if (thisloc == 65535)
					continue;
#endif
				// log this report
				dconf->reports[dconf->num_reports++] = thisloc;
			}
		}
	}

#if defined(SIMD_SIEVE_SCAN)
	// make it safe to perform floating point
	SCAN_CLEAN;
#endif

	if (dconf->num_reports >= MAX_SIEVE_REPORTS)
	{
		printf("error: too many sieve reports (found %d)\n",dconf->num_reports);
		exit(-1);
	}

	//remove small primes, and test if its worth continuing for each report
	filter_SPV(parity, dconf->sieve, dconf->numB-1,blocknum,sconf,dconf);
	filter_medprimes(parity, dconf->numB-1,blocknum,sconf,dconf);

	// factor all reports in this block
	for (j=0; j<dconf->num_reports; j++)
	{
		if (dconf->valid_Qs[j])
		{
			filter_LP(j, parity, blocknum, sconf, dconf);
			trial_divide_Q_siqs(j, parity, dconf->numB-1, blocknum,sconf,dconf);
		}
	}

	return 0;
}


