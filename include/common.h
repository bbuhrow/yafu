/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/24/08

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.  Much of the functionality in here
			is ignored, however the assembler defines perfected by
			Brian Gladman are particularly important, and several
			of the function prototypes are necessary for the rest
			of the code (savefile and block lanczos stuff).
--------------------------------------------------------------------*/

#ifndef _COMMON_H_
#define _COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#define _(x) #x
#define STRING(x) _(x)

	/* this byzantine complexity sets up the correct assembly
   language syntax based on the compiler, OS and word size 
   
   Where an inline assembler segment is provided in both 
   GCC and MSC format (i.e. alternative sections), the 
   Intel compiler is configured using guards with an A
   suffix to prefer the native version (GCC on Linux/Unix, 
   MSC on Windows). Where an inline assembler segment is 
   only provided in GCC or MSC format but not both (i.e. 
   exclusive sections) the guards have an X suffix.

   The Intel compiler on Windows appears to have some
   bugs in its processing of GCC inline assembler code.
   These are escaped with the _ICL_WIN_ define

   to print compiler predefined macros try:
   gcc -dM -E dummy.h

   */

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) || defined(__clang__)

	#define ASM_G __asm__
	#define ASM_M __asm

	/* for inline assembler on Unix/Linux */
	#if defined(__unix__)
		#if defined(__x86_64__)
			#define GCC_ASM64A
			#define GCC_ASM64X
			#define MSC_ASM64X
		#elif defined(__i386__)
			#define GCC_ASM32A
			#define GCC_ASM32X
			#define MSC_ASM32X
		#endif
	#endif

	/* for inline assembler on Windows */
	#if defined(_WIN32)
		#define _ICL_WIN_
		#if defined(_M_X64)
			#define MSC_ASM64A
			#define MSC_ASM64X
			#define GCC_ASM64X
		#elif defined(_M_IX86)
			#define MSC_ASM32A
			#define MSC_ASM32X
			#define GCC_ASM32X
		#endif
	#endif

#elif defined(__MINGW32__)

	#define ASM_M __asm__
	#define ASM_G __asm__

	#if defined(__x86_64__) 
		#define GCC_ASM64A
		#define GCC_ASM64X
	#elif defined(__i386__)
		#define GCC_ASM32A
		#define GCC_ASM32X
	#endif

#elif defined(__GNUC__)

	#define ASM_G __asm__

	#if defined(__x86_64__) 
		#define GCC_ASM64A
		#define GCC_ASM64X
	#elif defined(__i386__)
		#define GCC_ASM32A
		#define GCC_ASM32X
	#endif

#elif defined(_MSC_VER)

	#define ASM_M __asm

	#if defined(_M_IX86) && !defined(_WIN64)
		#define MSC_ASM32A
		#define MSC_ASM32X
	#elif defined(_WIN64)	
		

	#endif
#endif



/* loop alignment directives need to know whether
   we're using MSVC */

#ifndef _MSC_VER
	#define ALIGN_LOOP   ".p2align 4,,7 \n\t" 
#else
	#define ALIGN_LOOP /* nothing */
#endif


#ifdef __cplusplus
}
#endif

#endif /* _COMMON_H_ */

