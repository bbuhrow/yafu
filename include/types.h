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

#ifndef my_types
#define my_types

/* system-specific stuff ---------------------------------------*/

#ifdef WIN32

	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
	#include <process.h>
	#include <intrin.h>	

#else /* !WIN32 */

	#include <sys/types.h>
	#include <sys/stat.h>
	#include <fcntl.h>
	#include <unistd.h>
	#include <errno.h>
	#include <pthread.h>

#endif /* WIN32 */

/* system-independent header files ------------------------------------*/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/stat.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <signal.h>
#include <memory.h>
#include <malloc.h>


#ifndef _MSC_VER
#include <stdint.h>
#endif

/* basic types  -------------------------------------------------------*/

#if defined(_MSC_VER) || defined(NO_ASM) || defined(__MINGW32__)

	#if defined(__MINGW32__)
		//mingw32 builds using gcc need these types
		//check to see if mingw32 defines GNUC and __i386__
		//and we can combine this with the GNUC cases below

		typedef unsigned char uint8;
		typedef unsigned short uint16;
		typedef uint32_t uint32;
		typedef uint64_t uint64;
		typedef uint32_t fp_digit;
		typedef uint64_t fp_word;
		typedef int32_t fp_signdigit;
		typedef int64_t fp_signword;
		#define MAX_DIGIT 0xffffffff
		#define BITS_PER_DIGIT 32
		#define DEC_DIGIT_PER_WORD 9
		#define HEX_DIGIT_PER_WORD 8
		#define HIBITMASK 0x80000000
		#define MAX_HALF_DIGIT 0xffff
		#define MAX_DEC_WORD 0x3b9aca00
		#define ADDRESS_BITS 2

		#ifndef RS6K
		typedef char int8;
		typedef short int16;
		typedef int32_t int32;
		typedef int64_t int64;
		#endif

		/* portable 64-bit formatting */
		#define PRId64 "I64d"
		#define PRIu64 "I64u"
		#define PRIx64 "I64x"

		//sleep in milliseconds
		#define MySleep(x) usleep((x)*1000)

		

	#elif defined(_WIN64)
		//these types are for MSVC builds on a 64 bit compiler
		//until MSVC supports 64 bit inline assembly, the base 
		//type will continue to be 32 bit, and no assembly will be used.
		//if this changes in the future, then other parts of the code
		//will need to change, for instance the defines in lanczos.c

		#pragma intrinsic(_umul128)
		typedef __int8 int8;
		typedef __int16 int16;
		typedef __int32 int32;
		typedef __int64 int64;
		typedef unsigned __int8 uint8;
		typedef unsigned __int16 uint16;
		typedef unsigned __int32 uint32;
		typedef unsigned __int64 uint64;

		/*
		typedef unsigned __int64 fp_digit;
		typedef unsigned __int64 fp_word;
		typedef __int64 fp_signdigit;
		typedef __int64 fp_signword;
		#define MAX_DIGIT 0xffffffffffffffff
		#define BITS_PER_DIGIT 64
		#define DEC_DIGIT_PER_WORD 20
		#define HEX_DIGIT_PER_WORD 16
		#define HIBITMASK 0x8000000000000000
		#define MAX_HALF_DIGIT 0xffffffff
		#define MAX_DEC_WORD 0x8AC7230489E80000
		#define ADDRESS_BITS 3
		*/

		
		typedef unsigned __int32 fp_digit;
		typedef unsigned __int64 fp_word;
		typedef __int32 fp_signdigit;
		typedef __int64 fp_signword;
		#define MAX_DIGIT 0xffffffff
		#define BITS_PER_DIGIT 32
		#define DEC_DIGIT_PER_WORD 9
		#define HEX_DIGIT_PER_WORD 8
		#define HIBITMASK 0x80000000
		#define MAX_HALF_DIGIT 0xffff
		#define MAX_DEC_WORD 0x3b9aca00
		#define ADDRESS_BITS 2
		

		
		/* portable 64-bit formatting */
		#define PRId64 "I64d"
		#define PRIu64 "I64u"
		#define PRIx64 "I64x"

		//for gettimeofday
		//http://www.openasthra.com/c-tidbits/gettimeofday-function-for-windows/
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64

		struct timeval
		{
			long tv_sec;
			long tv_usec;
		};

		struct timezone 
		{
		  int  tz_minuteswest; /* minutes W of Greenwich */
		  int  tz_dsttime;     /* type of dst correction */
		};

		//sleep in milliseconds
		#define MySleep(x) Sleep((x))

	#else

		//and finally, MSVC builds using a 32 bit compiler
		typedef __int8 int8;
		typedef __int16 int16;
		typedef __int32 int32;
		typedef __int64 int64;
		typedef unsigned __int8 uint8;
		typedef unsigned __int16 uint16;
		typedef unsigned __int32 uint32;
		typedef unsigned __int64 uint64;
		typedef unsigned __int32 fp_digit;
		typedef unsigned __int64 fp_word;
		typedef __int32 fp_signdigit;
		typedef __int64 fp_signword;
		#define MAX_DIGIT 0xffffffff
		#define MAX_HALF_DIGIT 0xffff
		#define BITS_PER_DIGIT 32
		#define DEC_DIGIT_PER_WORD 9
		#define HEX_DIGIT_PER_WORD 8
		#define HIBITMASK 0x80000000
		#define MAX_DEC_WORD 0x3b9aca00
		#define ADDRESS_BITS 2
		
		/* portable 64-bit formatting */
		#define PRId64 "I64d"
		#define PRIu64 "I64u"
		#define PRIx64 "I64x"
		
		//for gettimeofday
		//http://www.openasthra.com/c-tidbits/gettimeofday-function-for-windows/
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64

		struct timeval
		{
			long tv_sec;
			long tv_usec;
		};

		struct timezone 
		{
		  int  tz_minuteswest; /* minutes W of Greenwich */
		  int  tz_dsttime;     /* type of dst correction */
		};

		//sleep in milliseconds
		#define MySleep(x) Sleep((x))

	#endif

#elif defined(__GNUC__) || defined(__INTEL_COMPILER)
	
	//for gettimeofday using gcc
	#include <sys/time.h>

	//sleep in milliseconds
	#define MySleep(x) usleep((x)*1000)

	#if defined(__x86_64__)

		typedef unsigned char uint8;
		typedef unsigned short uint16;
		typedef uint32_t uint32;
		typedef uint64_t uint64;
		typedef uint64_t fp_digit;
		typedef uint64_t fp_word;
		typedef int64_t fp_signdigit;
		typedef int64_t fp_signword;
		#define MAX_DIGIT 0xffffffffffffffff
		#define BITS_PER_DIGIT 64
		#define DEC_DIGIT_PER_WORD 20
		#define HEX_DIGIT_PER_WORD 16
		#define HIBITMASK 0x8000000000000000
		#define MAX_HALF_DIGIT 0xffffffff
		#define MAX_DEC_WORD 0x8AC7230489E80000
		#define ADDRESS_BITS 3

		#define PRId64 "ld"
		#define PRIu64 "lu"
		#define PRIx64 "lx"
		
		#ifndef RS6K
		typedef char int8;
		typedef short int16;
		typedef int32_t int32;
		typedef int64_t int64;
		#endif

	#elif defined(__i386__)

		typedef unsigned char uint8;
		typedef unsigned short uint16;
		typedef uint32_t uint32;
		typedef uint64_t uint64;
		typedef uint32_t fp_digit;
		typedef uint64_t fp_word;
		typedef int32_t fp_signdigit;
		typedef int64_t fp_signword;
		#define MAX_DIGIT 0xffffffff
		#define MAX_HALF_DIGIT 0xffff
		#define BITS_PER_DIGIT 32
		#define DEC_DIGIT_PER_WORD 9
		#define HEX_DIGIT_PER_WORD 8
		#define HIBITMASK 0x80000000
		#define MAX_DEC_WORD 0x3b9aca00
		#define ADDRESS_BITS 2

		#define PRId64 "lld"
		#define PRIu64 "llu"
		#define PRIx64 "llx"
		
		#ifndef RS6K
		typedef char int8;
		typedef short int16;
		typedef int32_t int32;
		typedef int64_t int64;
		#endif

	#else
	
		#error "unrecognized architecture in gcc"

	#endif

#else

	#error "unrecognized compiler"

#endif

#endif /* my types */


