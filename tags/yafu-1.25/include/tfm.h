/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@gmail.com

 * Modified:	Ben Buhrow
 * Date:		11/24/09
 * Purpose:		Port into Yafu-1.14.
 */

#ifndef TFM_H_
#define TFM_H_

#include "yafu.h"
#include "arith.h"
#include "common.h"

/* externally define this symbol to ignore the default settings, useful for changing the build from the make process */
#ifndef TFM_ALREADY_SET

/* do we want the large set of small multiplications ? 
   Enable these if you are going to be doing a lot of small (<= 16 digit) multiplications say in ECC
   Or if you're on a 64-bit machine doing RSA as a 1024-bit integer == 16 digits ;-)
 */
#define TFM_SMALL_SET
#define TFM_SMALL_MONT_SET

/* do we want huge code 
   Enable these if you are doing 20, 24, 28, 32, 48, 64 digit multiplications (useful for RSA)
   Less important on 64-bit machines as 32 digits == 2048 bits
 */
#if 0
#define TFM_MUL3
#define TFM_MUL4
#define TFM_MUL6
#define TFM_MUL7
#define TFM_MUL8
#define TFM_MUL9
#define TFM_MUL12
#define TFM_MUL17
#endif

#define TFM_MUL20
#define TFM_MUL24
#define TFM_MUL28
#define TFM_MUL32
#define TFM_MUL48
#define TFM_MUL64

#if 0
#define TFM_SQR3
#define TFM_SQR4
#define TFM_SQR6
#define TFM_SQR7
#define TFM_SQR8
#define TFM_SQR9
#define TFM_SQR12
#define TFM_SQR17
#endif

#define TFM_SQR20
#define TFM_SQR24
#define TFM_SQR28
#define TFM_SQR32
#define TFM_SQR48
#define TFM_SQR64

#define CHAR_BIT 8

/* do we want some overflow checks
   Not required if you make sure your numbers are within range (e.g. by default a modulus for fp_exptmod() can only be upto 2048 bits long)
 */
/* #define TFM_CHECK */

/* Is the target a P4 Prescott
 */
/* #define TFM_PRESCOTT */

/* Do we want timing resistant fp_exptmod() ?
 * This makes it slower but also timing invariant with respect to the exponent 
 */
/* #define TFM_TIMING_RESISTANT */

#endif /* TFM_ALREADY_SET */



/* Max size of any number in bits.  Basically the largest size you will be multiplying
 * should be half [or smaller] of FP_MAX_SIZE-four_digit
 *
 * You can externally define this or it defaults to 4096-bits [allowing multiplications upto 2048x2048 bits ]
 */
#ifndef FP_MAX_SIZE
   #define FP_MAX_SIZE           (4096+(8*DIGIT_BIT))
#endif

/* will this lib work? */
#if (CHAR_BIT & 7)
   #error CHAR_BIT must be a multiple of eight.
#endif
#if FP_MAX_SIZE % CHAR_BIT
   #error FP_MAX_SIZE must be a multiple of CHAR_BIT
#endif

/* autodetect x86-64 and make sure we are using 64-bit digits with x86-64 asm 
#if defined(__x86_64__)
   #if defined(TFM_X86) || defined(TFM_SSE2) || defined(TFM_ARM) 
       #error x86-64 detected, x86-32/SSE2/ARM optimizations are not valid!
   #endif
   #if !defined(TFM_X86_64) && !defined(TFM_NO_ASM)
      #define TFM_X86	//_64
   #endif
#endif
#if defined(TFM_X86_64)
    #if !defined(FP_64BIT)
       #define FP_64BIT
    #endif
#endif
*/
/* try to detect x86-32 */
/*
#if defined(__i386__) && !defined(TFM_SSE2)
   #if defined(TFM_X86_64) || defined(TFM_ARM) 
       #error x86-32 detected, x86-64/ARM optimizations are not valid!
   #endif
   #if !defined(TFM_X86) && !defined(TFM_NO_ASM)
      #define TFM_X86
   #endif
#endif
   */

/* make sure we're 32-bit for x86-32/sse/arm/ppc32 */
/*
#if (defined(TFM_X86) || defined(TFM_SSE2) || defined(TFM_ARM) || defined(TFM_PPC32)) && defined(FP_64BIT)
   #warning x86-32, SSE2 and ARM, PPC32 optimizations require 32-bit digits (undefining)
   #undef FP_64BIT
#endif
   */

/* multi asms? */
/*
#ifdef TFM_X86
   #define TFM_ASM
#endif
#ifdef TFM_X86_64
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
#ifdef TFM_SSE2
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
#ifdef TFM_ARM
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
#ifdef TFM_PPC32
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
#ifdef TFM_PPC64
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
#ifdef TFM_AVR32
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
   */

/* we want no asm? 
#ifdef TFM_NO_ASM
   #undef TFM_X86
   #undef TFM_X86_64
   #undef TFM_SSE2
   #undef TFM_ARM
   #undef TFM_PPC32
   #undef TFM_PPC64
   #undef TFM_AVR32
   #undef TFM_ASM   
#endif
*/
/* ECC helpers */
/*
#ifdef TFM_ECC192
   #ifdef FP_64BIT
       #define TFM_MUL3
       #define TFM_SQR3
   #else
       #define TFM_MUL6
       #define TFM_SQR6
   #endif
#endif

#ifdef TFM_ECC224
   #ifdef FP_64BIT
       #define TFM_MUL4
       #define TFM_SQR4
   #else
       #define TFM_MUL7
       #define TFM_SQR7
   #endif
#endif

#ifdef TFM_ECC256
   #ifdef FP_64BIT
       #define TFM_MUL4
       #define TFM_SQR4
   #else
       #define TFM_MUL8
       #define TFM_SQR8
   #endif
#endif

#ifdef TFM_ECC384
   #ifdef FP_64BIT
       #define TFM_MUL6
       #define TFM_SQR6
   #else
       #define TFM_MUL12
       #define TFM_SQR12
   #endif
#endif

#ifdef TFM_ECC521
   #ifdef FP_64BIT
       #define TFM_MUL9
       #define TFM_SQR9
   #else
       #define TFM_MUL17
       #define TFM_SQR17
   #endif
#endif
   */


/* some default configurations.
 */


/*
#if defined(FP_64BIT)
   // for GCC only on supported platforms 
#ifndef CRYPT
   typedef unsigned long ulong64;
#endif
   typedef ulong64            fp_digit;
   typedef unsigned long      fp_word __attribute__ ((mode(TI)));
#else
    //this is to make porting into LibTomCrypt easier :-) 
#ifndef CRYPT
   #if defined(_MSC_VER) || defined(__BORLANDC__) 
      typedef unsigned __int64   ulong64;
      typedef signed __int64     long64;
   #else
      typedef unsigned long long ulong64;
      typedef signed long long   long64;
   #endif
#endif
   typedef unsigned long      fp_digit;
   typedef ulong64            fp_word;
#endif
   */
/*
#ifdef _MSC_VER

	typedef unsigned __int32 fp_digit;
	typedef unsigned __int64 fp_word;

#else

	typedef uint32_t fp_digit;
	typedef uint64_t fp_word;

#endif
*/
/* # of digits this is */
#define DIGIT_BIT  (int)((CHAR_BIT) * sizeof(fp_digit))
#define FP_MASK    (fp_digit)(-1)
#define FP_SIZE    (FP_MAX_SIZE/DIGIT_BIT)

/* signs */
#define FP_ZPOS     0
#define FP_NEG      1

/* return codes */
#define FP_OKAY     0
#define FP_VAL      1
#define FP_MEM      2

/* equalities */
#define FP_LT        -1   /* less than */
#define FP_EQ         0   /* equal to */
#define FP_GT         1   /* greater than */

/* replies */
#define FP_YES        1   /* yes response */
#define FP_NO         0   /* no response */


/********************* Tom's fast math stuff **********************/
void fp_mul_comba(z *A, z *B, z *C);
void fp_mul_comba_small(z *A, z *B, z *C);
void fp_sqr_comba(z *A, z *B);
void fp_sqr_comba_small(z *A, z *B);
void fp_2expt(z *a, int b);
int fp_cmp_mag(z *a, z *b);
void fp_montgomery_calc_normalization(z *a, z *b);
void fp_mul_2(z * a, z * b);
int fp_montgomery_setup(z *a, fp_digit *rho);
void fp_montgomery_reduce(z *a, z *m, fp_digit mp);
void s_fp_sub(z *a, z *b, z *c);

#endif


/* $Source: /cvs/libtom/tomsfastmath/src/headers/tfm.h,v $ */
/* $Revision: 1.3 $ */
/* $Date: 2007/02/27 02:38:44 $ */
