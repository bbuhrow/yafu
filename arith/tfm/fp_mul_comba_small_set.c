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

#define TFM_DEFINES
#include "fp_mul_comba.c"
#include "yafu.h"

#if defined(TFM_SMALL_SET)

void fp_mul_comba_small(z *A, z *B, z *C)
{
	int sz = sizeof(fp_digit);
   fp_digit c0, c1, c2, at[32];

   switch (MAX(abs(A->size), abs(B->size))) { 

   case 1:
      memcpy(at, A->val, 1 * sz);
      memcpy(at+1, B->val, 1 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[1]); 
      COMBA_STORE(C->val[0]);
      COMBA_STORE2(C->val[1]);
      C->size = 2;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 2:
      memcpy(at, A->val, 2 * sz);
      memcpy(at+2, B->val, 2 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[2]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[3]);       MULADD(at[1], at[2]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[1], at[3]); 
      COMBA_STORE(C->val[2]);
      COMBA_STORE2(C->val[3]);
      C->size = 4;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 3:
      memcpy(at, A->val, 3 * sz);
      memcpy(at+3, B->val, 3 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[3]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[4]);       MULADD(at[1], at[3]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[5]);       MULADD(at[1], at[4]);       MULADD(at[2], at[3]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[1], at[5]);       MULADD(at[2], at[4]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[2], at[5]); 
      COMBA_STORE(C->val[4]);
      COMBA_STORE2(C->val[5]);
      C->size = 6;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 4:
      memcpy(at, A->val, 4 * sz);
      memcpy(at+4, B->val, 4 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[4]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[5]);       MULADD(at[1], at[4]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[6]);       MULADD(at[1], at[5]);       MULADD(at[2], at[4]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[7]);       MULADD(at[1], at[6]);       MULADD(at[2], at[5]);       MULADD(at[3], at[4]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[1], at[7]);       MULADD(at[2], at[6]);       MULADD(at[3], at[5]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[2], at[7]);       MULADD(at[3], at[6]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[3], at[7]); 
      COMBA_STORE(C->val[6]);
      COMBA_STORE2(C->val[7]);
      C->size = 8;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 5:
      memcpy(at, A->val, 5 * sz);
      memcpy(at+5, B->val, 5 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[5]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[6]);       MULADD(at[1], at[5]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[7]);       MULADD(at[1], at[6]);       MULADD(at[2], at[5]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[8]);       MULADD(at[1], at[7]);       MULADD(at[2], at[6]);       MULADD(at[3], at[5]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[9]);       MULADD(at[1], at[8]);       MULADD(at[2], at[7]);       MULADD(at[3], at[6]);       MULADD(at[4], at[5]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[1], at[9]);       MULADD(at[2], at[8]);       MULADD(at[3], at[7]);       MULADD(at[4], at[6]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[2], at[9]);       MULADD(at[3], at[8]);       MULADD(at[4], at[7]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[3], at[9]);       MULADD(at[4], at[8]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[4], at[9]); 
      COMBA_STORE(C->val[8]);
      COMBA_STORE2(C->val[9]);
      C->size = 10;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 6:
      memcpy(at, A->val, 6 * sz);
      memcpy(at+6, B->val, 6 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[6]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[7]);       MULADD(at[1], at[6]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[8]);       MULADD(at[1], at[7]);       MULADD(at[2], at[6]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[9]);       MULADD(at[1], at[8]);       MULADD(at[2], at[7]);       MULADD(at[3], at[6]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[10]);       MULADD(at[1], at[9]);       MULADD(at[2], at[8]);       MULADD(at[3], at[7]);       MULADD(at[4], at[6]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[11]);       MULADD(at[1], at[10]);       MULADD(at[2], at[9]);       MULADD(at[3], at[8]);       MULADD(at[4], at[7]);       MULADD(at[5], at[6]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[1], at[11]);       MULADD(at[2], at[10]);       MULADD(at[3], at[9]);       MULADD(at[4], at[8]);       MULADD(at[5], at[7]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[2], at[11]);       MULADD(at[3], at[10]);       MULADD(at[4], at[9]);       MULADD(at[5], at[8]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[3], at[11]);       MULADD(at[4], at[10]);       MULADD(at[5], at[9]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[4], at[11]);       MULADD(at[5], at[10]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[5], at[11]); 
      COMBA_STORE(C->val[10]);
      COMBA_STORE2(C->val[11]);
      C->size = 12;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 7:
      memcpy(at, A->val, 7 * sz);
      memcpy(at+7, B->val, 7 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[7]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[8]);       MULADD(at[1], at[7]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[9]);       MULADD(at[1], at[8]);       MULADD(at[2], at[7]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[10]);       MULADD(at[1], at[9]);       MULADD(at[2], at[8]);       MULADD(at[3], at[7]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[11]);       MULADD(at[1], at[10]);       MULADD(at[2], at[9]);       MULADD(at[3], at[8]);       MULADD(at[4], at[7]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[12]);       MULADD(at[1], at[11]);       MULADD(at[2], at[10]);       MULADD(at[3], at[9]);       MULADD(at[4], at[8]);       MULADD(at[5], at[7]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[13]);       MULADD(at[1], at[12]);       MULADD(at[2], at[11]);       MULADD(at[3], at[10]);       MULADD(at[4], at[9]);       MULADD(at[5], at[8]);       MULADD(at[6], at[7]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[1], at[13]);       MULADD(at[2], at[12]);       MULADD(at[3], at[11]);       MULADD(at[4], at[10]);       MULADD(at[5], at[9]);       MULADD(at[6], at[8]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[2], at[13]);       MULADD(at[3], at[12]);       MULADD(at[4], at[11]);       MULADD(at[5], at[10]);       MULADD(at[6], at[9]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[3], at[13]);       MULADD(at[4], at[12]);       MULADD(at[5], at[11]);       MULADD(at[6], at[10]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[4], at[13]);       MULADD(at[5], at[12]);       MULADD(at[6], at[11]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[5], at[13]);       MULADD(at[6], at[12]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[6], at[13]); 
      COMBA_STORE(C->val[12]);
      COMBA_STORE2(C->val[13]);
      C->size = 14;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 8:
      memcpy(at, A->val, 8 * sz);
      memcpy(at+8, B->val, 8 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[8]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[9]);       MULADD(at[1], at[8]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[10]);       MULADD(at[1], at[9]);       MULADD(at[2], at[8]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[11]);       MULADD(at[1], at[10]);       MULADD(at[2], at[9]);       MULADD(at[3], at[8]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[12]);       MULADD(at[1], at[11]);       MULADD(at[2], at[10]);       MULADD(at[3], at[9]);       MULADD(at[4], at[8]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[13]);       MULADD(at[1], at[12]);       MULADD(at[2], at[11]);       MULADD(at[3], at[10]);       MULADD(at[4], at[9]);       MULADD(at[5], at[8]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[14]);       MULADD(at[1], at[13]);       MULADD(at[2], at[12]);       MULADD(at[3], at[11]);       MULADD(at[4], at[10]);       MULADD(at[5], at[9]);       MULADD(at[6], at[8]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[15]);       MULADD(at[1], at[14]);       MULADD(at[2], at[13]);       MULADD(at[3], at[12]);       MULADD(at[4], at[11]);       MULADD(at[5], at[10]);       MULADD(at[6], at[9]);       MULADD(at[7], at[8]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[1], at[15]);       MULADD(at[2], at[14]);       MULADD(at[3], at[13]);       MULADD(at[4], at[12]);       MULADD(at[5], at[11]);       MULADD(at[6], at[10]);       MULADD(at[7], at[9]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[2], at[15]);       MULADD(at[3], at[14]);       MULADD(at[4], at[13]);       MULADD(at[5], at[12]);       MULADD(at[6], at[11]);       MULADD(at[7], at[10]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[3], at[15]);       MULADD(at[4], at[14]);       MULADD(at[5], at[13]);       MULADD(at[6], at[12]);       MULADD(at[7], at[11]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[4], at[15]);       MULADD(at[5], at[14]);       MULADD(at[6], at[13]);       MULADD(at[7], at[12]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[5], at[15]);       MULADD(at[6], at[14]);       MULADD(at[7], at[13]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[6], at[15]);       MULADD(at[7], at[14]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[7], at[15]); 
      COMBA_STORE(C->val[14]);
      COMBA_STORE2(C->val[15]);
      C->size = 16;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 9:
      memcpy(at, A->val, 9 * sz);
      memcpy(at+9, B->val, 9 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[9]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[10]);       MULADD(at[1], at[9]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[11]);       MULADD(at[1], at[10]);       MULADD(at[2], at[9]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[12]);       MULADD(at[1], at[11]);       MULADD(at[2], at[10]);       MULADD(at[3], at[9]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[13]);       MULADD(at[1], at[12]);       MULADD(at[2], at[11]);       MULADD(at[3], at[10]);       MULADD(at[4], at[9]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[14]);       MULADD(at[1], at[13]);       MULADD(at[2], at[12]);       MULADD(at[3], at[11]);       MULADD(at[4], at[10]);       MULADD(at[5], at[9]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[15]);       MULADD(at[1], at[14]);       MULADD(at[2], at[13]);       MULADD(at[3], at[12]);       MULADD(at[4], at[11]);       MULADD(at[5], at[10]);       MULADD(at[6], at[9]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[16]);       MULADD(at[1], at[15]);       MULADD(at[2], at[14]);       MULADD(at[3], at[13]);       MULADD(at[4], at[12]);       MULADD(at[5], at[11]);       MULADD(at[6], at[10]);       MULADD(at[7], at[9]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]);       MULADD(at[2], at[15]);       MULADD(at[3], at[14]);       MULADD(at[4], at[13]);       MULADD(at[5], at[12]);       MULADD(at[6], at[11]);       MULADD(at[7], at[10]);       MULADD(at[8], at[9]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[1], at[17]);       MULADD(at[2], at[16]);       MULADD(at[3], at[15]);       MULADD(at[4], at[14]);       MULADD(at[5], at[13]);       MULADD(at[6], at[12]);       MULADD(at[7], at[11]);       MULADD(at[8], at[10]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[2], at[17]);       MULADD(at[3], at[16]);       MULADD(at[4], at[15]);       MULADD(at[5], at[14]);       MULADD(at[6], at[13]);       MULADD(at[7], at[12]);       MULADD(at[8], at[11]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[3], at[17]);       MULADD(at[4], at[16]);       MULADD(at[5], at[15]);       MULADD(at[6], at[14]);       MULADD(at[7], at[13]);       MULADD(at[8], at[12]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[4], at[17]);       MULADD(at[5], at[16]);       MULADD(at[6], at[15]);       MULADD(at[7], at[14]);       MULADD(at[8], at[13]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[5], at[17]);       MULADD(at[6], at[16]);       MULADD(at[7], at[15]);       MULADD(at[8], at[14]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[6], at[17]);       MULADD(at[7], at[16]);       MULADD(at[8], at[15]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[7], at[17]);       MULADD(at[8], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[8], at[17]); 
      COMBA_STORE(C->val[16]);
      COMBA_STORE2(C->val[17]);
      C->size = 18;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 10:
      memcpy(at, A->val, 10 * sz);
      memcpy(at+10, B->val, 10 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[10]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[11]);       MULADD(at[1], at[10]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[12]);       MULADD(at[1], at[11]);       MULADD(at[2], at[10]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[13]);       MULADD(at[1], at[12]);       MULADD(at[2], at[11]);       MULADD(at[3], at[10]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[14]);       MULADD(at[1], at[13]);       MULADD(at[2], at[12]);       MULADD(at[3], at[11]);       MULADD(at[4], at[10]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[15]);       MULADD(at[1], at[14]);       MULADD(at[2], at[13]);       MULADD(at[3], at[12]);       MULADD(at[4], at[11]);       MULADD(at[5], at[10]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[16]);       MULADD(at[1], at[15]);       MULADD(at[2], at[14]);       MULADD(at[3], at[13]);       MULADD(at[4], at[12]);       MULADD(at[5], at[11]);       MULADD(at[6], at[10]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]);       MULADD(at[2], at[15]);       MULADD(at[3], at[14]);       MULADD(at[4], at[13]);       MULADD(at[5], at[12]);       MULADD(at[6], at[11]);       MULADD(at[7], at[10]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[18]);       MULADD(at[1], at[17]);       MULADD(at[2], at[16]);       MULADD(at[3], at[15]);       MULADD(at[4], at[14]);       MULADD(at[5], at[13]);       MULADD(at[6], at[12]);       MULADD(at[7], at[11]);       MULADD(at[8], at[10]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[0], at[19]);       MULADD(at[1], at[18]);       MULADD(at[2], at[17]);       MULADD(at[3], at[16]);       MULADD(at[4], at[15]);       MULADD(at[5], at[14]);       MULADD(at[6], at[13]);       MULADD(at[7], at[12]);       MULADD(at[8], at[11]);       MULADD(at[9], at[10]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[1], at[19]);       MULADD(at[2], at[18]);       MULADD(at[3], at[17]);       MULADD(at[4], at[16]);       MULADD(at[5], at[15]);       MULADD(at[6], at[14]);       MULADD(at[7], at[13]);       MULADD(at[8], at[12]);       MULADD(at[9], at[11]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[2], at[19]);       MULADD(at[3], at[18]);       MULADD(at[4], at[17]);       MULADD(at[5], at[16]);       MULADD(at[6], at[15]);       MULADD(at[7], at[14]);       MULADD(at[8], at[13]);       MULADD(at[9], at[12]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[3], at[19]);       MULADD(at[4], at[18]);       MULADD(at[5], at[17]);       MULADD(at[6], at[16]);       MULADD(at[7], at[15]);       MULADD(at[8], at[14]);       MULADD(at[9], at[13]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[4], at[19]);       MULADD(at[5], at[18]);       MULADD(at[6], at[17]);       MULADD(at[7], at[16]);       MULADD(at[8], at[15]);       MULADD(at[9], at[14]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[5], at[19]);       MULADD(at[6], at[18]);       MULADD(at[7], at[17]);       MULADD(at[8], at[16]);       MULADD(at[9], at[15]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[6], at[19]);       MULADD(at[7], at[18]);       MULADD(at[8], at[17]);       MULADD(at[9], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[7], at[19]);       MULADD(at[8], at[18]);       MULADD(at[9], at[17]); 
      COMBA_STORE(C->val[16]);
      /* 17 */
      COMBA_FORWARD;
      MULADD(at[8], at[19]);       MULADD(at[9], at[18]); 
      COMBA_STORE(C->val[17]);
      /* 18 */
      COMBA_FORWARD;
      MULADD(at[9], at[19]); 
      COMBA_STORE(C->val[18]);
      COMBA_STORE2(C->val[19]);
      C->size = 20;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 11:
      memcpy(at, A->val, 11 * sz);
      memcpy(at+11, B->val, 11 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[11]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[12]);       MULADD(at[1], at[11]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[13]);       MULADD(at[1], at[12]);       MULADD(at[2], at[11]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[14]);       MULADD(at[1], at[13]);       MULADD(at[2], at[12]);       MULADD(at[3], at[11]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[15]);       MULADD(at[1], at[14]);       MULADD(at[2], at[13]);       MULADD(at[3], at[12]);       MULADD(at[4], at[11]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[16]);       MULADD(at[1], at[15]);       MULADD(at[2], at[14]);       MULADD(at[3], at[13]);       MULADD(at[4], at[12]);       MULADD(at[5], at[11]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]);       MULADD(at[2], at[15]);       MULADD(at[3], at[14]);       MULADD(at[4], at[13]);       MULADD(at[5], at[12]);       MULADD(at[6], at[11]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[18]);       MULADD(at[1], at[17]);       MULADD(at[2], at[16]);       MULADD(at[3], at[15]);       MULADD(at[4], at[14]);       MULADD(at[5], at[13]);       MULADD(at[6], at[12]);       MULADD(at[7], at[11]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[19]);       MULADD(at[1], at[18]);       MULADD(at[2], at[17]);       MULADD(at[3], at[16]);       MULADD(at[4], at[15]);       MULADD(at[5], at[14]);       MULADD(at[6], at[13]);       MULADD(at[7], at[12]);       MULADD(at[8], at[11]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[0], at[20]);       MULADD(at[1], at[19]);       MULADD(at[2], at[18]);       MULADD(at[3], at[17]);       MULADD(at[4], at[16]);       MULADD(at[5], at[15]);       MULADD(at[6], at[14]);       MULADD(at[7], at[13]);       MULADD(at[8], at[12]);       MULADD(at[9], at[11]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[0], at[21]);       MULADD(at[1], at[20]);       MULADD(at[2], at[19]);       MULADD(at[3], at[18]);       MULADD(at[4], at[17]);       MULADD(at[5], at[16]);       MULADD(at[6], at[15]);       MULADD(at[7], at[14]);       MULADD(at[8], at[13]);       MULADD(at[9], at[12]);       MULADD(at[10], at[11]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[1], at[21]);       MULADD(at[2], at[20]);       MULADD(at[3], at[19]);       MULADD(at[4], at[18]);       MULADD(at[5], at[17]);       MULADD(at[6], at[16]);       MULADD(at[7], at[15]);       MULADD(at[8], at[14]);       MULADD(at[9], at[13]);       MULADD(at[10], at[12]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[2], at[21]);       MULADD(at[3], at[20]);       MULADD(at[4], at[19]);       MULADD(at[5], at[18]);       MULADD(at[6], at[17]);       MULADD(at[7], at[16]);       MULADD(at[8], at[15]);       MULADD(at[9], at[14]);       MULADD(at[10], at[13]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[3], at[21]);       MULADD(at[4], at[20]);       MULADD(at[5], at[19]);       MULADD(at[6], at[18]);       MULADD(at[7], at[17]);       MULADD(at[8], at[16]);       MULADD(at[9], at[15]);       MULADD(at[10], at[14]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[4], at[21]);       MULADD(at[5], at[20]);       MULADD(at[6], at[19]);       MULADD(at[7], at[18]);       MULADD(at[8], at[17]);       MULADD(at[9], at[16]);       MULADD(at[10], at[15]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[5], at[21]);       MULADD(at[6], at[20]);       MULADD(at[7], at[19]);       MULADD(at[8], at[18]);       MULADD(at[9], at[17]);       MULADD(at[10], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[6], at[21]);       MULADD(at[7], at[20]);       MULADD(at[8], at[19]);       MULADD(at[9], at[18]);       MULADD(at[10], at[17]); 
      COMBA_STORE(C->val[16]);
      /* 17 */
      COMBA_FORWARD;
      MULADD(at[7], at[21]);       MULADD(at[8], at[20]);       MULADD(at[9], at[19]);       MULADD(at[10], at[18]); 
      COMBA_STORE(C->val[17]);
      /* 18 */
      COMBA_FORWARD;
      MULADD(at[8], at[21]);       MULADD(at[9], at[20]);       MULADD(at[10], at[19]); 
      COMBA_STORE(C->val[18]);
      /* 19 */
      COMBA_FORWARD;
      MULADD(at[9], at[21]);       MULADD(at[10], at[20]); 
      COMBA_STORE(C->val[19]);
      /* 20 */
      COMBA_FORWARD;
      MULADD(at[10], at[21]); 
      COMBA_STORE(C->val[20]);
      COMBA_STORE2(C->val[21]);
      C->size = 22;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 12:
      memcpy(at, A->val, 12 * sz);
      memcpy(at+12, B->val, 12 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[12]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[13]);       MULADD(at[1], at[12]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[14]);       MULADD(at[1], at[13]);       MULADD(at[2], at[12]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[15]);       MULADD(at[1], at[14]);       MULADD(at[2], at[13]);       MULADD(at[3], at[12]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[16]);       MULADD(at[1], at[15]);       MULADD(at[2], at[14]);       MULADD(at[3], at[13]);       MULADD(at[4], at[12]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]);       MULADD(at[2], at[15]);       MULADD(at[3], at[14]);       MULADD(at[4], at[13]);       MULADD(at[5], at[12]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[18]);       MULADD(at[1], at[17]);       MULADD(at[2], at[16]);       MULADD(at[3], at[15]);       MULADD(at[4], at[14]);       MULADD(at[5], at[13]);       MULADD(at[6], at[12]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[19]);       MULADD(at[1], at[18]);       MULADD(at[2], at[17]);       MULADD(at[3], at[16]);       MULADD(at[4], at[15]);       MULADD(at[5], at[14]);       MULADD(at[6], at[13]);       MULADD(at[7], at[12]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[20]);       MULADD(at[1], at[19]);       MULADD(at[2], at[18]);       MULADD(at[3], at[17]);       MULADD(at[4], at[16]);       MULADD(at[5], at[15]);       MULADD(at[6], at[14]);       MULADD(at[7], at[13]);       MULADD(at[8], at[12]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[0], at[21]);       MULADD(at[1], at[20]);       MULADD(at[2], at[19]);       MULADD(at[3], at[18]);       MULADD(at[4], at[17]);       MULADD(at[5], at[16]);       MULADD(at[6], at[15]);       MULADD(at[7], at[14]);       MULADD(at[8], at[13]);       MULADD(at[9], at[12]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[0], at[22]);       MULADD(at[1], at[21]);       MULADD(at[2], at[20]);       MULADD(at[3], at[19]);       MULADD(at[4], at[18]);       MULADD(at[5], at[17]);       MULADD(at[6], at[16]);       MULADD(at[7], at[15]);       MULADD(at[8], at[14]);       MULADD(at[9], at[13]);       MULADD(at[10], at[12]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[0], at[23]);       MULADD(at[1], at[22]);       MULADD(at[2], at[21]);       MULADD(at[3], at[20]);       MULADD(at[4], at[19]);       MULADD(at[5], at[18]);       MULADD(at[6], at[17]);       MULADD(at[7], at[16]);       MULADD(at[8], at[15]);       MULADD(at[9], at[14]);       MULADD(at[10], at[13]);       MULADD(at[11], at[12]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[1], at[23]);       MULADD(at[2], at[22]);       MULADD(at[3], at[21]);       MULADD(at[4], at[20]);       MULADD(at[5], at[19]);       MULADD(at[6], at[18]);       MULADD(at[7], at[17]);       MULADD(at[8], at[16]);       MULADD(at[9], at[15]);       MULADD(at[10], at[14]);       MULADD(at[11], at[13]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[2], at[23]);       MULADD(at[3], at[22]);       MULADD(at[4], at[21]);       MULADD(at[5], at[20]);       MULADD(at[6], at[19]);       MULADD(at[7], at[18]);       MULADD(at[8], at[17]);       MULADD(at[9], at[16]);       MULADD(at[10], at[15]);       MULADD(at[11], at[14]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[3], at[23]);       MULADD(at[4], at[22]);       MULADD(at[5], at[21]);       MULADD(at[6], at[20]);       MULADD(at[7], at[19]);       MULADD(at[8], at[18]);       MULADD(at[9], at[17]);       MULADD(at[10], at[16]);       MULADD(at[11], at[15]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[4], at[23]);       MULADD(at[5], at[22]);       MULADD(at[6], at[21]);       MULADD(at[7], at[20]);       MULADD(at[8], at[19]);       MULADD(at[9], at[18]);       MULADD(at[10], at[17]);       MULADD(at[11], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[5], at[23]);       MULADD(at[6], at[22]);       MULADD(at[7], at[21]);       MULADD(at[8], at[20]);       MULADD(at[9], at[19]);       MULADD(at[10], at[18]);       MULADD(at[11], at[17]); 
      COMBA_STORE(C->val[16]);
      /* 17 */
      COMBA_FORWARD;
      MULADD(at[6], at[23]);       MULADD(at[7], at[22]);       MULADD(at[8], at[21]);       MULADD(at[9], at[20]);       MULADD(at[10], at[19]);       MULADD(at[11], at[18]); 
      COMBA_STORE(C->val[17]);
      /* 18 */
      COMBA_FORWARD;
      MULADD(at[7], at[23]);       MULADD(at[8], at[22]);       MULADD(at[9], at[21]);       MULADD(at[10], at[20]);       MULADD(at[11], at[19]); 
      COMBA_STORE(C->val[18]);
      /* 19 */
      COMBA_FORWARD;
      MULADD(at[8], at[23]);       MULADD(at[9], at[22]);       MULADD(at[10], at[21]);       MULADD(at[11], at[20]); 
      COMBA_STORE(C->val[19]);
      /* 20 */
      COMBA_FORWARD;
      MULADD(at[9], at[23]);       MULADD(at[10], at[22]);       MULADD(at[11], at[21]); 
      COMBA_STORE(C->val[20]);
      /* 21 */
      COMBA_FORWARD;
      MULADD(at[10], at[23]);       MULADD(at[11], at[22]); 
      COMBA_STORE(C->val[21]);
      /* 22 */
      COMBA_FORWARD;
      MULADD(at[11], at[23]); 
      COMBA_STORE(C->val[22]);
      COMBA_STORE2(C->val[23]);
      C->size = 24;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 13:
      memcpy(at, A->val, 13 * sz);
      memcpy(at+13, B->val, 13 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[13]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[14]);       MULADD(at[1], at[13]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[15]);       MULADD(at[1], at[14]);       MULADD(at[2], at[13]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[16]);       MULADD(at[1], at[15]);       MULADD(at[2], at[14]);       MULADD(at[3], at[13]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]);       MULADD(at[2], at[15]);       MULADD(at[3], at[14]);       MULADD(at[4], at[13]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[18]);       MULADD(at[1], at[17]);       MULADD(at[2], at[16]);       MULADD(at[3], at[15]);       MULADD(at[4], at[14]);       MULADD(at[5], at[13]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[19]);       MULADD(at[1], at[18]);       MULADD(at[2], at[17]);       MULADD(at[3], at[16]);       MULADD(at[4], at[15]);       MULADD(at[5], at[14]);       MULADD(at[6], at[13]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[20]);       MULADD(at[1], at[19]);       MULADD(at[2], at[18]);       MULADD(at[3], at[17]);       MULADD(at[4], at[16]);       MULADD(at[5], at[15]);       MULADD(at[6], at[14]);       MULADD(at[7], at[13]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[21]);       MULADD(at[1], at[20]);       MULADD(at[2], at[19]);       MULADD(at[3], at[18]);       MULADD(at[4], at[17]);       MULADD(at[5], at[16]);       MULADD(at[6], at[15]);       MULADD(at[7], at[14]);       MULADD(at[8], at[13]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[0], at[22]);       MULADD(at[1], at[21]);       MULADD(at[2], at[20]);       MULADD(at[3], at[19]);       MULADD(at[4], at[18]);       MULADD(at[5], at[17]);       MULADD(at[6], at[16]);       MULADD(at[7], at[15]);       MULADD(at[8], at[14]);       MULADD(at[9], at[13]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[0], at[23]);       MULADD(at[1], at[22]);       MULADD(at[2], at[21]);       MULADD(at[3], at[20]);       MULADD(at[4], at[19]);       MULADD(at[5], at[18]);       MULADD(at[6], at[17]);       MULADD(at[7], at[16]);       MULADD(at[8], at[15]);       MULADD(at[9], at[14]);       MULADD(at[10], at[13]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[0], at[24]);       MULADD(at[1], at[23]);       MULADD(at[2], at[22]);       MULADD(at[3], at[21]);       MULADD(at[4], at[20]);       MULADD(at[5], at[19]);       MULADD(at[6], at[18]);       MULADD(at[7], at[17]);       MULADD(at[8], at[16]);       MULADD(at[9], at[15]);       MULADD(at[10], at[14]);       MULADD(at[11], at[13]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[0], at[25]);       MULADD(at[1], at[24]);       MULADD(at[2], at[23]);       MULADD(at[3], at[22]);       MULADD(at[4], at[21]);       MULADD(at[5], at[20]);       MULADD(at[6], at[19]);       MULADD(at[7], at[18]);       MULADD(at[8], at[17]);       MULADD(at[9], at[16]);       MULADD(at[10], at[15]);       MULADD(at[11], at[14]);       MULADD(at[12], at[13]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[1], at[25]);       MULADD(at[2], at[24]);       MULADD(at[3], at[23]);       MULADD(at[4], at[22]);       MULADD(at[5], at[21]);       MULADD(at[6], at[20]);       MULADD(at[7], at[19]);       MULADD(at[8], at[18]);       MULADD(at[9], at[17]);       MULADD(at[10], at[16]);       MULADD(at[11], at[15]);       MULADD(at[12], at[14]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[2], at[25]);       MULADD(at[3], at[24]);       MULADD(at[4], at[23]);       MULADD(at[5], at[22]);       MULADD(at[6], at[21]);       MULADD(at[7], at[20]);       MULADD(at[8], at[19]);       MULADD(at[9], at[18]);       MULADD(at[10], at[17]);       MULADD(at[11], at[16]);       MULADD(at[12], at[15]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[3], at[25]);       MULADD(at[4], at[24]);       MULADD(at[5], at[23]);       MULADD(at[6], at[22]);       MULADD(at[7], at[21]);       MULADD(at[8], at[20]);       MULADD(at[9], at[19]);       MULADD(at[10], at[18]);       MULADD(at[11], at[17]);       MULADD(at[12], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[4], at[25]);       MULADD(at[5], at[24]);       MULADD(at[6], at[23]);       MULADD(at[7], at[22]);       MULADD(at[8], at[21]);       MULADD(at[9], at[20]);       MULADD(at[10], at[19]);       MULADD(at[11], at[18]);       MULADD(at[12], at[17]); 
      COMBA_STORE(C->val[16]);
      /* 17 */
      COMBA_FORWARD;
      MULADD(at[5], at[25]);       MULADD(at[6], at[24]);       MULADD(at[7], at[23]);       MULADD(at[8], at[22]);       MULADD(at[9], at[21]);       MULADD(at[10], at[20]);       MULADD(at[11], at[19]);       MULADD(at[12], at[18]); 
      COMBA_STORE(C->val[17]);
      /* 18 */
      COMBA_FORWARD;
      MULADD(at[6], at[25]);       MULADD(at[7], at[24]);       MULADD(at[8], at[23]);       MULADD(at[9], at[22]);       MULADD(at[10], at[21]);       MULADD(at[11], at[20]);       MULADD(at[12], at[19]); 
      COMBA_STORE(C->val[18]);
      /* 19 */
      COMBA_FORWARD;
      MULADD(at[7], at[25]);       MULADD(at[8], at[24]);       MULADD(at[9], at[23]);       MULADD(at[10], at[22]);       MULADD(at[11], at[21]);       MULADD(at[12], at[20]); 
      COMBA_STORE(C->val[19]);
      /* 20 */
      COMBA_FORWARD;
      MULADD(at[8], at[25]);       MULADD(at[9], at[24]);       MULADD(at[10], at[23]);       MULADD(at[11], at[22]);       MULADD(at[12], at[21]); 
      COMBA_STORE(C->val[20]);
      /* 21 */
      COMBA_FORWARD;
      MULADD(at[9], at[25]);       MULADD(at[10], at[24]);       MULADD(at[11], at[23]);       MULADD(at[12], at[22]); 
      COMBA_STORE(C->val[21]);
      /* 22 */
      COMBA_FORWARD;
      MULADD(at[10], at[25]);       MULADD(at[11], at[24]);       MULADD(at[12], at[23]); 
      COMBA_STORE(C->val[22]);
      /* 23 */
      COMBA_FORWARD;
      MULADD(at[11], at[25]);       MULADD(at[12], at[24]); 
      COMBA_STORE(C->val[23]);
      /* 24 */
      COMBA_FORWARD;
      MULADD(at[12], at[25]); 
      COMBA_STORE(C->val[24]);
      COMBA_STORE2(C->val[25]);
      C->size = 26;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 14:
      memcpy(at, A->val, 14 * sz);
      memcpy(at+14, B->val, 14 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[14]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[15]);       MULADD(at[1], at[14]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[16]);       MULADD(at[1], at[15]);       MULADD(at[2], at[14]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]);       MULADD(at[2], at[15]);       MULADD(at[3], at[14]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[18]);       MULADD(at[1], at[17]);       MULADD(at[2], at[16]);       MULADD(at[3], at[15]);       MULADD(at[4], at[14]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[19]);       MULADD(at[1], at[18]);       MULADD(at[2], at[17]);       MULADD(at[3], at[16]);       MULADD(at[4], at[15]);       MULADD(at[5], at[14]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[20]);       MULADD(at[1], at[19]);       MULADD(at[2], at[18]);       MULADD(at[3], at[17]);       MULADD(at[4], at[16]);       MULADD(at[5], at[15]);       MULADD(at[6], at[14]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[21]);       MULADD(at[1], at[20]);       MULADD(at[2], at[19]);       MULADD(at[3], at[18]);       MULADD(at[4], at[17]);       MULADD(at[5], at[16]);       MULADD(at[6], at[15]);       MULADD(at[7], at[14]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[22]);       MULADD(at[1], at[21]);       MULADD(at[2], at[20]);       MULADD(at[3], at[19]);       MULADD(at[4], at[18]);       MULADD(at[5], at[17]);       MULADD(at[6], at[16]);       MULADD(at[7], at[15]);       MULADD(at[8], at[14]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[0], at[23]);       MULADD(at[1], at[22]);       MULADD(at[2], at[21]);       MULADD(at[3], at[20]);       MULADD(at[4], at[19]);       MULADD(at[5], at[18]);       MULADD(at[6], at[17]);       MULADD(at[7], at[16]);       MULADD(at[8], at[15]);       MULADD(at[9], at[14]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[0], at[24]);       MULADD(at[1], at[23]);       MULADD(at[2], at[22]);       MULADD(at[3], at[21]);       MULADD(at[4], at[20]);       MULADD(at[5], at[19]);       MULADD(at[6], at[18]);       MULADD(at[7], at[17]);       MULADD(at[8], at[16]);       MULADD(at[9], at[15]);       MULADD(at[10], at[14]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[0], at[25]);       MULADD(at[1], at[24]);       MULADD(at[2], at[23]);       MULADD(at[3], at[22]);       MULADD(at[4], at[21]);       MULADD(at[5], at[20]);       MULADD(at[6], at[19]);       MULADD(at[7], at[18]);       MULADD(at[8], at[17]);       MULADD(at[9], at[16]);       MULADD(at[10], at[15]);       MULADD(at[11], at[14]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[0], at[26]);       MULADD(at[1], at[25]);       MULADD(at[2], at[24]);       MULADD(at[3], at[23]);       MULADD(at[4], at[22]);       MULADD(at[5], at[21]);       MULADD(at[6], at[20]);       MULADD(at[7], at[19]);       MULADD(at[8], at[18]);       MULADD(at[9], at[17]);       MULADD(at[10], at[16]);       MULADD(at[11], at[15]);       MULADD(at[12], at[14]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[0], at[27]);       MULADD(at[1], at[26]);       MULADD(at[2], at[25]);       MULADD(at[3], at[24]);       MULADD(at[4], at[23]);       MULADD(at[5], at[22]);       MULADD(at[6], at[21]);       MULADD(at[7], at[20]);       MULADD(at[8], at[19]);       MULADD(at[9], at[18]);       MULADD(at[10], at[17]);       MULADD(at[11], at[16]);       MULADD(at[12], at[15]);       MULADD(at[13], at[14]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[1], at[27]);       MULADD(at[2], at[26]);       MULADD(at[3], at[25]);       MULADD(at[4], at[24]);       MULADD(at[5], at[23]);       MULADD(at[6], at[22]);       MULADD(at[7], at[21]);       MULADD(at[8], at[20]);       MULADD(at[9], at[19]);       MULADD(at[10], at[18]);       MULADD(at[11], at[17]);       MULADD(at[12], at[16]);       MULADD(at[13], at[15]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[2], at[27]);       MULADD(at[3], at[26]);       MULADD(at[4], at[25]);       MULADD(at[5], at[24]);       MULADD(at[6], at[23]);       MULADD(at[7], at[22]);       MULADD(at[8], at[21]);       MULADD(at[9], at[20]);       MULADD(at[10], at[19]);       MULADD(at[11], at[18]);       MULADD(at[12], at[17]);       MULADD(at[13], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[3], at[27]);       MULADD(at[4], at[26]);       MULADD(at[5], at[25]);       MULADD(at[6], at[24]);       MULADD(at[7], at[23]);       MULADD(at[8], at[22]);       MULADD(at[9], at[21]);       MULADD(at[10], at[20]);       MULADD(at[11], at[19]);       MULADD(at[12], at[18]);       MULADD(at[13], at[17]); 
      COMBA_STORE(C->val[16]);
      /* 17 */
      COMBA_FORWARD;
      MULADD(at[4], at[27]);       MULADD(at[5], at[26]);       MULADD(at[6], at[25]);       MULADD(at[7], at[24]);       MULADD(at[8], at[23]);       MULADD(at[9], at[22]);       MULADD(at[10], at[21]);       MULADD(at[11], at[20]);       MULADD(at[12], at[19]);       MULADD(at[13], at[18]); 
      COMBA_STORE(C->val[17]);
      /* 18 */
      COMBA_FORWARD;
      MULADD(at[5], at[27]);       MULADD(at[6], at[26]);       MULADD(at[7], at[25]);       MULADD(at[8], at[24]);       MULADD(at[9], at[23]);       MULADD(at[10], at[22]);       MULADD(at[11], at[21]);       MULADD(at[12], at[20]);       MULADD(at[13], at[19]); 
      COMBA_STORE(C->val[18]);
      /* 19 */
      COMBA_FORWARD;
      MULADD(at[6], at[27]);       MULADD(at[7], at[26]);       MULADD(at[8], at[25]);       MULADD(at[9], at[24]);       MULADD(at[10], at[23]);       MULADD(at[11], at[22]);       MULADD(at[12], at[21]);       MULADD(at[13], at[20]); 
      COMBA_STORE(C->val[19]);
      /* 20 */
      COMBA_FORWARD;
      MULADD(at[7], at[27]);       MULADD(at[8], at[26]);       MULADD(at[9], at[25]);       MULADD(at[10], at[24]);       MULADD(at[11], at[23]);       MULADD(at[12], at[22]);       MULADD(at[13], at[21]); 
      COMBA_STORE(C->val[20]);
      /* 21 */
      COMBA_FORWARD;
      MULADD(at[8], at[27]);       MULADD(at[9], at[26]);       MULADD(at[10], at[25]);       MULADD(at[11], at[24]);       MULADD(at[12], at[23]);       MULADD(at[13], at[22]); 
      COMBA_STORE(C->val[21]);
      /* 22 */
      COMBA_FORWARD;
      MULADD(at[9], at[27]);       MULADD(at[10], at[26]);       MULADD(at[11], at[25]);       MULADD(at[12], at[24]);       MULADD(at[13], at[23]); 
      COMBA_STORE(C->val[22]);
      /* 23 */
      COMBA_FORWARD;
      MULADD(at[10], at[27]);       MULADD(at[11], at[26]);       MULADD(at[12], at[25]);       MULADD(at[13], at[24]); 
      COMBA_STORE(C->val[23]);
      /* 24 */
      COMBA_FORWARD;
      MULADD(at[11], at[27]);       MULADD(at[12], at[26]);       MULADD(at[13], at[25]); 
      COMBA_STORE(C->val[24]);
      /* 25 */
      COMBA_FORWARD;
      MULADD(at[12], at[27]);       MULADD(at[13], at[26]); 
      COMBA_STORE(C->val[25]);
      /* 26 */
      COMBA_FORWARD;
      MULADD(at[13], at[27]); 
      COMBA_STORE(C->val[26]);
      COMBA_STORE2(C->val[27]);
      C->size = 28;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 15:
      memcpy(at, A->val, 15 * sz);
      memcpy(at+15, B->val, 15 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[15]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[16]);       MULADD(at[1], at[15]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]);       MULADD(at[2], at[15]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[18]);       MULADD(at[1], at[17]);       MULADD(at[2], at[16]);       MULADD(at[3], at[15]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[19]);       MULADD(at[1], at[18]);       MULADD(at[2], at[17]);       MULADD(at[3], at[16]);       MULADD(at[4], at[15]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[20]);       MULADD(at[1], at[19]);       MULADD(at[2], at[18]);       MULADD(at[3], at[17]);       MULADD(at[4], at[16]);       MULADD(at[5], at[15]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[21]);       MULADD(at[1], at[20]);       MULADD(at[2], at[19]);       MULADD(at[3], at[18]);       MULADD(at[4], at[17]);       MULADD(at[5], at[16]);       MULADD(at[6], at[15]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[22]);       MULADD(at[1], at[21]);       MULADD(at[2], at[20]);       MULADD(at[3], at[19]);       MULADD(at[4], at[18]);       MULADD(at[5], at[17]);       MULADD(at[6], at[16]);       MULADD(at[7], at[15]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[23]);       MULADD(at[1], at[22]);       MULADD(at[2], at[21]);       MULADD(at[3], at[20]);       MULADD(at[4], at[19]);       MULADD(at[5], at[18]);       MULADD(at[6], at[17]);       MULADD(at[7], at[16]);       MULADD(at[8], at[15]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[0], at[24]);       MULADD(at[1], at[23]);       MULADD(at[2], at[22]);       MULADD(at[3], at[21]);       MULADD(at[4], at[20]);       MULADD(at[5], at[19]);       MULADD(at[6], at[18]);       MULADD(at[7], at[17]);       MULADD(at[8], at[16]);       MULADD(at[9], at[15]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[0], at[25]);       MULADD(at[1], at[24]);       MULADD(at[2], at[23]);       MULADD(at[3], at[22]);       MULADD(at[4], at[21]);       MULADD(at[5], at[20]);       MULADD(at[6], at[19]);       MULADD(at[7], at[18]);       MULADD(at[8], at[17]);       MULADD(at[9], at[16]);       MULADD(at[10], at[15]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[0], at[26]);       MULADD(at[1], at[25]);       MULADD(at[2], at[24]);       MULADD(at[3], at[23]);       MULADD(at[4], at[22]);       MULADD(at[5], at[21]);       MULADD(at[6], at[20]);       MULADD(at[7], at[19]);       MULADD(at[8], at[18]);       MULADD(at[9], at[17]);       MULADD(at[10], at[16]);       MULADD(at[11], at[15]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[0], at[27]);       MULADD(at[1], at[26]);       MULADD(at[2], at[25]);       MULADD(at[3], at[24]);       MULADD(at[4], at[23]);       MULADD(at[5], at[22]);       MULADD(at[6], at[21]);       MULADD(at[7], at[20]);       MULADD(at[8], at[19]);       MULADD(at[9], at[18]);       MULADD(at[10], at[17]);       MULADD(at[11], at[16]);       MULADD(at[12], at[15]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[0], at[28]);       MULADD(at[1], at[27]);       MULADD(at[2], at[26]);       MULADD(at[3], at[25]);       MULADD(at[4], at[24]);       MULADD(at[5], at[23]);       MULADD(at[6], at[22]);       MULADD(at[7], at[21]);       MULADD(at[8], at[20]);       MULADD(at[9], at[19]);       MULADD(at[10], at[18]);       MULADD(at[11], at[17]);       MULADD(at[12], at[16]);       MULADD(at[13], at[15]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[0], at[29]);       MULADD(at[1], at[28]);       MULADD(at[2], at[27]);       MULADD(at[3], at[26]);       MULADD(at[4], at[25]);       MULADD(at[5], at[24]);       MULADD(at[6], at[23]);       MULADD(at[7], at[22]);       MULADD(at[8], at[21]);       MULADD(at[9], at[20]);       MULADD(at[10], at[19]);       MULADD(at[11], at[18]);       MULADD(at[12], at[17]);       MULADD(at[13], at[16]);       MULADD(at[14], at[15]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[1], at[29]);       MULADD(at[2], at[28]);       MULADD(at[3], at[27]);       MULADD(at[4], at[26]);       MULADD(at[5], at[25]);       MULADD(at[6], at[24]);       MULADD(at[7], at[23]);       MULADD(at[8], at[22]);       MULADD(at[9], at[21]);       MULADD(at[10], at[20]);       MULADD(at[11], at[19]);       MULADD(at[12], at[18]);       MULADD(at[13], at[17]);       MULADD(at[14], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[2], at[29]);       MULADD(at[3], at[28]);       MULADD(at[4], at[27]);       MULADD(at[5], at[26]);       MULADD(at[6], at[25]);       MULADD(at[7], at[24]);       MULADD(at[8], at[23]);       MULADD(at[9], at[22]);       MULADD(at[10], at[21]);       MULADD(at[11], at[20]);       MULADD(at[12], at[19]);       MULADD(at[13], at[18]);       MULADD(at[14], at[17]); 
      COMBA_STORE(C->val[16]);
      /* 17 */
      COMBA_FORWARD;
      MULADD(at[3], at[29]);       MULADD(at[4], at[28]);       MULADD(at[5], at[27]);       MULADD(at[6], at[26]);       MULADD(at[7], at[25]);       MULADD(at[8], at[24]);       MULADD(at[9], at[23]);       MULADD(at[10], at[22]);       MULADD(at[11], at[21]);       MULADD(at[12], at[20]);       MULADD(at[13], at[19]);       MULADD(at[14], at[18]); 
      COMBA_STORE(C->val[17]);
      /* 18 */
      COMBA_FORWARD;
      MULADD(at[4], at[29]);       MULADD(at[5], at[28]);       MULADD(at[6], at[27]);       MULADD(at[7], at[26]);       MULADD(at[8], at[25]);       MULADD(at[9], at[24]);       MULADD(at[10], at[23]);       MULADD(at[11], at[22]);       MULADD(at[12], at[21]);       MULADD(at[13], at[20]);       MULADD(at[14], at[19]); 
      COMBA_STORE(C->val[18]);
      /* 19 */
      COMBA_FORWARD;
      MULADD(at[5], at[29]);       MULADD(at[6], at[28]);       MULADD(at[7], at[27]);       MULADD(at[8], at[26]);       MULADD(at[9], at[25]);       MULADD(at[10], at[24]);       MULADD(at[11], at[23]);       MULADD(at[12], at[22]);       MULADD(at[13], at[21]);       MULADD(at[14], at[20]); 
      COMBA_STORE(C->val[19]);
      /* 20 */
      COMBA_FORWARD;
      MULADD(at[6], at[29]);       MULADD(at[7], at[28]);       MULADD(at[8], at[27]);       MULADD(at[9], at[26]);       MULADD(at[10], at[25]);       MULADD(at[11], at[24]);       MULADD(at[12], at[23]);       MULADD(at[13], at[22]);       MULADD(at[14], at[21]); 
      COMBA_STORE(C->val[20]);
      /* 21 */
      COMBA_FORWARD;
      MULADD(at[7], at[29]);       MULADD(at[8], at[28]);       MULADD(at[9], at[27]);       MULADD(at[10], at[26]);       MULADD(at[11], at[25]);       MULADD(at[12], at[24]);       MULADD(at[13], at[23]);       MULADD(at[14], at[22]); 
      COMBA_STORE(C->val[21]);
      /* 22 */
      COMBA_FORWARD;
      MULADD(at[8], at[29]);       MULADD(at[9], at[28]);       MULADD(at[10], at[27]);       MULADD(at[11], at[26]);       MULADD(at[12], at[25]);       MULADD(at[13], at[24]);       MULADD(at[14], at[23]); 
      COMBA_STORE(C->val[22]);
      /* 23 */
      COMBA_FORWARD;
      MULADD(at[9], at[29]);       MULADD(at[10], at[28]);       MULADD(at[11], at[27]);       MULADD(at[12], at[26]);       MULADD(at[13], at[25]);       MULADD(at[14], at[24]); 
      COMBA_STORE(C->val[23]);
      /* 24 */
      COMBA_FORWARD;
      MULADD(at[10], at[29]);       MULADD(at[11], at[28]);       MULADD(at[12], at[27]);       MULADD(at[13], at[26]);       MULADD(at[14], at[25]); 
      COMBA_STORE(C->val[24]);
      /* 25 */
      COMBA_FORWARD;
      MULADD(at[11], at[29]);       MULADD(at[12], at[28]);       MULADD(at[13], at[27]);       MULADD(at[14], at[26]); 
      COMBA_STORE(C->val[25]);
      /* 26 */
      COMBA_FORWARD;
      MULADD(at[12], at[29]);       MULADD(at[13], at[28]);       MULADD(at[14], at[27]); 
      COMBA_STORE(C->val[26]);
      /* 27 */
      COMBA_FORWARD;
      MULADD(at[13], at[29]);       MULADD(at[14], at[28]); 
      COMBA_STORE(C->val[27]);
      /* 28 */
      COMBA_FORWARD;
      MULADD(at[14], at[29]); 
      COMBA_STORE(C->val[28]);
      COMBA_STORE2(C->val[29]);
      C->size = 30;
      fp_clamp(C);
      COMBA_FINI;
      break;

   case 16:
      memcpy(at, A->val, 16 * sz);
      memcpy(at+16, B->val, 16 * sz);
      COMBA_START;

      COMBA_CLEAR;
      /* 0 */
      MULADD(at[0], at[16]); 
      COMBA_STORE(C->val[0]);
      /* 1 */
      COMBA_FORWARD;
      MULADD(at[0], at[17]);       MULADD(at[1], at[16]); 
      COMBA_STORE(C->val[1]);
      /* 2 */
      COMBA_FORWARD;
      MULADD(at[0], at[18]);       MULADD(at[1], at[17]);       MULADD(at[2], at[16]); 
      COMBA_STORE(C->val[2]);
      /* 3 */
      COMBA_FORWARD;
      MULADD(at[0], at[19]);       MULADD(at[1], at[18]);       MULADD(at[2], at[17]);       MULADD(at[3], at[16]); 
      COMBA_STORE(C->val[3]);
      /* 4 */
      COMBA_FORWARD;
      MULADD(at[0], at[20]);       MULADD(at[1], at[19]);       MULADD(at[2], at[18]);       MULADD(at[3], at[17]);       MULADD(at[4], at[16]); 
      COMBA_STORE(C->val[4]);
      /* 5 */
      COMBA_FORWARD;
      MULADD(at[0], at[21]);       MULADD(at[1], at[20]);       MULADD(at[2], at[19]);       MULADD(at[3], at[18]);       MULADD(at[4], at[17]);       MULADD(at[5], at[16]); 
      COMBA_STORE(C->val[5]);
      /* 6 */
      COMBA_FORWARD;
      MULADD(at[0], at[22]);       MULADD(at[1], at[21]);       MULADD(at[2], at[20]);       MULADD(at[3], at[19]);       MULADD(at[4], at[18]);       MULADD(at[5], at[17]);       MULADD(at[6], at[16]); 
      COMBA_STORE(C->val[6]);
      /* 7 */
      COMBA_FORWARD;
      MULADD(at[0], at[23]);       MULADD(at[1], at[22]);       MULADD(at[2], at[21]);       MULADD(at[3], at[20]);       MULADD(at[4], at[19]);       MULADD(at[5], at[18]);       MULADD(at[6], at[17]);       MULADD(at[7], at[16]); 
      COMBA_STORE(C->val[7]);
      /* 8 */
      COMBA_FORWARD;
      MULADD(at[0], at[24]);       MULADD(at[1], at[23]);       MULADD(at[2], at[22]);       MULADD(at[3], at[21]);       MULADD(at[4], at[20]);       MULADD(at[5], at[19]);       MULADD(at[6], at[18]);       MULADD(at[7], at[17]);       MULADD(at[8], at[16]); 
      COMBA_STORE(C->val[8]);
      /* 9 */
      COMBA_FORWARD;
      MULADD(at[0], at[25]);       MULADD(at[1], at[24]);       MULADD(at[2], at[23]);       MULADD(at[3], at[22]);       MULADD(at[4], at[21]);       MULADD(at[5], at[20]);       MULADD(at[6], at[19]);       MULADD(at[7], at[18]);       MULADD(at[8], at[17]);       MULADD(at[9], at[16]); 
      COMBA_STORE(C->val[9]);
      /* 10 */
      COMBA_FORWARD;
      MULADD(at[0], at[26]);       MULADD(at[1], at[25]);       MULADD(at[2], at[24]);       MULADD(at[3], at[23]);       MULADD(at[4], at[22]);       MULADD(at[5], at[21]);       MULADD(at[6], at[20]);       MULADD(at[7], at[19]);       MULADD(at[8], at[18]);       MULADD(at[9], at[17]);       MULADD(at[10], at[16]); 
      COMBA_STORE(C->val[10]);
      /* 11 */
      COMBA_FORWARD;
      MULADD(at[0], at[27]);       MULADD(at[1], at[26]);       MULADD(at[2], at[25]);       MULADD(at[3], at[24]);       MULADD(at[4], at[23]);       MULADD(at[5], at[22]);       MULADD(at[6], at[21]);       MULADD(at[7], at[20]);       MULADD(at[8], at[19]);       MULADD(at[9], at[18]);       MULADD(at[10], at[17]);       MULADD(at[11], at[16]); 
      COMBA_STORE(C->val[11]);
      /* 12 */
      COMBA_FORWARD;
      MULADD(at[0], at[28]);       MULADD(at[1], at[27]);       MULADD(at[2], at[26]);       MULADD(at[3], at[25]);       MULADD(at[4], at[24]);       MULADD(at[5], at[23]);       MULADD(at[6], at[22]);       MULADD(at[7], at[21]);       MULADD(at[8], at[20]);       MULADD(at[9], at[19]);       MULADD(at[10], at[18]);       MULADD(at[11], at[17]);       MULADD(at[12], at[16]); 
      COMBA_STORE(C->val[12]);
      /* 13 */
      COMBA_FORWARD;
      MULADD(at[0], at[29]);       MULADD(at[1], at[28]);       MULADD(at[2], at[27]);       MULADD(at[3], at[26]);       MULADD(at[4], at[25]);       MULADD(at[5], at[24]);       MULADD(at[6], at[23]);       MULADD(at[7], at[22]);       MULADD(at[8], at[21]);       MULADD(at[9], at[20]);       MULADD(at[10], at[19]);       MULADD(at[11], at[18]);       MULADD(at[12], at[17]);       MULADD(at[13], at[16]); 
      COMBA_STORE(C->val[13]);
      /* 14 */
      COMBA_FORWARD;
      MULADD(at[0], at[30]);       MULADD(at[1], at[29]);       MULADD(at[2], at[28]);       MULADD(at[3], at[27]);       MULADD(at[4], at[26]);       MULADD(at[5], at[25]);       MULADD(at[6], at[24]);       MULADD(at[7], at[23]);       MULADD(at[8], at[22]);       MULADD(at[9], at[21]);       MULADD(at[10], at[20]);       MULADD(at[11], at[19]);       MULADD(at[12], at[18]);       MULADD(at[13], at[17]);       MULADD(at[14], at[16]); 
      COMBA_STORE(C->val[14]);
      /* 15 */
      COMBA_FORWARD;
      MULADD(at[0], at[31]);       MULADD(at[1], at[30]);       MULADD(at[2], at[29]);       MULADD(at[3], at[28]);       MULADD(at[4], at[27]);       MULADD(at[5], at[26]);       MULADD(at[6], at[25]);       MULADD(at[7], at[24]);       MULADD(at[8], at[23]);       MULADD(at[9], at[22]);       MULADD(at[10], at[21]);       MULADD(at[11], at[20]);       MULADD(at[12], at[19]);       MULADD(at[13], at[18]);       MULADD(at[14], at[17]);       MULADD(at[15], at[16]); 
      COMBA_STORE(C->val[15]);
      /* 16 */
      COMBA_FORWARD;
      MULADD(at[1], at[31]);       MULADD(at[2], at[30]);       MULADD(at[3], at[29]);       MULADD(at[4], at[28]);       MULADD(at[5], at[27]);       MULADD(at[6], at[26]);       MULADD(at[7], at[25]);       MULADD(at[8], at[24]);       MULADD(at[9], at[23]);       MULADD(at[10], at[22]);       MULADD(at[11], at[21]);       MULADD(at[12], at[20]);       MULADD(at[13], at[19]);       MULADD(at[14], at[18]);       MULADD(at[15], at[17]); 
      COMBA_STORE(C->val[16]);
      /* 17 */
      COMBA_FORWARD;
      MULADD(at[2], at[31]);       MULADD(at[3], at[30]);       MULADD(at[4], at[29]);       MULADD(at[5], at[28]);       MULADD(at[6], at[27]);       MULADD(at[7], at[26]);       MULADD(at[8], at[25]);       MULADD(at[9], at[24]);       MULADD(at[10], at[23]);       MULADD(at[11], at[22]);       MULADD(at[12], at[21]);       MULADD(at[13], at[20]);       MULADD(at[14], at[19]);       MULADD(at[15], at[18]); 
      COMBA_STORE(C->val[17]);
      /* 18 */
      COMBA_FORWARD;
      MULADD(at[3], at[31]);       MULADD(at[4], at[30]);       MULADD(at[5], at[29]);       MULADD(at[6], at[28]);       MULADD(at[7], at[27]);       MULADD(at[8], at[26]);       MULADD(at[9], at[25]);       MULADD(at[10], at[24]);       MULADD(at[11], at[23]);       MULADD(at[12], at[22]);       MULADD(at[13], at[21]);       MULADD(at[14], at[20]);       MULADD(at[15], at[19]); 
      COMBA_STORE(C->val[18]);
      /* 19 */
      COMBA_FORWARD;
      MULADD(at[4], at[31]);       MULADD(at[5], at[30]);       MULADD(at[6], at[29]);       MULADD(at[7], at[28]);       MULADD(at[8], at[27]);       MULADD(at[9], at[26]);       MULADD(at[10], at[25]);       MULADD(at[11], at[24]);       MULADD(at[12], at[23]);       MULADD(at[13], at[22]);       MULADD(at[14], at[21]);       MULADD(at[15], at[20]); 
      COMBA_STORE(C->val[19]);
      /* 20 */
      COMBA_FORWARD;
      MULADD(at[5], at[31]);       MULADD(at[6], at[30]);       MULADD(at[7], at[29]);       MULADD(at[8], at[28]);       MULADD(at[9], at[27]);       MULADD(at[10], at[26]);       MULADD(at[11], at[25]);       MULADD(at[12], at[24]);       MULADD(at[13], at[23]);       MULADD(at[14], at[22]);       MULADD(at[15], at[21]); 
      COMBA_STORE(C->val[20]);
      /* 21 */
      COMBA_FORWARD;
      MULADD(at[6], at[31]);       MULADD(at[7], at[30]);       MULADD(at[8], at[29]);       MULADD(at[9], at[28]);       MULADD(at[10], at[27]);       MULADD(at[11], at[26]);       MULADD(at[12], at[25]);       MULADD(at[13], at[24]);       MULADD(at[14], at[23]);       MULADD(at[15], at[22]); 
      COMBA_STORE(C->val[21]);
      /* 22 */
      COMBA_FORWARD;
      MULADD(at[7], at[31]);       MULADD(at[8], at[30]);       MULADD(at[9], at[29]);       MULADD(at[10], at[28]);       MULADD(at[11], at[27]);       MULADD(at[12], at[26]);       MULADD(at[13], at[25]);       MULADD(at[14], at[24]);       MULADD(at[15], at[23]); 
      COMBA_STORE(C->val[22]);
      /* 23 */
      COMBA_FORWARD;
      MULADD(at[8], at[31]);       MULADD(at[9], at[30]);       MULADD(at[10], at[29]);       MULADD(at[11], at[28]);       MULADD(at[12], at[27]);       MULADD(at[13], at[26]);       MULADD(at[14], at[25]);       MULADD(at[15], at[24]); 
      COMBA_STORE(C->val[23]);
      /* 24 */
      COMBA_FORWARD;
      MULADD(at[9], at[31]);       MULADD(at[10], at[30]);       MULADD(at[11], at[29]);       MULADD(at[12], at[28]);       MULADD(at[13], at[27]);       MULADD(at[14], at[26]);       MULADD(at[15], at[25]); 
      COMBA_STORE(C->val[24]);
      /* 25 */
      COMBA_FORWARD;
      MULADD(at[10], at[31]);       MULADD(at[11], at[30]);       MULADD(at[12], at[29]);       MULADD(at[13], at[28]);       MULADD(at[14], at[27]);       MULADD(at[15], at[26]); 
      COMBA_STORE(C->val[25]);
      /* 26 */
      COMBA_FORWARD;
      MULADD(at[11], at[31]);       MULADD(at[12], at[30]);       MULADD(at[13], at[29]);       MULADD(at[14], at[28]);       MULADD(at[15], at[27]); 
      COMBA_STORE(C->val[26]);
      /* 27 */
      COMBA_FORWARD;
      MULADD(at[12], at[31]);       MULADD(at[13], at[30]);       MULADD(at[14], at[29]);       MULADD(at[15], at[28]); 
      COMBA_STORE(C->val[27]);
      /* 28 */
      COMBA_FORWARD;
      MULADD(at[13], at[31]);       MULADD(at[14], at[30]);       MULADD(at[15], at[29]); 
      COMBA_STORE(C->val[28]);
      /* 29 */
      COMBA_FORWARD;
      MULADD(at[14], at[31]);       MULADD(at[15], at[30]); 
      COMBA_STORE(C->val[29]);
      /* 30 */
      COMBA_FORWARD;
      MULADD(at[15], at[31]); 
      COMBA_STORE(C->val[30]);
      COMBA_STORE2(C->val[31]);
      C->size = 32;
      fp_clamp(C);
      COMBA_FINI;
      break;
	default:
		fp_mul_comba(A,B,C);

   }
}

#endif
