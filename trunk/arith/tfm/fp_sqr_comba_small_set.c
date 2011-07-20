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
#include "fp_sqr_comba.c"
#include "yafu.h"

#if defined(TFM_SMALL_SET)

void fp_sqr_comba_small(z *A, z *B)
{
	int sz = sizeof(fp_digit);
   fp_digit *a, b[32], c0, c1, c2, sc0=0, sc1=0, sc2=0;

#ifdef TFM_ISO
   uint64   tt;   
#endif   
   switch (abs(A->size)) { 
   case 1:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);
      COMBA_STORE2(b[1]);
      COMBA_FINI;

      B->size = 2;
      memcpy(B->val, b, 2 * sz);
      fp_clamp(B);
      break;

   case 2:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);
      COMBA_STORE2(b[3]);
      COMBA_FINI;

      B->size = 4;
      memcpy(B->val, b, 4 * sz);
      fp_clamp(B);
      break;

   case 3:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);
      COMBA_STORE2(b[5]);
      COMBA_FINI;

      B->size = 6;
      memcpy(B->val, b, 6 * sz);
      fp_clamp(B);
      break;

   case 4:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
      SQRADD2(a[2], a[3]); 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
      SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);
      COMBA_STORE2(b[7]);
      COMBA_FINI;

      B->size = 8;
      memcpy(B->val, b, 8 * sz);
      fp_clamp(B);
      break;

   case 5:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
      SQRADD2(a[1], a[4]);    SQRADD2(a[2], a[3]); 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
      SQRADD2(a[2], a[4]);    SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
      SQRADD2(a[3], a[4]); 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
      SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);
      COMBA_STORE2(b[9]);
      COMBA_FINI;

      B->size = 10;
      memcpy(B->val, b, 10 * sz);
      fp_clamp(B);
      break;

   case 6:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
      SQRADD2(a[1], a[5]);    SQRADD2(a[2], a[4]);    SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
      SQRADD2(a[2], a[5]);    SQRADD2(a[3], a[4]); 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
      SQRADD2(a[3], a[5]);    SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
      SQRADD2(a[4], a[5]); 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
      SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);
      COMBA_STORE2(b[11]);
      COMBA_FINI;

      B->size = 12;
      memcpy(B->val, b, 12 * sz);
      fp_clamp(B);
      break;

   case 7:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
      SQRADD2(a[2], a[6]);    SQRADD2(a[3], a[5]);    SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
      SQRADD2(a[3], a[6]);    SQRADD2(a[4], a[5]); 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
      SQRADD2(a[4], a[6]);    SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
      SQRADD2(a[5], a[6]); 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
      SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);
      COMBA_STORE2(b[13]);
      COMBA_FINI;

      B->size = 14;
      memcpy(B->val, b, 14 * sz);
      fp_clamp(B);
      break;

   case 8:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
      SQRADD2(a[3], a[7]);    SQRADD2(a[4], a[6]);    SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
      SQRADD2(a[4], a[7]);    SQRADD2(a[5], a[6]); 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
      SQRADD2(a[5], a[7]);    SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
      SQRADD2(a[6], a[7]); 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
      SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);
      COMBA_STORE2(b[15]);
      COMBA_FINI;

      B->size = 16;
      memcpy(B->val, b, 16 * sz);
      fp_clamp(B);
      break;

   case 9:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
      SQRADD2(a[4], a[8]);    SQRADD2(a[5], a[7]);    SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
      SQRADD2(a[5], a[8]);    SQRADD2(a[6], a[7]); 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
      SQRADD2(a[6], a[8]);    SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
      SQRADD2(a[7], a[8]); 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
      SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);
      COMBA_STORE2(b[17]);
      COMBA_FINI;

      B->size = 18;
      memcpy(B->val, b, 18 * sz);
      fp_clamp(B);
      break;

   case 10:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[9]); SQRADDAC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[9]); SQRADDAC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[9]); SQRADDAC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[9]); SQRADDAC(a[4], a[8]); SQRADDAC(a[5], a[7]); SQRADDDB; SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
   SQRADDSC(a[4], a[9]); SQRADDAC(a[5], a[8]); SQRADDAC(a[6], a[7]); SQRADDDB; 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
      SQRADD2(a[5], a[9]);    SQRADD2(a[6], a[8]);    SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
      SQRADD2(a[6], a[9]);    SQRADD2(a[7], a[8]); 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
      SQRADD2(a[7], a[9]);    SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);

      /* output 17 */
      CARRY_FORWARD;
      SQRADD2(a[8], a[9]); 
      COMBA_STORE(b[17]);

      /* output 18 */
      CARRY_FORWARD;
      SQRADD(a[9], a[9]); 
      COMBA_STORE(b[18]);
      COMBA_STORE2(b[19]);
      COMBA_FINI;

      B->size = 20;
      memcpy(B->val, b, 20 * sz);
      fp_clamp(B);
      break;

   case 11:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[9]); SQRADDAC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[10]); SQRADDAC(a[1], a[9]); SQRADDAC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[10]); SQRADDAC(a[2], a[9]); SQRADDAC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[10]); SQRADDAC(a[3], a[9]); SQRADDAC(a[4], a[8]); SQRADDAC(a[5], a[7]); SQRADDDB; SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[10]); SQRADDAC(a[4], a[9]); SQRADDAC(a[5], a[8]); SQRADDAC(a[6], a[7]); SQRADDDB; 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
   SQRADDSC(a[4], a[10]); SQRADDAC(a[5], a[9]); SQRADDAC(a[6], a[8]); SQRADDDB; SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
   SQRADDSC(a[5], a[10]); SQRADDAC(a[6], a[9]); SQRADDAC(a[7], a[8]); SQRADDDB; 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
      SQRADD2(a[6], a[10]);    SQRADD2(a[7], a[9]);    SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);

      /* output 17 */
      CARRY_FORWARD;
      SQRADD2(a[7], a[10]);    SQRADD2(a[8], a[9]); 
      COMBA_STORE(b[17]);

      /* output 18 */
      CARRY_FORWARD;
      SQRADD2(a[8], a[10]);    SQRADD(a[9], a[9]); 
      COMBA_STORE(b[18]);

      /* output 19 */
      CARRY_FORWARD;
      SQRADD2(a[9], a[10]); 
      COMBA_STORE(b[19]);

      /* output 20 */
      CARRY_FORWARD;
      SQRADD(a[10], a[10]); 
      COMBA_STORE(b[20]);
      COMBA_STORE2(b[21]);
      COMBA_FINI;

      B->size = 22;
      memcpy(B->val, b, 22 * sz);
      fp_clamp(B);
      break;

   case 12:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[9]); SQRADDAC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[10]); SQRADDAC(a[1], a[9]); SQRADDAC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[11]); SQRADDAC(a[1], a[10]); SQRADDAC(a[2], a[9]); SQRADDAC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[11]); SQRADDAC(a[2], a[10]); SQRADDAC(a[3], a[9]); SQRADDAC(a[4], a[8]); SQRADDAC(a[5], a[7]); SQRADDDB; SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[11]); SQRADDAC(a[3], a[10]); SQRADDAC(a[4], a[9]); SQRADDAC(a[5], a[8]); SQRADDAC(a[6], a[7]); SQRADDDB; 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[11]); SQRADDAC(a[4], a[10]); SQRADDAC(a[5], a[9]); SQRADDAC(a[6], a[8]); SQRADDDB; SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
   SQRADDSC(a[4], a[11]); SQRADDAC(a[5], a[10]); SQRADDAC(a[6], a[9]); SQRADDAC(a[7], a[8]); SQRADDDB; 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
   SQRADDSC(a[5], a[11]); SQRADDAC(a[6], a[10]); SQRADDAC(a[7], a[9]); SQRADDDB; SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);

      /* output 17 */
      CARRY_FORWARD;
   SQRADDSC(a[6], a[11]); SQRADDAC(a[7], a[10]); SQRADDAC(a[8], a[9]); SQRADDDB; 
      COMBA_STORE(b[17]);

      /* output 18 */
      CARRY_FORWARD;
      SQRADD2(a[7], a[11]);    SQRADD2(a[8], a[10]);    SQRADD(a[9], a[9]); 
      COMBA_STORE(b[18]);

      /* output 19 */
      CARRY_FORWARD;
      SQRADD2(a[8], a[11]);    SQRADD2(a[9], a[10]); 
      COMBA_STORE(b[19]);

      /* output 20 */
      CARRY_FORWARD;
      SQRADD2(a[9], a[11]);    SQRADD(a[10], a[10]); 
      COMBA_STORE(b[20]);

      /* output 21 */
      CARRY_FORWARD;
      SQRADD2(a[10], a[11]); 
      COMBA_STORE(b[21]);

      /* output 22 */
      CARRY_FORWARD;
      SQRADD(a[11], a[11]); 
      COMBA_STORE(b[22]);
      COMBA_STORE2(b[23]);
      COMBA_FINI;

      B->size = 24;
      memcpy(B->val, b, 24 * sz);
      fp_clamp(B);
      break;

   case 13:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[9]); SQRADDAC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[10]); SQRADDAC(a[1], a[9]); SQRADDAC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[11]); SQRADDAC(a[1], a[10]); SQRADDAC(a[2], a[9]); SQRADDAC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[12]); SQRADDAC(a[1], a[11]); SQRADDAC(a[2], a[10]); SQRADDAC(a[3], a[9]); SQRADDAC(a[4], a[8]); SQRADDAC(a[5], a[7]); SQRADDDB; SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[12]); SQRADDAC(a[2], a[11]); SQRADDAC(a[3], a[10]); SQRADDAC(a[4], a[9]); SQRADDAC(a[5], a[8]); SQRADDAC(a[6], a[7]); SQRADDDB; 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[12]); SQRADDAC(a[3], a[11]); SQRADDAC(a[4], a[10]); SQRADDAC(a[5], a[9]); SQRADDAC(a[6], a[8]); SQRADDDB; SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[12]); SQRADDAC(a[4], a[11]); SQRADDAC(a[5], a[10]); SQRADDAC(a[6], a[9]); SQRADDAC(a[7], a[8]); SQRADDDB; 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
   SQRADDSC(a[4], a[12]); SQRADDAC(a[5], a[11]); SQRADDAC(a[6], a[10]); SQRADDAC(a[7], a[9]); SQRADDDB; SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);

      /* output 17 */
      CARRY_FORWARD;
   SQRADDSC(a[5], a[12]); SQRADDAC(a[6], a[11]); SQRADDAC(a[7], a[10]); SQRADDAC(a[8], a[9]); SQRADDDB; 
      COMBA_STORE(b[17]);

      /* output 18 */
      CARRY_FORWARD;
   SQRADDSC(a[6], a[12]); SQRADDAC(a[7], a[11]); SQRADDAC(a[8], a[10]); SQRADDDB; SQRADD(a[9], a[9]); 
      COMBA_STORE(b[18]);

      /* output 19 */
      CARRY_FORWARD;
   SQRADDSC(a[7], a[12]); SQRADDAC(a[8], a[11]); SQRADDAC(a[9], a[10]); SQRADDDB; 
      COMBA_STORE(b[19]);

      /* output 20 */
      CARRY_FORWARD;
      SQRADD2(a[8], a[12]);    SQRADD2(a[9], a[11]);    SQRADD(a[10], a[10]); 
      COMBA_STORE(b[20]);

      /* output 21 */
      CARRY_FORWARD;
      SQRADD2(a[9], a[12]);    SQRADD2(a[10], a[11]); 
      COMBA_STORE(b[21]);

      /* output 22 */
      CARRY_FORWARD;
      SQRADD2(a[10], a[12]);    SQRADD(a[11], a[11]); 
      COMBA_STORE(b[22]);

      /* output 23 */
      CARRY_FORWARD;
      SQRADD2(a[11], a[12]); 
      COMBA_STORE(b[23]);

      /* output 24 */
      CARRY_FORWARD;
      SQRADD(a[12], a[12]); 
      COMBA_STORE(b[24]);
      COMBA_STORE2(b[25]);
      COMBA_FINI;

      B->size = 26;
      memcpy(B->val, b, 26 * sz);
      fp_clamp(B);
      break;

   case 14:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[9]); SQRADDAC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[10]); SQRADDAC(a[1], a[9]); SQRADDAC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[11]); SQRADDAC(a[1], a[10]); SQRADDAC(a[2], a[9]); SQRADDAC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[12]); SQRADDAC(a[1], a[11]); SQRADDAC(a[2], a[10]); SQRADDAC(a[3], a[9]); SQRADDAC(a[4], a[8]); SQRADDAC(a[5], a[7]); SQRADDDB; SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[13]); SQRADDAC(a[1], a[12]); SQRADDAC(a[2], a[11]); SQRADDAC(a[3], a[10]); SQRADDAC(a[4], a[9]); SQRADDAC(a[5], a[8]); SQRADDAC(a[6], a[7]); SQRADDDB; 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[13]); SQRADDAC(a[2], a[12]); SQRADDAC(a[3], a[11]); SQRADDAC(a[4], a[10]); SQRADDAC(a[5], a[9]); SQRADDAC(a[6], a[8]); SQRADDDB; SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[13]); SQRADDAC(a[3], a[12]); SQRADDAC(a[4], a[11]); SQRADDAC(a[5], a[10]); SQRADDAC(a[6], a[9]); SQRADDAC(a[7], a[8]); SQRADDDB; 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[13]); SQRADDAC(a[4], a[12]); SQRADDAC(a[5], a[11]); SQRADDAC(a[6], a[10]); SQRADDAC(a[7], a[9]); SQRADDDB; SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);

      /* output 17 */
      CARRY_FORWARD;
   SQRADDSC(a[4], a[13]); SQRADDAC(a[5], a[12]); SQRADDAC(a[6], a[11]); SQRADDAC(a[7], a[10]); SQRADDAC(a[8], a[9]); SQRADDDB; 
      COMBA_STORE(b[17]);

      /* output 18 */
      CARRY_FORWARD;
   SQRADDSC(a[5], a[13]); SQRADDAC(a[6], a[12]); SQRADDAC(a[7], a[11]); SQRADDAC(a[8], a[10]); SQRADDDB; SQRADD(a[9], a[9]); 
      COMBA_STORE(b[18]);

      /* output 19 */
      CARRY_FORWARD;
   SQRADDSC(a[6], a[13]); SQRADDAC(a[7], a[12]); SQRADDAC(a[8], a[11]); SQRADDAC(a[9], a[10]); SQRADDDB; 
      COMBA_STORE(b[19]);

      /* output 20 */
      CARRY_FORWARD;
   SQRADDSC(a[7], a[13]); SQRADDAC(a[8], a[12]); SQRADDAC(a[9], a[11]); SQRADDDB; SQRADD(a[10], a[10]); 
      COMBA_STORE(b[20]);

      /* output 21 */
      CARRY_FORWARD;
   SQRADDSC(a[8], a[13]); SQRADDAC(a[9], a[12]); SQRADDAC(a[10], a[11]); SQRADDDB; 
      COMBA_STORE(b[21]);

      /* output 22 */
      CARRY_FORWARD;
      SQRADD2(a[9], a[13]);    SQRADD2(a[10], a[12]);    SQRADD(a[11], a[11]); 
      COMBA_STORE(b[22]);

      /* output 23 */
      CARRY_FORWARD;
      SQRADD2(a[10], a[13]);    SQRADD2(a[11], a[12]); 
      COMBA_STORE(b[23]);

      /* output 24 */
      CARRY_FORWARD;
      SQRADD2(a[11], a[13]);    SQRADD(a[12], a[12]); 
      COMBA_STORE(b[24]);

      /* output 25 */
      CARRY_FORWARD;
      SQRADD2(a[12], a[13]); 
      COMBA_STORE(b[25]);

      /* output 26 */
      CARRY_FORWARD;
      SQRADD(a[13], a[13]); 
      COMBA_STORE(b[26]);
      COMBA_STORE2(b[27]);
      COMBA_FINI;

      B->size = 28;
      memcpy(B->val, b, 28 * sz);
      fp_clamp(B);
      break;

   case 15:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[9]); SQRADDAC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[10]); SQRADDAC(a[1], a[9]); SQRADDAC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[11]); SQRADDAC(a[1], a[10]); SQRADDAC(a[2], a[9]); SQRADDAC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[12]); SQRADDAC(a[1], a[11]); SQRADDAC(a[2], a[10]); SQRADDAC(a[3], a[9]); SQRADDAC(a[4], a[8]); SQRADDAC(a[5], a[7]); SQRADDDB; SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[13]); SQRADDAC(a[1], a[12]); SQRADDAC(a[2], a[11]); SQRADDAC(a[3], a[10]); SQRADDAC(a[4], a[9]); SQRADDAC(a[5], a[8]); SQRADDAC(a[6], a[7]); SQRADDDB; 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[14]); SQRADDAC(a[1], a[13]); SQRADDAC(a[2], a[12]); SQRADDAC(a[3], a[11]); SQRADDAC(a[4], a[10]); SQRADDAC(a[5], a[9]); SQRADDAC(a[6], a[8]); SQRADDDB; SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[14]); SQRADDAC(a[2], a[13]); SQRADDAC(a[3], a[12]); SQRADDAC(a[4], a[11]); SQRADDAC(a[5], a[10]); SQRADDAC(a[6], a[9]); SQRADDAC(a[7], a[8]); SQRADDDB; 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[14]); SQRADDAC(a[3], a[13]); SQRADDAC(a[4], a[12]); SQRADDAC(a[5], a[11]); SQRADDAC(a[6], a[10]); SQRADDAC(a[7], a[9]); SQRADDDB; SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);

      /* output 17 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[14]); SQRADDAC(a[4], a[13]); SQRADDAC(a[5], a[12]); SQRADDAC(a[6], a[11]); SQRADDAC(a[7], a[10]); SQRADDAC(a[8], a[9]); SQRADDDB; 
      COMBA_STORE(b[17]);

      /* output 18 */
      CARRY_FORWARD;
   SQRADDSC(a[4], a[14]); SQRADDAC(a[5], a[13]); SQRADDAC(a[6], a[12]); SQRADDAC(a[7], a[11]); SQRADDAC(a[8], a[10]); SQRADDDB; SQRADD(a[9], a[9]); 
      COMBA_STORE(b[18]);

      /* output 19 */
      CARRY_FORWARD;
   SQRADDSC(a[5], a[14]); SQRADDAC(a[6], a[13]); SQRADDAC(a[7], a[12]); SQRADDAC(a[8], a[11]); SQRADDAC(a[9], a[10]); SQRADDDB; 
      COMBA_STORE(b[19]);

      /* output 20 */
      CARRY_FORWARD;
   SQRADDSC(a[6], a[14]); SQRADDAC(a[7], a[13]); SQRADDAC(a[8], a[12]); SQRADDAC(a[9], a[11]); SQRADDDB; SQRADD(a[10], a[10]); 
      COMBA_STORE(b[20]);

      /* output 21 */
      CARRY_FORWARD;
   SQRADDSC(a[7], a[14]); SQRADDAC(a[8], a[13]); SQRADDAC(a[9], a[12]); SQRADDAC(a[10], a[11]); SQRADDDB; 
      COMBA_STORE(b[21]);

      /* output 22 */
      CARRY_FORWARD;
   SQRADDSC(a[8], a[14]); SQRADDAC(a[9], a[13]); SQRADDAC(a[10], a[12]); SQRADDDB; SQRADD(a[11], a[11]); 
      COMBA_STORE(b[22]);

      /* output 23 */
      CARRY_FORWARD;
   SQRADDSC(a[9], a[14]); SQRADDAC(a[10], a[13]); SQRADDAC(a[11], a[12]); SQRADDDB; 
      COMBA_STORE(b[23]);

      /* output 24 */
      CARRY_FORWARD;
      SQRADD2(a[10], a[14]);    SQRADD2(a[11], a[13]);    SQRADD(a[12], a[12]); 
      COMBA_STORE(b[24]);

      /* output 25 */
      CARRY_FORWARD;
      SQRADD2(a[11], a[14]);    SQRADD2(a[12], a[13]); 
      COMBA_STORE(b[25]);

      /* output 26 */
      CARRY_FORWARD;
      SQRADD2(a[12], a[14]);    SQRADD(a[13], a[13]); 
      COMBA_STORE(b[26]);

      /* output 27 */
      CARRY_FORWARD;
      SQRADD2(a[13], a[14]); 
      COMBA_STORE(b[27]);

      /* output 28 */
      CARRY_FORWARD;
      SQRADD(a[14], a[14]); 
      COMBA_STORE(b[28]);
      COMBA_STORE2(b[29]);
      COMBA_FINI;

      B->size = 30;
      memcpy(B->val, b, 30 * sz);
      fp_clamp(B);
      break;

   case 16:
      a = A->val;
      COMBA_START; 

      /* clear carries */
      CLEAR_CARRY;

      /* output 0 */
      SQRADD(a[0],a[0]);
      COMBA_STORE(b[0]);

      /* output 1 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[1]); 
      COMBA_STORE(b[1]);

      /* output 2 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[2]);    SQRADD(a[1], a[1]); 
      COMBA_STORE(b[2]);

      /* output 3 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[3]);    SQRADD2(a[1], a[2]); 
      COMBA_STORE(b[3]);

      /* output 4 */
      CARRY_FORWARD;
      SQRADD2(a[0], a[4]);    SQRADD2(a[1], a[3]);    SQRADD(a[2], a[2]); 
      COMBA_STORE(b[4]);

      /* output 5 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB; 
      COMBA_STORE(b[5]);

      /* output 6 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]); 
      COMBA_STORE(b[6]);

      /* output 7 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB; 
      COMBA_STORE(b[7]);

      /* output 8 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[8]); SQRADDAC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]); 
      COMBA_STORE(b[8]);

      /* output 9 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[9]); SQRADDAC(a[1], a[8]); SQRADDAC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB; 
      COMBA_STORE(b[9]);

      /* output 10 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[10]); SQRADDAC(a[1], a[9]); SQRADDAC(a[2], a[8]); SQRADDAC(a[3], a[7]); SQRADDAC(a[4], a[6]); SQRADDDB; SQRADD(a[5], a[5]); 
      COMBA_STORE(b[10]);

      /* output 11 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[11]); SQRADDAC(a[1], a[10]); SQRADDAC(a[2], a[9]); SQRADDAC(a[3], a[8]); SQRADDAC(a[4], a[7]); SQRADDAC(a[5], a[6]); SQRADDDB; 
      COMBA_STORE(b[11]);

      /* output 12 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[12]); SQRADDAC(a[1], a[11]); SQRADDAC(a[2], a[10]); SQRADDAC(a[3], a[9]); SQRADDAC(a[4], a[8]); SQRADDAC(a[5], a[7]); SQRADDDB; SQRADD(a[6], a[6]); 
      COMBA_STORE(b[12]);

      /* output 13 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[13]); SQRADDAC(a[1], a[12]); SQRADDAC(a[2], a[11]); SQRADDAC(a[3], a[10]); SQRADDAC(a[4], a[9]); SQRADDAC(a[5], a[8]); SQRADDAC(a[6], a[7]); SQRADDDB; 
      COMBA_STORE(b[13]);

      /* output 14 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[14]); SQRADDAC(a[1], a[13]); SQRADDAC(a[2], a[12]); SQRADDAC(a[3], a[11]); SQRADDAC(a[4], a[10]); SQRADDAC(a[5], a[9]); SQRADDAC(a[6], a[8]); SQRADDDB; SQRADD(a[7], a[7]); 
      COMBA_STORE(b[14]);

      /* output 15 */
      CARRY_FORWARD;
   SQRADDSC(a[0], a[15]); SQRADDAC(a[1], a[14]); SQRADDAC(a[2], a[13]); SQRADDAC(a[3], a[12]); SQRADDAC(a[4], a[11]); SQRADDAC(a[5], a[10]); SQRADDAC(a[6], a[9]); SQRADDAC(a[7], a[8]); SQRADDDB; 
      COMBA_STORE(b[15]);

      /* output 16 */
      CARRY_FORWARD;
   SQRADDSC(a[1], a[15]); SQRADDAC(a[2], a[14]); SQRADDAC(a[3], a[13]); SQRADDAC(a[4], a[12]); SQRADDAC(a[5], a[11]); SQRADDAC(a[6], a[10]); SQRADDAC(a[7], a[9]); SQRADDDB; SQRADD(a[8], a[8]); 
      COMBA_STORE(b[16]);

      /* output 17 */
      CARRY_FORWARD;
   SQRADDSC(a[2], a[15]); SQRADDAC(a[3], a[14]); SQRADDAC(a[4], a[13]); SQRADDAC(a[5], a[12]); SQRADDAC(a[6], a[11]); SQRADDAC(a[7], a[10]); SQRADDAC(a[8], a[9]); SQRADDDB; 
      COMBA_STORE(b[17]);

      /* output 18 */
      CARRY_FORWARD;
   SQRADDSC(a[3], a[15]); SQRADDAC(a[4], a[14]); SQRADDAC(a[5], a[13]); SQRADDAC(a[6], a[12]); SQRADDAC(a[7], a[11]); SQRADDAC(a[8], a[10]); SQRADDDB; SQRADD(a[9], a[9]); 
      COMBA_STORE(b[18]);

      /* output 19 */
      CARRY_FORWARD;
   SQRADDSC(a[4], a[15]); SQRADDAC(a[5], a[14]); SQRADDAC(a[6], a[13]); SQRADDAC(a[7], a[12]); SQRADDAC(a[8], a[11]); SQRADDAC(a[9], a[10]); SQRADDDB; 
      COMBA_STORE(b[19]);

      /* output 20 */
      CARRY_FORWARD;
   SQRADDSC(a[5], a[15]); SQRADDAC(a[6], a[14]); SQRADDAC(a[7], a[13]); SQRADDAC(a[8], a[12]); SQRADDAC(a[9], a[11]); SQRADDDB; SQRADD(a[10], a[10]); 
      COMBA_STORE(b[20]);

      /* output 21 */
      CARRY_FORWARD;
   SQRADDSC(a[6], a[15]); SQRADDAC(a[7], a[14]); SQRADDAC(a[8], a[13]); SQRADDAC(a[9], a[12]); SQRADDAC(a[10], a[11]); SQRADDDB; 
      COMBA_STORE(b[21]);

      /* output 22 */
      CARRY_FORWARD;
   SQRADDSC(a[7], a[15]); SQRADDAC(a[8], a[14]); SQRADDAC(a[9], a[13]); SQRADDAC(a[10], a[12]); SQRADDDB; SQRADD(a[11], a[11]); 
      COMBA_STORE(b[22]);

      /* output 23 */
      CARRY_FORWARD;
   SQRADDSC(a[8], a[15]); SQRADDAC(a[9], a[14]); SQRADDAC(a[10], a[13]); SQRADDAC(a[11], a[12]); SQRADDDB; 
      COMBA_STORE(b[23]);

      /* output 24 */
      CARRY_FORWARD;
   SQRADDSC(a[9], a[15]); SQRADDAC(a[10], a[14]); SQRADDAC(a[11], a[13]); SQRADDDB; SQRADD(a[12], a[12]); 
      COMBA_STORE(b[24]);

      /* output 25 */
      CARRY_FORWARD;
   SQRADDSC(a[10], a[15]); SQRADDAC(a[11], a[14]); SQRADDAC(a[12], a[13]); SQRADDDB; 
      COMBA_STORE(b[25]);

      /* output 26 */
      CARRY_FORWARD;
      SQRADD2(a[11], a[15]);    SQRADD2(a[12], a[14]);    SQRADD(a[13], a[13]); 
      COMBA_STORE(b[26]);

      /* output 27 */
      CARRY_FORWARD;
      SQRADD2(a[12], a[15]);    SQRADD2(a[13], a[14]); 
      COMBA_STORE(b[27]);

      /* output 28 */
      CARRY_FORWARD;
      SQRADD2(a[13], a[15]);    SQRADD(a[14], a[14]); 
      COMBA_STORE(b[28]);

      /* output 29 */
      CARRY_FORWARD;
      SQRADD2(a[14], a[15]); 
      COMBA_STORE(b[29]);

      /* output 30 */
      CARRY_FORWARD;
      SQRADD(a[15], a[15]); 
      COMBA_STORE(b[30]);
      COMBA_STORE2(b[31]);
      COMBA_FINI;

      B->size = 32;
      memcpy(B->val, b, 32 * sz);
      fp_clamp(B);
      break;
	default:
		fp_sqr_comba(A,B);
}
}

#endif /* TFM_SMALL_SET */

/* $Source: /cvs/libtom/tomsfastmath/src/sqr/fp_sqr_comba_small_set.c,v $ */
/* $Revision: 1.1 $ */
/* $Date: 2007/02/15 00:31:32 $ */
