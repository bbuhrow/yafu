/* recurrence6.c
  By Jens Franke.
  6/13/04: Hacked up for use in GGNFS by Chris Monico.
  9/30/22: Vector AVX512 code contributed by Ben Buhrow

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/


#include <sys/types.h>
#include <limits.h>

#include "asm/siever-config.h"
#include "if.h"
#include "recurrence6.h"
#include <immintrin.h>

static u32_t A,A_bits,ub;

void rec_info_init(u32_t A1,u32_t ub1)
{
  u32_t aa;

  if(A1> USHRT_MAX)complain("Recurrence init: A=%u exceeds %u\n",A1,USHRT_MAX);
  A= A1;
  if(A%2!=0)complain("rec_info_init with odd range %u\n",A1);
  if(ub1> (USHRT_MAX / 4 + 1))
    complain("Recurrence init: ub=%u exceeds %u\n",ub1,(USHRT_MAX / 4 + 1));
  if(ub1<2)
    complain("Recurrence time %u too small\n",ub1);
  ub= ub1;
  for(aa= 1,A_bits= 0;A_bits<=CHAR_BIT*sizeof(u32_t);A_bits++,aa*= 2)
    if(aa==A)break;
  if(A_bits> CHAR_BIT*sizeof(u32_t))
    complain("rec_info_init: A=%u not a power of two\n",A);
}

u32_t get_recurrence_info(u32_t*res_ptr,u32_t p,u32_t r)
{
    u32_t b,c,s,t;

    if (r == 0) {
        b = 1;
        s = 2 * ub + 1;
        c = A - 3;
        t = 2 * ub;
    
        goto done;
    }
    if (r >= p) {
    
        b = 1;
        s = 0;
        c = A - 1;
        t = 2 * ub + 1;
    
        goto done;
    }

#if 0 //def HAVE_ASM_GETBC
    asm_getbc(r,p,A,&b,&s,&c,&t);
#else
    {
        b = r;
        s = 1;
        c = p;
        t = 0;
        for (;;) {
            u32_t k;

            if (b < A) {

                k = (c - A) / b + 1;
                t += k * s;
                c -= k * b;
                break;
            }
            k = c / b;
            c = c % b;
            t += k * s;
            if (c < A) {

                k = (b - A) / c + 1;
                s += k * t;
                b -= k * c;
                break;
            }
            k = b / c;
            b = b % c;
            s += k * t;
        }
    }
#endif

#if 0
{
  u32_t b1,c1,s1,t1;

  asm_getbc(r,p,A,&b1,&s1,&c1,&t1);
  if(b1!=b||c1!=c||s1!=s||t1!=t)Schlendrian("Asm ri\n");
}
#endif

{
    i32_t d;

    d = t - s;
    s = (s <= 2 * ub ? s : 2 * ub + 2 - (s & 1));
    if (t > 2 * ub) {
        if (d < 0 || d> 2 * ub)
            d = 2 * ub + 2 + (s & 1) - (t & 1);
        t = s + d;
    }
}

done:
    {
      res_ptr[0]= (s<<A_bits)+b;
      res_ptr[1]= (t<<A_bits)-c;
    }

#if 0
{
u32_t z[2];
u32_t have_it[3]= {0,0,0};

z[0]= A+b;
z[1]= s;
{
u32_t ot;
#ifdef ULONG_RI
u32_t*xx;
#else
u16_t*xx;
#endif

ot= (z[0]&1)+2*(z[1]&1);
if(ot==1)xx= x1;
else if(ot==2)xx= x2;
else xx= x3;
have_it[ot-1]= 1;
#ifdef ULONG_RI
*xx= ((z[1]/2)<<A_bits)|(z[0]/2);
#else
xx[0]= z[0]/2;
xx[1]= z[1]/2;
#endif
}

if(b+c<=A&&s<=t){
if(b+c<A){
z[0]= A;
z[1]= 1;
}else{
z[0]= 0;
z[1]= t-s;
}
}else{
z[0]= A+b-c;
z[1]= t+s;
}
#if 0
z[1]= (z[1]<=2*ub?z[1]:2*ub+2-(z[1]&1));
#else
z[1]= (z[1]<=USHRT_MAX?z[1]:2*ub+2-(z[1]&1));
#endif
{
u32_t ot;
#ifdef ULONG_RI
u32_t*xx;
#else
u16_t*xx;
#endif

ot= (z[0]&1)+2*(z[1]&1);
if(ot==1)xx= x1;
else if(ot==2)xx= x2;
else xx= x3;
have_it[ot-1]= 1;
#ifdef ULONG_RI
*xx= ((z[1]/2)<<A_bits)|(z[0]/2);
#else
xx[0]= z[0]/2;
xx[1]= z[1]/2;
#endif

z[0]= A-c;
z[1]= t;

{
u32_t ot;
#ifdef ULONG_RI
u32_t*xx;
#else
u16_t*xx;
#endif

ot= (z[0]&1)+2*(z[1]&1);
if(ot==1)xx= x1;
else if(ot==2)xx= x2;
else xx= x3;
have_it[ot-1]= 1;
#ifdef ULONG_RI
*xx= ((z[1]/2)<<A_bits)|(z[0]/2);
#else
xx[0]= z[0]/2;
xx[1]= z[1]/2;
#endif

if(have_it[0]==0)Schlendrian("???");
}

#endif

    return 2;
}


u32_t get_recurrence_info_16(u32_t* res_ptr, __m512i p, __m512i r)
{
    __m512i zero = _mm512_setzero_epi32();
    __m512i b = zero, c = zero, s = zero, t = zero;
    __m512i one = _mm512_set1_epi32(1);
    __m512i vA = _mm512_set1_epi32(A);
    __m512i ub2 = _mm512_slli_epi32(_mm512_set1_epi32(ub), 1);
    __mmask16 r0 = _mm512_cmpeq_epi32_mask(r, zero);
    __mmask16 rp = _mm512_cmpge_epu32_mask(r, p);

    //if (r == 0) {
    //	b = 1;
    //	s = 2 * ub + 1;
    //	c = A - 3;
    //	t = 2 * ub;
    //
    //	goto done;
    //}
    b = _mm512_mask_set1_epi32(b, r0, 1);
    t = _mm512_mask_mov_epi32(t, r0, ub2);
    s = _mm512_mask_add_epi32(s, r0, ub2, one);
    c = _mm512_mask_sub_epi32(c, r0, vA, _mm512_set1_epi32(3));

    //if (r >= p) {
    //
    //    b = 1;
    //    s = 0;
    //    c = A - 1;
    //    t = 2 * ub + 1;
    //
    //    goto done;
    //}

    b = _mm512_mask_set1_epi32(b, rp, 1);
    t = _mm512_mask_add_epi32(t, rp, ub2, one);
    s = _mm512_mask_mov_epi32(s, rp, zero);
    c = _mm512_mask_sub_epi32(c, rp, vA, one);
    
    __mmask16 skipped = r0 | rp;
    __mmask16 m = ~skipped;

    {
        //b = r;
        //s = 1;
        //c = p;
        //t = 0;
        b = _mm512_mask_mov_epi32(b, m, r);
        s = _mm512_mask_mov_epi32(s, m, one);
        c = _mm512_mask_mov_epi32(c, m, p);
        t = _mm512_mask_mov_epi32(t, m, zero);

        for (;;) {
            __m512i k;
            __mmask16 bA = _mm512_cmplt_epu32_mask(b, vA);

            //if (b < A) {
            //
            //    k = (c - A) / b + 1;
            //    t += k * s;
            //    c -= k * b;
            //    break;
            //}

            k = _mm512_mask_add_epi32(k, m & bA, _mm512_div_epu32(_mm512_sub_epi32(c, vA), b), one);
            t = _mm512_mask_add_epi32(t, m & bA, t, _mm512_mullo_epi32(k, s));
            c = _mm512_mask_sub_epi32(c, m & bA, c, _mm512_mullo_epi32(k, b));

            m = m & (~bA);
            if (m == 0)
                break;

            //k = c / b;
            //c = c % b;
            //t += k * s;

            k = _mm512_mask_div_epu32(k, m, c, b);
            c = _mm512_mask_rem_epu32(c, m, c, b);
            t = _mm512_mask_add_epi32(t, m, t, _mm512_mullo_epi32(k, s));

            __mmask16 cA = _mm512_cmplt_epu32_mask(c, vA);
            //if (c < A) {
            //
            //    k = (b - A) / c + 1;
            //    s += k * t;
            //    b -= k * c;
            //    break;
            //}

            k = _mm512_mask_add_epi32(k, m & cA, _mm512_div_epu32(_mm512_sub_epi32(b, vA), c), one);
            s = _mm512_mask_add_epi32(s, m & cA, s, _mm512_mullo_epi32(k, t));
            b = _mm512_mask_sub_epi32(b, m & cA, b, _mm512_mullo_epi32(k, c));

            m = m & (~cA);
            if (m == 0)
                break;

            //k = b / c;
            //b = b % c;
            //s += k * t;

            k = _mm512_mask_div_epu32(k, m, b, c);
            b = _mm512_mask_rem_epu32(b, m, b, c);
            s = _mm512_mask_add_epi32(s, m, s, _mm512_mullo_epi32(k, t));
        }
    }


    //{
    //    i32_t d;
    //
    //    d = t - s;
    //    s = (s <= 2 * ub ? s : 2 * ub + 2 - (s & 1));
    //    if (t > 2 * ub) {
    //        if (d < 0 || d> 2 * ub)
    //          d = 2 * ub + 2 + (s & 1) - (t & 1);
    //        t = s + d;
    //    }
    //}

    {
        m = ~skipped;
        __m512i d = _mm512_mask_sub_epi32(zero, m, t, s);
        __mmask16 s_ub2 = _mm512_cmpgt_epu32_mask(s, ub2);
        __m512i two = _mm512_set1_epi32(2);

        s = _mm512_mask_mov_epi32(s, m & s_ub2,
            _mm512_sub_epi32(_mm512_add_epi32(ub2, two), _mm512_and_epi32(s, one)));

        __mmask16 t_ub2 = _mm512_cmpgt_epu32_mask(t, ub2);
        __mmask16 d_0 = _mm512_cmplt_epi32_mask(d, zero);
        __mmask16 d_ub2 = _mm512_cmpgt_epi32_mask(d, ub2);

        d = _mm512_mask_sub_epi32(d, m & (d_0 | d_ub2),
            _mm512_add_epi32(_mm512_add_epi32(ub2, two), _mm512_and_epi32(s, one)), _mm512_and_epi32(t, one));

        t = _mm512_mask_add_epi32(t, m & t_ub2, s, d);
    }

//done:

    int i;
    u32_t sm[16], bm[16], tm[16], cm[16];
    _mm512_storeu_epi32(sm, s);
    _mm512_storeu_epi32(bm, b);
    _mm512_storeu_epi32(tm, t);
    _mm512_storeu_epi32(cm, c);


    for (i = 0; i < 16; i++)
    {
        res_ptr[2*i+0] = (sm[i] << A_bits) + bm[i];
        res_ptr[2*i+1] = (tm[i] << A_bits) - cm[i];
    }

    return 32;
}

