@* Scheduling function for the lattice siever.
Copyright (C) 2002 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
@(lasched.h@>=
u32_t *lasched(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t,u32_t,u32_t);
u32_t *lasched_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t,u32_t);

@
@c
#include <sys/types.h>
#include <limits.h>

#include "siever-config.h"
#include "lasched.h"
#include "../if.h"

#include <immintrin.h>
#include <stdlib.h>

#define L1_SIZE (1<<L1_BITS)
#define i_bits (I_bits-1)
#define n_i (1<<i_bits)

/* Amount by which to shift a larger number to provide space for a
   16-bit number. */

#define U16_SHIFT (CHAR_BIT*sizeof(u16_t))

#ifdef CONTIGUOUS_RI
#define RI_OFFSET FBsize
#define RI_OFFSET1 FBsize
#define RI_INCR 1
#else
#define RI_OFFSET 16
#define RI_OFFSET1 1
#define RI_INCR 2
#endif

#ifdef _MSC_VER
// so that I can read the code in MSVC without it being grayed out.
// It will not build in Visual studio.
#define AVX512_LASCHED
#endif


u32_t * ASM_ATTR lasched0(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched2(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched3(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched0nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched1nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched2nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched3nt(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);

u32_t * ASM_ATTR lasched0_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched1_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched2_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t * ASM_ATTR lasched3_1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);

u32_t *
lasched_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs,ot)
     u32_t *ri,*ij_ptr,*ij_ptr_ub,n1_j,**sched_ptr,fbi_offs,ot;
{
  u32_t ij,ij_ub;
  u32_t ot_mask,ot_tester;

  ij_ub=n1_j<<i_bits;

  switch(ot) {
  case 0:
    return lasched0_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  case 1:
    return lasched1_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  case 2:
    return lasched2_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  default:
    return lasched3_1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
  }
  return NULL;
}






u32_t *
lasched(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs,ot,FBsize)
     u32_t *ri,*ij_ptr,*ij_ptr_ub,n1_j,**sched_ptr,fbi_offs,ot,FBsize;
{
  u32_t ij,ij_ub;
  u32_t ot_mask,ot_tester;

  ij_ub=n1_j<<i_bits;

#ifdef SCHED_NT_BOUND
  if(ij_ub>SCHED_NT_BOUND)
    switch(ot) {
    case 0:
      return lasched0nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
    case 1:
      return lasched1nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
    case 2:
      return lasched2nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
    default:
      return lasched3nt(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
    }
#endif

#if !defined(AVX512_LASCHED)
	switch (ot) {
	case 0:
		return lasched0(ri, ij_ptr, ij_ptr_ub, n1_j, sched_ptr, fbi_offs);
	case 1:
		return lasched1(ri, ij_ptr, ij_ptr_ub, n1_j, sched_ptr, fbi_offs);
	case 2:
		return lasched2(ri, ij_ptr, ij_ptr_ub, n1_j, sched_ptr, fbi_offs);
	default:
		return lasched3(ri, ij_ptr, ij_ptr_ub, n1_j, sched_ptr, fbi_offs);
	}
#endif


	if (ot == 0)
	{

#if defined(AVX512_LASCHED)
		__m512i zero = _mm512_set1_epi32(0);
		__m512i vni = _mm512_set1_epi32(n_i);
		__m512i vni_m1 = _mm512_set1_epi32(n_i - 1);
		__m512i vL1m1 = _mm512_set1_epi32((L1_SIZE)-1);
		__m512i vfbi_offs = _mm512_set1_epi32(fbi_offs);
		__m512i vij_ub = _mm512_set1_epi32(ij_ub);
		__m512i vfbi_incr = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
		vfbi_offs = _mm512_add_epi32(vfbi_offs, vfbi_incr);
		vfbi_offs = _mm512_slli_epi32(vfbi_offs, U16_SHIFT);

		while (ij_ptr + 16 < ij_ptr_ub) {

			__m512i vri1 = _mm512_loadu_epi32(ri);
			__m512i vri2 = _mm512_loadu_epi32(ri + RI_OFFSET);
#ifndef CONTIGUOUS_RI
			// the compress/align is only necessary because ri is not
			// arranged well for this type of processing.  Next thing
			// to do might be to fix that.
			// After fixing it, it turns out not to matter.  The compresses/aligns
			// are apparently not a bottleneck.  In fact it is maybe slightly
			// slower this way... maybe because of the distance between ri
			// and ri+FBsize
			__m512i vri1a = _mm512_mask_compress_epi32(zero, 0x5555, vri1);
			__m512i vri1b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri1);
			__m512i vri2a = _mm512_mask_compress_epi32(zero, 0x5555, vri2);
			__m512i vri2b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri2);
			vri1 = _mm512_mask_alignr_epi32(vri1a, 0xff00, vri2a, zero, 8);
			vri2 = _mm512_mask_alignr_epi32(vri1b, 0xff00, vri2b, zero, 8);
#endif

			__m512i va = _mm512_sub_epi32(vni, _mm512_and_epi32(vri1, vni_m1));
			__m512i vb = _mm512_sub_epi32(vni, _mm512_and_epi32(vri2, vni_m1));
			__m512i vij = _mm512_load_epi32(ij_ptr);

			__mmask16 mij = _mm512_cmplt_epi32_mask(vij, vij_ub);
			u32_t memsched[16];
			u32_t memij[16];

#if 0
			while (mij == 0xffff) {

				__m512i vsched = _mm512_or_epi32(vfbi_offs, _mm512_and_epi32(vij, vL1m1));

				_mm512_store_epi32(memsched, vsched);
				_mm512_store_epi32(memij, vij);

				*(sched_ptr[memij[0] >> L1_BITS]++) = memsched[0];
				*(sched_ptr[memij[1] >> L1_BITS]++) = memsched[1];
				*(sched_ptr[memij[2] >> L1_BITS]++) = memsched[2];
				*(sched_ptr[memij[3] >> L1_BITS]++) = memsched[3];
				*(sched_ptr[memij[4] >> L1_BITS]++) = memsched[4];
				*(sched_ptr[memij[5] >> L1_BITS]++) = memsched[5];
				*(sched_ptr[memij[6] >> L1_BITS]++) = memsched[6];
				*(sched_ptr[memij[7] >> L1_BITS]++) = memsched[7];
				*(sched_ptr[memij[8] >> L1_BITS]++) = memsched[8];
				*(sched_ptr[memij[9] >> L1_BITS]++) = memsched[9];
				*(sched_ptr[memij[10] >> L1_BITS]++) = memsched[10];
				*(sched_ptr[memij[11] >> L1_BITS]++) = memsched[11];
				*(sched_ptr[memij[12] >> L1_BITS]++) = memsched[12];
				*(sched_ptr[memij[13] >> L1_BITS]++) = memsched[13];
				*(sched_ptr[memij[14] >> L1_BITS]++) = memsched[14];
				*(sched_ptr[memij[15] >> L1_BITS]++) = memsched[15];

				__m512i vi = _mm512_and_epi32(vij, vni_m1);
				__mmask16 mib = _mm512_mask_cmplt_epu32_mask(mij, vi, vb);
				__mmask16 mia = _mm512_mask_cmpge_epu32_mask(mij, vi, va);

				vij = _mm512_mask_add_epi32(vij, mib, vij, vri1);
				vij = _mm512_mask_add_epi32(vij, mia, vij, vri2);
				mij = _mm512_cmplt_epi32_mask(vij, _mm512_set1_epi32(ij_ub));
			}
#endif
			u32_t lastlane = 128;
			while (mij > 0) {
				__m512i vsched = _mm512_or_epi32(vfbi_offs, _mm512_and_epi32(vij, vL1m1));

				_mm512_store_epi32(memsched, vsched);
				_mm512_store_epi32(memij, vij); // _mm512_srli_epi32(vij, L1_BITS));

				// note: something to look into is how many buckets exist and
				// whether it might be more efficient to use compressstoreu
				// rather than this loop.  
				// The number of buckets varies, depending on input n1_j:
				// ij_ub = n1_j << i_bits.
				// Looking at 15e the number of buckets varied from 32 to well
				// over 1024, probably 4096.  As the number of buckets increases
				// obviously the density of hits decreases (fewer per bucket).
				// when the number of buckets is very small (maybe < 64) then
				// it might be better to use compressstoreu.  When it gets
				// very large (> 1024?) maybe scatter is an option.  In between
				// this routine is probably still best.
				// tried a bunch of things... compressstoreu doesn't help.
				// i also looked at iteration counts with different 
				// populations of mij.  A largeish fraction of the iterations
				// are done with a single set lane.  played with breaking out
				// when popcount(mij) == 1 and finishing with non-vector code
				// but that didn't give a speedup either (was close).

				u32_t m = mij;
				//bucket_hist[_mm_popcnt_u32(m)]++;
				int id;
				int c = 0;
				while (m > 0)
				{
					id = _tzcnt_u32(m);

					//if ((memij[id] >> L1_BITS) < 1024)
					//	bucket_hist[memij[id] >> L1_BITS]++;
					//else
					//	bucket_hist[1024]++;

					*(sched_ptr[memij[id] >> L1_BITS]++) = memsched[id];
					m = _blsr_u32(m);
					c++;
				}

				__m512i vi = _mm512_and_epi32(vij, vni_m1);
				__mmask16 mib = _mm512_mask_cmplt_epu32_mask(mij, vi, vb);
				__mmask16 mia = _mm512_mask_cmpge_epu32_mask(mij, vi, va);

				vij = _mm512_mask_add_epi32(vij, mib, vij, vri1);
				vij = _mm512_mask_add_epi32(vij, mia, vij, vri2);
				mij = _mm512_cmplt_epi32_mask(vij, vij_ub);

				if (c == 1)
				{
					lastlane = id;
					break;
				}

			}

			if (lastlane < 16)
			{
				u32_t fbioff = (fbi_offs + lastlane) << U16_SHIFT;
				_mm512_store_epi32(memij, vij);
				ij = memij[lastlane];
				u16_t a = n_i - (ri[lastlane*2 + 0] & (n_i - 1));
				u16_t b = n_i - (ri[lastlane*2 + 1] & (n_i - 1));
				u32_t ri0 = ri[lastlane * 2 + 0];
				u32_t ri1 = ri[lastlane * 2 + 1];
				while (ij < ij_ub) {
					u16_t i;
			
					*(sched_ptr[ij >> L1_BITS]++) = fbioff | (ij & (L1_SIZE - 1));
					i = ij & (n_i - 1);
					ij = (i < b) ? ij + ri0 : ij;
					ij = (i >= a) ? ij + ri1 : ij;
					//if (i < b)ij += ri[lastlane * 2 + 0];
					//if (i >= a)ij += ri[lastlane * 2 + 1];
				}
			
				memij[lastlane] = ij;
				vij = _mm512_load_epi32(memij);
			}

			ri += (RI_INCR * 16);
			_mm512_store_epi32(ij_ptr, _mm512_sub_epi32(vij, vij_ub));
			ij_ptr += 16;
			fbi_offs += 16;

			vfbi_offs = _mm512_set1_epi32(fbi_offs);
			vfbi_offs = _mm512_add_epi32(vfbi_offs, vfbi_incr);
			vfbi_offs = _mm512_slli_epi32(vfbi_offs, U16_SHIFT);
		}
#endif

		while (ij_ptr < ij_ptr_ub) {
			u16_t a, b;

			a = n_i - (ri[0] & (n_i - 1));
			b = n_i - (ri[RI_OFFSET1] & (n_i - 1));

			ij = *ij_ptr;

			while (ij < ij_ub) {
				u16_t i;

				*(sched_ptr[ij >> L1_BITS]++) = (fbi_offs << U16_SHIFT) | (ij & (L1_SIZE - 1));
				i = ij & (n_i - 1);
				if (i < b)ij += ri[0];
				if (i >= a)ij += ri[RI_OFFSET1];
			}
			ri += RI_INCR;
			*(ij_ptr++) = ij - ij_ub;
			fbi_offs++;
		}

	}
	else
	{
		ot_tester = (ot & 1) | ((ot & 2) << (i_bits - 1));
		ot_mask = n_i | 1;

#if defined(AVX512_LASCHED)
		__m512i vot_tester = _mm512_set1_epi32(ot_tester);
		__m512i vot_mask = _mm512_set1_epi32(ot_mask);
		__m512i zero = _mm512_set1_epi32(0);
		__m512i vni = _mm512_set1_epi32(n_i);
		__m512i vni_m1 = _mm512_set1_epi32(n_i - 1);
		__m512i vL1m1 = _mm512_set1_epi32(L1_SIZE - 1);
		__m512i vij_ub = _mm512_set1_epi32(ij_ub);
		__m512i vfbi_offs = _mm512_set1_epi32(fbi_offs);
		__m512i vfbi_incr = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
		vfbi_offs = _mm512_add_epi32(vfbi_offs, vfbi_incr);
		vfbi_offs = _mm512_slli_epi32(vfbi_offs, U16_SHIFT);

		while (ij_ptr + 16 < ij_ptr_ub) {

			__m512i vri1 = _mm512_loadu_epi32(ri);
			__m512i vri2 = _mm512_loadu_epi32(ri + RI_OFFSET);
#ifndef CONTIGUOUS_RI
			__m512i vri1a = _mm512_mask_compress_epi32(zero, 0x5555, vri1);
			__m512i vri1b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri1);
			__m512i vri2a = _mm512_mask_compress_epi32(zero, 0x5555, vri2);
			__m512i vri2b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri2);
			vri1 = _mm512_mask_alignr_epi32(vri1a, 0xff00, vri2a, zero, 8);
			vri2 = _mm512_mask_alignr_epi32(vri1b, 0xff00, vri2b, zero, 8);
#endif
			__m512i va = _mm512_sub_epi32(vni, _mm512_and_epi32(vri1, vni_m1));
			__m512i vb = _mm512_sub_epi32(vni, _mm512_and_epi32(vri2, vni_m1));
			u32_t memsched[16];
			u32_t memij[16];
			__m512i vij = zero;

			__mmask16 mri0 = _mm512_cmpeq_epu32_mask(vot_tester, _mm512_and_epi32(vri1, vot_mask));
			__mmask16 mri1 = _mm512_cmpeq_epu32_mask(
				_mm512_xor_epi32(vot_tester, vni), _mm512_and_epi32(vri2, vot_mask));
			__mmask16 mri2 = _mm512_cmple_epu32_mask(_mm512_and_epi32(vri1, vni_m1),
				_mm512_and_epi32(vri2, vni_m1));
			__mmask16 mri3 = _mm512_cmple_epu32_mask(vri1, vri2);
			__mmask16 mri4 = _mm512_cmpeq_epu32_mask(_mm512_and_epi32(vri1, vni_m1),
				_mm512_and_epi32(vri2, vni_m1));

			vij = _mm512_mask_mov_epi32(vij, mri0, vri1);
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & mri1, vri2);
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & (~mri1) & mri2 & mri3 & mri4, _mm512_sub_epi32(vri2, vri1));
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & (~mri1) & mri2 & mri3 & (~mri4), vni);
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & (~mri1) & ((~mri2) | (~mri3)), _mm512_add_epi32(vri2, vri1));
			vij = _mm512_srli_epi32(_mm512_add_epi32(vij, _mm512_andnot_epi32(vot_tester, vni)), 1);

			__mmask16 mij = _mm512_cmplt_epi32_mask(vij, vij_ub);

#if 0
			while (mij == 0xffff) {
				__m512i vsched = _mm512_or_epi32(vfbi_offs, _mm512_and_epi32(vij, vL1m1));

				_mm512_store_epi32(memsched, vsched);
				_mm512_store_epi32(memij, vij);

				*(sched_ptr[memij[0] >> L1_BITS]++) = memsched[0];
				*(sched_ptr[memij[1] >> L1_BITS]++) = memsched[1];
				*(sched_ptr[memij[2] >> L1_BITS]++) = memsched[2];
				*(sched_ptr[memij[3] >> L1_BITS]++) = memsched[3];
				*(sched_ptr[memij[4] >> L1_BITS]++) = memsched[4];
				*(sched_ptr[memij[5] >> L1_BITS]++) = memsched[5];
				*(sched_ptr[memij[6] >> L1_BITS]++) = memsched[6];
				*(sched_ptr[memij[7] >> L1_BITS]++) = memsched[7];
				*(sched_ptr[memij[8] >> L1_BITS]++) = memsched[8];
				*(sched_ptr[memij[9] >> L1_BITS]++) = memsched[9];
				*(sched_ptr[memij[10] >> L1_BITS]++) = memsched[10];
				*(sched_ptr[memij[11] >> L1_BITS]++) = memsched[11];
				*(sched_ptr[memij[12] >> L1_BITS]++) = memsched[12];
				*(sched_ptr[memij[13] >> L1_BITS]++) = memsched[13];
				*(sched_ptr[memij[14] >> L1_BITS]++) = memsched[14];
				*(sched_ptr[memij[15] >> L1_BITS]++) = memsched[15];

				__m512i vi = _mm512_and_epi32(vij, vni_m1);
				__mmask16 mib = _mm512_mask_cmplt_epu32_mask(mij, vi, vb);
				__mmask16 mia = _mm512_mask_cmpge_epu32_mask(mij, vi, va);

				vij = _mm512_mask_add_epi32(vij, mib, vij, vri1);
				vij = _mm512_mask_add_epi32(vij, mia, vij, vri2);
				mij = _mm512_cmplt_epi32_mask(vij, _mm512_set1_epi32(ij_ub));
			}
#endif
			u32_t lastlane = 128;
			while (mij > 0) {
				__m512i vsched = _mm512_or_epi32(vfbi_offs, _mm512_and_epi32(vij, vL1m1));

				_mm512_store_epi32(memsched, vsched);
				_mm512_store_epi32(memij, vij);

				u32_t m = mij;
				//bucket_hist[_mm_popcnt_u32(m)]++;

				int id;
				int c = 0;
				while (m > 0)
				{
					id = _tzcnt_u32(m);

					//if ((memij[id] >> L1_BITS) < 1024)
					//	bucket_hist[memij[id] >> L1_BITS]++;
					//else
					//	bucket_hist[1024]++;


					*(sched_ptr[memij[id] >> L1_BITS]++) = memsched[id];
					m = _blsr_u32(m);
					c++;
				}

				__m512i vi = _mm512_and_epi32(vij, vni_m1);
				__mmask16 mib = _mm512_mask_cmplt_epu32_mask(mij, vi, vb);
				__mmask16 mia = _mm512_mask_cmpge_epu32_mask(mij, vi, va);

				vij = _mm512_mask_add_epi32(vij, mib, vij, vri1);
				vij = _mm512_mask_add_epi32(vij, mia, vij, vri2);
				mij = _mm512_cmplt_epi32_mask(vij, vij_ub);

				if (c == 1)
				{
					lastlane = id;
					break;
				}
			}


			if (lastlane < 16)
			{
				u32_t fbioff = (fbi_offs + lastlane) << U16_SHIFT;
				_mm512_store_epi32(memij, vij);
				ij = memij[lastlane];
				u16_t a = n_i - (ri[lastlane * 2 + 0] & (n_i - 1));
				u16_t b = n_i - (ri[lastlane * 2 + 1] & (n_i - 1));
				u32_t ri0 = ri[lastlane * 2 + 0];
				u32_t ri1 = ri[lastlane * 2 + 1];
				while (ij < ij_ub) {
					u16_t i;

					*(sched_ptr[ij >> L1_BITS]++) = fbioff | (ij & (L1_SIZE - 1));
					i = ij & (n_i - 1);
					ij = (i < b) ? ij + ri0 : ij;
					ij = (i >= a) ? ij + ri1 : ij;
				}

				memij[lastlane] = ij;
				vij = _mm512_load_epi32(memij);
			}

			ri += (RI_INCR * 16);
			_mm512_store_epi32(ij_ptr, _mm512_sub_epi32(vij, vij_ub));
			ij_ptr += 16;
			fbi_offs += 16;

			vfbi_offs = _mm512_set1_epi32(fbi_offs);
			vfbi_offs = _mm512_add_epi32(vfbi_offs, vfbi_incr);
			vfbi_offs = _mm512_slli_epi32(vfbi_offs, U16_SHIFT);
		}
#endif

		while (ij_ptr < ij_ptr_ub) {
			u16_t a, b;

			a = n_i - (ri[0] & (n_i - 1));
			b = n_i - (ri[RI_OFFSET1] & (n_i - 1));

#line 102 "lasched.w"

			{
				ij = 0;
				if ((ri[0] & ot_mask) == ot_tester)ij = ri[0];
				else {
					if ((ri[RI_OFFSET1] & ot_mask) == (ot_tester ^ n_i))ij = ri[RI_OFFSET1];
					else {
						if ((ri[0] & (n_i - 1)) <= (ri[RI_OFFSET1] & (n_i - 1)) && ri[0] <= ri[RI_OFFSET1]) {


							if ((ri[0] & (n_i - 1)) == (ri[RI_OFFSET1] & (n_i - 1)))ij = ri[RI_OFFSET1] - ri[0];
							else ij = n_i;
							if (ot != 2)
								Schlendrian("Exceptional situation for oddness type %u ?\n",
									ot);
						}
						else ij = ri[0] + ri[RI_OFFSET1];
					}
				}
				ij = (ij + ((~ot_tester) & n_i)) / 2;
			}/*:3*/
#line 85 "lasched.w"

			while (ij < ij_ub) {
				u16_t i;

				*(sched_ptr[ij >> L1_BITS]++) = (fbi_offs << U16_SHIFT) | (ij & (L1_SIZE - 1));
				i = ij & (n_i - 1);
				if (i < b)ij += ri[0];
				if (i >= a)ij += ri[RI_OFFSET1];
			}
			ri += RI_INCR;
			*(ij_ptr++) = ij - ij_ub;
			fbi_offs++;
		}

	}

	//{
	//	int i;
	//	for (i = 0; i <= 16; i++)
	//	{
	//		if ((i & 15) == 0)
	//			printf("\n%d: %u ", i, bucket_hist[i]);
	//		else
	//			printf("%u ", bucket_hist[i]);
	//		bucket_hist[i] = 0;
	//	}
	//	printf("\n");
	//	
	//}

	return ri;
}

@
@<Calculate first sieving event from |ri|@>=
{
  ij=0;
  if( (ri[0]&ot_mask) == ot_tester ) ij=ri[0];
  else {
    if( (ri[1]&ot_mask) == (ot_tester^n_i) ) ij=ri[1];
    else {
      if((ri[0]&(n_i-1))<=(ri[1]&(n_i-1)) && ri[0]<=ri[1]) {
	/* This corresponds to the line
	   |if(b+c<=A && s<=t)| in recurrence6.w */
	if((ri[0]&(n_i-1))==(ri[1]&(n_i-1))) ij=ri[1]-ri[0];
	else ij=n_i;
	if(ot != 2)
	  Schlendrian("Exceptional situation for oddness type %u ?\n",
		      ot);
      }
      else ij=ri[0]+ri[1];
    }
  }
  ij=(ij+((~ot_tester)&n_i))/2;
}
