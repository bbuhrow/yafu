/*2:*/
#line 14 "medsched.w"

#include <sys/types.h> 
#include <limits.h> 

#include "siever-config.h"
#include "medsched.h"
#include "../if.h"

#include <immintrin.h>

#define L1_SIZE (1<<L1_BITS)
#define i_bits (I_bits-1)
#define n_i (1<<i_bits)

#define U16_SHIFT (CHAR_BIT*sizeof(u16_t))

#ifdef _MSC_VER
// so that I can read the code in MSVC without it being grayed out.
// It will not build in Visual studio.
#define AVX512_LASCHED
#endif

u32_t * medsched0(u32_t*,u32_t*,u32_t*,u32_t**,u32_t);

u32_t*
medsched(u32_t* ri, u32_t* ij_ptr, u32_t* ij_ptr_ub, u32_t**sched_ptr, u32_t fbi_offs, u32_t ot)
{
	u32_t ij;
	u32_t ot_mask, ot_tester;
	u32_t* sched;

	sched = *sched_ptr;

#ifndef AVX512_LASCHED
	if (ot != 0) {
		ot_tester = (ot & 1) | ((ot & 2) << (i_bits - 1));
		ot_mask = n_i | 1;
	}
	else return medsched0(ri, ij_ptr, ij_ptr_ub, sched_ptr, fbi_offs);
#endif

	if (ot == 0)
	{

#if defined(AVX512_LASCHED)
		__m512i zero = _mm512_set1_epi32(0);
		__m512i vni = _mm512_set1_epi32(n_i);
		__m512i vni_m1 = _mm512_set1_epi32(n_i - 1);
		__m512i vL1m1 = _mm512_set1_epi32((L1_SIZE)-1);
		__m512i vfbi_offs = _mm512_set1_epi32(fbi_offs);
		__m512i vfbi_incr = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
		vfbi_offs = _mm512_add_epi32(vfbi_offs, vfbi_incr);
		vfbi_offs = _mm512_slli_epi32(vfbi_offs, U16_SHIFT);

		while (ij_ptr + 16 < ij_ptr_ub) {

			__m512i vri1 = _mm512_loadu_epi32(ri);
			__m512i vri2 = _mm512_loadu_epi32(ri + 16);
			// the compress/align is only necessary because ri is not
			// arranged well for this type of processing.  Next thing
			// to do might be to fix that.
			__m512i vri1a = _mm512_mask_compress_epi32(zero, 0x5555, vri1);
			__m512i vri1b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri1);
			__m512i vri2a = _mm512_mask_compress_epi32(zero, 0x5555, vri2);
			__m512i vri2b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri2);
			__m512i vria = _mm512_mask_alignr_epi32(vri1a, 0xff00, vri2a, zero, 8);
			__m512i vrib = _mm512_mask_alignr_epi32(vri1b, 0xff00, vri2b, zero, 8);
			__m512i va = _mm512_sub_epi32(vni, _mm512_and_epi32(vria, vni_m1));
			__m512i vb = _mm512_sub_epi32(vni, _mm512_and_epi32(vrib, vni_m1));
			__m512i vij = _mm512_load_epi32(ij_ptr);

			__mmask16 mij = _mm512_cmplt_epi32_mask(vij, _mm512_set1_epi32(L1_SIZE));
			u32_t memsched[16];
			u32_t memij[16];

			while (mij > 0) {

				__m512i vsched = _mm512_or_epi32(vfbi_offs, vij);

				_mm512_store_epi32(memsched, vsched);

				u32_t m = mij;
				while (m > 0)
				{
					int id = _tzcnt_u32(m);
					*(sched++) = memsched[id];
					m = _blsr_u32(m);
				}

				__m512i vi = _mm512_and_epi32(vij, vni_m1);
				__mmask16 mib = _mm512_mask_cmplt_epu32_mask(mij, vi, vb);
				__mmask16 mia = _mm512_mask_cmpge_epu32_mask(mij, vi, va);

				vij = _mm512_mask_add_epi32(vij, mib, vij, vria);
				vij = _mm512_mask_add_epi32(vij, mia, vij, vrib);
				mij = _mm512_cmplt_epi32_mask(vij, _mm512_set1_epi32(L1_SIZE));
			}

			ri += 32;
			_mm512_store_epi32(ij_ptr, _mm512_sub_epi32(vij, _mm512_set1_epi32(L1_SIZE)));
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
			b = n_i - (ri[1] & (n_i - 1));
			ij = *ij_ptr;

			while (ij < L1_SIZE) {
				u16_t i;

				*(sched++) = (fbi_offs << U16_SHIFT) | ij;
				i = ij & (n_i - 1);
				if (i < b)ij += ri[0];
				if (i >= a)ij += ri[1];
			}
			ri += 2;
			*(ij_ptr++) = ij - L1_SIZE;
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
		__m512i vfbi_offs = _mm512_set1_epi32(fbi_offs);
		__m512i vfbi_incr = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
		vfbi_offs = _mm512_add_epi32(vfbi_offs, vfbi_incr);
		vfbi_offs = _mm512_slli_epi32(vfbi_offs, U16_SHIFT);

		while (ij_ptr + 16 < ij_ptr_ub) {

			__m512i vri1 = _mm512_loadu_epi32(ri);
			__m512i vri2 = _mm512_loadu_epi32(ri + 16);
			__m512i vri1a = _mm512_mask_compress_epi32(zero, 0x5555, vri1);
			__m512i vri1b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri1);
			__m512i vri2a = _mm512_mask_compress_epi32(zero, 0x5555, vri2);
			__m512i vri2b = _mm512_mask_compress_epi32(zero, 0xaaaa, vri2);
			__m512i vria = _mm512_mask_alignr_epi32(vri1a, 0xff00, vri2a, zero, 8);
			__m512i vrib = _mm512_mask_alignr_epi32(vri1b, 0xff00, vri2b, zero, 8);
			__m512i va = _mm512_sub_epi32(vni, _mm512_and_epi32(vria, vni_m1));
			__m512i vb = _mm512_sub_epi32(vni, _mm512_and_epi32(vrib, vni_m1));
			u32_t memsched[16];
			u32_t memij[16];
			__m512i vij = zero;

			__mmask16 mri0 = _mm512_cmpeq_epu32_mask(vot_tester, _mm512_and_epi32(vria, vot_mask));
			__mmask16 mri1 = _mm512_cmpeq_epu32_mask(
				_mm512_xor_epi32(vot_tester, vni), _mm512_and_epi32(vrib, vot_mask));
			__mmask16 mri2 = _mm512_cmple_epu32_mask(_mm512_and_epi32(vria, vni_m1),
				_mm512_and_epi32(vrib, vni_m1));
			__mmask16 mri3 = _mm512_cmple_epu32_mask(vria, vrib);
			__mmask16 mri4 = _mm512_cmpeq_epu32_mask(_mm512_and_epi32(vria, vni_m1),
				_mm512_and_epi32(vrib, vni_m1));

			vij = _mm512_mask_mov_epi32(vij, mri0, vria);
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & mri1, vrib);
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & (~mri1) & mri2 & mri3 & mri4, _mm512_sub_epi32(vrib, vria));
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & (~mri1) & mri2 & mri3 & (~mri4), vni);
			vij = _mm512_mask_mov_epi32(vij, (~mri0) & (~mri1) & ((~mri2) | (~mri3)), _mm512_add_epi32(vrib, vria));
			vij = _mm512_srli_epi32(_mm512_add_epi32(vij, _mm512_andnot_epi32(vot_tester, vni)), 1);

			__mmask16 mij = _mm512_cmplt_epi32_mask(vij, _mm512_set1_epi32(L1_SIZE));

			while (mij > 0) {
				__m512i vsched = _mm512_or_epi32(vfbi_offs, vij);

				_mm512_store_epi32(memsched, vsched);

				u32_t m = mij;
				while (m > 0)
				{
					int id = _tzcnt_u32(m);
					*(sched++) = memsched[id];
					m = _blsr_u32(m);
				}

				__m512i vi = _mm512_and_epi32(vij, vni_m1);
				__mmask16 mib = _mm512_mask_cmplt_epu32_mask(mij, vi, vb);
				__mmask16 mia = _mm512_mask_cmpge_epu32_mask(mij, vi, va);

				vij = _mm512_mask_add_epi32(vij, mib, vij, vria);
				vij = _mm512_mask_add_epi32(vij, mia, vij, vrib);
				mij = _mm512_cmplt_epi32_mask(vij, _mm512_set1_epi32(L1_SIZE));
			}

			ri += 32;
			_mm512_store_epi32(ij_ptr, _mm512_sub_epi32(vij, _mm512_set1_epi32(L1_SIZE)));
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
			b = n_i - (ri[1] & (n_i - 1));

			{
				ij = 0;
				if ((ri[0] & ot_mask) == ot_tester)ij = ri[0];
				else {
					if ((ri[1] & ot_mask) == (ot_tester ^ n_i))ij = ri[1];
					else {
						if ((ri[0] & (n_i - 1)) <= (ri[1] & (n_i - 1)) && ri[0] <= ri[1]) {


							if ((ri[0] & (n_i - 1)) == (ri[1] & (n_i - 1)))ij = ri[1] - ri[0];
							else ij = n_i;
							if (ot != 2)
								Schlendrian("Exceptional situation for oddness type %u ?\n",
									ot);
						}
						else ij = ri[0] + ri[1];
					}
				}
				ij = (ij + ((~ot_tester) & n_i)) / 2;
			}

			while (ij < L1_SIZE) {
				u16_t i;

				*(sched++) = (fbi_offs << U16_SHIFT) | ij;
				i = ij & (n_i - 1);
				if (i < b)ij += ri[0];
				if (i >= a)ij += ri[1];
			}
			ri += 2;
			*(ij_ptr++) = ij - L1_SIZE;
			fbi_offs++;
		}

	}
	
	*sched_ptr = sched;
	return ri;
}

/*:2*/
