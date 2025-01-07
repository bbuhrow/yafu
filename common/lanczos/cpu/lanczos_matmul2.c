/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#include "lanczos_cpu.h"

	/* code for handling transposed matrix multiplies 
	   when the matrix is in packed format */

/*-------------------------------------------------------------------*/
void mul_trans_one_med_block(packed_block_t *curr_block,
			v_t *curr_row, v_t *curr_b) {

	uint16 *entries = curr_block->d.med_entries;

	while (1) {
		v_t t;
#if defined(GCC_ASM64X)
		uint64 i = 0;
		uint64 row = entries[0];
		uint64 count = entries[1];
#else
		uint32 i = 0;
		uint32 row = entries[0];
		uint32 count = entries[1];
#endif

		if (count == 0)
			break;

		t = curr_row[row];

		/* Unlike the sparse blocks, medium-dense blocks
		   have enough entries that they can be stored in
		   row-major order, with many entries in each row.
		   One iteration of the while loop handles an entire
		   row at a time */

		/* curr_row and curr_b are both cached, so we have to
		   minimize the number of memory accesses and calculate
		   pointers as early as possible */

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG) && VWORDS == 1

	#define _txor(k)				\
		"movl 2*(2+2+(" #k "))(%2,%0,2), %%ecx	\n\t"	\
		"movzwl %%ax, %%edx      	\n\t"	\
		"shrl $16, %%eax         	\n\t"	\
		"movq (%1,%%edx,8), %%mm0	\n\t"	\
		"movq (%1,%%eax,8), %%mm1	\n\t"	\
		"pxor %4, %%mm0			\n\t"	\
		"pxor %4, %%mm1			\n\t"	\
		"movq %%mm0, (%1,%%edx,8)	\n\t"	\
		"movq %%mm1, (%1,%%eax,8)	\n\t"	\
		"movl 2*(2+4+(" #k "))(%2,%0,2), %%eax	\n\t"	\
		"movzwl %%cx, %%edx      	\n\t"	\
		"shrl $16, %%ecx         	\n\t"	\
		"movq (%1,%%edx,8), %%mm0	\n\t"	\
		"movq (%1,%%ecx,8), %%mm1	\n\t"	\
		"pxor %4, %%mm0			\n\t"	\
		"pxor %4, %%mm1			\n\t"	\
		"movq %%mm0, (%1,%%edx,8)	\n\t"	\
		"movq %%mm1, (%1,%%ecx,8)	\n\t"

	ASM_G volatile(
		"movl 2*(2+0)(%2,%0,2), %%eax	\n\t"
		"cmpl $0, %3			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addl $16, %0			\n\t"
		"cmpl %3, %0			\n\t"
		"jne 0b				\n\t"
		"1:				\n\t"
		:"+r"(i)
		:"r"(curr_b), "r"(entries),
		 "g"(count & (uint32)(~15)), "y"(t.w[0])
		:"%eax", "%ecx", "%edx", "%mm0", "%mm1", "memory", "cc");

	#undef _txor
	
		for (; i < count; i++)
			curr_b[entries[i+2]] = v_xor(t, curr_b[entries[i+2]]);

#elif defined(GCC_ASM64X) && VWORDS == 1

	#define _txor(k)				\
		"movzwq %%r8w, %%r9          	\n\t"	\
		"xorq %4, (%1,%%r9,8)         	\n\t"	\
		"shrq $16, %%r8              	\n\t"	\
		"xorq %4, (%1,%%r8,8)         	\n\t"	\
		"movl 2*(2+8+" #k ")(%2,%0,2), %%r8d	\n\t"	\
		"movzwq %%r10w, %%r11          	\n\t"	\
		"xorq %4, (%1,%%r11,8)         	\n\t"	\
		"shrq $16, %%r10              	\n\t"	\
		"xorq %4, (%1,%%r10,8)         	\n\t"	\
		"movl 2*(2+10+" #k ")(%2,%0,2), %%r10d	\n\t"	\
		"movzwq %%r12w, %%r13          	\n\t"	\
		"xorq %4, (%1,%%r13,8)         	\n\t"	\
		"shrq $16, %%r12              	\n\t"	\
		"xorq %4, (%1,%%r12,8)         	\n\t"	\
		"movl 2*(2+12+" #k ")(%2,%0,2), %%r12d	\n\t"	\
		"movzwq %%r14w, %%r15          	\n\t"	\
		"xorq %4, (%1,%%r15,8)         	\n\t"	\
		"shrq $16, %%r14              	\n\t"	\
		"xorq %4, (%1,%%r14,8)         	\n\t"	\
		"movl 2*(2+14+" #k ")(%2,%0,2), %%r14d	\n\t"

	ASM_G volatile(
		"movl 2*(2+0)(%2,%0,2), %%r8d	\n\t"
		"movl 2*(2+2)(%2,%0,2), %%r10d	\n\t"
		"movl 2*(2+4)(%2,%0,2), %%r12d	\n\t"
		"movl 2*(2+6)(%2,%0,2), %%r14d	\n\t"
		"cmpq $0, %3			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor(0) _txor(8)

		"addq $16, %0			\n\t"
		"cmpq %3, %0			\n\t"
		"jne 0b				\n\t"
		"1:				\n\t"
		:"+r"(i)
		:"r"(curr_b), "r"(entries),
		 "g"(count & (uint64)(~15)), "r"(t.w[0])
		:"%r8", "%r9", "%r10", "%r11", 
		 "%r12", "%r13", "%r14", "%r15", "memory", "cc");

	#undef _txor
	
		for (; i < count; i++)
			curr_b[entries[i+2]] = v_xor(t, curr_b[entries[i+2]]);

#elif defined(MSC_ASM32A) && defined(HAS_MMX) && VWORDS == 1

	#define _txor(k)				\
		ASM_M mov ecx, [2*(2+2+k)+ebx+esi*2]	\
		ASM_M movzx edx, ax			\
		ASM_M shr  eax, 16			\
		ASM_M movq mm0, [edi+edx*8]		\
		ASM_M movq mm1, [edi+eax*8]		\
		ASM_M pxor mm0, mm2			\
		ASM_M pxor mm1, mm2			\
		ASM_M movq [edi+edx*8], mm0		\
		ASM_M movq [edi+eax*8], mm1		\
		ASM_M mov eax, [2*(2+4+k)+ebx+esi*2]	\
		ASM_M movzx edx, cx			\
		ASM_M shr ecx, 16			\
		ASM_M movq mm0, [edi+edx*8]		\
		ASM_M movq mm1, [edi+ecx*8]		\
		ASM_M pxor mm0, mm2			\
		ASM_M pxor mm1, mm2			\
		ASM_M movq [edi+edx*8], mm0		\
		ASM_M movq [edi+ecx*8], mm1

	ASM_M
	{	
		push ebx
		mov esi, i
		mov edi, curr_b
		mov ebx, entries
		movq mm2, t
		mov ecx, count
		mov eax,[2*(2+0)+ebx+esi*2]
		and ecx, ~15
		je L1
		align 16
	L0:	push ecx
		_txor(0) _txor(4) _txor(8) _txor(12)
		pop ecx
		add esi, 16
		cmp esi, ecx
		jne L0
	L1: mov i, esi
		pop ebx
	}

	#undef _txor
	
		for (; i < count; i++)
			curr_b[entries[i+2]] = v_xor(t, curr_b[entries[i+2]]);
#else
		for (i = 0; i < count; i++)
			curr_b[entries[i+2]] = v_xor(t, curr_b[entries[i+2]]);
#endif
		entries += count + 2;
	}
}

/*-------------------------------------------------------------------*/
void mul_trans_one_block(const packed_block_t *curr_block,
				const v_t *curr_row, v_t * __restrict__ curr_b) {

	uint32 i = 0;
	uint32 num_entries = curr_block->num_entries;
	entry_idx_t *entries = curr_block->d.entries;

	/* unroll by 16, i.e. the number of matrix elements
	   in one cache line (usually). For 32-bit x86, we get
	   a huge performance boost by using either SSE or MMX
	   registers; not because they intrinsically are faster,
	   but because using them cuts the number of memory
	   operations in half, allowing the processor to buffer
	   more xor operations. Also convert two 16-bit reads into
	   a single 32-bit read with unpacking arithmetic */

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG) && VWORDS == 1

	#define _txor(x)				\
		"movl 4*" #x "(%1,%0,4), %%eax    \n\t"	\
		"movzwl %%ax, %%ecx               \n\t"	\
		"movq (%3,%%ecx,8), %%mm0         \n\t"	\
		"shrl $16, %%eax                  \n\t"	\
		"pxor (%2,%%eax,8), %%mm0         \n\t"	\
		"movq %%mm0, (%2,%%eax,8)         \n\t"

	ASM_G volatile(
		"cmpl $0, %4			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)

		"addl $16, %0			\n\t"
		"cmpl %4, %0			\n\t"
		"jne 0b				\n\t"
		"1:				\n\t"

		:"+r"(i)
		:"r"(entries), "r"(curr_b), "r"(curr_row), 
		 "g"(num_entries & (uint32)(~15))
		:"%eax", "%ecx", "%mm0", "memory", "cc");

	#undef _txor

	for (; i < num_entries; i++) {
		curr_b[entries[i].col_off] = v_xor(curr_b[entries[i].col_off],
						curr_row[entries[i].row_off]);
	}
#elif defined(MSC_ASM32A) && defined(HAS_MMX) && VWORDS == 1

	#define _txor(x)			\
		ASM_M mov eax,[4*x+edi+esi*4]	\
		ASM_M movzx ecx, ax		\
		ASM_M movq mm0, [edx+ecx*8]	\
		ASM_M shr eax, 16		\
		ASM_M pxor mm0, [ebx+eax*8]	\
		ASM_M movq [ebx+eax*8], mm0

	ASM_M
	{
		push ebx
		mov esi, i
		mov edi, entries
		mov ebx, curr_b
		mov ecx, num_entries
		mov edx, curr_row
		and ecx, ~15
		je L1
		align 16
	L0:	push ecx
		_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)
		pop ecx
		add esi, 16
		cmp esi, ecx
		jne L0
	L1:	mov i, esi
		pop ebx
	}

	#undef _txor

	for (; i < num_entries; i++) {
		curr_b[entries[i].col_off] = v_xor(curr_b[entries[i].col_off],
						curr_row[entries[i].row_off]);
	}
#elif defined(A64FX) && VWORDS == 8
        svbool_t __attribute__((arm_sve_vector_bits(512))) everything=svptrue_b64();
        for (i = 0; i < num_entries; i++) {
                #ifdef MANUAL_PREFETCH
                PREFETCH(entries + i + 48 / VWORDS);
                #endif
                uint64* b = (uint64*)(curr_b  + entries[i].col_off);
                const uint64* c = (const uint64*)(curr_row +entries[i].row_off);
                svstnt1_u64(everything,
                        b,
                        sveor_u64_x(everything,
                                svld1_u64(everything, b),
                                svldnt1_u64(everything, c)));
        }
#else
	for (i = 0; i < num_entries; i++) {
		#ifdef MANUAL_PREFETCH
		PREFETCH(entries + i + 48 / VWORDS);
		#endif
		curr_b[entries[i].col_off] = v_xor(curr_b[entries[i].col_off],
						curr_row[entries[i].row_off]);
	}
#endif
}
