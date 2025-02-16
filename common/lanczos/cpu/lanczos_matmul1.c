/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: lanczos_matmul1.c 1012 2017-06-11 17:42:14Z jasonp_sf $
--------------------------------------------------------------------*/

#include "lanczos_cpu.h"

	/* code for handling matrix multiplies when the
	   matrix is in packed format */

/*-------------------------------------------------------------------*/

static void mul_one_med_block(packed_block_t *curr_block,
			v_t *curr_col, v_t *curr_b) {

	uint16 *entries = curr_block->d.med_entries;

	while (1) {
		v_t accum;

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

		/* Unlike the sparse blocks, medium-dense blocks
		   have enough entries that they can be stored in
		   row-major order, with many entries in each row.
		   One iteration of the while loop handles an entire
		   row at a time */

		/* curr_col and curr_b are both cached, so we have to
		   minimize the number of memory accesses and calculate
		   pointers as early as possible */

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG) && VWORDS == 1

	#define _txor(k)				\
		"movzwl %%ax, %%edx		\n\t"	\
		"pxor (%2,%%edx,8), %0		\n\t"	\
		"shrl $16, %%eax		\n\t"	\
		"pxor (%2,%%eax,8), %%mm0	\n\t"	\
		"movl 2*(2+4+(" #k "))(%3,%1,2), %%eax \n\t"	\
		"movzwl %%cx, %%edx		\n\t"	\
		"pxor (%2,%%edx,8), %0		\n\t"	\
		"shrl $16, %%ecx		\n\t"	\
		"pxor (%2,%%ecx,8), %%mm0	\n\t"	\
		"movl 2*(2+6+(" #k "))(%3,%1,2), %%ecx \n\t"

	ASM_G volatile(
		"movl 2*(2+0)(%3,%1,2), %%eax	\n\t"
		"movl 2*(2+2)(%3,%1,2), %%ecx	\n\t"
		"pxor %0, %0			\n\t"
		"pxor %%mm0, %%mm0		\n\t"
		"cmpl $0, %4			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addl $16, %1			\n\t"
		"cmpl %4, %1			\n\t"
		"jne 0b				\n\t"
		"pxor %%mm0, %0			\n\t"
		"1:				\n\t"

		:"=y"(accum.w[0]), "+r"(i)
		:"r"(curr_col), "r"(entries),
		 "g"(count & (uint32)(~15))
		:"%eax", "%ecx", "%edx", "%mm0", "memory", "cc");

	#undef _txor

#elif defined(GCC_ASM64X) && VWORDS == 1

	#define _txor(k)				\
		"movzwq %%ax, %%rdx		\n\t"	\
		"xorq (%2,%%rdx,8), %0		\n\t"	\
		"shrq $16, %%rax		\n\t"	\
		"xorq (%2,%%rax,8), %%rsi	\n\t"	\
		"movl 2*(2+4+(" #k "))(%3,%1,2), %%eax \n\t"	\
		"movzwq %%cx, %%rdx		\n\t"	\
		"xorq (%2,%%rdx,8), %0		\n\t"	\
		"shrq $16, %%rcx		\n\t"	\
		"xorq (%2,%%rcx,8), %%rsi	\n\t"	\
		"movl 2*(2+6+(" #k "))(%3,%1,2), %%ecx \n\t"

	ASM_G volatile(
		"movl 2*(2+0)(%3,%1,2), %%eax	\n\t"
		"movl 2*(2+2)(%3,%1,2), %%ecx	\n\t"
		"xorq %0, %0			\n\t"
		"xorq %%rsi, %%rsi		\n\t"
		"cmpq $0, %4			\n\t"
		"je 1f				\n\t"
		ALIGN_LOOP
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addq $16, %1			\n\t"
		"cmpq %4, %1			\n\t"
		"jne 0b				\n\t"
		"xorq %%rsi, %0			\n\t"
		"1:				\n\t"

		:"=&r"(accum.w[0]), "+r"(i)
		:"r"(curr_col), "r"(entries), 
		 "g"(count & (uint64)(~15))
		:"%rax", "%rcx", "%rdx", "%rsi", "memory", "cc");

	#undef _txor

#elif defined(MSC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG) && VWORDS == 1

	#define _txor(k)				\
	    ASM_M movzx edx, ax				\
	    ASM_M pxor	mm1, [esi+edx*8]		\
	    ASM_M shr 	eax, 16				\
	    ASM_M pxor 	mm0, [esi+eax*8]		\
	    ASM_M mov 	eax, [2*(k+2+4)+ebx+edi*2]	\
	    ASM_M movzx edx, cx				\
	    ASM_M pxor 	mm1, [esi+edx*8]		\
	    ASM_M shr	ecx, 16				\
	    ASM_M pxor 	mm0, [esi+ecx*8]		\
	    ASM_M mov 	ecx, [2*(k+2+6)+ebx+edi*2]
			
	ASM_M
	{	
		push ebx
		mov esi, curr_col
		mov ebx, entries
		mov edi, i
		mov edx, count
		mov eax, [2*(2+0)+ebx+edi*2]
		mov ecx, [2*(2+2)+ebx+edi*2]
		pxor mm1, mm1
		pxor mm0, mm0
		and edx, ~15
		je L1
		align 16
	L0:	push edx
		_txor(0) _txor(4) _txor(8) _txor(12)
		pop edx
		add edi, 16
		cmp edi, edx
		jne L0
		pxor mm1, mm0
	L1:	movq accum, mm1
		mov i, edi
		pop ebx
	}

	#undef _txor

#else
	accum = v_zero;
	for (i = 0; i < (count & (uint32)(~15)); i += 16) {
		accum = v_xor(accum, curr_col[entries[i+2+0]]);
		accum = v_xor(accum, curr_col[entries[i+2+1]]);
		accum = v_xor(accum, curr_col[entries[i+2+2]]);
		accum = v_xor(accum, curr_col[entries[i+2+3]]);
		accum = v_xor(accum, curr_col[entries[i+2+4]]);
		accum = v_xor(accum, curr_col[entries[i+2+5]]);
		accum = v_xor(accum, curr_col[entries[i+2+6]]);
		accum = v_xor(accum, curr_col[entries[i+2+7]]);
		accum = v_xor(accum, curr_col[entries[i+2+8]]);
		accum = v_xor(accum, curr_col[entries[i+2+9]]);
		accum = v_xor(accum, curr_col[entries[i+2+10]]);
		accum = v_xor(accum, curr_col[entries[i+2+11]]);
		accum = v_xor(accum, curr_col[entries[i+2+12]]);
		accum = v_xor(accum, curr_col[entries[i+2+13]]);
		accum = v_xor(accum, curr_col[entries[i+2+14]]);
		accum = v_xor(accum, curr_col[entries[i+2+15]]);
	}

#endif
		for (; i < count; i++)
			accum = v_xor(accum, curr_col[entries[i+2]]);
		curr_b[row] = v_xor(curr_b[row], accum);
		entries += count + 2;
	}
}

/*-------------------------------------------------------------------*/
static void mul_one_block(packed_block_t *curr_block,
			v_t *curr_col, v_t *curr_b) {

	uint32 i = 0; 
	uint32 num_entries = curr_block->num_entries;
	entry_idx_t *entries = curr_block->d.entries;

	/* unroll by 16, i.e. the number of matrix elements
	   in one cache line (usually). For 32-bit x86, we get
	   a huge performance boost by using either SSE or MMX
	   registers; not because they intrinsically are faster,
	   but because using them cuts the number of memory
	   operations in half, allowing the processor to buffer
	   more xor operations. Also replace two 16-bit loads
	   with a single 32-bit load and extra arithmetic to
	   unpack the array indices */

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG) && VWORDS == 1

	#define _txor(x)				\
		"movl 4*" #x "(%1,%0,4), %%eax  \n\t"	\
		"movzwl %%ax, %%ecx             \n\t"	\
		"movq (%2,%%ecx,8), %%mm0       \n\t"	\
		"shrl $16, %%eax                \n\t"	\
		"pxor (%3,%%eax,8), %%mm0       \n\t"	\
		"movq %%mm0, (%2,%%ecx,8)       \n\t"

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
		:"r"(entries), "r"(curr_b), "r"(curr_col), 
		 "g"(num_entries & (uint32)(~15))
		:"%eax", "%ecx", "%mm0", "memory", "cc");

#elif defined(MSC_ASM32A) && defined(HAS_MMX) && VWORDS == 1

	#define _txor(x)				\
		ASM_M mov	eax, [4*x+edi+esi*4]	\
		ASM_M movzx ecx, ax			\
		ASM_M movq 	mm0, [ebx+ecx*8]	\
		ASM_M shr 	eax, 16			\
		ASM_M pxor 	mm0, [edx+eax*8]	\
		ASM_M movq 	[ebx+ecx*8], mm0
	
	ASM_M
	{
		push ebx
		mov esi, i
		mov edi, entries
		mov ebx, curr_b
		mov ecx, num_entries
		mov edx, curr_col
		and ecx, ~15
		je L1
		align 16
	L0:	push ecx
		_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)
		pop ecx
		add esi,16
		cmp esi,ecx
		jne L0
	L1:	mov i, esi
		pop ebx
	}

#else
	#define _txor(x) curr_b[entries[i+x].row_off] = v_xor( \
			      	curr_b[entries[i+x].row_off], \
				curr_col[entries[i+x].col_off])

	for (i = 0; i < (num_entries & (uint32)(~15)); i += 16) {
		#ifdef MANUAL_PREFETCH
		PREFETCH(entries + i + 48 / VWORDS);
		#endif

		_txor( 0); _txor( 1); _txor( 2); _txor( 3);
		_txor( 4); _txor( 5); _txor( 6); _txor( 7);
		_txor( 8); _txor( 9); _txor(10); _txor(11);
		_txor(12); _txor(13); _txor(14); _txor(15);
	}
#endif
	#undef _txor

	for (; i < num_entries; i++) {
		curr_b[entries[i].row_off] = v_xor(curr_b[entries[i].row_off],
						curr_col[entries[i].col_off]);
	}
}

/*-------------------------------------------------------------------*/
void mul_packed_core(void *data, int thread_num)
{
	/* we skip the first matrix row, since it is handled 
	   in the dense function below */

	la_task_t *task = (la_task_t *)data;
	packed_matrix_t *p = task->matrix;
	cpudata_t *c = (cpudata_t *)p->extra;

	uint32 start_block_c = task->block_num * c->superblock_size;
	uint32 num_blocks_c = MIN(c->superblock_size, 
				c->num_block_cols - start_block_c);

	packed_block_t *start_block = c->blocks + start_block_c +
					c->num_block_cols;
	v_t *x = c->x + start_block_c * c->block_size;
	uint32 i, j;

	for (i = task->task_num; i < c->num_block_rows - 1; 
					i += p->num_threads) {

		packed_block_t *curr_block = start_block + 
					i * c->num_block_cols;
		v_t *curr_x = x;
		uint32 b_off = i * c->block_size + c->first_block_size;
		v_t *b = c->b + b_off;

		if (start_block_c == 0)
			vv_clear(b, MIN(c->block_size, p->nrows - b_off));

		for (j = 0; j < num_blocks_c; j++) {
			mul_one_block(curr_block, curr_x, b);
			curr_block++;
			curr_x += c->block_size;
		}
	}
}

/*-------------------------------------------------------------------*/
void mul_packed_small_core(void *data, int thread_num)
{
	la_task_t *task = (la_task_t *)data;
	packed_matrix_t *p = task->matrix;
	cpudata_t *c = (cpudata_t *)p->extra;
	thread_data_t *t = c->thread_data + task->task_num;

	uint32 last_task = (task->task_num == p->num_threads - 1);
	uint32 num_blocks = c->num_block_cols / p->num_threads;
	uint32 block_off = num_blocks * task->task_num;
	uint32 off = c->block_size * block_off;
	uint32 vsize = num_blocks * c->block_size;
	v_t *x = c->x + off;
	v_t *b = t->tmp_b;
	packed_block_t *curr_block = c->blocks + block_off;
	uint32 i;

	vv_clear(b, MAX(c->first_block_size, VBITS * 
			(1 + (p->num_dense_rows + VBITS - 1) / VBITS)));

	if (p->num_threads == 1) {
		vsize = p->ncols;
	}
	else if (last_task) {
		num_blocks = c->num_block_cols - block_off;
		vsize = p->ncols - off;
	}

	for (i = 0; i < num_blocks; i++) {
		mul_one_med_block(curr_block, x, b);
		curr_block++;
		x += c->block_size;
	}

	/* multiply the densest few rows by x (in batches of VBITS rows) */

	for (i = 0; i < (p->num_dense_rows + VBITS - 1) / VBITS; i++)
		mul_BxN_NxB(c->dense_blocks[i] + off, 
				c->x + off, b + VBITS * i, vsize);
}
