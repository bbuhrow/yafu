/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: lanczos_vv.c 1015 2017-06-12 03:31:05Z jasonp_sf $
--------------------------------------------------------------------*/

#include "lanczos_cpu.h"

/*-------------------------------------------------------------------*/
void *vv_alloc(uint32 n, void *extra) {

	return aligned_malloc(n * sizeof(v_t), 64);
}

void vv_free(void *v) {

	aligned_free(v);
}

void vv_copyin(void *dest, v_t *src, uint32 n) {

	memcpy(dest, src, n * sizeof(v_t));
}

void vv_copy(void *dest, void *src, uint32 n) {

	memcpy(dest, src, n * sizeof(v_t));
}

void vv_copyout(v_t *dest, void *src, uint32 n) {

	memcpy(dest, src, n * sizeof(v_t));
}

void vv_clear(void *v, uint32 n) {

	memset(v, 0, n * sizeof(v_t));
}

void vv_xor(void *dest_in, void *src_in, uint32 n) {

	v_t *src = (v_t *)src_in;
	v_t *dest = (v_t *)dest_in;
	uint32 i;

	for (i = 0; i < (n & ~7); i += 8) {
		dest[i + 0] = v_xor(dest[i + 0], src[i + 0]);
		dest[i + 1] = v_xor(dest[i + 1], src[i + 1]);
		dest[i + 2] = v_xor(dest[i + 2], src[i + 2]);
		dest[i + 3] = v_xor(dest[i + 3], src[i + 3]);
		dest[i + 4] = v_xor(dest[i + 4], src[i + 4]);
		dest[i + 5] = v_xor(dest[i + 5], src[i + 5]);
		dest[i + 6] = v_xor(dest[i + 6], src[i + 6]);
		dest[i + 7] = v_xor(dest[i + 7], src[i + 7]);
	}
	for (; i < n; i++)
		dest[i] = v_xor(dest[i], src[i]);
}

void vv_mask(void *v_in, v_t mask, uint32 n) {

	v_t *v = (v_t *)v_in;
	uint32 i;

	for (i = 0; i < (n & ~7); i += 8) {
		v[i + 0] = v_and(v[i + 0], mask);
		v[i + 1] = v_and(v[i + 1], mask);
		v[i + 2] = v_and(v[i + 2], mask);
		v[i + 3] = v_and(v[i + 3], mask);
		v[i + 4] = v_and(v[i + 4], mask);
		v[i + 5] = v_and(v[i + 5], mask);
		v[i + 6] = v_and(v[i + 6], mask);
		v[i + 7] = v_and(v[i + 7], mask);
	}
	for (; i < n; i++)
		v[i] = v_and(v[i], mask);
}

/*-------------------------------------------------------------------*/
static void core_NxB_BxB_acc(v_t *v, v_t *c, v_t *y, uint32 n) {

	uint32 i;

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG) && VWORDS == 1
	i = 0;
	ASM_G volatile(
		     ALIGN_LOOP
		     "0:                                   \n\t"
		     "movq (%3,%0,8), %%mm0                \n\t"
		     "movl (%1,%0,8), %%eax                \n\t"
		     "incl %0                              \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "movq (%2,%%ecx,8), %%mm1             \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "pxor 1*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "shrl $16, %%eax                      \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "pxor 2*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "pxor 3*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movl 4-8(%1,%0,8), %%eax             \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "pxor 4*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "shrl $16, %%eax                      \n\t"
		     "cmpl %4, %0                          \n\t"
		     "pxor 5*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "pxor 6*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "pxor 7*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "pxor %%mm0, %%mm1                    \n\t"
		     "movq %%mm1, -8(%3,%0,8)              \n\t"
		     "jne 0b                               \n\t"
		     "emms                                 \n\t"
			:"+r"(i)
			:"r"(v.w), "r"(c.w), "r"(y.w), "g"(n)
			:"%eax", "%ecx", "%mm0", "%mm1", "memory");

#elif defined(MSC_ASM32A) && VWORDS == 1
	ASM_M
	{
		push    ebx
		mov	    edi,c
		mov	    esi,v
		mov     ebx,y
		xor	    ecx,ecx
		align 16
	L0:	movq	mm0,[ebx+ecx*8]
		mov	eax,[esi+ecx*8]
		inc	ecx
		movzx	edx, al
		movq	mm1,[edi+edx*8]
		movzx	edx,ah
		pxor	mm1,[1*256*8+edi+edx*8]
		shr	eax,16
		movzx	edx,al
		pxor	mm1,[2*256*8+edi+edx*8]
		movzx	edx,ah
		pxor	mm1,[3*256*8+edi+edx*8]
		mov	eax,[4-8+esi+ecx*8]
		movzx	edx,al
		pxor	mm1,[4*256*8+edi+edx*8]
		movzx	edx,ah
		shr	eax,16
		cmp	ecx,n
		pxor	mm1,[5*256*8+edi+edx*8]
		movzx	edx,al
		pxor	mm1,[6*256*8+edi+edx*8]
		movzx	edx,ah
		pxor	mm1,[7*256*8+edi+edx*8]
		pxor	mm1, mm0
		movq	[-8+ebx+ecx*8],mm1
		jne	L0
		pop	ebx
		emms
	}
#else
	for (i = 0; i < n; i++) {
		v_t vi = v[i];
		v_t accum =          c[ 0*256 + ((uint8)(vi.w[0] >>  0))];
		accum = v_xor(accum, c[ 1*256 + ((uint8)(vi.w[0] >>  8))]);
		accum = v_xor(accum, c[ 2*256 + ((uint8)(vi.w[0] >> 16))]);
		accum = v_xor(accum, c[ 3*256 + ((uint8)(vi.w[0] >> 24))]);
		accum = v_xor(accum, c[ 4*256 + ((uint8)(vi.w[0] >> 32))]);
		accum = v_xor(accum, c[ 5*256 + ((uint8)(vi.w[0] >> 40))]);
		accum = v_xor(accum, c[ 6*256 + ((uint8)(vi.w[0] >> 48))]);
		accum = v_xor(accum, c[ 7*256 + ((uint8)(vi.w[0] >> 56))]);
		#if VWORDS > 1
		accum = v_xor(accum, c[ 8*256 + ((uint8)(vi.w[1] >>  0))]);
		accum = v_xor(accum, c[ 9*256 + ((uint8)(vi.w[1] >>  8))]);
		accum = v_xor(accum, c[10*256 + ((uint8)(vi.w[1] >> 16))]);
		accum = v_xor(accum, c[11*256 + ((uint8)(vi.w[1] >> 24))]);
		accum = v_xor(accum, c[12*256 + ((uint8)(vi.w[1] >> 32))]);
		accum = v_xor(accum, c[13*256 + ((uint8)(vi.w[1] >> 40))]);
		accum = v_xor(accum, c[14*256 + ((uint8)(vi.w[1] >> 48))]);
		accum = v_xor(accum, c[15*256 + ((uint8)(vi.w[1] >> 56))]);
		#if VWORDS > 2
		accum = v_xor(accum, c[16*256 + ((uint8)(vi.w[2] >>  0))]);
		accum = v_xor(accum, c[17*256 + ((uint8)(vi.w[2] >>  8))]);
		accum = v_xor(accum, c[18*256 + ((uint8)(vi.w[2] >> 16))]);
		accum = v_xor(accum, c[19*256 + ((uint8)(vi.w[2] >> 24))]);
		accum = v_xor(accum, c[20*256 + ((uint8)(vi.w[2] >> 32))]);
		accum = v_xor(accum, c[21*256 + ((uint8)(vi.w[2] >> 40))]);
		accum = v_xor(accum, c[22*256 + ((uint8)(vi.w[2] >> 48))]);
		accum = v_xor(accum, c[23*256 + ((uint8)(vi.w[2] >> 56))]);
		#if VWORDS > 3
		accum = v_xor(accum, c[24*256 + ((uint8)(vi.w[3] >>  0))]);
		accum = v_xor(accum, c[25*256 + ((uint8)(vi.w[3] >>  8))]);
		accum = v_xor(accum, c[26*256 + ((uint8)(vi.w[3] >> 16))]);
		accum = v_xor(accum, c[27*256 + ((uint8)(vi.w[3] >> 24))]);
		accum = v_xor(accum, c[28*256 + ((uint8)(vi.w[3] >> 32))]);
		accum = v_xor(accum, c[29*256 + ((uint8)(vi.w[3] >> 40))]);
		accum = v_xor(accum, c[30*256 + ((uint8)(vi.w[3] >> 48))]);
		accum = v_xor(accum, c[31*256 + ((uint8)(vi.w[3] >> 56))]);
		#endif
		#endif
		#endif
		y[i] = v_xor(y[i], accum);
	}
#endif
}

/*-------------------------------------------------------------------*/
static const uint8 graycode[2 * 256] = {
   0, 0,    1, 0,    3, 1,    2, 0,    6, 2,    7, 0,    5, 1,    4, 0,   
  12, 3,   13, 0,   15, 1,   14, 0,   10, 2,   11, 0,    9, 1,    8, 0,   
  24, 4,   25, 0,   27, 1,   26, 0,   30, 2,   31, 0,   29, 1,   28, 0,   
  20, 3,   21, 0,   23, 1,   22, 0,   18, 2,   19, 0,   17, 1,   16, 0,   
  48, 5,   49, 0,   51, 1,   50, 0,   54, 2,   55, 0,   53, 1,   52, 0,  
  60, 3,   61, 0,   63, 1,   62, 0,   58, 2,   59, 0,   57, 1,   56, 0,   
  40, 4,   41, 0,   43, 1,   42, 0,   46, 2,   47, 0,   45, 1,   44, 0,  
  36, 3,   37, 0,   39, 1,   38, 0,   34, 2,   35, 0,   33, 1,   32, 0,   
  96, 6,   97, 0,   99, 1,   98, 0,  102, 2,  103, 0,  101, 1,  100, 0,  
 108, 3,  109, 0,  111, 1,  110, 0,  106, 2,  107, 0,  105, 1,  104, 0,  
 120, 4,  121, 0,  123, 1,  122, 0,  126, 2,  127, 0,  125, 1,  124, 0,  
 116, 3,  117, 0,  119, 1,  118, 0,  114, 2,  115, 0,  113, 1,  112, 0,  
  80, 5,   81, 0,   83, 1,   82, 0,   86, 2,   87, 0,   85, 1,   84, 0,   
  92, 3,   93, 0,   95, 1,   94, 0,   90, 2,   91, 0,   89, 1,   88, 0,   
  72, 4,   73, 0,   75, 1,   74, 0,   78, 2,   79, 0,   77, 1,   76, 0,   
  68, 3,   69, 0,   71, 1,   70, 0,   66, 2,   67, 0,   65, 1,   64, 0,  
 192, 7,  193, 0,  195, 1,  194, 0,  198, 2,  199, 0,  197, 1,  196, 0, 
 204, 3,  205, 0,  207, 1,  206, 0,  202, 2,  203, 0,  201, 1,  200, 0, 
 216, 4,  217, 0,  219, 1,  218, 0,  222, 2,  223, 0,  221, 1,  220, 0, 
 212, 3,  213, 0,  215, 1,  214, 0,  210, 2,  211, 0,  209, 1,  208, 0, 
 240, 5,  241, 0,  243, 1,  242, 0,  246, 2,  247, 0,  245, 1,  244, 0, 
 252, 3,  253, 0,  255, 1,  254, 0,  250, 2,  251, 0,  249, 1,  248, 0, 
 232, 4,  233, 0,  235, 1,  234, 0,  238, 2,  239, 0,  237, 1,  236, 0, 
 228, 3,  229, 0,  231, 1,  230, 0,  226, 2,  227, 0,  225, 1,  224, 0,  
 160, 6,  161, 0,  163, 1,  162, 0,  166, 2,  167, 0,  165, 1,  164, 0, 
 172, 3,  173, 0,  175, 1,  174, 0,  170, 2,  171, 0,  169, 1,  168, 0, 
 184, 4,  185, 0,  187, 1,  186, 0,  190, 2,  191, 0,  189, 1,  188, 0, 
 180, 3,  181, 0,  183, 1,  182, 0,  178, 2,  179, 0,  177, 1,  176, 0, 
 144, 5,  145, 0,  147, 1,  146, 0,  150, 2,  151, 0,  149, 1,  148, 0, 
 156, 3,  157, 0,  159, 1,  158, 0,  154, 2,  155, 0,  153, 1,  152, 0, 
 136, 4,  137, 0,  139, 1,  138, 0,  142, 2,  143, 0,  141, 1,  140, 0, 
 132, 3,  133, 0,  135, 1,  134, 0,  130, 2,  131, 0,  129, 1,  128, 0,  
};

static void mul_NxB_BxB_precomp(v_t *c, v_t *x) {

	/* fill c[][] with a bunch of "partial matrix multiplies". 
	   For 0<=j<256 and 0<=i<8*VWORDS, the i_th row of c[][] 
	   contains the matrix product

	   	( j << (8*i) ) * x[][]

	   where the quantity in parentheses is considered a 
	   1 x VBITS vector of elements in GF(2). The resulting
	   table will make matrix multiplies by x[][] about 8x faster
	 
	   We iterate through i in Gray code order to minimize overhead */

	uint32 i;
	v_t acc[8 * VWORDS];

	for (i = 0; i < 8 * VWORDS; i++)
		acc[i] = c[i * 256] = v_zero;

#define BXB_ACC(i) \
	acc[i] = v_xor(acc[i], x[i*8 + bit]); c[i*256 + word] = acc[i]

	for (i = 1; i < 256; i++) {

		uint32 word = graycode[2 * i];
		uint32 bit = graycode[2 * i + 1];

		BXB_ACC(0); BXB_ACC(1); BXB_ACC(2); BXB_ACC(3);
		BXB_ACC(4); BXB_ACC(5); BXB_ACC(6); BXB_ACC(7);
		#if VWORDS > 1
		BXB_ACC(8); BXB_ACC(9); BXB_ACC(10); BXB_ACC(11);
		BXB_ACC(12); BXB_ACC(13); BXB_ACC(14); BXB_ACC(15);
		#if VWORDS > 2
		BXB_ACC(16); BXB_ACC(17); BXB_ACC(18); BXB_ACC(19);
		BXB_ACC(20); BXB_ACC(21); BXB_ACC(22); BXB_ACC(23);
		#if VWORDS > 3
		BXB_ACC(24); BXB_ACC(25); BXB_ACC(26); BXB_ACC(27);
		BXB_ACC(28); BXB_ACC(29); BXB_ACC(30); BXB_ACC(31);
		#endif
		#endif
		#endif
	}
}

/*-------------------------------------------------------------------*/
void mul_NxB_BxB_acc(v_t *v, v_t *x, v_t *y, uint32 n) {

	/* let v[][] be a N x B matrix with elements in GF(2), 
	   represented as an array of n v_t structures. Let c[][]
	   be an (8*VWORDS) x 256 scratch matrix of v_t structures.
	   This code multiplies v[][] by the BxB matrix 
	   x[][], then XORs the n x 64 result into y[][] */

	v_t c[8 * VWORDS * 256];

	mul_NxB_BxB_precomp(c, x);

	core_NxB_BxB_acc(v, c, y, n);
}

/*-------------------------------------------------------------------*/
static void outer_thread_run(void *data, int thread_num)
{
	la_task_t *task = (la_task_t *)data;
	packed_matrix_t *p = task->matrix;
	cpudata_t *cpudata = (cpudata_t *)p->extra;
	thread_data_t *t = cpudata->thread_data + task->task_num;

	core_NxB_BxB_acc(t->x, t->b, t->y, t->vsize);
}

void vv_mul_NxB_BxB_acc(packed_matrix_t *matrix, 
			void *v_in, v_t *x,
			void *y_in, uint32 n) {

	cpudata_t *cpudata = (cpudata_t *)matrix->extra;
	v_t *v = (v_t *)v_in;
	v_t *y = (v_t *)y_in;
	v_t c[8 * VWORDS * 256];
	uint32 i;
	uint32 vsize = n / matrix->num_threads;
	uint32 off;
	task_control_t task = {NULL, NULL, NULL, NULL};

	mul_NxB_BxB_precomp(c, x);

	for (i = off = 0; i < matrix->num_threads; i++, off += vsize) {

		thread_data_t *t = cpudata->thread_data + i;

		t->x = v + off;
		t->b = c;
		t->y = y + off;
		if (i == matrix->num_threads - 1)
			t->vsize = n - off;
		else
			t->vsize = vsize;
	}

	task.run = outer_thread_run;

	for (i = 0; i < matrix->num_threads - 1; i++) {
		task.data = cpudata->tasks + i;
		threadpool_add_task(cpudata->threadpool, &task, 0);
	}
	outer_thread_run(cpudata->tasks + i, i);

	if (i > 0)
		threadpool_drain(cpudata->threadpool, 1);
}

/*-------------------------------------------------------------------*/
static void core_BxN_NxB(v_t *x, v_t *c, v_t *y, uint32 n) {

	uint32 i;

	memset(c, 0, 8 * VWORDS * 256 * sizeof(v_t));

#if defined(GCC_ASM32A) && defined(HAS_MMX) && defined(NDEBUG) && VWORDS == 1
	i = 0;
	ASM_G volatile(
		     ALIGN_LOOP
		     "0:                                   \n\t"
		     "movq (%3,%0,8), %%mm0                \n\t"
		     "movl (%1,%0,8), %%eax                \n\t"
		     "incl %0                              \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "movq %%mm0, %%mm1                    \n\t"
		     "pxor (%2,%%ecx,8), %%mm1             \n\t"
		     "movq %%mm1, (%2,%%ecx,8)             \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "movq %%mm0, %%mm1                    \n\t"
		     "pxor 1*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movq %%mm1, 1*256*8(%2,%%ecx,8)      \n\t"
		     "shrl $16, %%eax                      \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "movq %%mm0, %%mm1                    \n\t"
		     "pxor 2*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movq %%mm1, 2*256*8(%2,%%ecx,8)      \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "movq %%mm0, %%mm1                    \n\t"
		     "pxor 3*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movq %%mm1, 3*256*8(%2,%%ecx,8)      \n\t"
		     "movl 4-8(%1,%0,8), %%eax             \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "movq %%mm0, %%mm1                    \n\t"
		     "pxor 4*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movq %%mm1, 4*256*8(%2,%%ecx,8)      \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "shrl $16, %%eax                      \n\t"
		     "cmpl %4, %0                          \n\t"
		     "movq %%mm0, %%mm1                    \n\t"
		     "pxor 5*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movq %%mm1, 5*256*8(%2,%%ecx,8)      \n\t"
		     "movzbl %%al, %%ecx                   \n\t"
		     "movq %%mm0, %%mm1                    \n\t"
		     "pxor 6*256*8(%2,%%ecx,8), %%mm1      \n\t"
		     "movq %%mm1, 6*256*8(%2,%%ecx,8)      \n\t"
		     "movzbl %%ah, %%ecx                   \n\t"
		     "pxor 7*256*8(%2,%%ecx,8), %%mm0      \n\t"
		     "movq %%mm0, 7*256*8(%2,%%ecx,8)      \n\t"
		     "jne 0b                               \n\t"
		     "emms                                 \n\t"
			:"+r"(i)
			:"r"(x), "r"(c), "r"(y), "g"(n)
			:"%eax", "%ecx", "%mm0", "%mm1", "memory");

#elif defined(MSC_ASM32A) && VWORDS == 1
	ASM_M
	{
		push    ebx
		mov	    edi,c
		mov	    esi,x
		mov     ebx,y
		xor	    ecx,ecx
		align 16
    L0:	movq	mm0,[ebx+ecx*8]
		mov	    eax,[esi+ecx*8]
		inc	    ecx
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[edi+edx*8]
		movq	[edi+edx*8],mm1
		movzx	edx,ah
		movq	mm1, mm0
		pxor	mm1,[1*256*8+edi+edx*8]
		movq	[1*256*8+edi+edx*8],mm1
		shr	    eax,16
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[2*256*8+edi+edx*8]
		movq	[2*256*8+edi+edx*8],mm1
		movzx	edx,ah
		movq	mm1,mm0
		pxor	mm1,[3*256*8+edi+edx*8]
		movq	[3*256*8+edi+edx*8],mm1
		mov	    eax,[4-8+esi+ecx*8]
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[4*256*8+edi+edx*8]
		movq	[4*256*8+edi+edx*8],mm1
		movzx	edx,ah
		shr	    eax,16
		cmp	    ecx,n
		movq	mm1,mm0
		pxor	mm1,[5*256*8+edi+edx*8]
		movq	[5*256*8+edi+edx*8],mm1
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[6*256*8+edi+edx*8]
		movq	[6*256*8+edi+edx*8],mm1
		movzx	edx,ah
		pxor	mm0,[7*256*8+edi+edx*8]
		movq	[7*256*8+edi+edx*8],mm0
		jne	    L0
		emms 
		pop     ebx
	}
#else

	#define NXB_ACC(i) \
		c[i*256 + (uint8)(xi.w[i/8] >> (8*(i % 8)))] = \
		v_xor(c[i*256 + (uint8)(xi.w[i/8] >> (8*(i % 8)))], yi)

	for (i = 0; i < n; i++) {
		v_t xi = x[i];
		v_t yi = y[i];

		NXB_ACC(0); NXB_ACC(1); NXB_ACC(2); NXB_ACC(3);
		NXB_ACC(4); NXB_ACC(5); NXB_ACC(6); NXB_ACC(7);
		#if VWORDS > 1
		NXB_ACC(8); NXB_ACC(9); NXB_ACC(10); NXB_ACC(11);
		NXB_ACC(12); NXB_ACC(13); NXB_ACC(14); NXB_ACC(15);
		#if VWORDS > 2
		NXB_ACC(16); NXB_ACC(17); NXB_ACC(18); NXB_ACC(19);
		NXB_ACC(20); NXB_ACC(21); NXB_ACC(22); NXB_ACC(23);
		#if VWORDS > 3
		NXB_ACC(24); NXB_ACC(25); NXB_ACC(26); NXB_ACC(27);
		NXB_ACC(28); NXB_ACC(29); NXB_ACC(30); NXB_ACC(31);
		#endif
		#endif
		#endif
	}
#endif
}

/*-------------------------------------------------------------------*/
static void mul_BxN_NxB_postproc(v_t *c, v_t *xy) {

	uint32 i, j;

	#define NXB_POST(i) \
		a[i] = v_xor(a[i], c[i*256 + j])

	for (i = 0; i < 8; i++) {

		v_t a[8 * VWORDS];

		for (j = 0; j < 8 * VWORDS; j++)
			a[j] = v_zero;

		for (j = 0; j < 256; j++) {
			if ((j >> i) & 1) {
				NXB_POST(0); NXB_POST(1); NXB_POST(2); NXB_POST(3);
				NXB_POST(4); NXB_POST(5); NXB_POST(6); NXB_POST(7);
				#if VWORDS > 1
				NXB_POST(8); NXB_POST(9); NXB_POST(10); NXB_POST(11);
				NXB_POST(12); NXB_POST(13); NXB_POST(14); NXB_POST(15);
				#if VWORDS > 2
				NXB_POST(16); NXB_POST(17); NXB_POST(18); NXB_POST(19);
				NXB_POST(20); NXB_POST(21); NXB_POST(22); NXB_POST(23);
				#if VWORDS > 3
				NXB_POST(24); NXB_POST(25); NXB_POST(26); NXB_POST(27);
				NXB_POST(28); NXB_POST(29); NXB_POST(30); NXB_POST(31);
				#endif
				#endif
				#endif
			}
		}

		for (j = 0; j < 8 * VWORDS; j++)
			xy[8 * j] = a[j];
		xy++;
	}
}

/*-------------------------------------------------------------------*/
void mul_BxN_NxB(v_t *x, v_t *y, v_t *xy, uint32 n) {

	/* Let x and y be N x B matrices. This routine computes
	   the B x B matrix xy[][] given by transpose(x) * y */

	v_t c[8 * VWORDS * 256];

	core_BxN_NxB(x, c, y, n);

	mul_BxN_NxB_postproc(c, xy);
}

/*-------------------------------------------------------------------*/
static void inner_thread_run(void *data, int thread_num)
{
	la_task_t *task = (la_task_t *)data;
	packed_matrix_t *p = task->matrix;
	cpudata_t *cpudata = (cpudata_t *)p->extra;
	thread_data_t *t = cpudata->thread_data + task->task_num;

	mul_BxN_NxB(t->x, t->y, t->tmp_b, t->vsize);
}

void vv_mul_BxN_NxB(packed_matrix_t *matrix,
		   void *x_in, void *y_in,
		   v_t *xy, uint32 n) {


	cpudata_t *cpudata = (cpudata_t *)matrix->extra;
	v_t *x = (v_t *)x_in;
	v_t *y = (v_t *)y_in;
	uint32 i;
	uint32 vsize = n / matrix->num_threads;
	uint32 off;
	task_control_t task = {NULL, NULL, NULL, NULL};
#ifdef HAVE_MPI
	v_t xytmp[VBITS];
#endif

	for (i = off = 0; i < matrix->num_threads; i++, off += vsize) {
		thread_data_t *t = cpudata->thread_data + i;

		t->x = x + off;
		t->y = y + off;

		if (i == matrix->num_threads - 1)
			t->vsize = n - off;
		else
			t->vsize = vsize;
	}

	task.run = inner_thread_run;

	for (i = 0; i < matrix->num_threads - 1; i++) {
		task.data = cpudata->tasks + i;
		threadpool_add_task(cpudata->threadpool, &task, 0);
	}
	inner_thread_run(cpudata->tasks + i, i);

	/* All the scratch vectors used by threads get 
	   xor-ed into the final xy vector */

	vv_copy(xy, cpudata->thread_data[i].tmp_b, VBITS);

	if (i > 0) {
		threadpool_drain(cpudata->threadpool, 1);

		for (i = 0; i < matrix->num_threads - 1; i++) {
			thread_data_t *t = cpudata->thread_data + i;

			vv_xor(xy, t->tmp_b, VBITS);
		}
	}

#ifdef HAVE_MPI
	/* combine the results across an entire MPI row */

	global_xor(xy, xytmp, VBITS, matrix->mpi_ncols,
			matrix->mpi_la_col_rank,
			matrix->mpi_la_row_grid);

	/* combine the results across an entire MPI column */
    
	global_xor(xytmp, xy, VBITS, matrix->mpi_nrows,
			matrix->mpi_la_row_rank,
			matrix->mpi_la_col_grid);    
#endif
}
