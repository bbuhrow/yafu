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

void vv_xor(void * __restrict__ dest_in, void * __restrict__ src_in, uint32 n) {

	v_t *src = (v_t *)src_in;
	v_t *dest = (v_t *)dest_in;
	uint32 i;

#pragma omp parallel for
	for (i = 0; i < n; i++)
		dest[i] = v_xor(dest[i], src[i]);
}

void vv_mask(void *v_in, v_t mask, uint32 n) {

	v_t *v = (v_t *)v_in;
	uint32 i;

#pragma omp parallel for
	for (i = 0; i < n; i++)
		v[i] = v_and(v[i], mask);
}

/*-------------------------------------------------------------------*/
static void core_NxB_BxB_acc(const v_t *v, const v_t *c, v_t * __restrict__ y, uint32 n) {

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
#pragma omp parallel for
	for (i = 0; i < n; i++) {
		#ifdef MANUAL_PREFETCH
		PREFETCH(y+i+4);
		#endif
		uint32 j;
		v_t vi = v[i];
		v_t accum;
		for (j = 0; j < VWORDS; j++) accum.w[j] = 0;

		for (j = 0; j < 8 * VWORDS; j++) {
			uint32 k = j*256 + ((vi.w[(j >> 3)] >> (8*(j & 7))) & 255);
			accum = v_xor(accum, c[k]);
		}
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

	uint32 i, j;
	v_t acc[8 * VWORDS];

#pragma omp parallel for
	for (i = 0; i < 8 * VWORDS; i++)
		acc[i] = c[i * 256] = v_zero;

#define BXB_ACC(i) \
	acc[i] = v_xor(acc[i], x[i*8 + bit]); c[i*256 + word] = acc[i]

	for (i = 1; i < 256; i++) {

		uint32 word = graycode[2 * i];
		uint32 bit = graycode[2 * i + 1];

		for (j = 0; j < 8 * VWORDS; j++) { BXB_ACC(j); }
	}
}

/*-------------------------------------------------------------------*/
void mul_NxB_BxB_acc(v_t *v, v_t *x, v_t *y, uint32 n) {

	/* let v[][] be a N x B matrix with elements in GF(2), 
	   represented as an array of n v_t structures. Let c[][]
	   be an (8*VWORDS) x 256 scratch matrix of v_t structures.
	   This code multiplies v[][] by the BxB matrix 
	   x[][], then XORs the N x B result into y[][] */

	v_t c[8 * VWORDS * 256];

	mul_NxB_BxB_precomp(c, x);

	core_NxB_BxB_acc(v, c, y, n);
}

/*-------------------------------------------------------------------*/
void vv_mul_NxB_BxB_acc(packed_matrix_t *matrix, 
			void *v_in, v_t *x,
			void *y_in, uint32 n) {

	v_t *v = (v_t *)v_in;
	v_t *y = (v_t *)y_in;
	
	mul_NxB_BxB_acc(v, x, y, n);
}

/*-------------------------------------------------------------------*/
static void core_BxN_NxB(const v_t *x, v_t *c, const v_t *y, const uint32 n) {

	uint32 i, j;

	memset(c, 0, 8 * VWORDS * 256 * sizeof(v_t));
// #pragma omp parallel for
//	for (i = 0; i < 8 * VWORDS * 256; i++)
//		for (j = 0; j < VWORDS; j++) c[i].w[j] = 0;

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
		k = i*256 + ((xi.w[(i >> 3)] >> (8*(i & 7))) & 255); c[k] = v_xor(c[k], yi)

// #pragma omp parallel for private(j) reduction(^:c[0:8 * VWORDS * 256])
	for (i = 0; i < n; i++) {
		v_t xi = x[i];
		v_t yi = y[i];

		for (j = 0; j < 8 * VWORDS; j++) { 
			uint32 k;
			NXB_ACC(j); 
		}
	}
#endif
}

/*-------------------------------------------------------------------*/
static void mul_BxN_NxB_postproc(v_t *c, v_t *xy) {

	uint32 i, j, k;

	#define NXB_POST(i) \
		a[i] = v_xor(a[i], c[i*256 + j])

	for (i = 0; i < 8; i++) {

		v_t a[8 * VWORDS];

// #pragma omp parallel for
		for (j = 0; j < 8 * VWORDS; j++)
			a[j] = v_zero;

// #pragma omp parallel for private(k) reduction(^:a)
		for (j = 0; j < 256; j++) {
			if ((j >> i) & 1) {
				for (k = 0; k < 8 * VWORDS; k++) { NXB_POST(k); }
			}
		}

// #pragma omp parallel for
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
void mul_BxN_NxB_2(v_t *x, v_t *y, v_t *xy, uint32 n) {

	/* Let x and y be N x B matrices. This routine computes
	   the B x B matrix xy[][] given by transpose(x) * y 
	   This version is as fast as or faster than the previous
	   on most CPUs for VBITS >= 128 */

	/* use only for VWORDS = 2, 4, 6, or 8
	   processes 128 bits per loop */
	   
	uint32 w_x, w_y;

	for (w_x = 0; w_x < VWORDS; w_x += 2) {
		for (w_y = 0; w_y < VWORDS; w_y += 2) {

			int i, j, k;
			uint64 c[16][256][2];
			memset(c, 0, 16 * 256 * 2 * sizeof(uint64));

// #pragma omp parallel for private(j, k) reduction(^:c)
			for (i = 0; i < n; i++) {
				uint64 xi[2], yi[2];
				xi[0] = x[i].w[w_x];
				xi[1] = x[i].w[w_x + 1];
				yi[0] = y[i].w[w_y];
				yi[1] = y[i].w[w_y + 1];

				for (j = 0; j < 16; j++) { 
					k = (xi[(j >> 3)] >> (8*(j & 7))) & 255; 
					c[j][k][0] ^= yi[0];
					c[j][k][1] ^= yi[1];
				}
			}
			
			// Now combine the table entries

			for (i = 0; i < 8; i++) {
				uint64 a[16][2];
				memset(a, 0, 16 * 2 * sizeof(uint64));

// #pragma omp parallel for private(k) reduction(^:a)
				for (j = 0; j < 256; j++) {
					if ((j >> i) & 1) {
						for (k = 0; k < 16; k++) { 
							a[k][0] ^= c[k][j][0];
							a[k][1] ^= c[k][j][1];
						}
					}
				}

// #pragma omp parallel for 
				for (j = 0; j < 16; j++) {
					xy[64 * w_x + 8 * j + i].w[w_y] = a[j][0];
					xy[64 * w_x + 8 * j + i].w[w_y + 1] = a[j][1];
				}
			}
		}
	}
}

/*-------------------------------------------------------------------*/

void vv_mul_BxN_NxB(packed_matrix_t *matrix,
		   void *x_in, void *y_in,
		   v_t *xy, uint32 n) {

	v_t *x = (v_t *)x_in;
	v_t *y = (v_t *)y_in;
	int i, threads, size;
	v_t *xytmp;

#ifdef HAVE_OMP
	threads = omp_get_max_threads();
#else 
	threads = 1;
#endif
	xytmp = (v_t *) malloc(threads * VBITS * sizeof(v_t));
	size = n / threads + 1;

#pragma omp parallel for schedule(static, 1)
	for (i = 0; i < threads; i++) {
		uint32 my_n;
		if (i == threads-1) my_n = n + size - threads * size;
		else my_n = size;
#if VWORDS == 2 || VWORDS == 4 || VWORDS == 6 || VWORDS == 8
		mul_BxN_NxB_2(x + i * size, y + i * size, xytmp + i * VBITS, my_n);
#else
		mul_BxN_NxB(x + i * size, y + i * size, xytmp + i * VBITS, my_n);
#endif
	}

	vv_clear(xy, VBITS);
	for (i = 0; i < threads; i++) vv_xor(xy, xytmp + i * VBITS, VBITS);

#ifdef HAVE_MPI
	/* combine the results across an entire MPI row */

	global_xor(xy, xytmp, VBITS, matrix->mpi_ncols,
			matrix->mpi_la_col_rank,
			matrix->mpi_word, matrix->mpi_la_row_grid);

	/* combine the results across an entire MPI column */
    
	global_xor(xytmp, xy, VBITS, matrix->mpi_nrows,
			matrix->mpi_la_row_rank,
			matrix->mpi_word, matrix->mpi_la_col_grid);    
#endif
	free(xytmp);
}
