/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/24/08

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.
--------------------------------------------------------------------*/

#include "lanczos.h"
#include "util.h"
#include "qs.h"

#define BIT(x) ((uint64)(1) << (x))

static const uint64 bitmask[64] = {
	BIT( 0), BIT( 1), BIT( 2), BIT( 3), BIT( 4), BIT( 5), BIT( 6), BIT( 7),
	BIT( 8), BIT( 9), BIT(10), BIT(11), BIT(12), BIT(13), BIT(14), BIT(15),
	BIT(16), BIT(17), BIT(18), BIT(19), BIT(20), BIT(21), BIT(22), BIT(23),
	BIT(24), BIT(25), BIT(26), BIT(27), BIT(28), BIT(29), BIT(30), BIT(31),
	BIT(32), BIT(33), BIT(34), BIT(35), BIT(36), BIT(37), BIT(38), BIT(39),
	BIT(40), BIT(41), BIT(42), BIT(43), BIT(44), BIT(45), BIT(46), BIT(47),
	BIT(48), BIT(49), BIT(50), BIT(51), BIT(52), BIT(53), BIT(54), BIT(55),
	BIT(56), BIT(57), BIT(58), BIT(59), BIT(60), BIT(61), BIT(62), BIT(63),
};

/* for matrices of dimension exceeding MIN_POST_LANCZOS_DIM,
   the first POST_LANCZOS_ROWS rows are handled in a separate
   Gauss elimination phase after the Lanczos iteration
   completes. This means the lanczos code will produce about
   64 - POST_LANCZOS_ROWS dependencies on average. 
   
   The code will still work if POST_LANCZOS_ROWS is 0, but I 
   don't know why you would want to do that. The first rows are 
   essentially completely dense, and removing them from the main 
   Lanczos iteration greatly reduces the amount of arithmetic 
   in a matrix multiply, as well as the memory footprint of 
   the matrix */

#define POST_LANCZOS_ROWS 48
#define MIN_POST_LANCZOS_DIM 10000

/*-------------------------------------------------------------------*/
static uint64 * yafu_form_post_lanczos_matrix(fact_obj_t *obj, uint32 *nrows, 
				uint32 *dense_rows_out, uint32 ncols, 
				qs_la_col_t *cols) {

	uint32 i, j, k;
	uint32 num_dense_rows = *dense_rows_out;
	uint32 dense_row_words;
	uint32 new_dense_rows;
	uint32 new_dense_row_words;
	uint32 final_dense_row_words;
	uint64 mask;
	uint64 *submatrix;
#if BITS_PER_DIGIT == 32
	z tmp;
#else
	z32 tmp;
#endif

	/* if the matrix is going to have cache blocking applied,
	   proceed but do not form a post-Lanczos matrix if one
	   is not desired. We have to do this because the block
	   matrix multiply expects the number of dense rows to be
	   a multiple of 64.

	   Otherwise, don't do anything if the Lanczos iteration 
	   would finish quickly */

	submatrix = NULL;
	if (ncols >= MIN_NCOLS_TO_PACK ||
	    (POST_LANCZOS_ROWS > 0 && ncols >= MIN_POST_LANCZOS_DIM)) {

		if (POST_LANCZOS_ROWS > 0) {
			if (VFLAG > 0)
				printf("saving the first %u matrix rows "
					"for later\n", POST_LANCZOS_ROWS);
			logprint(obj->logfile, "saving the first %u matrix rows "
					"for later\n", POST_LANCZOS_ROWS);
			submatrix = (uint64 *)xmalloc(ncols * sizeof(uint64));
		}
	}
	else {
		return NULL;
	}

#if BITS_PER_DIGIT == 32
	zInit(&tmp);
#else
	zInit32(&tmp);
#endif

	mask = (uint64)(-1) >> (64 - POST_LANCZOS_ROWS);
	dense_row_words = (num_dense_rows + 31) / 32;
	/* we will be removing the first POST_LANCZOS_ROWS rows
	   from the matrix entirely, and packing together the
	   next few rows. The matrix may have dense rows already, 
	   or these rows may be partially or completely sparse, 
	   in which case we'll have to pack them manually. After
	   the post-lanczos rows are removed, the number of dense 
	   rows remaining is a multiple of 64 (minimum of 64) */

	new_dense_rows = MAX(num_dense_rows, POST_LANCZOS_ROWS);
	new_dense_rows += 64 - (new_dense_rows - POST_LANCZOS_ROWS) % 64;
	new_dense_row_words = (new_dense_rows + 31) / 32;
	final_dense_row_words = (new_dense_rows - POST_LANCZOS_ROWS) / 32;

	for (i = 0; i < ncols; i++) {
		uint32 curr_weight = cols[i].weight;
		uint32 *curr_row = cols[i].data;

		/* build up a bitfield of the rows that will be
		   stored in packed format. Start with the rows
		   that are already packed */

		for (j = 0; j < dense_row_words; j++)
			tmp.val[j] = curr_row[curr_weight + j];

		/* add in the rows from the sparse part of the matrix.
		   Entries from these rows are either added to the
		   new dense bitfield, or moved to fill the holes
		   created by packing the first few sparse rows. In 
		   the latter case, the row index must be biased to 
		   reflect the removed rows */

		for (; j < new_dense_row_words; j++)
			tmp.val[j] = 0;

		for (j = k = 0; j < curr_weight; j++) {
			uint32 curr_index = curr_row[j];

			if (curr_index < new_dense_rows)
				tmp.val[curr_index / 32] |= 
						bitmask[curr_index % 32];
			else
				curr_row[k++] = curr_index - POST_LANCZOS_ROWS;
		}

		//tmp.nwords = new_dense_row_words;
		tmp.size =  new_dense_row_words;

#if POST_LANCZOS_ROWS > 0
		/* remove the first POST_LANCZOS_ROWS bits from
		   the bitfield */
		submatrix[i] = ((uint64)tmp.val[0] |
				(uint64)tmp.val[1] << 32) & mask;
#endif

		/* move the rest of the bitfield and repack the (hopefully
		   shorter) current column in the heap */
		cols[i].weight = k;
		if (k + final_dense_row_words > 0) {
			cols[i].data = (uint32 *)xrealloc(curr_row, (k + 
						final_dense_row_words) * 
						sizeof(uint32));

#if BITS_PER_DIGIT == 32
		zShiftRight(&tmp,&tmp,POST_LANCZOS_ROWS);
#else
		zShiftRight32(&tmp,&tmp,POST_LANCZOS_ROWS);
#endif
			memcpy(cols[i].data + k, tmp.val, 
						final_dense_row_words * 
						sizeof(uint32));
		}
		else {
			free(cols[i].data);
			cols[i].data = NULL;
		}
	}

	*nrows -= POST_LANCZOS_ROWS;
	*dense_rows_out = new_dense_rows - POST_LANCZOS_ROWS;
	count_qs_matrix_nonzero(obj, *nrows, *dense_rows_out, ncols, cols);

#if BITS_PER_DIGIT == 32
	zFree(&tmp);
#else
	zFree32(&tmp);
#endif
	return submatrix;
}

/*-------------------------------------------------------------------*/
static void yafu_mul_64x64_64x64(uint64 *a, uint64 *b, uint64 *c ) {

	/* c[][] = x[][] * y[][], where all operands are 64 x 64
	   (i.e. contain 64 words of 64 bits each). The result
	   may overwrite a or b. */

	uint64 ai, bj, accum;
	uint64 tmp[64];
	uint32 i, j;

	for (i = 0; i < 64; i++) {
		j = 0;
		accum = 0;
		ai = a[i];

		while (ai) {
			bj = b[j];
			if (ai & 1)
				accum ^= bj;
			ai >>= 1;
			j++;
		}

		tmp[i] = accum;
	}
	memcpy(c, tmp, sizeof(tmp));
}

/*-----------------------------------------------------------------------*/
static void yafu_transpose_64x64(uint64 *a, uint64 *b) {

	uint32 i, j;
	uint64 tmp[64] = {0};

	for (i = 0; i < 64; i++) {
		uint64 word = a[i];
		uint64 mask = bitmask[i];
		for (j = 0; j < 64; j++) {
			if (word & bitmask[j])
				tmp[j] |= mask;
		}
	}
	memcpy(b, tmp, sizeof(tmp));
}

/*-------------------------------------------------------------------*/
void yafu_mul_Nx64_64x64_acc(uint64 *v, uint64 *x,
			uint64 *y, uint32 n) {

	/* let v[][] be a n x 64 matrix with elements in GF(2), 
	   represented as an array of n 64-bit words. Let c[][]
	   be an 8 x 256 scratch matrix of 64-bit words.
	   This code multiplies v[][] by the 64x64 matrix 
	   x[][], then XORs the n x 64 result into y[][] */

	uint32 i, j, k;
	uint64 c[8 * 256];

	/* fill c[][] with a bunch of "partial matrix multiplies". 
	   For 0<=i<256, the j_th row of c[][] contains the matrix 
	   product

	   	( i << (8*j) ) * x[][]

	   where the quantity in parentheses is considered a 
	   1 x 64 vector of elements in GF(2). The resulting
	   table will dramatically speed up matrix multiplies
	   by x[][]. */

	for (i = 0; i < 8; i++) {
		uint64 *xtmp = x + 8 * i;
		uint64 *ctmp = c + 256 * i;

		for (j = 0; j < 256; j++) {
			uint64 accum = 0;
			uint32 index = j;

			for (k = 0; k < 8; k++) {
				if (index & ((uint32)1 << k))
					accum ^= xtmp[k];
			}
			ctmp[j] = accum;
		}
	}

#if (defined(__GNUC__) || defined(__ICL)) && \
	defined(__i386__) && defined(HAS_MMX)
	i = 0;
	asm volatile(".p2align 4,,7                        \n\t"
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
			:"r"(v), "r"(c), "r"(y), "g"(n)
			:"%eax", "%ecx", "%mm0", "%mm1", "memory");

#elif defined(_MSC_VER) && !defined(_WIN64)
	i = 0;
	__asm
	{
		push	ebx
		mov	edi,y
		lea	ebx,c
		mov	esi,v
		mov	ecx,i
		align 16
	L0:	movq	mm0,[edi+ecx*8]
		mov	eax,[esi+ecx*8]
		inc	ecx
		movzx	edx, al
		movq	mm1,[ebx+edx*8]
		movzx	edx,ah
		pxor	mm1,[1*256*8+ebx+edx*8]
		shr	eax,16
		movzx	edx,al
		pxor	mm1,[2*256*8+ebx+edx*8]
		movzx	edx,ah
		pxor	mm1,[3*256*8+ebx+edx*8]
		mov	eax,[4-8+esi+ecx*8]
		movzx	edx,al
		pxor	mm1,[4*256*8+ebx+edx*8]
		movzx	edx,ah
		shr	eax,16
		cmp	ecx,n
		pxor	mm1,[5*256*8+ebx+edx*8]
		movzx	edx,al
		pxor	mm1,[6*256*8+ebx+edx*8]
		movzx	edx,ah
		pxor	mm1,[7*256*8+ebx+edx*8]
		pxor	mm1, mm0
		movq	[-8+edi+ecx*8],mm1
		jne	L0
		pop	ebx
		emms
	}
#else
	for (i = 0; i < n; i++) {
		uint64 word = v[i];
		y[i] ^=  c[ 0*256 + ((uint8)(word >>  0)) ]
		       ^ c[ 1*256 + ((uint8)(word >>  8)) ]
		       ^ c[ 2*256 + ((uint8)(word >> 16)) ]
		       ^ c[ 3*256 + ((uint8)(word >> 24)) ]
		       ^ c[ 4*256 + ((uint8)(word >> 32)) ]
		       ^ c[ 5*256 + ((uint8)(word >> 40)) ]
		       ^ c[ 6*256 + ((uint8)(word >> 48)) ]
		       ^ c[ 7*256 + ((uint8)(word >> 56)) ];
	}
#endif
}

/*-------------------------------------------------------------------*/
void yafu_mul_64xN_Nx64(uint64 *x, uint64 *y,
		   uint64 *xy, uint32 n) {

	/* Let x and y be n x 64 matrices. This routine computes
	   the 64 x 64 matrix xy[][] given by transpose(x) * y */

	uint32 i;
	uint64 c[8 * 256] = {0};

	memset(xy, 0, 64 * sizeof(uint64));

#if (defined(__GNUC__) || defined(__ICL)) && \
	defined(__i386__) && defined(HAS_MMX)
	i = 0;
	asm volatile(".p2align 4,,7                        \n\t"
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

#elif defined(_MSC_VER) && !defined(_WIN64)
	i = 0;
	__asm
	{
		push	ebx
		mov	edi,y
		lea	ebx,c
		mov	esi,x
		mov	ecx,i
		align 16
	L0:	movq	mm0,[edi+ecx*8]
		mov	eax,[esi+ecx*8]
		inc	ecx
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[ebx+edx*8]
		movq	[ebx+edx*8],mm1
		movzx	edx,ah
		movq	mm1, mm0
		pxor	mm1,[1*256*8+ebx+edx*8]
		movq	[1*256*8+ebx+edx*8],mm1
		shr	eax,16
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[2*256*8+ebx+edx*8]
		movq	[2*256*8+ebx+edx*8],mm1
		movzx	edx,ah
		movq	mm1,mm0
		pxor	mm1,[3*256*8+ebx+edx*8]
		movq	[3*256*8+ebx+edx*8],mm1
		mov	eax,[4-8+esi+ecx*8]
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[4*256*8+ebx+edx*8]
		movq	[4*256*8+ebx+edx*8],mm1
		movzx	edx,ah
		shr	eax,16
		cmp	ecx,n
		movq	mm1,mm0
		pxor	mm1,[5*256*8+ebx+edx*8]
		movq	[5*256*8+ebx+edx*8],mm1
		movzx	edx,al
		movq	mm1,mm0
		pxor	mm1,[6*256*8+ebx+edx*8]
		movq	[6*256*8+ebx+edx*8],mm1
		movzx	edx,ah
		pxor	mm0,[7*256*8+ebx+edx*8]
		movq	[7*256*8+ebx+edx*8],mm0
		jne	L0
		emms
		pop	ebx
	}
#else

	for (i = 0; i < n; i++) {
		uint64 xi = x[i];
		uint64 yi = y[i];
		c[ 0*256 + ((uint8) xi       ) ] ^= yi;
		c[ 1*256 + ((uint8)(xi >>  8)) ] ^= yi;
		c[ 2*256 + ((uint8)(xi >> 16)) ] ^= yi;
		c[ 3*256 + ((uint8)(xi >> 24)) ] ^= yi;
		c[ 4*256 + ((uint8)(xi >> 32)) ] ^= yi;
		c[ 5*256 + ((uint8)(xi >> 40)) ] ^= yi;
		c[ 6*256 + ((uint8)(xi >> 48)) ] ^= yi;
		c[ 7*256 + ((uint8)(xi >> 56)) ] ^= yi;
	}
#endif

	for(i = 0; i < 8; i++) {

		uint32 j;
		uint64 a0, a1, a2, a3, a4, a5, a6, a7;

		a0 = a1 = a2 = a3 = 0;
		a4 = a5 = a6 = a7 = 0;

		for (j = 0; j < 256; j++) {
			if ((j >> i) & 1) {
				a0 ^= c[0*256 + j];
				a1 ^= c[1*256 + j];
				a2 ^= c[2*256 + j];
				a3 ^= c[3*256 + j];
				a4 ^= c[4*256 + j];
				a5 ^= c[5*256 + j];
				a6 ^= c[6*256 + j];
				a7 ^= c[7*256 + j];
			}
		}

		xy[ 0] = a0; xy[ 8] = a1; xy[16] = a2; xy[24] = a3;
		xy[32] = a4; xy[40] = a5; xy[48] = a6; xy[56] = a7;
		xy++;
	}
}

/*-------------------------------------------------------------------*/
static uint32 yafu_find_nonsingular_sub(fact_obj_t *obj,
				uint64 *t, uint32 *s, 
				uint32 *last_s, uint32 last_dim, 
				uint64 *w) {

	/* given a 64x64 matrix t[][] (i.e. sixty-four
	   64-bit words) and a list of 'last_dim' column 
	   indices enumerated in last_s[]: 
	   
	     - find a submatrix of t that is invertible 
	     - invert it and copy to w[][]
	     - enumerate in s[] the columns represented in w[][] */

	uint32 i, j;
	uint32 dim;
	uint32 cols[64];
	uint64 M[64][2];
	uint64 mask, *row_i, *row_j;
	uint64 m0, m1;

	/* M = [t | I] for I the 64x64 identity matrix */

	for (i = 0; i < 64; i++) {
		M[i][0] = t[i]; 
		M[i][1] = bitmask[i];
	}

	/* put the column indices from last_s[] into the
	   back of cols[], and copy to the beginning of cols[]
	   any column indices not in last_s[] */

	mask = 0;
	for (i = 0; i < last_dim; i++) {
		cols[63 - i] = last_s[i];
		mask |= bitmask[last_s[i]];
	}
	for (i = j = 0; i < 64; i++) {
		if (!(mask & bitmask[i]))
			cols[j++] = i;
	}

	/* compute the inverse of t[][] */

	for (i = dim = 0; i < 64; i++) {
	
		/* find the next pivot row and put in row i */

		mask = bitmask[cols[i]];
		row_i = M[cols[i]];

		for (j = i; j < 64; j++) {
			row_j = M[cols[j]];
			if (row_j[0] & mask) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		/* if a pivot row was found, eliminate the pivot
		   column from all other rows */

		if (j < 64) {
			for (j = 0; j < 64; j++) {
				row_j = M[cols[j]];
				if ((row_i != row_j) && (row_j[0] & mask)) {
					row_j[0] ^= row_i[0];
					row_j[1] ^= row_i[1];
				}
			}

			/* add the pivot column to the list of 
			   accepted columns */

			s[dim++] = cols[i];
			continue;
		}

		/* otherwise, use the right-hand half of M[]
		   to compensate for the absence of a pivot column */

		for (j = i; j < 64; j++) {
			row_j = M[cols[j]];
			if (row_j[1] & mask) {
				m0 = row_j[0];
				m1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m0; 
				row_i[1] = m1;
				break;
			}
		}
				
		if (j == 64) {
			printf("lanczos error: submatrix "
					"is not invertible\n");
			logprint(obj->logfile, "lanczos error: submatrix "
					"is not invertible\n");
			return 0;
		}
			
		/* eliminate the pivot column from the other rows
		   of the inverse */

		for (j = 0; j < 64; j++) {
			row_j = M[cols[j]];
			if ((row_i != row_j) && (row_j[1] & mask)) {
				row_j[0] ^= row_i[0];
				row_j[1] ^= row_i[1];
			}
		}

		/* wipe out the pivot row */

		row_i[0] = row_i[1] = 0;
	}

	/* the right-hand half of M[] is the desired inverse */
	
	for (i = 0; i < 64; i++) 
		w[i] = M[i][1];

	return dim;
}

/*-----------------------------------------------------------------------*/
static void yafu_transpose_vector(uint32 ncols, uint64 *v, uint64 **trans) {

	/* Hideously inefficent routine to transpose a
	   vector v[] of 64-bit words into a 2-D array
	   trans[][] of 64-bit words */

	uint32 i, j;
	uint32 col;
	uint64 mask, word;

	for (i = 0; i < ncols; i++) {
		col = i / 64;
		mask = bitmask[i % 64];
		word = v[i];
		j = 0;
		while (word) {
			if (word & 1)
				trans[j][col] |= mask;
			word = word >> 1;
			j++;
		}
	}
}

/*-----------------------------------------------------------------------*/
static uint32 yafu_combine_cols(uint32 ncols, 
			uint64 *x, uint64 *v, 
			uint64 *ax, uint64 *av) {

	/* Once the block Lanczos iteration has finished, 
	   x[] and v[] will contain mostly nullspace vectors
	   between them, as well as possibly some columns
	   that are linear combinations of nullspace vectors.
	   Given vectors ax[] and av[] that are the result of
	   multiplying x[] and v[] by the matrix, this routine 
	   will use Gauss elimination on the columns of [ax | av] 
	   to find all of the linearly dependent columns. The
	   column operations needed to accomplish this are mir-
	   rored in [x | v] and the columns that are independent
	   are skipped. Finally, the dependent columns are copied
	   back into x[] and represent the nullspace vector output
	   of the block Lanczos code. */

	uint32 i, j, k, bitpos, col, col_words;
	uint64 mask;
	uint64 *matrix[128], *amatrix[128], *tmp;

	col_words = (ncols + 63) / 64;

	for (i = 0; i < 128; i++) {
		matrix[i] = (uint64 *)xcalloc((size_t)col_words, 
					     sizeof(uint64));
		amatrix[i] = (uint64 *)xcalloc((size_t)col_words, 
					      sizeof(uint64));
	}

	/* operations on columns can more conveniently become 
	   operations on rows if all the vectors are first
	   transposed */

	yafu_transpose_vector(ncols, x, matrix);
	yafu_transpose_vector(ncols, ax, amatrix);
	yafu_transpose_vector(ncols, v, matrix + 64);
	yafu_transpose_vector(ncols, av, amatrix + 64);

	/* Keep eliminating rows until the unprocessed part
	   of amatrix[][] is all zero. The rows where this
	   happens correspond to linearly dependent vectors
	   in the nullspace */

	for (i = bitpos = 0; i < 128 && bitpos < ncols; bitpos++) {

		/* find the next pivot row */

		mask = bitmask[bitpos % 64];
		col = bitpos / 64;
		for (j = i; j < 128; j++) {
			if (amatrix[j][col] & mask) {
				tmp = matrix[i];
				matrix[i] = matrix[j];
				matrix[j] = tmp;
				tmp = amatrix[i];
				amatrix[i] = amatrix[j];
				amatrix[j] = tmp;
				break;
			}
		}
		if (j == 128)
			continue;

		/* a pivot was found; eliminate it from the
		   remaining rows */

		for (j++; j < 128; j++) {
			if (amatrix[j][col] & mask) {

				/* Note that the entire row, *not*
				   just the nonzero part of it, must
				   be eliminated; this is because the
				   corresponding (dense) row of matrix[][]
				   must have the same operation applied */

				for (k = 0; k < col_words; k++) {
					amatrix[j][k] ^= amatrix[i][k];
					matrix[j][k] ^= matrix[i][k];
				}
			}
		}
		i++;
	}

	/* transpose rows i to 64 back into x[]. Pack the
	   dependencies into the low-order bits of x[] */

	for (j = 0; j < ncols; j++) {
		uint64 word = 0;

		col = j / 64;
		mask = bitmask[j % 64];

		for (k = i; k < 64; k++) {
			if (matrix[k][col] & mask)
				word |= bitmask[k - i];
		}
		x[j] = word;
	}

	for (j = 0; j < 128; j++) {
		free(matrix[j]);
		free(amatrix[j]);
	}

	if (i > 64)
		return 0;
	return 64 - i;
}

/*-----------------------------------------------------------------------*/
static void yafu_dump_lanczos_state(fact_obj_t *obj, 
			uint64 *x, uint64 **vt_v0, uint64 **v, 
			uint64 **vt_a_v, uint64 **vt_a2_v, uint64 **winv,
			uint32 n, uint32 dim_solved, uint32 iter,
			uint32 s[2][64], uint32 dim1) {

	char buf[256];
	FILE *dump_fp;

	sprintf(buf, "%s.chk", obj->savefile_name);
	dump_fp = fopen(buf, "wb");
	if (dump_fp == NULL) {
		printf("error: cannot open matrix checkpoint file\n");
		exit(-1);
	}

	fwrite(&n, sizeof(uint32), (size_t)1, dump_fp);
	fwrite(&dim_solved, sizeof(uint32), (size_t)1, dump_fp);
	fwrite(&iter, sizeof(uint32), (size_t)1, dump_fp);

	fwrite(vt_a_v[1], sizeof(uint64), (size_t)64, dump_fp);
	fwrite(vt_a2_v[1], sizeof(uint64), (size_t)64, dump_fp);
	fwrite(winv[1], sizeof(uint64), (size_t)64, dump_fp);
	fwrite(winv[2], sizeof(uint64), (size_t)64, dump_fp);
	fwrite(vt_v0[0], sizeof(uint64), (size_t)64, dump_fp);
	fwrite(vt_v0[1], sizeof(uint64), (size_t)64, dump_fp);
	fwrite(vt_v0[2], sizeof(uint64), (size_t)64, dump_fp);
	fwrite(s[1], sizeof(uint32), (size_t)64, dump_fp);
	fwrite(&dim1, sizeof(uint32), (size_t)1, dump_fp);

	fwrite(x, sizeof(uint64), (size_t)n, dump_fp);
	fwrite(v[0], sizeof(uint64), (size_t)n, dump_fp);
	fwrite(v[1], sizeof(uint64), (size_t)n, dump_fp);
	fwrite(v[2], sizeof(uint64), (size_t)n, dump_fp);
	fclose(dump_fp);
}

/*-----------------------------------------------------------------------*/
static void yafu_read_lanczos_state(fact_obj_t *obj, 
			uint64 *x, uint64 **vt_v0, uint64 **v, 
			uint64 **vt_a_v, uint64 **vt_a2_v, uint64 **winv,
			uint32 n, uint32 *dim_solved, uint32 *iter,
			uint32 s[2][64], uint32 *dim1) {

	uint32 read_n;
	uint32 status;
	char buf[256];
	FILE *dump_fp;

	sprintf(buf, "%s.chk", obj->savefile_name);
	dump_fp = fopen(buf, "rb");
	if (dump_fp == NULL) {
		printf("error: cannot open matrix checkpoint file\n");
		exit(-1);
	}

	status = 1;
	fread(&read_n, sizeof(uint32), (size_t)1, dump_fp);
	if (read_n != n) {
		printf("error: unexpected vector size\n");
		exit(-1);
	}
	status &= (fread(dim_solved, sizeof(uint32), (size_t)1, dump_fp) == 1);
	status &= (fread(iter, sizeof(uint32), (size_t)1, dump_fp) == 1);

	status &= (fread(vt_a_v[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_a2_v[1], sizeof(uint64), (size_t)64, dump_fp)==64);
	status &= (fread(winv[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(winv[2], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_v0[0], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_v0[1], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(vt_v0[2], sizeof(uint64), (size_t)64, dump_fp) == 64);
	status &= (fread(s[1], sizeof(uint32), (size_t)64, dump_fp) == 64);
	status &= (fread(dim1, sizeof(uint32), (size_t)1, dump_fp) == 1);

	status &= (fread(x, sizeof(uint64), (size_t)n, dump_fp) == n);
	status &= (fread(v[0], sizeof(uint64), (size_t)n, dump_fp) == n);
	status &= (fread(v[1], sizeof(uint64), (size_t)n, dump_fp) == n);
	status &= (fread(v[2], sizeof(uint64), (size_t)n, dump_fp) == n);

	fclose(dump_fp);
	if (status == 0) {
		printf("error: checkpoint recovery failed\n");
		exit(-1);
	}
}

/*-----------------------------------------------------------------------*/
static void yafu_init_lanczos_state(fact_obj_t *obj, 
			packed_matrix_t *packed_matrix,
			uint64 *x, uint64 *v0, uint64 **vt_v0, uint64 **v, 
			uint64 **vt_a_v, uint64 **vt_a2_v, uint64 **winv,
			uint32 n, uint32 s[2][64], uint32 *dim1) {

	uint32 i;

	/* The computed solution 'x' starts off random,
	   and v[0] starts off as B*x. This initial copy
	   of v[0] must be saved off separately */

	for (i = 0; i < n; i++) {
		x[i] = v[0][i] = 
			  (uint64)(get_rand(&obj->seed1, &obj->seed2)) << 32 |
		          (uint64)(get_rand(&obj->seed1, &obj->seed2));
	}

	yafu_mul_MxN_Nx64(packed_matrix, v[0], v[1]);
	yafu_mul_trans_MxN_Nx64(packed_matrix, v[1], v[0]);
	memcpy(v0, v[0], n * sizeof(uint64));

	/* Subscripts larger than zero represent past versions of 
	   these quantities, which start off empty (except for the 
	   past version of s[], which contains all the column 
	   indices) */
	   
	memset(v[1], 0, n * sizeof(uint64));
	memset(v[2], 0, n * sizeof(uint64));
	for (i = 0; i < 64; i++) {
		s[1][i] = i;
		vt_a_v[1][i] = 0;
		vt_a2_v[1][i] = 0;
		winv[1][i] = 0;
		winv[2][i] = 0;
		vt_v0[0][i] = 0;
		vt_v0[1][i] = 0;
		vt_v0[2][i] = 0;
	}
	*dim1 = 64;
}

/*-----------------------------------------------------------------------*/
static uint64 * yafu_block_lanczos_core(fact_obj_t *obj, 
				packed_matrix_t *packed_matrix,
				uint32 *num_deps_found,
				uint64 *post_lanczos_matrix,
				uint32 dump_interval) {
	
	/* Solve Bx = 0 for some nonzero x; the computed
	   solution, containing up to 64 of these nullspace
	   vectors, is returned */

	uint32 n = packed_matrix->ncols;
	uint64 *vnext, *v[3], *x, *v0;
	uint64 *winv[3], *vt_v0_next;
	uint64 *vt_a_v[2], *vt_a2_v[2], *vt_v0[3];
	uint64 *scratch;
	uint64 *tmp;
	uint32 s[2][64];
	uint64 d[64], e[64], f[64], f2[64];
	uint32 i, iter;
	uint32 dim0, dim1;
	uint64 mask0, mask1;

	uint32 dim_solved = 0;
	uint32 first_dim_solved = 0;
	uint32 report_interval = 0;
	uint32 next_report = 0;
	uint32 next_dump = 0;
	time_t first_time;

	if (packed_matrix->num_threads > 1)
	{	
		if (VFLAG > 0)
			printf("commencing Lanczos iteration (%u threads)\n",
					packed_matrix->num_threads);
		logprint(obj->logfile, "commencing Lanczos iteration (%u threads)\n",
					packed_matrix->num_threads);
	}
	else
	{
		if (VFLAG > 0)
			printf("commencing Lanczos iteration\n");
		logprint(obj->logfile, "commencing Lanczos iteration\n");
	}

	/* allocate all the 64x64 variables */

	winv[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	winv[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	winv[2] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a_v[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a_v[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a2_v[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_a2_v[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0[0] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0[1] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0[2] = (uint64 *)xmalloc(64 * sizeof(uint64));
	vt_v0_next = (uint64 *)xmalloc(64 * sizeof(uint64));

	/* allocate all of the size-n variables except v0,
	   which will be freed if it's not needed */

	v[0] = (uint64 *)xmalloc(n * sizeof(uint64));
	v[1] = (uint64 *)xmalloc(n * sizeof(uint64));
	v[2] = (uint64 *)xmalloc(n * sizeof(uint64));
	vnext = (uint64 *)xmalloc(n * sizeof(uint64));
	x = (uint64 *)xmalloc(n * sizeof(uint64));
	scratch = (uint64 *)xmalloc(n * sizeof(uint64));
	v0 = NULL;

	if (VFLAG > 0)
		printf("memory use: %.1f MB\n", (double)
			((6 * n * sizeof(uint64) +
			 yafu_packed_matrix_sizeof(packed_matrix))) / 1048576);
	logprint(obj->logfile, "memory use: %.1f MB\n", (double)
			((6 * n * sizeof(uint64) +
			 yafu_packed_matrix_sizeof(packed_matrix))) / 1048576);

	/* initialize */

	iter = 0;
	dim0 = 0;
	next_dump = dump_interval;

	if (obj->flags & MSIEVE_FLAG_NFS_LA_RESTART) {
		yafu_read_lanczos_state(obj, x, vt_v0, v, vt_a_v, vt_a2_v,
				winv, n, &dim_solved, &iter, s, &dim1);
		if (VFLAG > 0)
			printf("restarting at iteration %u (dim = %u)\n",
				iter, dim_solved);
		logprint(obj->logfile, "restarting at iteration %u (dim = %u)\n",
				iter, dim_solved);
	}
	else {
		v0 = (uint64 *)xmalloc(n * sizeof(uint64));
		yafu_init_lanczos_state(obj, packed_matrix, x, v0, vt_v0, v, 
				vt_a_v, vt_a2_v, winv, n, s, &dim1);
	}

	mask1 = 0;
	for (i = 0; i < dim1; i++)
		mask1 |= bitmask[s[1][i]];

	/* determine if the solver will run long enough that
	   it would be worthwhile to report progress */

	first_time = time(NULL);
	if (n > 60000 &&
	    obj->flags & (MSIEVE_FLAG_USE_LOGFILE |
	    		  MSIEVE_FLAG_LOG_TO_STDOUT)) {
		if (n > 1000000)
			report_interval = 200;
		else if (n > 500000)
			report_interval = 500;
		else if (n > 100000)
			report_interval = 2000;
		else
			report_interval = 8000;
		first_dim_solved = dim_solved;
		next_report = dim_solved + report_interval;
	}

	/* perform the iteration */

	while (1) {
		iter++;

		/* multiply the current v[0] by a symmetrized
		   version of B, or B'B (apostrophe means 
		   transpose). Use "A" to refer to B'B  */

		yafu_mul_MxN_Nx64(packed_matrix, v[0], scratch);
		yafu_mul_trans_MxN_Nx64(packed_matrix, scratch, vnext);

		/* compute v0'*A*v0 and (A*v0)'(A*v0) */

		yafu_mul_64xN_Nx64(v[0], vnext, vt_a_v[0], n);
		yafu_mul_64xN_Nx64(vnext, vnext, vt_a2_v[0], n);

		/* if the former is orthogonal to itself, then
		   the iteration has finished */

		for (i = 0; i < 64; i++) {
			if (vt_a_v[0][i] != 0)
				break;
		}
		if (i == 64)
			break;

		/* Find the size-'dim0' nonsingular submatrix
		   of v0'*A*v0, invert it, and list the column
		   indices present in the submatrix */

		dim0 = yafu_find_nonsingular_sub(obj, vt_a_v[0], s[0], 
					    s[1], dim1, winv[0]);
		if (dim0 == 0)
			break;

		/* mask0 contains one set bit for every column
		   that participates in the inverted submatrix
		   computed above */

		mask0 = 0;
		for (i = 0; i < dim0; i++)
			mask0 |= bitmask[s[0][i]];

		/* The block Lanczos recurrence depends on all columns
		   of v'Av appearing in the current and/or previous iteration. 
		   Verify that condition here
		  
		   Note that the test only applies if this is not the
		   last Lanczos iteration. I'm not sure that this is right, 
		   but the test fails on the last iteration much more often 
		   than would be expected by chance alone, yet ignoring
		   the failure still produces good dependencies 
		   
		   Note that the last iteration typically has dim_solved
		   slightly less than the number of rows, not the number
		   of columns (=n) */
	
		if (dim_solved < packed_matrix->nrows - 64) {
			if ((mask0 | mask1) != (uint64)(-1)) {
				printf("lanczos error (dim = %u): "
						"not all columns used\n",
						dim_solved);
				logprint(obj->logfile, "lanczos error (dim = %u): "
						"not all columns used\n",
						dim_solved);
				dim0 = 0;
				break;
			}
		}

		/* begin the computation of the next v. First mask
		   off the vectors that are included in this iteration */

		if (mask0 != (uint64)(-1)) {
			for (i = 0; i < n; i++)
				vnext[i] = vnext[i] & mask0;
		}

		/* begin the computation of the next v' * v0. For 
		   the first three iterations, this requires a full 
		   inner product. For all succeeding iterations, the 
		   next v' * v0 is the sum of three 64x64 products 
		   and is stored in vt_v0_next. */

		if (iter < 4) {
			yafu_mul_64xN_Nx64(v[0], v0, vt_v0[0], n);
		}
		else if (iter == 4) {
			free(v0);
			v0 = NULL;
		}

		/* compute d, fold it into vnext and update v'*v0 */

		for (i = 0; i < 64; i++)
			d[i] = (vt_a2_v[0][i] & mask0) ^ vt_a_v[0][i];

		yafu_mul_64x64_64x64(winv[0], d, d);

		for (i = 0; i < 64; i++)
			d[i] = d[i] ^ bitmask[i];

		yafu_mul_Nx64_64x64_acc(v[0], d, vnext, n);

		yafu_transpose_64x64(d, d);
		yafu_mul_64x64_64x64(d, vt_v0[0], vt_v0_next);

		/* compute e, fold it into vnext and update v'*v0 */

		yafu_mul_64x64_64x64(winv[1], vt_a_v[0], e);

		for (i = 0; i < 64; i++)
			e[i] = e[i] & mask0;

		yafu_mul_Nx64_64x64_acc(v[1], e, vnext, n);

		yafu_transpose_64x64(e, e);
		yafu_mul_64x64_64x64(e, vt_v0[1], e);
		for (i = 0; i < 64; i++)
			vt_v0_next[i] = vt_v0_next[i] ^ e[i];

		/* compute f, fold it in. Montgomery shows that 
		   this is unnecessary (f would be zero) if the 
		   previous value of v had full rank */

		if (mask1 != (uint64)(-1)) {
			yafu_mul_64x64_64x64(vt_a_v[1], winv[1], f);

			for (i = 0; i < 64; i++)
				f[i] = f[i] ^ bitmask[i];

			yafu_mul_64x64_64x64(winv[2], f, f);

			for (i = 0; i < 64; i++)
				f2[i] = ((vt_a2_v[1][i] & mask1) ^ 
					  vt_a_v[1][i]) & mask0;

			yafu_mul_64x64_64x64(f, f2, f);

			yafu_mul_Nx64_64x64_acc(v[2], f, vnext, n);

			yafu_transpose_64x64(f, f);
			yafu_mul_64x64_64x64(f, vt_v0[2], f);
			for (i = 0; i < 64; i++)
				vt_v0_next[i] = vt_v0_next[i] ^ f[i];
		}

		/* update the computed solution 'x' */

		yafu_mul_64x64_64x64(winv[0], vt_v0[0], d);
		yafu_mul_Nx64_64x64_acc(v[0], d, x, n);

		/* rotate all the variables */

		tmp = v[2]; 
		v[2] = v[1]; 
		v[1] = v[0]; 
		v[0] = vnext; 
		vnext = tmp;
		
		tmp = winv[2]; 
		winv[2] = winv[1]; 
		winv[1] = winv[0]; 
		winv[0] = tmp;
		
		tmp = vt_v0[2]; 
		vt_v0[2] = vt_v0[1]; 
		vt_v0[1] = vt_v0[0]; 
		vt_v0[0] = vt_v0_next; 
		vt_v0_next = tmp;
		
		tmp = vt_a_v[1]; vt_a_v[1] = vt_a_v[0]; vt_a_v[0] = tmp;
		
		tmp = vt_a2_v[1]; vt_a2_v[1] = vt_a2_v[0]; vt_a2_v[0] = tmp;

		memcpy(s[1], s[0], 64 * sizeof(uint32));
		mask1 = mask0;
		dim1 = dim0;

		/* possibly print a status update */

		dim_solved += dim0;
		if (report_interval) {
			if (dim_solved >= next_report) {
				time_t curr_time = time(NULL);
				double elapsed = difftime(curr_time, 
							first_time);
				uint32 eta = elapsed * (n - dim_solved) /
						(dim_solved - first_dim_solved);

				fprintf(NULL, "linear algebra completed %u "
					"of %u dimensions (%1.1f%%, ETA "
					"%dh%2dm)    \r",
					dim_solved, n, 100.0 * dim_solved / n,
					eta / 3600, (eta % 3600) / 60);
				next_report = dim_solved + report_interval;
				fflush(stderr);
			}
		}

		/* possibly dump a checkpoint file, check for interrupt.
		   Do not checkpoint or interrupt during the first few 
		   iterations, because the checkpoint file does not store 
		   v0 and we would be unable to restart in this case */

		if (dump_interval && iter >= 4) {
			if (dim_solved >= next_dump ||
			    obj->flags & MSIEVE_FLAG_STOP_SIEVING) {
				next_dump = dim_solved + dump_interval;
				yafu_dump_lanczos_state(obj, x, vt_v0, v, vt_a_v, 
						   vt_a2_v, winv, n, 
						   dim_solved, iter, s, dim1);
			}
			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;
		}
	}

	if (report_interval)
		fprintf(stderr, "\n");

	if (VFLAG > 0)
		printf("lanczos halted after %u iterations (dim = %u)\n", 
					iter, dim_solved);
	logprint(obj->logfile, "lanczos halted after %u iterations (dim = %u)\n", 
					iter, dim_solved);

	/* free unneeded storage */

	free(vnext);
	free(scratch);
	free(v0);
	free(vt_a_v[0]);
	free(vt_a_v[1]);
	free(vt_a2_v[0]);
	free(vt_a2_v[1]);
	free(winv[0]);
	free(winv[1]);
	free(winv[2]);
	free(vt_v0_next);
	free(vt_v0[0]);
	free(vt_v0[1]);
	free(vt_v0[2]);

	/* if a recoverable failure occurred, start everything
	   over again */

	if (dim0 == 0 || (obj->flags & MSIEVE_FLAG_STOP_SIEVING)) {
		free(x);
		free(v[0]);
		free(v[1]);
		free(v[2]);
		if (dim0 == 0)
		{
			printf("linear algebra failed; retrying...\n");
			logprint(obj->logfile,"linear algebra failed; retrying...\n");
			exit(1);
		}
		return NULL;
	}

	/* convert the output of the iteration to an actual
	   collection of nullspace vectors. Begin by multiplying
	   the output from the iteration by B */

	yafu_mul_MxN_Nx64(packed_matrix, x, v[1]);
	yafu_mul_MxN_Nx64(packed_matrix, v[0], v[2]);

	/* if necessary, add in the contribution of the
	   first few rows that were originally in B. We 
	   expect there to be about 64 - POST_LANCZOS_ROWS 
	   bit vectors that are in the nullspace of B and
	   post_lanczos_matrix simultaneously */

	if (post_lanczos_matrix) {
		for (i = 0; i < POST_LANCZOS_ROWS; i++) {
			uint64 accum0 = 0;
			uint64 accum1 = 0;
			uint32 j;
			mask0 = bitmask[i];
			for (j = 0; j < n; j++) {
				if (post_lanczos_matrix[j] & mask0) {
					accum0 ^= x[j];
					accum1 ^= v[0][j];
				}
			}
			v[1][i] ^= accum0;
			v[2][i] ^= accum1;
		}
	}

	*num_deps_found = yafu_combine_cols(n, x, v[0], v[1], v[2]);

	/* verify that these really are linear dependencies of B */

	yafu_mul_MxN_Nx64(packed_matrix, x, v[0]);

	for (i = 0; i < n; i++) {
		if (v[0][i] != 0)
			break;
	}
	if (i < n) {
		printf("lanczos error: dependencies don't work\n");
		logprint(obj->logfile, "lanczos error: dependencies don't work\n");
		exit(-1);
	}
	
	free(v[0]);
	free(v[1]);
	free(v[2]);
	//x isn't freed in this code path... but freeing it causes
	//segfaults sometimes, and I don't want to debug this...
	//free(x);

	if (*num_deps_found == 0)
	{
		printf("lanczos error: only trivial "
				"dependencies found\n");
		logprint(obj->logfile, "lanczos error: only trivial "
				"dependencies found\n");
	}
	else
	{
		if (VFLAG > 0)
			printf("recovered %u nontrivial dependencies\n", 
				*num_deps_found);
		logprint(obj->logfile, "recovered %u nontrivial dependencies\n", 
			*num_deps_found);
	}
	return x;
}

/*-----------------------------------------------------------------------*/
uint64 * qs_block_lanczos(fact_obj_t *obj, uint32 nrows, 
			uint32 num_dense_rows, uint32 ncols, 
			qs_la_col_t *B, uint32 *num_deps_found) {
	
	/* External interface to the linear algebra */

	uint64 *post_lanczos_matrix;
	uint64 *dependencies;
	packed_matrix_t packed_matrix;
	uint32 dump_interval;

	if (ncols <= nrows) {
		printf("matrix must have more columns than rows\n");
		logprint(obj->logfile, "matrix must have more columns than rows\n");
		exit(-1);
	}

	/* optionally remove the densest rows of the matrix, and
	   optionally pack a few more rows into dense format */

	post_lanczos_matrix = yafu_form_post_lanczos_matrix(obj, &nrows,
					&num_dense_rows, ncols, B);
	if (num_dense_rows) {
		if (VFLAG > 0)
			printf("matrix includes %u packed rows\n", 
					num_dense_rows);
		logprint(obj->logfile, "matrix includes %u packed rows\n", 
					num_dense_rows);
	}

	yafu_packed_matrix_init(obj, &packed_matrix, B, nrows, 
			   ncols, num_dense_rows);

	/* set up for writing checkpoint files. This only applies
	   to the largest matrices */

	dump_interval = 0;
	if (ncols > 1000000) {

		if (ncols > 4000000)
			dump_interval = 100000;
		else if (ncols > 2000000)
			dump_interval = 250000;
		else
			dump_interval = 500000;
		obj->flags |= MSIEVE_FLAG_SIEVING_IN_PROGRESS;
	}

	/* solve the matrix */

	do {
		dependencies = yafu_block_lanczos_core(obj, &packed_matrix,
						num_deps_found,
						post_lanczos_matrix,
						dump_interval);

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;

	} while (dependencies == NULL);

	if (dump_interval)
		obj->flags &= ~MSIEVE_FLAG_SIEVING_IN_PROGRESS;

	/* note that the following frees any auxiliary packed
	   matrix structures, and also frees the column entries from
	   the input matrix (whether packed or not) */

	yafu_packed_matrix_free(&packed_matrix);
	free(post_lanczos_matrix);
	return dependencies;
}
