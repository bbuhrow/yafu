/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"
#include "factor.h"
#include "util.h"

/* below, the routines for building and solving a qs matrix using
jasonp's block lanczos routines are implemented */

static void build_qs_matrix(uint32 ncols, qs_la_col_t *cols, 
		    	siqs_r *relation_list);


/*------------------------------------------------------------------*/
void qs_solve_linear_system(fact_obj_t *obj, uint32 fb_size, 
		    uint64 **bitfield, siqs_r *relation_list, 
		    qs_la_col_t *cycle_list, uint32 *num_cycles) {

	/* Generate linear dependencies among the relations
	   in full_relations and partial_relations */

	qs_la_col_t *cols;
	uint32 nrows, ncols;
	uint64 *dependencies;
	uint32 num_deps;

	ncols = *num_cycles;
	nrows = fb_size;
	cols = cycle_list;

	/* convert the list of relations from the sieving 
	   stage into a matrix. */

	build_qs_matrix(ncols, cols, relation_list);
	count_qs_matrix_nonzero(obj, fb_size, 0, ncols, cols);

	/* reduce the matrix dimensions to ignore almost empty rows */

	reduce_qs_matrix(obj, &nrows, 0, &ncols, cols, NUM_EXTRA_RELATIONS);

	if (ncols == 0) {
		printf("matrix is corrupt; skipping linear algebra\n");
		free(cols);
		*num_cycles = 0;
		return;
	}

	/* solve the linear system */

	dependencies = qs_block_lanczos(obj, nrows, 0, ncols, cols, &num_deps);

	if (num_deps == 0) {
		free(dependencies);
		return;
	}

	*bitfield = dependencies;
	*num_cycles = ncols;
}

/*------------------------------------------------------------------*/
uint32 qs_merge_relations(uint32 *merge_array,
		  uint32 *src1, uint32 n1,
		  uint32 *src2, uint32 n2) {

	/* Given two sorted lists of integers, merge
	   the lists into a single sorted list with
	   duplicate entries removed. If a particular
	   entry occurs an even number of times in the
	   two input lists, don't add it to the final list
	   at all. Returns the number of elements in the
	   resulting list */

	uint32 i1, i2, val1, val2, count1, count2;
	uint32 num_merge;

	i1 = i2 = 0;
	num_merge = 0;

	while (i1 < n1 && i2 < n2) {
		val1 = src1[i1];
		val2 = src2[i2];

		if (val1 < val2) {
			count1 = 0;
			do {
				i1++; count1++;
			} while (i1 < n1 && src1[i1] == val1);

			if (count1 & 1)
				merge_array[num_merge++] = val1;
		}
		else if (val1 > val2) {
			count2 = 0;
			do {
				i2++; count2++;
			} while (i2 < n2 && src2[i2] == val2);

			if (count2 & 1)
				merge_array[num_merge++] = val2;
		}
		else {
			count1 = count2 = 0;
			do {
				i1++; count1++;
			} while (i1 < n1 && src1[i1] == val1);
			do {
				i2++; count2++;
			} while (i2 < n2 && src2[i2] == val2);

			if ( (count1 + count2) & 1 )
				merge_array[num_merge++] = val1;
		}
	}

	if (i2 == n2) {
		src2 = src1;
		i2 = i1;
		n2 = n1;
	}

	while (i2 < n2) {
		count2 = 0; val2 = src2[i2];

		do {
			i2++; count2++;
		} while (i2 < n2 && src2[i2] == val2);

		if (count2 & 1)
			merge_array[num_merge++] = val2;
	}

	return num_merge;
}

/*------------------------------------------------------------------*/
#define MAX_COL_WEIGHT 1000

static void build_qs_matrix(uint32 ncols, qs_la_col_t *cols, 
			   siqs_r *relation_list) {

	/* Convert lists of relations from the sieving stage
	   into a sparse matrix. The matrix is stored by
	   columns, pointed to by 'cols'. The total number 
	   of nonzero entries in the matrix is returned */

	uint32 i, j;
	qs_la_col_t *col;

	/* Cycles are assumed to be sorted in order of increasing
	   number of relations, so that any cycles that are
	   not used would have created the heaviest matrix columns
	   anyway */
	
	for (i = 0; i < ncols; i++) {
		uint32 buf[MAX_COL_WEIGHT];
		uint32 accum[MAX_COL_WEIGHT];
		uint32 weight;

		/* merge each succeeding relation into the accumulated
		   matrix column */

		col = cols + i;

		for (j = weight = 0; j < col->cycle.num_relations; j++) {
			siqs_r *r = &relation_list[col->cycle.list[j]];
			weight = qs_merge_relations(accum, buf, weight,
						r->fb_offsets, r->num_factors);
			memcpy(buf, accum, weight * sizeof(uint32));
		}

		col->weight = weight;
		col->data = (uint32 *)malloc(weight * sizeof(uint32));
		memcpy(col->data, buf, weight * sizeof(uint32));
	}
}


/*
below, a crude block64 gaussian elimination routine is implemented,
integrated with sqrt and factor tracking code.
For small inputs, it is fast enough.

based on discussion in Factorization and Primalty Testing
by David M. Bressoud
*/

static uint64 bitValRead64(uint64 **m, int row, int col);

int BlockGauss_MPQS(mpqs_rlist *full, mpqs_rlist *partial, z *apoly, z *bpoly,
			fb_list_mpqs *fb, z *n, int mul, 
			FILE *sieve_log,z *factors,uint32 *num_factor,uint32 QSflag)
{
	int i,j,k,l,a,q,polynum;
	int *bl;
	uint8 **m;		//matrix of the powers of the prime decompositions of the relations over the factor base
	uint64 **m2_64;	//m mod 2, packed into 32 bit words
	uint64 **aug_64;	//matrix to store the permutations of the rows of m2, packed into 32 bit words
	uint32 largep, bool_val, B = fb->B;
	uint32 *partial_index;
	int num_f,num_p;
	int num_r,num_col,num_col_aug,set_continue;
	const int blocksz = 64;

	uint32 *pd;
	uint32 r;
	z zx, zy, tmp, tmp2, tmp3, tmp4, nn,poly_b,poly_d1,poly_d2,tmp_a;
	
	zInit(&zx);
	zInit(&zy);
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&tmp3);
	zInit(&tmp4);
	zInit(&nn);
	zInit(&poly_b);
	zInit(&poly_d1);
	zInit(&poly_d2);
	zInit(&tmp_a);

	num_f = full->num_r;
	num_p = partial->act_r;

	num_r = full->num_r + partial->act_r;
	num_col = (uint32)((B/blocksz)+1);
	num_col_aug = (uint32)(num_r/blocksz+1);

	//allocate storage based on total number of relations.
	pd = (uint32 *)malloc(B * sizeof(uint32));
	partial_index = (uint32 *)malloc(num_p * sizeof(uint32));

	aug_64 = (uint64 **)malloc(num_r * sizeof(uint64 *));
	for (i=0; i<num_r; i++)
		aug_64[i] = (uint64 *)malloc(num_col_aug * sizeof(uint64));

	m2_64 = (uint64 **)malloc(num_r * sizeof(uint64 *));
	for (i=0; i<num_r; i++)
		m2_64[i] = (uint64 *)malloc(num_col * sizeof(uint64));

	m = (uint8 **)malloc(num_r * sizeof(uint8 *));
	for (i=0; i<num_r; i++)
		m[i] = (uint8 *)malloc(B * sizeof(uint8));

	bl = (int *)malloc(num_r * sizeof(int));

	//write fulls to m
	for (i=0;i<num_f;i++)
	{
		//Initialize
		for (j=0;j<(int)B;j++)
			m[i][j] = 0;

		//copy the pd's of the fboffsets to the correct location in m
		//offset 0 is special - indicates the parity of the offset
		m[i][0] = (uint8)full->list[i]->fboffset[0];
		j=1;
		while (j<full->list[i]->num_factors)
		{
			m[i][full->list[i]->fboffset[j]]++;
			j++;
		}
	}

	//write fulls from partials to m, probably also redundant to do it this way?
	largep = partial->list[0]->largeprime;
	j=num_f;
	for (i=1;i<(int)partial->num_r;i++)
	{
		if (partial->list[i]->largeprime == largep)
		{
			//this partial's largep is the same as the one before, add the pd's and copy to m
			for (k=0;k<(int)B;k++)
				m[j][k] = 0;

			//do the factor of -1
			m[j][0] = (uint8)partial->list[i-1]->fboffset[0];
			//then the rest
			k=1;
			while (k<partial->list[i-1]->num_factors) 
			{
				m[j][partial->list[i-1]->fboffset[k]]++;
				k++;
			}
			//factor of -1
			m[j][0] += partial->list[i]->fboffset[0];
			//the rest
			k=1;
			while (k<partial->list[i]->num_factors)
			{
				m[j][partial->list[i]->fboffset[k]]++;
				k++;
			}
			//remember the index of the partial that made this full relation, we'll need it later
			partial_index[j-num_f]=i;
			//increment the relation counter
			j++;
		}
		largep = partial->list[i]->largeprime;
	}

	//construct the bit matrix
	for (i=0;i<num_r;i++)
	{
		for (j=0;j<num_col;j++)
		{
			m2_64[i][j] = 0;
			for (k=0;k<blocksz;k++)
			{
				if ((blocksz*j+k) < (int)B)
					m2_64[i][j] |= ((uint64)((uint64)m[i][blocksz*j+k]%2) << k);
			}
		}
	}

	//construct augmented matrix
	for (i=0;i<num_r;i++)
	{
		for (j=0;j<num_col_aug;j++)
		{
			aug_64[i][j] = 0;
			for (k=0;k<blocksz;k++)
			{
				if ((blocksz*j+k)==i)
					aug_64[i][j] = ((uint64)(1) << (uint64)(k));
			}
		}
	}

	*num_factor=0;
	if (VFLAG && !(QSflag == 1))
		printf("starting block64 gaussian elimination\n");
	//initialize blacklist
	for (i=0;i<num_r;i++) bl[i] = 0;
	//search over all columns, right to left (more sparse on the right side)
	for (i=B-1;i>=0;i--)
	{
		//and all rows
		for (j=0;j<num_r;j++)
		{
			//if the j'th row, i'th bit is 1 and not blacklisted, continue
			bool_val = (bitValRead64(m2_64,j,i) != 0) && (bl[j] == 0);
			if (bool_val)
			{
				//add the j'th row mod 2 to all rows after it with a 1 in the ith column
				for (k=j+1;k<num_r;k++)
				{
					bool_val = (bitValRead64(m2_64,k,i) != 0) && (bl[k] == 0);
					if (bool_val)
					{
						//found one in the k'th row.  add to the j'th row starting at column i.
						//record the addition in the augmented matrix
						for (l=(uint32)(i/blocksz);l>=0;l--)
							m2_64[k][l] = m2_64[j][l] ^ m2_64[k][l];
						for (l=0;l<num_col_aug;l++)
							aug_64[k][l] = aug_64[k][l] ^ aug_64[j][l];
						
						//then check if the row is all zeros
						a=0;
						for (l=(uint32)(i/blocksz);l>=0;l--)
							a = a || m2_64[k][l];

						if (a==0)
						{
							//initialize solution vector
							for (l=0;l<(int)B;l++) pd[l] = 0;

							//found a potential solution. check it.
							for (l=0;l<num_r;l++)
							{
								bool_val = bitValRead64(aug_64,k,l) != 0;
								if (bool_val)
								{
									//then the l'th row of m was involved
									for (q=0;q<(int)B;q++)
										pd[q] += m[l][q];
								}
							}

							//compute x mod n
							zCopy(&zOne,&zy);
							zCopy(&zOne,&zx);
							for (l=0;l<num_r;l++)
							{
								bool_val = bitValRead64(aug_64,k,l) != 0;
								if (bool_val)
								{
									if (QSflag == 1)
									{
										if (l >= num_f)
										{
											//if l >= num_f, then this row refers to a relation generated from two partials.
											//we'll need to go back to the two partial locations to find the two offsets to
											//multiply together
											//luckily, we've remembered the index in the complete list of partials that 
											//created this full relation
											
											//sqrt_n is stored in poly_a
											//polynum = partial->list[partial_index[l-num_f]]->polynum;
											//zCopy(&partial->list[partial_index[l-num_f]]->poly_a,&tmp_a);

											if (partial->list[partial_index[l-num_f]]->fboffset[0] == 0)
												zShortAdd(&apoly[0],partial->list[partial_index[l-num_f]]->offset,&nn);		//compute Q(x1)
											else
												zShortSub(&apoly[0],partial->list[partial_index[l-num_f]]->offset,&nn);

											//sqrt_n is stored in poly_a
											//zCopy(&partial->list[partial_index[l-num_f]-1]->poly_a,&tmp_a);
											//polynum = partial->list[partial_index[l-num_f]-1]->polynum;

											if (partial->list[partial_index[l-num_f]-1]->fboffset[0] == 0)
												zShortAdd(&apoly[0],partial->list[partial_index[l-num_f]-1]->offset,&tmp3);		//compute Q(x2)
											else
												zShortSub(&apoly[0],partial->list[partial_index[l-num_f]-1]->offset,&tmp3);

											zMul(&nn,&tmp3,&tmp4);	//compute Q(x1)*Q(x2)
											zMul(&zx,&tmp4,&tmp);	//accumulate with previous terms
											zDiv(&tmp,n,&tmp2,&zx);	//mod n
											
											//include the large prime in mp_y
											zShortMul(&zy,partial->list[partial_index[l-num_f]]->largeprime,&tmp2);
											zDiv(&tmp2,n,&tmp3,&zy);
										}
										else
										{
											//polynum = full->list[l]->polynum;
											if (full->list[l]->fboffset[0] == 0)
												zShortAdd(&apoly[0],full->list[l]->offset,&nn);
											else
												zShortSub(&apoly[0],full->list[l]->offset,&nn);

											zMul(&zx,&nn,&tmp);	//accumulate with previous terms
											zDiv(&tmp,n,&tmp2,&zx);	//mod n
										}
									}
									else
									{
										//printf("accumulating relation %d\n",l);
										//then the l'th relation is involved in the product of relations
										if (l >= num_f)
										{
											//if l >= num_f, then this row refers to a relation generated from two partials.
											//we'll need to go back to the two partial locations to find the two offsets to
											//multiply together
											//luckily, we've remembered the index in the complete list of partials that 
											//created this full relation
											
											//our relation is of the form (ax + b)^2 == a(ax^2 + 2bx + c) mod n
											//(ax^2 + 2bx + c) is what we trial divided, and we remembered
											//a and b, so we can form the left hand side easily
											
											if (partial->list[partial_index[l-num_f]]->largeprime != partial->list[partial_index[l-num_f]-1]->largeprime)
												printf("ERROR, large primes not equal\n");

											//recreate poly_b from poly_a (store instead??)
											polynum = partial->list[partial_index[l-num_f]]->polynum;
											zCopy(&bpoly[polynum],&poly_b);
											
											//compute Q1(x)
											zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]]->offset,&tmp);
											zNroot(&apoly[polynum],&poly_d1,2);
											if (partial->list[partial_index[l-num_f]]->parity)
												zSub(&tmp,&poly_b,&tmp2);
											else
												zAdd(&tmp,&poly_b,&tmp2);
											zCopy(&tmp2,&tmp);

											//include 'a'
											zMul(&zy,&poly_d1,&tmp2);
											zDiv(&tmp2,n,&tmp3,&zy);

											//compute Q2(x)
											polynum = partial->list[partial_index[l-num_f]-1]->polynum;
											zCopy(&bpoly[polynum],&poly_b);
											
											zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]-1]->offset,&tmp2);
											zNroot(&apoly[polynum],&poly_d2,2);
											if (partial->list[partial_index[l-num_f]-1]->parity)
												zSub(&tmp2,&poly_b,&tmp3);
											else
												zAdd(&tmp2,&poly_b,&tmp3);
											zCopy(&tmp3,&tmp2);

											//compute Q(x1)*Q(x2)
											zMul(&tmp,&tmp2,&tmp4);	
											zDiv(&tmp4,n,&tmp2,&tmp3);	//mod n
											zMul(&zx,&tmp3,&tmp);	//accumulate with previous terms
											zDiv(&tmp,n,&tmp2,&zx);	//mod n
											//include the large prime in mp_y
											sp2z(partial->list[partial_index[l-num_f]]->largeprime,&tmp);
											zMul(&tmp,&zy,&tmp2);
											zDiv(&tmp2,n,&tmp3,&zy);

											//include 'a'
											zMul(&zy,&poly_d2,&tmp2);	
											zDiv(&tmp2,n,&tmp3,&zy);
										}
										else
										{
											//recreate poly_b from poly_a (store instead??)
											polynum = full->list[l]->polynum;
											zCopy(&bpoly[polynum],&poly_b);
											
											//compute Q1(x)
											zShortMul(&apoly[polynum],full->list[l]->offset,&tmp);
											zNroot(&apoly[polynum],&poly_d1,2);
											if (full->list[l]->parity)
												zSub(&tmp,&poly_b,&nn);
											else
												zAdd(&tmp,&poly_b,&nn);

											zMul(&zx,&nn,&tmp);	//accumulate with previous terms
											zDiv(&tmp,n,&tmp2,&zx);	//mod n
											zMul(&zy,&poly_d1,&tmp2);	//sqrt(a) = d is part of mp_y
											zDiv(&tmp2,n,&tmp3,&zy);
										}
									}
								}
							}

							//printf("attempting to factor\n");
							//compute y mod n
							//ignore the factor of -1 in this operation
							for (l=1;l<(int)B;l++)
							{
								if (pd[l] > 0)
								{
									sp2z(fb->list[l].prime,&tmp);
									//pd tracks the exponents of the smooth factors.  we know they are all even
									//at this point.  we don't want to compute pd^2, so divide by 2.
									sp2z(pd[l]/2,&tmp2);
									zModExp(&tmp,&tmp2,n,&tmp3);
									zMul(&tmp3,&zy,&tmp4);
									zDiv(&tmp4,n,&tmp2,&zy);
								}
							}

							//split this off into a subroutine... also look for all non-trivial factors if one is composite
							//compute gcd(x-y,n)
							zSub(&zx,&zy,&tmp);
							zLEGCD(&tmp,n,&nn);
							
							if ((r = (uint32)zShortDiv(&nn,mul,&tmp)) == 0)
							{
								//mul divides this factor
								zCopy(&tmp,&nn);
							}

							if (!(nn.val[0] & 0x1))
							{
								//if it's not odd (factors of 2 creep in there, for some reason
								//remove the 2's
								while (!(nn.val[0] & 0x1))
									zShiftRight(&nn,&nn,1);
							}
							
							
							if ((zCompare(&nn,&zOne) > 0) && (zCompare(n,&nn) > 0))
							{
								//printf("non-trivial factor found = %s\n",z2decstr(&nn,&gstr1));
								zCopy(&nn,&tmp2);
								r = (uint32)zShortDiv(&tmp2,mul,&tmp3);
								if (r == 0)
									zCopy(&tmp3,&nn);

								if (isPrime(&nn))
								{
									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (zCompare(&nn,&factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									zCopy(&nn,&factors[*num_factor]);
									if (sieve_log != NULL)
										logprint(sieve_log,"prp%d = %s\n",ndigits(&nn),z2decstr(&nn,&gstr1));

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}
									//check if we're done by accumulating all factors and comparing to n
									zCopy(&factors[0],&nn);
									for (l=1;l<(int)*num_factor;l++)
									{
										zCopy(&factors[l],&tmp);
										zMul(&tmp,&nn,&tmp2);
										zCopy(&tmp2,&nn);
									}
									if (zBits(&nn) + 10 >= zBits(n))
									{
										//+ 10 accounts for the multiplier in n
										//found all factors, done
										goto free;
									}
								}

								//check the other factor
								zCopy(n,&tmp);
								zDiv(&tmp,&nn,&tmp2,&tmp3);

								//remove the multiplier if necessary
								zCopy(&tmp2,&tmp);
								r = zShortDiv(&tmp,mul,&tmp3);
								if (r == 0)
									zCopy(&tmp3,&tmp);
								else
									zCopy(&tmp2,&tmp);

								if (isPrime(&tmp))
								{
									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (zCompare(&tmp,&factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									zCopy(&tmp,&factors[*num_factor]);
									if (sieve_log != NULL)
										logprint(sieve_log,"prp%d = %s\n",ndigits(&tmp),z2decstr(&tmp,&gstr1));

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}
									//check if we're done by accumulating all factors and comparing to n
									zCopy(&factors[0],&nn);
									for (l=1;l<(int)*num_factor;l++)
									{
										zCopy(&factors[l],&tmp);
										zMul(&tmp,&nn,&tmp2);
										zCopy(&tmp2,&nn);
									}
									if (zBits(&nn) + 10 >= zBits(n))
									{
										//+ 10 accounts for the multiplier in n
										//found all factors, done
										goto free;
									}
								}
							} //if non-trivial factor
						} //if a == 0
					} //if found in k'th row
				} //add jth row mod 2 to all appropriate rows after it
				//blacklist the j'th row
				bl[j] = 1;
			} //if not blacklisted
		} //for all rows
	} //for all columns

	printf("matrix exhausted\n");
	r = (uint32)zShortDiv(n,mul,&tmp);
	for (i=0;(uint32)i<*num_factor;i++)
	{
		zCopy(&tmp,&nn);
		zDiv(&nn,&factors[i],&tmp,&tmp2);
	}
	if (sieve_log != NULL)
		logprint(sieve_log,"c%d = %s\n",ndigits(&tmp),z2decstr(&tmp,&gstr1));

free:
	free(pd);
	free(partial_index);
	for (i=0; i<num_r; i++)
		free(aug_64[i]);
	free(aug_64);
	for (i=0; i<num_r; i++)
		free(m2_64[i]);
	free(m2_64);
	for (i=0; i<num_r; i++)
		free(m[i]);
	free(m);
	
	/*
	free_lmatrix64(aug_64,0,num_r,0,num_col_aug-1);
	free_lmatrix64(m2_64,0,num_r,0,num_col-1);
	free_cmatrix(m,0,num_r,0,B-1);
	*/
	free(bl);
	zFree(&zx);
	zFree(&zy);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&tmp3);
	zFree(&tmp4);
	zFree(&nn);
	zFree(&poly_b);
	zFree(&poly_d1);
	zFree(&poly_d2);
	zFree(&tmp_a);
	return 0;
}

static uint64 masks64[64] = {0x1,0x2,0x4,0x8,
							0x10,0x20,0x40,0x80,
							0x100,0x200,0x400,0x800,
							0x1000,0x2000,0x4000,0x8000,
							0x10000,0x20000,0x40000,0x80000,
							0x100000,0x200000,0x400000,0x800000,
							0x1000000,0x2000000,0x4000000,0x8000000,
							0x10000000,0x20000000,0x40000000,0x80000000,
							0x100000000,0x200000000,0x400000000,0x800000000,
							0x1000000000,0x2000000000,0x4000000000,0x8000000000,
							0x10000000000,0x20000000000,0x40000000000,0x80000000000,
							0x100000000000,0x200000000000,0x400000000000,0x800000000000,
							0x1000000000000,0x2000000000000,0x4000000000000,0x8000000000000,
							0x10000000000000,0x20000000000000,0x40000000000000,0x80000000000000,
							0x100000000000000,0x200000000000000,0x400000000000000,0x800000000000000,
							0x1000000000000000,0x2000000000000000,0x4000000000000000,0x8000000000000000};

static uint64 bitValRead64(uint64 **m, int row, int col)
{
	//col is the column in 0 to B-1 representation
	//read the bit in the packed 64 bit representation of the appropriate row
	//don't bother to check the bounds of m w.r.t row and col, assume caller knows what it's doing
	//return 0 if bit not set, 1 << bit offset otherwize
	int offset, mcol;
	mcol = col/64;
	offset = (col%64);
	return (m[row][mcol] & masks64[offset]);
}
