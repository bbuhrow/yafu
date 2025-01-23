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

#ifdef HAVE_CUDA
#include "gpu/lanczos_gpu.h"
#else
#include "lanczos.h"
#endif

#ifdef HAVE_MPI

	/* The following is a functional replacement for
	   MPI_Allreduce(), but for large problems can be
	   configured to switch to a bucket accumulation
	   method that is asymptotically faster when the
	   vector is large and needs to be globally accumulated
	   and redistributed across a large number of nodes.

	   The algorithm uses the bucket strategy from the
	   paper "Global Combine on Mesh Architectures with 
	   Wormhole Routing". The implementation below is
	   based on code kindly contributed by Ilya Popovyan */

#if 1
#define GLOBAL_BREAKOVER 5000
#else
#define GLOBAL_BREAKOVER (uint32)(-1) /* turn off the fancy method */
#endif

#define MPI_NODE_0_START if (obj->mpi_la_row_rank + \
                                obj->mpi_la_col_rank == 0) {

#define MPI_NODE_0_END }

/*------------------------------------------------------------------*/
static void global_xor_async(v_t *send_buf, v_t *recv_buf, 
			uint32 total_size, uint32 num_nodes, 
			uint32 my_id, MPI_Datatype mpi_word, MPI_Comm comm, void *send_buf_in) {
	
	uint32 i, j;
	uint32 m, size, chunk, remainder;
	uint32 next_id, prev_id;
	MPI_Status mpi_status;
	MPI_Request mpi_req;
	v_t *curr_buf;
		
	/* split data */

	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;	
	
	/* we expect a circular topology here */

	next_id = mp_modadd_1(my_id, 1, num_nodes);
	prev_id = mp_modsub_1(my_id, 1, num_nodes);
			
	/* stage 1
	   P_m sends P_{m+1} the m-th chunk of data while receiving 
	   another chunk from P_{m-1}, and does the summation op 
	   on the received chunk and (m-1)-th own chunk */

	m = my_id;
	size = chunk;
	if (my_id == num_nodes - 1)
		size += remainder;

	curr_buf = send_buf;

	for (i = 0; i < num_nodes - 1; i++) {
				
		/* asynchronously send the current chunk */

		MPI_TRY(MPI_Isend(curr_buf + m * chunk, size, 
				mpi_word, next_id, 97, 
				comm, &mpi_req))

		/* switch to the recvbuf after the first send */

		curr_buf = recv_buf;
				
		size = chunk;
		if ((int32)(--m) < 0) {
			m += num_nodes;
			size += remainder;
		}

		/* don't wait for send to finish, start the recv 
		   from the previous node */

		MPI_TRY(MPI_Recv(curr_buf + m * chunk, size,
				mpi_word, prev_id, 97, 
				comm, &mpi_status))

		/* combine the new chunk with our own */

#ifdef HAVE_CUDAAWARE_MPI
		vv_xor_gpu(curr_buf + m * chunk, send_buf + m * chunk, size,
			((gpuvec_t *)send_buf_in)->gpudata);
#else
		for (j = 0; j < size; j++) 
			curr_buf[m * chunk + j] = v_xor(curr_buf[m * chunk + j], send_buf[m * chunk + j]);
#endif
		
		/* now wait for the send to end */

		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
#if defined(HAVE_CUDAAWARE_MPI)
		CUDA_TRY(cuCtxSynchronize())
#endif
	}	
		
	/* stage 2
	   P_m sends P_{m+1} m-th chunk of data, now containing 
	   a full summation of all m-th chunks in the comm,
	   while receiving another chunk from P_{m-1} and 
	   puts it to (m-1)-th own chunk */

	curr_buf = recv_buf + m * chunk;
	for (i = 0; i < num_nodes - 1; i++) {
		
		/* async send to chunk the next proc in circle */

		MPI_TRY(MPI_Isend(curr_buf, size, mpi_word, 
				next_id, 98, comm, &mpi_req))
		
		size = chunk;
		curr_buf -= chunk;
		if (curr_buf < recv_buf) {
			curr_buf += chunk * num_nodes;			
			size += remainder;
		}		
		
		/* don't wait for send to finish, start the recv 
		   from the previous proc in circle, put the new 
		   data just where it should be in recv_buf */

		MPI_TRY(MPI_Recv(curr_buf, size, mpi_word,
				prev_id, 98, comm, &mpi_status))
				
		/* now wait for the send to end */

		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
#if defined(HAVE_CUDAAWARE_MPI)
		CUDA_TRY(cuCtxSynchronize())
#endif
	}
}

/*------------------------------------------------------------------*/
void global_xor(void *send_buf_in, void *recv_buf_in, 
		uint32 total_size, uint32 num_nodes, 
		uint32 my_id, MPI_Datatype mpi_word, MPI_Comm comm) {
	
	v_t *send_buf, *recv_buf;

	if (num_nodes == 1) {
		vv_copy(recv_buf_in, send_buf_in, total_size);
		return;
	}

#ifdef HAVE_CUDA
#ifdef HAVE_CUDAAWARE_MPI
	send_buf = (v_t *)((gpuvec_t *)send_buf_in)->gpu_vec;
	recv_buf = (v_t *)((gpuvec_t *)recv_buf_in)->gpu_vec;
	CUDA_TRY(cuCtxSynchronize())
#else
	send_buf = (v_t *)((gpuvec_t *)send_buf_in)->host_vec;
	recv_buf = (v_t *)((gpuvec_t *)recv_buf_in)->host_vec;
	vv_copyout(send_buf, send_buf_in, total_size);
#endif
#else
	send_buf = (v_t *)send_buf_in;
	recv_buf = (v_t *)recv_buf_in;
#endif

	/* only get fancy for large buffers; even the
	   fancy method is only faster when many nodes 
	   are involved */
	
	if (total_size < GLOBAL_BREAKOVER || num_nodes < 2) {
		MPI_TRY(MPI_Allreduce(send_buf, 
				recv_buf, VWORDS * total_size,
				MPI_LONG_LONG, MPI_BXOR, comm))
	} else {
		global_xor_async(send_buf, recv_buf, 
			total_size, num_nodes, my_id, mpi_word, comm, send_buf_in);
	}
	
#if defined(HAVE_CUDA) && !defined(HAVE_CUDAAWARE_MPI)
	vv_copyin(recv_buf_in, recv_buf, total_size); 	
#endif
}

/*------------------------------------------------------------------*/
void global_chunk_info(uint32 total_size, uint32 num_nodes, 
			uint32 my_id, uint32 *chunk_size, 
			uint32 *chunk_start) 
{
	uint32 chunk, remainder;
	
	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;	
	
	if (chunk_start)
		*chunk_start = my_id*chunk;
	
	if (chunk_size) {
		if (my_id == num_nodes - 1)
			*chunk_size = chunk + remainder;
		else
			*chunk_size = chunk;
	}
}

/*------------------------------------------------------------------*/
void global_xor_scatter(void *send_buf_in, void *recv_buf_in, 
			void *scratch_in, uint32 total_size, 
			uint32 num_nodes, uint32 my_id, 
			MPI_Datatype mpi_word, MPI_Comm comm) {
	
	v_t *send_buf, *recv_buf, *scratch;
	uint32 i, j;
	uint32 m, size, chunk, remainder;
	uint32 next_id, prev_id;
	MPI_Status mpi_status;
	MPI_Request mpi_req;
	int *counts, *offsets;
    
	if (num_nodes == 1) {
		vv_copy(recv_buf_in, send_buf_in, total_size);
		return;
	}

#ifdef HAVE_CUDA
#ifdef HAVE_CUDAAWARE_MPI
	send_buf = (v_t *)((gpuvec_t *)send_buf_in)->gpu_vec;
	recv_buf = (v_t *)((gpuvec_t *)recv_buf_in)->gpu_vec;
	scratch  = (v_t *)((gpuvec_t *)scratch_in)->gpu_vec;
	CUDA_TRY(cuCtxSynchronize())
#else
	send_buf = (v_t *)((gpuvec_t *)send_buf_in)->host_vec;
	recv_buf = (v_t *)((gpuvec_t *)recv_buf_in)->host_vec;
	scratch  = (v_t *)((gpuvec_t *)scratch_in)->host_vec;
	vv_copyout(send_buf, send_buf_in, total_size);
#endif
#else
	send_buf = (v_t *)send_buf_in;
	recv_buf = (v_t *)recv_buf_in;
	scratch = (v_t *)scratch_in;
#endif
    
	/* split data */
    
	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;

#if 0 /* Use standard MPI collectives */

	counts = (int *)malloc(num_nodes * sizeof(int));
	offsets = (int *)malloc(num_nodes * sizeof(int));
	size = chunk;
	if (my_id == num_nodes - 1) size += remainder;

	for (i = 0; i < num_nodes; i++) {
		counts[i] = VWORDS * chunk;
		offsets[i] = i * VWORDS * chunk;
	}
	counts[num_nodes - 1] += VWORDS * remainder;

	MPI_TRY(MPI_Reduce_scatter(send_buf, recv_buf, counts,
		MPI_LONG_LONG, MPI_BXOR, comm))
	/* MPI_TRY(MPI_Reduce(send_buf, scratch, VWORDS * total_size, MPI_LONG_LONG,
               MPI_BXOR, 0, comm))
	MPI_TRY(MPI_Scatterv(scratch, counts, offsets,
                 MPI_LONG_LONG, recv_buf, VWORDS * size,
                 MPI_LONG_LONG, 0, comm)) */

	free(counts);
	free(offsets);
#else

	/* we expect a circular topology here */
    
	next_id = mp_modadd_1(my_id, 1, num_nodes);
	prev_id = mp_modsub_1(my_id, 1, num_nodes);
    
	/* P_m sends P_{m+1} the m-th chunk of data while receiving 
	   another chunk from P_{m-1}, and does the summation op 
	   on the received chunk and (m-1)-th own chunk */
    
	m = prev_id;
	size = chunk;
	if (m == num_nodes - 1)
		size += remainder;
        
	for (i = 0; i < num_nodes - 2; i++) {
        
		/* asynchroniously send the current chunk */
        
		MPI_TRY(MPI_Isend(send_buf + m * chunk, 
			size, mpi_word, next_id, 95, 
                        comm, &mpi_req))
        
		/* switch to the recvbuf after the first send */
                
		size = chunk;
		if ((int32)(--m) < 0) {
			m += num_nodes;
			size += remainder;
		}
        
		/* don't wait for send to finish, start the recv 
		   from the previous node */
        
		MPI_TRY(MPI_Recv(scratch, size,
                         mpi_word, prev_id, 95, 
                         comm, &mpi_status))
        
		/* combine the new chunk with our own */
        
#ifdef HAVE_CUDAAWARE_MPI
		vv_xor_gpu(send_buf + m * chunk, scratch, size,
			((gpuvec_t *)send_buf_in)->gpudata);
#else
		for (j = 0; j < size; j++) 
			send_buf[m * chunk + j] = v_xor(send_buf[m * chunk + j], scratch[j]);
#endif
		
		/* now wait for the send to end */
        
		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
#ifdef HAVE_CUDAAWARE_MPI
		CUDA_TRY(cuCtxSynchronize())
#endif
	}	
    
	/* asynchronously send the current chunk */
    
	MPI_TRY(MPI_Isend(send_buf + m * chunk, 
			size, mpi_word, 
			next_id, 95, comm, &mpi_req))
    
	/* switch to the recvbuf after the first send */
    
	size = chunk;
	if ((int32)(--m) < 0) {
		m += num_nodes;
	        size += remainder;
	}
    
	/* don't wait for send to finish, start the recv 
	   from the previous node */
    
	MPI_TRY(MPI_Recv(recv_buf, size,
                     mpi_word, prev_id, 95, 
                     comm, &mpi_status))
    
	/* combine the new chunk with our own */
    
#ifdef HAVE_CUDAAWARE_MPI
		vv_xor_gpu(recv_buf, send_buf + m * chunk, size,
			((gpuvec_t *)send_buf_in)->gpudata);
#else
		for (j = 0; j < size; j++) 
			recv_buf[j] = v_xor(recv_buf[j], send_buf[m * chunk + j]);
#endif

	/* now wait for the send to end */
    
	MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
#ifdef HAVE_CUDAAWARE_MPI
		CUDA_TRY(cuCtxSynchronize())
#endif
#endif

#if defined(HAVE_CUDA) && !defined(HAVE_CUDAAWARE_MPI)
	vv_copyin(recv_buf_in, recv_buf, size);
#endif
}

/*------------------------------------------------------------------*/
void global_allgather(void *send_buf_in, void *recv_buf_in, 
                        uint32 total_size, uint32 num_nodes, 
                        uint32 my_id, MPI_Datatype mpi_word, MPI_Comm comm) {
	

	uint32 i;
	uint32 size, chunk, remainder;
	uint32 next_id, prev_id;
	MPI_Status mpi_status;
	MPI_Request mpi_req;
	v_t *curr_buf, *send_buf, *recv_buf;
	int *counts, *offsets;

	if (num_nodes == 1) {
                vv_copy(recv_buf_in, send_buf_in, total_size);
                return;
        }
    
	/* split data */
    
	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;	
	
	/* we expect a circular topology here */
    
	next_id = mp_modadd_1(my_id, 1, num_nodes);
	prev_id = mp_modsub_1(my_id, 1, num_nodes);
    
	/* P_m sends P_{m+1} m-th chunk of data, now containing 
	   a full summation of all m-th chunks in the comm,
	   while receiving another chunk from P_{m-1} and 
	   puts it to (m-1)-th own chunk */
    
	size = chunk;
	if (my_id == num_nodes - 1)
		size += remainder;

#ifdef HAVE_CUDA
#ifdef HAVE_CUDAAWARE_MPI
	send_buf = (v_t *)((gpuvec_t *)send_buf_in)->gpu_vec;
	recv_buf = (v_t *)((gpuvec_t *)recv_buf_in)->gpu_vec;
	CUDA_TRY(cuCtxSynchronize())
#else
	send_buf = (v_t *)((gpuvec_t *)send_buf_in)->host_vec;
	recv_buf = (v_t *)((gpuvec_t *)recv_buf_in)->host_vec;
	vv_copyout(send_buf, send_buf_in, size);
#endif
#else
	send_buf = (v_t *)send_buf_in;
	recv_buf = (v_t *)recv_buf_in;
#endif

#if 0 /* Use standard MPI collectives */

	counts = (int *)malloc(num_nodes * sizeof(int));
	offsets = (int *)malloc(num_nodes * sizeof(int));
	for (i = 0; i < num_nodes; i++) {
		counts[i] = VWORDS * chunk;
		offsets[i] = i * VWORDS * chunk;
	}
	counts[num_nodes - 1] += VWORDS * remainder;

	MPI_TRY(MPI_Allgatherv(send_buf, VWORDS * size, MPI_LONG_LONG,
                   recv_buf, counts, offsets, MPI_LONG_LONG, comm))

	free(counts);
	free(offsets);
#else
	curr_buf = recv_buf + my_id * chunk;
    
	/* put own part in place first */
#if defined(HAVE_CUDAAWARE_MPI)
	CUDA_TRY(cuMemcpyDtoD(curr_buf, send_buf, size * sizeof(v_t)))
#else
	memcpy(curr_buf, send_buf, size * sizeof(v_t)); /* working on host vec */
#endif
    
	for (i = 0; i < num_nodes - 1; i++){
		
		/* async send to chunk the next proc in circle */
        
		MPI_TRY(MPI_Isend(curr_buf, size, mpi_word, 
				next_id, 96, comm, &mpi_req))
		
		size = chunk;
		curr_buf -= chunk;
		if (curr_buf < recv_buf) {
			curr_buf += chunk * num_nodes;			
			size += remainder;
		}		
		
		/* don't wait for send to finish, start the recv 
		   from the previous proc in circle, put the new 
		   data just where it should be in recv_buf */
        
		MPI_TRY(MPI_Recv(curr_buf, size, mpi_word,
				prev_id, 96, comm, &mpi_status))
        
		/* now wait for the send to end */
        
		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
#ifdef HAVE_CUDAAWARE_MPI
		CUDA_TRY(cuCtxSynchronize())
#endif
	}
#endif

#if defined(HAVE_CUDA) && !defined(HAVE_CUDAAWARE_MPI)
	vv_copyin(recv_buf_in, recv_buf, total_size); 	
#endif
}

/*-----------------------------------------------------------------------*/
v_t * gather_ncols(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *v_in, void *scratch_in, 
			v_t *out) {

#ifdef HAVE_CUDA 
	v_t *v = (v_t *)((gpuvec_t *)v_in)->host_vec;
	v_t *scratch = (v_t *)((gpuvec_t *)scratch_in)->host_vec;
	vv_copyout(v, v_in, packed_matrix->nsubcols); 
#else
	v_t *v = (v_t *)v_in;
	v_t *scratch = (v_t *)scratch_in;
#endif

	MPI_NODE_0_START
	if (out == NULL)
		out = (v_t *)aligned_malloc(packed_matrix->max_ncols * 
						sizeof(v_t), 64);
	MPI_NODE_0_END

	/* gather v into MPI row 0 */

	MPI_TRY(MPI_Gatherv(v,
			packed_matrix->nsubcols, 
			obj->mpi_word, scratch,
			packed_matrix->subcol_counts,
			packed_matrix->subcol_offsets,
			obj->mpi_word, 0, 
			obj->mpi_la_col_grid))

	/* gather row 0 into the root node */

	if (obj->mpi_la_row_rank == 0) {
		MPI_TRY(MPI_Gatherv(scratch,
				packed_matrix->ncols, 
				obj->mpi_word, out,
				packed_matrix->col_counts,
				packed_matrix->col_offsets,
				obj->mpi_word, 0, 
				obj->mpi_la_row_grid))
	}
	
	return out;
}

/*-----------------------------------------------------------------------*/
v_t * gather_nrows(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *scratch_in, v_t *out) {

#ifdef HAVE_CUDA 
	v_t *scratch = (v_t *)((gpuvec_t *)scratch_in)->host_vec;
	vv_copyout(scratch, scratch_in, packed_matrix->nrows); 
#else
	v_t *scratch = (v_t *)scratch_in;
#endif

	MPI_NODE_0_START
	if (out == NULL)
		out = (v_t *)aligned_malloc(packed_matrix->max_ncols * 
						sizeof(v_t), 64);
	MPI_NODE_0_END

	/* gather column 0 into the root node */

	if (obj->mpi_la_col_rank == 0) {
		MPI_TRY(MPI_Gatherv(scratch,
				packed_matrix->nrows, 
				obj->mpi_word, out,
				packed_matrix->row_counts,
				packed_matrix->row_offsets,
				obj->mpi_word, 0, 
				obj->mpi_la_col_grid))
	}
	
	return out;
}

/*-----------------------------------------------------------------------*/
void scatter_ncols(msieve_obj *obj,
			packed_matrix_t *packed_matrix,
			void *out_in, void *scratch_in, 
			v_t *in) {

	/* push out to the top MPI row */

#ifdef HAVE_CUDA 
	v_t *out = (v_t *)((gpuvec_t *)out_in)->host_vec;
	v_t *scratch = (v_t *)((gpuvec_t *)scratch_in)->host_vec;
#else
	v_t *out = (v_t *)out_in;
	v_t *scratch = (v_t *)scratch_in;
#endif

	if (obj->mpi_la_row_rank == 0)
		MPI_TRY(MPI_Scatterv(in, packed_matrix->col_counts,
				packed_matrix->col_offsets, 
				obj->mpi_word, scratch,
				packed_matrix->ncols,
				obj->mpi_word, 0, 
				obj->mpi_la_row_grid))

	/* push down each column */

	MPI_TRY(MPI_Scatterv(scratch, packed_matrix->subcol_counts,
	       			packed_matrix->subcol_offsets, 
	      			obj->mpi_word, out,
				packed_matrix->nsubcols,
	   			obj->mpi_word, 0, 
       				obj->mpi_la_col_grid))
#ifdef HAVE_CUDA
	vv_copyin(out_in, out, packed_matrix->nsubcols); 
#endif
}

#endif /* HAVE_MPI */
