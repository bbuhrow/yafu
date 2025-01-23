// Adapted from CUB's device/device_spmv.cuh

/******************************************************************************
 * Copyright (c) 2011, Duane Merrill.  All rights reserved.
 * Copyright (c) 2011-2018, NVIDIA CORPORATION.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/
 
#pragma once

#include <stdio.h>
#include <iterator>
#include <limits>
#include <cub/cub.cuh>
// #include <thrust/system/cuda/detail/core/triple_chevron_launch.h>

#include "dispatch_unaryspmv_orig.cuh"

struct DeviceUnarySpmv
{
    /* y = A * x + y */
    template <
        typename            ValueT>
    CUB_RUNTIME_FUNCTION
    static cudaError_t CsrMV(
        void*               d_temp_storage,                     ///< [in] Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
        size_t&             temp_storage_bytes,                 ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
        const int*                d_row_offsets,                      ///< [in] Pointer to the array of \p m + 1 offsets demarcating the start of every row in \p d_column_indices (with the final entry being equal to \p num_nonzeros)
        const int*                d_column_indices,                   ///< [in] Pointer to the array of \p num_nonzeros column-indices of the corresponding nonzero elements of matrix <b>A</b>.  (Indices are zero-valued.)
        const ValueT*             d_vector_x,                         ///< [in] Pointer to the array of \p num_cols values corresponding to the dense input vector <em>x</em>
        ValueT*             d_vector_y,                         ///< [out] Pointer to the array of \p num_rows values corresponding to the dense output vector <em>y</em>
        int                 num_rows,                           ///< [in] number of rows of matrix <b>A</b>.
        int                 num_cols,                           ///< [in] number of columns of matrix <b>A</b>.
        int                 num_nonzeros,                       ///< [in] number of nonzero elements of matrix <b>A</b>.
        ValueT              zero,                               ///< [in[ reduction identity element
        cudaStream_t        stream                  = 0,        ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
        bool                debug_synchronous       = false)    ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    {
        UnarySpmvParams<ValueT, int> spmv_params;
        spmv_params.d_row_end_offsets    = d_row_offsets + 1;
        spmv_params.d_column_indices     = d_column_indices;
        spmv_params.d_vector_x           = d_vector_x;
        spmv_params.d_vector_y           = d_vector_y;
        spmv_params.num_rows             = num_rows;
        spmv_params.num_cols             = num_cols;
        spmv_params.num_nonzeros         = num_nonzeros;
        spmv_params.zero                 = zero;

        return DispatchUnarySpmv<ValueT, int>::Dispatch(
            d_temp_storage,
            temp_storage_bytes,
            spmv_params,
            stream,
            debug_synchronous);
    }

};
