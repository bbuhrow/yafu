#include <stdio.h>
#include <vector>
#include <device_unaryspmv.cuh>

#include "spmv_engine.h"

#if defined(_WIN32) || defined (_WIN64)
	#define SPMV_ENGINE_DECL __declspec(dllexport)
#else
	#define SPMV_ENGINE_DECL __attribute__((visibility("default")))
#endif

#define CUDA_TRY(func) \
        {                                                               \
                cudaError_t status = func;                              \
                if (status != (cudaError_t) CUDA_SUCCESS) {             \
                        const char * str = cudaGetErrorString(status);  \
                        if (!str)                                       \
                                str = "Unknown";                        \
                        printf("error (%s:%d): %s\n", __FILE__, __LINE__, str);\
                        exit(-1);                                       \
                }                                                       \
        }

typedef unsigned int uint32;

struct spmv_engine
{
        spmv_engine()
          : temp_data(0), temp_size(0)
        {
        }

        ~spmv_engine()
        {
                if (temp_size)
                        CUDA_TRY(cudaFree(temp_data))
        }

        void * temp_data;
        size_t temp_size;
};

__device__ v_t operator+(const v_t& left, const v_t& right) {
	return v_xor(left, right);
};

extern "C"
{

SPMV_ENGINE_DECL void * 
spmv_engine_init(int * vbits)
{
	*vbits = VBITS;
	return new spmv_engine;	
}

SPMV_ENGINE_DECL void 
spmv_engine_free(void *e)
{
	delete (spmv_engine *)e;
}

SPMV_ENGINE_DECL void 
spmv_engine_run(void * e, spmv_data_t * data)
{
	spmv_engine *engine = (spmv_engine *)e;
	size_t temp_size;

	DeviceUnarySpmv::CsrMV(NULL, temp_size, 
		(int *)data->row_entries, (int *)data->col_entries, (v_t *)data->vector_in, (v_t *)data->vector_out,
		data->num_rows, data->num_cols, data->num_col_entries, v_zero);

	if (temp_size > engine->temp_size) {
		if (engine->temp_size) CUDA_TRY(cudaFree(engine->temp_data))
		CUDA_TRY(cudaMalloc(&engine->temp_data, temp_size))
		engine->temp_size = temp_size;
		printf("Allocated %0.1f MB for SpMV library\n", (double)temp_size / 1048576);
	}

	// Run SpMV: y = A x + y
	DeviceUnarySpmv::CsrMV(engine->temp_data, temp_size,
		(int *)data->row_entries, (int *)data->col_entries, (v_t *)data->vector_in, (v_t *)data->vector_out,
		data->num_rows, data->num_cols, data->num_col_entries, v_zero);
}

} // extern "C"
