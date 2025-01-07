#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cub/cub.cuh>

#include "sort_engine.h"

typedef unsigned int uint32;

#if defined(_WIN32) || defined (_WIN64)
	#define SORT_ENGINE_DECL __declspec(dllexport)
	typedef unsigned __int64 uint64;
#else
	#define SORT_ENGINE_DECL __attribute__((visibility("default")))
	typedef unsigned long long uint64;
#endif

#define CUDA_TRY(func) \
	{ 			 					\
		cudaError_t status = func;				\
		if (status != (cudaError_t) CUDA_SUCCESS) {		\
			const char * str = cudaGetErrorString(status);	\
			if (!str)					\
				str = "Unknown";			\
			printf("error (%s:%d): %s\n", __FILE__, __LINE__, str);\
			exit(-1);					\
		}							\
	}

struct sort_engine
{
	sort_engine()
	  : temp_data(0), temp_size(0)
	{
	}

	~sort_engine()
	{
		if (temp_size)
			CUDA_TRY(cudaFree(temp_data))
	}

	void * temp_data;
	size_t temp_size;
};

extern "C"
{

SORT_ENGINE_DECL void * 
sort_engine_init(void)
{
	return new sort_engine;
}

SORT_ENGINE_DECL void 
sort_engine_free(void * e)
{
	delete (sort_engine *)e;
}

SORT_ENGINE_DECL void 
sort_engine_run(void * e, sort_data_t * data)
{
	// arrays are assumed packed together; check
	// they would all start on a power-of-two boundary

	if (data->num_arrays > 1 && data->num_elements % 16) {
		printf("sort_engine: invalid array size\n");
		exit(-1);
	}

	sort_engine *engine = (sort_engine *)e;

	if (data->key_bits <= 32) {
		for (size_t i = 0; i < data->num_arrays; i++) {

			cub::DoubleBuffer<uint32> keys(
					(uint32 *)data->keys_in +
						i * data->num_elements,
					(uint32 *)data->keys_in_scratch +
						i * data->num_elements);

			cub::DoubleBuffer<uint32> values(
					(uint32 *)data->data_in +
						i * data->num_elements,
					(uint32 *)data->data_in_scratch +
						i * data->num_elements);

			// allocate device temp space (persistently)

			size_t temp_size;

			cub::DeviceRadixSort::SortPairs(
						0,
						temp_size,
						keys,
						values,
						data->num_elements,
						0,
						data->key_bits,
						data->stream);

			if (temp_size > engine->temp_size) {
				if (engine->temp_size)
					CUDA_TRY(cudaFree(engine->temp_data))

				CUDA_TRY(cudaMalloc(&engine->temp_data, temp_size))
				engine->temp_size = temp_size;
			}

			// sort for real

			cub::DeviceRadixSort::SortPairs(
						engine->temp_data,
						temp_size,
						keys,
						values,
						data->num_elements,
						0,
						data->key_bits,
						data->stream);

			if (keys.selector)
				std::swap(data->keys_in, data->keys_in_scratch);
			if (values.selector)
				std::swap(data->data_in, data->data_in_scratch);
		}
	}
	else {
		for (size_t i = 0; i < data->num_arrays; i++) {

			cub::DoubleBuffer<uint64> keys(
					(uint64 *)data->keys_in +
						i * data->num_elements,
					(uint64 *)data->keys_in_scratch +
						i * data->num_elements);

			cub::DoubleBuffer<uint32> values(
					(uint32 *)data->data_in +
						i * data->num_elements,
					(uint32 *)data->data_in_scratch +
						i * data->num_elements);

			// allocate device temp space (persistently)

			size_t temp_size;

			cub::DeviceRadixSort::SortPairs(
						0,
						temp_size,
						keys,
						values,
						data->num_elements,
						0,
						data->key_bits,
						data->stream);

			if (temp_size > engine->temp_size) {
				if (engine->temp_size)
					CUDA_TRY(cudaFree(engine->temp_data))

				CUDA_TRY(cudaMalloc(&engine->temp_data, temp_size))
				engine->temp_size = temp_size;
			}

			// sort for real

			cub::DeviceRadixSort::SortPairs(
						engine->temp_data,
						temp_size,
						keys,
						values,
						data->num_elements,
						0,
						data->key_bits,
						data->stream);

			if (keys.selector)
				std::swap(data->keys_in, data->keys_in_scratch);
			if (values.selector)
				std::swap(data->data_in, data->data_in_scratch);
		}
	}
}

} // extern "C"
