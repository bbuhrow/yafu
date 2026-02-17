#pragma once
#include <stdint.h>
#include "batch_factor.h"
#include "cuda_xface.h"

#ifdef _MSC_VER
// so I can browse the code in visual studio
#define HAVE_CUDA_BATCH_FACTOR
#endif

#ifdef HAVE_CUDA_BATCH_FACTOR
typedef struct {
	int gpunum;
	gpu_info_t* gpu_info;
} device_ctx_t;

typedef struct {
	// reference to the configured device
	device_ctx_t* dev;

	// config from calling code
	uint32_t lpba;			// lpb on the a-side
	uint32_t lpbr;			// lpb on the r-side
	uint32_t mfba;			// to determine if 3LP or 2LP
	uint32_t mfbr;			// to determine if 3LP or 2LP
	int verbose;

	// internal config
	int lpb_2lp;
	int lpb_3lp;
	int first_side;			// first side is the 2LP side, so we
							// can skip some of the more expensive 3LPs
	uint32_t b1_2lp;		
	uint32_t b2_2lp;
	uint32_t curves_2lp;
	uint32_t b1_3lp;
	uint32_t b2_3lp;
	uint32_t curves_3lp;
	uint32_t stop_nofactor;

	// threads need their own context and stream
	CUcontext gpu_context;
	CUmodule gpu_module;
	CUstream stream;
	gpu_launch_t* launch;
	CUevent start_event;
	CUevent end_event;

	// the data structure to work on
	relation_batch_t* rb;

	// gpu data:
	CUdeviceptr gpu_a_array;
	CUdeviceptr gpu_b_array;
	CUdeviceptr gpu_n_array;
	CUdeviceptr gpu_one_array;
	CUdeviceptr gpu_rsq_array;
	CUdeviceptr gpu_u32_array;
	CUdeviceptr gpu_res32_array;
	CUdeviceptr gpu_rho_array;

	uint32_t array_sz;
	uint32_t* u32_array;
	uint32_t flags;

	// host data
	uint64_t* a;
	uint64_t* b;

	// common
	uint32_t* rho;

	// for ecm64
	uint64_t* rsq;
	uint64_t* one;
	uint64_t* modulus_in;
	uint64_t* sigma_in;

	// for ecm96
	uint32_t* rsq96;
	uint32_t* one96;
	uint32_t* modulus96_in;
	uint32_t* factors96;

	// tracking factors to input residues
	uint32_t* residues_r_in;
	uint32_t* residues_a_in;
	uint32_t numres_r;
	uint32_t numres_a;
	uint32_t* rb_idx_r;
	uint32_t* rb_idx_a;
	uint32_t* rb_idx_2lp;
	uint32_t* rb_idx_3lp;
	uint32_t num_factors_2lp;
	uint32_t num_factors_3lp;

	int mode_2lp;			// 0: 2LP factorizations.  1: 3LP-cofactorizations

} device_thread_ctx_t;

// create and destroy gpu context info
device_ctx_t* gpu_device_init(int which_gpu);
device_thread_ctx_t* gpu_ctx_init(device_ctx_t* d, relation_batch_t* rb);
void gpu_dev_free(device_ctx_t* d);
void gpu_ctx_free(device_thread_ctx_t* t);

// do gpu cofactorization work
int do_gpu_cofactorization(device_thread_ctx_t* t, uint64_t* lcg,
	int b1_3lp_ovr, int b2_3lp_ovr, int b1_2lp_ovr, int b2_2lp_ovr,
	int curves_3lp_ovr, int curves_2lp_ovr);

#endif

