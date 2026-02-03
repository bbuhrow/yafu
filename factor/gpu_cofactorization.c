#include <stdio.h>
#include <string.h>
#include "gmp.h"
#include "microecm.h"
#include "batch_factor.h"

#ifdef HAVE_CUDA_BATCH_FACTOR

#include "gpu_cofactorization.h"

#define MAX_RESIDUE_WORDS 3

enum test_flags {
	DEFAULT_FLAGS = 0,				/* just a placeholder */
	FLAG_USE_LOGFILE = 0x01,	    /* append log info to a logfile */
	FLAG_LOG_TO_STDOUT = 0x02,		/* print log info to the screen */
	FLAG_STOP_GRACEFULLY = 0x04		/* tell library to stop */
};

// kernel function reference
enum {
	GPU_ECM_VEC = 0,
	GPU_ECM96_VEC,
	NUM_GPU_FUNCTIONS /* must be last */
};

// kernel function name (corresponding to an implemented function 
// in a .cu file)
static const char* gpu_kernel_names[] =
{
	"gbl_ecm",
	"gbl_ecm96",
};

// argument type lists for the kernels
static const gpu_arg_type_list_t gpu_kernel_args[] =
{
	/* ecm */
	{ 9,
		{
		  GPU_ARG_INT32,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_UINT32,
		  GPU_ARG_INT32,
		}
	 },
	/* ecm96 */
	{ 10,
		{
		  GPU_ARG_INT32,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_PTR,
		  GPU_ARG_UINT32,
		  GPU_ARG_UINT32,
		  GPU_ARG_INT32,
		}
	 },
};

#define BIGWORDS32 3

void mpz_to_bignum32(uint32_t* bignum, mpz_t gmp_in)
{
	int i;
	mpz_t t;
	mpz_init(t);
	mpz_set(t, gmp_in);

	for (i = 0; i < BIGWORDS32; i++)
	{
		bignum[i] = mpz_get_ui(t) & 0xffffffff;
		mpz_tdiv_q_2exp(t, t, 32);
	}

	mpz_clear(t);
	return;
}

void bignum32_to_mpz(mpz_t gmp_out, uint32_t* bignum)
{
	int i;

	mpz_set_ui(gmp_out, bignum[BIGWORDS32 - 1]);
	for (i = BIGWORDS32 - 2; i >= 0; i--)
	{
		mpz_mul_2exp(gmp_out, gmp_out, 32);
		mpz_add_ui(gmp_out, gmp_out, bignum[i]);
	}

	return;
}

static const double INV_2_POW_32 = 1.0 / (double)((uint64_t)(1) << 32);
static uint32_t uecm_lcg_rand_32B(uint32_t lower, uint32_t upper, uint64_t* ploc_lcg)
{
	*ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
	return lower + (uint32_t)(
		(double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}

uint32_t multiplicative_neg_inverse32(uint64_t a)
{
	uint32_t res = 2 + a;
	res = res * (2 + a * res);
	res = res * (2 + a * res);
	res = res * (2 + a * res);
	return res * (2 + a * res);
}

// the function we use to go and actually do work
// using the kernels, arguments, and GPU contexts/streams defined above.
uint32_t do_gpu_ecm64(device_thread_ctx_t* t)
{
	uint32_t quit = 0;

	gpu_arg_t gpu_args[GPU_MAX_KERNEL_ARGS];

	gpu_launch_t* launch;

	float elapsed_ms;
		
	int threads_per_block = 128;
	int num_blocks = t->array_sz / threads_per_block +
		((t->array_sz % threads_per_block) > 0);

	printf("commencing gpu 64-bit ecm in mode %d on %d inputs\n", 
		t->mode_2lp, t->array_sz);

	fflush(stdout);

	// copy sigma into device memory
	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_u32_array,
		t->u32_array,
		t->array_sz * sizeof(uint32_t),
		t->stream))

	// copy n into device memory
	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_n_array,
		t->modulus_in,
		t->array_sz * sizeof(uint64_t),
		t->stream))

	CUDA_TRY(cuEventRecord(t->start_event, t->stream))

	int curve = 0;
	int total_factors = 0;
	int i;

	// initialize on cpu
	// compute rho, one, and Rsq
	mpz_t rsq;
	mpz_t zn;
	mpz_init(rsq);
	mpz_init(zn);
	for (i = 0; i < t->array_sz; i++)
	{		
		t->rho[i] = multiplicative_neg_inverse32(t->modulus_in[i]);
		t->one[i] = ((uint64_t)0 - t->modulus_in[i]) % t->modulus_in[i];
		mpz_set_ui(rsq, 1);
		mpz_mul_2exp(rsq, rsq, 128);
		t->rsq[i] = mpz_tdiv_ui(rsq, t->modulus_in[i]);
	}

	// copy init values into device memory
	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rsq_array,
		t->rsq,
		t->array_sz * sizeof(uint64_t),
		t->stream))

	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rho_array,
		t->rho,
		t->array_sz * sizeof(uint32_t),
		t->stream))

	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_one_array,
		t->one,
		t->array_sz * sizeof(uint64_t),
		t->stream))

	// run 64-bit ecm curves
	launch = t->launch + GPU_ECM_VEC;

	int orig_size = t->array_sz;
	int nofactors = 0;

	while ((curve < t->curves_2lp) && (total_factors < orig_size)) {

		gpu_args[0].int32_arg = t->array_sz;
		gpu_args[1].ptr_arg = (void*)(t->gpu_n_array);		// n
		gpu_args[2].ptr_arg = (void*)(t->gpu_rho_array);	// rho
		gpu_args[3].ptr_arg = (void*)(t->gpu_one_array);	// unity
		gpu_args[4].ptr_arg = (void*)(t->gpu_rsq_array);	// Rsq
		gpu_args[5].ptr_arg = (void*)(t->gpu_u32_array);	// sigma
		gpu_args[6].ptr_arg = (void*)(t->gpu_a_array);		// f
		gpu_args[7].uint32_arg = t->b1_2lp;
		gpu_args[8].int32_arg = curve;

		gpu_launch_set(launch, gpu_args);

		// specify the x,y, and z dimensions of the thread blocks
		// that are created for the specified kernel function
		CUDA_TRY(cuFuncSetBlockShape(launch->kernel_func,
			threads_per_block, 1, 1))

		//printf("kernel %s, size <%d,%d>, ", gpu_kernel_names[GPU_ECM_VEC],
		//	num_blocks, threads_per_block); fflush(stdout);

		// launch the kernel with the size we just set and 
		// arguments configured by the gpu_launch_set command.
		CUDA_TRY(cuLaunchGridAsync(launch->kernel_func,
			num_blocks, 1, t->stream))

		// copy factors back to host
		CUDA_TRY(cuMemcpyDtoHAsync(t->a, t->gpu_a_array,
			t->array_sz * sizeof(uint64_t), t->stream))

		// swap factored inputs to the end of the list
		int n = t->array_sz;
		int c = 0;
		for (i = 0; i < n; i++)
		{
			mpz_set_ui(zn, t->modulus_in[i]);
			uint64_t factor = t->a[i];

			if ((factor > 1) && 
				(factor < t->modulus_in[i]))
			{
				mpz_set_ui(rsq, factor);

				int bits1 = mpz_sizeinbase(rsq, 2);
				mpz_tdiv_q(rsq, zn, rsq);
				int bits2 = mpz_sizeinbase(rsq, 2);

				if (t->rb_idx_2lp[i] < t->rb->num_relations)
				{
					//printf("recording factor %u in relation %d of %d\n",
					//	factor, t->rb_idx_2lp[i], t->rb->num_relations);
					if ((bits1 <= t->lpb_2lp) && (bits2 <= t->lpb_2lp))
					{
						// valid factorization, save it.
						if (t->mode_2lp == 0)
						{
							cofactor_t* c = t->rb->relations + t->rb_idx_2lp[i];

							if (t->first_side == 0)
							{
								c->lp_r[0] = factor;
								c->lp_r[1] = t->modulus_in[i] / factor;
							}
							else
							{
								c->lp_a[0] = factor;
								c->lp_a[1] = t->modulus_in[i] / factor;
							}
						}
						else
						{
							// in mode 1 these are a-side LPs whose indices
							// into the RB have been loaded into the r-side array.
							cofactor_t* c = t->rb->relations + t->rb_idx_2lp[i];
							uint8_t success = c->success;

							//printf("3LP cofactor %d has initial factor %u and "
							//	"r-side factors %u,%u (2LP success = %u)\n",
							//	t->rb_idx_2lp[i], c->lp_a[0], c->lp_r[0], 
							//	c->lp_r[1], success);

							// we end up here because the 3LP kernel already
							// found one valid factor and put it in position 0.
							// here we record the final two factors.
							if (t->first_side == 0)
							{
								c->lp_a[1] = factor;
								c->lp_a[2] = t->modulus_in[i] / factor;
							}
							else
							{
								c->lp_r[1] = factor;
								c->lp_r[2] = t->modulus_in[i] / factor;
							}
							c->success |= 0x0f;
						}
						t->num_factors_2lp++;
					}
					else
					{
						if (t->mode_2lp == 0)
						{
							//gmp_printf("2LP invalid factor sizes: %Zx = %x * %Zx\n",
							//	zn, factor, rsq);
							cofactor_t* c = t->rb->relations + t->rb_idx_2lp[i];
							c->success = 0;
						}
					}
				}
				else
				{
					printf("invalid relation index %d (of %d)\n", t->rb_idx_2lp[i],
						t->rb->num_relations);
				}

				// whether the factorization was valid or not, we are done
				// with the modulus after finding this factor.
				// load in a new input from the end of the list.
				// we do this so that the gpu continues to see a
				// continguous list of inputs.
				t->modulus_in[i] = t->modulus_in[n - 1];
				t->rsq[i] = t->rsq[n - 1];
				t->one[i] = t->one[n - 1];
				t->rho[i] = t->rho[n - 1];
				t->a[i] = t->a[n - 1];
				t->rb_idx_2lp[i] = t->rb_idx_2lp[n - 1];

				// shrink the list
				n--;
				c++;

				// visit this index again
				i--;
			}
		}

		int lastfactors = total_factors;
		total_factors += c;
		//printf("curve %d: %d of %d factored, %d of %d overall\n", 
		//	curve, c, t->array_sz, total_factors, orig_size);
		t->array_sz = n;

		if (lastfactors == total_factors)
			nofactors++;

		//if (nofactors > 2)
		//	break;

		num_blocks = t->array_sz / threads_per_block +
			((t->array_sz % threads_per_block) > 0);

		// copy new list of N to the gpu
		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_n_array,
			t->modulus_in,
			t->array_sz * sizeof(uint64_t),
			t->stream))

		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rsq_array,
			t->rsq,
			t->array_sz * sizeof(uint64_t),
			t->stream))

		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rho_array,
			t->rho,
			t->array_sz * sizeof(uint32_t),
			t->stream))

		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_one_array,
			t->one,
			t->array_sz * sizeof(uint64_t),
			t->stream))

		// new curves.  The gpu does 1 at a time.
		curve += 1;

	}
	//printf("\n");

	CUDA_TRY(cuEventRecord(t->end_event, t->stream))
	CUDA_TRY(cuEventSynchronize(t->end_event))
	CUDA_TRY(cuEventElapsedTime(&elapsed_ms,
		t->start_event, t->end_event))

	printf("found %d total factors (%d valid) in %1.4f ms\n", 
		total_factors, t->num_factors_2lp, elapsed_ms);

	mpz_clear(rsq);
	mpz_clear(zn);

	/* we have to synchronize now */
	CUDA_TRY(cuStreamSynchronize(t->stream))

	return quit;
}

uint32_t do_gpu_ecm96(device_thread_ctx_t* t)
{
	uint32_t quit = 0;

	gpu_arg_t gpu_args[GPU_MAX_KERNEL_ARGS];

	gpu_launch_t* launch;

	float elapsed_ms;

	int threads_per_block = 128;
	int num_blocks = t->array_sz / threads_per_block +
		((t->array_sz % threads_per_block) > 0);

	printf("commencing gpu 96-bit ecm on %d inputs\n", t->array_sz);
	fflush(stdout);

	// copy sigma into device memory
	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_u32_array,
		t->u32_array,
		t->array_sz * sizeof(uint32_t),
		t->stream))

	// copy n into device memory
	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_n_array,
		t->modulus96_in,
		t->array_sz * 3 * sizeof(uint32_t),
		t->stream))

	CUDA_TRY(cuEventRecord(t->start_event, t->stream))

	int total_factors = 0;
	int i;

	// initialize on cpu
	// compute rho, one, and Rsq
	mpz_t rsq, zn, zf, zc;
	mpz_init(rsq);
	mpz_init(zn);
	mpz_init(zf);
	mpz_init(zc);

	for (i = 0; i < t->array_sz; i++)
	{
		bignum32_to_mpz(zn, &t->modulus96_in[3 * i]);

		uint32_t n32 = t->modulus96_in[3 * i];
		t->rho[i] = multiplicative_neg_inverse32(n32);

		mpz_set_ui(rsq, 1);
		mpz_mul_2exp(rsq, rsq, 96);
		mpz_sub(rsq, rsq, zn);
		mpz_tdiv_r(rsq, rsq, zn);

		mpz_to_bignum32(&t->one96[3 * i], rsq);

		mpz_set_ui(rsq, 1);
		mpz_mul_2exp(rsq, rsq, 192);
		mpz_tdiv_r(rsq, rsq, zn);

		mpz_to_bignum32(&t->rsq96[3 * i], rsq);
	}

	// copy init values into device memory
	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rsq_array,
		t->rsq96,
		t->array_sz * 3 * sizeof(uint32_t),
		t->stream))

	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rho_array,
		t->rho,
		t->array_sz * sizeof(uint32_t),
		t->stream))

	CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_one_array,
		t->one96,
		t->array_sz * 3 * sizeof(uint32_t),
		t->stream))

	// run curves
	launch = t->launch + GPU_ECM96_VEC;

	uint64_t lcg = 0xbaddecafbaddecafull;

	int num2lp_retest = 0;
	int orig_size = t->array_sz;
	int no_factors = 0;
	int max_no_factors = 8;
	int max_curves = t->curves_3lp;
	int curve = 0;
	while ((curve < max_curves) && (total_factors < orig_size)) {

		gpu_args[0].int32_arg = t->array_sz;
		gpu_args[1].ptr_arg = (void*)(t->gpu_n_array);		// n
		gpu_args[2].ptr_arg = (void*)(t->gpu_rho_array);	// rho
		gpu_args[3].ptr_arg = (void*)(t->gpu_one_array);	// unity
		gpu_args[4].ptr_arg = (void*)(t->gpu_rsq_array);	// Rsq
		gpu_args[5].ptr_arg = (void*)(t->gpu_u32_array);	// sigma
		gpu_args[6].ptr_arg = (void*)(t->gpu_res32_array);	// f
		gpu_args[7].uint32_arg = t->b1_3lp;
		gpu_args[8].uint32_arg = 100 * t->b1_3lp;
		gpu_args[9].int32_arg = curve;

		gpu_launch_set(launch, gpu_args);

		int last_factors = num2lp_retest;

		// specify the x,y, and z dimensions of the thread blocks
		// that are created for the specified kernel function
		CUDA_TRY(cuFuncSetBlockShape(launch->kernel_func,
			threads_per_block, 1, 1))

		//printf("kernel %s, size <%d,%d>, ", gpu_kernel_names[GPU_ECM96_VEC],
		//	num_blocks, threads_per_block);

		// launch the kernel with the size we just set and 
		// arguments configured by the gpu_launch_set command.
		CUDA_TRY(cuLaunchGridAsync(launch->kernel_func,
			num_blocks, 1, t->stream))

		// copy factors back to host
		CUDA_TRY(cuMemcpyDtoHAsync(t->factors96, t->gpu_res32_array,
			t->array_sz * 3 * sizeof(uint32_t), t->stream))

		// swap factored inputs to the end of the list
		int n = t->array_sz;
		int c = 0;
		for (i = 0; i < n; i++)
		{
			bignum32_to_mpz(zf, &t->factors96[3 * i]);
			bignum32_to_mpz(zn, &t->modulus96_in[3 * i]);

			if ((mpz_cmp_ui(zf, 1) > 0) && (mpz_cmp(zf, zn) < 0))
			{
				mpz_tdiv_r(rsq, zn, zf);

				if (mpz_cmp_ui(rsq, 0) == 0)
				{
					int bits1 = mpz_sizeinbase(zf, 2);
					mpz_tdiv_q(zc, zn, zf);
					int bits2 = mpz_sizeinbase(zc, 2);

					if (bits1 <= t->lpb_3lp)
					{
						// the cofactor needs further gpu analysis.
						// build up a list on which we'll do 64-bit
						// factorizations as needed.

						// check cofactor
						uint64_t cofactor = mpz_get_ui(zc);

						if ((bits2 <= 64) && (prp_uecm(cofactor) == 0))
						{
							cofactor_t* c = t->rb->relations + t->rb_idx_3lp[i];

							// record the factor we found
							if (t->first_side == 0)
							{
								c->lp_a[0] = mpz_get_ui(zf);
							}
							else
							{
								c->lp_r[0] = mpz_get_ui(zf);
							}
							// and load the cofactor for further factorization
							t->modulus_in[num2lp_retest] = cofactor;
							t->rb_idx_2lp[num2lp_retest] = t->rb_idx_3lp[i];
							num2lp_retest++;
						}

					}
					else if (bits2 <= t->lpb_3lp)
					{
						// we found either an improbably large prime factor
						// or two smaller factors simultaneously.
						// submit this for further gpu analysis.
						// build up a list on which we'll do 64-bit
						// factorizations as needed.  possible 
						// to maybe also do the prp checks on gpu
						// but these are extremely cheap on cpu.

						// check cofactor
						uint64_t cofactor = mpz_get_ui(zf);

						if ((bits1 <= 64) && (prp_uecm(cofactor) == 0))
						{
							// record the factor we found
							cofactor_t* c = t->rb->relations + t->rb_idx_3lp[i];
							if (t->first_side == 0)
							{
								c->lp_a[0] = mpz_get_ui(zc);
							}
							else
							{
								c->lp_r[0] = mpz_get_ui(zc);
							}

							// and load the cofactor for further factorization
							t->modulus_in[num2lp_retest] = cofactor;
							t->rb_idx_2lp[num2lp_retest] = t->rb_idx_3lp[i];
							num2lp_retest++;
						}

					}

					// whether the factorization was valid or not, we are done
					// applying 3LP kernels to this modulus after factoring it.  
					// load in a new input from the end of the list.
					t->modulus96_in[3 * i + 0] = t->modulus96_in[3 * (n - 1) + 0];
					t->modulus96_in[3 * i + 1] = t->modulus96_in[3 * (n - 1) + 1];
					t->modulus96_in[3 * i + 2] = t->modulus96_in[3 * (n - 1) + 2];
					t->rsq96[3 * i + 0] = t->rsq96[3 * (n - 1) + 0];
					t->rsq96[3 * i + 1] = t->rsq96[3 * (n - 1) + 1];
					t->rsq96[3 * i + 2] = t->rsq96[3 * (n - 1) + 2];
					t->one96[3 * i + 0] = t->one96[3 * (n - 1) + 0];
					t->one96[3 * i + 1] = t->one96[3 * (n - 1) + 1];
					t->one96[3 * i + 2] = t->one96[3 * (n - 1) + 2];
					t->rho[i] = t->rho[n - 1];
					t->factors96[3 * i + 0] = t->factors96[3 * (n - 1) + 0];
					t->factors96[3 * i + 1] = t->factors96[3 * (n - 1) + 1];
					t->factors96[3 * i + 2] = t->factors96[3 * (n - 1) + 2];
					t->rb_idx_3lp[i] = t->rb_idx_3lp[n - 1];

					// shrink the list
					n--;
					c++;

					// visit this index again
					i--;
				}
			}
		}

		if (last_factors == num2lp_retest)
			no_factors++;

		total_factors += c;
		//printf("curve %d: %d of %d factored (%d valid), %d of %d overall\n",
		//	curve, c, t->array_sz, num2lp_retest - last_factors, total_factors, orig_size);
		t->array_sz = n;

		// run at least 40 curves and keep running curves until 
		// we get to max curves or we stop finding valid factors.
		if (no_factors >= max_no_factors)
		{
			curve = max_curves;
			printf("halting after not finding any valid factors in the "
				"last %d curves\n", max_no_factors);
		}
		else if (curve == (max_curves - 1))
		{
			curve = max_curves;
			printf("halting after running the maximum specified %d curves\n", 
				max_curves);
		}

		num_blocks = t->array_sz / threads_per_block +
			((t->array_sz % threads_per_block) > 0);

		// copy new list of N to the gpu
		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_n_array,
			t->modulus96_in,
			t->array_sz * sizeof(uint32_t),
			t->stream))

		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rsq_array,
			t->rsq96,
			t->array_sz * 3 * sizeof(uint32_t),
			t->stream))

		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_rho_array,
			t->rho,
			t->array_sz * sizeof(uint32_t),
			t->stream))

		CUDA_TRY(cuMemcpyHtoDAsync(t->gpu_one_array,
			t->one96,
			t->array_sz * 3 * sizeof(uint32_t),
			t->stream))

		// new curves.  The gpu does 1 at a time.
		curve += 1;

	}

	mpz_clear(rsq);
	mpz_clear(zn);
	mpz_clear(zf);
	mpz_clear(zc);

	CUDA_TRY(cuEventRecord(t->end_event, t->stream))
	CUDA_TRY(cuEventSynchronize(t->end_event))
	CUDA_TRY(cuEventElapsedTime(&elapsed_ms,
		t->start_event, t->end_event))

	printf("found %d total factors of 3LP inputs in %1.4f ms\n", 
		total_factors, elapsed_ms);
	t->array_sz = total_factors;

	/* we have to synchronize now */
	CUDA_TRY(cuStreamSynchronize(t->stream))


	printf("running 2LP kernel on %d 3LP-cofactors\n", num2lp_retest);
	t->array_sz = num2lp_retest;
	t->mode_2lp = 1;
	t->num_factors_2lp = 0;
	do_gpu_ecm64(t);

	printf("found %d valid factors\n", t->num_factors_2lp);
	t->num_factors_3lp = t->num_factors_2lp;

	return quit;
}

uint32_t gpu_cofactorization(device_thread_ctx_t* t)
{
	uint32_t quit = 0;
	int i;
	int j;

	// which list is 2LP?
	// load the residues into the 64-bit input array.
	if (t->first_side == 0)
	{
		// 2LP side has 2 words per input
		printf("setting up gpu to factor %d r-side 2LPs\n", t->numres_r);
		for (i = 0; i < t->numres_r; i++) {
			t->modulus_in[i] = ((uint64_t)t->residues_r_in[i * 2 + 1] << 32) |
				(uint64_t)t->residues_r_in[i * 2 + 0];
		}

		// the 2LP factorization code is agnostic to side, so
		// point it to which side it should be tracking.
		t->array_sz = t->numres_r;
		t->rb_idx_2lp = t->rb_idx_r;		
	}
	else
	{
		// 2LP side has 2 words per input
		printf("setting up gpu to factor %d a-side 2LPs\n", t->numres_a);
		for (i = 0; i < t->numres_a; i++) {
			t->modulus_in[i] = ((uint64_t)t->residues_a_in[i * 2 + 1] << 32) |
				(uint64_t)t->residues_a_in[i * 2 + 0];
		}
		// the 2LP factorization code is agnostic to side, so
		// point it to which side it should be tracking.
		t->array_sz = t->numres_a;
		t->rb_idx_2lp = t->rb_idx_a;
	}

	// try to completely factor the 2LP list.  The last handful
	// of curves don't typically make sense to run on the gpu (only
	// a few inputs left) but it also doesn't take much time, so
	// to be lazy we just finish it all here.
	t->mode_2lp = 0;
	t->num_factors_2lp = 0;
	do_gpu_ecm64(t);

	// sometimes the factors of a 2LP are not correctly sized.
	// when that happens, we can ignore the corresponding 3LP side cofactor.
	j = 0;
	for (i = 0; i < t->rb->num_relations; i++)
	{
		if (t->rb->relations[i].success == 0)
		{
			if (t->first_side == 0)
			{
				// if r-side was first, then we now take from a-residues
				t->numres_a--;
			}
			else
			{
				// if a-side was first, then we now take from r-residues
				t->numres_r--;
			}
		}
		else
		{
			if (t->first_side == 0)
			{
				// if r-side was first, then we now take from a-residues
				t->modulus96_in[3 * j + 0] = t->residues_a_in[3 * i + 0];
				t->modulus96_in[3 * j + 1] = t->residues_a_in[3 * i + 1];
				t->modulus96_in[3 * j + 2] = t->residues_a_in[3 * i + 2];
				t->rb_idx_a[j] = t->rb_idx_a[i];
			}
			else
			{
				// if a-side was first, then we now take from r-residues
				t->modulus96_in[3 * j + 0] = t->residues_r_in[3 * i + 0];
				t->modulus96_in[3 * j + 1] = t->residues_r_in[3 * i + 1];
				t->modulus96_in[3 * j + 2] = t->residues_r_in[3 * i + 2];
				t->rb_idx_r[j] = t->rb_idx_r[i];
			}

			// 2LP factorization was good
			j++;
		}
	}

	printf("ignoring %d 3LP-side cofactors due to invalid 2LP-side factorizations\n", 
		t->rb->num_relations - t->numres_a);
	
	// now run the 3LP kernels
	if (t->first_side == 0)
	{
		// if r-side was first, then we now take from a-residues
		// the 3LP factorization code is agnostic to side, so
		// point it to which side it should be tracking.
		t->array_sz = t->numres_a;
		t->rb_idx_3lp = t->rb_idx_a;
	}
	else
	{
		// if a-side was first, then we now take from r-residues
		// the 3LP factorization code is agnostic to side, so
		// point it to which side it should be tracking.
		t->array_sz = t->numres_r;
		t->rb_idx_3lp = t->rb_idx_r;
	}
	t->num_factors_3lp = 0;
	do_gpu_ecm96(t);

	// any survivors have now survived both (had factors found on both
	// 2LP and 3LP sides). 
	// printf("%d total inputs have correctly sized factors on both sides\n", 
	// 	t->num_factors_3lp);

	// double check the number that have success fully flagged
	t->num_factors_3lp = 0;
	for (i = 0; i < t->rb->num_relations; i++)
	{
		if (t->rb->relations[i].success == 0xff)
		{
			t->num_factors_3lp++;
			t->rb->relations[i].success = 1;
			t->rb->num_success++;
		}
		else
		{
			t->rb->relations[i].success = 0;
		}
	}
	printf("%d relations have been flagged as completely factored\n", t->num_factors_3lp);

	return quit;
}

/*------------------------------------------------------------------------*/
// definitions for ECM types/functions that use cuda_xface
/*------------------------------------------------------------------------*/

device_ctx_t* gpu_device_init(int which_gpu)
{
	gpu_config_t gpu_config;
	gpu_info_t* gpu_info;
	size_t gpu_mem;

	device_ctx_t* d = (device_ctx_t*)xcalloc(1, sizeof(device_ctx_t));

	gpu_init(&gpu_config);
	if (gpu_config.num_gpu == 0) {
		printf("error: no CUDA-enabled GPUs found\n");
		exit(-1);
	}

	d->gpunum = which_gpu;
	d->gpu_info = gpu_info = (gpu_info_t*)xmalloc(sizeof(gpu_info_t));
	memcpy(gpu_info, gpu_config.info + which_gpu,
		sizeof(gpu_info_t));

	printf("using GPU %u (%s)\n", which_gpu, gpu_info->name);
	printf("selected card has CUDA arch %d.%d\n",
		gpu_info->compute_version_major,
		gpu_info->compute_version_minor);
	printf("more GPU info:\n");
	printf("\tmax_grid_size: %d x %d x %d\n", gpu_info->max_grid_size[0],
		gpu_info->max_grid_size[1], gpu_info->max_grid_size[2]);
	printf("\tglobal_mem_size: %zd\n", gpu_info->global_mem_size);
	printf("\tconstant_mem_size: %d\n", gpu_info->constant_mem_size);
	printf("\tmax_threads_per_block: %d\n", gpu_info->max_threads_per_block);
	printf("\tmax_thread_dim: %d x %d x %d\n", gpu_info->max_thread_dim[0],
		gpu_info->max_thread_dim[1], gpu_info->max_thread_dim[2]);
	printf("\tnum_compute_units: %d\n", gpu_info->num_compute_units);
	printf("\tregisters_per_block: %d\n", gpu_info->registers_per_block);
	printf("\tshared_mem_size: %d\n", gpu_info->shared_mem_size);
	printf("\twarp_size: %d\n", gpu_info->warp_size);

	return d;
}

void gpu_dev_free(device_ctx_t* d)
{
	free(d->gpu_info);
	free(d);
}

device_thread_ctx_t* gpu_ctx_init(device_ctx_t* d, relation_batch_t *rb) {

	device_thread_ctx_t* t;

	t = (device_thread_ctx_t*)xcalloc(1, sizeof(device_thread_ctx_t));

	t->dev = d;
	t->rb = rb;

	/* every thread needs its own context; making all
	   threads share the same context causes problems
	   with the sort engine, because apparently it
	   changes the GPU cache size on the fly */
	CUDA_TRY(cuCtxCreate(&t->gpu_context,
		CU_CTX_BLOCKING_SYNC,
		d->gpu_info->device_handle))

	/* load GPU kernels */
	char ptxfile[80];
	if (d->gpu_info->compute_version_major == 2) {
		strcpy(ptxfile, "cuda_ecm20.ptx");
	}
	else if (d->gpu_info->compute_version_major == 3) {
		if (d->gpu_info->compute_version_minor < 5)
			strcpy(ptxfile, "cuda_ecm30.ptx");
		else
			strcpy(ptxfile, "cuda_ecm35.ptx");
	}
	else if (d->gpu_info->compute_version_major >= 9) {
		strcpy(ptxfile, "cuda_ecm90.ptx");
	}
	else if (d->gpu_info->compute_version_major >= 8) {
		strcpy(ptxfile, "cuda_ecm80.ptx");
	}
	else if (d->gpu_info->compute_version_major >= 5) {
		strcpy(ptxfile, "cuda_ecm50.ptx");
	}
	else
	{
		printf("sorry, Nvidia doesn't want to support your card\n");
		exit(-1);
	}

	CUDA_TRY(cuModuleLoad(&t->gpu_module, ptxfile))

	printf("successfully loaded kernel code from %s\n",
		ptxfile);

	t->launch = (gpu_launch_t*)xmalloc(NUM_GPU_FUNCTIONS *
		sizeof(gpu_launch_t));

	int i;
	for (i = 0; i < NUM_GPU_FUNCTIONS; i++) {
		gpu_launch_t* launch = t->launch + i;

		gpu_launch_init(t->gpu_module, gpu_kernel_names[i],
			gpu_kernel_args + i, launch);
	}

	/* threads each send a stream of kernel calls */
	CUDA_TRY(cuStreamCreate(&t->stream, 0))

	// for measuring elapsed time
	CUDA_TRY(cuEventCreate(&t->start_event, CU_EVENT_BLOCKING_SYNC))
	CUDA_TRY(cuEventCreate(&t->end_event, CU_EVENT_BLOCKING_SYNC))

	return t;
}

void gpu_ctx_free(device_thread_ctx_t* d)
{
	CUDA_TRY(cuEventDestroy(d->start_event))
	CUDA_TRY(cuEventDestroy(d->end_event))
	CUDA_TRY(cuStreamDestroy(d->stream))
	free(d->launch);
	CUDA_TRY(cuCtxDestroy(d->gpu_context))
}

/* external entry point */
int do_gpu_cofactorization(device_thread_ctx_t* t, uint64_t *lcg)
{
	uint32_t i;
	relation_batch_t* rb = t->rb;

	// the relation batch tracks all factors and metadata of a relation.
	// the gpu only cares about the large factors.  we need
	// to extract this info and put it in the data structures the gpu
	// code wants.
	t->array_sz = rb->num_relations;
	t->residues_r_in = (uint32_t*)xmalloc(sizeof(uint32_t) * MAX_RESIDUE_WORDS * rb->num_relations);
	t->residues_a_in = (uint32_t*)xmalloc(sizeof(uint32_t) * MAX_RESIDUE_WORDS * rb->num_relations);

	// the rb stores all factors, large and small, in one giant list.
	uint32_t* factors = rb->factors;

	// reference to relation_batch index
	t->rb_idx_r = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);
	t->rb_idx_a = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);

	uint32_t kr = 0;
	uint32_t ka = 0;
	t->numres_r = 0;
	t->numres_a = 0;
	t->first_side = -1;
	int j;
	for (i = 0; i < rb->num_relations; i++)
	{
		// metadata for this relation
		cofactor_t* c = rb->relations + i;

		// advance past the small factors for this relation
		factors += c->num_factors_r;
		factors += c->num_factors_a;

		// initialize success as a successful 2LP factorization.
		// this only won't be true if we have an actual 2LP that
		// fails to factor into correctly sized factors.  Many
		// 2LPs are already completely factored and the relation
		// just needs work on the 3LP side.
		c->success = 0xf0;

		// next in the list is the r-side residue
		if (c->lp_r_num_words > 0)
		{
			if (c->lp_r_num_words == 3)
			{
				if (t->first_side == 0)
				{
					// we've already picked the r-side as first side,
					// meaning that the a-side has 3lps as well.
					printf("cuda cofactorization can't handle 3lp on both sides yet\n");
					exit(1);
				}

				t->lpb_3lp = t->lpbr;		// 3lp's are on the r-side
				t->lpb_2lp = t->lpba;		// 2lp's are on the a-side
				t->first_side = 1;		// so do the a-side first.
			}
			for (j = 0; j < c->lp_r_num_words; j++)
			{
				t->residues_r_in[kr++] = factors[j];
			}
			factors += c->lp_r_num_words;
			t->rb_idx_r[t->numres_r] = i;
			t->numres_r++;
		}

		// and then the a-side
		if (c->lp_a_num_words > 0)
		{
			if (c->lp_a_num_words == 3)
			{
				if (t->first_side == 1)
				{
					// we've already picked the a-side as first side,
					// meaning that the r-side has 3lps as well.
					printf("cuda cofactorization can't handle 3lp on both sides yet\n");
					exit(1);
				}

				t->lpb_3lp = t->lpba;		// 3lp's are on the a-side
				t->lpb_2lp = t->lpbr;		// 2lp's are on the r-side
				t->first_side = 0;		// so do the r-side first.
			}
			for (j = 0; j < c->lp_a_num_words; j++)
			{
				t->residues_a_in[ka++] = factors[j];
			}
			factors += c->lp_a_num_words;
			t->rb_idx_a[t->numres_a] = i;
			t->numres_a++;
		}
	}

	if (t->first_side < 0)
	{
		printf("could not determine first side to factor\n");
		exit(1);
	}

	// determine ECM parameters from LPB sizes.
	if (t->lpb_2lp <= 26)
	{
		t->b1_2lp = 85;
		t->curves_2lp = 32;
	}
	else if (t->lpb_2lp <= 28)
	{
		t->b1_2lp = 125;
		t->curves_2lp = 32;
	}
	else if (t->lpb_2lp <= 30)
	{
		t->b1_2lp = 165;
		t->curves_2lp = 40;
	}
	else // <= 32
	{
		t->b1_2lp = 205;
		t->curves_2lp = 40;
	}

	if (t->lpb_3lp <= 26)
	{
		t->b1_3lp = 85;
		t->curves_3lp = 64;
	}
	else if (t->lpb_3lp <= 28)
	{
		t->b1_3lp = 125;
		t->curves_3lp = 64;
	}
	else if (t->lpb_3lp <= 30)
	{
		t->b1_3lp = 165;
		t->curves_3lp = 80;
	}
	else // <= 32
	{
		t->b1_3lp = 205;
		t->curves_3lp = 80;
	}

	/* set up device arrays */

	// ecm
	CUDA_TRY(cuMemAlloc(&t->gpu_u32_array, sizeof(uint32_t) * t->array_sz))
	CUDA_TRY(cuMemAlloc(&t->gpu_rsq_array, sizeof(uint32_t) * 3 * t->array_sz))
	CUDA_TRY(cuMemAlloc(&t->gpu_a_array, sizeof(uint64_t) * t->array_sz))
	CUDA_TRY(cuMemAlloc(&t->gpu_one_array, sizeof(uint32_t) * 3 * t->array_sz))
	CUDA_TRY(cuMemAlloc(&t->gpu_n_array, sizeof(uint32_t) * 3 * t->array_sz))
	CUDA_TRY(cuMemAlloc(&t->gpu_res32_array, sizeof(uint32_t) * 3 * t->array_sz))

	// common
	CUDA_TRY(cuMemAlloc(&t->gpu_rho_array, sizeof(uint32_t) * t->array_sz))

	// set up host arrays

	// ecm64
	t->u32_array = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);
	t->rsq = (uint64_t*)xmalloc(sizeof(uint64_t) * t->array_sz);
	t->modulus_in = (uint64_t*)xmalloc(sizeof(uint64_t) * t->array_sz);
	t->one = (uint64_t*)xmalloc(sizeof(uint64_t) * t->array_sz);
	t->a = (uint64_t*)xmalloc(sizeof(uint64_t) * t->array_sz);

	// ecm96
	t->rsq96 = (uint32_t*)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);
	t->modulus96_in = (uint32_t*)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);
	t->one96 = (uint32_t*)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);
	t->factors96 = (uint32_t*)xmalloc(sizeof(uint32_t) * 3 * t->array_sz);

	// common
	t->rho = (uint32_t*)xmalloc(sizeof(uint32_t) * t->array_sz);

	// generate a sigma for each input
	for (i = 0; i < t->array_sz; i++) {
		t->u32_array[i] = uecm_lcg_rand_32B(7, 0xffffffff, lcg);
	}

	// copy over the moduli to factor
	// printf("setting up gpu to factor %d r-side 2LPs\n", t->numres_r);
	// // 2LP
	// for (i = 0; i < t->numres_r; i++) {
	// 	t->modulus_in[i] = ((uint64_t)t->residues_r_in[i * 2 + 1] << 32) |
	// 		(uint64_t)t->residues_r_in[i * 2 + 0];
	// }
	// 
	// printf("setting up gpu to factor %d a-side 3LPs\n", t->numres_a);
	// // 3LP
	// for (i = 0; i < t->numres_a; i++) {
	// 	t->modulus96_in[3 * i + 0] = t->residues_a_in[i * 3 + 0];
	// 	t->modulus96_in[3 * i + 1] = t->residues_a_in[i * 3 + 1];
	// 	t->modulus96_in[3 * i + 2] = t->residues_a_in[i * 3 + 2];
	// }

	gpu_cofactorization(t);

	// clean up
	free(t->a);
	free(t->u32_array);
	free(t->modulus_in);
	free(t->rho);
	free(t->rsq);
	free(t->one);

	free(t->residues_r_in);
	free(t->residues_a_in);
	free(t->rb_idx_r);
	free(t->rb_idx_a);

	free(t->modulus96_in);
	free(t->rsq96);
	free(t->one96);
	free(t->factors96);

	CUDA_TRY(cuMemFree(t->gpu_a_array))
	CUDA_TRY(cuMemFree(t->gpu_u32_array))
	CUDA_TRY(cuMemFree(t->gpu_n_array))
	CUDA_TRY(cuMemFree(t->gpu_rho_array))
	CUDA_TRY(cuMemFree(t->gpu_rsq_array))
	CUDA_TRY(cuMemFree(t->gpu_one_array))
	CUDA_TRY(cuMemFree(t->gpu_res32_array))

	return 0;
}



#endif