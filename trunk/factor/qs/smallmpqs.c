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
#include "factor.h"
#include "soe.h"
#include "util.h"
#include "common.h"
#include <gmp.h>
#include "gmp_xface.h"

typedef struct
{
	uint32 prime_and_logp;
	uint32 roots;					//root1 is stored in the lower 16 bits, root2 in the upper 16
} smpqs_sieve_fb;

typedef struct
{
	uint32 small_limit;
	uint32 num_blocks;
	uint32 large_mult;
	double fudge_factor;
	uint32 num_extra_relations;
	uint32 dlp_lower;
	uint32 dlp_upper;
	int in_mem;
	int use_dlp;
} sm_mpqs_params;

//holds all the info for a factor base element
typedef struct
{
	uint32 *prime;
	uint32 *correction;
	uint32 *small_inv;
	uint8 *logprime;
} fb_element_sm_mpqs;

typedef struct
{
	uint32 B;
	uint32 small_B;
	uint32 med_B;
	fb_element_sm_mpqs *list;
} fb_list_sm_mpqs;

typedef struct
{
	mpz_t polyb;				//which polyb this relation uses
	uint32 largeprime;		//large prime in the pd.
	uint32 offset;			//offset specifying Q (the quadratic polynomial)
	uint32 polynum;			//which poly this relation uses
	uint32 parity;			//the sign of the offset (x) 0 is positive, 1 is negative
	uint16 *fboffset;		//offsets of factor base primes dividing Q(offset).  max # of fb primes < 2^16 with this choice
	uint8 num_factors;		//number of factor base factors in the factorization of Q
} sm_mpqs_r;

typedef struct
{
	uint32 num_r;
	uint32 act_r;
	uint32 allocated;
	sm_mpqs_r **list;
} sm_mpqs_rlist;

typedef struct
{
	uint64 poly_a;
	uint64 poly_b;
	mpz_t poly_c;
	uint32 poly_d;
	int poly_d_idp;
	int poly_d_idn;
	int side;
	int use_only_p;
} sm_mpqs_poly;

static void smpqs_sieve_block(uint8 *sieve, smpqs_sieve_fb *fb, uint32 start_prime, 
	uint8 s_init, fb_list_sm_mpqs *fullfb);

static int smpqs_check_relations(uint32 sieve_interval, uint32 blocknum, uint8 sieve[] ,mpz_t n, 
								sm_mpqs_poly *poly, uint8 closnuf,
								smpqs_sieve_fb *fb,fb_list_sm_mpqs *fullfb, sm_mpqs_rlist *full, 
								sm_mpqs_rlist *partial, uint32 cutoff,
								uint8 small_cutoff, uint32 start_prime, uint32 parity, 
								uint32 *num, int numpoly);

static void smpqs_trial_divide_Q(mpz_t Q, smpqs_sieve_fb *fb, sm_mpqs_rlist *full, sm_mpqs_rlist *partial,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, fb_list_sm_mpqs *fullfb, 
						  uint32 cutoff, uint8 small_cutoff, uint32 start_prime, int numpoly, 
						  uint32 parity, uint8 closnuf, uint8 bits);

static void smpqs_save_relation(sm_mpqs_rlist *list, uint32 offset, uint32 largeprime, uint32 num_factors, 
						  uint32 rnum, uint16 *fboffset, int numpoly, uint32 parity);

void smpqs_make_fb_mpqs(fb_list_sm_mpqs *fb, uint32 *modsqrt, mpz_t n);
void smpqs_computeB(sm_mpqs_poly *poly, mpz_t n);
void smpqs_nextD(sm_mpqs_poly *poly, mpz_t n);
void smpqs_computeRoots(sm_mpqs_poly *poly, fb_list_sm_mpqs *fb, uint32 *modsqrt, 
	smpqs_sieve_fb *fbp, smpqs_sieve_fb *fbn, uint32 start_prime);

void smpqs_get_more_primes(sm_mpqs_poly *poly);
uint8 smpqs_choose_multiplier(mpz_t n, uint32 fb_size);
int smpqs_BlockGauss(sm_mpqs_rlist *full, sm_mpqs_rlist *partial, uint64 *apoly, uint64 *bpoly,
			fb_list_sm_mpqs *fb, mpz_t n, int mul, 
			mpz_t *factors, uint32 *num_factor);
int sm_check_relation(mpz_t a, mpz_t b, sm_mpqs_r *r, fb_list_sm_mpqs *fb, mpz_t n);

__inline void sm_zcopy(mpz_t src, mpz_t dest)
{
	mpz_set(dest, src);
}

void sm_get_params(int bits, uint32 *B, uint32 *M, uint32 *BL);
int qcomp_smpqs(const void *x, const void *y);


#if defined(GCC_ASM64X) || defined(__MINGW64__)
	#define SM_SCAN_CLEAN asm volatile("emms");	

	#define SM_SIEVE_SCAN_64_VEC					\
		asm volatile (							\
			"movdqa (%1), %%xmm0   \n\t"		\
			"por 16(%1), %%xmm0    \n\t"		\
			"por 32(%1), %%xmm0    \n\t"		\
			"por 48(%1), %%xmm0    \n\t"		\
			"pmovmskb %%xmm0, %%r11   \n\t"		/* output results to 64 bit register */		\
			"testq %%r11, %%r11 \n\t"			/* AND, and set ZF */ \
			"jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
			"movdqa (%1), %%xmm0   \n\t"		/* else, we had hits, move sections of sieveblock back in */ \
			"movdqa 16(%1), %%xmm1   \n\t"		/* there are 16 bytes in each section */ \
			"movdqa 32(%1), %%xmm2   \n\t"		/* extract high bit masks from each byte */ \
			"movdqa 48(%1), %%xmm3   \n\t"		/* and combine into one 64 bit register */ \
			"pmovmskb %%xmm1, %%r9d   \n\t"		/*  */		\
			"pmovmskb %%xmm3, %%r11d   \n\t"	/*  */		\
			"salq $16, %%r9		\n\t"			/*  */ \
			"pmovmskb %%xmm2, %%r10d   \n\t"	/*  */		\
			"salq $48, %%r11		\n\t"		/*  */ \
			"pmovmskb %%xmm0, %%r8d   \n\t"		/*  */		\
			"salq $32, %%r10		\n\t"		/*  */ \
			"orq	%%r11,%%r9		\n\t"		/*  */ \
			"orq	%%r10,%%r8		\n\t"		/*  */ \
			"xorq	%%r11,%%r11		\n\t"		/* initialize count of set bits */ \
			"orq	%%r9,%%r8		\n\t"		/* r8 now holds 64 byte mask results, in order, from sieveblock */ \
			"xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
			"1:			\n\t"					/* top of bit scan loop */ \
			"bsfq	%%r8,%%rcx		\n\t"		/* put least significant set bit index into rcx */ \
			"addq	%%rcx,%%r10	\n\t"			/* add in the offset of this index */ \
			"movb	%%r10b, (%2, %%r11, 1) \n\t"		/* put the bit index into the output buffer */ \
			"shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
			"incq	%%r11		\n\t"			/* increment the count of set bits */ \
			"shrq	$1, %%r8 \n\t"				/* clear the bit */ \
			"testq	%%r8,%%r8	\n\t"			/* check if there are any more set bits */ \
			"jnz 1b		\n\t"					/* loop if so */ \
			"2:		\n\t"						/*  */ \
			"movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \
			: "=r"(result)						\
			: "r"(sieveblock + j), "r"(buffer)	\
			: "xmm0", "xmm1", "xmm2", "xmm3", "r8", "r9", "r10", "r11", "rcx", "cc", "memory");

#elif defined(GCC_ASM32X) || defined(__MINGW32__)
	#define SM_SCAN_CLEAN asm volatile("emms");	

	#define SM_SIEVE_SCAN_64		\
		asm volatile (							\
			"movdqa (%1), %%xmm0   \n\t"		\
			"orpd 16(%1), %%xmm0    \n\t"		\
			"orpd 32(%1), %%xmm0    \n\t"		\
			"orpd 48(%1), %%xmm0    \n\t"		\
			"pmovmskb %%xmm0, %0   \n\t"		\
			: "=r"(result)						\
			: "r"(sieveblock + j), "0"(result)	\
			: "%xmm0");


#elif defined(MSC_ASM32A)
	#define SM_SCAN_CLEAN ASM_M {emms};

	#define SM_SIEVE_SCAN_64	\
		do	{						\
			uint64 *localblock = sieveblock + j;	\
			ASM_M  {			\
				ASM_M mov edi, localblock			\
				ASM_M movdqa xmm0, XMMWORD PTR [edi]	\
				ASM_M por xmm0, XMMWORD PTR [edi + 16]	\
				ASM_M por xmm0, XMMWORD PTR [edi + 32]	\
				ASM_M por xmm0, XMMWORD PTR [edi + 48]	\
				ASM_M pmovmskb ecx, xmm0			\
				ASM_M mov result, ecx};			\
		} while (0);

#elif defined(_WIN64)
	#define SM_SCAN_CLEAN /*nothing*/

	#define SM_SIEVE_SCAN_64	\
		do	{				  		\
			__m128i local_block;	\
			__m128i local_block2;	\
			__m128i local_block3;	\
			__m128i local_block4;	\
			local_block = _mm_load_si128(sieveblock + j); \
			local_block2 = _mm_load_si128(sieveblock + j + 2); \
			local_block3 = _mm_load_si128(sieveblock + j + 4); \
			local_block = _mm_or_si128(local_block, local_block2); \
			local_block = _mm_or_si128(local_block, local_block3); \
			local_block4 = _mm_load_si128(sieveblock + j + 6); \
			local_block = _mm_or_si128(local_block, local_block4); \
			result = _mm_movemask_epi8(local_block); \
		} while (0);


#else	/* compiler not recognized*/

	#define SM_SCAN_CLEAN /*nothing*/
	#undef SM_SIMD_SIEVE_SCAN
	#undef SM_SIMD_SIEVE_SCAN_VEC

#endif
#define SM_SCAN_MASK 0x8080808080808080ULL

//uint64 total_locs;
//uint64 td_locs;

sm_mpqs_params sm_sieve_params;
#define SM_MAX_SMOOTH_PRIMES 100

void smpqs_make_fb_mpqs(fb_list_sm_mpqs *fb, uint32 *modsqrt, mpz_t n)
{
	//finds the factor base primes, and computes the solutions to the congruence x^2 = N mod p
	//for the QS, these are the starting positions of the sieve relative to the sqrt of N.
	//for the MPQS, additional work using the polynomial coefficents and these congruences 
	//needs to be done to compute the starting positions of the sieve.

	int i;
	uint32 b,j,r,k;
	uint32 prime, root1, root2;
	uint8 logp;
	mpz_t tmpr;
	int bits_n = mpz_sizeinbase(n,2);

	mpz_init(tmpr);

	//the 0th element in the fb is always  2, so start searching with 3
	j=2; i=1;
	while (j<fb->B)
	{
		r = mpz_tdiv_ui(n, (fp_digit)spSOEprimes[i]); //r = (uint32)zShortMod(n,(fp_digit)spSOEprimes[i]);
		if (r == 0)
		{
			//p divides n, which means it divides the multiplier.
			//we can still use it, but it only has one solution to x^2 == n mod p instead
			//of two.  just divide its logprime in half.
			//we also can't find the root using shanks-tonelli, but it will be very small
			//because the multiplier is very small, so just use brute force.
			prime = (uint32)spSOEprimes[i];
			b = mpz_tdiv_ui(n, prime); //b = (uint32)zShortMod(n,(fp_digit)prime);
			k=0;
			while (1)
			{
				if (((k*k) % prime) == b)
					break;
				k++;
			}
			root1 = k;
			root2 = prime - k;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5)/2;

			//fill in factor base
			fb->list->prime[j] = prime;
			modsqrt[j] = root1;
			fb->list->logprime[j] = logp;

			//store a couple things so we can replace single precision
			//mods with shifts and muls during trial division
			if (prime < 256)
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 32) / (uint64)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			else
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 40) / (uint64)prime);
				if (floor(256 * MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			
			j++;
			i++;
			continue;
		}

		mpz_set_ui(tmpr, r);
		b = mpz_kronecker_ui(tmpr, spSOEprimes[i]); //b = jacobi_1((fp_digit)r,(fp_digit)spSOEprimes[i]);
		if (b==1)
		{
			//this prime works
			prime = (uint32)spSOEprimes[i];
			if (bits_n > 80)
			{
				fp_digit kk;
				ShanksTonelli_1((fp_digit)r,(fp_digit)prime,&kk);
				k = (uint32)kk;
			}
			else
			{
				// with small n and small factor bases, its faster to brute force it.
				b = mpz_tdiv_ui(n, prime);
				k=0;
				while (1)
				{
					if (((k*k) % prime) == b)
						break;
					k++;
				}
			}
			root1 = (uint32)k;
			root2 = prime - root1;

			//compute logp
			logp = (uint8)(log((double)prime)/log(2.0) + .5);

			//fill in factor base
			fb->list->prime[j] = prime;
			modsqrt[j] = root1;
			fb->list->logprime[j] = logp;

			//store a couple things so we can replace single precision
			//mods with shifts and muls during trial division
			if (prime < 256)
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 32) / (uint64)prime);
				if (floor(MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			else
			{
				fb->list->small_inv[j] = (uint32)(((uint64)1 << 40) / (uint64)prime);
				if (floor(256 * MP_RADIX / (double)prime + 0.5) ==
								(double)fb->list->small_inv[j]) {
					fb->list->correction[j] = 1;
				}
				else {
					fb->list->correction[j] = 0;
					fb->list->small_inv[j]++;
				}
			}
			
			j++;
		}
		i++;
	}

	mpz_clear(tmpr);
	return;
}

#define SM_NUM_PARAM_ROWS 10
void sm_get_params(int bits, uint32 *B, uint32 *M, uint32 *BL)
{
	int i;
	double scale;

	//parameter table
	//bits, fb primes, lp mulitplier, 64k blocks
	//adjustment in v1.27 - more primes and less blocks for numbers > ~80 digits
	//also different scaling for numbers bigger than 100 digits (constant increase
	//of 20% per line)
	int param_table[SM_NUM_PARAM_ROWS][4] = {
		{50,	30,	30,	1},
		{60,	36,	50,	1},
		{70,	50,	50,	1},
		{80,	80,	50,	1},
		{90,	120,	50,	1},
		{100,	175,	50,	1},
		{110,	275,	50,	1},	
		{120,	375,	60,	1},
		{130,	550,	80,	1},
	};

	//linear interpolation according to bit size to determine
	//factor base bound.  use the closest parameter for lp multiplier
	//and number of blocks.

	*B = 0;
	if (bits <= param_table[0][0])
	{
		scale = (double)bits / (double)param_table[0][0];
		*B = (uint32)(scale * (double)(param_table[0][1]));		
		*M = 50;
		*BL = 1;
	}
	else
	{
		for (i=0;i<SM_NUM_PARAM_ROWS - 1;i++)
		{
			if (bits > param_table[i][0] && bits <= param_table[i+1][0])
			{
				scale = (double)(param_table[i+1][0] - bits) /
					(double)(param_table[i+1][0] - param_table[i][0]);
				*B = param_table[i+1][1] - 
					(uint32)(scale * (double)(param_table[i+1][1] - param_table[i][1]));
				
				//sconf->large_mult = (uint32)((double)param_table[i+1][2] - 
				//	(scale * (double)(param_table[i+1][2] - param_table[i][2])) + 0.5);
				*M = (uint32)((param_table[i+1][2] + param_table[i][2])/2.0 + 0.5);
				//sconf->num_blocks = (uint32)((double)param_table[i+1][3] - 
				//	(scale * (double)(param_table[i+1][3] - param_table[i][3])) + 0.5);
				*BL = (uint32)((param_table[i+1][3] + param_table[i][3])/2.0 + 0.5);
			}
		}
	}

	if (*B == 0)
	{
		//off the end of the table, extrapolate based on the slope of 
		//the last two

		scale = (double)(param_table[SM_NUM_PARAM_ROWS-1][1] - param_table[SM_NUM_PARAM_ROWS-2][1]) /
			(double)(param_table[SM_NUM_PARAM_ROWS-1][0] - param_table[SM_NUM_PARAM_ROWS-2][0]);
		*B = (uint32)(((double)bits - param_table[SM_NUM_PARAM_ROWS-1][0]) * 
			scale + param_table[SM_NUM_PARAM_ROWS-1][1]);
		*M = param_table[SM_NUM_PARAM_ROWS-1][2];	//reuse last one

		scale = (double)(param_table[SM_NUM_PARAM_ROWS-1][3] - param_table[SM_NUM_PARAM_ROWS-2][3]) /
			(double)(param_table[SM_NUM_PARAM_ROWS-1][0] - param_table[SM_NUM_PARAM_ROWS-2][0]);
		*BL = (uint32)(((double)bits - param_table[SM_NUM_PARAM_ROWS-1][0]) * 
			scale + param_table[SM_NUM_PARAM_ROWS-1][3]);
	}

	// minimum factor base - for use with really small inputs.
	// not efficient, but needed for decent poly selection
	if (*B < 25)
		*B = 25;

	return;
}


void smallmpqs(fact_obj_t *fobj)
{
	//input expected in fobj->qs_obj.gmp_n
	sm_mpqs_rlist *full, *partial;
	fb_list_sm_mpqs *fb;
	smpqs_sieve_fb *fb_sieve_p,*fb_sieve_n;
	sm_mpqs_poly *poly;
	uint32 *modsqrt;

	mpz_t n;
	mpz_t *factors, tmp, tmp2, tmp3, sqrt_n;
	uint64 *apoly, *bpoly;

	double t_time;
	struct timeval tstart, tend;
	TIME_DIFF *	difference;
	uint32 numpoly, polyalloc;
	uint32 mul,i,j, j2;
	uint32 pmax;							//largest prime in factor base
	uint32 cutoff;
	uint32 sieve_interval;
	uint32 start_prime;
	uint32 num, max_f;
	int digits_n, bits_n, charcount, pindex;
	uint32 num_factors;
	uint8 *sieve;							//sieve values
	uint8 s_init;							//initial sieve value
	uint8 closnuf, small_bits, max_bits;

	//total_locs = 0;
	//td_locs = 0;

	mpz_init(n);
	mpz_init(tmp);
	
	//copy to local variable
	mpz_set(n, fobj->qs_obj.gmp_n);

	if (mpz_cmp_ui(n,1) == 0)
	{
		mpz_clear(n);
		mpz_clear(tmp);
		return;
	}

	if (mpz_even_p(n))
	{
		mpz_clear(n);
		mpz_clear(tmp);
		gmp_printf("%Zu is not odd in smallmpqs\n",n);
		return;
	}

	// factor base bound
	bits_n = mpz_sizeinbase(n,2);

	if (bits_n > 130)
	{
		printf("input too big\n");
		mpz_clear(n);
		mpz_clear(tmp);
		return;
	}

	// don't use the logfile if we see this special flag
	if (fobj->qs_obj.flags != 12345)
	{
		if (fobj->logfile != NULL)
			logprint(fobj->logfile, "starting smallmpqs on C%d: %s\n",
				gmp_base10(n), mpz_conv2str(&gstr1.s, 10, n));
	}

	if (bits_n < 60)
	{
		mpz_t ztmp;
		mpz_init(ztmp);

		j = sp_shanks_loop(n, fobj);	
		if (j > 1)
		{
			mpz_set_64(ztmp, j);
			add_to_factor_list(fobj, ztmp);

			if (fobj->qs_obj.flags != 12345)
			{
				if (fobj->logfile != NULL)
					logprint(fobj->logfile,"prp%d = %u\n",ndigits_1(j), j);
			}

			mpz_tdiv_q_ui(n,n,j);
			add_to_factor_list(fobj, n);

			if (fobj->qs_obj.flags != 12345)
			{
				if (fobj->logfile != NULL)
					logprint(fobj->logfile,
						"prp%d = %s\n",gmp_base10(n),mpz_conv2str(&gstr1.s, 10, n));
			}

			mpz_set_ui(fobj->qs_obj.gmp_n, 1);

			mpz_clear(n);
			mpz_clear(tmp);
			return;
		}
	}

	//empircal tuning of sieve interval based on digits in n
	sm_get_params(bits_n,&j,&sm_sieve_params.large_mult,&sm_sieve_params.num_blocks);

	gettimeofday(&tstart, NULL);
	mpz_init(tmp2);
	mpz_init(tmp3);
	
	//default mpqs parameters
	sm_sieve_params.fudge_factor = 1.3;
	if (fobj->qs_obj.gbl_override_lpmult_flag != 0)
		sm_sieve_params.large_mult = fobj->qs_obj.gbl_override_lpmult;

	sm_sieve_params.num_extra_relations = 32;

	//allocate the space for the factor base
	fb = (fb_list_sm_mpqs *)malloc(sizeof(fb_list_sm_mpqs));
	
	//set fb size from above
	if (fobj->qs_obj.gbl_override_B_flag != 0)
		fb->B = fobj->qs_obj.gbl_override_B;
	else
		fb->B = j;

	if (fobj->qs_obj.gbl_override_blocks_flag != 0)
		sm_sieve_params.num_blocks = fobj->qs_obj.gbl_override_blocks;

	//compute the number of digits in n 
	digits_n = gmp_base10(n);

	//allocate storage for relations based on the factor base size
	max_f = fb->B + 3*sm_sieve_params.num_extra_relations;	
	full = (sm_mpqs_rlist *)malloc((size_t)(sizeof(sm_mpqs_rlist)));
	full->allocated = max_f;
	full->num_r = 0;
	full->act_r = 0;
	full->list = (sm_mpqs_r **)malloc((size_t) (max_f * sizeof(sm_mpqs_r *)));

	//we will typically also generate max_f/2 * 10 partials (empirically determined)
	partial = (sm_mpqs_rlist *)malloc((size_t)(sizeof(sm_mpqs_rlist)));
	partial->allocated = 10*fb->B;
	partial->num_r = 0;
	partial->act_r = 0;
	partial->list = (sm_mpqs_r **)malloc((size_t) (10*fb->B* sizeof(sm_mpqs_r *)));

	//set the sieve interval.  this depends on the size of n, but for now, just fix it.  as more data
	//is gathered, use some sort of table lookup.
	sieve_interval = 32768*sm_sieve_params.num_blocks;

	//allocate the space for the factor base
	modsqrt = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list = (fb_element_sm_mpqs *)malloc((size_t)(sizeof(fb_element_sm_mpqs)));
	fb->list->correction = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list->prime = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list->small_inv = (uint32 *)malloc(fb->B * sizeof(uint32));
	fb->list->logprime = (uint8 *)malloc(fb->B * sizeof(uint8));
	fb_sieve_p = (smpqs_sieve_fb *)malloc((size_t)(fb->B * sizeof(smpqs_sieve_fb)));
	fb_sieve_n = (smpqs_sieve_fb *)malloc((size_t)(fb->B * sizeof(smpqs_sieve_fb)));
	
	//allocate the sieve
	sieve = (uint8 *)xmalloc_align(32768 * sizeof(uint8));

	//allocate the current polynomial
	poly = (sm_mpqs_poly *)malloc(sizeof(sm_mpqs_poly));
	mpz_init(poly->poly_c);

	//allocate the polynomial lists
	polyalloc = 32;
	apoly = (uint64 *)malloc(polyalloc * sizeof(uint64));
	bpoly = (uint64 *)malloc(polyalloc * sizeof(uint64));

	//find multiplier
	mul = (uint32)smpqs_choose_multiplier(n,fb->B);
	mpz_mul_ui(n,n,mul);

	//find new sqrt_n
	mpz_init(sqrt_n);
	mpz_sqrt(sqrt_n,n);

	//find upper bound of Q values
	mpz_tdiv_q_2exp(tmp,n,1); //zShiftRight(&tmp,n,1);
	mpz_sqrt(tmp2,tmp); //zNroot(&tmp,&tmp2,2);
	mpz_mul_ui(tmp, tmp2, sieve_interval); //zShortMul(&tmp2,sieve_interval,&tmp);
	max_bits = mpz_sizeinbase(tmp,2); //(uint8)zBits(&tmp);
	
	//compute the first polynominal 'a' value.  we'll need it before creating the factor base in order
	//to find the first roots
	//'a' values should be as close as possible to sqrt(2n)/M, they should be a quadratic residue mod N (d/N) = 1,
	//and be a prime congruent to 3 mod 4.  this last requirement is so that b values can be computed without using the 
	//shanks-tonelli algorithm, and instead use faster methods.
	//since a = d^2, find a d value near to sqrt(sqrt(2n)/M)
	mpz_mul_2exp(tmp, n, 1); //zShiftLeft(&tmp,n,1);
	mpz_sqrt(tmp2, tmp); //zNroot(&tmp,&tmp2,2);
	mpz_tdiv_q_ui(tmp2, tmp2, sieve_interval); //zShortDiv(&tmp2,sieve_interval,&tmp2);

	mpz_sqrt(tmp, tmp2); //zNroot(&tmp2,&poly->poly_d,2);
	if (mpz_even_p(tmp))
		mpz_add_ui(tmp, tmp, 1); //zAdd(&poly->poly_d,&zOne,&poly->poly_d);
	
	mpz_nextprime(tmp, tmp); //zNextPrime_1(polyd, &fpt, &tmp, 1);
	poly->poly_d = (uint64)mpz_get_ui(tmp); //.val[0];

	if (spSOEprimes[szSOEp - 1] <= poly->poly_d)
		smpqs_get_more_primes(poly);

	pindex = bin_search_uint32(NUM_P, 0, poly->poly_d, spSOEprimes);
	if (pindex < 0)
	{
		printf("prime not found in binary search\n");
		exit(1);
	}
	poly->poly_d_idp = pindex;
	poly->poly_d_idn = pindex;
	poly->side = 1;
	poly->use_only_p = 0;

	smpqs_nextD(poly,n);
	smpqs_computeB(poly,n);	

	fb->list->prime[0] = 1;
	fb->list->prime[1] = 2;

	//construct the factor base, and copy to the sieve factor base
	smpqs_make_fb_mpqs(fb,modsqrt,n);

	for (i=2;i<fb->B;i++)
	{
		fb_sieve_p[i].prime_and_logp = (fb->list->prime[i] << 16) | (fb->list->logprime[i]);
		fb_sieve_n[i].prime_and_logp = (fb->list->prime[i] << 16) | (fb->list->logprime[i]);
	}

	//find the root locations of the factor base primes for this poly
	smpqs_computeRoots(poly,fb,modsqrt,fb_sieve_p,fb_sieve_n,2);

	pmax = fb->list->prime[fb->B-1];
	cutoff = pmax * sm_sieve_params.large_mult;

	//compute the number of bits in M/2*sqrt(N/2), the approximate value
	//of residues in the sieve interval
	//sieve locations greater than this are worthy of trial dividing
	closnuf = (uint8)(double)((bits_n - 1)/2);
	closnuf += (uint8)(log((double)sieve_interval/2)/log(2.0));
	closnuf -= (uint8)(sm_sieve_params.fudge_factor * log(cutoff) / log(2.0));
	
	closnuf += 6;

	//small prime variation -- hand tuned small_bits correction (likely could be better)
	small_bits = 7;
	closnuf -= small_bits;
	start_prime = 7;
	
	s_init = closnuf;

	//print some info to the screen and the log file
	if (VFLAG > 0)
	{
		gmp_printf("n = %Zd (%d digits and %d bits)\n",n,digits_n,bits_n);
		printf("==== sieve params ====\n");
		printf("factor base: %d primes (max prime = %u)\n",fb->B,pmax);
		printf("large prime cutoff: %u (%d * pmax)\n",cutoff,sm_sieve_params.large_mult);
		printf("sieve interval: %d blocks of size %d\n",sieve_interval/32768,32768);
		printf("multiplier is %u\n",mul);
		printf("trial factoring cutoff at %d bits\n",closnuf);
		printf("==== sieving in progress ====\n");
	}

	numpoly = 0;
	num = 0;
	charcount=0;
	while (1)
	{		
		//copy current poly into the poly lists
		if (numpoly < polyalloc)
		{
			apoly[numpoly] = poly->poly_a;
			bpoly[numpoly] = poly->poly_b;
		}
		else
		{
			// get more space for the polys, if needed
			polyalloc *= 2;
			apoly = (uint64 *)realloc(apoly, polyalloc * sizeof(uint64));
			bpoly = (uint64 *)realloc(bpoly, polyalloc * sizeof(uint64));

			apoly[numpoly] = poly->poly_a;
			bpoly[numpoly] = poly->poly_b;
		}

		for (j2=0; j2 < sm_sieve_params.num_blocks; j2++)
		{
			smpqs_sieve_block(sieve,fb_sieve_p,start_prime,s_init,fb);

			i = smpqs_check_relations(sieve_interval,j2,sieve,n,poly,
				s_init,fb_sieve_p,fb,full,partial,cutoff,small_bits,
				start_prime,0,&num,numpoly);

			smpqs_sieve_block(sieve,fb_sieve_n,start_prime,s_init,fb);
		
			i = smpqs_check_relations(sieve_interval,j2,sieve,n,poly,
				s_init,fb_sieve_n,fb,full,partial,cutoff,small_bits,
				start_prime,1,&num,numpoly);
		
			if (partial->num_r > 0)
			{
				//check the partials for full relations
				qsort(partial->list,partial->num_r,sizeof(sm_mpqs_r *),&qcomp_smpqs);
				j=0;
				for (i=0;i<partial->num_r-1;i++)
				{
					if (partial->list[i]->largeprime == partial->list[i+1]->largeprime)
						j++;
				}
				partial->act_r = j;
				
				if (j+(full->num_r) >= fb->B + sm_sieve_params.num_extra_relations) 
				{
					//we've got enough total relations to stop
					goto done;
				}
			}
		}

		//next polynomial
		smpqs_nextD(poly,n);
		smpqs_computeB(poly,n);
		smpqs_computeRoots(poly,fb,modsqrt,fb_sieve_p,fb_sieve_n,2);

		numpoly++;
	}

done:

	if (VFLAG > 0)
		printf("%d relations found: %d full + %d from %d partial, using %d polys\n",
			partial->act_r+full->num_r,full->num_r,partial->act_r,partial->num_r,numpoly);

	gettimeofday (&tend, NULL);
	difference = my_difftime (&tstart, &tend);

	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);

	if (VFLAG > 0)
		printf("QS elapsed time = %6.4f seconds.\n",t_time);

	//printf("%" PRIu64 " blocks scanned, %" PRIu64 " hit\n",total_locs, td_locs);

	//can free sieving structures now
	align_free(sieve);
	free(fb_sieve_p);
	free(fb_sieve_n);
	mpz_clear(poly->poly_c);
	free(poly);
	free(modsqrt);

	num_factors=0;
	factors = (mpz_t *)malloc(MAX_FACTORS * sizeof(mpz_t));
	for (i=0;i<MAX_FACTORS;i++)
		mpz_init(factors[i]);

	gettimeofday(&tstart,NULL);
	i = smpqs_BlockGauss(full,partial,apoly,bpoly,fb,n,mul,
		factors,&num_factors);

	gettimeofday (&tend, NULL);
	difference = my_difftime (&tstart, &tend);

	t_time = ((double)difference->secs + (double)difference->usecs / 1000000);
	free(difference);
	
	if (VFLAG > 0)
		printf("Gauss elapsed time = %6.4f seconds.\n",t_time);
		
	for(i=0;i<num_factors;i++)
	{
		add_to_factor_list(fobj, factors[i]);

		if (fobj->qs_obj.flags != 12345)
		{
			if (fobj->logfile != NULL)
				logprint(fobj->logfile,
					"prp%d = %s\n", gmp_base10(factors[i]),
					mpz_conv2str(&gstr1.s, 10, factors[i]));
		}
		
		mpz_tdiv_q(fobj->qs_obj.gmp_n, fobj->qs_obj.gmp_n, factors[i]);
	}

	mpz_clear(n);
	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
	mpz_clear(sqrt_n);

	for (i=0;i<MAX_FACTORS;i++)
		mpz_clear(factors[i]);
	free(factors);

	for (i=0;i<full->num_r;i++)
	{
		free(full->list[i]->fboffset);
		free(full->list[i]);
	}
	free(full->list);
	free(full);

	for (i=0;i<partial->num_r;i++)
	{
		free(partial->list[i]->fboffset);
		free(partial->list[i]);
	}
	free(partial->list);
	free(partial);

	free(apoly);
	free(bpoly);

	free(fb->list->correction);
	free(fb->list->prime);
	free(fb->list->logprime);
	free(fb->list->small_inv);
	free(fb->list);
	free(fb);

	return;
}

void smpqs_get_more_primes(sm_mpqs_poly *poly)
{
	uint64 num_p;
	int i;

	if (VFLAG > 1)
		printf("smallmpqs getting more primes: poly_d = %u\n",
		poly->poly_d);

	PRIMES = GetPRIMESRange(spSOEprimes, szSOEp, NULL, 0, 
		(uint64)((double)poly->poly_d * 1.25), &num_p);

	//save a batch of sieve primes too.
	spSOEprimes = (uint32 *)realloc(spSOEprimes, 
		(size_t) (num_p * sizeof(uint32)));

	for (i=0;i<num_p;i++)
		spSOEprimes[i] = (uint32)PRIMES[i];

	szSOEp = num_p;
	NUM_P = num_p;
	P_MIN = 0; 
	P_MAX = PRIMES[(uint32)NUM_P-1];

	if (VFLAG > 1)
		printf("prime finding complete, cached %u primes. pmax = %u\n",
		(uint32)NUM_P, (uint32)P_MAX);

	return;
}

static uint8 smpqs_mult_list[] =
	{1, 2, 3, 5, 7, 10, 11, 13, 15, 17, 19, 
	 23, 26, 29, 30, 31, 43, 59, 67, 73};

uint8 smpqs_choose_multiplier(mpz_t n, uint32 fb_size) 
{
	uint32 i, j;
	uint32 num_primes = MIN(2 * fb_size, 30);
	double best_score;
	uint8 best_mult;
	double scores[20];
	uint32 num_multipliers;
	mpz_t tmp;

	mpz_init(tmp);

	/* measure the contribution of 2 as a factor of sieve
	   values. The multiplier itself must also be taken into
	   account in the score. scores[i] is the correction that
	   is implicitly applied to the size of sieve values for
	   multiplier i; a negative score makes sieve values 
	   smaller, and so is better */

	for (i = 0; i < 20; i++) {
		uint8 curr_mult = smpqs_mult_list[i];
		uint8 knmod8 = (uint8)((curr_mult * mpz_get_ui(n)) % 8);
		double logmult = log((double)curr_mult);

		scores[i] = 0.5 * logmult;
		switch (knmod8) {
		case 1:
			scores[i] -= 2 * LN2;
			break;
		case 5:
			scores[i] -= LN2;
			break;
		case 3:
		case 7:
			scores[i] -= 0.5 * LN2;
			break;
		/* even multipliers start with a handicap */
		}
	}
	num_multipliers = i;

	/* for the rest of the small factor base primes */
	for (i = 1; i < num_primes; i++) {
		uint32 prime = (uint32)spSOEprimes[i];
		double contrib = log((double)prime) / (prime - 1);
		uint32 modp = (uint32)mpz_tdiv_ui(n,prime);

		for (j = 0; j < num_multipliers; j++) {
			uint8 curr_mult = smpqs_mult_list[j];
			uint32 knmodp = (modp * curr_mult) % prime;

			mpz_set_ui(tmp, knmodp);

			/* if prime i is actually in the factor base
			   for k * n ... */

			if (knmodp == 0 || mpz_kronecker_ui(tmp, prime) == 1) { //jacobi_1(knmodp, prime) == 1) {

				/* ...add its contribution. A prime p con-
				   tributes log(p) to 1 in p sieve values, plus
				   log(p) to 1 in p^2 sieve values, etc. The
				   average contribution of all multiples of p 
				   to a random sieve value is thus

				   log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
				   = (log(p) / p) * 1 / (1 - (1/p)) 
				   = log(p) / (p-1)

				   This contribution occurs once for each
				   square root used for sieving. There are two
				   roots for each factor base prime, unless
				   the prime divides k*n. In that case there 
				   is only one root */

				//printf("scores[%d] = %f\n",j,contrib);
				if (knmodp == 0)
					scores[j] -= contrib;
				else
					scores[j] -= 2 * contrib;
			}
		}

	}

	/* use the multiplier that generates the best score */
	best_score = 1000.0;
	best_mult = 1;
	for (i = 0; i < num_multipliers; i++) {
		
		double score = scores[i];
		if (score < best_score) {
			best_score = score;
			best_mult = smpqs_mult_list[i];
		}
	}

	mpz_clear(tmp);
	return best_mult;
}

static int smpqs_check_relations(uint32 sieve_interval, uint32 blocknum, uint8 *sieve, mpz_t n, sm_mpqs_poly *poly, uint8 closnuf,
						smpqs_sieve_fb *fb,fb_list_sm_mpqs *fullfb, sm_mpqs_rlist *full, sm_mpqs_rlist *partial, 
						uint32 cutoff, uint8 small_cutoff, uint32 start_prime, uint32 parity, uint32 *num, int numpoly)
{
	mpz_t Q,t1,t2,t3;
	uint32 offset,i,j;
	uint32 neg;
	uint64 *sieveblock;
	uint32 limit = 32768 >> 3;
		
	sieveblock = (uint64 *)sieve;
	mpz_init(Q);
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(t3);

	//total_locs += limit;

#ifdef SM_SIMD_SIEVE_SCAN_VEC

	for (j=0;j<limit;j+=8)	
	{		
		uint32 result;
		uint8 buffer[64];
		
		SM_SIEVE_SCAN_64_VEC;

		if (result == 0)
			continue;
		
		for (i=0; i<result; i++)
		{
			uint32 thisloc = (j << 3) + (uint32)buffer[i];

			//rarely, we get a result back that doesn't have the high bit set.
			//not sure why yet, but check for it here.
			//if ((sieve[thisloc] & 0x80) == 0)
			//	continue;

			(*num)++;

			offset = (blocknum<<15) + thisloc;
			mpz_set_64(t2, poly->poly_b);
			mpz_mul_2exp(t2, t2, 1); //zShiftLeft_1(&t2,&poly->poly_b);

			mpz_set_64(t1, poly->poly_a);
			mpz_mul_ui(t1, t1, offset); //zShortMul(&poly->poly_a,offset,&t1);
			if (parity)
				mpz_sub(t3, t1, t2); //zSub(&t1,&t2,&t3);
			else
				mpz_add(t3, t1, t2); //zAdd(&t1,&t2,&t3);

			mpz_mul_ui(t1, t3, offset); //zShortMul(&t3,offset,&t1);
			mpz_add(Q, t1, poly->poly_c); //zAdd(&t1,&poly->poly_c,&Q);
			if (mpz_sgn(Q) < 0)
			{
				neg = 1;
				mpz_neg(Q, Q);
			}
			else
				neg = 0;
				
			smpqs_trial_divide_Q(Q,fb,full,partial,sieve,offset,
				thisloc,neg,fullfb,cutoff,small_cutoff,start_prime,
				numpoly,parity,closnuf,sieve[thisloc]);
		}
	}

	SM_SCAN_CLEAN;

#elif defined(SM_SIMD_SIEVE_SCAN)

	for (j=0;j<limit;j+=8)	
	{
		uint32 result, k;

		SM_SIEVE_SCAN_64;

		if (result == 0)
			continue;

		for (i=0; i<8; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SM_SCAN_MASK) == (uint64)(0))
				continue;

			//at least one passed the check, find which one(s) and pass to 
			//trial division stage
			for (k=0;k<8;k++)
			{
				uint32 thisloc = ((j+i)<<3) + k;
				if ((sieve[thisloc] & 0x80) == 0)
					continue;

				(*num)++;

				offset = (blocknum<<15) + thisloc;
				mpz_set_64(t2, poly->poly_b);
				mpz_mul_2exp(t2, t2, 1); //zShiftLeft_1(&t2,&poly->poly_b);

				mpz_set_64(t1, poly->poly_a);
				mpz_mul_ui(t1, t1, offset); //zShortMul(&poly->poly_a,offset,&t1);
				if (parity)
					mpz_sub(t3, t1, t2); //zSub(&t1,&t2,&t3);
				else
					mpz_add(t3, t1, t2); //zAdd(&t1,&t2,&t3);

				mpz_mul_ui(t1, t3, offset); //zShortMul(&t3,offset,&t1);
				mpz_add(Q, t1, poly->poly_c); //zAdd(&t1,&poly->poly_c,&Q);
				if (mpz_sgn(Q) < 0)
				{
					neg = 1;
					mpz_neg(Q, Q);
				}
				else
					neg = 0;
				
				smpqs_trial_divide_Q(Q,fb,full,partial,sieve,offset,
					thisloc,neg,fullfb,cutoff,small_cutoff,start_prime,
					numpoly,parity,closnuf,sieve[thisloc]);
			}
		}
	}

	SM_SCAN_CLEAN;

#else


	for (j=0;j<limit;j+=8)	
	{
		uint32 k;

		if (((sieveblock[j] | sieveblock[j+1] | sieveblock[j+2] | sieveblock[j+3] |
		      sieveblock[j+4] | sieveblock[j+5] | sieveblock[j+6] | sieveblock[j+7]
			) & SM_SCAN_MASK) == (uint64)(0))
			continue;

		for (i=0; i<8; i++)
		{
			//check 8 locations simultaneously
			if ((sieveblock[j + i] & SM_SCAN_MASK) == (uint64)(0))
				continue;

			//at least one passed the check, find which one(s) and pass to 
			//trial division stage
			for (k=0;k<8;k++)
			{
				uint32 thisloc = ((j+i)<<3) + k;
				if ((sieve[thisloc] & 0x80) == 0)
					continue;

				(*num)++;

				offset = (blocknum<<15) + thisloc;
				mpz_set_64(t2, poly->poly_b);
				mpz_mul_2exp(t2, t2, 1); //zShiftLeft_1(&t2,&poly->poly_b);

				mpz_set_64(t1, poly->poly_a);
				mpz_mul_ui(t1, t1, offset); //zShortMul(&poly->poly_a,offset,&t1);
				if (parity)
					mpz_sub(t3, t1, t2); //zSub(&t1,&t2,&t3);
				else
					mpz_add(t3, t1, t2); //zAdd(&t1,&t2,&t3);

				mpz_mul_ui(t1, t3, offset); //zShortMul(&t3,offset,&t1);
				mpz_add(Q, t1, poly->poly_c); //zAdd(&t1,&poly->poly_c,&Q);
				if (mpz_sgn(Q) < 0)
				{
					neg = 1;
					mpz_neg(Q, Q);
				}
				else
					neg = 0;
				
				smpqs_trial_divide_Q(Q,fb,full,partial,sieve,offset,
					thisloc,neg,fullfb,cutoff,small_cutoff,start_prime,
					numpoly,parity,closnuf,sieve[thisloc]);
			}
		}
	}

#endif

	mpz_clear(Q);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(t3);
	return full->num_r;
}

#define DIVIDE_ONE_PRIME \
	while (mpz_tdiv_ui(Q, prime) == 0) \
	{						\
		fboffset[++smooth_num] = i;	\
		mpz_tdiv_q_ui(Q, Q, prime); 	\
	}

static void smpqs_trial_divide_Q(mpz_t Q, smpqs_sieve_fb *fb, sm_mpqs_rlist *full, sm_mpqs_rlist *partial,
						  uint8 *sieve, uint32 offset, uint32 j, uint32 sign, fb_list_sm_mpqs *fullfb, uint32 cutoff,
						  uint8 small_cutoff, uint32 start_prime, int numpoly, uint32 parity,uint8 closnuf,uint8 bits)
{
	smpqs_sieve_fb *fbptr;
	uint32 i,num_f,num_p;
	uint32 root1,root2,prime;
	uint32 r;
	int smooth_num;
	uint16 fboffset[SM_MAX_SMOOTH_PRIMES];
	uint8 logp;

	num_f = full->num_r;
	num_p = partial->num_r;
	
	//we have two signs to worry about.  the sign of the offset tells us how to calculate ax + b, while
	//the sign of Q(x) tells us how to factor Q(x) (with or without a factor of -1)
	//the square root phase will need to know both.  fboffset holds the sign of Q(x).  the sign of the 
	//offset is stored standalone in the relation structure.
	if (sign)
		fboffset[0] = 1;
	else
		fboffset[0] = 0;

	smooth_num=0;
	bits = (255-sieve[j]) + closnuf + 1;

	//take care of powers of two
	while (mpz_even_p(Q))
	{
		//right shift Q
		mpz_tdiv_q_2exp(Q,Q,1);
		fboffset[++smooth_num] = 1;
	}

	//completely unrolled trial division by primes that we have not 
	//sieved via the small prime variation
	{
		uint64 q64;
		uint32 tmp1;

		fbptr = fb + 2;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[2];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[2];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 2;
				mpz_tdiv_q_ui(Q,Q,prime); //zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (mpz_tdiv_ui(Q,prime) == 0); //zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 3;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[3];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[3];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 3;
				mpz_tdiv_q_ui(Q,Q,prime); //zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (mpz_tdiv_ui(Q,prime) == 0); //zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 4;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[4];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[4];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 4;
				mpz_tdiv_q_ui(Q,Q,prime); //zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (mpz_tdiv_ui(Q,prime) == 0); //zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 5;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[5];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[5];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 5;
				mpz_tdiv_q_ui(Q,Q,prime); //zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (mpz_tdiv_ui(Q,prime) == 0); //zShortMod32(&Q32,prime) == 0);
		}

		fbptr = fb + 6;
		root1 = fbptr->roots >> 16;	
		root2 = fbptr->roots & 0xffff;

		prime = fbptr->prime_and_logp >> 16;
		logp = fbptr->prime_and_logp & 0xff;

		tmp1 = offset + fullfb->list->correction[6];
		q64 = (uint64)tmp1 * (uint64)fullfb->list->small_inv[6];
		tmp1 = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp1 = offset - tmp1 * prime;

		if (tmp1 == root1 || tmp1 == root2)
		{
			do
			{
				fboffset[++smooth_num] = 6;
				mpz_tdiv_q_ui(Q,Q,prime); //zShortDiv32(&Q32,prime,&Q32);
				bits += logp;
			} while (mpz_tdiv_ui(Q,prime) == 0); //zShortMod32(&Q32,prime) == 0);
		}
	}

	if (bits < (closnuf + small_cutoff))
		return;

	i=start_prime;
	while (i < fullfb->B)
	{
		uint32 tmp;
		uint64 q64;

		fbptr = fb + i;
		root1 = (fbptr->roots >> 16); // + 32768 - j;	
		root2 = (fbptr->roots & 0xffff); // + 32768 - j;
		prime = fbptr->prime_and_logp >> 16;

		if (prime > 256)
			break;
		
		tmp = offset + fullfb->list->correction[i];
		q64 = (uint64)tmp * (uint64)fullfb->list->small_inv[i];
		tmp = q64 >> 32; 
		//at this point tmp1 is offset / prime
		tmp = offset - tmp * prime;

		if (tmp == root1 || tmp == root2)
			DIVIDE_ONE_PRIME;			

		i++;
	}

	while (i < fullfb->B)
	{
		uint32 tmp;
		uint64 q64;

		fbptr = fb + i;
		root1 = (fbptr->roots >> 16); // + 32768 - j;		
		root2 = (fbptr->roots & 0xffff); // + 32768 - j;
		prime = fbptr->prime_and_logp >> 16;

		tmp = offset + fullfb->list->correction[i];
		q64 = (uint64)tmp * (uint64)fullfb->list->small_inv[i];
		tmp = q64 >> 40; 
		//at this point tmp1 is offset / prime
		tmp = offset - tmp * prime;

		if (tmp == root1 || tmp == root2)
			DIVIDE_ONE_PRIME;		

		i++;
	}

	//check if it completely factored by looking at the unfactored portion in tmp
	if (mpz_cmp_ui(Q,1) == 0)
	{
		if (full->num_r == full->allocated) 
		{
			//printf("\nreallocating fulls\n");
			r = full->allocated;
			full->allocated *= 2;
			full->list = (sm_mpqs_r **)realloc(full->list, 
					full->allocated * sizeof(sm_mpqs_r *));
		}
		smpqs_save_relation(full,offset,1,smooth_num+1,num_f,fboffset,numpoly,parity);
	}
	else if (mpz_cmp_ui(Q,cutoff) < 0)
	{
		smpqs_save_relation(partial,offset,mpz_get_ui(Q),smooth_num+1,num_p,fboffset,numpoly,parity);

		if (partial->num_r == partial->allocated) 
		{
			//printf("\nreallocating partials\n");
			r = partial->allocated;
			partial->allocated *= 2;
			partial->list = (sm_mpqs_r **)realloc(partial->list, 
					partial->allocated * sizeof(sm_mpqs_r *));
		}
	}

	return;
}

static void smpqs_save_relation(sm_mpqs_rlist *list, uint32 offset, uint32 largeprime, uint32 num_factors, 
						  uint32 rnum, uint16 *fboffset, int numpoly, uint32 parity)
{
	uint32 i;
	list->list[rnum] = (sm_mpqs_r *)malloc(sizeof(sm_mpqs_r));
	list->list[rnum]->fboffset = (uint16 *)malloc(num_factors*sizeof(uint16));
	for (i=0;i<num_factors;i++)
		list->list[rnum]->fboffset[i] = fboffset[i];
	
	list->list[rnum]->offset = offset;
	list->list[rnum]->largeprime = largeprime;
	list->list[rnum]->parity = parity;
	list->list[rnum]->num_factors = (uint8)(num_factors);
	list->list[rnum]->polynum = numpoly;
	list->num_r++;
	return;
}

static void smpqs_sieve_block(uint8 *sieve, smpqs_sieve_fb *fb, uint32 start_prime, 
	uint8 s_init, fb_list_sm_mpqs *fullfb)
{
	uint32 prime, root1, root2, stop;
	uint32 B=fullfb->B;
	uint32 i, fourp;
	uint8 logp, *s2, *s3, *s4;
	smpqs_sieve_fb *fbptr;

	//initialize block
	memset(sieve,s_init,32768);
	
	//we've now filled the entire block with the small fb primes, proceed with the rest
	for (i=start_prime;i<B;i++)
	{	
		fbptr = fb + i;
		prime = fbptr->prime_and_logp >> 16;
		root1 = fbptr->roots >> 16;
		root2 = fbptr->roots & 0xffff;
		logp = fbptr->prime_and_logp & 0xff;

		fourp = prime << 2;
		stop = 32768 - fourp + prime; //stop = 32768 - prime;
		s2 = sieve + prime;
		s3 = s2 + prime;
		s4 = s3 + prime;

		while (root2 < stop)
		{
			sieve[root1] -= logp;
			sieve[root2] -= logp;
			s2[root1] -= logp;
			s2[root2] -= logp;
			s3[root1] -= logp;
			s3[root2] -= logp;
			s4[root1] -= logp;
			s4[root2] -= logp;
			root1 += fourp; //(prime << 1);
			root2 += fourp; //(prime << 1);
		}

		while (root2 < 32768)
		{
			sieve[root1] -= logp;
			sieve[root2] -= logp;
			root1 += prime;
			root2 += prime;
		}

	}

	return;
}

void smpqs_nextD(sm_mpqs_poly *poly, mpz_t n)
{
	uint32 r;

	if (poly->side || poly->use_only_p)
	{
		do 
		{
			poly->poly_d_idp++;
			if (poly->poly_d_idp >= szSOEp)
				smpqs_get_more_primes(poly);
			poly->poly_d = spSOEprimes[poly->poly_d_idp];
			r = mpz_tdiv_ui(n,poly->poly_d); //zShortMod(n,*polyd);
		} while ((jacobi_1(r,poly->poly_d) != 1) || ((poly->poly_d & 3) != 3));
	}
	else
	{
		do 
		{
			poly->poly_d_idn--;
			if (poly->poly_d_idn < 5)
			{
				poly->use_only_p = 1;
				smpqs_nextD(poly, n);
				return;
			}
			poly->poly_d = spSOEprimes[poly->poly_d_idn];			
			r = mpz_tdiv_ui(n,poly->poly_d); //zShortMod(n,*polyd);
		} while ((jacobi_1(r,poly->poly_d) != 1) || ((poly->poly_d & 3) != 3));
	}

	poly->side = !poly->side;

	return;
}


void gmpModExp_1(mpz_t a, uint32 b, uint32 m, uint32 *u)
{
	//computes a^b mod m = u using the binary method
	//see, for instance, the handbook of applied cryptography
	mpz_t aa,t;
	uint32 n,bb,ut;

	mpz_init(aa);
	mpz_init(t);

	n=1;
	mpz_set(aa,a);
	bb = b;
	while (bb != 0)
	{
		if (bb & 0x1)
		{
			mpz_mul_ui(t,aa,n); //zShortMul(&aa,n,&t);  //n*a
			n = mpz_tdiv_ui(t,m); //n = zShortMod(&t,m);   //n*a mod m
		}
		bb >>= 1;   
		//compute successive squares of a
		mpz_mul(t,aa,aa); //zSqr(&aa,&t);
		ut = mpz_tdiv_ui(t,m); //n = zShortMod(&t,m);
		mpz_set_ui(aa,ut); //sp2z(ut,&aa);
	}
	*u = n;

	mpz_clear(aa);
	mpz_clear(t);
	return;
}


void smpqs_computeB(sm_mpqs_poly *poly, mpz_t n)
{
	//using poly_d, compute poly_b and poly_a = poly_d^2
	//int i;
	mpz_t t1, t2, t3, t4, t5;
	uint32 ut1;
	uint32 polyd = poly->poly_d;

	mpz_init(t1);
	mpz_init(t2);
	mpz_init(t3);
	mpz_init(t4);
	mpz_init(t5);
	//poly_a = d^2.  we just found d.  also compute b using Hegel theorem and lifting.

	//t0 = n^(d-3)/4 mod d
	gmpModExp_1(n,(polyd-3) >> 2,polyd, &ut1);	

	//t1 = h1 = n*t0 mod d
	mpz_mul_ui(t1, n, ut1); //zShortMul(n,ut1,&t3);		
	mpz_tdiv_r_ui(t1, t1, polyd); //zShortMod(&t3,polyd);

	//t3 = n - h1^2
	mpz_mul(t2, t1, t1);
	mpz_sub(t3, n, t2); //zShortSub(n,ut1 * ut1,&t3);		

	//t4 = (n - h1^2)/d
	mpz_tdiv_q_ui(t4, t3, polyd); //zShortDiv(&t3,polyd,&t4);	

	//t4 = (n - h1^2)/d mod d
	mpz_tdiv_r_ui(t4, t4, polyd); //zShortMod(&t4,polyd);	

	//t5 = 2*h1;
	//zShiftLeft(&t5,&t1,1);		
	mpz_mul_2exp(t5, t1, 1);

	//compute t6 = (2*h1)^-1 mod d = (2*h1)^(d-2) mod d
	gmpModExp_1(t5,polyd-2,polyd,&ut1);

	//compute t3 = h2 = ((2*h1)^-1 * (n - h1^2)/d) mod d
	mpz_mul_ui(t5, t4, ut1); //zMul(&t6,&t4,&t5);		
	mpz_tdiv_r_ui(t3, t5, polyd); //zShortMod(&t5,polyd);

	//compute t5 = h1 + h2*D
	mpz_mul_ui(t4, t3, polyd); //zShortMul(&t3,polyd,&t4);	
	mpz_add(t5, t4, t1); //zAdd(&t4,&t1,&t5);

	//we're now done with d, so compute a = d^2
	poly->poly_a = (uint64)polyd * (uint64)polyd;

	//compute b = h1 + h2*d mod a
	mpz_set_64(t4, poly->poly_a);
	mpz_tdiv_r(t5, t5, t4);
	poly->poly_b = mpz_get_64(t5);

	//make sure b < a/2
	if (poly->poly_b > (poly->poly_a >> 1))
		poly->poly_b = poly->poly_a - poly->poly_b; //sp2z(pa - pb,&poly->poly_b);

	//now that we have b, compute c = (b*b - n)/a
	mpz_set_64(t2, poly->poly_b);
	mpz_mul(t1, t2, t2); //zSqr(&poly->poly_b,&t1);
	mpz_sub(t3, t1, n); //zSub(&t1,n,&t3);
	mpz_tdiv_q(poly->poly_c, t3, t4); //zShortDiv(&t3,pa,&poly->poly_c);

	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(t3);
	mpz_clear(t4);
	mpz_clear(t5);
	return;
}

int sm_check_relation(mpz_t a, mpz_t b, sm_mpqs_r *r, fb_list_sm_mpqs *fb, mpz_t n)
{
	int offset, lp, parity, num_factors;
	int j,retval;
	mpz_t Q, RHS,t1,t2;

	mpz_init(Q);
	mpz_init(RHS);
	mpz_init(t1);
	mpz_init(t2);

	offset = r->offset;
	lp = r->largeprime;
	parity = r->parity;
	num_factors = r->num_factors;

	mpz_set_ui(RHS, lp);
	for (j=1; j<num_factors; j++)
		mpz_mul_ui(RHS, RHS, fb->list->prime[r->fboffset[j]]);; //zShortMul(&RHS,fb->list->prime[r->fb_offsets[j]],&t1);

	//Q(x)/a = (ax + b)^2 - N, where x is the sieve index
	mpz_mul_ui(t1, a, offset); //zShortMul(a,offset,&t1);
	if (parity)
		mpz_sub(t2, t1, b); //zSub(&t1,b,&t2);
	else
		mpz_add(t2, t1, b); //zAdd(&t1,b,&t2);
	mpz_mul(t1, t2, t2); //zMul(&t2,&t2,&t1);
	mpz_sub(Q, t1, n); //zSub(&t1,n,&Q);
	mpz_tdiv_q(Q, Q, a);

	retval = 0;
	if (mpz_sgn(Q) < 0)
	{
		mpz_neg(Q,Q);
	}

	if (mpz_cmp(Q,RHS) != 0)
	{
		printf("failure to equate relation\n");
		gmp_printf("%Zd %Zd\n",Q,RHS);
		retval = 1;
	}

	mpz_clear(Q);
	mpz_clear(RHS);
	mpz_clear(t1);
	mpz_clear(t2);
	return retval;
}

void smpqs_computeRoots(sm_mpqs_poly *poly, fb_list_sm_mpqs *fb, uint32 *modsqrt, 
	smpqs_sieve_fb *fbp, smpqs_sieve_fb *fbn, uint32 start_prime)
{
	//the roots are computed using a and b as follows:
	//(+/-t - b)(a)^-1 mod p
	//assume b > t
	//uint32 root1, root2, prime, amodp;
	int root1, root2, prime, amodp, bmodp, x;
	uint32 i;
	
	for (i=2;i<7;i++)
	{
		uint64 q64, tmp, t2;

		prime = fb->list->prime[i]; 
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		bmodp = poly->poly_b % (uint64)prime;
		x = (int)root1 - bmodp;
		if (x < 0) x += prime; root1 = x;
		x = (int)root2 - bmodp;
		if (x < 0) x += prime; root2 = x;

		amodp = poly->poly_a % (uint64)prime;
		amodp = modinv_1(amodp,prime);

		t2 = (uint64)amodp * (uint64)root1;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root1 = t2 - tmp * prime;

		t2 = (uint64)amodp * (uint64)root2;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root2 = t2 - tmp * prime;

		// we don't sieve these primes, so ordering doesn't matter
		fbp[i].roots = (uint32)((root1 << 16) | root2);
		if (root2 == 0) fbn[i].roots = 0;
		else fbn[i].roots = (uint32)((prime - root2) << 16);
		if (root1 == 0) fbn[i].roots |= 0;
		else fbn[i].roots |= (uint32)(prime - root1);
	}

	for (i=7;i<fb->B;i++)
	{
		uint64 q64, tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		if (prime > 256)
			break;

		bmodp = poly->poly_b % (uint64)prime;
		x = (int)root1 - bmodp;
		if (x < 0) x += prime;
		root1 = x;
		x = (int)root2 - bmodp;
		if (x < 0) x += prime;
		root2 = x;
	
		//now (t - b) mod p is in root1 and (-t - b) mod p is in root2
		//find a^-1 mod p = inv(a mod p) mod p
		amodp = poly->poly_a % (uint64)prime;
		amodp = modinv_1(amodp,prime);

		t2 = (uint64)amodp * (uint64)root1;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root1 = t2 - tmp * prime;

		t2 = (uint64)amodp * (uint64)root2;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 32; 
		root2 = t2 - tmp * prime;

		if (root2 < root1)
		{
			fbp[i].roots = (uint32)((root2 << 16) | root1);
			fbn[i].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
		}
		else
		{
			fbp[i].roots = (uint32)((root1 << 16) | root2);
			fbn[i].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
		}
	}

	for ( ;i<fb->B;i++)
	{
		uint64 q64, tmp, t2;

		prime = fb->list->prime[i];
		root1 = modsqrt[i]; 
		root2 = prime - root1; 

		bmodp = poly->poly_b % (uint64)prime;
		x = (int)root1 - bmodp;
		if (x < 0) x += prime;
		root1 = x;
		x = (int)root2 - bmodp;
		if (x < 0) x += prime;
		root2 = x;
	
		//now (t - b) mod p is in root1 and (-t - b) mod p is in root2
		//find a^-1 mod p = inv(a mod p) mod p
		amodp = poly->poly_a % (uint64)prime;
		amodp = modinv_1(amodp,prime);

		t2 = (uint64)amodp * (uint64)root1;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 40; 
		root1 = t2 - tmp * prime;

		t2 = (uint64)amodp * (uint64)root2;
		tmp = t2 + (uint64)fb->list->correction[i];
		q64 = tmp * (uint64)fb->list->small_inv[i];
		tmp = q64 >> 40; 
		root2 = t2 - tmp * prime;

		if (root2 < root1)
		{
			fbp[i].roots = (uint32)((root2 << 16) | root1);
			fbn[i].roots = (uint32)(((prime - root1) << 16) | (prime - root2));
		}
		else
		{
			fbp[i].roots = (uint32)((root1 << 16) | root2);
			fbn[i].roots = (uint32)(((prime - root2) << 16) | (prime - root1));
		}
	}
	return;
}

int qcomp_smpqs(const void *x, const void *y)
{
	sm_mpqs_r **xx = (sm_mpqs_r **)x;
	sm_mpqs_r **yy = (sm_mpqs_r **)y;
	
	if (xx[0]->largeprime > yy[0]->largeprime)
		return 1;
	else if (xx[0]->largeprime == yy[0]->largeprime)
		return 0;
	else
		return -1;
}

static uint64 smpqs_bitValRead64(uint64 **m, int row, int col);

int smpqs_BlockGauss(sm_mpqs_rlist *full, sm_mpqs_rlist *partial, uint64 *apoly, uint64 *bpoly,
			fb_list_sm_mpqs *fb, mpz_t n, int mul, 
			mpz_t *factors,uint32 *num_factor)
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
	mpz_t zx, zy, tmp, tmp2, tmp3, tmp4, nn,tmp_a,input,zmul;
	
	mpz_init(zx);
	mpz_init(zy);
	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init(tmp3);
	mpz_init(tmp4);
	mpz_init(nn);
	mpz_init(input);
	mpz_init(tmp_a);
	mpz_init(zmul);

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

	// remove the multiplier from the input
	mpz_tdiv_q_ui(input, n, mul); //zShortDiv(n,mul,&input);
	mpz_set_ui(zmul, mul); //sp2z(mul,&zmul);

	//initialize blacklist
	for (i=0;i<num_r;i++) bl[i] = 0;
	//search over all columns, right to left (more sparse on the right side)
	for (i=B-1;i>=0;i--)
	{
		//and all rows
		for (j=0;j<num_r;j++)
		{
			//if the j'th row, i'th bit is 1 and not blacklisted, continue
			bool_val = (smpqs_bitValRead64(m2_64,j,i) != 0) && (bl[j] == 0);
			//bool_val = (((m2_64[(j)][(i >> 6)]) & (1ULL << ((uint64)i & 63ULL))) && (bl[j] == 0));
			if (bool_val)
			{
				//add the j'th row mod 2 to all rows after it with a 1 in the ith column
				for (k=j+1;k<num_r;k++)
				{
					bool_val = (smpqs_bitValRead64(m2_64,k,i) != 0) && (bl[k] == 0);
					//bool_val = (((m2_64[(k)][(i >> 6)]) & (1ULL << ((uint64)i & 63ULL))) && (bl[k] == 0));
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
								bool_val = smpqs_bitValRead64(aug_64,k,l) != 0;
								//bool_val = (((aug_64[(k)][(l >> 6)]) & (1ULL << ((uint64)l & 63ULL)))) != 0;
								if (bool_val)
								{
									//then the l'th row of m was involved
									for (q=0;q<(int)B;q++)
										pd[q] += m[l][q];
								}
							}

							//compute x mod n
							mpz_set_ui(zy, 1); //sm_zcopy(&zOne,&zy);
							mpz_set_ui(zx, 1); //sm_zcopy(&zOne,&zx);
							for (l=0;l<num_r;l++)
							{
								bool_val = smpqs_bitValRead64(aug_64,k,l) != 0;
								//bool_val = (((aug_64[(k)][(l >> 6)]) & (1ULL << ((uint64)l & 63ULL)))) != 0;
								if (bool_val)
								{
									//printf("accumulating relation %d\n",l);
									//then the l'th relation is involved in the product of relations
									if (l >= num_f)
									{
										uint64 pa;
										uint32 d1,d2;
										//if l >= num_f, then this row refers to a relation generated from two partials.
										//we'll need to go back to the two partial locations to find the two offsets to
										//multiply together
										//luckily, we've remembered the index in the complete list of partials that 
										//created this full relation
											
										//our relation is of the form (ax + b)^2 == a(ax^2 + 2bx + c) mod n
										//(ax^2 + 2bx + c) is what we trial divided, and we remembered
										//a and b, so we can form the left hand side easily
											
										// poly_a and b will fit in a uint64 in smallmpqs

										//recreate poly_b from poly_a (store instead??)
										polynum = partial->list[partial_index[l-num_f]]->polynum;
											
										//compute Q1(x)
										pa = apoly[polynum];
										d1 = (uint32)sqrt((int64)pa);
										
										mpz_set_64(tmp, apoly[polynum]);
										mpz_set_64(tmp3, bpoly[polynum]);
										mpz_mul_ui(tmp, tmp, partial->list[partial_index[l-num_f]]->offset); //zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]]->offset,&tmp);
										if (partial->list[partial_index[l-num_f]]->parity)
											mpz_sub(tmp, tmp, tmp3); //zShortSub(&tmp,pb,&tmp);
										else
											mpz_add(tmp, tmp, tmp3); //zShortAdd(&tmp,pb,&tmp);

										//include 'a'
										mpz_mul_ui(tmp2, zy, d1); //zShortMul(&zy,d1,&tmp2);
										mpz_tdiv_r(zy, tmp2, n); //zDiv(&tmp2,n,&tmp3,&zy);

										//compute Q2(x)
										polynum = partial->list[partial_index[l-num_f]-1]->polynum;
											
										//compute Q1(x)
										pa = apoly[polynum]; 
										d2 = (uint32)sqrt((int64)pa);
										
										mpz_set_64(tmp3, apoly[polynum]);
										mpz_set_64(tmp4, bpoly[polynum]);
										mpz_mul_ui(tmp3, tmp3, partial->list[partial_index[l-num_f]-1]->offset); //zShortMul(&apoly[polynum],partial->list[partial_index[l-num_f]-1]->offset,&tmp3);
										if (partial->list[partial_index[l-num_f]-1]->parity)
											mpz_sub(tmp2, tmp3, tmp4); //zShortSub(&tmp3,pb,&tmp2);
										else
											mpz_add(tmp2, tmp3, tmp4); //zShortAdd(&tmp3,pb,&tmp2);

										//compute Q(x1)*Q(x2)
										mpz_mul(tmp4, tmp, tmp2); //zMul(&tmp,&tmp2,&tmp4);	
										mpz_mul(tmp, zx, tmp4); //zMul(&zx,&tmp4,&tmp);		//accumulate with previous terms
										mpz_tdiv_r(zx, tmp, n); //zDiv(&tmp,n,&tmp2,&zx);		//mod n

										//include the large prime in mp_y
										mpz_mul_ui(tmp2, zy, partial->list[partial_index[l-num_f]]->largeprime); //zShortMul(&zy,partial->list[partial_index[l-num_f]]->largeprime,&tmp2);
										mpz_tdiv_r(zy, tmp2, n); //zDiv(&tmp2,n,&tmp3,&zy);

										//include 'a'
										mpz_mul_ui(tmp2, zy, d2); //zShortMul(&zy,d2,&tmp2);	
										mpz_tdiv_r(zy, tmp2, n); //zDiv(&tmp2,n,&tmp3,&zy);
									}
									else
									{
										uint64 pa,d1;

										//recreate poly_b from poly_a (store instead??)
										polynum = full->list[l]->polynum;
											
										//compute Q1(x)
										mpz_set_64(tmp, apoly[polynum]);
										mpz_set_64(tmp4, bpoly[polynum]);
										mpz_mul_ui(tmp, tmp, full->list[l]->offset); //zShortMul(&apoly[polynum],full->list[l]->offset,&tmp);
										pa = apoly[polynum]; //apoly[polynum].val[0];
										d1 = (uint32)sqrt((int64)pa);

										if (full->list[l]->parity)
											mpz_sub(nn, tmp, tmp4); //zShortSub(&tmp,pb,&nn);
										else
											mpz_add(nn, tmp, tmp4); //zShortAdd(&tmp,pb,&nn);

										mpz_mul(tmp, zx, nn); //zMul(&zx,&nn,&tmp);			//accumulate with previous terms
										mpz_tdiv_r(zx, tmp, n); //zDiv(&tmp,n,&tmp2,&zx);		//mod n

										mpz_mul_ui(tmp2, zy, d1); //zShortMul(&zy,d1,&tmp2);	//sqrt(a) = d is part of mp_y
										mpz_tdiv_r(zy, tmp2, n); //zDiv(&tmp2,n,&tmp3,&zy);
									}
								}
							}

							//compute y mod n
							//ignore the factor of -1 in this operation
							for (l=1;l<(int)B;l++)
							{
								if (pd[l] > 0)
								{
									mpz_set_ui(tmp, fb->list->prime[l]); //sp2z(fb->list->prime[l],&tmp);
									//pd tracks the exponents of the smooth factors.  we know they are all even
									//at this point.  we don't want to compute pd^2, so divide by 2.
									//computing the explicit exponentiation and then reducing is
									//slightly faster than doing modexp in smallmpqs.
									//zExp(pd[l]/2,&tmp,&tmp2);
									//zDiv(&tmp2,n,&tmp4,&tmp3);
									mpz_powm_ui(tmp3, tmp, pd[l] / 2, n);
									mpz_mul(tmp4, tmp3, zy); //zMul(&tmp3,&zy,&tmp4);
									mpz_tdiv_r(zy, tmp4, n); //zDiv(&tmp4,n,&tmp2,&zy);
								}
							}

							//split this off into a subroutine... also look for all non-trivial factors if one is composite
							//compute gcd(x-y,n)
							mpz_sub(tmp, zx, zy); //zSub(&zx,&zy,&tmp);
							mpz_gcd(nn, tmp, n); //zLEGCD(&tmp,n,&nn);

							//gmp_printf("gcd is %Zd\n",nn);
							
							/* remove any factors of the multiplier 
							   before saving tmp, and don't save at all
							   if tmp contains *only* multiplier factors */
							if (mul > 1) {
								uint32 ignore_me = spGCD(mul,
										mpz_tdiv_ui(nn, mul)); //zShortMod(&nn, mul));
								if (ignore_me > 1) {
									mpz_tdiv_q_ui(nn, nn, ignore_me); //zShortDiv(&nn, ignore_me, &tmp2);
									if (mpz_cmp_ui(nn, 1) == 0)
										continue;
								}								
							}								

							
							if ((mpz_cmp_ui(nn, 1) > 0) && (mpz_cmp(nn,input) < 0))
							{

								//gmp_printf("checking factor %Zd\n",nn);
								if (mpz_probab_prime_p(nn,5))
								{

									// sometime we find small primes that don't divide the input.
									// ignore these
									if (mpz_sizeinbase(nn,2) < 32)
									{
										if (mpz_tdiv_ui(input, mpz_get_ui(nn)) != 0)
											continue;

										//if (zShortMod(&input,nn.val[0]) != 0)
											//continue;
									}

									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (mpz_cmp(nn,factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									mpz_set(factors[*num_factor],nn);

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}

									//check if we're done by accumulating all factors and comparing to n
									mpz_set(nn,factors[0]);
									for (l=1;l<(int)*num_factor;l++)
										mpz_mul(nn,factors[l],nn); //,&nn);
									if (mpz_cmp(nn,input) == 0)
									{
										//found all factors, done
										goto free;
									}
								}

								//check the other factor
								//sm_zcopy(&input,&tmp);
								//zDiv(&tmp,&nn,&tmp2,&tmp3);
								//sm_zcopy(&tmp2,&tmp);
								
								mpz_tdiv_q(tmp, input, nn);								
	
								//gmp_printf("checking factor %Zd\n",tmp);
								if (mpz_probab_prime_p(tmp,5))
								{
									// sometime we find small primes that don't divide the input.
									// ignore these
									if (mpz_sizeinbase(nn,2) < 32)
									{
										if (mpz_tdiv_ui(input, mpz_get_ui(tmp)) != 0)
											continue;

										//if (zShortMod(&input,nn.val[0]) != 0)
											//continue;
									}

									//check that we havent' already found this one
									set_continue = 0;
									for (l=0;l<(int)*num_factor;l++)
									{
										if (mpz_cmp(tmp,factors[l]) == 0)
											set_continue = 1;
									}
									if (set_continue)
										continue;

									mpz_set(factors[*num_factor],tmp);

									(*num_factor)++;
									if (*num_factor > MAX_FACTORS)
									{
										printf("max number of factors found in block gauss\n");
										goto free;
									}

									//check if we're done by accumulating all factors and comparing to n
									mpz_set(tmp,factors[0]);
									for (l=1;l<(int)*num_factor;l++)
										mpz_mul(tmp,factors[l],tmp); //,&nn);
									if (mpz_cmp(tmp,input) == 0)
									{
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
	r = mpz_tdiv_q_ui(tmp, n, mul); //r = (uint32)zShortDiv(n,mul,&tmp);
	for (i=0;(uint32)i<*num_factor;i++)
	{
		//sm_zcopy(&tmp,&nn);
		//zDiv(&nn,&factors[i],&tmp,&tmp2);
		mpz_tdiv_q(tmp, tmp, factors[i]);
	}

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

	free(bl);
	mpz_clear(zx);
	mpz_clear(zy);
	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
	mpz_clear(tmp4);
	mpz_clear(nn);
	mpz_clear(tmp_a);
	mpz_clear(input);
	mpz_clear(zmul);
	return 0;
}

static uint64 smpqs_masks64[64] = {0x1,0x2,0x4,0x8,
							0x10,0x20,0x40,0x80,
							0x100,0x200,0x400,0x800,
							0x1000,0x2000,0x4000,0x8000,
							0x10000,0x20000,0x40000,0x80000,
							0x100000,0x200000,0x400000,0x800000,
							0x1000000,0x2000000,0x4000000,0x8000000,
							0x10000000,0x20000000,0x40000000,0x80000000ULL,
							0x100000000ULL,0x200000000ULL,0x400000000ULL,0x800000000ULL,
							0x1000000000ULL,0x2000000000ULL,0x4000000000ULL,0x8000000000ULL,
							0x10000000000ULL,0x20000000000ULL,0x40000000000ULL,0x80000000000ULL,
							0x100000000000ULL,0x200000000000ULL,0x400000000000ULL,0x800000000000ULL,
							0x1000000000000ULL,0x2000000000000ULL,0x4000000000000ULL,0x8000000000000ULL,
							0x10000000000000ULL,0x20000000000000ULL,0x40000000000000ULL,0x80000000000000ULL,
							0x100000000000000ULL,0x200000000000000ULL,0x400000000000000ULL,0x800000000000000ULL,
							0x1000000000000000ULL,0x2000000000000000ULL,0x4000000000000000ULL,0x8000000000000000ULL};



static uint64 smpqs_bitValRead64(uint64 **m, int row, int col)
{
	//col is the column in 0 to B-1 representation
	//read the bit in the packed 64 bit representation of the appropriate row
	//don't bother to check the bounds of m w.r.t row and col, assume caller knows what it's doing
	//return 0 if bit not set, 1 << bit offset otherwize
	int offset, mcol;
	mcol = col >> 6;
	offset = col & 63;
	return (m[row][mcol] & smpqs_masks64[offset]);
}





