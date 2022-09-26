/*
Copyright (c) 2014, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/

#include "gmp.h"
#include "soe.h"
#include "monty.h"
#include <stdint.h>

#define D 120

//#define DEBUG 1
#define base_t uint64_t
#define base_digits 2

typedef struct
{
	uint64_t base[2];
} u128_t;

typedef struct
{
	base_t X[base_digits];
	base_t Z[base_digits];
} tinyecm_pt;

typedef struct
{
	base_t sum1[base_digits];
	base_t diff1[base_digits];
	base_t sum2[base_digits];
	base_t diff2[base_digits];
	base_t tt1[base_digits];
	base_t tt2[base_digits];
	base_t tt3[base_digits];
	base_t tt4[base_digits];
	base_t tt5[base_digits];
	base_t s[base_digits];
	base_t n[base_digits];
	tinyecm_pt pt1;
	tinyecm_pt pt2;
	tinyecm_pt pt3;
	tinyecm_pt pt4;
	tinyecm_pt pt5;
	uint32_t sigma;

	tinyecm_pt Pa;
	tinyecm_pt Pd;
	tinyecm_pt Pad;
	tinyecm_pt Pb[20];
	base_t Paprod[base_digits];
	base_t Pbprod[20][base_digits];
	
	base_t stg2acc[base_digits];
	uint32_t stg1Add;
	uint32_t stg1Doub;
	uint32_t paired;
	uint32_t ptadds;
	uint64_t numprimes;
	uint64_t A;
	uint32_t last_pid;
	uint32_t amin;

	uint32_t stg1_max;
	uint32_t stg2_max;

} tinyecm_work;

static const uint32_t map[60] = {
	0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
	0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
	0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
	0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
	0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
	0, 0, 0, 16, 0, 0, 0, 0, 0, 17 };

static uint64_t* tecm_primes;
static uint64_t tecm_nump;
static uint64_t tecm_minp;
static uint64_t tecm_maxp;
static int tecm_primes_initialized = 0;

// local functions
void add(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P1, tinyecm_pt *P2, 
	tinyecm_pt *Pin, tinyecm_pt *Pout);
void duplicate(monty128_t *mdata, tinyecm_work *work, uint64_t * insum, uint64_t * indiff, tinyecm_pt *P);
void prac(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P, uint64_t c, double v);
int check_factor(uint64_t * Z, uint64_t * n, uint64_t * f);
void build_one_curve(tinyecm_pt *P, monty128_t *mdata, 
	tinyecm_work *work, uint32_t sigma, uint64_t * lcg_state, int verbose);

void ecm_stage1(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P);
void ecm_stage2(tinyecm_pt *P, monty128_t *mdata, tinyecm_work *work);

double bench_curves;
double bench_stg1;
double bench_stg2;

void u128_to_mpz(uint64_t *in, mpz_t out)
{
	mpz_set_ui(out, in[1]);
	mpz_mul_2exp(out, out, 64);
	mpz_add_ui(out, out, in[0]);
	return;
}

void mpz_to_u128(mpz_t in, uint64_t *out)
{
    out[0] = mpz_get_ui(in);
    mpz_tdiv_q_2exp(in, in, 64);
    out[1] = mpz_get_ui(in);

    // restore input
    mpz_mul_2exp(in, in, 64);
    mpz_add_ui(in, in, out[0]);

	return;
}

__inline void copy128(uint64_t *src, uint64_t *dest)
{
	dest[0] = src[0];
	dest[1] = src[1];
	return;
}

__inline void swap128(uint64_t *a, uint64_t *b)
{
	uint64_t tmp[2];
	tmp[0] = a[0];
	tmp[1] = a[1];
	a[0] = b[0];
	a[1] = b[1];
	b[0] = tmp[0];
	b[1] = tmp[1];
	return;
}

__inline void rot128(uint64_t *a, uint64_t *b, uint64_t *c)
{
	uint64_t tmp[2];
	tmp[0] = a[0];
	tmp[1] = a[1];
	a[0] = b[0];
	a[1] = b[1];
	b[0] = c[0];
	b[1] = c[1];
	c[0] = tmp[0];
	c[1] = tmp[1];
	return;
}

void tinyecm_work_init(tinyecm_work *work)
{
	work->stg1Add = 0;
	work->stg1Doub = 0;

	return;
}

void tinyecm_work_free(tinyecm_work *work)
{
	return;
}

void add(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P1, tinyecm_pt *P2, 
	tinyecm_pt *Pin, tinyecm_pt *Pout)
{
	// compute:
	//x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	//z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	// where:
	//x- = original x
	//z- = original z
	addmod128(P1->X, P1->Z, work->sum1, work->n);
	submod128(P1->X, P1->Z, work->diff1, work->n);
	addmod128(P2->X, P2->Z, work->sum2, work->n);
	submod128(P2->X, P2->Z, work->diff2, work->n);

	mulmod128(work->diff1, work->sum2, work->tt1, mdata);	//U
	mulmod128(work->sum1, work->diff2, work->tt2, mdata);	//V

	addmod128(work->tt1, work->tt2, work->tt3, work->n);
	submod128(work->tt1, work->tt2, work->tt4, work->n);
	sqrmod128(work->tt3, work->tt1, mdata);					//(U + V)^2
	sqrmod128(work->tt4, work->tt2, mdata);					//(U - V)^2

	// choosing the initial point Pz0 = 1 means that z_p-q = 1 and this mul isn't necessary...
	// but that involves a different way to initialize curves, so for now
	// we can't assume Z=1
	if (Pin->X == Pout->X)
	{
		mulmod128(work->tt1, Pin->Z, Pout->Z, mdata);		//Z * (U + V)^2
		mulmod128(work->tt2, Pin->X, Pout->X, mdata);		//x * (U - V)^2
		swap128(Pout->Z, Pout->X);
	}
	else
	{
		mulmod128(work->tt1, Pin->Z, Pout->X, mdata);		//Z * (U + V)^2
		mulmod128(work->tt2, Pin->X, Pout->Z, mdata);		//x * (U - V)^2
	}
	work->stg1Add++;
	return;
}

void duplicate(monty128_t *mdata, tinyecm_work *work, 
	uint64_t * insum, uint64_t * indiff, tinyecm_pt *P)
{
	sqrmod128(indiff, work->tt1, mdata);						// U=(x1 - z1)^2
	sqrmod128(insum,  work->tt2, mdata);						// V=(x1 + z1)^2
	mulmod128(work->tt1, work->tt2, P->X, mdata);				// x=U*V

	submod128(work->tt2, work->tt1, work->tt3, mdata->n);	    // w = V-U
	mulmod128(work->tt3, work->s, work->tt2, mdata);			// w = (A+2)/4 * w
	addmod128(work->tt2, work->tt1, work->tt2, mdata->n);       // w = w + U
	mulmod128(work->tt2, work->tt3, P->Z, mdata);				// Z = w*(V-U)
	work->stg1Doub++;
	return;
}

#define ADD 6.0
#define DUP 5.0
#define NV 10  

double getEcost(uint64_t d, uint64_t e)
{
	int doub = 0, add = 0;

	while (d > 0)
	{
		if ((e / 2) < d)
		{
			d = e - d;
		}
		else if ((d < (e / 4)) && ((e & 1) == 0))
		{
			e = e / 2;
			doub++;
			add++;
		}
		else
		{
			e = e - d;
			add++;
		}

	}
	return (doub + add) * 2 + add * 4 + doub * 3;
}

static double lucas_cost(uint64_t n, double v)
{
	uint64_t d, e, r;
	double c; /* cost */

	d = n;
	r = (uint64_t)((double)d * v + 0.5);
	if (r >= n)
		return (ADD * (double)n);
	d = n - r;
	e = 2 * r - n;
	c = DUP + ADD; /* initial duplicate and final addition */
	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
		}
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{ /* condition 1 */
			d = (2 * d - e) / 3;
			e = (e - d) / 2;
			c += 3.0 * ADD; /* 3 additions */
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{ /* condition 2 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		else if ((d + 3) / 4 <= e)
		{ /* condition 3 */
			d -= e;
			c += ADD; /* one addition */
		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		/* now d is odd and e is even */
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;
			c += 3.0 * ADD + DUP; /* three additions, one duplicate */
		}
		else /* necessarily e is even: catches all cases */
		{ /* condition 9 */
			e /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
	}

	if (d != 1)
	{
		c = 9999999.;
	}

	return c;
}

void lucas_opt(uint32_t B1, int num)
{
	uint64_t d, e, r, q, c;
	double cmin = 99999999., cost;
	int best_index = 0;
	int i, j;
	gmp_randstate_t gmprand;
	double *val;
	double trad_cost;
	uint64_t word;
	double *v;
	double *pcost;
	uint64_t f[128];
	int nump;

	gmp_randinit_default(gmprand);

	val = (double *)malloc((num + NV) * sizeof(double));
	val[0] = 0.61803398874989485;
	val[1] = 0.72360679774997897;
	val[2] = 0.58017872829546410;
	val[3] = 0.63283980608870629;
	val[4] = 0.61242994950949500;
	val[5] = 0.62018198080741576;
	val[6] = 0.61721461653440386;
	val[7] = 0.61834711965622806;
	val[8] = 0.61791440652881789;
	val[9] = 0.61807966846989581;

	for (i = NV; i < (num + NV); i++)
	{
		uint64_t ur = gmp_urandomb_ui(gmprand, 52);
		double r = 0.51 + ((double)ur * 2.2204460492503130808472633361816e-16) * 0.25;
		val[i] = r;
	}

	c = 1;
	e = 1;
	j = 0;
	i = 2;
	q = tecm_primes[i];
	trad_cost = 0.;

	while (q < 1000)
	{
		i++;
		q = tecm_primes[i];
	}
	nump = i;

	v = (double *)malloc(i * sizeof(double));
	pcost = (double *)malloc(i * sizeof(double));

	i = 2;
	q = tecm_primes[i];
	pcost[0] = 5.0;
	pcost[1] = 11.0;

	while (q < 1000)
	{
		printf("now optimizing prime %lu\n", q);
		for (d = 0, cmin = ADD * (double)q; d < (num + NV); d++)
		{
			cost = lucas_cost(q, val[d]);
			if (cost < cmin)
			{
				cmin = cost;
				best_index = d;
				printf("best cost is now %1.6f at index %d\n", cmin, d);
			}
		}

		printf("minimum prac cost of prime %lu is %1.6f at index %d, val %1.18f\n",
			q, cmin, best_index, val[best_index]);

		pcost[i] = cmin;
		v[i] = val[best_index];

		i++;
		q = tecm_primes[i];

	}

	printf("now optimizing composites\n");
	c = 6;
	while (c < 1000000)
	{
		int skip = 0;
		double fcost;

		j = 2;
		q = tecm_primes[j];
		while (q < c)
		{
			q = tecm_primes[j++];
			if (c == q)
			{
				skip = 1;
				break;
			}
		}

		if (skip)
		{
			c++;
			continue;
		}

		e = c;
		while ((e & 1) == 0)
			e /= 2;

		//printf("now optimizing composite %lu\n", e);
		for (d = 0, cmin = ADD * (double)e; d < (num + NV); d++)
		{
			cost = lucas_cost(e, val[d]);
			if (cost < cmin)
			{
				cmin = cost;
				best_index = d;
				//printf("best cost is now %1.6f at index %d\n", cmin, d);
			}
		}

		fcost = 0.0;
		d = e;
		j = 0;
		for (i = 0; i < nump; i++)
		{
			if (d % tecm_primes[i] == 0)
			{
				fcost += pcost[i];
				d /= tecm_primes[i];
				f[j++] = tecm_primes[i];
			}

			if (d == 1)
				break;
		}

		if (fcost > cmin)
		{
			printf("minimum prac cost of %lu is %1.1f at index %d, "
				"val %1.18f, factor cost of [", e, cmin, best_index, val[best_index]);
			for (i = 0; i < j; i++)
				printf("%lu,", f[i]);
			printf("\b] is %1.1f\n", fcost);
			fflush(stdout);
		}
		c++;

	}

	exit(0);

	while (q < B1)
	{
		//printf("now accumulating prime %lu, word is %lu\n", q, c);
		if ((e * q) < B1)
		{
			// add another of this prime to the word
			// as long as it still fits in the word.
			if ((0xFFFFFFFFFFFFFFFFULL / c) > q)
			{
				c *= q;
				e *= q;

				//printf("now optimizing traditional word %lu\n", q);
				for (d = 0, cmin = ADD * (double)q; d < (num + NV); d++)
				{
					cost = lucas_cost(q, val[d]);
					if (cost < cmin)
					{
						cmin = cost;
						best_index = d;
					}
				}

				trad_cost += cmin;
			}
			else
			{
				
				// word is full, proceed to optimize
				printf("now optimizing word %d = %lu\n", j, c);

				for (d = 0, cmin = ADD * (double)c; d < (num + NV); d++)
				{
					cost = lucas_cost(c, val[d]);
					if (cost < cmin)
					{
						cmin = cost;
						best_index = d;
						//printf("best cost is now %1.6f at index %d\n", cmin, d);
					}
				}

				printf("minimum prac cost is %1.6f at index %d, val %1.18f\n",
					cmin, best_index, val[best_index]);

				for (d = 0, cmin = ADD * (double)c; d < (num + NV); d++)
				{
					uint64_t dd = (uint64_t)((double)c * val[d]);
					uint64_t ee = c;

					if (spGCD(dd, ee) != 1)
					{
						continue;
					}

					cost = getEcost(dd, ee);
					if (cost < cmin)
					{
						cmin = cost;
						best_index = d;
					}

				}

				printf("minimum euclid cost is %1.6f at index %d, val %1.18f\n",
					cmin, best_index, val[best_index]);
				printf("traditional cost would have been: %1.6f\n", trad_cost);

				// next word
				trad_cost = 0.;
				c = 1;
				e = 1;
				j++;
			}
		}
		else
		{
			// next prime
			i++;
			e = 1;
			q = tecm_primes[i];
		}

	}

	printf("final word %d is %lu\n", j, c);
	for (d = 0, cmin = ADD * (double)c; d < (num + NV); d++)
	{
		cost = lucas_cost(c, val[d]);
		if (cost < cmin)
		{
			cmin = cost;
			best_index = d;
		}
	}

	printf("minimum prac cost is %1.6f at index %d, val %1.18f\n",
		cmin, best_index, val[best_index]);

	for (d = 0, cmin = ADD * (double)c; d < (num + NV); d++)
	{
		uint64_t dd = (uint64_t)((double)c * val[d]);
		uint64_t ee = c;

		if (spGCD(dd, ee) != 1)
		{
			continue;
		}

		cost = getEcost(dd, ee);
		if (cost < cmin)
		{
			cmin = cost;
			best_index = d;
		}

	}

	printf("minimum euclid cost is %1.6f at index %d, val %1.18f\n",
		cmin, best_index, val[best_index]);
	printf("traditional cost would have been: %1.6f\n", trad_cost);

	free(val);
	return;
}

void lucas_opt2(uint32_t B1, int num)
{
	uint64_t d, e, r, q, c;
	double cmin = 99999999., cost;
	int best_index = 0;
	int i, j, k, l;
	gmp_randstate_t gmprand;
	double *val;
	double trad_cost;
	uint64_t word;
	double *v;
	double *pcost;
	uint64_t f[128];
	int nump;

	gmp_randinit_default(gmprand);

	val = (double *)malloc((num + NV) * sizeof(double));
	val[0] = 0.61803398874989485;
	val[1] = 0.72360679774997897;
	val[2] = 0.58017872829546410;
	val[3] = 0.63283980608870629;
	val[4] = 0.61242994950949500;
	val[5] = 0.62018198080741576;
	val[6] = 0.61721461653440386;
	val[7] = 0.61834711965622806;
	val[8] = 0.61791440652881789;
	val[9] = 0.61807966846989581;

	for (i = NV; i < (num + NV); i++)
	{
		uint64_t ur = gmp_urandomb_ui(gmprand, 52);
		double r = 0.51 + ((double)ur * 2.2204460492503130808472633361816e-16) * 0.25;
		val[i] = r;
	}

	c = 1;
	e = 1;
	j = 0;
	i = 2;
	q = tecm_primes[i];
	trad_cost = 0.;

	while (q < 1000)
	{
		i++;
		q = tecm_primes[i];
	}
	nump = i;

	v = (double *)malloc(i * sizeof(double));
	pcost = (double *)malloc(i * sizeof(double));

	i = 2;
	q = tecm_primes[i];
	pcost[0] = 5.0;
	pcost[1] = 11.0;

	while (q < 1000)
	{
		printf("now optimizing prime %lu\n", q);
		for (d = 0, cmin = ADD * (double)q; d < (num + NV); d++)
		{
			cost = lucas_cost(q, val[d]);
			if (cost < cmin)
			{
				cmin = cost;
				best_index = d;
				printf("best cost is now %1.6f at index %d\n", cmin, d);
			}
		}

		printf("minimum prac cost of prime %lu is %1.6f at index %d, val %1.18f\n",
			q, cmin, best_index, val[best_index]);

		pcost[i] = cmin;
		v[i] = val[best_index];

		i++;
		q = tecm_primes[i];

	}

	printf("now optimizing composites composed of combinations of 4 primes\n");
	for (i = 3; i < 19; i++)
	{
		for (j = i + 1; j < 19; j++)
		{
			for (k = j + 1; k < 19; k++)
			{
				for (l = k + 1; l < 19; l++)
				{
					double fcost;
					e = tecm_primes[i] * tecm_primes[j] * tecm_primes[k] * tecm_primes[l];

					printf("now optimizing composite %lu = %u * %u * %u * %u\n", 
						e, tecm_primes[i], tecm_primes[j], tecm_primes[k], tecm_primes[l]);
					//printf("now optimizing composite %lu = %u * %u * %u\n",
				//		e, tecm_primes[i], tecm_primes[j], tecm_primes[k]); 

					for (d = 0, cmin = ADD * (double)e; d < (num + NV); d++)
					{
						cost = lucas_cost(e, val[d]);
						if (cost < cmin)
						{
							cmin = cost;
							best_index = d;
							//printf("best cost is now %1.6f at index %d\n", cmin, d);
						}
					}

					fcost = pcost[i] + pcost[j] + pcost[k] + pcost[l];
					f[0] = tecm_primes[i];
					f[1] = tecm_primes[j];
					f[2] = tecm_primes[k];
					f[3] = tecm_primes[l];

					if (fcost > cmin)
					{
						FILE *fid;
						fid = fopen("lucas_cost_comb.txt", "a");
						printf("minimum prac cost of %lu is %1.1f at index %d, "
							"val %1.18f, factor cost of [", e, cmin, best_index, val[best_index]);
						printf("%lu,", f[0]);
						printf("%lu,", f[1]);
						printf("%lu,", f[2]);
						printf("%lu,", f[3]);
						printf("\b] is %1.1f\n", fcost);
						fflush(stdout);

						fprintf(fid, "minimum prac cost of %lu is %1.1f at index %d, "
							"val %1.18f, factor cost of [", e, cmin, best_index, val[best_index]);
						fprintf(fid, "%lu,", f[0]);
						fprintf(fid, "%lu,", f[1]);
						fprintf(fid, "%lu,", f[2]);
						fprintf(fid, "%lu,", f[3]);
						fprintf(fid, "\b] is %1.1f\n", fcost);
						fclose(fid);
					}
				}
			}
		}
	}

	exit(0);
	free(val);
	return;
}

void prac70(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P)
{
	uint64_t *s1, *s2, *d1, *d2;
	int i;
	static const uint8_t steps[116] = {
		0,6,0,6,0,6,0,4,6,0,4,6,0,4,4,6,
		0,4,4,6,0,5,4,6,0,3,3,4,6,0,3,5,
		4,6,0,3,4,3,4,6,0,5,5,4,6,0,5,3,
		3,4,6,0,3,3,4,3,4,6,0,5,3,3,3,3,
		3,3,3,3,4,3,3,4,6,0,5,4,3,3,4,6,
		0,3,4,3,5,4,6,0,5,3,3,3,4,6,0,5,
		4,3,5,4,6,0,5,5,3,3,4,6,0,4,3,3,
		3,5,4,6 };

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	for (i = 0; i < 116; i++)
	{
		if (steps[i] == 0)
		{
			work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
			work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
			work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
			work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

			// point2 is [2]P
			duplicate(mdata, work, s1, d1, &work->pt1);
		}
		else if (steps[i] == 3)
		{
			// integrate step 4 followed by swap(1,2)
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}
		else if (steps[i] == 4)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if (steps[i] == 5)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
		}
		else if (steps[i] == 6)
		{
			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
		}

	}

	return;

}

void prac85(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P)
{
	uint64_t *s1, *s2, *d1, *d2;
	int i;
	static const uint8_t steps[146] = {
		0,6,0,6,0,6,0,6,0,4,
		6,0,4,6,0,4,4,6,0,4,
		4,6,0,5,4,6,0,3,3,4,
		6,0,3,5,4,6,0,3,4,3,
		4,6,0,5,5,4,6,0,5,3,
		3,4,6,0,3,3,4,3,4,6,
		0,4,3,4,3,5,3,3,3,3,
		3,3,3,3,4,6,0,3,3,3,
		3,3,3,3,3,3,4,3,4,3,
		4,6,0,3,4,3,5,4,6,0,
		5,3,3,3,4,6,0,5,4,3,
		5,4,6,0,4,3,3,3,5,4,
		6,0,4,3,5,3,3,4,6,0,
		3,3,3,3,5,4,6,0,3,3,
		3,4,3,3,4,6 };

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	for (i = 0; i < 146; i++)
	{
		if (steps[i] == 0)
		{
			work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
			work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
			work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
			work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

			// point2 is [2]P
			duplicate(mdata, work, s1, d1, &work->pt1);
		}
		else if (steps[i] == 3)
		{
			// integrate step 4 followed by swap(1,2)
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}
		else if (steps[i] == 4)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if (steps[i] == 5)
		{
			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
		}
		else if (steps[i] == 6)
		{
			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
		}
	}

	return;

}

void prac_good(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P, uint64_t c, double v_in)
{
	uint64_t d, e, r;
	double cmin, cost, v;
	int i;
	uint64_t *s1, *s2, *d1, *d2;
	uint32_t *sw_x, *sw_z;
	int shift = 0;

	while ((c & 1) == 0)
	{
		shift++;
		c >>= 1;
	}

	/* 1/val[0] = the golden ratio (1+sqrt(5))/2, and 1/val[i] for i>0
	   is the real number whose continued fraction expansion is all 1s
	   except for a 2 in i+1-st place */
	if (v_in < 0.0001)
	{
		static double val[NV] =
		{ 0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
		  0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
		  0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
		  0.61807966846989581 };

		/* chooses the best value of v */
		for (d = 0, cmin = ADD * (double)c; d < NV; d++)
		{
			cost = lucas_cost(c, val[d]);
			if (cost < cmin)
			{
				cmin = cost;
				i = d;
			}
		}
		v = val[i];
	}
	else
	{
		v = v_in;
	}

	d = c;
	r = (uint64_t)((double)d * v + 0.5);

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	/* first iteration always begins by Condition 3, then a swap */
	d = c - r;
	if ((c >= (2 * r)) || (r >= c))
	{
		printf("problem\n");
		exit(1);
	}
	e = 2 * r - c;

	// mpres_set(xB, xA, n);
	// mpres_set(zB, zA, n); /* B=A */
	// mpres_set(xC, xA, n);
	// mpres_set(zC, zA, n); /* C=A */
	// duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

	// the first one is always a doubling
	// point1 is [1]P
	work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
	work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
	work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
	work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

	submod128(work->pt1.X, work->pt1.Z, d1, work->n);
	addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

	// point2 is [2]P
	duplicate(mdata, work, s1, d1, &work->pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			//mpres_swap(xA, xB, n);
			//mpres_swap(zA, zB, n);
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}
		/* do the first line of Table 4 whose condition qualifies */
		if (d - e <= e / 4 && ((d + e) % 3) == 0)
		{ /* condition 1 */
			d = (2 * d - e) / 3;
			e = (e - d) / 2;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4); // T = A + B (C)
			add(mdata, work, &work->pt4, &work->pt1, &work->pt2, &work->pt5); // T2 = T + A (B)
			add(mdata, work, &work->pt2, &work->pt4, &work->pt1, &work->pt2); // B = B + T (A)

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xT2, zT2, xT, zT, xA, zA, xB, zB, n, u, v, w); /* T2 = f(T,A,B) */
			//add3(xB, zB, xB, zB, xT, zT, xA, zA, n, u, v, w); /* B = f(B,T,A) */
			//mpres_swap(xA, xT2, n);
			//mpres_swap(zA, zT2, n); /* swap A and T2 */
			swap128(work->pt1.X, work->pt5.X);
			swap128(work->pt1.Z, work->pt5.Z);
		}
		else if (d - e <= e / 4 && (d - e) % 6 == 0)
		{ /* condition 2 */
			d = (d - e) / 2;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt2);		// B = A + B (C)

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);
			duplicate(mdata, work, s1, d1, &work->pt1);		// A = 2A

			//add3(xB, zB, xA, zA, xB, zB, xC, zC, n, u, v, w); /* B = f(A,B,C) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */

		}
		else if ((d + 3) / 4 <= e)
		{ /* condition 3 */
			d -= e;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			//add3(xT, zT, xB, zB, xA, zA, xC, zC, n, u, v, w); /* T = f(B,A,C) */

			/* circular permutation (B,T,C) */
			//tmp = xB;
			//xB = xT;
			//xT = xC;
			//xC = tmp;
			//tmp = zB;
			//zB = zT;
			//zT = zC;
			//zC = tmp;
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if ((d + e) % 2 == 0)
		{ /* condition 4 */
			d = (d - e) / 2;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A

			//add3(xB, zB, xB, zB, xA, zA, xC, zC, n, u, v, w); /* B = f(B,A,C) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d+e is odd */
		else if (d % 2 == 0)
		{ /* condition 5 */
			d /= 2;

			add(mdata, work, &work->pt3, &work->pt1, &work->pt2, &work->pt3);		// C = C + A (B)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A

			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(C,A,B) */
			//duplicate(xA, zA, xA, zA, n, b, u, v, w); /* A = 2*A */
		}
		/* now d is odd, e is even */
		else if (d % 3 == 0)
		{ /* condition 6 */
			d = d / 3 - e;

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);
			duplicate(mdata, work, s1, d1, &work->pt4);		// T = 2A

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt5);		// T2 = A + B (C)
			add(mdata, work, &work->pt4, &work->pt1, &work->pt1, &work->pt1);		// A = T + A (A)
			add(mdata, work, &work->pt4, &work->pt5, &work->pt3, &work->pt4);		// T = T + T2 (C)

			//duplicate(xT, zT, xA, zA, n, b, u, v, w); /* T = 2*A */
			//add3(xT2, zT2, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T2 = f(A,B,C) */
			//add3(xA, zA, xT, zT, xA, zA, xA, zA, n, u, v, w); /* A = f(T,A,A) */
			//add3(xT, zT, xT, zT, xT2, zT2, xC, zC, n, u, v, w); /* T = f(T,T2,C) */

			/* circular permutation (C,B,T) */
			//tmp = xC;
			//xC = xB;
			//xB = xT;
			//xT = tmp;
			//tmp = zC;
			//zC = zB;
			//zB = zT;
			//zT = tmp;
			rot128(work->pt3.X, work->pt2.X, work->pt4.X);
			rot128(work->pt3.Z, work->pt2.Z, work->pt4.Z);
		}
		else if ((d + e) % 3 == 0)
		{ /* condition 7 */
			d = (d - 2 * e) / 3;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4);		// T = A + B (C)
			add(mdata, work, &work->pt4, &work->pt1, &work->pt2, &work->pt2);		// B = T + A (B)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt4);		// T = 2A
			add(mdata, work, &work->pt1, &work->pt4, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xB, zB, xT, zT, xA, zA, xB, zB, n, u, v, w); /* B = f(T,A,B) */
			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
		else if ((d - e) % 3 == 0)
		{ /* condition 8 */
			d = (d - e) / 3;

			add(mdata, work, &work->pt1, &work->pt2, &work->pt3, &work->pt4);		// T = A + B (C)
			add(mdata, work, &work->pt3, &work->pt1, &work->pt2, &work->pt3);		// C = C + A (B)

			//add3(xT, zT, xA, zA, xB, zB, xC, zC, n, u, v, w); /* T = f(A,B,C) */
			//add3(xC, zC, xC, zC, xA, zA, xB, zB, n, u, v, w); /* C = f(A,C,B) */
			//mpres_swap(xB, xT, n);
			//mpres_swap(zB, zT, n); /* swap B and T */
			swap128(work->pt2.X, work->pt4.X);
			swap128(work->pt2.Z, work->pt4.Z);

			submod128(work->pt1.X, work->pt1.Z, d1, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s1, work->n);
			duplicate(mdata, work, s1, d1, &work->pt4);		// T = 2A
			add(mdata, work, &work->pt1, &work->pt4, &work->pt1, &work->pt1);		// A = A + T (A) = 3A

			//duplicate(xT, zT, xA, zA, n, b, u, v, w);
			//add3(xA, zA, xA, zA, xT, zT, xA, zA, n, u, v, w); /* A = 3*A */
		}
		else /* necessarily e is even here */
		{ /* condition 9 */
			e /= 2;

			add(mdata, work, &work->pt3, &work->pt2, &work->pt1, &work->pt3);		// C = C + B (A)

			submod128(work->pt2.X, work->pt2.Z, d2, work->n);
			addmod128(work->pt2.X, work->pt2.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt2);		// B = 2B

			//add3(xC, zC, xC, zC, xB, zB, xA, zA, n, u, v, w); /* C = f(C,B,A) */
			//duplicate(xB, zB, xB, zB, n, b, u, v, w); /* B = 2*B */
		}
	}

	add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)
	//add3(xA, zA, xA, zA, xB, zB, xC, zC, n, u, v, w);

	for (i = 0; i < shift; i++)
	{
		submod128(P->X, P->Z, d1, work->n);
		addmod128(P->X, P->Z, s1, work->n);
		duplicate(mdata, work, s1, d1, P);		// P = 2P
	}

	if (d != 1)
	{
		printf("problem: d != 1\n");
	}

	return;

}

void prac(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P, uint64_t c, double v)
{
	uint64_t d, e, r;
	int i;
	uint64_t *s1, *s2, *d1, *d2;
	uint64_t swp;

	d = c;
	r = (uint64_t)((double)d * v + 0.5);

	s1 = work->sum1;
	s2 = work->sum2;
	d1 = work->diff1;
	d2 = work->diff2;

	d = c - r;
	e = 2 * r - c;

	// the first one is always a doubling
	// point1 is [1]P
	work->pt1.X[0] = work->pt2.X[0] = work->pt3.X[0] = P->X[0];
	work->pt1.X[1] = work->pt2.X[1] = work->pt3.X[1] = P->X[1];
	work->pt1.Z[0] = work->pt2.Z[0] = work->pt3.Z[0] = P->Z[0];
	work->pt1.Z[1] = work->pt2.Z[1] = work->pt3.Z[1] = P->Z[1];

	submod128(work->pt1.X, work->pt1.Z, d1, work->n);
	addmod128(work->pt1.X, work->pt1.Z, s1, work->n);

	// point2 is [2]P
	duplicate(mdata, work, s1, d1, &work->pt1);

	while (d != e)
	{
		if (d < e)
		{
			r = d;
			d = e;
			e = r;
			swap128(work->pt1.X, work->pt2.X);
			swap128(work->pt1.Z, work->pt2.Z);
		}

		if ((d + 3) / 4 <= e)
		{
			d -= e;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt4);		// T = B + A (C)
			rot128(work->pt2.X, work->pt4.X, work->pt3.X);
			rot128(work->pt2.Z, work->pt4.Z, work->pt3.Z);
		}
		else if ((d + e) % 2 == 0)
		{
			d = (d - e) / 2;

			add(mdata, work, &work->pt2, &work->pt1, &work->pt3, &work->pt2);		// B = B + A (C)

			submod128(work->pt1.X, work->pt1.Z, d2, work->n);
			addmod128(work->pt1.X, work->pt1.Z, s2, work->n);
			duplicate(mdata, work, s2, d2, &work->pt1);		// A = 2A
		}
		else
		{
			// empirically, tiny B1 values only need the above prac cases.
			// just in case, fall back on this.
			printf("unhandled case in prac\n");
			exit(1);
		}
	}

	add(mdata, work, &work->pt1, &work->pt2, &work->pt3, P);		// A = A + B (C)

	return;

}

void build_one_curve(tinyecm_pt *P, monty128_t *mdata, 
	tinyecm_work *work, uint32_t sigma, uint64_t * lcg_state, int verbose)
{
	base_t t1[2], t2[2], t3[2], t4[2], t5[2], s[3];
	base_t u[2], v[2], n[2];
	mpz_t gmpt, gmpn;
	mpz_init(gmpt);
	mpz_init(gmpn);

	n[0] = mdata->n[0];
	n[1] = mdata->n[1];

	if (verbose)
		printf("n = %016lx%016lx\n", n[1], n[0]);

	if (sigma == 0)
	{
		work->sigma = lcg_rand_32_range(7, (uint32_t)-1, lcg_state);
	}
	else
	{
		work->sigma = sigma;
	}
	sigma = work->sigma;

	u[0] = sigma;
	u[1] = 0;
	to_monty128(mdata, u);
	
	if (verbose)
		printf("monty(sigma) = %016lx%016lx\n", u[1], u[0]);

	t1[0] = 4;
	t1[1] = 0;
	to_monty128(mdata, t1);

	if (verbose)
		printf("monty(4) = %016lx%016lx\n", t1[1], t1[0]);

	mulmod128(u, t1, v, mdata);		// v = 4*sigma

	if (verbose)
		printf("v = 4*sigma = %016lx%016lx\n", v[1], v[0]);

	sqrmod128(u, u, mdata);
	t1[0] = 5;
	t1[1] = 0;
	to_monty128(mdata, t1);
	submod128(u, t1, u, mdata->n);		// u = sigma^2 - 5

	if (verbose)
		printf("u = sigma^2 - 5 = %016lx%016lx\n", u[1], u[0]);

	sqrmod128(u, t1, mdata);
	mulmod128(t1, u, P->X, mdata);	// x = u^3

	if (verbose)
		printf("x = u^3 = %016lx%016lx\n", P->X[1], P->X[0]);

	sqrmod128(v, t1, mdata);
	mulmod128(t1, v, P->Z, mdata);	// z = v^3

	if (verbose)
		printf("z = v^3 = %016lx%016lx\n", P->Z[1], P->Z[0]);

	//compute parameter A
	submod128(v, u, t1, mdata->n);		// (v - u)
	sqrmod128(t1, t2, mdata);
	mulmod128(t2, t1, t4, mdata);	// (v - u)^3

	if (verbose)
		printf("(v - u)^3 = %016lx%016lx\n", t4[1], t4[0]);

	t1[0] = 3;
	t1[1] = 0;
	to_monty128(mdata, t1);
	mulmod128(t1, u, t2, mdata);		// 3u
	addmod128(t2, v, t3, mdata->n);		// 3u + v

	if (verbose)
		printf("3u + v = %016lx%016lx\n", t3[1], t3[0]);

	mulmod128(t3, t4, t1, mdata);	// a = (v-u)^3 * (3u + v)

	if (verbose)
		printf("a = (v-u)^3 * (3u + v) = %016lx%016lx\n", t1[1], t1[0]);

	t2[0] = 16;
	t2[1] = 0;
	to_monty128(mdata, t2);
	mulmod128(P->X, t2, t3, mdata);	// 16*u^3
	mulmod128(t3, v, t4, mdata);		// 16*u^3*v

	if (verbose)
		printf("16*u^3*v = %016lx%016lx\n", t4[1], t4[0]);

	// u holds the denom, t1 holds the numer
	// accomplish the division by multiplying by the modular inverse
	t2[0] = 1;
	t2[1] = 0;
	mulmod128(t4, t2, t4, mdata);	// take t4 out of monty rep
	mulmod128(t1, t2, t1, mdata);	// take t1 out of monty rep

	u128_to_mpz(t4, gmpt);
	u128_to_mpz(n, gmpn);
	mpz_invert(gmpt, gmpt, gmpn);		// gmpt = t4^-1 mod n
	
	mpz_to_u128(gmpt, t3);
	
	if (verbose)
		printf("1/16*u^3*v = %016lx%016lx\n", t3[1], t3[0]);

	u128_to_mpz(t1, gmpn);
	mpz_mul(gmpt, gmpt, gmpn);			// gmpt = t1 * t4 
	u128_to_mpz(n, gmpn);
	mpz_tdiv_r(gmpt, gmpt, gmpn);		// gmpt = t1 * t4 % n
	mpz_to_u128(gmpt, work->s);
	to_monty128(mdata, work->s);

	u128_to_mpz(t4, gmpn);
	u128_to_mpz(t3, gmpt);
	mpz_mul(gmpt, gmpt, gmpn);
	u128_to_mpz(n, gmpn);
	mpz_tdiv_r(gmpt, gmpt, gmpn);

	//if (mpz_cmp_ui(gmpt, 1) != 0)
	//{
	//	gmp_printf("inversion produced result %Zd from %Zd\n", gmpt, gmpn);
	//	//exit(1);
	//}

	
	// t1 = b = (v - u)^3 * (3*u + v) / 16u^3v
	//mulmod128(t3, t1, work->s, s, mdata);

	mpz_clear(gmpt);
	mpz_clear(gmpn);
	return;
}

void tinyecm(mpz_t n, mpz_t f, uint32_t B1, uint32_t B2, uint32_t curves,
    uint64_t* lcg_state, int verbose)
{
	//attempt to factor n with the elliptic curve method
	//following brent and montgomery's papers, and CP's book
	base_t retval;
	base_t i, j;
	int curve;
	int tid;
	char *wstr;
	int found = 0;
	int result;
	uint64_t num_found;
	tinyecm_work work;
	tinyecm_pt P;
	monty128_t mdata;
	uint64_t n128[2];
	uint32_t sigma;

	mpz_to_u128(n, n128);
	tinyecm_work_init(&work);
	monty128_init(&mdata, n128);
	copy128(n128, work.n);
	work.stg1_max = B1;
	work.stg2_max = B2;

    if (!tecm_primes_initialized)
    {
        soe_staticdata_t* sdata = soe_init(0, 1, 32768);
        tecm_primes = soe_wrapper(sdata, 0, 1000000, 0, &tecm_nump, 0, 0);
        tecm_maxp = tecm_primes[tecm_nump - 1];
        soe_finalize(sdata);
        tecm_primes_initialized = 1;
    }

	if (0)
		lucas_opt2(B1, 1000000);

	if (0)
	{
		mpz_t gmp1, gmp2, gmp3, gmp4;
		int iterations = 10000;
        

		mpz_init(gmp1);
		mpz_init(gmp2);
		mpz_init(gmp3);
		mpz_init(gmp4);

		printf("commencing test of modular arithmetic mod %016lx%016lx\n", n128[1], n128[0]);

		// addmod
		for (i = 0; i < iterations; i++)
		{
			uint64_t a[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t b[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t c[2];

			u128_to_mpz(a, gmp1);
			u128_to_mpz(b, gmp2);
			mpz_mod(gmp1, gmp1, n);
			mpz_mod(gmp2, gmp2, n);
			mpz_to_u128(gmp1, a);
			mpz_to_u128(gmp2, b);
			mpz_add(gmp3, gmp1, gmp2);
			mpz_mod(gmp3, gmp3, n);
			addmod128(a, b, c, n128);
			u128_to_mpz(c, gmp4);
			if (mpz_cmp(gmp3, gmp4) != 0)
			{
				printf("addmod failure:\n");
				printf("a = %016lx%016lx\n", a[1], a[0]);
				printf("b = %016lx%016lx\n", b[1], b[0]);
				printf("c = %016lx%016lx\n", c[1], c[0]);
				printf("n = %016lx%016lx\n", n128[1], n128[0]);
				gmp_printf("gmp result = %Zx\n", gmp3);
			}

		}

		// submod
		for (i = 0; i < iterations; i++)
		{
			uint64_t a[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t b[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t c[2];

			u128_to_mpz(a, gmp1);
			u128_to_mpz(b, gmp2);
			mpz_mod(gmp1, gmp1, n);
			mpz_mod(gmp2, gmp2, n);
			mpz_to_u128(gmp1, a);
			mpz_to_u128(gmp2, b);
			mpz_sub(gmp3, gmp1, gmp2);
			if (mpz_sgn(gmp3) < 0)
				mpz_add(gmp3, gmp3, n);
			submod128(a, b, c, n128);
			u128_to_mpz(c, gmp4);
			if (mpz_cmp(gmp3, gmp4) != 0)
			{
				printf("submod failure:\n");
				printf("a = %016lx%016lx\n", a[1], a[0]);
				printf("b = %016lx%016lx\n", b[1], b[0]);
				printf("c = %016lx%016lx\n", c[1], c[0]);
				printf("n = %016lx%016lx\n", n128[1], n128[0]);
				gmp_printf("gmp result = %Zx\n", gmp3);
			}

		}

		// mulmod
		for (i = 0; i < iterations; i++)
		{
			uint64_t a[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t b[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t c[2];
			uint64_t s[3];

			u128_to_mpz(a, gmp1);
			u128_to_mpz(b, gmp2);
			mpz_mod(gmp1, gmp1, n);
			mpz_mod(gmp2, gmp2, n);
			mpz_to_u128(gmp1, a);
			mpz_to_u128(gmp2, b);
			mpz_mul(gmp3, gmp1, gmp2);
			mpz_mod(gmp3, gmp3, n);
			mpz_mul_2exp(gmp3, gmp3, 128);
			mpz_mod(gmp3, gmp3, n);
			to_monty128(&mdata, a);
			to_monty128(&mdata, b);
			mulmod128(a, b, c, &mdata);
			u128_to_mpz(c, gmp4);
			if (mpz_cmp(gmp3, gmp4) != 0)
			{
				printf("mulmod failure:\n");
				printf("a = %016lx%016lx\n", a[1], a[0]);
				printf("b = %016lx%016lx\n", b[1], b[0]);
				printf("c = %016lx%016lx\n", c[1], c[0]);
				printf("n = %016lx%016lx\n", n128[1], n128[0]);
				gmp_printf("gmp result = %Zx\n", gmp3);
			}

		}

		// sqrmod
		for (i = 0; i < iterations; i++)
		{
			uint64_t a[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t b[2] = { lcg_rand_64(lcg_state), lcg_rand_64(lcg_state) };
			uint64_t c[2];
			uint64_t s[3];

			u128_to_mpz(a, gmp1);
			u128_to_mpz(b, gmp2);
			mpz_mod(gmp1, gmp1, n);
			mpz_mod(gmp2, gmp2, n);
			mpz_to_u128(gmp1, a);
			mpz_to_u128(gmp2, b);
			mpz_mul(gmp3, gmp1, gmp1);
			mpz_mod(gmp3, gmp3, n);
			mpz_mul_2exp(gmp3, gmp3, 128);
			mpz_mod(gmp3, gmp3, n);
			to_monty128(&mdata, a);
			sqrmod128(a, c, &mdata);
			u128_to_mpz(c, gmp4);
			if (mpz_cmp(gmp3, gmp4) != 0)
			{
				printf("sqrmod failure:\n");
				printf("a = %016lx%016lx\n", a[1], a[0]);
				printf("c = %016lx%016lx\n", c[1], c[0]);
				printf("n = %016lx%016lx\n", n128[1], n128[0]);
				gmp_printf("gmp result = %Zx\n", gmp3);
			}

		}

		mpz_clear(gmp1);
		mpz_clear(gmp2);
		mpz_clear(gmp3);
		mpz_clear(gmp4);
	}

	mpz_set_ui(f, 1);
	for (curve = 0; curve < curves; curve++)
	{
		work.stg1Add = 0;
		work.stg1Doub = 0;
		work.last_pid = 0;

		if (verbose)
			printf("commencing curve %d of %u\n", curve, curves);

		sigma = 0;
		build_one_curve(&P, &mdata, &work, sigma, lcg_state, verbose);

		if (verbose)
		{
			printf("curve parameters:\n\tsigma = %u\n", work.sigma);
			printf("\tn = %016lx%016lx\n", work.n[1], work.n[0]);
			printf("\tx = %016lx%016lx\n", P.X[1], P.X[0]);
			printf("\tz = %016lx%016lx\n", P.Z[1], P.Z[0]);
			printf("\tb = %016lx%016lx\n", work.s[1], work.s[0]);
		}

		ecm_stage1(&mdata, &work, &P);

		result = check_factor(P.Z, mdata.n, mdata.mtmp1);
			
		if (result == 1)
		{
			if (verbose)
				printf("\nfound factor %016lx%016lx in stage 1 with sigma = %u\n",
					mdata.mtmp1[1], mdata.mtmp1[0], sigma);

			u128_to_mpz(mdata.mtmp1, f);
			break;
		}

		if (B2 > B1)
		{
			ecm_stage2(&P, &mdata, &work);
			result = check_factor(work.stg2acc, mdata.n, mdata.mtmp1);

			if (result == 1)
			{
				if (verbose)
					printf("\nfound factor %016lx%016lx in stage 2 with sigma = %u\n",
						mdata.mtmp1[1], mdata.mtmp1[0], sigma);

				u128_to_mpz(mdata.mtmp1, f);
				break;
			}
		}
	}

	return;
}

//#define TESTMUL
void ecm_stage1_good(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P, 
	base_t b1, base_t *primes, int verbose)
{
	int i;
	uint64_t q;
	uint64_t stg1 = (uint64_t)work->stg1_max;
	 
	// handle the only even case 
	q = 2;
	while (q < stg1)
	{
		submod128(P->X, P->Z, work->diff1, work->n);
		addmod128(P->X, P->Z, work->sum1, work->n);
		duplicate(mdata, work, work->sum1, work->diff1, P);
		q *= 2;
	}

	for (i = 1; (i < tecm_nump) && ((uint32_t)tecm_primes[i] < stg1); i++)
	{
		uint64_t c = 1;

		q = tecm_primes[i];
		do {
			prac(mdata, work, P, q, 0);
			c *= q;
		} while ((c * q) < stg1);
	}

	work->last_pid = i;

	if (verbose == 1)
	{
		printf("\nStage 1 completed at prime %lu with %u point-adds and %u point-doubles\n",
            tecm_primes[i - 1], work->stg1Add, work->stg1Doub);
		fflush(stdout);
	}
	return;
}

void ecm_stage1(monty128_t *mdata, tinyecm_work *work, tinyecm_pt *P)
{
	int i;
	uint64_t q;
	uint64_t stg1 = (uint64_t)work->stg1_max;

	// handle the only even case 
	q = 2;
	while (q < stg1 * 4)  // jeff: multiplying by 4 improves perf ~1%
	{
		submod128(P->X, P->Z, work->diff1, work->n);
		addmod128(P->X, P->Z, work->sum1, work->n);
		duplicate(mdata, work, work->sum1, work->diff1, P);
		q *= 2;
	}

	if (stg1 == 27)
	{
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 11, 0.580178728295464130);
		prac(mdata, work, P, 13, 0.618033988749894903);
		prac(mdata, work, P, 17, 0.618033988749894903);
		prac(mdata, work, P, 19, 0.618033988749894903);
		prac(mdata, work, P, 23, 0.522786351415446049);
	}
	else if (stg1 == 47)
	{
		// jeff: improved perf slightly by using one more uprac for 3,
		// and removing uprac for 47.
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 3, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 11, 0.580178728295464130);
		prac(mdata, work, P, 13, 0.618033988749894903);
		prac(mdata, work, P, 17, 0.618033988749894903);
		prac(mdata, work, P, 19, 0.618033988749894903);
		prac(mdata, work, P, 23, 0.522786351415446049);
		prac(mdata, work, P, 29, 0.548409048446403258);
		prac(mdata, work, P, 31, 0.618033988749894903);
		prac(mdata, work, P, 37, 0.580178728295464130);
		prac(mdata, work, P, 41, 0.548409048446403258);
		prac(mdata, work, P, 43, 0.618033988749894903);
		//        uecm_uprac(mdata, work, P, 47, 0.548409048446403258);
	}
	else if (stg1 == 59)
	{   // jeff: probably stg1 of 59 would benefit from similar changes
		// as stg1 of 47 above, but I didn't bother. Stg1 of 59 seems to
		// always perform worse than stg1 of 47, so there doesn't seem
		// to be any reason to ever use stg1 of 59.
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 3, 0.61803398874989485);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 5, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 7, 0.618033988749894903);
		prac(mdata, work, P, 11, 0.580178728295464130);
		prac(mdata, work, P, 13, 0.618033988749894903);
		prac(mdata, work, P, 17, 0.618033988749894903);
		prac(mdata, work, P, 19, 0.618033988749894903);
		prac(mdata, work, P, 23, 0.522786351415446049);
		prac(mdata, work, P, 29, 0.548409048446403258);
		prac(mdata, work, P, 31, 0.618033988749894903);
		prac(mdata, work, P, 1961, 0.552936068843375);   // 37 * 53
		prac(mdata, work, P, 41, 0.548409048446403258);
		prac(mdata, work, P, 43, 0.618033988749894903);
		prac(mdata, work, P, 47, 0.548409048446403258);
		prac(mdata, work, P, 59, 0.548409048446403258);
	}
	else if (stg1 == 70)
	{
		prac70(mdata, work, P);
		i = 19;
	}
	else // if (stg1 >= 85)
	{
		prac85(mdata, work, P);

		if (stg1 == 85)
		{
			prac(mdata, work, P, 61, 0.522786351415446049);
		}
		else
		{
			prac(mdata, work, P, 5, 0.618033988749894903);
			prac(mdata, work, P, 11, 0.580178728295464130);
			//            uecm_uprac(mdata, work, P, 61, 0.522786351415446049);
			prac(mdata, work, P, 89, 0.618033988749894903);
			prac(mdata, work, P, 97, 0.723606797749978936);
			prac(mdata, work, P, 101, 0.556250337855490828);
			prac(mdata, work, P, 107, 0.580178728295464130);
			prac(mdata, work, P, 109, 0.548409048446403258);
			prac(mdata, work, P, 113, 0.618033988749894903);

			if (stg1 == 125)
			{
				// jeff: moved 61 to here
				prac(mdata, work, P, 61, 0.522786351415446049);
				prac(mdata, work, P, 103, 0.632839806088706269);
			}
			else
			{
				prac(mdata, work, P, 7747, 0.552188778811121); // 61 x 127
				prac(mdata, work, P, 131, 0.618033988749894903);
				prac(mdata, work, P, 14111, 0.632839806088706);  // 103 x 137
				prac(mdata, work, P, 20989, 0.620181980807415);  // 139 x 151
				prac(mdata, work, P, 157, 0.640157392785047019);
				prac(mdata, work, P, 163, 0.551390822543526449);

				if (stg1 == 165)
				{
					prac(mdata, work, P, 149, 0.580178728295464130);
				}
				else
				{
					prac(mdata, work, P, 13, 0.618033988749894903);
					prac(mdata, work, P, 167, 0.580178728295464130);
					prac(mdata, work, P, 173, 0.612429949509495031);
					prac(mdata, work, P, 179, 0.618033988749894903);
					prac(mdata, work, P, 181, 0.551390822543526449);
					prac(mdata, work, P, 191, 0.618033988749894903);
					prac(mdata, work, P, 193, 0.618033988749894903);
					prac(mdata, work, P, 29353, 0.580178728295464);  // 149 x 197
					prac(mdata, work, P, 199, 0.551390822543526449);
				}
			}
		}
	}
	return;
}

// pre-paired sequences for various B1 and B2 = 25*B1
static const int numb1_27 = 76;
static const uint8_t b1_27[76] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7 };

static const int numb1_47 = 121;
static const uint8_t b1_47[121] = {
1,7,13,17,23,29,31,41,43,47,49,0,59,47,43,41,
37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,
23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,
11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,
11,1,47,0,59,49,41,31,23,17,11,7,1,19,37,47,
0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,
37,23,19,13,7,1,29,31,41,59,0,59,49,47,41,23,
19,17,13,7,1,43,53,0,59 };

static const int numb1_70 = 175;
static const uint8_t b1_70[175] = {
31,41,43,47,49,0,59,47,43,41,37,31,29,19,13,7,
1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,
0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,
41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,
31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,
29,31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,
53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,
53,49,47,43,31,29,23,11,17,0,47,43,41,37,31,23,
19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,
7,17,29,0,53,47,43,41,17,13,11,1,23,31,37 };

static const int numb1_85 = 225;
static const uint8_t b1_85[225] = {
	1,53,49,47,43,41,37,23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,
	31,23,17,11,7,1,19,29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,
	0,59,49,41,31,23,17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,
	31,41,59,0,59,49,47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,
	29,23,11,17,0,47,43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,
	17,13,11,1,23,31,37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,
	43,41,31,17,7,1,11,13,19,29 };

static const int numb1_125 = 319;
static const uint8_t b1_125[319] = {
	23,19,13,11,1,7,17,29,31,0,59,47,43,41,37,31,29,19,13,7,1,11,23,0,59,53,43,41,37,31,23,17,11,7,1,19,
	29,49,0,53,49,47,43,31,23,19,11,7,1,13,37,59,0,59,53,43,37,31,29,23,17,13,11,1,47,0,59,49,41,31,23,
	17,11,7,1,19,37,47,0,59,49,47,43,41,31,17,13,11,7,37,0,53,49,43,37,23,19,13,7,1,29,31,41,59,0,59,49,
	47,41,23,19,17,13,7,1,43,53,0,59,49,43,37,29,17,13,7,1,19,47,53,0,59,53,49,47,43,31,29,23,11,17,0,47,
	43,41,37,31,23,19,17,11,1,13,29,53,0,59,47,41,37,31,23,19,11,7,17,29,0,53,47,43,41,17,13,11,1,23,31,
	37,49,0,53,47,43,41,29,19,7,1,17,31,37,49,59,0,49,43,37,19,17,1,23,29,47,53,0,59,53,43,41,31,17,7,1,
	11,13,19,29,0,59,53,49,47,37,29,11,13,17,23,31,0,59,43,41,37,29,23,17,13,1,31,47,0,59,53,49,47,41,37,
	31,19,13,7,11,17,29,43,0,47,29,19,11,7,1,41,43,59,0,53,49,37,23,13,11,7,1,17,19,29,41,43,59,0,59,49,
	41,37,23,13,1,7,11,29,43,47,53,0,59,53,49,31,23,13,7,1,17,29,43,47,0,59,31,29,19,11,7,37,49,53 };

static const int numb1_165 = 425;
static const uint8_t b1_165[425] = {
	13,7,1,11,19,47,59,0,59,49,43,37,31,29,23,19,17,7,11,13,47,53,0,53,47,41,37,31,23,19,11,1,13,29,43,
	59,0,53,49,41,37,31,19,17,1,7,23,29,47,59,0,59,53,47,43,41,29,19,17,13,7,1,23,31,49,0,53,47,41,37,29,
	23,19,11,7,17,31,43,49,59,0,47,43,41,37,23,19,17,13,7,11,29,53,0,53,49,43,37,29,23,11,7,1,13,19,31,41,
	0,53,49,47,43,37,31,23,17,11,13,41,0,59,47,43,37,31,29,23,11,1,17,19,41,0,59,53,19,13,7,1,29,43,47,49,
	0,53,49,47,41,29,19,17,13,11,7,1,23,31,43,59,0,53,49,41,37,23,19,13,11,7,1,17,43,47,0,47,43,41,31,19,
	17,7,1,13,37,49,0,59,49,37,29,13,1,7,11,17,19,41,47,53,0,49,47,31,29,7,1,13,17,19,23,37,59,0,47,37,31,
	19,17,13,11,1,29,41,43,53,0,59,41,17,13,7,1,19,23,31,47,49,53,0,59,53,47,43,31,29,7,1,11,17,37,41,49,
	0,49,43,37,23,19,13,1,7,17,0,59,49,41,37,31,29,23,1,11,13,53,0,53,43,41,37,29,23,17,13,11,7,1,19,31,49,
	0,53,43,31,29,23,19,17,1,13,37,41,59,0,53,43,37,31,23,13,1,17,29,59,0,59,49,41,37,23,19,11,1,7,29,0,59,
	43,17,13,11,1,7,23,29,37,41,49,0,49,47,43,41,29,1,7,13,19,23,31,59,0,59,49,47,31,29,13,7,37,41,43,0,49,
	41,29,23,13,11,7,1,17,19,31,43,53,0,53,47,43,37,29,23,17,1,11,13,31,41,49,59,0,53,47,41,19,13,11,1,17,
	23,43,0,53,49,47,37,23,19,11,7,17,29,31,43,0,53,31,19,17,13,7,1,29,37,59 };

static const int numb1_205 = 511;
static const uint8_t b1_205[511] = {
	1,23,41,0,59,53,49,47,37,23,19,17,13,1,7,29,43,0,53,49,41,31,29,19,17,11,7,1,13,37,59,0,49,47,29,23,
	13,7,1,17,31,37,43,0,59,49,47,43,37,31,29,17,13,7,1,11,19,53,0,59,53,49,41,37,23,13,1,11,17,19,29,43,
	47,0,53,49,47,43,23,19,11,1,7,17,37,41,0,59,53,41,37,31,29,19,17,11,1,13,43,47,0,53,47,41,19,17,7,1,
	11,23,31,43,59,0,59,53,41,31,13,11,7,1,17,29,37,0,49,43,37,29,11,1,13,17,19,23,41,0,59,49,47,43,41,37,
	31,19,7,1,13,23,29,53,0,53,49,43,41,37,31,29,23,13,7,17,19,47,59,0,49,47,37,29,23,17,11,7,13,19,31,41,
	53,0,59,43,29,23,19,17,13,11,1,41,0,59,37,31,23,17,13,11,7,1,19,29,43,53,0,49,47,43,41,31,19,17,1,7,11,
	13,23,0,47,43,37,29,13,11,7,1,17,19,23,31,59,0,59,37,31,29,23,19,13,1,7,11,41,47,53,0,53,49,43,31,23,
	17,13,41,59,0,59,53,31,19,17,1,7,11,23,37,47,49,0,59,53,47,43,41,37,31,23,19,17,11,1,0,59,53,49,47,31,
	17,13,7,1,11,29,37,0,53,43,31,17,13,7,1,29,41,49,0,53,49,41,29,23,11,7,1,19,31,47,0,47,43,41,29,23,19,
	7,1,11,49,0,59,31,29,23,17,11,7,1,13,41,43,0,59,43,37,17,1,7,11,13,19,41,49,0,59,53,43,41,37,31,29,23,
	13,11,1,47,0,59,53,47,31,19,17,13,1,7,11,29,37,43,49,0,49,43,41,31,17,13,7,11,23,37,53,0,53,49,41,23,
	19,13,11,7,1,17,37,59,0,49,47,43,37,31,29,23,1,7,41,0,59,43,41,37,31,17,13,11,7,47,49,0,59,49,47,37,31,
	29,19,17,7,1,0,53,47,37,19,13,1,11,31,41,0,49,47,37,23,17,13,11,7,19,31,53,0,59,53,47,29,13,11,7,1,23,
	41,0,49,47,41,37,19,11,13,17,23,29,31,43,0,59,29,19,13,1,41,43,47,53,0,59,53,43,41,37,23,17,11,7,1,13,
	29,49 };


void ecm_stage2(tinyecm_pt* P, monty128_t* mdata, tinyecm_work* work )
{
	int b;
	int i, j, k;
	tinyecm_pt Pa1;
	tinyecm_pt* Pa = &Pa1;
	tinyecm_pt Pb[18];
	tinyecm_pt Pd1;
	tinyecm_pt* Pd = &Pd1;
	const uint8_t* barray = 0;
	int numb;
	uint32_t stg1_max = work->stg1_max;

#ifdef MICRO_ECM_VERBOSE_PRINTF
	ptadds = 0;
	stg1Doub = 0;
	stg1Add = 0;
#endif

	// this function has been written for MICRO_ECM_PARAM_D of 60, so you
	// probably don't want to change it.
	const int MICRO_ECM_PARAM_D = 60;

	//stage 2 init
	//Q = P = result of stage 1
	//compute [d]Q for 0 < d <= MICRO_ECM_PARAM_D

	uint64_t Pbprod[18][2];

	// [1]Q
	//Pb[1] = *P;
	copy128(P->X, Pb[1].X);
	copy128(P->Z, Pb[1].Z);
	mulmod128(Pb[1].X, Pb[1].Z, Pbprod[1], mdata);

	// [2]Q
	submod128(P->X, P->Z, work->diff1, work->n);
	addmod128(P->X, P->Z, work->sum1, work->n);
	duplicate(mdata, work, work->sum1, work->diff1, &Pb[2]);

		/*
		Let D = MICRO_ECM_PARAM_D.

		D is small in tinyecm, so it is straightforward to just enumerate the needed
		points.  We can do it efficiently with two progressions mod 6.
		Pb[0] = scratch
		Pb[1] = [1]Q;
		Pb[2] = [2]Q;
		Pb[3] = [7]Q;   prog2
		Pb[4] = [11]Q;  prog1
		Pb[5] = [13]Q;  prog2
		Pb[6] = [17]Q;  prog1
		Pb[7] = [19]Q;  prog2
		Pb[8] = [23]Q;  prog1
		Pb[9] = [29]Q;  prog1
		Pb[10] = [30 == D]Q;   // shouldn't this be [31]Q?
		Pb[11] = [31]Q; prog2   // shouldn't this be [37]Q?
		Pb[12] = [37]Q; prog2   // shouldn't this be [41]Q?
		Pb[13] = [41]Q; prog1   // [43]Q?
		Pb[14] = [43]Q; prog2   // [47]Q?
		Pb[15] = [47]Q; prog1   // [49]Q?
		Pb[16] = [49]Q; prog2   // [53]Q?
		Pb[17] = [53]Q; prog1   // [59]Q?
		// Pb[18] = [59]Q; prog1   // [60]Q?   note: we can remove this line I believe.  Pb[18] is never set, and never used.  Therefore I changed the definition of Pb above to have only 18 elements.

		two progressions with total of 17 adds to get 15 values of Pb.
		6 + 5(1) -> 11 + 6(5) -> 17 + 6(11) -> 23 + 6(17) -> 29 + 6(23) -> 35 + 6(29) -> 41 + 6(35) -> 47 + 6(41) -> 53 + 6(47) -> 59
		6 + 1(5) -> 7  + 6(1) -> 13 + 6(7)  -> 19 + 6(13) -> 25 + 6(19) -> 31 + 6(25) -> 37 + 6(31) -> 43 + 6(37) -> 49

		we also need [2D]Q = [60]Q
		to get [60]Q we just need one more add:
		compute [60]Q from [31]Q + [29]Q([2]Q), all of which we
		have after the above progressions are computed.

		we also need [A]Q = [((B1 + D) / (2D) * 2 D]Q
		which is equal to the following for various common B1:
		B1      [x]Q
		65      [120]Q
		85      [120]Q
		125     [180]Q      note: according to the A[Q] formula above, wouldn't this be [120]Q?  ( I suspect maybe the formula is supposed to be [((B1 + D) / D) * D]Q )
		165     [180]Q      note: same as above.
		205     [240]Q

		and we need [x-D]Q as well, for the above [x]Q.
		So far we are getting [x]Q and [x-2D]Q each from prac(x,Q).
		There is a better way using progressions of [2D]Q
		[120]Q = 2*[60]Q
		[180]Q = [120]Q + [60]Q([60]Q)
		[240]Q = 2*[120]Q
		[300]Q = [240]Q + [60]Q([180]Q)
		...
		etc.

		*/

	tinyecm_pt pt5, pt6;

	// Calculate all Pb: the following is specialized for MICRO_ECM_PARAM_D=60
	// [2]Q + [1]Q([1]Q) = [3]Q
	add(mdata, work, &Pb[1], &Pb[2], &Pb[1], &Pb[3]);        // <-- temporary

	// 2*[3]Q = [6]Q
	submod128(Pb[3].X, Pb[3].Z, work->diff1, work->n);
	addmod128(Pb[3].X, Pb[3].Z, work->sum1, work->n);
	duplicate(mdata, work, work->sum1, work->diff1, &pt6);   // pt6 = [6]Q

	// [3]Q + [2]Q([1]Q) = [5]Q
	add(mdata, work, &Pb[3], &Pb[2], &Pb[1], &pt5);    // <-- pt5 = [5]Q
	//Pb[3] = pt5;
	copy128(pt5.X, Pb[3].X);
	copy128(pt5.Z, Pb[3].Z);

	// [6]Q + [5]Q([1]Q) = [11]Q
	add(mdata, work, &pt6, &pt5, &Pb[1], &Pb[4]);    // <-- [11]Q

	i = 3;
	k = 4;
	j = 5;
	while ((j + 12) < MICRO_ECM_PARAM_D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		add(mdata, work, &pt6, &Pb[k], &Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// [6]Q + [1]Q([5]Q) = [7]Q
	add(mdata, work, &pt6, &Pb[1], &pt5, &Pb[3]);    // <-- [7]Q
	i = 1;
	k = 3;
	j = 1;
	while ((j + 12) < MICRO_ECM_PARAM_D)
	{
		// [j+6]Q + [6]Q([j]Q) = [j+12]Q
		add(mdata, work, &pt6, &Pb[k], &Pb[i], &Pb[map[j + 12]]);
		i = k;
		k = map[j + 12];
		j += 6;
	}

	// Pd = [2w]Q
	// [31]Q + [29]Q([2]Q) = [60]Q
	add(mdata, work, &Pb[9], &Pb[10], &Pb[2], Pd);   // <-- [60]Q

#ifdef MICRO_ECM_VERBOSE_PRINTF
	ptadds++;
#endif

	// temporary - make [4]Q
	tinyecm_pt pt4;
	submod128(Pb[2].X, Pb[2].Z, work->diff1, work->n);
	addmod128(Pb[2].X, Pb[2].Z, work->sum1, work->n);
	duplicate(mdata, work, work->sum1, work->diff1, &pt4);   // pt4 = [4]Q


	// make all of the Pbprod's
	for (i = 3; i < 18; i++)
	{
		 mulmod128(Pb[i].X, Pb[i].Z, Pbprod[i], mdata);
	}

	//initialize info needed for giant step
	tinyecm_pt Pad;

	// Pd = [w]Q
	// [17]Q + [13]Q([4]Q) = [30]Q
	add(mdata, work, &Pb[map[17]], &Pb[map[13]], &pt4, &Pad);    // <-- [30]Q

	// [60]Q + [30]Q([30]Q) = [90]Q
	add(mdata, work, Pd, &Pad, &Pad, Pa);

	tinyecm_pt pt90;   // set pt90 = [90]Q
	tinyecm_pt pt60;   // set pt60 = [60]Q

	copy128(Pa->X, pt90.X);
	copy128(Pd->X, pt60.X);
	copy128(Pa->Z, pt90.Z);
	copy128(Pd->Z, pt60.Z);

	// [90]Q + [30]Q([60]Q) = [120]Q
	add(mdata, work, Pa, &Pad, Pd, Pd);

	// [120]Q + [30]Q([90]Q) = [150]Q
	add(mdata, work, Pd, &Pad, Pa, Pa);

	//initialize accumulator
	uint64_t acc[2];
	copy128(mdata->one, acc);

	// adjustment of Pa and Pad for particular B1.
	// Currently we have Pa=150, Pd=120, Pad=30

	if (stg1_max < 70)
	{
		// first process the appropriate b's with A=90
		static const int steps27[16] = { 59,53,49,47,43,37,31,29,23,19,17,11,7,1,13,41 };
		static const int steps47[15] = { 43,37,31,29,23,19,17,11,7,1,13,41,47,49,59 };
		const int* steps;
		int numsteps;
		if (stg1_max == 27)
		{
			steps = steps27;
			numsteps = 16;
		}
		else // if (stg1_max == 47)
		{
			steps = steps47;
			numsteps = 15;
		}

		uint64_t pt90prod[2];
		mulmod128(pt90.X, pt90.Z, pt90prod, mdata);

		for (i = 0; i < numsteps; i++)
		{
			b = steps[i];
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
			uint64_t tt1[2];
			uint64_t tt2[2];
			submod128(pt90.X, Pb[map[b]].X, tt1, work->n);
			addmod128(pt90.Z, Pb[map[b]].Z, tt2, work->n);
			uint64_t tt3[2];
			mulmod128(tt1, tt2, tt3, mdata);
			addmod128(tt3, Pbprod[map[b]], tt1, work->n);
			submod128(tt1, pt90prod, tt2, work->n);

			uint64_t tmp[2];
			mulmod128(acc, tt2, tmp, mdata);
			if ((tmp[0] == 0) && (tmp[1] == 0))
				break;
			copy128(tmp, acc);
		}
	}
	else if (stg1_max == 70)
	{
		// first process these b's with A=120
		static const int steps[15] = { 49,47,41,37,31,23,19,17,13,11,7,29,43,53,59 };
		// we currently have Pd=120

		uint64_t pdprod[2];
		mulmod128(Pd->X, Pd->Z, pdprod, mdata);

		for (i = 0; i < 15; i++)
		{
			b = steps[i];
			// accumulate the cross product  (zimmerman syntax).
			// page 342 in C&P
			uint64_t tt1[2];
			uint64_t tt2[2];
			submod128(Pd->X, Pb[map[b]].X, tt1, work->n);
			addmod128(Pd->Z, Pb[map[b]].Z, tt2, work->n);
			uint64_t tt3[2];
			mulmod128(tt1, tt2, tt3, mdata);
			addmod128(tt3, Pbprod[map[b]], tt1, work->n);
			submod128(tt1, pdprod, tt2, work->n);

			uint64_t tmp[2];
			mulmod128(acc, tt2, tmp, mdata);
			if ((tmp[0] == 0) && (tmp[1] == 0))
				break;
			copy128(tmp, acc);
		}
	}
	else if (stg1_max == 165)
	{
		// Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
		// Need Pa = 180, Pd = 120, Pad = 60
		// either of these should be fine
		submod128(pt90.X, pt90.Z, work->diff1, work->n);
		addmod128(pt90.X, pt90.Z, work->sum1, work->n);
		duplicate(mdata, work, work->sum1, work->diff1, Pa);

		//Pad = pt60;
		copy128(pt60.X, Pad.X);
		copy128(pt60.Z, Pad.Z);
		// have pa = 180, pd = 120, pad = 60
	}
	else if (stg1_max == 205)
	{
		// Currently we have Pa=150, Pd=120, Pad=30,  and pt60=60, pt90=90
		// need Pa = 210, Pd = 120, Pad = 90

		// [120]Q + [90]Q([30]Q) = [210]Q
		add(mdata, work, Pd, &pt90, &Pad, Pa);

		//Pad = pt90;
		copy128(pt90.X, Pad.X);
		copy128(pt90.Z, Pad.Z);
	}

	//initialize Paprod
	uint64_t Paprod[2];
	mulmod128(Pa->X, Pa->Z, Paprod, mdata);

	if (stg1_max == 27)
	{
		barray = b1_27;
		numb = numb1_27;
	}
	else if (stg1_max == 47)
	{
		barray = b1_47;
		numb = numb1_47;
	}
	else if (stg1_max <= 70)
	{
		barray = b1_70;
		numb = numb1_70;
	}
	else if (stg1_max == 85)
	{
		barray = b1_85;
		numb = numb1_85;
	}
	else if (stg1_max == 125)
	{
		barray = b1_125;
		numb = numb1_125;
	}
	else if (stg1_max == 165)
	{
		barray = b1_165;
		numb = numb1_165;
	}
	else if (stg1_max == 205)
	{
		barray = b1_205;
		numb = numb1_205;
	}

	for (i = 0; i < numb; i++)
	{
		if (barray[i] == 0)
		{
			//giant step - use the addition formula for ECM
			tinyecm_pt point;
			copy128(Pa->X, point.X);
			copy128(Pa->Z, point.Z);

			//Pa + Pd
			add(mdata, work, Pa, Pd, &Pad, Pa);

			//Pad holds the previous Pa
			copy128(point.X, Pad.X);
			copy128(point.Z, Pad.Z);

			//and Paprod
			mulmod128(Pa->X, Pa->Z, Paprod, mdata);

			i++;
		}

		//we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
		//in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)

		b = barray[i];
		// accumulate the cross product  (zimmerman syntax).
		// page 342 in C&P
		uint64_t tt1[2];
		uint64_t tt2[2];
		submod128(Pa->X, Pb[map[b]].X, tt1, work->n);
		addmod128(Pa->Z, Pb[map[b]].Z, tt2, work->n);
		uint64_t tt3[2];
		mulmod128(tt1, tt2, tt3, mdata);
		addmod128(tt3, Pbprod[map[b]], tt1, work->n);
		submod128(tt1, Paprod, tt2, work->n);

		uint64_t tmp[2];
		mulmod128(acc, tt2, tmp, mdata);
		if ((tmp[0] == 0) && (tmp[1] == 0))
			break;
		copy128(tmp, acc);
	}

	copy128(acc, work->stg2acc);
	return;
}

int check_factor(uint64_t * Z, uint64_t * n, uint64_t * f)
{
	int status;
	mpz_t gmp_z;
	mpz_t gmp_n;
	mpz_t gmp_f;

	mpz_init(gmp_z);
	mpz_init(gmp_n);
	mpz_init(gmp_f);

	u128_to_mpz(Z, gmp_z);
	u128_to_mpz(n, gmp_n);

	mpz_gcd(gmp_f, gmp_z, gmp_n);

	status = 0;
	if (mpz_cmp_ui(gmp_f, 1) > 0)
	{
		if (mpz_cmp(gmp_f, gmp_n) == 0)
		{
			mpz_set_ui(gmp_f, 0);
			status = 0;
		}
		else
		{
			status = 1;
		}
	}
	
	mpz_to_u128(gmp_f, f);

	mpz_clear(gmp_f);
	mpz_clear(gmp_n);
	mpz_clear(gmp_z);

	return status;
}

void tinyecm_test(int sizeb, int num, int type)
{
	u128_t *inputs;
	int i;

	if (sizeb > 128)
	{
		printf("size must be <= 128 bits\n");
		return;
	}

	if (num > (1 << 30))
	{
		printf("please be reasonable\n");
		return;
	}

	inputs = (u128_t *)malloc(num * sizeof(u128_t));

	// build the requested input type, quantity and size
	for (i = 0; i < num; i++)
	{


	}

	// test them, building up timings in globals


	// report results

	return;
}


int tecm_dispatch(mpz_t n, mpz_t f, int targetBits, int arbitrary, uint64_t* ploc_lcg)
{
	int B1, curves;

	//uecm_stage2_pair(70, 180, 120, 120, 30, primes);
	//uecm_stage2_pair(180, 25 * 70, 150, 120, 30, primes);

	//uecm_stage2_pair(47, 150, 90, 120, 30, primes);
	//uecm_stage2_pair(150, 25 * 47, 150, 120, 30, primes);


	//uecm_stage2_pair(30, 150, 90, 120, 30, primes);
	//uecm_stage2_pair(150, 25 * 30, 150, 120, 30, primes);
	//exit(0);

	if (arbitrary)
	{
		// try fast attempts to find possible small factors.
		//{
		//	B1 = 47;
		//	curves = 1;
		//	tinyecm(n, f, B1, B1 * 25, curves, ploc_lcg, 0);
		//	if (mpz_get_ui(f) > 1)
		//		return 1;
		//}
		{
			B1 = 70;
			curves = 1;
			tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
			if (mpz_get_ui(f) > 1)
				return 1;
		}
		if (targetBits > 58)
		{
			B1 = 125;
			curves = 1;
			tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
			if (mpz_get_ui(f) > 1)
				return 1;
		}
	}

	//if (targetBits <= 40)
	//{
	//	B1 = 27;
	//	curves = 32;
	//	tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	//}
	//else if (targetBits <= 44)
	//{
	//	B1 = 47;
	//	curves = 32;
	//	tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	//}
	if (targetBits <= 48)
	{
		B1 = 70;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 52)
	{
		B1 = 85;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 58)
	{
		B1 = 125;
		curves = 32;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 62)
	{
		B1 = 165;
		curves = 42;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}
	else if (targetBits <= 64)
	{
		B1 = 205;
		curves = 42;
		tinyecm(n, f, B1, 25 * B1, curves, ploc_lcg, 0);
	}

	if (mpz_get_ui(f) > 1)
		return 1;
	else
		return 0;
}


static int tecm_get_bits(mpz_t n)
{
	return mpz_sizeinbase(n, 2);
}

// getfactor_tecm() returns 0 if unable to find a factor of n,
// Otherwise it returns 1 and a factor of n in argument f.
// 
// if the input is known to have no small factors, set is_arbitrary=0, 
// otherwise, set is_arbitrary=1 and a few curves targetting small factors
// will be run prior to the standard sequence of curves for the input size.
//  
// Prior to your first call of getfactor_tecm(), set *pran = 0  (or set it to
// some other arbitrary value); after that, don't change *pran.
// FYI: *pran is used within this file by a random number generator, and it
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to
// change *pran, since that would restart the sequence.
int getfactor_tecm(mpz_t n, mpz_t f, int is_arbitrary, uint64_t* pran)
{
	if (mpz_even_p(n))
	{
		mpz_set_ui(f, 2);
		return 1;
	}

	int bits = tecm_get_bits(n);
	return tecm_dispatch(n, f, bits, is_arbitrary, pran);
}

