/*3:*/
#line 33 "strategy.w"

#include <stdio.h> 
#include <sys/types.h> 
#include <math.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <limits.h> 
#include <string.h> 
#include <time.h> 
#include <gmp.h> 
#include "asm/siever-config.h"
#include "if.h"
#include "primgen32.h"
#include "asm/32bit.h"
#include "strategy.h"
#include "asm/montgomery_mul.h"
#include "ecm.h"
#include "pm1.h"

#include "mpqs.h"
#include "mpqs3.h"

#include <immintrin.h>

extern char*input_line;
extern size_t input_line_alloc;

static mpz_t ecm_f[2],ecm_aux;
static u32_t cf_maxcomp[2];

static u32_t nfecm= 0,nsecm= 0,nfpm1= 0,nspm1= 0;
static i64_t mpqsaux_clock;

#define CF_STAT

// for tinyecm/microecm
#include "../gmp-ecm/microecm.h"
static uint64_t pran;
static mpz_t uecm_factors[3];
static int uecm_initialized = 0;

#ifdef CF_STAT
static double**stat_cost;
static double**stat_yield;
static double cost,yield;
static u64_t cf_n= 0,cf_necm= 0,cf_naux= 0,cf_nauxmpqs= 0,cf_nauxmpqs3= 0;
static u64_t cf_nauxmpqstoobig= 0,cf_nauxecm= 0;
#endif
#line 74 "strategy.w"

#ifdef CF_STAT_EXACT
static u32_t**stat_cand;
static u32_t**stat_success;
static u32_t*stat_mpqsaux,*stat_aux;
#endif
#line 80 "strategy.w"


/*:3*//*4:*/
#line 83 "strategy.w"

void print_strategy(strat_t s)
{
	u32_t i, j, k, l, t;

	printf("Strategy:\n");
	for (i = 0; i <= s.mc[0]; i++)
		for (j = 0; j <= s.mc[1]; j++)
			if (t = s.stindex[i][j]) {
				printf("%u %u:", i, j);
				for (k = 0; l = s.stlist[t][k]; k++) {
					if (l & 1)printf(" Y"); else printf(" X");
					if ((l >> 2) == 1)printf("M");
					else {
						if ((l >> 2) & 1)printf("P%u:%u", (l >> 3) & 0x1fff, 4 * (l >> 16));
						else printf("E%u:%u", (l >> 3) & 0x7ff, l >> 14);
					}
					if (l & 2)printf("|");
				}
				printf("\n");
			}
	printf("\n");
	printf("bit:\n");
	for (i = 0; i <= s.mc[0]; i++) {
		printf("%3u: ", i);
		for (j = 0; j <= s.mc[1]; j++) {
			if (s.bit[i][j])printf("X"); else printf(".");
		}
		printf("\n");
	}
	printf("\n");
}


/*:4*//*5:*/
#line 118 "strategy.w"

static u32_t get_pm1_type(u32_t B1, u32_t B2)
{
	if (B1 < 2)return 0;
	if ((B1 >> 13))return 0;
	if (B2 < 2 * B1)return 0;
	if ((B2 >> 18))return 0;
	return 1 + 2 * B1 + ((B2 >> 2) << 14);
}


/*:5*//*6:*/
#line 130 "strategy.w"

static void get_pm1_param(u32_t* b1ptr, u32_t* b2ptr, u32_t type)
{
	*b1ptr = (type >> 3) & 0x1fff;
	*b2ptr = 4 * (type >> 16);
}


/*:6*//*7:*/
#line 139 "strategy.w"

static u32_t get_ecm_type(u32_t B1, u32_t B2)
{
	if (B1 < 2)return 0;
	if ((B1 >> 12))return 0;
	if (B2 < 2 * B1)return 0;
	if ((B2 >> 20))return 0;
	return 2 * ((B1 >> 1) << 1) + ((B2 >> 2) << 12);
}


/*:7*//*8:*/
#line 151 "strategy.w"

static void get_ecm_param(u32_t* b1ptr, u32_t* b2ptr, u32_t type)
{
	*b1ptr = (type >> 2) & 0xffe;
	*b2ptr = (type >> 14) << 2;
}


/*:8*//*9:*/
#line 160 "strategy.w"

u32_t get_fm_type(char** lptr)
{
	char* line, * tail;
	u32_t t, B1, B2;

	line = *lptr;
	if (*line == 'M') {
		t = 1; line++;
	}
	else if (*line == 'E') {
		line++;
		B1 = (u32_t)strtoul(line, &tail, 10);
		if ((tail == NULL) || (line == tail))
			complain("st file contains corrupt line:\n%s", line - 1);
		line = tail + 1;
		B2 = (u32_t)strtoul(line, &tail, 10);
		if ((tail == NULL) || (line == tail))
			complain("st file contains corrupt line:\n%s", line - 1);
		line = tail;
		t = get_ecm_type(B1, B2);
		if (!t)complain("bad ecm parameter in line:\n%s\n", *lptr);
	}
	else if (*line == 'P') {
		line++;
		B1 = (u32_t)strtoul(line, &tail, 10);
		if ((tail == NULL) || (line == tail))
			complain("st file contains corrupt line:\n%s", line - 1);
		line = tail + 1;
		B2 = (u32_t)strtoul(line, &tail, 10);
		if ((tail == NULL) || (line == tail))
			complain("st file contains corrupt line:\n%s", line - 1);
		line = tail;
		t = get_pm1_type(B1, B2);
		if (!t)complain("bad pm1 parameter in line:\n%s\n", *lptr);
	}
	else complain("unknown factorisation type\n");
	*lptr = line;
	return t;
}


/*:9*//*10:*/
#line 200 "strategy.w"

void read_strategy(strat_t* s, u16_t* maxcomp, char* basename, u16_t* maxpr)
{
	char* ifn;
	FILE* ifile;
	char* line, * tail;
	u32_t n0, n1, len;
	u32_t i, j, k, t, tt, b;
	size_t alloc = 0;
	u32_t* s_tmp;

	cf_maxcomp[0] = (u32_t)maxcomp[0];
	cf_maxcomp[1] = (u32_t)maxcomp[1];
	mpz_init(ecm_f[0]);
	mpz_init(ecm_f[1]);
	mpz_init(ecm_aux);
	init_montgomery_multiplication();
	memcpy(s->mc, maxcomp, 2 * sizeof(u16_t));
	s->stindex = (u32_t**)xmalloc((1 + maxcomp[0]) * sizeof(u32_t*));
	for (i = 0; i <= maxcomp[0]; i++)
		s->stindex[i] = (u32_t*)xcalloc((size_t)(1 + maxcomp[1]), sizeof(u32_t));
	s->stindex[0][0] = 1;
	adjust_bufsize((void**)&(s->stlist), &alloc, 1, 16, sizeof(u32_t*));
	s->stlist[0] = (u32_t*)xcalloc(1, sizeof(u32_t));
	s->stlist[0][0] = 0;
	len = 1;
#ifdef CF_STAT
	stat_cost = (double**)xmalloc((1 + maxcomp[0]) * sizeof(double*));
	for (i = 0; i <= maxcomp[0]; i++)
		stat_cost[i] = (double*)xcalloc((size_t)(1 + maxcomp[1]), sizeof(double));
	stat_yield = (double**)xmalloc((1 + maxcomp[0]) * sizeof(double*));
	for (i = 0; i <= maxcomp[0]; i++)
		stat_yield[i] = (double*)xcalloc((size_t)(1 + maxcomp[1]), sizeof(double));
	cost = 0.; yield = 0.;
#endif
#line 235 "strategy.w"
#ifdef CF_STAT_EXACT
	stat_cand = (u32_t**)xmalloc((1 + maxcomp[0]) * sizeof(u32_t*));
	for (i = 0; i <= maxcomp[0]; i++)
		stat_cand[i] = (u32_t*)xcalloc((size_t)(1 + maxcomp[1]), sizeof(u32_t));
	stat_success = (u32_t**)xmalloc((1 + maxcomp[0]) * sizeof(u32_t*));
	for (i = 0; i <= maxcomp[0]; i++)
		stat_success[i] = (u32_t*)xcalloc((size_t)(1 + maxcomp[1]), sizeof(u32_t));
	i = (maxcomp[0] < maxcomp[1] ? maxcomp[1] : maxcomp[0]);
	stat_aux = (u32_t*)xcalloc((size_t)(1 + i), sizeof(u32_t));
	stat_mpqsaux = (u32_t*)xcalloc((size_t)(1 + i), sizeof(u32_t));
#endif
#line 246 "strategy.w"

	asprintf(&ifn, "%s.st", basename);
	if ((ifile = fopen(ifn, "r")) != 0) {
		while (1) {
			if (skip_blanks_comments(&input_line, &input_line_alloc, ifile) == 0)break;
			if (input_line == NULL)break;
			line = input_line;
			n0 = (u32_t)strtol(line, &tail, 10);
			if ((tail == NULL) || (line == tail))
				complain("file %s contains corrupt line:\n%s", ifn, input_line);
			line = tail;
			n1 = (u32_t)strtol(line, &tail, 10);
			if ((tail == NULL) || (line == tail))
				complain("file %s contains corrupt line:\n%s", ifn, input_line);
			line = tail;
			if (n0 > maxcomp[0])continue;
			if (n1 > maxcomp[1])continue;
#ifdef CF_STAT
			line++;
			stat_cost[n0][n1] = strtod(line, &tail);
			if ((tail == NULL) || (line == tail))
				complain("file %s contains corrupt line:\n%s", ifn, input_line);
			line = tail + 1;
			stat_yield[n0][n1] = strtod(line, &tail);
			if ((tail == NULL) || (line == tail))
				complain("file %s contains corrupt line:\n%s", ifn, input_line);
			line = tail;
#endif
#line 274 "strategy.w"
			while (*line) {
				if (*line == 'X')break;
				if (*line == 'Y')break;
				line++;
			}
			s->stindex[n0][n1] = len;
			for (tail = line, i = 1; *tail; tail++)if (*tail == ',')i++;
			adjust_bufsize((void**)&(s->stlist), &alloc, (size_t)(len + 1),
				16, sizeof(u32_t*));
			s_tmp = (u32_t*)xcalloc((size_t)(i + 1), sizeof(u32_t));
			for (j = 0; j < i; j++) {
				t = 0;
				if (*line == 'Y')t++; else if (*line != 'X')complain("no XY\n");
				line++;
				tt = get_fm_type(&line);
				t += 4 * tt;
				s_tmp[j] = t;
				if (j + 1 < i)if (*line != ',')complain("no ,  %u %u\n", i, j);
				line++;
			}
			if (*line)complain("st-file corrupt\n");
			s_tmp[j++] = 0;

			for (i = j - 1; i; i--)if ((s_tmp[i - 1] & 1) == 0) { s_tmp[i - 1] ^= 2; break; }
			for (i = j - 1; i; i--)if ((s_tmp[i - 1] & 1) == 1) { s_tmp[i - 1] ^= 2; break; }

			for (i = 0; i < len; i++) {
				for (k = 0; k < j; k++)if (s_tmp[k] != s->stlist[i][k])break;
				if (k == j)break;
			}
			if (i < len)s->stindex[n0][n1] = i;
			else {
				s->stlist[len] = (u32_t*)xmalloc(j * sizeof(u32_t));
				memcpy(s->stlist[len], s_tmp, j * sizeof(u32_t));
				len++;
			}
			free(s_tmp);
		}
		fclose(ifile);
	}
	else {

		if ((3 * maxpr[0] + 3 < maxcomp[0]) || (3 * maxpr[1] + 3 < maxcomp[1])) {
			logbook(0, "Warning: >2LP but no .st file\n");
		}

		adjust_bufsize((void**)&(s->stlist), &alloc, (size_t)(len + 2),
			16, sizeof(u32_t*));
		s->stlist[len] = (u32_t*)xmalloc(3 * sizeof(u32_t));
		s->stlist[len][0] = 4 + 2; s->stlist[len][1] = 4 + 2 + 1; s->stlist[len][2] = 0;
		len++;
		s->stlist[len] = (u32_t*)xmalloc(3 * sizeof(u32_t));
		s->stlist[len][0] = 4 + 2 + 1; s->stlist[len][1] = 4 + 2; s->stlist[len][2] = 0;
		len++;
		for (n0 = 0; n0 <= maxcomp[0]; n0++) {
			for (n1 = 0; n1 <= maxcomp[1]; n1++) {
				if (n0 < n1)s->stindex[n0][n1] = len - 1;
				else s->stindex[n0][n1] = len - 2;


			}
		}
	}


	s->bit = (unsigned char**)xmalloc((maxcomp[0] + 1) * sizeof(unsigned char*));
	for (i = 0; i <= maxcomp[0]; i++)
		s->bit[i] = (unsigned char*)xcalloc((size_t)(maxcomp[1] + 1), sizeof(unsigned char));
	{
		u32_t i0, i1, j0, j1;

		for (i = 0; i <= maxcomp[0]; i++) {
			if (i > 3)i0 = i - 3; else i0 = 0;
			if (i == 0)i1 = maxpr[0]; else i1 = i;
			i1 += 1; if (i1 >= maxcomp[0])i1 = maxcomp[0];
			for (j = 0; j <= maxcomp[1]; j++) {
				if (j > 3)j0 = j - 2; else j0 = 0;
				if (j == 0)j1 = maxpr[1]; else j1 = j;
				j1 += 0; if (j1 >= maxcomp[1])j1 = maxcomp[1];
				if (s->stindex[i][j])
					for (k = i0; k <= i1; k++)
						for (t = j0; t <= j1; t++)
							s->bit[k][t] = 1;
			}
		}
	}
}


/*:10*//*11:*/
#line 363 "strategy.w"

int cofactorisation(strat_t* st, mpz_t** large_primes, mpz_t* large_factors,
	u16_t* max_primebits, u32_t* nlp, mpz_t* FBb_sq,
	mpz_t* FBb_cu)
{
	u32_t s, nb[2];
	u32_t m, * fm, t, j, done[2], B1, B2, pm1done[2];
	static ecm_t e[2];
	static int cf_ecm_init = 0;
	clock_t cl;

	if (!uecm_initialized)
	{
		// for tinyecm
		mpz_init(uecm_factors[0]);
		mpz_init(uecm_factors[1]);
		mpz_init(uecm_factors[2]);
		pran = 42;
		uecm_initialized = 1;
	}

	for (s = 0; s < 2; s++) {
		if (mpz_sgn(large_factors[s]) > 0) {
			if (mpz_cmp_ui(large_factors[s], 1) == 0)
				nlp[s] = 0;
			else {
				nlp[s] = 1;
				mpz_set(large_primes[s][0], large_factors[s]);
			}
			nb[s] = 0;
		}
		else {
			mpz_neg(large_factors[s], large_factors[s]);
			nb[s] = (u32_t)(mpz_sizeinbase(large_factors[s], 2));
			nlp[s] = 2;
		}
	}
	if ((nlp[0] < 2) && (nlp[1] < 2)) {
#ifdef CF_STAT
		yield += 1.;
#endif
#line 393 "strategy.w"
#ifdef CF_STAT_EXACT
		stat_cand[0][0]++;
		stat_success[0][0]++;
#endif
#line 397 "strategy.w"
		return 0;
	}


#ifdef CF_STAT
	yield += stat_yield[nb[0]][nb[1]];
	cost += stat_cost[nb[0]][nb[1]];
#endif
#line 404 "strategy.w"
#ifdef CF_STAT_EXACT
	stat_cand[nb[0]][nb[1]]++;
#endif
#line 407 "strategy.w"

	if (cf_ecm_init == 0) {
		ecm_curve_init(e[0]);
		ecm_curve_init(e[1]);
		cf_ecm_init = 1;
	}
	m = st->stindex[nb[0]][nb[1]];
	if (!m)return 1;
	fm = st->stlist[m];
#ifdef CF_STAT
	cf_n++;
#endif
#line 419 "strategy.w"
	for (s = 0; s < 2; s++)
		if (nlp[s] == 2) { nlp[s] = 0; done[s] = 0; pm1done[s] = 0; }
		else done[s] = 2;
	for (j = 0; t = fm[j]; j++) {
		i32_t nf, i;
		size_t sf[2];
		mpz_t* fac;

		s = t & 1;
		if (done[s] > 1)continue;
		if ((t >> 2) == 1) {

			// unless some kind of strategy file is read in, this
			// is the only cofactorization case that is used.
			// so we can drop in the tiny/micro ecm stuff here.


#ifdef GGNFS_MPQS

			if (mpz_sizeinbase(large_factors[s], 2) > 96)
				nf = mpqs3_factor(large_factors[s], max_primebits[s], &fac);
			else
				nf = mpqs_factor(large_factors[s], max_primebits[s], &fac);

#else
			nf = 0;
			fac = uecm_factors;

			if (mpz_sizeinbase(large_factors[s], 2) <= 64) {
				uint64_t n64 = mpz_get_ui(large_factors[s]);
				uint64_t f = getfactor_uecm(n64, 0, &pran);
				if (f > 1)
				{
					mpz_set_ui(fac[0], f);
					mpz_tdiv_q_ui(fac[1], large_factors[s], f);
					nf = 2;

					if (mpz_sizeinbase(fac[0], 2) > max_primebits[s]) 
					{
						nf = 0;
					}
					if (mpz_sizeinbase(fac[1], 2) > max_primebits[s]) 
					{
						nf = 0;
					}
					if (mpz_probab_prime_p(fac[0], 1) == 0)
					{
						nf = 0;
					}
					if (mpz_probab_prime_p(fac[1], 1) == 0)
					{
						nf = 0;
					}
				}
				else
				{
					nf = 0;
				}
			}
			else
			{
				if (getfactor_tecm(large_factors[s], fac[0],
					mpz_sizeinbase(large_factors[s], 2) / 3 - 2, &pran) > 0)
				{
					if (mpz_sizeinbase(fac[0], 2) <= max_primebits[s])
					{
						mpz_tdiv_q(fac[1], large_factors[s], fac[0]);

						// if the remaining residue is obviously too big, we're done.
						if (mpz_sizeinbase(fac[1], 2) > ((max_primebits[s] * 2) + 1))
						{
							nf = 0;
							goto done;
						}

						// check if the residue is prime.  could again use
						// a cheaper method.
						if (mpz_probab_prime_p(fac[1], 1) > 0)
						{
							if (mpz_sizeinbase(fac[1], 2) <= max_primebits[s])
							{
								// we just completed a DLP factorization involving
								// 2 primes whos product was > 64 bits.
								nf = 2;
								goto done;
							}
							nf = 0;
							goto done;
						}

						// ok, so we have extracted one suitable factor, and the 
						// cofactor is not prime and a suitable size.  Do more work to 
						// split the cofactor.
						// todo: target this better based on expected factor size.
						uint64_t q64;
						uint64_t f64;
						if (mpz_sizeinbase(fac[1], 2) <= 64)
						{
							q64 = mpz_get_ui(fac[1]);
							f64 = getfactor_uecm(q64, 0, &pran);
							mpz_set_ui(fac[2], f64);
						}
						else
						{
							getfactor_tecm(fac[1], fac[2], 32, &pran);
						}
						f64 = mpz_get_ui(fac[2]);

						if (f64 > 1)
						{
							mpz_tdiv_q_ui(fac[1], fac[1], f64);
							nf = 3;

							if (mpz_sizeinbase(fac[1], 2) > max_primebits[s]) {
								nf = 0;
							}
							if (mpz_sizeinbase(fac[2], 2) > max_primebits[s]) {
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[0], 1) == 0)
							{
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[1], 1) == 0)
							{
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[2], 1) == 0)
							{
								nf = 0;
							}
						}
						else
						{
							// what are the odds that the largefactor is p1 * p2^2?
							// and we've found p1 but ecm fails on the p2^2?
							// tinyecm checks for squares if it finds gcd==N,
							// but uecm doesn't.
							nf = 0;
						}
					}
					else
					{
						// check if the factor is prime.  could again use
						// a cheaper method.
						if (mpz_probab_prime_p(fac[0], 1) > 0)
						{
							// if the factor is obviously too big, give up.  This isn't a
							// failure since we haven't expended much effort yet.
							nf = 0;
						}
						else
						{
							// tecm found a composite first factor.
							// if it is obviously too big, we're done.
							if (mpz_sizeinbase(fac[0], 2) > ((max_primebits[s] * 2) + 1))
							{
								nf = 0;
								goto done;
							}

							// isolate the 2nd smaller factor, and check its size.
							mpz_tdiv_q(fac[1], large_factors[s], fac[0]);

							if (mpz_sizeinbase(fac[1], 2) > (max_primebits[s]))
							{
								nf = 0;
								goto done;
							}

							// todo: target this better based on expected factor size.
							uint64_t q64;
							uint64_t f64;
							if (mpz_sizeinbase(fac[0], 2) <= 64)
							{
								q64 = mpz_get_ui(fac[0]);
								f64 = getfactor_uecm(q64, 0, &pran);
								mpz_set_ui(fac[2], f64);
							}
							else
							{
								getfactor_tecm(fac[0], fac[2], 32, &pran);
							}
							f64 = mpz_get_ui(fac[2]);

							if (f64 > 1)
							{
								mpz_tdiv_q_ui(fac[0], fac[0], f64);
								nf = 3;

								if (mpz_sizeinbase(fac[0], 2) > max_primebits[s]) {
									nf = 0;
								}
								if (mpz_sizeinbase(fac[2], 2) > max_primebits[s]) {
									nf = 0;
								}
								if (mpz_probab_prime_p(fac[0], 1) == 0)
								{
									nf = 0;
								}
								if (mpz_probab_prime_p(fac[1], 1) == 0)
								{
									nf = 0;
								}
								if (mpz_probab_prime_p(fac[2], 1) == 0)
								{
									nf = 0;
								}
								
							}
							else
							{
								nf = 0;
							}
						}
					}
				}
				else
				{
					// if ecm can't find a factor, give up.  
					// unless this is a DLP with lpbr/a > 32... i.e., if the
					// large factor size is greater than 64 bits but less than
					// lpbr/a * 2.  In that case run mpqs... or tecm with
					// greater effort.

#if 0
					if (mpz_sizeinbase(large_factors[s1], 2) <= (max_primebits[s1] * 2))
					{
						if (getfactor_tecm(large_factors[s1], factor1, 33, &pran) > 0)
						{
							if (mpz_sizeinbase(factor1, 2) <= max_primebits[s1])
							{
								mpz_tdiv_q(factor2, large_factors[s1], factor1);

								// check if the residue is prime.  could again use
								// a cheaper method.
								if (mpz_probab_prime_p(factor2, 1) > 0)
								{
									if (mpz_sizeinbase(factor2, 2) <= max_primebits[s1])
									{
										// we just completed a DLP factorization involving
										// 2 primes whos product was > 64 bits.
										mpz_set(large_primes[s1][0], factor1);
										mpz_set(large_primes[s1][1], factor2);
										nlp[s1] = 2;
									}
									else
										break;
								}
								else
									break;
							}
							else
								break;
						}
						else
							break;
					}
					else
						break;
#else

					if (mpz_sizeinbase(large_factors[s], 2) <= (max_primebits[s] * 2))
					{
						nf = mpqs_factor(large_factors[s], max_primebits[s], &fac);
					}
					else
					{
						nf = 0;
					}
#endif
				}

			}

		done:

			_mm256_zeroupper();

		


#endif

			if (nf < 0)return-2;
			if (!nf)return 1;
			for (i = 0; i < nf; i++)
				mpz_set(large_primes[s][nlp[s] + i], fac[i]);
			nlp[s] += nf;
			done[s] = 2;


		}
		else {
			if ((t >> 2) & 1) {
				printf("trying pm1\n");
				get_pm1_param(&B1, &B2, t);
				nf = pm1_factor(large_factors[s], B1, B2, &fac);
				pm1done[s] = 1;
				nfpm1++; if (nf > 0)nspm1++;
			}
			else {
				printf("trying ecm\n");
				get_ecm_param(&B1, &B2, t);
				if (!done[s]) {
					if (ecm_curve_set(e[s], large_factors[s], B1, B2))return-3;
					done[s] = 1;
				}
				else {
					ecm_set_params(e[s], B1, B2);
				}
				nf = ecm(e[s], &fac);
				nfecm++; if (nf > 0)nsecm++;
			}
			if (nf < 0)return-3;
			if (nf) {
				/*12:*/
#line 480 "strategy.w"
				printf("trying big ecm\n");
				{
					int es, need_test[2], order[2], o;

#ifdef CF_STAT
					cf_necm++;
#endif
#line 487 "strategy.w"
					mpz_set(ecm_f[0], fac[0]);
					sf[0] = mpz_sizeinbase(ecm_f[0], 2);
					if (sf[0] <= 1)return-4;
					if (sf[0] < nb[s]) {
						mpz_fdiv_qr(ecm_f[1], ecm_aux, large_factors[s], ecm_f[0]);
						if (mpz_sgn(ecm_aux))Schlendrian("ecm found non-divisor\n");
						sf[1] = mpz_sizeinbase(ecm_f[1], 2);

						for (es = 0; es < 2; es++) {
							need_test[es] = 0;
							if (sf[es] <= max_primebits[s])continue;
							if (mpz_cmp(ecm_f[es], FBb_sq[s]) < 0)return 1;
							if (sf[es] <= 2 * max_primebits[s]) { need_test[es] = 1; continue; }
							if (mpz_cmp(ecm_f[es], FBb_cu[s]) < 0)return 1;
							need_test[es] = 1;
						}

						for (es = 0; es < 2; es++)
							if (need_test[es])

								if (psp(ecm_f[es]) == 1)return 1;

						order[0] = 0; order[1] = 1;

						if ((need_test[0] & need_test[1]) && (sf[0] > sf[1])) {
							order[0] = 1; order[1] = 0;
						}
						for (o = 0; o < 2; o++) {
							es = order[o];
							if (!need_test[es]) {
								mpz_set(large_primes[s][nlp[s]], ecm_f[es]); nlp[s]++;
							}
							else {
#ifdef CF_STAT_EXACT
								stat_aux[mpz_sizeinbase(ecm_f[es], 2)]++;
#endif
#line 522 "strategy.w"
#ifdef CF_STAT
								cf_naux++;
#endif
#line 525 "strategy.w"
								/*13:*/
#line 534 "strategy.w"

								{
									u32_t ne;
									size_t sz;

									cl = clock();
									ne = 0;
									while (1) {
										sz = mpz_sizeinbase(ecm_f[es], 2);
										if (sz <= 96) {
#ifdef CF_STAT_EXACT
											stat_mpqsaux[sz]++;
#endif
#line 547 "strategy.w"
#ifdef CF_STAT
											cf_nauxmpqs++;
#endif
#line 550 "strategy.w"
											nf = mpqs_factor(ecm_f[es], max_primebits[s], &fac);
										}
										else {
											if (5 * ne > (sz - 90)) {
#ifdef CF_STAT
												cf_nauxmpqs3++;
#endif
#line 556 "strategy.w"
#ifdef CF_STAT_EXACT
												stat_mpqsaux[sz]++;
#endif
#line 559 "strategy.w"
												if (sz > 128) {
#ifdef CF_STAT
													cf_nauxmpqstoobig++;
#endif
#line 563 "strategy.w"
													mpqsaux_clock += clock() - cl; return-1;
												}
												nf = mpqs3_factor(ecm_f[es], max_primebits[s], &fac);
											}
											else {
												B1 = 500; B2 = 36000;
												if (!done[s]) {
													if (ecm_curve_set(e[s], ecm_f[es], B1, B2)) {
														mpqsaux_clock += clock() - cl;
														return-3;
													}
													done[s] = 1;
												}
												else {
													if (!ne) {
														if (ecm_reset_n(e[s], ecm_f[es])) {
															mpqsaux_clock += clock() - cl;
															return-3;
														}
													}
													ecm_set_params(e[s], B1, B2);
												}
#ifdef CF_STAT
												cf_nauxecm++;
#endif
#line 586 "strategy.w"
												nf = ecm(e[s], &fac);
												nfecm++; if (nf > 0)nsecm++;
												ne++;
												if (nf < 0) { mpqsaux_clock += clock() - cl; return-3; }
												if (nf == 0)continue;

												if (mpz_sizeinbase(fac[0], 2) <= max_primebits[s]) {
													mpz_set(large_primes[s][nlp[s]], fac[0]); nlp[s]++;
													mpz_fdiv_qr(ecm_f[es], ecm_aux, ecm_f[es], fac[0]);
													if (mpz_sgn(ecm_aux)) {
														gmp_printf("%Zd %Zd\n%Zd\n", ecm_aux, fac[0], large_factors[s]);
														Schlendrian("ecm found non-divisor\n");
													}
													if (mpz_sizeinbase(ecm_f[es], 2) <= max_primebits[s]) {
														mpz_set(large_primes[s][nlp[s]], ecm_f[es]); nlp[s]++;
														break;

													}
													else if (psp(ecm_f[es]) == 1) {
														mpqsaux_clock += clock() - cl; return 1;
													}
													if (ecm_reset_n(e[s], ecm_f[es])) {
														mpqsaux_clock += clock() - cl;
														return-3;
													}
													continue;
												}

												if (psp(fac[0]) == 1) { mpqsaux_clock += clock() - cl; return 1; }
												continue;
											}
										}
										if (nf > 0) {
											for (i = 0; i < nf; i++)
												mpz_set(large_primes[s][nlp[s] + i], fac[i]);
											nlp[s] += nf;
											break;
										}
										mpqsaux_clock += clock() - cl;
										if (nf < 0)return-4;
										if (!nf)return 1;
									}
									mpqsaux_clock += clock() - cl;
								}

								/*:13*/
#line 525 "strategy.w"
								;
							}
						}
						done[s] = 2;
					}

				}

				/*:12*/
#line 458 "strategy.w"

			}
		}
		if (t & 2) {
			if (done[s] < 2)return 1;
		}
	}



	if ((done[0] != 2) || (done[1] != 2))return-1;

#ifdef CF_STAT_EXACT
	stat_success[nb[0]][nb[1]]++;
#endif
#line 470 "strategy.w"
	return 0;
}

/*:11*//*14:*/
#line 631 "strategy.w"

void print_strategy_stat()
{
#ifdef CF_STAT
logbook(0,"Expected yield/cost: %.4g  %.4g\n",yield,cost);
logbook(0,"p-1: %u tests, %u successes  ecm: %u tests, %u successes\n",
nfpm1,nspm1,nfecm,nsecm);
mpqsaux_clock= rint((1000.0*mpqsaux_clock)/CLOCKS_PER_SEC);
logbook(0,"MPQS-AUX %d\n",(int)mpqsaux_clock);
logbook(0,"COF: %llu tests, %llu ecm, %llu aux:\n",cf_n,cf_necm,cf_naux);
logbook(0,"       %llu mpqs, %llu mpqs3, %llu ecm, %llu too big\n",
cf_nauxmpqs,cf_nauxmpqs3,cf_nauxecm,cf_nauxmpqstoobig);
#endif
#line 644 "strategy.w"

#ifdef CF_STAT_EXACT
{
u32_t i,j;

logbook(0,"\nExact stat:\n");
for(i= 0;i<=cf_maxcomp[0];i++)
for(j= 0;j<=cf_maxcomp[1];j++)
if(stat_cand[i][j]){
logbook(0,"%u %u: %u -> %u (%g)\n",i,j,stat_cand[i][j],
stat_success[i][j],((double)stat_cand[i][j])*stat_yield[i][j]);
}
logbook(0,"\nmpqs-aux:\n");
j= (cf_maxcomp[0]<cf_maxcomp[1]?cf_maxcomp[1]:cf_maxcomp[0]);
for(i= 0;i<=j;i++)
if(stat_aux[i])
logbook(0,"%u: %u (%u)\n",i,stat_aux[i],stat_mpqsaux[i]);
}
#endif
#line 663 "strategy.w"
}/*:14*/
