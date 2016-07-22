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

#define LIMB_BLKSZ 10	
#define MAX_DIGITS 100
#define MP_RADIX 4294967296.0
#define LN2		0.69314718055994530942

//types of numbers
#define PRIME 0
#define PRP 1
#define COMPOSITE 2
#define UNKNOWN 3

#define fp_clamp(a)   { while ((a)->size && (a)->val[(a)->size-1] == 0) --((a)->size);}

uint32 mp_modadd_1(uint32 a, uint32 b, uint32 p);
uint32 mp_modsub_1(uint32 a, uint32 b, uint32 p);

//basic arbitrary precision arith routines
//contents of arith1.c
/********************* single precision arith **********************/
void spAdd(fp_digit u, fp_digit v, fp_digit *sum, fp_digit *carry);
void spAdd3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sum, fp_digit *carry);
void spSub3(fp_digit u, fp_digit v, fp_digit w, fp_digit *sub, fp_digit *borrow);
void spSub(fp_digit u, fp_digit v, fp_digit *sub, fp_digit *borrow);
void spMultiply(fp_digit u, fp_digit v, fp_digit *product, fp_digit *carry);
void spMulAdd(fp_digit u, fp_digit v, fp_digit w, fp_digit t, fp_digit *lower, fp_digit *carry);
void spMulMod(fp_digit u, fp_digit v, fp_digit m, fp_digit *w);
void spModExp(fp_digit a, fp_digit b, fp_digit m, fp_digit *u);
fp_digit spDivide(fp_digit *q, fp_digit *r, fp_digit u[2], fp_digit v);
fp_digit spBits(fp_digit n);
int bits64(uint64 n);
uint32 modinv_1(uint32 a, uint32 p);
uint32 modinv_1b(uint32 a, uint32 p);
uint32 modinv_1c(uint32 a, uint32 p);
int shortCompare(fp_digit p[2], fp_digit t[2]);
int shortSubtract(fp_digit p[2], fp_digit t[2], fp_digit w[2]);
uint64 spModExp_asm(uint64 b, uint64 e, uint64 m);
uint64 spPRP2(uint64 p);

/********************* arbitrary precision arith **********************/
//add
void zAdd(z *u, z *v, z *w);
void zAddb(z *u, z *v, z *y);
void zShortAdd(z *u, fp_digit v, z *w);
//sub
int zSub(z *u, z *v, z *w);
int zSubb(z *u, z *v, z *y);
void zShortSub(z *u, fp_digit v, z *w);

//mul
void zMul(z *u, z *v, z *w);
void zSqr(z *x, z *w);
void zShortMul(z *u, fp_digit v, z *w);
void zModMul(z *u, z *v, z *n, z *w);

//div
void zDiv(z *u, z *v, z *q, z *r);
fp_digit zShortDiv(z *u, fp_digit v, z *q);
fp_digit zShortMod(z *u, fp_digit v);
uint32 zShortEDiv32(z32 *u, uint32 v, uint32 inv);
uint32 zShortDiv32(z32 *u, uint32 v, z32 *q);
uint32 zShortMod32(z32 *u, uint32 v);
#ifdef _WIN64
	//uint32 mod_64(uint64 lo, uint64 hi, uint64 div);
#endif

//non-basic arbitrary precision arith routines
//contents of arith2.c
void zShiftLeft(z *a, z *b, int x);
void zShiftLeft_x(z *a, z *b, int x);
void zShiftLeft_1(z *a, z *b);
void zShiftLeft32(z32 *a, z32 *b, int x);
void zShiftRight(z *a, z *b, int x);
void zShiftRight_x(z *a, z *b, int x);
void zShiftRight_1(z *a, z *b);
void zShiftRight32(z32 *a, z32 *b, int x);
void zShiftRight32_x(z32 *a, z32 *b, int x);
int zNroot(z *u, z *w, int n);
int zFactorial(uint32 n, z *w);
int zPrimorial(uint32 n, z *w);
int zExp(uint32 k, z *u, z *w);
void zNeg(z *u);
void zModExp(z *a, z *b, z *m, z *u);
void zModExp_1(z *a, fp_digit b, fp_digit m, fp_digit *u);
int zBits(z *n);
void zShanksTonelli(z *n, fp_digit p, fp_digit *sq);
void ShanksTonelli_1(fp_digit a, fp_digit p, fp_digit *sq);
double zlog(mpz_t x);
void sim_mul_exp(uint16 *e, uint64 *g, mpz_t A, int k);
int getcolval(uint8 **em, int i, int k);

//applied arbitrary precision arith routines
//contents of arith3.c
fp_digit spGCD(fp_digit x, fp_digit y);
uint64 spBinGCD(uint64 x, uint64 y);
uint64 spBinGCD_odd(uint64 u, uint64 v);
uint64 gcd64(uint64 x, uint64 y);
void xGCD(z *a, z *b, z *x, z *y, z *g);
void dblGCD(double x, double y, double *w);
int dblFactorA(double *n, long p[], uint32 limit);
int isPrime(z *n);
int d_jacobi(double n, double p);
int zJacobi(z *n, z *p);
int rec_jacobi_1(uint32 n, uint32 p);
int jacobi_1(fp_digit n, fp_digit p);
int isSquare(z *n);
int llt(uint32 exp);
void lucas(uint32 n, z *L);
void fib(uint32 n, z *F);

//utilities for init/free/print/compare/convert 
//of arbitrary precision structures
//basically, the contents of arith0.c
void mp2gmp(z *src, mpz_t dest);
void gmp2mp(mpz_t src, z *dest);
void zInit(z *num);
void zInit32(z32 *num);
void zGrow(z *num, int newsz);
void zGrow32(z32 *num, int newsz);
void zCopy(z *src, z *dest);
void zCopy32(z32 *src, z32 *dest);
void zFree(z *num);
void zFree32(z32 *num);
void zClear(z *num);
void zClear32(z32 *num);
void str2hexz(char s[], z *u);
int isZero(z *n);
int isOne(z *n);
int isFive(z *n);
double z2dbl(z *a);
void sp642z(uint64 sp, z *mp);
void sp2z(fp_digit sp, z *mp);
uint64 z264(z *n);
int ndigits(z *n);
int ndigits_1(fp_digit n);
int gmp_base10(mpz_t x);
char *z2decstr(z *n, str_t *s);
char *z2hexstr(z *n, str_t *s);
int zCompare(z *u, z *v);
int zCompare32(z32 *u, z32 *v);
void zHex2Dec(z *u, z *v);
void zDec2Hex(z *u, z *v);
void swap(z *a, z *b);
double rint(double x);

// useful mpz utilities
int is_mpz_prp(mpz_t n);
uint64 mpz_get_64(mpz_t src);
void mpz_set_64(mpz_t dest, uint64 src);
void mpz_to_z32(mpz_t src, z32 *dest);
void z32_to_mpz(z32 *src, mpz_t dest);
char * mpz_conv2str(char **in, int base, mpz_t n);

// we need to convert between yafu bigints and msieve bigints occasionally
void mp_t2z(mp_t *src, z *dest);

char * mp_print(mp_t *a, uint32 base, FILE *f, char *scratch);
#define mp_sprintf(a, base, scratch) mp_print(a, base, NULL, scratch)

