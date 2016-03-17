/* Copyright 2011-2015 David Cleaver
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __MPZ_APRCL__
#define __MPZ_APRCL__

#ifndef HAVE_U64_T
#define HAVE_U64_T
typedef long long s64_t;
typedef unsigned long long u64_t;
#endif

#include "jacobi_sum.h"

/*
 * The PRP functions presented here are based on the paper:
 * Grantham, Jon. Frobenius Pseudoprimes. Math. Comp. 70 (2001), 873-891.
 */

/*************************************************************/
/*************************************************************/
/* These are the definitions for the probable prime routines */
/*************************************************************/
/*************************************************************/
#define PRP_ERROR -1
#define PRP_COMPOSITE 0
#define PRP_PRP 1
#define PRP_PRIME 2

/* ******************************************************************
 * mpz_prp: (also called a Fermat pseudoprime)
 * A "pseudoprime" to the base a is a composite number n such that,
 * (a,n)=1 and a^(n-1) = 1 mod n
 * ******************************************************************/
int mpz_prp(mpz_t n, mpz_t a);

/* *************************************************************************
 * mpz_euler_prp: (also called a Solovay-Strassen pseudoprime)
 * An "Euler pseudoprime" to the base a is an odd composite number n with,
 * (a,n)=1 such that a^((n-1)/2)=(a/n) mod n [(a/n) is the Jacobi symbol]
 * *************************************************************************/
int mpz_euler_prp(mpz_t n, mpz_t a);

/* *********************************************************************************************
 * mpz_sprp: (also called a Miller-Rabin pseudoprime)
 * A "strong pseudoprime" to the base a is an odd composite n = (2^r)*s+1 with s odd such that
 * either a^s == 1 mod n, or a^((2^t)*s) == -1 mod n, for some integer t, with 0 <= t < r.
 * *********************************************************************************************/
int mpz_sprp(mpz_t n, mpz_t a);

/* *************************************************************************
 * mpz_fibonacci_prp:
 * A "Fibonacci pseudoprime" with parameters (P,Q), P > 0, Q=+/-1, is a
 * composite n for which V_n == P mod n
 * [V is the Lucas V sequence with parameters P,Q]
 * *************************************************************************/
int mpz_fibonacci_prp(mpz_t n, long int p, long int q);

/* *******************************************************************************
 * mpz_lucas_prp:
 * A "Lucas pseudoprime" with parameters (P,Q) is a composite n with D=P^2-4Q,
 * (n,2QD)=1 such that U_(n-(D/n)) == 0 mod n [(D/n) is the Jacobi symbol]
 * *******************************************************************************/
int mpz_lucas_prp(mpz_t n, long int p, long int q);

/* *********************************************************************************************
 * mpz_stronglucas_prp:
 * A "strong Lucas pseudoprime" with parameters (P,Q) is a composite n = (2^r)*s+(D/n), where
 * s is odd, D=P^2-4Q, and (n,2QD)=1 such that either U_s == 0 mod n or V_((2^t)*s) == 0 mod n
 * for some t, 0 <= t < r. [(D/n) is the Jacobi symbol]
 * *********************************************************************************************/
int mpz_stronglucas_prp(mpz_t n, long int p, long int q);

/* *******************************************************************************************
 * mpz_extrastronglucas_prp:
 * Let U_n = LucasU(p,1), V_n = LucasV(p,1), and D=p^2-4.
 * An "extra strong Lucas pseudoprime" to the base p is a composite n = (2^r)*s+(D/n), where
 * s is odd and (n,2D)=1, such that either U_s == 0 mod n and V_s == +/-2 mod n, or
 * V_((2^t)*s) == 0 mod n for some t with 0 <= t < r-1 [(D/n) is the Jacobi symbol]
 * *******************************************************************************************/
int mpz_extrastronglucas_prp(mpz_t n, long int p);

/* ***********************************************************************************************
 * mpz_selfridge_prp:
 * A "Lucas-Selfridge pseudoprime" n is a "Lucas pseudoprime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the Lucas pseudoprime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * ***********************************************************************************************/
int mpz_selfridge_prp(mpz_t n);

/* *********************************************************************************************************
 * mpz_strongselfridge_prp:
 * A "strong Lucas-Selfridge pseudoprime" n is a "strong Lucas pseudoprime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the strong Lucase pseudoprime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * **********************************************************************************************************/
int mpz_strongselfridge_prp(mpz_t n);

/* **********************************************************************************
 * mpz_bpsw_prp:
 * A "Baillie-Pomerance-Selfridge-Wagstaff pseudoprime" is a composite n such that
 * n is a strong pseudoprime to the base 2 and
 * n is a Lucas pseudoprime using the Selfridge parameters.
 * **********************************************************************************/
int mpz_bpsw_prp(mpz_t n);

/* ****************************************************************************************
 * mpz_strongbpsw_prp:
 * A "strong Baillie-Pomerance-Selfridge-Wagstaff pseudoprime" is a composite n such that
 * n is a strong pseudoprime to the base 2 and
 * n is a strong Lucas pseudoprime using the Selfridge parameters.
 * ****************************************************************************************/
int mpz_strongbpsw_prp(mpz_t n);






/*******************************************************/
/*******************************************************/
/* These are the definitions for the APRT-CLE routines */
/*******************************************************/
/*******************************************************/
/* verbose = 0 means to only return the status  */
/*          it will not print out progress info */
/* verbose = 1 means to print out progress info */
/* verbose = 2 means to print out progress/failure/failover info */
#define APRTCLE_VERBOSE0 0
#define APRTCLE_VERBOSE1 1
#define APRTCLE_VERBOSE2 2

#define APRTCLE_ERROR -1
#define APRTCLE_COMPOSITE 0
#define APRTCLE_PRP 1
#define APRTCLE_PRIME 2

/* **********************************************************************************
 * APR-CL (also known as APRT-CLE) is a prime proving algorithm developed by:
 * L. Adleman, C. Pomerance, R. Rumely, H. Cohen, and H. W. Lenstra
 * APRT-CLE = Adleman-Pomerance-Rumely Test Cohen-Lenstra Extended version
 * You can find all the details of this implementation in the Cohen & Lenstra paper:
 *    H. Cohen and A. K. Lenstra, "Implementation of a new primality test",
 *    Math. Comp., 48 (1987) 103--121
 *
 * ----------------------------------------------------------------------------------
 *
 * This C/GMP version is a conversion of Dario Alpern's Java based APRT-CLE code
 * His code was based on Yuji Kida's UBASIC code
 *
 * Based on APRT-CLE Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
 * From: Last updated September 10th, 2011. See http://www.alpertron.com.ar/ECM.HTM
 *
 * On 2012/11/12 Dario Alpern has approved the conversion, from Java to C/GMP, of 
 * his implementation of the APR-CL algorithm, and that it be licensed under the LGPL.
 *
 * ----------------------------------------------------------------------------------
 *
 * With improvements based on Jason Moxham's APRCL v1.15 code, from 2003/01/01
 *
 * On 2013/04/14 Toby Moxham has approved the APR-CL code and data tables, 
 * originally written by his brother Jason Moxham on 2003/01/01, to be released 
 * under the LGPL.
 *
 * ----------------------------------------------------------------------------------
 *
 * v1.0 to v1.1 improvements:
 *
 * [The following fix was recommended by Dana Jacobsen and verified by Jon Grantham]
 *      - Bug fix: Removed unnecessary vl==0 check in mpz_extrastronglucas_prp
 * [The following improvements/fixes were recommended by Laurent Desnogues in 2013/08]
 *      - Speed improvement 1: Removed extraneous NormalizeJS calls in ARPCL
 *      - Speed improvement 2: Removed/consolidated calls to mpz_mod in APRCL
 *        (these improvements make the APRCL code about 1.5-2.2x faster)
 *      - Bug fix: Final test in APRCL routine is now correct
 *
 *
 * *********************************************************************************/

int mpz_aprcl(mpz_t N); /* Just return the status of the input, no progress is printed out */
int mpz_aprtcle(mpz_t N, int verbose);

#endif
