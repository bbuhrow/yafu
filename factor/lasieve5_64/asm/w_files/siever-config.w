@* Configuration file for GNFS lattice siever.

Copyright (C) 2002 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@(siever-config.h@>=
#ifndef __SIEVER_CONFIG_H__
#define __SIEVER_CONFIG_H__
#include <sys/types.h>
#include <gmp.h>

// SMJS Added
//#define NO_TD_CLOCK
// SMJS Seems a bit faster with HAVE_CMOV defined
#define HAVE_CMOV

#ifdef _WIN64
#define ASM_ATTR   __attribute__((sysv_abi))
#else
#define ASM_ATTR
#endif

#ifdef _WIN64
void bzero(void*, size_t);
#endif

int psp(mpz_t n);
// End added

#define L1_BITS 15
#define ULONG_RI
typedef unsigned u32_t;
typedef int i32_t;
typedef short int i16_t;
typedef unsigned short u16_t;

// Windows 64 has 32 bit long so 64 bit types are long long
#ifdef _WIN64
typedef unsigned long long u64_t;
typedef long long i64_t;
typedef unsigned short ushort;
typedef unsigned long long ulong;
#else
typedef unsigned long u64_t;
typedef long i64_t;
typedef unsigned short ushort;
typedef unsigned long ulong;
#endif

#define U32_MAX 0xffffffff
#define I32_MAX INT_MAX

// SMJS Not actually asm
int asm_cmp(ulong *a, ulong *b);

#define HAVE_ASM_GETBC
void ASM_ATTR asm_getbc(u32_t,u32_t,u32_t,u32_t*,u32_t*,u32_t*,u32_t*);
#define ASM_SCHEDSIEVE
void ASM_ATTR schedsieve(unsigned char,unsigned char*,u16_t*,u16_t*);
// schedsieve_1 only used in g version
void ASM_ATTR schedsieve_1(unsigned char,unsigned char*,u16_t*,u16_t*);
#define ASM_SCHEDTDSIEVE2
u16_t** ASM_ATTR tdsieve_sched2buf(u16_t**,u16_t*,unsigned char*,u16_t **,u16_t**);

#define ASM_MPZ_TD
#define PREINVERT
#if 1
void MMX_TdAllocate(int,size_t,size_t);
u16_t* MMX_TdInit(int,u16_t*,u16_t*,u32_t*,int);
void MMX_TdUpdate(int,int);
u32_t* MMX_Td(u32_t*,int,u16_t);
#define MMX_TD
#define MMX_REGW 8
#endif
#define ASM_LINESIEVER
u32_t* ASM_ATTR slinie(u16_t*,u16_t*,unsigned char*);
#define ASM_LINESIEVER3
u32_t* ASM_ATTR slinie3(u16_t*,u16_t*,unsigned char*);
#define ASM_LINESIEVER2
u32_t* ASM_ATTR slinie2(u16_t*,u16_t*,unsigned char*);
#define ASM_LINESIEVER1
u32_t* ASM_ATTR slinie1(u16_t*,u16_t*,unsigned char*);
#define ASM_TDSLINIE1
u32_t* ASM_ATTR tdslinie1(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_TDSLINIE2
u32_t* ASM_ATTR tdslinie2(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_TDSLINIE3
u32_t* ASM_ATTR tdslinie3(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_TDSLINIE
u32_t* ASM_ATTR tdslinie(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_SEARCH0
u32_t ASM_ATTR lasieve_search0(unsigned char*,unsigned char*,unsigned char*,
                               unsigned char*,unsigned char*,u16_t*,unsigned char*);
#define HAVE_MM_LASIEVE_SETUP
#define HAVE_MM_LASIEVE_SETUP2
#define HAVE_MM_LASIEVE_SETUP_64
#define  MAX_FB_PER_P 2
#define ASM_RESCALE

/* if coordinates of reduced vectors do not fit into 31 bit: */
#if 1
#define VERY_LARGE_Q
#endif


u32_t * ASM_ATTR asm_lasieve_mm_setup0(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup1(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup2(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup3(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup20(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup21(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup22(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup23(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);

// SMJS Added prototypes for 64 bit ones
u32_t * ASM_ATTR asm_lasieve_mm_setup0_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup1_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup2_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup3_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup20_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup21_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup22_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t * ASM_ATTR asm_lasieve_mm_setup23_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);

// More asm protos
void ASM_ATTR rescale_interval1(unsigned char *, u64_t);
void ASM_ATTR rescale_interval2(unsigned char *, u64_t);


u64_t ASM_ATTR asm_modadd64(u64_t, u64_t);
u64_t ASM_ATTR asm_modmul64(u64_t, u64_t);



@ Read-ahead safety for factor base, recurrence info.
@(siever-config.h@>=
#define FB_RAS 3

@ Machine specific initializations.
@c
static void
siever_init(void)
{
}

@ Bounds for the primes treated by the various types of schedule.
@(siever-config.h@>=
#define N_PRIMEBOUNDS 12

#endif 

@
@c
#ifndef MPQS_ONLY
/*
const ulong schedule_primebounds[N_PRIMEBOUNDS]=
	{0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,0x4000000,
	0x8000000,ULONG_MAX};
*/
const ulong schedule_primebounds[N_PRIMEBOUNDS]=
	{0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,0x4000000,
	0x8000000,0x10000000,0x20000000,0x40000000,0x80000000};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={20,21,22,23,24,25,26,27,28,29,30,32};
#endif
