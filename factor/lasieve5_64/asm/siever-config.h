/*1:*/
#line 12 "siever-config.w"

#ifndef __SIEVER_CONFIG_H__
#define __SIEVER_CONFIG_H__
#include <sys/types.h> 
#include <gmp.h> 




#define HAVE_CMOV

#ifdef _WIN64
#define ASM_ATTR   __attribute__((sysv_abi))
#else
#line 26 "siever-config.w"
#define ASM_ATTR
#endif
#line 28 "siever-config.w"

#ifdef _WIN64
void bzero(void*,size_t);
#endif
#line 32 "siever-config.w"

int psp(mpz_t n);


#define L1_BITS 15
#define ULONG_RI
typedef unsigned u32_t;
typedef int i32_t;
typedef short int i16_t;
typedef unsigned short u16_t;


#ifdef _WIN64
typedef unsigned long long u64_t;
typedef long long i64_t;
typedef unsigned short ushort;
typedef unsigned long long ulong;
#else
#line 50 "siever-config.w"
 typedef unsigned long u64_t;
typedef long i64_t;
typedef unsigned short ushort;
typedef unsigned long ulong;
#endif
#line 55 "siever-config.w"

#define U32_MAX 0xffffffff
#define I32_MAX INT_MAX


int asm_cmp(ulong*a,ulong*b);

#define HAVE_ASM_GETBC
void ASM_ATTR asm_getbc(u32_t,u32_t,u32_t,u32_t*,u32_t*,u32_t*,u32_t*);
#define ASM_SCHEDSIEVE
void ASM_ATTR schedsieve(unsigned char,unsigned char*,u16_t*,u16_t*);

void ASM_ATTR schedsieve_1(unsigned char,unsigned char*,u16_t*,u16_t*);
#define ASM_SCHEDTDSIEVE2
u16_t**ASM_ATTR tdsieve_sched2buf(u16_t**,u16_t*,unsigned char*,u16_t**,u16_t**);

#define ASM_MPZ_TD
#define PREINVERT
#if 1
void MMX_TdAllocate(int,size_t,size_t);
u16_t*MMX_TdInit(int,u16_t*,u16_t*,u32_t*,int);
void MMX_TdUpdate(int,int);
u32_t*MMX_Td(u32_t*,int,u16_t);
#define MMX_TD
#define MMX_REGW 8
#endif
#line 81 "siever-config.w"
#define ASM_LINESIEVER
u32_t*ASM_ATTR slinie(u16_t*,u16_t*,unsigned char*);
#define ASM_LINESIEVER3
u32_t*ASM_ATTR slinie3(u16_t*,u16_t*,unsigned char*);
#define ASM_LINESIEVER2
u32_t*ASM_ATTR slinie2(u16_t*,u16_t*,unsigned char*);
#define ASM_LINESIEVER1
u32_t*ASM_ATTR slinie1(u16_t*,u16_t*,unsigned char*);
#define ASM_TDSLINIE1
u32_t*ASM_ATTR tdslinie1(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_TDSLINIE2
u32_t*ASM_ATTR tdslinie2(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_TDSLINIE3
u32_t*ASM_ATTR tdslinie3(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_TDSLINIE
u32_t*ASM_ATTR tdslinie(u16_t*,u16_t*,unsigned char*,u32_t**);
#define ASM_SEARCH0
u32_t ASM_ATTR lasieve_search0(unsigned char*,unsigned char*,unsigned char*,
unsigned char*,unsigned char*,u16_t*,unsigned char*);
#define HAVE_MM_LASIEVE_SETUP
#define HAVE_MM_LASIEVE_SETUP2
#define HAVE_MM_LASIEVE_SETUP_64
#define  MAX_FB_PER_P 2
#define ASM_RESCALE


#if 1
#define VERY_LARGE_Q
#endif
#line 110 "siever-config.w"


u32_t*ASM_ATTR asm_lasieve_mm_setup0(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup1(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup2(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup3(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup20(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup21(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup22(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup23(u32_t*,u32_t*,size_t,u32_t,u32_t,u32_t,u32_t,u32_t*);


u32_t*ASM_ATTR asm_lasieve_mm_setup0_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup1_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup2_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup3_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup20_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup21_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup22_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);
u32_t*ASM_ATTR asm_lasieve_mm_setup23_64(u32_t*,u32_t*,size_t,u64_t,u64_t,u64_t,u64_t,u32_t*);


void ASM_ATTR rescale_interval1(unsigned char*,u64_t);
void ASM_ATTR rescale_interval2(unsigned char*,u64_t);


u64_t ASM_ATTR asm_modadd64(u64_t,u64_t);
u64_t ASM_ATTR asm_modmul64(u64_t,u64_t);



/*:1*//*2:*/
#line 142 "siever-config.w"

#define FB_RAS 3

/*:2*//*4:*/
#line 153 "siever-config.w"

#define N_PRIMEBOUNDS 12

#endif
#line 157 "siever-config.w"

/*:4*/
