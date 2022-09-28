#define USE_MEMCPY

/* have pmullw, pmulhw for 8 16bit integers */
#define HAVE_XMM_MUL

/* for mpqs <=96 Bit */

#define TINY

/* on athlon64 ASM_MPQS_TD needs HAVE_XMM_MUL */
#define ASM_MPQS_TD
#define A64_STYLE_TD
u32_t asm_td(u16_t*,u16_t,u64_t);

#define ASM_MPQS_SIEVE
void asm_sieve(void);

#define ASM_MPQS_EVAL
size_t asm_evaluate(unsigned *,unsigned *,u16_t*,unsigned);

#define ASM_MPQS_SIEVE_INIT
void asm_sieve_init(unsigned char*,u32_t,ushort*,u64_t*,unsigned char*,u32_t);
#define ASM_MPQS_SIEVE_INIT16
void asm_sieve_init16(unsigned char*,u32_t,ushort*,u64_t*,unsigned char*,u32_t);
#define ASM_MPQS_NEXT_POL

#define ASM_MPQS_GAUSS
void asm_gauss(void);

#ifdef TINY
static ushort mpqs_param[14][7]={
 { 30, 3, 2, 5, 7, 30, 768},      /* 44 */
 { 30, 3, 3, 5, 9, 30, 1024},     /* 48 */
 { 35, 3, 3, 5, 10, 40, 1280},    /* 52 */
 { 40, 3, 3, 6, 10, 50, 1536},    /* 56 */
 { 50, 3, 3, 6, 10, 50, 1536},    /* 60 */
 { 60, 3, 3, 6, 11, 60, 1536},    /* 64 */
 { 60, 3, 3, 7, 11, 60, 1792},    /* 68 */
 { 80, 3, 3, 7, 12, 80, 2048},    /* 72 */
 { 90, 3, 3, 7, 13, 90, 2048},    /* 76 */
 { 110, 3, 3, 7, 14, 110, 2560},   /* 80 */
 { 120, 3, 4, 7, 14, 120, 2560},   /* 84 */
 { 140, 3, 4, 7, 15, 140, 3072},   /* 88 */
 { 150, 3, 4, 7, 16, 150, 3072},   /* 92 */
 { 160, 4, 4, 8, 17, 160, 3072}    /* 96 */
};
/* second parameter useless in TINY-variant */
#else
static ushort mpqs_param[14][7]={
 { 40, 3, 2, 4, 11, 16, 16384},
 { 40, 3, 2, 4, 11, 16, 16384},
 { 40, 3, 2, 4, 11, 16, 16384},
 { 50, 3, 2, 4, 12, 16, 16384},
 { 60, 3, 3, 4, 15, 16, 16384},
 { 70, 3, 3, 5, 14, 16, 16384},
 { 80, 3, 3, 5, 14, 16, 16384},
 { 90, 3, 3, 5, 15, 20, 16384},
 { 110, 3, 3, 5, 17, 20, 16384},
 { 120, 3, 3, 5, 19, 20, 16384},
 { 140, 3, 4, 6, 18, 30, 16384},
 { 140, 3, 4, 6, 20, 40, 16384},
 { 160, 3, 4, 6, 21, 50, 16384},
 { 180, 4, 4, 6, 23, 70, 16384}
};
#endif


