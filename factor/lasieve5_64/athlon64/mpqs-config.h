
#define uchar unsigned char

#define USE_MEMCPY

/* have pmullw, pmulhw for 8 16bit integers */
#define HAVE_XMM_MUL

/* for mpqs <=96 Bit */

#define TINY

// SMJS Make all these ASM_ATTRs 

/* on athlon64 ASM_MPQS_TD needs HAVE_XMM_MUL */
#define ASM_MPQS_TD
#define A64_STYLE_TD
u32_t ASM_ATTR asm_td(u16_t*,u16_t,u64_t);

#define ASM_MPQS_SIEVE
void ASM_ATTR asm_sieve(void);

#define ASM_MPQS_EVAL
size_t ASM_ATTR asm_evaluate(unsigned *,unsigned *,u16_t*,unsigned);
size_t ASM_ATTR asm_evaluate_xmm(unsigned *,unsigned *,u16_t*,unsigned);

size_t ASM_ATTR asm_evaluate0(unsigned *,unsigned *,u16_t*,unsigned);
size_t ASM_ATTR asm_evaluate0_xmm(unsigned *,unsigned *,u16_t*,unsigned);

#define ASM_MPQS_SIEVE_INIT
void ASM_ATTR asm_sieve_init(unsigned char*,u32_t,ushort*,u64_t*,unsigned char*,u32_t);
#define ASM_MPQS_SIEVE_INIT16
void ASM_ATTR asm_sieve_init16(unsigned char*,u32_t,ushort*,u64_t*,unsigned char*,u32_t);
#define ASM_MPQS_NEXT_POL
void ASM_ATTR asm_next_pol3minus_xmm(u32_t,ushort *);
void ASM_ATTR asm_next_pol3minus(u32_t,ushort *);
void ASM_ATTR asm_next_pol3plus_xmm(u32_t,ushort *);
void ASM_ATTR asm_next_pol3plus(u32_t,ushort *);
void ASM_ATTR asm_next_pol10_xmm(u32_t,ushort *,ushort *,u32_t);
void ASM_ATTR asm_next_pol11_xmm(u32_t);


void ASM_ATTR asm3_next_pol3minus_xmm(u32_t,ushort *);
void ASM_ATTR asm3_next_pol3minus(u32_t,ushort *);
void ASM_ATTR asm3_next_pol3plus_xmm(u32_t,ushort *);
void ASM_ATTR asm3_next_pol3plus(u32_t,ushort *);
void ASM_ATTR asm3_next_pol10_xmm(u32_t,ushort *,ushort *,u32_t);
void ASM_ATTR asm3_next_pol11_xmm(u32_t);



#define ASM_MPQS_GAUSS
void ASM_ATTR asm_gauss(void);

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


/* for mpqs >=96 Bit */


#define TINY3


#define ASM_MPQS3_NEXT_POL
#define ASM_MPQS3_TD
u32_t ASM_ATTR asm3_td(u16_t*,u16_t,u32_t*);
/* to be completed
#define ASM_MPQS3_TDSIEVE
// SMJS Done in case
u32_t ASM_ATTR asm3_tdsieve(u16_t*,u16_t*,u16_t**,u16_t);

*/

#define ASM_MPQS3_SIEVE
void ASM_ATTR asm3_sieve(void);
void ASM_ATTR asm3_sievea(void);

#define ASM_MPQS3_SIEVE_INIT
void ASM_ATTR asm_sieve_init0(uchar *mpqs3_tinyarray,  u64_t *m64, u32_t mpqs3_tiny_prod);

#define ASM_MPQS3_EVAL

#define ASM_MPQS3_GAUSS
void ASM_ATTR asm_re_strip(u64_t *,u32_t,i16_t *,unsigned char*);


#ifdef TINY3
static ushort mpqs3_param[15][7]={ /* TODO: optimize parameters */
 { 160, 0, 4, 7, 17, 160, 4096},    /* 92 */
 { 170, 0, 4, 8, 17, 180, 4096},    /* 96 */
 { 200, 0, 4, 8, 17, 200, 5120},   /* 100 */
 { 220, 0, 4, 8, 19, 220, 5120},   /* 104 */
 { 240, 0, 4, 8, 19, 240, 6144},   /* 108 */
 { 270, 0, 5, 8, 21, 272, 7168},   /* 112 */
 { 310, 0, 5, 9, 21, 312, 8192},   /* 116 */
 { 350, 0, 5, 9, 21, 350, 10240},   /* 120 */
 { 380, 0, 5, 9, 22, 380, 11264},   /* 124 */
 { 410, 0, 5, 9, 22, 410, 12800},   /* 128 bit */

 { 440, 0, 5, 9, 22, 420, 12800},   /* experimental */
 { 500, 0, 6, 10, 22, 500, 16384},
 { 512, 0, 6, 10, 22, 512, 16384},
 { 512, 0, 6, 10, 22, 512, 16384},
 { 512, 0, 7, 10, 22, 512, 16384}
};
/* second parameter useless in TINY3-variant */
#else
static ushort mpqs3_param[15][7]={
 { 160, 4, 4, 6, 21, 50, 16384},
 { 180, 4, 4, 7, 22, 70, 16384},
 { 200, 4, 4, 7, 23, 70, 16384},
 { 220, 4, 5, 7, 23, 70, 16384},
 { 250, 4, 5, 8, 24, 70, 16384},
 { 290, 4, 5, 8, 25, 70, 16384},
 { 330, 5, 5, 9, 25, 70, 16384},
 { 370, 5, 5, 9, 25, 80, 16384},
 { 410, 4, 5, 9, 25, 80, 16384},
 { 440, 4, 5, 9, 25, 90, 16384},   /* 128 bit */

 { 480, 4, 5, 9, 25, 90, 16384},   /* experimental */
 { 540, 4, 6, 11, 25, 90, 16384},
 { 600, 4, 7, 11, 25, 90, 16384},
 { 660, 4, 7, 12, 25, 90, 16384},
 { 720, 4, 8, 12, 25, 90, 16384}
};
#endif

#define MOD3

