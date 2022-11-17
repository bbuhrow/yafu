#define NMAX_ULONGS   8
// SMJS #define ulong  unsigned long

// SMJS For ulong type
#include "siever-config.h"
/*
#ifdef _WIN64
typedef  unsigned long long ulong;
#else
typedef  unsigned long ulong;
#endif
*/

void ASM_ATTR (*asm_mulmod)(ulong *,ulong *,ulong *);
void ASM_ATTR (*asm_zero)(ulong *);
void ASM_ATTR (*asm_copy)(ulong *,ulong *);
void ASM_ATTR (*asm_half)(ulong *);
void ASM_ATTR (*asm_sub)(ulong *,ulong *,ulong *);
void ASM_ATTR (*asm_add2)(ulong *,ulong *);

void (*asm_sub_n)(ulong *,ulong *);

void (*asm_squmod)(ulong *,ulong *);
void (*asm_diff)(ulong *,ulong *,ulong *);
void (*asm_add2_ui)(ulong *,ulong);
int (*asm_inv)(ulong *,ulong *);

void init_montgomery_multiplication();
int set_montgomery_multiplication(mpz_t);

// SMJS Moved protos to here from c file so can be used in noasm64.c as well
// SMJS Removed externs
void ASM_ATTR asm_mulm64(ulong *,ulong *,ulong *);
void ASM_ATTR asm_zero64(ulong *);
void ASM_ATTR asm_sub64_3(ulong *,ulong *,ulong *);
void ASM_ATTR asm_copy64(ulong *,ulong *);
void ASM_ATTR asm_half64(ulong *);
void ASM_ATTR asm_add64(ulong *,ulong *);

void asm_sub_n64(ulong *,ulong *);
void asm_add64_ui(ulong *,ulong);
void asm_sqm64(ulong *,ulong *);
void asm_diff64(ulong *,ulong *,ulong *);
int asm_inv64(ulong *,ulong *);


void ASM_ATTR asm_sub128_3(ulong *,ulong *,ulong *);
void ASM_ATTR asm_zero128(ulong *);
void ASM_ATTR asm_copy128(ulong *,ulong *);
void ASM_ATTR asm_half128(ulong *);
void ASM_ATTR asm_mulm128(ulong *,ulong *,ulong *);
void ASM_ATTR asm_add128(ulong *,ulong *);

void asm_sub_n128(ulong *,ulong *);
void asm_sqm128(ulong *,ulong *);
void asm_diff128(ulong *,ulong *,ulong *);
void asm_add128_ui(ulong *,ulong);
int asm_inv128(ulong *,ulong *);


void ASM_ATTR asm_zero192(ulong *);
void ASM_ATTR asm_mulm192(ulong *,ulong *,ulong *);
void ASM_ATTR asm_sub192_3(ulong *,ulong *,ulong *);
void ASM_ATTR asm_add192(ulong *,ulong *);
void ASM_ATTR asm_copy192(ulong *,ulong *);
void ASM_ATTR asm_half192(ulong *);

void asm_sub_n192(ulong *,ulong *);
void asm_sqm192(ulong *,ulong *);
void asm_diff192(ulong *,ulong *,ulong *);
void asm_add192_ui(ulong *,ulong);
int asm_inv192(ulong *,ulong *);
