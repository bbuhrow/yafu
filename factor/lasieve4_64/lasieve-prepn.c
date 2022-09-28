/*4:*/
#line 19 "lasieve-prepn.w"

#include <sys/types.h> 
#include <math.h> 

#include "asm/siever-config.h"
#include "recurrence6.h"
#include "asm/32bit.h"

void
lasieve_setup(FB,proots,fbsz,a0,a1,b0,b1,ri_ptr)
u32_t*FB,*proots,fbsz,*ri_ptr;
i32_t a0,a1,b0,b1;
{
u32_t b0_ul,b1_ul,absa0,absa1;

if(b0<0||b1<0)
Schlendrian("lasieve_setup called with negative 2-nd coordinate (%d,%d)\n",
b0,b1);
if(fbsz<=0)return;
#ifdef HAVE_ASM_LASIEVE_SETUP
/*7:*/
#line 155 "lasieve-prepn.w"

if(FB[fbsz-1]<FLOAT_SETUP_BOUND1&&
fabs(a0)+FB[fbsz-1]*b0<FLOAT_SETUP_BOUND2&&
fabs(a1)+FB[fbsz-1]*b1<FLOAT_SETUP_BOUND2){
asm_lasieve_setup(FB,proots,fbsz,a0,a1,b0,b1,ri_ptr);
return;
}/*:7*/
#line 39 "lasieve-prepn.w"

#endif
b0_ul= (u32_t)b0;
b1_ul= (u32_t)b1;
if(a0>=0){
absa0= (u32_t)a0;
#define A0MOD0(p) (absa0%p)

#define A0MOD1(p) absa0
if(a1>=0){
absa1= (u32_t)a1;



#define A1MOD0(p) (absa1%p)
#define A1MOD1(p) absa1
/*5:*/
#line 95 "lasieve-prepn.w"

{
u32_t fbi,fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
fbp_bound= absa1<absa0?absa0:absa1;
if(fbp_bound<b0_ul)fbp_bound= b0_ul;
if(fbp_bound<b1_ul)fbp_bound= b1_ul;
for(fbi= 0;fbi<fbsz&&FB[fbi]<=fbp_bound;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)




for(;fbi<fbsz;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
}

/*:5*/
#line 55 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
}else{
u32_t aux;

absa1= (u32_t)(-a1);
#define A1MOD0(p) ((aux= absa1%p)> 0 ? p-aux : 0 )
#define A1MOD1(p) (p-absa1)
/*5:*/
#line 95 "lasieve-prepn.w"

{
u32_t fbi,fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
fbp_bound= absa1<absa0?absa0:absa1;
if(fbp_bound<b0_ul)fbp_bound= b0_ul;
if(fbp_bound<b1_ul)fbp_bound= b1_ul;
for(fbi= 0;fbi<fbsz&&FB[fbi]<=fbp_bound;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)




for(;fbi<fbsz;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
}

/*:5*/
#line 64 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
#undef A0MOD0
#undef A0MOD1
}
}else{
absa0= (u32_t)(-a0);
#define A0MOD0(p) ((aux= absa0%p)> 0 ? p-aux : 0 )
#define A0MOD1(p) (p-absa0)
if(a1>=0){
u32_t aux;
absa1= (u32_t)a1;
#define A1MOD0(p) (absa1%p)
#define A1MOD1(p) absa1
/*5:*/
#line 95 "lasieve-prepn.w"

{
u32_t fbi,fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
fbp_bound= absa1<absa0?absa0:absa1;
if(fbp_bound<b0_ul)fbp_bound= b0_ul;
if(fbp_bound<b1_ul)fbp_bound= b1_ul;
for(fbi= 0;fbi<fbsz&&FB[fbi]<=fbp_bound;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)




for(;fbi<fbsz;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
}

/*:5*/
#line 79 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
}else{
u32_t aux;
absa1= (u32_t)(-a1);
#define A1MOD0(p) ((aux= absa1%p)> 0 ? p-aux : 0 )
#define A1MOD1(p) (p-absa1)
/*5:*/
#line 95 "lasieve-prepn.w"

{
u32_t fbi,fbp_bound;



#define B0MOD(p) (b0_ul%p)
#define B1MOD(p) (b1_ul%p)
#define A0MOD(p) A0MOD0(p)
#define A1MOD(p) A1MOD0(p)
fbp_bound= absa1<absa0?absa0:absa1;
if(fbp_bound<b0_ul)fbp_bound= b0_ul;
if(fbp_bound<b1_ul)fbp_bound= b1_ul;
for(fbi= 0;fbi<fbsz&&FB[fbi]<=fbp_bound;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 108 "lasieve-prepn.w"

#undef A0MOD
#undef A1MOD
#undef B0MOD
#undef B1MOD
#define B0MOD(p) b0_ul
#define B1MOD(p) b1_ul
#define A0MOD(p) A0MOD1(p)
#define A1MOD(p) A1MOD1(p)




for(;fbi<fbsz;fbi++)/*6:*/
#line 129 "lasieve-prepn.w"

{
u32_t x;
modulo32= FB[fbi];

if(proots[fbi]==modulo32){
x= B0MOD(modulo32);
if(x==0)x= modulo32;
else{
x= modmul32(modinv32(x),B1MOD(modulo32));
if(x> 0)x= modulo32-x;
}
}else{
x= modsub32(A0MOD(modulo32),modmul32(proots[fbi],B0MOD(modulo32)));
if(x!=0){
x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[fbi],B1MOD(modulo32)),
A1MOD(modulo32)));
}else{

x= FB[fbi];
}
}
ri_ptr+= get_recurrence_info(ri_ptr,FB[fbi],x);
}

/*:6*/
#line 121 "lasieve-prepn.w"

#undef B0MOD
#undef B1MOD
#undef A0MOD
#undef A1MOD
}

/*:5*/
#line 87 "lasieve-prepn.w"

#undef A1MOD0
#undef A1MOD1
}
}
}

/*:4*/
