/*4:*/
#line 34 "redu2.w"

#include <math.h> 
#include <sys/types.h> 
#include "asm/siever-config.h"
#include "redu2.h"

#ifdef _MSC_VER
__inline double trunc(double x)
{
    return x < 0.0 ? ceil(x) : floor(x);
}

__inline double rint(double x) 
{ 
  double i, r = modf(x, &i); 
  if(r < 0.0) { 
    r += 1; i -= 1; 
  } 
  return (r > 0.5 || r == 0.5 && ((int)i & 1) ? i + 1.0 : i); 
} 
#endif
/*=== IMPORTANT: rint() with a sentinel counter is safe       */
/*    AND        yields more relations than using trunc().    */
/* A careful inspection is left as an excercise to the reader */

i32_t n_iter= 0;
void
reduce2(i32_t*a0_ptr,i32_t*b0_ptr,i32_t*a1_ptr,i32_t*b1_ptr,
i32_t a0,i32_t b0,i32_t a1,i32_t b1,float sigma)
{
float a0sq,a1sq,s;
i32_t j=0;

a0sq= ((float)a0)*a0;
a0sq+= sigma*((float)b0)*b0;
a1sq= ((float)a1)*a1;
a1sq+= sigma*((float)b1)*b1;

for(;;){
s= ((float)a0)*a1;
s+= sigma*((float)b0)*b1;

if(++j>1024) break; /* n_iter++; */
if(a0sq<a1sq){
i32_t k;

k= (i32_t)rint(s/a0sq);
if(k==0)break;
a1-= k*a0;
b1-= k*b0;
a1sq= ((float)a1)*a1;
a1sq+= sigma*((float)b1)*b1;
}else{
i32_t k;

k= (i32_t)rint(s/a1sq);
if(k==0)break;
a0-= k*a1;
b0-= k*b1;
a0sq= ((float)a0)*a0;
a0sq+= sigma*((float)b0)*b0;
}
}
n_iter += j;

if(b0<0){
b0= -b0;
a0= -a0;
}
if(b1<0){
b1= -b1;
a1= -a1;
}
if(a0sq> a1sq){
*a0_ptr= a0;
*b0_ptr= b0;
*a1_ptr= a1;
*b1_ptr= b1;
}else{
*a0_ptr= a1;
*b0_ptr= b1;
*a1_ptr= a0;
*b1_ptr= b0;
}
}/*:4*/
