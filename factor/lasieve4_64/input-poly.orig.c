/*4:*/
#line 17 "input-poly.w"

#include <stdio.h> 
#include <sys/types.h> 
#include <string.h> 
#include <gmp.h> 

#include "asm/siever-config.h"
#include "input-poly.h"
#include "if.h"

void
input_poly(mpz_t N,mpz_t**A,i32_t*adeg_ptr,mpz_t**B,i32_t*bdeg_ptr,mpz_t m,
FILE*input_file)
{
char*input_line= NULL;
size_t input_line_alloc= 0;
i32_t have_m= 0;

if(mpz_inp_str(N,input_file,10)==0)
complain("Cannot read number which is to be factored: %m\n");
/*5:*/
#line 59 "input-poly.w"

*adeg_ptr= -1;
*bdeg_ptr= -1;
while(have_m==0){
i32_t grad;
char*field;

if(skip_blanks_comments(&input_line,&input_line_alloc,input_file)<=0)
complain("Cannot read common root of NFS polynomials from input file\n");
switch(*input_line){
case'X':
if(sscanf(input_line+1,"%d",&grad)==0)
complain("Cannot understand input line %s\n",input_line);
if(grad> *adeg_ptr){
i32_t i;
if(*adeg_ptr>=0)*A= xrealloc(*A,(grad+1)*sizeof(**A));
else*A= xmalloc((grad+1)*sizeof(**A));
for(i= *adeg_ptr+1;i<=grad;i++)
mpz_init_set_ui((*A)[i],0);
*adeg_ptr= grad;
}
strtok(input_line," \t");
field= strtok(NULL," \t");
if(string2mpz((*A)[grad],field,10)!=0)
complain("Cannot understand number %s\n",field);
break;
case'Y':
if(sscanf(input_line+1,"%d",&grad)==0)
complain("Cannot understand input line %s\n",input_line);
if(grad> *bdeg_ptr){
i32_t i;
if(*bdeg_ptr>=0)*B= xrealloc(*B,(grad+1)*sizeof(**B));
else*B= xmalloc((grad+1)*sizeof(**B));
for(i= *bdeg_ptr+1;i<=grad;i++)
mpz_init_set_ui((*B)[i],0);
*bdeg_ptr= grad;
}
strtok(input_line," \t");
field= strtok(NULL," \t");
if(string2mpz((*B)[grad],field,10)!=0)
complain("Cannot understand number %s\n",field);
break;
case'M':
strtok(input_line," \t");
field= strtok(NULL," \t");
if(string2mpz(m,field,10)!=0)
complain("Cannot understand number %s\n",field);
have_m= 1;
break;
}
}

/*:5*/
#line 37 "input-poly.w"

if(*adeg_ptr==-1){
*adeg_ptr= 1;
*A= xmalloc(2*sizeof(**A));
mpz_init_set_ui((*A)[1],1);
mpz_init((*A)[0]);
mpz_neg((*A)[0],m);
}
if(*bdeg_ptr==-1){
*bdeg_ptr= 1;
*B= xmalloc(2*sizeof(**B));
mpz_init_set_ui((*B)[1],1);
mpz_init((*B)[0]);
mpz_neg((*B)[0],m);
}
if(*adeg_ptr==0||*bdeg_ptr==0)
complain("Polynomials of degree zero are not allowed\n");
/*6:*/
#line 112 "input-poly.w"

{
mpz_t x;
i32_t i;

if(mpz_sgn(*(*A+*adeg_ptr))==0){
complain("Leading coefficient (degree %u) vanishes\n",*adeg_ptr);
}
if(mpz_sgn(*(*B+*bdeg_ptr))==0){
complain("Leading coefficient (degree %u) vanishes\n",*bdeg_ptr);
}
for(i= 1,mpz_init_set(x,(*A)[*adeg_ptr]);i<=*adeg_ptr;i++){
mpz_mul(x,x,m);
mpz_add(x,x,(*A)[*adeg_ptr-i]);
}
mpz_fdiv_r(x,x,N);
if(mpz_sgn(x)!=0){
mpz_out_str(stderr,10,m);
complain(" not a root of first poly\n");
}
for(i= 1,mpz_set(x,(*B)[*bdeg_ptr]);i<=*bdeg_ptr;i++){
mpz_mul(x,x,m);
mpz_add(x,x,(*B)[*bdeg_ptr-i]);
}
mpz_fdiv_r(x,x,N);
if(mpz_sgn(x)!=0){
mpz_out_str(stderr,10,m);
complain(" not a root of second poly\n");
}
mpz_clear(x);
}/*:6*/
#line 54 "input-poly.w"

free(input_line);
}

/*:4*/
