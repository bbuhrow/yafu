/*2:*/
#line 111 "if.w"


#include <time.h> 
#include <unistd.h> 
#ifndef _WIN64 
#include <sys/times.h> 
#endif
#line 118 "if.w"
#include <stdio.h> 
#include <stdlib.h> 
#include <stdarg.h> 
#include <string.h> 
#include <gmp.h> 
#include <limits.h> 

#include "asm/siever-config.h"
#include "if.h"

int verbose= 0;
static unsigned int used_cols,ncol= 80;

/*:2*//*3:*/
#line 132 "if.w"

void*xmalloc(size_t size)
{
char*x;
if(size==0)return NULL;
if((x= malloc(size))==NULL)complain("xmalloc: %m\n");
return x;
}

/*:3*//*4:*/
#line 142 "if.w"

#if defined( _MSC_VER ) && defined( _WIN64 )
#error _MSC_VER
void*xvalloc(size_t size)
{
char*x;
static long g_pagesize= 0;
if(!g_pagesize)
{
SYSTEM_INFO system_info;
GetSystemInfo(&system_info);
g_pagesize= system_info.dwPageSize;
}
if(size==0)return NULL;
if((x= _aligned_malloc(size,g_pagesize))==NULL)complain("xvalloc: %m\n");
return x;
}

#elif defined (_WIN64)
#line 161 "if.w"

void*xvalloc(size_t size)
{
char*x;
if(size==0)return NULL;
if((x= _aligned_malloc(size,4096))==NULL)complain("xvalloc: %m\n");
return x;
}

#else
#line 171 "if.w"

void*xvalloc(size_t size)
{
char*x;
if(size==0)return NULL;
if((x= valloc(size))==NULL)complain("xvalloc: %m\n");
return x;
}
#endif
#line 180 "if.w"

/*:4*//*5:*/
#line 182 "if.w"

void*
xcalloc(size_t n,size_t s)
{
void*p;
if(n==0||s==0)return NULL;
if((p= calloc(n,s))==NULL)complain("calloc: %m\n");
return p;
}

/*:5*//*6:*/
#line 193 "if.w"

void*xrealloc(void*x,size_t size)
{
char*y;
if(size==0){
if(x!=NULL)free(x);
return NULL;
}
if((y= realloc(x,size))==NULL&&size!=0)complain("xrealloc: %m\n");
return y;
}

/*:6*//*7:*/
#line 207 "if.w"

FILE*logfile= NULL;



#define MAXMSGLEN 20000
void complain(char*fmt,...)
{
char msg[MAXMSGLEN];

va_list arglist;
va_start(arglist,fmt);

int qlen= vsnprintf(msg,MAXMSGLEN,fmt,arglist);

if(qlen<0||qlen> MAXMSGLEN){
fprintf(stderr,"ERROR: vsnprintf call failed during complain\n");
}else{
#if defined (__APPLE__) || defined (_WIN64)
int lenfmt= strlen(fmt);
if(lenfmt> 4&&fmt[lenfmt-3]=='%'&&fmt[lenfmt-2]=='m'&&fmt[lenfmt-1]=='\n'){

msg[qlen-2]= '\0';

snprintf(msg,MAXMSGLEN,"%s%s\n",msg,strerror(errno));
}
#endif
#line 234 "if.w"
}


fprintf(stderr,"%s",msg);

va_end(arglist);

if(logfile!=NULL){
fprintf(logfile,"%s",msg);
}
#ifdef HAVE_BOINC
boinc_finish(1);
#else
#line 247 "if.w"
 exit(1);
#endif
#line 249 "if.w"
}

/*:7*//*8:*/
#line 252 "if.w"

void Schlendrian(char*fmt,...)
{
va_list arglist;
va_start(arglist,fmt);
vfprintf(stderr,fmt,arglist);
if(logfile!=NULL)vfprintf(logfile,fmt,arglist);
#ifdef HAVE_BOINC
boinc_finish(1);
#else
#line 262 "if.w"
 abort();
#endif
#line 264 "if.w"
}

/*:8*//*9:*/
#line 267 "if.w"

#define NEW_NUMBER -0x10000
void logbook(int l,char*fmt,...)
{
if(l==NEW_NUMBER){
used_cols= 0;
return;
}
if(l<verbose){
va_list arglist;
char*output_str;
unsigned int sl;
va_start(arglist,fmt);
vasprintf(&output_str,fmt,arglist);
sl= strlen(output_str);
if(used_cols+sl> ncol){
fprintf(stderr,"\n");
if(logfile!=NULL)fprintf(logfile,"\n");
used_cols= 0;
}
fputs(output_str,stderr);
if(logfile!=NULL)fputs(output_str,logfile);
if(output_str[sl-1]=='\n')used_cols= 0;
else used_cols+= sl;
free(output_str);
}
}

/*:9*//*10:*/
#line 297 "if.w"

int
errprintf(char*fmt,...)
{
va_list arglist;
int res;

va_start(arglist,fmt);
if(logfile!=NULL)vfprintf(logfile,fmt,arglist);
res= vfprintf(stderr,fmt,arglist);
return res;
}

/*:10*//*11:*/
#line 311 "if.w"

void adjust_bufsize(void**buf,size_t*alloc,size_t req,
size_t incr,size_t item_size)
{
if(req> *alloc){
size_t new_alloc;
new_alloc= *alloc+incr*((req+incr-1-*alloc)/incr);
if(*alloc> 0)*buf= xrealloc(*buf,new_alloc*item_size);
else*buf= xmalloc(new_alloc*item_size);
*alloc= new_alloc;
}
}

/*:11*//*12:*/
#line 327 "if.w"

int yn_query(char*fmt,...)
{
va_list arglist;
char answer[10];

va_start(arglist,fmt);
if(logfile!=NULL)vfprintf(logfile,fmt,arglist);
vfprintf(stderr,fmt,arglist);
if(!isatty(STDIN_FILENO)||!isatty(STDERR_FILENO))return 0;

fflush(stderr);
while(scanf("%9s",answer)!=1||
(strcasecmp(answer,"yes")!=0&&strcasecmp(answer,"no")!=0)){
fprintf(stderr,"Please answer yes or no!\n");
vfprintf(stderr,fmt,arglist);
}
if(strcasecmp(answer,"yes")==0)return 1;
return 0;
}


/*:12*//*13:*/
#line 356 "if.w"

ssize_t
skip_blanks_comments(char**iline,size_t*iline_alloc,FILE*ifi)
{
while(getline(iline,iline_alloc,ifi)> 0){
if(**iline!='#'&&strspn(*iline,"\n\t ")<strlen(*iline))
return 1;
}
return 0;
}

/*:13*//*14:*/
#line 368 "if.w"

#ifdef BIGENDIAN
/*15:*/
#line 375 "if.w"

static u32_t
bswap_32(u32_t x)
{
return((x&0x000000ffUL)<<24)|((x&0x0000ff00UL)<<8)|
((x&0x00ff0000UL)>>8)|((x&0xff000000UL)>>24);
}

/*:15*//*16:*/
#line 384 "if.w"

static u64_t
bswap_64(u64_t x)
{
return((x&0xffULL)<<56)|((x&0xff00ULL)<<40)|((x&0xff0000ULL)<<24)|
((x&0xff000000ULL)<<8)|((x&0xff00000000ULL)>>8)|
((x&0xff0000000000ULL)>>24)|((x&0xff000000000000ULL)>>40)|
((x&0xff00000000000000ULL)>>56);
}

/*:16*//*17:*/
#line 395 "if.w"

int
write_i64(FILE*ofile,i64_t*buffer,size_t count)
{
size_t i;
int res;

for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:17*//*18:*/
#line 411 "if.w"

int
write_u64(FILE*ofile,u64_t*buffer,size_t count)
{
size_t i;
int res;

for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:18*//*19:*/
#line 427 "if.w"

int
write_i32(FILE*ofile,i32_t*buffer,size_t count)
{
size_t i;
int res;

for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}

/*:19*//*20:*/
#line 443 "if.w"

int
write_u32(FILE*ofile,u32_t*buffer,size_t count)
{
size_t i;
int res;

for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
res= fwrite(buffer,sizeof(*buffer),count,ofile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}

/*:20*//*21:*/
#line 459 "if.w"

int
read_i64(FILE*ifile,i64_t*buffer,size_t count)
{
size_t i;
int res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:21*//*22:*/
#line 473 "if.w"

int
read_u64(FILE*ifile,u64_t*buffer,size_t count)
{
size_t i;
int res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_64(buffer[i]);
return res;
}

/*:22*//*23:*/
#line 487 "if.w"

int
read_i32(FILE*ifile,i32_t*buffer,size_t count)
{
size_t i;
int res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}

/*:23*//*24:*/
#line 501 "if.w"

int
read_u32(FILE*ifile,u32_t*buffer,size_t count)
{
size_t i;
int res;

res= fread(buffer,sizeof(*buffer),count,ifile);
for(i= 0;i<count;i++)
buffer[i]= bswap_32(buffer[i]);
return res;
}

/*:24*/
#line 370 "if.w"

#endif 
#line 372 "if.w"


/*:14*//*25:*/
#line 516 "if.w"

#ifdef NEED_GETLINE
#define GETL_INCR 128
ssize_t
getline(char**lineptr,size_t*n,FILE*stream)
{
int rv;

if(*n==0){
*n= GETL_INCR;
*lineptr= xmalloc(*n);
}

rv= 0;
for(;;){
int m;
m= *n-rv;
if(fgets(*lineptr+rv,m-1,stream)==NULL)break;
rv= strlen(*lineptr);
if(rv==0||(*lineptr)[rv-1]=='\n')break;
*n+= GETL_INCR;
*lineptr= xrealloc(*lineptr,*n);
}
return rv;
}
#endif
#line 542 "if.w"

/*:25*//*26:*/
#line 544 "if.w"

#ifdef NEED_ASPRINTF
int
vasprintf(char**ptr,const char*template,va_list ap)
{

int n,size= 32;

while(1){
*ptr= xmalloc(size);

n= vsnprintf(*ptr,size,template,ap);
if(1+strlen(*ptr)<size)return n;
size*= 2;
free(*ptr);
}
}
#endif
#line 562 "if.w"

/*:26*//*27:*/
#line 564 "if.w"

#ifdef NEED_ASPRINTF
int asprintf(char**ptr,const char*template,...)
{
int rv;
va_list ap;

va_start(ap,template);
rv= vasprintf(ptr,template,ap);
va_end(ap);
return rv;
}
#endif
#line 577 "if.w"

/*:27*//*28:*/
#line 579 "if.w"

#ifdef NEED_FNMATCH
int fnmatch(char*s,char*fname,int dummy)
{
while(*fname)fname++;
if(*(fname--)!='z')return 1;
if(*(fname--)!='g')return 1;
if(*(fname--)!='.')return 1;
return 0;
}
#endif
#line 590 "if.w"


/*:28*//*29:*/
#line 594 "if.w"

int
u32_cmp012(const void*x,const void*y)
{
const u32_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return-1;
if(*xx> *yy)return 1;
return 0;
}

/*:29*//*30:*/
#line 607 "if.w"

int
u32_cmp210(const void*x,const void*y)
{
const u32_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return 1;
if(*xx> *yy)return-1;
return 0;
}

/*:30*//*31:*/
#line 620 "if.w"

int
u64_cmp012(const void*x,const void*y)
{
const u64_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return-1;
if(*xx> *yy)return 1;
return 0;
}

/*:31*//*32:*/
#line 633 "if.w"

int
u64_cmp210(const void*x,const void*y)
{
const u64_t*xx,*yy;
xx= x;
yy= y;
if(*xx<*yy)return 1;
if(*xx> *yy)return-1;
return 0;
}
/*:32*/
