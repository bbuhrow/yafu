/*1:*/
#line 18 "../if.w"

#include <stdarg.h> 
#include <stdio.h> 

void*xmalloc(size_t size);
void*xvalloc(size_t size);
void*xcalloc(size_t n,size_t s);
void*xrealloc(void*x,size_t size);
void complain(char*fmt,...);
void Schlendrian(char*fmt,...);
void logbook(int l,char*fmt,...);
int errprintf(char*fmt,...);
void adjust_bufsize(void**,size_t*,size_t,size_t,size_t);
extern int verbose;
extern FILE*logfile;
#ifdef BIGENDIAN
int write_i64(FILE*,i64_t*,size_t);
int write_u64(FILE*,u64_t*,size_t);
int write_i32(FILE*,i32_t*,size_t);
int write_u32(FILE*,u32_t*,size_t);
int read_i64(FILE*,i64_t*,size_t);
int read_u64(FILE*,u64_t*,size_t);
int read_i32(FILE*,i32_t*,size_t);
int read_u32(FILE*,u32_t*,size_t);
#else
#define write_i64(ofile,buffer,count) fwrite((void*)buffer,sizeof(i64_t),count,ofile)
#define write_u64(ofile,buffer,count) fwrite((void*)buffer,sizeof(u64_t),count,ofile)
#define write_u32(ofile,buffer,count) fwrite((void*)buffer,sizeof(u32_t),count,ofile)
#define write_i32(ofile,buffer,count) fwrite((void*)buffer,sizeof(i32_t),count,ofile)
#define read_i64(ofile,buffer,count) fread((void*)buffer,sizeof(i64_t),count,ofile)
#define read_u64(ofile,buffer,count) fread((void*)buffer,sizeof(u64_t),count,ofile)
#define read_u32(ofile,buffer,count) fread((void*)buffer,sizeof(u32_t),count,ofile)
#define read_i32(ofile,buffer,count) fread((void*)buffer,sizeof(i32_t),count,ofile)
#endif 
int yn_query(char*fmt,...);
ssize_t skip_blank_comments(char**,size_t*,FILE*);

#ifdef NEED_ASPRINTF
int asprintf(char**,const char*,...);
#endif

#ifdef NEED_GETLINE
ssize_t getline(char**,size_t*,FILE*);
#endif

#ifdef NEED_FNMATCH
int fnmatch(char*,char*,int)
#endif

typedef unsigned long long ullong;

/*:1*/
