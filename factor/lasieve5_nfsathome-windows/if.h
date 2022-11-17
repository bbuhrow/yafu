/*1:*/
#line 18 "if.w"

#ifdef _WIN64
#define NEED_ASPRINTF
#define NEED_FNMATCH
#define NEED_GETLINE
#endif
#line 24 "if.w"

#include <stdarg.h> 
#include <stdio.h> 
#include <unistd.h> 
#include <errno.h> 

#ifdef HAVE_BOINC
#include<stdarg.h> 
#ifdef _WIN32
#include"boinc_win.h"
#include"boinc_api.h"
#include"filesys.h"
#else
#include"boinc_api.h"
#include"filesys.h"
#endif
#endif
#line 41 "if.w"

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
#line 63 "if.w"
#define write_i64(ofile,buffer,count) fwrite((void*)buffer,sizeof(i64_t),count,ofile)
#define write_u64(ofile,buffer,count) fwrite((void*)buffer,sizeof(u64_t),count,ofile)
#define write_u32(ofile,buffer,count) fwrite((void*)buffer,sizeof(u32_t),count,ofile)
#define write_i32(ofile,buffer,count) fwrite((void*)buffer,sizeof(i32_t),count,ofile)
#define read_i64(ofile,buffer,count) fread((void*)buffer,sizeof(i64_t),count,ofile)
#define read_u64(ofile,buffer,count) fread((void*)buffer,sizeof(u64_t),count,ofile)
#define read_u32(ofile,buffer,count) fread((void*)buffer,sizeof(u32_t),count,ofile)
#define read_i32(ofile,buffer,count) fread((void*)buffer,sizeof(i32_t),count,ofile)
#endif 
#line 72 "if.w"
int yn_query(char*fmt,...);
ssize_t skip_blanks_comments(char**,size_t*,FILE*);
int u32_cmp012(const void*,const void*);
int u32_cmp210(const void*,const void*);
int u64_cmp012(const void*,const void*);
int u64_cmp210(const void*,const void*);

#ifdef NEED_ASPRINTF
int vasprintf(char**ptr,const char*template,va_list ap);
int asprintf(char**,const char*,...);
#endif
#line 83 "if.w"

#ifdef NEED_GETLINE
ssize_t getline(char**,size_t*,FILE*);
#endif
#line 87 "if.w"

#ifdef NEED_FNMATCH
int fnmatch(char*,char*,int);
#endif
#line 91 "if.w"

typedef unsigned long long ullong;


#ifdef _WIN64
#include <inttypes.h> 
#define UL_FMTSTR "%"PRIu64
#define DLL_FMTSTR "%"PRId64
#define UL_XFMTSTR "%"PRIX64
#define UL_xFMTSTR "%"PRIx64
#define SIZET_FMTSTR "%zu"
#else
#line 103 "if.w"
#define UL_FMTSTR "%lu"
#define DLL_FMTSTR "%ld"
#define UL_XFMTSTR "%lX"
#define UL_xFMTSTR "%lx"
#define SIZET_FMTSTR "%zu"
#endif
#line 109 "if.w"

/*:1*/
