@* Our interface to the operating system.

Copyright (C) 2000 Jens Franke.
This file is part of mpqs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.


The x-functions perform what their
normal counterpart does and exit in the case of an error.
We also offer two ways of exiting if a bug is found, and a way
to report how things are proceeding.

@(if.h@>=
#ifdef _WIN64
#define NEED_ASPRINTF
#define NEED_FNMATCH
#define NEED_GETLINE
#endif

#include <stdarg.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>

#ifdef HAVE_BOINC
        #include <stdarg.h>
        #ifdef _WIN32
                #include "boinc_win.h"
                #include "boinc_api.h"
                #include "filesys.h"
        #else
                #include "boinc_api.h"
                #include "filesys.h"
        #endif
#endif

void *xmalloc(size_t size);
void *xvalloc(size_t size);
void *xcalloc(size_t n,size_t s);
void *xrealloc(void *x,size_t size);
void complain(char *fmt,...);
void Schlendrian(char *fmt,...);
void logbook(int l,char *fmt,...);
int errprintf(char *fmt,...);
void adjust_bufsize(void **,size_t *,size_t,size_t,size_t);
extern int verbose;
extern FILE *logfile;
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
#endif /* littlebig end */
int yn_query(char *fmt,...);
ssize_t skip_blanks_comments(char**,size_t*,FILE*);
int u32_cmp012(const void*,const void*);
int u32_cmp210(const void*,const void*);
int u64_cmp012(const void*,const void*);
int u64_cmp210(const void*,const void*);

#ifdef NEED_ASPRINTF
int vasprintf (char **ptr, const char *template, va_list ap);
int asprintf(char**,const char *,...);
#endif

#ifdef NEED_GETLINE
ssize_t getline(char**,size_t*,FILE*);
#endif

#ifdef NEED_FNMATCH
int fnmatch (char*, char*, int);
#endif

typedef unsigned long long ullong;

// SMJS Added
#ifdef _WIN64
#include <inttypes.h>
#define UL_FMTSTR "%"PRIu64
#define DLL_FMTSTR "%"PRId64
#define UL_XFMTSTR "%"PRIX64
#define UL_xFMTSTR "%"PRIx64
#define SIZET_FMTSTR "%zu"
#else
#define UL_FMTSTR "%lu"
#define DLL_FMTSTR "%ld"
#define UL_XFMTSTR "%lX"
#define UL_xFMTSTR "%lx"
#define SIZET_FMTSTR "%zu"
#endif

@
@c

#include <time.h>
#include <unistd.h>
#ifndef _WIN64 // SMJS
#include <sys/times.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <gmp.h>
#include <limits.h>

#include "asm/siever-config.h"
#include "if.h"

int verbose=0;
static unsigned int used_cols,ncol=80;

@
@c
void *xmalloc(size_t size)
{
  char *x;
  if(size==0) return NULL;
  if((x=malloc(size))==NULL) complain("xmalloc: %m\n");
  return x;
}

@
@c
#if defined( _MSC_VER ) && defined( _WIN64 )
#error _MSC_VER
void*xvalloc(size_t size)
{
    char*x;
    static long g_pagesize = 0;
    if(!g_pagesize)
    {
        SYSTEM_INFO system_info;
        GetSystemInfo (&system_info);
        g_pagesize = system_info.dwPageSize;
    }
    if(size==0)return NULL;
    if((x = _aligned_malloc(size, g_pagesize)) == NULL) complain("xvalloc: %m\n");
    return x;
}

#elif defined (_WIN64)

void*xvalloc(size_t size)
{
  char*x;
  if (size==0) return NULL;
  if((x= _aligned_malloc(size,4096))==NULL) complain("xvalloc: %m\n");
  return x;
}

#else

void *xvalloc(size_t size)
{
  char *x;
  if(size==0) return NULL;
  if((x=valloc(size))==NULL) complain("xvalloc: %m\n");
  return x;
}
#endif

@
@c
void *
xcalloc(size_t n,size_t s)
{
  void *p;
  if(n==0 || s==0) return NULL;
  if((p=calloc(n,s))==NULL) complain("calloc: %m\n");
  return p;
}

@
@c
void *xrealloc(void *x,size_t size)
{
  char *y;
  if(size==0) {
    if(x!=NULL) free(x);
    return NULL;
  }
  if((y=realloc(x,size))==NULL && size!=0) complain("xrealloc: %m\n");
  return y;
}

@ The complain function also uses variable arguments. Complaints always go
to |stderr|.
@c
FILE *logfile=NULL;

// SMJS Hacked to allow machines without %m specifier to work
//      Assumes %m occurs at end of format str (is true for all current messages)
#define MAXMSGLEN 20000
void complain(char *fmt,...)
{
  char msg[MAXMSGLEN];

  va_list arglist;
  va_start(arglist,fmt);

  int qlen = vsnprintf(msg, MAXMSGLEN, fmt,arglist);

  if (qlen < 0 || qlen > MAXMSGLEN) {
    fprintf(stderr, "ERROR: vsnprintf call failed during complain\n");
  } else {
#if defined (__APPLE__) || defined (_WIN64)
    int lenfmt = strlen(fmt);
    if (lenfmt > 4 && fmt[lenfmt-3] == '%' && fmt[lenfmt-2] == 'm' && fmt[lenfmt-1] == '\n') {
// Remove m and newline
      msg[qlen-2] = '\0';
// Add sterror(errno)
      snprintf(msg,MAXMSGLEN,"%s%s\n",msg,strerror(errno));
    }
#endif
  }

  // What ever happened in vsnprintf try to print message
  fprintf(stderr,"%s",msg);

  va_end(arglist);

  if(logfile != NULL) {
    fprintf(logfile,"%s",msg);
  }
#ifdef HAVE_BOINC 
    boinc_finish(1);       
#else                
  exit(1);
#endif            
}

@ Use this function to report programming errors.
@c
void Schlendrian(char *fmt,...)
{
  va_list arglist;
  va_start(arglist,fmt);
  vfprintf(stderr,fmt,arglist);
  if(logfile != NULL) vfprintf(logfile,fmt,arglist);
#ifdef HAVE_BOINC 
    boinc_finish(1);       
#else                
  abort();
#endif          
}

@ Status logs are printed to |stderr| or written to a logbook in memory.
@c
#define NEW_NUMBER -0x10000
void logbook(int l, char *fmt,...)
{
  if(l==NEW_NUMBER) {
    used_cols=0;
    return;
  }
  if(l<verbose) {
    va_list arglist;
    char *output_str;
    unsigned int sl;
    va_start(arglist,fmt);
    vasprintf(&output_str,fmt,arglist);
    sl=strlen(output_str);
    if(used_cols+sl>ncol) {
      fprintf(stderr,"\n");
    if(logfile != NULL) fprintf(logfile,"\n");
      used_cols=0;
    }
    fputs(output_str,stderr);
    if(logfile != NULL) fputs(output_str,logfile);
    if(output_str[sl-1]=='\n') used_cols=0;@+
    else used_cols+=sl;
    free(output_str);
  }
}

@ This is almost the same as |fprintf(stderr,fmt,...)|, but in addition
writes the same output to the logfile, if a logfile exists.
@c
int
errprintf(char *fmt,...)
{
  va_list arglist;
  int res;

  va_start(arglist,fmt);
  if(logfile != NULL) vfprintf(logfile,fmt,arglist);
  res=vfprintf(stderr,fmt,arglist);
  return res;
}

@
@c
void adjust_bufsize(void **buf,size_t *alloc,size_t req,
			   size_t incr,size_t item_size)
{
  if(req>*alloc) {
    size_t new_alloc;
    new_alloc=*alloc+incr*((req+incr-1-*alloc)/incr);
    if(*alloc > 0) *buf=xrealloc(*buf,new_alloc*item_size);
    else *buf=xmalloc(new_alloc*item_size);
    *alloc=new_alloc;
  }
}

@ Return 1 if the user answers 'yes', -1 if he answers 'no', and 0 if
|stdin| is not a terminal. Always write the question asked to |stderr| and
possibly to the logfile.
@c
int yn_query(char *fmt,...)
{
  va_list arglist;
  char answer[10];

  va_start(arglist,fmt);
  if(logfile != NULL) vfprintf(logfile,fmt,arglist);
  vfprintf(stderr,fmt,arglist);
  if(!isatty(STDIN_FILENO) || !isatty(STDERR_FILENO)) return 0;
  /* Maybe we dont need this, but is safer so. */
  fflush(stderr);
  while(scanf("%9s",answer) != 1 ||
	(strcasecmp(answer,"yes") != 0 && strcasecmp(answer,"no") !=0)) {
    fprintf(stderr,"Please answer yes or no!\n");
    vfprintf(stderr,fmt,arglist);
  }
  if(strcasecmp(answer,"yes")==0) return 1;
  return 0;
}


@ A function which skips blank line and comment lines.
If it does not find a non-comment line, it returns |0|. Otherwise,
it return |1| and the first non-blank line is in |*iline| and the number
allocated bytes in |*iline_alloc|.

A line is skiped iff it starts with '#' or consists of blank and tab
characters.
@c
ssize_t
skip_blanks_comments(char **iline,size_t *iline_alloc,FILE *ifi)
{
  while(getline(iline,iline_alloc,ifi)>0) {
    if(**iline != '#' && strspn(*iline,"\n\t ")<strlen(*iline))
      return 1;
  }
  return 0;
}

@
@c
#ifdef BIGENDIAN
@<Big endian stuff@>@;
#endif /* littlebig end */


@
@<Big endian stuff@>=
static u32_t
bswap_32(u32_t x)
{
  return ((x&0x000000ffUL)<<24)|((x&0x0000ff00UL)<<8)|
    ((x&0x00ff0000UL)>>8)|((x&0xff000000UL)>>24);
}

@
@<Big endian stuff@>=
static u64_t
bswap_64(u64_t x)
{
  return ((x&0xffULL)<<56)|((x&0xff00ULL)<<40)|((x&0xff0000ULL)<<24)|
    ((x&0xff000000ULL)<<8)|((x&0xff00000000ULL)>>8)|
    ((x&0xff0000000000ULL)>>24)|((x&0xff000000000000ULL)>>40)|
    ((x&0xff00000000000000ULL)>>56);
}

@
@<Big endian stuff@>=
int
write_i64(FILE *ofile,i64_t *buffer,size_t count)
{
  size_t i;
  int res;

  for(i=0;i<count;i++)
    buffer[i]=bswap_64(buffer[i]);
  res=fwrite(buffer,sizeof(*buffer),count,ofile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_64(buffer[i]);
  return res;
}

@
@<Big endian stuff@>=
int
write_u64(FILE *ofile,u64_t *buffer,size_t count)
{
  size_t i;
  int res;

  for(i=0;i<count;i++)
    buffer[i]=bswap_64(buffer[i]);
  res=fwrite(buffer,sizeof(*buffer),count,ofile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_64(buffer[i]);
  return res;
}

@
@<Big endian stuff@>=
int
write_i32(FILE *ofile,i32_t *buffer,size_t count)
{
  size_t i;
  int res;

  for(i=0;i<count;i++)
    buffer[i]=bswap_32(buffer[i]);
  res=fwrite(buffer,sizeof(*buffer),count,ofile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_32(buffer[i]);
  return res;
}

@
@<Big endian stuff@>=
int
write_u32(FILE *ofile,u32_t *buffer,size_t count)
{
  size_t i;
  int res;

  for(i=0;i<count;i++)
    buffer[i]=bswap_32(buffer[i]);
  res=fwrite(buffer,sizeof(*buffer),count,ofile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_32(buffer[i]);
  return res;
}

@
@<Big endian stuff@>=
int
read_i64(FILE *ifile,i64_t *buffer,size_t count)
{
  size_t i;
  int res;

  res=fread(buffer,sizeof(*buffer),count,ifile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_64(buffer[i]);
  return res;
}

@
@<Big endian stuff@>=
int
read_u64(FILE *ifile,u64_t *buffer,size_t count)
{
  size_t i;
  int res;

  res=fread(buffer,sizeof(*buffer),count,ifile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_64(buffer[i]);
  return res;
}

@
@<Big endian stuff@>=
int
read_i32(FILE *ifile,i32_t *buffer,size_t count)
{
  size_t i;
  int res;

  res=fread(buffer,sizeof(*buffer),count,ifile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_32(buffer[i]);
  return res;
}

@
@<Big endian stuff@>=
int
read_u32(FILE *ifile,u32_t*buffer,size_t count)
{
  size_t i;
  int res;

  res=fread(buffer,sizeof(*buffer),count,ifile);
  for(i=0;i<count;i++)
    buffer[i]=bswap_32(buffer[i]);
  return res;
}

@ Note that, unlike the GNU version, this makeshift replacement is not able
to correctly deal with null characters in a line.
@c
#ifdef NEED_GETLINE
#define GETL_INCR 128
ssize_t
getline(char **lineptr, size_t *n, FILE *stream)
{
  int rv;

  if(*n==0) {
    *n=GETL_INCR;
    *lineptr=xmalloc(*n);
  }

  rv=0;
  for(;;) {
    int m;
    m=*n-rv;
    if(fgets(*lineptr+rv,m-1,stream)==NULL) break;
    rv=strlen(*lineptr);
    if(rv==0 || (*lineptr)[rv-1]=='\n') break;
    *n+=GETL_INCR;
    *lineptr=xrealloc(*lineptr,*n);
  }
  return rv;
}
#endif

@ This code is derived from the LINUX snprintf manual page.
@c
#ifdef NEED_ASPRINTF
int
vasprintf (char **ptr, const char *template, va_list ap)
{
  /* Guess we need no more than 100 bytes. */
  int n, size = 32;

  while (1) {
    *ptr = xmalloc (size);
    /* Try to print in the allocated space. */
    n = vsnprintf (*ptr, size, template, ap);
    if(1+strlen(*ptr)<size) return n;
    size*=2;
    free(*ptr);
  }
}
#endif

@
@c
#ifdef NEED_ASPRINTF
int asprintf (char **ptr, const char *template, ...)
{
  int rv;
  va_list ap;

  va_start(ap,template);
  rv=vasprintf(ptr,template,ap);
  va_end(ap);
  return rv;
}
#endif

@ Special purpose fnmatch replacement written by P. Leyland.
@c
#ifdef NEED_FNMATCH
int fnmatch (char *s, char *fname, int dummy)
{
  while (*fname) fname++;
  if (*(fname--) != 'z') return 1;
  if (*(fname--) != 'g') return 1;
  if (*(fname--) != '.') return 1;
  return 0;
}
#endif


@ Two comparison function for ulong numbers, mainly for use as the fourth
argument of |qsort|. Their names are self-explaining.
@c
int
u32_cmp012(const void *x,const void *y)
{
  const u32_t *xx,*yy;
  xx=x;
  yy=y;
  if(*xx<*yy) return -1;
  if(*xx>*yy) return 1;
  return 0;
}

@
@c
int
u32_cmp210(const void *x,const void *y)
{
  const u32_t *xx,*yy;
  xx=x;
  yy=y;
  if(*xx<*yy) return 1;
  if(*xx>*yy) return -1;
  return 0;
}

@
@c
int
u64_cmp012(const void *x,const void *y)
{
  const u64_t *xx,*yy;
  xx=x;
  yy=y;
  if(*xx<*yy) return -1;
  if(*xx>*yy) return 1;
  return 0;
}

@
@c
int
u64_cmp210(const void *x,const void *y)
{
  const u64_t *xx,*yy;
  xx=x;
  yy=y;
  if(*xx<*yy) return 1;
  if(*xx>*yy) return -1;
  return 0;
}

