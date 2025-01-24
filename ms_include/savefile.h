#pragma once



#ifndef _SAVEFILE_H_
#define _SAVEFILE_H_


#ifdef __cplusplus
extern "C" {
#endif


#include <stdint.h>
#include <stdio.h>

#define NO_ZLIB

#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
#include <windows.h>
#include <sys/types.h> 
#include <sys/stat.h>
#endif


#ifdef NO_ZLIB
#define gzFile   FILE
#define gzopen   fopen
#define gzclose  fclose
#define gzeof    feof
#define gzrewind rewind
#define gzprintf fprintf
#define gzputs(f,b)   fprintf(f, "%s", b)
#define gzgets(f,b,l) fgets(b,l,f)
#define gzflush(f,b)  fflush(f)
#else
#include <zlib.h>
#endif

/* structure encapsulating the savefile used in a factorization */

typedef struct {

#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
	HANDLE file_handle;
	uint32_t read_size;
	uint32_t eof;
#else

#ifdef NO_ZLIB
	FILE* fp;
#else
	gzFile* fp;
#endif
	char isCompressed;
	char is_a_FILE;
#endif
	char* name;
	char* buf;
	uint32_t buf_off;
} savefile_t;

/*---------------- SAVEFILE RELATED DECLARATIONS ---------------------*/

#define BIGNUM_BUF_SIZE 500
#define LINE_BUF_SIZE 300
#define SAVEFILE_READ 0x01
#define SAVEFILE_WRITE 0x02
#define SAVEFILE_APPEND 0x04

void savefile_init(savefile_t* s, char* filename);
void savefile_free(savefile_t* s);
void savefile_open(savefile_t* s, uint32_t flags);
void savefile_close(savefile_t* s);
uint32_t savefile_eof(savefile_t* s);
uint32_t savefile_exists(savefile_t* s);
void savefile_rewind(savefile_t* s);
void savefile_read_line(char* buf, size_t max_len, savefile_t* s);
void savefile_write_line(savefile_t* s, char* buf);
void savefile_flush(savefile_t* s);




#ifdef __cplusplus
}
#endif

#endif /* _SAVEFILE_H_ */