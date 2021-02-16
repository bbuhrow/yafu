/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 9/24/08

Modified:	Ben Buhrow
Date:		11/24/09
Purpose:	Port into Yafu-1.14.
--------------------------------------------------------------------*/

#include "qs.h"
#include "qs_impl.h"
#include "ytools.h"

/* we need a generic interface for reading and writing lines
   of data to the savefile while a factorization is in progress.
   This is necessary for two reasons: first, early msieve 
   versions would sometimes clobber their savefiles, and some
   users have several machines all write to the same savefile
   in a network directory. When output is manually buffered and 
   then explicitly flushed after writing to disk, most of the
   relations in the savefile will survive under these circumstances.

   The other reason is Windows-specific. Microsoft's C runtime 
   library has a bug that causes writes to files more than 4GB in
   size to fail. Thus, to deal with really large savefiles we have
   to call Win32 API functions directly, and these do not have
   any stdio-like stream functionality. Hence we need a homebrew
   implementation of some of stdio.h for the rest of the library
   to use */

#define SAVEFILE_BUF_SIZE 65536

/*--------------------------------------------------------------------*/
void qs_savefile_init(qs_savefile_t *s, char *savefile_name) {
	
	memset(s, 0, sizeof(qs_savefile_t));

	s->name = "siqs.dat";
	if (savefile_name)
	{
		s->name = savefile_name;
	}
	
	s->buf = (char *)xmalloc((size_t)SAVEFILE_BUF_SIZE);
}

/*--------------------------------------------------------------------*/
void qs_savefile_free(qs_savefile_t *s) {
	
	free(s->buf);
	memset(s, 0, sizeof(qs_savefile_t));
}

/*--------------------------------------------------------------------*/
void qs_savefile_open(qs_savefile_t *s, uint32_t flags) {
	
#if defined(WIN32) || defined(_WIN64)
	DWORD access_arg, open_arg;

	if (flags & SAVEFILE_READ)
		access_arg = GENERIC_READ;
	else
		access_arg = GENERIC_WRITE;

	if (flags & SAVEFILE_READ)
		open_arg = OPEN_EXISTING;
	else if (flags & SAVEFILE_APPEND)
		open_arg = OPEN_ALWAYS;
	else
		open_arg = CREATE_ALWAYS;

	s->file_handle = CreateFile(s->name, 
					access_arg,
					FILE_SHARE_READ |
					FILE_SHARE_WRITE, NULL,
					open_arg,
					FILE_FLAG_SEQUENTIAL_SCAN,
					NULL);

	if (s->file_handle == INVALID_HANDLE_VALUE) {
		printf("error: cannot open '%s'", s->name);
		exit(-1);
	}
	if (flags & SAVEFILE_APPEND) {
		LARGE_INTEGER fileptr;
		fileptr.QuadPart = 0;
		SetFilePointerEx(s->file_handle, 
				fileptr, NULL, FILE_END);
	}
	s->read_size = 0;
	s->eof = 0;

#else
	char *open_string;

	if (flags & SAVEFILE_APPEND)
		open_string = "a+";
	else if ((flags & SAVEFILE_READ) && (flags & SAVEFILE_WRITE))
		open_string = "r+w";
	else if (flags & SAVEFILE_READ)
		open_string = "r";
	else
		open_string = "w";

	s->fp = fopen(s->name, open_string);
	if (s->fp == NULL) {
		printf("error: cannot open '%s'", s->name);
		exit(-1);
	}
#endif

	s->buf_off = 0;
	s->buf[0] = 0;
}

/*--------------------------------------------------------------------*/
void qs_savefile_close(qs_savefile_t *s) {
	
#if defined(WIN32) || defined(_WIN64)
	CloseHandle(s->file_handle);
	s->file_handle = INVALID_HANDLE_VALUE;
#else
	fclose(s->fp);
	s->fp = NULL;
#endif
}

/*--------------------------------------------------------------------*/
uint32_t qs_savefile_eof(qs_savefile_t *s) {
	
#if defined(WIN32) || defined(_WIN64)
	return (s->buf_off == s->read_size && s->eof);
#else
	return feof(s->fp);
#endif
}

/*--------------------------------------------------------------------*/
uint32_t qs_savefile_exists(qs_savefile_t *s) {
	
#if defined(WIN32) || defined(_WIN64)
	struct _stat dummy;
	return (_stat(s->name, &dummy) == 0);
#else
	struct stat dummy;
	return (stat(s->name, &dummy) == 0);
#endif
}

/*--------------------------------------------------------------------*/
void qs_savefile_read_line(char *buf, size_t max_len, qs_savefile_t *s) {

#if defined(WIN32) || defined(_WIN64)
	size_t i, j;
	char *sbuf = s->buf;

	for (i = s->buf_off, j = 0; i < s->read_size && 
				j < max_len - 1; i++, j++) { /* read bytes */
		buf[j] = sbuf[i];
		if (buf[j] == '\n' || buf[j] == '\r') {
			buf[j+1] = 0;
			s->buf_off = i + 1;
			return;
		}
	}
	s->buf_off = i;
	if (i == s->read_size && !s->eof) {	/* sbuf ran out? */
		DWORD num_read;
		ReadFile(s->file_handle, sbuf, 
				SAVEFILE_BUF_SIZE, 
				&num_read, NULL);
		s->read_size = num_read;
		s->buf_off = 0;

		/* set EOF only if previous lines have exhausted sbuf
		   and there are no more bytes in the file */

		if (num_read == 0)
			s->eof = 1;
	}
	for (i = s->buf_off; i < s->read_size && 
				j < max_len - 1; i++, j++) { /* read more */
		buf[j] = sbuf[i];
		if (buf[j] == '\n' || buf[j] == '\r') {
			i++; j++;
			break;
		}
	}
	buf[j] = 0;
	s->buf_off = i;
#else
	fgets(buf, (int)max_len, s->fp);
#endif
}

/*--------------------------------------------------------------------*/
void qs_savefile_write_line(qs_savefile_t *s, char *buf) {

	if (s->buf_off + strlen(buf) + 1 >= SAVEFILE_BUF_SIZE)
		qs_savefile_flush(s);

	s->buf_off += sprintf(s->buf + s->buf_off, "%s", buf);
}

/*--------------------------------------------------------------------*/
void qs_savefile_flush(qs_savefile_t *s) {

#if defined(WIN32) || defined(_WIN64)
	if (s->buf_off) {
		DWORD num_write; /* required because of NULL arg below */
		WriteFile(s->file_handle, s->buf, 
				s->buf_off, &num_write, NULL);
	}
	FlushFileBuffers(s->file_handle);
#else
	fprintf(s->fp, "%s", s->buf);
	fflush(s->fp);
#endif

	s->buf_off = 0;
	s->buf[0] = 0;
}

/*--------------------------------------------------------------------*/
void qs_savefile_rewind(qs_savefile_t *s) {

#if defined(WIN32) || defined(_WIN64)
	LARGE_INTEGER fileptr;
	fileptr.QuadPart = 0;
	SetFilePointerEx(s->file_handle, fileptr, NULL, FILE_BEGIN);
	s->read_size = 0;   /* invalidate buffered data */
	s->buf_off = 0;
	s->eof = 0;
#else
	rewind(s->fp);
#endif
}

