/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: savefile.c 964 2014-05-03 03:30:03Z jasonp_sf $
--------------------------------------------------------------------*/

#include <ms_common.h>

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
void savefile_init(savefile_t *s, char *savefile_name) {
	
	memset(s, 0, sizeof(savefile_t));

	s->name = MSIEVE_DEFAULT_SAVEFILE;
	if (savefile_name)
		s->name = savefile_name;
	
	s->buf = (char *)xmalloc((size_t)SAVEFILE_BUF_SIZE);
}

/*--------------------------------------------------------------------*/
void savefile_free(savefile_t *s) {
	
	free(s->buf);
	memset(s, 0, sizeof(savefile_t));
}

/*--------------------------------------------------------------------*/
void savefile_open(savefile_t *s, uint32 flags) {
	
#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
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
#ifndef NO_ZLIB
	char name_gz[256];
	#if defined(WIN32) || defined(_WIN64)
	struct _stati64 dummy;
	#else
	struct stat dummy;
	#endif
#endif

	if (flags & SAVEFILE_APPEND)
		open_string = "a";
	else if ((flags & SAVEFILE_READ) && (flags & SAVEFILE_WRITE))
		open_string = "r+w";
	else if (flags & SAVEFILE_READ)
		open_string = "r";
	else
		open_string = "w";

	s->is_a_FILE = s->isCompressed = 0;

#ifndef NO_ZLIB
	sprintf(name_gz, "%s.gz", s->name);
	#if defined(WIN32) || defined(_WIN64)
	if (_stati64(name_gz, &dummy) == 0) {
		if (_stati64(s->name, &dummy) == 0) {
	#else
	if (stat(name_gz, &dummy) == 0) {
		if (stat(s->name, &dummy) == 0) {
	#endif
			printf("error: both '%s' and '%s' exist. "
			       "Remove the wrong one and restart\n",
				s->name, name_gz);
			exit(-1);
		}
		s->isCompressed = 1;
		s->fp = gzopen(name_gz, open_string);
		if (s->fp == NULL) {
			printf("error: cannot open '%s'\n", name_gz);
			exit(-1);
		}
		/* fprintf(stderr, "using compressed '%s'\n", name_gz); */
	} else if (flags & SAVEFILE_APPEND) {
		/* Unfortunately, append is not intuitive in zlib */
		/* Note: the .dat file may be a compressed file   */
		/*       we are using UNIX philosophy here:       */
		/*       it is the content, not filename, that matters */
		uint8 header[4];
		FILE *fp;
		int n;

		if((fp = fopen(s->name, "r"))) {
			if((n = fread(header, sizeof(uint8), 3, fp)) && 
		   	   (n != 3 || header[0]!=31 || header[1]!=139 || header[2]!=8))
				s->is_a_FILE = 1; 
			/* exists, non-empty and not gzipped,
			   so we will fopen a FILE to append plainly */
			fclose(fp);
		}
		if (s->is_a_FILE) {
			s->fp = (gzFile *)fopen(s->name, "a");
		} else {
			s->fp = gzopen(s->name, "a");
			s->isCompressed = 1;
		}
	} else
	{
		s->fp = gzopen(s->name, open_string);
	}

#else
	s->fp = fopen(s->name, open_string);
#endif
	
	if (s->fp == NULL) {
		printf("error: cannot open '%s'\n", s->name);
		exit(-1);
	}
#endif

	s->buf_off = 0;
	s->buf[0] = 0;
}

/*--------------------------------------------------------------------*/
void savefile_close(savefile_t *s) {
	
#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
	CloseHandle(s->file_handle);
	s->file_handle = INVALID_HANDLE_VALUE;
#else

#ifdef NO_ZLIB
	fclose(s->fp);
#else
	s->is_a_FILE ? fclose((FILE *)s->fp) : gzclose(s->fp);
#endif
	s->fp = NULL;
#endif
}

/*--------------------------------------------------------------------*/
uint32 savefile_eof(savefile_t *s) {
	
#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
	return (s->buf_off == s->read_size && s->eof);
#else

#if defined(NO_ZLIB)
	feof((FILE*)s->fp);
#else
	return (s->is_a_FILE ? feof((FILE *)s->fp) : gzeof(s->fp));
#endif
#endif
}

/*--------------------------------------------------------------------*/
uint32 savefile_exists(savefile_t *s) {
	
#if defined(WIN32) || defined(_WIN64)
	struct _stati64 dummy;
	return (_stati64(s->name, &dummy) == 0);
#else
	struct stat dummy;
	return (stat(s->name, &dummy) == 0);
#endif
}

/*--------------------------------------------------------------------*/
void savefile_read_line(char *buf, size_t max_len, savefile_t *s) {

#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
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

#if NO_ZLIB
	fgets(s->fp, buf, (int)max_len);
#else
	gzgets(s->fp, buf, (int)max_len);
#endif
#endif
}

/*--------------------------------------------------------------------*/
void savefile_write_line(savefile_t *s, char *buf) {

	if (s->buf_off + strlen(buf) + 1 >= SAVEFILE_BUF_SIZE)
		savefile_flush(s);

	s->buf_off += sprintf(s->buf + s->buf_off, "%s", buf);
}

/*--------------------------------------------------------------------*/
void savefile_flush(savefile_t *s) {

#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
	if (s->buf_off) {
		DWORD num_write; /* required because of NULL arg below */
		WriteFile(s->file_handle, s->buf, 
				s->buf_off, &num_write, NULL);
	}
	FlushFileBuffers(s->file_handle);
#else
	if (s->is_a_FILE) {
		fprintf((FILE *)s->fp, "%s", s->buf);
		fflush((FILE *)s->fp);
	} else {
#ifndef NO_ZLIB
		gzputs(s->fp, s->buf);
#endif
	}
#endif

	s->buf_off = 0;
	s->buf[0] = 0;
}

/*--------------------------------------------------------------------*/
void savefile_rewind(savefile_t *s) {

#if defined(NO_ZLIB) && (defined(WIN32) || defined(_WIN64))
	LARGE_INTEGER fileptr;
	fileptr.QuadPart = 0;
	SetFilePointerEx(s->file_handle, fileptr, NULL, FILE_BEGIN);
	s->read_size = 0;   /* invalidate buffered data */
	s->buf_off = 0;
	s->eof = 0;
#else

#ifdef NO_ZLIB
	rewind((FILE*)s->fp);
#else
	s->is_a_FILE ? rewind((FILE *)s->fp) : gzrewind(s->fp);
#endif
#endif
}

