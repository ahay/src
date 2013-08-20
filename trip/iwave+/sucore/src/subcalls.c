/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/***************************************************************************
SUBCALLS - routines for system functions with error checking
 
efopen		fopen with error check
efreopen	freopen with error check
efdopen		fdopen with error check
epopen		popen with error check
efclose		fclose with error check
epclose		pclose with error check
efflush		fflush with error check
eremove		remove with error check
erename		rename with error check
efseek		fseek with error check
efseeko		fseeko with error check
eftell		ftell with error check
eftello		ftello with error check
etmpfile	tmpfile with error check
erewind		rewind (dummy sub)
emkstemp		mkstemp with error check
emalloc		malloc with error check
erealloc	realloc with error check
ecalloc		calloc with error check
efgetpos	fgetpos with error check
efsetpos	fsetpos with error check
efread		fread with error check
efwrite		fwrite with error check

****************************************************************************
Function Prototypes:
FILE *efopen(const char *file, const char *mode);
FILE *efreopen(const char *file, const char *mode, FILE *stream1);
FILE *efdopen(int fd, const char *mode);
FILE *epopen(char *command, char *type);
int efclose(FILE *stream);
int epclose(FILE *stream);
int efflush(FILE *stream);
int eremove(const char *file);
int erename(const char *oldfile, const char* newfile);
int efseek(FILE *stream, off_t offset, int origin);
void erewind(FILE *stream);
long eftell(FILE *stream);
off_t eftello(FILE *streem);
int efseek(FILE *stream, off_t offset, int origin);
FILE *etmpfile(void);
void *emalloc(size_t size);
void *erealloc(void *memptr, size_t size);
void *ecalloc(size_t count, size_t size);
int efgetpos(FILE *stream, fpos_t *position);
int efsetpos(FILE *stream, const fpos_t *position);
size_t efread(void *bufptr, size_t size, size_t count, FILE *stream);
size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream);

****************************************************************************
Returns: All abort with message on error

efopen returns a FILE pointer
efreopen returns a FILE pointer
efdopen returns a FILE pointer
efclose returns 0
efflush returns 0
eremove returns 0
erename returns 0
efseek returns 0
efseeko returns 0
erewind is void
eftell returns file position indicator
eftello returns file position indicator
etmpfile returns a FILE pointer
emalloc returns void pointer to allocated memory
erealloc returns void pointer to allocated memory
ecalloc returns void pointer to allocated memory
efgetpos returns 0
efsetpos returns 0
efread returns number of items actually read
efwrite returns number of items actually written

****************************************************************************
Notes:
Getting less than the number of bytes asked for on a fread
is *not* a system subroutine error--usually it just means
end of file has been reached.  However, it *might* be an error
in some application.  Similarly coming up empty is not a system
error, but might be an application error.  It is left to the user
to trap these instances.  In particular, feof can be used to test
for end of file.

****************************************************************************
References: 
Rochkind, "Advanced UNIX Programming"
Kernighan and Pike, "The UNIX Programming Environment"
Kernighan and Ritchie, "The C Programming Language"
Mark Williams Company, "ANSI C--A Lexical Guide"

****************************************************************************
Authors: SEP: Rick Ottolini, Ron, Jon Claerbout, Stew Levin
CWP: Shuki Ronen, Jack Cohen
***************************************************************************/
/**************** end self doc ********************************/

#include "par.h"
#define ERROR	NULL

FILE *efopen(const char *file, const char *mode)
{
	FILE *stream;

	if (ERROR == (stream = fopen(file, mode)))
		suerr("%s: efopen: fopen failed", __FILE__);
	
	return stream;
}


FILE *efreopen(const char *file, const char *mode, FILE *stream1)
{
	FILE *stream2;

	if (ERROR == (stream2 = freopen(file, mode, stream1)))
			suerr("%s: efreopen: freopen failed", __FILE__);
	
	return stream2;
}

FILE *efdopen(int fd, const char *mode)
{
	FILE *stream;

	if (ERROR == (stream = fdopen(fd, mode)))
		      suerr("%s: efdopen: fdopen failed", __FILE__);
	
	return stream;
}


FILE *epopen(char *command, char *type)
{
	FILE *stream;

	if (ERROR == (stream = popen(command, type)))
		      suerr("%s: epopen: popen failed", __FILE__);
	
	return stream;
}


int efclose(FILE *stream)
{
	int status;

	if (EOF == (status = fclose(stream)))
		      suerr("%s: efclose: fclose failed", __FILE__);

	return status;
}


int epclose(FILE *stream)
{
	int status;

	if (EOF == (status = pclose(stream)))
		      suerr("%s: epclose: pclose failed", __FILE__);

	return status;
}


int efflush(FILE *stream)
{
	int status;

	if (EOF == (status = fflush(stream)))
		      suerr("%s: efflush: fflush failed", __FILE__);

	return status;
}


int eremove(const char *file)
{
	int status;

	if ((status = remove(file)))
		syssuerr("%s: eremove: remove failed", __FILE__);

	return status;
}


int erename(const char *oldfile, const char *newfile)
{
	int status;

	if ((status = rename(oldfile, newfile)))
		syssuerr("%s: erename: rename failed", __FILE__);

	return status;
}


int efseek(FILE *stream, off_t offset, int origin)
{
	if (fseek(stream, offset, origin))  /* non-zero => error */
		      suerr("%s: efseek: fseek failed", __FILE__);

	return 0;
}



void erewind(FILE *stream)	/* dummy function */
{
	rewind(stream);
	return;
}


long eftell(FILE *stream)
{
	long position;

	if (-1L == (position = ftell(stream)))
		syssuerr("%s: eftell: ftell failed", __FILE__);

	return position;
}

int efseeko(FILE *stream, off_t offset, int origin)
{

	/* non-zero => error */
	if (fseeko(stream, (off_t) offset, (int) origin))
		suerr("%s : efseeko: fseeko failed", __FILE__);

	return 0;
}

off_t eftello(FILE *streem)
{
	off_t eposition;
	off_t test=-1;

	eposition = ftello(streem);
	if (test == eposition) {
		fprintf(stderr,"sizeof(off_t)=%lu\n",
				(unsigned long) sizeof(eposition));
	}
	

	return eposition;
}

FILE *etmpfile(void)
{
	FILE *stream;

	if (ERROR == (stream = tmpfile()))
		      suerr("%s: etmpfile: tmpfile failed", __FILE__);
	
	return stream;
}


void *emalloc(size_t size)
{
	void *memptr;

	if (ERROR == (memptr = malloc(size)))
		suerr("%s : emalloc: malloc failed", __FILE__);
	
	return memptr;
}


void *erealloc(void *memptr, size_t size)
{
	void *newptr;

	if (ERROR == (newptr = realloc(memptr, size)))
		suerr("%s : erealloc: realloc failed", __FILE__);
	
	return newptr;
}


void *ecalloc(size_t count, size_t size)
{
	void *memptr;

	if (ERROR == (memptr = calloc(count, size)))
		suerr("%s : ecalloc: calloc failed", __FILE__);
	
	return memptr;
}

/* fgetpos and fsetpos may not exist on some systems */
/* if you get error messages about these then comment out the next two */
/* subroutine definitions */
/* beginning of fgetpos  and fsetpos block */


#ifndef SUN_A
	int efgetpos(FILE *stream, fpos_t *position)
	{
		int status;

		if ((status = fgetpos(stream, position)))
			syssuerr("%s: efgetpos: fgetpos failed", __FILE__);

		return status;
	}


	int efsetpos(FILE *stream, const fpos_t *position)
	{
		int status;

		if ((status = fsetpos(stream, position)))
			syssuerr("%s: efsetpos: fsetpos failed", __FILE__);

		return status;
	}
#endif /* end of SUN_A */
/* end of fgetpos, fsetpos block */

size_t efread(void *bufptr, size_t size, size_t count, FILE *stream)
{
	size_t nread;

	if (!size) suerr("%s: efread: fread given 0 item size", __FILE__);

	nread = fread(bufptr, size, count, stream);

	if (nread != count && ferror(stream))
		      suerr("%s: efread: fread only %d items of %d",
				__FILE__, nread, count);

	return nread;
}


size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream)
{
	size_t nwrite;

	nwrite = fwrite(bufptr, size, count, stream);

	if (nwrite != count)
		      suerr("%s: efwrite: fwrite only %d items of %d",
				__FILE__, nwrite, count);

	return nwrite;
}



