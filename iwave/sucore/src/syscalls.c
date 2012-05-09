/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*************************************************************************
SYSCALLS - routines for SYSTEM CALLs with error checking

ecreat		creat with error check
efork		fork with error check
eopen		open with error check
eclose		close with error check
eunlink		unlink with error check
elseek		lseek with error check
epipe		pipe with error check
eread		read with error check
ewrite		write with error check
 
**************************************************************************
Function Prototypes:
int ecreat(char *path, int perms);
int efork(void);
int eopen(char *path, int flags, int perms);
int eclose(int fd);
int eunlink(char *path);
off_t elseek(int fd, off_t offset, int origin);
int epipe(int fd[2]);
ssize_t eread(int fd, char *buf, size_t nbytes);
ssize_t ewrite(int fd, char *buf, size_t nbytes);

**************************************************************************
Returns:
ecreat returns a file descriptor
efork returns child pid to parent, 0 to child or -1 on error
eopen returns a file descriptor
eclose returns 0 on success or -1 on error
eunlink returns 0 on success or -1 on error
elseek returns number of bytes actually read
epipe returns 0 on success or -1 on error
eread returns number of bytes actually read
ewrite returns number of bytes actually written

**************************************************************************
Notes:
Last arg to read/write is unsigned int on some ANSI systems,
here we follow K&R, page 170.

Rochkind says creat is superfluous, K&R say its mandatory.
I think Rochkind is right--see TEST program below.

Getting less than the number of bytes asked for on a read
is *not* a system error--usually it just means end of file
has been reached.  However, it *might* be an error in some
application.  Similarly coming up empty is not a system error,
but might be an application error.  It is left to the user to
trap these instances.  Here is an example of the first situation
after Rochkind, page 43:

	#define SIZE sizeof(struct direct)
 	...
		while ((nread=eread(fd, &dlink, SIZE)) == SIZE) {
 			... (process dlink)
 		}
 		switch(nread) {
 		case 0:
 			return;
 		default:
 			suerr("partial read");
 		}
 
	In an application where end of file was an error, we could replace
 	case 0 by :
 		case 0:
 			suerr("EOF");

**************************************************************************
References: 
Rochkind, "Advanced UNIX Programming"
Kernighan and Pike, "The UNIX Programming Environment"
Kernighan and Ritchie, "The C Programming Language"
Mark Williams Company, "ANSI C--A Lexical Guide"

**************************************************************************
Authors: SEP: Rick Otolini, Ron, Jon Claerbout, Stew Levin
CWP: Shuki Ronen, Jack Cohen
*************************************************************************/
/**************** end self doc ********************************/

#include "par.h"
#define ERROR	-1

int ecreat(char *path, int perms)
{
	int fd;

	if (ERROR == (fd = creat(path, perms)))
		syssuerr("%s: ecreat: creat failed", __FILE__);
	
	return fd;
}

int efork(void)
{
	int fd;

	if (ERROR == (fd = fork()))
		syssuerr("%s: efork: fork failed", __FILE__);
	
	return fd;
}

int eopen(char *path, int flags, int perms)
{
	int fd;

	if (ERROR == (fd = open(path, flags, perms)))
		syssuerr("%s: eopen: open failed", __FILE__);
	
	return fd;
}

int eclose(int fd)
{
	int returnmode;

	if (ERROR == (returnmode = close(fd)))
		syssuerr("%s: eclose: close failed", __FILE__);
	
	return returnmode;
}

int eunlink(char *path)
{
	int returnmode;

	if (ERROR == (returnmode = unlink(path)))
		syssuerr("%s: eunlink: unlink failed", __FILE__);
	
	return returnmode;
}

off_t elseek(int fd, off_t offset, int origin)
{
	off_t position;

	if (ERROR == (position = lseek(fd, offset, origin)))
		syssuerr("%s: elseek: lseek failed", __FILE__);
	
	return position;
}


int epipe(int fds[2])
{
	int returnmode;

	if (ERROR == (returnmode = pipe(fds)))
		syssuerr("%s: epipe: pipe failed", __FILE__);
	
	return returnmode;
}


ssize_t eread(int fd, char *buf, size_t nbytes)
{
	
	ssize_t nread;	/* items actually read	*/

	if (ERROR == (nread = read(fd, buf, nbytes)))
		syssuerr("%s: eread: read failed", __FILE__);

	return nread;
}

#include <errno.h>

ssize_t ewrite(int fd, char *buf, size_t nbytes)
{
	ssize_t nwrite;	/* bytes actually written	*/

	if (nbytes != (nwrite = write(fd, buf, nbytes))) {
		if (nbytes != -1) errno = EFBIG;
		syssuerr("%s: ewrite: write failed", __FILE__);
	}

	return nwrite;
}


