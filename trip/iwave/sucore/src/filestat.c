/* Copyright (c) Colorado School of Mines, 2006.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/***************************************************************************
FILESTAT - Functions to determine and output the type of a file from file
		 descriptor

filestat - determine type of file from file descriptor
printstat - print the FileType as a string

****************************************************************************
Function Prototypes:
FileType filestat(int fd);
char *printstat(int fd);

****************************************************************************
filestat:
Input:
fd		file descriptor

Returns:	FileType filestat(fd)

****************************************************************************
printstat:
Input:
fd		file descriptor

Returns:	char *printstat(fd)

****************************************************************************
Notes:

Example:
	FileType ftype;
	...
	ftype = filestat(STDOUT)
	if (ftype = TTY) {
		...


BADFILETYPE is the error return and it is up to the calling program to
check for it.

FileType is defined in par.h as:
typedef enum {BADFILETYPE = -1,
   TTY, DISK, DIRECTORY, TAPE, PIPE, FIFO, SOCKET, SYMLINK} FileType;

SOCKETS and SYMLINKS are unlikely to be used.

Rochkind's Advanced Unix Programming assures us that the link count
for a pipe is 0.  But it's an old book.  It seems that PIPES are
sometimes implemented as FIFOs.  In most cases, the number of links for
a pipe is 0, even so.  But on NeXT (at least) the link count is 1.
In the code, I test for PIPE first and FIFO second, so for a PIPE
you'll actually get FIFO on machines that do it NeXT's way.

Portability: the code assumes that /dev/rmt0 and /dev/mt0
are tape devices on your system.  If not, make the obvious changes.
/dev/null is handled as a special case and classified as a disk file.

The check on tapes allows either the raw or buffered version.
This is moot and easily changed.

If new types are added, the typedef "FileType" in par.h must be extended.


****************************************************************************
Authors:
SEP: Einar Kjartansson , Stew Levin
CWP: Jack Cohen
***************************************************************************/
/**************** end self doc ********************************/

#ifndef __USE_BSD
#define __USE_BSD 1
#endif

#include <sys/stat.h>
#include "par.h"


#ifdef __S_IFMT
#define S_IFMT __S_IFMT
#endif
#ifdef __S_IFREG
#define S_IFREG __S_IFREG
#endif
#ifdef __S_IFDIR
#define S_IFDIR __S_IFDIR
#endif
#ifdef __S_IFIFO
#define S_IFIFO __S_IFIFO
#endif
#ifdef __S_IFSOCK
#define S_IFSOCK __S_IFSOCK
#endif
#ifdef __S_IFLNK
#define S_IFLNK __S_IFLNK
#endif

#ifndef __APPLE__
#include <sys/sysmacros.h>
#endif

/* determine type of file (DISK, PIPE, ...) */
FileType filestat(int fd)
{
	struct stat sfd;	/* for passed fd	*/
	struct stat sdn;	/* for /dev/null	*/
	struct stat smt;	/* for tape devices	*/


	if (-1 == fstat(fd, &sfd))  return BADFILETYPE;

	/* UNIX subroutine */
	if (isatty(fd))  return TTY;

	/* Standard stat test for regular file */
	if ((sfd.st_mode & S_IFMT) == S_IFREG) return DISK;

	/* Standard stat test for directory */
	if ((sfd.st_mode & S_IFMT) == S_IFDIR) return DIRECTORY;

	/* Only pipes have 0 links (might be FIFO too, so this test first) */
	if (!sfd.st_nlink) return PIPE;

	/* Standard stat test for FIFO */
	if ((sfd.st_mode & S_IFMT) == S_IFIFO) return FIFO;

	/* Standard stat test for socket */
	if ((sfd.st_mode & S_IFMT) == S_IFSOCK) return SOCKET;

	/* Standard stat test for symbolic link */
	if ((sfd.st_mode & S_IFMT) == S_IFLNK) return SYMLINK;

	/* Detect tape by comparing its major device number to	*/
	/* /dev/rmt0 (this is not quite portable).  		*/
	if (0 == stat("/dev/rmt0", &smt) &&
	     major(sfd.st_rdev) == major(smt.st_rdev)) return TAPE;
	if (0 == stat("/dev/mt0", &smt) &&
	     major(sfd.st_rdev) == major(smt.st_rdev)) return TAPE;

	/* Detect file as /dev/null by its device number and	*/
	/* classify it as a disk file.				*/
	if (0 == stat("/dev/null", &sdn) &&
		sfd.st_rdev == sdn.st_rdev) return DISK;

	/* error return */
	return BADFILETYPE;
}


/* Supply ascii string describing type of file */
char *printstat(int fd)
{
	switch (filestat(fd)) {
		case TTY:	return "TTY";
		case DISK:	return "DISK";
		case DIRECTORY:	return "DIRECTORY";
		case TAPE:	return "TAPE";
		case PIPE:	return "PIPE";
		case FIFO:	return "FIFO";
		case SOCKET:	return "SOCKET";
		case SYMLINK:	return "SYMLINK";
		default:	return "BADFILETYPE";
	}
}


