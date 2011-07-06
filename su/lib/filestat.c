/* Copyright (c) Colorado School of Mines, 2010.*/
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

#include <unistd.h>
#include <sys/stat.h>

#include "filestat.h"

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

#ifndef major	/* major() is a macro for getting the major device number */
#include <sys/sysmacros.h>
#endif

#ifndef _sf_filestat_h

typedef enum {BADFILETYPE = -1,
	TTY, DISK, DIRECTORY, TAPE, PIPE, FIFO, SOCKET, SYMLINK} FileType;
/*^*/

#endif

FileType filestat(int fd)
/*< determine type of file (DISK, PIPE, ...) >*/
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


#ifdef TEST

/* Test driver for function filestat
 *
 * Here are some tests using filestat to analyse STDIN, STDOUT and
 * the first command line argument along with the expected results:
 *	filestat filestat.c
 *		expect: TTY, TTY, DISK
 *	filestat <filestat.c /usr/local | cat
 *		expect: DISK, PIPE, DIRECTORY
 *	cat filestat.c | filestat filestat.c >/usr/tmp/junkxxx
 *		expect: PIPE, DISK, DISK
 *	filestat /dev/null
 *		expect: TTY, TTY, DISK
 *	filestat
 *		expect: TTY, TTY, "no filename given"
 */

main(int argc, char **argv)
{
	static FileType ftype;
	int fd;


	fprintf(stderr, "Checking filestat ...\n");
	/* Get FileType of stdin	*/
	switch (ftype = filestat(STDIN)) {
	case TTY:
		fprintf(stderr, "FileType(STDIN) = TTY\n");
	break;
	case DISK:
		fprintf(stderr, "FileType(STDIN) = DISK\n");
	break;
	case DIRECTORY:
		fprintf(stderr, "FileType(STDIN) = DIRECTORY\n");
	break;
	case TAPE:
		fprintf(stderr, "FileType(STDIN) = TAPE\n");
	break;
	case PIPE:
		fprintf(stderr, "FileType(STDIN) = PIPE\n");
	break;
	case FIFO:
		fprintf(stderr, "FileType(STDIN) = FIFO\n");
	break;
	case SOCKET:
		fprintf(stderr, "FileType(STDIN) = SOCKET\n");
	break;
	case SYMLINK:
		fprintf(stderr, "FileType(STDIN) = SYMLINK\n");
	break;
	case BADFILETYPE:
		fprintf(stderr, "FileType(STDIN) = BADFILETYPE\n");
	break;
	default:
	    fprintf(stderr, "filestat(stdin) failed: ftype = %d\n", ftype);
	}


	/* Get FileType of stdout	*/
	switch (ftype = filestat(STDOUT)) {
	case TTY:
		fprintf(stderr, "FileType(STDOUT) = TTY\n");
	break;
	case DISK:
		fprintf(stderr, "FileType(STDOUT) = DISK\n");
	break;
	case DIRECTORY:
		fprintf(stderr, "FileType(STDOUT) = DIRECTORY\n");
	break;
	case TAPE:
		fprintf(stderr, "FileType(STDOUT) = TAPE\n");
	break;
	case PIPE:
		fprintf(stderr, "FileType(STDOUT) = PIPE\n");
	break;
	case FIFO:
		fprintf(stderr, "FileType(STDIN) = FIFO\n");
	break;
	case SOCKET:
		fprintf(stderr, "FileType(STDIN) = SOCKET\n");
	break;
	case BADFILETYPE:
		fprintf(stderr, "FileType(STDOUT) = BADFILETYPE\n");
	break;
	default:
	    fprintf(stderr, "filestat(stdout) failed: ftype = %d\n", ftype);
	}

	/* Get FileType of argv[1]	*/
	if (argc == 1) {
		fprintf(stderr, "no filename given\n");
		exit(1);
	}
	if (-1 == (fd = open(argv[1], O_RDONLY))) {
		fprintf(stderr, "can't open %s", argv[1]);
		exit(2);
	}

	switch (ftype = filestat(fd)) {
	case TTY:
		fprintf(stderr, "FileType(fd) = TTY\n");
	break;
	case DISK:
		fprintf(stderr, "FileType(fd) = DISK\n");
	break;
	case DIRECTORY:
		fprintf(stderr, "FileType(fd) = DIRECTORY\n");
	break;
	case TAPE:
		fprintf(stderr, "FileType(fd) = TAPE\n");
	break;
	case PIPE:
		fprintf(stderr, "FileType(fd) = PIPE\n");
	break;
	case FIFO:
		fprintf(stderr, "FileType(STDIN) = FIFO\n");
	break;
	case SOCKET:
		fprintf(stderr, "FileType(STDIN) = SOCKET\n");
	break;
	case BADFILETYPE:
		fprintf(stderr, "FileType(argv[1]) = BADFILETYPE\n");
	break;
	default:
	    fprintf(stderr, "filestat(argv[1]) failed: ftype = %d\n", ftype);
	}


	fprintf(stderr, "Checking printstat ...\n");
	/* Print FileType of stdin	*/
	fprintf(stderr, "FileType(STDIN) = %s\n", printstat(STDIN));


	/* Print FileType of stdout	*/
	fprintf(stderr, "FileType(STDOUT) = %s\n", printstat(STDOUT));


	/* Print FileType of argv[1]	*/
	if (argc == 1) {
		fprintf(stderr, "no filename given\n");
		exit(1);
	}
	if (-1 == (fd = open(argv[1], O_RDONLY))) {
		fprintf(stderr, "can't open %s", argv[1]);
		exit(2);
	}

	fprintf(stderr, "FileType(fd) = %s\n", printstat(fd));

	exit(0);
}
#endif
