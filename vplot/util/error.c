/*
 * Copyright 1987 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
 */

/*
 *
 *  source file:   ./filters/utilities/error.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (SEP), June 23 1987
 *	Changed system("stty echo") to ioctl() to restore original tty
 *	settings after we've deliberately forced echo off.
 * Stewart A. Levin (SEP), July 5,1987
 *      Added reset of local mode word corresponding to frontend.c change.
 * Joe Dellinger, Nov 9, 1987
 *	Stew's changes make it necessary to explicitly put "Carriage-Return
 *	Line-Feed" instead of just "\n", since output stream may be
 *	uncooked. "CRLF" is defined in mesgcom.h.
 * Joe Dellinger, Jan 10, 1988
 *	Up to 8 arguments instead of just 3! This routine isn't very kosher.
 * Joe Dellinger, Feb 16, 1988
 *	Use vsprintf on machines that have them.
 * W. Bauske IBM 03-26-91
 *	Apply SysV fixes for RS/6000
 *	removed re-declare of sprintf and strcat for RS/6000
 &  Bob Clapp 10-98
 *   Changed LINUX to Posix1 from BSD
 */

/*
 * error reporting for VPLOT filters
 * To print text to the outside world, other routines must go through err,
 * which in turn must go through message.
 */

#include<sitedef.h>
#if defined (HAVE_TERMIO_H)
#include <termio.h>
#else /* USG */
#if defined(LINUX)
#include <bsd/sgtty.h>
#else
#include <sys/ioctl.h>
#include <sgtty.h>
#endif
#endif /* USG */

#include <stdio.h>
#include "../include/prototypes.h"
#include "../include/err.h"
#include "../include/closestat.h"
#include "../include/mesgcom.h"
#include "../include/extern.h"
extern int      device_open;

#if defined(HAVE_TERMIO_H)
extern struct termio tty_clean_state;
#else /* USG */
extern struct sgttyb tty_clean_state;
extern int      tty_clean_local_mode;
#endif /* USG */


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
#include <stdarg.h>
ERR (int type, char *filter, char *fmt,...)
{
va_list         apdum;

char            string[200];
char            text[150];

/*#if defined(HAVE_STRING_H)*/
/*#include <string.h>*/
/*#else*/
#include <strings.h>
/*#endif*/

    va_start (apdum,fmt);
/*
    type = va_arg (apdum, int);
    filter = va_arg (apdum, char *);
    fmt = va_arg (apdum, char *);
*/

    (void) vsprintf (text, fmt, apdum);
    (void) sprintf (string, "%s: ", filter);

    va_end (apdum);

_XFUNCPROTOEND
#else

ERR (type, filter, fmt, a1, a2, a3, a4, a5, a6, a7, a8)
    int             type;
    char           *filter, *fmt;
    double          a1, a2, a3, a4, a5, a6, a7, a8;
{
char            string[200];
char            text[150];

#if defined(__stdc__) || defined(__STDC__)
#include <string.h>
#else
#include <strings.h>
#endif

    (void) sprintf (text, fmt, a1, a2, a3, a4, a5, a6, a7, a8);
    (void) sprintf (string, "%s: ", filter);

#endif

    message (MESG_READY);

    switch (type)
    {
    case WARN:
	(void) strcat (string, "(warning) ");
	(void) strcat (string, text);
	(void) strcat (string, CRLF);
	message (MESG_TEXT, string);
	message (MESG_DONE);
	break;
    case COMMENT:
	(void) strcat (string, " ");
	(void) strcat (string, text);
	(void) strcat (string, CRLF);
	message (MESG_TEXT, string);
	message (MESG_DONE);
	break;
    case FATAL:
    default:
	(void) strcat (string, "(fatal) ");
	(void) strcat (string, text);
	(void) strcat (string, CRLF);
	message (MESG_TEXT, string);
	message (MESG_DONE);
	removtemp();
	if (device_open)
	{
	    dev.close (CLOSE_ERROR);
	    message (MESG_ON);
	    dev.close (CLOSE_DONE);
	    fflush (pltout);
	}
	if (!allowecho)		/* restore terminal to original tty state */
	{
#if defined(HAVE_TERMIO_H)
	    ioctl ((int) (fileno (pltout)), TCSETAW, &tty_clean_state);
#else /* USG */
	    (void) ioctl ((int) (fileno (pltout)), TIOCLSET,
			  (char *) (&tty_clean_local_mode));
	    (void) ioctl ((int) (fileno (pltout)), TIOCSETN,
			  (char *) (&tty_clean_state));
#endif /* USG */
	}
	exit (-1);
    }
}
