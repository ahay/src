/*
  Copyright (C) 1987 The Board of Trustees of Stanford University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

#include <termios.h>
#include <sys/ioctl.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../include/prototypes.h"
#include "../include/err.h"
#include "../include/closestat.h"
#include "../include/mesgcom.h"
#include "../include/extern.h"
extern int      device_open;

extern struct termio tty_clean_state;

extern void removtemp(void);

#include <stdarg.h>

int ERR (int type, char *filter, const char *fmt,...)
/*<  error reporting for VPLOT filters
 * To print text to the outside world, other routines must go through err,
 * which in turn must go through message. >*/
{
va_list         apdum;

char            string[200];
char            text[150];

    va_start (apdum,fmt);

    (void) vsprintf (text, fmt, apdum);
    (void) sprintf (string, "%s: ", filter);

    va_end (apdum);

    message (MESG_READY, NULL);

    switch (type)
    {
    case WARN:
	(void) strcat (string, "(warning) ");
	(void) strcat (string, text);
	(void) strcat (string, CRLF);
	message (MESG_TEXT, string);
	message (MESG_DONE, NULL);
	break;
    case COMMENT:
	(void) strcat (string, " ");
	(void) strcat (string, text);
	(void) strcat (string, CRLF);
	message (MESG_TEXT, string);
	message (MESG_DONE, NULL);
	break;
    case FATAL:
    default:
	(void) strcat (string, "(fatal) ");
	(void) strcat (string, text);
	(void) strcat (string, CRLF);
	message (MESG_TEXT, string);
	message (MESG_DONE, NULL);
	removtemp();
	if (device_open)
	{
	    dev.close (CLOSE_ERROR);
	    message (MESG_ON, NULL);
	    dev.close (CLOSE_DONE);
	    fflush (stdout);
	}
	exit (-1);
    }

    return 0;
}
