/* Handling warning and error messages. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#include "error.h"
#include "getpar.h"

void sf_error( char *format, ... )
/*< Outputs an error message to stderr and terminates the program. 

Format and variable arguments follow printf convention. Additionally, a ':' at
the end of format adds system information for system errors. >*/
{
    va_list args;

    (void) fflush(stdout);
    va_start(args,format);

    /* print out name of program causing error */
    fprintf (stderr,"%s: ",sf_getprog());

    /* print out remainder of message */
    (void) vfprintf( stderr, format, args );
    va_end(args);

    /* if format ends with ':', print system information */
    if (format[0] != '\0' && format[strlen(format)-1] == ':')
	fprintf (stderr, " %s", strerror(errno));

    fprintf (stderr, "\n");

    exit(EXIT_FAILURE); 
}

void sf_warning( char *format, ... )
/*< Outputs a warning message to stderr. 

Format and variable arguments follow printf convention. Additionally, a ':' at
the end of format adds system information for system errors. >*/
{
    va_list args;
    
    (void) fflush(stdout);
    va_start(args,format);

    /* print out name of program causing error */
    fprintf(stderr,"%s: ",sf_getprog()); 

    /* print out remainder of message */
    (void) vfprintf( stderr, format, args );
    va_end(args);

    /* if format ends with ':', print system information */
    if (format[0] != '\0' && format[strlen(format)-1] == ':')
	fprintf (stderr, " %s", strerror(errno));

    fprintf (stderr, "\n");
}

/* 	$Id$	 */

