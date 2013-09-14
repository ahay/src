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
#include <setjmp.h>

#include "_bool.h"
#include "error.h"

jmp_buf *python_except=NULL;

/* provided by getpar */
char* sf_getprog (void); 
bool sf_getbool (const char* key,/*@out@*/ bool* par);

void sf_error( const char *format, ... )
/*< Outputs an error message to stderr and terminates the program. 
---
Format and variable arguments follow printf convention. Additionally, a ':' at
the end of format adds system information for system errors. >*/
{
    va_list args;
    char *prog;
    int last;

    (void) fflush(stdout);
    va_start(args,format);

    last = strlen(format)-1;

    /* print out name of program causing error */
    prog = sf_getprog();
    fprintf (stderr,"%s: ",prog);

    /* print out remainder of message */
    (void) vfprintf( stderr, format, args );
    va_end(args);

    /* if format ends with ':', print system information */
    if (format[0] != '\0' && format[last] == ':')
	fprintf (stderr, " %s", strerror(errno));

    fprintf (stderr, "\n");
    (void) fflush(stderr);

    if (0==strcmp("python",prog)) longjmp(*python_except,1);
 
    exit(EXIT_FAILURE);
}

void sf_warning( const char *format, ... )
/*< Outputs a warning message to stderr. 
---
Format and variable arguments follow printf convention. Additionally, a ':' at
the end of format adds system information for system errors. >*/
{
    va_list args;
    char *prog;
    int last;
    bool dryrun;

    if (!sf_getbool("--dryrun",&dryrun)) dryrun=false;
    if (dryrun) return; /* silent if dryrun */
    
    (void) fflush(stdout);
    va_start(args,format);

    last = strlen(format)-1;

    /* if format ends with ';', carriage return */
    if (format[0] != '\0' && format[last] == ';')
	fprintf (stderr, "\r");

    /* if format starts with '.', new line */
    if (format[0] == '.') {
	fprintf (stderr, "\n");
	if (0 == last) {
	    (void) fflush(stderr);
	    return;
	}
    }

    /* print out name of program causing error */
    prog = sf_getprog();
    fprintf(stderr,"%s: ",prog); 

    /* print out remainder of message */
    (void) vfprintf( stderr, format, args );
    va_end(args);

    /* if format ends with ':', print system information */
    if (format[0] != '\0' && format[last] == ':')
	fprintf (stderr, " %s", strerror(errno));

    /* if format ends with ';', do not end line */
    if (format[0] == '\0' || format[last] != ';')
	fprintf (stderr, "\n");

    (void) fflush(stderr);
}

/* 	$Id$	 */

