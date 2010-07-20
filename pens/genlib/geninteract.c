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
 *  source file:   ./filters/genlib/geninteract.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger Feb 22 1988
 *	Created INT_PAUSE.
 * Joe Dellinger May 20 1989
 *	controltty may be NULL, meaning "/dev/tty" can't be read from.
 *	Handle this case gracefully.
 * Joe Dellinger April 2 1992
 *	Return 0 so dovplot won't get confused and exit accidentally due
 *	to a junk return value happening to be DOVPLOT_EXIT.
 */

/*
 * For most devices, this can be the interact routine.
 * This routine will work if the terminal control during plotting
 * is the same place as standard input.
 */

/*
 * Get a string from the user.
 */

#include <stdio.h>
#include "../include/extern.h"
#include "../include/intcom.h"
#include "../include/err.h"

#include "../utilities/util.h"

void geninteract (int what, FILE *controltty, char *string)
{
    switch (what)
    {
    case INT_PAUSE:
    case INT_F_PAUSE:
    case INT_GET_STRING:
	if (controltty == NULL)
	{
	    ERR (FATAL,name,"Sorry, can't read string from terminal.");
	}
	else
	{
	    if (NULL == fgets (string, 79, controltty)) 
		ERR (FATAL,name,"fgets error");
	}
	break;
    default:
	break;
    }
}
