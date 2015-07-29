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
 *  source file:   ./filters/proc_vplot.c
 *
 * Joe Dellinger (SEP), Feb 19 1988
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 *
 * Joe Dellinger Feb 22 1988
 *	Created INT_PAUSE to be separate from INT_GET_STRING.
 * Joe Dellinger Feb 27 1988
 *	Interact option turns of endpausing.
 * W. Bauske IBM 03-26-91 
 *	Apply SysV fixes for RS/6000
 * Stew Levin MOBIL 2-27-95
 *	Solaris mods
 & Bob Clapp
 *  Changed signals from BSD to POSIX1 for LINUX
 */

#include	<stdio.h>
#include	<math.h>
#define		GETPAR	getpar

#include	<sys/types.h>
#include	<sys/stat.h>
#include        <sys/ioctl.h>

#include	<ctype.h>
#include	<string.h>

#include	<rsfplot.h>

#include	"../include/params.h"	/* for machine dependencies */
#include	"../include/enum.h"
#include	"../include/err.h"
#include	"../include/attrcom.h"
#include	"../include/intcom.h"
#include	"../include/mesgcom.h"
#include	"../include/erasecom.h"
#include	"../include/closestat.h"
#include	"../include/pat.h"
#include	"../include/vertex.h"
#include	"../include/round.h"
#include	"../include/extern.h"

#include "../utilities/util.h"

#include "dovplot.h"

extern struct termio tty_clean_state;

extern bool      need_end_erase;
extern int      ever_called;
extern int      out_isatty;
extern int      nplots;
extern bool      endpause;
extern int      epause;
extern char     interact[];
extern FILE    *pltin;
extern FILE    *controltty;
extern char     pltname[];

void proc_vplot (int infileno, FILE *pltinarray[], char *pltinname[])
/*< This routine is responsible for processing the input files,
 * and performing the necessary pausing, etc, that may be needed
 * at the end before exiting. >*/
{
    char            *string=NULL;
    
/*
 * Finally, shove all the plot files off to be done!
 */
    
    dev.reader(infileno, pltinarray, pltinname);

/*
 * Normally, *genreader will be gen_do_dovplot, found in genlib
 */

    if (ever_called)
    {
	dev.close (CLOSE_FLUSH);
	if (epause > 0)
	{
	    sleep ((unsigned) epause);
	}
	if (dev.need_end_erase)
	{
	    dev.erase (ERASE_END);
	}
	/*
	 * Inquire point back from device. Skip endpause stuff if we do
	 * interact, Since that's a pause in itself. 
	 */
	if (interact[0] != '\0')
	{
	    getapoint ();
	}
	else
	{

/*
 * Pause at the end if the user specifically asks us to on the command line,
 * even if we don't think we should because we think it's a file.
 */
	    if (epause <= 0 &&
		(out_isatty || (NULL != (string = sf_getstring ("endpause"))))
		 && endpause)
	    {
		dev.close (CLOSE_PAUSE);
		dev.interact (INT_F_PAUSE, controltty, string);
	    }
	}
	message (MESG_ON,NULL);
	dev.close (CLOSE_NORMAL);
	dev.close (CLOSE_DONE);
	nplots++;
    }
    else
    {
	dev.close (CLOSE_NOTHING);
	ERR (COMMENT, name, "No input?");
	dev.close (CLOSE_DONE);
    }

}
