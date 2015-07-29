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
 *  source file:   ./filters/main_vplot.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stew Levin (SEP), December 7, 1993
 *      Dummy routine MAIN_ for shared library linkage
 * Stew Levin (Mobil) May 8, 1996
 *      Updated signal handling #ifdef's for LINUX
 *
 */

/*
 * generic pen -  VPLOT filter for whatever
 * Keyword: graphics vplot pen
 */

/*
 *  Edit History
 *
 *  "Pen" written 1979 for a PDP-11, Rob Clayton
 *  Eventually modified into "tekpen" by Glenn Kroeger early 1984
 *  Made device independent and extensively modified into "screenpen"
 *  by Joe Dellinger 1984 through 1985
 *  Reworked to be more GKS-like by Glenn and Michel Debiche, 1985
 *  Cleaned up, Joe 1986
 *  Added 'wstype' variable so one program can support multiple
 *		terminal types		-- Chuck Karish, Nov 1986
 *  Raster capability added -- Joe and Steve Cole 1987
 *  Cleaned up and documented for public release, Joe Dellinger 1987
 *  Whew!
 *
 *  Changed system("stty ...") calls to ioctls.  The system() calls
 *  are hanging for jon under 6.0.  Added sigblock() to cleanup
 *  error handler to avoid exit() problems.  Also shut down output
 *  translations to tty graphics terminals.  Stewart A. Levin  6-23-87
 *
 *  Shut down output translations with LLITOUT now that I've gotten
 *  the magic incantation that makes it stick on buggy BSD systems.
 *  Stewart A. Levin  7-5-87
 *
 *  Made "scale" a global so dovplot could use it.
 *  Joe Dellinger Oct 18 1987
 *
 *  Gather up all the incoming plot files, and then have "genreader"
 *  process them all. This allows someone to make a vplot-editor that
 *  can re-order plot files, repeat them, etc, etc. Added "buffer_input"
 *  and "allow_pipe" as part of this.
 *  Joe Dellinger Dec 16 1987
 *
 *  Added "hclose_done" to SEPlib version.
 *  Joe Dellinger Dec 19 1987
 *
 *  Added "reset_parameters" to frontend and dovplot
 *  Joe Dellinger Jan 8 1988
 *
 *  Added group_name, group_number, pltname
 *  Joe Dellinger Jan 20 1988
 *
 *  Inverted window should mean that everything gets clipped.
 *  Make it possible to override standard defaults for SEP.
 *  Joe Dellinger Feb 12 1988
 *
 *  Just use getpar, not getpar_.
 *  Joe Dellinger Feb 16 1988
 *
 *  Split "frontend.c" into 3 files,
 *  main_vplot, init_vplot, and proc_vplot.
 *  This so that other programs can use vplot and still keep their own
 *  mains.
 *  Joe Dellinger Feb 18 1988
 *
 *  Made some SEPlib warning messages a bit clearer.
 *  Joe Dellinger June 30 1988
 *
 * Dave Nichols (SEP), March 25 1989
 * Added cachepipe option to write piped input to a temporary file.
 * This allows them to be reread and replotted under X-window or Sunview.
 *
 * Dave Nichols (SEP), March 27 1989
 * Make sure the temporary file is cleaned up and that seplib input really
 * is coming from a pipe before it is copied.  
 * W. Bauske IBM, 03-27-91
 *	Applied SysV mods for RS/6000 
 *	Removed external re-declare of malloc for RS/6000 
 *
 * Dave Nichols (SEP), March 30 1992
 *  	Use sepxargv and sepxargc for traversing argument list to be consistent.
 *	This can be a problem is Xt initialisation deletes arguments.
 *
 * Dave Nichols (SEP), May 2 1992
 *	Check for environment variable VPLOTSPOOLDIR when creating temporary
 *	files.
 * Dave Nichols (SEP), Dec 11 1992
 *	Added check for zero length files on stdin (e.g. /dev/null )
 * Stew Levin 2-27-95
 *	Solaris mods
 *  Bob Clapp 10-98  Switched to POSIX (ala Sloaris) for LINUX signals
 */

#include <stdlib.h>

/* should be defined in stdlib */
extern int mkstemp (char *tmpl);

#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#define		GETPAR	getpar

#include	<sys/types.h>
#include	<sys/stat.h>
#include	<ctype.h>
#include	<signal.h>

#include <rsf.h>
#include <rsfplot.h>

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

#include "../genlib/genpen.h"
#include "../utilities/util.h"

#include "init_vplot.h"
#include "proc_vplot.h"
#include "main_vplot.h"

static FILE* tempcopy(FILE* infile,char* filename );

extern struct termio tty_clean_state;

extern int      allow_pipe;
extern char     callname[];
extern int      nplots;
extern bool      cachepipe;

/*
 * This routine is responsible for finding the input files,
 * setting up the input and output, and calling init_vplot
 * and proc_vplot. You can link to vplot without using this
 * routine, and so have your own main. See vplothacker.doc to
 * learn how to do this. (Especially you, Jon!)
 */

int main (int argc, char* argv[])
{
    int             in_isatty;
    int len;
    bool            docflag;
    int      infileno=0;
    FILE    *pltinarray[MAXIN];
    char    *pltinname[MAXIN];

    char           *cptr;
    char           *stringptr;
    int             ii;
    FILE           *temp;
    char            string[MAXFLEN + 1];
    
    sf_init(argc,argv);

    if ((stringptr = strrchr(argv[0], '/')))
	strncpy (callname, ++stringptr, 24);
    else
	strncpy (callname, argv[0], 24);

    in_isatty = isatty ((int) (fileno (stdin)));

    if (!sf_getbool ("selfdoc", &docflag)) docflag = (bool) (argc == 1);
    if (in_isatty && docflag)
    {
	for (ii = 0; ii < doclength; ii++)
	    printf ("%s\n", documentation[ii]);
	exit (0);
    }


/*
 ****************************************************************************
 * Set all global variables, open the device.
 ****************************************************************************
 */

    init_vplot (argc,argv);

/*
 ****************************************************************************
 * Start processing input files
 ****************************************************************************
 */

    /*
     * first process pipe input 
     */
    if (!in_isatty)
    {
	if (infileno >= MAXIN)
	{
	    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	}

  	if( dev.cachepipe )
        {
            if( (pltinarray[infileno] = tempcopy( stdin,string) ) == NULL )
            {
                ERR( FATAL, name, "copy of piped input failed");
	    }

            /* check for zero length input (e.g. /dev/null ) */
	    if( pltinarray[infileno] != (FILE*)-1 ) {

		len = strlen(string);
		pltinname[infileno] = sf_charalloc(len+1);

		strncpy( pltinname[infileno],string,len+1);

		infileno++;
	    }
        } else
        {
	    if (!allow_pipe)
	    {
	    	ERR (WARN, name, 
		     "cannot use pipes with this device, try cachepipe=y ");
	    }
	    else
	    {
		pltinname[infileno] = sf_charalloc(6);
	    	strncpy (pltinname[infileno], "stdin", 5);
	    	pltinarray[infileno] = stdin;
	    	infileno++;
	    }
	}
    }

    /*
     * finally process input line for non-getpar arguments and assume they
     * are also input files 
     */
    for (argc--, argv++; argc; argc--, argv++)
    {
	cptr = *argv;
	while (*cptr)
	{
	    if (*cptr == '=')
		break;
	    cptr++;
	}
	if (*cptr)
	    continue;
	cptr = *argv;
	if ((temp = fopen (cptr, "r")) != NULL)
	{
	    if (infileno >= MAXIN)
	    {
		ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	    }

	    len = strlen(cptr);
	    pltinname[infileno] = sf_charalloc(len+1);
	    
	    strncpy (pltinname[infileno], cptr,len+1);
	    pltinarray[infileno] = temp;
	    infileno++;
	}
	else
	{
	    ERR (WARN, name, "cannot open %s", cptr);
	}
    }

/*
 ****************************************************************************
 * Go do the plots
 ****************************************************************************
 */

    proc_vplot (infileno, pltinarray, pltinname);

    /*  delete the temporary copy of piped input if there is one*/
    removtemp();

    /* All done */
    exit (0);
}

/* routine to copy a file to a temporary file, used to copy stdin so that
 * it can be reread as required
 */

#define XFER_SIZE 1024

/* the name of the temporary file */
static char *tspoolnm=(char*)NULL;

void removtemp(void)
/*< remove temporary file >*/
{
    if( tspoolnm != (char*)NULL )
        if( unlink(tspoolnm) )
               ERR (WARN, name, "unable to delete temporary file");
    tspoolnm = (char*)NULL;
}

static FILE* tempcopy(FILE* infile,char* filename )
/* temp copy */
{
    FILE *temp;
    size_t len,total;
    char xfer_buf[ XFER_SIZE];

    temp = sf_tempfile(&tspoolnm,"w+");
 
    /* now copy everything from the input stream to our temporary file */
    total=0;
    while( (len=fread( xfer_buf, sizeof(char), XFER_SIZE, infile ) ) != 0 ){
	if( ferror(infile) ) {
	    ERR (WARN, name, "read error from pipe ");
	    return NULL;
	}
	    
        if( ( fwrite(xfer_buf,sizeof(char),len,temp) != len ) || ferror(temp) ) 	{
	    ERR (WARN, name, "write error to temp ");
	    return NULL;
	}
	total += len;
    }

    fclose(temp);

    /* We could unlink the file after we reopen it to make it really
     * temporary but we don't know for sure that a filter isn't going
     * to close the file and try to reopen it.
     */

    /* check for zero length input */
    if( total == 0 ) return (FILE*)-1;

    strcpy( filename, tspoolnm );
    return fopen( tspoolnm, "r");
}

