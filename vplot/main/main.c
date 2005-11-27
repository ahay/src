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

#include<sitedef.h>
#include <stdlib.h>

#ifdef SEP
extern void     sepwhere ();
extern char     sepoutwhere[];
extern char     sepheadwhere[];

#define		OUT	sepoutwhere
#define		HEAD	sepheadwhere
#if !defined(RS6000) && !defined(SOURCE)
#define		SOURCE  "\014Joe Dellinger, Stanford Exploration Project\014"
#endif
#include	<string.h>
#include	<sep.main>
#define		GETPAR	fetch

#else /* SEP */
#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#define		GETPAR	getpar
#endif /* SEP */


#if defined(HAVE_TERMIO_H)
#include	<termio.h>
#else
#if defined (HAVE_SGTTY_H)
#include	<sgtty.h>
#else
#include	<sys/ioctl.h>
#include	<sgtty.h>
#endif
#endif /* USG */
#include	<sys/types.h>
#include	<sys/stat.h>
#include	<ctype.h>
#include	<signal.h>

#include	<vplot.h>

#include	"./include/params.h"	/* for machine dependencies */
#include	"./include/enum.h"
#include	"./include/err.h"
#include	"./include/attrcom.h"
#include	"./include/intcom.h"
#include	"./include/mesgcom.h"
#include	"./include/erasecom.h"
#include	"./include/closestat.h"
#include	"./include/pat.h"
#include	"./include/vertex.h"
#include	"./include/round.h"
#include	"./include/extern.h"


#if defined(HAVE_TERMIO_H)
#else /* USG */
/*
 * signal catching
 */
#ifdef SIGFNC_RTN_VOID
void            cleanup ();
#else
int             cleanup ();
#endif
int             signum[] =
{
#ifdef LINUX
 SIGHUP, SIGINT, SIGQUIT, SIGIOT, SIGBUS, SIGPIPE, SIGTERM, SIGXCPU, SIGXFSZ
#else
 SIGHUP, SIGINT, SIGQUIT, SIGIOT, SIGEMT, SIGPIPE, SIGTERM, SIGXCPU, SIGXFSZ
#endif
};
#define NOSIG (sizeof (signum)/sizeof (int))	/* number of signals caught */
#if defined(SOLARIS ) || defined(LINUX)
struct sigaction   errhandler =
{
#ifdef LINUX
 cleanup, 0, 0
#else
 0, cleanup, 0
#endif
};
struct sigaction   ignored =
{
#ifdef LINUX
 SIG_IGN, 0, 0
#else
 0, SIG_IGN, 0
#endif
};
struct sigaction   oldvec;
#else /*SOLARIS*/
int             sigvec ();
struct sigvec   errhandler =
{
 cleanup, 0, 0
};
struct sigvec   ignored =
{
 SIG_IGN, 0, 0
};
struct sigvec   oldvec;
#endif /*SOLARIS*/
#endif /* USG */

#if defined(HAVE_TERMIO_H)
extern struct termio tty_clean_state;
#else /* USG */
extern struct sgttyb tty_clean_state;	/* external for utilities */
extern int      tty_clean_local_mode;
#endif /* USG */
extern int      allow_pipe;
extern char     callname[];
extern int      nplots;
extern int      allowecho;
extern int      cachepipe;

/*
 * file and terminal control variables
 */
extern int      genmessage ();
extern int      (*message) ();
extern FILE    *pltout;

FILE           *fopen ();
FILE           *fdopen ();
FILE	       *tempcopy();
extern int	unlink();

extern FILE    *pltinarray[MAXIN];
extern char     pltinname[MAXIN][MAXFLEN + 1];
extern int      infileno;
extern int      pltoutfd;

/*
 * This routine is responsible for finding the input files,
 * setting up the input and output, and calling init_vplot
 * and proc_vplot. You can link to vplot without using this
 * routine, and so have your own main. See vplothacker.doc to
 * learn how to do this. (Especially you, Jon!)
 */

#ifdef SEP
int             xsepxargc;
char          **xsepxargv;
char            scrap[MAXFLEN + 1];
int             hclose_done = NO;
int             fake_header = NO;
MAIN ()

#else /* SEP */
int             sepxargc;
char          **sepxargv;	/* for getpar */
MAIN_(){} /* dummy for shared linkage */
main (argc, argv)
    int             argc;
    char           *argv[];
#endif /* SEP */
{
#ifndef SEP
int             in_isatty, num_vplot, docflag;
char            instring[MAXFLEN + 1];
#endif /* SEP */

char           *cptr;
char           *stringptr;
int             ii;
FILE           *temp;
char            string[MAXFLEN + 1];
int		tempfileindex = -1;

    nulldev ();			/* Just to make sure it gets loaded */

#ifndef SEP
    if (stringptr = strrchr(argv[0], '/'))
	strncpy (callname, ++stringptr, 24);
    else
	strncpy (callname, argv[0], 24);
#else /* SEP */
    if (stringptr = strrchr(sepxargv[0], '/'))
	strncpy (callname, ++stringptr, 24);
    else
	strncpy (callname, sepxargv[0], 24);
#endif /* SEP */

#ifdef SEP
    pltout = outstream;
    if (redout ())
    {
	getch ("head", "s", scrap);
	if (strcmp (scrap, "/dev/null") == 0 &&
	    strcmp (sepheadwhere, "/dev/null") == 0)
	{
	    fake_header = YES;
	    getch ("out", "s", scrap);
	    if (strcmp (scrap, "stdout") != 0)
	    {
		/*
		 * They probably want the header output into the redirected
		 * output. (SEP only) 
		 */
		headstream = stdout;
		headfd = fileno (headstream);
		Puthead ("Fake header for special device-dependent data only.\n");
		Puthead ("To get input history passed along, over-ride default with head = stdout\n");
		/*
		 * If the output is going into a file, then put this
		 * information into the header file. 
		 */
		if (strcmp (scrap, "/dev/tty") != 0)
		{
		    fullnm (scrap, MAXFLEN + 1);
		    Puthead ("\tin=%s\n", scrap);
		}
		else
		{
		    Puthead ("\tin= nowhere\n");
		    Puthead ("\t(Sorry, the data went to the terminal; it's gone now!)\n");
		}

		/*
		 * Have to use my own puthead routine. Standard ones such as
		 * putch, etc, don't work with this. They remember where the
		 * header used to be. Puthead does this, checking for
		 * headstream not null. This is so this will work without
		 * sep, as well. 
		 */
	    }
	}
    }

    if (!fake_header)
    {
	Puthead ("\tn3 = unknown\n\t(Sorry, SEPlib requires the header be closed now.)\n");
	Puthead ("\t(Normally you won't care that n3 is unknown unless you're using Raspen!\n)");
	hclose ();
	hclose_done = YES;
    }
#else /* SEP */

    /*
     * If no arguments, and not in a pipeline, self document "wstype="
     * doesn't count as an argument for our purposes 
     */
    in_isatty = isatty ((int) (fileno (stdin)));
    sepxargc = argc;
    sepxargv = argv;
    docflag = 0;
    if (argc == 1)
	docflag = 1;
    if ((argc == 2) && !strncmp ("wstype=", argv[1], 7))
	docflag = 1;
    getpar ("selfdoc", "1", &docflag);
    if (in_isatty && docflag)
    {
	for (ii = 0; ii < doclength; ii++)
	    printf ("%s\n", documentation[ii]);
	exit (0);
    }

    pltout = stdout;
#endif /* SEP */

#if defined(HAVE_TERMIO_H)

#else /* USG */
    /*
     * This getpar for signal is only included for debugging purposes. By
     * using a signal option, one can stop any signals from being caught. 
     */
    if (getpar ("signal", "s", string) == 0)
    {
/*#ifdef SOLARIS*/
#if defined(SOLARIS) || defined(LINUX)
        sigfillset(&(errhandler.sa_mask));
#endif
	for (ii = 0; ii < NOSIG; ++ii)
	{
#if defined(SOLARIS) || defined(LINUX)
	    if (-1 == sigaction (signum[ii], &ignored, &oldvec))
	    {
		ERR (FATAL, name, "Bad sigvec call!");
	    }
	    if (oldvec.sa_handler == ignored.sa_handler)
		(void) sigaction (signum[ii], &oldvec, (struct sigaction *) NULL);
	    else
		(void) sigaction (signum[ii], &errhandler, (struct sigaction *) NULL);
#else
	    if (-1 == sigvec (signum[ii], &ignored, &oldvec))
	    {
		ERR (FATAL, name, "Bad sigvec call!");
	    }
	    if (oldvec.sv_handler == ignored.sv_handler)
		(void) sigvec (signum[ii], &oldvec, (struct sigvec *) NULL);
	    else
		(void) sigvec (signum[ii], &errhandler, (struct sigvec *) NULL);
#endif
	}
    }
#endif /* USG */

/*
 ****************************************************************************
 * Set all global variables, open the device.
 ****************************************************************************
 */

    init_vplot ();

/*
 ****************************************************************************
 * Start processing input files
 ****************************************************************************
 */


#ifdef SEP
    if (instream != NULL)
    {
	if (infileno >= MAXIN)
	{
	    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	}
	if( cachepipe && isapipe(fileno( instream)) )
	{
	    if( (pltinarray[infileno] = tempcopy(instream,string) ) == NULL )
	    {
	    	ERR( FATAL, name, "copy of piped input failed");
	    }
	 
            /* check for zero length input (e.g. /dev/null ) */
	    if( pltinarray[infileno] != (FILE*) -1 ) {

	    strcpy( pltinname[infileno], string );
	    /*remember what number this file is so we can delete it later*/
	    tempfileindex = infileno;
	    infileno++;
	    }
	}
	else
	{
	    if (!allow_pipe && isapipe (fileno (instream)))
	    {
	    	ERR (WARN, name, "cannot use pipes with this device, try cachepipe=y ");
	    }
	    else
	    {
	        strcpy (pltinname[infileno], "Pipe");
	    	pltinarray[infileno] = instream;
	    	infileno++;
	    }
	}
    }
    else
	ERR (WARN, name, "cannot read input pipe");

    xsepxargc = sepxargc;
    xsepxargv = sepxargv;

    for (xsepxargc--, xsepxargv++; xsepxargc; xsepxargc--, xsepxargv++)
    {
	cptr = *xsepxargv;
	while (*cptr)
	{
	    if (*cptr == '=')
		break;
	    cptr++;
	}
	if (*cptr)
	    continue;
	/* Ignore dummy arguments */
	if (strcmp (*xsepxargv, "dummy") == 0)
	    continue;
	if ((temp = fopen (*xsepxargv, "r")) == NULL)
	{
	    ERR (WARN, name, "cannot open header file %s", *xsepxargv);
	    continue;
	}
	fclose (temp);
	if (getch2 ("in", "s", string, *xsepxargv))
	{
	    if ((temp = fopen (string, "r")) != NULL)
	    {
		Puthead ("   +  %s --> in = %s\n", *xsepxargv, string);
		if (infileno >= MAXIN)
		{
		    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
		}
		strcpy (pltinname[infileno], string);
		pltinarray[infileno] = temp;
		infileno++;
	    }
	    else
	    {
		ERR (WARN, name, "cannot open input file %s", string);
	    }
	}
    }

#else /* SEP */
    /*
     * first process pipe input 
     */
    if (!in_isatty)
    {
	if (infileno >= MAXIN)
	{
	    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	}

  	if( cachepipe )
        {
            if( (pltinarray[infileno] = tempcopy( stdin,string) ) == NULL )
            {
                ERR( FATAL, name, "copy of piped input failed");
	    }

            /* check for zero length input (e.g. /dev/null ) */
	    if( pltinarray[infileno] != (FILE*)-1 ) {

	    strcpy( pltinname[infileno], string );
	    /* remember what number this file is so we can delete it later*/
	    tempfileindex = infileno;
            infileno++;
	    }
        }
        else
        {
	    if (!allow_pipe)
	    {
	    	ERR (WARN, name, "cannot use pipes with this device, try cachepipe=y ");
	    }
	    else
	    {
	    	strcpy (pltinname[infileno], "stdin");
	    	pltinarray[infileno] = stdin;
	    	infileno++;
	    }
	}
    }

    /*
     * next process in= inputfiles If they set num_vplot, also look for in1=
     * in2= etc 
     */

    num_vplot = 0;
    getpar ("numvplot", "d", &num_vplot);

    for (ii = 0; ii <= num_vplot; ii++)
    {
	if (ii == 0)
	    strcpy (instring, "in");
	else
	    sprintf (instring, "in%d", ii);

	if (getpar (instring, "s", string))
	{
	    if ((temp = fopen (string, "r")) != NULL)
	    {
		if (infileno >= MAXIN)
		{
		    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
		}
		strcpy (pltinname[infileno], string);
		pltinarray[infileno] = temp;
		infileno++;
	    }
	    else
	    {
		ERR (WARN, name, "cannot open %s", string);
	    }
	}
    }

    /*
     * finally process input line for non-getpar arguments and assume they
     * are also input files 
     */
    for (sepxargc--, sepxargv++; sepxargc; sepxargc--, sepxargv++)
    {
	cptr = *sepxargv;
	while (*cptr)
	{
	    if (*cptr == '=')
		break;
	    cptr++;
	}
	if (*cptr)
	    continue;
	cptr = *sepxargv;
	if ((temp = fopen (cptr, "r")) != NULL)
	{
	    if (infileno >= MAXIN)
	    {
		ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	    }
	    strcpy (pltinname[infileno], cptr);
	    pltinarray[infileno] = temp;
	    infileno++;
	}
	else
	{
	    ERR (WARN, name, "cannot open %s", cptr);
	}
    }
#endif /* SEP */

/*
 ****************************************************************************
 * Go do the plots
 ****************************************************************************
 */

    proc_vplot ();

#ifdef SEP
    if (!hclose_done)
    {
	Puthead ("\tn3=%d\n", nplots);
	hclose ();
	hclose_done = YES;
    }
#endif /* SEP */

    /*  delete the temporary copy of piped input if there is one*/
    removtemp();

    /* All done */
    exit (0);
    return 0;
}

#if defined(HAVE_TERMIO_H)
#else /* USG */
#ifdef SIGFNC_RTN_VOID
void cleanup ()
#else
cleanup ()
#endif
{
#ifndef SOLARIS
#ifndef LINUX
    sigblock (~(SIGKILL | SIGSTOP | SIGCONT));
#endif
#endif
    dev.close (CLOSE_INTERRUPT);
    message (MESG_ON);
    ERR (COMMENT, name, "Interrupted out.");
    dev.close (CLOSE_DONE);
    /*  delete the temporary copy of piped input if there is one*/
    removtemp();
    /*
     * Let them see what they are doing again 
     */
    if (!allowecho)
    {
/*#ifdef SOLARIS*/
#if defined(SOLARIS) || defined(LINUX)
	ioctl (pltoutfd, TCSETAW, (char *) (&tty_clean_state));
#else
	ioctl (pltoutfd, TIOCLSET, (char *) (&tty_clean_local_mode));
	ioctl (pltoutfd, TIOCSETN, (char *) (&tty_clean_state));
#endif
    }
    exit (0);
}
#endif /* USG */

/* routine to copy a file to a temporary file, used to copy stdin so that
 * it can be reread as required
 */

#define XFER_SIZE 1024

/* the name of the temporary file */
static char *tspoolnm=(char*)NULL;

removtemp(){
    if( tspoolnm != (char*)NULL )
        if( unlink(tspoolnm) )
               ERR (WARN, name, "unable to delete temporary file");
    tspoolnm = (char*)NULL;
}

FILE* tempcopy( infile, filename )
FILE* infile;
char* filename;
{
    FILE *temp;
    int len,total;
    char xfer_buf[ XFER_SIZE];
    char* spooldirnm;


    /* make up a temporary file name 
     * I know there are beter routines to make up the name 
     * but this one is more portable.
     * PEN_SPOOL is a directory we can put temporary files in
     * it is defined in params.h. It can be overridden by the
     * environment variable VPLOTSPOOLDIR.
     */
    if( (spooldirnm = getenv("VPLOTSPOOLDIR")) != NULL ){
	tspoolnm = (char *) malloc( strlen( spooldirnm ) +13 );
        strcpy( tspoolnm, spooldirnm );
    }else{
        tspoolnm = (char *) malloc( strlen( PEN_SPOOL ) + 13 );
        strcpy( tspoolnm, PEN_SPOOL );
    }
    tspoolnm = strcat( tspoolnm, "/vplotXXXXXX" );
    tspoolnm = mktemp( tspoolnm );

    if ( (temp = fopen(tspoolnm, "w")) == NULL) {
	ERR (WARN, name, "unable to create temporary file");
	return NULL;
    }

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
     * Hence the buisness with tempfileindex in the main routine.
     */

    /* check for zero length input */
    if( total == 0 ) return (FILE*)-1;

    strcpy( filename, tspoolnm );
    return fopen( tspoolnm, "r");
}

