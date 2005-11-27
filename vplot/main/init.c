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
 *  source file:   ./filters/init_vplot.c
 *
 * Joe Dellinger (SEP), Feb 18 1988
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger, Feb 20 1988
 *	#define GETPAR to be fetch if SEPlib version.
 * Joe Dellinger Feb 24 1988
 * 	Moved gen_do_dovplot to genlib, where it belongs.
 *	txfont, txprec, txovly, overlay defaults remembered.
 * Joe Dellinger Feb 28 1988
 *	Changed fatness to scale with size of screen and scale,
 *	like any regular geometric attribute.
 * Joe Dellinger Mar 4 1988
 *	Fixed slight bug with fatness scaling and style=old.
 * Joe Dellinger Mar 17 1988
 *	Added "colormask" option.
 * Joe Dellinger April 15 1988
 *	Check for illegal colormask values.
 * Joe Dellinger Oct 28 1988
 *	Put in "no_stretch_text" flag.
 * Steve Cole Feb 2 1989
 *      Added default x and y scale parameters for use by vppen.
 * Joe Dellinger May 7 1989
 *	Added break=ignore option.
 *	Made txsquare=yes the default.
 * Joe Dellinger May 20 1989
 *	"/dev/tty" may NOT be openable. Behave gracefully in this case.
 * Joe Dellinger August 2 1989
 *	Put in redpow, greenpow, bluepow, redmap, greenmap, bluemap
 * W. Bauske IBM 03-26-91 
 *	Apply the SysV fixes to the code for the RS/6000
 *	Removed re-declare of system alloc routines for RS/6000
 * Dave Nichols 9-17-92
 *      Include stdlib.h if it is available instead of defining malloc.  
 * Stew Levin  2/27/95 
 *	Solaris mods.
 * Bob Clapp 1/24/96
 *      Changed a whole bunch of fetches to getches
 * Bob Clapp 10/98 
 *      Changed signals from bsd to posix1 for linux
 */
#include<sitedef.h>
#include	<stdio.h>
#include	<math.h>
#if defined(HAVE_TERMIO_H)
#include	<termio.h>
#endif
#if defined(SYS_IOCTL_H)
#include	<sys/ioctl.h>
#endif
#if defined(HAVE_SGTTY_H)
#include	<sgtty.h>
#endif
#if  defined(HAVE_SYS_TYPES_H)
#include	<sys/types.h>
#endif

#if  defined(HAVE_SYS_STAT_H)
#include	<sys/stat.h>
#endif
#include	<ctype.h>
#include	<string.h>

#include	"./include/vplot.h"

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

#define		OPEN_ERROR	-1

#ifdef SEP
#define		GETPAR	fetch
#else
#define		GETPAR	getpar
#endif /* SEP */

#if defined(HAVE_TERMIO_H)
struct termio	tty_clean_state;
struct termio	tty_plot_state;
#else
struct sgttyb   tty_clean_state;
struct sgttyb   tty_plot_state;
int             tty_clean_local_mode;
int             tty_plot_local_mode;
#endif /* USG */

/* 
 * The following variables must ALWAYS
 * be set in the device open or device reset file.
 * All their defaults are set to absurd values so that
 * we'll crash quickly it they're not set as required.
 */
/* Screen dimensions */
int             dev_xmax = 0, dev_ymax = 0, dev_xmin = 0, dev_ymin = 0;
/* Number of pixels per inch in the horizontal (X) direction on the device */
float           pixels_per_inch = 0.;
/* vertical height (on screen) / horizontal width for a single pixel */
float           aspect_ratio = 0.;
/*
 * Number of SETTABLE COLORS that the device has.
 * Non settable colors don't count.
 */
int             num_col = -1;

/*
 * Other things that may need to be reset in dev.open
 * (Can't reset them some of them in dev.reset, because the user may
 * override from the command line.)
 */
/* Does device need erase at end? */
int             need_end_erase = NO;
/* should the output be buffered? */
int             buffer_output = YES;
/* should the input be buffered? */
int             buffer_input = YES;
/* should pipes be allowed for input? */
int             allow_pipe = YES;
/* Can the device do its own clipping? (of vectors and polygons.) */
int             smart_clip = NO;
/* Can the device stretch AND clip its own raster? */
int             smart_raster = NO;

/*
 * These may be reset to device-dependent defaults.
 */
float           fatmult = 1.;
float           patternmult = 1.;

/*
 * Look in extern.h for a better-organized list of the preceding.
 * Here I have left things which are also filter options with the
 * other filter options, even though they also belong in the above
 * categories.
 */

/*
 * flags and variables
 */
char            wstype[25];
char            callname[25];
int             xcenterflag, ycenterflag;
int             ever_called = NO;
int             first_time = YES;
int             device_open = NO;	/* used by ERR */
int             out_isatty = YES;	/* If NO, output is a file or pipe */
int             nplots = 0;	/* number of plots made */
int             no_stretch_text = YES;	/* Don't allow stretched text? */

/*
 * coordinate variables
 */
int             xwmax, xwmin, ywmax, ywmin;	/* window */
int             xnew, ynew;	/* new pen location */
int             xold, yold;	/* old pen location */
int             xorigin = 0, yorigin = 0;	/* global "origin" */
char           *txbuffer;
int             txbuflen, vxbuflen;
struct vertex  *vxbuffer;
int             xret, yret;
float           mxx, mxy, myx, myy;
/*
 * This is so the device can throw in a coordinate transformation
 * for its convenience ON TOP OF whatever transformation the user
 * wants.
 */
int             default_rotate = 0, default_hshift = 0, default_vshift = 0;
float           default_scale = 1., default_xscale = 1., default_yscale = 1.;

/*
 * attribute variables
 */
int             linestyle;
int             cur_color = DEFAULT_COLOR;
int             pat_color = DEFAULT_COLOR + 1;
int             next_color;
struct txalign  txalign;

int             txfont = DEFAULT_FONT;
int             txprec = DEFAULT_PREC;
int             txovly = OVLY_NORMAL;
int             default_txfont, default_txprec, default_txovly;
int             fat = 0;
int             afat = 0;
int             ipat = 0;	/* currently loaded pattern */
struct pat      pat[NPAT + 1];
float           dashsum = 0.;
int             dashon = NO;	/* Dashed lines? */
float           dashes[MAXDASH * 2];
float           dashpos = 0.;	/* Position in current dashing pattern */
int             color_set[MAX_COL + 1][_NUM_PRIM];
extern int      greycorr ();
int             num_col_8;


/*
 * filter options - flags
 */
/* Monochrome device? */
int             mono = NO;
/*
 * Invras determines the polarity of dithered raster.
 * invras=n means raster works the same way as vectors; what is normally
 * WHITE on most devices is BLACK on a hardcopy black and white device.
 * invras=y is the proper default to make dithered images not come out as negatives
 */
int             invras = YES;
int             window = YES;
int             shade = YES;
int             brake = BREAK_BREAK;
int             framewindows = NO;
int             endpause = NO;
int             cachepipe = NO;
/* Setting cachepipe = YES will copy any piped files to a temporary file,
 * this may get done in dev.open.
 * This is useful if the program may want to reverse seek. e.g. Xvpen, Sunpen
 */
int             allowecho = NEVER_SET;
/*
 * setting allowecho NO (which may get done
 * in dev.open) means that the output device
 * is the user's terminal and echoing chars
 * back at the user would insert nonsense
 * into the plot stream.  In this case we
 * explicitly shut down echoing and output
 * translation until the end of the job.
 * Arguably, output translation should be
 * set/unset in dev.open or dev.reset as we may
 * be plotting on a terminal not our own or
 * the device filter may have built in
 * workarounds for known output translation
 * effects. I use it here as a safety measure.
 */
int             wantras = YES;
int             colormask[5] = {1, 1, 1, 1, 1};

/*
 * filter options - enumerated
 */
int             style = NO_STYLE_YET;
int             default_style = STYLE;
int             rotate;
int             size = RELATIVE;
int             erase = FORCE_INITIAL | DO_LITERALS;

/*
 * filter options - valued
 */
int             xcenter, ycenter;	/* Vplot of point to force as center */
int             fatbase = 0;
int             epause = 0;	/* time to pause before erasing screen */
int             overlay = 0;	/* 1=overlay 0=replace */
int             default_overlay;

/*
 * 0 = none
 * 1 = random
 * 2 = Ordered
 * 3 = Floyd-Steinberg
 * 4 = Halftone
 */
int             dither = 1;	/* Dithering type */

int             xWmax, xWmin, yWmax, yWmin;
float           hshift, vshift;	/* Allow global translation of plot */
float           scale;	/* global scale */
float           xscale;	/* global x-scale */
float           yscale;	/* global y-scale */
float           hdevscale;	/* Vplot units to device units for x */
float           vdevscale;	/* Vplot units to device units for y */
float           txscale;/* global text scale */
float           mkscale;/* global marker scale */
float           dashscale;	/* global dashed line scale */
char            interact[MAXFLEN + 1] = "";	/* Where to store coordinate
						 * file */
float           greyc = 1.;	/* Nonlinear correction */
float           pixc = 1.;	/* Pixel overlap correction */
float           redpow = 1., redmap[4] = {1., 0., 0., 0.};
float           greenpow = 1., greenmap[4] = {0., 1., 0., 0.};
float           bluepow = 1., bluemap[4] = {0., 0., 1., 0.};

/* filter options - resettable between plots */
int             user_rotate;
float           user_txscale;
float           user_mkscale;
float           user_dashscale;
float           user_scale;
float           user_xscale;
float           user_yscale;
int             user_size;
float           user_hshift;
float           user_vshift;
int             user_xwmax_flag;
float           user_xwmax;
int             user_ywmax_flag;
float           user_ywmax;
int             user_xwmin_flag;
float           user_xwmin;
int             user_ywmin_flag;
float           user_ywmin;

int             ifat = 0;
float           fatmult_orig;

/*
 * file and terminal control variables
 */
int             pltoutfd, stderrfd, controlfd;
extern int      genmessage ();
int             (*message) () = genmessage;
struct stat     stderrstat;
struct stat     pltoutstat;
FILE           *pltout, *pltin;
FILE           *fopen ();
FILE           *fdopen ();
FILE           *controltty;
char            outbuf[BUFSIZ];
char           *getenv ();
#if defined(HAVE_STDLIB_H)
#include <stdlib.h>
#else
char           *malloc ();
char           *realloc ();
#endif
char            group_name[MAXFLEN + 1];
int             group_number = 0;
FILE           *pltinarray[MAXIN];
char            pltinname[MAXIN][MAXFLEN + 1];
char            pltname[MAXFLEN + 1] = "";
int             infileno = 0;
extern int      gen_do_dovplot ();
int             (*genreader) () = gen_do_dovplot;

/*
 * Initialize and declare global variables.
 * Use getpar and getenv to search for command-line options.
 */

init_vplot ()
{
char           *stringptr;
int             ii;
char            string[MAXFLEN + 1];
float           ftemp;


    txbuffer = (char *) malloc (TXBUFLEN);
    txbuflen = TXBUFLEN;
/*
 * Sun III lint complains about "pointer alignment problem" here,
 * although the documentation makes it sound like that shouldn't
 * be possible.
 */
    vxbuffer = (struct vertex *) malloc (sizeof (struct vertex) * VXBUFLEN);
    vxbuflen = VXBUFLEN;

    /*
     * If the device can't do color, it can set mono to YES If device is
     * mono, this overrides the user's insistence that it isn't. So GETPAR
     * before we call dev.open. 
     */
    GETPAR ("mono", "1", &mono);

    strcpy (wstype, "default");
    if ((stringptr = getenv ("WSTYPE")) != NULL)
	strcpy (wstype, stringptr);
    /*
     * Don't want to accidentally get this from a history file. Must come for
     * the command line, so getpar, not GETPAR 
     */
    getpar ("wstype", "s", wstype);

    /*
     * Initialize all patterns to be undefined. (device can create a
     * device-dependent default set if it wishes) 
     */
    for (ii = 0; ii <= NPAT; ii++)
    {
	pat[ii].patbits = NULL;
    }

/*
 * The device-independent code MAY NOT actually read from
 * the terminal itself. However, these variables are 
 * declared so that the device-dependent routines can use
 * them if they wish as a courtesy.
 */
    controlfd = open ("/dev/tty", 0);
    if (controlfd == OPEN_ERROR)
    {
	controltty = NULL;
    }
    else
    {
	controltty = fdopen (controlfd, "r");
    }

    /*
     * Call device open before doing anything else. this finalizes pltout 
     */
    dev.open ();
    device_open = YES;

/*
 * Beware trying to print out error messages until
 * we've figured out where to connect the "message" routine!
 */

    if (buffer_output)
	setbuf (pltout, outbuf);
    else
	setbuf (pltout, (char *) NULL);

    /*
     * If graphics output going to control terminal, disable echoing and
     * arrange for use of device's message routine 
     */
    pltoutfd = fileno (pltout);
    stderrfd = fileno (stderr);
    out_isatty = isatty (pltoutfd);

    if (allowecho == NEVER_SET)
    {
	allowecho = YES;
	if (out_isatty)
	{
	    fstat (pltoutfd, &pltoutstat);
	    fstat (stderrfd, &stderrstat);
#ifdef SEP
	    if ((pltoutstat.st_dev == stderrstat.st_dev))
		/*
		 * Something in seplib makes the next 2 tests not work, I
		 * don't know why. Unfortunately, this means the SEP versions
		 * aren't as good about detecting the case of sending the
		 * output to a different terminal than the one you're logged
		 * in at. We should probably get this fixed sometime. - Joe
		 * D. 
		 */
#else
	    if ((pltoutstat.st_dev == stderrstat.st_dev) &&
		(pltoutstat.st_ino == stderrstat.st_ino) &&
		(pltoutstat.st_rdev == stderrstat.st_rdev))
#endif /* SEP */
	    {
		allowecho = NO;
		message = dev.message;
	    }
	}
    }

    /*
     * process YES or NO option arguments. GETPAR means it could have come
     * from a seplib history file, getpar means it had to come from the
     * command line. Use "getpar" when having the right value of the
     * parameter is critical because setting the value overrides frontend's
     * judgement. 
     */
    getpar ("echo", "1", &allowecho);
    if (!allowecho)
    {
#if defined(HAVE_TERMIO_H)
	if (ioctl (pltoutfd, TCGETA, &tty_clean_state) == -1)
	{
		ERR (FATAL, name, "Bad ioctl call!");
	}
	tty_plot_state = tty_clean_state;
	tty_plot_state.c_lflag &= ~ECHO;
	if (ioctl (pltoutfd, TCSETAW, &tty_plot_state) == -1)
	{
		ERR (FATAL, name, "Bad ioctl call! (2)");
	}
#else /* USG */
	ioctl (pltoutfd, TIOCGETP, (char *) (&tty_clean_state));
	bcopy ((char *) (&tty_clean_state), (char *) (&tty_plot_state),
	       sizeof (struct sgttyb));
	ioctl (pltoutfd, TIOCLGET, (char *) (&tty_clean_local_mode));
	tty_plot_local_mode = tty_clean_local_mode | LLITOUT;
	ioctl (pltoutfd, TIOCLSET, (char *) (&tty_plot_local_mode));
	tty_plot_state.sg_flags &= (~ECHO);
	tty_plot_state.sg_flags &= (~LCASE);
	tty_plot_state.sg_flags |= (CBREAK);
	ioctl (pltoutfd, TIOCSETN, (char *) (&tty_plot_state));
#endif /* USG */
    }
    getpar ("endpause", "1", &endpause);
    GETPAR ("cachepipe", "1", &cachepipe);
    GETPAR ("shade", "1", &shade);
    GETPAR ("wantras", "1", &wantras);
    getch ("window", "1", &window);
    getch ("frame", "1", &framewindows);
    GETPAR ("overlay", "1", &overlay);
    default_overlay = overlay;
    GETPAR ("invras", "1", &invras);

    GETPAR ("txsquare", "1", &no_stretch_text);

    GETPAR ("colormask", "1", colormask);
    if (colormask[4] == NO)
    {
	if (colormask[0] + colormask[1] + colormask[2] + colormask[3] != 1)
	    ERR (FATAL, name, "Invalid colormask option.");
    }
    if (colormask[4] == NO && colormask[3] == YES)
    {
	colormask[4] = YES;
    }

    GETPAR ("redpow", "f", &redpow);
    GETPAR ("greenpow", "f", &greenpow);
    GETPAR ("bluepow", "f", &bluepow);

    if (redpow <= 0. || greenpow <= 0. || bluepow <= 0.)
	ERR (FATAL, name, "Invalid color pow option.");

    GETPAR ("red", "f", redmap);
    GETPAR ("green", "f", greenmap);
    GETPAR ("blue", "f", bluemap);

/*
 * Valued arguments
 */

    getpar ("dither", "d", &dither);
    getpar ("greyc", "f", &greyc);
    getpar ("pixc", "f", &pixc);

    GETPAR ("txfont", "d", &txfont);
    GETPAR ("txprec", "d", &txprec);
    GETPAR ("txovly", "d", &txovly);
    default_txfont = txfont;
    default_txprec = txprec;
    default_txovly = txovly;

    if (getch ("erase", "s", string))
    {
	if ((string[0] == 'n') || (string[0] == 'N'))
	    erase = NO;
	else
	if ((string[0] == 'o') || (string[0] == 'O'))
	    erase = FORCE_INITIAL;
	else
	if ((string[0] == 'l') || (string[0] == 'L'))
	    erase = DO_LITERALS;
	else
	    erase = FORCE_INITIAL | DO_LITERALS;
    }

    if (getch ("break", "s", string))
    {
	if (string[0] == 'b')
	    brake = BREAK_BREAK;
	else
	if (string[0] == 'e')
	    brake = BREAK_ERASE;
	else
	    brake = BREAK_IGNORE;
    }

    xcenter = 0;
    ycenter = 0;
    xcenterflag = NO;
    ycenterflag = NO;
    if (GETPAR ("xcenter", "f", &ftemp))
    {
	xcenterflag = YES;
	xcenter = ROUND (ftemp * RPERIN);
    }
    if (GETPAR ("ycenter", "f", &ftemp))
    {
	ycenterflag = YES;
	ycenter = ROUND (ftemp * RPERIN);
    }

    if ((stringptr = getenv ("PATTERNMULT")) != NULL)
    {
	sscanf (stringptr, "%f", &patternmult);
    }
    GETPAR ("patternmult", "f", &patternmult);

    GETPAR ("interact", "s", interact);
    GETPAR ("pause", "d", &epause);

    if (interact[0] != '\0')
    {
	epause = 0;		/* interact makes it own sort of pausing */
	endpause = NO;
    }

    /*
     * Find the default style 
     */
    stringptr = NULL;
    if (GETPAR ("style", "s", string))
	stringptr = string;
    else
	stringptr = getenv ("PLOTSTYLE");
    if (stringptr != NULL)
    {
	if ((stringptr[0] == 'r') || (stringptr[0] == 'R') ||
	    (stringptr[0] == 'm') || (stringptr[0] == 'M'))
	    default_style = ROTATED;
	else
	if ((stringptr[0] == 'o') || (stringptr[0] == 'O'))
	    default_style = OLD;
	else
	if ((stringptr[0] == 'a') || (stringptr[0] == 'A'))
	    default_style = ABSOLUTE;
	else
	    default_style = STANDARD;
    }

/*
 * Options changeable between calls to dovplot
 * (By calling reset_parameters or reset())
 * Dovplot calls reset every time it starts. It only calls
 * reset_parameters the first time.
 *
 * Things in this category:
 * user_*
 * fatmult_orig, fatbase
 */

    if ((stringptr = getenv ("FATMULT")) != NULL)
    {
	sscanf (stringptr, "%f", &fatmult);
    }
    GETPAR ("fatmult", "f", &fatmult);
    fatmult_orig = fatmult;

    user_rotate = 0;
    getch ("rotate", "d", &user_rotate);

    user_txscale = 1.0;
    user_mkscale = 1.0;
    user_dashscale = 1.0;
    GETPAR ("txscale", "f", &user_txscale);
    GETPAR ("mkscale", "f", &user_mkscale);
    GETPAR ("dashscale", "f", &user_dashscale);

    user_scale = 1.0;
    user_xscale = 1.0;
    user_yscale = 1.0;
    getch ("scale", "f", &user_scale);
    getch ("xscale", "f", &user_xscale);
    getch ("yscale", "f", &user_yscale);

    user_size = size;
    if (getch ("size", "s", string))
    {
	if ((string[0] == 'a') || (string[0] == 'A'))
	    user_size = ABSOLUTE;
	else
	    user_size = RELATIVE;
    }

    user_hshift = 0.;
    getch ("hshift xshift", "f", &user_hshift);
    user_vshift = 0.;
    getch ("vshift yshift", "f", &user_vshift);

    user_xwmax_flag = GETPAR ("xwmax", "f", &user_xwmax);
    user_ywmax_flag = GETPAR ("ywmax", "f", &user_ywmax);
    user_xwmin_flag = GETPAR ("xwmin", "f", &user_xwmin);
    user_ywmin_flag = GETPAR ("ywmin", "f", &user_ywmin);

    GETPAR ("fat", "d", &fatbase);



/*
 *  These parameters can simply be changed at any time with no
 *  need for re-initialization of anything else:
 *
 *  shade, wantras
 *
 */

}

init_colors ()
{
int             ii;

    /*
     * Set up the default color table 
     */

/*
 * First 7 colors are assumed to map to the standard 7 pen colors
 * on ALL devices, even those with no settable colors!
 */
    for (ii = 0; ii < 8; ii++)
    {
/*     if(color_set[ii][STATUS]!=SET){*/
	    color_set[ii][STATUS] = SET;
	    color_set[ii][_RED] = MAX_GUN * ((ii & 2) / 2);
	    color_set[ii][_GREEN] = MAX_GUN * ((ii & 4) / 4);
	    color_set[ii][_BLUE] = MAX_GUN * ((ii & 1) / 1);
	    color_set[ii][_GREY] = greycorr ((int) ((MAX_GUN * ii) / 7));
/*   }*/
    }

    if (mono)
    {
	/* Monochrome devices are assumed to have no settable colors */
	num_col = 0;
    }

    num_col_8 = (num_col > 8) ? num_col : 8;

    if (mono)
    {
	color_set[0][MAP] = 0;
	for (ii = 1; ii <= MAX_COL; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][MAP] = 7;
	}
	for (ii = num_col_8; ii < 256; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][_GREY] = color_set[((ii - 8) % 7) + 1][_GREY];
	}
	/* Put a grey scale in the upper half of the color table */
	for (ii = 256; ii <= MAX_COL; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][_GREY] = greycorr (ii - 256);
	}
    }
    else
    {
	/*
	 * Unmapped colors shouldn't be mapped; ie, they map to themselves 
	 */

	for (ii = 0; ii < num_col_8; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][MAP] = ii;
/*  fprintf(stderr,"MAPPING %d %d - %d \n",ii,MAP,color_set[ii][MAP]);*/
	 }
	/*
	 * Colors outside the range of this terminal map cyclically back into
	 * colors 1 through 7 
	 */
	for (ii = num_col_8; ii <= MAX_COL; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][MAP] = ((ii - 8) % 7) + 1;
	}
    }
}

setstyle (new_style)
    int             new_style;
{
    /*
     * Check to see if the style has changed 
     */
    if (new_style == style)
	return;

    style = new_style;
    reset_parameters ();
}

reset_parameters ()
{
float           inches;	/* scaling base for y axis */
float           screenheight, screenwidth;	/* true size of the screen */
int             ix, iy;

    xorigin = 0;
    yorigin = 0;

    rotate = default_rotate;
    rotate += user_rotate;

    txscale = user_txscale;
    mkscale = user_mkscale;
    dashscale = user_dashscale;

    scale = default_scale * user_scale;
    xscale = default_xscale * user_xscale;
    yscale = default_yscale * user_yscale;

    fatmult = fatmult_orig;

    size = user_size;

    switch (style)
    {
/*
 * The old standard on the machine erebus.
 * Pretty much dead now. Still useful for some old programs nobody's
 * wanted to update.
 */
    case OLD:
	txscale /= 3.;
	fatmult /= 3.;
	scale *= 3.;
	inches = STANDARD_HEIGHT;
	break;
/*
 * The old standard on the machine mazama. A useful coordinate system
 * for geological sorts of plots. 
 * The Y axis goes the long way, across the screen, which is the device's
 * horizontal axis.
 */
    case ROTATED:
	rotate += 90;
	inches = ROTATED_HEIGHT;
	break;
    case ABSOLUTE:
	size = ABSOLUTE;
    case STANDARD:
    default:
	inches = STANDARD_HEIGHT;
	break;
    }

    if (rotate >= 0)
	rotate = rotate % 360;
    else
	rotate = ((rotate % 360) + 360) % 360;

    mxx = cos (2. * 3.14159 * rotate / 360.);
    myy = cos (2. * 3.14159 * rotate / 360.);
    mxy = sin (2. * 3.14159 * rotate / 360.);
    myx = -sin (2. * 3.14159 * rotate / 360.);

    if (size == ABSOLUTE)
    {
	vdevscale = pixels_per_inch / (float) (RPERIN * aspect_ratio);
	hdevscale = pixels_per_inch / (float) RPERIN;
    }
    else
    {
	/*
	 * Fit the inches x inches unit square into a displayable box with
	 * aspect ratio SCREEN_RATIO 
	 */
	screenwidth = (dev_xmax - dev_xmin) * pixels_per_inch;
	screenheight =
	 (dev_ymax - dev_ymin) * pixels_per_inch * aspect_ratio;
	if ((screenheight / screenwidth) > SCREEN_RATIO)
	{
	    vdevscale = (SCREEN_RATIO * ((dev_xmax - dev_xmin) / aspect_ratio)) /
	     (inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	}
	else
	{
	    vdevscale = (dev_ymax - dev_ymin) / (inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	}
    }

    hshift = default_hshift;
    vshift = default_vshift;

    if (style == ROTATED)
    {
	vshift += dev_ymax - dev_ymin;
    }

    yscale *= scale;
    xscale *= scale;
    mkscale *= scale;

/*
 * Set up fatness multiplication factor
 */
    fatmult *= scale * hdevscale * RPERIN / FATPERIN;

    /*
     * The point (xcenter,ycenter) in vplot coordinates is to be centered in
     * the screen. 
     */
    if (xcenterflag || ycenterflag)
    {
	vptodevxy (xcenter, ycenter, &ix, &iy);
	if (xcenterflag)
	{
	    hshift += (dev_xmax + dev_xmin) / 2 - ix;
	}
	if (ycenterflag)
	{
	    vshift += (dev_ymax + dev_ymin) / 2 - iy;
	}
    }

    hshift += user_hshift * pixels_per_inch;
    vshift += user_vshift * pixels_per_inch / aspect_ratio;

/* plot window parameters defaulted */
/* to maximum size */

    devtovpw (dev_xmin, dev_ymin, dev_xmax, dev_ymax,
	      &xWmin, &yWmin, &xWmax, &yWmax);

    if (user_xwmax_flag)
	xWmax = ROUND (user_xwmax * RPERIN);
    if (user_ywmax_flag)
	yWmax = ROUND (user_ywmax * RPERIN);
    if (user_xwmin_flag)
	xWmin = ROUND (user_xwmin * RPERIN);
    if (user_ywmin_flag)
	yWmin = ROUND (user_ywmin * RPERIN);

    vptodevw (xWmin, yWmin, xWmax, yWmax, &xWmin, &yWmin, &xWmax, &yWmax);

    wlimit (dev_xmin, dev_xmax, &xWmin, &xWmax);
    wlimit (dev_ymin, dev_ymax, &yWmin, &yWmax);

    xwmax = xWmax;		/* plot window parameters defaulted */
    xwmin = xWmin;		/* to maximum size 		 */
    ywmax = yWmax;
    ywmin = yWmin;
    reset_windows ();
}

add_a_cor (filename, xcor, ycor)
    char           *filename;
    int             xcor, ycor;
{
static int      first_time = YES;
static FILE    *outfp;

    if (first_time == YES)
    {
	outfp = fopen (filename, "w");
	if (outfp == NULL)
	{
	    ERR (FATAL, name, "Can't open interact output file %s!", filename);
	}
	first_time = NO;
    }
    fprintf (outfp, "%f\t%f\n", (float) xcor / RPERIN, (float) ycor / RPERIN);
}

wlimit (min, max, xmin, xmax)
    int             min, max;
    int            *xmin, *xmax;
{
    if (*xmin < min)
	*xmin = min;

    if (*xmax > max)
	*xmax = max;
}
